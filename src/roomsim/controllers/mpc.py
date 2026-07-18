"""Economic model-predictive control, supervising a PI, with a Kalman filter.

Planner. Every Ts_pred (15 min) the controller solves, over a 24 h horizon with steps of Ts_pred,

    min   sum_j  c_j u_j  +  rho sum_j eps_j  +  lam sum_j |u_j - u_j-1|
    s.t.  T_lo,j - eps_j  <=  T_a,j+1  <=  T_hi,j + eps_j        (comfort band, soft)
          0 <= u_j <= U_max,j,   eps_j >= 0

where u_j is the sensible cooling (W), c_j the electrical energy per unit of sensible cooling
forecast for that step (the unit's COP depends on the outdoor temperature), times the price
relative to off-peak when a time-of-use tariff is given, and T_a,j+1 is
affine in u through the linear RC model, the forecast gains and the disturbance estimated by
the Kalman filter. The comfort band is narrow while the office is occupied and wide otherwise,
so the optimiser can let the room drift at night and pre-cool with the thermal mass when that is
cheaper. All terms are linear, so it is a linear programme (HiGHS, ~10 ms).

Supervision. The planned air temperature at the end of the first step becomes the set point of
the SIMC-tuned PI (with sigma-delta low-load handling), which runs every 60 s. Applying the
planned cooling directly, with the low-load handler converting it into cycling between solves,
was tried first: with no feedback between solves the cycling at 20-30 % demand produced
on/off periods of up to an hour and 3.75 K h outside the comfort band over three days. The PI
closes that loop; the planner only decides where the temperature should be. The occupied band
seen by the planner is narrowed by ``margin`` to leave room for the PI's ripple.

Estimation. The Kalman filter runs at every controller sample on the air-temperature sensor.
The cooling it assumes for the last sample is computed with the unit model at the measured air
state, so the disturbance state only has to absorb what the models really get wrong.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import linprog

from .. import psychro
from ..ac_unit import ACUnit
from ..control_model import RCModel
from ..estimation import DisturbanceKalman
from ..inputs import Scenario, hour_of_day
from ..plant import IDX_Y_RH, IDX_Y_TA
from ..tariff import TimeOfUseTariff
from .base import NO_VENTILATION, CompressorGuard, VentilationSchedule
from .pi_unit import PIUnitController, make_pi_controller

HOUR = 3600.0


@dataclass(frozen=True)
class ComfortSchedule:
    """Comfort band by time of day, and the matching set point for a PI.

    Occupied: T_set +- band. Unoccupied: up to T_hi_unocc (set back). ``early_start_h`` moves the
    PI set point back to T_set before occupancy (a fixed pre-cooling start).
    """

    T_set: float = 25.0
    band: float = 1.0
    T_hi_unocc: float = 28.0
    T_lo_unocc: float = 20.0
    start_h: float = 8.0
    end_h: float = 18.0
    early_start_h: float | None = None

    def occupied(self, t) -> np.ndarray:
        h = hour_of_day(t)
        return (h >= self.start_h) & (h < self.end_h)

    def occupied_at_instant(self, t) -> np.ndarray:
        """Occupied on either side of instant t: the state at 8:00 and at 18:00 both count."""
        t = np.asarray(t, dtype=float)
        return self.occupied(t) | self.occupied(t - 1.0)

    def bounds(self, t) -> tuple[np.ndarray, np.ndarray]:
        """Comfort band for the air temperature at instant t."""
        occ = self.occupied_at_instant(t)
        return (np.where(occ, self.T_set - self.band, self.T_lo_unocc),
                np.where(occ, self.T_set + self.band, self.T_hi_unocc))

    def setpoint(self, t: float) -> float:
        """Set point for a PI: T_set from the (early) start to the end of occupancy."""
        h = float(hour_of_day(t))
        start = self.start_h if self.early_start_h is None else self.early_start_h
        return self.T_set if start <= h < self.end_h else self.T_hi_unocc


class EconomicMPC:
    def __init__(self, model: RCModel, unit: ACUnit, forecast: Scenario,
                 comfort: ComfortSchedule = ComfortSchedule(),
                 ventilation: VentilationSchedule = NO_VENTILATION,
                 Ts: float = 60.0, Ts_pred: float = 900.0, horizon_h: float = 24.0,
                 rho: float = 5.0, lam: float = 2e-5, margin: float = 0.6,
                 W_ref: float = 0.0115, kalman: DisturbanceKalman | None = None,
                 tariff: TimeOfUseTariff | None = None,
                 inner: PIUnitController | None = None, **pi_kwargs):
        self.model, self.unit, self.forecast = model, unit, forecast
        self.comfort, self.ventilation = comfort, ventilation
        self.Ts, self.Ts_pred = Ts, Ts_pred
        self.N = int(round(horizon_h * HOUR / Ts_pred))
        self.n_pred = int(round(Ts_pred / Ts))
        self.rho, self.lam, self.margin, self.W_ref = rho, lam, margin, W_ref
        self.tariff = tariff
        self.kf = kalman or DisturbanceKalman(model, Ts)
        self.T_ref = comfort.T_set
        # the inner PI is tuned on the controller's own model parameters
        self.inner = inner or make_pi_controller("sigma-delta", model.p, unit, Ts=Ts,
                                                 T_set=lambda t: self.T_ref,
                                                 ventilation=ventilation, **pi_kwargs)
        self.k = 0
        self.u_plan = 0.0               # sensible cooling of the current planned move, W
        self.Q_c_prev = 0.0             # sensible cooling estimate over the last sample, W
        self.t_prev = None
        self.V_prev = 0.0
        self.last_solution = None
        self.solve_failures = 0

    @property
    def guard(self) -> CompressorGuard:
        return self.inner.guard

    # -- unit characteristics as the controller sees them -----------------------------------------
    def sensible_capacity(self, T_o: float) -> float:
        """Sensible cooling at full speed and full fan, at the comfort set point."""
        return self.unit.operate(1.0, 1.0, self.comfort.T_set, self.W_ref, T_o).Q_sens

    def energy_per_cooling(self, T_o: float) -> float:
        """Electrical energy per unit of sensible cooling at a typical part load, -."""
        st = self.unit.operate(0.5, 1.0, self.comfort.T_set, self.W_ref, T_o)
        return st.P_el / st.Q_sens

    def planning_band(self, t) -> tuple[np.ndarray, np.ndarray]:
        """Comfort band narrowed by ``margin`` while occupied, for the PI ripple."""
        lo, hi = self.comfort.bounds(t)
        occ = self.comfort.occupied_at_instant(t)
        return np.where(occ, lo + self.margin, lo), np.where(occ, hi - self.margin, hi)

    # -- optimisation ------------------------------------------------------------------------------
    def plan(self, t: float, z: np.ndarray) -> np.ndarray:
        """Solve the LP from estimated state z = [T_a, T_m, d] at time t; returns u (W) per step."""
        N, Ts = self.N, self.Ts_pred
        t_start = t + Ts * np.arange(N)
        t_mid = t_start + 0.5 * Ts
        d = self.forecast(t_mid)
        e = self.model.gains(d)
        T_o = np.asarray(d.T_o, dtype=float) * np.ones(N)
        V = np.array([self.ventilation(tt) for tt in t_start])
        lo, hi = self.planning_band(t_start + Ts)

        # affine prediction of the air temperature at the end of each step: T = f + S u
        S = np.zeros((N, N))
        f = np.zeros(N)
        X = np.zeros((2, N))
        x = z[:2].copy()
        for j in range(N):
            Ad, Bd, Ed = self.model.discrete(Ts, V[j])
            X = Ad @ X
            X[:, j] += Bd
            x = Ad @ x + Ed @ e[j] - Bd * z[2]           # d acts as heat into the air
            S[j] = X[0]
            f[j] = x[0]

        U_max = np.array([self.sensible_capacity(To) for To in T_o])
        c = np.array([self.energy_per_cooling(To) for To in T_o]) * Ts / 3.6e6   # kWh per W
        if self.tariff is not None:
            # price relative to the reference, so rho keeps its meaning in equivalent kWh
            c = c * self.tariff(t_mid) / self.tariff.reference
        # variables [u (N), eps (N), delta (N)]
        cost = np.concatenate([c, np.full(N, self.rho * Ts / HOUR), np.full(N, self.lam)])
        I, Z = np.eye(N), np.zeros((N, N))
        D = I - np.eye(N, k=-1)                              # u_j - u_j-1
        A_ub = np.block([[-S, -I, Z], [S, -I, Z], [D, Z, -I], [-D, Z, -I]])
        prev = np.zeros(N)
        prev[0] = self.u_plan
        b_ub = np.concatenate([f - lo, hi - f, prev, -prev])
        bounds = [(0.0, um) for um in U_max] + [(0.0, None)] * (2 * N)
        res = linprog(cost, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        if res.status != 0:
            self.solve_failures += 1
            return np.full(N, self.u_plan)
        u = res.x[:N]
        self.last_solution = dict(t=t_start, u=u, T=f + S @ u, lo=lo, hi=hi, eps=res.x[N:2 * N],
                                  U_max=U_max, cost=res.fun)
        return u

    # -- controller interface ----------------------------------------------------------------------
    def step(self, t: float, y: np.ndarray) -> np.ndarray:
        T_a = y[IDX_Y_TA]
        if self.t_prev is None:
            self.kf.initialise(T_a)
        else:
            self.kf.update(T_a, self.Q_c_prev, self.forecast(self.t_prev), self.V_prev)

        if self.k % self.n_pred == 0:
            u = self.plan(t, self.kf.z)
            self.u_plan = float(u[0])
            if self.last_solution is not None:
                # planned temperature at the end of the move, kept inside the band in force now
                lo, hi = self.planning_band(t)
                self.T_ref = float(np.clip(self.last_solution["T"][0], lo, hi))
        self.k += 1

        cmd = self.inner.step(t, y)
        # cooling actually requested, evaluated with the unit model at the measured air state
        s = cmd[0]
        if s > 0:
            W = float(psychro.W_from_RH(T_a, y[IDX_Y_RH]))
            st = self.unit.operate(s, 1.0, T_a, W, float(self.forecast(t).T_o))
            self.Q_c_prev = st.Q_sens - st.P_fan_in
        else:
            self.Q_c_prev = 0.0
        self.t_prev, self.V_prev = t, cmd[2]
        return cmd
