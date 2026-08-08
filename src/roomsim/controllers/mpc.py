"""Economic model-predictive control, supervising a PI, with a Kalman filter.

Planner. Every Ts_pred (15 min) the controller solves, over a 24 h horizon with steps of Ts_pred,

    min   sum_j  c_j u_j  +  rho sum_j eps_j  +  lam sum_j |u_j - u_j-1|
    s.t.  T_lo,j - eps_j  <=  T_a,j+1  <=  T_hi,j + eps_j        (comfort band, soft)
          0 <= u_j <= U_max,j,   eps_j >= 0

where u_j is the sensible cooling net of the indoor-fan heat (W), c_j the electrical energy per
unit of that cooling forecast for that step (the unit's COP depends on the outdoor temperature),
times the price relative to off-peak when a time-of-use tariff is given, and T_a,j+1 is
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
from typing import Any, TypedDict, cast

import numpy as np
from scipy.optimize import linprog

from .. import psychro
from .._typing import BoolArray, FloatArray
from ..ac_unit import ACUnit
from ..control_model import RCModel
from ..estimation import DisturbanceKalman
from ..inputs import Scenario, hour_of_day
from ..plant import IDX_Y_RH, IDX_Y_TA
from ..tariff import TimeOfUseTariff, check_daily_window, in_daily_window
from .base import NO_VENTILATION, CompressorGuard, VentilationSchedule, sampling_period
from .pi_unit import PIUnitController, make_pi_controller

HOUR = 3600.0

"""Array of booleans (alias used for annotations)."""


class PlanSolution(TypedDict):
    """Last solution of the planning LP, kept for inspection.

    Attributes
    ----------
    t : ndarray
        Start time of each planning step, s.
    u : ndarray
        Planned sensible cooling per step, W.
    T : ndarray
        Predicted air temperature at the end of each step, degC.
    lo, hi : ndarray
        Planning band at the end of each step, degC.
    eps : ndarray
        Comfort-band violation per step, K.
    U_max : ndarray
        Net sensible capacity per step, W.
    cost : float
        Optimal objective value, equivalent kWh.
    """

    t: FloatArray
    u: FloatArray
    T: FloatArray
    lo: FloatArray
    hi: FloatArray
    eps: FloatArray
    U_max: FloatArray
    cost: float


@dataclass(frozen=True)
class ComfortSchedule:
    """Comfort band by time of day, and the matching set point for a PI.

    Occupied: T_set +- band. Unoccupied: up to T_hi_unocc (set back). ``early_start_h`` moves the
    PI set point back to T_set before occupancy (a fixed pre-cooling start). A period whose start
    is later than its end crosses midnight.

    Attributes
    ----------
    T_set : float
        Occupied set point, degC.
    band : float
        Half-width of the occupied comfort band, K.
    T_hi_unocc, T_lo_unocc : float
        Upper and lower limits while unoccupied, degC.
    start_h, end_h : float
        Start and end of occupancy, clock hours (h).
    early_start_h : float or None
        Clock hour (h) from which the PI set point is T_set; None uses ``start_h``.

    Raises
    ------
    ValueError
        If a clock hour is outside [0, 24], or start_h or early_start_h equals end_h (see
        ``roomsim.tariff.check_daily_window``).
    """

    T_set: float = 25.0
    band: float = 1.0
    T_hi_unocc: float = 28.0
    T_lo_unocc: float = 20.0
    start_h: float = 8.0
    end_h: float = 18.0
    early_start_h: float | None = None

    def __post_init__(self) -> None:
        """Validate the occupied period and the early start."""
        check_daily_window("start_h", self.start_h, "end_h", self.end_h)
        if self.early_start_h is not None:
            check_daily_window("early_start_h", self.early_start_h, "end_h", self.end_h)

    def occupied(self, t: float | FloatArray) -> np.bool_ | BoolArray:
        """Return whether time t (s) lies in [start_h, end_h) of its day; elementwise."""
        return in_daily_window(hour_of_day(t), self.start_h, self.end_h)

    def occupied_at_instant(self, t: float | FloatArray) -> np.bool_ | BoolArray:
        """Return whether instant t (s) is occupied on either side; elementwise.

        The state at 8:00 and at 18:00 both count as occupied.
        """
        t = np.asarray(t, dtype=float)
        return self.occupied(t) | self.occupied(t - 1.0)

    def bounds(self, t: float | FloatArray) -> tuple[FloatArray, FloatArray]:
        """Return the comfort band (lower, upper) for the air temperature at instant t.

        Parameters
        ----------
        t : float or ndarray
            Time, s.

        Returns
        -------
        lo, hi : ndarray
            Lower and upper limits, degC.
        """
        occ = self.occupied_at_instant(t)
        return (np.where(occ, self.T_set - self.band, self.T_lo_unocc),
                np.where(occ, self.T_set + self.band, self.T_hi_unocc))

    def setpoint(self, t: float) -> float:
        """Return the PI set point at time t (s), degC.

        T_set from the (early) start to the end of occupancy, ``T_hi_unocc`` otherwise.
        """
        h = float(hour_of_day(t))
        start = self.start_h if self.early_start_h is None else self.early_start_h
        return self.T_set if in_daily_window(h, start, self.end_h) else self.T_hi_unocc


class EconomicMPC:
    """Economic MPC that plans the air temperature and supplies the set point of an inner PI.

    Parameters
    ----------
    model : RCModel
        Controller model (also used to tune the default inner PI).
    unit : ACUnit
        Air-conditioning unit model used for capacity, efficiency and the cooling estimate.
    forecast : Scenario
        Disturbance forecast as a function of time (s).
    comfort : ComfortSchedule, optional
        Comfort band and occupancy.
    ventilation : VentilationSchedule, optional
        Outdoor-air flow schedule, assumed known to the planner.
    Ts : float, optional
        Controller sampling period, s.
    Ts_pred : float, optional
        Planning step and re-planning interval, s.
    horizon_h : float, optional
        Planning horizon, h.
    rho : float, optional
        Weight of comfort violation, equivalent kWh per K h.
    lam : float, optional
        Weight of move size, per W of change in u.
    margin : float, optional
        Narrowing of the occupied band on each side, K.
    W_ref : float, optional
        Return-air humidity ratio assumed for capacity and efficiency, kg/kg.
    kalman : DisturbanceKalman or None, optional
        State estimator; None builds one from ``model`` and ``Ts``.
    tariff : TimeOfUseTariff or None, optional
        Electricity price; None weights energy uniformly.
    inner : PIUnitController or None, optional
        Inner PI; None builds a sigma-delta PI tuned on ``model.p``.
    **pi_kwargs
        Passed to ``make_pi_controller`` when ``inner`` is None.

    Attributes
    ----------
    N : int
        Number of planning steps.
    n_pred : int
        Controller samples per planning step.
    T_ref : float
        Current set point given to the inner PI, degC.
    u_plan : float
        Sensible cooling of the current planned move, W.
    Q_c_prev : float
        Sensible cooling estimate over the last sample, W.
    last_solution : PlanSolution or None
        Last successful LP solution.
    solve_failures : int
        Number of LP solves that did not succeed.

    Raises
    ------
    ValueError
        If ``Ts_pred`` is not a positive integer multiple of ``Ts``, or the Kalman filter or
        the inner PI assume a sampling period other than ``Ts``.
    """

    def __init__(self, model: RCModel, unit: ACUnit, forecast: Scenario,
                 comfort: ComfortSchedule = ComfortSchedule(),  # noqa: B008  (frozen dataclass)
                 ventilation: VentilationSchedule = NO_VENTILATION,
                 Ts: float = 60.0, Ts_pred: float = 900.0, horizon_h: float = 24.0,
                 rho: float = 5.0, lam: float = 2e-5, margin: float = 0.6,
                 W_ref: float = 0.0115, kalman: DisturbanceKalman | None = None,
                 tariff: TimeOfUseTariff | None = None,
                 inner: PIUnitController | None = None, **pi_kwargs: Any) -> None:
        self.model, self.unit, self.forecast = model, unit, forecast
        self.comfort, self.ventilation = comfort, ventilation
        self.Ts, self.Ts_pred = Ts, Ts_pred
        self.N = int(round(horizon_h * HOUR / Ts_pred))
        self.n_pred = int(round(Ts_pred / Ts))
        if self.n_pred < 1 or abs(self.n_pred * Ts - Ts_pred) > 1e-9 * Ts_pred:
            raise ValueError(f"Ts_pred = {Ts_pred!r} s is invalid; it must be a positive integer "
                             f"multiple of Ts = {Ts!r} s")
        for part, Ts_part in (("kalman", sampling_period(kalman)),
                              ("inner", sampling_period(inner))):
            if Ts_part is not None and abs(Ts_part - Ts) > 1e-9 * Ts:
                raise ValueError(f"{part}.Ts = {Ts_part!r} s does not match the MPC's sampling "
                                 f"period Ts = {Ts!r} s; they must be equal")
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
        self.t_prev: float | None = None
        self.V_prev = 0.0
        self.last_solution: PlanSolution | None = None
        self.solve_failures = 0

    @property
    def guard(self) -> CompressorGuard:
        """Compressor guard of the inner PI."""
        return self.inner.guard

    # -- unit characteristics as the controller sees them -----------------------------------------
    def sensible_capacity(self, T_o: float) -> float:
        """Return the net sensible cooling at full speed and full fan, at the comfort set point.

        Net of the indoor-fan heat, which stays in the room: the same quantity as the cooling
        fed to the Kalman filter, the input Q_c of the RC model.

        Parameters
        ----------
        T_o : float
            Outdoor temperature, degC.

        Returns
        -------
        float
            Net sensible cooling Q_sens - P_fan_in, W.
        """
        st = self.unit.operate(1.0, 1.0, self.comfort.T_set, self.W_ref, T_o)
        return st.Q_sens - st.P_fan_in

    def energy_per_cooling(self, T_o: float) -> float:
        """Return the electrical energy per unit of net sensible cooling at half speed, -.

        Parameters
        ----------
        T_o : float
            Outdoor temperature, degC.

        Returns
        -------
        float
            P_el / (Q_sens - P_fan_in), -.
        """
        st = self.unit.operate(0.5, 1.0, self.comfort.T_set, self.W_ref, T_o)
        return st.P_el / (st.Q_sens - st.P_fan_in)

    def planning_band(self, t: float | FloatArray) -> tuple[FloatArray, FloatArray]:
        """Return the comfort band narrowed by ``margin`` while occupied, for the PI ripple.

        Parameters
        ----------
        t : float or ndarray
            Time, s.

        Returns
        -------
        lo, hi : ndarray
            Lower and upper limits, degC.
        """
        lo, hi = self.comfort.bounds(t)
        occ = self.comfort.occupied_at_instant(t)
        return np.where(occ, lo + self.margin, lo), np.where(occ, hi - self.margin, hi)

    # -- optimisation ------------------------------------------------------------------------------
    def plan(self, t: float, z: FloatArray) -> FloatArray:
        """Solve the planning LP from the estimated state at time t.

        On solver failure ``solve_failures`` is incremented and the last successful plan is
        returned shifted to the current step (padded with its last move); without one, the
        current move is repeated.

        Parameters
        ----------
        t : float
            Current time, s.
        z : ndarray, shape (3,)
            Estimated state [T_a (degC), T_m (degC), d (W)].

        Returns
        -------
        ndarray, shape (N,)
            Planned sensible cooling per planning step, W.
        """
        N, Ts = self.N, self.Ts_pred
        t_start = cast("FloatArray", t + Ts * np.arange(N))   # float64; stubs say floating[Any]
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
        Id, Z = np.eye(N), np.zeros((N, N))
        D = Id - np.eye(N, k=-1)                             # u_j - u_j-1
        A_ub = np.block([[-S, -Id, Z], [S, -Id, Z], [D, Z, -Id], [-D, Z, -Id]])
        prev = np.zeros(N)
        prev[0] = self.u_plan
        b_ub = np.concatenate([f - lo, hi - f, prev, -prev])
        bounds = [(0.0, um) for um in U_max] + [(0.0, None)] * (2 * N)
        res = linprog(cost, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        if res.status != 0:
            self.solve_failures += 1
            j_now = self._plan_index(t)
            if j_now is None or self.last_solution is None:
                return np.full(N, self.u_plan)
            old = self.last_solution["u"][j_now:]
            return np.concatenate([old, np.full(N - len(old), old[-1])])
        u: FloatArray = res.x[:N]
        self.last_solution = {"t": t_start, "u": u, "T": f + S @ u, "lo": lo, "hi": hi,
                              "eps": res.x[N:2 * N], "U_max": U_max, "cost": res.fun}
        return u

    def _plan_index(self, t: float) -> int | None:
        """Return the step of the last successful plan that starts at time t (s).

        0 right after a successful solve, k after k failed solves; None without a plan or when
        t lies beyond its horizon.
        """
        if self.last_solution is None:
            return None
        j = int(round((t - float(self.last_solution["t"][0])) / self.Ts_pred))
        return j if 0 <= j < len(self.last_solution["T"]) else None

    # -- controller interface ----------------------------------------------------------------------
    def step(self, t: float, y: FloatArray) -> FloatArray:
        """Return the command for one controller sample.

        Updates the Kalman filter, re-plans every ``n_pred`` samples, runs the inner PI and
        estimates the sensible cooling requested, for the filter's next prediction.

        Parameters
        ----------
        t : float
            Current time, s.
        y : ndarray
            Measurements; uses the air temperature (degC) and relative humidity (-).

        Returns
        -------
        ndarray, shape (3,)
            Command [s (-), f (-), V_vent (m3/s)].
        """
        T_a = y[IDX_Y_TA]
        if self.t_prev is None:
            self.kf.initialise(T_a)
        else:
            self.kf.update(T_a, self.Q_c_prev, self.forecast(self.t_prev), self.V_prev)

        if self.k % self.n_pred == 0:
            u = self.plan(t, self.kf.state)
            self.u_plan = float(u[0])
            j = self._plan_index(t)
            if j is not None and self.last_solution is not None:
                # planned temperature at the end of the current move, kept inside the band that
                # applies at that instant (the band the LP used for it)
                lo, hi = self.planning_band(t + self.Ts_pred)
                self.T_ref = float(np.clip(self.last_solution["T"][j], lo, hi))
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
