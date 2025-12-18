"""Fixed-step simulation loop: continuous plant, discrete controller, logging.

Per step [t_k, t_k + dt):
  1. if t_k is a controller sampling instant, read the sensors and update the controller
     output (then hold it until the next sample);
  2. sample the piecewise-constant disturbances at t_k;
  3. advance the plant with Euler or RK4; continuous disturbances are evaluated at the stage times.

Without an air-conditioning unit the controller output is the plant input u directly.
With a unit, the controller output is a command [s, f, V_vent] (compressor speed, indoor fan,
ventilation) and the unit converts it into u at every Runge-Kutta stage, using the return-air
state of that stage. Holding Q_sens itself over the sampling interval instead would make the
closed loop only first-order accurate, because the coil output depends on the air state.

With a unit, the water held on its coil, M_c, is integrated as an additional state.
Cumulative boundary heat, water and CO2, and the electrical energy of the unit, are integrated
as extra states with the same method, so the balances in ``SimResult`` reflect exactly what the
integrator did.
"""

from __future__ import annotations

from typing import Protocol

import numpy as np

from .ac_unit import OFF_STATE, ACUnit, CoilState
from .inputs import PIECEWISE_CONSTANT, Disturbances, Scenario
from .integrators import STEPPERS
from .params import RoomParams
from .plant import (IDX_M_LAT, IDX_Q_SENS, IDX_TA, IDX_V_VENT, IDX_W, N_INPUTS, N_STATES,
                    STATE_SCALE, Fluxes, derivatives, fluxes, measure)
from .results import SimResult

# Command vector when a unit is present
IDX_CMD_S, IDX_CMD_FAN, IDX_CMD_VENT = 0, 1, 2
N_CMD = 3


class Controller(Protocol):
    def step(self, t: float, y: np.ndarray) -> np.ndarray:
        """Output from the measurement y = [T_a, RH_a, C], plus the coil temperature with a unit.

        Without a unit: u = [Q_hvac_sens, m_hvac_lat, V_vent].
        With a unit:    command = [s, f, V_vent].
        """
        ...


# Cumulative integrals: boundary heat (J), room water (kg), CO2 (m3), electrical energy (J),
# water entering from outside (kg) and condensate drained (kg)
N_AUX = 6
UNIT_LOG = ("Q_sens", "m_lat", "Q_coil", "P_el", "COP", "T_s", "T_out", "m_drain")


def _ratio(a: float, b: float, what: str) -> int:
    n = round(a / b)
    if n < 1 or abs(n * b - a) > 1e-9 * max(abs(a), 1.0):
        raise ValueError(f"{what} must be a positive integer multiple of the step")
    return n


def coil_sensor(st: CoilState, T_a: float) -> float:
    """Evaporator thermistor: coil surface while the unit runs, room air when it is idle."""
    return st.T_s if (st.Q_coil > 0.0 or st.P_fan_in > 0.0) and st.T_s == st.T_s else T_a


def unit_input(unit: ACUnit, cmd: np.ndarray, ramp: float, x: np.ndarray,
               d: Disturbances, M_c: float | None = None) -> tuple[np.ndarray, CoilState]:
    """Plant input u produced by the unit for command cmd at air state x."""
    st = unit.operate(cmd[IDX_CMD_S], cmd[IDX_CMD_FAN], x[IDX_TA], x[IDX_W], d.T_o, ramp, M_c)
    u = np.empty(N_INPUTS)
    u[IDX_Q_SENS] = -st.Q_sens + st.P_fan_in      # indoor fan heat stays in the room
    u[IDX_M_LAT] = st.m_lat
    u[IDX_V_VENT] = cmd[IDX_CMD_VENT]
    return u, st


def simulate(
    p: RoomParams,
    scenario: Scenario,
    x0: np.ndarray,
    t_end: float,
    dt: float = 60.0,
    t_start: float = 0.0,
    method: str = "rk4",
    controller: Controller | None = None,
    dt_ctrl: float = 60.0,
    unit: ACUnit | None = None,
    unit_running_at_start: bool = False,
    M_coil0: float = 0.0,
) -> SimResult:
    step = STEPPERS[method]
    n_steps = _ratio(t_end - t_start, dt, "t_end - t_start")
    n_ctrl = _ratio(dt_ctrl, dt, "dt_ctrl") if controller is not None else 1

    n_u = 1 if unit is not None else 0          # unit states: water on the coil
    i_aux = N_STATES + n_u
    x = np.concatenate([np.asarray(x0, dtype=float), [M_coil0] * n_u, np.zeros(N_AUX)])
    out = np.zeros(N_CMD if unit is not None else N_INPUTS)    # held controller output
    running = unit_running_at_start
    t_started = -np.inf if running else t_start

    t_log = t_start + dt * np.arange(n_steps + 1)
    x_log = np.empty((n_steps + 1, i_aux + N_AUX))
    u_log = np.empty((n_steps + 1, N_INPUTS))
    c_log = np.zeros((n_steps + 1, N_CMD)) if unit is not None else None
    d_log = np.empty((n_steps + 1, len(Disturbances._fields)))
    q_log = np.empty((n_steps + 1, len(Fluxes._fields)))
    unit_log = np.zeros((n_steps + 1, len(UNIT_LOG))) if unit is not None else None

    def plant_input(tt: float, xx: np.ndarray, d: Disturbances):
        if unit is None:
            return out, OFF_STATE
        ramp = unit.start_ramp(tt - t_started) if out[IDX_CMD_S] > 0 else 1.0
        return unit_input(unit, out, ramp, xx, d, M_c=xx[N_STATES])

    def log(k: int, t: float, x: np.ndarray, d: Disturbances) -> None:
        u, st = plant_input(t, x, d)
        x_log[k] = x
        u_log[k] = u
        d_log[k] = [float(v) for v in d]
        q_log[k] = fluxes(t, x[:N_STATES], u, d, p)
        if unit is not None:
            c_log[k] = out
            unit_log[k] = [getattr(st, name) for name in UNIT_LOG]

    for k in range(n_steps):
        t = t_log[k]
        d_k = scenario(t)
        if controller is not None and k % n_ctrl == 0:
            y = measure(x[:N_STATES], p)
            if unit is not None:
                y = np.append(y, coil_sensor(plant_input(t, x, d_k)[1], x[IDX_TA]))
            out = np.asarray(controller.step(t, y), dtype=float)
            if unit is not None:
                now_on = out[IDX_CMD_S] > 0
                if now_on and not running:
                    t_started = t
                running = now_on

        held = {name: getattr(d_k, name) for name in PIECEWISE_CONSTANT}

        def f(tt: float, xx: np.ndarray) -> np.ndarray:
            d = scenario(tt)._replace(**held)
            u, st = plant_input(tt, xx, d)
            fl = fluxes(tt, xx[:N_STATES], u, d, p)
            return np.concatenate([derivatives(fl, p), [st.dM_coil] * n_u, [fl.boundary_total,
                                   fl.water_boundary_total, fl.co2_total, st.P_el,
                                   fl.water_external, st.m_drain]])

        log(k, t, x, d_k)
        x = step(f, t, x, dt)

    log(n_steps, t_log[-1], x, scenario(t_log[-1]))

    return SimResult(
        t=t_log,
        x=x_log[:, :N_STATES],
        u=u_log,
        d=Disturbances(*d_log.T),
        q={name: q_log[:, i] for i, name in enumerate(Fluxes._fields)},
        E_in=x_log[:, i_aux],
        M_in=x_log[:, i_aux + 1],
        CO2_in=x_log[:, i_aux + 2],
        E_el=x_log[:, i_aux + 3],
        W_ext=x_log[:, i_aux + 4],
        M_drain=x_log[:, i_aux + 5],
        M_coil=x_log[:, N_STATES] if unit is not None else None,
        cmd=c_log,
        unit={name: unit_log[:, i] for i, name in enumerate(UNIT_LOG)} if unit is not None else None,
        p=p,
        method=method,
        dt=dt,
    )


def periodic_state(
    p: RoomParams,
    scenario: Scenario,
    x_guess: np.ndarray,
    period: float = 86_400.0,
    dt: float = 60.0,
    method: str = "rk4",
    controller_factory=None,
    dt_ctrl: float = 60.0,
    unit: ACUnit | None = None,
    max_iter: int = 10,
) -> np.ndarray:
    """Initial state of the periodic steady state, found by Newton shooting.

    Solves F(x) = Phi_T(x) - x = 0, where Phi_T maps the state at t = 0 to the state one
    period later. The mass time constant (~3 days) makes plain warm-up slow: after a week
    the daily maximum still drifts by ~0.15 K/day. The thermal part is linear, so Newton
    converges in one or two iterations; the Jacobian of Phi_T is built by finite differences.
    ``controller_factory`` returns a fresh controller per run, so controller state is reset.
    """
    scale = STATE_SCALE                           # finite-difference steps
    tol = 1e-3 * scale                            # periodicity error accepted after one period
    x = np.asarray(x_guess, dtype=float).copy()

    def phi(x0: np.ndarray) -> np.ndarray:
        ctrl = controller_factory() if controller_factory else None
        r = simulate(p, scenario, x0, period, dt=dt, method=method, controller=ctrl,
                     dt_ctrl=dt_ctrl, unit=unit)
        return r.x[-1]

    for _ in range(max_iter):
        F = phi(x) - x
        if np.all(np.abs(F) <= tol):
            return x
        Jf = np.empty((N_STATES, N_STATES))
        for j in range(N_STATES):
            e = np.zeros(N_STATES)
            e[j] = scale[j]
            Jf[:, j] = (phi(x + e) - phi(x - e)) / (2 * scale[j])
        x = x - np.linalg.solve(Jf - np.eye(N_STATES), F)
    raise RuntimeError("periodic shooting did not converge")
