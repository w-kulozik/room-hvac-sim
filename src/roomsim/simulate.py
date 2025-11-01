"""Fixed-step simulation loop: continuous plant, discrete controller, logging.

Per step [t_k, t_k + dt):
  1. if t_k is a controller sampling instant, read the sensors and update u (then hold it);
  2. sample the piecewise-constant disturbances at t_k;
  3. advance the plant with Euler or RK4; continuous disturbances are evaluated at the stage times.

Cumulative boundary heat, water and CO2 are integrated as three extra states with the same method,
so the balance closure in ``SimResult`` reflects exactly what the integrator did.
"""

from __future__ import annotations

from typing import Protocol

import numpy as np

from .inputs import PIECEWISE_CONSTANT, Disturbances, Scenario
from .integrators import STEPPERS
from .params import RoomParams
from .plant import N_INPUTS, N_STATES, STATE_SCALE, Fluxes, derivatives, fluxes, measure
from .results import SimResult


class Controller(Protocol):
    def step(self, t: float, y: np.ndarray) -> np.ndarray:
        """Return u = [Q_hvac_sens, m_hvac_lat, V_vent] from the measurement y = [T_a, RH_a, C]."""
        ...


N_AUX = 3   # cumulative boundary heat (J), water (kg), CO2 (m3)


def _ratio(a: float, b: float, what: str) -> int:
    n = round(a / b)
    if n < 1 or abs(n * b - a) > 1e-9 * max(abs(a), 1.0):
        raise ValueError(f"{what} must be a positive integer multiple of the step")
    return n


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
) -> SimResult:
    step = STEPPERS[method]
    n_steps = _ratio(t_end - t_start, dt, "t_end - t_start")
    n_ctrl = _ratio(dt_ctrl, dt, "dt_ctrl") if controller is not None else 1

    x = np.concatenate([np.asarray(x0, dtype=float), np.zeros(N_AUX)])   # + E_in, M_in, CO2_in
    u = np.zeros(N_INPUTS)

    t_log = t_start + dt * np.arange(n_steps + 1)
    x_log = np.empty((n_steps + 1, N_STATES + N_AUX))
    u_log = np.empty((n_steps + 1, N_INPUTS))
    d_log = np.empty((n_steps + 1, len(Disturbances._fields)))
    q_log = np.empty((n_steps + 1, len(Fluxes._fields)))

    def log(k: int, t: float, x: np.ndarray, d: Disturbances) -> None:
        x_log[k] = x
        u_log[k] = u
        d_log[k] = [float(v) for v in d]
        q_log[k] = fluxes(t, x[:N_STATES], u, d, p)

    for k in range(n_steps):
        t = t_log[k]
        if controller is not None and k % n_ctrl == 0:
            u = np.asarray(controller.step(t, measure(x[:N_STATES], p)), dtype=float)

        d_k = scenario(t)
        held = {name: getattr(d_k, name) for name in PIECEWISE_CONSTANT}

        def f(tt: float, xx: np.ndarray) -> np.ndarray:
            d = scenario(tt)._replace(**held)
            fl = fluxes(tt, xx[:N_STATES], u, d, p)
            return np.concatenate([derivatives(fl, p),
                                   [fl.boundary_total, fl.water_boundary_total, fl.co2_total]])

        log(k, t, x, d_k)
        x = step(f, t, x, dt)

    log(n_steps, t_log[-1], x, scenario(t_log[-1]))

    return SimResult(
        t=t_log,
        x=x_log[:, :N_STATES],
        u=u_log,
        d=Disturbances(*d_log.T),
        q={name: q_log[:, i] for i, name in enumerate(Fluxes._fields)},
        E_in=x_log[:, N_STATES],
        M_in=x_log[:, N_STATES + 1],
        CO2_in=x_log[:, N_STATES + 2],
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
        r = simulate(p, scenario, x0, period, dt=dt, method=method, controller=ctrl, dt_ctrl=dt_ctrl)
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
