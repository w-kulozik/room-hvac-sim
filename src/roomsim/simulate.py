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

Integrated state vector: [room state (5), water on the coil (with a unit), cumulative integrals].
The cumulative integrals (boundary heat, room water, CO2, external water, electrical energy,
drained condensate, net cooling delivered) use the same method as the state, so the balances in
``SimResult`` reflect exactly what the integrator did.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Protocol, cast

import numpy as np

from ._typing import FloatArray
from .ac_unit import OFF_STATE, ACUnit, CoilState, net_cooling
from .inputs import PIECEWISE_CONSTANT, Disturbances, Scenario
from .integrators import STEPPERS, euler_stability_limit, rk4_stability_limit
from .params import RoomParams
from .plant import (
    IDX_CO2,
    IDX_M_LAT,
    IDX_PHI_B,
    IDX_Q_SENS,
    IDX_TA,
    IDX_V_VENT,
    IDX_W,
    N_INPUTS,
    N_STATES,
    STATE_SCALE,
    Fluxes,
    derivatives,
    fluxes,
    jacobian,
    measure,
)
from .results import UNIT_FIELDS, SimResult, UnitLog

# Command vector when a unit is present
IDX_CMD_S, IDX_CMD_FAN, IDX_CMD_VENT = 0, 1, 2
N_CMD = 3

# Cumulative integrals appended to the state
AUX = ("E_in", "M_in", "CO2_in", "W_ext", "E_el", "M_drain", "E_cool")

STATE_NAMES = ("T_a", "T_m", "W", "phi_b", "C")
STABILITY_LIMITS = {"euler": euler_stability_limit, "rk4": rk4_stability_limit}


class Controller(Protocol):
    """Discrete-time controller called by the simulator at each sampling instant."""

    def step(self, t: float, y: FloatArray) -> FloatArray:
        """Compute the controller output from the measurement.

        Parameters
        ----------
        t : float
            Sampling time, s.
        y : ndarray
            Measurement [T_a degC, RH_a fraction, C ppm], plus the coil-sensor temperature
            (degC) when a unit is present.

        Returns
        -------
        ndarray, shape (3,)
            Without a unit: u = [Q_hvac_sens W, m_hvac_lat kg/s, V_vent m3/s].
            With a unit:    command = [s -, f -, V_vent m3/s].
        """
        ...


def _ratio(a: float, b: float, what: str) -> int:
    n = round(a / b)
    if n < 1 or abs(n * b - a) > 1e-9 * max(abs(a), 1.0):
        raise ValueError(f"{what} must be a positive integer multiple of the step")
    return n


def coil_sensor(st: CoilState, T_a: float) -> float:
    """Return the evaporator thermistor reading.

    Parameters
    ----------
    st : CoilState
        Unit operating point.
    T_a : float
        Indoor air temperature, degC.

    Returns
    -------
    float
        Coil surface temperature while the unit runs (coil load or indoor fan) and T_s is not
        NaN, otherwise the room air temperature, degC.
    """
    running = st.Q_coil > 0.0 or st.P_fan_in > 0.0
    return st.T_s if running and not np.isnan(st.T_s) else T_a


def unit_input(unit: ACUnit, cmd: FloatArray, ramp: float, x: FloatArray,
               d: Disturbances, M_c: float | None = None) -> tuple[FloatArray, CoilState]:
    """Compute the plant input u produced by the unit for command cmd at air state x.

    Parameters
    ----------
    unit : ACUnit
        Air-conditioning unit.
    cmd : ndarray, shape (3,)
        Command [s -, f -, V_vent m3/s].
    ramp : float
        Capacity fraction after a compressor start, -.
    x : ndarray
        State whose first entries are the room state [T_a degC, T_m degC, W kg/kg, ...].
    d : Disturbances
        Disturbances sampled at a single time.
    M_c : float, optional
        Water held on the coil, kg; passed to ``ACUnit.operate``.

    Returns
    -------
    u : ndarray, shape (3,)
        Plant input [Q_hvac_sens W, m_hvac_lat kg/s, V_vent m3/s]; Q_hvac_sens includes the
        indoor-fan heat.
    st : CoilState
        Unit operating point.
    """
    T_o = cast(float, d.T_o)    # d is sampled at a single time, so T_o is a scalar
    st = unit.operate(cmd[IDX_CMD_S], cmd[IDX_CMD_FAN], x[IDX_TA], x[IDX_W], T_o, ramp, M_c)
    u = np.empty(N_INPUTS)
    u[IDX_Q_SENS] = -st.Q_sens + st.P_fan_in      # indoor fan heat stays in the room
    u[IDX_M_LAT] = st.m_lat
    u[IDX_V_VENT] = cmd[IDX_CMD_VENT]
    return u, st


class _Actuation:
    """Turns the held controller output into plant inputs; tracks compressor starts."""

    def __init__(self, unit: ACUnit | None, t_start: float, running_at_start: bool):
        self.unit = unit
        self.out = np.zeros(N_CMD if unit is not None else N_INPUTS)
        self.running = running_at_start
        self.t_started = -np.inf if running_at_start else t_start

    @property
    def n_states(self) -> int:
        """Number of actuator states appended to the room state, -."""
        return 1 if self.unit is not None else 0       # water held on the coil

    def set_output(self, t: float, out: FloatArray) -> None:
        """Hold a new controller output from time t (s) and record compressor starts.

        Raises ValueError if the ventilation flow V_vent (m3/s) is negative or NaN: a balanced
        supply fan cannot deliver a negative flow, and a negative one would reverse the
        ventilation heat, moisture and CO2 terms.
        """
        self.out = np.asarray(out, dtype=float)
        V_vent = self.out[IDX_CMD_VENT if self.unit is not None else IDX_V_VENT]
        if not V_vent >= 0.0:
            raise ValueError(f"controller output V_vent = {V_vent!r} m3/s at t = {t!r} s is "
                             "invalid; it must be >= 0")
        if self.unit is not None:
            on = self.out[IDX_CMD_S] > 0
            if on and not self.running:
                self.t_started = t
            self.running = on

    def inputs(self, t: float, x: FloatArray, d: Disturbances) -> tuple[FloatArray, CoilState]:
        """Return the plant input and unit operating point at time t (s) and integrated state x."""
        if self.unit is None:
            return self.out, OFF_STATE
        ramp = self.unit.start_ramp(t - self.t_started) if self.out[IDX_CMD_S] > 0 else 1.0
        return unit_input(self.unit, self.out, ramp, x, d, M_c=x[N_STATES])


def simulate(
    p: RoomParams,
    scenario: Scenario,
    x0: FloatArray,
    duration: float,
    dt: float = 60.0,
    t_start: float = 0.0,
    method: str = "rk4",
    controller: Controller | None = None,
    dt_ctrl: float = 60.0,
    unit: ACUnit | None = None,
    unit_running_at_start: bool = False,
    M_coil0: float = 0.0,
) -> SimResult:
    """Simulate ``duration`` seconds from time ``t_start`` (seconds from midnight of day 0).

    Parameters
    ----------
    p : RoomParams
        Room parameters.
    scenario : Scenario
        Disturbances as a function of time.
    x0 : ndarray, shape (5,)
        Initial room state [T_a degC, T_m degC, W kg/kg, phi_b -, C ppm].
    duration : float
        Simulated time, s; a positive integer multiple of dt.
    dt : float, optional
        Integration step, s.
    t_start : float, optional
        Start time from midnight of day 0, s.
    method : str, optional
        Integrator, a key of ``STEPPERS`` ("euler" or "rk4").
    controller : Controller, optional
        Discrete controller; without one the plant input is held at zero.
    dt_ctrl : float, optional
        Controller sampling interval, s; a positive integer multiple of dt.
    unit : ACUnit, optional
        Air-conditioning unit; if given, the controller output is a unit command.
    unit_running_at_start : bool, optional
        Treat the compressor as already running (no start ramp) at t_start.
    M_coil0 : float, optional
        Initial water held on the coil, kg (used only with a unit).

    Returns
    -------
    SimResult
        Logged time series on the step grid and the cumulative balance integrals.

    Raises
    ------
    ValueError
        If duration or dt_ctrl is not a positive integer multiple of dt, or if the controller
        returns a negative ventilation flow.
    RuntimeError
        If the state becomes non-finite or W, C or (with an active buffer) phi_b becomes
        negative, which in practice means dt exceeds the stability limit of the integrator.
    KeyError
        If method is not a key of ``STEPPERS``.
    """
    step = STEPPERS[method]
    n_steps = _ratio(duration, dt, "duration")
    n_ctrl = _ratio(dt_ctrl, dt, "dt_ctrl") if controller is not None else 1
    # a controller whose discrete elements assume a sampling period (attribute Ts) must be
    # called at that period: at another one its integrators and predictions are scaled wrongly
    Ts_ctrl = getattr(controller, "Ts", None)
    if Ts_ctrl is not None and abs(Ts_ctrl - dt_ctrl) > 1e-9 * dt_ctrl:
        raise ValueError(f"dt_ctrl = {dt_ctrl!r} s does not match the controller's sampling "
                         f"period Ts = {Ts_ctrl!r} s; build the controller with Ts = dt_ctrl")
    act = _Actuation(unit, t_start, unit_running_at_start)
    i_aux = N_STATES + act.n_states
    x = np.concatenate([np.asarray(x0, dtype=float), [M_coil0] * act.n_states, np.zeros(len(AUX))])

    t_log = cast(FloatArray, t_start + dt * np.arange(n_steps + 1))   # float64 by promotion
    x_log = np.empty((n_steps + 1, len(x)))
    u_log = np.empty((n_steps + 1, N_INPUTS))
    d_log = np.empty((n_steps + 1, len(Disturbances._fields)))
    q_log = np.empty((n_steps + 1, len(Fluxes._fields)))
    c_log = np.zeros((n_steps + 1, N_CMD))
    op_log = np.zeros((n_steps + 1, len(UNIT_FIELDS)))

    def log(k: int, t: float, x: FloatArray, d: Disturbances) -> None:
        u, st = act.inputs(t, x, d)
        x_log[k], u_log[k] = x, u
        d_log[k] = [float(v) for v in d]
        q_log[k] = fluxes(t, x[:N_STATES], u, d, p)
        if unit is not None:
            c_log[k] = act.out
            op_log[k] = [getattr(st, name) for name in UNIT_FIELDS]

    buffer_active = p.buffer_capacity > 0.0     # otherwise phi_b is inert and may take any value

    for k in range(n_steps):
        t = t_log[k]
        d_k = scenario(t)
        if controller is not None and k % n_ctrl == 0:
            y = measure(x[:N_STATES], p)
            if unit is not None:
                y = np.append(y, coil_sensor(act.inputs(t, x, d_k)[1], x[IDX_TA]))
            act.set_output(t, controller.step(t, y))
        held = {name: getattr(d_k, name) for name in PIECEWISE_CONSTANT}

        def f(tt: float, xx: FloatArray) -> FloatArray:
            # f is called only by the step below, within this iteration
            d = scenario(tt)._replace(**held)  # noqa: B023
            u, st = act.inputs(tt, xx, d)
            fl = fluxes(tt, xx[:N_STATES], u, d, p)
            aux = [fl.boundary_total, fl.water_boundary_total, fl.co2_total, fl.water_external,
                   st.P_el, st.m_drain, net_cooling(st, xx[IDX_TA])]
            return np.concatenate([derivatives(fl, p), [st.dM_coil] * act.n_states, aux])

        log(k, t, x, d_k)
        x = step(f, t, x, dt)
        if not (np.all(np.isfinite(x)) and x[IDX_W] >= 0.0 and x[IDX_CO2] >= 0.0
                and (x[IDX_PHI_B] >= 0.0 or not buffer_active)):
            raise RuntimeError(_divergence_message(p, scenario, method, dt, t_start,
                                                   t_log[k + 1], x, x_log[0], u_log[0]))

    log(n_steps, t_log[-1], x, scenario(t_log[-1]))

    aux = {name: x_log[:, i_aux + i] for i, name in enumerate(AUX)}
    unit_log = None
    if unit is not None:
        unit_log = UnitLog(cmd=c_log, op={n: op_log[:, i] for i, n in enumerate(UNIT_FIELDS)},
                           M_coil=x_log[:, N_STATES], E_el=aux["E_el"], M_drain=aux["M_drain"],
                           E_cool=aux["E_cool"])
    return SimResult(
        t=t_log, x=x_log[:, :N_STATES], u=u_log, d=Disturbances(*d_log.T),
        q={name: q_log[:, i] for i, name in enumerate(Fluxes._fields)},
        E_in=aux["E_in"], M_in=aux["M_in"], CO2_in=aux["CO2_in"], W_ext=aux["W_ext"],
        p=p, method=method, dt=dt, unit=unit_log,
    )


def _divergence_message(p: RoomParams, scenario: Scenario, method: str, dt: float,
                        t_start: float, t: float, x: FloatArray, x0: FloatArray,
                        u0: FloatArray) -> str:
    """Describe a state that became non-finite or unphysical, with the likely cause.

    The stability limit is estimated from the plant Jacobian at the initial state and input
    (without the feedback of a controller or unit), so it is indicative only.
    """
    state = ", ".join(f"{n} = {v:.4g}" for n, v in zip(STATE_NAMES, x[:N_STATES], strict=True))
    msg = (f"simulation diverged at t = {t:.0f} s: the state became non-finite or unphysical "
           f"(negative W, C or active-buffer phi_b): {state}. The likely cause is a step "
           f"dt = {dt:g} s above the stability limit of the explicit '{method}' integrator")
    try:
        J = jacobian(t_start, x0[:N_STATES], u0, scenario(t_start), p)
        limit = STABILITY_LIMITS[method](J)
    except (KeyError, ValueError, np.linalg.LinAlgError):
        return msg + "; reduce dt."
    return msg + f" (about {limit:.0f} s at the initial state); reduce dt."


def periodic_state(
    p: RoomParams,
    scenario: Scenario,
    x_guess: FloatArray,
    period: float = 86_400.0,
    dt: float = 60.0,
    method: str = "rk4",
    controller_factory: Callable[[], Controller] | None = None,
    dt_ctrl: float = 60.0,
    unit: ACUnit | None = None,
    max_iter: int = 10,
) -> FloatArray:
    """Find the initial room state of the periodic steady state by Newton shooting.

    Solves F(x) = Phi_T(x) - x = 0, where Phi_T maps the room state at t = 0 to the state one
    period later. The mass time constant (~3 days) makes plain warm-up slow: after a week
    the daily maximum still drifts by ~0.15 K/day. The thermal part is linear, so Newton
    converges in one or two iterations; the Jacobian of Phi_T is built by finite differences.
    ``controller_factory`` returns a fresh controller per run. Controller states and the water
    on the coil restart at each evaluation, so with a controller the result is periodic only in
    the room state; closed-loop runs use warm-up days instead.

    Parameters
    ----------
    p : RoomParams
        Room parameters.
    scenario : Scenario
        Disturbances; periodic with ``period`` for a meaningful result.
    x_guess : ndarray, shape (5,)
        Initial guess of the room state [T_a degC, T_m degC, W kg/kg, phi_b -, C ppm].
    period : float, optional
        Period, s.
    dt : float, optional
        Integration step, s.
    method : str, optional
        Integrator, "euler" or "rk4".
    controller_factory : callable, optional
        Returns a new controller for each simulated period.
    dt_ctrl : float, optional
        Controller sampling interval, s.
    unit : ACUnit, optional
        Air-conditioning unit.
    max_iter : int, optional
        Maximum number of Newton iterations.

    Returns
    -------
    ndarray, shape (5,)
        Room state at t = 0 whose image after one period differs from it by at most
        ``1e-3 * STATE_SCALE`` in every component.

    Raises
    ------
    RuntimeError
        If the shooting does not converge within ``max_iter`` iterations, or if a state with
        no dynamics changes over one period, so that no periodic state exists.

    Notes
    -----
    A state whose value after one period does not depend on any state (phi_b with the
    buffer off, T_m with h_i = H_mo = 0) has no dynamics: its row of Jf - I is zero and
    the Newton matrix is singular. Such states keep their guess and Newton solves for the
    others; if one of them changes over the period, no periodic state exists.
    """
    scale = STATE_SCALE                           # finite-difference steps
    tol = 1e-3 * scale                            # periodicity error accepted after one period
    x = np.asarray(x_guess, dtype=float).copy()

    def phi(x0: FloatArray) -> FloatArray:
        ctrl = controller_factory() if controller_factory else None
        r = simulate(p, scenario, x0, period, dt=dt, method=method, controller=ctrl,
                     dt_ctrl=dt_ctrl, unit=unit)
        x_end: FloatArray = r.x[-1]
        return x_end

    for _ in range(max_iter):
        F = phi(x) - x
        if np.all(np.abs(F) <= tol):
            return x
        Jf = np.empty((N_STATES, N_STATES))
        for j in range(N_STATES):
            e = np.zeros(N_STATES)
            e[j] = scale[j]
            Jf[:, j] = (phi(x + e) - phi(x - e)) / (2 * scale[j])
        A = Jf - np.eye(N_STATES)
        # rows of A in units of STATE_SCALE; a zero row means the state has no dynamics
        frozen = np.all(np.abs(A) * scale[None, :] <= 1e-9 * scale[:, None], axis=1)
        drifting = frozen & (np.abs(F) > tol)
        if np.any(drifting):
            names = ", ".join(f"{STATE_NAMES[j]} (change {F[j]:.4g} per period)"
                              for j in np.flatnonzero(drifting))
            raise RuntimeError(f"no periodic state exists: {names} has no dynamics (its value "
                               "after one period does not depend on the state) but changes "
                               "over the period")
        if not np.any(frozen):
            x = x - np.linalg.solve(A, F)
        else:
            active = ~frozen
            x[active] -= np.linalg.solve(A[np.ix_(active, active)], F[active])
    raise RuntimeError("periodic shooting did not converge")
