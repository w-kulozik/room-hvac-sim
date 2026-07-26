"""Analytical reference solutions and convergence studies.

Shared by the test suite and by the scripts that draw the verification figures, so the
numbers in docs/verification.md are produced by exactly the code that the tests check.
"""

from __future__ import annotations

from dataclasses import replace
from typing import cast

import numpy as np
import numpy.typing as npt

from . import psychro
from ._typing import FloatArray
from .inputs import DAY, ConstantInputs, SummerScenario
from .params import RoomParams
from .plant import equilibrium_state
from .simulate import simulate

_DEFAULT_ROOM = RoomParams()

# ---- steady state ------------------------------------------------------------------------

def steady_state(p: RoomParams, d: ConstantInputs) -> FloatArray:
    """Return the analytical steady state under constant inputs and u = 0.

    Parameters
    ----------
    p : RoomParams
        Room parameters.
    d : ConstantInputs
        Constant boundary conditions (outdoor air, sun, occupancy).

    Returns
    -------
    FloatArray
        State [T_a (degC), T_m (degC), W (kg/kg), phi_b (-), C (ppm)], with the buffer in
        equilibrium with the air.
    """
    q_sol = p.g_win * p.A_win * d.G_sol
    q_occ = p.q_person_sens * d.n_occ
    A = np.array([[-(p.H_am + p.H_ao + p.H_inf), p.H_am],
                  [p.H_am, -(p.H_am + p.H_mo)]])
    b = -np.array([(p.H_ao + p.H_inf) * d.T_o + p.f_sol_air * q_sol + p.f_person_conv * q_occ,
                   p.H_mo * d.T_o + (1 - p.f_sol_air) * q_sol + (1 - p.f_person_conv) * q_occ])
    T_a, T_m = np.linalg.solve(A, b)
    W = d.W_o + d.n_occ * p.g_vap_person / p.m_dot_inf
    C = d.C_o + d.n_occ * p.co2_person / p.V_dot_inf * 1e6
    return equilibrium_state(T_a, T_m, W, C)


# ---- one-node limit ----------------------------------------------------------------------

def one_node_params(p: RoomParams = _DEFAULT_ROOM) -> RoomParams:
    """Return parameters with the air node decoupled from the mass and moisture buffer.

    The air node then behaves as a single first-order system.

    Parameters
    ----------
    p : RoomParams, optional
        Base room parameters.

    Returns
    -------
    RoomParams
        Copy of ``p`` with h_i = 0, H_mo = 0 and the moisture buffer off.
    """
    return replace(p, h_i=0.0, H_mo=0.0, moisture_buffer=False)


def one_node_error(method: str, dt: float, t_end: float = 48_000.0,
                   T_0: float = 20.0, T_o: float = 30.0) -> float:
    """Return the error of the one-node step response at t_end.

    The room starts in equilibrium at T_0 and the outdoor temperature is T_o from t = 0.

    Parameters
    ----------
    method : str
        Integration method passed to ``simulate``.
    dt : float
        Time step, s.
    t_end : float, optional
        End time, s.
    T_0 : float, optional
        Initial room temperature, degC.
    T_o : float, optional
        Outdoor temperature after the step, degC.

    Returns
    -------
    float
        |T_a,num - T_a,exact| at t_end, K.
    """
    p = one_node_params()
    tau = p.C_a / (p.H_ao + p.H_inf)
    r = simulate(p, ConstantInputs(T_o=T_o, W_o=0.01), equilibrium_state(T_0, T_0, 0.01),
                 t_end, dt=dt, method=method)
    exact = T_o + (T_0 - T_o) * np.exp(-t_end / tau)
    return cast(float, abs(r.T_a[-1] - exact))


# ---- energy-balance closure with an independent quadrature --------------------------------

def simpson(y: FloatArray, dt: float) -> float:
    """Integrate samples with the composite Simpson rule on a uniform grid.

    Parameters
    ----------
    y : FloatArray
        Samples at uniform spacing; len(y) - 1 must be even.
    dt : float
        Sample spacing, s.

    Returns
    -------
    float
        Integral of y over the grid, in units of y times s.

    Raises
    ------
    ValueError
        If the number of intervals is odd.
    """
    n = len(y) - 1
    if n % 2:
        raise ValueError("Simpson's rule needs an even number of intervals")
    return cast(float,
                dt / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum()))


def balance_params(p: RoomParams = _DEFAULT_ROOM) -> RoomParams:
    """Return parameters for the energy-balance study, with the moisture buffer off.

    The buffer does not affect the heat balance but would make explicit Euler unstable above
    about 790 s, and the study needs coarse steps.

    Parameters
    ----------
    p : RoomParams, optional
        Base room parameters.

    Returns
    -------
    RoomParams
        Copy of ``p`` with the moisture buffer off.
    """
    return replace(p, moisture_buffer=False)


def smooth_scenario() -> SummerScenario:
    """Return a summer scenario with only the sinusoidal outdoor temperature.

    Occupancy steps and the kinks of the half-sine sun would cap the convergence order of
    any fixed-grid quadrature, so both are switched off.

    Returns
    -------
    SummerScenario
        Scenario with G_max = 0 and no occupants.
    """
    return SummerScenario(G_max=0.0, n_people=0.0)


def energy_residual_simpson(method: str, dt: float, t_end: float = DAY) -> float:
    """Return the energy-balance residual against an independent Simpson quadrature.

    Parameters
    ----------
    method : str
        Integration method passed to ``simulate``.
    dt : float
        Time step, s.
    t_end : float, optional
        End time, s.

    Returns
    -------
    float
        |stored-energy change - Simpson integral of the logged boundary heat flows|, J.
    """
    p = balance_params()
    r = simulate(p, smooth_scenario(), equilibrium_state(20.0, 24.0, 0.01), t_end,
                 dt=dt, method=method)
    return cast(float, abs(r.energy_stored[-1] - simpson(r.boundary_heat_flow(), dt)))


def convergence_order(dts: npt.ArrayLike, errors: npt.ArrayLike) -> float:
    """Return the least-squares slope of log(error) against log(dt).

    Parameters
    ----------
    dts : array_like
        Time steps, s.
    errors : array_like
        Errors at those time steps, any consistent unit.

    Returns
    -------
    float
        Observed order of convergence, dimensionless.
    """
    return float(np.polyfit(np.log(dts), np.log(errors), 1)[0])


# ---- moisture decay ----------------------------------------------------------------------

def moisture_decay_exact_no_buffer(t: FloatArray, W_0: float, W_o: float,
                                   p: RoomParams) -> FloatArray:
    """Return the exact humidity decay of the room air without the moisture buffer.

    Parameters
    ----------
    t : FloatArray
        Times since the start of the decay, s.
    W_0 : float
        Initial room humidity ratio, kg/kg.
    W_o : float
        Outdoor humidity ratio, kg/kg.
    p : RoomParams
        Room parameters (only the air change rate is used).

    Returns
    -------
    FloatArray
        Room humidity ratio at ``t``, kg/kg.
    """
    return W_o + (W_0 - W_o) * np.exp(-t * p.ach / 3600.0)


def moisture_linear_system(p: RoomParams, T: float,
                           W_star: float) -> tuple[FloatArray, FloatArray]:
    """Linearise the isothermal air + buffer moisture system around equilibrium at W_star.

    Parameters
    ----------
    p : RoomParams
        Room parameters.
    T : float
        Air and buffer temperature, degC.
    W_star : float
        Equilibrium humidity ratio, kg/kg.

    Returns
    -------
    A : FloatArray
        2x2 system matrix, 1/s, for d/dt [W, phi_b] = A ([W, phi_b] - x_star).
    x_star : FloatArray
        Equilibrium state [W (kg/kg), phi_b (-)].
    """
    phi_star = float(psychro.RH_from_W(T, W_star))
    dpw_dW = psychro.P_ATM * psychro.EPS / (psychro.EPS + W_star) ** 2
    G, pws = p.buffer_conductance, float(psychro.p_ws(T))
    A = np.array([[(-p.m_dot_inf - G * dpw_dW) / p.m_air, G * pws / p.m_air],
                  [G * dpw_dW / p.buffer_capacity, -G * pws / p.buffer_capacity]])
    return A, np.array([W_star, phi_star])


def moisture_decay_exact_buffer(t: npt.ArrayLike, x0: npt.ArrayLike, A: FloatArray,
                                x_star: FloatArray) -> FloatArray:
    """Return the exact solution of the linear moisture system, a sum of two exponentials.

    x(t) = x* + V exp(L t) V^-1 (x0 - x*), with A = V L V^-1.

    Parameters
    ----------
    t : array_like
        Time or times, s.
    x0 : array_like
        Initial state [W (kg/kg), phi_b (-)].
    A : FloatArray
        System matrix from ``moisture_linear_system``, 1/s.
    x_star : FloatArray
        Equilibrium state [W (kg/kg), phi_b (-)].

    Returns
    -------
    FloatArray
        States, shape (len(t), 2): columns W (kg/kg) and phi_b (-).
    """
    lam, V = np.linalg.eig(A)
    c = np.linalg.solve(V, np.asarray(x0) - x_star)
    t = np.atleast_1d(t)
    return cast(FloatArray, x_star + (V @ (c[:, None] * np.exp(np.outer(lam, t)))).T.real)
