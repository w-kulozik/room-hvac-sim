"""Analytical reference solutions and convergence studies.

Shared by the test suite and by the scripts that draw the verification figures, so the
numbers in the README are produced by exactly the code that the tests check.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from . import psychro
from .inputs import DAY, ConstantInputs, SummerScenario
from .params import RoomParams
from .plant import IDX_PHI_B, IDX_W, U_OFF, buffer_flux, equilibrium_state
from .simulate import simulate


# ---- steady state ------------------------------------------------------------------------

def steady_state(p: RoomParams, d: ConstantInputs) -> np.ndarray:
    """Analytical steady state [T_a, T_m, W, phi_b] under constant inputs and u = 0."""
    q_sol = p.g_win * p.A_win * d.G_sol
    q_occ = p.q_person_sens * d.n_occ
    A = np.array([[-(p.H_am + p.H_ao + p.H_inf), p.H_am],
                  [p.H_am, -(p.H_am + p.H_mo)]])
    b = -np.array([(p.H_ao + p.H_inf) * d.T_o + p.f_sol_air * q_sol + p.f_person_conv * q_occ,
                   p.H_mo * d.T_o + (1 - p.f_sol_air) * q_sol + (1 - p.f_person_conv) * q_occ])
    T_a, T_m = np.linalg.solve(A, b)
    W = d.W_o + d.n_occ * p.g_vap_person / p.m_dot_inf
    return equilibrium_state(T_a, T_m, W)


# ---- one-node limit ----------------------------------------------------------------------

def one_node_params(p: RoomParams = RoomParams()) -> RoomParams:
    """Air node decoupled from the mass and moisture buffer: a single first-order system."""
    return replace(p, h_i=0.0, H_mo=0.0, moisture_buffer=False)


def one_node_error(method: str, dt: float, t_end: float = 48_000.0,
                   T_0: float = 20.0, T_o: float = 30.0) -> float:
    """|T_a,num - T_a,exact| at t_end after a step in outdoor temperature from T_0 to T_o."""
    p = one_node_params()
    tau = p.C_a / (p.H_ao + p.H_inf)
    r = simulate(p, ConstantInputs(T_o=T_o, W_o=0.01), equilibrium_state(T_0, T_0, 0.01),
                 t_end, dt=dt, method=method)
    exact = T_o + (T_0 - T_o) * np.exp(-t_end / tau)
    return abs(r.T_a[-1] - exact)


# ---- energy-balance closure with an independent quadrature --------------------------------

def simpson(y: np.ndarray, dt: float) -> float:
    """Composite Simpson rule on a uniform grid with an even number of intervals."""
    n = len(y) - 1
    if n % 2:
        raise ValueError("Simpson's rule needs an even number of intervals")
    return dt / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum())


def balance_params(p: RoomParams = RoomParams()) -> RoomParams:
    # Buffer off: it does not affect the heat balance but would make explicit Euler unstable
    # above ~790 s, and the study needs coarse steps.
    return replace(p, moisture_buffer=False)


def smooth_scenario() -> SummerScenario:
    # Only the sinusoidal outdoor temperature: occupancy steps and the kinks of the
    # half-sine sun would cap the convergence order of any fixed-grid quadrature.
    return SummerScenario(G_max=0.0, n_people=0.0)


def energy_residual_simpson(method: str, dt: float, t_end: float = DAY) -> float:
    """|stored-energy change - Simpson integral of the logged boundary heat flows|, J."""
    p = balance_params()
    r = simulate(p, smooth_scenario(), equilibrium_state(20.0, 24.0, 0.01), t_end,
                 dt=dt, method=method)
    return abs(r.energy_stored[-1] - simpson(r.boundary_heat_flow(), dt))


def convergence_order(dts, errors) -> float:
    """Least-squares slope of log(error) against log(dt)."""
    return float(np.polyfit(np.log(dts), np.log(errors), 1)[0])


# ---- moisture decay ----------------------------------------------------------------------

def moisture_decay_exact_no_buffer(t, W_0: float, W_o: float, p: RoomParams):
    return W_o + (W_0 - W_o) * np.exp(-t * p.ach / 3600.0)


def moisture_linear_system(p: RoomParams, T: float, W_star: float):
    """Linearised isothermal air + buffer moisture system around equilibrium at W_star.

    Returns (A, x_star) for d/dt [W, phi_b] = A ([W, phi_b] - x_star).
    """
    phi_star = float(psychro.RH_from_W(T, W_star))
    dpw_dW = psychro.P_ATM * psychro.EPS / (psychro.EPS + W_star) ** 2
    G, pws = p.buffer_conductance, float(psychro.p_ws(T))
    A = np.array([[(-p.m_dot_inf - G * dpw_dW) / p.m_air, G * pws / p.m_air],
                  [G * dpw_dW / p.buffer_capacity, -G * pws / p.buffer_capacity]])
    return A, np.array([W_star, phi_star])


def moisture_decay_exact_buffer(t, x0, A, x_star):
    """Sum of two exponentials: x(t) = x* + V exp(L t) V^-1 (x0 - x*)."""
    lam, V = np.linalg.eig(A)
    c = np.linalg.solve(V, np.asarray(x0) - x_star)
    t = np.atleast_1d(t)
    return x_star + (V @ (c[:, None] * np.exp(np.outer(lam, t)))).T.real
