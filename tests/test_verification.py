"""Verification against analytical solutions (specification, section 5)."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from roomsim import verification as v
from roomsim.inputs import DAY, HOUR, ConstantInputs, SummerScenario
from roomsim.params import RoomParams
from roomsim.plant import IDX_CO2, IDX_PHI_B, IDX_TA, IDX_TM, IDX_W, U_OFF, equilibrium_state, rhs
from roomsim.simulate import periodic_state, simulate

P = RoomParams()


# ---- 5.2 steady state --------------------------------------------------------------------

STEADY_CASES = [
    ConstantInputs(T_o=30.0, W_o=0.012),
    ConstantInputs(T_o=30.0, W_o=0.012, G_sol=300.0, n_occ=2),
    ConstantInputs(T_o=18.0, W_o=0.007, G_sol=0.0, n_occ=2),
]


@pytest.mark.parametrize("d", STEADY_CASES)
def test_steady_state_is_a_fixed_point(d):
    x_ss = v.steady_state(P, d)
    dx = rhs(0.0, x_ss, U_OFF, d(0.0), P)
    assert np.abs(dx[[IDX_TA, IDX_TM]]).max() < 1e-12      # K/s
    assert np.abs(dx[IDX_W]) < 1e-15                       # 1/s
    assert np.abs(dx[IDX_PHI_B]) < 1e-12
    assert np.abs(dx[IDX_CO2]) < 1e-10                     # ppm/s


@pytest.mark.parametrize("d", STEADY_CASES)
def test_simulation_converges_to_steady_state(d):
    # 60 days = 19 mass time constants, so the transient has decayed below 1e-7 K.
    # dt = 300 s: the air-buffer moisture mode stiffens with temperature (p_ws(T_m) grows),
    # and at the 48 degC steady state of case d1 the RK4 stability limit drops to ~530 s.
    x0 = equilibrium_state(24.0, 24.0, 0.009)
    r = simulate(P, d, x0, 60 * DAY, dt=300.0)
    x_ss = v.steady_state(P, d)
    assert r.x[-1, [IDX_TA, IDX_TM]] == pytest.approx(x_ss[[IDX_TA, IDX_TM]], abs=1e-6)
    assert r.x[-1, IDX_W] == pytest.approx(x_ss[IDX_W], abs=1e-10)
    assert r.x[-1, IDX_PHI_B] == pytest.approx(x_ss[IDX_PHI_B], abs=1e-7)
    assert r.x[-1, IDX_CO2] == pytest.approx(x_ss[IDX_CO2], abs=1e-6)


# ---- 5.3 one-node limit and convergence order ----------------------------------------------

DTS = np.array([4000.0, 2000.0, 1000.0, 500.0, 250.0])


def test_one_node_limit_is_first_order_exponential():
    p = v.one_node_params()
    tau = p.C_a / (p.H_ao + p.H_inf)
    assert tau == pytest.approx(24_970, rel=1e-3)
    r = simulate(p, ConstantInputs(T_o=30.0, W_o=0.01), equilibrium_state(20, 20, 0.01),
                 3 * DAY, dt=60.0)
    exact = 30.0 - 10.0 * np.exp(-r.t / tau)
    assert np.abs(r.T_a - exact).max() < 1e-9
    assert np.all(r.T_m == 20.0)          # mass fully decoupled


@pytest.mark.parametrize("method, order", [("rk4", 4.0), ("euler", 1.0)])
def test_convergence_order_one_node(method, order):
    errors = [v.one_node_error(method, dt) for dt in DTS]
    assert v.convergence_order(DTS, errors) == pytest.approx(order, abs=0.15)


# ---- 5.4 energy-balance closure -----------------------------------------------------------

@pytest.mark.parametrize("method", ["rk4", "euler"])
def test_energy_and_water_close_to_round_off(method):
    """(a) Boundary flows integrated as extra states: Runge-Kutta methods keep linear invariants."""
    S = SummerScenario()
    r = simulate(P, S, equilibrium_state(26.0, 30.0, 0.010), 3 * DAY, dt=60.0, method=method)
    E_abs = P.C_a * np.abs(r.T_a).max() + P.C_m * np.abs(r.T_m).max()   # stored energy, J
    assert np.abs(r.energy_stored).max() > 1e7      # the check is not trivially zero
    assert np.abs(r.water_stored).max() > 0.1
    # round-off: a few thousand steps at ~1e-16 relative precision of the stored totals
    assert np.abs(r.energy_residual).max() < 1e-11 * E_abs
    assert np.abs(r.water_residual).max() < 1e-11 * (P.m_air * 0.02 + P.buffer_capacity)


BALANCE_DTS = np.array([1200.0, 600.0, 300.0, 150.0])


@pytest.mark.parametrize("method, order", [("rk4", 4.0), ("euler", 1.0)])
def test_energy_residual_order_with_simpson(method, order):
    """(b) Independent quadrature of the logged flows: residual falls with the method order."""
    res = [v.energy_residual_simpson(method, dt) for dt in BALANCE_DTS]
    assert v.convergence_order(BALANCE_DTS, res) == pytest.approx(order, abs=0.2)
    # at the production step the residual is what the measured order predicts
    predicted = res[-1] * (60.0 / BALANCE_DTS[-1]) ** order
    assert v.energy_residual_simpson(method, 60.0) == pytest.approx(predicted, rel=0.1)


# ---- 5.5 moisture without sources --------------------------------------------------------

def test_moisture_decay_without_buffer():
    p = replace(P, moisture_buffer=False)
    W_0, W_o = 0.013, 0.008
    r = simulate(p, ConstantInputs(T_o=26.0, W_o=W_o), equilibrium_state(26, 26, W_0), 8 * HOUR)
    exact = v.moisture_decay_exact_no_buffer(r.t, W_0, W_o, p)
    assert np.abs(r.W - exact).max() < 1e-12
    # fitted time constant
    k = np.searchsorted(r.t, 2 * HOUR)
    tau_fit = -r.t[k] / np.log((r.W[k] - W_o) / (W_0 - W_o))
    assert tau_fit == pytest.approx(3600.0 / p.ach, rel=1e-8)


def test_moisture_decay_with_buffer_matches_linearised_solution():
    T, W_o, dW = 26.0, 0.0116, 2e-4
    x0 = equilibrium_state(T, T, W_o + dW)      # air and buffer both humid, outdoor drier
    r = simulate(P, ConstantInputs(T_o=T, W_o=W_o), x0, 24 * HOUR)
    A, x_star = v.moisture_linear_system(P, T, W_o)
    exact = v.moisture_decay_exact_buffer(r.t, x0[[IDX_W, IDX_PHI_B]], A, x_star)
    # linearisation error is second order in dW: about dW * dW/(eps + W) ~ 1e-7
    assert np.abs(r.W - exact[:, 0]).max() < 1e-3 * dW
    assert np.abs(r.phi_b - exact[:, 1]).max() < 1e-3 * (x0[IDX_PHI_B] - x_star[1])
    # one fast mode (air <-> buffer) and one slow mode (room <-> outdoors)
    tau = np.sort(-1.0 / np.linalg.eigvals(A).real)
    # coupled air-buffer mode; the buffer side alone, C_b / (G_b p_ws), would be about 820 s
    assert 350 < tau[0] < 450
    assert 12_000 < tau[1] < 18_000
    assert tau[1] > 2 * 3600.0 / P.ach     # buffering slows the room down


# ---- independent reference: implicit adaptive solver ------------------------------------

def test_matches_scipy_radau_on_full_scenario():
    S = SummerScenario()
    x0 = periodic_state(P, S, equilibrium_state(34.5, 35.3, 0.0119))
    r = simulate(P, S, x0, DAY, dt=60.0)

    # integrate segment by segment between occupancy switches, holding n_occ as the simulator does
    breaks = [0.0, 8 * HOUR, 18 * HOUR, DAY]
    x = x0.copy()
    for a, b in zip(breaks[:-1], breaks[1:]):
        held = S(a).n_occ
        sol = solve_ivp(lambda t, y: rhs(t, y, U_OFF, S(t)._replace(n_occ=held), P),
                        (a, b), x, method="Radau", rtol=1e-11, atol=[1e-9, 1e-9, 1e-13, 1e-11, 1e-7])
        x = sol.y[:, -1]
    assert r.x[-1, [IDX_TA, IDX_TM]] == pytest.approx(x[[IDX_TA, IDX_TM]], abs=1e-6)
    assert r.x[-1, IDX_W] == pytest.approx(x[IDX_W], abs=1e-10)
    assert r.x[-1, IDX_CO2] == pytest.approx(x[IDX_CO2], abs=1e-6)
