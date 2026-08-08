"""Regression tests for robustness fixes in the plant, inputs, psychrometrics and simulation loop.

None of these cases occurs in the published studies (default buffer, dt = 60 s, non-negative
flows), so the fixes leave every published number unchanged.
"""

import ast
import inspect
import warnings
from dataclasses import replace

import numpy as np
import pytest

from roomsim import psychro
from roomsim.inputs import DAY, HOUR, ConstantInputs, SummerScenario
from roomsim.params import RoomParams
from roomsim.plant import IDX_PHI_B, equilibrium_state
from roomsim.simulate import periodic_state, simulate
from roomsim.ventilation import VentilationFan
from roomsim.verification import one_node_params

P = RoomParams()


# ---- F1 / F5: buffer conductance consistent with its capacity ---------------------------

@pytest.mark.parametrize("field", ["buffer_slope", "buffer_depth", "buffer_density",
                                   "buffer_area"])
def test_zero_capacity_buffer_exchanges_no_water(field):
    """A buffer layer that cannot store water (C_b = 0) cannot exchange water with the air
    either: otherwise phi_b is frozen (dphi_b/dt = 0) while m_buf = G_b (phi_b p_ws - p_w)
    keeps feeding the air, an unlimited water source. Zero capacity must behave exactly like
    ``moisture_buffer=False``, and the room water balance must close to round-off."""
    p = replace(P, **{field: 0.0})
    assert p.buffer_capacity == 0.0
    assert p.buffer_conductance == 0.0
    x0 = equilibrium_state(26.0, 26.0, 0.0115)
    x0[IDX_PHI_B] = 0.9                       # inconsistent with the air; holds no water
    d = ConstantInputs(T_o=26.0, W_o=0.0115)
    r = simulate(p, d, x0, DAY)
    r_off = simulate(replace(P, moisture_buffer=False), d, x0, DAY)
    assert r.W[-1] == r_off.W[-1]
    assert np.max(np.abs(r.water_residual)) < 1e-12


@pytest.mark.parametrize("field", ["h_i", "buffer_permeability"])
def test_buffer_conductance_with_zero_transfer_coefficient(field):
    """G_b = A / (1/beta + d/(2 delta)) tends to 0 as beta -> 0 (h_i = 0, the documented way to
    decouple air and mass) or delta -> 0 (impermeable layer). The property must return that
    limit instead of raising ZeroDivisionError."""
    assert replace(P, **{field: 0.0}).buffer_conductance == 0.0


def test_default_buffer_conductance_unchanged():
    """The published runs use the default buffer; its conductance must keep every bit."""
    assert P.buffer_conductance == 5.200716408323909e-07


@pytest.mark.parametrize("field", ["buffer_slope", "buffer_depth", "buffer_density",
                                   "buffer_area", "buffer_permeability", "h_i", "ach"])
def test_negative_room_parameter_rejected(field):
    """Capacities, conductances and the infiltration rate are non-negative by definition; a
    negative value would, for example, give a negative buffer capacity that freezes phi_b while
    the buffer still exchanges water."""
    with pytest.raises(ValueError, match=field):
        replace(P, **{field: -1.0})


# ---- F2: periodic shooting with a state that has no dynamics ----------------------------

def test_periodic_state_with_frozen_buffer_state():
    """With the buffer off, phi_b has no dynamics (it passes through a period unchanged), so
    the Newton matrix Jf - I has a zero row and column. The shooting must solve for the other
    states and leave phi_b at its guess, for any guess including phi_b = 0 (where the
    finite-difference round-off that previously hid the singularity is exactly zero)."""
    p = replace(P, moisture_buffer=False)
    S = SummerScenario()
    results = []
    for phi in (0.0, 0.5):
        x_guess = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
        x_guess[IDX_PHI_B] = phi
        x_p = periodic_state(p, S, x_guess)
        assert x_p[IDX_PHI_B] == phi
        results.append(x_p)
    others = [0, 1, 2, 4]
    assert results[0][others] == pytest.approx(results[1][others], abs=1e-6)
    r = simulate(p, S, results[0], DAY)
    assert r.x[-1] == pytest.approx(results[0], abs=1e-6)


def test_periodic_state_open_integrator_raises_clearly():
    """With h_i = 0 and H_mo = 0 the mass node only accumulates the solar and radiant gains,
    so no periodic state exists. The shooting must say so instead of failing inside the
    linear solver or returning a meaningless state."""
    S = SummerScenario()
    x_guess = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
    with pytest.raises(RuntimeError, match="no periodic state"):
        periodic_state(one_node_params(), S, x_guess)


# ---- F3: divergence of the explicit integrator ------------------------------------------

def test_simulate_raises_on_unstable_step():
    """RK4 is unstable above about 1100 s for the default room (fastest mode: air-buffer
    moisture exchange). At dt = 1200 s the humidity ratio turns negative within hours; the
    simulator must stop with an error naming the time and the stability limit instead of
    returning W in [-12, 19] kg/kg."""
    S = SummerScenario()
    x0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
    with pytest.raises(RuntimeError, match=r"t = .* s.*stability limit"):
        simulate(P, S, x0, 2 * DAY, dt=1200.0)


def test_simulate_stable_long_step_still_runs():
    """Steps below the stability limit are unaffected by the check."""
    S = SummerScenario()
    x0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
    r = simulate(P, S, x0, 2 * DAY, dt=900.0)
    assert np.all(r.W > 0.0)


# ---- F4: negative ventilation flow ------------------------------------------------------

class _NegativeVent:
    def step(self, t, y):
        return np.array([0.0, 0.0, -0.02])


def test_negative_ventilation_flow_rejected():
    """A negative supply flow flips the sign of the ventilation heat, moisture and CO2 terms
    (35 degC outdoor air would cool the room) and makes the CO2 balance unstable once
    V_inf + V_vent < 0. A balanced supply fan cannot deliver a negative flow."""
    x0 = equilibrium_state(26.0, 26.0, 0.0115, 800.0)
    with pytest.raises(ValueError, match="V_vent"):
        simulate(P, ConstantInputs(n_occ=2, T_o=35.0), x0, HOUR, controller=_NegativeVent())


def test_fan_power_rejects_negative_flow():
    """The affinity law P = P_d (V/V_d)^3 would give a negative electrical power."""
    with pytest.raises(ValueError, match="V"):
        VentilationFan().power(np.array([0.01, -0.01]))


# ---- F6: limits and input validation ----------------------------------------------------

def test_dew_point_of_dry_air_is_minus_infinity():
    """p_w -> 0 gives T_dew -> -inf (inverse Magnus); dry air must return that limit, without
    the NaN and divide-by-zero warnings of log(0)/(B - log(0))."""
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert psychro.dew_point_from_W(0.0) == -np.inf
        T = psychro.dew_point_from_W(np.array([0.0, 0.01]))
    assert T[0] == -np.inf
    assert T[1] == psychro.dew_point_from_W(0.01)


def test_tau_moisture_sealed_room_is_infinite():
    """A sealed room (ach = 0) has no air exchange, so its moisture time constant is
    infinite, not a ZeroDivisionError."""
    assert replace(P, ach=0.0).tau_moisture == np.inf


@pytest.mark.parametrize("RH_o", [-0.1, 1.2])
def test_outdoor_rh_outside_unit_interval_rejected(RH_o):
    """RH is a fraction 0..1; 55 written as percent would give W_o far above saturation."""
    with pytest.raises(ValueError, match="RH_o"):
        SummerScenario(RH_o=RH_o)


@pytest.mark.parametrize("sunset_h", [6.0, 5.0])
def test_sunset_before_sunrise_rejected(sunset_h):
    """The half-sine divides by the day length sunset_h - sunrise_h, which must be positive."""
    with pytest.raises(ValueError, match="sunset_h"):
        SummerScenario(sunrise_h=6.0, sunset_h=sunset_h)


def test_occupancy_interval_across_midnight_wraps():
    """An interval with start > end (night shift 22 h to 6 h) is occupied from start to
    midnight and from midnight to end; previously it silently gave zero occupants."""
    S = SummerScenario(occupancy=((22.0, 6.0, 1.0), (9.0, 12.0, 2.0)))
    t = np.array([23.0, 3.0, 6.0, 10.0, 15.0, 22.0]) * HOUR
    assert list(S.n_occ(t)) == [1.0, 1.0, 0.0, 2.0, 0.0, 1.0]
    S2 = SummerScenario(occ_start_h=22.0, occ_end_h=6.0, n_people=1.0)
    assert list(S2.n_occ(t)) == [1.0, 1.0, 0.0, 0.0, 0.0, 1.0]


def test_psychro_has_no_orphan_string_literals():
    """Bare string expressions after the module docstring are dead code (unattached docs)."""
    tree = ast.parse(inspect.getsource(psychro))
    strings = [n for n in tree.body[1:]
               if isinstance(n, ast.Expr) and isinstance(n.value, ast.Constant)
               and isinstance(n.value.value, str)]
    assert strings == []
