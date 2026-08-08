"""Regression tests for the unit model (ac_unit.py): frost-limited cycle and fan-only drying."""

import math

import pytest
from scipy.optimize import brentq

from roomsim import psychro
from roomsim.ac_unit import ACUnit

UNIT = ACUnit()
P = UNIT.p
M_MAX = P.M_coil_max


def _W(T, RH):
    return float(psychro.W_from_RH(T, RH))


def _s_eff(s, f, T_in, W_in, T_o):
    """Speed whose compressor capacity the coil takes at T_evap_min (frost-limited point)."""
    T_r = P.T_evap_min
    return s * UNIT.coil_at(T_r, f, T_in, W_in).Q / UNIT.compressor_capacity(s, T_r, T_o)


# ---- F1: frost-limited operation ----------------------------------------------------------

def test_frost_limited_cop_uses_effective_speed():
    """With T_r held at T_evap_min the compressor runs at s_eff < s, so eta_rel(s_eff) applies.

    eta_rel describes the compressor at the speed it actually runs. At 20 degC / 10 % RH return
    air and -10 degC outdoors the coil takes only 37 % of the full-speed capacity, so a cycle
    efficiency evaluated at the commanded s = 1 overstated the COP by about 25 %.
    """
    s, f, T_in, T_o = 1.0, 0.4, 20.0, -10.0
    W_in = _W(T_in, 0.10)
    assert UNIT.operating_T_r(s, f, T_in, W_in, T_o) == P.T_evap_min
    s_eff = _s_eff(s, f, T_in, W_in, T_o)
    assert P.s_min < s_eff < s
    st = UNIT.operate(s, f, T_in, W_in, T_o)
    expected = UNIT._cop(UNIT.eta_full * UNIT.eta_rel(s_eff), st.Q_coil, P.T_evap_min, T_o)
    cop = st.COP
    assert cop == pytest.approx(expected, rel=1e-12)
    assert st.P_el == pytest.approx(st.Q_coil / expected + st.P_fan_in + P.P_fan_out, rel=1e-12)


def test_frost_limited_below_s_min_is_cycled_at_s_min():
    """Below s_min the unit cycles at s_min; the logged state is the cycle average.

    A compressor cannot run continuously below s_min. Cycling at s_min with duty
    d = s_eff / s_min delivers on average the coil capacity at T_evap_min (the logged Q_coil),
    while each on-period runs at s_min with the load Q_coil / d. So the cycle is evaluated at
    eta_rel(s_min) with the on-period condenser load, and P_comp = Q_coil / COP_on. With the
    commanded s = 1 the previous code used eta_rel(1) here (COP 25 % too high).
    """
    s, f, T_in, T_o = 1.0, 0.4, 12.0, -10.0
    W_in = _W(T_in, 0.10)
    s_eff = _s_eff(s, f, T_in, W_in, T_o)
    assert s_eff < P.s_min
    st = UNIT.operate(s, f, T_in, W_in, T_o)
    duty = s_eff / P.s_min
    cop_on = UNIT._cop(UNIT.eta_full * UNIT.eta_rel(P.s_min), st.Q_coil / duty, P.T_evap_min, T_o)
    cop = st.COP
    assert cop == pytest.approx(cop_on, rel=1e-12)
    assert st.Q_coil / st.COP == pytest.approx(st.P_el - st.P_fan_in - P.P_fan_out, rel=1e-12)


def test_frost_limited_power_is_continuous():
    """P_el and COP are continuous at frost onset and where s_eff crosses s_min.

    Both switches are regime changes of the same physical unit; a jump would be a modelling
    artefact (and a discontinuous right-hand side for the integrator).
    """
    f, T_o = 0.4, -10.0
    W_in = _W(20.0, 0.10)
    # s_eff = s_min crossing, swept through the return-air temperature at s = 1
    T_x = brentq(lambda T: _s_eff(1.0, f, T, _W(T, 0.10), T_o) - P.s_min, 10.0, 20.0, xtol=1e-13)
    d = 1e-7
    lo, hi = (UNIT.operate(1.0, f, T, _W(T, 0.10), T_o) for T in (T_x - d, T_x + d))
    assert abs(hi.P_el - lo.P_el) < 1e-3 and abs(hi.COP - lo.COP) < 1e-4
    # frost onset, swept through the compressor speed at fixed air
    def g(s):
        return (UNIT.coil_at(P.T_evap_min, f, 20.0, W_in).Q
                - UNIT.compressor_capacity(s, P.T_evap_min, T_o))
    s_x = brentq(g, P.s_min, 1.0, xtol=1e-14)
    lo, hi = (UNIT.operate(s, f, 20.0, W_in, T_o) for s in (s_x - d, s_x + d))
    assert abs(hi.P_el - lo.P_el) < 1e-3 and abs(hi.COP - lo.COP) < 1e-4


def test_unconstrained_cop_still_uses_commanded_speed():
    """Away from frost protection T_r is the balance root and s_eff = s: nothing changes."""
    for s, f, T_in, W_in, T_o in [(1.0, 1.0, 27.0, 0.0114, 35.0), (0.5, 0.6, 25.0, 0.010, 30.0)]:
        T_r = UNIT.operating_T_r(s, f, T_in, W_in, T_o)
        assert T_r > P.T_evap_min
        st = UNIT.operate(s, f, T_in, W_in, T_o)
        cop = st.COP
        assert cop == UNIT._cop(UNIT.eta_full * UNIT.eta_rel(s), st.Q_coil, T_r, T_o)


def test_frost_protection_also_active_at_full_fan():
    """The module docstring no longer claims frost protection needs low fan speed.

    Cool, dry return air with cold outdoor air is frost-limited even at full fan.
    """
    W_in = _W(15.0, 0.10)
    assert UNIT.operating_T_r(1.0, 1.0, 15.0, W_in, -10.0) == P.T_evap_min


# ---- F2: water held on a dry coil while the compressor runs (documented assumption) --------

def test_held_water_is_kept_on_a_dry_coil_with_compressor_on():
    """On-cycle assumption: with the compressor on, a dry coil neither drains nor evaporates.

    Documented in the module docstring; the retained water is bounded by M_coil_max and
    re-evaporates at the next fan-only period.
    """
    st = UNIT.operate(0.3, 1.0, 25.0, 0.0099, 32.0, M_c=0.2)
    assert st.wet_fraction == 0.0
    assert st.m_lat == 0.0 and st.dM_coil == 0.0 and st.m_drain == 0.0


# ---- F3: fan-only surface temperature and wet fraction as the film vanishes ---------------

def test_fan_only_flags_follow_the_evaporating_film():
    """T_s and wet_fraction go continuously to the dry values as M_c -> 0.

    Below film_ramp * M_coil_max the evaporation rate is proportional to M_c, so M_c decays
    exponentially and never reaches 0; the surface and wet fraction must follow the same
    weight e(M_c) instead of staying at the wet values (T_wb, 1) until M_c is exactly 0.
    """
    T_in, W_in = 25.0, 0.0099
    dry = UNIT.operate(0.0, 1.0, T_in, W_in, 32.0, M_c=0.0)
    tiny = UNIT.operate(0.0, 1.0, T_in, W_in, 32.0, M_c=1e-30)
    assert dry.T_s == T_in and dry.wet_fraction == 0.0
    assert tiny.T_s == pytest.approx(T_in, abs=1e-20)
    assert tiny.wet_fraction == pytest.approx(0.0, abs=1e-20)
    full = UNIT.operate(0.0, 1.0, T_in, W_in, 32.0, M_c=0.5 * M_MAX)
    w = P.film_ramp * M_MAX
    half = UNIT.operate(0.0, 1.0, T_in, W_in, 32.0, M_c=0.5 * w)
    assert full.wet_fraction == 1.0 and full.T_s < T_in
    assert half.wet_fraction == pytest.approx(0.5, rel=1e-12)
    assert half.T_s == pytest.approx(T_in + 0.5 * (full.T_s - T_in), rel=1e-12)


# ---- F4-F6: documented behaviour of operate() -------------------------------------------

def test_operate_documented_edge_behaviour():
    """M_c=None is dry in fan-only and full (all drains) with the compressor on; s > 0 with
    f <= 0 runs the fan at f_min; during the ramp COP refers to the requested point."""
    fan = UNIT.operate(0.0, 1.0, 25.0, 0.0099, 32.0)
    assert fan.m_lat == 0.0 and fan.wet_fraction == 0.0
    on = UNIT.operate(0.8, 1.0, 26.0, 0.012, 32.0)
    assert on.dM_coil == 0.0 and on.m_drain == on.m_lat > 0.0
    assert UNIT.operate(1.0, 0.0, 27.0, 0.0114, 35.0).P_fan_in == UNIT.fan_power(P.f_min)
    req = UNIT.operate(1.0, 1.0, 27.0, 0.0114, 35.0)
    ramped = UNIT.operate(1.0, 1.0, 27.0, 0.0114, 35.0, ramp=0.2)
    assert ramped.COP == req.COP and ramped.P_el == req.P_el
    assert math.isclose(req.Q_coil / req.COP, req.P_el - req.P_fan_in - P.P_fan_out,
                        rel_tol=1e-12)
