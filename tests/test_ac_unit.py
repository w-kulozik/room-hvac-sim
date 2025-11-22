"""Split air-conditioner model: energy identities, coil physics and performance trends."""

import math
import numpy as np
import pytest

from roomsim import psychro
from roomsim.ac_unit import (CP_V, CP_W, H_FG0, OFF_STATE, ACUnit, T_from_h_sat, _h,
                             _h_sat_and_slope)

U = ACUnit()
R = U.p   # rating conditions


def shr(st):
    return st.Q_sens / (st.Q_sens + st.m_lat * (H_FG0 + CP_V * R.T_in_rated))


# ---- numerics ---------------------------------------------------------------------------

def test_saturation_enthalpy_inversion():
    for T in np.linspace(-5.0, 30.0, 15):
        h, _ = _h_sat_and_slope(T)
        assert T_from_h_sat(h, T_guess=T + 8.0) == pytest.approx(T, abs=1e-9)


def test_saturation_enthalpy_matches_psychro():
    T = 12.0
    h, dh = _h_sat_and_slope(T)
    assert h / 1e3 == pytest.approx(float(psychro.enthalpy(T, psychro.W_from_RH(T, 1.0))), rel=1e-12)
    h2, _ = _h_sat_and_slope(T + 1e-5)
    h1, _ = _h_sat_and_slope(T - 1e-5)
    assert dh == pytest.approx((h2 - h1) / 2e-5, rel=1e-7)


# ---- rating point -----------------------------------------------------------------------

def test_rated_point():
    st = U.operate(1.0, 1.0, R.T_in_rated, R.W_in_rated, R.T_o_rated)
    assert st.wet
    assert st.Q_coil / st.P_el == pytest.approx(R.EER_rated, rel=1e-12)    # calibration
    assert 0.68 < shr(st) < 0.80                                           # typical split unit
    assert 8.0 < st.T_s < 14.0
    assert 12.0 < st.T_out < 17.0
    assert st.Q_sens + st.m_lat * (H_FG0 + CP_V * R.T_in_rated) == pytest.approx(R.Q_rated)


# ---- energy and mass identities ----------------------------------------------------------

CASES = [(1.0, 1.0, 27.0, 0.0114, 35.0), (0.5, 0.7, 24.0, 0.0095, 30.0),
         (0.3, 1.0, 26.0, 0.0060, 28.0), (0.8, 0.4, 28.0, 0.0140, 40.0)]


@pytest.mark.parametrize("s, f, T_in, W_in, T_o", CASES)
def test_air_side_energy_identity(s, f, T_in, W_in, T_o):
    st = U.operate(s, f, T_in, W_in, T_o)
    m = R.m_nom * min(max(f, R.f_min), 1.0)
    dh = m * (_h(T_in, W_in) - _h(st.T_out, st.W_out))
    Q_lat = st.m_lat * (H_FG0 + CP_V * T_in)
    assert st.Q_sens + Q_lat == pytest.approx(dh, rel=1e-12)
    assert st.Q_coil == pytest.approx(dh - st.m_lat * CP_W * st.T_s, rel=1e-12)
    assert st.m_lat == pytest.approx(m * (W_in - st.W_out), rel=1e-12)


@pytest.mark.parametrize("s, f, T_in, W_in, T_o", CASES)
def test_outlet_never_supersaturated(s, f, T_in, W_in, T_o):
    st = U.operate(s, f, T_in, W_in, T_o)
    assert psychro.RH_from_W(st.T_out, st.W_out) <= 1.0 + 1e-12
    assert st.W_out <= W_in


def test_outlet_saturation_guard_over_a_grid():
    for T_in in np.linspace(20.0, 32.0, 7):
        for RH in np.linspace(0.3, 0.95, 8):
            W_in = float(psychro.W_from_RH(T_in, RH))
            for f in (0.4, 0.7, 1.0):
                st = U.operate(1.0, f, T_in, W_in, 35.0)
                assert psychro.RH_from_W(st.T_out, st.W_out) <= 1.0 + 1e-12


# ---- coil physics -------------------------------------------------------------------------

def test_dry_coil_removes_no_water():
    st = U.operate(0.5, 1.0, 26.0, 0.0050, 30.0)        # dew point ~4 degC
    assert not st.wet
    assert st.m_lat == 0.0
    assert st.W_out == 0.0050


def test_dry_wet_transition_is_continuous():
    """Sweep inlet humidity through the point where the coil starts to condense."""
    W = np.linspace(0.0070, 0.0100, 3001)
    Qs = np.array([U.operate(0.6, 1.0, 26.0, w, 30.0).Q_sens for w in W])
    ml = np.array([U.operate(0.6, 1.0, 26.0, w, 30.0).m_lat for w in W])
    wet = np.array([U.operate(0.6, 1.0, 26.0, w, 30.0).wet for w in W])
    assert wet[0] == False and wet[-1] == True
    k = np.argmax(wet)
    # no jump at the switch: the step there is no larger than the neighbouring steps
    assert abs(Qs[k] - Qs[k - 1]) < 3 * abs(Qs[k + 2] - Qs[k + 1]) + 1e-9
    assert ml[k] < 1e-7


def test_bypass_factor_scaling():
    assert U.bypass_factor(1.0) == pytest.approx(R.BF_nom)
    assert U.bypass_factor(0.5) < R.BF_nom        # slower air: more contact, less bypass


# ---- performance trends -------------------------------------------------------------------

def test_lower_fan_speed_dehumidifies_more():
    st = [U.operate(1.0, f, 27.0, 0.0114, 35.0) for f in (1.0, 0.7, 0.5)]
    assert shr(st[0]) > shr(st[1]) > shr(st[2])
    assert st[0].T_s > st[1].T_s > st[2].T_s
    assert st[0].COP > st[1].COP > st[2].COP           # colder coil costs efficiency


def test_frost_protection_caps_evaporating_temperature():
    st = U.operate(1.0, 0.4, 27.0, 0.0114, 35.0)
    Q_air = st.Q_sens + st.m_lat * (H_FG0 + CP_V * 27.0)
    assert U.T_evap(Q_air, 0.4, 27.0, 0.0114) == pytest.approx(R.T_evap_min, abs=1e-3)
    assert Q_air < U.capacity(1.0, 35.0)
    # no cap when the fan is at full speed
    st_full = U.operate(1.0, 1.0, 27.0, 0.0114, 35.0)
    assert st_full.Q_sens + st_full.m_lat * (H_FG0 + CP_V * 27.0) == pytest.approx(U.capacity(1.0, 35.0))


def test_hot_outdoor_reduces_capacity_and_efficiency():
    a, b = (U.operate(1.0, 1.0, 27.0, 0.0114, T_o) for T_o in (30.0, 40.0))
    assert b.Q_coil < a.Q_coil
    assert b.Q_coil / b.P_el < a.Q_coil / a.P_el


def test_inverter_part_load_is_more_efficient():
    eer = [U.operate(s, 1.0, 27.0, 0.0114, 35.0) for s in (1.0, 0.6)]
    assert eer[1].Q_coil / eer[1].P_el > eer[0].Q_coil / eer[0].P_el


def test_commands_are_clipped():
    assert U.operate(0.0, 0.0, 27.0, 0.0114, 35.0) is OFF_STATE
    fan = U.operate(0.0, 1.0, 27.0, 0.0114, 35.0)              # fan only, dry coil
    assert fan.P_el == pytest.approx(U.p.P_fan_in_nom) and fan.m_lat == 0.0
    assert U.operate(0.1, 1.0, 27.0, 0.0114, 35.0) == U.operate(R.s_min, 1.0, 27.0, 0.0114, 35.0)
    assert U.operate(1.5, 1.0, 27.0, 0.0114, 35.0) == U.operate(1.0, 1.0, 27.0, 0.0114, 35.0)
    assert U.operate(1.0, 0.1, 27.0, 0.0114, 35.0) == U.operate(1.0, R.f_min, 27.0, 0.0114, 35.0)


# ---- start-up ramp -----------------------------------------------------------------------

def test_start_ramp_reduces_capacity_not_power():
    full = U.operate(1.0, 1.0, 27.0, 0.0114, 35.0)
    half = U.operate(1.0, 1.0, 27.0, 0.0114, 35.0, ramp=0.5)
    assert half.Q_coil < 0.55 * full.Q_coil
    assert half.P_el == pytest.approx(full.P_el)
    assert U.start_ramp(0.0) == 0.0
    assert U.start_ramp(R.tau_start) == pytest.approx(1 - math.exp(-1))
