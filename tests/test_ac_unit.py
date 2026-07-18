"""Split air-conditioner model: energy identities, coil physics and performance trends."""

import math

import numpy as np
import pytest

from roomsim import psychro
from roomsim.ac_unit import (
    CP_A,
    CP_V,
    CP_W,
    H_FG0,
    OFF_STATE,
    ACUnit,
    T_from_h_sat,
    _dew_point,
    _h,
    _h_sat_and_slope,
    _W_sat,
    net_cooling,
)

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
    assert st.wet_fraction == 1.0
    # calibration: net EER, the indoor fan heat (which stays in the room) subtracted from the capacity
    assert net_cooling(st, R.T_in_rated) / st.P_el == pytest.approx(R.EER_rated, rel=1e-12)
    assert 0.68 < shr(st) < 0.80                                           # typical split unit
    assert 8.0 < st.T_s < 14.0
    assert 12.0 < st.T_out < 17.0
    assert st.Q_sens + st.m_lat * (H_FG0 + CP_V * R.T_in_rated) == pytest.approx(R.Q_rated)


# ---- energy and mass identities ----------------------------------------------------------

CASES = [(1.0, 1.0, 27.0, 0.0114, 35.0), (0.5, 0.7, 24.0, 0.0095, 30.0),
         (0.3, 1.0, 26.0, 0.0060, 28.0), (0.8, 0.4, 28.0, 0.0140, 40.0),
         (0.3, 1.0, 25.0, 0.0127, 32.0)]                     # last: partially wet coil


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
    assert st.wet_fraction == 0.0
    assert st.m_lat == 0.0
    assert st.W_out == 0.0050


def test_partially_wet_coil_is_continuous_from_dry_to_wet():
    """Sweep inlet humidity: the wetted fraction rises continuously from 0 to 1."""
    W = np.linspace(0.0080, 0.0150, 2801)
    st = [U.operate(0.3, 1.0, 25.0, w, 32.0) for w in W]
    wet = np.array([s.wet_fraction for s in st])
    Qs = np.array([s.Q_sens for s in st])
    ml = np.array([s.m_lat for s in st])
    assert wet[0] == 0.0 and wet[-1] == 1.0
    partial = (wet > 0.0) & (wet < 1.0)
    assert partial.sum() > 100                           # a wide range of partially wet states
    assert np.all(np.diff(wet) >= 0.0) and np.all(np.diff(ml) >= -1e-15)
    # no jump anywhere: no step is much larger than both of its neighbours
    for y in (wet, Qs, ml):
        dy = np.abs(np.diff(y))
        assert np.all(dy[1:-1] <= 2.0 * np.maximum(dy[:-2], dy[2:]) + 1e-12 * np.abs(y).max())
    assert np.all(ml[wet == 0.0] == 0.0) and np.all(ml[partial] > 0.0)


def test_partially_wet_coil_meets_both_limits_at_its_boundaries():
    """Coil at a fixed evaporating temperature: no jump where the wet part appears (surface at
    the dew point at the air outlet) or where the dry part vanishes (surface at the dew point at
    the air inlet)."""
    T_in, W_in, f = 25.0, 0.0127, 1.0
    T_dp = _dew_point(W_in)
    m = R.m_nom * f
    UA_a = R.NTU_nom * m * (CP_A + CP_V * W_in)
    r = UA_a / (UA_a + R.UA_evap)
    NTU_o = R.NTU_nom * R.UA_evap / (UA_a + R.UA_evap)
    # wet part appears: T_r + r (T_in - T_r) exp(-NTU_o) = T_dp
    k = r * np.exp(-NTU_o)
    T_first = (T_dp - k * T_in) / (1.0 - k)
    # dry part vanishes: T_r + r (T_in - T_r) = T_dp
    T_all = (T_dp - r * T_in) / (1.0 - r)
    for T_r, expect in ((T_first, (0.0, 0.0)), (T_all, (1.0, 1.0))):
        warm, cold = (U.coil_at(T_r + d, f, T_in, W_in) for d in (1e-7, -1e-7))
        assert (warm.wet_fraction, cold.wet_fraction) == pytest.approx(expect, abs=1e-5)
        assert pytest.approx(warm.Q, rel=1e-6) == cold.Q
        assert cold.W_out == pytest.approx(warm.W_out, abs=1e-9)
        assert cold.T_out == pytest.approx(warm.T_out, abs=1e-5)
    assert U.coil_at(T_first + 1e-3, f, T_in, W_in).W_out == W_in       # still dry


@pytest.mark.parametrize("T_r", [5.0, 12.0, 15.5, 16.5, 17.0, 20.0])
def test_coil_outlet_is_physical(T_r):
    """Air-side enthalpy drop equals the capacity; outlet not supersaturated and no water added;
    surface and outlet temperatures between the refrigerant and the inlet air."""
    T_in, W_in = 25.0, 0.0127
    out = U.coil_at(T_r, 1.0, T_in, W_in)
    assert pytest.approx(R.m_nom * (_h(T_in, W_in) - _h(out.T_out, out.W_out)), rel=1e-12) == out.Q
    assert out.W_out <= W_in and out.W_out <= _W_sat(out.T_out) * (1 + 1e-12)
    assert T_r < out.T_s < T_in and T_r < out.T_out < T_in
    assert 0.0 <= out.wet_fraction <= 1.0


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
    T_in, W_in, T_o = 22.0, 0.0080, 30.0                  # cool, dry return air, minimum fan
    st = U.operate(1.0, 0.4, T_in, W_in, T_o)
    Q_air = st.Q_sens + st.m_lat * (H_FG0 + CP_V * T_in)
    assert U.T_evap(Q_air, 0.4, T_in, W_in) == pytest.approx(R.T_evap_min, abs=1e-9)
    assert Q_air < U.compressor_capacity(1.0, R.T_evap_min, T_o)     # the compressor slows
    # no cap when the fan is at full speed: compressor and coil capacities balance
    T_r = U.operating_T_r(1.0, 1.0, T_in, W_in, T_o)
    assert T_r > R.T_evap_min
    st_full = U.operate(1.0, 1.0, T_in, W_in, T_o)
    Q_full = st_full.Q_sens + st_full.m_lat * (H_FG0 + CP_V * T_in)
    assert Q_full == pytest.approx(U.compressor_capacity(1.0, T_r, T_o), rel=1e-9)


def test_lower_airflow_reduces_capacity():
    """Less air: colder evaporator, lower suction density, less capacity. At rating conditions
    the EnergyPlus example curve WindACCoolCapFFF (0.8 + 0.2 ff) gives 0.94 at 70 % flow."""
    def total(s, f, T_in, W_in, T_o):
        st = U.operate(s, f, T_in, W_in, T_o)
        return st.Q_sens + st.m_lat * (H_FG0 + CP_V * T_in)
    Q = [total(1.0, f, 27.0, 0.0114, 35.0) for f in (1.0, 0.85, 0.7, 0.55, 0.4)]
    assert np.all(np.diff(Q) < 0.0)
    assert Q[2] / Q[0] == pytest.approx(0.94, abs=0.02)
    # at minimum speed, room at 25 degC and 60 % RH: a few per cent
    W60 = float(psychro.W_from_RH(25.0, 0.60))
    Qm = [total(0.3, f, 25.0, W60, 32.0) for f in (1.0, 0.7)]
    assert 0.92 < Qm[1] / Qm[0] < 0.99


def test_capacity_sensitivities_match_the_curve():
    """cap_slope_evap and cap_slope were set so that the capacity at rating, coil included,
    rises with the entering wet bulb by 3.58 %/K and falls with the outdoor temperature by
    1.08 %/K, as in the EnergyPlus curve WindACCoolCapFT (decision D30)."""
    def W_at_wb(t):                                       # constant-enthalpy line at 27 degC
        return (_h(t, _W_sat(t)) - CP_A * 27.0) / (H_FG0 + CP_V * 27.0)

    def total(W_in, T_o):
        st = U.operate(1.0, 1.0, 27.0, W_in, T_o)
        return st.Q_sens + st.m_lat * (H_FG0 + CP_V * 27.0)
    wb = T_from_h_sat(_h(27.0, R.W_in_rated), T_guess=20.0)
    Q = [total(W_at_wb(wb + d), 35.0) for d in (-0.05, 0.05)]
    slope = (Q[1] - Q[0]) / 0.1 / R.Q_rated
    assert slope == pytest.approx(0.0358, abs=0.001)      # cap_slope_evap is rounded to 0.040
    Q = [total(R.W_in_rated, T_o) for T_o in (34.5, 35.5)]
    slope = (Q[1] - Q[0]) / R.Q_rated
    assert slope == pytest.approx(-0.0108, abs=0.0004)    # cap_slope is rounded to 0.016


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
    slow = U.operate(0.0, R.f_min, 27.0, 0.0114, 35.0)
    P_slow = R.P_fan_in_base + (R.P_fan_in_nom - R.P_fan_in_base) * R.f_min ** 3
    assert slow.P_el == pytest.approx(P_slow)
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


def test_cop_is_bounded_by_the_minimum_lift_and_continuous():
    """The cycle never lets the condensing temperature approach the evaporating one."""
    for T_o in np.linspace(0.0, 45.0, 46):
        for s in (0.3, 0.6, 1.0):
            for f in (0.4, 1.0):
                st = U.operate(s, f, 25.0, 0.0115, T_o)
                T_ev = U.operating_T_r(s, f, 25.0, 0.0115, T_o) + 273.15
                cop_max = U.eta_full * U.eta_rel(s) * T_ev / R.dT_lift_min
                assert 0.0 < st.COP <= cop_max * (1 + 1e-12)
                assert st.P_el > U.fan_power(f) + U.p.P_fan_out
    # no jump where the floor takes over
    T_o = np.linspace(15.0, 35.0, 2001)
    cop = np.array([U.operate(0.3, 1.0, 25.0, 0.0115, t).COP for t in T_o])
    assert np.abs(np.diff(cop)).max() < 0.01
