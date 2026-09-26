"""Heat pump in heating: cycle, frost, defrost, controllers and SCOP (decisions D50-D57)."""

import math
from itertools import pairwise

import numpy as np
import pytest
from scipy.optimize import brentq

from roomsim import psychro
from roomsim.ac_unit import ACUnit
from roomsim.controllers import (
    DefrostControl,
    OnOffThermostat,
    VentilationSchedule,
    make_heating_pi,
)
from roomsim.heat_pump import H_FUS, T0_K, HeatPump, W_sat_surface
from roomsim.inputs import DAY, HOUR, ConstantInputs
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.reference_unit import PUBLISHED
from roomsim.seasonal import (
    HEATING_BIN_HOURS,
    DeclaredPoint,
    heat_pump_declared_points,
    heating_part_load,
    scop,
)
from roomsim.simulate import simulate

HP = HeatPump()
P = RoomParams()
W_IN = 0.0073


def humid(T_o, RH):
    return float(psychro.W_from_RH(T_o, RH))


def test_saturation_over_ice():
    # Murphy and Koop (2005): 259.9 Pa at -10 degC and 611.2 Pa at 0 degC over ice; the Magnus
    # form of Alduchov and Eskridge is within 0.1 % there, and below water saturation
    assert float(psychro.p_wsi(-10.0)) == pytest.approx(259.9, rel=1e-3)
    assert float(psychro.p_wsi(0.0)) == pytest.approx(611.21)
    assert float(psychro.p_wsi(-10.0)) < float(psychro.p_ws(-10.0))
    # the surface saturation blends ice and water between -0.5 and 0.5 degC, continuously
    for T in (-0.5, 0.5):
        assert W_sat_surface(T - 1e-9) == pytest.approx(W_sat_surface(T + 1e-9), rel=1e-7)


@pytest.mark.parametrize("s,T_o,RH", [(0.3, 7.0, 0.8), (1.0, 2.0, 0.9), (1.68, -7.0, 0.8),
                                      (1.0, 12.0, 0.6), (0.6, -15.0, 0.7)])
def test_cycle_balances(s, T_o, RH):
    st = HP.operate(s, 1.0, 20.0, W_IN, T_o, W_o=humid(T_o, RH))
    G = HP.indoor_conductance(1.0, W_IN)
    coil, Q_h, cop = HP.outdoor_coil(s, G, 20.0, T_o, humid(T_o, RH), 0.0)
    P_comp = st.P_el - st.P_fan_in - HP.hp.P_fan_out
    # first law of the cycle and the heat delivered
    assert -st.Q_sens == pytest.approx(Q_h, rel=1e-12)
    assert Q_h == pytest.approx(coil.Q_e + P_comp, rel=1e-9)
    # the condenser temperature follows from the indoor coil, and the COP from it
    T_c = 20.0 + Q_h / G
    eta = HP.unit.eta_full * HP.eta_rel(s)
    lift = max(T_c - coil.T_e, HP.p.dT_lift_min)
    assert cop == pytest.approx(eta * (T_c + T0_K) / lift, rel=1e-9)
    # the outdoor coil gives what the compressor draws
    Q_air = HP.outdoor_air(coil.T_s, T_o, humid(T_o, RH), 0.0)[0]
    assert Q_air == pytest.approx(coil.Q_e, rel=1e-9)
    assert coil.T_e < coil.T_s < T_o


def test_capacity_and_cop_trends():
    # more speed, more heat; colder air, less heat and a lower COP at the same speed
    for T_o in (-10.0, 2.0, 7.0):
        q = [HP.capacity(s, T_o) for s in (0.3, 0.6, 1.0, 1.68)]
        assert all(b > a for a, b in pairwise(q))
    cops = []
    for T_o in (-15.0, -7.0, 2.0, 7.0, 12.0):
        st = HP.operate(1.0, 1.0, 20.0, W_IN, T_o)
        cops.append(-st.Q_sens / st.P_el)
    assert all(b > a for a, b in pairwise(cops))


def test_capacity_falls_with_condensing_temperature():
    # D57: at the cooling rating's condensing temperature the factor is 1
    G = HP.indoor_conductance(1.0, W_IN)
    T_in = HP.T_c_rated - 30.0
    Q_h, Q_e, _ = HP.cycle(1.0, -2.0, G, T_in)
    T_c = T_in + Q_h / G
    factor = 1.0 - HP.p.cap_slope * (T_c - HP.T_c_rated)
    assert Q_e == pytest.approx(HP.evaporator_capacity(1.0, -2.0) * factor, rel=1e-12)


def test_minimum_lift_switch_is_continuous():
    # find the outdoor temperature where the lift reaches its minimum at low speed and warm room
    G = HP.indoor_conductance(1.0, W_IN)

    def lift(T_e):
        Q_h = HP.cycle(0.3, T_e, G, 20.0)[0]
        return 20.0 + Q_h / G - T_e

    T_star = float(brentq(lambda T: lift(T) - HP.p.dT_lift_min, 0.0, 19.0))
    below, above = HP.cycle(0.3, T_star - 1e-7, G, 20.0), HP.cycle(0.3, T_star + 1e-7, G, 20.0)
    assert below[0] == pytest.approx(above[0], rel=1e-6)
    assert below[2] == pytest.approx(above[2], rel=1e-6)


def test_frost_grows_only_on_a_surface_below_zero_and_the_frost_point():
    # dry air: no deposit; warm humid air with the surface above 0 degC: condensate drains
    assert HP.operate(1.0, 1.0, 20.0, W_IN, 2.0, W_o=0.0).dM_coil == 0.0
    st = HP.operate(1.0, 1.0, 20.0, W_IN, 8.0, W_o=humid(8.0, 0.95))
    assert st.T_s > 0.5 and st.dM_coil == 0.0
    # humid air near 0 degC at maximum speed: frost at the mass-transfer rate
    W_o = humid(1.0, 0.9)
    st = HP.operate(1.68, 1.0, 20.0, W_IN, 1.0, W_o=W_o)
    assert st.T_s < -0.5
    expected = HP.hp.UA_air_out / 1006.0 * (W_o - W_sat_surface(st.T_s))
    assert st.dM_coil == pytest.approx(expected, rel=1e-9)


def test_frost_lowers_capacity_and_coil_temperature():
    W_o = humid(1.0, 0.9)
    clean = HP.operate(1.0, 1.0, 20.0, W_IN, 1.0, W_o=W_o, M_c=0.0)
    iced = HP.operate(1.0, 1.0, 20.0, W_IN, 1.0, W_o=W_o, M_c=HP.hp.M_half)
    assert -iced.Q_sens < -clean.Q_sens
    assert iced.T_s < clean.T_s


def test_defrost_and_off_states():
    h = HP.hp
    d = HP.operate(-1.0, 1.0, 20.0, W_IN, 0.0, M_c=0.5)
    assert d.P_el == h.P_defrost and d.P_fan_in == 0.0
    assert d.Q_sens == pytest.approx(h.Q_defrost - h.P_defrost)       # taken from the room
    assert d.dM_coil == pytest.approx(-h.eta_melt * h.Q_defrost / H_FUS)
    assert d.T_s == 0.0                                                  # melting plateau
    assert HP.operate(-1.0, 1.0, 20.0, W_IN, 0.0, M_c=0.0).T_s == h.T_coil_defrosted
    # compressor off: natural-convection melting in air above 0 degC, none below
    off = HP.operate(0.0, 0.0, 20.0, W_IN, 5.0, M_c=0.5)
    assert off.P_el == 0.0 and off.T_s == 5.0
    assert off.dM_coil == pytest.approx(-h.UA_natural_out * 5.0 / H_FUS)
    assert HP.operate(0.0, 0.0, 20.0, W_IN, -3.0, M_c=0.5).dM_coil == 0.0
    # melting stops smoothly as the frost vanishes
    assert HP.operate(-1.0, 1.0, 20.0, W_IN, 0.0, M_c=0.0).dM_coil == 0.0


def test_heating_pi_holds_the_set_point_and_the_balances_close():
    S = ConstantInputs(T_o=-5.0, W_o=humid(-5.0, 0.8), G_sol=0.0, n_occ=2)
    vent = VentilationSchedule()
    c = make_heating_pi(P, HP, ventilation=vent)
    r = simulate(P, S, equilibrium_state(20.0, 20.0, 0.005), 2 * DAY, controller=c, unit=HP)
    last = r.t >= DAY
    assert abs(float(r.T_a[last].mean()) - 20.0) < 0.1
    E_abs = P.C_a * 40 + P.C_m * 40
    assert np.abs(r.energy_residual).max() < 1e-11 * E_abs
    assert r.unit.M_coil.min() >= -1e-6                  # frost stays non-negative
    # the heat delivered is the integrated net cooling with the sign reversed
    assert -r.unit.E_cool[-1] > 0.0


def test_thermostat_in_heating():
    S = ConstantInputs(T_o=0.0, W_o=humid(0.0, 0.8), G_sol=0.0, n_occ=0)
    c = OnOffThermostat(T_set=20.0, heating=True, s_on=HP.hp.s_max)
    r = simulate(P, S, equilibrium_state(18.0, 18.0, 0.005), 6 * HOUR, controller=c, unit=HP)
    s = r.unit.cmd[:, 0]
    assert s[0] > 0                                      # 18 degC is below 20 - 0.5
    late = r.t >= 3 * HOUR
    # it cycles around the set point: starts and stops, mean within 0.3 K of 20 degC
    on = s[late] > 0
    assert on.any() and (~on).any()
    assert abs(float(r.T_a[late].mean()) - 20.0) < 0.3


class _Fixed:
    Ts = 60.0

    def __init__(self, s):
        self.s = s

    def step(self, t, y):
        return np.array([self.s, 1.0, 0.0])


def _supervised(strategy, s, T_o, RH, hours, **kw):
    S = ConstantInputs(T_o=T_o, W_o=humid(T_o, RH), G_sol=0.0, n_occ=0)
    c = DefrostControl(_Fixed(s), HP, strategy, ventilation=VentilationSchedule(V_occupied=0.0),
                       **kw)
    r = simulate(P, S, equilibrium_state(20.0, 20.0, 0.005), hours * HOUR, controller=c, unit=HP)
    return c, r


def test_timed_defrost_interval_and_cycle():
    # dry cold air: no frost, but the coil is cold, so the timed strategy defrosts every 40 min
    # of running; each defrost ends at once (no frost) and is followed by the drip time
    c, r = _supervised("timed", 1.0, -5.0, 0.0, 3)
    s = r.unit.cmd[:-1, 0]
    starts = np.flatnonzero((s < 0) & (np.r_[0.0, s[:-1]] >= 0))
    assert c.defrosts == len(starts) == 4
    assert np.all(np.diff(starts) == 40 + 1 + 1)          # 40 min running, 1 defrost, 1 drip
    assert r.unit.M_coil.max() == 0.0


def test_demand_defrost_waits_for_frost():
    # dry cold air: the coil matches its clean reference, so only the 3 h limit forces a defrost
    c, _ = _supervised("demand", 1.0, -5.0, 0.0, 3)
    assert c.defrosts == 0
    c, _ = _supervised("demand", 1.0, -5.0, 0.0, 3.2)
    assert c.defrosts == 1
    # humid air near 0 degC at full speed: frost builds and the coil falls below the clean
    # reference; the demand strategy defrosts and melts it
    c, r = _supervised("demand", 1.68, 1.0, 0.95, 4)
    assert c.defrosts >= 1
    s = r.unit.cmd[:-1, 0]
    first = int(np.flatnonzero(s < 0)[0])
    assert r.unit.M_coil[first] > 0.05
    after = first + int(np.flatnonzero(s[first:] >= 0)[0])
    assert r.unit.M_coil[after] < 0.02


def test_defrost_control_validation():
    with pytest.raises(ValueError, match="strategy"):
        DefrostControl(_Fixed(1.0), HP, "never")


def test_reproduces_the_declared_scop_of_a_real_unit():
    # FTXM25R + RXM25R, average season (D55): declared SCOP 5.10 and 659 kWh a year; the method
    # gives 5.12 and 656 kWh (residual 0.4 %, not explained), within the stated +-0.03
    d = PUBLISHED
    r = scop(d.heating_declared_points(), d.P_design_heating_kW, d.P_dh_design_kW,
             d.heating_low_power_modes())
    assert abs(r.SEER - d.SCOP) <= 0.03 and d.SCOP == 5.10
    assert abs(r.Q_CE / 659.0 - 1.0) < 0.005
    assert HEATING_BIN_HOURS.sum() == 4910


def test_scop_with_constant_cop_and_no_back_up():
    pts = [DeclaredPoint(T, float(heating_part_load(T, 3.0)), 3.0) for T in (-7.0, 2.0, 7.0,
                                                                             12.0)]
    r = scop(pts, 3.0, 3.0)
    assert abs(r.SEER_on - 3.0) < 1e-12 and abs(r.SEER - 3.0) < 1e-12
    with pytest.raises(ValueError, match="exactly one"):
        scop(pts[:3], 3.0, 3.0)


def test_declared_heating_points_of_the_model():
    P_d, pts = heat_pump_declared_points(HP)
    full = -HP.operate(HP.hp.s_max, 1.0, 20.0, W_IN, -10.0).Q_sens
    assert P_d * 1e3 == pytest.approx(full, rel=1e-3)    # rating humidity vs W_IN only
    for pt in pts:
        Ph = float(heating_part_load(pt.T_o, P_d))
        assert pt.Pdc == pytest.approx(Ph, rel=1e-9) or pt.Pdc > Ph
        assert math.isfinite(pt.EERd) and pt.EERd > 1.0


def test_cooling_unit_ignores_outdoor_humidity():
    u = ACUnit()
    a = u.operate(0.7, 1.0, 26.0, 0.011, 30.0)
    b = u.operate(0.7, 1.0, 26.0, 0.011, 30.0, W_o=0.02)
    assert a == b
