"""SEER by the bin method of Regulation (EU) No 206/2012 (decision D46), and the supporting
changes of D48 (thermostat schedule, integrated net cooling)."""

import numpy as np
import pytest

from roomsim import reference_unit
from roomsim.ac_unit import CP_V, H_FG0, ACUnit, net_cooling
from roomsim.controllers import OnOffThermostat
from roomsim.inputs import HOUR, ConstantInputs, SummerScenario
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.seasonal import (
    BIN_HOURS,
    BIN_T,
    CDC_DEFAULT,
    DeclaredPoint,
    W_at_wet_bulb,
    eer_bin_anchor,
    local_bin_hours,
    part_load,
    seer,
    unit_declared_points,
)
from roomsim.simulate import simulate

DATASHEET = reference_unit.PUBLISHED
PUBLISHED = DATASHEET.declared_points()
MODES = DATASHEET.low_power_modes()


def test_reference_season():
    assert BIN_HOURS.sum() == 2602                         # Regulation, Annex II, Table 1
    assert BIN_T[0] == 17 and BIN_T[-1] == 40
    assert float(part_load(35.0, 2.5)) == 2.5
    assert float(part_load(16.0, 2.5)) == 0.0


def test_reproduces_the_declared_seer_of_a_real_unit():
    # the declared SEER (8.65) and annual consumption (101 kWh) of the reference unit follow from
    # its declared points and low-power modes by the method of the Regulation (decision D46)
    r = seer(PUBLISHED, DATASHEET.P_design_kW, MODES)
    assert round(r.SEER, 2) == DATASHEET.SEER == 8.65
    assert round(r.Q_CE) == 101
    assert r.Q_C == 2.50 * 350
    # points A-C match the part load within 10 %; D (1.17 kW against 0.53 kW) is degraded
    assert r.EER_bin[:3] == pytest.approx([4.50, 6.52, 10.17])
    Pc_D = 2.50 * 4 / 19
    assert r.EER_bin[3] == pytest.approx(16.51 * (1 - CDC_DEFAULT * (1 - Pc_D / 1.17)))


def test_eer_bin_rules():
    P = 1.9                                                # part load at 30 degC: 1.4 kW
    assert eer_bin_anchor(DeclaredPoint(30.0, 1.45, 5.0), P) == 5.0
    assert eer_bin_anchor(DeclaredPoint(30.0, 2.0, 5.0), P) \
        == pytest.approx(5.0 * (1 - 0.25 * (1 - 1.4 / 2.0)))
    assert eer_bin_anchor(DeclaredPoint(30.0, 2.0, 5.0), P, Cdc=0.0) == 5.0
    with pytest.raises(ValueError, match="below the part load"):
        eer_bin_anchor(DeclaredPoint(30.0, 1.0, 5.0), P)


def test_constant_efficiency_gives_that_efficiency():
    pts = [DeclaredPoint(T, float(part_load(T, 3.0)), 4.0) for T in (35.0, 30.0, 25.0, 20.0)]
    r = seer(pts, 3.0)
    assert abs(r.SEER_on - 4.0) < 1e-13
    assert abs(r.SEER - 4.0) < 1e-13


def test_seer_validation():
    with pytest.raises(ValueError, match="exactly one"):
        seer(PUBLISHED[:3], 2.5)
    with pytest.raises(ValueError, match="exactly one"):
        seer((*PUBLISHED[:3], DeclaredPoint(25.0, 1.19, 10.17)), 2.5)
    with pytest.raises(ValueError, match="non-negative"):
        seer(PUBLISHED, 2.5, bin_hours=-BIN_HOURS)


def test_local_bins_round_to_the_nearest_degree():
    T = np.array([16.49, 16.5, 17.49, 39.6, 40.49, 40.5, 25.0, 25.0])
    h = local_bin_hours(T)
    assert h.sum() == 6                                    # 16.49 and 40.5 fall outside
    assert h[0] == 2 and h[1] == 0 and h[8] == 2 and h[-1] == 2


def test_declared_points_of_the_unit_model():
    unit = ACUnit()
    P_design, pts = unit_declared_points(unit)
    W_in = W_at_wet_bulb(27.0, 19.0)
    full = net_cooling(unit.operate(1.0, 1.0, 27.0, W_in, 35.0), 27.0) / 1e3
    assert P_design == pytest.approx(full, rel=1e-14)
    assert pts[0].Pdc == pytest.approx(P_design, rel=1e-14)
    for pt in pts[1:]:
        Pc = float(part_load(pt.T_o, P_design))
        s_min_cap = net_cooling(unit.operate(unit.p.s_min, 1.0, 27.0, W_in, pt.T_o), 27.0) / 1e3
        # either the part load is met, or the unit runs at its minimum capacity above it
        assert pt.Pdc == pytest.approx(max(Pc, s_min_cap), rel=1e-9)


def test_thermostat_follows_a_set_point_schedule():
    p, unit = RoomParams(), ACUnit()
    hot = ConstantInputs(T_o=32.0, W_o=0.012, G_sol=300.0, n_occ=2)
    x0 = equilibrium_state(26.0, 26.0, 0.011)

    def schedule(t):
        return 28.0 if t < 2 * HOUR else 24.0

    r = simulate(p, hot, x0, 4 * HOUR, controller=OnOffThermostat(T_set=schedule), unit=unit)
    s = r.unit.cmd[:, 0]
    first, second = r.t < 2 * HOUR, r.t >= 3 * HOUR
    assert not np.any(s[first] > 0)                       # 26 degC is inside 28 +- 0.5 K
    assert np.any(s[second] > 0)
    assert r.T_a[second].max() <= 24.5 + 0.2               # held near the new set point
    # a constant set point behaves exactly as before
    a = simulate(p, SummerScenario(), x0, 6 * HOUR, controller=OnOffThermostat(T_set=25.0),
                 unit=unit)
    b = simulate(p, SummerScenario(), x0, 6 * HOUR,
                 controller=OnOffThermostat(T_set=lambda t: 25.0), unit=unit)
    assert np.array_equal(a.x, b.x)


def test_net_cooling_is_integrated_with_the_state():
    # the integral agrees with the trapezoid of the logged net cooling on a smooth run, as the
    # electrical energy does (tests/test_unit_in_loop.py)
    p, unit = RoomParams(), ACUnit()
    hot = ConstantInputs(T_o=32.0, W_o=0.0140, G_sol=200.0, n_occ=2)

    class Fixed:
        def step(self, t, y):
            return np.array([0.7, 1.0, 0.0])

    r = simulate(p, hot, equilibrium_state(30.0, 31.0, 0.0135), 3 * HOUR, controller=Fixed(),
                 unit=unit, unit_running_at_start=True)
    op = r.unit.op
    q = op["Q_sens"] + op["m_lat"] * (H_FG0 + CP_V * r.T_a) - op["P_fan_in"]
    trap = np.sum(0.5 * (q[1:] + q[:-1])) * r.dt
    assert r.unit.E_cool[-1] == pytest.approx(trap, rel=1e-5)
    w = r.window(HOUR)
    assert w.unit.E_cool[0] == 0.0
    assert w.unit.E_cool[-1] == pytest.approx(r.unit.E_cool[-1] - r.unit.E_cool[60], rel=1e-12)
