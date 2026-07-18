"""Humidity-aware fan control: direct RH loop (H1) and coil-temperature cascade (H2)."""

import math

import numpy as np
import pytest

from criteria import KH_PER_DAY, RH_ACCURACY
from roomsim.ac_unit import OFF_STATE, ACUnit
from roomsim.controllers import OnOffThermostat, VentilationSchedule, make_pi_controller
from roomsim.controllers.humidity import CoilCascade, FanHumidityPI
from roomsim.inputs import DAY, SummerScenario
from roomsim.metrics import compute
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.simulate import coil_sensor, simulate

P = RoomParams()
UNIT = ACUnit()
F_MIN = UNIT.p.f_min


class FixedCompressor:
    """Stands in for the temperature controller: constant compressor speed, full fan."""

    def __init__(self, s):
        self.s = s

    def step(self, t, y):
        return np.array([self.s, 1.0, 0.0])


def y(T_a=25.0, RH=0.55, C=800.0, T_coil=15.0):
    return np.array([T_a, RH, C, T_coil])


# ---- sensor ----------------------------------------------------------------------------------

def test_coil_sensor_reads_the_coil_only_while_the_unit_runs():
    st = UNIT.operate(0.5, 1.0, 25.0, 0.011, 32.0)
    assert coil_sensor(st, 25.0) == st.T_s
    assert coil_sensor(OFF_STATE, 24.3) == 24.3


# ---- H1 --------------------------------------------------------------------------------------

def test_h1_leaves_the_fan_alone_when_dry_enough():
    c = FanHumidityPI(FixedCompressor(0.5), RH_set=0.55)
    for k in range(50):
        cmd = c.step(60.0 * k, y(RH=0.45))
    assert cmd[1] == 1.0


def test_h1_slows_the_fan_down_to_its_minimum_when_humid():
    c = FanHumidityPI(FixedCompressor(0.5), RH_set=0.55)
    fans = [c.step(60.0 * k, y(RH=0.70))[1] for k in range(300)]
    assert fans[0] < 1.0 and np.all(np.diff(fans) <= 1e-12)
    assert fans[-1] == pytest.approx(F_MIN)
    # under saturation the integrator settles at 1 - Kp e (1 - Tt/Ti) instead of winding up
    pi = c.pi
    assert pytest.approx(1.0 - pi.Kp * 15.0 * 0.5, rel=1e-6) == pi.I
    # once dry (e = -5 %) it unwinds at Kp Ts/Ti per %RH and sample; the fan is back at full
    # speed when Kp e + I reaches 0
    n_expected = (pi.I - pi.Kp * 5.0) / (pi.Kp * pi.Ts / pi.Ti * 5.0)
    recover = [c.step(60.0 * (300 + k), y(RH=0.50))[1] for k in range(int(n_expected) + 2)]
    assert recover[int(n_expected) - 2] < 1.0 and recover[-1] == 1.0


def test_h1_freezes_while_the_compressor_is_off():
    c = FanHumidityPI(FixedCompressor(0.0), RH_set=0.55)
    for k in range(20):
        c.step(60.0 * k, y(RH=0.70))
    assert c.pi.I == 0.0


# ---- H2 --------------------------------------------------------------------------------------

def test_h2_inner_loop_moves_the_fan_the_right_way():
    c = CoilCascade(FixedCompressor(0.5), RH_set=0.55, T_s_max=20.0)
    c.outer.I = 6.0                     # with RH at its set point the coil target stays at 14 degC
    c.f = 0.7
    for k in range(5):                                    # let the start-up settle
        c.step(60.0 * k, y(RH=0.55, T_coil=14.0))
    assert c.T_s_set == pytest.approx(14.0) and c.f == pytest.approx(0.7)
    f0 = c.f
    c.step(300.0, y(RH=0.55, T_coil=17.0))                # coil too warm: slow the fan
    assert c.f < f0
    f1 = c.f
    c.step(360.0, y(RH=0.55, T_coil=5.0))                 # coil too cold: speed it up
    assert c.f > f1


def test_h2_outer_loop_lowers_the_coil_set_point_when_humid():
    c = CoilCascade(FixedCompressor(0.5), RH_set=0.55, T_s_max=20.0, T_s_min=8.0)
    sets = [c.step(60.0 * k, y(RH=0.65, T_coil=15.0)) is not None and c.T_s_set for k in range(400)]
    assert sets[0] < 20.0 and sets[-1] == pytest.approx(8.0)
    assert np.all(np.diff(sets) <= 1e-12)


def test_h2_tracks_the_coil_set_point_on_the_real_coil():
    """Inner loop against the actual coil model at fixed compressor speed and return air."""
    c = CoilCascade(FixedCompressor(0.5), RH_set=0.55)
    T_a, W = 25.0, 0.0120
    st = UNIT.operate(0.5, 1.0, T_a, W, 32.0)
    for k in range(40):
        c.outer.I = 7.0                                   # RH at its set point: target 20 - 7 = 13
        cmd = c.step(60.0 * k, y(T_a=T_a, RH=0.55, T_coil=st.T_s))
        st = UNIT.operate(cmd[0], cmd[1], T_a, W, 32.0)
    assert c.T_s_set == pytest.approx(13.0)
    # the target is reachable: between the coil temperatures at minimum and at full fan
    T_s_min_fan = UNIT.operate(0.5, F_MIN, T_a, W, 32.0).T_s
    assert T_s_min_fan < 13.0 < UNIT.operate(0.5, 1.0, T_a, W, 32.0).T_s
    assert st.T_s == pytest.approx(c.T_s_set, abs=0.05)


# ---- closed loop -----------------------------------------------------------------------------

S = SummerScenario()
X0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
VENT = VentilationSchedule()


@pytest.fixture(scope="module")
def loops():
    def temp():
        return make_pi_controller("sigma-delta", P, UNIT, ventilation=VENT)

    ctrls = {"pi": temp(), "h1": FanHumidityPI(temp()), "h2": CoilCascade(temp()),
             "onoff": OnOffThermostat(ventilation=VENT)}
    return {k: simulate(P, S, X0, 5 * DAY, controller=c, unit=UNIT).window(3 * DAY)
            for k, c in ctrls.items()}


def test_humidity_controllers_hold_the_rh_set_point(loops):
    for name in ("h1", "h2"):
        m = compute(loops[name], 25.0)
        assert m.RH_mean_occ == pytest.approx(55.0, abs=RH_ACCURACY), name
        assert m.h_RH60_occ == 0.0, name
        assert m.Kh_outside_occ < KH_PER_DAY * 2, name


def test_dehumidification_costs_energy_but_less_than_on_off(loops):
    m = {k: compute(r, 25.0) for k, r in loops.items()}
    for name in ("h1", "h2"):
        assert m["pi"].energy_kWh < m[name].energy_kWh < m["onoff"].energy_kWh, name


def test_cascade_holds_humidity_more_steadily(loops):
    def rh_std(r):
        occ = np.asarray(r.d.n_occ) > 0
        return (100 * r.RH[occ]).std()
    assert rh_std(loops["h2"]) < rh_std(loops["h1"])


def test_fan_commands_within_limits_and_balances_close(loops):
    for name in ("h1", "h2"):
        r = loops[name]
        s, f = r.unit.cmd[:, 0], r.unit.cmd[:, 1]
        assert np.all(f[s > 0] >= F_MIN - 1e-12) and np.all(f <= 1.0)
        assert np.all(f[s == 0] == 0.0)                  # fan follows the compressor ('auto')
        assert np.abs(r.system_water_residual).max() < 1e-12
        assert math.isfinite(r.unit.E_el[-1])
