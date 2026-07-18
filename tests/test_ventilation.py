"""Ventilation fan model and ventilation strategies (decisions D23-D26)."""

import numpy as np
import pytest

from criteria import CO2_EXCEEDANCE_H_PER_DAY, CO2_LIMIT
from roomsim.ac_unit import ACUnit
from roomsim.controllers import (
    CoilCascade,
    DemandControlledPI,
    NightFlush,
    ScheduledVentilation,
    StagedVentilation,
    VentilationSchedule,
    WithVentilation,
    dcv_tuning,
    make_pi_controller,
)
from roomsim.inputs import DAY, HOUR, ConstantInputs, SummerScenario
from roomsim.params import RoomParams
from roomsim.plant import IDX_V_VENT, equilibrium_state
from roomsim.simulate import simulate
from roomsim.ventilation import VentilationFan

P = RoomParams()
UNIT = ACUnit()


def y(T_a=25.0, RH=0.55, C=800.0, T_coil=15.0):
    return np.array([T_a, RH, C, T_coil])


# ---- fan model -------------------------------------------------------------------------------

def test_fan_affinity_law_and_exact_energy():
    fan = VentilationFan()
    assert fan.power(0.03) == pytest.approx(30.0)
    assert fan.power(0.015) == pytest.approx(30.0 / 8)
    r = simulate(P, ConstantInputs(), equilibrium_state(25, 25, 0.01), 2 * HOUR,
                 controller=_FixedVent(0.02))
    assert fan.energy_kWh(r) == pytest.approx(fan.power(0.02) * 2 * HOUR / 3.6e6, rel=1e-12)


class _FixedVent:
    def __init__(self, V):
        self.V = V

    def step(self, t, y):
        u = np.zeros(3)
        u[IDX_V_VENT] = self.V
        return u


# ---- demand-controlled PI ------------------------------------------------------------------

def test_dcv_tuning_follows_the_linearised_co2_balance():
    Kp, Ti = dcv_tuning(P, C_set=900.0, n_people=2.0, tau_c=900.0, Ts=60.0)
    Q = 2 * P.co2_person * 1e6 / 480.0                       # total outdoor air at 900 ppm
    assert pytest.approx(0.02167, rel=1e-3) == Q
    tau, K = P.volume / Q, 480.0 / Q
    assert Ti == pytest.approx(tau)                           # tau < 4 (tau_c + theta)
    assert Kp == pytest.approx(tau / (K * (900.0 + 30.0)))


def test_dcv_is_off_when_unoccupied_and_bounded_when_occupied():
    v = DemandControlledPI(P)
    assert v.step(3 * HOUR, y(C=1500.0)) == 0.0
    assert v.step(9 * HOUR, y(C=500.0)) == pytest.approx(v.V_min)
    for k in range(200):
        V = v.step(9 * HOUR + 60.0 * k, y(C=1600.0))
    assert pytest.approx(v.V_max) == V


def test_dcv_holds_the_co2_set_point_in_closed_loop():
    d = ConstantInputs(T_o=25.0, W_o=0.010, n_occ=2)
    ctrl = WithVentilation(_FixedVent(0.0), DemandControlledPI(P, comfort=_AlwaysOccupied()))
    r = simulate(P, d, equilibrium_state(25, 25, 0.010, C=420.0), 12 * HOUR, controller=ctrl)
    assert r.CO2[-1] == pytest.approx(900.0, abs=1.0)
    Q = 2 * P.co2_person * 1e6 / 480.0
    assert r.u[-1, IDX_V_VENT] == pytest.approx(Q - P.V_dot_inf, rel=1e-3)


class _AlwaysOccupied:
    def occupied(self, t):
        return True


# ---- staged ---------------------------------------------------------------------------------

def test_staged_levels_hysteresis_and_dwell():
    s = StagedVentilation(dwell=300.0)
    t0 = 9 * HOUR
    assert s.step(t0, y(C=800.0)) == 0.006
    assert s.step(t0 + 60, y(C=860.0)) == 0.012            # up at 850
    assert s.step(t0 + 120, y(C=990.0)) == 0.012           # dwell: no second switch yet
    assert s.step(t0 + 360, y(C=990.0)) == 0.024           # up at 950 after 5 min
    assert s.step(t0 + 720, y(C=860.0)) == 0.024           # inside the hysteresis band
    assert s.step(t0 + 780, y(C=840.0)) == 0.012           # down at 850
    assert s.step(20 * HOUR, y(C=1200.0)) == 0.0           # unoccupied: off, reset
    with pytest.raises(ValueError):
        StagedVentilation(levels=(0.01, 0.02), up=(900.0, 1000.0), down=(800.0,))


# ---- night flush ----------------------------------------------------------------------------

def test_night_flush_conditions():
    base = ScheduledVentilation()
    nf = NightFlush(base, outdoor=lambda t: 21.0)
    assert nf.step(3 * HOUR, y(T_a=25.0)) == nf.V_flush        # empty, 4 K cooler outside
    assert nf.step(3 * HOUR, y(T_a=22.5)) == 0.0               # less than 2 K difference
    assert NightFlush(base, outdoor=lambda t: 15.0).step(3 * HOUR, y(T_a=21.5)) == 0.0  # room cold
    assert nf.step(10 * HOUR, y(T_a=27.0)) == base.step(10 * HOUR, y())  # occupied: inner rules


# ---- closed loop, against the criterion fixed beforehand (D24) -----------------------------

S = SummerScenario()
X0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))


def thermal():
    return CoilCascade(make_pi_controller("sigma-delta", P, UNIT))


@pytest.fixture(scope="module")
def runs():
    strategies = {
        "fixed": ScheduledVentilation(VentilationSchedule()),
        "dcv": DemandControlledPI(P),
        "staged": StagedVentilation(),
        "flush": NightFlush(DemandControlledPI(P), outdoor=S.T_o),
    }
    return {k: simulate(P, S, X0, 4 * DAY, controller=WithVentilation(thermal(), v),
                        unit=UNIT).window(2 * DAY) for k, v in strategies.items()}


def co2_exceedance_h(r):
    occ = np.asarray(r.d.n_occ)[:-1] > 0
    return float(np.sum(r.CO2[:-1][occ] > CO2_LIMIT) * r.dt / HOUR)


@pytest.mark.parametrize("name", ["dcv", "staged", "flush"])
def test_demand_control_meets_the_co2_requirement(runs, name):
    assert co2_exceedance_h(runs[name]) <= CO2_EXCEEDANCE_H_PER_DAY * 2


def test_strategies_respect_their_flow_limits(runs):
    for name, r in runs.items():
        V = r.u[:-1, IDX_V_VENT]
        occ = np.asarray(r.d.n_occ)[:-1] > 0
        assert np.all(V >= 0.0) and np.all(V <= 0.03 + 1e-12), name
        if name != "fixed":
            assert np.all(V[occ] >= 0.006 - 1e-12), name
        if name in ("dcv", "staged"):
            assert np.all(V[~occ] == 0.0), name
