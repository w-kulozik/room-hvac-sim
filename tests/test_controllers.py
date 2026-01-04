"""Discrete controllers: PI with anti-windup, compressor guard, low-load handlers, closed loop."""

import numpy as np
import pytest

from roomsim.ac_unit import ACUnit
from roomsim.controllers import (CompressorGuard, DiscretePI, HysteresisLowLoad,
                                 OnOffThermostat, SigmaDeltaLowLoad, VentilationSchedule,
                                 fan_command, make_pi_controller, room_pi_tuning, simc_pi)
from roomsim.inputs import DAY, HOUR, ConstantInputs, SummerScenario
from roomsim.metrics import compute
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.simulate import IDX_CMD_FAN, IDX_CMD_S, IDX_CMD_VENT, simulate

from criteria import KH_PER_DAY

P = RoomParams()
UNIT = ACUnit()
S_MIN = UNIT.p.s_min


# ---- PI ----------------------------------------------------------------------------------

def test_simc_rules():
    Kp, Ti = simc_pi(K=2.0, tau=100.0, theta=10.0, tau_c=40.0)
    assert Kp == pytest.approx(100.0 / (2.0 * 50.0))
    assert Ti == pytest.approx(100.0)
    assert simc_pi(2.0, 1000.0, 10.0, 40.0)[1] == pytest.approx(200.0)    # 4(tau_c + theta)


def test_room_tuning_values():
    Kp, Ti = room_pi_tuning(P, UNIT)
    assert Ti == pytest.approx(P.tau_air)
    assert 0.08 < Kp < 0.25          # demand per K: full capacity for an error of 4-12 K


def test_pi_integrates_error_when_unsaturated():
    pi = DiscretePI(Kp=0.2, Ti=600.0, Ts=60.0)
    for _ in range(10):
        v = pi.demand(1.0)
        pi.update(1.0, v, v)                      # actuator follows exactly: no tracking term
    assert pi.I == pytest.approx(10 * 0.2 * 60.0 / 600.0)


def _saturate_then_reverse(Tt, e_sat=5.0, e_rev=-0.5, n_sat=120):
    """Large error with the actuator stuck at 1, then the error reverses.

    Returns the integrator at the end of saturation and the samples until the demand reaches 0.
    """
    pi = DiscretePI(Kp=0.2, Ti=600.0, Ts=60.0, Tt=Tt)
    for _ in range(n_sat):
        v = pi.demand(e_sat)
        pi.update(e_sat, v, min(max(v, 0.0), 1.0))
    I_sat, n = pi.I, 0
    while pi.demand(e_rev) > 0.0:
        v = pi.demand(e_rev)
        pi.update(e_rev, v, min(max(v, 0.0), 1.0))
        n += 1
    return I_sat, n


def test_back_calculation_settles_at_the_analytical_value():
    Kp, Ti, e = 0.2, 600.0, 5.0
    for Tt in (150.0, 300.0, 600.0):
        I_sat, _ = _saturate_then_reverse(Tt, e_sat=e, n_sat=2000)
        assert I_sat == pytest.approx(1.0 - Kp * e * (1.0 - Tt / Ti), abs=1e-9)


def test_back_calculation_limits_windup():
    _, n_aw = _saturate_then_reverse(Tt=None)
    _, n_none = _saturate_then_reverse(Tt=np.inf)
    assert n_none > 20 * n_aw


def test_frozen_integrator():
    pi = DiscretePI(Kp=0.2, Ti=600.0, Ts=60.0, I=0.3)
    pi.update(2.0, pi.demand(2.0), 0.0, integrate=False)
    assert pi.I == 0.3


# ---- guard, fan, ventilation ---------------------------------------------------------------

def test_compressor_guard_minimum_times():
    g = CompressorGuard(min_on=180.0, min_off=300.0)
    assert g.can_start(0.0)
    g.set(0.0, True)
    assert not g.can_stop(120.0) and g.can_stop(180.0)
    g.set(200.0, False)
    assert not g.can_start(400.0) and g.can_start(500.0)
    g.set(500.0, True)
    assert g.starts == 2


def test_fan_modes_and_ventilation_schedule():
    assert fan_command(0.0, 0.8, "auto") == 0.0
    assert fan_command(0.5, 0.8, "auto") == 0.8
    assert fan_command(0.0, 0.8, "continuous") == 0.8
    with pytest.raises(ValueError):
        fan_command(0.5, 0.8, "sometimes")
    V = VentilationSchedule()
    assert V(9 * HOUR) == V.V_occupied and V(20 * HOUR) == 0.0


# ---- low-load handlers with synthetic demand ------------------------------------------------

def test_hysteresis_handler_cycles_and_freezes():
    h = HysteresisLowLoad(S_MIN, CompressorGuard(0.0, 0.0), dT_on=0.5, dT_off=0.5)
    assert h.step(0.0, 0.6, 26.0, 25.0) == (0.6, 0.6, True)          # modulating
    s, u, integrate = h.step(60.0, 0.1, 25.2, 25.0)                     # low load, still warm
    assert s == S_MIN and not integrate
    assert h.step(120.0, 0.1, 24.4, 25.0)[0] == 0.0                     # stops below 24.5
    assert h.step(180.0, 0.1, 25.2, 25.0)[0] == 0.0                     # waits until 25.5
    assert h.step(240.0, 0.1, 25.6, 25.0)[0] == S_MIN                   # restarts


@pytest.mark.parametrize("v", [0.05, 0.12, 0.2, 0.27])
def test_sigma_delta_meets_the_average_demand(v):
    Ts = 60.0
    sd = SigmaDeltaLowLoad(S_MIN, CompressorGuard(180.0, 180.0), Ts, T_ref=600.0)
    s = np.array([sd.step(k * Ts, v, 25.0, 25.0)[0] for k in range(int(24 * HOUR / Ts))])
    assert s.mean() == pytest.approx(v, abs=0.01)
    assert set(np.unique(s)) <= {0.0, S_MIN}


def test_sigma_delta_cycle_length_follows_the_design_formula():
    Ts, T_ref, v = 60.0, 600.0, 0.5 * S_MIN
    sd = SigmaDeltaLowLoad(S_MIN, CompressorGuard(0.0, 0.0), Ts, T_ref)
    for k in range(int(24 * HOUR / Ts)):
        sd.step(k * Ts, v, 25.0, 25.0)
    cycles_per_h = sd.guard.starts / 24.0
    assert cycles_per_h == pytest.approx(HOUR / (4 * T_ref), rel=0.1)


def test_sigma_delta_respects_minimum_times():
    Ts = 60.0
    sd = SigmaDeltaLowLoad(S_MIN, CompressorGuard(300.0, 300.0), Ts, T_ref=60.0)  # wants short cycles
    s = np.array([sd.step(k * Ts, 0.15, 25.0, 25.0)[0] for k in range(600)])
    on = s > 0
    edges = np.flatnonzero(np.diff(on.astype(int))) + 1
    runs = np.diff(edges)
    assert runs.min() * Ts >= 300.0


# ---- closed loop ---------------------------------------------------------------------------

S = SummerScenario()
X0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
VENT = VentilationSchedule()


def run(ctrl, days=4, keep=2):
    r = simulate(P, S, X0, days * DAY, controller=ctrl, unit=UNIT)
    return r.window((days - keep) * DAY)


@pytest.fixture(scope="module")
def loops():
    return {
        "onoff": run(OnOffThermostat(ventilation=VENT)),
        "pi_A": run(make_pi_controller("hysteresis", P, UNIT, ventilation=VENT)),
        "pi_C": run(make_pi_controller("sigma-delta", P, UNIT, ventilation=VENT)),
    }


def test_all_controllers_hold_the_set_point(loops):
    for name, r in loops.items():
        m = compute(r, 25.0)
        assert m.Kh_outside_occ < KH_PER_DAY * 2, name


def test_pi_tracks_better_and_starts_less_than_on_off(loops):
    m = {k: compute(r, 25.0) for k, r in loops.items()}
    assert m["pi_A"].T_rms_occ < m["onoff"].T_rms_occ
    assert m["pi_C"].starts < m["onoff"].starts
    assert m["pi_C"].starts_per_h_on < 0.5 * m["onoff"].starts_per_h_on   # gentler on the compressor
    assert m["pi_A"].energy_kWh < m["onoff"].energy_kWh


def test_commands_are_feasible(loops):
    for name, r in loops.items():
        s, f = r.cmd[:, IDX_CMD_S], r.cmd[:, IDX_CMD_FAN]
        assert np.all((s == 0) | ((s >= S_MIN - 1e-12) & (s <= 1.0))), name
        assert np.all((f == 0) == (s == 0)), name                       # fan 'auto'
        occupied = np.asarray(r.d.n_occ) > 0
        assert np.all(r.cmd[occupied, IDX_CMD_VENT] == VENT.V_occupied), name


def test_balances_close_in_closed_loop(loops):
    for name, r in loops.items():
        assert np.abs(r.energy_residual).max() < 1e-11 * (P.C_a + P.C_m) * 40, name
        assert np.abs(r.system_water_residual).max() < 1e-12, name


def test_minimum_times_hold_in_closed_loop(loops):
    for name, r in loops.items():
        on = r.cmd[:-1, IDX_CMD_S] > 0
        edges = np.flatnonzero(np.diff(on.astype(int))) + 1
        if len(edges) > 2:
            assert np.diff(edges).min() * r.dt >= 180.0, name


def test_pi_rejects_a_load_step_without_offset():
    d = ConstantInputs(T_o=32.0, W_o=0.0116, G_sol=300.0, n_occ=2)
    ctrl = make_pi_controller("sigma-delta", P, UNIT)
    r = simulate(P, d, equilibrium_state(25.0, 27.0, 0.0116), 12 * HOUR, controller=ctrl, unit=UNIT)
    w = r.window(8 * HOUR)
    assert abs(w.T_a.mean() - 25.0) < 0.05
    assert np.all(w.cmd[:, IDX_CMD_S] > S_MIN)            # high load: continuous modulation
