"""Air-conditioning unit coupled to the room through the simulation loop."""

import numpy as np
import pytest

from roomsim import psychro
from roomsim.ac_unit import ACUnit
from roomsim.controllers import make_pi_controller
from roomsim.inputs import HOUR, ConstantInputs, SummerScenario
from roomsim.params import RoomParams
from roomsim.plant import IDX_TA, IDX_W, equilibrium_state
from roomsim.simulate import simulate, unit_input
from roomsim.verification import convergence_order

P = RoomParams()
UNIT = ACUnit()
HOT = ConstantInputs(T_o=32.0, W_o=0.0140, G_sol=200.0, n_occ=2)
X_HOT = equilibrium_state(30.0, 31.0, 0.0135, C=900.0)


class Fixed:
    """Constant command [s, f, V_vent]."""

    def __init__(self, s, f=1.0, vent=0.0):
        self.cmd = np.array([s, f, vent])

    def step(self, t, y):
        return self.cmd


def test_unit_cools_and_dries_the_room():
    r = simulate(P, HOT, X_HOT, 6 * HOUR, controller=Fixed(0.8), unit=UNIT,
                 unit_running_at_start=True)
    assert r.T_a[-1] < r.T_a[0] - 4.0
    assert r.W[-1] < r.W[0]
    assert np.all(r.unit.op["Q_sens"] > 0) and np.all(r.unit.op["m_lat"] >= 0)
    assert r.unit.E_el[-1] > 0


def test_balances_close_with_the_unit_running():
    r = simulate(P, SummerScenario(), X_HOT, 24 * HOUR, controller=Fixed(0.6, 0.7, 0.01), unit=UNIT)
    E_abs = P.C_a * 40 + P.C_m * 40
    assert np.abs(r.energy_residual).max() < 1e-11 * E_abs
    assert np.abs(r.water_residual).max() < 1e-11
    assert np.abs(r.co2_residual).max() < 1e-13
    # the unit term in the balance is what the coil and fan delivered
    k = 500
    assert r.q["q_hvac"][k] == pytest.approx(-r.unit.op["Q_sens"][k] + UNIT.fan_power(0.7))
    assert r.q["m_hvac"][k] == pytest.approx(-r.unit.op["m_lat"][k])


def test_electrical_energy_matches_logged_power():
    r = simulate(P, HOT, X_HOT, 3 * HOUR, controller=Fixed(0.7), unit=UNIT,
                 unit_running_at_start=True)
    trap = np.sum(0.5 * (r.unit.op["P_el"][1:] + r.unit.op["P_el"][:-1])) * r.dt
    assert r.unit.E_el[-1] == pytest.approx(trap, rel=1e-5)


def test_start_up_ramp_after_a_start():
    r = simulate(P, HOT, X_HOT, 1 * HOUR, dt=10.0, controller=Fixed(1.0), unit=UNIT)
    assert r.unit.op["Q_coil"][0] == 0.0
    for k in (3, int(UNIT.p.tau_start / r.dt), 30):
        ramp = 1.0 - np.exp(-k * r.dt / UNIT.p.tau_start)
        # delivered capacity: the unit at the logged air state with the analytical ramp
        st = UNIT.operate(1.0, 1.0, r.T_a[k], r.W[k], HOT.T_o, ramp=ramp, M_c=r.unit.M_coil[k])
        assert r.unit.op["Q_coil"][k] == pytest.approx(st.Q_coil, rel=1e-12)
        # power: that of the full requested capacity at the same air state
        full = UNIT.operate(1.0, 1.0, r.T_a[k], r.W[k], HOT.T_o)
        assert r.unit.op["P_el"][k] == pytest.approx(full.P_el, rel=1e-12)
        assert st.Q_coil < full.Q_coil


# ---- accuracy of the coupling ---------------------------------------------------------

class HeldOutput:
    """Evaluates the unit once per sample and holds Q_sens and m_lat (the naive coupling)."""

    def __init__(self, s):
        self.s = s

    def step(self, t, y):
        T_a, RH = y[0], y[1]
        W = float(psychro.W_from_RH(T_a, RH))
        x = np.zeros(5); x[IDX_TA] = T_a; x[IDX_W] = W
        u, _ = unit_input(UNIT, np.array([self.s, 1.0, 0.0]), 1.0, x, HOT(t))
        return u


DTS = np.array([240.0, 120.0, 60.0, 30.0])
# dry outdoor air: the coil stays dry, so the right-hand side is smooth along the trajectory
DRY = ConstantInputs(T_o=32.0, W_o=0.0070, G_sol=200.0, n_occ=2)
X_DRY = equilibrium_state(30.0, 31.0, 0.0070, C=900.0)


def _final_error(dts, run, dt_ref=7.5):
    ref = run(dt_ref)
    return [abs(run(dt)[IDX_TA] - ref[IDX_TA]) for dt in dts]


def _wet_fractions(r, T_o):
    return np.array([UNIT.operate(c[0], c[1], r.T_a[k], r.W[k], T_o).wet_fraction
                     for k, c in enumerate(r.unit.cmd) if c[0] > 0])


def test_unit_inside_rk4_stages_keeps_fourth_order():
    def sim(dt):
        return simulate(P, DRY, X_DRY, 4 * HOUR, dt=dt, controller=Fixed(0.5), dt_ctrl=dt,
                        unit=UNIT, unit_running_at_start=True)
    assert np.all(_wet_fractions(sim(60.0), DRY.T_o) == 0.0)       # one coil regime throughout
    order = convergence_order(DTS, _final_error(DTS, lambda dt: sim(dt).x[-1]))
    assert order == pytest.approx(4.0, abs=0.3)


def test_crossing_into_the_partially_wet_regime_costs_order():
    """Where the wetted fraction leaves 1 the coil output is continuous but its derivative is not
    (a kink in the right-hand side). A fixed-step method then loses order on the step that
    contains the crossing: second order is expected and observed. The error stays small."""
    def sim(dt):
        return simulate(P, HOT, X_HOT, 4 * HOUR, dt=dt, controller=Fixed(0.8), dt_ctrl=dt,
                        unit=UNIT, unit_running_at_start=True)
    wet = _wet_fractions(sim(60.0), HOT.T_o)
    assert wet.max() == 1.0 and wet.min() < 0.9                    # the crossing happens
    dts = np.append(DTS, 15.0)
    err = _final_error(dts, lambda dt: sim(dt).x[-1], dt_ref=3.75)
    assert convergence_order(dts, err) == pytest.approx(2.0, abs=0.3)
    assert err[2] < 1e-8                                            # K, at the 60 s step


def test_closed_loop_with_sampled_pi_is_fourth_order():
    """The real closed loop: PI with sigma-delta low-load handling, sampled every 60 s, while the
    integration step is refined. Between samples the command is held, so the closed loop is a
    sampled-data system that RK4 solves to fourth order as long as no switching decision changes
    with the step (the compressor modulates throughout) and the coil stays in one regime."""
    def sim(dt):
        ctrl = make_pi_controller("sigma-delta", P, UNIT, Ts=60.0)
        return simulate(P, DRY, X_DRY, 6 * HOUR, dt=dt, controller=ctrl, dt_ctrl=60.0, unit=UNIT,
                        unit_running_at_start=True)
    r = sim(60.0)
    s = r.unit.cmd[:, 0]
    assert np.all(s > UNIT.p.s_min) and s.max() - s.min() > 0.2      # modulating, not switching
    assert np.all(_wet_fractions(r, DRY.T_o) == 0.0)
    dts = np.array([60.0, 30.0, 15.0])
    err = _final_error(dts, lambda dt: sim(dt).x[-1], dt_ref=3.75)
    assert convergence_order(dts, err) == pytest.approx(4.0, abs=0.3)


def test_holding_the_unit_output_is_only_first_order():
    def run(dt):
        return simulate(P, HOT, X_HOT, 4 * HOUR, dt=dt, controller=HeldOutput(0.8), dt_ctrl=dt).x[-1]
    order = convergence_order(DTS, _final_error(DTS, run))
    assert order == pytest.approx(1.0, abs=0.15)
