"""Air-conditioning unit coupled to the room through the simulation loop."""

import numpy as np
import pytest

from roomsim import psychro
from roomsim.ac_unit import ACUnit
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
    assert np.all(r.unit["Q_sens"] > 0) and np.all(r.unit["m_lat"] >= 0)
    assert r.E_el[-1] > 0


def test_balances_close_with_the_unit_running():
    r = simulate(P, SummerScenario(), X_HOT, 24 * HOUR, controller=Fixed(0.6, 0.7, 0.01), unit=UNIT)
    E_abs = P.C_a * 40 + P.C_m * 40
    assert np.abs(r.energy_residual).max() < 1e-11 * E_abs
    assert np.abs(r.water_residual).max() < 1e-11
    assert np.abs(r.co2_residual).max() < 1e-13
    # the unit term in the balance is what the coil and fan delivered
    k = 500
    assert r.q["q_hvac"][k] == pytest.approx(-r.unit["Q_sens"][k] + UNIT.p.P_fan_in_nom * 0.7**3)
    assert r.q["m_hvac"][k] == pytest.approx(-r.unit["m_lat"][k])


def test_electrical_energy_matches_logged_power():
    r = simulate(P, HOT, X_HOT, 3 * HOUR, controller=Fixed(0.7), unit=UNIT,
                 unit_running_at_start=True)
    trap = np.sum(0.5 * (r.unit["P_el"][1:] + r.unit["P_el"][:-1])) * r.dt
    assert r.E_el[-1] == pytest.approx(trap, rel=1e-5)


def test_start_up_ramp_after_a_start():
    r = simulate(P, HOT, X_HOT, 1 * HOUR, dt=10.0, controller=Fixed(1.0), unit=UNIT)
    Q = r.unit["Q_coil"]
    k_tau = int(UNIT.p.tau_start / r.dt)
    assert Q[0] == 0.0
    assert Q[k_tau] == pytest.approx((1 - np.exp(-1)) * Q[10 * k_tau], rel=0.05)
    # during the ramp the compressor draws the power of the full capacity at the same air state
    k = 3
    full = UNIT.operate(1.0, 1.0, r.T_a[k], r.W[k], HOT.T_o)
    assert r.unit["P_el"][k] == pytest.approx(full.P_el, rel=1e-12)
    assert r.unit["Q_coil"][k] < 0.25 * full.Q_coil


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


def _final_error(dts, run):
    ref = run(7.5)
    return [abs(run(dt)[IDX_TA] - ref[IDX_TA]) for dt in dts]


def test_unit_inside_rk4_stages_keeps_fourth_order():
    def run(dt):
        return simulate(P, HOT, X_HOT, 4 * HOUR, dt=dt, controller=Fixed(0.8), dt_ctrl=dt,
                        unit=UNIT, unit_running_at_start=True).x[-1]
    order = convergence_order(DTS, _final_error(DTS, run))
    assert order == pytest.approx(4.0, abs=0.3)


def test_holding_the_unit_output_is_only_first_order():
    def run(dt):
        return simulate(P, HOT, X_HOT, 4 * HOUR, dt=dt, controller=HeldOutput(0.8), dt_ctrl=dt).x[-1]
    order = convergence_order(DTS, _final_error(DTS, run))
    assert order == pytest.approx(1.0, abs=0.15)
