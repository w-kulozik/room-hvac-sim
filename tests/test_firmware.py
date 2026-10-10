"""The PI with sigma-delta handling as fixed-point C firmware, in the loop (decision D59).

The tests need a C compiler ($CC, or the zig toolchain from PyPI); they are skipped if the
firmware cannot be built.
"""

import ctypes
from fractions import Fraction

import numpy as np
import pytest

from roomsim.ac_unit import ACUnit
from roomsim.controllers import VentilationSchedule, make_pi_controller
from roomsim.controllers.firmware import (
    Q16,
    FirmwareController,
    FwState,
    firmware_params,
    load_library,
    to_q16,
)
from roomsim.inputs import DAY, SummerScenario
from roomsim.metrics import compute
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.simulate import simulate

try:
    LIB = load_library()
except Exception as exc:  # any build or load failure skips the module
    pytest.skip(f"firmware not available: {exc}", allow_module_level=True)

P, UNIT = RoomParams(), ACUnit()
RNG = np.random.default_rng(20261025)


def pi_controller():
    return make_pi_controller("sigma-delta", P, UNIT, ventilation=VentilationSchedule())


def round_half_away(x: Fraction) -> int:
    n = abs(x)
    r = int(n + Fraction(1, 2))
    return r if x >= 0 else -r


def test_rounding_primitives_against_exact_arithmetic():
    for _ in range(2000):
        x = int(RNG.integers(-2**52, 2**52))
        k = int(RNG.integers(1, 40))
        assert LIB.fw_rshift_round(x, k) == round_half_away(Fraction(x, 2**k))
    for _ in range(2000):
        g, v = int(RNG.integers(-2**30, 2**30)), int(RNG.integers(-2**24, 2**24))
        assert LIB.fw_mul_q24(g, v) == round_half_away(Fraction(g * v, 2**24))
    # ties go away from zero
    assert LIB.fw_rshift_round(3, 1) == 2 and LIB.fw_rshift_round(-3, 1) == -2


def test_parameters_are_the_python_tuning():
    c = pi_controller()
    fp = firmware_params(c)
    assert fp.kp_q24 / 2**24 == pytest.approx(c.pi.Kp, abs=2**-25)
    assert fp.ki_q24 / 2**24 == pytest.approx(c.pi.Kp * c.pi.Ts / c.pi.Ti, abs=2**-25)
    assert fp.threshold_q16 / Q16 == pytest.approx(c.lowload.threshold, abs=2**-17)
    assert (fp.ts_s, fp.min_on_s, fp.min_off_s) == (60, 180, 180)
    with pytest.raises(ValueError, match="cooling PI"):
        firmware_params(make_pi_controller("hysteresis", P, UNIT))


def _closed_loop(controller, days=6):
    S = SummerScenario()
    x0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
    return simulate(P, S, x0, days * DAY, controller=controller, unit=UNIT)


def test_open_loop_decisions_match_the_python_controller():
    # the measurement sequence of a Python closed-loop run, replayed into both controllers
    r = _closed_loop(pi_controller(), days=2)
    py, fw = pi_controller(), FirmwareController(pi_controller(), LIB)
    s_py, s_fw, v_py, v_fw = [], [], [], []
    for k in range(0, len(r.t) - 1):
        y = np.array([r.T_a[k], 0.5, 600.0, 15.0])
        s_py.append(py.step(float(r.t[k]), y)[0])
        v_py.append(py.last_demand)
        s_fw.append(fw.step(float(r.t[k]), y)[0])
        v_fw.append(fw.state.demand_q16 / Q16)
    s_py, s_fw = np.array(s_py), np.array(s_fw)
    assert np.array_equal(s_py > 0, s_fw > 0)                # every on/off decision
    assert np.abs(s_py - s_fw).max() < 1e-4
    assert np.abs(np.array(v_py) - np.array(v_fw)).max() < 1e-4
    assert fw.state.starts == py.guard.starts


def test_closed_loop_matches_the_python_controller():
    days, warm = 6, 3
    r_py = _closed_loop(pi_controller(), days).window(warm * DAY)
    r_fw = _closed_loop(FirmwareController(pi_controller(), LIB), days).window(warm * DAY)
    m_py, m_fw = compute(r_py, 25.0), compute(r_fw, 25.0)
    assert abs(m_fw.energy_kWh / m_py.energy_kWh - 1.0) < 0.005
    assert abs(m_fw.starts / m_py.starts - 1.0) < 0.03
    assert abs(m_fw.T_rms_occ - m_py.T_rms_occ) < 0.01


def test_state_layout_and_reset():
    assert LIB.fw_size_of_state() == ctypes.sizeof(FwState)
    fw = FirmwareController(pi_controller(), LIB)
    assert fw.state.on == 0 and fw.state.integrator_q32 == 0 and fw.Ts == 60.0
    # a warm room starts the compressor at once; a cold one keeps it off
    hot = fw.step(0.0, np.array([30.0, 0.5, 600.0, 15.0]))
    assert hot[0] > 0 and fw.state.on == 1
    cold = FirmwareController(pi_controller(), LIB).step(0.0, np.array([20.0, 0.5, 600.0, 15.0]))
    assert cold[0] == 0.0
    assert to_q16(25.0) == 25 * Q16
