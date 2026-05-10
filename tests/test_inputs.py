"""Disturbance scenarios reproduce the specified daily profiles."""

import numpy as np
import pytest

from roomsim import psychro
from roomsim.inputs import DAY, HOUR, ConstantInputs, SummerScenario

S = SummerScenario()
t_day = np.arange(0.0, DAY, 60.0)


def test_outdoor_temperature_extremes_and_timing():
    T = S.T_o(t_day)
    assert T.max() == pytest.approx(32.0, abs=1e-6)
    assert T.min() == pytest.approx(20.0, abs=1e-6)
    assert t_day[T.argmax()] / HOUR == pytest.approx(15.0)
    assert t_day[T.argmin()] / HOUR == pytest.approx(3.0)
    assert T.mean() == pytest.approx(26.0, abs=1e-9)


def test_default_outdoor_humidity_ratio_is_constant():
    t = np.linspace(0, DAY, 97)
    W = S.W_o(t)
    assert np.all(W == W[0])
    assert W[0] == pytest.approx(psychro.W_from_RH(26.0, 0.55))
    RH = psychro.RH_from_W(S.T_o(t), W)
    assert RH.max() == pytest.approx(0.55 * psychro.p_ws(26.0) / psychro.p_ws(20.0))


def test_constant_RH_option():
    t = np.linspace(0, DAY, 97)
    S_rh = SummerScenario(constant_W=False)
    assert psychro.RH_from_W(S_rh.T_o(t), S_rh.W_o(t)) == pytest.approx(0.55, rel=1e-12)


def test_solar_half_sine():
    G = S.G_sol(t_day)
    h = t_day / HOUR
    assert np.all(G[(h < 6) | (h >= 18)] == 0.0)
    assert G.max() == pytest.approx(500.0)
    assert h[G.argmax()] == pytest.approx(12.0)
    # Daily irradiation: integral of a half-sine = 2/pi * peak * duration
    daily = np.trapezoid(G, t_day)
    assert daily == pytest.approx(2 / np.pi * 500.0 * 12 * HOUR, rel=1e-4)


def test_occupancy_window_is_right_open():
    assert S.n_occ(8 * HOUR - 1) == 0
    assert S.n_occ(8 * HOUR) == 2
    assert S.n_occ(18 * HOUR - 1) == 2
    assert S.n_occ(18 * HOUR) == 0


def test_periodic_and_negative_time():
    t = np.array([-DAY + 15 * HOUR, 15 * HOUR, 2 * DAY + 15 * HOUR])
    d = S(t)
    for field in d:
        assert np.allclose(field, field[1])


def test_scalar_and_vector_calls_agree():
    t = np.array([0.0, 10 * HOUR, 13.5 * HOUR])
    vec = S(t)
    for i, ti in enumerate(t):
        sca = S(ti)
        for a, b in zip(sca, vec):
            assert float(a) == pytest.approx(b[i])


def test_constant_inputs_shape():
    d = ConstantInputs(T_o=30.0)(np.zeros(5))
    assert d.T_o.shape == (5,)
    assert np.all(d.T_o == 30.0)


def test_variable_occupancy_profile():
    from roomsim.inputs import VARIABLE_OCCUPANCY
    S_var = SummerScenario(occupancy=VARIABLE_OCCUPANCY)
    hours = np.array([7.5, 8.5, 10.0, 12.5, 14.0, 15.5, 17.0, 18.5])
    assert list(S_var.n_occ(hours * HOUR)) == [0, 1, 2, 0, 2, 3, 1, 0]
    t = np.arange(0.0, DAY, 60.0)
    assert S_var.n_occ(t).mean() * 24 / 10 == pytest.approx(1.6)     # mean over 8-18 h
