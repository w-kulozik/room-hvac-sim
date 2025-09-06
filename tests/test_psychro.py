"""Psychrometric reference values and internal consistency."""

import numpy as np
import pytest

from roomsim import psychro as ps


def test_saturation_pressure_20C():
    # Reference tables give 2.339 kPa; Magnus is within 0.3 %.
    assert ps.p_ws(20.0) == pytest.approx(2333.0, rel=1e-3)
    assert ps.p_ws(20.0) == pytest.approx(2339.0, rel=3e-3)


def test_humidity_ratio_25C_50pct():
    assert ps.W_from_RH(25.0, 0.50) == pytest.approx(0.0099, abs=1e-4)


def test_dew_point_20C_50pct():
    assert ps.dew_point(20.0, 0.50) == pytest.approx(9.3, abs=0.05)


def test_pw_W_round_trip():
    p_w = np.linspace(100.0, 7000.0, 50)
    assert ps.pw_from_W(ps.W_from_pw(p_w)) == pytest.approx(p_w, rel=1e-12)


def test_RH_W_round_trip():
    T = np.linspace(-10.0, 45.0, 12)
    RH = np.linspace(0.05, 1.0, 12)
    assert ps.RH_from_W(T, ps.W_from_RH(T, RH)) == pytest.approx(RH, rel=1e-12)


def test_dew_point_at_saturation_equals_dry_bulb():
    T = np.linspace(-20.0, 45.0, 14)
    assert ps.dew_point(T, 1.0) == pytest.approx(T, abs=1e-10)


def test_dew_point_from_W_consistent():
    T, RH = 26.0, 0.55
    assert ps.dew_point_from_W(ps.W_from_RH(T, RH)) == pytest.approx(ps.dew_point(T, RH), abs=1e-10)


def test_enthalpy_dry_air_and_latent_part():
    assert ps.enthalpy(20.0, 0.0) == pytest.approx(20.12)
    # 25 degC, 50 %: about 50.3 kJ/kg in standard charts
    assert ps.enthalpy(25.0, ps.W_from_RH(25.0, 0.5)) == pytest.approx(50.3, abs=0.3)
