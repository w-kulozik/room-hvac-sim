"""Derived parameters match the hand calculations in the specification."""

from dataclasses import replace

import pytest

from roomsim.params import RoomParams


def test_derived_values():
    p = RoomParams()
    assert p.volume == pytest.approx(54.0)
    assert p.C_a == pytest.approx(0.326e6, rel=1e-3)
    assert p.H_inf == pytest.approx(9.05, rel=1e-3)
    assert p.H_am == pytest.approx(180.0)
    assert p.g_vap_person == pytest.approx(0.05 / 3600)


def test_time_constants():
    p = RoomParams()
    assert p.tau_air == pytest.approx(1700.0, rel=0.02)
    assert p.tau_mass / 86400 == pytest.approx(3.1, rel=0.05)   # about three days
    assert p.tau_moisture == pytest.approx(7200.0)


def test_replace_keeps_derived_consistent():
    p = replace(RoomParams(), ach=1.0)
    assert p.H_inf == pytest.approx(2 * RoomParams().H_inf)
