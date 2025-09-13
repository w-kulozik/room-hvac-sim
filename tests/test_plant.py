"""Basic structural checks of the plant equations (signs, equilibria, bookkeeping)."""

import numpy as np
import pytest

from roomsim.inputs import ConstantInputs, SummerScenario
from roomsim.params import RoomParams
from roomsim.plant import IDX_TA, IDX_TM, IDX_W, U_OFF, fluxes, rhs

P = RoomParams()


def test_uniform_equilibrium_without_gains():
    d = ConstantInputs(T_o=26.0, W_o=0.011)(0.0)
    x = np.array([26.0, 26.0, 0.011])
    assert np.allclose(rhs(0.0, x, U_OFF, d, P), 0.0)


def test_warmer_outside_heats_both_nodes():
    d = ConstantInputs(T_o=30.0, W_o=0.011)(0.0)
    dx = rhs(0.0, np.array([20.0, 20.0, 0.011]), U_OFF, d, P)
    assert dx[IDX_TA] > 0 and dx[IDX_TM] > 0


def test_cooling_and_dehumidification_signs():
    d = ConstantInputs(T_o=26.0, W_o=0.011)(0.0)
    x = np.array([26.0, 26.0, 0.011])
    dx = rhs(0.0, x, np.array([-1000.0, 1e-4]), d, P)
    assert dx[IDX_TA] == pytest.approx(-1000.0 / P.C_a)
    assert dx[IDX_TM] == 0.0
    assert dx[IDX_W] == pytest.approx(-1e-4 / P.m_air)


def test_air_mass_exchange_cancels_in_boundary_total():
    S = SummerScenario()
    x = np.array([27.0, 25.0, 0.012])
    f = fluxes(0.0, x, U_OFF, S(13 * 3600.0), P)
    assert f.air_total + f.mass_total == pytest.approx(f.boundary_total, rel=1e-14)


def test_gains_split_conserves_totals():
    d = ConstantInputs(G_sol=400.0, n_occ=2)(0.0)
    f = fluxes(0.0, np.array([26.0, 26.0, 0.01]), U_OFF, d, P)
    assert f.q_solar_air + f.q_solar_mass == pytest.approx(P.g_win * P.A_win * 400.0)
    assert f.q_internal_conv + f.q_internal_rad == pytest.approx(2 * 70.0)
