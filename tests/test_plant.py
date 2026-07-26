"""Basic structural checks of the plant equations (signs, equilibria, bookkeeping)."""

from dataclasses import replace

import numpy as np
import pytest

from roomsim import psychro
from roomsim.inputs import ConstantInputs, SummerScenario
from roomsim.params import RoomParams
from roomsim.plant import IDX_PHI_B, IDX_TA, IDX_TM, IDX_W, U_OFF, equilibrium_state, fluxes, rhs

P = RoomParams()


def test_uniform_equilibrium_without_gains():
    d = ConstantInputs(T_o=26.0, W_o=0.011)(0.0)
    x = equilibrium_state(26.0, 26.0, 0.011)
    assert np.allclose(rhs(0.0, x, U_OFF, d, P), 0.0, atol=1e-15)


def test_warmer_outside_heats_both_nodes():
    d = ConstantInputs(T_o=30.0, W_o=0.011)(0.0)
    dx = rhs(0.0, equilibrium_state(20.0, 20.0, 0.011), U_OFF, d, P)
    assert dx[IDX_TA] > 0 and dx[IDX_TM] > 0


def test_cooling_and_dehumidification_signs():
    d = ConstantInputs(T_o=26.0, W_o=0.011)(0.0)
    x = equilibrium_state(26.0, 26.0, 0.011)
    dx = rhs(0.0, x, np.array([-1000.0, 1e-4, 0.0]), d, P)
    assert dx[IDX_TA] == pytest.approx(-1000.0 / P.C_a)
    assert dx[IDX_TM] == 0.0
    assert dx[IDX_W] == pytest.approx(-1e-4 / P.m_air)


def test_air_mass_exchange_cancels_in_boundary_total():
    S = SummerScenario()
    x = np.array([27.0, 25.0, 0.012, 0.5, 800.0])
    f = fluxes(0.0, x, U_OFF, S(13 * 3600.0), P)
    assert f.air_total + f.mass_total == pytest.approx(f.boundary_total, rel=1e-14)


def test_gains_split_conserves_totals():
    d = ConstantInputs(G_sol=400.0, n_occ=2)(0.0)
    f = fluxes(0.0, equilibrium_state(26.0, 26.0, 0.01), U_OFF, d, P)
    assert f.q_solar_air + f.q_solar_mass == pytest.approx(P.g_win * P.A_win * 400.0)
    assert f.q_internal_conv + f.q_internal_rad == pytest.approx(2 * 70.0)


# ---- moisture buffer ----------------------------------------------------------------

def test_buffer_parameters_are_plausible():
    assert P.beta_surface == pytest.approx(1.8e-8, rel=0.05)
    assert P.buffer_capacity == pytest.approx(1.44)           # kg water per unit RH
    # Seen from W, the buffer adds roughly the air's own moisture capacity
    dW_dphi = psychro.W_from_RH(25.0, 1.0)
    assert P.buffer_capacity / dW_dphi / P.m_air == pytest.approx(1.1, rel=0.1)


def test_water_exchange_air_buffer_is_conservative():
    d = ConstantInputs(T_o=30.0, W_o=0.014, n_occ=2)(0.0)
    x = np.array([25.0, 24.0, 0.009, 0.7, 600.0])
    u = np.array([0.0, 2e-5, 0.01])
    f = fluxes(0.0, x, u, d, P)
    dx = rhs(0.0, x, u, d, P)
    d_water = P.m_air * dx[IDX_W] + P.buffer_capacity * dx[IDX_PHI_B]
    assert d_water == pytest.approx(f.water_boundary_total, rel=1e-12)


def test_dry_buffer_absorbs_and_warm_wall_releases():
    d = ConstantInputs(T_o=26.0, W_o=0.011)(0.0)
    x = equilibrium_state(26.0, 26.0, 0.011)
    x_dry = x.copy()
    x_dry[IDX_PHI_B] -= 0.1
    assert fluxes(0.0, x_dry, U_OFF, d, P).m_buffer < 0           # buffer takes water
    x_warm = x.copy()
    x_warm[IDX_TM] += 2.0
    assert fluxes(0.0, x_warm, U_OFF, d, P).m_buffer > 0          # warmer wall releases


def test_buffer_can_be_disabled():
    p = replace(P, moisture_buffer=False)
    d = ConstantInputs(T_o=26.0, W_o=0.011)(0.0)
    dx = rhs(0.0, np.array([26.0, 26.0, 0.011, 0.2, 420.0]), U_OFF, d, p)
    assert dx[IDX_W] == 0.0 and dx[IDX_PHI_B] == 0.0
