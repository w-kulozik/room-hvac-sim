"""CO2 balance and mechanical ventilation."""

from dataclasses import replace

import numpy as np
import pytest

from roomsim import verification as v
from roomsim.inputs import DAY, HOUR, ConstantInputs, SummerScenario
from roomsim.params import RoomParams
from roomsim.plant import IDX_CO2, IDX_V_VENT, U_OFF, equilibrium_state, fluxes, rhs
from roomsim.simulate import simulate

P = RoomParams()


def test_occupied_steady_state_co2():
    # 2 people x 5.2e-6 m3/s into 54 m3 at 0.5 ach: rise of 1387 ppm above outdoors
    x_ss = v.steady_state(P, ConstantInputs(n_occ=2))
    assert x_ss[IDX_CO2] - 420.0 == pytest.approx(2 * 5.2e-6 / (54 * 0.5 / 3600) * 1e6)
    assert x_ss[IDX_CO2] == pytest.approx(1807, abs=1)


def test_co2_decay_time_constant_with_ventilation():
    V_vent = 0.01                                   # 36 m3/h, two thirds of an air change
    tau = P.volume / (P.V_dot_inf + V_vent)

    class Fixed:
        def step(self, t, y):
            u = np.zeros(3); u[IDX_V_VENT] = V_vent
            return u

    x0 = equilibrium_state(26.0, 26.0, 0.01, C=1500.0)
    r = simulate(P, ConstantInputs(T_o=26.0, W_o=0.01), x0, 4 * HOUR, controller=Fixed())
    exact = 420.0 + 1080.0 * np.exp(-r.t / tau)
    # RK4 at dt/tau = 0.02: truncation error ~1e-9 of the 1080 ppm amplitude
    assert np.abs(r.CO2 - exact).max() < 1e-9 * 1080.0


def test_ventilation_equals_extra_infiltration():
    """V_vent must act on heat, water and CO2 exactly like the same flow of infiltration."""
    d = SummerScenario()(14 * HOUR)
    x = np.array([28.0, 27.0, 0.011, 0.5, 900.0])
    V_vent = 0.02
    u = np.zeros(3); u[IDX_V_VENT] = V_vent
    p_more_inf = replace(P, ach=P.ach + V_vent * 3600.0 / P.volume)
    assert rhs(0.0, x, u, d, P) == pytest.approx(rhs(0.0, x, U_OFF, d, p_more_inf), rel=1e-12)


def test_co2_closes_to_round_off():
    r = simulate(P, SummerScenario(), equilibrium_state(26.0, 26.0, 0.0116), 2 * DAY)
    assert np.abs(r.co2_residual).max() < 1e-13 * P.volume * 2000e-6 * 1e3
    assert r.CO2.max() > 1500.0          # the office is poorly ventilated without a fan


def test_ventilation_does_not_touch_the_mass_node():
    d = ConstantInputs(T_o=30.0)(0.0)
    x = equilibrium_state(24.0, 24.0, 0.01)
    u = np.zeros(3); u[IDX_V_VENT] = 0.05
    f = fluxes(0.0, x, u, d, P)
    assert f.q_ventilation == pytest.approx(P.rho_air * 0.05 * P.cp_air * 6.0)
    assert rhs(0.0, x, u, d, P)[1] == pytest.approx(rhs(0.0, x, U_OFF, d, P)[1])
