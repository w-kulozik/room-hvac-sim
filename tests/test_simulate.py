"""Integrators, the simulation loop (control hold, input sampling) and periodic shooting."""

import numpy as np
import pytest

from roomsim.inputs import DAY, HOUR, ConstantInputs, SummerScenario
from roomsim.integrators import euler_stability_limit, euler_step, rk4_step
from roomsim.params import RoomParams
from roomsim.plant import IDX_Q_SENS, U_OFF, equilibrium_state, jacobian
from roomsim.simulate import periodic_state, simulate

P = RoomParams()


# ---- integrators on a scalar test equation ----------------------------------------------

@pytest.mark.parametrize("step, exact_one_step", [
    (euler_step, lambda z: 1 + z),
    (rk4_step, lambda z: 1 + z + z**2 / 2 + z**3 / 6 + z**4 / 24),
])
def test_amplification_factor(step, exact_one_step):
    lam, dt = -0.3, 0.7
    x1 = step(lambda t, x: lam * x, 0.0, np.array([1.0]), dt)
    assert x1[0] == pytest.approx(exact_one_step(lam * dt), rel=1e-15)


def test_rk4_is_exact_for_cubic_in_time():
    # dx/dt = 3 t^2 + 2 t  ->  x(t) = t^3 + t^2 ; Simpson weights make RK4 exact here
    x1 = rk4_step(lambda t, x: np.array([3 * t**2 + 2 * t]), 0.0, np.array([0.0]), 2.0)
    assert x1[0] == pytest.approx(12.0, rel=1e-15)


def test_euler_stability_limit_scalar_and_room():
    assert euler_stability_limit(np.array([[-0.01]])) == pytest.approx(200.0)
    S = SummerScenario()
    x0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
    dt_max = euler_stability_limit(jacobian(0.0, x0, U_OFF, S(0.0), P))
    # fastest mode is the air-buffer moisture exchange, not the air temperature
    assert 600.0 < dt_max < 1000.0
    assert dt_max > 60.0


# ---- simulation loop -------------------------------------------------------------------

class RecordingController:
    """Returns a different heat input at every call and records when it was called."""

    def __init__(self):
        self.calls = []

    def step(self, t, y):
        self.calls.append((t, y.copy()))
        return np.array([-100.0 * len(self.calls), 0.0, 0.0])


def test_controller_sampled_and_held():
    ctrl = RecordingController()
    x0 = equilibrium_state(26.0, 26.0, 0.011)
    r = simulate(P, ConstantInputs(), x0, 600.0, dt=20.0, controller=ctrl, dt_ctrl=120.0)
    assert [c[0] for c in ctrl.calls] == [0.0, 120.0, 240.0, 360.0, 480.0]
    # u is constant inside each 120 s sampling interval (6 steps of 20 s)
    q = r.u[:-1, IDX_Q_SENS].reshape(5, 6)
    assert np.all(q == q[:, :1])
    assert np.allclose(q[:, 0], [-100, -200, -300, -400, -500])
    # the controller sees only T_a, RH and CO2
    assert ctrl.calls[0][1].shape == (3,)


def test_controller_period_must_be_multiple_of_step():
    with pytest.raises(ValueError):
        simulate(P, ConstantInputs(), equilibrium_state(26, 26, 0.011), 600.0, dt=60.0,
                 controller=RecordingController(), dt_ctrl=90.0)


def test_occupancy_step_is_not_anticipated():
    """With the step ending exactly at 8:00, no vapour from occupants enters before 8:00."""
    S = SummerScenario()
    x0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
    r = simulate(P, S, x0, 8 * HOUR, dt=600.0, t_start=7 * HOUR)
    r_empty = simulate(P, SummerScenario(n_people=0.0), x0, 8 * HOUR, dt=600.0, t_start=7 * HOUR)
    assert np.array_equal(r.x, r_empty.x)


def test_logged_fluxes_match_state():
    S = SummerScenario()
    r = simulate(P, S, equilibrium_state(26.0, 26.0, 0.011), 13 * HOUR, t_start=11 * HOUR)
    k = 37
    assert r.q["q_mass_to_air"][k] == pytest.approx(P.H_am * (r.T_m[k] - r.T_a[k]))
    assert r.d.G_sol[k] == pytest.approx(float(S.G_sol(r.t[k])))


# ---- periodic steady state -------------------------------------------------------------

@pytest.fixture(scope="module")
def periodic():
    S = SummerScenario()
    x_p = periodic_state(P, S, equilibrium_state(26.0, 26.0, float(S.W_o(0.0))))
    return S, x_p, simulate(P, S, x_p, DAY)


def test_periodic_state_repeats(periodic):
    _, x_p, r = periodic
    assert r.x[-1] == pytest.approx(x_p, abs=1e-6)


def test_daily_mean_equals_steady_state_of_mean_inputs(periodic):
    """For a linear time-invariant system the mean of the periodic response equals the steady
    state under the mean inputs. This is an independent analytical check of the thermal part."""
    S, _, r = periodic
    t = r.t[:-1]
    T_o = S.T_o(t).mean()
    q_sol = P.g_win * P.A_win * S.G_sol(t).mean()
    q_occ = P.q_person_sens * S.n_occ(t).mean()
    A = np.array([[-(P.H_am + P.H_ao + P.H_inf), P.H_am],
                  [P.H_am, -(P.H_am + P.H_mo)]])
    b = -np.array([(P.H_ao + P.H_inf) * T_o + P.f_sol_air * q_sol + P.f_person_conv * q_occ,
                   P.H_mo * T_o + (1 - P.f_sol_air) * q_sol + (1 - P.f_person_conv) * q_occ])
    T_a_ss, T_m_ss = np.linalg.solve(A, b)
    # left-rectangle mean of a smooth periodic function is spectrally accurate
    assert r.T_a[:-1].mean() == pytest.approx(T_a_ss, abs=2e-3)
    assert r.T_m[:-1].mean() == pytest.approx(T_m_ss, abs=2e-3)
