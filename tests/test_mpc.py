"""Control model, Kalman filter and economic MPC."""

import numpy as np
import pytest
from scipy.linalg import expm

from criteria import KH_PER_DAY, ROBUST_ENERGY_REL
from roomsim.ac_unit import ACUnit
from roomsim.control_model import RCModel
from roomsim.controllers import (
    ComfortSchedule,
    EconomicMPC,
    VentilationSchedule,
    make_pi_controller,
)
from roomsim.estimation import DisturbanceKalman, augmented
from roomsim.inputs import DAY, HOUR, ConstantInputs, SummerScenario
from roomsim.metrics import compute, energy_cost
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.simulate import simulate
from roomsim.tariff import TimeOfUseTariff

P = RoomParams()
UNIT = ACUnit()
MODEL = RCModel(P)
S = SummerScenario()
VENT = VentilationSchedule()


# ---- control model -----------------------------------------------------------------------

def test_discretisation_is_exact():
    A, B, E = MODEL.continuous(V_vent=0.01)
    Ad, Bd, Ed = MODEL.discrete(900.0, V_vent=0.01)
    assert Ad == pytest.approx(expm(A * 900.0), rel=1e-12)
    # held inputs: the continuous steady state is a fixed point of the discrete map
    u, e = 500.0, np.array([30.0, 200.0, 300.0])
    x_ss = -np.linalg.solve(A, B * u + E @ e)
    assert Ad @ x_ss + Bd * u + Ed @ e == pytest.approx(x_ss, rel=1e-12)


def test_model_reproduces_the_plant_with_the_unit_off():
    """Same parameters, unit off: the RC model is the plant's thermal part, up to input sampling."""
    x0 = equilibrium_state(27.0, 28.0, 0.0116)
    r = simulate(P, S, x0, DAY, dt=60.0)
    z = x0[:2].copy()
    errs = []
    for k in range(len(r.t) - 1):
        Ad, Bd, Ed = MODEL.discrete(60.0)
        e = MODEL.gains(S(r.t[k] + 30.0))[0]          # inputs at the middle of the step
        z = Ad @ z + Ed @ e
        errs.append(abs(z[0] - r.T_a[k + 1]))
    # the only difference is sampling T_o and the sun at the middle of each 60 s step instead of
    # integrating them: an O(dt^2) error, measured at 2e-5 K
    assert max(errs) < 1e-4


def test_mismatched_model_scales_primary_parameters():
    m = RCModel.mismatched(P, {"C_m": 1.2, "furniture_factor": 0.8})
    assert m.p.C_m == pytest.approx(1.2 * P.C_m)
    assert m.p.C_a == pytest.approx(0.8 * P.C_a)


# ---- Kalman filter -----------------------------------------------------------------------

def test_augmented_system_is_observable_from_the_air_sensor():
    A, _, _ = augmented(MODEL, 60.0, 0.0)
    C = np.array([[1.0, 0.0, 0.0]])
    obs = np.vstack([C @ np.linalg.matrix_power(A, k) for k in range(3)])
    assert np.linalg.matrix_rank(obs) == 3


def test_kalman_recovers_mass_temperature_and_an_unknown_heat_input():
    """Truth: the controller's own model plus 250 W of heat into the air it does not know about."""
    Ad, Bd, Ed = MODEL.discrete(60.0)
    d_true, x = 250.0, np.array([25.0, 28.0])
    kf = DisturbanceKalman(MODEL, 60.0)
    d0 = ConstantInputs(T_o=30.0)(0.0)
    e = MODEL.gains(d0)[0]
    kf.initialise(x[0])                                  # starts with T_m = T_a: 3 K wrong
    err = []
    for k in range(int(6 * DAY / 60)):
        x = Ad @ x + Ed @ e - Bd * d_true
        kf.update(x[0], 0.0, d0, 0.0)
        if (k + 1) % int(DAY / 60) == 0:
            err.append(abs(kf.z[1] - x[1]))
    # unbiased, but slow: with only the air sensor, a wrong T_m and an unknown heat input are told
    # apart through the slow mass mode, so the error falls by a factor of ~6 per day (tau ~ 14 h)
    assert kf.z[1] == pytest.approx(x[1], abs=1e-3)
    assert kf.z[2] == pytest.approx(d_true, rel=1e-3)
    assert 3.0 < err[2] / err[3] < 10.0


# ---- comfort schedule ---------------------------------------------------------------------

def test_comfort_band_edges_belong_to_the_occupied_period():
    c = ComfortSchedule()
    assert c.occupied_at_instant(8 * HOUR) and c.occupied_at_instant(18 * HOUR)
    assert not c.occupied_at_instant(8 * HOUR - 60) and not c.occupied_at_instant(18 * HOUR + 60)
    assert c.setpoint(7 * HOUR) == c.T_hi_unocc
    assert ComfortSchedule(early_start_h=6.0).setpoint(7 * HOUR) == c.T_set


# ---- planner ------------------------------------------------------------------------------

Z0 = np.array([26.5, 27.0, 0.0])
T_PLAN = 2 * DAY + 6 * HOUR


def plan(**kw):
    c = EconomicMPC(MODEL, UNIT, S, ventilation=VENT, **kw)
    u = c.plan(T_PLAN, Z0)
    return c, u, c.last_solution


def test_plan_is_feasible_and_rides_the_upper_band_edge():
    c, u, sol = plan()
    assert np.all(u >= -1e-9) and np.all(u <= sol["U_max"] + 1e-6)
    assert sol["eps"].max() < 1e-6                       # no comfort violation planned
    lo, hi = c.planning_band(sol["t"] + c.Ts_pred)
    occ = c.comfort.occupied_at_instant(sol["t"] + c.Ts_pred)
    assert np.all(sol["T"][occ] <= hi[occ] + 1e-6) and np.all(sol["T"][occ] >= lo[occ] - 1e-6)
    # cooling costs energy, so the optimum sits at the warm edge of the band most of the day
    assert np.median(hi[occ] - sol["T"][occ]) < 0.05


def test_no_cooling_planned_while_the_band_is_wide():
    c = EconomicMPC(MODEL, UNIT, S, ventilation=VENT)
    u = c.plan(2 * DAY + 20 * HOUR, np.array([26.0, 26.5, 0.0]))
    assert np.all(u[: int(8 * HOUR / c.Ts_pred)] < 1.0)  # 20:00 - 04:00: room may float


def test_expensive_peak_shifts_cooling_out_of_it():
    def peak_cooling(ratio):
        _, u, sol = plan(tariff=TimeOfUseTariff(offpeak=0.15, peak=0.15 * ratio))
        h = (sol["t"] - sol["t"][0]) / HOUR + 6.0
        return u[(h >= 14) & (h < 20)].sum()
    assert peak_cooling(1.0) == pytest.approx(peak_cooling(2.0), rel=0.01)   # not worth it yet
    assert peak_cooling(6.0) < 0.9 * peak_cooling(1.0)


# ---- closed loop --------------------------------------------------------------------------

X0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
MISMATCH = {"C_m": 1.2, "furniture_factor": 0.8, "H_mo": 0.8, "h_i": 1.2, "g_win": 0.8}


@pytest.fixture(scope="module")
def loops():
    ctrls = {
        "pi": make_pi_controller("sigma-delta", P, UNIT, ventilation=VENT),
        "mpc": EconomicMPC(MODEL, UNIT, S, ventilation=VENT),
        "mismatch": EconomicMPC(RCModel.mismatched(P, MISMATCH), UNIT, S, ventilation=VENT),
    }
    out = {}
    for k, c in ctrls.items():
        r = simulate(P, S, X0, 4 * DAY, controller=c, unit=UNIT)
        out[k] = (c, r.window(2 * DAY))
    return out


def test_mpc_keeps_comfort_and_saves_energy(loops):
    m_pi, m_mpc = (compute(loops[k][1], 25.0) for k in ("pi", "mpc"))
    assert m_mpc.Kh_outside_occ < KH_PER_DAY * 2
    assert m_mpc.energy_kWh < m_pi.energy_kWh            # size of the saving: see docs/results.md


def test_mpc_is_robust_to_model_mismatch(loops):
    m = compute(loops["mismatch"][1], 25.0)
    m_nom = compute(loops["mpc"][1], 25.0)
    assert m.Kh_outside_occ < KH_PER_DAY * 2
    assert m.energy_kWh == pytest.approx(m_nom.energy_kWh, rel=ROBUST_ENERGY_REL)


def test_kalman_tracks_the_mass_in_closed_loop(loops):
    c, r = loops["mpc"]
    # the controller's thermal model is exact here, so the estimate must track the mass to
    # sensor-level accuracy; what the model misses (unit, humidity) goes into the disturbance
    assert c.kf.z[1] == pytest.approx(r.T_m[-1], abs=0.01)
    assert c.solve_failures == 0


def test_energy_cost_metric():
    tou = TimeOfUseTariff(offpeak=0.1, peak=0.3)
    r = simulate(P, ConstantInputs(T_o=30.0), X0, 2 * HOUR, t_start=13 * HOUR,
                 controller=make_pi_controller("sigma-delta", P, UNIT), unit=UNIT)
    P_el = r.unit.op["P_el"][:-1]
    h = r.t[:-1] / HOUR
    expected = (np.sum(P_el[h < 14]) * 0.1 + np.sum(P_el[h >= 14]) * 0.3) * r.dt / 3.6e6
    assert energy_cost(r, tou) == pytest.approx(expected, rel=1e-12)
