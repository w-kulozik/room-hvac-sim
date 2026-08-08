"""Regression tests for the control-area fixes (MPC, PI tuning, sampling period, daily windows)."""

from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from roomsim import psychro
from roomsim.ac_unit import ACUnit
from roomsim.control_model import RCModel
from roomsim.controllers import (
    CoilCascade,
    ComfortSchedule,
    DemandControlledPI,
    EconomicMPC,
    FanHumidityPI,
    OnOffThermostat,
    ScheduledVentilation,
    VentilationSchedule,
    WithVentilation,
    make_pi_controller,
)
from roomsim.controllers import mpc as mpc_module
from roomsim.estimation import DisturbanceKalman
from roomsim.inputs import HOUR, ConstantInputs, SummerScenario
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.simulate import simulate
from roomsim.tariff import TimeOfUseTariff

P = RoomParams()
UNIT = ACUnit()
MODEL = RCModel(P)
S = SummerScenario()
Y = np.array([25.0, 0.55, 600.0, 20.0])          # [T_a, RH, CO2, T_coil]


def _mpc(**kw: Any) -> EconomicMPC:
    return EconomicMPC(MODEL, UNIT, S, ventilation=VentilationSchedule(), **kw)


# ---- F1: band applied to the set point at the instant the planned temperature refers to -------

@pytest.mark.parametrize(("t_solve", "T_planned", "T_expected"), [
    (18 * HOUR, 27.0, 27.0),            # 18:15 is unoccupied: the band there goes up to 28 degC
    (7.75 * HOUR, 27.0, 25.4),          # 8:00 is occupied: 25 + 1 - 0.6 margin
])
def test_mpc_setpoint_clipped_with_band_at_end_of_first_step(
        monkeypatch: pytest.MonkeyPatch, t_solve: float, T_planned: float,
        T_expected: float) -> None:
    """The PI set point is the planned temperature at t + Ts_pred, T[0].

    It must therefore be kept in the planning band of that instant, the band the LP itself used
    for T[0]. Clipping with the band at the solve instant t held the PI at the occupied limit
    from 18:00 to 18:15, a quarter hour the planner lets float, and let it exceed the occupied
    limit from 7:45 to 8:00.
    """
    c = _mpc()

    def fake_plan(t: float, z: Any) -> Any:
        u = np.zeros(c.N)
        T = np.full(c.N, T_planned)
        c.last_solution = {"t": t + c.Ts_pred * np.arange(c.N), "u": u, "T": T,
                           "lo": T, "hi": T, "eps": u, "U_max": u, "cost": 0.0}
        return u

    monkeypatch.setattr(c, "plan", fake_plan)
    c.step(t_solve, Y)
    assert c.T_ref == pytest.approx(T_expected)


# ---- F3: LP capacity and cost on the same (net) cooling the Kalman filter is fed ----------------

def test_mpc_capacity_and_cost_net_of_indoor_fan_heat() -> None:
    """The LP variable u is the cooling of the RC model, the input of the Kalman filter.

    The filter is fed Q_sens - P_fan_in, the net heat removed from the air (the indoor-fan heat
    stays in the room). The capacity limit and the energy per unit of cooling must refer to that
    same quantity, otherwise U_max is too high by the fan heat and the cost too low.
    """
    c = _mpc()
    T_o = 32.0
    full = UNIT.operate(1.0, 1.0, c.comfort.T_set, c.W_ref, T_o)
    half = UNIT.operate(0.5, 1.0, c.comfort.T_set, c.W_ref, T_o)
    assert c.sensible_capacity(T_o) == pytest.approx(full.Q_sens - full.P_fan_in, rel=1e-12)
    assert c.energy_per_cooling(T_o) == pytest.approx(
        half.P_el / (half.Q_sens - half.P_fan_in), rel=1e-12)

    # the filter input at the same operating point is the LP capacity
    S_const = ConstantInputs(T_o=T_o)
    c2 = EconomicMPC(MODEL, UNIT, S_const)
    c2.inner.step = lambda t, y: np.array([1.0, 1.0, 0.0])        # type: ignore[method-assign]
    RH = float(psychro.RH_from_W(c2.comfort.T_set, c2.W_ref))
    c2.step(0.0, np.array([c2.comfort.T_set, RH, 420.0, 15.0]))
    assert c2.Q_c_prev == pytest.approx(c2.sensible_capacity(T_o), rel=1e-9)


# ---- F4: controller sampling period checked against the simulator's dt_ctrl --------------------

def _x0() -> Any:
    return equilibrium_state(26.0, 26.0, 0.0115)


@pytest.mark.parametrize("make", [
    lambda: make_pi_controller("sigma-delta", P, UNIT),
    lambda: make_pi_controller("hysteresis", P, UNIT),
    lambda: _mpc(),
    lambda: FanHumidityPI(make_pi_controller("sigma-delta", P, UNIT)),
    lambda: CoilCascade(make_pi_controller("sigma-delta", P, UNIT)),
    lambda: WithVentilation(OnOffThermostat(), DemandControlledPI(P)),
    lambda: WithVentilation(make_pi_controller("sigma-delta", P, UNIT), ScheduledVentilation()),
])
def test_simulate_rejects_controller_sampling_period_mismatch(make: Any) -> None:
    """Every discrete element uses its own Ts, not the elapsed time.

    Called every 30 s, a PI designed for 60 s integrates twice as fast and the Kalman filter
    predicts 60 s per 30 s sample. The simulator must refuse the combination.
    """
    with pytest.raises(ValueError, match="dt_ctrl"):
        simulate(P, S, _x0(), 600.0, dt=30.0, controller=make(), dt_ctrl=30.0, unit=UNIT)


def test_simulate_accepts_matching_or_unspecified_sampling_period() -> None:
    """A controller built for dt_ctrl, or one without a sampling period, runs."""
    simulate(P, S, _x0(), 600.0, dt=30.0, unit=UNIT, dt_ctrl=30.0,
             controller=make_pi_controller("sigma-delta", P, UNIT, Ts=30.0))
    simulate(P, S, _x0(), 600.0, dt=30.0, unit=UNIT, dt_ctrl=30.0, controller=OnOffThermostat())


def test_composite_controllers_reject_inconsistent_sampling_periods() -> None:
    """Parts of one controller that run at different periods are refused at construction."""
    with pytest.raises(ValueError, match="Ts"):
        FanHumidityPI(make_pi_controller("sigma-delta", P, UNIT), Ts=30.0)
    with pytest.raises(ValueError, match="Ts"):
        CoilCascade(make_pi_controller("sigma-delta", P, UNIT), Ts=30.0)
    with pytest.raises(ValueError, match="Ts"):
        WithVentilation(make_pi_controller("sigma-delta", P, UNIT), DemandControlledPI(P, Ts=30.0))
    with pytest.raises(ValueError, match="Ts"):
        _mpc(kalman=DisturbanceKalman(MODEL, Ts=30.0))
    with pytest.raises(ValueError, match="Ts"):
        _mpc(inner=make_pi_controller("sigma-delta", P, UNIT, Ts=30.0))


# ---- F5: daily windows that cross midnight -------------------------------------------------------

def test_tariff_peak_window_wraps_past_midnight() -> None:
    """A peak period from 20 h to 6 h covers 22 h and 2 h, not 12 h."""
    tou = TimeOfUseTariff(peak_start_h=20.0, peak_end_h=6.0)
    assert float(tou(22 * HOUR)) == tou.peak
    assert float(tou(2 * HOUR)) == tou.peak
    assert float(tou(12 * HOUR)) == tou.offpeak
    assert float(tou(6 * HOUR)) == tou.offpeak          # end excluded
    assert float(tou(20 * HOUR)) == tou.peak            # start included


def test_ventilation_schedule_wraps_past_midnight() -> None:
    """A night shift from 22 h to 6 h is ventilated at 23 h and 2 h, not at noon."""
    v = VentilationSchedule(start_h=22.0, end_h=6.0)
    assert v(23 * HOUR) == v.V_occupied
    assert v(2 * HOUR) == v.V_occupied
    assert v(12 * HOUR) == v.V_unoccupied


def test_comfort_schedule_wraps_past_midnight() -> None:
    """Occupancy, band edges and set point of a 22 h to 6 h schedule."""
    c = ComfortSchedule(start_h=22.0, end_h=6.0, early_start_h=21.0)
    assert bool(c.occupied(23 * HOUR)) and bool(c.occupied(2 * HOUR))
    assert not bool(c.occupied(12 * HOUR))
    assert bool(c.occupied_at_instant(6 * HOUR))         # the end instant belongs to occupancy
    lo, hi = c.bounds(np.array([2 * HOUR, 12 * HOUR]))
    assert list(hi) == [c.T_set + c.band, c.T_hi_unocc]
    assert c.setpoint(21.5 * HOUR) == c.T_set
    assert c.setpoint(3 * HOUR) == c.T_set
    assert c.setpoint(12 * HOUR) == c.T_hi_unocc


@pytest.mark.parametrize("build", [
    lambda: TimeOfUseTariff(peak_start_h=14.0, peak_end_h=14.0),
    lambda: TimeOfUseTariff(peak_start_h=-1.0, peak_end_h=6.0),
    lambda: VentilationSchedule(start_h=8.0, end_h=8.0),
    lambda: VentilationSchedule(start_h=8.0, end_h=25.0),
    lambda: ComfortSchedule(start_h=8.0, end_h=8.0),
    lambda: ComfortSchedule(early_start_h=18.0),
    lambda: ComfortSchedule(early_start_h=30.0),
])
def test_ill_defined_daily_windows_are_rejected(build: Any) -> None:
    """An empty window (start = end) or hours outside [0, 24] would silently never apply."""
    with pytest.raises(ValueError, match="_h"):
        build()


# ---- F6: plan after a failed LP solve ------------------------------------------------------------

def test_failed_solve_uses_old_plan_shifted_to_the_current_step(
        monkeypatch: pytest.MonkeyPatch) -> None:
    """When an LP solve fails, the previous plan is still the best information available.

    Its step 0 has just passed; the current step is its step 1. Both the set point (planned T at
    the end of the current step) and the planned move must come from index 1, not index 0.
    """
    c = _mpc()
    t0 = 10 * HOUR
    c.step(t0, Y)
    sol = c.last_solution
    assert sol is not None
    T_old, u_old = sol["T"].copy(), sol["u"].copy()

    def failing(*args: Any, **kwargs: Any) -> Any:
        return SimpleNamespace(status=2, x=None, fun=None)

    monkeypatch.setattr(mpc_module, "linprog", failing)
    for k in range(1, c.n_pred + 1):
        c.step(t0 + k * c.Ts, Y)
    assert c.solve_failures == 1
    t1 = t0 + c.Ts_pred
    lo, hi = c.planning_band(t1 + c.Ts_pred)
    assert c.T_ref == pytest.approx(float(np.clip(T_old[1], lo, hi)))
    assert c.u_plan == pytest.approx(float(u_old[1]))


# ---- F7: cached discretisation is read-only ------------------------------------------------------

def test_cached_discretisation_is_read_only() -> None:
    """``discrete`` returns cached arrays shared by every caller; writing to them must fail."""
    for arr in MODEL.discrete(60.0, 0.01):
        assert not arr.flags.writeable
        with pytest.raises(ValueError):
            arr[...] = 0.0
    Ad, Bd, Ed = MODEL.discrete(60.0, 0.01)
    assert np.all(np.isfinite(Ad)) and Ad[0, 0] != 0.0
