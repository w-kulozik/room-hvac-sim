"""Ventilation strategies, combined with any temperature/humidity controller.

A strategy has ``step(t, y) -> V_vent`` (m3/s). ``WithVentilation`` runs the inner controller
and replaces the ventilation element of its command, so the same compressor and fan control is
used with every strategy.

    ScheduledVentilation  fixed flow during occupied hours
    DemandControlledPI    PI on measured CO2 (demand-controlled ventilation, continuous)
    StagedVentilation     discrete fan levels switched on CO2 with hysteresis and dwell time
    NightFlush            adds free cooling with outdoor air when the office is empty
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Protocol

import numpy as np

from .._typing import FloatArray, FloatOrArray
from ..params import RoomParams
from ..plant import IDX_Y_CO2, IDX_Y_TA
from .base import Controller, VentilationSchedule
from .mpc import ComfortSchedule
from .pi import DiscretePI, simc_pi

IDX_VENT = 2
CO2_OUTDOOR = 420.0


class VentilationStrategy(Protocol):
    """Interface of a ventilation strategy."""

    def step(self, t: float, y: FloatArray) -> float:
        """Return the outdoor-air flow for measurement y at time t (s), m3/s."""
        ...


class WithVentilation:
    """Controller wrapper whose ventilation command comes from a separate strategy.

    Parameters
    ----------
    inner : Controller
        Controller providing the compressor and fan commands.
    strategy : VentilationStrategy
        Strategy providing the outdoor-air flow.
    """

    def __init__(self, inner: Controller, strategy: VentilationStrategy) -> None:
        self.inner, self.strategy = inner, strategy

    def step(self, t: float, y: FloatArray) -> FloatArray:
        """Return the inner command with its ventilation element replaced.

        Parameters
        ----------
        t : float
            Current time, s.
        y : ndarray
            Measurements, passed to the inner controller and the strategy.

        Returns
        -------
        ndarray, shape (3,)
            Command [s (-), f (-), V_vent (m3/s)].
        """
        cmd = np.array(self.inner.step(t, y), dtype=float)
        cmd[IDX_VENT] = self.strategy.step(t, y)
        return cmd


class ScheduledVentilation:
    """Fixed outdoor-air flow by time of day.

    Parameters
    ----------
    schedule : VentilationSchedule, optional
        Flow schedule.
    """

    def __init__(
        self, schedule: VentilationSchedule = VentilationSchedule(),  # noqa: B008  frozen dataclass
    ) -> None:
        self.schedule = schedule

    def step(self, t: float, y: FloatArray) -> float:
        """Return the scheduled flow at time t (s), m3/s; y is not used."""
        return self.schedule(t)


def dcv_tuning(p: RoomParams, C_set: float = 900.0, n_people: float = 2.0,
               tau_c: float = 900.0, Ts: float = 60.0) -> tuple[float, float]:
    """Return SIMC gains for the CO2 loop, linearised at the set point with n_people present.

    At the set point the total outdoor-air flow is Q = N G / (C_set - C_o); the process is first
    order with gain dC/dV = -(C_set - C_o)/Q (ppm per m3/s) and time constant V_room / Q.

    Parameters
    ----------
    p : RoomParams
        Room parameters (``co2_person`` in m3/s per person, ``volume`` in m3).
    C_set : float, optional
        CO2 set point, ppm.
    n_people : float, optional
        Number of occupants at the linearisation point, -.
    tau_c : float, optional
        Desired closed-loop time constant, s.
    Ts : float, optional
        Sampling period, s.

    Returns
    -------
    Kp : float
        Proportional gain, m3/s per ppm (positive: more CO2, more air).
    Ti : float
        Integral time, s.
    """
    Q = n_people * p.co2_person * 1e6 / (C_set - CO2_OUTDOOR)
    K = (C_set - CO2_OUTDOOR) / Q
    tau = p.volume / Q
    return simc_pi(K, tau, Ts / 2.0, tau_c)


class DemandControlledPI:
    """Continuous demand-controlled ventilation: PI on CO2 within [V_min, V_max] when occupied.

    Parameters
    ----------
    p : RoomParams
        Room parameters used for tuning.
    C_set : float, optional
        CO2 set point, ppm.
    V_min, V_max : float, optional
        Flow limits while occupied, m3/s.
    comfort : ComfortSchedule, optional
        Occupancy schedule; the flow is zero while unoccupied.
    Ts : float, optional
        Sampling period, s.
    """

    def __init__(self, p: RoomParams, C_set: float = 900.0, V_min: float = 0.006,
                 V_max: float = 0.03,
                 comfort: ComfortSchedule = ComfortSchedule(),  # noqa: B008  (frozen dataclass)
                 Ts: float = 60.0) -> None:
        Kp, Ti = dcv_tuning(p, C_set, Ts=Ts)
        self.pi = DiscretePI(Kp, Ti, Ts)
        self.C_set, self.V_min, self.V_max, self.comfort = C_set, V_min, V_max, comfort

    def step(self, t: float, y: FloatArray) -> float:
        """Return the outdoor-air flow for one sample and update the PI.

        Parameters
        ----------
        t : float
            Current time, s.
        y : ndarray
            Measurements; uses the CO2 concentration y[IDX_Y_CO2], ppm.

        Returns
        -------
        float
            Outdoor-air flow, m3/s.
        """
        if not self.comfort.occupied(t):
            return 0.0                                     # integrator kept for the next day
        e = y[IDX_Y_CO2] - self.C_set
        v = self.pi.demand(e)
        V = min(max(v, self.V_min), self.V_max)
        self.pi.update(e, v, V)
        return V


class StagedVentilation:
    """Discrete fan levels switched on CO2, with hysteresis and a minimum dwell time per level.

    Parameters
    ----------
    levels : sequence of float, optional
        Flow of each stage, m3/s, in increasing order.
    up : sequence of float, optional
        CO2 at or above which stage i moves to i + 1, ppm (one per pair of levels).
    down : sequence of float, optional
        CO2 at or below which stage i + 1 moves to i, ppm (one per pair of levels).
    dwell : float, optional
        Minimum time between stage changes, s.
    comfort : ComfortSchedule, optional
        Occupancy schedule; the flow is zero and the stage reset while unoccupied.

    Raises
    ------
    ValueError
        If ``up`` and ``down`` do not both have ``len(levels) - 1`` elements.
    """

    def __init__(self, levels: Sequence[float] = (0.006, 0.012, 0.024),
                 up: Sequence[float] = (850.0, 950.0), down: Sequence[float] = (750.0, 850.0),
                 dwell: float = 300.0,
                 comfort: ComfortSchedule = ComfortSchedule(),  # noqa: B008  (frozen dataclass)
                 ) -> None:
        if not (len(up) == len(down) == len(levels) - 1):
            raise ValueError("need one up and one down threshold between consecutive levels")
        self.levels, self.up, self.down = levels, up, down
        self.dwell, self.comfort = dwell, comfort
        self.stage, self.t_switch = 0, -np.inf

    def step(self, t: float, y: FloatArray) -> float:
        """Return the outdoor-air flow of the current stage after a possible switch.

        Parameters
        ----------
        t : float
            Current time, s.
        y : ndarray
            Measurements; uses the CO2 concentration y[IDX_Y_CO2], ppm.

        Returns
        -------
        float
            Outdoor-air flow, m3/s.
        """
        if not self.comfort.occupied(t):
            self.stage, self.t_switch = 0, -np.inf
            return 0.0
        C, i = y[IDX_Y_CO2], self.stage
        if t - self.t_switch >= self.dwell:
            if i < len(self.levels) - 1 and self.up[i] <= C:
                self.stage, self.t_switch = i + 1, t
            elif i > 0 and self.down[i - 1] >= C:
                self.stage, self.t_switch = i - 1, t
        return self.levels[self.stage]


class NightFlush:
    """Add free cooling with outdoor air while the office is empty.

    Full flow when the outdoor air is at least dT cooler than the room and the room is above
    T_min; otherwise the inner strategy's flow. ``outdoor`` is the outdoor sensor.

    Parameters
    ----------
    inner : VentilationStrategy
        Strategy used while occupied and when flushing does not apply.
    outdoor : callable
        Outdoor temperature as a function of time (s), degC.
    V_flush : float, optional
        Flushing flow, m3/s.
    dT : float, optional
        Minimum indoor-outdoor temperature difference for flushing, K.
    T_min : float, optional
        Air temperature at or below which flushing stops, degC.
    comfort : ComfortSchedule, optional
        Occupancy schedule.
    """

    def __init__(self, inner: VentilationStrategy, outdoor: Callable[[float], FloatOrArray],
                 V_flush: float = 0.03, dT: float = 2.0, T_min: float = 22.0,
                 comfort: ComfortSchedule = ComfortSchedule(),  # noqa: B008  (frozen dataclass)
                 ) -> None:
        self.inner, self.outdoor = inner, outdoor
        self.V_flush, self.dT, self.T_min, self.comfort = V_flush, dT, T_min, comfort

    def step(self, t: float, y: FloatArray) -> float:
        """Return the outdoor-air flow for one sample.

        Parameters
        ----------
        t : float
            Current time, s.
        y : ndarray
            Measurements; uses the air temperature y[IDX_Y_TA], degC (and whatever the inner
            strategy uses).

        Returns
        -------
        float
            Outdoor-air flow, m3/s.
        """
        V = self.inner.step(t, y)
        if self.comfort.occupied(t):
            return V
        T_a, T_o = y[IDX_Y_TA], float(self.outdoor(t))
        if T_o <= T_a - self.dT and T_a > self.T_min:
            return self.V_flush
        return V
