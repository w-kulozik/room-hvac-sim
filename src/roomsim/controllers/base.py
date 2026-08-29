"""Building blocks shared by the discrete controllers.

All controllers run every Ts seconds (the simulator's ``dt_ctrl``), read y = [T_a, RH_a, C]
and return the command [s, f, V_vent] for the air-conditioning unit. A controller whose
discrete elements assume a sampling period exposes it as the attribute ``Ts`` (s); the
simulator checks it against ``dt_ctrl``. Controllers that only use the time t have no ``Ts``.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Protocol

import numpy as np

from .._typing import FloatArray
from ..inputs import hour_of_day
from ..tariff import check_daily_window, in_daily_window

SetPoint = float | Callable[[float], float]
"""Constant set point (degC) or a schedule of time (s) returning degC."""


def set_point_at(T_set: SetPoint, t: float) -> float:
    """Return the set point at time t.

    Parameters
    ----------
    T_set : float or callable
        Set point, degC, or a function of time (s) returning it.
    t : float
        Time, s.

    Returns
    -------
    float
        Set point, degC.
    """
    return T_set(t) if callable(T_set) else T_set


class Controller(Protocol):
    """Interface of a discrete controller for the air-conditioning unit."""

    def step(self, t: float, y: FloatArray) -> FloatArray:
        """Return the command [s, f, V_vent] for measurement y at time t (s)."""
        ...


def sampling_period(obj: object) -> float | None:
    """Return the sampling period ``obj.Ts`` (s), or None if obj does not assume one.

    Parameters
    ----------
    obj : object
        Controller, controller part or ventilation strategy.

    Returns
    -------
    float or None
        Sampling period, s.
    """
    Ts = getattr(obj, "Ts", None)
    return None if Ts is None else float(Ts)


def common_sampling_period(owner: str, **periods: float | None) -> float | None:
    """Return the sampling period shared by the parts of a controller.

    Parameters
    ----------
    owner : str
        Name of the controller, used in the error message.
    **periods : float or None
        Sampling period of each part, s, by name; None for parts that do not assume one.

    Returns
    -------
    float or None
        The common sampling period, s; None if no part assumes one.

    Raises
    ------
    ValueError
        If two parts assume different sampling periods.
    """
    given = {name: Ts for name, Ts in periods.items() if Ts is not None}
    values = list(given.values())
    if any(not math.isclose(Ts, values[0], rel_tol=1e-9) for Ts in values[1:]):
        listing = ", ".join(f"{name} = {Ts!r} s" for name, Ts in given.items())
        raise ValueError(f"{owner}: the parts assume different sampling periods Ts ({listing}); "
                         "they must all be equal")
    return values[0] if values else None


@dataclass(frozen=True)
class VentilationSchedule:
    """Constant outdoor-air supply during occupied hours.

    Default: 0.0105 m3/s (38 m3/h). With 0.5 ach of infiltration this keeps the occupied
    steady-state CO2 at about 1000 ppm for two people. A period with start_h > end_h crosses
    midnight.

    Attributes
    ----------
    V_occupied : float
        Outdoor-air flow on [start_h, end_h), m3/s.
    V_unoccupied : float
        Outdoor-air flow outside those hours, m3/s.
    start_h, end_h : float
        Start and end of the occupied period, clock hours (h).

    Raises
    ------
    ValueError
        If a limit is outside [0, 24] or the two limits are equal (see
        ``roomsim.tariff.check_daily_window``).
    """

    V_occupied: float = 0.0105
    V_unoccupied: float = 0.0
    start_h: float = 8.0
    end_h: float = 18.0

    def __post_init__(self) -> None:
        """Validate the occupied period."""
        check_daily_window("start_h", self.start_h, "end_h", self.end_h)

    def __call__(self, t: float) -> float:
        """Return the scheduled outdoor-air flow at time t (s), m3/s."""
        h = float(hour_of_day(t))
        on = bool(in_daily_window(h, self.start_h, self.end_h))
        return self.V_occupied if on else self.V_unoccupied


NO_VENTILATION = VentilationSchedule(V_occupied=0.0)


class CompressorGuard:
    """Minimum run and stop times that protect the compressor from short cycling.

    Parameters
    ----------
    min_on : float, optional
        Minimum run time after a start, s.
    min_off : float, optional
        Minimum stop time after a stop, s.

    Attributes
    ----------
    on : bool
        Whether the compressor is running.
    t_switch : float
        Time of the last start or stop, s (-inf before the first switch).
    starts : int
        Number of starts so far.
    """

    def __init__(self, min_on: float = 180.0, min_off: float = 180.0) -> None:
        self.min_on, self.min_off = min_on, min_off
        self.on = False
        self.t_switch: float = -np.inf
        self.starts = 0

    def can_start(self, t: float) -> bool:
        """Return True if the compressor is off and has been off for at least ``min_off``.

        Parameters
        ----------
        t : float
            Current time, s.
        """
        return not self.on and t - self.t_switch >= self.min_off

    def can_stop(self, t: float) -> bool:
        """Return True if the compressor is on and has been on for at least ``min_on``.

        Parameters
        ----------
        t : float
            Current time, s.
        """
        return self.on and t - self.t_switch >= self.min_on

    def set(self, t: float, on: bool) -> None:
        """Record the compressor state at time t; a change of state resets the timer.

        Parameters
        ----------
        t : float
            Current time, s.
        on : bool
            New compressor state.
        """
        if on != self.on:
            self.on, self.t_switch = on, t
            self.starts += int(on)


def fan_command(s: float, f: float, mode: str) -> float:
    """Return the indoor-fan command for a compressor command and fan mode.

    Parameters
    ----------
    s : float
        Compressor command, -.
    f : float
        Requested fan speed, -.
    mode : str
        'auto' runs the fan only while the compressor runs; 'continuous' keeps it running.

    Returns
    -------
    float
        Fan command, -.

    Raises
    ------
    ValueError
        If ``mode`` is not 'auto' or 'continuous'.
    """
    if mode == "auto":
        return f if s > 0 else 0.0
    if mode == "continuous":
        return f
    raise ValueError(f"unknown fan mode {mode!r}")
