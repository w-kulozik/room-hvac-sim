"""Building blocks shared by the discrete controllers.

All controllers run every Ts seconds (the simulator's ``dt_ctrl``), read y = [T_a, RH_a, C]
and return the command [s, f, V_vent] for the air-conditioning unit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

from .._typing import FloatArray
from ..inputs import hour_of_day


class Controller(Protocol):
    """Interface of a discrete controller for the air-conditioning unit."""

    def step(self, t: float, y: FloatArray) -> FloatArray:
        """Return the command [s, f, V_vent] for measurement y at time t (s)."""
        ...


@dataclass(frozen=True)
class VentilationSchedule:
    """Constant outdoor-air supply during occupied hours.

    Default: 0.0105 m3/s (38 m3/h). With 0.5 ach of infiltration this keeps the occupied
    steady-state CO2 at about 1000 ppm for two people.

    Attributes
    ----------
    V_occupied : float
        Outdoor-air flow on [start_h, end_h), m3/s.
    V_unoccupied : float
        Outdoor-air flow outside those hours, m3/s.
    start_h, end_h : float
        Start and end of the occupied period, clock hours (h).
    """

    V_occupied: float = 0.0105
    V_unoccupied: float = 0.0
    start_h: float = 8.0
    end_h: float = 18.0

    def __call__(self, t: float) -> float:
        """Return the scheduled outdoor-air flow at time t (s), m3/s."""
        h = float(hour_of_day(t))
        return self.V_occupied if self.start_h <= h < self.end_h else self.V_unoccupied


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
