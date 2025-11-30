"""Building blocks shared by the discrete controllers.

All controllers run every Ts seconds (the simulator's ``dt_ctrl``), read y = [T_a, RH_a, C]
and return the command [s, f, V_vent] for the air-conditioning unit.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..inputs import hour_of_day


@dataclass(frozen=True)
class VentilationSchedule:
    """Constant outdoor-air supply during occupied hours.

    Default: 0.0105 m3/s (38 m3/h). With 0.5 ach of infiltration this keeps the occupied
    steady-state CO2 at about 1000 ppm for two people.
    """

    V_occupied: float = 0.0105
    V_unoccupied: float = 0.0
    start_h: float = 8.0
    end_h: float = 18.0

    def __call__(self, t: float) -> float:
        h = float(hour_of_day(t))
        return self.V_occupied if self.start_h <= h < self.end_h else self.V_unoccupied


NO_VENTILATION = VentilationSchedule(V_occupied=0.0)


class CompressorGuard:
    """Minimum run and stop times that protect the compressor from short cycling."""

    def __init__(self, min_on: float = 180.0, min_off: float = 180.0):
        self.min_on, self.min_off = min_on, min_off
        self.on = False
        self.t_switch = -np.inf
        self.starts = 0

    def can_start(self, t: float) -> bool:
        return not self.on and t - self.t_switch >= self.min_off

    def can_stop(self, t: float) -> bool:
        return self.on and t - self.t_switch >= self.min_on

    def set(self, t: float, on: bool) -> None:
        if on != self.on:
            self.on, self.t_switch = on, t
            self.starts += int(on)


def fan_command(s: float, f: float, mode: str) -> float:
    """Indoor fan: 'auto' runs it with the compressor, 'continuous' keeps it running."""
    if mode == "auto":
        return f if s > 0 else 0.0
    if mode == "continuous":
        return f
    raise ValueError(f"unknown fan mode {mode!r}")
