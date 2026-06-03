"""Humidity-aware control of an inverter unit through the indoor fan.

The compressor stays under the temperature PI (with its low-load handler). The indoor fan is
the dehumidification actuator: less airflow means a colder coil surface, more condensation and a
lower sensible heat ratio, at a small cost in efficiency. At minimum compressor speed, going from
full to 40 % fan takes the coil from about 19 degC (dry) to 15 degC (wet), SHR 1.00 -> 0.75
(room at 25 degC and 60 % RH, 32 degC outdoors).

H1  FanHumidityPI: one PI on relative humidity acts directly on the fan speed.
        d = PI(RH - RH_set) in [0, 1],   f = 1 - (1 - f_min) d

H2  CoilCascade: cascade with the evaporator thermistor as the intermediate measurement.
        outer  PI(RH - RH_set) -> coil temperature set point T_s* in [T_s_min, T_s_max]
        inner  integral action on the fan:  f <- f - Ki (T_coil - T_s*), clipped to [f_min, 1]
    The inner loop corrects anything that moves the coil temperature (compressor speed changes,
    return-air state) within a few samples, before it shows up in the room humidity.

Both freeze their integrators while the compressor is off, since the fan cannot dehumidify then,
and use back-calculation against their own output limits.
"""

from __future__ import annotations

import numpy as np

from ..plant import IDX_Y_RH, IDX_Y_TCOIL
from .pi import DiscretePI
from .pi_unit import PIUnitController

IDX_S, IDX_F = 0, 1


class FanHumidityPI:
    """H1: relative-humidity PI on the indoor fan speed."""

    def __init__(self, temperature: PIUnitController, RH_set: float = 0.55,
                 Kp: float = 0.05, Ti: float = 1800.0, f_min: float = 0.4, Ts: float = 60.0):
        self.temperature, self.RH_set, self.f_min = temperature, RH_set, f_min
        self.pi = DiscretePI(Kp, Ti, Ts)      # demand per %RH

    def step(self, t: float, y: np.ndarray) -> np.ndarray:
        cmd = self.temperature.step(t, y)
        e = 100.0 * (y[IDX_Y_RH] - self.RH_set)
        v = self.pi.demand(e)
        d = min(max(v, 0.0), 1.0)
        running = cmd[IDX_S] > 0
        self.pi.update(e, v, d, integrate=running)
        if running:
            cmd[IDX_F] = 1.0 - (1.0 - self.f_min) * d
        return cmd


class CoilCascade:
    """H2: RH (outer) -> coil temperature set point -> fan speed (inner)."""

    def __init__(self, temperature: PIUnitController, RH_set: float = 0.55,
                 Kp: float = 0.3, Ti: float = 1800.0, T_s_max: float = 20.0, T_s_min: float = 8.0,
                 Ki_fan: float = 0.04, f_min: float = 0.4, settle_samples: int = 2,
                 Ts: float = 60.0):
        self.temperature, self.RH_set = temperature, RH_set
        self.T_s_max, self.T_s_min = T_s_max, T_s_min
        self.Ki_fan, self.f_min, self.settle = Ki_fan, f_min, settle_samples
        self.outer = DiscretePI(Kp, Ti, Ts)   # coil set-point depression, K per %RH
        self.f = 1.0
        self.T_s_set = T_s_max
        self.samples_on = 0

    def step(self, t: float, y: np.ndarray) -> np.ndarray:
        cmd = self.temperature.step(t, y)
        running = cmd[IDX_S] > 0
        self.samples_on = self.samples_on + 1 if running else 0

        # outer loop: how far below T_s_max the coil should run
        e = 100.0 * (y[IDX_Y_RH] - self.RH_set)
        v = self.outer.demand(e)
        depression = min(max(v, 0.0), self.T_s_max - self.T_s_min)
        self.outer.update(e, v, depression, integrate=running)
        self.T_s_set = self.T_s_max - depression

        # inner loop: integral action on the fan, valid once the coil has cooled after a start
        if self.samples_on > self.settle:
            self.f -= self.Ki_fan * (y[IDX_Y_TCOIL] - self.T_s_set)   # coil too warm: slow the fan
            self.f = min(max(self.f, self.f_min), 1.0)
        if running:
            cmd[IDX_F] = self.f
        return cmd
