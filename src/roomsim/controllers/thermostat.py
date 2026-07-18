"""On/off thermostat with hysteresis: the reference every other strategy is compared with."""

from __future__ import annotations

import numpy as np

from ..plant import IDX_Y_TA
from .base import NO_VENTILATION, CompressorGuard, VentilationSchedule, fan_command


class OnOffThermostat:
    """Cooling thermostat: start at T_set + band/2, stop at T_set - band/2, fixed speed s_on.

    The compressor guard enforces minimum run and stop times.
    """

    def __init__(self, T_set: float = 25.0, band: float = 1.0, s_on: float = 1.0,
                 fan: float = 1.0, fan_mode: str = "auto",
                 ventilation: VentilationSchedule = NO_VENTILATION,
                 min_on: float = 180.0, min_off: float = 180.0):
        self.T_set, self.band, self.s_on = T_set, band, s_on
        self.fan, self.fan_mode, self.ventilation = fan, fan_mode, ventilation
        self.guard = CompressorGuard(min_on, min_off)

    def step(self, t: float, y: np.ndarray) -> np.ndarray:
        T, g = y[IDX_Y_TA], self.guard
        if not g.on and self.T_set + 0.5 * self.band <= T and g.can_start(t):
            g.set(t, True)
        elif g.on and self.T_set - 0.5 * self.band >= T and g.can_stop(t):
            g.set(t, False)
        s = self.s_on if g.on else 0.0
        return np.array([s, fan_command(s, self.fan, self.fan_mode), self.ventilation(t)])
