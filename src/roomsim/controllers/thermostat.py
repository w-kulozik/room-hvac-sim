"""On/off thermostat with hysteresis: the reference every other strategy is compared with."""

from __future__ import annotations

import numpy as np

from .._typing import FloatArray
from ..plant import IDX_Y_TA
from .base import NO_VENTILATION, CompressorGuard, VentilationSchedule, fan_command


class OnOffThermostat:
    """Cooling thermostat: start at T_set + band/2, stop at T_set - band/2, fixed speed s_on.

    The compressor guard enforces minimum run and stop times.

    Parameters
    ----------
    T_set : float, optional
        Set point, degC.
    band : float, optional
        Hysteresis band (total width), K.
    s_on : float, optional
        Compressor command while running, -.
    fan : float, optional
        Indoor-fan speed, -.
    fan_mode : str, optional
        'auto' or 'continuous' (see ``fan_command``).
    ventilation : VentilationSchedule, optional
        Outdoor-air flow schedule.
    min_on, min_off : float, optional
        Minimum run and stop times of the compressor, s.
    """

    def __init__(self, T_set: float = 25.0, band: float = 1.0, s_on: float = 1.0,
                 fan: float = 1.0, fan_mode: str = "auto",
                 ventilation: VentilationSchedule = NO_VENTILATION,
                 min_on: float = 180.0, min_off: float = 180.0) -> None:
        self.T_set, self.band, self.s_on = T_set, band, s_on
        self.fan, self.fan_mode, self.ventilation = fan, fan_mode, ventilation
        self.guard = CompressorGuard(min_on, min_off)

    def step(self, t: float, y: FloatArray) -> FloatArray:
        """Return the command for one controller sample.

        Parameters
        ----------
        t : float
            Current time, s.
        y : ndarray
            Measurements; uses the air temperature y[IDX_Y_TA], degC.

        Returns
        -------
        ndarray, shape (3,)
            Command [s (-), f (-), V_vent (m3/s)].
        """
        T, g = y[IDX_Y_TA], self.guard
        if not g.on and self.T_set + 0.5 * self.band <= T and g.can_start(t):
            g.set(t, True)
        elif g.on and self.T_set - 0.5 * self.band >= T and g.can_stop(t):
            g.set(t, False)
        s = self.s_on if g.on else 0.0
        return np.array([s, fan_command(s, self.fan, self.fan_mode), self.ventilation(t)])
