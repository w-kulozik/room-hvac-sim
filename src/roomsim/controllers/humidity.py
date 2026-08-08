"""Humidity-aware control of an inverter unit through the indoor fan.

The compressor stays under the temperature PI (with its low-load handler). The indoor fan is
the dehumidification actuator: less airflow means a colder coil surface, more condensation and a
lower sensible heat ratio, at a small cost in efficiency and capacity. At minimum compressor speed,
going from full to 40 % fan takes the mean coil surface from about 17.8 degC (dry) to 14.5 degC
(wet), SHR 1.00 -> 0.70, with 10 % less capacity (room at 25 degC and 60 % RH, 32 degC outdoors).

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

from .._typing import FloatArray
from ..plant import IDX_Y_RH, IDX_Y_TCOIL
from .base import common_sampling_period, sampling_period
from .pi import DiscretePI
from .pi_unit import PIUnitController

IDX_S, IDX_F = 0, 1


class FanHumidityPI:
    """H1: relative-humidity PI on the indoor fan speed.

    Parameters
    ----------
    temperature : PIUnitController
        Temperature controller that sets the compressor command.
    RH_set : float, optional
        Relative-humidity set point, - (fraction).
    Kp : float, optional
        Proportional gain, fan reduction demand per %RH.
    Ti : float, optional
        Integral time, s.
    f_min : float, optional
        Minimum fan speed, -.
    Ts : float, optional
        Sampling period, s.

    Raises
    ------
    ValueError
        If ``Ts`` differs from the sampling period of ``temperature``.
    """

    def __init__(self, temperature: PIUnitController, RH_set: float = 0.55,
                 Kp: float = 0.05, Ti: float = 1800.0, f_min: float = 0.4,
                 Ts: float = 60.0) -> None:
        common_sampling_period("FanHumidityPI", Ts=Ts,
                               temperature_Ts=sampling_period(temperature))
        self.Ts = Ts
        self.temperature, self.RH_set, self.f_min = temperature, RH_set, f_min
        self.pi = DiscretePI(Kp, Ti, Ts)      # demand per %RH

    def step(self, t: float, y: FloatArray) -> FloatArray:
        """Return the command for one sample, with the fan set by the humidity PI.

        Parameters
        ----------
        t : float
            Current time, s.
        y : ndarray
            Measurements; uses the air temperature (degC) and the relative humidity
            y[IDX_Y_RH] (-).

        Returns
        -------
        ndarray, shape (3,)
            Command [s (-), f (-), V_vent (m3/s)].
        """
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
    """H2: RH (outer) -> coil temperature set point -> fan speed (inner).

    Parameters
    ----------
    temperature : PIUnitController
        Temperature controller that sets the compressor command.
    RH_set : float, optional
        Relative-humidity set point, - (fraction).
    Kp : float, optional
        Outer proportional gain, K of coil set-point depression per %RH.
    Ti : float, optional
        Outer integral time, s.
    T_s_max, T_s_min : float, optional
        Limits of the coil temperature set point, degC.
    Ki_fan : float, optional
        Inner integral gain, fan speed change per K of coil error per sample, 1/K.
    f_min : float, optional
        Minimum fan speed, -.
    settle_samples : int, optional
        Samples after a compressor start before the inner loop acts.
    Ts : float, optional
        Sampling period, s (``Ki_fan`` and ``settle_samples`` are per sample).

    Raises
    ------
    ValueError
        If ``Ts`` differs from the sampling period of ``temperature``.

    Attributes
    ----------
    f : float
        Fan speed held by the inner loop, -.
    T_s_set : float
        Current coil temperature set point, degC.
    samples_on : int
        Consecutive samples with the compressor running.
    """

    def __init__(self, temperature: PIUnitController, RH_set: float = 0.55,
                 Kp: float = 0.3, Ti: float = 1800.0, T_s_max: float = 20.0, T_s_min: float = 8.0,
                 Ki_fan: float = 0.04, f_min: float = 0.4, settle_samples: int = 2,
                 Ts: float = 60.0) -> None:
        common_sampling_period("CoilCascade", Ts=Ts, temperature_Ts=sampling_period(temperature))
        self.Ts = Ts
        self.temperature, self.RH_set = temperature, RH_set
        self.T_s_max, self.T_s_min = T_s_max, T_s_min
        self.Ki_fan, self.f_min, self.settle = Ki_fan, f_min, settle_samples
        self.outer = DiscretePI(Kp, Ti, Ts)   # coil set-point depression, K per %RH
        self.f = 1.0
        self.T_s_set = T_s_max
        self.samples_on = 0

    def step(self, t: float, y: FloatArray) -> FloatArray:
        """Return the command for one sample, with the fan set by the cascade.

        Parameters
        ----------
        t : float
            Current time, s.
        y : ndarray
            Measurements; uses the air temperature (degC), the relative humidity y[IDX_Y_RH]
            (-) and the coil temperature y[IDX_Y_TCOIL] (degC).

        Returns
        -------
        ndarray, shape (3,)
            Command [s (-), f (-), V_vent (m3/s)].
        """
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
