"""PI temperature control of an inverter unit with a low-load handler.

Per sample: the PI computes a continuous demand v from the air-temperature error, the low-load
handler turns v into a compressor command the unit can follow, and the PI integrator is updated
with what was actually delivered (back-calculation) or frozen while the unit cycles.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np

from .._typing import FloatArray
from ..ac_unit import ACUnit
from ..params import RoomParams
from ..plant import IDX_Y_TA
from .base import NO_VENTILATION, CompressorGuard, VentilationSchedule, fan_command
from .lowload import HysteresisLowLoad, LowLoadHandler, SigmaDeltaLowLoad, ripple_T_ref
from .pi import DiscretePI, room_pi_tuning

SetPoint = float | Callable[[float], float]
"""Constant set point (degC) or a schedule of time (s) returning degC."""


class PIUnitController:
    """Air-temperature PI acting on the compressor through a low-load handler.

    Parameters
    ----------
    pi : DiscretePI
        Temperature PI (demand per K).
    lowload : LowLoadHandler
        Converts the PI demand into a compressor command; holds the compressor guard.
    T_set : float or callable, optional
        Set point, degC, or a function of time (s) returning it.
    fan : float, optional
        Indoor-fan speed, -.
    fan_mode : str, optional
        'auto' or 'continuous' (see ``fan_command``).
    ventilation : VentilationSchedule, optional
        Outdoor-air flow schedule.

    Attributes
    ----------
    last_demand : float
        Unsaturated PI demand of the last sample, -.
    """

    def __init__(self, pi: DiscretePI, lowload: LowLoadHandler, T_set: SetPoint = 25.0,
                 fan: float = 1.0, fan_mode: str = "auto",
                 ventilation: VentilationSchedule = NO_VENTILATION) -> None:
        self.pi, self.lowload = pi, lowload
        self.T_set, self.fan, self.fan_mode, self.ventilation = T_set, fan, fan_mode, ventilation
        self.last_demand = 0.0

    @property
    def guard(self) -> CompressorGuard:
        """Compressor guard of the low-load handler."""
        return self.lowload.guard

    def set_point(self, t: float) -> float:
        """Return the set point at time t (s), degC; T_set may be a constant or a schedule."""
        return self.T_set(t) if callable(self.T_set) else self.T_set

    def step(self, t: float, y: FloatArray) -> FloatArray:
        """Return the command for one controller sample and update the PI.

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
        T_set = self.set_point(t)
        e = y[IDX_Y_TA] - T_set
        v = self.pi.demand(e)
        s, u_applied, integrate = self.lowload.step(t, v, y[IDX_Y_TA], T_set)
        self.pi.update(e, v, u_applied, integrate)
        self.last_demand = v
        return np.array([s, fan_command(s, self.fan, self.fan_mode), self.ventilation(t)])


def make_pi_controller(kind: str, p: RoomParams, unit: ACUnit, Ts: float = 60.0,
                       T_set: SetPoint = 25.0, min_on: float = 180.0, min_off: float = 180.0,
                       ripple_K: float = 0.5, **kwargs: Any) -> PIUnitController:
    """Build a SIMC-tuned PI with low-load handler 'hysteresis' (A) or 'sigma-delta' (C).

    The sigma-delta time scale is chosen for an air-temperature ripple of ``ripple_K`` (peak to
    peak). While cycling the unit runs at s_min with a dry or partly wet coil, so the relevant
    capacity is the sensible output per unit of demand at minimum speed, not at full speed (where
    the coil is wet and part of the capacity is latent): using the latter underestimated the
    ripple by 1.5 times.

    Parameters
    ----------
    kind : str
        'hysteresis' or 'sigma-delta'.
    p : RoomParams
        Room parameters used for tuning.
    unit : ACUnit
        Air-conditioning unit.
    Ts : float, optional
        Sampling period, s.
    T_set : float or callable, optional
        Set point, degC, or a function of time (s) returning it.
    min_on, min_off : float, optional
        Minimum run and stop times of the compressor, s.
    ripple_K : float, optional
        Target air-temperature ripple for the sigma-delta handler, peak to peak, K.
    **kwargs
        Passed to ``PIUnitController`` (``fan``, ``fan_mode``, ``ventilation``).

    Returns
    -------
    PIUnitController
        The configured controller.

    Raises
    ------
    ValueError
        If ``kind`` is not recognised.
    """
    Kp, Ti = room_pi_tuning(p, unit, Ts)
    guard = CompressorGuard(min_on, min_off)
    lowload: LowLoadHandler
    if kind == "hysteresis":
        lowload = HysteresisLowLoad(unit.p.s_min, guard)
    elif kind == "sigma-delta":
        s_min = unit.p.s_min
        Q_per_demand = unit.operate(s_min, 1.0, 25.0, 0.0115, 32.0).Q_sens / s_min
        T_ref = ripple_T_ref(ripple_K, p.C_a, Q_per_demand, s_min)
        lowload = SigmaDeltaLowLoad(unit.p.s_min, guard, Ts, T_ref)
    else:
        raise ValueError(f"unknown low-load handler {kind!r}")
    return PIUnitController(DiscretePI(Kp, Ti, Ts), lowload, T_set=T_set, **kwargs)
