"""PI temperature control of an inverter unit with a low-load handler.

Per sample: the PI computes a continuous demand v from the air-temperature error, the low-load
handler turns v into a compressor command the unit can follow, and the PI integrator is updated
with what was actually delivered (back-calculation) or frozen while the unit cycles.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .._typing import FloatArray
from ..ac_unit import ACUnit
from ..params import RoomParams
from ..plant import IDX_Y_TA
from .base import (
    NO_VENTILATION,
    CompressorGuard,
    SetPoint,
    VentilationSchedule,
    common_sampling_period,
    fan_command,
    sampling_period,
    set_point_at,
)
from .lowload import HysteresisLowLoad, LowLoadHandler, SigmaDeltaLowLoad, ripple_T_ref
from .pi import DiscretePI, room_pi_tuning


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
    heating : bool, optional
        Heating mode: the error is T_set - T, so a room colder than the set point raises the
        demand (D53). Only the sigma-delta handler (C) supports it.

    Attributes
    ----------
    last_demand : float
        Unsaturated PI demand of the last sample, -.

    Raises
    ------
    ValueError
        If the PI and the low-load handler assume different sampling periods, or heating is
        requested with the hysteresis handler (its switching rule is written for cooling).
    """

    def __init__(self, pi: DiscretePI, lowload: LowLoadHandler, T_set: SetPoint = 25.0,
                 fan: float = 1.0, fan_mode: str = "auto",
                 ventilation: VentilationSchedule = NO_VENTILATION,
                 heating: bool = False) -> None:
        common_sampling_period("PIUnitController", pi_Ts=pi.Ts,
                               lowload_Ts=sampling_period(lowload))
        if heating and isinstance(lowload, HysteresisLowLoad):
            raise ValueError("heating = True needs the sigma-delta low-load handler; the "
                             "hysteresis handler switches on cooling thresholds")
        self.heating = heating
        self.pi, self.lowload = pi, lowload
        self.T_set, self.fan, self.fan_mode, self.ventilation = T_set, fan, fan_mode, ventilation
        self.last_demand = 0.0

    @property
    def guard(self) -> CompressorGuard:
        """Compressor guard of the low-load handler."""
        return self.lowload.guard

    @property
    def Ts(self) -> float:
        """Sampling period the PI (and the low-load handler) assume, s."""
        return self.pi.Ts

    def set_point(self, t: float) -> float:
        """Return the set point at time t (s), degC; T_set may be a constant or a schedule."""
        return set_point_at(self.T_set, t)

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
        e = T_set - y[IDX_Y_TA] if self.heating else y[IDX_Y_TA] - T_set
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
