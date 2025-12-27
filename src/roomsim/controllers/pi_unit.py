"""PI temperature control of an inverter unit with a low-load handler.

Per sample: the PI computes a continuous demand v from the air-temperature error, the low-load
handler turns v into a compressor command the unit can follow, and the PI integrator is updated
with what was actually delivered (back-calculation) or frozen while the unit cycles.
"""

from __future__ import annotations

import numpy as np

from ..ac_unit import ACUnit
from ..params import RoomParams
from ..plant import IDX_Y_TA
from .base import CompressorGuard, NO_VENTILATION, VentilationSchedule, fan_command
from .lowload import HysteresisLowLoad, SigmaDeltaLowLoad, ripple_T_ref
from .pi import DiscretePI, room_pi_tuning


class PIUnitController:
    def __init__(self, pi: DiscretePI, lowload, T_set=25.0, fan: float = 1.0,
                 fan_mode: str = "auto", ventilation: VentilationSchedule = NO_VENTILATION):
        self.pi, self.lowload = pi, lowload
        self.T_set, self.fan, self.fan_mode, self.ventilation = T_set, fan, fan_mode, ventilation
        self.last_demand = 0.0

    @property
    def guard(self) -> CompressorGuard:
        return self.lowload.guard

    def set_point(self, t: float) -> float:
        """T_set may be a constant or a schedule (callable of time)."""
        return self.T_set(t) if callable(self.T_set) else self.T_set

    def step(self, t: float, y: np.ndarray) -> np.ndarray:
        T_set = self.set_point(t)
        e = y[IDX_Y_TA] - T_set
        v = self.pi.demand(e)
        s, u_applied, integrate = self.lowload.step(t, v, y[IDX_Y_TA], T_set)
        self.pi.update(e, v, u_applied, integrate)
        self.last_demand = v
        return np.array([s, fan_command(s, self.fan, self.fan_mode), self.ventilation(t)])


def make_pi_controller(kind: str, p: RoomParams, unit: ACUnit, Ts: float = 60.0,
                       T_set=25.0, min_on: float = 180.0, min_off: float = 180.0,
                       ripple_K: float = 0.5, **kwargs) -> PIUnitController:
    """PI with low-load handler 'hysteresis' (A) or 'sigma-delta' (C), SIMC-tuned.

    The sigma-delta time scale is chosen for an air-temperature ripple of ``ripple_K`` (peak to
    peak) with the sensible capacity at full speed in typical summer conditions.
    """
    Kp, Ti = room_pi_tuning(p, unit, Ts)
    guard = CompressorGuard(min_on, min_off)
    if kind == "hysteresis":
        lowload = HysteresisLowLoad(unit.p.s_min, guard)
    elif kind == "sigma-delta":
        Q_max = unit.operate(1.0, 1.0, 25.0, 0.0115, 32.0).Q_sens
        T_ref = ripple_T_ref(ripple_K, p.C_a, Q_max, unit.p.s_min)
        lowload = SigmaDeltaLowLoad(unit.p.s_min, guard, Ts, T_ref)
    else:
        raise ValueError(f"unknown low-load handler {kind!r}")
    return PIUnitController(DiscretePI(Kp, Ti, Ts), lowload, T_set=T_set, **kwargs)
