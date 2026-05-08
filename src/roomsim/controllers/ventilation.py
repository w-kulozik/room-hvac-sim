"""Ventilation strategies, combined with any temperature/humidity controller.

A strategy has ``step(t, y) -> V_vent`` (m3/s). ``WithVentilation`` runs the inner controller
and replaces the ventilation element of its command, so the same compressor and fan control is
used with every strategy.

    ScheduledVentilation  fixed flow during occupied hours
    DemandControlledPI    PI on measured CO2 (demand-controlled ventilation, continuous)
    StagedVentilation     discrete fan levels switched on CO2 with hysteresis and dwell time
    NightFlush            adds free cooling with outdoor air when the office is empty
"""

from __future__ import annotations

from typing import Callable

import numpy as np

from ..params import RoomParams
from ..plant import IDX_Y_CO2, IDX_Y_TA
from .base import VentilationSchedule
from .pi import DiscretePI, simc_pi
from .mpc import ComfortSchedule

IDX_VENT = 2
CO2_OUTDOOR = 420.0


class WithVentilation:
    def __init__(self, inner, strategy):
        self.inner, self.strategy = inner, strategy

    def step(self, t: float, y: np.ndarray) -> np.ndarray:
        cmd = np.array(self.inner.step(t, y), dtype=float)
        cmd[IDX_VENT] = self.strategy.step(t, y)
        return cmd


class ScheduledVentilation:
    def __init__(self, schedule: VentilationSchedule = VentilationSchedule()):
        self.schedule = schedule

    def step(self, t: float, y: np.ndarray) -> float:
        return self.schedule(t)


def dcv_tuning(p: RoomParams, C_set: float = 900.0, n_people: float = 2.0,
               tau_c: float = 900.0, Ts: float = 60.0) -> tuple[float, float]:
    """SIMC gains for the CO2 loop, linearised at the set point with n_people present.

    At the set point the total outdoor-air flow is Q = N G / (C_set - C_o); the process is first
    order with gain dC/dV = -(C_set - C_o)/Q (ppm per m3/s) and time constant V_room / Q.
    Returns Kp in m3/s per ppm (positive: more CO2, more air) and Ti in s.
    """
    Q = n_people * p.co2_person * 1e6 / (C_set - CO2_OUTDOOR)
    K = (C_set - CO2_OUTDOOR) / Q
    tau = p.volume / Q
    return simc_pi(K, tau, Ts / 2.0, tau_c)


class DemandControlledPI:
    """Continuous demand-controlled ventilation: PI on CO2 within [V_min, V_max] when occupied."""

    def __init__(self, p: RoomParams, C_set: float = 900.0, V_min: float = 0.006,
                 V_max: float = 0.03, comfort: ComfortSchedule = ComfortSchedule(),
                 Ts: float = 60.0):
        Kp, Ti = dcv_tuning(p, C_set, Ts=Ts)
        self.pi = DiscretePI(Kp, Ti, Ts)
        self.C_set, self.V_min, self.V_max, self.comfort = C_set, V_min, V_max, comfort

    def step(self, t: float, y: np.ndarray) -> float:
        if not self.comfort.occupied(t):
            return 0.0                                     # integrator kept for the next day
        e = y[IDX_Y_CO2] - self.C_set
        v = self.pi.demand(e)
        V = min(max(v, self.V_min), self.V_max)
        self.pi.update(e, v, V)
        return V


class StagedVentilation:
    """Discrete fan levels switched on CO2, with hysteresis and a minimum dwell time per level."""

    def __init__(self, levels=(0.006, 0.012, 0.024), up=(850.0, 950.0), down=(750.0, 850.0),
                 dwell: float = 300.0, comfort: ComfortSchedule = ComfortSchedule()):
        if not (len(up) == len(down) == len(levels) - 1):
            raise ValueError("need one up and one down threshold between consecutive levels")
        self.levels, self.up, self.down = levels, up, down
        self.dwell, self.comfort = dwell, comfort
        self.stage, self.t_switch = 0, -np.inf

    def step(self, t: float, y: np.ndarray) -> float:
        if not self.comfort.occupied(t):
            self.stage, self.t_switch = 0, -np.inf
            return 0.0
        C, i = y[IDX_Y_CO2], self.stage
        if t - self.t_switch >= self.dwell:
            if i < len(self.levels) - 1 and C >= self.up[i]:
                self.stage, self.t_switch = i + 1, t
            elif i > 0 and C <= self.down[i - 1]:
                self.stage, self.t_switch = i - 1, t
        return self.levels[self.stage]


class NightFlush:
    """Adds free cooling while the office is empty: full flow when the outdoor air is at least
    dT cooler than the room and the room is above T_min. ``outdoor`` is the outdoor sensor."""

    def __init__(self, inner, outdoor: Callable[[float], float], V_flush: float = 0.03,
                 dT: float = 2.0, T_min: float = 22.0, comfort: ComfortSchedule = ComfortSchedule()):
        self.inner, self.outdoor = inner, outdoor
        self.V_flush, self.dT, self.T_min, self.comfort = V_flush, dT, T_min, comfort

    def step(self, t: float, y: np.ndarray) -> float:
        V = self.inner.step(t, y)
        if self.comfort.occupied(t):
            return V
        T_a, T_o = y[IDX_Y_TA], float(self.outdoor(t))
        if T_o <= T_a - self.dT and T_a > self.T_min:
            return self.V_flush
        return V
