"""Disturbance inputs d(t): outdoor weather, solar irradiance and occupancy.

Time t is in seconds from midnight of day 0. Negative times are allowed (warm-up days)
and every scenario is vectorised, so ``scenario(t_array)`` returns a forecast over a
horizon, which the predictive controllers will use later.

Solar input is the irradiance on the window plane, W/m2. The window properties
(area, g-value) live in ``RoomParams``, so the plant converts irradiance to a heat gain.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple, Protocol

import numpy as np

from . import psychro

DAY = 86_400.0
HOUR = 3_600.0


class Disturbances(NamedTuple):
    T_o: np.ndarray | float      # outdoor dry-bulb temperature, degC
    W_o: np.ndarray | float      # outdoor humidity ratio, kg/kg
    G_sol: np.ndarray | float    # solar irradiance on the window, W/m2
    n_occ: np.ndarray | float    # number of occupants, -


class Scenario(Protocol):
    def __call__(self, t) -> Disturbances: ...


def hour_of_day(t):
    """Clock hour in [0, 24) for time t in seconds."""
    return np.mod(np.asarray(t, dtype=float), DAY) / HOUR


@dataclass(frozen=True)
class SummerScenario:
    """Sinusoidal summer day repeated indefinitely.

    T_o   = T_mean + amp * sin(2*pi*(t/day - phase_h/24)); phase 9 h puts the maximum at 15 h.
    RH_o  constant, converted to W_o at the instantaneous T_o. Because p_ws is strongly
          nonlinear, W_o then doubles between 3 h and 15 h (0.008 -> 0.0165 kg/kg). Real
          air masses keep W roughly constant over a day while RH swings, so
          ``constant_W=True`` evaluates RH_o at T_mean instead (W_o = 0.0116 kg/kg).
    G_sol half-sine between sunrise and sunset, zero at night.
    n_occ step function, occupied on [occ_start, occ_end).
    """

    T_mean: float = 26.0
    T_amp: float = 6.0
    phase_h: float = 9.0
    RH_o: float = 0.55
    constant_W: bool = False
    G_max: float = 500.0
    sunrise_h: float = 6.0
    sunset_h: float = 18.0
    n_people: float = 2.0
    occ_start_h: float = 8.0
    occ_end_h: float = 18.0

    def T_o(self, t):
        t = np.asarray(t, dtype=float)
        return self.T_mean + self.T_amp * np.sin(2 * np.pi * (t / DAY - self.phase_h / 24.0))

    def W_o(self, t):
        if self.constant_W:
            ones = np.ones_like(np.asarray(t, dtype=float))
            return psychro.W_from_RH(self.T_mean, self.RH_o) * ones
        return psychro.W_from_RH(self.T_o(t), self.RH_o)

    def G_sol(self, t):
        h = hour_of_day(t)
        day_len = self.sunset_h - self.sunrise_h
        up = (h >= self.sunrise_h) & (h < self.sunset_h)
        return np.where(up, self.G_max * np.sin(np.pi * (h - self.sunrise_h) / day_len), 0.0)

    def n_occ(self, t):
        h = hour_of_day(t)
        occupied = (h >= self.occ_start_h) & (h < self.occ_end_h)
        return np.where(occupied, self.n_people, 0.0)

    def __call__(self, t) -> Disturbances:
        return Disturbances(self.T_o(t), self.W_o(t), self.G_sol(t), self.n_occ(t))


@dataclass(frozen=True)
class ConstantInputs:
    """Time-invariant disturbances, used by the verification tests.

    A step in outdoor conditions at t = 0 is modelled as constant inputs together with
    an initial state that differs from the new equilibrium.
    """

    T_o: float = 26.0
    W_o: float = 0.010
    G_sol: float = 0.0
    n_occ: float = 0.0

    def __call__(self, t) -> Disturbances:
        ones = np.ones_like(np.asarray(t, dtype=float))
        return Disturbances(self.T_o * ones, self.W_o * ones, self.G_sol * ones, self.n_occ * ones)
