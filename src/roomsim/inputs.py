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
from ._typing import BoolArray, FloatArray, FloatOrArray

DAY = 86_400.0
HOUR = 3_600.0

# Disturbances that are piecewise constant (they jump on grid points). The simulator samples
# them once per step at the step start instead of at every Runge-Kutta stage; otherwise the
# last RK4 stage of the step ending at 8:00 would already see the 8:00 value.
PIECEWISE_CONSTANT = ("n_occ",)


class Disturbances(NamedTuple):
    """Disturbance values at one time or over an array of times.

    Attributes
    ----------
    T_o : float or ndarray
        Outdoor dry-bulb temperature, degC.
    W_o : float or ndarray
        Outdoor humidity ratio, kg/kg dry air.
    G_sol : float or ndarray
        Solar irradiance on the window plane, W/m2.
    n_occ : float or ndarray
        Number of occupants, -.
    C_o : float or ndarray
        Outdoor CO2 concentration, ppm.
    """

    T_o: FloatArray | float      # outdoor dry-bulb temperature, degC
    W_o: FloatArray | float      # outdoor humidity ratio, kg/kg
    G_sol: FloatArray | float    # solar irradiance on the window, W/m2
    n_occ: FloatArray | float    # number of occupants, -
    C_o: FloatArray | float      # outdoor CO2 concentration, ppm


class Scenario(Protocol):
    """Callable that returns the disturbances at time t."""

    def __call__(self, t: FloatOrArray) -> Disturbances:
        """Evaluate the disturbances.

        Parameters
        ----------
        t : float or ndarray
            Time from midnight of day 0, s.

        Returns
        -------
        Disturbances
            Values at t, with the shape of t.
        """
        ...


def _in_interval(h: FloatArray, start: float, end: float) -> BoolArray:
    """Return where the clock hour lies in the interval [start, end).

    An interval with start > end wraps around midnight: [start, 24) together with [0, end).

    Parameters
    ----------
    h : ndarray
        Clock hour, h.
    start, end : float
        Interval limits, h.

    Returns
    -------
    ndarray of bool
        True inside the interval.
    """
    if start <= end:
        inside: BoolArray = (h >= start) & (h < end)
    else:
        inside = (h >= start) | (h < end)
    return inside


def hour_of_day(t: FloatOrArray) -> FloatArray:
    """Compute the clock hour in [0, 24).

    Parameters
    ----------
    t : float or ndarray
        Time from midnight of day 0, s (negative values allowed).

    Returns
    -------
    ndarray
        Clock hour, h (a numpy scalar for scalar input).
    """
    h: FloatArray = np.mod(np.asarray(t, dtype=float), DAY) / HOUR
    return h


def occupants(t: FloatOrArray, n_people: float, start_h: float, end_h: float,
              occupancy: tuple[tuple[float, float, float], ...] | None = None) -> FloatArray:
    """Compute the number of occupants of a daily schedule.

    Parameters
    ----------
    t : float or ndarray
        Time from midnight of day 0, s.
    n_people : float
        Occupants during [start_h, end_h), -.
    start_h, end_h : float
        Occupied period, clock hours (h); start_h > end_h crosses midnight.
    occupancy : tuple of (float, float, float), optional
        Intervals (start_h, end_h, people), h, h, -. If given, it replaces n_people and the
        occupied period; a later interval overrides an earlier one where they overlap.

    Returns
    -------
    ndarray
        Number of occupants, -.
    """
    h = hour_of_day(t)
    if occupancy is None:
        return np.where(_in_interval(h, start_h, end_h), n_people, 0.0)
    n = np.zeros_like(h)
    for start, end, people in occupancy:
        n = np.where(_in_interval(h, start, end), people, n)
    return n


@dataclass(frozen=True)
class SummerScenario:
    """Sinusoidal summer day repeated indefinitely.

    T_o   = T_mean + amp * sin(2*pi*(t/day - phase_h/24)); phase 9 h puts the maximum at 15 h.
    W_o   by default RH_o evaluated at T_mean and held constant (W_o = 0.0115 kg/kg), so the
          outdoor RH swings during the day as it does in a real air mass. With
          ``constant_W=False`` RH_o is held constant instead and converted at the instantaneous
          T_o; because p_ws is strongly nonlinear, W_o then doubles between 3 h and 15 h
          (0.008 -> 0.0165 kg/kg), which overstates the afternoon latent load.
    G_sol half-sine between sunrise and sunset, zero at night.
    n_occ step function, occupied on [occ_start, occ_end) with n_people; or, if ``occupancy`` is
          given, a profile of (start_h, end_h, people) intervals (zero outside them). An
          interval with start > end crosses midnight (e.g. 22 h to 6 h).

    Attributes
    ----------
    T_mean : float
        Daily mean outdoor temperature, degC.
    T_amp : float
        Amplitude of the outdoor-temperature sine, K.
    phase_h : float
        Phase of the outdoor-temperature sine, h.
    RH_o : float
        Outdoor relative humidity, fraction 0..1.
    constant_W : bool
        If True, W_o is held constant at its value for (T_mean, RH_o); otherwise RH_o is held.
    G_max : float
        Peak solar irradiance on the window, W/m2.
    sunrise_h, sunset_h : float
        Start and end of the solar half-sine, h.
    n_people : float
        Occupants during the occupied period, -.
    occ_start_h, occ_end_h : float
        Occupied period [occ_start_h, occ_end_h), h.
    C_o : float
        Outdoor CO2 concentration, ppm.
    occupancy : tuple of (float, float, float), optional
        Intervals (start_h, end_h, people), h, h, -. If given, it replaces n_people and the
        occupied period; a later interval overrides an earlier one where they overlap.

    Raises
    ------
    ValueError
        If RH_o is outside [0, 1] or sunset_h <= sunrise_h.
    """

    T_mean: float = 26.0
    T_amp: float = 6.0
    phase_h: float = 9.0
    RH_o: float = 0.55
    constant_W: bool = True
    G_max: float = 500.0
    sunrise_h: float = 6.0
    sunset_h: float = 18.0
    n_people: float = 2.0
    occ_start_h: float = 8.0
    occ_end_h: float = 18.0
    C_o: float = 420.0
    occupancy: tuple[tuple[float, float, float], ...] | None = None

    def __post_init__(self) -> None:
        """Reject a relative humidity outside 0..1 and a non-positive day length."""
        if not 0.0 <= self.RH_o <= 1.0:
            raise ValueError(f"RH_o = {self.RH_o!r} is invalid; it must be a fraction in [0, 1]")
        if not self.sunset_h > self.sunrise_h:
            raise ValueError(f"sunset_h = {self.sunset_h!r} is invalid; it must be greater than "
                             f"sunrise_h = {self.sunrise_h!r}")

    def T_o(self, t: FloatOrArray) -> FloatArray:
        """Compute the outdoor dry-bulb temperature.

        Parameters
        ----------
        t : float or ndarray
            Time from midnight of day 0, s.

        Returns
        -------
        ndarray
            Outdoor temperature, degC.
        """
        t = np.asarray(t, dtype=float)
        return self.T_mean + self.T_amp * np.sin(2 * np.pi * (t / DAY - self.phase_h / 24.0))

    def W_o(self, t: FloatOrArray) -> FloatArray:
        """Compute the outdoor humidity ratio.

        Parameters
        ----------
        t : float or ndarray
            Time from midnight of day 0, s.

        Returns
        -------
        ndarray
            Outdoor humidity ratio, kg/kg dry air.
        """
        if self.constant_W:
            ones = np.ones_like(np.asarray(t, dtype=float))
            W: FloatArray = psychro.W_from_RH(self.T_mean, self.RH_o) * ones
            return W
        return psychro.W_from_RH(self.T_o(t), self.RH_o)

    def G_sol(self, t: FloatOrArray) -> FloatArray:
        """Compute the solar irradiance on the window.

        Parameters
        ----------
        t : float or ndarray
            Time from midnight of day 0, s.

        Returns
        -------
        ndarray
            Irradiance, W/m2 (half-sine between sunrise and sunset, zero otherwise).
        """
        h = hour_of_day(t)
        day_len = self.sunset_h - self.sunrise_h
        up = (h >= self.sunrise_h) & (h < self.sunset_h)
        G: FloatArray = np.where(up, self.G_max * np.sin(np.pi * (h - self.sunrise_h) / day_len),
                                 0.0)
        return G

    def n_occ(self, t: FloatOrArray) -> FloatArray:
        """Compute the number of occupants.

        Parameters
        ----------
        t : float or ndarray
            Time from midnight of day 0, s.

        Returns
        -------
        ndarray
            Number of occupants, -.
        """
        return occupants(t, self.n_people, self.occ_start_h, self.occ_end_h, self.occupancy)

    def __call__(self, t: FloatOrArray) -> Disturbances:
        """Evaluate all disturbances.

        Parameters
        ----------
        t : float or ndarray
            Time from midnight of day 0, s.

        Returns
        -------
        Disturbances
            Values at t.
        """
        ones = np.ones_like(np.asarray(t, dtype=float))
        return Disturbances(self.T_o(t), self.W_o(t), self.G_sol(t), self.n_occ(t), self.C_o * ones)


@dataclass(frozen=True)
class ConstantInputs:
    """Time-invariant disturbances, used by the verification tests.

    A step in outdoor conditions at t = 0 is modelled as constant inputs together with
    an initial state that differs from the new equilibrium.

    Attributes
    ----------
    T_o : float
        Outdoor dry-bulb temperature, degC.
    W_o : float
        Outdoor humidity ratio, kg/kg dry air.
    G_sol : float
        Solar irradiance on the window, W/m2.
    n_occ : float
        Number of occupants, -.
    C_o : float
        Outdoor CO2 concentration, ppm.
    """

    T_o: float = 26.0
    W_o: float = 0.010
    G_sol: float = 0.0
    n_occ: float = 0.0
    C_o: float = 420.0

    def __call__(self, t: FloatOrArray) -> Disturbances:
        """Evaluate the disturbances.

        Parameters
        ----------
        t : float or ndarray
            Time from midnight of day 0, s.

        Returns
        -------
        Disturbances
            The constant values broadcast to the shape of t.
        """
        ones = np.ones_like(np.asarray(t, dtype=float))
        return Disturbances(self.T_o * ones, self.W_o * ones, self.G_sol * ones,
                            self.n_occ * ones, self.C_o * ones)


# Office day with variable occupancy (decision D27): late arrival, lunch, a meeting, early leave.
VARIABLE_OCCUPANCY = ((8.0, 9.0, 1.0), (9.0, 12.0, 2.0), (13.0, 15.0, 2.0), (15.0, 16.0, 3.0),
                      (16.0, 18.0, 1.0))
