"""Electricity price signals."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._typing import BoolArray, FloatArray, FloatOrArray
from .inputs import hour_of_day


def check_daily_window(start_name: str, start_h: float, end_name: str, end_h: float) -> None:
    """Check that [start_h, end_h) is a well-defined window of clock hours.

    A window with start_h > end_h crosses midnight. start_h = end_h is rejected because it is
    ambiguous (empty or the whole day).

    Parameters
    ----------
    start_name, end_name : str
        Names of the two limits, used in the error message.
    start_h, end_h : float
        Window limits, clock hours (h).

    Raises
    ------
    ValueError
        If a limit is outside [0, 24] or the two limits are equal.
    """
    for name, h in ((start_name, start_h), (end_name, end_h)):
        if not 0.0 <= h <= 24.0:
            raise ValueError(f"{name} = {h!r} is invalid; it must be a clock hour in [0, 24]")
    if start_h == end_h:
        raise ValueError(f"{start_name} = {start_h!r} equals {end_name} = {end_h!r}; the window "
                         f"[{start_name}, {end_name}) would be empty")


def in_daily_window(h: FloatOrArray, start_h: float, end_h: float) -> BoolArray:
    """Return where the clock hour h lies in [start_h, end_h); elementwise.

    With start_h > end_h the window crosses midnight: [start_h, 24) together with [0, end_h).

    Parameters
    ----------
    h : float or ndarray
        Clock hour in [0, 24), h.
    start_h, end_h : float
        Window limits, clock hours (h).

    Returns
    -------
    ndarray of bool
        True inside the window (a numpy scalar for scalar input).
    """
    h = np.asarray(h, dtype=float)
    if start_h <= end_h:
        inside: BoolArray = (h >= start_h) & (h < end_h)
        return inside
    wrapped: BoolArray = (h >= start_h) | (h < end_h)
    return wrapped


@dataclass(frozen=True)
class TimeOfUseTariff:
    """Two-level tariff, EUR/kWh: ``peak`` on [peak_start_h, peak_end_h), ``offpeak`` otherwise.

    A peak period with peak_start_h > peak_end_h crosses midnight.

    Attributes
    ----------
    offpeak : float
        Off-peak price, EUR/kWh.
    peak : float
        Peak price, EUR/kWh.
    peak_start_h, peak_end_h : float
        Peak period [peak_start_h, peak_end_h) in clock hours, h.

    Raises
    ------
    ValueError
        If a limit is outside [0, 24] or the two limits are equal (see
        ``check_daily_window``).
    """

    offpeak: float = 0.15
    peak: float = 0.40
    peak_start_h: float = 14.0
    peak_end_h: float = 20.0

    def __post_init__(self) -> None:
        """Validate the peak period."""
        check_daily_window("peak_start_h", self.peak_start_h, "peak_end_h", self.peak_end_h)

    def __call__(self, t: FloatOrArray) -> FloatArray:
        """Return the electricity price at time t.

        Parameters
        ----------
        t : float or ndarray
            Time from midnight of day 0, s.

        Returns
        -------
        ndarray
            Price, EUR/kWh.
        """
        h = hour_of_day(t)
        peak = in_daily_window(h, self.peak_start_h, self.peak_end_h)
        return np.where(peak, self.peak, self.offpeak)

    @property
    def reference(self) -> float:
        """Price used to express other costs (comfort) in equivalent kWh, EUR/kWh."""
        return self.offpeak


FLAT = TimeOfUseTariff(offpeak=0.20, peak=0.20)
