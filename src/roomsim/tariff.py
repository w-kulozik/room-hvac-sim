"""Electricity price signals."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._typing import FloatArray, FloatOrArray
from .inputs import hour_of_day


@dataclass(frozen=True)
class TimeOfUseTariff:
    """Two-level tariff, EUR/kWh: ``peak`` on [peak_start_h, peak_end_h), ``offpeak`` otherwise.

    Attributes
    ----------
    offpeak : float
        Off-peak price, EUR/kWh.
    peak : float
        Peak price, EUR/kWh.
    peak_start_h, peak_end_h : float
        Peak period [peak_start_h, peak_end_h) in clock hours, h.
    """

    offpeak: float = 0.15
    peak: float = 0.40
    peak_start_h: float = 14.0
    peak_end_h: float = 20.0

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
        return np.where((h >= self.peak_start_h) & (h < self.peak_end_h), self.peak, self.offpeak)

    @property
    def reference(self) -> float:
        """Price used to express other costs (comfort) in equivalent kWh, EUR/kWh."""
        return self.offpeak


FLAT = TimeOfUseTariff(offpeak=0.20, peak=0.20)
