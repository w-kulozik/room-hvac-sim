"""Electricity price signals."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .inputs import hour_of_day


@dataclass(frozen=True)
class TimeOfUseTariff:
    """Two-level tariff, EUR/kWh: ``peak`` on [peak_start_h, peak_end_h), ``offpeak`` otherwise."""

    offpeak: float = 0.15
    peak: float = 0.40
    peak_start_h: float = 14.0
    peak_end_h: float = 20.0

    def __call__(self, t) -> np.ndarray:
        h = hour_of_day(t)
        return np.where((h >= self.peak_start_h) & (h < self.peak_end_h), self.peak, self.offpeak)

    @property
    def reference(self) -> float:
        """Price used to express other costs (comfort) in equivalent kWh."""
        return self.offpeak


FLAT = TimeOfUseTariff(offpeak=0.20, peak=0.20)
