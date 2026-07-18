"""Mechanical ventilation fan: electrical power and energy."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .plant import IDX_V_VENT
from .results import SimResult


@dataclass(frozen=True)
class VentilationFan:
    """Fan affinity law around a design point: P = P_d (V / V_d)^3.

    Default: 0.03 m3/s (108 m3/h) at 30 W, a specific fan power of 1000 W/(m3/s) at design flow.
    """

    V_d: float = 0.03
    P_d: float = 30.0

    def power(self, V) -> np.ndarray:
        V = np.asarray(V, dtype=float)
        return self.P_d * (V / self.V_d) ** 3

    def energy_kWh(self, r: SimResult) -> float:
        """Fan energy over a result. The flow is held over each step (row k covers
        [t_k, t_k + dt)), so this sum is exact.
        """
        return float(np.sum(self.power(r.u[:-1, IDX_V_VENT])) * r.dt / 3.6e6)
