"""Mechanical ventilation fan: electrical power and energy."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ._typing import FloatArray, FloatOrArray
from .plant import IDX_V_VENT
from .results import SimResult


@dataclass(frozen=True)
class VentilationFan:
    """Fan affinity law around a design point: P = P_d (V / V_d)^3.

    Default: 0.03 m3/s (108 m3/h) at 30 W, a specific fan power of 1000 W/(m3/s) at design flow.

    Attributes
    ----------
    V_d : float
        Design volume flow, m3/s.
    P_d : float
        Electrical power at the design flow, W.
    """

    V_d: float = 0.03
    P_d: float = 30.0

    def power(self, V: FloatOrArray) -> FloatArray:
        """Compute the fan electrical power.

        Parameters
        ----------
        V : float or ndarray
            Ventilation volume flow, m3/s.

        Returns
        -------
        ndarray
            Electrical power, W.

        Raises
        ------
        ValueError
            If any flow is negative (the affinity law would give a negative power).
        """
        V = np.asarray(V, dtype=float)
        if np.any(V < 0.0):
            raise ValueError(f"V = {np.min(V)!r} m3/s is invalid; the flow must be >= 0")
        P: FloatArray = self.P_d * (V / self.V_d) ** 3
        return P

    def energy_kWh(self, r: SimResult) -> float:
        """Compute the fan energy over a simulation result.

        The flow is held over each step (row k covers [t_k, t_k + dt)), so this sum is exact.

        Parameters
        ----------
        r : SimResult
            Simulation result; the ventilation flow is read from its input log.

        Returns
        -------
        float
            Fan electrical energy, kWh.
        """
        return float(np.sum(self.power(r.u[:-1, IDX_V_VENT])) * r.dt / 3.6e6)
