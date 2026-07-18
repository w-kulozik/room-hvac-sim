"""Linear two-node thermal model used inside the controllers (not the simulated plant).

State   z = [T_a, T_m]
Input   Q_c   sensible heat removed from the air by the unit, W (>= 0)
Known   e = [T_o, Q_air, Q_mass]   outdoor temperature and the heat gains to each node, W
        (solar and occupant gains, split as in the plant)
Ventilation adds H_vent = rho cp V_vent to the air-outdoor conductance, so the model is
piecewise time-invariant: one matrix pair per ventilation level.

Continuous form  dz/dt = A z + B Q_c + E e
    A = [[-(H_am + H_ao + H_inf + H_vent), H_am], [H_am, -(H_am + H_mo)]] / diag(C_a, C_m)
    B = [-1/C_a, 0]
    E = [[(H_ao + H_inf + H_vent)/C_a, 1/C_a, 0], [H_mo/C_m, 0, 1/C_m]]
Discretised exactly for inputs held over the step (zero-order hold) through the matrix
exponential of the augmented matrix [[A, [B E]], [0, 0]] * Ts.

The controller's model can be built from mismatched parameters (``scale``) to test robustness.
It knows nothing about humidity, the coil or the unit's efficiency; the Kalman filter's
disturbance state absorbs what the model gets wrong.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import lru_cache

import numpy as np
from scipy.linalg import expm

from .inputs import Disturbances
from .params import RoomParams


@dataclass(frozen=True)
class RCModel:
    p: RoomParams = field(default_factory=RoomParams)

    @classmethod
    def mismatched(cls, p: RoomParams, scale: dict[str, float]) -> RCModel:
        """Model built from parameters multiplied by the given factors, e.g. {'C_m': 1.2}."""
        return cls(replace(p, **{k: getattr(p, k) * v for k, v in scale.items()}))

    def continuous(self, V_vent: float = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        p = self.p
        H_out = p.H_ao + p.H_inf + p.rho_air * p.cp_air * V_vent
        A = np.array([[-(p.H_am + H_out) / p.C_a, p.H_am / p.C_a],
                      [p.H_am / p.C_m, -(p.H_am + p.H_mo) / p.C_m]])
        B = np.array([-1.0 / p.C_a, 0.0])
        E = np.array([[H_out / p.C_a, 1.0 / p.C_a, 0.0],
                      [p.H_mo / p.C_m, 0.0, 1.0 / p.C_m]])
        return A, B, E

    def discrete(self, Ts: float, V_vent: float = 0.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return _discretise(self, float(Ts), float(V_vent))

    def gains(self, d: Disturbances) -> np.ndarray:
        """Known inputs e = [T_o, Q_air, Q_mass] from (forecast) disturbances; rows = time."""
        p = self.p
        q_sol = p.g_win * p.A_win * np.asarray(d.G_sol, dtype=float)
        q_occ = p.q_person_sens * np.asarray(d.n_occ, dtype=float)
        Q_air = p.f_sol_air * q_sol + p.f_person_conv * q_occ
        Q_mass = (1 - p.f_sol_air) * q_sol + (1 - p.f_person_conv) * q_occ
        return np.column_stack([np.asarray(d.T_o, dtype=float) * np.ones_like(Q_air), Q_air, Q_mass])


@lru_cache(maxsize=64)
def _discretise(model: RCModel, Ts: float, V_vent: float):
    A, B, E = model.continuous(V_vent)
    n, m = A.shape[0], 1 + E.shape[1]
    M = np.zeros((n + m, n + m))
    M[:n, :n] = A
    M[:n, n] = B
    M[:n, n + 1:] = E
    Phi = expm(M * Ts)
    return Phi[:n, :n], Phi[:n, n].copy(), Phi[:n, n + 1:].copy()
