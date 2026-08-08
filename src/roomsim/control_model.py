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

from ._typing import FloatArray
from .inputs import Disturbances
from .params import RoomParams

"""Array of float64 values (alias used for annotations in the control code)."""


@dataclass(frozen=True)
class RCModel:
    """Two-node (air, mass) linear thermal model of the room used for prediction and estimation.

    Attributes
    ----------
    p : RoomParams
        Room parameters the model is built from (capacities C_a, C_m in J/K, conductances
        H_am, H_ao, H_inf, H_mo in W/K, gain split factors).
    """

    p: RoomParams = field(default_factory=RoomParams)

    @classmethod
    def mismatched(cls, p: RoomParams, scale: dict[str, float]) -> RCModel:
        """Build a model from parameters multiplied by the given factors.

        Parameters
        ----------
        p : RoomParams
            Nominal room parameters.
        scale : dict of str to float
            Multiplicative factor per ``RoomParams`` field name, e.g. ``{'C_m': 1.2}``, -.

        Returns
        -------
        RCModel
            Model whose parameters are ``p`` with the listed fields scaled.
        """
        return cls(replace(p, **{k: getattr(p, k) * v for k, v in scale.items()}))

    def continuous(self, V_vent: float = 0.0) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Return the continuous-time matrices of dz/dt = A z + B Q_c + E e.

        Parameters
        ----------
        V_vent : float, optional
            Mechanical ventilation flow, m3/s; adds rho cp V_vent to the air-outdoor
            conductance.

        Returns
        -------
        A : ndarray, shape (2, 2)
            State matrix, 1/s.
        B : ndarray, shape (2,)
            Input vector for the sensible cooling Q_c, K/J.
        E : ndarray, shape (2, 3)
            Matrix for the known inputs e = [T_o, Q_air, Q_mass]; first column in 1/s, the
            others in K/J.
        """
        p = self.p
        H_out = p.H_ao + p.H_inf + p.rho_air * p.cp_air * V_vent
        A = np.array([[-(p.H_am + H_out) / p.C_a, p.H_am / p.C_a],
                      [p.H_am / p.C_m, -(p.H_am + p.H_mo) / p.C_m]])
        B = np.array([-1.0 / p.C_a, 0.0])
        E = np.array([[H_out / p.C_a, 1.0 / p.C_a, 0.0],
                      [p.H_mo / p.C_m, 0.0, 1.0 / p.C_m]])
        return A, B, E

    def discrete(self, Ts: float, V_vent: float = 0.0) -> tuple[FloatArray, FloatArray, FloatArray]:
        """Return the zero-order-hold discretisation (Ad, Bd, Ed) for the sampling period Ts.

        Results are cached per (model, Ts, V_vent) and shared by every caller, so the arrays
        are read-only; copy them before modifying.

        Parameters
        ----------
        Ts : float
            Sampling period, s.
        V_vent : float, optional
            Mechanical ventilation flow, m3/s.

        Returns
        -------
        Ad : ndarray, shape (2, 2)
            Discrete state matrix, -.
        Bd : ndarray, shape (2,)
            Discrete input vector for Q_c, K/W.
        Ed : ndarray, shape (2, 3)
            Discrete matrix for e = [T_o, Q_air, Q_mass]; first column -, the others K/W.
        """
        return _discretise(self, float(Ts), float(V_vent))

    def gains(self, d: Disturbances) -> FloatArray:
        """Compute the known inputs e = [T_o, Q_air, Q_mass] from (forecast) disturbances.

        Parameters
        ----------
        d : Disturbances
            Disturbances (scalars or arrays over time); uses T_o (degC), G_sol (W/m2) and
            n_occ (-).

        Returns
        -------
        ndarray, shape (n, 3)
            One row per time: outdoor temperature (degC), heat gain to the air node (W) and
            heat gain to the mass node (W).
        """
        p = self.p
        q_sol = p.g_win * p.A_win * np.asarray(d.G_sol, dtype=float)
        q_occ = p.q_person_sens * np.asarray(d.n_occ, dtype=float)
        Q_air = p.f_sol_air * q_sol + p.f_person_conv * q_occ
        Q_mass = (1 - p.f_sol_air) * q_sol + (1 - p.f_person_conv) * q_occ
        return np.column_stack([np.asarray(d.T_o, dtype=float) * np.ones_like(Q_air),
                                Q_air, Q_mass])


@lru_cache(maxsize=64)
def _discretise(model: RCModel, Ts: float,
                V_vent: float) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Discretise the model exactly for inputs held over Ts (see ``RCModel.discrete``)."""
    A, B, E = model.continuous(V_vent)
    n, m = A.shape[0], 1 + E.shape[1]
    M = np.zeros((n + m, n + m))
    M[:n, :n] = A
    M[:n, n] = B
    M[:n, n + 1:] = E
    Phi: FloatArray = expm(M * Ts)
    Ad, Bd, Ed = Phi[:n, :n], Phi[:n, n].copy(), Phi[:n, n + 1:].copy()
    for arr in (Ad, Bd, Ed):
        arr.setflags(write=False)      # cached and shared: callers must not modify them
    return Ad, Bd, Ed
