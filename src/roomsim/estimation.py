"""Kalman filter for the controller: estimates the unmeasured mass temperature and an
unknown heat input to the air, from the air-temperature sensor alone.

Augmented state  z = [T_a, T_m, d],  d in W, random walk.
    z_k+1 = [[Ad, Bd_d], [0, 1]] z_k + [Bd; 0] Q_c,k + [Ed; 0] e_k + w_k
    y_k   = T_a,k + v_k
where Bd_d = -Bd (d enters like a negative cooling input, i.e. as heat into the air).
The disturbance state gives offset-free predictions: any persistent error of the model
(wrong capacities or conductances, unmodelled latent/sensible split, fan heat, sensor bias)
shows up as a d that the MPC then compensates, as the integral term of a PI would.

The gain is the steady-state one from the discrete algebraic Riccati equation, computed for
each ventilation level (the model is piecewise time-invariant).
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_discrete_are

from .control_model import RCModel
from .inputs import Disturbances


def augmented(model: RCModel, Ts: float, V_vent: float):
    Ad, Bd, Ed = model.discrete(Ts, V_vent)
    A = np.zeros((3, 3))
    A[:2, :2] = Ad
    A[:2, 2] = -Bd
    A[2, 2] = 1.0
    B = np.append(Bd, 0.0)
    E = np.vstack([Ed, np.zeros(Ed.shape[1])])
    return A, B, E


C_MEAS = np.array([[1.0, 0.0, 0.0]])


class DisturbanceKalman:
    def __init__(self, model: RCModel, Ts: float = 60.0, sigma_Ta: float = 0.01,
                 sigma_Tm: float = 0.002, sigma_d: float = 20.0, sigma_y: float = 0.05):
        self.model, self.Ts = model, Ts
        self.Q = np.diag([sigma_Ta ** 2, sigma_Tm ** 2, sigma_d ** 2])
        self.R = np.array([[sigma_y ** 2]])
        self._gain = {}
        self.z = None

    def gain(self, V_vent: float) -> np.ndarray:
        if V_vent not in self._gain:
            A, _, _ = augmented(self.model, self.Ts, V_vent)
            P = solve_discrete_are(A.T, C_MEAS.T, self.Q, self.R)
            S = C_MEAS @ P @ C_MEAS.T + self.R
            self._gain[V_vent] = (P @ C_MEAS.T / S).ravel()   # gain on the predicted state
        return self._gain[V_vent]

    def initialise(self, T_a: float) -> None:
        self.z = np.array([T_a, T_a, 0.0])

    def update(self, y: float, Q_c_prev: float, d_prev: Disturbances, V_vent_prev: float) -> np.ndarray:
        """Predict over the last sample with the input that was applied, then correct with y."""
        if self.z is None:
            self.initialise(y)
            return self.z
        A, B, E = augmented(self.model, self.Ts, V_vent_prev)
        e = self.model.gains(d_prev)[0]
        z_pred = A @ self.z + B * Q_c_prev + E @ e
        K = self.gain(V_vent_prev)
        self.z = z_pred + K * (y - z_pred[0])
        return self.z
