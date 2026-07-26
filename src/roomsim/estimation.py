"""Kalman filter for the controller.

Estimates the unmeasured mass temperature and an unknown heat input to the air, from the
air-temperature sensor alone.

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

from ._typing import FloatArray
from .control_model import RCModel
from .inputs import Disturbances


def augmented(model: RCModel, Ts: float,
              V_vent: float) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Build the discrete model augmented with the disturbance state d.

    Parameters
    ----------
    model : RCModel
        Controller model.
    Ts : float
        Sampling period, s.
    V_vent : float
        Mechanical ventilation flow over the step, m3/s.

    Returns
    -------
    A : ndarray, shape (3, 3)
        Augmented state matrix for z = [T_a, T_m, d] (d in W).
    B : ndarray, shape (3,)
        Input vector for the sensible cooling Q_c, K/W (zero for d).
    E : ndarray, shape (3, 3)
        Matrix for the known inputs e = [T_o, Q_air, Q_mass] (zero row for d).
    """
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
    """Steady-state Kalman filter for z = [T_a, T_m, d] with the air temperature as output.

    Parameters
    ----------
    model : RCModel
        Controller model.
    Ts : float, optional
        Sampling period, s.
    sigma_Ta, sigma_Tm : float, optional
        Process-noise standard deviations of the air and mass temperatures per sample, K.
    sigma_d : float, optional
        Process-noise standard deviation of the disturbance random walk per sample, W.
    sigma_y : float, optional
        Measurement-noise standard deviation of the air-temperature sensor, K.

    Attributes
    ----------
    Q : ndarray, shape (3, 3)
        Process-noise covariance (diagonal), K2 and W2.
    R : ndarray, shape (1, 1)
        Measurement-noise covariance, K2.
    z : ndarray, shape (3,) or None
        Current state estimate [T_a (degC), T_m (degC), d (W)]; None before initialisation.
    """

    def __init__(self, model: RCModel, Ts: float = 60.0, sigma_Ta: float = 0.01,
                 sigma_Tm: float = 0.002, sigma_d: float = 20.0, sigma_y: float = 0.05) -> None:
        self.model, self.Ts = model, Ts
        self.Q = np.diag([sigma_Ta ** 2, sigma_Tm ** 2, sigma_d ** 2])
        self.R = np.array([[sigma_y ** 2]])
        self._gain: dict[float, FloatArray] = {}
        self.z: FloatArray | None = None

    @property
    def state(self) -> FloatArray:
        """Current state estimate [T_a (degC), T_m (degC), d (W)].

        Raises
        ------
        RuntimeError
            If the filter has not been initialised.
        """
        if self.z is None:
            raise RuntimeError("Kalman filter not initialised")
        return self.z

    def gain(self, V_vent: float) -> FloatArray:
        """Return the steady-state gain on the predicted state for a ventilation level.

        Computed from the discrete algebraic Riccati equation on first use and cached.

        Parameters
        ----------
        V_vent : float
            Mechanical ventilation flow, m3/s.

        Returns
        -------
        ndarray, shape (3,)
            Kalman gain; elements in K/K, K/K and W/K.
        """
        if V_vent not in self._gain:
            A, _, _ = augmented(self.model, self.Ts, V_vent)
            P = solve_discrete_are(A.T, C_MEAS.T, self.Q, self.R)
            S = C_MEAS @ P @ C_MEAS.T + self.R
            self._gain[V_vent] = (P @ C_MEAS.T / S).ravel()   # gain on the predicted state
        return self._gain[V_vent]

    def initialise(self, T_a: float) -> None:
        """Set the estimate to T_a for both temperatures and zero disturbance.

        Parameters
        ----------
        T_a : float
            Measured air temperature, degC.
        """
        self.z = np.array([T_a, T_a, 0.0])

    def update(self, y: float, Q_c_prev: float, d_prev: Disturbances,
               V_vent_prev: float) -> FloatArray:
        """Predict over the last sample with the input that was applied, then correct with y.

        On the first call (no estimate yet) the filter is initialised from y instead.

        Parameters
        ----------
        y : float
            Measured air temperature, degC.
        Q_c_prev : float
            Sensible cooling applied over the last sample, W.
        d_prev : Disturbances
            Disturbances at the start of the last sample.
        V_vent_prev : float
            Mechanical ventilation flow over the last sample, m3/s.

        Returns
        -------
        ndarray, shape (3,)
            Updated estimate [T_a (degC), T_m (degC), d (W)] (the array stored in ``z``).
        """
        if self.z is None:
            self.initialise(y)
            return self.state
        A, B, E = augmented(self.model, self.Ts, V_vent_prev)
        e = self.model.gains(d_prev)[0]
        z_pred = A @ self.z + B * Q_c_prev + E @ e
        K = self.gain(V_vent_prev)
        z: FloatArray = z_pred + K * (y - z_pred[0])
        self.z = z
        return z
