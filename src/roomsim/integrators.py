"""Fixed-step explicit integrators for dx/dt = f(t, x).

Fixed steps because the discrete controllers update on a regular grid and the
verification tests measure the order of convergence as the step is halved.
"""

from __future__ import annotations

from typing import Callable

import numpy as np

RHS = Callable[[float, np.ndarray], np.ndarray]


def euler_step(f: RHS, t: float, x: np.ndarray, dt: float) -> np.ndarray:
    """Explicit (forward) Euler, order 1."""
    return x + dt * f(t, x)


def rk4_step(f: RHS, t: float, x: np.ndarray, dt: float) -> np.ndarray:
    """Classical fourth-order Runge-Kutta."""
    k1 = f(t, x)
    k2 = f(t + 0.5 * dt, x + 0.5 * dt * k1)
    k3 = f(t + 0.5 * dt, x + 0.5 * dt * k2)
    k4 = f(t + dt, x + dt * k3)
    return x + dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


STEPPERS = {"euler": euler_step, "rk4": rk4_step}


def euler_stability_limit(J: np.ndarray) -> float:
    """Largest stable explicit-Euler step for the linear(ised) system dx/dt = J x, s.

    Euler is stable when |1 + dt*lambda| <= 1 for every eigenvalue lambda of J.
    For a real negative eigenvalue this is dt <= 2/|lambda|; the general form below
    also covers complex pairs: dt <= -2 Re(lambda) / |lambda|^2.
    """
    lam = np.linalg.eigvals(J)
    lam = lam[np.abs(lam) > 0]
    if np.any(lam.real >= 0):
        raise ValueError("system has non-decaying modes; no stability limit for decay")
    return float(np.min(-2.0 * lam.real / np.abs(lam) ** 2))
