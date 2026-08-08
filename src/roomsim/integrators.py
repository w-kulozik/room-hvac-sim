"""Fixed-step explicit integrators for dx/dt = f(t, x).

Fixed steps because the discrete controllers update on a regular grid and the
verification tests measure the order of convergence as the step is halved.
"""

from __future__ import annotations

from collections.abc import Callable

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


def _rk4_amplification(z: np.ndarray) -> np.ndarray:
    """Return the RK4 amplification factor R(z) = 1 + z + z^2/2 + z^3/6 + z^4/24."""
    return 1.0 + z * (1.0 + z * (0.5 + z * (1.0 / 6.0 + z / 24.0)))


def rk4_stability_limit(J: np.ndarray) -> float:
    """Largest stable classical RK4 step for the linear(ised) system dx/dt = J x, s.

    RK4 is stable when |R(dt*lambda)| <= 1 for every eigenvalue lambda of J, with
    R(z) = 1 + z + z^2/2 + z^3/6 + z^4/24. For a real negative eigenvalue this is
    dt <= 2.785/|lambda|. The limit is found by scanning dt along each eigenvalue's ray on a
    grid of 1e-4 * 4/|lambda|, so it is accurate to about 1e-4 relative.

    Parameters
    ----------
    J : ndarray, shape (n, n)
        Jacobian of the right-hand side, 1/s.

    Returns
    -------
    float
        Stability limit, s.

    Raises
    ------
    ValueError
        If J has an eigenvalue with non-negative real part other than zero.
    """
    lam = np.linalg.eigvals(J)
    lam = lam[np.abs(lam) > 0]
    if np.any(lam.real >= 0):
        raise ValueError("system has non-decaying modes; no stability limit for decay")
    s = np.linspace(1e-4, 4.0, 40_000)            # dt*|lambda|; |R| > 1 beyond 4 on every ray
    limits = []
    for mu in lam:
        unstable = np.abs(_rk4_amplification(s * mu / abs(mu))) > 1.0
        k = int(np.argmax(unstable)) if np.any(unstable) else len(s)
        limits.append(s[k - 1] / abs(mu) if k > 0 else 0.0)
    return float(min(limits))
