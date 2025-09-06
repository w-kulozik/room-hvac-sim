"""Psychrometric relations for moist air.

All functions accept scalars or numpy arrays. Units:
    T          dry-bulb temperature, degC
    p, p_w     total and water-vapour partial pressure, Pa
    W          humidity ratio, kg water / kg dry air
    RH         relative humidity, fraction 0..1 (not percent)
    h          specific enthalpy, kJ / kg dry air

Saturation pressure uses the Magnus form (Alduchov & Eskridge 1996 coefficients),
valid over liquid water for roughly -40..50 degC.
"""

from __future__ import annotations

import numpy as np

P_ATM = 101_325.0  # standard atmospheric pressure, Pa
EPS = 0.622        # ratio of molar masses M_w / M_da (rounded, as is usual in HVAC)

# Magnus coefficients
_A = 610.94   # Pa
_B = 17.625   # -
_C = 243.04   # degC


def p_ws(T):
    """Saturation vapour pressure over water, Pa."""
    T = np.asarray(T, dtype=float)
    return _A * np.exp(_B * T / (T + _C))


def W_from_pw(p_w, p=P_ATM):
    """Humidity ratio from vapour partial pressure."""
    p_w = np.asarray(p_w, dtype=float)
    return EPS * p_w / (p - p_w)


def pw_from_W(W, p=P_ATM):
    """Vapour partial pressure from humidity ratio (inverse of W_from_pw)."""
    W = np.asarray(W, dtype=float)
    return W * p / (EPS + W)


def W_from_RH(T, RH, p=P_ATM):
    """Humidity ratio at temperature T and relative humidity RH (fraction)."""
    return W_from_pw(np.asarray(RH, dtype=float) * p_ws(T), p)


def RH_from_W(T, W, p=P_ATM):
    """Relative humidity (fraction) at temperature T and humidity ratio W."""
    return pw_from_W(W, p) / p_ws(T)


def dew_point_from_pw(p_w):
    """Dew-point temperature, degC, from vapour partial pressure (inverse Magnus)."""
    g = np.log(np.asarray(p_w, dtype=float) / _A)
    return _C * g / (_B - g)


def dew_point(T, RH):
    """Dew-point temperature, degC, at temperature T and relative humidity RH."""
    return dew_point_from_pw(np.asarray(RH, dtype=float) * p_ws(T))


def dew_point_from_W(W, p=P_ATM):
    """Dew-point temperature, degC, from humidity ratio."""
    return dew_point_from_pw(pw_from_W(W, p))


def enthalpy(T, W):
    """Specific enthalpy of moist air, kJ/kg dry air (reference: 0 degC, liquid water)."""
    T = np.asarray(T, dtype=float)
    return 1.006 * T + np.asarray(W, dtype=float) * (2501.0 + 1.86 * T)
