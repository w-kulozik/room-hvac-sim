"""Psychrometric relations for moist air.

All functions accept scalars or numpy arrays. Units:
    T          dry-bulb temperature, degC
    p, p_w     total and water-vapour partial pressure, Pa
    W          humidity ratio, kg water / kg dry air
    RH         relative humidity, fraction 0..1 (not percent)
    h          specific enthalpy, kJ / kg dry air

Saturation pressure uses the Magnus form with the coefficients of Alduchov & Eskridge (1996),
recommended by them for liquid water between -40 and 50 degC; over ice (``p_wsi``), their
coefficients for ice, for -80 to 0 degC.
"""

from __future__ import annotations

import numpy as np

from ._typing import FloatArray, FloatOrArray

P_ATM = 101_325.0  # standard atmospheric pressure, Pa
EPS = 0.622        # ratio of molar masses M_w / M_da (rounded, as is usual in HVAC)

# Magnus coefficients
_A = 610.94   # Pa
_B = 17.625   # -
_C = 243.04   # degC
_A_ICE = 611.21   # Pa
_B_ICE = 22.587   # -
_C_ICE = 273.86   # degC


def p_ws(T: FloatOrArray) -> FloatArray:
    """Compute the saturation vapour pressure over liquid water (Magnus form).

    Parameters
    ----------
    T : float or ndarray
        Dry-bulb temperature, degC.

    Returns
    -------
    ndarray
        Saturation vapour pressure, Pa (a numpy scalar for scalar input).
    """
    T = np.asarray(T, dtype=float)
    out: FloatArray = _A * np.exp(_B * T / (T + _C))
    return out


def p_wsi(T: FloatOrArray) -> FloatArray:
    """Compute the saturation vapour pressure over ice (Magnus form).

    Parameters
    ----------
    T : float or ndarray
        Temperature, degC (intended for T <= 0).

    Returns
    -------
    ndarray
        Saturation vapour pressure over ice, Pa (a numpy scalar for scalar input).
    """
    T = np.asarray(T, dtype=float)
    out: FloatArray = _A_ICE * np.exp(_B_ICE * T / (T + _C_ICE))
    return out


def W_from_pw(p_w: FloatOrArray, p: float = P_ATM) -> FloatArray:
    """Compute the humidity ratio from the vapour partial pressure.

    Parameters
    ----------
    p_w : float or ndarray
        Water-vapour partial pressure, Pa.
    p : float, optional
        Total pressure, Pa.

    Returns
    -------
    ndarray
        Humidity ratio, kg/kg dry air.
    """
    p_w = np.asarray(p_w, dtype=float)
    out: FloatArray = EPS * p_w / (p - p_w)
    return out


def pw_from_W(W: FloatOrArray, p: float = P_ATM) -> FloatArray:
    """Compute the vapour partial pressure from the humidity ratio (inverse of `W_from_pw`).

    Parameters
    ----------
    W : float or ndarray
        Humidity ratio, kg/kg dry air.
    p : float, optional
        Total pressure, Pa.

    Returns
    -------
    ndarray
        Water-vapour partial pressure, Pa.
    """
    W = np.asarray(W, dtype=float)
    out: FloatArray = W * p / (EPS + W)
    return out


def W_from_RH(T: FloatOrArray, RH: FloatOrArray, p: float = P_ATM) -> FloatArray:
    """Compute the humidity ratio at a given temperature and relative humidity.

    Parameters
    ----------
    T : float or ndarray
        Dry-bulb temperature, degC.
    RH : float or ndarray
        Relative humidity, fraction 0..1.
    p : float, optional
        Total pressure, Pa.

    Returns
    -------
    ndarray
        Humidity ratio, kg/kg dry air.
    """
    return W_from_pw(np.asarray(RH, dtype=float) * p_ws(T), p)


def RH_from_W(T: FloatOrArray, W: FloatOrArray, p: float = P_ATM) -> FloatArray:
    """Compute the relative humidity at a given temperature and humidity ratio.

    Parameters
    ----------
    T : float or ndarray
        Dry-bulb temperature, degC.
    W : float or ndarray
        Humidity ratio, kg/kg dry air.
    p : float, optional
        Total pressure, Pa.

    Returns
    -------
    ndarray
        Relative humidity, fraction (not clipped to 0..1).
    """
    out: FloatArray = pw_from_W(W, p) / p_ws(T)
    return out


def dew_point_from_pw(p_w: FloatOrArray) -> FloatArray:
    """Compute the dew-point temperature from the vapour partial pressure (inverse Magnus).

    Parameters
    ----------
    p_w : float or ndarray
        Water-vapour partial pressure, Pa.

    Returns
    -------
    ndarray
        Dew-point temperature, degC; -inf for p_w = 0 (dry air), the limit of the inverse
        Magnus form.
    """
    p_w = np.asarray(p_w, dtype=float)
    dry = p_w == 0.0
    if np.any(dry):
        # log(0) = -inf would give -inf/inf = nan; return the limit instead
        g = np.log(np.where(dry, _A, p_w) / _A)
        limit: FloatArray = np.where(dry, -np.inf, _C * g / (_B - g))[()]
        return limit
    g = np.log(p_w / _A)
    out: FloatArray = _C * g / (_B - g)
    return out


def dew_point(T: FloatOrArray, RH: FloatOrArray) -> FloatArray:
    """Compute the dew-point temperature at a given temperature and relative humidity.

    Parameters
    ----------
    T : float or ndarray
        Dry-bulb temperature, degC.
    RH : float or ndarray
        Relative humidity, fraction 0..1.

    Returns
    -------
    ndarray
        Dew-point temperature, degC.
    """
    return dew_point_from_pw(np.asarray(RH, dtype=float) * p_ws(T))


def dew_point_from_W(W: FloatOrArray, p: float = P_ATM) -> FloatArray:
    """Compute the dew-point temperature from the humidity ratio.

    Parameters
    ----------
    W : float or ndarray
        Humidity ratio, kg/kg dry air.
    p : float, optional
        Total pressure, Pa.

    Returns
    -------
    ndarray
        Dew-point temperature, degC.
    """
    return dew_point_from_pw(pw_from_W(W, p))


def enthalpy(T: FloatOrArray, W: FloatOrArray) -> FloatArray:
    """Compute the specific enthalpy of moist air (reference: dry air and liquid water at 0 degC).

    Parameters
    ----------
    T : float or ndarray
        Dry-bulb temperature, degC.
    W : float or ndarray
        Humidity ratio, kg/kg dry air.

    Returns
    -------
    ndarray
        Specific enthalpy, kJ/kg dry air.
    """
    T = np.asarray(T, dtype=float)
    out: FloatArray = 1.006 * T + np.asarray(W, dtype=float) * (2501.0 + 1.86 * T)
    return out
