"""Solar position and irradiance on a tilted plane (decision D44).

Time is given as hours from 1 January, 00:00 UTC, of a 365-day year, as in a typical
meteorological year; it wraps around after 8 760 h.

Solar position: declination and equation of time from Spencer's Fourier series, evaluated at
the fractional year

    gamma = 2 pi / 365 (day - 1 + (hour - 12) / 24),

and the hour angle omega = 15 deg/h (solar time - 12 h), with
solar time = UTC + longitude / 15 deg/h + equation of time.

Angle of incidence on a plane of tilt beta and azimuth gamma_s (measured from south, positive
towards west):

    cos theta = sin d sin phi cos beta - sin d cos phi sin beta cos gamma_s
              + cos d cos phi cos beta cos omega + cos d sin phi sin beta cos gamma_s cos omega
              + cos d sin beta sin gamma_s sin omega

Irradiance on the plane, isotropic sky:

    G = Gb_n max(cos theta, 0) + Gd_h (1 + cos beta) / 2 + rho G_h (1 - cos beta) / 2
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from ._typing import FloatArray, FloatOrArray

HOURS_PER_YEAR = 8760.0


class SolarAngles(NamedTuple):
    """Declination, equation of time and hour angle at given instants.

    Attributes
    ----------
    declination : ndarray
        Solar declination, rad.
    equation_of_time : ndarray
        Equation of time, min (apparent minus mean solar time).
    hour_angle : ndarray
        Hour angle, rad (zero at solar noon, positive in the afternoon).
    """

    declination: FloatArray
    equation_of_time: FloatArray
    hour_angle: FloatArray


def solar_angles(t_utc_h: FloatOrArray, longitude: float) -> SolarAngles:
    """Compute the declination, the equation of time and the hour angle.

    Parameters
    ----------
    t_utc_h : float or ndarray
        Hours from 1 January 00:00 UTC of a 365-day year, h (wrapped modulo 8 760 h).
    longitude : float
        Longitude, degrees (east positive).

    Returns
    -------
    SolarAngles
        Declination (rad), equation of time (min) and hour angle (rad).
    """
    t = np.mod(np.asarray(t_utc_h, dtype=float), HOURS_PER_YEAR)
    hour = np.mod(t, 24.0)
    g = 2.0 * np.pi / 365.0 * (np.floor(t / 24.0) + (hour - 12.0) / 24.0)
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g)
            - 0.006758 * np.cos(2 * g) + 0.000907 * np.sin(2 * g)
            - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    eot = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g)
                    - 0.014615 * np.cos(2 * g) - 0.04089 * np.sin(2 * g))
    solar_time = hour + longitude / 15.0 + eot / 60.0
    omega = np.radians(15.0 * (solar_time - 12.0))
    return SolarAngles(decl, eot, omega)


def cos_incidence(t_utc_h: FloatOrArray, latitude: float, longitude: float, tilt: float,
                  azimuth: float) -> FloatArray:
    """Compute the cosine of the angle of incidence of the beam on a plane.

    Parameters
    ----------
    t_utc_h : float or ndarray
        Hours from 1 January 00:00 UTC of a 365-day year, h.
    latitude, longitude : float
        Site latitude (north positive) and longitude (east positive), degrees.
    tilt : float
        Tilt of the plane from the horizontal, degrees (90 for a wall or window).
    azimuth : float
        Azimuth of the plane's normal, degrees from south, positive towards west.

    Returns
    -------
    ndarray
        cos theta, -; negative when the sun is behind the plane. With ``tilt = 0`` it is the
        cosine of the solar zenith angle, negative when the sun is below the horizon.
    """
    a = solar_angles(t_utc_h, longitude)
    d, w = a.declination, a.hour_angle
    phi, beta, gam = np.radians(latitude), np.radians(tilt), np.radians(azimuth)
    c: FloatArray = (np.sin(d) * np.sin(phi) * np.cos(beta)
                     - np.sin(d) * np.cos(phi) * np.sin(beta) * np.cos(gam)
                     + np.cos(d) * np.cos(phi) * np.cos(beta) * np.cos(w)
                     + np.cos(d) * np.sin(phi) * np.sin(beta) * np.cos(gam) * np.cos(w)
                     + np.cos(d) * np.sin(beta) * np.sin(gam) * np.sin(w))
    return c


def plane_irradiance(Gb_n: FloatOrArray, Gd_h: FloatOrArray, G_h: FloatOrArray,
                     cos_theta: FloatOrArray, tilt: float, albedo: float = 0.2) -> FloatArray:
    """Compute the irradiance on a tilted plane with the isotropic-sky model.

    Parameters
    ----------
    Gb_n : float or ndarray
        Beam irradiance on a plane normal to the sun, W/m2.
    Gd_h : float or ndarray
        Diffuse irradiance on the horizontal plane, W/m2.
    G_h : float or ndarray
        Global irradiance on the horizontal plane, W/m2.
    cos_theta : float or ndarray
        Cosine of the angle of incidence on the plane (``cos_incidence``), -.
    tilt : float
        Tilt of the plane from the horizontal, degrees.
    albedo : float, optional
        Ground reflectance, -.

    Returns
    -------
    ndarray
        Irradiance on the plane: beam + sky diffuse + ground-reflected, W/m2.

    Raises
    ------
    ValueError
        If albedo is outside [0, 1].
    """
    if not 0.0 <= albedo <= 1.0:
        raise ValueError(f"albedo = {albedo!r} is invalid; it must be a fraction in [0, 1]")
    cb = np.cos(np.radians(tilt))
    G: FloatArray = (np.asarray(Gb_n, dtype=float) * np.maximum(cos_theta, 0.0)
                     + np.asarray(Gd_h, dtype=float) * (1.0 + cb) / 2.0
                     + albedo * np.asarray(G_h, dtype=float) * (1.0 - cb) / 2.0)
    return G
