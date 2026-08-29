"""Solar position and transposition (decision D44)."""

import math
from pathlib import Path

import numpy as np
import pytest

from roomsim import solar
from roomsim.weather import load_typical_year

DATA = Path(__file__).resolve().parents[1] / "data" / "weather"
T_YEAR = np.arange(0.0, 8760.0, 0.25)


def almanac(t_utc_h, year=2019):
    """Low-precision solar coordinates of the astronomical almanac, an independent formulation.

    Returns declination (rad) and equation of time (min) for hours from 1 January 00:00 UTC of
    the given year.
    """
    n = (np.datetime64(f"{year}-01-01T00:00") - np.datetime64("2000-01-01T12:00")) \
        / np.timedelta64(1, "D") + t_utc_h / 24.0
    L = np.radians((280.460 + 0.9856474 * n) % 360.0)
    g = np.radians((357.528 + 0.9856003 * n) % 360.0)
    lam = L + np.radians(1.915) * np.sin(g) + np.radians(0.020) * np.sin(2 * g)
    eps = np.radians(23.439 - 4e-7 * n)
    decl = np.arcsin(np.sin(eps) * np.sin(lam))
    ra = np.arctan2(np.cos(eps) * np.sin(lam), np.cos(lam))
    eot = 4.0 * np.degrees(np.angle(np.exp(1j * (L - ra))))
    return decl, eot


@pytest.mark.parametrize("year", [2019, 2020, 2021, 2022])
def test_solar_angles_agree_with_an_independent_formulation(year):
    # a typical year has no position in the leap-year cycle; over the cycle the date of the
    # equinox moves by about 18 h, i.e. 0.3 degrees of declination, so 0.6 degrees and 1.2 min
    # bound the difference of two correct formulations (measured: 0.3-0.5 deg, 0.7-0.9 min)
    a = solar.solar_angles(T_YEAR, 0.0)
    decl, eot = almanac(T_YEAR, year)
    assert np.degrees(np.abs(a.declination - decl)).max() < 0.6
    assert np.abs(a.equation_of_time - eot).max() < 1.2


def test_solar_angles_extremes():
    a = solar.solar_angles(T_YEAR, 0.0)
    assert np.degrees(a.declination.max()) == pytest.approx(23.44, abs=0.1)
    assert np.degrees(a.declination.min()) == pytest.approx(-23.44, abs=0.1)
    # equation of time: about -14 min in mid-February and +16 min early in November
    assert a.equation_of_time.min() == pytest.approx(-14.2, abs=0.5)
    assert a.equation_of_time.max() == pytest.approx(16.4, abs=0.5)


def test_hour_angle_is_zero_at_solar_noon():
    lon = 11.575
    t = 4000.0                               # any instant
    a = solar.solar_angles(t, lon)
    # solar time = UTC + lon/15 + EoT/60, so solar noon is at UTC 12 - lon/15 - EoT/60, with the
    # equation of time at noon itself (it changes by under a minute a day: two iterations)
    noon = math.floor(t / 24) * 24 + 12.0
    for _ in range(3):
        eot_noon = float(solar.solar_angles(noon, lon).equation_of_time)
        noon = math.floor(t / 24) * 24 + 12.0 - lon / 15.0 - eot_noon / 60.0
    assert float(solar.solar_angles(noon, lon).hour_angle) == pytest.approx(0.0, abs=1e-9)
    assert float(a.hour_angle) == pytest.approx(
        math.radians(15.0 * ((t % 24) + lon / 15.0 + float(a.equation_of_time) / 60.0 - 12.0)))


@pytest.mark.parametrize("tilt,azimuth", [(0, 0), (90, 0), (90, 90), (90, -90), (90, 180),
                                          (35, 20), (60, -135)])
def test_cos_incidence_equals_the_dot_product_of_sun_and_normal(tilt, azimuth):
    # independent construction in local east-north-up coordinates
    lat, lon = 40.417, -3.704
    a = solar.solar_angles(T_YEAR, lon)
    d, w, phi = a.declination, a.hour_angle, math.radians(lat)
    sun = np.stack([-np.cos(d) * np.sin(w),
                    np.cos(phi) * np.sin(d) - np.sin(phi) * np.cos(d) * np.cos(w),
                    np.sin(phi) * np.sin(d) + np.cos(phi) * np.cos(d) * np.cos(w)])
    b, g = math.radians(tilt), math.radians(azimuth)
    normal = np.array([-math.sin(b) * math.sin(g), -math.sin(b) * math.cos(g), math.cos(b)])
    expected = normal @ sun
    got = solar.cos_incidence(T_YEAR, lat, lon, tilt, azimuth)
    assert np.abs(got - expected).max() < 1e-12


def test_vertical_south_window_at_solar_noon():
    # at solar noon (hour angle zero) the general formula reduces to cos theta = sin(phi - delta)
    lat, lon = 48.137, 11.575
    for day in (100, 172, 250):
        noon = day * 24 + 12.0
        for _ in range(3):
            eot = float(solar.solar_angles(noon, lon).equation_of_time)
            noon = day * 24 + 12.0 - lon / 15.0 - eot / 60.0
        d = float(solar.solar_angles(noon, lon).declination)
        c = float(solar.cos_incidence(noon, lat, lon, 90.0, 0.0))
        assert c == pytest.approx(math.sin(math.radians(lat) - d), abs=1e-9)


def test_plane_irradiance_limits():
    # horizontal: beam on the horizontal plus diffuse, no ground reflection
    assert float(solar.plane_irradiance(800.0, 100.0, 500.0, 0.5, 0.0)) == pytest.approx(500.0)
    # vertical: half the sky and half the ground; no beam from behind
    assert float(solar.plane_irradiance(800.0, 100.0, 500.0, -0.3, 90.0, albedo=0.2)) \
        == pytest.approx(50.0 + 50.0)
    with pytest.raises(ValueError, match="albedo"):
        solar.plane_irradiance(0.0, 0.0, 0.0, 0.0, 90.0, albedo=1.5)


@pytest.mark.parametrize("site", ["munich", "madrid"])
def test_geometry_and_time_convention_reproduce_the_pvgis_horizontal_irradiance(site):
    # PVGIS gives G(h) = Gb(n) sin(altitude) + Gd(h) with its own solar geometry. Rebuilding it
    # with ours, at the stamped UTC hour plus the irradiance time offset, must leave only a
    # small residual (measured RMS 1.9 and 1.8 W/m2; with the offset ignored 8.3 and 12.4 W/m2,
    # decision D43). This fixes both the solar position and the time convention of the data.
    y = load_typical_year(DATA / f"{site}_tmy.csv")
    hours = np.arange(8760.0)
    m = y.G_h > 50.0

    def rms(offset):
        cz = solar.cos_incidence(hours + offset, y.latitude, y.longitude, 0.0, 0.0)
        rebuilt = y.Gb_n * np.maximum(cz, 0.0) + y.Gd_h
        return float(np.sqrt(np.mean((rebuilt[m] - y.G_h[m]) ** 2)))

    assert rms(y.irradiance_offset_h) < 3.0
    assert rms(0.0) > 2 * rms(y.irradiance_offset_h)
