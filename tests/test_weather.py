"""Typical meteorological year: reader and weather scenario (decision D43)."""

from pathlib import Path

import numpy as np
import pytest

from roomsim import psychro, solar
from roomsim.inputs import DAY, HOUR
from roomsim.weather import YEAR, WeatherScenario, day_index, load_typical_year, season_start

DATA = Path(__file__).resolve().parents[1] / "data" / "weather"
MUNICH = load_typical_year(DATA / "munich_tmy.csv")


def test_reader_returns_the_rows_of_the_file():
    y = MUNICH
    assert (y.latitude, y.longitude, y.elevation) == (48.137, 11.575, 524.0)
    assert y.irradiance_offset_h == pytest.approx(0.1795)
    assert y.T2m.shape == (8760,)
    # first row: 20070101:0000,7.47,65.61,0.0,0.0,0.0,95440.0
    assert (y.T2m[0], y.RH[0], y.SP[0]) == (7.47, pytest.approx(0.6561), 95440.0)
    assert y.header["radiation_db"] == "PVGIS-SARAH3"
    assert np.all((y.RH >= 0) & (y.RH <= 1))


def _rewrite(tmp_path, transform):
    text = (DATA / "munich_tmy.csv").read_text(encoding="utf-8")
    path = tmp_path / "bad.csv"
    path.write_text(transform(text), encoding="utf-8")
    return path


def test_reader_rejects_malformed_files(tmp_path):
    lines = (DATA / "munich_tmy.csv").read_text(encoding="utf-8").splitlines()
    header_end = next(i for i, line in enumerate(lines) if not line.startswith("#"))
    missing_row = _rewrite(tmp_path, lambda t: "\n".join(lines[:header_end + 5]
                                                         + lines[header_end + 6:]))
    with pytest.raises(ValueError, match="8760 consecutive hours"):
        load_typical_year(missing_row)
    with pytest.raises(ValueError, match="columns"):
        load_typical_year(_rewrite(tmp_path, lambda t: t.replace("Gd_h", "Gd")))
    with pytest.raises(ValueError, match="irradiance_time_offset_h"):
        load_typical_year(_rewrite(tmp_path, lambda t: t.replace("# irradiance_time", "# x")))
    with pytest.raises(ValueError, match="RH has values outside"):
        load_typical_year(_rewrite(tmp_path, lambda t: t.replace(",7.47,65.61,", ",7.47,165.61,")))


def test_day_index():
    assert day_index(1, 1) == 0
    assert day_index(5, 1) == 120
    assert day_index(12, 31) == 364
    assert season_start(5, 1) == 120 * DAY
    for bad in ((2, 29), (13, 1), (4, 31), (1, 0)):
        with pytest.raises(ValueError, match="365-day year"):
            day_index(*bad)


def test_values_at_the_time_stamps_and_between_them():
    S = WeatherScenario(MUNICH, utc_offset_h=2.0)
    k = 4500                                     # UTC hour of the year
    t = (k + 2.0) * HOUR                         # local time of that stamp
    assert float(S.T_o(t)) == MUNICH.T2m[k]
    assert float(S.T_o(t + 0.5 * HOUR)) == pytest.approx(0.5 * (MUNICH.T2m[k] + MUNICH.T2m[k + 1]))
    # the humidity ratio keeps the vapour pressure, hence the dew point, of the data
    p_w = MUNICH.RH[k] * float(psychro.p_ws(MUNICH.T2m[k]))
    assert float(psychro.pw_from_W(S.W_o(t))) == pytest.approx(p_w, rel=1e-12)
    # irradiance on the window at its own stamp (UTC hour + offset)
    t_sun = (k + MUNICH.irradiance_offset_h + 2.0) * HOUR
    cos_t = solar.cos_incidence(k + MUNICH.irradiance_offset_h, MUNICH.latitude,
                                MUNICH.longitude, 90.0, 0.0)
    expected = solar.plane_irradiance(MUNICH.Gb_n[k], MUNICH.Gd_h[k], MUNICH.G_h[k], cos_t, 90.0)
    assert float(S.G_sol(t_sun)) == pytest.approx(float(expected), rel=1e-12)


def test_year_is_periodic_and_continuous_across_new_year():
    S = WeatherScenario(MUNICH, utc_offset_h=1.0)
    last = (8759 + 1.0) * HOUR                   # local time of the last stamp
    assert float(S.T_o(last)) == MUNICH.T2m[-1]
    assert float(S.T_o(last + HOUR)) == MUNICH.T2m[0]
    assert float(S.T_o(last + 0.5 * HOUR)) == pytest.approx(0.5 * (MUNICH.T2m[-1] + MUNICH.T2m[0]))
    t = np.array([3.0e6, 1.7e7])
    assert np.allclose(S.T_o(t + YEAR), S.T_o(t), rtol=1e-12, atol=0.0)
    assert np.allclose(S.T_o(t - YEAR), S.T_o(t), rtol=1e-12, atol=0.0)


def test_vectorised_matches_scalar():
    S = WeatherScenario(MUNICH)
    t = season_start(7, 1) + np.linspace(0.0, 3 * DAY, 97)
    d = S(t)
    for i in (0, 13, 50, 96):
        di = S(float(t[i]))
        for name in ("T_o", "W_o", "G_sol", "n_occ", "C_o"):
            assert float(getattr(di, name)) == getattr(d, name)[i]


def test_sun_on_a_south_window_peaks_near_local_solar_noon():
    # solar noon in Munich (11.575 E) is at about 12:14 UTC + EoT; with UTC+2 near 14:14 minus
    # the equation of time. The hourly interpolation can move the maximum by up to an hour.
    S = WeatherScenario(MUNICH, utc_offset_h=2.0)
    for month in (5, 6, 7, 8, 9):
        t0 = season_start(month, 10)
        t = t0 + np.arange(0.0, DAY, 60.0)
        peak_h = float(t[np.argmax(S.G_sol(t))] - t0) / HOUR
        eot = float(solar.solar_angles((t0 / HOUR) + 12.0, MUNICH.longitude).equation_of_time)
        noon_local = 12.0 - MUNICH.longitude / 15.0 - eot / 60.0 + 2.0
        assert abs(peak_h - noon_local) < 1.0 + MUNICH.irradiance_offset_h
    night = season_start(7, 1) + 2.0 * HOUR
    assert float(S.G_sol(night)) == 0.0


def test_occupancy_follows_the_local_clock():
    S = WeatherScenario(MUNICH, utc_offset_h=2.0)
    t0 = season_start(6, 1)
    assert float(S.n_occ(t0 + 7.99 * HOUR)) == 0.0
    assert float(S.n_occ(t0 + 8.0 * HOUR)) == 2.0
    assert float(S.n_occ(t0 + 18.0 * HOUR)) == 0.0


@pytest.mark.parametrize("kwargs,match", [({"utc_offset_h": 15.0}, "utc_offset_h"),
                                          ({"window_tilt": -5.0}, "window_tilt"),
                                          ({"albedo": 2.0}, "albedo")])
def test_scenario_validation(kwargs, match):
    with pytest.raises(ValueError, match=match):
        WeatherScenario(MUNICH, **kwargs)


def test_cities_differ_as_their_climates():
    madrid = WeatherScenario(load_typical_year(DATA / "madrid_tmy.csv"))
    munich = WeatherScenario(MUNICH)
    t = season_start(5, 1) + np.arange(0.0, 153 * DAY, HOUR)
    # May-September means 23.4 and 15.9 degC (decision D43)
    assert float(madrid.T_o(t).mean()) == pytest.approx(23.4, abs=0.1)
    assert float(munich.T_o(t).mean()) == pytest.approx(15.9, abs=0.1)
