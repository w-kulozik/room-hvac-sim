"""Measured weather: a typical meteorological year as a disturbance scenario (decision D43).

The year is read from a CSV file written by ``scripts/fetch_weather.py`` (PVGIS 5.3 TMY,
hourly, UTC). Simulation time t is local clock time in seconds from midnight of 1 January,
local = UTC + ``utc_offset_h``, so the daily schedules (occupancy, set points, ventilation)
keep working on clock hours. The year has 365 days and is periodic.

Between the hourly values the outdoor temperature, the outdoor vapour pressure and the
irradiance on the window are interpolated linearly, each at its own time stamps: the
meteorological values at the stamped hour, the irradiance at the stamped hour plus the
irradiance time offset of the file. The humidity ratio is computed from the vapour pressure
at the model's standard pressure, which keeps the outdoor dew point.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import psychro, solar
from ._typing import FloatArray, FloatOrArray
from .inputs import DAY, HOUR, Disturbances, occupants

N_HOURS = 8760
YEAR = N_HOURS * HOUR
_DAYS_IN_MONTH = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
_COLUMNS = ("time_utc", "T2m", "RH", "G_h", "Gb_n", "Gd_h", "SP")


def day_index(month: int, day: int) -> int:
    """Return the number of whole days from 1 January to a date of a 365-day year.

    Parameters
    ----------
    month : int
        Month, 1-12.
    day : int
        Day of the month, from 1.

    Returns
    -------
    int
        Days before the date (0 for 1 January), so the date starts at t = index * DAY.

    Raises
    ------
    ValueError
        If the date does not exist in a 365-day year.
    """
    if not 1 <= month <= 12 or not 1 <= day <= _DAYS_IN_MONTH[month - 1]:
        raise ValueError(f"{month}/{day} (month/day) is not a date of a 365-day year")
    return sum(_DAYS_IN_MONTH[:month - 1]) + day - 1


@dataclass(frozen=True, eq=False)
class TypicalYear:
    """Hourly weather of a typical meteorological year.

    Attributes
    ----------
    latitude, longitude : float
        Site latitude (north positive) and longitude (east positive), degrees.
    elevation : float
        Site elevation, m.
    irradiance_offset_h : float
        Offset of the irradiance values from their time stamps, h.
    T2m : ndarray, shape (8760,)
        Air temperature at 2 m, degC.
    RH : ndarray, shape (8760,)
        Relative humidity, fraction 0..1.
    G_h, Gb_n, Gd_h : ndarray, shape (8760,)
        Global horizontal, beam normal and diffuse horizontal irradiance, W/m2.
    SP : ndarray, shape (8760,)
        Surface air pressure, Pa.
    header : dict of str to str
        The ``key=value`` header lines of the file (sources, months selected).

    Row k is stamped k hours after 1 January 00:00 UTC.
    """

    latitude: float
    longitude: float
    elevation: float
    irradiance_offset_h: float
    T2m: FloatArray
    RH: FloatArray
    G_h: FloatArray
    Gb_n: FloatArray
    Gd_h: FloatArray
    SP: FloatArray
    header: dict[str, str] = field(default_factory=dict)


def _expected_stamps() -> list[str]:
    """Return the time stamps MMDD:HHMM of the 8760 hours of a 365-day year."""
    return [f"{m:02d}{d:02d}:{h:02d}00" for m, n in enumerate(_DAYS_IN_MONTH, start=1)
            for d in range(1, n + 1) for h in range(24)]


def load_typical_year(path: str | Path) -> TypicalYear:
    """Read a typical meteorological year written by ``scripts/fetch_weather.py``.

    Parameters
    ----------
    path : str or Path
        CSV file: ``# key=value`` header lines, a column line, 8760 hourly rows.

    Returns
    -------
    TypicalYear
        The hourly data, with RH converted from % to a fraction.

    Raises
    ------
    ValueError
        If a column or header value is missing, the rows are not the 8760 consecutive hours of a
        365-day year starting on 1 January 00:00, or a value is outside its physical range.
    """
    header: dict[str, str] = {}
    rows: list[list[str]] = []
    columns: list[str] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            key, sep, value = line[1:].strip().partition("=")
            if sep:
                header[key.strip()] = value.strip()
        elif not columns:
            columns = line.split(",")
        elif line:
            rows.append(line.split(","))
    if tuple(columns) != _COLUMNS:
        raise ValueError(f"{path}: columns {columns} differ from the expected {list(_COLUMNS)}")
    for key in ("latitude", "longitude", "elevation_m", "irradiance_time_offset_h"):
        if key not in header:
            raise ValueError(f"{path}: header line '# {key}=...' is missing")
    stamps = [r[0][4:] for r in rows]   # drop the source year of each month
    if stamps != _expected_stamps():
        raise ValueError(f"{path}: the rows must be the 8760 consecutive hours of a 365-day "
                         f"year from 1 January 00:00 UTC ({len(rows)} rows found)")
    values = np.array([[float(v) for v in r[1:]] for r in rows])
    T2m, RH, G_h, Gb_n, Gd_h, SP = (values[:, i] for i in range(6))
    checks = {"T2m": (T2m, -90.0, 60.0), "RH": (RH, 0.0, 100.0), "G_h": (G_h, 0.0, 1500.0),
              "Gb_n": (Gb_n, 0.0, 1500.0), "Gd_h": (Gd_h, 0.0, 1500.0),
              "SP": (SP, 30_000.0, 110_000.0)}
    for name, (v, lo, hi) in checks.items():
        if not (np.all(np.isfinite(v)) and v.min() >= lo and v.max() <= hi):
            raise ValueError(f"{path}: {name} has values outside [{lo}, {hi}] "
                             f"(min {v.min()!r}, max {v.max()!r})")
    return TypicalYear(
        latitude=float(header["latitude"]), longitude=float(header["longitude"]),
        elevation=float(header["elevation_m"]),
        irradiance_offset_h=float(header["irradiance_time_offset_h"]),
        T2m=T2m, RH=RH / 100.0, G_h=G_h, Gb_n=Gb_n, Gd_h=Gd_h, SP=SP, header=header)


def _periodic_linear(t: FloatOrArray, t0: float, values: FloatArray) -> FloatArray:
    """Interpolate hourly values linearly, periodically over the year.

    Parameters
    ----------
    t : float or ndarray
        Time, s.
    t0 : float
        Time of values[0], s; values[k] belongs to t0 + k hours.
    values : ndarray, shape (8760,)
        Hourly values.

    Returns
    -------
    ndarray
        Interpolated values, with the shape of t.
    """
    x = (np.asarray(t, dtype=float) - t0) / HOUR
    k = np.floor(x)
    f = x - k
    i = np.mod(k, N_HOURS).astype(np.intp)
    j = np.where(i + 1 == N_HOURS, 0, i + 1)
    v: FloatArray = values[i] * (1.0 - f) + values[j] * f
    return v


@dataclass(frozen=True, eq=False)
class WeatherScenario:
    """Disturbances from a typical meteorological year, with an office occupancy schedule.

    Attributes
    ----------
    year : TypicalYear
        Hourly weather.
    utc_offset_h : float
        Local clock time minus UTC, h (2 for central European summer time).
    window_tilt : float
        Tilt of the window from the horizontal, degrees (90: vertical).
    window_azimuth : float
        Azimuth of the window's outward normal, degrees from south, positive towards west.
    albedo : float
        Ground reflectance, -.
    n_people : float
        Occupants during the occupied period, -.
    occ_start_h, occ_end_h : float
        Occupied period [occ_start_h, occ_end_h), clock hours (h).
    occupancy : tuple of (float, float, float), optional
        Intervals (start_h, end_h, people); replaces n_people and the occupied period (see
        ``roomsim.inputs.occupants``).
    C_o : float
        Outdoor CO2 concentration, ppm.

    Raises
    ------
    ValueError
        If utc_offset_h is outside [-12, 14], window_tilt outside [0, 180] or albedo outside
        [0, 1].
    """

    year: TypicalYear
    utc_offset_h: float = 2.0
    window_tilt: float = 90.0
    window_azimuth: float = 0.0
    albedo: float = 0.2
    n_people: float = 2.0
    occ_start_h: float = 8.0
    occ_end_h: float = 18.0
    occupancy: tuple[tuple[float, float, float], ...] | None = None
    C_o: float = 420.0
    _p_w: FloatArray = field(init=False, repr=False)
    _G_win: FloatArray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """Validate the settings and precompute the hourly vapour pressure and window sun."""
        if not -12.0 <= self.utc_offset_h <= 14.0:
            raise ValueError(f"utc_offset_h = {self.utc_offset_h!r} is invalid; it must be in "
                             "[-12, 14] h")
        if not 0.0 <= self.window_tilt <= 180.0:
            raise ValueError(f"window_tilt = {self.window_tilt!r} is invalid; it must be in "
                             "[0, 180] degrees")
        y = self.year
        hours = np.arange(N_HOURS, dtype=float) + y.irradiance_offset_h
        cos_t = solar.cos_incidence(hours, y.latitude, y.longitude, self.window_tilt,
                                    self.window_azimuth)
        G_win = solar.plane_irradiance(y.Gb_n, y.Gd_h, y.G_h, cos_t, self.window_tilt,
                                       self.albedo)
        object.__setattr__(self, "_p_w", y.RH * psychro.p_ws(y.T2m))
        object.__setattr__(self, "_G_win", G_win)

    @property
    def _t0_met(self) -> float:
        """Simulation time of the first meteorological value, s."""
        return self.utc_offset_h * HOUR

    @property
    def _t0_sun(self) -> float:
        """Simulation time of the first irradiance value, s."""
        return (self.utc_offset_h + self.year.irradiance_offset_h) * HOUR

    def T_o(self, t: FloatOrArray) -> FloatArray:
        """Compute the outdoor dry-bulb temperature.

        Parameters
        ----------
        t : float or ndarray
            Local time from midnight of 1 January, s.

        Returns
        -------
        ndarray
            Outdoor temperature, degC.
        """
        return _periodic_linear(t, self._t0_met, self.year.T2m)

    def W_o(self, t: FloatOrArray) -> FloatArray:
        """Compute the outdoor humidity ratio at standard pressure.

        Parameters
        ----------
        t : float or ndarray
            Local time from midnight of 1 January, s.

        Returns
        -------
        ndarray
            Outdoor humidity ratio, kg/kg dry air.
        """
        return psychro.W_from_pw(_periodic_linear(t, self._t0_met, self._p_w))

    def G_sol(self, t: FloatOrArray) -> FloatArray:
        """Compute the solar irradiance on the window.

        Parameters
        ----------
        t : float or ndarray
            Local time from midnight of 1 January, s.

        Returns
        -------
        ndarray
            Irradiance on the window plane, W/m2.
        """
        return _periodic_linear(t, self._t0_sun, self._G_win)

    def n_occ(self, t: FloatOrArray) -> FloatArray:
        """Compute the number of occupants.

        Parameters
        ----------
        t : float or ndarray
            Local time from midnight of 1 January, s.

        Returns
        -------
        ndarray
            Number of occupants, -.
        """
        return occupants(t, self.n_people, self.occ_start_h, self.occ_end_h, self.occupancy)

    def __call__(self, t: FloatOrArray) -> Disturbances:
        """Evaluate all disturbances.

        Parameters
        ----------
        t : float or ndarray
            Local time from midnight of 1 January, s.

        Returns
        -------
        Disturbances
            Values at t.
        """
        ones = np.ones_like(np.asarray(t, dtype=float))
        return Disturbances(self.T_o(t), self.W_o(t), self.G_sol(t), self.n_occ(t), self.C_o * ones)


def season_start(month: int, day: int) -> float:
    """Return the simulation time of local midnight at the start of a date.

    Parameters
    ----------
    month, day : int
        Date of a 365-day year.

    Returns
    -------
    float
        Local time from midnight of 1 January, s.
    """
    return day_index(month, day) * DAY
