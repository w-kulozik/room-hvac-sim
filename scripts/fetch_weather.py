"""Download the typical meteorological years used by the seasonal study and store them as CSV.

The data come from the PVGIS 5.3 web service of the European Commission's Joint Research Centre
(TMY tool, decision D43). The CSV files keep the hourly rows unchanged but only the columns the
simulator uses, with the location, the data sources and the months selected in a header. The
committed files in data/weather/ were written by this script; it is not part of reproduce.py,
because it needs the network and the service may be updated.

    python scripts/fetch_weather.py
"""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path
from typing import Any

DATA = Path(__file__).resolve().parents[1] / "data" / "weather"
API = "https://re.jrc.ec.europa.eu/api/v5_3/tmy"
SITES = {"munich": (48.137, 11.575), "madrid": (40.417, -3.704)}
START_YEAR, END_YEAR = 2005, 2023
COLUMNS = (("time(UTC)", "time_utc"), ("T2m", "T2m"), ("RH", "RH"), ("G(h)", "G_h"),
           ("Gb(n)", "Gb_n"), ("Gd(h)", "Gd_h"), ("SP", "SP"))


def fetch(lat: float, lon: float) -> dict[str, Any]:
    """Return the PVGIS TMY response for one location.

    Parameters
    ----------
    lat, lon : float
        Latitude and longitude, decimal degrees (east positive).

    Returns
    -------
    dict
        The decoded JSON response.
    """
    url = (f"{API}?lat={lat}&lon={lon}&startyear={START_YEAR}&endyear={END_YEAR}"
           "&outputformat=json")
    with urllib.request.urlopen(url, timeout=120) as response:
        data: dict[str, Any] = json.load(response)
    return data


def to_csv(data: dict[str, Any]) -> str:
    """Convert a PVGIS TMY response to the CSV format read by ``roomsim.weather``.

    Parameters
    ----------
    data : dict
        Decoded PVGIS TMY response.

    Returns
    -------
    str
        Header lines starting with '#', a column line and 8760 hourly rows.
    """
    loc, meteo = data["inputs"]["location"], data["inputs"]["meteo_data"]
    months = ",".join(f"{m['month']}:{m['year']}" for m in data["outputs"]["months_selected"])
    lines = [
        "# PVGIS 5.3 typical meteorological year (European Commission, Joint Research Centre)",
        "# (c) European Union, reused under CC BY 4.0; only the columns below are kept",
        f"# latitude={loc['latitude']}",
        f"# longitude={loc['longitude']}",
        f"# elevation_m={loc['elevation']}",
        f"# irradiance_time_offset_h={loc['irradiance_time_offset']}",
        f"# radiation_db={meteo['radiation_db']}",
        f"# meteo_db={meteo['meteo_db']}",
        f"# years={meteo['year_min']}-{meteo['year_max']}",
        f"# months_selected={months}",
        ",".join(name for _, name in COLUMNS),
    ]
    for row in data["outputs"]["tmy_hourly"]:
        lines.append(",".join(str(row[key]) for key, _ in COLUMNS))
    return "\n".join(lines) + "\n"


def main() -> None:
    """Download every site and write data/weather/<site>_tmy.csv."""
    DATA.mkdir(parents=True, exist_ok=True)
    for name, (lat, lon) in SITES.items():
        path = DATA / f"{name}_tmy.csv"
        path.write_text(to_csv(fetch(lat, lon)), encoding="utf-8", newline="\n")
        print(f"{path.name} written")


if __name__ == "__main__":
    main()
