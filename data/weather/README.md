# Weather data

Typical meteorological years (TMY) for the seasonal study (decision D43 in
[docs/decisions.md](../../docs/decisions.md)).

| File | Site | Latitude | Longitude | Elevation |
|---|---|---|---|---|
| `munich_tmy.csv` | Munich | 48.137 N | 11.575 E | 524 m |
| `madrid_tmy.csv` | Madrid | 40.417 N | 3.704 W | 665 m |

## Source

PVGIS 5.3, TMY tool, European Commission, Joint Research Centre
(<https://joint-research-centre.ec.europa.eu/photovoltaic-geographical-information-system-pvgis_en>).
Downloaded on 24 October 2026 with [scripts/fetch_weather.py](../../scripts/fetch_weather.py),
for the years 2005-2023. Irradiance from PVGIS-SARAH3, the other variables from ERA5. Each month
of a TMY comes from one year, chosen by the ISO 15927-4 procedure; the years are listed in the
file header.

## Licence and changes

© European Union. Reused under the Commission's reuse policy, which licenses the content under
the Creative Commons Attribution 4.0 International licence (CC BY 4.0).

Changes made: the JSON response was converted to CSV; the hourly rows are unchanged, but only
the columns below are kept (the infrared irradiance and the wind are dropped). The location,
the data sources and the months selected are written as `#` header lines.

## Columns

| Column | Meaning | Unit |
|---|---|---|
| `time_utc` | time stamp, `YYYYMMDD:HHMM`, UTC | - |
| `T2m` | air temperature at 2 m | °C |
| `RH` | relative humidity | % |
| `G_h` | global irradiance on the horizontal plane | W/m² |
| `Gb_n` | beam irradiance on a plane normal to the sun | W/m² |
| `Gd_h` | diffuse irradiance on the horizontal plane | W/m² |
| `SP` | surface air pressure | Pa |

The irradiance values belong to the time stamp plus the `irradiance_time_offset_h` given in
the header (the instant of the satellite observation); the other values to the time stamp
itself.
