"""Cooling season with measured weather, and the seasonal efficiency of the unit (D43-D47).

Runs the on/off thermostat, the PI with sigma-delta low-load handling (C) and the same PI with
the H2 humidity cascade through a cooling season (1 May to 30 September, after 7 warm-up days)
with the typical meteorological years of Munich and Madrid. Computes the SEER of the reference
unit (from its declared data) and of the unit model by the bin method of Regulation (EU)
No 206/2012, with the reference bins and with the bins of each city.

    python scripts/compare_season.py
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from roomsim.ac_unit import ACUnit
from roomsim.controllers import (
    CoilCascade,
    ComfortSchedule,
    OnOffThermostat,
    VentilationSchedule,
    make_pi_controller,
)
from roomsim.inputs import DAY, HOUR
from roomsim.metrics import compute
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.plotting import apply_style, plot_season_monthly
from roomsim.reference_unit import PUBLISHED
from roomsim.seasonal import (
    BIN_HOURS,
    local_bin_hours,
    seer,
    unit_declared_points,
)
from roomsim.simulate import Controller, simulate
from roomsim.weather import WeatherScenario, day_index, load_typical_year

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
DATA = ROOT / "data" / "weather"
SITES = {"Munich": "munich_tmy.csv", "Madrid": "madrid_tmy.csv"}
CONTROLLERS = ("On/off thermostat", "PI (C)", "PI (C) + H2 cascade")
WARMUP_DAYS = 7
MONTHS = ((5, "May"), (6, "Jun"), (7, "Jul"), (8, "Aug"), (9, "Sep"))
SEASON_START, SEASON_END = day_index(5, 1) * DAY, day_index(10, 1) * DAY
COMFORT = ComfortSchedule(early_start_h=6.0)   # 25 degC from 6 to 18 h, 28 degC otherwise (D45)
RH_SET = 0.55

@dataclass(frozen=True)
class SeasonRun:
    """Indicators of one seasonal run (1 May to 30 September)."""

    site: str
    controller: str
    electricity_kWh: float
    cooling_kWh: float
    runtime_h: float
    starts: int
    share_s_min: float           # share of the running time at minimum speed, -
    Kh_above_occ: float          # degree-hours above the occupied band, K h
    Kh_below_occ: float          # degree-hours below the occupied band, K h
    RH_mean_occ: float           # %
    h_RH60_occ: float            # h
    monthly_electricity_kWh: tuple[float, ...]
    monthly_cooling_kWh: tuple[float, ...]

    @property
    def eer_in_use(self) -> float:
        """Return the seasonal EER in use: net cooling over electricity, W/W."""
        return self.cooling_kWh / self.electricity_kWh


def build(name: str, p: RoomParams, unit: ACUnit) -> Controller:
    """Return a fresh controller of the given name, on the set-back schedule (D45)."""
    vent = VentilationSchedule()
    if name == "On/off thermostat":
        return OnOffThermostat(T_set=COMFORT.setpoint, ventilation=vent)
    pi = make_pi_controller("sigma-delta", p, unit, T_set=COMFORT.setpoint, ventilation=vent)
    if name == "PI (C)":
        return pi
    return CoilCascade(pi, RH_set=RH_SET)


def run_case(case: tuple[str, str]) -> SeasonRun:
    """Simulate one site with one controller through the season and reduce the result."""
    site, name = case
    p, unit = RoomParams(), ACUnit()
    S = WeatherScenario(load_typical_year(DATA / SITES[site]))
    t0 = SEASON_START - WARMUP_DAYS * DAY
    x0 = equilibrium_state(22.0, 22.0, float(S.W_o(t0)))
    r = simulate(p, S, x0, SEASON_END - t0, t_start=t0, controller=build(name, p, unit),
                 unit=unit).window(SEASON_START)
    m = compute(r, COMFORT.T_set, COMFORT.band)
    u = r.require_unit()
    s = u.cmd[:-1, 0]
    on = s > 0
    at_min = on & (s <= unit.p.s_min + 1e-9)
    occ = np.asarray(r.d.n_occ)[:-1] > 0
    T = r.T_a[:-1]
    above = np.maximum(T - (COMFORT.T_set + COMFORT.band), 0.0)
    below = np.maximum((COMFORT.T_set - COMFORT.band) - T, 0.0)
    edges = [day_index(month, 1) * DAY for month, _ in MONTHS] + [SEASON_END]
    k = np.searchsorted(r.t, edges)
    return SeasonRun(
        site=site, controller=name, electricity_kWh=m.energy_kWh,
        cooling_kWh=float(u.E_cool[-1] / 3.6e6), runtime_h=m.runtime_h, starts=m.starts,
        share_s_min=float(at_min.sum() / on.sum()) if on.any() else 0.0,
        Kh_above_occ=float(above[occ].sum() * r.dt / HOUR),
        Kh_below_occ=float(below[occ].sum() * r.dt / HOUR),
        RH_mean_occ=m.RH_mean_occ, h_RH60_occ=m.h_RH60_occ,
        monthly_electricity_kWh=tuple(float(v) for v in np.diff(u.E_el[k]) / 3.6e6),
        monthly_cooling_kWh=tuple(float(v) for v in np.diff(u.E_cool[k]) / 3.6e6),
    )


def seer_tables() -> tuple[str, str, dict[str, float]]:
    """Return the SEER table, the table of the model's declared points and SEERon by bins."""
    bins = {"reference season": BIN_HOURS}
    for site, fname in SITES.items():
        bins[site] = local_bin_hours(load_typical_year(DATA / fname).T2m)
    P_model, model_points = unit_declared_points(ACUnit())
    units = {"Reference unit (published data)": (PUBLISHED.declared_points(),
                                                  PUBLISHED.P_design_kW),
             "Unit model": (tuple(model_points), P_model)}
    modes = PUBLISHED.low_power_modes()   # the model has no standby electronics (D46)
    lines = ["| Bins | Bin hours (h) | Reference unit SEERon | Reference unit SEER | "
             "Model SEERon | Model SEER |", "|" + "---|" * 6]
    seer_on_model: dict[str, float] = {}
    for bname, hours in bins.items():
        ref, mod = (seer(points, P_design, modes, bin_hours=hours)
                    for points, P_design in units.values())
        lines.append(f"| {bname} | {hours.sum():.0f} | {ref.SEER_on:.2f} | {ref.SEER:.2f} | "
                     f"{mod.SEER_on:.2f} | {mod.SEER:.2f} |")
        seer_on_model[bname] = mod.SEER_on
    lines.append(f"| declared by the manufacturer | | | {PUBLISHED.SEER:.2f} | | |")
    res_model = seer(model_points, P_model, modes)
    pts = [f"Design load of the unit model: {P_model:.2f} kW", "",
           "| Point | Outdoor (°C) | Part load (kW) | Model capacity (kW) | Model EER | "
           "EERbin |", "|" + "---|" * 6]
    for name, pt, e_bin in zip("ABCD", model_points, res_model.EER_bin, strict=True):
        Pc = P_model * (pt.T_o - 16.0) / 19.0
        pts.append(f"| {name} | {pt.T_o:.0f} | {Pc:.2f} | {pt.Pdc:.2f} | {pt.EERd:.2f} | "
                   f"{e_bin:.2f} |")
    return "\n".join(lines), "\n".join(pts), seer_on_model


def main() -> None:
    """Run the season, write figures/season_comparison.md, figures/seer.md and the figure."""
    FIG.mkdir(exist_ok=True)
    apply_style()
    seer_table, points_table, seer_on_model = seer_tables()
    cases = [(site, name) for site in SITES for name in CONTROLLERS]
    with ProcessPoolExecutor(max_workers=len(cases)) as pool:
        runs = list(pool.map(run_case, cases))

    lines = ["| City, controller | Electricity (kWh) | Net cooling (kWh) | EER in use | "
             "Model SEERon, local bins | Runtime (h) | Starts | Running at minimum speed (%) | "
             "K·h above band occ. | K·h below band occ. | RH mean occ. (%) | "
             "h RH > 60 % occ. |", "|" + "---|" * 12]
    for run in runs:
        lines.append(
            f"| {run.site}, {run.controller} | {run.electricity_kWh:.1f} | "
            f"{run.cooling_kWh:.1f} | {run.eer_in_use:.2f} | {seer_on_model[run.site]:.2f} | "
            f"{run.runtime_h:.0f} | {run.starts} | {100 * run.share_s_min:.0f} | "
            f"{run.Kh_above_occ:.1f} | {run.Kh_below_occ:.1f} | {run.RH_mean_occ:.1f} | "
            f"{run.h_RH60_occ:.0f} |")
    table = "\n".join(lines)
    print("Cooling season 1 May - 30 September, typical meteorological years\n")
    print(table)
    print()
    print(seer_table)
    print()
    print(points_table)
    (FIG / "season_comparison.md").write_text(table + "\n", encoding="utf-8", newline="\n")
    (FIG / "seer.md").write_text(seer_table + "\n\n" + points_table + "\n", encoding="utf-8",
                                 newline="\n")

    monthly = {(run.site, run.controller): (run.monthly_electricity_kWh,
                                            run.monthly_cooling_kWh) for run in runs}
    plot_season_monthly(monthly, list(SITES), list(CONTROLLERS), [m for _, m in MONTHS],
                        FIG / "season_monthly.png")
    print(f"\nFigures written to {FIG}")


if __name__ == "__main__":
    main()
