"""Heating season with a heat pump: controllers, defrost strategies and SCOP (D50-D57).

Runs the on/off thermostat with timed defrost, the PI (C) with timed defrost and the PI (C) with
demand defrost through a heating season (1 October to 30 April, after 7 warm-up days) with the
typical meteorological years of Munich and Madrid, plus a sensitivity of the demand case in
Munich to the frost blockage parameter. Computes the SCOP of the reference unit (from its
declared data) and of the heat-pump model by the bin method of Regulation (EU) No 206/2012.

    python scripts/compare_heating.py
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from roomsim.controllers import (
    ComfortSchedule,
    DefrostControl,
    OnOffThermostat,
    VentilationSchedule,
    make_heating_pi,
)
from roomsim.heat_pump import HeatPump, HeatPumpParams
from roomsim.inputs import DAY, HOUR
from roomsim.metrics import compute
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.plotting import apply_style, plot_season_monthly
from roomsim.reference_unit import PUBLISHED
from roomsim.seasonal import (
    HEATING_BIN_HOURS,
    HEATING_BIN_T,
    heat_pump_declared_points,
    scop,
)
from roomsim.simulate import Controller, simulate
from roomsim.weather import YEAR, WeatherScenario, day_index, load_typical_year

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
DATA = ROOT / "data" / "weather"
SITES = {"Munich": "munich_tmy.csv", "Madrid": "madrid_tmy.csv"}
CASES = (("On/off, timed defrost", "on/off", "timed"), ("PI (C), timed defrost", "pi", "timed"),
         ("PI (C), demand defrost", "pi", "demand"))
SENSITIVITY = (("PI (C), demand defrost, M_half 0.25 kg", 0.25),
               ("PI (C), demand defrost, M_half 1.0 kg", 1.0))
WARMUP_DAYS = 7
UTC_OFFSET_H = 1.0                                        # central European time (D54)
SEASON_START = day_index(10, 1) * DAY
SEASON_END = YEAR + day_index(5, 1) * DAY                 # 1 May of the following year
MONTHS = ((10, "Oct"), (11, "Nov"), (12, "Dec"), (1, "Jan"), (2, "Feb"), (3, "Mar"), (4, "Apr"))
COMFORT = ComfortSchedule(T_set=20.0, band=1.0, T_hi_unocc=16.0, early_start_h=6.0)

@dataclass(frozen=True)
class HeatingRun:
    """Indicators of one heating season (1 October to 30 April)."""

    site: str
    case: str
    electricity_kWh: float
    heat_kWh: float               # heat delivered to the room, net of the defrost heat taken
    runtime_h: float
    starts: int
    defrosts: int
    defrost_h: float
    defrost_electricity_kWh: float
    defrost_room_kWh: float       # heat taken from the room during defrost
    share_s_max: float            # share of the running time at maximum speed, -
    Kh_below_occ: float           # degree-hours below the occupied band, K h
    frost_max_kg: float
    monthly_electricity_kWh: tuple[float, ...]
    monthly_heat_kWh: tuple[float, ...]

    @property
    def cop_in_use(self) -> float:
        """Return the seasonal COP in use: heat delivered over electricity, W/W."""
        return self.heat_kWh / self.electricity_kWh

    @property
    def defrost_penalty(self) -> float:
        """Return (defrost electricity + room heat taken in defrost) / electricity, -."""
        return (self.defrost_electricity_kWh + self.defrost_room_kWh) / self.electricity_kWh


def build(kind: str, strategy: str, p: RoomParams, hp: HeatPump) -> Controller:
    """Return a fresh heating controller wrapped in the defrost supervisor (D52-D53)."""
    vent = VentilationSchedule()
    room: Controller
    if kind == "on/off":
        room = OnOffThermostat(T_set=COMFORT.setpoint, ventilation=vent, heating=True,
                               s_on=hp.hp.s_max)
    else:
        room = make_heating_pi(p, hp, T_set=COMFORT.setpoint, ventilation=vent)
    return DefrostControl(room, hp, strategy, ventilation=vent)


def run_case(job: tuple[str, str, str, str, float]) -> HeatingRun:
    """Simulate one site and case through the heating season and reduce the result."""
    site, case, kind, strategy, M_half = job
    p = RoomParams()
    hp = HeatPump(hp=replace(HeatPumpParams(), M_half=M_half))
    S = WeatherScenario(load_typical_year(DATA / SITES[site]), utc_offset_h=UTC_OFFSET_H)
    t0 = SEASON_START - WARMUP_DAYS * DAY
    x0 = equilibrium_state(18.0, 18.0, float(S.W_o(t0)))
    ctrl = build(kind, strategy, p, hp)
    r = simulate(p, S, x0, SEASON_END - t0, t_start=t0, controller=ctrl, unit=hp)
    r = r.window(SEASON_START)
    m = compute(r, COMFORT.T_set, COMFORT.band)
    u = r.require_unit()
    s = u.cmd[:-1, 0]
    on, defrost = s > 0, s < 0
    starts_defrost = int(np.sum(defrost[1:] & ~defrost[:-1]) + defrost[0])
    occ = np.asarray(r.d.n_occ)[:-1] > 0
    below = np.maximum((COMFORT.T_set - COMFORT.band) - r.T_a[:-1], 0.0)
    dE = np.diff(u.E_el)
    edges = [day_index(mo, 1) * DAY + (0.0 if mo >= 10 else YEAR) for mo, _ in MONTHS]
    edges.append(SEASON_END)
    k = np.searchsorted(r.t, edges)
    return HeatingRun(
        site=site, case=case, electricity_kWh=m.energy_kWh,
        heat_kWh=float(-u.E_cool[-1] / 3.6e6), runtime_h=m.runtime_h, starts=m.starts,
        defrosts=starts_defrost, defrost_h=float(defrost.sum() * r.dt / HOUR),
        defrost_electricity_kWh=float(dE[defrost].sum() / 3.6e6),
        defrost_room_kWh=float(u.op["Q_sens"][:-1][defrost].sum() * r.dt / 3.6e6),
        share_s_max=float(np.sum(on & (s >= hp.hp.s_max - 1e-9)) / on.sum()) if on.any() else 0.0,
        Kh_below_occ=float(below[occ].sum() * r.dt / HOUR),
        frost_max_kg=float(u.M_coil.max()),
        monthly_electricity_kWh=tuple(float(v) for v in np.diff(u.E_el[k]) / 3.6e6),
        monthly_heat_kWh=tuple(float(v) for v in -np.diff(u.E_cool[k]) / 3.6e6),
    )


def local_heating_bins(T_hourly: np.ndarray) -> np.ndarray:
    """Count the hours of a weather record in the heating bins -10 to 15 degC (nearest degree).

    Hours below -10.5 degC, which the average season does not have, are not counted.
    """
    T_bin = np.floor(np.asarray(T_hourly, dtype=float) + 0.5)
    return np.array([np.sum(T_bin == Tj) for Tj in HEATING_BIN_T], dtype=float)


def scop_tables() -> tuple[str, str, dict[str, float]]:
    """Return the SCOP table, the table of the model's declared points and SCOPon by bins."""
    bins = {"average season": HEATING_BIN_HOURS}
    for site, fname in SITES.items():
        bins[site] = local_heating_bins(load_typical_year(DATA / fname).T2m)
    P_model, model_points = heat_pump_declared_points(HeatPump())
    pub_points, modes = PUBLISHED.heating_declared_points(), PUBLISHED.heating_low_power_modes()
    lines = ["| Bins | Bin hours (h) | Reference unit SCOPon | Reference unit SCOP | "
             "Model SCOPon | Model SCOP |", "|" + "---|" * 6]
    scop_on_model: dict[str, float] = {}
    for bname, hours in bins.items():
        ref = scop(pub_points, PUBLISHED.P_design_heating_kW, PUBLISHED.P_dh_design_kW, modes,
                   bin_hours=hours)
        mod = scop(model_points, P_model, P_model, modes, bin_hours=hours)
        lines.append(f"| {bname} | {hours.sum():.0f} | {ref.SEER_on:.2f} | {ref.SEER:.2f} | "
                     f"{mod.SEER_on:.2f} | {mod.SEER:.2f} |")
        scop_on_model[bname] = mod.SEER_on
    lines.append(f"| declared by the manufacturer | | | {PUBLISHED.SCOP:.2f} | | |")
    res_model = scop(model_points, P_model, P_model, modes)
    pts = [f"Design load of the heat-pump model: {P_model:.2f} kW (its output at s_max and "
           "-10 °C)", "",
           "| Point | Outdoor (°C) | Published capacity (kW) | Published COP | "
           "Model capacity (kW) | Model COP | Model COPbin |", "|" + "---|" * 7]
    for name, pub, pt, c_bin in zip("ABCD", pub_points, model_points, res_model.EER_bin,
                                    strict=True):
        pts.append(f"| {name} | {pt.T_o:.0f} | {pub.Pdc:.2f} | {pub.EERd:.2f} | {pt.Pdc:.2f} | "
                   f"{pt.EERd:.2f} | {c_bin:.2f} |")
    return "\n".join(lines), "\n".join(pts), scop_on_model


def main() -> None:
    """Run the heating season; write figures/heating_comparison.md, scop.md and the figure."""
    FIG.mkdir(exist_ok=True)
    apply_style()
    scop_table, points_table, scop_on_model = scop_tables()
    M_half = HeatPumpParams().M_half
    jobs = [(site, case, kind, strategy, M_half) for site in SITES
            for case, kind, strategy in CASES]
    jobs += [("Munich", case, "pi", "demand", mh) for case, mh in SENSITIVITY]
    with ProcessPoolExecutor(max_workers=len(jobs)) as pool:
        runs = list(pool.map(run_case, jobs))

    lines = ["| City, case | Electricity (kWh) | Heat delivered (kWh) | COP in use | "
             "Model SCOPon, local bins | Defrosts | Defrost time (h) | Defrost penalty (%) | "
             "Frost max (kg) | Running at maximum speed (%) | Starts | "
             "K·h below band occ. |", "|" + "---|" * 12]
    for run in runs:
        lines.append(
            f"| {run.site}, {run.case} | {run.electricity_kWh:.1f} | {run.heat_kWh:.1f} | "
            f"{run.cop_in_use:.2f} | {scop_on_model[run.site]:.2f} | {run.defrosts} | "
            f"{run.defrost_h:.1f} | {100 * run.defrost_penalty:.1f} | {run.frost_max_kg:.2f} | "
            f"{100 * run.share_s_max:.0f} | {run.starts} | {run.Kh_below_occ:.1f} |")
    table = "\n".join(lines)
    print("Heating season 1 October - 30 April, typical meteorological years\n")
    print(table)
    print()
    print(scop_table)
    print()
    print(points_table)
    (FIG / "heating_comparison.md").write_text(table + "\n", encoding="utf-8", newline="\n")
    (FIG / "scop.md").write_text(scop_table + "\n\n" + points_table + "\n", encoding="utf-8",
                                 newline="\n")
    main_cases = [case for case, _, _ in CASES]
    monthly = {(run.site, run.case): (run.monthly_electricity_kWh, run.monthly_heat_kWh)
               for run in runs if run.case in main_cases}
    plot_season_monthly(monthly, list(SITES), main_cases, [m for _, m in MONTHS],
                        FIG / "heating_monthly.png",
                        title="Heating season by month: electricity of the unit and heat "
                              "delivered per kWh (not shown below 1 kWh)",
                        ratio_label="COP in use (W/W)")
    print(f"\nFigures written to {FIG}")


if __name__ == "__main__":
    main()
