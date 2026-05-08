"""Ventilation strategies (decisions D23-D26, fixed before these runs).

Every strategy runs with the same temperature and humidity control (PI with sigma-delta low-load
handling and the coil-temperature cascade on the fan), so differences come from ventilation only.
Energy is reported for the air conditioner and for the ventilation fan separately.

    python scripts/compare_ventilation.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from roomsim.ac_unit import ACUnit
from roomsim.controllers import (CoilCascade, DemandControlledPI, NightFlush,
                                 ScheduledVentilation, StagedVentilation, VentilationSchedule,
                                 WithVentilation, make_pi_controller)
from roomsim.inputs import DAY, HOUR, SummerScenario
from roomsim.metrics import compute
from roomsim.params import RoomParams
from roomsim.plant import IDX_V_VENT, equilibrium_state
from roomsim.plotting import apply_style, plot_ventilation_day
from roomsim.simulate import simulate
from roomsim.ventilation import VentilationFan

FIG = Path(__file__).resolve().parents[1] / "figures"
WARMUP_DAYS, DAYS = 3, 3
CO2_LIMIT = 1000.0


def main() -> None:
    FIG.mkdir(exist_ok=True)
    apply_style()
    p, unit, S, fan = RoomParams(), ACUnit(), SummerScenario(), VentilationFan()
    x0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))

    def thermal():
        return CoilCascade(make_pi_controller("sigma-delta", p, unit), RH_set=0.55)

    strategies = {
        "No ventilation (infiltration only)": ScheduledVentilation(VentilationSchedule(V_occupied=0.0)),
        "V0 fixed schedule": ScheduledVentilation(VentilationSchedule()),
        "V1 demand-controlled PI": DemandControlledPI(p),
        "V2 staged levels": StagedVentilation(),
        "V3 V1 + night flush": NightFlush(DemandControlledPI(p), outdoor=S.T_o),
    }
    results = {}
    for name, strategy in strategies.items():
        r = simulate(p, S, x0, (WARMUP_DAYS + DAYS) * DAY,
                     controller=WithVentilation(thermal(), strategy), unit=unit)
        results[name] = r.window(WARMUP_DAYS * DAY)

    lines = ["| Strategy | AC (kWh) | Fan (kWh) | Total (kWh) | Air (m³) | CO₂ mean occ. (ppm) | "
             "CO₂ max occ. (ppm) | h > 1000 ppm | RH mean occ. (%) | K·h outside band |",
             "|" + "---|" * 10]
    for name, r in results.items():
        m = compute(r, 25.0)
        occ = np.asarray(r.d.n_occ)[:-1] > 0
        C = r.CO2[:-1]
        E_fan = fan.energy_kWh(r)
        air = float(np.sum(r.u[:-1, IDX_V_VENT]) * r.dt)
        lines.append(f"| {name} | {m.energy_kWh:.2f} | {E_fan:.3f} | {m.energy_kWh + E_fan:.2f} | "
                     f"{air:.0f} | {C[occ].mean():.0f} | {C[occ].max():.0f} | "
                     f"{np.sum(C[occ] > CO2_LIMIT) * r.dt / HOUR:.2f} | {m.RH_mean_occ:.1f} | "
                     f"{m.Kh_outside_occ:.3f} |")
    table = "\n".join(lines)
    print(f"Three days after {WARMUP_DAYS} warm-up days\n")
    print(table)
    (FIG / "ventilation_comparison.md").write_text(table + "\n", encoding="utf-8")

    shown = ("V0 fixed schedule", "V1 demand-controlled PI", "V2 staged levels", "V3 V1 + night flush")
    last_day = {n: results[n].window((WARMUP_DAYS + DAYS - 1) * DAY) for n in shown}
    plot_ventilation_day(last_day, CO2_LIMIT, FIG / "ventilation_day.png")
    print(f"\nFigures written to {FIG}")


if __name__ == "__main__":
    main()
