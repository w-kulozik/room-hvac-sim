"""Ventilation strategies (decisions D23-D27, fixed before these runs).

Every strategy runs with the same temperature and humidity control (PI with sigma-delta low-load
handling and the coil-temperature cascade on the fan), so differences come from ventilation only.
Energy is reported for the air conditioner and for the ventilation fan separately.

Two occupancy scenarios: two people throughout 8-18 h (the design case), and a variable office
day (D27). A last, exploratory table varies the night-flush flow.

    python scripts/compare_ventilation.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from roomsim.ac_unit import ACUnit
from roomsim.controllers import (
    CoilCascade,
    DemandControlledPI,
    NightFlush,
    ScheduledVentilation,
    StagedVentilation,
    VentilationSchedule,
    WithVentilation,
    make_pi_controller,
)
from roomsim.inputs import DAY, HOUR, VARIABLE_OCCUPANCY, SummerScenario
from roomsim.metrics import compute
from roomsim.params import RoomParams
from roomsim.plant import IDX_V_VENT, equilibrium_state
from roomsim.plotting import apply_style, plot_ventilation_day
from roomsim.simulate import simulate
from roomsim.ventilation import VentilationFan

FIG = Path(__file__).resolve().parents[1] / "figures"
WARMUP_DAYS, DAYS = 3, 3
CO2_LIMIT = 1000.0
P, UNIT, FAN = RoomParams(), ACUnit(), VentilationFan()


def thermal():
    return CoilCascade(make_pi_controller("sigma-delta", P, UNIT), RH_set=0.55)


def run(scenario, strategy):
    x0 = equilibrium_state(26.0, 26.0, float(scenario.W_o(0.0)))
    r = simulate(P, scenario, x0, (WARMUP_DAYS + DAYS) * DAY,
                 controller=WithVentilation(thermal(), strategy), unit=UNIT)
    return r.window(WARMUP_DAYS * DAY)


def row(name, r):
    m = compute(r, 25.0)
    occ = np.asarray(r.d.n_occ)[:-1] > 0
    C = r.CO2[:-1]
    E_fan = FAN.energy_kWh(r)
    air = float(np.sum(r.u[:-1, IDX_V_VENT]) * r.dt)
    return (f"| {name} | {m.energy_kWh:.2f} | {E_fan:.3f} | {m.energy_kWh + E_fan:.2f} | "
            f"{air:.0f} | {C[occ].mean():.0f} | {C[occ].max():.0f} | "
            f"{np.sum(C[occ] > CO2_LIMIT) * r.dt / HOUR:.2f} | {m.Kh_outside_occ:.3f} |")


HEADER = ["| Strategy | AC (kWh) | Fan (kWh) | Total (kWh) | Air (m³) | CO₂ mean occ. (ppm) | "
          "CO₂ max occ. (ppm) | h > 1000 ppm | K·h outside band |", "|" + "---|" * 9]


def strategies(S):
    return {
        "No ventilation (infiltration only)": ScheduledVentilation(VentilationSchedule(V_occupied=0.0)),
        "V0 fixed schedule": ScheduledVentilation(VentilationSchedule()),
        "V1 demand-controlled PI": DemandControlledPI(P),
        "V2 staged levels": StagedVentilation(),
        "V3 V1 + night flush": NightFlush(DemandControlledPI(P), outdoor=S.T_o),
    }


def main() -> None:
    FIG.mkdir(exist_ok=True)
    apply_style()
    out = []
    scenarios = {"Two people, 8-18 h (design case)": SummerScenario(),
                 "Variable occupancy (D27)": SummerScenario(occupancy=VARIABLE_OCCUPANCY)}
    last_results = None
    for label, S in scenarios.items():
        results = {name: run(S, strat) for name, strat in strategies(S).items()}
        out += [f"### {label}", ""] + HEADER + [row(n, r) for n, r in results.items()] + [""]
        last_results = results

    S = scenarios["Variable occupancy (D27)"]
    # exploratory: a fixed flow large enough for the 3-person peak (CO2 at 1000 ppm in steady state)
    V_peak = 3 * P.co2_person * 1e6 / (CO2_LIMIT - 420.0) - P.V_dot_inf
    out += ["### Fixed schedule sized for the peak (exploratory, variable occupancy)", ""] + HEADER
    out.append(row(f"V0 sized for 3 people ({V_peak * 3600:.0f} m³/h)",
                   run(S, ScheduledVentilation(VentilationSchedule(V_occupied=V_peak)))))
    out += ["", "### Night-flush flow (exploratory, design case)", ""] + HEADER
    for V in (0.01, 0.015, 0.02, 0.03):
        r = run(SummerScenario(), NightFlush(DemandControlledPI(P), outdoor=S.T_o, V_flush=V))
        out.append(row(f"V3 with flush at {V * 3600:.0f} m³/h", r))
    text = "\n".join(out)
    sys.stdout.reconfigure(encoding="utf-8")    # the table is not cp1252 (subscript 2 in CO2)
    print(text)
    (FIG / "ventilation_comparison.md").write_text(text + "\n", encoding="utf-8", newline="\n")

    shown = ("V0 fixed schedule", "V1 demand-controlled PI", "V2 staged levels", "V3 V1 + night flush")
    last_day = {n: last_results[n].window((WARMUP_DAYS + DAYS - 1) * DAY) for n in shown}
    plot_ventilation_day(last_day, CO2_LIMIT, FIG / "ventilation_day.png")
    print(f"\nFigures written to {FIG}")


if __name__ == "__main__":
    main()
