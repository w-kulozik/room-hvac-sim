"""Compare the standard controllers on three summer days.

Runs on/off and PI controllers (low-load handling by hysteresis, A, or sigma-delta, C), with the
indoor fan either following the compressor or left running between cycles. Each run has three
warm-up days that are discarded, so the controller and the thermal mass start the reported days
in their periodic regime.

    python scripts/compare_controllers.py
"""

from __future__ import annotations

from pathlib import Path

from roomsim.ac_unit import ACUnit
from roomsim.controllers import OnOffThermostat, VentilationSchedule, make_pi_controller
from roomsim.inputs import DAY, SummerScenario
from roomsim.metrics import Metrics, compute
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.plotting import apply_style, plot_controller_day, plot_controller_metrics
from roomsim.results import SimResult
from roomsim.simulate import Controller, simulate

FIG = Path(__file__).resolve().parents[1] / "figures"
T_SET = 25.0
WARMUP_DAYS, DAYS = 3, 3


def main() -> None:
    """Run the controllers, write figures/controller_comparison.md and the two figures."""
    FIG.mkdir(exist_ok=True)
    apply_style()
    p, unit, S, vent = RoomParams(), ACUnit(), SummerScenario(), VentilationSchedule()
    x0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))

    controllers: dict[str, Controller] = {
        "On/off thermostat": OnOffThermostat(T_set=T_SET, ventilation=vent),
        "On/off, fan continuous": OnOffThermostat(T_set=T_SET, ventilation=vent,
                                                  fan_mode="continuous"),
        "PI + hysteresis (A)": make_pi_controller("hysteresis", p, unit, T_set=T_SET,
                                                  ventilation=vent),
        "PI + sigma-delta (C)": make_pi_controller("sigma-delta", p, unit, T_set=T_SET,
                                                   ventilation=vent),
        "PI (C), fan continuous": make_pi_controller("sigma-delta", p, unit, T_set=T_SET,
                                                     ventilation=vent, fan_mode="continuous"),
    }
    results: dict[str, SimResult] = {}
    metrics: dict[str, Metrics] = {}
    for name, ctrl in controllers.items():
        r = simulate(p, S, x0, (WARMUP_DAYS + DAYS) * DAY, controller=ctrl, unit=unit)
        results[name] = r.window(WARMUP_DAYS * DAY)
        metrics[name] = compute(results[name], T_SET)

    header = ("| Controller | Energy (kWh) | Starts | Runtime (h) | T rms occ. (K) | "
              "K·h outside ±1 K | RH mean occ. (%) | h RH > 60 % occ. |")
    lines = [header, "|" + "---|" * 8]
    for name, m in metrics.items():
        lines.append(f"| {name} | {m.energy_kWh:.2f} | {m.starts} | {m.runtime_h:.1f} | "
                     f"{m.T_rms_occ:.2f} | {m.Kh_outside_occ:.2f} | {m.RH_mean_occ:.1f} | "
                     f"{m.h_RH60_occ:.1f} |")
    table = "\n".join(lines)
    print(f"Three summer days after {WARMUP_DAYS} warm-up days, set point {T_SET} °C\n")
    print(table)
    (FIG / "controller_comparison.md").write_text(table + "\n", encoding="utf-8", newline="\n")

    last_day = {n: r.window((WARMUP_DAYS + DAYS - 1) * DAY) for n, r in results.items()}
    plot_controller_day(last_day, T_SET, FIG / "controllers_day.png")
    plot_controller_metrics(metrics, FIG / "controllers_metrics.png")
    print(f"\nFigures written to {FIG}")


if __name__ == "__main__":
    main()
