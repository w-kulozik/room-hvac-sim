"""Economic MPC against PI control, with a time-of-use tariff and robustness cases.

All controllers keep the comfort band 25 +- 1 degC while the office is occupied (8-18 h); outside
those hours the room may float up to 28 degC. The PI references are the constant set point and a
night set-back with a fixed early start at 6:00. The MPC is run with its nominal model, with a
time-of-use tariff (peak 14-20 h at four times the off-peak price), with a model whose parameters
are off by 20 % and with a biased weather forecast.

    python scripts/compare_mpc.py
"""

from __future__ import annotations

from pathlib import Path

from roomsim.ac_unit import ACUnit
from roomsim.control_model import RCModel
from roomsim.controllers import (ComfortSchedule, EconomicMPC, VentilationSchedule,
                                 make_pi_controller)
from roomsim.inputs import DAY, SummerScenario
from roomsim.metrics import compute, energy_cost
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.plotting import apply_style, plot_mpc_day
from roomsim.simulate import simulate
from roomsim.tariff import TimeOfUseTariff

FIG = Path(__file__).resolve().parents[1] / "figures"
WARMUP_DAYS, DAYS = 3, 3
MISMATCH = {"C_m": 1.2, "furniture_factor": 0.8, "H_mo": 0.8, "h_i": 1.2, "g_win": 0.8}


def main() -> None:
    FIG.mkdir(exist_ok=True)
    apply_style()
    p, unit, S, vent = RoomParams(), ACUnit(), SummerScenario(), VentilationSchedule()
    comfort = ComfortSchedule()
    tou = TimeOfUseTariff(offpeak=0.15, peak=0.60)
    x0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
    model = RCModel(p)

    controllers = {
        "PI, 25 °C constant": make_pi_controller("sigma-delta", p, unit, ventilation=vent),
        "PI, set-back, start 6:00": make_pi_controller(
            "sigma-delta", p, unit, ventilation=vent,
            T_set=ComfortSchedule(early_start_h=6.0).setpoint),
        "MPC": EconomicMPC(model, unit, S, comfort, ventilation=vent),
        "MPC, time-of-use tariff": EconomicMPC(model, unit, S, comfort, ventilation=vent, tariff=tou),
        "MPC, model off by 20 %": EconomicMPC(RCModel.mismatched(p, MISMATCH), unit, S, comfort,
                                              ventilation=vent),
        "MPC, forecast +2 K, -30 % sun": EconomicMPC(
            model, unit, SummerScenario(T_mean=28.0, G_max=350.0), comfort, ventilation=vent),
    }
    results = {}
    for name, ctrl in controllers.items():
        r = simulate(p, S, x0, (WARMUP_DAYS + DAYS) * DAY, controller=ctrl, unit=unit)
        results[name] = r.window(WARMUP_DAYS * DAY)

    lines = ["| Controller | Energy (kWh) | Cost, ToU (EUR) | K·h outside band | T max occ. (°C) | "
             "Starts |", "|" + "---|" * 6]
    for name, r in results.items():
        m = compute(r, comfort.T_set)
        lines.append(f"| {name} | {m.energy_kWh:.2f} | {energy_cost(r, tou):.2f} | "
                     f"{m.Kh_outside_occ:.3f} | {m.T_max_occ:.2f} | {m.starts} |")
    table = "\n".join(lines)
    print(f"Three days after {WARMUP_DAYS} warm-up days; tariff {tou.offpeak:.2f} / {tou.peak:.2f} EUR/kWh\n")
    print(table)
    (FIG / "mpc_comparison.md").write_text(table + "\n", encoding="utf-8", newline="\n")

    shown = ("PI, set-back, start 6:00", "MPC", "MPC, time-of-use tariff")
    last_day = {n: results[n].window((WARMUP_DAYS + DAYS - 1) * DAY) for n in shown}
    plot_mpc_day(last_day, comfort, tou, FIG / "mpc_day.png")
    print(f"\nFigures written to {FIG}")


if __name__ == "__main__":
    main()
