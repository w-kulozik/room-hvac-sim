"""Humidity-aware fan control against the plain PI, and robustness to the actual unit.

All humidity controllers keep the compressor under the same PI temperature loop with sigma-delta
low-load handling; they differ only in how they drive the indoor fan:
    PI (C)          fan always at full speed
    H1 fan RH-PI    relative-humidity PI acting directly on the fan
    H2 cascade      RH -> coil temperature set point -> fan (evaporator thermistor)
The robustness runs keep the controllers tuned for the nominal unit and change the unit itself.

    python scripts/compare_humidity.py
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np

from roomsim.ac_unit import ACParams, ACUnit
from roomsim.controllers import (CoilCascade, FanHumidityPI, OnOffThermostat, VentilationSchedule,
                                 make_pi_controller)
from roomsim.inputs import DAY, SummerScenario
from roomsim.metrics import compute
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.plotting import apply_style, plot_humidity_day
from roomsim.simulate import simulate

FIG = Path(__file__).resolve().parents[1] / "figures"
T_SET, RH_SET = 25.0, 0.55
WARMUP_DAYS, DAYS = 3, 3


def rh_std_occ(r) -> float:
    occ = np.asarray(r.d.n_occ) > 0
    return float((100 * r.RH[occ]).std())


def fan_travel_per_h(r) -> float:
    """Total variation of the fan command per running hour: actuator wear and noise."""
    run = r.unit.cmd[:, 0] > 0
    f = r.unit.cmd[run, 1]
    return float(np.abs(np.diff(f)).sum() / (run.sum() * r.dt / 3600.0))


def main() -> None:
    FIG.mkdir(exist_ok=True)
    apply_style()
    p, S, vent = RoomParams(), SummerScenario(), VentilationSchedule()
    nominal = ACUnit()
    x0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))

    def temperature_loop():
        return make_pi_controller("sigma-delta", p, nominal, T_set=T_SET, ventilation=vent)

    builders = {
        "PI (C), full fan": temperature_loop,
        "H1 fan RH-PI": lambda: FanHumidityPI(temperature_loop(), RH_set=RH_SET),
        "H2 cascade": lambda: CoilCascade(temperature_loop(), RH_set=RH_SET),
    }

    def run(build, unit):
        r = simulate(p, S, x0, (WARMUP_DAYS + DAYS) * DAY, controller=build(), unit=unit)
        return r.window(WARMUP_DAYS * DAY)

    results = {name: run(b, nominal) for name, b in builders.items()}
    onoff = run(lambda: OnOffThermostat(T_set=T_SET, ventilation=vent), nominal)

    lines = ["| Controller | Energy (kWh) | Water removed (kg) | RH mean occ. (%) | RH std occ. (%) | "
             "h RH > 60 % occ. | T rms occ. (K) | Starts |", "|" + "---|" * 8]
    for name, r in {**results, "On/off (reference)": onoff}.items():
        m = compute(r, T_SET)
        water = float(r.unit.M_drain[-1])          # condensate drained in the window, kg
        lines.append(f"| {name} | {m.energy_kWh:.2f} | {water:.1f} | {m.RH_mean_occ:.1f} | "
                     f"{rh_std_occ(r):.2f} | {m.h_RH60_occ:.1f} | {m.T_rms_occ:.2f} | {m.starts} |")
    table = "\n".join(lines)

    units = {"nominal": nominal,
             "dirty coil (BF 0.30)": ACUnit(replace(ACParams(), BF_nom=0.30)),
             "larger unit (3.5 kW)": ACUnit(replace(ACParams(), Q_rated=3500.0))}
    rob = ["| Unit | Controller | RH mean occ. (%) | RH std occ. (%) | RH max occ. (%) | "
           "Fan travel per h |", "|" + "---|" * 6]
    for uname, unit in units.items():
        for cname in ("H1 fan RH-PI", "H2 cascade"):
            r = results[cname] if uname == "nominal" else run(builders[cname], unit)
            occ = np.asarray(r.d.n_occ) > 0
            rob.append(f"| {uname} | {cname} | {compute(r, T_SET).RH_mean_occ:.1f} | "
                       f"{rh_std_occ(r):.2f} | {100 * r.RH[occ].max():.1f} | "
                       f"{fan_travel_per_h(r):.2f} |")
    robustness = "\n".join(rob)

    print(f"Set points {T_SET} °C, {100 * RH_SET:.0f} % RH; three days after {WARMUP_DAYS} warm-up days\n")
    print(table)
    print("\nControllers tuned for the nominal unit, run on a different one\n")
    print(robustness)
    (FIG / "humidity_control.md").write_text(table + "\n\n" + robustness + "\n", encoding="utf-8",
                                             newline="\n")

    last_day = {n: r.window((WARMUP_DAYS + DAYS - 1) * DAY) for n, r in results.items()}
    plot_humidity_day(last_day, RH_SET, FIG / "humidity_control.png")
    print(f"\nFigures written to {FIG}")


if __name__ == "__main__":
    main()
