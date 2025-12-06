"""Day-1 deliverable: three summer days without air conditioning.

Produces figures/temperatures.png, figures/humidity.png, figures/energy_balance.png and
figures/convergence.png, and prints the free-floating baseline indicators.

    python scripts/run_baseline.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from roomsim import verification as v
from roomsim.inputs import DAY, HOUR, SummerScenario
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.plotting import (apply_style, plot_convergence, plot_energy_balance,
                              plot_humidity, plot_temperatures)
from roomsim.simulate import periodic_state, simulate

FIG = Path(__file__).resolve().parents[1] / "figures"


def daily_table(r, label: str) -> None:
    print(f"\n{label}")
    print("  day  T_a max (time)      T_a min   T_m mean   RH max")
    for k in range(int(round((r.t[-1] - r.t[0]) / DAY))):
        w = r.window(r.t[0] + k * DAY, r.t[0] + (k + 1) * DAY)
        i = int(np.argmax(w.T_a))
        print(f"  {k + 1:>3}  {w.T_a[i]:6.2f} °C ({(w.t[i] % DAY) / HOUR:5.2f} h)"
              f"  {w.T_a.min():6.2f}    {w.T_m.mean():6.2f}    {100 * w.RH.max():5.1f} %")


def main() -> None:
    FIG.mkdir(exist_ok=True)
    apply_style()
    p, S = RoomParams(), SummerScenario()
    x_guess = equilibrium_state(S.T_mean, S.T_mean, float(S.W_o(0.0)))

    # Reported run: starts on the periodic steady state (Newton shooting)
    x_p = periodic_state(p, S, x_guess)
    r = simulate(p, S, x_p, 3 * DAY)
    # For comparison: naive start with everything at the mean outdoor temperature
    r_cold = simulate(p, S, x_guess, 3 * DAY)

    daily_table(r, "Periodic steady state (reported baseline)")
    daily_table(r_cold, "Naive start at 26 °C (not periodic, for comparison)")
    k = int(np.argmax(r.T_a))
    print(f"\nFree-floating baseline: T_a,max = {r.T_a[k]:.2f} °C at {(r.t[k] % DAY) / HOUR:.2f} h")
    print(f"Energy residual max: {np.abs(r.energy_residual).max():.2e} J;"
          f" water residual max: {np.abs(r.water_residual).max():.2e} kg")

    plot_temperatures(r, FIG / "temperatures.png")
    plot_humidity(r, FIG / "humidity.png")
    plot_energy_balance(r, FIG / "energy_balance.png")

    dts1 = np.array([4000.0, 2000.0, 1000.0, 500.0, 250.0])
    dts2 = np.array([1200.0, 600.0, 300.0, 150.0])
    plot_convergence({
        "One-node step response": (dts1, {m: [v.one_node_error(m, dt) for dt in dts1]
                                          for m in ("rk4", "euler")}, "|T_a − exact| (K)"),
        "Energy balance vs Simpson": (dts2, {m: [v.energy_residual_simpson(m, dt) for dt in dts2]
                                             for m in ("rk4", "euler")}, "|residual| (J)"),
    }, FIG / "convergence.png")
    print(f"Figures written to {FIG}")


if __name__ == "__main__":
    main()
