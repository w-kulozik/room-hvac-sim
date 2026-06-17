"""Rated and part-load points of the unit model against published data of a real inverter split.

Reference: Daikin FTXM25R + RXM25R (2.5 kW wall-mounted inverter split, R-32), values from the
manufacturer's technical data books (full references in docs/references.md):
    RXM-R technical data book EEDEN20, pages 19-21 (combination FTXM25R + RXM25R)
    FTXM-R technical data EEDDE22, printed page 6, page 3 of the PDF file (indoor unit)
Ratings are at 27 degC dry bulb / 19 degC wet bulb indoors and 35 degC outdoors. The part-load
points B-D are the declared points of the seasonal (SEER) rating at 30, 25 and 20 degC outdoors.

The model is evaluated at the same indoor state with the indoor fan at full speed. Its part-load
capacities are taken in the same proportion to its own rated capacity as the published ones
(74 % at B, 47.6 % at C). At D the published capacity (1.17 kW) is close to the unit's minimum
(1.30 kW at rating conditions), so the model is evaluated at its minimum speed. No cycling
degradation is applied to the model.

    python scripts/compare_rated.py
"""

from __future__ import annotations

from pathlib import Path

from scipy.optimize import brentq

from roomsim.ac_unit import CP_A, CP_V, H_FG0, ACUnit, _h, _W_sat, net_cooling

FIG = Path(__file__).resolve().parents[1] / "figures"
T_IN, T_WB = 27.0, 19.0

# Published values, FTXM25R + RXM25R (see the module docstring)
PUBLISHED = {
    "cap_nom_kW": 2.50, "cap_min_kW": 1.30, "cap_max_kW": 3.20, "P_nom_kW": 0.56, "EER": 4.50,
    "SEER": 8.65, "airflow_m3h": (630, 456, 342, 246),      # high, medium, low, silent (cooling)
    "P_indoor_W": 25.0, "P_thermostat_off_W": 6.0,
    # declared part-load points (outdoor degC, Pdc kW, EERd)
    "points": (("A", 35.0, 2.50, 4.50), ("B", 30.0, 1.85, 6.52), ("C", 25.0, 1.19, 10.17),
               ("D", 20.0, 1.17, 16.51)),
}


def W_at_wet_bulb(T: float, T_wb: float) -> float:
    """Humidity ratio of air at dry bulb T on the constant-enthalpy line of wet bulb T_wb."""
    return (_h(T_wb, _W_sat(T_wb)) - CP_A * T) / (H_FG0 + CP_V * T)


def main() -> None:
    FIG.mkdir(exist_ok=True)
    unit = ACUnit()
    p = unit.p
    W_in = W_at_wet_bulb(T_IN, T_WB)

    def net(s: float, T_o: float) -> float:
        return net_cooling(unit.operate(s, 1.0, T_IN, W_in, T_o), T_IN)

    def eer(s: float, T_o: float) -> float:
        return net(s, T_o) / unit.operate(s, 1.0, T_IN, W_in, T_o).P_el

    rated = net(1.0, 35.0)
    rows = ["| Point | Outdoor (°C) | Published capacity (kW) | Published EER | "
            "Model capacity (kW) | Model EER |", "|" + "---|" * 6]
    for name, T_o, Pdc, EERd in PUBLISHED["points"]:
        if name == "D":
            s = p.s_min
        else:
            target = rated * Pdc / PUBLISHED["cap_nom_kW"]
            s = brentq(lambda x: net(x, T_o) - target, p.s_min, 1.0) if name != "A" else 1.0
        rows.append(f"| {name} | {T_o:.0f} | {Pdc:.2f} | {EERd:.2f} | {net(s, T_o) / 1e3:.2f} | "
                    f"{eer(s, T_o):.2f} |")
    points = "\n".join(rows)

    st = unit.operate(1.0, 1.0, T_IN, W_in, 35.0)
    f_low = p.f_min
    rating = [
        "| Quantity | Published | Model |", "|---|---|---|",
        f"| Net cooling capacity at rating (kW) | {PUBLISHED['cap_nom_kW']:.2f} | "
        f"{rated / 1e3:.2f} |",
        f"| Electrical power at rating (kW) | {PUBLISHED['P_nom_kW']:.2f} | {st.P_el / 1e3:.2f} |",
        f"| Net EER at rating | {PUBLISHED['EER']:.2f} | {rated / st.P_el:.2f} |",
        f"| Minimum capacity at rating conditions (kW) | {PUBLISHED['cap_min_kW']:.2f} | "
        f"{net(p.s_min, 35.0) / 1e3:.2f} |",
        f"| Maximum capacity (kW) | {PUBLISHED['cap_max_kW']:.2f} | {rated / 1e3:.2f} (no boost) |",
        f"| SEER | {PUBLISHED['SEER']:.2f} | not computed |",
        f"| Indoor airflow, highest / lowest step (m³/h) | {PUBLISHED['airflow_m3h'][0]} / "
        f"{PUBLISHED['airflow_m3h'][3]} | {3600 * p.V_nom:.0f} / {3600 * p.V_nom * f_low:.0f} |",
        f"| Indoor unit power, cooling nominal (W) | {PUBLISHED['P_indoor_W']:.0f} | "
        f"{unit.fan_power(1.0):.0f} (full fan) |",
        f"| Power with the compressor off, fan running (W) | {PUBLISHED['P_thermostat_off_W']:.0f} "
        f"(thermostat-off mode) | {unit.fan_power(f_low):.1f} (fan at minimum) |",
    ]
    table = "\n".join(rating)
    print(f"Indoor {T_IN:.0f} °C / {T_WB:.0f} °C wet bulb (W = {W_in:.4f}), full indoor fan\n")
    print(table)
    print()
    print(points)
    (FIG / "rated_comparison.md").write_text(table + "\n\n" + points + "\n", encoding="utf-8",
                                             newline="\n")


if __name__ == "__main__":
    main()
