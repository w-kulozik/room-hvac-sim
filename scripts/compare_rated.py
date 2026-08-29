"""Rated and part-load points of the unit model against published data of a real inverter split.

Reference: Daikin FTXM25R + RXM25R, a 2.5 kW inverter split; the published values and their
sources are in ``roomsim.reference_unit``. Ratings are at 27 degC dry bulb / 19 degC wet bulb
indoors and 35 degC outdoors. The part-load points B-D are the declared points of the seasonal
(SEER) rating at 30, 25 and 20 degC outdoors.

The model is evaluated at the same indoor state with the indoor fan at full speed. Its part-load
capacities are taken in the same proportion to its own rated capacity as the published ones
(74 % at B, 47.6 % at C). At D the published capacity (1.17 kW) is close to the unit's minimum
(1.30 kW at rating conditions), so the model is evaluated at its minimum speed. No cycling
degradation is applied to the model in this table; the model's SEER (decision D46) is computed
by the bin method, which applies it.

    python scripts/compare_rated.py
"""

from __future__ import annotations

from pathlib import Path

from scipy.optimize import brentq

from roomsim.ac_unit import ACUnit, net_cooling
from roomsim.reference_unit import PUBLISHED
from roomsim.seasonal import W_at_wet_bulb, seer, unit_declared_points

FIG = Path(__file__).resolve().parents[1] / "figures"
T_IN, T_WB = 27.0, 19.0

def main() -> None:
    """Evaluate the unit model at the published points; write figures/rated_comparison.md."""
    FIG.mkdir(exist_ok=True)
    unit = ACUnit()
    p = unit.p
    W_in = W_at_wet_bulb(T_IN, T_WB)

    def net(s: float, T_o: float) -> float:
        # net cooling at the indoor rating state, full fan, compressor speed s, outdoor T_o; W
        return net_cooling(unit.operate(s, 1.0, T_IN, W_in, T_o), T_IN)

    def eer(s: float, T_o: float) -> float:
        return net(s, T_o) / unit.operate(s, 1.0, T_IN, W_in, T_o).P_el

    rated = net(1.0, 35.0)
    rows = ["| Point | Outdoor (°C) | Published capacity (kW) | Published EER | "
            "Model capacity (kW) | Model EER |", "|" + "---|" * 6]
    for name, T_o, Pdc, EERd in PUBLISHED.points:
        if name == "D":
            s = p.s_min
        else:
            target = rated * Pdc / PUBLISHED.cap_nom_kW
            s = (brentq(lambda x, T_o=T_o, target=target: net(x, T_o) - target, p.s_min, 1.0)
                 if name != "A" else 1.0)
        rows.append(f"| {name} | {T_o:.0f} | {Pdc:.2f} | {EERd:.2f} | {net(s, T_o) / 1e3:.2f} | "
                    f"{eer(s, T_o):.2f} |")
    points = "\n".join(rows)

    st = unit.operate(1.0, 1.0, T_IN, W_in, 35.0)
    P_design, model_points = unit_declared_points(unit)
    model_seer = seer(model_points, P_design, PUBLISHED.low_power_modes()).SEER
    f_low = p.f_min
    rating = [
        "| Quantity | Published | Model |", "|---|---|---|",
        f"| Net cooling capacity at rating (kW) | {PUBLISHED.cap_nom_kW:.2f} | "
        f"{rated / 1e3:.2f} |",
        f"| Electrical power at rating (kW) | {PUBLISHED.P_nom_kW:.2f} | {st.P_el / 1e3:.2f} |",
        f"| Net EER at rating | {PUBLISHED.EER:.2f} | {rated / st.P_el:.2f} |",
        f"| Minimum capacity at rating conditions (kW) | {PUBLISHED.cap_min_kW:.2f} | "
        f"{net(p.s_min, 35.0) / 1e3:.2f} |",
        f"| Maximum capacity (kW) | {PUBLISHED.cap_max_kW:.2f} | {rated / 1e3:.2f} (no boost) |",
        f"| SEER | {PUBLISHED.SEER:.2f} | {model_seer:.2f} (D46) |",
        f"| Indoor airflow, highest / lowest step (m³/h) | {PUBLISHED.airflow_m3h[0]} / "
        f"{PUBLISHED.airflow_m3h[3]} | {3600 * p.V_nom:.0f} / {3600 * p.V_nom * f_low:.0f} |",
        f"| Indoor unit power, cooling nominal (W) | {PUBLISHED.P_indoor_W:.0f} | "
        f"{unit.fan_power(1.0):.0f} (full fan) |",
        f"| Power with the compressor off, fan running (W) | {PUBLISHED.P_thermostat_off_W:.0f} "
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
