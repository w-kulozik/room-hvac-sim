"""The PI with sigma-delta handling as fixed-point C firmware, against the Python controller (D59).

Open loop: the measurement sequence of a two-day Python closed-loop run is replayed into both
controllers. Closed loop: three standard summer days after three warm-up days with each
controller in the loop. Needs a C99 compiler ($CC, or the zig toolchain: pip install ziglang).

    python scripts/compare_firmware.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from roomsim.ac_unit import ACUnit
from roomsim.controllers import VentilationSchedule, make_pi_controller
from roomsim.controllers.firmware import Q16, FirmwareController, load_library
from roomsim.controllers.pi_unit import PIUnitController
from roomsim.inputs import DAY, SummerScenario
from roomsim.metrics import compute
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.results import SimResult
from roomsim.simulate import Controller, simulate

FIG = Path(__file__).resolve().parents[1] / "figures"
P, UNIT = RoomParams(), ACUnit()
WARMUP_DAYS, DAYS = 3, 3


def python_pi() -> PIUnitController:
    """Return the standard PI with sigma-delta handling (handler C)."""
    return make_pi_controller("sigma-delta", P, UNIT, ventilation=VentilationSchedule())


def run(controller: Controller, days: int) -> SimResult:
    """Simulate the standard summer days with a controller in the loop."""
    S = SummerScenario()
    x0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
    return simulate(P, S, x0, days * DAY, controller=controller, unit=UNIT)


def main() -> None:
    """Compare the controllers; write figures/firmware_sil.md."""
    lib = load_library()
    ref = run(python_pi(), 2)
    py, fw = python_pi(), FirmwareController(python_pi(), lib)
    ds, dv, decisions_differ, samples = 0.0, 0.0, 0, 0
    for k in range(len(ref.t) - 1):
        y = np.array([ref.T_a[k], 0.5, 600.0, 15.0])
        s_py, s_fw = py.step(float(ref.t[k]), y)[0], fw.step(float(ref.t[k]), y)[0]
        ds = max(ds, abs(s_py - s_fw))
        dv = max(dv, abs(py.last_demand - fw.state.demand_q16 / Q16))
        decisions_differ += int((s_py > 0) != (s_fw > 0))
        samples += 1
    open_loop = ["| Open loop, two days | Value |", "|---|---|",
                 f"| Samples | {samples} |",
                 f"| On/off decisions that differ | {decisions_differ} |",
                 f"| Largest speed difference (-) | {ds:.1e} |",
                 f"| Largest demand difference (-) | {dv:.1e} |"]

    total = WARMUP_DAYS + DAYS
    m_py = compute(run(python_pi(), total).window(WARMUP_DAYS * DAY), 25.0)
    m_fw = compute(run(FirmwareController(python_pi(), lib), total).window(WARMUP_DAYS * DAY),
                   25.0)
    closed = ["| Closed loop, three days | Python | Firmware | Difference |", "|---|---|---|---|",
              f"| Energy (kWh) | {m_py.energy_kWh:.3f} | {m_fw.energy_kWh:.3f} | "
              f"{100 * (m_fw.energy_kWh / m_py.energy_kWh - 1):+.3f} % |",
              f"| Starts | {m_py.starts} | {m_fw.starts} | "
              f"{100 * (m_fw.starts / m_py.starts - 1):+.1f} % |",
              f"| T rms occ. (K) | {m_py.T_rms_occ:.3f} | {m_fw.T_rms_occ:.3f} | "
              f"{m_fw.T_rms_occ - m_py.T_rms_occ:+.4f} K |",
              f"| RH mean occ. (%) | {m_py.RH_mean_occ:.2f} | {m_fw.RH_mean_occ:.2f} | "
              f"{m_fw.RH_mean_occ - m_py.RH_mean_occ:+.3f} |"]
    table = "\n".join(open_loop) + "\n\n" + "\n".join(closed)
    print(table)
    (FIG / "firmware_sil.md").write_text(table + "\n", encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
