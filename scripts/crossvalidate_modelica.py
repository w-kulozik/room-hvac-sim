"""Cross-validation of the room model against an independent Modelica implementation (D58).

Compares the Python model (RK4, 60 s) with the Modelica model in modelica/RoomSim.mo, solved by
OpenModelica (DASSL, relative tolerance 1e-10), on two three-day cases: A free-floating, and B
with ventilation and fixed cooling and dehumidification while occupied. The Modelica results
are committed in modelica/results/; with --run-modelica they are regenerated first (needs
OpenModelica; the omc executable is taken from $OPENMODELICAHOME/bin or the PATH).

    python scripts/crossvalidate_modelica.py [--run-modelica]
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
from collections.abc import Callable
from pathlib import Path

import numpy as np

from roomsim._typing import FloatArray
from roomsim.inputs import DAY, HOUR, SummerScenario
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.simulate import simulate

ROOT = Path(__file__).resolve().parents[1]
MODELICA = ROOT / "modelica"
RESULTS = MODELICA / "results"
FIG = ROOT / "figures"
DURATION, OUTPUT_STEP = 3 * DAY, 600.0
CASES = {"A free-floating": "FreeFloating", "B scheduled": "Scheduled"}
STATES = ("T_a", "T_m", "W", "phi_b", "C")
# acceptance bounds of D58, in the units of each state
BOUNDS = {"T_a": 1e-5, "T_m": 1e-5, "W": 1e-8, "phi_b": 1e-6, "C": 1e-4}
UNITS = {"T_a": "K", "T_m": "K", "W": "kg/kg", "phi_b": "-", "C": "ppm"}


class Scheduled:
    """Case B: ventilation, 800 W of cooling and 1e-4 kg/s of water removed from 8 to 18 h."""

    def step(self, t: float, y: FloatArray) -> FloatArray:
        """Return u = [Q_hvac_sens W, m_hvac_lat kg/s, V_vent m3/s] at time t (s)."""
        h = (t % DAY) / HOUR
        if 8.0 <= h < 18.0:
            return np.array([-800.0, 1e-4, 0.0105])
        return np.zeros(3)


def python_run(case: str) -> dict[str, FloatArray]:
    """Return the Python states of a case at the output times."""
    p, S = RoomParams(), SummerScenario()
    x0 = equilibrium_state(26.0, 26.0, float(S.W_o(0.0)))
    controller = Scheduled() if case == "B scheduled" else None
    r = simulate(p, S, x0, DURATION, controller=controller)
    keep = np.isclose(np.mod(r.t, OUTPUT_STEP), 0.0)
    return {"time": r.t[keep], **{name: r.x[keep, i] for i, name in enumerate(STATES)}}


def read_modelica(model: str) -> dict[str, FloatArray]:
    """Return the committed Modelica results of a model at the output times."""
    path = RESULTS / f"RoomSim.{model}_res.csv"
    header = path.read_text(encoding="utf-8").splitlines()[0]
    names = [n.strip('"') for n in header.split(",")]
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    columns = {n: data[:, i] for i, n in enumerate(names)}
    # OpenModelica writes the final time twice (an event at the stop time); keep the first
    t, first = np.unique(columns["time"], return_index=True)
    out = {"time": t}
    for name in STATES:
        out[name] = columns[f"room.{name}"][first]
    return out


def run_modelica() -> None:
    """Simulate both cases with OpenModelica and copy the results to modelica/results/."""
    home = os.environ.get("OPENMODELICAHOME")
    omc = shutil.which("omc", path=str(Path(home) / "bin") if home else None) or shutil.which("omc")
    if omc is None:
        raise SystemExit("omc not found: install OpenModelica or set OPENMODELICAHOME")
    RESULTS.mkdir(exist_ok=True)
    filt = "time|room\\\\.T_a|room\\\\.T_m|room\\\\.W|room\\\\.phi_b|room\\\\.C"
    lines = ["loadModel(Modelica); getErrorString();",
             f'loadFile("{(MODELICA / "RoomSim.mo").as_posix()}"); getErrorString();']
    for model in CASES.values():
        lines.append(f'simulate(RoomSim.{model}, stopTime={DURATION:.0f}, '
                     f'numberOfIntervals={DURATION / OUTPUT_STEP:.0f}, tolerance=1e-10, '
                     f'outputFormat="csv", variableFilter="{filt}"); getErrorString();')
    lines.append("getVersion();")
    with tempfile.TemporaryDirectory() as tmp:
        script = Path(tmp) / "run.mos"
        script.write_text("\n".join(lines) + "\n", encoding="utf-8")
        proc = subprocess.run([omc, str(script)], cwd=tmp, capture_output=True, text=True,
                              check=False)
        print(proc.stdout[-2000:])
        for model in CASES.values():
            src = Path(tmp) / f"RoomSim.{model}_res.csv"
            if not src.exists():
                raise SystemExit(f"OpenModelica did not produce {src.name}:\n{proc.stdout}")
            text = src.read_text(encoding="utf-8").replace("\r\n", "\n")
            (RESULTS / src.name).write_text(text, encoding="utf-8", newline="\n")
        version = proc.stdout.strip().splitlines()[-1].strip('"')
    (RESULTS / "VERSION.txt").write_text(version + "\n", encoding="utf-8", newline="\n")


def compare() -> tuple[str, dict[str, dict[str, float]]]:
    """Return the comparison table and the largest differences per case and state."""
    lines = ["| Case | " + " | ".join(f"max |Δ{n}| ({UNITS[n]})" for n in STATES) + " |",
             "|" + "---|" * (len(STATES) + 1)]
    diffs: dict[str, dict[str, float]] = {}
    fmt: Callable[[float], str] = lambda v: f"{v:.1e}"  # noqa: E731  one-line formatter
    for case, model in CASES.items():
        py, mo = python_run(case), read_modelica(model)
        if not np.array_equal(py["time"], mo["time"]):
            raise SystemExit(f"{case}: output times differ between Python and Modelica")
        diffs[case] = {n: float(np.max(np.abs(py[n] - mo[n]))) for n in STATES}
        lines.append(f"| {case} | " + " | ".join(fmt(diffs[case][n]) for n in STATES) + " |")
    lines.append("| Bound (D58) | " + " | ".join(fmt(BOUNDS[n]) for n in STATES) + " |")
    return "\n".join(lines), diffs


def main() -> None:
    """Optionally rerun Modelica; compare; write figures/modelica_crossvalidation.md."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run-modelica", action="store_true",
                        help="regenerate the Modelica results with OpenModelica first")
    args = parser.parse_args()
    if args.run_modelica:
        run_modelica()
    table, diffs = compare()
    version = (RESULTS / "VERSION.txt").read_text(encoding="utf-8").strip()
    print(f"Python RK4 at 60 s against {version}, DASSL, tolerance 1e-10\n")
    print(table)
    (FIG / "modelica_crossvalidation.md").write_text(table + "\n", encoding="utf-8",
                                                     newline="\n")
    failed = [f"{c}: {n}" for c, d in diffs.items() for n in STATES if d[n] >= BOUNDS[n]]
    if failed:
        raise SystemExit("outside the bounds of D58: " + ", ".join(failed))


if __name__ == "__main__":
    main()
