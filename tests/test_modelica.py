"""Cross-validation of the room model against the committed Modelica results (decision D58).

Uses the comparison of scripts/crossvalidate_modelica.py with the results in modelica/results/,
so OpenModelica is not needed to run it.
"""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("crossvalidate_modelica",
                                               ROOT / "scripts" / "crossvalidate_modelica.py")
assert _spec is not None and _spec.loader is not None
XV = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(XV)


def test_python_and_modelica_agree_within_the_bounds_of_D58():
    _, diffs = XV.compare()
    for case, d in diffs.items():
        for name, bound in XV.BOUNDS.items():
            assert d[name] < bound, f"{case}: {name} differs by {d[name]:.2e}"


def test_modelica_results_cover_the_three_days():
    for model in XV.CASES.values():
        mo = XV.read_modelica(model)
        assert mo["time"][0] == 0.0 and mo["time"][-1] == XV.DURATION
        assert len(mo["time"]) == XV.DURATION / XV.OUTPUT_STEP + 1
