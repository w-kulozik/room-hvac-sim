"""Regression tests for fixes in the indicators, result windows and plots.

The cases use small hand-built results so that the expected values can be counted by hand.
"""

import importlib.util
import math
import warnings
from pathlib import Path

import numpy as np
import pytest

from roomsim import metrics, plotting
from roomsim.inputs import DAY, HOUR, Disturbances
from roomsim.metrics import compute
from roomsim.params import RoomParams
from roomsim.plant import IDX_TA, IDX_W, N_STATES
from roomsim.results import UNIT_FIELDS, SimResult, UnitLog

DT = 60.0


def _result(s, f=None, n_occ=None, T_a=None):
    """Build a result on a 60 s grid with compressor speeds `s` and fan commands `f`."""
    s = np.asarray(s, dtype=float)
    n = len(s)
    f = np.zeros(n) if f is None else np.asarray(f, dtype=float)
    n_occ = np.ones(n) if n_occ is None else np.asarray(n_occ, dtype=float)
    x = np.zeros((n, N_STATES))
    x[:, IDX_TA] = 25.0 if T_a is None else np.asarray(T_a, dtype=float)
    x[:, IDX_W] = 0.010
    z = np.zeros(n)
    unit = UnitLog(cmd=np.column_stack([s, f, z]), op={k: z.copy() for k in UNIT_FIELDS},
                   M_coil=z.copy(), E_el=z.copy(), M_drain=z.copy(), E_cool=z.copy())
    return SimResult(t=np.arange(n) * DT, x=x, u=np.zeros((n, 3)),
                     d=Disturbances(z + 30.0, z + 0.012, z.copy(), n_occ, z + 420.0), q={},
                     E_in=z.copy(), M_in=z.copy(), CO2_in=z.copy(), W_ext=z.copy(),
                     p=RoomParams(), method="rk4", dt=DT, unit=unit)


# ---- F1: start at the first row of a window ---------------------------------------------

def test_window_does_not_count_a_start_for_a_compressor_already_running():
    """A start is an off -> on transition. When the compressor already runs in the row before
    the window, the first row of the window continues that cycle and is not a start; the
    window must carry the command of that row (``UnitLog.running_before``)."""
    r = _result([1, 1, 0, 0, 1, 1, 0, 0, 0])
    assert compute(r, 25.0).starts == 2          # first row of a run: no history, a start
    w = r.window(1 * DT)                         # rows 1.. : first row continues the cycle
    assert w.require_unit().running_before
    assert compute(w, 25.0).starts == 1
    assert compute(w.window(1 * DT), 25.0).starts == 1      # nested window keeps the flag
    assert compute(r.window(4 * DT), 25.0).starts == 1      # row 3 off: row 4 is a start
    assert compute(r.window(5 * DT), 25.0).starts == 0


# ---- F2 / F5: fan travel ------------------------------------------------------------------

def test_fan_travel_counts_only_steps_while_the_compressor_runs():
    """Fan travel is sum |f_k+1 - f_k| over consecutive rows both running, per running hour.
    The step from the last fan command of one cycle to the first of the next is not travel
    (the fan is off in between), and the last row is a stale command outside the result."""
    r = _result(s=[1, 1, 0, 0, 1, 1, 1], f=[0.5, 0.6, 0.0, 0.0, 0.9, 0.9, 0.2])
    # travel 0.1 (rows 0-1) over 4 running rows of 60 s
    assert metrics.fan_travel_per_h(r) == pytest.approx(0.1 / (4 * DT / HOUR))


def test_fan_travel_is_nan_when_the_compressor_never_runs():
    """With no running hours the indicator (travel per running hour) is undefined: NaN,
    without a division warning."""
    r = _result(s=np.zeros(5), f=np.zeros(5))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert math.isnan(metrics.fan_travel_per_h(r))


# ---- F3: no occupied rows -----------------------------------------------------------------

def test_compute_without_occupied_rows_returns_nan_for_occupied_means():
    """Means, maxima and RMS over an empty set are undefined (NaN); sums over it are zero.
    The energy and compressor indicators stay valid, so a warning rather than an error."""
    r = _result([1, 1, 0, 0, 1, 0], n_occ=np.zeros(6))
    with pytest.warns(RuntimeWarning, match="no occupied") as rec:
        m = compute(r, 25.0)
    assert len(rec) == 1                         # no numpy empty-slice warnings
    for name in ("T_rms_occ", "T_max_occ", "RH_mean_occ", "CO2_mean_occ"):
        assert math.isnan(getattr(m, name))
    assert m.Kh_outside_occ == 0.0 and m.h_RH60_occ == 0.0
    assert m.starts == 2 and m.runtime_h == pytest.approx(3 * DT / HOUR)


# ---- F4: window bounds --------------------------------------------------------------------

@pytest.mark.parametrize("t0, t1", [(10 * DAY, None), (0.0, 10 * DAY), (-DT, 3 * DT),
                                    (3 * DT, 1 * DT), (2 * DT, 2 * DT)])
def test_window_outside_the_run_or_reversed_raises(t0, t1):
    """A window must lie within the run and have t1 > t0; otherwise a clear ValueError
    instead of an IndexError on an empty selection or a silently truncated window."""
    r = _result(np.ones(6))
    with pytest.raises(ValueError, match="t0|t1"):
        r.window(t0, t1)


# ---- F6: masks and plots ------------------------------------------------------------------

def _load_script(name):
    path = Path(__file__).resolve().parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_rh_std_uses_the_same_rows_as_compute():
    """Row k covers [t_k, t_k + dt), so the last row lies outside the result; the RH spread
    must drop it like `compute` does."""
    T = np.full(6, 25.0)
    T[-1] = 15.0                                 # the last row would dominate the spread
    r = _result(np.ones(6), T_a=T)
    assert _load_script("compare_humidity").rh_std_occ(r) == pytest.approx(0.0, abs=1e-12)


def test_peak_spans_for_a_window_starting_or_ending_inside_the_peak():
    """Spans are index pairs [first peak row, first row after it); a peak at the start of the
    window starts at row 0 and a peak running to the end ends at the last row."""
    peak = np.array([1, 1, 0, 0, 1, 1, 0, 1], dtype=bool)
    assert plotting._peak_spans(peak) == [(0, 2), (4, 6), (7, 7)]
    assert plotting._peak_spans(np.zeros(4, dtype=bool)) == []
    assert plotting._peak_spans(np.array([0, 1, 1, 0], dtype=bool)) == [(1, 3)]


def test_step_values_hold_the_last_command_to_the_end():
    """The last row's command acts after the end of the result; a post step plot must end
    on the command of the previous row."""
    assert list(plotting._held(np.array([0.3, 1.0, 0.0]))) == [0.3, 1.0, 1.0]
