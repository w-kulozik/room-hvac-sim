"""Performance indicators for comparing controllers on the same scenario."""

from __future__ import annotations

import warnings
from collections.abc import Callable
from dataclasses import asdict, dataclass
from itertools import pairwise

import numpy as np

from ._typing import FloatArray
from .results import SimResult
from .simulate import IDX_CMD_FAN, IDX_CMD_S

HOUR = 3600.0


@dataclass(frozen=True)
class Metrics:
    """Indicators of one controlled run.

    Attributes
    ----------
    energy_kWh : float
        Electrical energy of the unit, kWh.
    starts : int
        Number of compressor starts.
    starts_per_h_on : float
        Starts per hour of compressor operation (0 when the compressor never runs), 1/h.
    T_rms_occ : float
        RMS deviation of the air temperature from the set point, occupied hours, K.
    Kh_outside_occ : float
        Degree-hours outside the set point +- band, occupied hours, K h.
    T_max_occ : float
        Highest air temperature while occupied, degC.
    h_RH60_occ : float
        Hours above 60 % RH while occupied, h.
    RH_mean_occ : float
        Mean relative humidity while occupied, %.
    CO2_mean_occ : float
        Mean CO2 concentration while occupied, ppm.
    runtime_h : float
        Compressor running hours, h.
    """

    energy_kWh: float             # electrical energy of the unit
    starts: int                   # compressor starts
    starts_per_h_on: float        # starts per hour of compressor operation
    T_rms_occ: float              # RMS deviation of T_a from the set point, occupied hours, K
    Kh_outside_occ: float         # degree-hours outside set point +- 1 K, occupied hours
    T_max_occ: float              # highest air temperature while occupied, degC
    h_RH60_occ: float             # hours above 60 % RH while occupied
    RH_mean_occ: float            # mean RH while occupied, %
    CO2_mean_occ: float           # mean CO2 while occupied, ppm
    runtime_h: float              # compressor running hours

    def as_dict(self) -> dict[str, float | int]:
        """Return the indicators as a dictionary keyed by field name.

        Returns
        -------
        dict of str to float or int
            Field names and values, in declaration order.
        """
        return asdict(self)


def compute(r: SimResult, T_set: float, band: float = 1.0) -> Metrics:
    """Compute the indicators over the whole result.

    Row k covers [t_k, t_k + dt), so the last row is dropped. Occupied rows are those with
    ``n_occ > 0``. A start is a row with the compressor running (s > 0) after a row with it
    off; the first row counts as a start unless the compressor was already running on the
    step before it (`UnitLog.running_before`, set when the result is a window of a longer
    run). Without occupied rows, `T_rms_occ`, `T_max_occ`, `RH_mean_occ` and `CO2_mean_occ`
    are NaN, `Kh_outside_occ` and `h_RH60_occ` are 0, and a RuntimeWarning is issued (the
    energy and compressor indicators remain valid).

    Parameters
    ----------
    r : SimResult
        Result of a run with an air-conditioning unit.
    T_set : float
        Temperature set point, degC.
    band : float, optional
        Half-width of the band around `T_set` used for `Kh_outside_occ`, K.

    Returns
    -------
    Metrics
        Indicators of the run.

    Raises
    ------
    ValueError
        If `r` has no unit record.

    Warns
    -----
    RuntimeWarning
        If no row is occupied.
    """
    unit = r.require_unit()
    dt = r.dt
    s = unit.cmd[:-1, IDX_CMD_S]
    on = s > 0
    starts = int(np.sum(on[1:] & ~on[:-1]) + (1 if on[0] and not unit.running_before else 0))
    occ = np.asarray(r.d.n_occ)[:-1] > 0
    T, RH, C = r.T_a[:-1], 100 * r.RH[:-1], r.CO2[:-1]
    dev = T - T_set
    outside = np.maximum(np.abs(dev) - band, 0.0)
    runtime_h = on.sum() * dt / HOUR
    if occ.any():
        T_rms = float(np.sqrt(np.mean(dev[occ] ** 2)))
        T_max, RH_mean, C_mean = float(T[occ].max()), float(RH[occ].mean()), float(C[occ].mean())
    else:
        warnings.warn("the result has no occupied rows (n_occ > 0); the occupied means, "
                      "maximum and RMS are NaN", RuntimeWarning, stacklevel=2)
        T_rms = T_max = RH_mean = C_mean = float("nan")
    return Metrics(
        energy_kWh=float((unit.E_el[-1] - unit.E_el[0]) / 3.6e6),
        starts=starts,
        starts_per_h_on=starts / runtime_h if runtime_h > 0 else 0.0,
        T_rms_occ=T_rms,
        Kh_outside_occ=float(outside[occ].sum() * dt / HOUR),
        T_max_occ=T_max,
        h_RH60_occ=float((RH[occ] > 60.0).sum() * dt / HOUR),
        RH_mean_occ=RH_mean,
        CO2_mean_occ=C_mean,
        runtime_h=float(runtime_h),
    )


def fan_travel_per_h(r: SimResult) -> float:
    """Return the indoor-fan travel per hour of compressor operation.

    The travel is the total variation of the fan command f while the compressor runs,
    ``sum |f[k+1] - f[k]|`` over the consecutive rows k, k+1 with s > 0 in both. The step
    between the last command of one compressor cycle and the first of the next is not
    counted (the fan is off in between). Row k covers [t_k, t_k + dt), so the last row is
    dropped, as in `compute`. The travel is divided by the running time (rows with s > 0
    times dt). It indicates actuator wear and noise.

    Parameters
    ----------
    r : SimResult
        Result of a run with an air-conditioning unit.

    Returns
    -------
    float
        Fan travel per running hour, 1/h (fan command is a fraction of full speed); NaN when
        the compressor never runs.

    Raises
    ------
    ValueError
        If `r` has no unit record.
    """
    cmd = r.require_unit().cmd[:-1]
    run = cmd[:, IDX_CMD_S] > 0
    runtime_h = run.sum() * r.dt / HOUR
    if runtime_h == 0:
        return float("nan")
    steps = np.abs(np.diff(cmd[:, IDX_CMD_FAN]))
    return float(steps[run[1:] & run[:-1]].sum() / runtime_h)


def low_load_ripples(r: SimResult, s_min: float) -> FloatArray:
    """Return the air-temperature ripple of each complete low-load cycle.

    A cycle runs from one compressor start to the next; it counts as low-load when the
    compressor speed is only `s_min` or zero over the cycle.

    Parameters
    ----------
    r : SimResult
        Result of a run with an air-conditioning unit.
    s_min : float
        Minimum compressor speed, -.

    Returns
    -------
    FloatArray
        Peak-to-peak air temperature over each low-load cycle, K.

    Raises
    ------
    ValueError
        If `r` has no unit record.
    """
    s = r.require_unit().cmd[:-1, IDX_CMD_S]
    on = s > 0
    starts = np.flatnonzero(on[1:] & ~on[:-1]) + 1
    ripples = []
    for a, b in pairwise(starts):
        seg = s[a:b]
        if np.all((seg == 0) | (np.abs(seg - s_min) < 1e-12)):
            ripples.append(np.ptp(r.T_a[a:b + 1]))
    return np.array(ripples)


def energy_cost(r: SimResult, tariff: Callable[[FloatArray], FloatArray]) -> float:
    """Return the electricity cost over the result.

    Power at each grid point is taken as constant over the following step (left rectangle;
    the 60 s step makes the error negligible).

    Parameters
    ----------
    r : SimResult
        Result of a run with an air-conditioning unit.
    tariff : callable
        Price as a function of time: maps times (s) to prices, EUR/kWh.

    Returns
    -------
    float
        Cost of the electrical energy of the unit, EUR.

    Raises
    ------
    ValueError
        If `r` has no unit record.
    """
    P = r.require_unit().op["P_el"][:-1]
    return float(np.sum(P * tariff(r.t[:-1])) * r.dt / 3.6e6)
