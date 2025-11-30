"""Performance indicators for comparing controllers on the same scenario."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .results import SimResult
from .simulate import IDX_CMD_S

HOUR = 3600.0


@dataclass(frozen=True)
class Metrics:
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

    def as_dict(self) -> dict:
        return asdict(self)


def compute(r: SimResult, T_set: float, band: float = 1.0) -> Metrics:
    """Indicators over the whole result. Row k covers [t_k, t_k + dt), so the last row is dropped."""
    dt = r.dt
    s = r.cmd[:-1, IDX_CMD_S]
    on = s > 0
    starts = int(np.sum(on[1:] & ~on[:-1]) + (1 if on[0] else 0))
    occ = np.asarray(r.d.n_occ)[:-1] > 0
    T, RH, C = r.T_a[:-1], 100 * r.RH[:-1], r.CO2[:-1]
    dev = T - T_set
    outside = np.maximum(np.abs(dev) - band, 0.0)
    runtime_h = on.sum() * dt / HOUR
    return Metrics(
        energy_kWh=float((r.E_el[-1] - r.E_el[0]) / 3.6e6),
        starts=starts,
        starts_per_h_on=starts / runtime_h if runtime_h > 0 else 0.0,
        T_rms_occ=float(np.sqrt(np.mean(dev[occ] ** 2))),
        Kh_outside_occ=float(outside[occ].sum() * dt / HOUR),
        T_max_occ=float(T[occ].max()),
        h_RH60_occ=float((RH[occ] > 60.0).sum() * dt / HOUR),
        RH_mean_occ=float(RH[occ].mean()),
        CO2_mean_occ=float(C[occ].mean()),
        runtime_h=float(runtime_h),
    )
