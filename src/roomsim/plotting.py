"""Shared matplotlib styling and the baseline figures.

Colours: first three slots of a categorical palette validated for colour-vision deficiency
on all pairs. Identity is never colour-only: every series also has a legend entry and,
where there is room, a direct label.
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING, cast

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ._typing import BoolArray, FloatArray
from .inputs import DAY, HOUR
from .results import SimResult

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from matplotlib.axes import Axes
    from numpy.typing import ArrayLike

    from .controllers import ComfortSchedule
    from .metrics import Metrics
    from .tariff import TimeOfUseTariff

StrPath = str | os.PathLike[str]    # output file accepted by Figure.savefig

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK_2, MUTED, GRID, SHADE = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1", "#f0efec"
REF = "#8a8984"   # reference lines (analytical slopes)


def apply_style() -> None:
    """Set the shared matplotlib rcParams: colours, fonts, grid, legend and output resolution."""
    plt.rcParams.update({
        "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb",
        "axes.edgecolor": MUTED, "axes.labelcolor": INK_2, "axes.titlecolor": INK,
        "axes.titlesize": 11, "axes.titleweight": "semibold", "axes.titlelocation": "left",
        "axes.labelsize": 9.5, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "xtick.color": INK_2, "ytick.color": INK_2, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "legend.frameon": False, "legend.fontsize": 8.5, "legend.labelcolor": INK_2,
        "lines.linewidth": 1.8, "font.family": "DejaVu Sans", "savefig.dpi": 150,
        "savefig.bbox": "tight",
    })


def _time_axis(ax: Axes, r: SimResult) -> FloatArray:
    """Set 12-hourly clock ticks on `ax`; return the time since the start of `r`, h."""
    th = cast(FloatArray, (r.t - r.t[0]) / HOUR)
    ax.set_xlim(th[0], th[-1])
    ax.set_xticks(np.arange(0, th[-1] + 1, 12))
    ax.set_xticklabels([f"{int(h % 24):02d}:00\nday {int(h // 24) + 1}" if h % 24 == 0
                        else f"{int(h % 24):02d}:00" for h in np.arange(0, th[-1] + 1, 12)])
    return th


def _shade_occupancy(ax: Axes, r: SimResult, th: FloatArray, label: bool = True) -> None:
    """Shade the occupied intervals of `r` on `ax`; `th` is the time axis, h."""
    occ = np.asarray(r.d.n_occ) > 0
    edges = np.flatnonzero(np.diff(occ.astype(int)))
    starts = list(th[edges[occ[edges + 1]] + 1]) if len(edges) else []
    ends = list(th[edges[~occ[edges + 1]] + 1]) if len(edges) else []
    if occ[0]:
        starts.insert(0, th[0])
    if occ[-1]:
        ends.append(th[-1])
    # rising and falling edges alternate, so every start has its end
    for i, (a, b) in enumerate(zip(starts, ends, strict=True)):
        ax.axvspan(a, b, color=SHADE, zorder=0, lw=0,
                   label="occupied" if (i == 0 and label) else None)


def _held(y: FloatArray) -> FloatArray:
    """Return `y` with its last value replaced by the previous one, for post step plots.

    The command of the last row acts after the end of the result; repeating the previous
    command ends the step line at t_N without drawing the stale value.
    """
    return np.append(y[:-1], y[-2]) if len(y) > 1 else y


def _peak_spans(peak: BoolArray) -> list[tuple[int, int]]:
    """Return the runs of True in `peak` as (first index, index after the run) pairs.

    A run reaching the end of the array ends at the last index.
    """
    d = np.diff(np.concatenate(([0], peak.astype(int), [0])))
    last = len(peak) - 1
    return [(int(a), min(int(b), last))
            for a, b in zip(np.flatnonzero(d > 0), np.flatnonzero(d < 0), strict=True)]


def _legend_top(ax: Axes, ncol: int) -> None:
    """Place the legend of `ax` above the plot area in `ncol` columns and raise the title."""
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=ncol, borderaxespad=0.2,
              handlelength=1.8, columnspacing=1.4)
    ax.set_title(ax.get_title(loc="left"), loc="left", pad=24)


def _direct_label(ax: Axes, x: float, y: float, text: str, color: str) -> None:
    """Write `text` next to the point (x, y) of `ax`, in a box outlined in `color`."""
    ax.annotate(text, (x, y), xytext=(4, 0), textcoords="offset points",
                va="center", fontsize=8.5, color=INK_2,
                bbox={"boxstyle": "round,pad=0.15", "fc": "#fcfcfb", "ec": color, "lw": 1.2})


def plot_temperatures(r: SimResult, path: StrPath) -> None:
    """Plot outdoor, thermal-mass and indoor air temperatures and save the figure.

    The maximum air temperature is marked and labelled; occupied hours are shaded.

    Parameters
    ----------
    r : SimResult
        Result to plot.
    path : str or os.PathLike
        Output file.
    """
    fig, ax = plt.subplots(figsize=(9, 3.8))
    th = _time_axis(ax, r)
    _shade_occupancy(ax, r, th)
    ax.plot(th, r.d.T_o, color=ORANGE, label="outdoor $T_o$")
    ax.plot(th, r.T_m, color=AQUA, label="thermal mass $T_m$")
    ax.plot(th, r.T_a, color=BLUE, label="indoor air $T_a$")
    k = int(np.argmax(r.T_a))
    ax.plot(th[k], r.T_a[k], "o", ms=6, color=BLUE, mec="#fcfcfb", mew=1.5)
    ax.annotate(f"max {r.T_a[k]:.1f} °C at {(r.t[k] % DAY) / HOUR:.1f} h", (th[k], r.T_a[k]),
                xytext=(-10, -6), textcoords="offset points", ha="right", va="top",
                fontsize=8.5, color=INK,
                bbox={"boxstyle": "round,pad=0.2", "fc": "#fcfcfb", "ec": "none", "alpha": 0.9})
    ax.set_ylabel("temperature (°C)")
    ax.set_title("Free-floating room temperatures, summer periodic steady state")
    _legend_top(ax, 4)
    fig.savefig(path)
    plt.close(fig)


def plot_humidity(r: SimResult, path: StrPath) -> None:
    """Plot relative humidity (top) and dew point (bottom), outdoor and indoor, and save.

    Parameters
    ----------
    r : SimResult
        Result to plot.
    path : str or os.PathLike
        Output file.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 5.6), sharex=True)
    th = _time_axis(ax2, r)
    for ax in (ax1, ax2):
        _shade_occupancy(ax, r, th, label=ax is ax1)
    from .psychro import RH_from_W
    ax1.plot(th, 100 * RH_from_W(r.d.T_o, r.d.W_o), color=ORANGE, label="outdoor")
    ax1.plot(th, 100 * r.RH, color=BLUE, label="indoor air")
    ax1.plot(th, 100 * r.phi_b, color=AQUA, lw=1.4, ls=(0, (4, 2)),
             label="wall buffer $\\varphi_b$")
    ax1.set_ylabel("relative humidity (%)")
    ax1.set_title("Relative humidity")
    _legend_top(ax1, 4)
    ax2.plot(th, r.T_dew_out, color=ORANGE, label="outdoor")
    ax2.plot(th, r.T_dew, color=BLUE, label="indoor air")
    ax2.set_ylabel("dew point (°C)")
    ax2.set_title("Dew point")
    _legend_top(ax2, 2)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_energy_balance(r: SimResult, path: StrPath) -> None:
    """Plot boundary heat against stored energy (top) and their residual (bottom), and save.

    Parameters
    ----------
    r : SimResult
        Result to plot.
    path : str or os.PathLike
        Output file.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 5.6), sharex=True)
    th = _time_axis(ax2, r)
    ax1.plot(th, r.E_in / 1e6, color=ORANGE, lw=3.2, label="∫ boundary heat flows dt")
    ax1.plot(th, r.energy_stored / 1e6, color=BLUE, lw=1.4,
             label="$C_a\\,\\Delta T_a + C_m\\,\\Delta T_m$")
    ax1.set_ylabel("energy since start (MJ)")
    ax1.set_title("Energy balance: heat in vs heat stored")
    _legend_top(ax1, 2)
    ax2.plot(th, r.energy_residual, color=INK_2, lw=1.2)
    ax2.axhline(0, color=MUTED, lw=0.8)
    lim = max(np.abs(r.energy_residual).max() * 1.3, 1e-9)
    ax2.set_ylim(-lim, lim)
    ax2.set_ylabel("residual (J)")
    E_scale = np.abs(r.energy_stored).max()
    ax2.set_title(f"Residual (max {np.abs(r.energy_residual).max():.1e} J, "
                  f"{np.abs(r.energy_residual).max() / E_scale:.0e} of the stored swing)")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_convergence(
    studies: Mapping[str, tuple[FloatArray, Mapping[str, ArrayLike], str]], path: StrPath
) -> None:
    """Plot error against time step for RK4 and Euler, one panel per study, and save.

    Each panel gives the fitted slope in the legend and a reference line of the theoretical
    order (4 for RK4, 1 for Euler) through the error at the smallest step.

    Parameters
    ----------
    studies : mapping
        ``{title: (dts, {method: errors}, ylabel)}``: time steps `dts` (s), one error per step
        for each of the methods ``"rk4"`` and ``"euler"`` (unit given by `ylabel`), and the
        y-axis label.
    path : str or os.PathLike
        Output file.
    """
    fig, axes = plt.subplots(1, len(studies), figsize=(9, 3.8))
    for ax, (title, (dts, errs, ylabel)) in zip(np.atleast_1d(axes), studies.items(),
                                                strict=True):
        from .verification import convergence_order
        for method, color, marker, order in (("rk4", BLUE, "o", 4), ("euler", ORANGE, "s", 1)):
            e = np.asarray(errs[method])
            slope = convergence_order(dts, e)
            ax.loglog(dts, e, marker=marker, ms=6, color=color, mec="#fcfcfb", mew=1.2,
                      label=f"{method.upper() if method == 'rk4' else 'Euler'}: slope {slope:.2f}")
            ref = e[-1] * (np.asarray(dts) / dts[-1]) ** order
            ax.loglog(dts, ref, color=REF, lw=1.0, ls=(0, (3, 3)), zorder=0)
            _direct_label(ax, dts[0], ref[0], f"$\\Delta t^{order}$", REF)
        ax.set_xticks(dts)
        ax.set_xticklabels([f"{d:g}" for d in dts])
        ax.minorticks_off()
        ax.set_xlabel("time step Δt (s)")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# ---- controller comparison ----------------------------------------------------------------

YELLOW, MAGENTA = "#eda100", "#e87ba4"
SERIES = (BLUE, ORANGE, AQUA, YELLOW, MAGENTA)    # fixed order: colour follows the controller


def plot_controller_day(results: Mapping[str, SimResult], T_set: float, path: StrPath,
                        band: float = 1.0) -> None:
    """Plot one day per controller: air temperature (top) and compressor speed (bottom).

    Parameters
    ----------
    results : mapping of str to SimResult
        One-day results with a unit record, keyed by controller name (at most five; colours
        follow `SERIES`).
    T_set : float
        Temperature set point, degC.
    path : str or os.PathLike
        Output file.
    band : float, optional
        Half-width of the shaded band around `T_set`, K.

    Raises
    ------
    ValueError
        If a result has no unit record.
    """
    n = len(results)
    fig = plt.figure(figsize=(9, 2.25 * n + 0.4))
    subfigs = fig.subfigures(n, 1, hspace=0.02)
    for i, (sf, (name, r), color) in enumerate(zip(subfigs, results.items(), SERIES,
                                                   strict=False)):
        axT, axS = sf.subplots(2, 1, sharex=True,
                               gridspec_kw={"height_ratios": [3, 1], "hspace": 0.12})
        th = (r.t - r.t[0]) / HOUR
        for ax in (axT, axS):
            _shade_occupancy(ax, r, th, label=False)
            ax.set_axisbelow(True)
        axT.axhspan(T_set - band, T_set + band, color=GRID, alpha=0.5, lw=0, zorder=0)
        axT.axhline(T_set, color=MUTED, lw=0.8, ls=(0, (3, 3)))
        axT.plot(th, r.T_a, color=color, lw=1.4)
        axT.set_ylim(T_set - band - 0.2, T_set + band + 0.2)
        axT.set_ylabel("$T_a$ (°C)")
        axT.set_title(name, loc="left", fontsize=10, pad=4)
        axS.step(th, _held(r.require_unit().cmd[:, 0]), where="post", color=color, lw=1.1)
        axS.set_ylim(-0.05, 1.08)
        axS.set_yticks([0, 0.3, 1])
        axS.set_ylabel("s (-)")
        axS.set_xlim(0, 24)
        axS.set_xticks(np.arange(0, 25, 3))
        axS.set_xticklabels([f"{h:02d}:00" for h in range(0, 25, 3)] if i == n - 1 else [])
    fig.suptitle("Air temperature and compressor speed over one day "
                 "(grey: occupied; band: set point ± 1 K)",
                 x=0.06, y=1.0, ha="left", va="bottom", fontsize=11, fontweight="semibold")
    fig.savefig(path)
    plt.close(fig)


def plot_controller_metrics(metrics: Mapping[str, Metrics], path: StrPath) -> None:
    """Plot small multiples of the main indicators, one bar per controller, and save.

    Panels: electrical energy (kWh), compressor starts, RMS temperature deviation while
    occupied (K) and hours above 60 % RH while occupied (h).

    Parameters
    ----------
    metrics : mapping of str to Metrics
        Indicators keyed by controller name (at most five; colours follow `SERIES`).
    path : str or os.PathLike
        Output file.
    """
    panels = [("energy_kWh", "Electrical energy (kWh)"), ("starts", "Compressor starts"),
              ("T_rms_occ", "RMS deviation, occupied (K)"),
              ("h_RH60_occ", "Hours above 60 % RH, occupied")]
    names = list(metrics)
    fig, axes = plt.subplots(1, len(panels), figsize=(11, 3.0), sharey=True)
    y = np.arange(len(names))[::-1]
    for ax, (key, title) in zip(axes, panels, strict=True):
        ax.set_axisbelow(True)
        vals = [getattr(metrics[n], key) for n in names]
        ax.barh(y, vals, color=SERIES[:len(names)], height=0.62, edgecolor="#fcfcfb", linewidth=2)
        vmax = max(vals) if max(vals) > 0 else 1.0
        for yi, v in zip(y, vals, strict=True):
            ax.text(v + 0.02 * vmax, yi, f"{v:.3g}" if isinstance(v, float) else str(v),
                    va="center", fontsize=8.5, color=INK_2)
        ax.set_xlim(0, vmax * 1.25)
        ax.set_title(title, fontsize=9.5)
        ax.grid(axis="y", visible=False)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(names)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def plot_humidity_day(results: Mapping[str, SimResult], RH_set: float, path: StrPath) -> None:
    """Plot one day per controller: indoor RH (top) and indoor fan speed while running (bottom).

    Parameters
    ----------
    results : mapping of str to SimResult
        One-day results with a unit record, keyed by controller name (at most five; colours
        follow `SERIES`).
    RH_set : float
        Relative-humidity set point, fraction.
    path : str or os.PathLike
        Output file.

    Raises
    ------
    ValueError
        If a result has no unit record.
    """
    n = len(results)
    fig = plt.figure(figsize=(9, 2.25 * n + 0.4))
    subfigs = fig.subfigures(n, 1, hspace=0.02)
    for i, (sf, (name, r), color) in enumerate(zip(subfigs, results.items(), SERIES,
                                                   strict=False)):
        axH, axF = sf.subplots(2, 1, sharex=True,
                               gridspec_kw={"height_ratios": [3, 1], "hspace": 0.12})
        th = (r.t - r.t[0]) / HOUR
        for ax in (axH, axF):
            _shade_occupancy(ax, r, th, label=False)
            ax.set_axisbelow(True)
        axH.axhline(100 * RH_set, color=MUTED, lw=0.8, ls=(0, (3, 3)))
        axH.axhline(60.0, color=MUTED, lw=0.8)
        axH.plot(th, 100 * r.RH, color=color, lw=1.4)
        axH.set_ylim(45, 70)
        axH.set_ylabel("RH (%)")
        axH.set_title(name, loc="left", fontsize=10, pad=4)
        cmd = r.require_unit().cmd
        fan = np.where(cmd[:, 0] > 0, cmd[:, 1], np.nan)
        axF.step(th, _held(fan), where="post", color=color, lw=1.1)
        axF.set_ylim(0.3, 1.05)
        axF.set_yticks([0.4, 1.0])
        axF.set_ylabel("fan (-)")
        axF.set_xlim(0, 24)
        axF.set_xticks(np.arange(0, 25, 3))
        axF.set_xticklabels([f"{h:02d}:00" for h in range(0, 25, 3)] if i == n - 1 else [])
    fig.suptitle(f"Indoor RH and fan speed over one day (dashed: set point {100 * RH_set:.0f} %; "
                 "solid: 60 %; grey: occupied)",
                 x=0.06, y=1.0, ha="left", va="bottom", fontsize=11, fontweight="semibold")
    fig.savefig(path)
    plt.close(fig)


def plot_mpc_day(results: Mapping[str, SimResult], comfort: ComfortSchedule,
                 tariff: TimeOfUseTariff, path: StrPath) -> None:
    """Plot air temperature against the comfort band (top) and hourly mean power (bottom).

    Peak-tariff hours are shaded. The grid of the first result is used for the band and the
    tariff.

    Parameters
    ----------
    results : mapping of str to SimResult
        One-day results with a unit record on the same grid, keyed by controller name (at
        most five; colours follow `SERIES`).
    comfort : ComfortSchedule
        Comfort band of the air temperature, degC.
    tariff : TimeOfUseTariff
        Electricity tariff; prices above `tariff.offpeak` mark the peak hours.
    path : str or os.PathLike
        Output file.

    Raises
    ------
    ValueError
        If a result has no unit record.
    """
    fig, (axT, axP) = plt.subplots(2, 1, figsize=(9, 6.0), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 2], "hspace": 0.12})
    r0 = next(iter(results.values()))
    th = (r0.t - r0.t[0]) / HOUR
    lo, hi = comfort.bounds(r0.t)
    axT.fill_between(th, lo, hi, color=GRID, alpha=0.6, lw=0, step="post", label="comfort band")
    peak = np.asarray(tariff(r0.t)) > tariff.offpeak
    for ax in (axT, axP):
        ax.set_axisbelow(True)
        for a, b in _peak_spans(peak):
            ax.axvspan(th[a], th[b], color=SHADE, lw=0, zorder=0)
    for (name, r), color in zip(results.items(), SERIES, strict=False):
        axT.plot(th, r.T_a, color=color, lw=1.3, label=name)
        hours = np.floor(th[:-1]).astype(int)
        P_h = np.bincount(hours, weights=r.require_unit().op["P_el"][:-1]) / np.bincount(hours)
        axP.step(np.arange(len(P_h) + 1), np.append(P_h, P_h[-1]), where="post", color=color,
                 lw=1.5)
    axT.set_ylim(min(19.5, min(r.T_a.min() for r in results.values()) - 0.3), 28.8)
    axT.set_ylabel("$T_a$ (°C)")
    axT.set_title("Air temperature and electrical power over one day "
                  "(grey band: comfort; shaded: peak tariff)", loc="left")
    _legend_top(axT, 4)
    axP.set_ylabel("mean power (W)")
    axP.set_xlim(0, 24)
    axP.set_xticks(np.arange(0, 25, 3))
    axP.set_xticklabels([f"{h:02d}:00" for h in range(0, 25, 3)])
    fig.savefig(path)
    plt.close(fig)


def plot_ventilation_day(results: Mapping[str, SimResult], co2_limit: float,
                         path: StrPath) -> None:
    """Plot one day per strategy: CO2 (top) and ventilation flow (bottom), and save.

    Parameters
    ----------
    results : mapping of str to SimResult
        One-day results keyed by strategy name (at most five; colours follow `SERIES`).
    co2_limit : float
        CO2 concentration drawn as a reference line, ppm.
    path : str or os.PathLike
        Output file.
    """
    from .plant import IDX_V_VENT
    n = len(results)
    fig = plt.figure(figsize=(9, 2.25 * n + 0.4))
    subfigs = fig.subfigures(n, 1, hspace=0.02)
    for i, (sf, (name, r), color) in enumerate(zip(subfigs, results.items(), SERIES,
                                                   strict=False)):
        axC, axV = sf.subplots(2, 1, sharex=True,
                               gridspec_kw={"height_ratios": [3, 1], "hspace": 0.12})
        th = (r.t - r.t[0]) / HOUR
        for ax in (axC, axV):
            _shade_occupancy(ax, r, th, label=False)
            ax.set_axisbelow(True)
        axC.axhline(co2_limit, color=MUTED, lw=0.8)
        axC.plot(th, r.CO2, color=color, lw=1.4)
        axC.set_ylim(400, 1300)
        axC.set_ylabel("CO$_2$ (ppm)")
        axC.set_title(name, loc="left", fontsize=10, pad=4)
        axV.step(th, _held(3600 * r.u[:, IDX_V_VENT]), where="post", color=color, lw=1.1)
        axV.set_ylim(-5, 115)
        axV.set_yticks([0, 50, 100])
        axV.set_ylabel("m$^3$/h")
        axV.set_xlim(0, 24)
        axV.set_xticks(np.arange(0, 25, 3))
        axV.set_xticklabels([f"{h:02d}:00" for h in range(0, 25, 3)] if i == n - 1 else [])
    fig.suptitle(f"CO$_2$ and ventilation flow over one day (line: {co2_limit:.0f} ppm; "
                 "grey: occupied)",
                 x=0.06, y=1.0, ha="left", va="bottom", fontsize=11, fontweight="semibold")
    fig.savefig(path)
    plt.close(fig)


def plot_season_monthly(monthly: Mapping[tuple[str, str], tuple[Sequence[float], Sequence[float]]],
                        sites: Sequence[str], controllers: Sequence[str], months: Sequence[str],
                        path: StrPath,
                        title: str = ("Cooling season by month: electricity of the unit and net "
                                      "cooling per kWh (not shown below 1 kWh)"),
                        ratio_label: str = "EER in use (W/W)") -> None:
    """Plot the monthly electricity and EER in use per site and controller, and save.

    Top row: electrical energy of the unit per month (kWh); bottom row: net cooling over
    electricity per month (W/W), left out for months below 1 kWh, where the ratio is dominated
    by start-up transients and fan-only operation. One column per site, one bar per controller.

    Parameters
    ----------
    monthly : mapping of (site, controller) to (electricity, delivered)
        Monthly electrical energy and net cooling (or heating) delivered, kWh, one value per
        month.
    sites : sequence of str
        Site names, one column each.
    controllers : sequence of str
        Controller names (at most five; colours follow `SERIES`).
    months : sequence of str
        Month labels.
    path : str or os.PathLike
        Output file.
    title : str, optional
        Figure title.
    ratio_label : str, optional
        Label of the bottom row (delivered over electricity).
    """
    n_c = len(controllers)
    width = 0.8 / n_c
    x = np.arange(len(months))
    fig, axes = plt.subplots(2, len(sites), figsize=(4.6 * len(sites), 5.4), sharex=True,
                             sharey="row", squeeze=False)
    for j, site in enumerate(sites):
        axE, axR = axes[0, j], axes[1, j]
        for i, (name, color) in enumerate(zip(controllers, SERIES, strict=False)):
            elec, cool = (np.asarray(v, dtype=float) for v in monthly[(site, name)])
            xi = x + (i - (n_c - 1) / 2) * width
            axE.bar(xi, elec, width, color=color, label=name, edgecolor="#fcfcfb", linewidth=1)
            eer = np.divide(cool, elec, out=np.full_like(cool, np.nan), where=elec >= 1.0)
            axR.bar(xi, eer, width, color=color, edgecolor="#fcfcfb", linewidth=1)
        axE.set_title(site)
        for ax in (axE, axR):
            ax.set_axisbelow(True)
            ax.grid(axis="x", visible=False)
        axR.set_xticks(x)
        axR.set_xticklabels(months)
    axes[0, 0].set_ylabel("Electricity (kWh)")
    axes[1, 0].set_ylabel(ratio_label)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.06, 1.0), ncol=n_c)
    fig.suptitle(title, x=0.06, y=1.06, ha="left", va="bottom", fontsize=11, fontweight="semibold")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
