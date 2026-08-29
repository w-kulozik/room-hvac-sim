"""Seasonal energy efficiency ratio (SEER) by the bin method of Regulation (EU) No 206/2012.

Decision D46. The bin hours, the operating hours of the low-power modes and the reference design
temperature are those of the Regulation, Annex II, Tables 1, 3 and 4. The calculation follows
the Commission communication 2012/C 172/01, section 3.1, for variable-capacity units:

    pl(Tj)    = (Tj - 16) / (Tdesignc - 16)                              part load ratio
    Pc(Tj)    = Pdesignc pl(Tj)                                          part load, kW
    EERbin    = EERd                         if |Pdc - Pc| < 10 % of Pc  (eq. 17)
              = EERd (1 - Cdc (1 - Pc/Pdc))  if Pc < Pdc (unit at its lowest capacity, eq. 15)
    SEERon    = sum hj Pc(Tj) / sum hj Pc(Tj) / EERbin(Tj)              (eq. 4)
    SEER      = QC / QCE,  QC = Pdesignc H_CE,
                QCE = QC / SEERon + H_TO P_TO + H_SB P_SB + H_OFF P_OFF + H_CK P_CK

with anchor points at 35, 30, 25 and 20 degC, linear interpolation between them and constant
values beyond them.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from ._typing import FloatArray
from .ac_unit import CP_A, CP_V, H_FG0, ACUnit, _h, _W_sat, net_cooling

# Regulation (EU) No 206/2012, Annex II, Table 1: reference cooling season
BIN_T: FloatArray = np.arange(17.0, 41.0, dtype=np.float64)       # bin temperature Tj, degC
BIN_HOURS = np.array([205, 227, 225, 225, 216, 215, 218, 197, 178, 158, 137, 109, 88, 63, 39,
                      31, 24, 17, 13, 9, 4, 3, 1, 0], dtype=float)  # hj, h per year (2 602 h)
T_DESIGN_C = 35.0        # reference design temperature for cooling, Table 3, degC
ANCHOR_T = (35.0, 30.0, 25.0, 20.0)   # outdoor temperatures of the declared points A-D, degC
# Table 4, cooling mode of a unit that also heats: equivalent active, thermostat-off, standby,
# off and crankcase-heater hours, h per year
H_CE, H_TO, H_SB, H_OFF, H_CK = 350.0, 221.0, 2142.0, 0.0, 2672.0
CDC_DEFAULT = 0.25       # default degradation coefficient for cooling, -
TOLERANCE = 0.10         # capacity tolerance of the variable-capacity rule, -
T_IN_RATING, T_WB_RATING = 27.0, 19.0   # indoor rating state, dry bulb and wet bulb, degC


@dataclass(frozen=True)
class DeclaredPoint:
    """Declared capacity and efficiency at one outdoor temperature.

    Attributes
    ----------
    T_o : float
        Outdoor temperature, degC.
    Pdc : float
        Declared cooling capacity, kW.
    EERd : float
        Declared energy efficiency ratio, W/W.
    """

    T_o: float
    Pdc: float
    EERd: float


@dataclass(frozen=True)
class LowPowerModes:
    """Power drawn in the modes other than active mode, kW.

    Attributes
    ----------
    P_TO : float
        Thermostat-off mode, kW.
    P_SB : float
        Standby mode, kW.
    P_OFF : float
        Off mode, kW.
    P_CK : float
        Crankcase-heater mode, kW.
    """

    P_TO: float = 0.0
    P_SB: float = 0.0
    P_OFF: float = 0.0
    P_CK: float = 0.0


@dataclass(frozen=True)
class SeasonalResult:
    """Result of the bin calculation.

    Attributes
    ----------
    SEER : float
        Seasonal energy efficiency ratio, W/W.
    SEER_on : float
        Active-mode seasonal energy efficiency ratio, W/W.
    Q_C : float
        Reference annual cooling demand, kWh.
    Q_CE : float
        Annual electricity consumption for cooling, kWh.
    EER_bin : ndarray
        EERbin at the anchor temperatures 35, 30, 25 and 20 degC, W/W.
    """

    SEER: float
    SEER_on: float
    Q_C: float
    Q_CE: float
    EER_bin: FloatArray


def part_load(T_j: float | FloatArray, P_design: float) -> FloatArray:
    """Return the part load Pc(Tj) = Pdesignc (Tj - 16) / (Tdesignc - 16).

    Parameters
    ----------
    T_j : float or ndarray
        Bin temperature, degC.
    P_design : float
        Design load for cooling, kW.

    Returns
    -------
    ndarray
        Part load, kW.
    """
    pl: FloatArray = P_design * (np.asarray(T_j, dtype=float) - 16.0) / (T_DESIGN_C - 16.0)
    return pl


def eer_bin_anchor(point: DeclaredPoint, P_design: float, Cdc: float = CDC_DEFAULT) -> float:
    """Return EERbin at a declared point (variable-capacity unit).

    Parameters
    ----------
    point : DeclaredPoint
        Declared capacity and EER at the anchor temperature.
    P_design : float
        Design load for cooling, kW.
    Cdc : float, optional
        Degradation coefficient for cooling, -.

    Returns
    -------
    float
        EERbin, W/W: EERd if the declared capacity matches the part load within the tolerance,
        EERd degraded by cycling if the declared capacity is larger.

    Raises
    ------
    ValueError
        If the declared capacity is below the part load by more than the tolerance (the rule
        would need a second, higher declared capacity).
    """
    Pc = float(part_load(point.T_o, P_design))
    if abs(point.Pdc - Pc) < TOLERANCE * Pc:
        return point.EERd
    if point.Pdc > Pc:
        return point.EERd * (1.0 - Cdc * (1.0 - Pc / point.Pdc))
    raise ValueError(f"declared capacity {point.Pdc!r} kW at {point.T_o!r} degC is below the "
                     f"part load {Pc:.3f} kW; the unit cannot meet the load")


def seer(points: Sequence[DeclaredPoint], P_design: float,
         modes: LowPowerModes = LowPowerModes(),  # noqa: B008  frozen dataclass
         Cdc: float = CDC_DEFAULT, bin_hours: FloatArray = BIN_HOURS) -> SeasonalResult:
    """Compute SEER and SEERon by the bin method.

    Parameters
    ----------
    points : sequence of DeclaredPoint
        Declared points at 35, 30, 25 and 20 degC (in any order).
    P_design : float
        Design load for cooling, kW (the declared capacity at 35 degC).
    modes : LowPowerModes, optional
        Power in the low-power modes, kW.
    Cdc : float, optional
        Degradation coefficient for cooling, -.
    bin_hours : ndarray, shape (24,), optional
        Hours of the bins 17-40 degC; the reference cooling season by default.

    Returns
    -------
    SeasonalResult
        SEER, SEERon, the annual demand and consumption, and the anchor EERbin values.

    Raises
    ------
    ValueError
        If the points are not exactly one at each of 35, 30, 25 and 20 degC, or bin_hours has
        the wrong shape or a negative entry.
    """
    by_T = {p.T_o: p for p in points}
    if sorted(by_T) != sorted(ANCHOR_T) or len(points) != len(ANCHOR_T):
        raise ValueError(f"declared points at {[p.T_o for p in points]} degC; exactly one at each "
                         f"of {list(ANCHOR_T)} is required")
    hours = np.asarray(bin_hours, dtype=float)
    if hours.shape != BIN_T.shape or np.any(hours < 0):
        raise ValueError(f"bin_hours must have {BIN_T.size} non-negative entries for 17-40 degC")
    anchors = np.array([eer_bin_anchor(by_T[T], P_design, Cdc) for T in ANCHOR_T])
    # np.interp needs increasing abscissae and holds the end values beyond them
    eer = np.interp(BIN_T, ANCHOR_T[::-1], anchors[::-1])
    Pc = part_load(BIN_T, P_design)
    seer_on = float(np.sum(hours * Pc) / np.sum(hours * Pc / eer))
    Q_C = P_design * H_CE
    Q_CE = (Q_C / seer_on + H_TO * modes.P_TO + H_SB * modes.P_SB + H_OFF * modes.P_OFF
            + H_CK * modes.P_CK)
    return SeasonalResult(SEER=Q_C / Q_CE, SEER_on=seer_on, Q_C=Q_C, Q_CE=Q_CE, EER_bin=anchors)


def local_bin_hours(T_hourly: FloatArray) -> FloatArray:
    """Count the hours of a weather record in the cooling bins 17-40 degC.

    Parameters
    ----------
    T_hourly : ndarray
        Hourly outdoor temperatures, degC.

    Returns
    -------
    ndarray, shape (24,)
        Hours per bin; an hour belongs to the bin of its temperature rounded to the nearest
        degree (halves upwards). Hours outside 16.5-40.5 degC are not counted.
    """
    T_bin = np.floor(np.asarray(T_hourly, dtype=float) + 0.5)
    counts: FloatArray = np.array([np.sum(T_bin == Tj) for Tj in BIN_T], dtype=float)
    return counts


def W_at_wet_bulb(T: float, T_wb: float) -> float:
    """Return the humidity ratio on the constant-enthalpy line of a wet-bulb temperature.

    Parameters
    ----------
    T : float
        Dry-bulb temperature, degC.
    T_wb : float
        Wet-bulb temperature, degC.

    Returns
    -------
    float
        Humidity ratio of air at dry bulb `T` with the enthalpy of saturated air at `T_wb`,
        kg/kg.
    """
    return (_h(T_wb, _W_sat(T_wb)) - CP_A * T) / (H_FG0 + CP_V * T)


def unit_declared_points(unit: ACUnit) -> tuple[float, list[DeclaredPoint]]:
    """Compute the design load and the declared points of the unit model.

    The unit is evaluated at the indoor rating state (27 degC, 19 degC wet bulb) with the indoor
    fan at full speed. The design load is its net capacity at 35 degC and full speed. At each
    anchor temperature the compressor speed is chosen to deliver the part load, or the minimum
    speed if the part load is below the minimum capacity.

    Parameters
    ----------
    unit : ACUnit
        Unit model.

    Returns
    -------
    P_design : float
        Design load for cooling, kW.
    points : list of DeclaredPoint
        Declared capacity (kW) and net EER at 35, 30, 25 and 20 degC.

    Raises
    ------
    ValueError
        If the part load at an anchor temperature exceeds the unit's full capacity there.
    """
    W_in = W_at_wet_bulb(T_IN_RATING, T_WB_RATING)

    def net(s: float, T_o: float) -> float:
        # net cooling at the rating indoor state, full fan, compressor speed s; W
        return net_cooling(unit.operate(s, 1.0, T_IN_RATING, W_in, T_o), T_IN_RATING)

    def eer(s: float, T_o: float) -> float:
        return net(s, T_o) / unit.operate(s, 1.0, T_IN_RATING, W_in, T_o).P_el

    P_design = net(1.0, T_DESIGN_C) / 1e3
    s_min = unit.p.s_min
    points = []
    for T_o in ANCHOR_T:
        Pc = float(part_load(T_o, P_design)) * 1e3
        if T_o == T_DESIGN_C:
            s = 1.0
        elif Pc <= net(s_min, T_o):
            s = s_min
        elif Pc <= net(1.0, T_o):
            s = brentq(lambda x, T_o=T_o, Pc=Pc: net(x, T_o) - Pc, s_min, 1.0, xtol=1e-12)
        else:
            raise ValueError(f"part load {Pc:.0f} W at {T_o} degC exceeds the full capacity")
        points.append(DeclaredPoint(T_o, net(s, T_o) / 1e3, eer(s, T_o)))
    return P_design, points
