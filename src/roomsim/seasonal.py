"""Seasonal efficiency (SEER, SCOP) by the bin method of Regulation (EU) No 206/2012.

The heating part (SCOP, average season, D55) is at the end of the module; it follows section
3.2 of the same communication and reuses the result type and the 10 % tolerance.

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
from .heat_pump import HeatPump

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



# ---- heating (D55) ------------------------------------------------------------------------
# Regulation (EU) No 206/2012, Annex II, Table 1: average heating season, bins -10 to 15 degC
HEATING_BIN_T: FloatArray = np.arange(-10.0, 16.0, dtype=np.float64)
HEATING_BIN_HOURS = np.array([1, 25, 23, 24, 27, 68, 91, 89, 165, 173, 240, 280, 320, 357, 356,
                              303, 330, 326, 348, 335, 315, 215, 169, 151, 105, 74],
                             dtype=float)                         # 4 910 h
T_DESIGN_H = -10.0       # reference design temperature, average season, degC
HEATING_ANCHOR_T = (-7.0, 2.0, 7.0, 12.0)   # declared points A-D, degC
# Table 4, heating mode (average season) of a unit that also cools, h per year
H_HE, H_TO_H, H_SB_H, H_OFF_H, H_CK_H = 1400.0, 179.0, 0.0, 0.0, 179.0
CDH_DEFAULT = 0.25
T_IN_HEATING, T_WB_HEATING = 20.0, 15.0     # indoor rating state, degC


def heating_part_load(T_j: float | FloatArray, P_design: float) -> FloatArray:
    """Return the heating part load Ph(Tj) = Pdesignh (Tj - 16) / (Tdesignh - 16), kW.

    Parameters
    ----------
    T_j : float or ndarray
        Bin temperature, degC.
    P_design : float
        Design load for heating, kW.

    Returns
    -------
    ndarray
        Part load, kW.
    """
    ph: FloatArray = P_design * (np.asarray(T_j, dtype=float) - 16.0) / (T_DESIGN_H - 16.0)
    return ph


def cop_bin_anchor(point: DeclaredPoint, P_design: float, Cdh: float = CDH_DEFAULT) -> float:
    """Return COPbin at a declared heating point (variable-capacity unit).

    Parameters
    ----------
    point : DeclaredPoint
        Declared capacity Pdh (field ``Pdc``, kW) and COPd (field ``EERd``) at the anchor
        temperature.
    P_design : float
        Design load for heating, kW.
    Cdh : float, optional
        Degradation coefficient for heating, -.

    Returns
    -------
    float
        COPd if Pdh matches the part load within the tolerance; COPd degraded by cycling if Pdh
        is larger; COPd if Pdh is smaller (a back-up heater covers the rest).
    """
    Ph = float(heating_part_load(point.T_o, P_design))
    if abs(point.Pdc - Ph) <= TOLERANCE * Ph:
        return point.EERd
    if point.Pdc > Ph:
        return point.EERd * (1.0 - Cdh * (1.0 - Ph / point.Pdc))
    return point.EERd


def scop(points: Sequence[DeclaredPoint], P_design: float, P_dh_design: float,
         modes: LowPowerModes = LowPowerModes(),  # noqa: B008  frozen dataclass
         Cdh: float = CDH_DEFAULT,
         bin_hours: FloatArray = HEATING_BIN_HOURS) -> SeasonalResult:
    """Compute SCOP and SCOPon for the average heating season, bivalent point -7 degC.

    Parameters
    ----------
    points : sequence of DeclaredPoint
        Declared points at -7, 2, 7 and 12 degC (capacity kW, COP in ``EERd``).
    P_design : float
        Design load for heating at -10 degC, kW.
    P_dh_design : float
        Declared capacity at -10 degC, kW (a back-up heater covers P_design - P_dh_design).
    modes : LowPowerModes, optional
        Power in the low-power modes, kW.
    Cdh : float, optional
        Degradation coefficient for heating, -.
    bin_hours : ndarray, shape (26,), optional
        Hours of the bins -10 to 15 degC; the average season by default.

    Returns
    -------
    SeasonalResult
        SCOP (field ``SEER``), SCOPon (``SEER_on``), the annual demand and consumption
        (``Q_C``, ``Q_CE``) and COPbin at -7, 2, 7 and 12 degC (``EER_bin``).

    Raises
    ------
    ValueError
        If the points are not exactly one at each anchor temperature, or bin_hours has the
        wrong shape or a negative entry.
    """
    by_T = {p.T_o: p for p in points}
    if sorted(by_T) != sorted(HEATING_ANCHOR_T) or len(points) != len(HEATING_ANCHOR_T):
        raise ValueError(f"declared points at {[p.T_o for p in points]} degC; exactly one at each "
                         f"of {list(HEATING_ANCHOR_T)} is required")
    hours = np.asarray(bin_hours, dtype=float)
    if hours.shape != HEATING_BIN_T.shape or np.any(hours < 0):
        raise ValueError(f"bin_hours must have {HEATING_BIN_T.size} non-negative entries for "
                         "-10 to 15 degC")
    anchors = np.array([cop_bin_anchor(by_T[T], P_design, Cdh) for T in HEATING_ANCHOR_T])
    # constant below Tbiv = -7 degC, linear between the anchors, extrapolated above 12 degC
    # from the points at 7 and 12 degC
    cop = np.interp(HEATING_BIN_T, HEATING_ANCHOR_T, anchors)
    above = HEATING_BIN_T > 12.0
    cop[above] = anchors[2] + (anchors[3] - anchors[2]) * (HEATING_BIN_T[above] - 7.0) / 5.0
    cap_T = (T_DESIGN_H, *HEATING_ANCHOR_T)
    cap = (P_dh_design, *(by_T[T].Pdc for T in HEATING_ANCHOR_T))
    Pdh = np.interp(HEATING_BIN_T, cap_T, cap)
    Ph = heating_part_load(HEATING_BIN_T, P_design)
    elbu = np.maximum(Ph - Pdh, 0.0)
    scop_on = float(np.sum(hours * Ph) / np.sum(hours * ((Ph - elbu) / cop + elbu)))
    Q_H = P_design * H_HE
    Q_HE = (Q_H / scop_on + H_TO_H * modes.P_TO + H_CK_H * modes.P_CK + H_OFF_H * modes.P_OFF
            + H_SB_H * modes.P_SB)
    return SeasonalResult(SEER=Q_H / Q_HE, SEER_on=scop_on, Q_C=Q_H, Q_CE=Q_HE, EER_bin=anchors)


def heat_pump_declared_points(hp: HeatPump) -> tuple[float, list[DeclaredPoint]]:
    """Compute the design load and the declared heating points of the heat-pump model (D55).

    The unit is evaluated at 20 degC indoors with full indoor fan, a clean outdoor coil and dry
    outdoor air. The design load is its full-speed (s_max) output at -10 degC. At each anchor
    temperature the compressor speed is chosen to deliver the part load, or s_min if the part
    load is below the minimum output.

    Parameters
    ----------
    hp : HeatPump
        Heat-pump model.

    Returns
    -------
    P_design : float
        Design load for heating, kW (also the declared capacity at -10 degC).
    points : list of DeclaredPoint
        Declared capacity (kW) and COP (in ``EERd``) at -7, 2, 7 and 12 degC.
    """
    W_in = W_at_wet_bulb(T_IN_HEATING, T_WB_HEATING)

    def operating(s: float, T_o: float) -> tuple[float, float]:
        st = hp.operate(s, 1.0, T_IN_HEATING, W_in, T_o)
        return -st.Q_sens, st.P_el                      # heat delivered, electrical power; W

    P_design = operating(hp.hp.s_max, T_DESIGN_H)[0] / 1e3
    s_min, s_max = hp.p.s_min, hp.hp.s_max
    points = []
    for T_o in HEATING_ANCHOR_T:
        Ph = float(heating_part_load(T_o, P_design)) * 1e3
        if Ph <= operating(s_min, T_o)[0]:
            s = s_min
        else:
            s = brentq(lambda x, T_o=T_o, Ph=Ph: operating(x, T_o)[0] - Ph, s_min, s_max,
                       xtol=1e-12)
        Q, P = operating(s, T_o)
        points.append(DeclaredPoint(T_o, Q / 1e3, Q / P))
    return P_design, points
