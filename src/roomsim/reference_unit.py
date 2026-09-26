"""Published data of the real unit the model is compared with (cooling and heating).

Daikin FTXM25R + RXM25R, a 2.5 kW wall-mounted inverter split (R-32). Values from the
manufacturer's technical data books (full references in docs/references.md):
    RXM-R technical data book EEDEN20, pages 19-21 (combination FTXM25R + RXM25R)
    FTXM-R technical data EEDDE22, printed page 6, page 3 of the PDF file (indoor unit)
Cooling ratings are at 27 degC dry bulb / 19 degC wet bulb indoors and 35 degC outdoors; the
part-load points B-D are the declared points of the seasonal (SEER) rating at 30, 25 and 20 degC.
Heating ratings are at 20 degC indoors and 7 degC dry bulb / 6 degC wet bulb outdoors; the
declared points of the seasonal (SCOP) rating, average season, are at -7, 2, 7 and 12 degC, with
the bivalent point at -7 degC. The data book prints the three sizes in interleaved columns; at
12 degC the capacity of this size is 0.99 or 1.08 kW and 1.08 is used (the SCOP changes by 0.004).
"""

from __future__ import annotations

from dataclasses import dataclass, replace

from .seasonal import DeclaredPoint, LowPowerModes


@dataclass(frozen=True)
class Datasheet:
    """Published values of the reference unit in cooling and heating.

    Attributes
    ----------
    cap_nom_kW, cap_min_kW, cap_max_kW : float
        Nominal, minimum and maximum cooling capacity, kW.
    P_nom_kW : float
        Power input at the nominal capacity, kW.
    EER : float
        Nominal energy efficiency ratio, W/W.
    SEER : float
        Seasonal energy efficiency ratio, W/W.
    P_design_kW : float
        Design load for cooling (Pdesignc), kW.
    airflow_m3h : tuple of int
        Indoor airflow at the high, medium, low and silent fan steps, m3/h.
    P_indoor_W : float
        Indoor unit power input in nominal cooling, W.
    P_thermostat_off_W, P_standby_W, P_off_W : float
        Power in thermostat-off, standby and off mode, W (no crankcase heater power is
        listed).
    points : tuple of (str, float, float, float)
        Declared part-load points: name, outdoor temperature (degC), capacity Pdc (kW), EERd.
    heat_nom_kW, heat_min_kW, heat_max_kW : float
        Nominal, minimum and maximum heating capacity, kW.
    COP : float
        Nominal coefficient of performance, W/W.
    SCOP : float
        Seasonal coefficient of performance, average season, W/W.
    P_design_heating_kW : float
        Design load for heating (Pdesignh) at -10 degC, kW.
    P_dh_design_kW : float
        Declared heating capacity at -10 degC, kW.
    P_thermostat_off_heating_W : float
        Power in thermostat-off mode in heating, W.
    heating_points : tuple of (str, float, float, float)
        Declared heating points: name, outdoor temperature (degC), capacity Pdh (kW), COPd.
    """

    cap_nom_kW: float = 2.50
    cap_min_kW: float = 1.30
    cap_max_kW: float = 3.20
    P_nom_kW: float = 0.56
    EER: float = 4.50
    SEER: float = 8.65
    P_design_kW: float = 2.50
    airflow_m3h: tuple[int, int, int, int] = (630, 456, 342, 246)  # high, medium, low, silent
    P_indoor_W: float = 25.0
    P_thermostat_off_W: float = 6.0
    P_standby_W: float = 1.0
    P_off_W: float = 1.0
    points: tuple[tuple[str, float, float, float], ...] = (
        ("A", 35.0, 2.50, 4.50), ("B", 30.0, 1.85, 6.52), ("C", 25.0, 1.19, 10.17),
        ("D", 20.0, 1.17, 16.51))
    heat_nom_kW: float = 2.80
    heat_min_kW: float = 1.30
    heat_max_kW: float = 4.70
    COP: float = 5.00
    SCOP: float = 5.10
    P_design_heating_kW: float = 2.40
    P_dh_design_kW: float = 2.30
    P_thermostat_off_heating_W: float = 7.0
    heating_points: tuple[tuple[str, float, float, float], ...] = (
        ("A", -7.0, 2.13, 3.60), ("B", 2.0, 1.29, 5.14), ("C", 7.0, 0.94, 6.26),
        ("D", 12.0, 1.08, 7.85))

    def declared_points(self) -> tuple[DeclaredPoint, ...]:
        """Return the declared points A-D for the bin calculation.

        Returns
        -------
        tuple of DeclaredPoint
            Outdoor temperature (degC), capacity (kW) and EER of each point.
        """
        return tuple(DeclaredPoint(T_o, Pdc, EERd) for _, T_o, Pdc, EERd in self.points)

    def low_power_modes(self) -> LowPowerModes:
        """Return the power in the low-power modes, kW (crankcase heater 0: none listed).

        Returns
        -------
        LowPowerModes
            Thermostat-off, standby, off and crankcase-heater power, kW.
        """
        return LowPowerModes(P_TO=self.P_thermostat_off_W / 1e3, P_SB=self.P_standby_W / 1e3,
                             P_OFF=self.P_off_W / 1e3, P_CK=0.0)

    def heating_declared_points(self) -> tuple[DeclaredPoint, ...]:
        """Return the declared heating points A-D (capacity in ``Pdc``, COP in ``EERd``).

        Returns
        -------
        tuple of DeclaredPoint
            Outdoor temperature (degC), capacity (kW) and COP of each point.
        """
        return tuple(DeclaredPoint(T_o, P, cop) for _, T_o, P, cop in self.heating_points)

    def heating_low_power_modes(self) -> LowPowerModes:
        """Return the power in the low-power modes in heating, kW.

        Returns
        -------
        LowPowerModes
            Thermostat-off, standby, off and crankcase-heater power, kW.
        """
        return replace(self.low_power_modes(), P_TO=self.P_thermostat_off_heating_W / 1e3)


PUBLISHED = Datasheet()
