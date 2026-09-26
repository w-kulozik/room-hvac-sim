"""Heat pump: the split unit run in heating, with frost on the outdoor coil and defrost.

Decisions D50-D51. The heating mode reuses the cooling calibration of ``ACUnit`` (Carnot
fraction, minimum lift, evaporating-temperature slope, fans, start-up ramp) instead of being
fitted to heating data.

Commands, as for the cooling unit: s compressor speed (0 off; clipped to [s_min, s_max] in
heating; s < 0 is a reverse-cycle defrost), f indoor fan speed.

Indoor coil (condenser), dry:
    Q_h = G (T_c - T_in),  G = eps m c_p,  eps = 1 - exp(-NTU_o),
    NTU_o = UA_o / (m c_p),  UA_o = (1/UA_a + 1/UA_r)^-1   (UA_a and UA_r as in cooling)
Outdoor coil (evaporator), air at T_o, surface T_s = T_e + Q_e / UA_r,o:
    Q_e = phi UA_a,o (T_o - T_s) + phi (UA_a,o / c_p) (W_o - W_sat(T_s))+ (h_fg + chi h_fus)
    phi = 1 / (1 + M_f / M_half)       frost blockage
    chi = fraction of the deposit that freezes: 1 below -0.5 degC, 0 above +0.5 degC, linear
          in between; W_sat blends the saturation over ice and over water with the same weight
Compressor (evaporator capacity):
    Q_e = s Q_rated exp[k_e (T_e - T_r,rated)] [1 - k_o (T_c - T_c,rated)]
with T_c,rated the condensing temperature of the cooling rating point (D57).
Cycle, temperatures in kelvin:
    COP_h = eta(s) T_c / max(T_c - T_e, dT_lift_min),   Q_h = Q_e + P_comp,  P_comp = Q_h / COP_h
With T_c = T_in + Q_h / G this is a quadratic in Q_h for each T_e, solved in closed form; T_e
follows from one root search on the outdoor-coil balance.
Frost: dM_f/dt = chi m_dep while heating; it melts with the compressor off in air above 0 degC
(natural convection) and during a defrost. Reverse-cycle defrost (s < 0): fixed electrical
power P_def, heat Q_def to the outdoor coil of which eta_melt melts frost; Q_def - P_def is
taken from the room.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import cached_property
from typing import NamedTuple

from scipy.optimize import brentq

from . import psychro
from .ac_unit import CP_A, CP_V, H_FG0, ACParams, ACUnit, CoilState

H_FUS = 334e3          # latent heat of fusion of ice, J/kg
T0_K = 273.15


@dataclass(frozen=True)
class HeatPumpParams:
    """Parameters of the heating mode that the cooling model does not have (D50-D51).

    Attributes
    ----------
    s_max : float
        Highest compressor speed in heating, -; the reference unit's ratio of maximum to
        nominal heating capacity.
    UA_air_out : float
        Air-side conductance of the clean outdoor coil, W/K.
    UA_ref_out : float
        Refrigerant-side conductance of the outdoor coil, W/K. In series with the air side
        it gives the cooling model's condenser conductance UA_cond.
    M_half : float
        Frost mass that halves the air-side conductance, kg.
    P_defrost : float
        Electrical power during a reverse-cycle defrost, W.
    Q_defrost : float
        Heat delivered to the outdoor coil during a defrost, W.
    eta_melt : float
        Fraction of Q_defrost that melts frost, -.
    T_coil_defrosted : float
        Outdoor-coil thermistor reading in defrost once the frost has melted, degC.
    UA_natural_out : float
        Conductance from the outdoor air to a frosted coil with the outdoor fan stopped
        (natural convection), W/K.
    frost_ramp : float
        Frost mass over which melting and the defrost thermistor reading go smoothly to their
        frost-free values, kg. It keeps the melting rate over this mass below the stability
        limit of RK4 at a 60 s step (D57).
    P_fan_out : float
        Outdoor fan power while the compressor runs, W.
    """

    s_max: float = 1.68
    UA_air_out: float = 600.0
    UA_ref_out: float = 600.0
    M_half: float = 0.5
    P_defrost: float = 700.0
    Q_defrost: float = 2000.0
    eta_melt: float = 0.5
    T_coil_defrosted: float = 10.0
    UA_natural_out: float = 60.0
    frost_ramp: float = 0.1
    P_fan_out: float = 40.0


class OutdoorCoil(NamedTuple):
    """Operating point of the outdoor coil in heating.

    Attributes
    ----------
    T_e : float
        Evaporating temperature, degC.
    T_s : float
        Coil surface temperature, degC.
    Q_e : float
        Heat taken from the outdoor air, W.
    m_dep : float
        Water deposited on the coil (frost and condensate), kg/s.
    chi : float
        Fraction of the deposit that freezes, -.
    """

    T_e: float
    T_s: float
    Q_e: float
    m_dep: float
    chi: float


def _frozen_fraction(T_s: float) -> float:
    """Return the fraction of the deposit that freezes: 1 below -0.5 degC, 0 above 0.5 degC."""
    return min(max(0.5 - T_s, 0.0), 1.0)


def W_sat_surface(T_s: float) -> float:
    """Return the humidity ratio at saturation over a coil surface, over ice or water.

    Parameters
    ----------
    T_s : float
        Surface temperature, degC.

    Returns
    -------
    float
        Saturation humidity ratio at standard pressure, kg/kg; over ice below -0.5 degC, over
        water above 0.5 degC, blended linearly in between.
    """
    chi = _frozen_fraction(T_s)
    p_w = chi * float(psychro.p_wsi(T_s)) + (1.0 - chi) * float(psychro.p_ws(T_s))
    return float(psychro.W_from_pw(p_w))


@dataclass(frozen=True)
class HeatPump:
    """Split unit in heating mode, built on the cooling model ``ACUnit``.

    Attributes
    ----------
    unit : ACUnit
        Cooling model whose calibration is reused.
    hp : HeatPumpParams
        Heating-only parameters.
    """

    unit: ACUnit = field(default_factory=ACUnit)
    hp: HeatPumpParams = field(default_factory=HeatPumpParams)

    @property
    def p(self) -> ACParams:
        """Parameters of the underlying unit (s_min, fans, start-up)."""
        return self.unit.p

    # -- components ---------------------------------------------------------------------------
    def fan_power(self, f: float) -> float:
        """Return the indoor fan power at speed f (as in cooling), W."""
        return self.unit.fan_power(f)

    def start_ramp(self, t_on: float) -> float:
        """Return the capacity fraction t_on seconds after a start (as in cooling), -."""
        return self.unit.start_ramp(t_on)

    def indoor_conductance(self, f: float, W_in: float) -> float:
        """Return G = eps m c_p of the indoor coil as a condenser, W/K.

        Parameters
        ----------
        f : float
            Indoor fan speed fraction, > 0.
        W_in : float
            Return-air humidity ratio, kg/kg.

        Returns
        -------
        float
            Heat delivered per kelvin of condensing temperature above the return air, W/K.
        """
        p = self.p
        m = p.m_nom * f
        cpm = CP_A + CP_V * W_in
        UA_a = p.NTU_nom * f ** -0.4 * m * cpm
        UA_o = 1.0 / (1.0 / UA_a + 1.0 / p.UA_evap)
        return (1.0 - math.exp(-UA_o / (m * cpm))) * m * cpm

    def evaporator_capacity(self, s: float, T_e: float) -> float:
        """Return the compressor's evaporator capacity, s Q_rated exp[k_e (T_e - T_r,rated)], W.

        Parameters
        ----------
        s : float
            Compressor speed fraction, -.
        T_e : float
            Evaporating temperature, degC.

        Returns
        -------
        float
            Heat the compressor draws from the evaporator, W.
        """
        p = self.p
        return s * p.Q_rated * math.exp(p.cap_slope_evap * (T_e - self.unit.T_r_rated))

    def outdoor_air(self, T_s: float, T_o: float, W_o: float, M_f: float) -> tuple[float, float,
                                                                                    float]:
        """Return the heat and water the outdoor air gives to a coil surface at T_s.

        Parameters
        ----------
        T_s : float
            Coil surface temperature, degC.
        T_o : float
            Outdoor air temperature, degC.
        W_o : float
            Outdoor humidity ratio, kg/kg.
        M_f : float
            Frost held on the coil, kg.

        Returns
        -------
        Q : float
            Sensible plus latent heat to the surface, W.
        m_dep : float
            Water deposited (frost and condensate), kg/s.
        chi : float
            Fraction of the deposit that freezes, -.
        """
        h = self.hp
        phi = 1.0 / (1.0 + max(M_f, 0.0) / h.M_half)
        chi = _frozen_fraction(T_s)
        m_dep = phi * h.UA_air_out / CP_A * max(W_o - W_sat_surface(T_s), 0.0)
        Q = phi * h.UA_air_out * (T_o - T_s) + m_dep * (H_FG0 + chi * H_FUS)
        return Q, m_dep, chi

    @cached_property
    def T_c_rated(self) -> float:
        """Condensing temperature of the cooling rating point, degC (reference of k_o)."""
        p, u = self.p, self.unit
        st = u.operate(1.0, 1.0, p.T_in_rated, p.W_in_rated, p.T_o_rated)
        return p.T_o_rated + st.Q_coil * (1.0 + 1.0 / st.COP) / p.UA_cond

    def eta_rel(self, s: float) -> float:
        """Return the relative cycle efficiency, the cooling law mirrored above full speed, -."""
        return self.unit.eta_rel(s if s <= 1.0 else 2.0 - s)

    def cycle(self, s: float, T_e: float, G: float, T_in: float) -> tuple[float, float, float]:
        """Return the heat delivered, the evaporator load and the COP at an evaporating temperature.

        The compressor's evaporator capacity falls with the condensing temperature,
        Q_e = A(T_e) [1 - k_o (T_c - T_c,rated)] = alpha - beta Q_h, and the cycle gives
        Q_h ((eta - 1) T_c + T_e) = eta T_c Q_e (or Q_h (eta T_c - lift) = eta T_c Q_e at the
        minimum lift), with T_c = T_in + Q_h / G: a quadratic in Q_h.

        Parameters
        ----------
        s : float
            Compressor speed fraction, -.
        T_e : float
            Evaporating temperature, degC.
        G : float
            Indoor-coil conductance (``indoor_conductance``), W/K.
        T_in : float
            Return-air temperature, degC.

        Returns
        -------
        Q_h : float
            Heat delivered by the condenser, W.
        Q_e : float
            Heat drawn from the evaporator, W.
        COP_h : float
            Heating COP of the compressor, Q_h / P_comp, -.
        """
        p = self.p
        eta = self.unit.eta_full * self.eta_rel(s)
        lift = p.dT_lift_min
        Ti, Te = T_in + T0_K, T_e + T0_K
        A = self.evaporator_capacity(s, T_e)
        alpha = A * (1.0 - p.cap_slope * (T_in - self.T_c_rated))
        beta = A * p.cap_slope / G
        a = ((eta - 1.0) + eta * beta) / G
        b = (eta - 1.0) * Ti + Te - eta * alpha / G + eta * beta * Ti
        c = -eta * alpha * Ti
        disc = b * b - 4.0 * a * c
        if disc < 0.0:
            raise ValueError(f"no heating cycle at T_e = {T_e!r} degC, T_in = {T_in!r} degC, "
                             f"s = {s!r}: the lift exceeds what the Carnot fraction allows")
        Q_h = (-b + math.sqrt(disc)) / (2.0 * a)
        T_c = Ti + Q_h / G
        if T_c - Te < lift:
            a = eta * (1.0 + beta) / G
            b = eta * Ti - lift - eta * alpha / G + eta * beta * Ti
            Q_h = (-b + math.sqrt(b * b - 4.0 * a * c)) / (2.0 * a)
            T_c = Ti + Q_h / G
            return Q_h, alpha - beta * Q_h, eta * T_c / lift
        return Q_h, alpha - beta * Q_h, eta * T_c / (T_c - Te)

    def outdoor_coil(self, s: float, G: float, T_in: float, T_o: float, W_o: float,
                     M_f: float) -> tuple[OutdoorCoil, float, float]:
        """Return the operating point where the outdoor coil and the cycle balance.

        Parameters
        ----------
        s : float
            Compressor speed fraction, -.
        G : float
            Indoor-coil conductance, W/K.
        T_in : float
            Return-air temperature, degC.
        T_o : float
            Outdoor air temperature, degC.
        W_o : float
            Outdoor humidity ratio, kg/kg.
        M_f : float
            Frost held on the coil, kg.

        Returns
        -------
        coil : OutdoorCoil
            Evaporating and surface temperatures, heat taken and water deposited.
        Q_h : float
            Heat delivered by the condenser, W.
        COP_h : float
            Heating COP of the compressor, -.
        """
        UA_r = self.hp.UA_ref_out

        def balance(T_e: float) -> float:
            Q_e = self.cycle(s, T_e, G, T_in)[1]
            return Q_e - self.outdoor_air(T_e + Q_e / UA_r, T_o, W_o, M_f)[0]

        # T_e lies a few kelvin below T_o (the coil's approach at full load is under 15 K); the
        # bracket stops 40 K below, short of temperatures at which the cycle has no solution
        T_e = float(brentq(balance, T_o - 40.0, T_o, xtol=1e-12, rtol=1e-14))
        Q_h, Q_e, cop = self.cycle(s, T_e, G, T_in)
        T_s = T_e + Q_e / UA_r
        _, m_dep, chi = self.outdoor_air(T_s, T_o, W_o, M_f)
        return OutdoorCoil(T_e, T_s, Q_e, m_dep, chi), Q_h, cop

    def _frost_weight(self, M_f: float) -> float:
        """Return 0 without frost, rising linearly to 1 at ``frost_ramp`` kg."""
        return min(max(M_f / self.hp.frost_ramp, 0.0), 1.0)

    # -- operating point ------------------------------------------------------------------------
    def operate(self, s: float, f: float, T_in: float, W_in: float, T_o: float,
                ramp: float = 1.0, M_c: float | None = None,
                W_o: float | None = None) -> CoilState:
        """Return the operating point in heating for commands (s, f).

        Signs follow ``CoilState`` (heat removed from the air): heating gives Q_sens < 0 and
        Q_coil < 0. T_s is the outdoor-coil surface temperature (the defrost thermistor) and
        dM_coil the rate of change of the frost on the outdoor coil.

        Parameters
        ----------
        s : float
            Compressor command: 0 off, > 0 heating (clipped to [s_min, s_max]), < 0 defrost.
        f : float
            Indoor fan speed fraction; clipped to [f_min, 1] while heating, ignored (off) in
            defrost.
        T_in : float
            Return-air temperature, degC.
        W_in : float
            Return-air humidity ratio, kg/kg.
        T_o : float
            Outdoor temperature, degC.
        ramp : float, optional
            Fraction of the requested heat delivered after a start, -.
        M_c : float or None, optional
            Frost on the outdoor coil, kg; None is a clean coil.
        W_o : float or None, optional
            Outdoor humidity ratio, kg/kg; None is dry outdoor air (no frost).

        Returns
        -------
        CoilState
            Operating point of the unit.
        """
        p, h = self.p, self.hp
        M_f = 0.0 if M_c is None else max(M_c, 0.0)
        W_o_ = 0.0 if W_o is None else W_o
        e = self._frost_weight(M_f)
        if s < 0.0:
            # reverse-cycle defrost, indoor fan off: frost melts, the room gives Q_def - P_def
            melt = h.eta_melt * h.Q_defrost / H_FUS * e
            Q_room = h.Q_defrost - h.P_defrost
            return CoilState(Q_room, 0.0, Q_room, h.P_defrost, 0.0, 0.0,
                             h.T_coil_defrosted * (1.0 - e), T_in, W_in, 0.0,
                             dM_coil=-melt, m_drain=0.0)
        if s == 0.0:
            # compressor off: frost melts in air above 0 degC; the fan may run (no coil effect)
            melt = h.UA_natural_out * max(T_o, 0.0) / H_FUS * e
            P_fan = self.fan_power(min(max(f, p.f_min), 1.0)) if f > 0.0 else 0.0
            return CoilState(0.0, 0.0, 0.0, P_fan, P_fan, 0.0, T_o, T_in, W_in, 0.0,
                             dM_coil=-melt, m_drain=0.0)
        s = min(max(s, p.s_min), h.s_max)
        f = min(max(f, p.f_min), 1.0)
        G = self.indoor_conductance(f, W_in)
        out, Q_h, cop = self.outdoor_coil(s, G, T_in, T_o, W_o_, M_f)
        P_fan_in = self.fan_power(f)
        P_el = Q_h / cop + P_fan_in + h.P_fan_out
        Q_del = Q_h * ramp
        m = p.m_nom * f
        T_out = T_in + Q_del / (m * (CP_A + CP_V * W_in))
        return CoilState(-Q_del, 0.0, -Q_del, P_el, P_fan_in, cop, out.T_s, T_out, W_in, 0.0,
                         dM_coil=out.chi * out.m_dep, m_drain=0.0)

    def capacity(self, s: float, T_o: float, T_in: float = 20.0) -> float:
        """Return the heat delivered at speed s with full fan, a clean coil and dry outdoor air.

        Parameters
        ----------
        s : float
            Compressor speed fraction, -.
        T_o : float
            Outdoor temperature, degC.
        T_in : float, optional
            Return-air temperature, degC.

        Returns
        -------
        float
            Heating capacity, W.
        """
        return -self.operate(s, 1.0, T_in, 0.0073, T_o).Q_sens

    def clean_coil_temperature(self, s: float, T_in: float, T_o: float) -> float:
        """Return what the defrost thermistor reads on a clean coil in dry outdoor air, degC.

        The reference of the demand defrost (D52): the unit's controller knows its own clean
        characteristic but has no outdoor humidity sensor.

        Parameters
        ----------
        s : float
            Compressor speed fraction, -.
        T_in : float
            Return-air temperature, degC.
        T_o : float
            Outdoor temperature, degC.
        """
        s = min(max(s, self.p.s_min), self.hp.s_max)
        return self.outdoor_coil(s, self.indoor_conductance(1.0, 0.0073), T_in, T_o, 0.0,
                                 0.0)[0].T_s

    def sensors(self, st: CoilState, T_a: float, T_o: float) -> list[float]:
        """Return the unit's own sensor readings appended to the room measurement.

        Parameters
        ----------
        st : CoilState
            Operating point.
        T_a : float
            Indoor air temperature, degC (not used).
        T_o : float
            Outdoor temperature, degC.

        Returns
        -------
        list of float
            [outdoor-coil thermistor degC, outdoor-air thermistor degC].
        """
        return [st.T_s, T_o]
