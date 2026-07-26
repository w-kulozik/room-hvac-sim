"""Inverter split air conditioner: cooling coil, capacity, power and start-up.

Commands
    s    compressor speed fraction; 0 = off, otherwise clipped to [s_min, 1]
    f    indoor fan speed fraction; 0 = off, otherwise clipped to [f_min, 1]
         (f > 0 with s = 0 is fan-only operation)

Coil (partially wet). The refrigerant evaporates at a uniform temperature T_r (two-phase; the
short superheated end is neglected). Air and refrigerant exchange heat through two conductances
in series:
    air side          UA_a = NTU_a m cp_m,  NTU_a = NTU_nom (m / m_nom)^(-0.4)   (h ~ Re^0.6)
                      BF = exp(-NTU_a) is the air-side bypass factor (0.15 at full fan)
    refrigerant side  UA_r = UA_evap
On a dry surface the local surface temperature lies between the air and the refrigerant,
    T_w = T_r + r (T - T_r),   r = UA_a / (UA_a + UA_r),
so the surface is warmest where the air enters and cools along the air path. The coil is dry
near the air inlet and wet downstream of the point where T_w reaches the dew point T_dp of the
entering air. This is the part-wet, part-dry structure of the detailed Coil:Cooling:Water model
of the EnergyPlus Engineering Reference (the wet/dry boundary sits where the surface reaches the
air dew point; the wet part uses an enthalpy-based conductance), applied here to a refrigerant
at uniform temperature, for which both parts have closed forms (decision D29). With x_d the
dry fraction of the surface and UA_o = (1/UA_a + 1/UA_r)^-1, NTU_o = UA_o / (m cp_m):
    dry part   T - T_r = (T_in - T_r) exp(-NTU_o x)
               boundary air temperature  T_b = T_r + (T_dp - T_r) / r
               x_d = ln[(T_in - T_r) / (T_b - T_r)] / NTU_o, clipped to [0, 1]
    wet part   (Lewis number 1)  UA_h = 1 / (c_s / UA_r + cp_m / UA_a),  kg/s
               h_out = h_sat(T_r) + (h_b - h_sat(T_r)) exp(-UA_h (1 - x_d) / m)
               c_s is the secant slope of h_sat between T_r and the surface temperature where
               the wet part starts (T_dp, or the inlet surface temperature if the coil is wet
               from the inlet). With this choice the two parts meet exactly at the boundary.
    outlet     effective condensing surface, as in the EnergyPlus model:
               eta = 1 - exp(-NTU_a (1 - x_d)),  h_c = h_b - (h_b - h_out) / eta,
               T_c = T_sat(h_c),  T_out = T_c + (T_b - T_c)(1 - eta),  W_out from (T_out, h_out)
As T_r falls, x_d moves continuously from 1 (dry) to 0 (wet from the inlet), so capacity, outlet
state and water removal are continuous: there is no dry/wet switch. The mean surface
temperature, read by the evaporator thermistor, is T_s = T_r + Q / UA_r.

Sensible/latent split, exact for h = cp_a T + W (h_fg0 + cp_v T):
    h_in - h_out = (cp_a + cp_v W_out)(T_in - T_out) + (W_in - W_out)(h_fg0 + cp_v T_in)
                 =        Q_sens / m              +          Q_lat / m
The refrigerant also takes the enthalpy of the condensate leaving at T_s:
    Q_coil = m (h_in - h_out) - m_lat cp_w T_s.

Compressor capacity. The speed sets the displaced volume; the mass flow, and so the capacity,
rises with the evaporating temperature (suction density) and falls with the outdoor temperature:
    Q_comp(T_r) = s Q_rated [1 + k_e (T_r - T_r,rated)] [1 - k_o (T_o - T_o,rated)]
The coil takes less heat the warmer T_r, so Q_air(T_r) = Q_comp(T_r) has a single root (Brent's
method). Less airflow, or cooler and drier return air, lowers T_r and with it the capacity: the
flow-fraction and entering-wet-bulb dependences that the EnergyPlus DX coil model writes as
curves (CAP-FFF, CAP-FT) follow here from the coil. k_o and k_e are set from the EnergyPlus
example curve WindACCoolCapFT (decision D30), both fitted through the coil at the rating point:
k_e so that the unit's capacity has the curve's entering-wet-bulb slope, then k_o so that the
capacity at fixed entering air has the curve's outdoor slope. k_o is larger than that slope
because less capacity lets the evaporator warm, which gives back part of the loss.

Refrigeration cycle (Carnot fraction with finite-UA condenser):
    T_cond = T_o + Q_cond / UA_cond,   Q_cond = Q_coil (1 + 1/COP)
    COP    = eta(s) T_r / (T_cond - T_r)
Substituting T_cond gives COP in closed form (no iteration):
    COP = (eta T_r - a) / (T_o - T_r + a),   a = Q_coil / UA_cond
eta is calibrated so that the unit meets its rated net EER at rating conditions: the net
capacity (coil capacity minus the indoor fan heat, which ends up in the room) over the total
electrical power, fans included.
Minimum lift: a real cycle cannot let the condensing temperature fall towards the evaporating
one. The compressor needs a minimum pressure ratio and the expansion device a minimum pressure
difference, and head-pressure control (outdoor fan speed) holds them. The cycle therefore keeps
T_cond - T_r >= dT_lift_min; when the formula above would go below it (mild outdoor air, low
load), T_cond = T_r + dT_lift_min and COP = eta T_r / dT_lift_min. The two expressions
meet continuously at the switch. Without this floor the COP exceeded 17 at 20 degC outdoors and
became negative below about 13 degC.

Frost protection: the controller of a real unit slows the compressor when the evaporating
temperature would fall below T_evap_min. Here T_r is held at T_evap_min and the capacity is what
the coil takes at that temperature (only active at low fan speed and high load).

Indoor fan: P = P_0 + (P_nom - P_0) f^3. The cubic part is the affinity law of a fan on a fixed
air path; P_0 is the speed-independent part (motor and drive losses, control electronics),
which dominates at low speed (decision D31).

Water held on the coil (mechanism as described by Henderson & Rengarajan, 1996).
State M_c, kg, capacity M_max:
    compressor on   condensate first wets the coil and drains only once the coil is full:
                    dM_c/dt = m_cond (1 - g),  drain = m_cond g,
                    g = 0 below 0.95 M_max -> 1 at M_max
    fan only        the film re-evaporates into the air stream at constant air enthalpy:
                    W_out = W_in + (1 - BF)(W_sat(T_wb) - W_in) e(M_c),  dM_c/dt = -m_evap
                    with e(M_c) = 1 above 0.05 M_max, falling linearly to 0 when the coil is dry
    both off        no change (natural drying is neglected)
The smooth ramps of g and e keep the right-hand side continuous for the integrator. The
enthalpy of the liquid film (cp_w T_wb ~ 80 kJ/kg vs h_fg ~ 2500 kJ/kg) is neglected in the
fan-only branch, so the supply air there has exactly the return-air enthalpy.
Re-evaporation after the compressor stops with the fan running is the main cause of the
part-load latent degradation of cycling units.

Start-up: after each start the delivered capacity ramps as 1 - exp(-t_on / tau_start) while
the compressor already draws the power of the requested capacity.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import cached_property
from typing import Literal, NamedTuple, overload

from scipy.optimize import brentq

from . import psychro

CP_A = 1006.0       # dry air, J/(kg K)
CP_V = 1860.0       # water vapour, J/(kg K)
CP_W = 4186.0       # liquid water, J/(kg K)
H_FG0 = 2501e3      # latent heat at 0 degC, J/kg
P = psychro.P_ATM
EPS = psychro.EPS
_A, _B, _C = psychro._A, psychro._B, psychro._C


# ---- fast scalar psychrometrics (same equations as psychro.py, J/kg units) ------------------

def _p_ws(T: float) -> float:
    return _A * math.exp(_B * T / (T + _C))


def _W_sat(T: float) -> float:
    pws = _p_ws(T)
    return EPS * pws / (P - pws)


def _h(T: float, W: float) -> float:
    return CP_A * T + W * (H_FG0 + CP_V * T)


def _h_sat_and_slope(T: float) -> tuple[float, float]:
    pws = _p_ws(T)
    Ws = EPS * pws / (P - pws)
    dpws = pws * _B * _C / (T + _C) ** 2
    dWs = EPS * P * dpws / (P - pws) ** 2
    h = CP_A * T + Ws * (H_FG0 + CP_V * T)
    dh = CP_A + CP_V * Ws + (H_FG0 + CP_V * T) * dWs
    return h, dh


def T_from_h_sat(h: float, T_guess: float = 10.0, tol: float = 1e-10) -> float:
    """Return the temperature at which saturated air has enthalpy ``h``.

    Newton iteration on h_sat(T); it usually converges in 3-5 steps.

    Parameters
    ----------
    h : float
        Moist-air enthalpy, J/kg dry air.
    T_guess : float, optional
        Starting temperature, degC.
    tol : float, optional
        Convergence tolerance on the Newton step, K.

    Returns
    -------
    float
        Saturation temperature, degC.

    Raises
    ------
    RuntimeError
        If the iteration has not converged after 30 steps.
    """
    T = T_guess
    for _ in range(30):
        hs, dh = _h_sat_and_slope(T)
        dT = (hs - h) / dh
        T -= dT
        if abs(dT) < tol:
            return T
    raise RuntimeError("saturation temperature did not converge")


def _dew_point(W: float) -> float:
    g = math.log(W * P / (EPS + W) / _A)
    return _C * g / (_B - g)


# ---- unit -------------------------------------------------------------------------------

@dataclass(frozen=True)
class ACParams:
    """Rating data and model parameters of the split unit.

    Attributes
    ----------
    Q_rated : float
        Total (air-side) cooling capacity at rating conditions, W.
    EER_rated : float
        Net EER at rating, (Q - P_fan_in) / P_el with fans included, W/W.
    T_o_rated : float
        Outdoor temperature at rating, degC.
    T_in_rated : float
        Return-air dry-bulb temperature at rating, degC.
    W_in_rated : float
        Return-air humidity ratio at rating (19.7 degC wet bulb), kg/kg.
    cap_slope : float
        Compressor capacity loss per K of outdoor temperature at fixed T_r, 1/K (D30).
    cap_slope_evap : float
        Compressor capacity gain per K of evaporating temperature, 1/K (D30).
    s_min : float
        Minimum compressor speed fraction, -.
    f_min : float
        Minimum indoor fan speed fraction, -.
    V_nom : float
        Indoor airflow at full fan, m3/s.
    BF_nom : float
        Air-side bypass factor at full fan, -.
    UA_evap : float
        Conductance refrigerant to coil surface, W/K.
    UA_cond : float
        Conductance refrigerant to outdoor air, W/K.
    eta_low_speed : float
        Relative loss of cycle efficiency at s_min with respect to full speed, -.
    dT_lift_min : float
        Minimum condensing minus evaporating temperature difference, K.
    P_fan_in_nom : float
        Indoor fan power at full speed, W.
    P_fan_in_base : float
        Speed-independent part of the indoor fan power, W (D31).
    P_fan_out : float
        Outdoor fan power while the compressor runs, W.
    T_evap_min : float
        Frost-protection limit on the evaporating temperature, degC.
    tau_start : float
        Time constant of the capacity ramp after a start, s (D9).
    M_coil_max : float
        Water the coil holds before it drains, kg.
    film_ramp : float
        Width of the full and dry ramps, fraction of M_coil_max, -.
    rho_air : float
        Air density used to convert airflow to mass flow, kg/m3.
    """

    Q_rated: float = 2500.0          # total (air-side) cooling capacity at rating conditions, W
    EER_rated: float = 3.5           # net EER at rating: (Q - P_fan_in) / P_el, fans included
    T_o_rated: float = 35.0          # outdoor, degC
    T_in_rated: float = 27.0         # return air dry bulb, degC
    W_in_rated: float = 0.0114       # return air humidity ratio (19.7 degC wet bulb)
    cap_slope: float = 0.016         # compressor capacity loss per K outdoor, at fixed T_r (D30)
    cap_slope_evap: float = 0.040    # capacity gain per K of evaporating temperature (D30)
    s_min: float = 0.3               # minimum compressor speed fraction
    f_min: float = 0.4               # minimum indoor fan fraction
    V_nom: float = 0.125             # indoor airflow at full fan, m3/s (450 m3/h)
    BF_nom: float = 0.15             # air-side bypass factor at full fan
    UA_evap: float = 600.0           # refrigerant <-> coil surface, W/K
    UA_cond: float = 300.0           # refrigerant <-> outdoor air, W/K
    eta_low_speed: float = 0.25      # relative loss of cycle efficiency at s_min vs full speed
    dT_lift_min: float = 15.0        # minimum condensing - evaporating temperature difference, K
    P_fan_in_nom: float = 25.0       # indoor fan at full speed, W
    P_fan_in_base: float = 5.0       # speed-independent part of the indoor fan power, W (D31)
    P_fan_out: float = 40.0          # outdoor fan while the compressor runs, W
    T_evap_min: float = 0.5          # frost protection limit on evaporating temperature, degC
    tau_start: float = 60.0          # capacity ramp after a start, s (decision D9)
    M_coil_max: float = 0.25         # water the coil holds before it drains, kg
    film_ramp: float = 0.05          # width of the full/dry ramps, fraction of M_coil_max
    rho_air: float = 1.2

    @property
    def m_nom(self) -> float:
        """Return the indoor air mass flow at full fan, kg/s."""
        return self.rho_air * self.V_nom

    @property
    def NTU_nom(self) -> float:
        """Return the air-side NTU at full fan, -ln(BF_nom), dimensionless."""
        return -math.log(self.BF_nom)


class CoilState(NamedTuple):
    """Operating point of the unit at one instant.

    Attributes
    ----------
    Q_sens : float
        Sensible heat removed from the air stream, W.
    m_lat : float
        Water removed from the air, kg/s (negative when the film re-evaporates).
    Q_coil : float
        Heat taken by the refrigerant, W.
    P_el : float
        Total electrical power, W.
    P_fan_in : float
        Indoor fan power, W (ends up as heat in the room).
    COP : float
        Q_coil / P_comp (compressor only), W/W.
    T_s : float
        Mean coil surface temperature, degC.
    T_out : float
        Supply-air temperature, degC.
    W_out : float
        Supply-air humidity ratio, kg/kg.
    wet_fraction : float
        Wetted fraction of the coil surface, 0 (dry) to 1.
    dM_coil : float
        Rate of change of the water held on the coil, kg/s.
    m_drain : float
        Condensate leaving through the drain, kg/s.
    """

    Q_sens: float      # sensible heat removed from the air stream, W
    m_lat: float       # water removed from the air, kg/s (negative when the film re-evaporates)
    Q_coil: float      # heat taken by the refrigerant, W
    P_el: float        # electrical power, W
    P_fan_in: float    # indoor fan power, W (ends up as heat in the room)
    COP: float         # Q_coil / P_comp (compressor only)
    T_s: float         # mean coil surface temperature, degC
    T_out: float       # supply air, degC
    W_out: float       # supply air, kg/kg
    wet_fraction: float  # wetted fraction of the coil surface, 0 (dry) to 1
    dM_coil: float     # rate of change of the water held on the coil, kg/s
    m_drain: float     # condensate leaving through the drain, kg/s


OFF_STATE = CoilState(Q_sens=0.0, m_lat=0.0, Q_coil=0.0, P_el=0.0, P_fan_in=0.0, COP=0.0,
                      T_s=math.nan, T_out=math.nan, W_out=math.nan, wet_fraction=0.0,
                      dM_coil=0.0, m_drain=0.0)


def net_cooling(st: CoilState, T_in: float) -> float:
    """Return the cooling delivered to the room: sensible + latent minus the indoor fan heat.

    Parameters
    ----------
    st : CoilState
        Operating point of the unit.
    T_in : float
        Return-air temperature, degC, at which the latent heat is evaluated.

    Returns
    -------
    float
        Net cooling, W.
    """
    return st.Q_sens + st.m_lat * (H_FG0 + CP_V * T_in) - st.P_fan_in


class CoilOutlet(NamedTuple):
    """Air-side result of the coil model at a given evaporating temperature.

    Attributes
    ----------
    Q : float
        Air-side total capacity m (h_in - h_out), W.
    T_out : float
        Outlet air temperature, degC.
    W_out : float
        Outlet air humidity ratio, kg/kg.
    T_s : float
        Mean coil surface temperature, degC.
    wet_fraction : float
        Wetted fraction of the coil surface, 0 (dry) to 1.
    """

    Q: float          # air-side total capacity m (h_in - h_out), W
    T_out: float
    W_out: float
    T_s: float         # mean surface temperature, degC
    wet_fraction: float


@dataclass(frozen=True)
class ACUnit:
    """Inverter split air conditioner built from the coil, compressor and cycle models.

    Attributes
    ----------
    p : ACParams
        Rating data and model parameters.
    """

    p: ACParams = field(default_factory=ACParams)

    # -- coil -------------------------------------------------------------------------------
    def bypass_factor(self, f: float) -> float:
        """Return the air-side bypass factor exp(-NTU_a).

        Parameters
        ----------
        f : float
            Indoor fan speed fraction, > 0.

        Returns
        -------
        float
            Bypass factor, dimensionless.
        """
        return math.exp(-self.p.NTU_nom * f ** -0.4)

    def fan_power(self, f: float) -> float:
        """Return the indoor fan electrical power.

        Parameters
        ----------
        f : float
            Indoor fan speed fraction, > 0.

        Returns
        -------
        float
            Fan power P_0 + (P_nom - P_0) f^3, W.
        """
        p = self.p
        return p.P_fan_in_base + (p.P_fan_in_nom - p.P_fan_in_base) * f ** 3

    @overload
    def _coil(self, T_r: float, f: float, T_in: float, W_in: float,
              outlet: Literal[False]) -> float: ...

    @overload
    def _coil(self, T_r: float, f: float, T_in: float, W_in: float,
              outlet: Literal[True]) -> CoilOutlet: ...

    def _coil(self, T_r: float, f: float, T_in: float, W_in: float,
              outlet: bool) -> float | CoilOutlet:
        """Coil with refrigerant at T_r. Returns Q (W) or, with ``outlet``, a CoilOutlet."""
        p = self.p
        m = p.m_nom * f
        cpm = CP_A + CP_V * W_in
        NTU_a: float = p.NTU_nom * f ** -0.4
        UA_a = NTU_a * m * cpm
        UA_r = p.UA_evap
        r = UA_a / (UA_a + UA_r)
        NTU_o = NTU_a * UA_r / (UA_a + UA_r)
        dT_in = T_in - T_r
        T_dp = _dew_point(W_in)
        # dry fraction x_d: the surface reaches the dew point where the air is at T_b
        if T_dp <= T_r or dT_in <= 0.0:
            x_d = 1.0
        else:
            T_b = T_r + (T_dp - T_r) / r
            x_d = 0.0 if T_b >= T_in else min(math.log(dT_in / (T_b - T_r)) / NTU_o, 1.0)
        T_b = T_r + dT_in * math.exp(-NTU_o * x_d)        # air at the end of the dry part
        wet = 1.0 - x_d
        if wet <= 0.0:
            Q = m * cpm * (T_in - T_b)
            if not outlet:
                return Q
            return CoilOutlet(Q, T_b, W_in, T_r + Q / UA_r, 0.0)
        # wet part: enthalpy-based conductance with the secant slope of h_sat
        T_up = min(T_dp, T_r + r * dT_in)                 # surface where the wet part starts
        h_r, dh_r = _h_sat_and_slope(T_r)
        c_s = (_h_sat_and_slope(T_up)[0] - h_r) / (T_up - T_r) if T_up - T_r > 1e-6 else dh_r
        UA_h = 1.0 / (c_s / UA_r + cpm / UA_a)
        h_b = _h(T_b, W_in)
        h_out = h_r + (h_b - h_r) * math.exp(-UA_h * wet / m)
        Q = m * (_h(T_in, W_in) - h_out)
        if not outlet:
            return Q
        eta = 1.0 - math.exp(-NTU_a * wet)
        if eta < 1e-9:
            W_out = W_in
            T_out = (h_out - H_FG0 * W_out) / (CP_A + CP_V * W_out)
        else:
            h_c = h_b - (h_b - h_out) / eta
            T_c = T_from_h_sat(h_c, T_guess=T_r)
            T_out = T_c + (T_b - T_c) * (1.0 - eta)
            W_out = (h_out - CP_A * T_out) / (H_FG0 + CP_V * T_out)
            if W_out > W_in:
                W_out = W_in
                T_out = (h_out - H_FG0 * W_out) / (CP_A + CP_V * W_out)
        if W_out > _W_sat(T_out):
            # outlet state beyond saturation: condense the excess at constant enthalpy
            T_out = T_from_h_sat(h_out, T_guess=T_out)
            W_out = _W_sat(T_out)
        return CoilOutlet(Q, T_out, W_out, T_r + Q / UA_r, wet)

    def coil_at(self, T_r: float, f: float, T_in: float, W_in: float) -> CoilOutlet:
        """Return the coil outlet with the refrigerant evaporating at T_r.

        Parameters
        ----------
        T_r : float
            Evaporating temperature, degC.
        f : float
            Indoor fan speed fraction, > 0.
        T_in : float
            Entering (return) air temperature, degC.
        W_in : float
            Entering air humidity ratio, kg/kg.

        Returns
        -------
        CoilOutlet
            Air-side capacity (W), outlet state and surface temperature.
        """
        return self._coil(T_r, f, T_in, W_in, True)

    def T_evap(self, Q: float, f: float, T_in: float, W_in: float) -> float:
        """Return the evaporating temperature at which the coil takes air-side capacity Q.

        Solved with Brent's method on [T_in - 80, T_in].

        Parameters
        ----------
        Q : float
            Air-side total capacity, W. For Q <= 0 the result is T_in.
        f : float
            Indoor fan speed fraction, > 0.
        T_in : float
            Entering air temperature, degC.
        W_in : float
            Entering air humidity ratio, kg/kg.

        Returns
        -------
        float
            Evaporating temperature, degC.
        """
        if Q <= 0.0:
            return T_in
        return float(brentq(lambda T: self._coil(T, f, T_in, W_in, False) - Q, T_in - 80.0,
                            T_in, xtol=1e-12, rtol=1e-14))

    # -- capacity and cycle -----------------------------------------------------------------
    def capacity(self, s: float, T_o: float) -> float:
        """Return the compressor capacity at the rated evaporating temperature.

        Parameters
        ----------
        s : float
            Compressor speed fraction.
        T_o : float
            Outdoor temperature, degC.

        Returns
        -------
        float
            Capacity, W.
        """
        return s * self.p.Q_rated * (1.0 - self.p.cap_slope * (T_o - self.p.T_o_rated))

    def compressor_capacity(self, s: float, T_r: float, T_o: float) -> float:
        """Return the compressor capacity with the refrigerant evaporating at T_r.

        Parameters
        ----------
        s : float
            Compressor speed fraction.
        T_r : float
            Evaporating temperature, degC.
        T_o : float
            Outdoor temperature, degC.

        Returns
        -------
        float
            Capacity Q_comp(T_r), W.
        """
        return self.capacity(s, T_o) * (1.0 + self.p.cap_slope_evap * (T_r - self.T_r_rated))

    @cached_property
    def T_r_rated(self) -> float:
        """Evaporating temperature at rating conditions (full speed and fan), degC."""
        p = self.p
        return self.T_evap(p.Q_rated, 1.0, p.T_in_rated, p.W_in_rated)

    def operating_T_r(self, s: float, f: float, T_in: float, W_in: float, T_o: float) -> float:
        """Return the evaporating temperature where compressor and coil capacities balance.

        Frost protection holds the result at T_r >= T_evap_min.

        Parameters
        ----------
        s : float
            Compressor speed fraction (already clipped to [s_min, 1]).
        f : float
            Indoor fan speed fraction, > 0.
        T_in : float
            Return-air temperature, degC.
        W_in : float
            Return-air humidity ratio, kg/kg.
        T_o : float
            Outdoor temperature, degC.

        Returns
        -------
        float
            Evaporating temperature, degC.
        """
        T_min = self.p.T_evap_min

        def g(T: float) -> float:
            return self._coil(T, f, T_in, W_in, False) - self.compressor_capacity(s, T, T_o)

        if T_min < T_in and g(T_min) <= 0.0:
            return T_min                     # the compressor slows to hold T_evap_min
        return float(brentq(g, T_in - 80.0, T_in, xtol=1e-12, rtol=1e-14))

    def eta_rel(self, s: float) -> float:
        """Return the cycle efficiency relative to full speed.

        Compressor losses grow quadratically as the speed falls towards s_min.

        Parameters
        ----------
        s : float
            Compressor speed fraction.

        Returns
        -------
        float
            Relative efficiency, dimensionless (1 at full speed).
        """
        x = (1.0 - s) / (1.0 - self.p.s_min)
        return 1.0 - self.p.eta_low_speed * x * x

    def _cop(self, eta: float, Q_coil: float, T_r: float, T_o: float) -> float:
        T_ev = T_r + 273.15
        a = Q_coil / self.p.UA_cond
        cop_floor = eta * T_ev / self.p.dT_lift_min          # COP with the minimum lift
        num, den = eta * T_ev - a, T_o + 273.15 - T_ev + a
        if num <= 0.0:
            raise ValueError("condenser too small for this capacity: no positive COP")
        if den <= 0.0:
            return cop_floor
        cop = num / den
        T_cond = T_o + 273.15 + a * (1.0 + 1.0 / cop)
        return cop if T_cond - T_ev >= self.p.dT_lift_min else cop_floor

    @cached_property
    def eta_full(self) -> float:
        """Carnot fraction at full speed, calibrated to the rated net EER (fans included)."""
        p = self.p
        T_r = self.T_r_rated
        out = self.coil_at(T_r, 1.0, p.T_in_rated, p.W_in_rated)
        m_lat = p.m_nom * (p.W_in_rated - out.W_out)
        Q_coil = out.Q - m_lat * CP_W * out.T_s
        P_fan = self.fan_power(1.0)
        P_comp = (out.Q - P_fan) / p.EER_rated - P_fan - p.P_fan_out
        cop = Q_coil / P_comp
        T_ev = T_r + 273.15
        a = Q_coil / p.UA_cond
        return (cop * (p.T_o_rated + 273.15 - T_ev + a) + a) / T_ev

    # -- full operating point ---------------------------------------------------------------
    def streams(
        self, T_r: float, f: float, T_in: float, W_in: float
    ) -> tuple[float, float, float, float, float, float, float]:
        """Return the heat and moisture streams with the refrigerant at T_r.

        Parameters
        ----------
        T_r : float
            Evaporating temperature, degC.
        f : float
            Indoor fan speed fraction, > 0.
        T_in : float
            Return-air temperature, degC.
        W_in : float
            Return-air humidity ratio, kg/kg.

        Returns
        -------
        tuple of float
            (Q_sens [W], m_lat [kg/s], Q_coil [W], T_s [degC], T_out [degC],
            W_out [kg/kg], wet_fraction [-]).
        """
        out = self.coil_at(T_r, f, T_in, W_in)
        m = self.p.m_nom * f
        m_lat = m * (W_in - out.W_out)
        Q_lat = m_lat * (H_FG0 + CP_V * T_in)
        Q_sens = m * (CP_A + CP_V * out.W_out) * (T_in - out.T_out)
        Q_coil = Q_sens + Q_lat - m_lat * CP_W * out.T_s
        return Q_sens, m_lat, Q_coil, out.T_s, out.T_out, out.W_out, out.wet_fraction

    def _drain_fraction(self, M_c: float) -> float:
        w = self.p.film_ramp * self.p.M_coil_max
        return min(max((M_c - (self.p.M_coil_max - w)) / w, 0.0), 1.0)

    def _evaporation_fraction(self, M_c: float) -> float:
        w = self.p.film_ramp * self.p.M_coil_max
        return min(max(M_c / w, 0.0), 1.0)

    def fan_only(self, f: float, T_in: float, W_in: float, M_c: float) -> CoilState:
        """Return the operating point with the fan running and the compressor off.

        The water film re-evaporates into the air stream at constant air enthalpy.

        Parameters
        ----------
        f : float
            Indoor fan speed fraction, clipped to [f_min, 1].
        T_in : float
            Return-air temperature, degC.
        W_in : float
            Return-air humidity ratio, kg/kg.
        M_c : float
            Water held on the coil, kg.

        Returns
        -------
        CoilState
            Operating point; m_lat <= 0 and dM_coil = m_lat.
        """
        p = self.p
        f = min(max(f, p.f_min), 1.0)
        m = p.m_nom * f
        P_fan_in = self.fan_power(f)
        h_in = _h(T_in, W_in)
        T_wb = T_from_h_sat(h_in, T_guess=T_in - 5.0)
        W_out = W_in + (1.0 - self.bypass_factor(f)) * (_W_sat(T_wb) - W_in) \
            * self._evaporation_fraction(M_c)
        T_out = (h_in - H_FG0 * W_out) / (CP_A + CP_V * W_out)
        m_lat = m * (W_in - W_out)                       # <= 0: water added to the air
        Q_sens = m * (CP_A + CP_V * W_out) * (T_in - T_out)
        # a wet coil sits at the wet-bulb temperature; a dry one at the return-air temperature
        T_surface = T_wb if M_c > 0.0 else T_in
        return CoilState(Q_sens, m_lat, 0.0, P_fan_in, P_fan_in, 0.0, T_surface, T_out, W_out,
                         1.0 if M_c > 0.0 else 0.0, dM_coil=m_lat, m_drain=0.0)

    def operate(self, s: float, f: float, T_in: float, W_in: float, T_o: float,
                ramp: float = 1.0, M_c: float | None = None) -> CoilState:
        """Return the operating point for commands (s, f).

        Parameters
        ----------
        s : float
            Compressor speed fraction; <= 0 is off, otherwise clipped to [s_min, 1].
        f : float
            Indoor fan speed fraction; <= 0 is off, otherwise clipped to [f_min, 1].
        T_in : float
            Return-air temperature, degC.
        W_in : float
            Return-air humidity ratio, kg/kg.
        T_o : float
            Outdoor temperature, degC.
        ramp : float, optional
            Fraction in [0, 1] of the requested capacity delivered after a start; the
            compressor already draws the power of the requested capacity.
        M_c : float or None, optional
            Water held on the coil, kg. None treats the coil as permanently full (steady
            wet-coil operation).

        Returns
        -------
        CoilState
            Operating point of the unit (OFF_STATE when both s and f are <= 0).
        """
        p = self.p
        if s <= 0.0:
            if f <= 0.0:
                return OFF_STATE
            return self.fan_only(f, T_in, W_in, 0.0 if M_c is None else M_c)
        s = min(max(s, p.s_min), 1.0)
        f = min(max(f, p.f_min), 1.0)
        T_r_req = self.operating_T_r(s, f, T_in, W_in, T_o)
        req = self.streams(T_r_req, f, T_in, W_in)
        if ramp >= 1.0:
            out = req
        else:
            Q_req = self._coil(T_r_req, f, T_in, W_in, False)
            out = self.streams(self.T_evap(Q_req * ramp, f, T_in, W_in), f, T_in, W_in)
        Q_sens, m_lat, Q_coil, T_s, T_out, W_out, wet = out
        cop = self._cop(self.eta_full * self.eta_rel(s), req[2], T_r_req, T_o)
        P_fan_in = self.fan_power(f)
        P_el = req[2] / cop + P_fan_in + p.P_fan_out
        g = 1.0 if M_c is None else self._drain_fraction(M_c)
        return CoilState(Q_sens, m_lat, Q_coil, P_el, P_fan_in, cop, T_s, T_out, W_out, wet,
                         dM_coil=m_lat * (1.0 - g), m_drain=m_lat * g)

    def start_ramp(self, t_on: float) -> float:
        """Return the capacity fraction after a compressor start.

        Parameters
        ----------
        t_on : float
            Time since the start, s (negative values count as 0).

        Returns
        -------
        float
            1 - exp(-t_on / tau_start), dimensionless.
        """
        return 1.0 - math.exp(-max(t_on, 0.0) / self.p.tau_start)
