"""Inverter split air conditioner: cooling coil, capacity, power and start-up.

Commands
    s    compressor speed fraction; 0 = off, otherwise clipped to [s_min, 1]
    f    indoor fan speed fraction; 0 = off, otherwise clipped to [f_min, 1]
         (f > 0 with s = 0 is fan-only operation)

Coil (bypass-factor / apparatus-dew-point model, equivalent to an effectiveness-NTU coil
with Lewis number 1 and a uniform surface temperature T_s):
    BF   = exp(-NTU),  NTU = NTU_nom * (m / m_nom)^(-0.4)   (air-side h ~ Re^0.6)
    dry  T_out = T_s + BF (T_in - T_s),          W_out = W_in
    wet  h_out = h_sat(T_s) + BF (h_in - h_sat(T_s)),
         W_out = W_sat(T_s) + BF (W_in - W_sat(T_s))
For a requested total capacity Q the surface temperature follows from
    dry  Q = m cp_m (1 - BF)(T_in - T_s)         (closed form)
    wet  Q = m (1 - BF)(h_in - h_sat(T_s))       (Newton on h_sat, analytic derivative)
The coil is dry when the dry solution gives T_s >= dew point of the inlet air. Both branches
coincide at T_s = T_dew,in, so capacity and outlet state are continuous across the switch.

Sensible/latent split, exact for h = cp_a T + W (h_fg0 + cp_v T):
    h_in - h_out = (cp_a + cp_v W_out)(T_in - T_out) + (W_in - W_out)(h_fg0 + cp_v T_in)
                 =        Q_sens / m              +          Q_lat / m
The refrigerant also takes the enthalpy of the condensate leaving at T_s:
    Q_coil = m (h_in - h_out) - m_lat cp_w T_s.

Refrigeration cycle (Carnot fraction with finite-UA heat exchangers):
    T_evap = T_s - Q_coil / UA_evap
    T_cond = T_o + Q_cond / UA_cond,   Q_cond = Q_coil (1 + 1/COP)
    COP    = eta(s) T_evap / (T_cond - T_evap)
Substituting T_cond gives COP in closed form (no iteration):
    COP = (eta T_evap - a) / (T_o - T_evap + a),   a = Q_coil / UA_cond
eta is calibrated so that the unit meets its rated EER (fans included) at rating conditions.

Frost protection: the controller of a real unit slows the compressor when the evaporating
temperature would fall below T_evap_min. Here the requested capacity is capped at the value
that gives T_evap = T_evap_min (bisection; only active at low fan speed and high load).

Water held on the coil (Henderson & Rengarajan, 1996). State M_c, kg, capacity M_max:
    compressor on   condensate first wets the coil and drains only once the coil is full:
                    dM_c/dt = m_cond (1 - g),  drain = m_cond g,  g = 0 below 0.95 M_max -> 1 at M_max
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
from typing import NamedTuple

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
    """Temperature at which saturated air has enthalpy h (J/kg). Newton, converges in 3-5 steps."""
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
    Q_rated: float = 2500.0          # total cooling capacity at rating conditions, W
    EER_rated: float = 3.5           # Q_coil / P_el at rating conditions (fans included)
    T_o_rated: float = 35.0          # outdoor, degC
    T_in_rated: float = 27.0         # return air dry bulb, degC
    W_in_rated: float = 0.0114       # return air humidity ratio (about 19 degC wet bulb)
    cap_slope: float = 0.01          # capacity loss per K of outdoor temperature above rating
    s_min: float = 0.3               # minimum compressor speed fraction
    f_min: float = 0.4               # minimum indoor fan fraction
    V_nom: float = 0.125             # indoor airflow at full fan, m3/s (450 m3/h)
    BF_nom: float = 0.15             # coil bypass factor at full fan
    UA_evap: float = 600.0           # refrigerant <-> coil surface, W/K
    UA_cond: float = 300.0           # refrigerant <-> outdoor air, W/K
    eta_low_speed: float = 0.25      # relative loss of cycle efficiency at s_min vs full speed
    P_fan_in_nom: float = 25.0       # indoor fan at full speed, W (cubic in speed)
    P_fan_out: float = 40.0          # outdoor fan while the compressor runs, W
    T_evap_min: float = 0.5          # frost protection limit on evaporating temperature, degC
    tau_start: float = 60.0          # capacity ramp after a start, s (40-90 s measured on small units)
    M_coil_max: float = 0.25         # water the coil holds before it drains, kg
    film_ramp: float = 0.05          # width of the full/dry ramps, fraction of M_coil_max
    rho_air: float = 1.2

    @property
    def m_nom(self) -> float:
        return self.rho_air * self.V_nom

    @property
    def NTU_nom(self) -> float:
        return -math.log(self.BF_nom)


class CoilState(NamedTuple):
    Q_sens: float      # sensible heat removed from the air stream, W
    m_lat: float       # water removed from the air, kg/s (negative when the film re-evaporates)
    Q_coil: float      # heat taken by the refrigerant, W
    P_el: float        # electrical power, W
    P_fan_in: float    # indoor fan power, W (ends up as heat in the room)
    COP: float         # Q_coil / P_comp (compressor only)
    T_s: float         # coil surface temperature, degC
    T_out: float       # supply air, degC
    W_out: float       # supply air, kg/kg
    wet: bool
    dM_coil: float     # rate of change of the water held on the coil, kg/s
    m_drain: float     # condensate leaving through the drain, kg/s


OFF_STATE = CoilState(Q_sens=0.0, m_lat=0.0, Q_coil=0.0, P_el=0.0, P_fan_in=0.0, COP=0.0,
                      T_s=math.nan, T_out=math.nan, W_out=math.nan, wet=False,
                      dM_coil=0.0, m_drain=0.0)


@dataclass(frozen=True)
class ACUnit:
    p: ACParams = field(default_factory=ACParams)

    # -- coil -------------------------------------------------------------------------------
    def bypass_factor(self, f: float) -> float:
        m = self.p.m_nom * f
        return math.exp(-self.p.NTU_nom * (m / self.p.m_nom) ** -0.4)

    def coil(self, Q: float, f: float, T_in: float, W_in: float) -> tuple:
        """Outlet state for total capacity Q (W) at fan fraction f.

        Returns (T_s, T_out, W_out, wet, m_air).
        """
        m = self.p.m_nom * f
        BF = self.bypass_factor(f)
        cpm = CP_A + CP_V * W_in
        T_s = T_in - Q / (m * cpm * (1.0 - BF))
        if T_s >= _dew_point(W_in):
            return T_s, T_s + BF * (T_in - T_s), W_in, False, m
        h_in = _h(T_in, W_in)
        h_s = h_in - Q / (m * (1.0 - BF))
        T_s = T_from_h_sat(h_s, T_guess=T_s)
        W_s = _W_sat(T_s)
        W_out = W_s + BF * (W_in - W_s)
        h_out = h_s + BF * (h_in - h_s)
        T_out = (h_out - H_FG0 * W_out) / (CP_A + CP_V * W_out)
        if W_out > _W_sat(T_out):
            # mixing line crosses the saturation curve near the surface: condense the excess
            T_out = T_from_h_sat(h_out, T_guess=T_out)
            W_out = _W_sat(T_out)
        return T_s, T_out, W_out, True, m

    # -- capacity and cycle -----------------------------------------------------------------
    def capacity(self, s: float, T_o: float) -> float:
        return s * self.p.Q_rated * (1.0 - self.p.cap_slope * (T_o - self.p.T_o_rated))

    def T_evap(self, Q: float, f: float, T_in: float, W_in: float) -> float:
        """Evaporating temperature at air-side capacity Q, degC."""
        return self.coil(Q, f, T_in, W_in)[0] - Q / self.p.UA_evap

    def frost_limited(self, Q: float, f: float, T_in: float, W_in: float) -> float:
        """Largest capacity <= Q with T_evap >= T_evap_min (T_evap falls monotonically with Q)."""
        T_min = self.p.T_evap_min
        if self.T_evap(Q, f, T_in, W_in) >= T_min:
            return Q
        lo, hi = 0.0, Q
        if self.T_evap(lo, f, T_in, W_in) <= T_min:
            return 0.0
        while hi - lo > 1e-6 * Q:
            mid = 0.5 * (lo + hi)
            if self.T_evap(mid, f, T_in, W_in) >= T_min:
                lo = mid
            else:
                hi = mid
        return lo

    def eta_rel(self, s: float) -> float:
        """Cycle efficiency relative to full speed; compressor losses grow at low speed."""
        x = (1.0 - s) / (1.0 - self.p.s_min)
        return 1.0 - self.p.eta_low_speed * x * x

    def _cop(self, eta: float, Q_coil: float, T_s: float, T_o: float) -> float:
        T_ev = T_s - Q_coil / self.p.UA_evap + 273.15
        a = Q_coil / self.p.UA_cond
        return (eta * T_ev - a) / (T_o + 273.15 - T_ev + a)

    @cached_property
    def eta_full(self) -> float:
        """Carnot fraction at full speed, calibrated to the rated EER (fans included)."""
        p = self.p
        T_s, T_out, W_out, _, m = self.coil(p.Q_rated, 1.0, p.T_in_rated, p.W_in_rated)
        m_lat = m * (p.W_in_rated - W_out)
        Q_coil = p.Q_rated - m_lat * CP_W * T_s
        P_comp = Q_coil / p.EER_rated - p.P_fan_in_nom - p.P_fan_out
        cop = Q_coil / P_comp
        T_ev = T_s - Q_coil / p.UA_evap + 273.15
        a = Q_coil / p.UA_cond
        return (cop * (p.T_o_rated + 273.15 - T_ev + a) + a) / T_ev

    # -- full operating point ---------------------------------------------------------------
    def streams(self, Q: float, f: float, T_in: float, W_in: float) -> tuple:
        """(Q_sens, m_lat, Q_coil, T_s, T_out, W_out, wet) for total air-side capacity Q."""
        T_s, T_out, W_out, wet, m = self.coil(Q, f, T_in, W_in)
        m_lat = m * (W_in - W_out)
        Q_lat = m_lat * (H_FG0 + CP_V * T_in)
        Q_sens = m * (CP_A + CP_V * W_out) * (T_in - T_out)
        Q_coil = Q_sens + Q_lat - m_lat * CP_W * T_s
        return Q_sens, m_lat, Q_coil, T_s, T_out, W_out, wet

    def _drain_fraction(self, M_c: float) -> float:
        w = self.p.film_ramp * self.p.M_coil_max
        return min(max((M_c - (self.p.M_coil_max - w)) / w, 0.0), 1.0)

    def _evaporation_fraction(self, M_c: float) -> float:
        w = self.p.film_ramp * self.p.M_coil_max
        return min(max(M_c / w, 0.0), 1.0)

    def fan_only(self, f: float, T_in: float, W_in: float, M_c: float) -> CoilState:
        """Fan running with the compressor off: the water film re-evaporates adiabatically."""
        p = self.p
        f = min(max(f, p.f_min), 1.0)
        m = p.m_nom * f
        P_fan_in = p.P_fan_in_nom * f ** 3
        h_in = _h(T_in, W_in)
        T_wb = T_from_h_sat(h_in, T_guess=T_in - 5.0)
        W_out = W_in + (1.0 - self.bypass_factor(f)) * (_W_sat(T_wb) - W_in) \
            * self._evaporation_fraction(M_c)
        T_out = (h_in - H_FG0 * W_out) / (CP_A + CP_V * W_out)
        m_lat = m * (W_in - W_out)                       # <= 0: water added to the air
        Q_sens = m * (CP_A + CP_V * W_out) * (T_in - T_out)
        return CoilState(Q_sens, m_lat, 0.0, P_fan_in, P_fan_in, 0.0, T_wb, T_out, W_out,
                         M_c > 0.0, dM_coil=m_lat, m_drain=0.0)

    def operate(self, s: float, f: float, T_in: float, W_in: float, T_o: float,
                ramp: float = 1.0, M_c: float | None = None) -> CoilState:
        """Operating point for commands (s, f) with return air (T_in, W_in) and outdoor T_o.

        ``ramp`` in [0, 1] scales the delivered capacity after a start; the compressor already
        draws the power of the requested capacity. ``M_c`` is the water held on the coil; when
        it is None the coil is treated as permanently full (steady wet-coil operation).
        """
        p = self.p
        if s <= 0.0:
            if f <= 0.0:
                return OFF_STATE
            return self.fan_only(f, T_in, W_in, 0.0 if M_c is None else M_c)
        s = min(max(s, p.s_min), 1.0)
        f = min(max(f, p.f_min), 1.0)
        Q_req = self.frost_limited(self.capacity(s, T_o), f, T_in, W_in)
        req = self.streams(Q_req, f, T_in, W_in)
        out = req if ramp >= 1.0 else self.streams(Q_req * ramp, f, T_in, W_in)
        Q_sens, m_lat, Q_coil, T_s, T_out, W_out, wet = out
        Q_coil_req, T_s_req = req[2], req[3]
        cop = self._cop(self.eta_full * self.eta_rel(s), Q_coil_req, T_s_req, T_o)
        P_fan_in = p.P_fan_in_nom * f ** 3
        P_el = Q_coil_req / cop + P_fan_in + p.P_fan_out
        g = 1.0 if M_c is None else self._drain_fraction(M_c)
        return CoilState(Q_sens, m_lat, Q_coil, P_el, P_fan_in, cop, T_s, T_out, W_out, wet,
                         dM_coil=m_lat * (1.0 - g), m_drain=m_lat * g)

    def start_ramp(self, t_on: float) -> float:
        """Capacity fraction t_on seconds after a compressor start."""
        return 1.0 - math.exp(-max(t_on, 0.0) / self.p.tau_start)
