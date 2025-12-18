"""Continuous-time room plant: two-node RC thermal model, two-node moisture balance and CO2.

State   x = [T_a, T_m, W, phi_b, C]
    T_a    indoor air temperature, degC
    T_m    thermal-mass temperature (walls, floor, furniture), degC
    W      indoor humidity ratio, kg/kg dry air
    phi_b  relative humidity in equilibrium with the moisture-buffer layer, fraction
    C      indoor CO2 concentration, ppm (volume)
Control u = [Q_hvac_sens, m_hvac_lat, V_vent]
    Q_hvac_sens  sensible heat delivered to the air, W (negative when cooling)
    m_hvac_lat   water removed from the air by the unit, kg/s (positive when dehumidifying)
    V_vent       mechanical outdoor-air supply, m3/s (balanced, no heat recovery)

Outdoor air enters by infiltration (H_inf, m_inf) and by ventilation (H_vent = rho cp V_vent,
m_vent = rho V_vent); both carry outdoor temperature, humidity and CO2.

Energy balances
    C_a dT_a/dt = H_am (T_m - T_a) + (H_ao + H_inf + H_vent)(T_o - T_a) + f_sol,a Q_sol
                  + Q_int,conv + Q_hvac_sens
    C_m dT_m/dt = H_am (T_a - T_m) + H_mo (T_o - T_m) + (1 - f_sol,a) Q_sol + Q_int,rad
Moisture balances (EMPD buffer layer on the interior surfaces)
    rho V dW/dt      = (m_inf + m_vent)(W_o - W) + N_occ g_vap - m_hvac_lat + m_buf
    C_b dphi_b/dt    = -m_buf
    m_buf            = G_b (phi_b p_ws(T_m) - p_w(W))      [buffer -> air, kg/s]
The buffer surface sits at the mass temperature, so a warming wall releases moisture.
The heat of sorption (~2.5 MJ/kg) is neglected; it is of order 10 W for daily cycles.
CO2 balance (well mixed; flows as volume of CO2, m3/s)
    V dC/dt * 1e-6   = (V n/3600 + V_vent)(C_o - C) * 1e-6 + N_occ G_CO2

``fluxes`` returns every term separately and ``rhs`` is assembled from it, so the energy
and water bookkeeping in the verification tests uses exactly the terms that drive the dynamics.
"""

from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np

from . import psychro
from .inputs import Disturbances
from .params import RoomParams

# State and input indices
IDX_TA, IDX_TM, IDX_W, IDX_PHI_B, IDX_CO2 = 0, 1, 2, 3, 4
N_STATES = 5
IDX_Q_SENS, IDX_M_LAT, IDX_V_VENT = 0, 1, 2
N_INPUTS = 3

# Size of a small change in each state (K, K, kg/kg, -, ppm), used for finite differences
STATE_SCALE = np.array([1e-3, 1e-3, 1e-6, 1e-4, 1e-1])

U_OFF = np.zeros(N_INPUTS)


class Fluxes(NamedTuple):
    """Heat flows in W and water flows in kg/s; sign convention: positive = into the node."""

    # into the air node
    q_mass_to_air: float     # H_am (T_m - T_a); the mass receives the opposite
    q_window: float          # H_ao (T_o - T_a)
    q_infiltration: float    # H_inf (T_o - T_a)
    q_ventilation: float     # rho cp V_vent (T_o - T_a)
    q_solar_air: float       # f_sol,a * g * A_win * G
    q_internal_conv: float   # convective part of occupant gains
    q_hvac: float            # Q_hvac_sens
    # into the mass node
    q_envelope: float        # H_mo (T_o - T_m)
    q_solar_mass: float      # (1 - f_sol,a) * g * A_win * G
    q_internal_rad: float    # radiant part of occupant gains
    # into the air moisture balance
    m_infiltration: float    # m_inf (W_o - W)
    m_ventilation: float     # rho V_vent (W_o - W)
    m_occupants: float       # N_occ * g_vap
    m_hvac: float            # -m_hvac_lat
    m_buffer: float          # released by the buffer layer; the buffer loses the same amount
    # into the room CO2 balance, m3/s of CO2
    c_outdoor_air: float     # (V n/3600 + V_vent)(C_o - C) * 1e-6
    c_occupants: float       # N_occ * G_CO2

    @property
    def air_total(self) -> float:
        return (self.q_mass_to_air + self.q_window + self.q_infiltration + self.q_ventilation
                + self.q_solar_air + self.q_internal_conv + self.q_hvac)

    @property
    def mass_total(self) -> float:
        return -self.q_mass_to_air + self.q_envelope + self.q_solar_mass + self.q_internal_rad

    @property
    def boundary_total(self) -> float:
        """Net heat entering the room from outside the two nodes (air-mass exchange cancels)."""
        return (self.q_window + self.q_infiltration + self.q_ventilation
                + self.q_solar_air + self.q_internal_conv
                + self.q_hvac + self.q_envelope + self.q_solar_mass + self.q_internal_rad)

    @property
    def moisture_total(self) -> float:
        """Net water flow into the air, kg/s."""
        return (self.m_infiltration + self.m_ventilation + self.m_occupants + self.m_hvac
                + self.m_buffer)

    @property
    def water_boundary_total(self) -> float:
        """Net water entering air + buffer together (air-buffer exchange cancels), kg/s."""
        return self.m_infiltration + self.m_ventilation + self.m_occupants + self.m_hvac

    @property
    def water_external(self) -> float:
        """Water entering from outside the room + unit system (outdoor air, occupants), kg/s."""
        return self.m_infiltration + self.m_ventilation + self.m_occupants

    @property
    def co2_total(self) -> float:
        """Net CO2 entering the room, m3/s."""
        return self.c_outdoor_air + self.c_occupants


def _p_ws(T: float) -> float:
    # scalar fast path of psychro.p_ws, same coefficients
    return psychro._A * math.exp(psychro._B * T / (T + psychro._C))


def buffer_flux(T_m: float, W: float, phi_b: float, p: RoomParams) -> float:
    """Water released by the buffer layer into the air, kg/s."""
    if p.buffer_conductance == 0.0:
        return 0.0
    return p.buffer_conductance * (phi_b * _p_ws(T_m) - W * psychro.P_ATM / (psychro.EPS + W))


def fluxes(t: float, x: np.ndarray, u: np.ndarray, d: Disturbances, p: RoomParams) -> Fluxes:
    """All heat and water flows at state x, control u and disturbances d."""
    T_a, T_m, W, phi_b, C = x[IDX_TA], x[IDX_TM], x[IDX_W], x[IDX_PHI_B], x[IDX_CO2]
    q_sol = p.g_win * p.A_win * d.G_sol
    q_occ = d.n_occ * p.q_person_sens
    V_vent = u[IDX_V_VENT]
    m_vent = p.rho_air * V_vent
    return Fluxes(
        q_mass_to_air=p.H_am * (T_m - T_a),
        q_window=p.H_ao * (d.T_o - T_a),
        q_infiltration=p.H_inf * (d.T_o - T_a),
        q_ventilation=m_vent * p.cp_air * (d.T_o - T_a),
        q_solar_air=p.f_sol_air * q_sol,
        q_internal_conv=p.f_person_conv * q_occ,
        q_hvac=u[IDX_Q_SENS],
        q_envelope=p.H_mo * (d.T_o - T_m),
        q_solar_mass=(1.0 - p.f_sol_air) * q_sol,
        q_internal_rad=(1.0 - p.f_person_conv) * q_occ,
        m_infiltration=p.m_dot_inf * (d.W_o - W),
        m_ventilation=m_vent * (d.W_o - W),
        m_occupants=d.n_occ * p.g_vap_person,
        m_hvac=-u[IDX_M_LAT],
        m_buffer=buffer_flux(T_m, W, phi_b, p),
        c_outdoor_air=(p.V_dot_inf + V_vent) * (d.C_o - C) * 1e-6,
        c_occupants=d.n_occ * p.co2_person,
    )


def derivatives(f: Fluxes, p: RoomParams) -> np.ndarray:
    """dx/dt from an already evaluated set of fluxes."""
    dx = np.empty(N_STATES)
    dx[IDX_TA] = f.air_total / p.C_a
    dx[IDX_TM] = f.mass_total / p.C_m
    dx[IDX_W] = f.moisture_total / p.m_air
    dx[IDX_PHI_B] = -f.m_buffer / p.buffer_capacity if p.buffer_capacity > 0 else 0.0
    dx[IDX_CO2] = f.co2_total / p.volume * 1e6
    return dx


def rhs(t: float, x: np.ndarray, u: np.ndarray, d: Disturbances, p: RoomParams) -> np.ndarray:
    """dx/dt of the room plant."""
    return derivatives(fluxes(t, x, u, d, p), p)


def jacobian(t: float, x: np.ndarray, u: np.ndarray, d: Disturbances, p: RoomParams) -> np.ndarray:
    """Jacobian d(rhs)/dx by central differences (used for stability limits and, later, MPC)."""
    x = np.asarray(x, dtype=float)
    J = np.empty((N_STATES, N_STATES))
    for j in range(N_STATES):
        e = np.zeros(N_STATES)
        e[j] = STATE_SCALE[j]
        J[:, j] = (rhs(t, x + e, u, d, p) - rhs(t, x - e, u, d, p)) / (2 * STATE_SCALE[j])
    return J


def measure(x: np.ndarray, p: RoomParams) -> np.ndarray:
    """Sensor outputs available to a controller: [T_a, RH_a, C]. T_m and phi_b are not measured."""
    T_a, W = x[IDX_TA], x[IDX_W]
    return np.array([T_a, float(psychro.RH_from_W(T_a, W)), x[IDX_CO2]])


IDX_Y_TA, IDX_Y_RH, IDX_Y_CO2 = 0, 1, 2
IDX_Y_TCOIL = 3     # appended by the simulator when an air-conditioning unit is present


def equilibrium_state(T_a: float, T_m: float, W: float, C: float = 420.0) -> np.ndarray:
    """State with the buffer layer in equilibrium with the air (no buffer flux)."""
    phi_b = float(psychro.pw_from_W(W) / psychro.p_ws(T_m))
    return np.array([T_a, T_m, W, phi_b, C])


def stored_energy(x: np.ndarray, p: RoomParams) -> np.ndarray:
    """Sensible energy stored relative to 0 degC, J. Works on a single state or on rows of states."""
    x = np.asarray(x)
    return p.C_a * x[..., IDX_TA] + p.C_m * x[..., IDX_TM]


def stored_water(x: np.ndarray, p: RoomParams) -> np.ndarray:
    """Water held in the room air plus the buffer layer, kg (buffer relative to phi_b = 0)."""
    x = np.asarray(x)
    return p.m_air * x[..., IDX_W] + p.buffer_capacity * x[..., IDX_PHI_B]


def stored_co2(x: np.ndarray, p: RoomParams) -> np.ndarray:
    """CO2 held in the room air, m3."""
    return p.volume * np.asarray(x)[..., IDX_CO2] * 1e-6
