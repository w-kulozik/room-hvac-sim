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
from typing import NamedTuple, cast

import numpy as np

from . import psychro
from ._typing import FloatArray
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
    """Heat flows in W and water flows in kg/s; sign convention: positive = into the node.

    Attributes
    ----------
    q_mass_to_air : float
        H_am (T_m - T_a), W; the mass node receives the opposite.
    q_window : float
        H_ao (T_o - T_a), W.
    q_infiltration : float
        H_inf (T_o - T_a), W.
    q_ventilation : float
        rho cp V_vent (T_o - T_a), W.
    q_solar_air : float
        Solar gain delivered to the air, f_sol,a g A_win G, W.
    q_internal_conv : float
        Convective part of the occupant gains, W.
    q_hvac : float
        Sensible heat delivered by the HVAC input, W.
    q_envelope : float
        H_mo (T_o - T_m), into the mass node, W.
    q_solar_mass : float
        Solar gain delivered to the mass, (1 - f_sol,a) g A_win G, W.
    q_internal_rad : float
        Radiant part of the occupant gains, into the mass node, W.
    m_infiltration : float
        m_inf (W_o - W), kg/s.
    m_ventilation : float
        rho V_vent (W_o - W), kg/s.
    m_occupants : float
        Vapour released by occupants, N_occ g_vap, kg/s.
    m_hvac : float
        Water added by the HVAC input, -m_hvac_lat, kg/s.
    m_buffer : float
        Water released by the buffer layer into the air, kg/s; the buffer loses the same amount.
    c_outdoor_air : float
        (V n/3600 + V_vent)(C_o - C) 1e-6, CO2 volume flow, m3/s.
    c_occupants : float
        CO2 released by occupants, N_occ G_CO2, m3/s.
    """

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
        """Net heat flow into the air node, W."""
        return (self.q_mass_to_air + self.q_window + self.q_infiltration + self.q_ventilation
                + self.q_solar_air + self.q_internal_conv + self.q_hvac)

    @property
    def mass_total(self) -> float:
        """Net heat flow into the mass node, W."""
        return -self.q_mass_to_air + self.q_envelope + self.q_solar_mass + self.q_internal_rad

    @property
    def boundary_total(self) -> float:
        """Net heat entering the room from outside the two nodes (air-mass exchange cancels), W."""
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
    """Compute the water released by the buffer layer into the air.

    Parameters
    ----------
    T_m : float
        Thermal-mass (buffer surface) temperature, degC.
    W : float
        Indoor humidity ratio, kg/kg dry air.
    phi_b : float
        Relative humidity in equilibrium with the buffer layer, fraction.
    p : RoomParams
        Room parameters.

    Returns
    -------
    float
        Water flow from the buffer to the air, kg/s (negative when the buffer absorbs).
    """
    if p.buffer_conductance == 0.0:
        return 0.0
    return p.buffer_conductance * (phi_b * _p_ws(T_m) - W * psychro.P_ATM / (psychro.EPS + W))


def fluxes(t: float, x: FloatArray, u: FloatArray, d: Disturbances, p: RoomParams) -> Fluxes:
    """Compute all heat and water flows at state x, control u and disturbances d.

    Parameters
    ----------
    t : float
        Time, s (not used; the plant is time-invariant given d).
    x : ndarray, shape (5,)
        Room state [T_a degC, T_m degC, W kg/kg, phi_b -, C ppm].
    u : ndarray, shape (3,)
        Plant input [Q_hvac_sens W, m_hvac_lat kg/s, V_vent m3/s].
    d : Disturbances
        Disturbances sampled at a single time (scalar fields).
    p : RoomParams
        Room parameters.

    Returns
    -------
    Fluxes
        Every heat flow (W), water flow (kg/s) and CO2 flow (m3/s) of the balances.
    """
    T_a, T_m, W, phi_b, C = x[IDX_TA], x[IDX_TM], x[IDX_W], x[IDX_PHI_B], x[IDX_CO2]
    # d is sampled at a single time here, so its fields are scalars
    G_sol, n_occ = cast(float, d.G_sol), cast(float, d.n_occ)
    q_sol = p.g_win * p.A_win * G_sol
    q_occ = n_occ * p.q_person_sens
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
        m_occupants=n_occ * p.g_vap_person,
        m_hvac=-u[IDX_M_LAT],
        m_buffer=buffer_flux(T_m, W, phi_b, p),
        c_outdoor_air=(p.V_dot_inf + V_vent) * (d.C_o - C) * 1e-6,
        c_occupants=n_occ * p.co2_person,
    )


def derivatives(f: Fluxes, p: RoomParams) -> FloatArray:
    """Compute dx/dt from an already evaluated set of fluxes.

    Parameters
    ----------
    f : Fluxes
        Flows at the current state.
    p : RoomParams
        Room parameters.

    Returns
    -------
    ndarray, shape (5,)
        State derivative [K/s, K/s, kg/(kg s), 1/s, ppm/s].
    """
    dx = np.empty(N_STATES)
    dx[IDX_TA] = f.air_total / p.C_a
    dx[IDX_TM] = f.mass_total / p.C_m
    dx[IDX_W] = f.moisture_total / p.m_air
    dx[IDX_PHI_B] = -f.m_buffer / p.buffer_capacity if p.buffer_capacity > 0 else 0.0
    dx[IDX_CO2] = f.co2_total / p.volume * 1e6
    return dx


def rhs(t: float, x: FloatArray, u: FloatArray, d: Disturbances, p: RoomParams) -> FloatArray:
    """Compute dx/dt of the room plant.

    Parameters
    ----------
    t : float
        Time, s.
    x : ndarray, shape (5,)
        Room state [T_a degC, T_m degC, W kg/kg, phi_b -, C ppm].
    u : ndarray, shape (3,)
        Plant input [Q_hvac_sens W, m_hvac_lat kg/s, V_vent m3/s].
    d : Disturbances
        Disturbances sampled at time t.
    p : RoomParams
        Room parameters.

    Returns
    -------
    ndarray, shape (5,)
        State derivative [K/s, K/s, kg/(kg s), 1/s, ppm/s].
    """
    return derivatives(fluxes(t, x, u, d, p), p)


def jacobian(t: float, x: FloatArray, u: FloatArray, d: Disturbances,
             p: RoomParams) -> FloatArray:
    """Compute the Jacobian d(rhs)/dx by central differences.

    Used for stability limits and by the predictive controller. The step for state j is
    ``STATE_SCALE[j]``.

    Parameters
    ----------
    t : float
        Time, s.
    x : ndarray, shape (5,)
        Room state [T_a degC, T_m degC, W kg/kg, phi_b -, C ppm].
    u : ndarray, shape (3,)
        Plant input [Q_hvac_sens W, m_hvac_lat kg/s, V_vent m3/s].
    d : Disturbances
        Disturbances sampled at time t.
    p : RoomParams
        Room parameters.

    Returns
    -------
    ndarray, shape (5, 5)
        Jacobian, entry (i, j) = d(dx_i/dt)/dx_j, in the units of state i per s per unit of
        state j.
    """
    x = np.asarray(x, dtype=float)
    J = np.empty((N_STATES, N_STATES))
    for j in range(N_STATES):
        e = np.zeros(N_STATES)
        e[j] = STATE_SCALE[j]
        J[:, j] = (rhs(t, x + e, u, d, p) - rhs(t, x - e, u, d, p)) / (2 * STATE_SCALE[j])
    return J


def measure(x: FloatArray, p: RoomParams) -> FloatArray:
    """Return the sensor outputs available to a controller; T_m and phi_b are not measured.

    Parameters
    ----------
    x : ndarray, shape (5,)
        Room state [T_a degC, T_m degC, W kg/kg, phi_b -, C ppm].
    p : RoomParams
        Room parameters (not used).

    Returns
    -------
    ndarray, shape (3,)
        Measurement [T_a degC, RH_a fraction, C ppm].
    """
    T_a, W = x[IDX_TA], x[IDX_W]
    return np.array([T_a, float(psychro.RH_from_W(T_a, W)), x[IDX_CO2]])


IDX_Y_TA, IDX_Y_RH, IDX_Y_CO2 = 0, 1, 2
IDX_Y_TCOIL = 3     # appended by the simulator when an air-conditioning unit is present


def equilibrium_state(T_a: float, T_m: float, W: float, C: float = 420.0) -> FloatArray:
    """Build a state with the buffer layer in equilibrium with the air (no buffer flux).

    Parameters
    ----------
    T_a : float
        Indoor air temperature, degC.
    T_m : float
        Thermal-mass temperature, degC.
    W : float
        Indoor humidity ratio, kg/kg dry air.
    C : float, optional
        Indoor CO2 concentration, ppm.

    Returns
    -------
    ndarray, shape (5,)
        State [T_a, T_m, W, phi_b, C] with phi_b = p_w(W) / p_ws(T_m).
    """
    phi_b = float(psychro.pw_from_W(W) / psychro.p_ws(T_m))
    return np.array([T_a, T_m, W, phi_b, C])


def stored_energy(x: FloatArray, p: RoomParams) -> FloatArray:
    """Compute the sensible energy stored relative to 0 degC.

    Parameters
    ----------
    x : ndarray, shape (..., 5)
        A single state or rows of states.
    p : RoomParams
        Room parameters.

    Returns
    -------
    ndarray, shape (...)
        C_a T_a + C_m T_m, J.
    """
    x = np.asarray(x)
    return p.C_a * x[..., IDX_TA] + p.C_m * x[..., IDX_TM]


def stored_water(x: FloatArray, p: RoomParams) -> FloatArray:
    """Compute the water held in the room air plus the buffer layer.

    Parameters
    ----------
    x : ndarray, shape (..., 5)
        A single state or rows of states.
    p : RoomParams
        Room parameters.

    Returns
    -------
    ndarray, shape (...)
        Water mass, kg (buffer relative to phi_b = 0).
    """
    x = np.asarray(x)
    return p.m_air * x[..., IDX_W] + p.buffer_capacity * x[..., IDX_PHI_B]


def stored_co2(x: FloatArray, p: RoomParams) -> FloatArray:
    """Compute the CO2 held in the room air.

    Parameters
    ----------
    x : ndarray, shape (..., 5)
        A single state or rows of states.
    p : RoomParams
        Room parameters.

    Returns
    -------
    ndarray, shape (...)
        CO2 volume, m3.
    """
    return p.volume * np.asarray(x)[..., IDX_CO2] * 1e-6
