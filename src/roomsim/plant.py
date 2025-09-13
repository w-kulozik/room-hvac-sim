"""Continuous-time room plant: two-node RC thermal model plus a moisture balance.

State   x = [T_a, T_m, W]
    T_a  indoor air temperature, degC
    T_m  thermal-mass temperature (walls, floor, furniture), degC
    W    indoor humidity ratio, kg/kg dry air
Control u = [Q_hvac_sens, m_hvac_lat]
    Q_hvac_sens  sensible heat delivered to the air, W (negative when cooling)
    m_hvac_lat   water removed from the air by the unit, kg/s (positive when dehumidifying)

Energy balances
    C_a dT_a/dt = H_am (T_m - T_a) + (H_ao + H_inf)(T_o - T_a) + f_sol,a Q_sol + Q_int,conv + Q_hvac_sens
    C_m dT_m/dt = H_am (T_a - T_m) + H_mo (T_o - T_m) + (1 - f_sol,a) Q_sol + Q_int,rad
Moisture balance (furnishings are assumed not to buffer moisture)
    rho V dW/dt = m_inf (W_o - W) + N_occ g_vap - m_hvac_lat

``fluxes`` returns every term separately and ``rhs`` is assembled from it, so the energy
bookkeeping in the verification tests uses exactly the terms that drive the dynamics.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

from .inputs import Disturbances
from .params import RoomParams

# State and input indices (C_CO2 will be appended to the state on Tuesday)
IDX_TA, IDX_TM, IDX_W = 0, 1, 2
N_STATES = 3
IDX_Q_SENS, IDX_M_LAT = 0, 1
N_INPUTS = 2

U_OFF = np.zeros(N_INPUTS)


class Fluxes(NamedTuple):
    """Heat flows in W and vapour flows in kg/s; sign convention: positive = into the node."""

    # into the air node
    q_mass_to_air: float     # H_am (T_m - T_a); the mass receives the opposite
    q_window: float          # H_ao (T_o - T_a)
    q_infiltration: float    # H_inf (T_o - T_a)
    q_solar_air: float       # f_sol,a * g * A_win * G
    q_internal_conv: float   # convective part of occupant gains
    q_hvac: float            # Q_hvac_sens
    # into the mass node
    q_envelope: float        # H_mo (T_o - T_m)
    q_solar_mass: float      # (1 - f_sol,a) * g * A_win * G
    q_internal_rad: float    # radiant part of occupant gains
    # into the air moisture balance
    m_infiltration: float    # m_inf (W_o - W)
    m_occupants: float       # N_occ * g_vap
    m_hvac: float            # -m_hvac_lat

    @property
    def air_total(self) -> float:
        return (self.q_mass_to_air + self.q_window + self.q_infiltration
                + self.q_solar_air + self.q_internal_conv + self.q_hvac)

    @property
    def mass_total(self) -> float:
        return -self.q_mass_to_air + self.q_envelope + self.q_solar_mass + self.q_internal_rad

    @property
    def boundary_total(self) -> float:
        """Net heat entering the room from outside the two nodes (air-mass exchange cancels)."""
        return (self.q_window + self.q_infiltration + self.q_solar_air + self.q_internal_conv
                + self.q_hvac + self.q_envelope + self.q_solar_mass + self.q_internal_rad)

    @property
    def moisture_total(self) -> float:
        return self.m_infiltration + self.m_occupants + self.m_hvac


def fluxes(t: float, x: np.ndarray, u: np.ndarray, d: Disturbances, p: RoomParams) -> Fluxes:
    """All heat and vapour flows at state x, control u and disturbances d."""
    T_a, T_m, W = x[IDX_TA], x[IDX_TM], x[IDX_W]
    q_sol = p.g_win * p.A_win * d.G_sol
    q_occ = d.n_occ * p.q_person_sens
    return Fluxes(
        q_mass_to_air=p.H_am * (T_m - T_a),
        q_window=p.H_ao * (d.T_o - T_a),
        q_infiltration=p.H_inf * (d.T_o - T_a),
        q_solar_air=p.f_sol_air * q_sol,
        q_internal_conv=p.f_person_conv * q_occ,
        q_hvac=u[IDX_Q_SENS],
        q_envelope=p.H_mo * (d.T_o - T_m),
        q_solar_mass=(1.0 - p.f_sol_air) * q_sol,
        q_internal_rad=(1.0 - p.f_person_conv) * q_occ,
        m_infiltration=p.m_dot_inf * (d.W_o - W),
        m_occupants=d.n_occ * p.g_vap_person,
        m_hvac=-u[IDX_M_LAT],
    )


def rhs(t: float, x: np.ndarray, u: np.ndarray, d: Disturbances, p: RoomParams) -> np.ndarray:
    """dx/dt of the room plant."""
    f = fluxes(t, x, u, d, p)
    dx = np.empty(N_STATES)
    dx[IDX_TA] = f.air_total / p.C_a
    dx[IDX_TM] = f.mass_total / p.C_m
    dx[IDX_W] = f.moisture_total / p.m_air
    return dx


def stored_energy(x: np.ndarray, p: RoomParams) -> np.ndarray:
    """Sensible energy stored relative to 0 degC, J. Works on a single state or on rows of states."""
    x = np.asarray(x)
    return p.C_a * x[..., IDX_TA] + p.C_m * x[..., IDX_TM]
