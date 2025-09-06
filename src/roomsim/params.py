"""Room parameters for the two-node RC thermal model and the moisture balance.

Defaults describe a 5 x 4 x 2.7 m office for two people with one 3 m2 window.
Derived quantities (capacities, conductances, time constants) are properties so that
changing a primary parameter with ``dataclasses.replace`` keeps everything consistent.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class RoomParams:
    # Geometry and air properties
    volume: float = 5.0 * 4.0 * 2.7          # V, m3
    rho_air: float = 1.2                     # kg/m3
    cp_air: float = 1006.0                   # J/(kg K)
    furniture_factor: float = 5.0            # multiplies air capacity to include light furnishings

    # Thermal mass node (walls, floor, heavy furniture)
    C_m: float = 10.0e6                      # J/K

    # Conductances
    h_i: float = 3.0                         # interior convective coefficient, W/(m2 K)
    A_int: float = 60.0                      # interior surface area exchanging with air, m2
    H_ao: float = 4.0                        # air <-> outdoor through the window, W/K (3 m2 x 1.3)
    H_mo: float = 25.0                       # mass <-> outdoor through opaque envelope, W/K
    ach: float = 0.5                         # infiltration, air changes per hour, 1/h

    # Window and solar split
    A_win: float = 3.0                       # m2
    g_win: float = 0.6                       # solar heat gain coefficient, -
    f_sol_air: float = 0.3                   # fraction of transmitted solar going directly to air

    # Occupants
    q_person_sens: float = 70.0              # sensible gain per person, W
    f_person_conv: float = 0.5               # convective fraction of sensible gain
    g_vap_person_kg_h: float = 0.05          # vapour released per person, kg/h

    # ---- derived quantities -------------------------------------------------------
    @property
    def m_air(self) -> float:
        """Mass of room air, kg."""
        return self.rho_air * self.volume

    @property
    def C_a(self) -> float:
        """Air-node heat capacity including light furnishings, J/K."""
        return self.rho_air * self.cp_air * self.volume * self.furniture_factor

    @property
    def H_inf(self) -> float:
        """Infiltration conductance rho*cp*V*n/3600, W/K."""
        return self.rho_air * self.cp_air * self.volume * self.ach / 3600.0

    @property
    def H_am(self) -> float:
        """Air <-> mass conductance 1/R_am = h_i*A_int, W/K."""
        return self.h_i * self.A_int

    @property
    def m_dot_inf(self) -> float:
        """Infiltration mass flow of dry air, kg/s."""
        return self.m_air * self.ach / 3600.0

    @property
    def g_vap_person(self) -> float:
        """Vapour released per person, kg/s."""
        return self.g_vap_person_kg_h / 3600.0

    @property
    def tau_air(self) -> float:
        """Air time constant with the mass node held fixed, s."""
        return self.C_a / (self.H_am + self.H_ao + self.H_inf)

    @property
    def tau_mass(self) -> float:
        """Mass time constant with the air node quasi-static, s.

        Seen from the mass, the air node is H_am in series with (H_ao + H_inf) to outdoors.
        """
        H_air_out = self.H_ao + self.H_inf
        H_series = self.H_am * H_air_out / (self.H_am + H_air_out) if H_air_out > 0 else 0.0
        return self.C_m / (self.H_mo + H_series)

    @property
    def tau_moisture(self) -> float:
        """Moisture time constant 3600/n, s."""
        return 3600.0 / self.ach
