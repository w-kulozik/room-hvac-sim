"""Container for simulation output with derived quantities and balance checks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import psychro
from .inputs import Disturbances
from .params import RoomParams
from .plant import (IDX_CO2, IDX_PHI_B, IDX_TA, IDX_TM, IDX_W, Fluxes, stored_co2,
                    stored_energy, stored_water)


@dataclass
class SimResult:
    """Time series on the step grid t_0 .. t_N.

    u[k], d[k] and the fluxes at row k are the values in force on [t_k, t_k+1)
    (zero-order hold for the control, step-start sample for piecewise-constant inputs).
    E_in, M_in and CO2_in are the cumulative boundary heat (J), water (kg) and CO2 (m3) since t_0,
    integrated alongside the state with the same method.
    """

    t: np.ndarray
    x: np.ndarray
    u: np.ndarray
    d: Disturbances
    q: dict[str, np.ndarray]
    E_in: np.ndarray
    M_in: np.ndarray
    CO2_in: np.ndarray
    p: RoomParams
    method: str
    dt: float
    E_el: np.ndarray | None = None              # cumulative electrical energy of the unit, J
    cmd: np.ndarray | None = None               # held commands [s, f, V_vent] (with a unit)
    unit: dict[str, np.ndarray] | None = None   # unit operating point at each grid point
    M_coil: np.ndarray | None = None            # water held on the coil, kg
    W_ext: np.ndarray | None = None             # cumulative water from outdoor air and occupants, kg
    M_drain: np.ndarray | None = None           # cumulative condensate drained, kg

    # ---- states and psychrometrics ----------------------------------------------------
    @property
    def T_a(self) -> np.ndarray:
        return self.x[:, IDX_TA]

    @property
    def T_m(self) -> np.ndarray:
        return self.x[:, IDX_TM]

    @property
    def W(self) -> np.ndarray:
        return self.x[:, IDX_W]

    @property
    def phi_b(self) -> np.ndarray:
        return self.x[:, IDX_PHI_B]

    @property
    def CO2(self) -> np.ndarray:
        return self.x[:, IDX_CO2]

    @property
    def RH(self) -> np.ndarray:
        return psychro.RH_from_W(self.T_a, self.W)

    @property
    def T_dew(self) -> np.ndarray:
        return psychro.dew_point_from_W(self.W)

    @property
    def T_dew_out(self) -> np.ndarray:
        return psychro.dew_point_from_W(self.d.W_o)

    # ---- balances ---------------------------------------------------------------------
    @property
    def energy_stored(self) -> np.ndarray:
        """Change in stored sensible energy since t_0, J."""
        E = stored_energy(self.x, self.p)
        return E - E[0]

    @property
    def energy_residual(self) -> np.ndarray:
        """Stored-energy change minus cumulative boundary heat, J (should be ~ round-off)."""
        return self.energy_stored - self.E_in

    @property
    def water_stored(self) -> np.ndarray:
        M = stored_water(self.x, self.p)
        return M - M[0]

    @property
    def water_residual(self) -> np.ndarray:
        return self.water_stored - self.M_in

    @property
    def co2_residual(self) -> np.ndarray:
        C = stored_co2(self.x, self.p)
        return (C - C[0]) - self.CO2_in

    @property
    def system_water_residual(self) -> np.ndarray:
        """Water balance of room air + buffer + coil film: storage change - (in - drained), kg."""
        M_c = np.zeros_like(self.t) if self.M_coil is None else self.M_coil - self.M_coil[0]
        drained = np.zeros_like(self.t) if self.M_drain is None else self.M_drain
        return self.water_stored + M_c - (self.W_ext - drained)

    def boundary_heat_flow(self) -> np.ndarray:
        """Net heat into the room at each grid point, W."""
        return sum(self.q[name] for name in BOUNDARY_HEAT_FIELDS)

    # ---- slicing ----------------------------------------------------------------------
    def window(self, t0: float, t1: float | None = None) -> "SimResult":
        """Sub-result on t0 <= t <= t1, with cumulative integrals restarted at t0."""
        t1 = self.t[-1] if t1 is None else t1
        m = (self.t >= t0 - 1e-9) & (self.t <= t1 + 1e-9)
        return SimResult(
            t=self.t[m], x=self.x[m], u=self.u[m],
            d=Disturbances(*(np.asarray(f)[m] for f in self.d)),
            q={k: v[m] for k, v in self.q.items()},
            E_in=self.E_in[m] - self.E_in[m][0], M_in=self.M_in[m] - self.M_in[m][0],
            CO2_in=self.CO2_in[m] - self.CO2_in[m][0],
            p=self.p, method=self.method, dt=self.dt,
            E_el=None if self.E_el is None else self.E_el[m] - self.E_el[m][0],
            cmd=None if self.cmd is None else self.cmd[m],
            unit=None if self.unit is None else {k: v[m] for k, v in self.unit.items()},
            M_coil=None if self.M_coil is None else self.M_coil[m],
            W_ext=None if self.W_ext is None else self.W_ext[m] - self.W_ext[m][0],
            M_drain=None if self.M_drain is None else self.M_drain[m] - self.M_drain[m][0],
        )


BOUNDARY_HEAT_FIELDS = ("q_window", "q_infiltration", "q_ventilation", "q_solar_air",
                        "q_internal_conv", "q_hvac", "q_envelope", "q_solar_mass", "q_internal_rad")
BOUNDARY_WATER_FIELDS = ("m_infiltration", "m_ventilation", "m_occupants", "m_hvac")
assert set(BOUNDARY_HEAT_FIELDS) | set(BOUNDARY_WATER_FIELDS) < set(Fluxes._fields)
