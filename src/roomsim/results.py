"""Container for simulation output with derived quantities and balance checks."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import psychro
from .inputs import Disturbances
from .params import RoomParams
from .plant import IDX_PHI_B, IDX_TA, IDX_TM, IDX_W, Fluxes, stored_energy, stored_water


@dataclass
class SimResult:
    """Time series on the step grid t_0 .. t_N.

    u[k], d[k] and the fluxes at row k are the values in force on [t_k, t_k+1)
    (zero-order hold for the control, step-start sample for piecewise-constant inputs).
    E_in and M_in are the cumulative boundary heat (J) and water (kg) since t_0,
    integrated alongside the state with the same method.
    """

    t: np.ndarray
    x: np.ndarray
    u: np.ndarray
    d: Disturbances
    q: dict[str, np.ndarray]
    E_in: np.ndarray
    M_in: np.ndarray
    p: RoomParams
    method: str
    dt: float

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
            p=self.p, method=self.method, dt=self.dt,
        )


BOUNDARY_HEAT_FIELDS = ("q_window", "q_infiltration", "q_solar_air", "q_internal_conv",
                        "q_hvac", "q_envelope", "q_solar_mass", "q_internal_rad")
BOUNDARY_WATER_FIELDS = ("m_infiltration", "m_occupants", "m_hvac")
assert set(BOUNDARY_HEAT_FIELDS) | set(BOUNDARY_WATER_FIELDS) < set(Fluxes._fields)
