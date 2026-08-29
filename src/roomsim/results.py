"""Container for simulation output with derived quantities and balance checks."""

from __future__ import annotations

from dataclasses import dataclass
from typing import cast

import numpy as np

from . import psychro
from ._typing import BoolArray, FloatArray
from .inputs import Disturbances
from .params import RoomParams
from .plant import (
    IDX_CO2,
    IDX_PHI_B,
    IDX_TA,
    IDX_TM,
    IDX_W,
    Fluxes,
    stored_co2,
    stored_energy,
    stored_water,
)

# Unit operating-point quantities logged at every grid point
UNIT_FIELDS = ("Q_sens", "m_lat", "Q_coil", "P_el", "P_fan_in", "COP", "T_s", "T_out",
               "wet_fraction", "m_drain")


@dataclass
class UnitLog:
    """Record of the air-conditioning unit, present when a unit is simulated.

    All arrays have one row per grid point of the owning `SimResult`.

    Attributes
    ----------
    cmd : FloatArray
        Held commands, one row per grid point, columns [compressor speed s (-), indoor fan
        speed f (-), ventilation flow V_vent (m3/s)].
    op : dict of str to FloatArray
        Operating point of the unit, keyed by the names in `UNIT_FIELDS`.
    M_coil : FloatArray
        Water held on the coil, kg.
    E_el : FloatArray
        Cumulative electrical energy since t_0, J.
    M_drain : FloatArray
        Cumulative condensate drained since t_0, kg.
    E_cool : FloatArray
        Cumulative net cooling delivered to the room since t_0 (sensible plus latent, minus
        the indoor-fan heat; ``roomsim.ac_unit.net_cooling``), J.
    running_before : bool
        True when the compressor was running (s > 0) on the step just before t_0, i.e. the
        record is a window of a longer run that starts inside a compressor cycle. False for
        the start of a simulation.
    """

    cmd: FloatArray              # held commands [s, f, V_vent]
    op: dict[str, FloatArray]    # operating point, keys UNIT_FIELDS
    M_coil: FloatArray           # water held on the coil, kg
    E_el: FloatArray             # cumulative electrical energy since t_0, J
    M_drain: FloatArray          # cumulative condensate drained since t_0, kg
    E_cool: FloatArray           # cumulative net cooling delivered since t_0, J
    running_before: bool = False  # compressor running on the step before t_0

    def window(self, m: BoolArray) -> UnitLog:
        """Return the contiguous rows selected by a mask, with the cumulative quantities restarted.

        Parameters
        ----------
        m : BoolArray
            Boolean mask over the grid points, selecting one contiguous, non-empty block.

        Returns
        -------
        UnitLog
            Selected rows; `E_el`, `M_drain` and `E_cool` start from zero at the first selected
            row, and `running_before` holds the compressor state of the row before it (inherited
            when the block starts at the first row).
        """
        i0 = int(np.argmax(m))
        running_before = (bool(self.cmd[i0 - 1, 0] > 0) if i0 > 0  # IDX_CMD_S of simulate
                          else self.running_before)
        return UnitLog(cmd=self.cmd[m], op={k: v[m] for k, v in self.op.items()},
                       M_coil=self.M_coil[m], E_el=self.E_el[m] - self.E_el[m][0],
                       M_drain=self.M_drain[m] - self.M_drain[m][0],
                       E_cool=self.E_cool[m] - self.E_cool[m][0],
                       running_before=running_before)


@dataclass
class SimResult:
    """Time series on the step grid t_0 .. t_N.

    u[k], d[k] and the fluxes at row k are the values in force on [t_k, t_k+1)
    (zero-order hold for the control, step-start sample for piecewise-constant inputs).
    The cumulative quantities since t_0 are integrated alongside the state with the same method:
    boundary heat E_in (J), room water M_in (kg), CO2 (m3), and water from outdoor air and
    occupants W_ext (kg).

    Attributes
    ----------
    t : FloatArray
        Grid times t_0 .. t_N, s.
    x : FloatArray
        States, one row per grid point, columns indexed by the ``IDX_*`` constants of `plant`.
    u : FloatArray
        Plant inputs held on each step, one row per grid point, columns [Q_hvac_sens (W),
        m_hvac_lat (kg/s), V_vent (m3/s)].
    d : Disturbances
        Disturbances sampled at the grid times.
    q : dict of str to FloatArray
        Heat flows (W) and water flows (kg/s) at the grid times, keyed by `plant.Fluxes` names.
    E_in : FloatArray
        Cumulative boundary heat into the room since t_0, J.
    M_in : FloatArray
        Cumulative water into the room since t_0, kg.
    CO2_in : FloatArray
        Cumulative CO2 into the room since t_0, m3.
    W_ext : FloatArray
        Cumulative water from outdoor air and occupants since t_0, kg.
    p : RoomParams
        Room parameters of the run.
    method : str
        Name of the integration method.
    dt : float
        Time step, s.
    unit : UnitLog or None
        Record of the air-conditioning unit, or None when no unit was simulated.
    """

    t: FloatArray
    x: FloatArray
    u: FloatArray
    d: Disturbances
    q: dict[str, FloatArray]
    E_in: FloatArray
    M_in: FloatArray
    CO2_in: FloatArray
    W_ext: FloatArray
    p: RoomParams
    method: str
    dt: float
    unit: UnitLog | None = None

    def require_unit(self) -> UnitLog:
        """Return the unit record of a run simulated with an air-conditioning unit.

        Returns
        -------
        UnitLog
            The `unit` attribute.

        Raises
        ------
        ValueError
            If the run was simulated without a unit.
        """
        if self.unit is None:
            raise ValueError("the result has no air-conditioning unit record")
        return self.unit

    # ---- states and psychrometrics ----------------------------------------------------
    @property
    def T_a(self) -> FloatArray:
        """Indoor air temperature, degC."""
        return self.x[:, IDX_TA]

    @property
    def T_m(self) -> FloatArray:
        """Thermal-mass temperature, degC."""
        return self.x[:, IDX_TM]

    @property
    def W(self) -> FloatArray:
        """Indoor humidity ratio, kg/kg."""
        return self.x[:, IDX_W]

    @property
    def phi_b(self) -> FloatArray:
        """Relative humidity of the wall buffer layer, fraction."""
        return self.x[:, IDX_PHI_B]

    @property
    def CO2(self) -> FloatArray:
        """Indoor CO2 concentration, ppm."""
        return self.x[:, IDX_CO2]

    @property
    def RH(self) -> FloatArray:
        """Indoor relative humidity, fraction."""
        return psychro.RH_from_W(self.T_a, self.W)

    @property
    def T_dew(self) -> FloatArray:
        """Indoor dew-point temperature, degC."""
        return psychro.dew_point_from_W(self.W)

    @property
    def T_dew_out(self) -> FloatArray:
        """Outdoor dew-point temperature, degC."""
        return psychro.dew_point_from_W(self.d.W_o)

    # ---- balances ---------------------------------------------------------------------
    @property
    def energy_stored(self) -> FloatArray:
        """Change in stored sensible energy since t_0, J."""
        E = stored_energy(self.x, self.p)
        return cast(FloatArray, E - E[0])

    @property
    def energy_residual(self) -> FloatArray:
        """Stored-energy change minus cumulative boundary heat, J (should be ~ round-off)."""
        return self.energy_stored - self.E_in

    @property
    def water_stored(self) -> FloatArray:
        """Change in water held by the room air and the buffer layer since t_0, kg."""
        M = stored_water(self.x, self.p)
        return cast(FloatArray, M - M[0])

    @property
    def water_residual(self) -> FloatArray:
        """Stored-water change minus cumulative water into the room, kg."""
        return self.water_stored - self.M_in

    @property
    def co2_residual(self) -> FloatArray:
        """Stored-CO2 change minus cumulative CO2 into the room, m3."""
        C = stored_co2(self.x, self.p)
        return cast(FloatArray, (C - C[0]) - self.CO2_in)

    @property
    def system_water_residual(self) -> FloatArray:
        """Water balance of room air + buffer + coil film: storage change - (in - drained), kg."""
        if self.unit is None:
            return self.water_stored - self.W_ext
        M_c = self.unit.M_coil - self.unit.M_coil[0]
        return cast(FloatArray, self.water_stored + M_c - (self.W_ext - self.unit.M_drain))

    def boundary_heat_flow(self) -> FloatArray:
        """Return the net heat into the room at each grid point.

        Returns
        -------
        FloatArray
            Sum of the flows named in `BOUNDARY_HEAT_FIELDS`, W.
        """
        return cast(FloatArray, sum(self.q[name] for name in BOUNDARY_HEAT_FIELDS))

    # ---- slicing ----------------------------------------------------------------------
    def window(self, t0: float, t1: float | None = None) -> SimResult:
        """Return the sub-result on t0 <= t <= t1, with cumulative integrals restarted at t0.

        Parameters
        ----------
        t0 : float
            Start time, s.
        t1 : float or None, optional
            End time, s. None selects up to the last grid time.

        Returns
        -------
        SimResult
            Rows with t0 - 1e-9 <= t <= t1 + 1e-9; cumulative quantities start from zero.

        Raises
        ------
        ValueError
            If t1 <= t0, if [t0, t1] does not lie within [t_0, t_N], or if fewer than two
            grid points fall in the window.
        """
        t1 = float(self.t[-1]) if t1 is None else t1
        tA, tB = float(self.t[0]), float(self.t[-1])
        if not t1 > t0:
            raise ValueError(f"window end t1 = {t1} s must be greater than its start t0 = {t0} s")
        if t0 < tA - 1e-9 or t1 > tB + 1e-9:
            raise ValueError(f"window [t0, t1] = [{t0}, {t1}] s must lie within the run "
                             f"[{tA}, {tB}] s")
        m = (self.t >= t0 - 1e-9) & (self.t <= t1 + 1e-9)
        if m.sum() < 2:
            raise ValueError(f"window [t0, t1] = [{t0}, {t1}] s holds {int(m.sum())} grid "
                             f"point(s); it needs at least two (dt = {self.dt} s)")
        return SimResult(
            t=self.t[m], x=self.x[m], u=self.u[m],
            d=Disturbances(*(np.asarray(f)[m] for f in self.d)),
            q={k: v[m] for k, v in self.q.items()},
            E_in=self.E_in[m] - self.E_in[m][0], M_in=self.M_in[m] - self.M_in[m][0],
            CO2_in=self.CO2_in[m] - self.CO2_in[m][0], W_ext=self.W_ext[m] - self.W_ext[m][0],
            p=self.p, method=self.method, dt=self.dt,
            unit=None if self.unit is None else self.unit.window(m),
        )


BOUNDARY_HEAT_FIELDS = ("q_window", "q_infiltration", "q_ventilation", "q_solar_air",
                        "q_internal_conv", "q_hvac", "q_envelope", "q_solar_mass", "q_internal_rad")
BOUNDARY_WATER_FIELDS = ("m_infiltration", "m_ventilation", "m_occupants", "m_hvac")
assert set(BOUNDARY_HEAT_FIELDS) | set(BOUNDARY_WATER_FIELDS) < set(Fluxes._fields)
