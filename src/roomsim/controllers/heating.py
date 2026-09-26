"""Heating control of the heat pump: the PI in heating mode and the defrost supervisor.

Decisions D52-D53. The defrost supervisor wraps a room controller. While it defrosts, and
during the drip time after a defrost, it overrides the command and does not call the room
controller, whose state is held, as in a real unit. Its measurement is that of
``HeatPump.sensors``: y = [T_a, RH, C, T_coil,out, T_o].

    T (timed)   start a defrost after t_timed of accumulated compressor running with the coil
                below T_cold since the end of the last defrost.
    D (demand)  start a defrost when the coil reads dT_demand below a clean coil at the same
                compressor speed, return air and outdoor temperature (dry outdoor air), after at
                least t_min of running since the last defrost; and in any case after t_max of
                running with the coil below T_cold.
Both end a defrost when the coil reads T_end or after t_end_max, then stop the compressor for
t_drip before handing back to the room controller.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from .._typing import FloatArray
from ..heat_pump import HeatPump
from ..params import RoomParams
from ..plant import IDX_Y_TA
from .base import CompressorGuard, Controller, SetPoint, VentilationSchedule
from .lowload import SigmaDeltaLowLoad, ripple_T_ref
from .pi import DiscretePI, simc_pi
from .pi_unit import PIUnitController

IDX_Y_TCOIL_OUT, IDX_Y_TOUT = 3, 4      # heat-pump sensors appended to the room measurement
T_O_TUNING = 2.0                        # outdoor temperature of the tuning point, degC (D53)


def make_heating_pi(p: RoomParams, hp: HeatPump, Ts: float = 60.0, T_set: SetPoint = 20.0,
                    min_on: float = 180.0, min_off: float = 180.0, ripple_K: float = 0.5,
                    tau_c: float = 900.0, **kwargs: Any) -> PIUnitController:
    """Build the PI with sigma-delta low-load handling for heating (D53).

    The SIMC rule of ``room_pi_tuning`` is applied with the heating capacity at full speed and
    2 degC outdoors as the process gain; the sigma-delta time scale gives an air-temperature
    ripple of ``ripple_K`` with the heating output at s_min.

    Parameters
    ----------
    p : RoomParams
        Room parameters used for tuning.
    hp : HeatPump
        Heat pump.
    Ts : float, optional
        Sampling period, s.
    T_set : float or callable, optional
        Set point, degC, or a function of time (s) returning it.
    min_on, min_off : float, optional
        Minimum run and stop times of the compressor, s.
    ripple_K : float, optional
        Target air-temperature ripple at low load, peak to peak, K.
    tau_c : float, optional
        Desired closed-loop time constant, s.
    **kwargs
        Passed to ``PIUnitController`` (``fan``, ``fan_mode``, ``ventilation``).

    Returns
    -------
    PIUnitController
        PI in heating mode, demand up to the heat pump's s_max.
    """
    K = hp.capacity(1.0, T_O_TUNING) / (p.H_am + p.H_ao + p.H_inf)
    Kp, Ti = simc_pi(K, p.tau_air, Ts / 2.0, tau_c)
    s_min = hp.p.s_min
    Q_per_demand = hp.capacity(s_min, T_O_TUNING) / s_min
    T_ref = ripple_T_ref(ripple_K, p.C_a, Q_per_demand, s_min)
    lowload = SigmaDeltaLowLoad(s_min, CompressorGuard(min_on, min_off), Ts, T_ref,
                                s_max=hp.hp.s_max)
    return PIUnitController(DiscretePI(Kp, Ti, Ts), lowload, T_set=T_set, heating=True,
                            **kwargs)


class DefrostControl:
    """Defrost supervisor around a heating controller (D52).

    Parameters
    ----------
    controller : Controller
        Room controller in heating mode.
    hp : HeatPump
        Heat pump; its clean-coil characteristic is the reference of the demand strategy.
    strategy : str
        'timed' or 'demand'.
    ventilation : VentilationSchedule, optional
        Outdoor-air flow schedule, held during defrost and drip time.
    Ts : float, optional
        Sampling period, s.
    t_timed : float, optional
        Running time with a cold coil that starts a timed defrost, s.
    T_cold : float, optional
        Coil temperature below which running time counts as frosting, degC.
    dT_demand : float, optional
        Coil temperature deficit against a clean coil that starts a demand defrost, K.
    t_min, t_max : float, optional
        Least running time before a demand defrost, and the cold running time that forces one,
        s.
    T_end : float, optional
        Coil temperature that ends a defrost, degC.
    t_end_max : float, optional
        Longest defrost, s.
    t_drip : float, optional
        Compressor stop after a defrost, s.

    Attributes
    ----------
    defrosts : int
        Number of defrosts started.
    mode : str
        'normal', 'defrost' or 'drip'.

    Raises
    ------
    ValueError
        If strategy is not 'timed' or 'demand'.
    """

    def __init__(self, controller: Controller, hp: HeatPump, strategy: str,
                 ventilation: VentilationSchedule = VentilationSchedule(),  # noqa: B008  frozen
                 Ts: float = 60.0, t_timed: float = 2400.0, T_cold: float = -1.0,
                 dT_demand: float = 2.0, t_min: float = 1200.0, t_max: float = 10800.0,
                 T_end: float = 8.0, t_end_max: float = 600.0, t_drip: float = 60.0) -> None:
        if strategy not in ("timed", "demand"):
            raise ValueError(f"strategy = {strategy!r} is invalid; it must be 'timed' or 'demand'")
        self.controller, self.hp, self.strategy, self.ventilation = (controller, hp, strategy,
                                                                     ventilation)
        self.Ts = Ts
        self.t_timed, self.T_cold, self.dT_demand = t_timed, T_cold, dT_demand
        self.t_min, self.t_max = t_min, t_max
        self.T_end, self.t_end_max, self.t_drip = T_end, t_end_max, t_drip
        self.mode = "normal"
        self.t_mode = -np.inf          # start of the current defrost or drip time, s
        self.run_since = 0.0           # running time since the last defrost, s
        self.run_cold = 0.0            # running time with a cold coil since the last defrost, s
        self.s_prev = 0.0              # compressor command of the previous sample, -
        self.defrosts = 0

    def _start_defrost(self, t: float) -> FloatArray:
        self.mode, self.t_mode = "defrost", t
        self.run_since = self.run_cold = 0.0
        self.defrosts += 1
        self.s_prev = -1.0
        return np.array([-1.0, 0.0, self.ventilation(t)])

    def _due(self, y: FloatArray) -> bool:
        """Count the previous sample's running time and return True if a defrost is due."""
        T_coil = y[IDX_Y_TCOIL_OUT]
        if self.s_prev > 0.0:
            # the previous command ran the compressor for one sample
            self.run_since += self.Ts
            if T_coil < self.T_cold:
                self.run_cold += self.Ts
        if self.strategy == "timed":
            return self.run_cold >= self.t_timed
        if self.run_cold >= self.t_max:
            return True
        if self.s_prev > 0.0 and self.run_since >= self.t_min:
            clean = self.hp.clean_coil_temperature(self.s_prev, y[IDX_Y_TA], y[IDX_Y_TOUT])
            return bool(clean - T_coil >= self.dT_demand)
        return False

    def step(self, t: float, y: FloatArray) -> FloatArray:
        """Return the command for one sample: the room controller's, or a defrost.

        Parameters
        ----------
        t : float
            Current time, s.
        y : ndarray
            Measurement [T_a degC, RH -, C ppm, T_coil,out degC, T_o degC].

        Returns
        -------
        ndarray, shape (3,)
            Command [s (-; -1 in defrost), f (-), V_vent (m3/s)].
        """
        T_coil = y[IDX_Y_TCOIL_OUT]
        if self.mode == "defrost":
            if T_coil >= self.T_end or t - self.t_mode >= self.t_end_max:
                self.mode, self.t_mode = "drip", t
            else:
                return np.array([-1.0, 0.0, self.ventilation(t)])
        if self.mode == "drip":
            if t - self.t_mode < self.t_drip:
                self.s_prev = 0.0
                return np.array([0.0, 0.0, self.ventilation(t)])
            self.mode = "normal"
        if self._due(y):
            return self._start_defrost(t)
        cmd = self.controller.step(t, y)
        self.s_prev = float(cmd[0])
        return cmd
