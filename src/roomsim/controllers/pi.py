"""Discrete PI controller with back-calculation anti-windup, and SIMC tuning for the room.

Positional form, sampled every Ts (cooling: a room warmer than the set point raises the demand):
    e_k  = y_k - r_k
    v_k  = Kp e_k + I_k                       unsaturated demand
    u_k  = what the actuator actually delivers (decided outside the PI)
    I_k+1 = I_k + Kp Ts/Ti e_k + Ts/Tt (u_k - v_k)
The last term (back-calculation) bleeds the integrator whenever the actuator cannot follow v,
so the integral never winds up beyond what the actuator can deliver. Under a constant error e
with the actuator saturated at u_max the integrator settles at
    I = u_max - Kp e (1 - Tt/Ti),
so a smaller Tt keeps it closer to the limit. Tt = Ti/2 is a choice made here for fast recovery
(decision D12); with a noisy sensor a longer Tt may be needed. With integrate=False the
integrator is frozen (conditional integration), used while the unit cycles at low load.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..ac_unit import ACUnit
from ..params import RoomParams


@dataclass
class DiscretePI:
    """Positional discrete PI with back-calculation anti-windup.

    Attributes
    ----------
    Kp : float
        Proportional gain, output units per unit of error (e.g. demand per K).
    Ti : float
        Integral time, s.
    Ts : float
        Sampling period, s.
    Tt : float or None
        Tracking time for anti-windup, s; None uses Ti/2, inf disables anti-windup.
    I : float
        Integrator state, output units.
    """

    Kp: float                 # demand per K
    Ti: float                 # integral time, s
    Ts: float = 60.0          # sampling period, s
    Tt: float | None = None   # tracking time for anti-windup, s (None: Ti/2; inf: no anti-windup)
    I: float = 0.0            # integrator state  # noqa: E741  (public attribute, PI notation)

    def demand(self, e: float) -> float:
        """Return the unsaturated demand v = Kp e + I for error e (no state change).

        Parameters
        ----------
        e : float
            Error y - r, in the units of the measurement.

        Returns
        -------
        float
            Unsaturated demand, output units.
        """
        return self.Kp * e + self.I

    def update(self, e: float, v: float, u_applied: float, integrate: bool = True) -> None:
        """Advance the integrator by one sample, with back-calculation.

        Parameters
        ----------
        e : float
            Error used for the demand, in the units of the measurement.
        v : float
            Unsaturated demand returned by ``demand``, output units.
        u_applied : float
            Output actually delivered, output units.
        integrate : bool, optional
            False freezes the integrator (conditional integration).
        """
        if not integrate:
            return
        Tt = 0.5 * self.Ti if self.Tt is None else self.Tt
        self.I += self.Kp * self.Ts / self.Ti * e + self.Ts / Tt * (u_applied - v)


def simc_pi(K: float, tau: float, theta: float, tau_c: float) -> tuple[float, float]:
    """Return SIMC PI settings (Skogestad 2003) for a process K e^(-theta s)/(tau s + 1).

    Parameters
    ----------
    K : float
        Process gain, output per unit of controller output.
    tau : float
        Time constant, s.
    theta : float
        Dead time, s.
    tau_c : float
        Desired closed-loop time constant, s.

    Returns
    -------
    Kp : float
        Proportional gain, in the inverse of the units of K.
    Ti : float
        Integral time, s.
    """
    Kp = tau / (K * (tau_c + theta))
    Ti = min(tau, 4.0 * (tau_c + theta))
    return Kp, Ti


def room_pi_tuning(p: RoomParams, unit: ACUnit, Ts: float = 60.0, tau_c: float = 900.0,
                   T_o: float = 32.0) -> tuple[float, float]:
    """Return PI gains for air-temperature control through the compressor demand.

    Process seen by the controller: demand in [0, 1] -> cooling -> T_a. With the mass node held
    fixed this is first order with tau = tau_air and gain Q_d / (H_am + H_ao + H_inf), where
    Q_d is the cooling per unit of demand. The SIMC rules are applied with Q_d = the total
    compressor capacity at full speed and the rated evaporating temperature (2620 W at 32 degC
    outdoors). This is a nominal value: T_a responds only to the net sensible cooling, and the
    unit is strongly nonlinear in the demand v. Below s_min (sigma-delta cycling, where the PI
    spends most of the time in the standard runs) the unit delivers on average about 3250 W per
    unit of demand; above s_min the coil is wet and the incremental sensible gain is
    700-1000 W. Measured in set-point steps of 0.5 K, the closed-loop 63 % time is 600-1140 s
    at low load (design 960 s, i.e. tau_c + theta plus the zero-order hold) and 2100-2340 s in
    the modulating range. Taking Q_d from the low-load sensible gain instead hardly changes
    the low-load response, because the sigma-delta modulation dominates it, and slows the
    modulating range further.

    The mass node adds a slow mode that the integral action handles. Dead time: half a sample
    for the zero-order hold; the sensor is read at the sampling instant without further delay,
    so theta = Ts/2. tau_c = 900 s trades speed for robustness to the unmodelled mass dynamics
    and the coil nonlinearity.

    Parameters
    ----------
    p : RoomParams
        Room parameters (conductances in W/K, ``tau_air`` in s).
    unit : ACUnit
        Air-conditioning unit; its full-speed capacity at T_o sets the process gain.
    Ts : float, optional
        Sampling period, s.
    tau_c : float, optional
        Desired closed-loop time constant, s.
    T_o : float, optional
        Outdoor temperature at which the capacity is evaluated, degC.

    Returns
    -------
    Kp : float
        Proportional gain, demand per K.
    Ti : float
        Integral time, s.
    """
    K = unit.capacity(1.0, T_o) / (p.H_am + p.H_ao + p.H_inf)
    return simc_pi(K, p.tau_air, Ts / 2.0, tau_c)
