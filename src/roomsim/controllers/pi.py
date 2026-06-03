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
    Kp: float                 # demand per K
    Ti: float                 # integral time, s
    Ts: float = 60.0          # sampling period, s
    Tt: float | None = None   # tracking time for anti-windup, s (None: Ti/2; inf: no anti-windup)
    I: float = 0.0            # integrator state

    def demand(self, e: float) -> float:
        """Unsaturated demand v for error e (no state change)."""
        return self.Kp * e + self.I

    def update(self, e: float, v: float, u_applied: float, integrate: bool = True) -> None:
        if not integrate:
            return
        Tt = 0.5 * self.Ti if self.Tt is None else self.Tt
        self.I += self.Kp * self.Ts / self.Ti * e + self.Ts / Tt * (u_applied - v)


def simc_pi(K: float, tau: float, theta: float, tau_c: float) -> tuple[float, float]:
    """SIMC rules (Skogestad 2003) for a first-order process K e^(-theta s)/(tau s + 1)."""
    Kp = tau / (K * (tau_c + theta))
    Ti = min(tau, 4.0 * (tau_c + theta))
    return Kp, Ti


def room_pi_tuning(p: RoomParams, unit: ACUnit, Ts: float = 60.0, tau_c: float = 900.0,
                   T_o: float = 32.0) -> tuple[float, float]:
    """PI gains for air temperature control through the compressor demand.

    Process seen by the controller: demand in [0, 1] -> cooling capacity -> T_a. With the mass
    node held fixed this is first order with tau = tau_air and gain capacity / (H_am + H_ao + H_inf).
    The mass node adds a slow mode that the integral action handles. Dead time: half a sample
    for the zero-order hold; the sensor is read at the sampling instant without further delay,
    so theta = Ts/2.
    tau_c = 900 s trades speed for robustness to the unmodelled mass dynamics and the coil
    nonlinearity.
    """
    K = unit.capacity(1.0, T_o) / (p.H_am + p.H_ao + p.H_inf)
    return simc_pi(K, p.tau_air, Ts / 2.0, tau_c)
