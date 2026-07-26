"""Low-load handling: turn a continuous demand into commands an inverter unit can follow.

The compressor cannot run below s_min, so a demand 0 < v < s_min has to be met by cycling.
Each handler returns (s, u_applied, integrate) for the PI:
    s          compressor command
    u_applied  the demand the unit actually meets, for back-calculation
    integrate  False to freeze the PI integrator

A. HysteresisLowLoad: what many commercial units do. Below s_min the unit runs at s_min until the
   room falls dT_off below the set point, stops, and restarts dT_on above it. The PI integrator is
   frozen meanwhile, otherwise the off periods would wind it down.

C. SigmaDeltaLowLoad: first-order sigma-delta modulation. The handler accumulates the difference
   between requested and delivered demand, starts when the accumulated deficit exceeds a
   threshold and stops when the surplus does. The average delivered demand equals the request,
   so the PI keeps integrating normally. With threshold A = s_min T_ref / 2 the on time is
   s_min T_ref / (s_min - v) and the off time s_min T_ref / v, i.e. at v = s_min/2 one cycle of
   4 T_ref. The accumulator is clamped to +-A: when the guard holds the compressor beyond a
   switching point, the extra deficit or surplus is dropped instead of being paid back later
   (paying it back kept the unit off for a long time while the room was already warm).

   Choosing T_ref. With Q_d the sensible cooling per unit of demand while cycling (at s_min),
   the unit delivers s_min Q_d between switches while the load takes v Q_d, and the energy
   mismatch A Q_d is stored in the air, so the air-temperature ripple is
       dT_pp = 2 A Q_d / C_a = s_min T_ref Q_d / C_a.
   ``ripple_T_ref`` inverts this for a target ripple.

Both respect the compressor guard (minimum run and stop times).
"""

from __future__ import annotations

from typing import Protocol

from .base import CompressorGuard


class LowLoadHandler(Protocol):
    """Interface of a low-load handler used by ``PIUnitController``."""

    guard: CompressorGuard

    def step(self, t: float, v: float, T: float, T_set: float) -> tuple[float, float, bool]:
        """Return (s, u_applied, integrate) for demand v at time t (s)."""
        ...


def ripple_T_ref(dT_pp: float, C_a: float, Q_d: float, s_min: float) -> float:
    """Return the sigma-delta time scale giving an air-temperature ripple dT_pp.

    Parameters
    ----------
    dT_pp : float
        Target air-temperature ripple, peak to peak, K.
    C_a : float
        Heat capacity of the air node, J/K.
    Q_d : float
        Sensible cooling per unit of demand while cycling at s_min, W.
    s_min : float
        Minimum compressor command, -.

    Returns
    -------
    float
        Time scale T_ref, s.
    """
    return dT_pp * C_a / (s_min * Q_d)


class HysteresisLowLoad:
    """Handler A: thermostatic cycling at s_min below the minimum speed, integrator frozen.

    Parameters
    ----------
    s_min : float
        Minimum compressor command, -.
    guard : CompressorGuard
        Minimum run and stop times.
    dT_on : float, optional
        Restart when the air is this far above the set point, K.
    dT_off : float, optional
        Stop when the air is this far below the set point, K.

    Attributes
    ----------
    cycling : bool
        Whether the last sample was handled by low-load cycling.
    """

    def __init__(self, s_min: float, guard: CompressorGuard, dT_on: float = 0.5,
                 dT_off: float = 0.5) -> None:
        self.s_min, self.guard = s_min, guard
        self.dT_on, self.dT_off = dT_on, dT_off
        self.cycling = False

    def step(self, t: float, v: float, T: float, T_set: float) -> tuple[float, float, bool]:
        """Return the compressor command for demand v at time t.

        Parameters
        ----------
        t : float
            Current time, s.
        v : float
            Unsaturated demand from the PI, -.
        T : float
            Measured air temperature, degC.
        T_set : float
            Set point, degC.

        Returns
        -------
        s : float
            Compressor command, -.
        u_applied : float
            Demand the unit meets, for back-calculation, -.
        integrate : bool
            False to freeze the PI integrator.
        """
        g = self.guard
        if v >= self.s_min:
            self.cycling = False
            if g.on or g.can_start(t):
                g.set(t, True)
                s = min(v, 1.0)
                return s, s, True
            return 0.0, 0.0, True                         # waiting for the minimum stop time
        # low load: thermostatic cycling at s_min, integrator frozen
        self.cycling = True
        if g.on and T_set - self.dT_off >= T and g.can_stop(t):
            g.set(t, False)
        elif not g.on and T_set + self.dT_on <= T and g.can_start(t):
            g.set(t, True)
        if v <= 0.0 and not g.on:
            # no demand at all: plain saturation at zero, let the integrator track it
            self.cycling = False
            return 0.0, 0.0, True
        return (self.s_min if g.on else 0.0), max(v, 0.0), False


class SigmaDeltaLowLoad:
    """Handler C: first-order sigma-delta modulation between 0 and s_min.

    Parameters
    ----------
    s_min : float
        Minimum compressor command, -.
    guard : CompressorGuard
        Minimum run and stop times.
    Ts : float, optional
        Sampling period, s.
    T_ref : float, optional
        Modulation time scale, s; the switching threshold is s_min T_ref / 2.

    Attributes
    ----------
    threshold : float
        Switching threshold of the accumulator, s.
    acc : float
        Integrated requested minus delivered demand, s.
    s_prev : float
        Compressor command of the previous sample, -.
    """

    def __init__(self, s_min: float, guard: CompressorGuard, Ts: float = 60.0,
                 T_ref: float = 600.0) -> None:
        self.s_min, self.guard, self.Ts = s_min, guard, Ts
        self.threshold = 0.5 * s_min * T_ref
        self.acc = 0.0              # requested minus delivered demand, integrated, s
        self.s_prev = 0.0

    def step(self, t: float, v: float, T: float, T_set: float) -> tuple[float, float, bool]:
        """Return the compressor command for demand v at time t.

        Parameters
        ----------
        t : float
            Current time, s.
        v : float
            Unsaturated demand from the PI, -.
        T : float
            Measured air temperature, degC (not used).
        T_set : float
            Set point, degC (not used).

        Returns
        -------
        s : float
            Compressor command, -.
        u_applied : float
            Demand the unit meets on average, for back-calculation, -.
        integrate : bool
            Always True.
        """
        g = self.guard
        if v >= self.s_min:
            # modulating range: follow the demand directly
            self.acc = 0.0
            if g.on or g.can_start(t):
                g.set(t, True)
            s = min(v, 1.0) if g.on else 0.0
            u_applied = s
        elif v <= 0.0:
            # no demand: stop as soon as the guard allows
            self.acc = 0.0
            if g.on and g.can_stop(t):
                g.set(t, False)
            s = self.s_min if g.on else 0.0
            u_applied = 0.0
        else:
            # low load: sigma-delta modulation between 0 and s_min
            self.acc += (v - self.s_prev) * self.Ts
            self.acc = min(max(self.acc, -self.threshold), self.threshold)
            if not g.on and self.acc >= self.threshold and g.can_start(t):
                g.set(t, True)
            elif g.on and self.acc <= -self.threshold and g.can_stop(t):
                g.set(t, False)
            s = self.s_min if g.on else 0.0
            u_applied = v                                 # met on average
        self.s_prev = s
        return s, u_applied, True
