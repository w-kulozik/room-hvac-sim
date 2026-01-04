"""Water held on the coil: filling, draining and re-evaporation (part-load latent degradation)."""

import numpy as np
import pytest

from roomsim import psychro
from roomsim.ac_unit import ACUnit, _h
from roomsim.inputs import HOUR, ConstantInputs
from roomsim.params import RoomParams
from roomsim.plant import equilibrium_state
from roomsim.simulate import simulate

from criteria import RH_RESOLUTION_W

P = RoomParams()
UNIT = ACUnit()
M_MAX = UNIT.p.M_coil_max
HUMID = ConstantInputs(T_o=30.0, W_o=0.0150, G_sol=150.0, n_occ=2)
X0 = equilibrium_state(26.0, 26.5, 0.0120)


class Schedule:
    """Compressor on/off with a fixed period; fan follows the compressor or runs continuously."""

    def __init__(self, on_s, off_s, fan_continuous, s=0.6):
        self.on_s, self.period, self.fan_cont, self.s = on_s, on_s + off_s, fan_continuous, s

    def step(self, t, y):
        on = (t % self.period) < self.on_s
        return np.array([self.s if on else 0.0, 1.0 if (on or self.fan_cont) else 0.0, 0.0])


# ---- unit level --------------------------------------------------------------------------

def test_coil_water_rates_are_consistent():
    for M_c in (0.0, 0.5 * M_MAX, 0.97 * M_MAX, M_MAX):
        st = UNIT.operate(0.8, 1.0, 26.0, 0.0120, 32.0, M_c=M_c)
        assert st.dM_coil + st.m_drain == pytest.approx(st.m_lat, rel=1e-14)
        assert st.m_drain >= 0.0 and st.dM_coil >= 0.0
    assert UNIT.operate(0.8, 1.0, 26.0, 0.0120, 32.0, M_c=0.5 * M_MAX).m_drain == 0.0
    assert UNIT.operate(0.8, 1.0, 26.0, 0.0120, 32.0, M_c=M_MAX).dM_coil == 0.0


def test_fan_only_re_evaporation_is_adiabatic():
    T_in, W_in = 25.0, 0.0100
    st = UNIT.operate(0.0, 1.0, T_in, W_in, 32.0, M_c=0.5 * M_MAX)
    assert st.m_lat < 0.0                                          # water returns to the air
    assert st.W_out > W_in and st.T_out < T_in
    assert _h(st.T_out, st.W_out) == pytest.approx(_h(T_in, W_in), rel=1e-13)
    assert psychro.RH_from_W(st.T_out, st.W_out) <= 1.0
    assert st.dM_coil == pytest.approx(st.m_lat)
    assert st.Q_coil == 0.0 and st.P_el == pytest.approx(UNIT.p.P_fan_in_nom)


def test_dry_coil_fan_only_changes_nothing_but_fan_heat():
    st = UNIT.operate(0.0, 1.0, 25.0, 0.0100, 32.0, M_c=0.0)
    assert st.m_lat == 0.0 and st.Q_sens == pytest.approx(0.0, abs=1e-9)


def test_rates_are_continuous_at_the_full_and_dry_limits():
    eps = 1e-9
    full = [UNIT.operate(0.8, 1.0, 26.0, 0.012, 32.0, M_c=M_MAX + d).m_drain for d in (-eps, eps)]
    assert full[0] == pytest.approx(full[1], rel=1e-6)
    dry = [UNIT.operate(0.0, 1.0, 25.0, 0.010, 32.0, M_c=d).m_lat for d in (eps, 0.0)]
    assert dry[0] == pytest.approx(dry[1], abs=1e-9)


# ---- in the loop -------------------------------------------------------------------------

def test_coil_fills_then_drains():
    r = simulate(P, HUMID, X0, 1 * HOUR, dt=10.0, controller=Schedule(4 * HOUR, 0, False, s=1.0),
                 unit=UNIT)
    M, drain = r.M_coil, r.unit["m_drain"]
    assert M[0] == 0.0 and np.all(np.diff(M) >= -1e-15)
    k_full = np.argmax(M > 0.95 * M_MAX)
    assert drain[:k_full].max() == 0.0                 # nothing drains while the coil wets
    assert drain[-1] > 0.9 * r.unit["m_lat"][-1]       # then condensate leaves at its own rate
    # time to wet the coil: M_max / condensation rate (rate ramps up after the start)
    t_wet = r.t[k_full] - r.t[0]
    assert 0.8 * M_MAX / r.unit["m_lat"][k_full] < t_wet < 2.0 * M_MAX / r.unit["m_lat"][k_full]


def test_water_of_room_and_coil_is_conserved():
    """Room air + buffer + coil film exchange water only with outdoor air, occupants and drain."""
    r = simulate(P, HUMID, X0, 3 * HOUR, controller=Schedule(20 * 60, 20 * 60, True), unit=UNIT)
    assert r.M_drain[-1] > 0.0 and r.M_coil.max() > 0.9 * M_MAX     # the check is not trivial
    assert np.abs(r.system_water_residual).max() < 1e-13
    assert np.abs(r.water_residual).max() < 1e-13                    # room alone, as before
    assert r.M_coil.min() >= -1e-12 and r.M_coil.max() <= M_MAX + 1e-12


def test_continuous_fan_returns_moisture_to_the_room():
    """Part-load latent degradation: same compressor cycling, fan left on between cycles."""
    runs = {fan: simulate(P, HUMID, X0, 6 * HOUR, controller=Schedule(15 * 60, 15 * 60, fan),
                          unit=UNIT) for fan in (False, True)}
    W_auto, W_cont = (runs[f].W[-int(2 * HOUR / 60):].mean() for f in (False, True))
    assert W_cont > W_auto + RH_RESOLUTION_W           # more humid by more than 1 %RH
    drained = {f: np.sum(r.unit["m_drain"]) * r.dt for f, r in runs.items()}
    assert drained[True] < 0.5 * drained[False]          # most condensate never reaches the drain
