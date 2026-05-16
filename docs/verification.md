# Verification

## Checks

All checks run in `pytest` (151 tests, about a minute and a half). The code was also audited
from scratch, without the development history; its findings and fixes are
in [docs/decisions.md](decisions.md) (D20-D22).

| Check | Result |
|---|---|
| Psychrometrics: p_ws(20 °C), W(25 °C, 50 %), T_d(20 °C, 50 %) | 2333 Pa, 0.00986 kg/kg, 9.26 °C |
| Steady state vs 2×2 linear solution (3 cases) | < 10⁻⁶ K |
| Daily mean of periodic solution vs steady state of mean inputs | < 2·10⁻³ K |
| One-node step vs exponential, convergence order | RK4 **4.04**, Euler **1.01** |
| Energy, water and CO2 closure, flows integrated as extra states | round-off (~10⁻¹³ relative) |
| Energy residual with independent Simpson quadrature, order | RK4 **4.16**, Euler **1.05** |
| Moisture decay without buffer, fitted τ vs 3600/n | 10⁻⁸ relative |
| Moisture decay with buffer vs linearised two-exponential solution | < 10⁻³ of the perturbation |
| CO2 decay with ventilation vs exponential | 10⁻⁹ of the amplitude |
| Full scenario vs `scipy` Radau (rtol 10⁻¹¹) | < 10⁻⁶ K |
| Coil: Q_sens + Q_lat vs ṁ Δh; outlet RH ≤ 100 % on 168 conditions | round-off; always |
| Room + coil water balance in closed loop (drain integrated) | round-off |
| Unit inside the integrator, closed-loop convergence order | **3.9** |
| PI anti-windup: integrator under saturation vs analytical value | 10⁻⁹ |
| Cycle COP over 0-45 °C, all speeds and fans: positive, bounded by the minimum lift, continuous | always |
| Sigma-delta ripple in closed loop vs design (0.5 K) | median within 10 % |
| Sigma-delta: mean delivered demand, cycle length vs design formula | ±0.01, ±10 % |
| Cascade inner loop on the real coil: coil temperature vs its set point | ±0.05 K |
| Control model: exact discretisation; fixed point = continuous steady state | 10⁻¹² |
| Control model vs plant with the unit off, one day | < 0.01 K |
| Kalman filter: mass temperature and unknown 250 W heat input after 6 days | 10⁻³ K, 0.1 % |
| MPC plan: feasible, inside the planning band, no cooling while the band is wide | always |

![Convergence](../figures/convergence.png)

## Notes on the checks

The energy balance is checked at two levels. Any Runge-Kutta method preserves
linear invariants, so integrating the boundary flows alongside the state closes the balance to
round-off whatever the step: that verifies signs and bookkeeping. Integrating the logged flows
with an independent quadrature (Simpson) gives a residual that falls with the order of the
method: that verifies accuracy.

**Coupling the unit.** The coil output depends on the return-air state, so it is evaluated at
every Runge-Kutta stage while only the commands are held between samples. Evaluating the unit
once per sample and holding Q_sens instead drops the closed loop to first order: at Δt = 60 s
the error grows from 10⁻¹⁰ K to 2·10⁻⁴ K.

**Stiffness.** The fastest mode is the air-buffer moisture exchange, and it stiffens with
temperature because the buffer flux scales with p_ws(T_m). The explicit-Euler stability limit is
≈ 790 s at 26 °C and ≈ 380 s at 48 °C; RK4 at 60 s has ample margin.

![Energy balance](../figures/energy_balance.png)

![Humidity](../figures/humidity.png)
