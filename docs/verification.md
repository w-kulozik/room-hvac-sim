# Verification

## Checks

All checks run in `pytest` (282 tests, about three and a half minutes on a laptop). The code and documents
were checked again at the end of the project; the defects found and their fixes are in
[decisions.md](decisions.md) (D20-D22). The unit model was revised after an external review
(D29-D32); the tests added then are marked "(D29-D32)" below. A final audit corrected further
defects and added robustness checks (D35-D42); the tests added then are marked "(D35-D42)".
The seasonal study added measured weather and the seasonal efficiency (D43-D48); its tests are
marked "(D43-D48)".

RK4 is the classical fourth-order Runge-Kutta method, PI a proportional-integral controller and
MPC model-predictive control. Convergence orders are fitted by least squares over a sequence of
halved time steps.

| Check | Result |
|---|---|
| Psychrometrics: p_ws(20 °C), W(25 °C, 50 %), T_d(20 °C, 50 %) | 2333 Pa, 0.00986 kg/kg, 9.26 °C |
| Steady state vs 2×2 linear solution (3 cases) | < 10⁻⁶ K |
| Daily mean of periodic solution vs steady state of mean inputs | < 2·10⁻³ K |
| One-node step vs exponential, convergence order | RK4 4.04, Euler 1.01 |
| Energy, water and CO2 closure, flows integrated as extra states | round-off (~10⁻¹³ relative) |
| Energy residual with independent Simpson quadrature, order | RK4 4.16, Euler 1.05 |
| Moisture decay without buffer, fitted τ vs 3600/n | 10⁻⁸ relative |
| Moisture decay with buffer vs linearised two-exponential solution | < 10⁻³ of the perturbation |
| CO2 decay with ventilation vs exponential | 10⁻⁹ of the amplitude |
| Full scenario vs `scipy` Radau (rtol 10⁻¹¹) | < 10⁻⁶ K |
| Coil: Q_sens + Q_lat vs ṁ Δh; outlet RH ≤ 100 % on 168 conditions | round-off; always |
| Partially wet coil: wetted fraction 0 → 1 over an inlet-humidity sweep, no jump in capacity, sensible heat or water removal; values on both sides of the two regime boundaries (D29-D32) | jump-free; 10⁻⁶ relative at the boundaries |
| Partially wet coil: air-side enthalpy drop vs capacity, no water added, no supersaturation, surface between refrigerant and air (D29-D32) | 10⁻¹² relative; always |
| Capacity falls with airflow; at 70 % flow vs the EnergyPlus curve WindACCoolCapFFF (0.94) (D29-D32) | monotone; 0.933 |
| Capacity slopes at rating with entering wet bulb and with outdoor temperature vs the EnergyPlus curve WindACCoolCapFT (+3.58 and −1.08 %/K) (D29-D32) | ±0.1 and ±0.04 %/K |
| Room + coil water balance in closed loop (drain integrated) | round-off |
| Unit inside the integrator, constant command, coil in one regime: convergence order | 4.05 |
| Closed loop, PI sampled every 60 s, integration step 60-15 s, coil in one regime: order (D29-D32) | 4.02 |
| Constant command crossing from a fully wet into a partially wet coil: order, error at 60 s (D29-D32) | 2.10, 4·10⁻¹⁰ K |
| Unit evaluated once per sample and Q_sens held: order | 1.13 |
| PI anti-windup: integrator under saturation vs analytical value | 10⁻⁹ |
| Cycle COP over 0-45 °C, all speeds and fans: positive, bounded by the minimum lift, continuous | always |
| Sigma-delta ripple in closed loop vs design (0.5 K) | median within 10 % |
| Sigma-delta: mean delivered demand, cycle length vs design formula | ±0.003, ±10 % |
| Cascade inner loop on the real coil: coil temperature vs its set point | ±0.05 K |
| Control model: exact discretisation; fixed point = continuous steady state | 10⁻¹² |
| Control model vs plant with the unit off, one day | < 10⁻⁴ K (measured 2·10⁻⁵ K) |
| Kalman filter: mass temperature and unknown 250 W heat input after 6 days | 10⁻³ K, 0.1 % |
| MPC plan: feasible, inside the planning band, no cooling while the band is wide | always |
| MPC set point clipped with the band at the end of the first prediction step: solve at 18:00 (end unoccupied) and at 7:45 (end occupied) (D35-D42) | 27.0 °C kept; clipped to 25.4 °C |
| MPC capacity and cost per unit of cooling vs the net cooling Q_sens − P_fan_in fed to the Kalman filter (D35-D42) | 10⁻¹² relative |
| MPC after a failed LP solve: set point and move taken from the previous plan shifted by one step (D35-D42) | exact |
| Compressor starts in a window that begins inside a running cycle; fan travel counts only steps between consecutive running samples, NaN without running time (D35-D42) | no start counted; exact on constructed records |
| Frost-limited unit: COP at the effective speed s_eff (20 °C, 10 % RH return air, −10 °C outdoors, s = 1, 40 % fan: s_eff = 0.37); cycling at s_min below it (D35-D42) | 10⁻¹² relative to the formula (COP 7.23 against 9.05 at the commanded speed) |
| Frost-limited unit: power and COP at frost onset and where s_eff crosses s_min (D35-D42) | jumps 7·10⁻⁶ W and 5·10⁻⁷; 1·10⁻⁶ W and 1·10⁻⁸ |
| Fan-only coil: surface temperature and wetted fraction as the film vanishes (D35-D42) | continuous, (T_in, 0) on a dry coil |
| Divergence detection: RK4 with Δt = 1 200 s on the summer scenario; Δt = 900 s still runs (D35-D42) | stops at t = 31 200 s, limit reported ≈ 1 094 s; runs |
| Zero-capacity moisture buffer with φ_b = 0.9, one day (D35-D42) | same W as with the buffer off; water residual 0 kg |
| Periodic state with the buffer off (φ_b has no dynamics), guesses φ_b = 0 and 0.5 (D35-D42) | φ_b kept; other states identical, periodic to 3·10⁻¹⁰; open mass integrator raises a clear error |
| Input validation: negative room parameters, negative ventilation flow, RH_o outside [0, 1], sunset before sunrise, empty or out-of-range daily windows, window outside the run, no occupied rows (D35-D42) | ValueError, or NaN with a warning for the occupied means |
| Daily windows crossing midnight (22-6 h): occupancy, ventilation schedule, comfort schedule, tariff peak (D35-D42) | active at 23 h and 2 h, not at 12 h |
| Controller sampling period vs `dt_ctrl` (seven controllers at 30 s against their 60 s), and parts of one controller with different periods (D35-D42) | ValueError |
| Cached exact discretisation of the control model is read-only (D35-D42) | write raises |
| Solar position vs the low-precision solar coordinates of the astronomical almanac, over four years of the leap cycle (D43-D48) | declination ±0.6°, equation of time ±1.2 min (measured 0.3-0.5°, 0.7-0.9 min) |
| Angle of incidence vs the dot product of sun and surface-normal vectors, seven orientations; vertical south wall at solar noon vs sin(φ − δ) (D43-D48) | 10⁻¹² |
| Time convention of the weather data: horizontal irradiance rebuilt from beam and diffuse with our solar geometry, Munich and Madrid (D43-D48) | RMS 1.9 and 1.8 W/m² (8.3 and 12.4 W/m² without the irradiance time offset) |
| Weather reader: rows as in the file, malformed files, periodic year, interpolation, dew point kept, local clock for the schedules (D43-D48) | exact; ValueError |
| SEER of the real unit from its published points and low-power modes (D43-D48) | 8.65 and 101 kWh, as published |
| SEER with the same EER at every bin; degradation rule; declared points of the unit model at the part load or at minimum speed (D43-D48) | 10⁻¹³; exact; 10⁻⁹ |
| Net cooling integrated with the state vs the trapezoid of the logged values (D43-D48) | 10⁻⁵ relative |
| On/off thermostat with a set-point schedule; with a constant set point identical to before (D43-D48) | follows the schedule; bit-identical |

COP is the coefficient of performance of the compressor.

![Error against time step on log-log axes for RK4 and Euler: one-node step response and energy-balance residual](../figures/convergence.png)

## Notes on the checks

**Energy balance.** It is checked at two levels. Any Runge-Kutta method preserves linear
invariants, so integrating the boundary flows alongside the state closes the balance to
round-off whatever the step: that verifies signs and bookkeeping. Integrating the logged flows
with an independent quadrature (Simpson) gives a residual that falls with the order of the
method: that verifies accuracy.

**Coupling the unit.** The coil output depends on the return-air state, so it is evaluated at
every Runge-Kutta stage while only the commands are held between samples. Evaluating the unit
once per sample and holding Q_sens instead drops the loop to first order: at Δt = 60 s the error
grows from 5·10⁻¹¹ K to 1·10⁻³ K.

**What the convergence tests show, and what they do not.** Fourth order is measured with a
constant command and in closed loop with the sampled PI, in runs where the coil stays in one
regime and the compressor modulates without switching. Two things lower the order:
- a trajectory that crosses from a fully wet into a partially wet coil, where the coil output is
  continuous but not differentiable: order 2 (2.10), with an error of 4·10⁻¹⁰ K at 60 s. Three
  exploratory runs crossing from a dry into a partially wet coil gave orders of 2.0, 4.0 and
  4.9, with irregular error sequences, so no order is stated for that boundary (the wet part
  starts there with zero length, so the kink is weaker);
- a switching decision (start, stop, fan step) that falls on a different sample for different
  integration steps. Then the runs differ by a whole switching event and there is no convergence
  order to measure. This belongs to the sampled controller, not to the integrator, and is not
  tested.
An earlier version of this page and the README said "fourth-order convergence, also in closed
loop"; the test behind it used a constant command (decision D8).

**Stiffness.** The fastest mode is the air-buffer moisture exchange, and it stiffens with
temperature because the buffer flux scales with p_ws(T_m). The explicit-Euler stability limit is
≈ 790 s at 26 °C and ≈ 380 s at 48 °C; the RK4 limit is ≈ 1 100 s at 26 °C, so RK4 at 60 s has
ample margin. With a step above the limit the simulator stops as soon as the state becomes
non-finite or unphysical (negative humidity ratio, CO2 or active-buffer φ_b) and reports the
limit estimated at the initial state, instead of returning diverged values.

**What the checks do not cover.** They verify that the equations are solved correctly. They do
not validate the equations against measurements of a real room or unit; see the known
limitations at the end of [decisions.md](decisions.md).

![Stored energy against integrated boundary heat flows over three days, with the residual of about 10⁻⁶ J](../figures/energy_balance.png)

![Outdoor, indoor and buffer relative humidity, and outdoor and indoor dew point, over three free-floating days](../figures/humidity.png)
