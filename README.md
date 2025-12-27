# roomsim: a dynamic room model for HVAC control

A small simulator of a single room with an inverter split air conditioner, verified against
analytical solutions, used to test and compare air-conditioning controllers.

The room has an air node, a thermal-mass node, a moisture balance with buffering by the
interior surfaces, a CO2 balance and proper psychrometrics. It is driven by outdoor weather,
sun through a window and occupants coming and going. The unit model resolves the cooling coil
(dry and wet), the water held on the coil, frost protection, part-load efficiency and start-up
losses. Controllers are discrete: they sample the sensors every Δt_c and hold their commands,
as a real unit does.

![Controllers over one day](figures/controllers_day.png)

## Model

State x = [T_a, T_m, W, φ_b, C] for the room, plus the water M_c held on the coil when a unit
is present. Plant input u = [Q_hvac_sens, ṁ_hvac_lat, V_vent].

```
C_a dT_a/dt  = H_am (T_m - T_a) + (H_ao + H_inf + H_vent)(T_o - T_a) + f_sol,a Q_sol + Q_int,conv + Q_hvac
C_m dT_m/dt  = H_am (T_a - T_m) + H_mo (T_o - T_m) + (1 - f_sol,a) Q_sol + Q_int,rad
ρV  dW/dt    = ρ(V n/3600 + V_vent)(W_o - W) + N g_vap - ṁ_hvac_lat + ṁ_buf
C_b dφ_b/dt  = -ṁ_buf,      ṁ_buf = G_b (φ_b p_ws(T_m) - p_w(W))
V   dC/dt    = (V n/3600 + V_vent)(C_o - C) + N G_CO2
```

| Symbol | Meaning | Value |
|---|---|---|
| C_a | air + light furnishings | 0.33 MJ/K |
| C_m | walls, floor, heavy furniture | 10 MJ/K |
| H_am | air ↔ mass, h_i·A_int | 180 W/K |
| H_ao, H_mo, H_inf | window, opaque envelope, infiltration (0.5 ach) | 4, 25, 9.05 W/K |
| C_b, G_b | moisture buffer (EMPD, painted gypsum) | 1.44 kg, 5.2·10⁻⁷ kg/(s·Pa) |
| G_CO2 | CO2 per person (1.2 met) | 0.0052 L/s |

Room 5 × 4 × 2.7 m, two occupants 8-18 h (70 W, 0.05 kg/h each), 3 m² window (g = 0.6),
sinusoidal summer day 20-32 °C. Full specification, including where and why it departs from
its first version: [docs/specification.md](docs/specification.md).

Time scales: air ≈ 1 700 s, mass ≈ 3.1 days, air-buffer moisture exchange ≈ 900 s,
apparent room moisture ≈ 15 000 s (twice the 7 200 s of the air alone).

### Air-conditioning unit

Inverter split unit, 2.5 kW at 35 °C outdoor, EER 3.5 at rating conditions. The controller sends
compressor speed s (0 or 0.3-1), indoor fan speed f and ventilation flow; the unit converts them
into Q_hvac_sens and ṁ_hvac_lat at every Runge-Kutta stage.

| Item | Model |
|---|---|
| Coil | bypass factor / apparatus dew point, BF = 0.15 at 450 m³/h, NTU ∝ ṁ^-0.4 |
| Dry / wet coil | closed form when dry; Newton on the saturation enthalpy when wet; continuous at the switch |
| Sensible / latent split | exact for the enthalpy definition: Q_sens + Q_lat = ṁ Δh to round-off |
| Cycle | Carnot fraction with finite-UA evaporator and condenser, closed form, calibrated to the rated EER |
| Frost protection | capacity capped so the evaporating temperature stays ≥ 0.5 °C |
| Start-up | capacity ramps with a 60 s time constant while the compressor draws full power |
| Coil water | condensate wets the coil (0.25 kg) before it drains; re-evaporates with the fan on and the compressor off |

At rating conditions the sensible heat ratio is 0.73 and the supply air 15 °C. Lower fan speed
gives a colder coil, more dehumidification and a lower COP; part load raises the EER to 4.0.

## Results

### Free-floating baseline

Without air conditioning, on the **periodic steady state** found by Newton shooting on the
one-day map, the room peaks at **T_a,max = 36.8 °C at 15:07** (minimum 33.4 °C, indoor RH
35-39 %). With 0.5 ach, an unshaded window and about 340 W of mean gains, the room sits ~9 K
above the mean outdoor temperature. Starting from a uniform 26 °C instead underestimates the
peak by almost 4 K on day 3, because the mass takes about three days to respond.

![Free-floating temperatures](figures/temperatures.png)

### Standard controllers

Three summer days after three warm-up days, set point 25 °C, ventilation 38 m³/h during
occupied hours (mean CO2 950 ppm):

| Controller | Energy (kWh) | Starts | T rms occ. (K) | RH mean occ. (%) | h RH > 60 % occ. |
|---|---|---|---|---|---|
| On/off thermostat (±0.5 K, full speed) | 11.38 | 174 | 0.35 | 49.2 | 0.0 |
| On/off, fan continuous | 9.92 | 117 | 0.40 | 64.1 | 26.2 |
| PI + hysteresis low-load handling (A) | 6.11 | 66 | 0.25 | 64.1 | 28.4 |
| PI + sigma-delta low-load handling (C) | 6.21 | 123 | 0.22 | 63.7 | 27.4 |
| PI (C), fan continuous | 7.21 | 126 | 0.22 | 63.7 | 27.4 |

![Controller metrics](figures/controllers_metrics.png)

**The PI halves the energy but loses humidity control.** The unit is large for this office, as
is common in practice, so the PI keeps it at minimum speed most of the time. At s = 0.3 and full
fan the coil surface sits near 19 °C, above the dew point of the room air (about 16 °C): the coil
runs dry and removes no water at all (sensible heat ratio 1.00 against 0.80 for on/off). The
mean EER rises from 2.6 to 4.5, but indoor RH climbs to 64 %. Recovering dehumidification at part
load is the subject of the next section.

**Leaving the fan on between cycles returns the condensate to the room.** With the on/off
thermostat, continuous fan raises the mean RH from 49 % to 64 %: the water held on the coil
re-evaporates during the off periods, and most of it never reaches the drain. This is the
part-load latent degradation reported for cycling units. Under PI control the unit modulates
rather than cycles during occupied hours, so the effect there is small.

**Low-load handling.** Below the 30 % minimum speed the compressor has to cycle. The PI computes
a continuous demand and a separate handler turns it into feasible commands, with minimum run and
stop times of 3 min. Handler A (thermostatic hysteresis, integrator frozen) is what many units
do; handler C (sigma-delta) meets the mean demand, so the PI keeps integrating. The sigma-delta
time scale follows from the ripple it causes: between switches the energy mismatch is stored in
the air, so the peak-to-peak ripple is s_min T_ref Q_max / C_a, and T_ref = 320 s gives 0.5 K.
Its accumulator is clamped so that a deficit held back by the minimum run or stop time is not
paid back later, which had kept the unit off while the room was already warm. Handler C ends up
with the tightest temperature control and more, but shorter, cycles than A (3.3 starts per
running hour against 11 for on/off).

### Humidity-aware fan control

The indoor fan is the dehumidification actuator: less airflow means a colder coil and a lower
sensible heat ratio. At minimum compressor speed, going from full to 40 % fan takes the coil from
19 °C (dry) to 15 °C (wet), SHR 1.00 to 0.76, for a 5 % loss of EER. The compressor stays under the
same PI temperature loop; two structures drive the fan:

- **H1**: a relative-humidity PI acting directly on the fan speed.
- **H2, cascade**: the RH loop sets a coil-temperature target, and an inner loop on the
  evaporator thermistor moves the fan to reach it.

Set points 25 °C and 55 % RH:

| Controller | Energy (kWh) | Water removed (kg) | RH mean occ. (%) | RH std occ. (%) | h RH > 60 % occ. |
|---|---|---|---|---|---|
| PI (C), full fan | 6.21 | 0.0 | 63.7 | 2.07 | 27.4 |
| H1 fan RH-PI | 7.27 | 5.2 | 55.3 | 0.81 | 0.0 |
| H2 cascade | 7.28 | 5.2 | 55.2 | 0.50 | 0.0 |
| On/off (reference) | 11.38 | 8.6 | 49.2 | 1.54 | 0.0 |

![Humidity control](figures/humidity_control.png)

Both bring the humidity back to its set point for 17 % more energy than the plain PI, still 36 %
less than on/off. The two structures have the same mean performance; the cascade holds humidity
more steadily and moves the fan about a third less. Keeping the controllers tuned for the nominal
unit and running them on a different one (dirty coil with bypass factor 0.30, or a 3.5 kW unit)
gives the same ordering: RH standard deviation 0.48-0.64 % for the cascade against 0.81-0.97 %
for H1. The advantage is real but small, as expected when the outer process (room humidity) is
much slower than anything the inner loop rejects.

### Model-predictive control

An economic MPC plans the air temperature over the next 24 h in 15-minute steps and re-plans
every 15 minutes. It minimises the forecast electrical energy (the unit's COP is better in the
cool hours), weighted by the electricity price when there is a tariff, plus a penalty for leaving
the comfort band (25 ± 1 °C while occupied, up to 28 °C otherwise). The prediction model is the
linear two-node RC model, discretised exactly; it knows nothing about humidity, the coil or the
unit's part-load behaviour. All terms are linear, so each plan is a linear programme solved by
HiGHS in about 10 ms.

A Kalman filter on the air-temperature sensor estimates the unmeasured mass temperature and an
unknown heat input to the air. The second state makes the predictions offset-free: whatever the
model gets wrong shows up as a disturbance that the planner then compensates. Started 3 K wrong
on the mass temperature, the estimate converges without bias but slowly (time constant ~14 h),
because with a single air sensor a wrong mass temperature and an unknown heat input can only be
told apart through the slow mass mode.

The MPC supervises the PI rather than driving the compressor directly. Applying the planned
cooling directly, with the low-load handler turning it into cycling between the 15-minute
solves, left the cycling without feedback and produced 3.75 K·h outside the band over three days.
In the final structure the plan sets the PI's set point, and the PI (60 s) does the tracking and
the cycling. The planner keeps 0.6 K inside the occupied band for the PI's ripple; the margin
was chosen from a sweep (0.3 K: 0.46 K·h outside the band, 0.6 K: 0.01 K·h, for 4 % more energy).

![MPC over one day](figures/mpc_day.png)

Three days, ventilation during occupied hours, tariff 0.15 / 0.60 EUR/kWh with the peak 14-20 h:

| Controller | Energy (kWh) | Cost (EUR) | K·h outside band | Starts |
|---|---|---|---|---|
| PI, 25 °C constant | 6.21 | 2.15 | 0.000 | 123 |
| PI, set-back, fixed start 6:00 | 5.62 | 1.91 | 0.000 | 39 |
| MPC | 5.45 | 1.95 | 0.009 | 61 |
| MPC, time-of-use tariff | 6.37 | **1.72** | 0.000 | 89 |
| MPC, model parameters off by 20 % | 5.46 | 1.94 | 0.001 | 64 |
| MPC, forecast +2 K and -30 % sun | 5.44 | 1.95 | 0.010 | 60 |

- **Against a well-configured PI the energy saving is small**: 3 % against set-back with a fixed
  early start, 12 % against a constant set point. Most of what the MPC finds (let the room float
  at night, start just in time, sit at the warm edge of the band) a PI with a good schedule
  already gets.
- **With a time-of-use tariff the MPC cuts the cost by 10 %** while using 13 % more energy: it
  cools the room down to 20 °C overnight, stores the cold in the thermal mass, and runs at lower
  power through the day and the peak. A PI has no way to do this.
- **Load shifting has a threshold.** Planning studies with different price ratios show no
  shifting below a peak/off-peak ratio of about 3; above it, at most 15-23 % of the cooling in the
  peak window moves out of it, depending on the width of the band. Pre-cooling the air pulls heat
  out of the strongly coupled mass (180 W/K) and much of it leaks back before the peak.
- **Robustness**: 20 % errors in five model parameters, or a forecast 2 K too warm with 30 % less
  sun, change energy and comfort by less than 1 %. The Kalman disturbance state absorbs both.

## Verification

All checks run in `pytest` (138 tests, about a minute and a half).

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
| Sigma-delta: mean delivered demand, cycle length vs design formula | ±0.01, ±10 % |
| Cascade inner loop on the real coil: coil temperature vs its set point | ±0.05 K |
| Control model: exact discretisation; fixed point = continuous steady state | 10⁻¹² |
| Control model vs plant with the unit off, one day | < 0.01 K |
| Kalman filter: mass temperature and unknown 250 W heat input after 6 days | 10⁻³ K, 0.1 % |
| MPC plan: feasible, inside the planning band, no cooling while the band is wide | always |

![Convergence](figures/convergence.png)

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

![Energy balance](figures/energy_balance.png)

![Humidity](figures/humidity.png)

## Layout

```
src/roomsim/
  psychro.py       psychrometric functions
  params.py        RoomParams dataclass, derived constants
  inputs.py        disturbance scenarios d(t), vectorised for forecasts
  plant.py         rhs(t, x, u, d, p), flux breakdown, Jacobian, sensors
  ac_unit.py       inverter split unit: coil, cycle, frost protection, start-up, coil water
  integrators.py   Euler, RK4, Euler stability limit
  simulate.py      fixed-step loop, discrete controller with zero-order hold, periodic shooting
  results.py       SimResult: states, psychrometrics, balance residuals
  metrics.py       energy, starts, comfort, humidity and CO2 indicators
  verification.py  analytical solutions and convergence studies (shared by tests and figures)
  plotting.py      figures
  control_model.py linear RC model for the controllers, exact discretisation
  estimation.py    Kalman filter with a disturbance state (offset-free)
  tariff.py        time-of-use electricity price
  controllers/     on/off thermostat, discrete PI with anti-windup, low-load handlers,
                   humidity-aware fan control (direct and cascade), economic MPC
scripts/run_baseline.py         free-floating baseline and verification figures
scripts/compare_controllers.py  controller comparison table and figures
scripts/compare_humidity.py     humidity-aware fan control and robustness to the unit
scripts/compare_mpc.py          MPC against PI, tariff and robustness cases
tests/
docs/specification.md
```

## Install and run

```bash
pip install -e .[dev]
pytest -v
python scripts/run_baseline.py
python scripts/compare_controllers.py
python scripts/compare_humidity.py
python scripts/compare_mpc.py
```

## Roadmap

1. Ventilation-level strategies driven by CO2.
2. MPC with humidity in the model, so that the planner can trade sensible against latent load.
