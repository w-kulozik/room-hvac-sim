# roomsim: a dynamic room model for HVAC control

A small simulator of a single room, verified against analytical solutions, used as the plant
for testing air-conditioning controllers (on/off, PI, cascade, model-predictive).

The room has an air node, a thermal-mass node, a moisture balance with buffering by the
interior surfaces, and proper psychrometrics. It is driven by outdoor weather, sun through
a window and occupants coming and going. Controllers are discrete: they sample the sensors
every Δt_c and hold their output, as a real unit does.

![Free-floating temperatures](figures/temperatures.png)

## Model

State x = [T_a, T_m, W, φ_b], control u = [Q_hvac_sens, ṁ_hvac_lat].

```
C_a dT_a/dt  = H_am (T_m − T_a) + (H_ao + H_inf)(T_o − T_a) + f_sol,a Q_sol + Q_int,conv + Q_hvac
C_m dT_m/dt  = H_am (T_a − T_m) + H_mo (T_o − T_m) + (1 − f_sol,a) Q_sol + Q_int,rad
ρV  dW/dt    = ṁ_inf (W_o − W) + N g_vap − ṁ_hvac_lat + ṁ_buf
C_b dφ_b/dt  = −ṁ_buf,      ṁ_buf = G_b (φ_b p_ws(T_m) − p_w(W))
```

| Symbol | Meaning | Value |
|---|---|---|
| C_a | air + light furnishings | 0.33 MJ/K |
| C_m | walls, floor, heavy furniture | 10 MJ/K |
| H_am | air ↔ mass, h_i·A_int | 180 W/K |
| H_ao, H_mo, H_inf | window, opaque envelope, infiltration (0.5 ach) | 4, 25, 9.05 W/K |
| C_b, G_b | moisture buffer (EMPD, painted gypsum) | 1.44 kg, 5.2·10⁻⁷ kg/(s·Pa) |

Room 5 × 4 × 2.7 m, two occupants 8-18 h (70 W, 0.05 kg/h each), 3 m² window (g = 0.6),
sinusoidal summer day 20-32 °C. Full specification, including where and why it departs from
its first version: [docs/specification.md](docs/specification.md).

Time scales: air ≈ 1 700 s, mass ≈ 3.1 days, air-buffer moisture exchange ≈ 900 s,
apparent room moisture ≈ 15 000 s (twice the 7 200 s of the air alone).

## Day-1 result: free-floating summer baseline

Three days on the **periodic steady state**, found by Newton shooting on the one-day map.
The ~3-day mass time constant makes plain warm-up slow, and starting from a uniform 26 °C
underestimates the peak by almost 4 K on day 3:

| Start | T_a max day 1 | day 2 | day 3 |
|---|---|---|---|
| periodic steady state | **36.76 °C** (15:07) | 36.76 °C | 36.76 °C |
| uniform 26 °C (naive) | 29.65 °C | 31.55 °C | 32.95 °C |

**Baseline for the controllers: T_a,max = 36.8 °C at 15:07**, T_a min 33.4 °C, mean mass
temperature 35.1 °C, indoor RH 35-39 %. With 0.5 ach, an unshaded window and about 340 W of
mean gains, the room sits ~9 K above the mean outdoor temperature.

![Humidity](figures/humidity.png)

The buffer layer (dashed) absorbs moisture while the room is occupied and releases it at night,
so the indoor dew point rises and falls more slowly than the air alone would allow.

## Verification

All checks run in `pytest` (54 tests, ~16 s).

| Check | Result |
|---|---|
| Psychrometrics: p_ws(20 °C), W(25 °C, 50 %), T_d(20 °C, 50 %) | 2333 Pa, 0.00986 kg/kg, 9.26 °C |
| Steady state vs 2×2 linear solution (3 cases) | < 10⁻⁶ K |
| Daily mean of periodic solution vs steady state of mean inputs | < 2·10⁻³ K |
| One-node step vs exponential, convergence order | RK4 **4.04**, Euler **1.01** |
| Energy and water closure, flows integrated as extra states | round-off (~10⁻¹³ relative) |
| Energy residual with independent Simpson quadrature, order | RK4 **4.16**, Euler **1.05** |
| Moisture decay without buffer, fitted τ vs 3600/n | 10⁻⁸ relative |
| Moisture decay with buffer vs linearised two-exponential solution | < 10⁻³ of the perturbation |
| Full scenario vs `scipy` Radau (rtol 10⁻¹¹) | < 10⁻⁶ K |

![Convergence](figures/convergence.png)

The energy balance is checked at two levels. Any Runge-Kutta method preserves
linear invariants, so integrating the boundary flows alongside the state closes the balance to
round-off whatever the step: that verifies signs and bookkeeping. Integrating the logged flows
with an independent quadrature (Simpson) gives a residual that falls with the order of the
method: that verifies accuracy.

![Energy balance](figures/energy_balance.png)

**Stiffness.** The fastest mode is the air-buffer moisture exchange, and it stiffens with
temperature because the buffer flux scales with p_ws(T_m). The explicit-Euler stability limit is
≈ 790 s at 26 °C and ≈ 380 s at 48 °C; RK4 at 60 s has ample margin.

## Layout

```
src/roomsim/
  psychro.py       psychrometric functions
  params.py        RoomParams dataclass, derived constants
  inputs.py        disturbance scenarios d(t), vectorised for forecasts
  plant.py         rhs(t, x, u, d, p), flux breakdown, Jacobian, sensors
  integrators.py   Euler, RK4, Euler stability limit
  simulate.py      fixed-step loop, discrete controller with zero-order hold, periodic shooting
  results.py       SimResult: states, psychrometrics, balance residuals
  verification.py  analytical solutions and convergence studies (shared by tests and figures)
  plotting.py      figures
  controllers/     discrete controllers (in progress)
scripts/run_baseline.py
tests/
docs/specification.md
```

## Install and run

```bash
pip install -e .[dev]
pytest -v
python scripts/run_baseline.py
```

## Roadmap

1. Air-conditioning unit model with capacity limits, sensible/latent split and COP.
2. On/off thermostat with hysteresis; discrete PI with anti-windup.
3. Cascade temperature and humidity control.
4. Linear MPC with a Kalman filter estimating the unmeasured mass temperature, tested under
   ±20 % model mismatch.
5. Common metrics: electrical energy, degree-hours outside the comfort band, hours above
   60 % RH, compressor starts.
