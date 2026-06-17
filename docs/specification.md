# Model specification

Specification of the room plant. Where the implementation showed that an earlier formulation
was physically inconsistent, it was revised; each revision is marked **[rev]** with the reason.

## 0. Structure

| Module | Content |
|---|---|
| `psychro.py` | psychrometric functions |
| `params.py` | room parameters in a frozen dataclass, derived quantities as properties |
| `inputs.py` | disturbances d(t): outdoor weather, sun, occupancy |
| `plant.py` | `rhs(t, x, u, d, p)` returning dx/dt (continuous plant) and `fluxes(...)` with every term |
| `integrators.py` | fixed-step explicit Euler and classical fourth-order Runge-Kutta (RK4) |
| `simulate.py` | fixed-step loop coordinating plant, controller and logging |
| `controllers/` | **discrete** controllers executed every Δt_c, output held (zero-order hold) |
| `tests/` | verification of section 5 |

Modules added later (air-conditioning unit, controllers, estimation, ventilation fan) are listed
in the layout section of the [README](../README.md#layout).

**State** x = [T_a, T_m, W, φ_b, C] **[rev]**
- T_a: indoor air temperature, °C
- T_m: thermal-mass temperature (walls, furniture), °C
- W: indoor humidity ratio, kg water / kg dry air
- φ_b: relative humidity in equilibrium with the moisture-buffer layer of the interior
  surfaces, fraction
- C: indoor CO2 concentration, ppm

**Control input** u = [Q_hvac_sens, ṁ_hvac_lat, V_vent]
- Q_hvac_sens: sensible heat delivered by the unit, W (negative when cooling)
- ṁ_hvac_lat: water removed by the unit, kg/s
- V_vent: mechanical outdoor-air supply, m³/s (balanced, no heat recovery)
- Zero in the baseline.

## 1. Thermal balance: two-node RC model

```
Air:   C_a·dT_a/dt = H_am·(T_m − T_a) + (H_ao + H_inf + H_vent)·(T_o − T_a)
                     + f_sol,a·Q_sol + Q_int,conv + Q_hvac_sens
Mass:  C_m·dT_m/dt = H_am·(T_a − T_m) + H_mo·(T_o − T_m) + (1 − f_sol,a)·Q_sol + Q_int,rad
```

with H_inf = ρ·c_p·V·n/3600, H_vent = ρ·c_p·V_vent, H_am = 1/R_am = h_i·A_int and
Q_sol = g·A_win·G_sol.

## 2. Moisture balance **[rev]**

```
ρ·V·dW/dt   = ρ·(V·n/3600 + V_vent)·(W_o − W) + N_pers·g_vap − ṁ_hvac_lat + ṁ_buf
C_b·dφ_b/dt = −ṁ_buf
ṁ_buf       = G_b·(φ_b·p_ws(T_m) − p_w(W))
```

*Reason:* the original balance stored moisture in the room air only. Interior surfaces
(painted gypsum) hold about as much exchangeable moisture over a day as the air itself, so the
air-only model makes humidity respond twice as fast as it does in a real room, which would flatter
any humidity controller. The buffer is the surface layer of the Effective Moisture Penetration
Depth (EMPD) model as documented in the EnergyPlus Engineering Reference, without its deep layer
and coating ([references.md](references.md)):

- C_b = ρ_mat·d_EMPD·A·(du/dφ): water stored per unit RH of the layer.
- G_b = A / (1/β + (d_EMPD/2)/δ): surface resistance in series with diffusion through half the
  layer.
- β from the Lewis analogy: β = h_i/(ρ·c_p)/(R_v·T).
- The layer is at the mass temperature T_m, so a warming wall releases moisture.
- The heat of sorption is neglected (order 10 W for daily cycles).
- Can be switched off with `RoomParams(moisture_buffer=False)`, which recovers the original
  model.

## 2b. CO2 balance

```
V·dC/dt = (V·n/3600 + V_vent)·(C_o − C) + N_pers·G_CO2·10⁶
```

with C in ppm, C_o = 420 ppm and G_CO2 = 0.0052 L/s = 5.2·10⁻⁶ m³/s per person (sedentary
office, 1.2 met); in the equation G_CO2 is in m³/s. The value is the one for an adult at 1.2 met
cited in ASHRAE Standard 62.1. Persily and de Jonge (2017) estimate a lower average for office
occupants, 0.0048 L/s, so the CO2 rise here is on the high side.
With infiltration alone the occupied steady state is 1 807 ppm, so the room needs mechanical
ventilation; V_vent is the actuator for ventilation-level strategies.

## 2c. Air-conditioning unit (inverter split, cooling)

Commands from the controller: compressor speed s (0 = off, otherwise [0.3, 1]), indoor fan f
([0.4, 1]) and ventilation V_vent. The unit turns them into the plant inputs
Q_hvac_sens = −Q_sens + P_fan,in and ṁ_hvac_lat at every Runge-Kutta stage, using the return-air
state of that stage.

Abbreviations: BF, bypass factor; NTU, number of transfer units; UA, heat-exchanger
conductance; COP, coefficient of performance of the compressor; EER, energy efficiency ratio of
the whole unit, net of the indoor fan heat; SHR, sensible heat ratio.

| Item | Model |
|---|---|
| Coil **[rev]** | refrigerant at uniform T_evap; air-side conductance (BF = exp(−NTU), NTU ∝ ṁ^−0.4, BF = 0.15 at 450 m³/h) in series with UA_evap = 600 W/K |
| Partially wet **[rev]** | dry from the air inlet to where the surface reaches the inlet dew point, wet after it; closed forms in both parts; wetted fraction continuous from 0 to 1 (D29) |
| Split | Q_sens = ṁ(c_pa + c_pv W_out)(T_in − T_out), Q_lat = ṁ_lat(h_fg0 + c_pv T_in): their sum equals ṁ Δh exactly |
| Outlet | brought to saturation at constant enthalpy if the mixing line would cross it |
| Capacity **[rev]** | s × 2.5 kW × [1 + 0.040 (T_evap − T_evap,rated)] × [1 − 0.016 (T_o − 35 °C)], solved with the coil for T_evap (D30) |
| Frost protection | evaporating temperature held at ≥ 0.5 °C; the capacity is what the coil takes at that temperature |
| Cycle | COP = η(s) T_evap/(T_cond − T_evap) with absolute temperatures, T_cond = max(T_o + Q_cond/UA_cond, T_evap + 15 K) with Q_cond = Q_coil (1 + 1/COP), closed form; η calibrated to a net EER of 3.5 (fans included, indoor fan heat subtracted from the capacity) at 35 °C / 27 °C / W = 0.0114 **[rev]** |
| Indoor fan **[rev]** | 5 W + 20 W × f³ (D31) |
| Start-up | delivered capacity × (1 − e^(−t_on/60 s)); compressor power that of the steady operating point at the commanded speed. Choice of 60 s: decision D9 in [decisions.md](decisions.md) |

At rating conditions: SHR 0.73, mean coil surface 13.4 °C, supply 15.1 °C. Lower fan speed gives
a colder coil, more dehumidification, less capacity and a lower COP; at 30 % speed under rating
conditions the net EER is 3.8.

*Reason for the revisions [rev]:* an external review pointed out that the uniform-surface coil
was either completely dry or completely wet, that the capacity did not depend on airflow, that
the fan power at low speed was too small, and that the EER counted the indoor fan heat as
cooling. Decisions D29-D32 in [decisions.md](decisions.md).

**Coupling accuracy [rev].** With the unit evaluated at the stage states, RK4 keeps fourth order:
measured 4.05 with a constant command and 4.02 in closed loop with the PI sampled every 60 s, as
long as the coil stays in one regime and no switching decision changes with the step. A
trajectory that crosses from a fully wet into a partially wet coil gives order 2.10 (the coil
output has a kink there), with an error of 4·10⁻¹⁰ K at Δt = 60 s. Evaluating the unit once per
sample and holding Q_sens gives first order (1.13): at Δt = 60 s the error is 1·10⁻³ K instead of
5·10⁻¹¹ K. The first version of this paragraph said "the closed loop keeps fourth order
(measured 3.9)"; that measurement used a constant command (decision D8).

## 3. Psychrometrics (p = 101 325 Pa)

- Saturation pressure, Magnus form with the coefficients of Alduchov and Eskridge, recommended
  by them for water between −40 and 50 °C (T in °C):
  p_ws = 610.94·exp(17.625·T/(T + 243.04)) Pa
- Humidity ratio: W = 0.622·p_w/(p − p_w); inverse p_w = W·p/(0.622 + W)
- Relative humidity: RH = p_w/p_ws(T)
- Dew point: T_d = 243.04·ln(p_w/610.94)/(17.625 − ln(p_w/610.94))
- Enthalpy: h = 1.006·T + W·(2501 + 1.86·T), kJ/kg dry air

## 4. Parameters

Room 5 × 4 × 2.7 m, office for two people.

| Quantity | Value |
|---|---|
| V | 54 m³ |
| C_a | ρ·c_p·V·5 ≈ 0.33 MJ/K (factor 5 for furnishings) |
| C_m | 10 MJ/K |
| h_i·A_int | 3 W/(m²·K) · 60 m² = 180 W/K |
| H_ao (window) | 3 m² · 1.3 W/(m²·K) ≈ 4 W/K |
| H_mo (walls) | 25 W/K |
| n | 0.5 h⁻¹ → H_inf ≈ 9 W/K |
| Person | 70 W sensible (half convective, half radiant), 0.05 kg/h vapour |
| Window | 3 m², g = 0.6, f_sol,a = 0.3 |
| Moisture buffer **[rev]** | gypsum, ρ = 800 kg/m³, d_EMPD = 3 mm, du/dφ = 0.01, δ = 2.5·10⁻¹¹ kg/(m·s·Pa), A = 60 m² → C_b = 1.44 kg, G_b ≈ 5.2·10⁻⁷ kg/(s·Pa) |

**Summer scenario, 3 days** (preceded by warm-up days, see below)
- T_o = 26 + 6·sin(2π(t/86400 − 9/24)) °C, maximum at 15 h.
- Outdoor humidity **[rev]**: W_o constant, equal to RH 55 % at the mean temperature 26 °C
  (W_o = 0.0115 kg/kg). *Reason:* holding RH_o = 55 % constant and converting at the
  instantaneous T_o makes W_o double between 3 h and 15 h (0.008 → 0.0165 kg/kg), which no real
  air mass does and which overstates the afternoon latent load. The original behaviour remains
  available with `SummerScenario(constant_W=False)`.
- Sun on the window: half-sine between 6 h and 18 h, peak 500 W/m².
- Occupancy: 2 people from 8 h to 18 h.

**Initial state and warm-up [rev]:** the mass time constant is about 3.1 days, as long as the
simulated period, so the result depends on the initial state. The reported run is preceded by
warm-up days (discarded) until the daily cycle is periodic.

**Time scales:** τ_air ≈ 1 700 s; τ_mass ≈ 3.1 days; air moisture alone 3 600/n = 7 200 s;
air-buffer exchange ≈ 400 s (coupled mode); apparent moisture time constant with buffering
≈ 15 000 s.
The system is stiff. RK4 with 60 s is accurate; explicit Euler is stable for steps below
2/|λ_max|, checked numerically from the Jacobian. **[rev]** The fastest mode is the air-buffer
moisture exchange, not the air temperature, and it stiffens with temperature because the
buffer flux scales with p_ws(T_m): the Euler limit is ≈ 790 s at 26 °C and ≈ 380 s at 48 °C
(RK4: ≈ 1 100 s and ≈ 530 s).

## 5. Verification

1. **Psychrometrics:** p_ws(20 °C) ≈ 2.33 kPa; W(25 °C, 50 %) ≈ 0.0099 kg/kg;
   T_d(20 °C, 50 %) ≈ 9.3 °C.
2. **Steady state:** with constant inputs, the simulation matches the analytical solution of the
   2×2 linear thermal system.
3. **One-node limit:** a step in outdoor temperature against the analytical exponential.
   Convergence order ≈ 4 with RK4 and ≈ 1 with Euler.
4. **Energy-balance closure [rev]:** residual C_a·ΔT_a + C_m·ΔT_m − ∫Σ boundary fluxes dt,
   checked two ways:
   - (a) integrating the fluxes as extra states with the same integrator: closes to round-off
     (every Runge-Kutta method preserves linear invariants), which verifies signs and bookkeeping;
   - (b) integrating the logged fluxes with composite Simpson: the residual falls as O(Δt) with
     Euler and O(Δt⁴) with RK4, which verifies the order.
   The same two checks apply to the water balance of air + buffer.
5. **Moisture without sources [rev]:**
   - buffer off: W tends to W_o with time constant 3 600/n;
   - buffer on, isothermal: W and φ_b follow the analytical solution of the linearised 2×2 system
     (one fast and one slow exponential).

## 6. Baseline deliverable

- 3-day simulation without air conditioning.
- Figures: temperatures (T_a, T_m, T_o); relative humidity with dew point; energy-balance closure.
- All tests passing.
- Maximum room temperature recorded: free-floating baseline for the controllers.
