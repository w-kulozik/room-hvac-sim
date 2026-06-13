# Model

The room model and the air-conditioning unit. The full specification, with each revision and
its reason, is in [specification.md](specification.md).

## Room

State x = [T_a, T_m, W, φ_b, C] for the room, plus the water M_c held on the coil when a unit
is present. Plant input u = [Q_hvac_sens, ṁ_hvac_lat, V_vent].

- T_a, T_m: air and thermal-mass temperatures, °C.
- W: humidity ratio of the room air, kg of water per kg of dry air.
- φ_b: relative humidity in equilibrium with the moisture-buffer layer of the interior
  surfaces.
- C: CO2 concentration, ppm.
- Q_hvac_sens: sensible heat delivered by the unit, W (negative when cooling, indoor fan heat
  included); ṁ_hvac_lat: water removed by the unit, kg/s; V_vent: mechanical outdoor-air
  supply, m³/s.

```
C_a dT_a/dt  = H_am (T_m - T_a) + (H_ao + H_inf + H_vent)(T_o - T_a)
               + f_sol,a Q_sol + Q_int,conv + Q_hvac_sens
C_m dT_m/dt  = H_am (T_a - T_m) + H_mo (T_o - T_m) + (1 - f_sol,a) Q_sol + Q_int,rad
ρV  dW/dt    = ρ(V n/3600 + V_vent)(W_o - W) + N g_vap - ṁ_hvac_lat + ṁ_buf
C_b dφ_b/dt  = -ṁ_buf,      ṁ_buf = G_b (φ_b p_ws(T_m) - p_w(W))
V   dC/dt    = (V n/3600 + V_vent)(C_o - C) + 10⁶ N G_CO2
```

Here n is the infiltration rate in air changes per hour (ach), H_vent = ρ c_p V_vent, N the
number of occupants, and G_CO2 is in m³/s per person in the CO2 balance.

| Symbol | Meaning | Value |
|---|---|---|
| C_a | air + light furnishings | 0.33 MJ/K |
| C_m | walls, floor, heavy furniture | 10 MJ/K |
| H_am | air to mass, h_i·A_int | 180 W/K |
| H_ao, H_mo, H_inf | window, opaque envelope, infiltration (0.5 ach) | 4, 25, 9.05 W/K |
| C_b, G_b | moisture buffer (EMPD, painted gypsum) | 1.44 kg, 5.2·10⁻⁷ kg/(s·Pa) |
| G_CO2 | CO2 per person (1.2 met) | 0.0052 L/s |

EMPD is the effective moisture penetration depth model of surface moisture buffering; the
buffer here is the surface layer of the EMPD model in the EnergyPlus Engineering Reference,
without its deep layer. G_CO2 is the value for an adult at 1.2 met cited in ASHRAE Standard
62.1. Both are referenced in [references.md](references.md).

Room 5 × 4 × 2.7 m, two occupants 8-18 h (70 W, 0.05 kg/h each), 3 m² window (g = 0.6),
sinusoidal summer day 20-32 °C, outdoor humidity ratio constant at 0.0115 kg/kg.

Time scales: air ≈ 1 700 s, mass ≈ 3.1 days, air-buffer moisture exchange ≈ 400 s,
apparent room moisture ≈ 15 000 s (twice the 7 200 s of the air alone).

## Air-conditioning unit

Inverter split unit, 2.5 kW at 35 °C outdoor, energy efficiency ratio (EER) 3.5 at rating
conditions. The controller sends compressor speed s (0 or 0.3-1), indoor fan speed f and
ventilation flow; the unit converts them into Q_hvac_sens and ṁ_hvac_lat at every Runge-Kutta
stage.

| Item | Model |
|---|---|
| Coil | bypass factor / apparatus dew point, BF = exp(−NTU) = 0.15 at 450 m³/h, NTU ∝ ṁ^-0.4 |
| Dry / wet coil | closed form when dry; Newton on the saturation enthalpy when wet; continuous at the switch |
| Sensible / latent split | exact for the enthalpy definition: Q_sens + Q_lat = ṁ Δh to round-off |
| Capacity | 2.5 kW at 35 °C outdoor, changing by 1 % per K of outdoor temperature (lower when warmer) |
| Cycle | Carnot fraction with finite-UA evaporator and condenser, closed form, calibrated to the rated EER; condensing temperature at least 15 K above evaporating (head-pressure control) |
| Frost protection | capacity capped so the evaporating temperature stays ≥ 0.5 °C |
| Start-up | capacity ramps with a 60 s time constant while the compressor already draws the power of the steady operating point |
| Coil water | condensate wets the coil (0.25 kg) before it drains; re-evaporates with the fan on and the compressor off |

NTU is the number of transfer units of the coil, UA a heat-exchanger conductance, and the
coefficient of performance (COP) the cooling per unit of compressor power.

At rating conditions the sensible heat ratio is 0.73, the coil surface 13 °C and the supply
air 15 °C. Lower fan speed gives a colder coil, more dehumidification and a lower COP; part load
raises the EER to 4.0 at minimum speed.

The ventilation fan is modelled separately with the affinity law P = P_d (V/V_d)³, 30 W at
108 m³/h (decision D23 in [decisions.md](decisions.md)).
