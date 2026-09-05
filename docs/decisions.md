# Decision log

Every modelling and design decision that departs from the original specification or that was
revised during the work, with the evidence behind it.

**Timing** says when the decision was taken:
- *a priori*: before seeing the result it affects;
- *after observation*: changed because a result looked wrong. These carry the larger risk of
  tuning to the answer, so the evidence and the alternatives are spelled out.

Where a test was changed after it failed, the entry says whether the model or the test was wrong.

Abbreviations: RH, relative humidity; EMPD, effective moisture penetration depth; COP,
coefficient of performance; EER, energy efficiency ratio; SHR, sensible heat ratio; PI,
proportional-integral controller; MPC, model-predictive control; SIMC, Skogestad's tuning rules
for PI controllers.

Publications named in this log are listed with full references in
[references.md](references.md).

---

## Plant

### D1. Outdoor humidity: constant humidity ratio, not constant relative humidity
- **Timing:** after observation (first run of the scenario).
- **Observation:** holding RH_o = 55 % and converting at the instantaneous T_o makes W_o go from
  0.0080 to 0.0165 kg/kg within a day. Real air masses keep W roughly constant over a day while
  RH swings.
- **Decision:** W_o = W(26 °C, 55 %) = 0.0115 kg/kg, constant. The original behaviour remains as
  `SummerScenario(constant_W=False)`.
- **Effect:** the afternoon latent load is no longer doubled.

### D2. Moisture buffering of interior surfaces (EMPD)
- **Timing:** a priori. Raised when the moisture balance was written; the specification stored
  moisture in the air only.
- **Evidence:** painted gypsum (800 kg/m³, 3 mm penetration depth, du/dφ = 0.01, 60 m²) holds
  about as much exchangeable water as the room air: 72 kg against 65 kg of air, in humidity-ratio
  units. The apparent moisture time constant goes from 7 200 s to about 15 000 s.
- **Alternatives:** a constant multiplier on the air's moisture capacity. Rejected because it
  misses the exchange dynamics (≈ 400 s; first estimated at 900 s, corrected in D22) and the
  release of moisture by warming walls.
- **Limits:** heat of sorption neglected (order 10 W in daily cycles).

### D3. Periodic steady state by Newton shooting
- **Timing:** after observation.
- **Observation:** starting at 26 °C, the daily maximum still rose by 0.15 K/day after seven
  warm-up days (mass time constant ≈ 3.1 days).
- **Decision:** solve x = Φ_24h(x) by Newton with a finite-difference Jacobian. It converges in
  two iterations because the thermal part is linear.
- **Check:** the daily mean of the periodic solution equals the steady state under mean inputs
  (exact for a linear time-invariant system), to < 2·10⁻³ K.

### D4. Piecewise-constant inputs sampled at the step start
- **Timing:** a priori.
- **Reason:** the last RK4 stage of the step ending at 8:00 would otherwise already see the
  occupants. A test checks that nothing from 8:00 leaks into [7:00, 8:00).

### D5. Mechanical ventilation as a third plant input
- **Timing:** a priori. Added together with CO2.
- **Evidence:** with infiltration alone, two occupants take the room to 1 807 ppm of CO2.

---

## Air-conditioning unit

### D6. Bypass-factor coil with an exact sensible/latent split
- **Timing:** a priori.
- **Decision:** the split follows from the enthalpy definition, so Q_sens + Q_lat = ṁ Δh to
  round-off; the dry and wet branches meet continuously at the dew point.
- **Later change:** the all-dry or all-wet coil was replaced by a partially wet one (D29). The
  exact split is kept.

### D7. Frost protection
- **Timing:** after observation.
- **Observation:** at 40 % fan and full speed the evaporating temperature reached −2.8 °C.
- **Decision:** cap the capacity so that T_evap ≥ 0.5 °C, as real units do.
- **Later change:** since the capacity falls with the evaporating temperature (D30), the same
  case gives 3.9 °C and the limit is no longer reached at rating conditions. It is reached with
  cool, dry return air (for example 22 °C and 0.008 kg/kg, full speed, 40 % fan); a test checks
  that case. The cycle efficiency in that regime was corrected in D40.

### D8. Unit evaluated at every Runge-Kutta stage
- **Timing:** a priori, then confirmed by a test.
- **Evidence (first version):** measured convergence order 3.9 with a constant compressor
  command. Evaluating once per sample and holding Q_sens gave 1.1, and at Δt = 60 s an error of
  2·10⁻⁴ K instead of 1·10⁻¹⁰ K. The earlier summary "fourth order, also in closed loop"
  overstated this: the test used a constant command, not a controller.
- **Evidence (after D29-D32):** four tests in `tests/test_unit_in_loop.py`, comparing the final
  air temperature with a reference run with a step of 7.5 s (3.75 s for the second and third):
  - constant command, coil dry throughout: order 4.05; error 5·10⁻¹¹ K at Δt = 60 s;
  - closed loop with the PI (sigma-delta handler) sampled every 60 s while the integration step
    goes from 60 to 15 s: order 4.02; error 1·10⁻¹⁰ K at Δt = 60 s. In that run the compressor
    modulates throughout and the coil stays dry, so no switching decision and no coil regime
    changes with the step;
  - constant command, trajectory crossing from a fully wet into a partially wet coil: order
    2.10. The coil output is continuous there but its derivative is not, and a fixed-step method
    loses order on the step that contains the crossing. Error 4·10⁻¹⁰ K at Δt = 60 s;
  - unit evaluated once per sample and Q_sens held: order 1.13; error 1·10⁻³ K at Δt = 60 s.
- **Not covered:** runs in which a compressor start or stop, or a fan step, happens at a
  different sample for different integration steps. Then the runs diverge by a whole switching
  event and no convergence order exists; this is a property of the sampled controller, not of
  the integrator.

### D9. Start-up time constant 60 s (initially 120 s)
- **Timing:** after observation.
- **Observation:** with 120 s, about 30 % of the on/off thermostat's energy went into start-up
  ramps (average on time 6.5 min).
- **Evidence for the new value:** Henderson, Huang and Parker (1999) give 40-60 s as typical of
  modern residential air conditioners and heat pumps, and use 60 s for their "typical" unit
  from literature values of 30-80 s. Those are ducted residential units, not small split units,
  so the value is a plausible order of magnitude, not a property of the modelled unit. Full
  citations are in [references.md](references.md).
- **Risk:** the value was revisited because a result looked large. It stays a parameter
  (`ACParams.tau_start`). On/off energy is the result most sensitive to it.

### D10. Water held on the coil
- **Timing:** a priori.
- **Mechanism:** the one described by Henderson and Rengarajan (1996): condensate held on the
  coil re-evaporates into the air stream when the compressor stops and the fan keeps running.
  Here the held water is a state of the simulation; their part-load equations are not used.
- **Parameters:** M_max = 0.25 kg, which fills in about 16 min at rated latent load. Not
  calibrated against a specific unit; the strength of the part-load latent degradation scales
  with it.
- **Limits:** film enthalpy neglected in fan-only operation; no natural drying with the fan off.

### D11. Drained and external water integrated as states
- **Timing:** after observation.
- **Observation:** a conservation test using the trapezoid rule on logged rates failed at every
  compressor switch, where the rates jump.
- **Verdict:** the test was wrong, not the model. Integrating these quantities alongside the
  state closes the water balance to round-off, consistent with the energy check.

---

## Controllers

### D12. Anti-windup tracking time T_t = T_i/2 (initially T_i)
- **Timing:** after observation.
- **Observation:** a test expected fast recovery after saturation and failed. Under saturation
  the integrator settles at I = u_max − K_p e (1 − T_t/T_i), so with T_t = T_i and a large error
  it sits at the limit and recovers only through integral action.
- **Verdict:** the test expectation was wrong.
- **Decision:** T_t = T_i/2. This is a choice made here, not a value from the literature; the
  test now checks the analytical settling value for three values of T_t.
- **Limits:** Åström and Murray (2008, section 10.4) suggest a reset gain that is a fraction of
  1/T_i, that is T_t longer than T_i, because a short tracking time lets measurement noise reset
  the integrator. The simulated sensor has no noise, so this risk does not appear here; on a
  real unit T_t would need checking against sensor noise.

### D13. Low-load handling separate from the PI
- **Timing:** a priori.
- **Reason:** the compressor cannot run below 30 % speed. The PI computes a continuous demand and
  a handler turns it into feasible commands:
  - A: thermostatic hysteresis with a frozen integrator;
  - C: sigma-delta, which meets the mean demand.

### D14. Sigma-delta: clamped accumulator and time scale from the ripple (initially T_ref = 600 s)
- **Timing:** after observation, while investigating the MPC's comfort violations.
- **Observations:**
  1. After a long run the accumulator owed off time and kept the unit off while the room was
     already above its set point.
  2. T_ref = 600 s had no physical basis; the resulting ripple was ≈ 0.9-1.1 K.
- **Decision:** clamp the accumulator to ± threshold, and set T_ref = ΔT_pp C_a / (s_min Q_max),
  with ΔT_pp = 0.5 K (T_ref ≈ 320 s). The capacity in this formula was later corrected (D21,
  T_ref = 210 s).
- **Effect on all PI-based results, at the time:** maximum occupied temperature 26.0 → 25.7 °C;
  starts per three days 84 → 123. Tables were regenerated. The current values, after D20 and
  D21, are in [results.md](results.md).

### D15. Humidity through the indoor fan; evaporator thermistor as a measurement
- **Timing:** a priori.
- **Evidence:** at minimum speed (room at 25 °C and 60 % RH, 32 °C outdoors), going from full to
  40 % fan takes the coil from 19 °C (dry) to 15 °C (wet), SHR 1.00 → 0.75, for a 5 % loss of
  EER.
- **After D29-D32:** the mean coil surface goes from 17.7 °C (dry) to 14.4 °C (fully wet), SHR
  1.00 → 0.70, with 10 % less capacity and a 3 % loss of net EER. At 70 % fan the coil is
  already 74 % wet.

### D16. MPC supervises the PI instead of driving the compressor
- **Timing:** after observation.
- **Observation:** applying the planned cooling directly, with the low-load handler cycling
  between the 15-minute solves, gave 3.75 K·h outside the comfort band in three days.
- **Decision:** the plan sets the PI's set point; the PI tracks it every 60 s with feedback.

### D17. Comfort band edges belong to the occupied period
- **Timing:** after observation. This was a bug.
- **Observation:** the band at 18:00 was evaluated as unoccupied, so the planner let the room
  warm during the last occupied quarter hour; the same at 8:00.

### D18. MPC comfort margin 0.6 K
- **Timing:** after observation. This is tuning, chosen by a sweep.
- **Evidence:**

  | Margin | Energy (kWh) | K·h outside the band |
  |---|---|---|
  | 0.3 K | 5.24 | 0.46 |
  | 0.4 K | 5.31 | 0.20 |
  | 0.5 K | 5.38 | 0.06 |
  | 0.6 K | 5.45 | 0.01 |

  This sweep was run before the corrections D20 and D21, so its energies are lower than the
  later ones. Repeated after D20 and D21 (same three days after three warm-up days):

  | Margin | Energy (kWh) | K·h outside the band |
  |---|---|---|
  | 0.3 K | 5.58 | 0.38 |
  | 0.4 K | 5.64 | 0.12 |
  | 0.5 K | 5.71 | 0.03 |
  | 0.6 K | 5.80 | 0.009 |

  Repeated again after the unit model changes D29-D32, with k_o as refitted in D30 (margin kept
  at 0.6 K). This sweep was run before the MPC corrections D35 and D36, after which the nominal
  MPC uses 6.27 kWh at 0.6 K, still with 0.000 K·h outside the band; it was not repeated:

  | Margin | Energy (kWh) | K·h outside the band |
  |---|---|---|
  | 0.3 K | 6.07 | 0.18 |
  | 0.4 K | 6.15 | 0.030 |
  | 0.5 K | 6.23 | 0.001 |
  | 0.6 K | 6.31 | 0.000 |

- **Reason for 0.6 K:** the comparison with the PI is made at equal comfort.

---

## Ventilation

Everything in D23-D26 was written down before any ventilation run, as a guard against tuning
to the answer. Changes made after seeing results are added as separate entries.

### D23. Ventilation fan power
- **Timing:** a priori.
- **Problem:** until now ventilation was free, which would favour any strategy that ventilates
  more.
- **Model:** fan affinity law around a design point. P = P_d (V / V_d)³, with V_d = 0.03 m³/s
  (108 m³/h) and P_d = 30 W, i.e. a specific fan power of 1 000 W/(m³/s) at design flow. This
  is an assumed value for a small ventilation unit, not taken from a product.
- **Computation:** ventilation flow is held between controller samples, so fan energy is the
  exact sum of P(V_k) dt.
- **Not modelled:** heat recovery (in summer it would mostly be bypassed).

### D24. Requirement and minimum flows
- **Timing:** a priori.
- **Requirement:** CO2 at or below 1 000 ppm during occupied hours. Reported as occupied hours
  above 1 000 ppm and the occupied maximum.
- **Minimum flow when occupied:** the area component of the ASHRAE Standard 62.1 rate for
  office space, 0.3 L/(s·m²) × 20 m² = 0.006 m³/s, to remove pollutants that do not scale
  with occupancy.
- **Unoccupied:** no flow, except when a strategy ventilates for cooling (D26).
- **Context:** the existing fixed schedule (0.0105 m³/s while occupied) is close to the 62.1
  breathing-zone rate for two people, people plus area component (2 × 2.5 L/s + 6 L/s =
  0.011 m³/s), before any correction for air distribution. Edition and source of these rates:
  [references.md](references.md).

### D25. Strategies compared
- **Timing:** a priori. All strategies use the same temperature and humidity control
  (PI + sigma-delta with the H2 cascade on the fan), so differences come from ventilation alone.
- **V0, fixed schedule:** 0.0105 m³/s from 8 to 18 h.
- **V1, demand-controlled, continuous:**
  - PI on measured CO2, set point 900 ppm (100 ppm of headroom below the requirement).
  - Output between the occupied minimum and 0.03 m³/s, with back-calculation anti-windup.
- **V2, demand-controlled, staged:** three fan levels (0.006, 0.012 and 0.024 m³/s), switched
  up at 850 / 950 ppm and down at 750 / 850 ppm, with a 5-minute minimum dwell time on each
  level. This is a generic form of staged demand-controlled ventilation.
- **V3, V1 plus night flush:**
  - Mechanical ventilation at 0.03 m³/s when unoccupied, outdoor air at least 2 K cooler than
    the room, and the room above 22 °C.
  - It tests whether the free cooling saves more compressor energy than the fan uses.

### D26. Tuning of the CO2 PI (V1)
- **Timing:** a priori.
- **Process:** dC/dt = (V_inf + V)(C_o − C)/V_room + N G/V_room. Linearised at 900 ppm with
  two people, the gain is ∂C/∂V = −(C − C_o)/(V_inf + V) and the time constant is
  V_room/(V_inf + V).
- **Tuning:** SIMC with τ_c = 900 s, as for the temperature loop, computed from the parameters
  (`dcv_tuning`). (For the temperature loop the achieved closed-loop time differs from the
  design because the unit is nonlinear; see D39. The CO2 loop has no such nonlinearity.)
- **Arithmetic correction:** the first version of this entry used V_inf + V ≈ 0.0175 m³/s.
  The steady state at 900 ppm requires N G / (C − C_o) = 0.0217 m³/s, so τ ≈ 2 500 s and
  ∂C/∂V ≈ −22 000 ppm per m³/s. The correction was committed together with the first
  ventilation runs. It changes only the numbers written here: the code has always computed
  the gains from the parameters (`dcv_tuning`), so every run used the corrected values.

### D27. Variable-occupancy scenario
- **Timing:** after observation. Written before running the new scenario; the strategies and
  their parameters stay as fixed in D25.
- **Observation:** with constant occupancy (two people, 8-18 h), demand-controlled ventilation
  used more air and energy than the fixed schedule. There is no variation for it to exploit,
  and its 900 ppm set point is stricter than where the fixed schedule ends up (≈ 1 000 ppm).
  This is a property of the scenario, not of the strategies.
- **Addition:** an office day with variable occupancy:

  | Hours | People |
  |---|---|
  | 8-9 | 1 |
  | 9-12 | 2 |
  | 12-13 | 0 (lunch) |
  | 13-15 | 2 |
  | 15-16 | 3 (meeting) |
  | 16-18 | 1 |

  Mean occupancy is 1.6 people, below the design value of 2, with a peak of 3.
- **Expectation stated beforehand:** the fixed schedule, sized for two people, exceeds 1 000 ppm
  during the meeting. Demand control meets the requirement with less air on average.
- **Exploratory, labelled as such:** a sensitivity of the night flush to its flow
  (0.01-0.03 m³/s), because fan power scales with the cube of the flow. V3 itself stays as
  defined.

### D28. Ventilation results against the expectations
- **Timing:** after the runs. This entry only records outcomes; no parameter was changed.
- **D27 expectation:** the fixed schedule exceeds 1 000 ppm during the meeting, and demand
  control meets the requirement with less air on average. Partly met:
  - the fixed schedule exceeded the limit for 3.9 h (peak 1 189 ppm), and demand control met
    the requirement;
  - "less air" holds only against a fixed flow sized for the peak (1 176 against 2 095 m³).
    Against the 2-person fixed flow, demand control used 4 % more air.
- **Unexpected:** night flush at the design flow costs more fan energy than it saves. It pays at
  roughly half that flow (exploratory sweep). V3 was not redefined; the sweep is reported
  separately.

## Corrections after the end-of-project check

At the end of the project the code and documents were checked again, looking for errors. The
check found two defects (D20, D21) and four inaccuracies (D22). The fixes follow.

### D20. Minimum condensing lift in the cycle model
- **Timing:** after observation; found in the end-of-project check.
- **Defect:** COP = (η T_ev − a)/(T_o − T_ev + a) had no floor on T_cond − T_ev. At minimum speed
  the COP reached 17.9 at 20 °C outdoors, 77 at 15 °C, and went negative below about 13 °C.
  36 % of the PI's running time in the summer scenario was at COP > 8.
- **Fix:** T_cond ≥ T_ev + 15 K, which represents head-pressure control holding the
  compressor's minimum pressure ratio and the expansion device's minimum pressure difference.
  Over 0-45 °C outdoors, all speeds and fan settings, COP now stays within 2.5-9.3 and is
  continuous where the floor takes over.
- **Effect on results:**

  | Result | Before | After |
  |---|---|---|
  | PI energy | 6.21 kWh | 6.85 kWh |
  | PI saving vs on/off | 45 % | 40 % |
  | MPC saving vs set-back PI | 3.0 % | 2.7 % |
  | Time-of-use MPC vs set-back PI, cost | −10 % | **+0.5 %** |

  The tariff conclusion was reversed: night pre-cooling only looked profitable because
  night-time cooling was nearly free in the defective model.
- **Uncertainty:** the 15 K value is a parameter (`ACParams.dT_lift_min`). Absolute part-load
  energies depend on it.
- **Later comparison with published data:** a real 2.5 kW inverter split declares a net EER of
  16.5 at 20 °C outdoors and part load ([results.md](results.md#comparison-with-a-real-unit)).
  A COP near 18 at 20 °C was therefore not implausible in itself; the defect was the unbounded
  and, below 13 °C, negative COP. With the 15 K floor the model reaches a net EER of about 5 at
  that point, so its night-time and part-load efficiency is conservative, and every conclusion
  that depends on the night-time COP (the tariff result above, D34) depends on this parameter.

### D21. Sigma-delta ripple design used the wrong capacity
- **Timing:** after observation; found in the end-of-project check.
- **Defect:** T_ref was computed with the sensible capacity at full speed (wet coil, 1 685 W).
  While cycling, the unit runs at minimum speed with a dry coil and delivers 1.5 times more
  sensible cooling per unit of demand. The measured ripple was therefore 0.74 K instead of the
  0.5 K intended.
- **Fix:** use the sensible output per unit of demand at minimum speed; T_ref = 210 s (163 s
  after D30, which raised the capacity at minimum speed).
- **Checks:**
  - Measured median ripple over about 100 low-load cycles: 0.50 K. A test now checks this.
  - The MPC margin sweep was repeated: 0.6 K still matches the PI's comfort (0.009 K·h).

### D22. Documentation and test corrections
- The coupled air-buffer moisture mode is about 400 s, not 900 s (≈ 820 s is the buffer side
  alone). The test bound was tightened from 300-1000 s to 350-450 s.
- The room dew point under PI control is about 17.7 °C while occupied, not about 16 °C. The
  coil at minimum speed (18.3-19.4 °C once the start-up ramp has decayed) is still dry, but
  with a margin of about 1 K, falling to 0.6 K at night. (Superseded by D29: the coil is now
  partially wet at these conditions.)
- The sigma-delta mean-demand test tolerance went from 0.01 to 0.003. Measured errors are
  ≤ 0.0017; at v = 0.05 the old tolerance allowed 20 %.
- The evaporator thermistor in fan-only mode now reads the wet-bulb temperature only if the coil
  is wet, and the return-air temperature if it is dry. (Refined in D42: the reading now moves
  continuously between the two as the film evaporates.)

## Changes after the external review

An external review of the results raised the points below. D29-D32 change the unit model; their
parameter values were fixed before any study was rerun, with one exception: the check at the end
of this round found that k_o had not been fitted like k_e, and it was refitted after the reruns
(D30). D33 records what the history shows about
a reference used in the MPC comparison. D34 lists the effect on the results.

### D29. Partially wet coil
- **Timing:** after observation (external review).
- **Observation:** the coil had a uniform surface temperature and was either completely dry or
  completely wet. Under PI control at minimum speed and full fan its surface (18.3-19.4 °C) was
  above the room dew point (about 17.7 °C), so it removed no water, while the evaporating
  temperature (16.9-18.0 °C) was below the dew point of the return air in about 70 % of those
  samples. A real coil is coldest where the air leaves
  it and would be partly wet. "No water removed" and the occupied RH of 64 % were therefore upper
  bounds within the model.
- **Decision:** refrigerant at a uniform evaporating temperature; air-side and refrigerant-side
  conductances in series, so the surface temperature varies along the air path. The coil is dry
  from the air inlet to the point where the surface reaches the dew point of the entering air,
  and wet after it. This is the structure of the part-wet, part-dry calculation of the
  detailed Coil:Cooling:Water model in the EnergyPlus Engineering Reference ([references.md](references.md)):
  wet/dry boundary where the surface reaches the air dew point, enthalpy-based conductance
  1/(c_s/UA_internal + c_p/UA_external) in the wet part, outlet state from an effective
  condensing surface. With a uniform refrigerant temperature both parts have closed forms; c_s
  is the secant slope of the saturation enthalpy over the wet part, which makes the two parts
  meet exactly at the boundary.
- **Checks:** the wetted fraction goes continuously from 0 to 1; capacity, outlet state and
  water removal have no jump at either boundary. At the rating point the coil is fully wet; SHR
  (0.73) and supply temperature (15.1 °C) are as before, and the mean surface is 13.4 °C against
  the 13.0 °C of the former uniform surface. Tests in `tests/test_ac_unit.py`.
- **Not modelled:** the superheated end of the evaporator (warmer, so drier), the temperature
  distribution across the fins (fin roots colder than the mean), and the effect of condensate on
  the air-side heat transfer. They act in opposite directions on the wetted area, so the new
  values are estimates, not bounds.
- **Numerical consequence:** where a trajectory crosses a regime boundary the right-hand side has
  a kink and RK4 loses order on that step (D8).

### D30. Capacity depends on evaporating and outdoor temperature
- **Timing:** after observation (external review).
- **Observation:** the capacity depended only on compressor speed and outdoor temperature, so a
  slower indoor fan did not reduce it. In a real unit less airflow lowers the evaporating
  temperature, the density of the suction gas and the refrigerant mass flow.
- **Decision:** Q = s Q_rated [1 + k_e (T_evap − T_evap,rated)] [1 − k_o (T_o − 35 °C)], solved
  together with the coil (one root, Brent's method). This is the role of the CAP-FT and CAP-FFF
  curves of the EnergyPlus DX coil model, which multiply the rated capacity by functions of the
  entering wet bulb and outdoor dry bulb, and of the flow fraction; here the flow and humidity
  dependence come out of the coil.
- **Values:** from the EnergyPlus example curve WindACCoolCapFT (test file WindACAuto.idf,
  EnergyPlus 9.6.0). At the rating point of this model (19.7 °C entering wet bulb, 35 °C) its
  slopes, relative to its value there, are +3.58 % per K of entering wet bulb and −1.08 % per K
  of outdoor temperature.
  - k_e = 0.040 per K: the value for which this model, with its coil, has the same wet-bulb
    slope at rating (fitted 0.0395, rounded). At 35 °C outdoors the outdoor factor is 1, so
    this fit does not depend on k_o.
  - k_o = 0.016 per K, fitted the same way once k_e is set: the value for which the unit's total
    capacity at fixed entering air (27 °C, 19.7 °C wet bulb, full speed and fan), coil included,
    falls with the outdoor temperature by the curve's 1.08 % per K (central difference over
    34.5-35.5 °C; fitted 0.0156, rounded; with 0.016 the slope is 1.10 % per K). k_o is larger
    than the curve's slope because less compressor capacity lets the evaporator warm, which
    gives back part of the loss. The rated evaporating temperature and the calibration of the
    cycle efficiency are at 35 °C and do not depend on k_o: the net EER at rating stays 3.50.
  - A test checks both slopes.
- **Correction after the end-of-round check (after observation):** the first version of this entry
  set k_o = 0.011 per K, the curve's outdoor slope applied to the compressor at a fixed
  evaporating temperature (it was 0.01 before D30), without fitting it through the coil. With the
  coil the unit's capacity at fixed entering air then fell by only 0.76 % per K. The check at the
  end of this round of changes found the mismatch after the studies had been rerun. k_o was
  refitted as above, because its fit was inconsistent with that of k_e, not because of a result,
  and every study was rerun. The three-day energies changed by −0.4 to +2.1 %. Effect on the
  results of D34, before the corrections D35-D37 (values with k_o = 0.011 in brackets): PI saving against on/off 30 % (29 %), PI
  occupied RH 62.4 % (62.6 %), humidity control +9 % against the PI (unchanged) and −23 to −24 %
  against on/off (−23 %), MPC against the set-back PI −1.6 % (−1.7 %), time-of-use MPC cost −6.6 %
  (unchanged) for 16 % more energy (18 %). The MPC with model parameters off by 20 % now uses
  1.6 % more than the nominal MPC (0.6 %).
- **Check against the flow curve of the same file,** WindACCoolCapFFF = 0.8 + 0.2 ff (stated
  range 0.5-1.5), at rating conditions:

  | Airflow | Model | Curve |
  |---|---|---|
  | 70 % | 0.933 | 0.94 |
  | 55 % | 0.876 | 0.91 |
  | 40 % | 0.788 | 0.88 (extrapolated) |

  At minimum compressor speed (room at 25 °C and 60 % RH, 32 °C outdoors) the losses are 4 % at
  70 % and 10 % at 40 % airflow.
- **Consequences:** at low load the evaporator is warmer, so minimum speed delivers more: 1.01 kW
  (0.99 kW net of the indoor fan heat) at 27 °C / 19 °C wet bulb and 35 °C; before the change it
  was 0.3 × 2.5 kW = 0.75 kW. The unit cycles more at low load. Frost protection is no longer
  reached at rating conditions (D7).
- **Limits:** one example curve, of a fixed-speed window unit, not data of an inverter split;
  linear dependence on both temperatures; the two slopes match the curve at the rating point
  only.

### D31. Indoor fan power at low speed
- **Timing:** after observation (external review). The value was fixed from published data before
  any rerun.
- **Observation:** with P = 25 W · f³ the fan drew 1.6 W at 40 % speed (32 W per m³/s). That
  leaves out the motor and drive losses and the control electronics, which do not fall with the
  cube of the speed.
- **Decision:** P = 5 W + 20 W · f³: 25 W at full speed as before, 6.3 W at 40 %.
- **Evidence:** the published data of a 2.5 kW wall-mounted inverter split (Daikin FTXM25R +
  RXM25R, [references.md](references.md)) give an indoor-unit power input of 25 W in nominal
  cooling with 630 m³/h at the highest fan step, that is 0.14 W per L/s (the model: 25 W at
  450 m³/h, 0.20 W per L/s), and 6 W in thermostat-off mode, with the compressor stopped. The
  5 W makes the model's fan-only power at minimum speed (6.3 W) close to those 6 W. That the
  thermostat-off power is mostly the indoor fan running slowly is an assumption; the data book
  does not say what it contains.
- **Not adopted:** the review suggested 0.3-0.6 W per L/s at nominal flow. The published data of
  the wall unit are lower (0.14 W per L/s), so the nominal 25 W was kept.

### D32. Net EER
- **Timing:** after observation (external review).
- **Observation:** the EER was the heat taken by the refrigerant over the electrical power. The
  indoor fan heat, which stays in the room, was counted as cooling.
- **Decision:** the rated value 3.5 is now a net EER: (coil capacity − indoor fan power) over the
  electrical power of compressor and both fans; the cycle efficiency is calibrated to it. Every
  EER in the documentation is net. The function `net_cooling` in `ac_unit.py` gives the net
  cooling of an operating point.
- **Effect:** at rating the gross air-side EER is 3.54 against the net 3.50.

### D33. The 6:00 start of the set-back PI used as the MPC reference
- **Timing:** not recorded. The set-back PI with a fixed start at 6:00 entered the code in the
  same commit as the MPC and its first results (a3fee3e, "Add economic MPC with Kalman filter,
  supervising the PI"). No other start time appears in the history, and no reason for 6:00 was
  written down. Whether 6:00 was chosen before or after the MPC results were seen cannot be told
  from the history, so it is treated as an unoptimised reference, not as one fixed a priori.
- **Sensitivity (after observation, exploratory, model after D29-D32 with k_o refitted):**
  set-back PI with other fixed start times, same three days:

  | Start | Energy (kWh) | Cost, time-of-use (EUR) | K·h outside band |
  |---|---|---|---|
  | 5:00 | 6.50 | 2.14 | 0.000 |
  | 6:00 | 6.41 | 2.13 | 0.000 |
  | 7:00 | 6.29 | 2.12 | 0.000 |
  | 7:30 | 6.23 | 2.12 | 0.000 |
  | MPC (for comparison) | 6.27 | 2.15 | 0.000 |

  The MPC row is after the corrections D35 and D36; before them it was 6.31 kWh and 2.19 EUR.
  The set-back PI rows were rerun after the final audit and are unchanged.
- **Consequence:** the MPC's energy saving against the set-back PI (2.2 % with the 6:00 start;
  1.6 % before D35-D36) is smaller than the effect of the PI's start time. With a start at 7:00
  the MPC saves 0.3 % (6.27 against 6.29 kWh; before D35-D36 the two matched, 6.31 against
  6.29 kWh), and a start at 7:30 uses 0.6 % less than the MPC. The comparison shows that the
  MPC finds a good start by itself, not that it saves energy against a well-tuned schedule. The script keeps the 6:00
  reference so that the table stays comparable with earlier versions.

### D34. Effect of D29-D32 on the results
- **Timing:** after the reruns. One parameter was changed after them: k_o, refitted because the
  end-of-round check found its fit inconsistent with that of k_e, not because of a result
  (D30). The table is after that change; "before" is the model before D29-D32. "After" is the
  state at the end of that round, before the corrections D35-D37 of the final audit, which
  changed some of these values again (current values in D35-D37 and [results.md](results.md)).

  | Result | Before | After |
  |---|---|---|
  | On/off thermostat energy | 11.38 kWh | 10.54 kWh |
  | PI (C) energy | 6.85 kWh | 7.39 kWh |
  | PI (C) saving against on/off | 40 % | 30 % |
  | PI (C) mean occupied RH | 63.7 % | 62.4 % |
  | PI (C) water drained in three days | 0.0 kg | 0.9 kg |
  | PI (C) compressor starts | 169 | 251 |
  | Humidity control (H1, H2) energy against the plain PI | +10 % | +9 % |
  | Humidity control energy against on/off | −34 % | −23 to −24 % |
  | MPC energy against set-back PI (6:00) | −2.7 % | −1.6 % |
  | Time-of-use MPC cost against set-back PI (6:00) | +0.5 % | −6.6 % |
  | Time-of-use MPC energy against set-back PI (6:00) | +7 % | +16 % |

- **Reversal:** with the revised unit the time-of-use tariff pays: the MPC that knows the tariff
  costs 6.6 % less than the set-back PI, for 16 % more energy. Before, it was level. The
  attribution to the individual changes is in [results.md](results.md#model-predictive-control).
- **Also changed:** the PI with handler C now starts the compressor more often in total than
  on/off (251 against 172; before, 169 against 174; 171 for on/off after D37), and the
  time-of-use MPC uses more energy than the PI at a constant set point (7.44 against 7.39 kWh;
  before, 6.39 against 6.85; 7.41 after D35-D36).
- **Unchanged:** the energy ordering of the five standard controllers, the humidity controllers
  reaching 55 % RH, the robustness of the MPC, and the ventilation findings.

## Corrections after the final audit

A final audit read the code, the scripts and the documents against each other and reran the
studies. It found two defects in the MPC (D35, D36), two in the indicators (D37, D38), a tuning
statement that did not hold (D39), a defect in the frost-limited unit (D40), missing input and
consistency checks (D41) and modelling assumptions that were not written down (D42). D35-D38
change published numbers; every table in `figures/` and in [results.md](results.md) was
regenerated. No parameter was changed.

### D35. MPC set point clipped with the band at the end of the first prediction step
- **Timing:** after observation; found in the final audit.
- **Defect:** the PI's set point is T[0], the planned air temperature at t + 15 min, the end of
  the step the PI then tracks. It was clipped to the planning band in force at the solve instant
  t instead of the band at t + 15 min, which the LP itself used for T[0]. At the 18:00 solve the
  instant 18:00 counts as occupied (D17), so the PI was held at ≤ 25.4 °C from 18:00 to 18:15,
  a quarter hour that is unoccupied and that the planner lets float. At the 7:45 solve the
  occupied limit was not applied to the step ending at 8:00; that had no effect, because the LP
  already kept T[0] inside the occupied band.
- **How found:** by checking the time each planned quantity refers to (forecasts at mid-step,
  ventilation at step start, band at step end) against the use of the plan in `step`, and
  rerunning the MPC cases with the clip at t + 15 min.
- **Fix:** clip with `planning_band(t + Ts_pred)`. A test checks both edges of the occupied
  period.
- **Effect on results (this correction alone):** MPC 6.31 → 6.26 kWh, 2.19 → 2.15 EUR; tariff MPC
  7.44 → 7.38 kWh, 1.99 → 1.96 EUR; comfort unchanged (0.000 K·h). The combined effect with D36
  is in D36.

### D36. MPC capacity and cost net of the indoor-fan heat
- **Timing:** after observation; found in the final audit.
- **Defect:** the LP variable is the cooling Q_c of the RC model, and the Kalman filter is fed
  the net cooling Q_sens − P_fan_in (the indoor-fan heat stays in the room). The capacity limit
  U_max and the energy per unit of cooling c were computed from the gross Q_sens: U_max was
  25 W (1.5 %) too high and c about 2.8 % too low at half speed.
- **How found:** by checking that the planner, the filter and the RC model use the same
  quantity for the cooling.
- **Fix:** subtract P_fan_in in both. A test checks that the LP capacity equals the filter input
  at the same operating point.
- **Effect on results (this correction alone):** MPC 6.31 → 6.30 kWh; tariff MPC cost 1.99 →
  1.98 EUR. Together with D35 (current values):

  | Result | Before D35-D36 | After |
  |---|---|---|
  | MPC energy | 6.31 kWh | 6.27 kWh |
  | MPC cost, time-of-use | 2.19 EUR | 2.15 EUR |
  | MPC starts | 149 | 148 |
  | Time-of-use MPC energy | 7.44 kWh | 7.41 kWh |
  | Time-of-use MPC cost | 1.99 EUR | 1.97 EUR |
  | Time-of-use MPC starts | 186 | 183 |
  | MPC, model off by 20 %: energy, starts | 6.41 kWh, 179 | 6.33 kWh, 168 |
  | MPC, biased forecast: energy, cost, starts | 6.28 kWh, 2.18 EUR, 147 | 6.25 kWh, 2.15 EUR, 145 |
  | MPC energy against set-back PI (6:00) | −1.6 % | −2.2 % |
  | MPC energy against set-back PI (7:00) | +0.3 % | −0.3 % |
  | Time-of-use MPC cost against set-back PI (6:00) | −6.6 % | −7.5 % |
  | Time-of-use MPC cost against MPC without tariff | −9 % | −8 % |
  | Time-of-use MPC energy against set-back PI (6:00) | +16 % | +16 % |
  | Effect of 20 % model error on MPC energy | +1.6 % | +1.0 % |
  | Effect of the biased forecast on MPC energy | −0.5 % | −0.3 % |

  Comfort is 0.000 K·h outside the band in every case, before and after. The conclusions stand:
  the MPC saves little against a well-scheduled PI, and the tariff pays at a price in energy.

### D37. Compressor starts counted at window boundaries
- **Timing:** after observation; found in the final audit.
- **Defect:** the indicators counted the first row of a result as a start whenever the
  compressor ran in it. The studies evaluate a window after three warm-up days; when the window
  begins inside a running cycle, that row is not a start. Of the published runs only the on/off
  thermostat was running at the window start.
- **How found:** by recounting the starts of every published run from the commands, including
  the row before each window.
- **Fix:** a windowed result records whether the compressor ran on the step before it
  (`UnitLog.running_before`), and the first row counts as a start only if it did not. A test
  checks a window that starts inside a cycle.
- **Effect on results:** on/off thermostat 172 → 171 starts, 10.4 → 10.3 starts per running hour
  (also the on/off reference of the humidity table). No other run changes.

### D38. Fan travel only between running samples
- **Timing:** after observation; found in the final audit.
- **Defect:** the fan travel of the humidity study (total variation of the fan command per
  running hour) concatenated the running samples, so the step between the last command of one
  compressor cycle and the first of the next was counted although the fan was off in between,
  and the stale command of the last row was included. Exactly half of H1's value was this
  cycle-to-cycle step: its proportional term moves while the compressor is off, whereas H2's
  integrator is held.
- **How found:** by recomputing the column from the commands with only consecutive running
  samples.
- **Fix:** `metrics.fan_travel_per_h` counts |f[k+1] − f[k]| only where both rows run, drops the
  last row like the other indicators, and returns NaN without running time. RH standard
  deviation and maximum in the same table now also drop the last row (no change in value).
- **Effect on results:** fan travel per running hour for H1, nominal / dirty coil / 3.5 kW unit:
  0.27 / 0.30 / 0.53 → 0.14 / 0.15 / 0.26; H2 unchanged (0.26 / 0.29 / 0.36). The ordering is
  reversed: H1 now moves the fan less than H2 for all three units, whereas the table suggested
  the opposite. The column had not been quoted in the documentation.

### D39. Process gain of the temperature PI
- **Timing:** after observation; found in the final audit. Analysis only; the tuning was kept.
- **Observation:** `room_pi_tuning` applies the SIMC rules with Q_d = 2620 W per unit of demand,
  the total compressor capacity (sensible plus latent) at full speed and the rated evaporating
  temperature, at 32 °C outdoors. The air temperature responds only to the net sensible cooling,
  and the unit is strongly nonlinear in the demand. Gains measured at the operating point:
  - below s_min (sigma-delta cycling, where the PI spends most of its time in the standard
    runs): about 3250 W of net sensible cooling per unit of demand;
  - above s_min, where the coil is wet: an incremental sensible gain of 700-1000 W.

  So the documented design (τ_c = 900 s) does not describe the achieved loop. Measured
  closed-loop 63 % times for set-point steps of 0.5 K: 600-1140 s at low load (design 960 s,
  τ_c + θ plus the zero-order hold) and 2100-2340 s in the modulating range.
- **Alternative tested:** taking Q_d as the net sensible output at s_min per unit of demand
  (≈ 3250 W). It hardly changed the low-load response (600-1200 s), because the sigma-delta
  modulation dominates it, and made the modulating range slower still. With it the PI (C) used
  7.36 instead of 7.39 kWh with 243 instead of 251 starts, but the MPC runs left the band
  (0.003-0.006 K·h outside it for the nominal, tariff and biased-forecast cases, against
  0.000), so the 0.6 K margin of D18 would have had to be retuned.
- **Decision:** keep the tuning and the D18 margin; document the gains and the measured
  closed-loop times in the docstring of `room_pi_tuning` and in the known limitations.
- **Effect on results:** none.

### D40. Frost-limited cycle efficiency at the effective speed
- **Timing:** after observation; found in the final audit.
- **Defect:** when frost protection holds the evaporating temperature at 0.5 °C, the compressor
  in effect slows to the speed s_eff = s Q_coil(0.5 °C) / Q_comp(s, 0.5 °C) whose capacity the
  coil takes. The cycle efficiency was still evaluated at the commanded speed s. At 20 °C and
  10 % RH return air, −10 °C outdoors, s = 1 and 40 % fan (s_eff = 0.37) this overstated the COP
  by 25 %; at 22 °C, 0.008 kg/kg and 30 °C by 0.3 %. Where s_eff fell below s_min the unit
  delivered that capacity continuously, which a compressor cannot do. The module docstring also
  said frost protection was active only at low fan speed and high load; it is also reached at
  full fan with cool, very dry return air and cold outdoor air.
- **How found:** by sweeping the unit model over the envelope (return air 15-35 °C and 10-95 %
  RH, outdoor −10 to 50 °C, all speeds and fan settings) and comparing the frost-limited points
  with the cycle evaluated at s_eff.
- **Fix:** the efficiency is evaluated at s_eff. If s_eff < s_min the compressor cycles at s_min
  with duty s_eff / s_min: the returned state is the cycle average and each on-period runs at
  s_min with the load Q_coil s_min / s_eff. Power and COP stay continuous at the onset of frost
  protection and at s_eff = s_min (tests).
- **Effect on results:** none. No published run is frost-limited; the evaporating temperature
  stays above 12.5 °C in the humidity runs. At the example point the COP goes from 9.05 to 7.23.

### D41. Robustness checks
- **Timing:** after observation; found in the final audit, by running the code with inputs
  outside those of the studies.
- **Defects and fixes** (behaviour before → after):
  - **Divergence:** RK4 above its stability limit (≈ 1 100 s for the default room) returned
    nonsense without warning, for example W between −12 and 19 kg/kg at Δt = 1 200 s. The
    simulator now stops as soon as the state is non-finite or W, C or an active buffer's φ_b
    is negative, and reports the stability limit estimated from the Jacobian at the initial
    state (`rk4_stability_limit` added next to the Euler limit).
  - **Input validation:** negative room parameters (h_i, ach, buffer), a negative ventilation
    flow (which reversed the ventilation terms and made the CO2 balance unstable: 6·10¹¹ ppm
    after a day), RH_o outside [0, 1], sunset before sunrise, a result window outside the run,
    and windows without occupied rows now raise a clear error, or give NaN with a warning for
    the occupied means, instead of silent nonsense or an obscure numpy error.
  - **Zero-capacity buffer:** with zero buffer density, depth or slope the buffer stored no
    water but still exchanged it, an unlimited source (water residual 5.7 kg in a day). A buffer
    that stores no water now exchanges none; the conductance also tends to 0 for h_i = 0 or zero
    permeability instead of raising ZeroDivisionError.
  - **Periodic state with frozen states:** with the buffer off, φ_b has no dynamics and the
    Newton matrix is singular; the shooting worked only through finite-difference round-off and
    failed for φ_b = 0. States without dynamics now keep their guess and Newton solves for the
    others; if such a state drifts (for example an open mass integrator), a clear error says that
    no periodic state exists.
  - **Daily windows across midnight:** occupancy intervals, ventilation and comfort schedules and
    tariff peaks with start > end (for example 22 h to 6 h) were silently never active. They now
    cross midnight; start = end and hours outside [0, 24] are rejected.
  - **Sampling period:** every discrete element (PI, low-load handler, Kalman filter, MPC, humidity
    and ventilation loops) uses its own Ts, so running a 60 s controller every 30 s silently
    doubled the PI integral gain and the Kalman prediction step. The simulator now rejects a
    controller whose Ts differs from `dt_ctrl`, and composite controllers reject parts with
    different periods.
  - **Read-only discretisation:** the cached exact discretisation of the control model returned
    the same arrays to every caller, so an in-place change would have corrupted the cache. They
    are now read-only.
  - **Failed LP:** after a failed solve the set point and the move were taken from step 0 of the
    previous plan, a step that had already passed. They are now taken from that plan shifted to
    the current step.
  - **Smaller items:** dew point of dry air −∞ instead of NaN; moisture time constant ∞ for a
    sealed room instead of ZeroDivisionError; peak shading of the MPC figure for a window that
    starts inside the peak; step plots no longer draw the command of the last row, which acts
    after the end of the result.
- **Effect on results:** none. The published runs use the default buffer, Δt = 60 s, controllers
  sampled at their own 60 s, non-negative flows, daily windows within one day and no failed LP
  solves; every table was regenerated and only the values listed in D35-D38 changed.

### D42. Documented modelling assumptions of the coil
- **Timing:** after observation; found in the final audit.
- **Held water on a dry coil while the compressor runs:** the held water changes only through
  condensation while the compressor runs. On a coil that is dry under the current conditions it
  neither drains nor re-evaporates, even with the surface above the dew point; it re-evaporates
  in the next fan-only period. This is the on-cycle assumption of Henderson and Rengarajan
  (1996). It was not written down. In the published PI runs a dry coil with the compressor on
  occurs 83 (handler C) and 34 (handler A) times per day, each time during a 2-3 minute start-up
  ramp; in steady dry running with water held the error in the drain and latent bookkeeping is
  at most the holdup of 0.25 kg. Kept as a documented assumption (module docstring, known
  limitations).
- **Fan-only coil state:** below 5 % of the holdup the film evaporates in proportion to the
  water held, so M_c decays exponentially and never reaches zero. The surface temperature and
  wetted fraction stayed at the wet values (wet bulb, 1) as long as any water was held, so the
  thermistor read the wet-bulb temperature indefinitely. They now follow the same evaporation
  weight as the outlet humidity and go continuously to the dry values (return-air temperature,
  0).
- **Docstrings corrected:** `operate(M_c=None)` treats the coil as full with the compressor on
  and as dry in fan-only operation; with the compressor on the fan runs at least at f_min even
  if f ≤ 0 is commanded; during the start-up ramp the reported COP is that of the requested
  operating point, not delivered capacity over compressor power.
- **Effect on results:** none. The humidity controllers freeze their integrators while the
  compressor is off, and no controller in the studies reads the thermistor in fan-only
  operation.

## Seasonal study with measured weather

Everything in D43-D48 was written down and committed before any seasonal run or any evaluation
of the model's seasonal efficiency. Outcomes are recorded afterwards in a separate entry; no
parameter fixed here is changed after seeing them.

### D43. Weather data
- **Timing:** a priori.
- **Source:** typical meteorological years (TMY) from PVGIS 5.3 of the European Commission's
  Joint Research Centre, built by the ISO 15927-4 procedure from hourly data of 2005-2023.
  Irradiance from PVGIS-SARAH3 (satellite), the other variables from the ERA5 reanalysis.
  Horizon shading from the terrain model is included.
- **Sites:**
  - Munich (48.137 N, 11.575 E, 524 m): moderate summer, May-September mean 15.9 °C, 208 h
    above 25 °C, maximum 33.1 °C.
  - Madrid (40.417 N, 3.704 W, 665 m): hot, dry summer, May-September mean 23.4 °C, 1 450 h
    above 25 °C, maximum 40.1 °C.
- **Licence:** © European Union, reused under the Commission's reuse policy (CC BY 4.0). The
  files in `data/weather/` keep the hourly rows unchanged but only the columns used (time,
  T2m, RH, G(h), Gb(n), Gd(h), SP), converted from JSON to CSV with the location and the
  months selected in a header; `data/weather/README.md` says so.
- **Time stamps:** the JSON labels the time column UTC. Temperature, humidity and pressure are
  taken at the stamped hour; irradiance at the stamped hour plus the `irradiance_time_offset`
  given in the file (0.180 h in Munich, 0.171 h in Madrid), as the PVGIS documentation
  instructs for solar geometry. Checked while designing: rebuilding the horizontal global
  irradiance from the beam and diffuse values with the solar geometry of D44 leaves an RMS
  residual of 1.9 W/m² (Munich) and 1.8 W/m² (Madrid) over hours above 50 W/m²; without the
  offset 8.3 and 12.4 W/m². This fixes the time convention; it becomes a test.
- **Local time:** the room's schedules run on local clock time. The cooling season lies inside
  daylight saving time, so local time is UTC + 2 h in both cities throughout.
- **Between the hourly values:** linear interpolation, applied to the outdoor temperature, the
  vapour pressure and the irradiance on the window (D44), each at its own time stamps. The year
  is periodic (31 December runs into 1 January); the months of a TMY come from different
  years, and PVGIS joins them.
- **Humidity:** the outdoor vapour pressure is taken from the data (RH times the saturation
  pressure at T2m) and converted to a humidity ratio at the model's standard pressure. This
  keeps the outdoor dew point, which is what decides condensation on the coil. The station is
  about 6 % below standard pressure; converting at the station pressure would raise the
  humidity ratio by about 6 % at the same dew point, inconsistently with the room model at
  standard pressure.

### D44. Sun on the window
- **Timing:** a priori.
- **Window:** vertical, facing due south, unshaded; area and g-value as in `RoomParams`
  (3 m², 0.6), without the fall of transmittance at large angles of incidence. The synthetic
  day used 500 W/m² of peak irradiance on the window; a vertical south window in Munich
  receives about that at noon in July.
- **Solar position:** declination and equation of time from Spencer's Fourier series in
  fractional-year form, hour angle from UTC, the longitude and the equation of time.
- **Transposition:** isotropic sky:
  G_win = Gb(n) max(cos θ, 0) + Gd(h) (1 + cos β)/2 + ρ G(h) (1 − cos β)/2,
  with θ the angle of incidence, β = 90° and ground reflectance ρ = 0.2. The isotropic model
  underestimates the diffuse light from the circumsolar region, so the sun on the window is
  on the low side.
- **Check:** the solar position is compared with a second, independent formulation (the
  low-precision solar coordinates of the astronomical almanac). The two agree to 0.3-0.5° in
  declination and 0.9 min in the equation of time; the difference is the phase of the leap
  year cycle, which a typical year does not define.

### D45. Season, schedules and controllers
- **Timing:** a priori.
- **Season:** 1 May to 30 September (153 days), after 7 warm-up days (24-30 April) that are
  discarded. The run starts from the equilibrium state of `equilibrium_state(22, 22, W_o)`.
- **Occupancy and ventilation:** as in the controller studies, two people from 8 to 18 h and
  the fixed ventilation schedule, on every day of the season. Weekends are not modelled; this
  overstates the season's demand of an office and does not affect the comparison between
  controllers.
- **Set point:** the set-back schedule of the MPC study (D33): 25 °C from 6 to 18 h, 28 °C
  otherwise, for every controller. The on/off thermostat is given the same schedule (it took
  a constant set point until now).
- **Controllers:** on/off thermostat; PI with sigma-delta low-load handling (C); the same PI
  with the H2 cascade on the indoor fan (RH set point 55 %). The MPC is not run: its forecast
  model would have to be extended to measured weather, which is a separate piece of work.
- **No heating:** the unit only cools. Hours with the room below the comfort band while
  occupied are reported, not corrected.
- **Indicators per city and controller:** electrical energy of the unit; net cooling delivered
  (sensible plus latent, minus the indoor-fan heat, integrated with the state like the
  electrical energy); seasonal EER in use (the ratio of the two); compressor running hours and
  starts; the share of running time at minimum speed; degree-hours above and below the band
  while occupied; mean occupied RH and occupied hours above 60 % RH.
- **Runtime:** about 1 s of computing time per simulated day per run, measured; six runs of
  160 days are run in parallel processes.

### D46. Seasonal energy efficiency ratio (SEER)
- **Timing:** a priori.
- **Method:** Commission Regulation (EU) No 206/2012, Annex II, with the calculation in the
  Commission communication 2012/C 172/01, section 3.1, for variable-capacity units:
  - bin hours of the reference cooling season, 17-40 °C, 2 602 h in total;
  - part load Pc(Tj) = Pdesignc (Tj − 16)/(35 − 16);
  - anchor points at 35, 30, 25 and 20 °C. Where the unit can deliver the part load within
    ± 10 %, EERbin is the EER at that load. Where the part load is below the minimum capacity,
    EERbin = EERd (1 − Cdc (1 − Pc/Pdc)) at the minimum capacity, with Cdc = 0.25 (the
    default value);
  - linear interpolation between the anchor points, constant beyond them;
  - SEERon = Σ hj Pc(Tj) / Σ hj Pc(Tj)/EERbin(Tj);
  - SEER = QC / QCE, with QC = Pdesignc × 350 h and QCE = QC/SEERon + 221 h PTO + 2 142 h PSB
    + 0 h POFF + 2 672 h PCK (a unit that also heats).
- **Check on the published data, done while designing:** with the declared points of the
  reference unit (D30, [references.md](references.md)), Pdesignc = 2.50 kW and the low-power
  modes of its data book (PTO 6 W, PSB 1 W, POFF 1 W; no crankcase heater power listed, so
  PCK = 0), the method gives SEERon = 8.96 and SEER = 8.654, against the published 8.65, and an
  annual consumption of 101 kWh, as published. The implementation becomes a test with this
  case.
- **The model's SEER:** the same method applied to the model, with its points computed at the
  rating indoor state (27 °C, 19 °C wet bulb) and full indoor fan, as in the rated comparison:
  - Pdesignc = the model's net capacity at 35 °C;
  - at each anchor temperature, the compressor speed that delivers Pc(Tj), or the minimum
    speed if Pc(Tj) is below the minimum capacity;
  - the model has no standby electronics or crankcase heater, so the low-power modes of the
    reference unit are used (6, 1, 1, 0 W).
- **Local bins:** SEERon is also computed with the bin hours of each city's typical year
  (hours rounded to the nearest degree, 17-40 °C), for the published data and for the model.
  This separates the effect of the climate on the rating from the effects that only a dynamic
  simulation shows.

### D47. Expectations stated beforehand
- **Timing:** a priori.
- **E1, the model's SEER:** between 4.0 and 5.0 (estimate 4.5 from the model's part-load EERs
  in the rated comparison, 3.45-5.08), roughly half the published 8.65. The model's EER rises
  much less at part load than the real unit's (D30-D32, Limits); the SEER gap is the
  part-load gap, not new information.
- **E2, local bins:** with the published data, SEERon is higher with Munich's bins and lower
  with Madrid's than with the reference season. For the model the same order, with smaller
  differences, because its EER varies less with load.
- **E3, in use against the rating:** the seasonal EER in use is 10-35 % below the model's own
  SEERon with the local bins. Reasons: the room is held at 25 °C, not 27 °C, which lowers the
  evaporating temperature; the solar and internal gains do not follow the outdoor temperature,
  as the bin method assumes; and the oversized unit cycles at minimum speed for much of the
  season.
- **E4, controllers:** in both cities the PI uses less electricity than the on/off thermostat,
  by 15-35 %; the H2 cascade uses 3-15 % more than the PI and keeps the occupied RH lower.
- **E5, cities:** Madrid needs more than 1.5 times Munich's electricity with every controller.
  Low confidence on one more: the seasonal EER in use is higher in Madrid, because more of the
  running time is above minimum speed.

### D48. Supporting changes
- **Timing:** a priori.
- `OnOffThermostat` accepts a set-point schedule, as the PI already does; with a constant set
  point it behaves as before.
- The unit record gains the integrated net cooling (D45). It is an extra integral appended to
  the state; the existing states, and therefore every published number, must stay bit for bit
  the same, which `reproduce.py --check` verifies.
- New modules `solar.py` (solar position and transposition), `weather.py` (the PVGIS reader
  and the weather scenario) and `seasonal.py` (bin method); new script
  `scripts/compare_season.py`, added to `reproduce.py`.

### D49. Seasonal results against the expectations
- **Timing:** after the runs. This entry only records outcomes; nothing fixed in D43-D48 was
  changed. Tables: [results.md](results.md#cooling-season-with-measured-weather).
- **E1, the model's SEER: met.** 4.43 (SEERon 4.51), against the published 8.65. The anchor
  EERs are 3.45, 4.26, 5.07 and 4.34 (D degraded for cycling); the real unit's rise from 4.50 to
  16.51, so the gap is the known part-load gap of D30-D32.
- **E2, local bins: met.** With the published data SEERon is 10.39 with Munich's bins, 8.96
  with the reference season and 7.46 with Madrid's; for the model 4.52, 4.51 and 4.27, the
  same order with much smaller differences.
- **E3, in use against the rating: met.** The seasonal EER in use is below the model's SEERon
  with the local bins by 11-13 % with the PI (3.95 against 4.52 in Munich, 3.81 against 4.27
  in Madrid), by 10 % with the H2 cascade in Munich and by 26-29 % with the on/off thermostat.
  The bin method charges cycling with Cdc = 0.25 at the 20 °C point only; the thermostat,
  which runs at full speed and cycles at every load, falls furthest below it.
- **E4, controllers: met in Munich, partly in Madrid.**
  - The PI uses less electricity than the on/off thermostat: 30 % less in Munich, 22 % less in
    Madrid (in the range 15-35 %).
  - The H2 cascade uses 6.5 % more than the PI in Munich and lowers the occupied RH slightly
    (56.8 against 57.4 %). In Madrid it is identical to the PI: the occupied RH there is 37 %,
    far below the 55 % set point, so the cascade keeps the fan at full speed all season. The
    expectation did not foresee a climate where the humidity loop never acts.
- **E5, cities: the first part met, the low-confidence part not.** Madrid uses 5.8-6.4 times
  Munich's electricity. The EER in use is lower in Madrid (3.81 against 3.95 with the PI), not
  higher: the hotter outdoor air raises the condensing temperature more than the longer
  modulating periods help. In Madrid the PI still runs at minimum speed for 75 % of its running
  time (89 % in Munich).
- **Not expected:** Munich needs little cooling (33-46 kWh, most of it in August, whose TMY
  month comes from the hot summer of 2015). With no heating, the room is below the occupied
  band for long periods in May, June and September (about 3 600 K·h over the season); this is
  the room free-floating in a cool climate, not a fault of the controllers.

## Structure

### D19. Refactor of the simulation interface
- **Timing:** after the controller work, as a clean-up.
- **Changes:**
  - `simulate` takes a `duration` instead of an end time. Passing a duration where the end
    time was expected caused two test errors.
  - Everything about the unit (commands, operating point, coil water, electrical energy,
    drained condensate) moves into one `UnitLog`, instead of six optional fields on
    `SimResult`.
  - Actuation (held output, compressor starts, unit evaluation) moves into one class inside
    `simulate`.
  - `CoilState` has no fields with defaults added late.
- **Check:** all tests pass, and every result table and figure is bit-identical before and
  after.

## Test changes after failures

| Test | What failed | Wrong side | Resolution |
|---|---|---|---|
| CO2 decay | tolerance 10⁻⁸ ppm | test | tolerance = RK4 truncation error ≈ 10⁻⁹ of the amplitude |
| start-up power | expected the same power as one hour later | test | compare with full capacity at the same air state |
| coil water conservation | trapezoid across switches | test | integrated states (D11) |
| PI recovery after saturation | expected < 5 samples | test | analytical settling value (D12) |
| H1 recovery, H2 inner loop | too few samples; set point overwritten by the outer loop | test | sample count derived analytically; outer integrator fixed |
| Kalman, 48 h | 0.09 K error on T_m | test | slow but unbiased convergence (τ ≈ 14 h); checked over 6 days and through the decay rate |
| starts PI vs on/off | ratio < 0.6 after D14 | test | compared per running hour, the quantity relevant to wear |
| T_max occupied < 26.0 | 26.03 °C; then relaxed to 26.5 without a reason | test (threshold was arbitrary) | removed; comfort is checked once, against the degree-hour criterion in `tests/criteria.py` |
| `dt_ctrl` without controller | ValueError | **code** | ratio checked only when a controller is given |
| PI (C) starts fewer than on/off | 166 against 114 after D29-D32 | test | the ordering of total starts was a result read off a run (in the three-day study 169 against 174 before, 251 against 172 after, 251 against 171 after D37), not a requirement; now checked per running hour (in the test run 8.0 against 10.4), plus fewer total starts for handler A |
| PI rejects a load step, "continuous modulation" | compressor cycled at s_min after D30 | test | the premise (high load) no longer held with the larger minimum capacity; sun on the window raised from 300 to 600 W/m² so that it holds again |
| H2 inner loop on the real coil | coil 11.9 °C against a target of 11.5 °C | test | the comment said the target was 13 °C but the RH input made it 11.5 °C, below what the coil reaches at 40 % fan after D30; RH input corrected to give 13 °C, and the test now checks that the target is reachable |
| Dry/wet transition continuous | no switch left after D29 | test replaced | continuity checked over the whole partially wet range and at both boundaries |
| Unit inside the RK4 stages, order 4 | 2.5 after D29 | test | the run crossed a coil regime boundary (kink); split into a run in one regime (order 4) and a crossing run (order 2), see D8 |

---

## Test tolerance audit

All 177 numerical assertions then in the suite (before D29-D34) were reviewed and sorted
into the four kinds described in [`tests/README.md`](../tests/README.md). Seven had tolerances looser than they needed to be, or no
stated basis:

| Test | Before | After |
|---|---|---|
| control model vs plant | < 0.01 K | < 1e-4 K; measured 2e-5 K, the O(Δt²) error of midpoint input sampling |
| Euler stability limit | between 600 and 1000 s | equal to 2/\|λ\| of the moisture subsystem to 1e-6 (the Jacobian is block triangular) |
| Simpson residual at 60 s | < 1 J | equal to the value predicted from the measured order, within 10 % |
| start-up ramp | ±5 % and an ad hoc +0.05 | exact: the unit at the logged state with the analytical ramp, to 1e-12 |
| Kalman in closed loop | ±0.2 K | ±0.01 K; measured 0.002 K with an exact thermal model |
| closed-loop comfort, RH, robustness | thresholds per test (0.05, 0.1, 26.5 °C, 5 %) | one set of criteria with rationale in `tests/criteria.py` |
| MPC saving | > 5 % | strict ordering; the size is reported, not asserted |

## Known limitations

- Water held on the coil neither evaporates nor drains while the compressor runs with a dry coil
  (surface above the dew point); it re-evaporates only in the next fan-only period (D42).
  Physically the film would evaporate. In the published runs this happens only during the
  start-up ramps of 2-3 minutes; the water so retained is at most the coil holdup of 0.25 kg.
  The effect on reported numbers has not been computed.
- The PI's process gain is a nominal value; the unit is strongly nonlinear in the demand, and
  the closed-loop time in the modulating range is 2.2-2.4 times the design value (D39).
- The coil has a uniform refrigerant temperature: no superheated end and no fin temperature
  distribution (D29). Capacity dependence on evaporating and outdoor temperature is linear and
  fitted to one example curve of a fixed-speed unit at the rating point only (D30).
- The nominal MPC has a perfect weather forecast and exact model parameters; the robustness runs
  relax each of these once.
- The MPC's Kalman filter estimates the last sample's cooling with full fan speed, no start-up
  ramp and the forecast outdoor temperature. The disturbance state absorbs the difference.
- Single zone, well-mixed air; no radiant temperature or operative-temperature comfort.
- The controller studies use one synthetic summer day (sinusoidal outdoor temperature,
  constant outdoor humidity ratio, half-sine solar irradiance), repeated.
- The seasonal study (D43-D49) uses typical meteorological years of two cities, cooling only:
  every day is a working day, the window has no shading and no fall of transmittance at large
  angles of incidence, the isotropic sky underestimates circumsolar diffuse light, and the MPC
  is not run on it.
- The model's SEER uses the low-power modes of the real unit (6, 1, 1 W), since the model has
  no standby electronics.
- The model is verified (the equations are solved correctly) but not validated against
  measurements of a real room or unit.
- The MPC's model has no humidity, so it cannot trade sensible against latent load.
- MPC disturbance forecasts are held over each 15-minute step.
- The Kalman filter needs about a day to separate a wrong mass temperature from a heat
  disturbance with a single air sensor.
- Unit parameters (UA values, coil holdup, start-up constant) are plausible for a 2.5 kW split
  unit but not fitted to a specific product. Against the published data of a real 2.5 kW
  inverter split the model's EER is 23 % lower at rating, half at 25 °C outdoors and less than
  a third at 20 °C, where its 15 K minimum lift bounds the COP ([results.md](results.md#comparison-with-a-real-unit)).
- The ventilation fan follows the cube law of a fixed duct system. Under constant-pressure
  control the fan power falls more slowly with flow, which would reduce the benefit of the
  reduced night-flush flows.
