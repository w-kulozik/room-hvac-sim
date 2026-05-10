# Decision log

Every modelling and design decision that departs from the original specification or that was
revised during the work, with the evidence behind it.

**Timing** says when the decision was taken:
- *a priori*: before seeing the result it affects;
- *after observation*: changed because a result looked wrong. These carry the larger risk of
  tuning to the answer, so the evidence and the alternatives are spelled out.

Where a test was changed after it failed, the entry says whether the model or the test was wrong.

---

## Plant

### D1. Outdoor humidity: constant humidity ratio, not constant relative humidity
- **Timing:** after observation (first run of the scenario).
- **Observation:** holding RH_o = 55 % and converting at the instantaneous T_o makes W_o go from
  0.0080 to 0.0165 kg/kg within a day. Real air masses keep W roughly constant over a day while
  RH swings.
- **Decision:** W_o = W(26 °C, 55 %) = 0.0116 kg/kg, constant. The original behaviour remains as
  `SummerScenario(constant_W=False)`.
- **Effect:** the afternoon latent load is no longer doubled.

### D2. Moisture buffering of interior surfaces (EMPD)
- **Timing:** a priori. Raised when the moisture balance was written; the specification stored
  moisture in the air only.
- **Evidence:** painted gypsum (800 kg/m³, 3 mm penetration depth, du/dφ = 0.01, 60 m²) holds
  about as much exchangeable water as the room air: 72 kg against 65 kg of air, in humidity-ratio
  units. The apparent moisture time constant goes from 7 200 s to about 15 000 s.
- **Alternatives:** a constant multiplier on the air's moisture capacity. Rejected because it
  misses the exchange dynamics (≈ 900 s) and the release of moisture by warming walls.
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

### D7. Frost protection
- **Timing:** after observation.
- **Observation:** at 40 % fan and full speed the evaporating temperature reached −2.8 °C.
- **Decision:** cap the capacity so that T_evap ≥ 0.5 °C, as real units do.

### D8. Unit evaluated at every Runge-Kutta stage
- **Timing:** a priori, then confirmed by a test.
- **Evidence:** measured closed-loop convergence order 3.9. Evaluating once per sample and holding
  Q_sens gives 1.1, and at Δt = 60 s an error of 2·10⁻⁴ K instead of 1·10⁻¹⁰ K.

### D9. Start-up time constant 60 s (initially 120 s)
- **Timing:** after observation.
- **Observation:** with 120 s, about 30 % of the on/off thermostat's energy went into start-up
  ramps (average on time 6.5 min).
- **Evidence for the new value:** reported capacity time constants of small units are 40-90 s
  (O'Neal and Katipamula; Henderson).
- **Risk:** the value was revisited because a result looked large. It stays a parameter
  (`ACParams.tau_start`). On/off energy is the result most sensitive to it.

### D10. Water held on the coil (Henderson)
- **Timing:** a priori.
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
- **Decision:** T_t = T_i/2, the usual compromise for a PI (Åström and Hägglund). The test now
  checks the analytical settling value for three values of T_t.

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
  with ΔT_pp = 0.5 K (T_ref ≈ 320 s).
- **Effect on all PI-based results:** maximum occupied temperature 26.0 → 25.7 °C; starts per
  three days 84 → 123. Tables were regenerated.

### D15. Humidity through the indoor fan; evaporator thermistor as a measurement
- **Timing:** a priori.
- **Evidence:** at minimum speed, going from full to 40 % fan takes the coil from 19 °C (dry) to
  15 °C (wet), SHR 1.00 → 0.76, for a 5 % loss of EER.

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
  (108 m³/h) and P_d = 30 W, i.e. a specific fan power of 1 000 W/(m³/s) at design flow, typical
  of a small ventilation unit.
- **Computation:** ventilation flow is held between controller samples, so fan energy is the
  exact sum of P(V_k) dt.
- **Not modelled:** heat recovery (in summer it would mostly be bypassed).

### D24. Requirement and minimum flows
- **Timing:** a priori.
- **Requirement:** CO2 at or below 1 000 ppm during occupied hours. Reported as occupied hours
  above 1 000 ppm and the occupied maximum.
- **Minimum flow when occupied:** the building component of ASHRAE 62.1 for offices,
  0.3 L/(s·m²) × 20 m² = 0.006 m³/s, to remove pollutants that do not scale with occupancy.
- **Unoccupied:** no flow, except when a strategy ventilates for cooling (D26).
- **Context:** the existing fixed schedule (0.0105 m³/s while occupied) is close to the full
  62.1 rate for two people (2 × 2.5 L/s + 6 L/s = 0.011 m³/s).

### D25. Strategies compared
- **Timing:** a priori. All strategies use the same temperature and humidity control
  (PI + sigma-delta with the H2 cascade on the fan), so differences come from ventilation alone.
- **V0, fixed schedule:** 0.0105 m³/s from 8 to 18 h.
- **V1, demand-controlled, continuous:**
  - PI on measured CO2, set point 900 ppm (100 ppm of headroom below the requirement).
  - Output between the occupied minimum and 0.03 m³/s, with back-calculation anti-windup.
- **V2, demand-controlled, staged:** three fan levels (0.006, 0.012 and 0.024 m³/s), switched
  up at 850 / 950 ppm and down at 750 / 850 ppm, with a 5-minute minimum dwell time on each
  level. This represents the generic form of staged demand-controlled ventilation.
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
  (`dcv_tuning`).
- **Arithmetic correction, made before any run:** the first version of this entry used
  V_inf + V ≈ 0.0175 m³/s. The steady state at 900 ppm requires N G / (C − C_o) = 0.0217 m³/s,
  so τ ≈ 2 500 s and ∂C/∂V ≈ −22 000 ppm per m³/s.

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

## Findings of the final audit

A reviewer with no access to the development history read the code and docs and tried to break
them. It confirmed two defects and three inaccuracies. The fixes follow.

### D20. Minimum condensing lift in the cycle model
- **Timing:** after observation; found in the final audit.
- **Defect:** COP = (η T_ev − a)/(T_o − T_ev + a) had no floor on T_cond − T_ev. At minimum speed
  the COP reached 17.9 at 20 °C outdoors, 77 at 15 °C, and went negative below about 13 °C.
  36 % of the PI's running time in the summer scenario was at COP > 8.
- **Fix:** T_cond ≥ T_ev + 15 K, which represents head-pressure control holding the compressor's
  minimum pressure ratio and the expansion device's minimum pressure difference. Over 0-45 °C outdoors,
  all speeds and fan settings, COP now stays within 2.5-9.3 and is continuous where the floor
  takes over.
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

### D21. Sigma-delta ripple design used the wrong capacity
- **Timing:** after observation; found in the final audit.
- **Defect:** T_ref was computed with the sensible capacity at full speed (wet coil, 1 685 W).
  While cycling, the unit runs at minimum speed with a dry coil and delivers 1.5 times more
  sensible cooling per unit of demand. The measured ripple was therefore 0.74 K instead of the
  0.5 K claimed.
- **Fix:** use the sensible output per unit of demand at minimum speed; T_ref = 210 s.
- **Checks:**
  - Measured median ripple over about 100 low-load cycles: 0.50 K. A test now checks this.
  - The MPC margin sweep was repeated: 0.6 K still matches the PI's comfort (0.009 K·h).

### D22. Documentation and test corrections
- The coupled air-buffer moisture mode is about 400 s, not 900 s (≈ 820 s is the buffer side
  alone). The test bound was tightened from 300-1000 s to 350-450 s.
- The room dew point under PI control is 17.6 °C, not about 16 °C. The coil at minimum speed
  (18.4-19.1 °C) is still dry, but with a margin of only about 1 K at night.
- The sigma-delta mean-demand test tolerance went from 0.01 to 0.003. Measured errors are
  ≤ 0.0017; at v = 0.05 the old tolerance allowed 20 %.
- The evaporator thermistor in fan-only mode now reads the wet-bulb temperature only if the coil
  is wet, and the return-air temperature if it is dry.

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

---

## Test tolerance audit

All 177 numerical assertions were reviewed and sorted into the four kinds described in
[`tests/README.md`](../tests/README.md). Seven had tolerances looser than they needed to be, or no
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
  (surface above the dew point). Physically the film would evaporate. Under PI control the coil
  is dry and receives no condensate, so the film is usually empty. Under on/off control the coil
  runs wet. The effect on reported numbers is expected to be small but has not been quantified.
- The MPC's Kalman filter estimates the last sample's cooling with full fan speed, no start-up
  ramp and the forecast outdoor temperature. The disturbance state absorbs the difference.

- Single zone, well-mixed air; no radiant temperature or operative-temperature comfort.
- The MPC's model has no humidity, so it cannot trade sensible against latent load.
- MPC disturbance forecasts are held over each 15-minute step.
- The Kalman filter needs about a day to separate a wrong mass temperature from a heat
  disturbance with a single air sensor.
- Unit parameters (UA values, coil holdup, start-up constant) are plausible for a 2.5 kW split
  unit but not fitted to a specific product.
