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

- Single zone, well-mixed air; no radiant temperature or operative-temperature comfort.
- The MPC's model has no humidity, so it cannot trade sensible against latent load.
- MPC disturbance forecasts are held over each 15-minute step.
- The Kalman filter needs about a day to separate a wrong mass temperature from a heat
  disturbance with a single air sensor.
- Unit parameters (UA values, coil holdup, start-up constant) are plausible for a 2.5 kW split
  unit but not fitted to a specific product.
