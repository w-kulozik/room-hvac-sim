# Questions a reviewer may ask

Short answers with the numbers behind them. The studies are in [results.md](results.md), the
decisions in [decisions.md](decisions.md) and the sources in [references.md](references.md).
Numbers marked "exploratory" come from runs that are not part of the scripts.

## 1. The coil at minimum speed: is it really dry, and how would the 64 % change?

It was dry in the first coil model, which gave the whole coil one surface temperature. Under PI
control that surface (18.3-19.4 °C) was above the room dew point (about 17.7 °C), so the coil
removed no water and the room reached 63.7 % RH while occupied. The evaporating temperature
(16.9-18.0 °C) was below the dew point of the return air most of the time, and a real coil is
coldest where the air leaves it, so the "no water" result was an upper bound on the humidity
within the model.

The coil is now partially wet (D29): dry where the air enters, wet downstream of the point where
the surface reaches the dew point. At minimum speed and full fan, once the start-up ramp has
decayed, 31-52 % of the surface is wet. With all the unit changes D29-D32 the PI removes 0.9 kg of
water in three days instead of none, and the occupied RH falls from 63.7 % to 62.4 %; the share of
D29 alone in that change was not computed. The main finding stands: at minimum speed with full
airflow the coil is too warm to hold 55 %.

The model has no product-specific "dry mode"; how commercial units implement theirs was not
examined. If such a mode holds the indoor fan at a fixed low speed, its model equivalent gives
(exploratory, PI with handler C, same three days):

| Indoor fan | Energy (kWh) | RH mean occ. (%) | Water drained (kg) |
|---|---|---|---|
| 100 % (plain PI) | 7.39 | 62.4 | 0.9 |
| 70 % fixed | 7.61 | 58.7 | 3.4 |
| 40 % fixed | 9.00 | 51.5 | 8.3 |
| H1 / H2, fan driven by an RH loop | 8.06 / 8.07 | 55.4 / 55.2 | 5.2 / 5.2 |

A fixed low speed overshoots the target and costs more; the humidity loops use only as much fan
reduction as needed.

## 2. Does slowing the fan reduce capacity, and what does that do to H1, H2 and frost protection?

Yes, since D30. Less air lowers the evaporating temperature, and the compressor capacity falls by
4 % per K of evaporating temperature. The coefficient was set so that the capacity rises with the
entering wet bulb as in the EnergyPlus example curve WindACCoolCapFT. The resulting flow
dependence, at rating conditions, against the EnergyPlus example flow curve WindACCoolCapFFF
(0.8 + 0.2 × flow fraction, stated for 0.5-1.5):

| Airflow | Model | EnergyPlus curve |
|---|---|---|
| 70 % | 0.933 | 0.94 |
| 55 % | 0.876 | 0.91 |
| 40 % | 0.788 | 0.88 (extrapolated) |

The loss depends on the operating point. At minimum compressor speed it is 6 % at 70 % and
13 % at 40 % airflow under rating conditions, and 4 % and 10 % with the room at 25 °C and 60 % RH
and 32 °C outdoors, where the coil is dry at full airflow. The EnergyPlus DX coil multiplies the
rated capacity by one flow curve at every operating point; here the loss follows from the coil
and the compressor.

- **H1 and H2:** both still reach 55 % RH. They cost 9 % more energy than the plain PI (10 %
  before the unit changes D29-D32). The temperature PI makes up for the lost capacity by running
  longer: 36.1 h of compressor operation against 31.6 h for the plain PI. While the compressor
  runs the fan averages 56 % and never reaches its minimum.
- **Frost protection:** the evaporating temperature stays above 12.5 °C in the H1 and H2 runs,
  far from the 0.5 °C limit. At full compressor speed and 40 % fan at rating conditions it is
  3.9 °C; before D30, with a capacity that did not fall, it was −2.8 °C and the protection had to
  act (D7). It is still reached with cool, dry return air (22 °C, 0.008 kg/kg, full speed, 40 %
  fan), and a test covers that case. The sections of the EnergyPlus Engineering Reference that
  were read do not cover frost protection of the indoor coil, so no comparison is made on that
  point.

## 3. Why T_t = T_i/2 in the back-calculation anti-windup, and why sigma-delta instead of hysteresis at low load?

**Tracking time.** All PIs in the code use T_t = T_i/2 (the default in `DiscretePI.update`,
`src/roomsim/controllers/pi.py`). Under a constant error e with the actuator saturated at u_max,
the integrator settles at I = u_max − K_p e (1 − T_t/T_i), so the unsaturated demand is
v = u_max + K_p e T_t/T_i. If the error then falls faster than the integrator can follow (T_t),
the output leaves saturation when v drops below u_max: with T_t = T_i the integrator sits at the limit
whatever the error, and that happens only once the error has changed sign; with T_t = T_i/2 it
happens as soon as the error has fallen below half its value under saturation. The value is a
choice made here for fast recovery (D12), checked by a test against the settling formula. Åström
and Murray (2008, section 10.4) recommend a longer tracking time when the sensor is noisy, because
a short one lets noise reset the integrator. The simulated sensor has no noise; on a real unit T_t
would have to be checked against it.

**Low-load handler.** Below 30 % speed the compressor must cycle. Hysteresis (handler A) switches
like a thermostat around the set point and freezes the integrator, so the mean temperature is
not controlled. Sigma-delta (handler C) delivers the PI's demand on average, so the integrator
keeps working. In the three-day runs:

| Handler | T rms occ. (K) | Starts | Starts per running hour | Energy (kWh) |
|---|---|---|---|---|
| A, hysteresis | 0.34 | 105 | 3.6 | 6.85 |
| C, sigma-delta | 0.20 | 251 | 8.0 | 7.39 |

C holds the temperature better and its ripple is designed (0.5 K, measured median 0.48 K), but
it starts the compressor more than twice as often and uses 8 % more energy. Which is preferable
depends on whether temperature accuracy or compressor starts matter more; the humidity and MPC
studies use C because they need the PI to keep integrating.

## 4. When would the MPC pay off, and how does the Kalman filter make it offset-free?

**Payoff.** Against a set-back PI the MPC saves 1.6 % energy, and only because the PI starts at
6:00; with a start at 7:00 the PI matches it (D33). The nominal MPC has a perfect forecast and
exact parameters; 20 % parameter errors raise its energy by 1.6 %, a biased forecast changes it
by less than 1 %.
The MPC pays where a fixed schedule cannot follow what changes:

- **Prices that vary in time.** With a peak/off-peak ratio of 4 the tariff-aware MPC costs 6.6 %
  less than the set-back PI, for 16 % more energy (after D29-D32; before, it was level). Planning
  runs show almost no shift below a ratio of about 3.
- **Start times that would otherwise need tuning** per season and per room: the MPC finds the
  start from the forecast; the PI needs its start time set by hand.
- **Not tested here:** variable occupancy known in advance, forecast weather with real errors,
  rooms with more thermal mass, and constraints such as a power limit. These are the cases a
  schedule handles worst, so they are where to look next.

**Offset-free prediction.** The filter estimates an augmented state [T_a, T_m, d], where d is an
unknown heat input to the air modelled as a random walk (`src/roomsim/estimation.py`). Any
persistent mismatch between the model and the room, such as wrong conductances, fan heat, latent
load or a biased forecast, makes the measured air temperature differ from the prediction, and the
filter moves d until the difference disappears. The planner predicts with that d, so in steady
state its prediction matches the measurement and the planned temperature is reached without
offset, as the integral term of a PI would ensure. In this project the PI underneath also has
integral action, so the room temperature itself is offset-free in any case; d matters for the
quality of the plan. With a single air sensor, a wrong mass temperature and d can only be told
apart through the slow mass mode: started 3 K wrong, the estimate converges with a time constant
of about 14 h.

## 5. Classification and prior art for a fan-speed / coil-temperature cascade

Groups checked in the CPC scheme for F24F, version 2026.08 (cooperativepatentclassification.org,
accessed 11 October 2026). Titles as in the scheme; the dots give the level under the parent.

| Group | Title (abridged where marked) | Relevance |
|---|---|---|
| F24F 11/77 | . . . . by controlling the speed of ventilators (under 11/74 "for controlling air flow rate or air velocity", 11/72, 11/70 "Control systems characterised by their outputs") | indoor fan speed as the actuator (H1, H2) |
| F24F 11/86 | . . . by controlling compressors within refrigeration or heat pump circuits (under 11/80 "for controlling the temperature of the supplied air") | compressor speed loop |
| F24F 11/46 | . . Improving electric energy efficiency or saving (under 11/30) | MPC; 11/47 ". . . Responding to energy costs" for the tariff; 11/48 ". . prior to normal operation, e.g. pre-heating or pre-cooling" for the pre-start |
| F24F 11/62 | . characterised by the type of control or by internal processing, e.g. using fuzzy logic, adaptive control or estimation of values | MPC, Kalman filter |
| F24F 11/41, 11/43 | . . Defrosting; Preventing freezing / . . . of indoor units | frost protection |
| F24F 3/1405 | {in which the humidity of the air is exclusively affected by contact with the evaporator of a closed-circuit cooling system or heat pump circuit} (under 3/14 "by humidification; by dehumidification") | dehumidification by the cooling coil |
| F24F 2110/10, 2110/20 | Temperature; Humidity (control inputs relating to air properties) | room temperature and RH as inputs |
| F24F 2140/20 | Heat-exchange fluid temperature (control inputs relating to system states) | evaporator thermistor |
| F24F 2140/30 | Condensation of water from cooled air | |
| F24F 2130/10 | Weather information or forecasts | MPC forecast |

The IPC counterparts were not checked.

**Search for a claim** of the form "a humidity controller sets a target coil (evaporator)
temperature, and an inner loop adjusts the indoor fan speed to reach it, while the compressor is
controlled on room temperature": combine F24F 11/77 with F24F 2140/20 and F24F 2110/20 (and
F24F 11/86 for the compressor loop), alternatively F24F 3/1405, with keywords such as indoor fan
speed, evaporator or coil temperature, target or set value, dehumidification, relative humidity,
cascade. No search has been carried out for this project, so nothing is claimed about novelty.
The results here suggest what a claim would have to rest on: in the model the cascade's
advantage over a direct RH loop on the fan is small (RH standard deviation 0.60 against 0.88 %
with the same mean and energy).
