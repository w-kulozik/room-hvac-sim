# roomsim: a dynamic room model for HVAC control

A small simulator of one office room with an inverter split air conditioner. It is verified
against analytical solutions and used to test and compare air-conditioning controllers.

- **The room:** an air node, a thermal-mass node, a moisture balance with buffering by the
  interior surfaces, a CO2 balance and proper psychrometrics.
- **What drives it:** outdoor weather, sun through a window, and occupants coming and going.
- **The unit model:** the cooling coil (dry and wet), the water held on the coil, frost
  protection, part-load efficiency and start-up losses.
- **The controllers:** discrete, as in a real unit. They sample the sensors every minute and
  hold their commands until the next sample.

![Air temperature and compressor speed over one day for the five standard controllers](figures/controllers_day.png)

## Main findings

All results are from simulation of one synthetic summer day, repeated (see Limits below).

1. **A proportional-integral (PI) controller uses 40 % less energy than an on/off thermostat
   but loses humidity control.**
   - The unit is oversized for the room, so the PI runs it at minimum speed with a dry coil.
   - Mean indoor relative humidity while occupied rises from 49 % to 64 %.
2. **Slowing the indoor fan restores dehumidification.** The humidity returns to its 55 % set
   point for 10 % more energy than the plain PI, and still 34 % less than on/off. A cascade
   through the coil temperature holds it a little more steadily.
3. **Model predictive control (MPC) saves little against a well-scheduled PI:** 3 % against a
   PI with set-back.
   - A time-of-use tariff does not pay off here once the cycle efficiency is modelled with a
     realistic minimum lift.
   - An earlier "10 % saving" was traced to an unbounded efficiency at night and withdrawn.
4. **Demand-controlled ventilation pays only when occupancy varies.** A fixed flow sized for
   the peak needs 78 % more air and 8 % more energy. Night flush pays only at about half the
   design flow.
5. **The numerics are verified.**
   - Fourth-order convergence, also in closed loop with the unit.
   - Energy, water and CO2 balances closed to round-off.
   - 151 automated tests.

## Limits

- One zone with well-mixed air, and one sinusoidal summer day repeated; no measured weather.
- The model has not been validated against measurements of a real room. Room and unit
  parameters are plausible for a small office and a 2.5 kW split unit, not fitted to a
  specific building or product.
- The findings are comparisons between controllers on the same model. Absolute energies
  depend on parameters such as the minimum condensing lift and the start-up time constant.
- The full list is at the end of [docs/decisions.md](docs/decisions.md).

## Documentation

Read in this order:

1. [docs/model.md](docs/model.md): the equations, the parameters and the air-conditioning unit.
2. [docs/results.md](docs/results.md): every study with its table and figure. The findings
   above come from here.
3. [docs/verification.md](docs/verification.md): what is checked and how accurately.
4. [docs/decisions.md](docs/decisions.md): the log of modelling and study decisions (D1-D28).
   It records which ones were fixed before running a study, the corrections made after the
   end-of-project check, and the known limitations.
5. [docs/specification.md](docs/specification.md): the full specification and its
   revisions.
6. [tests/README.md](tests/README.md): how test tolerances are chosen.
7. [docs/references.md](docs/references.md): the publications cited, what was read in each and
   what is taken from it.

## Install and run

Python 3.11 or later.

```bash
pip install -e ".[dev]"
pytest
python scripts/run_baseline.py
python scripts/compare_controllers.py
python scripts/compare_humidity.py
python scripts/compare_mpc.py
python scripts/compare_ventilation.py
```

Each script prints its table and writes its figures and its table to `figures/`.

## Layout

```
src/roomsim/
  psychro.py       psychrometric functions
  params.py        RoomParams dataclass, derived constants
  inputs.py        disturbance scenarios d(t), vectorised for forecasts
  plant.py         rhs(t, x, u, d, p), flux breakdown, Jacobian, sensors
  ac_unit.py       inverter split unit: coil, cycle, frost protection, start-up, coil water
  integrators.py   Euler, RK4 (classical fourth-order Runge-Kutta), Euler stability limit
  simulate.py      fixed-step loop, discrete controller with zero-order hold, periodic shooting
  results.py       SimResult: states, psychrometrics, balance residuals
  metrics.py       energy, starts, comfort, humidity and CO2 indicators
  verification.py  analytical solutions and convergence studies (shared by tests and figures)
  plotting.py      figures
  control_model.py linear RC model for the controllers, exact discretisation
  estimation.py    Kalman filter with a disturbance state (offset-free)
  tariff.py        time-of-use electricity price
  ventilation.py   ventilation fan power (affinity law)
  controllers/     on/off thermostat, discrete PI with anti-windup, low-load handlers,
                   humidity-aware fan control (direct and cascade), economic MPC,
                   ventilation strategies (scheduled, CO2 PI, staged, night flush)
scripts/           one script per study (see Install and run)
tests/             151 tests; tolerance rules in tests/README.md and tests/criteria.py
docs/              model, results, verification, decisions, specification, references
figures/           figures and result tables written by the scripts
```

## Roadmap

1. MPC with humidity in the model, so that the planner can trade sensible against latent load.
2. Ventilation inside the MPC, so that night flush and demand control are planned together.

## License

MIT, see [LICENSE](LICENSE).
