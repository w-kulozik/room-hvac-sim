# How the tests choose their tolerances

Every numerical comparison belongs to one of four kinds. Each kind has its own rule for the
tolerance, so that no threshold is read off a run and written back into the test.

1. **Identities.** Quantities equal by construction: balances closed with integrated flows,
   Q_sens + Q_lat = ṁ Δh, exact discretisation, round trips of the psychrometric functions.
   - Tolerance: round-off, 10⁻¹²-10⁻¹⁵ relative, scaled by the size of the stored quantity.

2. **Method errors.** Comparisons with analytical solutions or with an independent method.
   - Tolerance: the known truncation error of the method.
   - Convergence orders are fitted by least squares and must be within 0.15-0.3 of the
     theoretical order.
   - Where possible a value is *predicted* from the measured order instead of bounded.
     Example: the Simpson residual at Δt = 60 s is predicted from the one at 150 s.

3. **Physical plausibility.** Parameters and operating points of real equipment, such as the
   sensible heat ratio or the coil temperature at rating conditions.
   - Tolerance: the range found for that kind of equipment, written next to the assertion.

4. **Closed-loop requirements.** Comfort, humidity accuracy and robustness, defined once in
   [`criteria.py`](criteria.py) with their rationale.
   - Comparisons between controllers are strict orderings; the simulations are deterministic,
     so any difference is real. The size of the difference is reported in the README, not
     asserted.

Changes made to tests after they failed are listed in [`docs/decisions.md`](../docs/decisions.md),
with whether the test or the code was wrong.
