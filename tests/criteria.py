"""Acceptance criteria for closed-loop tests, defined once with their rationale.

Closed-loop results depend on tuning, so their tests check requirements stated here rather than
thresholds read off a run. Measured values are reported in the README, not encoded in the tests.
Orderings between controllers (A uses less energy than B) are tested as strict inequalities:
the simulations are deterministic, so any difference is real; its size is documented elsewhere.
"""

# Comfort: degree-hours outside the occupied band (25 +- 1 degC) per simulated day.
# 0.05 K h per day is 3 minutes at 1 K beyond the band, or 30 minutes at 0.1 K: below what an
# occupant would notice, and well below the hourly resolution of comfort standards.
KH_PER_DAY = 0.05

# Humidity control: mean occupied RH within 1 %RH of its set point. Room humidity sensors are
# typically accurate to +-2-3 %RH, so a 1 %RH offset is below what the controller can resolve.
RH_ACCURACY = 1.0

# Robustness of the MPC: with model parameters off by 20 %, energy within 5 % of the nominal run
# and the comfort criterion still met.
ROBUST_ENERGY_REL = 0.05

# A humidity difference counts as physically meaningful above 1 %RH, which at 25 degC is a
# humidity-ratio difference of about 2e-4 kg/kg.
RH_RESOLUTION_W = 2e-4
