# References

Publications named in the documentation and the code. For each one: what was read, and what is
taken from it. Nothing else in the project relies on these sources.

## Psychrometrics and room model

1. Alduchov, O. A. and Eskridge, R. E. (1996). Improved Magnus form approximation of saturation
   vapor pressure. *Journal of Applied Meteorology* 35(4), 601-609.
   doi:10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2
   - Read: the report version, US Department of Energy report DOE/ER/61011-T6 (1997),
     OSTI 548871, doi:10.2172/548871.
   - Used: the coefficients of the recommended Magnus form over water,
     e_w = 6.1094 hPa · exp(17.625 t/(243.04 + t)), equation (25) of the report, recommended
     for −40 to 50 °C with a maximum relative error below 0.4 %.
2. U.S. Department of Energy. *EnergyPlus Version 9.6 Engineering Reference*, section
   "Effective Moisture Penetration Depth (EMPD) Model".
   <https://bigladdersoftware.com/epx/docs/9-6/engineering-reference/effective-moisture-penetration-depth-empd.html>
   - Read: that section.
   - Used: the structure of the moisture buffer, a layer whose capacity is density × depth ×
     sorption slope, exchanging with the air through a boundary-layer resistance (from the
     Lewis relation) in series with diffusion through half the layer depth. The project keeps
     the surface layer only; the deep layer and the coating are not modelled.

## Occupants and ventilation

3. Persily, A. and de Jonge, L. (2017). Carbon dioxide generation rates for building occupants.
   *Indoor Air* 27(5), 868-879. doi:10.1111/ina.12383
   - Read: the full text (PubMed Central PMC5666301).
   - Used: that 0.0052 L/s of CO2 for an adult at 1.2 met is the value cited in ASHRAE
     Standard 62.1 and ASTM D6245, and their lower estimate for offices, 0.0048 L/s per person.
4. ASHRAE (2007). *ANSI/ASHRAE Standard 62.1-2007, Ventilation for Acceptable Indoor Air
   Quality*, Table 6-1, Minimum Ventilation Rates in Breathing Zone.
   - Read: not the standard itself. The office-space rates, 5 cfm per person plus 0.06 cfm per
     ft², are taken from ASHRAE Interpretation IC 62.1-2007-08 (2008), which applies them. In SI
     units they are 2.5 L/s per person and 0.3 L/(s·m²).
   - Used: the minimum occupied flow and the reference flow for two people (decision D24).

## Air-conditioning unit

5. Henderson, H. I. and Rengarajan, K. (1996). A model to predict the latent capacity of air
   conditioners and heat pumps at part-load conditions with constant fan operation. *ASHRAE
   Transactions* 102(1).
   - Read: the abstract only (OSTI 392453).
   - Used: the mechanism of part-load latent degradation, condensate held on the coil
     re-evaporating when the compressor stops and the fan keeps running (decision D10). Their
     equations are not used.
6. Henderson, H., Huang, Y. J. and Parker, D. (1999). *Residential equipment part load curves
   for use in DOE-2*. Report LBNL-42175, Lawrence Berkeley National Laboratory.
   - Read: the full report.
   - Used: start-up capacity time constants of 40-60 s as typical of modern residential units,
     and 60 s for their "typical" unit from literature values of 30-80 s (decision D9).

## Control

7. Skogestad, S. (2003). Simple analytic rules for model reduction and PID controller tuning.
   *Journal of Process Control* 13(4), 291-309. doi:10.1016/S0959-1524(02)00062-8
   - Read: the author's summary of the rules on his publication page; the paper itself was
     not read.
   - Used: the SIMC PI rules for a first-order process with time delay,
     K_p = τ/(k(τ_c + θ)) and T_i = min(τ, 4(τ_c + θ)) (`simc_pi` in
     `src/roomsim/controllers/pi.py`).
8. Åström, K. J. and Murray, R. M. (2008). *Feedback Systems: An Introduction for Scientists
   and Engineers*. Princeton University Press. Section 10.4, Integrator windup.
   - Read: section 10.4 of the online edition (version 2.11b).
   - Used: the back-calculation anti-windup scheme, and its advice that the reset gain be a
     fraction of 1/T_i so that measurement noise does not reset the integrator. The tracking
     time T_i/2 used here is not taken from it (decision D12).
