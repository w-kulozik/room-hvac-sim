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

## Coil, capacity and published data of a real unit

9. U.S. Department of Energy. *EnergyPlus Version 9.6 Engineering Reference*, chapter "Coils".
   <https://bigladdersoftware.com/epx/docs/9-6/engineering-reference/coils.html>, accessed
   11 October 2026.
   - Read: section "Single-Speed Electric DX Air Cooling Coil" (Model Description; Dry Coil
     Conditions) and section "Chilled-Water-Based Air Cooling Coil (Coil:Cooling:Water)"
     (Coil Completely Dry, Completely Wet and Part Wet Part Dry Calculations).
   - Used: for the partially wet coil (decision D29), the structure of the part-wet, part-dry
     calculation: the wet/dry boundary where the surface reaches the dew point of the air, the
     enthalpy-based conductance 1/(c_p,sat/UA_internal + c_p/UA_external) of the wet part, and the
     outlet state from an effective condensing surface with effectiveness 1 − exp(−UA_external/
     (ṁ c_p)). For the capacity (D30), the form of the DX coil model: rated capacity times a
     function of entering wet bulb and outdoor dry bulb and a function of the flow fraction. The
     EnergyPlus water coil iterates on a water temperature that varies along the coil; the
     closed forms here follow from the uniform refrigerant temperature and are this project's
     derivation.
10. EnergyPlus 9.6.0 example input file `testfiles/WindACAuto.idf`.
    <https://raw.githubusercontent.com/NREL/EnergyPlus/v9.6.0/testfiles/WindACAuto.idf>,
    accessed 11 October 2026. The same coefficients of WindACCoolCapFT appear as the
    Curve:Biquadratic example in the Input Output Reference, group "Performance Curves"
    (<https://bigladdersoftware.com/epx/docs/9-6/input-output-reference/group-performance-curves.html>).
    - Read: the curve objects WindACCoolCapFT, WindACCoolCapFFF, WindACEIRFT, WindACEIRFFF and
      WindACPLFFPLR, and the Coil:Cooling:DX:SingleSpeed objects that use them.
    - Used: WindACCoolCapFT (biquadratic in entering wet bulb and outdoor dry bulb,
      coefficients 0.942587793, 0.009543347, 0.000683770, −0.011042676, 0.000005249,
      −0.000009720) for the capacity slopes of D30; WindACCoolCapFFF (0.8 + 0.2 ff, stated for
      0.5-1.5) as the comparison for the flow dependence. These are example curves of a
      fixed-speed window unit.
11. Daikin Europe N.V. *RXM-R Technical data book*, EEDEN20. Read in a copy hosted by a reseller:
    <https://daichi.ru/upload/iblock/3c5/zzgxi8dx3mhawwwzlzisk5ek6rukjysw/RXM-R_EEDEN20_Data%20books_English.pdf>,
    accessed 11 October 2026.
    - Read: pages 19-21 (specifications of FTXM25R + RXM25R), page 27 (rating conditions) and
      page 28 (electrical data), as text extracted from the PDF.
    - Used: for FTXM25R + RXM25R in cooling, capacity 1.30 / 2.50 / 3.20 kW (minimum / nominal /
      maximum), nominal power input 0.56 kW, EER 4.50, SEER 8.65, power in thermostat-off mode
      6 W, and the declared part-load points A-D (capacity 2.50, 1.85, 1.19, 1.17 kW; EER 4.50,
      6.52, 10.17, 16.51 at 35, 30, 25, 20 °C outdoors). Rating conditions 27 °C dry bulb /
      19 °C wet bulb indoors, 35 °C outdoors, 5 m of piping. The extracted table interleaves two
      variants of each indoor unit; at points B and C the two variants differ in the second
      decimal (6.52 or 6.54; 10.17 or 10.14), and the first value is used.
12. Daikin Europe N.V. *Wandgerät Klimatisierung, Technische Daten, FTXM-R*, EEDDE22 (German).
    Read in a copy hosted by a reseller:
    <https://www.klimaanlagen-heizungen.de/media/pdf/91/75/08/Daikin-Perfera-Wandgeraet-FTXM-Datenblatt.pdf>,
    accessed 11 October 2026.
    - Read: page 6 (printed page number; page 3 of the PDF file), specifications of the indoor
      units, as text extracted from the PDF.
    - Used: for FTXM25R, power input in nominal cooling 0.025 kW; airflow in cooling 10.5 / 7.6 /
      5.7 / 4.1 m³/min (high / medium / low / quiet); fan motor output 22 W (decision D31 and the
      comparison in [results.md](results.md#comparison-with-a-real-unit)).

## Classification

13. Cooperative Patent Classification, scheme F24F, version 2026.08.
    <https://www.cooperativepatentclassification.org/sites/default/files/cpc/scheme/F/scheme-F24F.pdf>,
    accessed 11 October 2026.
    - Read: groups F24F 3/14-3/1405, 11/00-11/89 and the indexing codes F24F 2110/00-2140/60.
    - Used: the group titles in [faq.md](faq.md), question 5.
