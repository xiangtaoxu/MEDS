# Fast calibration: example_flux_tower_bci, variant interception_off

## Gates

| gate | pass | note |
|---|---|---|
| G3 | True | every fitted key has a non-zero, smooth gradient column (the triage fixes the others) |
| G4 | True |  |
| G5 | True | every key near a bound is listed with the target that pushed it |
| G7 | None | the full record with the slow tier on: run the calibrated configs (closed budgets; dry-season GPP and LE no worse than t |
| G10 | True | every declared alternative's shift is under 1 posterior sd, or its refit is reported |
| G12 | True | no scored output has a NaN: a trial with one fails (trials.finish) and never scores |
| G13 | False | a trait key more than 2 prior sd from its evidence: diagnose it (plan §5.1) or relabel it effective |

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | prior source |
|---|---|---|---|---|---|---|---|---|
| leaf_reflect_nir | trait | plant_type | 0.3204 | 0.3172-0.3238 | 0.45 | -3.09 | 0.11 | CLM5 broadleaf evergreen tropical 0.45; measured broadleaf NIR reflectance spans |
| leaf_angle_mean | trait | plant_type | 60.77 | 59.41-62.06 | 45 | +1.70 | 0.17 | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 30.6 | 29.79-31.44 | 41.01 | -0.59 | 0.06 | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 3.525 | 3.233-3.838 | 2.797 | +0.47 | 0.18 | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.002856 | 0.002184-0.00374 | 0.01 | -1.15 | 0.24 | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.05161 | 0.04974-0.05358 | 0.13 | -2.72 | 0.11 | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| wstress_sref_stomata | effective | plant_type | 1.153 | 0.9903-1.342 | 2 | -0.74 | 0.20 | the PFT file's value, a factor 2 either way (Sabot et al. 2022 fit 0.5-5 across  |
| stomata_psi_onset | trait | plant_type | -0.1559 | -0.2744--0.08768 | -0.8571 | +2.64 | 0.84 | half the PFT's turgor-loss point (the model's default), +-0.5 MPa |
| kappa | observation | observation | 0.7537 | 0.7164-0.7895 | 0.65 | +1.02 | 0.37 | best-practice plan §3.3: the tower's respiration is about the soil chambers' alo |

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| albedo | 0.98 | 1.30 | 0.978 | 1295 |
| lw_up | 0.97 | 1.31 | 0.990 | 3792 |
| le | 1.00 | 1.22 | 0.926 | 6691 |
| h | 1.40 | 3.00 | 0.862 | 1782 |
| gpp | 0.99 | 1.17 | 0.987 | 1472 |
| ustar | 0.30 | 1.00 | 0.938 | 1782 |

albedo, model/tower by local hour: 7 h 0.86, 8 h 0.93, 9 h 1.00, 10 h 0.99, 11 h 0.98, 12 h 0.98, 13 h 0.97, 14 h 0.99, 15 h 1.00, 16 h 0.98
albedo, by light: SW quartile 1 0.92, SW quartile 2 0.98, SW quartile 3 1.01, SW quartile 4 1.02

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 0.99, 3 h 0.99, 4 h 0.99, 5 h 0.99, 6 h 0.99, 7 h 0.98, 8 h 0.98, 9 h 0.99, 10 h 0.99, 11 h 0.99, 12 h 1.00, 13 h 1.00, 14 h 1.00, 15 h 1.00, 16 h 1.00, 17 h 0.99, 18 h 0.99, 19 h 0.99, 20 h 0.99, 21 h 0.99, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 0.99, SW quartile 4 0.99

le, model/tower by local hour: 6 h 1.28, 7 h 0.95, 8 h 0.95, 9 h 0.87, 10 h 0.88, 11 h 0.89, 12 h 0.89, 13 h 0.92, 14 h 0.95, 15 h 0.98, 16 h 1.03, 17 h 1.03, 18 h 0.78
le, by light: SW quartile 1 0.98, SW quartile 2 0.92, SW quartile 3 0.91, SW quartile 4 0.93

h, model/tower by local hour: 6 h 0.32, 7 h 0.73, 8 h 0.62, 9 h 0.65, 10 h 0.70, 11 h 0.78, 12 h 0.81, 13 h 0.85, 14 h 0.89, 15 h 0.97, 16 h 1.14, 17 h 2.39, 18 h -6.60
h, by light: SW quartile 1 7.79, SW quartile 2 0.98, SW quartile 3 0.83, SW quartile 4 0.76

gpp, model/tower by local hour: 6 h 0.70, 7 h 1.04, 8 h 0.96, 9 h 0.95, 10 h 0.97, 11 h 0.98, 12 h 0.96, 13 h 0.99, 14 h 1.00, 15 h 1.05, 16 h 1.06, 17 h 1.04, 18 h 0.63
gpp, by light: SW quartile 1 0.95, SW quartile 2 0.93, SW quartile 3 0.96, SW quartile 4 1.08

ustar, model/tower by local hour: 6 h 0.95, 7 h 0.91, 8 h 0.89, 9 h 0.87, 10 h 0.92, 11 h 0.94, 12 h 0.94, 13 h 0.95, 14 h 0.95, 15 h 0.95, 16 h 0.96, 17 h 0.98, 18 h 1.02
ustar, by light: SW quartile 1 1.01, SW quartile 2 0.96, SW quartile 3 0.91, SW quartile 4 0.90

## kappa

kappa = 0.754: the tower's respiration 3.28 -> 4.35 umol m-2 s-1; GPP 7.54 -> 8.61 (2.86 -> 3.26 kgC m-2 yr-1), means over the records with GPP and its respiration measured (daytime and night).

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar (GPP u* >= 0.325 (the fit used 0.4)): largest linear shift 0.36 sd
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 15.10 sd -- refitted: leaf_reflect_nir 0.3204->0.3207, leaf_angle_mean 60.77->60.65, vcmax25 30.6->26.39, stomatal_g1 3.525->4.74, stomatal_g0 0.002856->0.000575, z0m_ratio 0.05161->0.04996, wstress_sref_stomata 1.153->0.1552, stomata_psi_onset -0.1559->-0.1029, kappa 0.7537->0.7766
- alternative partitioning: the site TOML declares no daytime partitioning (GPP_DT, RECO_DT)

Validation cost: default 41156.5, MAP 17781.3.
