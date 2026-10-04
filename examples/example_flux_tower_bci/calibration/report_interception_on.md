# Fast calibration: example_flux_tower_bci, variant interception_on

## Gates

| gate | pass | note |
|---|---|---|
| G3 | True | every fitted key has a non-zero, smooth gradient column (the triage fixes the others) |
| G4 | True |  |
| G5 | True | every key near a bound is listed with the target that pushed it |
| G7 | None | the full record with the slow tier on: run the calibrated configs (closed budgets; dry-season GPP and LE no worse than t |
| G10 | True | every declared alternative's shift is under 1 posterior sd, or its refit is reported |
| G12 | True | no scored output has a NaN: a trial with one fails (trials.finish) and never scores |
| G13 | True | a trait key more than 2 prior sd from its evidence: diagnose it (plan §5.1) or relabel it effective |

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | prior source |
|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 54.02 | 49.92-57.88 | 45 | +0.92 | 0.43 | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 30.94 | 30.03-31.88 | 41.01 | -0.57 | 0.06 | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 2.992 | 2.879-3.108 | 2.797 | +0.13 | 0.08 | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.003154 | 0.002406-0.004136 | 0.01 | -1.07 | 0.24 | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.05165 | 0.04982-0.05356 | 0.13 | -2.71 | 0.11 | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| wstress_sref_stomata | effective | plant_type | 1.23 | 0.9714-1.552 | 2 | -0.66 | 0.31 | the PFT file's value, a factor 2 either way (Sabot et al. 2022 fit 0.5-5 across  |
| kappa | observation | observation | 0.7384 | 0.7014-0.7741 | 0.65 | +0.87 | 0.36 | best-practice plan §3.3: the tower's respiration is about the soil chambers' alo |

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 0.97 | 1.42 | 0.988 | 3792 |
| le | 0.98 | 1.27 | 0.912 | 6691 |
| h | 1.62 | 3.00 | 0.774 | 1782 |
| gpp | 0.98 | 1.18 | 0.987 | 1472 |
| ustar | 0.29 | 1.00 | 0.929 | 1782 |

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 0.99, 3 h 0.99, 4 h 0.99, 5 h 0.99, 6 h 0.98, 7 h 0.98, 8 h 0.98, 9 h 0.98, 10 h 0.99, 11 h 0.99, 12 h 0.99, 13 h 0.99, 14 h 1.00, 15 h 1.00, 16 h 0.99, 17 h 0.99, 18 h 0.99, 19 h 0.99, 20 h 0.99, 21 h 0.99, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 0.99, SW quartile 4 0.99

le, model/tower by local hour: 6 h 1.63, 7 h 1.03, 8 h 0.99, 9 h 0.87, 10 h 0.86, 11 h 0.86, 12 h 0.86, 13 h 0.89, 14 h 0.93, 15 h 0.96, 16 h 1.04, 17 h 1.07, 18 h 0.85
le, by light: SW quartile 1 1.10, SW quartile 2 0.94, SW quartile 3 0.89, SW quartile 4 0.89

h, model/tower by local hour: 6 h 0.41, 7 h 0.30, 8 h 0.47, 9 h 0.54, 10 h 0.62, 11 h 0.71, 12 h 0.75, 13 h 0.78, 14 h 0.82, 15 h 0.89, 16 h 1.04, 17 h 2.21, 18 h -6.19
h, by light: SW quartile 1 6.48, SW quartile 2 0.86, SW quartile 3 0.74, SW quartile 4 0.70

gpp, model/tower by local hour: 6 h 0.69, 7 h 1.04, 8 h 0.96, 9 h 0.95, 10 h 0.97, 11 h 0.97, 12 h 0.95, 13 h 0.98, 14 h 1.00, 15 h 1.05, 16 h 1.07, 17 h 1.04, 18 h 0.62
gpp, by light: SW quartile 1 0.95, SW quartile 2 0.94, SW quartile 3 0.96, SW quartile 4 1.07

ustar, model/tower by local hour: 6 h 0.94, 7 h 0.89, 8 h 0.86, 9 h 0.86, 10 h 0.90, 11 h 0.93, 12 h 0.93, 13 h 0.95, 14 h 0.94, 15 h 0.95, 16 h 0.95, 17 h 0.98, 18 h 1.02
ustar, by light: SW quartile 1 1.00, SW quartile 2 0.94, SW quartile 3 0.91, SW quartile 4 0.90

## kappa

kappa = 0.738: the tower's respiration 3.28 -> 4.44 umol m-2 s-1; GPP 7.54 -> 8.70 (2.86 -> 3.30 kgC m-2 yr-1), means over the records with GPP and its respiration measured (daytime and night).

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar (GPP u* >= 0.325 (the fit used 0.4)): largest linear shift 0.52 sd
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 16.72 sd -- refitted: leaf_angle_mean 54.02->51.01, vcmax25 30.94->27.86, stomatal_g1 2.992->4.507, stomatal_g0 0.003154->0.0007325, z0m_ratio 0.05165->0.04815, wstress_sref_stomata 1.23->0.1104, kappa 0.7384->0.756
- alternative partitioning: the site TOML declares no daytime partitioning (GPP_DT, RECO_DT)

Validation cost: default 26093.6, MAP 17288.4.
