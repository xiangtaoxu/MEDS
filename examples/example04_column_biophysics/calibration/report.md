# Fast calibration: example04_column_biophysics

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 3.82 | 3.51 | -8% |
| le | 1.63 | 1.60 | -1% |
| h | 4.18 | 4.32 | +3% |
| gpp | 2.34 | 1.59 | -32% |
| ustar | 1.75 | 0.67 | -62% |

RMSE in units of each target's sigma. Cost: default 63657.2, calibrated 52187.2.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 54.93 | 50.97-58.63 | 45 | +1.02 | 0.42 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 32.33 | 31.42-33.28 | 41.01 | -0.48 | 0.06 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 2.82 | 2.747-2.894 | 2.797 | +0.02 | 0.05 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.0003342 | 0.0002012-0.0007137 | 0.01 | -3.45 | 0.97 | beyond 2 prior sd, pushed by lw_up: diagnose it (plan §5.1) or relabel it effective | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.0528 | 0.05097-0.05473 | 0.13 | -2.65 | 0.11 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| wstress_sref_stomata | effective | plant_type | 2.028 | 0.9615-3.739 | 2 | +0.02 | 1.00 |  | the PFT file's value, a factor 2 either way (Sabot et al. 2022 fit 0.5-5 across  |
| kappa | observation | observation | 0.7336 | 0.6979-0.7683 | 0.65 | +0.82 | 0.35 |  | best-practice plan §3.3: the tower's respiration is about the soil chambers' alo |

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 0.97 | 1.39 | 0.988 | 3792 |
| le | 0.99 | 1.29 | 0.910 | 6691 |
| h | 1.49 | 3.00 | 0.803 | 1782 |
| gpp | 1.00 | 1.18 | 0.990 | 1472 |
| ustar | 0.28 | 1.00 | 0.941 | 1782 |

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 0.99, 3 h 0.99, 4 h 0.99, 5 h 0.99, 6 h 0.99, 7 h 0.98, 8 h 0.98, 9 h 0.98, 10 h 0.99, 11 h 0.99, 12 h 0.99, 13 h 1.00, 14 h 1.00, 15 h 1.00, 16 h 0.99, 17 h 0.99, 18 h 0.99, 19 h 0.99, 20 h 0.99, 21 h 0.99, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 0.99, SW quartile 4 0.99

le, model/tower by local hour: 6 h 1.80, 7 h 1.03, 8 h 0.96, 9 h 0.84, 10 h 0.83, 11 h 0.84, 12 h 0.85, 13 h 0.88, 14 h 0.94, 15 h 0.99, 16 h 1.10, 17 h 1.16, 18 h 0.97
le, by light: SW quartile 1 1.15, SW quartile 2 0.95, SW quartile 3 0.89, SW quartile 4 0.88

h, model/tower by local hour: 6 h 0.40, 7 h 0.42, 8 h 0.53, 9 h 0.59, 10 h 0.66, 11 h 0.74, 12 h 0.77, 13 h 0.81, 14 h 0.84, 15 h 0.90, 16 h 1.06, 17 h 2.23, 18 h -6.29
h, by light: SW quartile 1 6.64, SW quartile 2 0.89, SW quartile 3 0.77, SW quartile 4 0.72

gpp, model/tower by local hour: 6 h 0.69, 7 h 1.04, 8 h 0.98, 9 h 0.96, 10 h 0.97, 11 h 0.97, 12 h 0.95, 13 h 0.98, 14 h 1.01, 15 h 1.06, 16 h 1.07, 17 h 1.03, 18 h 0.63
gpp, by light: SW quartile 1 0.95, SW quartile 2 0.95, SW quartile 3 0.96, SW quartile 4 1.08

ustar, model/tower by local hour: 6 h 0.95, 7 h 0.91, 8 h 0.88, 9 h 0.87, 10 h 0.92, 11 h 0.94, 12 h 0.94, 13 h 0.96, 14 h 0.96, 15 h 0.96, 16 h 0.96, 17 h 0.99, 18 h 1.03
ustar, by light: SW quartile 1 1.01, SW quartile 2 0.96, SW quartile 3 0.92, SW quartile 4 0.91

## kappa

kappa = 0.734: the tower's respiration 3.28 -> 4.47 umol m-2 s-1; GPP 7.54 -> 8.25 (2.86 -> 3.13 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar (GPP u* >= 0.325 (the fit used 0.4)): largest linear shift 0.58 sd
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 20.94 sd -- refitted: leaf_angle_mean 54.93->50.53, vcmax25 32.33->28.97, stomatal_g1 2.82->4.763, stomatal_g0 0.0003342->0.0002139, z0m_ratio 0.0528->0.04967, wstress_sref_stomata 2.028->6.214, kappa 0.7336->0.7202
- alternative partitioning: the site TOML declares no daytime partitioning (GPP_DT, RECO_DT)
