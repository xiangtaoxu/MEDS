# Fast calibration: example04_column_biophysics

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 3.84 | 3.55 | -8% |
| le | 1.56 | 1.56 | +0% |
| h | 4.18 | 4.35 | +4% |
| gpp | 2.14 | 1.55 | -28% |
| ustar | 1.72 | 0.67 | -61% |

RMSE in units of each target's sigma. Cost: default 62022.3, calibrated 52501.3.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 57.99 | 54.1-61.54 | 45 | +1.36 | 0.43 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 31.97 | 31.06-32.91 | 41.01 | -0.50 | 0.06 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 2.83 | 2.735-2.927 | 2.797 | +0.02 | 0.07 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.001304 | 0.0007969-0.002224 | 0.01 | -1.87 | 0.49 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.05334 | 0.05149-0.05529 | 0.13 | -2.61 | 0.11 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| wstress_sref_stomata | effective | plant_type | 5.345 | 3.571-6.943 | 2 | +1.84 | 0.90 |  | the PFT file's value, a factor 2 either way (Sabot et al. 2022 fit 0.5-5 across  |
| kappa | observation | observation | 0.7254 | 0.6901-0.76 | 0.65 | +0.74 | 0.34 |  | best-practice plan §3.3: the tower's respiration is about the soil chambers' alo |

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 0.97 | 1.42 | 0.988 | 3792 |
| le | 0.99 | 1.26 | 0.934 | 6691 |
| h | 1.60 | 3.00 | 0.778 | 1782 |
| gpp | 0.98 | 1.20 | 0.998 | 1472 |
| ustar | 0.28 | 1.00 | 0.938 | 1782 |

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 0.99, 3 h 0.99, 4 h 0.99, 5 h 0.99, 6 h 0.99, 7 h 0.98, 8 h 0.98, 9 h 0.98, 10 h 0.99, 11 h 0.99, 12 h 0.99, 13 h 0.99, 14 h 1.00, 15 h 1.00, 16 h 0.99, 17 h 0.99, 18 h 0.99, 19 h 0.99, 20 h 0.99, 21 h 0.99, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 0.99, SW quartile 4 0.99

le, model/tower by local hour: 6 h 1.85, 7 h 1.09, 8 h 1.01, 9 h 0.88, 10 h 0.87, 11 h 0.87, 12 h 0.88, 13 h 0.91, 14 h 0.95, 15 h 0.99, 16 h 1.08, 17 h 1.13, 18 h 0.96
le, by light: SW quartile 1 1.15, SW quartile 2 0.96, SW quartile 3 0.91, SW quartile 4 0.91

h, model/tower by local hour: 6 h 0.41, 7 h 0.32, 8 h 0.49, 9 h 0.55, 10 h 0.62, 11 h 0.71, 12 h 0.75, 13 h 0.78, 14 h 0.82, 15 h 0.89, 16 h 1.05, 17 h 2.21, 18 h -6.21
h, by light: SW quartile 1 6.51, SW quartile 2 0.87, SW quartile 3 0.75, SW quartile 4 0.70

gpp, model/tower by local hour: 6 h 0.68, 7 h 1.04, 8 h 0.98, 9 h 0.96, 10 h 0.99, 11 h 0.99, 12 h 0.96, 13 h 1.00, 14 h 1.02, 15 h 1.06, 16 h 1.07, 17 h 1.03, 18 h 0.63
gpp, by light: SW quartile 1 0.95, SW quartile 2 0.95, SW quartile 3 0.97, SW quartile 4 1.09

ustar, model/tower by local hour: 6 h 0.95, 7 h 0.90, 8 h 0.87, 9 h 0.86, 10 h 0.91, 11 h 0.94, 12 h 0.94, 13 h 0.95, 14 h 0.95, 15 h 0.96, 16 h 0.96, 17 h 0.99, 18 h 1.03
ustar, by light: SW quartile 1 1.01, SW quartile 2 0.96, SW quartile 3 0.91, SW quartile 4 0.91

## kappa

kappa = 0.725: the tower's respiration 3.28 -> 4.52 umol m-2 s-1; GPP 7.54 -> 8.28 (2.86 -> 3.14 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar (GPP u* >= 0.325 (the fit used 0.4)): largest linear shift 0.50 sd
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 17.66 sd -- refitted: leaf_angle_mean 57.99->51.05, vcmax25 31.97->27.26, stomatal_g1 2.83->4.957, stomatal_g0 0.001304->0.0003557, z0m_ratio 0.05334->0.04894, wstress_sref_stomata 5.345->0.8407, kappa 0.7254->0.7393
- alternative partitioning: the site TOML declares no daytime partitioning (GPP_DT, RECO_DT)
