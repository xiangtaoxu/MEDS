# Fast calibration: example04_column_biophysics

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 3.81 | 3.54 | -7% |
| le | 1.65 | 1.58 | -4% |
| h | 4.16 | 4.28 | +3% |
| gpp | 2.05 | 1.44 | -30% |
| ustar | 1.73 | 0.67 | -61% |

RMSE in units of each target's sigma. Cost: default 61251.3, calibrated 51376.6.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 51.29 | 47.14-55.26 | 45 | +0.64 | 0.42 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 31.37 | 30.43-32.34 | 41.01 | -0.54 | 0.06 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 2.997 | 2.885-3.113 | 2.797 | +0.14 | 0.08 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.003018 | 0.00228-0.003998 | 0.01 | -1.11 | 0.25 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.05155 | 0.04971-0.05348 | 0.13 | -2.72 | 0.11 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| wstress_sref_stomata | effective | plant_type | 1.403 | 1.103-1.773 | 2 | -0.49 | 0.32 |  | the PFT file's value, a factor 2 either way (Sabot et al. 2022 fit 0.5-5 across  |
| kappa | observation | observation | 0.7238 | 0.6885-0.7584 | 0.65 | +0.72 | 0.34 |  | best-practice plan §3.3: the tower's respiration is about the soil chambers' alo |

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 0.97 | 1.40 | 0.988 | 3792 |
| le | 0.98 | 1.27 | 0.908 | 6691 |
| h | 1.62 | 3.00 | 0.772 | 1782 |
| gpp | 0.98 | 1.17 | 0.986 | 1472 |
| ustar | 0.29 | 1.00 | 0.926 | 1782 |

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 0.99, 3 h 0.99, 4 h 0.99, 5 h 0.99, 6 h 0.99, 7 h 0.98, 8 h 0.98, 9 h 0.98, 10 h 0.99, 11 h 0.99, 12 h 0.99, 13 h 0.99, 14 h 1.00, 15 h 1.00, 16 h 0.99, 17 h 0.99, 18 h 0.99, 19 h 0.99, 20 h 0.99, 21 h 0.99, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 0.99, SW quartile 4 0.99

le, model/tower by local hour: 6 h 1.61, 7 h 1.02, 8 h 0.99, 9 h 0.87, 10 h 0.86, 11 h 0.86, 12 h 0.85, 13 h 0.88, 14 h 0.92, 15 h 0.96, 16 h 1.04, 17 h 1.07, 18 h 0.84
le, by light: SW quartile 1 1.09, SW quartile 2 0.94, SW quartile 3 0.89, SW quartile 4 0.88

h, model/tower by local hour: 6 h 0.41, 7 h 0.31, 8 h 0.47, 9 h 0.54, 10 h 0.62, 11 h 0.71, 12 h 0.75, 13 h 0.78, 14 h 0.81, 15 h 0.88, 16 h 1.03, 17 h 2.19, 18 h -6.15
h, by light: SW quartile 1 6.45, SW quartile 2 0.85, SW quartile 3 0.74, SW quartile 4 0.70

gpp, model/tower by local hour: 6 h 0.68, 7 h 1.03, 8 h 0.97, 9 h 0.95, 10 h 0.97, 11 h 0.97, 12 h 0.95, 13 h 0.98, 14 h 1.00, 15 h 1.06, 16 h 1.07, 17 h 1.03, 18 h 0.61
gpp, by light: SW quartile 1 0.95, SW quartile 2 0.95, SW quartile 3 0.96, SW quartile 4 1.07

ustar, model/tower by local hour: 6 h 0.93, 7 h 0.89, 8 h 0.86, 9 h 0.85, 10 h 0.90, 11 h 0.92, 12 h 0.93, 13 h 0.94, 14 h 0.94, 15 h 0.94, 16 h 0.95, 17 h 0.97, 18 h 1.01
ustar, by light: SW quartile 1 0.99, SW quartile 2 0.94, SW quartile 3 0.90, SW quartile 4 0.89

## kappa

kappa = 0.724: the tower's respiration 3.28 -> 4.53 umol m-2 s-1; GPP 7.54 -> 8.29 (2.86 -> 3.14 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix.
- alternative gpp_ustar (GPP u* >= 0.325 (the fit used 0.4)): largest linear shift 0.53 sd
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 16.75 sd -- refitted: leaf_angle_mean 51.29->49.52, vcmax25 31.37->27.9, stomatal_g1 2.997->4.586, stomatal_g0 0.003018->0.0006371, z0m_ratio 0.05155->0.04811, wstress_sref_stomata 1.403->0.1117, kappa 0.7238->0.7539
- alternative partitioning: the site TOML declares no daytime partitioning (GPP_DT, RECO_DT)
