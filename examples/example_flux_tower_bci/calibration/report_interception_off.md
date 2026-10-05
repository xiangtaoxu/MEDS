# Fast calibration: example_flux_tower_bci, variant interception_off

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 3.75 | 3.44 | -8% |
| le | 1.50 | 1.45 | -3% |
| h | 4.05 | 4.26 | +5% |
| gpp | 2.09 | 1.44 | -31% |
| ustar | 1.76 | 0.68 | -61% |

RMSE in units of each target's sigma. Cost: default 60347.4, calibrated 50273.1.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 60.77 | 56.77-64.31 | 45 | +1.70 | 0.47 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 29.17 | 28.38-30 | 41.01 | -0.69 | 0.06 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 3.347 | 3.242-3.456 | 2.797 | +0.36 | 0.06 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.001199 | 0.0007018-0.002166 | 0.01 | -1.95 | 0.55 |  | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.05056 | 0.04869-0.05255 | 0.13 | -2.78 | 0.12 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| wstress_sref_stomata | effective | plant_type | 1.28 | 1.005-1.62 | 2 | -0.61 | 0.32 |  | the PFT file's value, a factor 2 either way (Sabot et al. 2022 fit 0.5-5 across  |
| kappa | observation | observation | 0.7801 | 0.7394-0.8179 | 0.65 | +1.29 | 0.41 |  | best-practice plan §3.3: the tower's respiration is about the soil chambers' alo |

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 0.95 | 1.40 | 0.989 | 3792 |
| le | 0.98 | 1.24 | 0.936 | 6691 |
| h | 1.67 | 3.00 | 0.778 | 1782 |
| gpp | 0.98 | 1.19 | 0.983 | 1472 |
| ustar | 0.31 | 1.00 | 0.924 | 1782 |

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 0.99, 3 h 0.99, 4 h 0.99, 5 h 0.99, 6 h 0.99, 7 h 0.98, 8 h 0.98, 9 h 0.98, 10 h 0.98, 11 h 0.99, 12 h 0.99, 13 h 0.99, 14 h 1.00, 15 h 1.00, 16 h 0.99, 17 h 0.99, 18 h 0.99, 19 h 0.99, 20 h 0.99, 21 h 0.99, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 0.99, SW quartile 4 0.99

le, model/tower by local hour: 6 h 1.24, 7 h 0.97, 8 h 0.96, 9 h 0.88, 10 h 0.89, 11 h 0.90, 12 h 0.90, 13 h 0.93, 14 h 0.97, 15 h 0.99, 16 h 1.06, 17 h 1.07, 18 h 0.75
le, by light: SW quartile 1 1.00, SW quartile 2 0.94, SW quartile 3 0.92, SW quartile 4 0.94

h, model/tower by local hour: 6 h 0.33, 7 h 0.48, 8 h 0.51, 9 h 0.55, 10 h 0.61, 11 h 0.70, 12 h 0.73, 13 h 0.77, 14 h 0.82, 15 h 0.90, 16 h 1.06, 17 h 2.29, 18 h -6.62
h, by light: SW quartile 1 7.39, SW quartile 2 0.89, SW quartile 3 0.74, SW quartile 4 0.68

gpp, model/tower by local hour: 6 h 0.71, 7 h 1.03, 8 h 0.95, 9 h 0.94, 10 h 0.97, 11 h 0.97, 12 h 0.95, 13 h 0.98, 14 h 1.00, 15 h 1.04, 16 h 1.06, 17 h 1.05, 18 h 0.65
gpp, by light: SW quartile 1 0.96, SW quartile 2 0.93, SW quartile 3 0.95, SW quartile 4 1.08

ustar, model/tower by local hour: 6 h 0.94, 7 h 0.90, 8 h 0.87, 9 h 0.85, 10 h 0.90, 11 h 0.92, 12 h 0.93, 13 h 0.94, 14 h 0.94, 15 h 0.94, 16 h 0.95, 17 h 0.97, 18 h 1.01
ustar, by light: SW quartile 1 1.00, SW quartile 2 0.94, SW quartile 3 0.90, SW quartile 4 0.89

## kappa

kappa = 0.780: the tower's respiration 3.28 -> 4.20 umol m-2 s-1; GPP 7.54 -> 8.47 (2.86 -> 3.21 kgC m-2 yr-1), means over the records with GPP and its respiration measured (daytime and night).

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar (GPP u* >= 0.325 (the fit used 0.4)): largest linear shift 0.44 sd
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 18.10 sd -- refitted: leaf_angle_mean 60.77->56.78, vcmax25 29.17->26.82, stomatal_g1 3.347->4.773, stomatal_g0 0.001199->0.0004139, z0m_ratio 0.05056->0.04921, wstress_sref_stomata 1.28->0.198, kappa 0.7801->0.7892
- alternative partitioning: the site TOML declares no daytime partitioning (GPP_DT, RECO_DT)
