# Fast calibration: example04_column_biophysics

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 3.87 | 3.54 | -9% |
| le | 1.70 | 1.59 | -6% |
| h | 4.34 | 4.35 | +0% |
| gpp | 2.50 | 1.58 | -37% |
| ustar | 1.72 | 0.67 | -61% |

RMSE in units of each target's sigma. Cost: default 66200.8, calibrated 52686.9.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 55.12 | 51.23-58.74 | 45 | +1.04 | 0.41 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 31.33 | 30.44-32.24 | 41.01 | -0.54 | 0.06 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 2.821 | 2.743-2.901 | 2.797 | +0.02 | 0.06 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.0004634 | 0.000254-0.001068 | 0.01 | -2.99 | 0.93 | beyond 2 prior sd, pushed by lw_up: diagnose it (plan §5.1) or relabel it effective | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.05357 | 0.05172-0.05552 | 0.13 | -2.60 | 0.11 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| kappa | observation | observation | 0.7669 | 0.7276-0.8039 | 0.65 | +1.15 | 0.39 |  | best-practice plan §3.3: the tower's respiration is about the soil chambers' alo |

Fixed at their priors: wstress_sref_stomata (dead: its gradient column is zero (a key the windows never use))

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 0.96 | 1.42 | 0.988 | 3792 |
| le | 0.98 | 1.29 | 0.922 | 6691 |
| h | 1.56 | 3.00 | 0.789 | 1782 |
| gpp | 0.99 | 1.18 | 0.986 | 1472 |
| ustar | 0.28 | 1.00 | 0.941 | 1782 |

lw_up, model/tower by local hour: 0 h 0.99, 1 h 0.99, 2 h 0.99, 3 h 0.99, 4 h 0.99, 5 h 0.99, 6 h 0.99, 7 h 0.98, 8 h 0.98, 9 h 0.98, 10 h 0.99, 11 h 0.99, 12 h 0.99, 13 h 0.99, 14 h 1.00, 15 h 1.00, 16 h 0.99, 17 h 0.99, 18 h 0.99, 19 h 0.99, 20 h 0.99, 21 h 0.99, 22 h 0.99, 23 h 0.99
lw_up, by light: SW quartile 2 0.99, SW quartile 3 0.99, SW quartile 4 0.99

le, model/tower by local hour: 6 h 1.96, 7 h 1.10, 8 h 1.01, 9 h 0.87, 10 h 0.86, 11 h 0.86, 12 h 0.86, 13 h 0.89, 14 h 0.94, 15 h 0.98, 16 h 1.08, 17 h 1.14, 18 h 1.00
le, by light: SW quartile 1 1.16, SW quartile 2 0.96, SW quartile 3 0.90, SW quartile 4 0.89

h, model/tower by local hour: 6 h 0.41, 7 h 0.34, 8 h 0.49, 9 h 0.56, 10 h 0.63, 11 h 0.72, 12 h 0.76, 13 h 0.80, 14 h 0.83, 15 h 0.90, 16 h 1.05, 17 h 2.21, 18 h -6.14
h, by light: SW quartile 1 6.51, SW quartile 2 0.88, SW quartile 3 0.76, SW quartile 4 0.71

gpp, model/tower by local hour: 6 h 0.71, 7 h 1.05, 8 h 0.97, 9 h 0.95, 10 h 0.97, 11 h 0.97, 12 h 0.94, 13 h 0.98, 14 h 1.00, 15 h 1.05, 16 h 1.07, 17 h 1.05, 18 h 0.66
gpp, by light: SW quartile 1 0.96, SW quartile 2 0.94, SW quartile 3 0.95, SW quartile 4 1.07

ustar, model/tower by local hour: 6 h 0.95, 7 h 0.91, 8 h 0.88, 9 h 0.87, 10 h 0.92, 11 h 0.94, 12 h 0.94, 13 h 0.96, 14 h 0.96, 15 h 0.96, 16 h 0.96, 17 h 0.99, 18 h 1.03
ustar, by light: SW quartile 1 1.01, SW quartile 2 0.96, SW quartile 3 0.92, SW quartile 4 0.91

## kappa

kappa = 0.767: the tower's respiration 3.28 -> 4.27 umol m-2 s-1; GPP 7.54 -> 8.13 (2.86 -> 3.08 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar (GPP u* >= 0.325 (the fit used 0.4)): largest linear shift 0.55 sd
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 19.87 sd -- refitted: leaf_angle_mean 55.12->53.56, vcmax25 31.33->28.23, stomatal_g1 2.821->4.442, stomatal_g0 0.0004634->0.0003678, z0m_ratio 0.05357->0.04945, kappa 0.7669->0.7422
- alternative partitioning: the site TOML declares no daytime partitioning (GPP_DT, RECO_DT)
