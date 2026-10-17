# Fast calibration: example04_column_biophysics

## Validation

**Passed**: the calibrated cost below the default's, and no target's RMSE more than 10% worse.

| target | default | calibrated | change |
|---|---|---|---|
| lw_up | 3.84 | 3.56 | -7% |
| le | 1.67 | 1.60 | -4% |
| h | 3.82 | 3.95 | +3% |
| gpp | 2.40 | 1.58 | -34% |
| ustar | 1.75 | 0.66 | -62% |

RMSE in units of each target's sigma. Cost: default 62356.2, calibrated 50904.9.

The full record with the slow tier on is the next check: run the calibrated configs (budgets closed; dry-season GPP and LE no worse than the default's).

## Keys

| key | kind | scope | MAP | 68 % | prior centre | prior z | sd ratio | note | prior source |
|---|---|---|---|---|---|---|---|---|---|
| leaf_angle_mean | trait | plant_type | 56.53 | 52.79-60 | 45 | +1.20 | 0.41 |  | the PFT file's value, +-10 degrees; spherical leaves average 57.3. The range is  |
| vcmax25 | trait | plant_type | 31.95 | 31.05-32.88 | 39.86 | -0.45 | 0.06 |  | EEO: the coordination of the Rubisco- and light-limited rates at the site's grow |
| stomatal_g1 | trait | plant_type | 2.785 | 2.715-2.857 | 2.798 | -0.01 | 0.05 |  | EEO: the least-cost optimum at the site's growing-season daytime climate (priors |
| stomatal_g0 | trait | plant_type | 0.00028 | 0.0001793-0.0005603 | 0.01 | -3.74 | 0.98 | beyond 2 prior sd, pushed by lw_up: diagnose it (plan §5.1) or relabel it effective | the PFT file's value, a factor e either way: the residual conductance is poorly  |
| z0m_ratio | effective | site | 0.05378 | 0.05202-0.05563 | 0.13 | -2.59 | 0.10 |  | the base value (ED2 0.13), +-0.04: closed canopies give z0/h of about 0.06-0.15  |
| kappa | observation | observation | 0.7434 | 0.7068-0.7788 | 0.65 | +0.92 | 0.36 |  | best-practice plan §3.3: the tower's respiration is about the soil chambers' alo |

Fixed at their priors: wstress_sref_stomata (uninformed: posterior/prior sd ratio 0.95 >= 0.9)

## Targets

| target | chi2/n | sigma scale | model/tower | n (cal) |
|---|---|---|---|---|
| lw_up | 0.98 | 1.44 | 0.988 | 3792 |
| le | 0.99 | 1.29 | 0.916 | 6691 |
| h | 1.04 | 2.94 | 0.861 | 1782 |
| gpp | 1.00 | 1.17 | 0.985 | 1472 |
| ustar | 0.27 | 1.00 | 0.952 | 1782 |

lw_up, model/tower by local hour: 0 h 0.98, 1 h 0.98, 2 h 0.98, 3 h 0.98, 4 h 0.98, 5 h 0.98, 6 h 0.98, 7 h 0.98, 8 h 0.98, 9 h 0.99, 10 h 0.99, 11 h 0.99, 12 h 0.99, 13 h 1.00, 14 h 1.00, 15 h 1.00, 16 h 0.99, 17 h 0.99, 18 h 0.99, 19 h 0.99, 20 h 0.98, 21 h 0.98, 22 h 0.98, 23 h 0.98
lw_up, by light: SW quartile 2 0.98, SW quartile 3 0.99, SW quartile 4 0.99

le, model/tower by local hour: 6 h 1.80, 7 h 1.04, 8 h 0.98, 9 h 0.85, 10 h 0.83, 11 h 0.85, 12 h 0.85, 13 h 0.89, 14 h 0.95, 15 h 0.99, 16 h 1.10, 17 h 1.16, 18 h 0.99
le, by light: SW quartile 1 1.15, SW quartile 2 0.96, SW quartile 3 0.89, SW quartile 4 0.88

h, model/tower by local hour: 6 h 0.49, 7 h 0.51, 8 h 0.65, 9 h 0.71, 10 h 0.77, 11 h 0.84, 12 h 0.85, 13 h 0.87, 14 h 0.87, 15 h 0.91, 16 h 1.00, 17 h 1.80, 18 h -4.10
h, by light: SW quartile 1 5.20, SW quartile 2 0.93, SW quartile 3 0.85, SW quartile 4 0.79

gpp, model/tower by local hour: 6 h 0.69, 7 h 1.04, 8 h 0.97, 9 h 0.95, 10 h 0.97, 11 h 0.97, 12 h 0.94, 13 h 0.98, 14 h 1.00, 15 h 1.05, 16 h 1.07, 17 h 1.04, 18 h 0.64
gpp, by light: SW quartile 1 0.95, SW quartile 2 0.94, SW quartile 3 0.96, SW quartile 4 1.07

ustar, model/tower by local hour: 6 h 0.93, 7 h 0.91, 8 h 0.90, 9 h 0.89, 10 h 0.94, 11 h 0.96, 12 h 0.96, 13 h 0.97, 14 h 0.97, 15 h 0.97, 16 h 0.97, 17 h 0.99, 18 h 1.02
ustar, by light: SW quartile 1 1.01, SW quartile 2 0.97, SW quartile 3 0.93, SW quartile 4 0.92

## kappa

kappa = 0.743: the tower's respiration 3.28 -> 4.41 umol m-2 s-1; GPP 7.54 -> 8.21 (2.86 -> 3.11 kgC m-2 yr-1), means over the records with GPP and its respiration measured, day and night; GPP gains the scaled respiration by day only.

## Uncertainty

Laplace, from the final gradient matrix -- LOCAL ONLY (the linearity check failed).
- alternative gpp_ustar (GPP u* >= 0.325 (the fit used 0.4)): largest linear shift 0.56 sd
- alternative closure (Bowen (H and LE scaled together) against the fit's attribution shares): largest linear shift 22.02 sd -- refitted: leaf_angle_mean 56.53->49.96, vcmax25 31.95->28.74, stomatal_g1 2.785->4.83, stomatal_g0 0.00028->0.0003499, z0m_ratio 0.05378->0.05034, kappa 0.7434->0.7235
- alternative partitioning: the site TOML declares no daytime partitioning (GPP_DT, RECO_DT)
