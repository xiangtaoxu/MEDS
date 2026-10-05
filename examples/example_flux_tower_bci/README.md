# Forcing from a flux tower: Barro Colorado Island

This example builds a MEDS forcing file from a flux tower's own meteorology, then drives MEDS with
it at the tower, starting from the forest the BCI 50-ha plot census measured. The site is Barro
Colorado Island (BCI), Panama: AmeriFlux PA-Bar, a 41 m tower above a seasonal tropical forest,
2012–2017.

What it shows:
- **declare, then validate.** A site TOML ([`bci_site.toml`](bci_site.toml)) declares what the
  file is. [`make_tower_forcing.py`](../../scripts/prepare_flux_tower/make_tower_forcing.py) checks
  each declaration against the sun and the data, and stops on a disagreement.
- **the model's conversions, not the provider's.** The file stores relative humidity as measured,
  and MEDS turns it into specific humidity with its own saturation curve.
- **explicit, flagged gap filling**, and a test of the longwave fill on observations it never saw.
- **the tower's height in the model.** Every sample is moved from 41 m above the ground to the top
  of each patch's canopy air space.
- **a start from a census, not a spin-up.** The 2010 census of the 50-ha plot is the stand, one patch
  per 20 m quadrat, and MEDS fuses it with its own restructuring before the first step.
- **a calibration of the fast parameters against the same tower**
  ([`scripts/calibrate_fast`](../../scripts/calibrate_fast/README.md)). It is scored on windows the fit
  never saw and shipped in [`calibration/`](calibration), and the figures show the default and the
  calibrated run side by side.

The design records are
[`docs/dev_plans/MEDS_FLUX_TOWER_FORCING_PLAN.md`](../../docs/dev_plans/MEDS_FLUX_TOWER_FORCING_PLAN.md),
[`docs/dev_plans/MEDS_BCI_CENSUS_INIT_PLAN.md`](../../docs/dev_plans/MEDS_BCI_CENSUS_INIT_PLAN.md)
and [`docs/dev_plans/MEDS_FAST_CALIBRATION_BEST_PRACTICE.md`](../../docs/dev_plans/MEDS_FAST_CALIBRATION_BEST_PRACTICE.md).

## The data

`BCI_v5.1.csv` is M. Detto's "Barro Colorado Island – eddy covariance flux data (2012–2017)",
[Zenodo 6456527](https://zenodo.org/records/6456527) (doi:10.5061/dryad.3tx95x6j5), released
under CC0. It is not in the repository. [`fetch_bci_data.py`](fetch_bci_data.py) downloads it into
`data/` and checks it against the published md5; `data/` and `output/` are gitignored. The data's
README asks that publications acknowledge the Center for Tropical Forest Science – Forest Global
Earth Observatory (CTFS-ForestGEO), which supported the tower.

**Note: the file's two longwave columns are swapped.** `Rl_up` is the downwelling longwave and
`Rl_dn` the upwelling, so [`bci_site.toml`](bci_site.toml) reads `LWdown` from `Rl_up`. Three
things show it:
- the provider's `Rnet` equals Rs − Rs_dn + Rl_up − Rl_dn to the last digit;
- at night `Rl_dn` averages 1.022 times the air's blackbody emission σT⁴, more than the sky can
  emit, and `Rl_up` 0.948;
- ERA5-Land's daily longwave follows `Rl_up` (r 0.86) and not `Rl_dn` (0.12).

**The census** is the BCI 50-ha plot's tree table for census 7 (2010), with census 6 (2005) for the
plot's biomass mortality: Condit R., Pérez R., Aguilar S., Lao S., Foster R., Hubbell S.P. 2019,
*Complete data from the Barro Colorado 50-ha plot: 423617 trees, 35 years*, Dryad
[doi:10.15146/5xcp-0d46](https://doi.org/10.15146/5xcp-0d46), CC0. The PIs ask to be told of papers
that use it. [`bci_census.toml`](bci_census.toml) names the copy to read and its checksum; the one it
names is the lab's CSV export, and the public tables are `bci.tree.zip` on the Dryad page, which
must be downloaded in a browser because Dryad refuses scripts.
## Running it

```bash
cd examples/example_flux_tower_bci
python run_example.py --forcing-only        # fetch, build the forcing, score the longwave fill, figures
python run_example.py                       # ... then the census file, the five tower years, the comparison
python run_example.py --calibrate --workers 40   # ... and redo the calibration first (about an hour on one 40-core node)
```

`run_example.py` needs numpy, pandas, netCDF4 and matplotlib, and for the model run a built
`meds_main` (`--meds-main`, default `../../build-ifx/meds_main`). The forcing build takes about
ten seconds, the census file about a minute, and each five-year run about 6 minutes on one core, or
1.5 minutes with `[run].n_threads = 8` in an OpenMP build (`docs/building.md`). It runs the five
years twice, with the default parameters and with the calibrated ones in `calibration/`.

## What the declarations are, and how each was checked

| Declared in `bci_site.toml` | Value | How it was established |
|---|---|---|
| clock | UTC−5, `stamp = "begin"` | V2: the shortwave envelope fits the model's sun best 5 min from this declaration (RMSE 36 W m⁻²; 61 if end-stamped). 0.006 % of the shortwave falls where the model sees night |
| VPD curve | Alduchov–Eskridge | V3: the provider's `vpd` is (1 − RH)·e_s(T) under it to 0.000 Pa. Bolton, the model's curve, misses by up to 5.2 Pa. The file stores RH, so this curve never reaches the model |
| heights | T/RH, wind and barometer at 41 m above the ground | the eddy-covariance height is 41 m. The mean pressure, 98.83 kPa, puts the barometer 180–210 m above sea level, not at the 150 m ground |
| location | 9.1568 N, −79.8486 E, 150 m | AmeriFlux PA-Bar |
| longwave | `LWdown` is the column `Rl_up` | the file labels its two longwave columns the wrong way round (the note under "The data") |

The forcing file is on a UTC clock and begin-stamped. It records the heights (`tq_height_m`,
`wind_height_m`, `height_above = "ground"`), and MEDS stops at open if `[forcing]` says otherwise.
Temperature, relative humidity, pressure, longwave and wind are values at the stamps, re-centred
from the tower's half-hour means, because MEDS interpolates states as instants. Rain and
shortwave stay means over each half hour.

## How each value was obtained

![Fill flags of the BCI forcing, per variable, 2012-2017](forcing_qc.png)

Temperature, humidity, pressure, rain and shortwave read as fully observed: the provider filled
their few gaps upstream, and the file does not say where. Two variables needed filling:
- **Wind:** 15 % is missing, including 212 negative values that the bounds screen removed. Gaps of
  up to 2 h are interpolated; longer ones take the mean diurnal variation.
- **Longwave:** 61 % is missing, because the radiometer arrived in 2015.

## Longwave: the fill, scored on observations it never saw

![The longwave fill scored on hidden observations](lw_comparison.png)

[`compare_longwave_fill.py`](../../scripts/prepare_flux_tower/compare_longwave_fill.py) hides 20 %
of the observed longwave in 10-day blocks. It fills the hidden records as the build does and scores
the fill, beside two references, against them (7,200 half hours):

| fill | bias | RMSE | r | diurnal-cycle RMSE |
|---|---|---|---|---|
| the model's synthesis, regressed onto the tower (the fill) | +2.3 | 14.5 | 0.79 | 4.3 |
| monthly day/night climatology | +4.4 | 18.6 | 0.63 | 7.1 |
| MEDS's `lwdown_source = "synthesize"` as it is | +20.9 | 29.7 | 0.61 | 22.1 |

(W m⁻²; the observed mean is 429.) Over other draws of the hidden blocks the fill scores RMSE
13.7–14.5. Filling from ERA5-Land or another source instead is left to the user: fill the tower
file before the build.

**The synthesis is regressed in two parts.** MEDS synthesizes longwave as εσT⁴[1 + a(1 − kt)],
with a = 0.22. Regressing on its clear-sky part εσT⁴ and its cloud part εσT⁴(1 − kt) separately lets
the site set the cloud coefficient: the pooled fit gives 0.10, and the two parts beat a single
regression on the whole synthesis on every draw (14.5 against 16.7 on the draw above). Used as it
is, the synthesis runs 21 W m⁻² above BCI's longwave, so a run here with no longwave at all would
receive that much too much.

## The model runs

[`meds_config_eval.toml`](meds_config_eval.toml) runs the five tower years, 2012-08-01 to
2017-08-01 UTC, with hourly output and CO₂ from the CMIP7 series. There is no spin-up.
- **The stand is the 2010 census.** [`make_census.py`](../../scripts/prepare_census/make_census.py)
  turns the 207,259 live trees of at least 1 cm into one patch per 20 m quadrat, 1,250 patches, and
  one row per distinct (quadrat, diameter), 84,937 rows. MEDS reads them and restructures the stand
  with its own cohort and patch fusion before the first step, down to 25 patches and 419 cohorts.
- **The soil starts at the site:** 298.65 K, the tower's mean air temperature, and 0.30 m³ m⁻³.
- **Soil carbon starts in steady state** with the census stand's litter: leaf and fine-root
  turnover from the PFT, and wood from the plot's biomass mortality, 1.90 % per year from 2005 to
  2010.
- **The PFT.** [`pft_parameters.toml`](pft_parameters.toml) is one evergreen broadleaf PFT for every
  tree, with MEDS's default allometry. Under it the census stand has LAI 5.6 and AGB 16.1 kgC m⁻²,
  against the census's own 15.1. The calibration below changes its fast parameters, not its
  allometry.
- **The leaf traits follow the canopy's light gradient** (`[trait_dynamics].trait_plasticity_on =
  true`). Each cohort's Vcmax25, Rd25, SLA and leaf lifespan are the PFT's top-of-canopy values
  scaled by the leaf area above it.
- **The canopy intercepts rain and dew** (`[fast].canopy_water_on = true`): leaves and wood hold a
  film of the water they catch, and it evaporates from them. The tower never measures a wet canopy,
  because its flag removes every raining half hour, so the calibration leaves the film capacities
  at the PFT's values.

`[forcing]` declares `tq_height = wind_height = 41`, `height_above = "ground"` and
`wind_exposure = "local"`. Each patch's forcing is moved from 41 m to its own canopy-air top.
- **Temperature** follows the dry adiabat: a 30 m canopy-air top over a 25 m stand is 0.11 K warmer
  than the tower.
- **Wind** follows the patch's log profile: that same patch gets 0.72 of the tower's wind.

The daily output carries each patch's `cas_depth_patch`, `air_temp_cas_top_patch` and
`wind_cas_top_patch`.

### Against the tower

![MEDS against the BCI tower: mean diurnal and seasonal cycles of carbon, water and energy](evaluation.png)

The default run takes 1.5 minutes with eight threads and 0.8 GB. MEDS reads the
census as 1,250 patches and 84,937 cohorts, with the stand's LAI 5.60 and AGB 16.12 kgC m⁻², exactly
as the census file states, and fuses it to 25 patches and 419 cohorts before the first step. The
count stays above `max_patch = 12` because `patch_light_tol_max` keeps dissimilar patches apart,
and the run says so at the end. The energy and water budgets close to machine precision.

**The stand over the five years:** LAI ends at 5.4 and AGB rises from 16.1 to 17.8 kgC m⁻², and
the patches fuse down to 15. Before the traits followed the light, every leaf had the canopy top's
traits, and LAI fell to 4.8.

**Against the tower**, over its measured hours (the site TOML's rule: FLAG = 1 for the turbulent
fluxes), 2012-08 to 2017-07, with the leaf's light use at #351's values (`phi_psii` 0.74, θ_J 0.7):

| | tower mean | MEDS mean | bias | r, hourly | r, mean seasonal cycle |
|---|---|---|---|---|---|
| GPP [µmol m⁻² s⁻¹] | 7.46 | 10.21 | +2.74 | 0.94 | 0.49 |
| NEE [µmol m⁻² s⁻¹] | −4.24 | −4.02 | +0.22 | 0.91 | 0.59 |
| latent heat [W m⁻²] | 75.5 | 71.5 | −4.0 | 0.92 | 0.74 |
| sensible heat [W m⁻²] | 32.4 | 65.0 | +32.5 | 0.94 | 0.96 |
| net radiation [W m⁻²] | 136.3 | 121.7 | −14.6 | 1.00 | 0.98 |

- **The diurnal cycles** are closely followed in shape (r ≥ 0.98 for every flux) and differ in size.
  Midday (11–14 h) GPP is 28.3 against the tower's 21.8 µmol m⁻² s⁻¹, and latent heat 207 against
  243 W m⁻².
- **Net radiation falls short by day because the canopy reflects too much:** the model's albedo,
  at incoming shortwave above 200 W m⁻², is 0.18 against the tower's 0.13.
- **The model puts too much of the day's energy into sensible heat.** At midday sensible heat is 238
  W m⁻² against the tower's 173: a Bowen ratio of 1.15 against 0.71. The tower's own fluxes close
  only 0.75 of its net radiation over whole days, and the gap behaves like missing sensible heat
  (the calibration's closure model), so part of the difference in H is the tower's. The model's
  friction velocity is about twice the tower's: 0.86 against 0.41 m s⁻¹ at night, and 1.14 against
  0.71 at midday.

**Before the longwave columns were corrected** the forcing carried the canopy's emission, 37 W m⁻²
more than the sky's. That hid the albedo error: net radiation looked right by day (bias +16.7
overall), and night-time net radiation was +6 against the tower's −33.

**Before v0.3.1 the reported sensible heat was about 105 W m⁻² too low.** It measured the canopy air
against the reference air's potential temperature referenced to the ground, 0.37 K above its actual
temperature at a 38 m canopy-air top, while the model exchanged heat on the actual temperatures. The
example then showed a sensible heat of −39 W m⁻², negative at every hour of the night, and the
reported fluxes left 105 W m⁻² of the net radiation unaccounted for. The same offset made every
stability solve 0.37 K more stable than the air was. Removing it raised the night-time friction
velocity from 0.78 to 0.87 m s⁻¹ and moved GPP, NEE, latent heat and net radiation by less than 1%.

## Calibrating the fast parameters

![The BCI calibration: each fitted key's prior z, and the validation scores](calibration.png)

[`calibration.toml`](calibration.toml) fits MEDS's sub-daily parameters to this tower with
[`scripts/calibrate_fast`](../../scripts/calibrate_fast/README.md), by the cross-site protocol of
[`MEDS_FAST_CALIBRATION_BEST_PRACTICE.md`](../../docs/dev_plans/MEDS_FAST_CALIBRATION_BEST_PRACTICE.md).
Every rule works at any tower; this tower's numbers come from its own data and forcing, and
`calibration.toml` holds only the choices that differ from the tool's defaults. The shipped set in
[`calibration/`](calibration) is the fit of 2026-10-04, with the canopy's interception on as in the
runs above.

### What the rules chose at BCI

- **The data.** The tower is read through [`bci_site.toml`](bci_site.toml), on its own half hours.
  The turbulent fluxes count where FLAG = 1, and every target only where the forcing's longwave and
  wind were observed.
- **u\* per target.** LE/Rnet is flat in u\*, and H/Rnet still rises at the top classes, so neither
  is filtered. GPP takes the provider's threshold, 0.4. Its own all-day diagnostic finds a plateau at
  0.33 (bootstrap 0.15–0.35): the declared alternative.
- **The closure.** The tower's daily H + LE is 0.75 of its net radiation, so f = 1.33. Within every
  VPD class H/Rnet rises 44–64 % from the calmest to the most turbulent third of the records, and
  LE/Rnet does not. The closure model therefore gives the whole gap to H, and LE stays as measured.
- **σ from the paired days:** LE 10.5 + 0.29 |LE|, H 8.9 + 0.14 |H|, and GPP's (NEE's) 2.2 +
  0.17 |NEE|, at the smoothed observation.
- **The windows.** Observed longwave confines every window to 2015–2017, the years the tower measured
  it. There are eight ten-day calibration windows, one per 1.5-month slot, and eight validation
  windows in the same slots of another of those years.
- **The seasonal runs** are the 2016 and 2017 dry seasons, 120 days each, ending at those years'
  deepest water deficits: 1,331 mm (the El Niño drought) and 698 mm. They are scored on LE only.
- **Fixed from coverage.** The flag removes every rain half hour, so the wet canopy is never sampled
  and the film capacities stay fixed.
- **The priors from the climate (EEO):**
  - `stomatal_g1`: the least-cost 2.80 kPa^0.5 (Lin et al.'s 3.77 for tropical rainforest is 0.9 sd
    away);
  - `vcmax25`: 41 µmol m⁻² s⁻¹, from the coordination of MEDS's own Rubisco- and light-limited rates;
  - Jmax/Vcmax 1.70, `ds_vcmax` 641.1 and `ds_jmax` 640.6, fixed at Kattge & Knorr for the growth
    temperature, 25.5 °C.
- **κ, the tower's respiration error.** The prior is 0.65 ± 0.10: the provider's respiration is
  about the soil chambers' alone, and its note says it "appeared underestimated".

### The fit

**Two choices by the owner after the first fit** (best-practice plan §0):
- **The albedo is off, and the NIR reflectance is fixed** at the PFT file's 0.45. With the albedo on,
  the fit took the leaf's NIR reflectance to 0.32 (prior z −3.1): reflectance plus transmittance
  0.57, against about 0.70 for tropical leaves. At realistic leaf optics the model's canopy reflects
  too much NIR, which points at its structure (clumping, wood optics) rather than its leaves. Both
  come back when that is revisited.
- **`stomata_psi_onset` is fixed** at half the PFT's turgor-loss point, −0.86 MPa. The first fit
  took it to its bound (z +2.6). The ten-day windows' longwave and LE set it, not the drought runs,
  and it traded against `stomatal_g1` (correlation 0.91). `wstress_sref_stomata` carries the
  drought response.

Levenberg–Marquardt on the seven remaining keys at once, from their priors' centres:
- the leaf angle;
- `vcmax25`, `stomatal_g1` and `g0`;
- the roughness length;
- `wstress_sref_stomata`;
- κ.

It ran 3 iterations, then one refresh of the chains, weights and σ, then 2 more. The cost went from
10,608 to 4,647. The screening's triage fixed no key. Without the albedo the leaf angle is informed
chiefly by GPP, then by LE and H, and less tightly (posterior/prior σ 0.43, against 0.17 with the
albedo). The fit took 616 trials (none failed), 16 minutes on three 40-core nodes. The revision's
staged fit took 4,804 trials and 31.8 core-hours.

| key | calibrated | 68 % | prior centre | prior z | posterior / prior σ |
|---|---|---|---|---|---|
| `vcmax25` [µmol m⁻² s⁻¹] | 30.9 | 30.0–31.9 | 41.0 (EEO) | −0.57 | 0.06 |
| `stomatal_g1` [kPa^0.5] | 2.99 | 2.88–3.11 | 2.80 (EEO) | +0.13 | 0.08 |
| `stomatal_g0` [mol m⁻² s⁻¹] | 0.0032 | 0.0024–0.0041 | 0.01 | −1.07 | 0.24 |
| `leaf_angle_mean` [°] | 54.0 | 49.9–57.9 | 45 | +0.92 | 0.43 |
| `z0m_ratio` (effective) | 0.052 | 0.050–0.054 | 0.13 | −2.71 | 0.11 |
| `wstress_sref_stomata` (effective) [MPa⁻¹] | 1.23 | 0.97–1.55 | 2.0 | −0.66 | 0.31 |
| κ (observation) | 0.738 | 0.701–0.774 | 0.65 | +0.87 | 0.36 |
| `leaf_reflect_nir` | 0.45, fixed | | | | |
| `stomata_psi_onset` [MPa] | −0.86, fixed | | | | |

- **κ = 0.74 makes the tower's respiration 4.44 µmol m⁻² s⁻¹ instead of 3.28**, about the soil
  chambers' ~4.3 alone. Its GPP becomes 3.30 kgC m⁻² yr⁻¹ instead of 2.86 over the measured records.
- **`vcmax25` lands at 30.9, inside its prior** (z −0.6), not at its floor as in every fit before
  this protocol. θ_J 0.7, `phi_psii` 0.74 and κ removed the level error the earlier fits pushed onto
  it.
- **The report flags no key.** No trait key is more than 2 prior sd from its evidence, and none is
  near a bound. `stomatal_g0` comes closest (z −1.07). The roughness length is 2.7 sd below ED2's
  0.13, but it is an effective key: it stands in for the canopy's structure.
- **The covariance is local only:** along the direction that is mostly `wstress_sref_stomata` the
  objective curves 0.45× as much as the quadratic predicts. `vcmax25` is correlated with
  `stomatal_g1` (−0.63) and κ (−0.66), and `stomatal_g1` with `g0` (−0.76).
- **The MAP is placed to about the width of its intervals.** The fit stops when an iteration lowers
  the cost by less than 0.1 %, about 5 here, while one posterior sd moves it by 1. The step a further
  iteration would take is 0.5–1.1 posterior sd on every key but `g0`, and 2.2 on `g0`.
- **The declared alternatives:**
  - **GPP u\* at the plateau's 0.33:** every key moves under 0.6 posterior sd.
  - **Bowen closure instead of the attribution's:** `wstress_sref_stomata` and `stomatal_g1` move
    about 16 sd in the linear estimate, so the fit was rerun under Bowen. There `vcmax25` is 27.9,
    `stomatal_g1` 4.51 and `wstress_sref_stomata` 0.11. The closure choice is still this site's
    largest uncertainty.
  - **Partitioning:** not applicable; the provider gives one partitioning.

**On the validation windows** (RMSE in units of each target's σ, never fitted), the validation
passes: the calibrated cost is lower, and no target is more than 10 % worse.

| target | default | calibrated |
|---|---|---|
| u\* | 1.74 | **0.67** |
| GPP (with κ) | 2.06 | **1.45** |
| upwelling longwave | 3.82 | **3.55** |
| LE (as measured) | 1.65 | 1.57 |
| H (closure-corrected) | 4.16 | 4.30 (+3 %) |
| objective | 61,507 | **51,650** |

The albedo is no longer scored (the first fit took its RMSE from 5.34 to 1.36 σ through the NIR
reflectance). LE improves a little, and H is a little worse. At the end of the fit H's χ² per row is
1.6, even with its σ tripled (the cap). By day the model's sensible heat is 0.77 of the
closure-corrected tower's, against 0.86 in the first fit: with the NIR reflectance back at 0.45 the
canopy absorbs less. The fit's report is [`calibration/report.md`](calibration/report.md).

### The calibrated run over the five years

`run_example.py` runs the five years again with
[`calibration/meds_config_calibrated.toml`](calibration/meds_config_calibrated.toml), which reads
[`calibration/pft_parameters_calibrated.toml`](calibration/pft_parameters_calibrated.toml) and writes
`output/cal-*`. Both files are the example's configs, comments kept, with the fitted values and the
Kattge & Knorr shape keys written in. Their header labels the effective keys. κ is not a model key,
so it is not in them. `evaluation.png` draws the calibrated run beside the default.

| over the tower's measured hours | tower | default | calibrated |
|---|---|---|---|
| GPP [µmol m⁻² s⁻¹] | 7.46 (8.61 with κ) | 10.21 | **8.54** |
| NEE [µmol m⁻² s⁻¹] | −4.24 | −4.02 | −3.31 |
| latent heat [W m⁻²] | 75.5 | 71.5 | 64.1 |
| sensible heat [W m⁻²] | 32.4 (as measured) | 65.0 | 70.4 |
| net radiation [W m⁻²] | 136.3 | 121.7 | 122.2 |
| albedo, at incoming shortwave above 200 W m⁻² | 0.129 | 0.182 | 0.173 |
| u\* at night / at midday [m s⁻¹] | 0.41 / 0.71 | 0.86 / 1.14 | **0.49 / 0.70** |
| midday GPP, LE, H | 21.8, 243, 173 | 28.3, 207, 238 | 23.4, 188, 243 |
| stand at the end: LAI, AGB [kgC m⁻²] | | 5.38, 17.8 | 5.50, 17.2 |

- **The full-record check passes.** The budgets close: whole-site energy and water have no failed
  check in 3.0 million. The dry season (January to April) is better than the default's: GPP RMSE
  3.92 against 6.33 µmol m⁻² s⁻¹, and LE RMSE 35.2 against 38.6 W m⁻².
- **GPP matches the tower corrected for κ** (8.54 against 8.61, over the same hours), and u\* is
  close to the tower's.
- **The albedo is 0.173 against the tower's 0.129.** The steeper leaves (54.0°) bring it down a little
  from the default's 0.182, and the rest is the canopy-structure error the albedo target is off for.
  Net radiation stays 14 W m⁻² short.
- **April's GPP in the drought year** is 7.2 against the tower's 6.3 in 2016 (the default's 6.5; the
  shipped set of v0.3 gave 4.1).
- **LE is lower than the default's** (64.1 against 71.5) and **H higher** (70.4 against 65.0). Over
  the record the closure-corrected tower's H is about 68 W m⁻². The model's H is above it at night
  (about 0 against the tower's −24 W m⁻²) and below it by day (0.77 of it in the fit's windows).
  That is the structural H misfit above.

### Redoing it

```bash
python run_example.py --calibrate --workers 40
```

This runs the default five years, runs the fit, and writes `calibration/` anew. On one 40-core node
the fit takes about an hour. Run the tool directly for a cluster: see
[`scripts/calibrate_fast/README.md`](../../scripts/calibrate_fast/README.md).

## Files

| File | What it is |
|---|---|
| [`bci_site.toml`](bci_site.toml) | the declaration of the tower data |
| [`fetch_bci_data.py`](fetch_bci_data.py) | downloads and verifies the data |
| [`run_example.py`](run_example.py) | runs every step |
| [`plot_forcing.py`](plot_forcing.py), [`plot_evaluation.py`](plot_evaluation.py), [`plot_calibration.py`](plot_calibration.py) | the figures |
| [`calibration.toml`](calibration.toml) | the calibration's choices: κ's prior and the fit's settings (the rules choose the rest) |
| [`calibration/`](calibration) | the fit's results and report (`fit.json`, `report.md`), the calibrated configs |
| [`bci_census.toml`](bci_census.toml) | the declaration of the census the run starts from |
| [`meds_config_eval.toml`](meds_config_eval.toml) | the MEDS run |
| [`output_variables.toml`](output_variables.toml), [`pft_parameters.toml`](pft_parameters.toml) | its output list and PFT |
