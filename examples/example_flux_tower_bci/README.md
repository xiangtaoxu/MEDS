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

The design records are
[`docs/dev_plans/MEDS_FLUX_TOWER_FORCING_PLAN.md`](../../docs/dev_plans/MEDS_FLUX_TOWER_FORCING_PLAN.md)
and [`docs/dev_plans/MEDS_BCI_CENSUS_INIT_PLAN.md`](../../docs/dev_plans/MEDS_BCI_CENSUS_INIT_PLAN.md).

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
```

`run_example.py` needs numpy, pandas, netCDF4 and matplotlib, and for the model run a built
`meds_main` (`--meds-main`, default `../../build-ifx/meds_main`). The forcing build takes about
ten seconds, the census file about a minute, and the five-year run about 7 minutes on one core, or
6 minutes with `[run].n_threads = 4` (more threads are slower for now, #325).

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
  tree, with MEDS's default allometry. The example is about the forcing and the start, and this is
  not a calibration for BCI. Under it the census stand has LAI 5.6 and AGB 16.1 kgC m⁻², against the
  census's own 15.1.

`[forcing]` declares `tq_height = wind_height = 41`, `height_above = "ground"` and
`wind_exposure = "local"`. Each patch's forcing is moved from 41 m to its own canopy-air top.
- **Temperature** follows the dry adiabat: a 30 m canopy-air top over a 25 m stand is 0.11 K warmer
  than the tower.
- **Wind** follows the patch's log profile: that same patch gets 0.72 of the tower's wind.

The daily output carries each patch's `cas_depth_patch`, `air_temp_cas_top_patch` and
`wind_cas_top_patch`.

### Against the tower

![MEDS against the BCI tower: mean diurnal and seasonal cycles of carbon, water and energy](evaluation.png)

The run takes 5.7 minutes with four threads, about 7 on one core, and 0.8 GB. MEDS reads the
census as 1,250 patches and 84,937 cohorts, with the stand's LAI 5.60 and AGB 16.12 kgC m⁻², exactly
as the census file states, and fuses it to 25 patches and 419 cohorts before the first step. The
count stays above `max_patch = 12` because `patch_light_tol_max` keeps dissimilar patches apart,
and the run says so at the end. The energy and water budgets close to machine precision.

**The stand over the five years:** LAI falls from 5.6 to 4.8 and AGB rises from 16.1 to 17.5 kgC m⁻²,
and the patches fuse down to 15. The large trees grow and the canopy thins; small trees do not grow
under this PFT with the default allometry.

**Against the tower**, over its measured hours (FLAG = 1 for the turbulent fluxes), 2012-08 to 2017-07:

| | tower mean | MEDS mean | bias | r, hourly | r, mean seasonal cycle |
|---|---|---|---|---|---|
| GPP [µmol m⁻² s⁻¹] | 7.46 | 11.07 | +3.61 | 0.94 | 0.43 |
| NEE [µmol m⁻² s⁻¹] | −4.24 | −4.04 | +0.20 | 0.90 | 0.45 |
| latent heat [W m⁻²] | 75.5 | 56.5 | −19.0 | 0.93 | 0.69 |
| sensible heat [W m⁻²] | 32.4 | −39.1 | −71.5 | 0.86 | −0.59 |
| net radiation [W m⁻²] | 136.3 | 120.4 | −16.0 | 1.00 | 0.98 |

- **The diurnal cycles** are closely followed in shape (r ≥ 0.99 for every flux) and differ in size.
  Midday GPP is 30 against the tower's 22 µmol m⁻² s⁻¹, and midday latent heat 164 against 237 W m⁻².
- **Net radiation falls short by day because the canopy reflects too much.** Over the 695 days the
  tower measured all four components, the model's albedo is 0.26 against the tower's 0.13: it
  absorbs 148 W m⁻² of shortwave where the tower's canopy absorbs 173. It emits 10 W m⁻² less
  longwave, which offsets part of that. At night the model reads −22 W m⁻² against the tower's −33.
- **Sensible heat is the largest miss:** −96 W m⁻² at night against the tower's −23, and 90 at
  midday against 163. The census canopy is taller than the tower: the canopy-air tops reach 51 m in
  the tallest patches, above the tower's 41 m, so those patches take forcing moved up from below
  their own top. Whether that explains the flux is not yet checked.
- **The seasonal cycles** are weaker: the model's GPP is highest in the dry season, the tower's early
  in the wet season, and the model's dry-season sensible heat falls where the tower's rises.

**Before the longwave columns were corrected** the forcing carried the canopy's emission, 37 W m⁻²
more than the sky's. That hid the albedo error: net radiation looked right by day (bias +16.7
overall), night-time net radiation was +6 against the tower's −33, and sensible heat was biased
−46 W m⁻².

## Files

| File | What it is |
|---|---|
| [`bci_site.toml`](bci_site.toml) | the declaration of the tower data |
| [`fetch_bci_data.py`](fetch_bci_data.py) | downloads and verifies the data |
| [`run_example.py`](run_example.py) | runs every step |
| [`plot_forcing.py`](plot_forcing.py), [`plot_evaluation.py`](plot_evaluation.py) | the figures |
| [`bci_census.toml`](bci_census.toml) | the declaration of the census the run starts from |
| [`meds_config_eval.toml`](meds_config_eval.toml) | the MEDS run |
| [`output_variables.toml`](output_variables.toml), [`pft_parameters.toml`](pft_parameters.toml) | its output list and PFT |
