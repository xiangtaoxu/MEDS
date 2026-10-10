# Example 04 — Column biophysics at a flux tower (Barro Colorado Island)

The coupled column of MEDS — canopy radiation, leaf gas exchange, the energy balances of leaves,
canopy air and soil, plant hydraulics, and the soil's water and heat — solved every 15 minutes in
each patch of a forest. The site is the Barro Colorado Island (BCI) flux tower in Panama (AmeriFlux
PA-Bar, 41 m above a seasonal tropical forest), over its five years, 2012–2017. How the column is
solved: [`docs/science/column_biophysics.md`](../../docs/science/column_biophysics.md).

The example:

1. builds a MEDS forcing file from the tower's own meteorology;
2. starts the forest from the 2010 census of the BCI 50-ha plot, with no spin-up;
3. calibrates the fast (sub-daily) parameters against the same tower;
4. compares the default and the calibrated runs with the tower's fluxes;
5. reruns ten days at half-hourly output to show what the column solves.

## Run

```bash
cd examples/example04_column_biophysics
python run_example.py --forcing-only            # the data, the forcing and its figures
python run_example.py                           # every step, with the shipped calibration
python run_example.py --calibrate --workers 40  # ... and redo the calibration first
```

It needs numpy, pandas, netCDF4, matplotlib and a built `meds_main` (`--meds-main`, default
`../../build-ifx/meds_main`). Every step together takes about 14 minutes on 8 threads. The
calibration takes about 5 hours on one 40-core node, less across a cluster with
[`scripts/calibrate_fast`](../../scripts/calibrate_fast/README.md).

## The data

- **Tower:** M. Detto, "Barro Colorado Island – eddy covariance flux data (2012–2017)",
  [Zenodo 6456527](https://zenodo.org/records/6456527), CC0. [`fetch_bci_data.py`](fetch_bci_data.py)
  downloads it into `data/` and checks its md5. Publications should acknowledge CTFS-ForestGEO,
  which supported the tower.
- **Census:** the BCI 50-ha plot's 2010 census (and 2005, for the plot's mortality). Condit et al.
  2019, Dryad [doi:10.15146/5xcp-0d46](https://doi.org/10.15146/5xcp-0d46), CC0; the PIs ask to be
  told of papers that use it. Dryad refuses scripts, so download `bci.tree.zip` in a browser;
  [`bci_census.toml`](bci_census.toml) names the copy to read.

## 1. The forcing

[`bci_site.toml`](bci_site.toml) declares what the tower file holds: its clock (UTC−5), its
measurement height (41 m), its location, and that the file swaps its two longwave labels.
[`make_tower_forcing.py`](../../scripts/prepare_flux_tower/make_tower_forcing.py) checks each
declaration against the sun and the data, and stops on a disagreement. MEDS moves each sample from
41 m to the top of each patch's canopy air.

![Fill flags of the BCI forcing, per variable, 2012-2017](forcing_qc.png)

Two variables have gaps: **wind** (15 %), filled from the neighbouring hours or the mean diurnal
cycle, and **longwave** (61 %: the radiometer arrived in 2015), filled by MEDS's longwave synthesis
corrected to the tower. On observations hidden from it, the longwave fill scores an RMSE of
14.5 W m⁻², against 18.6 for a monthly climatology.

![The longwave fill scored on hidden observations](lw_comparison.png)

## 2. The model

[`meds_config_eval.toml`](meds_config_eval.toml) runs 2012-08-01 to 2017-08-01.
- **The stand is the 2010 census:** 207,259 live trees of at least 1 cm, which MEDS groups into 130
  patches and 1,924 cohorts. One evergreen PFT ([`pft_parameters.toml`](pft_parameters.toml)) with
  MEDS's default allometry puts it at LAI 5.6 and 16.1 kgC m⁻² aboveground (the census's own
  estimate is 15.1).
- **Leaf traits follow the light** down the canopy, and the canopy holds rain and dew.
- **Each tree roots to a depth set by its height** (5 m for a 35 m tree), and water rises to its
  leaves through its sapwood, whose conductivity is set by wood density.
- **The soil is BCI's clay**, 6 m deep, its water retention and conductivity fitted to measurements
  on the island (Kupers et al. 2019; Godsey et al. 2004).

## 3. The calibration

![The BCI calibration: each fitted key's prior z, and the validation scores](calibration.png)

[`calibration.toml`](calibration.toml) fits the fast parameters with
[`scripts/calibrate_fast`](../../scripts/calibrate_fast/README.md), following
[a cross-site protocol](../../docs/dev_plans/MEDS_FAST_CALIBRATION_BEST_PRACTICE.md). The fit sees
the tower's latent heat (LE), sensible heat (H), GPP, upwelling longwave and u\* in eight 10-day
windows spread over the year, and LE over the 2016 and 2017 dry seasons; eight other windows
validate it. Two choices are BCI's own:
- the tower's energy balance closes to 0.75 of its net radiation, and the gap behaves like missing
  H, so the fit gives the whole gap to H;
- the tower's respiration is about the soil's alone, so its GPP is too low; the fit scales it by a
  factor κ.

The calibrated configs, the fit and its report are in [`calibration/`](calibration).

| key | calibrated (68 % range) | prior centre |
|---|---|---|
| `vcmax25` [µmol m⁻² s⁻¹] | 31.3 (30.4–32.2) | 41.0 (from the climate) |
| `stomatal_g1` [kPa^0.5] | 2.82 (2.74–2.90) | 2.80 (from the climate) |
| `stomatal_g0` [mol m⁻² s⁻¹] | 0.0005 (0.0003–0.0011) | 0.01 |
| `leaf_angle_mean` [°] | 55.1 (51.2–58.7) | 45 |
| `z0m_ratio` | 0.054 (0.052–0.056) | 0.13 |
| κ | 0.767 (0.728–0.804) | 0.65 |

On the validation windows (RMSE in units of each target's error):

| target | default | calibrated |
|---|---|---|
| u\* | 1.72 | 0.67 |
| GPP | 2.50 | 1.58 |
| upwelling longwave | 3.87 | 3.54 |
| LE | 1.70 | 1.59 |
| H | 4.34 | 4.35 |

The energy-balance gap is the largest uncertainty: sharing it between H and LE instead moves
`stomatal_g1` to 4.44 and `vcmax25` to 28.2.

## 4. Against the tower

![MEDS against the BCI tower: mean diurnal and seasonal cycles of carbon, water and energy](evaluation.png)

The dashed lines are the calibration's targets: GPP with κ, and H with the energy-balance gap added.
Means over the tower's measured hours, 2012-08 to 2017-07:

| | tower | target | default | calibrated |
|---|---|---|---|---|
| GPP [µmol m⁻² s⁻¹] | 7.46 | 8.04 | 10.54 | 8.54 |
| LE [W m⁻²] | 75.5 | 75.5 | 79.8 | 63.9 |
| H by day [W m⁻²] | 91.3 | 166.6 | 126.6 | 143.5 |
| NEE [µmol m⁻² s⁻¹] | −4.24 | | −3.58 | −2.68 |
| net radiation [W m⁻²] | 136.3 | | 122.0 | 122.6 |

- **Carbon:** the calibration brings GPP from 31 % above its target to 6 % above.
- **Energy:** the model's canopy reflects too much (albedo 0.17–0.18 against the tower's 0.13), so
  its net radiation is 14 W m⁻² short and its daytime H stays below the target.
- **Water:** the calibration lowers LE to 63.9 W m⁻², below the tower's 75.5, although it improves LE
  on the validation windows.
- Whole-site energy and water close in both runs.

## 5. Ten days inside the column

![Ten days in one patch at half-hourly resolution: forcing, temperatures, canopy-air CO2, leaf water potential, soil moisture and temperature](window.png)

The example runs the calibrated model to 20 April 2016, then restarts it for ten days with every
patch's and cohort's states written each half hour ([`output_window.toml`](output_window.toml)).
These are the last days of the 2016 El Niño dry season: four bright dry days, small showers, then
42 mm of rain on the 27th. The figure shows the patch that covers the most ground (8 % of the site,
LAI 5.9) and its tallest tree.
- **Temperatures.** Only the air at 41 m is forcing. The canopy air peaks 1–2 K above it, and the
  top tree's leaves 3–6 K above it at midday. The soil at 5 cm swings 1–3 K a day, peaking hours
  after the air.
- **Canopy-air CO₂** falls 4–10 µmol mol⁻¹ by day as the canopy draws it down, and builds up again
  at night, most on the calmest night.
- **Leaf water potential** refills to −0.41 MPa every night and falls to about −1.1 MPa at midday on
  the bright days: the tension that pulls the day's transpiration up through the sapwood. Panama's
  canopy trees measure −0.6 to −1.0 MPa before dawn and −1.4 to −2.0 at midday in a dry season, so
  the model's tree has more water than the measured ones: its roots reach soil that stays moist at
  depth.
- **Soil moisture.** The top 15 cm stays near 0.37 m³ m⁻³, wetter than the tower's sensor (0.28),
  because the trees draw on the whole 6 m column rather than on its top. The storm wets the top
  layers within hours.

## Files

| File | What it is |
|---|---|
| [`run_example.py`](run_example.py) | runs every step |
| [`fetch_bci_data.py`](fetch_bci_data.py) | downloads and checks the tower data |
| [`bci_site.toml`](bci_site.toml), [`bci_census.toml`](bci_census.toml) | the declarations of the tower data and the census |
| [`meds_config_eval.toml`](meds_config_eval.toml), [`pft_parameters.toml`](pft_parameters.toml), [`output_variables.toml`](output_variables.toml) | the MEDS run, its PFT and its output |
| [`calibration.toml`](calibration.toml) | the calibration's BCI choices |
| [`calibration/`](calibration) | the calibrated configs, the fit (`fit.json`) and its report |
| [`output_window.toml`](output_window.toml) | the ten-day window's output |
| [`plot_forcing.py`](plot_forcing.py), [`plot_evaluation.py`](plot_evaluation.py), [`plot_calibration.py`](plot_calibration.py), [`plot_window.py`](plot_window.py) | the figures |

Design records:
[forcing](../../docs/dev_plans/MEDS_FLUX_TOWER_FORCING_PLAN.md),
[census start](../../docs/dev_plans/MEDS_BCI_CENSUS_INIT_PLAN.md),
[calibration](../../docs/dev_plans/MEDS_FAST_CALIBRATION_BEST_PRACTICE.md).
