# Example 03 — Demography

The demography module of MEDS on its own: cohorts, patches, their fusion and fission, recruitment
and treefall disturbance, with no carbon, water or energy. **This example shows how the module runs
on vital rates handed to it from outside; it is not a calibration.** Here those rates are three
simple laws fitted to the censuses of the Barro Colorado Island (BCI) 50-ha plot, and a Python driver
passes them to the engine's law-free apply-primitives through
[`meds.demography`](../../python/meds/demography/_site.py) (`Site.apply_rates`).

## The laws

| Law | Rate | Form | Fitted by |
|---|---|---|---|
| growth | dbh growth g [cm/yr] per cohort | g = [g_min + (g_max − g_min) / (1 + (D/D₀)^−k)] · e^(−b·L) | Gamma quasi-likelihood |
| mortality | death rate m [1/yr] per cohort | m = γ + α · e^(−β·g) (Camac et al. 2018) | maximum likelihood, P(die in Δt) = 1 − e^(−m·Δt) |
| recruitment | new stems ≥ 1 cm [1/m²/yr] per PFT and patch | R = exp(c₀ + c₁·LAI + c₂·LAI_PFT) | Poisson GLM, weighted by area × Δt |

| PFT | wood density [g cm⁻³] | g_min [cm/yr] | g_max [cm/yr] | D₀ [cm] | k | b | γ [1/yr] | α [1/yr] | β [yr/cm] |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.36 | 0.36 | 1.03 | 17 | 1.3 | 0.21 | 0.023 | 0.23 | 11 |
| 2 | 0.50 | 0.17 | 0.59 | 8 | 2.7 | 0.16 | 0.019 | 0.08 | 14 |
| 3 | 0.68 | 0.08 | 0.50 | 23 | 1.7 | 0.07 | 0.019 | 0.11 | 38 |

- **Growth** rises with dbh D from g_min to g_max, half-way at D₀, and shade shrinks it by e^(−b·L).
  L is the overtopping LAI, the leaf area of taller trees per m² of ground: the competition index the
  engine computes for every cohort (`overtopping_lai`), which the driver reads as it is. For each
  census tree it is the leaf area of taller trees within 20 m, with height and leaf area from MEDS's
  pan-tropical allometry; trees of equal height share a layer, as in the engine. Zero and negative
  increments (the tape's millimetre, or a broken stem) are set to a tenth of the smallest positive
  one, 0.0017 cm/yr, which raises the census's mean growth from 0.085 to 0.094 cm/yr.
- **Mortality** depends on growth alone, as in Camac et al. (2018), who fitted the same form to the BCI
  censuses. Its growth is the growth law's, not the tree's measured growth: a cohort's growth is the
  law's too.
- **Recruitment** depends on the leaf area index of the patch, all of it (LAI) and the PFT's own
  (LAI_PFT) -- in the census, the leaf area within 20 m of a quadrat's centre. Trees and quadrats
  within 20 m of the plot's edge are not fitted.
- **PFTs are wood-density classes**: species' wood density from the Global Wood Density Database
  ([`bci_wood_density.csv`](bci_wood_density.csv): neotropical entries, species → genus → family),
  split at 0.42 and 0.58 g cm⁻³ so each class holds a third of the 1985 basal area.
- **A tree is one stem**, as in MEDS: a tree that loses its main stem has died, and returns as a
  recruit when a stem is measured again. Recruits enter at the census's 1 cm (`min_cohort_height`
  3.13 m); trees stop growing taller at the BCI canopy's 42 m (`hgt_max`, reached at dbh 99 cm) but
  keep thickening.
- **Treefall** (`patch_disturbance_rate` 0.014 yr⁻¹) kills canopy trees above 10 m on the disturbed
  area. The census counts those deaths too, so the driver takes the rate off the mortality of
  cohorts tall enough to die in a gap.

The coefficients are in [`vital_rates.json`](vital_rates.json); [`census_laws.py`](census_laws.py)
evaluates the three formulas each step.

## Data

The BCI 50-ha plot censuses 2–7 (1985–2010), Condit et al. 2019, Dryad
[doi:10.15146/5xcp-0d46](https://doi.org/10.15146/5xcp-0d46), CC0; the PIs ask to be told of papers
that use it. Dryad refuses scripts, so download `bci.tree.zip` in a browser and export the tree
tables `bci.tree2` … `bci.tree7` as `bci_1985.csv` … `bci_2010.csv`. The 1982 census is not used: it
rounded small stems to 5 mm. Only processed summaries are committed: the coefficients and
[`census_stand.csv`](census_stand.csv) (stems and basal area by PFT and size class, per census).

## Run

Build the Python library once (as in [example 01](../example01_leaf_gas_exchange/README.md#run)),
then from the repository root:

```bash
E=examples/example03_demography
python $E/prepare_census.py /path/to/bci_census   # tables in $E/data/, the 1985 stand, census_stand.csv
python $E/fit_vital_rates.py                      # scikit-learn, scipy; ~2 min -> vital_rates.json
PYTHONPATH=python python $E/run_demography.py --start census   # 1985 to 2100 -> output/census.nc
PYTHONPATH=python python $E/run_demography.py --start bare     # 300 years   -> output/bare.nc
python $E/plot_demography.py                      # demography.png
python post_proc/plot_forest_structure.py    $E/output/bare.nc -o $E/canopy_profile.gif
python post_proc/plot_landscape_3d.py        $E/output/bare.nc -o $E/landscape_3d.png
python post_proc/animate_landscape_growth.py $E/output/bare.nc -o $E/landscape_growth_3d.gif
```

The bare-ground run needs only the committed coefficients; the census run also needs the 1985 stand
that `prepare_census.py` writes.

## Results

![The census-trained laws and the stands they make](demography.png)

**Top row:** the three laws. **Bottom row:** a run from the 1985 census beside the next five
censuses, the size distribution in 2010, and a run from near-bare ground.

| Cross-validated on 1-ha blocks | trees (quadrats) | class means |
|---|---|---|
| growth, R² | 0.09 | 0.86 |
| mortality, AUC / R² | 0.60 | 0.74 |
| recruitment, R² | 0.29 | plot total 110 stems ha⁻¹ yr⁻¹, as observed |

| | stems ha⁻¹ | ≥ 10 cm ha⁻¹ | basal area m² ha⁻¹ | by PFT 1 / 2 / 3 |
|---|---|---|---|---|
| census 1985 | 4841 | 414 | 31.1 | 10.6 / 9.8 / 10.7 |
| census 2010 | 4145 | 416 | 30.5 | 10.7 / 9.4 / 10.4 |
| model 2010, from the 1985 census | 4038 | 411 | 31.5 | 9.9 / 10.5 / 11.1 |
| model 2100, from the 1985 census | 3675 | 543 | 34.9 | 9.2 / 14.3 / 11.4 |
| model, 100 years from bare ground | 3738 | 534 | 32.4 | 14.3 / 14.7 / 3.3 |
| model, 300 years from bare ground | 3639 | 512 | 36.1 | 11.4 / 16.6 / 8.1 |

- **From the census, the module tracks the plot for 25 years**: the stem decline, the trees above
  10 cm, the basal area and the size distribution class by class.
- **Run on to 2100, the stand gains basal area slowly** (0.03 m² ha⁻¹ yr⁻¹), most of it in PFT 2.
  These are BCI's present rates run forward, not a prediction.
- **From bare ground the stand goes through a succession.** Under an open canopy the light-wooded
  PFT 1 grows fastest and rises with PFT 2, reaching the census's basal area within about 90 years;
  as the canopy closes PFT 1 declines, PFT 2 holds the canopy, and the dense-wooded PFT 3 builds up
  slowly beneath. Basal area ends 18 % above the census's.

<p align="center">
  <img src="canopy_profile.gif" height="230" alt="Canopy-layer stand profile with vertical LAI">
  &nbsp;&nbsp;
  <img src="landscape_growth_3d.gif" height="230" alt="3D landscape, trees growing in place">
</p>

The near-bare-ground run as the leaf-area profile beside a stand cross-section (left), and as a
synthetic landscape whose trees keep their place as they grow (right); green PFT 1, blue PFT 2,
magenta PFT 3. [`landscape_3d.png`](landscape_3d.png) is its last year.

## Files

| File | |
|---|---|
| [`prepare_census.py`](prepare_census.py), [`fit_vital_rates.py`](fit_vital_rates.py) | the census tables, and the three laws fitted to them |
| [`census_laws.py`](census_laws.py), [`run_demography.py`](run_demography.py) | the laws the driver applies, and the driver |
| [`plot_demography.py`](plot_demography.py) | the figure |
| [`example_config_main.toml`](example_config_main.toml), [`example_config_pft.toml`](example_config_pft.toml) | the demographic settings, treefall, the census file, the PFTs |
| `_cadence.py`, `_write_nc.py` | the monthly and yearly steps; the netCDF the renders read |
