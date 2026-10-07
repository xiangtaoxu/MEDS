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
| growth | dbh growth g [cm/yr] per cohort | g = [g_min + (g_max − g_min) / (1 + (D/D₀)^−k)] · e^(−b·BAL) | Gamma quasi-likelihood |
| mortality | death rate m [1/yr] per cohort | m = γ + α · e^(−β·g) (Camac et al. 2018) | maximum likelihood, P(die in Δt) = 1 − e^(−m·Δt) |
| recruitment | new stems ≥ 1 cm [1/m²/yr] per PFT and patch | R = exp(c₀ + c₁·BA + c₂·BA_PFT) | Poisson GLM, weighted by area × Δt |

| PFT | wood density [g cm⁻³] | g_min [cm/yr] | g_max [cm/yr] | D₀ [cm] | k | b [ha m⁻²] | γ [1/yr] | α [1/yr] | β [yr/cm] |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.36 | 0.18 | 1.05 | 27 | 1.3 | 0.016 | 0.023 | 0.31 | 13 |
| 2 | 0.50 | 0.09 | 0.49 | 10 | 2.4 | 0.010 | 0.020 | 0.13 | 20 |
| 3 | 0.68 | 0.06 | 0.50 | 27 | 1.6 | 0.002 | 0.020 | 1.39 | 81 |

- **Growth** rises with dbh D from g_min to g_max, half-way at D₀, and shade shrinks it by e^(−b·BAL).
  BAL is the basal area of larger trees within 20 m [m² ha⁻¹]; a cohort counts the larger cohorts in
  its patch and half of its own size class, the average its trees have above them. Zero and negative
  increments (the tape's millimetre, or a broken stem) are set to a tenth of the smallest positive
  one, 0.0017 cm/yr, which raises the census's mean growth from 0.085 to 0.094 cm/yr.
- **Mortality** depends on growth alone, as in Camac et al. (2018), who fitted the same form to the BCI
  censuses. Its growth is the growth law's, not the tree's measured growth: a cohort's growth is the
  law's too.
- **Recruitment** depends on the basal area within 20 m of a quadrat's centre, all of it (BA) and the
  PFT's own (BA_PFT). Trees and quadrats within 20 m of the plot's edge are not fitted.
- **PFTs are wood-density classes**: species' wood density from the Global Wood Density Database
  ([`bci_wood_density.csv`](bci_wood_density.csv): neotropical entries, species → genus → family),
  split at 0.42 and 0.58 g cm⁻³ so each class holds a third of the 1985 basal area.
- **A tree is one stem**, as in MEDS: a tree that loses its main stem has died, and returns as a
  recruit when a stem is measured again. Recruits enter at the census's 1 cm (`min_cohort_height`
  3.13 m).
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
PYTHONPATH=python python $E/run_demography.py --start census --years 300                  # ~1 min
PYTHONPATH=python python $E/run_demography.py --start bare --years 300 --write-nc $E/output/bare_ground.nc
python $E/plot_demography.py                      # demography.png
python post_proc/plot_forest_structure.py    $E/output/bare_ground.nc -o $E/canopy_profile.gif
python post_proc/plot_landscape_3d.py        $E/output/bare_ground.nc -o $E/landscape_3d.png
python post_proc/animate_landscape_growth.py $E/output/bare_ground.nc -o $E/landscape_growth_3d.gif
```

The bare-ground run needs only the committed coefficients; the census run also needs the 1985 stand
that `prepare_census.py` writes.

## Results

![The census-trained laws and the stands they make](demography.png)

**Top row:** the three laws. **Bottom row:** a run from the 1985 census beside the next five
censuses, the size distribution in 2010, and a run from near-bare ground.

| Cross-validated on 1-ha blocks | trees (quadrats) | class means |
|---|---|---|
| growth, R² | 0.09 | 0.88 |
| mortality, AUC / R² | 0.60 | 0.78 |
| recruitment, R² | 0.29 | plot total 110 stems ha⁻¹ yr⁻¹, as observed |

| | stems ha⁻¹ | ≥ 10 cm ha⁻¹ | basal area m² ha⁻¹ | by PFT 1 / 2 / 3 |
|---|---|---|---|---|
| census 1985 | 4841 | 414 | 31.1 | 10.6 / 9.8 / 10.7 |
| census 2010 | 4145 | 416 | 30.5 | 10.7 / 9.4 / 10.4 |
| model 2010, from the 1985 census | 4071 | 406 | 30.9 | 9.4 / 10.3 / 11.2 |
| model 2285, from the 1985 census | 3676 | 480 | 33.3 | 6.8 / 14.0 / 12.5 |
| model, 300 years from bare ground | 3693 | 443 | 32.6 | 7.3 / 14.9 / 10.4 |

- **From the census, the module tracks the plot for 25 years**: the stem decline, the trees above
  10 cm, the basal area and the size distribution class by class.
- **Run on, the stand keeps its size and changes its make-up.** Basal area settles near 33 m² ha⁻¹
  as PFT 1 gives way to PFTs 2 and 3. The census shows the decline begun: PFT 1's stems fell by a
  third from 1985 to 2010. These are BCI's present rates run forward, not a prediction.
- **From bare ground** the stand reaches the census's basal area in about 200 years. PFT 2 leads,
  then PFT 3 takes over part of its place. There is no pioneer flush: on open ground PFT 1 recruits
  only 13 stems ha⁻¹ yr⁻¹, PFT 2 49 and PFT 3 69.

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
