# Demography example

A self-contained **demographic spin-up** of MEDS and the figures it produces. The configuration is a
pair of files: [`example_config_main.toml`](example_config_main.toml) (run/engine/IO settings, which
names the PFT file) and [`example_config_pft.toml`](example_config_pft.toml) (PFT traits + allometry +
mortality coefficients) — a **250-year** daily spin-up (2000-01-01 → 2250-01-01) from near-bare
ground, writing to `example_output/` with the prefix `example_output`. (For the standalone leaf-level
photosynthesis example, see [`../example_leaf_gas_exchange/`](../example_leaf_gas_exchange/).)

## Reproduce

Run from the **repository root** (the config's `output_dir = "examples/example_demography/example_output"`
is relative to where you launch `meds_main`):

```bash
LD_LIBRARY_PATH=$CONDA_PREFIX/lib ./build-ifx/meds_main examples/example_demography/example_config_main.toml
python post_proc/plot_site_timeseries.py  examples/example_demography/example_output/example_output-D-output.nc \
       -o examples/example_demography/example_output.png
python post_proc/plot_forest_structure.py examples/example_demography/example_output/example_output-D-output.nc \
       -o examples/example_demography/example_output_forest.gif            # every year, 3 fps
python post_proc/plot_landscape_3d.py     examples/example_demography/example_output/example_output-D-output.nc \
       -o examples/example_demography/forest3d_landscape.png               # 3D landscape (needs viz extra)
python post_proc/animate_landscape_growth.py examples/example_demography/example_output/example_output-D-output.nc \
       -o examples/example_demography/forest3d_growth.gif                  # 3D growth animation (needs viz extra)
```

## Model output (`example_output/`)

- **`example_output-D-output.nc`** — the diagnostic timeseries (one record per year): full
  cohort/patch/site state with derived diagnostics, each record stamped with its calendar date.
- **`example_output_pft_parameters.csv`** — the per-PFT parameter table actually used by the run
  (one row per PFT: wood density, allometry, growth, mortality-hazard and recruitment parameters), a
  provenance record written automatically alongside the output.

## Figures

- **`example_output.png`** — site-level totals over time: plant number, leaf area index, aboveground
  biomass, basal area, mean DBH, and cohort/patch structure counts.
- **`example_output_pft.png`** — per-PFT aboveground biomass over time, as lines in the classic ED /
  Moorcroft et al. (2001) colours (PFT 1 green, PFT 2 blue, PFT 3 magenta): the pioneer flush, then
  mid-successional dominance, then slow climax accumulation.
- **`example_output_forest.gif`** — the canopy-layer stand profile animated over the 250 years (one
  frame per year, 3 fps). The **left** panel is the site's vertical **LAI profile** at 2 m resolution
  (black stepped line) sharing the height axis with the **right** panel's stand cross-section: each
  cohort is a thin rectangle spanning its patch's full width (the flat canopy disk seen edge-on) at the
  cohort's height, thickness ∝ its LAI, colour = PFT. Patches tile oldest → youngest (width ∝ area) and
  keep stable slots via their persistent `global_patch_id`; the frame title shows the year since start
  and the panel header the total LAI.
- **`forest3d_landscape.png`** — a synthetic **3D landscape** of the whole site (last record):
  patches laid out as a contiguous, area-weighted Voronoi mosaic, each populated with allometric tree
  crowns (PFT 1 green, 2 blue, 3 magenta) shaded by Beer–Lambert light attenuation through the
  overtopping LAI (bright canopy top → dark understory). Needs the optional `viz` extra
  (`pyvista`, `scipy`, `netCDF4`); see `post_proc/plot_landscape_3d.py`.
- **`forest3d_growth.gif`** — the 3D landscape animated over the full 250-year spin-up (every 2nd
  year). Each cohort is tracked by its persistent `global_cohort_id`, so trees grow **in place**:
  positions are assigned in a backward pass (last record first, so the mature forest gets the cleanest
  layout) and recruits are scattered by a double-Poisson process, then frames play forward. The
  succession is unmistakable — bare ground → PFT-1 (green) pioneer flush → PFT-2 (blue) mid-
  successional canopy → PFT-3 (magenta) climax understory under blue emergents. Same optional `viz`
  extra; see `post_proc/animate_landscape_growth.py`.

## Census restart input

**`census_example.csv`** — a pseudo cohort census (one row per cohort:
`site_id,patch_id,cohort_id,dbh,height,pft,nplant`) used to start a run from existing stand
structure instead of bare ground. Point `[init].census_file` at it with `[init].init_mode = 1`. See
[`src/init/meds_init.f90`](../../src/init/meds_init.f90) (`init_from_census`).

## The empirical golden was recaptured (2026-09-10)

`test/golden/empirical_spinup_golden.csv` was originally taken by hand from the **original Fortran
empirical model**, which the reorg deleted — the empirical laws now live in `empirical_laws.py` and
drive the Fortran engine's law-free apply-primitives through the C-API. PR #137 then changed how the
recruit pool accrues: it is credited **every slow step** rather than as one monthly lump, so the pool
no longer receives a full month's recruits on day 1 before any time has passed. The first cohorts
therefore appear about a month later, and the whole early trajectory is offset.

The offset decays as the stand fills, and the trajectories converge:

| year | cohorts (P0) | cohorts (new) | AGB P0 | AGB new | rel. diff |
|---|---|---|---|---|---|
| 2  | 7   | 6   | 1.529e-03 | 4.607e-04 | 69.9 % |
| 3  | 12  | 6   | 7.468e-02 | 6.367e-02 | 14.7 % |
| 5  | 35  | 26  | 6.631e-01 | 6.208e-01 | 6.4 % |
| 7  | 83  | 70  | 1.669e+00 | 1.613e+00 | 3.4 % |
| 10 | 124 | 138 | 3.302e+00 | 3.250e+00 | 1.6 % |
| 20 | 229 | 197 | 6.961e+00 | 6.943e+00 | 0.3 % |
| 40 | 261 | 237 | 1.043e+01 | 1.045e+01 | **0.1 %** |

The golden is now captured from the current model, so the example reports zero error. **Recapture is
reproducible**: `python3 examples/example_demography/empirical_spinup.py --emit-golden`. It had no
such path before, which is why it went stale silently.
