# Demography example

A self-contained **demographic spin-up** of MEDS and the figures it produces. The configuration is a
pair of files: [`example_config_main.toml`](example_config_main.toml) (run/engine/IO settings, which
names the PFT file) and [`example_config_pft.toml`](example_config_pft.toml) (PFT traits + allometry +
mortality coefficients) — a **250-year** daily spin-up (2000-01-01 → 2250-01-01) from near-bare
ground, writing to `example_output/` with the prefix `example_output`. (For the standalone leaf-level
photosynthesis example, see [`../example_leaf_physiology/`](../example_leaf_physiology/).)

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

## Census restart input

**`census_example.csv`** — a pseudo cohort census (one row per cohort:
`site_id,patch_id,cohort_id,dbh,height,pft,nplant`) used to start a run from existing stand
structure instead of bare ground. Point `[init].census_file` at it with `[init].init_mode = 1`. See
[`src/init/meds_init.f90`](../../src/init/meds_init.f90) (`init_from_census`).
