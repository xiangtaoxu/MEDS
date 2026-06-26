# Examples

Reference outputs and inputs for MEDS.

## Reference figures (default 250-year spin-up)

From a default 250-year spin-up ([`meds_config.toml`](../meds_config.toml), daily time step),
produced by the netCDF output + post-processing pipeline:

```bash
LD_LIBRARY_PATH=$CONDA_PREFIX/lib ./build-io/meds_main meds_config.toml      # -> ./meds_output.nc
python post_proc/plot_site_timeseries.py meds_output.nc -o examples/meds_timeseries.png
```

- **`meds_timeseries.png`** — site-level totals over time: plant number, leaf area index,
  aboveground biomass, basal area, mean DBH, and cohort/patch structure counts.
- **`meds_timeseries_pft.png`** — per-PFT aboveground-biomass stack showing the successional
  composition (pioneer → mid → climax, the wood-density growth–mortality trade-off in action).

## Forest-structure animation (200-year spin-up)

**`meds_forest_structure.gif`** — a canopy-layer stand profile animated over a 200-year spin-up
(the full pioneer → mid → climax succession). Following MEDS's flat-canopy assumption, each cohort
is a thin horizontal rectangle spanning its patch's full width (the canopy disk seen edge-on) at the
cohort's height, with thickness ∝ its LAI and colour = PFT; the y-axis is height. Patches tile
left → right oldest → youngest (width ∝ area) and keep stable slots via their persistent
`global_patch_id`. Regenerate with (`--stride 2` plots every 2nd year to keep the GIF small):

```bash
# run200.toml: [run] years=200 ; [io] output_prefix="state200"
LD_LIBRARY_PATH=$CONDA_PREFIX/lib ./build-io/meds_main run200.toml          # -> ./state200.nc
python post_proc/plot_forest_structure.py state200.nc -o examples/meds_forest_structure.gif \
       --fps 8 --stride 2 --dpi 88
```

## Census restart input

**`census_example.csv`** — a pseudo cohort census (one row per cohort:
`site_id,patch_id,cohort_id,dbh,height,pft,nplant`) used to start a run from existing stand
structure instead of bare ground. Point `[init].census_file` at it in the config. See
[`src/init/meds_init.f90`](../src/init/meds_init.f90) (`init_from_census`).
