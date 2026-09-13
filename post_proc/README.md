# `post_proc/` — plotting the output

Four scripts that read a MEDS diagnostic netCDF file and draw it. They need **numpy**,
**matplotlib** and **netCDF4**; the two 3-D renderers also need the optional `viz` extra
(`pip install -e python/[viz]`, which adds **pyvista** and **scipy**).

All four take the diagnostic file and an output path:

```bash
python post_proc/plot_site_timeseries.py     out-D-output.nc -o timeseries.png
python post_proc/plot_forest_structure.py    out-D-output.nc -o forest.gif
python post_proc/plot_landscape_3d.py        out-D-output.nc -o landscape.png   # viz extra
python post_proc/animate_landscape_growth.py out-D-output.nc -o growth.gif      # viz extra
```

## What each one draws

**`plot_site_timeseries.py`** — the site totals over time: plant number, leaf area index,
aboveground biomass, basal area, mean diameter, and cohort and patch counts. Plus a per-PFT
aboveground-biomass plot showing successional composition, in the classic ED colours (PFT 1 green,
2 blue, 3 magenta, after Moorcroft et al. 2001).

**`plot_forest_structure.py`** — the stand structure, animated. Two panels sharing a height axis:
on the left the site's vertical LAI profile in 2 m layers, on the right a pseudo-spatial stand
cross-section. Following MEDS's flat-canopy assumption, each cohort is a thin horizontal rectangle
spanning its patch's full width — the canopy disk seen edge-on — at the cohort's height, with
thickness proportional to its LAI and colour by PFT. Patches tile oldest to youngest with width
proportional to area, and hold stable slots across frames via their persistent patch id.

**`plot_landscape_3d.py`** — a synthetic 3-D landscape of the whole site at one record. Patches are
laid out as a contiguous, area-weighted Voronoi mosaic, each populated with allometric tree crowns
shaded by Beer-Lambert attenuation through the overtopping leaf area: bright canopy top, dark
understory, and no cast-shadow artefacts.

**`animate_landscape_growth.py`** — that landscape over the whole run, as a GIF. Every cohort is
tracked by its persistent cohort id so **trees grow in place**. Positions are assigned in a
backward pass — last record first, so the mature forest gets the cleanest layout — with recruits
scattered into the gaps around it by a double-Poisson process; frames are then written forward.

## What makes cohort tracking possible

Every cohort and patch carries a `global_id`, stamped at creation and carried in lockstep through
every sort, fusion and compaction. Ids are never reused, and a fusion keeps the survivor's. That is
what lets a reader follow one plant across output records until it fuses away or is culled — and it
is why the growth animation can put a tree in the same place every frame.

The diagnostic file is ragged: one record per output interval, with `cohort_offset` and
`cohort_count` giving the patch-to-cohort map for each record.

## Reading the file yourself

The output is plain netCDF with CF-style metadata, so anything that reads netCDF will open it. Two
conventions worth knowing:

- **Which variables exist depends on the config.** Run `meds_main --dump-io-config` against your
  run's config to list every available name.
- **Read the coordinates, not your memory of the setup.** Soil depths come from the file's own
  `soil_z` coordinate and the reference CO₂ from `atm_co2_fast`, so a figure stays correct if a run
  changes its soil grid or its forcing. The biophysics example carries a cautionary note about what
  happens when a figure assumes a number instead.

The variable set, the axes and the aggregation rules are documented in
[`docs/science/diagnostics.md`](../docs/science/diagnostics.md).
