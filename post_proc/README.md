# `post_proc/` — plotting and evaluating the output

Five scripts and a notebook that read MEDS output and make sense of it. They need **numpy**,
**matplotlib** and **netCDF4**; the two 3-D renderers also need the optional `viz` extra
(`pip install -e python/[viz]`, which adds **pyvista** and **scipy**), and the notebook needs
**jupyter**.

The scripts take a file and an output path:

```bash
python post_proc/plot_site_timeseries.py     out-D-output.nc -o timeseries.png
python post_proc/plot_forest_structure.py    out-D-output.nc -o forest.gif
python post_proc/plot_pft_size.py            'out/*-M-*.nc'  -o pft_size.png
python post_proc/plot_landscape_3d.py        out-D-output.nc -o landscape.png   # viz extra
python post_proc/animate_landscape_growth.py out-D-output.nc -o growth.gif      # viz extra
```

The notebook takes an output **directory**:

```bash
MEDS_RUN=runs/ithaca_ark30/out jupyter notebook post_proc/evaluate_ithaca.ipynb
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

**`plot_pft_size.py`** — the two axes that describe *who* the stand is made of rather than how much
of it there is: composition by PFT, and structure by DBH class. Six panels — biomass and leaf area
by PFT over time, the size distribution first-against-last, biomass by class stacked over the run,
and the vital rates against size. It is the reference consumer for the PFT and size axes, so it also
prints the closure identities those axes promise (`Σ_pft agb_pft == agb_site`,
`Σ_class agb_size == agb_site`) before drawing anything: a figure built on axes that do not
partition the stand is a picture of a bug.

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

## The evaluation notebook

**`evaluate_ithaca.ipynb`** is a worked pass over one run, using each tier for the question that
tier can answer: diurnal energy and carbon from **FAST**, seasonal cycles from **DAILY**,
demographic trajectories and the closure identities from **MONTHLY/ANNUAL**.

It is not a benchmark. MEDS has never been scored against flux-tower or inventory data, and the
notebook does not pretend to. What it checks is what is knowable without observations: identities
the output must satisfy, and magnitudes and shapes that are either physically possible or are not.
That is still enough to catch a great deal — the first run of it turned up three defects (#245,
#246, #247) and one documented-behaviour surprise.

Two things it teaches that are easy to get wrong when reading these files by hand:

- **Turning a rate back into an amount needs the width of its own window**, `t[i+1] - t[i]`, and the
  trailing partial window has to be dropped because its width is unknown. Pairing a rate with the
  *previous* window's width is silently fine until the end of a run, where it turned one disturbance
  event into a 10.8 kgC/m² phantom.
- **`litter_*_site` is not all the carbon entering the soil.** The cull and disturbance mortality
  pathways add necromass straight to the patch soil-carbon pools, bypassing the litter accumulator
  the litter variables read. Closing the soil-carbon budget from file needs the
  `mort_carbon_*_site` variables as well; with them the budget closes to 0.3%, without them it
  misses by 29%.
