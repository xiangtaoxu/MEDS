# `forcing/` — prescribed external drivers

The home for **time-varying boundary conditions read from a file**, as opposed to state the model
evolves. Meteorology today; disturbance and land-use schedules and prescribed CO₂ or nitrogen
deposition later.

`libmeds_forcing` links the shared foundation and the netCDF C bindings **only** — never the
demography or state layer — so a prescribed driver stays low in the library graph.

## Modules

- **`meds_forcing_types`** — the runtime types. `met_forcing_t` is the instantaneous per-site
  atmospheric state the fast loop consumes: a read-only boundary-condition value, which carries the
  wind vector as well as the speed when the source supplies components. `met_record_t` is one raw
  file record. The reader state comes in two parts: `met_source_t`, one per run, holds the file,
  the time axis, the cell list and the loaded month, and is read-only while a step runs;
  `met_cursor_t`, one per polygon, holds the polygon's cell and location and the two records that
  bracket its model time. `met_cells_t` lists the archive cells a run reads (one for a site, the
  valid cells of a box), and `met_month_t` holds one month of them.
- **`meds_forcing_kernels`** — the `pure` and `elemental` math: per-variable temporal interpolation
  (linear or step) with an energy-conserving form for wind, the local apparent-solar-time transform
  (UTC plus longitude plus the equation of time), the **interval-mean-conserving** shortwave
  disaggregation, the total-to-four-stream shortwave partition (Erbs clearness index by default,
  Weiss-Norman available), humidity conversions over the shared saturation vapour pressure,
  precipitation phase, nearest-grid matching, and the wind-height and elevation lapse corrections.
- **`meds_met_driver`** — the reader. `met_open` opens a source (for a site, or for a region's
  cells), `met_cursor_init` places a polygon's cursor on it, `met_prefetch` loads what the next step
  reads, and `met_advance` / `met_instant` step and sample a cursor with no file access. Two file
  sources, chosen by `[forcing].format`: the MEDS multi-grid `(time, grid)` forcing NetCDF
  (`"netcdf"`), and the global ED_ERA5land archive (`"era5land"`), whose months it lays end to end
  as one hourly axis so bracketing and recycling are the same code for both. Dewpoint becomes
  specific humidity and the wind components become the speed at each stamp; shortwave is
  partitioned at ingest. Also the no-file constant-climate backend, used by the tests.
- **`meds_era5land_reader`** — the archive's files: path templates, the static file, site selection
  (nearest valid cell within `max_distance_km`) and box selection (across 180°), and a month of
  every domain cell read one chunk column at a time. It returns status codes, so each rejection is
  testable; the reader turns them into hard errors.

The `[forcing]` and `[site]` config type and all its selector codes live in
`src/config/meds_forcing_config.f90`, so `meds_config` — the root of the dependency graph — can carry
them with no back-edge into this library.

## Two rules worth knowing before you use it

**MEDS never gap-fills.** A missing or NaN required value is a hard error. Filling a gap silently
produces a run that looks healthy and is not, and there is no way for a downstream consumer to tell.

**The recycle window is declared, never inferred.** If `recycle = true`, then `recycle_start` and
`recycle_end` are required, and three things are checked rather than guessed: the span must be an
exact whole number of calendar years, the start must land exactly on a record stamp, and the file
must cover the window. The mapping is anchor-relative, so a cycle may begin anywhere in the
calendar, not only on 1 January, and hour-of-day is preserved exactly. A window that is not a whole
number of years drifts both hour-of-day and day-of-year on every wrap **while the daily mean stays
correct** — so nothing downstream complains, and a multi-decade run can end up reading the wrong
season at the wrong hour with a perfectly healthy-looking energy budget.

## How it reaches the model

Gated on `[forcing].forcing_on`. When on, `meds_main` opens the reader and threads it plus the
step-start time down to the fast loop, which refreshes a local context overlay **per sub-step** — so
the diurnal cycle lives inside the sub-step loop.

Forcing also drives the **canopy radiative transfer**: the met shortwave streams map onto the
two-stream's radiation record, the height-descending cohort gather order is reversed into the
two-stream's bottom-to-top contract, and the result is scattered back per cohort as absorbed
shortwave for the leaf energy balance plus an **incident-equivalent** photosynthetically active
radiation for gas exchange. The distinction matters: the leaf kernel re-applies leaf absorptance
internally, so feeding it raw absorbed visible radiation would count absorptance twice and cost
about 15 % of light-limited photosynthesis. Net longwave is wired the same way, with the canopy
emission temperature set to the canopy-air temperature so that leaf emission is counted exactly once
against the energy balance's own linearization.

The same bottom-to-top contract governs the aerodynamics call, so the in-canopy wind cascade runs in
the right direction for multi-cohort patches.

## Preparing forcing

**The ED_ERA5land archive** (`format = "era5land"`) is built once per installation with the tools in
`scripts/prepare_era5/` (`download_era5land_gdex.py` or `download_era5land_cds.py`, then
`build_era5land_static.py` and `build_era5land_archive.py`); a run then only names its folder in
`data_path`. See `MEDS_FORCING_DESIGN.md` §12–§14.

**A single forcing file** (`format = "netcdf"`) takes three steps. The first two are the ERA5-Land
tools in `scripts/prepare_era5/`, run in its `meds-era5` environment:

```bash
python scripts/prepare_era5/download_era5land_cds.py ...           # a box around the site, from the Copernicus data store
python scripts/prepare_era5/postprocess_era5land.py --split none ...  # decode GRIB: one box file per variable
python scripts/prep_era5land_forcing.py --in <box files> ...        # de-accumulate, convert, write the MEDS format
```

The full commands are in the header of `scripts/prep_era5land_forcing.py`.

The file format, the ERA5-Land de-accumulation recipe (including the hour-zero trap), and all the
disaggregation math are documented in [`docs/science/forcing.md`](../../docs/science/forcing.md).
The design record is [`docs/dev_plans/MEDS_FORCING_DESIGN.md`](../../docs/dev_plans/MEDS_FORCING_DESIGN.md).

**Tested** in `test/test_met_driver.f90`: the kernels, the constant backend, and a NetCDF round trip
that writes and reads a two-grid file, with and without the wind vector. `test/test_met_era5land.f90`
writes a small synthetic archive and covers the archive backend end to end.

## Not here yet

LWdown synthesis (the `"synthesize"` option is rejected by the config validator until it exists),
the full multi-polygon runtime, and a transient CO₂ stream. See
[`docs/ROADMAP.md`](../../docs/ROADMAP.md) §8.
