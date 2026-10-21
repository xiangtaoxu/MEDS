# `src/io/met/` — the atmosphere above the canopy

The home for **time-varying boundary conditions read from a file**, as opposed to state the model
evolves. Meteorology and prescribed CO₂ today; disturbance and land-use schedules and nitrogen
deposition later.

`libmeds_met` links the configuration library and the netCDF C bindings **only** — never the
demography or state layer — so a prescribed driver stays low in the library graph.

## Modules

- **`meds_met_types`** — the runtime types. `met_forcing_t` is the instantaneous atmospheric
  state the fast loop reads (a polygon's sample, and a site's or a patch's after the vertical
  moves): a read-only boundary-condition value, which carries the wind vector as well as the speed
  when the source supplies components; its one method, `swdown()`, sums the four shortwave streams.
  `met_record_t` is one ingested file record. The reader state comes in two parts: `met_source_t`,
  one per run, holds the file, the time axis, the cell list and the records in memory, and is
  read-only while a step runs; `met_cursor_t`, one per polygon, holds the polygon's cell and
  location and the two records that bracket its model time. `met_grid_t` describes a gridded
  product's grid (its static file, the names in it, and the chunk edge its data files are stored
  in); `met_cells_t` lists the cells of it a run reads (one for a site, the valid cells of a box),
  and `met_month_t` holds one month of them. `co2_series_t` is a prescribed CO₂ series, held on
  `met_source_t`.
- **`meds_forcing_kernels`** — the `pure` and `elemental` math: per-variable temporal interpolation
  (linear or step) with an energy-conserving form for wind, the solar zenith `met_solar_cosz` at
  local apparent solar time (UTC plus longitude plus the equation of time, a private step), the
  **interval-mean-conserving** shortwave disaggregation, the total-to-four-stream shortwave partition
  (Erbs clearness index by default, Weiss-Norman available), humidity conversions over the shared
  saturation vapour pressure, precipitation phase, the longwave synthesis and the sun on a slope.
- **`meds_lapse_rate`** — every vertical correction, in two steps (`docs/science/forcing.md` §8).
  The **terrain** lapse runs per record at ingest (`read_record`), from the forcing cell's elevation
  to the polygon's: temperature by a monthly lapse rate, pressure hydrostatically, humidity at
  constant relative humidity, and file longwave by the clear-sky ε·T⁴ ratio. Per sub-step, the same
  lapse takes the polygon's sample to each of its sites, with the beam on the site's slope
  (`met_at_site`). The **move to each patch's canopy-air top** runs per patch in the fast loop
  (`met_to_cas_top`), from the forcing's declared heights: potential temperature and humidity are
  conserved, and the wind follows the patch's own log profile after an open-terrain wind is returned
  to its blending height.
- **`meds_met`** — the reader. `met_open` opens a source (for a site, or for a region's
  cells), `met_cursor_init` places a polygon's cursor on it, `met_prefetch` loads what the next step
  reads, and `met_advance` / `met_instant` step and sample a cursor with no file access. Two file
  sources, chosen by `[forcing].format`, each an open and one record's values as stored:
  - `meds_met_default`: the MEDS multi-grid `(time, grid)` forcing NetCDF (`"ED_default"`);
  - `meds_met_era5land`: the global ED_ERA5land archive (`"ED_ERA5land"`), whose months it
    lays end to end as one hourly axis.

  Each is named after the `format` value that selects it. Both share `meds_met_share` (the
  status codes, the time axis, the recycle window and the check of a file against its config), so
  bracketing and recycling are the same code for both. A reader that refuses a file returns its
  status and its reason, and `met_open` prints the reason.
  `read_record` is the one ingest after that. Every clock is UTC. The humidity a file carries
  (dewpoint, relative humidity or specific humidity, exactly one) becomes specific humidity through
  the model's own saturation curve, and the wind components become the speed, at each stamp;
  shortwave is partitioned at ingest. Rain and shortwave come from the record whose interval
  contains the instant.
- **`meds_met_co2`** — the prescribed CO₂: reads a MEDS CO₂ file (format 1, a plain-text
  list of period means at a declared `timestep`) and looks it up at a model instant, linear between
  period middles. `met_open` reads it once and `met_instant` sets `met%co2` from it, or from
  `co2_const`, for every backend alike; the met file never carries CO₂. The shipped series is
  `data/co2/`, built by `scripts/prepare_co2/make_co2_file.py`.
- **`meds_met_grid`** — where on the Earth the met data are, for any gridded product described by
  a `met_grid_t`: a site's cell (the nearest valid cell within `max_distance_km`), a region's box of
  cells (also across 180°), the cells a mask keeps, and whether another file is on the grid
  (`same_grid`, which the maps rely on). It reads only the static file, never weather, and returns
  `GRID_*` status codes, so each rejection is testable. `meds_met_era5land` gives the archive's grid
  (`era5land_grid`) and reads its months one chunk column at a time.

The `[forcing]` and `[site]` config type and all its selector codes live in
`src/config/params/meds_forcing_params.f90`, so `meds_config` — the root of the dependency graph —
can carry them with no back-edge into this library.

## Two rules worth knowing before you use it

**MEDS never gap-fills.** A missing or NaN required value is a hard error. Filling a gap silently
produces a run that looks healthy and is not, and there is no way for a downstream consumer to tell.

**The recycle window is declared, never inferred.** If `recycle = true`, then `recycle_start` and
`recycle_end` are required, and three things are checked rather than guessed: the span must be an
exact whole number of calendar years, the start must land exactly on a record stamp, and the file
must cover the window. The mapping is anchor-relative, so a cycle may begin at any record stamp,
not only on 1 January (a region's at 00:00 or 01:00 on the 1st), and hour-of-day is preserved
exactly. A window that is not a whole number of years drifts both hour-of-day and day-of-year on
every wrap **while the daily mean stays correct** — so nothing downstream complains, and a
multi-decade run can end up reading the wrong season at the wrong hour with a perfectly
healthy-looking energy budget.

## How it reaches the model

Gated on `[forcing].forcing_on`. When on, the run (`meds_region`, one polygon or many) opens one
`met_source_t` (`met_open`) and gives each polygon a `met_cursor_t` (`met_cursor_init`). Before
each step the run loads what the step reads (`met_prefetch`), so no step reads a file. The fast
loop samples the forcing **once per sub-step** (`met_advance`, `met_instant`), outside the patch
loop, and moves the sample to each site (`met_at_site`) and then to each patch's canopy-air top
(`met_to_cas_top`) — so the diurnal cycle lives inside the sub-step loop.

How the sampled radiation enters the canopy two-stream (its bottom-to-top layer order, the
incident-equivalent PAR and the longwave emission temperature, in `apply_rt_forcing`) is described in
[`docs/science/canopy_radiation_transfer.md`](../../../../docs/science/canopy_radiation_transfer.md)
§3–§4.

## Preparing forcing

**The ED_ERA5land archive** (`format = "ED_ERA5land"`) is built once per installation with the tools in
`scripts/prepare_era5/` (`download_era5land_gdex.py` or `download_era5land_cds.py`, then
`build_era5land_static.py` and `build_era5land_archive.py`); a run then only names its folder in
`data_path`. The static file's `valid` mask comes from one raw GDEX file, so an archive built from
the CDS still needs one GDEX download. See `MEDS_FORCING_DESIGN.md` §12–§14.

**A single forcing file** (`format = "ED_default"`) comes from `scripts/prepare_era5/make_forcing_file.py`,
run in the `meds-era5` environment, with either input:

```bash
python scripts/prepare_era5/make_forcing_file.py --data-path <archive> --start ... --end ...  # cut a site from the archive
# or, without an archive, a small download first:
python scripts/prepare_era5/download_era5land_cds.py ...              # a box around the site, from the Copernicus data store
python scripts/prepare_era5/postprocess_era5land.py --split none ...  # decode GRIB: one box file per variable
python scripts/prepare_era5/make_forcing_file.py --box-dir <box files> ...   # de-accumulate, convert, write
```

The full commands are in the header of `make_forcing_file.py`.

**From flux-tower data**, `scripts/prepare_flux_tower/make_tower_forcing.py` builds the same kind of
file from AmeriFlux BASE, FLUXNET or a plain CSV, described by a site TOML that declares the source
clock, stamp convention, heights and units; the tool checks each declaration against the sun and the
data. `examples/example04_column_biophysics/` is the worked example.

The file format, the ERA5-Land de-accumulation recipe (including the hour-zero trap), and all the
disaggregation math are documented in [`docs/science/forcing.md`](../../../../docs/science/forcing.md).
The design record is [`docs/dev_plans/archive/MEDS_FORCING_DESIGN.md`](../../../../docs/dev_plans/archive/MEDS_FORCING_DESIGN.md).

**Tested** in `test/test_met_driver.f90`: the kernels, a NetCDF round trip
that writes and reads a two-grid file, with and without the wind vector, and the prescribed CO₂
(every format rule, and CO₂ that keeps rising under recycled met). `test/test_met_era5land.f90`
writes a small synthetic archive and covers the archive backend end to end. `test/test_met_tower.f90`
covers what a flux-tower file relies on: the humidity forms, the UTC and height checks, the interval
each flux is read from, and the round trip of the tower's relative humidity.

## Not here yet

Latitude-resolved CO₂: one global series drives every polygon. See
[`docs/ROADMAP.md`](../../../../docs/ROADMAP.md) §8 for what is planned.
