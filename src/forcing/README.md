# `forcing/` — prescribed external drivers

The home for **time-varying boundary conditions read from a file**, as opposed to state the model
evolves. Meteorology and prescribed CO₂ today; disturbance and land-use schedules and nitrogen
deposition later.

`libmeds_forcing` links the configuration library and the netCDF C bindings **only** — never the
demography or state layer — so a prescribed driver stays low in the library graph.

## Modules

- **`meds_forcing_types`** — the runtime types. `met_forcing_t` is the instantaneous per-site
  atmospheric state the fast loop consumes: a read-only boundary-condition value, which carries the
  wind vector as well as the speed when the source supplies components. `met_record_t` is one raw
  file record. The reader state comes in two parts: `met_source_t`, one per run, holds the file,
  the time axis, the cell list and the loaded month, and is read-only while a step runs;
  `met_cursor_t`, one per polygon, holds the polygon's cell and location and the two records that
  bracket its model time. `met_cells_t` lists the archive cells a run reads (one for a site, the
  valid cells of a box), and `met_month_t` holds one month of them. `co2_series_t` is a prescribed
  CO₂ series, held on `met_source_t`.
- **`meds_forcing_kernels`** — the `pure` and `elemental` math: per-variable temporal interpolation
  (linear or step) with an energy-conserving form for wind, the local apparent-solar-time transform
  (UTC plus longitude plus the equation of time), the **interval-mean-conserving** shortwave
  disaggregation, the total-to-four-stream shortwave partition (Erbs clearness index by default,
  Weiss-Norman available), humidity conversions over the shared saturation vapour pressure,
  precipitation phase, and nearest-grid matching.
- **`meds_lapse_rate`** — every vertical correction, in two steps (`docs/science/forcing.md` §8).
  The **terrain** lapse runs per record at ingest (`read_record`), from the forcing cell's elevation
  to the site's: temperature by a monthly lapse rate, pressure hydrostatically, humidity at constant
  relative humidity, and file longwave by the clear-sky ε·T⁴ ratio. The **move to each patch's
  canopy-air top** runs per patch in the fast loop (`met_to_cas_top`), from the forcing's declared
  heights: potential temperature and humidity are conserved, and the wind follows the patch's own
  log profile after an open-terrain wind is returned to its blending height.
- **`meds_met_driver`** — the reader. `met_open` opens a source (for a site, or for a region's
  cells), `met_cursor_init` places a polygon's cursor on it, `met_prefetch` loads what the next step
  reads, and `met_advance` / `met_instant` step and sample a cursor with no file access. Two file
  sources, chosen by `[forcing].format`: the MEDS multi-grid `(time, grid)` forcing NetCDF
  (`"netcdf"`), and the global ED_ERA5land archive (`"era5land"`), whose months it lays end to end
  as one hourly axis so bracketing and recycling are the same code for both. Dewpoint becomes
  specific humidity and the wind components become the speed at each stamp; shortwave is
  partitioned at ingest. Also the no-file constant-climate backend, used by the tests.
- **`meds_co2_series`** — the prescribed CO₂ (#184): reads a MEDS CO₂ file (format 1, a plain-text
  list of period means at a declared `timestep`) and looks it up at a model instant, linear between
  period middles. `met_open` reads it once and `met_instant` sets `met%co2` from it, or from
  `co2_const`, for every backend alike; the met file never carries CO₂. The shipped series is
  `data/co2/`, built by `scripts/prepare_co2/make_co2_file.py`.
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
must cover the window. The mapping is anchor-relative, so a cycle may begin at any record stamp,
not only on 1 January (a region's at 00:00 or 01:00 on the 1st), and hour-of-day is preserved
exactly. A window that is not a whole
number of years drifts both hour-of-day and day-of-year on every wrap **while the daily mean stays
correct** — so nothing downstream complains, and a multi-decade run can end up reading the wrong
season at the wrong hour with a perfectly healthy-looking energy budget.

## How it reaches the model

Gated on `[forcing].forcing_on`. When on, the run's driver (`meds_driver` for a site, `meds_region`
for a region) opens one `met_source_t` and gives each polygon a `met_cursor_t`, and each step's
forcing is loaded before the step (`met_prefetch`), so no step reads a file. The fast loop samples
the forcing **once per sub-step** (`met_advance`, `met_instant`), outside the patch loop, and moves
it to each patch's canopy-air top (`meds_lapse_rate`) — so the diurnal cycle lives inside the
sub-step loop.

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
`data_path`. The static file's `valid` mask comes from one raw GDEX file, so an archive built from
the CDS still needs one GDEX download. See `MEDS_FORCING_DESIGN.md` §12–§14.

**A single forcing file** (`format = "netcdf"`) comes from `scripts/prepare_era5/make_forcing_file.py`,
run in the `meds-era5` environment, with either input:

```bash
python scripts/prepare_era5/make_forcing_file.py --data-path <archive> --start ... --end ...  # cut a site from the archive
# or, without an archive, a small download first:
python scripts/prepare_era5/download_era5land_cds.py ...              # a box around the site, from the Copernicus data store
python scripts/prepare_era5/postprocess_era5land.py --split none ...  # decode GRIB: one box file per variable
python scripts/prepare_era5/make_forcing_file.py --box-dir <box files> ...   # de-accumulate, convert, write
```

The full commands are in the header of `make_forcing_file.py`.

The file format, the ERA5-Land de-accumulation recipe (including the hour-zero trap), and all the
disaggregation math are documented in [`docs/science/forcing.md`](../../docs/science/forcing.md).
The design record, now archived, is [`docs/dev_plans/archive/MEDS_FORCING_DESIGN.md`](../../docs/dev_plans/archive/MEDS_FORCING_DESIGN.md).

**Tested** in `test/test_met_driver.f90`: the kernels, the constant backend, a NetCDF round trip
that writes and reads a two-grid file, with and without the wind vector, and the prescribed CO₂
(every format rule, and CO₂ that keeps rising under recycled met). `test/test_met_era5land.f90`
writes a small synthetic archive and covers the archive backend end to end.

## Not here yet

The full multi-polygon runtime, and latitude-resolved CO₂ (one global series drives every polygon).
See [`docs/ROADMAP.md`](../../docs/ROADMAP.md) §8.
