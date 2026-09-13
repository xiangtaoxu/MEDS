# `forcing/` — prescribed external drivers

The home for **time-varying boundary conditions read from a file**, as opposed to state the model
evolves. Meteorology today; disturbance and land-use schedules and prescribed CO₂ or nitrogen
deposition later.

`libmeds_forcing` links the shared foundation and the netCDF C bindings **only** — never the
demography or state layer — so a prescribed driver stays low in the library graph.

## Modules

- **`meds_forcing_types`** — the runtime types. `met_forcing_t` is the instantaneous per-site
  atmospheric state the fast loop consumes: a read-only boundary-condition value. `met_record_t` is
  one raw file record, and `met_driver_t` is the mutable per-polygon reader buffer holding the two
  records that bracket the model time.
- **`meds_forcing_kernels`** — the `pure` and `elemental` math: per-variable temporal interpolation
  (linear or step) with an energy-conserving form for wind, the local apparent-solar-time transform
  (UTC plus longitude plus the equation of time), the **interval-mean-conserving** shortwave
  disaggregation, the total-to-four-stream shortwave partition (Erbs clearness index by default,
  Weiss-Norman available), humidity conversions over the shared saturation vapour pressure,
  precipitation phase, nearest-grid matching, and the wind-height and elevation lapse corrections.
- **`meds_met_driver`** — the reader. `met_open` / `met_advance` / `met_instant` / `met_close` over
  the MEDS multi-grid `(time, grid)` forcing NetCDF, with a per-polygon hyperslab read, the base time
  taken from the `time` variable's units attribute, and shortwave partitioned at ingest. Also the
  no-file constant-climate backend, used by the tests.

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

## Preparing a forcing file

Two scripts produce the file the reader consumes:

```bash
python scripts/download_era5land.py     # fetch from the Copernicus data store
python scripts/prep_era5land_forcing.py # de-accumulate, convert, write the MEDS format
```

The file format, the ERA5-Land de-accumulation recipe (including the hour-zero trap), and all the
disaggregation math are documented in [`docs/science/forcing.md`](../../docs/science/forcing.md).
The design record is [`docs/dev_plans/MEDS_FORCING_DESIGN.md`](../../docs/dev_plans/MEDS_FORCING_DESIGN.md).

**Tested** in `test/test_met_driver.f90`: the kernels, the constant backend, and a NetCDF round trip
that writes and reads a two-grid file. Green under ifx and nvfortran multicore.

## Not here yet

LWdown synthesis (the `"synthesize"` option is rejected by the config validator until it exists),
the full multi-polygon runtime, and a transient CO₂ stream. See
[`docs/ROADMAP.md`](../../docs/ROADMAP.md) §8.
