# forcing

The home for **prescribed external drivers** — time-varying boundary conditions read from a file, as
opposed to state the model evolves. Meteorology today; disturbance-event / land-use schedules and
prescribed CO₂ / N-deposition streams later. `libmeds_forcing` links `meds_shared` + the netCDF C
bindings (`meds_netcdf_c`) **only** — never the demography/state layer — so a prescribed driver stays
low in the library DAG. Design: `archive/MEDS_FORCING_DESIGN.md`.

**Implemented — P0 meteorological forcing** (single-site NetCDF reader, ERA5-Land at Ithaca NY):

- **`meds_forcing_types`** — the runtime types: `met_forcing_t` (the instantaneous per-site atmospheric
  state the fast loop consumes — a read-only boundary-condition value, the analogue of `rad_forcing_t` /
  `chydro_forcing_t`; defaults sum the four SW streams to 400 W/m² = the current CONST climate),
  `met_record_t` (one raw file record), and the **mutable** per-polygon reader buffer `met_driver_t`.
- **`meds_forcing_kernels`** — `pure`/`elemental` math: per-variable temporal `interpolate_forcing`
  (linear / step) + energy-conserving wind, the local **apparent-solar-time** transform (UTC + longitude
  + equation-of-time), the **interval-mean-conserving** shortwave disaggregation
  (`cosz_reconstruct_factor` returns `1/⟨cosz⟩_win`, *not* `⟨sec z⟩`), the total→(beam/diffuse)×(PAR/NIR)
  **`partition_shortwave`** (Erbs clearness-index), `dewpoint_to_specific_humidity` /
  `rh_to_specific_humidity` (reuse `meds_thermo`'s Bolton `esat`), and `precip_phase`.
- **`meds_met_driver`** — the reader: `met_open` / `met_advance` / `met_instant` / `met_close` over the
  **MEDS multi-grid `(time, grid)` forcing NetCDF** (per-polygon `grid_index` hyperslab read; base time
  from the `time:units` attribute; SW partitioned at ingest from total `SWdown`) + the no-file **CONST**
  reference-climate backend. **MEDS never gap-fills** — a missing/NaN required value is a hard error
  (`assert_finite`).

The `[forcing]`/`[site]` config type `forcing_config_t` + all selector codes live in **`src/shared`**
(`meds_forcing_config`) so `meds_config` (the DAG root) can carry it with no `shared → forcing` back-edge.
Tested in `test/test_met_driver.f90` (kernels, CONST backend, and a NetCDF round-trip that writes and reads
a `(time=25, grid=2)` file); green under ifx and nvfortran multicore.

**Wired into the fast loop** (the P0 driver integration, design §6; **opt-in**): the `[forcing]`/`[site]`
TOML block loads via `meds_config_io` — **gated on `forcing.forcing_on`** (a defaulted read, so a config
without the block runs the constant-forcing MVP unchanged). `meds_main` opens the reader when `forcing_on`
and threads it + the step-start time through `advance_one_step` → `run_fast_biophysics`, which refreshes a
local `fast_context_t` overlay **per sub-step** via `apply_met_to_ctx(met_instant(...))` — so the diurnal
cycle lives inside the sub-step loop, and the constant path stays bit-identical. `test/test_fast_loop.f90`
drives the loop from a real forcing NetCDF and asserts the diurnal signal (night SW=0 → GPP≈0, day → GPP>0).
The two ERA5-Land prep scripts (`scripts/download_era5land.py`, `scripts/prep_era5land_forcing.py`) produce
the file the reader consumes.

**P1/P2** (design §8): the canopy-RT join (per-cohort absorbed PAR from `meds_canopy_radiation` replaces the
LAI-share split of `rad_sw_top`), Weiss–Norman band-specific SW, multi-year cycling (solar-geometry / Feb-29
alignment), multi-polygon runtime (nearest-location `grid_index` match), and lapse / wind-height corrections.
