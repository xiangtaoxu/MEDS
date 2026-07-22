# forcing

The home for **prescribed external drivers** — time-varying boundary conditions read from a file, as
opposed to state the model evolves. Meteorology today; disturbance-event / land-use schedules and
prescribed CO₂ / N-deposition streams later. `libmeds_forcing` links `meds_shared` + the netCDF C
bindings (`meds_netcdf_c`) **only** — never the demography/state layer — so a prescribed driver stays
low in the library DAG. Design: `docs/dev_plans/MEDS_FORCING_DESIGN.md`.

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
  `rh_to_specific_humidity` (reuse `meds_therm_lib`'s Bolton `esat`), and `precip_phase`.
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
and threads it + the step-start time through `advance_one_step` → `fast_dynamics`, which refreshes a
local `fast_context_t` overlay **per sub-step** via `apply_met_to_ctx(met_instant(...))` — so the diurnal
cycle lives inside the sub-step loop, and the constant path stays bit-identical. `test/test_fast_loop.f90`
drives the loop from a real forcing NetCDF and asserts the diurnal signal (night SW=0 → GPP≈0, day → GPP>0).
The two ERA5-Land prep scripts (`scripts/download_era5land.py`, `scripts/prep_era5land_forcing.py`) produce
the file the reader consumes.

**Wired into the canopy radiative transfer** (the P1 RT join, design §6.3): when forcing is on, the fast
loop's per-sub-step `apply_rt_forcing` (`meds_fast_dynamics`) replaces the LAI-share shortwave split with the
real two-stream `meds_canopy_radiation.canopy_radiation`. It maps the met SW streams to `rad_forcing_t`
(`par_beam`/`par_diffuse`→VIS, `nir_beam`/`nir_diffuse`→NIR, `lwdown`→LW-diffuse), reverses the
height-DESCENDING cohort gather order into the two-stream's **BOTTOM(1)→TOP(n)** contract via an
ascending-height permutation, and inverse-scatters the result back per cohort: `abs_sw` = absorbed VIS+NIR
(leaf energy), plus below-canopy ground SW (the **net** `dn_ground − up_ground`, so soil albedo is
respected — a bare `ncoh=0` patch flows through the same `canopy_radiation` empty-canopy branch, so the
1→0 transition is continuous). For photosynthesis, `abs_par` is the two-stream absorbed VIS divided by the
leaf PAR absorptance (`par_per_w = 4.6` µmol W⁻¹) — i.e. an **incident-equivalent** PAR, because the leaf
gas-exchange kernel re-applies `leaf_absorptance` internally; feeding it raw absorbed VIS would count leaf
absorptance twice (~15 % low light-limited GPP). The per-patch forcing buffers resize per patch (a later
patch may hold more cohorts than the first).
PFT-uniform optics + soil albedo/emissivity are built once into `fast_context_t` by `build_fast_context`.
The same BOTTOM→TOP contract governs `canopy_aerodynamics`, so the sibling `aero_bottom_to_top`
(`meds_fast_ark`) now feeds it the reversed order too (fixing a latent wind-cascade inversion the
old direct call caused for multi-cohort patches). **Net longwave** is wired too: the two-stream's per-cohort
net leaf LW (`abs_leaf(RAD_LW,·)`) and net ground LW (`dn_ground − up_ground`) feed the leaf/ground energy
balance. The two-stream's canopy LW emission temperature is set to the **canopy-air temperature `tcas`**,
because the diagnostic leaf balance linearizes emission around `tcas` (`lw_slope·dtl`) and so needs
`abs_lw` = net-LW-*at-tcas* — this makes leaf emission count exactly once (no double-count with the
balance's own emission response). The ground balance carries no separate emission term, so the net
`dn_ground − up_ground` (soil emission baked in at `soil_temp`) is directly consistent.
Guarded by `test_canopy_radiation` (top-of-two-equal-LAI-cohorts absorbs more), `test_column_dynamics`
(aero order), and `test_fast_loop` (multi-cohort shading, hydraulics-neutralized so the light/gb ordering
drives GPP — with water stress ON, the taller cohort's more negative `psi_leaf` correctly wins the
GPP-per-leaf inversion; plus a night radiative-cooling check that the net-LW loss to a cold sky pulls the
canopy air below the atmosphere).

**P2 — implemented** (design §8): **Weiss–Norman** band-specific SW partition (`SWPART_WEISS_NORMAN`, ED2
`short_bdown_weissnorman` port; `partition_shortwave` gained a `psurf_pa` arg; Erbs stays the default and is
unchanged); **multi-year calendar recycling + Feb-29** (`file_lookup_sec` maps a whole-year Jan-1-aligned
file's records to the model's calendar year preserving day-of-year, `Feb-29 → Feb-28` when the file year is
non-leap; the SW reconstruction factor is anchored on the **model** window so the interval-mean-conserving
identity follows the model sun; a non-calendar file keeps the legacy absolute-seconds span-wrap);
**nearest-grid match** (`grid_match = "nearest"` binds the site to the file `grid` cell minimizing
great-circle distance from `[site]` lat/lon — the reusable atom of a future full multi-polygon runtime, which
stays deferred since MEDS is single-site); and **wind-height + elevation lapse** (opt-in `apply_wind_profile`
neutral-log wind lift to the reference height, `apply_elevation_lapse` hydrostatic T/P lapse from the
grid-cell to the site elevation), plus the ED2 `reference_height > hgt_max` config guard.

**Deferred**: the full multi-polygon runtime (a grid→polygon→site state hierarchy + array of `met_driver_t`
+ polygon loop + MPI — a large change orthogonal to forcing); LWdown synthesis (ERA5-Land ships `strd`); the
daily accumulator (`temp_day`/`daylength`/`doy`) for phenology.
