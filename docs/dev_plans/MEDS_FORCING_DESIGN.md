# MEDS Meteorological Forcing — Source & Wiring Design

The prerequisite that closes MEDS's largest coupling gap. Today MEDS has **no meteorological
forcing source**: the fast biophysics loop is OFF by default, and *every* atmospheric boundary
condition enters the fast kernels as a **constant, horizontally-uniform, in-type-default value**
hand-assembled by `fast_context_t` (`meds_fast_loop.f90`) and translated to the kernel inputs by two
private routines, `build_forcing` (→ `column_forcing_t`) and `fill_aenv` (→ `aero_env_t`). The one
time-derived environmental quantity in the whole model is the solar cosine `solar_cosz` in
`meds_time` — and it is not yet wired to a site latitude, because the TOML has no `[site]` block. The
CLAUDE.md "Reserved follow-ups" name this precisely: *"a meteorological forcing source"* is one of the
four missing pieces (with canopy RT, leaf energy balance, hydraulics) blocking real GPP.

This document designs (1) a per-site **forcing derived type** (`met_forcing_t`) and a small
**buffer/reader state** (`met_driver_t`), living at the polygon/site level — "a location sharing one
meteorological forcing" (CLAUDE.md §Architecture-1) — all now in a **new `src/forcing/` library**
(`libmeds_forcing`, §2), the home for every prescribed external driver (meteorology now; disturbance-event
and land-use schedules later); (2) a stateless **reader module** `src/forcing/meds_met_driver.f90` reading
a **MEDS forcing NetCDF** (the single file format, §7) that is **multi-grid from day one** — a `(time,
grid)` layout where each `grid` index is one polygon's location — aligned to the model calendar via
`meds_time`; (3) **per-variable temporal interpolation/disaggregation** from the forcing interval (hourly
for ERA5-Land) down to `dt_fast` (900 s), with a full ERA5-Land → MEDS-field mapping table and unit
conversions; (4) the **exact wiring** of forcing into the fast loop (the `build_forcing`/`fill_aenv` seam)
and, via daily aggregation, into the slow loop; (5) **two Python scripts** — a CDS-API **downloader** for
ERA5-Land over a requested grid cell and a **formatter** converting it into the MEDS multi-grid forcing
NetCDF — with the **Ithaca NY** grid cell as the reference test site; (6) a **phasing** P0→P1→P2 and a
**test plan** driving the leaf / RT kernels offline to reproduce a diurnal cycle.

The physics reference is **ED2**'s `ed_met_driver.f90` + `radiate_utils.f90` (HDF5 read →
`cgrid%metinput` → temporal interpolation + radiation breakdown + thermodynamic closure →
`cgrid%met` → per-site `cpoly%met`). The community interchange reference is the **ALMA** convention
(8 core drivers in SI, positive-downward flux sign) as realised by FLUXNET/PLUMBER2. As always: MEDS
ports ED2's *algorithms* (cosine-zenith SW disaggregation, Weiss–Norman partition, step-constant
precip) as **stateless `pure`/`elemental` kernels**, keeps the reader's mutable buffer in a
driver/IO-owned state type, and threads config as a passed derived type — never module globals.

> **DATA NOTE (read first).** The reference test forcing is now **ERA5-Land hourly reanalysis**
> (Copernicus CDS dataset `reanalysis-era5-land`, ~0.1° / ~9 km, hourly, UTC), sampled at the grid cell
> covering **Ithaca, NY (~42.44 °N, −76.50 °E; UTC−5 EST / −4 EDT)**. Two Python scripts prepare it
> (§7): `scripts/download_era5land.py` pulls the raw ERA5-Land NetCDF from the CDS API for a requested
> lat/lon (or box), and `scripts/prep_era5land_forcing.py` converts it into the **MEDS multi-grid forcing
> NetCDF** the Fortran reader consumes. The earlier **BCI tower CSV** path is **retired** — the reader now
> reads **only** the MEDS forcing NetCDF (plus the no-file CONST backend); BCI or any FLUXNET tower can
> re-enter later as *another formatter script* that emits the same NetCDF, so the Fortran reader never
> changes when the raw source does. ERA5-Land carries **no** CO₂, and only **total** downward shortwave
> (no direct/diffuse or PAR/NIR split) — so the shortwave-partition kernel (§5.6) and the CO₂ config
> constant (§5.2) move from "reserved" to **P0-required**.

---

## 0. Revision notes (2026-07-08) — three requested changes

This revision folds in three review comments; each is threaded through the affected sections below.

1. **New `src/forcing/` library — the home for all prescribed external drivers (§2).** All met-forcing
   code — the runtime types (`met_forcing_t`, `met_record_t`, `met_driver_t`), the disaggregation/partition
   kernels, and the reader (`meds_met_driver`) — moves out of the old `src/shared`+`src/io` split into one
   dedicated **`src/forcing/`** library (`libmeds_forcing`). This is the future home for **other forcing
   inputs too** (disturbance-event / fire / land-use schedules, prescribed-CO₂ or N-deposition streams).
   *Only* the plain-scalar `[forcing]`/`[site]` **config** stays in `src/shared` (so `meds_config`, the DAG
   root, can carry it without a backward `shared → forcing` edge — the same `soil_thermal_params_t` rule
   the energy design used). Consolidating the reader with its kernels **removes** the old design's central
   contortion — kernels were exiled to `src/shared` *only* to avoid an `io → biophysics` edge when the
   reader lived in `src/io`; with the reader now in its own library that reason evaporates (§2).
   **Folder-name recommendation: `forcing/`, not `input/`.** "Forcing" is the standard land-surface / ESM
   term for prescribed, time-varying *external boundary conditions* (meteorology, CO₂, N-deposition, and
   event schedules like fire / harvest / land-use) — precisely this folder's remit — and it cleanly
   **excludes** static config and initial state (which are `config` / `init`). `input/` is broader but
   muddier (config, restart files, and the PFT table are all "input" too); `drivers/` would collide with
   the existing `src/driver/` time-stepping engine. If you prefer `input/`, it is a pure rename — every
   path below swaps `forcing/ → input/` and nothing else changes.
   **netCDF + IO are now always compiled (comment 1):** the `-DMEDS_ENABLE_IO=OFF` netCDF-free build and
   its `meds_io_stub.f90` are deleted, so netCDF is a hard build dependency everywhere; the forcing library
   therefore needs **no stub reader** (`meds_met_driver_stub`/`meds_netcdf_c_stub` are gone) and always
   links the real `meds_netcdf_c`. The CONST backend remains as a legitimate no-file reader mode, not a
   netCDF-absence fallback (§2).

2. **ERA5-Land replaces the BCI flux file as the test source (§1, §5.2, §5.6, §7, §9).** A CDS-API
   downloader + a formatter script produce the MEDS forcing NetCDF from ERA5-Land for the Ithaca NY grid
   cell. ERA5-Land's specifics reshape three things: it is **hourly** (not 30-min), it is **UTC** (making
   the longitude-based solar-time path the *natural* one, §5.1), and it supplies **only total SWdown + no
   CO₂**, promoting the SW-partition kernel and the CO₂ constant to P0.

3. **NetCDF is THE forcing file format, multi-grid by construction (§2.1, §3, §7).** The MEDS forcing file
   is a single NetCDF with a `(time, grid)` layout: an **unstructured `grid` dimension** whose each index is
   one location (`latitude(grid)`, `longitude(grid)`), generalizing single-site (`grid=1`), a scattered set
   of towers, or a flattened regular grid. The reader selects a per-polygon `grid_index`. A single-site run
   uses `grid=1`; **multi-polygon runtime** (each polygon binding to a grid index by nearest-location match)
   stays P1/P2, but the **file format supports it from P0** so no re-cut is needed. The CSV backend is
   dropped.

4. **MEDS never gap-fills — a missing required value is a hard error (§5.5).** The `gap_policy` selector and
   all in-model gap-filling (persistence hold, climatology) are **removed**: if a required forcing field is
   missing/NaN at a needed step, the reader (and the formatter) **halt** with an error naming the field and
   time. Any gap-filling is the user's responsibility, entirely upstream/outside MEDS. This is orthogonal to
   the *necessary derivations* (SW partition, humidity-from-dewpoint), which compute — not fill — variables
   the source doesn't store directly.

---

## 1. Scope, target variables & reachability

### 1.1 What the forcing source must drive

**Primary target: the fast biophysics loop.** `run_fast_biophysics` (§`meds_fast_loop`) advances each
patch by `n_fast_per_slow = nint(dt_slow/dt_fast)` operator-split sweeps of `column_fast_step`. With
`dt_slow = 1 d` and `dt_fast = 900 s` that is **96 sub-steps spanning one day** — so the *diurnal
cycle itself lives inside the fast sub-step loop*, and the forcing must be re-evaluated **per
sub-step**, not held constant per slow day (the single most important consequence for the wiring,
§6). The kernels that consume forcing, and the field/type each reads, are (consolidated from the
consumption-seam brief):

| Forcing quantity | Consumed by (kernel) | Field / type | Units |
|---|---|---|---|
| **Air temperature** `Tair` | aerodynamics | `aero_env_t%theta_atm` | [K] air temp at `zref` (θ≈T near surface; no Poisson at P0, §6.2) |
| | CAS energy (`canopy_air_update`) | `enthalpy_atm` (folded via `cas_enthalpy_of_temp`) | [J/kg] |
| | (driver) | `fast_context_t%air_temp` | [K] |
| | phenology (slow) | `pheno_env_t%temp_day` | [K] daily mean |
| **Specific humidity** `Qair` | aerodynamics | `aero_env_t%shv_atm` | [kg/kg] |
| | CAS energy (humidity twin) | `enthalpy_atm` + `w_flux_ac` | [J/kg],[kg/m²/s] |
| | soil evap (`ground_evaporation`) | `chydro_forcing_t%q_air` (= CAS shv) | [kg/kg] |
| | (driver) | `column_forcing_t%shv_atm`, `fast_context_t%shv_atm` | [kg/kg] |
| **Pressure** `Psurf` | leaf gas exchange | `leaf_env_t%pressure` | [Pa] |
| | aerodynamics | `aero_env_t%press` | [Pa] |
| | veg/ground energy | `leaf_energy_env_t%press` | [Pa] |
| | (driver) | `fast_context_t%press` | [Pa] |
| **Precipitation** `Rainf` | soil hydrology | `chydro_forcing_t%precip_ground` | [kg/m²/s] post-interception |
| | per-cohort interception | `intercept_canopy_layer(rain_above)` | [kg/m²/s] |
| | (driver) | `column_forcing_t%precip`, `fast_context_t%precip` | [kg/m²/s] |
| **Shortwave** (direct+diffuse, PAR+NIR) | canopy radiation | `rad_forcing_t%incid_beam(band)`, `%incid_diff(band)` | [W/m²] canopy top |
| | | `rad_forcing_t%cosz` | [–] (derived, `solar_cosz`) |
| | leaf gas exchange | `leaf_env_t%par` (= abs_sw/lai·par_per_w) | [µmol/m²/s] |
| | veg energy | `leaf_energy_env_t%abs_sw` | [W/m²] |
| | (driver) | `column_forcing_t%abs_sw(:)`, `%abs_sw_ground`; `fast_context_t%rad_sw_top`, `%rad_sw_ground` | [W/m²] |
| **Longwave down** `LWdown` | canopy radiation | `rad_forcing_t%incid_diff(RAD_LW)` | [W/m²] |
| | veg energy | `leaf_energy_env_t%abs_lw` | [W/m²] net |
| | (driver) | `column_forcing_t%abs_lw(:)`, `%abs_lw_ground` | [W/m²] |
| **Wind** `Wind` | aerodynamics | `aero_env_t%u_ref`, `%zref` | [m/s], [m] |
| | (driver) | `fast_context_t%u_ref`, `%zref` | [m/s], [m] |
| **CO₂** `CO2air` | leaf gas exchange | `leaf_env_t%ca` (reads **prognostic** `bio%cas%can_co2`) | [µmol/mol] |
| | aerodynamics (`cstar`) | `aero_env_t%co2_atm` | [µmol/mol] free-atm |
| | CAS CO₂ (`column_co2_step`) | arg `co2_atm` | [µmol/mol] free-atm |
| | biogeochem fallback | `co2_opts_t%co2_atm_ref` | [µmol/mol] |
| | (driver) | `column_forcing_t%co2_atm`, `fast_context_t%co2_atm` | [µmol/mol] |
| **Air density** `ρ_air` (property, derived) | aero / CAS / leaf-energy / soil-evap / CO₂ | `%rho_air` (pervasive) | [kg/m³] |

**Two seam subtleties the design must respect (concerns from the brief):**

1. **Some "forcing-looking" fields are prognostic, not atmospheric.** `leaf_env_t%ca` reads the
   *prognostic* CAS `bio%cas%can_co2`; `leaf_env_t%leaf_temp` is the *prognostic* leaf-energy
   temperature; `ustar` is an **output** of aerodynamics (`aero_out_t%ustar`) consumed by the CAS
   energy and CO₂ kernels. The met source must fill only the **free-atmosphere reference** values
   (`co2_atm`, `theta_atm`, `enthalpy_atm`, `u_ref`) and let the prognostic CAS/leaf/`ustar`
   machinery evolve — it must **not** overwrite `can_co2`, `leaf_temp`, or supply `ustar`.
2. **Air temperature never enters the energy kernels as raw [K].** It is folded into
   `enthalpy_atm [J/kg]` (via `cas_enthalpy_of_temp(air_temp, shv)`, already in `meds_thermo`) for the
   CAS budget and into `theta_atm [K]` (potential temperature) for aerodynamics. The met source
   produces `air_temp` + `shv_atm`; the *existing* `build_forcing` already computes
   `enthalpy_atm = cas_enthalpy_of_temp(air_temp, shv_atm)`. So the seam is: **the met source writes
   the raw met scalars into `fast_context_t`'s met fields; `build_forcing`/`fill_aenv` keep their
   present logic unchanged.** This is the minimal-diff wiring (§6.2).

### 1.2 The canonical variable set (ALMA + ED2 union, what the reader produces)

The reader's job is to fill one instantaneous per-site record. The **canonical 12-field set** (ALMA
core + ED2's 4-stream SW split + derived `cosz`/`ρ_air`), with the interpolation policy each carries:

| Canonical field | Symbol | Units | Interp policy (§5) |
|---|---|---|---|
| `tair_k` | Tair | [K] | linear |
| `qair` | Qair | [kg/kg] specific humidity | linear |
| `psurf_pa` | Psurf | [Pa] | linear |
| `rainf` | Rainf | [kg/m²/s] | **step-constant** (conserve accumulation) |
| `wind` | Wind | [m/s] | **energy-form** (√ of variance-weighted squares, §5.3) |
| `lwdown` | LWdown | [W/m²] | linear |
| `par_beam`, `par_diffuse` | — | [W/m²] | **cosz-weighted** |
| `nir_beam`, `nir_diffuse` | — | [W/m²] | **cosz-weighted** |
| `co2` | CO2air | [µmol/mol] | linear (or config constant) |
| `cosz` | — | [–] | **recomputed** (`solar_cosz`), never stored |
| `rho_air` | ρ_air | [kg/m³] | **recomputed** from `Tair,Psurf,Qair` |

`SWdown = par_beam + par_diffuse + nir_beam + nir_diffuse` and
`rshort_diffuse = par_diffuse + nir_diffuse` recover the aggregate ED2 diagnostics. Snow (`Snowf`) is
folded into `Rainf` with an air-temperature phase split at ingest (§5.4); the tropical BCI record is
rain-only, so P0 treats all precip as liquid.

### 1.3 Reachability facts (from the ED2 & MEDS references)

- **ED2 stores forcing in derived types passed by argument**, never global scalars: raw time-buffered
  `met_driv_data` (14 arrays) → instantaneous `met_driv_state` (25 scalars) → per-site `cpoly%met`.
  MEDS mirrors this two-tier split (`met_driver_t` buffer → `met_forcing_t` instantaneous, §3).
- **ED2 disaggregates each SW stream by a mean daytime zenith factor** (`mean_daysecz`), forming a
  beam-perpendicular flux and rescaling to the instantaneous geometry; night (`cosz ≤ cosz_min`) is
  forced to zero. **MEDS deliberately does NOT port ED2's `⟨sec z⟩ = ⟨1/cosz⟩`** — that quantity does
  *not* conserve the interval-mean flux (by Jensen `⟨1/cosz⟩ ≥ 1/⟨cosz⟩`, badly so near
  sunrise/sunset). MEDS instead uses the **reciprocal of the window-mean cosine**, `1/⟨cosz⟩`, which is
  the unique factor that makes the reconstructed flux integrate back to `F_avg` (§5.1, the exact
  identity test 1d checks). MEDS already has `solar_cosz(t, t_sec, latitude_deg)` in `meds_time` — the
  analogue of `ed_zen`. It needs a **site latitude**, absent from the TOML today: adding `[site]`
  (§6.4) is a hard prerequisite for any radiation forcing.
- **ED2 interpolates precipitation step-constant** and wind as a squared (energy) quantity. MEDS
  ports both (§5.3–§5.4).
- **The `fast_context_t` comment is explicit**: *"the CALLER builds this from TOML in the production
  path; the MVP holds constant, horizontally-uniform boundary conditions."* Today only
  `test/test_fast_loop.f90` builds it. This design supplies the production caller.
- **CRITICAL WIRING GAP (unchanged by this doc, closed by §6):** `meds_main` line 124 calls
  `advance_one_step(site, cfg, is_new_month, is_new_year)` **without** the optional `fast_ctx`, and
  never calls `init_fast_reservoirs`, so `run_fast_biophysics` never fires in production even with
  `fast_biophysics_on = true`. Wiring the met source is done **together with** passing `fast_ctx` and
  seeding the reservoirs.

---

## 2. Where it lives — the new `src/forcing/` library (comment 1)

The DAG is `shared ← {allometry, plant} ← state ← demography ← aux ← main`, with `biophysics`/
`biogeochemistry` **stateless siblings** linking `shared` only. Meteorological forcing gets its **own
sibling library `libmeds_forcing` in `src/forcing/`** — the single home for every prescribed external
driver. All three forcing modules (runtime types, kernels, reader) live together; only the plain-scalar
config stays in `shared`:

```
shared ──┬─ meds_forcing_config  (in src/shared) — forcing_config_t: the [forcing]/[site] PLAIN-SCALAR
         │     config (path, format, dt, lat/lon, CO2, sw_partition, ...). Lives in shared so meds_config
         │     (the DAG ROOT) can carry it with NO backward shared→forcing edge (the soil_thermal_params_t
         │     rule from the energy design). Pure data; no methods.
         │
         ├─ meds_netcdf_c  (PROMOTED to its own low-level target — see below) — the netCDF ISO-C bindings,
         │     shared by libmeds_io (state serialization) AND libmeds_forcing (forcing read). IO-gated.
         ↓
   forcing ── libmeds_forcing  (NEW, src/forcing) — links meds_shared + meds_netcdf_c ONLY (NOT demography):
         ├─ meds_forcing_types    — met_forcing_t (instantaneous per-site record), met_record_t, and the
         │      MUTABLE reader buffer met_driver_t. Pure data; per-site scalars + 4 SW-band reals + a
         │      grid_index; GPU-mappable, no cohort-sized allocatables.
         ├─ meds_forcing_kernels  — PURE/ELEMENTAL disaggregate_shortwave, partition_shortwave,
         │      interpolate_forcing, cosz_reconstruct_factor, rh_to_specific_humidity, dewpoint_to_specific_
         │      humidity, precip_phase, synthesize_lwdown. Depends ONLY on meds_thermo + meds_time (shared);
         │      REUSES meds_thermo%air_density and %sat_vapor_pressure — does NOT re-invent air_density/esat.
         └─ meds_met_driver       — the READER: met_open/met_advance/met_instant/met_close over the MEDS
                forcing NetCDF (multi-grid, §7) + the no-file CONST backend; owns the only mutable forcing
                state (met_driver_t). Called by the driver (meds_fast_loop / meds_main); NEVER named by a
                physics library.
```

**Why one `src/forcing/` library — and why it is *simpler* than the old split.** The previous draft
scattered forcing across `src/shared` (types+kernels) and `src/io` (reader), *purely* to keep the reader —
then living in `src/io` — from forcing a `libmeds_io → libmeds_biophysics` edge when it called the physics
kernels. **That contortion is now unnecessary:** with the reader in its own `libmeds_forcing`, the kernels
sit **beside** it and both link `meds_shared` directly. The single library:

- **Owns the runtime forcing state cohesively.** `met_forcing_t` (the read-only boundary-condition value —
  to the fast loop what `rad_forcing_t` / `chydro_forcing_t` are to their kernels), `met_record_t` (one raw
  timestamped record), and the mutable `met_driver_t` buffer (ED2 `cgrid%metinput` analogue) all live in
  `meds_forcing_types`. No physics library names any of them — the acyclic DAG is preserved.
- **Reads NetCDF without dragging in demography.** The reader needs the netCDF C bindings but must **not**
  inherit `libmeds_io`'s `→ meds_demography` edge (io links demography only for *state serialization*, which
  forcing has no business touching). So **`meds_netcdf_c.f90` is promoted out of `libmeds_io` into its own
  minimal target** (links `netCDF::netcdf` + `meds_shared` for kinds); `libmeds_io` and `libmeds_forcing`
  both link it. This is a small, clean refactor (move one file's target membership) that *removes* an
  accidental dependency rather than adding one.
- **Keeps config in `shared`.** `forcing_config_t` is plain scalars on `meds_config_t`; putting the *type*
  in `shared` (`meds_forcing_config.f90`, or folded into `meds_config`) lets the DAG root carry it with no
  `shared → forcing` back-edge — exactly the `soil_thermal_params_t` precedent.
- **Reuses shared math.** `meds_thermo` already provides `elemental air_density(t_k,p_pa,shv)` and the
  Bolton-1980 `sat_vapor_pressure`/`sat_specific_humidity`; the reader calls them directly. Genuinely new
  math is only the SW disaggregation/partition, the RH- and **dewpoint-**to-specific-humidity conversions,
  and `precip_phase`.
- **Is the future home for other forcing.** Disturbance-event / fire / land-use schedules and prescribed
  CO₂ / N-deposition streams are the same *kind* of object — a time-indexed external driver read from a file
  — and will add `meds_*_driver` modules here without touching the physics libraries.

**netCDF + IO are ALWAYS compiled — no stubs (revision comment 1).** The `-DMEDS_ENABLE_IO=OFF`
netCDF-free build and its no-op `meds_io_stub.f90` are **removed**: netCDF is now a hard build dependency
for every configuration (`find_package(netCDF CONFIG REQUIRED)` unconditional; `enable_language(C)`
unconditional; `meds_io` always the real writer). Consequently the forcing library needs **no
`meds_met_driver_stub` and no `meds_netcdf_c_stub`** — `libmeds_forcing` always links the real
`meds_netcdf_c`, and `meds_aux` links `libmeds_forcing` unconditionally. The **CONST backend stays** — but
it is a *legitimate reader mode* (a no-file reference-climate box for tests / a run with no forcing file),
**not** a netCDF-absence fallback. `meds_netcdf_c` is still promoted to its own target so `libmeds_forcing`
gets the netCDF bindings without inheriting `libmeds_io`'s `→ meds_demography` (state-serialization) edge —
that DAG-cleanliness reason is independent of the (now-deleted) on/off gate. This build simplification is a
small, self-contained change to `CMakeLists.txt` (delete the `option`/`if`-`else` gate, delete
`meds_io_stub.f90`) that lands with the forcing library.

### 2.1 Files & CMake wiring

| File | Role | Analogue |
|---|---|---|
| `src/forcing/meds_forcing_types.f90` (new) | `met_forcing_t` (instantaneous record), `met_record_t`, `met_driver_t` (mutable buffer, incl. `grid_index`); `INTERP_*`/`SWPART_*`/`METAVG_*`/`MET_BACKEND_*` selector codes | `meds_biophysics_types` type block |
| `src/forcing/meds_forcing_kernels.f90` (new) | `pure`/`elemental` `disaggregate_shortwave`, `partition_shortwave`, `interpolate_forcing`, `cosz_reconstruct_factor`, `rh_to_specific_humidity`, `dewpoint_to_specific_humidity`, `precip_phase`, `synthesize_lwdown`. **Calls `meds_thermo%air_density` + `%sat_vapor_pressure` — no re-invented `air_density`/`esat`.** | `meds_thermo` |
| `src/forcing/meds_met_driver.f90` (new) | `met_open`, `met_advance`, `met_instant`, `met_close`; **NetCDF backend (multi-grid `(time,grid)`, §7) + CONST backend**; per-polygon `grid_index` selection. Links `meds_shared` + `meds_netcdf_c` | `meds_config_io` |
| `src/shared/meds_forcing_config.f90` (new) OR fold into `meds_config` | `forcing_config_t` plain-scalar `[forcing]`/`[site]` config on `meds_config_t` | `soil_thermal_params_t` rule |
| `src/shared/meds_time.f90` (**extend — prerequisite**) | `seconds_into_day(now)`, `time_advance_seconds(now, dsec)`, `seconds_between(a, b)` (second-accurate JDN + sec-of-day across day boundaries); `sec_of_day` is currently private (line 164) — the hourly→900 s interpolation + window slide need these | its `time_advance_days`/`days_between` block |
| `src/io/meds_netcdf_c.f90` (**move to own target**) + extend | promote out of `libmeds_io` into a minimal `meds_netcdf_c` library both `io`+`forcing` link (always built — netCDF is mandatory now); add `nc_inq_dimid`+`nc_inq_dimlen` bindings (the reader needs the `time`/`grid` dimension lengths) | — |
| `CMakeLists.txt` (**edit — comment 1**) | delete the `MEDS_ENABLE_IO` option + its `if`/`else` gate + `meds_io_stub`; `find_package(netCDF … REQUIRED)` + `enable_language(C)` unconditional; add `add_library(meds_forcing …)`; `meds_aux` links `meds_forcing` | — |
| `src/io/meds_io_stub.f90` (**delete — comment 1**) | removed; netCDF is always compiled | — |
| `src/io/meds_config_io.f90` (extend) | `[forcing]` + `[site]` key loaders (`req_*`), incl. `grid_index`, `sw_partition`, `latitude`/`longitude` | its `[fast]` block |
| `src/driver/meds_fast_loop.f90` (edit) | per-sub-step met refresh in `run_fast_biophysics`; `build_forcing`/`fill_aenv` read a `met_forcing_t` | — |
| `src/driver/meds_main.f90` (edit) | open the driver, seed reservoirs, thread `[site].reference_height` into `ctx%zref`, pass `fast_ctx` | — |
| `scripts/download_era5land.py` (new) | CDS-API download of ERA5-Land for a requested lat/lon (or box) → raw ERA5-Land NetCDF (§7) | — |
| `scripts/prep_era5land_forcing.py` (new) | raw ERA5-Land NetCDF → **MEDS multi-grid forcing NetCDF** (de-accumulate, unit-convert, humidity from dewpoint) (§7) | — |
| `test/test_met_driver.f90` (new) | CTest: read multi-grid NetCDF, interpolate, disaggregate, diurnal cycle | `test_fast_loop` |

CMake needs a **new `add_library(meds_forcing …)`** globbing `src/forcing/*.f90` (linking `meds_shared`
+ `meds_netcdf_c`), the **promotion** of `meds_netcdf_c` to its own always-built target, the **removal** of
the `MEDS_ENABLE_IO` gate + `meds_io_stub` (comment 1), and `meds_aux` gains `meds_forcing` on its link
line. Append the `test_met_driver` executable after
the `test_fast_loop` block. **Build nvfortran multicore on every new module** (a green ifx suite is not
sufficient — CLAUDE.md portability trap #7; the SW-partition kernels return arrays and must not be passed
straight into a call — bind to a named array first).

---

## 3. Derived types & state model

Two tiers (ED2's `met_driv_data` / `met_driv_state` split), plus the reader's config.

### 3.1 `met_forcing_t` — the instantaneous per-site record (`meds_forcing_types`, `src/forcing`)

```fortran
! SW band indices ALIAS the existing RT vocabulary (meds_biophysics_types: RAD_VIS=1, RAD_NIR=2,
! RAD_LW=3) so the forcing bands and the two-stream RT bands share ONE set of names at the seam:
integer(ik), parameter :: RAD_PAR = RAD_VIS   ! = 1_ik ; PAR ≡ the RT "visible" band (do NOT mint a new code)

type :: met_forcing_t                                        ! per-SITE instantaneous atmospheric state
   real(wp) :: tair_k       = 288.0_wp    ! [K]        air temperature at reference height
   real(wp) :: qair         = 0.008_wp    ! [kg/kg]    specific humidity
   real(wp) :: psurf_pa     = 101325.0_wp ! [Pa]       surface pressure
   real(wp) :: rainf        = 0.0_wp       ! [kg/m2/s] liquid precipitation rate (post phase-split)
   real(wp) :: snowf        = 0.0_wp       ! [kg/m2/s] frozen precip (0 for tropical BCI)
   real(wp) :: wind         = 2.0_wp       ! [m/s]     wind speed at reference height
   real(wp) :: lwdown       = 380.0_wp     ! [W/m2]    downwelling longwave (positive down)
   ! reference-climate SW: the four streams SUM to 400 W/m2 to reproduce fast_context_t%rad_sw_top=400
   ! (NOT 0 — a zeroed default gives zero canopy light and BREAKS the CONST no-op; see §6.2):
   real(wp) :: par_beam     = 180.0_wp     ! [W/m2]    direct-beam PAR at canopy top
   real(wp) :: par_diffuse  = 40.0_wp      ! [W/m2]    diffuse PAR
   real(wp) :: nir_beam     = 150.0_wp     ! [W/m2]    direct-beam NIR
   real(wp) :: nir_diffuse  = 30.0_wp      ! [W/m2]    diffuse NIR   (Σ = 400 W/m2)
   real(wp) :: co2          = 400.0_wp     ! [umol/mol] free-atmosphere CO2
   real(wp) :: cosz         = 0.0_wp       ! [-]       cos(solar zenith); DERIVED each substep, not stored
   real(wp) :: rho_air      = 1.2_wp       ! [kg/m3]   DERIVED from tair/psurf/qair
contains
   ! swdown() = par_beam+par_diffuse+nir_beam+nir_diffuse ; rshort_diffuse() = par_diffuse+nir_diffuse
end type met_forcing_t
```

All defaults are valid (clean under `-fpe0`), so `met_forcing_t()` is a usable "reference climate" box:
its four SW streams **sum to 400 W/m² = today's `fast_context_t%rad_sw_top`** (`meds_fast_loop.f90:41`;
`air_temp=288`, `shv_atm=0.008`, `co2=400`, `press=101325`, `u_ref=2.0`, `rho_air=1.2` likewise match).
The **ground SW default is not carried here** — today's ad-hoc `rad_sw_ground=60` (`:42`) is an
*independent* constant with no canopy-top relation; under CONST it is reproduced by threading a
reference `rad_sw_ground` (below) rather than by Beer's law, so the reference climate is an honest
reproduction of the current constants. See §6.2 for why this makes CONST a **near-no-op** (canopy-top
SW exact) rather than the bit-for-bit claim the earlier draft asserted.

### 3.2 `met_driver_t` — the reader buffer (MUTABLE, `meds_met_driver`, `src/forcing`)

The **only** mutable forcing state, **per polygon/site**. Holds the two bracketing records read from the
file at this site's `grid_index` plus the read cursor; it is the ED2 `cgrid%metinput` analogue, owned by
the driver and threaded like `meds_io_t`. For a multi-polygon run there is one `met_driver_t` per polygon,
each bound to its own `grid_index` (all reading the same file, different `grid` slices).

```fortran
type :: met_record_t                       ! one raw timestamped record as read (pre-interpolation)
   type(meds_time_t) :: when               ! timestamp of this record (interval label per METAVG_*)
   real(wp) :: tair_k, qair, psurf_pa, rainf, wind, lwdown, co2
   real(wp) :: par_beam, par_diffuse, nir_beam, nir_diffuse   ! 4-stream: split at INGEST from total SWdown
end type met_record_t                      !   (partition_shortwave, §5.6) using the record's interval-mean
                                           !   cosz — so a total-SW file (ERA5-Land) and a pre-split file
                                           !   (a future 4-stream tower) both fill these four the same way.

type :: met_driver_t                        ! MUTABLE per-POLYGON reader state
   type(forcing_config_t) :: fcfg           ! [forcing]/[site] config (path, format, dt, lat/lon, CO2, ...)
   integer(ik) :: backend    = MET_BACKEND_CONST   ! NETCDF | CONST(no file)   (CSV backend retired)
   integer(ik) :: ncid       = -1_ik        ! NetCDF handle (via meds_netcdf_c), -1 if const
   integer(ik) :: grid_index = 1_ik         ! which location (1..ngrid) on the (time,grid) file this site reads
   integer(ik) :: ngrid      = 1_ik         ! total locations in the file (from the `grid` dimension)
   integer(ik) :: nrec       = 0_ik         ! total time records in file (for cycling)
   integer(ik) :: irec_prev  = 0_ik         ! cursor: index of the "previous" bracketing record
   real(wp)    :: dt_forcing  = 3600.0_wp   ! [s] native forcing interval (hourly for ERA5-Land)
   type(met_record_t) :: rec_prev, rec_next ! the two records bracketing the model time (at grid_index)
   type(meds_time_t)  :: base_time          ! file time #1 timestamp (calendar anchor)
end type met_driver_t
```

### 3.3 `forcing_config_t` — the [forcing]/[site] block (`meds_forcing_config`, `src/shared`, plain scalars)

Lives on `meds_config_t` (a plain-scalar component — `meds_config` is the DAG root and must not
`use` a `forcing`/biophysics type, exactly the `soil_thermal_params_t` rule from the energy design).
Defaults are the **Ithaca NY / ERA5-Land** reference site:

```fortran
type :: forcing_config_t
   logical            :: forcing_on   = .false.    ! master gate (independent of fast_biophysics_on)
   integer(ik)        :: format       = MET_BACKEND_NETCDF  ! "netcdf" | "const"  (CSV retired)
   character(len=256) :: path         = ''         ! forcing NetCDF path
   integer(ik)        :: grid_index   = 1_ik        ! which (time,grid) location this polygon reads (§7)
   real(wp)           :: dt_forcing   = 3600.0_wp   ! [s] native interval (hourly ERA5-Land)
   integer(ik)        :: avg_convention = METAVG_END ! flux vars accumulate over the hour ENDING at the stamp
   integer(ik)        :: sw_partition = SWPART_WEISS_NORMAN ! ERA5-Land carries TOTAL SW -> partition required
   integer(ik)        :: lwdown_source = LW_FILE            ! file (ERA5-Land has strd) | synthesize
   ! (no gap_policy: MEDS never gap-fills — a missing required field is a hard error, §5.5)
   real(wp)           :: co2_const    = 420.0_wp   ! [umol/mol] ERA5-Land has NO CO2 -> single config authority (§5.2)
   real(wp)           :: rad_sw_ground_const = 60.0_wp ! [W/m2] CONST-backend ground SW (reproduces fast_context_t:42)
   logical            :: recycle      = .true.     ! cycle the record when the run outruns the file
   integer(ik)        :: start_clamp  = CLAMP_ERROR ! model start < base_time: error | hold record #1 (§4.2)
   ! ---- site geolocation ([site]) — prerequisites for solar geometry & lapse ----
   real(wp)           :: latitude_deg  = 42.44_wp  ! [deg +N]  Ithaca NY ~42.44 N
   real(wp)           :: longitude_deg = -76.50_wp ! [deg +E]  Ithaca NY ~-76.50 E (solar time, §5.1)
   real(wp)           :: utc_offset_h  = 0.0_wp    ! [h]       ERA5-Land is UTC -> file clock IS UTC (offset 0)
   logical            :: apply_solar_longitude = .true. ! UTC file -> longitude gives local solar time (§5.1)
   real(wp)           :: elevation_m   = 320.0_wp  ! [m]       Ithaca ~320 m; reserved for lapse (P2)
   real(wp)           :: reference_height = 40.0_wp! [m]       forcing reference height; must exceed max hgt_max
   real(wp)           :: wind_meas_height = 10.0_wp! [m]       ERA5-Land wind is 10 m (below tall canopy — §5.2/§10)
end type forcing_config_t
```

**Note the reversed solar-time default vs BCI.** Because ERA5-Land is **UTC-referenced**, the principled
path is `apply_solar_longitude = .true.` with `utc_offset_h = 0`: `solar_cosz` gets local apparent solar
time from the longitude term directly (§5.1), with no local-clock ambiguity. This is *cleaner* than the old
BCI local-standard-time default, which had a standing ~20–35 min phase error.

### 3.4 State-model summary

No global mutable forcing state. The mutable buffer `met_driver_t` is owned by the driver
(`meds_main` constructs it, threads it read-write into a refresh call each fast sub-step). The
instantaneous `met_forcing_t` is a read-only value produced by the reader and consumed by the fast
loop; the disaggregation/partition kernels are pure (no saved state); `cosz`/`ρ_air` are recomputed,
never stored. This keeps the physics libraries free of hidden met state and GPU-offload-safe — the
identical discipline as `rad_forcing_t` / `energy_forcing_t`.

---

## 4. The reader module — `src/forcing/meds_met_driver.f90`

Public API (mirrors `meds_io`'s open/write/close cadence and `ED2`'s `init/read/update` split):

```fortran
public :: met_driver_t, met_forcing_t
public :: met_open      ! (drv, fcfg)                 open file, read header, load records #1..2
public :: met_advance   ! (drv, now)                  slide the bracketing window to contain `now`
public :: met_instant   ! (drv, now, t_sec, lat) -> met_forcing_t   interpolate+disaggregate to `now`
public :: met_close     ! (drv)                       close handle / free
```

### 4.1 `met_open(drv, fcfg)` — header & first records

```
select case (fcfg%format)
case (MET_BACKEND_CONST)   ! no file: drv%rec_prev = rec_next = reference climate (met_forcing_t defaults)
case (MET_BACKEND_NETCDF)  ! nc_open via meds_netcdf_c; read the `grid` dim -> ngrid; validate
                           !   1 <= fcfg%grid_index <= ngrid; read time[] (nrec = len(time)); base_time from
                           !   the time-units attribute; load records #1..2 at THIS grid_index (§7 layout)
end select
drv%grid_index = fcfg%grid_index
drv%irec_prev  = 1
```

The NetCDF backend reuses the `meds_netcdf_c` ISO-C bindings (`nc_open`, `nc_inq_varid`,
`nc_get_vara_double`) — no netCDF-Fortran dependency (unavailable for ifx/nvfortran here). **It adds
`nc_inq_dimid` + `nc_inq_dimlen`** (§2.1): the reader needs both `nrec` = `len(time)` (cycling/EOF) **and**
`ngrid` = `len(grid)` (multi-grid validation); `meds_netcdf_c` binds no `nc_inq_dim*` today. **Multi-grid
read (comment 3).** The file is `(time, grid)` (§7); each variable is read as the **hyperslab at the fixed
`grid_index`** — `nc_get_vara_double(varid, start=[irec, grid_index], count=[2, 1])` pulls the two
bracketing time steps for this one location, so a multi-polygon run issues one such slabbed read per
polygon's `met_driver_t`, never loading the whole grid. `grid_index` comes from `[forcing].grid_index`
(P0: an explicit index; P1/P2: resolved by nearest-location match of the polygon's lat/lon against
`latitude(grid)`/`longitude(grid)`, ED2 `match_poly_grid`). A single-site file is simply `ngrid = 1`. The
canonical NetCDF layout (what the formatter writes, §7) is dims `time` (unlimited) × `grid`, coordinate
variables `time` (`units = "seconds since <base>"`), `latitude(grid)`, `longitude(grid)`, and the forcing
variables shaped `(time, grid)`, ALMA-named (`Tair`, `Qair`, `PSurf`, `Rainf`, `Wind`, `LWdown`, and either
`SWdown` total or the four `SWdown_*_beam/diffuse` streams — §7).

### 4.2 `met_advance(drv, now)` — window slide (calendar alignment via `meds_time`)

The bracketing records must satisfy `rec_prev%when ≤ now < rec_next%when`. Advance is a
`transfer`-style slide (ED2 `transfer_ol_month` analogue), reading forward until the window contains
`now`. **Start-before-`base_time` (spin-up) is handled explicitly** — BCI begins 2012-07-03, and a
plausible spin-up config sets `start_time` earlier, so `now < rec_prev%when` (record #1) is a real
case the earlier draft left unspecified (it only handled running off the *end*):

```
if ( time_lt(now, drv%base_time) ) then           ! model time precedes the first forcing record
   select case (drv%fcfg%start_clamp)
   case (CLAMP_HOLD)                                ! clamp to record #1 (persistence), warn once
      call warn_once('met_driver: start_time < base_time; holding record #1 until forcing begins')
      return                                        ! rec_prev = rec_next = record #1; w_next forced to 0 (§4.4)
   case default                                     ! CLAMP_ERROR (default): fail loud, symmetric with EOF
      error stop 'met_driver: model start precedes first forcing record and start_clamp=error'
   end select
end if
do while ( time_lt(drv%rec_next%when, now) .or. time_eq(drv%rec_next%when, now) )
   drv%rec_prev = drv%rec_next
   drv%irec_prev = drv%irec_prev + 1
   if (drv%irec_prev + 1 > drv%nrec) then          ! ran off the end
      if (drv%fcfg%recycle) then                    ! wrap to the first record (spin-up cycling)
         call rewind_to_first(drv) ; cycle
      else
         error stop 'met_driver: model time past end of forcing and recycle=off'
      end if
   end if
   call read_record(drv, drv%irec_prev + 1, drv%rec_next)
end do
```

Calendar alignment uses `meds_time` throughout: record timestamps are `meds_time_t`; `time_lt`/
`time_eq` bound the window; leap years are exact. **Feb-29 reconciliation** (ED2's leap quirk) is a P1
concern only when cycling a non-leap forcing year onto a leap model year — P0/BCI reads the record's
own real dates so there is no mismatch.

### 4.3 File backend: one NetCDF format, no CSV (comment 3)

The **CSV backend is retired.** The reader consumes exactly one file format — the **MEDS forcing NetCDF**
(§7) — plus the no-file CONST backend. The rationale that justified a CSV path ("drive BCI without a prep
step") is gone: the source is now ERA5-Land, which the CDS delivers as NetCDF, and the formatter script
(§7) always produces the MEDS NetCDF. Any future tower record (BCI, FLUXNET) re-enters as **another
formatter script emitting the same NetCDF** — so the Fortran reader is source-agnostic and never changes
when the raw data source does. All SI unit conversion, de-accumulation, humidity-from-dewpoint, and gap
handling move **into the formatter** (§7); the reader ingests already-canonical SI values (its only
per-read math is the SW partition of a total-SW file, §5.6, and the temporal interpolation, §5).

**One canonical saturation formula, shared reader ↔ formatter.** Any place the reader still computes
humidity — `rh_to_specific_humidity(RH,T,P)` and `dewpoint_to_specific_humidity(Td,P)` — is built on
`meds_thermo%sat_vapor_pressure` (the reader carries **no** private `esat`). The formatter script (§7) that
computes `Qair` from ERA5-Land dewpoint **must use the identical saturation formula and constants** (§7),
or the same `(Td, P)` yields a different `qair` on the two sides, breaking the Python↔Fortran agreement
test (§9.5). The formatter is written against the same `meds_thermo` form.

### 4.4 `met_instant(drv, now, t_sec, lat) -> met_forcing_t`

The hot path (ED2 `update_met_drivers` analogue). Computes the interpolation weight from the
bracketing timestamps, applies the per-variable policy (§5), recomputes `cosz` and `ρ_air`, and
returns the instantaneous record:

`seconds_between` / `time_advance_seconds` / `seconds_into_day` are the **new second-accurate
`meds_time` helpers** (§2.1 — prerequisites, not incidental), and `air_density` is `meds_thermo`'s
existing virtual-temperature routine, not a re-implementation:

```
dt_win = real(seconds_between(rec_prev%when, rec_next%when), wp)   ! = dt_forcing  (NEW meds_time helper)
elapsed = real(seconds_between(rec_prev%when, now), wp)
w_next  = clamp(elapsed / max(dt_win, tiny_num), 0, 1) ; w_prev = 1 - w_next
cosz    = solar_cosz(now, t_sec, lat)                              ! meds_time, floored >= 0

met%tair_k  = interpolate_forcing(INTERP_LINEAR, prev%tair_k, next%tair_k, w_next)
met%qair    = interpolate_forcing(INTERP_LINEAR, prev%qair,   next%qair,   w_next)
met%psurf_pa= interpolate_forcing(INTERP_LINEAR, prev%psurf_pa,next%psurf_pa,w_next)
met%wind    = interpolate_forcing(INTERP_LINEAR, prev%wind,   next%wind,   w_next)   ! §5.3 energy form
met%lwdown  = interpolate_forcing(INTERP_LINEAR, prev%lwdown, next%lwdown, w_next)
met%co2     = merge(fcfg%co2_const, interpolate_forcing(INTERP_LINEAR,prev%co2,next%co2,w_next), no_co2_in_file)
met%rainf   = interpolate_forcing(INTERP_STEP,   prev%rainf,  next%rainf,  w_next)   ! step-constant
call disaggregate_sw_bands(prev, next, cosz, cosz_bar, w_prev, w_next, met)          ! §5.1 (reciprocal-mean-cosz)
met%rho_air = air_density(met%tair_k, met%psurf_pa, met%qair)   ! REUSES meds_thermo%air_density (shared)
```

---

## 5. Temporal interpolation & disaggregation (hourly → 900 s)

The native ERA5-Land interval is **3600 s (hourly)**; `dt_fast = 900 s`. Each hourly record is split
into four 15-min fast sub-steps (`dt_forcing/dt_fast = 4`; the machinery handles any ratio). Per-variable
policy (community
practice; concerns from the literature brief):

### 5.1 Shortwave — solar-zenith-weighted (the load-bearing method)

Linear interpolation of SW smears the diurnal cycle and is badly wrong in the tropics where `cosz`
changes fast even at 30-min resolution. MEDS uses a **strictly interval-mean-conserving cosine-zenith
disaggregation**. The requirement is exact: reconstruct an instantaneous flux `F(t)` from the
interval-mean `F_avg` over `[rec_prev%when, rec_next%when]` such that its time-average over the SAME
window returns `F_avg`. If `F(t) = f_perp · cosz(t)` on the daytime part of the window (0 at night),
then

```
⟨F⟩_win = f_perp · (1/T_win) ∫_daytime cosz dt = f_perp · ⟨cosz⟩_win   ⟹   f_perp = F_avg / ⟨cosz⟩_win
F(now)  = F_avg · cosz(now) / ⟨cosz⟩_win        ! 0 at night (cosz(now) <= cosz_min)
```

where **`⟨cosz⟩_win` is the window-mean of `cosz` with night sub-samples clamped to 0** — NOT the mean
secant. This is the load-bearing correction over the earlier draft: ED2's `mean_daysecz` computes the
**mean of the secant `⟨sec z⟩ = ⟨1/cosz⟩`**, and by Jensen's inequality `⟨1/cosz⟩ ≥ 1/⟨cosz⟩` — the two
disagree, badly, near sunrise/sunset where `cosz → 0`. Only the **reciprocal of the mean cosine**
`1/⟨cosz⟩` conserves the interval mean (the identity test 1d checks); `⟨sec z⟩` reconstructs an SW that
does *not* integrate back to `F_avg` and is **biased high at low sun**. MEDS therefore deliberately
does **not** port `mean_daysecz`; the kernel is named `cosz_reconstruct_factor` (returning
`1/⟨cosz⟩_win`), not `mean_daytime_secant`.

```
cosz_bar = cosz_reconstruct_factor(lat, win_start, dt_sub, dt_win)   ! = 1 / ⟨cosz⟩_win ; ⟨⟩ over ALL subsamples
F(now)   = merge(F_avg * cosz(now) * cosz_bar, 0.0_wp, cosz_bar_valid .and. cosz(now) > cosz_min)
```

applied **independently to each of the four SW streams** (`par_beam, par_diffuse, nir_beam,
nir_diffuse`). `cosz_reconstruct_factor` subsamples `solar_cosz` over the window at `dt_sub` (midpoint
rule), forming `⟨cosz⟩_win = (1/N) Σ_i max(cosz_i, 0)` — night sub-samples contribute exactly 0, so the
normaliser is the mean over the **full window** (this is what makes the average close over a
sunrise/sunset window, where daytime seconds < window seconds). **The divide-by-small guard is on the
aggregate `⟨cosz⟩_win`, not per-step:** `cosz_bar_valid = ⟨cosz⟩_win > cosz_bar_min` (a fully-night
window has `⟨cosz⟩_win = 0` and `F_avg = 0`, so all SW routes to 0); there is no per-step secant to blow
up. **`avg_convention` (METAVG_*)** shifts the window reference (instant / end / begin / center) —
FLUXNET half-hourly is interval-ending-labelled, reanalysis differs; getting this wrong phase-shifts
the diurnal cycle by up to one window, on top of the local-time solar-phase offset (below).

**Solar-phase / timezone — ERA5-Land is UTC, so the longitude path is the *natural* one (comment 2).**
`solar_cosz(t, t_sec, latitude_deg)` (`meds_time.f90:79`) treats `t_sec` as **local apparent solar
seconds** — `frac_day = t_sec/86400`, so `t_sec` at 0.5·86400 is solar noon — with **no longitude / no
equation-of-time (EoT) term**. ERA5-Land timestamps are **UTC**, so feeding UTC-seconds straight into
`solar_cosz` would place solar noon at 00:00 UTC everywhere — wrong by the site's whole longitude offset
(Ithaca at −76.5 °E is ~5.1 h = ~76.5·240 s west of Greenwich). With a UTC file the fix is therefore
**mandatory, not optional**, and it is *cleaner* than BCI's local-clock case (no timezone-vs-meridian
ambiguity):
- **P0 default for ERA5-Land (`apply_solar_longitude=.true.`, `utc_offset_h=0`):** pass
  `t_sec_solar = t_sec_utc + longitude_deg·240 s/deg + eot(doy)` into `solar_cosz` (240 s per degree of
  longitude; `+`E, so Ithaca's negative longitude shifts solar noon *later* in UTC, correctly). The EoT
  term (±16 min) is a cheap `doy`-based series; include it or bound it in the test plan (test 2).
- A **local-clock file** (a future tower with local-standard timestamps) instead sets
  `apply_solar_longitude` per its convention with `utc_offset_h ≠ 0`; the general transform is
  `t_sec_solar = t_sec + (longitude_deg − 15·utc_offset_h)·240 s/deg + eot(doy)`, which reduces to the UTC
  case when `utc_offset_h = 0`. One code path, two configs.

If the file supplies only *total* SWdown (no PAR/diffuse split — **the ERA5-Land case**), the
`sw_partition` kernel first splits it (§5.6) before disaggregation. ERA5-Land supplies only total `ssrd`
(§5.2), so P0 **requires** the full total→(direct/diffuse)×(PAR/NIR) partition — the kernel that was
"reserved for reanalysis" in the earlier draft is now on the P0 critical path (§5.6).

### 5.2 ERA5-Land → MEDS-field mapping table (with unit conversions)

The **formatter** (`scripts/prep_era5land_forcing.py`, §7) does every conversion below; the reader ingests
already-canonical SI. Native ERA5-Land is **hourly, UTC**. Two variable classes need care: the
**instantaneous** states (t2m, d2m, sp, u10, v10 — valid *at* the timestamp) vs the **accumulated** fluxes
(tp, ssrd, strd — accumulated from 00 UTC; **de-accumulate** to a per-hour amount, then to a rate/mean;
§7). CO₂ and any direct/diffuse or PAR/NIR split are **absent** from ERA5-Land.

| ERA5-Land var (cdsapi name) | short | units (raw) | → MEDS field | conversion | interp policy |
|---|---|---|---|---|---|
| `2m_temperature` | `t2m` | K | `tair_k` | as-is | linear |
| `2m_dewpoint_temperature` | `d2m` | K | `qair` | `dewpoint_to_specific_humidity(d2m, sp)` (§7) | linear |
| `surface_pressure` | `sp` | Pa | `psurf_pa` | as-is | linear |
| `total_precipitation` | `tp` | m (accum) | `rainf` | de-accum→ hourly `m`; `× ρ_w(1000)/3600` → kg/m²/s | **step-constant** |
| `10m_u_component_of_wind` | `u10` | m/s | `wind` | with `v10` | **energy form** (§5.3) |
| `10m_v_component_of_wind` | `v10` | m/s | `wind` | `wind = √(u10² + v10²)`, floor `u_min` | (at **10 m** — §10) |
| `surface_solar_radiation_downwards` | `ssrd` | J/m² (accum) | total `SWdown` | de-accum→ hourly `J/m²`; `/3600` → W/m²; then **partition** (§5.6) | **cosz-weighted** per stream |
| `surface_thermal_radiation_downwards` | `strd` | J/m² (accum) | `lwdown` | de-accum→ hourly `J/m²`; `/3600` → W/m² | linear |
| *(none)* | — | — | `co2` | `fcfg%co2_const` (default 420 µmol/mol) | constant |

**De-accumulation (the load-bearing ERA5-Land step, §7).** `tp`/`ssrd`/`strd` are **accumulated from 00
UTC** and reset daily. The per-hour amount over the hour ending at valid time `H` is `accum(H) −
accum(H−1)`, with the first step of each UTC day handled per the ECMWF recipe (the 01 UTC value is the
00–01 total; the 00 UTC step boundary is the daily reset). Getting this off-by-one wrong double-counts or
zeroes the first hour of every day — so the formatter implements it explicitly and the test (§9) asserts a
known daily SW integral. The formatter writes **already-de-accumulated mean fluxes** (`Rainf` [kg/m²/s],
`SWdown` [W/m²], `LWdown` [W/m²]) into the MEDS NetCDF, timestamped at the hour end (`avg_convention=end`).

**Shortwave: no split in ERA5-Land → the partition kernel is P0-required (§5.6).** ERA5-Land gives only
total `ssrd`. The formatter (or the reader) splits total SWdown into **(direct beam, diffuse) × (PAR,
NIR)** via the clearness-index / Weiss–Norman method (§5.6), which needs the cosine of the solar zenith and
the top-of-atmosphere irradiance — both available from `solar_cosz` + the site lat/lon. **Design choice:
the split lives in the reader at ingest** (`partition_shortwave`, §5.6) using the record's interval-mean
`cosz`, so the MEDS NetCDF can carry a single `SWdown` variable (smaller, source-faithful) and the one
authoritative partition lives in Fortran; the formatter *may optionally* pre-split and write four streams
(the reader then passes them through, `sw_partition=passthrough`). Either way `met_record_t` ends up with
the four streams (§3.2). Standard constants: PAR ≈ 0.45–0.50 of total SW *energy*; `1/4.6 ≈ 0.2174` W per
µmol PAR (the 4.6 µmol J⁻¹ PAR conversion).

**CO₂ — ERA5-Land has no CO₂**; `co2 = fcfg%co2_const` (default **420 µmol/mol**, a present-day value),
the single free-atmosphere authority that fans out to `aero_env_t%co2_atm`, `column_forcing_t%co2_atm`, and
the CAS-CO₂ reference (the brief flags this value is duplicated across ~5 types — one config authority
writes it once). A transient/observed CO₂ stream (Mauna Loa / CAMS) is a P2 add-on (a CO₂-only variable on
the same file, non-cycling).

**Wind measurement height (§10).** ERA5-Land wind is diagnostic at **10 m**; a temperate forest at Ithaca is
~20–30 m tall, so 10 m can sit *below* canopy top. Like every reanalysis-forced DGVM (GSWP3/CRUNCEP/WFDE5),
P0 applies the near-surface reanalysis values at the model reference height and accepts the height
mismatch; `wind_meas_height` records the source height for a future log-profile adjustment (§10).

### 5.3 Wind — energy-conserving interpolation

Wind (and `ustar` if ever ingested) interpolate as **squared** quantities then square-rooted (ED2's
`vels` accumulates `u²`), so the kinetic-energy content is conserved across record boundaries:
`wind(now) = sqrt(w_prev·wind_prev² + w_next·wind_next²)`, floored at a `u_min` (numerical stability
of the M–O similarity in aerodynamics).

### 5.4 Precipitation — step-constant, phase-split at ingest

`Rainf` is held **step-constant** across the interval (never linearly interpolated — interpolation
smears intense events into physically wrong drizzle and breaks infiltration/runoff). Phase split
(`precip_phase`, ED2 Jin-1999 bands on `Tair`) partitions total precip into `rainf`/`snowf`; unlike
tropical BCI, **Ithaca has real sub-freezing precip**, so `precip_phase` fires on the actual ERA5-Land
file (tested, §9.8). Precip enters the fast loop as `column_forcing_t%precip` (ground-reaching,
post-interception) and `intercept_canopy_layer`'s `rain_above`.

### 5.5 Missing data — **MEDS never gap-fills; it errors (revision comment 2)**

**MEDS does no gap-filling, ever.** There is no `gap_policy` selector, no persistence hold, no
climatology fill. If any *required* forcing field at a needed time step is **missing or NaN**, the reader
**halts** (`error stop` / non-zero status) naming the field and timestamp. Gap-filling — if a dataset ever
needs it — is entirely **the user's responsibility, upstream and outside MEDS** (use a gap-filled product,
or fill it in your own preprocessing before the file reaches the MEDS formatter/reader). ERA5-Land is
gap-free over land, so this is a validation guard, not a routine path; making it a hard error keeps a
silently-degraded run from ever masquerading as a clean one.

Scope of "missing": the check is on genuine *data gaps* (a `_FillValue`/NaN in a variable that should have
a value). It does **not** touch the **necessary derivations** — the total→4-stream SW partition (§5.6) and
specific-humidity-from-dewpoint (§5.2) — which compute variables the source doesn't store directly from
variables it does; those are deterministic conversions, not gap-fills. The de-accumulation boundary value
(the single leading non-01Z sample, §7.3) is likewise not a data gap — the formatter drops it, so it is
never written to the MEDS file.

**The formatter does not gap-fill either** (`prep_era5land_forcing.py`, §7): it errors on an unexpected
NaN in a source variable rather than inventing a value, so the MEDS forcing NetCDF it writes is either
complete or the pipeline stops — consistent with the reader's contract.

### 5.6 SW partition seam (`partition_shortwave`) — **P0-required** for ERA5-Land (total SW)

A pluggable dispatch (ED2 `imetrad` analogue) over total SWdown + `psurf` + `cosz` → 4 streams:
`SWPART_PASSTHROUGH` (use the file's own split — for a *pre-split* file, e.g. a future 4-stream tower or
`prep_era5land_forcing.py --presplit`), `SWPART_WEISS_NORMAN` (WN85, band-specific PAR/NIR diffuse
fractions, needs pressure), `SWPART_SIB` (Sellers-86), `SWPART_CLEARIDX` (Boland/Tsubo/Erbs clearness
index). Ported as `pure` kernels from `radiate_utils.f90`. **Because ERA5-Land ships only total `ssrd`,
the P0 default is `SWPART_CLEARIDX`** (the simplest defensible total→direct/diffuse split, then a fixed
PAR fraction ≈ 0.45–0.50 of energy for the PAR/NIR split); `SWPART_WEISS_NORMAN` is the P1 upgrade
(band-specific diffuse fractions). `SWPART_PASSTHROUGH` is used only when the file already carries the
four streams (`sw_input_kind = "fourstream"`).

### 5.7 LWdown source / synthesis (`synthesize_lwdown`, `lwdown_source`)

When the file lacks LWdown or it is NaN (`LW_SYNTHESIZE`), reconstruct from `Tair`+`Qair` with a
clear-sky emissivity (Brutsaert 1975 / Idso-Jackson 1969) plus an optional cloud correction from the
SW clearness index: `LWdown = ε_clear·σ·Tair⁴·(1 + a·(1−kt))`. **ERA5-Land provides `strd`** (downward
longwave), so P0 uses `LW_FILE`; `LW_SYNTHESIZE` is the fallback for a future source without longwave.

---

## 6. Wiring — how forcing enters the loops

### 6.1 The injection points (exact lines)

Two private routines in `meds_fast_loop` translate the constant `fast_context_t` into the kernel
inputs, and are the **exact seam** a real forcing source replaces:

- **`build_forcing(forc, coh, ctx, sum_lai)`** → `column_forcing_t`: today reads `ctx%air_temp`,
  `ctx%shv_atm`, `ctx%co2_atm`, `ctx%rad_sw_ground`, `ctx%precip`, and LAI-splits `ctx%rad_sw_top`
  across cohorts. It already computes `enthalpy_atm = cas_enthalpy_of_temp(ctx%air_temp, ctx%shv_atm)`.
- **`fill_aenv(aenv, bio, ctx)`** → `aero_env_t`: reads `ctx%u_ref/zref/press/rho_air/air_temp/
  shv_atm/co2_atm` and the live CAS/ground state.

### 6.2 Minimal-diff wiring: refresh `fast_context_t`'s met fields per sub-step

The cleanest change that preserves `build_forcing`/`fill_aenv` verbatim: **make `fast_context_t`'s met
scalars time-varying**, refreshed from `met_instant` at the top of *each fast sub-step*. Because the
diurnal cycle lives across the `n_fast_per_slow` sub-steps (§1.1), the refresh must be **inside** the
`isub` loop of `run_fast_biophysics`, and the loop must know the sub-step wall-clock time.

Concretely, `run_fast_biophysics` gains the current date and the driver buffer, and the sub-step loop
rebuilds forcing each sweep:

```fortran
subroutine run_fast_biophysics(site, ctx, cfg, drv, now, worst_energy, worst_water, n_budget_fail)
   type(fast_context_t), intent(inout) :: ctx        ! met scalars refreshed each substep
   type(met_driver_t),   intent(inout) :: drv        ! reader buffer (mutable)
   type(meds_time_t),    intent(in)    :: now        ! start-of-slow-step date
   ...
   do isub = 1_ik, cfg%n_fast_per_slow
      t_sec  = seconds_into_day(now) + real(isub-1, wp)*cfg%dt_fast + 0.5_wp*cfg%dt_fast   ! substep midpoint
      t_sub  = time_advance_seconds(now, (isub-1)*nint(cfg%dt_fast))
      call met_advance(drv, t_sub)
      met = met_instant(drv, t_sub, t_sec, cfg%forcing%latitude_deg)     ! §4.4
      call apply_met_to_ctx(ctx, met)                                    ! write met scalars into ctx
      call build_forcing(forc, coh, ctx, sum_lai)                        ! UNCHANGED body
      call fill_aenv(aenv, bio, ctx)                                     ! UNCHANGED body
      call column_fast_step(cfg%dt_fast, cfg, ctx%ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh=gpp_coh)
      ... gpp_accum integration (unchanged) ...
   end do
```

`apply_met_to_ctx(ctx, met)` sets `ctx%air_temp = met%tair_k`, `ctx%shv_atm = met%qair`,
`ctx%press = met%psurf_pa`, `ctx%rho_air = met%rho_air`, `ctx%co2_atm = met%co2`,
`ctx%u_ref = met%wind`, `ctx%precip = met%rainf`, and `ctx%rad_sw_top = met%swdown()`. Four points the
earlier draft glossed:

- **Ground SW.** Today's `rad_sw_ground = 60` (`meds_fast_loop.f90:42`) is an *independent* ad-hoc
  constant, not a canopy-top derivative. To keep CONST a faithful reproduction, `apply_met_to_ctx`
  sets `ctx%rad_sw_ground = fcfg%rad_sw_ground_const` under the **CONST** backend (default 60), and
  `= met%swdown()·exp(−k·LAI)` (Beer's law) under CSV/NetCDF. At the P1 RT join (§6.3) this is
  replaced by the two-stream `abs_sw_ground`.
- **`theta_atm` is fed raw `T`, and the type label was misleading.** `fill_aenv` sets
  `aero_env_t%theta_atm = ctx%air_temp` (`meds_fast_loop.f90:220`) with **no Poisson
  `(p0/p)^(Rd/cp)` conversion**, so the "potential temp at zref" label in the §1.1 table is aspirational.
  MEDS uses the near-surface approximation **`θ ≈ T`** at P0 (the reference height is ~40 m, the factor
  is ≈1); the table label is corrected to "air temp at zref (θ≈T near surface)". Adding the true Poisson
  factor is a one-line P1 refinement in `fill_aenv`, not a forcing-reader concern.
- **`zref`.** `ctx%zref` defaults to **30 m** (`meds_fast_loop.f90:35`), but the configured BCI
  `[site].reference_height` is **40 m**. `meds_main` threads `cfg%site%reference_height → ctx%zref`
  when it builds the static context (§6.5) so the aerodynamics see the true measurement height; this is
  a build-time wiring line, not part of the per-substep refresh.
- **`lwdown` is carried but DEAD at P0.** `met%lwdown` (and `met%par_beam/par_diffuse/nir_*`,
  `met%cosz`) are handed to the RT seam **only once §6.3 lands**. At P0, `build_forcing` sets only
  `abs_sw`/`abs_sw_ground`; `column_forcing_t%abs_lw` is left unset until RT arrives, so **`lwdown`
  cannot affect any P0 result** and P0 tests must not assert on it.

Because the four SW streams sum to 400 and CONST holds them temporally flat (no cosz disaggregation —
a constant climate is diurnally flat by construction), the **canopy-top SW reproduces today's
`rad_sw_top = 400` exactly**, and air T / humidity / CO₂ / wind / pressure / ρ_air match the current
constants. This makes CONST a **near-no-op refactor** — but **not** the "bit-for-bit" claim the earlier
draft made: `cosz`/`rho_air` are recomputed each substep, and ground SW now flows through
`rad_sw_ground_const`. Test 4 (§9) is re-specified accordingly to compare against a documented CONST
baseline (canopy-top SW exact; budgets within round-off), **not** a byte-identical replay.

> **Design note — cleaner alternative (P1):** thread a `met_forcing_t` argument straight into
> `build_forcing`/`fill_aenv` (and later `column_forcing_t`), retiring the met scalars on
> `fast_context_t` entirely so `fast_context_t` carries only the *static* `column_config_t` +
> initial soil state. This removes the `apply_met_to_ctx` shim and the co2/temp duplication. P0 keeps
> the shim for minimal churn; P1 does the type-level cleanup once the RT-driven per-cohort SW lands.

### 6.3 Per-cohort shortwave — RT replaces the LAI-share split

Today `build_forcing` splits `rad_sw_top` across cohorts by **LAI share** (a placeholder). The real
per-cohort absorbed SW/PAR comes from **canopy RT** (`meds_canopy_radiation`), whose `rad_forcing_t`
consumes exactly `incid_beam(band)`, `incid_diff(band)`, `cosz` — which `met_forcing_t` now supplies
(`par_beam/par_diffuse/nir_beam/nir_diffuse` → the two bands, `cosz` from `solar_cosz`, `lwdown` →
`incid_diff(RAD_LW)`). The wiring: `run_fast_biophysics` calls `canopy_radiation(opt, rad_forcing_from
(met), coh, …)` to fill `column_forcing_t%abs_sw(:)`/`abs_lw(:)` per cohort and `abs_sw_ground`/
`abs_lw_ground`, **replacing** the LAI-share split. This is the RT↔forcing join; it is a P1 item
(P0 keeps the LAI split driven by the *time-varying* `rad_sw_top`, which already yields a diurnal GPP
cycle).

### 6.4 Slow-loop consumption — daily aggregation

The slow loop (`meds_vegetation_dynamics`) consumes forcing two ways:

1. **Via the fast→slow GPP handoff (primary).** `carbon_growth` reads `a_carbon =
   cohort%gpp_accum(j)` when `fast_biophysics_on`. `run_fast_biophysics` already integrates
   sub-daily GPP onto `gpp_accum`; with real time-varying SW forcing this becomes a **physically
   meaningful daily GPP** (today's constant SW gives a flat, unphysical value). No new slow-loop
   plumbing — the met source improves the *content* of `gpp_accum`, not the seam.
2. **Directly, for slow-loop drivers that need daily aggregates.** `pheno_env_t%temp_day` (daily-mean
   air T), `%avail_water`, `%daylength`, `%doy` are daily/seasonal cues. A small
   **daily accumulator** on the driver (sum `met%tair_k`, count, and `∫SWdown` over the day's fast
   sub-steps; `daylength` from `cosz > cosz_min` count; `doy = day_of_year(now)`) produces a
   `met_daily_t` the slow loop reads at day roll-over. P0 defers phenology forcing (the phenology
   module is a directional signal, already testable standalone); the accumulator is a P1 add when
   phenology is wired into the stepper.

### 6.5 Closing the production wiring gap (`meds_main`)

`meds_main` must (a) build a `met_driver_t` from `[forcing]`, (b) build a `fast_context_t` (static
`column_config_t` + initial soil from `[fast]`/`[site]`), **threading `cfg%site%reference_height` into
`ctx%zref`** so the aerodynamics use the configured 40 m rather than the 30 m field default
(`meds_fast_loop.f90:35`), (c) call `init_fast_reservoirs(site, ctx)` once, and (d) pass `fast_ctx`
**and** the driver/date into `advance_one_step`. Since
`advance_one_step(site, cfg, is_new_month, is_new_year, fast_ctx)` already has the optional `fast_ctx`
hook (and gates on `cfg%fast_biophysics_on .and. present(fast_ctx)`), it grows two optional args
`met_drv`/`now` threaded to `run_fast_biophysics`:

```fortran
! meds_main, replacing the line-124 call:
if (cfg%fast_biophysics_on) then
   call advance_one_step(site, cfg, is_new_month, is_new_year, fast_ctx=ctx, met_drv=drv, now=prev)
else
   call advance_one_step(site, cfg, is_new_month, is_new_year)
end if
```

**Latent-coupling hazard the wiring must avoid (brief concern):** if `fast_biophysics_on = true` but
`fast_ctx` is absent (today's `meds_main`), `carbon_growth` reads a `gpp_accum` that
`run_fast_biophysics` never filled/reset — stale/zero GPP. The fix is exactly this always-pass-when-on
wiring; `validate_config` additionally asserts `forcing_on` when `fast_biophysics_on` and
`growth_source == GS_CARBON` (no silent stub-GPP carbon run).

### 6.6 The `[forcing]` and `[site]` config blocks

New TOML sections, loaded by `meds_config_io` with the existing `req_*` helpers (`req_l`, `req_dur`
for durations, `req_r`, `req_i`, and a new `req_str` for paths / a mapped `req_enum`); every key
presence-mapped and hard-error on absence (the house rule — no defaults in the loader):

```toml
[site]
latitude         = 42.44       # [deg +N]  Ithaca NY reference cell
longitude        = -76.50      # [deg +E]  drives the solar-longitude/EoT correction (§5.1)
utc_offset       = 0.0         # [h]       ERA5-Land is UTC -> file clock IS UTC (offset 0)
elevation        = 320.0       # [m]       Ithaca ~320 m; reserved for lapse (P2)
reference_height = 40.0        # [m]       forcing reference height; MUST exceed every PFT hgt_max (asserted)
wind_meas_height = 10.0        # [m]       ERA5-Land wind is 10 m (below tall canopy — §5.2/§10)

[forcing]
forcing_on       = true        # master gate; independent of [fast].fast_biophysics_on
format           = "netcdf"    # "netcdf" (MEDS multi-grid forcing file) | "const" (reference climate)
path             = "data/forcing/ithaca_era5land_forcing.nc"
grid_index       = 1           # which (time,grid) location this polygon reads (§4, §7); 1 = single site
timestep         = "3600s"     # native forcing interval (ERA5-Land hourly)
avg_convention   = "end"       # flux vars are means over the hour ENDING at the stamp; states instantaneous
apply_solar_lon  = true        # UTC file -> longitude gives local solar time (MANDATORY for a UTC source, §5.1)
sw_partition     = "clearidx"  # ERA5-Land ships TOTAL SW -> partition required. "clearidx"(P0) | "weiss_norman"(P1) | "passthrough"(pre-split file)
lwdown_source    = "file"      # "file" (ERA5-Land strd) | "synthesize" (Brutsaert, for sources lacking LW)
# (no gap_policy key: MEDS never gap-fills — a missing required value hard-errors, §5.5)
co2_const        = 420.0       # [umol/mol] free-atm CO2 (ERA5-Land has no CO2 -> single config authority)
rad_sw_ground    = 60.0        # [W/m2] CONST-backend ground SW (reproduces today's constant)
start_clamp      = "error"     # model start < first record: "error" | "hold" (§4.2)
recycle          = true        # cycle the record when the run outruns the file (spin-up)
```

`req_dur("forcing.timestep")` reuses the duration parser; `format`/`avg_convention`/`sw_partition`/
`lwdown_source`/`start_clamp` map strings → the `MET_*`/`METAVG_*`/`SWPART_*`/`LW_*`/`CLAMP_*` codes via
small mappers (the `req_scheme` pattern). `validate_config` asserts:
`reference_height > max(hgt_max)` (ED2's abort rule), `latitude ∈ [−90,90]`, `dt_fast ≤ timestep`,
**`mod(timestep, dt_fast) == 0`** (the forcing interval must be an integer multiple of `dt_fast` for
clean sub-step windowing — 3600/900 = 4 for ERA5-Land), `1 ≤ grid_index ≤ ngrid` (against the file's
`grid` dimension, §4.1), and — for the NetCDF backend — that `start_time`/`end_time` lie **within the
forcing record grid** (else the `start_clamp`/`recycle` paths of §4.2 govern, and the assertion documents
which).

**UTC vs local-standard-time policy (now the P0 default, §5.1).** ERA5-Land timestamps are **UTC**
(`utc_offset = 0`), so `apply_solar_lon = true` is **mandatory** at P0 — `solar_cosz` needs **local
apparent solar** seconds, and the single conversion seam supplies them:
`t_sec_solar = sec_of_day(now) − utc_offset·3600 + (longitude·240 + eot(doy))`, gated by
`apply_solar_lon`. For a UTC source the `utc_offset` term is 0 and the longitude term carries the whole
local-time shift; a future **local-clock** tower file sets `utc_offset ≠ 0` and the same one-line seam
handles it — one code path, both conventions.

---

### 6.7 Forcing provenance — echo the assembled record to the diagnostic stream

Debugging the diurnal reconstruction (the cosz disaggregation, the avg-convention phase, the solar-time
offset) is nearly impossible from GPP alone. So the assembled instantaneous `met_forcing_t` is
**echoed to the diagnostic netCDF** (`meds_io`'s `-D-output.nc` stream) at the diagnostic cadence: the
13 scalar fields (`tair_k`, `qair`, `psurf_pa`, `rainf`, `wind`, `lwdown`, the four SW streams, `co2`,
`cosz`, `rho_air`) as site-level timeseries variables, plus `swdown()`/`rshort_diffuse()` diagnostics.
This is a small, additive `io_write_snapshot` extension (site-scalar variables, no new dimension) and
is the single most useful artefact for validating that the reconstructed SW tracks the measured `Rs`
and peaks at the right local time. It is **write-only provenance** — nothing consumes it — so it is
safe to add at P0 and costs one row of scalars per output record.

## 7. ERA5-Land prep — two scripts + the MEDS multi-grid forcing NetCDF (comments 2 & 3)

Two standalone scripts in `scripts/` (dependency-light: `cdsapi`, `xarray`/`netCDF4`, `numpy`), split by
concern so the slow network download is separate from the fast, re-runnable formatting:

1. **`scripts/download_era5land.py`** — pulls raw ERA5-Land hourly NetCDF from the CDS for a requested
   lat/lon (or bounding box) and date range.
2. **`scripts/prep_era5land_forcing.py`** — converts the raw ERA5-Land NetCDF into the **MEDS multi-grid
   forcing NetCDF** the Fortran reader consumes (de-accumulate fluxes, unit-convert, humidity from
   dewpoint, wind magnitude, optional SW pre-split).

### 7.1 The MEDS forcing NetCDF format (multi-grid, canonical — comment 3)

The single file format the reader reads (§4). **Dimensions:** `time` (UNLIMITED) × `grid` (number of
locations — `1` for a single site). **Coordinates:** `time(time)` (`units="seconds since <base>"`,
`calendar="proleptic_gregorian"`), `latitude(grid)` (`degrees_north`), `longitude(grid)` (`degrees_east`),
optional `elevation(grid)` (`m`). **Forcing variables, all shaped `(time, grid)`**, ALMA-named:

| Variable | Units | `cell_methods` | Notes |
|---|---|---|---|
| `Tair` | K | `time: point` | instantaneous state |
| `Qair` | kg/kg | `time: point` | from dewpoint (§7.3) |
| `PSurf` | Pa | `time: point` | |
| `Wind` | m/s | `time: point` | √(u10²+v10²) |
| `Rainf` | kg/m²/s | `time: mean` | de-accumulated hourly-mean rate |
| `LWdown` | W/m² | `time: mean` | de-accumulated hourly-mean flux |
| `SWdown` | W/m² | `time: mean` | **total** (default); reader partitions (§5.6). *Or* the four `SWdown_{par,nir}_{beam,diffuse}` if pre-split |
| `CO2air` | µmol/mol | `time: point` | optional; absent ⇒ reader uses `co2_const` |

**Global attributes:** `Conventions = "MEDS-forcing-1.0"` (ALMA-compatible names), `title`,
`source` (e.g. `"ERA5-Land hourly (reanalysis-era5-land)"`), `timestep_seconds` (3600),
`avg_convention = "end"` (flux vars are means over the hour ENDING at the stamp; state vars are point
values at the stamp), `sw_input_kind = "total" | "fourstream"` (tells the reader whether to partition),
`time_zone = "UTC"`, and provenance (download date range, cell indices, CO₂ constant + note, script
version). A single-site file is just `grid = 1`; a multi-cell file lists each location in `latitude(grid)`
/ `longitude(grid)` and the reader selects `grid_index` (§4).

> **Layout choice — unstructured `(time, grid)`, not `(time, lat, lon)`.** MEDS polygons are an arbitrary
> *list* of locations ("a location sharing one meteorological forcing"), so an unstructured `grid`
> dimension with `lat(grid)`/`lon(grid)` maps 1:1 to the polygon list, wastes no cells on a regular box,
> and makes nearest-match trivial (identity when the file is built from the polygon locations). A regular
> `(time, y, x)` reanalysis tile is representable by flattening `y×x → grid`. Single site = degenerate
> `grid = 1`.

### 7.2 `download_era5land.py` — CDS API (ERA5-Land hourly) — verified 2026-07-08

Uses `cdsapi >= 0.7.7` against the **current CDS** (post-2024 migration): dataset id
**`reanalysis-era5-land`** (hourly, 0.1°, UTC); `~/.cdsapirc` has exactly two lines
`url: https://cds.climate.copernicus.eu/api` (**no `/v2`**) and `key: <Personal-Access-Token>` (single
token, **no UID**); request keys `data_format:"netcdf"`, `download_format:"unarchived"`. Requests the eight
variables MEDS needs (`2m_temperature`, `2m_dewpoint_temperature`, `surface_pressure`,
`total_precipitation`, `10m_u_component_of_wind`, `10m_v_component_of_wind`,
`surface_solar_radiation_downwards`, `surface_thermal_radiation_downwards`) over a small `area`
**`[North, West, South, East]`** box around the requested lat/lon. Ithaca NY defaults: `lat=42.44,
lon=-76.50` ⇒ `area=[42.55, -76.60, 42.35, -76.40]`. **Two verified gotchas the script handles:**
- **ZIP fallback.** Since the 2024-11-26 CDS converter update, a NetCDF request that **mixes instantaneous
  (t2m/d2m/sp/u10/v10) and accumulated (tp/ssrd/strd) fields** is returned as a **`.zip`** (one `.nc` per
  GRIB `stepType`) *even with* `download_format:"unarchived"`. The script detects the zip magic bytes,
  extracts, and merges the member `.nc` files into one dataset (or, with `--split-requests`, issues two
  separate single-stepType requests). Do not assume a single `.nc` comes back.
- **De-accum needs a trailing day.** To recover the last day's `[23Z,00Z]` hour, the script **pads the
  request by one extra day** past the requested end (§7.3, §C of the fact sheet).

The downloaded NetCDF's variable names are `t2m`, `d2m`, `sp`, `u10`, `v10`, `tp`, `ssrd`, `strd`, on
`(valid_time|time, latitude, longitude)`; the formatter reads those names (not the GRIB shortnames).

### 7.3 `prep_era5land_forcing.py` — raw ERA5-Land → MEDS forcing NetCDF

The load-bearing conversions (§5.2). All operate per selected `grid` cell, then stack into `(time, grid)`:

```python
# --- humidity from dewpoint (must match meds_thermo saturation form, §4.3) ---
def dewpoint_to_q(Td_K, P_Pa):
    Tdc  = Td_K - 273.15
    esat = 611.2*np.exp(17.67*Tdc/(Tdc + 243.5))     # e(Td) = saturation vapor pressure AT the dewpoint
    return 0.622*esat/(P_Pa - 0.378*esat)            # specific humidity [kg/kg]

# --- de-accumulate an ERA5-Land accumulated field (accum since 00 UTC, hourly) to a per-hour amount ---
# ERA5-Land accumulations run from 00 UTC and reset daily (steps 1..24). THE 00:00 UTC STAMP IS THE
# WHOLE PREVIOUS DAY'S TOTAL (step 24), NOT zero and NOT a small hourly value — the #1 off-by-one trap.
# Correct recipe (verified CDS fact sheet, 2026-07-08): ordered by valid time, per-hour(H) = raw(H) -
# raw(H-1) EVERYWHERE (incl. the 00:00 stamp, which then yields the prior day's last hour) EXCEPT at
# 01:00 UTC, where per-hour = raw(01:00) AS-IS (the first step of a period, implicit 0 at 00:00).
def deaccumulate_hourly(accum, valid_time_utc):      # accum (time,), valid_time_utc sorted ascending
    per_hour = np.empty_like(accum)
    for i in range(len(accum)):
        if valid_time_utc[i].hour == 1:              # 01:00 UTC: first hour of the day, take as-is
            per_hour[i] = accum[i]
        elif i == 0:
            per_hour[i] = np.nan                      # cannot difference the very first non-01Z sample
        else:
            per_hour[i] = accum[i] - accum[i-1]
    return np.where(per_hour > NEG_EPS, per_hour, 0.0)  # clip GRIB-packing negatives (fact sheet)
# NOTE: to fully de-accumulate day D you need the 00:00 UTC stamp of day D+1 (it carries the [23Z,00Z]
# hour) — so download.py fetches ONE EXTRA trailing day (§7.2).

# --- per cell ---
Tair   = t2m                                         # K
PSurf  = sp                                          # Pa
Qair   = dewpoint_to_q(d2m, sp)                      # kg/kg
Wind   = np.hypot(u10, v10).clip(min=U_MIN)          # m/s  (10 m; height noted in attrs, §10)
Rainf  = deaccumulate_hourly(tp,   time_utc) * RHO_W / 3600.0   # m/hr -> kg/m2/s  (RHO_W = 1000)
SWdown = deaccumulate_hourly(ssrd, time_utc) / 3600.0          # J/m2/hr -> W/m2  (total; reader splits)
LWdown = deaccumulate_hourly(strd, time_utc) / 3600.0          # J/m2/hr -> W/m2
CO2air = np.full_like(Tair, CO2_CONST)               # ERA5-Land has none (attribute-flagged synthetic)

write_meds_forcing_nc(out_path,
    time_seconds, lat_grid, lon_grid, elev_grid,     # coords: time(time), lat/lon/elev(grid)
    Tair, Qair, PSurf, Wind, Rainf, LWdown, SWdown, CO2air,   # each (time, grid)
    attrs=dict(Conventions="MEDS-forcing-1.0", source="ERA5-Land hourly (reanalysis-era5-land)",
               timestep_seconds=3600, avg_convention="end", sw_input_kind="total",
               time_zone="UTC", co2_const=CO2_CONST))
```

The de-accumulation is the one step most easily gotten wrong (an off-by-one zeroes or doubles the first
hour of every UTC day) — §9 asserts a known daily SW/precip integral against the raw accumulated totals.
The **SW split** is deferred to the Fortran reader by default (`sw_input_kind="total"`, so the one
authoritative partition lives in `partition_shortwave`, §5.6); a `--presplit` flag can instead write the
four streams (`sw_input_kind="fourstream"`) for models that want them in-file. The formatter can stack
**multiple cells** into the `grid` dimension in one call (one MEDS file covering N locations), which is how
a future multi-polygon run gets its forcing.

---

## 8. Phasing

**P0 — single-site NetCDF reader wired to the fast loop, ERA5-Land at Ithaca (MVP).**
- New **`src/forcing/` library** (§2): `met_forcing_t`, `met_record_t`, `met_driver_t` (with
  `grid_index`); `forcing_config_t` in `shared`; `[site]` + `[forcing]` config; `req_*` loaders. Promote
  `meds_netcdf_c` to its own always-built target; **delete `MEDS_ENABLE_IO` + `meds_io_stub`** (netCDF is
  mandatory now, comment 1); `meds_aux` links `meds_forcing`.
- **`meds_time` second-level helpers** (`seconds_into_day`, `time_advance_seconds`, `seconds_between`)
  — a hard prerequisite; the interpolation and window slide do not compile without them (§2.1).
- `meds_met_driver`: **NetCDF backend** (the MEDS multi-grid `(time,grid)` file, reading the fixed
  `grid_index` slab) + **CONST** backend (CSV retired); `met_open`/`advance`/`instant`/`close`;
  **linear** interpolation for state vars, **step-constant** precip, **reciprocal-mean-cosz** SW
  reconstruction (`solar_cosz` + the new `cosz_reconstruct_factor` — NOT ED2's `⟨sec z⟩`),
  `meds_thermo%air_density` (reused), `dewpoint_to_specific_humidity`, and — **required now, not deferred**
  — the total→(direct/diffuse)×(PAR/NIR) **`partition_shortwave`** (§5.6), since ERA5-Land ships only total
  SW. UTC solar-time path (`apply_solar_longitude=.true.`, §5.1).
- Wiring: per-sub-step met refresh in `run_fast_biophysics` (§6.2 shim); `meds_main` builds the driver,
  seeds reservoirs, passes `fast_ctx` + date (closes the gap). Per-cohort SW stays the LAI-share split of
  the **time-varying** `rad_sw_top` (real RT deferred).
- **ERA5-Land prep scripts (§7):** `download_era5land.py` (CDS API) + `prep_era5land_forcing.py`
  (→ MEDS multi-grid NetCDF), with **Ithaca NY** as the reference cell. The **file format is multi-grid
  from P0** (`grid` dimension present); the reader reads `grid_index = 1`.
- Test: `test_met_driver` reproduces a diurnal GPP/temperature cycle offline from the Ithaca file (§9).

**P1 — physical fidelity.**
- **Canopy RT join** (§6.3): `rad_forcing_t` from `met_forcing_t`, per-cohort absorbed SW/PAR from
  `meds_canopy_radiation` replaces the LAI split; `cosz`/band SW/`lwdown` fully wired.
- **Weiss–Norman band-specific** SW partition (upgrade from the P0 clearness-index split),
  **LWdown synthesis** (Brutsaert, for sources lacking it — ERA5-Land has `strd`; this is variable
  *derivation*, not gap-filling, §5.5), **multi-year cycling** (`recycle`, Feb-29 reconciliation), the
  **`met_forcing_t`-argument cleanup** (retire the `apply_met_to_ctx` shim), the **daily accumulator**
  feeding phenology (`temp_day`, `daylength`, `doy`). (No gap-fill climatology — MEDS never gap-fills, §5.5.)

**P2 — multi-polygon runtime + gridded science.**
- **Multi-polygon runtime**: one `met_driver_t` per polygon, each binding to its `grid` slice by
  **nearest-location match** of the polygon lat/lon against `latitude(grid)`/`longitude(grid)` (ED2
  `match_poly_grid`) — the file format already supports this from P0 (§7), so this is reader/driver work
  only. Elevation **lapse-rate** correction (`lapse.f90` analogue) between grid cell and site, a
  **10 m→reference-height wind log-profile** adjustment (§5.2/§10), climate-change intercept/slope
  perturbations, other reanalysis products via the same format (CRUNCEP/GSWP3/WFDE5/full ERA5), and a
  **transient/observed CO₂** stream (a CO₂-only, non-cycling variable on the same file).

---

## 9. Test plan

Offline, driving the physics kernels with the **Ithaca NY ERA5-Land** file to reproduce a **diurnal
cycle** — the model-behaviour analogue of ED2's EDTS, and the first end-to-end use of the fast loop in
production. Reader-unit tests use small **synthetic NetCDF** fixtures (no CDS download in CI).

1. **Reader unit tests (`test_met_driver`).** (a) `met_open` on a synthetic **multi-grid** NetCDF
   (`grid=2`, 3 time records) loads records #1–2 **at the configured `grid_index`**, and reading
   `grid_index=2` returns that cell's (distinct) values — proving the `(time,grid)` slab read (§4.1);
   (b) `met_advance` slides the window correctly across a record boundary, **recycles** at EOF, and takes
   the **start-before-`base_time`** path (`start_clamp=error` stops; `=hold` clamps to record #1 with
   `w_next=0`); (c) `met_instant` at the exact record time returns the record value; at the midpoint
   returns the linear mean (state vars) / the previous value (precip, step-constant); (d) **conservation**
   — the reciprocal-mean-cosz SW reconstruction integrated over a window returns the interval mean to
   round-off: `(1/T_win) ∫ F_avg·cosz(t)/⟨cosz⟩_win dt = F_avg` (this closes **because** the factor is
   `1/⟨cosz⟩_win`, not `⟨sec z⟩`; a companion assert shows the `⟨sec z⟩` form does **not** close and is
   biased high on a sunrise window), and step-constant precip conserves accumulation.
2. **Solar-geometry sanity (UTC file).** `solar_cosz` with the **UTC + longitude** path
   (`apply_solar_longitude=.true.`) peaks near *local* noon at Ithaca (42.44 °N, −76.5 °E) — i.e. ~17:00
   UTC, **not** 12:00 UTC — proving the longitude term is applied; `cosz_reconstruct_factor` returns
   `1/⟨cosz⟩_win` over daytime sub-samples (night clamped to 0); night SW is exactly 0 (`cosz ≤ cosz_min`).
   Cross-check day-length against the astronomical value for the date/latitude, and **assert the modeled
   solar-noon time matches local apparent noon within one window** (bounding the EoT term, §5.1). A
   mid-latitude site also exercises the strong seasonal day-length swing BCI's tropics do not.
3. **Kernel offline drive (the headline).** Feed one day of the **Ithaca ERA5-Land** file (24 hourly
   records → 96 fast sub-steps) into `column_fast_step` for a single on-allometry cohort and assert a
   **physical diurnal cycle**: GPP tracks PAR (zero at night, peak near local noon), `leaf_temp` follows
   `Tair` + a radiative
   offset, transpiration tracks VPD, the CAS/soil budgets close each step (`budg%whole_*%worst` below
   tolerance). ERA5-Land has **no GPP** column, so GPP is checked for **physical plausibility** (sign,
   diurnal shape, midday magnitude in a reasonable µmol/m²/s range for a temperate broadleaf canopy) — not
   against an observed flux (that is a P1+ FLUXNET/AmeriFlux Ithaca-area comparison). Note **P0 leaf PAR is
   biased high**: `leaf_env%par = abs_sw/lai·par_per_w` converts the LAI-share of *total* SWdown, treating
   absorbed NIR as PAR-convertible, until the P1 RT join supplies true per-cohort absorbed PAR.
4. **De-accumulation correctness (formatter, the ERA5-Land-specific risk, §5.2/§7).** On a known 2-day
   synthetic ERA5-Land accumulation (`ssrd`/`tp`/`strd` ramping within each UTC day, resetting at 00 UTC),
   assert the formatter's de-accumulated **daily integral** equals the raw end-of-day accumulation
   (`Σ_hours SWdown·3600 == ssrd(23Z_end)` etc.), that the **first hour of each UTC day** is the accum
   value itself (not a cross-day difference), and that no hour is negative. This pins the off-by-one that
   would otherwise zero/double the first hour of every day.
5. **Multi-grid round-trip (comment 3).** The formatter writes a **2-cell** (`grid=2`) MEDS NetCDF from two
   distinct ERA5-Land cells; assert the file has a `grid` dimension of 2 with correct `latitude(grid)`/
   `longitude(grid)`, and that the Fortran reader at `grid_index=1` vs `2` returns each cell's series — a
   single-file, two-location forcing works end-to-end.
6. **CONST near-no-op refactor (re-specified — NOT bit-for-bit).** With `format = "const"`,
   `run_fast_biophysics` reproduces a **documented CONST baseline**: canopy-top SW is exactly 400 W/m²
   (the four reference streams sum to 400, held diurnally flat), and air T / humidity / CO₂ / wind /
   pressure match today's constants, so GPP and the CAS/soil budgets agree with the current
   `test_fast_loop` **within round-off**. It is **not** byte-identical (`cosz`/`rho_air` are recomputed
   each substep; ground SW flows through `rad_sw_ground_const=60`), and **`lwdown`/`abs_lw` is not
   asserted** (dead until the P1 RT join, §6.2).
7. **Python↔Fortran humidity agreement.** The formatter's `dewpoint_to_q` (§7.3) and the reader's
   `dewpoint_to_specific_humidity` must return the same `qair` for the same `(Td, P)` to round-off — the
   **single-saturation-formula** check (§4.3) guaranteeing the file's `Qair` and any reader-side humidity
   math agree.
8. **Precip-phase / frozen branch (Ithaca winter + synthetic).** Unlike tropical BCI, Ithaca has real
   sub-freezing precip, so `precip_phase` fires on the actual ERA5-Land file; a winter day (plus a
   synthetic `tair_k < t_3ple` record) asserts total precip splits conservatively into `rainf`+`snowf`
   (mass conserved, `snowf > 0` below the phase band).
9. **Fast→slow GPP handoff units.** With time-varying SW, assert the daily-integrated `gpp_accum` carries
   the correct units (`[kgC/plant/day]` at the seam `carbon_growth` reads) and is **reset once per slow
   day**, so a multi-day CONST run gives a constant daily GPP and a multi-day ERA5-Land run gives a
   day-varying one — pinning the reset/units the always-pass-`fast_ctx` wiring (§6.5) makes safe.
10. **Portability.** Build all new modules under **nvfortran multicore** (not just ifx): the
    `partition_shortwave` kernels return arrays — bind results to named arrays before any call (trap #7).
    Confirm the fast loop with real forcing runs on the multicore back end.

---

## 10. Open questions

1. **ERA5-Land de-accumulation edge (the sharpest risk — confirm against ECMWF docs, §5.2/§7).** The
   accumulation is since 00 UTC and resets daily; the exact treatment of the **00 UTC step** (is it 0, or
   does it carry the last hour of the previous day?) and the **01 UTC** first-difference determine whether
   the first hour of each UTC day is correct. Pinned by the verified CDS fact sheet (2026-07-08) and test 4;
   flagged here because a silent off-by-one biases the daily radiation/precip totals.
2. **10 m wind / 2 m T,q vs a tall canopy.** ERA5-Land near-surface diagnostics are at 10 m (wind) and 2 m
   (T,q), which can be **below** a 20–30 m forest canopy. P0 applies them at the model reference height as
   every reanalysis-forced DGVM does (GSWP3/CRUNCEP/WFDE5), recording `wind_meas_height`; a log-profile /
   Monin–Obukhov height adjustment from the measurement height to `reference_height` is a P2 refinement.
   Confirm `reference_height` still clears every PFT `hgt_max` (ED2 aborts if `zref ≤ hgt_max`).
3. **SW partition method at P0.** ERA5-Land gives only total SW, so a diffuse-fraction model is required
   from P0 (§5.6). The clearness-index (Erbs/Reindl/Spitters) split is the simplest defensible P0 choice;
   Weiss–Norman band-specific PAR/NIR diffuse fractions are the P1 upgrade. Both need `cosz` + TOA
   irradiance (have both); the choice affects the beam/diffuse ratio that canopy RT will consume at P1.
4. **CO₂ value / source.** ERA5-Land has no CO₂; P0 uses a single `co2_const` (420 µmol/mol). A
   transient observed CO₂ (Mauna Loa annual, or CAMS) is a P2 non-cycling stream on the same file.
5. **Where the daily phenology aggregation lives.** The `met_daily_t` accumulator (§6.4) could live on the
   driver (simplest) or on the site (visible to a future MPI/polygon layer). Deferred to P1 with phenology
   wiring; flag for the master-loop design.
6. **`met_forcing_t` on `fast_context_t` vs a first-class argument.** P0 uses the `apply_met_to_ctx` shim
   (minimal diff); P1 retires it. Decide at P1 whether `column_forcing_t`/`aero_env_t` should take
   `met_forcing_t` directly (cleaner, removes the co2/temp duplication the brief flags across ~5 types)
   versus keeping the `fast_context_t` carrier for backward compatibility with `test_fast_loop`.
7. **Multi-polygon grid binding.** The file is multi-grid from P0, but the P0 *reader* uses an explicit
   `grid_index`. When multi-polygon runtime lands (P2), the nearest-location match (ED2 `match_poly_grid`)
   and how a polygon with no covering land cell is handled (ERA5-Land is land-masked — coastal/urban cells
   can be missing) need specifying.
8. **Folder name (`forcing/` vs `input/`) — your call (§0).** Recommended `forcing/`; a pure rename if you
   prefer `input/`.
