# MEDS Diagnostic Aggregation & Output — Design

> **Revision (2026-07-09), against the current codebase.** Updated after the P3 fast-loop merge over two
> review passes. **Naming (ED2-faithful):** the per-step temporal reducers are the **`integrate_ / normalize_
> / reset`** family (`integrate` matches ED2's `integrate_ed_*mean_vars`; `integ_buffer_t`,
> `meds_output_integrate`, `output_integrate`) — was `accumulate/accum`; the per-step state reader is
> **`extract_variable`** / `extract_site` (pull the value *out of* the SoA and copy — was `harvest_var`);
> the registry helper is **`add_variable()`** (was `add`). **Every raw variable carries a scale suffix
> `_cohort` / `_patch` / `_site`** (`agb_cohort` vs `agb_site`; `_site` replaces the old `total_` prefix;
> fluxes drop `_rate` → `gpp_site`, `nee_site`, `transp_site`; identity/CSR/count fields exempt).
> **Other points:** (a) netCDF/IO is a hard dependency — the `MEDS_ENABLE_IO=OFF` **stub build was
> removed** (CLAUDE.md), so the `*_stub` twin is gone and the netCDF-free/serializer split is now DAG
> hygiene alone; (b) `mean_dbh`/`gpp_total`/`max_dbh` dropped; (c) the finest tier is the **fast** tier
> (`fmean`), a **multiple of `dt_fast`** like ED2, not "hourly"; (d) **cohort/patch structure changes only
> at the monthly and annual cadence**, so every cohort/patch output window ≤ 1 month sees a **fixed** slot
> set — the ragged-`global_id`-across-a-window machinery is deleted and **annual output is site-level
> only**; (e) the annual file drops the decade-bucket `0` (one whole-run file by default); (f) per-variable
> options move to an optional **`meds_io_config.toml`** named from the main config (mirroring the PFT
> config) that **overwrites** the main config's **high-level flux-group toggles** (`structure` /
> `carbon_fluxes` / `water_fluxes` / `energy_fluxes`); (g) the flux/state P1 additions report the
> **monthly mean only** — variance
> (`AGG_MEANSQ`) stays defined but deferred (no variable uses it yet).

A **stateful, streamed, multi-frequency** diagnostic-output subsystem for MEDS (`src/io/`) that turns
the instantaneous once-a-year snapshot the model writes today into a general **aggregation-and-output
engine**: every diagnostic variable declares itself once in a **registry** (name, long name, units,
dimensionality, aggregation operator, a variable **group**, the output streams it belongs to, and an
on/off flag), a set of per-frequency **integrators** fold it into a running dt-weighted reduction in
memory as the stepper ticks (fast → daily → monthly → annual, an ED2 `integrate_/normalize_/zero_`-style
`fmean → dmean → mmean` chain plus a yearly tier), and a
**flush** writes each closed period to a **time-chunked, per-stream netCDF file** through the existing
`meds_netcdf_c` C-binding layer. All of it is **user-mutable from TOML** — high-level flux-group toggles
and per-frequency switches in the main config, optional per-variable overrides in a separate
`meds_io_config.toml`, all with no source edit.

The reference architecture is **ED2's `var_table` + `average_utils` + `h5_output`** (registration ⊥
temporal averaging ⊥ serialization; a variable opts into *N* streams via *N* keyword flags, one
`integrate/normalize/zero` family per temporal tier). MEDS keeps that **decoupling** wholesale but
**rejects the mechanism**: ED2's `vt_info(maxvars,ngrids)` is a *global mutable registry of raw
pointers into the reallocatable SoA*, requiring a fragile `filltab_alltypes` re-hash after every
birth/death/fusion — exactly the hidden-global-mutable-state pattern MEDS/CLAUDE.md forbids. MEDS
instead carries an **explicit, passed** registry (`output_registry_t`), and the integrators **copy the
reduced value out of the live SoA each step** rather than caching pointers into it — so cohort/patch
count changes are a non-event (no re-hash, ever).

The design honours the three user goals directly: **PERFORMANCE** — in-memory integration + infrequent
buffered flush keeps I/O off the daily critical path (the load-bearing fact: daily/monthly/annual
flushes are rare and tiny); a P2 background-thread writer is an **optional, build-gated** extra for the
one dense case (the fast stream), *not* a load-bearing part of the perf story (§5.5); **FILE SIZE** —
one file **per stream per time-chunk** (a sim-year or sim-month), never one lump, all chunked +
deflated; **VERSATILITY** — user-mutable fast/daily/monthly/annual streams with high-level flux-group
toggles and an optional per-variable on/off override, from TOML.

Adding a variable is **two coupled edits, not one** — an honest accounting of the price of MEDS's
no-pointer registry: one `add_variable()` line in the registration list *and* one `case` in the
`extract_variable` switchboard that copies the value out of the SoA (§3.3). This is still far better than
today's three scattered edit sites (L5), and it is pointer-free.

Kinds `wp/ik`, `_wp` literals, `implicit none`, `error stop`, ≤132 columns, lowercase free-form,
nvfortran-safe (no array-valued function result passed straight into a call — issue #7). This is a
**design/plan document only**; no source is modified here.

---

## 1. Scope, current state & limitations

### 1.1 What MEDS writes today

There is exactly **one** diagnostic path (`src/io/meds_io.f90`, `libmeds_io.a`, **always compiled** — a
hard dependency) over a thin `iso_c_binding` wrapper of the netCDF **C** library (`meds_netcdf_c.f90`).
The `-DMEDS_ENABLE_IO=OFF` stub build and its no-op twin (`meds_io_stub.f90`) **were removed by design**
(CLAUDE.md), so there is no stub to keep in sync. It has two write streams and a restart reader, all
driven from `meds_main.f90`:

- **Diagnostic timeseries** `<prefix>-D-output.nc` — `io_create` defines a fixed schema (UNLIMITED
  `time`, FIXED `cohort` = `cfg%io_cohort_max`, FIXED `patch` = `cfg%io_patch_max`); `io_write_snapshot`
  appends **one time record** of the instantaneous live SoA prefix `[1:n]` (10 per-cohort vars, 6
  per-patch vars, 5 per-site totals, calendar); `io_close` closes. Cohort/patch slabs are `1×max`
  chunked + deflated (L4, shuffle) so the unused fill tail compresses. `error stop` if `n > cap`.
- **State checkpoint / restart** `<prefix>-S-<YYYYMMDDHHMMSS>.nc` — `io_write_state` writes
  **prognostic state only** (no diagnostics; geometry re-derived from `dbh` on read); `io_read_state`
  restarts. This stream is **orthogonal to diagnostics and stays exactly as-is** (§5.4).

### 1.2 The limitations this design removes

| # | Limitation (today) | Consequence |
|---|---|---|
| L1 | **Annual cadence only.** Writes fire only at `is_new_year` in `meds_main` (`mod(iyear,interval)==0`). | Finest frequency = 1 year, though the engine steps daily/monthly. No monthly/daily/fast output. |
| L2 | **Instantaneous only.** No integration/averaging anywhere; the module header states "no buffering beyond the live state is needed." | No time-mean, sum, min/max, or variance diagnostic is possible. |
| L3 | **One frequency for the whole stream.** A single `io_output_interval_years` governs every variable. | Cannot write different variables at different cadences. |
| L4 | **No per-variable toggle.** Every var is defined + written unconditionally; only knob is `io_write_output` on/off. | All-or-nothing; no config-selectable subset. |
| L5 | **Hard-coded variable list.** Adding one var = a 3-site Fortran edit (a var-id field in `meds_io_t`, a `def_*` in `io_create`, a `put_vara_*` in `io_write_snapshot`). | No data-driven registry. |
| L6 | **One lump file** grows for the whole run (`-D-output.nc`). | A century run is one ever-growing file; no per-year/-month split. |
| L7 | **Fixed cohort/patch dims** at config caps; live `n > cap` is a hard `error stop`. | Caps must be sized a priori; ragged growth across a window unhandled. |
| L8 | **TOML can't express per-variable config.** Flat `key=value`, real arrays only — no string arrays, no arrays-of-tables. | A variable-toggle list / per-variable frequency table is not representable. |

The single choke points are `io_create` (schema) + `io_write_snapshot` (record write) in `meds_io.f90`,
the `is_new_year`-gated block in `meds_main.f90` (~lines 126–137, the cadence owner), and the `[io]`
keys in `meds_config.f90` + `meds_config_io.f90`. This design replaces all three.

### 1.3 Design goals & non-goals

**Goals.** (G1) A **DRY variable registry**: one registration list feeds all streams; adding a variable
is a registry `add_variable()` line plus a `extract_variable` case (two adjacent, co-located edits —
§3.3 — replacing today's three scattered sites). (G2) **Multi-frequency in-memory aggregation**
(fast/daily/monthly/annual) with mean/sum/min/max/last/mean-of-squares operators **plus a dt-weighted
flux integrator** for rate quantities (§3.5, §4.3). (G3) **High-level flux-group toggles + per-variable
overrides** from TOML, with unknown-key validation (the fine-grained overrides live in an optional
`meds_io_config.toml`, §6). (G4) **Time-chunked per-stream files**, chunked+deflated, never one lump.
(G5) **I/O off the compute critical path** (in-memory integrate + buffered flush; optional P2 async on
OpenMP-enabled builds only). (G6) **CF-compliant averaged output**: `cell_methods` + `time_bnds` so a
reader can tell a mean from an instantaneous value.

**Non-goals (this design).** Multi-site/polygon output (single-site layout retained; §10 Q5 notes the
seam). Exposing the fast loop's per-sub-step diagnostics to the registry: the fast biophysics loop is
now wired into `meds_main` (opt-in `[fast].fast_biophysics_on`), but the **fast** output tier's
consumption of its sub-`dt_fast` CO₂/flux/column diagnostics is the tier's P1 wiring (§3.5). The
met-forcing reader (it exists — `src/forcing`). The diel-cycle (`qmean`) tier. On-disk format stays
netCDF-4/HDF5 via the existing C bindings (no zarr).

---

## 2. Where it lives (library DAG)

The DAG is `shared ← {allometry, plant} ← state ← demography ← aux ← main`. The aggregation engine
splits cleanly along the **state/process wall** already used everywhere else, and — crucially — along a
**netCDF wall** drawn *inside* the output family. netCDF is now a hard dependency (there is no stub
build to keep honest), so the wall is retained purely for **DAG hygiene**: it keeps the per-step
integrate edge off netCDF, so the stepper never links a C dependency it does not use.

- **The registry + integrators are pure data + pure kernels, and are netCDF-FREE** — they name only
  `site_t` (read-only) and plain arrays, exactly like `meds_demography_diagnostics` (which they call).
  `aux` (the stepper) references *only* these, for the per-step integrate tick. Because a Fortran static
  archive pulls only referenced objects, `aux` linking `libmeds_io.a` but calling only the integrator
  symbols **does not pull the serializer object, and so does not pull netCDF** — the one new edge
  (`aux → io-integ`) carries no C dependency.
- **The netCDF serialization is quarantined in `meds_output_stream`** (over `meds_netcdf_c`) —
  unchanged in *mechanism*, generalized in *content*. It is referenced **only from `meds_main`**, which
  drives every disk flush. The stepper never calls it.

**The stepper integrates; main serializes (the resolution of the §2↔§4.5 seam).** `advance_one_step`
calls the netCDF-free `output_integrate` each step (fold live state into buffers; at a period roll-over,
normalize + chain + reset + **stage** the closed record into a plain-array pending queue on the manager —
all buffer arithmetic, zero disk). After `advance_one_step` returns, `meds_main` calls
`output_serialize_pending(mgr)` — the **only** netCDF-touching entry. So the flush lives at `main`
(honouring "`io` is a leaf `main` links"), the stepper stays netCDF-free, and `output_manager_t` is a
**netCDF-free plain-data type** (it holds `stream_file_t` handles whose `ncid`s are just integers; only
the *operations* on them touch C, and those live in the serializer).

```
shared ─┬─ allometry ─ state ─ demography ─┬─ aux (stepper: ticks the netCDF-FREE integrators each step;
        │                                  │       normalize+chain+reset+STAGE at period roll-over; NO disk)
        │                                  │         │ links only meds_output_{types,registry,integrate}
        ├─ plant                           ├─ io  (meds_output_types     ← var_desc_t, integ_buffer_t, output_manager_t [netCDF-free])
        │                                  │      (meds_output_registry  ← the registration list + TOML overrides + freq index)
        └─ biophysics                      │      (meds_output_integrate  ← extract_variable + integrate/normalize/reset + the output_integrate tick [netCDF-free])
                                           │      (meds_output_stream     ← per-stream netCDF file lifecycle  ← the netCDF wall)
                                           │      (meds_output_manager    ← output_serialize_pending drains the stage → meds_output_stream [MAIN-only])
                                           │      (meds_netcdf_c          ← REUSED verbatim: the C bindings)
                                           │      (meds_io                ← state/restart, UNCHANGED)
                                           └─ main (owns output_manager_t; drives output_serialize_pending → the ONLY flush)
```

**Files & CMake.** `CMakeLists.txt` globs `src/io/*.f90` into `libmeds_io.a`, so new modules auto-add.
The whole family is **always compiled** — netCDF is a hard dependency, so there is no `MEDS_ENABLE_IO`
switch and no stub twin. The netCDF-free registry + integrator modules and the netCDF-quarantined
serializer are still separate objects (the DAG-hygiene split above), but both always build.

| File | Role | Analogue (ED2 / MEDS) |
|---|---|---|
| `src/io/meds_output_types.f90` (new) | `var_desc_t`, `stream_desc_t`, `integ_buffer_t`, `output_registry_t`, `output_manager_t`; the io-only `AGG_*` / `DIM_*` codes (`use`s `FREQ_*`/`GRP_*`/`N_FREQ`/`FC_*` from `src/shared/meds_output_config`, §6.4) | ED2 `var_table` + `idim_type`; MEDS `meds_biophysics_types` `SOIL_*`/`ENERGY_*` block |
| `src/shared/meds_output_config.f90` (new) | `output_config_t` (the `[output]` block) + the codes the config needs at the DAG root: `FREQ_*`, `N_FREQ`, `GRP_*`, `FC_*` | MEDS `meds_forcing_config` (the exact same anti-back-edge pattern) |
| `src/io/meds_output_registry.f90` (new) | `build_output_registry(cfg) → output_registry_t`: the ONE registration list; applies TOML enable flags | ED2 the `vtable_edio_*` call sites; MEDS `build_soil_params` |
| `src/io/meds_output_integrate.f90` (new) | `pure`/`elemental` integrate/normalize/reset kernels over `integ_buffer_t`; `extract_site(site, registry) → raw per-var slabs` (loops `extract_variable`); **`output_integrate(mgr, site, cadence)` — the per-step, netCDF-free tick the stepper calls** (folds extracted values into `mgr`'s buffers, normalize+chain+reset+stage at roll-over; touches NO netCDF, so `aux` linking it stays off netCDF) | ED2 `integrate_/normalize_/zero_ed_*mean_vars` (~8700 lines → one generic per-var loop) |
| `src/io/meds_output_stream.f90` (new) | per-stream netCDF file lifecycle: `stream_open_file`, `stream_write_record`, `stream_close_file`; time-chunk rollover | ED2 `h5_output`; MEDS `meds_io%io_create/io_write_snapshot/io_close` |
| `src/io/meds_output_manager.f90` (new) | **`output_serialize_pending(mgr)`** — drains `mgr`'s pending stage to disk via `meds_output_stream` (called by `main`; the serializer wall). Referenced **only from `main`**, so `aux` never pulls this object (and thus never pulls netCDF). The `output_manager_t` type itself is a netCDF-free plain-data record defined in `meds_output_types` (registry + integ buffers + `stream_file_t` handles + the pending stage) | ED2 `ed_output`; new |
| `src/io/meds_io.f90` (unchanged) | state/restart streams | — |
| `test/test_output_integrate.f90`, `test/test_output_registry.f90`, `test/test_output_roundtrip.f90` (new) | CTest | `test_column_hydrology` etc. |

**Build nvfortran multicore on every new module** — a green ifx suite is not sufficient (CLAUDE.md
issue #7). The integrate kernels are plain-array arithmetic (a natural fit for the same
`target teams distribute` discipline as `growth_step`, though P0 keeps them host-only — §4.5).

---

## 3. The variable registry (`var_desc_t`) — a modern `var_table`

### 3.1 The descriptor

Every diagnostic variable is one **`var_desc_t`** — a pure DATA descriptor, no pointers. The registry is
a fixed, source-defined **list** of these (the "registration list"), built once at start-up and then
**immutable**; the TOML surface only flips the `enabled` flag and the per-frequency membership mask.

```fortran
! io-only selector codes (in meds_output_types, beside the type defs) -----------------------!
! (FREQ_*, N_FREQ, GRP_*, FC_* live in src/shared/meds_output_config with output_config_t -- §6.4 --
!  so the DAG-root config type carries no io back-edge; meds_output_types `use`s them.)
integer(ik), parameter :: AGG_MEAN = 1_ik, AGG_SUM = 2_ik, AGG_MIN = 3_ik, AGG_MAX = 4_ik,   &
                          AGG_LAST = 5_ik, AGG_MEANSQ = 6_ik,    & ! MEANSQ pairs with MEAN -> variance
                          AGG_TMEAN = 7_ik, AGG_FLUXSUM = 8_ik     ! dt-weighted: state time-mean / flux integral

integer(ik), parameter :: DIM_SCALAR = 0_ik, DIM_COHORT = 1_ik, DIM_PATCH = 2_ik,            &
                          DIM_SOIL = 3_ik,  DIM_PFT = 4_ik        ! trailing axis
! demographic dims DIM_COHORT/DIM_PATCH change ONLY at the monthly (fiss/fuse/recruit/cull) and annual
! (disturbance) cadence, so within any output window (<= 1 month, §4.4) their slot set is FIXED -- they
! index directly like the stable axes DIM_SOIL/DIM_PFT/DIM_SCALAR. No id-keying is needed anywhere.

! variable GROUP + stream/frequency codes -> declared in src/shared/meds_output_config (§6.4), because
! output_config_t at the DAG root needs them; meds_output_types `use`s them. Shown here for reference.
integer(ik), parameter :: GRP_STRUCTURE = 1_ik, & ! demography/size: nplant, dbh, agb, height, ids, patch geometry
                          GRP_CARBON    = 2_ik, & ! carbon fluxes: gpp, npp, Rh, NEE
                          GRP_WATER     = 3_ik, & ! water fluxes: transpiration, soil moisture, evap
                          GRP_ENERGY    = 4_ik    ! energy fluxes/state: soil/leaf temp, sensible/latent

! stream/frequency bit positions -> a single membership bitmask per variable ---------------!
integer(ik), parameter :: FREQ_FAST = 1_ik, FREQ_DAILY = 2_ik, FREQ_MONTHLY = 4_ik,          &
                          FREQ_ANNUAL = 8_ik      ! ior() these; test with iand(). FAST = a multiple of dt_fast (§4.1)
integer(ik), parameter :: N_FREQ = 4_ik          ! number of temporal tiers

type :: var_desc_t
   character(len=32)  :: name       = ''          ! netCDF variable name (e.g. 'agb')
   character(len=96)  :: long_name  = ''          ! CF long_name attribute
   character(len=24)  :: units      = ''          ! CF units attribute
   integer(ik)        :: dim        = DIM_SCALAR   ! DIM_* : the trailing axis
   integer(ik)        :: agg        = AGG_MEAN     ! AGG_* : temporal reduction operator
   integer(ik)        :: group      = GRP_STRUCTURE! GRP_* : the flux-group a high-level toggle switches (§6)
   integer(ik)        :: streams    = 0_ik         ! ior(FREQ_*) membership bitmask
   logical            :: enabled    = .true.       ! master per-variable on/off (TOML-mutable)
   integer(ik)        :: xtype      = NC_DOUBLE     ! NC_DOUBLE | NC_INT (C-binding constraint)
   integer(ik)        :: source_id  = 0_ik          ! which SoA field / reduction supplies it (§3.3)
end type var_desc_t
```

`streams` is the direct analogue of ED2's per-variable `ihist/ianal/idail/imont/...` 0/1 flags,
compressed into one integer: a variable "belongs to the daily and annual streams" ⟺
`streams = ior(FREQ_DAILY, FREQ_ANNUAL)`; the daily flush writes variable `v` iff
`iand(v%streams, FREQ_DAILY) /= 0 .and. v%enabled`. One variable opts into *N* frequencies by *N* bits
— adding a frequency is one bit + one clause, exactly ED2's decoupling.

**One membership rule the registry enforces: a `DIM_COHORT`/`DIM_PATCH` variable may NOT declare
`FREQ_ANNUAL`.** Annual output is **site-level only** (`DIM_SCALAR` and the fixed site axes
`DIM_SOIL`/`DIM_PFT`), because a demographic window longer than a month would straddle the annual
disturbance restructuring (§4.4, §5.1). `build_output_registry` `error stop`s if a cohort/patch variable
carries the annual bit — this is what keeps the annual stream free of any ragged demographic dimension.

`dim` is MEDS's clean replacement for ED2's magic `idim_type` integer (a ~40-case header comment): an
enumerated axis code, not an arithmetic tuple. `geth5dims`'s job (idim_type → rank+shape) becomes a
trivial `select case(v%dim)` in the stream writer (§5.3).

**State vs. flux — the `agg` operator must respect physical units (ED2's integrated-vs-mean split).**
ED2 deliberately separates *integrated fluxes* from *mean states*, and MEDS must too. A **state**
(`agb_site`, `soil_temp_site`, `dbh_cohort`) is a stock: its period value is a **time-mean**, and if
steps are non-uniform the mean must be **dt-weighted** — `AGG_TMEAN` integrates `Σ(x·dt)` and normalizes
by `Σdt`, reducing to `AGG_MEAN` under uniform `dt`. A **flux** (`gpp_site` in `kgC/m²/s`, transpiration,
NEE) is a rate: a physically meaningful period *total* is the **integral** `Σ(x·dt)`, whose units are
the rate × time (`kgC/m²` over the period) — this is `AGG_FLUXSUM`. A plain `AGG_SUM` of per-step
samples is dt-unaware: under uniform `dt` it is off by exactly the step-count factor (wrong magnitude
*and* wrong units for a "total"), and across tiers with different sub-step counts (fast feeding daily)
a plain sum is not a coherent integral. `AGG_SUM` is therefore reserved for **count-like** integer
tallies (e.g. number of mortality events) where dt-weighting is meaningless. A physical **rate** may be
reported either way, and MEDS supports both: as a **mean rate** (`AGG_TMEAN`, dt-weighted, keeps the
`/s` unit — the P1 default for `gpp_site`/`nee_site`/`transp_site`, §3.5) *or* as a **period integral**
(`AGG_FLUXSUM`, drops the `/s`); a physical **stock** uses `AGG_TMEAN`/`AGG_MEAN`. What is *never*
allowed for a physical rate is plain `AGG_SUM` (dt-unaware). Each `var_desc_t` thus declares, via its
`agg`, how the quantity is reduced; `extract_variable`/`integrate` carry the step `dt` (§4.3). The
`units` string documents what the operator produces (a `_FLUXSUM` var's on-disk units drop the `/s`).

### 3.2 The registry container

```fortran
type :: output_registry_t
   type(var_desc_t), allocatable :: var(:)        ! the immutable registration list
   integer(ik)                   :: nvar = 0_ik
   ! per-frequency precomputed index lists (which vars are live in each tier) -> no re-scan per flush:
   integer(ik), allocatable :: idx_freq(:,:)      ! (max_vars_per_freq, N_FREQ): var indices enabled in tier
   integer(ik)              :: nidx(N_FREQ) = 0_ik
end type output_registry_t
```

`build_output_registry(cfg)` fills `var(:)` from the source registration list (§3.4), applies the TOML
`enabled` and frequency overrides (§6), then precomputes `idx_freq` so each flush walks only its own
live variables — the analogue of ED2's `varloop` bit-filter, done once up front instead of per write.

### 3.3 Binding a descriptor to live state — **no pointers**

ED2's fatal choice was `vt_vector(iptr)%var_rp => var` (a raw pointer into the SoA). MEDS **never caches
a pointer**. Instead each `var_desc_t` carries an integer `source_id` naming *which* reduction / SoA
field supplies its value, and a single **`extract_variable`** switchboard **copies** the current value out
of live state each step — the verb names exactly that (pull the value *out of* the SoA and copy it, a
snapshot, not a cached alias — the naming makes the load-bearing "read/copy, do not alias" divergence
from ED2 explicit). The integrators then fold that copied value into their running reductions.
`extract_site` loops the registry and calls `extract_variable` per variable:

```fortran
! meds_output_integrate.f90 -- pure extraction: current instantaneous value(s) of one variable.
subroutine extract_variable(site, v, scalar_out, slab_out, n_out)
   type(site_t),     intent(in)  :: site
   type(var_desc_t), intent(in)  :: v
   real(wp),         intent(out) :: scalar_out          ! for DIM_SCALAR
   real(wp),         intent(out) :: slab_out(:)         ! for DIM_COHORT/PATCH/... : the live [1:n] prefix
   integer(ik),      intent(out) :: n_out               ! valid length this step (the live cohort/patch count)
   ! select case (v%source_id): copy site%cohort%agb(1:n), or call total_agb(site), etc.
end subroutine
```

Because `extract_variable` **reads and copies** (never aliases), a cohort birth/death/fusion between steps
is invisible to the output layer — there is **no `filltab_alltypes` re-hash** and no
hidden-global-mutable-state. This is the load-bearing divergence from ED2 (Bundle-style takeaway: pass an
explicit descriptor, re-derive from state, do not cache pointers).

### 3.4 The DRY registration pattern

One source list — a single `contains`-local block in `build_output_registry` — declares every variable
once. A tiny `add_variable()` helper appends a `var_desc_t`; the group, frequency membership and default
operator are literal arguments. This is **one of the two** "add a variable" edit sites (the other is the
matching `extract_variable` case, §3.3): the `add_variable()` line names the variable and its metadata,
and the `SRC_*` case copies the right SoA field or calls the right reduction. The two are **deliberately
co-located** — the `SRC_*` parameter, the `add_variable()` line, and the `extract_variable` case are
grouped by dimension so a new cohort-dimensioned variable is a three-line diff in one region of two
adjacent files, never a hunt across `io_create` + `io_write_snapshot` + a struct field (today's L5). It
is **two coupled edits, not one** — the irreducible price of rejecting ED2's pointer registry — but
pointer-free and local:

```fortran
subroutine build_output_registry(reg, cfg)
   type(output_registry_t), intent(out) :: reg
   type(meds_config_t),     intent(in)  :: cfg
   integer(ik) :: k
   allocate(reg%var(MAX_OUTPUT_VARS)) ; reg%nvar = 0_ik
   ! name             long_name               units       dim         agg       group          streams(default) src
   call add_variable('agb_cohort','aboveground biomass', 'kgC/m2',  DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, DAY_MON, SRC_C_AGB)
   call add_variable('nplant_cohort','plant number dens','plant/m2',DIM_COHORT, AGG_LAST, GRP_STRUCTURE, DAY_MON, SRC_C_NPLANT)
   call add_variable('dbh_cohort','diameter at breast hgt','cm',    DIM_COHORT, AGG_MEAN, GRP_STRUCTURE, DAY_MON, SRC_C_DBH)
   call add_variable('agb_site', 'site aboveground biomass','kgC/m2',DIM_SCALAR,AGG_MEAN, GRP_STRUCTURE, DAY_MON_YR,SRC_S_AGB)
   call add_variable('lai_site', 'site leaf area index',  'm2/m2',  DIM_SCALAR, AGG_MEAN, GRP_STRUCTURE, DAY_MON_YR,SRC_S_LAI)
   call add_variable('gpp_site', 'site GPP',              'kgC/m2/s',DIM_SCALAR, AGG_TMEAN,GRP_CARBON,    MON,     SRC_S_GPP) ! P1, fast loop
   ! ... (full table in §3.5) ...
   call apply_group_toggles(reg, cfg)       ! disable whole GRP_* groups from the main config (§6)
   call apply_io_config(reg, cfg)           ! optional per-variable overrides from meds_io_config.toml (§6)
   call build_freq_index(reg)               ! precompute idx_freq (after the annual DIM_COHORT/PATCH guard)
contains
   subroutine add_variable(nm, ln, un, dm, ag, gr, st, src)
      ! ... reg%nvar = reg%nvar+1 ; reg%var(reg%nvar) = var_desc_t(nm, ln, un, dm, ag, gr, st, .true., ...)
   end subroutine
end subroutine
```

The `DAY_MON` (= D M), `MON` (= M), `DAY_MON_YR` (= D M Y), … are named `ior` combinations
(`integer(ik), parameter` — e.g. `DAY_MON_YR = ior(ior(FREQ_DAILY,FREQ_MONTHLY),FREQ_ANNUAL)`, which
deliberately **excludes** `FREQ_FAST` — the fast tier is opt-in) — a readable default membership per
variable, all overridable from TOML. Note the demographic (`DIM_COHORT`)
rows carry **no** `FREQ_ANNUAL` bit (`DAY_MON`, not `DAY_MON_YR`): annual output is site-level only
(§3.1), and the registry `error stop`s if a cohort/patch variable is given the annual bit.

### 3.5 Example variable-registry table (P0 set + P1 additions)

**Naming convention — every raw variable carries a scale suffix `_cohort` / `_patch` / `_site`** naming
the entity it describes (= its output dimension), so a quantity that exists at several scales is
unambiguous (`agb_cohort` = per-plant biomass on the cohort axis; `agb_site` = the stand total). Site
aggregates take `_site` (this **replaces** the old `total_` prefix). Fluxes drop any `_rate` suffix —
`gpp_site`, `nee_site`, `transp_site` — the unit string already says it is a rate. **Exempt**
(name already carries the entity, a scale suffix would be redundant/misleading): the identity, CSR and
count fields `global_cohort_id`, `global_patch_id`, `owner_patch`, `cohort_offset`, `cohort_count`,
`n_cohort`, `n_patch`.

Cohort/patch rows carry **no annual bit** (annual is site-level only, §3.1). The `group` column is what
the main config's high-level flux-group toggles switch on/off (§6).

| name | long_name | units | dim | agg | group | default streams | source |
|---|---|---|---|---|---|---|---|
| `nplant_cohort` | plant number density | plant/m² | cohort | last | structure | D M | `cohort%nplant` |
| `dbh_cohort` | diameter at breast height | cm | cohort | mean | structure | D M | `cohort%dbh` |
| `height_cohort` | height | m | cohort | mean | structure | D M | `cohort%height` |
| `basal_area_cohort` | basal area per plant | cm²/plant | cohort | mean | structure | M | `cohort%basal_area` |
| `agb_cohort` | aboveground biomass per plant | kgC/plant | cohort | mean | structure | D M | `cohort%agb` |
| `leaf_area_cohort` | leaf area per plant | m²/plant | cohort | mean | structure | M | `cohort%leaf_area` |
| `growth_avg_cohort` | moving-avg growth (mort. predictor) | cm/yr | cohort | mean | structure | M | `cohort%growth_avg` |
| `pft_cohort` | plant functional type index | – | cohort | last | structure | D M | `cohort%pft` (int) |
| `owner_patch` | owning patch (1-based) | – | cohort | last | structure | D M | `cohort%owner_patch` (int) |
| `global_cohort_id` | persistent cohort id | – | cohort | last | structure | D M | `cohort%global_id` (int) |
| `area_patch` | patch area fraction | – | patch | last | structure | D M | `patch%area` |
| `age_patch` | time since disturbance | yr | patch | last | structure | M | `patch%age` |
| `dist_type_patch` | disturbance type | – | patch | last | structure | M | `patch%dist_type` (int) |
| `cohort_offset` | first cohort of patch (CSR) | – | patch | last | structure | D M | `patch%cohort_offset` (int) |
| `cohort_count` | cohorts in patch (CSR) | – | patch | last | structure | D M | `patch%cohort_count` (int) |
| `global_patch_id` | persistent patch id | – | patch | last | structure | D M | `patch%global_id` (int) |
| `nplant_site` | site total plant number | plant/m² | scalar | mean | structure | D M Y | `total_nplant(site)` |
| `basal_area_site` | site total basal area | m²/m² | scalar | mean | structure | D M Y | `total_basal_area(site)` |
| `agb_site` | site aboveground biomass | kgC/m² | scalar | mean | structure | D M Y | `total_agb(site)` |
| `lai_site` | site leaf area index | m²/m² | scalar | mean | structure | D M Y | `total_lai(site)` |
| `n_cohort` | cohorts in use | – | scalar | last | structure | D M Y | `site%cohort%n` (int) |
| `n_patch` | patches in use | – | scalar | last | structure | D M Y | `site%patch%n` (int) |
| — *P1 additions (need the fast loop / new diagnostics)* — |||||||
| `gpp_site` | site gross primary productivity | kgC/m²/s | scalar | **tmean** | carbon | M | aggregate cohort GPP (dt-wtd) |
| `nee_site` | net ecosystem exchange | kgC/m²/s | scalar | **tmean** | carbon | M | site NEE (dt-wtd) |
| `transp_site` | site transpiration | mm/s | scalar | **tmean** | water | M | aggregate transpiration (dt-wtd) |
| `soil_temp_site` | soil temperature by layer | K | soil | tmean | energy | M | `patch%soil_e%soil_temp` |
| `soil_water_site` | volumetric soil moisture | m³/m³ | soil | tmean | water | M | `patch%soil_w%theta` |
| `mortality_cohort` | realized mortality intensity | 1/yr | cohort | tmean | structure | M | rate array |

(**last** for identity/index/CSR fields — an average of `global_id` or `cohort_offset` is meaningless;
they carry the *end-of-window* snapshot that interprets the fixed-within-window slabs and lets a reader
track a cohort *across* windows, §4.4. **tmean** is the dt-weighted state mean and is the default for
every physical stock — it collapses to `mean` under the uniform slow-step `dt` but is unit-correct once
the fast tier's `dt_fast` sub-steps feed it; **fluxsum** integrates a rate to a period total, dropping
the `/s` in the on-disk `units`; **sum** is reserved for integer count tallies, §3.5-note above.)

The P1 flux/state additions default to the **monthly mean** (`M`, `tmean`) — a single mean, **no
variance** (the `AGG_MEANSQ` variance-companion operator stays defined but no registry variable uses it
for now). A user who wants a flux on the fast or daily tier adds it in `meds_io_config.toml` (e.g.
`gpp_site = "F D M"`, §6.3), so the fast tier is populated on demand rather than by a heavy default.

### 3.6 Instantaneous escape hatch

A variable whose default streams include no averaging tier and whose `agg = AGG_LAST` on a per-step
frequency reproduces today's *values* exactly (instantaneous, written at flush). So the P0 registry with
the per-cohort fields at `AGG_LAST` on a **daily or monthly** stream carries the **same field values**
as the current `-D-output.nc` (today's per-cohort snapshot is likewise a demographic slab). The
migration is a strict generalization. (The *annual* stream can no longer host per-cohort slabs — it is
site-level only, §3.1 — so the today-equivalent per-cohort snapshot lives on the daily/monthly stream.)

**Row order now matches, because the within-window slot set is fixed.** Since cohort/patch structure
changes only at month/year boundaries (§4.4), an `AGG_LAST` cohort record is flushed in **live-SoA slot
order** — the *same* order `io_write_snapshot` writes today — so there is no first-seen-`global_id`
permutation to reconcile. The files still differ in **metadata** (the added per-variable
`cell_methods`/`_FillValue` attributes, §5.3), so the round-trip test (§8, test 5) compares values
**up to metadata** (the `global_id` slab remains the identity key for a reader), which is the honest and
testable claim.

---

## 4. The aggregation engine

### 4.1 Temporal tiers and the chaining rule

Four integration tiers, mirroring ED2's `fmean → dmean → mmean` chain (plus a yearly tier) but
**cadence-driven by the existing stepper flags** rather than second-counting. The finest tier is the
**fast** tier (ED2's `fmean`): like ED2, it is **not** an "hourly" tier — its integration window is a
**user-set multiple of `dt_fast`** (`[output.fast].interval_steps`, §6), so it integrates every `dt_fast`
sub-step and flushes every `interval_steps × dt_fast` (e.g. 4 × 15 min = one record/hour, or `1` for a
record every fast step). It exists only when the fast biophysics loop runs (`[fast].fast_biophysics_on`);
without it the daily tier is the finest.

| Tier | Fed from (per variable) | Integrated every | Normalized & flushed at | ED2 twin |
|---|---|---|---|---|
| **fast** (`fmean`) | raw state (if finest active) | `dt_fast` sub-step | fast-output boundary (`interval_steps × dt_fast`) | `fmean` |
| **daily** (`dmean`) | raw state, *or* the fast mean if the fast tier is active | slow step / fast boundary | day roll-over | `dmean` |
| **monthly** (`mmean`) | raw state, *or* the next-finer active tier | slow step / day roll-over | month roll-over | `mmean` |
| **annual** (`ymean`) | raw state, *or* the next-finer active tier | slow step / month roll-over | year roll-over | yearly vars |

**The feeder rule, stated concretely (per variable, not per tier).** For each variable, precompute its
ordered list of active tiers from `streams ∧ enabled` (call the finest `t₀`). **The finest active tier
`t₀` integrates from raw state every step** via `extract_variable` — this is the load-bearing correction:
*whatever* the finest active tier is, it reads state directly, so a **monthly-only** or **annual-only**
config works (the monthly buffer integrates every slow step, not off a non-existent daily feeder).
**Every coarser active tier chains from the next-finer active tier** at that finer tier's roll-over —
never from raw state, and never from an *inactive* intermediate tier. So daily+annual (no monthly)
chains annual off daily; monthly-only reduces raw state straight into the monthly buffer. The per-var
"who feeds whom" is stored once in `idx_freq`/a `feeder(tier)` map at registry build (§3.2), so
`output_integrate` folds raw state into the **precomputed finest tier per variable**, never a
hard-coded daily. **Configs to test explicitly: monthly-only, annual-only, and daily+annual-skip-monthly**
(§8, test 2).

**Chaining exactness.** For `AGG_MEAN`/`AGG_TMEAN` the chain is exact (a dt-weighted mean of dt-weighted
means, carrying `Σdt` as the chained weight, equals the direct dt-weighted mean over the coarse period).
For `AGG_FLUXSUM` the chain is a plain sum of sub-period integrals (integrals add). For
`AGG_SUM/MIN/MAX/LAST` the operator re-applies up the chain (sum of sums, min of mins, …), exact.
**`AGG_MEANSQ` chains on the lower tier's *mean*, not its meansq** — variance is
`mmean(x²)_from_finer_means − (mmean(x))²`, the diel/inter-day variance, matching ED2's `mmsqu`
second-moment tier. (The operator stays defined for this chaining, but variance is **deferred** — no
registry variable emits it yet, §3.5.) Each chained buffer carries its weight `Σdt` (not a bare sample
count) so the roll-up stays dt-correct even when sub-tiers have unequal step counts.

### 4.2 The integrator buffer

Per active `(variable, tier)` pair, one `integ_buffer_t` holds the running reduction, the sample count,
and the **dt weight**. There is a single **fixed-index** addressing mode for every dim: `DIM_SOIL`/
`DIM_PFT` are stable site axes, and — the key simplification (§4.4) — `DIM_COHORT`/`DIM_PATCH` are
**also fixed within any output window** (the slot set only changes at month/year boundaries, which are
also the window boundaries), so they integrate by direct slot index `[1:n]` too. **No `global_id`
keying, no per-row id map, no first-seen ordering.** The slab is dimensioned to the relevant cap
(`io_cohort_max`, `io_patch_max`, `n_soil_layer`, `n_pft`).

```fortran
type :: integ_buffer_t
   integer(ik) :: var_id  = 0_ik          ! index into registry%var(:)
   integer(ik) :: freq    = 0_ik          ! FREQ_* this buffer serves
   integer(ik) :: agg     = AGG_MEAN
   integer(ik) :: dim     = DIM_SCALAR
   real(wp)    :: scal    = 0.0_wp        ! scalar integrator (DIM_SCALAR)
   real(wp)    :: scal2   = 0.0_wp        ! scalar second moment  (AGG_MEANSQ, DIM_SCALAR)
   real(wp)    :: wsum    = 0.0_wp        ! scalar dt weight Sum(dt) (AGG_TMEAN/FLUXSUM/MEANSQ); nsamp for count-aggs
   real(wp),    allocatable :: slab(:)    ! per-index integrator, sized to cap (slot for cohort/patch, layer for soil, bin for pft)
   real(wp),    allocatable :: slab2(:)   ! second-moment integrator (AGG_MEANSQ only)
   real(wp),    allocatable :: wsum_slab(:)! per-index dt weight Sum(dt) (slab TMEAN/FLUXSUM)
   integer(ik), allocatable :: hits(:)    ! per-index sample count
   integer(ik) :: n_slab  = 0_ik          ! number of live indices this window ([1:n] cohort/patch count, or the fixed axis length)
   integer(ik) :: nsamp   = 0_ik          ! scalar sample count
   real(wp)    :: seed    = 0.0_wp        ! MIN=+huge, MAX=-huge, SUM/MEAN/FLUXSUM=0, LAST=undefined
end type integ_buffer_t
```

`scal2`/`slab2` carry the second moment for `AGG_MEANSQ`; `wsum`/`wsum_slab` carry the `Σdt` weight that
makes `AGG_TMEAN`/`AGG_FLUXSUM` and the chained roll-up dt-correct (§4.1, §4.3). Dropping the id map
(`id_of_row`/`n_active`/`id_keyed`) is the direct payoff of the fixed-within-window slot set (§4.4).

### 4.3 The three kernels (integrate / normalize / reset)

A **single generic per-variable loop** replaces ED2's ~8700 lines of hand-written per-field boilerplate.
`integrate_scalar` folds one extracted value into the buffer by operator; `normalize_buffer` closes the
period; `reset_buffer` re-seeds it. All `pure`/`elemental`-friendly (arithmetic + intrinsics):

```fortran
elemental subroutine integrate_scalar(buf, x, dt)      ! one entity, one step of length dt [s]
   type(integ_buffer_t), intent(inout) :: buf
   real(wp),             intent(in)    :: x, dt
   select case (buf%agg)
   case (AGG_MEAN);          buf%scal = buf%scal + x          ; buf%nsamp = buf%nsamp + 1_ik   ! equal-weight
   case (AGG_TMEAN);         buf%scal = buf%scal + x*dt       ; buf%wsum  = buf%wsum + dt       ! dt-weighted state
   case (AGG_FLUXSUM);       buf%scal = buf%scal + x*dt       ; buf%wsum  = buf%wsum + dt       ! integral of a rate
   case (AGG_SUM);           buf%scal = buf%scal + x          ; buf%nsamp = buf%nsamp + 1_ik   ! count-like tally
   case (AGG_MIN);           buf%scal = min(buf%scal, x)      ; buf%nsamp = buf%nsamp + 1_ik
   case (AGG_MAX);           buf%scal = max(buf%scal, x)      ; buf%nsamp = buf%nsamp + 1_ik
   case (AGG_LAST);          buf%scal = x                     ; buf%nsamp = buf%nsamp + 1_ik
   case (AGG_MEANSQ);        buf%scal = buf%scal + x*dt ; buf%scal2 = buf%scal2 + x*x*dt ; buf%wsum = buf%wsum + dt
   end select
end subroutine

subroutine normalize_buffer(buf, out, valid)            ! close the period; valid=.false. -> emit _FillValue
   ! empty-period guard (nsamp==0 .and. wsum==0): valid = .false., out = MISSING; NO divide-by-zero, NO +/-huge leak.
   ! AGG_MEAN            -> scal / nsamp                 (nsamp>0)
   ! AGG_TMEAN           -> scal / wsum                  (wsum >0)   dt-weighted state mean
   ! AGG_FLUXSUM         -> scal  (= Sum(x*dt))          as-is       period integral; units drop the /s
   ! AGG_SUM/MIN/MAX/LAST-> scal  as-is                              (MIN/MAX still-seeded -> valid=.false.)
   ! AGG_MEANSQ          -> mean = scal/wsum ; var = max(0, scal2/wsum - mean*mean)  (both to companion vars)
end subroutine

subroutine reset_buffer(buf)                            ! re-seed for the next period
   ! MEAN/TMEAN/FLUXSUM/SUM/MEANSQ: scal=scal2=0, wsum=0, nsamp=0 ; MIN: +huge ; MAX: -huge ; LAST: keep value
   ! slab: slab=slab2=0, hits=0, wsum_slab=0 (a fresh window over the current, fixed slot set)
end subroutine
```

**Empty-period normalization is defined, not left to chance.** `normalize_buffer` returns a `valid`
flag; when a tier saw **zero samples** in a period (the fast tier enabled but the fast loop didn't run,
or a slab index past the live count), `valid=.false.` and the serializer writes the CF `_FillValue` for
that record/slot — never a `0/0` NaN, and never the raw `±huge` MIN/MAX seed. The unused `[n+1:cap]`
slab fill tail is likewise written as `_FillValue`, not zero (§5.3). Test 1 covers the zero-sample
scalar and the never-occupied slab slot.

The **slab** variant (`integrate_slab`, every non-scalar dim) loops the live `[1:n]` prefix and folds
each index directly into `slab(i)`, incrementing that index's `hits` and `wsum_slab`; `normalize`
divides `AGG_MEAN` by the slot's `hits` and `AGG_TMEAN` by its `wsum_slab`. For `DIM_COHORT`/`DIM_PATCH`
the slot set is fixed for the whole window (§4.4), so every live slot has the full step count; for the
stable `DIM_SOIL`/`DIM_PFT` axes the same holds by construction. There is no id map, no first-seen
ordering, and no across-a-window cap-overflow path anywhere.

### 4.4 Cohort/patch dims are fixed within a window — the simplification

Cohorts and patches are born, fused, split, and culled — but **only at the monthly (recruit / cohort
fuse-terminate-split / sort) and annual (disturbance, patch fuse/terminate) cadence**, never mid-window:
between those boundaries `growth` and `mortality` change per-slot *values* (`dbh`, `nplant`) but **add or
remove no slot and never reorder** (CLAUDE.md order-of-operations). Combined with the rule that
cohort/patch output tiers never exceed monthly (§3.1) and the manager aligns each cohort/patch window to
the fiss/fuse boundary, **every output window lies entirely within one inter-restructuring interval, so
the SoA slot set `[1:n]` is constant across the window.** That single fact collapses the entire "ragged
demographic average" problem:

- **Integrate by direct slot index**, exactly like today's instantaneous snapshot — `slab(i) +=
  value(i)` for `i in [1:n]`. No `global_id` hash, no first-seen ordering, no per-id membership map.
- **The compact CSR (`cohort_offset`/`cohort_count`/`owner_patch`) written `AGG_LAST` describes the
  slabs exactly as today's reader expects** — there is no averaged-vs-instantaneous membership mismatch,
  because the set that was averaged *is* the set present at flush.
- **The persistent `global_id` slab is still written** (per §3.5), but only so a reader can track a
  cohort **across** windows (where the set legitimately changes at a boundary); it is **not** an
  integration key.
- **The cap check is the same simple instantaneous guard as today** (`n > cap`), never a separate
  "distinct ids across a window exceeded the cap" churn problem — that path is gone.

**Window/boundary alignment (the one wiring contract).** On `is_new_month`, the manager
**normalizes + stages + resets** the cohort/patch monthly buffers against the still-current slot set
**before** that month-boundary's fiss/fuse restructures it; the fresh window then integrates the new
month's (new) slot set. Equivalently: fiss/fuse defines the *start* of a window, never a point inside
one. (Finer tiers — fast, daily — sit strictly inside a month, so they never straddle a boundary at
all.) The old §4.4 machinery — id-keyed slabs, the averaged per-id `owner_patch` map, the
`ragged_avg_alive_only` flag, the mid-window fusion value-discontinuity, and the across-a-window
cap-overflow policy — is **deleted**: none of it can arise when the window cannot straddle a
restructuring event. `integrate_slab` folds each live index `i` in `[1:n]` directly into `slab(i)`,
`hits(i)`/`wsum_slab(i)` count that slot's contributions, and — because the slot set is fixed — every
live slot integrates the full window's step count.

**Fusion cleanliness comes for free.** Because fusion (which re-derives a survivor's `dbh` from conserved
AGB) happens only at the month boundary — a *window* boundary, never a point inside a window — no
`AGG_TMEAN` ever blends a pre- and post-fusion composition under one slot. The extensive/intensive
distinction that the old design had to flag as a limitation simply does not arise; a monthly mean is a
mean over one stable cohort composition, and the next month's fusion starts a fresh window.

**Fixed site axes are the same code path.** `DIM_SOIL`/`DIM_PFT` are stable by construction (no
birth/death), so they were always going to fold by direct index — and now `DIM_COHORT`/`DIM_PATCH` do
too. One `integrate_slab` serves every non-scalar dim; the only per-dim difference is which cap sizes
the slab and whether a CSR/`global_id` companion is written (cohort/patch: yes; soil/pft: no).

This "fixed slot set within a window, direct-index integrate, today-CSR at flush" scheme is what lets a
demographic state be time-averaged **without** ED2's per-birth/death pointer re-hash *or* the earlier
draft's id-keyed membership map — the monthly restructuring cadence is exactly the window cadence, so the
two never collide.

### 4.5 Coupling to the stepper cadence

The integrators tick from **`meds_stepper`** (the cadence owner), fed the same `is_new_month`/
`is_new_year` flags `meds_main` already computes, plus a new `is_new_day` (trivially `now%day /=
prev%day`) and (P1) `is_fast_out` (the fast-output boundary — the model crossed a multiple of `dt_fast`
equal to `interval_steps`, §4.1), and the step length `dt` (for dt-weighting, §4.3). The stepper's
per-step entry is **netCDF-free**: it integrates and *stages* closed records; it does **not** write.

```fortran
subroutine output_integrate(mgr, site, now, dt, is_fast_out, is_new_day, is_new_month, is_new_year)
   type(output_manager_t), intent(inout) :: mgr           ! netCDF-free plain data (buffers + pending stage)
   type(site_t),           intent(in)    :: site
   type(meds_time_t),      intent(in)    :: now
   real(wp),               intent(in)    :: dt            ! [s] length of the step just taken (dt_fast in the fast loop, dt_slow otherwise)
   logical,                intent(in)    :: is_fast_out, is_new_day, is_new_month, is_new_year
   ! 1. integrate: for each var, fold extract_variable(site,v) into its PRECOMPUTED FINEST active tier
   !                (mgr%feeder(:) from §4.1 — NOT hard-coded to daily; carries dt for TMEAN/FLUXSUM)
   ! 2. on is_fast_out  -> normalize+STAGE fast   records; chain fast-mean into the daily buffer; reset fast
   ! 3. on is_new_day   -> normalize+STAGE daily  records; chain into monthly (if monthly active) else annual; reset daily
   ! 4. on is_new_month -> normalize+STAGE monthly records BEFORE that boundary's cohort fiss/fuse (§4.4); chain into annual; reset monthly
   ! 5. on is_new_year  -> normalize+STAGE annual  records (site-level only, §3.1); reset annual
   ! all of the above is plain-array arithmetic on mgr; the STAGED records sit in mgr%pending until main drains them
end subroutine
```

Then, back in `meds_main`, once per step:

```fortran
call output_serialize_pending(mgr)   ! the ONLY netCDF-touching call
```

`output_serialize_pending` walks `mgr%pending`, opens/rolls the right per-stream file (§5.2), writes each
staged record, and clears the stage. Chaining feeds the **coarsest-below** *active* tier (the §4.1 feeder
map), so daily→annual-skip-monthly and monthly-only both route correctly. The roll-over tests fire **at
the boundary**, closing the period that just ended (ED2 names files by the *previous* period; §5.2), and
ordering is strict low→high tier so each chain reads a freshly-normalized feeder (the ED2 `ed_output`
ordering, distilled). P0 keeps the integrate step **host-only** (a cheap copy-reduce over the live
prefix, negligible beside the daily kernels); if profiling ever shows it matters, the slab integrate is
a `target teams distribute` sweep over plain arrays like `growth_step` — designed-in, not needed at P0.

**Where it wires:** `advance_one_step` gains an optional `mgr` argument. The **fast** tier integrates
inside the fast biophysics sub-step loop (every `dt_fast`, so `is_fast_out` and the fast `dt` come from
there); the coarser tiers integrate once per slow step. **The monthly cohort/patch flush must run
before that month boundary's fiss/fuse restructuring** (§4.4) — i.e. the monthly `normalize+stage+reset`
of the cohort/patch buffers precedes `vegetation_dynamics`' monthly block — so a window never straddles a
slot-set change. `meds_main` calls `output_serialize_pending` after `advance_one_step` returns — the
clean split of §2 that keeps `aux` netCDF-free and the flush at `main`.

---

## 5. File organization & performance

### 5.1 One file per stream per time-chunk (goal G4, FILE SIZE)

Never one lump. Each **stream** (frequency) writes to its own family of files, **time-chunked** so no
single file grows unbounded:

```
<output_dir>/<prefix>-F-<YYYYMMDD>.nc     fast    stream, one file per SIM-DAY   (interval-dependent record count)
<output_dir>/<prefix>-D-<YYYYMM>.nc       daily   stream, one file per SIM-MONTH (28-31 records)
<output_dir>/<prefix>-M-<YYYY>.nc         monthly stream, one file per SIM-YEAR  (12 records)
<output_dir>/<prefix>-Y.nc                annual  stream, ONE file for the WHOLE RUN (1 record/yr, site-level only)
<output_dir>/<prefix>-S-<YYYYMMDDHHMMSS>.nc   STATE/restart — UNCHANGED, one file per checkpoint
```

The stream letter (`F/D/M/Y`) mirrors ED2's fast-analysis-through-yearly convention; for the sub-annual
streams the date-stamp granularity is one tier **coarser** than the record cadence (a per-day-of-records
file for the fast stream, a per-month-of-records file for the daily stream, …), so each file holds a
tidy, bounded number of records. The **fast** stream's records-per-day follows `interval_steps` (e.g.
`4` fast steps → 24 records/day; `1` → one per `dt_fast`). The **annual** stream is different: it is now
**site-level only** (no cohort/patch slabs, §3.1) and holds one small record per year, so its default is
**one file for the whole run** — there is no decade bucket. (The old `-Y-<YYYY0>.nc` stamp zeroed the
last year digit to bucket a decade; that trick is dropped as needless — an annual file of a few hundred
scalar records is trivially small. If a user *does* want the annual stream split, `file_chunk` still
buckets it by a plain start year `-Y-<YYYY>.nc`, never a zero-padded decade.) The **time-chunk size is a
per-stream TOML knob** (`file_chunk`, §6). Each file has its own UNLIMITED `time` dim; records append
until the chunk boundary, then a new file opens. This is ED2's multi-record splitting, made a plain
rollover on the calendar flags.

### 5.2 Naming & rollover

`stream_write_record` compares the record's calendar bucket to the open file's bucket; on mismatch it
closes the current file and opens the next (`makefnam`-analogue: `<prefix>-<letter>[-<stamp>].nc`, stamp
granularity from `file_chunk`; the whole-run annual default carries no stamp at all). Files are named for
the **period whose data they hold** (the daily stream's `-D-202601.nc` holds January's daily means).
Rollover is driven purely by `now`, so it is exact across leap years (reusing `meds_time`).

### 5.3 netCDF layout per file

Each stream file is the **generalization of today's `io_create`**: UNLIMITED `time`; fixed ragged dims
`cohort`/`patch`/`soil`/`pft` sized to caps; per-variable definition driven by the registry loop instead
of the hard-coded list. `dim → (dimids, chunk)` is one `select case(v%dim)` (the clean `geth5dims`):

| `v%dim` | netCDF dims | chunk shape | deflate |
|---|---|---|---|
| `DIM_SCALAR` | `(time)` | default | L4+shuffle |
| `DIM_COHORT` | `(time, cohort)` | `(1, cohort_max)` | L4+shuffle |
| `DIM_PATCH` | `(time, patch)` | `(1, patch_max)` | L4+shuffle |
| `DIM_SOIL` | `(time, soil)` | `(1, n_soil)` | L4+shuffle |
| `DIM_PFT` | `(time, pft)` | `(1, n_pft)` | L4+shuffle |

The `1×max` chunk + deflate that already compresses the slab fill tail carries over verbatim (the
unused `[n+1:cap]` tail of every record is written as **`_FillValue`**, not zero — §4.3 — and deflates
near-free). (Deferred: a variance variable would get a `<name>_variance` companion, same dims, from the
`AGG_MEANSQ` buffer — no registry variable uses it yet, §3.5.) The per-record `global_id` slabs +
`n_cohort`/`n_patch` counts +
calendar are written for every cohort/patch-bearing stream. **Membership is always the today-CSR**
(`cohort_offset`/`cohort_count`/`owner_patch`, `AGG_LAST`) plus the `global_id` slab — the §4.4
fixed-within-window slot set means averaged and instantaneous records share one membership contract, so
each file is self-describing and independently readable by the existing reader. (The annual stream has no
cohort/patch dim at all, §3.1, so it carries no slab, no CSR, and no `global_id`.)

**CF metadata for aggregated data (goal G6) — a reader must be able to tell a mean from a snapshot.**
Every variable carries CF attributes emitted **from its `var_desc_t`**, closing the gap where averaged
records were previously indistinguishable from instantaneous ones:

| `agg` | `cell_methods` | on-disk units | notes |
|---|---|---|---|
| `AGG_MEAN`/`AGG_TMEAN` | `time: mean` | var units | dt-weighted for `TMEAN` |
| `AGG_FLUXSUM` | `time: sum` | rate units × time (drop `/s`) | period integral |
| `AGG_SUM` | `time: sum` | var units | count tally |
| `AGG_MIN` | `time: minimum` | var units | |
| `AGG_MAX` | `time: maximum` | var units | |
| `AGG_MEANSQ` | `time: mean` (+`_variance`: `time: variance`) | var units (var²) | mean+variance companion |
| `AGG_LAST` | `time: point` | var units | instantaneous snapshot |

Plus a `_FillValue` attribute (`NC_FILL_DOUBLE`/`NC_FILL_INT`) on **every** variable — the value the
serializer writes for an empty period (`valid=.false.`, §4.3) and for the unused slab fill tail — so a
reader never mistakes a fill for data.

**Time-coordinate semantics for a period (not a point).** An averaged record's `time` coordinate is the
period **midpoint**, and each file carries a **`time_bnds(2, time)`** bounds variable (`time:bounds =
"time_bnds"`) giving the period `[start, end)` in the same `time` units — the CF-standard way to stamp a
mean, so a reader integrates or plots against the correct interval rather than guessing. An `AGG_LAST`
stream stamps `time` at the instant and omits `time_bnds` (`cell_methods = time: point`). The calendar
companions (`year`/`month`/`day`) are the period-*start* stamp, matching the file name (§5.2).

**Provenance (paralleling the existing `pft_parameters.csv`).** Every diagnostic file carries global
attributes: `title`, `institution`/`source = "MEDS"`, the **run config echo** (`config_main`,
`config_pft` file names + a hash of their contents), the **git commit hash** + build compiler/flags
(threaded through from CMake as a generated `meds_build_info` constant), the `history` line, and CF
`Conventions = "CF-1.10"`. This makes a diagnostic file self-documenting to the same standard the run's
`pft_parameters.csv` already sets.

Only `NC_INT`/`NC_DOUBLE` and `NC_NETCDF4` are available in the C binding (§7); the registry's `xtype`
is constrained to those two (int for indices/ids/counts, double for everything else), which also fixes
the `_FillValue` fill constants to `NC_FILL_INT`/`NC_FILL_DOUBLE`. A lighter `NC_FLOAT` output is a
possible C-binding extension (§7.2), out of scope for P0.

### 5.4 State/restart stays separate

`io_write_state`/`io_read_state` (`<prefix>-S-<stamp>.nc`) are **untouched**: prognostic-only,
self-contained, no diagnostics, sized to `max(n,1)` (no caps), driven by `state_interval_years`. Keeping
the restart stream orthogonal to the diagnostic aggregation is deliberate — a checkpoint must be the raw
prognostic state at an instant, never a time-average. The manager owns both but never mixes them.

### 5.5 Performance: I/O off the compute path (goal G5)

The performance argument is layered, cheapest first; the recommended default needs only the first two:

1. **In-memory integration** — the daily kernels touch **no disk**; they only fold values into RAM
   buffers (a copy-reduce over the live `[1:n]` prefix, O(n) per step, ≪ the growth/mortality/fusion
   work). The current design already writes at most once a year; this keeps per-step disk touches at
   **zero**.
2. **Buffered, infrequent flush** — a flush writes one closed period's records (a day's / month's worth)
   in a burst, chunked+deflated, then returns. Flush frequency is the *output* cadence, not the *step*
   cadence: a daily stream flushes once per sim-day regardless of how many sub-steps ran. Between
   flushes the netCDF file handle stays open (append into the current time-chunk), so there is no
   per-record open/close cost. **This is the recommended P0/P1 default** — synchronous flush, because a
   flush is rare and small relative to the compute between flushes, so it never throttles.
3. **Double-buffering + asynchronous writer (P2, optional, build-gated)** — for the *fast* stream on a
   fast-loop run (where flushes are frequent and a synchronous HDF5 write could stall the step), the
   integrator can hold **two buffer sets**: the stepper fills buffer A while a **background writer**
   drains buffer B to disk; they swap at flush. **The async story has real platform constraints, stated
   honestly:**
   - **netCDF non-blocking I/O is NOT available on the HDF5 path.** MEDS's only format is
     `NC_NETCDF4`/HDF5, which has **no non-blocking put API** — the `ncmpi_i*` non-blocking puts are
     **PnetCDF / classic-format only**, and `nc_var_par_access`/`NC_INDEPENDENT` require an MPI-IO file
     opened with `nc_create_par`, which is neither in `meds_netcdf_c` nor meaningful single-process. So
     "netCDF's own async" is **scoped to a hypothetical future PnetCDF/MPI build**, not this one, and is
     dropped from the HDF5 design.
   - **The default ifx+netCDF build enables no host OpenMP runtime.** The `-mp` flag is `PUBLIC` only on
     the **nvfortran** demography target; on the flagship ifx build the `!$omp` lines are inert comments,
     so an `!$omp task` writer would run **inline/serially** and hide nothing. Async therefore requires
     an OpenMP-enabled host runtime — i.e. the **nvfortran `-mp` builds** — *or* a dedicated
     `pthread`/`c_f` std-thread spun via `iso_c_binding` (a portable alternative that does not depend on
     the compiler's OpenMP flag).
   - **HDF5 serializes all library calls under a global lock.** Even threadsafe HDF5 lets the background
     writer drain B while the stepper computes **only because compute does zero netCDF** — the design's
     "no per-step disk touch" is what makes async safe, and it **forecloses any future compute-side
     netCDF call** (noted so a later contributor does not break it).
   - Tradeoff: one extra buffer set (memory ≈ 2× the tiny integrator footprint), a lock/flag handshake,
     and a **run-end drain barrier**.

   **Recommendation: ship synchronous (options 1–2) at P0/P1, and treat async as an optional P2 that may
   well be dropped.** Because daily/monthly/annual flushes are already far off the critical path (§5.5
   intuition below), async only ever earns its keep for a *dense fast* stream on a fast-loop run — and
   even then the cheaper lever is to **coarsen the fast stream's `interval_steps`/file-chunk / variable
   set** rather than add threading. If P2 async is built, it lands on the nvfortran `-mp` build (or the
   pthread variant), never as "netCDF non-blocking on HDF5".

**Durability of long-held open handles (the cost of keeping the file open across a chunk).** Options 2–3
keep the HDF5 handle open across a whole time-chunk with no `nc_sync`, so a crash **loses the entire
in-flight chunk** — up to a sim-year for the annual stream, a sim-month for the daily. The design adds an
optional **`sync_every` knob** (`[output].sync_every = "flush" | "N records" | "never"`, default
`"flush"`): `nc_sync` after each flush (durable, one fsync per rare flush — negligible) or every N
records or never (fastest, crash-forfeits the chunk). A new `nc_sync` C binding is the only addition
needed (§7.2). On restart the diagnostic streams start **fresh files** anyway (open Q4), so a crash at
worst forfeits the current open chunk, never corrupts earlier ones.

Quantified intuition: the integrator footprint is `Σ_tiers (nvars × cap)` doubles — order
(4 tiers × ~25 vars × ~4096 cohorts × 8 B) ≈ a few MB, trivially RAM-resident. A daily-stream flush
writes ~25 vars × (1 record × cap, deflated) ≈ tens of KB on disk — microseconds of write beside a
sim-day of daily kernels. The synchronous default is comfortably off the critical path for everything
except a hypothetical dense fast stream, which is exactly what P2's async covers.

---

## 6. TOML config surface

### 6.1 Two-file design (main config + optional `meds_io_config.toml`)

Output control is split across **two files**, exactly mirroring how the main config already names a
separate **PFT config** file (`[init].pft_config`):

- **The MAIN config `[output]` block** carries the **high-level, low-cardinality** knobs a user reaches
  for most: the master switch, dir/prefix, caps, per-frequency enable + file-chunk, and — new —
  **flux-group toggles** (`structure` / `carbon_fluxes` / `water_fluxes` / `energy_fluxes` — `structure`
  is deliberately unsuffixed, the three flux groups carry `_fluxes`) that switch whole `GRP_*` groups
  of variables on or off. **This block alone is a complete, usable configuration** — with no
  `meds_io_config.toml` at all, the registry defaults (§3.5) plus these group/frequency switches fully
  determine what is written. This is the common case ("give me the daily structure + carbon fluxes,
  skip water and energy").
- **An OPTIONAL separate `meds_io_config.toml`**, named from the main config
  (`[output].io_config = "meds_io_config.toml"`), carries the **fine-grained per-variable** table — one
  entry per variable to force it on/off or re-target its streams. It exists only when a user wants
  variable-level control; if the key is absent, output is fully governed by the main block. This keeps
  the (potentially long) per-variable list out of the main config, just as the PFT trait table lives in
  its own file.

**Resolution order in `build_output_registry`** (later steps win):
1. Registry **source defaults** — each var's `group` + default `streams` (§3.5).
2. **Group toggles** (main config) — `apply_group_toggles`: if `structure` / `carbon_fluxes` /
   `water_fluxes` / `energy_fluxes` `= false`, clear every stream bit of every var in that `GRP_*` group.
3. **Per-frequency enable** (main config) — a tier with `enabled = false` has its `FREQ_*` bit cleared
   from every var's effective mask (that whole stream is off).
4. **Per-variable overrides** (`meds_io_config.toml`, if named) — `apply_io_config`: a `name = bool`
   forces the var on (at its registry-default streams) or off; a `name = "F D M"` **replaces** its
   stream mask. These are the finest granularity and **override** steps 2–3.
5. Enforce the **cohort/patch ⇒ no-annual** guard (§3.1) and build `idx_freq`.

**Reader mechanics.** Both files use the existing flat `section.key=value` parser (limitation L8 is
sidestepped, not lifted): the group toggles and per-frequency knobs are plain booleans/strings under
`[output]`; the per-variable table is a `[variables]` **sub-table of `name = bool` / `name = "F D M"`
scalars** in `meds_io_config.toml`, stored as fully-qualified `variables.<name>` keys — no string-array
or arrays-of-tables machinery needed. `meds_io_config.toml` is loaded with the same `meds_toml` reader
as the PFT file.

**Value grammar (specified, not sniffed loosely).** A `variables.<name>` value is one of: (a) a **bool**
— `true` force-enable at the registry default streams, `false` disable everywhere; or (b) a
**stream-string** — a space-separated subset of the tokens `{F, D, M, Y}` (e.g. `"M Y"`), parsed by
`parse_stream_mask(s) → ior(FREQ_*)` (unrecognized token, or `Y` on a cohort/patch var → `error stop`
naming the offender). The loader disambiguates by TOML value *type* (the flat parser retains each stored
value's scalar type). `file_chunk` is likewise a **named-bucket string** mapped by
`parse_file_chunk(s) → {FC_DAY, FC_MONTH, FC_YEAR, FC_RUN}` filling `file_chunk(N_FREQ)`.

**Unknown/typo'd key validation (closing the silent-ignore trap).** After building the registry the
loader **iterates the stored `variables.*` keys and `error stop`s (strict) or warns (lenient,
`[options]`-gated) on any key that does not match a registry variable name** — listing the offenders,
exactly as `load_meds_config` already lists missing mandatory keys. The same enumerate-and-check guards
unknown keys under `[output]` (a mistyped `file_chnk` or `carbn`). Versatility without a silent-failure
UX trap.

Backward-compat: `load_meds_config` keeps reading the legacy `[io]` keys if `[output]` is absent
(mapping `output_interval_years` → an annual stream with all site-level vars `AGG_LAST`), so existing
TOML runs unchanged. New keys use **optional-with-default** loaders (a new `opt_l`/`opt_s`/`opt_i`
alongside the existing `req_*`), because the group/frequency knobs must default sensibly rather than
force the user to name every variable — a deliberate softening of the "all `[io]` keys mandatory" rule
for the output block only.

### 6.2 Example MAIN-config `[output]` block (complete on its own)

```toml
[output]
enabled       = true            # master switch (replaces io.write_output)
dir           = "out"           # replaces io.output_dir
prefix        = "run01"         # replaces io.output_prefix
cohort_max    = 4096            # cohort/patch dim caps (replaces io.cohort_max / io.patch_max)
patch_max     = 256
strict_caps   = false           # false: warn+truncate on n>cap; true: error stop (regression mode)
sync_every    = "flush"         # nc_sync durability: "flush" | "<N> records" | "never"  (§5.5)
io_config     = "meds_io_config.toml"   # OPTIONAL per-variable override file. If PRESENT, its per-variable
                                        # entries OVERWRITE the high-level group/frequency options below
                                        # (finest granularity wins, §6.1). Omit -> group/freq settings alone rule.

  # --- high-level flux-group toggles: switch whole GRP_* groups without naming variables. ---
  # --- These are the coarse control; any variable named in meds_io_config.toml overrides them.  ---
  # --- Comments list the GROUP's membership (some entries are forward-looking; only vars in the ---
  # --- §3.5 registry are written today: structure*, gpp_site/nee_site, transp_site/soil_water_site, soil_temp_site). ---
  structure     = true          # demography/size group: nplant, dbh, agb, patch geometry, ids
  carbon_fluxes = true          # carbon group: gpp, nee (+ future npp/Rh)   (needs the fast loop to populate)
  water_fluxes  = false         # water group: transpiration, soil moisture (+ future evap)
  energy_fluxes = false         # energy group: soil temperature (+ future leaf temp/sensible/latent)

  # --- per-frequency streams: enable + records-per-file (the time-chunk) ---
  [output.fast]
  enabled       = false         # needs the fast loop; off by default
  interval_steps = 4            # flush a fast record every 4 * dt_fast (e.g. 4*15min = hourly); 1 = every dt_fast
  file_chunk    = "day"         # one file per sim-day

  [output.daily]
  enabled    = true
  file_chunk = "month"          # one file per sim-month

  [output.monthly]
  enabled    = true
  file_chunk = "year"           # one file per sim-year

  [output.annual]
  enabled    = true
  file_chunk = "run"            # ONE file for the whole run (site-level only, tiny; §5.1)

[state]                          # the restart stream stays its own block (unchanged semantics)
write_state          = true
state_interval_years = 10
```

### 6.3 Example `meds_io_config.toml` (optional, fine-grained)

```toml
# Per-variable output overrides. Named from the main config via [output].io_config.
# Each entry WINS over the main config's group/frequency defaults (§6.1 resolution order).
#   value = bool  -> force ON (at registry-default streams) / OFF everywhere
#   value = "F D M" -> replace the stream mask (tokens F/D/M/Y; Y forbidden on cohort/patch vars)
[variables]
growth_avg_cohort = false   # never write this variable, even though 'structure' is on
agb_cohort        = "M"     # narrow from the registry default (D M) to monthly only
soil_temp_site    = "F D"   # fast + daily, even though 'energy_fluxes' is off in the main config
gpp_site          = "F D M" # put the site GPP flux on the fast + daily + monthly tiers (default is M only)
```

### 6.4 Config type additions

`meds_config_t` gains an `[output]` sub-record `output_config_t`, a plain-data type. **The selector codes
it needs — `N_FREQ`, `FC_*`, `GRP_*`, `FREQ_*` — live in `src/shared` (a new `meds_output_config`
module, beside the type), NOT in `meds_output_types` (the io library).** This mirrors the forcing
precedent exactly (CLAUDE.md: `forcing_config_t` + its selector codes live in `src/shared/meds_forcing_config`
so `meds_config` carries the block **without a `shared → io` back-edge**); `meds_output_types` then
`use`s `meds_output_config` for those codes (and adds the io-only `AGG_*`/`DIM_*` + the `var_desc_t`/
`integ_buffer_t`/… types). So the config type at the DAG root pulls no io dependency:

```fortran
! codes below are declared in src/shared/meds_output_config (with output_config_t), NOT in the io layer.
type :: output_config_t
   logical            :: enabled     = .false.
   character(len=256) :: dir = '.', prefix = 'meds'
   character(len=256) :: io_config   = ''               ! path to meds_io_config.toml; '' -> group/freq defaults only
   integer(ik)        :: cohort_max = 4096_ik, patch_max = 256_ik
   logical            :: strict_caps = .false.
   integer(ik)        :: sync_every  = FC_FLUSH          ! nc_sync policy: FC_FLUSH | N | FC_NEVER (§5.5 durability)
   logical            :: grp_on(4)   = [.true.,.true.,.false.,.false.] ! GRP_STRUCTURE/CARBON/WATER/ENERGY toggles
   integer(ik)        :: fast_interval_steps = 4_ik      ! fast tier: flush every N*dt_fast (§4.1)
   logical            :: freq_on(N_FREQ)   = [.false.,.true.,.true.,.true.]  ! F D M Y defaults
   integer(ik)        :: file_chunk(N_FREQ) = [FC_DAY, FC_MONTH, FC_YEAR, FC_RUN]  ! records-per-file bucket per tier
   ! per-variable overrides are NOT stored here -- build_output_registry loads meds_io_config.toml (if
   ! io_config /= '') and queries variables.<name> by registry name, so this type stays fixed-size and the
   ! variable list stays source-defined (the registry is the single source of truth for WHICH vars exist).
end type
```

The group toggles + per-frequency knobs live here (fixed-size, in the main config); the per-variable map
is applied inside `build_output_registry` by loading `meds_io_config.toml` and looping the registry names
— so the config type stays fixed-size and the variable list stays source-defined (TOML only toggles what
the registry already declares).

---

## 7. netCDF layer

### 7.1 Reuse of `meds_netcdf_c`

The C bindings (`nc_create/def_dim/def_var/put_att_text/def_var_deflate/def_var_chunking/enddef/
put_vara_double/put_vara_int/put_var1_double/close/open/inq_varid/get_vara_int/get_vara_double/
strerror`, plus `cstr`/`nc_check` and the `NC_*` params) are **sufficient as-is** for the whole design:
the stream writer defines dims/vars from the registry, chunks+deflates each slab var, and appends
records exactly like `io_write_snapshot` does today — just in a loop over registry entries instead of a
hard-coded list, and to a rolling file instead of one. **No new C binding is required for P0/P1.**

### 7.2 Possible (optional) C-binding extensions

- `nc_sync` (`nf90_sync` / `nc_sync`) — **the one binding P1 wants**: the periodic-flush durability knob
  (§5.5), a one-line C wrapper, so a crash forfeits at most the current open chunk.
- `NC_FLOAT` (`nc_put_vara_float`) for half-size output if disk becomes the constraint — a
  single-precision diagnostic stream is defensible (diagnostics are not restarted from). Out of scope
  for P0.
- **No netCDF non-blocking/parallel binding is proposed for the HDF5 path.** `ncmpi_i*` non-blocking
  puts and `nc_var_par_access`/`nc_create_par` are PnetCDF/MPI-IO features absent from `NC_NETCDF4`
  single-process (§5.5); the optional P2 async writer uses a **host thread** (OpenMP `-mp` build or an
  `iso_c_binding` pthread), **not** a new netCDF binding. Should MEDS ever gain an MPI/PnetCDF build,
  the parallel bindings land with that work, not here.

### 7.3 Structuring the data-driven stream writer

`meds_output_stream` holds a `stream_file_t` (one per active frequency) carrying the open `ncid`, the
current time-chunk bucket, `nrec`, and — crucially — a **registry-indexed `var_id` array**
(`integer(c_int), allocatable :: vid(:)`) instead of today's fixed named fields (`v_c_pft`, `v_c_agb`,
…). This is what makes it data-driven: `stream_open_file` loops `reg%idx_freq(:,freq)` and defines +
stores each var-id by registry index, so `stream_write_record` writes by the same index. Adding a
variable never touches this module (it reads the registry). netCDF is a hard dependency, so this module
is always compiled — there is no stub twin to keep in sync, and the DAG-hygiene split (§2) keeps this
serializer object off the `aux` integrate edge without any build switch.

---

## 8. Test plan

Unit tests (CTest targets; the integrator/registry tests touch no netCDF but still link `libmeds_io.a`,
which is always built; build the **nvfortran multicore** back end too — issue #7):

1. **Integrator arithmetic** (`test_output_integrate`): feed a known scalar sequence; assert
   `AGG_MEAN` = arithmetic mean, `AGG_TMEAN` on **non-uniform `dt`** = the dt-weighted mean (≠ plain
   mean), `AGG_FLUXSUM` = `Σ(x·dt)` (units check: a constant rate over T seconds → `rate·T`), `SUM` =
   sum, `MIN`/`MAX` = extrema, `LAST` = final; assert `reset_buffer` re-seeds (min→+huge, max→−huge,
   mean/wsum→0). **Zero-sample guard:** normalize a never-touched buffer → `valid=.false.`, output =
   `_FillValue`, **no NaN and no ±huge leak**. (The `AGG_MEANSQ` variance path is a P1 test — variance is
   deferred, §3.5.)
2. **Chaining exactness & feeder routing**: (a) drive 30 daily values → monthly `AGG_TMEAN`
   equals the direct dt-weighted mean of the 30, annual = mean of 12 monthly means; `AGG_FLUXSUM` chain =
   sum of sub-period integrals. (b) **Feeder configs:** run
   **monthly-only**, **annual-only** (site-level vars), and **daily+annual-skip-monthly** and assert each
   variable's finest active tier integrates raw state directly and coarser tiers chain the right feeder
   — the monthly-only buffer is fed every step (not left empty), the daily+annual case chains annual off
   daily. Also: a **fast tier** feeding daily (`interval_steps` sub-steps per fast record, fast mean
   chains into daily).
3. **Fixed-within-window slab integrate** (`test_output_integrate`): drive a multi-step window over a
   **constant** cohort slot set; assert each slot's `AGG_TMEAN` is the dt-weighted mean over the window
   and the flushed slab is in slot order with the today-CSR intact. Assert the **boundary alignment**:
   the monthly cohort buffer is normalized+reset *before* a simulated fiss/fuse changes the slot set, so
   no window straddles the change. Assert the **cohort/patch ⇒ no-annual guard** `error stop`s when a
   `DIM_COHORT` var is given `FREQ_ANNUAL`. Assert `n > cap` hits the same strict/warn guard as today.
   Assert **`DIM_SOIL`/`DIM_PFT`** vars fold by direct index (same code path).
4. **Registry + group/variable toggles** (`test_output_registry`): build the registry; (a) flip a
   main-config **group toggle** (`water_fluxes=false`) and assert every `GRP_WATER` var is absent from
   `idx_freq`; (b) load a `meds_io_config.toml` with one `false` and one stream-string override and
   assert the override wins over the group/frequency defaults (§6.1 resolution order). **Unknown-key
   trap:** a typo'd `[variables]` name → `error stop` (strict) / warning (lenient) naming the offender,
   not a silent no-op (§6.1).
5. **Round-trip / back-compat** (`test_output_roundtrip`, needs netCDF): run a short spin-up with the
   **per-cohort P0 fields at `AGG_LAST` on a daily/monthly stream**; assert the produced `-D-*`/`-M-*`
   per-cohort slab is **equal to today's `-D-output.nc` up to metadata** — the within-window slot set is
   fixed (§4.4), so rows are in the **same slot order** as today (no permutation), and only the added
   `cell_methods`/`_FillValue` attributes differ. Then assert record counts match the calendar, rollover
   at month boundaries, `time`/`year`/`month`/`day` exact across a leap year, the annual `-Y.nc` is a
   single whole-run file of **site-level scalars only** (no cohort dim), and that each averaged variable
   carries the right **`cell_methods`**, a **`time_bnds`** interval, and a **`_FillValue`** attribute (§5.3).
6. **File organization**: assert one file per stream per time-chunk (glob the output dir), the annual
   stream is a single stampless `-Y.nc`, the state stream `-S-*.nc` is untouched, and that
   `strerror`-checked writes all return `NC_NOERR`.
7. **Performance smoke**: a multi-decade daily run; assert wall-time with output ON is within a few
   percent of output OFF (I/O off the critical path), and peak RSS grows by only the integrator
   footprint.

---

## 9. Phasing

**P0 — registry + daily/annual means + group/variable toggles + split files (the core value).**
- `meds_output_types` (`var_desc_t` incl. `group`, `integ_buffer_t` incl. `scal2`/`slab2`/`wsum`/
  `wsum_slab`, `output_registry_t`, netCDF-free `output_manager_t`, `AGG_*`/`FREQ_*`/`DIM_*`/`GRP_*`/
  `FC_*` codes).
- `build_output_registry` with the P0 table (§3.5), `extract_variable`/`extract_site`, the per-var
  **feeder map** (§4.1), the group toggles + resolution order (§6.1), and the cohort/patch⇒no-annual guard.
- Daily + annual tiers; `AGG_MEAN/TMEAN/SUM/LAST` operators (dt carried, empty-period `_FillValue`
  guard); scalar + **fixed-slot-set slab** integrate (§4.4, no id map); `DIM_SOIL`/`DIM_PFT` share the
  slab path.
- `meds_output_stream`: per-stream, per-time-chunk files, data-driven var-id array (always compiled — no
  stub), **`cell_methods` + `time_bnds` + `_FillValue`** CF metadata, and **provenance** globals (config
  echo + git hash, §5.3).
- Main-config `[output]` block (group toggles + per-frequency knobs) + optional `meds_io_config.toml`
  per-variable map (value grammar + **unknown-key validation**); `opt_*` loaders; legacy `[io]` back-compat.
- `output_integrate` (netCDF-free) wired into `meds_stepper` (monthly cohort flush **before** fiss/fuse,
  §4.5); `output_serialize_pending` driven from `meds_main` (the flush wall) — replaces the `is_new_year`
  write block.
- Tests 1, 3, 4, 5, 6. **Deliverable:** user-mutable daily+annual output in split files, time-averaged,
  CF-tagged, the current per-cohort stream reproducible **up to metadata** (§3.6).

**P1 — fast + monthly tiers, min/max, flux integrals, durability.**
- Monthly tier + full chaining; `AGG_MIN/MAX`; `AGG_FLUXSUM` flux integrals (fed by the fast loop).
  (`AGG_MEANSQ`/variance stays defined but is **deferred** — no registry variable uses it yet, §3.5.)
- **Fast** tier (`interval_steps × dt_fast`, §4.1) integrating inside the fast biophysics sub-step loop
  and consuming its CO₂/flux/column diagnostics (`gpp_site`, `nee_site`, `soil_temp_site`, … — put on the
  fast tier via `meds_io_config.toml`, default `M`); the P1 registry rows (§3.5).
- `nc_sync` durability knob (`sync_every`, §5.5) + its one C binding.
- Tests 2 (fast-feeds-daily), and the flux paths of 1.

**P2 — optional asynchronous writer (may be dropped).**
- Only if a dense fast stream proves to stall the step: double buffer set + **host-thread** background
  flush (nvfortran `-mp` OpenMP task **or** an `iso_c_binding` pthread — **not** netCDF non-blocking,
  which is unavailable on the HDF5 path), run-end drain barrier. Gated on an OpenMP-enabled build.
- Optional `NC_FLOAT` C binding for lighter diagnostic output.
- A perf test (async variant: assert no compute stall on the fast stream).

Each phase is independently shippable and leaves the model runnable; P0 alone removes limitations
L1–L6, L8 and gives the user everything in the stated goals except sub-daily/flux/variance (P1) and the
last slice of write-latency hiding (P2, optional).

---

## 10. Open questions

1. **Area-weighted spatial roll-up for scalar means.** The site totals are already area-weighted patch
   reductions (`meds_demography_diagnostics`); should any *cohort-dimensioned* variable also be offered
   as an nplant- or area-weighted site aggregate in the same registry (a `DIM_SCALAR` derived twin), or
   is that left to post-processing? Leaning: expose the existing site totals as the canonical scalars;
   defer arbitrary weighted roll-ups to the reader.
2. **Weighting inside `AGG_MEAN` for cohorts that change `nplant` within a window.** A simple per-step
   mean weights every step equally; an nplant- or agb-weighted temporal mean may be more ecologically
   meaningful for per-plant quantities. Proposed: default equal-weight (matches ED2's `dmean`); add an
   optional `AGG_WMEAN` (weighted by a companion source, e.g. nplant) in P1 if wanted.
3. **`file_chunk = "run"` for high-cadence streams.** A whole-run single file for the annual stream is
   the default (few, tiny records); for daily it defeats FILE-SIZE. Should the loader **reject** a
   `file_chunk` coarser than N records for sub-monthly streams, or just warn? Leaning: warn, trust the user.
4. **Restart interaction with open diagnostic files.** On restart mid-run (`INIT_RESTART`), the
   diagnostic streams start **fresh files** at the restart date (no append into a pre-restart file). Is
   that acceptable, or should the manager detect and append to the matching time-chunk file? Leaning:
   fresh files (simpler, and the calendar stamp disambiguates); document it. (The whole-run annual `-Y.nc`
   is the one stream that would want appending — note it.)
5. **Multi-site / polygon dimension.** The whole layout is single-site (no polygon axis), matching the
   engine. When multi-site lands, does a `site`/`polygon` dim slot in as a *leading* fixed dim on every
   var (cheap, uniform) — and does that interact with the per-time-chunk file split (one file per
   stream per chunk **per site**, or a site dim inside one file)? Deferred, but the registry/`dim`
   enum should reserve a `DIM_*` leading-axis convention now so it is not a breaking change later.
6. **Keep P2 async at all?** With daily/monthly/annual flushes already microseconds beside a sim-day of
   kernels, async only ever pays on a *dense fast* stream — and even there, coarsening `file_chunk` /
   `interval_steps` or the fast variable set is the cheaper lever. Leaning: **do not build async until a
   measured fast-stream stall demands it**, and if built, on the nvfortran `-mp` build (or a pthread),
   never as netCDF non-blocking. Revisit only after an hourly-scale fast stream is exercised in anger.
7. **Time stamp for a period mean.** Chosen (§5.3): `time` = period **midpoint** with a `time_bnds`
   `[start, end)` interval (CF-standard). Alternative conventions (stamp at start, or at end as ED2's
   "previous period" file naming implies) are equally decodable; the midpoint+bounds pair is picked so a
   plotter without bounds-awareness still places the mean sensibly. Confirm against the post-processing
   readers (`post_proc/`) before P0 freeze.
8. **`AGG_*` classification per variable.** The registry must tag each rate as either a **mean rate**
   (`TMEAN`, keeps `/s`) or a **period integral** (`FLUXSUM`, drops `/s`) — never plain `SUM`, which is
   reserved for integer count tallies (§3.3, §3.5). The P1 default is `TMEAN` (a mean rate); a `FLUXSUM`
   twin is offered if a user wants a period total. Is a hand audit of the P0 vars enough, or should the
   registry assert units consistency (a `/s` in `units` ⟹ `FLUXSUM`/`TMEAN`, never `SUM`)? Leaning: a
   compile-time-style audit plus a units-sanity assertion in the registry test.
