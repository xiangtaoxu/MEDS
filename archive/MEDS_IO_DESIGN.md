# MEDS Diagnostic Aggregation & Output — Design

A **stateful, streamed, multi-frequency** diagnostic-output subsystem for MEDS (`src/io/`) that turns
the instantaneous once-a-year snapshot the model writes today into a general **aggregation-and-output
engine**: every diagnostic variable declares itself once in a **registry** (name, long name, units,
dimensionality, aggregation operator, the output streams it belongs to, and an on/off flag), a set of
per-frequency **accumulators** integrate it in memory as the stepper ticks (hourly → daily → monthly →
annual, `fmean → dmean → mmean → qmean`-style chaining), and a **flush** writes each closed period to a
**time-chunked, per-stream netCDF file** through the existing `meds_netcdf_c` C-binding layer. All of it
is **user-mutable from TOML** — enable/disable a whole frequency, toggle any single variable, all with
no source edit.

The reference architecture is **ED2's `var_table` + `average_utils` + `h5_output`** (registration ⊥
temporal averaging ⊥ serialization; a variable opts into *N* streams via *N* keyword flags, one
`integrate/normalize/zero` family per temporal tier). MEDS keeps that **decoupling** wholesale but
**rejects the mechanism**: ED2's `vt_info(maxvars,ngrids)` is a *global mutable registry of raw
pointers into the reallocatable SoA*, requiring a fragile `filltab_alltypes` re-hash after every
birth/death/fusion — exactly the hidden-global-mutable-state pattern MEDS/CLAUDE.md forbids. MEDS
instead carries an **explicit, passed** registry (`output_registry_t`), and the accumulators **copy the
reduced value out of the live SoA each step** rather than caching pointers into it — so cohort/patch
count changes are a non-event (no re-hash, ever).

The design honours the three user goals directly: **PERFORMANCE** — in-memory accumulation + infrequent
buffered flush keeps I/O off the daily critical path (the load-bearing fact: daily/monthly/annual
flushes are rare and tiny); a P2 background-thread writer is an **optional, build-gated** extra for the
one dense case (an hourly stream), *not* a load-bearing part of the perf story (§5.5); **FILE SIZE** —
one file **per stream per time-chunk** (a sim-year or sim-month), never one lump, all chunked +
deflated; **VERSATILITY** — user-mutable hourly/daily/monthly/annual streams with a per-variable on/off
toggle, from TOML.

Adding a variable is **two coupled edits, not one** — an honest accounting of the price of MEDS's
no-pointer registry: one `add()` line in the registration list *and* one `case` in the `harvest_var`
switchboard that copies the value out of the SoA (§3.3). This is still far better than today's three
scattered edit sites (L5), and it is pointer-free.

Kinds `wp/ik`, `_wp` literals, `implicit none`, `error stop`, ≤132 columns, lowercase free-form,
nvfortran-safe (no array-valued function result passed straight into a call — issue #7). This is a
**design/plan document only**; no source is modified here.

---

## 1. Scope, current state & limitations

### 1.1 What MEDS writes today

There is exactly **one** diagnostic path (`src/io/meds_io.f90`, `libmeds_io.a`, on by default) over a
thin `iso_c_binding` wrapper of the netCDF **C** library (`meds_netcdf_c.f90`), with a no-op twin
(`meds_io_stub.f90`, same module name) linked when `-DMEDS_ENABLE_IO=OFF`. It has two write streams and
a restart reader, all driven from `meds_main.f90`:

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
| L1 | **Annual cadence only.** Writes fire only at `is_new_year` in `meds_main` (`mod(iyear,interval)==0`). | Finest frequency = 1 year, though the engine steps daily/monthly. No monthly/daily/hourly output. |
| L2 | **Instantaneous only.** No accumulation/averaging anywhere; the module header states "no buffering beyond the live state is needed." | No time-mean, sum, min/max, or variance diagnostic is possible. |
| L3 | **One frequency for the whole stream.** A single `io_output_interval_years` governs every variable. | Cannot write different variables at different cadences. |
| L4 | **No per-variable toggle.** Every var is defined + written unconditionally; only knob is `io_write_output` on/off. | All-or-nothing; no config-selectable subset. |
| L5 | **Hard-coded variable list.** Adding one var = a 3-site Fortran edit (a var-id field in `meds_io_t`, a `def_*` in `io_create`, a `put_vara_*` in `io_write_snapshot`). | No data-driven registry. |
| L6 | **One lump file** grows for the whole run (`-D-output.nc`). | A century run is one ever-growing file; no per-year/-month split. |
| L7 | **Fixed cohort/patch dims** at config caps; live `n > cap` is a hard `error stop`. | Caps must be sized a priori; ragged growth across a window unhandled. |
| L8 | **TOML can't express per-variable config.** Flat `key=value`, real arrays only — no string arrays, no arrays-of-tables. | A variable-toggle list / per-variable frequency table is not representable. |
| L9 | **Stub drift risk.** `meds_io_stub.f90` must mirror the real public API by hand. | Every new public routine must land in both. |

The single choke points are `io_create` (schema) + `io_write_snapshot` (record write) in `meds_io.f90`,
the `is_new_year`-gated block in `meds_main.f90` (~lines 126–137, the cadence owner), and the `[io]`
keys in `meds_config.f90` + `meds_config_io.f90`. This design replaces all three.

### 1.3 Design goals & non-goals

**Goals.** (G1) A **DRY variable registry**: one registration list feeds all streams; adding a variable
is a registry `add()` line plus a `harvest_var` case (two adjacent, co-located edits — §3.3 — replacing
today's three scattered sites). (G2) **Multi-frequency in-memory aggregation** (hourly/daily/monthly/
annual) with mean/sum/min/max/last/mean-of-squares operators **plus a dt-weighted flux integrator** for
rate quantities (§3.5, §4.3). (G3) **Per-frequency + per-variable user toggles** from TOML, with
unknown-key validation. (G4) **Time-chunked per-stream files**, chunked+deflated, never one lump. (G5)
**I/O off the compute critical path** (in-memory accumulate + buffered flush; optional P2 async on
OpenMP-enabled builds only). (G6) Stub stays interface-identical when IO is OFF. (G7) **CF-compliant
averaged output**: `cell_methods` + `time_bnds` so a reader can tell a mean from an instantaneous value.

**Non-goals (this design).** Multi-site/polygon output (single-site layout retained; §12 notes the
seam). CO₂/flux tower sub-`dt_fast` output (the fast loop is not yet wired into `meds_main`; the hourly
tier is designed to consume fast-loop diagnostics *when* they land, §3.5). The met-forcing reader. The
diel-cycle (`qmean`) tier is **P1**, not P0. On-disk format stays netCDF-4/HDF5 via the existing C
bindings (no zarr).

---

## 2. Where it lives (library DAG)

The DAG is `shared ← {allometry, plant} ← state ← demography ← aux ← main`. The aggregation engine
splits cleanly along the **state/process wall** already used everywhere else, and — crucially — along a
**netCDF wall** drawn *inside* the output family so the DAG and stub boundary stay honest:

- **The registry + accumulators are pure data + pure kernels, and are netCDF-FREE** — they name only
  `site_t` (read-only) and plain arrays, exactly like `meds_demography_diagnostics` (which they call).
  `aux` (the stepper) references *only* these, for the per-step accumulate tick. Because a Fortran static
  archive pulls only referenced objects, `aux` linking `libmeds_io.a` but calling only the accumulator
  symbols **does not pull the serializer object, and so does not pull netCDF** — the one new edge
  (`aux → io-accum`) carries no C dependency.
- **The netCDF serialization is quarantined in `meds_output_stream`** (over `meds_netcdf_c`), gated by
  `MEDS_ENABLE_IO` with the stub twin — unchanged in *mechanism*, generalized in *content*. It is
  referenced **only from `meds_main`**, which drives every disk flush. The stepper never calls it.

**The stepper accumulates; main serializes (the resolution of the §2↔§4.5 seam).** `advance_one_step`
calls the netCDF-free `output_accumulate` each step (fold live state into buffers; at a period roll-over,
normalize + chain + reset + **stage** the closed record into a plain-array pending queue on the manager —
all buffer arithmetic, zero disk). After `advance_one_step` returns, `meds_main` calls
`output_serialize_pending(mgr)` — the **only** netCDF-touching entry, stubbed to a no-op when IO is OFF.
So the flush lives at `main` (honouring "`io` is a leaf `main` links"), the stepper stays netCDF-free,
and `output_manager_t` is a **netCDF-free plain-data type** (it holds `stream_file_t` handles whose
`ncid`s are just integers; only the *operations* on them touch C, and those live in the serializer).

```
shared ─┬─ allometry ─ state ─ demography ─┬─ aux (stepper: ticks the netCDF-FREE accumulators each step;
        │                                  │       normalize+chain+reset+STAGE at period roll-over; NO disk)
        │                                  │         │ links only meds_output_{types,registry,accum}
        ├─ plant                           ├─ io  (meds_output_types     ← var_desc_t, accum_buffer_t, output_manager_t [netCDF-free])
        │                                  │      (meds_output_registry  ← the registration list + TOML overrides + freq index)
        └─ biophysics                      │      (meds_output_accum     ← harvest_var + accumulate/normalize/reset kernels)
                                           │      (meds_output_stream    ← per-stream netCDF file lifecycle  ← the netCDF wall)
                                           │      (meds_output_stream_stub← no-op twin, same module name)
                                           │      (meds_netcdf_c         ← REUSED verbatim: the C bindings)
                                           │      (meds_io               ← state/restart, UNCHANGED)
                                           └─ main (owns output_manager_t; drives output_serialize_pending → the ONLY flush)
```

**Files & CMake.** `CMakeLists.txt` globs `src/io/*.f90` into `libmeds_io.a`, so new modules auto-add.
`MEDS_ENABLE_IO=OFF` must select a **stub twin for the whole new family** (§7.4): the serializer modules
(`meds_output_stream`) get a no-op twin; the registry + accumulator modules are **netCDF-free and always
compiled** (they only touch `site_t` and plain arrays — no C bindings), so the stub boundary is drawn at
the serializer, not the accumulator. This keeps the in-memory aggregation testable with `-DMEDS_ENABLE_IO=OFF`.

| File | Role | Analogue (ED2 / MEDS) |
|---|---|---|
| `src/io/meds_output_types.f90` (new) | `var_desc_t`, `stream_desc_t`, `accum_buffer_t`, `output_registry_t`, `output_manager_t`; the `AGG_*` / `FREQ_*` / `DIM_*` selector codes | ED2 `var_table` + `idim_type`; MEDS `meds_biophysics_types` `SOIL_*`/`ENERGY_*` block |
| `src/io/meds_output_registry.f90` (new) | `build_output_registry(cfg) → output_registry_t`: the ONE registration list; applies TOML enable flags | ED2 the `vtable_edio_*` call sites; MEDS `build_soil_params` |
| `src/io/meds_output_accum.f90` (new) | `pure`/`elemental` accumulate/normalize/reset kernels over `accum_buffer_t`; `harvest_site(site, registry) → raw per-var slabs` | ED2 `integrate_/normalize_/zero_ed_*mean_vars` (~8700 lines → one generic per-var loop) |
| `src/io/meds_output_stream.f90` (new) | per-stream netCDF file lifecycle: `stream_open_file`, `stream_write_record`, `stream_close_file`; time-chunk rollover | ED2 `h5_output`; MEDS `meds_io%io_create/io_write_snapshot/io_close` |
| `src/io/meds_output_stream_stub.f90` (new) | no-op twin of `meds_output_stream` (same module name) | `meds_io_stub.f90` |
| `src/io/meds_output_manager.f90` (new) | `output_manager_t` glue (netCDF-free plain data: registry + accum buffers + `stream_file_t` handles + the pending-flush stage). `output_accumulate(mgr, site, cadence)` — the per-step, netCDF-free tick (called by the stepper). `output_serialize_pending(mgr)` — drain the stage to disk (called by `main`; the serializer wall) | ED2 `ed_output`; new |
| `src/io/meds_io.f90` (unchanged) | state/restart streams | — |
| `test/test_output_accum.f90`, `test/test_output_registry.f90`, `test/test_output_roundtrip.f90` (new) | CTest | `test_column_hydrology` etc. |

**Build nvfortran multicore on every new module** — a green ifx suite is not sufficient (CLAUDE.md
issue #7). The accumulate kernels are plain-array arithmetic (a natural fit for the same
`target teams distribute` discipline as `growth_step`, though P0 keeps them host-only — §4.5).

---

## 3. The variable registry (`var_desc_t`) — a modern `var_table`

### 3.1 The descriptor

Every diagnostic variable is one **`var_desc_t`** — a pure DATA descriptor, no pointers. The registry is
a fixed, source-defined **list** of these (the "registration list"), built once at start-up and then
**immutable**; the TOML surface only flips the `enabled` flag and the per-frequency membership mask.

```fortran
! selector codes (in meds_output_types, beside the type defs) ------------------------------!
integer(ik), parameter :: AGG_MEAN = 1_ik, AGG_SUM = 2_ik, AGG_MIN = 3_ik, AGG_MAX = 4_ik,   &
                          AGG_LAST = 5_ik, AGG_MEANSQ = 6_ik,    & ! MEANSQ pairs with MEAN -> variance
                          AGG_TMEAN = 7_ik, AGG_FLUXSUM = 8_ik     ! dt-weighted: state time-mean / flux integral

integer(ik), parameter :: DIM_SCALAR = 0_ik, DIM_COHORT = 1_ik, DIM_PATCH = 2_ik,            &
                          DIM_SOIL = 3_ik,  DIM_PFT = 4_ik        ! trailing axis
! ragged (id-keyed, born/die/fuse mid-window): DIM_COHORT, DIM_PATCH.  fixed-index (stable axis): DIM_SOIL, DIM_PFT, DIM_SCALAR.

! stream/frequency bit positions -> a single membership bitmask per variable ---------------!
integer(ik), parameter :: FREQ_HOURLY = 1_ik, FREQ_DAILY = 2_ik, FREQ_MONTHLY = 4_ik,        &
                          FREQ_ANNUAL = 8_ik      ! ior() these; test with iand()
integer(ik), parameter :: N_FREQ = 4_ik          ! number of temporal tiers

type :: var_desc_t
   character(len=32)  :: name       = ''          ! netCDF variable name (e.g. 'agb')
   character(len=96)  :: long_name  = ''          ! CF long_name attribute
   character(len=24)  :: units      = ''          ! CF units attribute
   integer(ik)        :: dim        = DIM_SCALAR   ! DIM_* : the trailing (ragged) axis
   integer(ik)        :: agg        = AGG_MEAN     ! AGG_* : temporal reduction operator
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

`dim` is MEDS's clean replacement for ED2's magic `idim_type` integer (a ~40-case header comment): an
enumerated axis code, not an arithmetic tuple. `geth5dims`'s job (idim_type → rank+shape) becomes a
trivial `select case(v%dim)` in the stream writer (§5.3).

**State vs. flux — the `agg` operator must respect physical units (ED2's integrated-vs-mean split).**
ED2 deliberately separates *integrated fluxes* from *mean states*, and MEDS must too. A **state**
(`agb`, `soil_temp`, `dbh`) is a stock: its period value is a **time-mean**, and if steps are
non-uniform the mean must be **dt-weighted** — `AGG_TMEAN` accumulates `Σ(x·dt)` and normalizes by
`Σdt`, reducing to `AGG_MEAN` under uniform `dt`. A **flux** (`total_gpp` in `kgC/m²/s`, transpiration,
NEE) is a rate: a physically meaningful period *total* is the **integral** `Σ(x·dt)`, whose units are
the rate × time (`kgC/m²` over the period) — this is `AGG_FLUXSUM`. A plain `AGG_SUM` of per-step
samples is dt-unaware: under uniform `dt` it is off by exactly the step-count factor (wrong magnitude
*and* wrong units for a "total"), and across tiers with different sub-step counts (hourly feeding daily)
a plain sum is not a coherent integral. `AGG_SUM` is therefore reserved for **count-like** integer
tallies (e.g. number of mortality events) where dt-weighting is meaningless; every physical rate uses
`AGG_FLUXSUM` and every physical stock uses `AGG_TMEAN`/`AGG_MEAN`. Each `var_desc_t` thus declares, via
its `agg`, whether it is a stock or a flux; `harvest`/`accumulate` carry the step `dt` (§4.3). The
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
field supplies its value, and a single **`harvest_site`** switchboard copies the current value out each
step:

```fortran
! meds_output_accum.f90 -- pure extraction: current instantaneous value(s) of one variable.
subroutine harvest_var(site, v, scalar_out, slab_out, n_out)
   type(site_t),     intent(in)  :: site
   type(var_desc_t), intent(in)  :: v
   real(wp),         intent(out) :: scalar_out          ! for DIM_SCALAR
   real(wp),         intent(out) :: slab_out(:)         ! for DIM_COHORT/PATCH/... : the live [1:n] prefix
   integer(ik),      intent(out) :: n_out               ! valid length this step (ragged, §4.4)
   ! select case (v%source_id): copy site%cohort%agb(1:n), or call total_agb(site), etc.
end subroutine
```

Because `harvest_var` **reads and copies** (never aliases), a cohort birth/death/fusion between steps is
invisible to the output layer — there is **no `filltab_alltypes` re-hash** and no
hidden-global-mutable-state. This is the load-bearing divergence from ED2 (Bundle-style takeaway: pass an
explicit descriptor, re-derive from state, do not cache pointers).

### 3.4 The DRY registration pattern

One source list — a single `contains`-local block in `build_output_registry` — declares every variable
once. A tiny `add()` helper appends a `var_desc_t`; the frequency membership and default operator are
literal arguments. This is **one of the two** "add a variable" edit sites (the other is the matching
`harvest_var` case, §3.3): the `add()` line names the variable and its metadata, and the `SRC_*` case
copies the right SoA field or calls the right reduction. The two are **deliberately co-located** — the
`SRC_*` parameter, the `add()` line, and the `harvest_var` case are grouped by dimension so a new
cohort-dimensioned variable is a three-line diff in one region of two adjacent files, never a hunt
across `io_create` + `io_write_snapshot` + a struct field (today's L5). It is **two coupled edits, not
one** — the irreducible price of rejecting ED2's pointer registry — but pointer-free and local:

```fortran
subroutine build_output_registry(reg, cfg)
   type(output_registry_t), intent(out) :: reg
   type(meds_config_t),     intent(in)  :: cfg
   integer(ik) :: k
   allocate(reg%var(MAX_OUTPUT_VARS)) ; reg%nvar = 0_ik
   ! name           long_name                  units       dim         agg       streams(default)          src
   call add('agb',      'aboveground biomass',     'kgC/m2',  DIM_COHORT, AGG_MEAN, DAY_MON_YR,  SRC_C_AGB)
   call add('nplant',   'plant number density',    'plant/m2',DIM_COHORT, AGG_LAST, DAY_MON_YR,  SRC_C_NPLANT)
   call add('dbh',      'diameter at breast hgt',  'cm',      DIM_COHORT, AGG_MEAN, MON_YR,      SRC_C_DBH)
   call add('total_agb','site aboveground biomass','kgC/m2',  DIM_SCALAR, AGG_MEAN, ALL_FREQ,    SRC_S_AGB)
   call add('total_lai','site leaf area index',    'm2/m2',   DIM_SCALAR, AGG_MEAN, ALL_FREQ,    SRC_S_LAI)
   call add('total_gpp','site GPP',                'kgC/m2/s',DIM_SCALAR, AGG_MEAN, HOUR_DAY_MON,SRC_S_GPP)  ! P1, fast loop
   ! ... (full table in §3.6) ...
   call apply_toml_overrides(reg, cfg)      ! flip enabled / edit streams from [output.variables]
   call build_freq_index(reg)               ! precompute idx_freq
contains
   subroutine add(nm, ln, un, dm, ag, st, src)
      ! ... reg%nvar = reg%nvar+1 ; reg%var(reg%nvar) = var_desc_t(nm, ln, un, dm, ag, st, .true., ...)
   end subroutine
end subroutine
```

The `DAY_MON_YR`, `ALL_FREQ`, … are named `ior` combinations (`integer(ik), parameter`) — a readable
default membership per variable, all overridable from TOML.

### 3.5 Sample variable-registry table (P0 set + P1 additions)

| name | long_name | units | dim | agg | default streams | source |
|---|---|---|---|---|---|---|
| `nplant` | plant number density | plant/m² | cohort | last | D M Y | `cohort%nplant` |
| `dbh` | diameter at breast height | cm | cohort | mean | M Y | `cohort%dbh` |
| `height` | height | m | cohort | mean | M Y | `cohort%height` |
| `basal_area` | basal area per plant | cm²/plant | cohort | mean | Y | `cohort%basal_area` |
| `agb` | aboveground biomass per plant | kgC/plant | cohort | mean | D M Y | `cohort%agb` |
| `leaf_area` | leaf area per plant | m²/plant | cohort | mean | M Y | `cohort%leaf_area` |
| `growth_avg` | moving-avg growth (mort. predictor) | cm/yr | cohort | mean | M Y | `cohort%growth_avg` |
| `pft` | plant functional type index | – | cohort | last | D M Y | `cohort%pft` (int) |
| `owner_patch` | owning patch (1-based) | – | cohort | last | D M Y | `cohort%owner_patch` (int) |
| `global_cohort_id` | persistent cohort id | – | cohort | last | D M Y | `cohort%global_id` (int) |
| `patch_area` | patch area fraction | – | patch | last | D M Y | `patch%area` |
| `patch_age` | time since disturbance | yr | patch | last | M Y | `patch%age` |
| `dist_type` | disturbance type | – | patch | last | M Y | `patch%dist_type` (int) |
| `cohort_offset` | first cohort of patch (CSR) | – | patch | last | D M Y | `patch%cohort_offset` (int) |
| `cohort_count` | cohorts in patch (CSR) | – | patch | last | D M Y | `patch%cohort_count` (int) |
| `global_patch_id` | persistent patch id | – | patch | last | D M Y | `patch%global_id` (int) |
| `total_nplant` | site total plant number | plant/m² | scalar | mean | D M Y | `total_nplant(site)` |
| `total_basal_area` | site total basal area | m²/m² | scalar | mean | D M Y | `total_basal_area(site)` |
| `total_agb` | site aboveground biomass | kgC/m² | scalar | mean | D M Y | `total_agb(site)` |
| `total_lai` | site leaf area index | m²/m² | scalar | mean | D M Y | `total_lai(site)` |
| `mean_dbh` | basal-area-weighted mean DBH | cm | scalar | mean | D M Y | `mean_dbh(site)` |
| `n_cohort` | cohorts in use | – | scalar | last | D M Y | `site%cohort%n` (int) |
| `n_patch` | patches in use | – | scalar | last | D M Y | `site%patch%n` (int) |
| — *P1 additions (need the fast loop / new diagnostics)* — |||||||
| `gpp_total` | site GPP time-integral | kgC/m² (per period) | scalar | **fluxsum** | H D M | `Σ gpp·dt` |
| `gpp_rate` | site GPP mean rate | kgC/m²/s | scalar | **tmean** | H D M | `gpp` (flux, dt-wtd) |
| `gpp_rate` (variance twin) | — | — | scalar | **meansq** | M | — |
| `soil_temp` | soil temperature by layer | K | soil | tmean | H D M | `patch%soil_e%soil_temp` |
| `soil_water` | volumetric soil moisture | m³/m³ | soil | tmean | D M | `patch%soil_w%theta` |
| `mort_rate` | realized mortality intensity | 1/yr | cohort | tmean | M Y | rate array |
| `max_dbh` | largest DBH this period | cm | scalar | **max** | M Y | `maxval(dbh)` |

(**last** for identity/index/CSR fields — an average of `global_id` or `cohort_offset` is meaningless;
they carry the *end-of-period* snapshot needed to interpret the ragged slabs, §4.4. **tmean** is the
dt-weighted state mean and is the default for every physical stock — it collapses to `mean` under the
uniform slow-step `dt` but is unit-correct once the hourly tier's `dt_fast` sub-steps feed it;
**fluxsum** integrates a rate to a period total, dropping the `/s` in the on-disk `units`; **sum** is
reserved for integer count tallies, §3.5-note above.)

### 3.6 Instantaneous escape hatch

A variable whose default streams include no averaging tier and whose `agg = AGG_LAST` on a per-step
frequency reproduces today's *values* exactly (instantaneous, written at flush). So the P0 registry with
everything at `AGG_LAST` on the annual stream carries the **same field values** as the current
`-D-output.nc` — the migration is a strict generalization.

**But not bit-for-bit, and the doc does not claim it.** The accumulate path flushes ragged slabs in
**first-seen-`global_id` order** (§4.4), whereas `io_write_snapshot` today writes **live-SoA slot
order** — i.e. cohorts *height-sorted* by the demography reorder. So even at `AGG_LAST` the rows are a
**permutation** of today's, and the two `global_id` slabs are the key that re-pairs them. Add
per-variable `cell_methods`/`_FillValue` attributes (§5.3) and the files differ in metadata too. The
round-trip test (§8, test 5) therefore compares as an **id-keyed set** (join on `global_id`, assert
equal values), *not* byte-for-byte — "equivalent up to row permutation and metadata", which is the
honest and testable claim.

---

## 4. The aggregation engine

### 4.1 Temporal tiers and the chaining rule

Four accumulation tiers, mirroring ED2's `fmean → dmean → mmean → qmean` but **cadence-driven by the
existing stepper flags** rather than second-counting:

| Tier | Fed from (per variable) | Integrated every | Normalized & flushed at | ED2 twin |
|---|---|---|---|---|
| **hourly** (`hmean`) | raw state (if finest active) | `dt_fast` sub-step (P1) | hour roll-over | `fmean` |
| **daily** (`dmean`) | raw state, *or* the hourly mean if hourly is active | slow step / hour roll-over | day roll-over | `dmean` |
| **monthly** (`mmean`) | raw state, *or* the next-finer active tier | slow step / day roll-over | month roll-over | `mmean` |
| **annual** (`ymean`) | raw state, *or* the next-finer active tier | slow step / month roll-over | year roll-over | yearly vars |

**The feeder rule, stated concretely (per variable, not per tier).** For each variable, precompute its
ordered list of active tiers from `streams ∧ enabled` (call the finest `t₀`). **The finest active tier
`t₀` accumulates from raw state every step** via `harvest_var` — this is the load-bearing correction:
*whatever* the finest active tier is, it reads state directly, so a **monthly-only** or **annual-only**
config works (the monthly buffer accumulates every slow step, not off a non-existent daily feeder).
**Every coarser active tier chains from the next-finer active tier** at that finer tier's roll-over —
never from raw state, and never from an *inactive* intermediate tier. So daily+annual (no monthly)
chains annual off daily; monthly-only reduces raw state straight into the monthly buffer. The per-var
"who feeds whom" is stored once in `idx_freq`/a `feeder(tier)` map at registry build (§3.2), so
`output_accumulate` folds raw state into the **precomputed finest tier per variable**, never a
hard-coded daily. **Configs to test explicitly: monthly-only, annual-only, and daily+annual-skip-monthly**
(§8, test 2).

**Chaining exactness.** For `AGG_MEAN`/`AGG_TMEAN` the chain is exact (a dt-weighted mean of dt-weighted
means, carrying `Σdt` as the chained weight, equals the direct dt-weighted mean over the coarse period).
For `AGG_FLUXSUM` the chain is a plain sum of sub-period integrals (integrals add). For
`AGG_SUM/MIN/MAX/LAST` the operator re-applies up the chain (sum of sums, min of mins, …), exact.
**`AGG_MEANSQ` chains on the lower tier's *mean*, not its meansq** — variance is
`mmean(x²)_from_finer_means − (mmean(x))²`, the diel/inter-day variance, matching ED2's `mmsqu`
second-moment tier. Each chained buffer carries its weight `Σdt` (not a bare sample count) so the
roll-up stays dt-correct even when sub-tiers have unequal step counts.

### 4.2 The accumulator buffer

Per active `(variable, tier)` pair, one `accum_buffer_t` holds the running reduction, the sample count,
and the **dt weight**. Two addressing modes share the type, keyed by `dim`: **id-keyed** ragged
(`DIM_COHORT`/`DIM_PATCH` — born/die/fuse mid-window, indexed by persistent `global_id`, §4.4) and
**fixed-index** (`DIM_SOIL`/`DIM_PFT` — a stable axis with no birth/death, indexed directly, §4.4). Both
are dimensioned to the relevant cap (`io_cohort_max`, `io_patch_max`, `n_soil_layer`, `n_pft`).

```fortran
type :: accum_buffer_t
   integer(ik) :: var_id  = 0_ik          ! index into registry%var(:)
   integer(ik) :: freq    = 0_ik          ! FREQ_* this buffer serves
   integer(ik) :: agg     = AGG_MEAN
   integer(ik) :: dim     = DIM_SCALAR
   logical     :: id_keyed = .false.      ! .true. for DIM_COHORT/PATCH (global_id map); .false. fixed-index
   real(wp)    :: scal    = 0.0_wp        ! scalar accumulator (DIM_SCALAR)
   real(wp)    :: scal2   = 0.0_wp        ! scalar second moment  (AGG_MEANSQ, DIM_SCALAR)  <-- was missing (bug fix)
   real(wp)    :: wsum    = 0.0_wp        ! scalar dt weight Sum(dt) (AGG_TMEAN/FLUXSUM/MEANSQ); nsamp for count-aggs
   real(wp),    allocatable :: slab(:)    ! ragged/axis accumulator, sized to cap
   real(wp),    allocatable :: slab2(:)   ! second-moment accumulator (AGG_MEANSQ only)
   real(wp),    allocatable :: wsum_slab(:)! per-entity dt weight Sum(dt) (ragged TMEAN/FLUXSUM)
   integer(ik), allocatable :: hits(:)    ! per-entity sample count (ragged: entities appear/vanish mid-window)
   integer(ik), allocatable :: id_of_row(:)! global_id occupying each slab row (id-keyed only); 0 = empty
   integer(ik) :: n_active = 0_ik         ! number of distinct ids seen this window (id-keyed)
   integer(ik) :: nsamp   = 0_ik          ! scalar sample count
   real(wp)    :: seed    = 0.0_wp        ! MIN=+huge, MAX=-huge, SUM/MEAN/FLUXSUM=0, LAST=undefined
end type accum_buffer_t
```

The scalar `scal2` is a real field (it was referenced by the `AGG_MEANSQ` kernel in an earlier draft but
never declared); `wsum`/`wsum_slab` carry the `Σdt` weight that makes `AGG_TMEAN`/`AGG_FLUXSUM` and the
chained roll-up dt-correct (§4.1, §4.3).

### 4.3 The three kernels (integrate / normalize / reset)

A **single generic per-variable loop** replaces ED2's ~8700 lines of hand-written per-field boilerplate.
`accumulate_step` folds one harvested value into the buffer by operator; `normalize_buffer` closes the
period; `reset_buffer` re-seeds it. All `pure`/`elemental`-friendly (arithmetic + intrinsics):

```fortran
elemental subroutine accumulate_scalar(buf, x, dt)      ! one entity, one step of length dt [s]
   type(accum_buffer_t), intent(inout) :: buf
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
   ! ragged: hits=0, wsum_slab=0, id_of_row=0, n_active=0 (a fresh id-window each period)
end subroutine
```

**Empty-period normalization is defined, not left to chance.** `normalize_buffer` returns a `valid`
flag; when a tier saw **zero samples** in a period (an hourly tier enabled but the fast loop didn't run,
or a ragged slot no entity ever occupied), `valid=.false.` and the serializer writes the CF
`_FillValue` for that record/slot — never a `0/0` NaN, and never the raw `±huge` MIN/MAX seed. The
unused `[n+1:cap]` ragged fill tail is likewise written as `_FillValue`, not zero (§5.3). Test 1 covers
the zero-sample scalar and the never-occupied ragged slot.

The **id-keyed** ragged variant (`accumulate_slab`, `DIM_COHORT`/`DIM_PATCH`) loops over the live
`[1:n]` prefix and folds each entity into its `global_id`-addressed row (§4.4), incrementing that row's
`hits` and `wsum_slab`; `normalize` divides `AGG_MEAN` by the row's **own** `hits` and `AGG_TMEAN` by
its **own** `wsum_slab` (an entity present for only part of the window is averaged over the steps/seconds
it existed, matching a demographic mean-while-alive). The **fixed-index** variant (`DIM_SOIL`/`DIM_PFT`)
skips the id map entirely — soil layers and PFT bins are a stable axis, folded directly by index `k`
(§4.4) — so no hash, no first-seen ordering, no cap-overflow path.

### 4.4 Ragged cohort/patch dims across a window — the hard part

Cohorts and patches are **born, fused, split, and culled within a window**; their SoA slot index is not
stable, so accumulating "cohort slot 7" across a month is meaningless. MEDS already has the fix in
hand: the **persistent `global_id`** (monotonic, never reused, carried in lockstep through every
sort/fusion). The window accumulator is **keyed on `global_id`, not slot**:

- The ragged `slab` buffer is a **dense array sized to the cap**, plus a companion
  `id_of_slot(cap)` filled at flush and an `slot_of_id` hash (or a sorted-id binary search) built once
  per step from the live `global_id(1:n)`. `accumulate_slab` maps each live cohort to its stable buffer
  row; `hits` counts how many steps that id contributed.
- **Birth mid-window** → a new `global_id` claims a fresh buffer row (seeded on first sight);
  `hits < nsteps`, so its mean is over its lifetime. **Death/fusion mid-window** → the id stops
  contributing; its partial mean is still valid and is flushed. **Fusion** keeps the survivor id (its
  row keeps accumulating); the absorbed id's row closes — a faithful record.
- **Cap overflow** is handled more gracefully than today's hard `error stop`: if the count of distinct
  ids seen in a window exceeds the cap, the flush writes the cap-many highest-`agb` (or first-seen) ids
  and logs a one-line warning, rather than aborting a long run. (The `error stop` remains available via
  a strict-mode config flag for regression runs.)

**The flush membership contract — instantaneous and averaged are NOT the same, and the doc says so.**
The naïve "serialize exactly like today's ragged prefix + CSR" is **wrong for any averaged window with
churn**: the averaged slab spans `[1:n_active]` = *all* ids seen (a superset that includes mid-window
dead/absorbed cohorts), whereas a `cohort_offset`/`cohort_count`/`owner_patch` CSR written `AGG_LAST`
describes only the *surviving, end-of-period compact* ordering. Those two orderings differ, so a reader
that slices the slab by the end-of-period CSR reads the wrong rows, and dead-id rows have no patch
mapping at all. MEDS resolves this by **membership semantics that depend on the stream's `agg`**:

- **Instantaneous ragged streams** (all vars `AGG_LAST`, the today-compatible mode): the slab **is** the
  live compact `[1:n]` set in slot order, and the **CSR contract is retained unchanged** —
  `cohort_offset`/`cohort_count`/`owner_patch` (all `AGG_LAST`) plus the `global_id` slab reconstruct
  membership exactly as today's reader does. Nothing changes for this path.
- **Averaged ragged streams** (any tier-mean/sum/flux var): the flush writes a **full per-id membership
  map** — an `owner_patch` slab **one entry per slab row, in first-seen-`global_id` slab order**
  (`AGG_LAST` per id = the last patch that id belonged to) — and **drops the compact CSR
  `offset`/`count`** for that stream. The reader keys membership directly off the paired
  `(global_id, owner_patch)` slabs, so a mid-window-dead cohort's partial mean is fully interpretable
  (it carries its own id and its own last patch). The CSR prefix-sum contract is an
  *instantaneous-snapshot* device and is not reused for time-averaged records.
- **Config `ragged_avg_alive_only`** (default `false`): if set, an averaged ragged stream instead emits
  **only ids alive at flush**, in which case the compact CSR *is* snapshot-consistent and the simpler
  today-reader applies — at the cost of losing mid-window-dead cohorts' partial means. The default keeps
  the fuller per-id-map record; the flag is the "snapshot-consistent CSR" opt-in.

**Fixed-index dims skip all of this.** `DIM_SOIL`/`DIM_PFT` have a **stable, birth/death-free** axis, so
they accumulate by direct index `k` and serialize the full `[1:n_soil]`/`[1:n_pft]` axis with no
`global_id` slab, no CSR, no cap-overflow path — the id-keyed machinery above is **restricted to
`DIM_COHORT`/`DIM_PATCH`**.

**Known semantic — fusion value-discontinuity inside a mean window (documented, not fixed).** When a
surviving-id cohort absorbs another via fusion mid-window, its `dbh` (and derived geometry) is
**re-derived from the conserved total AGB** (CLAUDE.md invariant), so the survivor's per-plant series
takes a step *inside* the averaging window — the resulting `AGG_TMEAN` mixes two physical compositions
(pre- and post-fusion) under one id's mean. This is inherent to averaging a demographically mutating
state and is left as a **flagged limitation**: conserved *extensive* quantities (AGB, plant number) are
still exact; only per-plant *intensive* means (`dbh`, `height`) carry the blend. A user needing
fusion-clean series uses a finer tier or the instantaneous stream.

This "accumulate by persistent id, flush the per-id membership map" scheme is the key algorithm that
lets a *ragged* demographic state be *time-averaged* — ED2 sidesteps it by re-hashing pointers every
birth/death; MEDS uses the id it already maintains.

### 4.5 Coupling to the stepper cadence

The accumulators tick from **`meds_stepper`** (the cadence owner), fed the same `is_new_month`/
`is_new_year` flags `meds_main` already computes, plus a new `is_new_day` (trivially `now%day /=
prev%day`) and (P1) `is_new_hour`, and the step length `dt` (for dt-weighting, §4.3). The stepper's
per-step entry is **netCDF-free**: it accumulates and *stages* closed records; it does **not** write.

```fortran
subroutine output_accumulate(mgr, site, now, dt, is_new_hour, is_new_day, is_new_month, is_new_year)
   type(output_manager_t), intent(inout) :: mgr           ! netCDF-free plain data (buffers + pending stage)
   type(site_t),           intent(in)    :: site
   type(meds_time_t),      intent(in)    :: now
   real(wp),               intent(in)    :: dt            ! [s] length of the step just taken
   logical,                intent(in)    :: is_new_hour, is_new_day, is_new_month, is_new_year
   ! 1. accumulate: for each var, fold harvest_var(site,v) into its PRECOMPUTED FINEST active tier
   !                (mgr%feeder(:) from §4.1 — NOT hard-coded to daily; carries dt for TMEAN/FLUXSUM)
   ! 2. on is_new_hour  -> normalize+STAGE hourly records; chain hourly-mean into the daily buffer; reset hourly
   ! 3. on is_new_day   -> normalize+STAGE daily  records; chain into monthly (if monthly active) else annual; reset daily
   ! 4. on is_new_month -> normalize+STAGE monthly records; chain into annual; reset monthly
   ! 5. on is_new_year  -> normalize+STAGE annual  records; reset annual
   ! all of the above is plain-array arithmetic on mgr; the STAGED records sit in mgr%pending until main drains them
end subroutine
```

Then, back in `meds_main`, once per step:

```fortran
call output_serialize_pending(mgr)   ! the ONLY netCDF-touching call; no-op stub when MEDS_ENABLE_IO=OFF
```

`output_serialize_pending` walks `mgr%pending`, opens/rolls the right per-stream file (§5.2), writes each
staged record, and clears the stage. Chaining feeds the **coarsest-below** *active* tier (the §4.1 feeder
map), so daily→annual-skip-monthly and monthly-only both route correctly. The roll-over tests fire **at
the boundary**, closing the period that just ended (ED2 names files by the *previous* period; §5.2), and
ordering is strict low→high tier so each chain reads a freshly-normalized feeder (the ED2 `ed_output`
ordering, distilled). P0 keeps the accumulate step **host-only** (a cheap copy-reduce over the live
prefix, negligible beside the daily kernels); if profiling ever shows it matters, the slab accumulate is
a `target teams distribute` sweep over plain arrays like `growth_step` — designed-in, not needed at P0.

**Where it wires:** `advance_one_step` gains an optional `mgr` argument; after the slow loop mutates
state it calls the netCDF-free `output_accumulate` (the analogue of how it already optionally runs the
fast loop before `vegetation_dynamics`), and `meds_main` calls `output_serialize_pending` after
`advance_one_step` returns — the clean split of §2 that keeps `aux` netCDF-free and the flush at `main`.

---

## 5. File organization & performance

### 5.1 One file per stream per time-chunk (goal G4, FILE SIZE)

Never one lump. Each **stream** (frequency) writes to its own family of files, **time-chunked** so no
single file grows unbounded:

```
<output_dir>/<prefix>-H-<YYYYMMDD>.nc     hourly   stream, one file per SIM-DAY   (24 records)
<output_dir>/<prefix>-D-<YYYYMM>.nc       daily    stream, one file per SIM-MONTH (28-31 records)
<output_dir>/<prefix>-M-<YYYY>.nc         monthly  stream, one file per SIM-YEAR  (12 records)
<output_dir>/<prefix>-Y-<YYYY0>.nc        annual   stream, one file per SIM-DECADE(10 records) [or whole run]
<output_dir>/<prefix>-S-<YYYYMMDDHHMMSS>.nc   STATE/restart — UNCHANGED, one file per checkpoint
```

The stream letter (`H/D/M/Y`) mirrors ED2's `I/D/E/Q/Y` convention; the date stamp granularity is one
tier **coarser** than the record cadence (a per-day-of-records file for the hourly stream, a
per-month-of-records file for the daily stream, …), so each file holds a tidy, bounded number of
records (24 / ~30 / 12 / 10). The **time-chunk size is a per-stream TOML knob** (`file_chunk`, §6) so a
user can pick "one daily file per year" (365 records) instead of per month. Each file has its own
UNLIMITED `time` dim; records append until the chunk boundary, then a new file opens. This is ED2's
multi-record `out_time_fast`/`nrec`/`irec` splitting, made a plain rollover on the calendar flags.

### 5.2 Naming & rollover

`stream_write_record` compares the record's calendar bucket to the open file's bucket; on mismatch it
closes the current file and opens the next (`makefnam`-analogue: `<prefix>-<letter>-<stamp>.nc`, stamp
granularity from `file_chunk`). Files are named for the **period whose data they hold** (the daily
stream's `-D-202601.nc` holds January's daily means). Rollover is driven purely by `now`, so it is
exact across leap years (reusing `meds_time`).

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

The `1×max` chunk + deflate that already compresses the ragged fill tail carries over verbatim (the
unused `[n+1:cap]` tail of every record is written as **`_FillValue`**, not zero — §4.3 — and deflates
near-free). Variance variables get a `<name>_variance` companion (same dims); `AGG_MEANSQ` writes both
mean and variance from one buffer. The per-record `global_id` slabs + `n_cohort`/`n_patch` counts +
calendar are written for every stream. Membership is reconstructed per the §4.4 contract: the
instantaneous stream from the CSR, an averaged stream from the paired `(global_id, owner_patch)` per-id
map — each file is self-describing and independently readable.

**CF metadata for aggregated data (goal G7) — a reader must be able to tell a mean from a snapshot.**
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
serializer writes for an empty period (`valid=.false.`, §4.3) and for the unused ragged fill tail — so a
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

1. **In-memory accumulation** — the daily kernels touch **no disk**; they only fold values into RAM
   buffers (a copy-reduce over the live `[1:n]` prefix, O(n) per step, ≪ the growth/mortality/fusion
   work). The current design already writes at most once a year; this keeps per-step disk touches at
   **zero**.
2. **Buffered, infrequent flush** — a flush writes one closed period's records (a day's / month's worth)
   in a burst, chunked+deflated, then returns. Flush frequency is the *output* cadence, not the *step*
   cadence: a daily stream flushes once per sim-day regardless of how many sub-steps ran. Between
   flushes the netCDF file handle stays open (append into the current time-chunk), so there is no
   per-record open/close cost. **This is the recommended P0/P1 default** — synchronous flush, because a
   flush is rare and small relative to the compute between flushes, so it never throttles.
3. **Double-buffering + asynchronous writer (P2, optional, build-gated)** — for the *hourly* stream on a
   fast-loop run (where flushes are frequent and a synchronous HDF5 write could stall the step), the
   accumulator can hold **two buffer sets**: the stepper fills buffer A while a **background writer**
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
   - Tradeoff: one extra buffer set (memory ≈ 2× the tiny accumulator footprint), a lock/flag handshake,
     and a **run-end drain barrier**.

   **Recommendation: ship synchronous (options 1–2) at P0/P1, and treat async as an optional P2 that may
   well be dropped.** Because daily/monthly/annual flushes are already far off the critical path (§5.5
   intuition below), async only ever earns its keep for a *dense hourly* stream on a fast-loop run — and
   even then the cheaper lever is to **coarsen the hourly stream's file-chunk / variable set** rather
   than add threading. If P2 async is built, it lands on the nvfortran `-mp` build (or the pthread
   variant), never as "netCDF non-blocking on HDF5".

**Durability of long-held open handles (the cost of keeping the file open across a chunk).** Options 2–3
keep the HDF5 handle open across a whole time-chunk with no `nc_sync`, so a crash **loses the entire
in-flight chunk** — up to a sim-year for the annual stream, a sim-month for the daily. The design adds an
optional **`sync_every` knob** (`[output].sync_every = "flush" | "N records" | "never"`, default
`"flush"`): `nc_sync` after each flush (durable, one fsync per rare flush — negligible) or every N
records or never (fastest, crash-forfeits the chunk). A new `nc_sync` C binding is the only addition
needed (§7.2). On restart the diagnostic streams start **fresh files** anyway (open Q4), so a crash at
worst forfeits the current open chunk, never corrupts earlier ones.

Quantified intuition: the accumulator footprint is `Σ_tiers (nvars × cap)` doubles — order
(4 tiers × ~25 vars × ~4096 cohorts × 8 B) ≈ a few MB, trivially RAM-resident. A daily-stream flush
writes ~25 vars × (1 record × cap, deflated) ≈ tens of KB on disk — microseconds of write beside a
sim-day of daily kernels. The synchronous default is comfortably off the critical path for everything
except a hypothetical dense hourly stream, which is exactly what P2's async covers.

---

## 6. TOML config surface

### 6.1 Design

The current `[io]` block (8 mandatory presence-mapped keys) is **superseded** by an `[output]` block
with **per-frequency enable + file-chunk** knobs and an **`[output.variables]` per-variable toggle
map**. Two mechanisms are needed from the TOML reader, which today only does flat `key=value` + real
arrays (limitation L8):

- **Per-frequency enable** — plain booleans/strings under `[output]`, already expressible.
- **Per-variable toggle** — a `[output.variables]` **sub-table of `name = bool` (or `name = "D M Y"`
  stream-string)** entries. The flat `section.key=value` parser stores these as
  `output.variables.<name> = <value>` keys directly (no nested-table machinery needed — the existing
  parser already stores fully-qualified `section.key`; a `[output.variables]` header just sets the
  section prefix). So the **only reader change is to accept an arbitrary key set under one section**,
  which `toml_has`/`toml_logical`/`toml_string` already support by name. A per-variable *frequency
  string* ("D M Y") needs the string form; `toml_string` already exists. **No string-array or
  arrays-of-tables support is required** — the sub-table-of-scalars form sidesteps L8 entirely.

**Value grammar (specified, not sniffed loosely).** A `[output.variables].<name>` value is one of:
(a) a **bool** — `true` force-enable at the registry default streams, `false` disable everywhere; or
(b) a **stream-string** — a space-separated subset of the tokens `{H, D, M, Y}` (e.g. `"M Y"`), parsed
by `parse_stream_mask(s) → ior(FREQ_*)` (unrecognized token → `error stop` naming the offender). The
loader disambiguates by TOML value *type* (the flat parser retains each stored value's scalar type), so
a bare `true`/`false` reads as (a) and a quoted string as (b); a value that is neither (e.g. an integer)
is a config error. `file_chunk` is likewise a **named-bucket string** mapped by
`parse_file_chunk(s) → {FC_DAY, FC_MONTH, FC_YEAR, FC_DECADE, FC_RUN}` filling the
`file_chunk(N_FREQ)` integer array — the enum the config type comments as `...` today is spelled out
here.

**Unknown/typo'd key validation (closing the silent-ignore trap).** Because `build_output_registry`
queries `output.variables.<name>` only for names it already knows, a typo (`grwoth_avg = false`) would
today be **silently ignored**. The flat parser *does* store every key under the section as an
enumerable array, so after building the registry the loader **iterates the stored
`output.variables.*` keys and `error stop`s (strict) or warns (lenient, `[options]`-gated) on any key
that does not match a registry variable name** — listing the offenders, exactly as `load_meds_config`
already lists missing mandatory keys. The same enumerate-and-check guards unknown keys under `[output]`
itself (a mistyped `file_chnk`). Versatility without a silent-failure UX trap.

Backward-compat: `load_meds_config` keeps reading the legacy `[io]` keys if `[output]` is absent
(mapping `output_interval_years` → an annual stream with all vars `AGG_LAST`), so existing TOML runs
unchanged. New keys use **optional-with-default** loaders (a new `opt_l`/`opt_s`/`opt_i` alongside the
existing `req_*`), because per-frequency/per-variable knobs must default sensibly rather than force the
user to name every variable — a deliberate softening of the "all `[io]` keys mandatory" rule for the
output block only.

### 6.2 Example `[output]` block

```toml
[output]
enabled       = true            # master switch (replaces io.write_output)
dir           = "out"           # replaces io.output_dir
prefix        = "run01"         # replaces io.output_prefix
cohort_max    = 4096            # ragged-dim caps (replaces io.cohort_max / io.patch_max)
patch_max     = 256
strict_caps   = false           # false: warn+truncate on overflow; true: error stop (regression mode)
sync_every    = "flush"         # nc_sync durability: "flush" | "<N> records" | "never"  (§5.5)
ragged_avg_alive_only = false   # averaged ragged: false=full per-id membership map; true=alive-at-flush CSR (§4.4)

  # --- per-frequency streams: enable + how many records per file (the time-chunk) ---
  [output.hourly]
  enabled    = false            # needs the fast loop; off by default
  file_chunk = "day"            # one file per sim-day (24 records)

  [output.daily]
  enabled    = true
  file_chunk = "month"          # one file per sim-month

  [output.monthly]
  enabled    = true
  file_chunk = "year"           # one file per sim-year

  [output.annual]
  enabled    = true
  file_chunk = "run"            # one file for the whole run (few records)

  # --- per-variable overrides: turn a variable OFF, or re-target its streams ---
  #     absent -> the registry default (§3.5). value = bool OFF/ON, or a stream string.
  [output.variables]
  growth_avg       = false      # never write this variable
  basal_area       = "M Y"      # override default membership: monthly + annual only
  soil_temp        = "H D"      # hourly + daily (P1)
  total_gpp        = true       # keep registry default streams, force-enabled

[state]                          # the restart stream stays its own block (unchanged semantics)
write_state          = true
state_interval_years = 10
```

### 6.3 Config type additions

`meds_config_t` gains an `[output]` sub-record (a plain-data type, no io dependency at the DAG root):

```fortran
type :: output_config_t
   logical            :: enabled     = .false.
   character(len=256) :: dir = '.', prefix = 'meds'
   integer(ik)        :: cohort_max = 4096_ik, patch_max = 256_ik
   logical            :: strict_caps = .false.
   logical            :: ragged_avg_alive_only = .false. ! averaged ragged: false=per-id map, true=alive-at-flush CSR (§4.4)
   integer(ik)        :: sync_every  = FC_FLUSH          ! nc_sync policy: FC_FLUSH | N | FC_NEVER (§5.5 durability)
   logical            :: freq_on(N_FREQ)   = [.false.,.true.,.true.,.true.]  ! H D M Y defaults
   integer(ik)        :: file_chunk(N_FREQ) = [FC_DAY, FC_MONTH, FC_YEAR, FC_RUN]  ! records-per-file bucket per tier
   ! per-variable overrides are read straight from the toml table by name in build_output_registry,
   ! so they are NOT stored here as a fixed struct (avoids a name/length limit) -- the registry
   ! re-queries cfg%toml (or a retained toml_table_t handle) for output.variables.<name>.
end type
```

The per-variable map is applied inside `build_output_registry` by looping the registry names and
querying `output.variables.<name>` — so the config type stays fixed-size and the variable list stays
source-defined (the registry is the single source of truth for *which* variables exist; TOML only
toggles them).

---

## 7. netCDF layer & the stub

### 7.1 Reuse of `meds_netcdf_c`

The C bindings (`nc_create/def_dim/def_var/put_att_text/def_var_deflate/def_var_chunking/enddef/
put_vara_double/put_vara_int/put_var1_double/close/open/inq_varid/get_vara_int/get_vara_double/
strerror`, plus `cstr`/`nc_check` and the `NC_*` params) are **sufficient as-is** for the whole design:
the stream writer defines dims/vars from the registry, chunks+deflates each ragged var, and appends
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

### 7.3 Structuring the stream writer for the stub

`meds_output_stream` holds a `stream_file_t` (one per active frequency) carrying the open `ncid`, the
current time-chunk bucket, `nrec`, and — crucially — a **registry-indexed `var_id` array**
(`integer(c_int), allocatable :: vid(:)`) instead of today's fixed named fields (`v_c_pft`, `v_c_agb`,
…). This is what makes it data-driven: `stream_open_file` loops `reg%idx_freq(:,freq)` and defines +
stores each var-id by registry index, so `stream_write_record` writes by the same index. Adding a
variable never touches this module (it reads the registry).

### 7.4 The stub twin (goal G6)

Following the `meds_io`/`meds_io_stub` pattern exactly: `meds_output_stream_stub.f90` provides a no-op
module of the **same name** (`meds_output_stream`), same public interface (`stream_open_file`,
`stream_write_record`, `stream_close_file`) — CMake compiles exactly one, selected by
`MEDS_ENABLE_IO`. Because the **registry and accumulator modules are netCDF-free** (they touch only
`site_t`, plain arrays, and the reductions — no C bindings), they are **always compiled**; only the
serializer is stubbed. The manager (`output_manager_t`) calls the serializer through the stubbable
interface, so with IO OFF the accumulators still run (testable, and cheap) and the flush is a no-op.
Keeping the two serializer twins interface-identical is the one hand-sync discipline (mitigated by a
CTest that compiles both and checks the public symbol set — §8), and the surface is smaller than today
because the whole variable list no longer lives in the interface.

---

## 8. Test plan

Unit tests (CTest targets, links `meds_io` or `-DMEDS_ENABLE_IO=OFF` for the netCDF-free ones; build the
**nvfortran multicore** back end too — issue #7):

1. **Accumulator arithmetic** (`test_output_accum`, netCDF-free): feed a known scalar sequence; assert
   `AGG_MEAN` = arithmetic mean, `AGG_TMEAN` on **non-uniform `dt`** = the dt-weighted mean (≠ plain
   mean), `AGG_FLUXSUM` = `Σ(x·dt)` (units check: a constant rate over T seconds → `rate·T`), `SUM` =
   sum, `MIN`/`MAX` = extrema, `LAST` = final, `MEANSQ`-derived variance = the textbook value; assert
   `reset_buffer` re-seeds (min→+huge, max→−huge, mean/wsum→0). **Zero-sample guard:** normalize a
   never-touched buffer → `valid=.false.`, output = `_FillValue`, **no NaN and no ±huge leak**.
2. **Chaining exactness & feeder routing** (netCDF-free): (a) drive 30 daily values → monthly `AGG_TMEAN`
   equals the direct dt-weighted mean of the 30, monthly variance = inter-day variance, annual = mean of
   12 monthly means; `AGG_FLUXSUM` chain = sum of sub-period integrals. (b) **Feeder configs:** run
   **monthly-only**, **annual-only**, and **daily+annual-skip-monthly** and assert each variable's
   finest active tier accumulates raw state directly and coarser tiers chain the right feeder — the
   monthly-only buffer is fed every step (not left empty), the daily+annual case chains annual off daily.
3. **Ragged window by `global_id`**: synthesize a window with a cohort **born**, one **culled**, and a
   **fusion** (survivor id continues); assert each id's mean is over *its* `hits`/`wsum_slab`, the
   absorbed id's partial mean flushes, and the **averaged per-id membership map** (`global_id` +
   per-row `owner_patch`) reconstructs each row's patch — including a mid-window-dead row (the
   contract that the instantaneous-CSR path would get wrong, §4.4). Assert the instantaneous
   `AGG_LAST` stream still reconstructs from the CSR unchanged, and `ragged_avg_alive_only=true` emits
   the alive-only compact CSR. Assert **cap overflow** warns+truncates (non-strict) and `error stop`s
   (strict). Assert **`DIM_SOIL`/`DIM_PFT` fixed-index** vars accumulate by direct index with no id map.
4. **Registry + TOML toggles** (`test_output_registry`): build the registry, apply an
   `[output.variables]` map (one `false`, one stream-string override), assert the disabled var is absent
   from every `idx_freq` and the retargeted var appears only in its named tiers. **Unknown-key trap:** a
   typo'd `[output.variables]` name → `error stop` (strict) / warning (lenient) naming the offender, not
   a silent no-op (§6.1).
5. **Round-trip / back-compat** (`test_output_roundtrip`, needs netCDF): run a short spin-up with the
   **P0 registry set to `AGG_LAST` on the annual stream**; assert the produced `-Y-*.nc` is
   **equal to today's `-D-output.nc` up to row permutation and metadata** — compare as an **id-keyed
   set** (join the two `global_id` slabs, assert equal values), *not* byte-for-byte, since the
   accumulate path flushes in first-seen-id order vs today's height-sorted slot order (§3.6). Then a
   daily-stream run: re-open, assert record counts match the calendar, rollover at month boundaries,
   `time`/`year`/`month`/`day` exact across a leap year, and that each averaged variable carries the
   right **`cell_methods`**, a **`time_bnds`** interval, and a **`_FillValue`** attribute (§5.3).
6. **File organization**: assert one file per stream per time-chunk (glob the output dir), that the
   state stream `-S-*.nc` is untouched, and that `strerror`-checked writes all return `NC_NOERR`.
7. **Stub parity**: compile with `-DMEDS_ENABLE_IO=OFF`; assert `output_accumulate` still runs (buffers
   hold the right values) and `output_serialize_pending` is a silent no-op; a CTest greps both serializer twins'
   `public ::` lists for equality (the drift guard).
8. **Performance smoke**: a multi-decade daily run; assert wall-time with output ON is within a few
   percent of output OFF (I/O off the critical path), and peak RSS grows by only the accumulator
   footprint.

---

## 9. Phasing

**P0 — registry + daily/annual means + per-variable toggle + split files (the core value).**
- `meds_output_types` (`var_desc_t`, `accum_buffer_t` incl. `scal2`/`wsum`/`wsum_slab`/`id_of_row`,
  `output_registry_t`, netCDF-free `output_manager_t`, `AGG_*`/`FREQ_*`/`DIM_*`/`FC_*` codes).
- `build_output_registry` with the P0 table (§3.5), `harvest_var`, the per-var **feeder map** (§4.1).
- Daily + annual tiers; `AGG_MEAN/TMEAN/SUM/LAST` operators (dt carried, empty-period `_FillValue`
  guard); scalar + ragged-by-`global_id` accumulate with the **averaged per-id membership map** (§4.4);
  `DIM_SOIL`/`DIM_PFT` fixed-index accumulate.
- `meds_output_stream` (+ stub twin): per-stream, per-time-chunk files, data-driven var-id array,
  **`cell_methods` + `time_bnds` + `_FillValue`** CF metadata, and **provenance** globals (config echo +
  git hash, §5.3).
- `[output]` TOML block + `[output.variables]` toggle map (value grammar + **unknown-key validation**);
  `opt_*` loaders; legacy `[io]` back-compat.
- `output_accumulate` (netCDF-free) wired into `meds_stepper`; `output_serialize_pending` driven from
  `meds_main` (the flush wall) — replaces the `is_new_year` write block.
- Tests 1, 3, 4, 5, 7. **Deliverable:** user-mutable daily+annual per-variable output in split files,
  time-averaged, CF-tagged, the current stream reproducible **up to row permutation + metadata** (§3.6).

**P1 — hourly + monthly tiers, min/max/variance, flux integrals, durability.**
- Monthly tier + full chaining; `AGG_MIN/MAX/MEANSQ` (variance companion vars); `AGG_FLUXSUM` flux
  integrals (fed by the fast loop when it lands).
- Hourly tier scaffolding consuming fast-loop diagnostics (`gpp`, CAS/soil column state) **when** the
  fast loop is wired into `meds_main`; the P1 registry rows (§3.5).
- `nc_sync` durability knob (`sync_every`, §5.5) + its one C binding.
- Tests 2, 6, and the variance/flux paths of 1.

**P2 — optional asynchronous writer (may be dropped).**
- Only if a dense hourly stream proves to stall the step: double buffer set + **host-thread** background
  flush (nvfortran `-mp` OpenMP task **or** an `iso_c_binding` pthread — **not** netCDF non-blocking,
  which is unavailable on the HDF5 path), run-end drain barrier. Gated on an OpenMP-enabled build.
- Optional `NC_FLOAT` C binding for lighter diagnostic output.
- Test 8 (async variant: assert no compute stall on the hourly stream).

Each phase is independently shippable and leaves the model runnable; P0 alone removes limitations
L1–L6, L8 and gives the user everything in the stated goals except sub-daily/flux/variance (P1) and the
last slice of write-latency hiding (P2, optional).

---

## 10. Open questions

1. **Area-weighted spatial roll-up for scalar means.** The site totals are already area-weighted patch
   reductions (`meds_demography_diagnostics`); should any *cohort-dimensioned* variable also be offered
   as an nplant- or area-weighted site aggregate in the same registry (a `DIM_SCALAR` derived twin), or
   is that left to post-processing? Leaning: expose the existing five site totals as the canonical
   scalars; defer arbitrary weighted roll-ups to the reader.
2. **Weighting inside `AGG_MEAN` for cohorts that change `nplant` mid-window.** A simple per-step mean
   weights every step equally; an nplant- or agb-weighted temporal mean may be more ecologically
   meaningful for per-plant quantities. Proposed: default equal-weight (matches ED2's `dmean`); add an
   optional `AGG_WMEAN` (weighted by a companion source, e.g. nplant) in P1 if wanted.
3. **`file_chunk = "run"` for high-cadence streams.** A whole-run single file for the annual stream is
   fine (few records); for daily it defeats FILE-SIZE. Should the loader **reject** a `file_chunk`
   coarser than N records for sub-monthly streams, or just warn? Leaning: warn, trust the user.
4. **Restart interaction with open diagnostic files.** On restart mid-run (`INIT_RESTART`), the
   diagnostic streams start **fresh files** at the restart date (no append into a pre-restart file). Is
   that acceptable, or should the manager detect and append to the matching time-chunk file? Leaning:
   fresh files (simpler, and the calendar stamp disambiguates); document it.
5. **`global_id` cap vs. the accumulator cap.** The window may see more distinct ids than the flush cap
   even if the *instantaneous* count never exceeds it (churn). The truncation policy (§4.4) handles it,
   but the "right" cap for a long monthly window with heavy fusion needs a heuristic — is `2×io_cohort_max`
   a safe accumulator cap? Needs a spin-up measurement.
6. **Multi-site / polygon dimension.** The whole layout is single-site (no polygon axis), matching the
   engine. When multi-site lands, does a `site`/`polygon` dim slot in as a *leading* fixed dim on every
   var (cheap, uniform) — and does that interact with the per-time-chunk file split (one file per
   stream per chunk **per site**, or a site dim inside one file)? Deferred, but the registry/`dim`
   enum should reserve a `DIM_*` leading-axis convention now so it is not a breaking change later.
7. **Keep P2 async at all?** With daily/monthly/annual flushes already microseconds beside a sim-day of
   kernels, async only ever pays on a *dense hourly* stream — and even there, coarsening `file_chunk` or
   the hourly variable set is the cheaper lever. Leaning: **do not build async until a measured hourly
   stall demands it**, and if built, on the nvfortran `-mp` build (or a pthread), never as netCDF
   non-blocking. Revisit only after the fast loop + met forcing land and an hourly stream exists.
8. **Time stamp for a period mean.** Chosen (§5.3): `time` = period **midpoint** with a `time_bnds`
   `[start, end)` interval (CF-standard). Alternative conventions (stamp at start, or at end as ED2's
   "previous period" file naming implies) are equally decodable; the midpoint+bounds pair is picked so a
   plotter without bounds-awareness still places the mean sensibly. Confirm against the post-processing
   readers (`post_proc/`) before P0 freeze.
9. **`AGG_SUM` vs `AGG_FLUXSUM` classification per variable.** The registry must tag each rate as a flux
   (integrate, `FLUXSUM`) and each stock as a state (`TMEAN`), and reserve plain `SUM` for integer count
   tallies (§3.5). Is a hand audit of the ~25 P0 vars enough, or should `harvest_var` assert units
   consistency (e.g. a `/s` in `units` ⟹ `FLUXSUM`/`TMEAN`, never `SUM`)? Leaning: a compile-time-style
   audit plus a units-sanity assertion in the registry test.
