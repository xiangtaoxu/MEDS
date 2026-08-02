# MEDS v0.1 — Diagnostic I/O: audit, target architecture, and variable plan

**Status: IMPLEMENTED (2026-08-02), branch `feature/io-v01`.** P0–P5 all landed. This document is
kept as the design record; `docs/science/diagnostics.md` is the user-facing page.
**Date:** 2026-08-02. Design decisions D1–D6 (§8) settled the same day; nothing in this plan is open.

> **What shipped, against what was planned.** ~203 registered variables (planned ~229; the shortfall is
> the disturbance-area and mortality-by-pathway rows, and per-band albedo, which need physics-side
> bookkeeping that does not exist yet — see the deferred list below). All 8 groups and all 7 axes are
> live. Verified: ifx Release + Debug 38/38, nvfortran multicore 38/38, and output **byte-identical at
> 1 vs 4 threads** across all 75 files of a 3-year run with `ecophys` + `radiation` on.
>
> **Two findings from the implementation, both recorded in the sections below:**
> 1. The P0 reduction port is **byte-identical at `-O0`** across all 75 output files but differs by
>    1–2 ULP on 6 variables at `-O2`. That is FMA contraction on mathematically equivalent expressions,
>    not a logic change — so "bit-identical" is the right gate only without contraction, and §7.1 now
>    says so.
> 2. The tendency variables (`dbh_growth`, `agb_growth`, `mort_rate`) were first sourced from
>    `site%deriv`. That bundle is deliberately **not** lockstep-reordered, so at the output tick — which
>    runs after the monthly fiss/fuse — index `i` referred to a different plant. It passed at one thread
>    and failed at four, because thread count perturbs which cohorts fuse. Fixed by recording them into
>    `cohort%sdiag`, which does ride the lockstep. **The thread-invariance test is what caught it**, and
>    that is the strongest argument in this whole document for keeping it in the suite.
>
> **Deferred, and honestly out rather than quietly dropped:** mortality carbon by pathway (background /
> cull / disturbance) and the disturbance area flux, which need new bookkeeping in the core engine's
> kill paths; per-band albedo and the up-welling SW/LW terms, which need the two-stream `rad_flux_t`
> surfaced the way the leaf/hydraulic fluxes now are; `AGG_MEANSQ`/variance, still defined and still
> unused; and the P2 asynchronous writer, which §8 D-open already recommended not building until a
> measured fast-stream stall demands it.
**Scope:** the last major component before a MEDS v0.1 release — make the model able to *report itself*.
Concretely: (a) a principled module layer that **computes** diagnostic quantities and **aggregates them
across the demographic scales** (cohort → patch → site, and → PFT / → size class), (b) a much larger,
ED2-informed **variable set** covering fast-scale fluxes/states, slow-scale demography, and
individual-level ecophysiology, and (c) a **TOML surface** that turns any variable on or off at any
temporal tier.

Companion documents (read first, do not duplicate): `MEDS_IO_DESIGN.md` (the aggregation subsystem, P0/P1
shipped), `MEDS_FAST_OUTPUT_DESIGN.md` (the FAST tier, shipped). This plan is the **next layer up**: it
takes the shipped machinery as given and fixes the *source* side.

---

## 0. Executive summary

The output **sink** — registry, temporal integrators, per-tier time-chunked netCDF files, CF metadata,
TOML group/frequency toggles, the FAST tier — is **built and working**. It is not the thing standing
between MEDS and a usable v0.1.

The **source** side is the bottleneck, and it fails in four specific ways:

1. **There is no diagnostics *module*, only a bag of totals.** `meds_output_diagnostics` is ~20
   hand-written `total_*` functions, each of which re-walks the patch/CSR loop with its own weighting
   inlined. Every new site-level number is a new bespoke function. There is no reusable
   cohort→patch→site reduction, no PFT axis, no size-class axis, and no place where a *derived*
   quantity (soil matric potential from θ, VPD from CAS twins, LAI from `nplant·leaf_area`, PLC from
   ψ) is computed once and shared.

2. **The fast loop computes ~40 per-cohort and per-patch quantities every `dt_fast` and throws away
   all but 3 + 9.** `A_net`, `g_sw`, `C_i`, transpiration, `ψ_leaf`, `ψ_wood`, PLC, sapflow, root
   uptake per layer, boundary-layer conductances, in-canopy wind, absorbed PAR/SW/LW per cohort,
   albedo, ground SW/LW, runoff, drainage, infiltration, throughfall, snow state, per-patch CAS twins
   — none of these survive the sub-step. **This is the single blocker for everything the user asked
   for under "individual-level ecophys and biophysical variables."**

3. **The FAST tier is a second, parallel extraction path.** Adding one fast-capable per-cohort
   variable today costs **four** coupled edits (an `SRC_*` code, an `add_variable()` line, a slow-path
   `extract_variable` case, and a fast-path `fast_sample_t` field + staging in `fast_dynamics` +
   `extract_fast_scalar`/`output_integrate_fast` cases). The design doc's honest "two coupled edits"
   price has doubled.

4. **Whole classes of output do not exist at all**: no per-patch physical variables, no PFT-resolved
   output (ED2's most-used analysis axis), no size-class/DBH-class output, no slow-loop demographic
   *rates* (growth, mortality, recruitment, litterfall by flux), no NPP allocation components, and no
   water/energy budget residuals.

The plan below fixes these in six phases (`P0 → {P1, P2} → P3 → P4 → P5`), ending at ~229 registered
variables across 8 groups and 8 axes — cohort, patch, site, PFT, DBH size class, soil layer,
**soil layer × patch**, and radiation band — all TOML-selectable per tier and per variable. §4 is the
deliverable variable list; §4.0 is the short "minimum set to judge whether a run is sane," which is
what the request asked for first. Scope is **single-site only** (§8 D5).

---

## 1. Audit — what exists today

### 1.1 The aggregation subsystem (shipped, sound, keep)

| Module | Lines | Role | Verdict |
|---|---|---|---|
| `src/shared/config/meds_output_config.f90` | 100 | `output_config_t`; `FREQ_*` (4 tiers), `GRP_*` (5), `FC_*`, `SYNC_*` | keep; extend |
| `src/io/meds_output_types.f90` | 207 | `var_desc_t`, `integ_buffer_t`, `pending_record_t`, `stream_file_t`, `output_manager_t`, `fast_sample_t`; `AGG_*` (8), `DIM_*` (5) | keep; extend |
| `src/io/meds_output_registry.f90` | 424 | the registration list (56 vars) + group/tier/override resolution + buffer allocation | keep; restructure the list |
| `src/io/meds_output_integrate.f90` | 519 | `extract_variable` switchboard, `integrate/normalize/reset`, `output_integrate` tick, `output_integrate_fast` | **split** (see §3) |
| `src/io/meds_output_stream.f90` | 347 | per-tier per-chunk netCDF files, CF `cell_methods`/`_FillValue`, F/D/M/Y rollover | keep; add axes |
| `src/io/meds_output_manager.f90` | 54 | `output_serialize_pending` — the only netCDF-touching drain | keep |
| `src/io/meds_output_diagnostics.f90` | 298 | the `total_*` reductions | **replace** (see §3.3) |

CMake targets: `meds_output_{diagnostics,types,integrate,registry}` build into **`meds_io_prep`**
(links `meds_core` only, no netCDF — the DAG-hygiene wall); `meds_io` + `meds_output_{stream,manager}`
build into **`meds_io_stream`** (adds `meds_netcdf_c`). The two new modules of §3.2/§3.3 join
`meds_io_prep`, so **every diagnostic module stays inside `src/io/`**.

The **one** deliberate exception is the existing `src/shared/config/meds_output_config.f90`, which
holds `output_config_t` and the `FREQ_*`/`GRP_*`/`FC_*` codes. It cannot move into `src/io/`: it is
what lets `meds_config` — the DAG **root** — carry the `[output]` block with no `shared → io`
back-edge (the same anti-back-edge pattern `meds_forcing_config` uses). This plan's config additions
(§5.1: the `axes_*` toggles, the three new `GRP_*` codes, `dbh_class_edges`) extend that file in place;
nothing else is added under `src/shared/`.

The load-bearing invariants this machinery already guarantees, which the plan must not break:

- **No pointers into the SoA.** Values are *copied* out of live state each step; a cohort birth/death
  is a non-event for output.
- **Cohort/patch slot sets are fixed within an output window** (restructuring only at month/year
  boundaries), so slabs fold by direct slot index. `close_tier` fires **before** the new step is
  folded, so a month-boundary step closes the old window on the old slot set and opens the new window
  on the new one. *Verified consistent — do not "fix" this.*
- **The netCDF wall:** the stepper edge (`meds_output_{config,types,registry,integrate}`) is
  netCDF-free; only `meds_output_{stream,manager}`, referenced from `meds_main` alone, touch C.
- **Serial-order reduction** in the threaded fast loop, so output is byte-identical at any thread count.

### 1.2 Current registry inventory — 56 variables

| Group | Axis | Count | Variables |
|---|---|---|---|
| STRUCTURE | cohort | 10 | `nplant_cohort` `dbh_cohort` `height_cohort` `basal_area_cohort` `agb_cohort` `leaf_area_cohort` `growth_avg_cohort` `pft_cohort` `owner_patch` `global_cohort_id` |
| STRUCTURE | patch | 6 | `area_patch` `age_patch` `dist_type_patch` `cohort_offset` `cohort_count` `global_patch_id` |
| STRUCTURE | site | 4 | `nplant_site` `basal_area_site` `agb_site` `lai_site` |
| CARBON | site | 10 | `gpp_site` `npp_site` `rh_site` + 7 `soilc_*_site` |
| ENERGY | site | 2 | `cas_temp_site` `soil_temp_top_site` |
| ENERGY/WATER | soil | 2 | `soil_temp_site` `soil_water_site` |
| WATER | site | 1 | `et_site` |
| NUMERICS | site | 11 | `work_integ_steps/rej`, `work_soil_nsub`, `work_hydro_nsub`, `work_nonconv`, `work_hydro_thrash`, `work_rk45_rescue`, `work_clamp_{stage,commit,mass,energy}` |
| CARBON/ENERGY/WATER | site (FAST-only) | 7 | `gpp_rate_fast` `le_flux_fast` `h_flux_fast` `rnet_fast` `sw_in_fast` `ustar_fast` `air_temp_fast` |
| CARBON/ENERGY | cohort (FAST-only) | 3 | `leaf_temp_cohort_fast` `gpp_cohort_fast` `height_cohort_fast` |

`MAX_OUTPUT_VARS = 128`. Dead code found: `SRC_S_N_COHORT` / `SRC_S_N_PATCH` have source codes and
`extract_scalar_source` cases but **no `add_variable()` line** — the `n_cohort`/`n_patch` variables the
design table promises are never registered.

### 1.3 The legacy `[io]` snapshot stream (to retire)

`meds_io.f90` still carries `io_create` / `io_write_snapshot` / `io_close` writing
`<prefix>-D-output.nc`: an annual-cadence, instantaneous, hard-coded 21-variable schema, driven from a
second `is_new_year` block in `meds_main`. It duplicates what `[output]` does strictly better and it
**collides on the `-D-` filename prefix** with the aggregation subsystem's daily stream. It exists only
for back-compat.

`io_write_state` / `io_read_state` (the `-S-<stamp>.nc` restart stream) are **orthogonal and stay** —
a checkpoint must be raw prognostic state at an instant, never a time-average.

---

## 2. Findings — the gaps, ranked

**F1 — the sink is done; the source is the bottleneck.** No further work is needed on temporal
aggregation, file chunking, CF metadata, or the TOML resolution order. All remaining effort belongs on
*producing* the numbers.

**F2 — `meds_output_diagnostics` is not a module, it is a bag of totals.** Twenty near-identical
functions, each re-deriving the patch loop and the CSR slice, each hard-wiring its own weight
(`area`, `area·nplant`, `area·nplant·1e-4`). Consequences: (a) adding a site variable is ~10 lines of
boilerplate; (b) there is no patch-level output at all, because the patch-level intermediate is never
materialized; (c) there is no PFT or size-class axis, because those need a *grouped* reduction that
does not exist; (d) derived quantities (soil ψ from θ, VPD from CAS enthalpy+shv, LAI from
`nplant·leaf_area`, PLC from ψ) have nowhere to live and are computed nowhere.

**F3 — the FAST tier is a parallel hand-wired extraction path.** `output_integrate` loops
`do t = 2, N_FREQ`, hard-skipping FAST; the FAST tier is fed instead by `output_integrate_fast` reading
`fast_sample_t` (9 hand-coded scalar fields) plus three 2-D cohort/soil staging arrays filled inline in
`fast_dynamics`. A registry variable that carries `FREQ_FAST` but is not one of those 12 hand-staged
sources silently produces nothing. **Adding a fast-capable per-cohort variable is four edit sites.**

**F4 — the fast loop discards nearly everything it computes.** Per `dt_fast`, per patch, MEDS
currently computes and then drops:

- *per cohort:* `leaf_flux_t{A_net, A_gross, gs, ci, cs, transpiration, rd, limitation}`;
  `hydro_flux_t{sapflow, root_uptake, root_uptake_layer(:), psi_leaf, psi_wood, plc, nsub}`;
  `forc%{abs_sw, abs_par, abs_lw, abs_sw_wood, abs_lw_wood}`; `aero%{wind, leaf_gbh, leaf_gbw,
  wood_gbh, wood_gbw}`; `bio%{wood_temp, leaf_water_mass, wood_water_mass, leaf_surf_water,
  wood_surf_water}`; the β\_stomata / β\_nonstomata stress factors; intercepted water.
- *per patch:* `rad_flux_t{albedo(band), dn_ground(band), up_ground(band)}`; `frozen%{infiltration,
  drainage, runoff_surf, precip_ground, uptake}`; `aero%{tstar, qstar, cstar, zeta, obu, ggnet, rough,
  displace, uh}`; snow column; ponding; the ground surface energy balance terms;
  `budg%{whole_energy, whole_water, cas_co2}` residuals; `budg%nee_last`.

Only `leaf_temp`, per-cohort `gpp`, `height` (cohort) and 9 site scalars are staged. **Everything the
request asks for under "individual-level ecophys and biophysical variables (Anet, gsw, gsc)" lives in
this discarded set.**

**F5 — no patch-level physical output; `DIM_PFT` is declared but never implemented.** The `DIM_PATCH`
axis carries only geometry/identity (`area`, `age`, `dist_type`, CSR, id). `DIM_PFT` exists in
`meds_output_types` and in the `select case` of `meds_output_stream`, but no variable uses it and
`stream_open_file` never defines a `pft` dimension.

**F6 — no PFT-resolved and no size-class-resolved output.** ED2's most-used analysis products are
`agb(dbh-class, pft)`, `nplant(dbh-class, pft)`, `basal_area(dbh-class, pft)` — the forest-inventory
view. MEDS can only emit the raw cohort slab and make the reader do the binning, which is fragile
(cohort identity churns) and puts model knowledge (PFT count, DBH class edges) in post-processing.

**F7 — slow-loop rates are never surfaced.** `cohort_deriv_block` (`d_dbh_dt`, `d_agb_dt`,
`d_wood_carbon_dt`, `dln_nplant_dt`, …) is computed every step and is still resident on `site%deriv` at
the output tick, but nothing reads it. NPP allocation components (`g_leaf`, `g_fineroot`, `g_wood`,
`npp_store`, `g_repro`, `growth_resp`) are locals inside `compute_carbon_allocation`. `litter_input_t`
is built per patch and consumed by biogeochemistry without ever being reported. Recruitment counts,
mortality-by-cause, and disturbance area fluxes are likewise invisible.

**F8 — the config surface is too coarse for the target variable count.** Five groups (`structure`,
`carbon_fluxes`, `water_fluxes`, `energy_fluxes`, plus `numerics`) over ~180 variables is a blunt
instrument; the per-variable `meds_io_config.toml` path is implemented but has **no shipped example
file anywhere in the repo**, and the `[output]` block in `meds_config_main.toml` is entirely commented
out — so the feature is effectively undiscoverable.

**F9 — hygiene.** `MAX_OUTPUT_VARS = 128` is below the target count. `n_cohort`/`n_patch` are dead.
The `_fast` name suffix on `gpp_rate_fast` etc. conflicts with the design's own scale-suffix convention
(`_site`/`_patch`/`_cohort`) — a *tier* is not a *scale*, and the same quantity should not change name
because it is reported hourly. Integrator slab buffers are sized to `cohort_max` (4096) rather than the
live cap, which at ~60 cohort-dimensioned variables × 4 tiers × 4 arrays is ~30 MB of mostly-fill.

**F10 — legacy duplication.** Two diagnostic writers, two `[io]`/`[output]` config blocks, colliding
`-D-` filenames.

---

## 3. Target architecture

### 3.1 The five-stage diagnostic wall

Today there are effectively two stages (raw state → temporal integrate). The target inserts three:

```
  raw prognostic state (site_t, patch reservoirs, cohort SoA)
      │
  [1] DERIVE          meds_diagnostic_kernels     pure per-entity derived quantities
      │                                            (soil psi(theta), VPD(cas), LAI(nplant,leaf_area),
      │                                             PLC(psi), wetness, albedo, C balance closure)
      │
  [2] CAPTURE         cohort_diag_t / patch_diag_t  dt-weighted sub-step accumulators for
      │                                             everything the fast loop computes and drops
      │
  [3] REDUCE          meds_diagnostic_reduce       generic weighted scale aggregation
      │                                            cohort -> {patch, site, pft, size-class}
      │                                            patch  -> site
      │
  [4] INTEGRATE       meds_output_integrate        (UNCHANGED) temporal folding per (var, tier)
      │
  [5] SERIALIZE       meds_output_stream/manager   (UNCHANGED) per-tier per-chunk netCDF
```

Stages 4–5 are shipped. Stages 1–3 are the work.

**Every new module in this plan lives in `src/io/`.** Nothing diagnostic goes into
`src/shared/functions/`. The rule and its consequences, stated once so the phases can just follow it:

- **Both new modules join the `meds_io_prep` CMake target** (`src/io/meds_output_diagnostics.f90`,
  `..._types`, `..._integrate`, `..._registry`), which links `meds_core` **only** — no netCDF. So they
  sit inside the existing DAG-hygiene wall: the stepper edge (`meds_aux → meds_io_prep`) still pulls
  no C dependency, and `meds_output_diagnostics` — the module §3.3 replaces — is already there, so the
  replacement is an in-place swap rather than a new library edge.
- **`meds_io_prep` is an explicit file list, not a `file(GLOB)`** ([CMakeLists.txt:131-135](../../CMakeLists.txt#L131-L135)).
  Adding a source there is a required manual edit; it will not auto-add the way `src/plant` and
  `src/biophysics` do.
- **Accepted consequence:** diagnostic kernels in `src/io/` cannot be called from anything *below* io
  in the DAG — `plant`, `biophysics`, `biogeochemistry`, `core`, `shared` — without a back-edge.
  That is correct rather than merely tolerable: a *diagnostic* is by definition something no physics
  kernel needs, and putting it below the wall would be the first step toward physics depending on its
  own reporting layer. If a quantity ever turns out to be needed by physics, that is the signal it was
  never a diagnostic and belongs in the owning physics library instead.
- Tests link `meds_io_prep` already, so the kernels stay directly unit-testable. The one thing this
  placement costs is reuse from the `libmeds_plant_c` Python path, which links `meds_shared` only —
  acceptable, since that API exposes leaf gas exchange and phenology, not site diagnostics.

### 3.2 New: `src/io/meds_diagnostic_kernels.f90`

`pure`/`elemental` functions taking plain arguments, returning derived quantities. No `site_t`, no
config aggregator — so they stay directly unit-testable. Examples:

```fortran
pure elemental real(wp) function cohort_lai(nplant, leaf_area)          ! [m2/m2]
pure elemental real(wp) function cohort_gsc(gsw)                        ! gsw / 1.6  [mol CO2/m2/s]
pure elemental real(wp) function cohort_wue(a_net, transpiration)       ! [umol CO2 / mmol H2O]
pure elemental real(wp) function soil_wetness(theta, theta_res, theta_sat)
pure elemental real(wp) function cas_vpd(can_temp, can_shv, can_prss)   ! [Pa]
pure elemental real(wp) function plc_from_psi(psi, p50, a_vuln)         ! [-]
```

Content rule: anything already living in a physics library (`soil_psi_from_theta` in `meds_hydr_lib`,
`uext_to_temp` in `meds_therm_lib`, `psi_from_water_content`) is **called**, not re-implemented — those
live in `meds_shared`, which `meds_io_prep` already links transitively, so calling down is free. This
module holds only diagnostics that no physics kernel owns.

### 3.3 New: `src/io/meds_diagnostic_reduce.f90` — replaces `meds_output_diagnostics`

One generic reduction family, parameterized by an explicit **weight kind**, replacing 20 bespoke
`total_*` functions:

```fortran
integer(ik), parameter :: W_NONE = 0, W_NPLANT = 1, W_LEAF_AREA = 2, W_BASAL_AREA = 3,   &
                          W_AGB = 4, W_AREA = 5, W_NPLANT_AREA = 6

! cohort field -> per-patch value (weighted mean if `mean`, weighted sum otherwise)
pure subroutine reduce_cohort_to_patch(site, field, wkind, mean, out, n_out)
! cohort field -> site scalar (patch-area-weighted roll-up of the above)
pure real(wp) function reduce_cohort_to_site(site, field, wkind, mean)
! cohort field -> per-PFT slab                              [DIM_PFT]
pure subroutine reduce_cohort_to_pft(site, field, wkind, mean, n_pft, out)
! cohort field -> per-DBH-class slab                        [DIM_SIZE]
pure subroutine reduce_cohort_to_size(site, field, wkind, mean, edges, out)
! patch field  -> site scalar (area-weighted)
pure real(wp) function reduce_patch_to_site(site, field, mean)
! per-layer patch column -> site column                     [DIM_SOIL / DIM_SNOW]
pure subroutine reduce_patch_column_to_site(site, col, n, out)
! per-layer patch column -> the FLATTENED (layer, patch) slab, no reduction  [DIM_SOIL_PATCH]
pure subroutine gather_patch_columns(site, col, nlayer, out, n_out)
```

Payoff, stated plainly: **every per-cohort variable gets its patch, site, PFT and size-class twins for
free**, from one registry declaration. That is what turns 56 variables into ~180 without 180 hand-written
reductions, and it is the direct answer to "some of them are extensive and can be aggregated to
patch/site levels."

Two rules the reducer must enforce, because getting them wrong is the classic diagnostics bug:

- **Extensive vs intensive.** An extensive per-plant quantity (`agb` [kgC/plant], `leaf_area`
  [m²/plant]) aggregates as `Σ area·nplant·x` and lands in per-ground-area units. An intensive
  quantity (`dbh`, `leaf_temp`, `psi_leaf`, `gsw`) aggregates as a **weighted mean**, and the weight
  must be declared (basal area for `dbh`, leaf area for `leaf_temp`/`gsw`/`A_net`, `nplant` for
  demographic rates). The registry descriptor gains a `weight` field so this is data, not code.
- **Empty-set guard.** A patch with no cohorts, or a PFT/size bin with no members, emits
  `_FillValue`, never `0/0`. The existing `normalize_*` `valid` flag already carries this downstream.

### 3.4 New state: the fast diagnostic accumulator blocks

This is the largest and most invasive change, and the one that unlocks F4. Model it on ED2's
`fmean_*` fields on `patchtype`, but keep MEDS's discipline:

```fortran
! src/core/meds_core_state_types.f90 (or a sibling meds_core_diag_types)
type :: cohort_diag_block          ! SoA, one entry per live cohort slot
   integer(ik) :: n = 0, cap = 0
   real(wp) :: wsum = 0.0_wp                  ! Sum(dt) — the dt weight for every field below
   real(wp), allocatable :: a_net(:), a_gross(:), gsw(:), ci(:), cs(:), rd(:)
   real(wp), allocatable :: transp(:), sapflow(:), root_uptake(:)
   real(wp), allocatable :: psi_leaf(:), psi_wood(:), plc(:)
   real(wp), allocatable :: leaf_temp(:), wood_temp(:), leaf_vpd(:)
   real(wp), allocatable :: abs_par(:), abs_sw(:), abs_lw(:), light_level(:)
   real(wp), allocatable :: gbw(:), wind(:)
   real(wp), allocatable :: beta_stomata(:), beta_nonstomata(:)
   real(wp), allocatable :: intercept(:), leaf_surf_water(:)
end type

type :: patch_diag_block           ! SoA, one entry per live patch slot
   integer(ik) :: n = 0, cap = 0
   real(wp) :: wsum = 0.0_wp
   real(wp), allocatable :: le(:), h(:), rnet(:), ground_heat(:)
   real(wp), allocatable :: sw_up(:), lw_up(:), albedo(:), par_ground(:), sw_ground(:)
   real(wp), allocatable :: transp(:), evap_ground(:), evap_canopy(:)
   real(wp), allocatable :: infiltration(:), drainage(:), runoff(:), throughfall(:)
   real(wp), allocatable :: snowmelt(:), swe(:), snow_depth(:), snow_cover(:)
   real(wp), allocatable :: nee(:), gpp(:), reco(:)
   real(wp), allocatable :: cas_shv(:), cas_co2(:), cas_vpd(:), ground_temp(:)
   real(wp), allocatable :: resid_energy(:), resid_water(:), resid_co2(:)
   real(wp), allocatable :: soil_psi(:,:), soil_fliq(:,:), root_uptake_layer(:,:)   ! (layer, patch)
end type
```

**Six design decisions, with the reasoning, because this block is where the plan can go wrong:**

1. **dt-weighted sums, normalized on read.** Every field accumulates `x·dt_fast`; `wsum` accumulates
   `dt_fast`. `extract_variable` divides. This makes the block correct under the adaptive integrator
   (unequal sub-steps) and under `AGG_TMEAN` chaining, for free.
2. **Reset once per slow step**, exactly like `gpp_accum` / `et_accum` / `xi_accum` — the established
   fast→slow bridge lifecycle. Do not invent a second reset cadence.
3. **TRANSIENT — never written to the restart file.** These are diagnostics; a checkpoint restores
   prognostic state only. This is what keeps the change out of `io_write_state`/`io_read_state`.
4. **They must still ride the cohort/patch lockstep reorder.** The reset happens once per slow step,
   but restructuring (fuse/split/cull/recruit/disturb) happens *inside* the slow step, after the fast
   loop has filled the block and possibly before the monthly window closes. So `cohort_reorder` /
   `copy_cohort_slot` / `cohort_compact` / `cohort_ensure_capacity` and the patch `sort_patches` /
   `patch_compact` sites **must** be extended. This is the real cost of the change and is exactly the
   "add a per-cohort field → update *these*" invariant in CLAUDE.md. Fusion rule: the same
   extensive/intensive split as §3.3 (leaf-area-weight the intensive fields, `nplant`-sum the
   extensive ones) so a fused cohort's diagnostics are the physically correct blend.
5. **Allocated only when at least one variable in the block is registered live.** A run with
   `[output].carbon_fluxes = false` and no cohort ecophys variables pays nothing. Gate on a
   `mgr`-derived `need_cohort_diag` / `need_patch_diag` flag computed once from the registry.
6. **Filled in the same serial-order fold as the existing `red_fast` staging**, so threading stays
   byte-identical.

**This block also collapses F3.** Once per-cohort and per-patch diagnostics live in site state, the
FAST tier no longer needs its own extraction path: `output_integrate` can fold *any* registry variable
at *any* tier through the *one* `extract_variable` switchboard. The FAST tier's remaining special case
is only *when* it closes (every `fast_interval_steps × dt_fast`), which is cadence, not extraction.

**Alternative considered and rejected:** keep the staging arrays on `output_manager_t` and grow them.
That avoids touching the lockstep, but (a) it puts model state in the IO manager, (b) it cannot survive
a mid-slow-step restructuring, and (c) it entrenches the two-path extraction (F3). The state-hub
placement is more invasive once and cheaper forever.

### 3.5 New: slow-loop rate capture

Small and cheap, unlike §3.4. Two additions:

- **Read `site%deriv` directly.** It is live at the output tick and already carries `d_dbh_dt`,
  `d_agb_dt`, `d_wood_carbon_dt`, `d_leaf_carbon_dt`, `dln_nplant_dt`. Add `SRC_D_*` cases in
  `extract_variable`. Zero new state. Covers growth and mortality rates at cohort resolution.
- **A per-patch `slow_diag_t`** for the quantities that are currently locals: NPP allocation
  components (`npp_leaf`, `npp_fineroot`, `npp_wood`, `npp_storage`, `npp_repro`, `growth_resp`,
  `storage_deficit`), the `litter_input_t` fluxes, recruitment count/carbon, and the mortality carbon
  flux by pathway (background / cull / disturbance). Filled once per slow step by
  `meds_vegetation_dynamics` and `meds_core`, read at the tick, reset after. Same lockstep obligations
  as §3.4 but only on the patch axis (6 permute sites, not the cohort machinery).

### 3.6 New axes

| Code | Axis | Length | Serializer work |
|---|---|---|---|
| `DIM_PFT` | plant functional type | `cfg%pft%n_pft` (run-time; §8 D2) | define the `pft` dim (**declared but never built today**) |
| `DIM_SIZE` | DBH size class | `n_dbh_class` from TOML edges | new dim + a `dbh_class_edges` coordinate variable |
| `DIM_SNOW` | snow layer | `n_snow_layer_max` | new dim, mirrors `DIM_SOIL` |
| `DIM_BAND` | radiation band | `n_band` (VIS/NIR/TIR) | new dim; only for albedo/RT diagnostics |
| `DIM_SOIL_PATCH` | soil layer × patch | `n_soil_layer_max × patch_max` | **included** (§3.6.1, its own phase P1) |

`DIM_PFT` and `DIM_SIZE` are the highest-value additions in the whole plan per line of code, because
they turn the cohort slab (which a reader must interpret through churning ids) into stable,
directly-plottable model output.

#### 3.6.1 `DIM_SOIL_PATCH` — the 2-D (layer, patch) axis

Per-patch soil columns diverge strongly — a young disturbance gap and a closed canopy dry down at
different rates and to different depths — and the area-weighted site column hides exactly that. This
axis emits the per-patch column directly.

**The implementation is a flattening, not a second buffer type.** `integ_buffer_t%slab(:)` stays 1-D;
a `DIM_SOIL_PATCH` variable folds into it at

```
   index = (ip - 1) * n_soil_layer_max + k        ! k = 1..n_soil_layer_max, ip = 1..n_patch
```

striding by the **compile-time constant** `n_soil_layer_max = 20`, never by the live patch or layer
count, so the mapping is stable regardless of how many layers are active. `n_slab = n_patch ·
n_soil_layer_max`. Consequences, all favourable:

- **`integrate_slab` / `normalize_slab` / `reset_buffer` are untouched.** They already fold an
  arbitrary-length 1-D slab by direct index; the fixed-slot-set argument of §4.4 of `MEDS_IO_DESIGN.md`
  carries over verbatim because the patch slot set is fixed within a window for the same reason the
  cohort one is.
- **No new C binding.** `nc_put_vara_double` takes assumed-size `startp(*)`/`countp(*)`
  ([meds_netcdf_c.f90:108](../../src/io/meds_netcdf_c.f90#L108)), so a rank-3 write is
  `start = [nrec, 0, 0]`, `count = [1, n_patch, n_soil_layer_max]` with no interface change.
- **Serializer work is confined to two places:** `stream_open_file` defines the variable with dims
  `(time, patch, soil)` and chunk `(1, patch_max, n_soil_layer_max)`; `stream_write_record` selects
  the 3-D `start`/`count` in the existing `select case (v%dim)`.
- **Cost:** `max_slab` grows to `max(cohort_max, patch_max, patch_max · n_soil_layer_max)` — at
  `patch_max = 256` that is 5120 doubles, comparable to `cohort_max = 4096`, so the sizing story is
  unchanged. The `_FillValue` tail for inactive layers/patches deflates to near-nothing.

Gated by its own `[output].axes_soil_patch` toggle (default **off**) because it is the highest-volume
non-cohort axis.

### 3.7 Registry ergonomics

With ~180 variables the current one-`call add_variable(...)`-per-line block becomes unreadable. Two
changes:

- **Group the registration list into `contains`-local `register_<group>_<axis>` subroutines**
  (`register_carbon_cohort`, `register_water_site`, …), one per §4 sub-table, so the source order
  matches the documentation order and a reviewer can diff one domain.
- **Auto-generate scale twins.** A helper `add_variable_family(base, ..., scales = 'CPSF')` emits the
  cohort/patch/site/PFT variants of one quantity from one line, with the correct weight kind and the
  correct default streams per scale. This is what keeps the "two coupled edits" price honest at scale:
  one family line + one `extract_variable` case for the *cohort-level* source, and the reducer supplies
  every aggregate.
- Raise `MAX_OUTPUT_VARS` to 512; size integrator slabs to the live cap rather than `cohort_max`.

---

## 4. The variable plan

### 4.0 The minimum set — "is this run sane?"

Before the exhaustive list: **these 32 variables are what you actually look at first.** They are the
answer to "key variables we need to check the simulated fast-scale and slow-scale dynamics," and they
should be the shipped default `[output]` configuration.

**Fast-scale (diurnal, FAST + DAILY) — 14**

| Variable | Why it is the first thing to check |
|---|---|
| `sw_in_site`, `rnet_site` | is the radiation forcing/absorption physically shaped? |
| `le_site`, `h_site` | the Bowen ratio — the single most diagnostic energy number |
| `gpp_site` | does the light response have the right shape and magnitude? |
| `nee_site` | night = +Reco, day = −uptake; the sign flip is the integration test |
| `et_site` | closes against `le_site`; catches a unit/latent-heat error instantly |
| `cas_temp_site`, `cas_vpd_site` | canopy-air coupling; the first thing to oscillate if `dt_fast` is too big |
| `cas_co2_site` | drawdown amplitude; a stuck value means the CO₂ twin is not coupled |
| `soil_temp_site` (layer) | the diurnal thermal wave, and the annual wave's damping depth |
| `soil_water_site` (layer) | drydown/rewetting shape; catches a broken Richards BC |
| `ustar_site` | is the surface layer turbulent at all? (a floored `ustar` fakes decoupling) |
| `resid_energy_site`, `resid_water_site` | budget closure — **a nonzero value invalidates everything above** |

**Slow-scale (DAILY/MONTHLY/ANNUAL) — 18**

| Variable | Why |
|---|---|
| `agb_site`, `lai_site`, `basal_area_site`, `nplant_site` | the four stand-structure scalars |
| `agb_pft`, `lai_pft` | is PFT composition doing anything, or is one PFT winning trivially? |
| `agb_size`, `nplant_size` | the size distribution — the demographic core's actual output |
| `npp_site`, `rh_site`, `nep_site` | the carbon balance; `nep` drift is the long-run sanity check |
| `ra_site` | autotrophic fraction — should be ~0.5 of GPP, a strong parameter check |
| `agb_growth_site`, `agb_mort_site`, `agb_recruit_site` | the demographic rates that *make* the AGB trajectory |
| `soilc_total_site` | is soil carbon spinning up, equilibrating, or running away? |
| `n_cohort`, `n_patch` | fuse/fission health — a monotonically climbing count is a config bug |
| `canopy_height_site` | stand development, and the CAS-depth driver |

---

### 4.1 Conventions

- **Scale suffix, always:** `_cohort` / `_patch` / `_site` / `_pft` / `_size`. Identity/CSR/count
  fields are exempt (`global_cohort_id`, `cohort_offset`, `n_cohort`).
- **No tier suffix.** Retire `_fast` from `gpp_rate_fast` etc. (F9); the same quantity keeps one name
  at every tier, and `cell_methods` + the file's stream letter say how it was reduced.
- **Rates keep `/s` and use `AGG_TMEAN`** (a dt-weighted mean rate) unless the variable is explicitly a
  period total, in which case `AGG_FLUXSUM` and the units drop the `/s`. Plain `AGG_SUM` stays reserved
  for integer count tallies. *Current `gpp_site`/`npp_site`/`et_site`/`rh_site` use `AGG_SUM` on a
  slow-step accumulator — these need re-classifying to `FLUXSUM` when the fast diag block lands, since
  the accumulator already contains the dt weighting.*
- **Streams column:** `F` fast, `D` daily, `M` monthly, `Y` annual. Cohort/patch-dimensioned variables
  may not carry `Y` (registry-enforced).
- **Status column:** `READY` = registered or one `extract_variable` case away; `PLUMB` = the number is
  computed today but discarded, needs §3.4/§3.5 capture; `DERIVE` = needs a §3.2 kernel; `NEW` = needs
  new physics-side bookkeeping.

---

### 4.2 Group taxonomy (revised from 5 to 8)

| Group | Covers | Default |
|---|---|---|
| `GRP_STRUCTURE` | demography, size, identity, patch geometry | on |
| `GRP_CARBON` | GPP/NPP/Ra/Rh/NEE, allocation, litter, carbon stocks | on |
| `GRP_WATER` | ET, transpiration, evaporation, runoff/drainage, soil & plant water | on |
| `GRP_ENERGY` | LE/H/Rnet/G, temperatures, CAS state | on |
| `GRP_RADIATION` | albedo, absorbed SW/PAR/LW, light levels | off |
| `GRP_ECOPHYS` | **new** — individual-level leaf gas exchange + hydraulics (cohort axis) | off |
| `GRP_BIOGEOCHEM` | soil carbon pools, decomposition scalars, (future) N | on |
| `GRP_NUMERICS` | integrator work counters, budget residuals, clamp telemetry | off |

`GRP_ECOPHYS` is split out of `GRP_CARBON` deliberately: it is the highest-volume group (cohort axis ×
20 variables) and the one a production run will most often want off.

---

### 4.3 Carbon fluxes

| name | long_name | units | dim | agg | streams | source | ED2 analogue | status |
|---|---|---|---|---|---|---|---|---|
| `gpp_site` | gross primary productivity | kgC/m²/s | scalar | tmean | F D M Y | `cohort%gpp_accum` → reduce | `fmean_gpp` | READY (reclassify) |
| `npp_site` | net primary productivity | kgC/m²/s | scalar | tmean | F D M Y | GPP − Ra | `fmean_npp` | READY |
| `ra_site` | autotrophic respiration | kgC/m²/s | scalar | tmean | F D M Y | leaf+stem+root maint + growth resp | `fmean_plresp` | PLUMB (growth resp) |
| `rh_site` | heterotrophic respiration | kgC/m²/s | scalar | tmean | F D M Y | `xi_accum%rh_fast_accum` | `fmean_rh` | READY |
| `reco_site` | ecosystem respiration | kgC/m²/s | scalar | tmean | F D M Y | Ra + Rh | — | DERIVE |
| `nee_site` | net ecosystem exchange | kgC/m²/s | scalar | tmean | F D M Y | `budg%nee_last` | `fmean_nep` (sign) | PLUMB |
| `nep_site` | net ecosystem production | kgC/m²/s | scalar | tmean | D M Y | GPP − Ra − Rh | `fmean_nep` | DERIVE |
| `leaf_resp_site` | leaf dark respiration | kgC/m²/s | scalar | tmean | D M Y | `cohort%leaf_resp_accum` | `fmean_leaf_resp` | READY |
| `stem_resp_site` | stem maintenance respiration | kgC/m²/s | scalar | tmean | D M Y | `cohort%stem_resp_accum` | `fmean_stem_resp` | READY |
| `root_resp_site` | fine-root maintenance respiration | kgC/m²/s | scalar | tmean | D M Y | `cohort%root_resp_accum` | `fmean_root_resp` | READY |
| `growth_resp_site` | growth respiration | kgC/m²/s | scalar | tmean | D M Y | `compute_carbon_allocation` local | `fmean_*_growth_resp` | PLUMB (§3.5) |
| `npp_leaf_site` | NPP allocated to leaf | kgC/m²/s | scalar | tmean | M Y | allocation local | `mmean_nppleaf` | PLUMB |
| `npp_fineroot_site` | NPP allocated to fine root | kgC/m²/s | scalar | tmean | M Y | allocation local | `mmean_nppfroot` | PLUMB |
| `npp_wood_site` | NPP allocated to wood | kgC/m²/s | scalar | tmean | M Y | allocation local | `mmean_nppwood` | PLUMB |
| `npp_storage_site` | NPP to non-structural storage | kgC/m²/s | scalar | tmean | M Y | allocation local | — | PLUMB |
| `npp_repro_site` | NPP to reproduction | kgC/m²/s | scalar | tmean | M Y | allocation local | `mmean_nppseeds` | PLUMB |
| `litter_leaf_site` | leaf litterfall carbon | kgC/m²/s | scalar | tmean | M Y | `litter_input_t` | `mmean_fgc_in` | PLUMB |
| `litter_fineroot_site` | fine-root litter carbon | kgC/m²/s | scalar | tmean | M Y | `litter_input_t` | `mmean_fsc_in` | PLUMB |
| `litter_wood_site` | coarse woody debris input | kgC/m²/s | scalar | tmean | M Y | `litter_input_t` | `mmean_stgc_in` | PLUMB |
| `gpp_patch`, `npp_patch`, `nee_patch`, `reco_patch` | patch twins | kgC/m²/s | patch | tmean | D M | reduce | — | PLUMB |
| `gpp_pft`, `npp_pft` | PFT twins | kgC/m²/s | pft | tmean | M Y | reduce | `mmean_gpp_py` | NEW |
| `gpp_cohort`, `npp_cohort` | per-plant carbon flux | kgC/plant/s | cohort | tmean | F D M | `gpp_accum`, resp accums | `fmean_gpp` (co) | READY/PLUMB |

**Carbon stocks** (`AGG_TMEAN`, `GRP_CARBON`):

| name | units | dim | streams | source | status |
|---|---|---|---|---|---|
| `agb_site` `agb_patch` `agb_pft` `agb_size` | kgC/m² | scalar/patch/pft/size | D M Y | `cohort%agb` → reduce | READY / NEW |
| `bgb_site` | kgC/m² | scalar | D M Y | `wood_carbon·(1−aboveground_frac)` + fineroot | DERIVE |
| `leaf_carbon_site` `_patch` `_pft` | kgC/m² | — | D M Y | `cohort%leaf_carbon` | READY/NEW |
| `fineroot_carbon_site` | kgC/m² | scalar | D M Y | `cohort%fineroot_carbon` | READY |
| `wood_carbon_site` | kgC/m² | scalar | D M Y | `cohort%wood_carbon` | READY |
| `storage_carbon_site` | kgC/m² | scalar | D M Y | `cohort%nonstructural_carbon` | READY |
| `veg_carbon_site` | kgC/m² | scalar | D M Y | sum of the four pools | DERIVE |
| `soilc_total_site` | kgC/m² | scalar | D M Y | sum of the 7 CENTURY pools | DERIVE |
| `soilc_*_site` (7 pools) | kgC/m² | scalar | D M Y | `patch%soil_carbon` | READY |
| `soilc_*_patch` (7 pools) | kgC/m² | patch | M | reduce | NEW |
| `carbon_resid_site` | kgC/m² | scalar | D M Y | carbon budget closure | NEW |

---

### 4.4 Water fluxes and states

| name | long_name | units | dim | agg | streams | source | ED2 | status |
|---|---|---|---|---|---|---|---|---|
| `et_site` | evapotranspiration | kg/m²/s | scalar | tmean | F D M Y | `site%et_accum` | — | READY (reclassify) |
| `transp_site` | canopy transpiration | kg/m²/s | scalar | tmean | F D M Y | `leaf_flux%transpiration` | `fmean_transp` | PLUMB |
| `evap_ground_site` | ground/soil evaporation | kg/m²/s | scalar | tmean | F D M Y | `ground_evaporation` | `fmean_vapor_gc` | PLUMB |
| `evap_canopy_site` | canopy-interception evaporation | kg/m²/s | scalar | tmean | F D M | surface-film evaporation | `fmean_vapor_lc/wc` | PLUMB |
| `root_uptake_site` | root water uptake | kg/m²/s | scalar | tmean | F D M | `frozen%uptake` | `fmean_water_supply` | PLUMB |
| `precip_site` | precipitation | kg/m²/s | scalar | tmean | F D M Y | `forc%precip + forc%snowf` | `fmean_pcpg` | PLUMB |
| `snowfall_site` | snowfall | kg/m²/s | scalar | tmean | D M Y | `forc%snowf` | — | PLUMB |
| `throughfall_site` | canopy throughfall | kg/m²/s | scalar | tmean | D M | interception residual | `fmean_throughfall` | PLUMB |
| `infiltration_site` | soil infiltration | kg/m²/s | scalar | tmean | D M | `frozen%infiltration` | — | PLUMB |
| `runoff_site` | surface runoff | kg/m²/s | scalar | tmean | D M Y | `frozen%runoff_surf` | `fmean_runoff` | PLUMB |
| `drainage_site` | bottom drainage | kg/m²/s | scalar | tmean | D M Y | `frozen%drainage` | `fmean_drainage` | PLUMB |
| `snowmelt_site` | snowmelt water flux | kg/m²/s | scalar | tmean | D M | snow kernels | — | PLUMB |
| `resid_water_site` | column water-budget residual | kg/m² | scalar | fluxsum | F D M Y | `budg%whole_water` | `water_residual` | PLUMB |

**Water states:**

| name | units | dim | agg | streams | source | ED2 | status |
|---|---|---|---|---|---|---|---|
| `soil_water_site` | m³/m³ | soil | tmean | F D M Y | `soil_w%theta` | `fmean_soil_water` | READY |
| `soil_psi_site` | MPa | soil | tmean | F D M Y | `soil_psi_from_theta` | `fmean_soil_mstpot` | DERIVE |
| `soil_wetness_site` | – | soil | tmean | D M Y | `(θ−θ_res)/(θ_sat−θ_res)` | `fmean_soil_wetness` | DERIVE |
| `soil_fliq_site` | – | soil | tmean | D M Y | `soil_e%soil_fliq` | `fmean_soil_fliq` | READY |
| `root_uptake_layer_site` | kg/m²/s | soil | tmean | D M | `hydro_flux%root_uptake_layer` | `fmean_wflux_gw_layer` | PLUMB |
| `soil_water_patch` | m³/m³ | patch | tmean | M | column mean, reduce | — | NEW |

**2-D per-patch soil columns** (`DIM_SOIL_PATCH`, §3.6.1, `axes_soil_patch = true`, streams `D M`):

| name | long_name | units | agg | source | status |
|---|---|---|---|---|---|
| `soil_water_layer_patch` | volumetric soil moisture by layer and patch | m³/m³ | tmean | `patch%soil_w%theta` | NEW (axis only) |
| `soil_psi_layer_patch` | soil matric potential by layer and patch | MPa | tmean | `soil_psi_from_theta` | DERIVE |
| `soil_temp_layer_patch` | soil temperature by layer and patch | K | tmean | `patch%soil_e%soil_temp` | NEW (axis only) |
| `soil_fliq_layer_patch` | soil liquid fraction by layer and patch | – | tmean | `patch%soil_e%soil_fliq` | NEW (axis only) |
| `soil_wetness_layer_patch` | relative saturation by layer and patch | – | tmean | derived | DERIVE |
| `root_uptake_layer_patch` | root water uptake by layer and patch | kg/m²/s | tmean | `hydro_flux%root_uptake_layer` | PLUMB (needs P2) |

The first four are `NEW (axis only)`: the state is already per-patch and prognostic
(`site%patch%soil_w(ip)%theta(:)`), so nothing physics-side has to change — **only the addressing mode
of §3.6.1**. That is why this lands as its own small phase before the fast capture. Only
`root_uptake_layer_patch` waits for the per-cohort/per-patch fast accumulator.
| `w_surface_site` | kg/m² | scalar | tmean | D M | `soil_w%w_surface` | `fmean_sfcw_mass` | PLUMB |
| `swe_site` | kg/m² | scalar | tmean | D M Y | `snow%swe` | `fmean_sfcw_mass` | PLUMB |
| `snow_depth_site` | m | scalar | tmean | D M Y | `snow%snow_depth` | `fmean_sfcw_depth` | PLUMB |
| `snow_temp_site` | K | snow | tmean | D M | `snow%snow_temp` | `fmean_sfcw_temp` | PLUMB |
| `snow_cover_site` | – | scalar | tmean | D M Y | Niu–Yang cover fraction | — | PLUMB |
| `canopy_water_site` | kg/m² | scalar | tmean | D M | `leaf_surf_water + wood_surf_water` | `fmean_leaf_water` | PLUMB |

---

### 4.5 Energy and radiation

| name | long_name | units | dim | agg | streams | source | ED2 | status |
|---|---|---|---|---|---|---|---|---|
| `le_site` | latent heat flux | W/m² | scalar | tmean | F D M Y | `atm_fluxes` | `fmean_vapor_ac`·λ | READY (rename) |
| `h_site` | sensible heat flux | W/m² | scalar | tmean | F D M Y | `atm_fluxes` | `fmean_sensible_ac` | READY (rename) |
| `rnet_site` | net all-wave radiation | W/m² | scalar | tmean | F D M Y | absorbed SW+LW sum | `fmean_rnet` | READY (rename) |
| `ground_heat_site` | ground heat flux (G) | W/m² | scalar | tmean | F D M | soil top-face flux | `fmean_sensible_gg` | PLUMB |
| `sw_in_site` | incident shortwave | W/m² | scalar | tmean | F D M Y | `ctx%rad_sw_top` | `fmean_atm_rshort` | READY (rename) |
| `sw_up_site` | reflected shortwave | W/m² | scalar | tmean | F D M | `rad_flux%albedo`·SW | `fmean_rshortup` | PLUMB |
| `lw_in_site` | incident longwave | W/m² | scalar | tmean | F D M | met forcing | `fmean_atm_rlong` | PLUMB |
| `lw_up_site` | outgoing longwave | W/m² | scalar | tmean | F D M | two-stream upward TIR | `fmean_rlongup` | PLUMB |
| `albedo_site` | all-wave shortwave albedo | – | scalar | tmean | D M Y | `rad_flux%albedo` | `fmean_albedo` | PLUMB |
| `albedo_band_site` | band albedo (VIS/NIR) | – | band | tmean | D M | `rad_flux%albedo(:)` | `fmean_albedo_par/nir` | PLUMB |
| `sw_ground_site` | shortwave reaching ground | W/m² | scalar | tmean | D M | `forc%abs_sw_ground` | `fmean_rshort_gnd` | PLUMB |
| `par_ground_site` | PAR reaching ground | W/m² | scalar | tmean | D M | `rad_flux%dn_ground(VIS)` | `fmean_par_gnd` | PLUMB |
| `ustar_site` | friction velocity | m/s | scalar | tmean | F D M | `aero%ustar` | `fmean_ustar` | READY (rename) |
| `tstar_site` `qstar_site` `cstar_site` | turbulent scales | K, kg/kg, µmol/mol | scalar | tmean | F D | `aero%` | `fmean_tstar` etc. | PLUMB |
| `zeta_site` | M–O stability parameter | – | scalar | tmean | F D | `aero%zeta` | — | PLUMB |
| `ggnet_site` | ground conductance | m/s | scalar | tmean | F D | `aero%ggnet` | `fmean_can_ggnd` | PLUMB |
| `rough_site` `displace_site` | roughness, displacement | m | scalar | tmean | D M | `aero%` | `fmean_rough`, `fmean_veg_displace` | PLUMB |
| `resid_energy_site` | column energy-budget residual | J/m² | scalar | fluxsum | F D M Y | `budg%whole_energy` | `energy_residual` | PLUMB |

**Energy states:**

| name | units | dim | agg | streams | source | ED2 | status |
|---|---|---|---|---|---|---|---|
| `cas_temp_site` | K | scalar | tmean | F D M Y | `cas%can_temp` | `fmean_can_temp` | READY |
| `cas_shv_site` | kg/kg | scalar | tmean | F D M | `cas%can_shv` | `fmean_can_shv` | PLUMB |
| `cas_co2_site` | µmol/mol | scalar | tmean | F D M Y | `cas%can_co2` | `fmean_can_co2` | PLUMB |
| `cas_vpd_site` | Pa | scalar | tmean | F D M | derived from twins | `fmean_can_vpdef` | DERIVE |
| `cas_depth_site` | m | scalar | tmean | M Y | `cas%can_depth` | — | PLUMB |
| `soil_temp_site` | K | soil | tmean | F D M Y | `soil_e%soil_temp` | `fmean_soil_temp` | READY |
| `soil_temp_top_site` | K | scalar | tmean | F D M Y | layer 1 | `fmean_gnd_temp` | READY |
| `ground_temp_site` | K | scalar | tmean | F D M | ground-skin balance | `fmean_gnd_temp` | PLUMB |
| `leaf_temp_site` | K | scalar | tmean | F D M | leaf-area-weighted reduce | `fmean_leaf_temp` | PLUMB |
| `wood_temp_site` | K | scalar | tmean | F D M | wood-area-weighted reduce | `fmean_wood_temp` | PLUMB |
| `air_temp_site` | K | scalar | tmean | F D M Y | met forcing | `fmean_atm_temp` | READY (rename) |
| `atm_shv_site` `atm_prss_site` `atm_vels_site` `atm_co2_site` | forcing echo | scalar | tmean | D M | met forcing | `fmean_atm_*` | PLUMB |

Echoing the forcing into the diagnostic file is not redundant: it makes each output file
self-contained for flux-tower comparison, which is the primary v0.1 evaluation workflow.

---

### 4.6 Demography and structure

**Cohort axis** (`GRP_STRUCTURE`, streams `D M`):

| name | units | agg | source | status |
|---|---|---|---|---|
| `nplant_cohort` `dbh_cohort` `height_cohort` `basal_area_cohort` `agb_cohort` `leaf_area_cohort` | — | last/mean | SoA | READY |
| `growth_avg_cohort` `pft_cohort` `owner_patch` `global_cohort_id` | — | mean/last | SoA | READY |
| `lai_cohort` | m²/m² | tmean | `nplant·leaf_area` | DERIVE |
| `leaf_carbon_cohort` `fineroot_carbon_cohort` `wood_carbon_cohort` `storage_carbon_cohort` | kgC/plant | tmean | SoA | READY |
| `sla_cohort` `vcmax25_cohort` `rd25_cohort` `llspan_cohort` | traits | tmean | SoA (dynamic traits) | READY |
| `overtopping_lai_cohort` | m²/m² | tmean | SoA | READY |
| `dbh_growth_cohort` | cm/yr | tmean | `deriv%d_dbh_dt` | READY (§3.5) |
| `agb_growth_cohort` | kgC/plant/yr | tmean | `deriv%d_agb_dt` | READY (§3.5) |
| `mort_rate_cohort` | 1/yr | tmean | `−deriv%dln_nplant_dt` | READY (§3.5) |
| `dmax_psi_leaf_cohort` | MPa | tmean | SoA | READY |
| `pheno_flush_cohort` `pheno_shed_cohort` | – | tmean | SoA | READY |

**Patch axis** (streams `D M`):

| name | units | source | status |
|---|---|---|---|
| `area_patch` `age_patch` `dist_type_patch` `cohort_offset` `cohort_count` `global_patch_id` | — | SoA | READY |
| `lai_patch` `agb_patch` `nplant_patch` `basal_area_patch` `n_cohort_patch` | — | reduce | NEW |
| `canopy_height_patch` | m | max cohort height | NEW |

**Site axis:**

| name | units | streams | source | status |
|---|---|---|---|---|
| `nplant_site` `basal_area_site` `agb_site` `lai_site` | — | D M Y | reduce | READY |
| `wai_site` | m²/m² | D M Y | wood area index | NEW |
| `n_cohort` `n_patch` | – | D M Y | counts | **dead code — register** |
| `canopy_height_site` | m | D M Y | tallest cohort | DERIVE |
| `mean_dbh_site` | cm | M Y | BA-weighted | READY (unregistered) |
| `agb_growth_site` `agb_mort_site` `agb_recruit_site` | kgC/m²/yr | M Y | `deriv` + recruit/cull hooks | NEW |
| `ba_growth_site` `ba_mort_site` `ba_recruit_site` | m²/m²/yr | M Y | ditto | NEW |
| `nplant_mort_site` `nplant_recruit_site` | plant/m²/yr | M Y | ditto | NEW |
| `disturbance_area_site` | 1/yr | Y | treefall area flux | NEW |

**PFT axis** (`DIM_PFT`, streams `M Y`): `agb_pft`, `lai_pft`, `nplant_pft`, `basal_area_pft`,
`gpp_pft`, `npp_pft`, `agb_growth_pft`, `agb_mort_pft`, `agb_recruit_pft`. All NEW; all one
`reduce_cohort_to_pft` call each.

**Size-class axis** (`DIM_SIZE`, streams `M Y`): `nplant_size`, `agb_size`, `basal_area_size`,
`lai_size`, `agb_growth_size`, `agb_mort_size`, `agb_recruit_size`. Class edges configurable
(`[output].dbh_class_edges = [0, 10, 20, 30, 50, 70, 100]`), emitted as a coordinate variable so the
file is self-describing. All NEW.

**Binning follows ED2: a size class of *plants*, not of cohorts** (§8 D4). Each cohort is assigned to
the single bin containing its mean `dbh` and contributes with weight `area · nplant`; a cohort is never
split across bins. So `nplant_size` is a genuine stem-density distribution directly comparable to a
forest inventory, and `Σ_class nplant_size == nplant_site` exactly. The same
`reduce_cohort_to_size(..., wkind = W_NPLANT_AREA)` call serves every row.

---

### 4.7 Individual-level ecophysiology — `GRP_ECOPHYS`, cohort axis

**This is the group the request singles out and the group that does not exist at all today.** Most rows
are `PLUMB` via the §3.4 `cohort_diag_block`; default streams `F D M` (off unless `ecophys = true`).

Four rows are marked `NEW†` because the quantity is a **local inside `leaf_gas_exchange`** and is not on
`leaf_flux_t` at all: `beta_stomata` / `beta_nonstomata` (computed at
`src/plant/meds_leaf_gas_exchange.f90:216-221` and discarded) and the three potential rates behind
`limitation`. Surfacing them is a small, self-contained addition of four fields to `leaf_flux_t` — worth
doing, because a water-stress closure whose two limbs are never observable is exactly the kind of thing
that stays silently inert (it already did once: `env%psi` defaulted to 0 and `beta_stomata` was
identically 1 until issue #95).

| name | long_name | units | agg | source | ED2 |
|---|---|---|---|---|---|
| `anet_cohort` | net leaf assimilation | µmol CO₂/m²leaf/s | tmean | `leaf_flux%A_net` | `fmean_a_net` |
| `agross_cohort` | gross leaf assimilation | µmol CO₂/m²leaf/s | tmean | `leaf_flux%A_gross` | — |
| `gsw_cohort` | stomatal conductance to H₂O | mol H₂O/m²leaf/s | tmean | `leaf_flux%gs` | `fmean_leaf_gsw` |
| `gsc_cohort` | stomatal conductance to CO₂ | mol CO₂/m²leaf/s | tmean | `gsw/1.6` | — |
| `gbw_cohort` | leaf boundary-layer conductance | mol H₂O/m²leaf/s | tmean | `aero%leaf_gbw` | `fmean_leaf_gbw` |
| `ci_cohort` | intercellular CO₂ | µmol/mol | tmean | `leaf_flux%ci` | — |
| `cs_cohort` | leaf-surface CO₂ | µmol/mol | tmean | `leaf_flux%cs` | — |
| `ci_ca_cohort` | Cᵢ/Cₐ ratio | – | tmean | derived | — |
| `rd_cohort` | leaf dark respiration | µmol/m²leaf/s | tmean | `leaf_flux%rd` | `fmean_leaf_resp` |
| `transp_cohort` | leaf transpiration | mol H₂O/m²leaf/s | tmean | `leaf_flux%transpiration` | `fmean_transp` |
| `wue_cohort` | water-use efficiency | µmol CO₂/mmol H₂O | tmean | `A_net/transp` | — |
| `limitation_cohort` | binding photosynthesis limit | – (enum) | last | `leaf_flux%limitation` | `fmean_a_rubp/light/co2` |
| `a_rubp_cohort` `a_light_cohort` `a_co2_cohort` | the three potential rates | µmol/m²/s | tmean | FvCB internals — **NEW†** | `fmean_a_*` |
| `leaf_temp_cohort` | leaf temperature | K | tmean | `bio%leaf_temp` | `fmean_leaf_temp` |
| `wood_temp_cohort` | wood temperature | K | tmean | `bio%wood_temp` | `fmean_wood_temp` |
| `leaf_vpd_cohort` | leaf-to-air VPD | Pa | tmean | leaf env | `fmean_leaf_vpdef` |
| `psi_leaf_cohort` | leaf water potential | MPa | tmean | `hydro_flux%psi_leaf` | `fmean_leaf_psi` |
| `psi_wood_cohort` | wood water potential | MPa | tmean | `hydro_flux%psi_wood` | `fmean_wood_psi` |
| `plc_cohort` | percent loss of conductance | – | tmean | `hydro_flux%plc` | — |
| `sapflow_cohort` | wood→leaf sapflow | kg/plant/s | tmean | `hydro_flux%sapflow` | `fmean_wflux_wl` |
| `root_uptake_cohort` | soil→root uptake | kg/plant/s | tmean | `hydro_flux%root_uptake` | `fmean_wflux_gw` |
| `leaf_water_cohort` `wood_water_cohort` | internal tissue water | kg/plant | tmean | SoA | `fmean_leaf_water_int` |
| `beta_stomata_cohort` | stomatal water-stress factor | – | tmean | leaf kernel local — **NEW†** | `fmean_fs_open` |
| `beta_nonstomata_cohort` | non-stomatal stress factor | – | tmean | leaf kernel local — **NEW†** | — |
| `abs_par_cohort` | absorbed PAR | W/m² ground | tmean | `forc%abs_par` | `fmean_par_l` |
| `abs_sw_cohort` | absorbed shortwave | W/m² ground | tmean | `forc%abs_sw` | `fmean_rshort_l` |
| `abs_lw_cohort` | net longwave | W/m² ground | tmean | `forc%abs_lw` | `fmean_rlong_l` |
| `light_level_cohort` | relative light level | – | tmean | RT profile | `fmean_light_level` |
| `wind_cohort` | in-canopy wind speed | m/s | tmean | `aero%wind` | — |

Site/patch/PFT twins of `anet`, `gsw`, `psi_leaf`, `plc` (leaf-area-weighted) come free from the
reducer and are worth registering — a leaf-area-weighted canopy `gsw` and canopy-mean `ψ_leaf` are
standard flux-tower comparison quantities.

---

### 4.8 Numerics and budget health — `GRP_NUMERICS`

Existing 11 work counters stay. Add: `resid_energy_site`, `resid_water_site`, `resid_co2_site`,
`carbon_resid_site` (moved here from their physical groups if `numerics` is on),
`theta_ood_max_site`, `adapt_dt_site` (mean accepted controller step), `hydro_nonconv_site`.
Recommendation: **flip `GRP_NUMERICS` and the three residuals ON by default** for v0.1. A budget
residual that nobody looks at is the exact failure mode
`feedback_conservation_is_not_stability` warns about — but its inverse, a run whose closure was never
recorded at all, is worse.

### 4.9 Counts

| Group | New | Total after |
|---|---|---|
| STRUCTURE | ~35 | ~55 |
| CARBON | ~28 | ~38 |
| WATER | ~32 | ~35 |
| ENERGY | ~24 | ~28 |
| RADIATION | ~8 | ~8 |
| ECOPHYS | ~30 | ~30 |
| BIOGEOCHEM | ~9 | ~17 |
| NUMERICS | ~7 | ~18 |
| **Total** | **~173** | **~229** |

`MAX_OUTPUT_VARS` → 512.

---

## 5. TOML configuration surface

### 5.1 Main config `[output]` — extended

```toml
[output]
enabled     = true
dir         = "out"
prefix      = "run01"
io_config   = "meds_io_config.toml"   # optional per-variable overrides (wins over everything below)
cohort_max  = 4096
patch_max   = 256
strict_caps = false
sync_every  = "flush"
dbh_class_edges = [0.0, 10.0, 20.0, 30.0, 50.0, 70.0, 100.0]   # NEW: DIM_SIZE bin edges [cm]

  # --- variable-group toggles (8 groups) ---
  structure  = true
  carbon     = true
  water      = true
  energy     = true
  radiation  = false
  ecophys    = false      # per-cohort leaf gas exchange + hydraulics (high volume)
  biogeochem = true
  numerics   = true       # includes the budget residuals -- recommended ON

  # --- axis toggles: suppress a whole dimension without naming variables ---
  axes_cohort     = true  # per-cohort slabs (the largest output by volume)
  axes_patch      = true
  axes_pft        = true
  axes_size       = true
  axes_soil_patch = false # 2-D (soil layer x patch) columns -- section 3.6.1, high volume

  [output.fast]
  enabled = false ; interval_steps = 4 ; file_chunk = "day"
  [output.daily]
  enabled = true  ; file_chunk = "month"
  [output.monthly]
  enabled = true  ; file_chunk = "year"
  [output.annual]
  enabled = true  ; file_chunk = "run"
```

**New:** the `axes_*` toggles. With ~55 cohort-dimensioned variables, "give me everything at site level
and nothing per-cohort" is the single most common request and is currently inexpressible. Resolution
order becomes: registry defaults → **axis toggles** → group toggles → per-tier enable → per-variable
overrides.

### 5.2 `meds_io_config.toml` — ship a real one

The per-variable override loader is implemented but **no example file exists in the repo**. Ship a
fully-commented `meds_io_config.toml` at the repo root listing **every registered variable** with its
default stream mask commented out, generated by a `--dump-io-config` flag on `meds_main`. That flag is
also the discoverability fix for F8 and the regression guard that the registry and the documentation
cannot drift apart.

```toml
[variables]
# value = true/false        -> force on (at registry defaults) / off everywhere
# value = "F D M Y"         -> replace the stream mask (Y forbidden on cohort/patch vars)
anet_cohort   = "F D"       # sub-daily leaf physiology for one diagnostic run
gsw_cohort    = "F D"
agb_size      = "M Y"
growth_avg_cohort = false
```

### 5.3 Retire `[io]`

Move `write_state` / `state_interval_years` into a `[state]` block (as the design doc already
specifies), delete `write_output` / `output_interval_years` / `output_dir` / `output_prefix` /
`cohort_max` / `patch_max`, delete `io_create`/`io_write_snapshot`/`io_close` from `meds_io.f90` and
the `is_new_year` write block from `meds_main`. **Then rename the aggregation subsystem's daily stream
files only after `[io]` is gone**, so the `-D-` collision never exists in a released version.

---

## 6. Phasing

Each phase leaves the model runnable and is independently shippable. `nvfortran` multicore must be
built for every phase (a green ifx suite is not sufficient — issue #7), and the threaded fast loop
must stay byte-identical at 1/2/4/8 threads.

### P0 — the reduction layer (no new state, no new physics)
- New `src/io/meds_diagnostic_kernels.f90` (§3.2) and `src/io/meds_diagnostic_reduce.f90` (§3.3),
  **added to the `meds_io_prep` explicit source list** (it is not a GLOB); delete
  `meds_output_diagnostics.f90` from that list once ported. Port every existing `total_*` to the
  generic reducer and **assert bit-identical output** against the current registry as the acceptance
  test.
- Wire `DIM_PFT` and `DIM_SIZE` in `meds_output_stream` + `[output].dbh_class_edges`.
- Register everything reachable from existing state: the derived cohort/patch/site/PFT/size structure
  variables, the carbon stocks, `n_cohort`/`n_patch`, `mean_dbh_site`, `canopy_height_site`,
  `soil_psi_site`, `soil_wetness_site`, `cas_vpd_site`, the `deriv`-sourced growth/mortality rates.
- Registry restructure (§3.7), `MAX_OUTPUT_VARS` → 512, slab sizing fix.
- **Deliverable:** ~110 variables, PFT and size-class axes live, zero physics-side changes.
  This alone covers most of the *slow-scale* checking list in §4.0.

The phase graph is `P0 → {P1, P2} → P3 → P4 → P5`: P1 and P2 are independent of each other, so if the
2-D soil output is wanted early it does not have to queue behind the invasive fast-capture work.

### P1 — the 2-D (soil layer × patch) axis
Small, self-contained, and **independent of P2** — it can be built in parallel with the fast capture
because it touches only the addressing mode and the serializer, no physics and no state.
- `DIM_SOIL_PATCH` + the flattened index (§3.6.1); `max_slab` grows to include
  `patch_max · n_soil_layer_max`; the rank-3 `start`/`count` case in `stream_write_record`; the
  `(time, patch, soil)` definition + chunking in `stream_open_file`.
- `gather_patch_columns` in the reducer; the `[output].axes_soil_patch` toggle (default off).
- Register `soil_water_layer_patch`, `soil_psi_layer_patch`, `soil_temp_layer_patch`,
  `soil_fliq_layer_patch`, `soil_wetness_layer_patch` (§4.4). `root_uptake_layer_patch` is registered
  here but stays `_FillValue` until P2 supplies it.
- **Deliverable:** per-patch soil drydown and thermal profiles — the gap-vs-closed-canopy contrast the
  area-weighted site column hides. **Acceptance test:** the area-weighted mean over the patch axis of
  `soil_water_layer_patch` reproduces `soil_water_site` to roundoff, which is what proves the flattened
  index is right.

### P2 — the fast diagnostic capture (the big one)
- `cohort_diag_block` + `patch_diag_block` (§3.4): types, allocation gated on the registry, dt-weighted
  fill in `fast_dynamics`/`column_fast_step`, reset per slow step, **lockstep reorder + fusion rules**,
  and the extensive/intensive fusion blend.
- Unify the extraction path: `output_integrate` folds every tier including FAST through one
  `extract_variable`; delete `fast_sample_t`, `extract_fast_scalar`, `output_integrate_fast` and the
  bespoke staging arrays. **F3 closes here.**
- Register `GRP_ECOPHYS` (§4.7), the water/energy flux partitions (§4.4/§4.5), the budget residuals.
- Surface the four `NEW†` leaf-kernel locals (§4.7) as fields on `leaf_flux_t`.
- Fill `root_uptake_layer_patch` (the one P1 row left as fill).
- **Deliverable:** the full fast-scale checking list; ~190 variables; the four-edit-site tax is back to
  two.

### P3 — the slow-loop rates
- `slow_diag_t` per patch (§3.5): NPP allocation components, litter fluxes, recruitment,
  mortality carbon by pathway, disturbance area flux.
- Register the census-style demographic rates (`agb_growth/mort/recruit`, `ba_*`, `nplant_*`) on the
  site, PFT and size axes.
- **Deliverable:** the carbon and demographic budgets are closeable *from the output file alone* —
  `d(agb)/dt ≈ growth − mortality + recruitment` becomes a testable assertion, not a hope.

### P4 — config, cleanup, provenance
- The extended `[output]` block with axis toggles; `--dump-io-config`; the shipped
  `meds_io_config.toml`. **This is what delivers the per-variable debugging toggle (§8 D6): the
  override mechanism already exists and works — what is missing is a generated file that lists every
  variable, so a user can find the one they want to switch on.**
- Delete the legacy `[io]` diagnostic writer, move `[state]` out, resolve the `-D-` collision.
- Retire the `_fast` name suffixes.
- Provenance globals the design specifies but that should be re-verified: git hash, config echo +
  hash, compiler/flags, CF `Conventions`.
- Update `post_proc/` readers (`plot_site_timeseries.py`, `plot_forest_structure.py`) to the new
  stream layout and add a PFT/size-class plotter.

### P5 — documentation and evaluation
- `docs/science/diagnostics.md`: the full variable table, the extensive/intensive aggregation rules,
  the `cell_methods` semantics, and the ED2 name mapping (so an ED2 user can find their variable).
- A worked evaluation notebook against the Ithaca 30-yr test bed: diurnal energy/carbon cycles from the
  FAST tier, seasonal cycles from DAILY, demographic trajectories from MONTHLY/ANNUAL.

---

## 7. Test plan

Extending `test_output_integrate` / `test_output_registry` / `test_output_roundtrip`:

1. **Reducer equivalence (P0 gate) — RESULT.** Byte-identical across all 75 output files of a 3-year
   run at `-O0`. At `-O2` six variables differ by 1–2 ULP (max 3.8e-16 relative). The gate is therefore
   stated as: **exact without FMA contraction, ≤4 ULP with it.** Demanding bit-identity at `-O2` across
   a refactor is not achievable in general — the compiler is free to contract mathematically equivalent
   expressions differently — and pretending otherwise would just have meant loosening the tolerance
   later without saying why.
2. **Extensive/intensive correctness.** A two-cohort patch with known `nplant`/`leaf_area`: assert
   `agb_site = Σ area·nplant·agb` (extensive) and `leaf_temp_site` = leaf-area-weighted mean
   (intensive), and that swapping the weight kind changes the answer — so a wrong weight cannot pass.
3. **PFT/size closure.** `Σ_pft agb_pft == agb_site` and `Σ_class agb_size == agb_site` to roundoff;
   a cohort exactly on a class edge lands in **exactly one** bin (assert the total stem count is
   conserved, which is what catches a half-open/closed interval error).
3b. **2-D soil closure (P1 gate).** The patch-area-weighted mean over the patch axis of
   `soil_water_layer_patch` reproduces `soil_water_site` layer-by-layer to roundoff. Separately, drive
   a fixture with **different layer counts and a changed patch count between windows** and assert the
   flattened index still addresses `(k, ip)` correctly — the `n_soil_layer_max` stride, not a live
   count, is the thing being tested.
4. **Empty-set guard.** A patch with zero cohorts, a PFT with zero members, a size class with zero
   members → `_FillValue`, never NaN, never 0.
5. **Fast-capture conservation (P1 gate).** Over one slow step,
   `Σ_substep anet_cohort·leaf_area·dt` reconciles with `gpp_accum − leaf_resp_accum` to roundoff.
   This is the test that catches a dt-weighting or unit error in the whole `cohort_diag_block`.
6. **Fusion blend.** Fuse two cohorts with different `psi_leaf`/`gsw` mid-window; assert the survivor's
   diagnostic is the leaf-area-weighted blend, and that `Σ nplant·anet` is conserved.
7. **Lockstep integrity.** Sort/fuse/split/cull/recruit/disturb with the diag blocks allocated; assert
   every diag field permutes in lockstep with `global_id` (the standard MEDS reorder test, extended).
8. **Thread invariance.** Full output byte-identical at 1/2/4/8 threads with `ecophys = true`.
9. **Budget closure from file.** P2: read the monthly file and assert
   `Δagb ≈ (growth − mort + recruit)·Δt` within the fusion tolerance, and
   `Δsoilc ≈ litter − Rh`.
10. **Config resolution order.** Axis toggle vs group toggle vs tier enable vs per-variable override,
    all four in one fixture, asserting the documented precedence; typo'd key → `error stop` naming it.
11. **Performance.** A multi-decade daily run: wall time with `structure+carbon+water+energy` on within
    a few percent of output-off; with `ecophys = true` + FAST tier, quantify the cost and **report it**
    (do not assume it is free — ~30 cohort variables × 4096 slots × 4 tiers is the one configuration
    that can matter).
12. **Round-trip metadata.** Every variable carries `cell_methods`, `_FillValue`, `units`,
    `long_name`; averaged records carry `time_bnds`; the size-class axis carries its edges.

---

## 8. Risks and resolved decisions

### Risks

**R1 — the `cohort_diag_block` lockstep is the highest-risk change in this plan.** It adds ~25 arrays
to the machinery that CLAUDE.md explicitly names as the fix for ED2's "forgot to reallocate" class of
bug. Mitigation: one `contains`-local helper that lists the diag arrays once, called from every reorder
site, plus test 7. Do not spread the array list across `copy_cohort_slot`, `cohort_compact`,
`cohort_ensure_capacity` and `move_alloc_block` by hand.

**R2 — memory.** ~55 cohort-dimensioned variables × 4 tiers × 4 buffer arrays × `cohort_max` doubles is
~115 MB at `cohort_max = 4096`, and P1 raises `max_slab` further to `patch_max · n_soil_layer_max`
(5120 at `patch_max = 256`) for the 2-D soil axis. The slab-sizing fix (live cap, not `cohort_max`) is
therefore **not optional at P0** — it is a prerequisite for both P1 and P2.

**R3 — `AGG_SUM` → `AGG_FLUXSUM` reclassification changes existing output values.** `gpp_site`,
`npp_site`, `et_site`, `rh_site` are currently `AGG_SUM` over slow-step accumulators. Changing them is
right but breaks any existing analysis. Do it once, at P2, and note it in the release notes.

**R4 — the flattened `DIM_SOIL_PATCH` index is a silent-corruption risk if the stride is ever taken
from a live count.** Striding by `n_soil_layer_max` (a compile-time constant) rather than the active
layer count is what makes the mapping stable; getting it wrong produces plausible, shifted profiles
rather than an error. Test §7.3b exists specifically to catch this and should vary both the layer count
and the patch count.

### Decisions (settled 2026-08-02; no longer open)

**D1 — the fast diag block is per-cohort, NOT per-(cohort, sub-step).** One dt-weighted accumulator
value per cohort per slow step; a per-sub-step array would be `ncohort × nsub × nvar` and only the FAST
tier could consume it. **Documented consequence: a FAST-tier per-cohort variable is a mean over the
fast window, not an instantaneous sub-step value** — identical for `fast_interval_steps = 1`, a
4-sub-step mean for `= 4`. This must be stated in the variable's `cell_methods` and in
`docs/science/diagnostics.md`, because "hourly `A_net`" that is silently a 4-step mean is the kind of
thing that misleads a reader comparing against a flux tower.

**D2 — the `pft` axis length is the run-time PFT count** (`cfg%pft%n_pft`), not a compile-time ceiling.
Because that makes the dim run-dependent, each file must carry a `pft` **coordinate variable** with the
PFT indices plus the per-PFT wood density (the model's PFT contrast axis) as an attribute, so a file
remains self-describing when compared across runs with different PFT tables.

**D3 — 2-D `(soil layer, patch)` output is IN, as its own phase P1.** Design in §3.6.1 (flattened
index, no new buffer type, no new C binding), variables in §4.4, acceptance test in §7.3b. Gated by
`[output].axes_soil_patch`, default off.

**D4 — `DIM_SIZE` is an ED2-style size class of *plants*.** Bin by the cohort's mean `dbh`, weight by
`area · nplant`, never split a cohort across bins. Detailed in §4.6.

**D5 — single-site only for v0.1.** Multi-site/polygon output is explicitly out of scope; no
placeholder axis, no reserved leading dimension, no speculative generality. When multi-site lands it
brings its own file-layout decision.

**D6 — v0.1 default configuration** is as proposed: `structure + carbon + water + energy + biogeochem +
numerics` on; `radiation + ecophys` off; `axes_cohort/patch/pft/size` on, `axes_soil_patch` off; FAST
tier off; daily/monthly/annual on. That emits the §4.0 checking list out of the box at modest file
size, and each heavy group is one boolean away.

**Per-variable toggling for debugging is a stated requirement.** It is worth being precise about what
already exists: the `meds_io_config.toml` mechanism — `name = true|false` or `name = "F D M Y"`, with
an unknown-key `error stop` that names typos — **is implemented and working today**
(`apply_io_overrides`, resolution step 4). What is missing is purely discoverability: no example file
ships, and the `[output]` block in `meds_config_main.toml` is entirely commented out. P4's
`--dump-io-config` generator closes that, and it doubles as the regression guard that the registry and
its documentation cannot drift apart. So this is a **documentation-and-ergonomics** task, not new
machinery.
