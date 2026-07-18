# MEDS design: extract the hydraulics constitutive curves into `shared/functions/meds_hydro_curve`

**Status:** design-only (review before any code lands). **Date:** 2026-07-18.
**Scope:** a behaviour-preserving relocation (Path A) plus the enabling seam for consolidated
parameter derivation. **No numerical change** to any result is intended or expected.

## 1. Motivation

Two threads converge here:

1. **Consolidated parameter derivation.** We want a single `derive_parameters(cfg)` routine that
   installs *every* derived/precomputed parameter — allometry coefficients, per-PFT rate
   derivations, and (new) the per-PFT hydraulic vulnerability-curve lookup tables — callable from
   both `load_meds_config` and `build_test_config`. The DAG blocks this today only because the
   table build needs `flux_potential`, which lives in `plant`, while the natural home for the
   orchestrator is `shared/config` (which cannot `use` `plant`).

2. **The functions that block it were mis-filed.** `flux_potential` and the whole vulnerability /
   pressure-volume family are **stateless, `elemental pure`, scalar-in/scalar-out constitutive
   relations** — material properties of tissue, structurally identical to the allometry functions
   that already live in `src/shared/functions/`. They are not "plant the process"; only the
   *solver* that assembles them into a network ODE is.

**Path A** dissolves the blocking dependency: move the constitutive curves down into
`src/shared/functions/meds_hydro_curve.f90` (next to `meds_allometry`), leave the solver in
`src/plant/meds_plant_hydraulics.f90`, and then `derive_parameters` in `shared/config` can build
the tables directly — with no procedure-pointer injection and no dependency inversion.

Precedent: the soil column already separates constitutive curves (`meds_soil_parameters`) from the
integrator (`meds_soil_solver`). This is the same split one tier lower, pushed into
`shared/functions` specifically so config-load can install the tables — exactly as `set_allometry`
is installed today.

## 2. What moves, what stays (verified against the current file)

### 2a. → NEW `src/shared/functions/meds_hydro_curve.f90` (constitutive; `elemental pure`)

**Vulnerability / conductance family**
| symbol | current line | notes |
|---|---|---|
| `plc_retained` | 59 | 1 − PLC |
| `dplc_dpsi` | 67 | d(plc)/dψ (public; no current external caller) |
| `flux_potential` | 77 | Kirchhoff Φ (closed form kexp∈{1,2}; 7-pt GL else) |
| `phi_inverse` | 100 | bisection inverse of Φ (public; no current external caller) |
| `kirchhoff_edge` | 116 | edge conductance from ΔΦ/Δψ |

**Pressure-volume family**
| symbol | current line |
|---|---|
| `pv_psi_tlp` | 147 |
| `pv_rwc_tlp` | 153 |
| `psi_from_rwc` | 159 |
| `rwc_from_psi` | 172 |
| `water_content` | 190 |
| `capacitance` | 199 |
| `pv_water_cap_from_traits` | 215 |

**Module parameters that move** (used only by the above): `dpsi_eps` (21), `NG` (22),
`gl_x` (25), `gl_w` (29), `rwc_floor` (39).

**New members** (the precompute this enables): `hydro_table_t` (uniform-grid lookup: `r_max`,
`inv_h`, `n`, `g(:)`), `build_hydro_table(kexp, r_max, n) -> hydro_table_t` (uses `flux_potential`
at `psi50=1` to fill `G(r)`), and the runtime lookup `flux_potential_lin(psi, psi50, tab)`
(linear interp — chosen over cubic per the CPU-gather benchmark). These are additive and unused on
the hot path until a general-exponent PFT is configured; today's kexp=2 path is unaffected.

**Dependencies of the new module:** `use meds_kinds, only : wp, ik` — nothing else. No
`meds_constants`, no `meds_plant_types`. It is a pure leaf of the DAG, cleaner than `meds_allometry`
(which also needs `meds_constants`). **No cycle is possible.**

### 2b. → STAYS in `src/plant/meds_plant_hydraulics.f90` (the solver; touches `hydro_*_t`)

- `solve_plant_water` (227) + nested `freeze_coeffs`, `apply_expm`, `expm_step`, `cond_max`
- `plant_water_tendency` (375) — the IMEX-ARK RHS
- `rhizosphere_cond` (131) — kept here per decision (a soil↔root boundary-condition helper, not a
  tissue constitutive curve; still uses `pi`/`tiny_num` from `meds_constants`)
- Module parameters that stay (solver/stepper): `c_floor` (45), `k_floor` (46), `sinhc_eps` (47),
  `safety` (48), `fmax` (49), `fmin` (50), `h_floor` (51)
- The module keeps its name `meds_plant_hydraulics` and its `use meds_plant_types` /
  `use meds_constants` (for `pi, grav_head, safe_exp, tiny_num`).

**One added dependency:** the solver now imports the curves it consumes —
`use meds_hydro_curve, only : kirchhoff_edge, capacitance, water_content, plc_retained`
(`kirchhoff_edge` in `freeze_coeffs`:334 and `plant_water_tendency`:413; `capacitance` in
`freeze_coeffs`:332-333 **and** `plant_water_tendency`:411-412; `water_content` in
`solve_plant_water`:299-302; `plc_retained` in the `flux%plc` diagnostic:308). `plant → shared` is an existing allowed edge. The solver does **not**
call `flux_potential` directly (only via `kirchhoff_edge`, which moved), so no other imports needed.

## 3. Blast radius (verified)

`grep` for every moving symbol across `src/` + `test/`:

- **Production code:** *none outside the solver*. No `plant`/`biophysics`/`biogeochemistry`/driver
  file calls a PV or vulnerability-curve function today. (Prospective consumers — carbon dynamics'
  sapwood water, a future leaf/wood energy balance's tissue heat capacity — are why the family
  belongs in `shared`, but none exist yet, so the move is low-risk.)
- **C-API / Python:** no references. No binding changes.
- **Tests:** exactly **one** file changes its `use` — `test/test_plant_hydraulics.f90`
  (lines 20–22 import `pv_psi_tlp, rwc_from_psi, psi_from_rwc, water_content, capacitance,
  plc_retained, flux_potential, kirchhoff_edge`). Repoint those to `use meds_hydro_curve`. Its
  `solve_plant_water`/`plant_water_flux` imports (via the interface) are unchanged.

Every other `use meds_plant_hydraulics` site imports only solver symbols and is **untouched**:
`meds_plant_interface.f90:35`, `meds_ark_stepper.f90:26`, `meds_column_derivs.f90:34`,
`test/test_column_derivs.f90:32`.

## 4. The `derive_parameters` seam (the payoff)

A single orchestrator in `shared/config` (renamed from the earlier `init_derived_parameters` per
decision). **Corrected after adversarial review:** the current loader hand-rolls **four**
derivation calls, not three, and `set_allometry` takes **9 global scalars**, not `cfg%pft`:

```fortran
! src/shared/config/meds_config.f90  (or a small sibling meds_config_derive)
subroutine derive_parameters(cfg)
   type(meds_config_t), intent(inout) :: cfg
   ! set_allometry installs 9 GLOBAL coefficients (b1Ht..light_ext), NOT per-PFT. Today they are
   ! config_io LOCALS (from TOML, meds_config_io.f90:446/638-646) and hardcoded LITERALS in
   ! build_test_config — never stored in cfg. To make derive_parameters cfg-only, promote them to
   ! cfg fields (see §7); then:
   call set_allometry(cfg%allom%b1Ht, …, cfg%allom%light_ext)  ! use meds_allometry (new shared edge)
   call derive_config(cfg)                  ! existing
   call derive_pft_rates(cfg%pft)           ! existing — mort_gamma/alpha/beta from wood density
   call derive_leaf_params(cfg%pft)         ! existing — jmax25/tpu25/rd25 (MUST NOT be dropped)
   ! call set_hydro_tables(cfg)             ! Phase B (§8): BLOCKED — wood_kexp is not in cfg yet
end subroutine
```

- **The four `derive_*`/`set_*` routines are mutually order-independent** (each reads only primary
  loaded fields; `set_allometry` installs a global none of the others consume). That independence —
  not mere "duplication" — is what licenses collapsing the two current sequences, which in fact
  differ in order: production runs `derive_config` first (`meds_config_io.f90:658-661`), the test
  runs it last (`meds_test_support.f90:117-121`).
- **`set_hydro_tables` is Phase B, not Phase A.** It would loop PFTs building
  `build_hydro_table(cfg%pft(i)%wood_kexp, …)` — but **`wood_kexp` is not in `cfg`**: `pft_table_t`
  carries no hydraulics traits, `meds_config_t` has no hydro block, and `hydro_params_t` is a
  hardcoded stub in `meds_fast_loop.f90:101-107` populated after config load. So there is nothing
  per-PFT in `cfg` to build a table from until hydraulics traits are plumbed into `cfg` from TOML
  (Phase B, §8). When it lands, `meds_config` will `use meds_hydro_curve` + `meds_pft_params` (both
  in `libmeds_shared`; no cycle) and `set_hydro_tables` runs last, after the `[derived]` override.
- **Loader order is unchanged:** `read → derive_parameters → [derived] override → validate` (validate
  LAST, as today at `meds_config_io.f90:673`). The `[derived]` override block
  (`meds_config_io.f90:663-671`, which pins `mort_gamma/alpha/beta`) **stays in `config_io` after
  `derive_parameters`** — it needs the TOML parser handle, so it cannot move into a shared/config
  routine. §4b explains why `set_hydro_tables` is nonetheless safe there today.
- **Call-path coverage (confirmed):** `load_meds_config` runs on *every* production path — restart
  included (`meds_main.f90:69`, before the `init_mode` select-case) — and the C-API path
  (`meds_demography_capi.f90:66`) also routes through it; `build_test_config` covers tests. No path
  bypasses `derive_parameters`. No `init/` involvement is needed, because the curves now live in
  `shared` (the whole point of Path A).

### 4a. Table storage — DECISION: Option 2 (explicit, in `hydro_params_t`)

**Locked: Option 2.** The table lives on `hydro_params_t` (in `plant`), co-located with the
`wood_kexp`/`wood_psi50` it is built from. Rationale:
- `hydro_params_t` already carries `wood_kexp` and is already handed to `solve_plant_water`, so the
  solver receives the table with no new plumbing, and the staleness guard is a one-line invariant on
  a single struct: `assert tab%kexp == p%wood_kexp` (§4b).
- No hidden global mutable state (unlike Option 1's `protected` module array); matches the repo rule
  and the earlier GPU-residency discussion.
- `meds_plant_types` gains `use meds_hydro_curve, only : hydro_table_t` (`plant → shared`, allowed).

**Refinement (removes the only con):** make `hydro_table_t` a **fixed-size POD value type** —
`g(0:NTAB)` with `NTAB` a module parameter (512, from the benchmark) + scalars `kexp, r_max, inv_h,
n` — **not** an allocatable component. A POD `hydro_table_t` is trivially copyable
(`fro%hydro_p = ccfg%hydro_p` at `meds_column_dynamics.f90:966`) and GPU-mappable (no deep copy /
`declare target` headache), at a fixed ~4 KB per PFT. `build_hydro_table` is a **subroutine**
(`intent(out)` table), not an array-returning function, to dodge the nvfortran array-temp trap.

### 4b. Override & staleness safety (the "overwritten derived hydro parameter" case)

*What if a derived hydro parameter is overridden in the PFT TOML — is the table still correct?*

- **Today the case does not arise for hydraulics.** The table's only input, `wood_kexp`, is a
  *primary* trait (read directly, never computed — it is only ever assigned, e.g.
  `meds_fast_loop.f90:105`), and the `[derived]`/`override_derived` mechanism currently pins **only**
  the mortality Camac params (`mort_gamma/alpha/beta`, `meds_config_io.f90:663-671`). No hydro
  quantity is derived-and-stored (`pv_water_cap_from_traits`'s `water_cap` has *no* production
  consumer), so there is nothing overridable that the table depends on.
- **The table is override-robust by construction.** `G(r)` depends **only on `kexp`**; `psi50`
  enters at lookup time as a runtime multiply (`flux_potential_lin = psi50·G(r)`). So a `psi50`
  override needs **no** rebuild, and the *only* override that could invalidate a table is a `kexp`
  override.
- **Ordering invariant (future-proofing).** If `kexp` ever becomes a *derived* trait (e.g. a
  function of wood density in `derive_pft_rates`) and thus `[derived]`-overridable, then
  `set_hydro_tables` **must run after** that override so it builds from the final `kexp`. Because
  today's mortality override touches nothing the table reads, `set_hydro_tables`' position relative to
  the override is currently immaterial — the invariant is written down so it is honored the moment a
  table input first becomes derivable.
- **Guard against post-derive mutation.** Tests and default-setters poke `wood_kexp` *directly* after
  the config is built (`meds_fast_loop.f90:105`, several tests) — a path no override-ordering rule
  covers. To make any desync fail-fast instead of silently using a stale table, **store the source
  `kexp` inside `hydro_table_t`** and have the solver assert `tab%kexp == p%wood_kexp`. Cheap, and it
  catches both the override-ordering case and direct field mutation. (Dormant while `kexp=2` keeps the
  table unused, but this is the concrete "yes, the scenario is handled.")
  - **As implemented** (`solve_plant_water`, hardened after adversarial review): the guard is
    `inv_h <= 0 .or. |tab%kexp - wood_kexp| > 1e-9`. The `inv_h<=0` clause closes a blind spot where a
    never-built table (`kexp=0`) coincides with a defaulted `wood_kexp=0` (both 0 ⇒ the kexp test alone
    passes); a built table always has `inv_h = NTAB/r_max > 0`. The `pure` ARK RHS
    (`plant_water_tendency`) cannot error-stop, so it carries the precondition as a documented invariant
    enforced at hydro-params assembly (`build_fast_context` builds the table, which rides the `hydro_p`
    copy into the ARK path).
- **Out of scope (pre-existing).** If `water_cap` ever becomes derived *and* overridable while the
  nonlinear `capacitance(psi,…)` keeps computing independently from `pi0/eps/af`, an override could
  desync the linear proxy from the curve. This relocation is behaviour-preserving and neither creates
  nor fixes that; track it separately.

## 5. Behaviour preservation

- All moved functions are copied verbatim (same bodies, same `elemental pure` attributes, same
  parameter values). The only source edits are `use` lines and the module split.
- The new table members are **additive and dormant**: nothing calls `flux_potential_lin` until a
  PFT with `kexp∉{1,2}` is configured and the solver is switched to consult the table. Today's
  `kexp=2` runs take the identical `atan` closed form. Golden anchors (e.g. `tc_split` values,
  `test_plant_hydraulics` assertions) must remain bit-identical; that is the acceptance test.
- Compiler coverage: build **ifx Debug** (`-stand f18 -check all`) *and* **nvfortran multicore**
  (the array-temp portability trap is compiler-specific; a green ifx run is not sufficient).

## 6. Migration steps (for the implementation PR, not this doc)

1. Create `src/shared/functions/meds_hydro_curve.f90`; move the §2a functions + params verbatim;
   add `hydro_table_t` + `build_hydro_table` + `flux_potential_lin`.
2. Add the file to `SHARED_SOURCES` in `CMakeLists.txt` (no new library; no new link edges).
3. In `meds_plant_hydraulics.f90`: delete the moved code; add
   `use meds_hydro_curve, only : kirchhoff_edge, capacitance, water_content, plc_retained`;
   drop the now-unused moved params; keep the solver + `rhizosphere_cond` + solver params.
4. Repoint `test/test_plant_hydraulics.f90` lines 20–22 to `use meds_hydro_curve`.
5. Add `derive_parameters` in `shared/config` wrapping the **four** existing calls (`set_allometry`,
   `derive_config`, `derive_pft_rates`, `derive_leaf_params`) + `set_hydro_tables` last; replace the
   two derivation sequences (`meds_config_io.f90:658-661`, `meds_test_support.f90:117-121`) with a
   call to it, leaving the `[derived]` override block (`meds_config_io.f90:663-671`) and
   `validate_config` (`:673`) in place *after* it. Either promote the 9 allometry coefficients to
   `cfg` fields (so the call is cfg-only) or keep `set_allometry` at each site (§7).
6. Build ifx Debug + nvfortran multicore; run the full CTest suite; confirm bit-identical goldens.

## 7. Risks & out-of-scope

- **Risk (low): intra-`shared` module ordering.** `derive_parameters` sits above `meds_hydro_curve`
  and `meds_pft_params` within the same library; CMake's auto-dependency handles `.mod` order, but
  confirm no accidental *use*-cycle inside `shared` (functions must not `use` config).
- **Risk (low): `set_hydro_tables` signature vs storage choice** (§4a) — resolve Option 1 vs 2
  before writing it; it changes whether the tables are module state or `cfg` data.
- **DECISION: add a `cfg%allom` block.** `set_allometry`'s 9 globals are currently `config_io`
  locals / test literals, never stored in `cfg`. Locked choice: promote them to a `cfg%allom`
  derived-type field (`b1Ht, b2Ht, agb_c1, agb_c2, ca_b1, ca_b2, lai_b1, lai_b2, light_ext`) that
  both call sites populate (config_io from TOML at `:638-646`; `build_test_config` from the literals
  at `:117-118`). `derive_parameters(cfg)` then installs via `call set_allometry(cfg%allom%b1Ht, …)`,
  becoming genuinely cfg-only. This makes `cfg` the complete run record (consistent with "config is
  the source of truth") and removes the config_io-locals-never-in-cfg smell. `meds_config` gains
  `use meds_allometry, only : set_allometry` (new, non-cyclic `shared → shared/functions` edge).
- **Must-not-drop: `derive_leaf_params`.** The consolidated routine wraps *four* calls; omitting the
  leaf-param derivation leaves `jmax25/tpu25/rd25` unset and silently breaks bit-identicality. The
  golden-anchor + `pft_parameters.csv` provenance checks are the backstop.
- **Out of scope:** the linear-table *adoption* on the hot path (switching the solver to consult
  `flux_potential_lin` for general kexp), the Padé alternative, and any numerical retuning. This doc
  only relocates code and adds the dormant table machinery + the derivation seam.

## 8. Implementation phasing (scoping correction)

Discovered during implementation planning: **hydraulics traits are not in `cfg`** — `pft_table_t`
has no hydraulics fields, `meds_config_t` has no hydro block, and the per-plant `hydro_params_t` is a
hardcoded stub (`meds_fast_loop.f90:101-107`) filled after config load. So `set_hydro_tables(cfg)`
cannot be written yet (nothing per-PFT to build from). The work therefore splits:

### Phase A — this PR (behaviour-preserving; ships now)
1. **`meds_hydro_curve.f90`** in `shared/functions` (`use meds_kinds` only): move the 12 constitutive
   functions + moving params verbatim; add the POD `hydro_table_t` (fixed `g(0:NTAB)`, `NTAB=512`,
   `r_max=8`, plus stored `kexp`), the `build_hydro_table` **subroutine**, and `flux_potential_lin`.
2. Add to `SHARED_SOURCES` in CMake (no new library/edge).
3. **`meds_plant_hydraulics.f90`**: delete moved code; add
   `use meds_hydro_curve, only : kirchhoff_edge, capacitance, water_content, plc_retained`; keep the
   solver, `rhizosphere_cond`, and the solver params.
4. **`meds_plant_types.f90`**: `use meds_hydro_curve, only : hydro_table_t`; add a dormant
   `vuln_table` field to `hydro_params_t` (Option 2 storage, unpopulated in Phase A).
5. **`test_plant_hydraulics.f90`**: repoint the two `use` lines to `meds_hydro_curve`; add a unit
   test for `build_hydro_table` + `flux_potential_lin` vs the quadrature (proves the dormant machinery
   + the kexp field).
6. **`cfg%allom`** block on `meds_config_t` (§7 decision); config_io and `build_test_config` populate
   it; add `derive_parameters(cfg)` = `set_allometry(cfg%allom%…)` + `derive_config` +
   `derive_pft_rates` + `derive_leaf_params`.
7. Replace the two duplicated sequences (`meds_config_io.f90:658-661`,
   `meds_test_support.f90:117-121`) with `call derive_parameters(cfg)`; leave `[derived]` override
   (`:663-671`) and `validate_config` (`:673`) after it.
8. Build ifx Debug + nvfortran multicore; full ctest; **bit-identical goldens** (e.g. `tc_split`).

### Phase B — DONE (integrated into this PR)
- **`hydraulics_config_t`** (shared, in `meds_config`) holds the PFT-uniform hydraulics parameters
  (13 PV/vulnerability/conductance scalars + `rhizo_cond`) with field defaults = the former stub
  values; `cfg%hydraulics` carries it. **Design choice:** PFT-uniform (matching the current model,
  which has one hydraulics parameterization, not per-PFT), and expressed as a *shared* type — `cfg`
  cannot hold the plant `hydro_params_t` without a `shared → plant` cycle.
- **`[hydraulics]` TOML block** (`meds_config_io`): 14 opt-in *defaulted* reads from the main file
  (each key defaults to the `hydraulics_config_t` value), so a config without the block is
  byte-identical to the old hardcoded stub. Documented block added to `meds_config_main.toml`.
- **`apply_hydraulics_config`** (`meds_column_dynamics`) is the single flatten seam:
  `hydraulics_config_t → hydro_params_t + rhizo_cond`, and it builds the `vuln_table` from
  `wood_kexp` (mirroring the leaf seam that flattens PFT photosynthesis traits). The
  `meds_fast_loop` hardcoded stub is **removed** — `build_fast_context` calls the seam — and the four
  tests that duplicated the stub (`test_column_dynamics/picard_coupling/column_ark/fast_loop`) now
  call it too (5 copies of the stub collapsed to one definition).
- Hot-path adoption of `flux_potential_lin` for `kexp∉{1,2}` + the `tab%kexp == p%wood_kexp`
  (`inv_h > 0`) staleness guard were delivered in Phase B's first cut (§4a/§4b) and are unchanged.
- Verified: ifx Debug 32/32 + nvfortran multicore 32/32 (bit-identical), `[hydraulics]` key names
  match the reader 1:1, and the general-kexp path is exercised end-to-end.

  *Not done (genuine future work):* per-PFT hydraulics (each PFT its own curve + table) — needs
  `hydro_params_t` selected per cohort's PFT in the fast loop, a deeper change than this MVP.
```
