# MEDS Plant Ecophysiology — Integrated Design

**Two things in one plan:** (I) a **behaviour-preserving refactor** of `src/plant/` into a single flat
plant-ecophysiology library plus a self-contained `demography/`; and (II) the **first new ecophysiology
process** landing into that module — **woody (stem) + fine-root maintenance respiration**. Part I is pure
reorganization (no science change); Part II is the feature. Do Part I first, green, then Part II.

Supersedes and merges the former `MEDS_PLANT_REFACTOR_DESIGN.md` + `MEDS_WOOD_ROOT_DESIGN.md`.

---

## 0. Overview & goals

MEDS's `src/plant/` currently holds four separate libraries behind per-domain subfolders (`leaf/`,
`hydraulics/`, `phenology/`, `vitality/`) plus `meds_allometry.f90`. Two structural cuts fix what the code
*is*:

1. **`src/plant/` → ONE flat library `libmeds_plant` = plant ECOPHYSIOLOGY only** (leaf gas exchange,
   respiration, hydraulics, phenology). These are tightly integrated (the leaf solver couples
   photosynthesis + stomata + Rd; hydraulics couples through `psi_leaf`; respiration shares the
   temperature response), individually small, and there is no reason to compile them apart. Every file is
   already `meds_<domain>_*`-prefixed, so flattening is pure `git mv` — **no renames, no collisions**.

2. **The empirical vital rates move OUT of `plant/` into `demography/`.** `growth`/`mortality`/
   `recruitment`/`vital_rates` are **purely empirical demographic functions** (structure-only relations in
   dbh / wood density / overtopping LAI), not ecophysiology. Moving them makes `demography/` a
   **self-contained demography-only model** — it produces *and* applies its own rates and runs the full
   spin-up with no plant-ecophysiology dependency.

The pay-off of cut #2 on cut #1: with vitality gone, **nothing left in `plant/` touches the cohort state
(`site_t`)**, so `libmeds_plant` links **`meds_shared` only** — standalone, orthogonal to state,
unit-testable (as `leaf/` is today), and the optional Python shared library does **not** drag the
demographic SoA in.

### The two-halves architecture

```
  EMPIRICAL DEMOGRAPHY  ──(growth[], mortality[], recruitment[] data arrays)──►  demographic engine
  (self-contained: demography/ owns the phenomenological rates AND applies them)

  MECHANISTIC ECOPHYSIOLOGY  (plant/: leaf gas exchange, respiration, hydraulics, phenology)
  (standalone kernels; NOT yet wired; will later feed the SAME update_demography rate seam,
   replacing the empirical rates with carbon/water-balance-derived ones)
```

The `update_demography(site, growth, mortality, recruitment, …)` **data-array seam is unchanged.** Today
demography's own empirical provider fills those arrays; tomorrow a plant-ecophysiology provider fills them
— no engine change, no new dependency edge, both halves independently testable.

---

# PART I — The refactored architecture (behaviour-preserving)

## 1. Target source layout

```
src/
  shared/
    meds_kinds.f90  meds_constants.f90  meds_pft_params.f90  meds_config.f90  meds_time.f90
    meds_temp_response.f90        ← MOVED from plant/leaf/meds_leaf_temp_response.f90 (module renamed)
  allometry/
    meds_allometry.f90            ← MOVED from plant/ (own library; a shared foundation below state — §3)
  state/
    meds_demography_types.f90
  plant/                          ← FLAT; ONE library libmeds_plant; ECOPHYSIOLOGY only; links meds_shared ONLY
    meds_plant_types.f90          ← ALL derived types (leaf / hydro / pheno / respiration) in ONE module (§8a)
    meds_plant_interface.f90      ← THE public seams (one door): leaf_gas_exchange, plant_water_flux, update_phenology
    meds_leaf_gas_exchange.f90    ← leaf COMPUTE: photosynthesis + stomata + Ci solver (one module)
    meds_plant_hydraulics.f90     ← hydraulics COMPUTE: conductance + pressure-volume + matrix-exp solver (one module)
    meds_phenology.f90            ← phenology COMPUTE: the cue engine
    meds_plant_respiration.f90    ← NEW, Part II (stem + fine-root maint. resp + growth_respiration; leaf Rd in solver)
    meds_plant_capi.f90           ← REMOVE_ITEM'd from the static lib; compiled into the pylib only
  demography/                     ← now a SELF-CONTAINED demography-only model
    meds_demography_interface.f90  meds_demography_dynamics.f90
    meds_demography_structure.f90  meds_demography_diagnostics.f90
    meds_growth.f90  meds_mortality.f90  meds_recruitment.f90  meds_vital_rates.f90   ← MOVED from plant/vitality/
  biophys/  biogeochem/  io/  driver/  init/       (unchanged by this refactor)
```

Notes: `src/biophys/` (canopy radiation, a sibling stateless-kernel library) is **out of scope** — it may
adopt the same flattening later, independently. Moved vitality files keep their module names
(`meds_growth`, `meds_mortality`, `meds_recruitment`, `meds_phenomenological_vital_rates`); files gain the
`meds_` prefix to match demography.

## 2. Libraries — before → after

| Library (before) | Fate | Library (after) | Links |
|---|---|---|---|
| `meds_allometry` (src/plant/) | move dir | `meds_allometry` (src/allometry/) | `meds_shared` |
| `meds_leaf_physiology` | **dissolve →** | `meds_plant` | `meds_shared` **only** |
| `meds_hydraulics` | **dissolve →** | `meds_plant` | — |
| `meds_phenology` | **dissolve →** | `meds_plant` | — |
| `meds_vitality` | **dissolve →** | `meds_demography` | (already links state+allometry) |
| `meds_leaf_c` (pylib) | rename | `meds_plant_c` | `meds_plant` |
| `meds_state`,`meds_demography`,`meds_aux`,`meds_io`,`meds_config_io` | keep (link-list edits) | same | — |

Four plant libraries collapse to **one** (`meds_plant`); the empirical rates fold into `meds_demography`;
`meds_allometry` relocates. Blast radius is tiny: the dissolved physiology libraries' only consumers were
their own tests + the capi; `meds_vitality`'s only consumer was `meds_aux`.

## 3. The dependency DAG — allometry is a foundation, not ecophysiology

```
shared
  ├── allometry ─────────────┐
  │                          ▼
  │                        state ──── demography ──── aux ──── main
  └── plant (ecophysiology) ┄┄┄┄┄┄┄ (standalone; links shared only; NOT yet wired to demography)
```

- `plant` links **`shared` only** and is orthogonal to `state` — consumed today only by its tests and the
  Python pylib. Allometry-derived geometry (`dbh`, `height`, `wai`, `sap_area`, `broot`) reaches the
  ecophysiology kernels through their `*_env_t` structs, computed by the caller, so the kernels never
  `use meds_allometry` (`hydro_env_t` already carries `sap_area/bleaf/bsap` this way).
- **Why allometry can't stay in the flat `plant/` library.** It *can*, acyclically (plant no longer
  depends on state, so there is no cycle) — but it *must not*, because **`state` depends on allometry**
  (`set_cohort_size` → `dbh_to_height`; `agb_to_dbh` in fusion). If allometry were compiled into
  `libmeds_plant`, then `state` — and transitively the whole demographic core — would depend on the
  *entire* ecophysiology library (leaf solver included). That breaks `cmake --build --target
  meds_demography` standalone **and directly contradicts the "self-contained demography" goal** of cut #2.
  Allometry is a *shared structural-geometry foundation* used by **both** the demographic core and
  ecophysiology, so it belongs in its own library below `state`. It is already a separate library today;
  the refactor just moves the file to `src/allometry/` to match its DAG position.
  *(Alternative: fold `allometry` into `meds_shared` — one fewer library, but bloats the foundation.
  Recommended: keep it its own library — it has real domain logic and its own `test_allometry`.)*

## 4. CMake — the rewrite

```cmake
add_library(meds_shared STATIC ${SHARED_SOURCES})            # SHARED_SOURCES gains meds_temp_response.f90

add_library(meds_allometry STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/allometry/meds_allometry.f90)
target_link_libraries(meds_allometry PUBLIC meds_shared)
meds_fortran_flags(meds_allometry)

add_library(meds_state STATIC ${CMAKE_CURRENT_SOURCE_DIR}/src/state/meds_demography_types.f90)
target_link_libraries(meds_state PUBLIC meds_shared meds_allometry)

# ONE flat plant-ecophysiology library (links shared ONLY)
file(GLOB PLANT_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/plant/*.f90)
file(GLOB PLANT_CAPI    CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/plant/*_capi.f90)
list(REMOVE_ITEM PLANT_SOURCES ${PLANT_CAPI})               # capi -> pylib only, never the static core
add_library(meds_plant STATIC ${PLANT_SOURCES})
target_link_libraries(meds_plant PUBLIC meds_shared)
meds_fortran_flags(meds_plant)

# demography absorbs the empirical vital rates (its GLOB already covers src/demography/*.f90)
add_library(meds_demography STATIC ${DEMOG_SOURCES})
target_link_libraries(meds_demography PUBLIC meds_state meds_allometry)

# physiology tests all link the one library
foreach(t leaf_physiology plant_hydraulics plant_phenology plant_respiration)
  add_executable(test_${t} test/test_${t}.f90)
  target_link_libraries(test_${t} PRIVATE meds_plant)
  meds_fortran_flags(test_${t})
endforeach()

# meds_aux: drop meds_vitality (folded into demography)
#   target_link_libraries(meds_aux PUBLIC meds_demography meds_config_io)   # was: + meds_vitality

option(MEDS_BUILD_PYLIB "Build the plant C-API shared library for Python (ctypes)" OFF)
if(MEDS_BUILD_PYLIB)
  set_target_properties(meds_shared meds_plant PROPERTIES POSITION_INDEPENDENT_CODE ON)
  add_library(meds_plant_c SHARED ${PLANT_CAPI})            # GLOB: leaf capi today; +respiration/hydro later
  target_link_libraries(meds_plant_c PRIVATE meds_plant)
  meds_fortran_flags(meds_plant_c)
endif()
```

Because `libmeds_plant` links `meds_shared` only, the pylib is PIC over just `{shared, plant}` — the cohort
SoA is **not** pulled into the leaf-curve shared library.

## 5. C-API — symbol-stable, GLOB-extensible

`meds_leaf_capi.f90` → `src/plant/meds_plant_capi.f90` (module `meds_leaf_capi` → `meds_plant_capi`):

- **Every exported C symbol and C-mirror struct is preserved** — `meds_leaf_solve`, `meds_assim_demand_c3`,
  `meds_electron_transport_j`, `meds_peaked_arrhenius`, `meds_arrhenius`, and the `leaf_env_c` /
  `leaf_params_c` / `leaf_flux_c` / `leaf_c3_demand_c` structs. `python/meds/leaf/_ffi.py` bindings do not
  change. The only internal edit is `use meds_leaf_temp_response` → `use meds_temp_response`.
- **`*_capi.f90` is GLOBbed into `libmeds_plant_c`**, so each domain adds its own capi file into the same
  shared library as it gains a Python need (e.g. `meds_respiration_capi.f90` → `meds_stem_resp`,
  `meds_fine_root_resp`; or a `plant_hydrodynamics` shim).
- **The Python surface becomes process-oriented** (designed in the package, decoupled from Fortran layout):
  ```python
  from meds import leaf                 # meds.leaf.gas_exchange(env, pft)  -> A-Ci / A-PAR / A-T / gs-VPD
  from meds import respiration          # meds.respiration.stem(...) / .fine_root(...)   (future capi)
  # plant_hydrodynamics(env) integrating soil -> root -> stem -> leaf is one more capi shim + wrapper.
  ```

## 6. Interface / compute split — one `meds_plant_interface`

The flat module is organized as **one public-seam module over per-domain compute modules** (a readability
reorg on top of the flatten; each domain's small modules collapse into one compute file, and all the seams
live together):

- **`meds_plant_interface`** — THE single public door. Hosts the thin seams `leaf_gas_exchange(env, cfg,
  ipft, flux)` (flattens `cfg%pft` → params), `plant_water_flux(...)`, and `update_phenology(...)`, and
  re-exports the public types. Callers `use meds_plant_interface` and nothing else. (Future:
  `plant_hydrodynamics` integrating soil → root → stem → leaf, and the respiration seams, land here too.)
- **`meds_leaf_gas_exchange`** — leaf compute: FvCB C3 + Collatz C4 demand + electron transport + the
  Leuning/Medlyn/Katul stomata + the bracketed Ci solver, one module. Its raw kernels
  (`solve_leaf_gas_exchange`, `assim_demand_c3`, `electron_transport_j`) are also called directly by
  `meds_plant_capi` — so leaf compute stays a separately-`use`-able module, NOT folded into the interface.
- **`meds_plant_hydraulics`** — hydraulics compute (pressure-volume + Kirchhoff conductance + matrix-exp
  solver), one module. (Symmetric with leaf: compute here, seam in the interface — the interface is NOT
  bloated with hydraulics math.)
- **`meds_phenology`** — phenology compute (the cue engine).

The process-oriented Python functions wrap the interface seams. Leaf Rd is computed inside the leaf solver
(the Ci solve needs it) and returned as `leaf_flux_t%rd`; `meds_plant_respiration` (Part II) adds the
non-leaf respiration compute, with its seams exposed through `meds_plant_interface`.

## 7. What moves — migration map

| From | To | Change |
|---|---|---|
| `plant/leaf/meds_leaf_temp_response.f90` | `shared/meds_temp_response.f90` | rename module → `meds_temp_response`; repoint every `use` |
| `plant/meds_allometry.f90` | `allometry/meds_allometry.f90` | own library; no code change |
| `plant/leaf/*.f90` (minus temp_response, capi) | `plant/*.f90` | `git mv`; no rename (already `meds_leaf_*`) |
| `plant/hydraulics/*.f90` | `plant/*.f90` | `git mv`; no rename (already `meds_hydro_*`) |
| `plant/phenology/*.f90` | `plant/*.f90` | `git mv`; no rename (already `meds_pheno_*`) |
| `plant/vitality/{growth,mortality,recruitment,vital_rates}.f90` | `demography/meds_*.f90` | `git mv` + `meds_` prefix; module names unchanged |
| `plant/leaf/meds_leaf_capi.f90` | `plant/meds_plant_capi.f90` | module → `meds_plant_capi`; C symbols unchanged; `use meds_temp_response` |

Delete the emptied `plant/{leaf,hydraulics,phenology,vitality}/` dirs; fold their READMEs into one
`src/plant/README.md`.

## 7a. Consolidated types — one `meds_plant_types` module

All the per-domain `*_types` modules merge into a **single `meds_plant_types.f90`** (one module
`meds_plant_types`), organized in clearly-commented sections:

```fortran
module meds_plant_types
   use meds_kinds, only : wp, ik
   implicit none
   ! ==== leaf ====   PATH_C3/C4, LIM_*, leaf_env_t, leaf_photo_params_t, leaf_flux_t
   ! ==== hydraulics = hydro_env_t, hydro_flux_t, node topology params, ...
   ! ==== phenology == CUE_*, PHEN_*, pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t
   ! ==== respiration = wood_env_t, wood_params_t, wood_flux_t, root_env_t, root_params_t, root_flux_t
end module
```

Every procedure module (`meds_leaf_gas_exchange`, `meds_plant_hydraulics`, `meds_phenology`,
`meds_plant_respiration`, `meds_plant_interface`, and `meds_plant_capi`) does
`use meds_plant_types, only: …`.

**Consequences of the merge (the explicit question):**
- **Computational: NONE.** Module organization is compile-time only; a derived type generates identical
  layout/code regardless of which module defines it. Zero runtime effect.
- **Logical: none beyond mechanics.** The types are mutually independent records (no leaf↔hydro↔pheno type
  references, so no ordering or circularity issue), and every name is distinct (`leaf_env_t`/`hydro_env_t`/
  `pheno_env_t`/`wood_env_t`/`root_env_t`; `PATH_*`/`LIM_*`/`CUE_*`/`PHEN_*`) — no collisions. The only edits
  are mechanical: `use meds_leaf_types` → `use meds_plant_types` at each call site (kept minimal by
  `only:`). The C-mirror `bind(c)` types stay in `meds_plant_capi` (unchanged).
- **Recompile granularity:** editing any type recompiles every type-consumer. But `libmeds_plant` already
  compiles as a whole, so that is the existing granularity — no practical loss.
- This consciously overrides CLAUDE.md's "one module per file" *for the types module only* (a single
  cohesive data module for a single-library domain); the procedure modules keep one-module-per-file.

*(Lower-churn alternative, if preferred: keep `meds_leaf_types`/`meds_hydro_types`/`meds_pheno_types` as
separate modules inside the one `meds_plant_types.f90` file — then no `use` lines change. Recommended is
the single module above, matching "a single plant_types".)*

## 8. Migration order + acceptance (each step verified EXACT before the next)

1. **Promote temp response** → `src/shared/meds_temp_response`; repoint `use`s; delete old file.
2. **Relocate allometry** → `src/allometry/`.
3. **Fold vitality into demography** (`git mv` + `meds_` prefix; drop `meds_vitality` from CMake + `aux`).
4. **Flatten plant** (`git mv` leaf/hydraulics/phenology up; delete subfolders; one `meds_plant` lib; fix
   test links).
5. **Consolidate types** → merge `meds_{leaf,hydro,pheno}_types` into one `meds_plant_types` (§7a); repoint
   `use` lines. Verify: all physiology ctests + leaf example EXACT.
6. **Rename the capi** → `meds_plant_capi`; capi GLOB → `libmeds_plant_c`; bump lib name in `_ffi.py`.

**Acceptance (behaviour-preserving):**
- ifx `-stand f18 -check all -fpe0` and nvfortran multicore `-O3`: full ctest suite green (a green ifx run
  is not sufficient — CLAUDE.md issue #7; the array-temporary trap applies to any new glue).
- 250-yr spin-up **bit-reproduces** LAI = 7.3783 / AGB = 16.535 (336 cohorts, 12 patches).
- `-DMEDS_BUILD_PYLIB=ON` builds `libmeds_plant_c`; after the one lib-name bump in `_ffi.py`,
  `examples/example_leaf_physiology/reproduce_slot2017.py` reproduces `slot2017.png`; `pytest
  python/tests/test_leaf.py` green.
- `demography/` builds and runs standalone (no `meds_plant` link) — the self-containment check.
- `plant/` builds and tests standalone (`cmake --build … --target meds_plant`) — links `meds_shared` only.

Do steps 1–6 as one refactor PR; the Part-II feature is a separate PR.

---

# PART II — First ecophysiology feature: wood + root maintenance respiration

> **IMPLEMENTED 2026-07-04 (PR after #12, branch `feature/plant-respiration`).** First cut follows the
> hydraulics/phenology precedent: **module-local params** (`wood_params_t`/`root_params_t` in
> `meds_plant_types`, built by the caller/test) — **not** `pft_table_t`/config — so it is non-invasive and
> standalone-tested. **All maintenance factors are 25 °C-referenced** (MEDS's single model-wide
> reference); the kernels just consume `stem_resp_factor25`/`root_resp_factor25`. ED2 references stem/root
> respiration at **15 °C**, so seeding a MEDS default from an ED2/Chambers number is a **one-time
> conversion done at PARAMETER INIT** (where PFT params are chosen) — **not** a runtime helper and **not**
> in `meds_pft_params` (`init` sits high in the DAG and can use `meds_temp_response` with no cycle). §13/§14
> below (the `pft_table_t` + `[respiration]` config path) are the **deferred config-driven seam**, exactly
> as for hydraulics/phenology. Verified: ifx `-check all` 14/14 ctests, nvfortran multicore 14/14;
> `test_plant_respiration` = 10/10 checks.

Lands as a single **compute** module `src/plant/meds_plant_respiration.f90` (the 7th plant file), with its
public seams exposed through `meds_plant_interface` — symmetric with leaf/hydraulics/phenology (compute in
its own module, seam in the interface; §6). Its types (`wood_env_t`/`root_env_t`/…) go in
`meds_plant_types`. Leaf dark respiration (Rd) stays inside the leaf gas-exchange solver — the Ci root-find
needs `A_net = A_gross − Rd`, and it is already exposed as `leaf_flux_t%rd`. So this module covers the
**non-leaf** maintenance respiration (stem + fine root) — exactly FATES's `maintresp_nonleaf` grouping.

## 9. Scope & framing — non-leaf maintenance respiration + the growth-respiration formula

MEDS's carbon flow will be `A = GPP − Rleaf − Rm_wood − Rm_root` (after-maintenance carbon, ED2's
"metabolic NPP"), then `Rg = growth_resp_factor·max(0, A)`, and `NPP = A − Rg`.

The module hosts **three respiration terms** (leaf Rd is the exception — it stays in the leaf gas-exchange
solver because the Ci root-find needs `A_net = A_gross − Rd`, and it is already exposed as
`leaf_flux_t%rd`):

- **Stem maintenance respiration `Rm_wood`** and **fine-root maintenance respiration `Rm_root`** —
  instantaneous, temperature-driven, per-plant fluxes (§11–§12).
- **Growth respiration `Rg`** — a **pure formula, a fraction of the input after-maintenance carbon**
  (ED2 `growth_resp = growth_resp_factor·(GPP − Rleaf − Rroot − Rstem)`, clamped ≥ 0). It takes `A` as an
  **argument** — it does *not* compute GPP or maintenance, and it owns no carbon pool — so it fits the
  stateless module cleanly (§12, §13). The future allocation/NPP module computes `A` (summing GPP − leaf Rd
  − wood Rm − root Rm), calls `growth_respiration(A, factor)`, and applies `NPP = A − Rg` to allocation.

So the split is: **the respiration module owns the respiration *formulas* (including `Rg`'s fraction rule);
the future allocation module owns the carbon *budget* (computing `A`, applying `NPP`).** This keeps all
autotrophic-respiration arithmetic except leaf Rd in one place, without giving the module a carbon pool.

### What the reference survey found (drives the decisions)

- **ED2** (the model we port). Stem respiration (`stem_resp_driv.f90`) is **surface-area based** (Chambers
  et al. 2004): `stem_area = (dbh·height·nplant + WAI)·π/agf_bs`, size-dependent baseline
  `10^(log10(SRF)+size_scaler·dbh)`, Arrhenius/Collatz response ÷ low/high-T inhibition, driven by explicit
  **wood temperature**; grasses = 0. Fine-root respiration (`soil_respiration.f90:root_resp_norm`) is
  **fine-root-biomass based** (`factor·broot·nplant`), integrated over soil layers by **soil temperature**.
  Both factors are referenced to **15 °C**, and their T-params **default to the leaf Rd params**.
- **FATES / CLM5 / ELM** use per-tissue-**nitrogen** maintenance (`Rm = N·2.525e-6 gC gN⁻¹ s⁻¹ ·
  Q10^((T−20)/10)`, Ryan 1991; only *live/sapwood* tissue respires; `grperc≈0.11`). **MEDS cannot use
  this** — no tissue N, no biomass pools — and FATES's *own default sapwood N:C = 1e-8 zeros live-wood
  respiration*. So **this is an ED2 port, not a FATES/CLM port**; FATES/CLM inform only the deferred
  growth-fraction and acclimation term.
- **Literature.** Both fluxes are large and badly constrained; form matters more than tuning. The #1
  improvement over ED2 is **thermal acclimation** (stems: Ziegler 2024 *Science*; roots: Reich 2016
  *Nature* — ~80% of the fixed-Q10 warming response disappears). Measured stem CO₂ **efflux ≠ local
  respiration** (xylem-transported CO₂; Trumbore 2013) — the Chambers coefficient is an apparent-efflux
  parameter (a code comment, not a mechanism). Wood-density ρ is a legitimate inverse proxy for
  sapwood-parenchyma N — a MEDS-native modulation of the stem coefficient.

## 10. Design decisions

1. **ED2 port, maintenance-only, one stateless module.** Wood = ED2 Chambers surface-area form; root = ED2
   per-`broot` form. Lives in the flat `plant/` (links `meds_shared` only), driven by env structs, not
   wired to demography.
2. **Reuse the promoted `meds_temp_response` (Part I §1).** Wood/root use `temp_response(TRESP_PEAKED, …,
   ea_rd, hd_rd, ds_rd, …)` — **peaked** (bounded at high T) and consistent with leaf Rd. MEDS
   **deliberately drops ED2's `tlow/thigh` inhibition denominators** (ED2 itself just copies them from leaf
   Rd; "harmless but unjustified") — a documented simplification.
3. **One reference temperature (25 °C); convert ED2's 15 °C factors once, at ingest** (`derive_resp_params`,
   next to `derive_leaf_params`). `pft_table_t` stores only 25 °C-referenced factors. *(The single most
   likely bug otherwise — §16 Defect ①.)*
4. **Per-plant flux unit** (`[µmol CO₂ plant⁻¹ s⁻¹]`) — the `×stem_area`/`×broot` bookkeeping is inside the
   kernel — so leaf Rd + wood Rm + root Rm compose with the per-plant AGB pool.
5. **Ship stem first; root has a hard prerequisite.** Stem runs off `dbh/height/nplant` (WAI=0 default,
   `wood_temp` an env input). Root is 0 without `broot`, so the `broot = (leaf_area/SLA)·q` allometry is a
   stated prerequisite (reused by hydraulics), and root uses a single root-weighted mean `soil_temp` (the
   layered integral stays *outside* the kernel).
6. **Reserve the acclimation seam now, off by default** — a `t_acclim` field in the env struct + an
   `acclimation` config enum, so enabling Atkin/CLM5 acclimation later is a parameter flip, not a re-cut.
7. **N-ready internally** — the root coefficient is written `base_rate_per_N × n_conc` (`n_conc` a per-PFT
   constant today), so a future stoichiometry state drops in with no interface change.
8. **`is_woody` guard** — MEDS has only trees (wood-density axis, no grass flag); the guard gives grasses 0
   the day a grass PFT is added.

## 11. The module — `meds_plant_respiration.f90`

Public procedures (stateless; the two maintenance kernels are pure functions of `(env, params)`; the
growth-respiration term is a pure `elemental` function of scalars — the reserved `t_acclim` is read-only):

```fortran
subroutine stem_maintenance_respiration(env, params, out)      ! wood
   type(wood_env_t),    intent(in)  :: env      ! wood_temp (+ reserved t_acclim), dbh, height, wai, nplant
   type(wood_params_t), intent(in)  :: params   ! flat per-PFT
   type(wood_flux_t),   intent(out) :: out       ! stem_resp [umol CO2 / plant / s]
end subroutine

subroutine fine_root_maintenance_respiration(env, params, out) ! root
   type(root_env_t),    intent(in)  :: env      ! soil_temp (+ reserved t_acclim), broot
   type(root_params_t), intent(in)  :: params
   type(root_flux_t),   intent(out) :: out       ! root_resp [umol CO2 / plant / s]
end subroutine

! Growth (construction) respiration as a fraction of the input after-maintenance carbon (ED2 form).
! npp_in = GPP - Rleaf - Rm_wood - Rm_root (ED2's "metabolic NPP"), supplied by the caller; NO carbon
! pool is owned here. rg carries npp_in's units (per plant). true NPP = npp_in - rg.
elemental pure function growth_respiration(npp_in, growth_resp_factor) result(rg)
   real(wp), intent(in) :: npp_in              ! [carbon flux / plant] after-maintenance carbon (A)
   real(wp), intent(in) :: growth_resp_factor  ! [--] per-PFT construction-cost fraction
   real(wp)             :: rg                   ! [same units as npp_in] growth respiration
   rg = growth_resp_factor * max(0.0_wp, npp_in)
end function growth_respiration
```

Types (declared in the consolidated `meds_plant_types`, §7a):

```fortran
type :: wood_env_t
   real(wp) :: wood_temp, dbh, height, wai, nplant, t_acclim   ! t_acclim RESERVED (unused v1)
end type
type :: wood_params_t
   logical  :: is_woody
   real(wp) :: stem_resp_factor25, stem_resp_size_scaler, agf_bs, ea, hd, ds
end type
type :: wood_flux_t ; real(wp) :: stem_resp ; end type   ! [umol CO2 / plant / s]

type :: root_env_t   ; real(wp) :: soil_temp, broot, t_acclim ; end type
type :: root_params_t; real(wp) :: root_resp_factor25, ea, hd, ds ; end type   ! factor25 = base_rate_per_N*n_conc
type :: root_flux_t  ; real(wp) :: root_resp ; end type  ! [umol CO2 / plant / s]
```

(`growth_respiration` is a bare scalar function — no env/flux struct — because it needs only `npp_in` and
the per-PFT factor; it is unit-agnostic, returning whatever flux unit the caller passes in.)

## 12. Governing equations

`s(T) ≡ temp_response(TRESP_PEAKED, 1, ea, hd, ds, T)` is the dimensionless peaked-Arrhenius scale (=1 at
298.15 K), from the shared `meds_temp_response`.

**Stem (per plant):**
```
if (.not. is_woody) stem_resp = 0
srf25(dbh)         = stem_resp_factor25 * 10**( stem_resp_size_scaler * dbh )   ! scaler=0 => flat (default)
stem_area_perplant = ( pi*(dbh*1e-2)*height + pi*wai/nplant ) / agf_bs          ! [m2 stem / plant]
stem_resp          = srf25(dbh) * s(wood_temp) * stem_area_perplant             ! [umol CO2 / plant / s]
```
The size scaler ships **off** (ED2's positive DBH exponent is noisy and fights the shrinking-sapwood
argument); the MEDS-native alternative is ρ-modulation of `stem_resp_factor25` in `derive_resp_params`.

**Fine root (per plant), single effective soil temperature:**
```
root_resp = root_resp_factor25 * s(soil_temp) * broot                           ! [umol CO2 / plant / s]
```
`soil_temp` is a root-biomass-weighted mean over an assumed exponential root profile (the weighting lives
in the caller; the kernel sees one scalar, so a later per-layer sum never touches the seam). `broot=0 ⇒ 0`.

**Growth (construction) respiration — a fraction of input after-maintenance carbon (ED2 form):**
```
rg = growth_resp_factor * max(0, npp_in)     ! npp_in = GPP - Rleaf - Rm_wood - Rm_root ; NPP = npp_in - rg
```
Temperature-independent (a construction-cost overhead, McCree/Thornley), so no `s(T)` factor. `npp_in` and
`rg` share the caller's carbon-flux unit (per plant). `growth_resp_factor` is a per-PFT trait (§13). This
is the *formula only*; the caller supplies `npp_in` and applies `NPP`.

## 13. PFT traits + the 15 °C → 25 °C ingest conversion  *(DEFERRED config-driven path — the shipped first cut uses module-local, 25 °C-based params; any 15→25 °C conversion is a one-time parameter-init step, not a runtime helper. See the Part II banner)*

New `pft_table_t` columns: `is_woody`, `stem_resp_factor_15`, `stem_resp_size_scaler` (=0),
`agf_bs`, `root_resp_factor_15`, `root_n_conc` (optional), **`growth_resp_factor`** (for §12's `Rg`), and
DERIVED `stem_resp_factor25` / `root_resp_factor25`. Conversion in a new `derive_resp_params` (beside
`derive_leaf_params`):

```fortran
subroutine derive_resp_params(pft)
   type(pft_table_t), intent(inout) :: pft
   real(wp), parameter :: t_ref_ed2 = 288.15_wp                 ! 15 degC (ED2 stem/root reference)
   real(wp), allocatable :: s15(:)
   ! rate whose 25 degC value is 1, evaluated at 15 degC  =>  factor25 = factor15 / s15
   s15 = peaked_arrhenius_scale(1.0_wp, pft%ea_rd, pft%hd_rd, pft%ds_rd, t_ref_ed2)
   pft%stem_resp_factor25 = pft%stem_resp_factor_15 / s15
   pft%root_resp_factor25 = pft%root_resp_factor_15 / s15
end subroutine
```
(Bind `s15` to a named array — nvfortran array-temporary trap.) ED2 defaults to seed the PFT TOML:
`stem_resp_factor_15 = 10^(-0.672-0.19)/2.2 ≈ 0.0587` µmol/m²stem/s@15 °C; `stem_resp_size_scaler = 0.0041`
(ship 0); `agf_bs ≈ 0.7`; `root_resp_factor_15 = 0.2455` (tropical) / `0.528` (temperate) µmol/kgC/s@15 °C;
`ea_rd/hd_rd/ds_rd = 46390/200000/490`; `growth_resp_factor ≈ 0.33` (ED2 tropical broadleaf; ED2 conifer
0.45, grass ⅓; FATES/CLM `grperc = 0.11` — tunable per PFT).

## 14. Config  *(DEFERRED — first cut has no config block; params are module-local, like hydraulics/phenology)*

```toml
[respiration]
# Non-leaf (stem + fine-root) MAINTENANCE respiration (ED2 forms). Not yet coupled to demographic carbon.
temp_response_form = "peaked"     # "arrhenius" | "peaked" (shared with leaf Rd; peaked default)
stem_size_scaling  = false        # ED2's 10^(scaler*dbh) size effect (off => flat baseline)
acclimation        = "none"       # "none" | "linear" (RESERVED: Atkin/CLM5 running-mean-T; unused v1)
```

## 15. Units

Both `*_flux_t` fields are `[µmol CO₂ plant⁻¹ s⁻¹]` (× nplant → ED2's per-ground units); stated in the
field comment (as `leaf_flux_t%rd` / `hydro_flux_t` do). The future carbon module forms `NPP = GPP − Rleaf
− Rm_wood − Rm_root` per plant, then `Rg = grperc·max(0,NPP)`, converting to kgC at its own boundary.

## 16. Defects this design pre-empts (from the adversarial review)

| # | Sev | Defect | Guard |
|---|-----|--------|-------|
| ① | CRIT | ED2 15 °C factor into the 25 °C-locked `temp_response` → **~1.8–2× silent low bias** | §13 ingest conversion; only 25 °C factors reach kernels |
| ② | HIGH | "Reuse leaf machinery" silently changes the T-form (Medlyn ≠ ED2 tlow/thigh; plain Arrhenius unbounded) | §10.2: peaked (bounded), inhibition dropped on purpose, documented |
| ③ | HIGH | Root kernel is 0 without `broot`; stem/root not symmetric | §10.5: `broot` allometry prerequisite; stem ships first |
| ④ | HIGH | Jumping to FATES/CLM N-basis (no N; FATES default zeros live-wood MR) | ED2 forms; N-ready internal factor only (§10.7) |
| ⑤ | MED | Leaf/wood/root fluxes won't compose (per-leaf-area vs per-ground vs per-plant) | §15 per-plant unit |
| ⑥ | MED | Acclimation retrofit breaks the sealed seam | §10.6 `t_acclim` + config enum reserved now |
| ⑦ | LOW | ED2 DBH scaler is a noisy fudge | shipped off; ρ-modulation reserved |
| ⑧ | LOW | Three WAI notions (biophys 0.1·LAI, ED2 allometry, new) | single `size→wai` allometry (§17); MVP `wai=0` |

---

# PART III — Cross-cutting

## 17. Reserved follow-ups (documented; not in the first cuts)

Each deferred with an explicit seam:
- **Growth-respiration *carbon budget*** — the `Rg` *formula* now ships in the module (§11–§12); what is
  deferred is the *caller* that computes `A = GPP − Rleaf − Rm_wood − Rm_root` (summing all maintenance
  terms) and applies `NPP = A − Rg` to allocation — the future allocation/NPP module.
- **Storage / turnover respiration** — needs `bstorage`.
- **Tissue turnover → litter** (`root/bark_turnover_rate` → soil-C) — the natural next wood/root process,
  but it feeds biogeochemistry, not `plant/`.
- **Thermal acclimation** (Atkin 2015 / CLM5 `rootstem_acc`) — highest science upgrade; `t_acclim` seam
  reserved (§10.6); the running-mean-T *updater* + its state land when a forcing source exists (template:
  phenology's caller-owned `pheno_state_t`).
- **Per-N maintenance basis** — the root factor is already `base_rate_per_N × n_conc`; add a root-N state.
- **`broot` allometry** — `broot = (leaf_area/SLA)·q` in `meds_allometry` (reused by hydraulics).
- **WAI allometry** — one `size→wai` function in `meds_allometry` (`wai = nplant·b1WAI·size^b2WAI`), shared
  by `biophys` (which currently fabricates `wai ≈ 0.1·LAI`) and respiration; MVP uses `wai=0`.
- **Layered soil respiration** — ED2's per-layer `Σ root_resp_norm(T_soil(k))·Δz(k)`, when a soil thermal
  column + `krdepth` exist; kernel already takes one effective `soil_temp`.
- **`plant_hydrodynamics` seam** — integrate soil → root → stem → leaf hydraulics (the whole-plant
  water-transport seam) into the flat module.
- **Wiring ecophysiology into demography** — a mechanistic rate provider feeding the `update_demography`
  data-array seam, replacing the empirical vital rates.

## 18. Validation & milestones

- **Part I (refactor)** — §8 acceptance: EXACT reproduction (spin-up bit-identical; leaf example + pytest;
  ifx + nvfortran) plus the two standalone-build checks (`demography` without plant; `plant` shared-only).
- **Part II (respiration)** — `test_plant_respiration` (links `meds_plant`): `is_woody=.false. ⇒ 0`;
  `s(298.15 K)=1` identity at 25 °C; 15→25 °C conversion round-trip; hand-computed reference for known
  `(dbh,height,nplant,wai,T)`; `wai=0` vs `wai>0`; `broot=0 ⇒ 0`; monotone rise with `soil_temp` to the
  peaked optimum; **growth respiration** `growth_respiration(npp,f) = f·npp` for `npp>0`, `= 0` for
  `npp≤0`. **Oracle = numerical agreement with a reference implementation of the chosen form** (no
  ecosystem consumer yet, so "it runs" is not validation). Green on ifx (`-check all -fpe0`) + nvfortran
  multicore.
- **Milestone order:** Part I steps 1–5 (one PR) → `broot` allometry → `meds_plant_respiration` (stem, then
  root) (a second PR).

## 19. How this differs from ED2 and FATES/CLM (respiration)

| Aspect | ED2 | FATES / CLM5 | **MEDS (this design)** |
|--------|-----|--------------|------------------------|
| Stem basis | surface area (cylinder + WAI), Chambers | sapwood **N** × base rate | **surface area** (ED2), WAI optional |
| Root basis | per **kgC** fine root, layer-integrated | per **gN** fine root, layer-integrated | **per kgC** (ED2), single effective soil T; N-ready internally |
| Reference T | **15 °C** | **20 °C** | **25 °C** (one reference; converted at ingest) |
| T-response | Arrhenius/Collatz **÷ tlow·thigh** | plain **Q10 = 1.5** | **peaked Arrhenius** (bounded), inhibition dropped |
| Growth resp | `grf·(P−Rleaf−Rroot)` | `grperc·max(0,GPP−ΣRm)` | **formula in module** (`growth_respiration(npp_in, grf)`, ED2 form); budget deferred to allocation |
| Acclimation | none | CLM5 `rootstem_acc` | **seam reserved**, off in v1 |
| Flux unit | µmol / m² ground / s | gC / m² ground / s | **µmol / plant / s** |
| Grass | stem resp = 0 | `woody(ft)` gate | `is_woody` guard |

## 20. References

Stem: Chambers et al. 2004 (*Ecol. Appl.*); Rowland et al. 2018; Ryan 1990/1995, Pruyn 2002, Spicer &
Holbrook 2007 (sapwood-volume basis); Salomón et al. 2020 (TReSpire). Xylem-CO₂ confound: Trumbore et al.
2013, Angert et al. 2012, Bloemen et al. 2014. Fine root: Burton et al. 2002 (R∝root N, Q10 2.4–3.1);
Atkinson et al. 2007. Temperature / acclimation: Atkin & Tjoelker 2003; Atkin et al. 2015 (GlobResp); Reich
et al. 2016 (*Nature*); Ziegler et al. 2024 (*Science* — stem acclimation). Paradigm: Amthor 2000; Ryan
1991 (MR∝N). Cross-model: CLM5 Plant-Respiration Tech Note; FATES `FatesPlantRespPhotosynthMod.F90` +
`fates_params_default.json`; ELM `CNMRespMod.F90`. ED2: `dynamics/stem_resp_driv.f90`,
`dynamics/soil_respiration.f90`, `dynamics/growth_balive.f90`, `utils/allometry.f90`, `memory/pft_coms.f90`
+ `init/ed_params.f90`.

## 21. Open questions / decisions during implementation

1. **`allometry` home — DECIDED: its own `src/allometry/` library** (not folded into `meds_shared`).
   Rationale recorded in GitHub issue xiangtaoxu/MEDS#11 (allometry is a shared structural-geometry
   foundation that `state` depends on; keeping it inside `libmeds_plant` would make the demographic core
   depend on the whole ecophysiology library, breaking self-contained demography — see §3).
2. **Types — DECIDED: one `meds_plant_types` module** (§7a), overriding one-module-per-file for the types
   module only. No computational or logical consequence (compile-time only, no name collisions).
3. **Vital-rate file naming in `demography/`** — `meds_growth.f90` … (module names unchanged) vs a
   `meds_demography_rates_*` scheme.
4. **`src/plant/README.md`** — fold the four per-domain READMEs into one ecophysiology-module README.
5. **`wood_temp` proxy** — drive `stem_resp_norm` with canopy/air temperature until a wood energy balance
   exists; document at the call site.
6. **ρ-modulation of the stem factor** — enable in the first cut or reserve (reserved by default: the
   ρ↔parenchyma-N link is directionally sound but quantitatively scattered).
7. **Sequencing** — land Part I (refactor) green first, then Part II (respiration); never mix the reorg
   with new science.
