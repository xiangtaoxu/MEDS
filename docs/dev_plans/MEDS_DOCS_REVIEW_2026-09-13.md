# MEDS documentation review, 2026-09-13

**Status:** review record, read-only. No source or documentation was changed by this review; this
file is its only output. Reviewer: Claude (Fable 5.1). Baseline: `main @ b9596c5` (PR #147,
2026-09-11), 422 commits, 118 merged PRs since 2026-07-01, 14 open issues.

**Scope.** Every tracked `.md` file (70 files, ~41 k lines) read against the current `src/` tree
(87 modules, 30.3 k lines, 17 CMake libraries, 45 ctests), `CMakeLists.txt`, `git log`,
`gh pr list`, `gh issue list`, and the Claude Code documentation for memory files, rules and skills.
The 41 design plans were each checked three ways: their own status header, the modules/flags/keys
they said they would create (grepped in `src/` and `test/`), and the merged PR that shipped them.
Source comments in `src/` and `test/` were audited for change-narrative text.

**The five requests, one line each.**

| # | Request | Verdict |
|---|---|---|
| 1 | Dev plans: covered / out-of-date / invalidated | Archive 30 of 41, keep 11. Three "design-only" headers describe work that shipped the same day. Nine files share a banner whose item 3 was overturned by PR #90. One plan is gitignored yet cited. |
| 2 | README too verbose | 310 lines → about 110. Half the file is figure captions and build/config/output detail that lives elsewhere. Its carbon numbers are already stale against the regenerated example. |
| 3 | `src/` structure figure | Proposed `src/README.md` in §3: tree with line counts, four placement rules, the library DAG, a "where does a new file go" lookup. |
| 4 | CLAUDE.md evaluation | 777 lines, ~17 k tokens per session; the Claude Code docs target under 200 lines. Split into a ~150-line CLAUDE.md, five path-scoped `.claude/rules/*.md`, a gitignored `CLAUDE.local.md`, `docs/ROADMAP.md`, `CHANGELOG.md`. Skip AGENTS.md for now. |
| 5 | CHANGELOG | 434 history comment lines in 84 files, zero `TODO`/`FIXME` anywhere. Seed from the 118 PR titles, then sweep comments and eight reader-facing pages. |

---

## 0. Four things found on the way that are not documentation

1. **PR #144 did not flip `soil_carbon_on`.** The type default at `src/config/meds_config.f90:389`
   is `.true.`, but the loader at `src/config/meds_config_io.f90:814` reads the key with default
   `.false.`, and `toml_logical` returns that default when the key is absent (verified in
   `src/config/meds_toml.f90:150-170`). The root `meds_config_main.toml` has no `soil_carbon_on`
   key, so it runs with soil carbon OFF. Only the two `example_biophysics` TOMLs set it explicitly.
   One-line fix; also update the two comments that still say "opt-in, default .false."
   (`meds_config_io.f90:810`, `meds_demography_patch_fusefiss.f90:498`).
2. **`lwdown_source = "synthesize"` is a silent no-op.** The config parses it
   (`meds_config_io.f90:440-452`, `LW_SYNTHESIZE`) but `meds_met_driver.f90:320,493` always reads
   longwave from the file. Implement it or reject the value in `validate_config`.
3. **`MEDS_CAPI_REORG_DESIGN.md` is gitignored** (`.gitignore:69`, "kept on disk, never
   committed") yet the tracked `MEDS_CODE_STRUCTURE_DESIGN.md:534` cites it. Every other clone has a
   dangling reference. Decide: track it (then archive it) or drop the citation.
4. **The README's carbon numbers are stale.** README: GPP 437.2, R_eco 410.3, net 26.9 gC m⁻²,
   sink 372 of 744 h, CAS +25.9 ppm at night. `examples/example_biophysics/README.md`, regenerated
   2026-09-10/11 with soil carbon on and the Rh units fixed (#140): 434.6 / 288.9 / 145.7, 55 %,
   +17.2 ppm. The PNGs were regenerated (2026-09-11); the README prose was not (last touched
   2026-08-02).

---

## 1. Dev plans: archive, update, or keep

### 1.1 The convention has to change first

`docs/dev_plans/README.md` (whose title still reads "# archive") says superseded documents stay in
place and `archive/` is only for actively misleading ones. That rule is what left 41 files in one
flat directory, three of them with headers saying "design-only" for work that merged within a day.

Recommended rule: `dev_plans/` holds only plans with open items, or reference sections still cited
by section number from code. Everything else moves to `archive/` with a five-line tombstone:
status, shipping PR(s) with dates, what replaced it, where the live description now is (science
page or module), and "line references are the pre-2026-09 tree".

Moving is safe for source citations. Of the 166 comment lines in `src/` + `test/` that name a plan,
142 use the bare filename and 24 carry the `docs/dev_plans/` prefix; only the 24 need a path edit.
Three scripts (`scripts/numerics_sweep.py`, `parity_fidelity.py`, `parity_scenarios.py`) already
carry a doubled `docs/dev_plans/docs/dev_plans/archive/` path from the last move.

Verdict key: **ARCHIVE** = complete or superseded, tombstone + move; **ARCHIVE·INVALID** = the
conclusions are wrong, tombstone must say "do not cite"; **KEEP·LIVE** = un-built items still
intended; **KEEP·REF** = still cited by section number from code, needs a header rewrite.

### 1.2 Integrator and numerics family (12)

| Document | Lines | Verdict | Evidence, and what the tombstone must say |
|---|---|---|---|
| `MEDS_IMEX_ARK_DESIGN` | 362 | ARCHIVE | P0–P5 shipped 2026-07-11 by direct merge (e77f7f0, no PR); default since PR #88. Shipped tableau is ARS(2,2,2), not the planned ARK2(1)3L. P6 rejected by PRODUCTION §8. Cited from 5 src/test files by bare name. |
| `MEDS_INTEGRATOR_PHYSICS_PARITY_PLAN` | 639 | ARCHIVE | Phases 0–4 in PR #88; Phase 5 realised by the exact-exponential tissue store. Its own §7 still says Phases 4/5 "NOT implemented". Residue: `scripts/numerics_sweep.py` `--parity` preset pins `fast.leaf_energy_model`/`wood_energy_model`, keys that no longer exist. |
| `MEDS_INTEGRATOR_TEST` | 307 | ARCHIVE·INVALID | 24-run matrix scored against the retired split scheme, inside the frozen-conductance oscillation regime. The 900 s headline coincidentally matches today; the accuracy ranking and cost ratios do not. |
| `MEDS_NUMERICS_SCOPING` | 1225 | KEEP·REF | §5.1 process mask, §11 bare-array convention, §12.6 ED2 catalogue cited from 20 files. Scheme roadmap (Strang, ARK4, BDF2) superseded by PRODUCTION §8; §7 still says BB2/BB3 "COMMITTED" though GPU evaluation refuted them. |
| `MEDS_P3_COUPLED_SURFACE_DESIGN` | 421 | ARCHIVE | Not design-only: P3a–d shipped PR #37, replaced by the Newton arrowhead (2026-07-11), deleted with split (PR #88) and the `picard_*` keys (PR #105). Residues: `meds_fast_types.f90:50`, `meds_config_io.f90:725`. |
| `MEDS_PRODUCTION_INTEGRATOR_PLAN` | 2173 | KEEP·LIVE | The active numerics roadmap. Open: N5 adaptive freeze cadence; RK45 production warning (recommended, never added); E1 ψ-clamp artefact (#104); E5 RK45 rescue snapshot; §8 soil water into the tableau (#93 Phase 1); ψ_leaf still the one non-converged state at 900 s. §9 "where the code is" has 13 stale paths. Split the measurement record from the queue. |
| `MEDS_VEG_ENERGY_INTEGRATION_PLAN` | 551 | KEEP·REF | §1–8, 12–13 are the tissue-store design and still true. §9–11 and §14 ("no single coefficient is the culprit", "keep dt_fast ≤ 150 s") were overturned the same day by PR #90. Nine other files point here as "current state". Add a correction header. |
| `MEDS_ED2_RK45_DESIGN` | 1329 | ARCHIVE | Header says "not yet PR'd"; P0–P4 merged PR #67, P5/P6 PR #68. 53 source citations, all bare-name. Rescue target is now ARK; `INTEG_RK4` shipped as `INTEG_RK45`. |
| `MEDS_COLUMN_DYNAMICS_DESIGN` | 975 | ARCHIVE | Parts I/III/IV live as `meds_canopy_aerodynamics`, `meds_numerics`, `meds_budget_check` (PR #28). Part II operator-split sweep retired. 20 stale paths. |
| `MEDS_HIGH_PRIORITY_BUGFIX_PLAN` | 978 | ARCHIVE | All seven fixes in PR #29 (2026-07-07); bug 7's Dunne fix later mooted when Dunne was deleted. Cited by nothing. |
| `MEDS_GPU_EVALUATION` | 530 | KEEP·LIVE | Measurement stands (PR #110). Five of seven recommendations open: mark BB2/BB3 refuted in SCOPING §7; fix the GPU overselling in `CMakeLists.txt:6` and `CLAUDE.md:76`; allocator self-time in `build_column_frozen`; cohort-axis threading; `wp = real32` experiment. Not indexed anywhere. |
| `archive/MEDS_INTEGRATOR_PARITY` | 1309 | already archived | Tombstone adequate. Fix: its "where to look" points at VEG_ENERGY §9–14 (itself overturned); `examples/example_biophysics/meds_config_july.toml:29` quotes its invalidated "~21× more accurate" number without a [RETIRED] tag. |

### 1.3 Plant and output family (9)

| Document | Lines | Verdict | Evidence, and what the tombstone must say |
|---|---|---|---|
| `MEDS_PLANT_ECOPHYSIOLOGY_DESIGN` | 592 | ARCHIVE | Part I layout (flat `src/plant/`, façade, `src/allometry/`) shipped 2026-07-04 and was replaced by the structure plan. Part II respiration shipped PR #13. Caveat: §11–12 is the only written record of the stem/root maintenance-respiration equations (see §1.6). 18 stale paths. |
| `MEDS_PLANT_CARBON_ALLOCATION_REFACTOR_DESIGN` | 278 | ARCHIVE | PR #52; `docs/science/plant_carbon_allocation.md` covers it. Cited by path from the kernel header and two science pages; repoint. |
| `MEDS_PLANT_CARBON_DYNAMICS_DESIGN` | 332 | ARCHIVE | Superseded banner already present. Five comments still name the deleted `meds_plant_carbon_dynamics` module. |
| `MEDS_PLANT_TRAIT_DYNAMICS_DESIGN` | 194 | ARCHIVE | PR #54; `plant_traits.md` covers it. Carry §5 thermal acclimation to the roadmap. |
| `MEDS_PHENOLOGY_DESIGN` | 339 | ARCHIVE | Tri-state design retired by the rate refactor; no banner today, so its §4 equations read as current. §10 cites the removed `-DMEDS_ENABLE_IO=OFF` build. |
| `MEDS_PHENOLOGY_RATE_REFACTOR_DESIGN` | 665 | KEEP·LIVE | P3 verified un-wired: the driver hard-codes `avail_water = 0`, `dmax_leaf_psi = 0`, `rad = 0` (`meds_vegetation_dynamics.f90:936-943`) and `validate_config` still rejects WATER/HYDRO/LIGHT cues (`meds_config.f90:662`), though the kernel side is coded and the fast loop now produces all three drivers. P4 `retained_carbon_fraction` and P5 declination unification also open. |
| `MEDS_IO_DESIGN` | 1164 | KEEP·REF | Cited by section from CMake, 3 config modules, 2 io modules, 3 tests. Banner needs: FAST tier shipped 2026-07-13 (its §9 deferral is obsolete); §3.5 table superseded by IO_V01 §4; async writer shelved; variance still deferred. |
| `MEDS_IO_V01_PLAN` | 1099 | KEEP·REF | PR #111; 208 registered variables (doc says ~203). Cited from 12 files + CLAUDE.md. Open: mortality by pathway; disturbance area (`PD_DISTURB_AREA` declared, never written); per-band albedo; FAST staging unification (`fast_sample_t`, `extract_fast_scalar`, `output_integrate_fast` were slated for deletion and still exist); `[io]` → `[state]` rename. |
| `MEDS_FAST_OUTPUT_DESIGN` | 299 | ARCHIVE | Shipped 2026-07-13 by direct merge (01347bb, no PR). Its deferred P2 was absorbed by IO_V01 D1. Six hard line-number citations are stale. |

### 1.4 Structure, reorg and review family (9)

| Document | Lines | Verdict | Evidence, and what the tombstone must say |
|---|---|---|---|
| `MEDS_CODE_STRUCTURE_DESIGN` | 1894 | KEEP·LIVE, split | Reference (13 %): §1 decisions, §4 target tree, §5 straddlers, §6 placement rules, §7.1–7.2. Migration log (76 %): §2–3, §8–14, §15.0/15.3/15.5; §10.2 alone is 549 lines of conservation measurements and could be its own record. Roadmap (11 %): §15.1, 15.2, 15.4, 15.6–15.8. §4 tree has drifted (nine module renames, seven modules missing). Decision #13 (test/ mirrors tree) never decided. Stale statements: §12.6 "`meds_biophysics_interface` stays" (gone), §13.5 `reconcile_tissue_water_capacity` "open" (exists), §13.8 `agf_bs` "open" (closed #130), §14.2 "three shims" (four). |
| `MEDS_CODE_REVIEW_2026-07-06` | 483 | ARCHIVE | Self-marked fully addressed (PRs #29–#34). Line references are the pre-September tree. |
| `MEDS_CODE_REVIEW_2026-09-08` | 562 | ARCHIVE | Executed by PRs #119–#124; leftovers re-homed into STRUCTURE §15.6/15.7. Header still says "items 4–10 not covered". |
| `MEDS_REORG_DESIGN` (v6) | 290 | ARCHIVE·INVALID header | Header: "DESIGN-ONLY, not implemented, no branch". Reality: PR #43 merged P0–P4 the same day (2026-07-16). Body is complete; the header gets the history backwards. |
| `MEDS_DRIVER_REORG_DESIGN` | 403 | ARCHIVE | Header: "not executed". Reality: PR #62 shipped Phases 1–5 (2026-07-21). Cited by `meds_fast_rk4_oracle.f90:5`. One undecided item: §7 retire the IMEX-Euler oracle tier. |
| `MEDS_CORE_MODULE_REORG_DESIGN` | 309 | ARCHIVE | Shipped PR #44/#45, then every `meds_core_*` name was reversed on 2026-09-09. 83 stale module names, the most of any file. |
| `MEDS_BIOPHYSICS_DEDUP_DESIGN` | 267 | ARCHIVE | PR #58; both deferred items closed later (#59, #88). Cited by nothing. |
| `MEDS_BIOPHYSICS_REORG_DESIGN` | 234 | ARCHIVE | PR #57; every module it created was moved or renamed again by the structure plan. |
| `MEDS_CAPI_REORG_DESIGN` | 315 | ARCHIVE, untracked | Gitignored. One-`.so` principle adopted (PR #138), 9-file split rejected (STRUCTURE §14.2). Residue: `meds_capi_demography.f90:165-169` still inlines allometry; `examples/example_demography/empirical_laws.py:18` hard-codes B1HT/B2HT. |

### 1.5 Biophysics and column family (12)

| Document | Lines | Verdict | Evidence, and what the tombstone must say |
|---|---|---|---|
| `MEDS_BIOGEOCHEMISTRY_DESIGN` | 1127 | KEEP·LIVE | P0 (PR #35) and P3 (PR #64: per-patch state, `[soil_carbon]`, restart, litter seam) landed. P1 half: SASU and EXPM done; the DAMM kernel exists with no TOML key and no caller. N cycle is parsed but no kernel reads it. P2 (per-layer pools, N limitation, fire) unbuilt. §9 "Ra is still 0" is stale. The only written record of the CENTURY matrix; no science page. |
| `MEDS_COLUMN_CO2_BALANCE_DESIGN` | 1580 | ARCHIVE | Twin shipped PR #26 under different names (`cas_column_*` in `meds_cas_biophysics`); Rh authority moved to the CENTURY matrix. §2 placement invalidated by the reorg. Q10 and DAMM Rh kernels both have zero callers. 74 stale references, the most of any file. Carry §3.5 multi-layer CAS to the roadmap. |
| `MEDS_COLUMN_HYDROLOGY_DESIGN` | 748 | ARCHIVE·INVALID | Its own "P2 IMPLEMENTED: retention-integral Zeng–Decker" (line 714) is now false: ZD, the lumped aquifer store, baseflow and `z_wt` were deleted 2026-07-30 (55ffca0); the aquifer is a head-driven boundary. Everything else (Richards solve, retention curves, BCs, ponding) landed. |
| `MEDS_ENERGY_BALANCE_DESIGN` | 752 | ARCHIVE | P0–P2a shipped PRs #24/#25. Its P3 "coupled fixed point" is realised as the ESDIRK2 Newton arrowhead, not the Picard it planned. Five science pages link to it by path; retarget. Carry §14 P2b skin layer and P5 multi-layer CAS. |
| `MEDS_LEAF_WOOD_ENERGY_DESIGN` | 353 | ARCHIVE | P0–P3 PR #41; P4 bordered arrowhead abandoned for the exact-exponential store. Already banner-marked superseded by VEG_ENERGY; extend to a tombstone. |
| `MEDS_SNOW_DESIGN` | 485 | KEEP·LIVE | P0 PR #42. Header and §6 still say "ARK stays snow-free (split-path only)": resolved, ARK/RK45 run the shared snow stage since PR #77/#80. Live: P1 multi-layer snow, compaction, aging albedo; P2 canopy snow interception. Kernels live in `meds_ground_biophysics`, not the named modules. |
| `MEDS_HYDRAULICS_DESIGN` | 580 | ARCHIVE | PR #9 then PR #49. P5 leaf↔hydraulics Picard superseded by ARK plus the PR #91 corrector. `plant_hydraulics.md` covers it. |
| `MEDS_HYDRO_CURVE_EXTRACTION_DESIGN` | 308 | ARCHIVE | Header "design-only"; `src/functions/meds_hydr_lib.f90` is exactly this extraction, shipped PR #49 the same day. Cited by nothing. Per-PFT hydraulics remains optional roadmap. |
| `MEDS_MULTILAYER_ROOTS_DESIGN` | 161 | ARCHIVE | Phase A shipped PR #49; the opt-in flag deleted 2026-07-30 so the body's "bit-identical when off" goal is gone. Carry Phase B per-layer root nodes and hydraulic redistribution. |
| `radiative_transfer_design` | 433 | ARCHIVE | PR #5. Line 14 "the code does not exist yet" and the "biophys:" title prefix are stale; six planned modules collapsed to two. `dev_plans/README.md:24` still says `src/biophys/`. |
| `MEDS_FORCING_DESIGN` | 1247 | KEEP·LIVE | P0–P3 shipped (PRs #36, #69). Also the only specification of the forcing NetCDF format (§7.1) and the ERA5 de-accumulation recipe (§7.3), cited from two scripts and CMake. Live: LWdown synthesis (see §0), multi-polygon runtime, transient CO₂ stream. `src/forcing/README.md:100` lists the daily accumulator as deferred; it exists as `site%pheno_tair_sum/pheno_tair_n`. No science page. |
| `MEDS_SLOW_DYNAMICS_DESIGN` | 408 | ARCHIVE | Header "design-only" is stale: Part I PR #63, Part II PR #64. Cited from 18 source comments by bare name (safe to move); 16 dead `src/driver/…#L` anchors inside the doc. |

### 1.6 Three science pages are missing, and they gate three archive moves

Archiving is clean only when the live equations have another home. Three do not:

- **Soil carbon**: the CENTURY matrix ODE, SASU, litter partition and the `xi_int` seam exist only
  in `MEDS_BIOGEOCHEMISTRY_DESIGN`. `docs/README.md` already lists `science/soil_carbon.md` as
  planned.
- **Meteorological forcing**: cosz reconstruction, shortwave partition, precipitation phase,
  recycling and the file format exist only in `MEDS_FORCING_DESIGN`.
- **Stem and root maintenance respiration**: only in `MEDS_PLANT_ECOPHYSIOLOGY_DESIGN` §11–12.

Write `docs/science/soil_carbon.md`, `forcing.md` and `plant_respiration.md` first, or keep those
three plans as references until then.

### 1.7 The shared banner

Nine files carry the identical "SUPERSEDED IN PART — 2026-07-31" banner (IMEX_ARK, PHYSICS_PARITY,
INTEGRATOR_TEST, NUMERICS_SCOPING, P3_COUPLED_SURFACE, ED2_RK45, COLUMN_DYNAMICS, DRIVER_REORG,
LEAF_WOOD_ENERGY). Items 1 and 2 (split retired; "IMEX-ARK" is ESDIRK2) are right. Item 3 says
`dt_fast` is stability-limited with a 150 s default and sends readers to VEG_ENERGY §9–14. PR #90
merged the same day, removed the oscillation and set 900 s; `meds_config.f90:535` now says "an
ACCURACY parameter". Rewrite the banner once and propagate to all nine.

### 1.8 Proposed `docs/dev_plans/` after the move

**Stays (11 + index).**

- `README.md`, retitled (it currently says "# archive"), indexing every file with one line and a
  verdict.
- Live (7): `MEDS_CODE_STRUCTURE_DESIGN` (trimmed to §1/4/5/6/7 + §15),
  `MEDS_PRODUCTION_INTEGRATOR_PLAN`, `MEDS_PHENOLOGY_RATE_REFACTOR_DESIGN`, `MEDS_GPU_EVALUATION`,
  `MEDS_BIOGEOCHEMISTRY_DESIGN`, `MEDS_SNOW_DESIGN`, `MEDS_FORCING_DESIGN`.
- Reference (4): `MEDS_NUMERICS_SCOPING`, `MEDS_IO_DESIGN`, `MEDS_IO_V01_PLAN`,
  `MEDS_VEG_ENERGY_INTEGRATION_PLAN`. Each gets a rewritten status header dated 2026-09.

**Moves to `archive/` (30, of which 3 are "invalid").** Each gets the five-line tombstone. Three
headers are corrected in the same commit (REORG v6, DRIVER_REORG, CORE_MODULE: "design-only" →
shipped). The open items below are extracted into `docs/ROADMAP.md` *before* the move, so nothing
intended is lost.

### 1.9 Items to extract into `docs/ROADMAP.md`

Collected from the plans above, `CLAUDE.md` "Reserved follow-ups", and the ~40 "deferred / MVP /
placeholder / not yet" comment markers in `src/` and `test/`.

- **Phenology** (RATE_REFACTOR §9): P3 wire WATER/HYDRO/LIGHT cues from the fast loop (cohort
  columns `pheno_water_avg`, `pheno_low_psi_days`, `pheno_high_psi_days`, `pheno_light_avg`; daily
  max leaf ψ; lift the `validate_config` rejection); P4 `retained_carbon_fraction`; P5 one
  `solar_declination(doy)` for `solar_cosz` and `daylength`.
- **Plant physiology**: thermal acclimation (TRAIT §5, ECOPHYS §17); storage respiration; layered
  root respiration; per-PFT hydraulics (HYDRO_CURVE §8); Phase B per-layer root nodes and hydraulic
  redistribution (ROOTS §5, HYDRAULICS §16).
- **Output** (IO_V01 §4.5–4.6, §6): mortality carbon by pathway; disturbance area flux
  (`PD_DISTURB_AREA` declared, unwritten); per-band albedo and up-welling SW/LW; unify FAST
  extraction through `extract_variable` and delete `fast_sample_t`/`extract_fast_scalar`/
  `output_integrate_fast`; rename `[io]` → `[state]` and purge "legacy [io]" comments; variance
  companion (`AGG_MEANSQ`, IO_DESIGN §3.5).
- **Production numerics** (PRODUCTION §5–8): N5 adaptive freeze cadence with a real error
  estimator; RK45 production-cadence warning; E1 `rwc_floor` ψ-clamp artefact (#104); E5 RK45
  rescue snapshot; fold soil water into the ARK tableau (#93 Phase 1); ψ_leaf non-convergence at
  900 s.
- **Vegetation energy** (VEG_ENERGY §6–7): separate canopy film store with phase change; retire
  `veg_energy_step_implicit` (unify on `veg_energy_diagnostic`); 1.25·h free-convection slope;
  honest wood sizing (`bsap` placeholder, `meds_vegetation_dynamics.f90:617`).
- **Numerics scoping** (SCOPING §7, §11.3): mark BB2/BB3 refuted; MB2 soil-energy substepping
  (adaptive knobs in `energy_opts_t` never read); bare-array forms of `cas_column_step_implicit`,
  `soil_energy_step_implicit`, `soil_carbon_step`, snow kernels.
- **GPU evaluation** (§12): fix the GPU overselling in `CMakeLists.txt:6` and `CLAUDE.md:76`;
  allocator self-time in `build_column_frozen` (~24 %); cohort-axis threading; `wp = real32`
  experiment.
- **Structure** (CODE_STRUCTURE §15): pass `column_params_t` via `column_config_t` instead of
  copying into `column_frozen_t`; per-layer face `budget_imbalance` on the committed path; delete
  `column_cohort_t` (38 refs in 11 files) for `*_fast_slice_t`; consolidate the 14 per-test `check`
  routines into `meds_test_support`; frozen-seam contract note (Phase 3, optional); packed
  `column_state_t` (#146); decide the year-rollover `seam[soil_carbon_rh]` residual (§15.2);
  decision #13 (test/ mirrors tree) undecided.
- **Biogeochemistry** (BIOGEOCHEM §7): DAMM as a selectable moisture response (add `hr_model`
  key or delete the dead kernel); nitrogen twin; vertically resolved pools; N limitation of NPP;
  fire consumption of ground pools; explicit CWD pool; multi-layer CAS (CO2 §3.5, ENERGY §14 P5).
- **Snow** (SNOW §7): P1 multi-layer snow, compaction, aging albedo, density-dependent `k_snow`;
  P2 canopy snow interception.
- **Forcing** (FORCING §5.7, §8): LWdown synthesis (or reject the value); multi-polygon runtime
  (`nearest_grid_index` is the atom); transient/observed CO₂ stream.
- **Small residues**: Phase 6 remove the dead `--parity` preset (PHYSICS_PARITY §7); retire the
  IMEX-Euler oracle tier or record the decision (DRIVER_REORG §7); C-API shim inlines allometry
  (CAPI residue); soil bottom thermal BC (#145).
- **Open GitHub issues**: #1, #6, #7, #47, #74, #89, #96, #104, #114, #117, #118, #145, #146,
  #148.

Several current markers are themselves stale and should be deleted rather than moved: "optics
PFT-UNIFORM" (`meds_fast_dynamics.f90:188`, fixed by #131); "SPLIT PATH ONLY … error-stops under
INTEG_ARK" (`meds_fast_types.f90:164`, no such stop exists); "LAI-share SW split MVP until RT join"
(`meds_fast_dynamics.f90:930`, the join is done; this is the no-forcing fallback); the `gpp_ref`
"fast-off fallback" premise (`meds_config.f90:341`, the fast loop is always on); "rhizosphere
conductance single-layer MVP" (`meds_config.f90:160`, per-layer path is unconditional).

---

## 2. README: what an ecologist should get in two minutes

### 2.1 Diagnosis

The README is 310 lines and 22 KB. Its shape today:

| Block | Lines | Problem |
|---|---|---|
| Opening paragraph + demography GIFs | ~20 | Good. Keep. |
| Fast-loop figure + caption | ~35 | Caption carries ~15 pinned numbers and a paragraph of history ("removed the old period-2 oscillation"). |
| Carbon + soil figures + captions | ~35 | ~25 more pinned numbers, several now wrong (§0 item 4). Both figures are already in `examples/example_biophysics/README.md` with current numbers. |
| Design goals | 10 | Generic; says nothing an ED2 user cannot guess. |
| Status + Highlights | ~60 | A changelog in disguise ("always on", "no longer", "three interchangeable integrators", "some 200 variables"). |
| Building + Installing dependencies | ~50 | Duplicates CLAUDE.md and the example READMEs. |
| Configuration | ~40 | The only place the two-file layout, init modes and parameter philosophy are written. Must survive, but not here. |
| Dependencies & environment | ~20 | Overlaps Building. |
| Output & post-processing | ~55 | Four script descriptions that belong beside the scripts; the output system is documented in `docs/science/diagnostics.md`. |
| Scientific reference + pointer to CLAUDE.md | ~8 | Keep; point to `src/README.md` instead of CLAUDE.md. |

The novelty (a demographic model whose sub-daily surface physics is solved rather than
prescribed, in modern testable Fortran, drivable from Python) is there, but a reader has to find it
between captions.

### 2.2 Proposed outline (~110 lines)

1. Title, one paragraph: ED2 lineage, Fortran 2018, site model, two coupled timescales, what is
   different (solved surface physics, conservation ledgers, Python-drivable).
2. Two figures: the demography pair as now; one biophysics figure (`biophysics_july.png`) with a
   three-sentence caption and **no pinned numbers**. Link to the example for the numbers.
3. "What the model does": a 10–12 row table, process · scheme · where documented. This replaces
   Status, Highlights and Design goals. Rows: canopy radiation (ED2 two-stream), aerodynamics
   (CLM5 MO + ED2 Nusselt), canopy air space (3 prognostic twins), leaf gas exchange (FvCB/Collatz,
   3 stomatal models), plant hydraulics (matrix-exp, multi-layer roots), soil water (implicit
   Richards, vG/Campbell), soil and tissue energy (internal energy, freeze/thaw by inversion), snow,
   phenology (rate-based), carbon allocation (PARTEH-H1), soil carbon (CENTURY matrix),
   demography (cohort/patch fuse-fission, treefall), integrators (ESDIRK2 default, RK45),
   output (~200 variables × 4 timescales), parallelism (patch threading, byte-identical).
4. Quick start: six lines (install netCDF via the script, configure, build, ctest, run the
   demography example, plot).
5. Where next: `examples/`, `docs/science/`, `docs/ed2_comparison.md` ("coming from ED2?"),
   `src/README.md` (code structure), `python/README.md`, `docs/ROADMAP.md`, `CHANGELOG.md`.
6. Scientific references.

### 2.3 What moves where

| From README | To |
|---|---|
| Carbon and soil figure blocks | Drop; already in `examples/example_biophysics/README.md` |
| Building, compiler activation, installer scripts | New `docs/building.md` (short), or a `scripts/README.md` |
| Configuration (two files, init modes, parameter philosophy, derived block) | New `docs/configuration.md` |
| Output & post-processing (four scripts) | New `post_proc/README.md`; output system → `docs/science/diagnostics.md` |
| Highlights that describe history | `CHANGELOG.md` |

Also fix, in the same PR: `examples/README.md:11` and `examples/example_leaf_gas_exchange/README.md:4`
still say the leaf kernel is "not yet coupled to the demographic spin-up" and cite `src/plant/`;
`examples/example_phenology/README.md:3` cites `src/plant/meds_phenology.f90`.

---

## 3. The `src/` figure for a new `src/README.md`

There is no `src/README.md` today; the tree lives as an eight-line sketch inside CLAUDE.md and a
drifted §4 inside the structure plan. The proposal below is drawn from the current tree and
`CMakeLists.txt` at b9596c5, so line counts and library names are exact. It is meant to be the one
place a new file's home is looked up, and the page to cite when introducing the code.

```
src/                                    30.3 k lines · 87 modules · 17 CMake libraries
│
│  ── FOUNDATION ─────────────────────── libmeds_shared: no model state, no process
├── base/          kinds, physical and calendar constants                      109
├── functions/     stateless constitutive laws: allometry, thermodynamics,
│                  retention/PV curves, canopy optics, temperature response   1 085
├── util/          calendar time, numerics (expm, root finders), budget checks  882
├── config/        TOML reader and loader, PFT trait table, meds_config_t,
│                  the per-domain *_opts leaves                               2 907
│
│  ── STATE, a layer in two halves ───────────────────────────────────────────
├── state/
│   ├── column/    ONE patch seen vertically: soil water, soil energy, snow,
│   │              canopy air, soil carbon, and their parameter bundles         558
│   └── site/      ALL patches: cohort SoA, patch CSR, lockstep reorder,
│                  site_t, diagnostic accumulators                            1 674
│
│  ── PROCESSES, timescale first ──────────────────────────────────────────────
├── fast_dynamics/       sub-daily, dt_fast                                 ≈10 900
│   ├── canopy/    the medium: two-stream radiation, aerodynamics, CAS box      983
│   ├── plant/     the organisms: leaf gas exchange, hydraulics,
│   │              maintenance respiration, tissue energy                      1 606
│   ├── soil/      the ground: soil water, soil energy, ground skin, snow      1 472
│   ├── numerics/  integrators (ARK default, RK45), state vector, frozen
│   │              record, pre-pass, error control, RK4 test oracle           5 378
│   └── driver/    walks one slow step in dt_fast sub-steps over patches      1 538
│
├── slow_dynamics/       daily to annual                                     ≈4 600
│   ├── plant/     phenology, carbon allocation, trait plasticity                587
│   ├── soil/      CENTURY soil-carbon matrix, litter partition                  912
│   ├── demography/ vital-rate LAWS and the operators that apply them:
│   │              update, cohort/patch fuse-fission, recruitment, treefall   1 308
│   └── driver/    slow coordinator, vegetation and biogeochem drivers,
│                  the slow conservation ledger                               1 816
│
│  ── EDGES ────────────────────────────────────────────────────────────────────
├── forcing/       prescribed drivers: met reader, disaggregation kernels      1 067
├── io/            netCDF C bindings, restart stream, the diagnostic wall
│                  (derive → capture → reduce → integrate → serialize)        4 563
├── init/          initial community: bare ground, census                       181
├── capi/          bind(c) shims → one libmeds.so (leaf, phenology,
│                  demography, full run)                                        887
└── main/          meds_stepper (cadence), meds_driver (open/step/finalize),
                   meds_main (the PROGRAM)                                       740
```

**Four rules, in the order you apply them.**

1. **Kernel folders never see `site_t`.** `fast_dynamics/{canopy,plant,soil}`,
   `slow_dynamics/{plant,soil}` and `state/column` link `state_column` + `config` only, so each
   kernel library builds standalone and stays OpenMP-target eligible. If a new routine needs
   `site_t`, it is driver code.
2. **A kernel goes where its caller's timescale is**, then in its domain folder. A kernel called
   from both tiers is a documented seam, never a folder. There is one:
   `heterotrophic_respiration_matrix`, so sub-daily Rh respires the same CENTURY matrix the daily
   step debits.
3. **A derived type lives with whoever mutates it.** Two mutators means boundary state, which is
   `state/column`. Parameters are derived once and never integrated, so `meds_column_params` is a
   separate module from `meds_column_state_types`.
4. **Laws and operators do not touch.** `demography/` holds the vital-rate laws
   (`meds_demography_rates`) and the appliers (`update_*`, `*_fusefiss`, `apply_recruitment`,
   `apply_patch_disturbance`); the appliers take rate arrays and never import the laws. The slow
   driver is the one place they meet, and the Python `apply_rates` path is the standing test.

**Library DAG** (acyclic; from `CMakeLists.txt`):

```
shared ─→ state_column ─→ config ─→ state_site ─→ demography ─→ io_prep ─┐
             ├─→ fast_kernels ──────────────────────────────────┐        │
             └─→ slow_kernels ──────────────────────────┐       │        │
netcdf_c ─→ forcing ────────────────────────────────────┴───────┴────→ fast ─┐
                                                       demography ────→ slow ─┼→ stepper ─→ model
                                                       demography ────→ init ─┘   (INTERFACE)
io_stream (meds_io + output_stream/manager) ──────────────┐
model ────────────────────────────────────────────────────┴→ driver ─→ meds_main │ meds_py (libmeds.so)
```

**Where does a new file go?** A short lookup to end the README:

| New thing | Home |
|---|---|
| a soil or plant process | `fast_dynamics/<domain>` or `slow_dynamics/<domain>` by cadence |
| a reservoir two subsystems mutate | `state/column` |
| a per-cohort field | `state/site` + the lockstep reorder + fusion policy + creation sites |
| a constitutive curve | `functions/` |
| a TOML block | `config/` as an `*_opts` leaf; `meds_config` carries it |
| an output variable | one `add_variable` line in `io/meds_output_registry` |
| a Python entry point | `capi/` + its mandatory ctest target |

Worth stating in that README because code comments currently disagree with it: `meds_stepper.f90:3`
says it "lives in `src/driver/`", `meds_vegetation_dynamics.f90:3` says "compiled into meds_aux",
`meds_config_io.f90:3` says "lives in `src/io/`", `meds_capi_demography.f90:3` says `libmeds_c`.
Seventeen module headers cite paths or libraries that no longer exist.

---

## 4. CLAUDE.md: too much, in the wrong place

### 4.1 Diagnosis

CLAUDE.md is 777 lines and 67 KB, roughly 17 000 tokens loaded into every session. The Claude Code
documentation targets under 200 lines per file and says adherence drops as the file grows. Section
sizes: "Source layout & libraries" 380 lines (49 % of the file), "Invariants" 87, "Repository
status" 62, "Build" 62, "Architecture (inherited from ED2)" 57, "Modernization guidelines" 38,
"Writing docs — GitHub math" 25, "Reserved follow-ups" 26. The 380-line section is subsystem prose
that duplicates the module headers, the science pages and the four folder READMEs (zero verbatim
overlap, so it is a fourth independent copy to keep in sync). It also carries this machine (conda
prefix, oneAPI path, HPC SDK path, the RTX card) in a tracked file.

Accuracy has drifted with the restructure. Verified stale claims:

| Line | Claim | Reality |
|---|---|---|
| 34 | `app/` | does not exist |
| 76 | nvfortran "is the parallel/GPU path" | GPU evaluation: not viable as scoped |
| 112 | "38 tests" | 45 |
| 141 | libraries `meds_biophysics`/`meds_plant` | do not exist |
| 244 | seam `meds_leaf_physiology%leaf_gas_exchange` | module is `meds_leaf_gas_exchange` |
| 254 | `meds_plant_phenology` + `meds_pheno_engine` | module is `meds_phenology` |
| 255 | `meds_plant_carbon_dynamics` | `meds_plant_carbon_allocation` |
| 264 | plant kernels "NOT yet wired into the demographic stepper" | contradicts the run-model paragraph at line 44 |
| 294, 312 | `meds_soil_solver` | no such module (Thomas sweep is in `meds_numerics`) |
| 318 | "coupled fixed point deferred to P3" | history; P3 shipped, then replaced by ESDIRK2 |
| 369 | "`src/utils/` remains an empty placeholder" | does not exist |
| 385, 439 | `libmeds_io` | target is `meds_io_stream` |
| 412 | "legacy span-wrap kept for non-calendar files" | PR #69 made it a hard error |
| 631–640 | C-API for the demographic engine listed as future | exists (`meds_capi_demography`, `meds_capi_run`) |

### 4.2 Proposed layout

Two facts from the Claude Code docs shape this. `@import` lines still load at launch, so splitting
CLAUDE.md into imports saves nothing; **path-scoped `.claude/rules/*.md`** (YAML frontmatter
`paths: [...]`) is the mechanism that keeps text out of context until a matching file is touched.
And Claude Code does not read `AGENTS.md`; if another tool (Codex, Cursor) is ever used on MEDS,
make CLAUDE.md a one-line `@AGENTS.md` import. Until then a separate AGENTS.md is a second copy to
keep in sync.

| File | Loads | Contents (from today's CLAUDE.md) |
|---|---|---|
| `CLAUDE.md` (~150 lines) | always | What MEDS is (5 lines). Build and test commands (three blocks, `$CONDA_PREFIX` not a hard path). The four placement rules and "kernels never see `site_t`". The six engine invariants as one line each with a pointer: lockstep reorder, persistent ids, conservation asserts, tendencies as data, fusion policy declared once, never read `site%deriv`. Naming convention. "Test as you port; build nvfortran too." Portability traps (array-valued function results, `-auto`, BLOCK in parallel regions, static mold). Pointers to `src/README.md`, `docs/ROADMAP.md`, `CHANGELOG.md`, `docs/dev_plans/README.md`. Git conventions. A four-line PR checklist (§5.5). |
| `.claude/rules/fast-loop.md` | `src/fast_dynamics/**`, `test/test_column_*.f90`, `test/test_fast_loop.f90` | Frozen record and dt_fast semantics, the two integrators and the rescue, the `column_view` filler, fast-loop fusion policy, ledger expectations, the per-thread scratch pool rationale. |
| `.claude/rules/state-demography.md` | `src/state/**`, `src/slow_dynamics/**`, `src/init/**` | SoA + CSR, lockstep machinery and creation sites, conserved invariants, order of operations, laws vs operators, the slow ledger. |
| `.claude/rules/output.md` | `src/io/**`, `src/config/meds_output_config.f90` | Five-stage wall, weight/mean/scale contract, axes, range-partitioned source ids, explicit CMake list, `[io]` is restart-only. |
| `.claude/rules/config.md` | `src/config/**`, `**/*.toml` | No hard-coded parameters, two-file layout, presence map, `[soil_column]` vs `[soil]`, defaults live only in `build_test_config`. |
| `.claude/rules/docs-math.md` | `docs/**/*.md` | The GitHub math section, verbatim. |
| `CLAUDE.local.md` (gitignored) | always, this machine | oneAPI and HPC SDK activation paths, conda prefix, the GPU present, the untracked `runs/` test beds. |
| `docs/ROADMAP.md` | on demand | Replaces "Reserved follow-ups" and every deferred marker (§1.9). |
| `docs/ed2_comparison.md` §6 | on demand | Absorbs the 57-line "Architecture (inherited from ED2)" file map. It is reference for a porter, not a rule for every session. |

A `TODO.md` inside `.claude/` is not recommended: the roadmap is a project document and belongs in
`docs/`, referenced from CLAUDE.md. Optional and lazy-loaded (cost is one description line until
invoked): `.claude/skills/build-test/SKILL.md` holding the three configure/build/ctest recipes and
the "a green ifx run is not sufficient" rule; `.claude/skills/ledger-check/SKILL.md` for the
conservation-run harness.

---

## 5. CHANGELOG: one home for "what changed, when"

### 5.1 Diagnosis

There is no CHANGELOG, ROADMAP or TODO anywhere in the repo, and no `TODO`/`FIXME` token in any
source file. History and intent are instead spread over 434 comment lines in 84 Fortran files
(291 in `src/`, 143 in `test/`), 40 passages in CLAUDE.md, and change stories in eight
reader-facing pages. The pattern is consistent: a comment says what the code *used to* do and which
PR changed it, where a reader only needs why the code is the way it is now.

History-type lines concentrate in the fast loop: `test_column_dynamics` 30, `meds_fast_frozen` 26,
`meds_fast_rk45` 25, `test_column_rk45` 23, `meds_fast_types` 21, then `meds_site_state_types`,
`meds_soil_water`, `meds_fast_ark`, `meds_fast_dynamics`, `meds_config` at 14 each.

### 5.2 Shape of the file

- **Location.** `CHANGELOG.md` at the repository root is where GitHub and release tooling look;
  `docs/CHANGELOG.md` works if the README links it. Keep a Changelog format: `Unreleased`, then
  `0.1.0 — 2026-08-02`, then a `pre-0.1` section grouped by subsystem for 2026-06-23 → 2026-08-02.
- **Entry shape.** One line, PR link, and for anything that changes numbers the before → after
  with the magnitude ("CENTURY Rh reached the atmosphere 964× too small; fixed, #140").
- **Seed.** `gh pr list --state merged --limit 200` returns all 118 titles with merge dates. Two
  features shipped without a PR number and need a hash: IMEX-ARK (e77f7f0, 2026-07-11) and the
  FAST output tier (01347bb, 2026-07-13). The per-plan audits in §1 list every PR, issue, hash and
  date cited inside the plan headers.

### 5.3 What moves out of comments

Eight worked examples of the edit, each with the sentence that survives:

| Where | Today (history) | Stays (rationale) |
|---|---|---|
| `meds_config.f90:61-64` | "the operator-split third was RETIRED 2026-07-31 (it converged to a different limit…). The `integration_scheme` selector went with it." | "`[fast].time_integrator` selects ARK (default) or RK45; see numerical_scheme.md §3 for why there is no operator-split scheme." |
| `meds_fast_step.f90:7-8` | "The operator-split integrator that used to live here was retired 2026-07-31…" | nothing; lines 4–6 already state the module's role. |
| `meds_column_state_types.f90:34-39` | "The pond used to be a MASS buffer with no thermal state… The books closed but the column lost energy (#78 item 4)." | "Prognostic pond enthalpy so every pond seam is a paired (mass, enthalpy) transfer; a mass-only pond lets energy leave the ledger while the water still holds it." |
| `meds_site_state_types.f90:305-311` | "This used to live only on the per-patch SCRATCH, which BB1 phase 1 hoisted OUT of the patch loop, so it was loop-carried…" | "Per-patch controller state stored as a column reservoir, so the answer is independent of patch order and thread scheduling." |
| `meds_fast_types.f90:169-172` | "The picard_* mirrors that used to sit here were DELETED (plan E4)… The comment claiming … was false." | delete, or one line: "ARK's Newton cap is `NEWT_MAX` in meds_fast_ark, not a config field." |
| `meds_io.f90:5-9` | "This module used to carry a second, DIAGNOSTIC writer … retired at v0.1 … COLLIDED on the -D- prefix." | "Restart only. Diagnostics are the `[output]` registry subsystem; a checkpoint is raw prognostic state at an instant, never a time average." |
| `meds_fast_ark.f90:10-17` | "What used to live here and no longer does … moved by PR #120." | "ARK only. Shared machinery: meds_fast_frozen, meds_fast_be_stage, meds_column_state_ops." |
| `meds_soil_types.f90:35-41` | "Keeping mass here and enthalpy in the callers is what produced the two defects fixed in PR #81…" | "Soil and rainfall temperatures are inputs because this kernel owns the pond's enthalpy as well as its mass; split ownership lets the store drift." |

**What must not move.** Comments that explain a present mechanism stay even when they cite an
issue: the `theta_atm` default that floors `ustar` 44× (`meds_canopy_types.f90:212`, #97); the
transpiration corrector argument (`meds_fast_be_stage.f90:353`); why `site%deriv` must not be
read by the output layer (`meds_site_diag_types.f90:115`); why RK45 does not clamp the committed
state (`meds_fast_rk45.f90:231`). Only the "that was the first implementation and it was wrong"
clauses go.

**Stale facts inside comments, to fix in the same sweep:** `meds_leaf_gas_exchange.f90:6,457`
present `meds_plant_capi` as a current caller (it is `meds_capi_leaf`); `meds_fast_split` is
treated as current code in `meds_config_io.f90:265` and four tests; `meds_fast_types.f90:4,7` cite
`meds_column_dynamics`/`meds_column_derivs`; eight `[RETIRED]`-tagged citations of
`MEDS_INTEGRATOR_PARITY.md` omit its `archive/` path; ghost references to retired symbols
(`multilayer_roots` 11, `plant_water_tendency` 6, `snow_on` 6, `picard_*` 4, `GRP_PSI` 3) and 75
mentions of "split path" / "operator-split" across `src/` and `test/`.

### 5.4 What moves out of reader-facing pages

- `docs/science/numerical_scheme.md` §7 items 4–9 (lines 484–527): a defect-history essay
  (debug_error had no TOML reader; three RK45 instances fixed; pond now holds heat; snowmelt
  rerouted). Largest single candidate. Also §5 lines 229–234 ("tables that used to fill this
  section were retired") and §5a lines 252–256.
- `examples/example_biophysics/README.md` lines 297–370 "Three bugs this example found" (#139,
  #140, soil-energy time-level split) and lines 166–177 on why dt_fast changed. Keep the
  diagnostic pattern in one paragraph; the stories go to the changelog.
- `examples/example_demography/README.md` lines 68–78: the golden recapture and PR #137 offset
  table.
- `src/forcing/README.md` lines 86–92: the silently-falling-back recycle classifier that PR #69
  replaced (366 d 22 h span, ~10 h offset after 29 years).
- `src/fast_dynamics/plant/README.md`: façade deletion story, stale `src/allometry/` and
  `src/demography/` paths, "not wired yet", phenology described as tri-state.
- `src/slow_dynamics/soil/README.md`: "Note (module reorg)" cites `src/biophysics/`; its "reserved
  follow-ups" list P3 items that shipped in PR #64.
- `src/fast_dynamics/README.md`: still describes "the split sweep" (lines 12, 55, 66), cites
  `test_picard_coupling.f90` (does not exist, line 100), has a garbled coupling sentence at lines
  94–95, and `shared/functions`/`shared/state` paths at 79–85.
- `python/README.md` struck-through roadmap rungs (lines 70–77); `docs/README.md:26` "(formerly
  the top-level archive/)"; `python/meds/_libmeds.py:3,25-26` docstring history.
- CLAUDE.md: the 40 passages (default flips with dates, rename histories, "deleted in step 9", the
  `build_fast_context` cleanup story at lines 516–523, the step-10 rename list at 733–738).

### 5.5 The rule going forward

Add four lines to CLAUDE.md as a PR checklist:

1. A source comment states present-tense rationale and may cite an issue number or a science
   section. It does not say what the code used to do.
2. What changed and when goes in `CHANGELOG.md`, in the same PR.
3. What is deferred goes in `docs/ROADMAP.md` with an issue number, not in a "deferred / MVP /
   placeholder" comment.
4. A plan document whose last item ships gets its tombstone and moves to `docs/dev_plans/archive/`
   in the same PR.

---

## 6. Suggested sequence

| PR | Scope | Check |
|---|---|---|
| 0 | Fix the `soil_carbon_on` loader default and its two comments; decide `lwdown_source = "synthesize"`. | Root TOML run now reports the soil-carbon seam line. |
| 1 | Create `CHANGELOG.md` (seeded from the 118 PR titles + §1 facts) and `docs/ROADMAP.md` (§1.9). Retitle and complete `docs/dev_plans/README.md`. Decide the gitignored CAPI doc. | Every open item in §1 has a roadmap line. |
| 2 | Archive move: tombstones on 30 plans; three corrected headers; the shared banner rewritten in nine files; 24 path-prefixed citations and three doubled script paths fixed; `meds_config_july.toml:29` tagged. Trim `MEDS_CODE_STRUCTURE_DESIGN` to reference + roadmap, §10.2 as its own record. Write the three missing science pages, or defer those three moves. | `grep -rn "docs/dev_plans/" src test scripts` resolves every path. |
| 3 | README rewrite (§2); new `src/README.md` (§3), `post_proc/README.md`, `docs/configuration.md`, `docs/building.md`; refresh the four folder READMEs and the three stale example READMEs; re-pin or remove the README's carbon numbers. | No retired path (`src/plant/`, `src/biophysics/`, `src/driver/`, `src/shared/`, `src/demography/`, `src/allometry/`) in any tracked `.md` outside `archive/`. |
| 4 | CLAUDE.md split (§4.2): ~150-line CLAUDE.md, five `.claude/rules/*.md`, `CLAUDE.local.md`; delete the stale claims in §4.1; move the ED2 map to the comparison page. | `wc -l CLAUDE.md` < 200; every module name in it exists. |
| 5 | Comment sweep (§5.3), one PR per directory group: `fast_dynamics/numerics`; rest of `fast_dynamics`; `state` + `slow_dynamics`; `io` + `config`; `test`. | Comments only, so each PR is byte-identical on both compilers by construction, which is also the check. |
