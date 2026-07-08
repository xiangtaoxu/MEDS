# MEDS High-Priority Bugfix Implementation Plan

**Date:** 2026-07-07
**Status:** DESIGN / PLAN — report-only. No source is modified by this document.
**Source:** Assembled from the 7 confirmed fix designs answering the 8 critical/high
physical-process findings of the MEDS code review (`archive/MEDS_CODE_REVIEW_2026-07-06.md`, §2.1).
**Author:** Lead engineer (MEDS).

---

## 1. Scope & Bug Status

This plan covers the 8 critical/high physical-process bugs from the §2.1 review. Two of the
findings were double-checked against the ED2 reference and the running code; their verdicts are
recorded below. Every item here is a **design**: implementers must still write the code, tests,
and run both back ends.

| Bug | Title (abbrev.) | Severity | Status | Effort |
|-----|-----------------|----------|--------|--------|
| **1+3** | Fast biophysics loop never wired into `meds_main` (`gpp_accum=0` ⇒ zero carbon growth) + maintenance respiration omitted from carbon-mode NPP | Critical | **CONFIRMED** (coupled, one fix) | Medium-Large |
| **2** | Cohort fusion/fission reset carbon-mode prognostic pools (storage buffer destroyed) | High | **CONFIRMED** | Medium |
| **4** | C4 PEP/CO2-limited slope `kp_eff` has no temperature response | High | **CONFIRMED** | Low |
| **5** | TOML parsers silently swallow `iostat` errors; malformed required values pass as zeros/defaults | High | **CONFIRMED** | Small |
| **6** | Soil interior heat advection upwind/sign at `meds_column_energy.f90:82` | none | **REFUTED** — code is correct (docs only) | ~none |
| **7** | Dunne saturation-excess runoff spuriously large under free-drain/bedrock (`z_wt` pinned at column bottom) | High | **CONFIRMED** | Low |
| **8** | atm↔CAS scalar conductance omits the Monin-Obukhov `temp1/temp2` (c3) profile factor in the reusable CAS-twin kernels | High | **CONFIRMED** (kernels + docs; live driver already correct) | Small-Medium |

**Net:** 7 real fixes (bug1+3 counts once) + 1 documentation-only clarification (bug6, refuted).

---

## 2. Recommended Sequencing

The only hard cross-fix dependency is **bug2 depends on bug1+3**: bug1+3 introduces three new
per-cohort **extensive** accumulator fields (`leaf_resp_accum`, `stem_resp_accum`,
`root_resp_accum`) that the fusion routine `fuse_2_cohorts` enumerates *by hand*; if bug2's
carbon-mode fuse/split conservation lands first, those fields do not yet exist, and if bug1+3
lands after bug2 ships, the new fields will be silently dropped on fusion (an `nplant`-weighted
merge that was never added). So **bug1+3 lands first, bug2 immediately after**, sharing the same
`cohort_block` field additions and the same fusefiss patch sites.

Everything else is independent and can proceed in parallel: bug4 (single-kernel leaf physics),
bug5 (self-contained IO library), bug7 (single hydrology branch), bug8 (kernel signatures + docs),
and bug6 (docs only).

**Dependency list**

```
bug1+3  ──▶ bug2        (bug2 must conserve bug1+3's new accumulators on fuse/split)
bug4     (independent)
bug5     (independent)
bug6     (independent — docs only)
bug7     (independent)
bug8     (independent)
```

**Suggested order**

1. **bug1+3** — fast-slow carbon coupling + maintenance respiration. Adds the shared `cohort_block`
   accumulators (see §3). Largest fix; unblocks bug2.
2. **bug2** — carbon-mode fuse/fission conservation. Reuses §3's field list; extends the fuse
   `nplant`-weighted merge to the new accumulators and the four carbon pools.
3. **bug4 / bug5 / bug7 / bug8** — parallel, no shared state. Land in any order.
4. **bug6** — documentation clarification only; land whenever convenient.

**Rationale for landing bug1+3 first:** it is the only fix that (a) changes the state SoA
(`cohort_block`) and (b) creates new extensive per-ground fields that a second fix (bug2) must be
aware of. Sequencing it first collapses the two field-plumbing efforts into one lockstep audit,
avoiding the ED2 "forgot to reallocate a new cohort field" class of bug twice.

---

## 3. Shared Infrastructure Changes (do these once)

### 3.1 New per-cohort accumulator fields on `cohort_block` (bug1+3, also touched by bug2)

Three new accumulators added to `cohort_block` in
`/home/xiangtao/projects/MEDS/src/state/meds_demography_types.f90` next to `gpp_accum` (~line 84).
All are `[kgC/plant]`, reset each slow step, written by the fast loop, mirroring `gpp_accum`. Keep
`gpp_accum` **gross**; track loss terms separately (matches ED2).

| Field | ED2 analogue | Units |
|-------|--------------|-------|
| `leaf_resp_accum` | `today_leaf_resp` (leaf dark Rd) | kgC/plant |
| `stem_resp_accum` | `today_stem_resp` | kgC/plant |
| `root_resp_accum` | `today_root_resp` | kgC/plant |

**Lockstep routines that MUST be updated for EACH field (single-place invariant — the ED2
"forgot to reallocate" hazard).** Every one of these is a hand-enumerated site in
`meds_demography_types.f90` unless noted:

1. Type declaration (add the 3 `real(wp), allocatable :: X(:)`).
2. `cohort_alloc` — allocate + init `0.0` (beside line 179-180).
3. `site_free` — deallocate list (line 160).
4. `cohort_ensure_capacity` — copy block (after line 242).
5. `move_alloc_block` — (after line 277).
6. `cohort_reorder` — permute (after line 344).
7. `copy_cohort_slot` — (after line 394); **also covers `split_cohorts`**.
8. **Fusion weighting** in `/home/xiangtao/projects/MEDS/src/demography/meds_demography_fusefiss.f90`
   at the cohort merge (~line 242-243, beside the existing `gpp_accum` `nplant`-weighted blend,
   **before** `nplant(recc)=ntot`): add the identical `nplant`-weighted merge for all 3:
   `X(recc)=(nr*X(recc)+nd*X(donc))/ntot`.

**Do NOT** add these to `gather_pft_params`, `set_cohort_size`, or `set_cohort_size_from_carbon`
(they are fluxes, not PFT params or re-derived geometry). No per-patch field is added, so
`sort_patches`/`patch_compact` are untouched.

> **Invariant:** these accumulators are **extensive per-ground fluxes**, so fusion weights them by
> `nplant` (never leaf-area). Leaf-area weighting would silently break carbon conservation across
> fusion. Bug2's new carbon-pool conservation guard error-stops if a newly added extensive field is
> left unweighted, surfacing an omission at runtime.

### 3.2 New `cas_atm_forcing_t.temp1` / `.temp2` fields (bug8)

Add `real(wp) :: temp1` and `real(wp) :: temp2` (dimensionless Monin-Obukhov scalar-transfer
coefficients for heat and vapour/scalar) to `cas_atm_forcing_t` in
`/home/xiangtao/projects/MEDS/src/biophysics/meds_biophysics_types.f90` (~line 244-250). The
producer already computes them: `aero_out_t%temp1/temp2` (line 365). Prefer **no safe default**
(or a sentinel) so a caller that forgets to fill them fails loudly rather than silently
reinstating `c3=1`. If a default is required for standalone tests, default `1.0` but document that
the caller **must** copy `aero_out%temp1/temp2` each fast substep.

### 3.3 New PFT traits / TOML keys

**None required by the recommended paths.** Specifically:

- **bug4** reuses the existing Vcmax activation-energy set (`ea_vcmax/hd_vcmax/ds_vcmax`, already
  in `leaf_photo_params_t`). No new trait, no TOML change. (An optional dedicated `ea_kp` path is
  documented but **not recommended** — see §4.4.)
- **bug1+3** needs a `build_fast_context` builder assembling soil/resp/hydro/co2 params from
  `cfg%pft` or documented constants (MVP constant forcing; no `column_config` TOML loader exists in
  production yet). This is new *code*, not new persistent state.
- **bug5** adds only an optional private helper `toml_parse_error`; no new state/type/DAG change.

---

## 4. Per-Bug Fix Sections

---

### 4.1 Bug 1+3 — Fast biophysics loop never wired; maintenance respiration omitted from NPP

- **Severity:** Critical · **Effort:** Medium-Large · **Status:** CONFIRMED (both bugs, one coupled fix)
- **Files:** `src/driver/meds_main.f90`, `src/driver/meds_stepper.f90`, `src/driver/meds_fast_loop.f90`,
  `src/driver/meds_column_dynamics.f90`, `src/driver/meds_vegetation_dynamics.f90`,
  `src/state/meds_demography_types.f90`, `src/demography/meds_demography_fusefiss.f90`,
  `test/test_fast_loop.f90`

**Root cause (confirmed).**

- **BUG1 (fast loop never invoked):** `meds_main.f90:124` calls
  `advance_one_step(site,cfg,is_new_month,is_new_year)` **without** the optional `fast_ctx`.
  `meds_stepper.f90:37` runs `run_fast_biophysics` only when
  `(cfg%fast_biophysics_on .AND. present(fast_ctx))`; with no ctx the fast loop is **silently
  skipped**. `run_fast_biophysics` (`meds_fast_loop.f90:164`) is the ONLY writer of
  `cohort%gpp_accum`. So with `growth_source=GS_CARBON` and `fast_biophysics_on=.true.`,
  `carbon_growth` (`meds_vegetation_dynamics.f90:98-99`) reads `a_carbon = cohort%gpp_accum(j) = 0`
  ⇒ `net_carbon<=0` ⇒ `plant_carbon_allocation` grows nothing. Secondary: `init_fast_reservoirs`
  (`meds_fast_loop.f90:53`) is called ONLY in tests, so even a passed ctx would have zero cas/soil
  reservoirs.
- **BUG3 (maintenance resp dropped):** `carbon_growth` (`meds_vegetation_dynamics.f90:103`) sets
  `env%net_carbon = a_carbon - growth_respiration(a_carbon, growth_resp_factor)` — only GROWTH
  respiration is deducted. Leaf dark resp (`lf%rd`), stem maintenance (`wf%stem_resp`), fine-root
  maintenance (`rf%root_resp`) ARE computed inside `column_fast_step`
  (`meds_column_dynamics.f90:227/258/260` as `ra_leaf/ra_stem/ra_root`) but only for the CAS CO2
  twin (line 268); never accumulated per-cohort nor handed to the slow seam. NPP overestimated by
  the maintenance term. ED2 `growth_balive.f90:857-878` confirms the intended math.

**Fix (numbered).**

**(a) State: add the 3 accumulators** to `cohort_block` and update every lockstep site — see §3.1
(items 1-8, including the fusefiss `nplant`-weighted merge).

**(b) Producer: emit per-plant resp fluxes.** In `column_fast_step`
(`meds_column_dynamics.f90`, cohort loop 216-263) add 3 OPTIONAL `intent(out)` args mirroring
`gpp_coh(:)`: `leaf_resp_coh(:)`, `stem_resp_coh(:)`, `root_resp_coh(:)`, each `[µmol CO2/plant/s]`,
populated per-plant, guarded by `present()`:

```fortran
if (present(leaf_resp_coh)) leaf_resp_coh(i) = lf%rd * coh%leaf_area(i)   ! per plant
if (present(stem_resp_coh)) stem_resp_coh(i) = wf%stem_resp               ! already per-plant
if (present(root_resp_coh)) root_resp_coh(i) = rf%root_resp               ! already per-plant
```

> `ra_leaf` at line 232 is `*nplant` (per-ground) — do **NOT** reuse it. Multiply `lf%rd` by
> `leaf_area` to get per-plant. Stem/root kernels already return per-plant.

**(c) Integrator: accumulate.** In `run_fast_biophysics` (`meds_fast_loop.f90`): extend the
pre-window reset (line 104) to zero the 3 new accumulators; allocate the 3 `*_resp_coh` named
arrays alongside `gpp_coh` (line 113); pass them into `column_fast_step` (line 159); in the substep
integrator (161-165) integrate each with the SAME conversion as `gpp`:

```fortran
X_accum(i) = X_accum(i) + X_coh(j) * cfg%dt_fast * umol_2_kgC   ! [kgC/plant]
```

Keep `gpp_accum` gross; do **not** net here.

**(d) Slow seam: subtract maintenance.** Rewrite `carbon_growth`
(`meds_vegetation_dynamics.f90:95-103`):

```fortran
if (cfg%fast_biophysics_on) then
   gross_gpp  = cohort%gpp_accum(j)
   resp_maint = leaf_resp_accum(j) + stem_resp_accum(j) + root_resp_accum(j)
else                                   ! fast-off stub (no mechanistic maintenance source)
   gross_gpp  = cfg%gpp_ref * leaf_area * dt_yr
   resp_maint = 0.0_wp
end if
npp_after_maint = gross_gpp - resp_maint
env%net_carbon  = npp_after_maint - growth_respiration(npp_after_maint, pft%growth_resp_factor(pf))
```

Downstream `get_plant_flux_slow → npp/npp_repro → carbon_vital_rates →
update_demography/apply_carbon_npp` is UNCHANGED. `net_carbon` may go negative ⇒ already handled by
`plant_carbon_allocation`'s deficit/starving path.

**(e) Harden the gate + wire `meds_main`.** In `meds_stepper.f90` (37-39) replace the silent gate:

```fortran
if (cfg%fast_biophysics_on) then
   if (.not. present(fast_ctx)) &
      error stop 'advance_one_step: fast_biophysics_on but no fast_context supplied'
   call run_fast_biophysics(site, fast_ctx, cfg)
end if
```

In `meds_main.f90`: declare `type(fast_context_t) :: fast_ctx`; if `cfg%fast_biophysics_on`, call
`build_fast_context(cfg,fast_ctx)` + `init_fast_reservoirs(site,fast_ctx)` **once** after the
community build (after line 85), and pass `fast_ctx` into `advance_one_step` (line 124). The ctx
MUST be owned by `meds_main` and passed as an arg — the stepper cannot cache it (hidden global state
forbidden).

**(f) `build_fast_context` (new).** Small builder (in `meds_fast_loop.f90` or a `meds_main`
contains-routine) assembling `ctx%ccfg` (soil via `build_soil_params/build_soil_thermal`, wood/root
resp params from `cfg%pft`, hydro params, co2 opts) + reference met from `cfg`. MVP constant forcing
since production has no `column_config` TOML loader yet.

**(g) Fusion coupling (bug2 seam).** The 3 accumulators are extensive per-ground fluxes ⇒ add them
to the fuse `nplant`-weighted merge list (§3.1 item 8). Split needs no change — `copy_cohort_slot`
already carries every per-plant field.

**(h) Empirical bit-identical guarantee.** Default `cfg` keeps `fast_biophysics_on=.false.` +
`growth_source=GS_EMPIRICAL`; no ctx is built, the guard is never reached, and the 3 new arrays are
allocated + zeroed but never read. Verify empirical spin-up output unchanged.

**(i) Test update.** `test/test_fast_loop.f90` gate-off sub-case (86-93) currently sets
`fast_biophysics_on=.true.` then calls `advance_one_step` with no ctx expecting a silent skip — now
error-stops under the guard. Set `fast_biophysics_on=.false.` for that sub-case (or assert the
abort). Gate-ON (95-99) and end-to-end handoff (106-126) stay valid.

**Corrected math / units.** Per-plant resp accumulation mirrors `gpp_accum`: `gpp_coh`,
`leaf_resp_coh`, `stem_resp_coh`, `root_resp_coh` are `[µmol CO2/plant/s]`; each substep
`accum += flux * cfg%dt_fast[s] * umol_2_kgC[kgC/µmol]`, over `n_fast_per_slow` substeps = one slow
step ⇒ `[kgC/plant]`.

```
net_carbon = gpp_accum
           - (leaf_resp_accum + stem_resp_accum + root_resp_accum)
           - Rg,      Rg = growth_resp_factor * max(0, gpp_accum - resp_maint)
```

ED2 (`growth_balive.f90:857-878`): `metab_resp = today_leaf_resp + today_root_resp + today_stem_resp;
metnpp = today_gpp - metab_resp`. ED2's `f_unitconv = umol_2_kgC*day_sec/nplant` because ED2
accumulates per-GROUND µmol then divides by `nplant`; MEDS accumulates **per-PLANT** directly
(`leaf_resp_coh` carries `leaf_area`; stem/root already per-plant), so **no `/nplant`**. Fusion
blend is `nplant`-weighted: `X_new = (nr*Xr + nd*Xd)/(nr+nd)`.

**ED2 reference.** `dynamics/growth_balive.f90:857-878` (`metab_resp`, `metnpp`, `f_unitconv` at
:863); `update_carbon_balances :1036-1161`; `stem_resp_driv.f90` for `today_stem_resp`.

**Test plan.**
- Empirical bit-identical: default spin-up before/after ⇒ identical `-D-output.nc` + `print_summary`.
- Carbon+fast end-to-end (extend `test_fast_loop.f90` case 4): fast on, `GS_CARBON`, `gpp_ref=0` ⇒
  assert `gpp_accum>0` AND `leaf/stem/root_resp_accum>0`, and net gain = `gpp_accum - resp_maint - Rg`
  (maintenance visibly reduces NPP vs pre-fix).
- Guard trap: `fast_biophysics_on=.true.` + `advance_one_step` with NO `fast_ctx` must error stop.
- Lockstep integrity: force fusion and split on a carbon-mode multi-cohort site; assert the 3
  accumulators are `nplant`-weighted on fusion / copied on split, site carbon still conserves.
- Reset semantics: two consecutive slow steps; assert each accumulator resets at
  `run_fast_biophysics` entry (no cross-step leakage).
- Build `ifx -stand f18 -check all` AND `nvfortran` multicore (green ifx not sufficient; new optional
  array out-args risk the nvfortran array-temp miscompile — bind to **named** arrays).

**Risks.**
- `build_fast_context` is new work with no production `column_config` TOML loader; resp factors
  (`stem_resp_factor25`, `root_resp_factor25`, `agf_bs`, `ea/hd/ds`) and soil params must come from
  `cfg%pft` or documented constants, else maintenance magnitude is arbitrary.
- Hard guard changes `test_fast_loop.f90` gate-off semantics (must update or it error-stops).
- `net_carbon` can now go negative ⇒ exercises `plant_carbon_allocation` storage-drawdown/starving
  path more often (physically correct; tissue destruction is a future PR).
- nvfortran portability: new optional array out-args + `*_resp_coh` temporaries must be **named**
  arrays, never function-result temporaries.
- Fusion weighting MUST use `nplant` (extensive per-ground); leaf-area weights break conservation.
- Fast-off stub keeps `resp_maint=0` by design; document so the stub is not mistaken for
  maintenance-inclusive NPP.

---

### 4.2 Bug 2 — Cohort fusion & fission destroy carbon-mode prognostic pools

- **Severity:** High · **Effort:** Medium · **Status:** CONFIRMED · **Depends on:** bug1+3
- **File:** `src/demography/meds_demography_fusefiss.f90`

**Root cause (confirmed).** In `meds_demography_fusefiss.f90`, `fuse_2_cohorts` (222-256) conserves
`nplant` + total AGB, re-derives `dbh` from the EMPIRICAL `agb_to_dbh` (:249), then calls
`set_cohort_size` (:250). `set_cohort_size` (`meds_demography_types.f90:429-432`) OVERWRITES all
four carbon pools (`leaf_carbon`, `fineroot_carbon`, `wood_carbon`, `nonstructural_carbon`) with
on-allometry values. `split_cohorts` does the same at :291. The pools are never `nplant`-weight
merged the way agb/gpp_accum/leaf_temp/psi are.

In `GS_CARBON` mode `wood_carbon` is the prognostic size anchor: because
`agb = p_aboveground_frac*wood_carbon` and `dbh=wood_to_dbh` is the exact inverse, the
`agb_to_dbh → set_cohort_size` round-trip recovers `wood_carbon` and AGB within tolerance, so the
AGB guard passes and `wood_carbon` survives. The DAMAGE is to the OTHER THREE pools:
`set_cohort_size` resets `leaf_carbon = leaf_area/sla` (destroys phenological/NPP deviations),
`fineroot_carbon = root_to_leaf_ratio*leaf_carbon`, and — most seriously —
`nonstructural_carbon = storage_cushion*leaf_carbon`, obliterating the prognostic storage buffer
whose depletion is the carbon-starvation signal feeding growth-dependent (Camac) mortality. The
carbon-aware inverse `set_cohort_size_from_carbon` (`meds_demography_types.f90:441-450`) takes the
pools as INPUTS (does not overwrite) but is unused by fuse/split.

> **Energy/hcap is NOT part of this bug.** MEDS has no explicit `hcap`/`leaf_energy` field. Leaf
> heat content is conserved *implicitly* by the existing leaf-area-weighted `leaf_temp`/`psi` merge
> (:237-241) — already correct; preserve it.

**Fix (numbered).**

1. **Module wiring:** add `set_cohort_size_from_carbon` to the `meds_demography_types` `only`
   import (26-28); add `GS_CARBON` to the `meds_config` import (already `used`). No new state/config
   fields — both already exist.
2. **Thread `growth_source` into `fuse_2_cohorts`:** signature becomes
   `fuse_2_cohorts(site, recc, donc, conservation_tol, growth_source)`; the single caller `fuse_pass`
   (:202) passes `cfg%growth_source`. `split_cohorts`/`new_fuse_cohorts` already receive `cfg`.
3. **FUSE Step 1 (conserve FIRST, ~228-247):** keep `nr,nd,ntot`, the leaf-area-weighted
   `leaf_temp`/`psi` merge (:237-241), and the `nplant`-weighted `gpp_accum` merge (:243) EXACTLY.
   **Add** an `nplant`-weighted merge of the four carbon pools using the OLD `nplant` weights:
   ```fortran
   leaf_carbon(recc)          = (nr*leaf_carbon(recc)          + nd*leaf_carbon(donc))          / ntot
   fineroot_carbon(recc)      = (nr*fineroot_carbon(recc)      + nd*fineroot_carbon(donc))      / ntot
   wood_carbon(recc)          = (nr*wood_carbon(recc)          + nd*wood_carbon(donc))          / ntot
   nonstructural_carbon(recc) = (nr*nonstructural_carbon(recc) + nd*nonstructural_carbon(donc)) / ntot
   ```
   Then `nplant(recc)=ntot` and (empirical only) `agb(recc)=agb_tot/ntot` as today.
4. **FUSE Step 2 (derive geometry SECOND, branch):**
   ```fortran
   if (growth_source == GS_CARBON) then
      call set_cohort_size_from_carbon(site%cohort, recc)   ! dbh=wood_to_dbh(wood_carbon); pools untouched
   else                                                     ! GS_EMPIRICAL — unchanged
      dbh(recc) = agb_to_dbh(agb(recc), ...)
      call set_cohort_size(site%cohort, recc)               ! legitimately re-derives on-allometry diagnostics
   end if
   ```
5. **FUSE Step 3 (guards):** keep the existing AGB-density guard (:252-254) for BOTH modes. **Add**,
   guarded by `GS_CARBON`, a carbon-pool conservation check: capture the four pre-merge
   donor+survivor totals BEFORE overwriting, then per pool verify
   `|ntot*pool(recc) - total_before| <= conservation_tol*max(total_before, tiny_num)`, else
   `error stop 'fuse_2_cohorts: carbon-pool conservation violated'`.
6. **SPLIT — EMPIRICAL branch:** when `cfg%growth_source /= GS_CARBON` leave current arithmetic
   bit-identical (renorm from `p_agb=agb_c2*(2+b2Ht)`, `dbh=d0*(1±eps)*renorm`, `set_cohort_size`,
   AGB guard :306-307). Add a comment documenting the `hgt_max` interaction (see risks).
7. **SPLIT — CARBON branch (new, exactly conservative, `hgt_max`-immune):** do NOT perturb `dbh`.
   Perturb the conserved anchor `wood_carbon` symmetrically:
   ```fortran
   wc0 = wood_carbon(i)
   wood_carbon(i) = wc0 * (1.0_wp + eps)          ! '+eps' daughter kept in slot i
   call copy_cohort_slot(cohort, m, i)            ! copies ALL per-plant pools + accumulators
   wood_carbon(m) = wc0 * (1.0_wp - eps)          ! '-eps' daughter
   ! leaf_carbon / fineroot_carbon / nonstructural_carbon UNCHANGED in both daughters
   call set_cohort_size_from_carbon(cohort, i)
   call set_cohort_size_from_carbon(cohort, m)
   ```
   `renorm` is exactly 1 for `wood_carbon` (`0.5n·wc0(1+eps)+0.5n·wc0(1-eps)=n·wc0`), independent of
   the `hgt_max` cap ⇒ conserves ALL four pools AND agb exactly.
8. **SPLIT guards (carbon):** extend the post-loop check to assert total `wood_carbon` (hence AGB)
   conserved via the sum-before/sum-after pattern already at :281/:299; optionally assert the other
   three pool totals. Keep fresh-global-id stamping (:303-305) and `rebuild_csr/sort_cohorts` as-is.
9. **bug1+3 coupling:** any NEW extensive per-plant accumulator from bug1+3 MUST be added to the
   FUSE Step-1 `nplant`-weighted merge list (fuse enumerates fields by hand). SPLIT needs no change
   (`copy_cohort_slot` is centralized).

**Corrected math / units.** Carbon pools are `[kgC/plant]` (intensive per-plant); conserved
quantity is the per-ground total `nplant*pool [kgC/m²]`; merge is `pool_f=(nr*pool_r+nd*pool_d)/(nr+nd)`
— exactly ED2's `rnplant/dnplant` weighting of `bleaf/broot/bsapwood/bdead/bstorage`. ED2 derives
`dbh/hite` from the STRUCTURAL pool (`bd2dbh(bdead)`) AFTER pool fusion; MEDS's carbon-mode analogue
is `set_cohort_size_from_carbon` deriving `dbh=wood_to_dbh(wood_carbon)`. Since
`agb=p_aboveground_frac*wood_carbon`, conserving `nplant*wood_carbon` conserves `nplant*agb`
automatically. The leaf-area-weighted `leaf_temp` merge conserves `Σ(hcap·T)` (constant-cp,
no-phase-change equivalent of ED2's add-then-invert). The pool weighting MUST run before the
geometry refresh (it uses pre-fusion `leaf_area`).

**ED2 reference.** `utils/fuse_fiss_utils.f90` `fuse_2_cohorts` (1986-4428): `nplant`-weight fusion
of `bdead*/bleaf/broot/bsapwood*/bstorage` via `rnplant/dnplant` (:2060-2159); `dbh/hite` from the
structural pool AFTER (`bd2dbh/dbh2h`, :2178-2189); leaf/wood energy/water/hcap added extensive then
temperature recovered via `uextcm2tl` (:2229-2234, 2655-2662); LAI-weighted psi (:2439-2455). Split
analogue ~:1689-1960.

**Test plan.**
- Carbon-mode fusion conservation: 2-cohort same-PFT patch (`GS_CARBON`) with `nonstructural_carbon`
  deliberately OFF allometry, force fuse; assert `Σ(nplant*pool)` for all four pools conserved to
  `conservation_tol` AND `nonstructural_carbon(recc) != storage_cushion*leaf_carbon` (storage
  carried, not reset). Pre-fix fails this.
- Carbon-mode fission conservation: cohort with `LAI>cohort_lai_cap` + off-allometry storage, call
  `split_cohorts`; assert total `nplant*wood_carbon`, total agb, and the other three pool totals
  conserved, and both daughters' per-plant storage equals parent's (unperturbed).
- `hgt_max`-immunity: cohort at/above `hgt_max` with high LAI; carbon-mode split conserves
  `wood_carbon/agb` exactly where the empirical `dbh`-renorm would drift.
- Empirical regression: run fusion/spin-up under `GS_EMPIRICAL` before/after ⇒ bit-identical
  `dbh/agb/nplant` trajectories.
- Guard trip: inject imbalance ⇒ new carbon-pool conservation error stop fires; existing AGB guard
  still fires on an AGB violation.
- Build both back ends: `ifx -stand f18 -check all` (full ctest) + nvfortran multicore.

**Risks.**
- Behaviour change confined to `GS_CARBON`; empirical path untouched (verify bit-identical).
- `fuse_2_cohorts` enumerates every merged field by hand — a second maintenance site alongside
  `cohort_reorder/copy_cohort_slot`; the new conservation guard surfaces omissions at runtime.
- `hgt_max` interaction in the EMPIRICAL split (approximate AGB conservation absorbed by the 1%
  guard); carbon-mode split is exact and immune. Follow-up: split empirical daughters on agb
  directly.
- Carbon-mode split leaves `leaf_carbon` unperturbed while `wood_carbon` is `±eps` (daughters
  slightly off each other's allometry) — intentional; pools are prognostic.
- `conservation_tol` applied to totals via `max(total,tiny_num)`; capture pre-merge totals BEFORE
  the survivor slot is overwritten (in-place overwrite trap).

---

### 4.3 Bug 4 — C4 PEP/CO2-limited slope `kp_eff` has no temperature response

- **Severity:** High · **Effort:** Low (one-line kernel change; no new trait/TOML; add one test)
  · **Status:** CONFIRMED
- **Files:** `src/plant/meds_leaf_gas_exchange.f90`, `test/test_leaf_physiology.f90`

**Root cause (confirmed).** In `meds_leaf_gas_exchange.f90`, `solve_leaf_gas_exchange`
temperature-scales every C4 biochemical term except the PEP/CO2 slope: `vcmax/jmax/tpu/rd` go
through `temp_response`/`arrhenius_scale` (209-212), but line 227 builds
`kp_eff = p%kp25 * pressure / p_std` — a pressure correction ONLY, frozen at 25 °C.
`assim_demand_c4` (:88) uses `ap = kp_eff*ci` as the C4 CO2/PEP-limited rate, so the PEPcase initial
slope never warms up or rolls off. ED2 (`farq_leuning.f90:1606`) sets `rho = klowco28 * vm`, i.e.
`kp` is tied to the FULLY temperature-corrected Vcmax; `klowco2` (`ed_params.f90:1222`) is a
temperature-INDEPENDENT ratio, so ED2's `kp` inherits 100% of Vcmax's T-dependence and the high-T
deactivation, keeping `kp/Vcmax` constant across temperature.

**Fix (recommended — reuse Vcmax's activation-energy set; ED2-faithful, zero new params).**

1. Replace `meds_leaf_gas_exchange.f90` line 227
   `kp_eff = p%kp25 * pressure / p_std` with:
   ```fortran
   kp_eff = temp_response(tresp, p%kp25, p%ea_vcmax, p%hd_vcmax, p%ds_vcmax, &
                          t_leaf) * pressure / p_std
   ```
2. No new imports/args: `temp_response` already imported (line 19); `tresp`, `t_leaf`, and
   `p%ea_vcmax/hd_vcmax/ds_vcmax` are all in scope. No change to `meds_pft_params`,
   `meds_plant_types`, or either TOML.
3. Keep the `pressure/p_std` correction exactly. `temp_response` returns `1.0` at
   `t_leaf=t_ref_photo` (25 °C), so `kp_eff` is byte-for-byte unchanged at reference temperature
   (existing 25 °C C4 tests stay green — pure regression-safe extension).
4. This gives `kp` the SAME treatment as Vcmax under whichever `tresp` is configured:
   `TRESP_ARRHENIUS` ⇒ pure Arrhenius rise (`kp/vcmax=const(T)`, the ED2 invariant);
   `TRESP_PEAKED` ⇒ the same high-T deactivation envelope (`hd_vcmax/ds_vcmax`). Do **NOT** add a
   separate envelope for `kp`.
5. Add a unit test (see test plan).

> **Alternative (NOT recommended):** decouple `kp`'s Q10 from Vcmax by adding a shared `ea_kp`
> (~52950 J/mol for an exact Collatz-1992 Q10=2) to `[leaf_physiology]` in `meds_config_main.toml`,
> a field in `leaf_photo_params_t`, and its loader (interface + capi + python `_ffi.py`). Breaks the
> ED2 `kp/Vcmax=const` relationship, adds a param, touches 4-5 files for a ~0.35 Q10 difference.

**Corrected math / units.**
```
kp_eff = temp_response(tresp, kp25, ea_vcmax, hd_vcmax, ds_vcmax, t_leaf) * pressure / p_std
Arrhenius:  f(T) = exp[ Ea/(R·Tref) · (1 - Tref/T) ],   Tref = 298.15 K
Peaked   :  × fH(Tref)/fH(T),   fH(T) = 1 + exp[(dS·T - Hd)/(R·T)]
```
Units: `kp25 [mol m⁻² s⁻¹] × f(T)[-] × pressure/p_std[-] = kp_eff [mol m⁻² s⁻¹]`; then
`ap = kp_eff·ci` (`ci [µmol mol⁻¹]`) gives `ap [µmol m⁻² s⁻¹]`, consistent with `vcmax/ac`. Effective
Q10 of `ea_vcmax=65330 J/mol` near 25 °C ≈ 2.35 ⇒ `kp_eff/vcmax` constant with T (the ED2 invariant).

**ED2 reference.** `dynamics/farq_leuning.f90:1606` (`rho = klowco28 * vm`); Vcmax temp response
:514-573; `klowco2` T-independent, `init/ed_params.f90:1222`. Physics: Collatz, Ribas-Carbo & Berry
1992.

**Test plan.**
- Regression: C4 PFT at `t_leaf=25 °C` ⇒ `kp_eff`/CO2-limited outcome unchanged (`temp_response(Tref)=1`).
- Warming (`TRESP_ARRHENIUS`): C4 leaf under CO2/PEP-limiting `Ci` at 15/25/35 °C ⇒ PEP-limited rate
  `ap` (limitation flag `LIM_C4_PEP`) strictly increasing with temperature.
- Invariant: under `TRESP_ARRHENIUS`, `kp_eff/vcmax` temperature-invariant.
- Deactivation: under `TRESP_PEAKED`, `kp_eff(45 °C) < kp_eff(35 °C)`.
- Build ifx (`-DMEDS_ENABLE_IO=OFF` Debug `-check all`) + nvfortran multicore; `ctest -R
  leaf_physiology`. Re-baseline any C4 assertion pinned at a non-25 °C temperature.

**Risks.**
- Reusing `ea_vcmax/hd_vcmax/ds_vcmax` couples `kp`'s response (incl. peaked roll-off) to Vcmax; a
  future C4 PFT wanting a distinct `kp` Q10 needs the dedicated-trait path.
- MEDS's peaked form models only high-T deactivation (no ED2 `tlow_fun`); `kp` inherits that
  MEDS-wide omission at cold leaf temperatures.
- Any test asserting absolute C4 assimilation at a non-25 °C temperature must be re-baselined.
- Effective Q10~2.35 slightly above Collatz's canonical 2; acceptable (keeps the ED2 invariant).

---

### 4.4 Bug 5 — TOML parsers silently swallow `iostat` errors

- **Severity:** High · **Effort:** Small (one self-contained IO library; no state/type/DAG change)
  · **Status:** CONFIRMED
- **Files:** `src/io/meds_toml.f90`, `src/io/meds_config_io.f90`

**Root cause (confirmed).** The getters (`meds_toml.f90` header 9-10) return a default when a key is
absent OR unparseable — contradicting the `meds_config_io` contract (every parameter REQUIRED, no
defaults, missing/unparseable = HARD ERROR). The presence-map (`keymiss_t`) only catches ABSENT
keys; a PRESENT-but-malformed key never reaches `note_missing`, so garbage flows through. Specifics:

1. `toml_real_array` (:187-198) does `read(...,iostat=ios)` but NEVER checks `ios`; a malformed
   token makes list-directed read stop early leaving trailing `tmp` at `0.0`, and
   `nout=count_tokens(s)` counts ALL tokens ⇒ OVER-reports length. `req_pa` (config_io:216-226) then
   sees `nout==npft` and accepts a vector with silent zeros.
2. `toml_logical` (:147-148) fixed-position `s(1:4)=='true'`/`s(1:5)=='false'`; any other spelling
   silently keeps default, and `'trueish'` wrongly matches.
3. `req_dur` (config_io:112-117) reads with `iostat` but only assigns `secs` when `ios==0`; a
   present-but-unparseable duration silently keeps the `86400 s` initializer.
4. PFT-count inference (config_io:347-351) sets `npft=nout` from `pft.wood_density` with no `MAXPFT`
   bound; `buf` is `MAXPFT(64)`-sized ⇒ >64 tokens is an OOB read/write inside `toml_real_array`.

No external callers of the getters exist outside `meds_toml.f90`/`meds_config_io.f90`, so an
internal error-stop policy is signature-safe.

**Fix (numbered).**

1. **Policy:** keep public getter signatures UNCHANGED (they still take `default`, still return it
   for key-ABSENT, which `req_*` never reaches). Change ONLY the present-but-unparseable path to a
   HARD ERROR. Scalar/array getters in `meds_toml` error stop INTERNALLY; config_io helpers that own
   their parsing and have `m` in scope (`req_dur`, `req_date`-style) route to `note_missing`. Do NOT
   add an `iostat`/`err` out-arg.
2. **`meds_toml.f90`:** add private helper near `find_key`:
   ```fortran
   subroutine toml_parse_error(key, raw, what)   ! what = 'integer'|'real'|'logical'|'real array'
      write(*,'(a)') ' meds_toml: cannot parse '//what//' for key "'//trim(key)//'"'
      write(*,'(a)') '   raw value: "'//trim(raw)//'"'
      error stop 'meds_toml: malformed configuration value'
   end subroutine
   ```
3. **`toml_int` (115-124) / `toml_real` (126-135):** replace `if (idx>0 .and. ios/=0) v=default`
   with:
   ```fortran
   idx = find_key(t, key);  if (idx == 0_ik) return          ! absent: keep v=default
   read(t%val(idx), *, iostat=ios) v
   if (ios /= 0) call toml_parse_error(key, t%val(idx), 'integer')   ! or 'real'
   ```
4. **`toml_logical` (137-149):** drop fixed-position compares. After `s = adjustl(t%val(idx))`:
   ```fortran
   select case (trim(s))
   case ('true','.true.','True','TRUE');   v = .true.
   case ('false','.false.','False','FALSE'); v = .false.
   case default; call toml_parse_error(key, t%val(idx), 'logical')
   end select
   ```
   Exact-match kills the `'trueish'` false-positive; unknown spelling is a hard error.
5. **`toml_real_array` (170-198):** rewrite so absent/present-malformed/bad-token/too-many are
   distinct:
   ```fortran
   idx = find_key(t, key);  if (idx == 0_ik) then; nout = 0_ik; return; end if   ! (a) ABSENT
   lb = index(s,'[');  rb = index(s,']', back=.true.)
   if (lb == 0 .or. rb <= lb) call toml_parse_error(key, t%val(idx), 'real array') ! (b) not a [..] literal
   ! commas -> spaces
   ntok = count_tokens(s)
   if (ntok == 0_ik) then; nout = 0_ik; return; end if                            ! (c) empty []
   if (ntok > size(out, kind=ik)) call toml_parse_error(key, t%val(idx), 'real array') ! (d) buffer/OOB bound
   read(s, *, iostat=ios) (tmp(k), k = 1_ik, ntok)
   if (ios /= 0) call toml_parse_error(key, t%val(idx), 'real array')              ! (e) every token must parse
   out(1:ntok) = tmp(1:ntok);  nout = ntok                                         ! (f) TRUE parsed count
   ```
   Now `nout` is never over-reported; `req_pa`'s `nout==npft` remains the length authority. Delete
   the wrong inline comment at line 192.
6. **`req_dur` (config_io 88-118):** change initializer `secs=86400.0_wp` → `secs=0.0_wp` (dead
   fallback). After the read: `if (ios==0) then; secs=v*mult; else; call note_missing(m,key); end
   if`. Treat `n==0` (empty string) as `note_missing`. Present-but-unparseable duration becomes an
   aggregated hard error, consistent with `req_date`.
7. **PFT-count inference (config_io 347-351):** after the existing `if (nout<1_ik) error stop ...`,
   add `if (nout > MAXPFT) error stop 'meds_config: pft.wood_density lists more than MAXPFT PFTs in
   '//trim(cfg%pft_config)` before `npft=nout`. (The new `toml_real_array` `ntok>size(out)` guard
   fires first for a `MAXPFT`-sized buf; this gives the clear domain message.)
8. **Docs:** update the `meds_toml.f90` header (9-10) to state a getter returns the default ONLY
   when the key is ABSENT; a present-but-unparseable value is a hard error.
9. **Flag only (out of scope, do NOT change):** `toml_string` (:162-164) computes
   `v = s(2:index(s(2:),'"'))` — `index()` is relative to `s(2:)` but used as an absolute upper
   bound, mis-slicing quoted strings when the opening quote is not at position 1. Recommend a
   follow-up ticket; it is not an `iostat`-swallow issue.

**Corrected math / units.** No numeric math changes (parsing/validation only). The `req_dur`
multiplier table is unchanged and authoritative: `s/S=1, m/M=60, h/H=3600, d/D=86400, w/W=604800 s`;
a bare number = seconds. `iostat` variables may keep `integer(ik)` kind.

**ED2 reference.** N/A — MEDS-only config I/O (`meds_toml` is a bespoke minimal TOML reader with no
ED2 analogue). The relevant contract is the MEDS "no hard-coded parameters / every parameter
required from TOML" rule (`CLAUDE.md` + `meds_config_io` header).

**Test plan.**
- Build the netCDF-free strict-bounds Debug tree (ifx, `-DMEDS_ENABLE_IO=OFF`, `-check all`); full
  suite stays green (well-formed configs are a no-op).
- Also build nvfortran multicore (`-DMEDS_GPU=multicore`).
- Malformed PFT trait token: `pft.sla = [9.0, x, 11.0]` ⇒ aborts via `toml_parse_error` naming
  `pft.sla` and printing the raw text, not a zeros run.
- Short/over-reported array: `pft.vcmax25` fewer entries than `npft` ⇒ aggregated missing-key +
  `load_meds_config` error stops; trailing garbage token ⇒ hard parse error.
- Logical spellings: `.true.`/`True` parse; `yes`/`ture` error stop with key+raw.
- Duration: `run.dt_slow='oned'` ⇒ reported as missing/invalid required key and aborts; `'1800s'`,
  `'12h'`, `'7d'`, `'3600'` still parse.
- PFT bound: `pft.wood_density` with >64 entries ⇒ `MAXPFT` error stop, no OOB write.
- Regression: real `meds_config_main.toml` + PFT file end-to-end ⇒ byte-identical PFT-parameters CSV.

**Risks.**
- Silent-default → error stop may surface EXISTING malformed-but-tolerated entries in shipped
  configs (a bool as `.true.`, a suffix-less duration); verify shipped files parse before shipping.
- `toml_logical` rejects spellings outside the accepted set; document the accepted set.
- Two coexisting report styles remain by design (ABSENT keys aggregate; present-but-malformed
  scalars abort immediately).
- A config that wrote a scalar where an array is expected now says "malformed" not "missing" — a
  wording change for any test asserting on the missing-key list.

---

### 4.5 Bug 6 — Soil interior heat advection upwind/sign — REFUTED (code is correct)

- **Severity:** none (refuted) · **Effort:** ~none (optional ~5 min doc clarifications)
  · **Status:** REFUTED
- **Files:** `src/biophysics/meds_column_energy.f90`, `src/driver/meds_column_dynamics.f90`,
  `src/biophysics/meds_biophysics_types.f90`

**Root cause — REFUTED.** The report assumes `soil_energy_flux` receives a DOWNWARD-positive
`w_flux`; it does not. The kernel's face-flux convention is UPWARD-positive, and under that
convention the upwind selection and sign at `meds_column_energy.f90:82-86` are both correct, mirroring
ED2 `rk4_derivs.f90:691-697`.

Evidence:
- **(a)** `hf` is upward-positive: `hf(0)=-forcing%g_top` (:77), `hf(k)=-kf*(t_new(k)-t_new(k+1))/dz_node`
  (:81), `hf(n)=forcing%geothermal` (:88). The divergence `div=(hf(k)-hf(k-1))+(qwf(k)-qwf(k-1))`
  (:94) enters ONE update, so `qwf` MUST share `hf`'s upward-positive convention.
- **(b)** The ONLY producer is `meds_column_dynamics.f90:341`:
  `eforc%w_flux(1:nsl) = -hflux%w_flux(1:nsl)` (down-positive hydro → up-positive energy). The
  hydrology source is documented DOWNWARD Darcy flux; negating yields UPWARD-positive.
  `test_column_energy.f90:140` sets `w_flux=-1.0e-7` ("gentle downward flow").
- **(c)** MEDS `k=1` top, `k=n` bottom (negative-z); face `k` between layer `k` (above) and `k+1`
  (below). ED2 uses the OPPOSITE ordering.
- **(d)** At :82-86, `w_flux(k)<=0` (downward) uses `t_new(k)` (layer above = source), `w_flux(k)>0`
  (upward) uses `t_new(k+1)` (layer below = source) — exactly ED2's scheme transposed to the flipped
  ordering. Physical check (downward infiltration of warm surface water) confirms the interior WARMS.
  `qwf(0)=qwf(n)=0` so the interior term telescopes to zero and does not perturb `energy_resid`.

**Fix.** No code change required. Optional doc-only hardening (behaviour-preserving):
1. In `meds_biophysics_types.f90:239` annotate `energy_forcing_t%w_flux` as
   `UPWARD-positive (caller flips the downward-positive hydrology flux, see
   meds_column_dynamics.f90:341)`.
2. Expand the terse comment at `meds_column_energy.f90:82` to spell out the upwind convention.

**Corrected math / units (documenting the correct current expression).**
```
qwf(k) = w_flux(k) * rho_h2o * internal_energy_liquid(T_source)      ! [W/m2], upward-positive
T_source = t_new(k)    if w_flux(k) <= 0   (downward, source = layer k above face k)
         = t_new(k+1)  if w_flux(k) >  0   (upward,   source = layer k+1 below face k)
```
Sign carried entirely by `w_flux` (no extra minus) since both `w_flux` and `qwf` are
upward-positive — identical to ED2 `rk4_derivs.f90:691-697`
(`qw_flux_g = w_flux_g*wdns8*tl2uint8(T_source)`). Because `qwf(0)=qwf(n)=0`, interior advection is
conservative and leaves `energy_resid` identically zero.

**ED2 reference.** `dynamics/rk4_derivs.f90:686-706` (source-temperature upwind, :691-697) and
:732-740 (conservative telescoping update). Producer of the MEDS up-positive flux:
`meds_column_dynamics.f90:341`.

**Test plan.** No fix to verify. To lock in behaviour and guard the latent risk, optionally add a
directional unit test to `test_column_energy.f90`: isothermal column, `w_flux(1:n-1)=-1.0e-6`
(downward) with a warm top layer, step, assert a lower interior layer's `soil_temp` rises while
`energy_resid` stays ~0. Mirror case: `w_flux>0` (upward) with a warm bottom warms an upper layer.
Existing suite remains green.

**Risks.**
- **Latent trap:** if a future caller populates `w_flux` from the DOWNWARD-positive hydrology flux
  without the sign flip, advection reverses. The undocumented direction on the type field makes this
  a trap; the optional doc fix removes it.
- The advection path is opt-in and OFF by default (`ccfg%advect_soil_heat=.false.`) and the existing
  test checks only budget closure, not direction — a directional regression here would not be caught.
  Consider adding a directional assertion if this path is promoted from opt-in.

---

### 4.6 Bug 7 — Dunne saturation-excess runoff spuriously large under free-drain/bedrock

- **Severity:** High · **Effort:** Low · **Status:** CONFIRMED
- **Files:** `src/biophysics/meds_column_hydrology.f90`, `test/test_column_hydrology.f90`

**Root cause (confirmed).** In `column_hydrology_flux`
(`meds_column_hydrology.f90:100-104`), `z_wt` is set — for every bottom BC except
`SOIL_BC_AQUIFER` — to `z_wt = z_bottom = params%soil_layer_z(n+1)`, a FIXED geometric constant, NOT
a water table tied to soil moisture. Line 122 then unconditionally computes
`f_sat = opts%f_max*exp(0.5*opts%f_over*z_wt)` and removes `q_over = f_sat*forcing%precip_ground`
from throughfall on EVERY step, regardless of actual saturation. Because `z_wt` is a fixed constant
under free-drain/bedrock, this yields a constant runoff fraction depending only on column depth:
with default `soil_depth=8.5 m`, `f_sat=0.4*exp(0.5*0.5*(-8.5))≈0.048` (~4.8% of ALL rain shed even
from a bone-dry column, forever); at `z_bottom=-2 m`, `f_sat≈0.24` (~24%); at `-1 m`, ~31%. This
contradicts the module design: `MEDS_COLUMN_HYDROLOGY_DESIGN.md` sec 3e ("`z_wt` is inert" under
FREE_DRAIN/BEDROCK), sec 1.2 (MVP runoff = infiltration-excess only; Dunne deferred to P2), sec 3d
(Dunne tagged P2/SIMTOP). The exponent sign is correct (`z_wt<=0`); the bug is substituting a finite
geometric bottom for the water-table depth. ED2 `lsm_hyd.f90` computes runoff from a genuinely
diagnosed water-table depth and gates on `soil_sat_water > 1.e-6`.

**Fix (recommended — option c: make Dunne a no-op unless a genuine water table exists).**

1. In `column_hydrology_flux`, replace the unconditional line 122 with a branch that only computes a
   nonzero `f_sat` under `SOIL_BC_AQUIFER` (where `col%z_wt` is a real prognostic water table updated
   at line 143):
   ```fortran
   if (opts%bottom_bc == SOIL_BC_AQUIFER) then
      f_sat = opts%f_max * exp(0.5_wp * opts%f_over * z_wt)
   else
      f_sat = 0.0_wp     ! no diagnosed water table under FREE_DRAIN/BEDROCK (design 3e/1.2)
   end if
   q_over = f_sat * forcing%precip_ground
   q_liq  = forcing%precip_ground - q_over
   ```
2. Add a one-line comment citing design sec 3e/1.2: only Horton infiltration-excess (ponding
   overflow, :173-175) generates surface runoff in the MVP.
3. Leave the `z_wt` assignment (100-104) as-is: `z_wt=z_bottom` is still a harmless datum for the
   default-off Zeng-Decker `psi_e` path (:108) where `z_wt` cancels; it is now no longer read by the
   Dunne term.
4. **OPTIONAL future upgrade (P2, do NOT do now):** diagnose a real water table by scanning from the
   bottom for the top of the contiguous saturated (`theta>=theta_sat`) zone; treat "no saturated
   layer" as effectively-infinite depth (`f_sat->0`). Deferred because under a unit-gradient
   free-drain BC the column essentially never saturates from below, so `f_sat=0` is the correct MVP.
5. Add a regression test asserting Dunne runoff is ZERO under free-drain and bedrock for an
   unsaturated column; the existing aquifer Dunne test (`test_column_hydrology.f90:395-413`,
   `col%z_wt=-0.1`) is unaffected.

**Corrected math / units.**
```
SIMTOP/TOPMODEL:  f_sat = f_max * exp(-0.5 * f_over * d_wt),  d_wt = water-table DEPTH >= 0 [m]
MEDS elevation :  z_wt <= 0,  d_wt = -z_wt  =>  f_sat = f_max*exp(0.5*f_over*z_wt)   (sign correct)
q_over = f_sat*precip_ground [kg/m2/s];  q_liq = precip_ground - q_over
Defaults: f_max=0.4, f_over=0.5 [1/m]
```
Spurious magnitude vs depth (dry column): `-8.5 m → 0.048`, `-2 m → 0.243`, `-1 m → 0.314`,
`-0.5 m → 0.353`. Setting `f_sat=0` ⇒ `q_over=0` ⇒ `q_liq=precip_ground`; runoff (:174) becomes
ponding-overflow only; `mass_resid` uses the same code path with `f_sat=0`, closure preserved. Under
aquifer BC nothing changes (`z_wt=col%z_wt` prognostic).

**ED2 reference.** `dynamics/lsm_hyd.f90`: `moist_zi` diagnosed from the site water balance (:161,
:381), saturation runoff gated on `soil_sat_water > 1.e-6` (:356), whole scheme gated by
`useTOPMODEL` (:122,:254). Bottom-BC enum `ISOILBC` in `memory/soil_coms.F90:87`. Design spec:
`MEDS_COLUMN_HYDROLOGY_DESIGN.md` sec 1.2 / 3d / 3e.

**Test plan.**
- `test_dunne_free_drain_zero`: dry free-drain column, sub-capacity precip (no ponding), assert
  `flux%runoff_surf == 0` exactly.
- Same with `SOIL_BC_BEDROCK`.
- Shallow-column variant: `z_bottom≈-1 m` under free-drain, sub-capacity precip, assert
  `runoff_surf == 0` (pre-fix shed ~31%).
- Aquifer regression: existing test (`col%z_wt=-0.1`, `SOIL_BC_AQUIFER`) still holds
  `runoff_surf ≈ 0.4*exp(0.5*0.5*(-0.1))*precip`.
- Heavy-downpour free-drain test (`:238-253`) still passes (Horton overflow independent of `q_over`).
- `flux%mass_resid` within `opts%atol` in all cases (Debug `debug_error=.true.`).
- Build + run nvfortran multicore in addition to ifx.

**Risks.**
- Behaviour change: downstream calibration/tests relying on baseline free-drain Dunne runoff see it
  drop to zero (intended correction).
- A user wanting saturation-excess runoff without an aquifer must switch `bottom_bc` to aquifer (or
  wait for P2); document that `f_max/f_over` are aquifer-only knobs in the MVP.
- Very low numerics risk: the fix only zeroes `q_over` on two BC branches; `mass_resid` closure
  structurally unchanged.

---

### 4.7 Bug 8 — atm↔CAS scalar conductance omits the M-O `temp1/temp2` (c3) profile factor

- **Severity:** High · **Effort:** Small-Medium (kernel signatures + one forcing-type edit + two
  test updates + doc corrections) · **Status:** CONFIRMED (kernels + docs; live driver already correct)
- **Files:** `src/biophysics/meds_column_energy.f90`, `src/biogeochemistry/meds_column_co2.f90`,
  `src/biophysics/meds_biophysics_types.f90`, `test/test_surface_energy.f90`,
  `test/test_column_co2.f90`, `archive/MEDS_COLUMN_CO2_BALANCE_DESIGN.md`,
  `src/biogeochemistry/README.md`, `src/driver/meds_column_dynamics.f90`

**Root cause (confirmed, with scoping correction).** The reusable stateless CAS-twin kernels form
the atm↔CAS scalar conductance as `rho*ustar` (enthalpy: `meds_column_energy.f90:267`
`gatm = rho_air*ustar`) and `can_dmol*ustar` (CO2: `meds_column_co2.f90:60`
`gatm_co2 = can_dmol*ustar`), both dropping the dimensionless M-O scalar-transfer coefficient `c3`.
This over-couples the CAS to the free atmosphere by `1/c3` (~5-7× near-neutral), suppressing diurnal
decoupling. HOWEVER: (1) the producer already exists — `mo_surface_layer` computes
`temp1 = vonk/d_heat(...)` (`meds_canopy_aerodynamics.f90:151,157`), which IS `c3`, stored in
`aero_out_t%temp1/temp2` (:61-63); (2) the **live** fast-loop driver `column_fast_step` ALREADY
consumes it correctly — `meds_column_dynamics.f90:314-316` sets `gah=rho*ustar*temp1`,
`gaw=rho*ustar*temp2`, `gac=can_dmol*ustar*temp2` and advances all three twins with those (:321-323).
So the defect is NOT in the running model; it is that the canonical named kernels
(`canopy_air_update`, `canopy_air_co2_update`) and `cas_atm_forcing_t` hardcode `c3=1`, leaving a
divergent wrong second code path plus stale docs.

**Fix (numbered).**

1. **Producer: NO CHANGE.** `temp1 = vonk/d_heat` (`meds_canopy_aerodynamics.f90:151,157`),
   `aero_out_t%temp1/temp2` written at :61-63 (`temp2=temp1` since `z0q=z0h`). This `temp1` IS `c3`.
2. **Forcing fields:** add `temp1` (heat) + `temp2` (vapour/scalar) to `cas_atm_forcing_t` (§3.2).
   Prefer a sentinel/non-default so an unfilled forcing fails a Debug guard instead of silently
   coupling at `c3=1`.
3. **Enthalpy twin:** in `canopy_air_update` (`meds_column_energy.f90:253-277`) add `temp1` (and
   `temp2`) args and change `gatm = rho_air*ustar` → `gatm = rho_air*ustar*temp1` (value used at :269
   and the :273 residual). Byte-matches the driver's `gah`.
4. **Humidity twin (RECOMMENDED, matches driver):** make the vapour twin IMPLICIT like the enthalpy
   twin and like `column_fast_step`. Take `shv_atm`, set `gatm_w = rho_air*ustar*temp2`:
   ```fortran
   shv_new = (wcapcan*cas_shv + dt*(coh_w_flux + coh_transp + ground_w_flux - dew &
             + gatm_w*shv_atm)) / (wcapcan + dt*gatm_w)
   ```
   (MINIMAL alternative: keep the explicit `w_flux_ac` but require callers to build it as
   `rho_air*ustar*temp2*(shv_atm - can_shv) = rho_air*ustar*qstar`.)
5. **CO2 twin:** in `canopy_air_co2_update` (`meds_column_co2.f90:39-70`) add a `temp2` arg and change
   `gatm_co2 = can_dmol*ustar` → `gatm_co2 = can_dmol*ustar*temp2` (:60); flows into `co2_new` (:62),
   `resid` (:64), `loss2atm` (:67). Thread `temp2` through `column_co2_step` (:184-214).
6. **Tests:** update `test_surface_energy.f90:144` and every `canopy_air_co2_update` call site in
   `test_column_co2.f90` (:95,108,120,125,142,155,160,176,234) for the new signatures + add
   regression assertions (see test plan).
7. **Recommended consolidation:** refactor `column_fast_step`
   (`meds_column_dynamics.f90:310-326`) to CALL the now-fixed kernels instead of re-inlining
   `enth1/shv1/co21`, so model and reusable kernels are ONE profile-factored path. At minimum keep
   them numerically identical to prevent re-divergence.
8. **Docs:** retire the `c3→1` caveat in `MEDS_COLUMN_CO2_BALANCE_DESIGN.md` (exec summary :71-75,
   sec 1.5 :165-174, sec 4a :765, sec 3.5.7 :609, decisions table :1339, limitations :1522-1523),
   `src/biogeochemistry/README.md`, and the matching note in the energy design doc. Reframe to
   "`c3 = aero temp1/temp2` threaded from the M-O solver, consistent across all three twins". Grep
   for `c3` across `archive/` and `src/biogeochemistry` before closing.
9. Rebuild + validate on BOTH ifx (full ctest) and nvfortran multicore (scalars only here; the
   array-temp trap does not apply, but green ifx is not sufficient).

**Corrected math / units.** See the Appendix for the full `temp1` explainer. In brief:
`temp1 = c3 = vonk / [ ln((zref-d)/z0h) - psih(zeta) + psih(zeta0h) ]`, dimensionless, `~0.05-0.4`
(~0.15-0.20 near neutral). Correct atm↔CAS scalar conductance is `g = rho*ustar*c3`; `rho*ustar`
alone (`c3=1`) is wrong. `g_heat = rho*ustar*temp1`, `g_vap = rho*ustar*temp2` `[kg/m²/s]`; CO2 uses
dry-air molar density: `gatm_co2 = can_dmol*ustar*temp2 [mol_dryair/m²/s]`. Omitting `c3` inflates
the conductance by `1/c3` (~5-7× near neutral, up to ~20× stable), damping diurnal decoupling and
under-predicting nocturnal sub-canopy CO2 accumulation.

**ED2 reference.** `ed_stars8` in `dynamics/canopy_struct_dynamics.f90:1798-1800`
(`phistar = c3*(phi_atm - phi_can); F_phi = rho*ustar*phistar`), `c3 = vonk/(Prandtl*(ln((z-d)/z0h) -
psih(zeta) + psih(zeta0h)))`; CAS CO2 ODE `dynamics/rk4_derivs.f90:2163`
(`cflxac = can_dmol*ustar*cstar`). MEDS's `temp1 = vonk/d_heat` is exactly this `c3` with CLM4/CLM5
denominators (Prandtl folded to 1).

**Test plan.**
- `test_surface_energy.f90` (`test_cas`): pass `temp1=temp2<1` (e.g. 0.18) into `canopy_air_update`.
  Assert (a) CAS energy residual ~0 over 30 steps; (b) CAS approaches the atmosphere STRICTLY SLOWER
  than with `temp1=1` (run `temp1=1` vs `temp1=0.18`, assert larger final gap for smaller `temp1`).
- `test_column_co2.f90`: pass `temp2<1` at each call site. Assert `budget%resid ~ 0`; halving `temp2`
  slows equilibration of `can_co2` toward `co2_atm`; expected time constant
  `tau = can_depth/(ustar*temp2)`.
- New parity/anti-regression test: surfaces zeroed (pure atmospheric relaxation), assert the
  standalone `canopy_air_update`/`canopy_air_co2_update`, fed the same `ustar` + aero `temp1/temp2`,
  reproduce `column_fast_step`'s `enth1/shv1/co21` (`:321-323`) to machine precision — LOCKS the
  kernels to the live driver.
- If the humidity twin becomes implicit, add a test that the vapour residual closes and matches the
  driver's `shv1`.
- Full ctest on ifx + nvfortran multicore build.

**Risks.**
- Divergent code paths: the live model already applies `temp1/temp2`, so runtime NEE/energy does not
  change — only the reusable kernels + docs do. The recommended consolidation removes the two-path
  risk; if deferred, add the parity test.
- Making the humidity twin implicit changes its discretization (explicit → L-stable implicit) — more
  correct and matches the driver, but re-verify the residual and update tolerances.
- Silent-regression hazard: if `temp1/temp2` default to 1.0 and a caller forgets to fill them, the
  bug returns silently — mitigate with a sentinel + Debug guard, or make the aero→forcing copy the
  only construction site.
- Signature churn touches `test_surface_energy.f90` + ~9 `test_column_co2.f90` sites (compile-time,
  low silent-miss risk).
- Doc consistency: the `c3→1` note recurs in ~6 places + README + energy doc; grep before closing.

---

## 5. Combined Test & Validation Strategy

**Per-fix unit tests** (owned by each section above):

| Fix | New/updated unit tests |
|-----|------------------------|
| bug1+3 | `test_fast_loop.f90`: carbon+fast end-to-end (case 4), guard-trap, reset-semantics, lockstep fuse/split integrity |
| bug2 | `meds_demography_fusefiss` carbon-mode fuse + fission conservation, `hgt_max`-immunity, guard-trip, empirical bit-identical |
| bug4 | `test_leaf_physiology.f90`: 25 °C regression, warming monotonicity, `kp/vcmax` invariant, peaked deactivation |
| bug5 | negative-config tests (malformed token, short/over-reported array, logical spellings, duration, PFT bound) + valid-config regression |
| bug6 | (optional) directional advection test in `test_column_energy.f90` |
| bug7 | `test_dunne_free_drain_zero` (free-drain + bedrock + shallow), aquifer regression, mass-resid closure |
| bug8 | CAS-twin slowdown assertions + driver/kernel parity test to machine precision |

**Cross-cutting regression suites:**

1. **Empirical bit-identical guarantee (bug1+3, bug2).** Run the default spin-up (`GS_EMPIRICAL`,
   fast off) before and after the demography changes; assert byte-identical `-D-output.nc` and
   `print_summary`. This is the primary guard that the state-SoA and fusefiss edits do not perturb
   the shipped default path.

2. **Carbon-conservation regression across fuse/split (bug1+3 × bug2).** On a carbon-mode
   (`GS_CARBON`) multi-cohort site, force both fusion and fission, and assert conservation of every
   extensive per-ground quantity: `Σ(nplant·wood_carbon)`, `Σ(nplant·leaf_carbon)`,
   `Σ(nplant·fineroot_carbon)`, `Σ(nplant·nonstructural_carbon)`, `Σ(nplant·agb)`, and the three new
   `Σ(nplant·*_resp_accum)` accumulators. This is the single test that proves the coupled bug1+3/bug2
   field plumbing is complete — a forgotten `nplant`-weighted merge site trips it (or trips bug2's
   runtime conservation guard).

3. **Diurnal-cycle offline test (bug8, and downstream bug1+3), once forcing exists.** No production
   forcing loader exists yet (fast loop OFF by default, MVP constant forcing in
   `build_fast_context`). When an offline forcing driver lands, run a single-column day/night cycle
   and assert: (a) the CAS decouples from the free atmosphere during stable nights (small `c3`),
   showing nocturnal sub-canopy CO2 accumulation; (b) daytime canopy-air warm/moist/CO2-drawdown
   departure is present (not damped as under `c3=1`); (c) the three CAS twins share one `ustar` and
   one `c3`. Until then, the driver/kernel parity test (bug8) is the stand-in that locks kernel
   physics to the already-correct live driver.

4. **nvfortran multicore build check (ALL fixes).** Per `CLAUDE.md`, a green ifx run is not
   sufficient. Every fix must build and pass ctest under both `ifx -stand f18 -check all` (Debug,
   `-DMEDS_ENABLE_IO=OFF`) and `nvfortran -DMEDS_GPU=multicore`. Special attention for bug1+3: the new
   optional array out-args (`*_resp_coh`) must be bound to **named** arrays, never
   function-result temporaries, to avoid the nvfortran array-temporary miscompile (`CLAUDE.md`
   issue #7). bug4/5/7/8 are scalar-only, so the array-temp trap does not apply there.

---

## Appendix A — What `temp1` (the M-O scalar-transfer coefficient c3) Is

*Self-contained note, cleaned from bug8's `math_or_units`. Paste-ready for the design docs.*

**What `temp1` is.** `temp1` is the dimensionless Monin-Obukhov **bulk scalar-transfer coefficient**,
ED2's `c3`. It is the proportionality between a turbulent scalar flux and the atmosphere-minus-canopy
difference of that scalar: for any scalar φ (potential temperature, specific humidity, CO2), the
eddy flux into the canopy air space is

```
F_phi   = rho * ustar * phistar
phistar = c3 * (phi_atm - phi_can)
```

So the correct atm↔CAS scalar **conductance** is `g = rho * ustar * c3`, i.e. `temp1` supplies the
missing `c3`; `rho*ustar` alone (`c3=1`) is wrong.

**Formula (as MEDS computes it, `meds_canopy_aerodynamics.f90` `d_heat` + `mo_surface_layer`).**

```
temp1 = c3 = vonk / D_h,   with (neutral/unstable branch)
D_h        = ln( (zref - d) / z0h ) - psih(zeta) + psih(zeta0h)
=>  c3     = vonk / [ ln((zref-d)/z0h) - psih(zeta) + psih(zeta0h) ]
```

Inputs: `vonk` = von Karman (0.4); `z0h` = scalar roughness (= `z0m` = `z0m_ratio*veg_height` here);
`zldis = zref - d` (reference height above displacement `d = d_ratio*veg_height`);
`zeta = zldis/L` the M-O stability parameter (from the fixed-iteration `ustar/L` solve);
`zeta0h = z0h*zeta/zldis = z0h/L`; `psih` = integrated heat/scalar stability correction (CLM
`StabilityFunc2`). The neutral Prandtl number is folded into the CLM four-range denominators
(Prandtl = 1), so `temp1 = vonk/d_heat` matches ED2's `c3 = vonk/(Prandtl*(...))`. `temp2` is the same
coefficient for vapour/CO2 (`z0q`); it equals `temp1` today because `z0q = z0h`. `temp1` is finite and
strictly `> 0` (`D_h > 0` in every branch), so restoring it introduces no division-by-zero and is
`-fpe0` safe.

**Units.** `c3` (`temp1/temp2`) dimensionless `[-]`. `ustar [m/s]`. `rho_air [kg/m3]`. So
`g_heat = rho*ustar*temp1` and `g_vap = rho*ustar*temp2` are mass conductances `[kg/m2/s]`; the CO2
twin uses the DRY-AIR MOLAR density, `gatm_co2 = can_dmol*ustar*temp2 [mol_dryair/m2/s]`, with
`can_dmol = rho_air*(1-can_shv)/mmdry [mol/m3]`.

**Magnitude.** `c3 ~ 0.05-0.4`, roughly `0.15-0.20` near neutral (e.g. `veg_height 20 m`,
`zref 30 m`: `d~12.6 m`, `z0h~2.6 m`, `zldis~17.4 m`, `ln(zldis/z0h)~1.9`, `c3~0.4/1.9~0.21`).
Unstable (daytime, `psih<0`) shrinks `D_h` ⇒ `c3` rises toward `~0.3-0.4`; stable (nocturnal,
`psih>0`) grows `D_h` ⇒ `c3` falls to `~0.05-0.1`.

**Physical consequence of omitting it (`c3=1`).** The atm↔CAS conductance is inflated by `1/c3`
(~5-7× near neutral, up to ~20× under stable stratification), so the CAS enthalpy/humidity/CO2 are
relaxed toward the free-atmosphere values far too fast. The canopy air space then hugs the
atmosphere, and the model loses diurnal decoupling: the daytime canopy-air warm/moist/CO2-drawdown
departure is damped and, critically, nocturnal sub-canopy CO2 accumulation (which physically depends
on WEAK stable-night coupling, small `c3`) is strongly under-predicted. Restoring `temp1/temp2` makes
all three CAS twins share one turbulence basis (`ustar`) AND one consistent profile factor (`c3`),
exactly as ED2's `ed_stars8` does.
