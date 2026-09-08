# MEDS code review, 2026-09-08 — items 1–3 (budgets, state ownership, main inventory)

**Status:** review record, read-only. No source changed by this document. Reviewer: Claude (Fable 5.1)
with three parallel audit passes, every top-ranked finding re-verified by hand against the tree at
`2b1c015` (= tag `v0.1.0`). Line numbers refer to that commit.

Items 4–10 of the review plan (interfaces, helper consolidation, naming, numerics, flags, efficiency,
ED2 fidelity) are NOT covered here.

---

## Item 3 — what is on `main`

- `main` == `v0.1.0` (`2b1c015`). Working tree clean. Every feature branch in the project history
  (leaf-wood-energy #41, snow #42, ed2-rk45-water-mass #67, integrator-parity #77/#88, transp
  corrector #91, io-v01 #111) is MERGED. Only `origin/docs/readme-figures` remains and it is merged.
- Two stashes (2026-07-30, branch `integrator-physics-parity`, "option B") are DEAD: stash@{0}
  touches `src/driver/meds_fast_split.f90` and `test/test_picard_coupling.f90`, neither of which exists;
  stash@{1} adds `leaf_temp/wood_temp` to `column_state_t`, which PR #88 delivered differently. Safe
  to `git stash drop` both.
- Design-doc status headers are stale in ~8 files (say "design-only"/"branch X" for shipped work:
  SLOW_DYNAMICS, BIOGEOCHEMISTRY, SNOW, PHENOLOGY_RATE_REFACTOR, IO_V01, LEAF_WOOD_ENERGY ...).
  `docs/dev_plans/README.md` declares them point-in-time, so this is hygiene, not error.
- Physics gated OFF by default: `canopy_water_on=.false.` (no interception film; condensate goes to
  soil layer 1, issue #74), `soil_carbon_on=.false.`, `trait_plasticity_on=.false.`. The first two are
  "flag whose off path is known-incomplete physics" candidates for item 8.
- 28 `build-*` directories at repo root (gitignored). Open issues: #118 #117 #104 #96 #89 #74 #47 #7 #6 #1.

---

## Item 1 — budget/conservation checks

### 1A. Fast-loop ledgers (`meds_budget_check` + callers)

Ranked, all verified:

1. **Store-scaled tolerance + no signed accumulation ⇒ blind to sustained ~1 W/m² leaks.**
   `meds_fast_ark.f90:1385-1402`, `meds_fast_rk45.f90:825-837`: `scale = abs(e_soil1 + wcap*enth1)`
   ≈ 1–1.5e9 J/m² (liquid datum `tsupercool_liq` ≈ 56 K ⇒ u_liq ≈ 1 MJ/kg), `rtol = 1e-6` ⇒ tol
   ≈ 1e3 J/m²/step ≈ 1.1–1.7 W/m² at 900 s. `budget_t` keeps `resid` (last), `worst` (max |resid|),
   counters — no signed sum (`meds_budget_check.f90:24-33,63-71`). A one-signed 0.5 W/m² bias passes
   every step and drifts soil ~3 K/yr with `n_fail = 0`.
   Water: `max(w_soil1+…,1)`·1e-6 ≈ 6–9e-4 kg/m²/step ≈ 0.06–0.09 mm/day — marginal.
2. **`atol = 5e6 J/m²` whenever `canopy_water_on`** (`ark:1400-1402`, `rk45:835-837`) = 5.6 kW/m²
   per step: the whole-energy check is effectively disabled on canopy-water runs. Comment admits it
   hides the film-valuation mismatch (film valued at `u_liq(rain_temp)`, evaporation credited to CAS
   at `enthalpy_vapor(tl)`, leaf pays only `latent_heat_vap`).
3. **`n_fail`/`worst` never read at run end** (`meds_stepper.f90:50,52` pass no `worst_*`,
   `n_budget_fail`). The one production channel `PD_RESID_ENERGY` (`meds_fast_dynamics.f90:1144-1145`)
   accumulates `worst*dt` and is time-averaged and labelled W/m² (`meds_core_diag_types.f90:167`): it
   is an unsigned running max in J/m², off by a factor dt from its label.
4. **Hard stop dead in every configured run**: `halt_budgets = debug_error .and. mask_is_full`;
   `debug_error` defaults false and no integration test sets it.
5. **ARK books the theta_res floor's energy (`e_floor`, `ark:299`) but not its mass** in
   `whole_water` (`ark:1379`; `acc%whole_wat_in ≡ 0` at :302). RK45 books both (`rk45:790,797`).
6. **Condensate deposit valued at end-of-step T_cas, debited per stage** (`ark:1280`, `rk45:604-606`
   vs `meds_fast_time_derivs.f90:251`); RK45's correct `bw_cond_enth` (`rk45:272`) is never used.
7. **Unbooked mass edits on the committed path**: `max(·,tiny_num)` floor on tissue water
   (`ark:860-863`, fires under issue #104); slow→fast seam `clamp_water_to_capacity` discards excess
   (`meds_fast_dynamics.f90:583-589`); film `surf_overflow` booked as OUT but not routed to ground.
8. **RK45 asserts ledgers on a step the dispatcher may roll back** (`rk45:815-837` vs
   `meds_fast_step.f90:112-117`).
9. **Per-layer / per-cohort checks exist only in the scratch hydrology under debug_error**
   (`meds_soil_water.f90:299,305`); RK45's own theta trajectory has none; no tissue water/energy
   residual anywhere. Vertical redistribution errors are invisible to every whole-column ledger.
10. **PHYSICS, closes the ledger but is wrong (verified): the soil root heat sink is charged twice
    for the liquid enthalpy of transpired water.** `coh_qsoil = Σ transp·(h_vap(tl) − L) = transp·u_liq(tl)`
    (`time_derivs:182`) AND `qloss = uptake·u_liq(soil)` (`ark:2035`) both enter
    `root_heat_sink` (`ark:233`, `time_derivs:361`); the leaf receives `qwflux_wl = sapflow·u_liq(wood)`
    (`time_derivs:175`) but pays only `L` for evaporation (`:163-164`) while the CAS receives the full
    `transp·h_vap(tl)` (`:182`). Steady state: soil −2x, leaf +x, CAS +x, wood 0 ⇒ ledger closes, but a
    spurious soil→leaf transfer of `transp·u_liq ≈ 30 W/m²` at 3 mm/day. The pre-P2 proxy (soil pays
    for the vapour's liquid part) was never removed when P2 added the explicit soil→wood→leaf chain.
    Tell-tale: the result depends on the energy datum, which a correct formulation cannot.
    Fix: leaf pays `transp·u_liq(tl)` in addition to `L` (i.e. loses `h_vap(tl)` to the CAS) and
    `coh_qsoil` leaves the soil sink. Check against ED2 `rk4_derivs` leaf/CAS enthalpy terms.

Proposed utility changes (ranked): (i) `budget_t` gains `resid_sum` (signed) and `flux_gross`;
accumulator NOT reset per slow step (`fast_dynamics:607`) — keep a per-site run accumulator;
(ii) `scale = influx + outflux`, `atol = rate_floor·dt` (0.01 W/m², 1e-8 kg/m²/s); delete the 5e6
merge and fix the film valuation; (iii) `budget_report` at end of run from `meds_stepper`;
(iv) `PD_RESID_*` accumulate signed `resid` (J/m²) and `|resid|/dt` (W/m²); (v) add
`sum(floor_layer)*dt` to ARK `w_in`, deposit condensate at b-weighted enthalpy; (vi) per-layer
`budget_imbalance` over faces on the committed path, per-cohort tissue water/energy; (vii) run RK45
ledgers after the rail decision.

### 1B. Slow loop, fusion/fission, seams (cohort → patch → site)

Ranked, top four verified by hand:

1. **Cohort fission duplicates canopy film water.** `split_cohorts` halves `nplant`
   (`meds_core_cohort_fusefiss.f90:281`) then `copy_cohort_slot` copies `leaf/wood_surf_water`
   verbatim (`meds_core_state_types.f90:757-758`); these are [kg/m² ground] (`:172-173`) and fusion
   SUMS them (`fusefiss:214-215`). Every split doubles the film. Only AGB asserted (`:302`). Same in
   `apply_patch_disturbance` survivor copy (`meds_core_patch_fusefiss.f90:457-458`). Diag/sdiag
   accumulators are also copied, not halved.
2. **Patch fusion does not rescale ground-referenced film water** (`patch_fusefiss:243,248` scale
   only `nplant`). `test_patch.f90` checks N, area, theta, CAS enthalpy, soil energy only.
3. **Dormant deciduous cohorts get wood water reset to PSI_INIT daily.** Lazy-init sentinel
   `leaf_water_mass <= 0` (`meds_fast_dynamics.f90:567`) re-seeds BOTH leaf and wood. Snap-to-bare
   shed (`vegetation_dynamics.f90:456-458`) drives `leaf_water_mass = 0`; every dormant day `bleaf=0`
   ⇒ leaf seed 0 ⇒ sentinel trips again ⇒ yesterday's integrated wood water discarded. "0 =
   uninitialised" is ambiguous with the reachable state "bare canopy".
4. **Recruitment carbon ≠ reproduction debit.** `npp_to_recruitment` divides by `carbon_min` =
   AGB only (`meds_plant_vital_rates.f90:56`), but `init_cohort` gives each recruit wood + leaf +
   fineroot + storage (`state_types:815-855`). Reproduction is never a carbon store (only a density
   survives in `patch%recruit_pool`); with `repro_carbon_efficiency = 1e-3` 99.9 % of repro carbon
   vanishes while the recruits that appear carry more than the debit. `seed_rain_recruits` is an
   undeclared external source.
5. **Starvation `deficit` and `growth_resp` never reach the CAS carbon balance** (verified:
   `meds_plant_carbon_allocation.f90:98-101`, `vegetation_dynamics.f90:390-397`; zero hits for
   `growth_resp` in fast/biophysics). CAS-side NEE ≠ pool-side ΔC by (growth_resp − deficit) daily.
6. **No site-level carbon ledger at any timescale.** Grid-scale checks are area (`meds_main.f90:314`)
   and a NaN scan. `rh_seam_gap` is computed (`meds_biogeochem_dynamics.f90:149-151`) but dropped by
   the production stepper. `soilc_audit_t%resid` is a tautology (`rh_today := litter_in − dc_pool`).
7. **Mortality, cull, disturbance-kill discard tissue water, film water, tissue heat**; cull and
   disturbance carbon reach soil only when `soil_carbon_on`, else lost; their litter bypasses
   `audit%litter_in`.
8. **CENTURY matrix conservation asserted nowhere** (`assemble_transfer_matrix`,
   `meds_soil_biogeochem.f90:162-223`; no check `0 ≤ er_j ≤ 1`).
9. **Second-order**: `terminate_patches` rescales survivors by `1/s` (all stores); disturbance gap
   gets `recruit_pool = 0`; `blend_cas` blends specific quantities AND `can_depth` (non-conserving
   when depths differ); mortality litter uses unfloored fraction when the `negligible` floor binds;
   `fuse_2_cohorts` leaf/wood temperature is leaf-area-weighted (conserved quantity is `hcap·T`).
10. **`conservation_tol` (0.001 in every shipped TOML, required key) is applied to merges whose
    error is round-off (~1e-16)** — a 0.1 %/fusion leak passes. `size_tol` (`meds_constants.f90:78`)
    is unused.

Proposed: (i) `site_ledger_t` snapshot + daily/annual assert in `advance_slow_dynamics` for C,
water, energy with seam events (seeds, clamps, mortality water) as DECLARED boundary terms;
(ii) fix extensive/ground-referenced field policy in `copy_cohort_slot`/split/patch-fuse/disturbance
and add `Σ surf_water` to the asserts; (iii) replace the `≤ 0` sentinel with an explicit seeded flag
(seed leaf and wood independently); (iv) fuse/split asserts at ~1e-12 relative, keep
`conservation_tol` for the site ledger, delete `size_tol`; (v) assert matrix columns and `rh_seam_gap`
in production; make `audit%litter_in` include cull/disturbance adds.

---

## Item 2 — state ownership / frozen-flux defect class

Ranked; top four verified by hand.

1. **RK45 never fills `psi_leaf_coh`.** `meds_fast_step.f90:112-114` returns on RK45 success before
   the fill at `:126-132`; `psi_leaf_pool` (`meds_fast_dynamics.f90:452-454`) is uninitialised and
   consumed by `max()` at `:772`. In `time_integrator="rk45"` runs the predawn-ψ accumulator holds the
   previous patch's (or garbage) value ⇒ `beta_stomata → 1`, stress removed, thread-count dependent.
   RK45 is not a valid oracle for any water-stressed comparison.
2. **Thread race: `npp_patch` is shared.** Declared at `fast_dynamics:309`, missing from the
   `private` clause at `:465-466`, written `:723`, read `:726` ⇒ FAST-tier `npp_rate` is
   thread-count dependent. (Post-dates the PR #109 byte-identity test.)
3. **RK45 applies the ψ-wilting factor twice at the wood↔soil interface.** `fro%uptake =
   hflux%uptake_total` already contains `fwilt` (`meds_soil_water.f90:498,581`); RK45 feeds it back
   through `soil_water_time_deriv`, which applies `fwilt` again (`:541-543,581`). Wood is credited the
   full `uptake_frozen` ⇒ water fabricated ∝ `uptake·(1−fwilt)·dt` on dry soil. ARK does not call
   `soil_water_time_deriv` (commits scratch θ1) so it is unaffected.
4. **RK45 lacks the PR #91 transpiration corrector** (`column_derivs:407`: `sapflow_frozen − transp_i`)
   ⇒ the "reference" now carries the exact `dt·(transp_pp − transp_realised)` inconsistency the
   corrector removed from ARK. ARK-vs-RK45 ψ_leaf comparisons at production dt measure this term.
5. **Reported LE/H/ET use state-n conductance on end-of-step state** (`atm_fluxes`,
   `meds_fast_step.f90:148-151`, `aero` from the pre-pass) while the CAS exchanged at per-stage
   re-solved `gah/gaw`. Headline ET can differ from the ledger's b-weighted `acc%whole_wat_out`,
   which is not exported.
6. **Tissue-water `tiny_num` floor fabricates water with no ledger entry and is known to fire**
   (`ark:860-862`; #104 pinned-wood case).
7. **Corrector asymmetry**: committed mass moves at `sapflow_c` (`ark:846`) while advected enthalpy
   (`qwflux_wl`, `ark:2034`) and `CD_SAPFLOW` (`ark:1947`) use `sapflow_b`.
8. **Per-cohort diagnostics are state-n pre-pass values**, not committed fluxes (`CD_ROOT_UPTAKE`
   pre-scale, `CD_TRANSP` from the state-n leaf kernel, `CD_*_TEMP` previous-step averages).
9. **Leaf LW linearisation base moves with the stage while `abs_lw` is frozen at `tcas^n`**
   (`time_derivs:162,173-174`); ground has no LW slope at all. Phantom radiative input booked into
   `coh_rnet`. Undocumented.
10. **`ggnet` still frozen** (`types:343`, `time_derivs:223`, set once `ark:1847`); no comment states
    the freeze. Known: 0.70 K soil-T error at 900 s.
11. ARK stage-3 `clamp_cas` shift is committed if the step is accepted (`ark:520-523`), contrary to
    the comment at `:497-500`. Low probability.
12. Warm-start seed shared across schemes after an RK45→ARK rescue (`fast_step:116`, `ark:1118`,
    `rk45:509`).
13. Stale text: RK45 "C5 rejects SOIL_BC_AQUIFER" (no such guard); `column_frozen_t` wood/leaf
    fields (`types:505-517`) documented for routines that no longer exist; `ark:1019-1027` header
    contradicts `:24-30`.

Clean: `surface_derivs`/`column_derivs` are `pure`; every kernel they call is pure/elemental; no
module-scope mutable state in fast-loop or biophysics modules; the only SAVE is `write_fast_probe`
(single-thread guarded). Output layer reads no `site%deriv`. Per-patch warm start is clean except #1.
The ARK frozen-field discipline is largely deliberate and documented (see the audit's frozen-field
table in the session transcript; key rows reproduced by classification: (a) documented: h_coeff,
g_tr, a_leaf/a_wood, film conductances, snow stage, scratch hydrology authority, uptake_frozen;
(b) undocumented: abs_lw base, ggnet, qwflux at sapflow_b; (c) defect: items 1–4 above).

---

## Fix log

- **2026-09-08, commit 1 (item 2 #1-#3):** RK45 psi_leaf report, `npp_patch` private, double
  wilting factor -- shipped with regression tests; ARK serial byte-identical to v0.1.0.
- **2026-09-08, commit 2 (item 1A utility):** `budget_t` gains signed `resid_sum`, `abs_sum`,
  `flux_gross`, `elapsed`; `budget_check` (flux-scaled tolerance + rate floor), `budget_merge`,
  `budget_report`; run-level ledgers reported by `meds_main`; `resid_*_site` are now signed mean
  rates with honest units; ARK books the theta_res floor mass and the clip mass; condensate deposit
  carries the b-weighted stage enthalpy on both schemes (item 1A #5, #6). The `canopy_water_on`
  slack survives as a NAMED `atol_extra` (item 1A #2) until the film valuation is fixed.
  **Open, found by the tighter tolerance:** under a melting/refreezing snow pack the ARK whole-
  column energy ledger carries an unattributed residual of up to ~0.7 J/m2 per 150 s step
  (5e-3 W/m2), oscillating in sign (-14 J/m2 net over a day). The CAS, soil, pack and pond ledgers
  each close to 1e-9 on the same steps and the condensate enthalpy is not it; it correlates loosely
  with the tissue heat-store change. The energy rate floor is set to 1e-2 W/m2 to leave it visible
  but non-fatal.
- **2026-09-08, commit 3 (item 1A #10, #2 -- shared physics):** the leaf pays the FULL vapour
  enthalpy of what it transpires (`enthalpy_vapor(t_cas)`, the linearization reference, on both the
  leaf side and the CAS credit); the soil root heat sink is `qloss` only (advected enthalpy of the
  extracted water, once), matching ED2 `rk4_derivs` (`qtransp = transp*tq2enthalpy(T_leaf)`,
  `qloss = wloss*uint_water`). The canopy film pays `enthalpy_vapor - u_liq(rain_temp)`, so film +
  tissue + CAS close exactly and the 5e6 J/m2 `atol_extra` slack is deleted. New test: the canopy is
  energy-neutral (absorbed == sensible + full vapour enthalpy handed to the CAS); it fails on the old
  source by exactly the old soil proxy (~108 W/m2 in the fixture). **Result change (July example,
  new - main):** soil-top temperature +4.0 K mean / +6.3 K max, leaf +0.25 K, CAS +0.20 K, sensible
  heat -8 W/m2 mean, GPP +0.24 umol/m2/s mean. `docs/ed2_comparison.md` soil-temperature statements
  predate this and need re-running.
- **2026-09-08, commit 4 (item 1B #1-#3 -- fusion and seam):** new `scale_cohort_ground_fields`
  applies the nplant factor to the ground-referenced interception films on cohort fission (0.5),
  patch fusion (area weights) and the disturbance survivor copy; `split_cohorts` and
  `fuse_2_patches` now assert film conservation at 1e-12. The fast-loop lazy seed tests leaf and
  wood water INDEPENDENTLY, so a bare (dormant) cohort no longer has its wood water re-seeded at
  PSI_INIT every day; the leaf-out seed remains an undeclared source until a slow ledger exists.
  Tests: split, two-patch fusion, two-donor disturbance (a one-donor case cannot discriminate --
  its dilution factor is 1), and a masked-hydraulics one-step dormant run; all four fail on the
  previous source.
- **2026-09-08, commit 5 (golden regeneration):** `examples/example_biophysics` re-spun from bare
  ground for 50 years on the fixed code (10 min, 4 threads), July re-run, the three README figures
  regenerated and every figure-derived number in `README.md` and the example README refreshed
  (soil surface 23.4 C mean, was 18.3; 1.7 m depth 16.1 C, was 9.9; net July uptake 26.9 gC/m2, was
  92.5, because the warmer soil respires more; leaf-air excess +3.8 K, was +4.3). The old-physics
  spin-up states are kept beside the new ones as `spinup-old-physics-S-*.nc` (gitignored).
  **Run-level ledger over the 50-year spin-up:** water closes to -7e-11 kg/m2; energy accumulates
  +6.75 MJ/m2 (mean +4.3e-3 W/m2, worst step 544 J/m2, 9% of checks breach the 1e-2 W/m2 floor at
  dt_fast = 900 s). A one-year daily diagnostic (`resid_energy_site`) localises it: **Dec-Mar only**,
  +1.6 to +3.2e-3 W/m2 while the soil surface sits at 273 K; Apr-Nov close to <1e-6 W/m2. It is the
  same class as the sub-1 J/m2 melting-pack residual seen at 150 s in `test_column_dynamics` (the
  CAS, soil, pack and pond sub-ledgers each close). Candidates: the soil freeze/thaw plateau
  bookkeeping across the ARK stage/commit seam, or the frozen-surface film/condensate valuation.
  **RESOLVED, commit 6 (below).**
- **2026-09-08, commit 6 (the winter residual -- TWO mechanisms, both found by a per-step probe):**
  1. *Pond enthalpy valuation.* Every large step (up to 11 J/m2 at 150 s) had NO pack, sub-freezing
     precipitation routed to the ground as liquid at the canopy-air temperature (~270 K: rain, or
     sub-threshold snowfall), an EMPTY pond and a soil top at 273 K. In `column_hydrology_flux` the
     pond received the water at 270 K, `uext_to_temp` put it on the melt plateau (T = t_3ple, fliq<1),
     infiltration was valued at `u_liq(t_3ple)` -- more than the water carried by L_f*(1-fliq) -- the
     pond drained negative and the empty-pond reset zeroed the deficit. Predicted
     m*cp_liq*(t_3ple - t_precip) = 11.0 J/m2 on the probe step, observed 11.3. Fix: infiltration and
     overflow are valued at the pond's MEAN specific enthalpy e/w via `temp_of_liquid_enthalpy(e/w)`
     (exact inverse of `internal_energy_liquid`), in the shared kernel and RK45's own pond overflow.
     `test_pond_subfreezing_inflow` fails on the old source by exactly this amount.
  2. *Shed water under a pack.* The remaining bias was a per-patch CONSTANT residual on every step
     from November on, independent of snow, melt or precipitation, and equal to
     `shed_water_rate*dt*u_liq(t_melt)` (1.97e-6 kg/m2 x 9.06e5 J/kg = 1.78 J/m2, observed 1.78). The
     daily leaf/root-turnover shed water (P4) is routed into the pond together with meltwater at
     `t_precip = t_melt` whenever a pack exists, but both whole-column ledgers ZEROED the ground-
     inflow enthalpy term under a pack (`merge(0, precip_ground*dt*u_liq(rain_temp), snowfac>0)`,
     with `rain_temp` pinned at `tsupercool_liq`), so the shed water's enthalpy entered the pond
     unbooked. A restart clears `shed_water_rate` (it is not in the state file), which is why a
     January restarted from its own state was clean while the continuous run was not. Fix: the
     boundary term is `(precip_ground - snow_melt_rate)*dt*u_liq(t_precip)` on both schemes -- the
     non-melt inflow at the temperature the pond receives it; identical to the old term without a
     pack. New `surface_frozen_t%snow_melt_rate`.
  Threading was checked and exonerated (bit-identical ledgers at 1/2/4 threads). Results: January
  2074 cumulative +35.5 -> -0.34 J/m2 (mechanism 1); July 2074 -> January 2075 continuous run
  +4127 -> -2.1 J/m2, worst step 2.75 -> 0.057 J/m2, 0/1165824 breaches (both). One-year daily
  diagnostic (Jul 2074 - Jul 2075, 150 s, 4 threads): Dec-Mar mean residual 1.5-3.1e-3 W/m2 ->
  |<= 5e-7| W/m2, year cumulative -4.5 J/m2, 0/2312640 breaches; water -4.5e-12 kg/m2.
  Follow-ups noted, not done: `shed_water_rate` is not persisted in the restart file (a restart
  loses at most one day of it); sub-threshold snowfall onto bare ground is deposited as LIQUID at the
  canopy-air temperature (the fusion enthalpy is created; ledger-consistent but physically wrong).

## Suggested order of fixes

1. Item-2 #1, #2, #3 (RK45 ψ fill, `npp_patch` private, double `fwilt`) — small, unambiguous,
   and they restore RK45 as an oracle and the threading byte-identity guarantee.
2. Item-1A #10 (double-charged root heat sink) — physics; verify against ED2 first, then a
   datum-invariance test.
3. Item-1B #1–#3 (film-water duplication on split/fuse/disturbance; dormancy sentinel).
4. Ledger utility: signed cumulative residual, flux-scaled tolerance, end-of-run report,
   `PD_RESID_*` units; delete the 5e6 atol.
5. Site-level carbon/water/energy ledger with declared seam sources; recruitment carbon debit.
6. Everything else in ranked order.
