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
- **2026-09-08, commit 7 (the remaining energy-budget open items):**
  1. *Melting-pack residual (~0.7 J/m2 per 150 s step in `test_column_dynamics` RUN 8) -- a THIRD
     condensate time level.* `column_be_stage` debited the CAS at the Newton-converged canopy-air
     temperature inside `surface_derivs`, but valued the deposit at `t_cas1` from the flux-form commit,
     a few mK apart; with 0.04 kg/m2 of dew per step under a cold pack that is ~0.7 J/m2. Fix:
     `surface_tend_t%cond_enth` carries exactly what `surface_derivs` debited and both schemes deposit
     that number. RUN 8 worst: 0.698 -> 9.6e-7 J/m2.
  2. *Sub-threshold snowfall onto bare ground.* Now valued as ICE at min(t_3ple, tair) -- the same
     valuation `snow_accumulate` gives snow that forms a pack -- mixed with rain at the canopy-air
     temperature into one effective inflow temperature (`temp_of_liquid_enthalpy` of the mixture
     enthalpy), used for the pond inflow, the ledger and the canopy film alike. The soil then pays the
     fusion enthalpy through the plateau instead of receiving it from nowhere. New RUN 9b: the same
     cold day as snow vs as rain leaves the soil+pond poorer by ~0.56 L_f per kg (the rest is fed back
     through reduced surface losses); the old source gives exactly 0.
  3. *Restart.* `shed_water_rate` is now written and (optionally) read. Found while there: the state
     writer put `fw(:,4)` -- a column nothing ever assigned -- as `soil_w_surface_enth`, so every
     state file carried an UNDEFINED pond enthalpy (only harmful when a pond exists at checkpoint
     time); it now writes `fw(:,2)`. `test_state_roundtrip` now round-trips a ponded store with
     enthalpy and the shed rate (the pond-enthalpy check fails on the old writer).
  **Result:** one-year daily diagnostic (Jul 2074 - Jul 2075, 150 s, 4 threads): whole-column energy
  mean leak 3.6e-11 W/m2 (was 7.3e-4 before commit 6, 1.4e-7 after), worst step 1.35e-6 J/m2 (was
  87 then 0.057), every month at round-off, 0/2312640 breaches; water -4.9e-12 kg/m2. The energy
  ledger is now closed to machine precision across a full seasonal cycle including snow, freeze/thaw,
  dew, and sub-freezing precipitation.

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

---

# Items 4-6 (interfaces, duplicated helpers, naming) -- findings, 2026-09-08

Read-only audit of `main` at `10dc661` (after PR #119), three parallel passes, every ranked finding
below re-verified by hand. No code changed by this section.

## Item 4 -- module interfaces (physics-based seams)

The kernels themselves are largely clean: `veg_energy_diagnostic`, `ground_surface_fluxes`,
`snow_energy_step`, `soil_energy_*`, `soil_water_*`, `cas_column_*`, `rhizosphere_cond` and
`solve_plant_water_batch` all take physical quantities with units. The problems are one layer up.

1. **The frozen structs are the de-facto interface and mix seven seams.** `surface_frozen_t` has
   61 fields and `column_frozen_t` 50; together they carry radiation, aerodynamic conductances,
   tissue coefficients, hydrology boundary, snow outcome, root zone, plant geometry, six parameter
   structs copied out of `column_config_t`, two protocol flags (`mo_live`, `cas_condensation`) and
   ten dead or duplicated fields. Every stage kernel takes the whole thing. Proposed decomposition
   (blast radius: 6 files, incl. `rk4_oracle` and `test_column_derivs` which hand-builds the struct):
   `cas_boundary_t`, `tissue_coefficients_t`, `canopy_film_frozen_t`, `ground_boundary_t`,
   `snow_stage_t` (already exists; today copied field by field), `soil_hydrology_frozen_t`,
   `root_zone_t`; parameters passed as `ccfg`, not copied.
2. **Dead selectors and ghost fields (verified).** `leaf_energy_model` / `wood_energy_model` are
   parsed from TOML, copied into `column_config_t`, and read by NO physics routine; `LEAFEN_*`/
   `WOODEN_*` constants exist only in `use` lists; `test_column_ark` and `test_column_rk45` toggle
   `wood_energy_model` believing it switches the wood model ("RK45 PROG-WOOD" asserts on a model that
   does not change). `veg_energy_step_implicit` and `le_conductance_flux`/`lw_emission_slope` have
   zero callers in src (`le_conductance_flux` now encodes the OLD latent-only leaf payment, i.e. wrong
   physics). Frozen fields written but never read: `wood_gbh`, `wood_abs_sw/lw`, `wood_area`,
   `snow_melt_enth`, `snow_t_melt`; `src_frac` is hard-wired 1.0 and multiplied in 4 places.
   -> Delete all of it (memory rule: delete flags that gate wrong or absent physics).
3. **ARK and RK45 duplicate eight post-march blocks** (mask restore, film clamp, unpack, condensate
   deposit, tissue commit, store totals, boundary amounts, final `surface_derivs`); the RK45 mask
   restore already diverged (pond fields handled elsewhere) and the RK45 final `surface_derivs`
   (rk45:714-718) is dead work. -> shared pure helpers, extracted verbatim (see item 5 #1-#5).
4. **`column_prepass` fuses five processes** (aerodynamics, leaf gas exchange, autotrophic
   respiration, heterotrophic Rh, CAS capacities) behind 9 aggregate arguments and 18 outputs.
   -> split into physically named routines; the caller assembles `nee_biotic`.
5. **`cas_column_step_implicit` is exported but unused**; `column_be_stage` re-implements the same
   BE box inline (ark:196-201). -> call the kernel.
6. **`t_ground` is a per-evaluation input smuggled through the frozen struct**: every RHS evaluation
   deep-copies `surface_frozen_t` (~20 allocatable arrays) to overwrite one scalar (6x per RK45 step).
   -> explicit argument of `surface_derivs`.
7. **`advance_snow_stage`, `aero_bottom_to_top`, `apply_rt_forcing`, `atm_fluxes`,
   `accumulate_patch_diag`** take 4-6 aggregates for a dozen scalars each. -> scalar signatures.
8. **`leaf_gas_exchange` re-flattens ~45 PFT fields from `meds_config_t` for every leaf every
   dt_fast.** -> build `leaf_photo_params_t(n_pft)` once at config time.
9. **Three encodings of the inflow temperature** (`t_precip`, `rain_temp`, `film_u_ref`) and two
   pond implementations (scratch hydrology vs RK45 hand-composition, rk45:685-706).
   -> one `t_inflow`; `pond_update` kernel in `meds_soil_water` used by both.
10. **Soil-energy forcing assembled by hand twice** with different face conventions (ark:211-243,
    time_derivs:369-421). -> `soil_energy_forcing(...)` assembler.
11. **Reported LE/H (`atm_fluxes`) are not the ledger's fluxes**: `rho*ustar*temp2*(q-q_atm)*L_v`
    vs the CAS's `gaw*(shv1-shv_atm)` with `gah = rho*ustar*temp1` and enthalpy. The headline ET
    output and the conserved vapour export are different numbers. -> report `gaw*(...)`.
12. **Fast->slow handoff** is 12 cohort + 8 patch SoA fields, each repeated in four lockstep lists in
    core plus the restart writer; the fast GATHER also writes `site%cohort%*_water_mass` (lazy seed +
    capacity clamp). SoA is the right design (fusion needs per-field rules); the repetition is not.
    -> `cohort_fast_slice_t`/`patch_fast_slice_t` components; move seed/clamp to a slow-loop
    `reconcile_tissue_water_capacity`. Large blast radius (core, fusefiss, io); do last.
13. **Core facade incomplete**: drivers/io reach past `meds_core_interface` for `site_alloc`,
    `site_free`, `rebuild_csr`, `set_cohort_size`, `gather_pft_params`, the `DMAX_PSI_LEAF_*`
    sentinels and all of `meds_core_diag_types`; `meds_fast_types` depends on core for one sentinel.

## Item 5 -- duplicated and misplaced helpers (all verified by grep)

| # | Concept | Copies | Proposed home |
|---|---|---|---|
| 1 | Whole-column store totals (soil water, soil energy, plant water, film, tissue) | ARK 1298-1323 + RK45 754-774 verbatim; partials in soil_water/soil_energy/tests | `meds_column_stores` (pure functions, same loop order => bit-identical) |
| 2 | Post-march unpack + soil-T diagnosis + final surface_derivs | ARK/RK45 | `unpack_column_state`, `diagnose_soil_temps` in `meds_fast_types` |
| 3 | Process-mask restore | ARK 11 fields / RK45 9 fields (diverged) | `apply_process_mask` |
| 4 | Film capacity clamp + overflow/deficit | ARK/RK45 verbatim | `clamp_canopy_film` in `meds_vegetation_biophysics` |
| 5 | Condensate deposit | ARK/RK45 | `deposit_condensate` |
| 6 | LW emission slope `4*eps*sigma*T^3*A` | 5 inline + a dead helper | call `lw_emission_slope` |
| 7 | Soil-top temperature diagnosis `uext_to_temp(e(1), theta(1)*rho, hcap(1))` | 6 | `soil_layer_temp` |
| 8 | CAS<->atm conductance assembly, THREE formulations (temp1 vs temp2) | time_derivs / prepass / atm_fluxes | `cas_atm_conductances` in aerodynamics |
| 9 | Dry-air molar density `rho*(1-q)/mmdry` (tests hard-code 0.0289655) | 4 src + 5 test | `cas_molar_density` in therm_lib |
| 10 | Root-weighted mean soil T: explicit loop AND `root_weighted_psi` (misnamed, 1 caller, for temperature) | 2 | `weighted_mean` in `meds_numerics` |
| 11 | theta bounds relief with mass bookkeeping | 3 implementations | `relieve_theta_bounds` in soil_water |
| 12 | `clamp01` re-implemented inline | 14 sites, helper exists with 4 callers | call `clamp01` |
| 13 | Uniform soil (theta,T) seeding loop | driver + 4 tests | `seed_soil_column` |
| 14 | PSI_INIT tissue-water seed | driver + 3 tests | `seed_plant_water` in hydr_lib |
| 15 | Constants declared twice: `MAX_RECYCLE_YEARS` (config + met_driver), `safety = 0.9` (soil_water, control, hydraulics), `lnexp_min` (-38 in constants, shadowed by -30 in pft_params) | | one definition each |
| 16 | `debug_error` declared in 3 option types; only `energy.debug_error` is read from TOML, so the soil-water hard stops are UNREACHABLE and `co2_opts_t%debug_error` has no reader | | one `debug_error` on the column config |
| 17 | 21 test programs each define their own `check`/`check_true`; `meds_test_support` has a different signature nobody uses for numerics; three column tests hand-build the same fixture (~18 lines x3, `reset_state` x3, diurnal forcing x3) | | `check_abs`/`check_true` + `build_test_column`/`seed_column_state`/`set_diurnal_forcing` in `meds_test_support` |

Silent-omission matrix (state fields x combinators): `zero_like` (rk45) never allocates
`leaf_surf_water`/`wood_surf_water` -- safe only because the error norm excludes films; `state_err_diff`
and `zero_like` leave the pond fields default-initialised while `state_sub` subtracts them; the
test-local `copy_state` zeroes the pond on every copy. Every state field is enumerated in ~12 places
(combinators, mask restores, pack/unpack, error norm); adding a field fails nowhere at compile time.
Also: the generic `column_state_t` algebra (`state_*`, `bflux_*`, `clamp_*`) lives in `meds_fast_ark`
and RK45 imports it from there -- an RK45->ARK dependency for non-ARK code; move to the type owner.

## Item 6 -- naming and stale documentation

Top renames (all mechanical, byte-identical; Fortran is case-insensitive so use `sed -I -w`;
`bf`/`acc` have unrelated homonyms in `meds_optics_lib`):
`fro`->`frozen` (532 occ), `bio`->`patch_biophys` (299), `coh`->`column_cohort` (254; clashes with the
core `cohort` alias), `ccfg`->`column_config` (216; one letter from `cfg`), `ys`->`stage_state`,
`budg`->`budget`, `sf`->`surface_tendency`, `bf`/`acc`->`boundary_flux`/`boundary_flux_total`,
`hforc/hflux/eforc/eflux` + `chydro_*_t` -> `soil_water_forcing/flux`, `soil_energy_forcing/flux`;
fields `wcap/ccap`->`cas_mass_capacity/cas_molar_capacity`, `gah/gaw/gac`->`g_cas_atm_heat/vapour/co2`,
`hydro`(soil opts) vs `hydro_p/hydro_o`(plant) -> `soil_water_opts`/`hydraulics_params/opts`,
`h_coeff_f/g_tr_f/g_film_f` -> `_leaf`, `a_leaf/a_wood/a_store`->`*_hcap_per_dt` (`a_store` also
means carbon-to-storage in the allocator), `enth_atm`->`enthalpy_atm`, `snowf/tair/precip`->
`snowfall/air_temp/rainfall` (`forc%precip` is DOCUMENTED as "ground-reaching rainfall" but is
assigned the met rainfall), `src_frac`, `mo_live`, `coh_rnet/coh_transp`->`canopy_*`.
Mechanism-changed names: `INTEG_RK4` selects Cash-Karp RK45 (rename `INTEG_RK45`, keep accepting the
TOML string); `veg_energy_diagnostic` is the exact prognostic relaxation (rename `veg_energy_balance`);
`column_hydrology_flux` commits state (`advance_soil_water_column`); `column_prepass`
(`column_gas_exchange_prepass`); `uext_to_temp`/`temp_to_uext` (ED2 token; `internal_energy_to_temp`).
Concept-consistency: `theta` = soil moisture AND potential temperature (`theta_atm`, `can_theta`);
`*_energy` (store) vs `*_enth` (advected) is followed except `w_surface_enth`/`snow_enth0/1`;
`rain_temp` is not a rain temperature (pinned to `tsupercool_liq` under a pack); plant library uses
ED2's `bleaf/bsap/broot/bwood` against core's `*_carbon` spelling; `wood_area`/`sai` duplicate `wai`.
Keep (established, file-format bound): `agb`, `lai`, `wai`, `swe`, `ustar`, `ggnet`, `can_*`,
`gbh/gbw/gsw`, `dbh`, `nplant`, `pft`, `ipft`.
Do NOT rename netCDF registry strings or TOML keys without a compatibility note.

Stale documentation (counts in src comments): `meds_fast_split` 26 (several present it as LIVE:
rk45:4, ark:76-77, ark:1409, biophysics README:85); "split path" 38; `snow_on` 7 (3 read as if the flag
exists); `MEDS_INTEGRATOR_PARITY.md` 7 (+ a doubled path in `docs/science/numerical_scheme.md:11`);
undefined routines cited as real: `advance_hydraulics_full`, `state_wrms`, `advance_wood_energy_full`,
`veg_store_correction`; `plant_water_tendency` listed as a live kernel in two science pages;
`phenology_on` in `plant_phenology.md:177`; "C5 rejects SOIL_BC_AQUIFER" (rk45:650,686) contradicted
at rk45:497; CLAUDE.md says "8 tests"/"7/7" (21 add_test entries) and "dt_fast default 150 s" (no code
default; shipped 900 s; accuracy-limited, warn > 225 s). 153 comment lines narrate history
("USED TO", "RETIRED", "PR #"); the ten longest blocks to rewrite are listed in the audit output
(ark:1-33, fast_step:1-26, fast_snow:1-31, fast_types:299-331 and 182-212, ark:1700-1730,
rk45:622-653, control:172-198, config:445-473, veg_biophysics:107-143).

## Recommended sequence (each step verifiable byte-identical with the existing harness)

1. Deletions only: dead selectors/constants/fields/helpers (item 4 #2, item 5 #15-16), the dead RK45
   `surface_derivs`, stale comments and the ten history blocks; fix the two misdocumented facts
   (`forc%precip`, `t_ground`), CLAUDE.md counts/defaults. Zero arithmetic change.
2. Extract the eight ARK/RK45 blocks into shared pure helpers, verbatim (item 5 #1-#7, #12-#14);
   move the `column_state_t` algebra to its type owner; complete `zero_like`.
3. `t_ground` explicit; call `cas_column_step_implicit`; `soil_energy_forcing` and `pond_update`
   assemblers; `cas_atm_conductances` + fix `atm_fluxes` to report the ledger flux (this one CHANGES
   an output diagnostic, not the state).
4. Scalar signatures for the small drivers (`advance_snow_stage`, `aero_bottom_to_top`, ...).
5. Split `column_prepass`; `leaf_photo_params_t` table.
6. Mechanical renames (item 6 group 1-2) in one commit per group, `sed -I -w`, byte-identical.
7. Frozen-struct decomposition; `integrator_opts_t`; core facade completion.
8. Fast/slow slices in core (touches restart I/O) -- last.

## Execution log for items 4-6 (2026-09-08)

- **PR #120 (`refactor/review-steps-1-5`), steps 1, 1b, 2, 3 done:** dead selectors/kernels/fields/
  constants deleted; stale comments and pages corrected; `meds_column_state_ops` owns the column-state
  algebra and the shared post-march helpers; `t_ground` is an explicit argument of `surface_derivs`;
  the CAS box commit uses `cas_column_step_implicit`. Steps 4-5 NOT done (see PR for the list).
- **PR (stacked, `refactor/review-step-6-renames`), step 6 partial:** mechanical renames of the
  fast-loop argument vocabulary -- `fro`->`frozen`, `bio`->`biophys`, `coh`->`col_cohort`,
  `ccfg`->`col_config`, `budg`->`budget`, `sf`->`surf_tend`, `ys`->`y_stage`, `INTEG_RK4`->`INTEG_RK45`
  (TOML string unchanged), `veg_energy_diagnostic`->`veg_energy_balance`; 41 over-long lines wrapped
  (the 132-column rule). Shorter than the audit's proposals (`patch_biophys`, `column_config`,
  `surface_tendency`) because those pushed ~100 lines past 132 columns. Field renames (`wcap/ccap`,
  `gah/gaw/gac`, `hydro*`, `_f/_w` suffixes, `snowf/tair/precip`, `enth_atm`) and the
  mechanism-changed routine names (`column_hydrology_flux`, `column_prepass`, `uext_to_temp`) are NOT
  done.
- **Verification protocol learned:** deleting type fields or moving procedures between modules changes
  ifx code generation; outputs are then not byte-identical but hour-1 differences sit at round-off
  (2e-16 relative on state) and grow chaotically to ~1e-6 over a month. Pure renames ARE data-identical.
  Compare netCDF DATA (not bytes; headers carry a timestamp) with the first records at round-off as
  the acceptance criterion for type-touching steps.

- **nvfortran multicore built on both branches (38/38) before merging.** One trap surfaced: the moved
  `bflux_zero` assigned `acc = column_bflux_t()` to its intent(out) dummy, and nvfortran 25.11 rejects
  that in `meds_column_state_ops` ("Empty structure constructor", F-0155) although it compiled the
  identical line in `meds_fast_ark` on `main`. The line was redundant (intent(out) default-initialises,
  F2018 8.5.10) and was removed. Rule of thumb: a green ifx build does not cover a module move; build
  the NVHPC back end whenever a procedure changes module.

## Decisions after items 4-6 (2026-09-08, with the author)

- **Frozen record (step 7).** ONE container `column_frozen_t` of physically named sub-records
  (`cas_boundary_t`, `tissue_coefficients_t`, `canopy_film_capacity_t`, `ground_boundary_t`,
  `snow_stage_t`, `soil_hydraulics_t`, `root_zone_t`); parameters passed, not copied;
  `surface_frozen_t` dissolves into them. No `_frozen_` in the piece names -- "frozen" is a
  statement about lifetime, carried by the container and by `intent(in)`, and the pieces are
  reusable live (a per-stage conductance refresh would recompute `cas_boundary_t`). The record
  stays in `src/driver`: `shared/state` is for state the CORE must carry across fusion; the frozen
  record is a driver work record that depends upward on biophysics/plant/core. `meds_fast_types.f90`
  is NOT split.
- **`column_cohort_t` goes (step 8, widened).** It is a read-only view of the demographic inputs
  plus derived geometry, not state. The hand-built test views are allometrically inconsistent and
  never set `bwood` (the allocator zeroes every array but that one, so the wood heat capacity in
  those tests runs on unset memory floored to the minimum); a site built through the core's own
  allometry cannot drift like that, and the per-thread scratch exists only because the gather
  writes into it. Derived geometry (`lai`, `wai`, sapwood carbon and area, total wood carbon)
  becomes cohort-block fields refreshed after growth, consistent with `basal_area`/`agb`/`leaf_area`;
  the three hard-coded constants (leaf width 0.04 m, branch diameter 0.02 m, crown fraction 1.0)
  become PFT parameters; the driver reads contiguous CSR cohort sections. Kernels in
  `src/biophysics` never take the core type. Do this together with the fast/slow slices so the
  lockstep lists are collapsed once, and after step 7.
- **Order of the remaining work:** (1) split `column_prepass` + PFT leaf table; (2) scalar
  signatures; (3) `soil_energy_forcing` assembler, `pond_overflow` kernel, LE/H = ledger flux;
  (4) frozen decomposition + `integrator_opts_t` + core facade; (5) step 8 widened as above;
  (6) remaining field and routine renames. Extensibility note: after (5) a new per-cohort INPUT is
  one field on the cohort block; a new per-cohort PROGNOSTIC field still touches the ~12 state
  enumeration sites (silent-omission matrix) -- unaddressed.

## Execution log, PR 3 (steps 1-3 of the order above; 2026-09-08)

- `meds_fast_prepass` (new): `column_prepass` is a thin orchestrator over
  `refresh_canopy_aerodynamics`, `root_zone_environment`, `canopy_leaf_gas_exchange`,
  `canopy_maintenance_respiration`, `patch_heterotrophic_respiration`,
  `cas_capacities_and_conductances`; `aero_bottom_to_top` moved with it. Patch totals keep their
  i = 1..n order. `leaf_photo_table_t` (per-PFT `leaf_photo_params_t` + capacity ratios + solver
  selectors) is built once per run by `build_leaf_photo_table` into `column_config_t%leaf_photo`;
  `leaf_gas_exchange_batch` takes the table, `leaf_photo_params_for_pft` is the one flattening.
  `cas_atm_conductances` (aerodynamics) is the one gah/gaw/gac formula. Hour-1..3 differences vs
  main 3e-11 relative (module-move round-off), 1e-6 after a month.
- Scalar signatures: `advance_snow_stage`, `aero_bottom_to_top`, `apply_rt_forcing`,
  `accumulate_patch_diag`. Data-identical.
- `assemble_soil_energy_forcing` (state ops) for the ARK stage and the whole-column RHS; the
  provenance of faces/drainage/clip stays the caller's explicit choice. Data-identical.
- `pond_overflow` (soil water) for the scratch hydrology and the RK45 commit. ARK data-identical;
  RK45 round-off (exact inverse pair instead of e/w; POND_TINY threshold).
- LE/H reported from the ledger's b-weighted CAS->atmosphere export
  (`stage_bflux_t%atm_vap_out/atm_heat_out` -> `column_budget_t%atm_vap_export/atm_heat_export`).
  H is an explicit sensible export, NOT enthalpy-minus-latent: the CAS enthalpy values vapour at
  ~3.4 MJ/kg (liquid datum), so that difference put ~40% of LE into H (measured 5x). Output change
  confined to le/h and their derivatives: July LE mean 85.29 -> 85.69 W/m2, H 15.50 -> 15.40.
- NOT done here: `atm_fluxes` still takes the budget (fine); the `t_inflow` unification of
  `t_precip`/`rain_temp`/`film_u_ref` (item 4 #9) and `apply_process_mask` (item 5 #3) remain.
