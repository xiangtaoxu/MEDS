# MEDS P3 — The Coupled Surface Fixed Point (leaf ↔ CAS ↔ ground ↔ soil)

**Status (2026-07-08):** **P3a–P3d IMPLEMENTED** on `feature/p3-coupled-surface` (the Picard fixed point in
`column_fast_step`: leaf↔CAS, ground+soil-thermal in the loop, RT leaf-LW emission re-based to `leaf_temp`,
`[fast]` TOML config surface). Conserving, split path byte-unchanged, green on ifx Release + ifx Debug
(`-check all`) + nvfortran multicore (27/27); guarded by `test_picard_coupling`. **P3e (prognostic-leaf
option) and the P3f lagged/frozen soil-water optimization are DEFERRED** — `leaf_energy_model="prognostic"`
errors clearly (never a silent fall-back); `soil_water_coupling` LAGGED and COUPLED both re-solve the soil
water each pass (required for conservation). Two deviations from the plan below, forced by testing: the
coupled Picard **re-solves soil water + hydraulics from `state^n` every pass** (freezing them while the leaf
demand iterates leaks water/enthalpy), and `picard_relax` defaults to **0.5** (the CAS↔ground sensible
coupling makes the fixed-point map oscillatory, slope ≈ −1).

This document is the plan for P3 of the MEDS fast-biophysics loop: the simultaneous, within-sub-step solve of
the surface energy budget that today MEDS advances by a single operator-split (Gauss–Seidel) sweep with
lagged sibling temperatures.

**Scope decisions taken (user, 2026-07-08).** Two design forks are resolved as **config-selectable options,
both implemented** (not one-or-the-other):

1. **Leaf thermal model** — `leaf_energy_model = "diagnostic" | "prognostic"`. Diagnostic (zero-heat-capacity
   steady-state leaf, today's inline `tl = tcas + dtl`) is the default; a prognostic `leaf_energy` store
   (ED2/ClimaLand transient leaf via `veg_energy_balance`) is a selectable upgrade.
2. **Soil-water coupling inside the Picard loop** — `soil_water_coupling = "lagged" | "coupled"`. Lagged
   (soil water / hydraulics / `soil_evap` / `src_frac` solved once per sub-step and frozen across passes) is
   the default; coupled (re-solved every pass from the `state^n` snapshot) is the selectable tighter option.

Everything else follows the recommended mechanism below. The recommended mechanism and defaults give a
correct, ED2-faithful coupled surface with the smallest risk; the two optional paths are additive and land in
their own phases so `main` stays green throughout.

---

## 1. The problem P3 solves

The fast loop advances four families of *stores* that continuously exchange energy, water, and CO₂ every
`dt_fast` (~900 s) sub-step:

- **leaf / wood** — per cohort — `leaf_temp` (`site%cohort%leaf_temp`, `meds_demography_types.f90:82`);
- **canopy air space (CAS)** — per patch — prognostic `can_enthalpy`/`can_shv`/`can_co2`, derived `tcas`;
- **ground skin** — `t_ground = soil_temp(1)`;
- **soil thermal column** — per layer — `soil_energy → soil_temp`.

These are circularly dependent: leaf temp needs `tcas` (sensible) and net radiation, but `tcas` is set by the
leaf+ground fluxes; the ground needs `tcas` and conducts to the soil; photosynthesis sets `gs → transpiration
→ latent heat` off the leaf, and that needs the leaf temp; and **radiation ↔ energy** is circular (the
two-stream needs the leaf/ground temps to emit longwave; the energy balance needs the two-stream's absorbed
fluxes as its source — seams #8/#9 of `docs/dev_plans/MEDS_ENERGY_BALANCE_DESIGN.md` §9).

**Today (P0–P2):** `column_fast_step` (`src/driver/meds_column_dynamics.f90`) makes **one sequential sweep**
per sub-step — leaf balance (diagnostic, linearized around a *held-fixed* `tcas`) → CAS twin update → ground →
soil — with each store reading the neighbours' *previous* values. No iteration to convergence. The net-longwave
wiring feeds the two-stream `tcan_bt = tcas` (lagged) as the leaf emission temperature — this *is* the
operator-split placeholder contract.

**P3 goal:** solve `leaf_temp`, `CAS(enth/shv/co2)`, `t_ground`, and `soil_temp` **simultaneously and
implicitly** within each sub-step, so every inter-store flux is mutually consistent at the end of the step.
The physically important payoffs that only close under simultaneity:

- the **evaporative-cooling feedback** — evaporation cools the CAS and raises its humidity → the vapour
  gradient driving further evaporation shrinks → transpiration self-throttles (requires `can_enthalpy` and
  `can_shv` to advance *together* each iteration);
- the **radiation ↔ energy** loop closing to a consistent leaf/ground temperature.

The four energy kernels in `src/biophysics/meds_column_energy.f90` (`veg_energy_balance`,
`ground_surface_balance`, `canopy_air_update`, `soil_energy_flux`) were written **stateless-first** — each
takes the other stores' temperatures/fluxes as forced inputs — *specifically* so that P3 wraps them in the
coupled solve with **"only the driver stops lagging"** (`docs/dev_plans/MEDS_ENERGY_BALANCE_DESIGN.md` §6.3).

---

## 2. Recommended mechanism — an outer Picard fixed point inside `column_fast_step`

Wrap the existing sweep in an outer **Picard (successive-substitution) iteration** that re-solves the coupled
block until the store temperatures stop moving, dispatched on `cfg%integration_scheme`
(`SCHEME_SPLIT_SEQUENTIAL` | `SCHEME_PICARD_COUPLED`, both already defined in `meds_config.f90:58`).

**Why Picard over an ED2-style RK/hybrid coupled integrator.** The stateless-first kernel design *mandates*
"only the driver stops lagging." An RK/hybrid coupled solve would need a new raw-tendency layer, a hand-built
coupled (arrow/Schur) linear solver, **and** promotion of the leaf to a prognostic `leaf_energy` threaded
through the entire cohort lockstep — materially more code and a whole new "forgot-to-reallocate" bug class,
for a benefit that is second-order at `dt_fast ≈ 900 s` (the leaf thermal time constant is minutes, so the
zero-heat-capacity steady-state leaf is an *accurate* root — exactly CLM's leaf shape, ideal for Picard). We
therefore make Picard the P3 mechanism and keep an RK/hybrid `SCHEME_RK_COUPLED` only as a documented future
cross-validation reference against ED2 `integration_scheme=3` (§12), not on the P3 path.

**The load-bearing numerical idea — snapshot `state^n`, re-solve the same BE step each pass.** A *true* fixed
point must, on every pass, re-solve the *same* backward-Euler CAS/soil step **from the sub-step's initial
state `state^n`**, not accumulate. So `column_fast_step` snapshots `enth0/shv0/co20` and `soil_energy^n`
**once** before the loop and **resets** the prognostic soil store to that snapshot at the top of every pass;
only the *forcing* of the BE step (the source terms, evaluated at the current iterate's `tcas/qcas/t_ground`)
changes between passes, and the committed state at convergence is the exact implicit solution. This is why the
loop must live **inside** `column_fast_step` — wrapping externally in `run_fast_biophysics` would advance N
distinct sub-steps, or force a signature change plus snapshot/restore of the whole `bio` bundle.

**The bit-identity safety property.** `SCHEME_SPLIT_SEQUENTIAL` is exactly `niter = 1` of the same loop — one
pass, no convergence test — so the default path is **byte-identical** to today's sweep. This is the invariant
that keeps `main` green while the restructure lands, and it is a first-class test (§10).

**What is iterated vs frozen.** The iterated block is `{ leaf energy (per cohort) → ground balance → soil
thermal BE step → implicit CAS three-twin update }`. A **pre-pass** computes, once per sub-step, everything
that ED2 also freezes per DTLSM and feeds it in frozen: the aerodynamic conductances
(`aero_bottom_to_top`), leaf gas exchange + `gs` + GPP/Rd + the transpiration conductance `g_tr` and sensible
coefficient `h_coeff`, stem/root/heterotrophic respiration and `nee_biotic`, and — in the *default* lagged
soil-water mode — the soil-water column, `soil_evap`, the supply limiter `src_frac`, and plant hydraulics.
Freezing `gs` and `src_frac` (the non-smooth `min()` limiter) keeps the Picard residual smooth and the
iterated block small; the CO₂ twin is thermally passive (`nee_biotic` frozen) and is solved once after
convergence.

---

## 3. Config surface

New selectors and knobs (all defaulted so absent config = today's behaviour). Placed on `column_config_t`
(`meds_column_dynamics.f90` ~:71) for the MVP, then promoted to a `[fast]` TOML block with the presence-map
+ `validate_config` wiring in **P3d**:

| Key | Type / default | Meaning |
|---|---|---|
| `integration_scheme` | `"split"` (default) \| `"picard"` | already exists (`meds_config.f90:58`, `req_scheme`); the `error stop` at `meds_fast_loop.f90:184-186` is removed so `picard` is live |
| `leaf_energy_model` | `"diagnostic"` (default) \| `"prognostic"` | steady-state inline leaf vs prognostic `leaf_energy` store (§5) |
| `soil_water_coupling` | `"lagged"` (default) \| `"coupled"` | freeze soil-water/hydraulics per sub-step vs re-solve inside the loop (§6) |
| `picard_max_iter` | int, 20 | outer-iteration cap |
| `picard_tol_temp` | real [K], 1e-3 | temperature convergence tolerance |
| `picard_tol_shv` | real [kg/kg], 1e-6 | CAS specific-humidity convergence tolerance |
| `picard_relax` | real, 1.0 | optional under-relaxation of the next-pass seed (1.0 = none) |
| `picard_fixed_iter` | logical, false | GPU warp-uniform fixed pass count (no early exit); host-only today |

Reporting-only diagnostics on `column_budget_t`: `picard_iters` (int), `picard_worst_resid` (real),
`picard_nonconv` (int); surfaced up through `run_fast_biophysics` as an `n_picard_nonconv` out-arg (mirrors
`worst_energy`/`n_budget_fail`).

The two new selector-code families (`LEAFEN_DIAGNOSTIC/LEAFEN_PROGNOSTIC`, `SOILH2O_LAGGED/SOILH2O_COUPLED`)
live in `meds_biophysics_types` next to the existing `SOIL_*`/`ENERGY_*` codes (self-contained convention).

---

## 4. The Picard loop — structure

Inside `column_fast_step`, after the pre-pass:

```
snapshot state^n:  enth0 = cas%can_enthalpy ; shv0 = cas%can_shv ; co20 = cas%can_co2
                   soil_e_n(:) = soil_e%soil_energy(:)          ! only touched when soil is in the loop (P3b+)
niter = (scheme == PICARD) ? picard_max_iter : 1

do iter = 1, niter
   tcas_in = tcas ; leaf_in(:) = leaf_temp(:) ; tground_in = t_ground      ! previous-iterate snapshot (norm)
   reset:      soil_e%soil_energy(:) = soil_e_n(:)                          ! P3b+ (soil in the loop)
   [soil water re-solve from state^n]                                      ! P3f only (soil_water_coupling=coupled)
   leaf:       per cohort, evaluate the leaf energy balance at the CURRENT tcas/qcas   (§5: diagnostic or prognostic)
   ground:     h_ground, le_ground, g_top  at the current tcas/t_ground    (P3b+)
   soil:       one BE thermal step from soil_e_n given g_top               (P3b+)
   CAS:        src_enth = coh_h + coh_qw + h_ground + le_ground
               src_vap  = coh_transp + soil_evap
               enth1 = (wcap*enth0 + dt*(src_enth + gah*enth_atm)) / (wcap + dt*gah)     ! implicit, FROM enth0
               shv1  = (wcap*shv0  + dt*(src_vap  + gaw*shv_atm )) / (wcap + dt*gaw)
               tcas  = cas_temp_of_enthalpy(enth1, shv1) ; qcas = shv1
   converge?:  resid_T = max(|tcas-tcas_in|, max|leaf_temp-leaf_in|, |t_ground-tground_in|)   ! last two per phase
               resid_q = |shv1 - shv0_prev|
               if (scheme==PICARD .and. resid_T < tol_temp .and. resid_q < tol_shv) exit
end do

after convergence (ONCE): co21 (CO₂ twin) ; commit cas ; accumulate per-kernel + whole-column budgets ;
                          fill gpp_coh / *_resp_coh out-args ; record picard_iters/nonconv.
```

**Convergence norm.** Temperature change of the prognostic thermal stores + the CAS humidity twin. CO₂ is
excluded (thermally passive, solved once). Near-zero-LAI cohorts are **slaved** (`leaf_temp = tcas`) and
excluded from `resid_T` — the ED2 resolvable-cohort device (`rk4_integ_utils.f90:769-778` analogue) — so a
near-singular leaf denominator can never stall the loop.

**Relaxation.** `picard_relax = 1.0` default. The dominant feedback (evaporative cooling) is a *stabilizing*
negative feedback, the diagnostic leaf is an already-damped linearized root, and the CAS/soil BE steps are
L-stable with large heat capacities, so the coupled map is expected to be a contraction converging in ~2–6
passes. Optional under-relaxation `ω < 1` on the *next-pass seed* is available if oscillation appears; because
the committed `enth1` is always the exact BE solution given the converged source, relaxation never leaks a
budget (at convergence relaxed == unrelaxed).

**Non-convergence contract.** On reaching `picard_max_iter`: **clamp** `bio` to the last-converged iterate
(never write half-solved state), set `converged = .false.`, increment `picard_nonconv`, and under
`ccfg%energy%debug_error` (the existing `soil_energy_flux` discipline) `error stop`; otherwise accept the clamp
and surface the count. Fast→slow accumulators and per-patch budgets update **once** per converged sub-step,
never per pass.

---

## 5. Leaf thermal model — both options (`leaf_energy_model`)

### 5.1 `LEAFEN_DIAGNOSTIC` (default) — zero-heat-capacity steady-state leaf

Today's inline balance, promoted into the Picard kernel unchanged in arithmetic. Given the current iterate's
`tcas`/`qcas`, per cohort:

```
lw_slope = 4·leaf_emiss·σ·T_emit³·lai ; le_slope, le_ref from g_tr, qsat_c, qcas
dtl = (abs_sw + abs_lw - le_ref) / (h_coeff + le_slope + lw_slope) ; tl = tcas + dtl ; leaf_temp = tl
```

The steady-state root is exact for the leaf's minutes-scale time constant at `dt_fast ≈ 900 s`. This is the
recommended default and the P3a–P3d baseline.

### 5.2 `LEAFEN_PROGNOSTIC` — transient leaf via `veg_energy_balance` (P3e)

Uses the existing prognostic kernel `veg_energy_balance(store_energy, env, tparams, dt, is_leaf, flux)`
(`meds_column_energy.f90:160`), which takes one **L-stable linearized backward-Euler step** on a prognostic
`leaf_energy` (the design's `−8·ε·σ·T³·lai` two-sided LW damping and `ρ_air`-corrected latent term). In the
Picard loop, each pass re-steps `leaf_energy` from the snapshot `leaf_energy^n` given the current `tcas`; at
convergence `leaf_energy^{n+1}` is committed and `leaf_temp` is re-derived from it.

**New prognostic state (the one genuinely new SoA field the design adds), threaded through the centralized
cohort lockstep** (`meds_demography_types.f90`) exactly like `leaf_temp`:

- `cohort%leaf_energy(:)` (and `wood_energy(:)` if wood is prognostic) `[J/m²]` on `cohort_block`, plus the
  `patch_biophys_t` working copy.
- Lockstep sites to update (the "single place that touches every array"): `cohort_ensure_capacity`
  (allocate + init, ~:185), the `tmp%…` copy (~:248), `move_alloc_block` (~:286), `cohort_reorder` permute
  (~:356), `copy_cohort_slot` (~:409).
- **Creation stamps** — every cohort-birth site initializes `leaf_energy` from `leaf_temp` via the shared
  `temp_to_uext` inverter: `add_cohort`, `apply_recruitment`, `split_cohorts`, `apply_patch_disturbance`.
- **Patch fusion** area-blends `leaf_energy` (the `recruit_pool`/`soil_water` blend precedent).
- The leaf/wood store enters the **whole-column energy budget** (a new `Δ(leaf_energy)` term), so the budget
  test also guards the prognostic path.

Because the diagnostic and prognostic leaves share the same `env`/forcing and both feed `coh_h`/`coh_qw`/
`coh_transp`/`coh_rnet`, the Picard loop body dispatches on `leaf_energy_model` at the per-cohort leaf step and
nothing else in the loop changes.

---

## 6. Soil-water coupling — both options (`soil_water_coupling`)

### 6.1 `SOILH2O_LAGGED` (default)

Solve the soil-water column (`column_hydrology_flux`), plant hydraulics, `soil_evap`, and the supply limiter
`src_frac` **once** in the pre-pass and **freeze** them across all Picard passes. The CAS gains only the
realized water each pass. This lags the `soil_evap ↔ t_ground ↔ can_shv` coupling by the sub-step — exactly as
ED2 freezes `gsw`/hydraulics per DTLSM — and keeps the residual smooth (the non-smooth `min()` limiter is out
of the loop). Recommended default.

### 6.2 `SOILH2O_COUPLED` (P3f)

Re-solve the soil-water column + plant hydraulics + `soil_evap` + `src_frac` **inside** the Picard loop, from
a `soil_water^n` snapshot reset each pass (mirroring the `soil_energy^n` reset), so the
moisture ↔ evaporation ↔ `t_ground` ↔ `can_shv` coupling is simultaneous. Two design consequences:

- **Snapshot/reset `theta^n` each pass** alongside `soil_energy^n` — forgetting the reset double-advances the
  water column and breaks the water budget (the single highest-severity mechanical risk, caught by the
  conservation test).
- **The non-smooth supply limiter `src_frac = min(…)`** can make the residual piecewise and slow/stall Picard.
  Mitigation, in order of preference: (a) hold `src_frac`'s *active branch* fixed within a converged Picard
  cluster and only re-evaluate it between outer restarts; (b) a smoothed (softmin) limiter behind the coupled
  option; (c) rely on under-relaxation + the `max_iter` clamp. The design ships (a) first and flags (b) as a
  tunable.

Both options share the same loop skeleton; the coupled path simply moves the (already existing)
`column_hydrology_flux`/hydraulics calls from the pre-pass to the top of the pass, guarded by the selector.

---

## 7. Radiation ↔ energy coupling — lag one sub-step (ED2-faithful)

**Decision: lag the two-stream one `dt_fast`, do NOT re-run it inside the Picard loop.** `apply_rt_forcing`
(`meds_fast_loop.f90:268`) stays once per sub-step, before `column_fast_step`; the resulting
`abs_sw/abs_par/abs_lw/abs_sw_ground/abs_lw_ground` are **frozen** across all passes. This matches ED2
(`canopy_radiation` runs once per DTLSM before `odeint`; `rshort_l`/`rlong_l` stay frozen across every internal
substep). Re-running the two-stream per pass would force moving the loop up a level and re-basing the emission
temperature to the iterated `leaf_temp`, requiring removal of the `lw_slope` linearization to avoid a
double-count — a bigger, riskier reformulation for the expensive per-pass RT cost we deliberately avoid.

**Within-sub-step LW consistency** comes from the `lw_slope` Jacobian already in the leaf/ground balance: as
the temperature moves off the RT emission base during the iteration, the linearized net-LW response
`abs_lw − lw_slope·(T − T_emit)` supplies the restoring longwave without a fresh two-stream call. **Headline
invariant:** MEDS `abs_lw` is *already net* of each cohort's own emission (the two-stream folds
`emission = σ·canopy_temp⁴` into the scattering source, `meds_canopy_radiation.f90:84`), so the kernels must
**never re-subtract `σT⁴`** into the conservative tendency — only the slope correction — and the whole-column
budget's `coh_rnet` uses the same net-LW-at-`T` expression (the double-count guard the budget test enforces).

**Closing the circular seam (P3c).** Today the emission temp is substituted with `tcas` (`tcan_bt = tcas`),
which is why the inline leaf uses `lw_slope·dtl` with `dtl = tl − tcas`. P3c re-bases the emission temp to the
*previous sub-step's converged* per-cohort `leaf_temp` (`tcan_bt(j) = leaf_temp(perm(j))`) and `soil_temp(1)`,
and generalizes the leaf LW linearization base to `T_emit = leaf_temp_prev` with
`lw_slope = 4·leaf_emiss·σ·T_emit³·lai`. Because the generalized balance **reduces exactly to the current form
when `T_emit == tcas`**, the split path (which still feeds `tcas`) is provably bit-identical — which is why the
seam re-base is safe to stage as its own phase. Across sub-steps the one-step lag ties `leaf_temp(n) → RT
emission(n+1)`; at `dt_fast ≈ 900 s` and a minutes leaf time constant this is tight. If the RT-vs-energy
`leaf_temp` mismatch ever exceeds tolerance, the deferred fallback is a cheap linearized LW correction
`Δabs_lw = −4·ε·σ·T³·ΔT` inside the loop while keeping the full two-stream lagged.

---

## 8. Numerics, conservation, stability (summary)

- **Every pass restores the prognostic soil (and, in coupled mode, water) store to `state^n` and re-applies
  boundary advection from that snapshot**, and the CAS is re-solved from `enth0` each pass, so the committed
  state at convergence is the mutually-consistent implicit solution and the per-kernel budgets close by
  construction each pass. Whole-column energy/water budgets close to round-off, accumulated **once** after
  convergence — including the `advect_soil_heat = .true.` re-run.
- **Stiffness:** the diagnostic leaf denominator `h_coeff + le_slope + lw_slope` is the damped steady-state
  root; the prognostic leaf uses the L-stable linearized BE step; CAS/soil BE steps are L-stable. Near-zero-LAI
  slaving removes the only near-singular denominator.
- **GPU:** host-only today, so the adaptive `picard_iters` count is fine; `picard_fixed_iter` runs a uniform
  pass count for warp uniformity when `column_fast_step` is later offloaded.
- **nvfortran issue #7:** the new snapshot/reset code must bind array-valued function results to named
  temporaries before passing them into calls (green ifx is not sufficient — build nvfortran multicore).

---

## 9. Phasing (each phase independently mergeable + tested; `split` stays the default throughout)

- **P3a — Picard scaffold + leaf↔CAS fixed point.** Remove the `error stop`; dispatch on `integration_scheme`
  in `column_fast_step`; split the cohort loop into pre-pass (frozen) + inner leaf-energy re-eval; snapshot
  `enth0/shv0/co20`; iterate leaf(diagnostic)↔CAS to convergence; near-zero-LAI slaving; non-convergence clamp.
  Ground/soil stay once-after-convergence. **Test:** split bit-identity, Picard convergence (≤6 iters),
  evaporative-cooling feedback (Picard vs split@1).
- **P3b — ground + soil thermal inside the loop.** Move `ground_surface_balance` + the soil thermal BE step
  inside; `soil_energy^n` snapshot/reset each pass. **Test:** ground/soil self-consistency, conservation under
  coupling (incl. `advect_soil_heat`).
- **P3c — radiation↔energy seam re-base (lagged one sub-step).** `tcan_bt = leaf_temp_prev`,
  `T_emit`-generalized LW. **Test:** one-step-lag identity, `|leaf_temp_RT − leaf_temp_energy|` bound, budget
  still closes (double-count guard), bare `ncoh=0` patch continuous.
- **P3d — hardening + config surface.** `[fast]` TOML block (presence-map + `validate_config`), non-convergence
  contract surfaced, `picard_fixed_iter`. **Test:** config round-trip, non-convergence contract, `test_fast_loop`
  site-level Picard path.
- **P3e — `LEAFEN_PROGNOSTIC` option.** `leaf_energy`/`wood_energy` SoA + full cohort lockstep + creation
  stamps + patch-fusion blend + whole-column budget term; dispatch the leaf step on `leaf_energy_model`.
  **Test:** prognostic-vs-diagnostic agree at quasi-steady state, prognostic leaf lags diagnostic under a fast
  radiation ramp, budget closes with the new store, lockstep survives fuse/split/sort/disturbance.
- **P3f — `SOILH2O_COUPLED` option.** Re-solve soil water + hydraulics inside the loop from `theta^n`;
  `src_frac` active-branch freeze within a cluster. **Test:** coupled-vs-lagged agree in well-watered
  conditions and differ (sign-correct) under strong drought stress; water budget closes; convergence not
  degraded past `max_iter`.

Recommended merge order: P3a → P3b → P3c → P3d (the correct, ED2-faithful coupled surface with the
recommended defaults), then P3e and P3f as additive selectable options in either order.

---

## 10. Test plan (consolidated)

1. **Split bit-identity** — `test_column_dynamics` diurnal case under `split` vs `picard@max_iter=1`: `cas`,
   `soil_e`, `leaf_temp` trajectories bit-identical every step; plus a golden compare against pre-P3 `main` to
   catch round-off drift from the loop restructure.
2. **Picard convergence** — `picard_nonconv == 0` over 96 steps, `picard_iters ≤ ~6` (contraction).
3. **Evaporative-cooling feedback** (the simultaneity value proof) — a high-VPD / water-stressed step yields
   lower transpiration and a measurably different `can_shv`/`tcas` than split@1, in the correct direction; at
   midday `|tl − (tcas + dtl(tcas))| < tol` under Picard but not under split.
4. **Ground/soil self-consistency** — `|Δsoil_temp(1)| < tol` at convergence; soil-surface diurnal temp differs
   from split by a bounded, sign-correct amount.
5. **Conservation under coupling** — CAS energy/water/CO₂ + soil-thermal + soil-water + whole-column budgets
   `n_fail == 0` on the Picard path (incl. `advect_soil_heat`).
6. **RT↔energy consistency** — one-step-lag identity; `|leaf_temp_RT − leaf_temp_energy|` within a few tenths K
   at diurnal quasi-steady; whole-column energy budget still closes with `T_emit = leaf_temp` (LW double-count
   guard); bare `ncoh=0` patch continuous.
7. **Non-convergence contract** — pathological forcing → `picard_nonconv` increments, state clamps to the
   last-converged iterate with budgets closed in production, `error stop` under `debug_error`.
8. **Stiff-case + relaxation** — large LAI + small `can_depth` converges within `max_iter`; per-pass `max|ΔT|`
   monotone for default `ω`; `picard_relax < 1` restores convergence if it stalls; below-`veg_hcap_min` cohort
   slaved.
9. **Prognostic leaf (P3e)** — diagnostic vs prognostic agree at quasi-steady; prognostic lags under a fast
   ramp; budget closes with the `leaf_energy` term; lockstep survives fuse/split/sort/disturbance.
10. **Coupled soil-water (P3f)** — coupled vs lagged agree well-watered, differ sign-correctly under drought;
    water budget closes; convergence preserved.
11. **Config round-trip (P3d)** — `[fast]` keys parse / presence-map / `validate`; `split` bit-identical.
12. **Cross-compiler (all phases)** — nvfortran multicore green on every touched module, CPU-vs-multicore
    identical; issue-#7 array-temp discipline in the snapshot/reset code.

---

## 11. Risks

- **Bit-identity regression on the default `split` path** from restructuring the sweep into a `niter=1` loop.
  *Mitigation:* `iter==1` guards reproduce today's order exactly; split bit-identity + golden compare.
- **Snapshot/reset correctness of the prognostic soil store (P3b) and water column (P3f)** — the single
  highest-severity mechanical risk; forgetting the reset double-advances and breaks conservation. *Caught by
  the conservation test.*
- **LW double-count when re-basing emission temp `tcas → leaf_temp` (P3c)** — `abs_lw` is already net; only the
  slope correction may be applied. *Guard:* exact reduction at `T_emit == tcas` + budget closes.
- **Prognostic-leaf lockstep omission (P3e)** — a new SoA field missed at one of the reorder/copy/alloc/creation
  sites (the ED2 "forgot-to-reallocate" class). *Mitigation:* the lockstep checklist in §5.2 + fuse/split/sort
  test.
- **Non-smooth `src_frac` stalling coupled Picard (P3f).** *Mitigation:* active-branch freeze / softmin /
  under-relaxation + clamp.
- **Oscillation / non-convergence on stiff steps** (small `can_depth`, huge LAI, dawn radiation sign flips).
  *Mitigation:* under-relaxation, `max_iter` cap, near-zero-LAI slaving, clamp-or-error contract.
- **RT one-sub-step lag** is an approximation to the true radiation↔energy fixed point; adequate at diurnal
  cadence, flagged, with the per-iter linearized-LW correction as fallback.
- **GPU:** host-adaptive iteration count breaks warp uniformity; switch to `picard_fixed_iter` at offload.

---

## 12. Open questions (remaining, for the domain expert)

- **`src_frac` handling in `SOILH2O_COUPLED`** — active-branch freeze vs a smoothed softmin limiter: pin the
  choice against ED2 with real BCI/ERA5-Land drought forcing (the coupled path's whole reason to exist).
- **Convergence norm scope** — temperature stores + `|Δcan_shv|`, or is the temp norm alone sufficient? Confirm
  via the evaporative-feedback test; keep `picard_tol_shv` active until then.
- **Relaxation default** `ω = 1.0` (assume contraction) vs a conservative `ω < 1` — determine empirically from
  the stiff-case test.
- **Aerodynamic conductances** — keep lagged (recommended, ED2-like) or add a CLM-style outer Monin–Obukhov
  pass around the temperature iteration? Only if stability oscillation is observed.
- **GPU `picard_fixed_iter` count** that reproduces the host-adaptive result to tolerance across the diurnal
  cycle, and whether it scales with `dt_fast` — pin before the offload path is trusted.
- **Retire `SCHEME_SPLIT_SEQUENTIAL`** in favour of `picard@max_iter=1` (shared code path), or keep both
  selectors? *Recommendation:* keep both selectors, one code path.
- **`SCHEME_RK_COUPLED` cross-validation integrator** (ED2 `integration_scheme=3` analogue) — build later as a
  reference/escape hatch, or drop? Out of scope for P3.
