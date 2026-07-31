# MEDS Fast-Loop Numerics Overhaul: IMEX-ARK Implementation Plan

> **SUPERSEDED IN PART — 2026-07-31.** Three things in this document no longer hold:
> 1. **The `split` integrator is RETIRED** and `[fast].integration_scheme` is deleted. `ark` is the
>    default; `rk45` is the accuracy baseline and its stiff rescue now redoes the step on `ark`.
>    Anything here that treats split as the default, the reference, or a comparison anchor is history.
> 2. **"IMEX-ARK" is a misnomer.** The biotic CO₂ source is folded implicit, so `f_E == 0` and the
>    scheme is a 2-solve **ESDIRK2** (γ = 1 − 1/√2). The config string stays `"ark"`.
> 3. **`dt_fast` is STABILITY-limited, not accuracy-limited.** Above ~150–225 s the frozen surface
>    coupling drives a sustained period-2 canopy-air oscillation (~8 K at 900 s) that **no
>    conservation ledger detects**. Default is now 150 s. Every measurement in this document taken at
>    900 s or 1800 s is inside that regime.
>
> Current state: `docs/science/numerical_scheme.md` §2–3 and
> `docs/dev_plans/MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §9–14.

*Prepared for the MEDS lead developer. This synthesizes the tendency-coupling map, the removal/plumbing inventory, the fast-slow-seam-and-test map, the IMEX-ARK design memo, and the adversarial review into one buildable plan. Where the review found a real hole in the memo, the plan adopts the review's correction — three memo inversions are fixed below (§2, §5). No code is changed here.*

---

## 0. Goal & scope

**Brief:** overhaul and simplify the numerical scheme of the MEDS *fast* processes — remove the operator-split + Picard machinery and implement one coupled IMEX additive Runge–Kutta (ARK) integrator. Slow processes stay as simple daily updates.

| Layer | Status | Action |
|---|---|---|
| **Fast loop** (`column_fast_step`, sub-`dt_fast`, stiff enth/shv/CO₂/soil-thermal/soil-water/ψ) | operator-split (`niter=1`) **or** under-relaxed Picard (`niter>1`), plus **3 internal adaptive substeppers** | **REPLACED** by one IMEX-ARK stage sequence + one column-level controller |
| **Fast→slow seam** (four `*_accum` carbon totals `[kgC/plant/day]`, `meds_demography_types.f90:90-96`) | plain left-Riemann `Σ flux·Δt·unit` over sub-steps (`meds_fast_loop.f90:287-293`) | **UNTOUCHED** — confirmed the slow loop consumes only daily scalars (`meds_vegetation_dynamics.f90:100-108`); no integrator sophistication wanted here |
| **Slow loop** (vegetation dynamics, daily) | simple daily Euler/exact updates | **UNTOUCHED** — explicitly out of scope |
| **RK4** | does not exist | **ADDED** as a non-production reference oracle (`INTEG_RK4`) for order/stiffness tests only |

**Scope boundary (load-bearing):** the rewrite lives *entirely inside* `column_fast_step`. Its output contract — four optional flux out-args (`gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh`, `meds_column_dynamics.f90:177-180`), the `column_budget_t budg` struct, and the `intent(inout) bio` advance — is preserved verbatim. `run_fast_biophysics` (`meds_fast_loop.f90:282-293`) and the accumulator seam require **zero edits** beyond dropping one `SCHEME_PICARD_COUPLED` branch (`:427`).

**Honest headline (see §8):** this net-removes *methods* and *config knobs* (fewer integrators, no outer Picard, no 3 substeppers, 7 config fields × 3 sites deleted) but net-adds *numerical-linear-algebra machinery* (an arrowhead/Schur surface solve that does not exist today). Bill it as **"one coupled solver + one controller,"** not "less code."

---

## 1. Target architecture

### 1.1 The four pieces

1. **`column_derivs`** — a `pure` stateless RHS `f(y)` returning `dy/dt` for every reservoir, committing nothing, evaluated at an arbitrary stage state `Y_i` with **current-state (explicit) face/flux coefficients**. This is the central refactor: every current kernel is *advance-and-commit* (form flux → hide in a BE denominator → mutate `intent(inout)`); the ARK needs the pre-denominator numerator.
2. **The implicit-stage arrowhead solve** — `(I − γ·dt·J)·Δ = r`, one factorization structure per stage, reusing `thomas_solve` (soil columns) + a small dense Schur border (CAS + soil-top + condensed leaf/ground).
3. **The ARK tableau** — `ARK2(1)3L[2]SA` (Kennedy–Carpenter): 3 stages, 2nd-order solution + embedded 1st-order estimate, single `γ = 1 − 1/√2 ≈ 0.29289322`, stiffly accurate, L-stable.
4. **One column-level adaptive controller** — a single WRMS error norm over the whole coupled state + one accept/reject, reusing `adaptive_step_update` (`test_numerics.f90:85-93`). This is what *subsumes* the three deleted internal substeppers.

### 1.2 New fast step vs old split (text diagram)

```
OLD (meds_column_dynamics.f90:323-458):
  freeze pre-pass (:237-314)          [KEEP]
  snapshot y^n (:316-321)             [KEEP enth0/shv0/co20; DELETE store snapshots]
  do iter = 1, niter                  ── niter=1 => Lie split;  niter>1 => Picard
     leaf-T diagnostic @ old CAS  (:336-366)
     reset stores to y^n; soil-water BE (:368-400)
     ground skin              (:402-405)
     CAS 3-twin BE-in-atm     (:407-412)
     reset soil_e; soil-thermal BE (:414-434)
     convergence test / under-relax (:436-457)
  end do                              [DELETE ENTIRE SHELL]
  commit CAS/CO2 (:461-464); budgets (:480-507)   [KEEP]

NEW:
  freeze pre-pass (:237-314)          [REUSE VERBATIM  -> fast_frozen_t]
  snapshot y^n                        [REUSE mechanism; RK-stage base]
  ARK loop over dt_ark sub-steps of dt_fast:
     for stage i in 1..3:
        Y_i  = y^n + dt*Σ_j (a^E_ij*f_E(Y_j) + a^I_ij*f_I(Y_j))   [explicit part known]
        call column_derivs(Y_i, fro, ...) -> f, diag  [PURE, Step A/B/C below]
        SOLVE implicit block for stage increment:
           - CAS{enth,shv} + soil-top-E + soil-water-supply : dense Schur border
             with leaf-T(n cohorts) + ground skin STATIC-CONDENSED
           - soil-heat interior : Thomas column (top row joins border)
           - soil-water interior : Thomas column (Celia, capped inner Picard)
           - hydraulics 2x2 per cohort : plain ESDIRK (I - γ dt M)^-1   [ONE-WAY, after block]
           - CO2 : scalar divide
        active-set pass on src_frac (R1); single-flux-per-interface commit (R2)
     accept/reject on WRMS(y_2 - ŷ_1); adaptive_step_update
  commit; b-weighted boundary-flux budget ledger (§6)   [KEEP budget code]
```

### 1.3 `column_derivs` signature and interior

```fortran
pure subroutine column_derivs(y, fro, ccfg, geom, f, diag)
   type(column_state_t),  intent(in)  :: y     ! enth,shv,co2,soil_E(:),theta(:),psi(:, :)
   type(fast_frozen_t),   intent(in)  :: fro   ! §1.2 frozen pre-pass + forcing
   type(column_config_t), intent(in)  :: ccfg
   type(column_geom_t),   intent(in)  :: geom  ! dz, dz_node, root_frac, ncoh, nsl
   type(column_tend_t),   intent(out) :: f     ! dy/dt for every component; commits nothing
   type(column_bflux_t),  intent(out) :: diag  ! boundary fluxes for the ledger (§6)
end subroutine
```

- **Step A — diagnose from `y`:** `tcas = cas_temp_of_enthalpy(enth,shv)`; per soil layer `uext_to_temp → soil_temp(k),fliq(k)`, `soil_psi_of_theta → psi_soil(k)`, and `K(θ),C(ψ)` from the constitutive curves. `t_ground=soil_temp(1)`, `t_bot=soil_temp(n)`.
- **Step B — algebraic (capacity-free) constraints** solved *at this stage state* (index-1 DAE rows; τ_leaf≈7 s, ground algebraic ≪ `dt_fast`): the diagnostic leaf-T `dtl_i` (`meds_column_dynamics.f90:352-354`) and the ground skin `g_top = abs_sw_ground + abs_lw_ground − h_ground − le_ground` (`meds_column_energy.f90:240-250`).
- **Step C — assemble tendencies** in flux-divergence / capacity form (see §2 table for the seven rows).

---

## 2. The IMEX partition

**Reviewer-corrected from the memo.** Two inversions are fixed: **(i) soil-water joins the coupled surface block** (partly) because the `src_frac` supply limiter makes soil-water↔CAS a *two-way* coupling; **(ii) hydraulics leaves the block** — `gs`/`ψ_leaf` is frozen/lagged one step (`meds_column_dynamics.f90:274`), so hydraulics is a **one-way downstream** 2×2, not part of the dense arrowhead.

### 2.1 Partition table

| Reservoir / coupling | State (units) | τ | ARK class | Solver (reused) |
|---|---|---|---|---|
| **R1 CAS enthalpy** | `can_enthalpy` [J/kg] (`state_types.f90:47`) | atm 150–240 s | **implicit — coupled surface block** | dense Schur border |
| **R2 CAS shv** | `can_shv` [kg/kg] (`:48`) | leaf/DSL | **implicit — coupled surface block** | dense Schur border |
| **R7 leaf-T** (per cohort) | `leaf_temp` [K] (`:388` working) | ≈7 s | **implicit — static-condensed into block** | closed-form elim (slopes `h_coeff_f,le_slope,lw_slope` `:353`) |
| ground skin | — (no capacity) | algebraic | **implicit — static-condensed into block** | closed-form elim (`ggnet·ρ·cp` `:403`) |
| **R4 soil-top-E (layer 1)** | `soil_energy(1)` [J/m³] | ~4000 s seam | **implicit — joins block** (arrowhead top row via `g_top`) | Thomas top row ↔ border |
| **R4 soil-heat interior** | `soil_energy(2:n)` | 4e3–1e6 s | **implicit — tridiagonal** | `soil_heat_be_step` (`energy:137-148`) |
| **R5 soil-water supply/evap** | (coupling to R1/R2) | 100–540 s wet | **implicit — in block via src_frac active-set** | Schur + active-set (§5) |
| **R5 soil-water interior** | `theta(1:n)` [m³/m³] (`:32`) | 100–540 s | **implicit — tridiagonal (Celia, capped inner Picard)** | `soil_be_single_step` (`hydro:304-380`) |
| **R6 hydraulics 2×2** (per cohort) | `psi(NODE_LEAF/WOOD,i)` [MPa] | 17 s | **implicit — own 2×2, ONE-WAY after block** | plain ESDIRK `(I−γdt·M)⁻¹` (NOT ETD — see §5/R6) |
| **R3 CAS CO₂** | `can_co2` [µmol/mol] (`:49`) | atm | **implicit — standalone scalar** | scalar divide (`:462`) |
| biotic GPP/Ra/Rh → CO₂ | frozen (`:276-305`) | — | **explicit** (once/step, zero Jacobian) | frozen source |
| advective soil-heat faces `qwf` | opt-in (`:428`) | — | **explicit** | frozen face |
| ponding/aquifer (R8), leaf-film (R9) | `w_surface,w_aquifer,leaf_water` | non-stiff | **explicit** buckets | forward update |
| κ(θ,fliq)/C(ψ), K(θ) coefficients | — | — | **explicit** (recompute per stage) | constitutive curves |

### 2.2 Stability proof — no stiff coupling left explicit

Explicit RK real-axis cap ≈ **19 s**. Every coupling with τ at or below that is inside an implicit block:

- **leaf↔CAS↔ground↔soil-top surface energy-balance fixed point** (τ 7–4000 s): fully implicit, static-condensed into the border. The review derived the CAS self-feedback through algebraic slaves `d(coh_h)/dtcas ≈ −h_coeff_f·(le_slope+lw_slope)/denom` — a large negative feedback with λ·Δt ≈ 900/7…900/80 ≫ 19 s. **Any explicit treatment of `src_enth`'s enth/shv dependence is unconditionally unstable** (R3). The memo's "or ARK-explicit" option for this seam is **rejected**.
- **CAS↔atm self-relaxation** (150–240 s): implicit, L-stable (the BE-in-atm form, generalized to the γ-stage form per R9).
- **Hydraulics stiff modes** (leaf-node 17 s; rhizosphere `rhz·(ψ_soil−ψ_w)`): the internal 2×2 is implicit (eigenvalues ≤ 0 ⇒ L-stable); the rhizosphere term is implicit *on the hydraulics side* (`m22,c2`, `hyd:337,339`). Only the transpiration *magnitude* passes between hydraulics and soil, and that is bounded by `src_frac ∈ [0,1]` — a contractive coupling, amplification ≤ 1 regardless of dt.
- **Left explicit, provably safe:** (i) biotic sources — frozen, zero Jacobian, no stability constraint; (ii) soil-evap→shv — DSL resistance makes `∂evap/∂qcas` small, τ ≫ dt; (iii) leaf-demand→soil-sink direction — a per-stage source whose back-direction (`src_frac`) is contractive.

**Gauss–Seidel-within-step subsumes Picard:** the ESDIRK's 3 stages, evaluated in block order (surface+demand → soil water with that demand → soil ψ feeds next stage's demand), exchange information ≈3 times per step — the same information-exchange count the old outer Picard used, now embedded in a high-order tableau.

---

## 3. What gets REMOVED

Every guarantee the deleted code carried is re-listed with where the ARK now provides it.

### 3.1 `src/driver/meds_column_dynamics.f90` — Picard outer iteration + split sweep

| Lines | What | ARK preserves it via |
|---|---|---|
| `:44` | `use … SCHEME_SPLIT_SEQUENTIAL, SCHEME_PICARD_COUPLED` | import removed |
| `:87-96` | `column_config_t` Picard knobs (`picard_max_iter/tol_temp/tol_shv/relax/fixed_iter`, `leaf_energy_model`, `soil_water_coupling`) | replaced by `ark_rtol/ark_atol/ark_adaptive/ark_fixed_substep` (§4) |
| `:99-103` | `LEAFEN_*`, `SOILH2O_*` param constants + `public` | leaf stays diagnostic (sole model); soil-water always re-solved in-stage |
| `:137-140` | `column_budget_t` Picard diagnostics (`picard_iters/nonconv/worst_resid`) | rename → `ark_substeps/ark_rejects` (or keep names as opaque counters to minimize test churn) |
| `:181-182` | `converged/iters` optional out-args | **KEEP** — natural home for ARK substep count / accept flag (preserves `test_picard_coupling`/`test_fast_loop` reporting) |
| `:198-209` | Picard-only locals (`t_emit,te,soil_w_n,soil_e_n,psi_n,resid_T,iter,niter,…`) | deleted |
| `:219-220` | `picard = (scheme==SCHEME_PICARD_COUPLED)`; `niter` setup | scheme dispatch gone (replaced by `integrator` dispatch) |
| `:224-225` | `LEAFEN_PROGNOSTIC` guard `error stop` | diagnostic leaf is sole path |
| `:316-321` | state^n snapshot | **KEEP** `enth0/shv0/co20` (RK-stage base); **DELETE** `soil_w_n/soil_e_n/psi_n` store snapshots (the ARK re-snapshots stores at stage entry — the *mechanism* survives, the Picard *reset-each-pass* semantics do not) |
| `:323-458` | **entire `do iter` fixed point** (leaf 3a `:336-366`, LAI-slave `:342-344`, LW base `:348`, soil-water reset+resolve `:368-400`, ground `:402-405`, CAS 3-twin `:407-412`, soil-thermal reset+step `:414-434`, convergence `:436-448`, under-relax `:449-457`) | **Guarantees carried:** (i) CAS L-stable-in-atm → generalized γ-stage form (§5/R9); (ii) realized transp = root uptake via `src_frac` → active-set in-stage (§5/R1); (iii) water/enthalpy conservation across leaf↔soil↔CAS seam → single-flux-per-interface commit (§5/R2, §6); (iv) leaf LW counted once at leaf_temp → in Step B |
| `:459` | `if (picard .and. picard_fixed_iter) nconv=.true.` | deleted |
| `:466-474` | Picard diagnostics writeback + non-convergence contract | replaced by ARK substep-controller non-convergence contract (flag, no device `error stop`) |

**Net:** keep the pre-pass (`:237-314`), `enth0/shv0/co20` capture, CAS/CO₂ commit (`:461-464`), all budgets (`:480-507`); excise the `do iter` shell, collapse its body into one ARK advance.

### 3.2 `src/driver/meds_fast_loop.f90` — Picard LW-emission branch

| Lines | What | Action |
|---|---|---|
| `:18` | `use meds_config, only: …, SCHEME_PICARD_COUPLED` | drop `SCHEME_PICARD_COUPLED` |
| `:427-431` | `apply_rt_forcing`: `if (scheme==SCHEME_PICARD_COUPLED) tcan_bt(j)=bio%leaf_temp(imin) else tcan_bt(j)=tcas` | collapse to `tcan_bt(j)=tcas` (diagnostic leaf linearizes LW around `tcas`). This is the only Picard branch in the site driver; the `column_fast_step` call `:266-294` and accumulators `:287-293` are unchanged |

### 3.3 Internal adaptive substeppers → folded into the column controller

| File / routine | Lines to delete | Guarantee → ARK |
|---|---|---|
| `meds_column_hydrology.f90` `soil_water_advance` | `:231-297` wrapper: `SOIL_SUBSTEP_FIXED` loop `:253-268` + adaptive step-doubling `:270-296` (call site `:143-144`) | L-stable BE + step-doubling + mass-conservative flux-form → **KEEP `soil_be_single_step` `:304-380`** (called once/stage, with capped 2–3-iter inner Picard, §5/R8); accuracy control → one column WRMS norm incl. θ (§5/R8). Dead after fold: `soil_opts_t%{nsub,h_init,max_substep,substep}` (`biophysics_types.f90:157,161-162`), `SOIL_SUBSTEP_ADAPTIVE/FIXED` (`:119-120,:38`) |
| `meds_plant_hydraulics.f90` `solve_plant_water` substep | `:256-296`: `HYDRO_SUBSTEP_FIXED` loop `:258-266` + adaptive `:267-296` | exact/unconditionally-stable → **KEEP `freeze_coeffs`/`apply_expm` `:328-371`** but as the *reference* path only; production hydraulics uses plain ESDIRK 2×2 (§5/R6); ψ enters column norm (§5/R7). Dead after fold: `hydro_opts_t%{substep_mode,h_init,max_substep}` (`plant_types.f90:149,153-154`), `HYDRO_SUBSTEP_ADAPTIVE/FIXED` (`:113-114`) |

### 3.4 Config plumbing

| File | Lines | Action |
|---|---|---|
| `meds_config.f90` | `:28` public list; `:56-59` `SCHEME_*` constants + comment; `:72` field `integration_scheme`; `:74-82` Picard/leaf/soil-water block; `:228-230` validation | replace with `INTEG_IMEX_ARK`/`INTEG_RK4` constants, field `integrator`, `ark_rtol/ark_atol/ark_adaptive/ark_fixed_substep`, integrator-range + `ark_rtol>0`/`ark_atol>0` validation (inside `if (cfg%fast_biophysics_on)`) |
| `meds_config_io.f90` | `:21` import; `:166-180` **`req_scheme` subroutine**; `:477` call; `:478-488` defaulted Picard reads | replace `req_scheme`→`req_integrator` (`"imex_ark"→INTEG_IMEX_ARK` default, `"rk4"→INTEG_RK4`); call `req_integrator(tm,'fast.integrator',cfg%integrator,miss)`; add defaulted `ark_*` reads (all defaulted → existing configs keep running) |
| `meds_test_support.f90` | `:13` import; `:33` `cfg%integration_scheme = SCHEME_SPLIT_SEQUENTIAL` | set new integrator default (`INTEG_IMEX_ARK`) |
| `test/test_picard_coupling.f90` | entire program | **repurpose** as RK4-vs-ARK reference cross-check (it is the only consumer of `SCHEME_PICARD_COUPLED`, `picard_max_iter`, `picard_fixed_iter`, `budg%picard_*`); repoint golden anchor (see §7) |

**TOML surface:** `[fast].integration_scheme = "split"|"picard"` → `[fast].integrator = "imex_ark"|"rk4"`; `[fast].picard_*`/`leaf_energy_model`/`soil_water_coupling` retired; `[fast].ark_rtol/ark_atol/ark_adaptive/ark_fixed_substep` added.

---

## 4. What gets ADDED / REFACTORED

File-by-file, MEDS naming conventions (spell out unconventional acronyms; keep `dbh/nplant/pft`; no ED2-style contractions).

| File (new / edited) | Content | Reuses |
|---|---|---|
| **NEW** `src/driver/meds_column_derivs.f90` | `pure subroutine column_derivs` (§1.3); types `column_state_t`, `column_tend_t`, `column_bflux_t`, `column_geom_t`, `fast_frozen_t`. Step-A diagnose, Step-B algebraic leaf/ground, Step-C tendency assembly. **No commits, no BE denominators, current-state coefficients.** | pulls the algebra out of `meds_column_dynamics.f90:337-434` and `meds_column_energy.f90:240-250` |
| **NEW** `src/driver/meds_ark_stepper.f90` | ARK2(1)3L[2]SA tableau (`a^E,a^I,b,b̂,c,γ`); the sub-step loop; per-reservoir stage storage `k_i`; WRMS error norm; accept/reject; fixed-substep lockstep branch (§5 GPU). Calls `column_derivs` + the implicit stage solve. | `adaptive_step_update` (`meds_numerics`) |
| **NEW** `src/driver/meds_ark_stage_solve.f90` | The implicit stage: assemble `(I−γdt·J)`; Thomas-factor soil-heat + soil-water columns; assemble + closed-form-invert the dense Schur border (CAS{enth,shv} + soil-top-E, with leaf-T(n) + ground condensed); `src_frac` active-set pass; per-cohort hydraulics 2×2 ESDIRK; CO₂ scalar divide; single-flux-per-interface commit. | `thomas_solve`, `soil_be_single_step`, `soil_heat_be_step`, the CAS BE form (γ-generalized) |
| **EDIT** `src/driver/meds_column_dynamics.f90` | `column_fast_step` keeps pre-pass + snapshot + commit + budgets; body = pack `y` → `call ark_stepper` → unpack → b-weighted ledger. Delete the `do iter` shell (§3.1). | — |
| **EDIT** `src/driver/meds_fast_loop.f90` | `:427-431` collapse; `build_fast_context` `:105-112` swap 7 Picard copies for `ark_*` copies. | — |
| **EDIT** `src/shared/meds_config.f90`, `src/io/meds_config_io.f90` | §3.4 config surface. | — |
| **NEW (reference path)** RK4 dispatch | `INTEG_RK4` case in `column_fast_step` dispatch (where split/picard dispatched, `:219`): a fully-explicit adaptive RK4 over `column_derivs` (hydraulics via `apply_expm` reference path). Non-production oracle for §7. | `column_derivs`, `apply_expm` |

**State vector pack/unpack (§ MAP part 3):** per-patch `bio%cas%{can_enthalpy,can_shv,can_co2}` + `bio%soil_e%soil_energy(1:n_active)` + `bio%soil_w%theta(1:n_active)`; per-cohort `bio%psi(NODE_LEAF:NODE_WOOD, 1:ncoh)`. Invariants: `N_HYDRO_NODE (state_types.f90:26) == N_HYDRO (plant_types.f90:94)`; `n_active ≤ n_soil_layer_max=20`. `can_temp`, `soil_temp/soil_fliq`, `leaf_temp` are **diagnosed**, not integrated.

**PRESERVE / REUSE VERBATIM** (the ARK's building blocks): `thomas_solve` (`hydro:348`), `soil_be_single_step` (`hydro:304-380`), `soil_heat_be_step` (`energy:44/:126`), `apply_expm`/`freeze_coeffs`/`expm_step` (`hyd:328-371`, reference path), the frozen pre-pass (`:256-314`), the budget checks (`:480-507`, `track_resid :558-565`), `aero_bottom_to_top` (`:516-555`), the RT-forcing contract (`:386-462`, minus the `tcan_bt` collapse), `src_frac` (`:384-400`).

---

## 5. Handling the hard parts (reviewer mitigations folded in)

### R1 (P0-FATAL) — `src_frac` complementarity → conservation break

`src_frac = min(1, uptake_total/coh_transp)` (`:385`), then `coh_qw *= src_frac; coh_qsoil *= src_frac; coh_transp *= src_frac` (`:398-400`). The old code re-solved soil water from state^n each Picard pass *precisely so realized uptake = realized CAS latent* (comment `:368-372`). A single frozen-Jacobian solve can land on the wrong branch of the `min` → CAS gains latent the soil never released → `whole_water%n_fail`/`we<1e-3` fails.
**Mitigation (mandatory):** stage-internal **active-set** — solve assuming `src_frac=1`; if `uptake_total < coh_transp`, clamp transpiration to supply, recompute leaf-latent/`le_ground`/`coh_qw` with the clamped value, re-solve the block once (bounded 1–2 pass, *not* the old outer Picard). This is why **soil-water supply/evap is IN the coupled block** (§2), correcting the memo.

### R2 (P0-FATAL) — machine-precision conservation is not automatic under a monolithic state solve

Budgets close today because each interface flux is *computed once and applied to both reservoirs* (e.g. `coh_qw` into CAS = water soil released, both × the same `src_frac`). A `(I−γdt·J)⁻¹` *state* solve where each reservoir then reads off its own linearized equation closes budgets only to O(Newton-tol).
**Mitigation:** after the stage state solve, **evaluate each interface flux once from the stage state and drive reservoir increments in flux-divergence form** (leaf↔CAS, CAS↔ground, CAS↔atm, root↔soil-water, soil-top↔CAS). The budget ledger (§6) uses those *same* fluxes. Never "each reservoir solves its own BE and reads a state."

### R3 (P0-FATAL) — soil-top↔CAS / leaf↔CAS must not be explicit

Covered in §2.2. Leaf-T and ground skin are **static-condensed** (Schur/rank-1 using the already-assembled slopes `h_coeff_f,le_slope,lw_slope` `:349-353`, `ggnet·ρ·cp` `:403`) **into** the implicit CAS+soil-top block. The soil-top-energy row (layer 1) **joins the arrowhead**, it is not a decoupled Thomas column. Reject any partition with `src_enth` as a purely explicit frozen source.

### R6 (P1) — hydraulics via plain ESDIRK, NOT ETD `apply_expm`

`apply_expm` is the exact frozen-2×2 exponential — but exponential integrators have their *own* φ-function order conditions; dropping an exact-exp block into an additive-RK tableau does **not** inherit the tableau order (the shared coupling terms `e_transp`, `ψ_soil` are only ESDIRK-order → reintroduces splitting error).
**Mitigation:** advance hydraulics with the **same ESDIRK stage** — another 2×2 diagonal block `(I−γdt·M)⁻¹`, M eigenvalues ≤ 0 ⇒ L-stable. One solver technology, uniform order. This also *simplifies* vs the memo. Keep `apply_expm` only as the RK4-reference hydraulics path. Because `gs`/`ψ_leaf` is lagged one step (`:274`), hydraulics is **one-way downstream** — solved *after* the surface block, not inside the arrowhead.

### R9 (P2) — generalize the CAS BE form to the γ-stage form

The verbatim BE form `enth1=(wcap·enth0+Δt(src+gah·atm))/(wcap+Δt·gah)` (`:410`) is BE = ARK1 only. Copy-pasting pins you to 1st order. **Generalize** to `(wcap·Y_base + γdt(src_impl + gah·atm))/(wcap + γdt·gah)` — atm term + condensed leaf/ground in the implicit tableau, biotic sources in the explicit tableau. The implicit Jacobian diagonal becomes `1 + γdt·wci·gah`; soil tridiagonal diagonals become `c_eff·dz/(γdt)` and `cc·dz/(γdt)+dsk·dz`.

### R7 (P1) — ψ must enter the column error controller

`keff = kirchhoff_edge(plc(ψ_wood))`, `cl,cw = capacitance(ψ)` are nonlinear (`hyd:332-334`). One freeze/stage over 900 s is coarser than the deleted substepper. **Mitigation:** the single column WRMS norm **must include ψ** (per-cohort, own `atol_ψ/rtol_ψ`) so a cavitation front forces whole-step rejection/shrink. Folding the hydraulics substepper is safe *only* if ψ is in the norm.

### R8 (P1) — Richards nonlinearity: frozen-Jacobian vs in-stage Newton

At infiltration into dry soil, K(θ)/C(ψ) span orders across a layer; a single frozen-coefficient BE at 900 s can push θ past saturation (mass conserved by flux form, but *state* wrong). **Mitigation (minimum-safe, do not ship without it):** (a) θ in the column norm with tight `atol_θ` so fronts trigger shrink; (b) retain a **capped modified-Picard (2–3 iters)** *inside* the soil-water implicit stage (reuse `soil_be_single_step`'s `SOIL_LIN_PICARD` loop `:328,331-360`, node-local, cheap). This recovers the old Celia robustness without the deleted *outer* split. Optional opt-in `ccfg%stage_newton` allows 1–2 Newton iters on the qsat latent border term (off by default, off on GPU).

### DAE order honesty (R4/R5/R10)

- **R10:** an IMEX-Euler MVP is **1st order = same order as the split it replaces** — it gains the coupling win but not the accuracy win, so the "≥2nd order" headline would be unverifiable at MVP. **Therefore the MVP is ARK2(1)3L[2]SA, not IMEX-Euler** — minimal extra cost, delivers *both* the coupling win and a demonstrable order bump.
- **R5:** ARK has **stage order 2**; on stiff/algebraic components a general ARK reduces toward stage order. Stiff accuracy (L[2]SA) buys *order-2 robustness on the algebraic variables* — it does **not** recover classical order 3. State the target as **≥2nd order**; calibrate tests to p≈2; require the stage-internal algebraic (leaf/ground/`src_frac`) solve to be *consistent at each stage* (else you drop below 2).
- **R4:** the algebraic leaf constraint is *not exact* — `qsat` linearized around `tcas`, LW around `te`, `gs` frozen (`:337-353`) — one Newton step leaving an O((ΔT_leaf)²) residual that does **not** vanish as Δt→0. **Do not order-test leaf T.** Order-test a *differential* variable (CAS enthalpy or soil-top temperature) via self-convergence against the same code at dt/64 so the identical constraint linearization cancels (§7).

### GPU warp-uniformity / fixed-substep mode (R12 + §7 map)

The adaptive accept/reject branch is per-patch data-dependent → warp divergence. Provide `ccfg%ark_adaptive=.false.` + `ark_fixed_substep = N` → every patch takes exactly `N` equal ARK steps of `dt_fast/N`, no error norm, no accept/reject, no per-lane branch (mirrors the folded-away `SOIL_SUBSTEP_FIXED`/`HYDRO_SUBSTEP_FIXED`). Stage kernels are already branch-free/fixed-size: closed-form 2×2/3×3 border inverse, fixed Thomas sweep, uniform `1..ncoh` hydraulics loop, scalar CO₂. No device `error stop` (convergence failure sets a host-checked flag). Nonlinear treatment on the offload path stays at single frozen-Jacobian (no in-stage Newton). Adaptive = CPU/ifx production path; fixed-substep = nvfortran-GPU path; both share identical `column_derivs` + stage solvers.

---

## 6. Conservation ledger under stage accumulation

Committed state `y^{n+1} = y^n + dt·Σ_i b_i·f(Y_i)`. Because every tendency in `f` is in flux-divergence / capacity form (§2 table), interior fluxes telescope to zero per column and only **boundary faces** survive:
`Δ(store) = dt·Σ_i b_i·(boundary_in_i − boundary_out_i)` — exact to machine precision.

Implementation:
- `column_derivs` returns boundary fluxes in `diag`: `g_top`, `gah·(enth−enth_atm)`, `gaw·(shv−shv_atm)`, `gac·(co2−co2_atm)`, `drainage`, `runoff_surf`, `infiltration`, `coh_rnet`.
- During the stage loop, accumulate each boundary flux with the **same `b_i` weight**: `bflux_accum += b_i·diag_i`.
- At commit, call the **unchanged** `budget_accumulate`/`track_resid` (`:481-507`) with the b-weighted boundary fluxes as `in`/`out`. The seven budgets (`cas_energy/water/co2`, `soil_energy/water`, `whole_water/whole_energy`) close by construction. The stiffly-accurate tableau (last stage = solution) makes endpoint and weighted-sum fluxes coincide at the final stage → tight ledger. **Zero edits to the budget code.**

---

## 7. Verification & test plan

New `test/test_ark_integrator.f90` (ctest `ark_integrator`, links `meds_aux meds_testsupport`, mirroring `test_picard_coupling`'s CMake block `:298-301`), plus targeted edits.

| Test | Requirement | Concrete assertion |
|---|---|---|
| **(a) Regression gate — small-dt agreement** | ARK@refined dt reproduces the split's continuous limit (both are consistent; they *legitimately differ at 900 s* — that is the point, cf. `test_picard_coupling.f90:107-114` asserts `max|dTcas|>0.3`) | Run the 24h diurnal setup from `test_column_dynamics` under ARK at `dt_fast=900/16≈56 s`, sub-sample to the 900 s grid, assert `max|T_ark − T_split_refined|` below budget `atol`; gap shrinks under further refinement (O(dt) agreement, not bit-identical) |
| **(b) Order-of-accuracy** (primary scientific justification; no analogue today) | ARK ≥ 2nd order on a **differential** variable via **self-convergence** (R4/R5) | Integrate a smooth diurnal window at `dt∈{dt0,dt0/2,dt0/4}`, reference = same code at `dt0/64` (Richardson); measure `p=log2(e(dt)/e(dt/2)) ≥ 1.9` on **CAS enthalpy or soil-top T** (NOT leaf T); hold `gs` frozen consistently across refinement |
| **(c) Budget/conservation to machine precision** | all 7 budgets close (structural, order-independent) | Copy `test_column_dynamics.f90:93-99` + `test_fast_loop.f90:82-84`: `whole_energy%n_fail==0`, `whole_water%n_fail==0`, `we<1e-3`, `ww<1e-8`, all `cas_*`/`soil_*` `n_fail==0`. Preserve the `column_budget_t` residual-tracking calls |
| **(d) RK4 independent oracle** | ARK and fully-explicit RK4 agree (shares no code with the implicit assembly → rules out shared-bug false pass) | `INTEG_RK4` dispatch; at small `dt_fast` (RK4 stable within ~19 s cap) assert ARK↔RK4 agree to tight tol on a non-stiff-dominated window |
| **(e) Stiffness robustness** | one ARK step at full `dt_fast` stays bounded (stiffness ratio ~6.5e5) | Run 24h at `dt_fast=900 s` **and** `1800 s`: assert `minval(T)>280 .and. maxval(T)<315` (cf. `:108-110`) + budgets close. Optional non-gating diagnostic: RK4@900 s violates a bound (documents the stability separation; keep non-gating to avoid platform-fragile blow-up asserts) |
| **(f) Folded-controller observability** | the single column controller's adaptivity is regression-testable (replaces per-kernel `nsub` reporting) | Assert `column_fast_step` reports a substep/iteration count via the `iters` out-arg (`:182`) |
| **Golden anchor** | bit-level regression pin survives the overhaul | Keep `test_picard_coupling.f90:86-89` (`tc_split(54)==292.450065`, `ss_split(54)==296.218258`) while split exists in RK4-reference form; once split is removed, **repoint the anchor to the ARK trajectory** at the same `(54)` indices (regenerate once in Debug/ifx, pin) |

**Cross-cutting gates (every new test/module):**
- **Dual-compiler.** Must build+pass under **both** `ifx` (`-stand f18 -check all -fpe0`) **and** `nvfortran` multicore (`-DMEDS_GPU=multicore -Mbounds -Ktrap=fp`). Green ifx alone is **not sufficient** (CLAUDE.md portability trap). Build both `build-ifx` and `build-mc` on the new modules before believing any green result.
- **nvfortran array-temp trap (issue #7) — highest-risk pitfall here.** ARK stage assembly (`k_i` stage vectors, arrowhead RHS, per-cohort ψ stages) is exactly the array-temp-heavy code that triggers it. **Never** `call thomas_solve(a,b,c, assemble_rhs(...), …)` — bind `rhs = assemble_rhs(...)` to a named array first. nvfortran silently miscompiles the temp descriptor at `-O2`; ifx only emits an `arg_temp_created` remark, so the ifx suite hides it.
- **Kernel tests retained.** `column_energy`, `surface_energy`, `column_hydrology`, `column_co2`, `plant_hydraulics` still exercise the *kernels* the ARK couples — unchanged.

---

## 8. Complexity ledger & risks

### Honest ledger

| REMOVED | ~lines | ADDED | machinery |
|---|---|---|---|
| outer Picard/split shell (`MCD:323-458`) | ~135 | ARK tableau + stage storage `k_i` | new |
| `soil_water_advance` substepper | ~65 | **arrowhead/Schur assembly** (off-diag Jacobian, Thomas-eliminate soil columns, dense border solve, back-sub) | **does not exist today** |
| `solve_plant_water` substepper | ~40 | per-cohort leaf static-condensation (n rank-1 updates) | new |
| 7 config fields × 3 sites + `req_scheme` + validation | ~30 | `src_frac` complementarity active-set | new |
| — | — | single-flux-per-interface post-solve commit | new |
| — | — | one column error controller | reuses `adaptive_step_update` |

**Verdict — the "simplify" claim is half true.** Net-negative on *methods* and *config knobs* (real: fewer integrators, no outer iteration, one controller instead of three substeppers, net-negative config line count). Net-**positive** on *numerical-linear-algebra machinery* and a strictly harder conservation-debugging story: today a budget failure localizes to one independent block; under the coupled arrowhead a single wrong Jacobian off-diagonal sign breaks conservation globally. Trading "several simple weakly-coupled BE solvers + an outer iteration" for "one strongly-coupled solver" is the right call for order and Lie-error removal, but **bill it as "one coupled solver + one controller," not "less code."**

### Top risks → mitigations (prioritized)

| # | Risk | Severity | Mitigation | §ref |
|---|---|---|---|---|
| R1 | `src_frac` complementarity → conservation break | P0-FATAL | stage-internal active-set (1–2 pass), soil-water supply IN block | §5/R1, §2 |
| R2 | monolithic state solve ≠ machine-precision budgets | P0-FATAL | single-flux-per-interface commit feeding existing budgets | §5/R2, §6 |
| R3 | soil-top↔CAS / leaf↔CAS left explicit → unstable at 900 s | P0-FATAL | static-condense leaf/ground; soil-top row in arrowhead; reject "or explicit" | §2.2, §5/R3 |
| R7/R8 | folded substepper misses cavitation/wetting fronts | P1 | ψ AND θ in the one column norm; capped inner soil-water Picard (2–3) | §5/R7,R8 |
| R6 | ETD `apply_expm` in additive ARK is order-fragile | P1 | plain ESDIRK 2×2; `apply_expm` reference-only | §5/R6 |
| R4/R5/R10 | order test false-fails; MVP no order gain; over-sold order 3 | P1/P2 | order-test differential var via self-convergence; MVP = ARK2(1); advertise ≥2nd | §5, §7 |
| R12 | nvfortran array-temp miscompile in stage assembly | P1 | named array temps; build both compilers | §7 |

### De-risk-first order (prototype before committing to the arrowhead)

1. Extract `column_derivs` and verify it reproduces the split trajectory when driven by IMEX-Euler (γ=1, 1 stage) at small dt — bit-anchor bridge. **Cheapest proof the pure-RHS refactor is faithful.**
2. Prototype the `src_frac` active-set + single-flux commit on the *existing* split solver first (R1/R2 are the fatal ones and are independent of the ARK tableau) — prove budgets still close to machine precision *before* adding stage machinery.

### Phased roadmap

| Phase | Deliverable | Effort / risk | Acceptance test (green on ifx **and** nvfortran) |
|---|---|---|---|
| **P0** | Extract `pure column_derivs` (Step A/B/C) from pre-pass + `:337-434`; leave the split solver calling it | M / low | split trajectory bit-reproduced; existing `column_dynamics`/`fast_loop`/`picard_coupling` unchanged-green |
| **P1** | `src_frac` active-set + single-flux-per-interface commit on the existing solver | M / **high (R1/R2)** | 7 budgets close to machine precision at 900 s across the diurnal cycle |
| **P2** | Implicit stage solvers (reuse `thomas_solve`, `soil_be_single_step` capped-Picard, ESDIRK 2×2); assemble dense Schur border + soil-top arrowhead row | L / **high (R3, arrowhead)** | coupled surface block closes budgets; matches P1 fluxes |
| **P3** | Wire **ARK2(1)3L[2]SA** + `adaptive_step_update` controller (WRMS incl. θ and ψ); b-weighted ledger; `iters` reporting | M / med | order test `p≥1.9` on CAS enthalpy; stiffness bounded at 900 s & 1800 s; substep count observable |
| **P4** | Config surface swap (`integrator`, `ark_*`); RK4 reference dispatch; repurpose `test_picard_coupling`; new `test_ark_integrator`; retire `SCHEME_*` | M / low | full §7 suite green; RK4-oracle agreement; golden anchor repointed |
| **P5** | GPU fixed-substep lockstep mode (`ark_fixed_substep`, branch-free, no device `error stop`) | M / med | nvfortran `-DMEDS_GPU=multicore` green; warp-uniform stage sequence; no array-temp miscompile |
| **P6 (follow-on)** | Option A: fold soil-water top node + (if ever un-lagged) hydraulics into the arrowhead; promote to ARK3(2)4L[2]SA | L / med | order/stiffness gates still green at higher order |

Each phase is independently landable; the regression gate at every phase is "reproduce current split at small `dt_fast` + all 7 budgets close." P0–P3 already remove the dominant split (Lie-Trotter) error — that is the smallest correct MVP delivering both the coupling and accuracy wins.

---

## 9. Vetting notes & clarifications (post-review Q&A)

These points were verified directly against source during review and clarify three questions that came up; they refine — do not overturn — the design above.

### 9.1 Hydraulics is solved *by* the integrator (it "leaves the block," not the integrator)

"Hydraulics leaves the block" (§2) means it leaves the **dense coupled surface arrowhead**, not the IMEX-ARK. There are three implicit tiers: (a) the coupled surface *block* (CAS enth/shv + soil-top energy, leaf-T & ground condensed in) solved as one dense system; (b) **own separate implicit blocks** — soil-heat column (Thomas), soil-water column (Thomas), **hydraulics 2×2 per cohort (ESDIRK)**, CAS CO₂ (scalar); (c) explicit sources. Hydraulics is tier (b): fully implicit, L-stable, its ~17 s stiffness properly absorbed — just solved as a standalone per-cohort 2×2 sequenced *after* the surface block, not folded into the arrowhead.

Why it can be standalone: the within-step coupling to the surface is **one-way**, because (i) `gs`/`ψ_leaf` is lagged one step (`meds_column_dynamics.f90:274` — stomatal conductance uses start-of-step ψ, frozen for the whole step), so ψ does not feed back into transpiration demand within a step; and (ii) the only back-coupling, the `src_frac` supply limiter (`:385`), is contractive (∈[0,1]). Information flows surface+soil-water → hydraulics; nothing hydraulics produces re-enters the surface energy/vapor balance within the step, so there is no two-way coupling to split and no accuracy lost by sequencing. **Contingency (P6):** if MEDS ever un-lags `gs↔ψ` (within-step stomatal response to leaf water potential — a real drought feedback), hydraulics becomes two-way and must fold into the arrowhead.

### 9.2 Radiation: short-wave is a frozen source; only long-wave contributes an implicit term

The two-stream RT solver (**both** bands) runs once per sub-step, **lagged**, in `apply_rt_forcing` (`meds_fast_loop.f90:449`), *outside* the integrator — it is frozen forcing, not re-solved inside the ARK stages. Inside the coupled energy balance:
- **Short-wave** (`abs_sw`, `abs_par`, `abs_sw_ground`) — pure constant source, zero temperature dependence → **explicit/frozen**, never in the implicit solve.
- **Long-wave** enters as (a) the two-stream **net LW at the base temperature** (`abs_lw`, `abs_lw_ground`) — a frozen constant, like SW — **plus** (b) the **local Stefan–Boltzmann self-emission linearization** `lw_slope = 4·ε·σ·T³·LAI` (`meds_column_dynamics.f90:349`), which sits in the implicit denominator (`:353`) → **this is the only radiative term in the coupled solve.**

So "only LW is coupled" is correct — with the precision that what is implicit is the `4εσT³` self-emission slope, **not** the two-stream LW field (which is lagged along with SW). The ARK plan preserves this: SW + `abs_lw` stay explicit frozen sources in `column_derivs`; the leaf `lw_slope` is carried into the surface border via the static-condensed leaf-T row.

**Gap to tidy (free):** today the implicit `lw_slope` is applied to the **leaf only**. The ground surface balance (`meds_column_energy.f90:245`) uses `rn = abs_sw + abs_lw` fully frozen — no `4εσT_ground³` term — so the ground's own LW-emission feedback is lagged, and only ground sensible + latent couple its temperature implicitly. When the ARK static-condenses the ground skin into the surface block, add a symmetric ground `lw_slope` (numerically ~5 W/m²/K; an accuracy consistency fix, not a stability issue).

### 9.3 ARK tableau attribution — pin the coefficients before coding

The plan names the tableau "ARK2(1)3L[2]SA" with `γ = 1 − 1/√2 ≈ 0.29289322` and attributes it to Kennedy–Carpenter. That γ value is actually the ARS(2,2,2) / Giraldo–Kelly–Constantinescu "ARK2" scheme; Kennedy–Carpenter's smallest additive pair is ARK3(2)4L[2]SA (γ ≈ 0.435866). The **scheme class** specified (3-stage, order-2 with embedded order-1 estimate, single γ, L-stable, stiffly-accurate ESDIRK-IMEX) is correct and appropriate regardless of attribution — but **pin the exact `a^E`, `a^I`, `b`, `b̂`, `c`, `γ` coefficients and the citation when implementing P3** (candidates: ARS(2,2,2); Giraldo et al. 2013 MWR "ARK2"; Kennedy–Carpenter 2003 ARK3(2)4L[2]SA if going to order 3).

### 9.4 Facts verified directly against source during review

Removal target = the single `do iter = 1, niter` loop (`meds_column_dynamics.f90:323-458`), with the `state^n` store snapshot/reset (`:374`) present *only* because `column_hydrology_flux`/`plant_water_flux` mutate state in place. All three CAS twins (enth `:410`, shv `:411`, CO₂ `:462`) are implicit-in-atm BE (identical `(wcap+dt·ga)` denominator) — there is **no** vapor-twin stability gap. Leaf-T is the diagnostic algebraic solve (`:352-354`, index-1 DAE). `src_frac = min(1, uptake/demand)` at `:385`. Config plumbing: `SCHEME_*` at `meds_config.f90:58-59`, field `:72`, picard knobs `:76-80`, validation `:228-230`; `req_scheme` at `meds_config_io.f90:166-180`; `build_fast_context` copy `meds_fast_loop.f90:106`, Picard branch `:427`. Fast→slow seam (four `*_accum` daily carbon totals) is untouched; slow loop is confirmed simple daily Euler/exact.

---

*Design doc for the MEDS fast-loop numerics overhaul. Implementation tracked on branch `feature/imex-ark-numerics`, starting with P0 (extract the pure `column_derivs` RHS + faithfulness test).*