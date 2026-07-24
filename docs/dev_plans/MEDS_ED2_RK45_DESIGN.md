# MEDS ED2-Faithful RK45 + Internal-Water-Mass Design

**Status:** **P0+P1 IMPLEMENTED** (2026-07-23, on branch `feature/ed2-rk45-water-mass` — P0 committed, P1
not yet; not yet PR'd): internal water MASS is the persisted plant-hydraulics state (both split and ARK
paths), the split path's transpiration↔uptake ledger closes to machine precision (§8 gate 2), and the
canopy-surface water store (interception/film-evap/dew, §3.4) is wired on the split path (see "P0
implementation notes" and "P1 implementation notes" at the end of §8). P2 (the RK45 stepper itself) and P3
(biomass seam/clamp hardening) remain design-only. Companion to `MEDS_NUMERICS_SCOPING.md` §12 (the
baseline-schemes round) and §12.6 (the verified ED2 DTLSM frozen/integrated/refreshed catalogue). This doc
specifies the FIRST of the four baseline schemes — a faithful reproduction of ED2's adaptive Cash–Karp
**RK45** — and the state-vector change it is paired with: **leaf/wood internal water mass become
prognostic states so the column water budget closes within each `dt_fast`.**

Every ED2 equation below was read directly from `ED2/ED/src/dynamics/rk4_derivs.f90` and cross-checked; the
MEDS-side facts from a source sweep; and the design was hardened by an adversarial pass that **refuted an
earlier naive version** (see §9, Pitfalls — it is the most important section).

---

## 1. Purpose and the one decision that drives everything

Two coupled deliverables:

1. **ED2-faithful RK45** — the adaptive embedded Cash–Karp 5(4) integrator (ED2's RKQS), over the ED2
   partition (§12.6): physiology/radiation/**hydraulic fluxes** frozen once per `dt_fast`, the physical
   stores integrated, temperatures/aerodynamics refreshed per stage.
2. **Internal water mass as state** — replace MEDS's current *ψ-as-prognostic* (advanced by an exact
   matrix exponential, operator-split out of the tableau) with *water-mass-as-prognostic* (ψ diagnosed
   from mass via the existing PV curve), so soil → wood → leaf → CAS → atmosphere is a chain of mass
   stores, each `d(store)/dt = in_flux − out_flux`, and the column water budget closes by construction.

**The single decision that makes both work: FREEZE the hydraulic fluxes for the RK march**, computed as the
**time-averaged flux from an exact ψ solve over the step** — exactly as ED2 does. ED2 computes
`wflux_wl`/`wflux_gw_layer` **once per DTLSM (~600 s)** in `plant_hydro_driver`: an *analytical
exponential-decay projection* of the end-of-step ψ (`calc_plant_water_flux`, `exp(a·dt)`), giving
`wflux = dW/dt + transp` — the average flux, held constant across every RK stage. MEDS reproduces this with
its existing matrix-exponential ψ solve as the pre-pass, freezing the same average flux (§5). Freezing the
average is what makes the whole thing tractable:

- **Non-stiff for an explicit RK45.** If the fluxes were recomputed from the evolving mass each stage,
  `flux = conductance·Δψ(mass)` reintroduces the ~17 s hydraulic RC eigenvalue into the explicit stability
  region — the exact stiffness MEDS removed by splitting ψ out. With the fluxes frozen, `d(wood_water)/dt`
  is a constant and `d(leaf_water)/dt = frozen_sapflow − transp` is driven only by the CAS-humidity
  timescale (~130 s), so the explicit RK45 stays stable at its CAS-limited step (~360 s → ~3 substeps).
  The stiffness is resolved *inside* the pre-pass (the matrix exponential is exact/unconditionally stable),
  not in the RK45.
- **Single-definition closure.** The *same* frozen uptake value debits the soil and credits the wood; the
  *same* frozen sapflow debits the wood and credits the leaf. Exactly one number per interface, so closure
  is arithmetic.
- **Same method as ED2 — the average flux is not an improvement, it is the reproduction.** ED2 does *not*
  use an instantaneous flux: `calc_plant_water_flux` (`plant_hydro.f90:769-776, 885-892`) projects the
  end-of-step ψ **analytically** (`exp(a·dt)`, the exact RC decay over the ~600 s DTLSM) and forms
  `wflux = (proj_ψ − ψ)·c/dt + transp = dW/dt + transp` — the time-**averaged** flux. MEDS's matrix
  exponential is the same construction. So freezing the exact-solve average IS ED2-faithful; the earlier
  claim that it "improves on ED2's instantaneous flux" was wrong (corrected §5.3).

An earlier design draft recomputed the fluxes from evolving mass every stage. That version is **both stiff
and broken on closure** (two incompatible uptake definitions). Do not build it. §9 records why.

---

## 2. The ED2 equations (verified against source)

All quantities per **m² ground** (ED2 nplant-weights the per-plant fluxes at the interface — this is the
units discipline the port must copy). From `rk4_derivs.f90`:

```
wloss    = wflux_gw_layer(k1,ico) * nplant        ! [kg/m2g/s] FROZEN root uptake, soil layer k1 -> wood   (:879)
wflux_wl = cpatch%wflux_wl(ico)   * nplant        ! [kg/m2g/s] FROZEN sapflow, wood -> leaf                 (:2096)
transp                                            ! [kg/m2g/s] REFRESHED, frozen gsw x evolving CAS gradient

d(wood_water_im2)/dt =   wloss  - wflux_wl        ! uptake in, sapflow out                          (:895, :2099)
d(leaf_water_im2)/dt = wflux_wl - transp          ! sapflow in, transpiration out                          (:2098)
```

Energy travels **with** the water (advective enthalpy — ED2 comment: *"we will do it in two steps so we
ensure energy is conserved"*), at the **upwind** temperature, flux frozen but temperature refreshed:

```
qloss     = wloss    * u_liq(T_up)   T_up = soil-water temp if wloss>=0    else wood_temp   (:881-885)
qwflux_wl = wflux_wl * u_liq(T_up)   T_up = wood_temp        if wflux_wl>=0 else leaf_temp  (:2109-2112)

d(wood_energy)/dt += qloss - qwflux_wl            ! uptake enthalpy in, sapflow enthalpy out   (:897, :2123)
d(leaf_energy)/dt += qwflux_wl                    ! sapflow enthalpy in                              (:2122)
d(soil_energy(k1))/dt -= qloss                    ! soil loses the uptake water's enthalpy
```

The transpiration then carries the leaf water (and its latent heat) to the CAS vapour store. ψ is **not**
an ED2 state; it is diagnosed from `leaf_water_im2`/`wood_water_im2` via the PV curve for the frozen
gs (Category-0) and for diagnostics.

---

## 3. The column water-budget closure (the heart of the request)

The CAS vapour store is fed by **four** parallel pathways, not just transpiration. This plan makes the
**transpiration** pathway close (via internal water mass); **soil evaporation** and **snow sublimation**
already close in the current MEDS ledger; the **canopy-surface (interception/dew)** pathway is a separate
store that is only partly built (§3.3). A *complete* column budget sums all of them — the plan must not
close the transpiration pathway while silently leaving another open.

### 3.1 The full flux graph (per m² ground, one consistent flux per interface)

```
                         ┌─ leaf_int_water ──transp──┐
 soil ──wloss──> wood_int_water ──wflux_wl──┘         │
 soil ─────────────────────── soil_evap ─────────────┤
 snow ─────────────────────── sublimation ───────────┼──> CAS_vapour ──ET──> atmosphere
 leaf/wood_surface_water ──film_evap / −dew ──────────┘
 precip ─> {interception→surface_water, throughfall→soil, snowfall→snow}
```

### 3.2 The transpiration pathway (what THIS plan adds — internal water mass)

```
d(soil_water·ρ)/dt   −= wloss            (root uptake; + infiltration/drainage boundary terms, unchanged)
d(wood_int_water)/dt  =  wloss − wflux_wl
d(leaf_int_water)/dt  =  wflux_wl − transp
d(CAS_vapour)/dt     +=  transp
```

Sub-sum telescopes: `d/dt(soil_uptake_part + wood_int + leaf_int) = −transp`, exactly offset by the CAS
credit — `wloss` and `wflux_wl` cancel between the store they leave and the store they enter, **because
each is one number used on both sides**. This is the closure the current MEDS ARK lacks: today the soil
sink is fixed from the state-ⁿ demand while the stages re-evaluate transpiration and the difference is
dropped (`meds_fast_ark.f90` inflates the `whole_water` tolerance to `max(1e-3, |fro%uptake|·dt_fast)` to
hide it). With internal water mass as a state the plant capacitance **absorbs** the uptake↔transpiration
mismatch explicitly, so that inflation is removed for a real machine-precision check.

**Why frozen-uptake / refreshed-transpiration still closes** (adversarially confirmed): closure needs only
that each interface use *one* flux on *both* sides within the step — not that the fluxes match each other.
`wloss`/`wflux_wl` are frozen constants (trivially identical on both sides); `transp` is refreshed but the
*same* refreshed value debits the leaf store and credits the CAS.

### 3.3 The other three pathways (must be in the SAME ledger)

- **Soil / ground evaporation** — a direct soil→CAS vapour flux, parallel to transpiration, NOT through the
  plant. **Already wired and closing:** `soil_evap = hflux%soil_evap` is the single ground-latent authority
  (`meds_fast_split.f90:431`), enters `src_vap = coh_transp + soil_evap + subl_mass/dt`
  (`:499`), and the `whole_water` ledger already carries the soil-moisture debit against it
  (`w_soil0/1 = Σ θ·dz·ρ + w_surface`, `:618`). The plan **reuses this unchanged**; it must simply be
  shown as a term (`d(soil)/dt −= soil_evap`, `d(CAS_vapour)/dt += soil_evap`) so the combined ledger is
  complete.
- **Snow sublimation + melt** — `subl_mass` sublimes snow→CAS vapour; melt drains snow→soil.
  **Already in the ledger** (`swe0/swe1` summed into `whole_water`, `:618`). Reuse unchanged; show the
  terms.
- **Canopy-surface (interception film + dew/condensation)** — precip intercepted onto leaf/wood surfaces,
  re-evaporating to the CAS, and the reverse (CAS supersaturation condensing as **dew** onto the surface).
  **This plan WIRES it in** (revised per review — it is required for ED2 column-physics parity, not
  deferred). Today `intercept_canopy_layer` exists (`meds_vegetation_biophysics.f90:192`) but is never
  called in the fast loop, there is no persistent surface-water state (`leaf_water` lives only on the
  working `leaf_energy_env_t`), and the §8g condensate **leaks out of the column**
  (`bf%whole_wat_out += sf%cond`, `meds_fast_ark.f90:194`). §3.4 specifies the fix.

### 3.4 The canopy-surface-water store (wired, ED2 `leaf_water` analogue)

Add a persistent **surface** water store per cohort — `leaf_surf_water`, `wood_surf_water` [kg/m² ground]
— DISTINCT from the internal (xylem/symplast) water of §3.2. ED2 carries both as separate RK4 states
(`leaf_water` external vs `leaf_water_im2` internal, §12.6); MEDS must too, because they are different
physics: the surface film evaporates directly through the leaf boundary layer (no stomata), the internal
water feeds transpiration through stomata.

**Water tendencies** (per m² ground):

```
d(leaf_surf_water)/dt =  intercept_leaf − film_evap_leaf   (film_evap<0 ⇒ dew/condensation IN)
d(wood_surf_water)/dt =  intercept_wood − film_evap_wood
d(soil_water)/dt      += throughfall + (precip − Σ intercept)   (canopy drip + uncaught precip)
d(CAS_vapour)/dt      += film_evap_leaf + film_evap_wood        (evap OUT of / dew IN to the CAS)
```

- **Interception** — `intercept_canopy_layer` (already written), swept top→bottom over height-sorted
  cohorts, a capacity-limited Beer-fraction bucket; overflow + uncaught precip is throughfall to the soil.
- **Film evaporation / dew** — a single signed flux from the **wetted fraction** `f_wet =
  min(1, surf_water/capacity)` at the leaf/wood boundary conductance (no stomatal resistance):
  `film_evap = f_wet · gbw · ρ · (q_sat(T_surf) − q_cas)`. Positive → evaporation; negative → **dew**
  (the §8g condensate), which now lands on `surf_water` instead of leaking. This replaces the §8g
  condensation sink and **fixes that leak** — the condensate is a real deposit, budget-tracked.
- **Transpiration competes with film evaporation for the leaf.** The stomatal (transpiration) flux draws
  from the **dry fraction** `(1 − f_wet)`; the wet fraction evaporates at the potential rate. This
  partition is ED2-faithful (wetted-fraction weighting) and must be applied so a fully wet leaf transpires
  ~0 and evaporates its film.

**Energy coupling (for ED2 energy parity):** the surface water is part of the leaf/wood **heat capacity**
(`wmass` on `leaf_energy_env_t` already carries total water for `hcap`), so `d(leaf_energy)/dt` gains the
film's storage term; **film evaporation removes latent heat** from the leaf (→ CAS as vapour enthalpy),
**dew condensation releases latent heat** onto the leaf; and interception/throughfall carry the rain's
liquid enthalpy. All at the (refreshed) surface temperature. Omitting any of these breaks the *energy*
budget even if water closes.

### 3.5 The complete ledger + scope

```
d/dt( Σ soil_water + snow + leaf_int_water + wood_int_water + leaf_surf_water + wood_surf_water + CAS_vapour )
   = precip − ET − drainage − runoff
```

**Scope:** the plan makes MEDS's column water AND energy physics match ED2 — it adds the **internal** water
(§3.2, transpiration pathway) AND wires the **surface** water (§3.4, interception/film/dew); `soil_water`,
`snow`, `CAS_vapour` and the soil-evap/sublimation pathways already close and are reused. The acceptance
gate (§8) is that the `whole_water` and `whole_energy` ledgers sum **every** store above and close to
machine precision, with the §8g condensate leak eliminated (routed to dew).

---

## 4. State-vector change (mass replaces ψ; ψ becomes diagnostic)

**Add** per-cohort `leaf_water_mass`, `wood_water_mass` [**kg/m² ground**] at the three state homes ψ
already lives in, and **remove ψ as a prognostic**:

| home | current | change |
|---|---|---|
| `cohort_block` SoA (`meds_core_state_types`) | `psi(N_HYDRO_NODE, ncoh)` | add `leaf_water_mass(ncoh)`, `wood_water_mass(ncoh)`; update **every** lockstep site ψ touches (alloc/dealloc, grow-copy, permute, `copy_cohort_slot`, fresh-slot reset, **fusion**) |
| `patch_biophys_t` (`meds_biophysics_types`) | `psi(:,:)` | add the two mass arrays beside it |
| `column_state_t` / `column_tend_t` (`meds_fast_types`) | `psi` / `dpsi_dt` | replace with `leaf_water_mass`/`wood_water_mass` and `dleaf_water_dt`/`dwood_water_dt`; update the RK state algebra (`state_init`/`axpy`/`sub`/`err_diff`) |

**Diagnose ψ from mass** — one new elemental function (all constituents already exist in `meds_hydr_lib`):

```fortran
elemental pure function psi_from_water_content(w, pi0, elastic_mod, apoplast_frac, water_sat, biomass) result(psi)
   w_sat = water_sat * biomass                      ! saturated tissue water   [kg/plant]
   rwc   = (w - apoplast_frac*w_sat) / ((1.0 - apoplast_frac)*w_sat)
   psi   = psi_from_rwc(rwc, pi0, elastic_mod)      ! existing exact inverse   (meds_hydr_lib.f90:198)
end function
```

The forward map `water_content(ψ)` (`meds_hydr_lib.f90:229`) and capacitance `dW/dψ` already exist; the
Kirchhoff conductance/vulnerability curves and `rhizosphere_cond` are reused unchanged.

**The gs seam stays FROZEN.** `column_prepass` gathers `leaf_water_mass`ⁿ, diagnoses ψ_leaf **once** per
`dt_fast`, and feeds `leaf_gas_exchange_batch` (which uses ψ_leaf in `beta_nonstomata`). This preserves
ED2's Category-0 frozen gs — do **not** refresh ψ_leaf for gs per stage. Only the transpiration *flux*
and the enthalpy tag refresh (via the evolving CAS gradient and leaf/wood temperature).

---

## 5. What the frozen flux does to the plant-hydraulics solver

This is the largest structural change, so in detail.

### 5.1 What MEDS does TODAY — and that it already computes the average flux

`solve_plant_water` / `solve_plant_water_batch` (`meds_plant_hydraulics.f90:116-392`) is a genuine **ODE
integrator for ψ**: a 2-node RC network `dψ/dt = M·ψ + c` (M, c frozen for the step via `freeze_coeffs`),
advanced **exactly** over `dt` with a matrix exponential `e^{M·dt}` (`advance_exact_linear`, underflow-safe
`sinhc` form), under adaptive step-doubling. Crucially, the fluxes it returns are already the
**time-averaged** fluxes over that exact solve — from the *converged storage change*, not an instantaneous
rate:

```
dW_leaf = water_content(ψ_leaf_end) − water_content(ψ_leaf_start)     ! :203-207
flux%sapflow     = dW_leaf/dt + e_transp                              ! = <sapflow>  over dt   :208
flux%root_uptake = (dW_leaf + dW_wood)/dt + e_transp                  ! = <uptake>   over dt   :209
```

So MEDS **already has proposal 2's ingredient.** The only defect is downstream: today ψ is the persistent
state, mass is a diagnostic, and (in the ARK) `flux%root_uptake` fed to the soil differs from the
transpiration the stages re-evaluate — the transp↔uptake gap.

### 5.2 The design (revised per proposal 2): average-flux pre-pass → freeze → integrate mass

**Keep the matrix exponential as the pre-pass; freeze its average flux; integrate mass in the RK45.** Once
per `dt_fast`, before the RK march:

1. Diagnose ψ_leafⁿ, ψ_woodⁿ from `leaf_water_mass`ⁿ, `wood_water_mass`ⁿ (`psi_from_water_content`, §4).
2. Run the exact ψ solve (existing `solve_plant_water`) over the full `dt_fast` with the state-ⁿ transp
   demand → obtain the **time-averaged** `<sapflow>`, `<uptake>` it already returns (§5.1).
3. **Freeze** `<sapflow>`, `<uptake>(k)`, nplant-weighted, for the whole march.
4. Hand the SAME frozen `<uptake>(k)` to the soil sink and to the wood store; the SAME `<sapflow>` to the
   wood and leaf stores. The RK45 then integrates the two mass states with these frozen fluxes:
   `d(wood_water)/dt = <uptake> − <sapflow>` (constant), `d(leaf_water)/dt = <sapflow> − transp(t)`
   (transp refreshed).

This is a minimal change to the existing code: the pre-pass IS today's `solve_plant_water` call; the new
work is (a) making mass the persistent state and ψ the diagnostic (§4), (b) *freezing* the returned fluxes
and using them as the single per-interface number, and (c) the two mass derivatives in `column_derivs`.
The exponential is **kept**, not retired.

### 5.3 The flux IS the exact-solve average — same as ED2 (corrected)

Both ED2 and MEDS compute the frozen flux as the **exact-solve time average**, not an instantaneous rate.
ED2 projects the end-of-step ψ analytically (`proj_leaf_psi = f(exp(a·dt))`, `plant_hydro.f90:769-776`)
and sets `wflux_wl = (proj_ψ − ψ)·c/dt + transp` = `dW/dt + transp`; MEDS's matrix exponential produces the
identical `flux%sapflow = dW/dt + e_transp`. This is the exponential-integrator effective flux — it
accounts for the flux decaying as ψ relaxes over the ~600 s step, so it transfers the correct *total* mass
and cannot overshoot. So the design choice is **not** average-vs-instantaneous (both are average); it is
**freeze the average for the march (non-stiff, exact closure) vs recompute the flux from evolving mass each
stage (stiff, §9)**. Freezing wins. The stiffness is absorbed inside the exact pre-pass (unconditionally
stable); the RK45 sees a constant flux.

### 5.4 Which transpiration drives the pre-pass flux (the review question)

**First, the mechanics — transp is NOT frozen everywhere; the RK45 refreshes it every stage.** In the code,
`transp = g_tr_f · (q_sat(T_leaf) − q_cas)` (`meds_fast_split.f90:332`, `surface_derivs`): `g_tr_f` (the
series conductance carrying the **frozen** gs) is constant, but the leaf→CAS gradient `q_sat(T_leaf) − q_cas`
is made of **state variables** (`T_leaf` from leaf energy, `q_cas` = CAS humidity) that evolve within the
step. So the RK stage RHS re-evaluates `transp(stage)` at each stage, and `d(leaf_water)/dt =
<sapflow>_frozen − transp(stage)`, `d(CAS_vapour)/dt += transp(stage)` both use that refreshed value. **The
only place a single "state-ⁿ transp" appears is the pre-pass ψ projection that sets the frozen `<sapflow>`/
`<uptake>` MAGNITUDE.**

Why the asymmetry (freeze sapflow, refresh transp): transp's driver (`T_leaf`, `q_cas`) has **no fast
feedback** — gs is frozen — so it moves at the ~130 s CAS timescale and refreshing it per explicit stage is
non-stiff. sapflow's driver `ψ_leaf = f(leaf_water)` carries the **~17 s RC feedback** (drain → ψ drops →
flux changes), so refreshing *it* per explicit stage puts that eigenvalue in the RK45 stability region and
blows up. Refresh the feedback-free flux; freeze the one with the fast feedback.

Consequently the transp-BC choice below affects **only the frozen `<sapflow>`/`<uptake>` magnitude** (the
soil-drawdown), never the leaf-water trajectory or the CAS (those use the refreshed transp), and never
closure. The options for that single projection transp:

  | transp for the projection | lag | notes |
  |---|---|---|
  | **last-step average** (ED2 `psi_open`/`psi_closed`, `:207`) | ~1 `dt_fast` (600–900 s) | a genuine average, but lagged a full step |
  | **state-ⁿ** (MEDS today, `transp_pp`) | ~½ `dt_fast` | fresher — a start-of-step snapshot; **recommended MVP** |
  | **step-midpoint** (predicted) | ~0 | least biased estimate of the step-average; needs a cheap state predictor |
  | **self-consistent** (Picard on `<transp>`) | 0 | exact but ~2–3× cost; over-engineering here |

  **Recommendation: state-ⁿ (MEDS's current choice) for the MVP — it is strictly fresher than ED2's
  last-step lag at no cost**, and the lag matters more at 600–900 s than it did historically. Crucially,
  **this is an ACCURACY knob, not a closure knob:** the frozen `<sapflow>`/`<uptake>` are single numbers
  used on both sides of each interface regardless of which transp produced them, and the leaf/wood storage
  absorbs any mismatch between the projection's transp and the RK45's realized transp. So a wrong transp
  here biases the *soil↔plant partitioning and the ψ_leaf lag*, never the total water. The natural upgrade
  is the **step-midpoint** transp, obtainable from the very same midpoint-predictor the scoping-doc §8f
  wants for the coefficient freeze — so improving this and improving the freeze are one piece of
  infrastructure, not two. Defer the Picard-consistent option until a rapid-change window (sunrise,
  cloud-edge, dry-down) shows the state-ⁿ lag actually matters.

  **Resolving it WITH the RK45 (predictor-corrector).** The RK45 hands back the exact step-average
  `<transp>_RK = (1/dt)∫transp(stage) dt` for free — it is the accumulated CAS-vapour credit. So the clean
  one-pass upgrade is: (1) freeze `<sapflow>`/`<uptake>` from transp-ⁿ and run the RK45; (2) read off
  `<transp>_RK`; (3) re-freeze the hydraulic flux with `<transp>_RK` (optionally re-run the march — full
  Picard). This replaces the state-ⁿ lag with the exact realized step-average at the cost of one extra
  (cheap) exponential solve, using the RK45's own output rather than a separate predictor. Recommended only
  if the measurement in the previous paragraph shows the lag matters.

### 5.5 The remaining approximations (bounded)

- **Frozen-coefficient linearization** (M frozen at ψⁿ conductances/PLC) — the same approximation ED2 and
  MEDS already make; the exponential is exact *given* M.
- **2-node topology only** (leaf+wood; 3-node error-stops) — inherited; fine for the MVP.

### 5.6 What changes in the code

| kept | changed / added |
|---|---|
| `solve_plant_water` matrix-exp + step-doubling (now the **averaging pre-pass**) | its returned fluxes are **frozen** and used as the single per-interface number (not re-diagnosed downstream) |
| the constitutive layer (`edge_cond`, `rhizosphere_cond`, `water_content`, capacitance, `psi_from_rwc`) | + `psi_from_water_content` (one elemental fn, §4) |
| — | `leaf_water_mass`/`wood_water_mass` persistent states (3 homes, §4); ψ becomes diagnostic |
| — | the two mass derivatives in `column_derivs`; `GRP_LEAF_W`/`GRP_WOOD_W` tol-groups (§6) |
| — | the ARK `plant_water_tendency` stage RHS is dropped (ψ no longer integrated in the tableau) |

Because the exponential is retained, there is **no ψ-accuracy loss** relative to today — the earlier draft's
"within-step ψ accuracy for closure" trade is **gone**; we get the accurate average flux AND exact closure.
The only cost is the mass bookkeeping. The exponential stays the hydraulics engine for **all** paths
(split/ARK/RK45), so no config split of hydraulics representation is needed.

---

## 6. The RK45 stepper

- **Tableau:** Cash–Karp embedded 5(4), 6 stages, fully explicit (`k_i = f(y + Σ_{j<i} a_ij k_j)`). Two
  b-vectors give the 5th- and 4th-order solutions; their difference is the embedded error.
- **Adaptive (RKQS):** substep within `dt_fast` for stability + accuracy; the embedded error drives step
  size through the **existing** `meds_fast_control` controller (`state_wrms_grouped` + `step_control_factor`
  + the warm-start `adapt_dt_last` already landed). No new controller.
- **State + error norm:** integrate CAS enthalpy/shv/co2, soil energy, soil water, leaf/wood energy, and
  the new leaf/wood water mass. The WRMS gets two new groups — `GRP_LEAF_W`, `GRP_WOOD_W` (atol ~ a small
  fraction of saturated tissue water). ψ is dropped from the norm (it is no longer a state).
- **Reuse:** `meds_fast_rk4_oracle` already does a full-column RK4 over `column_derivs` — lift its stage
  structure. New: the Cash–Karp coefficients, the embedded 4th-order b-vector, the RKQS accept/reject
  driver (or reuse `adaptive_ark_march`'s control loop with an RK stage in place of the ARK stage), and
  the internal-water derivative + enthalpy terms in `column_derivs`.
- **Home:** a new `meds_fast_rk45` module (peer of `meds_fast_ark`), dispatched by a new `INTEG_RK4`
  selector in `meds_fast_split`'s two-way gate.

**Stability estimate (the explicit-scheme risk, resolved by §1):** with fluxes frozen, the stiffest
integrated mode is CAS (τ ≈ 130 s) → dt ≲ 360 s → ~3 substeps at `dt_fast = 900 s`. The internal-water
mass does **not** add a stiff mode because its inflow is frozen (constant) and its outflow (transp) moves
at the CAS timescale. **This holds only while the fluxes are frozen** — see §9.

---

## 7. The two divergences from current MEDS

1. **Aerodynamics (decide before coding).** ED2 REFRESHES `canopy_turbulence8` every RK stage (§12.6);
   MEDS FREEZES aero in `column_prepass`. A *faithful* RK45 calls `canopy_aerodynamics` inside
   `column_derivs` (per stage). Cost: one aero solve per stage (~3×). Make it a config toggle
   (`[fast].rk45_refresh_aero`, default faithful=true) so frozen-aero is measurable against it. Recommend
   measuring the accuracy delta before committing — it may not be worth the cost.
2. **Hydraulics representation.** ED2 has no ψ state; MEDS had ψ as the state. This design moves MEDS to a
   mass state with ψ diagnostic and the flux frozen — but frozen as the **exact-solve average** (proposal
   2, §5.2-5.3), keeping MEDS's matrix exponential as the averaging pre-pass. So it is ED2's *structure*
   with a *better-than-ED2* flux value. This applies to ALL fast paths (split/ARK/RK45), not just the
   RK45, since it is the closure mechanism — so split/ARK are no longer byte-identical (they gain closure).

---

## 8. Config, verification, phasing

**Config (`[fast]`):** `time_integrator = "rk45"` (new `INTEG_RK4`); `rk45_refresh_aero` (default true);
the RK45 reuses `error_level`/`step_controller`/tolerance knobs. All default-off for existing configs.

**Acceptance gates** (per §12.3 — NOT production RMSE at 900 s, which is freeze-limited):
1. **Order-of-accuracy** on a frozen-forcing / manufactured-solution test — the RK45 must show ~5th order
   (embedded 4th) where the coefficient freeze is zero by construction.
2. **Machine-precision column water closure** — the `whole_water` ledger must sum **every store present**
   (soil + snow + wood_int + leaf_int + CAS vapour, **plus canopy-surface film if/once wired**, §3.4) vs
   the boundary terms (precip − ET − drainage − runoff), and close to ~1e-13 with the inflated ARK
   tolerance removed. It must include the soil-evap and sublimation pathways (already present) as terms,
   and must NOT be declared passing while the §8g condensate still leaks. This is the headline deliverable
   — gate on it hard.
3. **Machine-precision column energy closure** — with the `qloss`/`qwflux_wl` enthalpy terms carried.
4. **Cost/stability** via the work counters (`MEDS_NUMERICS_SCOPING.md` §5.3 — substeps, rejections;
   already emitted) — bounded, low rejection.
5. ifx 36/36 + nvfortran multicore; existing paths (split/ARK) NOT byte-identical (the hydraulics
   representation changes for all paths — but they gain closure). Regression vs a fine-`dt` reference.

**Phasing:**
- **P0 — internal water mass + averaged-frozen-flux closure (transpiration pathway).** Add the internal
  mass fields + `psi_from_water_content`; keep `solve_plant_water` as the averaging pre-pass; freeze its
  returned `<sapflow>`/`<uptake>` and use the same value on both sides; prove the transpiration-pathway
  ledger on the split path first (simplest). Delivers the transp↔uptake closure independent of RK45.
- **P1 — canopy-surface water (interception / film-evap / dew), §3.4.** Wire `intercept_canopy_layer` as a
  persistent state; add film evaporation + dew (route the §8g condensate onto it, fixing that leak); the
  wet/dry-fraction transpiration partition; the surface-water enthalpy + latent coupling. Extends the
  ledger to the full store set — this is what makes the column water AND energy physics match ED2.
- **P2 — the Cash–Karp RK45** over the P0+P1 state, with the enthalpy coupling and the aero toggle.
- **P3 — the biomass fast/slow seam** (§9) and cavitation/saturation clamp bookkeeping hardened.

### P0 implementation notes (2026-07-23)

State-vector change landed on **both** paths, not just split: `patch_biophys_t`/`cohort_block` carry
`leaf_water_mass`/`wood_water_mass` [kg/plant] (a `<=0` sentinel means "not yet lazily seeded"; the true
seed `water_content(PSI_INIT,...)` is written on first gather in `meds_fast_dynamics.f90`, since the core/
biophysics layers link `shared` only and cannot reach the plant PFT trait types needed to compute it
themselves). `meds_hydr_lib` gained `psi_from_water_content` — the exact inverse of `water_content`,
composing the existing `psi_from_rwc`. Cohort fusion now nplant-weights the mass fields (extensive, like
AGB) instead of leaf-area-weighting them (§9's pitfall, applied).

**ARK path:** kept the tableau's own `column_state_t%psi` machinery completely unchanged (still a psi
state, still advanced by the operator-split exact matrix-exponential exactly as before) and added only a
thin diagnose-in/convert-out boundary — `psi_from_water_content` when packing `y%psi` in
`build_column_frozen`, `water_content` when committing `y_out%psi` back to `bio%*_water_mass` in
`column_fast_step_ark`. `diagnose(convert(x)) == x` for exact inverses, so this unifies the persisted
representation with minimal risk to ARK's own (still out-of-scope) water-closure behavior.

**Split path — the actual closure fix — was reordered**, not just re-typed: plant hydraulics now runs
**before** the soil Richards solve each Picard pass (previously: soil solve, then hydraulics fed the
soil's post-solve `psi_soil`/theta). Boundary psi (`psi_soil_from_theta`, new elemental use of the
existing soil retention curve) and rhizosphere conductance are diagnosed from **state^n theta** (pre-solve),
matching the frozen-BC spirit of §5. Hydraulics runs at the **full, unthrottled** transpiration demand —
the old instantaneous `src_frac` clamp (`min(1, uptake_total/coh_transp)` scaling transpiration itself) is
gone; the plant's own leaf/wood storage now absorbs any step-to-step soil-supply/demand mismatch instead.
The solver's own aggregate `root_uptake` (`(dw_leaf+dw_wood)/dt + transp`, i.e. the demand INCLUDING any
storage-refill term) becomes the soil's root-sink forcing; if the soil's own `fwilt` limiting can't honour
it in full, a rescale `scale = uptake_total/total_uptake_b` is applied **only** to the wood-side mass
credit (sapflow — the internal wood→leaf transfer — is untouched), so the whole-column ledger closes to
the soil's *true* realized supply. The mass update itself is the explicit Euler form
(`mass_new = mass_old + dt·(flux_in − flux_out)`) using the solver's own returned `sapflow`/`root_uptake`,
not a `water_content(psi_out)` conversion — algebraically identical when `scale==1` (the common case; this
is exactly the identity `solve_plant_water` uses internally to derive those fluxes from its own `dw`), and
the form that stays closure-consistent when `scale<1`. `whole_water`'s ledger gained the plant's own
`nplant·(leaf_water_mass+wood_water_mass)` storage term (start- and end-of-step) — without it, the new
store's real change would read as a leak.

**A real bug the verification pass caught (not anticipated in §9):** `solve_plant_water`'s aggregate
`root_uptake` carries no floor, and — unlike the per-layer multilayer distribution, which the codebase
already floors to `max(...,0)` because hydraulic redistribution (HR) is intentionally not enabled anywhere
in this model — nothing previously fed that aggregate to the soil directly, so the missing floor was never
exercised. Once it became a live soil-forcing input, a freshly-seeded cohort (`PSI_INIT`, a "well-hydrated"
constant) sitting in soil drier than that constant implies produces a genuinely *negative* aggregate uptake
(the solver's own psi trajectory implying net efflux while it re-equilibrates) — which broke
`column_hydrology_flux`'s own internal mass closure when fed straight through. Fixed by flooring
`root_uptake_b`/`root_uptake_layer_b` to `max(...,0)` right after the batch call, extending the project's
existing HR-disabled convention to this new aggregate pathway (both ifx and nvfortran runs below reflect
the fix).

`multilayer_roots`' existing `root_sink_share` one-step-lag mechanism (distributing the requested total
across layers by last step's realized shares) was **left as-is**, just fed the new `total_uptake_b` instead
of raw transpiration demand as the thing it distributes — the reorder makes the lag unnecessary in
principle (this step's own per-layer breakdown is available before the soil call now), but removing it is
an orthogonal simplification to a feature that is off by default, out of scope for this pass.

**Known deferred imprecision (not a P0 regression, and explicitly not gated on):** `eforc%root_heat_sink`
(the soil energy sink for extracted water) is still keyed to full leaf transpiration at leaf temperature,
not the soil's actual realized uptake — getting that right needs the upwind-temperature `qloss`/
`qwflux_wl` advective treatment (§2), deferred with the rest of the energy coupling. `whole_energy` closure
is out of this pass's scope by design (only `whole_water` is gated).

**Verification:** ifx Debug 36/36 + nvfortran multicore 36/36. `whole_water` worst residual ~1e-13 kg/m²
(machine precision) on both the default-coupling and inter-layer-advection runs in `test_column_dynamics`,
and ~1.6e-13 in `test_fast_loop` — gate 2 (§8) met on the split path. Five tests needed updates for the
mass-not-psi representation (`test_fast_loop`, `test_fusion_cohort`, `test_column_ark`,
`test_column_dynamics`, `test_picard_coupling`); one hardcoded golden-value regression anchor in
`test_picard_coupling` (noon CAS/soil-surface temperature under the split scheme) was re-pinned to the new
values — an explicit, documented departure from byte-identical (the user-authorized scope of this pass),
not a silent behavior change.

### P1 implementation notes (2026-07-23)

Canopy-surface water landed on the **split path only** (`column_fast_step` in `meds_fast_split.f90` is the
unified dispatch entry point for both integrators, so its own top-of-routine `error stop` guard intercepts
`canopy_water_on .and. time_integrator==INTEG_ARK` before any dispatch happens — same deferred-scope
pattern already used there for `wood_energy_model`/`leaf_energy_model` PROGNOSTIC-under-ARK). New opt-in `[fast].canopy_water_on` TOML flag (default `.false.`,
threaded `meds_config_t` → `column_config_t` alongside `snow_on`) gates the whole feature; off keeps every
path byte-identical.

**State:** `leaf_surf_water`/`wood_surf_water` [kg/m² **ground**] added to `cohort_block` and
`patch_biophys_t`, seeded to a true `0.0` at every lockstep site (no lazy-init needed — a new cohort's
canopy starts bone dry). This is the OPPOSITE fusion convention from §3.2's internal water mass: the
surface store is already ground-area-referenced, so cohort fusion **sums** it (not nplant-weighted) —
weighting it would double-count the area normalization already baked into each cohort's contribution.

**Mechanism:** one combined `intercept_canopy_layer` bucket per cohort (its own Beer-law screen + capacity
already combine `lai+wai`), swept top→bottom over the height-descending gather order, called with
`e_canopy=0` (capture/capacity only — this call never evaporates); the combined result is split back into
the two persisted stores by area-index share (`lai/pai`, `wai/pai`) since leaf and wood are separate
energy-balance stores. The ACTUAL film evaporation is a separate signed flux computed from the
Picard-converged leaf/wood energy balance: `veg_energy_diagnostic` (`meds_vegetation_biophysics.f90`)
gained four trailing OPTIONAL args (`f_wet`, `le_slope_wet`, `le_ref_wet`, `film_evap`) that fold the
wetted-fraction partition into the SAME linearized solve as the dry/stomatal pathway, rather than bolting
evaporation on afterward as an independent estimate — proven algebraically (and unit-tested in
`test_surface_energy.f90`) to reduce EXACTLY to the pre-P1 behavior when the new args are absent, and to
the fully-wet case's `le_slope`/`le_ref` swap when `f_wet=1`. `leaf_film_coeff` (new, alongside the
existing `leaf_transp_coeff`) gives the boundary-layer-only (no stomatal resistance) conductance for the
wet pathway. The wetted fraction itself (`f_wet_c`) is frozen once per `dt_fast` from the interception
sweep (not re-solved per Picard pass), mirroring how `gs`/hydraulics are already frozen. Film evaporation
is clamped to the water on hand (`min(film_evap, store/dt_fast)`) at the point it is computed — the CAS
credit and the post-Picard-loop store debit therefore use the identical number by construction, the same
discipline as P0's `sapflow_b`/`root_uptake_b` last-pass-commit. Dew (`film_evap<0`) is deliberately left
**unclamped** at that site (see below).

**A cascading set of three real bugs**, each caught by the whole-column budget checks and fixed in turn:

1. **Missing CAS vapour-mass credit.** Film evaporation's *energy* was folded into `coh_qw` but its *mass*
   was never added to `src_vap` — a new `coh_film_evap` accumulator (mirroring `coh_transp`) closes this;
   `src_vap` now reads `coh_transp + coh_film_evap + soil_evap + subl_mass/dt_fast`.
2. **Capacity-overflow silently discarded.** After fix 1, `whole_water` still failed on days with strong
   dew: `intercept_canopy_layer`'s own `w_max` ceiling was being hit and (until this point) any store
   surplus above it just vanished with no offsetting boundary term. The first attempted fix — clamping DEW
   itself so the store could never exceed capacity — is the third bug below. The surviving fix tracks the
   post-commit capacity overflow explicitly (`surf_overflow`) and books it out through `w_out`/`e_out` at
   the same `rain_temp` reference as the rest of the interception ledger, rather than clamping the flux.
3. **Clamp-induced energy-kernel desync.** Clamping dew at the `veg_energy_diagnostic` call site (to keep
   the store within capacity) breaks a subtler invariant: `dh`/`drnet`/`tl` are already solved *assuming
   the full unclamped wet-pathway conductance was active*, so retroactively clamping the credited
   `film_evap` strands the leaf's entire latent capacity into sensible heat — a multi-hundred-W/m² error at
   `f_wet=1`, far worse than the capacity overflow it was meant to prevent. Resolution: dew is never clamped
   at the kernel; only the post-hoc store-capacity overflow (bug 2's fix) bookkeeps the excess. This is the
   general lesson of §9 applied one level deeper — the closure-safe clamp point is a store's *boundary*,
   not a flux feeding an already-self-consistent solve.

**Known deferred imprecision (mirrors the P0 `root_heat_sink` note):** the surface water's enthalpy
(`surf_enth0`/`surf_enth1`) is booked at a fixed `rain_temp` reference rather than a real prognostic
surface-water temperature — unlike transpiration, which closes exactly because the SAME `tl` appears on
both the store and CAS sides of its interface, the surface-water energy ledger uses `rain_temp` on the
store side and `tl` on the CAS-credit side, so a `film_evap · cp_liq · (tl − rain_temp)` mismatch survives,
amplified by `internal_energy_liquid`'s deliberately low reference temperature (`tsupercool_liq≈55.75 K`
makes `u_liq` values ~10⁶ J/kg at normal temperatures, so a small mismatch reads as a large-looking but
still bounded residual). Closing this exactly needs real thermal-inertia infrastructure for the surface
film — genuinely P2/P3 scope per this doc's own "`d(leaf_energy)/dt` gains the film's storage term"
language (§3.4) — so `whole_energy` is asserted against an explicit bound (`< 5e6` J/m²) rather than
machine precision when `canopy_water_on` is on; `whole_water` still closes to machine precision
unconditionally (the bound is energy-only).

**nvfortran gotcha (extends the CLAUDE.md portability trap, a new flavor):** the whole-program optimizer's
sensitivity is not limited to array-valued function results fed straight into a call — a block of
**unconditionally-executed arithmetic that algebraically evaluates to exactly zero** (the surface-water
ledger terms computed once per sub-step even when `canopy_water_on` is off, since every input is seeded
0.0 and never touched) was enough to perturb code generation for unrelated, pre-existing computation
elsewhere in the same routine: `test_fast_loop` failed on nvfortran multicore only (ifx Debug stayed clean
36/36 throughout), on the zero-cohort bare patch specifically, with a large `whole_water` residual that
bore no algebraic relationship to any new term (all confirmed exactly 0.0 via targeted, then-removed debug
prints) — and, tellingly, the failure disappeared whenever an unrelated diagnostic `write` was added
nearby (the classic signature of an optimizer miscompilation, not a logic bug). Binding the `sum()`
array-section argument to a named temporary first (the documented fix for the *other* flavor of this trap)
did **not** resolve it. What did: gating the entire surf_water0/1·surf_enth0/1·intercepted_total block
behind `if (ccfg%canopy_water_on)` — i.e., skipping the dead arithmetic entirely rather than computing then
discarding it — which is both a legitimate performance win (avoids needless per-step FLOPs on the
byte-identical default path, consistent with how every other opt-in fast-loop feature in this codebase
already behaves) and, empirically, the fix. Lesson for future opt-in fast-loop features: gate a new
feature's ledger arithmetic behind its own flag from the start, rather than computing it unconditionally
on the theory that it telescopes to a harmless no-op — on nvfortran, "harmless no-op" and "invisible to the
optimizer" are not the same thing.

**Verification:** ifx Debug 36/36 + nvfortran multicore 36/36. With `canopy_water_on` on
(`test_column_dynamics` RUN 6): peak canopy film water reaches its `dewmx·(lai+wai)` capacity (~0.35 kg/m²
for that test's cohort, confirming the capacity clamp engages), `whole_water` closes with `n_fail==0`, and
`whole_energy`'s worst residual is ~1.7e6 J/m² — within the documented bound above. With the flag off
(`test_fast_loop`, RUN 1–5 of `test_column_dynamics`), residuals are unchanged from the P0 baseline
(`whole_water` ~1e-13–2e-13 kg/m², `whole_energy` ~4e-7–6e-7 J/m²), confirming the feature is a true no-op
on the default path.

---

## 9. Pitfalls (the adversarial pass — read this before coding)

An earlier design draft was refuted here; each item is a real trap with its resolution.

- **FATAL if the flux is refreshed from evolving mass in the RK stages.** Recomputing `<sapflow>`/`<uptake>`
  from the evolving mass-diagnosed ψ each stage reintroduces the ~17 s hydraulic RC eigenvalue into the
  explicit RK45 stability region → tiny substeps or blow-up. **Resolution: freeze the flux for the march**
  (§1, §5) — computed as the exact-solve *average* (§5.3) so it is accurate and non-overshooting, then held
  constant across the stages. The stiffness is resolved inside the exact pre-pass, never in the RK45. This
  is the load-bearing decision of the whole design.
- **FATAL if uptake has two definitions.** If the soil sink uses a scratch-hydrology uptake while the
  plant store uses `rhiz·(ψ_soil − ψ_wood)`, the interface has two numbers and closure is false.
  **Resolution: compute ONE frozen `wloss` in the pre-pass and use it on both the soil and the wood side.**
  This also makes the current soil-before-plant call ordering irrelevant.
- **Units.** The plant store and fluxes are **per plant** (`water_content` → kg/plant, `sapflow`/
  `root_uptake` → kg/plant/s, `rhizosphere_cond` divides by nplant); the soil and CAS are **per m² ground**.
  ED2 multiplies by `nplant` at the interface (`wflux_wl·nplant`). **Copy that exactly** — an unconverted
  per-plant flux crossing into a per-ground store is a silent leak.
- **Biomass fast/slow seam (new leak the ψ-state was immune to).** `W = capacity(biomass)·rwc(ψ)`; biomass
  changes on the slow step, so a persistent conserved `leaf_water_mass` carried across a biomass update has
  the wrong capacity. **Resolution: carry ψ continuous across the slow seam** — diagnose ψ from mass before
  the slow step, re-derive `mass = capacity(new biomass)·rwc(ψ)` after. ψ is the seam-continuous quantity;
  mass is the within-fast-loop prognostic. (Within a fast day biomass is constant, so no in-loop leak.)
- **Cohort fusion semantics change.** ψ is intensive → fusion **leaf-area-weights** it
  (`meds_core_cohort_fusefiss.f90:184`). Water mass is **extensive per plant** → fusion must **nplant-weight
  and conserve the total** (like AGB), NOT leaf-area-weight. Getting this wrong silently violates column
  water conservation across every fuse.
- **Enthalpy coupling must be carried.** MEDS now has separate prognostic leaf and wood energy; the sapflow
  and uptake advect enthalpy between them (`qwflux_wl`, `qloss`, upwind temperature). Omitting it breaks the
  *energy* budget even if water closes. Frozen flux × refreshed temperature (§2).
- **Closure-killers at the stores' bounds.** Internal water hitting 0 (cavitation) or saturation, and the
  soil sink floored at `max(0,…)`: each must be **bookkept** (route the clamped excess to a tracked term),
  not silently applied. The soil solver already bookkeeps its below-wilting give-back into `uptake_total`
  (`meds_soil_water.f90:134-148`) — mirror that discipline on the plant stores.

---

## 10. What is verified vs assumed

- **Verified by direct read:** the ED2 derivative equations and enthalpy coupling (§2, `rk4_derivs.f90`
  :879-897, :2096-2123); that **ED2's `calc_plant_water_flux` is an analytical exponential-decay
  projection giving the time-AVERAGED flux** `dW/dt + transp` (`plant_hydro.f90:769-776, 885-892`), NOT an
  instantaneous rate — using **last-step average transp** (`:34, :207`); that MEDS's `solve_plant_water`
  produces the identical average flux (`meds_plant_hydraulics.f90:208-209`); the nplant-weighting; that
  MEDS has every constitutive piece and needs only `psi_from_water_content`; the three state homes and
  every lockstep touch point; that `intercept_canopy_layer` exists but is unwired and the §8g condensate
  leaks; the fusion weighting difference.
- **Adversarially confirmed:** the closure telescoping; that frozen-`<uptake>`/refreshed-transp does not
  leak; that a per-stage refreshed flux is fatal on stability.
- **Design choice (transp for the projection):** state-ⁿ (MEDS today) over ED2's last-step lag — an
  accuracy knob, not a closure knob (§5.4); step-midpoint via the §8f predictor is the upgrade path.
  **To measure** on a rapid-change window (sunrise / cloud-edge / dry-down) whether the state-ⁿ lag matters
  before adding the Picard-consistent option.
- **Assumed / to measure:** the ~3-substep stability estimate (wet + dry forced windows); the aero
  refresh-vs-freeze accuracy delta (§7.1); the RK45's 5th-order gate (frozen-forcing test); the
  wet/dry-fraction transpiration partition against ED2. The RK45 stepper module and the surface-water
  wiring have not been prototyped.

---

## 11. The fast loop after this development (summary)

### 11.1 End-state architecture — one `dt_fast`

After P0–P3 land, every fast integrator (split, ARK, RK45) evaluates the **same** pure RHS
`meds_fast_time_derivs%column_derivs` over the same state vector, and the whole column water AND energy
budgets close to machine precision. One `dt_fast` runs in three acts:

**Act 1 — the frozen pre-pass (once, at tⁿ, `column_prepass` + the hydraulic averaging solve).**
- Diagnose ψ_leafⁿ, ψ_woodⁿ from the prognostic `leaf_water_mass`ⁿ / `wood_water_mass`ⁿ
  (`psi_from_water_content`).
- Freeze the **Category-0** physiology: gs / GPP / Rd / leaf+stem+root respiration / radiation absorption,
  with gs reading ψ_leafⁿ (ED2-faithful — §12.6).
- Aerodynamics: frozen here by default, or refreshed per stage inside `column_derivs` when
  `[fast].rk45_refresh_aero` (the ED2-faithful choice).
- Run the exact matrix-exponential ψ solve (`solve_plant_water`) once with the state-ⁿ transp demand, and
  **freeze its time-averaged `<sapflow>`, `<uptake>(k)`** (nplant-weighted). One number per interface,
  handed to both the soil sink and the plant stores.

**Act 2 — the RK45 march (Cash–Karp 5(4), adaptive, ~3 substeps for CAS stability).** The state vector and
its treatment:

| state (per m² ground) | in the march | driver |
|---|---|---|
| CAS enthalpy / shv / CO₂ | integrated (in WRMS) | atmospheric exchange + surface sources |
| soil energy (per layer) | integrated (in WRMS) | BE-diffusion terms in `column_derivs` |
| soil water θ (per layer) | integrated (in WRMS) | Richards `soil_water_time_deriv` |
| **leaf / wood internal water mass** | integrated (in WRMS, new `GRP_LEAF_W`/`GRP_WOOD_W`) | `d/dt = <sapflow> − transp(stage)`, `<uptake> − <sapflow>` |
| **leaf / wood surface (film) water** | integrated | interception in, film-evap/dew out |
| leaf / wood energy | integrated (prognostic modes) | radiation + sensible + **advective enthalpy** (`qwflux_wl`, `qloss`) |
| snow (swe + energy) | operator-split (existing) | accumulation / sublimation / melt |
| ψ_leaf, ψ_wood | **diagnostic** (from mass) | — |
| gs, GPP, Rd, `<sapflow>`, `<uptake>` | **frozen** (Act 1) | — |
| transp, temperatures, (aero if faithful) | **refreshed each stage** | evolving CAS/leaf state |

The RK45 reuses the existing `meds_fast_control` controller (WRMS + PI/I + warm start) unchanged.

**Act 3 — commit + close.** The `whole_water` ledger sums **every** store — soil + snow + leaf/wood
internal water + leaf/wood surface water + CAS vapour — against the true boundary (precip − ET − drainage
− runoff), and closes to ~1e-13; `whole_energy` likewise with the advective-enthalpy terms carried.
Condensation is a **dew deposit** on the surface store, not a boundary leak (the §8g fix), so the inflated
ARK water tolerance is gone. The transp↔uptake gap is gone because `<uptake>` is one frozen number on both
sides.

**The result:** MEDS's column water + energy physics matches ED2 (internal + surface water, advective
enthalpy, all four CAS-vapour pathways), on a validated ED2-faithful adaptive RK45 — with a
machine-precision closed budget the current split/ARK lack.

### 11.2 Targeted file changes

**P0 status (2026-07-23):** the `meds_core_state_types.f90`/`meds_core_cohort_fusefiss.f90`/
`meds_biophysics_types.f90` (mass fields only, not the P1 surface store)/`meds_hydr_lib.f90` rows below are
DONE. `meds_fast_split.f90` is DONE for the transpiration↔uptake ledger (not the `INTEG_RK4` dispatch arm,
which is P2). `meds_fast_ark.f90` is DONE for the mass↔psi gather/commit boundary only — the tableau's own
`psi`-in-`column_state_t` representation, `plant_water_tendency`, and `column_prepass`'s other duties are
all UNCHANGED (deliberately out of scope; see the P0 implementation notes after §8). Every other row
(`meds_fast_rk45.f90`, `meds_fast_types.f90`'s `column_state_t.psi` retirement, `meds_fast_time_derivs.f90`,
`meds_fast_control.f90`'s `GRP_LEAF_W`/`GRP_WOOD_W`, `meds_fast_rk4_oracle.f90`, `meds_plant_hydraulics.f90`'s
`plant_water_tendency` retirement, the surface-water wiring, and the config plumbing) is still P1/P2 design-only.

**New artifacts (do not exist today):**
- `src/driver/meds_fast_rk45.f90` — the Cash–Karp 5(4) module (peer of `meds_fast_ark`): coefficients +
  embedded 4th-order b-vector + RKQS driver (or reuse of `adaptive_ark_march`'s control loop with a
  Cash–Karp stage), over `column_derivs`.
- `INTEG_RK4` selector (`meds_config`) + a third dispatch arm in `column_fast_step` (`meds_fast_split.f90:185`).
- `[fast].rk45_refresh_aero` config field + reader parse.
- `psi_from_water_content` — one `elemental pure` fn in `meds_hydr_lib` (the only constitutive addition).
- `GRP_LEAF_W` / `GRP_WOOD_W` tol groups in `meds_fast_control`.
- `leaf_water_mass`/`wood_water_mass` (internal) + `leaf_surf_water`/`wood_surf_water` (surface) states at
  the three state homes.

**Modified files:**

| file | current role | change |
|---|---|---|
| `src/core/meds_core_state_types.f90` | persistent cohort SoA (`psi(:,:)` `:102`) | add the 2 (P0) / 4 (P1) mass arrays at every lockstep site ψ touches: `cohort_alloc:305`, `site_free:281`, `cohort_ensure_capacity:390`, `move_alloc_block:437`, `cohort_reorder:529`, `copy_cohort_slot:591`, `init_cohort:677` |
| `src/core/meds_core_cohort_fusefiss.f90` | cohort fusion (`fuse_2_cohorts:164`; ψ leaf-area-weighted `:184`) | mass fuses **extensive** (nplant-weighted, conserve total — the carbon pattern `:207-214`), NOT leaf-area-weighted; ψ dropped/re-derived |
| `src/biophysics/meds_biophysics_types.f90` | `patch_biophys_t` (`psi:348`); `leaf_energy_env_t.leaf_water:226` (transient film) | add mass + surface arrays beside `psi`, thread `alloc_patch_biophys:448` / `ensure_patch_biophys_capacity`; the film becomes a persistent surface store |
| `src/driver/meds_fast_types.f90` | flat RK `column_state_t.psi:252` / `column_tend_t.dpsi_dt:298`; `column_budget_t` | replace ψ with `leaf/wood_water_mass` + `leaf/wood_surf_water`; `dpsi_dt` → mass derivatives; update `state_init/axpy/sub/err_diff` |
| `src/driver/meds_fast_time_derivs.f90` | `column_derivs` (`plant_water_tendency:239` → `f%dpsi_dt`); `surface_derivs` (`TAU_COND` sink `:141-156`) | drop `dpsi_dt`; add mass derivatives + `qloss`/`qwflux_wl` enthalpy + surface-water/film terms; **route `sf%cond` to dew, not out** |
| `src/driver/meds_fast_split.f90` | split stepper + dispatch gate (`:185`) | add `INTEG_RK4` arm; `solve_plant_water_batch:459` becomes the averaging pre-pass with frozen `<uptake>` feeding soil+wood; no longer byte-identical (gains closure); P0 proves the ledger here first |
| `src/driver/meds_fast_ark.f90` | ARK stepper; `column_prepass:847`; `build_column_frozen`; §8g leak `:192-194` | `column_prepass` gathers mass → ψ_leaf for gs; `build_column_frozen` freezes `<sapflow>`/`<uptake>`; drop `plant_water_tendency`/ψ-in-tableau; remove the `sf%cond` ledger terms; aero-refresh toggle |
| `src/driver/meds_fast_control.f90` | controller; 7 tol groups incl. `GRP_PSI:47`; `with_psi:182` | add `GRP_LEAF_W`/`GRP_WOOD_W`; retire `GRP_PSI` + `with_psi`; controller otherwise reused |
| `src/driver/meds_fast_rk4_oracle.f90` | test-only full-column RK4 over `column_derivs` | structural template the RK45 lifts; its `psi` stage-algebra terms → mass terms; stays an oracle |
| `src/plant/meds_plant_hydraulics.f90` | `solve_plant_water` matrix-exp (`<flux>:208-209`); `plant_water_tendency:407` | `solve_plant_water` **kept** as the averaging pre-pass (fluxes now frozen downstream); `plant_water_tendency` **retired** |
| `src/shared/functions/meds_hydr_lib.f90` | PV curves: `water_content:229`, `psi_from_rwc:198`, capacitance | add `psi_from_water_content` (composes existing pieces); rest reused |
| `src/biophysics/meds_vegetation_biophysics.f90` | `intercept_canopy_layer:192` (**unwired**); `veg_energy_diagnostic`; `veg_surface_fluxes` | wire interception as a persistent state; add wet/dry-fraction transp partition; enable film-evap/dew + its enthalpy |
| `src/shared/config/meds_config.f90` + `src/io/meds_config_io.f90` | `INTEG_SPLIT/ARK` (`:65-66`); `[fast]` fields | add `INTEG_RK4` + `rk45_refresh_aero`; extend `time_integrator` parse (`config_io:665`) + validation |

**Reused unchanged** (shown as ledger terms only): `meds_soil_water` `ground_evaporation`/`soil_evap` +
the below-wilting give-back template (`:134-148`); `meds_ground_biophysics` snow accumulate/sublimation/melt;
`meds_output_registry` `GRP_NUMERICS` work counters (the RK45 feeds the existing `work_integ_steps`/`_rej`).
Note: neither design doc adds output-registry rows for the new water stores — they are ledger-internal;
emitting them as diagnostics would be a separate, additive change.
