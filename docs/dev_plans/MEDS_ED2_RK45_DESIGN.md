# MEDS ED2-Faithful RK45 + Internal-Water-Mass Design

**Status:** **P0+P1+P2+P3+P4 IMPLEMENTED** (2026-07-24, on branch `feature/ed2-rk45-water-mass`, all
phases committed, not yet PR'd/merged to `main`): internal water MASS is the persisted plant-hydraulics
state on every path (split, ARK, and the genuinely-integrated RK45), every column water/energy ledger
closes to machine precision, canopy-surface water (interception/film-evap/dew) is wired on all three
integrators, the Cash–Karp RK45 stepper empirically holds ~5th order (including on a variable ARK's own
operator split degrades), the biomass fast/slow seam conserves water MASS (not ψ) across a slow-loop
biomass update with saturation-ceiling clamp bookkeeping for the one hazard that creates, and P4 extends
that same conservation discipline across the two OTHER slow-step events that move mass across the seam —
tree mortality (audited: already correct by the existing per-plant/`nplant` architecture, no code needed)
and leaf/root turnover (a real gap: net tissue loss now sheds water to the ground exactly like throughfall,
its own tracked variable, with energy riding the existing infiltration treatment for free) — see "P0",
"P1", "P2 gates 1 and 4", "P2c", "P3", and "P4" implementation notes at the end of §8 for the full
phase-by-phase history. Companion to `MEDS_NUMERICS_SCOPING.md` §12 (the
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

### P2 implementation notes (2026-07-23) — the RK45 stepper, in three sub-parts

**P2 part 1 — retire ψ from the shared ARK/RK45 tableau.** `column_state_t`'s `psi` field is gone; the
tableau now carries `leaf_water_mass`/`wood_water_mass` directly (§4), matching what P0 already made the
persisted representation. ARK gained a frozen-flux water closure to go with it: `column_be_stage` grew an
optional `sf_out`, `ark2_column_step` b-weights the two ESDIRK stages' own `transp_c` (`transp_bw =
(1-γ)·sf2%transp_c + γ·sf3%transp_c`) into a rewritten `advance_water_mass_full` that takes this array as an
explicit input rather than re-deriving an endpoint-only value that didn't match what the CAS itself
integrated — the "one flux, both sides of an interface" principle (§9) applied to the wood↔leaf mass
transfer specifically. Verified: ifx 36/36, nvfortran 36/36, `whole_water` closing to machine precision.

**P2 part 2 — advective enthalpy (ED2's `qwflux_wl`/`qloss`, §2).** Rather than adding genuine prognostic
leaf/wood energy state to the shared tableau (the literal reading of §2's "closing this needs
`d(leaf_energy)/dt`" language), a simpler mechanism was used: since leaf/wood stay diagnostic
(zero-capacitance) on the ARK/RK45 path, the advective terms fold in as *extra source terms* on the
existing algebraic energy balance instead. `build_column_frozen` computes a root-frac-weighted mean soil
temperature and, per cohort, an upwind temperature (`t_up_wl`, wood or leaf depending on sapflow's sign),
freezing **both** the mass flux (already frozen per P0) **and** the reference temperature at state^n —
avoiding circularity with `surface_derivs` (which produces those temperatures), at the cost of some
thermal-fidelity approximation (a documented, bounded imprecision, mirroring P0/P1's own deferred-precision
notes). `fro%surf%qwflux_wl`/`fro%surf%q_wood_net`/`fro%qloss_frozen` [W/m² ground] carry the result;
`q_wood_net := qloss − qwflux_wl` by construction, so the two together always sum to the full soil-side
debit regardless of how sapflow/uptake individually split it.

**P2 part 3 — the Cash–Karp stepper itself**, `meds_fast_rk45.f90`: the tableau of §6 exactly (6 stages,
`A21..A65`, 5th-order `B1/B3/B4/B6` committed via local extrapolation, embedded 4th-order `BS1..BS6`),
generalizing `adaptive_step_update` (`meds_numerics.f90`) with an optional `p_order` (default 1, preserving
every existing caller) so the step-size controller's `-1/(p+1)` exponent matches whichever embedded pair is
driving it — ARK's ARS(2,2,2) estimate is 1st-order, Cash–Karp's is 4th. **Diverges from §7.1's "faithful"
default**: aero is reused frozen from `build_column_frozen` (shared with ARK), not refreshed per stage —
`rk45_refresh_aero` was never implemented; the faithful/frozen accuracy-vs-cost comparison §7.1 recommends
measuring first remains undone. `column_fast_step_rk45` checks the whole-column water/energy budgets only
(§8 gates 2/3) — RK45 has no operator split at all, so there is no separate per-kernel "soil_water (rk45)"
check the way ARK needs one; a per-kernel `cas_co2` closure is also not yet tracked (would need per-stage
CO2-exchange accumulation `rk45_column_step` doesn't currently carry).

**Three bugs, found and fixed while validating against the whole-column ledger** (the first is RK45-local;
the other two are pre-existing, shared with ARK, and were only exposed — not introduced — by RK45's
stricter/newly-added tests):

1. **`state_sub` argument aliasing** (RK45-local, a genuine Fortran standard violation, not a compiler
   quirk): the embedded-error computation called `state_sub(y_out, y_err, n, nsl, y_err)`, aliasing `y_err`
   as both the `intent(in) b` and `intent(out) out` arguments of the same call — undefined behavior, since
   the compiler is permitted to assume `out` doesn't alias its other arguments. Symptom: the embedded error
   estimate did not shrink with `dt` at all across ~9 orders of magnitude (mathematically impossible for a
   correct embedded pair, whose difference must vanish as `dt→0`), so the adaptive controller never met
   tolerance except at a floor-clamped `dt`, driving 4094 substeps instead of the §6 estimate of ~3. Fixed
   by accumulating the 4th-order sum into a separate `y_4th` temporary before subtracting.
2. **`veg_energy_diagnostic`'s `drnet` double-counting the P2-part-2 advective term.** `drnet`'s formula
   (`abs_sw + abs_lw − lw_slope·(...)`) simply re-includes whatever was passed as `abs_sw` — so feeding it
   `abs_sw + qwflux_wl` (the original P2-part-2 call convention) made `coh_rnet` (and hence the
   `e_in` ledger) count the soil→leaf transfer as if it were *external* radiative input, on top of the
   genuine soil-side debit (`root_heat_sink += qloss_total`) — a real double-count. Fixed by giving
   `veg_energy_diagnostic` a new optional `q_extra` input that still shifts the `dt_temp` equilibrium (and
   hence `dh`/`transp`, which correctly reach the CAS) but is excluded from `drnet`'s formula, which must
   stay radiative-only for the boundary-flux identity to close.
3. **ARK's `column_be_stage` never added `qloss_total` to its own `root_heat_sink`** (only `column_derivs`,
   meds_fast_time_derivs.f90, had it) — invisible until bug 2 was fixed, because the *old*, double-counting
   `drnet` happened to inflate ARK's `e_in` by exactly the amount its soil was never debited, i.e. two bugs
   canceling by coincidence on that one path. Fixing bug 2 alone regressed ARK's whole-energy closure from
   machine precision to ~7.5e4 J; fixing this third bug (mirroring `column_derivs`'s treatment exactly, incl.
   the per-kernel `soil_enth_out` ledger) restored it.
4. **`column_derivs` never added infiltration/runoff/drainage's advected enthalpy to `root_heat_sink`**
   the way `column_be_stage` does (§3.3's boundary water-enthalpy advection) — so any nonzero background
   drainage left the soil-energy *state* unaware of an amount the `e_in`/`e_out` ledger still counted as
   crossing the boundary. Invisible to every prior ARK/split test (all effectively zero-drainage scenarios)
   and to RK45's own first pass (same reason); caught only once a test exercised a small but genuinely
   nonzero background drainage — `internal_energy_liquid`'s low reference temperature (see the P1 note
   above) turns even a tiny mass flux into a large-looking energy term, so the omission was not subtle in
   magnitude once triggered. Fixed by mirroring `column_be_stage`'s `e_infil`/`e_runof`/`e_drain` treatment
   in `column_derivs`. (Counted as a 4th, ARK-adjacent finding rather than folded into bug 3's list above,
   since it lives in the RK45-only RHS, not the shared `veg_energy_diagnostic` kernel.)

Diagnosing bugs 2–4 required verifying the whole-column identity term-by-term — per-stage CAS-ODE and
soil-ODE self-consistency (`k_i%d_cas_enthalpy·wcap == src_enth_i − atm_enth_i`, `Σ_k dedt(k)·dz(k) ==
g_top+geothermal−coh_qsoil−qloss_total`) both held to machine precision at every stage checked throughout,
which by itself ruled out a stage-computation bug and correctly pointed at the boundary-flux *reconciliation*
(the `e_in`/`e_out` construction) instead — the actual defect was two frozen boundary terms
(`e_infil`/`e_runof`/`e_drain`, `qloss_total`) that the ledger assumed but the state's own RHS never
received, not an error in any single per-stage rate.

**Verification:** ifx Debug 37/37 + nvfortran multicore 37/37 (`test_column_rk45.f90`, new: GPP parity vs.
the split, a 24-step dry-window physical/bounded march, and a 96-step dry diurnal `whole_water`/
`whole_energy` closure + substep-count check). `whole_energy` closes to ~5e-7 J/m² (vs. an initial ~5e4 J
before the fixes above); adaptive substeps land at 2 per `dt_fast` — consistent with §6's ~3-substep
stability estimate (down from 4094 before bug 1 was fixed). `whole_water` closes to ~1e-6 kg/m².

### P2 gates 1 and 4 (2026-07-23) — wet-precip test, order-of-accuracy test, and a 4th bug

Closing out the two gates the P2-part-3 pass above left unexercised.

**A 4th bug, found by the new wet-precipitation test under REAL (not just incidental background)
drainage:** `column_fast_step_rk45`'s outer `e_in` added `forc%precip*dt_fast*u_liq(cas_temp)` ON TOP of
`e_in_acc` (the b-weighted accumulation of `rk45_column_step`'s own per-substep `e_in`, which already
carries `fro%infiltration*u_liq(rain_temp)` — the SAME frozen quantity feeding `column_derivs`'
`root_heat_sink(1)`, i.e. what the soil state actually receives). Since `infiltration ≈ precip` whenever
runoff is small (the common case) and `rain_temp ≈ cas_temp` (both are the same Act-1 reference), this
double-counted nearly the FULL infiltrating share — a residual of ~7e4 J/m² per step at a modest
continuous rain rate (8e-5 kg/m²/s), invisible in the P2-part-3 pass because its only precip-adjacent
scenario was a tiny incidental background drainage where `infiltration=0` made the redundant term a
no-op. Mirrors ARK's own `whole_energy` ledger exactly: `acc%whole_enth_in` (with `e_infil` folded into
`bf%whole_enth_in`) is used DIRECTLY, with no further outer precip addition. Fix: `e_in = e_in_acc`, full
stop — `w_in` keeps its `forc%precip*dt_fast` term unchanged (there is no matching double on the water
side, since `w_out_acc` has no infiltration-side counterpart to double against).

**Wet-precipitation test** (`test_column_rk45.f90`, mirrors `test_column_ark`'s Test F exactly:
`forc%precip = 8.0e-5` continuous over a 96-step diurnal march): unlike ARK's wet test (gated at a
"lagged-ponding operator-split" WATER tolerance, per that test's own docstring), RK45's `whole_water`
closes to the SAME tight, non-inflated tolerance as its dry test — no split-vs-continuous mismatch to
tolerate, since soil water is genuinely integrated. After the double-count fix, `whole_energy` closes to
~4e-7 J/m² even under sustained rain.

**Order-of-accuracy test** (§8 gate 1; `test_column_derivs.f90`'s new `test_rk45_order`, alongside the
existing `test_ark2` order test it mirrors): self-convergence of the committed 5th-order solution against
step size, via a fine reference (`dt=2400/1024 s`, 1024 steps) vs. three halvings (`dt=300/150/75 s`,
comfortably under the ~360 s single-step stability estimate). A 5th-order method's error collapses to
double-precision noise (~1e-13) at far coarser `dt` than a 2nd-order one, so — unlike `test_ark2`'s
`dt=200/8..32 s` range — the step sizes had to be widened substantially before the truncation term cleared
the roundoff floor. Two variables checked: (a) CAS CO2 (a decoupled affine ODE, mirroring `test_ark2`'s own
"clean tableau order" baseline) — observed order ≈4.9; (b) the soil-top temperature — the SAME variable
`test_ark2`'s part (a2) shows degraded to ~1.2 order under ARK's operator split — observed order ≈5.8 to
≈6.1 for RK45, i.e. **no order reduction on a genuinely coupled variable**, empirically confirming §6's "no
operator split at all" design claim rather than just asserting it.

**Verification:** ifx Debug 37/37 + nvfortran multicore 37/37 (both `test_column_rk45` and
`test_column_derivs`).

### P2c (2026-07-24) — canopy-surface water wired into the shared ARK/RK45 tableau

Closes the last deferred item: interception/film-evap/dew (§3.4/P1, previously split-path-only) now runs
under ARK and RK45 too, via `column_state_t`'s new `leaf_surf_water`/`wood_surf_water` [kg/m² ground] —
the surface-film analogue of §4's internal water mass, added at the same three homes (`column_tend_t` gets
matching `d_*_surf_water`; `surface_frozen_t` gets the frozen wetted fraction `f_wet_c` and film-evap
conductances `g_film_f`/`g_film_w`; `column_frozen_t` gets the frozen interception rates `intercept_leaf`/
`intercept_wood`; `surface_tend_t` gets `film_evap_leaf`/`film_evap_wood` outputs, mirroring `transp_c`'s
role for the internal-mass ODE).

**Design, in three pieces:**
1. **Interception is frozen ONCE per `dt_fast`** in `build_column_frozen`, mirroring `sapflow_frozen`/
   `uptake_frozen`'s "one frozen number, no per-stage re-solve" convention (§6's stability argument) rather
   than the split path's own per-Picard-pass re-run — a single height-sorted `intercept_canopy_layer` sweep
   (no `aero` dependency, so it runs *before* `column_prepass`) converts the split path's one-shot bucket
   update into an equivalent frozen RATE (`(post-interception − pre-interception) / dt_fast`; exact under
   explicit Euler over the full step). Its own `sigma_w` output becomes the frozen `f_wet_c`. Critically,
   **the soil hydrology's `hforc%precip_ground` now uses `throughfall_total`** (precip minus what the
   sweep caught), not raw `forc%precip` — feeding the unreduced total while also crediting interception
   would create water from nothing.
2. **Film evaporation is genuinely per-stage**, via `veg_energy_diagnostic`'s existing wetted-fraction
   extension (§3.4/P1): `surface_derivs` now computes `le_slope_wet`/`le_ref_wet` from the frozen
   conductance but state-dependent `dqdt`/`qsat_c−qcas`, mirroring exactly how the dry/stomatal pathway
   splits `g_tr_f` (frozen) from its live-state terms. `film_evap_leaf`/`wood` feed `coh_qw`/`src_vap`
   (mass+energy into the CAS) the same way `transp_c` does.
3. **RK45 integrates `leaf_surf_water`/`wood_surf_water` genuinely** (`column_derivs`: `d_*_surf_water =
   intercept_* − film_evap_*`, no operator split, matching how it treats internal mass). **ARK
   operator-splits it out** exactly like internal mass: passed through unchanged across the ESDIRK stages
   (`state_init`/`state_axpy`/`state_accum`/`state_sub`/`state_extrap`/`state_err_diff` all extended),
   then committed once over the full `dt` by a new `advance_surf_water_full` (mirrors
   `advance_water_mass_full`) using the b-weighted per-stage `film_evap_leaf`/`wood` (captured via the
   existing `sf_out` mechanism `advance_water_mass_full`'s `transp_bw` already relies on). No new
   `process_mask_t` field — surface water rides the existing `hydraulics` mask entry (no test scenario
   needs it reduced independently of internal mass yet).

**A real bug, found by the new canopy-water tests (not anticipated in design): mass fabrication from
over-evaporation.** `film_evap` is driven by a state^n-frozen conductance (`g_film_f`/`g_film_w`), oblivious
to the store depleting mid-step — nothing prevented it from draining more than was actually present,
driving `leaf_surf_water`/`wood_surf_water` **negative**. The consequence wasn't just an unphysical sign:
`intercept_canopy_layer`'s own internal `max(·,0)` floor on next step's *starting* bucket silently
"fixed" the negative carry-in, materializing water that was never really lost — a genuine mass leak, not a
cosmetic one, caught by `whole_water%n_fail` (not `whole_energy`, which already tolerates a bounded gap
here). Two-part fix, both needed:
1. **Rescale the frozen film-evap conductance by availability** in `build_column_frozen` — the same
   pattern `uptake_frozen` already uses when the soil can't honour the full request (`scale =
   realized/requested`). A preliminary `surface_derivs` call (`sf0`, already computed there for the
   hydraulics pre-pass) gives a state^n potential `film_evap`; if `potential·dt_fast` would exceed
   `current store + this step's frozen interception`, `g_film_f`/`g_film_w` are scaled down accordingly
   and `sf0` is re-evaluated so its own `transp_c` (feeding the hydraulics batch solve) sees the final
   conductance. This is an *approximation*, not a hard guarantee — evaporative demand can still grow
   through the step as `tcas`/leaf temperature evolve past their state^n values.
2. **A `surf_deficit` backstop**, symmetric to the pre-existing `surf_overflow` (capacity-ceiling)
   bookkeeping: floor the committed `leaf_surf_water`/`wood_surf_water` at 0 and credit the shortfall as a
   *negative* addition to the whole-column ledger's outflow (flooring a negative store *up* to 0 makes the
   store appear to gain, so the ledger's outflow must shrink to match — the exact sign-mirror of
   `surf_overflow`, which must *grow* the outflow when a store is capped *down*). This guarantees exact
   `whole_water` closure regardless of how good the piece-1 approximation is; a residual gap of the same
   kind only shows up in the (already bounded, not machine-precision) `whole_energy` ledger. Needed for ARK
   too: `advance_surf_water_full` originally floored at 0 *internally*, silently discarding the same
   information the caller's bookkeeping needs — the floor was moved out to the caller.

**Verification:** ifx Debug 37/37 + nvfortran multicore 37/37. New tests: `test_column_ark.f90`'s Test G
and `test_column_rk45.f90`'s Test E (both mirror `test_column_dynamics.f90`'s own RUN 6 for the split path)
— a 96-step diurnal march with a 5-step morning rain pulse (istep 20–24, `precip=5e-5` kg/m²/s) under
`canopy_water_on=.true.`: the pulse is genuinely intercepted (peak film ≈0.107 kg/m² for both integrators,
matching each other closely as expected since they share `build_column_frozen`), `whole_water` closes
exactly (`n_fail==0`), `whole_energy` stays within the same documented bound (`<5e6` J/m²) the split path's
own RUN 6 uses (≈1.2e3 J/m² observed for both ARK and RK45 — comfortably inside it). With the flag off,
every existing test is unaffected (all new fields/arithmetic are zero and gated behind `canopy_water_on`
per the P1 nvfortran lesson), confirming the default path stays byte-identical.

**Deferred, not part of this pass:** the ARK §8g `TAU_COND` condensation sink is a separate, pre-existing
mechanism, left untouched rather than unified with the new dew pathway (`film_evap<0`) — avoid enabling
both `canopy_water_on` and `cas_condensation` simultaneously for now. A dedicated `process_mask_t` field
for surface water (vs. riding the `hydraulics` entry). Order-of-accuracy / self-convergence coverage for
the surf_water ODE specifically (the existing `test_rk45_order` in `test_column_derivs.f90` predates this
feature and does not exercise it).

---

### P3 (2026-07-24) — the biomass fast/slow seam, mass-conserving (user-directed revision of §9)

**Design decision (user, overriding this doc's own earlier draft resolution in §9):** the seam conserves
water MASS, not water potential. `leaf_water_mass`/`wood_water_mass` are left completely UNCHANGED across
a slow-loop biomass update — no seam code touches them in the growth direction at all. Since capacity
`W_sat = water_sat·biomass` grows with biomass while mass stays fixed, the next fast-loop touch diagnoses
a lower (drier) rwc/ψ automatically, via the same `psi_from_water_content` every dt_fast already calls.
That lower ψ is the physically-correct signal that draws more water from the soil over subsequent fast
steps — "redistribute the water mass to the new tissue biomass ... plants will absorb water from soil in
next timestep to make it up," in the user's framing. This is also the mass-conservation-consistent choice:
P0's entire premise was replacing ψ (intensive, capacity-blind) with mass (extensive, conserved) as the
fast-loop prognostic specifically so the column ledger closes; making ψ the seam-continuous quantity
instead (this doc's original §9 resolution) would have silently created/destroyed mass at exactly the
seam, undermining that premise. **In the growth direction this needed zero new code** — nothing in
`meds_vegetation_dynamics`'s carbon growth path ever touched `leaf_water_mass`/`wood_water_mass`, so mass
was already conserved across the seam by omission; the fix was recognizing that the doc's own drafted
"resolution" would have been a regression, not a gap.

**The one real hazard mass conservation creates: a capacity SHRINK.** Biomass is not monotonic — the
phenology dormant-canopy leaf snap-to-bare (`meds_vegetation_dynamics%leaf_shed_amount`) can shed the
*entire* current `leaf_carbon` pool in one slow step (`shed = pool`, not a gradual decay), so `bleaf` can
collapse far faster than "small every day" in that one case. A carried-forward mass that was reasonable
against yesterday's (larger) capacity can then exceed today's (collapsed) ceiling — `rwc > 1`, a tissue
state with no valid inverse on the turgid PV branch (not a NaN, since the turgid formula is linear+1/r
and does not blow up, but a non-physical extrapolation past full saturation). Two-part fix, mirroring the
project's established clamp/bookkeep discipline (P2c's `surf_overflow`/`surf_deficit`):
1. `clamp_water_to_capacity(w, water_sat, biomass)` — one new `elemental` function in `meds_hydr_lib.f90`
   (beside `psi_from_water_content`): `w_capped = min(w, water_sat·biomass)`. A strict no-op whenever
   biomass grows (the common case); caps exactly at the new ceiling whenever biomass has shrunk enough
   that the old mass no longer fits, including collapsing cleanly to exactly 0 for a fully-bare cohort.
2. Wired into `meds_fast_dynamics.f90`'s daily per-cohort gather, in the `else` branch of the existing
   lazy-init check (the branch that runs for every already-seeded cohort, i.e. every day after the first):
   `leaf_water_mass`/`wood_water_mass` are clamped against that day's freshly-gathered `coh%bleaf`/
   `coh%bsap+coh%broot` before being copied into the fast loop's working bundle. The released excess is
   **not** carried forward into any ledger — there is no slow-timescale water budget in this model to
   bookkeep it into, and the fast loop's own `whole_water` ledger is unaffected either way (it spans one
   `dt_fast`, entirely *after* this gather, so the clamp is invisible to it by construction — confirmed by
   the driver-level test below, which asserts `n_fail==0` on the very call that triggers the clamp).

**A pre-existing, branch-unrelated nvfortran bug found by the new tests (not a P3 regression — confirmed
by reverting to the pre-P3 commit and reproducing the identical crash).** `phi_inverse`'s bisection
(`meds_hydr_lib.f90`) hardcodes its upper bracket at `0.0_wp`; for any `wood_kexp` outside the closed-form
`{1,2}` set, that endpoint evaluation reaches `flux_potential`'s quadrature branch with `r=0` on its very
first call, every time, computing `kirchhoff_integrand(0.0) = 1/(1+0**kexp)`. Under nvfortran's strict
`-Ktrap=fp`, `0.0**kexp` for a general (non-integer-recognized) real `kexp` apparently lowers through a
`log`-based path and traps on `log(0)` — even though the overall mathematical result (`0**kexp=0` for
`kexp>0`) is perfectly well-defined and ifx's own `-fpe0` does not trip on it. Invisible on every other
existing test because `wood_kexp=2.0` (the closed-form, no-power-operator branch) is this project's
universal default, and the one test that already exercised a general `kexp` (`test_hydro_table`'s sweep)
happens to never land exactly on `r=0`; `phi_inverse`'s hardcoded bracket is the only call site that
lands there deterministically. Fixed by guarding the `r<=0`/`u<=0` case explicitly (`f=1.0`/`y=1.0`, the
correct limit) in `plc_retained` and the `flux_potential`-local `kirchhoff_integrand`, avoiding the power
operator entirely at the trap-prone point; `dplc_dpsi` shares the same latent hazard (plus a further `0**0`
case at `kexp=1`) but is currently uncalled anywhere, so it is left as-is with a comment flagging it rather
than guarded blind. Diagnosed via the same bisection-with-stderr-markers discipline this design doc's
earlier phases established (stdout is fully block-buffered across an unhandled SIGFPE, so `write(0,...)`+
`flush(0)` between calls was needed to localize it at all).

**New tests:** `test_plant_hydraulics.f90` gained `test_biomass_seam` (library-level: growth at fixed mass
lowers ψ_wood; feeding that lower ψ into `solve_plant_water` against an identical fixed soil boundary
condition draws strictly more `root_uptake` than the pre-growth ψ does — the end-to-end physical claim)
and `test_seam_capacity_clamp` (`clamp_water_to_capacity` is an exact no-op on growth; caps exactly at the
new ceiling on a synthetic 10x shrink with a non-negative released excess; collapses a fully-bare tissue
to exactly 0). `test_fast_loop.f90` gained a driver-level integration test: a forced 100x `leaf_carbon`
collapse between two `fast_dynamics` calls (standing in for one slow-loop day of the phenology snap-to-
bare, without needing a full `vegetation_dynamics` call) drops `leaf_water_mass` to well under half its
pre-collapse value — proving the clamp is actually wired into the real daily gather, not just correct in
isolation — while the whole-column budgets stay closed on that same call.

**Verification:** ifx Debug 37/37 + nvfortran multicore 37/37 (both suites unchanged in count — P3 added
assertions to two existing test binaries rather than new ones). Flag-off/default-path behavior is
byte-identical (the seam was already a no-op by omission; the clamp only ever activates on a capacity
shrink, which no existing test before this one constructed). P3 was the last item scoped by this design
doc; only the separately-tracked P3 cavitation/saturation clamp for *within-fast-loop* dynamics (as
opposed to the seam-specific ceiling clamp above) was not part of this pass — no test in the existing
suite has been observed to need it, and this pass did not go looking for one.

---

### P4 (2026-07-24) — conserving carbon/water/energy across mortality and leaf/root turnover

P3 closed the seam for *biomass growth* (mass-conserving, with a ceiling clamp for the one hazard that
creates). The user then raised two further slow-step events that also move carbon/water/energy across the
fast/slow boundary and asked for the same treatment: **tree mortality** (a cohort's `nplant` drops) and
**leaf/root turnover** (a cohort's tissue biomass drops without a whole individual dying). The two turned
out to need very different responses — one already worked, the other needed real new wiring.

**1. Mortality — audited, no code needed.** `update_cohort_states_kernel` (`meds_core_state_types.f90`)
applies `nplant(i) *= exp(dln*dt)` and touches **only** `dbh/height/basal_area/agb/leaf_area/leaf_carbon/
fineroot_carbon/wood_carbon/nonstructural_carbon/nplant` — every per-plant field (`leaf_water_mass`,
`wood_water_mass`, `leaf_temp`, `wood_temp`, `leaf_surf_water`, `wood_surf_water`) is left completely
alone. Since every column-wide total is `Σ nplant(i)·field(i)`, a shrinking `nplant` alone already reduces
water and energy out of the tracked system exactly in proportion to how many individuals died — with zero
risk of double-touching a per-plant field that has no death-specific meaning to begin with (a *surviving*
individual's own water content does not change just because a neighbor died). Carbon is handled too, and
already had its own explicit pathway: `accumulate_mortality_litter` computes `died_nplant · pool(j)` for
every carbon pool and routes it into the same per-patch `litter_input_t` the turnover pathway below
feeds, when `[soil_carbon].soil_carbon_on`. **Conclusion: both requirements were already met by the
existing architecture** — the per-plant/`nplant`-multiplier convention this whole codebase already uses
IS the mass-conservation mechanism for mortality; nothing needed adding. (An explicit multi-day, whole-site
water+energy ledger that could *assert* this was considered and rejected as unnecessary scope: no such
ledger exists at this timescale, per-plant fields are provably untouched by inspection above, and P3 already
established that the fast loop's own `whole_water`/`whole_energy` checks are — by design — blind to
anything that happens between two `fast_dynamics` calls.)

**2. Leaf/root turnover — real gap, fixed.** Unlike mortality, a turnover event does NOT touch `nplant` —
individuals stay alive, but the phenology/turnover pathway (`meds_vegetation_dynamics%cohort_carbon_demand`
→ `update_biomass_turnover`) can still shed *tissue* biomass (leaf, fine root) faster than growth replaces
it on a given day. P3's ceiling clamp alone was too crude for this case: it only intervened once the
carried-forward mass exceeded the NEW (smaller) capacity, silently discarding the excess with no
destination. The user's refinement: **check the NET pool change at the end of the slow step.** A net GAIN
needs nothing further (P3's own mass-conserving seam already handles it — capacity grows, mass sits still,
ψ reads lower next touch). A net LOSS should shed water in the SAME proportion as the carbon lost — keeping
the REMAINING tissue's rwc unchanged, exactly like a leaf carrying its own share of water away when it
abscises — and that shed water should reach the ground **the same way throughfall does, but as its own
distinctly-tracked variable**, so nothing is silently lost from the whole-column ledger.

**Design.** A new subroutine `shed_turnover_water(site, cfg, npp)` (`meds_vegetation_dynamics.f90`), called
right after `compute_carbon_allocation` returns this step's `npp` (and, critically, BEFORE
`update_cohort_states` commits the new pools — the shed *fraction* needs the PRE-step `leaf_carbon`/
`fineroot_carbon` as its denominator):
```fortran
if (npp%leaf(j) < 0.0_wp) then
   shed_frac   = min(1.0_wp, -npp%leaf(j) / max(cohort%leaf_carbon(j), tiny_num))
   shed_amount = cohort%leaf_water_mass(j) * shed_frac
   cohort%leaf_water_mass(j) = cohort%leaf_water_mass(j) - shed_amount
   patch_total(ip) = patch_total(ip) + cohort%nplant(j) * shed_amount
end if
```
(fine roots mirror this into `wood_water_mass` — the lumped leaf+wood 2-node model has no separate root
water pool, so `wood_water_mass`'s shed fraction uses fine-root carbon's own fractional loss as the proxy;
`bsap` doesn't shed today, and it's itself still an MVP placeholder — see `meds_fast_dynamics.f90`'s
`0.10*wood_carbon` stub — so a more exact split isn't yet meaningful). The nplant-weighted per-patch total
is divided by `cfg%dt_slow` and SET (not accumulated — it's "today's rate," fully replaced every step, so a
patch with no shedding this step correctly reads exactly 0) onto a new persistent field, `site%patch%
shed_water_rate` [kg/m² ground/s] — threaded through the FULL patch lockstep (`patch_alloc`/
`patch_ensure_capacity`/`site_free` in `meds_core_state_types.f90`; `sort_patches`/`fuse_2_patches`/
`patch_compact`/`apply_patch_disturbance` in `meds_core_patch_fusefiss.f90`), **area-weighted on fusion
like `age`** (not extensive-summed like `soil_carbon`/`xi_accum` — it's a per-area rate, and both fusing
patches' rates are already real, freshly-computed numbers for THIS step by the time patch fusion runs
later in the same `vegetation_dynamics` call) and **reset to a fresh 0 on a new disturbance gap, also like
`age`** (a brand-new gap patch has no cohorts of its own that contributed to today's rate — unlike
`soil_carbon`/`xi_accum`, which are genuine material inherited from the disturbed area).

The frozen rate is bridged into the fast loop exactly like the soil-carbon pool: `patch_biophys_t` and
`column_forcing_t` each gain a `shed_water_rate` scalar; `meds_fast_dynamics.f90`'s daily gather copies
`site%patch%shed_water_rate(ip)` into `bio%shed_water_rate` once per patch (beside the existing `bio%
soil_carbon` copy), and `fill_forcing` (now taking `bio` as an extra argument) copies it into `forc%
shed_water_rate`, constant across every `dt_fast` sub-step of that day. From there, **one line each** in
`meds_fast_split.f90` and (shared by RK45) `meds_fast_ark.f90%build_column_frozen`:
```fortran
hforc%precip_ground = hforc%precip_ground + forc%shed_water_rate
```
puts the water through the SAME infiltration/ponding/runoff solve throughfall gets — not a hand-rolled
parallel mechanism — and one line each in the three integrators' whole-column ledger assembly
(`w_in`/`acc%whole_wat_in`) recognizes it as a boundary input distinct from `forc%precip`, so the ledger
stays exact rather than seeing water "appear from nowhere."

**No separate energy wiring was needed — a deliberate simplification, not an oversight.** Once the shed
water is mixed into `hforc%precip_ground`, its energy automatically rides the SAME `e_infil = fro%
infiltration · internal_energy_liquid(fro%rain_temp)` term every OTHER infiltrating input already gets —
exactly "similarly as precipitation," not a more precise mechanism than precipitation itself receives (this
codebase already treats ALL infiltrating water at one reference temperature; see P1/P2c's own documented
`rain_temp` simplification for canopy-surface water). Tracking the shed water's OWN temperature (leaf/wood
temp at shed time) all the way through a parallel energy pathway was considered and rejected: it would only
matter for a term this small, would require either splitting `column_hydrology_flux`'s single infiltration
input into two temperature-tagged streams or accepting a second, harder-to-justify approximation elsewhere,
and the existing rain_temp treatment is already the established, documented precedent for this exact class
of imprecision. **Empirically confirmed, not just argued**: the new tests below show whole_energy closing
to the SAME machine-precision bound as every other test in these files (worst residual `<1e-6` J), with
`forc%shed_water_rate` active and zero energy-specific code added.

**Tests.** `test_carbon_growth.f90` gained a `shed_turnover_water` unit test (net loss sheds water exactly
proportional to the carbon loss on both leaf and fine-root paths; net gain touches nothing and resets the
patch rate to exactly 0; a full net-loss-of-the-whole-pool sheds ALL the water). `test_disturbance.f90`
gained one assertion that a fresh treefall gap's `shed_water_rate` is exactly 0. `test_column_ark.f90`
gained `test_ark_shed_water` AND `test_split_shed_water` (both integrators share that file's fixture);
`test_column_rk45.f90` gained `test_rk45_shed_water` — all three run a 96-step diurnal march with a constant
`forc%shed_water_rate` (precip held at 0) and assert: the soil column wets from shed water ALONE, and both
`whole_water` and `whole_energy` still close with `n_fail==0`.

**Verification:** ifx Debug 37/37 + nvfortran multicore 37/37 (test COUNT unchanged — P4 extended four
existing binaries). Default path (`shed_water_rate` always 0 unless a cohort actually loses tissue biomass
net this step) is unaffected in every pre-existing test — confirmed by the full suite passing unchanged.

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
- **Biomass fast/slow seam (new leak the ψ-state was immune to) — REVISED, see P3 implementation notes.**
  `W = capacity(biomass)·rwc(ψ)`; biomass changes on the slow step, so a persistent `leaf_water_mass`
  carried across a biomass update reads against a *different* capacity next touch. An earlier draft of
  this doc resolved this by making ψ the seam-continuous quantity (diagnose ψ from mass before the slow
  step, re-derive `mass = capacity(new biomass)·rwc(ψ)` after) — but that re-derivation silently creates
  or destroys mass at exactly the seam, which is backwards for a design whose entire P0 motivation was
  making mass (not ψ) the thing that conserves. **Resolution (final): mass is the seam-continuous
  quantity, not ψ** — `leaf_water_mass`/`wood_water_mass` simply carry forward UNCHANGED across a
  biomass update; a capacity GROWTH is then read as a slightly lower rwc/ψ next touch (drier), which is
  the physically-correct signal that draws more water from the soil over subsequent fast steps, not a
  defect to patch over. The only guard needed is the opposite direction: a capacity SHRINK (e.g. the
  phenology dormant-canopy leaf snap-to-bare) can leave more mass than the new ceiling admits — a
  non-reachable, supersaturated tissue state — so that direction alone is clamped, with the excess
  bookkept rather than silently retained (see `clamp_water_to_capacity`, P3 notes). (Within a fast day
  biomass is constant, so no in-loop leak either way.)
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
