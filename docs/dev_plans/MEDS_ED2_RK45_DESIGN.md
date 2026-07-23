# MEDS ED2-Faithful RK45 + Internal-Water-Mass Design

**Status:** design-only (2026-07-23). Companion to `MEDS_NUMERICS_SCOPING.md` §12 (the baseline-schemes
round) and §12.6 (the verified ED2 DTLSM frozen/integrated/refreshed catalogue). This doc specifies the
FIRST of the four baseline schemes — a faithful reproduction of ED2's adaptive Cash–Karp **RK45** — and
the state-vector change it is paired with: **leaf/wood internal water mass become prognostic states so the
column water budget closes within each `dt_fast`.**

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

**The single decision that makes both work: FREEZE the hydraulic fluxes.** ED2 computes the sapflow
(`wflux_wl`) and root uptake (`wflux_gw_layer`) **once per DTLSM** in `plant_hydro_driver` and holds them
constant across every RK stage. This is not incidental — it is what makes the whole thing tractable:

- **ED2-faithful** by definition (§12.6).
- **Non-stiff for an explicit RK45.** If the fluxes were recomputed from the evolving mass each stage,
  `flux = conductance·Δψ(mass)` reintroduces the ~17 s hydraulic RC eigenvalue into the explicit stability
  region — the exact stiffness MEDS removed by splitting ψ out. With the fluxes frozen, `d(wood_water)/dt`
  is a constant and `d(leaf_water)/dt = frozen_sapflow − transp` is driven only by the CAS-humidity
  timescale (~130 s), so the explicit RK45 stays stable at its CAS-limited step (~360 s → ~3 substeps).
- **Single-definition closure.** The *same* frozen uptake value debits the soil and credits the wood; the
  *same* frozen sapflow debits the wood and credits the leaf. There is exactly one number per interface,
  so closure is arithmetic, not approximate.

An earlier design draft recomputed the fluxes from evolving mass. That version is **both non-faithful and
stiff, and its closure is broken** (two incompatible uptake definitions). Do not build it. §9 records why.

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

**Flux graph**, per m² ground, one consistent flux per interface:

```
soil(k1) --wloss--> wood_water --wflux_wl--> leaf_water --transp--> CAS_vapour --ET--> atmosphere
```

**Store tendencies:**

```
d(soil_water·ρ)/dt   = −wloss            (+ infiltration/drainage boundary terms, unchanged)
d(wood_water)/dt     =  wloss − wflux_wl
d(leaf_water)/dt     =  wflux_wl − transp
d(CAS_vapour)/dt     =  transp − ET      (ET = CAS↔atm vapour flux)
```

**Closure — sum telescopes:**

```
d/dt(soil + wood_water + leaf_water + CAS_vapour)
   = −wloss + (wloss − wflux_wl) + (wflux_wl − transp) + (transp − ET)
   = −ET
```

Only the atmospheric exchange `ET` (and the soil boundary terms) leaves the column. `wloss`, `wflux_wl`
and `transp` each cancel between the store they leave and the store they enter — **because each is a
single number used on both sides.** This is the closure the current MEDS ARK lacks: today the soil sink is
fixed from the state-ⁿ demand while the stages re-evaluate transpiration and the difference is dropped
(`meds_fast_ark.f90` inflates the `whole_water` tolerance to `max(1e-3, |fro%uptake|·dt_fast)` to hide
it). With internal water mass as a state, the plant capacitance **absorbs** the uptake↔transpiration
mismatch explicitly, so that tolerance inflation is removed and replaced by a real machine-precision check.

**Why frozen-uptake / refreshed-transpiration still closes** (adversarially confirmed): closure needs only
that each interface use *one* flux on *both* sides within the step — not that the fluxes match each other.
`wloss` and `wflux_wl` are frozen constants (identical on both sides trivially); `transp` is refreshed but
the *same* refreshed value debits the leaf store and credits the CAS. The internal stores carry whatever
imbalance results. Closure is exact regardless of the freeze.

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

## 5. The frozen-flux pre-pass (MEDS's `plant_hydro_driver`)

Once per `dt_fast`, before the RK march, from state-ⁿ potentials:

1. Diagnose ψ_leafⁿ, ψ_woodⁿ from `leaf_water_mass`ⁿ, `wood_water_mass`ⁿ.
2. Diagnose ψ_soil from the soil column (root-fraction-weighted, as today).
3. Compute the **frozen** fluxes: `wflux_wl = edge_cond(ψ_wood,ψ_leaf)·(ψ_wood − ψ_leaf)` (sapflow);
   `wloss = rhizosphere_cond·(ψ_soil − ψ_wood)` per root layer (uptake). nplant-weight both to per-m²
   ground. Floor uptake at ≥ 0 per layer (no hydraulic redistribution, as today).
4. Hand the frozen `wflux_wl`, `wloss(k)` to the stage RHS and to the soil sink — **the same values**.

This is the seam that guarantees §3's single-flux-per-interface property. It replaces `solve_plant_water`'s
ψ-exponential; the exact-exponential machinery is retired for this scheme (it was the *integrator* for ψ,
now unnecessary because ψ is diagnostic and the flux is frozen).

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
2. **Hydraulics representation.** ED2 has no ψ state; MEDS had ψ as the state. This design moves MEDS to
   ED2's representation (mass state, ψ diagnostic, flux frozen) — the change this whole doc specifies.

---

## 8. Config, verification, phasing

**Config (`[fast]`):** `time_integrator = "rk45"` (new `INTEG_RK4`); `rk45_refresh_aero` (default true);
the RK45 reuses `error_level`/`step_controller`/tolerance knobs. All default-off for existing configs.

**Acceptance gates** (per §12.3 — NOT production RMSE at 900 s, which is freeze-limited):
1. **Order-of-accuracy** on a frozen-forcing / manufactured-solution test — the RK45 must show ~5th order
   (embedded 4th) where the coefficient freeze is zero by construction.
2. **Machine-precision column water closure** — the new `whole_water` ledger (soil + wood_water +
   leaf_water + CAS vapour vs boundary ET/precip) closes to ~1e-13, and the inflated ARK tolerance is
   removed. This is the headline deliverable — gate on it hard.
3. **Machine-precision column energy closure** — with the `qloss`/`qwflux_wl` enthalpy terms carried.
4. **Cost/stability** via the §5.3 work counters (substeps, rejections) — bounded, low rejection.
5. ifx 36/36 + nvfortran multicore; existing paths (split/ARK) byte-identical (RK45 is purely additive).

**Phasing:**
- **P0 — mass state, ψ diagnostic, no new integrator.** Add the mass fields + `psi_from_water_content`;
  convert the *existing* split/ARK hydraulics to mass-prognostic with frozen fluxes; prove the closure
  ledger on the split path first (simplest). This delivers the water-closure win independent of RK45.
- **P1 — the Cash–Karp RK45** over the P0 state, with the enthalpy coupling and the aero toggle.
- **P2 — the biomass fast/slow seam** (§9) and cavitation/saturation clamp bookkeeping hardened.

---

## 9. Pitfalls (the adversarial pass — read this before coding)

An earlier design draft was refuted here; each item is a real trap with its resolution.

- **FATAL if fluxes are refreshed, not frozen.** Recomputing `wflux_wl`/`wloss` from the evolving
  mass-diagnosed ψ each stage (a) diverges from ED2 (which freezes them) and (b) reintroduces the ~17 s
  hydraulic RC eigenvalue into the explicit RK45 stability region → tiny substeps or blow-up.
  **Resolution: freeze the fluxes** (§1, §5). This is the load-bearing decision of the whole design.
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
  :879-897, :2096-2123); the nplant-weighting; that MEDS has every constitutive piece and needs only
  `psi_from_water_content`; the three state homes and every lockstep touch point; the current non-closure
  and the inflated tolerance that hides it; the fusion weighting difference.
- **Adversarially confirmed:** the closure telescoping; that frozen-uptake/refreshed-transp does not leak;
  that the naive refreshed-flux design is fatal on both fidelity and stability.
- **Assumed / to measure:** the ~3-substep stability estimate (measure on a wet + a dry forced window);
  the aero refresh-vs-freeze accuracy delta (§7.1 — measure before committing); the RK45's 5th-order gate
  (on a frozen-forcing test). The three failed research scouts (ED2-derivatives, closure, stepper) were
  re-covered by direct reading + the verifier pass, so no gap remains — but the RK45 stepper module itself
  has not been prototyped.
