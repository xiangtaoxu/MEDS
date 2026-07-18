# MEDS Slow-Timescale Biogeochemistry — Design Document

*Definitive implementation plan for the slow (daily) carbon — and optional nitrogen — cycle: litter and
coarse-woody-debris inputs from mortality and tissue turnover feed a prognostic multi-pool soil-carbon
column; ED2-faithful CENTURY decomposition returns the pools' carbon to the atmosphere as heterotrophic
respiration, which drives the existing canopy-air-space CO2 twin (`meds_column_co2`) and closes NEE. The
whole decomposition network is organized as one **matrix ODE** `dX/dt = B·I(t) + A·ξ(t)·K·X(t)`, which
buys ED2's science plus a semi-analytic accelerated spin-up and traceability diagnostics for free.*

Design-only (no source committed by this document). Companion to
`MEDS_COLUMN_CO2_BALANCE_DESIGN.md` (the fast CAS-CO2 twin this feeds, and the biogeochemistry charter it
amends), `MEDS_PLANT_CARBON_DYNAMICS_DESIGN.md` (the PARTEH allocation whose turnover streams become
litter here), and `MEDS_ENERGY_BALANCE_DESIGN.md` / `MEDS_COLUMN_HYDROLOGY_DESIGN.md` (the soil
temperature and moisture that drive the environmental scalar ξ). ED2 physics reference:
`soil_respiration.f90`, `decomp_coms.f90`, `ed_params.f90`, `structural_growth.f90`, `growth_balive.f90`,
`vegetation_dynamics.f90`. Matrix reference: Luo et al. 2022 (JAMES synthesis), Xia et al. 2012/2013,
Huang et al. 2018, Lu et al. 2020.

---

## 0. Executive summary

**What / where / why (one paragraph).** MEDS today closes only the *fast* half of the carbon cycle: a
prognostic canopy-air CO2 twin advanced from cohort GPP/respiration and a heterotrophic-respiration flux
read off a **frozen, prescribed** soil-carbon scalar (`meds_column_co2.f90`). The slow half is entirely
unbuilt — litter and coarse woody debris from tissue turnover, starvation, mortality and treefall are
silently **dropped**; there are no prognostic soil-C pools, no decomposition transfers, and no nitrogen.
This document specifies the module family that closes it: a new stateless kernel module
**`meds_soil_biogeochem.f90`** in `src/biogeochemistry/` that (i) collects the demography litter/CWD
streams into a **CENTURY-family multi-pool soil-carbon state** (`soil_carbon_t`, expanded from its
one-pool MVP shape), (ii) advances those pools one **daily** step by ED2's decomposition rate laws, and
(iii) returns the respired carbon as the `heterotrophic_respiration` the fast CAS-CO2 kernel already
consumes — making `fast_soil_carbon` a pool that is finally *written* rather than a dead input. The whole
decomposition network — pools, decay rates, inter-pool transfers, temperature/moisture limitation — is
expressed once as the **carbon matrix ODE** `dX/dt = B·I(t) + A·ξ(t)·K·X(t)`; this is not a different
model from ED2's CENTURY, it is ED2's CENTURY *reorganized* so that a semi-analytic accelerated spin-up
and residence-time / storage-capacity traceability diagnostics fall out of the same assembled matrices.

**The five decisive calls:**

1. **The matrix form is the organizing abstraction, ED2 is the physics.** Every ED2 decomposition
   coefficient (decay rate, respired fraction, inter-pool transfer, lignin brake, N immobilization) maps
   onto exactly one entry of `(B, A, ξ, K)` (§2, §3). Nothing in the represented process changes; the
   payoff is that spin-up becomes a linear solve and model diagnosis becomes a decomposition of the same
   matrices (§4). **Scalar placement is `A·ξ·K`, not `ξ·(A·K)`** — because ξ is *per-pool* and `A` is
   non-diagonal, the environmental scalar must multiply each donor's decay *before* the transfer split,
   or mass balance and the capacity formulas break (§3.3, the single most common matrix-form error).

2. **Pools are ED2's CENTURY set, per-patch, above/below split, lignin-tracked.** Seven-pool ceiling —
   metabolic (fast) and structural litter each split above-ground (flammable "grnd") vs below-ground
   ("soil"), plus microbial, slow (humified) and passive SOM (§2.1). Default scheme runs the three-active
   (5-pool) form (microbial+passive inert, `K=0`); the full five-active CENTURY (Bolker 1998 / Koven 2013
   scheme) is a config switch that only turns on two decay rates and the clay-controlled transfer
   off-diagonals. Coarse woody debris routes into structural litter via `f_labile_stem`, faithful to ED2
   (an explicit CWD pool is a reserved additive slot, §2.2).

3. **State lives per-patch in `state/`, threaded through the patch lockstep at P3.** The prognostic
   `soil_carbon_t` (now a pool *vector*) joins `patch_index` as a new per-patch SoA field beside the
   existing fast reservoirs `cas` / `soil_e` / `soil_w` — the same home, the same lockstep obligation:
   `sort_patches` **and** `patch_compact`, `patch_ensure_capacity`, the `fuse_2_patches` area-fraction
   blend, and the disturbance donor-copy (§5.2, the ED2 "forgot to reallocate" trap). Kernels stay
   stateless/`pure`; the pool vector is passed in as a derived-type arg.

4. **Two timescales, one carbon-mass contract, one Rh authority.** The pools are advanced **once per day**
   by the slow `soil_carbon_step`; within a day they are **frozen and read-only** to the fast loop, whose
   CAS-CO2 kernel is the **sole source** of the heterotrophic-respiration flux to the atmosphere — it draws
   an *instantaneous* Rh off the frozen pool with the fast temperature/moisture modifier (exactly ED2's
   exploited fast/slow seam) **and accumulates it sub-daily into a `today_rh` accumulator** (ED2's
   `today_rh`/`today_*_loss`). At the daily slot, `soil_carbon_step` decrements the pools by **that
   accumulated day's flux** (not by an Rh recomputed from daily-mean T/θ), so the pool decrement equals the
   integrated atmospheric flux **by construction** — closing the Jensen gap that recomputing `ξ(⟨T⟩,⟨θ⟩)`
   would leave open. The daily matrix step *never* feeds a second Rh into the CAS twin; it only does
   pool bookkeeping + inter-pool transfers driven by the same accumulated loss. A **daily carbon-mass audit**
   `ΔC_pool = litter_in − rh_today − Δlignin_export` and an explicit **fast-integrated-Rh ≡
   daily-pool-respiration** check guard the boundary — the checks the fast kernel's self-certifying
   `resid ≡ 0` cannot see (§5.4, §8). (NEE = GPP − Ra − Rh closes only its Rh term here; committed
   autotrophic Ra is still 0 in `cohort_co2_flux_t` and NEE stays incomplete until Ra is wired — §9.)

5. **Nitrogen is an optional parallel matrix, deferred but shaped-in.** The same `(A, K)` topology carries
   a companion N pool vector with C:N stoichiometry; N immobilization down-regulates structural
   decomposition through the ξ-family scalar `f_decomp`, and N limitation of NPP is a P2 feedback into the
   plant carbon seam. C-only is the default and the whole MVP (§2.3, §7-P2). PARTEH-H1 (the MEDS allocation
   model) is validated pure-carbon, so nothing upstream needs N to run.

**Honesty scope of the MVP.** P0 ships the pools + Q10 decomposition + the matrix bookkeeping wired to the
*existing* fast Rh, with litter input still a prescribed scalar and pools advanced by a plain daily Euler
step — a structurally complete, carbon-closed slow loop, not yet a validated soil-C product. The matrix
spin-up, DAMM, and C:N land at P1; vertically-resolved pools and N limitation of NPP at P2. Each deferred
piece is named, each carries a test, and none is claimed working before it is built.

---

## 1. Scope, motivation & what exists today

### 1.1 The gap this closes

MEDS's biogeochemistry is a P0, stateless, shared-only sibling library implementing **only** fast carbon
exchange (design `MEDS_COLUMN_CO2_BALANCE_DESIGN.md`). Concretely, in `meds_column_co2.f90`:

- `canopy_air_co2_update` advances the prognostic CAS CO2 twin from `f_bio = Rh + Ra − GPP`;
- `heterotrophic_respiration_flux` computes `Rh = pool · rh_k_base · f_temp(T) · f_water(θ)` off a
  **`soil_carbon_t` that holds one lumped `fast_soil_carbon [kgC/m2]`, prescribed and `intent(in)`** —
  nothing writes it;
- `aggregate_cohort_co2_fluxes` / `column_co2_step` reduce cohort fluxes and assemble the column step.

Everything slow is missing, and the missing pieces are exactly the ecosystem carbon-conservation path:

1. **Litter / CWD input is dropped.** `tissue_turnover_rates` + `plant_carbon_allocation`
   (`meds_plant_carbon_dynamics.f90`) compute `turnover_leaf` / `turnover_fineroot` and the starvation
   `deficit`, and report them out of the kernel — but the driver routes them nowhere. `mortality_step`
   and `apply_patch_disturbance` (treefall) destroy plants and tissue with **no necromass routing**. No
   module collects these streams into litter.
2. **No prognostic soil-C pools, no decomposition transfers.** The slow file `meds_soil_carbon.f90` does
   not exist. `soil_carbon_t` has a single pool; the CENTURY expansion (fast/structural/slow/microbial/
   passive, above/below split, lignin, transfer fractions, respired fractions) is unbuilt. The fast Rh
   draws down a pool that never refills.
3. **No nitrogen.** No N pools, no C:N, no immobilization/mineralization.
4. **No slow-loop plumbing.** No daily carbon-mass audit, no per-layer root-weighted soil T/θ, no
   committed growth/storage respiration (still 0 in `cohort_co2_flux_t`).

The `soil_carbon_t` and `co2_opts_t` types (`meds_biogeochem_types.f90`) were **deliberately shaped** so
these additions are additive — this document fills that reserved space.

### 1.2 Why the matrix form (and not just "port ED2's daily Euler loop")

A direct port of ED2's `update_C_and_N_pools` (a daily forward-Euler advance of seven scalars) would
already close the cycle. We adopt the **matrix organization** on top because it delivers three things the
scalar port cannot, at essentially no extra cost since the coefficients are identical:

- **Semi-analytic accelerated spin-up (SASU).** MEDS has *no meteorological forcing source yet* and the
  soil/passive pools have millennial turnover — a naive bare-ground integration to equilibrium is
  thousands of model-years. The matrix form makes the steady state a single linear solve over the
  **active pool sub-block** `X_ss = −(A·ξ̄·K)^{-1}_active·B·Ī` over climatological-mean drivers (§4.1;
  the full 7×7 solve is *singular* because inert pools have `K=0` — the reduction to active indices is
  mandatory, not cosmetic). One caveat is load-bearing for the current repo state: the climatological
  means `ξ̄`, `Ī` are **not yet producible** from a live run — `Ī` derives from litter, hence GPP and
  vegetation, hence the missing met forcing. So at **P0/P1 SASU runs on *prescribed* climatological means**
  and is validated against a brute-force integration on the *same* prescribed drivers (§4.1, test 5);
  it becomes self-consistent only once a forcing source and a spun-up canopy feed it real means.
- **Traceability diagnostics for free.** Carbon residence time, the storage capacity/potential split, and
  the attribution tree (input × allocation × transfer × turnover × environment) are read directly off the
  assembled `(B, A, ξ, K)` (§4.2-4.3) — the standard tool for diagnosing why two runs (or MEDS vs ED2)
  store different amounts of carbon.
- **A clean, GPU-friendly, per-patch-parallel kernel.** The assembly of `(A, ξ, K)` and the frozen-matrix
  daily advance are stateless pure arithmetic over small fixed-size arrays, embarrassingly parallel across
  patches — consistent with MEDS's OpenMP-target posture and the stateless-kernel discipline of every
  sibling module.

The matrix form is **exactly** ED2's CENTURY: `update_C_and_N_pools` already separates
`loss = decay_rate·scalar·pool` (that is `ξ·K·X`), the respired fraction `er_i`, and the transfer
fractions `ex_i_j` (that is `A`) — it is a matrix model written out element-by-element. We give it its
matrix name so the diagnostics and the spin-up become available.

### 1.3 Reachability facts from the ED2 reference

- ED2 pools are **per-patch scalars** (`sitetype csite`, indexed `ipa`), *not* vertically resolved; only
  the climate scalars `A_decomp`/`B_decomp` use soil-layer T/θ down to `rh_active_depth`
  (`soil_respiration.f90`). MEDS mirrors this at MVP (per-patch pools; vertical resolution is P2).
- ED2 runs decomposition on **two timescales**: sub-daily flux accumulation in `soil_respiration_driver`
  (every `DTLSM`, accumulating `today_rh` and `today_*_loss`), and a **once-per-day** Euler pool advance
  in `update_C_and_N_pools`, called from `vegetation_dynamics` after phenology + allocation, before
  disturbance, with a fatal carbon-budget check (`check_budget_soilc`). MEDS adopts this exact cadence
  (§5.3).
- ED2's default `DECOMP_SCHEME` (0-4) collapses to **three active pools** (microbial + passive decay rates
  are literally 0); the full five-active CENTURY is `DECOMP_SCHEME=5` (Bolker/Koven). MEDS exposes the
  same enum (§6).
- ED2 has **no explicit coarse-woody-debris pool** — dead wood (`bdead`/`bsapwood`/`bbark`) enters
  structural litter through the `f_labile_stem(ipft)` split. MEDS follows this (§2.2).
- N is **prognostic but minimal** in ED2 (fast + structural + mineralized pools only; microbial/slow/
  passive carry no N; no nitrification/denitrification/gaseous loss/deposition/fixation; limitation only
  via `f_decomp` when `N_DECOMP_LIM=1`, default off). MEDS's optional N cycle inherits this scope (§2.3).

---

## 2. ED2-faithful pool structure

### 2.1 The pool vector

The prognostic state is a per-patch carbon vector `X [kgC/m2]`, ordered litter → SOM → passive so that the
**default-scheme active block of `A` is lower-triangular** (back-substitution solvable, §4.1); scheme 5's
microbial⇄slow / slow⇄passive back-transfers make its active block non-triangular (a small dense solve,
§3.2/§4.1) — and the *full* matrix is singular in every scheme (`K=0` inert pools), so the solve is always
on the active sub-block, never the full 7×7. The
seven-pool ceiling matches ED2 exactly; the above-/below-ground split (ED2 "grnd" flammable vs "soil") is
carried so fire (a later disturbance class) can consume the above-ground fraction:

| idx | pool | field | `K` [1/yr] (repr.) | receives litter (B) | env scalar ξ |
|---|---|---|---|---|---|
| 1 | metabolic litter, above | `fast_grnd_carbon` | `k_fast ≈ 11` | ✔ labile·above | `A_decomp` (surface) |
| 2 | metabolic litter, below | `fast_soil_carbon` | `k_fast ≈ 11` | ✔ labile·below | `B_decomp` (active layer) |
| 3 | structural litter, above | `struct_grnd_carbon` | `k_struct ≈ 4.5` | ✔ structural·above + CWD·above | `A_decomp·L_g·f_decomp` |
| 4 | structural litter, below | `struct_soil_carbon` | `k_struct ≈ 4.5` | ✔ structural·below + CWD·below | `B_decomp·L_s·f_decomp` |
| 5 | microbial SOM | `microbial_carbon` | `k_micr` (0 default; ≈6 sch.5) | ✗ (transfers only) | `B_decomp` |
| 6 | slow / humified SOM | `slow_carbon` | `k_slow ≈ 0.2` | ✗ | `B_decomp` |
| 7 | passive SOM | `passive_carbon` | `k_pass` (0 default; ≈0.012 sch.5) | ✗ | `B_decomp` |

Structural pools additionally carry a **lignin sub-state** `struct_grnd_lignin` / `struct_soil_lignin`
`[kgC/m2]` (ED2 `structural_grnd_L`) so the lignin fraction `f_lignin = L/C` can brake decomposition and
weight the respired fraction (§2.4). **Lignin decrement law (ED2-faithful, made explicit here).** ED2
removes structural lignin in **strict proportion** to the structural carbon it decays — the same
first-order rate acts on `L` as on `C`, so `L` and `C` decay together and `f_lignin = L/C` is
**conserved through pure decomposition** (ED2 does *not* preferentially spare lignin within a step; the
lignin brake instead makes high-`f_lignin` litter decay slower *as a whole*). Concretely, per structural
pool `s`: `L_s ← L_s·(1 − ξ_s·K_s·dt)` in lockstep with `C_s ← C_s + (u_s − ξ_s·K_s·C_s)·dt`, and
incoming litter adds `lignin_in` (§2.2) so `f_lignin` drifts only via *inputs* of differing chemistry,
never via decay. Lignin is a **passive tracer bounded `0 ≤ L_s ≤ C_s`**; it is outside the
`n_soil_pool` carbon matrix (it carries no independent mass — its carbon is already counted in
`struct_*_carbon`) and therefore is **not** double-counted in the `soilc_audit_t` carbon balance, but it
gets its **own** bounds check `0 ≤ L_s ≤ C_s` and its own balance `ΔL_s = lignin_in_s − ξ_s·K_s·L_s·dt`
(test 7). A preferential-loss variant (lignin enriches as non-lignin decays faster) is a reserved config
switch, not the ED2 default. The decay rates are the ED2 `ed_params.f90` values (converted
`/yr → /day` by `yr_day` at the seam); the representative numbers above are for orientation —
**the exact per-scheme values are required from TOML and verified against `ed_params.f90` at load** (§6).
`k_slow` in particular is scheme-dependent (the ED2 memory-brief's "100.2" concatenation must be resolved
to the `ed_params.f90` value before wiring; treat as a load-time provenance check).

`soil_carbon_t` grows from its one-pool MVP shape (which stays index 2, `fast_soil_carbon`, so the fast Rh
keeps compiling) to the full vector additively:

```fortran
!----- Slow, stateful per-patch soil-carbon pools (written DAILY, read-only in the fast loop). --------!
!      Expanded from the P0 one-pool shape; fast_soil_carbon KEEPS its name/role so meds_column_co2      !
!      compiles unchanged. The active-pool count is set by decomp_scheme (default 3-active: idx 5,7 = 0). !
type :: soil_carbon_t
   ! carbon pools [kgC/m2] -- the matrix state vector X, ordered litter -> SOM -> passive
   real(wp) :: fast_grnd_carbon    = 0.0_wp   !< X(1) metabolic litter, above-ground (flammable)
   real(wp) :: fast_soil_carbon    = 0.0_wp   !< X(2) metabolic litter, below-ground  (the P0 pool)
   real(wp) :: struct_grnd_carbon  = 0.0_wp   !< X(3) structural litter + CWD, above
   real(wp) :: struct_soil_carbon  = 0.0_wp   !< X(4) structural litter + CWD, below
   real(wp) :: microbial_carbon    = 0.0_wp   !< X(5) microbial SOM        (scheme-5 only)
   real(wp) :: slow_carbon         = 0.0_wp   !< X(6) slow / humified SOM
   real(wp) :: passive_carbon      = 0.0_wp   !< X(7) passive SOM          (scheme-5 only)
   ! lignin sub-state of the structural pools [kgC/m2] (fraction f_lignin = L/C brakes decomposition)
   real(wp) :: struct_grnd_lignin  = 0.0_wp
   real(wp) :: struct_soil_lignin  = 0.0_wp
   ! optional nitrogen twin [kgN/m2] -- present only when opts%n_cycle_on (P1/P2); C-only default
   real(wp) :: fast_grnd_n = 0.0_wp, fast_soil_n = 0.0_wp
   real(wp) :: struct_grnd_n = 0.0_wp, struct_soil_n = 0.0_wp
   real(wp) :: mineralized_n = 0.0_wp
end type soil_carbon_t
```

A parallel fixed-size parameter `n_soil_pool = 7_ik` and index parameters (`IP_FAST_GRND=1`, …,
`IP_PASSIVE=7`) live in `meds_biogeochem_types`; a pair of thin `pack_pool_vector` / `unpack_pool_vector`
helpers marshal between the readable named type and the length-`n_soil_pool` array the matrix kernels
operate on, so the kernels stay array-based (GPU-friendly) while state stays self-documenting.

### 2.2 Litter and coarse woody debris inputs — the demography seam

The input vector `u = B·I [kgC/m2/day]` is assembled from the three currently-dropped demography streams,
split labile-vs-structural by `f_labile_leaf(ipft)` / `f_labile_stem(ipft)` and above-vs-below by
`agf_bs(ipft)`, exactly as ED2's `update_litter_inputs` / structural-growth litter staging:

- **Tissue turnover + phenological drop** — `turnover_leaf`, `turnover_fineroot` from
  `plant_carbon_allocation` (`carbon_npp_t`), currently reported and dropped. Leaf turnover splits
  labile/structural and above (`agf`) / below; fine-root turnover is below-ground.
- **Starvation / negative-NPP deficit** — the `deficit` the plant kernel reports when storage is
  exhausted; the shed leaf/fineroot carbon becomes litter.
- **Mortality + treefall necromass** — the leaf/fineroot/wood/storage carbon of killed plants.
  `mortality_step` reduces `nplant` by `Δn = nplant·(1−exp(−mort·dt))`; the lost biomass
  `Δn·(leaf+fineroot+wood+nonstructural)` is necromass. Treefall (`apply_patch_disturbance`) kills the
  canopy; its biomass is necromass on the source patch. **Coarse woody debris = the wood component of this
  necromass**, routed to the structural pools via `f_labile_stem` (mostly structural, high lignin) — ED2's
  no-separate-CWD-pool convention. An explicit CWD pool (a slower-turnover structural sub-pool, index 8) is
  a reserved additive slot for fire/fragmentation fidelity; MVP folds CWD into structural.

The lignin input to the structural pools is `struct_litter · (l2n_stem/c2n_stem)`-weighted (ED2), so
`f_lignin` of the incoming structural litter is set by the stem lignin:C ratio.

**Ownership and the library wall (critical).** Litter production is a *demography/plant* quantity;
decomposition is *biogeochemistry*. Neither library may name the other's types. The streams therefore
cross the wall as **plain data arrays**, assembled by the driver (`meds_aux`, the one layer above both) —
the exact mirror of `update_demography(growth[], npp[], mortality[], …)`. The new seam is a per-patch
litter-input record built by the driver from the demography turnover/mortality/treefall outputs and handed
to the biogeochem kernel as `intent(in)`. **Demography gains no biogeochem edge and biogeochem gains no
demography edge** (§5.1). Concretely, `apply_carbon_npp` / `mortality_step` / `apply_patch_disturbance`
must additionally *emit* their necromass into a per-patch accumulator (a new `carbon_flux_block`-style
plain array on the state hub, zeroed daily) that the driver reads — this is the demography-side change
(pure data, no plant/biogeochem dependency).

### 2.3 Optional nitrogen (deferred, shaped-in)

Faithful to ED2's minimal N: prognostic N only in the fast + structural + mineralized pools; microbial/
slow/passive carry no N mass. C:N ratios (`decomp_coms`): `c2n_structural = 150`, `c2n_slow = 10` (used
for microbial+slow+passive), `c2n_fast` initialized `≈ 30`, per-PFT `c2n_leaf` (SLA-dependent, ~86 for
tropical), `c2n_stem = 150`. N enters through the same litter seam (C flux ÷ tissue C:N). N immobilization
down-regulates **structural** decomposition via `f_decomp = N_supply/(N_demand + N_supply)` (§2.4), the
one ξ-family scalar that couples the C and N matrices. No nitrification/denitrification/gaseous loss/
deposition/fixation. **N limitation of NPP** (feeding `f_decomp`-analogous down-regulation back into the
plant carbon seam) is the P2 feedback. Default is C-only (`n_cycle_on = .false.`), for which every N field
and `f_decomp` collapse to inert (`f_decomp ≡ 1`).

### 2.4 The ED2 rate laws (the exact physics the matrix encodes)

Per pool `i`, the daily carbon decomposition loss is `loss_i = ξ_i · K_i · X_i` (units `kgC/m2/day` with
`K` in `/day`). The pieces, straight from `soil_respiration.f90` / `decomp_coms.f90`:

**Environmental scalar `ξ` = temperature × moisture × oxygen** (ED2 `het_resp_weight`,
`rh_weight = temperature_scale · water_limitation · oxygen_limitation`), evaluated on the **surface** layer
for above-ground ("grnd") pools (`A_decomp`) and over the **active layer** for below-ground pools
(`B_decomp`). Six selectable functional forms by `decomp_scheme`; MEDS ships the two that matter for the
MVP and reserves the rest:

- **Temperature.** Scheme 0 (default): `min(1, exp(resp_temp_increase·(T − 318.15)))`,
  `resp_temp_increase = 0.0757` — ED2's capped exponential (identical to the fast Rh's `HR_EXP_ED2`).
  Scheme 5: Collatz **Q10** `rh0·rh_q10^((T−288.15)/10)`, `rh0=0.700`, `rh_q10=1.500`, ref 15 °C
  (identical to the fast Rh's `HR_Q10`). The matched functional form is deliberate — but note the
  reconciliation is **not** "the daily ξ is the time-average of the fast `f_temp`": because ξ is the
  *convex product* `f_temp(T(t))·f_water(θ(t))`, Jensen's inequality gives
  `⟨f_temp·f_water⟩ ≠ f_temp(⟨T⟩)·f_water(⟨θ⟩)` under any diurnal cycle, so recomputing ξ from
  daily-mean state would leave the pool decrement and the integrated atmospheric flux **mismatched**.
  MEDS instead follows ED2 exactly: the fast loop **accumulates the true sub-daily flux** into `today_rh`
  and the daily step decrements the pool by *that* integral (§5.4). Matched chemistry then guarantees the
  fast Rh (flux to atmosphere) and the daily pool loss are the **same number**, not merely close.
  Scheme 5's *uncapped* Q10 exceeds 1 above 15 °C (≈2.4 at 45 °C), so the correct invariant is
  **`ξ > 0`**, not `ξ ∈ (0,1]`; only scheme 0's `min(1,·)`-capped temperature term is bounded by 1
  (the capacity/residence formulas need only `ξ > 0`, §3.1, §4).
- **Moisture.** Scheme 0: one-sided about `resp_opt_water = 0.8938` — dry side
  `exp((rel_m − resp_opt_water)·5.0786)`, wet side 1. Scheme 2 (CENTURY): logistic
  `1/(1+exp(0.24·(18°C − T)))·1/(1+exp(0.60·(T − 45°C)))` for temperature and a paired dry/wet logistic
  for moisture. `rel_m = (θ − θ_dry)/(θ_sat − θ_dry)`.
- **Oxygen (anoxia).** Scheme 0: wet-side `exp((resp_opt_water − rel_m)·4.5139)`. Scheme 5:
  `rel_oxygen^0.6` off air-filled porosity.

**Air-dry floor `θ_dry` is a correctness dependency, not a free choice** (inherited caveat from the fast
Rh): it must resolve to a real `soilcp`/retention-curve air-dry moisture, never silently to the van
Genuchten residual `θ_res` (§9). The DAMM option (already in `meds_column_co2`) *removes* this floor
problem via mechanistic dual-substrate diffusion and is the preferred P1 form (§6).

**Lignin brake (structural pools only).** `f_lignin = L/C` (=1 if C=0);
`L_g = exp(−e_lignin·f_lignin)`, `e_lignin = 3.0` slows structural decay; `L_s` analogous. The structural
respired fraction is lignin-weighted: `er_struct = (1−f_lignin)·r_struct_o + f_lignin·r_struct_l`.

**N immobilization brake (structural pools, N cycle on).** `f_decomp = N_supply/(N_demand + N_supply)`
with `N_demand = ξ·K·((1−er_struct)/c2n_slow − 1/c2n_structural)·X_struct` and
`N_supply = N_immobil_supply_scale·mineralized_n`. Multiplies the structural ξ. `≡ 1` when N is off.

**Respired fraction and transfer split.** Each pool `j` sends a fraction `er_j` of its loss to CO2
(heterotrophic respiration) and the complement `(1 − er_j)` to downstream pools via the transfer fractions
`ex_j→i`. Default (scheme 0-4): `er_fast = 1.0` (metabolic litter respires 100%, transfers nothing);
`er_struct = 0.3` (30% respired, 70% → slow); `er_slow = 1.0` (slow is the sole SOM outlet, respires
100%); microbial/passive inert. Scheme 5: `er_fast = 0.55`, `er_struct ∈ {0.20 lignin, 0.50 non-lignin}`,
`er_micr = 0.60 + 0.17·xsand` (sand-controlled), `er_slow = er_pass = 0.55`, with clay-controlled
transfers to passive `ex_·→pass = min(1−er, 0.003 + slope·xclay)`. These `er` and `ex` are the entries of
`A` (§3.2).

---

## 3. The matrix formulation

### 3.1 The canonical equation

The per-patch soil-carbon column obeys the standard matrix ODE (Luo et al. 2022):

```
dX/dt = B·I(t)  +  A·ξ(t)·K·X(t)                                 [kgC/m2/day]
```

element-wise (the mass balance the matrix encodes, and the exact form MEDS integrates):

```
dX_i/dt = u_i(t)  +  Σ_{j≠i} a_ij·ξ_j·K_j·X_j  −  ξ_i·K_i·X_i
```

with

- **`X`** — the length-`n_soil_pool` carbon vector `[kgC/m2]` (§2.1);
- **`u = B·I`** — the litter/CWD input vector `[kgC/m2/day]` (§2.2); `B` is the allocation/partition
  operator, nonzero only in the litter pools (1-4), `I` the total litter carbon; MEDS assembles `u`
  directly (multiple chemically-distinct streams), which is `B·I` in canonical form;
- **`K`** — diagonal baseline (potential) decay rates `[1/day]` (§2.1), strictly positive on active
  pools;
- **`ξ(t)`** — diagonal per-pool environmental scalar (temperature × moisture × oxygen, with the lignin
  `L` and N `f_decomp` brakes folded into the structural entries) `[-]`, **`> 0`** (only the capped
  scheme-0 form is additionally `≤ 1`; scheme 5's Q10 is not — §2.4). All the matrix algebra below needs
  only `ξ > 0`;
- **`A`** — the donor/transfer matrix `[-]`: `a_jj = −1` **uniformly** for every pool `j` (all of a
  donor's decayed carbon leaves it), off-diagonal `a_ij ∈ [0,1]` is the fraction of pool `j`'s loss
  entering pool `i`; the **respired fraction** of `j` is `er_j = 1 − Σ_{i≠j} a_ij` (the CO2 that leaves
  the network entirely). The convention is uniform on purpose: for an **inert** pool (`K_j=0`) the entry
  `a_jj=−1` is *harmless* — `ξ_j·K_j·X_j = 0` zeroes its whole row/column contribution — and it keeps the
  column-sum bookkeeping `Σ_i a_ij = −er_j` (test 2) valid for every pool. Inert pools are simply
  **excluded from the active sub-block** the solver operates on (§4.1); they are never given `a_jj=0`
  (which would falsely read as `er_j=1` in the column check).

**Heterotrophic respiration** — the flux that closes NEE — is the respired complement, the carbon that
`A` does *not* transfer:

```
Rh(t) = Σ_j er_j · ξ_j(t) · K_j · X_j(t)  =  −1ᵀ·A·ξ(t)·K·X(t)    [kgC/m2/day]
```

(the column sums of `A` are `−er_j`, so `−1ᵀA` picks out the respired fractions). This **instantaneous**
form is what the **fast loop** evaluates each substep off the frozen pool `X` — converted to `[umol/m2/s]`
by `kgCday_2_umols` and handed to `canopy_air_co2_update` as `heterotrophic_respiration` — and it is the
**single** Rh authority into the CAS-CO2 twin. The daily `soil_carbon_step` **accumulates** the same flux
(`today_rh`, via `xi_int`) and uses it *only* to decrement the pools and drive inter-pool transfers; it
does **not** hand a second, separately-computed Rh to CAS (§5.4 — resolving the double-count hazard).

### 3.2 The concrete `A` and `K` for the MEDS pool set

Ordering `X = [fast_grnd, fast_soil, struct_grnd, struct_soil, microbial, slow, passive]ᵀ`.

**Default (`decomp_scheme` 0-4, three-active).** Microbial (5) and passive (7) are inert (`K_5=K_7=0`),
so their rows/columns are null; metabolic litter respires fully; structural sends 70% to slow; slow
respires fully:

```
        col1    col2    col3        col4        col5  col6   col7
       fast_g  fast_s  struct_g    struct_s    micr  slow   pass
row1  [ -1      0       0           0            0     0      0   ]   fast_grnd
row2  [  0     -1       0           0            0     0      0   ]   fast_soil
row3  [  0      0      -1           0            0     0      0   ]   struct_grnd
row4  [  0      0       0          -1            0     0      0   ]   struct_soil
row5  [  0      0       0           0           -1     0      0   ]   microbial (INERT, K_5=0 -> excluded from active solve)
row6  [  0      0      (1−er_stg)  (1−er_sts)    0    -1      0   ]   slow      <- 0.7 from each structural
row7  [  0      0       0           0            0     0     -1   ]   passive   (INERT, K_7=0 -> excluded from active solve)

K = diag(k_fast, k_fast, k_struct, k_struct, 0, k_slow, 0)          [1/day]
er = (1.0, 1.0, er_stg, er_sts, 1.0*, 1.0, 1.0*)  ;  er_stg=er_sts=0.3 (lignin-weighted, §2.4)
                                                  (* inert: a_jj=−1 kept uniformly (§3.1), but K_j=0 zeros
                                                     the row and the pool is dropped from the active solve, §4.1)
```

`Rh = ξ1·k_fast·X1 + ξ2·k_fast·X2 + er_stg·ξ3·k_struct·X3 + er_sts·ξ4·k_struct·X4 + ξ6·k_slow·X6`.

**Full CENTURY (`decomp_scheme = 5`, five-active, Bolker/Koven).** Microbial and passive activate; fast
litter → microbial; structural non-lignin → microbial, lignin → slow; microbial ⇄ slow with a
clay-controlled tail to passive; slow ⇄ passive:

```
row5 (microbial):  a_5,1=(1−er_f), a_5,2=(1−er_f), a_5,3=(1−er_stg)(1−f_lig), a_5,4=(1−er_sts)(1−f_lig),
                   a_5,6 = (1−er_slow) − ex_slow→pass                     (slow's non-passive tail returns to micr)
row6 (slow):       a_6,3=(1−er_stg)·f_lig, a_6,4=(1−er_sts)·f_lig,
                   a_6,5 = (1−er_micr) − ex_micr→pass,  a_6,7 = (1−er_pass)
row7 (passive):    a_7,5 = ex_micr→pass = min(1−er_micr, 0.003+0.032·xclay),
                   a_7,6 = ex_slow→pass = min(1−er_slow, 0.003+0.009·xclay)

K = diag(12, 12, 1.5, 1.5, 6.0, 0.15, 0.012)/yr   (repr.; from ed_params.f90, verify at load)
```

The transfer topology is **texture-controlled** exactly as CENTURY: sand controls the microbial respired
fraction `er_micr`, clay controls the fractions routed to the passive pool. `A` is lower-triangular in the
litter → microbial → slow → passive order **only up to the microbial⇄slow (`a_5,6`, `a_6,5`) and
slow⇄passive (`a_6,7`, `a_7,6`) back-transfers, which place nonzero entries *above* the diagonal — so in
scheme 5 `A` is genuinely non-triangular** and `det A ≠ ±1`. Pure back-substitution therefore does **not**
apply to scheme 5; its (≤5×5) active block is solved by a small dense LU / Gaussian elimination with
partial pivoting (§4.1). The system is still non-singular — the active-block diagonal `−ξ·K` dominates the
weak off-diagonal transfers — but invertibility is an argument about the *active block*, not a triangular
determinant.

### 3.3 Scalar placement — the load-bearing correctness rule

The environmental scalar sits **inside** the transfer, `A·ξ·K`, **not outside** as `ξ·(A·K)`. Because ξ is
per-pool (surface `A_decomp` for grnd pools, `B_decomp` for below; the structural entries additionally
carry `L` and `f_decomp`) and `A` is non-diagonal, the two orderings differ. The physically correct
statement is: **each donor `j` decomposes at its own rate `ξ_j·K_j`, and only then is the resulting loss
split** (fraction `er_j` respired, fractions `a_ij` transferred). Writing `ξ·A·K` would incorrectly apply
the *receiving* pool's scalar to a donor's transferred carbon. Since ξ and K are diagonal they commute
(`A·ξ·K = A·K·ξ`), but neither commutes with `A`. Getting this wrong silently breaks mass balance and
every capacity/residence-time formula (§4) — it is the single most common error in matrix-form ports and
is called out here so the assembly kernel is written `A·(ξ·K·X)` from the start.

### 3.4 Input convention

`I` is **litter carbon** (post-mortem / post-turnover biomass), so `B·I` enters the pools directly with no
respiration term — the autotrophic respiration that would appear if `I` were GPP has already been removed
upstream in the plant carbon budget (`net_carbon = GPP − Ra`). Do not mix conventions: the biogeochem
matrix input is litter, never GPP.

---

## 4. Semi-analytic spin-up & traceability diagnostics

### 4.1 Accelerated spin-up (SASU)

The steady state under climatological-mean drivers is a direct linear solve, not a multi-millennium
integration. Setting `dX/dt = 0` with time-averaged input `Ī` and scalar `ξ̄`:

```
0 = B·Ī + A·ξ̄·K·X_ss    ⟹    X_ss = −(A·ξ̄·K)^{-1}·B·Ī      [kgC/m2]   (over ACTIVE indices only)
```

The minus sign and the use of **time-averaged fluxes** (arithmetic mean of `ξ` and `I` over the forcing
cycle — *not* averaged states) are both essential.

**The solve runs on the ACTIVE sub-block only — the full `n_soil_pool` system is singular.** In the
default scheme `K_5=K_7=0`, so columns 5 and 7 of `A·ξ̄·K` are identically zero: `det(A·ξ̄·K)=0`, the
inverse **does not exist**, and a blind triangular back-substitution would **divide by zero** on the inert
diagonal. This is not a corner case — it is the *default* scheme. So the solver first selects the active
index set `𝒜 = { j : K_j > 0 }` (driven by `decomp_scheme`: `{1,2,3,4,6}` for schemes 0-4;
`{1,2,3,4,5,6,7}` for scheme 5), forms the reduced square block `M = (A·ξ̄·K)[𝒜,𝒜]` and RHS
`b = −(B·Ī)[𝒜]`, solves `M·X_ss[𝒜] = b`, and **holds every inert pool at its conserved value** (0 —
it has no litter input and no transfer *into* it in the default scheme, so its equilibrium is 0). All
residence-time and capacity diagnostics (§4.2-4.3) operate on the **same** active block; the full-matrix
inverse is never formed or referenced.

**How the active block is solved — scheme-dependent.**
- **Default (schemes 0-4), 3-active + grnd/below = the `{1,2,3,4,6}` block:** with `a_jj=−1` and no
  above-diagonal entries, `M` *is* lower-triangular after the litter → slow ordering, so
  **back-substitution** solves it (`O(|𝒜|²)`, a few dozen flops, no LAPACK, GPU-safe).
- **Scheme 5, 5-active block:** the microbial⇄slow and slow⇄passive back-transfers make `M`
  **non-triangular** (§3.2), so back-substitution is invalid — MEDS uses a small **dense LU / Gaussian
  elimination with partial pivoting** on the `≤5×5` block (still allocation-free, LAPACK-free,
  branch-light, GPU-safe). Diagonal dominance of `−ξ̄·K` over the weak transfers guarantees a
  well-conditioned pivot. A single `solve_active_block(M, b)` helper dispatches triangular vs dense on the
  block's structure, so the caller is oblivious.

**Iterated quasi-linear fixed point (Xia 2012 SASU).** Because ξ (moisture/temperature responses) and, in
scheme 5, some transfers are mildly state- or driver-dependent, spin-up is an outer fixed-point loop, not
a one-shot inverse: (1) run the model one climatological cycle to refresh `ξ̄`, `Ī` (and, with N on, the
C:N-dependent `A`, `f_decomp`); (2) solve the analytic `X_ss` on the active block; (3) reset the pools;
(4) repeat. A handful of outer iterations replaces thousands of bare-ground years. For the C-only,
fixed-ξ̄ MVP the inner solve is a single active-block back-substitution.

**Two honest caveats on SASU's value at the current repo state.** (a) **Prescribed means, not live means.**
SASU consumes climatological-mean litter `Ī` and scalar `ξ̄`, but `Ī` derives from GPP/vegetation, which
needs the still-absent met forcing — so at **P0/P1 the means are *prescribed*** (a supplied litter/scalar
climatology), and SASU is validated by agreeing with a brute-force integration on the *same* prescribed
drivers (test 5), not yet by matching a live spun-up canopy. It becomes fully self-consistent only when a
forcing source closes the litter→ξ̄/Ī loop. (b) **The mean-driver equilibrium is an approximation to the
true periodic attractor.** Replacing `ξ(t),I(t)` by their means discards the `⟨ξ·X⟩ − ξ̄·X̄` covariance of
periodic forcing (the fundamental Xia-2012 bias); the outer fixed-point loop **reduces but does not
eliminate** it. Scheme-5/N state-dependence is an *additional* source of iteration, not the only one — the
covariance bias is present even in the linear, fixed-topology default.

### 4.2 Carbon residence time

```
τ    = −(A·ξ̄·K)^{-1}·B        (per-pool residence-time contributions, X_ss = τ·Ī)          [active block]
τ_E  = −1ᵀ·(A·ξ̄·K)^{-1}·B     (scalar ecosystem soil-carbon residence time)
τ'_E = −1ᵀ·(A·K)^{-1}·B        (baseline, climate-removed)   ⟹   τ_E = τ'_E / ξ̄   (scalar-ξ limit)
```

**Unit convention (load-bearing — the kernel runs `K` in `/day`).** Since `K` in the kernel is `[1/day]`
(§2.4, §3.1), `−(A·ξ̄·K)^{-1}` comes out in **days**; to report residence times in **years** either
(i) assemble the diagnostic with the *pre-conversion* `/yr` `K` (recommended — `K` is loaded in `/yr` and
only converted at the step seam), or (ii) multiply the `/day` result by `yr_day`. The `[yr]` labels above
assume convention (i); the diagnostic routine documents which it used. Both `τ` solves reuse the **active
sub-block** machinery of §4.1 (inert pools carry no residence time — they are dropped, not inverted).

### 4.3 Storage capacity, potential, and the traceability tree

The transient state chases a moving attractor (Luo 2017):

```
X_c(t) = −(A·ξ(t)·K)^{-1}·B·I(t)     [kgC/m2]  instantaneous carbon storage CAPACITY (the equilibrium the
                                                system chases under current climate + litter input)
X_p(t) = X_c(t) − X(t)               [kgC/m2]  storage POTENTIAL (remaining sink; <0 = a source)
dX/dt  = A·ξ·K·(X − X_c) = −A·ξ·K·X_p           the transient relaxes toward X_c at rate set by A·ξ·K
```

and the diagnostic attribution `X_c = I · τ'_E · (1/ξ_temp) · (1/ξ_water)` decomposes any MEDS-vs-ED2 (or
scenario) difference in stored soil carbon hierarchically into **input** (litter production), **baseline
residence time** (allocation `B`, transfer network `A`, turnover `K`), and **environmental scalars**
(temperature, moisture, and — N on — nutrient `1/ξ_N`). These are **pure post-processing** off the
already-assembled matrices (§5.1), no new state — the standard tool for diagnosing the soil-C module
against ED2 (§8).

---

## 5. MEDS module design

### 5.1 Library-DAG placement

`biogeochemistry` is a **shared-only sibling** of `biophysics` and `plant` (charter already widened by the
CO2 design §2.4 to own *both* fast CO2 exchange and slow carbon pools). The new slow module keeps that
wall: it links `meds_shared` **only** — no `site_t`, no `meds_demography_types`, no biophysics types.
Compute is stateless; the pool vector and the litter-input vector arrive as derived-type / plain-array
args; the driver (`meds_aux`) weaves demography ⟷ biogeochem.

```
shared ─┬─ allometry ─ state ─ demography ─┐
        ├─ plant  (carbon allocation: emits turnover/deficit litter as DATA)   │
        ├─ biophysics (soil T, θ  ->  drive ξ)                                  ├─ aux (driver:
        └─ biogeochemistry ──────────────────────────────────────────────────┘    collect litter[] from
             meds_biogeochem_types   (extend: soil_carbon_t vector, decomp_opts_t, litter_input_t)          demography/plant,
             meds_column_co2         (REUSED: fast CAS-CO2 twin + fast Rh read off the pools)                soil T/θ from
             meds_soil_biogeochem    (NEW: assemble A/ξ/K, build litter u, daily matrix step, spin-up,       biophysics, call
                                      residence/capacity diagnostics) -- links shared ONLY                    soil_carbon_step,
                                                                                                             feed Rh to column_co2_step)
```

No illegal edges: demography emits litter as plain data (no biogeochem dependency); biogeochem consumes
plain data (no demography/plant/biophysics dependency); the driver is the only layer that sees both.

### 5.2 Where the pool state lives — per-patch `state/`, justified

**Decision: `soil_carbon_t` becomes a per-patch SoA field on `patch_index` (`state/`), beside the existing
fast reservoirs `cas` / `soil_e` / `soil_w`.** ED2 holds the pools as per-patch scalars (`csite`), and
MEDS already has the exact home and pattern: `patch_index` carries `cas(ip)`, `soil_e(ip)`, `soil_w(ip)`
as per-patch reservoir arrays that thread through fusion/disturbance. `soil_carbon_t(ip)` joins them.

The rejected alternative — a free-standing "column-state" type outside `state/` — is worse: it would
duplicate the patch-lockstep machinery and break the single-source-of-truth invariant that *one* place
knows every per-patch field.

**Lockstep obligation (the ED2 "forgot to reallocate" trap).** Patches have **no single reorder routine**,
so a new per-patch field must be threaded through *every* permute/pack/blend site (CLAUDE.md; identical to
the P3 obligation the CO2 twin's `can_co2` carries):

- `patch_alloc` / `patch_ensure_capacity` (`meds_demography_types`) — allocate the new array;
- `sort_patches` **and** `patch_compact` (`meds_demography_fusefiss`) — the two permute/pack sites;
- `fuse_2_patches` — **area-fraction blend**: pools are extensive per unit ground area, so the fused pool
  is the area-weighted mean `(a₁·X₁ + a₂·X₂)/(a₁+a₂)` (conserves patch-total carbon, like `recruit_pool`);
- `apply_patch_disturbance` donor-copy — the age-0 gap **inherits area-weighted donor pools** (the
  gap's soil carbon is the mix of its donors', exactly as it inherits `cas`/`soil_e`/`soil_w`).

Miss one and the pools silently desync on the next fusion — the class of bug the centralized lockstep
exists to prevent. This is P3 work, landing together with the CO2 twin's patch-state lockstep (§7).

**Restart/checkpoint I/O — the expanded state MUST be serialized (P3).** `soil_carbon_t` grows from one
scalar to **7 carbon pools + 2 lignin sub-states + 5 N fields**, all prognostic per-patch. `meds_io`'s
state stream (`io_write_state` → `-S-<date>.nc`, `io_read_state`) currently serializes only the demographic
state and the fast reservoirs; **if the new pools are not added to the netCDF state schema, a
checkpoint/restart silently zeroes all soil carbon (and lignin, and N)** — a soil-C spin-up erased on the
first restart. So the P3 lockstep list extends to **`meds_io`**: add the 7+2+5 fields to
`io_write_state`/`io_read_state` and the state schema, beside `can_co2`/`soil_e`/`soil_w`. As with the rest
of the state stream, only the *prognostic* pools + lignin + N are written; the assembled matrices and
diagnostics are re-derived on read.

**Cold-start initial values (a fresh, pre-spin-up bare-ground run).** On a cold start every carbon pool,
lignin sub-state, and N field is initialized to **0** (a bare mineral surface with no accumulated litter
or SOM); `f_lignin` is then undefined and read as 0 (the `C=0 ⟹ f_lignin=0` guard, §2.4). A cold run is
therefore *defined* but *un-spun-up*; the intended path to a realistic initial state is **SASU (§4.1)** —
solve `X_ss` on prescribed climatological means and stamp the pools with it — not a literal
multi-millennium bare-ground integration. `init_bare_ground` sets the zeros; a `[soil_carbon.init]` option
may instead prescribe non-zero starting pools (e.g. an observed SOC profile) or request an
`initialize_from_steady_state` SASU stamp at init.

### 5.3 The kernel module `meds_soil_biogeochem.f90`

Stateless, `pure` where the arithmetic allows; all state passed as args. Public procedures:

```fortran
!---------------------------------------------------------------------------------------------------!
! Assemble the diagonal INSTANTANEOUS environmental scalar xi at the given soil T/theta (+ lignin, N   !
! brakes). Per-pool: A_decomp (surface) for grnd pools, B_decomp (active layer) for below; scheme-      !
! selectable. This is the SAME functional form the fast CAS-CO2 kernel uses per substep; the fast loop  !
! ACCUMULATES its sub-daily integral into xi_int (fed to soil_carbon_step) rather than passing a single !
! daily-mean xi -- so no Jensen bias (§2.4, §5.4). For spin-up it is evaluated at climatological-mean    !
! T/theta to build xi_bar (§4.1). xi > 0 (scheme-0 additionally <= 1; scheme-5 Q10 is not, §3.1).       !
!---------------------------------------------------------------------------------------------------!
pure subroutine assemble_env_scalar(soil_temp_surf, soil_temp_active, theta_active, theta_dry,      &
                                    theta_sat, pools, opts, xi)
   real(wp),             intent(in)  :: soil_temp_surf, soil_temp_active     !< [K]
   real(wp),             intent(in)  :: theta_active, theta_dry, theta_sat   !< [m3/m3]
   type(soil_carbon_t),  intent(in)  :: pools                                !< for f_lignin, f_decomp (N)
   type(decomp_opts_t),  intent(in)  :: opts
   real(wp),             intent(out) :: xi(n_soil_pool)                      !< [-] diagonal scalar, in (0,1]
end subroutine

!---------------------------------------------------------------------------------------------------!
! Assemble the (mostly static) transfer matrix A and decay diagonal K from config + texture + lignin.!
! A is rebuilt when f_lignin changes (structural respired-fraction blend); K is scheme-fixed.        !
!---------------------------------------------------------------------------------------------------!
pure subroutine assemble_transfer_matrix(pools, opts, a_mat, k_diag, er)
   type(soil_carbon_t), intent(in)  :: pools
   type(decomp_opts_t), intent(in)  :: opts
   real(wp),            intent(out) :: a_mat(n_soil_pool, n_soil_pool)   !< [-] a_ii=-1, a_ij transfer frac
   real(wp),            intent(out) :: k_diag(n_soil_pool)               !< [1/day] baseline decay rates
   real(wp),            intent(out) :: er(n_soil_pool)                   !< [-] respired fraction = 1 - sum_i a_ij
end subroutine

!---------------------------------------------------------------------------------------------------!
! Map the ALREADY-PARTITIONED, pool-destined litter fluxes onto the input vector u = B*I.            !
! The per-PFT labile/structural (f_labile) and above/below (agf) split is applied PER-COHORT by the  !
! driver BEFORE summing to the patch (§2.2, §5.4): each cohort's tissue carbon is split by ITS PFT's !
! f_labile_leaf/f_labile_stem and agf, then accumulated into the four litter destinations. This      !
! kernel therefore does NO f_labile arithmetic (it has no per-PFT data) — it just places the four     !
! pre-split streams + lignin into the pool vector. Routing CWD -> structural is part of the driver    !
! split (wood uses f_labile_stem). This resolves the earlier ambiguity: the DRIVER splits per-cohort; !
! the KERNEL maps the patch-summed, pool-destined fluxes.                                             !
!---------------------------------------------------------------------------------------------------!
pure subroutine build_litter_input(litter, u, lignin_in)
   type(litter_input_t), intent(in)  :: litter        !< pre-partitioned per-pool litter fluxes (driver-built)
   real(wp),             intent(out) :: u(n_soil_pool) !< [kgC/m2/day] input to each pool (1-4 nonzero)
   real(wp),             intent(out) :: lignin_in(2)   !< [kgC/m2/day] lignin flux to struct_grnd/struct_soil
end subroutine

!---------------------------------------------------------------------------------------------------!
! Advance the pools ONE DAILY step: dX/dt = u + A*(xi_int/dt)*K*X, decrementing each donor by the      !
! TRUE sub-daily loss integral. The step takes the fast loop's ACCUMULATED per-pool scalar integral    !
! xi_int_j = INT_day xi_j(t) dt [day] (ED2 today_* accumulators), NOT xi recomputed from daily-mean     !
! T/theta -- so the pool decrement equals the day's integrated atmospheric Rh by construction (closes   !
! the Jensen gap, §5.4). Structural lignin decays in proportion (L_s *= (1 - xi_int_s*K_s)); lignin_in  !
! adds fresh lignin (§2.1). Two solvers over the FROZEN A*(xi_int*K):                                   !
!   DECOMP_STEP_EULER -- ED2-faithful forward Euler (one step; the daily k*dt << 1 so already stable).  !
!   DECOMP_STEP_EXPM  -- exact matrix exponential of the frozen operator; its value is EXACTNESS for    !
!                        LARGE accelerated/spin-up steps, NOT daily-step stability (forward Euler is     !
!                        already stable at dt=1 day: even k_fast~11/yr gives k*dt~0.03).                 !
! rh_today is reported SOLVER-CONSISTENTLY as the ACTUAL pool loss to atmosphere, rh_today =            !
! litter_in + lignin-neutral - dC_pool (inter-pool transfers cancel in the total), so it matches the    !
! pool decrement for BOTH solvers and equals the fast loop's accumulated today_rh. It is an AUDIT/      !
! diagnostic value -- it is NOT re-fed to the CAS-CO2 twin (the fast loop is the sole CAS Rh authority, !
! §5.4). Both solvers mass-conserving to round-off.                                                     !
!---------------------------------------------------------------------------------------------------!
subroutine soil_carbon_step(pools, u, lignin_in, xi_int, opts, rh_today, audit)
   type(soil_carbon_t),      intent(inout) :: pools
   real(wp),                 intent(in)    :: u(n_soil_pool), lignin_in(2)
   real(wp),                 intent(in)    :: xi_int(n_soil_pool)             !< [day] INT_day xi_j dt (fast accumulator)
   type(decomp_opts_t),      intent(in)    :: opts
   real(wp),                 intent(out)   :: rh_today                        !< [kgC/m2/day] actual pool loss = fast today_rh
   type(soilc_audit_t),      intent(out)   :: audit                           !< dC = litter_in - rh_today; resid ~ 0
end subroutine

!---------------------------------------------------------------------------------------------------!
! Semi-analytic spin-up over the ACTIVE sub-block (indices with K>0; the full matrix is SINGULAR when   !
! inert pools have K=0, §4.1): select A = {j : K_j>0}, solve M*X_ss[A] = -(B*I_bar)[A] and hold inert   !
! pools at 0. Dispatches on structure: triangular back-substitution for the (triangular) default block, !
! small dense LU / Gaussian elimination w/ partial pivoting for the (non-triangular) scheme-5 block --   !
! both LAPACK-free, allocation-free, GPU-safe. No full inverse ever formed.                             !
!---------------------------------------------------------------------------------------------------!
pure subroutine solve_soil_carbon_steady_state(a_mat, k_diag, xi_bar, u_bar, opts, pools_ss)
   real(wp),            intent(in)  :: a_mat(n_soil_pool, n_soil_pool), k_diag(n_soil_pool)
   real(wp),            intent(in)  :: xi_bar(n_soil_pool), u_bar(n_soil_pool)   !< climatological means (prescribed at P0/P1)
   type(decomp_opts_t), intent(in)  :: opts                                      !< decomp_scheme -> active index set
   type(soil_carbon_t), intent(out) :: pools_ss
end subroutine

!---------------------------------------------------------------------------------------------------!
! Diagnostics (pure post-processing off the assembled matrices): residence time, storage capacity,    !
! storage potential, and the traceability decomposition. No state mutated.                            !
!---------------------------------------------------------------------------------------------------!
pure subroutine soil_carbon_diagnostics(a_mat, k_diag, xi, u, pools, diag)
   real(wp),               intent(in)  :: a_mat(n_soil_pool,n_soil_pool), k_diag(n_soil_pool)
   real(wp),               intent(in)  :: xi(n_soil_pool), u(n_soil_pool)
   type(soil_carbon_t),    intent(in)  :: pools
   type(soilc_diag_t),     intent(out) :: diag   !< tau(n), tau_E, X_c(n), X_p(n), trace tree components
end subroutine
```

Supporting types (in `meds_biogeochem_types`):

```fortran
type :: litter_input_t                    ! per-patch, per-day litter -- ALREADY PARTITIONED to pool destinations.
   ! The driver sums each cohort's leaf/fineroot/wood/storage necromass into these four labile/structural x
   ! above/below bins USING THAT COHORT'S PFT f_labile_leaf/f_labile_stem and agf (§2.2, §5.4). This type
   ! deliberately carries NO per-tissue, per-PFT breakdown -- the split is done upstream, per-cohort, so the
   ! kernel has everything it needs to map straight onto u (the earlier per-tissue form could NOT split by
   ! the per-PFT f_labile it never saw). Units [kgC/m2/day].
   real(wp) :: labile_grnd  = 0.0_wp      !< -> X(1) fast_grnd  (labile leaf/storage, above)
   real(wp) :: labile_soil  = 0.0_wp      !< -> X(2) fast_soil  (labile fineroot/storage, below)
   real(wp) :: struct_grnd  = 0.0_wp      !< -> X(3) struct_grnd (structural leaf + CWD, above)
   real(wp) :: struct_soil  = 0.0_wp      !< -> X(4) struct_soil (structural fineroot + belowground CWD)
   real(wp) :: lignin_grnd  = 0.0_wp      !< lignin flux to struct_grnd (sets f_lignin of incoming litter)
   real(wp) :: lignin_soil  = 0.0_wp      !< lignin flux to struct_soil
   ! optional N twin (driver-split by tissue C:N), present only when n_cycle_on:
   real(wp) :: n_labile_grnd = 0.0_wp, n_labile_soil = 0.0_wp
   real(wp) :: n_struct_grnd = 0.0_wp, n_struct_soil = 0.0_wp
end type litter_input_t

type :: soilc_audit_t                      ! daily carbon-mass conservation guard (the fast/slow contract)
   real(wp) :: litter_in    = 0.0_wp       !< [kgC/m2/day] sum(u)
   real(wp) :: rh_out       = 0.0_wp       !< [kgC/m2/day] Rh reported by soil_carbon_step (= litter_in - dC_pool)
   real(wp) :: rh_fast_accum= 0.0_wp       !< [kgC/m2/day] fast loop's accumulated today_rh (the CAS-fed flux)
   real(wp) :: dC_pool      = 0.0_wp       !< [kgC/m2/day] net pool change
   real(wp) :: resid        = 0.0_wp       !< [kgC/m2/day] dC_pool - (litter_in - rh_out); ~0 by construction (both solvers)
   real(wp) :: rh_seam_gap  = 0.0_wp       !< [kgC/m2/day] rh_out - rh_fast_accum; the fast/slow reconciliation check (~0)
   real(wp) :: lignin_resid = 0.0_wp       !< [kgC/m2/day] max_s |dL_s - (lignin_in_s - xi_int_s*K_s*L_s)|; passive-tracer check
end type soilc_audit_t

type :: decomp_opts_t                      ! pre-extracted decomposition selectors + params (like co2_opts_t)
   integer(ik) :: decomp_scheme = 0_ik     !< {0..5}; 0-4 = 3-active, 5 = full 5-pool CENTURY
   integer(ik) :: step_solver   = DECOMP_STEP_EULER   !< {EULER, EXPM}
   logical     :: n_cycle_on    = .false.  !< optional nitrogen (P1/P2)
   real(wp)    :: k_fast, k_struct, k_micr, k_slow, k_pass          !< [1/yr] baseline decay (verify vs ed_params)
   real(wp)    :: er_fast, er_struct_lig, er_struct_nonlig, er_micr, er_slow, er_pass
   real(wp)    :: e_lignin = 3.0_wp
   real(wp)    :: xsand = 0.0_wp, xclay = 0.0_wp                    !< texture -> er_micr, transfers to passive
   real(wp)    :: agf_fast = 0.5_wp, agf_struct = 0.7_wp           !< above/below split FALLBACKS (ED2 agf_fsc/agf_stsc);
                                                                   !< the driver splits per-COHORT with per-PFT agf/f_labile
                                                                   !< (PFT-file traits, §6) -- these are not used by the kernel
   real(wp)    :: c2n_structural=150.0_wp, c2n_slow=10.0_wp, c2n_fast=30.0_wp   !< N stoichiometry (N on)
   real(wp)    :: n_immobil_supply_scale = 0.0_wp
   ! + the environmental-scalar params (resp_temp_increase, rh_q10, resp_opt_water, ...) shared with co2_opts_t
end type decomp_opts_t
```

### 5.4 The slow-loop driver seam & the fast/slow carbon contract

The daily step hooks `meds_vegetation_dynamics` at ED2's `update_C_and_N_pools` slot — **after** phenology
+ carbon allocation (so the day's turnover/deficit litter exists) and mortality/treefall (so necromass
exists), **before** disturbance-driven patch restructuring, once per day. In the driver (`meds_aux`, the
only layer above demography + biogeochem + biophysics), per patch:

```
DURING the day (each fast substep, per patch) -- pools FROZEN & intent(in):
  a. xi_inst = assemble_env_scalar(T(t), theta(t), pools, opts)          [instantaneous, per-pool]
  b. Rh_inst = heterotrophic_respiration_flux(pools, xi_inst)  -> fed to canopy_air_co2_update (CAS)
  c. ACCUMULATE  xi_int(:) += xi_inst(:)*dt_sub ;  today_rh += Rh_inst*dt_sub     (ED2 today_* accumulators)
     (fast loop off / no forcing: xi_int = xi_bar*1day from a prescribed daily mean; today_rh from same)

AT the daily slot (meds_vegetation_dynamics, update_C_and_N_pools position), per patch:
  1. collect litter[]  <- per-cohort turnover/deficit (apply_carbon_npp's carbon_npp_t) + necromass
        (mortality_step + apply_patch_disturbance accumulator), SPLIT per-cohort by PFT f_labile/agf,
        summed into the pool-destined litter_input_t                                        [DATA]
  2. call build_litter_input(litter) -> u, lignin_in
  3. call assemble_transfer_matrix(pools, opts) -> A, K, er
  4. call soil_carbon_step(pools(ip), u, lignin_in, xi_int, opts) -> pools updated, rh_today, audit
       (decrements pools by the ACCUMULATED xi_int loss; rh_today = litter_in - dC_pool)
  5. audit guards (else error stop -- carbon not conserved across the day):
       |audit%resid| < tol                       (soil-side closure, both solvers)
       |audit%rh_seam_gap| < tol                  (rh_today == today_rh: pool decrement == CAS-fed flux)
       0 <= L_s <= C_s ;  |audit%lignin_resid| < tol
  6. zero xi_int, today_rh for the next day.
```

**The fast/slow contract (the whole correctness hinge): one Rh authority, the accumulated integral.** The
**fast loop is the sole source of the atmospheric Rh flux** — its per-substep `heterotrophic_respiration_flux`
off the frozen pool feeds `canopy_air_co2_update`, and *only* that flux reaches the CAS-CO2 twin. The daily
`soil_carbon_step` **never re-feeds a second Rh to CAS**; it does pool bookkeeping + inter-pool transfers,
and it decrements the donors by the fast loop's **accumulated** sub-daily loss `xi_int` (not by an Rh
recomputed from daily-mean T/θ). Because the two use identical chemistry and the *same* accumulated
integral, the pool decrement equals the day's integrated atmospheric flux **by construction** — no Jensen
gap (§2.4). The slow step is the pools' **sole writer**, once per day; the pools are `intent(in)` (enforced
by `pure`) in every fast kernel. Three separate guards catch the failure modes the fast twin's self-certifying
`resid ≡ 0` is blind to: the soil-side closure `ΔC_pool = litter_in − rh_today`, the **fast/slow seam check**
`rh_today ≡ today_rh` (the pool loss equals the CAS-fed flux — the check that actually catches a
double-applied or mismatched decomposition), and the lignin passive-tracer bound (§9, the CO2 design §9
risk 1 generalized to the slow loop).

**Which solver runs when (removes the EXPM/seam tension).** In **coupled daily operation** (fast loop on)
the pools are frozen over the day, so the correct daily advance is the **accumulated-loss forward step**
(`DECOMP_STEP_EULER` with `xi_int`): the pool loss is exactly the fast loop's frozen-pool integral, so
`rh_today ≡ today_rh` and the seam check is tight. `DECOMP_STEP_EXPM` is reserved for **offline spin-up /
large accelerated steps**, where there is *no* concurrent fast loop and hence no CAS flux to reconcile
against (`today_rh` is undefined); there `rh_today` is still reported solver-consistently as the actual
pool loss (`litter_in − dC_pool`) so the soil-side audit closes, but the `rh_seam_gap` check is simply not
exercised. The two never mix in one step.

Once the pools are prognostic, the fast Rh finally draws down a pool that is refilled by litter and
decremented by decomposition — **ecosystem carbon is conserved** for the first time (today the fast Rh
draws down a dead prescribed input; §1.1). **NEE is still only partly closed:** `NEE = GPP − Ra − Rh`, and
this module closes the **Rh** term, but committed autotrophic/growth/storage respiration **Ra is still 0**
in `cohort_co2_flux_t` (§1.1) — so NEE remains incomplete until Ra is wired (a plant-carbon/CO2 task, not
this module's; §9).

---

## 6. Config: `[soil_carbon]` / `[biogeochem]` TOML blocks

Two-part idiom matching `co2_opts_t` / `soil_opts_t`: `decomp_opts_t` (§5.3) carries in-type defaults so
the standalone kernels/tests compile pre-P3; run values are **REQUIRED from TOML** and presence-mapped
(no hard-coded model parameters — the MEDS provenance rule). The decay-rate and transfer parameters are
**loaded and verified against `ed_params.f90`** at config time (§2.1 — resolve the scheme-dependent
`k_slow` provenance here).

```toml
[soil_carbon]
decomp_scheme        = 0            # -> req_decomp_scheme mapper: 0..5 (0-4 = 3-active, 5 = full CENTURY)
step_solver          = "euler"      # "euler" (ED2-faithful) | "expm" (matrix exponential, stiff-stable)
n_cycle_on           = false        # optional nitrogen (P1/P2)
# baseline decay rates [1/yr] -- verified against ED2 ed_params.f90 at load
k_fast               = 11.0
k_struct             = 4.5
k_slow               = 0.2
k_micr               = 0.0          # 6.0 when decomp_scheme = 5
k_pass               = 0.0          # 0.012 when decomp_scheme = 5
# respired fractions [-]
er_fast              = 1.0          # 0.55 (scheme 5)
er_struct_nonlig     = 0.3          # 0.50 (scheme 5)
er_struct_lig        = 0.3          # 0.20 (scheme 5)
er_slow              = 1.0          # 0.55 (scheme 5)
er_pass              = 1.0          # 0.55 (scheme 5)
e_lignin             = 3.0
agf_fast             = 0.5          # above-ground init split (ED2 agf_fsc)
agf_struct           = 0.7          # (ED2 agf_stsc)
xsand                = 0.35         # texture -> er_micr, clay -> transfers to passive
xclay                = 0.25
# environmental-scalar params are SHARED with [co2] (resp_temp_increase, rh_q10, resp_opt_water, ...)

[soil_carbon.nitrogen]              # presence-mapped only when n_cycle_on = true
c2n_structural       = 150.0
c2n_slow             = 10.0
c2n_fast             = 30.0
n_immobil_supply_scale = 40.0       # [1/yr] converted by yr_day
```

Per-PFT litter-partition traits (`f_labile_leaf`, `f_labile_stem`, `c2n_leaf`, `l2n_stem`) go in the PFT
file via `req_pa('pft.<trait>', …, npft, miss)`. `validate_config` guards: all `k_* ≥ 0`,
`0 ≤ er_* ≤ 1`, `0 ≤ agf_* ≤ 1`, `0 ≤ xsand,xclay ≤ 1`, `decomp_scheme ∈ {0..5}`, and (scheme 5)
`k_micr > 0 ∧ k_pass > 0`, and (N on) all `c2n_* > 0`. A `req_decomp_scheme` mapper cloned from
`req_temp_response` (`select case`, `case default → note_missing` so an unknown value is a hard error).

---

## 7. Phasing P0 → P3

Consistent with the biophysics/biogeochem P0→P3 discipline (RT, hydrology, energy, CO2 all followed it):

- **P0 — pools + Q10 + matrix bookkeeping, wired to the existing fast Rh (the bulk of the value,
  ~standalone).** Expand `soil_carbon_t` to the 7-pool vector (index 2 keeps `fast_soil_carbon`, so
  `meds_column_co2` compiles unchanged); add `decomp_opts_t`, `litter_input_t`, `soilc_audit_t`,
  `n_soil_pool`/index params, `pack/unpack_pool_vector`. Create `meds_soil_biogeochem.f90` with
  `assemble_env_scalar`, `assemble_transfer_matrix`, `build_litter_input`, `soil_carbon_step`
  (`DECOMP_STEP_EULER`, ED2-faithful daily forward Euler), and the `Rh = −1ᵀAξKX` respired-complement.
  Litter input is a **prescribed scalar/array input** (no demography wiring yet); pools advanced with a
  plain daily Euler step; the fast Rh reads the (now-written) pool vector. CMake adds the module to the
  `meds_biogeochemistry` glob (auto). Tests 1-7 (§8). Everything bare-array, `pure` where possible,
  GPU-eligible. **No per-patch state, no config wiring, no driver seam** — exactly how soil hydrology/
  energy/CO2 landed their P0.
- **P1 — matrix spin-up + DAMM + C:N + diagnostics.** Add `solve_soil_carbon_steady_state` (SASU on the
  **active sub-block** — triangular back-substitution for the default block, small dense LU for the
  scheme-5 block, §4.1) and `soil_carbon_diagnostics` (residence time, capacity/potential, traceability;
  active block, unit convention §4.2). SASU runs on **prescribed climatological means** at this phase (no
  live met forcing yet) and is validated against brute-force on the same prescribed drivers (§4.1, test 5).
  Switch the preferred decomposition moisture response to **DAMM** (already in `meds_column_co2`; removes
  the `θ_dry` floor problem). Turn on the optional **N cycle** (parallel N pools, C:N stoichiometry,
  `f_decomp` structural brake). Add `DECOMP_STEP_EXPM` — the exact matrix exponential of the frozen
  `A·ξ·K`; its value is **exactness for the LARGE steps of accelerated/spin-up integration**, not daily-step
  stability (forward Euler is already stable at `dt=1 day`, where every `k·dt ≪ 1`). Still stateless
  (pools in, pools out).
- **P2 — vertically-resolved pools + N limitation of NPP.** Promote the per-patch scalar pools to
  **per-layer** pools (CLM-style depth-resolved SOC): `X` becomes `n_soil_pool × n_soil_layer`, with
  vertical transport (cryo-/bio-turbation, advection) as additional off-diagonal `A` entries on a
  block-tridiagonal system (reuses the promoted `meds_column_solver` Thomas sweep, mirroring the CO2
  multi-layer forward-design). Per-layer root-weighted soil T/θ drive per-layer ξ (replacing the
  single-`T_soil` collapse). **N limitation of NPP** feeds an `f_decomp`-analogous down-regulation back
  into the plant carbon seam (`get_plant_flux_slow`) — the demand-vs-supply nutrient gate. Fire as a new
  disturbance class consumes the above-ground ("grnd") pools.
- **P3 — per-patch state + config + driver coupling (the deferred integration all biophysics stores
  share).**
  - **State lockstep:** `soil_carbon_t(ip)` joins `patch_index`; thread through `patch_alloc` /
    `patch_ensure_capacity`, `sort_patches` **and** `patch_compact`, the `fuse_2_patches` area-fraction
    blend, and the disturbance donor-copy (§5.2) — in lockstep beside `can_co2` / `cas` / `soil_e` /
    `soil_w`. **Serialize the expanded state in `meds_io`:** the 7 carbon pools + 2 lignin + 5 N fields
    join `io_write_state`/`io_read_state` and the netCDF state schema (§5.2) — else a restart silently
    zeroes soil carbon. Cold-start zeros the pools; `init_bare_ground` / the optional SASU stamp set the
    initial state (§5.2).
  - **Config:** `[soil_carbon]` into `meds_config_t` / `meds_config_io` / `validate_config` /
    `derive_config`; add `meds_biogeochemistry` to the `meds_aux` link line.
  - **Driver seam (single Rh authority):** demography emits necromass into a per-patch daily accumulator
    (a `carbon_flux_block` sibling, plain data — no plant/biogeochem edge); the driver splits litter[]
    per-cohort by PFT `f_labile`/`agf` into the pool-destined `litter_input_t`, reads soil T/θ from the
    biophysics reservoirs, and calls `soil_carbon_step` once per day at the `update_C_and_N_pools` slot in
    `meds_vegetation_dynamics`, passing the fast loop's **accumulated `xi_int`** and running the audit
    guards. The **fast loop remains the sole Rh authority into `column_co2_step`** (its per-substep Rh off
    the frozen pool) — the daily step does **not** feed a second Rh to the CAS twin; it only decrements the
    pools by the accumulated loss so the pool change matches the CAS-fed flux (§5.4). The fast/slow contract
    is enforced: pools `intent(in)` in the fast loop, written once per day.

---

## 8. Test plan & validation vs ED2

Program `test_soil_biogeochem` in the house style (`check` / `check_true`, `nfail`, `error stop 1`),
linking `meds_biogeochemistry`. **Build under both ifx *and* nvfortran multicore** — a green ifx run is
not sufficient (CLAUDE.md issue #7; the array-temporary miscompile is silent at `-O2`, and this module is
array-heavy). Never pass an array-valued function result straight into a call — bind to a named array
first.

1. **`soilc_mass_closure`** — randomized pools, litter, `xi_int`: assert `|audit%resid| < 1e-12·max(ΣX,1)`
   for one `soil_carbon_step`. **Because `rh_out` is reported solver-consistently as the ACTUAL pool loss
   `rh_out = litter_in − dC_pool`** (inter-pool transfers cancel in the total; §5.3), `resid ≡ 0` holds by
   construction for **both** `EULER` and `EXPM` — the earlier instantaneous-`−1ᵀAξKX` definition would
   *not* close for `EXPM` (the exact integral ≠ the start-of-day rate). The core conservation guarantee.
2. **`rh_respired_complement`** — assert `Rh = Σ_j er_j·ξ_j·K_j·X_j = −1ᵀ·A·ξ·K·X` recomputed
   independently; and that the transferred + respired fractions sum to the total loss per pool
   (`Σ_i a_ij + er_j = 0` column check, `a_jj = −1`).
3. **`scalar_placement`** — a per-pool ξ with distinct grnd/below values: assert `A·ξ·K·X` differs from
   the wrong `ξ·A·K·X`, and that the *correct* form conserves mass while the wrong one does not (guards
   §3.3, the load-bearing error).
4. **`decomp_scheme_topology`** — scheme 0 gives the 3-active matrix (microbial/passive inert, structural
   → slow, 100% respired fast/slow); scheme 5 gives the 5-active clay/sand-controlled topology; both
   reproduce hand-computed `A`.
5. **`spinup_steady_state`** — `solve_soil_carbon_steady_state` on the **active sub-block** with constant
   prescribed `(ξ̄, ū)`: assert `soil_carbon_step` from `X_ss` returns `dX/dt ≈ 0` to round-off (the analytic
   equilibrium is a fixed point of the integrator) and inert pools stay 0; and that a long forward
   integration from bare ground converges to the same `X_ss` (SASU vs brute force agree). Run **both** the
   default 3-active block (triangular back-substitution) and the scheme-5 5-active block (dense LU) — the
   scheme-5 block is non-triangular (§3.2). A guard asserts that attempting the **full 7×7** solve is
   rejected/short-circuited (it is singular, `K_5=K_7=0`, §4.1).
6. **`residence_capacity`** — `X_ss = τ·ū` (residence time reproduces the steady state); `X_c = X_ss` at
   equilibrium and `X_p = 0`; a perturbed pool relaxes toward `X_c` (`dX/dt = −A·ξ·K·X_p`).
7. **`litter_partition` + lignin tracer** — two parts. (a) The **driver-side** per-cohort split: a
   two-PFT cohort set with distinct `f_labile_leaf`/`f_labile_stem`/`agf` sums into a pool-destined
   `litter_input_t` whose `labile_*`/`struct_*` bins match a hand computation; CWD (wood) lands in
   structural via `f_labile_stem`. (b) The **kernel** `build_litter_input` maps that record straight onto
   `u`/`lignin_in` with no further splitting (it has no per-PFT data). Then assert the **lignin tracer
   law** over a step: `L_s` decays in proportion to `C_s` so `f_lignin` is conserved under pure decay,
   `lignin_in` raises it, and `0 ≤ L_s ≤ C_s` with `|audit%lignin_resid| < tol` (§2.1).
8. **`temperature_moisture_response`** — daily ξ ≈ doubles per 10 K at `rh_q10 = 2` (Q10) and saturates at
   45 °C (`HR_EXP_ED2`); moisture peaks at `resp_opt_water`; `Rh = 0` at `X = 0`; a nonzero `(X, k)` gives
   nonzero `Rh` (the stub-zeroing guard). **Reproduces the fast Rh's `f_temp`/`f_water` at matched inputs**
   (the fast/slow reconciliation, §5.4).
9. **`daily_carbon_audit` (fast/slow contract, the seam)** — drive a synthetic **diurnal** T/θ cycle;
   accumulate the fast per-substep Rh + `xi_int` over the day at frozen pools; step the pools with that
   `xi_int`. Assert **`|audit%rh_seam_gap| = |rh_today − today_rh| < tol`** (the pool decrement equals the
   CAS-fed flux, by construction). **Counter-check the Jensen gap:** recomputing ξ from the daily-mean
   `(⟨T⟩,⟨θ⟩)` instead yields a *nonzero* seam gap under the diurnal cycle — the test asserts the
   accumulator path closes and the daily-mean path does **not** (guarding §2.4/§5.4). Assert the pool
   `intent(in)` contract is never violated (unchanged by any fast kernel).
10. **N-cycle closure** (P1) — with `n_cycle_on`, C:N stoichiometry holds per pool; `f_decomp` brakes
    structural decay under N shortage and `≡ 1` under surplus; mineralization = fast+structural N losses.
11. **`nvfortran_multicore_build`** — the module compiles and tests pass under `nvfortran -mp`.

**Validation vs ED2 (the EDTS analogue).** Beyond unit tests, the science acceptance is a
whole-module comparison against ED2's `soil_respiration` + `update_C_and_N_pools` on a controlled patch:
identical initial pools, litter inputs, soil T/θ time series, and `DECOMP_SCHEME` → assert MEDS pools and
daily Rh track ED2's within tolerance over a multi-year run (both the 3-active default and scheme-5).
The **traceability tree** (§4.3) is the diagnostic when they diverge — it attributes any steady-state
soil-C difference to input, `τ'_E` (allocation/transfer/turnover), or the environmental scalars, isolating
the offending coefficient. Spin-up is validated by SASU-vs-brute-force agreement (test 5) and by the
equilibrium `X_ss` matching an ED2 long spin-up on the same climatological drivers.

---

## 9. Open questions, risks & deferred work

**The biggest risk — carbon-mass conservation across the fast/slow seam.** The fast CAS-CO2 kernel's
`resid ≡ 0` certifies only the twin's internal algebra; it is blind to whether the pool was correctly
refilled and decremented. If a contributor lets a fast kernel mutate the pool, or double-applies
decomposition (once in the daily matrix step *and* again through the fast Rh), the per-substep CO2
residual still reads ~0 while ecosystem carbon silently drifts. **Mitigation (three-legged):** (1) **one Rh
authority** — the fast loop is the sole flux into CAS; the daily step feeds *no* second Rh, so there is no
double-count path to begin with (§5.4). (2) The daily step decrements the pools by the fast loop's
**accumulated `xi_int`**, not by a recomputed daily-mean ξ, so the pool loss equals the integrated
atmospheric flux by construction — no Jensen drift (§2.4). (3) Pools are `intent(in)` in every fast kernel
(enforced by `pure`), and the daily audit runs a *separate* mass check plus the explicit
**`rh_today ≡ today_rh` seam check** (test 9) that a single-side soil audit would miss.

**Air-dry moisture floor `θ_dry`** — a correctness dependency inherited from the fast Rh: resolve to a
real `soilcp`/retention air-dry value, never silently to `θ_res`. The DAMM path (P1) removes the problem
mechanistically and is the recommended default once available.

**`k_slow` provenance** — the ED2 memory-brief carried an ambiguous concatenated value ("100.2");
the exact per-scheme decay rates must be read from `ed_params.f90` and asserted at config load (§6),
not transcribed from the brief.

**Scheme-5 quasi-linearity + non-triangular back-transfers** — the clay/sand-controlled transfers and (N
on) the C:N-dependent `f_decomp` make `A` mildly state-dependent, so SASU is an iterated fixed point, not
a one-shot inverse (§4.1); treat the analytic solve as the inner step of a few outer cycles. Separately,
the microbial⇄slow and slow⇄passive **back-transfer cycles** put entries above the diagonal, so scheme
5's active block is **genuinely non-triangular** — it needs a real dense solve (LU with partial pivoting,
§3.2/§4.1), *not* the back-substitution that suffices for the default block. Getting this wrong (assuming
triangularity) silently corrupts the scheme-5 equilibrium.

**SASU is a mean-driver approximation, not the true attractor** — replacing `ξ(t),I(t)` by their means
discards the `⟨ξ·X⟩ − ξ̄·X̄` covariance of periodic forcing (Xia-2012 bias). The outer fixed-point loop
reduces but does not remove it; it is present even in the linear default scheme, independent of scheme-5/N
state-dependence (§4.1). And at P0/P1 the means themselves are **prescribed**, since litter→ξ̄/Ī needs the
absent met forcing.

**NEE is only partly closed here** — this module closes the **Rh** term of `NEE = GPP − Ra − Rh`, but the
committed autotrophic/growth/storage respiration **Ra is still 0** in `cohort_co2_flux_t` (§1.1, §5.4). NEE
remains incomplete until Ra is wired (a plant-carbon/CO2-twin task); the soil-C module must not be read as
"NEE closed".

**Restart/checkpoint & cold-start** — the expanded prognostic state (7 C pools + 2 lignin + 5 N) MUST be
added to `io_write_state`/`io_read_state` and the netCDF state schema at P3, or a restart silently zeroes
soil carbon; cold-start initializes the pools/lignin/N to 0 (bare ground), with SASU (§4.1) the intended
route to a realistic pre-run state (§5.2).

**Vertical resolution deferred (P2)** — MVP pools are per-patch scalars like ED2, with only the climate
scalars using layer-resolved T/θ; depth-resolved SOC (and its block-tridiagonal vertical transport) is a
P2 upgrade sharing the promoted `meds_column_solver`.

**Deferred work.** Explicit coarse-woody-debris pool (MVP folds CWD into structural, ED2-faithful); fire
consumption of the above-ground pools (needs the fire disturbance class); full N cycle beyond ED2's
minimal scope (nitrification/denitrification, gaseous loss, deposition, fixation — all absent in ED2 too);
phosphorus (the CABLE/CASA-CNP extension the matrix framework bolts on additively via a parallel P matrix);
and C-N-P coupled data assimilation (the matrix form's analytic spin-up is what makes Bayesian/ML
calibration of the soil-C parameters tractable — Hararuk 2014, Tao 2020 — a research direction, not a
deliverable).

---

## 10. File manifest (all absolute)

**New.**
- `/home/xiangtao/projects/MEDS/src/biogeochemistry/meds_soil_biogeochem.f90` — the slow kernels
  (`assemble_env_scalar`, `assemble_transfer_matrix`, `build_litter_input`, `soil_carbon_step`,
  `solve_soil_carbon_steady_state`, `soil_carbon_diagnostics`); links `meds_shared` only.
- `/home/xiangtao/projects/MEDS/test/test_soil_biogeochem.f90` — tests 1-11 (§8).

**Edit.**
- `/home/xiangtao/projects/MEDS/src/biogeochemistry/meds_biogeochem_types.f90` — expand `soil_carbon_t`
  to the 7-pool vector + lignin + optional N (§2.1); add `litter_input_t`, `soilc_audit_t`,
  `soilc_diag_t`, `decomp_opts_t`, `n_soil_pool`/`IP_*` index params, `DECOMP_STEP_*`/`decomp_scheme`
  selectors, `pack/unpack_pool_vector`.
- `/home/xiangtao/projects/MEDS/src/state/meds_demography_types.f90` (P3) — `soil_carbon_t(ip)` on
  `patch_index`; `patch_alloc` / `patch_ensure_capacity`.
- `/home/xiangtao/projects/MEDS/src/demography/meds_demography_fusefiss.f90` (P3) — thread the pools
  through `sort_patches`, `patch_compact`, the `fuse_2_patches` area-fraction blend.
- `/home/xiangtao/projects/MEDS/src/demography/meds_demography_dynamics.f90` (P3) — necromass emission
  from `mortality_step` / `apply_patch_disturbance` into a per-patch daily accumulator (plain data).
- `/home/xiangtao/projects/MEDS/src/driver/meds_vegetation_dynamics.f90` (P3) — collect litter[], read
  soil T/θ, call `soil_carbon_step` at the `update_C_and_N_pools` slot, run the audit, feed Rh to
  `column_co2_step`.
- `/home/xiangtao/projects/MEDS/src/shared/meds_config.f90` + `.../io/meds_config_io.f90` (P3) —
  `[soil_carbon]` block, `req_decomp_scheme` mapper, `validate_config` guards, `ed_params.f90` provenance
  check.
- `/home/xiangtao/projects/MEDS/src/io/meds_io.f90` (P3) — serialize the expanded prognostic
  `soil_carbon_t` (7 C pools + 2 lignin + 5 N fields) in `io_write_state` / `io_read_state` and the netCDF
  state schema, beside `can_co2` / `soil_e` / `soil_w`; keep `meds_io_stub.f90` in sync (§5.2). Without
  this a checkpoint/restart silently loses all soil carbon.
- `/home/xiangtao/projects/MEDS/src/shared/meds_constants.f90` — `yr_day` (if absent) for `/yr → /day`
  decay-rate conversion (`kgCday_2_umols` already present from the CO2 design).
- `/home/xiangtao/projects/MEDS/CMakeLists.txt` — `test_soil_biogeochem` after `test_column_co2`; the
  `meds_biogeochemistry` glob auto-picks the new module.
- `/home/xiangtao/projects/MEDS/src/biogeochemistry/README.md` + `/home/xiangtao/projects/MEDS/CLAUDE.md`
  — note the slow soil-C matrix module now realizes the "owns fast *and* slow carbon" charter.

**Reference kernels to mirror.**
- `/home/xiangtao/projects/MEDS/src/biogeochemistry/meds_column_co2.f90` — the fast Rh + CAS-CO2 twin this
  feeds (`heterotrophic_respiration_flux`, `column_co2_step`).
- `/home/xiangtao/projects/MEDS/src/plant/meds_plant_carbon_dynamics.f90` — the turnover/deficit litter
  streams (`carbon_npp_t`).

**ED2 references.**
- `/home/xiangtao/projects/ED2/ED/src/dynamics/soil_respiration.f90` (`het_resp_weight`, `resp_rh`,
  `resp_f_decomp`, `root_resp_norm`, the sub-daily flux accumulation).
- `/home/xiangtao/projects/ED2/ED/src/dynamics/vegetation_dynamics.f90` (`update_C_and_N_pools` — the
  daily pool advance + `check_budget_soilc`).
- `/home/xiangtao/projects/ED2/ED/src/memory/decomp_coms.f90`, `.../init/ed_params.f90` (pools, decay
  rates, transfer/respired fractions, C:N ratios — the provenance source for §6).
- `/home/xiangtao/projects/ED2/ED/src/dynamics/structural_growth.f90`, `.../growth_balive.f90`
  (`update_litter_inputs` — the litter staging MEDS's `build_litter_input` mirrors).

**Matrix references.** Luo et al. 2022 (JAMES `10.1029/2022MS003008`, synthesis); Xia et al. 2012 (GMD
`10.5194/gmd-5-1259-2012`, SASU); Xia et al. 2013 (GCB `10.1111/gcb.12172`, traceability); Luo et al. 2017
(BG `10.5194/bg-14-145-2017`, capacity/potential); Huang et al. 2018 (GCB `10.1111/gcb.13948`,
matrix-CLM4.5); Lu et al. 2020 (JAMES `10.1029/2020MS002105`, matrix-CLM5 C+N); Sierra et al. 2012 (GMD
`10.5194/gmd-5-1045-2012`, SoilR compartmental); Bolker/Pacala/Parton 1998 (Ecol. Appl. 8:425, CENTURY
scheme-5 structure); Koven et al. 2013 (BG 10:7109, scheme-5 Q10 + moisture/oxygen); Parton et al. 1987
(SSSAJ 51:1173, CENTURY 3-pool origin).
