# MEDS Plant Hydraulics — Module Design

A **stateless**, per-individual plant-hydraulics compute library for MEDS (`src/plant/hydraulics/`).
It takes ED2 **X16** (Xu et al. 2016) as the physical reference but revises it to be *more accurate*
(coupled nodes, nonlinear pressure–volume, Kirchhoff-integrated conductance, adaptive sub-stepping)
and *more efficient* (fixed small dimension, closed-form linear algebra, branch-light smooth
limiters — GPU/SIMD-friendly). It mirrors `src/plant/leaf/` exactly: links `meds_shared` only, no
`site_t` dependency, compiles and tests standalone, one public seam. The prognostic state (node
water potentials `psi`) lives in the cohort SoA and is passed by argument — the FATES `*Mem`/compute
split. Soil hydrology is out of scope: the rhizosphere boundary `(soil_psi, rhizo_cond)` is a
*given* scalar BC.

The goal is to reproduce ED2-hydro **qualitatively** (diurnal ψ depression, nocturnal recovery,
drought down-regulation) with better physics and numerics — not to match X16 bit-for-bit. The fine
explicit reference integrator (§12) is the accuracy ground truth.

**Units are MPa throughout** (never meters of head): matches `leaf_env_t%psi_leaf` and the
`[pft].wstress_psi_*` TOML. Kinds `wp/ik`, `_wp` literals, `implicit none`, `pure`/`elemental`
kernels, `error stop`, ≤132 cols.

---

## 1. Module layout

Under `src/plant/hydraulics/`, mirroring `src/plant/leaf/`. One static lib `libmeds_hydraulics.a`
(`PUBLIC meds_shared`), one public seam.

| File | Role | Analogue in `leaf/` |
|------|------|---------------------|
| `meds_hydro_types.f90` | `hydro_env_t` (BCs), `hydro_flux_t` (outputs), `hydro_params_t` (flat per-PFT), node/solver enums, `N_HYDRO` parameter | `meds_leaf_types.f90` |
| `meds_hydro_pv.f90` | PV kernels: `water_content(ψ)`, `rwc_from_psi(ψ)`, `capacitance(ψ)` (Bartlett/Christoffersen nonlinear, closed-form inverse), `pv_water_cap_from_traits` (legacy-linear mapping) — `pure`/`elemental` | (part of `photosynthesis`) |
| `meds_hydro_conductance.f90` | `plc_retained`, `flux_potential`(Φ, closed-form `kexp∈{1,2}` or quadrature) + `phi_inverse`, `kirchhoff_edge` (`K_eff=ΔΦ/Δψ` with `Δψ→0` guard), `k_plant_max_from_allometry` (segment mode), `stem_psi_at_height` diagnostic, optional `rhizosphere_cond` (test helper) — `pure`/`elemental` | `meds_leaf_stomata.f90` |
| `meds_hydro_solver.f90` | coupled adaptive stepper `solve_plant_water`; internal `pure` `CONTAINS` helpers; no allocatables | `meds_leaf_solver.f90` |
| `meds_plant_hydraulics.f90` | **THE seam** `plant_water_flux(env, cfg, ipft, dt, psi, flux)`; flattens `cfg%pft`→`hydro_params_t`, reads `cfg` selectors, calls the solver | `meds_leaf_physiology.f90` |
| `meds_hydro_capi.f90` | optional `bind(c)` shim (`-DMEDS_BUILD_PYLIB`), excluded from core lib | `meds_leaf_capi.f90` |
| `test/test_plant_hydraulics.f90` | CTest target (links `meds_testsupport` + `meds_hydraulics`) | `test/test_leaf_physiology.f90` |

CMake: mirror the leaf-physiology block (`add_library(meds_hydraulics STATIC …)`,
`target_link_libraries(meds_hydraulics PUBLIC meds_shared)`, `meds_fortran_flags(…)`, exclude the
`_capi.f90` from the core glob).

---

## 2. The stateless contract

```fortran
subroutine plant_water_flux(env, cfg, ipft, dt, psi, flux)
   type(hydro_env_t),   intent(in)    :: env          ! boundary conditions (read-only)
   type(meds_config_t), intent(in)    :: cfg          ! selectors + shared constants
   integer(ik),         intent(in)    :: ipft
   real(wp),            intent(in)    :: dt           ! [s] fast step
   real(wp),            intent(inout) :: psi(N_HYDRO) ! [MPa] node potentials — SoA slice, owned OUTSIDE
   type(hydro_flux_t),  intent(out)   :: flux         ! sapflow, root_uptake, psi_leaf, plc, nsub, converged
end subroutine
```

* `psi` is `intent(inout)`, a **fixed-size** `real(wp)` array of length `N_HYDRO` (compile-time
  parameter) — a slice of the cohort SoA handed in by the caller, the *only* mutable thing and not
  module state. Two threads advancing two cohorts touch disjoint `psi(:)` ⇒ pure/reentrant.
* Everything in `hydro_env_t` is a boundary condition: transpiration `E`, the aggregated rhizosphere
  `soil_psi`/`rhizo_cond`, and the per-plant geometry/biomass (frozen over the step).
* `hydro_params_t` is flattened from `cfg%pft(ipft)` inside the seam (like `leaf_photo_params_t`);
  `solve_plant_water` never reaches into `meds_config`.
* No `save`, no module variables, no allocatables in the hot path. The solver is `pure` end to end.

| Quantity | Lives in | Passed how |
|----------|----------|-----------|
| `psi(N_HYDRO)` node potentials | cohort SoA (`src/state/`) | `intent(inout)` |
| `E`, `soil_psi`, `rhizo_cond` | env / other modules | `hydro_env_t` (in) |
| geometry `bleaf`,`bsap`,`broot`,`sap_area`,`height` | cohort SoA (cached) | `hydro_env_t` (in) |
| PFT hydraulic traits | `pft_table_t` | flattened → `hydro_params_t` |

---

## 3. Node topology

Fixed compile-time `integer(ik), parameter :: N_HYDRO = 3`. The active count is a config enum; the
solver loops `1:n_active`, `n_active ∈ {2,3}`, so the array shape never changes.

| Topology | Nodes | Edges | When |
|----------|-------|-------|------|
| **`HYDRO_NODES_2`** (**default**) | Leaf `L`, Wood `W` (stem+root lumped) | `L–W` (whole-plant xylem), `soil–W` (rhizosphere) | production first cut; the internal edge is the whole-plant conductance (§4) |
| **`HYDRO_NODES_3`** (opt-in) | Leaf `L`, Stem-base `T`, Root `R` | `L–T`, `T–R`, `soil–R` (chain) | separate root storage/vulnerability, emergent hydraulic redistribution, ψ(z) diagnostic (§6.4) |

Both are chains ⇒ tridiagonal `A` ⇒ 2-node is a hand-unrolled 2×2, 3-node a 3-unknown Thomas sweep
(no LAPACK). `NODE_LEAF=1, NODE_STEM=2, NODE_ROOT=3`; in the 2-node run `psi(3)` is unused and
`W ≡ node 2`.

The **2-node default** is the cleanest home for the whole-plant-conductance parameterization (§4):
a single internal edge `L–W` = whole-plant xylem conductance, plus the given `soil–W` rhizosphere
edge — i.e. Bartlett's Ohm's-law plant with two capacitive pools. **Node heights:** `z_leaf =
height`, `z_wood = 0`. `N_HYDRO` stays 3 so the 3-node upgrade is a config flip, not an ABI change;
when it lands, split `W` into stem-base `T` (`z=0`) and root `R` (`z=0`) with `k_plant_max`
partitioned into series segment conductances (§4), and the soil edge moves to `R`.

**Soil is the root boundary condition** — the caller supplies the aggregated scalar pair
`(soil_psi, rhizo_cond)`; the `root_beta`/RAI aggregation belongs to a future soil module. A
test-only pure helper `rhizosphere_cond(soil_cond, broot, sra, root_frac, dz, nplant)` reproduces
the Katul-2003 `soil_cond·√RAI/(π·dz)/n` form. (Per-layer root nodes are a planned extension, §16.)

---

## 4. Governing equations (MPa)

Total head `Ψ_i = ψ_i + grav_head·z_i`, `grav_head = 9.804e-3` MPa/m (`meds_constants`). Flux `j→i`:
`q = K_ij(ψ_j − ψ_i + g_ij)`, `g_ij = grav_head·(z_j − z_i) = −g_ji`. Per-node mass balance
(kg H₂O per plant):

```
C_i(ψ_i) dψ_i/dt = Σ_{j~i} K_ij(ψ) (ψ_j − ψ_i + g_ij) − S_i          (1)
```

* `C_i` **capacitance** [kg MPa⁻¹] (§6). Unit box: if porting an X16 per-meter value,
  `cap[kg/kgC/MPa] = cap_X16[kg/kgC/m] × 101.97`. Define TOML traits natively in per-MPa.
* `K_ij` **internal xylem conductance** [kg MPa⁻¹ s⁻¹] — parameterized by **whole-plant
  conductance** (default, `HYDRO_COND_KPLANT`): `k_plant_max` is the leaf-area-specific maximum
  whole-plant conductance [kg s⁻¹ MPa⁻¹ m⁻²_leaf], the measured quantity and exactly the `K_plant`
  in Bartlett's Ohm's law `ψ_leaf = ψ_soil − E/K_plant`. **Vulnerability enters via the Kirchhoff
  flux potential only — never pointwise** (§6.2): the edge conductance is the ψ-averaged retained
  fraction over the drop, `K_LW = k_plant_max·leaf_area·⟨plc⟩`, `⟨plc⟩ = ΔΦ/(k_plant_max·Δψ)`.
  (Alternative `HYDRO_COND_SEGMENT`: derive `k_plant_max = wood_kmax·sap_area/(height·vessel_curl·
  leaf_area) = wood_kmax·H_v/(height·vessel_curl)` from stem allometry, `H_v = sap_area/leaf_area`
  the Huber value — same edge law, size-dependent conductance capturing tall-tree hydraulic
  limitation.) The soil edge is always the given `rhizo_cond`.
* `g_ij` gravity offset [MPa], a constant in the RHS.
* `S_L = E` (transpiration sink, kg s⁻¹); `S_R = 0`; soil enters node `R`/`W` as
  `rhizo_cond(soil_psi − ψ + g)`.

**Whole-plant vs edge k are one series network.** For a multi-node chain
`1/K_plant = Σ_edges 1/K_edge`. In 2-node the single internal edge *is* `K_plant` (nothing to
partition). For 3-node, partition `k_plant_max` into series segment conductances via resistance
fractions `f_leaf+f_stem+f_root = 1` (`1/K_edge_i = f_i/(k_plant_max·leaf_area)`), so whole-plant k
stays the top-level trait and the edge k's are derived — the solver (§5) never changes, only how the
`A`-matrix entries are filled.

**2-node system** (`g_WL = −grav_head·height`): `C ψ' = A ψ + b`, `ψ' = M ψ + c`, `M = C⁻¹A`. `A` is
the grounded Laplacian (symmetric negative-definite, grounded through `rhizo_cond > 0`), so `M` is
similar to symmetric `Ĥ = C^{-1/2}AC^{-1/2}` ⇒ **all eigenvalues of `M` are real ≤ 0 and `M` is
diagonalizable** — no complex/defective branch; the §5 closed forms are always valid. This is X16's
two pools but **coupled** (X16's leaf-then-wood split is the first-order error removed in §5.3).

---

## 5. The solver

**Default: coupled, adaptively sub-stepped, frozen-coefficient — exact on each sub-step.**

| Enum | Method | Use |
|------|--------|-----|
| `HYDRO_SOLVER_EXPM` (**default**, pairs with 2-node) | closed-form 2×2 **matrix exponential** of the frozen system | production 2-node; exact frozen step, cheapest |
| `HYDRO_SOLVER_BE` (3-node / stiff) | linearly-implicit **backward Euler** (Rosenbrock-1), 3×3 Thomas / 2×2 Cramer, L-stable | 3-node chain; near-cavitation / leafless-floored |

`hydro_enable = 0` is the off path (ψ held saturated at 0). On `[t,t+h]` freeze `M,c`; exact frozen
solution `ψ(t+h) = ψ* + e^{Mh}(ψ(t) − ψ*)`, `ψ* = −M⁻¹c`.

### 5.1 Underflow-safe 2×2 matrix exponential

With `N = Mh`, `μ = ½tr N`, `Δ = √max(μ² − det N, 0)`:

```
e_p = safe_exp(μ+Δ) ; e_m = safe_exp(μ−Δ)          ! both ≤ 1: no overflow
CH  = 0.5*(e_p + e_m)
if (Δ > sinhc_eps) then
   SH = 0.5*(e_p − e_m)/Δ
else
   SH = CH                                          ! smooth Δ→0 limit; NO 0/0
end if
e^{Mh} = [ CH+SH*(N11−μ)  SH*N12 ;  SH*N21  CH+SH*(N22−μ) ]
ψ_new  = ψ* + e^{Mh} (ψ − ψ*)
```

L-stable for the fast (embolism) mode. `ψ*` uses the closed 2×2 inverse with a `det` floor. Do
**not** use `merge` for the `sinhc` limit (its dead branch's `0/0` trips `-fpe0`/`-Ktrap=fp`).

### 5.2 Backward Euler / Rosenbrock-1 (3-node & stiff 2-node)

`(C/h − J) Δψ = F(ψⁿ)`, one linear solve (3×3 Thomas / 2×2 Cramer). `J` uses the finite `dplc/dψ`
(§6.2) on the `C`-side and the bounded Kirchhoff edge Jacobian. L-stable: `λ → −∞`
(leafless-floored `C_L`) ⇒ amplification `1/(1−hλ) → 0` with no transcendental. Auto-selected when a
node is at its `C`-floor or `|λ_min|h` is large.

### 5.3 Coupling vs operator splitting (accuracy)

X16 = Lie–Trotter split `e^{M_W h}e^{M_L h}` as one `h = DTLSM` step. BCH gives local error
`−½h²[M_L,M_W]+O(h³)` (first order), `‖·‖ ≈ (h²/2)·K_LW²/(C_L C_W)·|ψ_L−ψ_W+g_WL|`, maximal when
`C_L ~ C_W` with large `K_LW`, and X16 never sub-cycles. The coupled `expm` evaluates
`e^{(M_L+M_W)h}` — the commutator term never appears — and adaptive sub-stepping (§5.4) drives the
residual frozen-coefficient error to tolerance. Strictly more accurate on two counts, branch-free.

### 5.4 Adaptive sub-stepping (step doubling)

Re-freeze `M,c` (plc, K, C) at each sub-step start. Estimate error one step `h` vs two `h/2`;
`err = max_i |ψ_{h/2}−ψ_h|/(atol+rtol|ψ|)`; accept extrapolated `ψ_{h/2}` if `err ≤ 1`;
`h ← h·clip(safety·err^{−1/(p+1)}, fmin, fmax)`, `p=1`; cap `n_sub ≤ hydro_max_substep`.
Well-watered midday ⇒ ~1 sub-step; near-`ψ50` ⇒ auto-refine to `O(10¹–10²)`. **GPU:** variable
`n_sub` per cohort is warp divergence; `HYDRO_SUBSTEP_FIXED` runs a fixed count for the offload path.

---

## 6. Constitutive curves (`pure`/`elemental`)

### 6.1 Pressure–volume / capacitance — nonlinear (Bartlett/Tyree–Hammel), only

Symplastic relative water content `R ∈ (0,1]`, `R = 1` at full turgor. Fundamental PV traits (per
tissue, leaf/wood): `pi0 = π_o` osmotic potential at full turgor [MPa, <0], `eps = ε` bulk modulus
of elasticity [MPa, >0], `af` apoplastic fraction [0–1]. **Derived in `derive_hydro_params`
(Bartlett et al. 2012, Ecol. Lett. 15:393, eqns 1–2):**

```
psi_tlp = pi0*eps/(pi0 + eps)                  ! turgor loss point  [MPa]   (Bartlett eqn 1)
rwc_tlp = (pi0 + eps)/eps                       ! symplastic RWC at TLP      (Bartlett eqn 2)
```

Symplast water potential (Tyree & Hammel 1972; the model analysed by Bartlett 2012):

```
turgid  (rwc_tlp ≤ R ≤ 1):  psi(R) = eps*(R − rwc_tlp) + pi0/R      ! turgor ψ_p + osmotic ψ_s
flaccid (0 < R < rwc_tlp):  psi(R) = pi0/R                          ! ψ_p = 0
```

At `R = rwc_tlp`, `ψ_p = 0` and `ψ = pi0/rwc_tlp = pi0*eps/(pi0+eps) = psi_tlp` — consistent with
Bartlett eqn 1. **Closed-form inverse** (no Newton): flaccid `R = pi0/ψ`; turgid solves the
quadratic `eps·R² − (ψ+eps+pi0)·R + pi0 = 0`:

```
R = ( (psi+eps+pi0) + sqrt((psi+eps+pi0)**2 − 4*eps*pi0) ) / (2*eps)      ! + root ∈ [rwc_tlp, 1]
```

(discriminant `> 0` since `−4·eps·pi0 > 0`). Switch branch at `ψ = psi_tlp`.

Water content and capacitance (`W_sat_sym = (1−af)·water_sat·biomass`; apoplast is a constant
reservoir `af·water_sat·biomass`; Christoffersen et al. 2016, GMD 9:4227 adds the wood capillary
term):

```
W(psi) = W_sat_sym*R(psi) + af*water_sat*biomass
C(psi) = dW/dpsi = W_sat_sym * cR(psi),   cR = dR/dpsi:
   turgid :  cR = 1 / (eps − pi0/R**2)          ! → 1/(eps+|pi0|) at R=1
   flaccid:  cR = -R**2/pi0 = |pi0|/psi**2      ! rises sharply below TLP
```

`C(ψ)` has a real C⁰ kink at `psi_tlp` (capacitance jumps up as turgor is lost — physical); since
`C` is re-frozen each sub-step, adaptive stepping refines through it. This is the **only** PV curve
(π_o/ε/a_f are the best-constrained PV traits — Bartlett's 317-species database, π_o the dominant
driver — and the closed-form inverse makes it cheap). Cite both papers in the `meds_hydro_pv.f90`
header.

**Reproducing legacy linear-PV (X16) runs.** X16's constant capacitance is a pure-elastic reservoir
with no turgor-loss point; the nonlinear curve can't be made exactly linear over a wide range, but
matching the full-turgor capacitance `C(ψ→0) = W_sat_sym/(ε+|π_o|)` reproduces linear results over
the turgid range. Documented helper (both directions):

```fortran
pure real(wp) function pv_water_cap_from_traits(pi0, eps, af, water_sat) result(water_cap)
   water_cap = (1.0_wp - af)*water_sat / (eps + abs(pi0))         ! linear water_cap ← nonlinear
end function
! inverse: given a legacy water_cap and a Bartlett pi0, back out eps = (1-af)*water_sat/water_cap − |pi0|
```

To make a nonlinear run track an old linear one that dipped below the (nonexistent) TLP, pick `|π_o|`
so `psi_tlp` sits below the run's minimum ψ; the curves then stay close throughout.

### 6.2 Vulnerability & the Kirchhoff conductance law (the ONLY conductance path)

`plc(ψ)` is the *pointwise* retained-conductance fraction; it is used **only** through the Kirchhoff
flux potential `Φ`, **never** sampled at a single ψ to form an edge conductance. This is more
physically exact for a storage-free edge (`q = leaf_area·ΔΦ` is the exact flux, not an approximation)
and monotone in ψ ⇒ **M-matrix Jacobian** ⇒ better-conditioned near cavitation than the pointwise
`K(ψ)·Δψ` form (which can go non-monotone and break Newton). There is no harmonic-mean/pointwise edge
mode.

```fortran
elemental real(wp) function plc_retained(psi, psi50, kexp) result(f)     ! clamp ratio for psi>0
   real(wp) :: r ;  r = max(psi/psi50, 0.0_wp) ;  f = 1.0_wp/(1.0_wp + r**kexp)
end function
pure real(wp) function flux_potential(psi, kmax, psi50, kexp) result(phi) ! Φ = kmax·∫plc dψ
   ! closed form for kexp=1 (kmax·psi50·ln(1+r)) or kexp=2 (kmax·psi50·atan(r)); else fixed 5–7 pt Gauss
end function
```

**Edge flux (exact, storage-free edge):** `q = leaf_area·(Φ(ψ_up) − Φ(ψ_down))` for the matric part,
plus the additive gravity head. Two solver entry points on the *same* `Φ`:

- **`EXPM`** (frozen-linear): freeze the **secant conductance** at sub-step start,
  `K_eff = leaf_area·ΔΦ/Δψ`, so `q = K_eff·(ψ_up − ψ_down + g)` stays linear. **`Δψ→0` guard**
  (nights / equilibrated): `K_eff → k_plant_max·leaf_area·plc(½(ψ_up+ψ_down))` (the L'Hôpital limit —
  the pointwise formula survives *only* here, as a numerical fallback, not a user mode).
- **`BE`/Rosenbrock** (implicit): use the exact nonlinear `q = leaf_area·ΔΦ` with the **bounded,
  correct-sign Jacobian** `∂q/∂ψ_up = +leaf_area·k_plant_max·plc(ψ_up)`,
  `∂q/∂ψ_down = −leaf_area·k_plant_max·plc(ψ_down)` (`plc∈[0,1]` ⇒ M-matrix; no near-cavitation
  non-monotonicity). `dplc_dpsi` is then only needed for the `C`-side of the Jacobian, not the edge.

Floors on `C`/`K` go **only** in the linearization/Jacobian/reported flux, never inside `W(ψ)` or
`Φ`. Restrict `wood_kexp ∈ {1,2}` for closed-form `Φ`, or accept the fixed-quadrature cost.

### 6.3 Conservation (boundary-only, from converged ΔW)

```
dW_L = water_content(psiL_new) − water_content(psiL_old)     ! per PV curve
dW_W = water_content(psiW_new) − water_content(psiW_old)
flux%sapflow     = dW_L/dt + E                                ! W→L sap = leaf storage Δ + transp
flux%root_uptake = (dW_L + dW_W)/dt + E                       ! closes ΔW_stored = (uptake−E)·dt
```

Only `E` and `root_uptake` feed a future water budget; internal edge fluxes never do.

### 6.4 Diagnostic vertical stem ψ profile (3-node)

Reconstruct ψ along the bole (the `L–T` edge, from stem base `z_stem=0` to leaf `z_leaf=height`)
from the converged nodes + stem sapflow `Q` — purely diagnostic (state unchanged). A constant-`k`
profile would be linear in `z` but **underestimates the top tension**, because `k(ψ)=wood_kmax·
plc(ψ)` *falls* as ψ drops, steepening the gradient toward the top (top-down embolism). The same
Kirchhoff flux potential used for the edge law (§6.2) linearizes this exactly:

```
Phi(psi) = wood_kmax · ∫ plc(psi') dpsi'                                    (matric flux potential)
```

Steady mass balance ⇒ (gravity neglected) `dΦ/dz = −Q/sap_area = const`, so **Φ is linear in z for
any nonlinear k(ψ)**. Reconstruct by two-point interpolation on Φ between the bole endpoints, then
invert:

```
Phi(z) = Phi(psi_stem) + ( Phi(psi_leaf) − Phi(psi_stem) ) · (z − z_stem)/height
psi(z) = phi_inverse( Phi(z) )                     ! 1-D monotone bisection (plc>0 ⇒ Φ monotone)
```

**Gravity:** the exact balance is `dΦ/dz = −Q/sap_area − wood_kmax·plc(ψ)·grav_head`, no longer
exactly linear — either superpose the hydrostatic head `−grav_head·(z−z_stem)` on the gravity-free
Φ solution (good when the flow gradient dominates), or RK-integrate the 1-D ODE
`dψ/dz = −Q/(sap_area·wood_kmax·plc(ψ)) − grav_head` with a handful of steps (exact, cheap for a
diagnostic).

```fortran
pure real(wp) function stem_psi_at_height(z, psi_stem, psi_leaf, height, ...)   ! Φ-interp + invert
```

**One Φ kernel, everywhere:** (1) the dynamic edge-conductance law (§6.2), (2) this diagnostic
profile, and (3) the hydraulic supply function `E(ψ)=∫k dψ` (Sperry) for a future gain-risk stomatal
optimization. Uses: locate `z*` where `psi(z*)=wood_psi50` (bole embolism height, tall trees fail
top-down); height-resolved ψ for a future multi-layer canopy coupling; a continuous profile from the
lumped node. Assumes no taper.

---

## 7. State additions (cohort SoA)

Add one 2-D prognostic field to `cohort_block` (pattern: `growth_hist(nwin,cap)`):

```fortran
real(wp), allocatable :: psi_node(:,:)     ! (N_HYDRO, cap) [MPa] node water potentials, PROGNOSTIC
```

Touch the **7 centralized lockstep routines** in `meds_demography_types.f90`, in lockstep:
`cohort_alloc` (allocate + init `0.0_wp`), `site_free`, `cohort_ensure_capacity`
(`tmp%psi_node(:,1:m)=…`), `move_alloc_block` (`move_alloc`), `cohort_reorder`
(`psi_node(:,1:m)=psi_node(:,perm(1:m))`), `cohort_compact` (via reorder), `copy_cohort_slot`. Not
geometry-derived ⇒ **not** in `set_cohort_size`. **Initialize `= soil_psi` at every creation site**
(`add_cohort`/`init_bare_ground`, `apply_recruitment`, `split_cohorts`, `apply_patch_disturbance`) —
outside the 7 routines; grep `assign_cohort_id` sites as the checklist. **Fusion/fission** blends
node ψ by conserved water mass (`Σ nplant·W(ψ)` → invert to survivor ψ via `meds_hydro_pv`), at the
`meds_demography_structure` sites. **Biomass-change conservation:** within the kernel biomass is a
frozen BC (no issue); the fixed-ψ/changed-`C` `ΔW` from growth/allocation is booked as a water
source/sink at the demography/fast-loop boundary (future coupling).

Future GPU sweep: add `p_hydro_*` gathered arrays in `gather_pft_params` (not needed for the
standalone seam).

---

## 8. PFT traits

Add to `pft_table_t`, allocate in `alloc_pft_table`, load via `req_pa` presence-map, one row each in
`meds_config_pft.toml`:

| Field | Units | Notes |
|-------|-------|-------|
| `k_plant_max(:)` | kg s⁻¹ MPa⁻¹ m⁻²_leaf | **leaf-area-specific whole-plant conductance** (default `HYDRO_COND_KPLANT`; measurable, = Bartlett `K_plant`) |
| `wood_psi50(:)`, `wood_kexp(:)` | MPa(<0), – | xylem vulnerability (modulates `k_plant_max`) |
| `wood_kmax(:)`, `vessel_curl(:)` | kg m⁻¹ s⁻¹ MPa⁻¹, – | sapwood conductivity + path-length — only for `HYDRO_COND_SEGMENT` |
| `leaf_water_sat(:)`, `wood_water_sat(:)` | kg H₂O kgC⁻¹ | saturated water / biomass |
| `leaf_pi0(:)`, `leaf_eps(:)`, `leaf_af(:)` | MPa(<0), MPa(>0), – | leaf PV (Bartlett 2012) |
| `wood_pi0(:)`, `wood_eps(:)`, `wood_af(:)` | MPa(<0), MPa(>0), – | wood PV (+ capillary, Christoffersen 2016) |
| `root_sra(:)`, `root_beta(:)` | m² kgC⁻¹, – | rhizosphere test-helper / future per-layer roots (§16) |

(`leaf_water_cap`/`wood_water_cap` are dropped — no linear PV path; legacy `water_cap` maps to
`pi0,eps` via `pv_water_cap_from_traits`, §6.1. For 3-node, add resistance partition fractions
`f_leaf,f_stem,f_root` to split `k_plant_max`.)

`derive_hydro_params(pft)` (new, beside `derive_leaf_params`): compute `psi_tlp`, `rwc_tlp` per
tissue from `pi0,eps` (Bartlett eqns 1–2); precompute `W_sat_sym`, per-PFT `psi_min` bounds. Pure
trait-table transform, no state. `validate_config`: `wood_psi50<0`, `wood_kexp>0`, `pi0<0`, `eps>0`,
`0≤af<1`, and `psi_tlp` above `wstress_psi_close`. Example `[pft]` rows (pioneer/mid/climax):

```toml
k_plant_max = [ 6.0e-4, 4.0e-4, 2.5e-4]  # [kg/s/MPa/m2_leaf] whole-plant conductance (default cond mode)
wood_psi50  = [-1.5, -2.5, -3.5]   # [MPa] climax most cavitation-resistant
wood_kexp   = [ 2.0,  3.0,  4.0]
leaf_pi0    = [-1.2, -1.8, -2.2]   # [MPa] osmotic at full turgor  (Bartlett 2012)
leaf_eps    = [ 8.0, 15.0, 25.0]   # [MPa] elastic modulus; note the eps–pi0 negative covariance
leaf_af     = [ 0.30, 0.25, 0.20]
# wood_kmax / vessel_curl only needed for HYDRO_COND_SEGMENT
```

---

## 9. Config `[hydraulics]` block

```fortran
integer(ik), parameter :: HYDRO_NODES_2 = 2_ik, HYDRO_NODES_3 = 3_ik           ! default _2
integer(ik), parameter :: HYDRO_SOLVER_EXPM = 1_ik, HYDRO_SOLVER_BE = 2_ik      ! default _EXPM
integer(ik), parameter :: HYDRO_COND_KPLANT = 1_ik, HYDRO_COND_SEGMENT = 2_ik   ! default _KPLANT (whole-plant k)
integer(ik), parameter :: HYDRO_SUBSTEP_ADAPTIVE = 1_ik, HYDRO_SUBSTEP_FIXED = 2_ik
```

(No `HYDRO_PV_*` — nonlinear Bartlett curve is the only PV, §6.1. No `HYDRO_EDGE_*` — the Kirchhoff
flux potential is the only conductance law, §6.2.) `meds_config_t` scalars (from `[hydraulics]` in
`meds_config_main.toml`): `hydro_enable` (0=off, ψ saturated), `hydro_topology` (default 2),
`hydro_solver` (default EXPM), `hydro_cond_mode` (default whole-plant K), `hydro_substep_mode`,
`hydro_gravity_on`, `hydro_rtol`, `hydro_atol`, `hydro_h_init`, `hydro_max_substep`.

---

## 10. Coupling contract with the leaf module

Leaf and hydraulics are **duals**; neither owns the fixed point.

```
leaf:        psi_leaf ─► beta = clamp((psi_leaf−psi_close)/(psi_open−psi_close),0,1) ─► E
hydraulics:  E (+ soil BC) ─► advance psi over dt ─► psi_leaf = psi(NODE_LEAF)
```

* **In:** `hydro_env_t%transp = E`; caller bridges units `E_kg = E[mol m⁻² s⁻¹]·leaf_area·0.018`.
* **Out:** `flux%psi_leaf = psi(NODE_LEAF)`, consistent with the `E` used.
* **Who iterates:** the *future* fast loop, not this kernel — it freezes `E` over `dt` (an
  explicit-coupling lag the loop's 1–2 Picard passes close; DTLSM is short). Both modules stay
  sealed and `meds_shared`-only. **No fast loop exists today** — hydraulics ships standalone
  (CTest + optional Python), as `leaf` and `biophys` did before wiring.

---

## 11. GPU / nvfortran

* Kernels `pure`/`elemental`, arithmetic + intrinsics only; flat params; fixed `N_HYDRO` ⇒ no
  allocatables/runtime shapes. **Branch-light:** X16 `zero_flow`/`is_small` traps → smooth limiters
  (Kirchhoff edge, `C`/`K` floors, `plc` clamp, `sinhc` `if`); only adaptive accept/reject and
  `n_sub` are data-dependent (mitigated by `HYDRO_SUBSTEP_FIXED`).
* **Array-valued-function-result trap:** the 2-node kernel is scalar (trap-free); the 3-node array
  path binds every array result to a named array before any call.
* **FPE-safe in Debug** (`-fpe0`/`-Ktrap=fp`): the `sinhc` guard, the `plc` clamp, the `Δψ→0`
  conductance guard, and the PV closed-form inverse all avoid `0/0` and `(neg)**real`.
* No `site_t` on the device. A green **ifx** suite is not sufficient — build **nvfortran multicore**
  on the new module.

---

## 12. Validation & milestones

**CTest `test_plant_hydraulics`:**

1. **Single-cohort diurnal** — sinusoidal `E`, well-watered soil: `psi_leaf` tracks `E`, midday
   minimum, nocturnal recovery toward `soil_psi`; bounded `≤ 0`.
2. **Conservation** — `E=0 ⇒ dψ=0` to `atol`; with `E`, `ΔW_stored = (uptake−E)·dt` to machine
   precision (ΔW-based fluxes), both topologies and both conductance modes.
3. **PV curve vs Bartlett** — nonlinear `ψ(R)` reproduces `psi_tlp`/`rwc_tlp` (Bartlett eqns 1–2);
   `water_content(rwc_from_psi(ψ)) == ψ` round-trip (closed-form inverse) to round-off; `C(ψ)`
   matches finite-difference `dW/dψ`; the `C`-jump sits exactly at `psi_tlp`.
4. **vs fine explicit reference** — `EXPM`/`BE` adaptive matches a 1000-substep explicit RK4 to
   `rtol` in well-watered **and** near-`ψ50` regimes; sub-step count drops when wet; the coupled
   scheme beats an X16-style operator split under `C_L~C_W`. This reference — not X16 — is the
   accuracy ground truth.
5. **Qualitative vs X16** — the diurnal `psi_leaf` trajectory shows the same midday depression /
   nocturnal recovery as ED2-hydro (pattern check, not bit-for-bit).
6. **Diagnostic profile** — 3-node `stem_psi_at_height` (Kirchhoff) is monotone in `z`, reduces to
   hydrostatic at `Q=0`, and its `ψ=wood_psi50` crossing height moves up-canopy as drought deepens.
7. **FPE + degenerate** — `Δ→0`, nightly `ψ≈0`, leafless `bleaf→0`: no NaN/trap, continuous limits.
8. **nvfortran multicore build green** (portability gate).

**Phased milestones:**

| Phase | Deliverable | Tests |
|-------|-------------|-------|
| **P0** | Types + CMake lib + seam skeleton; `meds_hydro_pv` nonlinear (Bartlett/Christoffersen) + closed-form inverse + legacy-linear mapping; `flux_potential`/`phi_inverse` Φ kernels | 3 |
| **P1** | **Production default:** 2-node `EXPM` + whole-plant `k_plant_max` + Kirchhoff edge + adaptive sub-stepping + ΔW-based fluxes | 1,2,4,5,7 |
| **P2** | *(opt-in)* 3-node `BE` (Thomas) + `k_plant_max` partition + diagnostic `stem_psi_at_height` | 1,2,4,6 (3-node) |
| **P3** | SoA `psi_node` + 7 lockstep routines + creation-site init + fusion water-mass blend; TOML traits + `derive_hydro_params` + `validate_config` | full engine build |
| **P4** | C-API + Python pkg; nvfortran GPU parity; `HYDRO_SUBSTEP_FIXED` | 8 |
| **P5** | *(future)* DTLSM fast loop; leaf↔hydraulics Picard; demography `ΔW` booking | — |

---

## 13. References

- Xu, Medvigy, Powers, Becknell, Guan (2016) *New Phytol.* 212:80 — ED2 X16 plant hydrodynamics
  (physical structure referenced; numerics superseded).
- Bartlett, Scoffoni, Sack (2012) *Ecol. Lett.* 15:393 — leaf turgor-loss-point / p–v theory;
  `psi_tlp = π_o ε/(π_o+ε)`, `RWC_tlp = (π_o+ε)/ε`; π_o dominant, ε–π_o covary negatively.
- Christoffersen et al. (2016) *Geosci. Model Dev.* 9:4227 — FATES-HYDRO / TFS v.1-Hydro;
  apoplast/capillary reservoir, tissue PV partitioning, per-layer rhizosphere.
- Tyree & Hammel (1972) *J. Exp. Bot.* 23:267 — the symplast turgor+osmotic p–v model.

---

## 14. Open questions

1. **`k_plant_max` values & vulnerability sharing** — the leaf-area-specific whole-plant conductance
   per PFT, and whether one `wood_psi50/kexp` governs the whole plant in 2-node (yes for now).
2. **Hydraulic redistribution** — bidirectional root/soil flow (branch-free) vs a smooth ramp.
   *Recommend:* bidirectional (also the natural behavior of the per-layer extension, §16).
3. **Fusion ψ blend** — water-mass-conserving vs reset to `soil_psi`. *Recommend:* water-mass-conserving.
4. **`hydro_atol/rtol`** defaults (`atol=1e-3` MPa, `rtol=1e-3`?) and `hydro_max_substep` (200?).
5. *(deferred to the 3-node upgrade)* root-vs-stem trait split, stem node height, gravity treatment
   in the ψ(z) diagnostic.

---

## 15. How this differs from X16 and FATES-HYDRO

MEDS sits deliberately **between** X16 (cheap, analytic, lumped) and FATES-HYDRO (rich, fully
numerical, high-dimensional): it keeps X16's low dimension and per-step analytic exactness, adopts
FATES's nonlinear PV and full node coupling, and adds two things neither has — a **Kirchhoff-
integrated conductance** and a **stateless, GPU-ready** kernel.

| Aspect | **X16** (ED2, Xu 2016) | **FATES-HYDRO** (Christoffersen 2016) | **MEDS** (this design) |
|--------|------------------------|----------------------------------------|------------------------|
| Nodes | 2 pools (leaf, stem+root lumped) | ~10+: rhizosphere shells + absorbing root (per layer) + transporting root + N stem + leaf | fixed small: 2 (default) / 3; per-layer roots planned (§16) |
| Pressure–volume | linear, constant capacitance | nonlinear (elastic + capillary + apoplast) | nonlinear Bartlett/Tyree–Hammel + Christoffersen apoplast; **closed-form + closed-form inverse** |
| Vulnerability → conductance | PLC **pointwise** at ψ_wood, frozen over `dtlsm` | PLC pointwise per segment | **Kirchhoff flux potential Φ=∫k dψ** — integrated over the potential drop; monotone / M-matrix |
| Conductance parameterization | stem conductivity × sapwood area / length | per-segment conductivities | **whole-plant, leaf-area-specific** `k_plant_max` (measurable; Bartlett Ohm's law); segment mode optional |
| Node coupling | **operator-split** leaf-then-wood (1st-order) | fully coupled (implicit) | **fully coupled**, no split |
| Time integration | one analytic step per `dtlsm`, no sub-cycling | implicit, internal iterations | frozen-coefficient **matrix-exponential** (exact per sub-step) or L-stable BE; **adaptive sub-stepping** |
| Degenerate cases | ~5 `zero_flow` traps + `is_small` lumping | iteration safeguards | **smooth limiters / floors** (branch-light, no case traps) |
| Soil–root interface | Katul aggregate; HR disabled | explicit rhizosphere shells per layer; HR emergent | given BC now; per-layer roots + emergent HR planned (§16) |
| Architecture | stateful, in ED2's RK4 tracers | stateful; `*Mem`/compute split | **stateless** kernel; ψ in cohort SoA (`*Mem`/compute split) |
| Units / platform | meters of head; CPU | SI; CPU | MPa; **CPU + GPU** (`pure`/`elemental`, OpenMP target) |

**In one line:** MEDS = X16's small, analytic, lumped model — but *coupled* instead of split,
*nonlinear-PV* instead of linear, *Kirchhoff-integrated* conductance instead of pointwise,
*adaptively sub-cycled*, and *stateless / GPU-ready* — approaching FATES's physical fidelity at a
fraction of the dimensionality and cost.

---

## 16. Future extension: per-soil-layer fine-root water nodes (compatible — no core rework)

**Planned (not now):** replace the single aggregated `(soil_psi, rhizo_cond)` BC with a fine-root
water node `R_k` in **each soil layer the roots reach** — each with its own capacitance (fine-root
storage in that layer), connected to the common stem/collar node and to its soil layer. Topology
becomes a chain with a **fan at the collar**: `Leaf – Stem/collar – {R_1 … R_nlayer}`, each `R_k –
soil layer k`. This aligns MEDS's below-ground structure with FATES-HYDRO's per-layer absorbing
roots.

This design was built to accommodate it; confirming compatibility point by point:

1. **Governing physics unchanged.** Eq. (1) is written for an *arbitrary* node graph (`Σ_{j~i}`);
   adding nodes/edges only extends the sums. The Kirchhoff conductance law (§6.2) and the nonlinear
   PV (§6.1) apply per edge / per node exactly as-is.
2. **State field generalizes directly.** `psi_node(N_HYDRO, cap)` already carries per-node ψ; bump
   the compile-time `N_HYDRO` to `2 + max_root_layers` (still a fixed shape ⇒ GPU-friendly). The 7
   lockstep routines (§7) handle the 2-D field generically — **no new bookkeeping**.
3. **Solver structure unchanged.** Still `C ψ' = A ψ + b`; only the sparsity changes from tridiagonal
   (chain) to **arrow / bordered** (collar coupled to all `R_k` in parallel). `A` stays a symmetric
   negative-definite grounded Laplacian, so the eigen/stability guarantees (§4) hold, and the arrow
   system still solves in **O(N)** by eliminating the leaf and the per-layer root nodes into the
   collar (Schur complement). `BE` (§5.2) becomes the natural solver; the matrix-exp generalizes via
   the symmetric-eigen route if wanted.
4. **Two capabilities the current single node cannot give — for free.** (a) **Per-layer soil-water
   uptake** `rhizo_cond_k·(soil_psi_k − ψ_Rk)` (from converged ψ), the sink a future soil module
   needs; (b) **emergent hydraulic redistribution** — a root in a dry layer with `soil_psi_k < ψ_Rk`
   effluxes water through the collar to wetter layers, purely from the sign of the flux (bidirectional,
   no HR switch).
5. **Conductance partitioning is already specified.** The collar→`R_k` edges partition the below-
   ground whole-plant conductance by layer via the `root_beta` root-fraction profile (§3) — the same
   resistance-fraction mechanism used for the 3-node split (§4).

**Localized changes when it lands:** `N_HYDRO` max bump; `hydro_env_t` carries per-layer
`soil_psi(:)`/`rhizo_cond(:)` arrays; `hydro_flux_t` returns per-layer `root_uptake(:)`; the A-matrix
assembly handles the arrow structure (`BE` default); conservation (§6.3) sums the boundary terms over
layers (structure unchanged). Nothing in the 2-/3-node core blocks it — the equation, the Kirchhoff
conductance, the nonlinear PV, the stateless contract, and the SoA state field all extend directly.
