# biophys: canopy radiative transfer — design & implementation plan

A faithful reimplementation of ED2's **two-stream** multi-layer canopy radiative transfer (the
`icanrad=2` solver in `../ED2/ED/src/dynamics/twostream_rad.f90`; Liou 2002 ch.6, Longo et al. 2019),
re-expressed in the MEDS idiom — a **sealed, mostly-stateless, array-interface calculator** in the mould
of `src/plant/leaf/`. Two deliberate departures from a literal ED2 port, decided with the project owner:

1. **One unified multi-band solver** (default bands VIS / NIR / LW) instead of separate SW and LW
   routines. Every band carries a thermal-emission source that is simply **zero for VIS/NIR**.
2. **SCOPE / 4SAIL leaf-angle scattering coefficients** (van der Tol 2009; Verhoef 1984; exact
   leaf-inclination integrals) in place of ED2/CLM's `φ₁/φ₂/μ̄` + `(1+χ_L)²` backscatter approximation,
   which breaks for near-horizontal leaves.

This document is the plan; the code does not exist yet.

---

## 1. Scope, target scheme, and reachability facts

ED2 ships three canopy-RT schemes selected by `icanrad`, plus a dead one:

| `icanrad` | routine | status |
|-----------|---------|--------|
| 0 | `old_sw/lw_two_stream` | **dead** — `ed_opspec.F90:2563` rejects `icanrad<1` |
| 1 | `sw/lw_multiple_scatter` | Zhao & Qualls multiple-scattering |
| 2 | `sw/lw_two_stream` (`twostream_rad.f90`) | **the two-stream MEDS targets** (Longo 2019) |

Two constraints from `ed_opspec.F90` shape the port: `icanrad` is namelist-restricted to `{1,2}` (the
"old" two-stream is unreachable), and `icanrad=2` **forbids `crown_mod=1`**, so under the two-stream the
crown-area index `cai ≡ 1`. We therefore implement the **`cai=1` (closed-crown) path** — this is ED2's
own supported configuration, not a MEDS simplification. Finite crown area is a documented non-goal (§10).

What we keep from ED2: the **multi-layer vertical structure** (one radiation layer per cohort, a
`2N+2` interface-continuity two-stream solve). What we replace: the SW/LW duplication (→ one band
solver, §3) and the leaf-angle optics (→ SCOPE, §6).

---

## 2. Where it lives (library DAG)

Mirroring the leaf precedent, the **solver core is `meds_state`-free**: it takes plain arrays
(`pft(:)`, `lai(:)`, `wai(:)`, per-cohort optics) plus a forcing struct — the RT analogue of
`leaf_env_t`. A thin orchestration layer reads `site_t` and calls the core per patch. This keeps the
core unit-testable in isolation, device-eligible, and reusable.

```
shared ─┬─ allometry ─ state ─┐
        │                     ├─ biophys(core: shared-only)  ← band solver + optics, array interface
        └─ pft_params ────────┘
                                    biophys(orchestration) reads state, walks patches ─────────────┘
```

`add_library(meds_biophys ...)`; the solver-core TUs link **`meds_shared` only**; the one orchestration
TU also links `meds_state`. CMake globs `src/biophys/*.f90`.

| file | module | role |
|------|--------|------|
| `meds_rad_types.f90` | `meds_rad_types` | data-only: forcing, band spec, per-cohort optics, outputs |
| `meds_leaf_angle.f90` | `meds_leaf_angle` | **pure** SCOPE LIDF + `volscatt` → `bf`, `G(μ)` (§6) |
| `meds_canopy_optics.f90` | `meds_canopy_optics` | **pure** ρ,τ + `bf`,`G` → ω, β, β0; leaf/wood blend (§6) |
| `meds_surface_optics.f90` | `meds_surface_optics` | **pure** ground reflectance + emission (bare-soil stub, §5) |
| `meds_rad_linsys.f90` | `meds_rad_linsys` | banded / block-tridiagonal solve of the 2N+2 system (§8) |
| `meds_twostream_band.f90` | `meds_twostream_band` | **the unified single-band solver** (§4) |
| `meds_canopy_radiation.f90` | `meds_canopy_radiation` | **THE public seam** + per-patch orchestration over `site_t` |

`meds_canopy_radiation` is the analogue of `meds_leaf_physiology` — the one module production callers
`use`. Everything else is internal.

---

## 3. The unified multi-band idea (feedback #1)

The diffuse two-stream operator is **identical across bands**; only source terms and boundary conditions
differ. So we solve one band at a time with a `band_spec` and loop over a **config-defined band list**
(default VIS, NIR, LW):

| band | direct beam? | emission? | ω (per-cohort) | top BC | bottom BC |
|------|:---:|:---:|---|---|---|
| VIS  | yes | (0) | ρ+τ ≈ 0.15 | incident diffuse sky | ground albedo·refl + beam refl |
| NIR  | yes | (0) | ≈ 0.8 | incident diffuse sky | ground albedo·refl + beam refl |
| LW   | no  | yes | 1−ε ≈ 0.05 | atmospheric downwelling `rlong` | (1−ε_g)·refl + ε_g·σT_g⁴ |

The emission source is **universal but identically zero for VIS/NIR** — exactly the requested framing.
This is validated by SCOPE's own architecture: its thermal module `RTMt_sb` reuses the *same* diffuse
scattering coefficients (`sigb, sigf, att, m`) as the solar module `RTMo`, adding only a per-layer
Planck source `(1−ρ−τ)·σT⁴`.

Consequences (both simplifications):
- **Work in absolute W/m² per band**, not ED2's normalize-to-1 fractions. Broadband albedo is then just
  `up(top)/incident`, and ED2's three-regime renormalization dance disappears.
- **Day/night gating is automatic**: a band with zero incident *and* zero emission is skipped (VIS/NIR
  at night); LW always runs (it always has an emission source).
- Kills ED2's ~600-line SW/LW near-duplication → one `meds_twostream_band`. (MEDS "never reproduce a big
  file" rule.)

---

## 4. The public seam & data types

**Core (shared-only, array interface):**
```fortran
subroutine solve_twostream_band(ncoh, etai, band, cohort_optics, forcing, out)
   integer(ik),          intent(in)  :: ncoh
   real(wp),             intent(in)  :: etai(ncoh)        ! effective area index per cohort (bottom→top)
   type(band_spec_t),    intent(in)  :: band              ! has_beam/has_emission, incident, ground optics
   type(band_optics_t),  intent(in)  :: cohort_optics     ! per-cohort ω, β, β0, k, emission for THIS band
   type(rad_forcing_t),  intent(in)  :: forcing           ! cosz, band incident fluxes
   type(band_out_t),     intent(out) :: out               ! per-cohort absorbed W/m2 + boundary fluxes
end subroutine
```

**Orchestration seam over the demographic state:**
```fortran
subroutine canopy_radiation(site, cfg, forcing, surf, rad)   ! walks patches; fills a diagnostic field
```
`rad` is the per-cohort diagnostic SoA the biophys README promises ("absorbed-PAR field consumed by
`plant/leaf`") — laid out in lockstep with the flat cohort block, **not** stored in `cohort_block` (RT
is a re-derivable fast flux, not prognostic state — keep it out of the fusion/reorder lockstep).

```fortran
type :: rad_forcing_t                       ! the "leaf_env_t" of canopy RT — explicit, no globals
   real(wp) :: cosz                          ! cosine solar zenith / angle of incidence (floored)
   real(wp) :: par_beam, par_diff, nir_beam, nir_diff   ! [W/m2] incident, absolute
   real(wp) :: rlong                         ! [W/m2] downwelling atmospheric LW
   real(wp) :: canopy_temp                   ! [K] leaf/wood temperature (until an energy balance exists)
end type

type :: band_spec_t
   logical  :: has_beam, has_emission
   real(wp) :: incident_beam, incident_diff  ! [W/m2] top-of-canopy for this band
   real(wp) :: ground_refl                   ! albedo (SW) or 1-emiss (LW)
   real(wp) :: ground_emission               ! [W/m2] ε_g·σT_g⁴ (LW); 0 for SW
end type

type :: band_optics_t                        ! per-cohort, for one band (gathered — no PFT-table indexing in kernels)
   real(wp), allocatable :: omega(:)         ! single-scattering albedo ω = ρ+τ (leaf/wood blended)
   real(wp), allocatable :: beta(:)          ! diffuse backscatter fraction  (SCOPE)
   real(wp), allocatable :: beta0(:)         ! direct-beam upscatter fraction (SCOPE; beam bands only)
   real(wp), allocatable :: kdir(:)          ! direct-beam extinction G(μ)/μ (beam bands only)
   real(wp), allocatable :: emission(:)      ! [W/m2] per-layer (1-ω)σT⁴ (emission bands only; else 0)
   real(wp), allocatable :: leaf_frac(:)     ! (1-ω_leaf)·elai / [(1-ω_leaf)·elai + (1-ω_wood)·ewai]  (leaf/wood split)
end type
```

---

## 5. Surface reflectance (feedback #2)

A `meds_surface_optics` module with a **pure** `surface_band_optics(surf, band) → (reflectance,
emission)`. `surf` is a `surface_state_t` stub carrying (eventually) soil moisture, litter, standing
water, snow — but **for now returns bare soil**:

```
SW bands: reflectance = soil_albedo_band                       (configured; later f(moisture,colour))
LW band : reflectance = 1 - soil_emiss ;  emission = soil_emiss·σ·T_soil⁴
```
This is exactly where ED2's bug **B6** lives (`radiate_driver.f90:639`, the `albedo_damp_nir` self-
assignment) — write it correctly from the start. It slots in as the bottom boundary condition of the
band solver. The full version (soil-moisture/colour albedo, snow, surface water, the McCumber–Pielke
peat branch) lands when soil state exists.

---

## 6. Leaf-angle scattering — the SCOPE scheme (feedback #4)

Replace ED2/CLM's `φ₁, φ₂, μ̄` extinction and the `(1+χ_L)²` backscatter with the **exact
leaf-inclination integrals** of SCOPE/4SAIL. This fixes both CLM approximations (linear G-function
*and* the backscatter surrogate) with one internally-consistent leaf-angle distribution (LIDF).

### 6a. LIDF — a two-parameter Beta distribution (feedback #2)
Represent the leaf-inclination distribution as a **Beta(p, q) distribution** on the normalized
inclination `t = θ_leaf/(π/2) ∈ [0,1]` (Goel & Strebel 1984, *Agron. J.* — the standard statistical LAD
model). This is a single, generic, continuous two-parameter family that **subsumes both Campbell's
ellipsoidal (a one-parameter curve in Beta space) and Verhoef's archetypes**: planophile, erectophile,
plagiophile, extremophile, uniform (Beta(1,1)) and near-spherical are all specific `(p,q)` — Goel &
Strebel tabulate them. It is the most flexible and future-proof choice, and differentiable (nice for
future calibration/data assimilation). Expose the two PFT traits as **mean leaf angle + concentration**
(more interpretable, the successor to ED2's single `orient_factor`), mapped internally to `(p,q)` via
`p = t̄(t̄(1−t̄)/σ_t² − 1)`, `q = (1−t̄)(…)`.

Discretize onto the SCOPE 13-class grid `litab = [5,15,…,75, 81,83,85,87,89]°` (fine near vertical);
the class weights are `lidf(i) = I_{t_{i+1}}(p,q) − I_{t_i}(p,q)` where `I_x(p,q)` is the regularized
incomplete beta function (a ~30-line pure `betai`/`betacf` continued-fraction helper — the one special
function to add). Everything downstream (`bf`, `G(μ)`, `β`, `β0`) is unchanged — only the `lidf`
generator differs.

### 6b. Angle integrals
```
bf   = Σ_i cos²(θ_i)·lidf_i                       ! 2nd moment of leaf inclination   [μ-INDEPENDENT → precompute]
G(μ) = Σ_i chi_s(θ_i, μ)·lidf_i                   ! exact Ross G-function via volscatt [μ-dependent → per step]
k    = G(μ)/μ                                      ! direct-beam extinction coefficient
```
`chi_s` is the azimuth-averaged leaf-area projection toward the sun (`volscatt`, Verhoef). **MEDS needs
fluxes only, not directional reflectance**, so the view/bidirectional half of `volscatt` (`chi_o, frho,
ftau, sob…`) is *not needed* — only `chi_s` (for `G`) and `bf`.

### 6c. Per-band scattering coefficients (per tissue, then blended)
For band single-tissue reflectance ρ and transmittance τ:
```
ω  = ρ + τ                                                     ! single-scattering albedo
β  = ½·[ 1 + bf·(ρ−τ)/ω ]        (= sigb/ω)                     ! DIFFUSE backscatter fraction
β0 = ½·[ 1 + (bf/k)·(ρ−τ)/ω ]    (= sb/(k·ω))                   ! DIRECT-beam upscatter fraction
```
> **Note the `β0` normalizer is `k·ω`, not `ω`.** `sb = sdb·ρ + sdf·τ` with `sdb,sdf = ½(k±bf)`, and the
> total direct single-scatter is `sb+sf = k·ω`, so `β0 = sb/(kω)`. Dividing by ω alone is a known slip.

**Leaf/wood blend** (feedback #3): each area is clumping-corrected — `elai = clumping_leaf·lai`,
**`ewai = clumping_wood·wai`** — and `etai = elai + ewai`. Per cohort the area weights are `leaf_weight =
elai/etai`, `wood_weight = ewai/etai`, and each optical coefficient is the area-weighted mix of the leaf
value (leaf ρ,τ) and the wood value (wood ρ,τ), both using the same `bf`, `k`. (`clumping_wood` is a
separate PFT trait; stems are more clumped/vertical than leaves — a separate wood LIDF is a later
refinement.)

### 6d. Feeding the two-stream (Sellers structure retained, μ̄ ≡ 1)
The SCOPE diffuse eigenstructure `att = 1−sigf`, `m = √(att²−sigb²)` is *algebraically identical* to
Sellers' `λ = √((1−εω)(1−ω))/μ̄` with **`μ̄ ≡ 1`** (one can show `att−sigb = 1−ω`, `att+sigb = 1−εω`
with `ε = 1−2β`). So we keep ED2's Sellers solver form verbatim and just set:
```
μ̄  = 1                                   (SAIL convention: diffuse geometric extinction = 1 per unit area)
λ  = m = √((1−ω)(1−εω)) ,  ε = 1−2β
γ± = ½·(1 ± √((1−ω)/(1−εω)))
μ0 = 1/k ,  expm0 = exp(−k·etai)          (direct-beam transmittance; k is the EXACT G(μ)/μ)
expl± = exp(±λ·etai)
```
This removes the `μ̄` machinery (and its extreme-angle singularities) entirely — leaf angle now enters
only through `bf` (scattering asymmetry) and `G(μ)` (direct extinction), the SCOPE way. It differs
slightly from ED2/CLM's angle-dependent diffuse `μ̄`; that is the intended SCOPE-consistency trade.

### 6e. Cost split
`bf`, `ω`, `β`, `λ`, `γ±` are **μ-independent → precompute once per (PFT, band)** at setup. Only `G(μ)`,
`k`, `β0` depend on the sun angle → recompute per step per PFT (a 13-element dot product; cheap).

### 6f. Why CLM breaks (the bug feedback #4 flags)
CLM's `β = ½[1 + ((1+χ_L)/2)²·(ρ−τ)/ω]` substitutes `((1+χ_L)/2)²` for the exact `bf`. With `χ_L`
clamped to `[−0.4, 0.6]`, the surrogate maxes at `0.64 < 1`, so **planophile backscatter is structurally
under-estimated** (it can never reach the horizontal-leaf limit `β = ρ/(ρ+τ)`), it is already wrong at
spherical (`0.25` vs the exact `bf = 1/3`), and it cannot represent bimodal LADs at all. SCOPE evaluates
`bf` exactly by quadrature, recovering every limit.

---

## 7. Wood area & the leaf/wood split (feedback #3)

`waiᵢ = wai_lai_ratio · laiᵢ` (config'd, default 0.1; from a wood/branch-area allometry later). Both
areas are clumping-corrected (feedback #3): `elai = clumping_leaf·lai`, **`ewai = clumping_wood·wai`**,
`etaiᵢ = elaiᵢ + ewaiᵢ`. This gives per-cohort leaf+wood-blended optics (§6c) and the **leaf/wood
absorption split**: the fraction of a cohort's per-band absorbed radiation that lands on leaves is
```
leaf_frac_band = (1−ω_leaf)·elai / [ (1−ω_leaf)·elai + (1−ω_wood)·ewai ]
```
derived from the **same** clumping-corrected areas `(elai, ewai)` and absorptivity `(1−ω)` the solver used. This *is*
the fix for ED2 bug **B7** (their TIR split dropped `clumping`, diverging from the solver's own
weighting) — here it is one formula, consistent across all bands. Needs wood optical traits
(`wood_reflect/trans_vis/nir`, `wood_emiss_tir`) in the PFT table. Wood temperature for the LW emission
aliases `canopy_temp` until a wood energy balance exists (`par_l = par_v · leaf_frac_vis` for
photosynthesis).

---

## 8. Solving the linear system (feedback #5, benchmarked)

**Not LINPACK** (obsolete, superseded by LAPACK), and **not a library dependency in the core.** This was
checked empirically, not asserted: a micro-benchmark of the exact pattern — 100,000 independent
variable-size (N∈[5,80]) block-tridiagonal 2×2 systems (= the two-stream 2N+2 system, one per patch),
one source compiled serial / `-mp` / `-mp=gpu`, on the actual RTX 3050 Ti + 8-core CPU (nvfortran 25.11).
Time to solve all 100k systems (residual ~1e-16 everywhere):

| storage layout | serial | multicore (8) | GPU |
|---|---|---|---|
| **contiguous** (per-patch) | 0.174 s | **0.059 s** | 0.110 s |
| **interleaved** (padded to Nmax) | 0.649 s | 0.363 s | **0.045 s** |

Findings:
1. **The custom in-kernel block-Thomas/adding solver is correct (residual 1e-16) and offloads across all
   three backends from ONE source** — confirmed. It is O(N) vs ED2's dense Gauss `lisys_solver8` at O(N³).
2. **CPU and GPU want opposite memory layouts.** CPU is fastest with contiguous per-patch storage
   (cache-friendly); GPU is fastest with interleaved+padded storage (coalesced). Best CPU (0.059 s,
   multicore/contiguous) and best GPU (0.045 s, interleaved) are within ~1.3× on this laptop — the GPU is
   **not** a decisive win for this small, sequential-per-system, memory-bound pattern.
3. **The solve is cheap.** 100k patches in ~60 ms on 8 CPU cores → a realistic single-site MEDS run
   (~10²–10³ patches) is sub-millisecond. GPU offload only pays off at very large batch (regional/global,
   10⁶+ patches); below that it is kernel-launch-latency-bound and loses to the CPU.
4. **BLAS/LAPACK doesn't help and costs more.** It cannot be called inside an OpenMP `target` region
   (device-side), so it can't serve the GPU path at all; per-system `dgbsv` on CPU carries call/pivot
   overhead for tiny systems and forks a second code path. cuSPARSE's batched penta-diagonal
   `gpsvInterleavedBatch` needs the *same* interleaved+padded, uniform-size layout the custom GPU kernel
   uses **and** is CUDA-only + host-launched — i.e. the custom interleaved kernel already does what it
   would, portably and dependency-free.

**Plan:**
- **Custom banded/adding solver in `meds_rad_linsys`, self-contained** (no external dependency, preserving
  MEDS's dependency-free engine and the single serial/multicore/GPU source).
- **Primary path = CPU (serial + multicore) with contiguous per-patch layout** — this is where RT runs
  for typical MEDS runs, and its cost is negligible.
- **GPU offload = secondary**, worthwhile only for very large batches; it needs the interleaved layout
  (or a transpose before the kernel). Keep the solver kernel **layout-parameterizable** so this switch is
  a local change, not a rewrite.
- Use the transmittance-scaled **adding form** (`T = exp(−λ·etai) ∈ (0,1]`) — no `exp(+λ·etai)` overflow,
  and **no pivoting needed** for the physical RT matrix (the benchmark's unpivoted elimination on the
  diagonally-dominant system gives 1e-16 residuals; the adding method is the standard stable layered-RT
  solver).
- Keep **LAPACK `dgbsv`** (band, `kl=ku=2`) and ED2's dense solver as **CPU-only test oracles**, never in
  the hot/device path.

---

## 9. The math, end to end (per band)

```
elai = clumping_leaf·lai ; ewai = clumping_wood·wai ; etai = elai + ewai
Optics (§6): ω, β, β0, k, emission  (leaf/wood blended)
Direct beam (beam bands):  down0(i) = down0(i+1)·exp(−k·etai(i))          top→bottom
Diffuse 2N+2 system: interior continuity with sources
    = beam-scattering particular solution (β0-driven; beam bands)
    + emission particular solution B_i=(1−ω)σT⁴ (emission bands)
  bottom BC = ground_refl·(down diffuse + transmitted beam) + ground_emission
  top BC    = incident diffuse
Solve (banded, §8) → up(i), down(i) at interfaces
Absorbed(i) = [down0(i+1)−down0(i)] + [down(i+1)−down(i) + up(i)−up(i+1)]   (flux divergence)
Split: leaf gets Absorbed(i)·leaf_frac(i) (§7)
Albedo_band = up(ncoh+1)/(incident_beam+incident_diff)
```
Band SW split of `rshort` into PAR/NIR beam/diffuse uses `fvis_beam=0.43`, `fvis_diff=0.52` (the ED2
constants), applied in the orchestration layer before the band loop. LW stays absolute; its incident is
`rlong`, its emission from `canopy_temp` and the ground.

---

## 10. Non-goals / deferred (documented, not forgotten)

- **Finite crown area (`cai<1`)** — ED2 forbids it with the two-stream; add only with a decision on the
  crown-only ED2 bugs B2/B3 (§12), using the LW-consistent forms.
- **Wood/branch-area allometry** — replaces the `wai = 0.1·lai` placeholder.
- **Soil/snow/water ground optics** — replaces the bare-soil `surface_optics` stub (B6-correct).
- **Directional reflectance / BRDF** — the `volscatt` view half; only needed for remote-sensing output.
- **Horizontal shading (`ihrzrad`, CCI), multiple-scattering (`icanrad=1`), big-leaf (`ibigleaf=1`)** —
  out of scope; the band-solver seam should not preclude a second scheme.

---

## 11. Test strategy (CTest targets, the EDTS analogue)

1. **Energy conservation** (per band, machine precision): `Σ absorbed + transmitted_to_ground +
   reflected_out_top + (emitted−absorbed for LW) = incident`. The strongest invariant.
2. **Analytic limits:** zero LAI → albedo = ground albedo; single cohort direct beam → Beer's law
   `exp(−k·L)`; `ω→0` → pure absorption; **leaf-angle limits** — horizontal LAD → `β = ρ/(ρ+τ)`,
   vertical/spherical → `β = 0.5` (the SCOPE win over CLM, and a direct unit test of §6c).
3. **Cross-check the optics** against the SCOPE reference values (`bf`, `k`, `sigb`, `β`) for canonical
   LIDFs (planophile/erectophile/spherical) — a pure-function test needing no solver.
4. **Regression vs ED2** — drive ED2's `sw_two_stream`/`lw_two_stream` on a fixed single-PFT stack with
   `cai=1` and **`orient_factor=0` (→ bf=1/3, μ̄=1)** so the CLM and SCOPE angle treatments coincide;
   match within tolerance. This isolates the multi-layer solver from the (intended) optics change.
5. **The B1 witness test:** a two-PFT stack with contrasting leaf reflectance. MEDS (correct per-cohort
   optics, §12) differs from unpatched ED2 (which reuses `pft(ncoh)` for all diffuse optics) — encode the
   physically-correct expectation.
6. **Determinism / kind-independence:** identical serial vs multicore (and CPU↔GPU once offloaded), per
   the existing MEDS 7/7 cross-backend discipline.

---

## 12. Bugs found in the ED2 two-stream scheme

Verified by adversarial multi-agent review against the LW sibling, `old_twostream_rad.f90`,
`multiple_scatter.f90`, and the CLM/Sellers/SCOPE references.

**B1 — Stale PFT index in the SW diffuse loop (real; affects production runs). ★**
`twostream_rad.f90:652-665`: `diffuseloop` computes each cohort's diffuse ω/β via `leaf_scatter(ipft)`
etc., but `ipft` is never assigned in the loop — it retains `pft(ncoh)` from `directloop` (`:556`). So
**every cohort uses the canopy-top cohort's PFT scattering** for all diffuse optics (the area weights are
correctly per-cohort; the direct-beam `beta0/epsil0` are fine). **Manifests** under `icanrad=2`,
shortwave, in any multi-PFT patch — corrupts diffuse PAR/NIR, sub-canopy light, and albedo. The LW
solver, `old_sw_two_stream`, and `sw_multiple_scatter` all index `pft(i)` correctly. **Fix:** `ipft =
pft(i)` at the loop top. **MEDS avoids it by construction** — per-cohort optics are *gathered into arrays*
(like `gather_pft_params`), so no reusable scalar index can go stale. *Worth reporting upstream.*

**B6 — `albedo_damp_nir` self-assignment (real; narrow).** `radiate_driver.f90:639`, bedrock branch:
`albedo_damp_nir = albedo_damp_nir` (uninitialised local); should be `albedo_soil_nir`. Bites only for
bedrock top-soil with snow/surface water. Fixed when `meds_surface_optics` is written (§5).

**B7 — TIR leaf/wood split omits `clumping_factor` (real; split-only).** `radiate_driver.f90:1103-1109`
(and big-leaf copy `:1215-1221`): the longwave leaf/wood split lacks the `clumping` factor that the SW
split and the solver's own `elai/etai` weighting use. Total LW conserved; only the leaf/wood partition
shifts. **MEDS avoids it** by the unified §7 split (same `elai`, all bands).

**B5 — Night-time fast-mean albedo diagnostic (real; diagnostic-only).**
`radiate_driver.f90:1653-1664` accumulates `fmean_albedo` with no daytime guard (the daily-mean sibling
has one), biasing the fast-mean with meaningless night-time soil albedo. In MEDS, gate albedo diagnostics
on daytime.

**B2 / B3 — crown-area-only, latent (moot at `cai=1`).** SW diffuse `μ` drops the `cai` factor that LW
`μ` and direct `μ0` carry (`:590` vs `:151`); and a contested `cai`-weighted LW blackbody source
(`:164`). Both require `cai<1`, which `icanrad=2` forbids → **zero impact in valid ED2 runs**. Relevant
only if MEDS ever adds finite crowns; use the LW-consistent forms then. (In the SCOPE formulation these
particular expressions are superseded anyway.)

*Correctly rejected as non-bugs:* `etai=0` division (resolvable cohorts have positive area); `rlong`
zero-guard (downwelling LW > 0 at the surface); unused `fvis_*_def` (documented fallback); backscatter/μ̄
singularities (only via out-of-range XML overrides). MEDS still validates inputs at construction to keep
these impossible.

---

## 13. Config & parameters (no hard-coded values, per MEDS rule)

**PFT file `[pft]`** (SoA in `pft_table_t`, derived optics via a pure `derive_rad_optics`):
`leaf_reflect_vis/nir`, `leaf_trans_vis/nir`, `leaf_emiss_tir`, `wood_reflect_vis/nir`,
`wood_trans_vis/nir`, `wood_emiss_tir`, `clumping_leaf`, `clumping_wood`, and the **Beta-LAD parameters**
(`mean_leaf_angle`, `leaf_angle_conc` → `(p,q)`; §6a).

**Main file `[canopy_radiation]`:** the band list (default VIS/NIR/LW with wavelength ranges +
`fvis_beam/fvis_diff`), `wai_lai_ratio` (0.1), `cosz_min`, `rshort_twilight_min`, `rt_method` selector,
and (until soil state) bare-soil ground optics (`soil_albedo_vis/nir`, `soil_emiss`, `soil_temp` source).
All required + presence-mapped like the rest of the config; derived optics computed by a pure routine so
overrides never leave them stale.

---

## 14. Phased plan

1. **Optics foundation** — `meds_rad_types`; `meds_leaf_angle` (LIDF + `volscatt` `chi_s`/`bf`/`G(μ)`);
   `meds_canopy_optics` (ω, β, β0, leaf/wood blend). Unit tests vs SCOPE reference values + the leaf-angle
   limits. Wire PFT/config fields.
2. **Band solver (leaf+wood, `cai=1`)** — `meds_twostream_band` + `meds_rad_linsys` (dense oracle first,
   then the banded/adding solver). Energy-conservation & analytic-limit tests; the `orient=0` regression
   vs ED2.
3. **Surface optics + orchestration** — `meds_surface_optics` (bare soil, B6-correct); `meds_canopy_
   radiation`: band split/regimes, per-patch walk over `site_t` (bottom→top index over each CSR slice,
   as ED2's reversed loop), leaf/wood split, diagnostic absorbed-PAR field. Regression + B1 witness test.
4. **Longwave band** — reuse the same solver with the emission source and `canopy_temp`.
5. **Wire into the stepper / leaf coupling** — feed per-cohort absorbed PAR into `meds_leaf_physiology`
   as the fast-loop light input; first real step of the reserved "couple leaf physiology into growth".
6. **Later** — wood-area allometry (replace `wai=0.1·lai`); soil/snow ground optics; finite crowns
   (B2/B3); directional reflectance; GPU offload of the banded solve.
