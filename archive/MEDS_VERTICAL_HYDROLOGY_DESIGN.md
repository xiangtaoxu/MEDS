# MEDS Vertical Hydrology — Module Design

A **stateless** vertical soil-water compute library for MEDS (`src/biophysics/meds_vertical_hydrology.f90`). It owns the one-dimensional ground-water column: **canopy rain interception**, **surface infiltration + ponding**, the **multi-layer soil-water balance** (inter-layer Darcy/Richards flux, drainage, runoff, root-uptake sink, soil evaporation). It is the below-ground continuation of the same 1-D grounded-Laplacian chain the plant-hydraulics kernel already solves — and it exists chiefly to **close the `hydro_env_t%(soil_psi, rhizo_cond)` boundary condition** that hydraulics today receives as a hand-set scalar.

The physics closure is **van Genuchten–Mualem** `ψ(θ), K(θ)` by default (the direction the ED2 community is moving; ED2 carries it as the `vG80` scheme), with **Campbell/Clapp–Hornberger a config option** for ED2-BC64 reproducibility. Because the hydraulics coupling is through **potential `ψ` — curve-independent at the interface (§3b) — not water content**, the retention curve is a free choice with no cross-module constraint. Both forms are closed-form invertible. The **numerics are CLM's, not ED2's**: a linearly-implicit backward-Euler tridiagonal (Thomas) solve on the flux-form mixed Richards equation, adaptive step-doubling, and machine-precision finite-volume conservation — the plant-hydraulics `HYDRO_SOLVER_BE` branch scaled from 3 nodes to `n_soil_layer_max`. This escapes the wet-soil stiffness that pins ED2's explicit RK4 at its floor step, and maps cleanly onto OpenMP-target (columns are the parallel axis). The MVP interior flux is the **plain gravitational (unit-gradient) form**; the Zeng–Decker (2009) equilibrium correction is **deferred to P2 with the aquifer BC**, where — done properly (retention-integral `ψ_E`, interior faces only) — it earns its coarse-layer-drift fix. As written with a linear-midpoint `ψ_E`, ZD is algebraically identical to plain gravity and buys nothing, so the MVP does not carry it (§3e).

It mirrors `src/plant/hydraulics/` and `src/biophysics/` (canopy RT) exactly: links `meds_shared` only, no `site_t` dependency, compiles and tests standalone against synthetic forcing, one public seam. The prognostic state — per-patch soil moisture and ponding, plus the **per-cohort canopy-interception film** — lives in the `meds_state` SoA and is passed by argument (the FATES `*Mem`/compute split). **Units are SI head in meters internally** (ED2 `slpots` convention), converted to **MPa** only at the plant-hydraulics seam via `grav_head`. Kinds `wp/ik`, `_wp` literals, `implicit none`, `pure`/`elemental` constitutive kernels, `error stop`, ≤132 cols, nvfortran-safe.

**Deferred siblings (out of scope here).** The canopy-air-space water/vapor balance and the leaf boundary-layer conductance kernel are *separate* biophysics modules, to be designed later; this document is the soil column only. Two forward seams are left for them: soil evaporation exposes its `T_ground`/`q_air`/`r_aero` forcing as the hook where the future canopy-air-space + energy balance plugs in (§3g), and per-cohort interception exposes a wetted-fraction `σ_w` as the hook the future canopy-air-space evaporation will consume (§3c).

---

## 1. Scope, target processes & reachability facts

### 1.1 Coordinate & sign conventions (fixed once, used throughout)

**ED2-style vertical coordinate: `z` = elevation [m], positive UP, `z = 0` at the ground surface** (MEDS's `slz`). Belowground the soil nodes are **negative** (`z_k < 0`, deeper = more negative, matching ED2 `slz`); aboveground canopy heights are positive — so a **single axis cleanly separates above- and below-ground distances**. Depth is `d = −z ≥ 0`. MEDS indexes soil layers **`k = 1` = top (shallowest, nearest surface) → `k = n_soil_layer_active` = bottom (most negative)**, the natural infiltration→drainage order (`z_node(1) ≈ −0.01 m`, `z_node(N) ≈ −soil_depth`). Layer thickness `dz_k > 0` and internode spacing `dz_node(k) = z_node(k) − z_node(k+1) > 0` are **positive magnitudes** (the upper node is less negative than the lower). Water-table elevation `z_wt ≤ 0`.

Matric potential `ψ ≤ 0` [m of head]; volumetric water content `θ` [m³ m⁻³]; **total hydraulic head `H = ψ + z`** (elevation head is the actual elevation `z`). Fluxes `q` are **positive downward** (the natural sense for infiltration/drainage): `q = −(−K·∂H/∂z) = K·∂ψ/∂z + K`, where the `+K` is gravity. The discrete forms in §3d/§3e use `dz_node, dz > 0` and `q` positive-down, so the **stencils are sign-agnostic to the elevation convention** — only `H = ψ + z`, the negative soil `z`, and the equilibrium potential `ψ_E = z_wt − z` (§3e) carry the ED2 sign.

### 1.2 Target processes and MVP vs deferred

| Process | Governing choice | MVP (P1) | Deferred |
|---|---|---|---|
| **Canopy interception** | **per-cohort** capacity-limited bucket, top→bottom throughfall cascade (§3c) | per-cohort `leaf_water` SoA, capacity `dewmx·PAI`, cascade sweep; exposes `σ_w` | Deardorff wetted-fraction *evaporation* + `(1−σ_w)` transpiration split (needs canopy-air-space + leaf energy balance) |
| **Throughfall + drip** | per-cohort cascade: rain_above − Δstore − evap → cohort below → soil (§3c) | ✔ | phase (snow) splitting |
| **Surface infiltration + ponding** | Neumann-flux top BC, Ksat-capped; excess → ponding `w_surface`, then Dunne runoff (§3d) | capped-infiltration + instantaneous surface runoff, small ponding store | Green–Ampt sharp-front; Dirichlet-head saturation switch |
| **Inter-layer soil water balance** | flux-form mixed Richards, implicit BE-Thomas; **plain gravity flux** (§3e, §5) | ✔ full interior solve (plain gravity) | Zeng–Decker equilibrium correction (retention-integral `ψ_E`), with the aquifer BC |
| **Bottom BC / drainage** | free drainage (unit-gradient) MVP; SIMTOP baseflow + unconfined aquifer + water table (§3e) | free drainage | aquifer/`z_wt`, TOPMODEL lateral (`lsm_hyd`) |
| **Surface + subsurface runoff** | Dunne saturation-excess `f_sat` + Horton infiltration-excess (§3e) | infiltration-excess only | `f_sat`/`f_max` CTI map, exponential baseflow |
| **Root-uptake sink** | **FATES delegation**: consume per-plant `hydro_flux_t%root_uptake`, ψ-limited `S_k` on the diagonal (§3f) | ✔ (prescribed-profile sink for standalone tests; hydraulics-driven when wired) | per-soil-layer fine-root nodes (hydraulics §16) |
| **Soil evaporation** | Swenson–Lawrence (2014) dry-surface-layer resistance + `α_soil` from `ψ_1` (§3g) | ✔ | McCumber–Pielke peat, litter resistance |
| **Retention/conductivity** | **van Genuchten–Mualem default** (closed-form invertible); Campbell/BC64 option (§3b) | ✔ | vertical texture variation |

Water-only in v1: **no soil energy/thermal, no freezing, no snow/frozen precip** — the advective/conductive heat fluxes ED2 co-integrates (`qw_flux_g`, `slcpd`, `fzcorr`) are out of scope until an energy balance exists. `fliq` is treated as 1 (all liquid); the frozen-soil throttle is a hook, not implemented. Relative to the CLM5/ED2 references, "vertical water balance" here means **liquid-water only** — snow, frozen precip, and the sfcwater phase stack are deferred (§14 P5).

### 1.3 Reachability facts from the ED2 reference

- ED2's **only** soil-top water flux in the RK4 derivative is ground evaporation (`w_flux_g(mzg+1)=wflxgc·wdnsi`). **Rain never infiltrates directly**: it buffers in a virtual/TSW pool and is dumped into the top layer up to saturation room *between* RK4 steps (`adjust_sfcw_properties`). MEDS replaces this with a proper conductivity-limited top Darcy flux (§3d) — an explicit improvement, not a port.
- ED2's alternate conductivity-based infiltration (`infiltration_method/=0`) is **dead code** — it calls `fatal_error(...)` immediately. Ignore.
- `icanrad=2` / bedrock class 13 / `ISOILBC` branches are all reachable ED2 configurations; the bedrock hydraulics-disable is a valid lower-BC mode we mirror as a config option, not a hard-coded `nsoil==13` switch.

---

## 2. Where it lives (library DAG)

### 2.1 The state/process wall

The DAG is `shared ← {allometry, plant} ← state ← demography ← aux ← main`, and `biophysics` is a **sibling of `plant`** that links `meds_shared` **only** (`CMakeLists.txt` `target_link_libraries(meds_biophysics PUBLIC meds_shared)`). So the kernel library **cannot name `site_t`** and **cannot `use meds_demography_types`**. This forces the identical split the plant-hydraulics kernel already uses: the compute is **stateless**, and the prognostic soil column is owned by `meds_state` (the one layer every science lib can see, which owns all prognostic SoA + the lockstep machinery).

**The controlling precedent: soil water is to `patch_index` what `psi` is to `cohort_block`** — prognostic state in `meds_state`, a stateless kernel in `meds_biophysics` that receives it as an `intent(inout)` value type, woven together only in `meds_aux`.

**DAG ownership of typed vs plain data (see §16 blockers).** `meds_config` is the DAG root in `src/shared` (it `use`s only `meds_kinds`/`meds_constants`/`meds_pft_params`/`meds_time`); it therefore **must not** `use meds_biophysics_types`. Consequently `soil_params_t` is a **biophysics** type and **never appears as a `meds_config_t` component**. `meds_config_t` carries only **plain `real`/`integer` allocatables** (`soil_layer_z`, `dz`, `dz_node`, `z_node`, and the per-layer texture arrays), derived in `derive_config` exactly like `height_edges`. The typed `soil_params_t` is assembled by a biophysics routine `meds_soil_parameters%derive_soil_params(cfg) → soil_params_t`, called from the aux/init layer that sees both `meds_config` and `meds_biophysics` — the exact analogue of `meds_optics%derive_rad_optics` filling `rad_pft_optics_t` (a biophysics routine, not `derive_config`). Likewise the **BC/solver enums are `integer(ik), parameter` codes in `meds_config` (shared)**, re-exported alongside `TS_*/BK_*/DIST_*`, so `meds_config_io` maps the TOML strings exactly like every existing `TS_MONTHLY`/`SM_LEUNING`/`GS_CARBON` selector; biophysics `use`s `meds_config` to read them.

```
shared ─┬─ allometry ─ state ─ demography ─┐
        │  (meds_config owns: soil geom +  │
        │   texture as PLAIN real arrays;  ├─ aux (fast-loop driver: RT → leaf → hydraulics
        │   BC/solver enum PARAMETERs)     │      → THIS kernel; the ⊥ weave, host-only;
        ├─ plant  (hydraulics: closes ─────┤      also calls derive_soil_params → soil_params_t)
        │          soil_psi BC)            │
        └─ biophysics(core: shared-only) ──┘
             meds_vertical_hydrology  ← Richards BE-Thomas kernel, array interface
             meds_soil_parameters   ← vG/Campbell curves + derive_soil_params (assembles soil_params_t)
             meds_biophysics_types    ← soil_column_t, vhydro_forcing_t, vhydro_flux_t, soil_params_t
```

### 2.2 Files & CMake wiring

| File | Role | Analogue |
|---|---|---|
| `src/biophysics/meds_biophysics_types.f90` (extend) | `soil_column_t`, `vhydro_forcing_t`, `vhydro_flux_t`, `soil_params_t`, `soil_opts_t`; `n_soil_layer_max` parameter | `meds_rad_types` / `meds_hydro_types` |
| `src/biophysics/meds_soil_parameters.f90` (new) | `pure`/`elemental` `psi_of_theta`, `theta_of_psi`, `hydr_cond`, `moist_capacity` (`C=dθ/dψ`) for **van Genuchten (default) + Campbell (option)**, closed-form inverses; **`derive_soil_params(cfg)`** assembling `soil_params_t` | `meds_soil_coms` constitutive fns + `derive_rad_optics` |
| `src/biophysics/meds_soil_solver.f90` (new) | hand-rolled tridiagonal `thomas_solve(a,b,c,d,x)` — **subroutine, `intent(out) x(n_soil_layer_max)`** | `meds_rad_linsys` |
| `src/biophysics/meds_vertical_hydrology.f90` (new) | **THE seam** `vertical_hydrology_flux`; inner device-eligible `soil_water_be_step`; interception/infiltration/evap sub-steps; the `psi_soil` aggregator | `meds_canopy_radiation` |
| `test/test_vertical_hydrology.f90` (new) | CTest target, links `meds_biophysics` + `meds_testsupport` | `test_plant_hydraulics` |

(No `bind(c)` C-API file is created now — the Python/C-API shim is deferred to the **end of biophysics development**, §14 P5.)

CMake: sources are GLOB'd with `CONFIGURE_DEPENDS` into `libmeds_biophysics`, so the new module files auto-add. Append a `test_vertical_hydrology` executable following the `test_canopy_radiation` block (`target_link_libraries(... PRIVATE meds_biophysics)`; `add_test(...)`). To wire the fast loop later, add `meds_biophysics` to `meds_aux`'s deps (`target_link_libraries(meds_aux PUBLIC …)` — it does not yet link biophysics). Build **nvfortran multicore** on the new module: a green ifx suite is not sufficient (§10).

### 2.3 The fast-loop / met-forcing gap and standalone testability

The master stepper `advance_one_step` calls only the slow loop `vegetation_dynamics`; there is **no sub-daily loop and no meteorological forcing** (transpiration/GPP are stubbed, `STUB_TISSUE_TEMP = 298.15`). So — exactly like RT, leaf, and hydraulics — this module **ships standalone on its unit test**: the seam takes **forcing passed in as a value type** (`vhydro_forcing_t`), never read from a global. When a met-forcing reader (`src/io`, precip + SW/LW + air T + VPD) and a fast-loop driver land, `advance_one_step` will drive it (§9); until then it stands alone. The still-absent piece is that met reader, not this kernel.

---

## 3. Vertical structure & governing equations

### 3a. Soil-layer discretization (`slz`/`dz`/`n_soil_layer_max`)

Fixed compile-time `integer(ik), parameter :: n_soil_layer_max = 20` (the array-shape ceiling; the active count `n_soil_layer_active ≤ n_soil_layer_max` from config, so the array shape never changes — GPU-friendly, no hot-path allocatables). The parameter is spelled out per the MEDS naming convention (unconventional ED2 acronyms like `nzg` are expanded; `dbh`/`nplant`/`pft`-style domain tokens are kept). **Two grid-construction modes, both from config (§8):** (i) an **exponential generator** — fine at surface, coarse at depth, to `soil_depth` (default 8.5 m, CLM-like) via the closed form below; or (ii) an **explicit user-supplied `soil_layer_z` array** of interface **elevations** (`z ≤ 0`, descending from `0`, e.g. `[0, −0.02, −0.06, …]`, ED2 `slz` style), used verbatim — so an ED2 coarse 12-layer grid, a CLM 20-layer grid, or any bespoke profile reproduce exactly. When `soil_layer_z` is present it **overrides** the generator (and sets `n_soil_layer_active = size(soil_layer_z)−1`); otherwise `soil_depth`+`grid_growth` build the exponential grid. Either way the same derived `soil_layer_z/z_node/dz/dz_node` arrays result, computed once in `derive_config` (the `height_edges` pattern) as **plain real arrays on `meds_config_t`**, then packed into `soil_params_t` by `derive_soil_params`:

```
soil_layer_z(k) = −soil_depth · (exp(f_grid·(k-1)/n_active) − 1)/(exp(f_grid) − 1)  k = 1..n_active+1  [m, elevation ≤ 0]
soil_layer_z(1) = 0                                                                   (surface)
dz(k)     = soil_layer_z(k) − soil_layer_z(k+1)                                             layer thickness > 0  [m]
z_node(k) = 0.5·(soil_layer_z(k) + soil_layer_z(k+1))                                       node (T-point) elevation ≤ 0 [m]
dz_node(k)= z_node(k) − z_node(k+1)                                            internode spacing > 0 (M-point) [m]
```

The gradient uses internode spacing `dz_node(k)`; the divergence uses layer thickness `dz(k)` — correct staggered finite volume (the same ED2 `dslzti`/`dslzi` distinction). The **surface-to-first-node** distance is `−z_node(1)` (the first-node **depth**, positive; used by the infiltration-capacity gradient, §3d), distinct from `dz_node(1)` (node1↔node2 spacing). `f_grid` is the geometric growth factor (config, default `≈ 4` giving a ~2–4 cm top layer). Below `n_soil_layer_active` (bedrock / depth-to-bedrock cutoff) `∂θ/∂t = 0`.

### 3b. Soil hydraulic constitutive curves (van Genuchten–Mualem default; Campbell option, `pure`/`elemental`)

**Decision: van Genuchten–Mualem is the MVP default; Campbell/BC64 is a config option.** The ED2 community is moving to van Genuchten as the default closure (ED2 already carries it as the `vG80` scheme); it fits measured retention data better than Clapp–Hornberger and is smooth across the air-entry region (no kink). The curve is a **free choice with no cross-module constraint**: the hydraulics kernel embeds **no soil curve** — it consumes only `env%(soil_psi [MPa], rhizo_cond [kg/s/MPa])` and couples through **potential `ψ`, curve-independent at the interface** (confirmed at `meds_plant_types.f90:119`, `meds_plant_hydraulics.f90:250,338`). Both vG and Campbell have **closed-form `θ(ψ)` and `ψ(θ)` inverses** (no Newton), so both are `pure`/`elemental` and cheap; vG's only extra care is its steeper near-saturation `K(θ)` and the specific capacity `C→0` at saturation.

Effective saturation `Se = clamp((θ−θr)/(θsat−θr), 0, 1)`; `m = 1 − 1/n`:

```
θ(ψ) = θr + (θsat−θr)·[1 + (α·|ψ|)^n]^(−m)      (ψ<0; θ=θsat for ψ≥0)   [m³/m³]   van Genuchten retention
ψ(θ) = −(1/α)·( Se^(−1/m) − 1 )^(1/n)                                   [m, ≤0]    closed-form inverse (no Newton)
K(θ) = max(K_min, K_sat · Se^l · [1 − (1 − Se^(1/m))^m]^2 )            [m/s]      Mualem (l = 0.5 pore-connectivity)
C(ψ) = dθ/dψ = α·m·n·(θsat−θr)·(α·|ψ|)^(n−1)·[1 + (α·|ψ|)^n]^(−m−1)    [1/m]      specific moisture capacity (≥0)
```

`α` [1/m] is the inverse air-entry scale, `n > 1` the pore-size-distribution index (`m = 1−1/n`), `l = 0.5` Mualem pore-connectivity. `K_min = 1.16e-13 m/s` (ED2 `hydcond_min`) floors the conductivity; `Se` is clamped so `ψ, K, C` stay finite — the **residual/air-dry floor**, required unconditionally for well-posedness (contrast the *wilting-point* uptake cutoff, §3f). A small `C`-floor near saturation keeps the BE diagonal non-singular where vG's `C→0`. The frozen-soil throttle `K *= exp(max(lnexp_min, −freezecoef·(1−fliq)))` is a **hook only** (`fliq≡1` in v1).

The **Campbell/BC64 option** (config `retention = "campbell"`, ED2-faithful, the §11-test-8 baseline) swaps in `ψ = ψ_sat·Se^(−b)`, `K = K_sat·Se^(2b+3)`, and the clamped closed-form `θ(ψ)` (`min(1,·)` so `θ ≤ θ_sat` at/above air-entry, no supersaturation). Interface conductivity by **upstream weighting** (§5) for the gravity/advective term (harmonic-mean smears wetting fronts); both curve families are `elemental`, branch-light, FPE-safe (clamp `Se`, floor `K`, guard `C`) under `-fpe0`/`-Ktrap=fp`. The round-trip + FPE tests (§11 test 2) assert both.

**Texture-class retention params** — vG: `θ_sat`, `θ_res`, `α`, `n`, `K_sat`; Campbell: `θ_sat`, `θ_res`, `ψ_sat`, `b`, `K_sat` — are REQUIRED from TOML (§7, §8), never hard-coded; the **threshold moistures** `θ_fc`, `θ_wp`, `θ_cp` are **derived** per class from the retention curve (§7). The one genuine cross-module contract is mass-bookkeeping (§9.2), independent of curve: the soil module debits `col%theta` with the same `hydro_flux_t%root_uptake` flux/units the hydraulics kernel reports. **Unit discipline:** always convert `ψ[m] = ψ[MPa]·1e6/(ρ_w·g)`.

**Texture-class retention params** (`θ_sat`, `θ_res`, `ψ_sat`, `b`, `K_sat`) are REQUIRED from TOML (§7, §8), never hard-coded; the **threshold moistures** `θ_fc`, `θ_wp`, `θ_cp` are **derived** per class from the retention curve (§7). **Unit discipline (ED2 quirk #3):** always convert `ψ[m] = ψ[MPa]·1e6/(ρ_w·g)` — never rely on ED2's `slpotcp`/`slpotwp` coincidence that holds only because `wdns=1000`.

### 3c. Canopy interception (per-cohort storage, capacity, throughfall cascade, drip)

**Decision: per-cohort interception with a prognostic `leaf_water` film on the cohort SoA, computed by a stateless per-layer kernel the orchestrator sweeps top→bottom.** This is ED2-faithful (ED2 carries per-cohort `leaf_water` with a wetted fraction `σ_w`, `rk4_derivs.f90:1463`) and it is the foundation the future canopy-air-space module needs: the `σ_w` this kernel exports is exactly the wetted fraction that will scale interception evaporation and the `(1−σ_w)` transpiration split once a canopy-air-space + leaf energy balance exist. Prognostic state: one per-cohort scalar `leaf_water` [kg m⁻²(ground)] on `cohort_block` (§6), carried through the cohort lockstep and summed on cohort fusion.

**The stateless kernel operates on one canopy layer** (one cohort, or the whole canopy lumped):

```
pai      = lai + sai                                          [m2/m2]  cohort plant-area index (leaf+stem)
f_pi     = alpha_pi · (1 − exp(−k_int·pai))                   [–]      interception fraction (Beer, per-layer; α_pi=1)
w_max    = dewmx · pai                                        [kg/m2]  holding capacity (dewmx = 0.1 kg/m2 per unit PAI)
q_grab   = f_pi · rain_above                                  [kg/m2/s]  potential interception from the flux above
```

Capacity-limited bucket update over `dt` — this is where MEDS **fixes ED2 quirk B-3** (ED2 lets `leaf_water` overshoot mid-step then caps it post-hoc; MEDS caps at grab time and drips the remainder as a proper flux):

```
room     = max(0, w_max − leaf_water)/dt                      [kg/m2/s]  free-capacity rate
q_intr   = min(q_grab, room + E_canopy)                       actually intercepted (bounded by capacity + evap headroom)
leaf_water = clamp(leaf_water + (q_intr − E_canopy)·dt, 0, w_max)
q_drip   = max(0, q_grab − q_intr)                            [kg/m2/s]  overflow drips to the layer below
throughfall_below = (rain_above − q_grab) + q_drip            [kg/m2/s]  gap-throughfall + drip → next cohort / soil
sigma_w  = min(1, (leaf_water/w_max)^(2/3))                   [–]        Deardorff wetted fraction (EXPORTED, §4)
```

**Vertical cascade (the orchestrator, host).** Cohorts are already height-sorted in the SoA; the orchestrator sweeps **top→bottom**, feeding each cohort's `throughfall_below` as the next cohort's `rain_above`, starting from `precip_rain` at the canopy top. The bottom cohort's `throughfall_below` is the **ground-reaching liquid `q_liq,g`** handed to infiltration (§3d). This is a *real* vertical cascade — an improvement over ED2's TAI-weighted split of one patch-total, which has no inter-cohort shading (quirk B-2, now superseded rather than reproduced). The degenerate call (`n_cohort = 0`, or a caller wanting a big-leaf closure) invokes the kernel **once** with the patch-total PAI and a single lumped `leaf_water` — "a single cohort mimicking the whole canopy."

**Evaporation is a forcing for now.** `E_canopy` [kg/m²/s] per cohort arrives on the forcing type (`vhydro_forcing_t`, capped at `leaf_water/dt + q_grab`). Until the canopy-air-space module exists there is **no humidity feedback**: the classic interception-loss *suppression* of transpiration (an evaporating film raises `can_shv`, shrinking the transpiration gradient) is **deferred** — it cannot be reproduced by a split PET partition and requires the shared prognostic `can_shv`, which is exactly why the canopy-air-space balance is its own future module. This kernel owns the *water* budget (capacity, throughfall, drip, storage) and **exports `σ_w`** as the seam that module will consume; the `(1−σ_w)` transpiration correction (ED2 bug B-1) also lands there, not here.

### 3d. Surface infiltration + ponding / surface water

The ground-reaching flux `q_liq,g` is partitioned into infiltration (the interior solve's **top face `q_½`**), ponding, and surface runoff — decided **before** the interior solve (operator split step 3, §5). Prognostic per-patch ponding store `w_surface` [kg m⁻²]. The infiltration capacity is the **top-face Darcy velocity from a ponded surface** (`ψ_0 = 0` at `z = 0`, gradient over the surface→first-node **depth** `−z_node(1) > 0`), so suction *increases* capacity as the soil dries:

```
q_pond_in = w_surface/dt + q_liq,g                            available surface water rate  [kg/m2/s]
q_inf,max = K(θ_1) · ( 1 + (0 − psi_1)/(−z_node(1)) )
          = K(θ_1) · ( 1 + |psi_1|/(−z_node(1)) )   (−z_node(1) = first-node depth)  ⇒ q_inf,max ≥ K(θ_1) > 0
q_infl    = min( q_pond_in, q_inf,max·rho_h2o )               MVP top-face MASS flux (Neumann)         [kg/m2/s]
w_surface = w_surface + (q_liq,g − q_infl/rho_h2o... )·dt      ponding accumulates the excess
q_runoff  = max(0, w_surface − w_pond,max)/dt                 ponding overflow → surface runoff        [kg/m2/s]
w_surface = min(w_surface, w_pond,max)
```

`q_inf,max` is a **Darcy velocity [m/s]** (surface-minus-node suction `(ψ_0 − ψ_1) = −ψ_1 ≥ 0` over the first-node depth, plus the `+1` gravity term, per §1.1's `q = K∂ψ/∂z + K`); multiplying by `rho_h2o` once gives the mass flux — there is no double-`ρ_w`. Using the **first-node depth `−z_node(1)`** (surface→first node, positive), **not** `dz_node(1)`, is the correct length scale.

**Why this differs from ED2 under heavy rain.** ED2 has no infiltration flux in its derivative at all: it buffers rain in a virtual/TSW pool and dumps it into the top layer up to *available pore space* `(θ_sat − θ_1)·dz` between RK4 steps (quirk A1, S1). Infiltration is therefore **limited by storage room, not by conductivity** — so under intense rain on a dry, low-K soil (e.g. clay) ED2 soaks in nearly everything until the top layer saturates and **under-generates Hortonian (infiltration-excess) runoff**. The MEDS cap `q_inf,max ≈ K(θ_1)·(1 + |ψ_1|/z_node(1))` instead routes rain beyond the soil's hydraulic capacity to ponding→runoff (correct infiltration-excess); and because the suction term `|ψ_1|` is large on a dry soil and collapses as it wets, it reproduces the **declining infiltration capacity** of a wetting front (Green–Ampt-like) that ED2 cannot represent. The single-node `q_inf,max` is itself a crude Green–Ampt proxy (it reads the current top-node `ψ_1`, not a sharp front); the sharp-front form is the P2 upgrade.

`w_pond,max` is a small config depth (default ~5 mm; ED2 `water_stab_thresh`). When the surface saturates (`θ_1 → θ_sat`) the P2 upgrade switches the top face from **Neumann flux to Dirichlet head** (`ψ_0 = w_surface/ρ_w`) — the standard infiltration/ponding switch. **Dunne saturation-excess runoff** (P2): `q_over = f_sat·q_liq,g` with `f_sat = f_max·exp(−0.5·f_over·z_wt)` (SIMTOP), split off before the infiltration-capacity cap. Snow/frozen precip and the ED2 sfcwater multi-layer stack are deferred (water-only v1).

### 3e. Inter-layer flux + bottom BC + drainage + runoff

Flux-form (finite-volume) mixed Richards, `q` positive down (elevation `z ≤ 0`, §1.1), `S_k` root sink [m³/m³/s]:

```
dz_k · dθ_k/dt = q_{k−½} − q_{k+½} − S_k·dz_k                                        (1)  [m/s]
q_{k+½} = −K_{k+½} · [ (psi_{k+1} − psi_k)/dz_node(k) − 1 ]   =  −K·∂psi/∂z + K       (2, MVP: plain gravity)
```

**MVP interior flux (P1): the plain gravitational form**, eq. (2). The `−1` (equivalently `+K` after the `−K` factor) is the gravity term of §1.1; `K_{k+½}` is **upstream-weighted** (§3b — upstream keeps wetting fronts monotone; harmonic-mean smears them; not ED2's log-linear interpolation). This converges to the unit-gradient drainage profile under the free-drainage bottom BC, and — with a no-flow bottom and zero sink — to a hydrostatic `q→0` state that exhibits the coarse-layer drainage drift ED2 also suffers (S9); removing that drift is exactly what the **P2** Zeng–Decker correction buys.

**Why Zeng–Decker is a P2 addition, not part of the MVP.** ZD (2009) references the gradient to the hydrostatic-equilibrium profile `ψ_E`:

```
q_{k+½} = −K_{k+½} · [ (psi_{k+1} − psi_k) − (psi_E,k+1 − psi_E,k) ] / dz_node(k)     (2′, Zeng–Decker, P2)
psi_E(z) = z_wt − z         (hydrostatic line; ∂psi_E/∂z = −1)                        [m]
```

Its value hinges **entirely on how `ψ_E,k` is averaged over a layer.** With the *linear-midpoint* mean `ψ_E,k = z_wt − z_node(k)`, `(ψ_E,k+1 − ψ_E,k) = z_node(k) − z_node(k+1) = dz_node(k)` and eq. (2′) collapses **algebraically to the plain form (2)** — toggling ZD then changes *nothing* and `z_wt` cancels from every interior face. The genuine coarse-layer-drift fix comes **only** from computing `ψ_E,k` as the layer-average of `ψ(θ(z_wt − z))` through the **nonlinear retention curve**, so that over a thick layer the mean `ψ_E` departs from the linear midpoint. **Real ZD (P2) therefore means: (i) `ψ_E` from the retention integral, (ii) applied on interior faces only, (iii) with the bottom-flux BC kept independent** (below). The MVP carries neither a `zeng_decker` flag nor the retention-integral machinery — it uses (2) and lets the free-drainage BC set the bottom.

**What `ψ_E` actually carries (and what it does not).** `ψ_E = z_wt − z` is the matric potential in hydrostatic equilibrium with a water table at elevation `z_wt`: the gravity part is `−z` (slope `−1`) and the water-table datum is the `+z_wt` offset. Operationally the two act very differently — the `−z` part is a *gradient* (`∂ψ_E/∂z = −1`), so it survives differentiation and **is** the gravity term, whereas `+z_wt` is a *constant offset* whose gradient is zero, so **`z_wt` cancels identically from every interior flux**. Above the table (`z > z_wt`, shallower) `ψ_E < 0` (suction), at it `ψ_E = 0`, below (`z < z_wt`, deeper) `ψ_E > 0` (positive pressure). So the water table couples in **not through the interior gradient** but only through (i) the equilibrium *target* the column relaxes toward (raising `z_wt` shifts the whole target wetter) and (ii) the **bottom boundary condition** (`SOIL_BC_AQUIFER` baseflow `∝ exp(−f_drai·z_wt)`). Under the MVP `SOIL_BC_FREE_DRAIN`/`SOIL_BC_BEDROCK` bottoms, `ψ_E` is never referenced absolutely, so `z_wt` is **inert**.

**Bottom boundary condition** (config `bottom_bc`, mirroring ED2 `ISOILBC`):

| mode | bottom face `q_{N+½}` | when |
|---|---|---|
| `SOIL_BC_FREE_DRAIN` (**MVP default**, ED2 ISOILBC=1) | unit-gradient `q = K(θ_N)` (pure gravity, positive-down) | idealized / benchmark |
| `SOIL_BC_AQUIFER` (P2, SIMTOP, Niu 2005/07) | recharge to lumped unconfined aquifer; `dW_a/dt = q_recharge − q_drai`, baseflow `q_drai = q_drai,max·exp(−f_drai·z_wt)` | physical water-table dynamics |
| `SOIL_BC_BEDROCK` (ED2 ISOILBC=0) | `q = 0` (no drainage) | shallow soil / bedrock cutoff; hydrostatic `q→0` test |
| `SOIL_BC_SLOPE` (ED2 ISOILBC=2) | slope-reduced `q = K(θ_N)·sin(sldrain)` | lateral drainage |

**Keep the bottom-flux BC independent of the interior correction.** As *continuous physics* the hydrostatic (`q→0`) and unit-gradient free-drainage (`q = K`) steady states are genuinely different attractors and cannot both hold strictly — that is the real content of "separate them." The clean way to honor it (P2): apply ZD's equilibrium reference to **interior** faces only, and let the **bottom face** be whatever the chosen BC prescribes (free-drainage `K(θ_N)`, aquifer recharge, or no-flow), so the two never algebraically collide. **There is no interior conflict to resolve in the MVP** — the MVP has no ZD — so the two idealized limits are simply validated under **distinct BCs** (§11 test 3): the **unit-gradient free-drainage steady state** (`dψ/dz = 0`, `q = K`) under `SOIL_BC_FREE_DRAIN`; the **hydrostatic `q→0` limit** under `SOIL_BC_BEDROCK` (no-flow) with `S = 0` and an initial `ψ = ψ_E` profile (this is also the case that exhibits coarse-layer drift under plain gravity and is *fixed* once real ZD lands at P2). Drainage and total runoff are diagnostics from the **same converged `ψ^{n+1}`**: `drainage = q_{N+½}(ψ^{n+1})·ρ_w` [kg/m²/s], `runoff_total = q_runoff + q_over + q_drai`. The **anti-overshoot limiter** (ED2 zeroes flux into a full / out of a dry layer) is *not* ported as a hard flux-kill (quirk A5 — it strands water and breaks conservation); the implicit solver + `θ`-clip with excess routed to ponding handles overshoot conservatively instead.

### 3f. Root water-uptake sink

**Decision: keep the FATES-HYDRO / ED2-hydraulics delegation seam, NOT CLM's empirical `β·E_demand` partition** (bundle D §6.1, §8). This module's job is only to (a) expose per-layer `ψ_soil,k` and per-layer soil→root conductance, and (b) **accept the plant-hydraulics kernel's per-plant uptake as the sink** `S_k`. This is the exact FATES "host owns water, vegetation owns uptake" contract, and it is what MEDS's already-planned stateless hydraulics kernel is built to consume/produce.

The patch's cohorts deposit their per-plant `hydro_flux_t%root_uptake` [kg/s] (the **plant** hydraulics type — distinct from this module's `vhydro_flux_t`, §4), scaled `×nplant` [plant/m²] and summed, then partitioned across layers by an exponential root-fraction profile `root_frac(k)` (per PFT, §7) weighted by layer conductance. Crucially, the demanded uptake is **made ψ-dependent by a smooth wilting cutoff** so the semi-implicit diagonal genuinely self-limits (the prescribed flux alone has `dS/dψ = 0` and would not):

```
U_total    = Σ_cohort nplant · flux%root_uptake                                  [kg/m2/s]  patch demand
w_k        = g_rhizo_k / Σ_j g_rhizo_j          (conductance weights, §9.1)        [–]
f_wilt(psi_k) = clamp( (psi_k − psi_wilt)/(psi_open − psi_wilt), 0, 1 )   (smoothstep-ramped) [–]
S_k        = U_total · w_k · f_wilt(psi_k) / (rho_h2o · dz_k)                     [m³/m³/s]  volumetric sink
dS_k/dpsi_k = U_total · w_k / (rho_h2o · dz_k) · f_wilt'(psi_k)                    (nonzero near wilting)
```

**On the two limiters — what is physics and what is scaffolding.** In the *coupled* path the sink is the hydraulics kernel's `root_uptake`, drawn across the rhizosphere edge as `q ∝ (ψ_soil − ψ_root)`; this **self-limits** as the layer dries (`ψ_soil → ψ_root ⇒ q → 0`) and the plant's own PLC/stomatal closure cuts `E` further. So neither `f_wilt` nor the `θ_wp` cap is *independent physics* — the drought limitation already lives in `Δψ = ψ_soil − ψ_root` and the plant vulnerability curve. They are kept for **numerical**, not physical, reasons:

- **`f_wilt` — a Jacobian/standalone term, not a second drought law.** A *prescribed* flux has `dS/dψ = 0`, contributing nothing to the implicit diagonal, so the semi-implicit step could over-extract within one solve. `f_wilt` supplies the stabilizing `dS_k/dψ_k = U_total·w_k/(ρ_w·dz_k)·f_wilt'(ψ_k) ≠ 0` (§5.2) and is the self-limiter for the **standalone prescribed-sink test**. **In coupled mode set `ψ_wilt` *below* the plant's PLC/stomatal cutoff** so `f_wilt ≈ 1` whenever the plant still transpires — otherwise it *double-counts* the plant's down-regulation and breaks the soil↔plant `ΔW` budget (§9.2). `ψ_wilt = −1.5 MPa` is the standalone default, not a coupled-mode physical threshold.
- **The hard `θ_wp` cap — a debug/standalone backstop, not operative physics.** `θ_wp = θ(ψ=−1.5 MPa)` is an agronomic convention, not a soil-water barrier (water still moves below it; `K > 0`, `ψ` finite until the *residual* floor). Its role is to stop the **standalone prescribed-sink** (blind to soil state) from over-drawing a layer to unphysical `θ`: post-solve extraction is capped at `max(0,(θ_k − θ_wp))·dz_k/dt` and the shortfall bookkept as `uptake_deficit` (§5.4). In coupled mode the hydraulic self-limit + `f_wilt` keep `θ` above `θ_wp`, so the cap should **rarely bind**; when it does, treat `uptake_deficit` as a **diagnostic of soil↔plant coupling inconsistency**, not as physics.

Distinct from §3b: the **residual/air-dry `θ_res` floor is the *curve* regularizer** (keeps `ψ, K, C` finite as `θ→0`) and is required **unconditionally**, uptake or not; `θ_wp` is a *plant* convention layered on top — do not conflate the two. ED2's single-`kroot`-deposit-then-redistribute-upward scheme (quirk A10) and its wilting-factor availability weighting are the reference for the profile shape but are superseded by the explicit hydraulic uptake.

### 3g. Soil evaporation

**Decision: adopt Swenson & Lawrence (2014) dry-surface-layer (DSL) resistance + Philip `α_soil`** (bundle D §5; directly targets the pan-tropical savanna/dry-forest soil-evaporation overestimation, better than ED2's cruder `ed_grndvap8`). Evaporation removes mass from the **top layer only** (top-face sink):

```
alpha_soil = exp( psi_1 · grav / (r_wv · T_g) )             soil-surface RH (Philip 1957)         [–]
L_dsl      = D_max · (θ_init − θ_1)/θ_init   for θ_1 < θ_init ;  L_dsl = 0 otherwise                [m]
r_soil     = L_dsl / (D_v · tau)                            DSL vapor resistance                  [s/m]
E_soil     = rho_atm · (alpha_soil·q_sat(T_g) − q_air) / (r_aw + r_soil)                            [kg/m2/s]
```

`θ_init = 0.8·θ_sat` (DSL initiation), `D_max ≈ 15 mm`, `D_v = 2.12e-5·(T_g/273.15)^1.75` (molecular vapor diffusivity), `τ = φ_air²·(φ_air/φ)^(3/b)` (air-filled-porosity tortuosity). This is the CLM5 form (Swenson & Lawrence 2014, GRACE/FLUXNET-validated). **ClimaLand** (CliMA, Julia) uses the same `α_soil` pore-space RH but the Sakaguchi–Zeng (2009) *exponential* DSL and full Monin–Obukhov coupling — a modern cross-check, and the reason the *linear* CLM5 DSL is preferred here (one knob `D_max` vs three fitted `d_ds, α, p`). It replaces ED2's `β = 0.5(1−cos(π·smterm))` Lee–Pielke throttle, which acts on `q_sat` (no pore-RH), has no diffusion length, and over-evaporates at intermediate moisture — the exact bias Swenson–Lawrence fixed; do **not** stack a `β` on top of `α_soil` (double-counts moisture limitation).

**Prescribed `T_ground` — the seam for a future soil-energy module.** `α_soil` and `q_sat(T_g)` need the ground skin temperature `T_g`, properly the solution of a surface *energy* balance MEDS does not yet have. So for v1 `T_ground` (with `q_air`, `r_aero`) is a **forcing** on `vhydro_forcing_t`; the cleanest default is `T_ground = T_air` (or a damped/lagged air temperature). Consequence: MEDS soil evaporation is **supply/diffusion-limited only, not energy-limited** — no net-radiation cap on `LE`, and **no evaporative-cooling feedback** on `T_g` (which would depress the skin temperature and self-limit `E`). `E_soil` is capped by top-layer available water. `forcing%t_ground` is the explicit drop-in seam: when a soil-energy balance lands (the natural sequel that solves `T_g` and unlocks `LE`-cooling and frozen-soil handling), only the *source* of `T_g` changes — the flux kernel is untouched. ED2's `ggsoil0`/`kksoil` MH91 forms remain the fallback for a future `IED_GRNDVAP`-style selector.

---

## 4. Public seam & data types (the RT analogue)

Declared in `meds_biophysics_types` (its docstring already reserves it as "the intended home for future energy-balance and hydrology types"). The soil types are **renamed away from the plant `hydro_flux_t`/`hydro_forcing_t`** (which `meds_plant_types` already defines for per-plant hydraulics, and which §9.2 consumes) so the aux fast loop can `use` both without aliasing: `vhydro_forcing_t` / `vhydro_flux_t` (per-ground-area soil fluxes) vs the plant `hydro_flux_t` (per-plant, `%root_uptake`). All **pure DATA** value types — the kernel never sees `site_t`. **Each derived-type component carries its own default initializer** (Fortran initializes only the last name on a shared-initializer line, so every real is written out) to stay clean under `-fpe0`/`-Ktrap=fp`.

```fortran
integer(ik), parameter :: n_soil_layer_max = 20_ik              ! compile-time max column depth

type :: soil_column_t                                          ! the MUTABLE prognostic slice (one patch)
   real(wp) :: theta(n_soil_layer_max) = 0.0_wp                 ! [m3/m3] volumetric soil moisture (PROGNOSTIC)
   real(wp) :: w_surface = 0.0_wp                               ! [kg/m2] ponded surface water
   real(wp) :: w_aquifer = 0.0_wp                               ! [kg/m2] lumped aquifer store (SOIL_BC_AQUIFER, P2)
   real(wp) :: z_wt      = 0.0_wp                               ! [m] diagnosed water-table depth (P2; feeds psi_E)
end type
! Canopy interception is per-COHORT: the leaf_water film lives on cohort_block (§6), not here.

type :: vhydro_forcing_t                                       ! soil-column boundary conditions (read-only)
   real(wp) :: precip_ground = 0.0_wp                          ! [kg/m2/s] GROUND-reaching liquid q_liq,g (post interception cascade, §3c)
   real(wp) :: root_uptake(n_soil_layer_max) = 0.0_wp          ! [kg/m2/s] per-layer transpiration DEMAND (×nplant summed)
   real(wp) :: t_ground = 298.15_wp                            ! [K] ground skin temp — FORCED (=T_air) until soil energy (§3g)
   real(wp) :: q_air = 0.0_wp, rho_air = 0.0_wp                 ! [kg/kg],[kg/m3] canopy-air humidity & density (soil evap)
   real(wp) :: r_aero = 0.0_wp                                  ! [s/m] aerodynamic resistance (soil evap series)
end type

type :: soil_params_t                                         ! per-column geometry + texture (biophysics type; per site)
   integer(ik) :: nzg_active = n_soil_layer_max               ! active layer count (n_soil_layer_active)
   real(wp) :: soil_layer_z(n_soil_layer_max+1) = 0.0_wp
   real(wp) :: z_node(n_soil_layer_max)   = 0.0_wp
   real(wp) :: dz(n_soil_layer_max)       = 0.0_wp
   real(wp) :: dz_node(n_soil_layer_max)  = 0.0_wp             ! [m] geometry (derived once from cfg plain arrays)
   real(wp) :: theta_sat(n_soil_layer_max) = 0.0_wp
   real(wp) :: theta_res(n_soil_layer_max) = 0.0_wp           ! [m3/m3] porosity + residual
   real(wp) :: ksat(n_soil_layer_max)      = 0.0_wp           ! [m/s]   saturated conductivity
   real(wp) :: vg_alpha(n_soil_layer_max)  = 0.0_wp           ! [1/m]   van Genuchten inverse air-entry (DEFAULT)
   real(wp) :: vg_n(n_soil_layer_max)      = 0.0_wp           ! [–]     van Genuchten pore-size index (>1)
   real(wp) :: psi_sat(n_soil_layer_max)   = 0.0_wp           ! [m]     Campbell air-entry potential (option)
   real(wp) :: b_camp(n_soil_layer_max)    = 0.0_wp           ! [–]     Campbell/Clapp–Hornberger exponent (option)
   real(wp) :: theta_fc(n_soil_layer_max)  = 0.0_wp
   real(wp) :: theta_wp(n_soil_layer_max)  = 0.0_wp
   real(wp) :: theta_cp(n_soil_layer_max)  = 0.0_wp           ! [m3/m3] DERIVED thresholds (from curve, §7)
   real(wp) :: root_frac(n_soil_layer_max) = 0.0_wp           ! [–] normalized root fraction (Σ=1)
end type

type :: soil_opts_t                                           ! pre-extracted selectors + tolerances (NOT whole cfg)
   integer(ik) :: solver      = SOIL_SOLVER_BE                ! codes live in meds_config (shared)
   integer(ik) :: retention   = SOIL_RETENTION_VG            ! SOIL_RETENTION_VG (default) | SOIL_RETENTION_CAMPBELL
   integer(ik) :: bottom_bc   = SOIL_BC_FREE_DRAIN           ! SOIL_BC_FREE_DRAIN|AQUIFER|BEDROCK|SLOPE
   integer(ik) :: linearize   = SOIL_LIN_FROZEN              ! SOIL_LIN_PICARD | SOIL_LIN_FROZEN
   integer(ik) :: substep     = SOIL_SUBSTEP_ADAPTIVE        ! | SOIL_SUBSTEP_FIXED
   logical     :: zeng_decker = .false.                      ! P2 only; MVP uses plain-gravity flux (§3e)
   real(wp)    :: rtol = 1.0e-3_wp, atol = 1.0e-4_wp
   real(wp)    :: h_init = 900.0_wp
   integer(ik) :: max_substep = 200_ik, max_picard = 5_ik
   real(wp)    :: w_pond_max = 5.0_wp, dewmx = 0.1_wp, intercept_alpha = 1.0_wp, intercept_k = 0.5_wp
   real(wp)    :: dsl_dmax = 0.015_wp, dsl_theta_init = 0.8_wp
   real(wp)    :: f_drai = 2.5_wp, q_drai_max = 5.5e-6_wp, f_over = 0.5_wp, f_max = 0.4_wp
   logical     :: debug_error = .false.                      ! error stop on cap-hit / mass_resid in Debug
end type

type :: vhydro_flux_t                                        ! outputs + diagnostics (the "band_out_t")
   real(wp) :: infiltration = 0.0_wp, drainage = 0.0_wp
   real(wp) :: runoff_surf  = 0.0_wp, runoff_sub = 0.0_wp     ! [kg/m2/s] boundary fluxes
   real(wp) :: soil_evap = 0.0_wp, canopy_evap = 0.0_wp
   real(wp) :: throughfall = 0.0_wp, drip = 0.0_wp            ! [kg/m2/s]
   real(wp) :: uptake_total = 0.0_wp                          ! [kg/m2/s] Σ S_k·dz_k (budget term)
   real(wp) :: uptake_deficit = 0.0_wp                        ! [kg/m2/s] capped (unmet) sink (budget term)
   real(wp) :: clip_excess = 0.0_wp                           ! [kg/m2/s] θ-clip water routed to ponding (budget term)
   real(wp) :: psi_soil(n_soil_layer_max) = 0.0_wp            ! [MPa] per-layer matric potential (EXPORTED to hydraulics)
   real(wp) :: soil_cond(n_soil_layer_max) = 0.0_wp           ! [kg/s/MPa per plant] soil→root conductance (§9.1)
   real(wp) :: mass_resid = 0.0_wp                            ! [kg/m2] closed-budget residual (should be ~0)
   integer(ik) :: nsub = 0_ik                                 ! adaptive sub-steps taken
   logical  :: converged = .true.                             ! .false. on any cap-hit (§5.3)
end type
```

**The public seam** (host, derived-type interface). It takes a **small `soil_opts_t` selectors record**, not the whole `meds_config_t` — matching the RT precedent (`canopy_radiation(opt, forcing, …)` deliberately takes no `cfg`), so the biophysics seam does not couple to the entire run-config type. `soil_opts_t` is filled once from `cfg` at the aux layer.

```fortran
subroutine vertical_hydrology_flux(col, forcing, params, opts, dt, flux)
   type(soil_column_t),    intent(inout) :: col        ! theta(:)+stores — SoA slice, owned OUTSIDE, only mutable thing
   type(vhydro_forcing_t), intent(in)    :: forcing     ! BCs (read-only)
   type(soil_params_t),    intent(in)    :: params      ! geometry + texture (read-only)
   type(soil_opts_t),      intent(in)    :: opts        ! pre-extracted selectors + tolerances (no whole cfg)
   real(wp),               intent(in)    :: dt          ! [s] fast step
   type(vhydro_flux_t),    intent(out)   :: flux        ! boundary fluxes + psi_soil + diagnostics
end subroutine
```

**Second entry point: the per-cohort interception kernel** (§3c). Stateless per canopy layer; the orchestrator loops the height-sorted cohorts **top→bottom**, feeding each `throughfall` as the next cohort's `rain_above` (the top cohort's `rain_above` is canopy-top precip; the bottom cohort's `throughfall` is `forcing%precip_ground`). Not `elemental` — the cascade is a sequential recurrence, one call per cohort:

```fortran
pure subroutine intercept_canopy_layer(leaf_water, rain_above, lai, sai, e_canopy, dt,        &
                                       dewmx, k_int, alpha_pi, throughfall, drip, sigma_w)
   real(wp), intent(inout) :: leaf_water                ! [kg/m2] per-cohort film (PROGNOSTIC — cohort_block, §6)
   real(wp), intent(in)    :: rain_above, lai, sai      ! [kg/m2/s], [m2/m2]  flux from above + cohort area
   real(wp), intent(in)    :: e_canopy, dt              ! [kg/m2/s] wet-film evaporation demand (forcing), [s]
   real(wp), intent(in)    :: dewmx, k_int, alpha_pi    ! capacity + Beer-fraction params (from soil_opts_t)
   real(wp), intent(out)   :: throughfall, drip, sigma_w ! [kg/m2/s],[kg/m2/s],[–] to layer below + wetted fraction
end subroutine
```

Two threads advancing two patches touch disjoint `col%theta(:)` ⇒ pure/reentrant. No `save`, no module variables, no hot-path allocatables. The interior arithmetic is delegated to a **bare-array, device-eligible** inner routine (§5, §10) — the `growth_step` precedent — and `soil_opts_t`/derived types stay strictly out of that inner routine.

---

## 5. The solver / numerics

**Decision: linearly-implicit backward-Euler (Rosenbrock-1) on the flux-form mixed Richards equation, hand-rolled Thomas sweep over the fixed `n_soil_layer_max` column, wrapped in the hydraulics stateless + adaptive-step-doubling idiom** (bundles D, E). Explicitly **reject ED2's explicit adaptive RK4** — it is exactly the wet-soil stiffness that pins ED2 at its floor step (K(θ) spans 6–10 orders of magnitude; thin top layers give a diffusion CFL `Δt ≤ dz²/(2·max D)` collapsing to sub-seconds when wet). Backward Euler is **A-stable and L-stable** ⇒ one long `dt` regardless of stiffness. Drop the hydraulics `EXPM` branch — a closed-form matrix exponential is the wrong tool at `n_soil_layer_max` dimension; the soil column is the `HYDRO_SOLVER_BE` branch only, scaled from 3 nodes to `n_soil_layer_max` (Thomas instead of 3×3 Cramer).

### 5.1 Operator-splitting order

Sequential boundary/source split (Lie–Trotter, 1st-order — matches the 1st-order BE interior, nothing wasted; each is a cheap non-stiff bucket) then one coupled implicit interior solve:

1. **Interception cascade** (§3c): sweep the height-sorted cohorts top→bottom, each updating its per-cohort `leaf_water` film (capacity `dewmx·PAI`, minus `E_canopy`) and passing `throughfall` down; the bottom cohort's `throughfall` is the ground-reaching `forcing%precip_ground`.
2. **Infiltration / ponding partition** (§3d): cap `q_liq,g` at top-layer capacity; excess → `w_surface` then surface runoff. Sets the interior **top face** (Neumann `q_½ = infiltration`; P2 Neumann→Dirichlet switch on surface saturation).
3. **Soil evaporation** (§3g): DSL flux as an additional top-layer sink.
4. **Coupled implicit interior solve** — inter-layer redistribution + drainage + root sink, **ONE tridiagonal Thomas solve**. Do **not** split redistribution from drainage; do **not** split the root sink out — fold the ψ-limited `S_k` on the diagonal (semi-implicit, §3f) so extraction backs off as a layer dries.
5. **Post-solve clip + book-keeping**: clamp `θ ∈ [θ_res, θ_sat]`, routing any **saturation excess upward to ponding and counting it** as `flux%clip_excess`; **cap per-layer extraction at `max(0,(θ_k−θ_wp))·dz_k/dt`** and report the capped shortfall as `flux%uptake_deficit`. Both clip terms enter the mass budget (§5.4).

### 5.2 Linearization (Celia modified-Picard, production) / frozen-coefficient (MVP)

Discretize (1)+(2) in flux form; evaluate `q^{n+1}` by Taylor/Newton in the layer increments. The **conservative Celia et al. (1990) modified-Picard** iterates on `ψ`, updating `θ` through `θ(ψ)` so mass is conserved (a pure ψ-based Richards leaks mass):

```
(C^m/dt)·(psi^{m+1} − psi^m)  =  ∇·[K^m ∇(psi^{m+1} − psi_E)]  −  (θ^m − θ^n)/dt  −  S(psi^{m+1})
```

`C^m = dθ/dψ` and face-`K^m` frozen at iterate `m` ⇒ the matrix stays linear and tridiagonal, one Thomas solve per iterate (cap `max_picard ~ 3–5`). In the **MVP `ψ_E ≡ 0`** (the plain-gravity flux of §3e eq. 2; the retention-integral `ψ_E` enters only at P2). **MVP fallback: drop the Picard loop** — one frozen-coefficient (`K,C` at `θ^n`) linearly-implicit BE step per substep (Rosenbrock-1, identical to hydraulics §5.2's single linear solve), leaning on adaptive step-doubling for the nonlinearity.

The tridiagonal rows (interior `k`), with `K` and `ψ_E` frozen:

```
a_k · Δψ_{k−1} + d_k · Δψ_k + c_k · Δψ_{k+1} = r_k
a_k = −K_{k−½}/dz_node(k−1)
c_k = −K_{k+½}/dz_node(k)
d_k = C_k·dz_k/dt − a_k − c_k + dS_k/dψ_k·dz_k        (ψ-limited root sink on the diagonal, §3f)
r_k = q_{k−½} − q_{k+½} − S_k·dz_k − (θ^m−θ^n)·dz_k/dt
```

Top row uses `q_½ = infiltration − E_soil` (Neumann); bottom row uses the §3e bottom-BC flux.

**The `θ` update is the conservative flux-divergence update — stated authoritatively.** After the solve, `θ^{n+1}` is defined by the finite-volume balance using boundary/interface fluxes diagnosed from the converged `ψ^{n+1}`:

```
θ_k^{n+1} = θ_k^n + (dt/dz_k)·( q_{k−½}(ψ^{n+1}) − q_{k+½}(ψ^{n+1}) − S_k(ψ^{n+1})·dz_k )
```

In the **production Picard path** this `θ^{n+1}` and `θ(ψ^{n+1})` reconcile to the iteration tolerance (Celia), so the stored `θ` lies on the retention curve. In the **frozen/non-Picard MVP path** the stored `θ` is off-curve by the (bounded) linearization residual — step-doubling (§5.3) controls this residual. In **both** cases the **exported `ψ_soil` is re-derived from the stored `θ` via `psi_of_theta`**, never from the linear-solve `ψ` iterate, so `ψ_soil` is always consistent with the conserved `θ` handed to the hydraulics kernel.

### 5.3 Stateless contract, substep / step-doubling, cap-hit contract

Fixed `n_soil_layer_max`; loop `1:nzg_active`; all scratch as fixed-length **automatic** arrays (`a,b,c,d,cp,dp(n_soil_layer_max)`), stack-friendly on device. Reuse the hydraulics step-doubling controller verbatim: freeze coefficients at substep start, take one step `h` and two `h/2`,

```
err = max_k |θ_{h/2,k} − θ_{h,k}| / (atol + rtol·|θ_k|)
accept extrapolated θ_{h/2} if err ≤ 1 ;  h ← h·clip(safety·err^(−1/(p+1)), fmin, fmax),  p = 1
```

BE is unconditionally stable, so substepping buys **accuracy at the wetting front and Picard convergence**, not stability: quiescent/well-drained ⇒ ~1 substep; active infiltration front ⇒ auto-refine to `O(10¹)`. A `SOIL_SUBSTEP_FIXED` path runs a warp-uniform fixed count for the GPU offload (data-dependent `nsub` is warp divergence).

**Cap-hit contract (defined, not silent).** Two caps can be hit:

- **`max_substep` (adaptive step-doubling):** if the wetting front still needs `err ≤ 1` after `nsub = max_substep`, the last extrapolated step is **accepted** with `flux%converged = .false.` and a logged diagnostic.
- **`max_picard`:** if the Celia iteration has not met `rtol/atol` after `max_picard` iterates, the last iterate is **accepted** with `flux%converged = .false.` and a logged diagnostic.

In **Debug** (`opts%debug_error = .true.`) either cap-hit is escalated to `error stop` — the same discipline as the `mass_resid` check (§5.4). In production the run continues with the flag and diagnostic set, so an under-resolved step is never silently indistinguishable from a converged one. A validation case (§11 test 4b) deliberately drives a sharp front into the cap to exercise this path.

### 5.4 Mass-conservation budget check

Flux form guarantees conservation **structurally**: summing the conservative update (§5.2) over `k`, every interior face `q_{k+½}` telescopes to zero, leaving only the two boundary faces and the total sink. The **top face carries both infiltration and soil evaporation** (`q_½ = Infiltration − E_soil`, §5.2):

```
Σ_k (θ_k^{n+1} − θ_k^n)·dz_k  =  dt·( (Infiltration − E_soil) − Drainage − Σ_k S_k·dz_k )        (layer store)
```

**Diagnose the boundary fluxes from the same converged `ψ^{n+1}` used in the update** (`Infiltration = q_½(ψ^{n+1})`, `E_soil` from step 3, `Drainage = q_{N+½}(ψ^{n+1})`, `Uptake_k = S_k(ψ^{n+1})`) — never recompute from an independent formula — so this layer-store balance closes to machine precision for the conservative flux-divergence `θ` update of §5.2 (the direct analog of hydraulics §6.3's `ΔW_stored = (uptake−E)·dt`). The **total site water balance across all stores** closes to machine precision **only after the two clip terms of §5.1 step 5 are added as explicit budget terms** — the `θ`-clip excess (an internal transfer into `w_surface`, hence conserved) and the capped unmet-sink deficit (reported, not silently dropped):

```
Δ(Σ_cohort leaf_water + w_surface + Σ_k θ_k·dz_k·ρ_w + w_aquifer)
   = dt·( precip − canopy_evap − soil_evap − uptake_total − runoff_total − drainage )
   with  uptake_total = Σ_k S_k·dz_k·ρ_w  (actual, after the θ_wp cap)  and  the capped deficit
         uptake_deficit  and clip_excess  carried as explicit diagnostics so the residual is exactly the
         floating-point round-off, not a hidden clip.
```

Ship `flux%mass_resid` as a diagnostic and `error stop` in Debug if `|resid| > atol` — the demography-invariant discipline.

### 5.5 nvfortran-safe Thomas solve

Hand-rolled as a **subroutine with `intent(out) x(n_soil_layer_max)`**, never an array-valued function fed into a call (issue #7 — silently wrong at `-O2`, segfault at `-O0`; a green ifx run hides it). Never `call apply(thomas(a,b,c,d))`; bind first: `call thomas_solve(a,b,c,d,x); … use x`. Forward-elimination scratch (`cp,dp`) are local automatics, never returned into a call. Diagonally dominant (`d_k` dominates via `C·dz/dt`) ⇒ **no pivoting needed**.

---

## 6. State additions (per-patch soil column + per-cohort interception film)

Two prognostic homes. Per bundle C, **geometry + texture stay shared/immutable per site** (same discretization, same `soil_params_t`), and the **prognostic water lives where its process lives**: soil moisture + ponding per PATCH (each patch has a distinct transpiration sink → distinct θ, ED2-faithful — reconciling "site shares the soil column" with per-patch state), the interception film per COHORT (§3c).

**Per-patch fields** on `patch_index` in `meds_demography_types.f90`:

```fortran
real(wp), allocatable :: soil_water(:,:)      ! (n_soil_layer_max, patch) [m3/m3]  PROGNOSTIC — theta(:) SoA
real(wp), allocatable :: surface_water(:)     ! (patch)                   [kg/m2]  ponding store
```

**Per-cohort field** on `cohort_block`:

```fortran
real(wp), allocatable :: leaf_water(:)        ! (cohort)                  [kg/m2]  intercepted canopy film (§3c)
```

**Patch lockstep** — patches have no single reorder routine; carry the two patch fields through **six** sites:

1. `patch_alloc` — allocate; init `soil_water = θ_fc`, `surface_water = 0` **for a true cold-start column only**.
2. `patch_ensure_capacity` — `tmp%…(:,1:m)=…` grow-copy.
3. `sort_patches` (`meds_demography_fusefiss`) — permute with the patch order.
4. `patch_compact` (`meds_demography_fusefiss`) — pack live patches.
5. **`fuse_2_patches` (`meds_demography_fusefiss:435`)** — the **area-weighted blend** site: `soil_water`, `surface_water` combined by area-fraction average exactly where `recruit_pool` is already area-weight-combined (line 447). Omitting this permutes/packs the columns but never blends them on fusion — the intensive-quantity fusion invariant.
6. `site_free` — patch `deallocate`.

**Cohort lockstep** — `leaf_water` rides the **single centralized cohort reorder** (`cohort_reorder`/`copy_cohort_slot`/`rebuild_csr`/`cohort_ensure_capacity` — the CLAUDE.md "add a per-cohort field → update *these*" rule) plus the fuse/split/recruit/terminate sites. Conservation (the conserved quantity is total canopy water `Σ_cohort leaf_water`, §5.4): **fusion** sums the two films (kg/m² of ground is extensive per cohort); **fission** splits by daughter area share; **recruitment** inits `leaf_water = 0`; **termination** sheds the film to `w_surface` (conserved, a small explicit budget term) rather than dropping it. Unlike AGB, `leaf_water` is *not* geometry-derived, so it is **not** in any `set_cohort_size`.

**Creation-site conservation (disturbance-gap init).** A treefall/disturbance **gap fragment is carved from an existing donor patch whose soil is physically unchanged by canopy loss**, so the gap column **copies the donor's `soil_water(:,donor)` and `surface_water`** — it must **not** be reset to `θ_fc` (resetting would inject/remove `(θ_fc − θ_donor)·dz·area_gap ≠ 0` of water and break the §5.4 site budget). Only **`init_bare_ground` / true cold-start columns** initialize to `θ_fc`, stores `= 0`. So: `apply_patch_disturbance` gap fragment → **copy donor**; `init_bare_ground` → `θ_fc`.

The orchestration layer (§9) copies one patch's SoA slice → `soil_column_t` value → kernel (`intent(inout)`) → writes back into `patch%soil_water(:,ip)`; the interception sweep updates `cohort%leaf_water(:)` in place.

The `soil_params_t` geometry/texture is assembled once per site by the biophysics routine `derive_soil_params(cfg)` (§2.1) from the plain geometry/texture arrays that `derive_config` placed on `meds_config_t`, and cached at the site/aux level — **never a `meds_config_t` component and never per patch**.

---

## 7. Soil & PFT/root parameters

**No hard-coded model parameters** (CLAUDE.md): every texture/retention/root parameter is REQUIRED from TOML, presence-mapped, `error stop` if missing. Only genuine universal constants go in `meds_constants` (§8.2). The retention params — **van Genuchten** `θ_sat`, `θ_res`, `α`, `n`, `K_sat` (default) or **Campbell** `θ_sat`, `θ_res`, `ψ_sat`, `b`, `K_sat` (option) — are supplied **directly per class from TOML**; MEDS runs **no in-model pedotransfer** on the hot path. Offline default-generators seed the `[soil]` table below: **Carsel & Parrish (1988) / Wösten (1999)** texture-class vG parameters for the default, and the ED2.2 **Cosby-1984 pedotransfer** for the Campbell option — neither embedded as source literals. If a texture→param pedotransfer is ever wanted at run time, its coefficients must be loaded at run start via a `set_soil_pedotransfer([soil_pedotransfer])` installer into `protected` module variables — the `meds_allometry`/`set_allometry` precedent — never as in-source literals:

```
(reference / offline only — Cosby 1984, seeds the Campbell-OPTION defaults, NOT in-model literals)
ksat     = 10^(−0.60 + 1.26·sand − 0.64·clay) · 0.0254/3600      [m/s]
b_camp   = 3.10 + 15.7·clay − 0.3·sand                            [–]
psi_sat  = −(10^(2.17 − 0.63·clay − 1.58·sand)) · 0.01            [m]
theta_sat= 0.01·(50.5 − 14.2·sand − 3.7·clay)                     [m3/m3]
```

**Threshold moistures are DERIVED, not required.** `θ_fc`, `θ_wp`, `θ_cp` are computed per class in `derive_config` from the supplied retention params at fixed potentials — they are removed from the required-key list (resolving the §8.1 required-vs-derived ambiguity):

```
theta_wp = theta_of_psi(psi = −1.5 MPa)      theta_cp = theta_of_psi(psi = −3.1 MPa)
theta_fc = theta_of_psi(psi = psi_fc)        psi_fc a single documented field-capacity head (config default, S8)
```

Reference **van Genuchten** default table (Carsel & Parrish 1988 texture-class parameters; loaded verbatim from TOML). Representative rows for defaults/tests (`θ_fc`/`θ_wp` are *derived* from the curve, illustrative — they depend on the `psi_fc` definition, §8):

| # | key | θ_sat | θ_res | α [1/m] | n | Ksat [m/s] | θ_fc(derived) | θ_wp(derived) |
|---|---|---|---|---|---|---|---|---|
| 1 | Sand | 0.43 | 0.045 | 14.5 | 2.68 | 8.25e−5 | ~0.06 | ~0.03 |
| 5 | Loam | 0.43 | 0.078 | 3.6 | 1.56 | 2.89e−6 | ~0.20 | ~0.11 |
| 8 | ClayLoam | 0.41 | 0.095 | 1.9 | 1.31 | 7.2e−7 | ~0.32 | ~0.19 |
| 11 | Clay | 0.38 | 0.068 | 0.8 | 1.09 | 5.6e−7 | ~0.36 | ~0.27 |

(The Campbell-option table — `ψ_sat`, `b` from Cosby-1984 — remains available for ED2-BC64 reproducibility runs.)

**Root vertical profile** — per-PFT exponential root fraction mapped onto the fixed layers (FATES delegation), normalized `Σ_k root_frac(k) = 1`:

```
root_frac(k) ∝ exp(root_beta·z_node(k))·dz(k)      (z_node ≤ 0 ⇒ decays with depth; root_beta from pft_table_t)
```

`root_beta` (or ED2's `root_depth`) and already-adjacent hydraulics traits (`root_sra`, `root_to_leaf_ratio`) live in `pft_table_t`; a per-cohort deepest-rooting-layer `kroot` derives from `root_depth` × allometric height. Peat (ED2 class 12) and bedrock (class 13, hydraulics-disabled) are config-selectable classes, **not** hard-coded `nsoil==` switches (ED2 quirk A6/A7).

---

## 8. Config `[hydrology]` / `[soil]` block

### 8.1 Keys (all REQUIRED unless noted DERIVED, presence-mapped, MAIN file — non-PFT settings)

```toml
[soil]
n_soil_layer   = 20            # n_soil_layer_active (≤ n_soil_layer_max); ignored if soil_layer_z given
soil_depth     = 8.5           # [m] column bottom (exponential-generator mode)
grid_growth    = 4.0           # exponential layer-spacing factor (exponential-generator mode)
# soil_layer_z  = [0.0, -0.02, -0.06, -0.12, ...] # OPTIONAL explicit interface ELEVATIONS [m, z≤0, descending];
#                                                # if present, OVERRIDES soil_depth+grid_growth (§3a)
texture_class  = "Loam"        # or per-layer array; selects the retention row
retention      = "van_genuchten"  # "van_genuchten" (default) | "campbell"  → SOIL_RETENTION_VG | _CAMPBELL
porosity       = 0.43          # θ_sat  (or per-layer)     REQUIRED
theta_res      = 0.078         # residual water content    REQUIRED
ksat           = 2.89e-6       # [m/s]                     REQUIRED
# --- van Genuchten (retention = "van_genuchten", DEFAULT) ---
vg_alpha       = 3.6           # [1/m] inverse air-entry   REQUIRED (vG)
vg_n           = 1.56          # pore-size index (>1)      REQUIRED (vG)
# --- Campbell/BC64 (retention = "campbell") ---
# psi_sat      = -0.26         # [m] air-entry potential   REQUIRED (Campbell)
# campbell_b   = 5.65          # Clapp–Hornberger exponent REQUIRED (Campbell)
psi_fc_mpa     = -0.033        # [MPa] field-capacity head (fixes θ_fc definition, S8)
# theta_fc / theta_wp / theta_cp are DERIVED in derive_config (§7) — NOT keys

[hydrology]
dtlsm_sec        = 900.0       # fast-loop sub-daily step (ED2 DTLSM)
solver           = "backward_euler"    # SOIL_SOLVER_BE (only supported; error stop otherwise)
linearization    = "frozen"    # "picard" (production) | "frozen" (MVP)
zeng_decker      = false       # P2 ONLY (retention-integral ψ_E); MVP uses plain-gravity flux (§3e)
bottom_bc        = "free_drain" # free_drain | aquifer | bedrock | slope → SOIL_BC_*
f_drai           = 2.5         # [1/m] SIMTOP baseflow decay (SOIL_BC_AQUIFER)
q_drai_max       = 5.5e-6      # [kg/m2/s]
f_over           = 0.5         # [1/m] saturated-fraction decay (Dunne runoff)
f_max            = 0.4         # max saturated fraction (CTI-derived)
w_pond_max       = 5.0         # [kg/m2] ponding capacity before surface runoff
dewmx            = 0.1         # [kg/m2 per unit PAI] per-cohort canopy storage capacity (§3c)
intercept_alpha  = 1.0         # Beer interception efficiency (per cohort)
intercept_k      = 0.5         # [1/PAI] Beer extinction for the interception fraction (§3c)
dsl_dmax         = 0.015       # [m] DSL max thickness (soil evap)
dsl_theta_init   = 0.8         # DSL initiation θ_1/θ_sat
substep_mode     = "adaptive"  # adaptive | fixed (GPU)
rtol             = 1.0e-3
atol             = 1.0e-4      # [m3/m3]
h_init           = 900.0       # [s] initial substep
max_substep      = 200
max_picard       = 5
debug_error      = false       # error stop on cap-hit / mass_resid in Debug builds
```

### 8.2 Wiring (mirrors the RT/hydraulics config precedent)

Declare the scalar/geometry fields (incl. the optional `soil_layer_z` explicit-elevation array) and **plain per-layer texture arrays** on `meds_config_t` (flat-block precedent of `[carbon]`/`[leaf]`). The **string→enum selectors are `integer(ik), parameter` codes declared in `meds_config` (shared)** — `SOIL_SOLVER_BE`, `SOIL_RETENTION_VG`/`SOIL_RETENTION_CAMPBELL`, `SOIL_BC_FREE_DRAIN`/`SOIL_BC_AQUIFER`/`SOIL_BC_BEDROCK`/`SOIL_BC_SLOPE`, `SOIL_LIN_PICARD`/`SOIL_LIN_FROZEN`, `SOIL_SUBSTEP_ADAPTIVE`/`SOIL_SUBSTEP_FIXED` (all `SOIL_`-prefixed, distinct from the plant `HYDRO_SOLVER_BE`) — re-exported alongside the existing `TS_*/BK_*/DIST_*` codes, so `load_meds_config` (in `meds_config_io`, which links shared+allometry, **not** biophysics) maps the TOML strings exactly like every existing `TS_MONTHLY`/`SM_LEUNING` selector (§2.1 blocker fix). Register keys via the typed `req_r/req_i/req_l/req_s` helpers under the MAIN-file section; the presence-check block then `error stop`s listing any missing key with no code change. Derive `soil_layer_z/z_node/dz/dz_node` and the per-layer texture arrays in `derive_config` (the `height_edges` pattern) as **plain real allocatables on `meds_config_t`**; the typed `soil_params_t` and the derived `θ_fc/θ_wp/θ_cp` are assembled later by the biophysics `derive_soil_params` (§2.1, §6). Validate in `validate_config` with the `tag//'…'` idiom: `if (cfg%soil_porosity <= 0.0_wp) error stop tag//'porosity <= 0'`; **per curve** — vG `vg_alpha > 0`, `vg_n > 1` (default) *or* Campbell `psi_sat < 0`, `campbell_b > 0`; `0 ≤ θ_res ≤ θ_wp ≤ θ_fc ≤ θ_sat` (checked on the derived thresholds), `n_soil_layer_active ≤ n_soil_layer_max`.

**Constants to ADD to `meds_constants`** (genuine universal constants, so allowed there — never TOML):

| constant | value | why |
|---|---|---|
| `rho_h2o` | 1000 kg/m³ | kg/m² ↔ m ↔ m³/m³ conversions; unbundle `grav_head` |
| `grav` | 9.80665 m/s² | standalone g for head & `α_soil` (only `ρ_w·g` fused as `grav_head` exists today) |
| `latent_heat_vap` | 2.501e6 J/kg | transpiration/evap mass ↔ latent-heat flux (future energy balance) |
| `r_wv` | 461.5 J/kg/K | water-vapor gas constant for Philip `α_soil` |

Existing `grav_head = 9.804e-3 MPa/m` (= `ρ_w·g`) is the exact soil-head→MPa converter for the plant seam (§9).

---

## 9. Coupling contracts

### 9.1 Closes the plant-hydraulics `psi_soil` boundary condition (the headline)

`hydro_env_t` (read-only per-plant BC into `solve_plant_water`) has the stubs this module fills:

| field | units | role | supplied by |
|---|---|---|---|
| `env%soil_psi` | **MPa** | aggregated rhizosphere water potential — **the BC to close** | §9.2 |
| `env%rhizo_cond` | **kg/s/MPa per plant** | soil→root conductance per plant | §9.2 |
| `env%transp` | kg/s | transpiration demand E (sink into hydraulics) | leaf/fast loop |

The soil module exports `flux%psi_soil(k)` [MPa] and `flux%soil_cond(k)` **[kg/s/MPa per plant]** (unambiguously per-plant, matching `g_rhizo_k` below — not a per-RAI basis), and the orchestration layer (in `meds_aux`, the `meds_vegetation_dynamics` analogue) aggregates them into the per-plant scalar pair using the existing `rhizosphere_cond(soil_cond, broot, sra, root_frac, dz, nplant)` hook (labelled "Production fills `rhizo_cond` from a future soil module" — **that future module is this one**):

```
psi_soil_k = grav_head · psi_head_k[m]                          [MPa]     retention curve (vG/Campbell) → MPa
g_rhizo_k  = rhizosphere_cond(K_k, broot, sra, root_frac_k, dz_k, nplant)  [kg/s/MPa per plant]
env%rhizo_cond = Σ_k g_rhizo_k                                              (parallel conductances, per plant)
env%soil_psi   = Σ_k g_rhizo_k · psi_soil_k / Σ_k g_rhizo_k                 (Thévenin-equivalent BC)
```

The **conductance-weighted mean** is the correct single-BC reduction of the parallel per-layer network: a lumped 2-node hydraulics solve driven by `(env%soil_psi, env%rhizo_cond)` reproduces the network's total uptake to first order, and the per-layer weights `w_k = g_rhizo_k/Σg_rhizo` are reused to redistribute the returned uptake (§3f) — so the BC and the sink partition are **mutually consistent**. When the per-soil-layer fine-root extension lands (hydraulics §16), `hydro_env_t` carries `soil_psi(:)`/`rhizo_cond(:)` arrays directly and the aggregation drops out — this design accommodates that with no core rework.

### 9.2 Consumes the transpiration sink

The hydraulics output **`hydro_flux_t%root_uptake`** [kg/s per plant] — the **plant** hydraulics type (distinct from this module's `vhydro_flux_t`), returned as `(dw_l + dw_w)/dt + e_transp`, "the budget term" — is scaled `×nplant`, summed over the patch's cohorts, and partitioned by `w_k` into the volumetric *demand* `U_total·w_k`, then ψ-limited and diagonal-folded into `S_k` (§3f). **Unit discipline:** hydraulics works in per-plant kg/s, MPa; the soil column in per-ground-area kg/m²/s, m³/m³, m-head. The orchestration layer does the `×nplant` and the m³/m³ ↔ kg/m² ↔ MPa conversions (via `rho_h2o`, `grav_head`) — the kernel never crosses unit systems.

### 9.3 Feeds surface optics, phenology, drought demography

- **Soil albedo** (RT §5, `meds_surface_optics`): `col%theta(1)` → moisture-darkening of the bare-soil reflectance (`surface_state_t` currently returns dry bare soil); the `w_surface`/snow branch fills the RT B6 stub. Data handoff: the aux layer writes `col%theta(1)` into `surface_state_t%theta_top` before the RT call.
- **Phenology WATER cue** (`meds_pheno_engine`): the per-PFT TEMP/WATER/HYDRO cue bitmask consumes **`ψ_soil` (MPa) at rooting depth** as the WATER limiting signal — `ψ` (not `θ`) is chosen because it is texture-normalized and already the currency of the HYDRO cue, so the two cues share units. **Seam:** the aux layer writes a per-PFT `pheno_env_t%psi_root [MPa]` (the `root_frac`-weighted mean of `flux%psi_soil`) that `meds_pheno_engine%update_water_cue(psi_root, …)` consumes — currently stubbed.
- **Future drought demography**: `ψ_soil` / integrated water stress feeds the Camac mortality hazard and a fire-threshold `θ_fr` diagnostic — the mechanistic replacement for the structure-only competition path.

### 9.4 The fast-loop drive (reserved)

A new fast-loop driver in `src/driver` (compiled into `meds_aux`), called from `advance_one_step` **before** `vegetation_dynamics`, runs an inner sub-daily loop over `cfg%hydro_dtlsm_sec`. Per patch, per sub-step: `canopy_radiation` → absorbed PAR → `leaf_gas_exchange` → per-cohort transpiration → `×nplant` summed to `forcing%root_uptake`; **`intercept_canopy_layer` swept top→bottom** over the height-sorted cohorts (met precip → `cohort%leaf_water`, bottom throughfall → `forcing%precip_ground`) → `vertical_hydrology_flux` updates `patch%soil_water(:,ip)` → aggregate to `env%(soil_psi, rhizo_cond)` (§9.1) → `plant_water_flux` advances `psi`. In this weave both the plant `hydro_flux_t` and this module's `vhydro_flux_t` are in scope; the rename (§4) lets the aux layer `use` both without aliasing. The **still-absent piece is a met forcing reader** (precip, SW/LW, air T, VPD → a new `src/io` source feeding `rad_forcing_t` + `vhydro_forcing_t`); until it exists the fast loop is disabled and this module stands on its unit test alone (like RT, leaf, hydraulics before wiring). This module lands alongside the `[hydraulics]` flatten-from-`cfg%pft` wrapper the interface header already flags as pending.

---

## 10. GPU / nvfortran portability

- **Two-layer structure:** the seam `vertical_hydrology_flux` takes derived types (host); the interior arithmetic `soil_water_be_step(theta, dz, dz_node, ksat, psi_sat, b, theta_sat, s_k, …, nzg)` takes **bare fixed-size arrays + `firstprivate` scalar params**, no derived types (`soil_opts_t` stays out of it) — the `growth_step` precedent (the device cannot read host module vars, so texture/geometry flow in as arguments, not module globals). The per-column sweep is a natural `!$omp target` candidate; **the parallel axis is columns (patches/sites), never within a column** — Thomas is a sequential recurrence, so parallelize across the many columns, one sweep per thread. Per-patch orchestration stays **host-only** (all restructuring is host-only).
- **Fixed `n_soil_layer_max`** ⇒ no allocatables/runtime shapes; scratch (`a,b,c,d,cp,dp`) as fixed-length automatics, stack-friendly on device.
- **Issue #7 trap:** `thomas_solve` is a subroutine with `intent(out) x(n_soil_layer_max)` — never `call foo(thomas(...))`. Bind every array result to a named array before any call. **Build nvfortran multicore on the new module** — a green ifx suite (only an `arg_temp_created` remark) is *not* sufficient; issue #7 was found in `src/biophysics` itself.
- **FPE-safe in Debug** (`-fpe0`/`-Ktrap=fp`): `Θ` clamp, `K` floor, `C`-near-saturation guard, clamped closed-form retention inverse — no `0/0`, no `(neg)**real`. Constitutive kernels `pure`/`elemental`, branch-light smooth limiters instead of `case`/`zero_flow` traps.
- **Do NOT use `-stdpar=gpu`** (managed-allocator double-free on the allocatable-component `site`); OpenMP `target` + `-gpu=mem:separate`, `MEDS_GPU=multicore|gpu`.
- **`SOIL_SUBSTEP_FIXED`** for the offload path (variable `nsub` is warp divergence).

---

## 11. Validation & milestones

CTest `test_vertical_hydrology` (drives synthetic `vhydro_forcing_t`/`soil_params_t` from `defaults()`, marches the kernel, asserts):

1. **Mass conservation (the strongest invariant):** the total-store balance `Δ(all stores) = dt·(precip − evap − uptake − runoff − drainage)` closes to **machine precision for the conservative flux-divergence `θ` update (§5.2)** once `clip_excess` and `uptake_deficit` are included as budget terms; `flux%mass_resid ≈ 0` in wet and dry columns, in **both** linearizations. In the **frozen** path additionally assert the stored `θ` sits off the retention curve by no more than the step-doubling-bounded linearization residual (and that `ψ_soil` is re-derived from that `θ`); in the **Picard** path assert `θ` reconciles to `θ(ψ)` within iteration tolerance.
2. **Constitutive round-trip (both curves):** `theta_of_psi(psi_of_theta(θ)) == θ` to round-off for **van Genuchten and Campbell**; vG is smooth through air-entry, the Campbell `min(1,·)` clamp holds `θ ≤ θ_sat` for `ψ ≥ ψ_sat` (no supersaturation); `C(ψ)` matches finite-difference `dθ/dψ`; `K,ψ` monotone, FPE-clean at `θ_res` and `θ_sat` (incl. vG's `C→0` at saturation via the `C`-floor).
3. **Analytic / steady-state limits (§3e):** (a) **unit-gradient free drainage** under `SOIL_BC_FREE_DRAIN` → `dψ/dz = 0`, `q = K` steady state (MVP plain gravity); (b) **hydrostatic `q→0`** under `SOIL_BC_BEDROCK` (no-flow), `S=0`, initial `ψ=ψ_E` → the MVP plain-gravity flux exhibits the coarse-layer drainage **drift** (S9, the baseline), and the **P2 retention-integral ZD drives `q ≡ 0` for any `dz`** (the ZD acceptance test); (c) zero-flux column with `S=0` → `dθ = 0`.
4. **vs a fine explicit RK4 reference:** (a) BE-adaptive matches a 1000-substep explicit reference to `rtol` in wet **and** dry columns; substep count drops when quiescent and refines at an infiltration front; no NaN/trap at saturation. (b) **cap-hit case:** a sharp front forced against `max_substep` (and, in Picard, `max_picard`) sets `flux%converged=.false.`, logs a diagnostic, and `error stop`s under `debug_error=.true.` (§5.3).
5. **The coupling check (closes the loop):** the column's aggregated `env%soil_psi`/`rhizo_cond` (§9.1) reproduce the hydraulics steady state `psi_wood → soil_psi − E/rhizo_cond` (the `test_plant_hydraulics` assertion) — the soil↔plant seam verified end-to-end.
6. **Interception cascade:** per cohort `throughfall + drip + Δleaf_water + E_canopy·dt = rain_above·dt`, and the patch total `Σ_cohort` closes to canopy-top precip; capacity `dewmx·PAI` respected; drip only above capacity; the top→bottom sweep gives a real cascade (taller cohorts shade shorter).
7. **Infiltration/ponding partition:** infiltration ≤ Ksat capacity; `q_inf,max ≥ K(θ_1) > 0` and *increases* as the surface dries (suction sign); excess accumulates in `w_surface` then runs off; total closes.
8. **ED2 comparison (regression):** drive ED2's soil column on a single texture class with free drainage and a prescribed root sink; MEDS matches the equilibrium θ-profile and drainage flux **within the calibrated tolerance of test 4 (`rtol = 1e-3` on θ-profile, 2 % on integrated drainage flux over the run)** — isolating the intended numerics change from the shared CH closure. (Open Q6 tracks re-calibrating these once a fine reference exists.)
9. **nvfortran multicore build green** (portability gate); serial↔multicore↔GPU determinism (the 7/7 cross-backend discipline).

**Phased milestones:**

| Phase | Deliverable | Tests |
|---|---|---|
| **P0** | Types + CMake lib + seam skeleton; `meds_soil_parameters` **vG + Campbell** `ψ/θ/K/C` + closed-form inverses + `derive_soil_params`; hand-rolled `thomas_solve` (subroutine) | 2 |
| **P1 (MVP / production default)** | flux-form BE + frozen-coefficient single-Newton + upstream `K` + **plain-gravity flux** + ψ-limited diagonal root sink + capped-infiltration/free-drainage BCs + **per-cohort interception cascade** + DSL soil evap + adaptive step-doubling + conservative-θ ΔW budget + cap-hit contract | 1,3,4,6,7 |
| **P2** | Celia modified-Picard + **retention-integral Zeng–Decker (interior faces)** + Neumann→Dirichlet ponding switch + SIMTOP aquifer/water-table bottom BC + Dunne `f_sat` runoff + van Genuchten option | 1,3,4,8 |
| **P3** | SoA `soil_water(:,:)`/stores through the **6** patch lockstep sites + `leaf_water` through the **cohort lockstep** (incl. `fuse_2_patches` blend + donor-copy disturbance init) + TOML `[soil]`/`[hydrology]` + `derive_config` + `derive_soil_params` + `validate_config`; export `psi_soil`, close the hydraulics BC | 5, full build |
| **P4** | fast-loop driver in `meds_aux` (once met forcing lands); leaf→hydraulics→soil weave | 4b + **P4 smoke: closed site water balance over a driven synthetic day** |
| **P5** | nvfortran GPU parity; `SOIL_SUBSTEP_FIXED`; per-soil-layer fine-root nodes (hydraulics §16); **C-API + Python (end of biophysics dev)** | 9 |

---

## 12. How this differs from ED2 / FATES / CLM

MEDS **owns the soil column** (unlike FATES, which delegates to a host — MEDS *is* its own host), **adopts CLM's numerics and boundary conditions**, **keeps ED2's CH closure**, and **connects root uptake via the FATES delegation seam** to its own hydraulics kernel. **Note on the columns below:** CLM5 and FATES are split deliberately because they contribute *different kinds* of thing — **CLM5 supplies the solver + boundary conditions** (it has a native soil-hydrology solver), while **FATES has no native soil solver at all** (it delegates soil water to a host land model) and instead supplies the **host-delegation / per-layer uptake contract** (`Q_j`, `*Mem`/compute split). Conflating them would hide that MEDS borrows the *numerics* from CLM and the *coupling seam* from FATES.

| Aspect | **ED2** | **CLM5** (solver + BCs) | **FATES** (delegation/uptake contract) | **MEDS (this design)** |
|---|---|---|---|---|
| Time integration | explicit adaptive RK4, co-integrated with canopy energy (stiff, GPU-hostile) | implicit backward-Euler tridiagonal, one long `dt` | *(delegates to host — none native)* | **implicit BE-Thomas**, operator-split from canopy, stateless over fixed `n_soil_layer_max` |
| Equilibrium correction | none (coarse-layer drainage drift → needs fine substeps) | Zeng–Decker (2009) | *(host)* | **Zeng–Decker at P2** (retention-integral `ψ_E`, interior faces); **MVP plain gravity** |
| Retention/conductivity | Campbell/BC64 (default) + `vG80` option | Clapp–Hornberger | *(host)* | **van Genuchten–Mualem default; Campbell/BC64 option** — ψ-coupled, curve-independent |
| Surface infiltration | **none direct** — buffered in virtual/TSW pool, dumped to sat-room between steps | `q_infl` = throughfall − Dunne runoff, top BC | *(host)* | **conductivity-limited top Darcy flux** + ponding + Dirichlet switch |
| Bottom BC / drainage | free-drain / bedrock / aquifer (`ISOILBC`) | SIMTOP baseflow + unconfined aquifer + water table | *(host)* | **free-drain MVP → SIMTOP aquifer** (config) |
| Runoff | e-folding pond removal + TOPMODEL lateral (off by default) | Dunne `f_sat` saturation-excess + Horton | *(host)* | **Horton MVP → Dunne `f_sat`** |
| Root uptake | single-`kroot` deposit + upward availability redistribution | CLM `β·E` empirical partition | **per-layer sink `Q_j`, host owns water / veg owns uptake** | **FATES delegation** — consume hydraulics per-plant uptake, ψ-limited `S_k` on diagonal |
| Soil evaporation | Lee–Pielke `β` + `ed_grndvap8` MH91 forms | Swenson–Lawrence DSL resistance | *(host)* | **Swenson–Lawrence DSL adopted** |
| Interception | per-cohort TAI split, no cascade, Deardorff `sigmaw`, post-step drip | `tanh(L+S)` fraction + `dewmx·(L+S)` bucket | *(host)* | **per-cohort capacity-limited bucket + real top→bottom cascade**; `σ_w` exported, evap → canopy-air-space module |
| Mem/compute split | monolithic stateful | grid state + physics | **`*Mem` state / stateless compute** | **stateless kernel; state in `meds_state` SoA** (FATES-style) |
| Platform | stateful, CPU | host land model, CPU | host land model, CPU | **stateless kernel; CPU + GPU** (`pure`/`elemental`, OpenMP target) |

**One line:** MEDS = CLM's implicit solver + boundary conditions (BE-Thomas + SIMTOP + Swenson–Lawrence; Zeng–Decker at P2) wrapped around **van Genuchten** soil physics (Campbell/BC64 optional for ED2 reproducibility), connected to MEDS's own FATES-HYDRO-like hydraulics kernel via FATES's host-delegation seam — GPU-friendly and numerically robust, ψ-coupled (curve-independent).

---

## 13. Bugs / quirks found in the ED2 reference

Curated from the vertical-soil-water and canopy-interception extractions; flagged to fix-by-construction or carry deliberately.

**S1 — No direct rain infiltration flux (fix).** Precip/dew/drip buffer in the virtual/TSW pool and only enter the soil (up to saturation room) in `adjust_sfcw_properties` *between* RK4 steps; the derivative has no conductivity-limited surface flux (only ground evaporation `w_flux_g(mzg+1)`). **MEDS uses a proper top Darcy/Neumann flux** (§3d).

**S2 — `infiltration_method/=0` is dead code.** The alternate conductivity-based infiltration immediately calls `fatal_error('Running alt infiltation when we shouldn''t be')`. Ignored.

**S3 — `slpotcp` vs `slpotwp` fragile unit coincidence.** Written with different formulas that agree *only because* `wdns=1000` (`1e6/wdns = 1000 = wdns`). **MEDS always uses `ψ[m] = ψ[MPa]·1e6/(ρ_w·g)`** with explicit `rho_h2o`, `grav`.

**S4 — BC64 `soilre=0` ⇒ `ψ→−∞` as θ→0.** Bounded only by `soilcp`/`soilwp` clamps + the `drysoil` limiter. **MEDS makes the clamp explicit** (`Se` clamped in the elemental curve, `min(1,·)` on the Campbell inverse) rather than relying on limiters elsewhere. The **vG default uses `θ_res > 0`**, so `ψ` stays finite as `θ→θ_res` — this divergence is a Campbell-option concern only.

**S5 — Anti-overshoot limiter strands water (do not port as-is).** Zeroing `w_flux_g` on dry/sat neighbors is not flux-conservative and can freeze exchange between a saturated and an adjacent full layer. **MEDS handles overshoot via the implicit solve + `θ`-clip with excess routed to ponding and counted** (`clip_excess`, conservative).

**S6 — Bedrock hard-coded `nsoil==13` in multiple places** (flux zero, transpiration skip, `slcpd`), parallel to the `'BDRK'` method tag — two switches. **MEDS uses one config `bottom_bc = bedrock` / a per-class `hydraulics_off` flag**, no magic index.

**S7 — Peat (class 12) fully hard-coded to unverified LEAF-3 values** (source comment warns ED-2.2 legacy differs). **MEDS drives every class from the TOML retention table**, no hard-coded overrides; any pedotransfer coefficients are run-loaded, not literals.

**S8 — Field-capacity definition is scheme-inconsistent** (K=0.1 mm/day vs ψ=−0.01 MPa vs ψ=−0.033 MPa across schemes). **MEDS fixes one definition per config** (`psi_fc_mpa`, default −0.033 MPa) and derives `θ_fc` from the curve; the definition is documented and validated (§7).

**S9 — Free-drainage over-drains in long dry spells** (ghost layer mirrors `klsl` moisture, so there is always downward gravitational loss). **MEDS removes the spurious drainage with the P2 retention-integral Zeng–Decker correction** (interior faces only, §3e); the MVP plain-gravity flux still exhibits this coarse-layer drift (the baseline validated in §11 test 3b), and the fix is deferred to P2 with the aquifer BC — *not* claimed for the MVP. (The earlier "ZD with `z_wt` at the column bottom" framing was dropped: a linear-midpoint `ψ_E` makes ZD algebraically identical to plain gravity, and `z_wt` cancels from every interior flux.)

**S10 — `freezecoef = 7·ln(10)` is an extremely steep frozen throttle** (~7 orders from `fliq=1→0`). Carried as a **hook only** (v1 is `fliq≡1`); to be re-derived when soil thermal lands.

**B(interception)-1 — Wet-fraction / transpiration double surface.** ED2 scales film evaporation by `sigmaw` but does **not** scale transpiration by `(1−sigmaw)` — fully-wet leaves transpire *and* evaporate at full rate over the same LAI, bounded only by the `can_rhv≥1` limiter. MEDS already carries the per-cohort `leaf_water` film and exports `σ_w` (§3c); the **`(1−σ_w)` transpiration factor lands with the canopy-air-space module** (which owns interception evaporation and the transpiration gradient), not in this soil column.

**B(interception)-2 — No inter-cohort vertical rain cascade.** ED2 interception is a TAI-weighted split of one patch-total `intercepted_max`; taller cohorts do not shade shorter ones. **MEDS supersedes this with a real top→bottom cascade** (§3c): each cohort intercepts from the throughfall of those above, so taller cohorts shade shorter ones — an improvement over ED2, not a faithful reproduction of its limitation.

**B(interception)-3 — Post-step instantaneous shedding, not a rate.** In-derivative `wshed≡0` ("TURNING OFF SHEDDING FOR NOW"); `leaf_water` legally overshoots capacity mid-step and is capped in `adjust_veg_properties`. **MEDS's bucket drip is a proper flux** (`q_drip`), so the solver sees it.

---

## 14. Phased implementation plan (MVP → full)

- **P0 — foundation.** `meds_biophysics_types` extended with `soil_column_t`, `vhydro_forcing_t`, `vhydro_flux_t`, `soil_params_t`, `soil_opts_t`, `n_soil_layer_max` (enums live in `meds_config`); `meds_soil_parameters` (`pure`/`elemental` **van Genuchten + Campbell** `ψ/θ/K/C` + closed-form inverses + `derive_soil_params`); `meds_soil_solver` (subroutine, `intent(out)`). Tests: constitutive round-trip (both curves), `dθ/dψ` vs finite difference. CMake lib builds; nvfortran multicore green.
- **P1 — MVP / production default.** `vertical_hydrology_flux` seam + bare-array `soil_water_be_step`: flux-form BE, frozen-coefficient single-Newton, upstream `K`, **plain-gravity flux**, ψ-limited diagonal root sink (prescribed-profile forcing), capped-infiltration + ponding + free-drainage bottom, **per-cohort capacity-limited interception cascade** (`intercept_canopy_layer` + top→bottom sweep), DSL soil evaporation, adaptive step-doubling, conservative-θ ΔW budget, cap-hit contract. Tests 1,3,4,6,7. Standalone on its CTest.
- **P2 — richer physics.** Celia modified-Picard; **retention-integral Zeng–Decker equilibrium correction** (interior faces only, bottom-flux BC kept independent); van Genuchten option; Neumann→Dirichlet ponding switch; SIMTOP aquifer + water-table `z_wt` bottom BC; Dunne `f_sat` runoff. Tests 1,3,4,8.
- **P3 — state + config + coupling.** `soil_water(:,:)`/`surface_water` through the **6** patch lockstep sites + `leaf_water` through the **cohort lockstep** (fuse/split/recruit/terminate) + area-weighted fusion blend (`fuse_2_patches`) + donor-copy disturbance init; `[soil]`/`[hydrology]` TOML (incl. optional `soil_layer_z`) + presence map + `derive_config` (plain arrays) + `derive_soil_params` (typed) + `validate_config`; add `rho_h2o`/`grav`/`latent_heat_vap`/`r_wv` to `meds_constants`; declare `SOIL_*` enum parameters in `meds_config`; export `flux%psi_soil` and close `hydro_env_t%(soil_psi, rhizo_cond)` via the §9.1 aggregator. Test 5 (the coupling check), full engine build.
- **P4 — fast-loop wire (blocked on met forcing).** Fast-loop driver in `meds_aux`, called from `advance_one_step` before `vegetation_dynamics`; RT→leaf→soil→hydraulics weave over `dtlsm_sec`. Needs the met reader (`src/io`). **Validation:** the cap-hit case (test 4b) plus a **P4 smoke/integration assertion** — a closed site water balance over a driven synthetic day (all stores + boundary fluxes reconcile to `atol`), so the milestone is not test-free.
- **P5 — platform + extension.** nvfortran GPU parity + `SOIL_SUBSTEP_FIXED`; per-soil-layer fine-root nodes (drops the §9.1 aggregation — hydraulics §16); soil thermal/freezing; snow/sfcwater stack (lifting the liquid-water-only restriction of §1.2); **C-API + Python (`bind(c)` shim + `-DMEDS_BUILD_PYLIB`) — deferred to the end of biophysics development, not created now.**

---

## 15. Open questions

1. **`n_soil_layer_max` and grid.** Fix `n_soil_layer_max = 20` with exponential spacing to 8.5 m (CLM-like), or keep an ED2 coarse-12 default for regression? *Recommend:* `n_soil_layer_max = 20` compile-time ceiling, `n_soil_layer_active` + `soil_depth` + `grid_growth` from config so both regimes reproduce.
2. **Zeng–Decker in the MVP.** *Resolved (§3e):* the MVP carries **no** ZD. As written with a linear-midpoint `ψ_E` it is algebraically identical to plain gravity and buys nothing, and its `z_wt`-at-bottom rationale was mechanically wrong (`z_wt` cancels from interior fluxes). The MVP uses the plain-gravity flux; **real ZD — retention-integral `ψ_E`, interior faces only, bottom-flux BC independent — is deferred to P2** with the aquifer BC, where `z_wt` becomes physically active. Remaining choice: retention-integral quadrature order for `ψ_E,k` (midpoint vs 2-point) — calibrate against test 3b coarse-layer drift.
3. **Interception granularity.** *Resolved (§3c):* per-cohort `leaf_water` film with a top→bottom cascade is the MVP. The wetted-fraction *evaporation* and the `(1−σ_w)` transpiration split (bug B-1) await the canopy-air-space module; this kernel already exports `σ_w`.
4. **`ψ_soil` aggregation.** Conductance-weighted Thévenin mean (§9.1) vs availability-weighted vs deepest-layer. *Recommend:* conductance-weighted (consistent with the sink partition), superseded by per-layer `soil_psi(:)` when the fine-root extension lands.
5. **van Genuchten vs Campbell default.** *Resolved (§3b):* **van Genuchten–Mualem is the MVP default** (the direction ED2 is moving; better fit to data; smooth at air-entry). Campbell/BC64 is a config option for ED2-BC64 reproducibility. The coupling is ψ-based (curve-independent), so this is a free choice. Watch: vG's steep near-saturation `K` and `C→0` at saturation — the `C`-floor + adaptive substep handle it, and the P2 Picard path helps if a case stiffens.
6. **Numerics tolerances/defaults:** `atol = 1e-4 m³/m³`, `rtol = 1e-3`, `max_substep = 200`, `max_picard = 5`, `dtlsm_sec = 900`? To calibrate against the fine explicit reference (test 4); the calibrated `rtol`/drainage tolerance then pins the test-8 ED2 regression bound.
7. **Surface runoff timescale.** Instantaneous Horton (MVP) vs ED2's `runoff_time = 3600 s` e-folding pond removal vs Dunne `f_sat`. *Recommend:* instantaneous Horton P1 → Dunne P2; the ED2 e-folding form is a config fallback.
8. **Aquifer `f_max` / CTI map** source for `SOIL_BC_AQUIFER` (needs a topographic-index parameter) — deferred with the aquifer BC; free-drainage needs none.

---

## 16. References

- **Zeng & Decker (2009)** *J. Hydrometeorol.* 10:308 — equilibrium-potential-corrected Richards; the coarse-layer drainage-drift fix. **Deferred to P2 (§3e)** — the fix requires a *retention-integral* `ψ_E`; with a linear-midpoint `ψ_E` it is algebraically identical to plain gravity, so the MVP does not carry it.
- **Swenson, Lawrence & Lee (2012)** *JGR* 117, D21107 & **Sakaguchi & Zeng (2009)** *JGR* 114, D01107 — dry-surface-layer / soil-resistance evaporation lineage (the exponential DSL ClimaLand uses). **Soil-evaporation alternative (§3g).**
- **Deck, Kohler et al. (2026)** *J. Adv. Model. Earth Syst.* — ClimaLand (CliMA Land, Julia): `α_soil` pore-space RH + SZ09 exponential DSL + Monin–Obukhov coupling; a modern cross-check for the soil-evaporation kernel. **Cross-reference (§3g).**
- **Celia, Bouloutas & Zarba (1990)** *Water Resour. Res.* 26:1483 — mass-conservative modified-Picard mixed-form Richards. **Production linearization (§5.2).**
- **Lawrence et al. (2019, 2020)** CLM5 Description / Technical Note (NCAR) — implicit tridiagonal soil hydrology, SIMTOP aquifer, `tanh(L+S)` interception, DSL evaporation; the ZD discrete-flux form is verified against this tech note (§3e). **Numerics + BCs reference (§3, §5).**
- **Niu et al. (2005)** *JGR* 110, D21106 (SIMTOP `f_sat`/baseflow) & **Niu et al. (2007)** *JGR* 112, D07103 (unconfined aquifer). **Bottom BC / runoff (§3d, §3e).**
- **Swenson & Lawrence (2014)** *JGR* 119:10299 — dry-surface-layer soil-resistance evaporation. **Soil evaporation (§3g).**
- **van Genuchten (1980)** *Soil Sci. Soc. Am. J.* 44:892 & **Mualem (1976)** *Water Resour. Res.* 12:513 — the van Genuchten retention curve + Mualem unsaturated-conductivity model, the **default closure**; **Carsel & Parrish (1988)** *Water Resour. Res.* 24:755 — texture-class vG parameters (the default table, §7). **Default closure + parameters (§3b, §7).**
- **Clapp & Hornberger (1978)** *Water Resour. Res.* 14:601; **Campbell (1974)** *Soil Sci.* 117:311; **Cosby et al. (1984)** *Water Resour. Res.* 20:682 — CH/Campbell retention & the sand/clay pedotransfer, the **Campbell-option** generator (§7). **Option closure + parameters (§3b, §7).**
- **Christoffersen et al. (2016)** *Geosci. Model Dev.* 9:4227 & **FATES-HYDRO V1.0, Xu/Christoffersen et al. (2023)** *GMD* 16:6267 — host-delegates-water, vegetation-owns-uptake seam; per-layer `Q_j`; FATES has no native soil solver (host-delegated). **Root-uptake contract (§3f, §9, §12).**
- **Longo et al. (2019)** *GMD* 12:4309 — ED-2.2 technical description; the soil-water scheme this module reimplements.
- **ED2 source** (`../ED2/ED/src`): `memory/soil_coms.F90`, `init/ed_params.f90` (`init_soil_coms`), `dynamics/rk4_derivs.f90` (`leaftw_derivs`), `dynamics/rk4_misc.f90` (`adjust_sfcw_properties`, `adjust_veg_properties`), `utils/ed_therm_lib.f90` (`ed_grndvap8`), `dynamics/lsm_hyd.f90` — the extracted reference (§3, §13).
- **MEDS internal:** `archive/MEDS_HYDRAULICS_DESIGN.md` (§2/§5/§6.3/§11/§16 — the stateless + adaptive-step-doubling + boundary-only-conservation idiom this module reuses); `archive/radiative_transfer_design.md` (the `meds_biophysics` seam/DAG precedent, incl. `derive_rad_optics` filling `rad_pft_optics_t`); `CLAUDE.md` (state/process wall, no-hard-coded-parameters, naming convention, issue #7 nvfortran trap).
