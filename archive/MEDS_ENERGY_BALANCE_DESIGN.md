# MEDS Energy Balance — Module Family Design

A **stateless** biophysics module family for MEDS (`src/biophysics/`) owning the land-surface **thermal energy budget**: the internal-energy balance of the four coupled stores that carry temperature in an ecosystem-demography column — **leaf and wood** (per cohort), the **canopy air space (CAS)** (per patch), the **soil thermal column** (per patch, per layer), and the **ground/surface skin** (the top soil layer, no separate skin state). It is the thermal twin of the column-hydrology module: same negative-z soil geometry, same `meds_soil_solver` Thomas sweep, same stateless-kernel + state/process-wall discipline, same "forcing passed as a value type, prognostic state in `meds_state`" split. Where the hydrology module moves **water** through the column and closes the plant-hydraulics `soil_psi` BC, this family moves **energy** through the same column and closes every **forced-temperature seam** the rest of MEDS currently hard-wires: `leaf_env_t%leaf_temp`, `wood_env_t%wood_temp`, `chydro_forcing_t%t_ground`, the two respiration/phenology `soil_temp` stubs, the two-stream longwave **lower** BC `surface_state_t%soil_temp`, and the two-stream **canopy** leaf/wood emission temperature `canopy_radiation(…canopy_temp…)` — **nine seams in all** (§9).

The physics reference is **ED2** (`leaftw_derivs`, the RK4 soil/CAS/vegetation energy budgets) revised toward **CLM5 numerics** (implicit backward-Euler Thomas soil heat, analytic Clausius linearization) and the **ClimaLand internal-energy state design** (prognostic volumetric internal energy, closed-form temperature read-off, phase change as a state inversion rather than a solve-then-adjust). The single load-bearing state decision — **prognose internal energy / enthalpy, never temperature, wherever water can change phase** — is taken directly from ED2's `uextcm2tl` pattern and ClimaLand's `ρe_int`, and is what makes freeze/thaw a smooth invertible read-off and makes leaf-flush / water-mass changes a re-partition rather than a re-solve.

**Per the user directive, this family is STATELESS-FIRST.** Each of the four per-store kernels is designed and unit-tested standalone now, taking the *other* stores' temperatures as **forced inputs** (exactly the seams in §9). The coupled fixed point — the simultaneous leaf ↔ CAS ↔ ground ↔ soil solve that ED2 and CLM iterate every sub-step — is **deferred to P3 orchestration**, precisely as the column-hydrology module deferred its `psi_soil` coupling and fast-loop wiring. Nothing in the kernels changes when the coupling lands; only the driver stops lagging.

Units are SI throughout: **internal energy `[J/m³]` (soil, volumetric) and `[J/m²ground]` (leaf/wood, extensive), specific enthalpy `[J/kg]` (CAS), temperature `[K]`, heat flux `[W/m²]`**. Kinds `wp/ik`, `_wp` literals, `implicit none`, `pure`/`elemental` constitutive + inversion kernels, `error stop`, ≤132 columns, nvfortran-safe.

**Deferred siblings (out of scope here).** The Monin–Obukhov turbulence/roughness closure that produces the aerodynamic conductances (`ggnet`, `ustar`, `r_aero`, leaf `gbh`/`gbw`) is a *separate* biophysics module; this family consumes those conductances as forcing. Snow / temporary-surface-water thermal stack is deferred to P2b (mirrors hydrology's water-only-v1). The coupled surface fixed point is P3.

---

## 1. Scope, target processes & reachability

### 1.1 The four energy stores and their prognostic representation

| Store | Level | Prognostic variable | Recover | Rationale |
|---|---|---|---|---|
| **Leaf** | cohort | `leaf_energy` [J/m²ground] | `(leaf_temp, leaf_fliq)` | leaf-held water (hydraulic + interception film) can freeze; ED2 `leaf_energy` |
| **Wood** | cohort | `wood_energy` [J/m²ground] | `(wood_temp, wood_fliq)` | sapwood/bark water freezes; ED2 `wood_energy` |
| **Canopy air space** | patch | `can_enthalpy` [J/kg] (+ `can_shv` [kg/kg]) | `can_temp` | air has no phase change, but enthalpy mirrors ED2, is condensation-safe, and fuses the sensible+latent split into one budget |
| **Soil thermal column** | patch, per layer | `soil_energy(k)` [J/m³] | `(soil_temp(k), soil_fliq(k))` | soil water freezes; latent plateau is the whole point; ED2 `soil_energy` |
| **Ground / surface skin** | patch | *(none — is `soil_energy(1)`)* | `t_ground = soil_temp(1)` | ED2 has no separate skin; the surface energy balance is a **flux BC on the top soil layer** |

**The governing decision (Bundle E, ED2, ClimaLand): internal energy, not temperature, wherever phase change lives.** Internal energy is the conserved quantity; every flux (radiation, sensible, latent, sap-flow, advected soil water) is an additive `[W/m²]` source/sink, so the budget is a clean sum with no products of state variables. Crucially it makes freeze/thaw a **smooth, monotone, invertible** map: as energy accumulates through the mixed-phase plateau, temperature pins at the triple point while the *liquid fraction* absorbs the latent heat — impossible if temperature were prognostic (the map `T(u)` is multivalued at fixed `T` on the plateau, single-valued in `u`). It also makes heat-capacity changes (leaf flush/drop, water mass changing, soil moisture changing) a **non-event for the state variable** — you re-partition energy, you do not re-solve temperature. This is ED2's `leaf_energy`/`wood_energy`/`soil_energy` design and ClimaLand's `ρe_int`, and it is adopted wholesale. CLM5's temperature-prognostic-plus-two-step-phase-adjust is **rejected** in favor of the enthalpy state (Bundle D rec #1, #7): phase change becomes a state read-off, energy conserves by construction, and the branch-free algebraic inversion is GPU-friendly.

**CAS is the one enthalpy-not-energy store.** The CAS carries **specific enthalpy per kg of moist air `[J/kg]`** (ED2 `can_enthalpy`), companioned by specific humidity `can_shv`, with `can_temp` *diagnosed*. This choice fuses sensible and latent into one linear state function and makes the CAS heat capacity identically equal to its water-vapour capacity (§4b) — an elegance lost if the CAS carried temperature. There is no soil-air phase change, so temperature *alone* would be defensible, but enthalpy mirrors ED2 exactly and stays robust to condensation and to the sensible/latent partition (Bundle B).

### 1.2 Target processes and MVP vs deferred

| Process | Governing choice | MVP (P1) | Deferred |
|---|---|---|---|
| **Leaf energy balance** | prognostic `leaf_energy`, `Rn−H−LE`, linearized BE step (§4a) | ✔ liquid-only (`fliq≡1`) | leaf-water freeze/thaw (P2a) |
| **Wood energy balance** | prognostic `wood_energy`, `Rn−H` (no transpiration), `π·WAI` area (§4a) | ✔ | sap-flow enthalpy advection (needs dynamic hydraulics, P2b) |
| **CAS enthalpy budget** | prognostic `can_enthalpy`+`can_shv`, vapour-carried enthalpy twins (§4b) | ✔ single-layer, fixed `can_depth` | multi-layer `K_diff(z)`; variable-depth work term |
| **Soil heat diffusion** | implicit BE-Thomas on `soil_temp`, conservative `soil_energy` update (§4c, §6) | ✔ conduction + advective heat, liquid-only | freeze/thaw plateau (P2a, automatic in inverter) |
| **Ground/surface skin** | surface energy balance as **top-layer flux BC**; `t_ground = soil_temp(1)` (§4d) | ✔ | explicit thin skin layer for diurnal-T fidelity (P2b) |
| **Soil thermal conductivity** | de Vries / Johansen `κ(θ, fliq, texture)` (§4c, §8) | ✔ (ice-aware form; `fliq≡1` in P1) | organic/quartz refinement |
| **Coupled surface fixed point** | leaf↔CAS↔ground↔soil simultaneous solve | **DEFERRED (P3)** | the whole coupling |

**Thermal-only in v1, and phase-change-off in P1.** Exactly as the hydrology module shipped water-only with `fliq≡1` and freezing a hook, this family ships **liquid-only in P1** (`fliq≡1` everywhere; the freezing plateau is exercised only by the inverter's unit test) and turns on the freeze/thaw plateau at **P2a**. Because the state is internal energy and the inverter reads `fliq` off `u`, **P2a freeze/thaw requires no change to the solver or the state layout** — only the inverter's plateau branch and the ice-aware conductivity/heat-capacity terms activate. Snow / temporary-surface-water thermal layers are deferred to P2b with the hydrology's sfcwater stack.

### 1.3 Reachability facts from the ED2 reference

- ED2 **does not integrate temperature** for any store — leaf, wood, and every soil layer carry internal energy in the RK4 tracer vector; temperature and liquid fraction are recovered each sub-step by `uextcm2tl`/`uint2tl` (Bundle A §1, Bundle C §1). MEDS mirrors this.
- ED2 **has no separate skin temperature**: `ground_temp = topsoil_temp` (`ed_therm_lib.f90:677`); the surface energy balance is imposed as the top-face flux of the top soil layer (`h_flux_g(mzg+1)`, Bundle C §4). MEDS adopts this; a distinct skin layer is an optional P2b fidelity upgrade, not the MVP.
- ED2's CAS thermodynamic state is a **single specific enthalpy** `can_enthalpy [J/kg]`; there is **no `H + λE` formulation** anywhere — every water flux carries its enthalpy twin and the atmosphere exchange is one `e_flux_ac`, with the pure sensible `h_flux_ac` reconstructed only as a diagnostic (Bundle B §2, §4). MEDS ports the enthalpy-flux form directly.
- ED2's bottom soil heat BC is **zero conductive flux** (adiabatic, no geothermal): the conduction loop starts at `klsl+1`, leaving `h_flux_g(klsl)` at its per-step-zeroed value; heat leaves the bottom **only advectively** with drained water (Bundle C §5). MEDS makes this BC **explicit** (`κ_face(N)=0` or a config geothermal value), not implicit in the loop bound.
- ED2 folds `ρ_air·c_p` into `leaf_gbh`/`wood_gbh` and `ρ_air·c_p·ggnet` into the ground exchange (Bundle A bug #1, Bundle B/C). MEDS keeps aerodynamic conductances **pure** (`[m/s]`) and multiplies `ρ_air·c_p` explicitly at the flux site — a deliberate divergence to avoid double-counting when porting.

---

## 2. Where it lives (library DAG)

### 2.1 The state/process wall

The DAG is `shared ← {allometry, plant} ← state ← demography ← aux ← main`, and `biophysics` is a **sibling of `plant`** that links `meds_shared` **only** (`target_link_libraries(meds_biophysics PUBLIC meds_shared)`). So the energy kernels **cannot name `site_t`**, **cannot `use meds_demography_types`**, and **cannot read a global config** — the compute is stateless, the prognostic thermal state is owned by `meds_state`, and the two are woven only in `meds_aux`. This is the identical split the RT, hydraulics, and column-hydrology kernels already use.

**The controlling precedent: leaf/wood energy is to `cohort_block` what `psi_node` is to `cohort_block`; soil energy is to `patch_index` what `soil_water` will be to `patch_index`; CAS enthalpy is to `patch_index` what `surface_water` will be to `patch_index`** — prognostic state in `meds_state`, stateless kernels in `meds_biophysics` receiving it as `intent(inout)` value types, woven together in `meds_aux`.

**Typed-vs-plain ownership (the hydrology blocker, re-applied).** `meds_config` is the DAG root in `src/shared`; it **must not** `use meds_biophysics_types`. Consequently the thermal parameter type `soil_thermal_params_t` is a **biophysics** type and **never a `meds_config_t` component**; `meds_config_t` carries only plain `real`/`integer` texture/thermal arrays (`soil_quartz`, `soil_dry_conductivity`, `soil_dry_heat_capacity`, …). The typed `soil_thermal_params_t` is assembled by a biophysics routine `meds_soil_thermal%build_soil_thermal(cfg) → soil_thermal_params_t`, called from the aux/init layer — the exact analogue of **`build_soil_params`** (`meds_soil_parameters.f90:133`, the routine that assembles the typed `soil_params_t`) and `derive_rad_optics` (`meds_optics.f90:206`).

**The `ENERGY_*` selector enums live in `meds_biophysics_types`, next to the existing `SOIL_*` codes** (`meds_biophysics_types.f90:30–116` — `SOIL_SOLVER_BE`, `SOIL_BC_*`, `SOIL_SUBSTEP_*`, …), following the current self-contained convention (no new `use meds_config` edge). If and when the whole family of selector codes migrates to `meds_config` at P3, `SOIL_*` and `ENERGY_*` migrate *together* and `meds_biophysics_types` would then `use meds_config` for the enum params; the MVP does not do this. TOML strings are mapped to the codes by `meds_config_io`.

**Soil geometry/state is P3-deferred, not existing infrastructure.** Per CLAUDE.md the whole hydrology config+state — soil geometry (`soil_layer_z`/`dz`/`nzg`), `patch_index%soil_water`, and the `psi_soil` coupling — is itself **P3-deferred**; `meds_config` has no soil geometry today and `build_soil_params` is called only from `test_column_hydrology.f90`. The thermal arrays are therefore derived in `derive_config` **as the P3 end-state that lands together with the hydrology config/state**, not as something already present. The only blend precedent that exists *today* to point at is `recruit_pool` in `fuse_2_patches` (`meds_demography_fusefiss.f90:435`); `soil_water` is the sibling blend that lands with the same P3 work.

```
shared ─┬─ allometry ─ state ─ demography ─┐
        │  (meds_config owns, AT P3, with  │
        │   hydrology config: soil geom +  ├─ aux (fast-loop driver: RT → leaf-gas-exchange
        │   texture+thermal as PLAIN real  │      → leaf/wood energy → CAS → hydrology → soil heat;
        │   arrays)                        │      the ⊥ weave, host-only; P3 coupled fixed point;
        ├─ plant  (leaf gas exchange,      │      also calls build_soil_thermal → soil_thermal_params_t)
        │          respiration, phenology, │
        │          hydraulics — CONSUMES   │
        │          the closed temperatures)│
        └─ biophysics(core: shared-only) ──┘
             meds_column_energy    ← ALL stateless energy kernels for the whole soil-veg-air
                                      column: soil-heat BE-Thomas + leaf/wood/ground/CAS surface
                                      (sibling of meds_column_hydrology)
             meds_soil_thermal     ← κ(θ,fliq,texture), C_vol(θ,fliq) constitutive + build_soil_thermal
             meds_thermo (extend)  ← enthalpy↔T inverter, Clausius slope, air_density, cp_moist
             meds_soil_solver      ← REUSED verbatim (thomas_solve, nvfortran-safe subroutine)
             meds_biophysics_types ← thermal_state / energy_forcing / energy_flux types; ENERGY_* codes
                                     (next to the existing SOIL_* codes)
```

### 2.2 Files & CMake wiring

| File | Role | Analogue |
|---|---|---|
| `src/biophysics/meds_biophysics_types.f90` (extend) | `soil_energy_column_t`, `cas_state_t`, `leaf_energy_env_t`, `leaf_energy_flux_t`, `energy_forcing_t`, `energy_flux_t`, `soil_thermal_params_t`, `veg_thermal_params_t`, `cas_atm_forcing_t`, `energy_opts_t`; `ENERGY_*` codes | the `soil_column_t`/`chydro_*`/`SOIL_*` block |
| `src/biophysics/meds_soil_thermal.f90` (new) | `pure`/`elemental` `soil_thermal_cond(θ,fliq,…)`, `soil_heat_cap_vol(θ,fliq,…)`; `build_soil_thermal(cfg)` | `meds_soil_parameters` (`build_soil_params`) |
| `src/biophysics/meds_column_energy.f90` (new) | **ALL energy kernels for the whole soil-veg-air column.** Soil seam `soil_energy_flux` + inner device-eligible `soil_heat_be_step` (conservative `soil_energy`/`soil_temp` update); leaf/wood cohort `veg_energy_balance`; per-patch `ground_surface_balance`; per-patch `canopy_air_update` | `meds_column_hydrology` |
| `src/shared/meds_thermo.f90` (extend) | enthalpy↔(T,fliq) inverter `uext_to_temp`, `d_sat_vapor_pressure_dt`, `air_density`, `cp_moist` | its own sat-vapour block |
| `src/shared/meds_constants.f90` (extend) | `cp_air`, `cp_vap`, `cp_liq`, `cp_ice`, `latent_heat_fusion`, `k_water`, `k_ice`, `k_air`, `t_3ple`, `tsupercool_liq`, `tsupercool_vap` | its `stefan`/`latent_heat_vap` block |
| `src/biophysics/meds_soil_solver.f90` | **REUSED unchanged** | — |
| `test/test_column_energy.f90` (new) | CTest, links `meds_biophysics` + `meds_testsupport` | `test_column_hydrology` |
| `test/test_surface_energy.f90` (new) | CTest, leaf/wood/CAS/ground kernels | `test_canopy_radiation` |

CMake globs `src/biophysics/*.f90` into `libmeds_biophysics` with `CONFIGURE_DEPENDS`, so the new modules auto-add. Append `test_column_energy`/`test_surface_energy` executables following the `test_column_hydrology` block. **Build nvfortran multicore on every new module** — a green ifx suite is not sufficient (§10).

### 2.3 The fast-loop gap and standalone testability

There is still no sub-daily loop and no meteorological forcing (`advance_one_step` calls only the slow loop; `STUB_TISSUE_TEMP = 298.15`). So — exactly like RT, leaf, hydraulics, and column-hydrology — this family **ships standalone on its unit tests**: every kernel takes forcing as a value type (`energy_forcing_t`, `leaf_energy_env_t`, `cas_atm_forcing_t`), never from a global. The still-absent piece is the met reader (`src/io`: SW/LW down, air T, VPD, wind, precip), not these kernels.

---

## 3. State variable & representation

### 3.1 The shared inverter (the `uextcm2tl` analogue)

Every phase-changing store recovers `(temperature, liquid_fraction)` from prognostic internal energy through **one shared `elemental` inverter in `meds_thermo`**, `uext_to_temp(uext, wmass, dry_hcap, temp, fliq)` — the direct port of ED2 `uextcm2tl`/`uint2tl` (Bundle A §1, Bundle C §1). Given extensive internal energy `uext`, total water mass on the store `wmass`, and the dry (water-free) heat capacity `dry_hcap`:

```
u_freeze = (dry_hcap + wmass·cp_ice)·t_3ple                          ! energy if all-ice at triple point
u_melt   = u_freeze + wmass·latent_heat_fusion                       ! energy if all-liquid at triple point

uext <  u_freeze :  fliq = 0 ;  temp = uext / (dry_hcap + wmass·cp_ice)                             ! all ice
uext >  u_melt   :  fliq = 1 ;  temp = (uext + wmass·cp_liq·tsupercool_liq)/(dry_hcap + wmass·cp_liq) ! all liquid
u_freeze==u_melt :  fliq = 0.5 ; temp = t_3ple                       ! negligible-water singularity guard
otherwise        :  temp = t_3ple ; fliq = clamp((uext−u_freeze)/(wmass·latent_heat_fusion),0,1)     ! plateau
```

with `t_3ple = 273.16 K`, `cp_ice ≈ 2093`, `cp_liq ≈ 4186 J/kg/K`, `latent_heat_fusion (L_f) ≈ 3.34e5 J/kg`, and the **ED2-consistent liquid reference**

```
tsupercool_liq = t_3ple − (cp_ice·t_3ple + latent_heat_fusion)/cp_liq ≈ 56.8 K.
```

This is the value implied by ED2's `cmtl2uext`/`uextcm2tl`, satisfying `cp_liq·(t_3ple − tsupercool_liq) − cp_ice·t_3ple = L_f`. It is chosen precisely so the **liquid branch returns `t_3ple` at `uext = u_melt` for any `dry_hcap`** — the ice branch is referenced at 0 K (`temp = uext/(dry_hcap+wmass·cp_ice)`), so the liquid branch must be offset by exactly this amount to make `T(u)` **continuous at both `u_freeze` and `u_melt`**. (The naïve `t_3ple − L_f/cp_liq ≈ 193.36 K` is *wrong* here: it would make the liquid branch jump to ≈409.7 K at the top of the plateau for `dry_hcap=0, wmass=1`, a >130 K discontinuity — it is only correct if the ice branch is itself re-referenced at `t_3ple`, which MEDS's 0-K ice branch is not.) A round-trip **continuity test** asserts `temp(u_freeze) = temp(u_melt) = t_3ple` to round-off (test 2).

The forward map `temp_to_uext(dry_hcap, wmass, temp, fliq)` (ED2 `cmtl2uext`) rebuilds `uext` whenever heat capacity or water mass changes, keeping temperature continuous across the phase discontinuity — this is how creation/fusion sites stamp a fresh cohort's `leaf_energy` from an initial temperature (§7).

The `u_freeze == u_melt` fallback (`fliq=0.5, temp=t_3ple`) is the deliberate **negligible-water singularity guard**; the plateau branch's `wmass·L_f` denominator would otherwise blow up as `wmass→0`. MEDS uses a banded `abs(u_melt − u_freeze) < eps` test (not exact-float equality) — Bundle B/C flag ED2's exact `== t3ple` comparisons as brittle.

**The `dry_hcap`/`wmass` split — one convention, consistently.** The dry heat capacity holds only the water-free tissue/soil matrix; the water contribution enters dynamically through `wmass·(cp_ice | cp_liq)` inside the inverter, **not** baked into `dry_hcap`. This is ED2's convention and it correctly tracks changing water mass and phase (Bundle A §3, Bundle C §3). The port must pick this convention consistently — ED2's static-hydraulics path folds bound water into the tissue specific heats instead, and mixing the two **double-counts internal-water heat capacity** (Bundle A bug #5). MEDS puts **all** store water in `wmass`.

### 3.2 Per-store heat capacities

```
leaf: dry_hcap = nplant·C2B·bleaf·c_leaf(pft)                                  [J/m²/K]   (Bundle A §3)
      wmass    = leaf_water + leaf_water_hydraulic                            [kg/m²]    interception film + xylem water
wood: dry_hcap = nplant·C2B·(bsap·c_sapw(pft) + bdead·c_dead(pft) + bbark·c_bark(pft))
      wmass    = wood_water + wood_water_hydraulic
soil: dry_hcap = soil_dry_heat_capacity(k)   [J/m³/K]  (texture param, §8)
      C_eff(k) = soil_dry_heat_capacity(k) + soil_water(k)·rho_h2o·(fliq·cp_liq + (1−fliq)·cp_ice)
CAS:  hcapcan  = air_density·can_depth   [kg/m²]  (≡ water capacity, §4b — specific-enthalpy convention)
```

`c_leaf/c_sapw/c_dead/c_bark` are per-tissue specific heats `[J/kg/K]` (PFT params, §8); `C2B ≈ 2` (carbon→dry-biomass); `bleaf/bsap/bdead/bbark` `[kgC/plant]` from allometry. Wood uses **sapwood + heartwood + bark**, not `c_leaf` — keep them distinct (Bundle A §4).

---

## 4. Governing equations per store

Sign convention: heat fluxes positive **out of** the store into the CAS / atmosphere for leaf/wood/ground; radiation absorbed is a positive source; all leaf/wood/ground/CAS turbulent fluxes are **per unit ground area** `[W/m²]`, extensive over the cohort's LAI/WAI or the patch. Soil fluxes are `[W/m²]` on the negative-z column, **positive upward** (out of the ground).

### 4a. Leaf / wood cohort energy balance

**Leaf tendency** (per cohort, ED2 `rk4_derivs.f90:1660`; prognostic `leaf_energy [J/m²]`):

```
d(leaf_energy)/dt = rshort_l + rlong_l − h_flux_lc − qw_flux_lc − q_transp + leaf_qintercepted    [W/m²]
```

- **Absorbed shortwave `rshort_l`** `[W/m²]` — from the canopy RT solver (`rad_flux_t%abs_leaf(VIS)+abs_leaf(NIR)`). Pure source; no local gradient computed here. `par_l` is carried separately for photosynthesis; the *energy* budget uses full `rshort_l`.
- **Net longwave `rlong_l`** `[W/m²]` — `rad_flux_t%abs_leaf(LW)`, which is **already NET of each cohort's own thermal emission**: the two-stream solver folds `emission(coh) = stefan·canopy_temp⁴` into the absorbed-flux scattering source (`meds_canopy_radiation.f90:76`, `meds_twostream_band.f90:71–73`), driven by the per-cohort emission temperature `canopy_temp` fed IN via seam #9 (`canopy_temp = flux%temp`, §9). **The leaf kernel therefore does NOT subtract `2·ε_leaf·σ·leaf_temp⁴·LAI` again** — doing so would double-count LW loss and make the leaf budget wrong. (The local emission response is retained *only* as a linearization term in the BE Jacobian for stiff-mode damping, §6.2 — it never appears in the conservative tendency `R`.)
- **Sensible heat to CAS** `h_flux_lc = effarea_heat·LAI·gbh·ρ_air·cp_air·(leaf_temp − can_temp)`, `effarea_heat = 2` (both leaf sides). **MEDS keeps `gbh` a pure aerodynamic conductance `[m/s]` and multiplies `ρ_air·cp_air` explicitly** — ED2 bakes them into `leaf_gbh`; do not double-count (Bundle A bug #1).
- **Interception-film evaporation** `qw_flux_lc = w_flux_lc · enthalpy_vapor(leaf_temp)`, with the wetted-fraction throttle `w_flux_lc = effarea_evap·LAI·gbw·ρ_air·(qsat(leaf_temp) − can_shv)·σ_w`, `effarea_evap = 1`, `σ_w = min(1,(leaf_water/w_max)^{2/3})` (the same Deardorff wetted fraction the hydrology module already exports). `enthalpy_vapor(T)` is the specific enthalpy of the vapour leaving `[J/kg]` (ED2 `tq2enthalpy`, §4b) — the energy carried away as latent heat.
- **Transpiration** `q_transp = transp · enthalpy_vapor(leaf_temp)`, with

  ```
  transp = ρ_air·LAI·(fs_open·g_open + (1−fs_open)·g_closed)·(qsat(leaf_temp) − can_shv)·effarea_transp   [kg/m²/s]
  ```

  where the conductance is boundary layer **in series with** stomata `g = gbw·gsw/(gbw+gsw)`, and `effarea_transp` is the **per-PFT** sidedness (2 amphistomatous, else 1 — multiplies transpiration only, not sensible; Bundle A bug #3). **`ρ_air` is mandatory** — `qsat`/`can_shv` are specific humidities `[kg/kg]` and `g` is `[m/s]`, so without `ρ_air [kg/m³]` the driver would carry the wrong units; this matches the sibling `w_flux_lc`, which also multiplies `ρ_air`. Transpired water is always liquid (`fliq=1` in the enthalpy) and its *mass* comes from the soil/xylem, not `leaf_water`.
- **Interception source `leaf_qintercepted`** — enthalpy of precipitation caught on the leaf (mass twin owned by the hydrology `intercept_canopy_layer`). Shedding is handled in the slow adjust step, **not** the tendency (ED2 `wshed≡0` in-derivative; Bundle A bug #6).

**Wood tendency** (ED2 `rk4_derivs.f90:1942`) = leaf budget **minus transpiration**, on the **`π·WAI`** area basis (branches as cylinders — circumference, not `2·LAI` flat plates; Bundle A bug #2):

```
d(wood_energy)/dt = rshort_w + rlong_w − h_flux_wc − qw_flux_wc + wood_qintercepted             [W/m²]
h_flux_wc = pi·WAI·gbh_wood·ρ_air·cp_air·(wood_temp − can_temp)
```

Wood has no stomata (no `q_transp`, no soil extraction). Under **dynamic hydraulics (P2b)** an internal sap-flow enthalpy `qw_flux_wl = w_flux_wl·internal_energy_liquid(T_source)` moves with xylem water from wood to leaf (upwind on the source temperature, using **liquid internal energy `cp_liq·(T−tsupercool_liq)`, NOT vapour enthalpy** — no phase change on internal advection; Bundle A §4, bug #4); `d(leaf_energy) += qw_flux_wl`, `d(wood_energy) −= qw_flux_wl`.

**Sap-stream energy is a BIOPHYSICS term that consumes hydraulics WATER fluxes (ownership).** The whole transpiration path — soil layer `k` (`soil_temp(k)`) → root → wood (`wood_temp`) → leaf (`leaf_temp`), then evaporation — carries *sensible* liquid enthalpy up the column, and this energy accounting lives **here in the energy balance**, not in `plant/`. Plant hydraulics stays purely **water + potential**: it supplies the fluxes (per-layer `root_uptake`, the per-segment xylem/sap flow `w_flux_wl`); the energy family multiplies each by `internal_energy_liquid(T_source)` upwind and books it as a source/sink on the stores it connects (soil loses `uptake·u_liq(soil_temp(k))` in the §4c advective term; wood/leaf exchange `qw_flux_wl`). This is the exact mirror of the soil advective-heat term (§4c), which already consumes the hydrology `w_flux`. Conservation along the chain is by construction — every segment removes `u_liq(T_source)` from the upwind store and adds it to the downwind store; the vapour **latent** heat departs separately at the leaf as `q_transp`. (P1 runs with these advective terms off — `fliq≡1`, static hydraulics — and turns them on at P2b when dynamic hydraulics provides the segment fluxes.)

**Dew, frost & condensation are automatic — with two caveats.** Every latent flux above uses the gradient `(qsat(T_surface) − can_shv)`, so when the air is more humid than the surface saturates (`can_shv > qsat`) the flux is **negative = condensation**: negative `w_flux_lc`/`w_flux_gc` deposits water on the leaf/ground and its enthalpy twin `qw_flux < 0` **releases latent heat** into the store — dew forms with no special branch. Because phase is read off the store's internal energy, deposition when the surface is below freezing is **frost/deposition** automatically (the added mass takes the surface `fliq`) — a free benefit of the enthalpy state. Two conditions make it work end-to-end: **(i)** the leaf wetted-fraction `σ_w` must throttle **only evaporation**, never condensation, or dew can never nucleate on a dry leaf (`σ_w=0`); apply it as `σ_eff = (w_flux_lc ≥ 0) ? σ_w : 1` (ED2 keeps a separate dew branch). **(ii)** the hydrology's ground-evaporation `E_soil ≥ 0` clamp ("no dew in v1") must be **lifted** so `w_flux_gc` can go negative — a coupling-time change, since the condensed mass must land in the soil/pond water store (the mass twin of the released latent heat).

**Numerics: prognostic, L-stable linearized BE step (§6.2).** Recover `(leaf_temp, fliq)` from `leaf_energy` via the inverter; linearize the tendency about `leaf_temp^n`; take one implicit step for the tiny-heat-capacity stiff mode; then advance `leaf_energy` conservatively and re-invert. This is ED2's transient leaf (prognostic `leaf_energy`) made **unconditionally stable** — CLM's zero-heat-capacity steady-state Newton is rejected because ED2 in fact prognoses leaf energy (Bundle A §1 is authoritative; Bundle D mislabels ED2 as steady-state), and the implicit step removes the stiff-mode objection Bundle D raised against a prognostic leaf.

**Resolvability floor (Bundle A §5).** A cohort's leaf/wood energy is integrated only if resolvable: `leaf_resolvable = exposed .and. (leaf_hcap > veg_hcap_min(pft))`, where `veg_hcap_min` is a per-PFT heat-capacity floor (from a minimum resolvable plant density at the smallest leaf size). Below it, the cohort is slaved to `can_temp` and its tendency is zeroed — this prevents the `1/leaf_hcap` temperature inversion from exploding for near-zero-mass cohorts. The transition is booked **conservatively** by a `temp_to_uext` re-derivation (energy re-computed from retained temperature) so no energy is created/destroyed on the resolvable↔unresolvable switch.

### 4b. Canopy-air-space energy budget

**Prognostic specific enthalpy `can_enthalpy [J/kg]`** (Bundle B). The enthalpy⇄temperature relation is a **linearised, phase-free** state function (constant specific heats, linear-in-T latent heat collapse the phase discontinuity into one reference offset):

```
forward:  can_enthalpy = (1−shv)·cp_air·T + shv·cp_vap·(T − tsupercool_vap)                 [J/kg]
inverse:  can_temp = (can_enthalpy + shv·cp_vap·tsupercool_vap) / ((1−shv)·cp_air + shv·cp_vap)   [K]
```

with `cp_vap = 1859 J/kg/K` (water **vapour** cp — **not** `cp_liq=4186`; Bundle B bug #1) and `tsupercool_vap = t_3ple − (u_ice_t3 + latent_heat_vap_sub)·(1/cp_vap) ≈ −1.56e3 K`, the extrapolated reference at which vapour specific enthalpy is zero. This single offset bakes the latent heat of vaporisation/sublimation baseline into `cp_vap·(T − tsupercool_vap)` = the *total* (thermal + phase) vapour specific enthalpy — which is **why any vapour flux automatically transports its latent heat** and why no separate `λE` term appears anywhere in the CAS budget.

**Full tendency** (ED2 `rk4_derivs.f90:2156`; `hcapcani = 1/hcapcan`):

```
d(can_enthalpy)/dt = hcapcani·(  Σ_coh h_flux_lc + Σ_coh h_flux_wc + h_flux_gc          ! sensible: leaf, wood, ground
                               + Σ_coh qw_flux_lc + Σ_coh qw_flux_wc + Σ_coh q_transp   ! vapour-enthalpy twins: leaf/wood/transp
                               + qw_flux_gc − qdew_gnd_flux                           ! vapour-enthalpy: soil evap, dew
                               + e_flux_ac )                                       ! TOTAL enthalpy exchange atm↔CAS
```

Each sensible term is `g·A·ρ_air·cp_air·(T_source − can_temp)`; each **vapour-enthalpy twin** is `qw_flux = w_flux · enthalpy_vapor(T_source)` with the source temperature evaluated at **the surface the water leaves** (`leaf_temp`, `wood_temp`, `t_ground`), never `can_temp` (Bundle B §2, quirk). The **coupling rule for the port is one line: whenever you move water mass `m` at temperature `T`, move enthalpy `m·enthalpy_vapor(T) = m·cp_vap·(T − tsupercool_vap)` in the same step**, using the *same* `T_source` in both the mass driver and its enthalpy twin. **This rule binds the ground twin as tightly as the leaf twin** — the ground/soil surface loses `w_flux_gc·enthalpy_vapor(t_ground)` in `G_top` (§4d) and the CAS receives the identical `qw_flux_gc = w_flux_gc·enthalpy_vapor(t_ground)`, so the closed soil↔CAS budget balances to round-off.

**Atmosphere exchange is the whole story in `e_flux_ac`:** `e_flux_ac = ρ_air·ustar·estar`, `estar = c3·(enthalpy_atm − enthalpy_can)`, the Monin–Obukhov enthalpy scale. `e_flux_ac` already contains **both** the dry sensible flux and the full vapour enthalpy carried by `w_flux_ac`. The pure sensible `h_flux_ac = e_flux_ac + w_flux_ac·cp_vap·tsupercool_vap` is **diagnostic-only** — never used to drive the state. **Do not attempt an `H + λE` formulation** (Bundle B bug #2): port the single-`e_flux_ac` enthalpy-flux form or you will double-count / drop the latent baseline.

**Heat capacity ≡ water capacity** (Bundle B §3, `can_whccap8`): because the state is *specific* enthalpy per kg of air, the capacity converting a `[W/m²]` flux to a `[J/kg/s]` tendency is the **column air mass per ground area** `hcapcan = air_density·can_depth [kg/m²]` — identically the water-vapour capacity `wcapcan`. The same `ρ·Δz` divides the enthalpy sum and the `can_shv` water-vapour sum (`d(can_shv)/dt` is the mass twin, §9). MVP fixes `can_depth`; the variable-depth 1st-law work term (`ddens_dt_effect`) is a known small non-conservation flagged for P2.

**Numerics: prognostic single step, no root-find** (Bundle E C.3) — the CAS is the mixing node. `can_enthalpy^{n+1} = can_enthalpy^n + dt·(tendency)`; recover `can_temp`. Multi-layer extension is clean (§12): replace the single `ggnet`/`ustar` exchange with layer-to-layer `F_{k→k+1} = K_diff·ρ·(h_k − h_{k+1})/Δz`, per-layer `hcapcan(k)`; the same `K_diff` drives water, enthalpy, and CO₂ layer exchanges so the twin structure carries over unchanged.

### 4c. Soil thermal column

**Prognostic per-layer volumetric internal energy `soil_energy(k) [J/m³]`** on the negative-z column (`k=1` top/surface, `k=n_soil_layer_active` bottom; geometry reused verbatim from `soil_params_t`). Recover `(soil_temp(k), soil_fliq(k))` via the inverter with `wmass = soil_water(k)·rho_h2o` and `dry_hcap = soil_dry_heat_capacity(k)`.

**Volumetric heat diffusion** (ED2 `rk4_derivs.f90`, CLM5 §2.6). With faces indexed so that face `k` sits between layer `k` and layer `k+1` below it, and all fluxes **positive up**:

```
d(soil_energy(k))/dt = dz_i(k)·[ (h_flux(k) − h_flux(k−1))            ! conduction divergence (into k from below − out to above)
                               + (qw_flux(k) − qw_flux(k−1)) ] + Q_k  ! advective heat + root/absorbed heat sink
h_flux(k)  = −κ_face(k)·(soil_temp(k) − soil_temp(k+1))/dz_node(k)     [W/m²]   (positive up)
Q_k        = root_heat_sink(k)/dz(k)                                   [W/m³]   (units convention, §6.4)
```

`h_flux(0) = −G_top` is the top surface face (upward flux = negative of the net-into-ground `G_top`), so the layer-1 divergence `h_flux(1) − h_flux(0) = h_flux(1) + G_top` correctly **adds** `+G_top` to layer 1 — mutually consistent with the matrix RHS `r(1)` (§6.1). This is the correct sign: a hot surface layer warms the layer below it, and a positive `G_top` heats layer 1.

**Thermal conductivity vs moisture** — the de Vries / Johansen `κ(θ, fliq, texture)` (Bundle C §2a, Bundle D §3):

```
κ(θ, fliq) = K_e·κ_sat + (1 − K_e)·κ_dry                                          [W/m/K]
K_e        = max(0, min(1, log10(max(S_r, S_r_floor)) + 1))  (unfrozen)  |  max(0, min(1, S_r))  (frozen)
S_r        = θ/θ_sat ,   S_r_floor ≈ 0.05
κ_sat      = κ_solid^{1−θ_sat}·k_water^{θ_sat·fliq}·k_ice^{θ_sat·(1−fliq)}         geometric mean over phases
κ_dry, κ_solid  ← soil_dry_conductivity(k), soil_solid_conductivity(k)  (texture params, §8)
```

The **Kersten number `K_e` is clamped to `[0,1]`** with a saturation floor on the `log10` argument: Johansen's `K_e = log10(S_r)+1` is only valid above `S_r ≈ 0.05–0.1` and is unbounded below (as the column dries `K_e → −∞`, giving `κ < κ_dry` or even negative — non-physical and destabilizing — and `log10(0)` traps under `-fpe0`). The `max(S_r, S_r_floor)` and outer `max(0,min(1,·))` make the dry end well-behaved and FPE-safe (§10).

**Interface conductivity is the dz-weighted HARMONIC mean** (thermal resistors in series; Bundle D §3, adopting CLM over ED2) — ED2's geometric/log-linear interpolation overestimates the effective conductivity when adjacent layers differ strongly and lacks the series-resistor justification (Bundle C bug #3):

```
κ_face(k) = (dz(k) + dz(k+1)) / ( dz(k)/κ(k) + dz(k+1)/κ(k+1) )                    [W/m/K]
```

The dz-weighting is the correct series-resistor face conductance for the **non-uniform** negative-z grid (thin top layers); the unweighted `2·κ₁κ₂/(κ₁+κ₂)` equal-thickness form is dropped because it biases exactly the surface conduction (skin fidelity) the design cares about most.

**Advective heat** carried by the hydrology's inter-layer water flux `w_flux(k) [m/s]` (positive up), **upwind** on the source-layer temperature, liquid assumed (Bundle C §2c): `qw_flux(k) = w_flux(k)·rho_h2o·internal_energy_liquid(T_source)`, `T_source = soil_temp(k)` if `w_flux(k) ≤ 0` (downward) else `soil_temp(k+1)`. This is the ClimaLand `ρe_int,l·K∇h` term — energy travels with infiltrating/draining water, which CLM's temperature-only diffusion omits — and it is **nearly free once internal energy is the state** (Bundle D rec #6). `w_flux` is a forcing from the column-hydrology module (`chydro_flux_t` inter-layer fluxes exposed at the aux seam). Advection is **lagged/explicit** (not in the BE tridiagonal, §6.1), so its own first-order-upwind CFL bounds the substep independently of BE's L-stability.

**Freeze/thaw is automatic in the inverter** — no explicit latent term in the tendency (Bundle C §5, ClimaLand). Energy simply moves across `u_freeze`/`u_melt` and `soil_fliq` responds; the P1 build runs `fliq≡1` and P2a turns the plateau on with zero solver change. **Conductivity is ice-aware from the start** (`κ_sat` uses `fliq`), fixing ED2's bug of using liquid-water conductivity regardless of phase (Bundle C bug #2).

### 4d. Ground / surface skin balance

**MEDS adopts ED2's design: no separate skin temperature.** The surface energy balance is imposed as the **top-face flux BC on the top soil layer**, whose temperature `soil_temp(1)` *is* the prognostic skin; `t_ground = soil_temp(1)` (Bundle C §4). Closure is automatic once `soil_energy(1)` is advanced.

```
G_top = Rn_ground − H_ground − LE_ground                                          [W/m²]  (net into surface)
Rn_ground = rshort_g + rlong_g                        ! absorbed SW+LW at ground (rad_flux_t%dn_ground/up_ground)
H_ground  = (1−snowfac)·ggnet·ρ_air·cp_air·(t_ground − can_temp)                   ! sensible to CAS
LE_ground = w_flux_gc · enthalpy_vapor(t_ground)         ! soil evaporation — SAME enthalpy twin as qw_flux_gc (CAS side)
```

**`LE_ground` uses the identical enthalpy twin `w_flux_gc·enthalpy_vapor(t_ground)` that the CAS receives as `qw_flux_gc`** (§4b) — same water mass, same `T_source`, same specific enthalpy — so the closed leaf+CAS+ground+soil budget balances to round-off (test 1). The soil-water pool returns `internal_energy_liquid(t_ground)` with the removed mass. Using an ad-hoc `w_flux_gc·latent_heat_vap` (a *constant* specific energy) on the surface while the CAS receives `w_flux_gc·enthalpy_vapor(t_ground)` would leak ≈`w_flux_gc·(enthalpy_vapor − L_v) ≈ w_flux_gc·0.9e6 W/m²` and violate the E4/E6/E8 rule — MEDS does not.

`G_top` enters the soil BE as the top-layer Neumann flux (`h_flux(0) = −G_top`); the top layer integrates `Rn − H − LE − G_conduction` as a storage change, **not** an algebraic skin solve (Bundle C §4). This diagnosed `t_ground = soil_temp(1)` is exactly what fills **seam #2** (`chydro_forcing_t%t_ground`, unlocking energy-limited soil evaporation and evaporative cooling in the hydrology) and **seam #8** (`surface_state_t%soil_temp`, the two-stream longwave *lower* BC `grnd_emiss = soil_emiss·stefan·soil_temp⁴`). A distinct thin skin layer for higher diurnal-surface-T fidelity is a P2b upgrade (Bundle C bug #7) — the MVP keeps `dz(1)` small instead.

**Bottom BC = explicit zero-flux geothermal** (`κ_face(n_active) = 0`), or a config `geothermal_flux`/prescribed deep temperature. MEDS **asserts** this in the derivative rather than relying on the loop bound, fixing ED2's stale-flux hazard (Bundle C bug #1). Heat still leaves the bottom advectively with drained water (`qw_flux(n_active)`).

**Standing water & litter surfaces (deferred, P2b).** The surface balance above lands `G_top` directly on `soil_energy(1)` — the top BC is the **mineral-soil** surface. A persistent **standing-water pool** (a flooded forest / wetland — the hydrology `w_surface` ponding grown thermally significant) or a **litter / organic horizon** is physically a thin surface energy store *between* the atmosphere/CAS and `soil_energy(1)`: its own heat capacity (water `ρ_w·cp_liq·depth`, or litter organic `C_vol`), its own `(temp, fliq)` inversion, its own radiation / sensible / latent exchange with the CAS, and it conducts `G_top` down into the soil. This is the **same machinery** as the deferred snow / temporary-surface-water (sfcwater) thermal stack (§1.2, P2b): snow, standing water, and litter are all **stacked surface layers** sharing the inverter, the BE-Thomas conduction, and the surface energy balance — a small per-patch `surface_layer(:)` array, the analogue of ED2's `sfcwater_energy(nzs)`. **P1 lumps them into the soil top-layer BC** (no separate thermal mass — a flooded surface is treated as wet mineral soil); flooded-forest / wetland / litter fidelity, and the associated surface-water evaporative buffering, arrive with that stack. Wetland-specific effects (anoxia, methane, a near-constant water-table skin temperature) then follow naturally from a resolved standing-water layer.

---

## 5. Public seams & data types

Declared in `meds_biophysics_types` (its docstring reserves it as "the intended home for future energy-balance types"), next to the existing `SOIL_*` codes. Renamed away from any plant-type collision so the aux fast loop can `use` both. All **pure DATA** value types; each component carries its own default initializer (clean under `-fpe0`).

```fortran
integer(ik), parameter :: n_soil_layer_max = 20_ik            ! shared with hydrology
! ENERGY_* selector codes live HERE, beside SOIL_SOLVER_BE / SOIL_BC_* / SOIL_SUBSTEP_*:
integer(ik), parameter :: ENERGY_SOLVER_BE = 1_ik
integer(ik), parameter :: ENERGY_BC_GEOTHERMAL = 1_ik, ENERGY_BC_PRESCRIBED_T = 2_ik
integer(ik), parameter :: ENERGY_PHASE_OFF = 0_ik, ENERGY_PHASE_ON = 1_ik
integer(ik), parameter :: ENERGY_SUBSTEP_ADAPTIVE = 1_ik, ENERGY_SUBSTEP_FIXED = 2_ik

type :: soil_energy_column_t                                  ! MUTABLE prognostic slice (one patch)
   real(wp) :: soil_energy(n_soil_layer_max) = 0.0_wp         ! [J/m3] volumetric internal energy (PROGNOSTIC)
   real(wp) :: soil_temp(n_soil_layer_max)   = 0.0_wp         ! [K]    diagnosed each step
   real(wp) :: soil_fliq(n_soil_layer_max)   = 1.0_wp         ! [–]    diagnosed liquid fraction
end type

type :: cas_state_t                                           ! MUTABLE prognostic slice (one patch)
   real(wp) :: can_enthalpy = 0.0_wp                          ! [J/kg] specific enthalpy (PROGNOSTIC)
   real(wp) :: can_shv      = 0.0_wp                          ! [kg/kg] specific humidity (PROGNOSTIC, twin)
   real(wp) :: can_temp     = 0.0_wp                          ! [K]    diagnosed
   real(wp) :: can_depth    = 0.0_wp                          ! [m]    CAS depth (from canopy height; forcing)
end type

type :: soil_thermal_params_t                                 ! per-column THERMAL TEXTURE (biophysics type)
   integer(ik) :: nzg_active = n_soil_layer_max
   real(wp) :: soil_solid_conductivity(n_soil_layer_max) = 0.0_wp   ! [W/m/K] κ_solid
   real(wp) :: soil_dry_conductivity(n_soil_layer_max)   = 0.0_wp   ! [W/m/K] κ_dry
   real(wp) :: soil_dry_heat_capacity(n_soil_layer_max)  = 0.0_wp   ! [J/m3/K] slcpd (dry matrix)
   ! NOTE: porosity (theta_sat) and geometry come via soil_params_t (already passed) — NOT duplicated here.
end type

type :: veg_thermal_params_t                                  ! per-PFT vegetation THERMAL params (biophysics type)
   real(wp) :: leaf_emiss     = 0.95_wp                       ! [–]  LW emissivity (Jacobian −8εσT³ term, §6.2)
   real(wp) :: effarea_heat   = 2.0_wp                        ! [–]  sensible sidedness (both leaf sides)
   real(wp) :: effarea_evap   = 1.0_wp                        ! [–]  film-evap sidedness
   real(wp) :: effarea_transp = 1.0_wp                        ! [–]  transpiration sidedness (per-PFT; E3)
   real(wp) :: veg_hcap_min   = 20.0_wp                       ! [J/m2/K] resolvability floor (E17)
   real(wp) :: leaf_width     = 0.05_wp                       ! [m]  boundary-layer length scale
   real(wp) :: c_leaf = 0.0_wp, c_sapw = 0.0_wp               ! [J/kg/K] tissue specific heats
   real(wp) :: c_dead = 0.0_wp, c_bark = 0.0_wp
end type

type :: energy_forcing_t                                      ! soil-column thermal BCs (read-only)
   real(wp) :: g_top          = 0.0_wp                        ! [W/m2] net ground heat flux (Rn−H−LE), top Neumann
   real(wp) :: geothermal     = 0.0_wp                        ! [W/m2] bottom flux (default 0)
   real(wp) :: soil_water(n_soil_layer_max) = 0.0_wp          ! [m3/m3] θ from hydrology (κ, C_eff, wmass)
   real(wp) :: w_flux(n_soil_layer_max)     = 0.0_wp          ! [m/s]  inter-layer water flux (advective heat)
   real(wp) :: root_heat_sink(n_soil_layer_max) = 0.0_wp      ! [W/m2] enthalpy removed with root uptake
end type

type :: cas_atm_forcing_t                                     ! atmospheric forcing feeding the CAS (read-only)
   real(wp) :: ustar        = 0.0_wp                          ! [m/s]  friction velocity (M-O)
   real(wp) :: enthalpy_atm = 0.0_wp                          ! [J/kg] reference-level specific enthalpy
   real(wp) :: w_flux_ac       = 0.0_wp                          ! [kg/m2/s] atm↔CAS water-vapour mass flux
   real(wp) :: rho_air      = 0.0_wp                          ! [kg/m3] air density
end type

type :: energy_flux_t                                         ! soil outputs + diagnostics
   real(wp) :: ground_heat = 0.0_wp                           ! [W/m2] conductive flux into layer 1
   real(wp) :: bottom_heat = 0.0_wp                           ! [W/m2] advective+geothermal bottom loss
   real(wp) :: energy_resid = 0.0_wp                          ! [J/m2] closed-budget residual (~0)
   integer(ik) :: nsub = 0_ik
   logical  :: converged = .true.
end type

type :: leaf_energy_env_t                                     ! per-cohort surface BCs (leaf OR wood)
   real(wp) :: abs_sw = 0.0_wp, abs_lw = 0.0_wp               ! [W/m2] absorbed SW, NET LW (rad_flux_t; §4a)
   real(wp) :: can_temp = 0.0_wp, can_shv = 0.0_wp            ! [K],[kg/kg] CAS state (FORCED sibling)
   real(wp) :: gbh = 0.0_wp, gbw = 0.0_wp                     ! [m/s] boundary-layer heat/vapour conductance
   real(wp) :: gsw = 0.0_wp, fs_open = 1.0_wp                 ! [m/s] stomatal (leaf only), open fraction
   real(wp) :: area_index = 0.0_wp                            ! [m2/m2] LAI (leaf) or WAI (wood)
   real(wp) :: leaf_water = 0.0_wp, wmass = 0.0_wp            ! [kg/m2] film (σ_w) and total water (heat cap)
   real(wp) :: dry_hcap = 0.0_wp                              ! [J/m2/K] tissue heat capacity
   real(wp) :: rho_air = 0.0_wp, press = 0.0_wp               ! [kg/m3],[Pa] CAS air
end type

type :: leaf_energy_flux_t                                    ! per-cohort outputs
   real(wp) :: temp = 0.0_wp, fliq = 1.0_wp                   ! [K],[–] diagnosed store temperature
   real(wp) :: h_flux = 0.0_wp, qw_flux = 0.0_wp, q_transp = 0.0_wp ! [W/m2] sensible, film-evap, transp (→ CAS twins)
   real(wp) :: w_flux = 0.0_wp, transp = 0.0_wp                 ! [kg/m2/s] mass twins (→ CAS can_shv)
   real(wp) :: energy_resid = 0.0_wp                          ! [J/m2] closed-budget residual (~0; Test 1)
end type

type :: energy_opts_t                                         ! selectors + tolerances (codes in meds_biophysics_types)
   integer(ik) :: soil_solver   = ENERGY_SOLVER_BE
   integer(ik) :: bottom_bc     = ENERGY_BC_GEOTHERMAL        ! _GEOTHERMAL | _PRESCRIBED_T
   integer(ik) :: phase_change  = ENERGY_PHASE_OFF            ! _OFF (P1) | _ON (P2a)
   integer(ik) :: substep       = ENERGY_SUBSTEP_ADAPTIVE     ! | _FIXED (GPU)
   real(wp)    :: rtol = 1.0e-3_wp, atol = 1.0e-2_wp          ! atol in [K]
   real(wp)    :: h_init = 900.0_wp
   integer(ik) :: max_substep = 200_ik
   logical     :: debug_error = .false.                       ! error stop on energy_resid in Debug
end type
```

**The public soil seam** (host, derived-type interface; mirrors `column_hydrology_flux`):

```fortran
subroutine soil_energy_flux(col, forcing, therm, soil, opts, dt, flux)
   type(soil_energy_column_t), intent(inout) :: col      ! soil_energy(:) — SoA slice, only mutable thing
   type(energy_forcing_t),     intent(in)    :: forcing  ! G_top, θ, w_flux, root sink (read-only)
   type(soil_thermal_params_t),intent(in)    :: therm    ! κ texture (read-only)
   type(soil_params_t),        intent(in)    :: soil     ! geometry + porosity (theta_sat) — reused, not duplicated
   type(energy_opts_t),        intent(in)    :: opts
   real(wp),                   intent(in)    :: dt        ! [s]
   type(energy_flux_t),        intent(out)   :: flux
end subroutine
```

**Non-convergence contract.** On `nsub` reaching `opts%max_substep` the step sets `flux%converged = .false.`. Under `opts%debug_error` it logs and `error stop`s. **In production (`debug_error=.false.`) the contract is not silent:** `soil_energy_flux` returns without writing back a capped/unconverged column; the orchestrator (P3) must, on `converged=.false.`, **flag the patch and either clamp to the last converged sub-step state or propagate a hard error** — a non-converged energy step may never write state and proceed unchecked. The same contract binds the surface kernels and the CAS update below.

**The per-store surface kernels** (stateless, one call per cohort/patch):

```fortran
pure subroutine veg_energy_balance(store_energy, env, tparams, dt, is_leaf, flux)  ! leaf OR wood, per cohort
   real(wp),                intent(inout) :: store_energy    ! [J/m2] leaf_energy or wood_energy (PROGNOSTIC)
   type(leaf_energy_env_t), intent(in)    :: env             ! forced sibling temps + conductances
   type(veg_thermal_params_t), intent(in) :: tparams         ! per-PFT thermal (emiss, effarea_*, hcap_min, …)
   real(wp),                intent(in)    :: dt
   logical,                 intent(in)    :: is_leaf         ! toggles transpiration + area factor (2·LAI vs π·WAI)
   type(leaf_energy_flux_t),intent(out)   :: flux            ! temp, fluxes, mass twins, energy_resid
end subroutine

pure subroutine ground_surface_balance(t_ground, env, g_top, h_ground, le_ground)  ! per patch — INSTANTANEOUS
   ! returns t_ground = soil_temp(1) skin diagnostics + G_top for the soil seam.
   ! NO dt: ground has no prognostic variable (§1.1); it emits instantaneous fluxes only.
   real(wp),                intent(in)  :: t_ground          ! [K] = soil_temp(1) (from soil kernel)
   type(leaf_energy_env_t), intent(in)  :: env               ! can_temp, rho_air, conductances, Rn_ground
   real(wp),                intent(out) :: g_top, h_ground, le_ground

subroutine canopy_air_update(cas, coh_h_flux, coh_qw_flux, coh_w_flux, coh_transp,      &  ! per patch
                             ground_h_flux, ground_qw_flux, ground_w_flux, dew,          &
                             atm, dt, resid)
   type(cas_state_t),       intent(inout) :: cas             ! can_enthalpy/can_shv step + can_temp recover
   real(wp),                intent(in)    :: coh_h_flux, coh_qw_flux      ! [W/m2] Σ sensible, Σ vapour-enthalpy (leaf+wood+transp)
   real(wp),                intent(in)    :: coh_w_flux, coh_transp     ! [kg/m2/s] Σ MASS twins → can_shv budget
   real(wp),                intent(in)    :: ground_h_flux, ground_qw_flux! [W/m2]
   real(wp),                intent(in)    :: ground_w_flux, dew         ! [kg/m2/s] soil-evap + dew MASS twins
   type(cas_atm_forcing_t), intent(in)    :: atm             ! ustar, enthalpy_atm, w_flux_ac, rho_air
   real(wp),                intent(in)    :: dt
   real(wp),                intent(out)   :: resid           ! [J/m2] closed-budget residual (~0; Test 1)
end subroutine
```

The CAS update advances **both twins in one stateless call**: the enthalpy budget from the `h_flux`/`qw_flux` sums and `e_flux_ac(atm)`, and the co-prognostic `can_shv` from `d(can_shv)/dt = wcapcani·(coh_w_flux + coh_transp + ground_w_flux − dew + atm%w_flux_ac)` — so both the enthalpy and mass budgets close in isolation. The mass fluxes (`coh_w_flux`, `coh_transp`, `ground_w_flux`, `dew`, `atm%w_flux_ac`) are the missing water drivers the twin needs.

Two threads advancing two patches touch disjoint `col%soil_energy(:)` / `cas%…` ⇒ pure/reentrant. The soil interior arithmetic is delegated to a **bare-array device-eligible** inner routine `soil_heat_be_step(soil_temp, dz, dz_node, kappa, c_eff, g_top, geo, qw, q, nzg)` — the **`growth_step` precedent** (`meds_demography_dynamics.f90:52`, a genuine bare-array OpenMP-target routine). Note column-hydrology's inner step `soil_be_single_step` (`meds_column_hydrology.f90:276`) is **not yet** bare-array (it takes `soil_params_t`/`soil_opts_t`; device-eligibility deferred per CLAUDE.md), so `soil_heat_be_step` *establishes* the bare-array pattern rather than following the hydrology sibling; `energy_opts_t` stays out of it.

---

## 6. The solver / numerics

### 6.1 Soil heat — implicit backward-Euler Thomas (reuse `meds_soil_solver`)

The 1-D heat equation on the negative-z grid is structurally identical to the Richards solve. **Reject ED2's explicit adaptive RK4** (co-integrated with canopy energy, stiff, GPU-hostile — thin surface layers give a diffusion CFL `dt ≤ dz²/(2·max D)` collapsing to sub-seconds). Adopt **backward Euler (`α=1`), not Crank–Nicolson** (Bundle D rec #2): BE is monotone (no CN ringing at sharp freezing fronts), L-stable, and simpler to couple with the phase-change read-off; expose `α` as a parameter if CN is ever wanted.

Recover `soil_temp^n` from `soil_energy^n`; build the tridiagonal in `soil_temp^{n+1}` with dz-weighted harmonic-mean face `κ`, mirroring the hydrology rows exactly (`meds_column_hydrology.f90:311`):

```
a(k) = −κ_face(k−1)/dz_node(k−1)                     (sub-diagonal)
c(k) = −κ_face(k)/dz_node(k)                         (super-diagonal)
b(k) = C_eff(k)·dz(k)/h − a(k) − c(k)                (diagonal; C_eff = current effective vol. heat capacity)
r(k) = C_eff(k)·dz(k)/h·soil_temp(k)^n + dz(k)·Q_k   (source; TOP row r(1) ADDS +G_top, bottom the BC)
call thomas_solve(a, b, c, r, soil_temp_new, nzg)    ! nvfortran-safe subroutine, intent(out) result
```

with `Q_k = root_heat_sink(k)/dz(k) [W/m³]` (§6.4). The top row's `+G_top` in `r(1)` is consistent with the conservative update below, where `h_flux(0) = −G_top` makes the top face into layer 1 equal to `+G_top`.

**Conduction is implicit; advection is lagged/explicit (operator split).** The advective `qw_flux` is **excluded** from the tridiagonal (matrix carries conduction + `Q_k` only); it is applied explicitly in the energy update using upwind first-order fluxes. BE L-stability covers only the conductive stiff mode, so the substep must additionally respect the **advective CFL** of the explicit upwind term whenever `w_flux ≠ 0` (bounded in the adaptive controller; tests 3 and 5 include a nonzero-`w_flux` case to catch splitting error and CFL violations). If a fully consistent semi-implicit solve is later wanted, the upwind coefficients fold into `a/b/c` — deferred.

**Phase change lives entirely in the state inversion, NOT in the matrix** (Bundle E C.2, "more robust, ED2's choice"). `C_eff` uses the *current* phase-appropriate water heat capacity with **no latent spike** — the matrix stays well-conditioned. After the BE solve for `soil_temp^{n+1}`, advance the prognostic energy by the **conservative flux-divergence update** (the direct analogue of the hydrology θ-update, `meds_column_hydrology.f90:387`):

```
soil_energy(k)^{n+1} = soil_energy(k)^n
   + (dt/dz(k))·( h_flux(k) − h_flux(k−1) + qw_flux(k) − qw_flux(k−1) ) + dt·Q_k
```

with the conductive faces diagnosed from `soil_temp^{n+1}`, the advective faces from the lagged/upwind `w_flux`, `h_flux(0) = −G_top`, and `Q_k = root_heat_sink(k)/dz(k)`. Then invert `soil_energy^{n+1} → (soil_temp, soil_fliq)`. The sign `h_flux(k) − h_flux(k−1)` (not the reverse) is what makes a positive `G_top` **add** energy to layer 1, in agreement with the matrix RHS. Freeze/thaw is **automatic**: energy crossing `[u_freeze, u_melt]` pins `soil_temp = t_3ple` while `soil_fliq` runs 0→1, absorbing `wmass·L_f` of latent heat — no apparent-capacity spike, no CLM two-step "solve-then-adjust." This is decisively cleaner than CLM's temperature-based excess-energy bookkeeping and than ClimaLand's τ-relaxation (which introduces a non-physical thermal-time tuning parameter, Bundle D rec #7). The **exported `soil_temp` is re-derived from the stored energy**, never from the transport-variable iterate — always consistent with the conserved energy.

**Boundary conditions:** top = Neumann `G_top` from the ground surface balance (§4d); bottom = zero-flux geothermal (`κ_face(nzg)=0`) or config `geothermal`/prescribed deep T. Reuse the hydrology's **adaptive step-doubling / fixed-substep wrapper** (`ENERGY_SUBSTEP_ADAPTIVE`/`_FIXED`) so heat and water share substep control and the GPU warp-uniform path. BE is L-stable, so *conductive* substepping is **accuracy-only** (wetting/freezing front sharpness); the advective CFL adds a stability floor when `w_flux ≠ 0`.

### 6.2 Surface stores — small linearized single BE step

Leaf, wood, and ground are advanced by **one L-stable linearized step per sub-step** on the prognostic energy (quasi-steady heat capacity is small but nonzero — ED2/ClimaLand transient, not CLM's zero-capacity Newton). Recover `T^n`, linearize the tendency `R(T)` about `T^n`:

```
dR/dT = −8·ε·σ·T^3·LAI − effarea_heat·LAI·gbh·ρ·cp − L_v·ρ_air·g_total·(dqsat/dT)   (leaf; wood drops the g_total term)
T*   = T^n + R(T^n)/(C/h − dR/dT)                                                     (implicit increment)
store_energy^{n+1} = store_energy^n + dt·R_lin(T*)      ! conservative; then re-invert → (temp, fliq)
```

Two corrections vs a naïve form: (i) the **longwave damping is `−8·ε·σ·T³·LAI`** (the emission `2·ε·σ·T⁴·LAI` is two-sided, so its derivative is `8·ε·σ·T³·LAI`, not `4`) — this term is a *linearization for stiff-mode damping only* (the net LW itself comes from the RT via `abs_lw`, §4a; it is not re-added to the conservative tendency `R`); (ii) the **latent term carries `ρ_air`** (`−L_v·ρ_air·g_total·(dqsat/dT)`), matching the `ρ_air`-corrected transpiration driver (§4a). `dqsat/dT` is the **Clausius–Clapeyron slope** — the load-bearing nonlinear term — supplied by a new `meds_thermo` helper `d_sat_vapor_pressure_dt`, `dq_sat/dT = 0.622·P/(P − 0.378·e_sat)²·de_sat/dT`. Because photosynthesis needs `leaf_temp` and `leaf_temp` needs transpiration (which needs `gsw` from photosynthesis), the leaf balance and gas exchange are a small **inner** coupling — the natural first fixed point, lagging `gsw` from the previous sub-iterate. The CAS is advanced by a **prognostic single step, not a root-find** (§4b) — it is the mixing node.

### 6.3 The coupled fixed point — DEFERRED to P3

The user directive: coupling deferred, kernels stateless now. The design that makes that clean — **each per-store kernel takes the other stores' temperatures as forced inputs (the §9 seams)** — is honored by the signatures in §5: `veg_energy_balance` takes `env%can_temp`; `canopy_air_update` takes cohort/ground fluxes; `ground_surface_balance` takes `env%can_temp` and returns `G_top`; `soil_energy_flux` takes `forcing%g_top`. Until the coupling lands, each kernel runs with the sibling temperatures **lagged from the previous sub-step** — precisely today's placeholder contract. The **P3 orchestration** wraps them in the simultaneous leaf ↔ CAS ↔ ground ↔ soil solve (ED2's nested RK4 audit, or CLM's outer M-O / inner Newton, or a Picard fixed point over the sub-step); nothing in the kernels changes, only the driver stops lagging. This is the identical discipline the hydrology follows (`column_hydrology_flux` is "THE stateless seam," device-eligibility and coupling deferred).

### 6.4 Energy-conservation budget & nvfortran-safe solve

Flux form guarantees conservation **structurally**: summing the conservative update over `k`, every interior face telescopes, leaving only the two boundary faces and the sinks. With the **`Q`-units convention fixed once** — `Q_k = root_heat_sink(k)/dz(k)` is `[W/m³]` *inside* the kernel, while `root_heat_sink` and the residual are `[W/m²]` — the residual is

```
flux%energy_resid = Σ_k (soil_energy^{n+1}−soil_energy^n)·dz(k)
                    − dt·( G_top − bottom_heat − Σ_k root_heat_sink(k) )
```

so `dz(k)·Q_k = root_heat_sink(k)` and the layer-thickness factor cancels exactly (asserted in test 1). This ships as a diagnostic and `error stop`s in Debug (`opts%debug_error`) if `|resid| > atol` — the demography-invariant discipline. The **leaf/wood** (`leaf_energy_flux_t%energy_resid`) and **CAS** (`canopy_air_update`'s `resid` out-arg) paths carry the same instrumentation, so all four stores' round-off conservation checks (test 1) are testable through the public types, uniformly with the soil path.

`thomas_solve` is reused as a **subroutine with `intent(out) x(n_soil_layer_max)`**, never an array-valued function fed into a call (issue #7 — silently wrong at `-O2`, segfault at `-O0`); the surface kernels that return per-cohort/per-layer temperatures bind every array result to a named array first.

---

## 7. State additions & lockstep obligations

Geometry + texture stay **shared/immutable per site** (`soil_params_t` geometry + porosity, `soil_thermal_params_t` κ texture, `veg_thermal_params_t` per-PFT). The **prognostic thermal state lives where its process lives**: leaf/wood energy per COHORT, soil energy + CAS state per PATCH. All of this lands **at P3, together with the hydrology config/state** (soil geometry and `patch_index%soil_water` do not exist today).

**Per-cohort fields** on `cohort_block` (`meds_demography_types.f90`):

```fortran
real(wp), allocatable :: leaf_energy(:)      ! (cohort) [J/m2 ground]  PROGNOSTIC
real(wp), allocatable :: wood_energy(:)      ! (cohort) [J/m2 ground]  PROGNOSTIC
```

**Per-patch fields** on `patch_index`:

```fortran
real(wp), allocatable :: soil_energy(:,:)    ! (n_soil_layer_max, patch) [J/m3]  PROGNOSTIC — lands with soil_water
real(wp), allocatable :: can_enthalpy(:)     ! (patch) [J/kg]  PROGNOSTIC
real(wp), allocatable :: can_shv(:)          ! (patch) [kg/kg] PROGNOSTIC — CAS water twin
```

**Cohort lockstep** — `leaf_energy`/`wood_energy` ride the **single centralized cohort reorder** (`cohort_alloc` allocate+init, `site_free` dealloc, `cohort_ensure_capacity` copy, `move_alloc_block`, `cohort_reorder` permute, `copy_cohort_slot`; `cohort_compact`/`rebuild_csr` route through `cohort_reorder` — the CLAUDE.md "add a per-cohort field → update *these*" rule). **Not geometry-derived** ⇒ **not** in `set_cohort_size`. Because energy is extensive `[J/m²ground]`, **fusion sums** the two cohorts' `leaf_energy`/`wood_energy`; **fission splits** by daughter number/area share; conservation is total energy. Crucially, at every creation site (`add_cohort`/`init_bare_ground`, `apply_recruitment`, `split_cohorts`, `apply_patch_disturbance`) a fresh cohort's energy must be **stamped on-allometry from an initial temperature** via `temp_to_uext(dry_hcap, wmass, T_init, fliq_init)` — an energy-analogue call alongside `set_cohort_size`, exactly as `set_cohort_size` re-derives `agb`/`leaf_carbon`. A recruit inits at `T_init = can_temp` (thermal equilibrium with the CAS).

**Patch lockstep** — patches have no single reorder routine; carry `soil_energy(:,:)`, `can_enthalpy(:)`, `can_shv(:)` through the **six** sites: `patch_alloc` (init `soil_energy` from `T_init` per layer, `can_enthalpy` from `T_init`/`can_shv`), `patch_ensure_capacity`, `sort_patches`, `patch_compact`, **`fuse_2_patches`** (the **area-weighted blend** — combine `soil_energy`, `can_enthalpy`, `can_shv` by area fraction **exactly where `recruit_pool` is blended today**, `meds_demography_fusefiss.f90:435`; `soil_water` becomes the sibling blend when the hydrology state lands), `site_free`. **Disturbance-gap init copies the donor** (`apply_patch_disturbance` gap fragment copies `soil_energy(:,donor)`/`can_enthalpy(donor)` — the soil column and CAS are physically unchanged by canopy loss; resetting to `T_init` would inject spurious energy and break the budget). Only `init_bare_ground`/true cold start initializes from a prescribed profile.

The orchestration layer copies one patch's SoA slice → `soil_energy_column_t`/`cas_state_t` value → kernel (`intent(inout)`) → writes back into `patch%soil_energy(:,ip)`/`can_enthalpy(ip)`; the cohort sweep updates `cohort%leaf_energy(:)`/`wood_energy(:)` in place.

---

## 8. Parameters & config (no hard-coded values)

**No hard-coded model parameters** (CLAUDE.md): every texture-thermal, tissue-heat-capacity, and emissivity parameter is REQUIRED from TOML, presence-mapped, `error stop` if missing. Only genuine universal constants go in `meds_constants`.

**Constants to ADD to `meds_constants`** (true physical constants — allowed there, never TOML):

| constant | value | why |
|---|---|---|
| `cp_air` | 1004.6 J/kg/K | dry-air cp; sensible fluxes, CAS enthalpy |
| `cp_vap` | 1859 J/kg/K | water-**vapour** cp; ALL CAS enthalpy + vapour-enthalpy twins (NOT cp_liq) |
| `cp_liq` | 4186 J/kg/K | liquid-water cp; store water heat capacity, advected soil-water enthalpy |
| `cp_ice` | 2093 J/kg/K | ice cp; frozen-store heat capacity |
| `latent_heat_fusion` | 3.34e5 J/kg | freezing-band width; inverter plateau |
| `k_water` | 0.57 W/m/K | soil thermal conductivity constituent |
| `k_ice` | 2.29 W/m/K | frozen-soil conductivity (ED2 uses k_water regardless — bug #2) |
| `k_air` | 0.025 W/m/K | air-filled-porosity conductivity |
| `t_3ple` | 273.16 K | triple point (distinguished from `t_kelvin=273.15`) |
| `tsupercool_liq` | `t_3ple − (cp_ice·t_3ple + latent_heat_fusion)/cp_liq` ≈ 56.8 K | liquid internal-energy zero reference (ED2 `cmtl2uext`-consistent; §3.1) |
| `tsupercool_vap` | `t_3ple − (u_ice_t3+L_sub)/cp_vap` | vapour specific-enthalpy zero reference |

Present and reused: `stefan`, `latent_heat_vap = 2.501e6` (tagged "future energy balance"), `t_kelvin`, `rho_h2o`, `grav`, `r_wv`, `r_gas`, `p_std`.

**Thermo helpers to ADD to `meds_thermo`** (`pure`/`elemental`): `uext_to_temp`/`temp_to_uext` (the phase-change inverter/forward, §3.1); `d_sat_vapor_pressure_dt` (Clausius slope, §6.2); `air_density(t_k, p_pa, shv)` and `cp_moist(shv)` for the CAS terms; `enthalpy_vapor(t_k)` = `cp_vap·(t_k − tsupercool_vap)` and `internal_energy_liquid(t_k)` = `cp_liq·(t_k − tsupercool_liq)`. Keep the Bolton rule: universal empirical coefficients in the formula, tunable model parameters not.

**TOML parameters** (per-PFT tissue heat capacities and emissivity on `pft_table_t`; texture-thermal per soil class on the `[soil]` block):

```toml
[soil]                                   # thermal (extends the hydrology [soil] block)
soil_quartz          = 0.35              # quartz fraction (κ_solid geometric mean)   REQUIRED
soil_solid_cond      = 3.0               # [W/m/K] mineral solid conductivity          REQUIRED (or derived from quartz)
soil_dry_cond        = 0.15              # [W/m/K] dry-soil conductivity               REQUIRED
soil_dry_heat_cap    = 2.0e6             # [J/m3/K] dry-matrix volumetric heat cap     REQUIRED
soil_emissivity      = 0.96              # ground LW emissivity (two-stream lower BC)  REQUIRED
geothermal_flux      = 0.0               # [W/m2] bottom BC (default adiabatic)

[energy]
solver               = "backward_euler"  # ENERGY_SOLVER_BE (only supported)
bottom_bc            = "geothermal"      # geothermal | prescribed_t → ENERGY_BC_*
phase_change         = false             # P1 off (fliq≡1) | P2a on
substep_mode         = "adaptive"        # adaptive | fixed (GPU)
rtol                 = 1.0e-3
atol                 = 1.0e-2            # [K]
max_substep          = 200
debug_error          = false

[pft]                                    # per-PFT (pft_table_t → veg_thermal_params_t), one row each
c_leaf   = [3.2e3, 3.2e3, 3.2e3]         # [J/kg/K] leaf tissue specific heat
c_sapw   = [2.7e3, 2.7e3, 2.7e3]         # sapwood
c_dead   = [2.3e3, 2.3e3, 2.3e3]         # heartwood
c_bark   = [2.0e3, 2.0e3, 2.0e3]         # bark
leaf_emiss = [0.95, 0.95, 0.95]          # leaf LW emissivity
leaf_width = [0.05, 0.05, 0.05]          # [m] characteristic leaf dimension (boundary-layer gbh)
effarea_transp = [1.0, 1.0, 2.0]         # per-PFT transpiration sidedness (2 amphistomatous; E3)
veg_hcap_min = [20.0, 20.0, 20.0]        # [J/m2/K] resolvability floor (Bundle A §5)
```

Snow thermal parameters (density-dependent conductivity, Bundle D §3) are deferred to P2b with the sfcwater stack. Wiring mirrors the RT/hydrology precedent: scalar/thermal fields as plain arrays on `meds_config_t`, `ENERGY_*` enum `parameter`s in `meds_biophysics_types` (beside `SOIL_*`), mapped by `meds_config_io`, derived in `derive_config`, typed `soil_thermal_params_t` assembled by biophysics **`build_soil_thermal`** (modeled on `build_soil_params`). `validate_config`: `soil_dry_cond>0`, `soil_dry_heat_cap>0`, `0<soil_emissivity≤1`, `0<leaf_emiss≤1`, `c_*>0`.

---

## 9. Coupling contracts (DEFERRED wiring)

The headline: this family **closes every forced-temperature seam** MEDS currently hard-wires — **nine in all**. Exact fields, units, current forced value, and closer:

| # | Field (identifier) | Units | Location | Current forced value | Closed by |
|---|---|---|---|---|---|
| 1 | `leaf_env_t%leaf_temp` | K | `meds_plant_types.f90:53` | consumed as `t_leaf` (`meds_leaf_gas_exchange.f90:198`); drives all biochemistry T-scaling + `Kc/Ko/Γ*` | **leaf energy balance** → `flux%temp` |
| 2 | `chydro_forcing_t%t_ground` | K | `meds_biophysics_types.f90:130` | `= 298.15`, "FORCED = T_air until soil energy"; drives `α_soil`, `qsat(t_ground)` in `ground_evaporation` | **ground/top-soil skin** → `soil_temp(1)` |
| 3 | `wood_env_t%wood_temp` | K | `meds_plant_types.f90:246` | drives peaked-Arrhenius stem respiration | **wood energy balance** → `flux%temp` |
| 4 | `root_env_t%soil_temp` | K | `meds_plant_types.f90:271` | root-weighted soil T for fine-root respiration | **soil energy** → `root_frac`-weighted `soil_temp(:)` |
| 5 | `pheno_env_t%soil_temp` | K | `meds_plant_types.f90:186` | cold-drop phenology trigger (`cold_drop_soiltemp1/2`) | **soil energy** → shallow `soil_temp(k)` |
| 6 | `pheno_env_t%temp_day` | K | `meds_plant_types.f90:185` | daily-mean canopy/air T for GDD/chill sums | **CAS enthalpy** → `can_temp` (daily mean) |
| 7 | `carbon_env_t%tissue_temp` / `turnover_env_t%tissue_temp` | K | `meds_plant_types.f90:354` / `:302` | evergreen cold-suppression of turnover | **leaf/wood energy** → `flux%temp` |
| 8 | `surface_state_t%soil_temp` | K | `meds_optics.f90:50` | `= 298.0`; ground LW emission `grnd_emiss = soil_emiss·stefan·soil_temp⁴` (two-stream **lower** BC) | **ground/soil-surface** → `soil_temp(1)` |
| 9 | `canopy_radiation(…canopy_temp…)` | K | `meds_canopy_radiation.f90:35/40` | currently unwired (tests only); drives per-cohort thermal emission `emission(1:ncoh)=stefan·canopy_temp⁴` (`:76`) into `solve_band` (`meds_twostream_band.f90:71–73`) | **leaf/wood energy** → `flux%temp` |

**Seams #8 and #9 are the two halves of the physically important radiation↔energy circular coupling (the deferred fixed point).** `meds_optics`/`meds_canopy_radiation` need `soil_temp(1)` (#8) to emit longwave *down* from the ground **and** each cohort's `canopy_temp = flux%temp` (#9) to emit longwave from the leaf/wood — the latter is why `rad_flux_t%abs_leaf(LW)` is **already net** of a cohort's own emission (§4a). Symmetrically, the soil surface balance and the leaf/wood balance need the two-stream `up_ground`/`dn_ground`/`abs_leaf`/`abs_wood` (`rad_flux_t`) as their net-radiation *source*. In P1–P2 each kernel runs with the sibling lagged; the **P3 fast loop iterates radiation ↔ energy** (or lags one step): the leaf/wood/ground temperatures are fed INTO the RT (#8, #9), and the RT's net absorbed fluxes come back OUT as sources. The single orchestration seam that consumes all of this is the `TODO` at `meds_plant_interface.f90:184` (`call leaf_gas_exchange(env, cfg, ipft, flux) ! TODO: + hydraulics + leaf energy balance`); `get_plant_flux_fast` (`:179`) is where a diagnosed `leaf_temp` gets written into `leaf_env_t` before the gas-exchange call.

**The water↔energy CAS twin.** The CAS energy budget (`can_enthalpy`) and the CAS water budget (`can_shv`) are **structurally coupled** — every water flux in `d(can_shv)/dt = wcapcani·(Σ w_flux_lc + Σ w_flux_wc + Σ transp + w_flux_gc − dewgnd + w_flux_ac)` has an enthalpy twin `qw_flux = w_flux·enthalpy_vapor(T_source)` in the enthalpy budget, with `hcapcan ≡ wcapcan = air_density·can_depth` the shared divisor (Bundle B §3, §4). **The forced-`t_ground` hydrology soil-evaporation seam (#2) and this energy family are two halves of the same CAS control volume**: the hydrology owns `w_flux_gc` (soil evap mass, keyed on `t_ground`), this family owns `qw_flux_gc` (its enthalpy twin, keyed on the *same* `t_ground`, §4d) and diagnoses `t_ground = soil_temp(1)`. When both are wired at P3, the **evaporative-cooling feedback** closes — evaporation adds vapour enthalpy while sensible fluxes cool the air, `can_temp` drops via `hq2temp`, `can_rhv` rises, throttling further evaporation (Bundle B §5). Forcing `can_temp` (or evolving it without its `can_shv` partner) **severs this loop** and mis-fires the saturation guard — which is exactly why `can_enthalpy`+`can_shv` are prognostic and `can_temp` is derived, and why `canopy_air_update` takes the mass twins (§5).

**The fast-loop drive (reserved P3).** A fast-loop driver in `src/driver` (compiled into `meds_aux`), called from `advance_one_step` before `vegetation_dynamics`: per patch, per sub-step — `canopy_radiation` (now fed `canopy_temp` #9 and ground `soil_temp` #8) → `abs_leaf`/`up_ground` → `veg_energy_balance` (per cohort, lagged CAS) → `leaf_gas_exchange` (now with diagnosed `leaf_temp`) → `ground_surface_balance` → `G_top` → `soil_energy_flux` → `soil_temp(:)` → `canopy_air_update` (CAS enthalpy + `can_shv` from all cohort/ground mass+enthalpy fluxes) → the P3 Picard pass tightens the leaf/CAS/ground/soil/radiation fixed point. In this weave both `chydro_flux_t` (water) and `energy_flux_t` (heat) are in scope. The still-absent piece is the **met forcing reader** (SW/LW down, air T, VPD, wind, precip → `rad_forcing_t` + `energy_forcing_t` + `chydro_forcing_t` + `cas_atm_forcing_t`); until it lands the fast loop is disabled and this family stands on its unit tests.

---

## 10. GPU / nvfortran portability

- **Two-layer structure:** the seams `soil_energy_flux`/`veg_energy_balance`/`canopy_air_update` take derived types (host); the interior arithmetic `soil_heat_be_step(soil_temp, dz, dz_node, kappa, c_eff, g_top, geo, qw, q, nzg)` takes **bare fixed-size arrays + `firstprivate` scalar params**, no derived types (`energy_opts_t` stays out) — the **`growth_step` precedent** (`meds_demography_dynamics.f90:52`, a real bare-array OpenMP-target routine; the device cannot read host module vars, so texture/geometry flow in as arguments). Column-hydrology's `soil_be_single_step` is **not yet** bare-array (derived-type signature, device-eligibility deferred), so `soil_heat_be_step` establishes the pattern rather than following that sibling. **The parallel axis is columns (patches/sites), never within a column** — Thomas is a sequential recurrence, one sweep per thread. Per-patch/per-cohort orchestration stays host-only.
- **Fixed `n_soil_layer_max`** ⇒ no allocatables/runtime shapes; scratch (`a,b,c,r,cp,dp`) as fixed-length automatics, stack-friendly on device.
- **Issue #7 trap:** `thomas_solve` is a subroutine with `intent(out) x(n_soil_layer_max)` — never `call foo(thomas(...))`. The surface kernels bind every array result to a named array before any call. **Build nvfortran multicore on every new module** — a green ifx suite (only an `arg_temp_created` remark) is *not* sufficient; issue #7 was found in `src/biophysics` itself.
- **FPE-safe in Debug** (`-fpe0`/`-Ktrap=fp`): the inverter's `u_freeze==u_melt` banded guard; the Clausius slope's `(P − 0.378·e_sat)²` positive denominator; the dz-weighted harmonic-mean `κ_face` (both `κ>0`); the **Kersten-number `log10(max(S_r, S_r_floor))` floor + `[0,1]` clamp** (no `log10(0)` at a dry layer, no negative `κ`); `qsat`/`cp_moist` — no `0/0`, no `(neg)**real`. Constitutive + inversion kernels `pure`/`elemental`, branch-light smooth limiters instead of `case`/`zero_flow` traps (ED2's `if (can_temp == t3ple)` exact-float tests are replaced by banded checks; Bundle B/C).
- **Do NOT use `-stdpar=gpu`** (managed-allocator double-free on the allocatable-component `site`); OpenMP `target` + `-gpu=mem:separate`, `MEDS_GPU=multicore|gpu`.
- **`ENERGY_SUBSTEP_FIXED`** for the offload path (variable `nsub` is warp divergence); heat and water share the substep controller so a coupled offload sweep stays warp-uniform.

---

## 11. Validation & milestones

CTest `test_column_energy` + `test_surface_energy` (drive synthetic forcing from `defaults()`, march the kernels, assert):

1. **Energy conservation (the strongest invariant), instrumented on ALL four stores:** `Δ(Σ_k soil_energy·dz) = dt·(G_top − bottom_heat − Σ root_heat_sink)` closes to **machine precision** for the conservative flux-divergence update (§6.1), in wet and dry columns; `flux%energy_resid ≈ 0`. Leaf/wood: `Δstore_energy = dt·(Rn − H − LE + intercept)`, checked via `leaf_energy_flux_t%energy_resid ≈ 0`. CAS: `Δ(hcapcan·can_enthalpy) = dt·Σ fluxes`, checked via the `canopy_air_update` `resid` out-arg ≈ 0. The **ground+CAS latent twin** balances to round-off (identical `w_flux_gc·enthalpy_vapor(t_ground)` on both sides, §4d/§4b).
2. **Enthalpy↔temperature round-trip + continuity:** `uext_to_temp(temp_to_uext(dry_hcap, wmass, T, fliq)) == (T, fliq)` to round-off across all-ice, plateau, all-liquid; **the map `T(u)` is continuous at both `u_freeze` and `u_melt`** (`temp(u_freeze)=temp(u_melt)=t_3ple` for several `dry_hcap`, catching a wrong `tsupercool_liq`); the plateau pins `temp = t_3ple` and reads `fliq` linearly; the `wmass→0` singularity guard returns `fliq=0.5` with no trap. `d_sat_vapor_pressure_dt` matches finite-difference `de_sat/dT`.
3. **Steady state (incl. nonzero G_top and nonzero w_flux):** constant `G_top`, adiabatic bottom, `Q=0` → the column relaxes to a **linear `soil_temp(z)` profile** with `κ_face·dT/dz = G_top` at every face (Fourier steady state); a **nonzero-`w_flux`** variant checks the advective sign/splitting; a nonzero-`G_top` case pins the surface-flux sign. Leaf with `H=LE=0` relaxes to radiative equilibrium `Rn_net = 0`.
4. **Phase-change latent heat (P2a):** apply a constant cooling `G_top<0` to a wet layer at `t_3ple+ε`; assert `soil_temp` **pins at `t_3ple`** while `soil_fliq` falls 0-ward, absorbing exactly `wmass·L_f` before the temperature resumes dropping (the zero-curtain); total energy conserved through the plateau.
5. **vs a fine explicit reference:** BE-adaptive matches a 1000-substep explicit RK4 to `rtol` in a diurnally-forced column, **including a nonzero-`w_flux` case** (exposes advective splitting error / CFL); substep count drops when quiescent, refines at a sharp thermal/freezing front; no NaN/trap. Cap-hit case sets `flux%converged=.false.`, logs, and `error stop`s under `debug_error` (production contract per §5).
6. **The coupling check (closes seams, P3):** a lagged Picard weave of `veg_energy_balance` ↔ `canopy_air_update` ↔ `ground_surface_balance` ↔ `soil_energy_flux` ↔ radiation converges to a consistent `(leaf_temp, can_temp, t_ground, soil_temp(1))` fixed point; the diagnosed `t_ground` reproduces the hydrology's `ground_evaporation` behavior, the two-stream LW lower BC (#8), and the canopy emission (#9).
7. **vs ED2 (regression):** drive an ED2 single-texture column with prescribed net radiation and a root heat sink; MEDS matches the equilibrium temperature profile and ground heat flux within `rtol=1e-3` on the profile / 2% on integrated `G`. Leaf: reproduce ED2's diurnal `leaf_temp` depression/recovery pattern (qualitative, not bit-for-bit).
8. **nvfortran multicore build green** (portability gate); serial↔multicore↔GPU determinism (the 7/7 cross-backend discipline).

**Phased milestones:**

| Phase | Deliverable | Tests |
|---|---|---|
| **P0** | Types + `ENERGY_*` codes (in `meds_biophysics_types`) + CMake lib; `meds_constants` additions; `meds_thermo` inverter + Clausius slope + `air_density`/`cp_moist`; `meds_soil_thermal` (`κ`, `C_vol`, `build_soil_thermal`); reuse `thomas_solve` | 2 |
| **P1 (MVP)** | four **stateless per-store kernels** standalone: `soil_energy_flux` (BE-Thomas, conservative-energy update, `fliq≡1`) + `veg_energy_balance` (leaf/wood linearized BE) + `ground_surface_balance` + `canopy_air_update`; adaptive step-doubling; energy-conservation budget on all four stores + cap-hit contract | 1,3,5 |
| **P2a ✅** | **phase-change plateau** (inverter `_ON`, ice-aware `κ_sat(fliq)`/`C_eff(fliq)`) — self-contained, zero solver change, **unblocked**; delivered as tests only (freeze/thaw zero-curtain + ice-aware κ) | 1,3,4,5 |
| **P2b** | richer thermal, each **gated on its external dependency**: snow/sfcwater thermal stack (gated on hydrology sfcwater layers) + optional thin skin layer + sap-flow enthalpy advection (gated on dynamic plant hydraulics) | 1,4,5 |
| **P3** | **COUPLING (deferred, after this):** SoA state (`leaf_energy`/`wood_energy` cohort lockstep + `soil_energy`/`can_enthalpy`/`can_shv` 6 patch sites + creation-site `temp_to_uext` stamps + `fuse_2_patches` blend + donor-copy disturbance init); TOML `[soil]`/`[energy]`/`[pft]` + `derive_config` + `build_soil_thermal` + `validate_config` (lands with hydrology config/state); the leaf↔CAS↔ground↔soil↔radiation fixed-point orchestration in `meds_aux`; **close all 9 seams** | 6,7, full build |
| **P4** | fast-loop wire (once met forcing lands); RT→leaf-energy→gas-exchange→ground→soil→CAS weave over `dtlsm_sec`; evaporative-cooling feedback closed with the hydrology twin | 6 + driven-day energy balance |
| **P5** | nvfortran GPU parity + `ENERGY_SUBSTEP_FIXED`; multi-layer CAS (`K_diff(z)`); C-API + Python (end of biophysics dev) | 8 |

---

## 12. How this differs from ED2 / CLM / ClimaLand

MEDS **owns all four stores** (it is its own host, unlike FATES), **adopts the ClimaLand internal-energy state**, **adopts CLM's implicit BE-Thomas numerics and harmonic-mean interface**, and **keeps ED2's enthalpy-flux CAS, no-skin ground, two-source partition, and hard-pinned phase change**.

| Aspect | **ED2** | **CLM5** | **ClimaLand** | **MEDS (this design)** |
|---|---|---|---|---|
| Soil state | internal energy (RK4 tracer) | **temperature** + 2-step phase adjust | volumetric internal energy `ρe_int` | **internal energy** `soil_energy` (ClimaLand-style) |
| Leaf state | prognostic `leaf_energy` (RK4) | zero-heat-capacity **steady-state Newton** | prognostic `T_c` (heat capacity) | **prognostic `leaf_energy`**, L-stable linearized BE (ED2/ClimaLand transient) |
| CAS | specific enthalpy `can_enthalpy` (+`can_shv`) | diagnostic weighted-mean `T_s`/`q_s` (zero-capacity node) | — | **specific enthalpy** (ED2), `hcapcan ≡ wcapcan` |
| Soil integration | explicit adaptive RK4, co-integrated w/ canopy (stiff) | implicit BE/CN tridiagonal | implicit IMEX Jacobian | **implicit BE-Thomas**, operator-split (advection lagged), stateless over fixed `n_soil_layer_max` |
| Interface κ | geometric/log-linear (overestimates) | **harmonic** (series resistors) | harmonic | **dz-weighted harmonic** (CLM) |
| Phase change | enthalpy read-off (`uextcm2tl`) | temperature excess-energy 2-step | τ-relaxation (tuning param) | **enthalpy read-off** (ED2/ClimaLand) — automatic in inverter, no adjust step, no τ |
| Frozen κ | uses liquid κ regardless (bug) | ice-aware Kersten | ice-aware Balland-Arp | **ice-aware** `κ_sat(fliq)`, Kersten clamped `[0,1]` |
| Advective heat | yes (`qw_flux`, liquid) | **omitted** (T-only diffusion) | yes (`ρe_int,l K∇h`) | **yes** (ClimaLand; nearly free with energy state; explicit/CFL-bounded) |
| Ground skin | none — `t_ground = topsoil_temp` | none — top-layer flux BC | Neumann on `ρe_int` | **none — `t_ground = soil_temp(1)`** (ED2), optional thin skin P2b |
| Ground LE twin | vapour enthalpy at `t_ground` | — | liquid-enthalpy consistent | **`w_flux_gc·enthalpy_vapor(t_ground)` on BOTH surface & CAS** (no `λE` mismatch) |
| Bottom heat BC | zero-flux (implicit in loop bound — stale hazard) | zero-flux / geothermal | Neumann | **explicit** zero-flux geothermal / prescribed T |
| `ρ cp`/`ρ` in flux drivers | baked into `gbh`/`ggnet` | explicit | explicit | **explicit** (avoids double-count; ED2 bug #1; `ρ_air` on transp too) |
| Coupling | nested RK4 sub-step iteration | outer M-O / inner Newton | coupled IMEX Jacobian | **deferred P3** Picard fixed point; kernels stateless w/ lagged siblings now |
| Platform | stateful, CPU | host land model, CPU | Julia, CPU/GPU | **stateless kernel; CPU + GPU** (`pure`/`elemental`, OpenMP target) |

**One line:** MEDS = ClimaLand's internal-energy state (closed-form phase read-off) executed with CLM5's implicit BE-Thomas numerics (dz-weighted harmonic interface, analytic Clausius linearization), wrapped around ED2's enthalpy-flux CAS, no-separate-skin ground, and hard-pinned phase change — stateless, GPU-friendly, energy-conserving to round-off, closing all nine forced-temperature seams, with the coupled fixed point deferred to P3.

---

## 13. Bugs / quirks found in the ED2 reference

Curated from the leaf/wood, CAS, and soil-thermal extractions; flagged to fix-by-construction or carry deliberately.

**E1 — `ρ_air·cp` baked into `gbh`/`ggnet` (fix by explicit multiply).** ED2's `leaf_gbh`/`wood_gbh`/`ggnet·can_rhos·can_cp` already carry air density and specific heat (`h_flux_lc = 2·LAI·gbh·ΔT`, no explicit `cp`). **MEDS keeps aerodynamic conductances pure `[m/s]` and multiplies `ρ_air·cp_air` at the flux site** — the conceptual form. Likewise the transpiration/evap mass drivers multiply `ρ_air` explicitly (§4a). Do not double-count (Bundle A bug #1).

**E2 — Wood area `π·WAI` vs leaf `2·LAI` (carry both).** Leaves are two flat plates (`effarea_heat=2·LAI`); branches are cylinders (`π·WAI`). Easy to accidentally unify; keep distinct (Bundle A bug #2).

**E3 — `effarea_transp` is per-PFT and transpiration-only.** 2 for amphistomatous conifers, else 1; multiplies transpiration only, while sensible heat always uses 2. A single "leaf-sidedness" constant is wrong — carried on `veg_thermal_params_t` (Bundle A bug #3).

**E4 — Enthalpy vs internal-energy asymmetry (mis-books latent heat).** Evaporation/transpiration use vapour **enthalpy** `enthalpy_vapor` (flow work + latent); sap flow uses liquid **internal energy** `internal_energy_liquid` (no phase change on internal advection). Mixing them mis-books latent heat (Bundle A bug #4). MEDS provides two named helpers and drives **ground LE with the same `enthalpy_vapor(t_ground)` twin on both surface and CAS** (§4d) — no `λE` shortcut.

**E5 — Heat capacity convention must be picked once.** ED2's dry `hcap` excludes prognostic water (added in `uextcm2tl` via `wmass`) under dynamic hydraulics, but folds bound water into tissue specific heats under static hydraulics. Mixing double-counts internal-water heat capacity. **MEDS puts ALL store water in `wmass`, consistently** (Bundle A bug #5).

**E6 — `cph2o = 1859` is vapour cp, not liquid (`4186`) — in ALL CAS terms.** Every `qw_flux`/`enthalpy_vapor`/`h_flux_ac` uses the vapour reference line with the single `(cp_vap, tsupercool_vap)` pair; a naive re-derivation plugging in `cp_liq` corrupts both state recovery and the sensible/latent split. **Keep the single pair everywhere** (Bundle B bug #1).

**E7 — Sensible and latent are inseparable in the CAS state.** `h_flux_ac` is diagnostic-only, reconstructed from `e_flux_ac − (−w_flux_ac·cp_vap·tsupercool_vap)`. **Do not attempt an `H+λE` formulation** — port the single-`e_flux_ac` enthalpy-flux form or double-count/drop the latent baseline (Bundle B bug #2).

**E8 — Vapour-enthalpy twins use SOURCE temperature, not `can_temp`.** `qw_flux*` evaluates enthalpy at `leaf_temp`/`t_ground`/`sfcwater_tempk`, matching the donor pool's energy equation which subtracts the same twin. Using `can_temp` violates the donor's energy balance (Bundle B quirk). **The twin must use identical `T_source` on both sides of the interface** — enforced for the ground evap twin (§4d/§4b).

**E9 — No geothermal / deep-T heat BC; bottom flux zero implicitly.** ED2's bottom conductive flux is zero only because the loop starts at `klsl+1` plus per-step zeroing — a stale non-zero flux would leak if the derivative ran without re-zeroing. **MEDS sets `κ_face(nzg)=0` (or a config geothermal value) explicitly** (Bundle C bug #1).

**E10 — Thermal conductivity ignores ice (uses liquid κ regardless of phase).** ED2's `k_th(θ)` uses `k_water=0.57` even when frozen; ice conducts ~4× better (`k_ice=2.29`). `soil_fliq` is available but unused. **MEDS's `κ_sat(fliq)` is ice-aware from the start**, with the Kersten number clamped to `[0,1]` (Bundle C bug #2). Fine for tropics; correct for freezing soils.

**E11 — Interface κ by geometric interpolation, not harmonic.** Conduction through series layers should combine as the thickness-weighted harmonic mean; ED2's geometric interpolant overestimates when adjacent layers differ strongly. **MEDS uses the dz-weighted harmonic mean** `(dz_k+dz_{k+1})/(dz_k/κ_k+dz_{k+1}/κ_{k+1})` (CLM; Bundle C bug #3).

**E12 — `slcpd` air term is a frozen midpoint approximation.** ED2 assumes moisture permanently halfway between `soilcp` and `slmsts` for the air heat-capacity contribution. Small term (`air_hcapv≈1.2e3` vs mineral `~2e6`), so minor — MEDS computes `C_eff` from the actual current `θ, fliq` (Bundle C bug #4).

**E13 — Advected water enthalpy assumes liquid (`fliq=1`), upwind-only.** Consistent with liquid Darcy flow but off if ice mobilizes; first-order upwind is numerically diffusive but stable, and carries its own CFL (bounded explicitly, §6.1). **MEDS keeps liquid-upwind** (matches liquid-only P1; revisit if frozen transport is modeled) (Bundle C bug #5).

**E14 — `t3ple` exact-float equality tests are brittle.** ED2's `if (can_temp == t3ple)` / `if (soil_temp == t3ple)` are exact comparisons at the triple point. **MEDS uses banded `abs(·) < eps` / `abs(u_melt−u_freeze) < eps` tests** (Bundle B/C).

**E15 — Skin = top-layer temperature limits diurnal surface-T amplitude.** No sub-grid skin solve means surface-T amplitude is capped by the top-layer heat capacity; a thin top layer is needed for realistic `t_ground`/LW emission. **MEDS keeps `dz(1)` small; an explicit skin layer is an optional P2b upgrade** (Bundle C bug #7).

**E16 — Shedding disabled in-derivative.** ED2's `wshed≡0` in the tendency ("TURNING OFF SHEDDING"); real shedding is in `adjust_veg_properties`. **MEDS replicates the split** (interception enthalpy in the fast step, shedding in the slow adjust), not the literal dead `qwshed` (Bundle A bug #6).

**E17 — `veg_hcap_min` resolvability floor prevents `1/hcap` blowup.** Below the per-PFT heat-capacity floor a cohort is slaved to `can_temp`; the resolvable↔unresolvable switch is booked conservatively by `temp_to_uext`. **MEDS carries the floor (on `veg_thermal_params_t`) and the conservative re-derivation** (Bundle A §5).

---

## 14. Phased implementation plan

- **P0 — foundation.** `meds_biophysics_types` extended (`soil_energy_column_t`, `cas_state_t`, `leaf_energy_env_t`/`leaf_energy_flux_t`, `energy_forcing_t`, `energy_flux_t`, `soil_thermal_params_t`, `veg_thermal_params_t`, `cas_atm_forcing_t`, `energy_opts_t`; `ENERGY_*` codes **beside the existing `SOIL_*` codes**); `meds_constants` additions (D-table, incl. the corrected `tsupercool_liq`); `meds_thermo` extended (`uext_to_temp`/`temp_to_uext` inverter, `d_sat_vapor_pressure_dt`, `air_density`, `cp_moist`, `enthalpy_vapor`, `internal_energy_liquid`); `meds_soil_thermal` (`pure`/`elemental` `κ`, `C_vol`, `build_soil_thermal`); reuse `thomas_solve`. Tests: inverter round-trip + **continuity at `u_freeze`/`u_melt`** (all-ice/plateau/all-liquid), Clausius slope vs finite difference. CMake lib builds; **nvfortran multicore green**.
- **P1 — MVP, four stateless per-store kernels standalone.** `soil_energy_flux` seam + bare-array `soil_heat_be_step` (BE-Thomas, dz-weighted harmonic κ, advective heat, conservative-energy update with correct flux-divergence sign, `fliq≡1`, adaptive step-doubling, energy budget on all four stores + cap-hit contract); `veg_energy_balance` (leaf/wood linearized BE, explicit `ρ cp` and `ρ_air`, `−8εσT³` LW damping, `π·WAI` vs `2·LAI`, `veg_hcap_min` floor); `ground_surface_balance` (`t_ground = soil_temp(1)` skin, `G_top` with `enthalpy_vapor(t_ground)` LE twin, no `dt`); `canopy_air_update` (CAS enthalpy+`can_shv` single step from mass+enthalpy fluxes, `hcapcan ≡ wcapcan`, source-T twins). Tests 1,3,5. Each on its CTest, siblings lagged/forced.
- **P2a — phase change (self-contained, unblocked). ✅ IMPLEMENTED.** Inverter plateau `_ON` + ice-aware `κ_sat(fliq)`/`C_eff(fliq)` — freeze/thaw automatic, **zero solver change**: the design's promise held — *no kernel edit was needed*, the P1 code already gated `fliq_use = fl_n(k)` on `opts%phase_change` and the ice-aware constitutive laws + 3-branch inverter were in place from P0/P1. P2a delivered as **tests only** (`test_column_energy`): the zero-curtain (single-layer wet control volume, sealed bottom, steady `G_top` cooling → energy removed per step is exactly `|G_top|·dt`), asserting `soil_temp` pins at `t_3ple` to round-off while `soil_fliq` falls 1→0 absorbing exactly `wmass·L_f` (0.3% of expected), plus the thaw mirror and an ice-aware-conductivity unit check (`κ(fliq=0) > κ(fliq=1)`). Green ifx + nvfortran mc. Tests 1,3,4,5.
- **P2b — richer thermal (each gated on an external dependency).** Snow/sfcwater thermal stack (**gated on the hydrology module's sfcwater layers**); optional explicit thin skin layer; sap-flow enthalpy advection (**gated on dynamic plant hydraulics**). Tests 1,4,5.
- **P3 — COUPLING (deferred, after this).** SoA state through the cohort lockstep (`leaf_energy`/`wood_energy`) + the **6** patch sites (`soil_energy`/`can_enthalpy`/`can_shv`, incl. `fuse_2_patches` area-blend at the `recruit_pool` precedent + donor-copy disturbance init) + creation-site `temp_to_uext` stamps; `[soil]`/`[energy]`/`[pft]` TOML + presence map + `derive_config` (plain arrays) + `build_soil_thermal` (typed) + `validate_config` — landing **together with the hydrology config/state**; the **leaf↔CAS↔ground↔soil↔radiation fixed-point orchestration** in `meds_aux` closing **all 9 seams** (radiation↔energy iteration for seams #8 and #9). Tests 6,7, full engine build.
- **P4 — fast-loop wire (blocked on met forcing).** Fast-loop driver in `meds_aux`, called from `advance_one_step` before `vegetation_dynamics`; RT→leaf-energy→gas-exchange→ground→soil→CAS weave over `dtlsm_sec`; evaporative-cooling feedback closed with the hydrology CAS-water twin. Needs the met reader (`src/io`). Validation: driven-day closed energy balance.
- **P5 — platform + extension.** nvfortran GPU parity + `ENERGY_SUBSTEP_FIXED`; multi-layer CAS (`K_diff(z)`, per-layer `hcapcan`); C-API + Python (end of biophysics dev). Test 8.

---

## 15. Open questions

1. **Prognostic-leaf vs steady-state-leaf.** *Resolved (§4a):* prognostic `leaf_energy` (ED2-faithful, phase-safe) with an L-stable linearized BE step removing the stiff-mode objection. Revisit a fully-coupled IMEX leaf only if MEDS moves to a monolithic implicit land solver.
2. **Skin layer.** Keep ED2's `t_ground = soil_temp(1)` (MVP) or add an explicit thin skin for diurnal surface-T / LW-emission fidelity? *Recommend:* keep `dz(1)` small at P1; add a skin layer at P2b if the two-stream LW lower BC (seam #8) proves sensitive.
3. **Phase-change mechanism.** *Resolved (§6.1):* enthalpy read-off in the inverter (ED2/ClimaLand), not CLM's temperature 2-step nor ClimaLand's τ-relaxation. The one choice left: freezing-point depression (supercooled `θ_liq,max`, `tsupercool`) — adopt the smooth Dall'Amico form at P2a or keep the sharp `t_3ple` pin? *Recommend:* sharp pin at P2a (parameter-lean), depression later if needed.
4. **CAS depth.** Fixed `can_depth` (MVP) vs variable-depth control volume with the `ddens_dt_effect` work term. *Recommend:* fixed at P1 (flag the small diurnal-density non-conservation); variable-depth at P2b with the multi-layer CAS.
5. **Soil conductivity model.** de Vries rational (ED2) vs Johansen/Farouki Kersten (CLM) vs Balland-Arp (ClimaLand). *Recommend:* Johansen/Farouki (CLM), ice-aware, Kersten clamped `[0,1]`, dz-weighted harmonic interface — best-documented and matches the BE-Thomas numerics choice.
6. **Numerics tolerances/defaults:** `atol=1e-2 K`, `rtol=1e-3`, `max_substep=200`, `dtlsm_sec=900`? Calibrate against the fine explicit reference (test 5); heat and water should share the substep controller (with the advective-CFL floor when `w_flux≠0`).
7. **Advective-heat coupling to hydrology.** The soil-heat `qw_flux` needs the hydrology's inter-layer `w_flux(k)` as forcing. Does the P3 orchestration run hydrology-then-heat (heat sees the just-updated `w_flux`) or lag one sub-step? *Recommend:* hydrology-then-heat within a sub-step (heat is downstream of water), lagging only across the outer coupling Picard; either way the explicit advective term respects its CFL.
8. **Tissue specific heats and `veg_hcap_min` values** per PFT — the least-constrained new parameters; seed from ED2 `ed_params.f90` and Longo 2019, expose in TOML.

---

## 16. References

- **Lawrence et al. (2019, 2020)** *CLM5.0 Technical Note* (NCAR) — Ch. 2.5 (leaf/two-source flux Newton solve, Clausius linearization eqs. 2.5.128–2.5.157) and Ch. 2.6 (implicit tridiagonal soil/snow heat, Crank–Nicolson/BE, Farouki/Johansen conductivity, phase change). **Numerics + BCs (§4c, §6).**
- **Johansen (1975)** & **Farouki (1981)** — Kersten-number soil thermal-conductivity scheme (valid `S_r ≳ 0.05–0.1`; MEDS clamps). **Soil κ (§4c).**
- **Balland & Arp (2003)** *J. Environ. Eng. Sci.* — thermal conductivity used by ClimaLand soil. **Cross-reference (§4c).**
- **Deck, Kohler et al. (2026)** *J. Adv. Model. Earth Syst.* 10.1029/2025MS005118 — ClimaLand: prognostic `ρe_int` soil energy, closed-form `T(ρe_int, θ_i)`, advective liquid-enthalpy term, prognostic canopy `T_c`, τ-relaxation phase change. **State-design source (§1, §3, §6).**
- **ClimaLand.jl source** (`main`): `src/standalone/Soil/energy_hydrology.jl`, `soil_heat_parameterizations.jl` (`volumetric_heat_capacity`, `temperature_from_ρe_int`, `thermal_conductivity`, `phase_change_source`), `src/standalone/Vegetation/canopy_energy.jl`. **Verified closures (§3, §4).**
- **Longo et al. (2019)** *GMD* 12:4309 — ED-2.2 technical description; the vegetation/CAS/soil energy budgets this family reimplements. **Reference physics (§4, §13).**
- **ED2 source** (`../ED2/ED/src`): `dynamics/rk4_derivs.f90` (`leaftw_derivs` — leaf ~1411–1770, wood ~1790–2050, sap-flow/veg-sum 2108–2145, CAS closure 2156, soil conduction 447/527, advection 691/737, surface BC 513); `utils/therm_lib.f90`/`therm_lib8.f90` (`uextcm2tl`/`uint2tl`/`tl2uint`/`cmtl2uext`/`tq2enthalpy`/`hq2temp`); `utils/ed_therm_lib.f90` (`calc_veg_hcap` 67, `update_veg_energy_cweh` 156, `ground_temp = topsoil_temp` 677); `utils/stable_cohorts.f90` (`is_resolvable`); `init/ed_params.f90` (`effarea_*`, `veg_hcap_min`, `slcpd`/`thcond0..3` 2260–2295); `dynamics/canopy_struct_dynamics.f90` (`ed_stars8`, `can_whccap8`, `leaf_gbh` fallback ~729); `memory/consts_coms.F90` (`cliq/cice/cph2o/cpdry/alli/t3ple/tsupercool_*` — pull numeric values when porting, especially `tsupercool_liq` from `cmtl2uext`, absent in the extracted checkout). **The extracted reference (§4, §13).**
- **MEDS internal:** `archive/MEDS_COLUMN_HYDROLOGY_DESIGN.md` (the direct sibling — stateless-kernel + state/process-wall pattern, `meds_soil_solver` Thomas reuse, negative-z geometry, adaptive step-doubling, conservative flux-divergence update, `build_soil_params`, the forced `t_ground` seam this family closes); `archive/MEDS_HYDRAULICS_DESIGN.md` (the stateless + L-stable BE + boundary-only-conservation idiom, the `hydro_env_t%leaf_temp` seam); `src/biophysics/meds_biophysics_types.f90:30–116` (the existing `SOIL_*` codes the `ENERGY_*` codes sit beside), `meds_demography_dynamics.f90:52` (`growth_step`, the real bare-array OpenMP-target precedent), `meds_demography_fusefiss.f90:435` (`recruit_pool` blend precedent), `meds_canopy_radiation.f90:35/76` + `meds_twostream_band.f90:71–73` (the `canopy_temp` emission seam #9); `CLAUDE.md` (state/process wall, no-hard-coded-parameters, naming convention, issue #7 nvfortran array-result trap, `-stdpar=gpu` prohibition, the reserved master-loop follow-up this family's P3 coupling belongs to).
