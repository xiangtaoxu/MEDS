# MEDS biophysics ⇄ driver de-duplication — design (design-only, no code changes)

Status: **PARTIALLY IMPLEMENTED** (branch `refactor/biophysics-dedup`). Goal: **each physical
process lives ONCE in `biophysics`/`biogeochemistry` as pure flux + two-form (tendency /
implicit-step) kernels; the drivers only orchestrate — and BOTH numeric schemes (operator-split and
IMEX-ARK) are kept.** The de-dup shares the *physics*, not the *scheme*.

**Done (bit-identical, ifx 34/34 + nvfortran 34/34):**
- Heterotrophic respiration moved to `meds_soil_biogeochem` (§4); dead `column_co2_step` deleted.
- **CAS box two-form** (`cas_column_time_deriv` / `cas_column_step_implicit`) extracted to
  `meds_cas_biophysics` and called by BOTH integrators — the CAS enthalpy/humidity/CO2 box math no
  longer inlined in `surface_derivs` (ARK) or `column_fast_step` (split); the scheme-specific source
  adjustments (ARK condensation sink, split snow sublimation) stay in the drivers' source assembly.

**Deferred follow-ups:**
- `leaf_energy_diagnostic` + `ground_surface_fluxes` extraction (the leaf/wood + ground skin are
  still inlined; the split path's snow-blend orchestration makes these thinner-value + higher-care).
- Removing `canopy_air_update` / `canopy_air_co2_update`: retained as isolated **budget-closure**
  reference kernels (they also return the closed-budget residual + NEE/NEP that the bare box kernels
  don't); removing them means migrating those unit-test assertions.

The remainder of this doc is the original plan / target structure.

## 1. The overlap — surface/CAS physics is triplicated

`meds_column_derivs.surface_derivs` (the ARK RHS) and `meds_column_dynamics.column_fast_step` (the
split step) each **inline** the leaf/wood/ground/CAS physics, and `src/biophysics` *also* ships
stand-alone kernels for the same physics. The soil-heat / soil-water / plant-hydraulics tendencies
are the ONE part done right — the ARK RHS delegates them to the kernels. Everything above the soil
is duplicated:

| physics | biophysics kernel | `surface_derivs` (ARK, live) | `column_fast_step` (split, live) |
|---|---|---|---|
| leaf energy (diagnostic) | `veg_energy_step_implicit` | inline (`:225-233`) | inline |
| wood energy (diagnostic) | `veg_energy_step_implicit(is_leaf=F)` | inline (`:243-248`) | calls kernel (`:550`) |
| ground skin | `ground_surface_balance` | inline (`:254-256`) | inline (`:644-646`) |
| CAS enthalpy + humidity | `canopy_air_update` | inline (`:274-275`) | inline |
| CAS CO₂ twin | `canopy_air_co2_update` | inline (`:276`) | inline |
| **soil heat** | `soil_energy_time_deriv` | **calls kernel** (`:327`) ✓ | `soil_energy_step_implicit` ✓ |
| **soil water** | `soil_water_time_deriv` | **calls kernel** (`:333`) ✓ | `column_hydrology_flux` ✓ |
| **plant hydraulics** | `plant_water_tendency` | **calls kernel** (`:343`) ✓ | (calls kernel) ✓ |

## 2. Principle — physics once, both schemes call it

The soil block shows the target pattern: a pure `*_time_deriv` (RHS) for ARK and a `*_step_implicit`
(backward-Euler advance) for the split, over a **shared flux core**, all in biophysics. Extend that
**two-form API** to the surface block. The key realization for keeping both schemes:

- The **fluxes** (leaf/wood diagnostic temps + `coh_h`/`coh_qw`/`transp`, the ground skin `h`/`le`/
  `g_top`) are **scheme-agnostic** — identical for split and ARK.
- Only the **CAS box update** differs by scheme, and only in the *stepping*: ARK wants the tendency
  `d_cas = (src + gatm·(atm−state))/cap`; the split wants the implicit box step
  `state^{n+1} = (state + dt·(src+gatm·atm)/cap)/(1 + dt·gatm/cap)`. Both consume the **same** `src`.

So a shared flux kernel + a two-form CAS box lets each integrator keep its own stepping while the
physics lives once. The scheme stays a driver toggle (`integration_scheme = split | ARK`).

## 3. Target code structure

**Interface convention (see §3d):** kernels take **POD derived-type bundles** (flat: scalars +
fixed-size arrays, NO `allocatable`/`pointer` components), matching every existing biophysics kernel;
field/argument names are spelled out for readability (established domain tokens — `cas_enthalpy`,
`can_shv`, `dbh`, `nplant`, `pft`, loop `i`/`ip` — are kept).

### 3a. `src/biophysics` — pure physics kernels

**`meds_vegetation_biophysics`** — leaf/wood surface energy (diagnostic + prognostic two-form):

```fortran
! Shared flux basis at a given tissue temperature (already exists; keep).
pure subroutine veg_surface_fluxes(tissue_temp, env, tparams, is_leaf, sensible_flux, ...)

! DIAGNOSTIC mode: the linearized instantaneous leaf/wood balance (extracted verbatim from
! surface_derivs:225-248). Pure; both drivers call it per cohort.
type :: leaf_diag_env_t          ! per-cohort frozen inputs (from surface_frozen_t)
   real(wp) :: sensible_coeff, transp_conductance, absorbed_shortwave, absorbed_net_longwave
   real(wp) :: area_index, emissivity
end type
type :: leaf_diag_out_t
   real(wp) :: tissue_temp, sensible_flux, latent_flux, soil_shed_enthalpy
   real(wp) :: transpiration, net_radiation
end type
pure subroutine leaf_energy_diagnostic(env, cas_temp, cas_shv, cas_qsat, cas_dqsat_dtemp,  &
                                       air_density, is_leaf, out)

! PROGNOSTIC two-form store (step exists; add the tendency twin for symmetry with soil).
pure subroutine veg_energy_time_deriv  (store, env, tparams, is_leaf, energy_tendency)
pure subroutine veg_energy_step_implicit(store, env, tparams, dt, is_leaf, flux)   ! exists
```

**`meds_ground_biophysics`** — ground skin (a pure flux; scheme-agnostic):

```fortran
! Takes the FROZEN soil_evaporation (the single hydrology authority), NOT a re-computed
! evaporation — this is why it supersedes the dead ground_surface_balance, which recomputed le.
pure subroutine ground_surface_fluxes(ground_temp, cas_temp, ground_conductance, air_density, &
                                      soil_evaporation, absorbed_shortwave_ground,            &
                                      absorbed_net_longwave_ground,                           &
                                      ground_sensible_flux, ground_latent_flux, soil_top_heat_flux)
! snow kernels unchanged (snow_cover_fraction / snow_accumulate / snow_drain_meltwater /
!                        snow_energy_step + private surface_fluxes/base_conductance)
```

**`meds_cas_biophysics`** — the canopy-air-space **column** in two forms (surface analogue of
`soil_energy_*`). NOTE: a single well-mixed layer today, but `cas_column_t` is shaped so a
**multi-layer canopy-air column** (K-theory Thomas solve, `n>1`) is the general case and the box is
the `n=1` degenerate one (design `MEDS_COLUMN_CO2_BALANCE_DESIGN.md §4c′`); field names read
per-layer-ready:

```fortran
type :: cas_column_t             ! frozen canopy-air-space params (per patch, per substep)
   real(wp) :: air_mass_capacity, air_molar_capacity          ! wcap, ccap  [kg/m2],[mol/m2]
   real(wp) :: atm_conductance_enthalpy, atm_conductance_vapor, atm_conductance_co2
   real(wp) :: atm_enthalpy, atm_specific_humidity, atm_co2, pressure, condensation_timescale
end type
type :: cas_source_t             ! summed surface sources into the CAS
   real(wp) :: surface_enthalpy_source, surface_vapor_source, biotic_co2_source
end type

! Private shared core: applies the smooth condensation sink + forms the column rates.
pure subroutine cas_column_rates(state, source, column, d_enthalpy, d_shv, d_co2, condensation)

pure subroutine cas_column_time_deriv  (state, source, column, tendency)          ! ARK  (:270-276)
pure subroutine cas_column_step_implicit(state, source, column, dt, updated, residual)  ! split
                                          ! (replaces canopy_air_update + canopy_air_co2_update)
! aggregate_cohort_co2_fluxes: KEEP as the cohort GPP/Ra -> patch summing helper the driver calls.
```

### 3b. `src/biogeochemistry` — carbon fluxes by domain

**`meds_soil_biogeochem`** gains the fast heterotrophic-respiration kernels (moved from
`meds_cas_biophysics`), beside the slow `heterotrophic_respiration_matrix`:

```fortran
pure function heterotrophic_respiration_flux(fast_soil_carbon, soil_temp, theta, theta_dry, &
                                             theta_sat, opts) result(rh)      ! moved here
pure function heterotrophic_respiration_damm(...) result(rh)                  ! moved here
! co2_opts_t + damm_params_t move to meds_biogeochem_types (from meds_biophysics_types).
```

No `biophysics→biogeochemistry` edge appears: the ONLY remaining caller is the driver (which links
both libraries) — see §4. The driver is then the single Rh authority (Rh from biogeochem + GPP/Ra
from plant → `nee_biotic`).

### 3c. `src/driver` — orchestration only

The two integrators keep their frozen-coefficient assembly + their stepping, but delegate all surface
physics. `surface_derivs` (ARK RHS) becomes thin:

```fortran
pure subroutine surface_derivs(state, frozen, n, tend):
   cas_temp = cas_temp_of_enthalpy(state%cas_enthalpy, state%cas_shv)
   cas_qsat = sat_specific_humidity(cas_temp, frozen%pressure)
   cas_dqsat_dtemp = sat_specific_humidity_temp_deriv(cas_temp, frozen%pressure)
   source = 0
   do i = 1, n:
      call leaf_energy_diagnostic(frozen%leaf(i), cas_temp, state%cas_shv, cas_qsat,        &
                                  cas_dqsat_dtemp, frozen%air_density, .true. , leaf_out)
      call leaf_energy_diagnostic(frozen%wood(i), cas_temp, state%cas_shv, cas_qsat,        &
                                  cas_dqsat_dtemp, frozen%air_density, .false., wood_out)
      tend%leaf_temp(i) = leaf_out%tissue_temp;  tend%transp_c(i) = leaf_out%transpiration
      source%surface_enthalpy_source += leaf_out%sensible_flux + leaf_out%latent_flux       &
                                      + wood_out%sensible_flux
      source%surface_vapor_source    += leaf_out%transpiration
   call ground_surface_fluxes(frozen%ground_temp, cas_temp, frozen%ground_conductance,      &
                              frozen%air_density, frozen%soil_evaporation,                   &
                              frozen%absorbed_shortwave_ground, frozen%absorbed_net_longwave_ground, &
                              ground_sensible_flux, ground_latent_flux, tend%soil_top_heat_flux)
   source%surface_enthalpy_source += ground_sensible_flux + ground_latent_flux
   source%surface_vapor_source    += frozen%soil_evaporation
   source%biotic_co2_source        = frozen%biotic_co2_source
   call cas_column_time_deriv(state, source, frozen%cas, tend)   ! d_cas_enthalpy/shv/co2 (+ cond)
```

`column_fast_step` (split) uses the SAME kernels inside its Picard loop, then the implicit step:

```fortran
Picard iterate k = 1..niter:
   do i: call leaf_energy_diagnostic(...) ; accumulate source     ! same kernels as ARK
   call ground_surface_fluxes(...)                                ! same kernels as ARK
   call cas_column_step_implicit(cas_state, source, cas_column, dt_fast, cas_updated, residual)
```

**Stays in the driver (orchestration):** `build_column_frozen` / the per-substep pre-pass that fills
`cas_column_t` + the per-cohort `leaf_diag_env_t` by calling `canopy_aerodynamics`, `canopy_radiation`,
the plant `leaf_gas_exchange`, `column_hydrology_flux` (the `soil_evaporation` authority), and Rh; the
ARK tableau (`meds_ark_stepper`); the split sequencing + Picard loop; the closed-budget ledgers.

### 3d. Interface convention — POD bundles + a per-domain Python capi

- **Fortran kernels** take POD derived-type bundles (readable, grouped, matching all existing
  biophysics kernels). Keep them **flat** (scalars + fixed-size arrays only) — no `allocatable`/
  `pointer` components — so the two boundaries that need contiguous data can marshal them
  mechanically. The bare-array conversion for OpenMP-`target` offload stays deferred (P5, CLAUDE.md);
  these surface kernels are not on the offloaded hot path today.
- **Python exposure** is a per-domain `bind(c)` capi (a thin flattening wrapper), exactly like
  `meds_plant_capi` flattens `cfg%pft` → `leaf_photo_params_t` for `meds_leaf_solve`. Because the
  domains are stateless pure-kernel seams, **radiative transfer, CAS-column dynamics, and soil
  dynamics can each get an independent capi entry point and be called separately from Python** — the
  caller supplies the inputs the Fortran driver would otherwise assemble, and gets fluxes/tendencies
  back. The POD bundles do not block this; the capi wrapper marshals bare arrays into them. (This is
  the reason for the "flat, no allocatables" rule above.)

## 4. Repurpose vs delete (corrected from the first draft)

Most "dead" kernels are dead *because* the drivers inline instead of calling them — the fix is to
make them the **live shared kernels**, not to remove them:

| kernel | disposition |
|---|---|
| `canopy_air_update` | **repurpose** → `cas_column_step_implicit`; wire the split path to call it (live) |
| `canopy_air_co2_update` | **fold** into the CO₂ branch of the two `cas_column_*` (live in both) |
| `ground_surface_balance` | **replace** with `ground_surface_fluxes` (frozen `soil_evap`; live in both) |
| leaf/wood inline (both drivers) | **extract** → `leaf_energy_diagnostic` (live in both) |
| `aggregate_cohort_co2_fluxes` | **keep** as the driver's cohort-CO₂ summing helper |
| `heterotrophic_respiration_flux`/`_damm` | **move** to `meds_soil_biogeochem` (§3b) |
| `column_co2_step` | **delete** — genuinely redundant all-in-one assembler (0 callers) |
| `intercept_canopy_layer` | **separate call**: wire into the throughfall cascade, or drop (test-only, `sigma_w` unused) |

## 5. Ordering — extract → rewire → delete, one store at a time

Deleting first only removes reference code while the physics stays inline (no progress). Per surface
store: (1) **extract** the driver-inline physics into the new pure kernel, asserting it reproduces
the inline result; (2) **rewire** `surface_derivs` AND `column_fast_step` to call it (delete both
inline copies); (3) **delete/repurpose** the old advance-commit kernel. Recommended order:
**ground → CAS box → leaf/wood → Rh-to-biogeochem**. Rebuild + assert the golden anchor
(`tc_split(54)=292.450065`) + the per-store closed budgets **bit-identical on ifx AND nvfortran**
after each store (issue #7).

**Tests:** coverage already lives in the coupled tests (`test_column_derivs` validates `surface_derivs`
bit-for-bit vs the split; `test_column_dynamics` validates the split). Replace the deleted unit tests
(`test_surface_energy` CAS/ground cases, `test_column_co2`) with small unit tests on the three new
kernels (`leaf_energy_diagnostic`, `ground_surface_fluxes`, `cas_column_*`) so isolated coverage is not lost.

## 6. Appendix — the gated ARK-only end-state (NOT this refactor)

The deepest possible de-dup is a *single* integration path: retire the split `column_fast_step` and
keep only IMEX-ARK, removing the second inline copy entirely. This is **explicitly out of scope here**
and stays gated on ARK reaching feature parity (aquifer / Zeng–Decker, which `column_fast_step_ark`
hard-errors on) — the split remains the feature-complete, validated reference. §2–§5 keep **both**
schemes; this appendix is recorded only so the two-form kernels above are chosen to make an eventual
ARK-only collapse a driver-side deletion, not a physics rewrite.

## 7. Risk

- **Bit-identity.** §5 is behaviour-preserving only if each extracted kernel reproduces the inline
  arithmetic exactly (same operation order). Golden anchors + closed-budget residuals are the guard;
  build both compilers.
- **Scheme fidelity.** The shared flux kernels must feed both the ARK tendency and the split implicit
  box without changing either — verified by the existing per-reservoir `test_column_derivs` equalities.
- **Scope.** Touches `meds_column_derivs`, `meds_column_dynamics`, `meds_ark_stepper` (frozen-type
  fields), `meds_cas_biophysics`, `meds_ground_biophysics`, `meds_vegetation_biophysics`,
  `meds_soil_biogeochem`, `meds_biogeochem_types`, and ~3 tests — a focused refactor, its own PR.
