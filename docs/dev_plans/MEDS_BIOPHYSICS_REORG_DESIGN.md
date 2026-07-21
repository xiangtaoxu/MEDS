# MEDS biophysics module reorganization — design (design-only, no code changes)

Status: **IMPLEMENTED** (branch `refactor/biophysics-surface-subsystems`) — the surface-subsystem
split (P0–P3) + the two-form rename API + the `meds_biophysics_interface` façade landed as
**bit-identical** moves/renames (ifx Debug 34/34, nvfortran multicore 34/34). The operator-**split**
integrator is **retained** (see §7). Deferred to follow-ups: the §4 shared-helper *dedup* (each
`*_time_deriv` / `*_step_implicit` pair currently keeps its own faces/κ code — renamed, moved, but not
yet factored into the shared `*_flux_divergence` helper, to guarantee bit-identity this pass), the
physics/numerics *collapse* (§7), and a coupled `meds_soil_dynamics`.

Original goal: reorganize `src/biophysics/` from process-lumped `meds_column_*` files to
**surface-subsystem** modules, with a shared optics library and a library **interface façade**.

## 1. Current state (what we are refactoring)

`src/biophysics/` today (post the 2026-07 reorg):

| file | contents | kernel style |
|---|---|---|
| `meds_optics` | leaf-angle Beta LIDF, `G(μ)`, `ω/g` scattering, ground optics | pure |
| `meds_twostream_band` | one-band adding solver | pure |
| `meds_canopy_radiation` | **RT seam** `canopy_radiation` (loops bands, calls `solve_band`) | seam + logic |
| `meds_canopy_aerodynamics` | MO surface layer + Nusselt BL + wind extinction + ground conductance | pure |
| `meds_column_energy` | **4 stores**: soil heat (`soil_energy_flux`+`soil_energy_tendency`+`soil_heat_be_step`); veg (`veg_energy_balance`+`veg_surface_fluxes`); ground (`ground_surface_balance`); CAS (`canopy_air_update`, enthalpy+shv) | mixed |
| `meds_column_hydrology` | soil water (`column_hydrology_flux`+`soil_water_advance`+`soil_be_single_step`+`soil_water_tendency`+`face_and_sink`+`compute_psi_e`+`ground_evaporation`); `intercept_canopy_layer` | mixed |
| `meds_column_co2` | CAS CO₂ (`canopy_air_co2_update`+`column_co2_step`+`aggregate_cohort_co2_fluxes`); `heterotrophic_respiration_flux`/`_damm` | mixed |
| `meds_snow` | mass + energy snow store (`snow_*`) | advance-commit |
| `meds_biophysics_types` | all biophysics derived types + `SOIL_*`/`ENERGY_*`/`HR_*` selectors | data |

### 1.1 Two integration styles coexist (the context)

Biophysics carries **two forms** of the same physics:

- **Advance-and-commit** (`*_flux` / `*_update` / `*_balance`): form a flux, hide it in a
  backward-Euler denominator, mutate the `intent(inout)` store. Consumed by the **operator-split**
  driver `meds_column_dynamics`.
- **Pure tendency / RHS** (`*_tendency`, `surface_derivs`): side-effect-free `f(y)=dy/dt`. Consumed
  by the **IMEX-ARK** driver `meds_column_derivs` + `meds_ark_stepper`.

`canopy_air_update`, `canopy_air_co2_update`, `column_co2_step` are **test-only** today — the split
driver *inlines* the CAS twin updates; the ARK driver forms `d_cas_*` in `surface_derivs`.

## 2. Resolved decisions

| # | question | resolution | rationale |
|---|---|---|---|
| **D1** | integration path | **Keep both** (retain the split integrator; keep ARK) | The split path is the robust fallback and INTEG_ARK still hard-errors on aquifer / Zeng–Decker. Retiring it now costs capability; §7. |
| **D2** | subsystem granularity | one module per surface subsystem; **split only soil** (two large independent solves) | §3. |
| **D3** | ground skin + snow | one module `meds_ground_biophysics` (ground skin balance + snow store) | The two are the bare/snow-covered **modes of one interface**, both producing the soil-top BC + surface→CAS fluxes. |
| **D4** | vegetation | one module `meds_vegetation_biophysics` (leaf/wood energy + canopy interception) | Coupled through the interception film (`leaf_water` → `sigma_w` → leaf latent flux); also disambiguates from `meds_vegetation_dynamics` (slow). |
| **D5** | CAS | one module `meds_cas_biophysics` (enthalpy/shv/CO₂ twins + flux aggregation + heterotrophic Rh) | The three twins are one well-mixed box (shared capacity + `ρ·ustar·temp` conductance). `cas` token matches `cas_state_t`. |
| **D6** | interface module | **create** pure façade `meds_biophysics_interface`; **do not** rename `meds_canopy_radiation` | The plant/core façades are logic-free re-exports; `meds_canopy_radiation` holds real RT logic (§6). |
| **D7** | canopy RT files | **split `meds_optics`**: pure optical-property kernels → shared `meds_optics_lib` (`shared/functions/`); RT assembly + solver + seam → **one** biophysics module `meds_canopy_radiation` (absorbing `meds_twostream_band`); `surface_state_t` → `meds_biophysics_types` | The leaf-angle/scattering kernels are pure semi-empirical constitutive properties (bare scalars, kinds/constants only, already used standalone by setup + tests) — the same properties-in-shared vs assembly-in-domain split as `meds_hydr_lib` vs the soil solver. `solve_band` has a single consumer (the seam). |
| **D8** | `ground_evaporation` home | **stays in `meds_soil_water`** (not the ground module) | It is woven into the Richards substep as the θ-dependent, availability-capped top water-flux BC (`meds_column_hydrology.f90:117,118,139`), and hydrology is the established single authority for `soil_evap`. The ground/energy modules consume it as a frozen input. |

Naming: subsystem process modules take the **`_biophysics`** suffix (`meds_cas_biophysics`,
`meds_vegetation_biophysics`, `meds_ground_biophysics`) — it reads clearly at the `use` site
and separates fast-biophysics from the slow `*_dynamics` drivers. Soil keeps store-based names
(`meds_soil_energy` / `meds_soil_water`) because it is the one split subsystem. RT/aero keep their
established process names.

Scope realized here: **surface-subsystem split (Points 1) + per-store dedup & two-form API (§4) +
façade (Point 3)**. The full physics/numerics *collapse* (moving steppers driver-side) is deferred
(§7).

## 3. Granularity principle — one module per surface subsystem

> **One module per surface subsystem; split a subsystem into multiple files only where it holds two
> large, independently-solved stores.**

- **CAS → 1 file** (`meds_cas_biophysics`): three scalar twins in one box, shared capacity +
  conductance — one coupled solve.
- **Vegetation → 1 file** (`meds_vegetation_biophysics`): leaf/wood energy + interception, coupled
  through the wetted-fraction film — one surface, tightly related.
- **Ground surface → 1 file** (`meds_ground_biophysics`): bare-ground skin + snow, the two
  mutually-exclusive **modes** of one interface.
- **Soil → 2 files** (`meds_soil_energy`, `meds_soil_water`): the sole split — two *large,
  independent* implicit solves (heat diffusion vs Richards) that merely share the negative-z grid +
  the Thomas solver, and are operator-split (weakly coupled: `θ → κ`, lagged). They would merge into
  a single `meds_soil_dynamics` **only** if a *coupled* heat–moisture solve is adopted (§7).

**Corollary — the boundary follows the solve coupling, not the physical label.** A flux lives with
the store whose *solve* it is woven into, even when its physical category suggests elsewhere. Case in
point: `ground_evaporation` is physically a surface latent flux, but it is computed inside the
Richards substep from the evolving top-layer moisture and capped against available water (D8), so it
stays in `meds_soil_water`; the ground/energy modules consume the resulting `soil_evap` as a frozen
input.

## 4. Per-store shape — one physics core, two public forms

Each prognostic store exposes the **same physics** in two forms, over a shared set of private
helpers evaluated at a *supplied* state field. Naming convention:

- `<store>_time_deriv` — **pure** RHS `d/dt` at the current state (explicit, `dt`-independent
  tangent; consumed by the ARK integrator). Renames the old `*_tendency`.
- `<store>_step_implicit` — the **implicit BE advance** over `dt` (mutates the store, re-diagnoses,
  closes the budget; a `dt`-dependent secant — a *step*, not a derivative). Renames the old
  advance-commit `*_flux` / `*_be_single_step`.

Shared private helpers (evaluated at a supplied `T`/`ψ`), so the two forms can never drift:

- `<store>_properties` — state → diagnosed vars + constitutive coefficients (κ, `C_eff` / `K`, `C(ψ)`).
- `<store>_face_cond` — face conductivity (harmonic-mean κ for heat; upstream-weighted `K` for water).
- `<store>_flux_divergence` — fluxes at the supplied field → `d/dt` divergence (+ source/sink).

The symmetry this makes structural (and provably `dt→0`-consistent, vs. today's two hand-kept copies):

```
rate(T) := <store>_flux_divergence(T, ...)                     ! the shared physics
time_deriv     :  dXdt   = rate(X^n)                            ! explicit tangent
step_implicit  :  X^{n+1} solves  X^{n+1} = X^n + dt·rate(X^{n+1})   ! backward Euler
```

### 4.1 Rename map

| today | becomes | module |
|---|---|---|
| `soil_energy_tendency` | `soil_energy_time_deriv` | `meds_soil_energy` |
| `soil_energy_flux` | `soil_energy_step_implicit` | `meds_soil_energy` |
| `soil_heat_be_step` | `soil_heat_be_solve` (private tridiagonal + `thomas_solve`) | `meds_soil_energy` |
| `soil_water_tendency` | `soil_water_time_deriv` | `meds_soil_water` |
| `soil_be_single_step` | `soil_water_step_implicit` | `meds_soil_water` |
| `veg_energy_balance` | `veg_energy_step_implicit` (+ `veg_energy_time_deriv` if a prognostic-leaf RHS is needed) | `meds_vegetation_biophysics` |

### 4.2 The water store's extra layer

Soil water has a **third** level that heat lacks: the adaptive substep controller + surface BCs
(`soil_water_advance`, ponding/runoff/infiltration/DSL in `column_hydrology_flux`) that *orchestrate*
`soil_water_step_implicit`. Under keep-both it stays in `meds_soil_water` as the seam
`column_hydrology_flux`; it is the natural first thing to move driver-side under a future collapse
(§7). Heat's "column advance" is trivial in the MVP (one step), so `soil_energy_step_implicit` *is*
its seam.

## 5. Target module layout

| module | from | holds |
|---|---|---|
| `meds_canopy_radiation` | RT-assembly half of `meds_optics` + `meds_twostream_band` + the seam | `derive_rad_optics` / `blend_cohort_optics` / `ground_optics` + the two-stream solver (`solve_band`/`layer_rt`) + `canopy_radiation` |
| `meds_canopy_aerodynamics` | unchanged | turbulence + conductances |
| `meds_soil_energy` | soil-heat half of `meds_column_energy` | `soil_energy_time_deriv` + `soil_energy_step_implicit` + physics helpers + `soil_heat_be_solve` |
| `meds_soil_water` | `meds_column_hydrology` (soil part) | `soil_water_time_deriv` + `soil_water_step_implicit` + `face_and_sink`/`compute_psi_e`/curve wrappers + `ground_evaporation`; the `column_hydrology_flux` substep/BC seam |
| `meds_vegetation_biophysics` | veg half of `meds_column_energy` + interception from `meds_column_hydrology` | `veg_surface_fluxes` + `veg_energy_step_implicit` (+ `veg_energy_time_deriv`); `intercept_canopy_layer` |
| `meds_ground_biophysics` | `ground_surface_balance` (from `meds_column_energy`) + `meds_snow` | ground skin balance + the snow mass/energy store (`snow_*`) |
| `meds_cas_biophysics` | `meds_column_co2` + CAS half of `meds_column_energy` | CAS enthalpy/shv/CO₂ update + flux aggregation + `heterotrophic_respiration_*` |
| `meds_biophysics_types` | unchanged (not split here) | data |
| `meds_biophysics_interface` | **new** | pure façade re-exporting all of the above |

Removed by consolidation: `meds_column_energy`, `meds_column_hydrology`, `meds_column_co2`,
`meds_snow` (into `meds_ground_biophysics`), `meds_twostream_band` + the RT-assembly half of
`meds_optics` (into `meds_canopy_radiation`).

Moves to **`shared/functions/`**: the pure optical-property kernels (`leaf_class_angle`,
`beta_pdf_kernel`, `beta_lidf`, `beta_params_from_mean`, `leaf_bf`, `gfun_direct`, `scatter_pair` +
the `N_LEAF_CLASS`/`leaf_*_deg` grid) become a shared **`meds_optics_lib`** library (kinds/constants
only), the optical analogue of `meds_hydr_lib`/`meds_therm_lib`. The type `surface_state_t` moves to
`meds_biophysics_types` beside the other `rad_*` types.

Naming convention (in `meds_optics_lib` + `meds_canopy_radiation`): **do not** name the
`rad_pft_optics_t` variable/argument `opt`/`opts` — that token is reserved for the options structs
(`soil_opts_t`, `energy_opts_t`, `co2_opts_t` all bind to `opts`). Use `optics` (or `optic_params`).

## 6. The interface façade (`meds_biophysics_interface`)

A **logic-free re-export** module (the pattern of `meds_plant_interface` / `meds_core_interface`):
`use` + `public` for every seam (`canopy_radiation`, `canopy_aerodynamics`, the soil/veg/CAS/
ground-surface `*_step_implicit` + `*_time_deriv`) and the shared types. Per the plant precedent, it
is a *convenience* for black-box callers (e.g. `meds_fast_loop`); the white-box integrators
(`meds_ark_stepper`, `meds_column_derivs`) keep `use`-ing the kernel modules directly for the
`*_time_deriv` RHS at each stage. `meds_canopy_radiation` stays a seam-with-logic that the façade
re-exports (it is **not** renamed into the façade).

## 7. Retained split solver & the deferred collapse

**Decision: keep the operator-split integrator.** Rationale (from the coupled-vs-split analysis):

- Split keeps each store a self-contained kernel with the right method per store (implicit Thomas
  for stiff heat; Richards/Picard for water), testable in isolation, GPU-friendly, and cheap when
  coupling is weak.
- Its main weakness — splitting error / lagged stiff feedback — is already mitigated in MEDS: the
  **enthalpy** formulation makes freeze/thaw a read-off (no separate stiff phase equation), and
  `SCHEME_PICARD_COUPLED` + IMEX-ARK recover the coupled fixed point when needed, without a joint
  Jacobian.
- A **fully coupled** solve wins only in a specific regime — rapidly advancing freeze fronts in wet
  soil, where ice↔κ↔flow is stiff enough that a lagged split mis-times the front. That regime, not
  aesthetics, is the trigger to (a) merge `meds_soil_energy` + `meds_soil_water` into a coupled
  `meds_soil_dynamics`, and (b) reconsider collapse.

**Deferred collapse (future, gated on ARK feature parity — aquifer, Zeng–Decker):** decompose each
`*_step_implicit` into physics (stays: `*_face_cond`, `C_eff`/`C(ψ)`, `*_flux_divergence`) + numerics
(moves driver-side: tridiagonal assembly + `thomas_solve` + commit + the water substep controller),
then retire the advance-commit forms and the split branch. Because §4 already makes each
`*_step_implicit` a thin wrapper over shared helpers, this later lift is largely mechanical.

## 8. Phased plan (keep-both; each phase behaviour-preserving, builds + tests green ifx & nvfortran)

- **P0 — façade + RT consolidation.** Add `meds_biophysics_interface` re-exporting the current seams +
  types; repoint black-box callers (`meds_fast_loop`) to it (leave the white-box integrators on the
  kernel modules). Split `meds_optics`: pure kernels → shared `meds_optics_lib` (`shared/functions/`),
  `surface_state_t` → `meds_biophysics_types`, RT-assembly + `meds_twostream_band` → the one
  `meds_canopy_radiation` module. Update `use` lists (`meds_fast_loop`, `test_canopy_radiation`). Pure
  addition/relocation; zero behaviour change.
- **P1 — split `meds_column_energy`.** soil-heat kernels → `meds_soil_energy` (with the §4 dedup +
  `soil_energy_time_deriv` / `soil_energy_step_implicit` / `soil_heat_be_solve` renames); veg →
  `meds_vegetation_biophysics`; `ground_surface_balance` → `meds_ground_biophysics`;
  `canopy_air_update` → `meds_cas_biophysics` (created in P3). Golden anchors (`tc_split(54)`,
  per-store budgets) bit-identical.
- **P2 — split `meds_column_hydrology`.** soil-water kernels → `meds_soil_water` (with the §4 dedup +
  `soil_water_time_deriv` / `soil_water_step_implicit` renames); `intercept_canopy_layer` →
  `meds_vegetation_biophysics`. Keep `column_hydrology_flux` intact as the seam (§4.2).
- **P3 — consolidate CAS + ground/snow.** Fold `meds_column_co2` + the CAS half of
  `meds_column_energy` into `meds_cas_biophysics` (move `heterotrophic_respiration_*` with it); fold
  `meds_snow` + `ground_surface_balance` into `meds_ground_biophysics`. Update tests + façade.

Each phase is a relocation/rename + helper factoring — no algorithm change. Verify the standing
discipline: a green ifx run is **not** sufficient; build the nvfortran multicore back end too
(issue #7), and assert the machine-precision budget residuals unchanged.

## 9. Deferred / out of scope

- **Collapse to tendency-only** (move steppers + water substep controller driver-side) — §7, gated
  on ARK aquifer/Zeng–Decker parity.
- **Coupled `meds_soil_dynamics`** (joint heat–moisture implicit solve) — only on a freeze-front
  accuracy need; §7.
- **Split `meds_biophysics_types`** (per-store `energy_*` / `chydro_*` type modules) — a separate,
  larger change.

## 10. Risks & notes

- **Bit-identity.** P0–P3 are pure relocations/renames + a physics dedup that must reproduce the
  existing golden anchors and closed budgets on both compilers.
- **Rename churn** touches the two drivers + the tests (mechanical, module-name + symbol updates);
  the snow tests move with `meds_ground_biophysics`.
- **Docs to update on implementation:** CLAUDE.md biophysics bullet, `src/biophysics/README.md`, and
  the three `docs/science/` pages (module/function names only).
