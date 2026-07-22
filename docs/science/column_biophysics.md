# Column biophysics

The integrative fast-loop (sub-daily) surface model: how MEDS advances the coupled internal-energy,
water, and CO₂ budgets of the whole soil–vegetation–atmosphere column each `dt_fast`. This page ties
together the per-store kernels — each documented on its own page — and describes the two integrators
(`meds_column_dynamics`) that weave the stores together:

- [canopy_air_space_biophysics](canopy_air_space_biophysics.md) — the CAS coupling reservoir and its
  three prognostic twins (enthalpy, humidity, CO₂).
- [soil_biophysics](soil_biophysics.md) — the coupled soil water (Richards) and soil thermal
  (internal-energy heat) columns, canopy interception, ground evaporation, heterotrophic respiration.
- [vegetation_energy_dynamics](vegetation_energy_dynamics.md) — the leaf/wood tissue energy balance
  (diagnostic and prognostic modes).
- [snow_biophysics](snow_biophysics.md) — the opt-in snow / temporary-surface-water store.
- [canopy_aerodynamics](canopy_aerodynamics.md) — the Monin-Obukhov surface layer and boundary-layer
  conductances that set the shared $`u_*`$ atm↔CAS transfer.
- [canopy_radiation_transfer](canopy_radiation_transfer.md) — the two-stream absorbed SW / net LW that
  drive the surface energy stores.
- [plant_hydraulics](plant_hydraulics.md) — the soil→plant→stomata water pathway coupled through
  `psi_soil` and the realized transpiration.

---

## Two cross-cutting principles

1. **Prognostic internal energy, not temperature.** Every thermal store carries specific enthalpy or
   internal energy as its prognostic variable and diagnoses temperature (and liquid fraction) by
   inverting the thermodynamics (`meds_therm_lib%uext_to_temp`). Freeze/thaw is then a *read-off* of the
   inverter — cooling a wet layer pins its temperature at the triple point while the internal energy
   keeps falling and the liquid fraction absorbs $`w\,L_f`$. Phase change comes for free, with zero
   solver change.
2. **Closed-budget discipline.** Every store returns a residual that is the change in storage minus the
   time-integral of its boundary fluxes, algebraically zero by construction and asserted to round-off in
   Debug builds. A whole-column ledger (`budg%whole_energy`/`whole_water`) catches cross-seam leaks the
   per-kernel budgets would miss.

The **canopy air space (CAS)** is the shared coupling reservoir. It carries three implicit prognostic
twins — enthalpy, specific humidity, CO₂ mixing ratio — each advanced from its surface sources and an
implicit atmospheric-exchange term that shares the $`u_*\cdot\mathrm{temp1/temp2}`$ conductance from the
aerodynamics kernel. See [canopy_air_space_biophysics](canopy_air_space_biophysics.md).

---

## Weaving the stores: operator split vs IMEX-ARK

`column_fast_step` (`meds_column_dynamics`) advances one patch one `dt_fast`. It offers two integrators
(`cfg%time_integrator`).

### Operator split (default) with leaf↔CAS Picard coupling

Per sub-step the driver runs a **pre-pass** (leaf gas exchange → GPP/$g_s$/$R_d$; stem+root respiration;
CAS capacities; aerodynamics), then a sequence: aerodynamics → {diagnostic leaf balance, ground balance}
→ soil-water column (infiltration / DSL evaporation / root uptake / drainage → `src_frac`, the
supply-limited fraction of transpiration demand) → CAS three-twin update (implicit atm; see
[canopy_air_space_biophysics](canopy_air_space_biophysics.md)) → soil thermal column (implicit BE,
reading the just-updated moisture). Plant hydraulics then advances each cohort's node potentials from the
realized transpiration and the root-weighted `psi_soil`, feeding *next* step's leaf gas exchange (the
soil→plant→stomata drought feedback, lagged one `dt_fast`).

Under `integration_scheme="picard"` (`SCHEME_PICARD_COUPLED`) the sequence
{leaf energy → soil-water/`src_frac` → CAS twins → soil thermal} is wrapped in an **outer Picard fixed
point**. State$`^n`$ is snapshotted once; each pass re-solves the *same* backward-Euler steps from the
snapshot, under-relaxing the CAS seed (`picard_relax`, ~0.5, to tame the slope-$-1$ oscillation) while
the committed BE solutions stay exact. `niter=1` reproduces the pure operator split **bit-identically**
(`SCHEME_SPLIT_SEQUENTIAL`), so the coupled path is a strict generalization. Convergence tests the
inter-iterate CAS+leaf+ground temperature and CAS humidity; a non-converged run clamps to the last
iterate (never a partial state).

### IMEX-ARK (opt-in)

`column_fast_step_ark` drives the coupled column with an additive Runge-Kutta stepper
(`meds_ark_stepper`, `docs/dev_plans/MEDS_IMEX_ARK_DESIGN.md`), which needs a side-effect-free
$`f(y)=dy/dt`$ for the whole column — supplied by `meds_fast_time_derivs` (`column_derivs`). Each
reservoir's tendency is the explicit RHS whose BE advance reproduces its split kernel as
$`\Delta t\to0`$ (soil heat/water, via `soil_energy_time_deriv`/`soil_water_time_deriv`) or exactly (the CAS
implicit-in-atm twins; the plant-hydraulics $2\times2$ matrix-exponential). The surface hydrology BCs
(`q_top`, `psi_e`, `soil_psi_root`) and the leaf-gas-exchange pre-pass are the **frozen, explicit** part
of the additive split, held constant across the macro-step; soil water is currently operator-split out
(the robust ponding/runoff scratch solve commits `theta1`, a lagged split). An ARK **conservation
ledger** accumulates the $b$-weighted boundary-flux amounts so the same budgets the split closes also
close on the ARK path.

Whichever integrator runs, **every store closes its residual** and the whole-column ledger closes to
machine precision — the invariant that makes the coupled surface model trustworthy.

---

## Prognostic state (per patch, `patch_biophys_t` / `meds_column_state_types`)

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Canopy air space | `cas_state_t` | `can_enthalpy`, `can_shv`, `can_co2` | `can_temp` |
| Soil water | `soil_column_t` | `theta(:)`, `w_surface`, `w_aquifer` | `z_wt`, `psi_soil` |
| Soil thermal | `soil_energy_column_t` | `soil_energy(:)` [J m⁻³] | `soil_temp`, `soil_fliq` |
| Snow / surface water | `snow_column_t` | `swe`, `snow_energy` [J m⁻²], `snow_depth` | `snow_temp`, `snow_fliq` |
| Leaf / wood tissue | `patch_biophys_t` | `leaf_temp`, `wood_temp` (internal energy) | temp / fliq read-offs |
| Plant hydraulics | `patch_biophys_t%psi` | node water potentials (leaf/wood/root) | — |

Area-weighted `blend_*` mixers conserve these intensively when patches fuse or a disturbance gap is
carved (temp/fliq are re-diagnosed, never blended).

## References
- ED2 `../ED2/ED/src/dynamics/`: `rk4_derivs.f90`, `canopy_struct_dynamics.f90`, `soil_respiration.f90`,
  `lsm_hyd.f90` — the fast-loop stores in the reference model.
- Kennedy & Carpenter (2003) — additive Runge-Kutta (IMEX) methods.
- Per-store references (retention curves, DSL evaporation, snow-cover fraction) are on the domain pages
  linked above.
- Design docs: `MEDS_COLUMN_DYNAMICS_DESIGN.md`, `MEDS_ENERGY_BALANCE_DESIGN.md`,
  `MEDS_COLUMN_HYDROLOGY_DESIGN.md`, `MEDS_COLUMN_CO2_BALANCE_DESIGN.md`, `MEDS_IMEX_ARK_DESIGN.md`,
  `MEDS_SNOW_DESIGN.md`, `MEDS_LEAF_WOOD_ENERGY_DESIGN.md`, `MEDS_P3_COUPLED_SURFACE_DESIGN.md`.

## Code map (whole column)

| Concept | Routine |
|---|---|
| CAS enthalpy + vapour + CO₂ twins | `meds_cas_biophysics`: `cas_column_step_implicit`, `cas_column_time_deriv` |
| soil water (implicit Richards) | `meds_soil_water`: `column_hydrology_flux`, `soil_water_step_implicit`, `soil_water_advance` |
| canopy interception | `meds_vegetation_biophysics`: `intercept_canopy_layer` |
| soil thermal (implicit BE heat) | `meds_soil_energy`: `soil_energy_step_implicit`, `soil_heat_be_solve` |
| leaf/wood energy | `meds_vegetation_biophysics`: `veg_energy_diagnostic` (shared diagnostic solve), `veg_energy_step_implicit` (prognostic store) |
| ground skin fluxes | `meds_ground_biophysics`: `ground_surface_fluxes` |
| soil heterotrophic Rh | `meds_soil_biogeochem`: `heterotrophic_respiration_flux`, `heterotrophic_respiration_damm` |
| snow energy / base conductance | `meds_ground_biophysics`: `snow_energy_step`, `snow_base_conductance` |
| snow mass / cover / melt | `meds_ground_biophysics`: `snow_accumulate`, `snow_cover_fraction`, `snow_drain_meltwater` |
| operator-split + Picard step | `meds_column_dynamics`: `column_fast_step` |
| IMEX-ARK step | `meds_column_dynamics`: `column_fast_step_ark`; `meds_ark_stepper` |
| whole-column tendency RHS | `meds_fast_time_derivs`: `column_derivs`, `surface_derivs` |
| explicit tendency siblings | `soil_energy_time_deriv`, `soil_water_time_deriv`, `plant_water_tendency` |
| prognostic column types | `meds_column_state_types`: `cas_state_t`, `soil_column_t`, `soil_energy_column_t`, `snow_column_t` |
