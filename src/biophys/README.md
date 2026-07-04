# biophys

**Fast, sub-daily, mostly-stateless physical flux calculators.** Device-eligible and netCDF-free
(forcing enters via passed-in value types, never a direct `use netcdf`). Sealed and orthogonal to the
demographic engine: the modules here link `meds_shared` only (no `site_t`), so they compile and test
standalone, exactly like `src/plant/leaf/`.

## Canopy radiative transfer (`meds_canopy_radiation`)

A faithful reimplementation of ED2's two-stream canopy RT (the `icanrad=2` solver in
`../ED2/ED/src/dynamics/twostream_rad.f90`; Liou 2002, Longo et al. 2019), modernized per
`archive/radiative_transfer_design.md`:

- **One unified multi-band solver** (`meds_twostream_band`) for the default bands VIS / NIR / LW,
  working in absolute W m⁻². Every band carries a thermal-emission source, identically zero for VIS/NIR.
- **SCOPE / 4SAIL leaf-angle scattering** (`meds_leaf_angle`, `meds_canopy_optics`): a two-parameter
  **Beta** leaf-inclination distribution → `bf = <cos^2(theta)>` and the exact Ross `G(mu)`; per band
  `omega = rho+tau`, diffuse backscatter `beta = 0.5(1+bf(rho-tau)/omega)`, beam upscatter
  `beta0 = 0.5(1+(bf/k)(rho-tau)/omega)`. Replaces ED2/CLM's `phi1/phi2/mu_bar` + `(1+chi)^2`
  backscatter, which breaks for near-horizontal leaves.
- **O(N) adding-method solve** (no external BLAS/LAPACK): the block-tridiagonal two-stream solved by a
  bottom-up / top-down flux recursion, energy-conserving by construction.
- **Leaf + wood** with clumping-corrected `elai`, `ewai` (WAI defaults to `0.1*LAI`), and a leaf/wood
  absorption split consistent with the solver's own weighting.
- **Surface reflectance** (`meds_surface_optics`) is a bare-soil placeholder (per-band albedo + thermal
  emission); it grows a full soil/snow/water model when soil state exists.

The public seam is `meds_canopy_radiation%canopy_radiation` -- the RT analogue of
`meds_leaf_physiology%leaf_gas_exchange`. It returns per-cohort absorbed radiation (leaf & wood, per
band) plus patch albedo and below-canopy fluxes; the absorbed leaf PAR (`RAD_VIS`) is the field the
leaf-physiology module will consume. Exercised by `test/test_canopy_radiation.f90`.

**Not yet wired into the stepper** (no `site_t` orchestration, meteorological forcing, or leaf/energy
coupling) -- a standalone calculator, like the leaf module was at its first cut. Reserved follow-ups:
per-patch `site_t` orchestration + stepper coupling; TOML config wiring of the optical traits; wood-area
allometry; a real soil-optics model; GPU offload of the solve (interleaved layout). Also here (future):
leaf/canopy energy balance and soil hydrology.
