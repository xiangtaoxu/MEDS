# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What MEDS is

MEDS (**Modular Ecosystem Dynamics Simulator**, also read as **Modern Ecosystem Demography
Simulator** to flag its descent from the ED / ED2 lineage) is a ground-up reimplementation, in
**Fortran 2018**, of the Ecosystem Demography model **ED2** (https://github.com/EDmodel/ED2). It
keeps ED2's science — a size- and age-structured (demographic) terrestrial biosphere model — while
replacing the legacy code structure with modular, testable, standards-conformant Fortran.

The goal is *not* a line-by-line translation. It is to re-express ED2's proven algorithms with
modern language features (modules + submodules, derived-type encapsulation, explicit interfaces,
`allocatable` over `pointer`, parameterized real kinds) so the model is easier to extend, test, and
reason about. When in doubt about *what* a process should compute, the ED2 source is the reference;
when deciding *how* to structure it, follow the modernization guidelines below — do not copy ED2's
structure.

### Reference implementation
The original ED2 tree is checked out as a **sibling directory**: `../ED2`. The core model lives in
`../ED2/ED/src` (driver, dynamics, init, io, memory, mpi, utils). Treat it as read-only reference.
Key entry points to read when porting a feature:
- `../ED2/ED/src/memory/ed_state_vars.F90` — the entire 34k-line state hierarchy (see below).
- `../ED2/ED/src/driver/ed_model.F90` — the main time-integration loop.
- `../ED2/ED/src/driver/edmain.F90` → `ed_driver.F90` — program entry and dispatch.

Background papers: Moorcroft et al. 2001 (Ecol. Monogr.); Medvigy et al. 2009 (JGR);
**Longo et al. 2019 (GMD 12:4309)** — the definitive ED-2.2 technical description (read the
Supporting Information). ED2 wiki: https://github.com/EDmodel/ED2/wiki.

## Repository status

**Implemented: the standalone demographic core** (`src/`, `app/`, `test/`, `CMakeLists.txt`). It runs
cohort & patch dynamics — individual-tree growth, mortality, recruitment, cohort/patch fusion/fission,
and treefall **patch disturbance** — driven by demographic *rates supplied from outside* the engine as
three plain arrays. Size follows the pan-tropical (ED2 `iallom==3`) allometry (`meds_allometry`); each
cohort carries **AGB (carbon)** and **leaf area**, and cohort fusion/fission conserve total AGB. There
is deliberately **no** radiative transfer or full biogeochemistry yet. The standalone Fortran model is the
**carbon-driven** path, driven by `meds_vegetation_dynamics`: the plant carbon seam `get_plant_flux_slow`
turns a stub GPP into per-pool NPP, `wood_carbon` is the prognostic size anchor, the driver's
`compute_slow_derivatives` runs the `wood_carbon→dbh` flip to back out the per-cohort tendency bundle, and
the core engine's `update_cohort_states` applies it (light competition still enters through the stored
`overtopping_lai`). The former **empirical** growth/mortality/recruitment LAWS were moved OUT of Fortran to
the Python example (`examples/example_demography/`, driven through the opt-in `libmeds_c` C-API). See
"Demographic core" below.

Toolchain on this machine (installed, but **off the default PATH** — activate before building):
- **Intel `ifx` 2026** — `source /opt/intel/oneapi/setvars.sh`. Strict-standards CPU compiler; runs the
  full test suite. (`scripts/install_ifx.sh` installs it elsewhere.)
- **NVIDIA `nvfortran` 25.11** (HPC SDK) — add `…/hpc_sdk/Linux_x86_64/25.11/compilers/bin` to PATH.
  This is the parallel/GPU path. The hot kernel (`meds_core`: `update_cohort_states`) carry explicit OpenMP
  `target` regions over plain arrays, so the build picks the device via the NVHPC `-mp` flag
  (`MEDS_GPU=gpu` → `-mp=gpu -gpu=mem:separate`; `MEDS_GPU=multicore` → `-mp`; no flag → serial). All
  three back ends are validated on the RTX 3050 Ti: ifx 7/7, nvfortran multicore 7/7, nvfortran GPU
  7/7 (CPU↔GPU results identical). **Do NOT use `-stdpar=gpu`**:
  it forces the global CUDA-managed allocator, whose deep-copy/finalize of the allocatable-component
  `site` type double-frees on the host. OpenMP `target` + `-gpu=mem:separate` keeps all state in
  normal host memory and moves only the mapped arrays.
- **CMake 4.3.4** (conda, on PATH) — recent enough for both IntelLLVM and NVHPC compiler IDs.
- `gfortran` is not installed (`scripts/install_gfortran.sh` adds it); the build supports it too.
- **Portability trap (nvfortran):** never pass an array-valued *function result* straight into a call
  (`call foo(bar(x))` where `bar` returns an array). nvfortran's whole-program optimizer miscompiles the
  temporary descriptor — **silently wrong values at `-O2`/`-fast`, segfault at `-O0`** — while
  `ifx -stand f18 -check all` tolerates it (only an `arg_temp_created` remark), so a green ifx suite
  hides it. Bind to a named array first (`tmp = bar(x); call foo(tmp)`); this also clears the ifx remark.
  Corollary: a green ifx run is **not** sufficient — build the nvfortran multicore back end on new
  modules too. (See issue #7; found porting `src/biophysics`.)

## Build (CMake)

CMake auto-resolves Fortran `.mod` dependencies — this is the deliberate fix for ED2's "run `make` six
times" hack and per-platform `include.mk` files. Never reintroduce manual object lists or repeated builds.

netCDF output + I/O are **always compiled** (a hard dependency — the `-DMEDS_ENABLE_IO=OFF` netCDF-free
stub build was removed by design), so every build needs the netCDF C library (its CMake config also pulls
HDF5 -> the build enables the C language); `scripts/install_netcdf.sh` installs it (conda or apt) and
prints the `-DCMAKE_PREFIX_PATH` to pass. Every configure must therefore point at the netCDF-C install via
`-DCMAKE_PREFIX_PATH=<prefix>`. Release builds the C-binding IO layer cleanly (Debug `-check all` emits a
harmless arg_temp_created remark per netCDF call), so prefer Release for production output.

```bash
# Default build (Intel ifx). Point at the netCDF-C install via CMAKE_PREFIX_PATH (REQUIRED — netCDF is mandatory).
source /opt/intel/oneapi/setvars.sh
cmake -S . -B build-ifx -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/common
cmake --build build-ifx -j
ctest --test-dir build-ifx --output-on-failure          # 8 tests
LD_LIBRARY_PATH=$HOME/miniforge3/envs/common/lib ./build-ifx/meds_main meds_config_main.toml  # -> <prefix>-D-output.nc (+ -S-<ts>.nc if write_state)

# Debug build (strict -check all) — still needs netCDF (CMAKE_PREFIX_PATH REQUIRED); use Debug for engine checks.
cmake -S . -B build-debug -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/common
cmake --build build-debug -j && LD_LIBRARY_PATH=$HOME/miniforge3/envs/common/lib ctest --test-dir build-debug --output-on-failure

# Multicore CPU / GPU via OpenMP target (NVIDIA nvfortran); netCDF still required (pass CMAKE_PREFIX_PATH)
export PATH=…/hpc_sdk/Linux_x86_64/25.11/compilers/bin:$PATH
cmake -S . -B build-mc  -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_BUILD_TYPE=Release -DMEDS_GPU=multicore \
      -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/common
cmake -S . -B build-gpu -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_BUILD_TYPE=Release -DMEDS_GPU=gpu \
      -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/common
cmake --build build-gpu -j

# Run a single test by regex
ctest --test-dir build-debug -R fusion_cohort --output-on-failure
```

`MEDS_GPU` (`none|multicore|gpu`) only affects NVHPC builds; the `-mp` flags are `PUBLIC` on the
`meds_core` target so `meds_main`/tests inherit the offload compile+link flags. ifx/gfortran ignore
it — the `!$omp` lines are comments without an OpenMP flag, so those builds stay serial. Per-compiler
Debug flags live in the `meds_fortran_flags()` function in `CMakeLists.txt` (ifx `-stand f18 -check all
-fpe0`; gfortran `-std=f2018 -fcheck=all -ffpe-trap`; nvfortran `-Mbounds -Ktrap=fp`).

## Demographic core

### Source layout & libraries
The source tree is organised by process domain over a **state/process wall** (an acyclic library
DAG `shared ← {allometry, plant} ← state ← demography ← aux ← main`, where `plant` links `shared`
ONLY and is orthogonal to `state`); moving a file changes only CMake wiring because Fortran `use` is
by module name and all `.mod`s share one directory. **The 2026-07-04 plant refactor** flattened
`src/plant/` into one ecophysiology library and moved the empirical vital rates into `demography`
(making it self-contained); design: `docs/dev_plans/MEDS_PLANT_ECOPHYSIOLOGY_DESIGN.md`.
- **`src/shared/`** → `libmeds_shared.a` — the foundation, NOT tied to any process: `meds_kinds`
  (precision), `meds_constants`, `meds_pft_params` (the PFT trait table, incl. per-PFT `hgt_max`),
  `meds_time` (calendar + leap-year-aware Gregorian arithmetic), `meds_config` (run config), and
  `meds_temp_response` (Arrhenius / peaked deactivation — promoted here from the leaf module so leaf,
  respiration and any tissue share one code path without a plant→plant library edge). Root of the DAG.
- **`src/allometry/meds_allometry.f90`** → `libmeds_allometry.a` — pan-tropical (`iallom==3`)
  size↔height↔AGB↔leaf-area relations. A shared structural-geometry foundation used by BOTH `state`
  (cohort geometry caching / fusion via `set_cohort_size`/`agb_to_dbh`) and the plant ecophysiology
  library, so it is its OWN library BELOW `state` — it cannot live in `libmeds_plant` without making
  the demographic core depend on all of ecophysiology (see issue #11). `hgt_max` is a per-PFT argument
  to `dbh_to_height`/`agb_to_dbh`.
- **`src/core/`** → `libmeds_core.a` — the CORE ecosystem-structure engine: the cohort/patch STATE
  ontology PLUS the SELF-CONTAINED apply-PRIMITIVES (links `shared` only — `allometry` now lives in
  `shared/functions/`; NO plant-ecophysiology dependency, so the engine mutates state on its own).
  **Five files:** `meds_core_state_types` (the flat site-wide Structure-of-Arrays `cohort_block` + patch
  CSR, the ONE centralized lockstep `cohort_reorder`/`rebuild_csr`/`copy_cohort_slot`/`set_cohort_size`
  machinery, cohort birth `init_cohort`, and the transient tendency bundle `cohort_deriv_block`);
  `meds_core_state_update` (the pure appliers — the OpenMP-target `update_cohort_states` that advances the
  cohort SoA by a supplied per-cohort tendency bundle, the patch-level `update_patch_states`, and the
  `update_overtopping_lai` competition sweep); `meds_core_cohort_fusefiss` (sort + cohort fuse/fission +
  `apply_recruitment`); `meds_core_patch_fusefiss` (sort + patch fuse/fission + treefall
  `apply_patch_disturbance`; depends on the cohort sibling); and `meds_core_interface` (the one-`use`
  public façade). The engine NEVER computes a rate — it APPLIES the tendencies/arrays it is handed; the
  vegetation-dynamics DRIVER computes them. The empirical growth/mortality/recruitment LAWS were moved to
  the Python example, and the carbon vital-rate kernels to `plant` — so `src/core/` no longer hosts a rate
  provider (the former `meds_demography_rates` is deleted). **Naming:** the library and its modules were
  renamed `demography → core` (2026-07-16); design `docs/dev_plans/MEDS_CORE_MODULE_REORG_DESIGN.md`.
- **`src/plant/`** → `libmeds_plant.a` — ONE flat, self-contained plant-PHYSIOLOGY kernel library (links
  `meds_shared` only; NO `site_t`, compiles standalone via `cmake --build … --target meds_plant`).
  Mechanistic per-plant PHYSICAL fluxes only (demographic rate laws live in `demography`, by domain).
  All derived types are consolidated in **`meds_plant_types`**. It holds: **leaf gas
  exchange** — the seam `meds_leaf_physiology%leaf_gas_exchange(env, cfg, ipft, flux)` over
  `meds_leaf_photosynthesis` (FvCB C3 + Collatz C4), `meds_leaf_stomata` (Leuning / Medlyn / Katul),
  `meds_leaf_solver` (bracketed Ci root-find); **hydraulics** (`meds_plant_hydraulics` +
  `meds_hydr_lib`; constitutive PV/Kirchhoff curves in shared, the matrix-exp solver in plant).
  The multi-layer root boundary + soil↔hydraulics coupling are opt-in
  (`[hydraulics].multilayer_roots`, default off ⇒ single root-fraction-weighted BC, bit-identical);
  **hydraulic redistribution (HR) is intentionally NOT enabled** — per-layer root efflux is floored to
  0 in both the plant solver and the soil sink, so uptake is non-negative and conserved. HR is deferred
  to a future version (see `docs/dev_plans/MEDS_MULTILAYER_ROOTS_DESIGN.md`);
  **phenology** (`meds_plant_phenology` + `meds_pheno_engine`); **respiration** (`meds_plant_respiration`);
  and **carbon dynamics** (`meds_plant_carbon_dynamics`). The optional
  Python C-API (`meds_plant_capi.f90`, `-DMEDS_BUILD_PYLIB=ON` → `libmeds_plant_c`, GLOB
  `src/plant/*_capi.f90`) is compiled only into the shared lib and exposes BOTH leaf gas exchange
  (`meds_leaf_solve`) and the phenology kernel (`meds_phenology_step`), through the `meds.plant.leaf` Python
  package + its `meds.plant.pheno` submodule (reproduces Slot & Winter 2017 in
  `examples/example_leaf_gas_exchange/`; the four phenology strategies in `examples/example_phenology/`).
  NOT yet wired into the demographic stepper.
- **`src/biophysics/`** → `libmeds_biophysics.a` — self-contained fast (sub-daily) stateless physical
  kernels, links `shared` only; a sibling stateless-kernel library to `plant`. Modules are grouped **by
  surface subsystem** (one per thermal/chemical store), with a logic-free re-export façade
  **`meds_biophysics_interface`** (the analogue of `meds_plant_interface`) exposing every seam through one
  `use`. **(1) Canopy radiative transfer** (ED2 two-stream `icanrad=2`): the pure optical-property kernels
  (leaf-angle + canopy `scatter_pair` + the `beta_*`/`leaf_bf`/`gfun_direct` family) live in the shared
  **`meds_optics_lib`** (`src/shared/functions/`); the RT assembly (`derive_rad_optics`/
  `blend_cohort_optics`/`ground_optics`), the two-stream solver (`solve_band`/`layer_rt`), and the sealed
  seam `canopy_radiation` all live together in **`meds_canopy_radiation`**.
  **(2) Soil water** (P0/P1/P2; design `docs/dev_plans/MEDS_COLUMN_HYDROLOGY_DESIGN.md`): the 1-D
  soil-water column seam **`meds_soil_water%column_hydrology_flux`** (step `soil_water_step_implicit`,
  explicit sibling `soil_water_time_deriv`, `ground_evaporation`) — implicit backward-Euler Thomas
  Richards with **Celia modified-Picard** or frozen-coefficient linearization, **upstream-weighted K**,
  **retention-integral Zeng–Decker** equilibrium correction, **adaptive step-doubling** substepping,
  conductivity-limited infiltration + ponding, **Dunne (`f_sat`) runoff**, DSL soil evaporation, ψ-limited
  root sink, and a **free-drain / bedrock / SIMTOP-aquifer** bottom BC (diagnosed water table `z_wt`).
  Per-cohort interception (`intercept_canopy_layer`) now lives in `meds_vegetation_biophysics` (below).
  Over the van Genuchten (default) / Campbell
  soil retention curves (`soil_theta_from_psi` / `soil_psi_from_theta` / `soil_hydr_cond_from_theta` /
  `soil_moist_cap_from_psi` + `SOIL_RETENTION_*`), which live in **`meds_hydr_lib`** (`src/shared/
  functions/`) as the soil-water analogue of the tissue PV curves, and the tridiagonal
  **`meds_soil_solver`**; every step closes a machine-precision water budget (`flux%mass_resid`). The
  per-column `soil_params_t` bundle + its `pure` assembler `build_soil_hydr_params` live in
  **`meds_column_state_types`** (beside the prognostic soil columns they describe).
  **(3) Energy balance** (P0/P1/P2a; design `docs/dev_plans/MEDS_ENERGY_BALANCE_DESIGN.md`): four stateless per-store
  kernels solving the land-surface thermal budget, now split **by store** across the surface-subsystem
  modules — leaf/wood (the diagnostic `veg_energy_diagnostic` + prognostic `veg_energy_step_implicit`,
  in **`meds_vegetation_biophysics`**), ground surface (`ground_surface_fluxes`, in
  **`meds_ground_biophysics`**), canopy air space (the two-form box `cas_column_step_implicit` /
  `cas_column_time_deriv`, in **`meds_cas_biophysics`**), and the soil thermal column (`soil_energy_step_implicit` +
  `soil_heat_be_solve`, in **`meds_soil_energy`**, implicit BE-Thomas heat diffusion
  **reusing `meds_soil_solver` + the negative-z geometry**). Prognostic **internal energy / enthalpy (not
  temperature)**, so **freeze/thaw** is a read-off of the shared `meds_therm_lib` inverter (`uext_to_temp`) —
  **P2a turns the plateau on with zero solver change** (`energy_opts_t%phase_change = ENERGY_PHASE_ON`, ice-aware
  `κ_sat(fliq)`/`C_eff(fliq)`), so cooling a wet layer pins `soil_temp` at `t_3ple` while `soil_fliq` absorbs
  `wmass·L_f` (zero-curtain, tested). Closes the forced-temperature seams (`leaf_temp`, `t_ground`, `soil_temp`, RT surface temp).
  Every step closes a machine-precision energy budget. The coupled leaf↔CAS↔ground↔soil fixed point is
  deferred to P3 (kernels take sibling temps as forced inputs). Shared thermal constants live in
  `meds_constants`, and **`meds_therm_lib`** holds both the moist-air thermodynamics (incl.
  `sat_specific_humidity_temp_deriv`, the one shared dqsat/dT slope used by every implicit latent-flux
  linearization) AND the soil thermal constitutive kernels `soil_thermal_cond`/`soil_heat_cap_vol` (the
  thermal twin of the retention curves); the `soil_thermal_params_t` type + `pure` `build_soil_therm_params`
  live in **`meds_column_state_types`**.
  **(4) Canopy aerodynamics** (**`meds_canopy_aerodynamics`**): CLM5 Monin-Obukhov surface layer + ED2
  Nusselt leaf/wood boundary layers + per-cohort wind extinction + CLM ground conductance, emitting the
  `temp1`/`temp2` scalar-transfer factors that set the shared `ustar`-conductance for all three CAS twins.
  **(5) Snow / temporary-surface-water** (the `snow_*` kernels — Niu-Yang cover fraction, accumulation,
  meltwater percolation, snow-surface energy balance + snow-base→soil-top conductance — now live in
  **`meds_ground_biophysics`** alongside the ground-skin balance). **(6) Canopy-air-space CO2 balance**
  (**`meds_cas_biophysics`**; design
  `docs/dev_plans/MEDS_COLUMN_CO2_BALANCE_DESIGN.md`): `can_co2` is the **third prognostic CAS twin**,
  advanced by the shared `cas_column_*` box (the driver assembles the biotic source `Reco − GPP` and
  emits `budg%nee_last`; `heterotrophic_respiration_flux` incl. `HR_DAMM` lives in `meds_soil_biogeochem`)
  — a fast diffusion/venting exchange, so it lives here, NOT in biogeochemistry. Shared derived types live
  in **`meds_biophysics_types`**, which re-exports: the run-config bundles (`soil_opts_t`/`energy_opts_t`/
  `snow_params_t`/`aero_cfg_t` + the `SOIL_*`/`ENERGY_*` selector codes) from **`meds_biophysics_opts`**
  (a low-level `shared/config` leaf, not the `meds_config` aggregator — so the sealed kernels stay
  device-eligible), the soil `*_params_t` types from `meds_column_state_types`, and `SOIL_RETENTION_*`
  from `meds_hydr_lib`. Science pages:
  `docs/science/{canopy_radiation_transfer,canopy_aerodynamics,column_biophysics}.md` (the last with
  per-store pages `{canopy_air_space,soil,snow}_biophysics.md` + `vegetation_energy_dynamics.md`). State-free like RT
  — the per-patch STATE + TOML config + the `psi_soil` and cross-store coupling land at P3 (to couple the
  whole fast loop). The hydrology Neumann→Dirichlet ponded-surface switch and the energy freeze/thaw
  plateau are deferred (P2).
- **`src/biogeochemistry/`** → `libmeds_biogeochemistry.a` — the ecosystem column's **slow soil-carbon /
  nutrient cycle**, a stateless shared-only sibling of `biophysics`. Links `shared` only; kernels are
  `pure`-where-possible / GPU-eligible. (The **fast** canopy-air CO2 exchange `meds_cas_biophysics` + its
  `co2_opts_t`/`damm_params_t`/`HR_*` types are a sub-daily biophysical process and now live in
  `src/biophysics/`; the fast/slow-seam test still re-uses `heterotrophic_respiration_flux` from there.)
  **Slow soil-carbon matrix**
  (P0; design `docs/dev_plans/MEDS_BIOGEOCHEMISTRY_DESIGN.md`): **`meds_soil_biogeochem`** is ED2's CENTURY
  decomposition reorganized as the carbon matrix ODE `dX/dt = B·I + A·ξ·K·X` — a 7-pool `soil_carbon_t`
  (metabolic/structural litter × above/below, microbial, slow, passive; lignin sub-tracer; optional N),
  `assemble_env_scalar`/`assemble_transfer_matrix` (scheme-0 3-active / scheme-5 5-active CENTURY; scalar
  placement `A·(ξ·K·X)`, `a_jj=-1`, `er_j=1-Σaᵢⱼ`), the daily `soil_carbon_step` (forward **EULER** on the
  fast loop's *accumulated* `xi_int`, exact augmented **EXPM** for spin-up; carbon + lignin audits), the
  respired complement `heterotrophic_respiration_matrix`, **SASU** `solve_soil_carbon_steady_state`
  (active-block `{K_j>0}` Gaussian solve — the full 7×7 is singular), and `soil_carbon_diagnostics`
  (capacity/potential, residence time). Shared types + `HR_*`/`DECOMP_*`/`IP_*` selectors in
  **`meds_biogeochem_types`**. Tests: `test/test_column_co2.f90`, `test/test_soil_biogeochem.f90`. State-free
  like the biophysics stores — per-patch `soil_carbon_t` in `state/`, `[soil_carbon]` TOML +
  `ed_params.f90` provenance, netCDF restart, and the demography→litter→Rh driver seam land at P3;
  optional N cycle, DAMM decomposition moisture, and vertically-resolved pools are P1/P2. `src/utils/`
  remains an empty placeholder.
- **`src/forcing/`** → `libmeds_forcing.a` — the home for **prescribed external drivers** (time-varying
  boundary conditions read from a file: meteorology now; disturbance/land-use schedules, prescribed CO₂/N
  later). Links `meds_shared` + the netCDF bindings `meds_netcdf_c` **only** (NOT demography). Design:
  `docs/dev_plans/MEDS_FORCING_DESIGN.md`. **P0 met forcing:** **`meds_forcing_types`** (`met_forcing_t` — the
  instantaneous per-site atmospheric state, a read-only boundary-condition value with defaults summing the
  4 SW streams to 400 W/m²; `met_record_t`; the mutable per-polygon buffer `met_driver_t`),
  **`meds_forcing_kernels`** (`pure`/`elemental`: `interpolate_forcing`, energy-conserving wind, the
  UTC+longitude+EoT apparent-solar transform, the **interval-mean-conserving** shortwave disaggregation
  `cosz_reconstruct_factor` = `1/⟨cosz⟩` **not** `⟨sec z⟩`, the total→(beam/diffuse)×(PAR/NIR) Erbs
  `partition_shortwave`, `dewpoint_to_specific_humidity`/`rh_to_specific_humidity` over `meds_therm_lib`'s
  Bolton esat, `precip_phase`), and **`meds_met_driver`** (`met_open`/`advance`/`instant`/`close` over the
  MEDS multi-grid **`(time,grid)`** forcing NetCDF with per-polygon `grid_index` hyperslab reads + the
  no-file **CONST** backend; **MEDS never gap-fills → a NaN in a required field is a hard error**). The
  `[forcing]`/`[site]` config type `forcing_config_t` + selector codes live in **`src/shared`**
  (`meds_forcing_config`) so `meds_config` carries it without a `shared→forcing` back-edge. `meds_netcdf_c`
  was promoted out of `libmeds_io` into its own always-built target both `io` + `forcing` link. Tested in
  `test/test_met_driver.f90` (ifx + nvfortran). The ERA5-Land prep scripts (`scripts/download_era5land.py`,
  `scripts/prep_era5land_forcing.py`) produce the file it reads. **WIRED into the fast loop** (design §6,
  opt-in): the `[forcing]`/`[site]` TOML block (gated on `forcing.forcing_on`, default false, so configs
  without it are unchanged) loads via `meds_config_io`; `meds_main` opens the reader when `forcing_on` and
  threads it + the step-start time through `advance_one_step` → `fast_dynamics`, which refreshes a
  local `fast_context_t` overlay per sub-step via `apply_met_to_ctx(met_instant(...))` (the diurnal cycle
  lives in the sub-step loop; the constant-forcing path is bit-identical). `test_fast_loop` asserts the
  diurnal signal (night GPP≈0, day GPP>0). **Also wired into the canopy RT (P1 join, design §6.3):** when
  forcing is on, `apply_rt_forcing` (`meds_fast_dynamics`) replaces the LAI-share SW split with the real
  two-stream `meds_canopy_radiation` — maps the met SW streams to `rad_forcing_t`, reverses the
  height-DESCENDING cohort gather order into the two-stream's **BOTTOM(1)→TOP(n)** contract, and
  inverse-scatters per-cohort `abs_sw` (absorbed VIS+NIR, leaf energy) + `abs_par` (VIS ÷ `leaf_absorptance`
  = **incident-equivalent** PAR, since the leaf kernel re-applies absorptance — else it's double-counted,
  ~15% low) + **net** ground SW (`dn_ground−up_ground`, albedo respected; a bare `ncoh=0` patch uses the
  same empty-canopy branch so 1→0 is continuous); PFT-uniform optics live on `fast_context_t`, and the
  per-patch forcing buffers resize per patch. The SAME contract fixed a
  latent bug in `meds_column_dynamics`: `column_fast_step` was calling `canopy_aerodynamics` with the raw
  gather order (inverting its top→bottom wind cascade for multi-cohort patches); the new
  `aero_bottom_to_top` reverses into BOTTOM→TOP and scatters the per-cohort wind/`gb` outputs back (both
  reversals are identity for n≤1, so single-cohort/const paths are unchanged). **Net longwave** is wired:
  per-cohort net leaf LW + net ground LW from the two-stream feed the leaf/ground energy balance, with the
  two-stream's canopy LW emission temperature set to `tcas` (the canopy-air temp) so `abs_lw` = net-LW-at-`tcas`
  matches the diagnostic leaf balance's linearization base (`lw_slope·dtl`) — leaf emission counted once.
  **P2 done:** Weiss–Norman band-specific SW (`SWPART_WEISS_NORMAN`, ED2 `short_bdown_weissnorman`; Erbs stays
  default); multi-year **calendar** recycling + Feb-29 reconciliation (`file_lookup_sec` maps a whole-year
  Jan-1 file to the model calendar year, day-of-year exact, SW reconstruction anchored on the model sun; legacy
  span-wrap kept for non-calendar files); **nearest-grid** match (`grid_match="nearest"`, great-circle argmin —
  the atom of a still-deferred full multi-polygon runtime); wind-height log-profile + hydrostatic elevation
  lapse (opt-in) + the ED2 `reference_height > hgt_max` guard. Deferred: full multi-polygon runtime, LWdown
  synthesis, the phenology daily accumulator.
- **`src/driver/`, `src/init/`** → all part of `libmeds_aux.a` — the top-level utilities that wire the
  process modules together: `meds_stepper` (the thin master stepper / cadence owner, `src/driver`; seed
  of a future all-process **master loop**, ED2-`ed_model` analogue), `meds_vegetation_dynamics` (the
  slow-loop **vegetation-dynamics driver**, ED2-`veg_dynamics_driver` analogue — assembles the carbon NPP
  via the plant seam, computes the per-cohort tendency bundle in `compute_slow_derivatives`, and applies it
  via the core engine's `update_cohort_states`/`update_patch_states`), and
  `meds_init` (`src/init` — the initial-community builders:
  `init_bare_ground`, `add_cohort`, and `init_from_census`). The `meds_aux` target globs `src/driver/*.f90`
  + `src/init/*.f90`, EXCLUDING `src/driver/meds_main.f90`, and links `meds_core` + `meds_plant` +
  `meds_config_io` (the one layer above BOTH the engine and the plant kernels).
- **`src/driver/meds_main.f90`** → the executable `meds_main`, the single entry point (read config →
  build community → run → save output → exit; merged from the former `app/meds_demo` + `app/meds_io_demo`).
  Excluded from `meds_aux` (it is a PROGRAM) and built as the `meds_main` target, linking `meds_aux` +
  the I/O library (real or stub, below).
- **`src/io/`** — netCDF I/O via the netCDF **C** library through `iso_c_binding` (`meds_netcdf_c`
  bindings + `meds_io`, `libmeds_io.a`), built **by default**. `meds_io` has two output streams and a
  restart reader: a DIAGNOSTIC timeseries (`io_create`/`io_write_snapshot`/`io_close` ->
  `<prefix>-D-output.nc`, with derived diagnostics), instantaneous STATE checkpoints
  (`io_write_state` -> `<prefix>-S-<YYYYMMDDHHMMSS>.nc`, prognostic state only — no diagnostics, cached
  geometry re-derived on read), and `io_read_state` (restart). **Diagnostic AGGREGATION subsystem**
  (opt-in `[output]`, design `docs/dev_plans/MEDS_IO_DESIGN.md`, P0 done): a **registry** of scale-suffixed
  variables (`var_desc_t` in `meds_output_types`; the P0 table + §6.1 group/tier/override resolution in
  `meds_output_registry`) feeds per-`(variable,tier)` **integrators** (`meds_output_integrate`:
  `extract_variable` switchboard + `integrate/normalize/reset` + the netCDF-free per-step `output_integrate`
  tick that stages closed periods) that a netCDF **serializer** (`meds_output_stream` per-tier per-time-chunk
  files with CF `cell_methods`/`_FillValue` + `F/D/M/Y` rollover; `meds_output_manager` drains the stage).
  The netCDF-free half (`meds_output_config` in `src/shared` + types/integrate/registry in the
  `meds_io_prep` CMake target) keeps the stepper edge off netCDF (the §2 DAG-hygiene wall); the
  serializer half rides `meds_io`. Wired into `meds_main` (opt-in, coexists with the legacy `[io]` snapshot);
  each active tier integrates raw state independently (exact for MEAN/TMEAN/SUM/MIN/MAX/LAST — chaining +
  variance + the FAST tier are P1). `enable_language(C)`
  fires here so netCDF-C's config can resolve HDF5/Threads. netCDF-Fortran is unavailable for
  ifx/nvfortran here (its `.mod` is gfortran-only); the C API is compiler-agnostic. Point CMake at the
  netCDF-C install with `-DCMAKE_PREFIX_PATH=<prefix>` (here the
  conda env `~/miniforge3/envs/common`). netCDF is a **hard dependency** — `meds_io` is always built (the
  former `-DMEDS_ENABLE_IO=OFF` stub `meds_io_stub.f90` was removed), so every configure must pass
  `-DCMAKE_PREFIX_PATH`. Also here: the always-on TOML config reader
  (`meds_toml` + `meds_config_io`, `libmeds_config_io.a`, no external deps), which also writes the
  per-PFT parameter table to `<output_dir>/<prefix>_pft_parameters.csv` (`write_pft_params_csv`,
  netCDF-free, one row per PFT — a provenance record of the run's trait values).
  - **No hard-coded model parameters.** The source defines only true constants (`meds_constants`,
    plus `meds_allometry`'s coefficients which are `protected` module vars *installed at load* via
    `set_allometry`). Every parameter is REQUIRED from TOML across TWO files — a MAIN file (all non-PFT
    settings, names the PFT file via `[init].pft_config`) and a PFT file (`[pft]` traits + `[camac]`
    mortality coefficients + `[allometry]` coefficients). `load_meds_config` builds a **presence map**
    while reading and `error stop`s listing every missing key; a missing file is also a hard error.
    There is no `build_config`/defaults: derived quantities come from `derive_config` + `derive_pft_rates`
    (overridable via `[options].override_derived` + a `[derived]` block). Tests get a complete config
    from `build_test_config()` in `test/meds_test_support.f90` (the only place "default" values live in
    code). The offloaded appliers take their scalars/arrays as **plain arguments** (they can't read host
    module vars on the device); host allometry functions read them from the module.
- Build order via link deps: `meds_shared → meds_core → meds_aux` (+ optional `meds_io_stream`). The
  core engine compiles standalone (`cmake --build <dir> --target meds_core`).

### Invariants to build on when extending the engine
- **State = flat site-wide Structure-of-Arrays** (`meds_core_state_types`): all cohorts of the whole
  site in one contiguous set of 1-D arrays (`cohort_block`), patch membership as a CSR map
  (`cohort_offset`/`cohort_count` + `owner_patch`). The dominant daily kernels are a single
  unit-stride sweep.
- **OpenMP `target` over plain arrays for the hot kernel** (`meds_core_state_update`:
  `update_cohort_states`) — it takes bare arrays (no `site_t`, no derived types), so the `map` clauses are
  clean and the host keeps all state in normal memory. Keep it arithmetic-only (intrinsics only).
  All restructuring (sort/fuse/split/terminate/recruit) and the tendency COMPUTATION are **host-only**
  (the carbon `compute_slow_derivatives` in the driver calls the branchy allometric flip; the empirical
  one lives in the capi).
- **Tendencies arrive as plain DATA** (the `cohort_deriv_block` bundle) that the vegetation-dynamics
  DRIVER computes and `update_cohort_states` APPLIES (`state += rate·dt`; `nplant *= exp(dln·dt)`,
  floored). The bundle carries every field's time-derivative + the log-space mortality rate, so the engine
  needs no allometry and no empirical/carbon branch — that distinction lives entirely in the driver/capi.
  The bundle is TRANSIENT (per step, not lockstep-reordered). The engine never computes a rate; the
  driver does — **carbon** (the standalone Fortran path): NPP via the plant seam → the `wood_carbon→dbh`
  flip backs out the tendencies. The empirical growth/mortality/recruitment LAWS live in the Python
  example — **growth** = intrinsic (a capped log-linear function of dbh) × competition
  suppression (`exp(growth_lai_slope·overtopping LAI)`) × reproductive-allocation suppression;
  **mortality** = the Camac-2018 additive hazard `mort_gamma + mort_alpha·exp(−mort_beta·growth_avg)`,
  where `growth_avg` is the cohort's tracked simple moving-average growth (window `growth_memory_days`,
  a prognostic per-cohort field); and **recruitment** = baseline `seed_rain_recruits` plus a
  reproduction flux (the carbon diverted from growth to reproduction, per PFT and patch, converted to
  min-size recruits). The engine applies recruitment via `apply_recruitment` (pool + cohort spawn).
  There is
  **no environmental driver** (no temperature) in the current model. A mechanistic module produces the
  same arrays with no engine change. (There is no `rate_provider_t` — data arrays, not a class seam, so
  the engine stays callable from Python.)
- **Conserved invariants** (every fuse/split asserts within 1%, else `error stop`): cohort fusion/split
  conserve **total aboveground biomass (carbon)** and plant number — `dbh` is *re-derived* from the
  conserved per-plant AGB (`agb_to_dbh`), never averaged; patch fusion conserves **site-level plant
  number** via area-fraction rescaling, and patch area always renormalizes to 1. Patch **disturbance**
  conserves area (donors shed a fraction into a new age-0 gap) but intentionally does NOT conserve
  plant number — the killed canopy is the disturbance.
- **One centralized lockstep reorder** (`meds_core_state_types`: `cohort_reorder`/`cohort_compact`/
  `copy_cohort_slot`/`rebuild_csr`/`cohort_ensure_capacity`/`move_alloc_block`, plus `set_cohort_size`
  which fills the cached height/BA/AGB/leaf-area of one slot). When you add a per-cohort field, update
  *these* — the single place that touches every array (the fix for ED2's "forgot to reallocate" class).
  (Patch arrays have no single reorder routine; their permute/pack sites are `sort_patches` and
  `patch_compact` in `meds_core_patch_fusefiss` — update both when adding a per-patch field.)
- **Persistent identity** (`global_id` on `cohort_block` and `patch_block`, monotonic `next_*_id`
  counters on `site_t`): every cohort/patch is stamped at creation via `assign_cohort_id`/
  `assign_patch_id` and carries that id, in lockstep with all other fields, through every
  sort/fusion/compaction; ids are never reused. Fusion keeps the *survivor's* id (the absorbed one
  disappears); a split daughter, a recruit, and a disturbance-gap fragment each get a *fresh* id. This
  lets an external reader (e.g. `post_proc/plot_forest_structure.py`) track one cohort/patch across
  output records until it fuses away or is culled. Creation sites that must stamp ids: `add_cohort`/
  `init_bare_ground` (`meds_init`), `apply_recruitment`, `split_cohorts`, `apply_patch_disturbance`.
- **Order of operations** (`meds_vegetation_dynamics%vegetation_dynamics`, driven by
  `meds_stepper%advance_one_step` in `src/driver/`; cadence flags from the
  caller's calendar): every step `growth → mortality → patch-age`; monthly (`do_cohort_fissfuse`)
  `recruit → cohort fuse/terminate/split → sort`; annual `disturbance` (`do_patch_disturbance`) then
  patch restructuring (`do_patch_fissfuse`: `patch sort/fuse/terminate → cohort consolidate`) — the two
  are INDEPENDENT triggers. Disturbance integrates its yearly rate over the 1-year patch-dynamics
  interval, NOT the per-step `dt_yr`. Same routine serves daily/weekly/monthly steps.
- Cohort fusion/fission keys on **height & LAI** (size anchor is height): similarity is
  `|Δheight| < hgt_max·tol` (the cohort's per-PFT height cap) with geometric tolerance relaxation up to `max_cohort`, and a cohort
  splits (or is blocked from fusing) when its LAI exceeds `cohort_lai_cap`. **Patch fusion** compares a
  cumulative-LAI **light profile** (ED2 scheme, `patch_light_profile`): per-height-layer
  `light = exp(-light_ext·LAI_above)`, fused when the mean AND max layer light difference are within
  `patch_light_tol` / `patch_light_maxdev_factor`, same disturbance class only.
- **PFT contrast is wood density** (`meds_pft_params`). Low ρ ⇒ higher mortality hazard (the
  Camac-2018 `mort_gamma`/`mort_alpha`/`mort_beta` are a power law in ρ, `param_0·(ρ/0.6)^exp`, always
  positive). The growth,
  competition and reproduction parameters (`growth_dbh_*`, `growth_lai_slope`,
  `reproduction_investment_fraction`, `repro_carbon_efficiency`, `seed_rain_recruits`) are per-PFT but
  default to uniform values; the two height thresholds (`min_cohort_height` = recruit/birth size,
  `min_reproduction_height` = reproductive maturity) are shared. Add traits here; derive ρ-dependent
  rates in `derive_pft_rates`.

### Reserved follow-ups (not yet implemented)
A top-level master loop over all processes; an `!$omp target data` region keeping the cohort arrays
device-resident across the daily loop (cuts the per-step map overhead that currently makes the GPU
spin-up migration-bound); a `bind(c)` C-API + shared library for Python (`ctypes`/`cffi`) for the
DEMOGRAPHIC engine — `f2py` will not handle the derived-type/allocatable design, and the data-array
interface (not a Fortran class) is the intended foreign-call layer (the PLANT module already has this:
`src/plant/meds_plant_capi.f90` + `-DMEDS_BUILD_PYLIB=ON` → `libmeds_plant_c`, exposed through
the `meds.plant` Python package (`python/meds/plant`: `meds.plant.leaf` gas exchange + `meds.plant.pheno`
phenology, a clean ctypes-free API installed with `pip install -e python/`), exercised by
`examples/example_leaf_gas_exchange/reproduce_slot2017.py` and `examples/example_phenology/run_phenology.py`);
and **coupling the leaf-physiology module into the demographic
growth** — `src/plant/` exists as a standalone plant-ecophysiology library (FvCB C3 + Collatz
C4, Leuning/Medlyn/Katul stomata, Arrhenius/peaked temperature response), but wiring its assimilation
into the rate provider needs the still-missing canopy radiative transfer (per-cohort absorbed PAR),
leaf energy balance (leaf temperature), a meteorological forcing source, and plant hydraulics
(`psi_leaf`). (Done since the first cut: a single `meds_main` entry point
(`src/driver`); netCDF output — `src/io`, always compiled (a hard dependency) — split into a diagnostic
timeseries and instantaneous STATE checkpoints; initialization from a
cohort census or a restart state file; recruitment moved to the physiology rate layer; temperature
(environmental control) removed; the cumulative-LAI light patch-fusion; the wood-density PFT axis;
treefall disturbance; a real leap-year calendar (`meds_time`) — the run is bounded by `start_time`/
`end_time` dates (not a year count), `meds_main` walks the Gregorian calendar to set the monthly/annual
cadence, both output streams stamp the simulated date, and state files are named `-S-<sim-date>.nc`.)

## Architecture (inherited from ED2)

The science is organized around one central idea: **ecosystem state is a nested demographic
hierarchy**, and the model advances it on **two coupled timescales**. Understanding these two things
is the key to the whole codebase.

### 1. The demographic state hierarchy
State is a tree of nested, ragged arrays. In ED2 this is the monolithic `ed_state_vars.F90`
(`edtype` → `polygontype` → `sitetype` → `patchtype`, ~34k lines). In MEDS this becomes **one focused
module per level** with explicit allocate/deallocate/copy routines (ideally type-bound):

- **grid** — bookkeeping container linking polygons.
- **polygon** — a location sharing one meteorological forcing.
- **site** — shares a soil column and ground hydrology.
- **patch** — shares a disturbance history and stand age (a successional cohort of land).
- **cohort** — plants of one PFT in one height bin (the demographic atom).

Height and age are continuous but **binned**; "identical" means same bin. Bins are defined
dynamically, which is why **fusion/fission** (merging/splitting cohorts and patches to control
resolution) is a first-class operation — see `../ED2/ED/src/utils/fuse_fiss_utils.f90` and
`fusion_fission_coms.f90`. Plant functional types (**PFTs**) parameterize cohort behavior; PFT traits
live in `pft_coms.f90`.

### 2. The two-timescale integration loop
`ed_model.F90` drives the run. Per polygon:

- **Fast loop (sub-hourly, `DTLSM`)** — biophysics: radiation (two-stream), canopy air space,
  photosynthesis + stomatal conductance (Farquhar–Leuning and Katul variants), respiration,
  energy/water/CO₂ budgets, soil/leaf thermodynamics, plant hydraulics. These are stiff ODEs solved
  by a **selectable integrator**: Euler / RK4 / Heun / hybrid (`integration_scheme`; see
  `rk4_*`, `euler_driver.f90`, `heun_driver.f90`, `hybrid_driver.f90`, `bdf2_solver.f90`).
- **Slow loop (daily/monthly/yearly)** — vegetation dynamics (`veg_dynamics_driver`): carbon
  allocation and structural growth, phenology, mortality, reproduction/recruitment, and disturbance
  (treefall, fire, land use / forestry). This is where the demographic structure actually changes,
  followed by fusion/fission to keep cohort/patch counts bounded.

Process families and their ED2 source (under `../ED2/ED/src/dynamics` unless noted):
- Photosynthesis: `farq_leuning.f90`, `farq_katul.f90`, `photosyn_driv.f90`
- Radiation: `twostream_rad.f90`, `multiple_scatter.f90`, `radiate_driver.f90`
- Canopy/turbulence: `canopy_struct_dynamics.f90`
- Phenology: `phenology_driv.f90`, `phenology_aux.f90`
- Growth/allocation: `growth_balive.f90`, `structural_growth.f90`, `stem_resp_driv.f90`
- Mortality / reproduction: `mortality.f90`, `reproduction.f90`
- Disturbance: `disturbance.f90`, `fire.f90`, `forestry.f90`
- Soil C/N + hydrology: `soil_respiration.f90`, `plant_hydro.f90`, `lsm_hyd.f90`, `decomp_coms.f90`

### 3. Parameters, I/O, and parallelism
- **Parameters** — ED2 holds tunable parameters/constants in many global `*_coms` modules
  (`consts_coms.F90`, `pft_coms.f90`, `soil_coms.F90`, `physiology_coms.f90`, …). In MEDS these
  become **derived-type configuration objects passed explicitly** (see modernization guidelines) so
  state is not hidden global mutable data.
- **I/O** — ED2 reads an `ED2IN` Fortran namelist plus optional XML config, and writes HDF5 history
  (restart) and analysis output (`../ED2/ED/src/io`). Keep I/O behind a dedicated layer; isolate the
  HDF5/netCDF dependency so the science modules never touch file formats directly.
- **Parallelism** — MPI distributes polygons across ranks (master/slave, guarded by `-DRAMS_MPI`),
  with OpenMP threading inside a rank. Both are optional and must build cleanly when disabled.

## Modernization guidelines (ED2 → MEDS)

These are the rules that make MEDS "modern and modular." They override ED2's structure wherever they
conflict.

- **One responsibility per module; one module per file.** Never reproduce a 34k-line file. Put type
  definitions in `*_types` modules and heavy procedures in **submodules** to cut compile coupling.
- **`implicit none` everywhere; explicit interfaces via modules.** Drop ED2's `USE_INTERF` flag that
  *bypassed* interface checking — the compiler should check every call.
- **Parameterized kinds.** A single `shared`/`kinds` module defines `sp`, `dp` (and integer kinds);
  use `real(dp)` etc. Never write bare `real*8`, `kind=8`, or `1.0d0` literals scattered in code.
- **No hidden global mutable state.** Replace mutable `*_coms` globals with derived-type config/state
  passed as arguments. This is what makes routines unit-testable, thread-safe, and reentrant.
- **`allocatable` over `pointer`** for ownership; use `move_alloc` and automatic deallocation instead
  of manual bookkeeping. Reserve `pointer` for genuine linked structures.
- **`pure`/`elemental`** for leaf math kernels (photosynthesis, allometry, thermodynamics) — they are
  the easiest to test and the most reusable.
- **Errors via `error stop` / status codes**, not ED2's `fatal_error` halt pattern, so failures are
  catchable in tests.
- **≤132 columns, free-form, lowercase keywords.** ED2 is intentionally 132-column compliant; keep
  that — do not rely on `-ffree-line-length-none` to mask over-long lines.
- **Readability over succinctness when performance is equal**, especially for ecologically/physically
  meaningful names. Spell out `site`, `cohort`, `patch`, `dbh_critical`, `hgt_max` (the site-level
  type is `site_t`, the instance `site`; containers are `site%cohort` / `site%patch`); keep terse
  *loop indices* (`i`, `ip`, `pf`/`ipft`) and the established domain tokens `dbh`, `nplant`, `pft`.
  Inside `associate`, alias to the full word (`cohort => site%cohort`, `patch => site%patch`,
  `pft => cfg%pft`), not single letters.
- **Test as you port.** Each ported kernel gets a unit test (CTest target); validate whole-model
  behavior against ED2 outputs (the analog of ED2's `EDTS` regression suite). A port is not done
  until it reproduces the reference within tolerance.

## Git conventions

Track Fortran source (`*.f90 *.F90`), CMake files (`CMakeLists.txt`, `*.cmake`), namelists, and docs.
Build artifacts (`*.o *.mod *.smod *.a`, the `build*/` trees) are ignored — see `.gitignore`. Line
endings are normalized to LF via `.gitattributes`. Keep generated HDF5/netCDF output and large input
datasets out of the repo.

## Writing docs — GitHub Markdown math

The science pages (`docs/science/*.md`) render on GitHub, whose Markdown sanitizer runs *inside*
`$…$` / `$$…$$` **before** MathJax. It unescapes backslash-punctuation and applies emphasis to the
LaTeX, so `\,`→`,`, `\;`→`;`, `\!`→`!`, `\\[2pt]` breaks every `\begin{cases}` row, and paired `_`
get eaten by italics. **Escaping does not help — GitHub is the thing that unescapes.** The only fix
is to *protect* the math from Markdown, using GitHub's two protected-math syntaxes:

- **Display equations → ` ```math ` fenced code blocks, never `$$…$$`.** A code fence isn't
  Markdown-processed, so `\,` spacing, `\\` cases-row separators, `\begin{cases}`, and `\_` all reach
  MathJax intact. (Fences need a blank line before/after; put the number inline as `\qquad(1)`.)
- **Inline math with any Markdown-conflicting char → `` $`…`$ `` (dollar-backtick), not `$…$`.**
  "Conflicting" = contains `\,` `\;` `\!` `\_` `\{` `\}` `\*`, **two or more** brace-subscripts
  (`X_{a}…Y_{b}`, which pair into emphasis), or a `}_` subscript-after-brace.
- **Simple inline stays plain `$…$`** — intraword subscripts like `$C_i$`, `$g_s$`, `$A_g$`,
  `$\psi_{50}$` (a single brace-subscript) are never mangled.
- **Blocked regardless of protection:** `\operatorname` (use `\mathrm{…}`) and `\tag{}` (renders as a
  vertical jumble; number manually with `\qquad(1)`).
- **Prose underscores** outside code spans that could pair into italics (e.g. `β_stomata / β_nonstomata`)
  must be escaped (`β\_stomata`) or wrapped in backticks.

Inside a protected block/span the LaTeX is honored verbatim, so write it cleanly: `\Gamma^*` (not
`\Gamma^\*`), real config-key underscores (`\text{vessel\_curl}`), proper `\,` spacing. There is **no
local MathJax preview** — validate structurally (no `$$` left, `` ```math `` fences balanced, every
conflicting construct inside a fence or `` $`…`$ `` span) and eyeball the rendered file on a branch/PR.
