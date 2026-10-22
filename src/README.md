# `src/` — the MEDS source tree

**37.2 k lines, 105 modules, 21 CMake libraries.** This page is the map: the layout, the four rules
that decide where a new file goes, and the library graph the build enforces. It is the page to
read before adding code, and the one to hand someone who asks how MEDS is organized.

The science is in [`docs/science/`](../docs/science/); the design record is in
[`docs/dev_plans/`](../docs/dev_plans/).

## The shape, in one sentence

MEDS is an **operator-split fast/slow model**, so the tree splits **by timescale** first and by
**domain** second, over a **state layer** that both halves share.

```
src/                                    37.2 k lines · 106 modules · 21 CMake libraries
│
│  ── FOUNDATION ─────────────────────── no model state, no process: libmeds_share
├── share/         uses nothing outside itself; every layer may use it: kinds,
│   │              constants, numerics (clamps, means, the conservation test),
│   │              linear algebra (Thomas sweep, dense LU, arrowhead solve,
│   │              matrix exponential), calendar time and the sun's
│   │              position, budget checks                                   1 249
│   └── functions/ stateless constitutive laws, one library per domain:
│                  allometry, plant hydraulics, soil hydraulics,
│                  thermodynamics (moist air and its lapse, with the
│                  temperature response of enzyme rates), canopy optics,
│                  the radiation above the canopy                            2 129
│
│  ── CONFIGURATION ──────────────────── the run's inputs: libmeds_config
├── config/        meds_config_t (one field per key) and the routines that
│   │              derive what the parameters imply (the PFT table's arrays;
│   │              a site's soil profile, meds_soil_derive), the TOML reader,
│   │              the loader, and the checks of the keys and their values   3 214
│   └── params/    the run's parameters, one leaf per TOML block: the types
│                  and codes a kernel imports without meds_config (the PFT
│                  table, the soil, forcing, output, the fast loop, snow and
│                  aerodynamics, leaf, biogeochemistry, sites, region, maps)   832
│
│  ── STATE, one file per level of the hierarchy, each in sections ────────────
├── state/         meds_cohort_state_vars: the flat cohort structure-of-arrays
│                  and its lockstep; meds_patch_state_vars: one patch's stores
│                  (canopy air, snow, soil water and heat, soil carbon), their
│                  blends, and the patch block with its CSR map;
│                  meds_site_state_vars: the site block;
│                  meds_polygon_state_vars: polygon_t; meds_diag_vars: the
│                  diagnostic accumulators; meds_state_defaults: the values
│                  a field holds before it has a real one                    2 756
│
│  ── PROCESSES, timescale first ──────────────────────────────────────────────
├── fast_dynamics/       sub-daily, one dt_fast per step                      9 311
│   ├── canopy/    the medium: two-stream radiation, aerodynamics, the
│   │              canopy-air box                                            1 028
│   ├── plant/     the organisms: leaf gas exchange, hydraulics, maintenance
│   │              respiration, tissue energy                                1 063
│   ├── soil/      the ground: soil water, soil energy, the ground skin,
│   │              the snow store                                              749
│   ├── integrator/ the integrator machinery: the state vector, the frozen
│   │              work record, ARK and RK45, error control, the pre-pass    5 034
│   └── driver/    the fast loop's setup from the configuration, and the
│                  walk of one slow step in dt_fast sub-steps over the
│                  patches                                                   1 437
│
├── slow_dynamics/       daily to annual                                      4 826
│   ├── plant/     phenology, carbon allocation, trait plasticity               601
│   ├── soil/      the CENTURY soil-carbon matrix, the litter partition         826
│   ├── demography/ the vital-rate LAWS, and the operators that apply them:
│   │              state update, cohort and patch fuse-fission, recruitment,
│   │              treefall disturbance                                      1 382
│   └── driver/    the slow coordinator, the vegetation and biogeochemistry
│                  drivers, and the slow conservation ledger                 2 017
│
│  ── EDGES ────────────────────────────────────────────────────────────────────
├── io/            every data file a run reads or writes                     10 013
│   ├── (root)     the netCDF C bindings, and the helpers both writers share    615
│   ├── init/      the initial state, by [init].init_mode: bare ground, a
│   │              cohort census, or a restart file, which the checkpoints
│   │              write; init_polygon starts a polygon from one and seeds
│   │              the fast loop's stores (init_fast_state)                   1 693
│   ├── met/       the atmosphere above the canopy: the met reader (a MEDS
│   │              forcing file, or the ED_ERA5land archive) and its sample at
│   │              a site and at a patch's canopy-air top (meds_met), the grid
│   │              its cells are on (meds_met_grid) and the CO2 series
│   │              (meds_met_co2)                                             2 905
│   ├── maps/      static per-cell maps on the forcing grid: the one place
│   │              a map is read (meds_maps)                                    350
│   └── output/    the diagnostic wall: derive → capture → reduce →
│                  integrate → serialize                                      4 450
├── c_api/         bind(c) shims → one libmeds.so: leaf, phenology,
│                  demography, and the full coupled run                        938
└── main/          meds_stepper (one coupled slow step, and the boundary
                   restructuring), meds_polygon (one polygon's run
                   machinery and step), meds_region (a run of one polygon
                   or many: open / step / finalize), meds_main
                   (the PROGRAM)                                             1 403
```

**`share/` is one library, and it admits only what depends on nothing outside it.** Every module
there `use`s other `share/` modules and compiler intrinsics and nothing else: no model state, no
configuration, no external library. That is what lets every other layer, and the tests that link
one narrow library, use it freely. A file that needs `meds_config`, a state type or netCDF is not
foundation however small or widely used it is, and goes to its caller's layer instead; `config/`
is a library of its own because configuration is not foundation. The test is what keeps `share/`
from becoming the place a file goes when no other place fits.

## Four rules, in the order you apply them

**1. Kernel folders never see `polygon_t`.** `fast_dynamics/{canopy,plant,soil}` and
`slow_dynamics/{plant,soil}` take from the state layer only the per-patch store types and their
blends, through `use meds_patch_state_vars, only : ...`, and link only `state` and `config`. That is
what keeps each kernel library building standalone, keeps the kernels eligible for OpenMP `target`
offload (a portability property worth holding, though offload is not currently a speedup — see
[`dev_plans/archive/MEDS_GPU_EVALUATION.md`](../docs/dev_plans/archive/MEDS_GPU_EVALUATION.md)), and lets a
test link one kernel library alone. It is checked on every ctest run, by `kernels_no_site_state`
(`scripts/lint/check_kernel_state.py`), because the library boundary cannot hold it: every module
file goes to one shared directory, so any file can use any built module. **If a new routine needs
`polygon_t`, it is driver code**, and it belongs in a `driver/` folder.

**2. A kernel goes where its caller's timescale is**, and within that, in its domain folder. A
kernel called from *both* tiers is a documented seam, never a folder of its own. There is exactly
one: `heterotrophic_respiration_matrix`, so that the sub-daily respiration debits the same CENTURY
matrix the daily step does. Co-locating it with the daily step is what closes the seam gap to
machine precision.

**3. A derived type lives with whoever mutates it.** If two subsystems mutate it, it is state, and
belongs in the `state/` file of the level that holds it: a patch's stores in
`meds_patch_state_vars`, a cohort's fields in `meds_cohort_state_vars`. **Parameters are not state**: they are derived once and never
integrated, so they live with the configuration they are derived from: their types in
`config/params/` (`soil_params_t` and `soil_profile_t` in `meds_soil_params`, `pft_table_t` in
`meds_pft_params`), and the routines that derive them in `config/` (a site's soil profile in
`meds_soil_derive`, the PFT table's leaf capacities in `meds_config`). A kernel's own argument types are not parameters in this
sense: `*_args_t` (`leaf_photo_args_t`, `pheno_args_t`, …) carry the values a kernel reads, a driver
fills them from the parameters, and they live with the kernel.

**4. Laws and operators do not touch.** `demography/` holds both the vital-rate laws
(`meds_demography_rates`) and the operators that apply them (`update_*`, `*_fusefiss`,
`terminate_*`, `apply_recruitment`, `apply_patch_disturbance`) — and the operators take rate
**arrays** as arguments and never import the laws. The slow driver is the one place a rate meets
its application. The Python `apply_rates` path, which feeds externally computed rates through the
same operators, is the standing test that the separation holds.

## The library graph

Acyclic by construction, and the build enforces it. Every arrow is a `target_link_libraries` edge
in [`CMakeLists.txt`](../CMakeLists.txt).

```
  share  ──→ config ──→ state ──────────────────────→ demography ──→ output
     │                    │                              ↑             │
     │                    ├──→ fast_kernels ──────┐      │             │
     │                    └──→ slow_kernels ──────┼──────┘             │
     │                                            │                    │
  netcdf_c, config ──→ met ─────────────────────→ fast  ←──────────────┘
                                                  slow
                                                  init  ←── slow, netcdf_c
                                                    └──→ stepper ──→ model
                                                                     (INTERFACE)
  met, config ──→ maps ──────────────────────────┐
  output_stream ─────────────────────────────────┤
  model ─────────────────────────────────────────┴──→ driver ──→ meds_main
                                                               └→ meds_py (libmeds.so)
```

One edge is worth knowing because it is easy to get wrong:

- **`output` carries no netCDF.** The netCDF-free half of the output subsystem (the registry, the
  file set's setup, the extract, the reductions, the temporal integrators) is a separate target
  from the writer, `output_stream`, so the stepper can tick the integrators without pulling a C
  dependency. It globs `io/output/` less the writer. The fast loop links it because it stages the
  FAST tier's sub-step samples in a polygon's output buffers.

## Where does a new file go?

| You are adding | It goes in |
|---|---|
| a soil, canopy or plant process | `fast_dynamics/<domain>/` or `slow_dynamics/<domain>/`, by the **cadence of its caller** |
| a per-patch store the kernels and the patch block share | `state/meds_patch_state_vars.f90`: its type, its blend, and its section of `patch_block` |
| a per-cohort field | `state/meds_cohort_state_vars.f90`: its section of `cohort_block` and of each lockstep list — and every creation site, the fusion policy, and, if it persists, the restart file (`io/init/meds_init_restart.f90`: three adjacent lines in its cohort section) |
| a per-patch field | `state/meds_patch_state_vars.f90`: its section of `patch_block` and of each patch list — and `fuse_2_patches`, the new gap in `apply_patch_disturbance`, and, if it persists, the restart file's patch section |
| a stateless constitutive curve: a physical or biological law with parameters | `share/functions/` |
| a model-independent tool every layer needs (numerics, calendar, checks) | `share/`, if it uses nothing outside `share/` |
| a TOML block | a leaf in `config/params/`, so `meds_config` can carry it with no back-edge; its keys in `meds_config_main.toml`, its checks in `config/meds_config_check.f90` |
| a prescribed time-varying input (weather, CO2) | `io/met/`: a reader named after its `[forcing].format` value; a gridded product, a `met_grid_t` for `meds_met_grid` |
| a way to start a run (`[init].init_mode`) | `io/init/`: a `meds_init_<start>` module, dispatched by `init_polygon` |
| a static per-cell input map (soil, a cell mask) | `io/maps/meds_maps.f90` |
| an output variable | one `add_variable` line in `io/output/meds_output_registry.f90`; a new field also needs its source id and accessor case in `io/output/meds_output_extract.f90` |
| a Python entry point | `c_api/`, plus its mandatory ctest target |
| a top-level orchestration step | `main/meds_stepper.f90` or a `driver/` folder |

## Invariants the tree is built on

These are the things that break quietly if you do not know them.

- **State is a flat, site-wide structure of arrays.** Every cohort of the whole site lives in one
  contiguous set of 1-D arrays, with patch membership as a CSR map. The dominant daily kernels are
  a single unit-stride sweep.
- **One centralized lockstep reorder.** When you add a per-cohort field, add it to its section of
  `cohort_block` and of each lockstep routine (allocate, free, grow, move, reorder, copy) in
  `state/meds_cohort_state_vars.f90` — the one file that lists every cohort field, every list in
  the same section order. A field that persists across a restart is defined, written and read in
  the same section of the restart file (`io/init/meds_init_restart.f90`). This is the structural
  fix for ED2's "forgot to reallocate one array" class of bug. Patch fields work the same way in `state/meds_patch_state_vars.f90`, whose
  `reorder_patch_block` is what `sort_patches` and `patch_compact` call.
- **Persistent identity.** Every cohort and patch is stamped with a monotonic `global_id` at
  creation and carries it, in lockstep, through every sort, fusion and compaction. Ids are never
  reused; fusion keeps the survivor's. That is what lets a reader track one cohort across output
  records.
- **Tendencies arrive as plain data.** The driver computes a per-cohort tendency bundle; the engine
  applies it. The engine never computes a rate, which is why it is callable from Python with no
  class seam.
- **Conservation is asserted, not assumed.** Every fuse and split conserves aboveground biomass and
  plant number or stops the run. Every fast step closes its store's budget; the slow step closes a
  site-wide ledger. **A closed budget proves bookkeeping, not plausibility** — several real leaks
  hid behind ledgers that balanced, because the ledger declared the same quantity it consumed.
- **Declare a field's fusion kind once.** In the cohort diagnostic block, intensive fields blend
  leaf-area-weighted, extensive ones by plant number, ground-referenced ones sum; in the prognostic
  state every fast field is a per-plant amount (the tissue energies among them) and blends by plant
  number, and an intensive trait blends on the weight that conserves its extensive product (SLA on
  leaf carbon). Getting this wrong at a call site is invisible: the biomass assertion passes and both
  ledgers close.
- **Do not read `site%deriv` from the output layer.** The tendency bundle is transient and
  deliberately not lockstep-reordered, which is correct for its own consumer but wrong for the
  output tick: the step re-sorts the cohorts after filling it. Use the diagnostic block.

## Building one piece

Each kernel library compiles standalone, which is the fastest way to check that rule 1 still holds:

```bash
cmake --build build-ifx --target meds_fast_kernels
cmake --build build-ifx --target meds_slow_kernels
cmake --build build-ifx --target meds_demography
```

Full build, test and run instructions are in [`docs/building.md`](../docs/building.md).

## `test/` is flat, deliberately

`src/` is two levels deep; `test/` is one flat directory. That asymmetry is a **decision, not an
oversight**.

Flat wins because the suite's files have globally unique names, so a mirror buys no lookup that
`test_<subject>.f90` does not already give — while costing a test rename on every source move. The
thing that actually needs to stay honest is the **link line**, not the directory: many tests
deliberately link one narrow library each (`meds_share`, `meds_config`,
`meds_fast_kernels`, `meds_met`) so those layers stay standalone-buildable, and that discipline
is visible in `CMakeLists.txt` whatever directory the file sits in. It is also why the assertion
helpers live in `meds_test_assert`, which depends on nothing but a kind.

Revisit if the suite outgrows a single `ls`, or if two tests ever want the same basename.

## Per-folder notes

Five folders carry their own README with subsystem detail:
[`fast_dynamics/`](fast_dynamics/README.md), [`fast_dynamics/plant/`](fast_dynamics/plant/README.md),
[`io/met/`](io/met/README.md), [`io/maps/`](io/maps/README.md),
[`slow_dynamics/soil/`](slow_dynamics/soil/README.md).
