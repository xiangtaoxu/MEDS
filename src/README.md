# `src/` — the MEDS source tree

**30 k lines, 86 modules, 19 CMake libraries.** This page is the map: the layout, the four rules
that decide where a new file goes, and the library graph the build enforces. It is the page to
read before adding code, and the one to hand someone who asks how MEDS is organized.

The science is in [`docs/science/`](../docs/science/); the design record is in
[`docs/dev_plans/`](../docs/dev_plans/).

## The shape, in one sentence

MEDS is an **operator-split fast/slow model**, so the tree splits **by timescale** first and by
**domain** second, over a **state layer** that both halves share.

```
src/                                    30.3 k lines · 86 modules · 19 CMake libraries
│
│  ── FOUNDATION ─────────────────────── no model state, no process: libmeds_shared
├── shared/        uses nothing outside itself; every layer may use it
│   ├── base/      kinds, physical and calendar constants                      109
│   ├── functions/ stateless constitutive laws: allometry, thermodynamics,
│   │              retention and pressure-volume curves, canopy optics,
│   │              temperature response                                      1 085
│   └── util/      calendar time, numerics (matrix exponential, Thomas sweep,
│                  root finders), budget checks                                882
│
│  ── CONFIGURATION ──────────────────── the run's inputs: libmeds_config
├── config/        TOML reader and loader, PFT trait table, meds_config_t,
│                  and the per-domain *_opts leaves                          2 920
│
│  ── STATE, a layer in two halves ───────────────────────────────────────────
├── state/
│   ├── column/    ONE patch, seen vertically: the soil-water, soil-energy,
│   │              snow, canopy-air and soil-carbon reservoirs, their fusion
│   │              blends, and the parameter bundles that describe them         558
│   └── site/      ALL the patches: the flat cohort structure-of-arrays, the
│                  patch CSR map, the lockstep reorder machinery, site_t,
│                  and the diagnostic accumulators                           1 674
│
│  ── PROCESSES, timescale first ──────────────────────────────────────────────
├── fast_dynamics/       sub-daily, one dt_fast per step                     10 977
│   ├── canopy/    the medium: two-stream radiation, aerodynamics, the
│   │              canopy-air box                                              983
│   ├── plant/     the organisms: leaf gas exchange, hydraulics, maintenance
│   │              respiration, tissue energy                                1 606
│   ├── soil/      the ground: soil water, soil energy, the ground skin,
│   │              the snow store                                            1 472
│   ├── numerics/  the integrator machinery: the state vector, the frozen
│   │              work record, ARK and RK45, error control, the pre-pass,
│   │              and a test-only RK4 oracle                                5 378
│   └── driver/    walks one slow step in dt_fast sub-steps over the patches 1 538
│
├── slow_dynamics/       daily to annual                                      4 623
│   ├── plant/     phenology, carbon allocation, trait plasticity               587
│   ├── soil/      the CENTURY soil-carbon matrix, the litter partition         912
│   ├── demography/ the vital-rate LAWS, and the operators that apply them:
│   │              state update, cohort and patch fuse-fission, recruitment,
│   │              treefall disturbance                                      1 308
│   └── driver/    the slow coordinator, the vegetation and biogeochemistry
│                  drivers, and the slow conservation ledger                 1 816
│
│  ── EDGES ────────────────────────────────────────────────────────────────────
├── forcing/       prescribed drivers: the met reader and its disaggregation
│                  kernels                                                   1 067
├── io/            netCDF C bindings, the restart stream, and the diagnostic
│                  wall: derive → capture → reduce → integrate → serialize   4 563
├── init/          the initial community: bare ground, or a cohort census       181
├── c_api/         bind(c) shims → one libmeds.so: leaf, phenology,
│                  demography, and the full coupled run                        887
└── main/          meds_stepper (the cadence owner), meds_driver
                   (open / step / finalize), meds_main (the PROGRAM)           740
```

**`shared/` is one library, and it admits only what depends on nothing outside it.** Every module
there `use`s other `shared/` modules and compiler intrinsics and nothing else: no model state, no
configuration, no external library. That is what lets every other layer, and the tests that link
one narrow library, use it freely. A file that needs `meds_config`, a state type or netCDF is not
foundation however small or widely used it is, and goes to its caller's layer instead; `config/`
itself is a separate library for exactly that reason, since it links `state_column`. The test is
what keeps `shared/` from becoming the place a file goes when no other place fits.

## Four rules, in the order you apply them

**1. Kernel folders never see `site_t`.** `fast_dynamics/{canopy,plant,soil}`,
`slow_dynamics/{plant,soil}` and `state/column` link only `state_column` and `config`. That is what
keeps each kernel library building standalone, keeps the kernels eligible for OpenMP `target`
offload (a portability property worth holding, though offload is not currently a speedup — see
[`dev_plans/MEDS_GPU_EVALUATION.md`](../docs/dev_plans/MEDS_GPU_EVALUATION.md)), and lets twelve of
the tests link one kernel library alone. It is checked, not assumed:
there is no occurrence of `site_t` in any of those folders. **If a new routine needs `site_t`, it
is driver code**, and it belongs in a `driver/` folder.

**2. A kernel goes where its caller's timescale is**, and within that, in its domain folder. A
kernel called from *both* tiers is a documented seam, never a folder of its own. There is exactly
one: `heterotrophic_respiration_matrix`, so that the sub-daily respiration debits the same CENTURY
matrix the daily step does. Co-locating it with the daily step is what closes the seam gap to
machine precision.

**3. A derived type lives with whoever mutates it.** If two subsystems mutate it, it is boundary
state and belongs in `state/column`. **Parameters are not state**: they are derived once and never
integrated, so `meds_column_params` is a separate module from `meds_column_state_types`.

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
  shared ──→ state_column ──→ config ──→ state_site ──→ demography ──→ io_prep
     │            │                          ↑              │             │
     │            ├──→ fast_kernels ─────────┤              │             │
     │            └──→ slow_kernels ─────────┘              │             │
     │                                                      │             │
  netcdf_c ──→ forcing ───────────────────────────────→  fast  ←──────────┤
                                                          slow  ←─────────┘
                                                          init
                                                            └──→ stepper ──→ model
                                                                             (INTERFACE)
  io_stream ─────────────────────────────────────┐
  model ─────────────────────────────────────────┴──→ driver ──→ meds_main
                                                               └→ meds_py (libmeds.so)
```

Two edges are worth knowing because they are easy to get wrong:

- **`config` links `state_column`**, not the other way round: `meds_config` takes
  `n_soil_layer_max` from `meds_column_params` for the `[soil_column]` bounds check. The edge was
  missing until a Ninja build exposed it — the Makefile generator happened to order the two module
  files correctly, so a parallel build with a different scheduler failed where the everyday one
  passed.
- **`io_prep` carries no netCDF.** The netCDF-free half of the output subsystem (the reductions,
  the registry, the temporal integrators) is a separate target from the serializer, so the stepper
  can tick the integrators without pulling a C dependency. It is an explicit file list, not a glob:
  adding a diagnostic module there is a deliberate CMake edit.

## Where does a new file go?

| You are adding | It goes in |
|---|---|
| a soil, canopy or plant process | `fast_dynamics/<domain>/` or `slow_dynamics/<domain>/`, by the **cadence of its caller** |
| a reservoir two subsystems mutate | `state/column/` |
| a per-cohort field | `state/site/` — and the lockstep reorder, every creation site, and the fusion policy |
| a stateless constitutive curve | `shared/functions/` |
| a numerical, calendar or checking helper every layer needs | `shared/util/`, if it uses nothing outside `shared/` |
| a TOML block | `config/`, as an `*_opts` leaf so `meds_config` can carry it with no back-edge |
| an output variable | one `add_variable` line in `io/meds_output_registry.f90` |
| a Python entry point | `c_api/`, plus its mandatory ctest target |
| a top-level orchestration step | `main/meds_stepper.f90` or a `driver/` folder |

## Invariants the tree is built on

These are the things that break quietly if you do not know them.

- **State is a flat, site-wide structure of arrays.** Every cohort of the whole site lives in one
  contiguous set of 1-D arrays, with patch membership as a CSR map. The dominant daily kernels are
  a single unit-stride sweep.
- **One centralized lockstep reorder.** When you add a per-cohort field, update the reorder,
  compaction, slot-copy and capacity routines in `state/site/meds_site_state_types.f90` — the
  single place that touches every array. This is the structural fix for ED2's "forgot to reallocate
  one array" class of bug. Patch arrays have no single reorder routine; their permute and pack
  sites are `sort_patches` and `patch_compact`.
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
- **Declare a field's fusion kind once.** Intensive fields blend leaf-area-weighted, extensive ones
  by plant number, ground-referenced ones sum. Getting this wrong at a call site is invisible: the
  biomass assertion passes and both ledgers close.
- **Do not read `site%deriv` from the output layer.** The tendency bundle is transient and
  deliberately not lockstep-reordered, which is correct for its own consumer but wrong for the
  output tick, which runs after the monthly fuse-fission. Use the diagnostic block.

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
oversight** (#193, structure-plan decision #13, taken 2026-09-13).

Flat wins because the suite is ~47 files with globally unique names, so a mirror buys no lookup that
`test_<subject>.f90` does not already give — while costing a test rename on every source move, in a
tree that has moved twice. The thing that actually needs to stay honest is the **link line**, not the
directory: seventeen tests deliberately link one narrow library each (`meds_shared`, `meds_config`,
`meds_fast_kernels`, `meds_forcing`) so those layers stay standalone-buildable, and that discipline
is visible in `CMakeLists.txt` whatever directory the file sits in. It is also why the assertion
helpers live in `meds_test_assert`, which depends on nothing but a kind.

Revisit if the suite outgrows a single `ls`, or if two tests ever want the same basename.

## Per-folder notes

Four folders carry their own README with subsystem detail:
[`fast_dynamics/`](fast_dynamics/README.md), [`fast_dynamics/plant/`](fast_dynamics/plant/README.md),
[`forcing/`](forcing/README.md), [`slow_dynamics/soil/`](slow_dynamics/soil/README.md).
