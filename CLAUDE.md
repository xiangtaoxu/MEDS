# CLAUDE.md

Guidance for Claude Code working in this repository.

## What MEDS is

A ground-up reimplementation in **Fortran 2018** of the Ecosystem Demography model **ED2**
(https://github.com/EDmodel/ED2): a size- and age-structured terrestrial biosphere model, advanced
on two coupled timescales. The goal is *not* a line-by-line translation — it is to re-express ED2's
proven algorithms with modern language features so the model is easier to extend, test and reason
about. **When in doubt about *what* a process should compute, ED2 is the reference; when deciding
*how* to structure it, follow the rules here and do not copy ED2's structure.**

ED2 is checked out as a sibling directory, `../ED2`, read-only. Its core is `../ED2/ED/src`, with
the state hierarchy in `memory/ed_state_vars.F90` and the time loop in `driver/ed_model.F90`. The
process-by-process mapping, including which ED2 options MEDS keeps and which it deletes, is in
[`docs/ed2_comparison.md`](docs/ed2_comparison.md).

Background: Moorcroft et al. 2001 (*Ecol. Monogr.*); Medvigy et al. 2009 (*JGR*); **Longo et al.
2019 (*GMD* 12:4309)**, the definitive ED-2.2 technical description.

## Orientation

| To find | Read |
|---|---|
| the source layout, placement rules, library graph | [`src/README.md`](src/README.md) |
| the equations | [`docs/science/`](docs/science/) |
| building, testing, the compiler traps | [`docs/building.md`](docs/building.md) |
| the config surface | [`docs/configuration.md`](docs/configuration.md) |
| what changed and when | [`CHANGELOG.md`](CHANGELOG.md) |
| what is deferred | [`docs/ROADMAP.md`](docs/ROADMAP.md) |
| design records | [`docs/dev_plans/`](docs/dev_plans/) — its README says which are live |

Subsystem detail loads automatically when you touch the relevant files: see `.claude/rules/` for the
fast loop, state and demography, output, config, and the docs math rules.

## Build and test

netCDF is a hard dependency; every configure needs `-DCMAKE_PREFIX_PATH`, which the presets take
from `$CONDA_PREFIX`. Toolchain activation for this machine is in `CLAUDE.local.md`.

```bash
FC=ifx cmake --preset dev               # the dev kind (Release, -O2) in build-dev
cmake --build build-dev -j
ctest --test-dir build-dev --output-on-failure          # 86 tests, ~90 s; -j 5: 33 s
ctest --test-dir build-dev -L unit --output-on-failure  # the 67 unit tests, ~15 s

# debug (-stand f18 -check all -fpe0) for engine work; the whole suite passes here too:
FC=ifx cmake --preset debug
ctest --test-dir build-debug -R slow_update --output-on-failure   # one test by regex
```

The tests mirror `src/`, one per module (ctest `slow_update` is `test/slow/test_slow_update.f90`),
and each has one label: `unit` (67), `integration` (5), `lint` (10) or `python` (4).
`-DMEDS_TEST_LONG=ON` adds the two `long` Palo Verde runs, to make before a physics change merges
(`docs/building.md`, "Running the tests").

Four build kinds, one preset each (`docs/building.md`, "Build kinds"): `debug`, `dev`, `profile`
(dev with `-g`) and `production` (IPO and `-xCORE-AVX512`, 26–36 % faster, for long runs; it builds
`meds_main` and `libmeds.so`, not the tests, whose bit-exact checks it breaks at round-off; the
failures and production's output against dev's are in `docs/building.md`, "Running the tests"). A
configure that names no build type builds dev. dev is the reference: compare outputs between commits
on dev builds, never on production ones. dev builds `meds_main`, `libmeds.so` and the tests, and
`build-dev` is the default: the example scripts run `build-dev/meds_main` and the Python package loads
`build-dev/libmeds.so`, unless `--meds-main` (or `MEDS_MAIN`) or `MEDS_LIB` names another build.

CMake auto-resolves Fortran module dependencies. This is the deliberate fix for ED2's "run `make`
six times" hack — **never reintroduce manual object lists or repeated builds.**

**A green ifx run is not sufficient.** Build the nvfortran multicore back end on new modules too,
and gfortran. Seven traps, each invisible to ifx and each of which has bitten at least once:

- **Never pass an array-valued function result straight into a call.** nvfortran's optimizer
  miscompiles the temporary descriptor — silently wrong at `-O2`, segfault at `-O0`. Bind to a named
  array first.
- **nvfortran rejects `BLOCK` inside a parallel region.**
- **ifx builds `private` copies of a derived type through a static mold** every thread writes, which
  is why per-thread scratch is an explicit pool, not a data-sharing clause.
- **Never pass a component section `a(:)%c` as an actual argument when `c`'s type has allocatable
  components.** gfortran copies it through a temporary whose copy-out dangles them. Keep such
  objects in a contiguous array of their own.
- **Never nest an implied-do inside a `pack` inside an array constructor.** nvfortran returns wrong
  elements at `-O2` and segfaults at `-O0`. Pack a named array.
- **Keep `findloc` off LOGICAL arrays and off shorter character values.** nvfortran aborts on the
  first ("unimplemented") and returns 0 for `findloc(names, trim(s))`. Search `merge(1, 0, mask)`;
  pass the untrimmed string.
- **Never reset with `x = t()` when `t` has a fixed-size array component of a type with
  allocatable components.** nvfortran allocates a garbage size. Assign a default-initialised local.

## Where a new file goes

Four rules, applied in order. The long form, with the library graph and a lookup table, is in
[`src/README.md`](src/README.md).

1. **Kernels never see `polygon_t`.** The kernel folders take only the per-patch store types from the
   state layer (`meds_patch_state_vars`, through an only-list), which is what keeps them
   device-eligible and standalone-buildable; ctest `kernels_patch_only` checks it. A routine that
   needs `polygon_t` is driver code. The integrator also takes the patch's site (`site_surface_t`) and
   nothing past it; ctest `integrator_boundary` checks it.
2. **A kernel goes where its caller's timescale is**, then in its process module (fast:
   `fast_dynamics/processes/`; slow: `slow_dynamics/processes/`). A kernel called
   from both tiers is a documented seam, never a folder — there is exactly one.
3. **A derived type lives with whoever mutates it.** Two mutators means state, in the `state/` file
   of the level that holds it (`meds_cohort_state_vars`, `meds_patch_state_vars`,
   `meds_site_state_vars`, `meds_polygon_state_vars`). Parameters are not state.
4. **Laws and operators do not touch.** The appliers (`meds_slow_update`) take rate and decision
   arrays and never import the laws (`meds_demography`). `meds_slow_dynamics` is the one place they meet.

## Invariants that break quietly

- **State is a flat site-wide structure of arrays**; adding a per-cohort field means adding it to its
  section of `cohort_block` and of **each lockstep list** in `meds_cohort_state_vars` (one section
  order throughout), every creation site, the fusion policy, and, if it persists, the same section
  of the restart file (`meds_init_restart`: define, write, read). Patch fields: the same, in
  `meds_patch_state_vars`.
- **Persistent ids** are stamped at creation and carried through every sort, fusion and compaction.
  Never reused.
- **Tendencies arrive as plain data.** The engine applies; the driver computes. The bundle is
  transient and deliberately not lockstep-reordered — **never read it from the output layer.**
- **A field's fusion kind is declared once**, not at the call site. Getting it wrong is invisible.
- **Only prognostic state has tendencies.** A cohort's state is `nplant` and its four carbon pools.
  Its geometry (dbh, height, basal area, AGB, leaf area, wood caches) has one writer,
  `set_cohort_size_from_carbon`, called after every operator that moves carbon or SLA; it is never
  integrated, blended or read back from a file on its own. Integrating leaf area beside leaf carbon
  drifted it to -1e-20 on a leafless cohort and crashed a run 50 years later.
- **Bounds hold by construction and are checked every step.** An operator that moves an amount
  limits the outflow to what is there; a clamp's bound comes from a quantity that cannot go
  negative. `cohort_bad_field` and `soil_carbon_bad_pool` stop the run on the step a bound breaks,
  naming the field: a new bounded field goes there.
- **An intensive quantity blends on the weight that conserves its extensive product**: SLA on leaf
  carbon, vcmax25 on leaf area, tissue temperature on heat capacity. Blend with `pair_mean`, which
  has no size threshold: a leafless cohort's leaf carbon of 1e-193 is still a weight, and a 1e-30
  cutoff there broke the leaf area fusion conserves.
- **A conservation tolerance scales with the terms, not the net.** `conserved(before, after,
  magnitude, rtol)` with magnitude = sum of |terms|: a net near zero turned a 1e-12 relative
  tolerance into 1e-42 absolute. Its only absolute floor is the format's, tiny/epsilon
  (~1e-292), below which an optimized build flushes the arithmetic to zero.
- **Conservation is asserted, not assumed** — but **a closed budget proves bookkeeping, not
  plausibility.** Several real leaks sat behind ledgers that balanced, because the ledger declared
  the same quantity it consumed. Watch consecutive-step traces.
- **Measure a filed issue's premise before building its fix.** Two items have dissolved or inverted
  under measurement.
- **A green suite proves nothing** if the code is outside every build path, or the assertions are
  regime-blind, or a permissive mock is standing in for the real library.

## Conventions

**Names.** Spell out ecologically or physically meaningful words — `site`, `cohort`, `patch`,
`dbh_critical`, `hgt_max`. Keep terse loop indices and the established domain tokens `dbh`,
`nplant`, `pft`. Inside `associate`, alias to the full word. Readability beats succinctness when
performance is equal. Name a helper by its operation, not by its caller's domain.

**Constants.** A named constant is spelled by its kind, not by who uses it. The physical,
mathematical and unit-conversion constants and the kinds, in `meds_constants` and `meds_kinds`, are lower
case: they read as symbols in a formula (`grav`, `rho_h2o`, `day_sec`, `wp`). Every other named constant
is UPPER case, public or private, module-level or local: a selector code, an index or a flag, with one
prefix per family (`SOIL_RETENTION_VG`, `INTEG_ESDIRK2`), and a number chosen in the code, a tolerance, a
floor, a table size or a fitted coefficient (`RWC_FLOOR`, `HYDRO_TABLE_NTAB`). A constant keeps its name
when another module starts to use it. ctest `constants_named` (`scripts/lint/check_constant_case.py`)
enforces it.

**Style.** `implicit none` everywhere; explicit interfaces via modules; parameterized kinds (never
bare `real*8` or `1.0d0`); `allocatable` over `pointer`; `pure`/`elemental` for leaf math;
`error stop` rather than a halt routine, so failures are catchable in tests; one responsibility per
module, one module per file, at most 132 columns.

**License header.** Every source file starts with an `SPDX-License-Identifier: Apache-2.0` comment,
on the line after the `#!` line if there is one. A new contributor adds themselves to `AUTHORS`; see
`CONTRIBUTING.md`.

**No hidden global mutable state.** Config and state are derived types passed as arguments. This is
what makes routines unit-testable, thread-safe and reentrant.

**No bare thresholds.** A guard against dividing by zero or leaving a function's domain uses
`tiny_num` (or `safe_exp`). A physical threshold is a named setting or parameter, with its source.
A literal clamp that is right where it is (a solver's step size, a function's domain) says why on the
line: `! clamp-ok: <reason>`. ctest `clamps_named` (`scripts/lint/check_clamps.py`) enforces it:
`max(lai, 0.1)` once spread a thin cohort's light over leaf area it did not have.

**Test as you port.** Every ported kernel gets a CTest target. A port is not done until it
reproduces its reference within tolerance.

## Documentation rules

These are the rules this repository most often broke, so they are stated rather than implied.

1. **A source comment states present-tense rationale.** It says what the code does and why, and
   cites the science page it implements (`docs/science/<page>.md §4.1`, `eq. (7)`). It never cites
   an issue or pull-request number and does not say what the code used to do: that is history.
   ctest `no_issue_refs` (`scripts/lint/check_issue_refs.py`) keeps issue and pull-request numbers
   out of every file except `CHANGELOG.md`, `docs/ROADMAP.md` and `docs/dev_plans/`.
2. **What changed and when goes in [`CHANGELOG.md`](CHANGELOG.md)**, in the same pull request.
3. **What is deferred goes in [`docs/ROADMAP.md`](docs/ROADMAP.md)** with an issue number — not in a
   `deferred` / `MVP` / `placeholder` comment, which is findable only by grep and goes stale
   silently.
4. **When a design plan's last open item ships**, give it a tombstone and move it to
   `docs/dev_plans/archive/` in the same pull request.

Design-plan **section numbers are never renumbered**: about 50 source comments still cite a plan by
section, where the plan is the only written record of a rationale; everything a science page covers is
cited there instead. Most cite by bare filename, so moving a file is safe but renumbering it is not.
**Science-page equation numbers are cited from the code too**: inserting or removing an equation
renumbers the rest of its page, so update every citation in the same change
(`.claude/rules/docs-math.md`).

## Git

Track Fortran source, CMake files, configs and docs; build artifacts and generated output are
ignored. Line endings are LF via `.gitattributes`. Keep generated netCDF and large input datasets
out of the repository. Commit or push only when asked; if on the default branch, branch first.

**`beta` is the integration branch.** Branch from `beta` and open pull requests against it, not
`main`: `beta` collects merged work that is not yet in a release, and merges into `main` when a
release is cut. A closing keyword (`Fixes #N`) in a pull request into `beta` does not close the
issue, because GitHub acts on those only for the default branch; close it when the release lands.
