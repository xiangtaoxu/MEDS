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

netCDF is a hard dependency; every configure needs `-DCMAKE_PREFIX_PATH`. Toolchain activation for
this machine is in `CLAUDE.local.md`.

```bash
cmake -S . -B build-ifx -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
cmake --build build-ifx -j
ctest --test-dir build-ifx --output-on-failure          # 45 tests

# Debug (-stand f18 -check all -fpe0) for engine work:
cmake -S . -B build-debug -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
ctest --test-dir build-debug -R fusion_cohort --output-on-failure   # one test by regex
```

CMake auto-resolves Fortran module dependencies. This is the deliberate fix for ED2's "run `make`
six times" hack — **never reintroduce manual object lists or repeated builds.**

**A green ifx run is not sufficient.** Build the nvfortran multicore back end on new modules too.
Three traps, each invisible to ifx and each of which has bitten once:

- **Never pass an array-valued function result straight into a call.** nvfortran's optimizer
  miscompiles the temporary descriptor — silently wrong at `-O2`, segfault at `-O0`. Bind to a named
  array first. (Issue #7.)
- **nvfortran rejects `BLOCK` inside a parallel region.**
- **ifx builds `private` copies of a derived type through a static mold** every thread writes, which
  is why per-thread scratch is an explicit pool, not a data-sharing clause.

## Where a new file goes

Four rules, applied in order. The long form, with the library graph and a lookup table, is in
[`src/README.md`](src/README.md).

1. **Kernels never see `site_t`.** The kernel folders link `state_column` and `config` only, which
   is what keeps them device-eligible and standalone-buildable. A routine that needs `site_t` is
   driver code.
2. **A kernel goes where its caller's timescale is**, then in its domain folder. A kernel called
   from both tiers is a documented seam, never a folder — there is exactly one.
3. **A derived type lives with whoever mutates it.** Two mutators means boundary state, which is
   `state/column`. Parameters are not state.
4. **Laws and operators do not touch.** The appliers take rate arrays and never import the rate
   laws. The slow driver is the one place they meet.

## Invariants that break quietly

- **State is a flat site-wide structure of arrays**; adding a per-cohort field means updating the
  **one centralized lockstep reorder**, every creation site, and the fusion policy.
- **Persistent ids** are stamped at creation and carried through every sort, fusion and compaction.
  Never reused.
- **Tendencies arrive as plain data.** The engine applies; the driver computes. The bundle is
  transient and deliberately not lockstep-reordered — **never read it from the output layer.**
- **A field's fusion kind is declared once**, not at the call site. Getting it wrong is invisible.
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

**Style.** `implicit none` everywhere; explicit interfaces via modules; parameterized kinds (never
bare `real*8` or `1.0d0`); `allocatable` over `pointer`; `pure`/`elemental` for leaf math;
`error stop` rather than a halt routine, so failures are catchable in tests; one responsibility per
module, one module per file, at most 132 columns.

**No hidden global mutable state.** Config and state are derived types passed as arguments. This is
what makes routines unit-testable, thread-safe and reentrant.

**Test as you port.** Every ported kernel gets a CTest target. A port is not done until it
reproduces its reference within tolerance.

## Documentation rules

These are the rules this repository most often broke, so they are stated rather than implied.

1. **A source comment states present-tense rationale.** It may cite an issue number or a science
   section. It does not say what the code used to do.
2. **What changed and when goes in [`CHANGELOG.md`](CHANGELOG.md)**, in the same pull request.
3. **What is deferred goes in [`docs/ROADMAP.md`](docs/ROADMAP.md)** with an issue number — not in a
   `deferred` / `MVP` / `placeholder` comment, which is findable only by grep and goes stale
   silently.
4. **When a design plan's last open item ships**, give it a tombstone and move it to
   `docs/dev_plans/archive/` in the same pull request.

Design-plan **section numbers are never renumbered**: roughly 170 source comments cite them by
number. Most cite by bare filename, so moving a file is safe but renumbering it is not.

## Git

Track Fortran source, CMake files, configs and docs; build artifacts and generated output are
ignored. Line endings are LF via `.gitattributes`. Keep generated netCDF and large input datasets
out of the repository. Commit or push only when asked; if on the default branch, branch first.
