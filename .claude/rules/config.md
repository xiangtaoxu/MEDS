---
paths:
  - "src/config/**"
  - "**/*.toml"
---

# Configuration

User-facing description: [`docs/configuration.md`](../../docs/configuration.md).

## No hard-coded model parameters

The source defines only true **constants** — numerical, geometric, calendar. Every model parameter
is required from TOML across two files: a main file with all non-PFT settings, which names a PFT
file carrying the trait table, the mortality-hazard coefficients and the allometry coefficients.

`load_meds_config` builds a **presence map** while reading and stops listing every missing key. A
missing file is also a hard error. There is no `build_config` and no defaults path: derived
quantities come from the derive routines, overridable via an options flag plus a derived block.

**Tests get their config from `build_test_config()` in `test/meds_test_support.f90`** — the only
place "default" values live in code.

## The absent-key trap

A few blocks read each key with a fallback to its in-type default so the feature needs no TOML
edits. **When you do that, the fallback you pass to the reader must equal the in-type default.** The
reader returns the supplied fallback whenever the key is absent, so if the two disagree the loader
wins for every config that omits the key, including the shipped one. A subsystem switch the type
declares on, read with an off fallback, turns the whole subsystem off in silence, and hides every
defect in it from every code path anyone runs.

If a value is genuinely required, make it required and let the presence map catch it.

## `[soil]` is the ground; `[fast]` holds the solvers' settings

Every property of the ground is one block, whichever solver reads it: the column, its surface and
both bottom boundaries, so one Fortran type does not map to one block (`soil_column_config_t`,
`soil_opts_t` and `energy_opts_t` are each filled from `[soil]`, the last two also from `[fast]`).
The physical column — depth and layer bottoms, hydraulic texture, retention family, thermal
properties — is validated at load, because a layer count over the compile-time ceiling, layer
bottoms that do not increase or a saturated water content below the residual produce a silently
wrong column rather than a crash. The root profile is a plant trait: it lives in `[hydraulics]`,
with optional per-PFT values in the PFT file.

## Do not add a flag whose off path is known-wrong physics

That is why there is no switch for snow, the water balance, soil advection, multi-layer roots or
phase change. A flag whose off path is the cruder physics gets plumbed through every call site
forever and eventually reads as "this cannot happen". Delete the branch instead.

The corollary: a config value that parses but does nothing is worse than one that is absent. If a
selector is unimplemented, reject it in `validate_config`.

## Where a config type lives

A per-domain `*_opts` leaf in `src/config/`, so `meds_config` — the root of the dependency graph —
can carry it with no back-edge into a kernel library, and so the sealed kernels stay
device-eligible. Never put a config type in a kernel module.

## Renaming

netCDF registry strings and TOML keys are **not** renamed with Fortran identifiers. A rename pass
that touches a string literal changes output files and breaks user configs; one that touches only
identifiers cannot.
