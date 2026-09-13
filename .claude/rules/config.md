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

A few late-added blocks read each key with a fallback to its in-type default so the feature needs no
TOML edits. **When you do that, the fallback you pass to the reader must equal the in-type default.**
They disagreed once for the soil-carbon switch — the type said on, the loader passed off — and since
the reader returns the supplied default whenever the key is absent, the loader won for every config
that omitted the key, including the shipped one. An entire subsystem ran off for weeks, and it hid
two real defects from every code path anyone ran.

If a value is genuinely required, make it required and let the presence map catch it.

## `[soil_column]` is the ground; `[soil]` is the solver over it

The physical column — layer count, depth, grid growth, hydraulic texture, retention family, root
profile, thermal properties — is validated at load, because a layer count over the compile-time
ceiling or a saturated water content below the residual produces a silently wrong column rather than
a crash.

## Do not add a flag whose off path is known-wrong physics

Several have been deleted for exactly this: the snow switch, two water-balance switches, the soil
advection switch, the multi-layer root switch, the phase-change switch. A flag whose off path is the
cruder physics gets plumbed through every call site forever and eventually reads as "this cannot
happen". Delete the branch instead.

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
