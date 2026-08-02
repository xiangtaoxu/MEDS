# MEDS v0.1.0, for people who use ED2

MEDS is a ground-up reimplementation of [ED2](https://github.com/EDmodel/ED2) in Fortran 2018. If you
already run ED2, this page tells you what carried over, what changed, what is missing, and what a MEDS
run costs you to set up.

**Reference points.** MEDS at tag `v0.1.0`; ED2 mainline at commit `125f814d` (ED-2.2). ED2
options are named by their `ED2IN` namelist key, and "ED2 default" means the value in the shipped
`ED/run/ED2IN`.

**One thing to establish up front: MEDS has not been benchmarked against ED2.** No EDTS-equivalent
regression suite has been run, no site has been compared flux-for-flux, and no output of any kind has
been scored against observations. What has been verified is internal — 38 unit tests on two compilers,
per-step conservation ledgers, and thread-invariant output. Everything below describes *what the code
does*, not *how well it does it*. Treat the numbers a MEDS run produces as a working model's numbers.

---

## 1. Summary

**What is the same.** MEDS keeps ED2's central idea intact: an ecosystem is a size- and
age-structured demographic hierarchy, advanced on two coupled timescales. Cohorts are plants of one
PFT in one size class; patches are land sharing a disturbance history; fusion and fission are
first-class operations that keep resolution bounded. Underneath runs a sub-daily biophysics loop —
two-stream canopy radiative transfer, a prognostic canopy air space, Monin–Obukhov aerodynamics, leaf
photosynthesis with stomatal coupling, and soil columns carried as internal energy — feeding a
daily/monthly/annual vegetation-dynamics loop that does growth, mortality, recruitment, and
disturbance. Most individual algorithms are recognisably ED2's: the `icanrad = 2` two-stream solver,
the Nusselt leaf boundary layers, the CENTURY decomposition cascade, the `iallom = 3` pan-tropical
allometry, ED2's negative-*z* soil geometry, Chambers-2004 stem respiration.

**What is different, in one sentence each.**

- **Scope is much narrower.** One site, one soil column, no fire, no land use or harvest, no nitrogen,
  no MPI, no gridded/regional runs. MEDS v0.1 is a site model.
- **Options are decisions, not switches.** ED2 exposes ~40 scheme selectors (`ICANRAD`, `IPHEN_SCHEME`,
  `IALLOM`, `DECOMP_SCHEME`, …), most with legacy and beta branches. MEDS mostly implements one path —
  usually ED2's default or its best-supported alternative — and deletes the rest. Where MEDS keeps a
  choice, it is because both branches are defensible.
- **Some defaults are more aggressive than ED2's.** Plant hydraulics is always on and unconditional
  (ED2 default `PLANT_HYDRO_SCHEME = 0`, no hydraulics). The soil-water solver is implicit Richards
  with van Genuchten retention by default (ED2 default is Campbell/Cosby inside its RK4). Phenology is
  unconditional and rate-based.
- **Some are more conservative.** Soil freeze/thaw is opt-in (`[energy].phase_change`, default off).
  Soil carbon is opt-in. Snow is a single bulk layer against ED2's multi-layer temporary surface water.
- **Configuration and output are wholly rebuilt.** TOML instead of `ED2IN` + XML; netCDF instead of
  HDF5; and roughly 200 output variables individually switchable per timescale rather than a fixed
  schema gated by frequency flags.
- **The code is unit-tested and thread-deterministic.** 38 CTest targets; output byte-identical at any
  OpenMP thread count.

**Who should look at it.** If you need fire, land use, a nitrogen cycle, regional runs, or the ED2
PFT parameterisations validated across biomes, stay on ED2 — MEDS does not have them. If you are doing
single-site process work on plant hydraulics, canopy biophysics, sub-daily carbon–water coupling, or
numerics, MEDS is a substantially smaller and more tractable codebase (~75 source files against ED2's
`ed_state_vars.F90` alone being 34k lines), with the conservation and determinism properties spelled
out.

---

## 2. Processes

Read the notes column — that is where the qualifications live.

### 2.1 Structure and demography

| | MEDS v0.1.0 | ED2 (ED-2.2) | Notes |
|---|---|---|---|
| **State hierarchy** | site → patch → cohort | grid → polygon → site → patch → cohort | MEDS drops the grid/polygon levels and supports one site. State is a flat site-wide Structure-of-Arrays with a CSR patch map, not nested ragged arrays. |
| **PFTs** | Run-time count; every trait supplied from a TOML file. No built-in table. | 17 hard-coded PFTs with defaults in `ed_params.f90`, overridable by XML. | MEDS has **no PFT defaults at all** — a missing key is a hard error. That is a deliberate trade: no hidden parameterisation, but no curated ED2 PFT set to inherit either. The shipped example uses three PFTs contrasted along wood density. |
| **Allometry** | Pan-tropical, ED2 `IALLOM = 3` | `IALLOM` 0–4, default 3 | Same family. MEDS additionally inverts AGB→DBH for fusion. |
| **Cohort fusion / fission** | Keys on height similarity + an LAI cap; conserves per-plant AGB (DBH re-derived, never averaged) and plant number | `fuse_fiss_utils.f90`, `MAXCOHORT` target with tolerance relaxation | Same idea, different similarity metric. MEDS asserts the conservation invariants at 1% and `error stop`s otherwise. |
| **Patch fusion** | ED2-style cumulative-LAI light profile over height layers, same disturbance class only | `IFUSION` 0/1, `MAXPATCH` | MEDS ports the light-profile criterion. |
| **Treefall disturbance** | Annual; a fraction of every patch peels into one age-0 gap; tall die, short understory survives | `TREEFALL_DISTURBANCE_RATE`, plus small-tree background mortality | Directly ported, including the tall/short split. |
| **Fire** | **Absent** | `INCLUDE_FIRE` 0–3, `SM_FIRE` | |
| **Land use / harvest / plantation** | **Absent** | `IANTH_DISTURB`, the `SL_*` / `CL_*` logging block, agriculture and pasture stock PFTs | |
| **Cohort/patch identity** | Persistent `global_id`, never reused, survives every sort/fuse/compaction | Not tracked | Lets an external reader follow one cohort through the output until it fuses away. |

### 2.2 Plant physiology

| | MEDS v0.1.0 | ED2 (ED-2.2) | Notes |
|---|---|---|---|
| **Photosynthesis** | FvCB (C3) with Rubisco / RuBP-regeneration / TPU limitation; Collatz (C4). Arrhenius or peaked temperature response | `IPHYSIOL` 0–3; default 2 (Collatz Q10 with Moorcroft high/low-T corrections). Options 1 and 3 add Jmax and TPU | MEDS is closest to `IPHYSIOL = 3` in structure (explicit Jmax + TPU) but exposes the temperature response as a per-trait choice rather than a bundled scheme number. |
| **Stomatal conductance** | Leuning 1995, Medlyn 2011 (USO), Katul 2010 optimality | `ISTOMATA_SCHEME` 0 (Leuning) or 1 (Katul) | Medlyn USO is the addition. |
| **A–Ci coupling** | One nonlinear equation in Ci, solved by bracketing + bisection; uniform across all three gs models and both pathways | Analytical solution per case in `farq_leuning.f90`; separate solver in `farq_katul.f90` | MEDS trades ED2's analytic speed for one code path that every model shares. |
| **Water stress on carbon** | Two limbs. Stomatal limb (Sabot 2022, driven by soil water potential) **on**; non-stomatal capacity limb (ramp on Vcmax/Jmax/TPU in leaf water potential) **off by default** | `H2O_PLANT_LIM` 0–5, default 2 (supply/demand FSW from rooting-zone water) | Different formulation entirely. The capacity limb is off because its two parameters are weakly constrained and it is rarely measured directly, and because leaf water potential is not converged at production `dt_fast` (§3.3) — wiring an unconverged potential into carbon through a poorly constrained ramp is not a trade worth making. The limb's formulation is itself under review. |
| **Plant hydraulics** | **Always on, unconditional.** 2- or 3-node network (leaf / wood / root), nonlinear Bartlett pressure–volume capacitance, Kirchhoff-integrated xylem conductance, per-layer rhizosphere conductances, solved exactly by matrix exponential over the step | `PLANT_HYDRO_SCHEME` 0 (default, off — tissues always saturated), 1 or 2 (Xu 2016 / Christoffersen 2016) | The single biggest default difference. MEDS has no "no hydraulics" path. Hydraulic redistribution is deliberately **not** enabled (per-layer root efflux floored at zero). |
| **Trait plasticity** | Light-driven SLA / Vcmax25 / Rd25 / leaf lifespan, turnover-limited; opt-in, default off | `TRAIT_PLASTICITY_SCHEME` 0–3 and −1/−2, default 0 | MEDS's single scheme corresponds most closely to ED2's option 3 (change constrained by leaf turnover). |
| **Phenology** | One generic rate engine: daily cues → two governor drives → a flush rate and a shed rate, both per-day. Per-PFT cue masks select which cues drive which governor. Unconditional | `IPHEN_SCHEME` −1 … 4 (evergreen / drought-deciduous old and new / prescribed / light / hydraulic) | This is a genuine restructuring, not a port. ED2's scheme numbers become per-PFT cue masks, so cold-deciduous, drought-deciduous and light-driven habits coexist in one run without a global switch. Prescribed phenology from files (`IPHEN_SCHEME = 1`, `PHENPATH`) has **no MEDS equivalent**. |
| **Carbon allocation** | FATES PARTEH-H1: four pools (leaf, fine root, wood, non-structural), allometric targets, daily priority ladder, wood as the residual sink and the prognostic size anchor | `growth_balive` (daily) + `structural_growth` (monthly), `ISTRUCT_GROWTH_SCHEME` 0/1 | Unified into one daily step. Growth respiration is charged on realized growth only. |
| **Maintenance respiration** | Leaf (from the gas-exchange kernel), stem (Chambers 2004, surface-area based), fine root | `GROWTH_RESP_SCHEME`, `STORAGE_RESP_SCHEME`, `ISTEM_RESPIRATION_SCHEME` (default 1 = Chambers) | Same stem formulation as ED2's current default. |
| **Mortality** | One hazard: Camac 2018 additive `gamma + alpha·exp(−beta·growth)`, on a tracked moving-average carbon growth rate, with per-PFT coefficients as power laws in wood density. Plus treefall. | Six additive pathways: ageing, carbon-starvation (`CARBON_MORTALITY_SCHEME`, `IDDMORT_SCHEME`, `CBR_SCHEME`), treefall background, cold/frost, hydraulic failure (`HYDRAULIC_MORTALITY_SCHEME`), and disturbance | **The largest science gap.** MEDS has no frost mortality, no explicit hydraulic-failure mortality, and no separate carbon-starvation term — growth-dependence carries all of it. ED2's Camac option (`CARBON_MORTALITY_SCHEME = 2`) is the closest analogue. |
| **Recruitment** | Reproduction carbon → recruits **within the parent patch**, plus a per-PFT external seed rain; pooled until a minimum size, then spawned | `REPRO_SCHEME` 0–3, default 2 (seeds exchanged among all patches of a polygon) | MEDS does **not** disperse seeds between patches. In a strongly heterogeneous stand this matters. |

### 2.3 Biophysics

| | MEDS v0.1.0 | ED2 (ED-2.2) | Notes |
|---|---|---|---|
| **Canopy radiative transfer** | Two-stream, multi-layer, one solver run per band (VIS, NIR, thermal LW); Beta leaf-angle distribution; absolute W m⁻² throughout | `ICANRAD` 0/1/2, default 2 (Liou two-stream) | Ported from `twostream_rad.f90`. MEDS generalises the band loop and drops ED2's normalise-to-unit-incidence convention. No horizontal shading (`IHRZRAD`), no finite crown radius (`CROWN_MOD`). |
| **Canopy turbulence** | CLM5 Monin–Obukhov surface layer + per-cohort exponential wind extinction + ED2 Nusselt leaf/wood boundary layers + CLM ground conductance | `ICANTURB` 0–4, default 2 (Massman 1997); `ISFCLYRM` 1–4, default 3 (Beljaars–Holtslag) | MEDS's combination is closest to ED2's `ICANTURB = 4` / `ISFCLYRM = 4` (the CLM-based branches), not to the ED2 defaults. |
| **Canopy air space** | Prognostic temperature, humidity **and CO₂** — three twins in one box model | Prognostic in the RK4 state, same three | Equivalent. MEDS's CAS depth is per-patch state derived from the tallest cohort plus freeboard. |
| **Leaf / wood energy** | Separate leaf and wood temperatures, relaxed **exactly** over the step (the tissue ODE is linear under the frozen coefficients, so a closed-form exponential with an endpoint weight and a step-average weight) | `IBRANCH_THERMO` 0/1/2, default 1 (branch and leaf as one pool in the biophysics) | MEDS is closest to `IBRANCH_THERMO = 2`. Caveat carried in the code: the wood area index and sapwood placeholders make the modelled wood time constant several times too small. |
| **Soil water** | Implicit backward-Euler Thomas solve of Richards, Celia modified-Picard or frozen-coefficient linearisation, upstream-weighted K, adaptive step doubling. van Genuchten (default) or Campbell retention | Integrated inside the RK4 patch state; Campbell(-Mualem) with Cosby PTF by default (`SOIL_HYDRO_SCHEME` 0), Tomasella-Hodnett or van Genuchten as beta options 1/2 | Different numerical treatment (see §3) and a different default retention curve. |
| **Soil hydraulic parameters** | Given directly per column (saturated water content, Ksat, curve parameters); uniform texture broadcast over layers; no pedotransfer function, no texture classes | 17 texture classes, PTFs on sand/silt/clay (+ SOC, pH, CEC, bulk density for scheme 2), soil-texture and soil-depth map databases | **A real setup difference.** ED2 will build a soil column from a texture map; MEDS wants numbers. |
| **Soil bottom BC** | free drain / bedrock / aquifer, where the aquifer is a head-driven two-way boundary against a saturated zone at the column base | `ISOILBC` 0–3 (bedrock / free drainage / lateral drainage with `SLDRAIN` / aquifer) | No lateral-drainage-by-slope option in MEDS. |
| **Soil thermal** | Implicit backward-Euler heat diffusion on prognostic internal energy; freeze/thaw plateau implemented but **opt-in, default off**; bottom boundary adiabatic | Prognostic internal energy in the RK4 state, phase change always active | Two caveats on the record: MEDS's default column is shallower than the annual damping depth, so the annual wave reflects off a zero-flux base; and running with phase change off in a seasonally frozen site is wrong physics that the config permits. |
| **Snow / temporary surface water** | **Single bulk layer**; Niu–Yang cover fraction, accumulation, meltwater percolation, surface energy balance, snow-base conductance | Up to `NZS` layers (default 4), `IPERCOL` 0–2 | The clearest resolution downgrade in MEDS. |
| **Ground evaporation** | Dry-surface-layer resistance formulation | `IED_GRNDVAP` 0–5, default 0 (modified Lee–Pielke) | Different formulation. |
| **Interception** | Per-cohort Beer-law interception with a per-PAI storage capacity | `LEAF_MAXWHC` | Comparable. |
| **Soil biogeochemistry** | CENTURY reorganised as an explicit carbon matrix ODE: 7 pools (metabolic and structural litter × above/below, microbial, slow, passive) with a lignin sub-tracer; daily forward Euler on the fast loop's accumulated environmental scalar; exact matrix-exponential and a SASU steady-state solve for spin-up. **Opt-in, default off** | `DECOMP_SCHEME` 0–5, default 2 (CENTURY-like, 3 active pools); scheme 5 is the 5-pool Bolker CENTURY | MEDS's pool structure corresponds to ED2's `DECOMP_SCHEME = 5`. The matrix form is what makes the fast-loop heterotrophic respiration and the daily pool debit agree **by construction** — the seam closes to ~1e-13. The SASU steady-state solve has no ED2 equivalent and is a real spin-up accelerator. |
| **Nitrogen** | Scaffolded in the types, **not active** (carbon-only) | `N_PLANT_LIM`, `N_DECOMP_LIM` | |

---

## 3. Numerics

This is where MEDS diverges most deliberately, and where its documentation is most specific. The
reference page is [`docs/science/numerical_scheme.md`](science/numerical_scheme.md).

### 3.1 What is frozen — the part that is *the same*

Before any integrator runs, both models make the same semi-discretisation choice: photosynthesis,
stomatal conductance, dark respiration, canopy radiative transfer, plant hydraulic fluxes and
aerodynamic conductances are computed **once per fast step** and held fixed while the column is
advanced. ED2 calls `plant_hydro_driver` and `canopy_photosynthesis` immediately before entering its
integrator, and the RK stages only read the results; MEDS's `column_prepass` does the same thing for
the same reason. Freezing the ~20 s hydraulic mode is what makes any affordable step size possible.

MEDS states the consequence explicitly, and ED2 users should carry it too: **refining the time step
refines the freeze and the discretisation together**, and making the integrator more accurate at fixed
step size eventually stops helping, because the freeze error is the floor. On MEDS at 150 s, ~97% of
the remaining error is the freeze, not the integrator tolerance.

**One coefficient is deliberately not frozen.** The canopy-air ↔ atmosphere conductances are re-solved
at *every* integrator stage against that stage's own canopy-air state. ED2 does this too
(`update_diagnostic_vars` → `canopy_turbulence8` per stage). In MEDS it is load-bearing: freezing them
produced a sustained period-2 oscillation in canopy-air temperature of up to ~8 K step-to-step at
900 s **while every conservation budget closed to ~1e-6 J**. That is the most useful warning on this
page — *conservation is not stability*, and it applies to any model with this structure.

### 3.2 Integrators

| | MEDS v0.1.0 | ED2 (ED-2.2) |
|---|---|---|
| **Default** | `ark` — a 2-stage, second-order, **L-stable ESDIRK2** implicit scheme with an embedded error estimate and adaptive sub-stepping inside each `dt_fast` | `INTEGRATION_SCHEME = 1` — **fourth-order Runge–Kutta**, adaptive, `RK4_TOLERANCE = 0.01` |
| **Alternative** | `rk45` — adaptive explicit **Cash–Karp 5(4)** over the whole column state, kept as the accuracy baseline and the closest analogue of ED2's RK4 | `0` Euler (`NSUB_EULER`), `2` Heun, `3` hybrid (BDF2 implicit canopy + explicit rest) |
| **Retired** | An operator-split backward-Euler scheme was the historical default and was **removed** in 2026-07 | — |
| **Soil water** | Outside the ARK tableau (advanced once per step by the implicit Richards solver); inside the RK45 tableau | Inside the RK4 state |
| **Leaf ↔ canopy-air coupling** | Solved implicitly within every ARK stage by a direct 2×2 Newton | Explicit within the RK stages (except `INTEGRATION_SCHEME = 3`) |
| **Fallback** | When the explicit RK45 march cannot resolve a step (cold night, closed canopy), the driver discards it and redoes that `dt_fast` on `ark`. A counter reports how often | — |
| **Fast step** | `dt_fast`, required in TOML; shipped configs use **900 s** | `DTLSM`, default **600 s** |

The naming caveat, since it propagates: despite the config string `ark`, the scheme is **not IMEX** —
the biotic CO₂ source is folded implicit, the explicit tableau is empty, and it reduces to a clean
2-solve ESDIRK2.

### 3.3 What the step size costs, measured

Against a 12.5 s reference on a high-LAI sunlit stand over 24 h:

| `dt_fast` | 150 s | 300 s | 900 s |
|---|---|---|---|
| canopy-air T RMSE | 0.017 K | 0.055 K | 0.16 K |
| soil T RMSE | 0.06 K | 0.22 K | 0.70 K |
| GPP error (default config) | −0.0% | — | **−0.05%** |

Two things an ED2 user should take from this. First, `dt_fast` in MEDS is an **accuracy** parameter,
not a stability boundary — that changed when the conductances stopped being frozen. In the default
configuration the 900 s production step is cheap: daily GPP within 0.05% of a resolved reference, ET
within ~1%, canopy-air temperature within 0.16 K, with soil-surface temperature the loosest at 0.70 K.

Second, and independently of carbon, **leaf water potential itself does not converge in `dt_fast`** —
daytime mean −0.23 MPa at 12.5 s against −1.19 MPa at 900 s, because the plant water-mass update is an
explicit step with frozen sapflow and root uptake, so the per-step excursion grows with the step. Any
study keyed to leaf water potential — hydraulic stress, potential-driven mortality — should run at
≤ 150 s regardless of what the carbon budget looks like.

### 3.4 Conservation

Both models check budgets. ED2 verifies energy, water, CO₂ and carbon every photosynthesis step
against relative tolerances (`tol_subday_budget`, `tol_carbon_budget` in `utils/budget_utils.f90`).
MEDS closes whole-column water and energy ledgers to machine precision (~1e-13 kg m⁻² and ~1e-7 J m⁻²)
in the unsaturated regime, with snow and condensation active, and asserts this in the test suite; an
opt-in `[energy].debug_error` turns a non-closing budget into a hard stop.

Two hard-won qualifications, both stated in the MEDS docs and worth repeating to anyone porting
analysis code:

- A **whole-column** ledger cannot detect enthalpy placed in the wrong *layer* — the error is purely
  vertical and still sums correctly against the boundary. Three such defects were found in MEDS only
  by per-layer face checks or by a temperature that stopped being plausible.
- A closed budget says the bookkeeping is right and nothing at all about whether the trajectory is
  physical (see the 8 K oscillation in §3.1).

### 3.5 Parallelism and determinism

MEDS threads the patch axis with OpenMP (opt-in at both build and run time), and **the output is
byte-identical at any thread count** — site accumulators are staged per (sub-step, patch) and folded
back in patch order rather than by `reduction(+:)`, so the last bits do not drift with thread
scheduling. The hot demographic kernel carries OpenMP `target` regions and runs on GPU under
nvfortran. There is **no MPI**; ED2's polygon-parallel decomposition has no counterpart because MEDS
runs one site.

Measured on a 50-year forced spin-up: ~2.0× on 4 cores, which is ~67% of that machine's *measured*
4-core aggregate throughput (3.03×, not 4× — turbo and shared cache). Speedup is bounded by patch
count, so a bare-ground spin-up is the worst case.

### 3.6 Known numerical limitations

Carried verbatim from the MEDS docs, because they are the things that would bite you:

1. Leaf water potential does not converge in `dt_fast` (§3.3). Unfixed.
2. The stability threshold has been mapped on one synthetic forcing at one wind speed; the feedback
   involved is wind-dependent and light-to-moderate wind is the dangerous band.
3. The ARK error norm carries structurally-zero terms (soil moisture and plant water are split out of
   its tableau), so it runs slightly looser than its stated tolerance.
4. The whole evidence base in the numerics page is one site, one month per cell, one forcing year, no
   confidence intervals. Differences under ~2× should not be read as real.

---

## 4. Running them

### 4.1 Configuration

| | MEDS | ED2 |
|---|---|---|
| Format | **Two TOML files** — a main file (everything non-PFT) and a PFT file it names | `ED2IN` Fortran namelist + optional XML (`IEDCNFGF`) for parameters, + `EVENT_FILE` |
| Missing keys | **Hard error**, listing every missing key at once. There are no built-in defaults for model parameters | Defaults in `ed_params.f90`; XML overrides |
| Scheme selectors | Few; most processes have one implementation | ~40 (`I*` flags), many with legacy/beta branches |
| Derived parameters | `derive_config` / `derive_pft_rates`, overridable via an explicit `[derived]` block | `ed_params.f90` |
| Parameter provenance | The per-PFT table actually used is written to `<prefix>_pft_parameters.csv` every run | `ATTACH_METADATA` attaches metadata to HDF5 datasets |

The "no defaults" rule is the biggest practical difference in setting up a run. It buys you a config
that fully determines the simulation and cannot silently inherit a parameterisation you did not read.
It costs you the ability to start from a curated ED2-style PFT set — you must supply every trait.

### 4.2 Initialization

| | MEDS | ED2 |
|---|---|---|
| Modes | Bare ground, a cohort **census** (inventory-style), or a netCDF restart | `IED_INIT_MODE` −1 … 7: bare/near-bare, ED-1.0, ED-2.0, ED-2.1 HDF5 histories, multi-file nearest-neighbour matching, lidar-compatible variants |
| Soil state | From config | `ISOILSTATEINIT`, `SOILSTATE_DB`, or from the history file |
| Site properties | From config (latitude, elevation, reference height, soil column) | Map databases: `VEG_DATABASE`, `SOIL_DATABASE`, `SLCOL_DATABASE`, `SOILDEPTH_DB`, `LU_DATABASE`, `THSUMS_DATABASE` |

MEDS has no database-driven initialization at all. For a single site with known properties this is a
simplification; for anything regional it is a wall.

### 4.3 Meteorological forcing

| | MEDS | ED2 |
|---|---|---|
| Format | One netCDF file on a `(time, grid)` layout, or a constant-forcing backend with no file | `ED_MET_DRIVER_HEADER` pointing at HDF5 files, `IMETTYPE` |
| Gap filling | **None — a NaN in a required field is a hard error** | — |
| Shortwave partitioning | Erbs (default), Weiss–Norman, or a clearness-index scheme | `IMETRAD` 0–5 (as-is, SiB, Weiss–Norman, all-diffuse, all-direct, clearness index) |
| Sub-daily reconstruction | Interval-mean-conserving disaggregation anchored on the model sun (UTC + longitude + equation of time) | `IMETAVG` 0–3 declares the averaging convention |
| Multi-year cycling | Calendar recycling with Feb-29 reconciliation, day-of-year exact | `METCYC1` / `METCYCF`, `ISHUFFLE` |
| Preparation | `scripts/download_era5land.py` + `scripts/prep_era5land_forcing.py` | Community drivers in ED2 format |
| CO₂ | Constant from config, or from the forcing file; echoed into the output | `INITIAL_CO2` or from the driver |

The no-gap-filling policy is deliberate and will reject forcing files ED2 would happily run.

**A known bug, since you would hit it:** multi-decade runs on recycled met can get anti-phased
sub-daily shortwave. Single-year evaluation is unaffected. It is root-caused and unfixed.

### 4.4 Output

This is the part v0.1 rebuilt, and the largest usability difference.

| | MEDS | ED2 |
|---|---|---|
| Format | netCDF, written through the C library (no netCDF-Fortran dependency, so it builds under ifx and nvfortran) | HDF5 |
| Variable selection | **~200 variables, each switchable individually, per timescale.** Resolution order: registry defaults → axis toggles → group toggles → per-tier → per-variable overrides | A fixed schema per file type; the controls are the frequency flags (`IFOUTPUT`, `IDOUTPUT`, `IMOUTPUT`, `IQOUTPUT`, `IYOUTPUT`, `ITOUTPUT`, `IOOUTPUT`) and `IADD_{SITE,PATCH,COHORT}_MEANS` |
| Timescales | Four tiers: sub-daily (F), daily (D), monthly (M), annual (Y), each independently enabled and chunked | Per-file-type frequencies plus `FRQFAST` / `FRQSTATE` |
| Axes | scalar, cohort, patch, soil layer, **PFT**, **DBH size class**, and 2-D (patch × soil layer) | polygon / site / patch / cohort levels |
| Aggregation | Declared as data on each variable — a weight (none / plant density / leaf area / basal area / AGB) and mean-vs-sum. One registry line emits a field's patch, site, PFT and size-class rollups | Hard-coded per variable |
| Discovering variables | `meds_main --dump-io-config` writes a complete override file listing every variable | Read the source or the wiki |
| Empty sets | A mean over an empty patch/PFT/class is `_FillValue`; a sum is a true 0, so the PFT and size-class sums still equal the site total when a PFT goes locally extinct | — |
| Restart | A separate, deliberately orthogonal stream: `<prefix>-S-<YYYYMMDDHHMMSS>.nc`, raw prognostic state at an instant, never a time average | `ISOUTPUT`, `FRQSTATE`, HDF5 history files |

Two conventions worth knowing before you write analysis code. **Sub-daily variables are separate
registry entries from their coarse namesakes** (`cas_temp_fast` reads staging captured inside the fast
loop; `cas_temp_site` reads end-of-step state) — putting either on the wrong stream is rejected at
start-up rather than silently writing fill. And **cohort and patch axes are forbidden on the annual
stream**, because a window longer than a month would straddle the disturbance restructuring.

### 4.5 Build, test, scale

| | MEDS | ED2 |
|---|---|---|
| Build | CMake ≥ 3.20, automatic Fortran module dependency resolution | `make` with per-platform `include.mk`; the well-known "run make six times" pattern |
| Compilers | ifx and nvfortran verified each release; gfortran supported | gfortran, ifort, others |
| Dependencies | netCDF-C (mandatory) | HDF5, optionally MPI |
| Tests | 38 CTest targets, run on both compilers, Release and Debug | `EDTS` regression suite comparing whole-model output |
| Parallelism | OpenMP over patches (bit-identical), OpenMP `target` GPU offload | MPI over polygons + OpenMP within a rank |
| Scale | One site | Single points, regional grids, or coupled to BRAMS |

MEDS has no equivalent of EDTS. That is the gap behind the "not benchmarked" caveat at the top of
this page, and closing it is the obvious next piece of work.

---

## 5. Quick checklist: what an ED2 user will find missing

- Fire (`INCLUDE_FIRE`), land use and logging (`IANTH_DISTURB`, `SL_*`, `CL_*`)
- Nitrogen limitation (`N_PLANT_LIM`, `N_DECOMP_LIM`)
- Multiple sites per polygon, multiple polygons, regional grids, MPI
- Map-database initialization and soil texture classes / pedotransfer functions
- Prescribed phenology from files (`IPHEN_SCHEME = 1`)
- Frost mortality and explicit hydraulic-failure mortality
- Seed dispersal between patches (`REPRO_SCHEME = 2`)
- Big-leaf mode (`IBIGLEAF`), horizontal shading (`IHRZRAD`), finite crown radius (`CROWN_MOD`)
- Multi-layer snow (`NZS` > 1)
- The 17-PFT parameter set — MEDS ships no PFT defaults
- Any published benchmarking

## 6. Where to look

| topic | MEDS | ED2 |
|---|---|---|
| Time integration | `src/driver/meds_fast_{step,ark,rk45,time_derivs}.f90`; [`docs/science/numerical_scheme.md`](science/numerical_scheme.md) | `dynamics/rk4_*.f90`, `euler_driver.f90`, `heun_driver.f90`, `hybrid_driver.f90`, `bdf2_solver.f90` |
| Photosynthesis | `src/plant/meds_leaf_gas_exchange.f90`; [`leaf_gas_exchange.md`](science/leaf_gas_exchange.md) | `farq_leuning.f90`, `farq_katul.f90`, `photosyn_driv.f90` |
| Radiation | `src/biophysics/meds_canopy_radiation.f90`; [`canopy_radiation_transfer.md`](science/canopy_radiation_transfer.md) | `twostream_rad.f90`, `radiate_driver.f90` |
| Turbulence | `src/biophysics/meds_canopy_aerodynamics.f90`; [`canopy_aerodynamics.md`](science/canopy_aerodynamics.md) | `canopy_struct_dynamics.f90` |
| Hydraulics | `src/plant/meds_plant_hydraulics.f90`; [`plant_hydraulics.md`](science/plant_hydraulics.md) | `plant_hydro.f90` |
| Soil | `src/biophysics/meds_soil_{water,energy}.f90`; [`soil_biophysics.md`](science/soil_biophysics.md) | `rk4_misc.f90`, `lsm_hyd.f90`, `soil_coms.F90` |
| Allocation / growth | `src/plant/meds_plant_carbon_allocation.f90`; [`plant_carbon_allocation.md`](science/plant_carbon_allocation.md) | `growth_balive.f90`, `structural_growth.f90` |
| Decomposition | `src/biogeochemistry/meds_soil_biogeochem.f90` | `soil_respiration.f90`, `decomp_coms.f90` |
| Fusion / fission | `src/core/meds_core_{cohort,patch}_fusefiss.f90` | `fuse_fiss_utils.f90` |
| State | `src/core/meds_core_state_types.f90` | `memory/ed_state_vars.F90` |
| Output | `src/io/`; [`docs/science/diagnostics.md`](science/diagnostics.md) | `io/` |

Background reading for both: Moorcroft et al. 2001 (*Ecol. Monogr.*), Medvigy et al. 2009 (*JGR*),
and **Longo et al. 2019** (*GMD* 12:4309) — the ED-2.2 technical description, whose Supporting
Information is the reference MEDS was written against.
