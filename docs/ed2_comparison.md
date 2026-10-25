# MEDS v0.3.0, for people who use ED2

> **This page describes MEDS v0.3.0, and a rewrite against the current code is planned
> ([`ROADMAP.md`](ROADMAP.md) §12).** Until then, where the page disagrees with this list, the list
> is right:
>
> - **Water stress and resorption.** The non-stomatal capacity limb is always on
>   ([`leaf_gas_exchange.md`](science/leaf_gas_exchange.md) §4.1), and leaf and fine-root resorption
>   default to 0.3 ([`plant_carbon_allocation.md`](science/plant_carbon_allocation.md) §1).
> - **Leaf water potential converges at 900 s.** Soil and plant water are solved together, implicitly
>   ([`soil_plant_water.md`](science/soil_plant_water.md)); the midday leaf potential agrees between
>   75 and 900 s ([`numerical_scheme.md`](science/numerical_scheme.md) §5a), so a study of it needs
>   no shorter `dt_fast`. A drought season does: 900 s runs Palo Verde's driest months' GPP and ET
>   about 6 % high, and 300 s is the advice for dry-season leaf-habit timing (§5a). The canopy air is
>   stable at 900 s in light wind as in strong (§5a′, §7).
> - **Roots and sapwood.** Each cohort roots to a depth set by its height, and its sapwood's
>   conductivity follows its wood density ([`plant_hydraulics.md`](science/plant_hydraulics.md) §3–§4).
> - **Setup.** Height can follow Chave et al. (2014) ([`configuration.md`](configuration.md), "The PFT
>   file"); a soil column can come from texture through a pedotransfer function, and from the soil
>   maps ([`soil_biophysics.md`](science/soil_biophysics.md), "Hydraulic parameters from soil texture"
>   and "The column from the soil maps"); and a region runs its polygons on OpenMP threads, one polygon
>   per thread ([`configuration.md`](configuration.md), "Regional runs").
> - **Counts.** The suite has 86 tests (84 in a debug build) and the output registry 366 variables.

MEDS is a ground-up reimplementation of [ED2](https://github.com/EDmodel/ED2) in Fortran 2018. If you
already run ED2, this page tells you what carried over, what changed, what is missing, and what a MEDS
run costs you to set up.

**Reference points.** MEDS at tag `v0.3.0`; ED2 mainline at commit `125f814d` (ED-2.2). ED2
options are named by their `ED2IN` namelist key, and "ED2 default" means the value in the shipped
`ED/run/ED2IN`.

**What changed between releases**, and which numbers each change moved, is in
[`CHANGELOG.md`](../CHANGELOG.md); this page describes the code as it is.

**One thing to establish up front: MEDS has not been benchmarked against ED2.** No EDTS-equivalent
regression suite has been run, no site has been compared flux-for-flux, and no output of any kind has
been scored against observations. What has been verified is internal — 53 unit tests on two compilers,
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

- **Scope is much narrower.** One soil column per patch, no fire, no land use or harvest, no nitrogen,
  no MPI. A run is one site, or a region of independent polygons on the ERA5-Land grid (one site each,
  in one process, from bare ground).
- **Options are decisions, not switches.** ED2 exposes ~40 scheme selectors (`ICANRAD`, `IPHEN_SCHEME`,
  `IALLOM`, `DECOMP_SCHEME`, …), most with legacy and beta branches. MEDS mostly implements one path —
  usually ED2's default or its best-supported alternative — and deletes the rest. Where MEDS keeps a
  choice, it is because both branches are defensible.
- **Some defaults are more aggressive than ED2's.** Plant hydraulics is always on and unconditional
  (ED2 default `PLANT_HYDRO_SCHEME = 0`, no hydraulics). Soil and plant water are solved together as
  one implicit network, with van Genuchten retention by default (ED2 default is Campbell/Cosby inside
  its RK4). Phenology is
  unconditional and rate-based.
- **Some are more conservative.** Snow is a single bulk layer against ED2's multi-layer temporary
  surface water. (Soil freeze/thaw and soil biogeochemistry are not opt-in: the freeze/thaw plateau
  is inherent to a column prognostic in internal energy, with ice-aware conductivity and heat capacity
  always active, and soil biogeochemistry is on by default.)
- **Configuration and output are wholly rebuilt.** TOML instead of `ED2IN` + XML; netCDF instead of
  HDF5; and roughly 200 output variables individually switchable per timescale rather than a fixed
  schema gated by frequency flags.
- **The code is unit-tested and thread-deterministic.** 53 CTest targets on two compilers; output
  byte-identical at any OpenMP thread count.

**Who should look at it.** If you need fire, land use, a nitrogen cycle, regional runs, or the ED2
PFT parameterisations validated across biomes, stay on ED2 — MEDS does not have them. If you are doing
single-site process work on plant hydraulics, canopy biophysics, sub-daily carbon–water coupling, or
numerics, MEDS is a substantially smaller and more tractable codebase (86 source files against ED2's
`ed_state_vars.F90` alone being 34k lines), with the conservation and determinism properties spelled
out.

---

## 2. Processes

Read the notes column — that is where the qualifications live.

### 2.1 Structure and demography

| | MEDS v0.3.0 | ED2 (ED-2.2) | Notes |
|---|---|---|---|
| **State hierarchy** | polygon → site → patch → cohort, one site per polygon | grid → polygon → site → patch → cohort | MEDS drops the grid level: a region is N polygons on the ERA5-Land grid, each stepped on its own, and a site run is one. State is a flat site-wide Structure-of-Arrays with a CSR patch map, not nested ragged arrays. |
| **PFTs** | Run-time count; every trait supplied from a TOML file. No built-in table. | 17 hard-coded PFTs with defaults in `ed_params.f90`, overridable by XML. | MEDS has **no PFT defaults at all** — a missing key is a hard error. That is a deliberate trade: no hidden parameterisation, but no curated ED2 PFT set to inherit either. The shipped example uses three PFTs contrasted along wood density. |
| **Allometry** | Pan-tropical, ED2 `IALLOM = 3` | `IALLOM` 0–4, default 3 | Same family. Biomass is Chave et al. (2014) eq. 4 in carbon, where ED2 uses its own refit of it; leaf area is ED2's BAAD fit divided by `C2B`, as ED2 applies it. MEDS additionally inverts AGB→DBH for fusion. |
| **Cohort fusion / fission** | Keys on height similarity + an LAI cap; conserves per-plant AGB (DBH re-derived, never averaged) and plant number | `fuse_fiss_utils.f90`, `MAXCOHORT` target with tolerance relaxation | Same idea, different similarity metric. MEDS asserts the conservation invariants at 1% and `error stop`s otherwise. |
| **Patch fusion** | ED2-style cumulative-LAI light profile over height layers, same disturbance class only | `IFUSION` 0/1, `MAXPATCH` | MEDS ports the light-profile criterion. |
| **Treefall disturbance** | Annual; a fraction of every patch peels into one age-0 gap; tall die, short understory survives | `TREEFALL_DISTURBANCE_RATE`, plus small-tree background mortality | Directly ported, including the tall/short split. |
| **Fire** | **Absent** | `INCLUDE_FIRE` 0–3, `SM_FIRE` | |
| **Land use / harvest / plantation** | **Absent** | `IANTH_DISTURB`, the `SL_*` / `CL_*` logging block, agriculture and pasture stock PFTs | |
| **Cohort/patch identity** | Persistent `global_id`, never reused, survives every sort/fuse/compaction | Not tracked | Lets an external reader follow one cohort through the output until it fuses away. |

### 2.2 Plant physiology

| | MEDS v0.3.0 | ED2 (ED-2.2) | Notes |
|---|---|---|---|
| **Photosynthesis** | FvCB (C3) with Rubisco / RuBP-regeneration / TPU limitation; Collatz (C4). Arrhenius or peaked temperature response. C3 and C4 carry **separate co-limitation curvatures** (`theta_cj_c3` 0.98 / `theta_ip_c3` 0.95) | `IPHYSIOL` 0–3; default 2 (Collatz Q10 with Moorcroft high/low-T corrections). Options 1 and 3 add Jmax and TPU | MEDS is closest to `IPHYSIOL = 3` in structure (explicit Jmax + TPU) but exposes the temperature response as a per-trait choice rather than a bundled scheme number. C3's own curvatures raise leaf A **22.4%** and coupled GPP **18.4%** over sharing C4's. |
| **Thermal acclimation** | Kattge & Knorr 2007: entropy terms and the Jmax:Vcmax ratio track a 30-day running mean growth temperature. **Opt-in, default off** | `IPHYSIOL` has no acclimation; `TRAIT_PLASTICITY_SCHEME` acclimates to *light*, not temperature | Vcmax optimum moves 28.2 → 38.2 °C over a 20 K growth range. Coupled GPP **+42.8%** when enabled — and two thirds of that is the **Jmax:Vcmax ratio**, which K&K puts at 2.24 at Ithaca's growth temperature against the PFT file's fixed 1.7. Requires the peaked temperature response. |
| **Stomatal conductance** | Leuning 1995, Medlyn 2011 (USO), Katul 2010 optimality | `ISTOMATA_SCHEME` 0 (Leuning) or 1 (Katul) | Medlyn USO is the addition. |
| **A–Ci coupling** | One nonlinear equation in Ci, solved by bracketing + bisection; uniform across all three gs models and both pathways | Analytical solution per case in `farq_leuning.f90`; separate solver in `farq_katul.f90` | MEDS trades ED2's analytic speed for one code path that every model shares. |
| **Water stress on carbon** | Two limbs. Stomatal limb (Sabot 2022, driven by soil water potential) **on**; non-stomatal capacity limb (ramp on Vcmax/Jmax/TPU in leaf water potential) **off by default** | `H2O_PLANT_LIM` 0–5, default 2 (supply/demand FSW from rooting-zone water) | Different formulation entirely. The capacity limb is off because its two parameters are weakly constrained and it is rarely measured directly, and because leaf water potential is not converged at production `dt_fast` (§3.3) — wiring an unconverged potential into carbon through a poorly constrained ramp is not a trade worth making. The limb's formulation is itself under review. |
| **Plant hydraulics** | **Always on, unconditional.** Two stores per cohort (leaf, wood), nonlinear Bartlett pressure–volume capacitance, Kirchhoff-integrated xylem conductance, per-layer Gardner rhizosphere conductances, solved implicitly with the soil as one network | `PLANT_HYDRO_SCHEME` 0 (default, off — tissues always saturated), 1 or 2 (Xu 2016 / Christoffersen 2016) | The single biggest default difference. MEDS has no "no hydraulics" path. Hydraulic redistribution is deliberately **not** enabled (no uptake from a layer drier than the wood). Thirteen hydraulic traits are **per-PFT**; they are selectable but **uncalibrated**, so MEDS ships hydraulically identical PFTs. |
| **Trait plasticity** | Light-driven SLA / Vcmax25 / Rd25 / leaf turnover rate, turnover-limited; opt-in, default off | `TRAIT_PLASTICITY_SCHEME` 0–3 and −1/−2, default 0 | MEDS's single scheme corresponds most closely to ED2's option 3 (change constrained by leaf turnover). |
| **Phenology** | One generic cue kernel: daily cues → a smoothed flush tendency and a smoothed senescence tendency. Per-PFT **flush** and **shed** cue masks select which of three cues drive each side: temperature (warmth and cold degree-day sums), the hours of light (a running mean of the hours the PAR reaching each cohort's top exceeds the PFT's `par_min`; a low `par_min` counts the photoperiod) and water (predawn leaf potential summed against the turgor-loss point). Evergreen is a leaf-cover floor, not a flag | `IPHEN_SCHEME` −1 … 4 (evergreen / drought-deciduous old and new / prescribed / light / hydraulic) | A genuine restructuring, not a port. ED2's scheme numbers become per-PFT cue masks, so cold-deciduous, drought-deciduous and light-driven habits coexist in one run without a global switch. Prescribed phenology from files (`IPHEN_SCHEME = 1`, `PHENPATH`) has **no MEDS equivalent**. **Read the caveat in §2.2a before using any of it.** |
| **Carbon allocation** | FATES PARTEH-H1: four pools (leaf, fine root, wood, non-structural), allometric targets, daily priority ladder, wood as the residual sink and the prognostic size anchor | `growth_balive` (daily) + `structural_growth` (monthly), `ISTRUCT_GROWTH_SCHEME` 0/1 | Unified into one daily step. Growth respiration is charged on realized growth only. Leaf **resorption** on shed is available (a retained fraction returns to storage instead of becoming litter), **default 0** — the shipped `leaf_turnover_rate` is calibrated on observed litterfall, which already has resorption in it. |
| **Maintenance respiration** | Leaf (from the gas-exchange kernel), stem (Chambers 2004, surface-area based), fine root (summed **over layers**, so the temperature response is applied per layer rather than to a mean) | `GROWTH_RESP_SCHEME`, `STORAGE_RESP_SCHEME`, `ISTEM_RESPIRATION_SCHEME` (default 1 = Chambers) | Same stem formulation as ED2's current default. MEDS has no storage-pool maintenance (ED2's `growth_balive` turnover, `STORAGE_RESP_SCHEME`): at ED2's temperate-broadleaf 0.6243/yr an Ithaca run loses 35% GPP and 43% AGB, and no shipped PFT used it. |
| **Mortality** | One hazard: Camac 2018 additive `gamma + alpha·exp(−beta·growth)`, on a tracked moving-average carbon growth rate, with per-PFT coefficients (shipped as Camac's power laws in wood density, computed by `scripts/pft_traits`). Plus treefall. | Six additive pathways: ageing, carbon-starvation (`CARBON_MORTALITY_SCHEME`, `IDDMORT_SCHEME`, `CBR_SCHEME`), treefall background, cold/frost, hydraulic failure (`HYDRAULIC_MORTALITY_SCHEME`), and disturbance | **The largest science gap.** MEDS has no frost mortality, no explicit hydraulic-failure mortality, and no separate carbon-starvation term — growth-dependence carries all of it. ED2's Camac option (`CARBON_MORTALITY_SCHEME = 2`) is the closest analogue. The **carbon** leaving by each pathway (background hazard / cull at the tracking floor / treefall) is reported separately, which ED2 does not do. |
| **Recruitment** | Reproduction carbon → recruits **within the parent patch**, plus a per-PFT external seed rain; pooled until a minimum size, then spawned | `REPRO_SCHEME` 0–3, default 2 (seeds exchanged among all patches of a polygon) | MEDS does **not** disperse seeds between patches. In a strongly heterogeneous stand this matters. |

### 2.2a Two caveats on what is selectable

**Phenology is selectable and self-consistent. The coupled model is not validated.** The cue
kernel has been fitted offline at four sites (Harvard Forest, Hyytiälä, Barro Colorado Island, Palo
Verde; `examples/example02_canopy_phenology`), but no coupled MEDS run's leaf-area cycle has yet been
scored against an observation, and the water cue has been driven only by a soil-moisture surrogate
and by a MEDS run's leaf water potential without the shed-to-water feedback. Thresholds are site- and
biome-specific: a bright-hours senescence trigger set for a tropical dry season, used at a temperate
site, strips the canopy every summer.

**Thermal acclimation, storage-pool maintenance respiration, leaf resorption and the non-stomatal
water-stress limb all ship off.** Each moves numbers substantially; turning any of them on is a
decision to leave the shipped baseline.

### 2.3 Biophysics

| | MEDS v0.3.0 | ED2 (ED-2.2) | Notes |
|---|---|---|---|
| **Canopy radiative transfer** | Two-stream, multi-layer, one solver run per band (VIS, NIR, thermal LW); Beta leaf-angle distribution; absolute W m⁻² throughout | `ICANRAD` 0/1/2, default 2 (Liou two-stream) | Ported from `twostream_rad.f90`. MEDS generalises the band loop and drops ED2's normalise-to-unit-incidence convention. No horizontal shading (`IHRZRAD`), no finite crown radius (`CROWN_MOD`). |
| **Canopy turbulence** | CLM5 Monin–Obukhov surface layer + per-cohort exponential wind extinction + ED2 Nusselt leaf/wood boundary layers + CLM ground conductance | `ICANTURB` 0–4, default 2 (Massman 1997); `ISFCLYRM` 1–4, default 3 (Beljaars–Holtslag) | MEDS's combination is closest to ED2's `ICANTURB = 4` / `ISFCLYRM = 4` (the CLM-based branches), not to the ED2 defaults. |
| **Canopy air space** | Prognostic temperature, humidity **and CO₂** — three twins in one box model | Prognostic in the RK4 state, same three | Equivalent. MEDS's CAS depth is per-patch state derived from the tallest cohort plus freeboard. |
| **Leaf / wood energy** | Each cohort's leaf and wood energy [J per plant] is its state, and its temperatures are derived from it. Each step relaxes the temperatures **exactly** (the tissue ODE is linear under the frozen coefficients, so a closed-form exponential with an endpoint weight and a step-average weight) and commits the energy | `IBRANCH_THERMO` 0/1/2, default 1 (branch and leaf as one pool in the biophysics) | MEDS is closest to `IBRANCH_THERMO = 2`. As in ED2, the state is the energy (`leaf_energy`, `wood_energy`), and a slow step that changes a tissue's mass keeps its temperature. Unlike ED2, the energy is per plant, not per m² of ground; the canopy films are valued apart from it; and the tissue water stays liquid below 0 °C, so there is no liquid fraction. The wood's heat capacity is sized from the ED2 allometry (all the wood dry, the whole wood store's water wet), which gives it a time constant of ~20–30 min. |
| **Soil water** | Mixed-form Richards, solved by implicit Euler with Newton and adaptive sub-steps together with the pond and the plants' water (one network); thickness-weighted geometric-mean face K; no saturation clip (a specific storage above saturation). van Genuchten (default) or Campbell retention | Integrated inside the RK4 patch state; Campbell(-Mualem) with Cosby PTF by default (`SOIL_HYDRO_SCHEME` 0), Tomasella-Hodnett or van Genuchten as beta options 1/2 | Different numerical treatment (see §3) and a different default retention curve. |
| **Soil hydraulic parameters** | Given directly per column (saturated water content, Ksat, curve parameters); uniform texture broadcast over layers; no pedotransfer function, no texture classes | 17 texture classes, PTFs on sand/silt/clay (+ SOC, pH, CEC, bulk density for scheme 2), soil-texture and soil-depth map databases | **A real setup difference.** ED2 will build a soil column from a texture map; MEDS wants numbers. |
| **Soil bottom BC** | drainage $`q = f\,K`$ (free drain $f = 1$, a soil map's leaky base at its `f_leak`, bedrock $f = 0$) / aquifer, where the aquifer is a head-driven two-way boundary against a saturated zone at the column base | `ISOILBC` 0–3 (bedrock / free drainage / lateral drainage with `SLDRAIN` / aquifer) | No lateral-drainage-by-slope option in MEDS; ED2's `SLDRAIN` scales free drainage by the slope, MEDS's leaky base by the map's factor. |
| **Soil thermal** | Implicit backward-Euler heat diffusion on prognostic internal energy; the freeze/thaw plateau and ice-aware conductivity/heat capacity are **always active**; bottom boundary geothermal (default) or an optional **Dirichlet anchor** below the column | Prognostic internal energy in the RK4 state, phase change always active | Equivalent on the plateau. The caveat that stands: MEDS's default 2 m column is shallower than the annual damping depth, so the annual wave reflects off a near-zero-flux base — at Ithaca the annual swing at −1.73 m is **18.6 K in a 2 m column against 13.6 K in a 3 m one**, and 3 m is not converged either. The alternative to deepening is to anchor the base at a fixed temperature at a derived depth (`d/√2` below the deepest node, 3.12 m for the shipped column). Validated against the analytic semi-infinite solution `exp(−z/d)` — the adiabatic base holds **+82%** too much amplitude at the base node, the anchor **−2%**. Off by default, so the shipped configuration is unchanged. |
| **Snow / temporary surface water** | **Single bulk layer**; Niu–Yang cover fraction, accumulation, meltwater percolation, surface energy balance, snow-base conductance | Up to `NZS` layers (default 4), `IPERCOL` 0–2 | The clearest resolution downgrade in MEDS. |
| **Saturation vapour pressure** | Bolton 1980 over liquid, Magnus over ice, blended on a liquid fraction that crosses at exactly 0 °C | Murphy–Koop family in `therm_lib8` | The liquid branch agrees with Murphy–Koop to 0.2%. The ice branch matters: over months with snow it lowers latent heat **4.26%** against the liquid curve alone. |
| **Ground evaporation** | Dry-surface-layer resistance formulation | `IED_GRNDVAP` 0–5, default 0 (modified Lee–Pielke) | Different formulation. |
| **Interception** | Per-cohort Beer-law interception with a per-PAI storage capacity | `LEAF_MAXWHC` | Comparable. |
| **Soil biogeochemistry** | CENTURY reorganised as an explicit carbon matrix ODE: 7 pools (metabolic and structural litter × above/below, microbial, slow, passive) with a lignin sub-tracer; daily forward Euler on the fast loop's accumulated environmental scalar; exact matrix-exponential and a SASU steady-state solve for spin-up. **On by default** | `DECOMP_SCHEME` 0–5, default 2 (CENTURY-like, 3 active pools); scheme 5 is the 5-pool Bolker CENTURY | MEDS's pool structure corresponds to ED2's `DECOMP_SCHEME = 5`. The matrix form is what makes the fast-loop heterotrophic respiration and the daily pool debit agree **by construction** — the seam closes to ~1e-13. The SASU steady-state solve has no ED2 equivalent and is a real spin-up accelerator. |
| **Nitrogen** | **Not built** (carbon-only); planned after the move to ED-MORE | `N_PLANT_LIM`, `N_DECOMP_LIM` | |

---

## 3. Numerics

This is where MEDS diverges most deliberately, and where its documentation is most specific. The
reference page is [`docs/science/numerical_scheme.md`](science/numerical_scheme.md).

### 3.1 What is frozen — the part that is *the same*

Before any integrator runs, both models make the same semi-discretisation choice: photosynthesis,
stomatal conductance, dark respiration, canopy radiative transfer, plant hydraulic fluxes and
aerodynamic conductances are computed **once per fast step** and held fixed while the column is
advanced. ED2 calls `plant_hydro_driver` and `canopy_photosynthesis` immediately before entering its
integrator, and the RK stages only read the results; MEDS's `build_frozen` does the same thing for
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

| | MEDS v0.3.0 | ED2 (ED-2.2) |
|---|---|---|
| **Default** | `esdirk2` — a 2-stage, second-order, **L-stable ESDIRK2** implicit scheme with an embedded error estimate and adaptive sub-stepping inside each `dt_fast` | `INTEGRATION_SCHEME = 1` — **fourth-order Runge–Kutta**, adaptive, `RK4_TOLERANCE = 0.01` |
| **Alternative** | `rk45` — adaptive explicit **Cash–Karp 5(4)** over the whole column state, a run's explicit alternative and the closest analogue of ED2's RK4 | `0` Euler (`NSUB_EULER`), `2` Heun, `3` hybrid (BDF2 implicit canopy + explicit rest) |
| **Soil water** | Advanced once per step by the soil-plant water network, on both schemes; the stages move it at the network's time-mean rates | Inside the RK4 state |
| **Leaf ↔ canopy-air coupling** | Solved implicitly within every ESDIRK2 stage by a direct 2×2 Newton | Explicit within the RK stages (except `INTEGRATION_SCHEME = 3`) |
| **Fallback** | When the explicit RK45 march cannot resolve a step (cold night, closed canopy), the driver discards it and redoes that `dt_fast` on `esdirk2`. A counter reports how often | — |
| **Fast step** | `dt_fast`, required in TOML; shipped configs use **900 s** | `DTLSM`, default **600 s** |

The naming caveat, since it propagates: the ARS(2,2,2) tableau is usually presented as an IMEX pair,
but the scheme is **not IMEX** — the biotic CO₂ source is folded implicit, the explicit tableau is
empty, and it reduces to a clean 2-solve ESDIRK2, the config string's name for it.

### 3.3 What the step size costs, measured

Against a 12.5 s reference on a high-LAI sunlit stand over 24 h. **These convergence errors were
measured at v0.1.0**, whose absolute GPP differs from the current model's by ~18%; they are about how
far a coarse step sits from a resolved one, which that does not obviously change. Read them as the
shape of the `dt_fast` trade, not as current measurements:

| `dt_fast` | 150 s | 300 s | 900 s |
|---|---|---|---|
| canopy-air T RMSE | 0.017 K | 0.055 K | 0.16 K |
| soil T RMSE | 0.06 K | 0.22 K | 0.70 K |
| GPP error (default config) | −0.0% | — | **−0.05%** |

Two things an ED2 user should take from this. First, `dt_fast` in MEDS is an **accuracy** parameter,
not a stability boundary, because the conductances are re-solved at every stage (§3.1). In the default
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
opt-in `[fast].debug_error` turns a non-closing budget into a hard stop.

Two hard-won qualifications, both stated in the MEDS docs and worth repeating to anyone porting
analysis code:

- A **whole-column** ledger cannot detect enthalpy placed in the wrong *layer* — the error is purely
  vertical and still sums correctly against the boundary, and a temperature that stops being
  plausible is not a detector. `faces[soil_layer_mass]` is: it reports, per interior layer, the mass
  the soil-energy equation was charged for minus the mass the committed water actually moved.
  Reintroducing a wrong-layer defect leaves **both** whole-column ledgers at zero failures over
  103 680 checks while this number moves from 2.3e-13 to 1.1e-01 kg m⁻². If you maintain a model with
  this structure, that comparison is the most useful thing on this page.
- A closed budget says the bookkeeping is right and nothing at all about whether the trajectory is
  physical (see the 8 K oscillation in §3.1).

### 3.5 Parallelism and determinism

MEDS threads the patch axis with OpenMP (compiled in by default, asked for at run time), and **the output is
byte-identical at any thread count** — site accumulators are staged per (sub-step, patch) and folded
back in patch order rather than by `reduction(+:)`, so the last bits do not drift with thread
scheduling. There is no GPU build: OpenMP `target` offload measured 1.4× slower than the CPU
(`docs/dev_plans/archive/MEDS_GPU_EVALUATION.md`). There is **no MPI**: a region's polygons run one
after another in one process, and a larger domain is split into boxes run as separate jobs.

Measured on a 50-year forced spin-up: ~2.0× on 4 cores, which is ~67% of that machine's *measured*
4-core aggregate throughput (3.03×, not 4× — turbo and shared cache). Speedup is bounded by patch
count, so a bare-ground spin-up is the worst case.

### 3.6 Known numerical limitations

Carried verbatim from the MEDS docs, because they are the things that would bite you:

1. Leaf water potential converges in `dt_fast`: it is solved implicitly with the soil by the water
   network, and the midday probe agrees to 0.005 MPa between 75 and 900 s
   (`science/numerical_scheme.md` §5a). The wood potential and a restart's first day at 900 s are
   not yet scored on the network (see [`ROADMAP.md`](ROADMAP.md)).
2. The stability threshold has been mapped on one synthetic forcing at one wind speed; the feedback
   involved is wind-dependent and light-to-moderate wind is the dangerous band.
3. The ESDIRK2 error norm carries structurally-zero terms (soil moisture and plant water are split out of
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
| Derived parameters | `derive_config`; no PFT trait is derived from another in the model (fitted traits come from `scripts/pft_traits`) | `ed_params.f90` |
| Parameter provenance | The per-PFT table actually used is written to `<prefix>_pft_parameters.csv` every run | `ATTACH_METADATA` attaches metadata to HDF5 datasets |

The "no defaults" rule is the biggest practical difference in setting up a run. It buys you a config
that fully determines the simulation and cannot silently inherit a parameterisation you did not read.
It costs you the ability to start from a curated ED2-style PFT set — you must supply every trait.

### 4.2 Initialization

| | MEDS | ED2 |
|---|---|---|
| Modes | Bare ground, a cohort **census** (inventory-style), or a netCDF restart | `IED_INIT_MODE` −1 … 7: bare/near-bare, ED-1.0, ED-2.0, ED-2.1 HDF5 histories, multi-file nearest-neighbour matching, lidar-compatible variants |
| Soil state | From config | `ISOILSTATEINIT`, `SOILSTATE_DB`, or from the history file |
| Site properties | From config (latitude, elevation, reference height, soil column), or for each cell of a region from the `[maps]` maps on the forcing grid: soil depth, texture and chemistry, the bottom boundary, the deep temperature and Chave's E | Map databases: `VEG_DATABASE`, `SOIL_DATABASE`, `SLCOL_DATABASE`, `SOILDEPTH_DB`, `LU_DATABASE`, `THSUMS_DATABASE` |

MEDS reads no vegetation database: a run starts from bare ground, a census or a restart. A bare-ground
start's `[pft].nbg_init_density` is ED2's `init_density`: every patch gets one cohort of each PFT whose
density is above zero, that many plants per m², at `min_cohort_height`.

### 4.3 Meteorological forcing

| | MEDS | ED2 |
|---|---|---|
| Format | A MEDS forcing file of time series at a list of points (`"default"`), the global ED_ERA5land archive of one netCDF file per variable per month on the 0.1° grid (`"ERA5land"`, which a region requires), or a constant-forcing backend with no file | `ED_MET_DRIVER_HEADER` pointing at HDF5 files, `IMETTYPE` |
| Gap filling | **None — a NaN in a required field is a hard error** | — |
| Shortwave partitioning | Erbs (default), Weiss–Norman, or a clearness-index scheme. MEDS's Weiss–Norman uses the paper's 600 / 720 W m⁻² visible / NIR beam, where ED2 uses 0.43 / 0.57 of the solar constant | `IMETRAD` 0–5 (as-is, SiB, Weiss–Norman, all-diffuse, all-direct, clearness index) |
| Sub-daily reconstruction | Interval-mean-conserving disaggregation anchored on the model sun (UTC + longitude + equation of time) | `IMETAVG` 0–3 declares the averaging convention |
| Multi-year cycling | Calendar recycling with Feb-29 reconciliation, day-of-year exact | `METCYC1` / `METCYCF`, `ISHUFFLE` |
| Preparation | ERA5-Land tools in `scripts/prepare_era5/` | Community drivers in ED2 format |
| Downwelling longwave | From the file, or **synthesized** (Brutsaert clear-sky + an Idso-style cloud correction driven by a clearness index) | Required in the driver |
| CO₂ | A constant (`co2_const`), or a prescribed series (`co2_file`; the CMIP7 global annual means ship) read on model time, so it keeps rising under recycled met; never from the met file. Echoed into the output | `INITIAL_CO2` or from the driver |

The no-gap-filling policy is deliberate and will reject forcing files ED2 would happily run.

The synthesized longwave is a **fallback, not a substitute**: against ERA5-Land `strd` the clear-sky
term alone is −29.9 W m⁻² biased (RMSE 43.6 on a 312.3 mean), so the cloud correction is not optional,
and a run on synthesized LW sits 1.37 K cooler at the soil surface than the same run on the file's.
Use it when your forcing has no longwave, not in preference to one that does.

Recycled met keeps its sub-daily phase over multi-decade runs: the recycle window is declared in
config and validated at open, and `test_met_file` checks recycle phase over 29 years.

What the forcing code in `src/io/met` takes from ED2:

- **The sample a step reads** is ED2's `cgrid%metinput`, held in `met_forcing_t` and passed by the
  driver rather than kept as a global. A region's polygon keeps its cell's latitude and longitude, as
  an ED2 polygon does.
- **Wind** interpolates as a squared speed, conserving kinetic energy, as ED2's `vels` does.
- **The Weiss–Norman partition** follows ED2's `short_bdown_weissnorman` (`radiate_utils.f90`). MEDS
  keeps ED2's visible diffuse term, which scales the beam it subtracts by cos z once more than other
  implementations of WN85 eq. 3, and ED2's behaviour just above the grazing cutoff (cos z ≈
  0.0175–0.029), where the NIR diffuse stream can go slightly negative while the four streams still
  sum to the total ([forcing](science/forcing.md) §6).
- **The rain/snow split** is a simplified form of ED2's Jin (1999) ramp.
- **The grid match**: `great_circle_distance` ports ED2's `dist_gc` (`great_circle.f90`), and
  `nearest_grid_index` the argmin loop of `match_poly_grid`.
- **Recycling** maps a model year into the window as ED2's `year_use` substitution does, anchored on
  the window's declared start rather than on 1 January; 29 February becomes 28 February in a non-leap
  file year, as ED2's `read_ol_file` repeats 28 February.
- **The temperature lapse** is ED2's `calc_met_lapse`.

### 4.4 Output

This is the largest usability difference.

| | MEDS | ED2 |
|---|---|---|
| Format | netCDF, written through the C library (no netCDF-Fortran dependency, so it builds under ifx and nvfortran) | HDF5 |
| Statistics | Every variable's time **mean** over its record's period (`mean_<x>`), and a state's value at the period's end (`<x>`) | Means; `IQOUTPUT` adds a mean diurnal cycle |
| Variable selection | **~370 variables, each selectable per timescale, by name or pattern** (`"*_co" = false`, `"mean_*_py" = "M Y"`): registry defaults, then the `io_config` lines, the most specific winning | A fixed schema per file type; the controls are the frequency flags (`IFOUTPUT`, `IDOUTPUT`, `IMOUTPUT`, `IQOUTPUT`, `IYOUTPUT`, `ITOUTPUT`, `IOOUTPUT`) and `IADD_{SITE,PATCH,COHORT}_MEANS` |
| Timescales | Four tiers: sub-daily (F), daily (D), monthly (M), annual (Y), each independently enabled and chunked | Per-file-type frequencies plus `FRQFAST` / `FRQSTATE` |
| Axes | scalar, cohort, patch, soil layer, **PFT**, **DBH size class**, and 2-D (patch × soil layer) | polygon / site / patch / cohort levels |
| Aggregation | Declared as data on each variable — a weight (none / plant density / leaf area / basal area / AGB) and mean-vs-sum. One registry line emits a field's patch, site, PFT and size-class rollups | Hard-coded per variable |
| Discovering variables | `meds_main --dump-io-config` lists every variable by group, with its default tiers and the tiers it can be on | Read the source or the wiki |
| Empty sets | A mean over an empty patch/PFT/class is `_FillValue`; a sum is a true 0, so the PFT and size-class sums still equal the site total when a PFT goes locally extinct | — |
| Restart | A separate, deliberately orthogonal stream, `[output]`'s `state_*` keys, at any interval of days, months or years: `<prefix>-S-<YYYYMMDDHHMMSS>.nc`, raw prognostic state at an instant, never a time average | `ISOUTPUT`, `FRQSTATE`, HDF5 history files |
| Post-processing | `post_proc/` — site time series, stand structure, a 3-D landscape, a per-PFT / per-size-class plotter, and a worked evaluation notebook that checks the conservation identities **from the file** | Community scripts |

Two conventions worth knowing before you write analysis code. **Sub-daily variables are separate
registry entries from their coarse namesakes** (`mean_cas_temp_py` reads staging captured inside the fast
loop; `mean_cas_temp_py` reads end-of-step state) — putting either on the wrong stream is rejected at
start-up rather than silently writing fill. And **cohort and patch axes are forbidden on the annual
stream**, because a window longer than a month would straddle the disturbance restructuring.

What the output code in `src/io/output` takes from ED2:

- **The temporal reduction** follows ED2's `integrate_` / `normalize_` / `zero_` family: a variable is
  folded into its window each step, normalized when the window closes, and re-seeded
  (`integrate_scalar` / `integrate_slab`, `normalize_scalar` / `normalize_slab`, `reset_buffer`).
- **The DBH size classes** are ED2's: a class of plants, and a cohort lands whole in the class of its
  mean diameter, never split across classes, which makes `mean_nplant_size_py` a stem-density distribution
  comparable to a forest inventory.
- **The relative saturation** `soil_wetness`, (θ − θr)/(θs − θr), carries ED2's name.

### 4.5 Build, test, scale

| | MEDS | ED2 |
|---|---|---|
| Build | CMake ≥ 3.20, automatic Fortran module dependency resolution | `make` with per-platform `include.mk`; the well-known "run make six times" pattern |
| Compilers | ifx and gfortran verified for v0.3.0; nvfortran is the multicore back end | gfortran, ifort, others |
| Dependencies | netCDF-C (mandatory) | HDF5, optionally MPI |
| Tests | 53 CTest targets, run on both compilers, Release and Debug | `EDTS` regression suite comparing whole-model output |
| Parallelism | OpenMP over patches, bit-identical at any thread count. (No GPU build: offload measured 1.4× *slower* than the CPU.) | MPI over polygons + OpenMP within a rank |
| Scale | One site, or a region of polygons in one process | Single points, regional grids, or coupled to BRAMS |

MEDS has no equivalent of EDTS. That is the gap behind the "not benchmarked" caveat at the top of
this page, and closing it is the obvious next piece of work.

### 4.6 Model state

MEDS's state files, `src/state/meds_{cohort,patch,site}_state_vars.f90`, do the job of ED2's
`memory/ed_state_vars.F90` and are named after it: they define the cohort, patch and site types and
the routines that allocate, copy and reorder them. Three differences:

- **No global state.** ED2 holds the model's state in the module-level `edgrid_g`. MEDS's modules
  define the types and the routines on them only; the state is passed as an argument.
- **One list per operation.** ED2 copies or packs each cohort array in every routine that moves
  cohorts, which is where its "forgot to reallocate one array" bugs come from. MEDS lists every
  cohort field once in each of a few lockstep routines (allocate, free, grow, move, reorder, copy),
  all in one section order, and every sort, fusion and compaction goes through them.
- **Rates are not state.** ED2 keeps `ddbh_dt`, `dagb_dt` and `mort_rate` in the cohort type, so every
  copy carries them. MEDS computes the tendencies into a separate bundle each step and applies them
  at once; the output reads them from the slow diagnostic block.

The fast-to-slow maintenance-respiration accumulators (`leaf_resp_accum`, `stem_resp_accum`,
`root_resp_accum`) mirror ED2's `today_leaf_resp`, `today_stem_resp` and `today_root_resp`, kept
apart from gross GPP as ED2 keeps them.

### 4.7 Where the shared functions come from

The curves in `src/share/functions` state MEDS's own rationale in their comments and their science
pages; what they take from ED2 is listed here.

- **Allometry** (`meds_allometry_lib`). The pan-tropical curves are ED2's `iallom = 3` (`ed_params.f90`,
  `utils/allometry.f90`). MEDS ships one coefficient set for every PFT; a PFT file may give each PFT its
  own curve and coefficients. The leaf-area coefficient
  `lai_b1` = 0.4677/2 comes from ED2's `size2bl` leaf fit,
  bleaf = `c14f15_bl_xx(1)`/(SLA·C2B)·(D²H)^0.6410495: SLA·bleaf, the leaf area, does not depend on
  SLA, so MEDS uses it directly. The wood area index is ED2's `cpatch%wai` form, b1·D^b2, and the
  sapwood fraction ED2's `dbh2sf`.
- **Plant hydraulics** (`meds_plant_hydr_lib`). The sapwood conductivity against wood density is Xu et al.
  (2016) as ED2 codes it (`plant_hydro_scheme = 2`), with wood density held to ED2's 0.35–0.95 g cm⁻³; it is
  computed outside the model, by `scripts/pft_traits`.
  The root profile is ED2's cumulative exponential; over the rooting depth it sums to 1 − β, and MEDS
  divides by that so the fractions sum to one.
- **Soil water** (`meds_soil_hydr_lib`). The van Genuchten parameters follow
  Hodnett & Tomasella (2002) and the saturated conductivity Ottoni et al. (2019), as ED2's
  `soil_hydro_scheme = 2` uses them; the residual water content is floored as ED2 floors it; the pore
  connectivity −1 is ED2's ([soil_biophysics](science/soil_biophysics.md) eq. (2)).
- **Thermodynamics** (`meds_therm_lib`). The energy-to-temperature inverter is ED2's
  `uextcm2tl`/`uint2tl`. ED2's saturated soil conductivity counts liquid water only; MEDS's counts the
  ice in a frozen layer.
- **Calendar** (`meds_time`). ED2's daylength takes the polar-day branch at arg ≤ 1, which is wrong;
  MEDS takes it at arg ≤ −1. A month step clamps the day into the target month, as ED2's does.
- **Temperature response** (`meds_therm_lib`). The peaked form is Medlyn et al. (2002)'s, as FATES
  also uses it.

### 4.8 Where the configuration defaults come from

The reference config, `meds_config_main.toml`, gives each key's default and what it means; what the
defaults take from ED2 and the models ED2 builds on is listed here.

- **Snow** (`[snow]`). The aged-snow albedos, 0.518 visible and 0.435 near infrared, are ED2's, and
  `min_new_snow_mass` is ED2's `tiny_sfcwater_mass`. The liquid-holding fraction 0.10 is LEAF-3's 1:9
  ratio of liquid to ice.
- **Canopy aerodynamics** (`[aerodynamics]`). The roughness and displacement ratios are ED2's `vh2vr` and
  `vh2dh`; the Monin–Obukhov constants and the convective velocity scale `wc` are CLM5's, and the
  dense-canopy ground conductance (`cs_dense`, `gamma_g`) CLM4's
  ([canopy_aerodynamics](science/canopy_aerodynamics.md)). The flat-plate (leaf) and cylinder (wood)
  Nusselt coefficients keep ED2's and CLM's values in `aero_params_t`, with no key.
- **Soil carbon** (`[soil_carbon]`). Every default is ED2's (`decomp_coms.f90`, `ed_params.f90`), from
  its schemes 0–4: the decay rates `k_fast` … `k_pass` are `decay_rate_fsc`, `_stsc`, `_msc`, `_ssc`
  and `_psc` (ED2's `decay_rate_ssc` depends on the scheme; MEDS takes 0.2), the respired fractions
  `er_*` are ED2's `r_fsc`, `r_stsc_o`, `r_stsc_l`, `r_ssc`, `r_psc` and `r_msc_int`/`_slp` (respired is
  1 − transferred), and the capped temperature and water responses are ED2's. `"ed2"` and `"century5"`
  carry ED2's `DECOMP_SCHEME` codes, 0 and 5, and the default `"euler"` step is ED2's forward Euler
  ([soil_carbon](science/soil_carbon.md)).
- **PFT traits** (`meds_config_pft.toml`). `root_to_leaf_ratio` is ED2's `q`, and `nbg_init_density`
  ED2's `init_density`.
- **The soil column.** Depths are negative downward, interface elevations ≤ 0, as ED2's `slz`.

### 4.9 Where the fast loop's forms come from

The fast-loop kernels in `src/fast_dynamics` state MEDS's own rationale in their comments and their
science pages. What they take from ED2, beyond the frozen coefficients, the per-stage turbulence and
the integrators of §3.1–3.2 and the biophysics of §2.3, is listed here.

- **Leaf gas exchange** (`meds_leaf_gas_exchange`). ED2 sets the C4 initial slope `kp = klowco2·vm`, so
  it inherits Vcmax's temperature response, high-temperature deactivation included; MEDS gives `kp25`
  Vcmax's activation and deactivation set, as it does TPU. The leaf water potential that drives the
  stomata is diagnosed once per `dt_fast` from the step-start leaf water and held through the stages,
  as ED2 holds it through `DTLSM`.
- **Roots** (`meds_plant_hydraulics`). The rooting depth `root_depth_b1·H^root_depth_b2` is ED2's
  height allometry `IALLOM = 1`: 5 m for a 35 m tree.
- **Maintenance respiration** (`meds_plant_respiration`). The stem term is Chambers et al. (2004)'s
  surface-area form on the stem (a cylinder plus the branch area from the WAI), scaled by the
  aboveground structural fraction, as ED2 codes it; the fine-root term is ED2's per-`broot` form, a
  per-kgC baseline times the layer-summed temperature scale times the fine-root carbon. ED2 references
  both at 15 °C and MEDS at 25 °C, so a value taken from ED2 or Chambers is converted once, when the PFT
  parameters are set.
- **Leaf and wood boundary layers** (`meds_canopy_aerodynamics`). The Nusselt kernel and its
  coefficients are ED2's, as are the air's viscosity and thermal diffusivity, linear in temperature
  about 0 °C. ED2 forms the Reynolds number with the thermal diffusivity, $`u\ell/\alpha`$; MEDS
  uses the kinematic viscosity, $`u\ell/\nu`$, so its Re is 1.42 times ED2's and a laminar leaf's
  forced conductance 1.19 times.
- **The ground's coupling to the canopy air** (`meds_canopy_aerodynamics`). ED2's bare-ground conductance
  is $`u_*\,\mathrm{temp1}`$, the canopy air's own factor to the atmosphere, so bare ground crosses the
  surface layer twice. MEDS takes CLM5's soil coefficient, $`\kappa u_*/(0.13\,(z_{0g}u_*/\nu)^{0.45})`$,
  which spans only the layer between the soil and its roughness height, so bare ground's resistance to
  the reference height is CLM5's. A patch without cohorts has ED2's bare roughness, `soil_rough` 0.01 m.
- **Tissue and soil energy** (`meds_fast_integ_derivs`, `meds_fast_integ_frozen`, `meds_soil_heat`). The leaf
  pays the full vapour enthalpy of what it transpires and the soil the liquid internal energy of what
  the roots draw, each once, as ED2's `rk4_derivs` charges `tq2enthalpy(T_leaf)` and `uint_water`. The
  advected-enthalpy names are ED2's: `qwflux_wl` is the sap flow's enthalpy into the leaf,
  `qloss_frozen` is ED2's `qloss` (the root water's enthalpy leaving the soil), and `q_wood_net` is
  `qloss − qwflux_wl`. The soil water-enthalpy faces take the upward-positive sign of ED2's
  `qw_flux_g`.
- **Soil water** (`meds_soil_water`). The face conductivity is the thickness-weighted geometric mean,
  ED2's rule for soil layers.
- **Canopy air** (`meds_fast_integ_backward_euler`). Neither model clamps the canopy-air humidity at saturation:
  transient supersaturation is relaxed by a condensation sink, not cut off.
- **Carbon bookkeeping** (`meds_fast_state`). Gross GPP and the maintenance losses are summed apart
  over the day and netted by the daily allocator, as ED2 does.

### 4.10 Where the slow loop's forms come from

The slow-loop code in `src/slow_dynamics` states MEDS's own rationale in its comments and its science
pages. What it takes from ED2 is listed here.

- **The slow driver.** `meds_slow_dynamics` plays the part of ED2's `veg_dynamics_driver`. Leaf
  phenology runs inside the slow loop, as ED2's `phenology_driv` runs inside its vegetation dynamics.
- **Allocation** (`meds_plant_carbon_allocation`).
  - It covers ED2's `growth_balive` and `structural_growth` in one daily step.
  - It follows FATES's PARTEH hypothesis 1 ("allometrically guided, carbon only"): every pool has an
    allometric target, and the residual builds wood.
- **Soil carbon** (`meds_soil_carbon`).
  - The seven pools are ED2's CENTURY set, and the network is ED2's `update_C_and_N_pools`,
    reorganized as a carbon matrix.
  - The `euler` step is ED2's forward Euler.
  - The above-ground pools take the surface temperature and the below-ground ones the active layer's,
    as ED2's `A_decomp` and `B_decomp`.
- **Cohort and patch fusion and fission** (`meds_slow_update`).
  - They are ED2's `fuse_fiss_utils`; `fuse_cohorts` and `fuse_patches` keep ED2's names.
  - Patch fusion compares ED2's cumulative-LAI light profiles (`patch_pft_size_profile`).
- **Treefall disturbance** (`apply_patch_disturbance`). The survivorship is ED2's treefall: cohorts
  at or above the survival height die in the gap, and the understory survives into it.
- **Light plasticity** (`meds_plant_traits`). It ports ED2's `update_cohort_plastic_trait`
  (Lloyd et al. 2010).

### 4.11 Where the run loop's form comes from

- **ED2.** `edmain` calls `ed_driver`, which initializes the run and calls `ed_model`, one time loop
  at the integration step (`DTLSM`): each step the canopy radiation and the integrator advance every
  site, a new day calls `veg_dynamics_driver`, a new month reads the meteorological drivers, and the
  output is written at its intervals.
- **MEDS.** `meds_main` drives `meds_region`: `region_open` builds the run, `region_step_month`
  steps it a month at a time, and `region_finalize` closes it. Each slow step of each polygon is
  `polygon_step`, and `meds_stepper`'s `advance_one_step` runs the fast loop's sub-steps inside it,
  then the slow dynamics ([order_of_processes](science/order_of_processes.md) §1–§2). The loop is
  the slow step, with the fast loop nested in it, where ED2's is the fast step with the dynamics
  called on the calendar. The forcing is loaded, and the output written, between the steps of
  every polygon rather than inside them.

---

## 5. Quick checklist: what an ED2 user will find missing

- Fire (`INCLUDE_FIRE`), land use and logging (`IANTH_DISTURB`, `SL_*`, `CL_*`)
- Nitrogen limitation (`N_PLANT_LIM`, `N_DECOMP_LIM`)
- Multiple sites per polygon, MPI, restarts of regional runs
- Map-database initialization and soil texture classes / pedotransfer functions
- Prescribed phenology from files (`IPHEN_SCHEME = 1`)
- Frost mortality and explicit hydraulic-failure mortality
- Seed dispersal between patches (`REPRO_SCHEME = 2`)
- Big-leaf mode (`IBIGLEAF`), horizontal shading (`IHRZRAD`), finite crown radius (`CROWN_MOD`)
- Multi-layer snow (`NZS` > 1)
- The 17-PFT parameter set — MEDS ships no PFT defaults
- Any published benchmarking

## 5a. ED2 two-stream defects found during the port

Six defects in ED2's `icanrad = 2` canopy radiative transfer, found while porting it and verified
against the longwave sibling routine, `old_twostream_rad.f90`, `multiple_scatter.f90`, and the
CLM / Sellers 1985 / SCOPE references. **All six are reported upstream.** They are recorded here
because two of them are the reason MEDS's RT is structured the way it is.

Which ED2 scheme MEDS targets: `icanrad` is namelist-restricted to `{1,2}`, so the *old* two-stream
is unreachable dead code; `icanrad = 2` is what MEDS reimplements. That option forbids
`crown_mod = 1`, so the crown area index is identically 1 — MEDS implements that closed-crown path,
which is ED2's own supported configuration rather than a simplification.

| # | ED2 location | Severity | What | MEDS |
|---|---|---|---|---|
| **B1** | `twostream_rad.f90:652-665` | **high** | In `sw_two_stream`'s `diffuseloop`, `ipft` is never assigned in the loop — it keeps `pft(ncoh)` from the earlier `directloop`, so **every cohort's diffuse ω/β use the canopy-top cohort's PFT**. Area weights are per-cohort and the direct beam is fine; the LW solver and both siblings index correctly. Wrong diffuse PAR/NIR, sub-canopy light and albedo on any patch with ≥ 2 PFTs of differing optics. | **Impossible by construction.** MEDS gathers per-cohort optics into arrays, so there is no reusable scalar index to go stale. Regression-guarded by `test_canopy_radiation`'s multi-PFT case. |
| **B6** | `radiate_driver.f90:639` | high (narrow) | Bedrock branch self-assigns an uninitialised local (`albedo_damp_nir = albedo_damp_nir`); should be `albedo_soil_nir`. Bites on bedrock top-soil with snow or surface water. | Implemented correctly in `ground_optics`. |
| **B7** | `radiate_driver.f90:1103-1109, 1215-1221` | low | The longwave leaf/wood absorption split omits the `clumping_factor` that the SW split and the solver's own `elai`/`etai` weighting use. Total LW is conserved; only the leaf/wood partition shifts. | **Impossible by construction.** MEDS derives the split from the same clumping-corrected `elai`/`ewai` and `(1-ω)` the solver used — one formula, all bands. |
| **B5** | `radiate_driver.f90:1653-1664` | low (diagnostic) | Fast-mean albedo is accumulated with no daytime guard, though its daily-mean sibling has one, so night-time soil albedo biases it. | Avoided. MEDS emits top-of-canopy incident and upwelling fluxes per band rather than a time-averaged albedo, so the reader forms `Σ up / Σ down` over whatever period they choose — which is what a satellite product is, and has no night-time 0/0. |
| **B2/B3** | `twostream_rad.f90:590, 164` | low (latent) | SW diffuse `μ` drops the `cai` factor that LW `μ` and the beam `μ0` carry; and a contested `cai`-weighted LW blackbody source. Reachable only at `cai < 1`, which `icanrad = 2` forbids — **zero impact in any valid ED2 run**. | Moot at `cai = 1`; the SCOPE formulation supersedes both expressions. Revisit only if MEDS adds finite crowns. |

**Correctly rejected as non-bugs** during verification, recorded so they are not re-investigated:
the `etai = 0` division (resolvable cohorts have positive area), the `rlong` zero-guard (downwelling
LW is always > 0 at the surface), the "unused" `fvis_*_def` fallback fractions, and the
backscatter / `μ̄` singularities (reachable only through out-of-range XML overrides). MEDS validates
inputs at construction to keep these unreachable.

### How the MEDS RT maps onto ED2's

MEDS keeps ED2's multi-layer vertical structure — one radiation layer per cohort, the
interface-continuity two-stream — and replaces the SW/LW duplication and the CLM leaf-angle optics.

| MEDS | ED2 | Mapping |
|---|---|---|
| `meds_canopy_radiation` (`canopy_radiation`) | `radiate_driver.f90 :: sfcrad_ed` | Per-patch orchestration: per-cohort optics, band loop, leaf/wood split, absorbed radiation + albedo. State-free — arrays in, `rad_out_t` out — rather than mutating `csite`. |
| `solve_band` / `layer_rt` | `sw_two_stream` **and** `lw_two_stream` | **Unified.** ED2's two ~600-line routines collapse to one: SW is beam with no emission, LW is emission with no beam, and the diffuse operator is identical. Absolute W m⁻² throughout, dropping ED2's normalise-to-unit-incidence. O(N) adding method in place of the dense `lisys_solver8`. |
| `meds_pft_derive` (`derive_rad_optics`) | `init_can_rad_params` + the scatter/backscatter block | ED2/CLM's `(1+χ)²` backscatter fit → **SCOPE/4SAIL**: Beta leaf-angle distribution → `bf`, exact Ross `G(μ)`, `β = ½(1 + bf(ρ−τ)/ω)`. Fixes the near-horizontal-leaf inaccuracy of the fit. |
| `ground_optics` | `sfcrad_ed` ground-albedo block | Per-band albedo + thermal emission, implemented B6-correct. |
| `meds_canopy_radiation` (its argument records) | `canopy_radiation_coms.f90` | ED2's mutable module globals → explicit derived types passed as arguments. |

---

## 6. Where to look

| topic | MEDS | ED2 |
|---|---|---|
| Time integration | `src/fast_dynamics/integrator/meds_fast_integ_{stepper,esdirk2,rk45,derivs}.f90`; [`docs/science/numerical_scheme.md`](science/numerical_scheme.md) | `dynamics/rk4_*.f90`, `euler_driver.f90`, `heun_driver.f90`, `hybrid_driver.f90`, `bdf2_solver.f90` |
| Photosynthesis | `src/fast_dynamics/processes/meds_leaf_gas_exchange.f90`; [`leaf_gas_exchange.md`](science/leaf_gas_exchange.md) | `farq_leuning.f90`, `farq_katul.f90`, `photosyn_driv.f90` |
| Radiation | `src/fast_dynamics/processes/meds_canopy_radiation.f90`; [`canopy_radiation_transfer.md`](science/canopy_radiation_transfer.md) | `twostream_rad.f90`, `radiate_driver.f90` |
| Turbulence | `src/fast_dynamics/processes/meds_canopy_aerodynamics.f90`; [`canopy_aerodynamics.md`](science/canopy_aerodynamics.md) | `canopy_struct_dynamics.f90` |
| Hydraulics | `src/fast_dynamics/processes/meds_plant_hydraulics.f90`; [`plant_hydraulics.md`](science/plant_hydraulics.md) | `plant_hydro.f90` |
| Soil | `src/fast_dynamics/meds_soil_{water,energy}.f90`; [`soil_biophysics.md`](science/soil_biophysics.md) | `rk4_misc.f90`, `lsm_hyd.f90`, `soil_coms.F90` |
| Allocation / growth | `src/slow_dynamics/processes/meds_plant_carbon_allocation.f90`; [`plant_carbon_allocation.md`](science/plant_carbon_allocation.md) | `growth_balive.f90`, `structural_growth.f90` |
| Decomposition | `src/slow_dynamics/processes/meds_soil_carbon.f90`; [`soil_carbon.md`](science/soil_carbon.md) | `soil_respiration.f90`, `decomp_coms.f90` |
| Respiration | `src/fast_dynamics/processes/meds_plant_respiration.f90`; [`plant_respiration.md`](science/plant_respiration.md) | `stem_resp_driv.f90` |
| Phenology | `src/slow_dynamics/processes/meds_plant_phenology.f90`; [`plant_phenology.md`](science/plant_phenology.md) | `phenology_driv.f90`, `phenology_aux.f90` |
| Mortality / recruitment | `src/slow_dynamics/processes/meds_demography.f90`; [`demography.md`](science/demography.md) | `mortality.f90`, `reproduction.f90` |
| Disturbance | `src/slow_dynamics/driver/meds_slow_update.f90` (`apply_patch_disturbance`) | `disturbance.f90`, `fire.f90`, `forestry.f90` |
| Forcing | `src/io/met/`; [`forcing.md`](science/forcing.md) | `ed_met_driver.f90` |
| Fusion / fission | `src/slow_dynamics/driver/meds_slow_update.f90` | `fuse_fiss_utils.f90` |
| Parameters | `src/config/` (TOML, all required); [`docs/configuration.md`](configuration.md) | `ed_params.f90`, the `*_coms` modules, `ED2IN` |
| State | `src/state/meds_{cohort,patch,site}_state_vars.f90`, named after ED2's; MEDS has no global like `edgrid_g` | `memory/ed_state_vars.F90` |
| Output | `src/io/output/`; [`docs/science/diagnostics.md`](science/diagnostics.md) | `io/` |

Background reading for both: Moorcroft et al. 2001 (*Ecol. Monogr.*), Medvigy et al. 2009 (*JGR*),
and **Longo et al. 2019** (*GMD* 12:4309) — the ED-2.2 technical description, whose Supporting
Information is the reference MEDS was written against.
