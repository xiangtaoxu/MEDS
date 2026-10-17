# MEDS roadmap

**What this file is.** Every piece of work MEDS has deliberately deferred, in one place. It
replaces the "reserved follow-ups" sections and the scattered `deferred` / `MVP` / `placeholder`
comments that used to be the only record of intent — those were findable only by grepping, and
several of them had quietly become true while nobody updated the comment.

**The rule.** A source comment states present-tense rationale. What changed and when goes in
[`CHANGELOG.md`](../CHANGELOG.md). What is deferred goes here, **and every item carries a GitHub
issue number** — so a deferred decision is findable, assignable and closable rather than a sentence
in a document nobody greps.

Items #150–#200 were filed on 2026-09-13 from this file. When you add an item here, open the issue
in the same change.

**Status vocabulary.** *Planned* — intended, with a design. *Candidate* — worth doing, no design
yet. *Open question* — the decision itself is undecided. Nothing here is a commitment to a date.

---

## 1. Known defects and open issues

The issues that predate this roadmap. Everything in sections 2–12 has an issue too (#150–#201,
filed 2026-09-13); these are listed separately because they were filed from measurement or from a
failure, not from a plan.

| # | Title | Note |
|---|---|---|
| [#74](https://github.com/xiangtaoxu/MEDS/issues/74) | Condensate is deposited into soil layer 1, not onto leaf/wood surface water | |
| [#96](https://github.com/xiangtaoxu/MEDS/issues/96) | Dynamic vapour pressure for leaf transpiration (Kelvin $`e_i`$) | Built, measured, removed; likely route to foliar water uptake |
| [#254](https://github.com/xiangtaoxu/MEDS/issues/254) | Passive deep **thermal** layers below the hydrologically active column | The Dirichlet anchor shipped and cut the base-layer amplitude error from +82 % to −2 %, but a purely resistive termination cannot reflect less than 0.41 — closing the rest needs heat *capacity* below the column, i.e. a thermal grid that extends past the water grid |
| [#265](https://github.com/xiangtaoxu/MEDS/issues/265) | Sub-canopy conductance: MEDS ports ED2's non-default `icanturb = 4`, giving an 8–16× too-stiff ground resistance | Matters in gaps |
| [#268](https://github.com/xiangtaoxu/MEDS/issues/268) | No litter layer: no surface organic horizon for the ground energy balance or soil evaporation to act on | |
| [#269](https://github.com/xiangtaoxu/MEDS/issues/269) | The canopy air space is one well-mixed slab: 1.2 K warmer than the free air at midday, where a real sub-canopy is cooler and steadier | |
| [#393](https://github.com/xiangtaoxu/MEDS/issues/393) | Leafless cohorts with no storage never die (9 % of BCI cohorts) | Nothing in the mortality law sees a cohort's carbon state, so a cohort that has lost its leaves and its storage persists while its leaf carbon decays geometrically. A carbon deficit destroys no tissue either: it is netted from the CO₂ the cohort owes (`docs/science/plant_carbon_allocation.md` §2). |
| [#394](https://github.com/xiangtaoxu/MEDS/issues/394) | Pool floors clamp overdraws instead of limiting the outflow | The slow applier (`update_cohort_states_kernel`) and `reacclimate_plant_traits` keep a pool non-negative with `max(…, 0)`, which creates the carbon an overdraw would have taken; the per-step bounds guard cannot see it. |
| [#395](https://github.com/xiangtaoxu/MEDS/issues/395) | `reconcile_tissue_water_capacity` creates and destroys tissue water without booking it | It computes the seeded and the discarded water (`seeded`, `discarded`), but no caller books them. |
| [#416](https://github.com/xiangtaoxu/MEDS/issues/416) | `grav_head` (9.804×10⁻³ MPa/m) is not `rho_h2o·grav` (9.80665×10⁻³) | The plant hydraulics and the pedotransfer functions use one value, the soil water the other; 0.03 %. |
| [#423](https://github.com/xiangtaoxu/MEDS/issues/423) | The `post_proc` readers use example 03's stand-layout names, but the README says to run them on `meds_main` output | Four scripts; deferred by the owner. Meanwhile `post_proc/README.md` says which layout each script reads. |
| [#390](https://github.com/xiangtaoxu/MEDS/issues/390) | Dry-soil conductivity: a film-flow term for the van Genuchten–Mualem curve | Enhancement. Until then the pore connectivity L = −1 (Schaap & Leij 2000) holds K above ED2's floor of 1.16×10⁻¹³ m s⁻¹ to about −2 MPa for a soil of n near 1.4 |

---

## 2. Phenology

Source: `docs/dev_plans/MEDS_PHENOLOGY_SENESCENCE_PLAN.md`. Science page:
[`science/plant_phenology.md`](science/plant_phenology.md).

One generic cue model (temperature, hours of light, water) with evergreen as a leaf-cover floor.
The kernel is fitted offline at three sites (`examples/example02_canopy_phenology`); no coupled MEDS
run's leaf-area cycle has been scored against an observation yet. That is the benchmarking gap in
§1, of which this is one instance.

- **The water cue on the model's own predawn potential** — *Planned.* A dry-forest run (Guanacaste)
  with the cue on, since the cue needs the shed → ψ feedback; example02 drives it with a
  soil-moisture surrogate.
- **`leaf_lifespan_toc` → `leaf_turnover_rate_toc`** — *Planned.* [#379](https://github.com/xiangtaoxu/MEDS/issues/379) The background turnover rate is
  the model's parameter; leaf lifespan is emergent (plan D2).

- **Drought-driven fine-root dieback** — *Candidate.* Fine roots age on the leaves' turnover gate and
  share their resorption (#387), but they are not shed with the leaves: root length sets the
  root–soil conductance, so a drop would cut a leafless plant's water uptake. A root dieback would act
  per soil layer, on that layer's water, not on the leaves.

---

## 3. Soil biogeochemistry

Source: `docs/dev_plans/archive/MEDS_BIOGEOCHEMISTRY_DESIGN.md` §7. Science page:
[`science/soil_carbon.md`](science/soil_carbon.md).

- **The nitrogen twin.** *Planned.* [#154](https://github.com/xiangtaoxu/MEDS/issues/154) Shaped in already: `n_cycle_on` is parsed and the N fields
  exist on `soil_carbon_t`, but no kernel reads the flag and the restart skips the fields.
  Needs the decomposition N limitation (`f_decomp`), the mineralization/immobilization flux, and
  N limitation of NPP on the plant side.
- **Vertically resolved soil-carbon pools.** *Candidate.* [#155](https://github.com/xiangtaoxu/MEDS/issues/155) The `n = 1` layer is the degenerate
  case of the intended profile.
- **An explicit coarse woody debris pool.** *Candidate.* [#156](https://github.com/xiangtaoxu/MEDS/issues/156) Today treefall necromass enters the
  structural pools directly.
- **Fire consumption of the ground pools.** *Candidate.* [#157](https://github.com/xiangtaoxu/MEDS/issues/157) No fire anywhere in MEDS.
- **`[soil_carbon]` keys that are read but drive nothing.** *Candidate.* [#418](https://github.com/xiangtaoxu/MEDS/issues/418) `agf_fast` and
  `agf_struct` are never used (the litter is partitioned by tissue), and the nitrogen keys wait for
  #154; refuse them or document them as inert.

---

## 4. Fast-loop numerics

Source: `docs/dev_plans/MEDS_PRODUCTION_INTEGRATOR_PLAN.md` §5–§8. Science page:
[`science/numerical_scheme.md`](science/numerical_scheme.md).

- **Soil and plant water as one network.** *Done 2026-10-12.* `MEDS_SPAC_WATER_NETWORK_PLAN.md`;
  science page [`science/soil_plant_water.md`](science/soil_plant_water.md). The pond, the soil layers
  and every cohort's wood and leaves are solved together, implicitly, on both schemes; the soil-water
  solver, the plant-hydraulics solver, the transpiration corrector, the tissue-water floor and RK45's
  saturation clip are gone. What is left:
  - **Newton work.** A step takes 3.5 Newton iterations per sub-step on BCI, mostly the small
    cohorts' wood and leaves at each change in demand, and partly root uptake switching on and off
    across layers (widening its 10⁻³ MPa ramp to 5×10⁻² cuts the iterations by 14 %). A warm start
    (persisted unknowns per patch) or eliminating each cohort nonlinearly are the candidates.
  - **The C API.** The network's process functions and solver take prescribed-potential nodes and
    fixed-flux boundaries, so a plant-only or soil-only solve is a configuration; exposing it through
    the C API is open.
  - **Sensible heat.** On BCI's dry season at 150 s the network's H is 0.8 W m⁻² above the old
    scheme's while LE agrees to 0.1 W m⁻²; not yet attributed.

- **N5 — an adaptive freeze cadence with a real error estimator.** *Closed 2026-09-30; revisit with
  the numerical scheme.* [#158](https://github.com/xiangtaoxu/MEDS/issues/158) It was proposed when the
  stability limit forced `dt_fast = 150 s` for the whole day. N2a removed that limit and the shipped
  `dt_fast` is 900 s, so what is left is a bounded efficiency gain: the coefficient pre-pass is
  4.4–11 % of a step. Any cadence would have to stay site-wide, because the threaded patch loop shares
  the sub-step's forcing samples and output staging.
- **Fold soil water into the ARK tableau.** *Closed 2026-09-30; revisit with the numerical scheme.*
  [#159](https://github.com/xiangtaoxu/MEDS/issues/159) — successor to #93. In-stage soil water costs
  +14–26 % and removes a split error of 0.011–0.018 % of column water. The committed soil state and
  the fluxes the soil-heat stages use come from one solve, the water network, on both schemes, so the
  borrowed-flux defect class is absent. The drought result behind closing #93 (2.3 % at 900 s) was measured under
  `ARREST_GS_CLAMP`, since replaced (#335) and then removed in favour of the always-on capacity
  limb on the current leaf potential; re-measure it in that round.
- **The `rwc_floor` clamp artefact** ([#104](https://github.com/xiangtaoxu/MEDS/issues/104)). *Closed
  2026-09-30, no change.* A wood store below its apoplastic floor maps to a potential of about
  −10⁴ MPa, and the hydraulics then take many sub-steps to refill it. Measured on the current scheme
  (`dev_plans/MEDS_EFFICIENCY_SWEEP_PLAN.md`, Phase 3 status), it hardly happens any more:
  - five calibrated BCI years never reach the floor;
  - the line-searched parameter set reaches it 3 times in a million cohort-steps and never empties
    a store;
  - only a stand with no rain at all collapses.

  Every entry happens on a step whose uptake the soil throttled, the path to revisit if it ever
  matters. Leaf and wood already share one pressure-volume curve (`water_curve_t`), so a
  physiology change would add a special case rather than remove one.
- **`psi_wood` at 900 s, and a restart's first day.** *Re-measure on the water network.*
  [#162](https://github.com/xiangtaoxu/MEDS/issues/162) `psi_leaf` converged at 900 s on ARK with the
  transpiration corrector (#91): over July on the established Ithaca stand, daily means at 900 s
  matched a 12.5 s run to 0.001 MPa on every day after the first,
  with a start-up transient on the first day after a restart (−1.83 against −0.28 MPa) and 0.17 MPa
  left in `psi_wood` on the 3-hour midday probe through the frozen uptake seam. The water network
  (2026-10) now solves leaf and wood water implicitly with the soil, which removes that seam; the
  midday probe's `psi_leaf` agrees to 0.005 MPa between 75 and 900 s. `psi_wood` and the restart
  transient are still to be re-measured.
- **Is `clamp_soil_energy` still needed?** *Open question.* [#420](https://github.com/xiangtaoxu/MEDS/issues/420) Its reason, ground evaporation's
  `pow()` inside the RK45 stages, is gone (evaporation is frozen from the water network); it moves energy
  with no ledger entry when it bites. Measure how often it fires, then remove it or restate why.

---

## 5. Vegetation energy

Source: `docs/dev_plans/archive/MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` §6–§7. Science page:
[`science/vegetation_energy_dynamics.md`](science/vegetation_energy_dynamics.md).

- **A separate canopy film store with phase change.** *Planned.* [#165](https://github.com/xiangtaoxu/MEDS/issues/165) Intercepted water currently
  has no independent thermal state and cannot freeze.
- **The free-convection slope.** *Closed 2026-09-30; revisit with the numerical scheme.*
  [#167](https://github.com/xiangtaoxu/MEDS/issues/167) The design note says the true sensible-heat
  slope is `1.25·h` and the solved `ΔT_leaf` is overstated ~20 % in calm conditions. Measured, it is
  not. `1.25` is the **pure free-convection** limit: `H ∝ ΔT^{1+m}` holds only for the Grashof part
  of the Nusselt number, so the real slope factor is `1 + m·f_free` with
  `f_free = Nu_free/(Nu_forced+Nu_free)`. At the shipped `leaf_width = 0.04 m`, and at the
  `ugbmin = 0.25 m/s` in-canopy wind **floor** — the most free-convection-favourable state the model
  can reach — that factor is **1.02–1.07** over the leaf-minus-CAS temperature range an Ithaca run
  actually produces, and only **1.15** at an extreme 15 K. It reaches 1.25 only at exactly zero
  wind, which `ugbmin` makes unreachable.

  The effect on `ΔT_leaf` is smaller again, because `h_coeff` is one of four terms in
  `veg_energy_balance`'s denominator (`h_coeff + le_slope + le_slope_wet + lw_slope`).

  And the fix is not the one the note describes. A flat factor on the denominator alone would break
  the kernel's exact conservation identity (`store_hcap_per_dt*(dt_end − dt_prev) + denom*dt_avg ==
  numer`). The correct form is a **Newton linearization about the frozen `ΔT₀`**, which changes the
  numerator, the denominator, the reported flux `dh`, *and* the exponential relaxation time constant
  `τ = cap/denom` together — in the kernel with the most delicate conservation invariant in the
  model, whose adaptive controller has already been broken once by a discontinuity here. That is not
  a proportionate trade for a few percent of one denominator term.
- **`veg_coupling_floor` is in ground-area units.** *Candidate.* [#349](https://github.com/xiangtaoxu/MEDS/issues/349) `veg_energy_balance`
  floors a tissue's coupling conductance at 1 W m⁻² K⁻¹ of ground, while every term it floors scales
  with the tissue's own area index, so a thin cohort's leaf or wood temperature is slaved to the
  canopy air. The deficit goes to the canopy air, so energy closes.
- **The wood's heat capacity uses a fixed sapwood moisture.** *Candidate.* [#419](https://github.com/xiangtaoxu/MEDS/issues/419) The water
  network carries the wood's own water, and the leaf's capacity already uses its prognostic water; the
  wood's does not.

---

## 6. Output and diagnostics

Source: `docs/dev_plans/MEDS_IO_V01_PLAN.md` §4.5–§4.6, §6; `MEDS_IO_DESIGN.md` §3.5. Science
page: [`science/diagnostics.md`](science/diagnostics.md).

- **A spectrally resolved surface radiative record.** *Candidate.* [#255](https://github.com/xiangtaoxu/MEDS/issues/255)
  Follow-up to #171: Per-band incident and upwelling fluxes
  shipped in v0.2.0 on the VIS/NIR/LW three-band grid the two-stream solves. Comparing against a
  multispectral product (MODIS bands, Sentinel-2) needs finer bands, which is a change to the RT's
  band structure rather than to its output.
- **A patch axis on the FAST tier.** *Done in beta.* [#270](https://github.com/xiangtaoxu/MEDS/issues/270) Each FAST quantity of the patch
  block has a `*_patch_fast` twin, staged from the patch block's own row.
- **Variances.** *Removed; add back if needed.* [#275](https://github.com/xiangtaoxu/MEDS/issues/275)
  - **Done.** The skin temperature ships as `skin_temp_*`, the black-body temperature of the
    longwave up. MEDS has no separate ground skin to plumb: snow-free, its ground surface is the top
    soil layer.
  - **Removed.** The four variances (`cas_temp_var_site`, `soil_temp_top_var_site`,
    `cas_vpd_var_site`, `leaf_temp_var_site`) and their operator are gone. Nothing used them, and
    they squared one end-of-step sample a day, which is not the within-step spread their names
    suggested.
  - **To add one back**, first choose its definition, then accumulate a sum of squares inside the
    step.
    - **The variance of the site-mean series** is what a tower sees. It needs each sub-step's site
      mean, which the FAST staging forms only when that tier is on.
    - **The variance pooled over patches and time** comes directly from a sum-of-squares row in the
      patch block.
- **Fast-only rows in a slow-only run read 0.** *Candidate.* [#299](https://github.com/xiangtaoxu/MEDS/issues/299) v0.3.0 weights the patch block
  with the slow step when the fast loop is off, so the slow rows report their rates. The rows only
  the fast loop fills (fluxes, the forcing echo) still read 0 there, where `_FillValue` would say
  "not simulated"; that needs a per-variable flag in the registry.
- **One serializer for site and region files.** *Candidate.* [#312](https://github.com/xiangtaoxu/MEDS/issues/312) The region writer duplicates
  the site writer, and the forcing echo is enumerated in three places.

---

## 7. Plant physiology

- **PER-COHORT thermal acclimation.** *Candidate.* [#256](https://github.com/xiangtaoxu/MEDS/issues/256) Follow-up to #176:
  Photosynthetic acclimation shipped in v0.2.0 as a **site-level** growth temperature (Kattge &
  Knorr 2007), which is what that fit is calibrated on — the mean air temperature of the preceding
  weeks. Letting a shaded understory cohort acclimate differently from a sunlit canopy one needs a
  relation calibrated on **leaf** rather than air temperature; applying the air-temperature fit to a
  per-cohort leaf temperature would use it well outside its range. Acclimation of **respiration**
  (as opposed to photosynthetic capacity) is also still open.
- **Whether storage maintenance should be ON by default.** *Open question.* The mechanism
  shipped in v0.2.0 (#177) with `storage_turnover_rate` defaulting to **0**, which reproduces the
  earlier behaviour. ED2's temperate-broadleaf value of 0.6243 yr⁻¹ costs an Ithaca run 35 % of GPP
  and 43 % of AGB over five years, so switching the default on is a rebaseline decision, not a
  parameter tweak — and it lands with the same question #118 and #176 raised about presets that were
  never climate- or biome-specific.
- **Vertically resolved root respiration** — *shipped in v0.2.0* (#178). Kept here only as the
  pointer that the remaining vertical-resolution gap is **root biomass**, not temperature: the
  response is now summed per layer against the root profile, but `broot` itself is still a single
  per-plant pool distributed by the static `root_beta` profile rather than a prognostic per-layer
  one.
- **Per-PFT hydraulic traits are shipped** (#179, v0.2.0); what remains is that nothing has
  *calibrated* them. The `[hydraulics]` traits are selectable per PFT and default to the shared
  `[hydraulics]` block, so a run that does not set them is unchanged — which means MEDS ships with
  hydraulically identical PFTs until somebody puts real trait values in. Wood density is the obvious
  axis to derive them from, as the Camac mortality coefficients already are, and `wood_kmax` now
  is (Xu et al. 2016); `wood_psi50` (denser wood ⇒ more negative) is next.
- **Phase B — per-layer root nodes.** *Planned.* [#180](https://github.com/xiangtaoxu/MEDS/issues/180) Each soil layer is now its own path into a cohort's
  wood in the water network, with no root store between. A root node would be the first unknown of
  each cohort's block in the arrowhead solve, which takes it as it stands. Sources:
  `MEDS_MULTILAYER_ROOTS_DESIGN.md` §5, `MEDS_HYDRAULICS_DESIGN.md` §16,
  `MEDS_SPAC_WATER_NETWORK_PLAN.md`.
- **Hydraulic redistribution.** *Planned.* [#181](https://github.com/xiangtaoxu/MEDS/issues/181) Root uptake is zero from a layer drier than the
  wood (a smoothed ramp in `root_water_uptake`), so it is non-negative by construction. Turning
  redistribution on is a change to that one flow law, plus deciding how the ledger treats water moving
  between layers through the plant. Same source.
- **A root radial conductance in series.** *Planned.* [#356](https://github.com/xiangtaoxu/MEDS/issues/356) The soil-to-root path is the
  rhizosphere's Gardner conductance alone; a radial root conductance in series is a different flow
  law for that path in `meds_hydr_lib`.
- **Γ* with surface pressure.** *Candidate.* [#342](https://github.com/xiangtaoxu/MEDS/issues/342) `gstar25` is a partial pressure at 100 kPa, so
  Γ* as a mole fraction rises as the pressure falls, where it should not depend on pressure. Negligible
  at BCI (Γ* +0.7 %, A −0.1 %); to fix before high-elevation sites.
- **Cohort-level sunlit and shaded leaves.** *Open question.* [#343](https://github.com/xiangtaoxu/MEDS/issues/343) Every leaf of a cohort
  photosynthesises at the cohort's mean light, as in ED2; FATES and CLM5 split each leaf layer into
  sunlit and shaded fractions.
- **Whether the soil-to-root conductance grows in proportion to fine-root length.** *Open question.*
  [#377](https://github.com/xiangtaoxu/MEDS/issues/377) The Gardner form (`docs/science/plant_hydraulics.md` eq. (14)) makes it linear in root
  length up to a logarithmic crowding term; real uptake may grow less than linearly.
- **Size-dependent physiological parameters.** *Open question.* [#381](https://github.com/xiangtaoxu/MEDS/issues/381) Turnover rates, the
  root-to-shoot ratio and leaf photosynthetic traits vary with plant size; MEDS holds them per PFT.

---

## 8. Forcing

Source: `docs/dev_plans/archive/MEDS_FORCING_DESIGN.md` §5.7, §8. Science page:
[`science/forcing.md`](science/forcing.md).

- **A better cloud term for the LWdown synthesis.** *Candidate.* [#257](https://github.com/xiangtaoxu/MEDS/issues/257)
  Follow-up to #182: The Brutsaert/Idso synthesis shipped in
  v0.2.0, so a source lacking longwave can now drive MEDS. Its cloud correction is one empirical
  coefficient on the SW clearness index, and it holds the last daytime index through the night —
  adequate for a fallback, but the residual is real: driving Ithaca from synthesis leaves the soil
  surface 1.37 K cooler than the file's `strd`. A cloud-fraction formulation (Crawford & Duchon
  1999) or a nocturnal index carried from a longer window would close more of it.
- **Forcing adapters for NLDAS-3, Daymet and CHIRPS** (forcing phase F6). *Candidate.*
  [#302](https://github.com/xiangtaoxu/MEDS/issues/302) Deferred 2026-09-27; ERA5-Land stays the
  only product. Each would follow the ED_ERA5land pattern: a downloader, the same flat monthly
  archive, and a `met_source` entry (`MEDS_FORCING_DESIGN.md` §16). NLDAS-3 (~1 km, hourly) could
  drive MEDS alone but needs regional archives and its own chunk size (OD3); Daymet and CHIRPS are
  daily corrections to an hourly base.
- **ED_ERA5land archive years before June 2002.** *Planned.*
  [#303](https://github.com/xiangtaoxu/MEDS/issues/303) Deferred until a study needs them. GDEX
  starts in June 2002, so they come from the CDS, with the tools that already built 2026-04 and
  2026-06.
- **The multi-polygon runtime.** *In beta, for one node.* [#183](https://github.com/xiangtaoxu/MEDS/issues/183) Region runs exist since R2 (#289):
  `[run].mode = "region"` runs every selected ED_ERA5land cell of a `[region]` box as its own polygon, in
  one process, all sharing one forcing reader. `MEDS_POLYGON_RUNTIME_PLAN.md` §14 (RP0–RP7) made a
  region fill one node, replacing R4–R6 of its §10:
  - a cell mask, each polygon's soil and Chave E from the `[maps]` files, and a seeded
    bare-ground start (`pft.nbg_init_density`);
  - region state files and restarts;
  - every variable and tier in the region's one file set (`detail_polygons` retired).

  Not planned: tiles, launchers and MPI; a larger area runs as several regions, one node each.
  Still open:
  - patch threads inside a region's polygons ([#400](https://github.com/xiangtaoxu/MEDS/issues/400));
  - a leaky soil bottom ([#401](https://github.com/xiangtaoxu/MEDS/issues/401));
  - the rest of the runtime consolidation ([#310](https://github.com/xiangtaoxu/MEDS/issues/310));
  - example 06 as region runs, which awaits the owner's go for its cost.
- **A source interface for the met reader.** *Done in beta.* [#311](https://github.com/xiangtaoxu/MEDS/issues/311) The two file sources
  are modules of their own (`meds_met_file_source`, `meds_met_archive_source`), each an open and one
  record's values as stored, over what they share (`meds_met_source_common`). The reader keeps
  opening, the cursor, the stepping and the one ingest. F7-F12 went in with the efficiency sweep.
- **The shortwave-to-PAR partition.** *Open question.* [#369](https://github.com/xiangtaoxu/MEDS/issues/369) The PAR share, the beam/diffuse
  split and the Weiss–Norman form are uncertain by about ±10 %: against towers' PAR sensors the
  `clearidx` fractions give PPFD 3–11 % high, mostly from the diffuse 0.52, and WN85's visible band is
  the 320–700 nm share of the solar spectrum, near-UV included (`docs/science/forcing.md` §6).
- **The neutral move to each patch's canopy-air top.** *Open question.* [#350](https://github.com/xiangtaoxu/MEDS/issues/350) The forcing is moved
  from its own heights to each patch's canopy-air top assuming a neutral surface layer
  (`docs/science/forcing.md` §8, Step 2): small at BCI, one-sided where the sensor is far above the
  canopy-air tops.

---

## 9. Snow

Source: `docs/dev_plans/archive/MEDS_SNOW_DESIGN.md` §7. Science page:
[`science/snow_biophysics.md`](science/snow_biophysics.md).

- **P1 — multi-layer snow.** *Planned.* [#186](https://github.com/xiangtaoxu/MEDS/issues/186) Compaction and densification, an aging albedo, and a
  density-dependent thermal conductivity. The single-layer store is the degenerate case.
- **P2 — canopy snow interception and unloading.** *Planned.* [#187](https://github.com/xiangtaoxu/MEDS/issues/187) Snow currently reaches the ground
  through the canopy unimpeded.

---

## 10. Source structure and testing

Source: `docs/dev_plans/MEDS_CODE_STRUCTURE_DESIGN.md` §15.

- **Pass `column_params_t` through `column_config_t`** rather than copying it into the frozen
  record every step. *Done* (2026-09-30, `MEDS_EFFICIENCY_SWEEP_PLAN.md` Phase 2): the routines take
  `col_config`, and `column_params_t` is gone. [#188](https://github.com/xiangtaoxu/MEDS/issues/188)
- **Per-layer face budget imbalance on the committed path**, per-cohort tissue residuals, and
  RK45 ledgers asserted after the rail decision. *Planned; the RK45 item is done* (2026-09-30:
  `rk45_ledgers_stop` stops only for a step it keeps). [#189](https://github.com/xiangtaoxu/MEDS/issues/189)
- **Delete `column_cohort_t`** in favour of per-field policy tables. *Closed* (2026-09-30).
  [#190](https://github.com/xiangtaoxu/MEDS/issues/190) Four of the five benefits the design claimed
  landed piecemeal: `column_cohort_init` gives the test fixtures allometric consistency, the three
  hard-coded constants are PFT parameters, the derived geometry is on the cohort block, and
  `reconcile_tissue_water_capacity` took the seed and clamp out of the gather. The fifth, a fusion
  policy that cannot omit a field, a table would not deliver either: the compiler cannot check that
  a table lists every field. The policy stays in `fuse_cohort_fast_state` and
  `scale_cohort_ground_fields`, with a witness per kind in `test_fusion_cohort`.
- **One field list for the fast-loop state** ([#146](https://github.com/xiangtaoxu/MEDS/issues/146)).
  *Done* (2026-09-30, option B of `MEDS_EFFICIENCY_SWEEP_PLAN.md` Appendix B): `state_to_array`,
  `array_to_state` and `tend_to_array` list the fields once, `state_entry_rules` says how each
  enters the error, and the combinators work on the flat array. The compiler still cannot check
  the list; `test_state_combinators` runs every field through every combinator.
- **Dead code.** *Candidate.* [#417](https://github.com/xiangtaoxu/MEDS/issues/417) `size_tol`/`area_tol`, `surface_tend_t%coh_transp`,
  `apply_rt_forcing`'s `tcas` dummy, `veg_energy_balance`'s `t_emit` generality, and unused imports.
- **Two test residues.** *Candidate.* [#421](https://github.com/xiangtaoxu/MEDS/issues/421) A duplicate `test_rk45_shed_water` in
  `test_column_ark.f90`, and `test_mass_correction_neutrality`, which guards a removed correction.

---

## 11. Performance

Source: `docs/dev_plans/archive/MEDS_GPU_EVALUATION.md` §12.

- **Attack the allocator traffic.** *Done in part* (2026-09-30). [#195](https://github.com/xiangtaoxu/MEDS/issues/195)
  With #188, #146 and the per-thread ARK storage, a patch-step allocates 104 times instead of 304,
  and the allocator's share of a serial BCI run fell from 14.9 % to 7.7 %. What remains is mostly the
  frozen record, built afresh each step (4 % of the serial fast loop): it relies on a new record's
  defaults for its scalars, so reusing it needs a reset that lists them.
- **The water network's cost.** *Within its gate.* On the BCI profile the soil-plant water network
  runs at +13.5 % against the solvers it replaced at 12 threads (+19 % serial), inside the plan's
  +15 % (`MEDS_SPAC_WATER_NETWORK_PLAN.md` §12.4). Most of it is Newton iterations (§4 above) and the
  per-cohort × per-layer root-uptake evaluation.
- **Thread and vectorise the cohort axis on the CPU.** *Planned.* [#196](https://github.com/xiangtaoxu/MEDS/issues/196) This is the evaluation's
  headline recommendation. Patch-axis threading already ships; the cohort axis is untouched.
- **A single-precision experiment** (`wp = real32`). *Decided: no.* [#197](https://github.com/xiangtaoxu/MEDS/issues/197) Measured 56× on the device for
  FP32, which is why it was worth asking what MEDS needs; the answer (2026-09-13) is that MEDS stays
  `real64` throughout. Kept here as the record; the issue closes with v0.3.0.
- **Region runs: read the next month's forcing while the current month steps.** *Planned.*
  [#407](https://github.com/xiangtaoxu/MEDS/issues/407) Part of [#183](https://github.com/xiangtaoxu/MEDS/issues/183). A region loads each month's forcing on one thread while the others
  wait; reading and compute never overlap.

---

## 12. Smaller residues

Gaps and residues found in the 2026-10-17 annotation review.

- **The snow-burial and open-canopy blends are inactive.** *Candidate.* [#410](https://github.com/xiangtaoxu/MEDS/issues/410) The fast-loop driver
  passes `snowfac = 0` and `opencan_frac = 0` to the aerodynamics, so the snow-burial roughness blend
  and the open-canopy ground-conductance blend (`docs/science/canopy_aerodynamics.md` §1, §3) never act.
- **Ground optics are bare soil plus snow.** *Candidate.* [#411](https://github.com/xiangtaoxu/MEDS/issues/411) Litter, standing water and soil moisture
  do not change the ground's albedo (`docs/science/canopy_radiation_transfer.md` §4; see also #268).
- **Maintenance respiration's missing limbs.** *Candidate.* [#412](https://github.com/xiangtaoxu/MEDS/issues/412) No maintenance respiration on
  reproductive tissue, no substrate or moisture limitation of maintenance respiration, and the stem
  size scaler is σ = 0 (`docs/science/plant_respiration.md` §7).
- **An unattributed energy residual under a melting snow pack.** *Open question.* [#413](https://github.com/xiangtaoxu/MEDS/issues/413) The whole-column
  ledger carries a sign-oscillating residual of up to about 5×10⁻³ W m⁻² (0.7 J m⁻² per 150 s step,
  about −14 J m⁻² over a day) while the canopy-air, soil, pack and pond ledgers each close to 10⁻⁹; it
  is why the energy floor in `meds_budget_check` is 10⁻² W m⁻².
- **The canopy reflects too much near-infrared at realistic leaf optics.** *Open question.* [#414](https://github.com/xiangtaoxu/MEDS/issues/414) The fast
  calibration's albedo target stays off (`scripts/calibrate_fast/calibration_reference.toml`) because
  meeting it would move a leaf's optics past its evidence; the canopy structure has to change instead.
- **`docs/ed2_comparison.md` describes v0.3.0.** *Planned.* [#422](https://github.com/xiangtaoxu/MEDS/issues/422) Deferred by the owner. Until the
  rewrite, a banner at the top of the page lists the statements that are now wrong.
- **Example 05's Palo Verde PFTs** label values "example 04's calibration" that differ from it. *Candidate.*
  [#424](https://github.com/xiangtaoxu/MEDS/issues/424)
- **`leaf_gas_exchange.md` §2 gives the θ_J borrowing cost two ways** (−25 % and 30–32 %). *Candidate.*
  [#425](https://github.com/xiangtaoxu/MEDS/issues/425)
- **Numerics checks not yet scored on the water network.** *Planned.* [#415](https://github.com/xiangtaoxu/MEDS/issues/415) The stability threshold on
  a forced multi-season run and across wind speeds; a consistent midpoint freeze to replace the
  `forcing_sample_frac` trade; a drought season against `dt_fast` (`docs/science/numerical_scheme.md`
  §7; #159 mentions the last).
