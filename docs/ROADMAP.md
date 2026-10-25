# MEDS roadmap

**What this file is.** The work MEDS has deliberately deferred, in one place. Each item names its
GitHub issue, so a deferred decision can be found, assigned and closed. An issue that closes leaves
this list.

**The rule.** A source comment states the present rationale. What changed and when goes in
[`CHANGELOG.md`](../CHANGELOG.md). What is deferred goes here, with an issue number: when you add an
item, open the issue in the same change.

**Status words.** *Planned*: intended, with a design. *Candidate*: worth doing, no design yet.
*Open question*: the decision itself is open. *Closed, to revisit*: the issue is closed because the
gain is small today, and the item stays here for the round that reopens it. Nothing here is a
commitment to a date.

Design records for the larger items are in [`dev_plans/`](dev_plans/README.md).

---

## 1. Open questions from measurement

- **Dynamic vapour pressure for leaf transpiration** (Kelvin $`e_i`$). *Candidate.* [#96](https://github.com/xiangtaoxu/MEDS/issues/96)
  Built, measured and removed. The likely route back is foliar water uptake.
- **Passive deep thermal layers below the hydrologically active soil column.** *Candidate.* [#254](https://github.com/xiangtaoxu/MEDS/issues/254)
  The `dirichlet` option of `[soil].thermal_bottom_bc` anchors the base layer to a deep temperature,
  but a purely resistive termination cannot reflect less than 0.41. Closing the rest needs heat
  capacity below the water grid: a thermal grid that extends past it.
- **Intercepted-water evaporation under ERA5-Land forcing.** *Open question.* [#433](https://github.com/xiangtaoxu/MEDS/issues/433)
  At BCI the canopy films evaporate about 34 % of the rain and 75 % of ET on ERA5-Land, against about
  11 % and 20 % on the tower's forcing. The reanalysis's drizzle, the film capacities and the form of
  the wet fraction are the candidates.
- **A film-flow term for the van Genuchten–Mualem conductivity in dry soil.** *Candidate.* [#390](https://github.com/xiangtaoxu/MEDS/issues/390)
  Until then the pore connectivity L = −1 (Schaap & Leij 2000) keeps K above ED2's floor of
  1.16×10⁻¹³ m s⁻¹ down to about −2 MPa, for a soil with n near 1.4.
- **Patch canopy height against a fixed sensor height.** *Open question.* [#481](https://github.com/xiangtaoxu/MEDS/issues/481)
  Canopy-air fluxes use each patch's own height, while a tower's observations are at one fixed height.
  See also the forcing's move to the canopy-air top (§8).

---

## 2. Phenology

Science page: [`science/plant_phenology.md`](science/plant_phenology.md). Plan:
[`dev_plans/MEDS_PHENOLOGY_SENESCENCE_PLAN.md`](dev_plans/MEDS_PHENOLOGY_SENESCENCE_PLAN.md).

- **Leaf shedding as an up-regulation of the PFT's baseline turnover.** *Planned.* [#379](https://github.com/xiangtaoxu/MEDS/issues/379)
  The turnover is a rate (`leaf_turnover_rate_toc`). The new form is $`(1 + a \cdot \text{shed tendency})`$
  times the baseline, in place of the absolute `shed_rate_max`, with a light cue relative to the
  cohort's own light. It changes examples 02 and 06, so it waits for their next refit and rerun.
- **Drought-driven fine-root dieback.** *Candidate.* [#546](https://github.com/xiangtaoxu/MEDS/issues/546) Fine roots age with the leaves'
  turnover but are not shed with them, because root length sets the root–soil conductance. A dieback
  would act per soil layer, on that layer's water.

---

## 3. Soil biogeochemistry

Science page: [`science/soil_carbon.md`](science/soil_carbon.md).

- **The nitrogen cycle.** *Planned, after the move to ED-MORE.* [#154](https://github.com/xiangtaoxu/MEDS/issues/154)
  MEDS is carbon-only. Nitrogen pools, the decomposition N limitation, mineralization and
  immobilization, and N limitation of NPP come with the move.
- **Vertically resolved soil-carbon pools.** *Candidate.* [#155](https://github.com/xiangtaoxu/MEDS/issues/155)
  The single layer is the degenerate case of the intended profile.
- **An explicit coarse woody debris pool.** *Candidate.* [#156](https://github.com/xiangtaoxu/MEDS/issues/156)
  Treefall necromass enters the structural pools directly.
- **Fire consumption of the ground pools.** *Candidate.* [#157](https://github.com/xiangtaoxu/MEDS/issues/157)
  There is no fire in MEDS.
- **Soil carbon does not decompose when the fast loop is off.** *Candidate.* [#545](https://github.com/xiangtaoxu/MEDS/issues/545)
  The soil step decays the pools by the environmental scalar the fast loop integrates, so with
  `[fast].fast_biophysics_on = false` nothing decomposes, silently. Refuse the combination, give the
  step a fixed scalar, or warn.

---

## 4. Fast-loop numerics

Science page: [`science/numerical_scheme.md`](science/numerical_scheme.md). Plan:
[`dev_plans/MEDS_PRODUCTION_INTEGRATOR_PLAN.md`](dev_plans/MEDS_PRODUCTION_INTEGRATOR_PLAN.md).

- **The `dt_fast` checks on the water network.** *Planned.* [#415](https://github.com/xiangtaoxu/MEDS/issues/415)
  At 900 s Palo Verde's February–April GPP and ET run 5.5–6 % high, linearly in `dt_fast`, because
  the gas exchange and the transpiration demand use the step-start leaf potential
  (`science/numerical_scheme.md` §5a, §7). A predictor–corrector that re-evaluates the gas exchange at
  the mid-step leaf potential would remove it, for about one extra pre-pass per step. Until then 300 s
  is the advice for dry-season leaf-habit studies. A consistent midpoint freeze to replace the
  `forcing_sample_frac` trade (§7) is the other open check.
- **The step controller's error norm is not invariant to splitting a cohort.** *Candidate.* [#444](https://github.com/xiangtaoxu/MEDS/issues/444)
  One cohort and two at half density sub-step differently (1.6×10⁻³ K on a thin stand) though the
  physics agrees to 2×10⁻¹¹ K.
- **An analytic 2×2 Jacobian for the ESDIRK2 canopy-air Newton solve.** *Candidate.* [#494](https://github.com/xiangtaoxu/MEDS/issues/494)
  The forward-difference Jacobian was 8.3 % of the fast loop when measured; re-measure first.
- **The water sub-step's error estimate.** *Candidate.* [#491](https://github.com/xiangtaoxu/MEDS/issues/491)
  It pairs the last factorisation with the latest coupling rows. The effect is at round-off level; the
  consistent form costs about 5 % of one factorisation.
- **Water network leftovers.** *Candidate.* [#547](https://github.com/xiangtaoxu/MEDS/issues/547) Newton iterations per sub-step (a warm start
  or a nonlinear elimination of each cohort), a C API entry for a plant-only or soil-only solve, and
  an unattributed 0.8 W m⁻² sensible-heat difference on BCI's dry season against the earlier scheme.
- **Adaptive freeze cadence** ([#158](https://github.com/xiangtaoxu/MEDS/issues/158)) **and soil water in
  the ARK tableau** ([#159](https://github.com/xiangtaoxu/MEDS/issues/159)). *Closed, to revisit.* The
  gain is bounded (the coefficient pre-pass is 4–11 % of a step; in-stage soil water costs 14–26 %),
  and the water network removed the borrowed-flux defect. Re-measure the drought case when the scheme
  is revisited.
- **The free-convection slope** in the tissue energy linearization. *Closed, to revisit.* [#167](https://github.com/xiangtaoxu/MEDS/issues/167)
  Measured, the missing factor is 1.02–1.07 over the leaf-minus-canopy-air range the model produces
  (1.15 at 15 K). The correct fix is a Newton linearization in the kernel with the most delicate
  conservation invariant.

---

## 5. Vegetation energy and the ground

Science page: [`science/vegetation_energy_dynamics.md`](science/vegetation_energy_dynamics.md).

- **A freezing canopy film.** *Planned.* [#165](https://github.com/xiangtaoxu/MEDS/issues/165)
  Intercepted water carries its tissue's energy and temperature but stays liquid below 0 °C. Rime and
  glaze wait for tissue freezing (below).
- **Tissue water never freezes.** *Candidate.* [#456](https://github.com/xiangtaoxu/MEDS/issues/456)
  Leaf and wood water stay supercooled liquid, and nothing defoliates or kills a plant for cold.
  Cold-driven mortality (a PFT cold tolerance, frost leaf loss, freeze–thaw embolism) would build on
  the tissue energy state.
- **The canopy-air capacity changes between steps without booking the exchange** (ED2's "zcan
  effect"). *Candidate.* [#435](https://github.com/xiangtaoxu/MEDS/issues/435)
  The slab's depth follows the tallest cohort each slow step; the heat and water it gains or loses
  are not booked. Small.
- **The sub-canopy ground conductance.** *Open question.* [#265](https://github.com/xiangtaoxu/MEDS/issues/265)
  MEDS ports ED2's non-default `icanturb = 4`, CLM4's dense-canopy form, whose ground resistance is
  8–16× that of the Sellers et al. (1986) form ED2 ships. Raising `cs_dense` from 0.004 to 0.04 raises
  BCI's LE by 6 % and its ground evaporation by 43 %. To be decided with a later refit of example 04.
- **A litter layer.** *Open question.* [#268](https://github.com/xiangtaoxu/MEDS/issues/268)
  There is no surface organic horizon for the ground energy balance or soil evaporation. To be
  measured per patch: the case behind it was a site mean dominated by gap patches.
- **One well-mixed canopy air space.** *Open question.* [#269](https://github.com/xiangtaoxu/MEDS/issues/269)
  A real sub-canopy is cooler and steadier than the air above the crowns; a layered canopy air would
  resolve it. To be measured per patch, as the litter layer.
- **Snow does not open the canopy.** *Candidate.* [#410](https://github.com/xiangtaoxu/MEDS/issues/410)
  A pack buries the canopy's roughness and displacement, but never raises the ground conductance's open
  fraction (`science/canopy_aerodynamics.md` §1, §3). Trigger: the first snowy site.
- **Ground optics.** *Candidate.* [#411](https://github.com/xiangtaoxu/MEDS/issues/411)
  Litter, standing water and soil moisture do not change the ground's albedo
  (`science/canopy_radiation_transfer.md` §4).
- **The canopy reflects too much near-infrared at realistic leaf optics.** *Open question.* [#414](https://github.com/xiangtaoxu/MEDS/issues/414)
  The fast calibration's albedo target stays off (`scripts/calibrate_fast/calibration_reference.toml`)
  because meeting it would move a leaf's optics past its evidence. The canopy structure has to change
  instead.

---

## 6. Output and diagnostics

Science page: [`science/diagnostics.md`](science/diagnostics.md).

- **A spectrally resolved surface radiative record.** *Candidate.* [#255](https://github.com/xiangtaoxu/MEDS/issues/255)
  Per-band fluxes are on the VIS/NIR/LW grid the two-stream solves. Comparing against a multispectral
  product (MODIS bands, Sentinel-2) needs finer bands: a change to the radiative transfer's band
  structure, not to its output.
- **Monthly average diurnal cycle output.** *Candidate.* [#467](https://github.com/xiangtaoxu/MEDS/issues/467)

---

## 7. Plant physiology

- **Per-cohort thermal acclimation, and acclimation of respiration.** *Candidate.* [#256](https://github.com/xiangtaoxu/MEDS/issues/256)
  Photosynthetic acclimation uses a site-level growth temperature (Kattge & Knorr 2007), calibrated on
  the preceding weeks' air temperature. A per-cohort version needs a relation calibrated on leaf
  temperature. Respiration does not acclimate.
- **Root uptake in wet soil has no root conductance in series.** *Planned.* [#356](https://github.com/xiangtaoxu/MEDS/issues/356)
  The soil-to-root path is the rhizosphere's Gardner conductance alone. A radial root conductance in
  series is a different flow law in `meds_plant_hydraulics`.
- **Hydraulic redistribution.** *Planned.* [#181](https://github.com/xiangtaoxu/MEDS/issues/181)
  Root uptake is zero from a layer drier than the wood (a smoothed ramp in
  `meds_plant_hydraulics.f90`). Allowing redistribution changes that one flow law, and the ledger must
  decide how it books water moving between layers through the plant.
- **Cohort-level sunlit and shaded leaves.** *Open question.* [#343](https://github.com/xiangtaoxu/MEDS/issues/343)
  Every leaf of a cohort photosynthesises at the cohort's mean light, as in ED2. FATES and CLM5 split
  each leaf layer into sunlit and shaded fractions.
- **Size-dependent physiological parameters.** *Open question.* [#381](https://github.com/xiangtaoxu/MEDS/issues/381)
  Turnover rates, the root-to-shoot ratio and leaf photosynthetic traits vary with plant size; MEDS
  holds them per PFT.
- **Maintenance respiration's missing limbs.** *Candidate.* [#412](https://github.com/xiangtaoxu/MEDS/issues/412)
  No maintenance respiration on reproductive tissue, and no substrate or moisture limitation
  (`science/plant_respiration.md` §7).
- **Calibrated hydraulic traits.** *Candidate.* [#548](https://github.com/xiangtaoxu/MEDS/issues/548) The shipped PFTs differ hydraulically
  only in `wood_kmax` (Xu et al. 2016's fit to wood density); `wood_psi50` is the same for all. A
  wood-density fit for it belongs in `scripts/pft_traits`.

---

## 8. Forcing and sites

Science page: [`science/forcing.md`](science/forcing.md).

- **A better cloud term for the LWdown synthesis.** *Candidate.* [#257](https://github.com/xiangtaoxu/MEDS/issues/257)
  The Brutsaert/Idso synthesis has one empirical coefficient on the shortwave clearness index and holds
  the last daytime index through the night. Driving Ithaca from it leaves the soil surface 1.37 K
  cooler than the file's `strd`. A cloud-fraction form (Crawford & Duchon 1999) or a nocturnal index
  from a longer window would close more of it.
- **Forcing adapters for NLDAS-3, Daymet and CHIRPS.** *Candidate.* [#302](https://github.com/xiangtaoxu/MEDS/issues/302)
  ERA5-Land is the only product. Each adapter would follow the ED_ERA5land pattern: a downloader, the
  same flat monthly archive, and a `met_source` entry.
- **The shortwave-to-PAR partition.** *Open question.* [#369](https://github.com/xiangtaoxu/MEDS/issues/369)
  The PAR share, the beam/diffuse split and the Weiss–Norman form are uncertain by about ±10 %.
  Against towers' PAR sensors, the `clearidx` fractions give PPFD 3–11 % high, mostly from the diffuse
  0.52 (`science/forcing.md` §6).
- **The neutral move to each patch's canopy-air top.** *Open question.* [#350](https://github.com/xiangtaoxu/MEDS/issues/350)
  The forcing moves from its own heights to each patch's canopy-air top assuming a neutral surface
  layer (`science/forcing.md` §8). The move is small at BCI and one-sided where the sensor is far above
  the canopy-air tops.
- **The multi-polygon runtime.** *In beta, for one node.* [#183](https://github.com/xiangtaoxu/MEDS/issues/183)
  `[run].mode = "region"` runs every selected ED_ERA5land cell of a `[region]` box as its own polygon
  in one process, over one forcing reader. Tiles, launchers and MPI are not planned: a larger area runs
  as several regions, one node each. Still open:
  - patch threads inside a region's polygons ([#400](https://github.com/xiangtaoxu/MEDS/issues/400));
  - one run container, one loop and one output-file-set sequence ([#310](https://github.com/xiangtaoxu/MEDS/issues/310)).
- **Several sites per polygon: what is not built.** *Planned.* [#455](https://github.com/xiangtaoxu/MEDS/issues/455)
  A polygon's sites (`[[sites]]` in a site run) each take the polygon's forcing to their elevation and
  their beam onto their slope (`science/forcing.md` §8). Issue #455 asks for an example. Not built:
  - the diffuse light a slope receives, and the light terrain reflects onto it;
  - a site map for a region (each cell's sites, shares, soils and terrain); a region runs one site per
    cell;
  - lateral water from upslope to downslope sites, and seed dispersal between sites;
  - an optional `site` axis for the stand variables (`agb`, `lai`).
- **The default forcing reader takes a numeric `_FillValue` as data.** *Candidate.* [#544](https://github.com/xiangtaoxu/MEDS/issues/544)
  `read_series` does not turn a file's `_FillValue` or `missing_value` into NaN, so a gap marked
  1e20 runs as data. The archive and map readers already refuse or mask fills.

---

## 9. Snow

Science page: [`science/snow_biophysics.md`](science/snow_biophysics.md).

- **Multi-layer snow.** *Planned.* [#186](https://github.com/xiangtaoxu/MEDS/issues/186)
  Compaction and densification, an aging albedo and a density-dependent thermal conductivity. The
  single-layer store is the degenerate case.
- **Canopy snow interception and unloading.** *Planned.* [#187](https://github.com/xiangtaoxu/MEDS/issues/187)
  Snow reaches the ground through the canopy unimpeded.

---

## 10. Code, tests and documentation

- **Dead code.** *Candidate.* [#417](https://github.com/xiangtaoxu/MEDS/issues/417)
  Unused tolerances, a never-assigned tendency field, two unused dummy arguments and unused module
  imports.
- **The `post_proc` readers use example 03's stand-layout names.** *Candidate.* [#423](https://github.com/xiangtaoxu/MEDS/issues/423)
  Four scripts read the layout example 03 writes, not `meds_main`'s names. Either teach them
  `meds_main`'s names or keep saying which script reads which layout in `post_proc/README.md`.
- **A redistributable Python wheel.** *Candidate.* [#549](https://github.com/xiangtaoxu/MEDS/issues/549)
  The wheel `python/pyproject.toml` builds is machine-local. A PyPI wheel needs a gfortran build,
  `cibuildwheel` and `auditwheel`, and netCDF's dependency closure vendored.

---

## 11. Performance

- **Thread and vectorise the cohort axis on the CPU.** *Planned.* [#196](https://github.com/xiangtaoxu/MEDS/issues/196)
  Patch-axis threading ships; the cohort axis is untouched.
- **Region runs: read the next month's forcing while the current month steps.** *Planned.* [#407](https://github.com/xiangtaoxu/MEDS/issues/407)
  A region reads each month's forcing on one thread while the others wait.
- **Decided against.** Single precision (`wp = real32`): MEDS stays `real64` throughout
  ([#197](https://github.com/xiangtaoxu/MEDS/issues/197), closed). GPU offload: OpenMP `target` offload
  measured 1.4× slower than the CPU, so the build option is removed.
