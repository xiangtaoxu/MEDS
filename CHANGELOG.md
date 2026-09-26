# Changelog

All notable changes to MEDS. Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
MEDS does not yet promise semantic versioning, because v0.x is explicitly pre-benchmark.

**This file is the only place change history lives.** Source comments state present-tense
rationale; what changed and when goes here; what is deferred goes in
[`docs/ROADMAP.md`](docs/ROADMAP.md). See `CLAUDE.md` for the rule.

Entries cite the pull request that shipped them. Anything that moves a number states the
before and after.

---

## [Unreleased]

### Added

- **ERA5-Land forcing tools in `scripts/prepare_era5/`** (#279). Downloading and post-processing are
  separate tools:
  - `download_era5land_gdex.py` fetches NSF NCAR GDEX d633008's global 5-day files unchanged into a
    raw pool, in parallel (capped at GDEX's per-user limit of 10 streams), verified and resumable;
  - `download_era5land_cds.py` fetches the Copernicus CDS, in GRIB by default: one variable ×
    12 months per request, which is half the requests NetCDF needs;
  - `postprocess_era5land.py` turns either source's raw files into NetCDF box files.

  They run in their own `meds-era5` environment (`scripts/prepare_era5/environment.yml`). On a New
  York State box, the CDS and GDEX outputs agree to within 0.00024 K. The old
  `scripts/prep_era5land_forcing.py` stays until the forcing reader upgrade lands.
- **Forcing-data design and a polygon runtime plan** (#279).
  - `MEDS_FORCING_DESIGN.md` gains Part II (§11–§19): a global per-variable monthly `ED_ERA5land_`
    archive, and a reader upgrade that reads a site or a box from it one month at a time and
    converts dewpoint to specific humidity inside the model. Its header now marks what is already
    implemented.
  - `MEDS_POLYGON_RUNTIME_PLAN.md` designs regional runs as an OpenMP loop over polygons, without
    MPI.
- **The global monthly `ED_ERA5land` archive builder** (#280).
  - `build_era5land_static.py` writes the static file: a valid-data mask defined from the data,
    ERA5-Land orography and land fraction.
  - `build_era5land_archive.py` writes one global file per variable per month,
    `ED_ERA5land_<Var>_<YYYYMM>.nc`, flat in the archive folder:
    - `Tair`, `Tdew`, `PSurf`, `u10`, `v10` as delivered; `Rainf`, `SWdown`, `LWdown`
      de-accumulated;
    - month-long chunks, and quantized;
    - it checks every hour against the static mask and plausibility bounds, records each file in a
      manifest with its checksum, and deletes raw files once verified.
  - It reads GDEX raw files or global CDS GRIB (`--source cds`). The CDS path indexes the GRIB
    headers and decodes only the valid cells, one field at a time.
  - The July 2022 pilot (GDEX) built in 457 s on 8 cores (16.0 GB). Its New York values match the
    independent box output to 0.0039 K, the quantization.
  - June 2022 was built from CDS in 411 s on 8 cores; its `Tair` and `Rainf` are bit-identical to a
    GDEX build of the same month.
- **`download_era5land_cds.py --bbox global`** requests the native global grid, and `--parallel`
  (default 3) keeps several requests in the CDS queue at once (#280).

### Changed

- **`build_era5land_archive.py --work-dir`** writes each output on another disk, such as a compute
  node's local drive, then copies the finished file into the archive in one sequential pass that
  also computes its checksum. Writing HDF5 chunks directly over a network filesystem made builds
  3–4 times slower once several ran at once.

### Removed

- **`scripts/download_era5land.py`** (#280), replaced by `scripts/prepare_era5/download_era5land_cds.py`
  and `postprocess_era5land.py`. `scripts/prep_era5land_forcing.py` now reads their box files
  (`--in` takes several files), and the forcing README, science doc, config comment and the
  `example_biophysics` instructions show the new commands.

## [0.2.2] — 2026-09-25

An **open-source and layout** release. MEDS is now licensed under the Apache License 2.0, and the
source tree is laid out so that each folder's purpose is clearer from its name. **Nothing in the
model changes.** The biophysics example's output is byte-identical to v0.2.1, and `libmeds.so`
exports the same C interface.

It is also the first release made through the **`beta` integration branch**. Pull requests now
target `beta`, which collects merged work until a release merges it into `main`; see
[`CONTRIBUTING.md`](CONTRIBUTING.md).

### Added

- **MEDS is open source under the Apache License 2.0** (#277). Before this, the repository had no
  license at all, so nobody could legally reuse the code, while `python/pyproject.toml` declared
  MIT. MEDS now ships:
  - `LICENSE`;
  - a `NOTICE` holding the copyright line, *The MEDS Authors*, and the attribution to ED2, whose
    CC BY 4.0 terms cover the portions MEDS adapts from it;
  - an `AUTHORS` file listing the copyright holders;
  - a `CONTRIBUTING.md`.

  The license also covers v0.1.0 through v0.2.1. Every source file now starts with an
  `SPDX-License-Identifier: Apache-2.0` line. The wheel declares `License-Expression: Apache-2.0`
  and carries `LICENSE` and `NOTICE`, which raises the build floor to scikit-build-core 0.11, the
  first release that reads a license expression.

### Changed

- **The foundation library is one folder again: `src/shared/`** (#276). `base/`, `functions/` and
  `util/` together are `libmeds_shared`, and they were the only library whose sources spanned
  three top-level folders. They now sit under `src/shared/{base,functions,util}`, while `config/`
  and `state/` stay top-level. `src/README.md` states the admission test: a module belongs in
  `shared/` only if it uses nothing outside it. That test keeps the regrouped folder from becoming
  the catch-all that the #125 restructure dissolved. No module, library target or Python name
  changes.
- **`src/capi/` is now `src/c_api/`, and its four shims are `meds_c_api_{leaf,phenology,demography,run}`**
  (#276). Their ctest targets are now `test_c_api_*` (select them with `ctest -R c_api`). The
  `bind(c)` names are unchanged, so `libmeds.so` exports the same C symbols and the Python package
  is unaffected. The ctypes-mirror banner in `python/meds/plant/pheno.py` now points at the
  phenology shim; it named a file an earlier restructure had removed.

## [0.2.1] — 2026-09-15

A diagnostics-and-boundaries release. Nothing here changes the demographic core; what it changes is
**what the model reports and what its lower boundaries do**, and two of the four move numbers a
reader would otherwise have trusted.

The thread that produced it is worth stating, because it is a methodological caution as much as a
set of fixes. A question about why the biophysics example's soil ran warmer than the air turned out
to have a **site-mean artefact** as its main answer: that stand carries a disturbance gap, patch LAI
spanned 0.64 to 5.36, and the gap supplied ~70 % of the apparent anomaly from 21 % of the area. Per
patch, a closed canopy was already behaving correctly — 2 cm soil sitting at the daily-mean air
temperature and running *below* air at midday. **In a demography model a site mean of a nonlinear
surface quantity is not a coarse answer, it is a wrong one**, and the example now plots one patch
rather than the stand. Chasing that question nevertheless turned up three real defects, all fixed
here, and one output gap that hid it (#270, open).

### Changed

- **The shipped examples now use the absorbing thermal bottom boundary** (#267). #145 built
  `[energy].bottom_bc = "dirichlet"`, `deep_temp`, the derived `deep_depth = 3.12 m` optimum and
  `test_soil_annual_damping`, and then nothing switched to it: the default is still `geothermal`,
  which holds the bottom heat flux at zero and makes the base an adiabatic wall that **reflects**
  the annual temperature wave. Measured in `examples/example_biophysics`, switching it moves the
  July soil temperature at the 1.73 m base node from **20.5 °C to 16.2 °C** and its annual range
  from **3.1 K to 1.2 K**; the analytic damping for that depth is 0.42 of the surface amplitude,
  against 0.80 for the adiabatic base and 0.46 for the Dirichlet one. The surface moves only
  0.29 K, so this matters for deep soil temperature — soil-carbon Q10, root-zone temperature, the
  timing of spring thaw — not for the canopy energy balance. `deep_temp` is a SITE constant with no
  defensible global value, so the default is unchanged and `meds_config_main.toml` now documents
  all three keys instead of omitting them.
- **`examples/example_biophysics` regenerated** against that boundary and against the #266 ground-
  evaporation fix. The stand ends at 14 cohorts / 2 patches, peak LAI 4.16, AGB 9.6 kgC m⁻², mean
  dbh 24.7 cm, 16.2 kgC m⁻² soil carbon (was 20 / 3, 4.06, 9.0, 23.5, 15.4). July gross uptake
  442.6 gC m⁻² against 408.5, net 206.6 against 177.7.
- **The example README is trimmed to what a reader needs.** The `dt_fast` 150 s-vs-900 s
  convergence table and the Python-driver-vs-executable agreement study were development evidence,
  not reader-facing; the substance of the latter is already recorded here under #139. The
  reader-facing rules both carried — `dt_fast` is not the output cadence, and long runs compare
  through site aggregates rather than cohort by cohort — are kept. Also corrected: the opening
  sentence said the energy balance is solved every 30 minutes; `dt_fast` is 900 s, so it is 15.

### Fixed

- **Ground evaporation used the wrong air-filled porosity, which inverted its moisture response**
  (#266). `ground_evaporation` referenced CLM5's air-filled pore space to the **bulk** top-layer
  moisture, `phi_air = phi - theta1`. CLM5 eq 5.80 references it to `theta_air`, the **air-dry**
  water content of eq 5.78 — a texture constant obtained by inverting the retention curve at
  `psi = -1e4 m`, not a state. With `phi_air` tied to `theta1` the tortuosity falls as
  `(phi - theta1)^(10/3)`, faster than the dry-surface-layer thickness falls linearly, so `r_soil`
  **rises** with wetness and **a wetter soil evaporates less**: for the shipped loam the ground
  latent flux had a *minimum* at `theta = 0.31`, and the approach to `theta_init` was a 60x cliff
  rather than a limit. At the biophysics example's actual state — `theta` 0.246, `psi` −0.0106 MPa,
  soil-surface relative humidity 0.99992, i.e. water not limiting — `r_soil` was **9167 s/m**
  against an aerodynamic 1389 s/m, pinning ground LE at 3.0 W m⁻² (3.6 % of site ET). Corrected,
  the response is monotone, joins the DSL-free limit continuously, and the closed-canopy floor's
  daily-mean warm bias drops from +1.42 K to +0.99 K. The DSL thickness denominator is corrected to
  `(theta_init - theta_air)` (eq 5.77) at the same time; the tortuosity keeps the Millington–Quirk
  exponent rather than CLM5's `phi_air^2 (phi_air/phi)^(3/B1)`, because `B1` is Clapp–Hornberger and
  the shipped curve is van Genuchten — with `phi_air` constant the two differ only by a constant
  factor, so the response *shape* is unaffected. **Note that simply disabling the DSL also removes
  the bias but desiccates the top layer to a daily-mean −17.8 MPa**: the DSL does real work and was
  referenced to the wrong porosity. `test_column_hydrology` now asserts the physics rather than the
  formula — evaporation monotone in soil moisture, and continuous into the DSL-free limit — and
  both assertions fail on the previous code. Nothing in the suite had constrained this: the only
  soil-evaporation assertions were snow-cover area weighting, which is why a 9x error in ground
  latent flux sat behind a green suite.
- **`examples/example_biophysics` plots the closed-canopy patch, not a site mean across a gap.**
  The stand's patch LAI spans 0.64 to 5.36 and the gap — 21 % of the area, passing ~70 % of incident
  shortwave straight to the ground — behaves like bare soil and runs more than 10 K above air at
  midday, supplying ~70 % of the site-mean soil warm anomaly on its own. The figure now selects the
  highest-LAI patch, and the soil series is relabelled **"Surface soil layer"**: it is the top soil
  layer (node ~1.8 cm, 0–4 cm thick), not a skin or litter temperature, and MEDS has no surface
  organic horizon to confuse it with. Per-patch sub-daily data comes from the opt-in
  `[fast].fast_probe`; the carbon and soil figures remain site means (#270).
- **Sub-daily state variables declared `AGG_TMEAN` were once-per-window snapshots, not time
  means** (#264). A family of `FLD_P_*` output sources read `site` at the output tick — once per
  `dt_slow` — and were then folded by `AGG_TMEAN` as though they had been time-integrated. For a
  quantity with a diurnal cycle that makes the "mean" one instantaneous sample per window, taken at
  whatever local time the boundary falls on, **so the bias is a function of the site's longitude**.
  Measured over one July at Ithaca: `soil_temp_top_site` read **27.79 °C** against the identical
  field accumulated properly at **26.39 °C**, and `cas_temp_site` **23.71 °C** against **22.80 °C**.
  `cas_temp`, `cas_shv`, `cas_co2`, `cas_vpd`, `soil_temp_top` and `w_surface` (site and patch
  axes, nine registry entries) now read the dt-weighted `PD_*` accumulators. After the fix
  `soil_temp_top_site` agrees with the sub-step reference to **0.005 K** and `cas_temp_site` to
  0.03 K. Stocks that are correct as instants — area, age, the soil-carbon pools, SWE, snow depth —
  are unchanged. `PD_CAS_VPD` and `PD_W_SURFACE` are new accumulator rows; VPD is nonlinear in
  temperature, so it is formed per sub-step rather than derived from the averaged twins.
  `test_output_registry` now asserts over the **whole registry** that no `AGG_TMEAN` variable is
  sourced from a once-per-tick state read, so a new variable cannot reintroduce it.
- `air_vpd` and `specific_humidity_to_vpd` moved from `src/io/meds_diagnostic_kernels` to
  `meds_therm_lib`, where the thermodynamics belongs: the fast loop now needs the same formula the
  output layer uses, and a kernel may not depend on the output layer. `meds_diagnostic_kernels`
  re-exports them, so no caller changes. The io layer's private `PRSS_REF` duplicate of `p_std` is
  retired in the same change.

- **Exactly-tied cohorts are now co-dominant in the two-stream, not stacked** (#207).
  `update_overtopping_lai` already treated equal-height cohorts as co-dominant, but
  `canopy_radiation` builds one discrete RT layer per cohort and its input is merely sorted, so
  tied cohorts **stacked** and whichever sorted first was placed above. The tie is exact and
  routine, not hypothetical: `apply_recruitment` uses **one scalar** `recruit_dbh` for every PFT,
  height is re-derived from it, and the sort is stable — so at every recruitment event the winner
  was the **lower PFT index**. Nothing about the ecology chose that. Measured on two optically
  identical tied cohorts, the lower slot absorbed **12.6 % less shortwave at zero overstory,
  rising to 15.1 % under LAI 6** — worst exactly where a whole-canopy albedo or GPP check is least
  likely to notice it, because the total is roughly conserved and only the split between PFTs
  moves. `canopy_radiation` now takes `height` and merges exactly-tied cohorts into one RT layer,
  weighting the blended optics by area index — the same weight `blend_cohort_optics` already uses
  to mix leaf against wood inside a cohort — and scattering the layer's absorbed flux back by
  absorptivity-weighted area, which is identically the quantity `leaf_frac` is built from. A
  one-member layer takes its cohort's values verbatim, so a canopy with no ties is bit-identical.
  After the fix the positional penalty is **0.0 % at every overstory LAI**.
  `test_canopy_radiation` asserts both that optically identical twins absorb equally and that
  permuting the tied pair's PFT order moves neither one; all five assertions fail on the previous
  code. **#1 stays open** — breaking the tie *ecologically* is its standing question, and this
  change does not answer it.

## [0.2.0] — 2026-09-14

**Read this before comparing a v0.2.0 run against a v0.1.0 one.** The release moved real numbers,
and the largest change is that **phenology now runs at all**. See
[`docs/ed2_comparison.md` §0](docs/ed2_comparison.md) for the before/after table.

### Known limitations, stated

Three things this release does **not** fix, carried here because they would otherwise be found by
surprise:

- **Leaf water potential does not converge at the shipped `dt_fast`** (#162, open). Daytime mean
  −0.23 MPa at 12.5 s against −1.19 MPa at 900 s: the plant water-mass update is an explicit step
  with frozen sapflow and root uptake, so the per-step excursion grows with the step. The error is
  **inherited from the canopy air and amplified about 4×**, and the residual relocates to
  `psi_wood` through the frozen uptake seam. Every other state and flux converges. Any study keyed
  to leaf water potential — hydraulic stress, potential-driven mortality — should run at ≤ 150 s
  regardless of what the carbon budget looks like.
- **Phenology is selectable and self-consistent. It is not validated.** No MEDS leaf-area cycle has
  been scored against an observation, at any site, under any strategy. The thresholds are literature
  values for the biome each strategy describes, not site calibrations — selecting `CUE_LIGHT` with
  its default 200 W/m² onset at Ithaca strips the canopy every summer, because temperate summer
  insolation sits above a threshold chosen for a tropical dry season. And it did not run at all
  before this release (#245), so **no MEDS result published before v0.2.0 had a leaf-area cycle**.
- **MEDS has still never been benchmarked.** No EDTS-equivalent regression suite, no site compared
  flux-for-flux, no output scored against observations. What is verified is internal: 49 CTest
  targets on two compilers, per-step conservation ledgers that close to machine precision, a
  per-layer face-closure check, and output byte-identical at any thread count.

Four features that move numbers substantially ship **off**, because the rebaseline window closed
before they landed: Kattge–Knorr thermal acclimation, storage-pool maintenance respiration, leaf
resorption on shed, and the non-stomatal water-stress limb. Per-PFT hydraulic traits are selectable
but uncalibrated, so MEDS ships hydraulically identical PFTs.

### Added

- **An evaluation notebook and a PFT / size-class plotter** (#175). The PFT and DBH-class output
  axes shipped in v0.1 with no reference consumer, and the v0.1 IO plan's worked evaluation against
  the Ithaca test bed was never written. Both now exist in [`post_proc/`](post_proc/).

  `plot_pft_size.py` draws composition by PFT and structure by size class — six panels — and prints
  the closure identities those axes promise (`Σ_pft agb_pft == agb_site`,
  `Σ_class agb_size == agb_site`, and the same for stem density) **before** drawing anything. All
  four hold to roundoff (worst 6.5e-16). A figure built on axes that do not partition the stand is a
  picture of a bug.

  `evaluate_ithaca.ipynb` is a worked pass over one run, using each tier for what that tier can
  answer: diurnal energy and carbon from FAST, seasonal cycles from DAILY, demographic trajectories
  and the conservation identities from MONTHLY/ANNUAL. It is explicitly **not** a benchmark — MEDS
  has never been scored against flux-tower or inventory data — but it checks what is knowable
  without observations, and that turned out to be plenty.

  **Its first run found three defects, all now filed:**

  - **#245 — the phenology cue never reaches the model.** Every PFT has been running evergreen
    whatever it declares, because `load_phenology_pft` skips its whole block unless a key the
    shipped `meds_config_pft.toml` does not document (`flush_cue_mask`) is present, which bypasses
    the presence map and leaves eleven required keys silently missing. The Ithaca reference stand,
    declared cold-deciduous, holds LAI 5.28–5.66 through every January of a 50-year run. Supplying
    the real key names produces a textbook cycle (LAI 0.34 → 5.28 → 2.89), so the phenology works —
    it has never been switched on.
  - **#247 — an out-of-bounds write in the output path.** `cohort_diag_reorder` permutes the
    diagnostic rows but never updates the block's `n`, so after a cull the block's count exceeds the
    live cohort count while `extract_variable` sizes its scratch buffer from the live one. Memory
    corruption in Release. Invisible until now because the reference stand only ever grows.
  - **#246 — the soil axis is padded to `n_soil_layer_max` with `0` / `NaN`** instead of
    `_FillValue`, so reducing over it gives a 136 K soil column, and `soil_psi_site` divides by zero
    on the padding.

  It also documents two things that are easy to get wrong reading these files by hand: turning a
  rate into an amount needs the width of **its own** window (pairing it with the previous window's
  width turned one disturbance event into a 10.8 kgC/m² phantom), and **`litter_*_site` is not all
  the carbon entering the soil** — the cull and disturbance pathways bypass that accumulator, so the
  soil-carbon budget closes from file to **0.28%** with the `mort_carbon_*_site` variables of #169
  and misses by **29%** without them.

- **A per-layer face-closure check, for the one defect class both whole-column ledgers are blind to**
  (#189, item 1 of 3). The fast loop's energy and water ledgers close to machine precision against
  the column boundary — and a purely *vertical* error survives that untouched, because enthalpy put
  in the wrong **layer** still sums correctly. This repository has paid for that three times (PRs
  #77, #85, #86): a scheme advecting soil enthalpy on a mass flux borrowed from the frozen scratch
  solve while committing θ from its own stages, so heat moves for water that never did. What found
  the last instance was an implausible temperature, which is not a detector.

  `faces[soil_layer_mass]` is: per interior layer, |the mass the soil-energy equation was charged
  for − the mass the committed θ actually moved|. It is measured where each scheme commits — on ARK
  from the scratch solve whose θ it commits verbatim, on RK45 at its own b-weighted commit, before
  the clip/floor guards edit θ. Reported next to the two budgets and asserted in the suite.

  **Measured, and mutation-verified.** Clean, it is machine-zero: **2.29e-13 kg/m²** worst over
  420,480 checks (ARK, one Ithaca year) and **2.29e-13** over 207,360 (RK45). Reintroducing the #85
  defect — handing the energy forcing `w_flux_frozen` instead of the stage's own faces:

  | | clean | with the defect |
  |---|---|---|
  | `whole_energy` | 0 fails / 103,680 | **0 fails / 103,680** |
  | `whole_water` | 0 fails / 103,680 | **0 fails / 103,680** |
  | `faces[soil_layer_mass]` | 2.29e-13 | **1.15e-01 kg/m²** |

  Both ledgers are *exactly* as clean with the defect present as without it. The new check moves
  twelve orders of magnitude. And this is at **Ithaca**, where the original temperature symptom was
  invisible (the historical note records 276.49 → 276.50 K) — so the check is strictly more
  sensitive than the accident that first caught the bug.

  Also surfaces `face_mass_resid`, which `advance_soil_water_column` had been computing and
  publishing to a field nothing read.

  **Items 2 and 3 of #189 are not in this change** and the issue stays open for them: per-cohort
  tissue residuals, and moving the RK45 ledger assertion after the rail decision. On the second,
  the concern is confirmed real — `budget_check(..., halt_budgets)` runs inside
  `column_fast_step_rk45`, *before* `column_fast_step` decides whether to roll the step back and
  rescue it on ARK, so under `[energy].debug_error` a step that was about to be discarded can abort
  the run on a ledger breach it would never have committed. The `n_fail` counting is already safe
  (the rollback restores the whole budget); only the hard stop is exposed.

- **The shipped `meds_io_config.toml` example is regenerated, and a test now keeps it that way**
  (#241). That file is generated by `meds_main --dump-io-config` but is also *tracked*, because it
  is the shipped example of the per-variable override surface that
  [`docs/science/diagnostics.md`](docs/science/diagnostics.md) points users at. Nobody regenerated
  it for five pull requests, so it had fallen **twenty variables** behind the registry — the four
  variance variables (#174), `storage_resp_site` (#177), `disturb_area_site` (#170), the five
  top-of-canopy radiative fluxes (#171), the whole five-variable fast-tier CO₂ family, and #169's
  four — plus two stale long names.

  It degrades *quietly*, which is why it drifted: an unknown name in the file is a hard error, but a
  **missing** one is not an error at all — that variable just keeps its registry-default streams. A
  user copying the example as a starting point and switching everything off would silently fail to
  switch off the variables the file forgot.

  New `io_config_example` test asserts every registered variable appears in the tracked file, and
  names the ones that do not. It checks name *presence*, not a byte diff: regeneration rewrites
  every comment and mask line, so a byte comparison would fail on cosmetic churn and end up
  suppressed.

- **Mortality carbon split by pathway** (#169). The output carried mortality *rates* and the total
  litter flux, but nothing separated the three ways a MEDS plant can die, so a stand thinning
  continuously and a stand being knocked over looked alike. Three variables now carry it:
  `mort_carbon_background_site` (the continuous hazard), `mort_carbon_cull_site` (cohorts dropping
  below the tracked size floor), `mort_carbon_disturb_site` (canopy killed when treefall opens a
  gap), plus `mort_carbon_background_patch`. Each is the whole individual — leaf, fine root, wood
  and storage — because that is what a death removes.

  **Emitted whatever `[soil_carbon].soil_carbon_on` says.** Two of the three litter sites sit behind
  that switch, and how much biomass died is a demographic question that does not stop being asked
  when the soil pools are off. The write is placed on the demography side of each guard.

  Measured on a spun-up Ithaca stand (AGB 17.5 kgC/m², `veg_carbon_site` 26.2 kgC/m²), annual means:
  **background 0.463, disturbance 0.348, cull 0.000 kgC/m²/yr** — so treefall is **42.8%** of the
  stand's mortality carbon.

  That headline decomposes rather than standing on its own. The disturbance term is the 1.39%/yr
  treefall hazard applied to 25.07 kgC/m², which is **95.8% of the stand's plant carbon** — i.e. to
  essentially the whole canopy, since only small recruits sit below the 10 m
  `disturbance_survive_height`. And the background term is an effective **1.771%/yr** of the live
  pool against `agb_mort_site`'s **1.758%/yr**, two numbers computed on entirely different code
  paths (the `PD_` patch block and the `CS_` cohort block) agreeing to 0.75%.

  **The cull reading a hard zero is physical, not a silent zero** — the failure mode #170 was. At
  the shipped `negligible_nplant = 1e-8` no cohort ever reaches the floor, because cohort fusion
  consolidates small cohorts long before they decay that far. Confirmed by raising the floor to
  5e-3 in a scratch run, where the variable reports and the stand duly collapses (AGB 16.7 → 1.3).

  The disturbance term is valued on the donor patch's *surviving* ground, so it carries a
  `1/(1-frac)` factor: the carbon died on the `frac` of the donor that became the gap, and the site
  aggregation weights by area at read time, after the shrink. Without it the site total is one-signed
  low by 1.4%/yr. Writing it on the gap instead is not an option — that slot is cleared, and carries
  no weight, for the remainder of the step.

- **Top-of-canopy radiative fluxes on the diagnostic path** (#171). The two-stream forms the upward
  flux every step and the output discarded it, so a run could not be compared against a radiometer or
  a satellite product without re-deriving it offline. Five variables now carry it: `sw_in_vis_site`,
  `sw_in_nir_site`, `sw_up_vis_site`, `sw_up_nir_site`, `lw_up_site` (surface emission included).

  **Fluxes, not a time-averaged albedo — deliberately.** A period-mean albedo is the mean of a
  *ratio*, which is not the ratio of the means, and at night the shortwave ratio is 0/0. The albedo a
  reader wants, and what a satellite product *is*, is `Σ up / Σ down` over whatever period they
  choose. Emitting a mean albedo would have produced a number that looks right and is wrong.

  Measured on a spun-up Ithaca stand, annual sums: **VIS 0.083, NIR 0.320, broadband 0.208**. The
  VIS/NIR contrast is the vegetation red-edge that every vegetation index is built on. Winter months
  rise to 0.23 / 0.45 as the snow albedo ramps in, and `lw_up_site` runs 305–468 W/m², consistent
  with σT⁴ across the annual surface-temperature range — so the numbers are checkable against
  something external, which is the point of the issue.

  No second RT solve: `canopy_radiation` already forms `albedo(b)`, so the upwelling is that ratio
  times the incident it was formed from.

- **Variance output** (#174). `AGG_MEANSQ` had existed since the IO design with no consumer:
  `normalize_scalar` computed a variance into an optional argument that nothing passed, and
  `normalize_slab` had no variance path at all. It is now `AGG_VARIANCE`, emitting the dt-weighted
  $\langle x^2\rangle - \langle x\rangle^2$ directly, with the dead `out2`/`has2` plumbing removed.

  **Registered as ordinary variables sharing their partner's source id**, not as a companion slot
  bolted to the mean. Two consequences, both deliberate: the existing per-variable buffer / normalize
  / serialize path carries them with **no new machinery** — no change to the pending record, the
  serializer, or `close_tier` — and each variance is **independently switchable** through the
  `[variables]` override, exactly like every other output. A bolted-on slot would have been neither.

  Four ship, monthly and annual only, off unless requested: `cas_temp_var_site`,
  `leaf_temp_var_site`, `soil_temp_top_var_site`, `cas_vpd_var_site`. Units are the partner's
  **squared**; take the square root for a standard deviation. Emitting the standard deviation
  directly would have lost the additivity that lets a variance combine across periods.

  Measured on a spun-up Ithaca stand: canopy-air standard deviation **2.3 K in July against 5.4 K in
  December**, and VPD 341 Pa in July against 70 Pa in January. A monthly mean cannot tell you that,
  and it is most of what a sub-daily process actually sees.

  The variance is floored at zero — the two moments accumulate independently, so round-off can put
  the difference a hair below zero for a near-constant series, and a negative variance in an output
  file is worse than a zero. `test_output_integrate` covers the dt-weighting (18.75, not the
  equal-weight 25) and the constant-series case.

### Changed

- **`[io]` is now `[state]`** (#173), with a deprecation path. The block was named for a legacy
  diagnostic writer retired at v0.1; what remained was the restart stream, so `io` named the one
  output path it did *not* cover — every diagnostic goes through `[output]`.

  | old | new |
  |---|---|
  | `[io].output_dir` | `[state].output_dir` |
  | `[io].output_prefix` | `[state].output_prefix` |
  | `[io].write_state` | `[state].write_state` |
  | `[io].state_interval_years` | `[state].interval_years` |

  The last one loses its `state_` prefix, which was stuttering once the block itself is called
  `state`.

  **`[io]` still loads** in v0.2.x, printing **one** deprecation warning for the block rather than
  one per key, and will be removed in a later release. A config that silently stopped being read
  would fall back to whatever the missing-key report produced rather than to the values the user
  wrote, which is why this is a shim and not a straight edit. Verified byte-identical: a run under
  each spelling produced 4 identical netCDF files.

  `test_capi_run` now derives a legacy-spelling config **from the shipped one** and asserts the two
  parse to the same values — so the deprecation path cannot rot into a block that reads nothing, and
  the test cannot drift from what the shipped config actually contains.

  Internally `cfg%io_*` became `cfg%state_*`. One incidental fix: `write_derived_config` in
  `test_capi_run` emitted `[io] write_state = false` and then appended the shipped config, so once
  the shipped file moved to `[state]` the new key would have won and silently flipped it back on.

### Added

- **Per-PFT hydraulic traits** (#179). All thirteen — the leaf and wood pressure–volume curves, the
  xylem vulnerability, and the conductance parameterization — can now be given per-PFT arrays in the
  `[pft]` table under the **same key names** the `[hydraulics]` block uses. Until now every PFT shared
  one parameter set, so wood density was the only axis on which PFTs could differ hydraulically, in a
  model whose point is that plant strategies differ.

  **All optional**: an absent key falls back to the `[hydraulics]` scalar, so a config can make *one*
  trait per-PFT without restating the other twelve. A run that sets none is **byte-identical** — 13
  netCDF files compared against a binary built from the parent commit.

  **Each table entry carries its own Kirchhoff lookup**, rebuilt from that PFT's `wood_kexp`. Sharing
  one across PFTs would have given every PFT the first one's vulnerability *shape* while every budget
  still closed — the defect class this repository keeps finding.

  **`solve_plant_water` itself is untouched.** The batch selects the cohort's table entry and hands
  the solver the same `hydro_params_t` it always took, so the numerically delicate part of the
  hydraulics did not change at all; only the parameter selection moved.

  **The first cut segfaulted five tests**, because the fixtures build `col_config` by hand and did not
  populate the new table — the fifth instance of the fixture-does-not-mirror-the-driver trap this
  release. Rather than patch five fixtures, `apply_hydraulics_config` now *is* the table builder and
  there is deliberately **no** PFT-uniform companion field: a fixture that forgets it fails to
  **compile**. That converts a runtime segfault into a build error, which is the only durable answer
  to a trap that has recurred five times.

- **Longwave synthesis for a forcing source without `LWdown`** (#182). `lwdown_source =
  "synthesize"` was a declared selector that parsed and then took longwave from the file anyway;
  PR #149 made `validate_config` reject it rather than let it lie. It is now implemented and the
  rejection is gone:

  ```
  LW↓ = ε_clear · σ · Tₐ⁴ · [1 + a·(1−kt)]
  ```

  `lw_clear_form` selects Brutsaert (1975), a power law in screen-level vapour pressure, or Idso &
  Jackson (1969), temperature alone — the fallback when a source's humidity is not trustworthy. `kt`
  is the **same** clearness index the shortwave partition uses, exposed rather than recomputed, so
  there is one definition of "how clear is the sky" in the forcing path. `lw_cloud_a` defaults to
  0.22.

  **The cloud term is not cosmetic.** Against the Ithaca ERA5-Land file, clear-sky Brutsaert alone
  underestimates `strd` by a mean of **−29.9 W/m²** (RMSE 43.6) on a file mean of 312.3 — clear-sky
  formulations miss the cloud enhancement and a temperate site is cloudy most of the time.

  **At night `kt` is undefined**, and `clearness_index` returns a **negative sentinel** rather than
  zero: zero would read as fully overcast and give every night the maximum cloud correction. The
  driver holds the last *daytime* index through the dark. That scalar is updated in `met_advance`,
  not `met_instant` — `met_instant` is `intent(in)` and runs per sub-step, while the advance sits
  outside the patch loop, so the state cannot become a data race when the patch axis is threaded.

  **It is a fallback, not a substitute**: driving Ithaca from synthesis instead of the file's `strd`
  leaves the soil surface **1.37 K cooler** in the annual mean. Use the file's longwave when the file
  has it.

  `LWdown` is now read *optionally* when synthesizing — demanding a variable the feature exists to
  do without would have defeated it. The two selector codes live in `meds_forcing_config`, beside
  `LW_FILE`/`LW_SYNTHESIZE`, not with the kernel that evaluates them: `meds_forcing` links
  `meds_config`, so putting them in the kernel would have pointed the config layer at the forcing
  layer and closed a dependency cycle.

- **Leaf resorption on shed** (#151), per-PFT `retained_carbon_fraction`. A fraction of the **active**
  (senescence) leaf shed returns to the non-structural pool instead of entering litter; until now
  every gram of shed leaf carbon became litter, which overstates litter input and understates the
  plant's retained reserve.

  **The leaf pool loses the full shed either way** — the design note warns that crediting storage
  while removing only the litter share *creates* carbon, so `npp%leaf` keeps the complete removal and
  the split happens between the two destinations:

  ```
  leaf   -= S_base + S_act          storage += f · S_act          litter = S_base + (1−f)·S_act
  ```

  **Only the active excess is resorbed, not the baseline.** The shed rate is
  `max(k_shed·drive, k_turn)` — a *max*, not a sum — so the active excess is exactly
  `shed_rate − base_rate`, which decomposes the max correctly in both regimes. Baseline turnover is
  excluded because `leaf_turnover_rate` is calibrated against observed **litterfall**, which already
  has resorption in it; resorbing it again would double-count. `pheno_drives_to_rates` now reports the
  baseline share so the carbon layer can make that split.

  **Default 0**, reproducing the earlier behaviour. Most measured resorption is of N and P rather than
  C, so carbon fractions are modest — 0.1–0.2 is defensible. On a **temperate-deciduous** stand at
  `f = 0.35`, Ithaca, 5 years:

  | | f = 0 | f = 0.35 | |
  |---|---|---|---|
  | storage pool | 0.000424 | 0.000688 | **+62.2 %** |
  | soil carbon | 0.001516 | 0.001449 | **−4.4 %** |
  | NPP to storage | 0.000375 | 0.000286 | −23.9 % |
  | GPP / LAI | | | +22.0 % / +12.2 % |

  The storage/soil-carbon pair is the direct signature — carbon moved from the litter path to the
  plant reserve — and less NPP then has to be diverted to refill storage. On the **evergreen** shipped
  config the effect is *exactly zero*, because that stand never actively sheds; that is the intended
  behaviour and `test_plant_phenology` asserts it.

  The slow-loop carbon ledger closes at **−4.9E-17** against a declared 7.41E-02 with `f = 0.35`.

- **Storage-pool maintenance respiration** (#177), per-PFT `storage_turnover_rate` [yr⁻¹]. The
  non-structural pool was the one live carbon store that cost nothing to hold — a cohort could carry
  an arbitrarily large reserve for free. It is now charged a fractional turnover each step, ED2's
  `growth_balive.f90` form:

  ```
  M_s = C_s · min(1, storage_turnover_rate · dt)
  ```

  **No temperature dependence**, because ED2 has none here: it sets `maintenance_temp_dep = 1.0` for
  storage and leaves the temperature form commented out as "experimental and arbitrary". Adopting one
  would go beyond the reference rather than follow it.

  **Default 0**, which reproduces the earlier behaviour — verified not by a byte compare (a new
  output variable changes the file layout) but by comparing **1677 variable instances** across 13
  monthly files against a binary built from the parent commit: worst absolute difference **exactly
  0.000e+00**, with `storage_resp_site` the only addition. ED2's own values are temperate broadleaf
  0.6243, temperate grass and conifer 0, tropical non-grass 1/6, tropical grass 1/3 — so zero is a
  legitimate ED2 setting, not an absence of physics. Where it is turned on the charge is large: at
  ED2's temperate-broadleaf rate an Ithaca run loses **35 % of GPP and 43 % of AGB** over five years,
  because the drain compounds through stand development. Whether it should be on by default is a
  v0.3.0 rebaseline question, recorded on `docs/ROADMAP.md`.

  **The ledger caught the first implementation.** Decrementing `nonstructural_carbon` in place — the
  obvious reading of ED2, which does exactly that — left the slow-loop carbon ledger with
  **−2.4339E-03 kgC** undeclared in the `allocate` phase and **+2.4339E-03** over-declared in
  `grow+mortality`, equal and opposite. The phase that *declares* the CO₂ efflux has to be the phase
  where the pool drops. Routing the charge through `npp%nonstructural` as a **tendency** — the
  standing "driver computes, engine applies" rule — closes it: the residual is now **−6.0E-17**
  against a declared 5.13E-02. The allocator still sees the post-maintenance reserve, so ED2's
  maintenance-before-growth ordering is preserved even though the pool moves one phase later.

  Reported as `storage_resp_site` [kgC/m²/yr], and it rides `co2_owed` — the same per-patch channel
  growth respiration already uses to reach the fast loop's `nee_biotic`.

### Changed

- **Fine-root maintenance respiration is summed over soil layers, not taken at a mean temperature**
  (#178). The model resolves a soil temperature per layer; the root response collapsed it to a
  root-weighted mean *before* the temperature function, throwing that resolution away. It now
  evaluates `Σ_k root_frac_k · f(T_k)` instead of `f(Σ_k root_frac_k · T_k)`.

  This is the Jensen argument #145 makes for Rh, with one extra turn worth stating: the **peaked**
  Arrhenius form is convex below its optimum and **concave** near it, so the error changes sign with
  season instead of biasing one way. Measured at Ithaca on the default 2 m column with `root_beta =
  2`:

  | | Dec | Mar | Jun | Sep | annual mean |
  |---|---|---|---|---|---|
  | layered vs mean-temperature | **+4.8 %** | +1.8 % | **−8.3 %** | −0.3 % | **−0.17 %** |

  So it is a **seasonal** correction to root respiration, not an annual-budget one — reporting only
  the annual figure would understate it by a factor of thirty.

  The weighted scale is patch-uniform, so it is evaluated once per patch and broadcast over the
  cohort array; `fine_root_maintenance_respiration` accordingly takes a temperature *scale* rather
  than a temperature, which keeps it `elemental` over cohorts. `test_plant_respiration` asserts that
  a **uniform** profile reproduces the old answer exactly (so the change is inert without a
  gradient) and that a **split** profile with the identical weighted mean does not.

### Added

- **Thermal acclimation of photosynthetic capacity** (#176), opt-in via
  `[leaf_physiology].thermal_acclimation`. Kattge & Knorr (2007): the peaked form's entropy term and
  the capacity ratio follow a running-mean growth temperature $T_g$ (°C),

  ```
  dS_v = 668.39 − 1.07·Tg     dS_j = 659.70 − 0.75·Tg     Jmax25/Vcmax25 = 2.59 − 0.035·Tg
  ```

  so a warm-grown and a cold-grown stand of the same PFT no longer share one temperature response.
  Measured across a 20 K range of growth temperature, the Vcmax optimum moves **28.2 °C → 35.0 °C**.

  **Measured in the coupled model, and decomposed** — because the headline number is misleading on
  its own. Ithaca, 5 years, growth temperature ≈ 10 °C:

  | | off | dS only | dS + ratio | dS | ratio | total |
  |---|---|---|---|---|---|---|
  | GPP | 0.004286 | 0.004823 | 0.006118 | +12.5 % | +26.8 % | **+42.8 %** |
  | NPP | 0.003675 | 0.004154 | 0.005332 | +13.0 % | +28.4 % | +45.1 % |
  | LAI | 0.048615 | 0.053505 | 0.062691 | +10.1 % | +17.2 % | +29.0 % |

  **Two thirds of the effect is the capacity ratio, not the optimum shift.** At Ithaca's ~10 °C
  growth temperature Kattge & Knorr put Jmax25/Vcmax25 at **2.24** against the PFT file's fixed
  **1.7** — a 32 % higher Jmax. The fixed value is a single global number being applied at a cold
  site; acclimation replaces it with a climate-dependent one. That is the same class of finding as
  #118's Vcmax review: a preset that was never climate-specific.

  Three deliberate choices, each recorded in `docs/science/leaf_gas_exchange.md`:

  - **The ratio acclimates too.** Shifting dS alone would acclimate the *shape* of each response
    while pinning the two branches in fixed proportion, which is not what the study measured.
  - **$T_g$ is the growth temperature, not the leaf temperature** — the fit is calibrated against the
    mean **air** temperature of the preceding weeks, so MEDS tracks an exponential running mean of
    the daily mean (window `acclim_window_days`, default 30 d).
  - **Site-level, not per cohort**, for the same reason. Per-cohort acclimation would need a
    relation calibrated on leaf rather than air temperature; it is on the ROADMAP, together with
    acclimation of *respiration*, which this does not touch.

  Applied by refreshing the per-PFT leaf table once per slow step from the **config** reference
  coefficients, which are never overwritten — so no kernel signature changes, the law lives in one
  place, and the refresh is idempotent. The running mean is prognostic and slow, so it is written to
  the state file; rebuilding a month of memory on every restart would make the optimum jump.
  `thermal_acclimation` **requires** `temp_response_form = "peaked"` — the Arrhenius form has no dS
  to shift, and the loader refuses rather than letting the flag be silently inert. Off by default,
  so the shipped path is unchanged.

- **All five phenology cues are wired; the drought-deciduous and light-exchanging strategies now run
  from configuration alone** (#150). `validate_config` used to reject the `WATER(2)`, `HYDRO(4)` and
  `LIGHT(16)` cue bits in either mask, because the kernel computed all five cues while the driver fed
  three of them zeros. Both halves are closed.

  **Four cue accumulators became per-cohort state.** `water_avg`, `low_psi_days`, `high_psi_days` and
  `light_avg` lived only inside `advance_leaf_phenology`, where `state = pheno_state_t()` re-zeroed
  them every slow step. A 10-day running mean reset daily is just its own instantaneous input, and a
  *consecutive*-dry-day counter reset daily never exceeds one — so those cues could not have worked
  even with their drivers present. They now carry the lockstep reorder, every creation site, and a
  **declared** fusion policy (survivor-keeps, stated explicitly rather than left implicit).

  **The drivers**, each a fast-loop daily reduction unless noted:

  | cue | driver |
  |---|---|
  | `CUE_WATER` | root-weighted fraction of extractable water, `Σ f_root,k · clamp01[(θ−θ_wp)/(θ_fc−θ_wp)]` |
  | `CUE_HYDRO` | the cohort's published predawn leaf potential (`dmax_psi_leaf`, already reduced for #95) vs a turgor-loss point **derived from the same pressure–volume curve** the leaf stress arrestor uses |
  | `CUE_LIGHT` | daily-mean incident shortwave at the canopy top (ED2 `rad_avg`) |
  | cold-drop | **top-layer soil temperature**, replacing the air-temperature proxy |

  The soil-temperature swap is a **number-mover for the already-shipped temperate-deciduous
  strategy**, not only an unlock: soil lags and damps air, so an air proxy crosses the 284.3 K and
  275.15 K cold-drop thresholds earlier in autumn than the soil does.

  **Nine `[phenology]` cue parameters became configurable**, through a new optional per-PFT array
  loader (`opt_pa`). Five had no table entry at all and four had a table entry but no loader, so the
  WATER/HYDRO/LIGHT cues could not have been tuned even once their drivers landed. They are optional
  — absent keys take the table defaults — because those cues are opt-in and a temperature-strategy
  config should not have to supply five numbers it never reads. A wrong-length array is still an
  error.

  **Acceptance, measured.** Both strategies run from a TOML config on the Ithaca driver, and
  reproduce the design's patterns 3 and 4 end to end: drought-deciduous holds a full canopy when
  watered, sheds under sustained drought (shed drive 0.996) and **reflushes on rewet** (0.996) — the
  reflush being the part that needs the persisted counters; light-exchanging sheds with light (0.9999
  bright vs 0.004 dim) while its flush stays permissive (1.000).

  **Selectable is not validated, and the release notes say so.** The kernel is unit-tested for all
  four strategies, the four drivers each have a hand-computed unit test, and both new strategies are
  asserted end to end — but **no MEDS run's leaf-area cycle has ever been scored against a phenology
  observation, under any strategy**, including the two that have shipped since v0.1.0. The thresholds
  are literature values for the biome each strategy describes, not site calibrations: selecting
  `CUE_LIGHT` with its default 200 W/m² onset at Ithaca strips the canopy every summer, because
  temperate summer insolation sits above a threshold chosen for a tropical dry season.

### Fixed

- **Soil-dimensioned output wrote its inactive padding as data** (#246). Every soil variable was
  emitted over all `n_soil_layer_max` (20) slots of what is normally a 10-layer column, because
  `layer_source_field` set `nlayer = n_soil_layer_max` against a comment that already said "the
  ACTIVE layer count". The `DIM_SOIL` branch marks every slot it is handed as valid, so the padding
  was written as data rather than left for netCDF to fill.

  Two consequences. Reducing over the soil axis gave nonsense — `soil_temp_site` averaged over its
  own soil dimension read **136 K**, ten real layers averaged with ten zeros. And
  `soil_matric_potential` was evaluated on padding whose `theta_sat` is 0: a divide by zero that a
  Debug build (`-fpe0`) **aborts** on and a Release build turns into NaN in the file. That abort was
  only reachable with `[output].water_fluxes = true`, which the shipped Ithaca config has off — so
  the one configuration that crashes was not one the suite or the test bed exercised.

  `nlayer` is now `dp%n_soil`, and `soil_z` is written over the same extent. The padding is
  therefore never written and comes back as `_FillValue`, masked by any CF-aware reader — which is
  exactly how the cohort and patch axes have always handled their own unused capacity. Measured
  after: `soil_temp_site` over the unmasked layers reads **273.8 K**, no unmasked non-finite value
  anywhere, and the Debug build completes with `soil_psi_site` enabled.

  `test_output_roundtrip` now runs a 6-layer column against the 20-slot ceiling and asserts the tail
  is fill for both the variable and the coordinate. Its source array is 290 K in *every* slot
  including the padding, so a writer that emitted the ceiling would read 290 and pass for the wrong
  reason; it has to read `_FillValue`.

- **Phenology was silently disabled in every run** (#245). **This is the largest behavioural change
  in v0.2.0: any PFT with a `[phenology]` block now actually has a leaf-area cycle.**

  `load_phenology_pft` gated the whole section on one key inside it —
  `toml_has(t, 'phenology.flush_cue_mask')` — and that is a key the shipped `meds_config_pft.toml`
  never documented. It documented `cue_mask`, a name no reader consumes. So a config that wrote a
  full, deliberate `[phenology]` block was **skipped in silence**: the presence map never saw
  twenty-three required keys go missing, and every PFT fell back to `CUE_NONE` / `CUE_NONE` —
  `flush = 1`, `shed = 0`, the evergreen fixed point — whatever leaf habit it declared.

  The Ithaca reference stand is declared cold-deciduous, has sane thresholds, and reaches 270.5 K
  soil temperature in winter. It held **LAI 5.28–5.66 through every January of a 50-year run**, with
  `flush = 1.0000, shed = 0.0000` in every month. **No MEDS run has ever had a leaf-area cycle.**

  Three changes:

  - **The gate is the section, not a key inside it.** New `toml_has_section` asks whether the author
    meant to configure phenology at all, rather than whether they spelled one particular key the way
    the loader wanted. A `[phenology]` block missing keys is now a hard error naming all of them.
  - **`cue_mask` is rejected by name**, with a migration message giving both replacements and the cue
    bits. It was retired when the flush and shed cues became independently selectable; ignoring it
    silently changed a PFT's leaf habit, which is the one outcome worse than stopping.
  - **The shipped `meds_config_pft.toml` block is corrected and completed** (all twenty-three keys,
    and a note that `evergreen = [0]` alone does not give a PFT a season), as is
    `examples/example_biophysics/pft_parameters.toml`, which used the stale name live along with two
    more keys — `on_threshold` / `off_threshold` — left over from the retired status/deadband design.

  With the real key names, the Ithaca stand runs as it always should have: **LAI 0.001 (Feb) → 4.72
  (Sep) → 0.70 (Dec)**, shed governor 0.97 in January and 0.00 through summer. The phenology kernel
  was never broken — it was never switched on.

  This is the **second** occurrence of the absent-key trap that `.claude/rules/config.md` already
  describes, in a different block. The rule text described it exactly.

  **Migration:** a config carrying `cue_mask` now stops with instructions. A config with no
  `[phenology]` section is unchanged (evergreen defaults), which is the documented fallback. A config
  that had a complete block under the real key names was already working and is unaffected.

- **An out-of-bounds write in the output path: the diagnostic blocks kept a stale count after a
  cull** (#247). `cohort_diag_reorder` and `patch_diag_reorder` permuted the diagnostic rows but
  never updated the block's own `n`, while `cohort_reorder` set `cohort%n = m` right after calling
  them. So after any operator that **shrinks** the array, the block claimed more live slots than the
  array had — and `extract_variable` sizes its scratch buffer from the live count while
  `cohort_diag_value` filled `1..d%n` into it. A write past the end: silent in Release, `subscript
  97 > upper bound 96` under Debug, a SIGSEGV some steps later.

  Invisible until now because the reference stand only ever **grows** — Ithaca goes 114 → 125
  cohorts over three years, so the block's count was never the larger one. It becomes routine the
  moment phenology works (#245): shedding drives cohorts below the tracking floor,
  `terminate_cohorts` culls them, and a two-year run dies in its first December.

  Fixed at both ends. The count now rides the lockstep like every other field, and
  `cohort_diag_value` / `patch_diag_value` bound their write by `size(x)` as well as by the block —
  the second is what turns the next count mismatch into a short read instead of memory corruption.
  `test_containers` asserts both after its cull, and fails on the unfixed code.

- **Slow per-patch diagnostics emitted a per-SECOND rate under a per-YEAR label** (#239). **Six
  shipped variables read a factor `yr_sec = 3.1557e7` too small**: `litter_leaf_site`,
  `litter_leaf_patch`, `litter_fineroot_site`, `litter_struct_site`, `nplant_recruit_site` and
  `disturb_area_site`. Anyone who plotted litterfall or recruitment from a MEDS run before this
  should re-read those series.

  The patch diagnostic block is `(value, weight)` and the reader returns `value/weight`. The weight
  is `Σ dt` in **seconds** — `accumulate_patch_diag` adds `dt_fast` every fast step — so every row
  must contribute *(its rate, in the units the registry declares)* × *(dt in seconds)*. The fast rows
  always did, with `flux * dt_fast`. The five slow rows, added later, contributed a bare per-step
  **amount** (and, for recruitment, `rate * dt_years`) — dimensionally a per-second rate.

  On a spun-up Ithaca stand (114 cohorts, AGB 16.7 kgC/m², LAI 5.56), annual means before → after:
  `litter_leaf_site` **9.80e-09 → 0.309 kgC/m²/yr**, `litter_fineroot_site` **1.21e-08 → 0.380**,
  `litter_struct_site` **2.07e-08 → 0.654**, `nplant_recruit_site` **1.11e-09 → 0.0349 plant/m²/yr**.
  The corrected leaf number is the independent check: LAI 5.56 at SLA ≈ 20 m²/kgC is ≈ 0.28 kgC/m² of
  leaf carbon turned over annually.

  **Why nothing caught it.** A diagnostic is downstream of every conservation ledger, and the slow
  ledger's own litter term reads `patch%litter_in` directly rather than this slot, so the two never
  had to agree. The one test that did touch a slow slot — `test_disturbance`, added with #170 —
  asserted the **raw accumulator** rather than the number `patch_diag_value` emits, so it certified
  the contribution and was blind to the convention it was written against. It now asserts through the
  reader, and a new `slow_diag_units` test asserts, for all five slow slots, the identity the
  per-year label *means*: `reported rate × elapsed years == the amount that actually flowed`, against
  oracles outside the diagnostic path (`patch%litter_in`, the recruit pool, the configured hazard).

- **The phenology memory now survives a restart** (#150). No phenology state was written to the state
  file at all — not the two governors, not the GDD and chilling sums. A restart resurrected every
  cohort at its **birth** values (`flush_drive = 1`, `shed_drive = 0`, `gdd = chill = 0`, the
  evergreen fixed point), so a temperate-deciduous stand restarted in January came back with flushing
  permitted and no chilling accumulated — it leafed out in midwinter and then had to rebuild weeks of
  thermal memory. This affected the two strategies that have shipped **since v0.1.0**, not only the
  two #150 adds. None of it is derivable from the instantaneous state: these are time integrals over
  the preceding weeks. All eight columns are now written, and optional on read, so an older state
  file still restarts exactly as it did.

### Examples

- **Every example regenerated against v0.2.0**, and two of them were broken.
  `example_leaf_gas_exchange` raised `TypeError` on every run: it passed `theta=p.theta_j` to
  `assimilation_demand_c3`, which #118 renamed when it split the C3 co-limitation curvatures out —
  literally the confusion #118 was about, sitting in the example. `example_demography` is the entry
  below. `example_phenology` now genuinely exercises all four strategies; two could not be selected
  before #150.

  **`example_biophysics` is re-spun and genuinely cold-deciduous for the first time** (#245). It
  ends at 20 cohorts / 3 patches, peak LAI 4.06, AGB 9.00 kgC/m², mean dbh 23.5 cm, 15.4 kgC/m² soil
  carbon — against 128 / 12 / 5.32 / 15.75 / 35.2 / 24.2 when it ran evergreen. Establishment is
  about twice as slow: a deciduous seedling rebuilds its canopy from storage every spring and spends
  its early growing seasons near break-even, so little happens until year 15. Its trajectory figure
  now plots the **annual maximum** LAI — it sampled at the year boundary, which is 1 January, and for
  a deciduous stand that plotted a bare canopy topping out at 0.28 for a stand whose July LAI is 4.1.

  **Its July stage drops from `dt_fast = 150 s` to the 900 s production default.** The short step was
  justified on the argument that a diel figure needs sub-daily fidelity and a long step smears the
  energy partitioning. Measured over that exact July, it does not: canopy air 22.69 vs 22.67 °C mean
  and 9.34 vs 9.37 K diel amplitude, leaf 23.32 vs 23.31 and 12.92 vs 12.91, soil surface 23.42 vs
  23.39, leaf−air +4.04 K by day either way, `dmax_psi_leaf` −0.2951 vs −0.2954 MPa. Every number the
  figures report agrees to 0.03 K, for 14.3 s of wall time against 4.0 s. `dt_fast` is an **accuracy**
  parameter and not the output cadence — the hourly records come from `fast_interval_steps` — and the
  config and README now say so, along with what 900 s *is* wrong for (the sub-daily `psi_leaf`
  excursion, #162, for which MEDS emits no diagnostic).

- **The demography example is slow-scale demography only** (#260), which is what it always claimed
  to be: cohort and patch dynamics, fusion and fission, growth, mortality, recruitment and treefall.
  **No carbon dynamics and no soil carbon.**

  It had stopped being that without anyone noticing. Its vital-rate laws are LAI-driven and
  empirical, and the reorg moved them out of Fortran into `empirical_laws.py`; the Fortran model has
  only the *carbon* path. So the example's config, brought up to schema in September and never
  re-run, was silently pointing `meds_main` at **a different model** — stub GPP
  (`gross_gpp = gpp_ref * leaf_area`) with no light competition, hence no negative feedback on leaf
  area. It diverged to LAI 501 and AGB 907 kgC/m² and died in `cohort_reorder`. The committed output
  beside it was older still, from the deleted *Fortran* empirical model.

  The example's driver is now the Python one that actually implements its laws, and it writes its own
  output: new `meds_site_get_patch_real` / `meds_site_get_patch_int` in the demography C-API expose
  the patch areas, ages and the cohort→patch CSR map, and `_write_nc.py` writes the ragged
  cohort/patch netCDF `post_proc/` already reads.

  **The stand now equilibrates and is physically sensible**: 250 years settles at ~0.97 stems m⁻²,
  **AGB 16.5 kgC m⁻², LAI 7.4**, 266–355 cohorts over 12 patches, with a textbook inverse-J size
  distribution — 0.43 stems m⁻² below 1 cm DBH down to 0.0008 above 50 cm, and **76% of the biomass
  in stems over 20 cm**. The previous committed output had AGB 121 kgC m⁻², about four times the
  densest forest on Earth. The run takes **22 seconds**, against 2 min 18 s for the carbon run that
  crashed, and `empirical_spinup.py` still reproduces its golden exactly.

  The stale `example_output_pft_parameters.csv` is deleted: it described the deleted Fortran model
  (`growth_lai_slope`, `mort_gamma/alpha/beta`, no carbon traits), and a provenance record for a
  model that does not exist is worse than none.

### Documentation

- **The leaf water-stress divergence from ED2 is now recorded as a decision, not an open question**
  (#47), in `docs/science/leaf_gas_exchange.md` §4.3. MEDS keeps the Sabot et al. (2022) two-limb
  scheme; ED2's `farq_katul.f90` uses a Manzoni-style single `stoma_beta` plus a turgor-loss
  capacity term. The three divergences and the reason each one stands: the capacity limb's shape
  (linear $\psi_{leaf}$ ramp vs ED2's 6th-power turgor-loss) and its targets ($V_{cmax}$, $J_{max}$
  and TPU vs ED2's $J_{max}$ and $\alpha$ only, which acts only on the light-limited branch); the
  stomatal limb's scope (Leuning/Medlyn $g_1$ **and** Katul $\lambda$ vs Katul only, which would
  leave the two explicit schemes with no drought response at all); and the parameter split, where
  MEDS's $s_{ref}\,e$ is exactly ED2's `stoma_beta`, so a calibrated value transfers as
  `stoma_beta = -sref_stomata * lambda_psi_exp`.

  **The `psi_soil` wiring defect the issue also filed was already closed by #95** and needed no work
  here: the field is `psi`, the driver passes `dmax_psi_leaf` on both call paths, and unset cohorts
  are seeded from the surface-layer soil potential. Verified rather than assumed.

### Fixed

- **Saturation over a frozen surface now uses the ice curve** (#89). `sat_vapor_pressure`,
  `sat_specific_humidity` and both their temperature derivatives take an optional `fliq` (liquid
  fraction). Absent means pure liquid, which is exactly what every caller did before, so omitting it
  is **bit-identical**. Supplied, they blend:

  ```
  e_sat = fliq * 611.2 exp(17.67 Tc/(Tc+243.5))  +  (1-fliq) * 611.2 exp(21.87 Tc/(Tc+265.5))
  ```

  Both branches carry the same 611.2 Pa constant, so they cross **exactly** at `Tc = 0` and the
  blend is continuous in temperature *and* in `fliq` — a pack that freezes or melts slides between
  the curves instead of stepping, which matters because a step here lands in the right-hand side an
  adaptive controller integrates.

  **What was wrong.** `snow_surface_fluxes` drove sublimation with `sat_specific_humidity` on the
  liquid curve, while the enthalpy side of the same routine was already ice-aware — removing
  `enthalpy_vapor` from an ice-referenced layer debits sublimation (vaporization + fusion)
  automatically. So the model treated the snow surface as ice for *energy* and as liquid for *vapour
  pressure*. Over ice `e_sat` is **10 % lower at −10 °C, 22 % at −20 °C and 34 % at −30 °C**, so the
  driving gradient was overstated by those factors. The same applied to ground evaporation from
  frozen soil, which now uses the top layer's `soil_fliq`.

  **Measured** (Ithaca, 5 years, `dt_fast = 900 s`), over the 18 months with snow on the ground:

  | | liquid curve | ice branch | change |
  |---|---|---|---|
  | latent heat flux | 0.96083 | 0.91988 | **−4.26 %** |
  | snow water equivalent | 4.85938 | 4.87776 | **+0.38 %** |
  | sensible heat flux | −2.69802 | −2.69890 | −0.03 % |
  | soil surface temperature | 276.595 | 276.596 | +0.000 % |

  Peak SWE rises 0.18–0.26 % in each of the five winters; the February shoulder gains 2.1 %. The
  flux responds by less than `e_sat` does because what drives it is the *gradient*
  `q_sat(T_s) − q_CAS`, and the winter canopy air is often near saturation — and because the Ithaca
  pack spends much of its life near 0 °C, where the two curves nearly coincide. At a cold
  continental site holding −20 °C for months the correction is a much larger fraction.

  **Not applied to dewpoint conversion or diagnostic VPD.** Dewpoint is *defined* over liquid, so an
  ice branch there would mis-convert the forcing. Also not applied to the canopy, which has no ice
  state at all — intercepted snowfall is held as liquid with no fusion debit, which is a separate
  known gap.

  The ED2 proposal this issue points at (EDmodel/ED2#442) concerns the liquid formula. Measured
  against Murphy & Koop (2005), MEDS's existing Bolton form is within **0.2 %** from −30 to +40 °C,
  so the liquid curve was never the problem and is left alone. The ice form added here is within
  0.1 % at −10 °C and 0.9 % at −30 °C.

### Changed

- **One solar declination for the whole model** (#152). `solar_cosz` used Cooper (1969),
  `daylength` used White (1997); both now call `solar_declination(doy)`, which is the Cooper form.
  The two were the same function offset by 1.25 days — $-\cos x = \sin(x-\pi/2)$ puts White's
  ascending zero crossing at doy 82.25 against Cooper's 81 — and Cooper's is the closer to the true
  vernal equinox near doy 79–80. It is also the form the radiation path already used every
  `dt_fast`, so `daylength` moved rather than `solar_cosz`.

  At Ithaca (42.44 °N) daylength changes by at most **3.7 minutes**, at the equinoxes, and by
  essentially nothing at the solstices; the 10.5 h autumn phenology cue fires **2 days earlier**
  (doy 297 → 295). `solar_cosz` and therefore the radiation are unchanged.

  `test_time` now asserts the shared source by **inverting each consumer back to a declination**
  rather than re-implementing the formula, so a future divergence is caught even if the formula
  itself changes. Verified to fail on the pre-fix code.

### Added

- **C3 gets its own co-limitation curvatures** (#118): `theta_cj_c3` (default 0.98, the
  $A_c$/$A_j$ transition) and `theta_ip_c3` (0.95, the transition with $A_p$). Both are new
  **required** keys in the `[pft]` table. C3 previously passed `theta_j` — the curvature of the
  electron-transport hyperbola — into both of its smoothings, while C4 already had a dedicated pair.
  A co-limitation curvature and a light-saturation curvature share units and a functional form and
  nothing else; CLM and FATES put the C3 $A_c$/$A_j$ curvature near 0.98 against the 0.7–0.9 that
  fits the $J$ hyperbola.

  **At the leaf**, PFT-1 kinetics, $C_i = 280$ µmol mol⁻¹, saturating light
  ($A_c = 14.40$, $A_j = 17.53$, $A_p = 45.0$, so $\min = 14.40$):

  | curvatures | after $A_c$/$A_j$ | after $A_p$ | total vs min() |
  |---|---|---|---|
  | `theta_j` = 0.85 (before) | 11.31 (−21.4 %) | 10.80 (−4.5 %) | **−25.0 %** |
  | 0.98 / 0.95 (now) | 13.49 (−6.3 %) | 13.22 (−2.0 %) | **−8.2 %** |

  i.e. **+22.4 % gross assimilation** at that point. The shortfall was nearly independent of
  $V_{cmax}$ (30, 31, 31, 32 % at $V_{cmax,25}$ = 60, 90, 120, 150), so it was a systematic offset,
  not a regime effect — no measured $V_{cmax}$ reproduced a measured rate.

  **In the coupled model**, Ithaca, one year from a common spun-up stand with **demography frozen**,
  so stand structure cannot diverge and the difference is physiology:

  | | `theta_j` | `theta_*_c3` | change |
  |---|---|---|---|
  | GPP | 0.07633 | 0.09039 | **+18.4 %** |
  | NPP | 0.06496 | 0.07869 | **+21.1 %** |
  | NEE | −1.4888 | −1.8274 | −22.8 % |
  | LAI (allocation still runs) | 0.9796 | 1.0128 | +3.4 % |

  With demography **live**, 10 years from cold start, the same change reads GPP +66 % and AGB +89 %
  — that is the *compounded* figure, because a persistent rate increase accelerates stand
  development, and it is not an equilibrium sensitivity. The controlled number above is the one to
  quote.

  **The $V_{cmax,25}$ presets were re-examined and deliberately left alone.** The concern was that
  they might have been tuned against the over-smoothing, in which case the compensation would have
  been hidden inside a parameter named for a different process. They were not: `vcmax25 = [60, 45,
  40]` entered in the commit that first added the leaf module and has never been revised, the values
  sit mid-range for their PFT descriptions, and the only other GPP-facing knob (`gpp_ref`) is the
  stub used when the fast loop is off. So this is a correction, not the unwinding of a calibration,
  and the presets should not be lowered to absorb it.

  `leaf_photo_params_t%theta_cj`/`theta_ic` are renamed `theta_cj_c4`/`theta_ic_c4` so the pathway is
  visible at the point of use — which is the whole defect — and the C API mirror, `meds/plant/_ffi.py`
  and `LeafParams` follow. `meds_assimilation_demand_c3` now takes two curvatures instead of one.

- **A Dirichlet bottom thermal boundary for the soil column** (#145), selected by
  `[energy].bottom_bc = "dirichlet"` with `deep_temp` [K] and `deep_depth` [m]. The bottom node
  conducts to a plane held at `deep_temp`, a distance `deep_depth - |z_node(n)|` below it. The
  default stays `geothermal`, and the Neumann path is bit-identical to before.

  **The premise, re-measured against an exact oracle.** A homogeneous column at constant `theta` has
  constant `kappa` and `C`, so the analytic semi-infinite solution `exp(-z/d)` with
  `d = sqrt(2*alpha/omega)` holds exactly — an oracle independent of the model's own machinery. The
  harness is validated by a 12 m column with the adiabatic base, which reproduces that profile to
  0.1 % over the top three damping depths. At the default 2 m column, amplitude relative to the
  surface at the bottom node (−1.73 m):

  | bottom BC | amplitude ratio | vs analytic 0.421 | RMS error over the profile |
  |---|---|---|---|
  | `geothermal` (adiabatic) | 0.764 | **+82 %** | 0.145 |
  | `dirichlet`, `deep_depth = 3.12` | 0.412 | **−2 %** | 0.015 |

  **The anchor depth is derived, not fitted.** A resistive termination reflects least when its
  impedance `kappa/l` matches the magnitude of the half-space impedance `kappa*(1+i)/d`, i.e. when
  `l = d/sqrt(2)`. For the default column that puts the anchor at `1.727 + 1.973/sqrt(2) = 3.12` m; a
  measured sweep puts the optimum at 3.0–3.1 m, so the derivation is right to within one sweep step.

  **In the coupled model** (Ithaca, 10 years from cold start, `dt_fast = 900 s`, `deep_temp = 283.24 K`
  = the forcing's mean annual air temperature; final-year monthly means):

  | | `geothermal` | `dirichlet` | change |
  |---|---|---|---|
  | annual swing at −1.73 m | 26.25 K | 15.39 K | **−41 %** |
  | annual mean at −1.73 m | 290.47 K | 286.41 K | **−4.06 K** |
  | Rh, annual mean | 0.00240 | 0.00230 | **−4.2 %** |
  | GPP / NPP / AGB / LAI | — | — | +0.6 % / +0.8 % / +1.0 % / +0.7 % |

  The adiabatic profile *flattens* below −0.45 m (swing 30.0, 28.6, 27.9, 27.1, 26.3 K) — the
  reflection signature — while the anchored one keeps decaying (29.0, 26.7, 23.7, 20.2, 15.4 K). The
  Rh change is seasonal in shape, not a level shift: the August peak falls 0.0057 → 0.0051 while
  April rises 0.0016 → 0.0018, which is the Jensen argument in #145 made visible. The soil-carbon
  difference (+3.6 %) is a 10-year transient-accumulation difference, not an equilibrium one — the
  pool is still climbing at the end of the run.

  **What it does not fix, stated plainly.** A purely resistive termination cannot reflect less than
  0.41 in amplitude at any `l`: matching a complex impedance with a real one leaves the phase wrong by
  45°. Closing the rest needs heat *capacity* below the column — passive deep thermal layers under the
  hydrologically active one — which is on `docs/ROADMAP.md` as the #145 follow-up. The anchor buys
  about a factor of ten in this metric, not exactness.

  `deep_temp` is **required** when the Dirichlet BC is selected. An error in it is a steady flux
  `kappa/l * error` into the column base, so a silent default would reintroduce a mean-annual
  deep-soil bias — the very thing this boundary condition removes.

### Fixed

- **The PFT-parameter CSV dump lost its last two columns and printed a bit pattern** (#118). The
  `write_pft_params_csv` format had 44 item slots against a 46-item output list. A Fortran format
  shorter than its list does not fail — it reverts to the last repeat group and keeps going — so an
  integer edit descriptor silently received a real (`fineroot_turnover_rate` printed as
  `4605380978949069210`) and `f_labile_stem` / `struct_lignin_frac` vanished. Introduced by adding
  the two curvature columns in this same release and caught by a new column-count assertion in
  `test_pft_optics_config`, which is verified to fail on the broken format.

- **The forcing file and the config are now checked against each other** (#185). Five items, decided
  individually:
  - **`dt_forcing` is validated against the file's actual record spacing**, and the axis against
    itself for uniformity. This was the real hazard and it is *not* the one the issue named: the
    value was read straight from the config and **never compared with the file at all**, while
    placing interval midpoints, disaggregating shortwave and bracketing the recycle seam. A config
    saying 3600 s against a half-hourly file mis-timed all three, plausibly. Checking the spacing is
    strictly stronger than checking the `timestep_seconds` attribute, which stays provenance.
  - **`avg_convention` and `sw_input_kind` are read and validated** when present, skipped when
    absent so older files still load. `sw_input_kind = "total"` against `sw_partition =
    "passthrough"` used to crash on a missing variable; it is now a clear rejection.
  - **`avg_convention = "instant"` and `"center"` are rejected.** Both parsed and then ran the
    end-of-interval path anyway, because only `METAVG_BEGIN` has a branch — so selecting either got
    a scheme the user did not ask for, silently. Same precedent as `lwdown_source`.
  - **`SWPART_SIB` deleted.** A reserved code with no implementation and no TOML spelling: it could
    never be selected, and `partition_shortwave`'s default branch would have routed it to Erbs.
  - **`elevation(grid)` stays unread**, deliberately: site elevation is a `[site]` property and the
    variable is provenance about the source grid.

  Reported through `met_open`'s existing status channel rather than `error stop`, which is what
  makes each rejection testable — the module's own comment says that is what the channel is for.
  The synthetic test fixture now writes the two global attributes the prep script writes; without
  that the new checks would have been skipped in the test while firing in production.
- **The tissue-water floor now reports the water it creates** (#148). `advance_water_mass_full`
  floors `leaf_water_mass` and `wood_water_mass` at a tiny positive value so the linear mass Euler
  step cannot go negative — and creates water doing so. The code's own comment called the case
  "unobserved in this pass's test scenarios", which was **a belief, not a measurement**: the floor
  fires per cohort per tissue, while the whole-column water ledger sums leaf + wood over all
  cohorts, so water created in one cohort's wood is indistinguishable from a redistribution between
  cohorts. The budget closed to ~4×10⁻¹² kg m⁻² with the floor entirely unmonitored. It is now
  reported through the existing `budget%clamp_mass` / `clamp_commit_n` commit-clamp channel, which
  already reduces to `site%work_clamp_mass` and already has an output variable — so it surfaces
  end-to-end with no new reporting surface, which is what the issue asked for. Measured on a forced
  fixture: 2 activations creating 2.92×10⁻² kg m⁻², and **zero on an ordinary step**. The claim in
  `ark2_column_step` that its commit counter "stays 0 by construction" was corrected — that claim
  was the whole of the defect.
- **`PD_DISTURB_AREA` is written and emitted** (#170). The slot was declared in the patch
  diagnostic block with no writer and no registry row, so the disturbed-area flux read as a **silent
  zero** rather than a missing variable — the harder failure to notice, and one no conservation
  check can see, because zero disturbed area is a perfectly conservative answer. Written on the
  donor patches inside `apply_patch_disturbance` *before* the gap is appended, which is where the
  diag slots still line up with the donors; writing after the append is the out-of-bounds trap this
  file has already paid for once. New output variable `disturb_area_site`.
- **A run selecting `time_integrator = "rk45"` above `dt_fast` = 300 s is now warned** (#160). The
  transpiration corrector that cut a ~1 MPa `psi_leaf` error by 314× (PR #91) lives in
  `advance_water_mass_full`, which ARK calls and RK45 does not, so an RK45 production run silently
  carried an error the default path does not. A warning rather than an error, because RK45 is the
  accuracy baseline and is meant for a fine step.
- **Γ\* now responds to `o2_mol_frac`, so the O₂ knob propagates to both places oxygen enters the
  C3 demand** (#117). The compensation point is set by Rubisco's CO₂/O₂ specificity and is
  proportional to the O₂ partial pressure, but it was computed with no O₂ dependence at all — so
  raising O₂ correctly inhibited carboxylation through `Kc(1 + O/Ko)` while leaving the entire
  photorespiratory penalty on `Aj` untouched. Scaled by `o2_mol_frac / 0.209`, the O₂ the shipped
  `gstar25` was measured at (Bernacchi et al. 2001). **The factor is exactly 1 at the default, so
  every shipped configuration is bit-identical** (10.18931 µmol m⁻² s⁻¹ before and after on the test
  fixture). Away from it, measured on a midday tropical leaf: halving O₂ raises A_net by 22.0 %
  against 8.2 % before, and 35 % O₂ cuts it by 23.7 % against 10.5 %.
- **Two roadmap items described code that no longer exists**, found by re-measuring every filed
  premise before scheduling it (#168, #166). `bsap` stopped being a placeholder in PR #125 —
  `set_cohort_wood_geometry` derives it from ED2's real `b1SA`/`b2SA` sapwood-area allometry. The
  old `0.10 * wood_carbon` placeholder made the wood thermal time constant **6.5–10× too short**
  across the whole size range (`f_sap` runs 1.00 at dbh ≤ 19 cm to 0.655 at 117 cm, against 0.10).
  `veg_energy_step_implicit` was deleted in PR #120, and `veg_energy_diagnostic` does not exist
  either; `veg_energy_balance` is the single closure. Both entries corrected in `docs/ROADMAP.md`
  and their four stale source/doc references removed.
- **`scripts/numerics_sweep.py` could not run a cross-scheme comparison.** The `--parity` preset
  pinned three config keys that no longer exist (`fast.integration_scheme`,
  `fast.leaf_energy_model`, `fast.wood_energy_model` — #199 named two), and the `SCHEMES` table
  still offered `split` and `picard`, which are now a hard error. The preset is removed rather than
  repointed: every difference it pinned has since been closed by making the schemes agree.
- **The C-API demography shim inlined the allometry** instead of calling `meds_allometry`, so a
  coefficient change updated the model and not the shim (#200). It now calls `dbh_to_height`,
  `dbh_to_agb`, `dbh_to_leaf_area`, `size2leaf_carbon` and `size2wood_carbon`.
  `examples/example_demography/empirical_laws.py` reads the four allometry coefficients from the
  `[allometry]` block of its shipped PFT config instead of hard-coding them; the values are
  unchanged, so the example's behaviour is unchanged.

### Documentation

- **`fuse_cohort_fast_state` says what it is: the one declaration of the per-cohort fast-state
  policy, centralised but not enforced** (#190, deferred). Four of the five benefits that justified
  deleting `column_cohort_t` have already landed piecemeal — `column_cohort_init` gives the test
  fixtures allometric consistency (the `bwood`-on-uninitialized-memory bug is fixed), the three
  hard-coded canopy constants are PFT parameters, the derived geometry is on the cohort block, and
  `reconcile_tissue_water_capacity` took the seed and clamp out of the fast gather. The fusion and
  scaling policy is centralised too, in `fuse_cohort_fast_state` and `scale_cohort_ground_fields`.
  What is left is **completeness you cannot forget** — a table the blend *iterates* cannot omit a
  field that a hand-written routine can — and that is #146's silent-omission class. The two are now
  **paired for v0.3.0**: one packed, policy-carrying layout should serve the fast state vector and
  the cohort slice together, because separately each is a large refactor buying a fraction of one
  property.
- **The RK45 stiff rescue keeps its whole-step rollback** (#161), measured rather than optimised.
  E5 proposed snapshotting at the point of failure to avoid redoing accepted sub-steps. Over a full
  simulated year at Ithaca on `rk45` at `dt_fast` = 900 s the rescue fires **zero times**
  (`work_rk45_rescue_site` = 0, against 70 577 integrator sub-steps on the same run — the counter is
  live and the zero is real), so there is no work to save on the reference workload. It is also not
  free to build: retaining part of RK45's boundary-flux accumulation while ARK finishes the interval
  means one `dt_fast`'s ledger summing two schemes' contributions, and it is well-defined only for
  the `stiff_bail` trigger — `rk45_state_railed` tests the *final* state, so on that path no
  sub-step is identified to resume from. Recorded at the site.
- **The FAST output tier's staging path is not a duplicate switchboard, and is kept** (#172). It was
  slated for deletion on the belief that it duplicated the general extraction path. It does not:
  the tier is *already* on the general machinery — it shares the registry, the buffers, `close_tier`
  and the serializer — and the only bespoke part is the extraction *source*, which is forced by when
  the values exist. Sub-daily quantities are sampled **during** the fast loop; by the time `main`
  folds them, `site` holds the post-fast-loop, end-of-slow-step snapshot, so resolving them against
  live site state would silently emit end-of-day values on a sub-daily axis. The two functions
  cannot overlap either: `SRC_S_*` occupy 4000–4999 and `SRC_F_*` 5000–5999, disjoint by
  construction. Deleting the staging would have deleted sub-daily sampling. Recorded at the site.
- **The frozen-seam contract is written down** (#201):
  [`docs/dev_plans/MEDS_FROZEN_SEAM_CONTRACT.md`](docs/dev_plans/MEDS_FROZEN_SEAM_CONTRACT.md). The
  Λ = F·dt/S criterion and its case split (linear-in-store is scale-free and sound; a prescribed
  flux against a prognostic store is unsound, with the +30.7 µmol m⁻² s⁻¹ NEE scar to prove it), the
  four seams classified against it, why Λ is meaningless for a *rate* seam and what replaces it
  (debit-before-credit), and the arbitration rule — scale all demands by `min(1, S/D)` — for a store
  with several consumers, which is order-independent where per-process clamping is not.
- **The ED2 two-stream defects found during the port are folded into
  [`docs/ed2_comparison.md`](docs/ed2_comparison.md) §5a** (#6), with the MEDS ↔ ED2 RT structure
  mapping. All six are reported upstream. Two of them — the stale-PFT diffuse index and the missing
  clumping factor in the longwave split — are *impossible by construction* in MEDS, and that is why
  its RT is shaped the way it is.
- **`test/` stays flat, deliberately** (#193, structure-plan decision #13), recorded in
  `src/README.md` with the reason: the discipline that matters is the link line, not the directory.
- **The year-rollover `rh_seam_gap` residual is attributed** (#192). 8.370×10⁻⁴ kgC m⁻² at a year
  boundary is the annual patch cadence, not a leak: the fast loop accumulates against one patch
  composition and the daily step debits another, blended through a matrix nonlinear in the lignin
  fraction. Documented at the field, with the instruction not to widen a tolerance to absorb it.
- **The build files no longer describe GPU offload as the parallel path** (#194). Measured, the
  offload build runs **1.4× slower** than the CPU (49.8 s against 36.2 s), one kernel sits at 0.4 %
  occupancy, and the device treated as 20 slow cores is 26× slower than 4 CPU cores. Patch-axis CPU
  threading is the parallel path. `MEDS_GPU=gpu` is kept as a reproducible experiment.
  `CMakeLists.txt`, `docs/building.md`, `src/README.md` and `docs/ed2_comparison.md` corrected.
- **A v0.2.0 release plan** in [`docs/dev_plans/MEDS_V02_RELEASE_PLAN.md`](docs/dev_plans/MEDS_V02_RELEASE_PLAN.md):
  all 66 open issues triaged into six phases plus release, 48 in scope and 18 deferred to v0.3+,
  with twelve decisions recorded.
- **The documentation was reorganized against the restructured source tree.** Thirty design
  plans moved to `docs/dev_plans/archive/` with tombstones; eleven stay live or as reference.
  The README dropped from 310 lines to a reader's entry point, with building, configuration
  and post-processing moved to their own pages. A new [`src/README.md`](src/README.md)
  documents the source layout, the placement rules and the library graph.
  `CLAUDE.md` split into a short always-loaded file plus path-scoped rules. This file and
  `docs/ROADMAP.md` were created. Three missing science pages were written:
  `docs/science/soil_carbon.md`, `forcing.md`, `plant_respiration.md`.

### Changed

- **The soil-energy adaptive-substep surface is deleted** (#163), after measuring rather than
  assuming. `soil_energy_step_implicit` hard-codes `flux%nsub = 1` and backward Euler is
  unconditionally stable, so an inner substepper could only ever buy *accuracy* — and the accuracy
  of the integrated state is already owned by the outer adaptive march through `GRP_SE` (soil
  internal energy, which is what the state vector carries). A second controller on a quantity the
  first one already controls is not worth wiring.
  The dead surface was wider than filed: besides `substep`, `h_init` and `max_substep`, `[energy].rtol`
  fed only `GRP_SOIL_T` — **a tolerance group with no member in `state_wrms_grouped`**, which sums
  `GRP_SE` and `GRP_THETA` and never `GRP_SOIL_T`. The group, its `ATOL_SOIL_T_DEF` default and the
  `ENERGY_SOLVER_BE` / `ENERGY_SUBSTEP_ADAPTIVE` single-valued enums went too; `N_TOL_GROUP` is 7.
  `[energy].atol` and `debug_error` stay — both are live, `atol` as the budget-closure threshold for
  the debug halt, not a step tolerance. `bottom_bc` stays because #145 wires the Dirichlet thermal
  anchor onto it in this same release.
  One coupling was preserved deliberately: `atol_scale` used to reach `[energy].atol` by
  round-tripping config → `tols(GRP_SOIL_T)` → config, so deleting the dead group would have
  silently dropped it. It is now applied directly, and a test asserts it.
  **Verified byte-identical**: all 14 output files of a 1-year Ithaca ARK run — restart state and
  every monthly diagnostic — are unchanged.
- **The two unreachable heterotrophic-respiration kernels are deleted** (#153).
  `heterotrophic_respiration_damm` (Davidson 2012) and `heterotrophic_respiration_flux` (Q10 / ED2
  capped exponential) were tested but **could not be selected**: there was no `hr_model` TOML key
  anywhere and `co2_opts_t` was never carried by `meds_config`, so the selector could not be set
  and the kernels could not be called. A tested-but-unreachable kernel is the worst of both worlds —
  maintenance and review weight for nothing, while reading to a newcomer as an available option.
  There is **one** production Rh authority: the CENTURY matrix.
  The dead surface was wider than the issue recorded: the `HR_*` selector codes, `co2_opts_t`,
  `damm_params_t`, the shared `water_modifier` helper, three DAMM-only constants in
  `meds_constants`, and unused `co2_opts_t` imports in `meds_fast_types` and `meds_fast_prepass`.
  Net **−237 lines**. The implementations are preserved on branch `archive/damm-hr` (at `2fcb647`).
  Five subtests in `test_column_co2` went with the kernels they exercised; in
  `test_soil_biogeochem` only part (a) of the fast/slow seam test went — part (b), the diurnal
  accumulation and the Jensen counter-check, builds its own factors inline and is untouched.
- **The IMEX-Euler oracle tier is retired** (#198). It was not an oracle: it returned no reference
  trajectory, and `imex_euler_column_step` was a two-line wrapper around `column_be_stage` plus
  `advance_water_mass_full` — the production scheme's own kernels — so its independence was in the
  tableau, not in the machinery it was supposed to check. All five of its consumers kept their
  coverage: four now call a test-local `be_euler_step` (the same two-line composition, living in the
  code that degrades the scheme), and `test_adaptive_march` was **ported to `adaptive_ark_march`**,
  so it now exercises the production controller instead of one that existed only to serve the tier
  (8 sub-steps at `rtol` 1e-3 against 25 at 1e-6). `meds_fast_rk4_oracle` holds one oracle and its
  name is accurate again.
- **The RK4 oracle's independence is now structural.** With the tier gone, the module no longer
  imports `meds_fast_be_stage` at all — no `column_be_stage`, no `newton_surface_solve`, no
  `advance_water_mass_full`. It sees the pure right-hand side and the state algebra and nothing
  else, which is what makes agreement between it and an implicit scheme rule out a shared-bug false
  pass. Before, the module imported the BE machinery for the tier's benefit while the oracle itself
  never touched it.
- **One implementation of each test assertion helper** (#191). Twenty-three local copies across
  nineteen files, consolidated behind generic interfaces so all ~1000 call sites compile unchanged;
  net −403 lines. The two families (fatal condition-first, accumulating name-first) are kept
  deliberately — they differ in failure behaviour, not just signature. A new `meds_test_assert`
  holds the assertions and depends on nothing but a kind, because seventeen tests deliberately link
  one narrow library each. Found and closed one coverage hole: `test_biogeochem_dynamics`' local
  helper error-stopped while the shared one accumulates, which turned four assertions into no-ops
  until a verdict call was added.

### Fixed

- **`soil_carbon_on` now actually defaults to on.** PR #144 flipped the in-type default to
  `.true.` but the TOML loader still passed `.false.` as its absent-key default, and
  `toml_logical` returns the supplied default when the key is absent. Every config that
  omitted the key, including the shipped `meds_config_main.toml`, therefore ran with soil
  carbon off. Off is not a coarser soil-carbon model, it is no soil carbon: litter is
  discarded and `rh = 0`.
- **`forcing.lwdown_source = "synthesize"` is rejected rather than silently ignored.** The
  value parsed to `LW_SYNTHESIZE` but the met driver always read LWdown from the file, so it
  selected the file path under a name promising Brutsaert/Idso synthesis. `validate_config`
  now stops on it until the synthesis exists.

### Changed

- Test coverage for the three untested fast-loop state combinators, with two silent
  exclusions named (#147).
- `soil_carbon_on` default flipped to `.true.` in the type (#144). Measured over a year at
  Ithaca, off reports annual-mean NEE at −4.657 against −2.464 µmol m⁻² s⁻¹, an 89 % stronger
  apparent sink, with Rh identically zero against 0.833 kgC m⁻² yr⁻¹, for no wall-clock saving.
- `[energy].phase_change` deleted; ice-aware soil thermal properties are unconditional (#143).
  The flag only ever gated ice-aware `κ_sat(f_liq)` and `C_eff(f_liq)`; the difference was
  ≤ 0.61 K and 35 692 substeps either way.

### Added

- **The full coupled model is drivable from Python** (`meds.model.Run`): open, step, finalize
  through a C-API shim (#139). The Python and executable paths differ only by libm-versus-libimf
  interposition: worst relative difference ~1×10⁻¹² over a simulated day.
- **One `libmeds.so`** built by scikit-build-core, with a mandatory ctest target per C-API shim
  so an ABI change is a build failure in a default build rather than a silent break in an
  optional one (#138).
- **A slow-loop conservation ledger**, with birth and death as paired transfers (#132).

### Fixed (slow-loop conservation, #132–#137)

The ledger closed a series of real carbon, energy and water leaks that a green test suite and
the per-store fast budgets had both missed:

- **Growth respiration was 23 % of growth carbon that was never exhaled.** The allocator's
  outputs now all have a destination (#133).
- Tissue thermal mass is an exchange, not an appearance (#134).
- The time-averaged diagnostics now close for energy and water everywhere (#135).
- Mortality is valued on what the applier actually removed, not on what was requested (#136).
- Reproduction carbon became a flow rather than a disappearance, which closed the ledger
  (#137). Side effect: the recruit pool now accrues every slow step instead of as a monthly
  lump, so the first cohorts appear about a month later and the early trajectory is offset.
  The offset decays as the stand fills: 70 % in AGB at year 2, 1.6 % at year 10, 0.1 % at
  year 40.
- **Soil respiration reached the atmosphere 964× too small** — a CENTURY Rh unit error, hidden
  because no shipped config turned soil carbon on. The seam check that would have caught it was
  computed but never reported (#140).
- An out-of-bounds litter read: a per-patch quantity indexed outside the patch block when
  disturbance added a patch mid-step. The ledger declared the same garbage it consumed, so
  conservation balanced on it (#139).

### Changed (source-tree restructure, #125–#127, #141)

The tree is now **timescale-first for processes** over a **state layer in two halves**. Moving a
file changes only CMake wiring, because Fortran `use` is by module name, so every step was
verified byte-identical on both back ends.

- Steps 1–6: state as a layer (`state/column`, `state/site`), processes timescale-first
  (`fast_dynamics/`, `slow_dynamics/`) (#125). `src/shared/` dissolved into `base/`,
  `functions/`, `config/`, `state/`.
- Steps 9–10: every façade module dissolved; the fast-loop argument vocabulary renamed so
  names tell the truth — `wcap`/`ccap` → `cas_mass_capacity`/`cas_molar_capacity`,
  `gah`/`gaw`/`gac` → `g_atm_heat`/`g_atm_vapour`/`g_atm_co2`,
  `uext_to_temp`/`temp_to_uext` → `internal_energy_to_temp`/`temp_to_internal_energy`, and the
  three different things called `hydro` split into `soil_water_opts` versus
  `hydraulics_params`/`hydraulics_opts` (#126). netCDF registry strings and TOML keys were
  deliberately not renamed, so output files and configs are unchanged.
- Steps 0 and 7: the fast-loop state vector (#127).
- The library name `demography` was renamed to `core` in July and back to `demography` here:
  `core` stopped carrying information once its state half became `state/site`.
- Structure plan §15 phased remainder, Phase 1 executed (#141).

### Changed (configuration, #129, #130)

- **The soil column comes from config, not from literals in the driver** (#129): the new
  `[soil_column]` block carries layer count, depth, grid growth, hydraulic texture, retention
  family, root profile and the three thermal properties, validated at load. `depth` is the knob
  for the known too-shallow-column defect (2.0 m against a ~2.5 m annual damping depth).
- **The fast driver's hard-coded parameters found their homes** (#130). `agf_bs` duplicated the
  per-PFT `aboveground_frac` with a hard-coded 0.7, so a run whose PFTs differed in allocation
  used their values everywhere except stem respiration (#128). Canopy optics — leaf and wood
  reflectance and transmittance per band, clumping, leaf-angle mean and standard deviation —
  became `[pft]` traits, so two PFTs can finally differ in how they intercept light (#131).
  Longwave is configured as emissivity, with reflectance derived as `1 − emissivity` and
  transmittance zero, because a leaf is opaque at thermal wavelengths.

### Fixed (2026-09 review, #119–#124)

- Conservation fixes: a closed whole-column energy ledger, RK45 and threading defects,
  film-water and seam conservation (#119). A 50-year spin-up's energy imbalance went from
  +6.75 MJ m⁻² to 3.6×10⁻¹¹ W m⁻².
- Dead code removed, shared column-state algebra extracted, `t_ground` made explicit, one CAS
  box kernel instead of two (#120).
- The fast-loop argument vocabulary renamed; over-long lines wrapped (#122).
- `column_prepass` split into the five processes it fused; a per-PFT leaf photosynthesis table;
  scalar signatures; shared assemblers; ledger-consistent latent and sensible heat (#123).
- The frozen work record decomposed by physical content; `integrator_opts_t`;
  `apply_process_mask` (#124).

---

## [0.1.0] — 2026-08-02

First tagged release. **Unbenchmarked**: no EDTS-equivalent regression suite has been run, no
site compared flux-for-flux, no output scored against observations. What is verified is
internal — the test suite on two compilers, per-step conservation ledgers, and thread-invariant
output.

### Added

- **The v0.1 diagnostic output layer** (#111): per-variable, per-timescale output control.
  ~208 registered variables across 8 groups and 7 axes (cohort, patch, site, soil, PFT,
  DBH size class, and the 2-D soil×patch slab), each switchable individually per timescale
  from TOML. Extensive quantities carry their own aggregation weight, so one registry line
  emits a cohort field's patch, site, PFT and size-class rollups. `meds_main --dump-io-config`
  lists every variable. Verified byte-identical at 1 versus 4 threads across all 75 files of a
  3-year run.
- **The ED2-to-MEDS comparison** for people who already run ED2 (#113, #115),
  [`docs/ed2_comparison.md`](docs/ed2_comparison.md).
- **Patch-axis threading** (#109): byte-identical output at 1, 2, 4 and 8 threads. 2.03× at
  4 threads against a measured 3.03× hardware ceiling. Three silent compiler traps were found
  and worked around: Intel's `-auto-scalar` default placing local arrays in static storage,
  nvfortran rejecting `BLOCK` inside a parallel region, and ifx building `private` copies of a
  derived type through a compiler-generated static mold.

### Changed

- **`dt_fast` became an accuracy parameter, not a stability one** (#90). The period-2
  canopy-air oscillation was traced to one frozen coefficient — the canopy-air-to-atmosphere
  conductance — and re-solving it at every integrator stage removed the oscillation while
  *reducing* integrator work. The production default went from 150 s to **900 s**. The
  non-stomatal water-stress limb was gated off by default in the same change.
- The legacy `[io]` diagnostic writer was retired at v0.1: an annual-cadence, instantaneous,
  hard-coded 21-variable schema that duplicated `[output]` and collided with it on the `-D-`
  filename prefix. `meds_io` now carries the state (restart) stream only.

### Fixed

- **The transpiration seam in the plant water update** (#91): a pure flux inconsistency,
  exactly `dt·(transp_pp − transp_bw)`, corrected by a transpiration corrector. Improved
  `psi_leaf` convergence 314×. ARK only; RK45 does not carry the corrector.
- Stomatal water-stress closure (#98, issue #95): `beta_stomata` was identically 1 until fixed.
  Restart persistence for the daily-maximum leaf water potential.
- Pathological hydraulics sub-stepping is now detected (#105): a collapsed, floored wood store
  burned 13× the wall clock, silently (issue #104).
- The adaptive warm start got per-patch storage (#108, issue #106) — it had been loop-carried,
  so patch 1 cold-started and patches 2..N inherited a neighbour's step size.

### Documentation

- **GPU offload evaluated against measurement and found not viable as scoped** (#110). The GPU
  build ran 1.4× *slower* than the CPU (49.8 s against 36.2 s), with one kernel at 0.4 %
  occupancy; the GPU treated as 20 slow cores was 26× slower than 4 CPU cores. The
  recommendation is to thread and vectorise the cohort axis on the CPU instead.

---

## Pre-0.1 development — 2026-06-23 to 2026-08-02

MEDS began on 2026-06-23. The sections below group the first six weeks by subsystem rather than
by date, because the work proceeded as a dozen parallel subsystem builds.

### Demographic core

- Source tree organized by process domain; per-PFT `hgt_max` (#3).
- The engine became **carbon-driven**: carbon pools on the cohort structure-of-arrays (#17),
  carbon-driven growth wired in (#19), the empirical vital-rate laws moved out to the Python
  example and the engine reduced to law-free apply-primitives (#16, #18, #43).
- The demography engine reorganized into `core` with a tendency-seam growth/mortality
  interface (#44, #45).
- **Run-model decision** (#63): the fast biophysics loop is always on, and phenology is
  unconditional. A slow-only, empirical-demography run is the Python C-API path.

### Fast-loop biophysics

- **Canopy radiative transfer** (#5): a faithful modernized ED2 two-stream (`icanrad = 2`) with
  a unified multi-band solver, SCOPE/4SAIL leaf-angle scattering over a Beta leaf-angle
  distribution, and an in-house block-tridiagonal solve. ED2 bugs found during the port are
  recorded in issue #6.
- **Plant hydraulics** (#9, #49): a stateless matrix-exponential network solver, shared
  constitutive curves (pressure-volume, Kirchhoff conductance, xylem vulnerability) extracted
  into `meds_hydr_lib`, and multi-layer root water uptake.
- **Leaf gas exchange** (#2, #46): FvCB C3 and Collatz C4 demand, Leuning / Medlyn / Katul
  stomatal models, a bracketed C_i solver, and Sabot two-limb water stress.
- **Column soil hydrology** (#22, #23): implicit backward-Euler Thomas Richards with Celia
  modified-Picard linearization, upstream-weighted conductivity, adaptive step-doubling,
  infiltration and ponding, and a free-drain / bedrock / aquifer bottom boundary.
- **Energy balance** (#24, #25): four stateless per-store thermal kernels carrying prognostic
  **internal energy rather than temperature**, so the freeze/thaw plateau is a read-off of the
  enthalpy inverter rather than a special case.
- **Canopy air space CO₂** (#26, #27): `can_co2` as the third prognostic CAS twin, with a DAMM
  heterotrophic-respiration option.
- **Fast-loop coupling capstone** (#28): the sub-daily loop coupled and owning its own state.
- **Prognostic leaf and wood energy** (#41) with a separate wood temperature; FAST diurnal
  diagnostics in the same change.
- **Snow and temporary surface water** (#42): stateless kernels and a conserving fast-loop
  coupling that closes whole-column mass and energy budgets to machine precision through
  accumulation, sublimation, melt into infiltration, and the snow-albedo ramp.

### Numerics

- **IMEX-ARK** (2026-07-11, e77f7f0, no PR number): an L-stable ESDIRK2 on the ARS(2,2,2)
  tableau with an arrowhead Newton surface solve. The config string stays `"ark"`; the explicit
  part is empty, so despite the name it is a diagonally implicit scheme.
- **Error-control infrastructure** (#65) and tolerance unification, a process mask, a sweep
  harness (#66).
- **ED2-faithful RK45** (#67, #68): an adaptive Cash-Karp 5(4) march, internal water carried as
  **mass** rather than potential, and a mass-conserving demographic seam.
- **Integrator parity** (#77, #80–#87): a long sequence making the schemes solve the same
  physics. Along the way: the reference used to score them was mis-timed; ARK and RK45 dropped
  snowfall all winter; `veg_coupling_floor` destroyed energy; the split scheme's soil-energy
  budget was a tautology; the pond held no enthalpy so water crossing into it shed its heat
  while the books still closed.
- **The operator-split integrator was retired** (#88). It converged to a different limit
  (~0.45 K in canopy-air temperature) that refinement never removed and nobody ever attributed,
  and it could not carry the coupled tissue heat store. `ark` became the default and `rk45` the
  accuracy baseline, with the RK45 stiff rescue redoing the step on `ark`. The tissue heat
  store was turned on in the same change, integrated by an exact exponential with two weights
  (endpoint and step-average) because the tissue ODE is linear under the frozen coefficients.

### Slow loop

- **Phenology**: a stateless signal module (#10), wired into the run loop (#39), then refactored
  to a **rate-based** signal-only kernel emitting two relative tendencies instead of a
  directional tri-state (#51).
- **Plant carbon**: pure carbon-dynamics kernels (#14), carbon-pool allometry and PFT traits
  (#15), then an **elemental growth-allocation kernel** following FATES PARTEH Hypothesis-1
  (#52), with growth respiration charged inside the kernel on realized growth.
- **Trait plasticity** (#54): light-driven acclimation of specific leaf area, V_cmax, dark
  respiration and leaf lifespan.
- **Non-leaf maintenance respiration** (#13): stem and fine-root, an ED2 Chambers-2004 port.
- **Slow soil-carbon biogeochemistry** (#35): ED2's CENTURY decomposition reorganized as the
  carbon matrix ODE `dX/dt = B·I + A·ξ·K·X`, with a 7-pool state, a lignin sub-tracer, an exact
  augmented matrix exponential for accelerated steps, and a SASU steady-state solve. Wired into
  the slow loop (#64): the fast loop's heterotrophic respiration respires the same frozen pool
  the daily step debits, so the day's total fast Rh equals the daily debit by construction.

### Forcing

- **Meteorological forcing** (#36): a single-site NetCDF reader over a multi-grid `(time, grid)`
  file, `pure`/`elemental` disaggregation kernels including an interval-mean-conserving
  shortwave reconstruction, the canopy-RT join, net longwave, Weiss-Norman band-specific
  shortwave, multi-year calendar recycling, nearest-grid matching, and wind-height and
  elevation lapse. MEDS never gap-fills: a missing or NaN required value is a hard error.
- **The recycle window became declared rather than inferred** (#69). The previous classifier
  accepted only Jan-1 00:00 files and silently fell back to an absolute-seconds span wrap
  otherwise. Real ERA5-Land records are stamped at the *end* of each interval, so their first
  record is 01:00:00 and they always took that fallback: the Ithaca file spans 366 d 22 h, and
  a 29-year run ended up reading late May at a ~10 h offset. Nothing caught it, because the
  cosz reconstruction is mean-conserving, so daily-mean shortwave stayed correct and the slow
  demography looked healthy.

### Output and I/O

- **Diagnostic aggregation and output subsystem** (#38, #40): the registry, the per-tier
  temporal integrators, and the netCDF serializer, written through the netCDF **C** library via
  `iso_c_binding` so the output layer builds under ifx and nvfortran (netCDF-Fortran's module
  format is gfortran-only).
- A **FAST (sub-daily) output tier** (2026-07-13, 01347bb, no PR number) for diurnal-cycle
  analysis.

### Code review and bug fixes

- **Adversarial code review, 2026-07-06**, over ~10.4k lines and 47 modules, targeting
  physical-process bugs, numerical defects and organization. All sections addressed across
  #29–#34: 8 critical/high physical-process bugs, then medium, then low-priority guards, then
  performance and solver issues, then organization.

### Tooling and infrastructure

- CMake with automatic Fortran module-dependency resolution — the deliberate fix for ED2's
  "run `make` six times" hack and per-platform `include.mk` files.
- **nvfortran portability trap documented** (#8, issue #7): never pass an array-valued function
  result straight into a call. nvfortran's whole-program optimizer miscompiles the temporary
  descriptor — silently wrong values at `-O2`, segfault at `-O0` — while
  `ifx -stand f18 -check all` tolerates it. A green ifx run is not sufficient.
- 3D visualization of forest structure (#4); the biophysics example (#70).
- Design plans relocated to `docs/dev_plans/` (#48); GitHub math rendering fixed in the science
  pages (#50).

---

[Unreleased]: https://github.com/xiangtaoxu/MEDS/compare/v0.2.2...beta
[0.2.2]: https://github.com/xiangtaoxu/MEDS/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/xiangtaoxu/MEDS/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/xiangtaoxu/MEDS/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/xiangtaoxu/MEDS/releases/tag/v0.1.0
