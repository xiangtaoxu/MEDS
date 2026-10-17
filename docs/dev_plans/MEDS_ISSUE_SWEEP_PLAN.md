# MEDS open-issue sweep — plan

**Revision 3, 2026-10-17. Every decision is taken (§6); work under way, decisions made during it are logged in §8.**
- Revision 1 read every open issue and proposed three phases and decisions D1–D11.
- Revision 2 recorded the owner's decisions (§6).
- Revision 3 records the last four:
  - D9: a leafless leaf takes its wood's temperature;
  - D12: a note in the READMEs of the runs not repeated;
  - #379's new turnover form is deferred;
  - the wind kernel's crown argument goes.

**Against:** `beta` at `38c715b` (after the annotation review, #409). Every open issue was read with
its comments, and each issue's premise was checked against that commit's code. Where a quick
read-only measurement settled a claim, it was made: ex06's region state files, the ERA5-Land
orography, a Python copy of the leaf kernel. Numbers marked *estimate* are analytic, not model runs.

**Scope.** All 53 open issues (#74–#425), plus 11 defects the review found that are not filed (§5).
The owner's priorities: **bug fixes and documentation fixes first**; for enhancements, a recommendation
of what to add now and what to defer. The sweep is meant to leave MEDS clean before the move to
ED-MORE (`MEDS_TO_ED_MORE_TRANSITION_PLAN.md`, draft); §7 says where this plan departs from that
draft's issue dispositions. **No long-term run is repeated in this sweep** (D6, D7): example 05's
regeneration and regrowth and example 06 are a separate job after it.

---

## 1. Summary

An issue can sit in two rows (a doc line now, its fix later).

| Verdict | Count | Issues |
|---|---|---|
| **Fix now: bugs and cleanups that leave results unchanged** (Phase 1) | 7 + 6 new | #394, #412 (part 3), #416 text, #417, #418 (nitrogen removed, D5), #420, #421; new: Campbell seed, missing-file fallback, `meds.canopy` reopen (D4), region resume without a failed cell (D8), crown area retired (D1), example paths and hygiene |
| **Fix now: documentation** (Phase 1, in parallel) | 6 | #422 (banner only), #423 (stop-gap only), #424 (relabel, D6), #425, #377 (answer, then close), #350 (one claim); plus the doc lines the other fixes need |
| **Measure first** (Phase 2) | 7 | #395, #413, #415 (parts 1 and 3), #419, #356, #265 (D10), #74 (size only) |
| **Fix now: bugs that move results** (Phase 3, one batch, contemporary runs only) | 7 + 1 new | #342, #349 (D9), #393 (D2), #395, #410 (CLM form, D1), #416, #419, plus new: the cohort stuck at the plant-number floor |
| **Enhancements to add now** | 2 | #401 (leaky bottom), #379 (the rename now, D3) |
| **Close** | 3 | #183 (shipped), #377 (answered), #412 part 1 (nothing to respire) |
| **Defer, with a trigger** | 27, and parts of 5 more | §7; #369 joins #343 (D11) |

In time: about **4 days** for Phase 1 with its docs, **3 days** of measurements, and **4–5 days**
for Phase 3 with its reruns. That fits the transition draft's bug-sweep window (Oct 26 – Nov 6) and
could start on Oct 19. The release that closes the sweep is v0.3.5.

---

## 2. Phase 1 — fixes that leave results unchanged, and the documentation (≈ 4 days)

Each PR here must be byte-identical on every shipped config and on the `rr` region case, on ifx
Release, ifx Debug and gfortran. The doc PRs can be reviewed in parallel with the code PRs.

### PR 1a — correctness that a valid config never sees

| Item | What is wrong | Fix |
|---|---|---|
| **New: Campbell seed** | Under Campbell retention, `meds_fast_prepass.f90:247-248` passes `vg_alpha(1)`/`vg_n(1)`, which Campbell never sets (both 0). The predawn-potential seed of a new cohort reads 0 MPa (a fully wet soil), and `dse_dhead` divides by zero, which should trap at once in a Debug build (inferred from `-fpe0`, not run). The water network itself uses the right accessors. | Give `soil_params_t` one parameter pair and drop the van Genuchten–named fields, so no caller can pick the wrong pair. A Campbell column test in `test_spac_solver`. |
| **New: missing census or restart file** | A site run whose `census_file` or `restart_file` cannot be opened silently starts from bare ground; the only notice needs `verbose`, and Python runs (`verbose=False`) print nothing. The fast calibration restarts every trial from a state file, so a missing file becomes a plausible wrong trial. `docs/configuration.md` documents the fallback. A region stops, but says "not a region state file" for a file that does not exist. | Stop with a message naming the mode and the file; bare ground only for `init_mode = 0`. Give the region its own "cannot open" message. Refuse a region state file in a site run (rank-2 variables read as rank 1, not checked). A ctest per mode. |
| **New: `meds.canopy` reopen** (D4: keep `meds.canopy`) | `Canopy(path)` skips reopening a path it already has open, so a config rewritten at the same path is silently ignored, which is what a calibration loop does. | Cache by instance, not by path; a test that rewrites a config at the same path. |
| **New: resume a region without a failed cell** (D8) | `meds_main` stops a whole region at the first failed polygon, and a region restart must hold exactly the region's polygons. So a 600-cell run whose cell fails in year 120 cannot mask that cell and resume. | Let a region restart from a state file whose polygons include the region's: the region's polygons are a **subset** of the file's, matched by cell. Masking the failed cell (`[region].cell_mask`) and resuming then works. Byte-identical when the two sets are equal. A test in `test_region`: resume with one cell masked out. |

### PR 1b — clean-ups (byte-identical)

| Item | What |
|---|---|
| **#417** dead code | `size_tol` (unused); `area_tol` replaces the literal `1.0e-5` in the two finalize checks; the `coh_transp` field (written, never read); `apply_rt_forcing`'s unused `tcas` and the call that computes it; `veg_energy_balance`'s `t_emit` and its `(t_cas − t_emit)` terms, always zero; `meds_run_t%is_open`. **164 unused imported names in 35 modules**, removed by a script; a lint like `clamps_named` is optional. |
| **#418 + D5** nitrogen and inert keys | Retire `[soil_carbon].agf_fast`/`agf_struct` and the five nitrogen keys through the retired-key mechanism, with a message naming what replaces them (nitrogen: "not built; planned after the move to ED-MORE"). **Remove the nitrogen from the code:** the always-zero nitrogen twin of the soil pools and of the litter input, its area blends, `n_cycle_on`, the `f_decomp` placeholder, and the science page's nitrogen paragraph (`soil_carbon.md` §10 says what is not built). ROADMAP #154 says nitrogen and C–N coupling come after the move. |
| **D1** crown area retired | The crown-area allometry is never called: retire `[allometry].ca_b1`/`ca_b2` (required keys today) and `dbh_to_crown_area`. Retire the PFT constant `crown_area_frac` too, since it is 1.0 in every shipped config and its only use is the in-canopy wind's crown cover. The wind kernel's crown argument goes with it, and so does `test_aerodynamics`'s "small crown" test; eq. (8) of `canopy_aerodynamics.md` becomes the closed-crown form u·exp(−LAI/2). |
| **#420** `clamp_soil_energy` | Its reason (ground evaporation's `pow()` inside the RK45 stages) is gone, it never runs on ARK, and its `denergy` output is never used. Count how often it fires on the suite and on a forced RK45 year; delete it with its unused imports if zero, or name the stage term that needs it. (It does not create energy off the books, as the issue feared: it clamps a throwaway stage input.) |
| **#421** tests | Delete the duplicate `test_rk45_shed_water` in `test_column_ark.f90` (the RK45 file's copy is the superset); rename `test_mass_correction_neutrality` and relabel its checks (it still guards the liquid-enthalpy datum). |
| **#394** pool floors | The applier's `max(pool + rate·dt, 0)` becomes a round-off snap with anything larger reported by the bounds guard; `reacclimate_plant_traits` takes extra leaf carbon only up to the storage there is. The floors never bind today (worst 2.5×10⁻¹³ kgC m⁻²), so results stay bit for bit. |
| **#412** part 3 | `stem_resp_size_scaler` is a module default no key can set, applied as `10**(0·dbh)` every step: delete it, drop the factor from eq. (3), and correct `plant_respiration.md` §7, which says shipped runs have no size scaling (ex05 and ex06 scale stem respiration through the sapwood basis). |
| small | Rename `store_hcap_per_dt` (the storage refill in the allocation kernel, misnamed in an old rename pass) to `storage_refill`. Default `snow_stage_t%t_melt` to the triple point, not 0 K (N8). |

### PR 1c — examples and repository hygiene (byte-identical)

- Every example runs without `/ibstorage`. Affected: `example04/bci_census.toml` (two census paths),
  `example05/palo_verde/meds_config_inventory.toml` (`data_path`), `example02/make_era5_par_hours.py`
  (default archive). Each takes the path from an option or the environment and stops with a
  message when it is missing, as example 06 does.
- Delete the root `meds_output_pft_parameters.csv` (untouched since 2026-07-31, read by nothing).
- `CLAUDE.md` says "53 tests"; it is 72 (73 with the Python library). `README.md` says 252
  variables; it is 284.

### PR 1d — documentation (no code)

| Issue | What |
|---|---|
| **#422** | Deferred by the owner. Meanwhile, a short **banner** at the top of `ed2_comparison.md`, the page `docs/README.md` sends every ED2 user to first. It says the page describes v0.3.0 and names the statements that are now wrong (the capacity limb is always on; leaf water potential converges at 900 s; per-cohort roots and sapwood conductivity; Chave height, soil maps, threaded regions), with links. |
| **#423** | Deferred by the owner. Meanwhile, `post_proc/README.md` says which scripts read `meds_main` files (`plot_pft_size.py`, the notebook) and which read only example 03's single-file layout (the other four), with command lines that work. **README.md's Quick start is broken the same way:** it runs example 03's config, which writes no diagnostics, then plots a file no run writes. Point it at a config that writes output, or end it at `ctest`. |
| **#424** (D6: relabel only) | Palo Verde's PFT file labels three values "example 04's calibration": they are from example 04's fit on the 1.4 m column (`5566b85`); example 04 was refitted on 2 m after that. Relabel them as the Palo Verde set's own values. The `g1` comment says Lin's slope −3.97 per g cm⁻³ where the values follow it scaled (−2.70). The `g0` comment cites Duursma et al. (2019) for a value ten times below theirs. |
| **#425** | The θ_J table in `leaf_gas_exchange.md` §2 cannot be reproduced on today's kinetics. A Python copy of the kernel gives, at saturating light and Ci 280: −27.4 % (θ = 0.85 on both transitions), −23.7 % (on Ac/Aj only), −9.8 % (the shipped 0.98/0.95), independent of Vcmax; the "30–32 %" is the net rate at finite light (−31 to −34 %). Re-measure with `meds.plant.leaf`, state each number with its conditions, and align the two code comments. §4.4 quotes a stale `g0` as well. |
| **#377** | Answer, then close: a constant absorbing fraction of root length is exactly a rescaled specific root length (both L and L_v scale), so it adds no freedom. One sentence at `plant_hydraulics.md` eq. (14); the sub-linear limit belongs to #356. |
| **#350** | `MEDS_FAST_CALIBRATION_BEST_PRACTICE.md:122` says the data report prints the share of area whose canopy-air top is above the sensor; it does not. Add the print or correct the line. |
| others | `vegetation_energy_dynamics.md:137-138` says the sapwood ring is the hydraulic store (it is sapwood plus fine roots; #419). The example 04/05 configs say the canopy intercepts "rain and dew" (dew goes to the soil; #74). `numerical_scheme.md` §5a′'s table (Φ′ = −1.16 at 150 s for the 30 m stand) contradicts §5a's "stable at every `dt_fast`" (resolved by #415's measurement). `configuration.md`: a failed polygon stops the region after its month. ROADMAP: #269's row repeats a retracted figure; #265/#268/#269 sit under "Known defects" but are modelling choices; #416's row text; the stale example 06 bullet; the archived forcing-design path; #342's "before high-elevation sites" is already past (below). The polygon plan header still says "implemented on feat/region-runtime". |

---

## 3. Phase 2 — measurements that decide the rest (≈ 3 days)

Read-only runs: no code is merged. Each result goes into its issue as a comment.

| Issue | Question | Measurement |
|---|---|---|
| **#395** | How much water does `reconcile_tissue_water_capacity` create (seeding recruits) and destroy (clamping a store above today's capacity), off both ledgers? | Return its two totals per year at BCI (ex04, ex05's census check) and Palo Verde's inventory run. |
| **#413** | Where is the whole-column energy residual under a melting pack (≤ 5×10⁻³ W m⁻²)? | In the snow fixture, give the tissue store and the canopy film per-step ledgers (neither has one) and print each interface's enthalpy on both sides. A second lead: the film is valued at the current step's inflow temperature, so a stored film is silently revalued between steps (N7). |
| **#415** (1), (3) | Is the production `dt_fast` = 900 s safe in light wind and through a drought season, on the water network? | From ex04 (BCI) and Palo Verde inventory restarts (a Dec–May dry season): 150/300/450/900 s, scoring canopy-air temperature variance at the 2·dt period binned by wind, GPP, ET, leaf and wood ψ. Re-measure Φ′ on today's ARK at 2 and 10 m/s and fix §5a/§5a′. Part (2), the midpoint freeze, is new numerics: defer. |
| **#419** | How much does the wood's heat capacity change if its water is the sapwood share of the store the water network carries, instead of a fixed `bsap × 2 × 1.0`? (At default traits the fixed value is twice the saturated sapwood water; Palo Verde's PFTs carry 0.99–4.41.) Also: is `c_sapw` = 2700 J kg⁻¹ K⁻¹ a fresh-wood value (then the water is counted twice)? ED2 uses 1159 for oven-dry wood. | Wood temperature amplitude, sensible heat and closure on the ex04 window and a Palo Verde dry season. |
| **#356** | Is root uptake still unlimited in wet soil on the network, where the Gardner conductance is 30–230× the old Katul one? | ex04 evaluation run over 2017-04-15…20 (the storm front): the largest cohort's peak uptake over its peak sapflow. Above ~2×, a series root conductance is worth designing (defer the build either way: it needs refits). |
| **#265** (D10) | After #410, does the CLM4 dense-canopy ground conductance still matter? (It is a faithful port of a published scheme, not a bug; ED2's default is the Sellers integration, about 10× less resistive at BCI midday by *estimate*.) | ex04's 10-day window and the Palo Verde inventory run at `cs_dense` 0.004 and 0.04: H, LE, soil temperature and water, LAI. Within noise: close as a documented choice; otherwise port ED2's default with a refit, as a later item. |
| **#74** | Is dew worth routing to the leaves? | The yearly condensate total at BCI against the yearly film evaporation. |

---

## 4. Phase 3 — bugs that move results (≈ 4–5 days with the reruns)

One PR per item, each with its before and after numbers, merged in this order so the contemporary
runs are repeated **once** at the end.

| Order | Item | Fix | What moves |
|---|---|---|---|
| 1 | **#416** `grav_head` | `grav_head = rho_h2o · grav · 10⁻⁶` (the soil side uses `grav_head` too, so the issue's text is corrected). | Heads and potentials, ~3×10⁻⁴ relative. |
| 2 | **#342** Γ* | Divide by a named reference pressure (100 kPa, Bernacchi's), not the local one; one kinetics helper replaces the three copies of the conversion (kernel, a Python test, example 01). A test that Γ* in ppm does not depend on pressure. Example 06's domain needs it: 29 % of its cells are above 1000 m, and Γ* is 23 % too high at the 90th-percentile cell. | GPP: BCI +0.1 %, highland cells several %. |
| 3 | **#349** coupling floor (D9) | Replace the 1 W m⁻² K⁻¹-of-ground floor by a divide guard and drop the slaving conductance: the always-on heat store already bounds the step. Make `veg_hcap_min` per unit of tissue area (N9), so a tissue's heat capacity and all its fluxes scale with its area. A tissue with no area keeps no heat and takes the temperature of what it grows from (D9): a leafless leaf its wood's, a cohort's wood with no area the canopy air's. A split-invariance test: one cohort against two at half density. | Young and sparse stands (GPP −2.6 %, LAI −2 % over 5 years from bare ground); BCI 0.00–0.04 %. |
| 4 | **#410** open canopy (D1: CLM) | Each patch's open fraction W = exp(−(LAI + WAI)) blends the bare-soil and dense-canopy ground conductances, as CLM does, so a sparse or empty patch couples its ground like bare soil. The snow-burial half waits for a snowy site. | Ground H, LE, soil temperature and water in gaps and sparse patches; tiny at closed BCI. An empty patch is ~28× too stiff today (*estimate*). |
| 5 | **#393** starving cohorts (D2), **N1** the floor cull | No new hazard. **A special case:** a cohort with no leaves and no storage dies within one month, through the existing mortality path, so its carbon goes to litter. Its mortality rate is set so its density reaches the cull floor in 30 days. Example 06's final states hold 7.7 % (region 1) and 9.6 % (region 2) of cohorts in that state. ROADMAP's #393 entry says this is interim, to be revisited with a growth- or storage-based mortality. Fix the cull at the same time: the applier floors plant density at exactly `negligible_nplant` while the cull keeps `>=` it, so a cohort that decays to the floor lives forever (24 such cohorts in example 06, dbh ≥ 19 cm). | Runs with leafless, storage-less cohorts (Palo Verde, example 06's domain); BCI unchanged. |
| 6 | **#395** tissue water | If Phase 2 shows it matters: the clamp's excess goes to the ground through `shed_water_rate` and is declared, as turnover's already is; recruits are born with their starting water, declared in the recruitment phase, which removes the zero-water sentinel. A ledger test across a recruitment and a snap-to-bare. | Soil water slightly; slow-ledger water terms. |
| 7 | **#419** wood heat capacity | If Phase 2 confirms: the wood's thermal water is the sapwood share of the water-network store; one function for the fast loop and the fusion weights; the duplicate constants go. | Wood heat capacity −10 to −20 % at default traits, up to ~2× for Palo Verde; small flux changes. |
| 8 | **#401** leaky bottom *(enhancement, added now)* | One drainage factor f in [0, 1] replaces the three cases: free drain f = 1, leaky f = `f_leak` from the map, bedrock f = 0; the aquifer stays its own case. Region files carry `f_leak`. Tests: f = 1 and f = 0 reproduce today's free drain and bedrock bit for bit; a short map-driven region run (`test_region`). | Site runs unchanged; map-driven runs (82 % of example 06's cells are leaky and drain ~20× too fast today), seen in the separate long-run job. |
| 9 | **#379** turnover as a rate (D3: rename now) | Every lifespan name goes, and `leaf_turnover_rate` is the only one. Today the code still has `leaf_lifespan_toc` (the PFT key, 31 lines) and the cohort trait `llspan` (83 lines in 23 files, written to the state file); `leaf_turnover_rate` exists only in the phenology kernel. Rename the key to `leaf_turnover_rate_toc` [1/yr] (the old key is refused with the conversion in its message) and the cohort trait to a rate; the light plasticity becomes rate = rate_toc·exp(−k·LAI); the trait fuses on leaf carbon, which conserves the turnover flux (a leaf-area mean of a lifespan does not); a restart reads an old `llspan` as 1/x; the related keys (`kplastic_llspan`, `fineroot_lifespan_plastic`, the `[derived]` override) follow. The issue's new turnover form (Q1) is deferred: it would change example 02's and example 06's seasonality, and neither is refit or rerun in this sweep. | The rename: round-off to small (the fusion weight). |

**The rerun (D7: contemporary runs only).** After the batch: ctest on three compilers, then:

| Example | Rerun? |
|---|---|
| Ex01–ex04 | Every run, as they stand. |
| Ex05 | The BCI census check and the Palo Verde 2009 inventory run. |
| Ex05 regeneration (1600–2020) and regrowth (150 years) | No. |
| Ex06 | No. |

The outputs of the reruns are archived as v0.3.5. Example 04 is refitted only if BCI moves beyond the fit's noise; the expected moves there are ≤ 0.1 %. **The examples not rerun keep figures made with the code before the sweep.** Their READMEs say which release made them, one line each (D12), until the separate long-run job repeats them.

---

## 5. Defects found by this review that are not filed

Each would get a self-contained issue when its work starts, closed by its PR.

| | Defect | Where it goes |
|---|---|---|
| N1 | A cohort at exactly the plant-density floor is never culled. | Phase 3, with #393 |
| N2 | Campbell seed reads van Genuchten parameters. | Phase 1, PR 1a |
| N3 | A missing census or restart file falls back to bare ground silently; a region state file in a site run is not refused. | Phase 1, PR 1a |
| N4 | `meds.canopy` ignores a config rewritten at the same path. | Phase 1, PR 1a |
| N5 | Examples need `/ibstorage`; stale root CSV; wrong test and variable counts; README Quick start. | Phase 1, PR 1c/1d |
| N6 | A region cannot resume without a failed cell. | Phase 1, PR 1a (D8) |
| N7 | The canopy film is revalued between steps at the new inflow temperature, unbooked (≈ 0.34 MJ kg⁻¹ when a pack forms on a wet canopy). Code reading only. | Phase 2, with #413 |
| N8 | `snow_stage_t%t_melt` defaults to 0 K; a trace pack in a restart would value its shed water at 0 K. | Phase 1, PR 1b |
| N9 | `veg_hcap_min` (20 J m⁻² K⁻¹) is per ground area, the same split-dependence as #349. | Phase 3, with #349 |
| N10 | The crown-area allometry is required config but never called. | Phase 1, PR 1b (D1) |
| N11 | `[soil_carbon].xsand`/`xclay` (fixed) drive CENTURY5's texture terms even when soil maps give the texture. Only `decomp_scheme = "century5"`. | Defer, or Phase 1 if trivial |

Also worth asking: example 06's region 1 stopped logging after 1946-01-01 and was resumed from its
state file (`meds_runs/ex06_region/region01_stalled_1946`); if that was a MEDS hang, not a node
fault, it is unfiled.

---

## 6. Decisions

### Taken by the owner, 2026-10-17

| | Decision | Owner's answer |
|---|---|---|
| D1 | #410's open fraction | CLM's W = exp(−(LAI + WAI)); retire the crown area. |
| D2 | #393 | No separate starvation hazard: growth already sets mortality. A cohort with no leaves and no storage dies within one month (a special case), with a ROADMAP note to revisit. |
| D3 | #379 | Rename now. |
| D4 | `meds.canopy` | Keep it (and fix N4); a radiative-transfer example comes later. |
| D5 | Nitrogen | Remove it now; nitrogen and C–N coupling come after the move to ED-MORE. |
| D6 | #424 | Relabel only. No example 06 rerun in any phase. |
| D7 | Rerun scope | Contemporary runs only; example 05's and 06's long-term runs are a separate job after the sweep. |
| D8 | Region resume | Include it if small: it is (PR 1a). A region resumes from a state file holding its polygons and more, its polygons being a subset of the file's. |
| D10 | #265 | Measure, then decide. |
| D11 | #369 | Defer the PAR-share investigation and change; test it across flux towers together with sunlit/shaded leaves (#343). |

### Taken after revision 2

| | Decision | Owner's answer |
|---|---|---|
| D9 | A leafless leaf's temperature | The wood's: a bud sits on its twig. |
| D12 | Note the runs not repeated | Yes: one line in each affected README. |
| — | #379's new turnover form | Defer. |
| — | The wind kernel's crown argument | Remove it with the crown area. |

**D9, worked out.** NaN was considered and rejected. Every reader of a leaf temperature weights it,
and NaN times a zero weight is NaN:
- the two-stream's emission temperature, which would turn the cohort's wood emission into NaN;
- the output means and the cohort fusion;
- the restart file;
- the run's NaN guard, which stops on a NaN tissue temperature.

Each would need a special case. Once #349's floor goes and N9 makes the heat-capacity minimum per unit
of tissue area, a tissue's heat capacity and every flux scale with its area. A small leaf is then simply
a leaf, and only exactly zero area needs a definition. The rule is that **a tissue with no area takes
the temperature of what it grows from**:
- a leafless leaf takes its wood's (the owner's choice: a bud sits on its twig);
- wood with no area (a grass, or a seedling without wood area) takes the canopy air's.

Nothing weights that value, because its area and heat capacity are zero. New leaves reach their own
balance within a step, so the choice barely moves results. Zero heat capacity holds zero heat, and
the thermal mass regrowth adds is declared as growth's is today.

---

## 7. Enhancements: add now or defer

**Add now:** #401 (small; makes the soil maps mean what they say; consolidates three bottom
conditions into one), #379's rename (D3).

**Defer**, with what would trigger each. All are new science or capability, none is needed for a
correct, clean MEDS, and each would be re-filed in ED-MORE:

| Issue | Trigger |
|---|---|
| #74 dew on leaves | #165, and Phase 2 showing dew matters |
| #96 Kelvin leaf humidity | a foliar-uptake study |
| #154 nitrogen, #155 layered soil carbon | after the move to ED-MORE (D5); a deep-carbon or radiocarbon application |
| #156 coarse woody debris, #157 fire | scoring a post-disturbance carbon budget; a savanna or dry-forest fire study (example 06 is the first candidate) |
| #165 canopy film heat, #186/#187 snow | the first boreal or alpine site; #413 located first |
| #181 hydraulic redistribution, #356 root conductance, #377 | designed together (one soil-to-root edge); #356's Phase 2 number |
| #196 cohort threads, #400 idle threads, #407 forcing read-ahead | a site run limited by its patch count; small regions becoming common; the first long transient regional run |
| #254 deep thermal layers | a soil-temperature-at-depth target |
| #255 spectral record, #411 soil albedo, #414 NIR albedo | a band-reflectance target; leafless or sparse canopies scored; the radiative-transfer example (D4) for #414 |
| #256 acclimation | a temperature-gradient study (the per-cohort half has no leaf-temperature fit to use) |
| #257 longwave cloud term | the first run on a source without LWdown |
| #268 litter layer, #269 layered canopy air | re-measured per patch after #410 (both cases were gap-patch site means) |
| #302 forcing adapters | a study needing finer forcing |
| #310 (R3) one run container | a region Python entry (none planned). Its R10 has dissolved: close that part. |
| #343 sunlit/shaded leaves, #369 PAR share | one study across flux towers (D11). **Save the #343 prototype first**: it exists only as uncommitted changes in the `sunlit-shaded` worktree. |
| #350 neutral move to the canopy-air top | calibrating at a short-vegetation tower |
| #379's new turnover form | the next refit of example 02 and rerun of example 06 |
| #381 size-dependent parameters | ex05's growth residual by size surviving light plasticity |
| #390 film flow in dry soil | a dry-season ψ target |
| #393's full mortality | revisit the one-month special case with a growth- or storage-based mortality (D2) |
| #412 part 2 (substrate-limited respiration) | after #393's revisit (they pull in opposite directions) |
| #415 part 2 (midpoint freeze) | new numerics |

**Close:** #183 (the region runtime shipped in #408; its tails are #400, #401, #407 and #310),
#377 (answered, §2), #412 part 1 (reproductive carbon is never held in a pool, so there is nothing
to respire; the issue stays open for part 2).

**Where this departs from the transition draft (Appendix D).** The draft keeps #393 and #401 out of
the bug sweep because they change the examples. This plan takes both in: #393 is a bug that keeps
dead plants in example 06's reported density (the one-month special case), and #401 is small and
makes half of example 06's cells drain as their map says. Neither forces a rerun here (D6, D7). The
draft's "comment rewrite" and its lint test are done (#409), and the region runtime is merged (#408).

---

## 8. Log — decisions made during the work (the owner was away from 2026-10-17)

The owner asked for the work to continue through the phases without waiting, and for every decision
that differs from this plan to be recorded here.

| | When | Decision | Why |
|---|---|---|---|
| L1 | Phase 1 start | Nothing is merged into `beta` while the owner is away. Each piece of work is a PR into `beta`; later phases build on an integration branch that merges the earlier PRs, so they can be reviewed and merged in order. | Merging into `beta` was always the owner's call. |
| L2 | Phase 1 start | Phase 1 runs as three parallel PRs: 1a (correctness), 1b (clean-ups) and 1c+1d (docs and example hygiene, one PR). The tree-wide unused-import removal (the rest of #417) comes after them, to avoid conflicts in 35 modules. | Parallel work in separate worktrees; one broad mechanical change last. |
| L3 | Phase 1 start | New issues filed: #427 (N2, Campbell seed), #428 (N3, missing file), #429 (N4, `meds.canopy`), #430 (N5, example paths and hygiene), #431 (N6, region resume), #432 (N10, crown area). N8 (`t_melt` default) has no issue of its own; it rides with #417's PR. | §5: each defect gets an issue when its work starts. |
| L4 | Phase 2 start | Phase 2 runs alongside Phase 1, except #265, which needs #410 and is measured on Phase 3's branch. | The measurements read beta's behaviour, which Phase 1 does not change. |
| L5 | Phase 2, M1 | **#395 is negligible in size** (after the census fill, 1e-5–1e-3 kg m⁻² yr⁻¹, under 1e-5 of transpiration), but it stays in Phase 3 as a bookkeeping clean-up with results unchanged to ~1e-5. The routine's own outputs also lack the patch-area weight (they overstate by roughly the patch count); that is fixed in the same PR. | The plan said "if Phase 2 shows it matters"; ledger cleanliness is worth the 20–40 lines. |
| L6 | Phase 2, M1 | **#74 closed** as measured and negligible: dew is 0.01–0.5 % of ET at BCI (≤ 5 mm yr⁻¹). The plan had it deferred. | Nothing to carry to ED-MORE. |
| L7 | Phase 2, M1 | **#356 measured:** after the 2017-04-17 storm, big trees take up 13–14× their peak sapflow (≈ 1 t in 15 min) from a −0.5 MPa drive; site ratio 4.1. The build stays deferred as planned (it needs refits); the design is recorded on the issue. | As planned. |
| L8 | Phase 2, M1 | **New issue #433:** under ERA5-Land forcing at BCI, intercepted-water evaporation is ~34 % of rain and ~75 % of ET (11 % of rain on the tower forcing). A question for later, not this sweep. | Found while measuring #74. |
| L9 | Phase 2, M4 | **#419 confirmed, and widened:** `c_sapw` = 2700 J kg⁻¹ K⁻¹ is a fresh-wood value (100 % moisture), applied to dry mass with the water added on top, so the water is counted twice. Phase 3 takes the wood's thermal water from the water network's store (sapwood share), sets `c_sapw` to ED2's dry 1217, and counts heartwood dry. Midday H at BCI +20–27 W m⁻², night H −4 to −11, daily means unchanged. **Example 04 is refitted** after the batch, since its diurnal H is in the calibration. | The plan refitted example 04 only if BCI moved beyond noise; this does. |
| L10 | Phase 2, M4 | **New bug #434, fixed in Phase 3 with #419:** fine-root turnover drains ~10× too much water from the wood store (≈ 3.6 kg m⁻² per 9 months at BCI against ≈ 0.4 held by the fine roots). | Same code and same wood-water accounting as #419. |
| L11 | Phase 2, M3 | **#413 was already fixed** (2026-09-08, in #119: the condensate's enthalpy is valued at one temperature on both sides); the "open" note came from a squash message. With the old valuation put back the residual reappears exactly (0.63 J m⁻² per step) and the suite's RUN 9b fails. Phase 3 opens with a byte-identical commit: the energy floor back to 1e-3 W m⁻² (every test passes with it), the floor's comment and the ROADMAP item corrected; #413 closes with it. | Measured; nothing to hunt. |
| L12 | Phase 2, M3 | **N7 measured and moved to #165:** the film's between-step revaluation is real but nets small over a season (−3.2e4 J m⁻² in 90 days at Ithaca); M3 also found the film stranded under a snow pack all winter. Both are fixed by #165's paired film mass and enthalpy, deferred to the first snowy site. The plan's 0.34 MJ kg⁻¹ was wrong: 0.91 on MEDS's datum. | Same mechanism as #165. |
| L13 | Phase 2, M3 | **New issue #435:** the canopy-air capacity changes between steps without booking the exchange (ED2's "zcan effect"); low priority, not in this sweep. The tissue store's between-step jump (fixed wood water; the leaf's 0 K datum) is recorded on #419 and handled in its Phase 3 PR. | Found while measuring #413. |
| L14 | Phase 2, M3 | **Test bug folded into #421:** RUN 9's "SNOWF r45" leg never rebuilds the integrator options, so it runs ARK; the clean-up PR fixes it. | Test-only, same PR as the other test residues. |
| L15 | Phase 2, M2 | **#415 (1), (3) measured:** 900 s is stable in every wind band (no period-2 signal; BCI means within 0.1 % of 150 s), and §5a′'s Φ′ = −1.16 is obsolete (|Φ′| ≤ 0.11 at 2 and 10 m s⁻¹ on today's ARK). At Palo Verde, 900 s overestimates Feb–Apr GPP and ET by 5.5–6 %, linear in dt and caused by the per-step freeze, not by the water solver. **900 s stays the default**; Phase 3's docs state the drought bias, recommend 300 s for dry-season leaf-habit studies, and replace §5a′'s table. #415 stays open (deferred) for a predictor–corrector fix and its part (2). | The plan expected either a pass or a dt change; the answer is a documented bias. |
| L16 | Phase 2, M2 | **New bug #436, fixed in Phase 3 with #395:** re-flushing leaves start nearly empty of water, so the first step after each daily update reads ψ_leaf ≈ −10⁴ MPa (pollutes `psi_leaf_site`, closes the leaf for a step). The new capacity is filled at the plant's potential from the cohort's wood store, conserving the plant's water. | Same routine as #395. |
| L17 | Phase 1 | PRs: #437 (1c+1d, docs and examples), #438 (1a, correctness). In 1a the soil's parameter pair is named `curve_par_a`/`curve_par_n` after the config keys, with the accessor functions removed. The Campbell seed check sits in `test_column_rk45`, which has the column fixture, and the region state-file cases in `test_region`. 1c also ignores `*_pft_parameters.csv`, so a run's table is not committed again. | Implementation details; the plan named neither. |
| L18 | Phase 1 | **New bug #439:** `region_free` does not close an opened region's met source and output files, so reopening the same variable crashes. Production runs are unaffected. It is fixed in Phase 3 with #401 (same region code). | Found while building #431. |
| L19 | Phase 1 | The byte-identity harness had to be fixed mid-phase: it had compared an empty folder. Every Phase 1 result quoted in the PRs comes from the fixed harness or the docs agent's own comparator, which agree. | Correctness of the verification itself. |

## Appendix — every open issue

Type as found (B bug, D doc, E enhancement, Q question, C clean-up); verdict and phase.

| # | Title (short) | Type | Verdict |
|---|---|---|---|
| 74 | Dew to soil layer 1 | E | Defer; doc line now; measure size |
| 96 | Kelvin leaf humidity | E | Defer |
| 154 | Nitrogen | E | Defer to ED-MORE; the code's nitrogen removed now (D5) |
| 155 | Layered soil carbon | E | Defer |
| 156 | Coarse woody debris | E | Defer |
| 157 | Fire | E | Defer |
| 165 | Canopy film heat | E | Defer |
| 181 | Hydraulic redistribution | E | Defer |
| 183 | Multi-polygon runtime | E | **Close** (shipped) |
| 186 | Snow layers | E | Defer |
| 187 | Canopy snow | E | Defer |
| 196 | Cohort threads | E | Defer |
| 254 | Deep thermal layers | E | Defer |
| 255 | Spectral record | E | Defer |
| 256 | Acclimation | E | Defer |
| 257 | Longwave cloud term | E | Defer |
| 265 | Sub-canopy conductance | Q | Measure after #410, then decide (D10) |
| 268 | Litter layer | E | Defer |
| 269 | Layered canopy air | E | Defer; ROADMAP row fix |
| 302 | Forcing adapters | E | Defer; ROADMAP path fix |
| 310 | One run container | C | Close R10; defer R3 |
| 342 | Γ* and pressure | B | **Phase 3** |
| 343 | Sunlit/shaded | E | Defer, with #369 (D11); save prototype |
| 349 | Coupling floor | B | **Phase 3** (D9) |
| 350 | Neutral move | Q | Defer; doc claim now |
| 356 | Root conductance | E | Measure; defer build |
| 369 | PAR share | Q | Defer, with #343 (D11) |
| 377 | Root length question | Q | **Close** (answered) |
| 379 | Turnover as a rate | E | **Rename now** (D3); new form deferred |
| 381 | Size-dependent parameters | Q | Defer |
| 390 | Dry-soil film flow | E | Defer |
| 393 | Starving cohorts | B | **Phase 3**, one-month special case (D2); full mortality deferred |
| 394 | Pool floors | C | **Phase 1** |
| 395 | Tissue water booking | B | Measure; **Phase 3** |
| 400 | Idle threads | E | Defer |
| 401 | Leaky bottom | E | **Add now** (Phase 3) |
| 407 | Forcing read-ahead | E | Defer |
| 410 | Open canopy, snow burial | B + E | **Phase 3** (CLM form, D1); defer snow |
| 411 | Soil albedo | E | Defer |
| 412 | Respiration limbs | E | **Phase 1** part 3; close part 1; defer part 2 |
| 413 | Melting-pack energy residual | B | Measure (Phase 2) |
| 414 | NIR albedo | Q | Defer (with the RT example) |
| 415 | `dt_fast` checks | Q | Measure parts 1, 3; defer part 2 |
| 416 | `grav_head` | B | **Phase 3**; text fix Phase 1 |
| 417 | Dead code | C | **Phase 1** |
| 418 | Inert soil-carbon keys | C | **Phase 1**, with the nitrogen removed (D5) |
| 419 | Wood heat capacity | B | Measure; **Phase 3**; doc line Phase 1 |
| 420 | `clamp_soil_energy` | Q | **Phase 1** after a count |
| 421 | Test residues | C | **Phase 1** |
| 422 | `ed2_comparison.md` | D | Banner now; rewrite deferred |
| 423 | `post_proc` names | D | Stop-gap now; fix deferred |
| 424 | Palo Verde values | D | **Phase 1**, relabel only (D6) |
| 425 | θ_J cost | D | **Phase 1** |
