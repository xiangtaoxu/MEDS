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

[Unreleased]: https://github.com/xiangtaoxu/MEDS/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/xiangtaoxu/MEDS/releases/tag/v0.1.0
