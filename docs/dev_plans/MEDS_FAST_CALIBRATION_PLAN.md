# MEDS fast-parameter calibration plan

**Status:** written 2026-09-29 against `beta` at `522def6`, and revised the same day after the
owner's decisions (§12). P0d (#327) and P0h (#328) are merged. P0a, P0b, P0c and P0g are #329
(`feat/fast-calibration-p0`). P1–P3 are #330 (`feat/fast-calibration-tool`), stacked on it, and the
BCI fit is done. §13 records what the first fit showed and what to change next.

**Goal:** a crude, fast, repeatable refinement of the parameters that govern MEDS's sub-daily
physics, against eddy-covariance data, with the vegetation structure held at its initial state. It
also yields the parameters' covariance, as the starting point of a longer-term calibration. Barro
Colorado Island (BCI), started from the 2010 census, is the first site, and the calibration becomes
part of that example. The products are:
- a calibrated parameter set: a PFT file plus the fast blocks of the main TOML;
- its covariance;
- scores for the default and the calibrated sets on windows the fit never saw;
- a tool that repeats the procedure at another tower.

Speed is preferred over comprehensiveness at this stage. The plan fits short windows, and it uses
gradients, not space-filling samples (§6). It does not tune slow processes (allocation, growth,
mortality, recruitment, phenology, soil-carbon turnover).

## 1. Decisions

| # | Decision |
|---|---|
| D1 | **Frozen structure.** Every trial runs with `[run].slow_on = false`. `demography_on`, `do_cohort_fissfuse` and `do_patch_fissfuse` stay `true`, because they gate the census fusion at open; with them off, all 1,250 patches and 84,937 cohorts would reach the fast loop. Patches, cohorts, nplant, dbh, height, carbon pools, LAI and WAI are then constant (§3.1). |
| D2 | **Only fast parameters** (§4), for the one BCI PFT. The allometry, and therefore LAI (`[allometry] lai_b1/lai_b2`), is structure and stays fixed. |
| D3 | **Targets** are the tower's own fluxes with FLAG = 1 and its radiation components: reflected shortwave (albedo), upwelling longwave (surface temperature), net radiation, LE, H and GPP. NEE is a check, not a target, because frozen runs carry no growth respiration (§3.3). |
| D4 | **Closure: the Bowen-ratio-preserving correction** (owner, 2026-09-29). H and LE are scaled together so that H + LE closes on Rnet, with their ratio kept (§5.2). The evaporative fraction is scored uncorrected. |
| D5 | **Short windows** (owner). The fit uses a set of short windows, 10 days each, spread over the seasons of the calibration period (§5.1), not whole years. |
| D6 | **Trait plasticity stays on** (owner). The per-cohort Vcmax25, Rd25, SLA and leaf lifespan follow the canopy light gradient, and the calibrated values are the PFT's top-of-canopy values. |
| D7 | **Every trial restarts from a shared state** (owner). Each window has one state file, from a frozen spin-up with the current parameters. A restart re-acclimates the plastic traits to the trial's PFT file (P0b), so the traits follow the trial while the plastic canopy profile is kept. |
| D8 | **Gradient-based estimation** (§6, owner's question): Levenberg–Marquardt on the stacked, weighted residuals of every target in every window. The Jacobian comes from central finite differences run in parallel, and Gaussian priors come from each parameter's literature range. This gives the maximum a posteriori (MAP) set and its Laplace covariance. The first Jacobian doubles as the screening. |
| D9 | **Canopy interception both ways** (owner). The fit is run with `[fast].canopy_water_on` false and true, and the two are compared on the validation windows. |
| D10 | **The example shows both sets** (owner). The BCI README and figures show the default and the calibrated parameters side by side, and the calibration is a step of `run_example.py`. |
| D11 | **Fix first, then fit.** Code issues a calibration would absorb are fixed before fitting (P0, §7). The P0d physics fixes went in together as #327 (owner). |
| D12 | **Physical bounds.** Each parameter has a literature range, used both as the prior and as a hard bound, through a transform (§6.1). A target the model cannot reach inside the ranges is reported as a structural deficit, not forced. The albedo is the likely first case (§4.1). |

## 2. What the BCI tower offers (measured 2026-09-29)

From `BCI_v5.1.csv`, with the longwave columns read the right way round (`MEDS_FLUX_TOWER_FORCING_PLAN.md` D9):

| year | H, LE, NEE, GPP with FLAG = 1 | Rnet, reflected SW | both longwave components |
|---|---|---|---|
| 2012 (Jul–Dec) | 37 % | 0 | 0 |
| 2013 | 64 % | 48 % | 0 |
| 2014 | 65 % | 99 % | 0 |
| 2015 | 70 % | 88 % | 46 % |
| 2016 | 72 % | 88 % | 88 % |
| 2017 (Jan–Aug) | 80 % | 100 % | 100 % |

- **Closure.** On the 230 days with all 48 half hours good and Rnet present, (H + LE)/Rnet = 0.77. By day it is 0.72, and the slope of H + LE on Rnet is 0.73. The file has no ground heat flux.
- **Partition.** The daytime evaporative fraction is 0.61.
- **Night.** 17 % of night half hours have u* < 0.2 m s⁻¹; the median is 0.40.
- **Seasonal gaps.** FLAG = 1 coverage is lowest in February and March (50 % and 38 %), the dry season.
- **SWC** ranges 0.20–0.81. The top of that range exceeds any plausible porosity, so SWC is a relative index of wetting and drying only.
- **GPP** is the provider's partitioned product, with the partitioning's uncertainty.

## 3. What MEDS offers (read and measured 2026-09-29)

### 3.1 Frozen runs and shared states

`[run].slow_on` (`meds_config.f90:185-190`, read at `meds_config_io.f90:957`) gates the whole slow
step and the calendar-boundary restructuring (`meds_stepper.f90:86-124`); the fast loop still runs.
A full restart (`init_mode = 2`) restores the cohorts, patches, plastic traits and every fast
reservoir, and the run reports an exact restart of the canopy-air, soil and snow states. A state
is written at the end of any run with `[state].write_state = true`.

Measured at BCI on one core, all frozen:

| run | wall time | note |
|---|---|---|
| census start, one year | 127 s | 0.83 GB; stand identical at start, mid-year and end; budgets closed |
| census start, one day, writing a state (plasticity on) | 18.3 s | the census read, fusion and instant trait acclimation are most of it; a 138 kB state file |
| restart from that state, 10 days | 6.6 s | one trial of one window |

**Smoothness.** Four 10-day restarts with `stomatal_g1` at −1 %, 0, +1 % and +10 %:

| flux (10-day mean) | +1 % change / −1 % change | +10 % change / +1 % change | hourly correlation of the ±1 % changes |
|---|---|---|---|
| GPP | 1.006 | 9.53 | 0.999 |
| LE | 1.085 | 9.50 | 0.938 |
| H | 1.003 | 9.88 | 0.931 |
| NEE | 1.008 | 9.53 | 0.999 |

The model responds smoothly and nearly linearly at the percent level, so finite-difference
gradients are usable. LE's asymmetry is why the Jacobian is central (§6.2). Every other key's
smoothness is checked by the screening Jacobian.

**One restart mismatch.** The state run ends with LAI 5.6380 and the restart reports 5.6389. The
state file held dbh but none of the carbon pools, and a restart re-derived them on the allometry.
Cohort fusion keeps the pools and leaves a fused cohort below the allometry for its dbh (leaf area
is concave in AGB), so the restart's leaf area came out higher; with plasticity on it also reset the
fine-root carbon, by up to 2× in the most shaded cohorts. P0g stores the pools, the geometry and the
LAI above each cohort, and a run split at a restart now matches the unsplit run bit for bit.

### 3.2 How parameters reach a run

- **TOML only.** There are no command-line or environment overrides (`meds_main.f90:13-58`). A trial is a main TOML that points (`[init].pft_config`) at its own PFT TOML, with absolute paths. `scripts/numerics_sweep.py` already does this for numerics: it parses the TOML, sets keys, writes one config per cell and runs `meds_main` with a timeout.
- **Unknown keys are silently ignored** in the optional blocks, so a misspelt key keeps its default.
- **The provenance CSV is incomplete** (`meds_config_io.f90:1347-1358`), so it cannot prove what a trial ran with.
- **No in-process ensemble** (`MEDS_POLYGON_RUNTIME_PLAN.md` OR2). A trial is one `meds_main` process, and the fit runs many of them in parallel.
- **Restart traits.** A restart takes each cohort's `sla`, `vcmax25`, `rd25` and `llspan` from the state file (`meds_io.f90:459-470`), so a trial's PFT file would not reach them. With plasticity on, a census start sets each trait to its top-of-canopy value times exp(k × LAI above the cohort), instantly (`meds_driver.f90:184-187`; `meds_plant_trait_dynamics.f90`). In a frozen stand that profile is fixed, so repeating the instant acclimation at a restart reproduces exactly what a census start with the trial's PFT file would give (P0b).

### 3.3 What frozen runs lack

- **Growth respiration** reaches the fast NEE only as `slow_co2_rate`, which only the slow allocator sets. So `nee_fast` and `reco_fast` lack it.
- **Turnover water** to the ground (`shed_water_rate`) is also 0; it is negligible for the energy balance.
- **Phenology.** LAI stays at 5.6. BCI's canopy loses leaf area in the dry season, so dry-season windows carry that structural error.

### 3.4 Outputs

- **Hourly** (`meds_output_registry.f90:761-801`): `gpp_rate_fast`, `nee_fast`, `reco_fast`, `le_flux_fast`, `h_flux_fast`, `rnet_fast`, `ustar_fast`, `cas_temp_fast`, `soil_temp_top_fast`, `soil_water_site_fast`, and the forcing.
- **Daily only:** reflected shortwave (`sw_up_vis_site`, `sw_up_nir_site`), `lw_up_site` and the site leaf temperature. The hourly albedo and radiometric surface temperature need P0a.

### 3.5 Code issues a calibration would absorb

**Fixed in #327 (P0d):**
- the energy balance now uses the PFT's leaf and wood emissivities, as the radiation solver does;
- the wood boundary layer uses the wood temperature;
- the root profile comes from `[hydraulics].root_beta` and `root_depth`, and `[soil_column].root_beta` is refused;
- `[hydraulics].conductance = "segment"` makes `wood_kmax` and `vessel_curl` take effect;
- the stomatal closure and the drought-phenology cue use each PFT's own `leaf_pi0` and `leaf_elastic_mod`;
- the ground albedo and emissivity are settable in `[soil]`;
- the leaf-angle comment is corrected and the dead `phenology_on` key is removed.

On the five-year BCI census run these fixes move sensible heat from −39.01 to −39.08 W m⁻². Nothing
else changes at the README's precision.

**Fixed in #328 (P0h):** the reported sensible heat was about 100 W m⁻² too low, and every
surface-layer solve too stable. The reference air's potential temperature was referenced to the
ground, `T + (g/cp)·zref`, and the canopy air's was not, while the canopy-air energy budget exchanged
heat on the actual temperatures. So the reported H sat `g_ah·g·zref` below the flux the budget booked,
and the stability solve saw 0.37 K of stable stratification that was not there. Both potential
temperatures are now referenced to the canopy-air top. On BCI the mean H goes from −39.1 to 77.3 W m⁻²
(tower 32.4), its seasonal correlation from −0.59 to 0.95, and Rnet − H − LE over five years from
+104.6 to −11.4 W m⁻²; GPP, NEE, LE and Rnet move by less than 1%. H is a target (D3), so a fit before
this fix would have spent the aerodynamic keys on a reporting error.

**Still open:**

| issue | where | consequence |
|---|---|---|
| `leaf_absorptance` cancels when forcing is on | `meds_fast_dynamics.f90:1161`, `meds_leaf_gas_exchange.f90:61` | not identifiable: excluded |
| stem and root respiration always use the peaked form with fixed Ea, Hd and ΔS | `meds_plant_types.f90:276-291` | not calibratable; ignores `temp_response_form` |
| canopy interception is off by default (`[fast].canopy_water_on`) | `meds_config.f90:218` | no wet-canopy evaporation; fitted both ways (D9) |

### 3.6 The baseline the fit starts from (after #327 and #328)

The five-year census run, over the tower's measured hours. The corrected tower column applies §5.2.1's
daytime closure correction with one whole-period factor, f = 1.385 (the fit uses the sliding window):

| | MEDS | tower | tower, corrected |
|---|---|---|---|
| H, daytime (06–18 h) [W m⁻²] | 156 | 92 | 128 |
| LE, daytime [W m⁻²] | 107 | 145 | 201 |
| evaporative fraction, daytime | 0.41 | 0.61 | 0.61 |
| GPP, midday (10–14 h) [µmol m⁻² s⁻¹] | 30.0 | 21.7 | |
| H, night (19–05 h) [W m⁻²] | −0.1 | −23.6 | (not corrected) |
| u*, night / day [m s⁻¹] | 0.87 / 1.00 | 0.41 / 0.59 | |
| albedo | 0.26 | 0.13 | |

What this means for the fit:
- **Latent heat is the largest miss:** half the corrected tower's by day. GPP is too high at the same
  time, so the model's water-use efficiency is too high. The joint fit (§6) should raise
  `stomatal_g1` and lower `vcmax25` together. `stomatal_g1` may reach its 6.0 bound, which G5 would
  report.
- **Canopy interception matters more than it looked (D9).** The forcing carries 2.04 m of rain a
  year and the model evaporates 0.68 m. Interception loss in lowland tropical forest is commonly
  10–20% of rainfall: 0.6–1.1 mm d⁻¹ here, 16–32 W m⁻², a third to two-thirds of the deficit. The
  interception-on variant is expected to fit better.
- **Part of the deficit is energy.** The model absorbs about 25 W m⁻² less shortwave (§4.1), and the
  optics fit can recover only part of it (D12).
- **u\* is about twice the tower's**, and it is not a target (D3). The aerodynamic keys (§4.3) set both
  u* and H, so without u* they can trade against the stomatal keys. The tower measures u* at 41 m,
  near the 38 m mean canopy-air top, inside the roughness sublayer, where the model's Monin–Obukhov
  u* is not strictly comparable. A loose σ allows for that (§12, row 9).
- **Night H** is near zero against the tower's −24. It follows the night radiative cooling (Rnet −22
  against −33) and the canopy-air coupling, which `canopy_freeboard` and `ustmin` touch.

## 4. The parameters

Per-PFT keys are in the PFT TOML; the rest are in the main TOML. The defaults are the BCI values.
Each range is the prior's ±2σ and a hard bound (§6.1); P1's registry holds them with their
literature sources. The screening Jacobian (§6.3) decides which enter the fit.

### 4.1 Radiation

| key | default | range | constrained by |
|---|---|---|---|
| `leaf_reflect_nir` | 0.45 | 0.30–0.55 | albedo |
| `leaf_transmit_nir` | 0.25 | 0.15–0.45 | albedo |
| `leaf_reflect_vis` | 0.10 | 0.05–0.12 | albedo, absorbed PAR |
| `leaf_transmit_vis` | 0.05 | 0.02–0.08 | albedo, absorbed PAR |
| `leaf_clumping` | 0.80 | 0.50–1.00 | albedo, light profile |
| `leaf_angle_mean` | 45° | 45–65° | albedo, beam extinction |

The model's albedo is 0.26 and the tower's 0.13. A one-layer reconstruction of the model's
two-stream formulas (not a model run) gives:

| optics | albedo |
|---|---|
| shipped | 0.24 |
| NIR ρ 0.35, τ 0.15 | 0.18 |
| NIR ρ 0.30, τ 0.15; VIS ρ 0.08, τ 0.04 | 0.155 |

A leaf angle of 57.3° instead of 45° lowers any of these by 0.015. Reaching 0.13 on optics alone
needs a leaf NIR scattering (ρ + τ) near 0.45, well below measured leaves. The fit therefore stays
inside the ranges, and the albedo error left over is reported as structural: crown-scale clumping
and shading that a one-dimensional canopy does not see. The VIS optics also set absorbed PAR, which
is one reason the fit is joint (§6).

### 4.2 Photosynthesis and stomata

| key | default | range | constrained by |
|---|---|---|---|
| `vcmax25` (top of canopy) | 45 µmol m⁻² s⁻¹ | 25–80 | GPP, midday |
| `stomatal_g1` | 3.0 kPa^0.5 | 1.5–6.0 | LE and GPP together |
| `jmax_vcmax_ratio` | 1.7 | 1.4–2.2 | GPP light response |
| `stomatal_g0` | 0.01 mol m⁻² s⁻¹ | 0.001–0.05 | night LE |
| `theta_j` | 0.90 | 0.70–0.95 | GPP light response |
| `ds_vcmax` | 650 J mol⁻¹ K⁻¹ | 630–670 | GPP on the hottest afternoons |

### 4.3 Energy partition and aerodynamics

| key | default | range | constrained by |
|---|---|---|---|
| `[aerodynamics] z0m_ratio` | 0.13 | 0.05–0.20 | H, u* |
| `[aerodynamics] d_ratio` | 0.63 | 0.50–0.80 | H, u* |
| `leaf_width` | 0.04 m | 0.02–0.15 | H, surface temperature |
| `[aerodynamics] canopy_freeboard` | 5 m | 2–10 | canopy-air storage, night H |
| `[aerodynamics] ustmin` | 0.10 m s⁻¹ | 0.05–0.20 | night H |
| `[soil] dsl_dmax` | 0.015 m | 0.005–0.05 | LE after rain |
| `[soil] dewmx`, `intercept_k` (interception on only, D9) | 0.1, 0.5 | 0.05–0.3, 0.3–0.8 | LE on wet days |

### 4.4 Water stress

| key | default | range | constrained by |
|---|---|---|---|
| `wstress_sref_stomata` | 2.0 MPa⁻¹ | 0.5–5.0 | the dry windows' LE and GPP |
| `k_plant_max` | 6e-4 kg s⁻¹ MPa⁻¹ m⁻² | 2e-4–2e-3 | midday LE in the dry windows (`conductance = "whole_plant"`; the fit does not also take `wood_kmax`) |
| `wood_psi50` | −2.0 MPa | −3.5 to −1.0 | the same |
| `[hydraulics] root_beta` (log) | e⁻⁴ ≈ 0.018, with `root_depth` 2 m | 3.4e-4–0.37 (an e-folding of 0.25–2 m; ED2's 0.001 inside) | access to deep water |
| `leaf_pi0` (per PFT) | −1.5 MPa (ψ_tlp −1.71) | −2.5 to −1.0 | the dry windows' closure at 2 ψ_tlp |
| `[soil] psi_wilt` | −152.96 m | −300 to −100 | the same |

Ten days are too short for soil water to respond to these keys, so the dry windows constrain them
through the soil water the shared states carry (§5.1). That is the weakest part of a
short-window fit; the covariance (§6.4) shows it, and a longer-term calibration refines it.

### 4.5 Respiration

| key | default | range | constrained by |
|---|---|---|---|
| `rd_vcmax_ratio` | 0.015 | 0.008–0.025 | night NEE |
| `stem_resp_factor25` | 0.06 µmol m⁻² s⁻¹ | 0.02–0.15 | night NEE |
| `root_resp_factor25` | 0.30 µmol kgC⁻¹ s⁻¹ | 0.1–0.6 | night NEE |
| `[soil_carbon] resp_temp_increase` | 0.0757 K⁻¹ | 0.04–0.10 | night NEE against temperature |

Night NEE enters with FLAG = 1 and u* ≥ 0.2 m s⁻¹, with growth respiration added offline (§5.3).
The soil-carbon pools are frozen at their steady state, so heterotrophic respiration's size is the
spin-up's, not these keys'.

### 4.6 Not calibrated

- **Structure:** allometry, LAI, WAI, sapwood area, `crown_area_frac`.
- **Slow:** everything slow, including the plasticity slopes (not settable today) and `sla`, which changes leaf carbon but not LAI.
- **Not identifiable:** `leaf_absorptance` (§3.5). The leaf and wood emissivities: across 0.95–0.99 they move the upwelling longwave by roughly 0.04 (σT⁴ − L↓) ≈ 1 W m⁻², against the model's 10 W m⁻² deficit, which is canopy temperature.
- **Ground optics** (`[soil].ground_albedo_vis`, `ground_albedo_nir`, `ground_emissivity`, settable since #327): under LAI 5.6 they barely reach the tower. A sparser site could let the screening take them.
- **Rooting depth:** `root_depth` is held at the column depth, and `root_beta` alone sets the profile.
- **Biochemistry:** `kc25`, `ko25`, `gstar25` and their activation energies (Bernacchi 2001).
- **Soil properties:** `[soil_column]` retention and `ksat`, which are site measurements.
- **Numerics:** `dt_fast` and the tolerances, fixed at the run's values, because hydraulic parameters absorb time-step error (`meds_config.f90:639-687`).

## 5. Windows, targets and residuals

### 5.1 Windows and states

**Calibration windows.** Eight windows of 10 days in 2015-08 → 2017-07, the period with every
radiation component observed:
- two in the wet season;
- two in the wet-to-dry transition;
- two in the dry season, one of them in the 2016 El Niño drought;
- two in the dry-to-wet transition.

Each window is the 10-day span in its season with the best FLAG = 1 coverage and observed
longwave. **Validation windows** sit at the same calendar places in 2013–2014, and in 2015–2017
between the calibration windows.

**States.** A frozen spin-up chain starts from the census at 2015-04-01, four months before the
first window. It runs to each window's start in turn, writing a state there and restarting from it
to the next. That is one frozen pass of about two years, about 5 minutes. The soil water and
temperature at each window's start are then those of a spun-up run, not the census's
initialization. The chain runs with the current parameters, so it is re-run once with the MAP set
and the fit refined from there (§6.2).

### 5.2 Targets

For each window, the hours used are the tower's FLAG = 1 hours for the turbulent fluxes, and hours
whose forcing is observed (qc 0 for `LWdown` and `Wind`). The targets and their observation
errors σ:

| target | model | tower | σ |
|---|---|---|---|
| albedo (hours with SW↓ > 200 W m⁻²) | `sw_up_fast` / `sw_in_fast` (P0a) | Rs_dn / Rs | 0.01 |
| upwelling longwave | `lw_up_fast` (P0a) | the file's `Rl_dn` column | 5 W m⁻² |
| net radiation | `rnet_fast` | Rnet | 10 + 0.05·\|Rnet\| W m⁻² |
| LE, H | `le_flux_fast`, `h_flux_fast` | closure-corrected (§5.2.1) | 10 + 0.15·\|flux\| W m⁻² |
| evaporative fraction (daytime means) | LE/(H + LE) | uncorrected | 0.05 |
| GPP (daytime) | `gpp_rate_fast` | gpp | 1.5 + 0.15·GPP µmol m⁻² s⁻¹ |
| night NEE (u* ≥ 0.2) | `nee_fast` + growth respiration | NEE | 2 µmol m⁻² s⁻¹ |
| u* | `ustar_fast` | ustar | 0.1 + 0.2·u* m s⁻¹ |

The σ forms follow the random-error scaling of eddy-covariance fluxes (Richardson et al. 2006). The
coefficients are starting values; the fit reports each target's normalized residual so that σ can
be revisited.

#### 5.2.1 Closure

Following FLUXNET2015's energy-balance-closure correction, in a sliding ±15-day window of daytime
hours:

```
f = Σ Rnet / Σ (H + LE),   H_c = f · H,   LE_c = f · LE
```

G is taken as 0, since the file has none. Night hours are not corrected.

### 5.3 Growth respiration for night NEE

The five-year census run with the slow loop on gives the stand's growth respiration
(`growth_resp_site`, daily) by calendar month. Before comparison, the frozen model's NEE is
increased by that monthly climatology.

### 5.4 The residual vector

For parameters u in the transformed space (§6.1), the residual vector is

```
r(u) = [ (y_model(u) − y_obs) / σ ]   for every target, window and hour used
       [ (u − u_prior) / σ_prior ]    one row per parameter
```

and the fit minimizes ½‖r(u)‖². Hourly residuals are autocorrelated, so the covariance (§6.4)
inflates each target's weight by its effective sample size, from the lag-1 autocorrelation of its
residuals at the MAP.

## 6. The method

This borrows two ideas from differentiable modeling: gradients of the model output with respect to
its parameters, and the curvature they give. The gradients here come from finite differences, not
automatic differentiation (§11).

### 6.1 Transform and prior

- **Transform.** A parameter bounded to [a, b] is fitted as u = logit((θ − a)/(b − a)). A positive scale parameter (conductances, `k_plant_max`, the respiration factors) uses log first, then the same bound. The fit is unconstrained in u, and every θ stays inside its range.
- **Prior.** Gaussian in u, centred on the default, with σ_prior such that the range is ±2σ.

### 6.2 Levenberg–Marquardt

Each iteration:
1. **Jacobian.** J = ∂r/∂u by central differences, one perturbation of ±h per parameter, with h = 0.04 in u (about 1 % of the range near the centre). For k parameters and W windows that is 2kW trials of 10 days, all independent. The base run for W windows comes with them.
2. **Step.** Solve (JᵀJ + λ·diag(JᵀJ)) δ = −Jᵀr. Run the W windows at u + δ; accept the step and lower λ if ‖r‖ fell, otherwise raise λ and try again.
3. **Stop** when the relative drop in ‖r‖² is below 10⁻³, or after 15 iterations.

**Robustness:**
- **Multiple starts.** The fit starts from the default and from two points drawn from the prior; the best MAP is kept, and disagreement between the three is reported.
- **Failed trials.** A failed or timed-out trial (twice the median time; the hydraulic cost cliff, ROADMAP #104) makes its step a rejection, and a failed Jacobian column is retried at h/2.
- **Rough keys.** A key whose central difference is not smooth (§6.3) is left out of the gradient fit and set by a 1-D line search at the end.

**State refresh.** After the first MAP, the state chain (§5.1) is re-run with the MAP set, and the
fit continues from the MAP for at most five more iterations. This removes most of the mismatch
between the default-parameter soil water and the trial's parameters.

### 6.3 Screening from the first Jacobian

The Jacobian at the default serves as the screening. It gives:
- **the sensitivity of every target to every key**: normalized, J scaled by σ_prior;
- **identifiability**: the singular values of the weighted Jacobian. A key that loads only on directions with small singular values is fixed at its prior mean;
- **collinearity**: pairs of keys the tower cannot separate, which then share one direction;
- **harness checks**: a key with an exact-zero column is a bug (misspelt, dead, or a restarted trait), not an insensitive parameter;
- **smoothness**: the ratio of the + and − differences, as in §3.1. A ratio far from 1 marks a rough key (§6.2).

### 6.4 Covariance

At the MAP the Laplace approximation gives

```
Σ_post ≈ (JᵀWJ + Σ_prior⁻¹)⁻¹
```

with W the effective-sample-size weights of §5.4. The deliverables are the MAP, Σ_post in θ and u,
the correlation matrix, and each key's posterior-to-prior σ ratio: how much the tower taught.

A linearity check comes with it: for each of the three leading principal directions, runs at ±1σ
compare the actual change in ‖r‖² with the quadratic prediction. Where they differ by more than a
factor of 2, the covariance is reported as local only. The covariance is the prior or the proposal
of a longer-term calibration (§11).

### 6.5 Validation

With the MAP set:
- the validation windows (§5.1), scored as in §5.2;
- the full five-year census run with the slow loop on, beside the default run, in the example's own figures (D10).

The calibrated set is accepted only if it passes the gates of §8.

## 7. Phases

**P0 — model changes the fit needs**, each with a test:
- **P0a (done: `sw_up_fast`, `lw_up_fast`):** hourly `sw_up_fast` (VIS + NIR) and `lw_up_fast`, staged beside `rnet_fast` (`meds_fast_dynamics.f90:765-800`).
- **P0b (done: `reacclimate_plant_traits`):** `[init].reacclimate_traits`, default false. On a restart with plasticity on, it repeats the census start's instant acclimation (`meds_driver.f90:184-187`) with the current PFT file. With plasticity off, it keeps the PFT file's traits instead of the state's. Test: a census start and a restart with a changed `vcmax25` give identical per-cohort traits in a frozen stand.
- **P0c (done: the parameter record `<prefix>_parameters.csv`, every key the loader read):** complete the parameter provenance, so every §4 key is in the run's provenance file. The tool checks each trial's provenance against what it set, which also catches a misspelt key.
- **P0d — #327, open against `beta`:** the physics fixes of §3.5, in one pull request (owner), with a before-and-after on the BCI example.
- **P0h — #328, stacked on #327:** the reported sensible heat and the stability solve on the budget's temperature basis (§3.5, §3.6).
- **P0g (done: the state stores the pools, geometry, LAI above, film water, growth buffer and `slow_co2_rate`):** a restart reproduces the state it was written from. Explain or remove the LAI mismatch of §3.1, with a test that a run split at a restart matches the unsplit run.

**P1 (done) — `scripts/calibrate_fast/`**, modelled on `scripts/numerics_sweep.py` and kept small:
- **registry:** a parameter registry (TOML) giving each key's file, section, default, range, transform and source.
- **site declaration:** windows, targets, σ and the tower file. It lives in the example (`examples/example_flux_tower_bci/calibration.toml`).
- **states:** the state chain of §5.1.
- **trials:** a trial writer that deep-sets keys in parsed copies of the base TOMLs, and a runner (`meds_main`, one thread, timeout, provenance check).
- **residuals:** the reader and the residual builder of §5.
- **fit:** Levenberg–Marquardt with the parallel Jacobian, the screening report, the covariance and the linearity check.
- **workers:** one Slurm allocation holds a pool of single-thread workers for the whole fit, and the driver dispatches trials to them, so no trial waits in the queue. On a workstation the same pool is local processes.

**P2 (done; §13) — BCI** (D9, D10):
- the fit with interception off and on;
- the validation;
- `pft_parameters_calibrated.toml` and the calibrated main-TOML blocks in the example;
- `run_example.py --calibrate` to repeat the fit, and the default and calibrated runs side by side in the README and figures.

**P3 (done: 17 unit tests and ctest `calibrate_fast`) — tests:**
- pytest for the registry and transforms, the trial writer (keys land where the loader reads them), the closure correction, the residuals, and the fit on a synthetic linear model with a known answer and covariance;
- a CTest smoke test: one Jacobian column on a 3-day window that moves the output, and a repeated trial that reproduces it byte for byte.

## 8. Gates

| gate | requirement |
|---|---|
| G1 | the stand is identical at the start and end of every trial (true today, §3.1) |
| G2 | the same parameters give byte-identical output (one thread; the OpenMP build is thread-invariant) |
| G3 | every fitted key has a non-zero, smooth Jacobian column (§6.3) |
| G4 | on the validation windows the MAP lowers ‖r‖², and no target's normalized RMSE rises by more than 10 % |
| G5 | every MAP value lies inside its range; one within 5 % of a bound is reported with the target that pushed it there |
| G6 | the three starts reach MAPs whose ‖r‖² agree within 5 %, or the disagreement is reported |
| G7 | the full five-year run with the slow loop on completes with the calibrated set and closed budgets |

## 9. Cost

A 10-day trial takes 6.6 s on one core (§3.1). For k = 20 keys and W = 8 windows:

| step | trials | core-hours | wall time on 160 cores |
|---|---|---|---|
| state chain (per refresh) | 1 chain | 0.08 | ~5 min, sequential |
| one Jacobian (2kW + W) | 328 | 0.6 | ~20 s |
| one LM iteration with 1–3 step tries | ~340–350 | 0.65 | ~30 s |
| a fit: 3 starts × 15 iterations + refresh | ~17,000 | ~30 | ~30 min |
| both interception variants | ×2 | ~60 | ~1 h |
| validation | ~100 trials + 2 full runs | ~1 | ~15 min |

The worker pool keeps queue time out of these numbers. On a 4-core workstation one Jacobian takes
about 9 minutes and one start's fit a few hours.

## 10. What limits it

- **Short windows.** The soil-water memory and the slow hydraulic response are carried by the shared states, not fitted (§4.4). The covariance says how poorly the water-stress keys are known, and the longer-term calibration refines them.
- **Local covariance.** The Laplace covariance is valid near the MAP; the linearity check (§6.4) marks where it is not.
- **Frozen LAI.** BCI's dry-season leaf loss is absent (§3.3).
- **One PFT.** A single PFT stands for a diverse canopy, so the fitted values are effective canopy values.
- **Footprint.** The tower sees the island's forest; the stand is the 50-ha plot's census.
- **The tower's own uncertainty:** closure, the GPP partitioning, and night-time flux under low turbulence.
- **Structural error.** The albedo (§4.1) and the friction velocity at twice the tower's (§3.6) may be structural, so the gates and bounds are there to show it.

## 11. Toward differentiable calibration

The owner asked whether ideas from differentiable modeling apply. What this plan takes, and what
would come later:

1. **Now: gradients by parallel finite differences.** MEDS is Fortran 2018 with OpenMP and an adaptive integrator, and cannot be differentiated automatically today. Its response is smooth at the percent level (§3.1), so a central-difference Jacobian costs 2k short trials, all in parallel. With it come the gradient-based fit, the identifiability analysis and the curvature, that is, the covariance. The same machinery would take an automatic-differentiation Jacobian unchanged.
2. **Next: a differentiable surrogate.** The fit's ~17,000 trials sample the parameter space densely near the MAP. They can train an emulator of each window's outputs (a neural network or a Gaussian process) whose gradients are exact and cheap. An emulator allows MCMC for a proper posterior, and fits over many more windows, as the longer-term calibration.
3. **Later: a differentiable fast loop.** Differentiating the fast loop itself, by source-to-source or LLVM-level AD (Tapenade, Enzyme) or by a port of the column physics to an AD framework, is a project of its own. It would allow differentiable parameter learning across towers: a network that maps site attributes to parameters, trained through the model. This plan's registry, residuals and gates carry over to it.

## 12. Owner decisions

| # | question | decision (2026-09-29) |
|---|---|---|
| 1 | closure handling | the Bowen-ratio-preserving correction (D4) |
| 2 | periods | short windows, speed first (D5): 8 × 10 days in 2015-08 → 2017-07, validated elsewhere (§5.1) |
| 3 | canopy interception | fit both ways (D9) |
| 4 | shared state | yes, with trait plasticity kept (D6, D7, P0b) |
| 5 | order of the P0d fixes | all five in one pull request: #327 |
| 6 | the example | shows the default and calibrated sets; the calibration is part of it (D10) |
| 7 | where the root profile lives | `[hydraulics]`, as a plant trait, not `[soil_column]` (#327) |
| 8 | the sensible-heat bias | a reporting error, fixed before the fit (#328) |
| 9 | u* as a target | yes, with σ = 0.1 + 0.2·u* m s⁻¹ for the roughness-sublayer caveat (§3.6, §5.2) |
| 10 | order of the remaining P0 work | all of P0a, P0b, P0c and P0g in one pull request, then the tool and the fit |
| — | method | gradient-based, with covariance (D8, §6, §11) |

## 13. What the first fit showed (2026-09-29)

The BCI fit ran on the calib-tool build, with trait plasticity on in the example
(`examples/example_flux_tower_bci/README.md`, "Calibrating the fast parameters", has the tables).

### 13.1 Results

| | interception off (shipped) | interception on |
|---|---|---|
| keys screened / fitted / rough | 28 / 20 / 2 (`wood_psi50`, `leaf_pi0`) | 30 / 20 / 0 |
| trials, failed | 19,096, 2 | 14,800, 2 |
| wall time, cores | 37 min on 320 (8 × 40) | 56 min on 256 |
| median trial, loaded node | 18.6 s (8.5 s on an idle node) | 47 s |
| objective, calibration windows | 70,094 → 25,971 | 69,732 → 25,925 |
| objective, validation windows | 68,018 → 29,790 | 66,985 → 30,093 |
| start spread (G6) | 0.3 % | 0.3 % |
| keys at a bound (G5) | 8 | 8 |
| G7, five years with the slow loop | pass | fails: 53 water-budget breaches |

- **G1–G6 pass for both variants.** G1 first failed on `leaf_carbon`: re-acclimation rebuilt leaf
  carbon from leaf area and SLA. It now scales the carbon by the SLA's change (#329).
- **Every validation target improves**, by 7 % (upwelling longwave) to 70 % (albedo).
- **Over five years the calibrated set fixes the level of GPP (6.75 against the tower's 7.46),
  net radiation and u\*.** LE lands on the tower's closure-corrected value.
- **Structure shows at the bounds.** The eight at-bound keys are `vcmax25`, `jmax_vcmax_ratio`,
  `rd_vcmax_ratio`, `stomatal_g1`, `leaf_clumping`, `canopy_freeboard`, `z0m_ratio` and
  `leaf_transmit_vis`. H stays 36 W m⁻² high with `g1` at its ceiling, so the gap is the
  energy-partition structure (§3.6), not a parameter.

### 13.2 What went wrong, and the rule each gave

- **Rough keys: default, not line search.** The line search set `wood_psi50` = −1.24 and
  `leaf_pi0` = −2.43. Those values fit the windows better but broke the five-year water budget:
  659 breaches. Bisection found that reverting the two keys closes it.
  - Rule: `[fit].rough_keys = "default"`.
  - The cost: the late dry season is too stressed (13.3).
- **A trial fails on any whole-site budget breach**, not only on a crash. A breach means the
  candidate is off the model's valid domain.
- **The timeout floor is the configured timeout, and the median counts only completed trials.**
  The first interception-on fit timed out every trial. Its timeout was set from the first trials,
  which ran on an idle node 2–3× faster than a full one.
- **HDF5 is not thread-safe.** The driver reads trial output from threads, which segfaulted until
  every netCDF read took one lock.
- **The trial and chain digests must be canonical and must include the base configs.** The chain
  key lacked them, so a changed base config reused stale states.
- **Validation windows need their own forcing mask.** The 2013–14 windows had no hour with observed
  longwave, so the validation mask requires observed wind only.
- **Trait plasticity exposed an undeclared tissue heat change** in the slow ledger's allocate
  phase: −5,903 J m⁻² over five years. It is now declared.
- **The per-layer soil-water face check** (`faces[soil_layer_mass]`) worsens under every calibrated
  set. The worst residual is 1.1–2.6 kg m⁻², against 0.0011 at the default, while the whole-site
  water budget closes. With `debug_error`, `advance_soil_water_column` stops in October 2012 on
  "per-face mass budget did not close". This is a solver issue to fix in the model, not the
  calibration's.

### 13.3 The late dry season

In the five-year run with the shipped set, April GPP collapses in 2014, 2016 and 2017. In April
2016 it is 2.9 against the tower's 6.3 µmol m⁻² s⁻¹. The fit's 2016-04-13 window shows the same:
daytime GPP 6.4 at the MAP against the tower's 14.3.

- **Every Jacobian neighbour of the MAP gives 6.3–6.4.** Only the rough keys at their line-searched
  values give 9.5.
- **The chain states are not the cause.** The MAP chain's soil water at the window start matches
  the five-year run's: 0.054 against 0.053 m³ m⁻³ averaged over the rooted layers.
- **The dry-season response is carried by the two keys the method cannot use.** A fixed-state
  Jacobian also never sees that more transpiration in January dries April.

So the next fit needs two things:
- **A smooth hydraulic response,** or a derivative-free step for the hydraulic keys under a budget
  constraint.
- **Windows long enough, or chained, for the dry-season drawdown to enter the gradient.**

### 13.4 Efficiency, measured

- **Utilization was about 50 %.** The three starts converged to the same point, so two of them were
  wasted work.
- **The serial chains took about 32 % of the wall time.** Each chain is one long run, and each
  refresh reruns it.
- **A trial on a full node runs 2.2× slower** than on an idle one, for a cause not yet measured.
  Fewer slots per node may give more trials per hour.
- **The LM hit its iteration cap (15)** while still improving by 0.02–1 % per iteration.
- **Eight of the 20 keys sat at a bound** yet were differenced every iteration: 16 of the 41 trials
  per window per Jacobian.
