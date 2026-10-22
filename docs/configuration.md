# Configuring a MEDS run

A run is described by **two [TOML](https://toml.io) files**, and the first is the only command-line
argument:

```bash
./build/meds_main path/to/run_main.toml
```

With no argument, `meds_main` reads `./meds_config_main.toml`. Relative paths inside the config
resolve against the directory you launch from, not the config's own location.

## The parameter philosophy

**MEDS hard-codes no model parameters.** The source defines only true *constants* — numerical,
geometric and calendar. Every parameter is user-mutable in TOML and **required**: a missing
parameter, or a missing file, is a **hard error**. The reader builds a presence map while loading
and aborts listing every absent key, rather than falling back to a default nobody chose.

This is deliberate and it is not free — a MEDS config is long. The reason is that a silent default
is indistinguishable from a considered value once a run is a year old, and the model has already
been bitten by exactly that: a flag whose loader default disagreed with its type default ran an
entire subsystem off for every config that omitted the key.

Two exceptions, both explicit:

- **Derived quantities** are computed from the primary ones — the step in years, the height class
  edges, the leaf capacities from their ratios to Vcmax25. No PFT trait is derived from another: a
  trait that follows a fitted relation (the sapwood conductivity or the mortality hazard from wood
  density, the plasticity slopes from Vcmax25) is computed outside MEDS, by
  [`scripts/pft_traits`](../scripts/pft_traits/README.md), and written into the PFT file.
- **Optional keys.** Solver settings, process switches and a few late-added blocks (`[soil]`,
  `[snow]`, `[soil_carbon]`) read each key with a fallback to its default, so a feature needs no
  TOML edits. The reference file lists each one commented out at that default.

**Every key MEDS reads is listed in the two reference configs**, `meds_config_main.toml` and
`meds_config_pft.toml`, set or commented out at its default. They are the complete list. A key that
is not listed in its file's reference stops the run before it starts: the error names every such key
at once, together with any missing required ones. Each line gives a suggestion:
- a key that moved section gets its new home (`state.cohort_max` → `output.cohort_max`);
- a misspelling gets the key it most likely meant.

Renamed and retired keys are listed in [`CHANGELOG.md`](../CHANGELOG.md). A listed key that the run's
other settings give nothing to do is refused the same way, with the setting that makes it inapplicable:
`forcing.path` with `format = "ERA5land"`, a `[region]` block in a site run, a soil key the
`[maps]` provide.

A key that parses and does nothing is worse than one that is absent: `[forcing] dt_forcing = 7200`
(the key is `timestep`) would otherwise run silently at the file's own spacing. The test
`config_keys_listed` holds the two references equal to the keys the loader reads, so a new key is
listed when it is added.

Every run that writes a state writes `<output_dir>/<prefix>_pft_parameters.csv` — one row per PFT
with every per-PFT trait actually used. It is a provenance record: what the run really did, not what
you meant.

Every run also writes the **parameter record**, `<prefix>_parameters.csv`, beside its diagnostic
output (`[output].dir`, `[output].prefix`) and beside its state files (`[output].state_file_dir`,
`state_file_prefix`), whichever it writes. It
has one row for every key the loader read, from every file it read (the main file, the PFT file,
the output list): `source,key,index,present,value`. `index` is 0 for a scalar and the element for an
array; `present` says whether the key was in the file or took its compiled default; `value` is what
the run used, to the last digit.

## The two files

### The main file

All non-PFT settings. Named on the command line; it names the PFT file via `[init].pft_config`.

| Block | What it sets |
|---|---|
| `[run]` | The slow timestep, and the run span as **calendar dates** (`start_time`, `end_time`, leap-year-aware Gregorian, `"YYYY-MM-DD[ HH:MM:SS]"`). Thread count. The run `mode`: one site, or a region. |
| `[fast]` | The sub-daily loop: `dt_fast`, the integrator, tolerances, error control; the relative tolerance `rtol` of every error-controlled store, the soil-plant water network's absolute tolerance `water_atol`, the budget checks (`soil_energy_budget_atol`, `debug_error`) and the canopy interception (`intercept_alpha`, `intercept_k`). |
| `[init]` | How the run starts, and the path to the PFT file. |
| `[demography]` | Cohort and patch fusion/fission, the cadence switches. |
| `[disturbance]`, `[recruitment]` | Treefall rates; the seed-rain and reproduction settings. |
| `[soil]` | The **ground**: depth and layers, hydraulic texture, retention family and thermal properties; its surface (the pond, the dry surface layer, and the optics `ground_albedo_vis`, `ground_albedo_nir` and `ground_emissivity`); and its two bottom boundaries, `water_bottom_bc` and `thermal_bottom_bc` with the thermal anchor's `deep_temp` and `deep_depth`. |
| `[snow]`, `[aerodynamics]` | The snow store's and the canopy aerodynamics' parameters. |
| `[soil_carbon]` | The CENTURY decomposition: selectors, rate parameters, cold-start spin-up. |
| `[forcing]`, `[site]` | The meteorological driver, and where the site is. |
| `[[sites]]` | A site run's polygon as several sites, one table each: its share, terrain and ground ([Several sites](#several-sites)). |
| `[region]` | For `mode = "region"` only: the box of forcing cells and which of them to simulate. |
| `[output]` | What a run writes: the state files (`state_*`) and the diagnostics (`io_config`, the tiers). |

**`[soil]` is the ground; `[fast]` holds the solvers' settings.** Every property of the ground is in
one block, whichever solver reads it: the column, its surface and both bottom boundaries. The column
is validated at load, because layer bottoms that do not increase, or a saturated water content below
the residual, produce a silently wrong column rather than a crash.

**The soil layers.**
- **Default:** every column cuts one list of layer bottoms at its `depth`: SoilGrids 2.0's depths
  (0.05, 0.15, 0.30, 0.60 and 1.00 m), then 1.5, 2, 3, 4.5, 6.5, 9, 12, 16, 21 and 27 m. A remainder
  thinner than half of the layer it falls in joins the layer above.
- **What that gives:** the top layer is 5 cm everywhere, shallow soils get few layers, and every column
  shares one depth axis.
- **Override:** `layer_bottoms = [...]` [m] replaces the list for a measured profile, and the column
  ends at its last bottom unless `depth` cuts it.
- **The science:** [`science/soil_biophysics.md`](science/soil_biophysics.md), "The soil layers".

**The bottom water boundary — `[soil].water_bottom_bc`.** `free_drain` (default), `bedrock` or `aquifer`.
The first two are the ends of one drainage boundary, $`q = f\,K(\theta_n)`$ at $f = 1$ and $f = 0$; a
leaky base, $0 < f < 1$, comes from the soil maps' `f_leak` (Input maps, below). `aquifer` is a water
table at the column base, for wet sites only ([`science/soil_biophysics.md`](science/soil_biophysics.md),
"Boundary conditions").

### The bottom thermal boundary — `[soil].thermal_bottom_bc`

The default 2.0 m column sits against an annual thermal damping depth of about 2 m, so the annual
temperature wave has barely attenuated by the time it reaches the base. `thermal_bottom_bc` decides what
happens to it there.

| key | values | meaning |
|---|---|---|
| `thermal_bottom_bc` | `geothermal` (default) \| `dirichlet` | prescribed bottom flux (held at zero — an adiabatic wall that **reflects** the annual wave) vs conduction to a fixed deep temperature |
| `deep_temp` | [K] | the anchor temperature. **Required** when `thermal_bottom_bc = "dirichlet"`; there is no default |
| `deep_depth` | [m], default 3.145 | depth of the anchor plane below the surface; must lie below the bottom soil node |

`deep_temp` is the mean annual soil temperature below the damping depth — close to the mean annual
air temperature of the forcing that drives the run, and a site property like latitude. It is required
rather than defaulted because an error in it is a steady flux into the column base, and so a
mean-annual bias in deep-soil temperature; a silent default would put back, somewhere else, the bias
this boundary condition exists to remove.

`deep_depth` is a physical choice, not a numerical one, and the default is derived: a resistive
termination reflects least when the anchor sits `d/sqrt(2)` below the bottom node, `d` being the
annual damping depth. For the default column that is 1.75 + 1.973/√2 = 3.145 m. **Recompute it if you
change `[soil].depth` or the thermal texture** — the derivation, the measured sweep and what
the anchor does and does not fix are in [`docs/science/soil_biophysics.md`](science/soil_biophysics.md).

Deepening the column with `[soil].depth` remains available and is the other way out, at the
cost of layers; it also moves the root profile and the drainage, so it is not a thermal-only change.

### Threads and the slow tier — `[run]`

- **Threads.** `[run].n_threads` (default 1) runs a site's patches, or a region's polygons, on that
  many OpenMP threads. Patches and polygons are independent within a step and the site's sums are
  folded in patch order, so the output is byte-identical at any thread count: it is a wall-clock
  setting only, bounded by the patch count and by how unevenly the patches cost. The default of 1
  keeps a build compiled with OpenMP serial until asked; a build without OpenMP ignores the key. MEDS
  does not read `OMP_NUM_THREADS`, so a run's config records its thread count
  ([`science/numerical_scheme.md`](science/numerical_scheme.md) §6a).
- **`[run].slow_on = false`** skips the whole slow tier every step: growth, mortality, phenology,
  traits, demography and soil carbon. The stand stays as it is while the fast loop runs. It is broader
  than `[demography].demography_on`, which freezes only the restructuring.
- **The slow ledger** (`[run].slow_ledger_on`, default on) snapshots the site's carbon, water and
  energy stores either side of every slow-step operator, and prints at the end of the run what each
  operator could not account for, judged against the flux it declared
  ([`science/order_of_processes.md`](science/order_of_processes.md) §4). It writes no state and changes
  no result, and costs about 0.3 % of the wall clock. It is the one check that sees a store jump across
  the slow step; the fast loop's ledgers sum per-sub-step residuals and cannot.

### Cohort and patch fusion — `[demography]`

Fusion merges similar cohorts within a patch every month, and similar patches every year. Both run
in passes. Each pass merges every pair that passes a similarity test at one tolerance, and the
tolerance grows geometrically from pass to pass, from a minimum to a **ceiling**. Passes stop early
once the count is within its **target**.

| | Cohorts | Patches |
|---|---|---|
| similar when | their heights differ by less than `tol · h_top`, with `h_top` the height of the patch's tallest cohort, and their combined LAI is below `cohort_lai_cap` | their light profiles differ by at most `tol` on average, and by at most `tol · patch_light_maxdev_factor` in any layer |
| tolerance, first pass | `cohort_size_tol_min` | `patch_light_tol` |
| tolerance ceiling | `cohort_size_tol_max` | `patch_light_tol_max` (default 0.15) |
| passes, at most | `n_cohort_fusion_iter` | `n_patch_fusion_iter` |
| target | `max_cohort`, in the most crowded patch | `max_patch` |

The cohort tolerance scales with the stand: in a 3 m regrowth patch the first pass merges cohorts
within about 6 cm, so a cohort that has pulled ahead keeps its light and its suppressed neighbours can
thin, while in a 30 m forest it merges within about 0.6 m.

A patch's light profile is the fraction of full sunlight left under each of `n_height_layers` equal
layers up to the tallest `hgt_max`, `exp(−k × LAI at and above the layer)` with k =
`light_ext_coeff_fusefiss` (default 0.5), which nothing else uses. Layers fully lit in both patches
are not compared.

**The ceiling is hard and the target is not.** Two cohorts or patches more different than the
ceiling stay apart even if the count is then above its target, because merging them would lose
structure the model needs. A run in which that happened says so at the end, for example
`NOTE: patch fusion left 24 patches on 2013-01-01 (max_patch = 12): ...`. A negative `max_cohort` or
`max_patch` forces fusion instead, ignoring the test.

`patch_light_tol_max` is optional. Absent, it is 0.15, or `patch_light_tol` if that is larger. Set
it to `patch_light_tol × 1.5^(n_patch_fusion_iter − 1)`, 0.759375 for the shipped 0.10 and six
passes, to reproduce the schedule before the ceiling existed.

### The PFT file

Everything PFT-specific: the `[pft]` trait table, each PFT's size allometry, `[allometry]`, and each
PFT's plant hydraulics, `[hydraulics]`. **The number of PFTs is the length of `[pft].wood_density`** —
every other trait array must match it. Every `[allometry]` and `[hydraulics]` key takes one value, for
every PFT, or an array with one value per PFT.

`[hydraulics]` holds plant water transport (`science/plant_hydraulics.md`): the leaf and wood
pressure–volume curves, the xylem's vulnerability (`wood_psi50`, `wood_kexp`), its sapwood specific
conductivity `wood_kmax` and path factor `vessel_curl`, and the roots (`root_beta`,
`root_depth_b1`/`b2`, `specific_root_length`, `fine_root_radius`). The root profile is a plant trait,
so these keys live in the PFT file only.

Wood density is the primary PFT axis. The mortality hazard's parameters, `[pft].mort_gamma`,
`mort_alpha` and `mort_beta`, ship as Camac et al. (2018)'s power laws in it, so low-density PFTs get a
higher growth-independent hazard and a steeper low-growth penalty; `scripts/pft_traits` computes them,
and a site's own census fit replaces them (example 05 does). The one exception to the growth-driven
hazard is a cohort with no leaves and no storage, which dies in the step it is found
(`science/plant_carbon_allocation.md` §2). Growth, competition and reproduction parameters are per-PFT
but ship uniform.

**Recruitment.** A PFT's recruits are a baseline external seed rain, `[pft].seed_rain_recruits`
[plants m⁻² yr⁻¹], plus those its reproduction buys. From `min_repro_height` up, a cohort diverts
`reproduction_investment_fraction` of its growth carbon to reproduction
([plant_carbon_allocation](science/plant_carbon_allocation.md) §2). `recruit_carbon_efficiency` of
that carbon (seed set, germination and seedling survival together) becomes recruits at
`min_cohort_height`, a share that falls as exp(−`recruit_shade_decay` × patch LAI). `include_pft = 0`
keeps a PFT from recruiting.

Each PFT's height follows the curve `[allometry].height_allometry` names for it (one name for every
PFT, or one per PFT, `["chave2014", "power", "power"]`). `"power"`, the default, is the pan-tropical
power law `exp(b1Ht + b2Ht ln dbh)`. `"chave2014"` is the log-quadratic ln H = `chave_a` − E +
`chave_b` ln dbh + `chave_c` (ln dbh)² + `chave_rse`²/2, taken at its mean as BIOMASS's `retrieveH` does;
its coefficients default to Chave et al. (2014)'s pantropical fit (0.893, 0.760, −0.0340, 0.243), and a
site's own heights give a local fit (examples 03 and 05 use one for Barro Colorado Island); the
`chave_*` keys are refused when no PFT uses `"chave2014"`. E
(`[allometry].chave_e`, one value) is the site's stress index from its climate (Chave's eq. 6b), the same
for every PFT: the drier, more seasonal or hotter the site, the larger E and the shorter its trees. With
`[maps].chave_e_map`, each polygon takes E from its own cell instead, and `chave_e` is refused
([Input maps](#input-maps)). Every curve is capped at the PFT's `[allometry].hgt_max`, and biomass
(`agb_c1`, `agb_c2`) and leaf area (`lai_b1`, `lai_b2`) take the height it gives; the wood area
(`wai_b1`, `wai_b2`), the sapwood (`sapwood_area_b1`, `sapwood_area_b2`) and the share of the wood above
ground (`aboveground_frac`) are per PFT too. A cohort reproduces once
its height reaches its PFT's `[pft].min_repro_height` [m], a height set outside MEDS.

Leaves are parameterized by their turnover: `[pft].leaf_turnover_rate_toc` [yr⁻¹, > 0] is the share
of a top-of-canopy cohort's leaf carbon its background turnover replaces a year, and with
`[trait_dynamics].trait_plasticity_on` it slows in shade (`docs/science/plant_traits.md` eq. (1)).
The light-plasticity slopes, `[pft].kplastic_vcmax`, `kplastic_rd`, `kplastic_sla` and
`kplastic_leaf_turnover`, ship as ED2's fits to Vcmax25 and the turnover rate (Lloyd et al. 2010), from
`scripts/pft_traits`.

## Several sites

A polygon is one site unless `[[sites]]` divides it, in a site run. Each `[[sites]]` table is one site:

| Key | Meaning |
|---|---|
| `area_fraction` | its share of the polygon (required; the shares sum to 1) |
| `elevation` [m] | default: the polygon's, `[site].elevation` |
| `slope`, `aspect` [degrees] | from horizontal, and clockwise from north the way the slope faces; default flat (0, 0) |
| `census_site_id` | the census rows it reads (`init_mode = 1`, where every site needs one) |
| any of `[soil]`'s column keys | `depth`, `layer_bottoms`, `retention`, the hydraulic and thermal texture, `psi_fc`, `water_bottom_bc`, `deep_temp`, `deep_depth`, under the same names |

- **The ground.** A table that sets any of `[soil]`'s column keys gives its site its own column,
  taking `[soil]`'s value for each key it leaves out; a table that sets none takes `[soil]`'s
  column. Each site's column is checked as `[soil]`'s is. The ground's surface (pond, dry surface
  layer, optics) and its thermal boundary are the run's. With the soil maps, every site takes the
  cell's soil, and the column keys are refused.
- **The terrain.** Each sub-step the polygon's forcing is moved to each site: to its elevation, by the
  run's lapse (`[site].apply_elevation_lapse`), and its direct beam onto its slope
  ([forcing](science/forcing.md) §8, Each site of a polygon). Each site's daily-mean air temperature
  drives its phenology and growth temperature, and its forcing is reported as `mean_*_si`.
- **What a site is.** Each site has its own patches, made at the start (`init_mode = 0`, the same
  number per site) or from its census rows; the patch operators work within a site
  ([order of processes](science/order_of_processes.md)), so no area moves between sites. The output
  reports each site's soil column and properties on a `site` axis, and every `_py` variable is the
  polygon's ([diagnostics](science/diagnostics.md) §2–3).
- **A restart** reads the sites back from the state file, which must hold as many as the config and,
  for each, as many soil layers.
- **A region refuses `[[sites]]`**: each cell is one site until a site map gives its sites.

```toml
[[sites]]                 # a quarter of the polygon on a shallow, south-east-facing slope
area_fraction = 0.25
depth         = 0.6
elevation     = 400.0
slope         = 12.0
aspect        = 135.0
[[sites]]                 # the rest on [soil]
area_fraction = 0.75
```

The `[site]` block is the polygon's place; the name `[[sites]]` differs because a TOML file cannot
hold both a table and an array of tables of one name.

## How a run starts

`[init].init_mode` selects one of three. The file a mode starts from is required by that mode and
optional otherwise; the Python census start (`meds.demography.Site(census=True)`) reads `census_file`
whatever the mode:

| Mode | Start from | Needs |
|---|---|---|
| `0` | **near-bare ground** (the default) | nothing; `[pft].nbg_init_density` seeds it (below) |
| `1` | a **cohort census** | `[init].census_file` — a CSV with one row per cohort, columns matched by name from its header: `patch_id`, `dbh` [cm], `pft` and `nplant` [plants per m² of the patch] are required; `patch_area`, `site_id`, `cohort_id` and `height` are optional. `dbh` drives the allometry. |
| `2` | a **state checkpoint** | `[init].restart_file` — a `<prefix>-S-*.nc` written by a previous run. Continues the exact instantaneous state: a run split at a restart ends where the unsplit run ends, bit for bit. |

Each distinct `patch_id` is a patch of age 0. With `patch_area`, in any unit, the patches take their
areas normalized to the site; without it they share it equally. A file without a header whose rows
are seven numbers is read in the old positional order, `site_id, patch_id, cohort_id, dbh, height,
pft, nplant`. An unknown or repeated column, a missing required one, or a patch whose rows disagree
on `patch_area` stops the run, naming it.

**A census stand is restructured before the first step**, by the slow step's own operators: its
monthly cohort fusion, fission and cull, then its yearly patch fusion, without recruitment or
disturbance, and under the same `[demography]` switches. So a census can carry one row per measured
tree size and one patch per plot cell, and the model starts from a stand within `max_cohort` and
`max_patch`, or above them where `patch_light_tol_max` keeps dissimilar patches apart. The run log
prints the counts before and after.

**A restart can take this run's leaf traits** instead of the state file's:
`[init].reacclimate_traits = true` (restart only; default false). The plastic traits (`sla`,
`vcmax25`, `rd25`, `leaf_turnover_rate`) are then set as a census start sets them from this run's PFT file:
acclimated to each cohort's leaf area above it, as the state file holds it, with
`[trait_dynamics].trait_plasticity_on`, and the PFT's top-of-canopy values without it. Leaf area stays
as read, and leaf carbon scales by the SLA's change, with storage taking up the difference; a lower
SLA takes at most the storage a cohort holds, and its leaf area then falls short by the rest. A restart
whose PFT file gives the traits the state already has changes nothing. A calibration
trial restarts from a shared state this way, so a changed `vcmax25` reaches the cohorts
(`scripts/calibrate_fast`).

**Near-bare ground is empty unless seeded.** Each patch starts without cohorts, and the stand
grows from recruitment. `[pft].nbg_init_density` [plants m⁻²], optional and 0 by default, gives
each patch one cohort of every PFT with a positive value. That many plants per m² start at
`min_cohort_height`, on the polygon's own height curve, where recruits are born. This is ED2's
`init_density`; a regional spin-up uses it to start every polygon with the same seedlings.

A census is how you start from a field inventory; see
[`data/census_example.csv`](../data/census_example.csv)
and `init_from_census` in [`../src/io/init/meds_init_census.f90`](../src/io/init/meds_init_census.f90).

**A census or state file that cannot start the run stops it**, naming the mode, the file and what
is wrong with it: one that cannot be opened, a census with no row for the site or a malformed
column or row, or a state file an older MEDS wrote, without a row for each of the run's polygons
(a region's, given to a site run) or with a PFT, site or soil-layer count unlike the run's. A state
file must hold every variable this MEDS reads: one from an older MEDS stops the run at the first it
lacks, so start that run again from its census or bare ground. Only `init_mode = 0` starts from near-bare ground; a run that asks for a stand never starts
from bare ground instead. `meds_main` stops, and `meds.model.Run` raises `RuntimeError`.

**The initial soil state** is `[init].soil_temp` [K] and `[init].soil_theta` [m³ m⁻³], both optional:
every soil layer of every patch starts there, unless a state checkpoint restores the soil. They
default to 288 K and 0.30. A run with no spin-up starts from them, so set them to the site: the
mean annual air temperature, for instance, keeps a warm site's column from starting cold. The
loader refuses a temperature outside 233–333 K and a water content outside
(`[soil].theta_res`, `theta_sat`]. A column from the soil maps has a texture per layer, so
each layer starts at its own field capacity instead, and `soil_theta` is refused.

**A new cohort starts in its patch's state**, at setup and at every recruitment: at the canopy-air
temperature, holding the water its tissues hold at the potential of the patch's top soil layer. At
setup the fast loop has not seeded the patches yet, so a cohort starts at `soil_temp` and the starting
soil water's potential, and is re-seeded from its patch once the patch is seeded.

## Output

Two streams, both set in `[output]`:

- **Diagnostic timeseries** — the `[output]` subsystem. 373 variables: each a mean over its
  record's period (`mean_<x>`), and a state also its value at the period's end (`<x>`), at four timescales
  (sub-daily, daily, monthly, annual). `[output].io_config` chooses them by name or pattern
  (`"*_co" = false`, `"mean_*_py" = "M Y"`); `meds_main --dump-io-config` lists every variable. See
  [`science/diagnostics.md`](science/diagnostics.md) §5.
- **State files** — `<state_file_dir>/<state_file_prefix>-S-<YYYYMMDDHHMMSS>.nc`, written with
  `state_write = true`: the model's whole state at an instant, no diagnostics, to restart from. One
  every `state_interval` (`"N days"`, `"N months"` or `"N years"`: every N days from the run's start,
  whole slow steps, or at every N-th month or year boundary) and one at the run's end. **The
  timestamp is the simulated date**, so pointing `[init].restart_file` at one resumes from exactly
  that date, and the run continues as the one that wrote it would have, bit for bit.

Every file of both streams has a `polygon` dimension, one row per polygon; a site run's files have one
row ([Regional runs](#regional-runs) gives the layout).

A checkpoint is raw prognostic state at an instant and a diagnostic record is a time average; the
two streams are deliberately separate because conflating them is how a restart silently starts from
a mean.

## The forcing file

`[forcing].forcing_on` turns the meteorological driver on; `[site]` says where the site is.
Without it the fast loop runs against a constant reference climate, which is useful for tests and
useless for science.

`[forcing].format` picks the source:

- **`"ERA5land"`** reads the global ED_ERA5land archive: one file per variable per month, built
  with the tools in `scripts/prepare_era5/`. Give it `data_path` (the archive folder) and
  `max_distance_km`. The reader finds the site's cell itself and takes that cell's elevation from
  the archive, so `path`, `grid_index`, `grid_match` and `[site].grid_elevation` do not apply and
  are rejected if present. Every month file the run needs must exist when it starts, and the run
  must step daily from midnight (`dt_slow = "1d"`, `start_time` at 00:00:00): the reader loads a
  month at a time before each step.
- **`"default"`** reads one forcing file, named by `path`: cut from ERA5-Land by
  `scripts/prepare_era5/make_forcing_file.py`, or built from flux-tower data by
  `scripts/prepare_flux_tower/make_tower_forcing.py`.

The file formats, the ERA5-Land and flux-tower preparation, and the recycling rules are documented in
[`science/forcing.md`](science/forcing.md). Three things to know before writing a config:

- **Every forcing clock is UTC.** A file must carry `time_zone = "UTC"`; the solar geometry takes
  local solar time from `[site].longitude`. Convert a local-time source when you build the file,
  and shift the output to local time afterwards.

- **The recycle window is declared, never inferred.** If `recycle = true`, then `recycle_start` and
  `recycle_end` are required, and the span must be an exact whole number of calendar years. A
  window of any other length drifts both hour-of-day and day-of-year on every wrap while the daily
  mean stays correct, so nothing downstream complains. The window may start at any record stamp,
  and the seam between its last record and its first may fall anywhere in a day. On an end-stamped
  file such as ERA5-Land, a calendar year's first record is at 01:00.
- **MEDS never gap-fills.** A missing or NaN required value is a hard error, not an interpolation.

**The forcing is moved to each patch's canopy-air top.** There is no fixed reference height:
every sample is taken from the forcing's own heights to the top of each patch's canopy air space, which
grows with the stand, and the aerodynamics runs from there ([`science/forcing.md`](science/forcing.md) §8).
- **Declare the forcing's own heights in `[forcing]`:** `tq_height` and `wind_height`; `height_above`
  (`"zero_plane"` for a reanalysis, `"ground"` for a flux tower); and `wind_exposure` (`"open_terrain"`
  for ERA5's 10 m wind, with `wind_exposure_z0 = 0.03` and `wind_blending_height = 40`, or `"local"`).
  `meds_config_main.toml` shows ERA5-Land's values.
- **The terrain lapse** (`[site].apply_elevation_lapse`) moves temperature, pressure, humidity (at constant
  relative humidity) and file longwave from the forcing cell's elevation to the site's.
  `[site].lapse_rate_tair` takes one rate or twelve monthly rates, and `[site].grid_elevation` gives
  the cell's elevation for an `"default"` file. Both are required with the lapse on and refused
  with it off; a flux tower measures at the site and runs with it off.
- **A file may state its heights.** When it carries `tq_height_m`, `wind_height_m` or
  `height_above`, they must agree with `[forcing]`, or the run stops at open.

**CO₂ is set in `[forcing]`, never by the met file.**
- **`co2_source = "const"`**, the default, holds `co2_const` for the whole run.
- **`co2_source = "file"`** reads `co2_file`, a MEDS CO₂ file looked up on model time, so the CO₂
  keeps rising while recycled met repeats.
  - The repository ships `data/co2/co2_cmip7_global_annual_1000-2022.txt`, CMIP7 global annual
    means for 1000–2022.
  - The path is relative to where MEDS runs, like `path`.
  - A run the file does not cover is refused at startup.
- **Rejected combinations.** The other mode's key is rejected, and so is a met file that carries
  `CO2air`.

The format (a `timestep` line, a `units` line, then one `<period start> <value>` row per period) and
the shipped series are in [`science/forcing.md`](science/forcing.md) §12.

## Regional runs

`[run].mode = "region"` simulates every selected cell of a box of the ED_ERA5land archive as its own
polygon, in one process. Each polygon is an independent stand at its cell's centre, at the cell's
orography, in UTC, and it computes exactly what a site run at that point would: a site run is a
region of one polygon, stepped the same way. Scattered sites are not a region; run them as separate
site runs (a job array, or one allocation filled with GNU parallel).

```toml
[run]
mode = "region"

[maps]                                               # required: each polygon's soil (see below)
soil_depth_map     = "ED_ERA5land_soil_depth_v2.nc"
soil_texture_map   = "ED_ERA5land_soil_texture_v2.nc"
soil_chemistry_map = "ED_ERA5land_soil_biogeochemical_properties_v2.nc"

[region]
box_nwse          = [42.95, -76.95, 41.95, -75.95]   # [N, W, S, E]; may cross 0 or 180 degrees
land_fraction_min = 0.5                              # skip lakes and fractional coastal cells
cell_mask         = "natural_mask.nc"                # optional: a map on the forcing grid ...
cell_mask_var     = "mask"                           # ... its (lat, lon) variable ...
cell_mask_min     = 0.5                              # ... and the value a cell needs to run
```

- **Selection.** Every cell of the box that the archive has data for and whose static land fraction
  is at least `land_fraction_min` (default 0.5) becomes a polygon, in row-major order. A polygon's
  id is its cell's index on the global grid, `row * 3600 + col` at 0.1 degrees, so ids agree across
  regions.
- **Cell mask.** With `cell_mask`, a cell of the box runs only where the mask's `cell_mask_var` is at
  least `cell_mask_min` (default 0.5), and is skipped where it is below or the map holds no value. A
  0/1 mask works as it is; a fraction, such as a natural-vegetation share, takes a threshold. The mask
  must be on the forcing grid: the file's `lat` and `lon` must be the ED_ERA5land static file's, or the
  run stops. MEDS never regrids a map; the run log says how many cells the mask kept.
- **Output.** One file per tier per time chunk for the whole region, with a `polygon` dimension and
  the coordinates `polygon_id`, `lat`, `lon`, `row` and `col`, as every run's files have (a site
  run's polygon is id 0). Region files hold every variable and tier a site run writes, the sub-daily one included; `[output].io_config` chooses them, and with
  them the volume. A cohort or patch variable is `(time, polygon, cohort)` or `(time, polygon,
  patch)`, its axis as long as the largest live count among the polygons when the file opens;
  `n_cohort` and `n_patch` `(time, polygon)` give each polygon's own count, and the entries past it
  are the fill value. Once a polygon has failed, the end of
  the run writes nothing more, so a closed file is never rewritten.
- **Each polygon's sites, their soils and map inputs in its files.** Every file has a `site` axis
  behind `polygon`, as long as the largest polygon's site count; `n_site (polygon)` gives each
  polygon's own, and the entries past it are the fill value. A region runs one site per cell. Sites
  have their own layers, so the files' soil coordinate is per site: `soil_z (polygon, site, soil)`
  holds each site's layer node depths, and its CF bounds `soil_z_bnds (polygon, site, soil, nv)` each
  layer's top and bottom, so the bottom layer ends at that site's own soil depth. Both are negative
  downward, and the fill value below the site's bottom layer, as its soil variables are. Each file
  also carries each site's area and terrain, `area_si` (its fraction of the polygon),
  `elevation_si` [m], `slope_si` and `aspect_si` [degrees], and the inputs the cell's maps gave
  it, as it runs them, all `(polygon, site)`: `n_soil_si` (its layer count), `soil_depth_si` [m],
  `bottom_bc_si` (`flag_values` 1, 2 for drain, aquifer), `drain_factor_si` (the drainage factor
  $f$ of the bottom flux $`q = f\,K`$: 1 on a free-draining base, the cell's `f_leak` on a leaky one,
  0 on a sealed one; the fill value on an aquifer) and `deep_temp_si` [K], under the same names as
  in the restart file; with `height_allometry = "chave2014"`, `chave_e_py (polygon)`.
- **Soil.** Every polygon takes its soil column from its own cell, through the `[maps]` soil
  maps, which a region requires. With `height_allometry = "chave2014"` it also needs
  `[maps].chave_e_map`. Each polygon is then a site run with that cell's soil, and its
  configuration is checked on its own.
- **Rules.** A region needs `[forcing].format = "ERA5land"` with forcing and the fast loop on. It
  takes its locations from its cells, so `[site].latitude`, `longitude`, `elevation`
  and `[forcing].max_distance_km` are refused; the rest of `[site]`, the terrain-lapse switch and
  rates, still applies. A region starts from bare ground (`init_mode = 0`) or from a region state
  file (`init_mode = 2`), never a census, and it runs without the fast probe. A region loads each
  month's forcing once, before the month, so a recycle window must start at 00:00 or 01:00 on the
  1st of a month.
- **Threads.** A region steps its polygons on `[run].n_threads` OpenMP threads, each polygon on one
  thread. A worker thread's stack is `OMP_STACKSIZE`, 4 MB by default with ifx, not the process's,
  and a grown stand's automatic arrays overflow it after a few decades of regrowth: the run dies with
  a segmentation fault. Set `OMP_STACKSIZE=512M` (it reserves address space; only what a thread
  touches is used). A region with more than one thread warns at open when it is not set.
- **State files.** With `[output].state_write`, a region writes one state file every
  `state_interval` and at the end, `<prefix>-S-<date>.nc` like a site's. It holds every polygon: a
  `polygon` axis leads each variable, `polygon_id` names the rows, the cohort and patch axes are as
  long as the largest polygon's, and `n_cohort`, `n_patch` and `n_window` mark each polygon's live
  part (the rest is the fill value); `lat` and `lon` are its cell's centre. A region restarted from
  it (`init_mode = 2`, `restart_file`) finds each polygon's row by its id, the cell's index on the
  forcing grid. The file must hold every polygon of the region and may hold more, so a region whose
  cell failed resumes from its last state file with that cell masked out (`[region].cell_mask`); the
  other polygons run on as they would have. A restart file that cannot be opened, one an older MEDS
  wrote, or one without a row for each polygon stops the region, naming the file. A restart
  continues exactly: a region run in two halves, across a recycle wrap, writes the second half as
  the whole run does. A site run writes the same layout with one row.
- **Threads.** In a region, `[run].n_threads` runs that many polygons side by side, each with a
  single-threaded patch loop; the results are the same at any thread count. A polygon that fails
  (a NaN, a cohort or soil-carbon state out of bounds) stops at the failing step, and the others
  finish the month. The month's records are then written, with the fill value where the failed
  polygon has none, and so is a state file that falls due at the month's end (the failed polygon as
  it stood); the log names the polygon and the date, and `meds_main` stops with an error. So the
  region ends after the month in which a polygon failed: no later month runs, and no end-of-run state
  file is written.
- **Cost.** Work and memory grow with the polygon count. Output is written between months, so a
  crash loses at most the current month. The 100 cells of a 1° box around Ithaca take 13 minutes a
  year on one thread and 81 s on 40 (measured with a detail polygon, before the region files held
  every variable). The hourly tier's cost grows with the polygons; choose its variables with
  `[output].io_config`. See `docs/dev_plans/MEDS_POLYGON_RUNTIME_PLAN.md` §14 for the production
  plan.

## Input maps

`[maps]` names the static maps a run takes its soil from, and Chave's site index E. Each map is
on the ED_ERA5land grid: its `lat` and `lon` must be the archive's static file's, so a cell's value
sits at the cell's own row and column, and MEDS never regrids. `src/io/maps/meds_maps.f90`
reads every map and is the only module that does.

| Key | File (ED-SoilMaps v2) | Variables |
|---|---|---|
| `soil_depth_map` | `ED_ERA5land_soil_depth_v2.nc` | `soil_depth` [m], `bc_class`, `f_leak`, `deep_temp` [K] |
| `soil_texture_map` | `ED_ERA5land_soil_texture_v2.nc` | `sand`, `silt`, `clay` [%], `bulk_density` [g cm⁻³], `coarse_fragments` [% volume], on SoilGrids' six layers (`depth_bnds` [m]) |
| `soil_chemistry_map` | `ED_ERA5land_soil_biogeochemical_properties_v2.nc` | `organic_carbon` [g kg⁻¹], `cec` [cmol(c) kg⁻¹], `ph_h2o`, on the same layers |
| `chave_e_map` | `ED_ERA5land_Chave2014_E.nc` | `chave_e` |

The three soil maps go together. A region needs them. A site may use them as well, read at the
forcing cell its `[site]` location reads (so with `format = "ERA5land"`), or it sets the
`[soil]` column by hand. From the maps, each cell's column is built as follows:

- **Depth and layers.** The column ends at `soil_depth`, on the default layers
  ([`[soil]`](#the-main-file)).
- **Each layer's soil.** A layer takes the map layers it overlaps, weighted by the overlap; below the
  maps' 2 m, it takes the deepest map layer. Its retention curve and Ks come from MEDS's
  pedotransfer function, and its thermal properties from the same texture by CLM5's mineral-soil
  formulas ([`science/soil_biophysics.md`](science/soil_biophysics.md)).
- **The bottom water boundary** follows `bc_class`. The first three drain through the base at
  $`q = f\,K(\theta_n)`$ ([`science/soil_biophysics.md`](science/soil_biophysics.md), Boundary
  conditions):
  - 1 is `free_drain`, $f = 1$;
  - 2 (leaky) drains at the cell's `f_leak` (0.05 in ED-SoilMaps v2), which must be in (0, 1];
  - 3 (sealed) is `bedrock`, $f = 0$;
  - 4 is `aquifer`.
- **The thermal anchor** sits at `deep_temp`, annual damping depth / √2 below the bottom node, where it
  reflects the annual wave least. It is used when `[soil].thermal_bottom_bc = "dirichlet"`.

So that a run never holds two values for one quantity, every key the soil maps provide is refused
while they are set:

- in `[soil]`: `depth`, `layer_bottoms`, `retention` and the texture keys, `water_bottom_bc`, and
  `deep_temp` and `deep_depth`;
- `[init].soil_theta`.

`[soil].psi_fc` still defines field capacity, where each layer starts. With `chave_e_map`,
`[allometry].chave_e` is refused in the PFT file, and the map needs a PFT on `height_allometry = "chave2014"`.

**MEDS fills no gaps.** A run cell a map leaves empty stops the run, naming the map, the variable and
the cell. A product is completed in its own pipeline, not in the model.

## Worked examples

The example configs are the fastest way in, and each is a complete pair:

- [`examples/example03_demography/`](../examples/example03_demography/) — the demography engine
  driven from Python by growth, mortality and recruitment trained on the BCI censuses, from the
  1985 census and from near-bare ground.
- [`examples/example04_column_biophysics/`](../examples/example04_column_biophysics/) — the
  coupled model at the Barro Colorado Island flux tower: forcing from the tower, a start from the
  plot's census, five years at hourly output, and ten days restarted at half-hourly output.
