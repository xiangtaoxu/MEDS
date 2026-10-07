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
  edges, the mortality-hazard coefficients from wood density. Set `[options].override_derived =
  true` and supply a `[derived]` block to pin them instead.
- **Optional keys.** Solver settings, process switches and a few late-added blocks (`[soil_column]`,
  `[snow]`, `[soil_carbon]`) read each key with a fallback to its default, so a feature needs no
  TOML edits. The reference file lists each one commented out at that default.

**Every key MEDS reads is listed in the two reference configs**, `meds_config_main.toml` and
`meds_config_pft.toml`, set or commented out at its default. They are the complete list. A key that
is not listed in its file's reference stops the run before it starts: the error names every such key
at once, together with any missing required ones. Each line gives a suggestion:
- a retired key gets what replaced it;
- a key that moved section gets its new home (`state.cohort_max` → `output.cohort_max`);
- a misspelling gets the key it most likely meant.

A key that parses and does nothing is worse than one that is absent: `[forcing] dt_forcing = 7200`
(the key is `timestep`) used to run silently at the file's own spacing. The test
`config_keys_listed` holds the two references equal to the keys the loader reads, so a new key is
listed when it is added.

Every run that writes a state writes `<output_dir>/<prefix>_pft_parameters.csv` — one row per PFT
with every per-PFT trait actually used. It is a provenance record: what the run really did, not what
you meant.

Every run also writes the **parameter record**, `<prefix>_parameters.csv`, beside its diagnostic
output (`[output].dir`, `[output].prefix`) and beside its state (`[state]`), whichever it writes. It
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
| `[fast]` | The sub-daily loop: `dt_fast`, the integrator, tolerances, error control. |
| `[init]` | How the run starts, and the path to the PFT file. |
| `[demography]` | Cohort and patch fusion/fission, the cadence switches. |
| `[disturbance]`, `[recruitment]` | Treefall rates; the seed-rain and reproduction settings. |
| `[soil_column]` | The **physical ground**: layer count, depth, grid growth, hydraulic texture, retention family, thermal properties. |
| `[soil]`, `[energy]`, `[snow]`, `[aerodynamics]` | The **solvers over it**: selectors and tolerances. `[soil]` also holds the bare ground's optics: `ground_albedo_vis`, `ground_albedo_nir` and `ground_emissivity`. |
| `[soil_carbon]` | The CENTURY decomposition: selectors, rate parameters, cold-start spin-up. |
| `[hydraulics]` | Plant water transport: the pressure–volume and vulnerability traits, the conductance form (`conductance = "whole_plant"` or `"segment"`), and the root profile (`root_beta`, `root_depth`). The root profile is a plant trait, so `[soil_column].root_beta` is refused. |
| `[forcing]`, `[site]` | The meteorological driver, and where the site is. |
| `[region]` | For `mode = "region"` only: the box of forcing cells and which of them to simulate. |
| `[output]` | Which diagnostics are written, on which axes, at which timescales. |
| `[state]` | Restart checkpointing: output directory, prefix, interval. *(Called `[io]` before v0.3.0; that spelling is now refused.)* |
| `[options]` | `override_derived` and other run switches. |

**`[soil_column]` is the ground; `[soil]` is the solver over it.** The split matters: a layer count
above the compile-time ceiling or a saturated water content below the residual produces a silently
wrong column rather than a crash, so the physical block is validated at load.

### The bottom thermal boundary — `[energy].bottom_bc`

The default 2.0 m column sits against an annual thermal damping depth of about 2 m, so the annual
temperature wave has barely attenuated by the time it reaches the base. `bottom_bc` decides what
happens to it there.

| key | values | meaning |
|---|---|---|
| `bottom_bc` | `geothermal` (default) \| `dirichlet` | prescribed bottom flux (held at zero — an adiabatic wall that **reflects** the annual wave) vs conduction to a fixed deep temperature |
| `deep_temp` | [K] | the anchor temperature. **Required** when `bottom_bc = "dirichlet"`; there is no default |
| `deep_depth` | [m], default 3.12 | depth of the anchor plane below the surface; must lie below the bottom soil node |

`deep_temp` is the mean annual soil temperature below the damping depth — close to the mean annual
air temperature of the forcing that drives the run, and a site property like latitude. It is required
rather than defaulted because an error in it is a steady flux into the column base, and so a
mean-annual bias in deep-soil temperature; a silent default would put back, somewhere else, the bias
this boundary condition exists to remove.

`deep_depth` is a physical choice, not a numerical one, and the default is derived: a resistive
termination reflects least when the anchor sits `d/sqrt(2)` below the bottom node, `d` being the
annual damping depth. For the default column that is 1.727 + 1.973/√2 = 3.12 m. **Recompute it if you
change `[soil_column].depth` or the thermal texture** — the derivation, the measured sweep and what
the anchor does and does not fix are in [`docs/science/soil_biophysics.md`](science/soil_biophysics.md).

Deepening the column with `[soil_column].depth` remains available and is the other way out, at the
cost of layers; it also moves the root profile and the drainage, so it is not a thermal-only change.

### Cohort and patch fusion — `[demography]`

Fusion merges similar cohorts within a patch every month, and similar patches every year. Both run
in passes. Each pass merges every pair that passes a similarity test at one tolerance, and the
tolerance grows geometrically from pass to pass, from a minimum to a **ceiling**. Passes stop early
once the count is within its **target**.

| | Cohorts | Patches |
|---|---|---|
| similar when | their heights differ by less than `tol · hgt_max`, and their combined LAI is below `cohort_lai_cap` | their light profiles differ by at most `tol` on average, and by at most `tol · patch_light_maxdev_factor` in any layer |
| tolerance, first pass | `cohort_size_tol_min` | `patch_light_tol` |
| tolerance ceiling | `cohort_size_tol_max` | `patch_light_tol_max` (default 0.15) |
| passes, at most | `n_cohort_fusion_iter` | `n_patch_fusion_iter` |
| target | `max_cohort`, in the most crowded patch | `max_patch` |

A patch's light profile is the fraction of full sunlight left under each of `n_height_layers` equal
layers up to the tallest `hgt_max`, `exp(−light_ext × LAI at and above the layer)`. Layers fully lit
in both patches are not compared.

**The ceiling is hard and the target is not.** Two cohorts or patches more different than the
ceiling stay apart even if the count is then above its target, because merging them would lose
structure the model needs. A run in which that happened says so at the end, for example
`NOTE: patch fusion left 24 patches on 2013-01-01 (max_patch = 12): ...`. A negative `max_cohort` or
`max_patch` forces fusion instead, ignoring the test.

`patch_light_tol_max` is optional. Absent, it is 0.15, or `patch_light_tol` if that is larger. Set
it to `patch_light_tol × 1.5^(n_patch_fusion_iter − 1)`, 0.759375 for the shipped 0.10 and six
passes, to reproduce the schedule before the ceiling existed.

### The PFT file

Everything PFT-specific: the `[pft]` trait table, the `[camac]` mortality-hazard derivation
coefficients, and the `[allometry]` coefficients. **The number of PFTs is the length of
`[pft].wood_density`** — every other trait array must match it.

Wood density is the primary PFT axis. The Camac-2018 mortality coefficients are derived from it as
power laws, so low-density PFTs get a higher growth-independent hazard and a steeper low-growth
penalty. Growth, competition and reproduction parameters are per-PFT but ship uniform.

## How a run starts

`[init].init_mode` selects one of three, and the files for the unselected modes are ignored rather
than removed:

| Mode | Start from | Needs |
|---|---|---|
| `0` | **near-bare ground** (the default) | nothing |
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
`vcmax25`, `rd25`, leaf lifespan) are then set as a census start sets them from this run's PFT file:
acclimated to each cohort's leaf area above it, as the state file holds it, with
`[trait_dynamics].trait_plasticity_on`, and the PFT's top-of-canopy values without it. Leaf area stays
as read, and leaf carbon scales by the SLA's change, with storage taking up the difference; a restart
whose PFT file gives the traits the state already has changes nothing. A calibration
trial restarts from a shared state this way, so a changed `vcmax25` reaches the cohorts
(`scripts/calibrate_fast`).

A census is how you start from a field inventory; see
[`examples/example_demography/census_example.csv`](../examples/example_demography/census_example.csv)
and `init_from_census` in [`../src/init/meds_init.f90`](../src/init/meds_init.f90). Unusable input
falls back to near-bare ground with a warning.

**The initial soil state** is `[init].soil_temp` [K] and `[init].soil_theta` [m³ m⁻³], both optional:
every soil layer of every patch starts there, unless a state checkpoint restores the soil. They
default to 288 K and 0.30. A run with no spin-up starts from them, so set them to the site: the
mean annual air temperature, for instance, keeps a warm site's column from starting cold. The
loader refuses a temperature outside 233–333 K and a water content outside
(`[soil_column].theta_res`, `theta_sat`].

## Output

Two streams, both with the stem `<output_dir>/<output_prefix>` from `[state]`:

- **Diagnostic timeseries** — the `[output]` subsystem. 252 variables across 9 groups and 7 axes,
  each switchable individually per timescale (sub-daily, daily, monthly, annual). Run
  `meds_main --dump-io-config` to generate a file listing every available variable name; the
  override mechanism always worked, what was missing was any way to learn what exists.
  See [`science/diagnostics.md`](science/diagnostics.md).
- **State checkpoints** — `<prefix>-S-<YYYYMMDDHHMMSS>.nc`, enabled with `[state].write_state`. The
  instantaneous prognostic state only, no diagnostics, written every `[state].interval_years` and at
  run end. **The timestamp is the simulated date**, so pointing `[init].restart_file` at one
  resumes from exactly that date.

A checkpoint is raw prognostic state at an instant and a diagnostic record is a time average; the
two streams are deliberately separate because conflating them is how a restart silently starts from
a mean.

> **`[io]` was renamed to `[state]`.** The block was named for a legacy diagnostic writer that was
> retired at v0.1; what remained was the restart stream, so `io` named the one output path it did
> *not* cover. Since v0.3.0 a config that still spells it `[io]` is refused at load, with a message
> naming the new keys: rename the block, and `io.state_interval_years` to `state.interval_years`
> (the `state_` prefix was stuttering once the block itself was called `state`).

## The forcing file

`[forcing].forcing_on` turns the meteorological driver on; `[site]` says where the site is.
Without it the fast loop runs against a constant reference climate, which is useful for tests and
useless for science.

`[forcing].format` picks the source:

- **`"ED_ERA5land"`** reads the global ED_ERA5land archive: one file per variable per month, built
  with the tools in `scripts/prepare_era5/`. Give it `data_path` (the archive folder) and
  `max_distance_km`. The reader finds the site's cell itself and takes that cell's elevation from
  the archive, so `path`, `grid_index`, `grid_match` and `[site].grid_elevation` do not apply and
  are rejected if present. Every month file the run needs must exist when it starts, and the run
  must step daily from midnight (`dt_slow = "1d"`, `start_time` at 00:00:00): the reader loads a
  month at a time before each step.
- **`"ED_default"`** reads one forcing file, named by `path`: cut from ERA5-Land by
  `scripts/prepare_era5/make_forcing_file.py`, or built from flux-tower data by
  `scripts/prepare_flux_tower/make_tower_forcing.py`.

The earlier names, `"netcdf"` and `"era5land"`, are refused with a message naming the new one.
The file formats, the ERA5-Land and flux-tower preparation, and the recycling rules are documented in
[`science/forcing.md`](science/forcing.md). Three things to know before writing a config:

- **Every forcing clock is UTC.** A file must carry `time_zone = "UTC"`; the solar geometry takes
  local solar time from `[site].longitude`. Convert a local-time source when you build the file,
  and shift the output to local time afterwards. There is no `[site].utc_offset` or
  `apply_solar_longitude`, and both are refused.

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
  the cell's elevation for an `"ED_default"` file. Both are required with the lapse on and refused
  with it off; a flux tower measures at the site and runs with it off.
- **A file may state its heights.** When it carries `tq_height_m`, `wind_height_m` or
  `height_above`, they must agree with `[forcing]`, or the run stops at open.
- The old `[site].reference_height`, `wind_meas_height`, `apply_wind_profile` and `wind_roughness_z0`
  are rejected, with a message naming these keys.

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
orography, in UTC, and it computes exactly what a site run at that point would: the region steps
every polygon with the site run's own step. Scattered sites are not a region; run them as separate
site runs (a job array, or one allocation filled with GNU parallel).

```toml
[run]
mode = "region"

[region]
box_nwse          = [42.95, -76.95, 41.95, -75.95]   # [N, W, S, E]; may cross 0 or 180 degrees
land_fraction_min = 0.5                              # skip lakes and fractional coastal cells
detail_polygons   = [1714634]                        # optional: these also write single-site files
```

- **Selection.** Every cell of the box that the archive has data for and whose static land fraction
  is at least `land_fraction_min` (default 0.5) becomes a polygon, in row-major order. A polygon's
  id is its cell's index on the global grid, `row * 3600 + col` at 0.1 degrees, so ids agree across
  regions.
- **Output.** One file per tier per time chunk for the whole region, with a `polygon` dimension and
  the coordinates `polygon_id`, `lat`, `lon`, `row` and `col`. Region files hold the variables that
  have the same shape everywhere: site totals, per PFT, per size class and per soil layer. Cohort
  and patch variables, and the sub-daily tier, are written only for `detail_polygons`, as ordinary
  single-site files named `<prefix>-p<polygon id>-...`.
- **Rules.** A region needs `[forcing].format = "ED_ERA5land"` with forcing and the fast loop on. It
  takes its locations from its cells, so `[site].latitude`, `longitude`, `elevation`
  and `[forcing].max_distance_km` are refused; the rest of `[site]`, the terrain-lapse switch and
  rates, still applies. Until restarts of regions exist, a region starts from bare ground
  (`init_mode = 0`) and writes no checkpoints (`[state].write_state = false`), and it runs without
  the fast probe. A region loads each month's forcing once, before the month, so a recycle window
  must start at 00:00 or 01:00 on the 1st of a month.
- **Threads.** In a region, `[run].n_threads` runs that many polygons side by side, each with a
  single-threaded patch loop; the results are the same at any thread count. A polygon that fails
  (a NaN, an impossible soil-carbon pool) is reported and stops, and the others finish the month.
- **Cost.** Work and memory grow with the polygon count. Output is written between months, so a
  crash loses at most the current month. The 100 cells of a 1° box around Ithaca take 13 minutes a
  year on one thread and 81 s on 40. A detail polygon costs about four times an ordinary one (its
  hourly site files and per-cohort diagnostics), so on many threads it sets the month's pace: the
  same run without one takes 53 s on 40 threads. See `docs/dev_plans/MEDS_POLYGON_RUNTIME_PLAN.md`
  for the roadmap to restarts and tiles.

## Worked examples

The example configs are the fastest way in, and each is a complete pair:

- [`examples/example_demography/`](../examples/example_demography/) — a 250-year demographic
  spin-up from near-bare ground.
- [`examples/example04_column_biophysics/`](../examples/example04_column_biophysics/) — the
  coupled model at the Barro Colorado Island flux tower: forcing from the tower, a start from the
  plot's census, five years at hourly output, and ten days restarted at half-hourly output.
