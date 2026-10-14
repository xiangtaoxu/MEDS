# MEDS BCI census initialization plan

**Status:** written 2026-09-28 against `beta` at `73ab654`, revised the same day after review.
**First trial implemented 2026-09-29** on `feat/bci-census-init`: P0a–P0c, P1 and the P2 example,
with the five-year run from the census finishing within the flux period and compared with the tower
(`examples/example04_column_biophysics/README.md`). G2 passes (the model's first stand equals the census
file: LAI 5.60, AGB 16.12 kgC m⁻²), G4 matches the simulation (25 patches against 24), and G8 is 7.2
minutes. `plot_census.py` and the remaining gates are left for the next round. It builds on the
allometry fix merged into `beta` as #321 (§5.2): every number below uses those defaults. The 2010
census was measured from the local copy in `/ibstorage/xiangtao/bci_census/` (§4), and MEDS's patch
fusion was simulated on it (§6.3). The user took the open decisions on 2026-09-28 (§11). Gates G2
and G4–G9 need model runs. The leaf-area half of #321 stops the biophysics example's spin-up from growing
any stand (§5.2). The user accepted that on 2026-09-29: the biophysics example is retired once this
one passes its gates (§11). The same PFT carbon balance is this example's main risk (G9).

**Goal:** the BCI flux-tower example stops spinning up from bare ground. It starts on 2012-08-01
from the stand the public Barro Colorado Island 50-ha plot census measured, and runs the five tower
years once. Every tree is the example's one evergreen PFT: functional diversity is out of scope.

## 1. Decisions

| # | Decision |
|---|---|
| D1 | **Initialize from census 7 (2010)**, the last census before the run starts. The run window stays 2012-08-01 → 2017-08-01 UTC. There is no comparison with a later census in this example. |
| D2 | **Every live tree with a measured diameter of at least 1 cm**, from the census-7 tree table (one row per tree). The diameter is the table's `dbh`, which is the main stem only: the local copy has no `ba` column, and 12.5 % of trees have more than one stem. Every tree is PFT 1. |
| D3 | **The allometry is MEDS's default**, after the fix that makes the default biomass law Chave et al. (2014) in carbon and divides the default leaf-area scale by `C2B` as ED2 does (§5.2). The example changes no allometry relationship. |
| D4 | **A patch is one plot cell, 20 m by default and any size allowed.** The file gives each patch its area, so a grid that does not tile the plot keeps its edge cells at their true area. |
| D5 | **The preparation only maps census data to patches: `patch_id`, `patch_area` and `nplant`.** It does no aggregation or fusion. Trees in one cell with the same measured diameter share a row, with `nplant` their count over the patch area, because they are the same cohort. |
| D6 | **MEDS fuses the census stand before the model runs** (P0c). After reading a census it applies the slow step's own cohort and patch restructuring, so the first step starts from a stand within `max_cohort` and `max_patch`. |
| D7 | **One MEDS run.** `meds_config_spinup.toml`, its restart file and the spin-up step of `run_example.py` are deleted. |
| D8 | **Three Fortran changes, none of which alters another run:** `[init]` keys for the initial soil state, off by default (P0a); a header-matched census reader with an optional `patch_area` column, which still reads today's files (P0b); and restructuring at census initialization (P0c). |
| D9 | **Soil carbon starts from the steady-state solve** that already exists (`[soil_carbon].spinup_steady`). The litter input is estimated from the census stand by the tool and checked against the model's first-year litter output. |
| D10 | **No data in the repository.** The tool reads the census in place, from the path `bci_census.toml` declares (§4), and checks it against the recorded checksum. The README carries the citation and acknowledgement the data asks for. |
| D11 | **The converter is general, the settings are BCI's.** `scripts/prepare_census/` turns a ForestGEO tree table into a MEDS census CSV; `examples/example04_column_biophysics/bci_census.toml` declares the BCI source, plot geometry and choices, as `bci_site.toml` does for the tower. |

## 2. What the spin-up supplies today, and what replaces it

| The spin-up supplied | Replaced by |
|---|---|
| vegetation structure: LAI 4.8, AGB 15.3 kgC m⁻², 114 cohorts after 50 years (under the uncorrected allometry) | the census stand (D1–D5) |
| a patch mosaic from 50 years of disturbance | MEDS's own fusion of the census cells, about 9 patches of unequal area (§6.3). Every patch has age 0 (§3, F1) |
| soil water and temperature in equilibrium with the forcing | `[init].soil_temp` and `[init].soil_theta` (P0a) |
| soil-carbon pools after 50 years of litter, not in steady state | the steady-state solve with a census-based litter input (D9) |
| cohort memories: growth average, phenology drives, leaf-water extremes | birth defaults. Each is already handled for new cohorts (§3, F7) |

The example gets faster: the 11-minute spin-up disappears. The evaluation run carries more cohorts
than the spin-up's 114, up to `max_cohort = 60` in each patch (gate G8).

## 3. Findings in MEDS that shape the plan (read 2026-09-28)

- **F1. The census path exists, with equal areas only.** `init_mode = 1` reads a CSV
  (`src/init/meds_init.f90`, `init_from_census`). It reads seven columns by position,
  `site_id, patch_id, cohort_id, dbh, height, pft, nplant`, with no area column. Each distinct
  `patch_id` becomes an equal-area primary patch of age 0. `dbh` [cm] drives the allometry;
  `height` and `cohort_id` are read for provenance only. `nplant` is plants per m². A test exists
  (`test/test_init_census.f90`).
- **F2. Nothing fuses at initialization.** Cohort restructuring first runs at the first month
  boundary, and patch restructuring at the first new year. A census read as it stands would run
  its first month with every row a cohort and its first five months with every cell a patch. P0c
  closes that.
- **F3. Cohort fusion compares heights, not diameters.** Two cohorts fuse when their height
  difference is below `hgt_max × tol`, starting at `cohort_size_tol_min = 0.02`, which is 0.92 m for
  `hgt_max = 46 m`. They never fuse if their combined LAI reaches `cohort_lai_cap = 1.0`. Fission
  splits any cohort above that cap. Trees above 116.7 cm share the 46 m cap, so within a patch they
  merge into one cohort unless the LAI cap keeps them apart.
- **F4. The initial soil state is hard-coded.** `init_fast_reservoirs` sets every layer to
  θ = 0.30 m³ m⁻³ and 288 K (`fast_context_t` in `src/fast_dynamics/driver/meds_fast_dynamics.f90`),
  with no config key. The example's Dirichlet bottom is 298.65 K, so the column would start
  10.65 K cold. Warming 2 m of soil by that much takes about 70 MJ m⁻², on the order of 13 W m⁻²
  for two months, which would bias the first months' surface fluxes.
- **F5. Soil carbon starts at zero** without `spinup_steady`: every `soil_carbon_t` pool defaults to
  0, so heterotrophic respiration and NEE would start at zero.
- **F6. The default allometry at BCI sizes** (defaults after #321, wood density 0.60):

  | dbh [cm] | height [m] | leaf area [m²] | AGB [kgC] |
  |---|---|---|---|
  | 1 | 3.1 | 0.5 | 0.06 |
  | 10 | 11.5 | 21 | 20 |
  | 30 | 21.4 | 130 | 310 |
  | 60 | 31.6 | 407 | 1,758 |
  | 100 | 42.2 | 944 | 6,315 |
  | 200 | 46.0 | 2,427 | 26,603 |

  A 1 cm stem is taller than `min_cohort_height = 2 m`, so no census tree is below the recruit size.
- **F7. New cohorts are already handled.** Camac mortality uses only its instantaneous term until a
  cohort has a growth average. Leaf water, hydraulic memory and phenology are seeded on first
  touch. Trait plasticity is off in this example, so the census path's instantaneous acclimation
  does not apply.
- **F8. The litter outputs exist.** `litter_leaf_site`, `litter_fineroot_site` and
  `litter_struct_site` check the tool's litter estimate against the model's first year (D9).
- **F9. The census stand under the default allometry** (measured on the 2010 census, 207,259 trees):

  | Allometry | stand AGB [kgC m⁻²] | stand LAI |
  |---|---|---|
  | the defaults before #321 | 39.8 | 11.2 |
  | **the defaults after #321 (D3)** | **16.1** | **5.6** |
  | the census's own `agb`, at a carbon fraction of 0.5 | 15.1 | |

  The stand's AGB is 7 % above the census's own estimate, which comes from the default heights and
  one wood density for every species. The LAI is close to the 4.8 the old spin-up settled at. Leaf
  area by size class:

  | dbh [cm] | stems ha⁻¹ | LAI |
  |---|---|---|
  | 1–10 | 3,729 | 1.34 |
  | 10–30 | 335 | 1.59 |
  | 30–100 | 78 | 2.19 |
  | ≥ 100 | 3.0 | 0.47 |
- **F10. The restructuring operators work on any site.** `new_fuse_cohorts`, `terminate_cohorts`,
  `split_cohorts`, `sort_cohorts`, `sort_patches`, `new_fuse_patches` and `terminate_patches` take a
  `site_t` and the config, and `src/init` already calls `sort_cohorts`. Cohort storage grows by half
  when full, so 85,000 rows load without quadratic copying. A cohort carries 59 per-cohort arrays and
  about 45 diagnostic fields, roughly 1.5 kB, so 85,000 cohorts take about 130 MB before fusion.
  Patches are sorted by age, and every census patch has age 0, so fusion visits them in file order.

## 4. The data

**The copy this example uses** is `/ibstorage/xiangtao/bci_census/bci_<year>.csv`: seven censuses,
1982–2010, each with 394,658 rows, one per tree, in the same tree order in every file. The example
needs census 6 (2005, for the mortality rate in §6.5) and census 7 (2010).

| File | SHA-256 |
|---|---|
| `bci_2010.csv` | `eb5b8373de32118250ec5955b91df122223867c9046af5e91f2085836158bf7b` |
| `bci_2005.csv` | `4794db4da9017d0a13b4b0acc2540e91aae5adfe2e43f62fae16c366696473ae` |

The columns the tool uses:
- `gx`, `gy`: metres from the west and south borders, in [0, 1000) and [0, 500);
- `dbh`: **mm**, of one stem only, the main stem;
- `agb`: aboveground biomass of the tree, **Mg** dry mass, 0 for dead trees;
- `status`: A alive, D dead, P not yet recruited, M missed; `nostems`, `hom`, `ExactDate`.

There is no `ba` column. The tool keeps status A, counts everything it drops, and reports the census
dates from `ExactDate`.

**What the 2010 file holds** (measured 2026-09-28):

| | |
|---|---|
| measured | 2010-01-18 to 2011-03-28, median 2010-05-17 |
| live trees | 221,758, of which 14,499 have no diameter and zero `agb`, and are dropped |
| used | 207,259 trees: 4,145 stems ha⁻¹, 416 ha⁻¹ at ≥ 10 cm, main-stem basal area 30.5 m² ha⁻¹ |
| multi-stem trees | 25,985 of the 207,259 |
| largest | 248 cm; 92 trees at or above 116.7 cm, where the default height reaches its 46 m cap |
| AGB, the table's own | 302.4 Mg ha⁻¹ (15.1 kgC m⁻²) |
| 2005 → 2010 | AGB 305.5 → 302.4 Mg ha⁻¹; 1.90 % of standing biomass and 2.67 % of stems die per year |

**Which release this is.** The canonical public release is Condit et al. (2019), *Complete data from
the Barro Colorado 50-ha plot: 423617 trees, 35 years*, Dryad doi:10.15146/5xcp-0d46, CC0, with eight
censuses to 2015 and a `ba` column. The local copy has seven censuses, 394,658 trees, no `ba`, and the
extra columns `pom`, `MeasureID` and `CensusID`. That matches the 2012 Smithsonian release the Dryad
version replaced; the match is inferred, not confirmed. The README states which release the figures
come from and cites Condit et al. The dataset asks that the PIs, Condit and Hubbell, be told of papers
using it.

**Getting it elsewhere.** A user without the local copy points `bci_census.toml` at their own. Dryad
blocks scripted downloads: its API answers 401, "must have current bearer token", and its web file
link returns a bot-challenge page to a script. The public tables come as `bci.tree.zip` (48,394,285
bytes, SHA-256 `07c0fce4528c18d11916c2dd2e8e9f069f08bf62dfee38a67b18bda5c3b987ec`), downloaded in a
browser. Its tables are R `.rdata`, read with `pyreadr`, an optional dependency of this example only.

**The plot and the tower.** The plan assumes the 50-ha plot's forest stands for the tower's
footprint. That is not checked here.

## 5. P0 — the model changes

### 5.1 P0a: the initial soil state, in `[init]`

Two optional keys, read with defaults so that no existing config changes:

| Key | Units | Default | BCI value | Why |
|---|---|---|---|---|
| `[init].soil_temp` | K | 288.0 | 298.65 | the tower's mean air temperature, already the Dirichlet `deep_temp` |
| `[init].soil_theta` | m³ m⁻³ | 0.30 | 0.30 | head −0.51 m on the default column: wetter than field capacity (0.164 at −3.37 m) and below saturation (0.43). August 1 is three months into the wet season |

- `build_fast_context` copies them into `theta_init` and `soil_temp_init`; `init_fast_reservoirs`
  is unchanged.
- The loader refuses `soil_theta` outside (`theta_res`, `theta_sat`] and `soil_temp` outside
  233–333 K, naming the key.
- A state restart that restores the fast reservoirs ignores them, as it ignores `census_file`.
- The canopy air, leaves and wood still start at about 288 K and reach the forcing within the first
  day. That is left alone.
- Docs: a paragraph in `docs/configuration.md`, "How a run starts". CHANGELOG entry.

### 5.2 The allometry fix this plan builds on

Merged into `beta` as #321 on 2026-09-28 (merge commit `32ccf9a`), as its own change:
- MEDS's default `agb_c1, agb_c2` were ED2's `IALLOM = 3` refit of Chave (`c14f15_bs_tf`) and
  `lai_b1` ED2's BAAD leaf fit (`c14f15_bl_xx`). Both are dry mass or its equivalent, and ED2
  divides both by `C2B = 2` before use (`size2bd`, `size2bl`). MEDS used them undivided, so its AGB
  was dry mass labelled as carbon and its leaf area twice ED2's.
- The fix makes the biomass law Chave et al. (2014) eq. 4 itself, in carbon: `agb_c1 = 0.0673/2 =
  0.03365`, `agb_c2 = 0.976`. It halves `lai_b1` to 0.23384770. Both land in `meds_config_pft.toml`,
  every example's PFT file (this example's included) and the `meds_allometry` initializers.
- A new test, `allometry_defaults`, checks that the shipped `[allometry]` block equals the
  initializers, that the biomass law is Chave's in carbon, and that the leaf scale is ED2's over
  `C2B`. It fails on the old values.
- ifx Release and gfortran Release: 55 of 55 tests pass, the new one included.

**Found after the merge (2026-09-29): the halved leaf area stops small trees growing.** The
biophysics example's 50-year spin-up (`meds_config_spinup.toml`, Ithaca) grows a stand before
#321 and none after it. Each half of #321 was then applied alone:

| Allometry | Final LAI | Final AGB [kgC m⁻²] | Mean dbh [cm] |
|---|---|---|---|
| before #321 | 4.10 | 9.63 | 24.7 |
| Chave biomass only | 4.94 | 6.77 | 31.6 |
| halved leaf area only | 0.00 | 0.002 | 0.45 |
| #321 (both) | 0.00 | 0.001 | 0.45 |

The recruits appear at 0.45 cm and never grow. The example PFT's carbon balance was tuned with the
doubled leaf area, and at half the leaf area a small tree cannot stay positive. The test suite
passed because no test runs a stand long enough to grow.

ED2 was read from `~/ED2_tutorial/ED2`, commit `f8d8228a` (2020-08-12). The sibling `../ED2` that
`CLAUDE.md` names does not exist on this machine. ED2's source carries a comment next to the leaf
line asking whether the `C2B` belongs there; the fix follows ED2 as it runs.

### 5.3 P0b: a header-matched census reader, with `patch_area`

`init_from_census` reads its columns by name from the header line instead of by position:

| Column | Required | Meaning |
|---|---|---|
| `patch_id` | yes | integer; each distinct value is one patch |
| `dbh` | yes | [cm]; drives the allometry |
| `pft` | yes | PFT index |
| `nplant` | yes | [plants per m² of the patch] |
| `patch_area` | no | [m²]; the same on every row of a patch. Absent, the patches have equal areas, as today |
| `site_id` | no | as today: the first site's rows are read. Absent, every row is one site |
| `cohort_id`, `height` | no | provenance only, as today; `height`, if present, must be positive |

- The areas are normalized by their sum, so the patches' fractions of the site add to 1 and the
  column's unit cancels. A patch whose rows disagree on `patch_area`, a non-positive area, a
  missing required column or an unknown column name stops the run, naming it.
- Today's seven-column files, `examples/example_demography/census_example.csv` and the test's,
  load unchanged.
- Docs: the census row of the `docs/configuration.md` "How a run starts" table. CHANGELOG entry.

### 5.4 P0c: restructuring at census initialization

After `init_from_census` succeeds, and before `polygon_prepare` seeds the fast reservoirs, the driver
applies the slow step's own restructuring to the census stand, minus recruitment and disturbance:
1. the monthly cohort block: `new_fuse_cohorts`, `terminate_cohorts`, `split_cohorts`,
   `sort_cohorts`;
2. the annual patch block: `sort_patches`, `new_fuse_patches`, `terminate_patches`,
   `new_fuse_cohorts`, `terminate_cohorts`, `sort_cohorts`.

- **The same routines, not a copy.** One routine in `src/init`, which already links the
  demography operators (F10), calls them in that order. The two blocks honour
  `do_cohort_fissfuse` and `do_patch_fissfuse`, so a config with fusion off gets none at start
  either.
- **Census starts only.** Bare ground has nothing to fuse, and a state restart already carries a
  restructured stand.
- **What it keeps.** Stems and carbon are conserved by fusion. No census tree is culled: a single
  1 cm tree in a 400 m² patch holds 1.5e-4 kgC m⁻² of AGB, far above `min_cohort_agb = 1e-6`.
- **The fast reservoirs.** Fusion blends each patch's canopy-air and soil stores, which hold their
  default values at this point; `init_fast_reservoirs` then seeds every surviving patch. The Debug
  build (`-check all -fpe0`) must pass this path.
- **Reporting.** With `verbose`, the run prints the patch and cohort counts and the site's stems,
  AGB and LAI before and after, so G2 and G4 can be read off the log.
- Docs: the census paragraph of "How a run starts". CHANGELOG entry.

## 6. P1 — `scripts/prepare_census/make_census.py`

### 6.1 Inputs

- the census tree table, as CSV or `.rdata` (§4);
- `bci_census.toml`: the paths of the source files (default the local copy of §4) and their
  checksums, which the tool verifies before reading, the census year, the plot extent
  (1000 m × 500 m), `cell_size` (default 20 m) and the minimum diameter (10 mm);
- for the litter estimate and the summary only, the MEDS PFT file: the `[allometry]` coefficients,
  `wood_density`, `hgt_max`, `sla`, `aboveground_frac` and the turnover parameters. The census file
  itself needs none of them.

### 6.2 What the tool does

1. Keep status A with a `dbh` of at least 10 mm and coordinates inside the plot. Count every
   exclusion, by status, in the summary, including the live trees without a diameter.
2. Convert `dbh` from mm to cm, and record the multi-stem count from `nostems`.
3. Assign each tree to its square cell of side `cell_size`, numbering cells west to east and south
   to north. That is the order fusion visits them (F10). A grid that does not tile the plot, such
   as 40 m on the 500 m side, ends in a row of partial cells; each keeps its true area, here
   40 × 20 m.
4. Write one row per distinct (cell, `dbh`): `nplant` is the number of trees over the cell's area.

Nothing else: no ranking, no binning, no fusion. On the 2010 census at 20 m this gives 1,250
patches of 400 m² and 84,937 rows, a median of 68 per patch and at most 95, and the rows add back
to all 207,259 stems.

### 6.3 What MEDS's fusion makes of it

Simulated on the 2010 census with the defaults and MEDS's patch-fusion test: light profiles on 16
layers of 46 m, `patch_light_tol = 0.10` relaxed by 1.5 per pass, `max_patch = 12`
(`~/claude_workspace/meds_runs/bci_census/simulate_init_fusion.py`). The simulation takes patch
profiles from single trees, before cohort fusion, so the model's own result will differ in detail
(G4).

| Cell | Patches written | After each fusion pass | Final patches: area fraction and LAI, largest first |
|---|---|---|---|
| 20 m | 1,250 | 89, 26, 9 | 0.18/3.8, 0.18/5.0, 0.16/5.7, 0.14/4.6, 0.12/6.1, 0.10/9.0, 0.08/6.5, 0.04/7.1, 0.002/6.8 |
| 25 m | 800 | 88, 30, 9 | 0.33/4.6, 0.16/5.7, 0.16/7.7, 0.14/5.5, 0.10/5.9, 0.07/6.4, 0.05/3.6, 0.005/6.5, 0.001/7.8 |
| 40 m | 300 | 58, 14, 5 | 0.66/5.2, 0.30/6.5, 0.02/5.9, 0.01/4.6, 0.007/6.6 |
| 50 m | 200 | 25, 7 | 0.46/5.3, 0.43/6.0, 0.08/4.6, 0.015/5.5, 0.01/7.2, 0.005/4.8, 0.005/6.6 |

(The 40 m row was simulated with the partial strip dropped, before step 3 kept it.)

- From 20 m cells MEDS reaches `max_patch` in three passes and ends below it, at 9 patches with
  LAI 3.8 to 9.0. The areas are the model's own grouping, not a ranking the tool imposed.
- Larger cells average away more of the gap structure before fusion sees it: 50 m cells end with
  two patches holding 89 % of the area.
- Fusion stops at the first pass that reaches `max_patch`, so the third pass's larger tolerance
  overshoots, to 9 rather than 12. A larger `max_patch` keeps more patches.
- **Under the tolerance ceiling** (`patch_light_tol_max = 0.15`, merged
  as #322), the tolerance steps from 0.10 to 0.15 over six passes and stops.
  From 20 m cells the passes leave 89, 59, 45, 38, 29 and 24 patches, so the run keeps 24, above
  `max_patch`, with LAI from 2.9 to 9.0. From 50 m cells they reach 9 in the fourth pass.

### 6.4 Outputs

- `data/bci_census2010_meds.csv`, the file `init_from_census` reads (§5.3). A `#` header records
  the source file and checksum, the census dates, the filters, the cell size and the git commit.
  From the 2010 census at 20 m it begins:

  ```
  # MEDS census: Barro Colorado Island 50-ha plot, census 7 (2010-01-18 to 2011-03-28)
  # ...
  site_id,patch_id,patch_area,dbh,pft,nplant
  1,1,400.0,53.3,1,2.500000e-03
  1,1,400.0,37.2,1,2.500000e-03
  1,1,400.0,27.2,1,2.500000e-03
  ...
  1,1250,400.0,1.1,1,2.000000e-02
  1,1250,400.0,1.0,1,5.000000e-03
  ```

  A density of 2.5e-3 is one tree in 400 m²; 2.0e-2 is eight trees of 1.1 cm in the same cell.
- `data/bci_census2010_summary.json`: the plot's stems per ha, basal area, AGB and LAI under the PFT
  allometry, the raw table's own totals and `agb` for comparison (G1, G3), the exclusion counts,
  and the litter estimate (§6.5).

### 6.5 The litter estimate for the steady-state soil carbon

The steady state needs four litter inputs [kgC m⁻² day⁻¹]. The tool derives them from the census
stand and the PFT:
- leaf litter: leaf carbon over `leaf_lifespan_toc`;
- fine-root litter: fine-root carbon times `fineroot_turnover_rate`;
- wood: stand wood carbon times the plot's observed mortality rate by biomass from 2005 to 2010,
  measured at 1.90 % per year;
- the labile and structural split from `f_labile_leaf` and `f_labile_stem`.

The model's first-year `litter_*_site` outputs are then compared with this estimate (G7). The
climatological decomposition scalar `spinup_xi` stays at 1.0 in this attempt.

## 7. P2 — the example

| File | Change |
|---|---|
| `meds_config_spinup.toml` | deleted |
| `meds_config_eval.toml` | `init_mode = 1`, `census_file = "data/bci_census2010_meds.csv"`, the two `[init]` soil keys, and `[soil_carbon] spinup_steady = true` with the tool's litter input. `max_patch` stays 12 |
| `pft_parameters.toml` | nothing beyond #321, which gives it the corrected defaults |
| `bci_census.toml` | new: the census declaration (§6.1), with `cell_size = 20` |
| `run_example.py` | the spin-up step is replaced by building the census CSV. `--forcing-only` is unchanged |
| `plot_census.py` | new: the census stand against MEDS's after restructuring: stems per ha by size class, and each patch's area and LAI |
| `plot_evaluation.py` | unchanged. NEE becomes meaningful with a steady-state soil, so a NEE panel is worth adding, but it is optional |
| `README.md` | "The model runs" is rewritten around the census stand. "The PFT" reports the starting LAI and AGB against the census's own (F9). The data section cites both datasets. The evaluation table and figures are regenerated, and the runtime re-measured. The paragraph on the soil-water fix moves to the CHANGELOG entry it already has |

`data/` and `output/` are already ignored.

## 8. Validation gates

Stop and report at a failed gate. Do not tune past it.

| Gate | Check | Pass |
|---|---|---|
| G1 | conversion | the rows add back to every kept tree, and the patch areas to the plot's 50 ha. **Measured:** 84,937 rows, 207,259 stems, 1,250 patches of 400 m² at 20 m |
| G2 | MEDS reads the file as meant | before restructuring, MEDS's patch count, area fractions, stems, AGB and LAI equal the tool's summary to 1e-6 relative |
| G3 | plausibility of the initial stand | AGB against the census's own, and LAI against a measured BCI value the README cites. **Measured:** AGB 16.1 against 15.1 kgC m⁻², 7 % high; LAI 5.6 (F9). Report both; the allometry is MEDS's default by decision |
| G4 | restructuring at initialization | stems and carbon conserved to 1e-9 relative; at most `max_patch` patches and `max_cohort` cohorts per patch. **Predicted** (§6.3): 9 patches from 20 m cells. Report the model's patch areas and LAIs beside the simulation's |
| G5 | first cohort restructuring | the cohort count and site LAI change little across 2012-09-01, apart from growth and mortality, because the stand already went through the same operators |
| G6 | first patch restructuring | how many patches fuse on 2013-01-01. Few should, for the same reason |
| G7 | soil state | soil temperature does not drift by more than 1 K over the first 90 days. First-year litter is within 25 % of the §6.5 estimate |
| G8 | runtime | the restructuring at initialization and the five-year run are timed, and the README states both. More than 30 minutes on one core is a reason to revisit D4 |
| G9 | stand stability | the change in stand AGB and LAI over the first year, against BCI's observed −0.2 % per year for AGB (2005–2010). A loss above 10 % means the PFT's carbon balance cannot hold the census stand, and the README says so rather than hiding it |

## 9. P3 — tests

- **`scripts/prepare_census/tests/`** (pytest, registered in CTest as `prepare_census` beside
  `prepare_flux_tower`), on a synthetic ForestGEO table:
  - dead, prior, missing, diameter-less and sub-10 mm records are dropped;
  - mm convert to cm, and CSV and `.rdata` inputs give the same stand;
  - cells of 20, 25 and 50 m tile the plot; 40 m gives a row of partial cells with their true
    areas; the areas always add to the plot's;
  - trees with the same cell and diameter share one row, and the rows add back to every tree.
- **Fortran, the reader (`test/test_init_census.f90`):** columns in any order; the optional columns
  absent and present; `patch_area` normalized to fractions; a disagreeing area, a missing required
  column and an unknown column each refused; today's seven-column file giving the same site as
  before.
- **Fortran, restructuring at initialization:** a synthetic census of many small patches conserves
  stems and carbon, ends within `max_patch` and `max_cohort`, and culls nothing. It runs for a
  census start only, and not with `do_cohort_fissfuse` and `do_patch_fissfuse` off. Also run it in
  the Debug build.
- **Fortran, config:** the two `[init]` keys are read, defaulted and bounds-checked, and
  `init_fast_reservoirs` uses them. Extend the existing config test rather than adding a target.
- **Suite:** ifx Release and Debug, and gfortran. nvfortran is not installed on this machine
  (`CLAUDE.local.md`), so the nvfortran check waits for the GPU platform, and the PR says so.

## 10. Not in scope

- functional diversity: species to PFTs by wood density or other traits;
- changing any allometry relationship for this example. Chave et al. (2014)'s height model was
  measured as an alternative on 2026-09-28: with BCI's E = 0.0465 it gives the census AGB
  15.3 kgC m⁻², but it needs a (ln D)² term the MEDS height law lacks. The user chose the default
  heights;
- a comparison with a later census (decided 2026-09-28);
- patch ages from the census, which would need an optional `patch_age` column (a ROADMAP item);
- how well the 50-ha plot stands for the tower footprint (§4);
- lianas, and palm-specific allometry;
- buttress and height-of-measurement corrections to the diameter, and the other stems of multi-stem
  trees, which need the stem table;
- projecting census 7 forward from 2010 to the 2012 start.

## 11. Decisions taken (2026-09-28)

1. **Allometry:** MEDS's defaults, after the fix that makes the default biomass law Chave et al.
   (2014) in carbon and divides the leaf-area scale by `C2B` as ED2 does (D3, §5.2).
2. **No later-census check** in this example (D1).
3. **Cells:** 20 m quadrats to start, with `cell_size` configurable to any size (D4, §6.2).
4. **The preparation only maps census data to patches** (`patch_id`, `patch_area`, `nplant`), and
   MEDS does all fusion itself, before the model runs (D5, D6, §5.3–5.4).
5. **Patch fusion gets a tolerance ceiling**, `[demography].patch_light_tol_max = 0.15`, so
   `max_patch` is a target (2026-09-29, merged as #322). From 20 m cells the
   census start keeps 24 patches (§6.3).
6. **The biophysics example is retired** once this example passes its gates, and its #321
   regression is not fixed (2026-09-29).
