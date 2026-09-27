# MEDS polygon runtime plan — regional runs as an OpenMP loop over polygons, without MPI

> # 📐 DESIGN — written 2026-09-26, revised 2026-09-27. No code yet; R0–R2 are planned in detail (§10.1–§10.3).
>
> **What this plan does:** lets one MEDS process simulate a contiguous **region** of independent
> **polygons** (one polygon = one forcing grid cell with its own `site_t`), with the polygons
> advanced in parallel by an **OpenMP loop**. It uses no MPI. Large regions run in batches, or as
> HPC job arrays over spatial tiles. Scattered site networks are **not** run in one process: they
> run as separate processes (§8.1).
>
> **What it covers:** the non-MPI part of ROADMAP #183 ("a grid → polygon → site state hierarchy, an
> array of readers, a polygon loop and MPI").
>
> **What it depends on:** the forcing reader upgrade (F4, `MEDS_FORCING_DESIGN.md` §15), ✅ merged
> 2026-09-26 (#282): monthly block reads of the ED_ERA5land archive, and box selection in the library.
>
> **Revision 2026-09-27:**
> - the collection of polygons a run covers is a **region**, and it is always contiguous (§1.1);
> - site networks run as separate processes (§8.1), so the `polygons_csv` option is dropped;
> - a polygon's location lives in the polygon, not in the shared config or `site_t` (B12);
> - each output frequency needs a record queue, not one slot (B11);
> - output performance is analysed in §6.1;
> - R0–R2 are planned step by step in §10.1–§10.3.
>
> Baseline: `beta` at `b65c4b8`. The `file:line` citations are to that commit.

---

## 1. Goals and non-goals

**Goals:**

1. **Run the N polygons of a region in one process.** Each polygon has its own state (`site_t`, the
   fast context, budgets, the ledger) and its own forcing cell. All share one configuration and one
   forcing reader.
2. **Parallelise with OpenMP over polygons.** Results must be **byte-identical at any thread
   count**, as the patch loop already guarantees within one site (`test_fast_loop.f90:246-253`).
3. **Keep every single-site run working and unchanged:** a site run goes through the same code as a
   region of one polygon.
4. **Scale to continents and the globe** by batching polygons and by running spatial tiles as
   independent jobs, not by MPI.

**Non-goals:**

- **MPI.** Tiles map naturally onto ranks later if ever needed.
- **Scattered site networks in one process.** They share no forcing reads, so running them together
  saves nothing; separate processes are simpler and isolate failures (§8.1).
- **Interactions between polygons** (seed dispersal, fire spread, lateral water). ED2's seed
  exchange is among the patches *within* a polygon (`ed2_comparison.md:126`), which MEDS already
  handles inside `site_t`.
- **Several sites per polygon, and spatially varying parameters and initial conditions.** Soil
  texture, PFT maps and census data are a later plan (§12, OR2). Here every polygon holds one site,
  shares the parameters and the initialisation recipe, and differs only in forcing and location.
- **GPU offload.**

### 1.1 Vocabulary

| Term | Meaning in MEDS | ED2 counterpart |
|---|---|---|
| **site** | the ecological state of one location: soil column, patches, cohorts (`site_t`) | site (`sitetype`): a soil or terrain subdivision of a polygon |
| **polygon** | one location with one meteorological forcing (an ED_ERA5land cell) and its own state: location, forcing cursor, output accumulators, budgets, **one** `site_t` (`meds_polygon_t`, §10.3) | polygon (`polygontype`), which may hold several sites |
| **region** | the polygons of one run: always **contiguous**; a lat/lon box now, a mask later (`meds_region_t`, §10.3) | a regional grid (`edtype`, one per `N_ED_REGION` box) |
| **patch**, **cohort** | as today | patch, cohort |

- **Why ED2's order holds.** In ED2 a *grid* is a mesh (the atmospheric model's, with optional
  nests) and a polygon is one cell of it: `edtype` "contains arrays of polygons that populate the
  current grid" (ED2 `ed_state_vars.F90:2591`). MEDS keeps that order. It calls the collection a
  region because it has one mesh (the forcing archive) and no nests.
- **A region is not a state level.** No process acts at region scale (interactions between polygons
  are a non-goal), so a region is a *selection* of polygons plus the unit of I/O. Regional means over
  sub-areas (ecoregions, watersheds) are an output aggregation, not a hierarchy level.
- **One polygon, one forcing.** An area spanning several forcing cells is several polygons.

## 2. What the code already supports

From the 2026-09-26 survey, refreshed at `b65c4b8`:

- **One self-contained run object.**
  - `meds_run_t` (`src/main/meds_driver.f90:79-105`) holds the config, one `site_t`, the output
    manager, the fast context, one `met_driver_t`, the budgets, the slow ledger, the clocks and the
    counters.
  - It is a local of the program, and the code says so: "no state hides in module scope, so two runs
    can be open at once" (`meds_driver.f90:75-78`).
  - The main loop is thin: `driver_open` → `driver_step` loop → `driver_finalize`
    (`meds_main.f90:57-68`). `driver_step` advances one slow step (one day by default).
- **State is passed as arguments.** `site_t` is a flat structure of arrays
  (`src/state/site/meds_site_state_types.f90:342`). There is no RNG, and every file unit uses
  `newunit`.
- **The patch loop is already OpenMP-parallel and deterministic:**
  - the directive is at `meds_fast_dynamics.f90:514-818`;
  - per-thread scratch pools are *local* allocatables sized to the thread count, allocated on each
    call, so each polygon thread naturally gets its own;
  - results are folded back in patch order.
- **The met samples of a slow step are taken outside the patch loop** (`meds_fast_dynamics.f90:471-485`):
  `met_advance` and `met_instant` run once per sub-step, in increasing time, before the patches.
- **The build already turns on recursive locals** (`MEDS_AUTO_LOCALS` in `CMakeLists.txt`), which is
  mandatory when threading and removes races on static locals.
- **The forcing side is ready (F4).** The ED_ERA5land reader loads a month of any list of cells with
  one chunk-column read per 16 × 16 block, and selects a box's valid cells in row-major order,
  across 180° (`meds_era5land_reader`: `era5land_select_box`, `era5land_load_month`).

## 3. What blocks a parallel polygon loop today, and how each is handled

| # | Blocker (citation at `b65c4b8`) | Resolution |
|---|---|---|
| B1 | **netCDF inside the time step.** `driver_step` flushes output twice per day (`output_serialize_pending`, `meds_driver.f90:367` for the fast tier and `:379` for the others) and writes yearly checkpoints (`:420-424`). Forcing is read inside `met_advance` (`meds_fast_dynamics.f90:482`): the archive backend loads a month on demand (`ensure_month`, `meds_met_driver.f90:707`), and the single-file backend reads every value with a 1 × 1 hyperslab (`read_scalar`, `:860`). HDF5 serialises every call under a global lock, and netCDF-C is not thread-safe (`MEDS_IO_DESIGN.md` §5.5). | Split each step into a **compute phase** (no netCDF at all) and an **I/O phase** that runs serially at month boundaries (§4, R1). Forcing is loaded in the I/O phase; the single-file backend holds its cell's whole series in memory. |
| B2 | **Module-level allometry.** Nine `protected` reals (`src/shared/functions/meds_allometry.f90:50-58`) are written by `set_allometry` from every config load (`meds_config.f90:473`). | Load config **once, serially, before** the parallel loop. All polygons share one parameter set (non-goal: per-polygon parameters). A later refactor moves these into the config type (OR2). |
| B3 | **Saved state in the fast-loop probe.** `write_fast_probe` keeps `save` variables (`meds_fast_dynamics.f90:919-920`). | Reject `fast_probe` in region mode, as `n_threads > 1` already does (`meds_config.f90:691`). |
| B4 | **Nested OpenMP.** A polygon `parallel do` around the step nests the patch-loop region. | In region mode force `cfg%n_threads = 1` inside polygons (a `validate_config` rule), and set `omp_set_max_active_levels(1)`. The NVHPC `target` loop (`meds_demography_update.f90:128`) must also stay single-level; region mode requires `MEDS_GPU = none` for now. |
| B5 | **Shared file names.** Output streams `<dir>/<prefix>-<tier>…nc` (`meds_output_stream.f90:176-178`), restarts `<prefix>-S-…nc` and `_pft_parameters.csv` (`meds_driver.f90:254`) would collide between polygons. | Region mode writes **one file per tier with a `polygon` dimension**, and one ragged restart file per checkpoint (§6). There are no per-polygon files except for `detail_polygons`. `ensure_output_dir` (`meds_driver.f90:603`) runs once, serially. |
| B6 | **`error stop` kills the process.** There are 12 sites in the fast kernels, plus `nc_check`, the reader and census parsing. One bad polygon ends the whole run. | R3 fails fast, with the polygon id in the message. R5 converts the common failures to a status return, marks the polygon failed, carries on, and reports the failed polygons at the end (§7). |
| B7 | **Interleaved stdout.** `print_summary`, budget reports and reader or output lines would interleave between threads. | In region mode per-polygon printing is off (`verbose` forced low). A per-month progress line and an end-of-run per-polygon status table are written in the serial phase. |
| B8 | **Allocator contention.** About 26 `allocate` calls per patch per `dt_fast` in `build_column_frozen` take about 24% of allocator self-time (#195); more threads mean more heap contention. | The R5 performance item: fix #195 by pre-allocation, or link a thread-scalable allocator (jemalloc or tcmalloc via `LD_PRELOAD`), whichever measures better. R0 records the baseline. |
| B9 | **`BLOCK` constructs.** nvfortran rejects `BLOCK` inside a parallel region (`meds_fast_dynamics.f90:678-680`), and `driver_step` has one (`meds_driver.f90:401`). | The rule concerns constructs lexically inside the region, and `driver_step` is *called* from it. **R0 checks this with nvfortran** where it is installed (it is not on the development cluster). If nvfortran rejects it, R1 hoists the `BLOCK`. |
| B10 | **The ifx `private` derived-type mold.** Compiler-generated static molds are shared between threads (`docs/building.md`). | The polygon loop declares **no** `private` or `firstprivate` derived types. Each iteration works on its own element `poly(p)` of a shared array, and scratch lives in called routines. |
| B11 | **One pending record per output frequency.** A closed period is staged in `mgr%pending(t)` (`meds_output_integrate.f90:802`) and must be written before the next one closes. With writes deferred to month end, the daily tier closes about 30 records in between, and the fast tier up to 744. | Each frequency gets a **record queue**, filled in the compute phase and drained in the I/O phase (R1). |
| B12 | **Site location in the shared config.** Latitude, longitude and elevation live in `cfg%forcing`, and the model reads latitude outside the reader, for the day length that drives phenology (`meds_vegetation_dynamics.f90:1018-1019`). All polygons share one config. | Location belongs to the **polygon**, as in ED2. The reader's per-polygon part carries it for solar geometry, and latitude is passed down to the day-length call (R2). `site_t` stays ecological state only, so it can later hold several sites per polygon. |

## 4. Execution model: month-synchronous, compute and I/O phases

```
serial:   load config once (B2) → build the region: polygons from the box + selection rules (§5)
          allocate poly(1:N): site_t, fast_ctx, budgets, ledger, output part, forcing cursor, location
          initialise every polygon (bare ground in R2; restart in R4), no I/O inside
loop over months:
  serial:   forcing reader loads the month block for ALL polygons, plus the one record before it
  parallel: !$omp parallel do schedule(dynamic, chunk)                       (R3; a plain loop in R2)
              do p = 1, N:  advance poly(p) through the whole month (fast + slow steps),
                            forcing from the shared read-only month buffer at poly(p)%cell,
                            output records queued in poly(p)'s output part
  serial:   pack and write every polygon's queued records to the region files (§6)
            yearly: write the checkpoint (ragged restart, R4); flush logs and status
```

- **Why a month is the synchronisation unit:**
  - the forcing archive is stored by month (forcing design FD7);
  - a month of one polygon is long enough that the per-region cost of OpenMP is negligible;
  - `schedule(dynamic)` balances polygons that cost different amounts, for example many-cohort
    tropical forest against tundra.
- **Why compute and I/O are split** (the reasoning behind B1, B5 and B11):
  - netCDF cannot be called from several threads at once, and a region file is shared by all
    polygons, so all file work is gathered where one thread does it;
  - bulk I/O is far faster than piecemeal I/O: a month of all cells in one chunk-column pass instead
    of thousands of tiny reads, and one hyperslab per variable per record instead of many small
    writes;
  - record order is the polygon order, and every polygon is at the same model time at a write or a
    checkpoint, so results cannot depend on thread timing;
  - with no netCDF in the compute phase, a writer thread can later overlap the I/O phase with the
    next month (§6.1).
- **A site run is a region of one polygon** through the same code, with I/O between months instead
  of every step. Output content must be byte-identical to today's; R1 checks this.
- **Costs accepted.** Output records wait in memory until the month ends: about 100 MB at worst for
  one site's fast tier with per-cohort variables. A crash loses at most the current month's
  unwritten output, and the yearly checkpoint still bounds the rerun.

## 5. Region and polygon selection

- **The region** is a contiguous set of forcing cells: a lat/lon box (`[region].box_nwse`, crossing
  0° or 180° correctly) now, and a mask or ID raster later. Each valid forcing cell inside it is a
  candidate polygon, in row-major order (`era5land_select_box`).
- **Selection rules** (config), since not every valid forcing cell should carry vegetation:
  - **`land_fraction_min`.** The valid mask includes large lakes, for example the whole Caspian
    (6,866 valid cells have `lsm = 0`), and coastal fractional cells. The default is 0.5, read from
    the static file's `land_fraction`.
  - **`exclude_ice`** (later). The Antarctic and Greenland ice sheets are about 0.74 M of the 2.2 M
    valid cells. The rule needs a glacier mask in the static file (forcing design §14.3), which
    ERA5-Land does not supply directly; until then, boxes simply avoid the ice sheets.
- **Per-polygon attributes** come from the static file: latitude and longitude (the cell centre),
  `grid_elevation`, and the polygon id. `utc_offset` is 0, because the archive is UTC.
- **Polygon id** (OR6): the cell's global row-major index on the forcing grid, `row × nlon + col`.
  An id is then stable across regions and tiles, and outputs from different jobs line up.

## 6. Output and restart in region mode

- **Output streams:** one file per tier per time chunk for the whole region, with a leading
  **`polygon`** dimension. `MEDS_IO_DESIGN.md` §10 Q5 already anticipates a leading-axis
  convention.
  - **Coordinates:** `polygon_id(polygon)`, `lat(polygon)`, `lon(polygon)`, and the grid indices
    `row(polygon)` and `col(polygon)`, so a region maps back onto the lat/lon grid in one step (OR7).
  - **Variables:** only **fixed-shape** ones: site totals `(time, polygon)`, per-PFT
    `(time, polygon, pft)`, per-size-class `(time, polygon, dbh_class)`, and per-soil-layer
    `(time, polygon, soil)`.
  - **Cohort-level and patch-level variables** (ragged across polygons) are not written for the
    region. Full detail goes only to a short `detail_polygons` list, as ordinary single-site files
    named by polygon id.
  - **Fast-tier (`F`) output** is limited to `detail_polygons`. Daily, monthly and yearly tiers are
    written for all polygons.
  - **Writing** happens in the serial phase, one hyperslab per variable per record over all
    polygons, into the configured output folder. On a cluster that folder should be on node-local
    disk, with each file copied to shared storage when its chunk closes. HDF5's many small writes are
    slow over NFS: builds of the forcing archive ran 3–4× slower when several wrote over NFS at once.
- **Restart (R4):** one netCDF file per checkpoint for the whole region, using CF **contiguous ragged
  arrays**:
  - per-polygon counts of patches, and per-patch counts of cohorts, plus offsets;
  - the site, patch and cohort fields concatenated.

  Restoring a region restores each polygon exactly. A single-site restart stays in today's format.
- **Scale check** for 305,000 North American polygons: a monthly record of 50 site-level variables is
  305k × 50 × 4 B, about 61 MB, which is trivial. Restart size depends on cohort counts; R0 measures
  `site_t` size after spin-up.

### 6.1 Output performance

Writing and computing both grow with the number of polygons, so the region size cancels out. The
share of wall time spent in the serial write is about

> (bytes written per polygon-month ÷ write speed) × cores ÷ compute time per polygon-month.

With about 100 MB/s of compressed writing, 40 cores and about 1 s of compute per polygon-month (the
real figure is an R0 measurement):

| Output for all polygons | Bytes per polygon-month | Write share |
|---|---|---|
| Monthly, about 50 site variables | about 200 B | about 0.01% |
| Daily, with PFT and soil variables | about 25 KB | about 1% |
| Fast tier (hourly), about 200 variables | about 600 KB | about 25% |

The last row is why the fast tier is limited to `detail_polygons`. The yearly restart (roughly 10 GB
for 300,000 established polygons) is written once a year and does not register.

Ways to use more cores for output, cheapest first:

1. **Pack in the parallel phase.** Each polygon copies its closed records into its own slice of
   region-wide record arrays. No netCDF is involved, so it is safe. The serial phase then issues one
   write per variable per record. Done from R2.
2. **Overlap the I/O phase with the next month.** Load month m+1's forcing, start computing it on
   N−1 cores, and let one thread write month m. This is safe because only that thread calls netCDF
   while the others compute; it needs double record buffers. Added only if R3 measures a write share
   above about 5% for the output actually requested.
3. **Offload to other processes.** Write fast, uncompressed, write-optimised files during the run,
   and let a follow-on Slurm job compress and rechunk them for time-series reading, in parallel across
   variables and years.
4. **Tiles as a job array** (§8). Every job writes its own files at the same time, and readers open
   them together along `polygon`.

Not used: netCDF's own asynchronous or parallel I/O (not available on the HDF5 path without MPI and
PnetCDF, `MEDS_IO_DESIGN.md` §5.5), or several threads writing netCDF at once.

## 7. Failure handling, logging, determinism

- **Determinism.**
  - Polygons are independent, so each polygon's trajectory depends only on its inputs, never on the
    thread count or the schedule.
  - Output order is the polygon order, fixed by the region, not by thread completion.
  - A CTest asserts byte identity between 1 and 4 polygon threads (R3).
- **Failures:**
  - **R3 fails fast:** the first failure stops the run, naming the polygon id and its lat/lon.
  - **R5 isolates failures:**
    - the step path returns a status for the common failures: the NaN check, conservation-ledger
      failures and forcing gaps;
    - a failed polygon is frozen and flagged `failed` in the outputs, with the step and reason, and
      the others continue;
    - `error stop` stays only for programming errors.
- **Logging:** each polygon keeps a small in-memory log (a ring buffer). The serial phase writes the
  logs of failed or flagged polygons, plus a per-month progress line (polygons done, failed,
  wall time).

## 8. Memory, batching and very large regions

- **Memory per polygon** is `site_t` plus the fast context, the output part and the forcing cursor.
  Cohort counts drive it; **R0 measures it** after spin-up on representative sites.
- **Batching:** if N × memory per polygon exceeds `memory_limit_gb`, the region runs in **batches**
  of `batch_size` polygons. Each batch runs its full simulation period (every month, with I/O as in
  §4) before the next starts.
  - Batches are built from **spatially compact tiles**, aligned to the forcing chunks (16 × 16
    cells), so each batch's monthly forcing read touches few chunks.
  - The outputs of all batches go into the same region files, as disjoint polygon ranges.
- **Job arrays:** because polygons are independent, a continental or global region can also be split
  into tiles submitted as an HPC **job array**. Each job runs one tile in region mode and writes its
  own files, and readers open the tile files together along `polygon` (polygon ids are stable, §5).
  This scales like MPI with none of its machinery.
- **The forcing buffer per batch** is 8 variables × cells × 744 h × 4 B: about 0.5 GB for 20,000
  cells, and about 7.3 GB for all North American land.
- **Loading it** with chunk-column reads took **5.0 s per variable-month** for all of North America,
  so about 40 s for all 8 variables in the serial phase. That is small next to a month of compute for
  300,000 polygons.

### 8.1 Site networks: separate processes

A network of scattered sites (flux towers, plots) runs as **separate `meds_main` processes**, one
config per site generated from a template (latitude, longitude, output prefix). Scattered sites share
no forcing chunks, so one process would save nothing, and separate processes isolate failures. R6
documents two recipes:

- **A job array:** one site per task (`sbatch --array=1-N%K`), as the archive build does.
- **One allocation filled with a task farm,** for schedulers that allocate whole nodes or cap jobs per
  user:
  ```bash
  #SBATCH --nodes=1 --ntasks=1 --cpus-per-task=40
  export OMP_NUM_THREADS=1
  parallel -j $SLURM_CPUS_PER_TASK --joblog sites.joblog \
      'meds_main configs/{}.toml > logs/{}.log 2>&1' :::: sites.txt
  ```
  GNU parallel starts the next site as one finishes, so uneven sites balance, and `--resume-failed`
  reruns the failures. Across nodes, each run is an `srun --exclusive -N1 -n1 -c1` job step.

A small merge script assembles a network-wide file from the per-site outputs when one is wanted.
Many copies of one site with different parameters (calibration ensembles) are a different case: they
need per-polygon parameters (OR2) and would then run in a region-like container.

## 9. Configuration sketch

```toml
[run]
mode = "region"                   # "site" (default; today's behaviour) | "region"

[forcing]                         # as implemented in F4 (MEDS_FORCING_DESIGN.md §15)
format          = "era5land"
data_path       = "<installation-specific path to the ED_ERA5land archive>"

[region]
box_nwse          = [45.1, -79.8, 40.4, -71.8]   # contiguous box, may cross 0 or 180 degrees
land_fraction_min = 0.5           # skip lakes and fractional coastal cells
detail_polygons   = []            # polygon ids that also get full cohort/patch and fast-tier output
threads           = 32            # R3: OpenMP threads over polygons; patch threading forced to 1 (B4)
schedule_chunk    = 4             # R3
memory_limit_gb   = 64            # R5: beyond this, the region runs in batches of compact tiles
```

`validate_config` enforces the following in region mode:
- `[forcing].format = "era5land"`, because a single MEDS forcing file holds only a few locations;
- `n_threads = 1` (patch threads) and no `fast_probe`;
- `MEDS_GPU = none`;
- until R4: bare-ground initialisation only, and no checkpoints (`[state].write_state = false`);
- the `[site]` location keys and `[forcing].max_distance_km` are rejected, since each polygon takes
  its location from its cell.

## 10. Phases

| Phase | Content | Acceptance |
|---|---|---|
| **R0** measure and verify | Cost of a simulated year and of a polygon-month from a spun-up forest; `site_t` memory; the allocator profile (#195); nvfortran on a `BLOCK` in a routine called from a parallel region (B9) | Numbers recorded in §10.1. B9 answered where nvfortran exists. |
| **R1** compute/I-O split | The step split into a compute phase (no netCDF) and a month-boundary I/O phase; per-frequency record queues (B11); forcing loaded only in the I/O phase; a site run otherwise unchanged (§10.2) | The CTest suite green; single-site outputs **byte-identical** to before. |
| **R2** region and polygon container, serial | The reader split into a shared source and per-polygon cursors; location in the polygon (B12); the output manager split into shared and per-polygon parts; `meds_region_t` and `meds_polygon_t`; the month-synchronous loop without OpenMP; region-dimension output; `detail_polygons` (§10.3) | N polygons run as one region produce outputs **byte-identical** to N separate single-site runs (a 4-polygon synthetic test). |
| **R3** OpenMP polygon loop | `!$omp parallel do schedule(dynamic)` over polygons; the B3, B4 and B7 rules in `validate_config`; fail-fast messages with the polygon id; the write share measured (§6.1) | Byte identity between 1 and 4 threads; scaling measured to the core count; nvfortran build green where available. |
| **R4** ragged restart | Region checkpoint and restart with CF contiguous ragged arrays | A restart round trip is bit-identical to an uninterrupted run. |
| **R5** robustness and scale | Status-based failure isolation (B6), batching by tiles, the tile job-array recipe, allocator work (B8), per-polygon logs | A failure-injection test: one polygon fails and the rest match the reference. Throughput and memory recorded for a 20,000-polygon box. |
| **R6** interfaces and docs | Optional C API and Python entry point for region runs (the C API registry holds at most 4 runs and is unsynchronised, `meds_c_api_run.f90:46-47`); docs, including the site-network recipes (§8.1); an example | A documented example runs a small box end to end. |

**Order:** R0 anytime (it opens the R1 work). R1 needs F4, which is done. R2 needs R1, and R3 needs
R2. R4 and R5 follow R3.

### 10.1 R0 — measure and verify (about half a day)

1. **Cost of a simulated year and of a polygon-month.** Run the `example_biophysics` spin-up to its
   restart state, then time the July stage and a full year from that state, from the archive
   (`format = "era5land"`), in Release, on one thread. Record wall time per polygon-month for a
   bare-ground start and for the established stand.
2. **`site_t` memory after spin-up.** Patch and cohort counts, and bytes per polygon for `site_t`,
   the fast context and the output accumulators. These set the R2 memory budget and the R5 batch
   size.
3. **Allocator profile (#195).** The share of time in `allocate` and `deallocate` per sub-step: the
   baseline for B8.
4. **B9 on nvfortran**, where it is installed: a routine containing `BLOCK`, called from inside a
   `parallel do`. This is not possible on the development cluster, which has no nvfortran; it is
   flagged, not skipped.
5. **The inputs to §6.1:** bytes per polygon-month for the daily and monthly tiers of the default
   output set.

Deliverable: the numbers written into this section; no code change.

### 10.2 R1 — compute/I-O split (single site; outputs unchanged)

**Goal:** a step makes no netCDF call, and all file work happens in an I/O phase at month boundaries
and at the end of the run. Single-site outputs are byte-identical to before.

1. **Split the step** (`meds_driver.f90`).
   - `driver_step` becomes `driver_compute_step` (today's body, lines 318-418, without the two
     flushes and the checkpoint) and `driver_io_phase` (drain the record queues, write the yearly
     checkpoint, load the next month's forcing).
   - The I/O phase runs when the step just taken closes a calendar month, and in `driver_finalize`
     before the streams close. Yearly checkpoints (`:420-424`) already fall on a month boundary.
   - `driver_step` keeps its signature (compute, then the I/O phase if a month closed), so the C API
     (`meds_c_api_run.f90:99`) and Python do not change.
2. **Record queues for output** (B11; `meds_output_types.f90`, `meds_output_integrate.f90`,
   `meds_output_manager.f90`).
   - `pending(t)` becomes a queue of closed records per frequency, grown on demand. `close_tier`
     appends, and `output_serialize_pending` drains the queue in order.
   - The fast tier still closes every `fast_interval_steps` sub-steps (`meds_driver.f90:361-371`),
     but no longer writes there; its records wait in the queue.
   - Memory: a month of one site's fast-tier records, sized by the live cohort count rather than
     `cohort_max`.
3. **Forcing in memory** (`meds_met_driver.f90`).
   - **Single file:** `met_open` reads the chosen cell's whole series once (at most a few years of
     hourly values), and `read_record` takes values from memory instead of `read_scalar`.
   - **Archive:** a new `met_prefetch(drv, year, month)`, called from the I/O phase and from
     `met_open`, loads the archive month that the coming model month maps to (through the recycle
     mapping when recycling). It also keeps the one record before that month: 00:00 on the 1st from
     the previous month's file, or the window's last record at the recycle wrap.
   - `ensure_month` loses its load path, and a record that is not in memory becomes a programming
     error (`error stop`). The reader counts its loads (`drv%n_loads`), so a test can assert that
     none happen inside the compute phase.
4. **Nothing else touches netCDF in the step.** The checkpoint moves into the I/O phase; the census
   and restart reads stay in `driver_open`; `fast_probe` (formatted text, not netCDF) is unchanged.
5. **Tests.**
   - **Queued against immediate writes:** a unit test integrates the same records with a drain after
     each one and with one drain at the end, and requires identical files.
   - **Month seam without loads:** the synthetic archive of `test_met_era5land`, across a month
     boundary and the recycle wrap, asserting that `n_loads` changes only in `met_prefetch`.
   - **Byte identity on real runs:** a script run before and after R1 compares every variable and
     attribute of the outputs of the July 2024 Ithaca runs (archive and single file, as in the F4
     acceptance run) and of `example_demography`. The result goes in the PR.
6. **Acceptance:** the CTest suite green (ifx Release, and Debug for the touched tests); the
   byte-identical outputs of step 5; nvfortran where available.
7. **Size:** 2–3 working sessions, one PR.

### 10.3 R2 — region and polygon container, serial

**Goal:** the N polygons of a region run in one process, one after another, each identical to its
own single-site run. Three PRs, each keeping single-site output unchanged.

**PR 1 — the reader split and the polygon's location.**

1. **Split `met_driver_t`** (`meds_forcing_types.f90`, `meds_met_driver.f90`).
   - `met_source_t`, shared: the forcing config, the time axis and recycle window, the cell list
     (`met_domain_t` renamed `met_cells_t`), and the month buffer with its carried record.
   - `met_cursor_t`, per polygon: the cell index, latitude, longitude and grid elevation, the two
     bracketing records, the cursor, the seam flag, and the last daytime clearness index.
   - `met_advance(src, cur, t)` and `met_instant(src, cur, t)`; `met_prefetch(src, ...)` loads every
     cell. A site run is one source and one cursor.
   - `fast_dynamics` and `advance_one_step` take the source and the cursor instead of `met_drv`
     (`meds_fast_dynamics.f90:271-277`, `meds_stepper.f90:36-43`).
2. **Location in the polygon** (B12). The cursor carries latitude and longitude for solar geometry,
   and latitude is passed through `advance_one_step` to the day-length call
   (`meds_vegetation_dynamics.f90:1018-1019`). In site mode the values come from `[site]`, exactly
   as today.
3. **Tests:** `test_met_driver` and `test_met_era5land` on the new interfaces, including two cursors
   on one source reading two cells correctly; the R1 byte-identity script.

**PR 2 — the output manager split.**

4. **Split `output_manager_t`** (`meds_output_types.f90`, `meds_output_registry.f90`,
   `meds_output_integrate.f90`, `meds_output_manager.f90`, `meds_output_stream.f90`).
   - Shared (`output_shared_t`): the registry, the diagnostic parameters, the file settings, the
     streams and the provenance attributes.
   - Per polygon (`output_part_t`): the integration buffers, `t_open` and `has_data`, the record
     queues and the fast-tier staging.
   - A site run is one shared part and one polygon part, and writes today's files.
5. **Tests:** the output unit tests on the new types; the R1 byte-identity script. This is the
   largest and riskiest PR of R2.

**PR 3 — the region driver, region files and the equivalence test.**

6. **Config** (`meds_config.f90`, `meds_config_io.f90`, and a new `meds_region_opts` leaf in
   `src/config/`): `[run].mode = "site" | "region"`; the `[region]` block (`box_nwse`,
   `land_fraction_min`, `detail_polygons`); the region-mode rules of §9.
7. **Types and driver** (a new `src/main/meds_region.f90`).
   - `meds_polygon_t`: id, cell, location, `site_t`, fast context, budgets, slow ledger, seam
     statistics, `met_cursor_t`, `output_part_t` and status.
   - `meds_region_t`: the config, `met_source_t`, `output_shared_t` and `poly(1:N)`.
   - `region_open`: select the box's valid cells (`era5land_select_box`, extended to return the
     static `land_fraction`), apply `land_fraction_min`, build the polygons, and initialise each from
     bare ground.
   - The month loop of §4 without OpenMP: prefetch, advance each polygon through the month with the
     same compute step as a site run, then the I/O phase.
   - `region_finalize`: the end-of-run per-polygon status table (budget verdicts, conservation) and
     closing the files.
   - `meds_main` dispatches on `[run].mode`.
8. **Region files** (a region stream beside the site stream in `meds_output_stream.f90`).
   - Dimensions `(time, polygon[, pft | dbh_class | soil])`; coordinates `polygon_id`, `lat`, `lon`,
     `row` and `col`; only fixed-shape variables (§6, OR1).
   - Packing: in the I/O phase, each record's per-polygon values are gathered into region-wide
     arrays (§6.1 item 1), then written with one hyperslab per variable per record.
   - `detail_polygons` write ordinary single-site files named `<prefix>-p<polygon_id>-…`.
9. **Tests.**
   - **Serial equivalence:** a 4-polygon region on the synthetic archive against 4 single-site runs at
     the same cells. Every region variable's `polygon = p` slice must equal the site run's series, and
     a detail polygon's files must equal the site run's files.
   - **Selection:** a box across 180°, `land_fraction_min`, and the polygon id convention.
   - **Config rejections:** the region-mode rules of §9.
10. **Demo:** a 1° box around Ithaca (about 100 polygons) for one year from the archive, run
    serially, with wall time and memory recorded here as the baseline for R3.
11. **Size:** 4–6 working sessions across the three PRs.

## 11. Tests (CTest)

- **Queued writes:** drained at once or per record, identical files (R1).
- **Month seam:** forcing continuity across a month boundary and the recycle wrap, with no loads in
  the compute phase (R1).
- **Serial equivalence:** 4 synthetic polygons as a region produce byte-identical outputs to 4
  single-site runs (R2).
- **Thread invariance:** 1 against 4 polygon threads, byte-identical (R3).
- **Restart round trip** of a region with different cohort counts per polygon (R4).
- **Failure isolation:** inject a NaN into one polygon's forcing; that polygon is flagged and the
  others are unchanged (R5).

## 12. Open decisions

| # | Decision | Recommendation |
|---|---|---|
| OR1 | The variable set for region output | Start from the site-level, per-PFT and per-size-class variables enabled in `meds_io_config.toml`, as a `region_output` group. |
| OR2 | Spatially varying inputs (soil, PFT composition, parameters), several sites per polygon, calibration ensembles | A separate plan after R3. It needs the allometry refactor (B2), which moves parameters from module state into the config. |
| OR3 | Failure policy default: abort or skip-and-flag | Abort during development (R3). Skip-and-flag is available from R5, selectable by config. |
| OR4 | Batching inside one process, or job arrays across tiles | Both: batching for a workstation or a single node, job arrays for continental and global runs on HPC. |
| OR5 | Whether to simulate lake and ice cells at all | No by default (`land_fraction_min = 0.5`, and boxes that avoid the ice sheets). Revisit with a lake or glacier surface scheme. |
| OR6 | Polygon id | The cell's global row-major index on the forcing grid (§5), stable across regions and tiles. |
| OR7 | Region output layout | A 1-D `polygon` dimension with `row` and `col` indices (§6), rather than a gridded `(time, lat, lon)` layout; a gridded view is a one-line reindex. |
| OR8 | Region mode before R4 | Bare-ground starts only, and no checkpoints until the ragged restart exists; census starts need spatial inputs (OR2). |

**Decided 2026-09-27:** the vocabulary of §1.1 (region, always contiguous; polygon; site), and site
networks as separate processes (§8.1) instead of a `polygons_csv` option.

## 13. Relationship to other documents

- **ROADMAP #183** (multi-polygon runtime with MPI): this plan delivers everything except MPI. The
  ROADMAP entry is updated when R3 ships.
- **`MEDS_FORCING_DESIGN.md`:** its §3.2 and §8 P2 ("one `met_driver_t` per polygon") are superseded
  by one shared reader with per-polygon cursors (R2). The `domain`/`box_nwse` keys of its §15.2
  arrive here as `[run].mode = "region"` and `[region].box_nwse` (§9). Its §10 Q7 (a polygon with no
  covering land cell) is answered by its §15.5 and the selection rules here.
- **`MEDS_GPU_EVALUATION.md`** §10 ("MEDS has no site or ensemble axis at all"): the polygon axis is
  exactly that axis. Its §7 and §12.4 (allocator traffic, #195) become more pressing with polygon
  threads (B8).
- **`MEDS_IO_DESIGN.md`** §5.5 (the HDF5 global lock, the shelved asynchronous writer) and §10 Q5 (the
  leading-axis convention) are resolved by the serial I/O phase, §6.1 and the `polygon` dimension.
