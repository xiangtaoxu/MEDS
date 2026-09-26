# MEDS polygon runtime plan — regional runs as an OpenMP loop over polygons, without MPI

> # 📐 DESIGN — written 2026-09-26. No code yet.
>
> **What this plan does:** lets one MEDS process simulate many independent **polygons** (one polygon
> = one forcing grid cell with its own `site_t`), with the polygons advanced in parallel by an
> **OpenMP loop**. It uses no MPI. Large domains run in batches, or as HPC job arrays over spatial
> tiles.
>
> **What it covers:** the non-MPI part of ROADMAP #183 ("a grid → polygon → site state hierarchy, an
> array of readers, a polygon loop and MPI").
>
> **What it depends on:** the forcing reader upgrade (F4) in `MEDS_FORCING_DESIGN.md` §15, for
> `domain = "box"` and its monthly block reads.
>
> Baseline: `main` = `beta` at `92cad44`. The `file:line` citations are to that commit, from a code
> survey on 2026-09-26.

---

## 1. Goals and non-goals

**Goals:**

1. **Run N polygons in one process.** Each polygon has its own state (`site_t`, the fast context,
   budgets, the ledger) and its own forcing cell. All share one configuration and one forcing
   reader.
2. **Parallelise with OpenMP over polygons.** Results must be **byte-identical at any thread
   count**, as the patch loop already guarantees within one site (`test_fast_loop.f90:249-277`).
3. **Keep every single-site run working and unchanged:** today's mode is simply N = 1.
4. **Scale to continents and the globe** by batching polygons and by running spatial tiles as
   independent jobs, not by MPI.

**Non-goals:**

- **MPI.** Tiles map naturally onto ranks later if ever needed.
- **Interactions between polygons** (seed dispersal, fire spread, lateral water). ED2's seed
  exchange is among the patches *within* a polygon (`ed2_comparison.md:126`), which MEDS already
  handles inside `site_t`.
- **Spatially varying parameters and initial conditions.** Soil texture, PFT maps and census data
  are a later plan (§9, OR2). Here every polygon shares the parameters and the initialisation
  recipe, and differs only in forcing and location.
- **GPU offload.**

## 2. What the code already supports

From the 2026-09-26 survey:

- **One self-contained run object.**
  - `meds_run_t` (`src/main/meds_driver.f90:78-104`) holds the config, one `site_t`, the output
    manager, the fast context, one `met_driver_t`, the budgets, the slow ledger, the clocks and the
    counters.
  - It is a local of the program (`meds_main.f90:31`), and the code says so: "no state hides in
    module scope, so two runs can be open at once" (`meds_driver.f90:75-77`).
  - The main loop is thin: `driver_open` → `driver_step` loop → `driver_finalize`
    (`meds_main.f90:56-72`).
- **State is passed as arguments.** `site_t` is a flat structure of arrays
  (`src/state/site/meds_site_state_types.f90:342`). There is no RNG, and every file unit uses
  `newunit`.
- **The patch loop is already OpenMP-parallel and deterministic:**
  - the directive is at `meds_fast_dynamics.f90:514-818`;
  - per-thread scratch pools are *local* allocatables sized to the thread count, allocated on each
    call (:319-351, :499-512), so each polygon thread naturally gets its own;
  - results are folded back in patch order (:820-856).
- **The build already turns on recursive locals** (`MEDS_AUTO_LOCALS`, `CMakeLists.txt:62-107`,
  mandatory when threading), which removes races on static locals.
- **The forcing format is already multi-cell** (`(time, grid)`, with `grid_index` and nearest
  matching), and the new archive (`MEDS_FORCING_DESIGN.md` §14) is a regular global grid.

## 3. What blocks a parallel polygon loop today, and how each is handled

| # | Blocker (survey citation) | Resolution |
|---|---|---|
| B1 | **netCDF inside the time step.** `driver_step` writes output (`output_serialize_pending`, `meds_driver.f90:359,371`), yearly checkpoints (:412-416), and forcing reads through `met_advance` (`meds_fast_dynamics.f90:482`). HDF5 serialises every call under a global lock and is not thread-safe (`MEDS_IO_DESIGN.md` §5.5). | Split each step into a **compute phase** (no netCDF at all) and an **I/O phase** that runs serially (§4, R1). Forcing becomes an in-memory month buffer (`MEDS_FORCING_DESIGN.md` §15.3), so no forcing reads happen inside the step. |
| B2 | **Module-level allometry.** Nine `protected` reals (`meds_allometry.f90:50-58`) are written by `set_allometry` from every `load_meds_config` (`meds_config.f90:472`). | Load config **once, serially, before** the parallel loop. All polygons share one parameter set (non-goal: per-polygon parameters). A later refactor moves these into the config type (OR2). |
| B3 | **Saved state in the fast-loop probe.** `write_fast_probe` keeps `save` variables (`meds_fast_dynamics.f90:919-920`). | Reject `fast_probe` in polygon mode, as `n_threads > 1` already does (`meds_config.f90:690`). |
| B4 | **Nested OpenMP.** A polygon `parallel do` around `driver_step` nests the patch-loop region. | In polygon mode force `cfg%n_threads = 1` inside polygons (`validate_config` rule), and set `omp_set_max_active_levels(1)`. The patch loop then runs serially with a one-slot pool. The NVHPC `target` loop (`meds_demography_update.f90:128`) must also stay single-level; polygon mode requires `MEDS_GPU = none` for now. |
| B5 | **Shared file names.** Output streams `<dir>/<prefix>-<tier>…nc` (`meds_output_stream.f90:174-176`), restarts `<prefix>-S-…nc` (`meds_io.f90:80`) and `_pft_parameters.csv` (`meds_driver.f90:245`) would collide between polygons. | Polygon mode writes **one file per tier with a `polygon` dimension**, and one ragged restart file per checkpoint (§6). No per-polygon files. `ensure_output_dir` (`execute_command_line`, :599) runs once, serially. |
| B6 | **`error stop` kills the process.** There are 12 sites in the fast kernels, plus `nc_check`, `met_advance` and census parsing. One bad polygon ends the whole run. | R3 fails fast, with the polygon id in the message. R5 converts the common failures to a status return, marks the polygon failed, carries on, and reports the failed polygons at the end (§7). |
| B7 | **Interleaved stdout.** `print_summary`, budget reports and `met_open`/`output`/`state` lines would interleave between threads. | In polygon mode per-polygon printing is off (`verbose` forced low). A per-run summary and a per-polygon status table are written in the serial phase. |
| B8 | **Allocator contention.** About 26 `allocate` calls per patch per `dt_fast` in `build_column_frozen` take about 24% of allocator self-time (#195); more threads mean more heap contention. | The R5 performance item: fix #195 by pre-allocation, or link a thread-scalable allocator (jemalloc or tcmalloc via `LD_PRELOAD`), whichever measures better. |
| B9 | **`BLOCK` constructs.** nvfortran rejects `BLOCK` inside a parallel region (`meds_fast_dynamics.f90:678-680`), and `driver_step` has one (`meds_driver.f90:393`). | The rule concerns constructs lexically inside the region. `driver_step`'s `BLOCK` sits in a routine *called* from the region, which should be allowed. **R0 verifies this with nvfortran;** if it is rejected, the `BLOCK` is hoisted. |
| B10 | **The ifx `private` derived-type mold.** Compiler-generated static molds are shared between threads (`meds_fast_dynamics.f90:324-334`; `docs/building.md:51-54`). | The polygon loop declares **no** `private` or `firstprivate` derived types. Each iteration works on its own element `poly(p)` of a shared array, and scratch lives in called routines. |

## 4. Execution model: month-synchronous, compute and I/O phases

```
serial:   load config once (B2) → build domain: polygon list from forcing box + selection rules (§5)
          allocate poly(1:N): site_t, fast_ctx, budgets, ledger, output accumulators, forcing view
          initialise every polygon (bare ground or restart), serially or in parallel (no I/O inside)
loop over months:
  serial:   forcing reader loads the month block for ALL polygons (chunk-column reads, forcing design §15.3)
  parallel: !$omp parallel do schedule(dynamic, chunk)
              do p = 1, N:  advance poly(p) through the whole month (fast + slow steps),
                            forcing from the shared read-only month buffer at poly(p)%cell,
                            output accumulated in memory in poly(p)'s accumulators
  serial:   drain every polygon's completed output records to the polygon-dimension files (§6)
            yearly: write the checkpoint (ragged restart, §6); flush logs and status
```

- **Why a month is the synchronisation unit:**
  - the forcing reader already works in months (forcing design FD7);
  - a month of one polygon is long enough that the per-region cost of OpenMP is negligible;
  - `schedule(dynamic)` balances polygons that cost different amounts, for example many-cohort
    tropical forest against tundra.
- **Single-site mode is N = 1** through the same code path, with I/O between months instead of
  every step. Output content must be byte-identical to today's; R1 checks this with the existing
  single-site outputs.
- **Output records** that close within the month (daily, and fast-tier records if enabled) are held
  in each polygon's in-memory queue until the serial phase. Fast-tier output is therefore restricted
  in polygon mode (§6).

## 5. Domain and polygon selection

- **The domain** comes from the forcing reader's `domain = "box"` (`box_nwse`, forcing design
  §15.5). Each valid forcing cell inside it is a candidate polygon, in row-major order.
- **Selection rules** (config), since not every valid forcing cell should carry vegetation:
  - **`land_fraction_min`:** the valid mask includes large lakes, for example the whole Caspian
    (6,866 valid cells have `lsm = 0`), and coastal fractional cells. The suggested default is 0.5.
  - **`exclude_ice`:** the Antarctic and Greenland ice sheets are about 0.74 M of the 2.2 M valid
    cells. A glacier mask is added to the static file for this (forcing design §14.3).
  - **An explicit polygon list** (`polygons_csv`, with lat, lon, id) for site networks. Each entry
    maps to its nearest valid cell (forcing design §15.5), and duplicates share a cell but keep their
    own state.
- **Per-polygon attributes** set from the static file: latitude and longitude (cell centre),
  `grid_elevation`, `utc_offset` (from longitude, or zero), and the polygon id (stable: row-major
  index into the global grid).

## 6. Output and restart in polygon mode

- **Output streams:** one file per tier per time chunk for the whole domain, with a leading
  **`polygon`** dimension. `MEDS_IO_DESIGN.md` §10 Q5 already anticipates a leading-axis
  convention.
  - **Coordinates:** `polygon_id(polygon)`, `lat(polygon)`, `lon(polygon)`.
  - **Variables:** only **fixed-shape** ones: site totals `(time, polygon)`, per-PFT
    `(time, polygon, pft)`, per-size-class `(time, polygon, dbh_class)`, and per-soil-layer
    `(time, polygon, soil)`.
  - **Cohort-level and patch-level variables** (ragged across polygons) are not written for the
    domain. Full detail goes only to a short `detail_polygons` list, as ordinary single-site style
    files named by polygon id.
  - **Fast-tier (`F`) output** is limited to `detail_polygons`. Daily, monthly and yearly tiers are
    written for all polygons.
  - **Writing** happens in the serial phase, one hyperslab per variable per record over all
    polygons.
- **Restart:** one netCDF file per checkpoint for the whole domain, using CF **contiguous ragged
  arrays**:
  - per-polygon counts of patches, and per-patch counts of cohorts, plus offsets;
  - the site, patch and cohort fields concatenated.

  Restoring a domain restores each polygon exactly. A single-site restart stays in today's format.
- **Scale check** for 305,000 North American polygons: a monthly record of 50 site-level
  variables is 305k × 50 × 4 B, about 61 MB, trivial. Restart size depends on cohort counts; R0
  measures `site_t` size after spin-up.

## 7. Failure handling, logging, determinism

- **Determinism.**
  - Polygons are independent, so each polygon's trajectory depends only on its inputs, never on the
    thread count or the schedule.
  - Output order is the polygon order, fixed by the domain, not by thread completion.
  - A CTest asserts byte identity between 1 and 4 polygon threads.
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

## 8. Memory, batching and very large domains

- **Memory per polygon** is `site_t` plus the fast context, the accumulators and the forcing view.
  Cohort counts drive it; **R0 measures it** after spin-up on representative sites.
- **Batching:** if N × memory per polygon exceeds `memory_limit_gb`, the domain runs in **batches**
  of `batch_size` polygons. Each batch runs its full simulation period (every month, with I/O as in
  §4) before the next starts.
  - Batches are built from **spatially compact tiles**, aligned to the forcing chunks (16 × 16
    cells), so each batch's monthly forcing read touches few chunks.
  - The outputs of all batches go into the same domain files, as disjoint polygon ranges.
- **Job arrays:** because polygons are independent, a continental or global run can also be split
  into tile groups submitted as an HPC **job array**. Each job runs one tile group in polygon mode
  and writes its own files, and a small merge step (or just file lists) assembles them. This scales
  like MPI with none of its machinery.
- **The forcing buffer per batch** is 8 variables × cells × 744 h × 4 B: about 0.5 GB for 20,000
  cells, and about 7.3 GB for all North American land.
- **Loading it** with chunk-column reads (forcing design §15.3) took **5.0 s per variable-month**
  for all of North America, so about 40 s for all 8 variables in the serial phase. That is small next
  to a month of compute for 300,000 polygons.

## 9. Configuration sketch

```toml
[run]
mode = "polygons"                 # "site" (default; today's behaviour) | "polygons"

[forcing]                         # forcing design §15.2
met_source = "era5land"
data_path  = "<installation-specific path to the ED_ERA5land archive>"
domain     = "box"
box_nwse   = [45.1, -79.8, 40.4, -71.8]

[polygons]
land_fraction_min = 0.5           # skip lakes and fractional coastal cells
exclude_ice       = true          # skip ice-sheet cells (static glacier mask)
# polygons_csv    = "sites.csv"   # alternative to the box: explicit site network (id, lat, lon)
threads           = 32            # OpenMP threads over polygons; patch threading is forced to 1 inside (B4)
schedule_chunk    = 4
memory_limit_gb   = 64            # beyond this, the domain runs in batches of spatially compact tiles
detail_polygons   = []            # ids that also get full cohort/patch and fast-tier output
```

`validate_config` enforces the following in polygon mode:
- `n_threads = 1`;
- no `fast_probe`;
- `MEDS_GPU = none`;
- `met_source` must not be `legacy_file` with more than one polygon (legacy files hold only a few
  cells).

## 10. Phases

| Phase | Content | Acceptance |
|---|---|---|
| **R0** measure and verify | The wall time of one single-site simulated year, in the default configuration; `site_t` memory after spin-up; the allocator profile (#195); nvfortran on a `BLOCK` in a routine called from a parallel region (B9) | Numbers recorded here. B9 answered. |
| **R1** compute/I-O split | `driver_step` split into a compute phase (no netCDF) and an I/O phase (serialise pending records, checkpoints). Forcing reads only at month boundaries (needs forcing design F4, §15). Single-site behaviour unchanged. | The existing CTest suite is green on ifx and nvfortran; single-site outputs are **byte-identical** to before. |
| **R2** domain and polygon container, serial | `meds_domain_t` (shared config, shared month buffer, polygon list) and `meds_polygon_t` (per-polygon state and accumulators). The month-synchronous loop without OpenMP. Polygon-dimension output for fixed-shape variables. `detail_polygons`. | N polygons run serially produce outputs **byte-identical** to N separate single-site runs (a 4-polygon synthetic test). |
| **R3** OpenMP polygon loop | `!$omp parallel do schedule(dynamic)` over polygons; the B3, B4 and B7 rules in `validate_config`; fail-fast messages with the polygon id | Byte identity between 1 and 4 threads; scaling measured to the core count; nvfortran build green. |
| **R4** ragged restart | Domain checkpoint and restart with CF contiguous ragged arrays | A restart round trip is bit-identical against an uninterrupted run. |
| **R5** robustness and scale | Status-based failure isolation (B6), batching by tiles, the job-array recipe, allocator work (B8), per-polygon logs | A failure-injection test: one polygon fails and the rest match the reference. Throughput and memory recorded for a 20,000-polygon box. |
| **R6** interfaces and docs | Optional C API and Python entry point for polygon runs (the current C API registries hold at most 4 runs and are unsynchronised, `meds_c_api_run.f90:46-48`); docs; an example | A documented example runs a small box end to end. |

**Order:** R0 anytime. R1 needs forcing design F4 (§15, monthly block reads). R2 and R3 need R1. R4 and R5
follow R3.

## 11. Tests (CTest)

- **Serial equivalence:** 4 synthetic polygons as a domain produce byte-identical outputs to 4
  single-site runs (R2).
- **Thread invariance:** 1 against 4 polygon threads, byte-identical (R3).
- **Month seam:** forcing continuity across a month boundary with the shared buffer, reusing the
  forcing design §15.6 files (R1).
- **Restart round trip** of a domain with different cohort counts per polygon (R4).
- **Failure isolation:** inject a NaN into one polygon's forcing; that polygon is flagged and the
  others are unchanged (R5).

## 12. Open decisions

| # | Decision | Recommendation |
|---|---|---|
| OR1 | The variable set for polygon-dimension output | Start from the site-level and PFT/size-class variables that exist in `meds_io_config.toml`, as a `polygon_output` group. |
| OR2 | Spatially varying inputs (soil, PFT composition, parameters) | A separate plan after R3. It needs the allometry refactor (B2), moving parameters from module state into the config. |
| OR3 | Failure policy default: abort or skip-and-flag | Abort during development (R3). Skip-and-flag is available from R5, selectable by config. |
| OR4 | Batching inside one process or job arrays across tiles | Both. Batching for a workstation or a single node; job arrays for continental and global runs on HPC. |
| OR5 | Whether to simulate lake and ice cells at all | No by default (`land_fraction_min = 0.5`, `exclude_ice = true`). Revisit with a lake or glacier surface scheme. |

## 13. Relationship to other documents

- **ROADMAP #183** (multi-polygon runtime with MPI): this plan delivers everything except MPI. The
  ROADMAP entry is updated when R3 ships.
- **`MEDS_FORCING_DESIGN.md`:** its §3.2 and §8 P2 ("one `met_driver_t` per polygon") are
  superseded by one shared reader with per-polygon views of the month buffer. Its §10 Q7 (a polygon
  with no covering land cell) is answered by its §15.5 and the selection rules here.
- **`MEDS_GPU_EVALUATION.md`** §10 ("MEDS has no site or ensemble axis at all"): the polygon axis is
  exactly that axis. Its §7 and §12.4 (allocator traffic, #195) become more pressing with polygon
  threads (B8).
- **`MEDS_IO_DESIGN.md`** §5.5 (the HDF5 global lock) and §10 Q5 (the leading-axis convention) are
  resolved by the serial I/O phase and the `polygon` dimension.
