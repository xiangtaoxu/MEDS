# MEDS numerics and efficiency sweep — plan

**Revision 4, 2026-09-30.**
- Revision 2 took in the owner's comments and decisions D1–D7 (§5).
- Revision 3 took in the second round: the forcing-writer layout (§5.1), the flattened pool deferred
  and patch threads compared with polygon threads (§2.7), the ecophysiology of the fixed wood water
  (Appendix A), and how option B works (Appendix B).
- Revision 4 records the owner's final choices:
  - threads over polygons for region runs (Phase 5);
  - #104 by option (b) in this sweep (Phase 3);
  - #146/#190 by option B (Phase 2).

  Every decision in §5 is now taken.
- **Phase 3 finding, 2026-09-30 (during implementation).** Steps 1 and 2 of Phase 3 are done. Step
  3, the physiology change, is **paused for the owner**: the instrument and a closer reading of the
  sources change the case for option (b). See "Phase 3 status" under Phase 3.

**Against:** `beta` = `main` = v0.3.1 (`eea42ec`). **Status:** being implemented on branch
`dev/efficiency-sweep` (§0.1). Every number in §2 was measured on an unmodified `git archive eea42ec`
build, or read from the code (marked *code reading*). The measurement harness is kept outside the
repository, with the profiles and run scripts: `~/claude_workspace/meds_efficiency_plan/` on the
development cluster (`harness/`, `profiles/`).

**Scope.** All 40 open issues were read. This sweep takes 16 of them that are numerics, compute
cost, or the code structure those depend on (§1), plus six new items the measurements turned up
(§3.7). Three numerics issues were closed on 2026-09-30 (D4). The rest are science or feature work
and stay out of the sweep.

---

## 0. Summary

1. **#325 comes first.** At 8 threads, 83% of all CPU time goes to ifx's allocator for "bound
   procedure values": 232 of 280 CPU-seconds on the 60-day BCI case. Of that, 87% comes from
   `flux_potential`/`kirchhoff_edge` and 13% from `solve_leaf_gas_exchange`. Run serially, the same
   cost is only 7% of the fast loop. The prototype patch still applies (offset by one line).
2. **Forcing files are written with 1×1 chunks** (new). The shared writer
   `scripts/forcing_common/meds_forcing_file.py` makes `time` unlimited and takes the default
   chunking, so opening the BCI file reads 90,528 one-value chunks per variable. Rechunking that
   file to (8760, 1) has three effects:
   - each run starts 3.3 s sooner;
   - the file shrinks from 63 MB to 4.2 MB;
   - every output variable stays bit-identical.

   That is roughly a third of each 10-day calibration trial (median 8.5 s on an idle node).
3. **Allocator traffic is real, but not where #195 says** (#195). It is 15% of serial fast-loop CPU.
   The largest site is `surface_derivs` (5%), not `build_column_frozen` (2%), and the state
   combinators and the ARK stage locals make up most of the rest. #188's copy is now a heap deep
   copy of the hydraulics table every step.
4. **Numerics.**
   - #158, #159 and #167 are closed (D4); the integrator is revisited in a later round (§3, Phase 4).
   - #104's cost halves with #325; the rest waits on the owner's physics choice (Appendix A).
   - #189 item 3 is a small, bit-identical fix.
5. **Structure.**
   - **Region threading.** The pool of all patches across polygons is deferred. Of the two
     simple shapes (§2.7), threads over polygons give 17–25× on a 40-core node. Patch threads inside a
     plain polygon loop give 1.4–3.3×, because Ithaca stands have only 2–6 patches. Neither needs a
     compilation flag (§2.6, §2.7, Phase 5).
   - **#146 no longer has a speed case of its own**; Appendix B gives the options.
   - **#311 F8 is housekeeping, not performance** (milliseconds per simulated year).
   - **#312 is cleanup,** and it must not collide with the region work.

**Phase order:**
- **0:** baselines.
- **1:** quick wins: #325, forcing chunking, #189 item 3, count the corrector's work, docs.
- **2:** one reusable per-thread workspace: #188, then #195, with #146 (option B) and #190.
- **3:** #104 by option (b): the wood's apoplastic water drains as its conduits embolise (a
  physiology PR).
- **5:** runtime consolidation and polygon threads for region runs: #310, #183, F10, then #196.
- **6:** output and forcing cleanup: #312, #311, #270, #299, #275.

Phase 4 (integrator re-measure) moves to the later numerics round.

---

## 0.1 Working rules for every change in this sweep

**Code (owner's comment 2).**
- **Readability and modularity.** Every revision keeps the code at least as readable and modular as
  it was, and improves it where it can.
- **Choosing a fix.** Prefer the solution that consolidates, simplifies and removes special cases.
  When two fixes perform alike, take the one that deletes more code, or more special cases.
- **The PR states the effect.** Each PR says what it consolidated or removed: routines merged,
  branches deleted, lines net.
- **Comments are plain language.** A comment says what the code does and why, in words a modeller
  reads without a numerical-analysis or compiler background.
  - Avoid unexplained jargon and acronyms; for example, say "a procedure passed as an argument"
    rather than "BPV thunk".
  - Name the physical quantity rather than the implementation trick.
  - Keep CLAUDE.md's rules: present-tense rationale; history goes in the CHANGELOG, deferrals in
    the ROADMAP.
  - Shorten long banner comments that repeat a design document. Cite the document's section instead.

**Cluster (owner's comment 4).** Every build, test and timing runs on partition `R128C40` only
(`cbsuxu01`–`08`), with `--partition=R128C40 --exclusive --mem=0` for timings, one run per node.
Never `R256C128`.

**Numbers (D1).** Rounding-level movement of outputs is acceptable in every PR. The PR states its
size: the first day's largest relative difference, and the five-year tower statistics, budgets and
stand summary. A change bigger than rounding must be explained. Pure refactors should still come out
bit-identical on short windows. If one does not, say so, rather than hiding it under "rounding".

**Git (D6, D7).**
- No releases and no version bumps.
- One development branch, `dev/efficiency-sweep`, off `beta`: one commit per item, and pull requests
  into `beta`. CHANGELOG entries go under `[Unreleased]`.
- No new GitHub issues. New findings (§3.7) are noted in the commit and PR messages.

---

## 1. Triage of the open issues

### 1.1 In the sweep

| # | Title (short) | Phase | Verdict |
|---|---|---|---|
| 325 | OpenMP >4 threads slower (bound-procedure allocations) | 1 | Do now. Measured at 83% of CPU at 8 threads |
| 188 | column_params_t copied into the frozen record | 2 | First commit of #195. It is now a heap deep copy (#179) |
| 195 | Allocator traffic in build_column_frozen | 2 | Do it, re-scoped: 15% of the fast loop, spread over the ARK march |
| 146 / 190 | Packed state / cohort field table | 2 | Option B: one field list and one rules table (Appendix B) |
| 104 | Hydraulics thrash on a collapsed wood store | 3 | #325 halves its cost. Then option (b): wood apoplastic water drains with embolism (Appendix A) |
| 189 | Per-cohort tissue residuals; RK45 ledger order | 1 (item 3), 3 (item 2, water half) | Item 3 is tiny and bit-identical. Item 2 is #104's instrument |
| 310 | Runtime consolidation (R2–R12, O11) | 5 | The prep PR is bit-identical; then one run container |
| 183 | Multi-polygon runtime, R3–R6 | 5 | R3 by threads over polygons, with no compilation flag; the flattened pool deferred (§2.7) |
| 196 | Thread and vectorise the cohort axis | 5 | After R3. Vectorise `surface_derivs`; cohort tasks only when patches < cores |
| 311 | Forcing source interface (F6–F12) | 5 (F10), 6 | F10 before R3; the rest is cleanup |
| 312 | Output serializer (O6–O12) | 6 | O7+O8, then O6+#270, then O9+#299; O10 and O12 separately |
| 270 | FAST tier has no patch axis | 6 | Goes with O6 |
| 299 | Slow-only fast rows read 0, not missing | 6 | Remaining half goes with O9 |
| 275 | Skin temperature; within-step variance | 6 | Remaining half follows O6/#270 |

### 1.2 Closed on 2026-09-30 (D4)

#158 (N5 adaptive freeze cadence), #159 (soil water into the ARK tableau) and #167 (free-convection
slope) are closed. Each got a comment giving the measured reason and the plan to revisit the
numerical scheme in a later round. ROADMAP §4 and §5 still list them; Phase 1e updates those
entries.

### 1.3 Settled earlier; do not reopen

Decided by measurement in `MEDS_PRODUCTION_INTEGRATOR_PLAN.md` and the closed issues:

- the split integrator retired;
- N2a's ground analogue, N2b/N2c-gah, N2e and all of N6 refuted;
- N3 and N4 not to be built;
- E1 as first written refuted, as is loosening the hydraulics tolerance (it flips W_leaf between
  clusters);
- E2 as first written, E5 (#161) and the RK45 changes (#160);
- no higher-order tableau, and no hydraulics in the tableau;
- FP32 no (#197); MB2 (#163); the bare-array conversions (#164); #93's Phases 2–4.

### 1.4 Out of this sweep (science and features)

#1, #74, #96, #154, #155, #156, #157, #165, #180, #181, #186, #187, #254, #255, #256, #257, #258,
#265, #268, #269, #302, #316. #316 (the necromass ledger declaration) is a small correctness fix
worth doing on its own, but it is not efficiency.

---

## 2. Baseline: what v0.3.1 costs today

The build is ifx 2026.1 Release with `-g`, one run per idle `--exclusive` node on `R128C40`. The
case is BCI from the census, 2012-08-01 to 10-01, starting at 25 patches and 419 cohorts. The
harness is in `harness/` and the reports in `profiles/` of the directory named in the header.

### 2.1 Wall time

| case | serial build | OpenMP 1 thread | 4 threads | 8 threads | 16 threads |
|---|---|---|---|---|---|
| 60 days | 34.6 s | 36.2 s | **30.0 s** | 50.9 s | 67.9 s |
| 5 years, default parameters | | | 5:17 | | |
| 5 years, calibrated parameters | | | 6:02 | | |

### 2.2 Where the serial time goes

VTune, 60-day case, 26.4 CPU-seconds in all:

| block | share | notes |
|---|---|---|
| **start-up (driver_open)** | **30.6%** | once per run |
| – census restructuring | 14.9% | every patch fusion reorders every cohort array, about 1,200 times |
| – forcing read | 12.3% | per-chunk overhead of the 1×1 chunks, not disk |
| **fast loop** | **67.5%** | 17.8 s |
| – ARK march | 34.8% | Newton 17.1% (`surface_derivs` 11.6% + Jacobian 5.4%); **the hydraulics re-solve in the corrector 9.2%** |
| – build_column_frozen | 24.9% | leaf gas exchange 8.0%; the scratch soil-water solve 6.4%; the hydraulics pre-pass 3.9% |
| – canopy radiation | 4.2% | |
| output | 1.5% | |

Cross-cutting, as a share of the fast loop:

- **allocator: 15.3%.** By caller: `surface_derivs` 0.88 s, `build_column_frozen` 0.39,
  `state_init` 0.35, `column_fast_step_ark` 0.30, `ark2_column_step` 0.16, `column_be_stage` 0.14,
  the other combinators 0.22, `alloc_rad_flux` 0.08.
- **bound-procedure allocations: 7.4%.**

Separately, `pow` in the soil hydraulic functions is 6.5% of all CPU.

### 2.3 At 8 threads

280 CPU-seconds for 51.7 s of wall time: 5.2 of the node's 40 cores are effectively used.

- **Bound-procedure allocations: 232 s (83%).** 200.8 s from `flux_potential` and `kirchhoff_edge`,
  31.2 s from `solve_leaf_gas_exchange`.
- **OpenMP imbalance: 10.2 s.** The allocations hide it today. It becomes the next limit, because
  14–15 patches cannot feed 16 threads.

### 2.4 Solver work over five years

Daily `work_*` counters, per patch, area-weighted, run at 4 threads:

| counter | default: mean/day | default: max | calibrated: mean/day | calibrated: max |
|---|---|---|---|---|
| accepted ARK sub-steps | 201 (2.09 per dt_fast) | 233 | 195.5 | 216 |
| rejected ARK steps | 4.5 (2.2%) | 31 | 1.6 | 15 |
| soil-water sub-steps | 110 (1.14 per dt_fast) | 557 | 109 | 554 |
| hydraulics sub-steps (pre-pass only) | 1,880 | 1,989 | 1,882 | 2,408 |
| hydraulics thrash | 0 | 0 | 0 except on 6 days | 0.86 |
| whole-site energy/water fails | 0 / 3.0 M checks | | 0 / 3.05 M | |

The stand ends at 14 patches and 232 cohorts. #104's regime is rare with both shipped parameter
sets, but common for calibration candidates.

### 2.5 Findings from reading the code

- **The hydraulics solve runs 1 + (ARK trials) times per dt_fast.** The corrector re-solves it on
  every trial, including rejected ones (`meds_fast_be_stage.f90:429`). Its sub-step and convergence
  counts (`nsub_c`, `conv_c`, lines 394 and 435; checked) are discarded, so the work counters and
  the thrash detector see the pre-pass only.
- **Only three routines pass an internal procedure as an argument.** They are `flux_potential`,
  `solve_leaf_gas_exchange` and `phi_inverse`; the last is used only by a test.
- **Heap allocations per patch per dt_fast.** A code count gives about 260; the plan's LD_PRELOAD
  measurement gave 106. Phase 0 re-counts.
  - `surface_derivs`: 5 per call, 10 or more calls per stage.
  - `column_state_t` constructions: about 17.
  - `build_column_frozen`: 46, where #195 says 26.
  - `apply_rt_forcing`: 10.
- **`apply_rt_forcing` runs a selection sort (cost ∝ cohorts²) every sub-step** (checked).
- **#162's old figures survive in four places:**
  - ROADMAP §4;
  - `numerical_scheme.md` §5a/§6/§7.1;
  - the dt_fast > 900 s warning at `meds_config.f90:733-738`;
  - `examples/example_biophysics/meds_config_july.toml:40-45`.

  #91 fixed psi_leaf, and `test_psi_dt_convergence` guards it (checked).
- **Local, not in the repo:** the main checkout's `build-ifx` caches `MEDS_OPENMP=OFF`, so
  `run_example.py`'s default binary ignores `n_threads`.

### 2.6 What threading a region could buy (owner's comment 3)

The measurements cover the 100-cell Ithaca box, one year (2016) from bare ground, the v0.3.1 serial
build, one `R128C40` node per case. The harness is `harness/region_exp/` (see the header).

| measurement | result |
|---|---|
| each of the 100 cells run alone, one after another | 1,224 s in all; mean 12.2 s per cell, range 10.3–33.1 s (spread 32%) |
| the same 100 runs, 20 at a time on one node | 95.9 s (12.8× faster); each run 1.55× slower than alone |
| the same 100 runs, 40 at a time on one node | 71.1 s (17.2× faster); each run **1.93× slower** than alone |
| the region run, serial, as today (with its daily and monthly output) | **1,012 s**: 17% less than the 100 separate runs, because the polygons share one forcing read |

**What that means:**

1. **The node itself caps any parallel shape at about 20× on 40 cores** for this workload. With 40
   runs at once, each is 1.93× slower: memory bandwidth, shared cache and clock speed are shared.
   Part of it may be each process reading the forcing archive for itself, which a threaded region
   run would share. So a threaded region run could come in a little above 20×, but not near 40×.
2. **A thread per polygon** (R3 as first tabled: `!$omp parallel do` over polygons once per month)
   loses time only in the tail, when the last expensive polygons finish while other cores are idle.
   The flattened pool (one task per patch-day, over all polygons) has a much smaller tail, but it
   can use only as many cores as there are patches.

   A schedule simulation with these per-cell costs, 40 threads, costliest tasks first. It assumes a
   polygon's cost is spread evenly over its months and patches:

   | polygons | node use, thread per polygon | flattened, 2 patches per polygon (bare ground) | flattened, 6 patches per polygon |
   |---|---|---|---|
   | 10 | 25% (only 10 cores can work) | 49% | 75% |
   | 40 | 43% | 82% | 94% |
   | 100 (this box) | 86% (59% if polygons start in random order) | 90% | 97% |
   | 400 | 94% | 97% | 99% |
   | 1,000 | 98% | 100% | 100% |

3. **So the flattened pool buys, over a thread per polygon:**
   - **5–15% for a region of about 100 polygons,** if the costliest polygons start first;
   - **about 5% or less** at 400 or more;
   - **2–3× only for small regions** with fewer polygons than cores;
   - **nothing for a single site**: a one-polygon pool is today's patch loop.

   Regions larger than a node already run as tiles in a job array (polygon plan §8), so a node
   usually holds hundreds to thousands of polygons. That is where the two shapes tie.
4. **Against what users can already do,** threading a region is a smaller win than against serial
   region mode.
   - Region mode serial today: 1,012 s. Threaded, either shape: roughly 55–70 s (about 15–18×).
   - The 100 cells as separate runs, 40 at a time, with no code change: 71 s.

   So threading region mode is worth about 1.0–1.3× over running the cells as separate processes.
   Its real gains are elsewhere: one shared forcing read (17% of the work here), one output file per
   tier instead of one per cell, and no hand-made tile configs.

**Readability.**
- **A thread per polygon** is about 30 lines around `polygon_step`, plus a rule that turns patch
  threads off inside it. The fast loop stays as it is.
- **The flattened pool** splits the 690-line `fast_dynamics` into three parts (forcing samples,
  patch body, sums). Each part is clearer than today's single routine, but the day's state then
  has to cross three routine boundaries: the per-polygon staging arrays and forcing samples, and a
  day-by-day lock-step of all polygons. It is the bigger, riskier change, and it adds indirection:
  a task number that maps to (polygon, patch).

**Decision (owner, 2026-09-30):** the flattened pool is **deferred**.

### 2.7 A plain loop over polygons: patch threads inside, or polygons side by side?

This is the comparison the owner asked for, between two regions made of the same polygons.

- **(A)** A plain, serial loop over polygons. Each polygon's fast loop threads its patches, as a site
  run does today.
- **(B)** Threads over polygons, one polygon per thread, with no patch threads inside.

**What each would take.**
- **No new compilation flag for either.** The build has one OpenMP switch (`MEDS_OPENMP`, on by
  default), and which loop is threaded is decided when the model runs, from the mode:
  - **(A)** needs no new code. Remove the region-mode `n_threads = 1` refusal
    (`meds_config.f90:799`), whose own message says it is a placeholder until R3, and test that
    region runs with patch threads give the same results.
  - **(B)** needs one `!$omp parallel do` around the polygon loop (about 30 lines with its
    scheduling). It also needs a rule that runs the patch loop with one thread inside: MEDS already
    passes the patch loop its thread count explicitly (`num_threads(n_thread)`), so the rule sets
    that to 1.
- **The real work in (B)** is checking that nothing a polygon step touches is shared between polygons.
  The polygon plan (§3, blockers B1–B12) lists those things, and R1/R2 already removed most of them:
  no file I/O inside the step, per-polygon output buffers, no saved variables.
- **One setting, `n_threads`,** would mean patch threads in a site run and polygon threads in a
  region run.

**Measured** on one `R128C40` node, v0.3.1, each polygon one year:

| polygon | alone, 1 thread | (A) patch threads | (B) 40 polygons at once, 1 thread each |
|---|---|---|---|
| established Ithaca stand (2 patches, 14 cohorts) | 8.4 s | 6.0 s with 2 threads (1.4×), 6.0 s with 4 | 40 polygons in 13.6 s: **25×** the serial throughput (each run 1.57× slower than alone) |
| bare-ground Ithaca cell (6 patches, no cohorts) | 12.0 s | 6.3 s with 2, 5.7 s with 3, **3.6 s with 6** (3.3×) | 100 polygons in 71 s: **17×** (each 1.93× slower) |
| BCI stand (14–25 patches; 60 days) | 34.6 s | 30.0 s with 4 (1.15×) today, limited by #325; the #325 prototype gave about 2× at 8 threads | not measured |

**Why they differ.**
- **(A) can never use more threads than a polygon has patches,** and patches differ in cost.
  Temperate stands have 2–6 patches, so (A) buys 1.4–3.3×. A many-patch tropical stand can do
  better once #325 is fixed.
- **(B) is capped only by the node.** Running many polygons at once slows each by 1.6–1.9×, so 40
  cores give 17–25×.
- **For any region of more than a handful of polygons, (B) is 5–15× faster than (A).** (A) is
  worth having only for regions with fewer polygons than cores that are made of many-patch stands.

**Complexity.**
- (A) adds nothing.
- (B) adds one parallel loop and one rule, and site runs keep their patch loop as it is.
- Neither needs a compilation flag or a second build.

**Decided (owner, 2026-09-30): (B), threads over polygons,** for region runs (Phase 5, step 4).
Option (A) is not pursued, and the flattened pool stays deferred.

---

## 3. The plan, by phase

### Phase 0 — Baselines and harness (no source change)

1. **Regression references at `eea42ec`.**
   - The r1 case files predate v0.3.1's renames. They still say `format = "era5land"` and carry
     `[site].utc_offset`/`apply_solar_longitude`, which v0.3.1 refuses. Update them first; the
     region demo configs needed the same edits.
   - Then record default and `-fp-model consistent` references.
2. **Fill the r1 gaps.** Add a threaded case (1/4/8), the BCI census start (60 days), an RK45 case,
   and a threaded region case (for Phase 5).
3. **Benchmark set,** every run alone on an `R128C40` node:

   | id | case |
   |---|---|
   | B1 | BCI 60 days at 1/4/8/16 threads |
   | B2 | BCI 5 years, default and calibrated |
   | B3 | one calibration trial (10-day restart) |
   | B4 | the 100-cell region box (§2.6) |
   | B5 | `test_fast_loop` thread invariance |

4. **Re-count heap allocations per patch-step,** to settle 106 against 260.

### Phase 1 — Quick wins (each a small PR into `beta`)

**1a. #325, both halves.**
- **`flux_potential`:** apply the prototype, which moves the quadrature into its own function, or
  write the 7-point sum inline. Either way no procedure is passed on the fast path.
- **`solve_leaf_gas_exchange`:**
  - move its two residuals to module level with one small, explicitly passed record of what they
    need (22 reals, 3 integers, 3 logicals);
  - give that residual shape a bisection;
  - convert `phi_inverse` the same way.

  **This consolidates:** one residual shape and one bisection for the three root-finds.
- **Stop it coming back:** a ctest (or a grep in CI) that fails if a contained procedure is passed
  as an argument anywhere in `src/`.
- **Expected:** 8 threads from 50.9 s to about 12–15 s on B1. #104's collapsed-state cost roughly
  halves.
- **Gate:**
  - B2 faster at 8 and 16 threads than at 4;
  - `test_fast_loop` invariance holds;
  - the CHANGELOG states the rounding movement;
  - the "at most four threads" advice is removed from `docs/building.md` and the BCI README.

**1b. Forcing chunking.**
- **Where the writer lives (the owner's decision, §5.1).**
  - The `forcing_common/` folder goes.
  - The writer becomes one script, `scripts/meds_forcing_file.py`: `write_forcing_file`,
    `check_complete` and the constants they use.
  - The tower-only copies of the model's conversions (humidity, pressure, solar geometry, the
    longwave synthesis) move into `scripts/prepare_flux_tower/`.
  - Update the imports in `prepare_era5/make_forcing_file.py`, the tower scripts and
    `calibrate_fast/tests/test_smoke.py`, and the references in `docs/science/forcing.md`, the tower
    README and the CHANGELOG.
- **The chunking.** The writer writes a fixed `time` dimension, or passes
  `chunksizes=(min(n, 8760), ngrid)`. One fix covers every tool that writes a forcing file.
- **The note (D2).** `met_open` prints a one-line note when a file's time chunk is a single record,
  and the note says why it is printed. Proposed wording:

  > `note: <file> stores one time record per chunk, which makes reading it slow (about 3 s per
  > 90,000 records). Rewriting it once with "nccopy -c time/8760,grid/1 in.nc out.nc" makes it
  > about 100x faster to open, with identical values.`

- **Measured:** −3.3 s per run; outputs bit-identical.

**1c. #189 item 3.** The RK45 ledger check decides whether to stop only after the step is accepted
or rescued. About 30 lines, bit-identical. Add a test proving that a committed breach still stops.

**1d. Count the corrector's hydraulics work** into `work_hydro_*` and `work_nonconv`. Physics is
unchanged; only those counters change.

**1e. Docs.**
- Correct the stale #162 figures in the four places listed in §2.5.
- Update ROADMAP §4/§5 for the three closed issues, per the ROADMAP rule that every deferred item
  names an open issue: move them to "revisit with the numerical scheme".
- Correct #325's "per stage" wording.
- Fix three out-of-date comments: the `meds_fast_frozen.f90:55-77` banner, `numerical_scheme.md`
  §3's `ark_niter`, and the `column_fast_step_ark` header.

### Phase 2 — One reusable per-thread workspace (#188, then #195)

1. **#188 first.** The frozen record points at `col_config`'s static parameters instead of copying
   them. Most of all, stop the deep copy of the hydraulics table (n_pft × 4.2 KB) every step.
2. **#195.** Add a `fast_workspace_t` per thread, grown only when a patch has more cohorts than
   before. It holds:
   - the frozen record;
   - the ARK stage states;
   - the `surface_tend_t` records;
   - the tissue flux integrals;
   - the radiation records.

   `build_column_frozen`, `state_init`, `surface_derivs` and `canopy_radiation` then fill it in
   place, and the per-stage deep copy `sf_out = surf_tend` goes. **This consolidates:** one
   workspace instead of about 260 allocations per patch-step, scattered over ten routines.
3. **#146 by option B (decided, owner 2026-09-30; Appendix B).**
   - One pair of routines beside `column_state_t` copies its fields to and from a flat array, and one
     rules table says, per entry, whether it steps, counts in the error, or is clamped.
   - The ~15 bookkeeping routines become array operations on the workspace's arrays.
   - It rewrites the same routines as #195, so doing both here avoids touching them twice.
   - **#190 in the same phase:** the cohort gather's per-field fusion and scaling rules
     (`fuse_cohort_fast_state`, `scale_cohort_ground_fields`) move into the same kind of table, if
     that deletes code. Otherwise #190 closes in that PR, citing the four of its five benefits that
     have already landed (ROADMAP §10).
4. **`apply_rt_forcing`'s sort.** Use the height order the cohort block already keeps, and drop the
   per-sub-step selection sort. The order must stay identical.

**Expected:** 10–15% off serial fast-loop CPU, and better thread scaling, since the allocator also
contends. **Gate:** allocations per patch-step close to 0 in steady state.

### Phase 3 — Hydraulics: #104

**Decided (owner, 2026-09-30): option (b), in this sweep.** The wood's apoplastic water drains as
its conduits embolise. This is a physiology change, so it goes in its own PR.

**Phase 3 status (2026-09-30): steps 1 and 2 done; step 3 paused for the owner.**

- **Step 2 is in** (`7c7b481`): `water_curve_t` and the split into `meds_water_retention.f90`,
  bit-identical. A related fix came out of it (N-7, `10d6482`): the daily tissue-water reconcile now
  uses each PFT's own curve.
- **Step 1, the instrument** (a temporary build, not committed). Per patch-step it recorded the
  soil-throttle scale, and per cohort whether the committed wood water ended below the curve's
  apoplastic floor (`apoplast_frac × W_sat`), and whether it had been above it at the step's start:

  | run | patch-steps | throttled | wood below the floor (cohort-steps) | entries | on a throttled step | emptied to 0 |
  |---|---|---|---|---|---|---|
  | BCI five years, calibrated | 3.05 M | 18% (scale down to 0.11) | 0 | 0 | – | 0 |
  | BCI five years, line-searched (`wood_psi50` −1.24, `leaf_pi0` −2.43) | 3.05 M | 18% (to 0.085) | 173 of 55 M | 7 | 7 | 0 (lowest 0.49 × floor) |
  | Ithaca stand, 400 days, no rain (constant forcing) | 0.10 M | 17% (to 0) | 1,047 | 64 | 64 | 25 |

  - **The drain path in Appendix A is confirmed.** Every one of the 71 entries below the floor
    happened on a step whose uptake the soil throttled; the leaf never went below its floor.
  - **The cost has mostly gone.** Since #335's stomatal closure, the calibrated BCI stand never
    reaches the floor, and the line-searched set does so 3 times in a million cohort-steps. On those
    patch-steps the hydraulics take 47 sub-steps per cohort, against 3.1 otherwise, so the total cost
    is negligible. Only a stand with no water at all collapses.
- **The literature, read more closely, does not support the bound fraction step 3 planned:**
  - **TFS-Hydro's sapwood residual fraction is not a small bound part.** Christoffersen et al.
    (2016) count the water released by embolised conduits in the drainage *above* the residual
    fraction, which is where the potential reaches −∞, exactly as MEDS's curve does at its
    apoplastic fraction. Their Table 2 formulas give that residual as about 0.30 of saturated water
    at wood density 0.4 and 0.45 at 0.6, larger than MEDS's whole 0.20. It cannot be the bound part
    of MEDS's apoplastic water.
  - **Wood science puts more water in the bound part, not less.** Cell-wall water at fibre
    saturation is about 0.28 g per g of dry wood. MEDS's saturated wood water, 1.0 kg per kg of
    carbon, is about 0.48 g per g of dry wood. So cell-wall water is roughly 58% of MEDS's
    saturated wood water, more than the 20% MEDS already holds fixed.
  - **What that means.** MEDS's "fixed" 20% is not conduit water waiting to drain. The conduit water
    is already in the part of the curve that drains. Appendix A's point 4 ("the wood's 20% lumps
    drainable conduit and capillary water with the small amount that is truly bound") was wrong
    for MEDS's parameters.
  - **Only SurEau** (Cochard et al. 2021) drains all of the apoplastic water with embolism, a bound
    fraction of 0.
- **Options for the owner:**
  - **(b) as decided,** with a bound fraction of 0 (SurEau). This releases water that wood science
    counts as bound, and it moves ordinary drought behaviour, including the calibrated BCI runs,
    which never reach the floor today.
  - **(c) stop the drain where it starts.** The instrument points at the throttle, which cuts the
    wood's uptake but not its sap flow. This is the only option that fixes the cause.
  - **(a) bound the potential.** A finite floor in place of the −10⁴ MPa marker. It is small and
    changes no run that stays off the floor, which is now every realistic run seen.
  - **Close #104 as resolved by #335,** with this table, and file nothing new.
- **Recommendation.** Do not change the wood's physiology for #104: the evidence for (b) is weak,
  and the cost it was meant to remove is gone. Either close #104 with the table above, or take (a)
  now as a guard for dry runs and leave (c) to the later numerics round (Phase 4's list).

1. **Instrument first** (state unchanged). Use #189 item 2's water half: per-cohort committed mass
   against the kernel's endpoint. Log per step the wood water against its apoplastic minimum and the
   soil-throttle scale, on the calibrated and the line-searched sets. This records the before
   picture and confirms or refutes the drain path in Appendix A.
2. **One record per tissue for the water curve** (consolidation first; bit-identical).
   - Today about ten files call `water_content`, `capacitance` and `psi_from_water_content` with the
     same five loose traits (π₀, ε, apoplastic fraction, saturated water, biomass).
   - Give leaf and wood each one small record of their curve traits, and pass that record instead of
     the loose traits.
   - Step 3 then adds the embolism traits in one place, not at every call.
   - **Where it lives (owner's question, 2026-09-30).** The curves are already in `shared/functions`
     (`meds_hydr_lib.f90`), and no copy of that math exists elsewhere (checked).
     - The new record is defined there too, not in `fast_dynamics/plant/meds_plant_types.f90`, so
       `shared/functions` keeps depending only on `shared/base`. `hydro_params_t` then holds one
       record for the leaf and one for the wood.
     - **Recommended:** split `meds_hydr_lib.f90`, which holds four jobs in 428 lines. The storage
       curves (tissue pressure–volume and soil retention) move into
       `shared/functions/meds_water_retention.f90`. Conductance (vulnerability curve, Kirchhoff flux)
       and the root profile stay in `meds_hydr_lib.f90`.
     - The new module uses the vulnerability function from `meds_hydr_lib`, a one-way dependency,
       because embolism release needs it.
     - The split is a pure code move with no output change, in this commit.
     - **The owner confirms the split.**
3. **The physiology change** (`meds_hydr_lib`).
   - **The new curve.** Wood water W(ψ) = living-cell water (unchanged) + the drainable apoplastic
     water × the fraction of conductance retained at ψ + a small bound residual. The retained
     fraction comes from the vulnerability curve the model already has (`wood_psi50`,
     `wood_kexp`).
   - **Capacitance** gains the matching term: the drainable water × d(retained)/dψ.
   - **The leaf curve is unchanged,** since the leaf's apoplastic water is bound: its record simply
     has no drainable part.
   - **Mass → potential.** With two terms there is no closed form. Use a safeguarded Newton
     iteration started from today's closed-form answer. Pass the curve record explicitly, not an
     internal procedure (Phase 1a's pattern), so #325's cost does not come back.
   - **The `rwc_floor` clamp.** It stays only as a guard below the bound residual.
   - **New trait: the wood's bound fraction.** Its default is taken from TFS-Hydro's sapwood residual
     fraction (Christoffersen et al. 2016) for the example PFTs, and checked in the PR.
   - **Memory.** MEDS's loss of conductance has none, so released conduit water refills when the
     potential recovers. Keep that, consistent with how conductance is treated today. Say it plainly
     in `docs/science/plant_hydraulics.md`, and leave a memory for loss of conductance to a later
     change.
4. **Tests.**
   - mass → ψ → mass exact over the whole range, both tissues;
   - capacitance continuous and positive;
   - a collapsed-store fixture whose sub-step count is bounded;
   - `test_plant_hydraulics`'s existing cases, whose expected values move only where the wood is
     dry.
5. **Evaluation** (the outputs change by design):
   - **Runs.** BCI five years, default and calibrated, before against after:
     - the tower statistics (GPP, LE, H, Rnet);
     - predawn and midday ψ;
     - the thrash and sub-step counters;
     - whole-site water at 0 fails;
     - the stand summary.
   - **Line-searched set** (`wood_psi50` −1.24, `leaf_pi0` −2.43): thrash near 0, and no uptake
     bursts.
   - **Calibration.** Re-fitting (about 100 core-hours) is not part of this PR unless the owner asks.
     The PR reports how far the shipped calibration's scores move.
6. **Docs.**
   - `docs/science/plant_hydraulics.md` §2: the new curve and its sources (Appendix A).
   - The CHANGELOG.
   - ROADMAP §4's `rwc_floor` entry.
   - Close #104 with a comment.

### Phase 4 — Moved to a later numerics round (D4)

To pick up then:

- the T8 convergence checklist on the current code (24 h, and the drought test under
  `linear_decline`);
- whether psi_wood's 0.174 MPa at 900 s matters;
- the reopen conditions of #158 and #159;
- #189 item 2's energy half.

### Phase 5 — Runtime consolidation and region threading (#310, #183, F10, #196)

**Decided (owner, 2026-09-30): region runs thread over polygons (option B of §2.7).** The flattened
pool and patch threads inside a region (option A) are not pursued. It needs #325 (Phase 1a) and the
per-thread workspace (Phase 2) first: without #325, polygon threads would queue on the same allocator
lock that slows patch threads today.

**Phase 5 status (2026-09-30).**

- **Done:**
  - F10 (`e77bf62`);
  - polygon threads with R4, R5 and R12 (`b3bd130`);
  - R2, R6, R9 and R11 (`45a16a4`);
  - the O11 comment (`6301093`).

  Every r1 case is bit-identical. The region's 403 files are identical at 1, 10, 20 and 40
  threads.
- **Measured.** 100 cells, one year:

  | threads | wall time |
  |---|---|
  | 1 | 785 s (v0.3.1: 1,012 s) |
  | 10 | 131 s |
  | 20 | 92 s |
  | 40 | 81 s |

  Below the 17–25× expected in §2.7, because the one detail polygon costs about four times an
  ordinary one and sets each month's pace: without it, 40 threads take 53 s (15×). The region's own
  files cost about 2 s. **Next lever:** give a detail polygon the spare threads for its patch loop,
  or schedule it first with a smaller chunk.
- **R10, one shared fast context: deferred, with the reason.** The only part of a polygon's context
  that differs is its acclimated leaf table, and that sits inside `col_config`, which the whole fast
  loop reads. Sharing the context would mean either:
  - copying `col_config` into each polygon step, which #188 removed; or
  - passing the leaf table through four more routines.

  The saving is a few tens of kB per polygon. The stand-in that R10 names (polygon 1's soil curve
  for the region's files) is harmless, because every polygon's curve comes from the same config.
- **Item 5, one run container (rv-R3): proposed, not built.** It changes what callers see, so the
  owner decides two things first:
  - **The step.** One `driver_step` advances every polygon by one slow step: prefetch, a parallel
    loop over polygons, and the I/O phase at a month's end. This keeps the C-API's per-day stepping,
    which the Python examples use. Stepping a region by the day costs nothing measurable (365
    fork-joins a year, against 21 ms per polygon-step).
  - **Failures.** A failed polygon is flagged and skipped from then on, and `driver_step` returns
    the failure. Does `meds_main` stop there, writing the output so far, or finish the run and
    report the failed polygons at the end (the polygon plan's R5)?

  With those answered, the two containers, the two `meds_main` loops and the two step routines
  become one each. Open and finalize keep a site branch (census, restart, checkpoints) and a
  region branch (cells, detail polygons).
- **#196, first step measured: no gain from compiler vectorisation.**
  - `surface_derivs`'s cohort loop does not vectorise: it calls `veg_energy_balance` in another
    module once per cohort, and ifx does not inline across modules without `-ipo`.
  - Targeting the nodes' AVX-512 (`-xHost`, Xeon Gold 6230) gives 24.4 s against 23.9 s on 60 days
    of BCI, serial: no gain. The fast loop is branchy per-cohort work, and libm already chooses its
    AVX-512 routines at run time.
  - `-ipo` needs the LLVM archiver wired into CMake and was not measured.
  - What is left of #196 is the cohort-level tasks, for a site with fewer patches than cores.
- **Also done in this phase:**
  - #312 O7+O8, one serializer (`e50fbaf`, below Phase 6);
  - a stale-value read in slow-only runs (`11e00bf`).

1. **#310 prep PR,** bit-identical:
   - one `open_output_files` in place of three copies (R2);
   - fold the stepper's optional-argument branches (R6), which also mends the dropped
     `latitude_deg`;
   - the location on the polygon (R9);
   - R11 and R12;
   - the O11 comment.

   **This consolidates** three copies of the output set-up sequence and six call branches.
2. **The month loop, simplified before it is threaded.** A failed month is consumed (R4), and a
   month's steps are listed once (R5).
3. **F10:** `met_instant` becomes `pure`, a one-word change, so the compiler checks that the
   forcing sample is safe to call from threads.
4. **Polygon threads (#183 R3).**
   - **The loop.** One `!$omp parallel do schedule(dynamic)` over the region's polygons, once per
     month, in the month's compute phase. The I/O phase stays serial, as it is now. Costliest
     polygons start first, by last month's cost.
   - **Inside it,** each polygon runs its patch loop with one thread. MEDS already passes the patch
     loop its thread count, so this is a single rule.
   - **Settings and build.** In region mode, `[run].n_threads` means polygon threads. The placeholder
     refusal goes. No compilation flag: the build keeps its one OpenMP switch.
   - **Blockers.** Before threading, confirm the polygon plan's blockers B1–B12 against the current
     code: no file I/O and no saved variables inside a polygon's step, and per-polygon output
     buffers.
   - **Gate:** a region run gives identical results at 1, 4 and 8 threads, and a new threaded
     ctest case asserts it. Node use is measured on the 100-cell box at 10/20/40 threads (expected
     17–25× at 40, §2.7).
5. **One run container.** A site run is a region of one polygon (rv-R3). The thread rule then
   becomes one sentence: threads go to the polygons when there are several, and to the patches when
   there is one. **This consolidates** two containers and two program loops into one.
6. **One shared fast context,** with the leaf table per polygon (R10).
7. **The flattened pool is deferred** (owner, 2026-09-30).
8. **#196, after the polygon threads.**
   - First vectorise the cohort loop in `surface_derivs`; `-qopt-report` may show it already is.
   - Then cohort-level tasks for the per-dt_fast leaf and hydraulics solves, when a site has fewer
     patches than cores (BCI: 14 against 40).

   Measure the fusion-off cost law first (GPU evaluation §12.5).

### Phase 6 — Output, forcing and state-layout cleanup

**Phase 6 status (2026-09-30).**

- **Done.** Every r1 case is bit-identical unless noted.
  - O7+O8, one serializer (`6519ef5`).
  - O9: `strict_caps` retired, and the caps are checked at the step (`14a6550`).
  - O10, O12, F9, F11, F12 (`27e1245`).
  - F7+F8 (`06a6b2c`).
  - F10 (`e77bf62`, with Phase 5).
  - The #299 fix (`385344d`): a slow-only run reports its slow rates and reads its fast-loop
    variables as missing. It moves `demography_30yr` only.
  - A stale-value read found on the way (`11e00bf`).
- **Not started, and why:**
  - **O6 + #270:** one forcing-echo table and a patch axis on the FAST tier. #270 is a new
    output feature, and it touches the same staging lines as the fast loop.
  - **#275:** the skin temperature, a science-output change.
  - **F6:** splitting the 1,250-line reader into sources.

  Each is a PR of its own size; this one is large enough to review as it stands.
- **#299 caveat:** litterfall stays 0 in a run with soil carbon off, which does not accumulate
  litter at all. Whether it should (it relates to #316) is a separate decision.

- **#312.**
  - O7+O8: one serializer for site and region files. It can go at any time.
  - Then O6 + #270 together: one table of forcing echoes, and a patch axis on the FAST tier.
  - Then #275's skin temperature.
  - Then the #299 half with O9: a "needs the fast loop" flag gives missing values, and
    `cohort_max` is enforced early.
  - O10 and O12 separately.

  **This consolidates:** two serializers into one, and three lists of the forcing echo into one.
  O6/#270/#275 edit the same staging lines as Phase 5, so they land after it, never in parallel.
- **#311.**
  - F6: split the reader into sources.
  - F7: one hourly-cadence constant and one "seconds since" parser.
  - F8: resolve field indices at open, and reuse the next bracket record when the bracket slides,
    which halves the reads.
  - F9: delete dead code; `wind_log_profile`'s test goes with it.
  - F11 and F12.

### 3.7 New items found by this review

Per D6 these are not filed as issues; each is named in the commit and PR message that handles it.

| id | item | where | phase |
|---|---|---|---|
| N-1 | Forcing files written with 1×1 chunks | `scripts/forcing_common/meds_forcing_file.py` | 1b |
| N-2 | The corrector's hydraulics counts are discarded | `meds_fast_be_stage.f90:394, 435` | 1d |
| N-3 | Stale #162 figures and warning text | four files (§2.5) | 1e |
| N-4 | A selection sort (cost ∝ cohorts²) every sub-step | `apply_rt_forcing` | 2 |
| N-5 | Census fusion reorders all cohorts after every fused pair | `patch_fuse_pass` → `rebuild_csr` | optional: 3.9 s once per census start; reorder once per pass |
| N-6 | `pow` in the soil hydraulic functions | `meds_hydr_lib` | optional: 6.5% of CPU; precomputed exponents move rounding everywhere |
| N-7 | The daily tissue-water reconcile ignores a PFT's own curve | `meds_fast_reconcile` | 3 (fixed, `10d6482`) |
| N-8 | A diagnostic reader leaves its slots unset when its block has no entries, so a slow-only run read stale values | `meds_site_diag_types` | 6 (fixed, `11e00bf`) |
| N-9 | A detail polygon (hourly files, per-cohort diagnostics) costs about four times an ordinary one and sets a threaded region's pace | `meds_region` | 5 (measured; the lever is noted in the Phase 5 status) |

---

## 4. Gates and verification, per PR

- The ctest suite on ifx Release, ifx Debug and gfortran Release; the Python tests with `MEDS_LIB`.
  One gfortran profile for #325, because the mechanism is ifx-specific.
- r1 `compare_runs.py` against the Phase 0 references, with the movement stated (D1).
- B1–B5 timings for every performance PR, each run alone on an `R128C40` node, in the PR body.
- `test_fast_loop` thread invariance on an OpenMP build.
- nvfortran is not on cbsuxu; each PR says the portability check was not run.
- Each PR states what it consolidated or removed (§0.1).

## 5. Decisions

| # | question | answer |
|---|---|---|
| D1 | rounding-level movement | **accepted** in every PR, stated with its size |
| D2 | forcing chunking | **fix the shared writer**, and print a note that says why it is printed (Phase 1b) |
| D3 | #104 option | **(b), in this sweep** (owner, 2026-09-30): the wood's apoplastic water drains as its conduits embolise; its own physiology PR (Phase 3) |
| D4 | #158, #159, #167 | **closed 2026-09-30**; the numerical scheme is revisited later |
| D5 | #146/#190 | **option B** (owner, 2026-09-30): one field list and one rules table, done inside Phase 2 |
| D6 | new issues | **none filed**; new items noted in commit and PR messages |
| D7 | releases | **none**; every commit and PR updates `beta` |
| C1 | where the forcing writer lives | **one script, `scripts/meds_forcing_file.py`**; tower-only conversions into `prepare_flux_tower/` (§5.1) |
| C3 | region threading shape | **threads over polygons** (option B of §2.7; owner, 2026-09-30). The flattened pool is deferred, and patch threads inside a region are not pursued. No compilation flag |

### 5.1 Where the forcing writer lives (owner's comment 1)

`meds_forcing_file.py` has three users:
- `prepare_era5/make_forcing_file.py`, which cuts site files from the ERA5-Land archive, for
  example `data/forcing/ithaca_forcing.nc`;
- `prepare_flux_tower/make_tower_forcing.py`;
- `calibrate_fast/tests/test_smoke.py`.

It was made shared in #320 (flux-tower plan §6) so the two tools write one format. Folding it into
`prepare_flux_tower` would make the ERA5 tool import from the flux-tower example's folder.

The file holds two things:
- **the writer** (`write_forcing_file`, `check_complete`, about 55 lines), which all three use;
- **Python copies of the model's conversions** (humidity, pressure, solar geometry, the longwave
  synthesis, about 190 lines), which only the tower tool uses.

Three layouts are possible:

| option | what moves | effect |
|---|---|---|
| (i) keep as is | nothing | the chunking fix is made once and every tool gets it |
| (ii) the conversions move into `prepare_flux_tower` | about 190 lines | `forcing_common` keeps only the writer; each folder holds only what it uses |
| (iii) everything moves into one tool's folder | the whole file | the other tool imports across folders |

**Decided by the owner (2026-09-30):** no subfolder. The writer is one script,
`scripts/meds_forcing_file.py`, and the tower-only conversions move into the tower example's folder,
`scripts/prepare_flux_tower/`. Phase 1b carries it out.

## 6. Expected gains (estimates; Phase 1–2 measurements replace them)

| workload | today | after Phase 1 | after Phase 2 | after Phase 5 |
|---|---|---|---|---|
| BCI 60 d, best thread count | 30.0 s (4 threads) | about 12–15 s (8 threads) | about 10–13 s | same |
| BCI 5 yr, best thread count | 5:17 (4 threads) | about 2 min (8–14 threads; 14 patches cap it) | −10–15% | same (cohort tasks may add some) |
| calibration trial (10 d, 1 thread) | 8.5 s median | about 5 s | −10–15% | same |
| 100-cell region, 1 year | 1,012 s serial (71 s as 100 separate runs, 40 at a time) | serial | serial | about 55–70 s threaded on one node |

---

## Appendix A — #104: the three options, explained (D3)

**What goes wrong.**
- **The two stores.** Each cohort's leaf and wood hold water, and the model tracks the mass `W` of
  each. The water potential ψ, which drives flow between soil, wood, leaf and air, is read off a
  pressure–volume curve from `W`.
- **What the curve can represent.** It splits a store's water into two parts:
  - a fixed apoplastic part, 20% of the saturated water for wood (`wood_apoplast_frac`), which the
    curve treats as never draining;
  - the living-cell (symplastic) part, whose relative content R sets ψ. Below turgor loss,
    ψ = π₀/R.

  The curve therefore cannot describe a store holding less water than its apoplastic part.
- **What happens below it.** `psi_from_water_content` computes a negative R, clamps it to
  `rwc_floor` = 1e-4, and returns ψ = π₀/1e-4, about −10⁴ MPa (`meds_hydr_lib.f90:221-231`, `:288`).
  That value is a marker, not a potential: xylem is fully cavitated by −10 MPa.
- **How a store gets there** (*code reading*; Phase 3 confirms or refutes it). When the soil cannot
  supply what the plant asks for, the pre-pass cuts only the wood's uptake. Sap flow from wood to
  leaf is not cut (`meds_fast_frozen.f90:536-539`). The stage then updates the wood's mass by a
  straight line, floored at zero (`meds_fast_be_stage.f90:470-490`), so one step can drain the wood
  below its apoplastic part.
- **What the next step does.** The solver sees a −10⁴ MPa wood against, say, −2 MPa soil. It asks
  for a huge refill, sub-steps up to about 137 times per cohort, and delivers the 3–8 kg m⁻² uptake
  bursts in single steps. It also prices the refill from the curve (`W(ψ)` differences), and the
  curve believes the apoplastic water is still there. So the refill falls short by that gap, and the
  store can stay below the curve for several steps. That is the consecutive-step thrash.

**What the "fixed" water is, and whether it drains in nature** (literature search, 2026-09-30):

1. **Where the fraction comes from.** A pressure–volume curve plots 1/ψ against the water a drying
   tissue has lost. The straight part below turgor loss is extended to 1/ψ = 0, that is, to an
   infinitely negative ψ, and the relative water content left there is called the **apoplastic
   fraction**. It is water outside the living cells: in the cell walls, held by matric forces, and in
   the xylem ([Bartlett et al. 2012](https://adinet.ahacentre.org/repository/detail/the_determinants_of_leaf_turgor_loss_point_and_prediction_of_drought_tolerance_of_species_and_biomes_a_global_meta_analysis_20120322);
   [PrometheusWiki protocol](https://prometheusprotocols.net/function/water-relations/pressure-volume-curves/leaf-pressure-volume-curve-parameters/);
   [Neufeld, water-potential components](https://appstate.edu/~neufeldhs/pltphys/waterpotentialcomponents.htm)).
   MEDS uses this standard curve: ψ = π₀/R goes to −∞ as the living-cell water R goes to 0, leaving
   the apoplastic water behind.
2. **In leaves it is mostly bound water, and treating it as fixed is sound.**
   - A leaf's apoplastic water is mostly cell-wall water. It leaves only when tissue air-dries, far
     below any potential a living leaf survives.
   - TFS-Hydro uses Bartlett's relation RWC_r = 0.01·ε + 0.17 for the leaf's residual fraction
     ([Christoffersen et al. 2016, GMD](https://gmd.copernicus.org/articles/9/4227/2016/)). MEDS's
     leaf default ε = 12 MPa gives 0.29; MEDS uses 0.30.
   - What real leaves keep losing after their stomata shut is water through the cuticle and leaky
     stomata (the minimum conductance g_min). In a long drought that loss drives them to lethal
     dehydration ([Duursma et al. 2019, New Phytologist](https://hal.inrae.fr/hal-02627387v1)). That
     water comes from the living cells, which the curve already lets drain.
3. **In wood, most of the "fixed" water does drain.** Sapwood water-release curves have three phases
   ([Tyree & Yang 1990](https://link.springer.com/article/10.1007/BF02411394), as summarised in
   [PMC6446748](https://pmc.ncbi.nlm.nih.gov/articles/PMC6446748)):
   - (I) capillary water in fibre and vessel cavities and intercellular spaces, released between 0
     and about −0.5 MPa;
   - (II) elastic release from the living cells;
   - (III) water freed from the conduits as they cavitate, at moderate to low potentials.

   **Cavitation release is a real drought buffer.** It hands water to the transpiration stream and
   briefly raises the xylem potential, which matters most in large trees and in drought
   ([Hölttä et al. 2009](https://hal.archives-ouvertes.fr/hal-01189354)).

   **Models already treat it as drainable:**
   - Christoffersen et al. 2016 renamed the fixed part the **residual fraction**, "in light of the
     considerable amount of water released when vessels embolize in stems", and gave sapwood a
     capillary phase.
   - SurEau ([Cochard et al. 2021](https://hal.inrae.fr/hal-03269429)) describes each organ's
     apoplasm by its vulnerability curve, so that water is released as embolism proceeds. It follows
     tissue desiccation beyond stomatal closure.

   **Wood science draws the same line.** Water bound in the cell walls, up to the fibre-saturation
   point of about 30% of dry mass, is held; the free water in the cell cavities is not
   ([Oregon State Extension EM 8600](https://extension.oregonstate.edu/catalog/pub/em-8600-wood-moisture-relationships)).
4. **For MEDS** (*points 2 and 3 of this list were corrected on 2026-09-30, during Phase 3: with
   MEDS's saturated wood water, the cell-wall water alone exceeds the 20% held fixed, and TFS-Hydro's
   residual is larger still. See "Phase 3 status".*):
   - The leaf's 30% fixed water is well founded.
   - The wood's 20% is not: it lumps drainable conduit and capillary water with the small amount
     that is truly bound.
   - MEDS's wood curve has neither the capillary phase nor cavitation release.

**The three options.**

**(a) Bound the potential.**
- **What changes.** A store below its apoplastic part reports a finite floor, ψ_min (about −10 to
  −20 MPa), instead of −10⁴ MPa. Alternatively, raise `rwc_floor` so that π₀/rwc_floor equals
  ψ_min.
- **Effect.** The refill starts from a realistic gradient, so the stiff burst and most of the
  sub-stepping go away.
- **Cost.** About 10 lines.
- **What it leaves.**
  - It adds a parameter (ψ_min) and keeps a clamp: one special case replaces another.
  - The mismatch between the mass and the curve stays, so the refill is still priced wrong, only
    less violently.
- **Who it touches.** Every reader of ψ in that state: stomata (`linear_decline` already shuts them
  fully below 2ψ_tlp), conductance loss, and any mortality rule reading ψ. The default BCI run is
  unchanged: it never gets there.

**(b) Let the wood's apoplastic water drain as its conduits embolise** (revised after the literature
search above).
- **What changes.**
  - Wood water becomes: living-cell water (unchanged) + apoplastic water × (1 − loss of
    conductance at ψ) + a small bound residual. The loss of conductance comes from the vulnerability
    curve the model already has (`wood_psi50`, `wood_kexp`). This is SurEau's approach, and the
    capillary phase of TFS-Hydro is the same idea near saturation.
  - The leaf curve is unchanged.
- **Effect.**
  - Every wood mass down to the small bound residual gets a finite potential: roughly where the
    conduits are embolised, a few times ψ50 (about −5 to −10 MPa when ψ50 is −2.5 MPa). So for wood
    the −10⁴ MPa marker and the short
    refills go away, together with most of the reason for the clamp.
  - Embolism becomes a source of stored water in drought, as observed.
- **What it needs.**
  - One new trait: the bound fraction, which is small for wood.
  - A decision on memory. MEDS's loss of conductance is a function of today's potential, with no
    memory, so the released conduit water would refill as soon as the potential recovers. Real
    embolised conduits refill slowly or not at all. Either accept that simplification, or give loss
    of conductance a memory, a separate and larger change.
- **Who it touches.** Every run whose wood dries toward ψ50, not only collapsed stores. With the
  smooth vulnerability curve, some conduit water is already released at modest potentials (at
  kexp = 2, about 14% of it at ψ = 0.4·ψ50). So it moves ordinary drought behaviour, including the
  calibrated BCI runs.
- **Kind of change.** A science change: evaluate it on the BCI five-year runs and check the
  calibration, not as a performance fix.
- **Cost.** About 50–100 lines in `meds_hydr_lib` plus round-trip tests over the whole range.

**(c) Stop the drain where it starts.**
- **What changes.** When the soil throttles uptake, the hydraulics solver takes the soil's
  deliverable uptake as its root boundary (a specified-flux, "Neumann", root boundary), instead of
  the soil potential. Sap flow then adapts to the water that really arrives, and the wood is not
  over-drained.
- **Blend.** To stay continuous, blend the two boundaries by how hard the soil throttles (weight
  1 − scale), so an unthrottled step (scale = 1) is exactly today's.
- **History.** This boundary was built and measured once. In ordinary conditions it moved the
  error from the wood to the leaf (psi_leaf error 4.6e-3 → 0.16 MPa), so it was not landed. The
  patch file (`neumann_root_bc.patch`) was in an old scratch area and is gone, so it would have to
  be rebuilt.
- **Effect.** It removes the cause on this path, but not the marker. A store that reaches the floor
  another way, such as a large swing in canopy-air humidity within a step, still returns −10⁴ MPa.
- **Cost and risk.**
  - About 60–120 lines, and two boundary conditions inside one solver: a special case added.
  - The leaf, which drives stomata, takes more of the error when the soil throttles.
- **Default run.** Unchanged where the soil does not throttle.

**Decided (owner, 2026-09-30): (b), in this sweep, as its own physiology PR (Phase 3).** The
recommendation it was chosen against is kept below for the record.

**Recommendation (before the decision).**
- **The literature supports (b)** for wood. It removes the marker where it arises and adds a buffer
  the real plant has.
- **But (b) changes ordinary drought behaviour,** so it is a physiology change to evaluate on its
  own, arguably in the later numerics and physiology round, not a performance fix.
- **For this sweep, (a)** stops the cost and the bursts without changing any run that never reaches
  collapse, and (b) can later replace it.
- **(c)** adds the most complexity for the least coverage.
- **Either way,** Phase 3's instrument goes first, to confirm the drain path.

## Appendix B — #146 and #190: what the options mean (D5)

**The problem.**
- **The type.** The fast integrator's state is one derived type, `column_state_t`, with 11 fields:
  canopy-air enthalpy, humidity and CO₂; soil energy and water per layer; the pond and its
  enthalpy; and leaf and wood water and surface films per cohort.
- **The hand-written routines.** About 15 routines each list those fields by hand: copy
  (`state_init`), `y + a·k` (`state_axpy`), accumulate, extrapolate, the error difference, subtract,
  zero, the clamps, the commit (`unpack_column_state`), plus the pack in `build_column_frozen` and
  the error norm.
- **Per-field exceptions.** Several routines treat some fields specially, each written out where it
  applies:
  - the pond has no tendency, so it is copied, not stepped;
  - films are left out of the error norm;
  - masses are zeroed in the error.
- **The hazard.** Add a 12th field and forget it in one routine: the code compiles and the tests
  pass. The budgets still close, because a field left at its default reads as a consistent zero,
  and the step controller simply sees a smaller error.
- **Today's guard.** `test_state_combinators` catches an existing field dropped from a routine,
  which is the failure that has actually happened. It does not catch a new field wired nowhere.
  There are no live defects of this kind.
- **Performance.** Phase 2 removes the combinators' allocations either way, so none of the options
  below is about speed.

**The options.**

| | what it is | what a new field costs | omission caught? | readability | size |
|---|---|---|---|---|---|
| **A. Leave it** (close #146) | as today | edit about 15 routines | existing fields only (the test) | unchanged | 0 |
| **B. One list, in one place** | Beside the type, one pair of routines copies the fields into and out of a plain array, and one small table names each field's rules (steps? in the error norm? clamped?). Every combinator becomes a one-line array operation: `y%v = x%v + a*k%v`. | edit the type and the two adjacent routines | not by the compiler; a test checks that the pair is exact and that the array length matches the type | **better**: 15 hand lists become 2, and the per-field exceptions become one table instead of comments scattered across routines | about 300–400 lines, net shorter |
| **C. A generator** | A small table file (field names plus their rules) and a Python script that writes the type's body and every combinator as a Fortran include file. The generated file is committed, and a ctest fails if it no longer matches the table (as `io_config_example` does for `meds_io_config.toml`). | edit the table, re-run the script | yes: generated code cannot forget a field | **worse** for reading and grepping the combinators: they live in a generated file, and a contributor must know to edit the table | about 500 lines plus a 200-line script |
| **D. Fully packed state** | Every field becomes an index into one array (`y%v(I_CAS_ENTHALPY)`) throughout the physics code | nothing to forget | yes | **worse**: about 470 references in 9 source files lose their names | 1,200–1,600 lines |

**#190** (the cohort gather type `column_cohort_t`) has the same problem: its per-field fusion and
scaling rules are already in two routines (`fuse_cohort_fast_state`, `scale_cohort_ground_fields`).
The table from option B would be the natural home for them if it is chosen, and otherwise #190 can
close with #146.

**Does B mean wrapping the fast-loop state in a structure and passing all of it to every function?**
No.
- **It is already one structure,** `column_state_t`, and only the integrator's bookkeeping uses it:
  8 files, all in `src/fast_dynamics/numerics/`.
- **The physics routines never receive it, and B keeps it that way.**
  - `surface_derivs` gets a smaller canopy-air record (`surface_state_t`) and its coefficient
    records.
  - Hydraulics, gas exchange, soil water and radiation take plain arrays and numbers.
- **B changes only the inside of the ~15 bookkeeping routines.** Instead of each listing the 11
  fields by hand, they call one pair of routines written once beside the type, and do their arithmetic
  on the flat array:

  ```fortran
  ! beside column_state_t: the one place that lists the fields, and the one table of their rules
  pure subroutine state_to_array(y, n, nsl, values)     ! fields -> one array
  pure subroutine array_to_state(values, n, nsl, y)     ! one array -> fields
  !   rules%steps(:)     does this entry change within a stage? (the pond does not)
  !   rules%in_error(:)  does it count in the step-size error? (the films do not)

  ! a bookkeeping routine, e.g. y_stage = y + a*k, then reads:
  call state_to_array(y, n, nsl, y_values)
  call tend_to_array(k, n, nsl, k_values)
  y_values = y_values + a * k_values * rules%steps
  call array_to_state(y_values, n, nsl, y_stage)
  ```
- **No cost from passing the whole structure.** Fortran passes a derived type by reference, so
  handing a routine the whole structure copies nothing. `use`-ing the module only makes the type's
  definition visible. The flat arrays would live in the Phase 2 workspace, so B adds no allocations.
- **Rounding.** The step-size error today adds the fields in one particular order. B keeps that
  order, or accepts a rounding-level change (D1).

**Decided (owner, 2026-09-30): B, inside Phase 2.** The reasoning: #195 rewrites the same routines anyway, and it is the
option that consolidates and removes special cases (§0.1), with no extra tool. C is the only one
that gives a compile-time guarantee, at a readability cost. If neither is wanted, close both (A).
