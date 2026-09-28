# MEDS code review, 2026-09-27 — `beta` → `main` for v0.3.0

**Status:** review record, read-only. No source changed by this document. Reviewer: Claude (Fable 5.1)
with seven parallel passes (runtime, forcing, output + slow dynamics, Python tools, documentation,
issue triage, build + test); every finding ranked here was re-verified by hand against the tree at
`e33cb56` (`beta`, 62 commits and 82 files ahead of `main` = `v0.2.2`, `92cad44`). Line numbers
refer to `e33cb56`.

**Executed 2026-09-28 in #314:** phases 0–5 of §9, and the archive audit of §10. R2 (Phase 1.7) and
P3 (Phase 4.7) are deferred to #310 and #313; the integrator plan stays live until #162 and E5 are
re-homed. At the owner's request the same PR also removes the `[io]` shim (#309), makes the Debug
suite green (#308) and does the `prepare_era5` dedup (#313, which includes P3). The after-merge
items of §8 item 6 that remain are #310–#312.

The three questions asked: (1) can changes be consolidated under a generic principle and special
cases removed; (2) is the documentation consistent with the code; (3) are the new modules
modularized and readable. In parallel: which open issues does `beta` close, and which should be
fixed before it merges.

---

## 0. Verdict and the recommended pre-merge set

**The architecture is sound and the release is close.** The R1/R2 split is clean where it matters:
one shared `met_source_t` with per-polygon `met_cursor_t`; `output_files_t` shared and
`output_buffers_t` per polygon, with netCDF ids touched only by the serializer; `meds_polygon` owns
one polygon's state and step, and no module-level mutable state was introduced anywhere in the diff
(the only `save` in the fast loop, `meds_fast_dynamics.f90:967-968`, predates `beta`). The `forcing`
output group rides the ordinary registry. The Python ↔ Fortran archive contract (names, units,
epoch, chunking, fill, latitude order, humidity formula) checks out byte for byte.
`meds_io_config.toml` lists all 253 registered variables (#241 has not regressed).
`order_of_processes.md` §2/§5/§6 match `polygon_step` line for line. On a clean export, ifx Release
passes **53/53** and the 4-cell regional smoke run completes in 5 s with the expected file shape; the
three ifx Debug failures are identical on `main` (§6).

**Two latent crashes and one wrong-ledger bug should not ship**, and a handful of small fixes belong
at a release boundary. Ranked:

| # | Item | Where | Effort |
|---|---|---|---|
| 1 | Archive + recycling **crashes at the year wrap for any `recycle_start` not at 01:00** (§2 F1). Every config in the tree uses `01:00`, which is why it never fired. Fix the prefetch, or reject other anchors in `validate_config` with a message. | `meds_met_driver.f90:300-316, 434-439, 710-717, 905-918` | S |
| 2 | **#290, the −0.3…−4 W/m² energy leak, is a ledger omission, not physics**: under `bottom_bc = "dirichlet"` the solve applies `−g_deep·(T_bot − deep_temp)` but every ledger books the bottom face as `geothermal ≡ 0`. On `main` since v0.2.1 (#267), so not a regression — but v0.3.0's own examples print a budget warning at every step. Bookkeeping-only fix (§7 B1). | `meds_soil_energy.f90:126`, `meds_fast_frozen.f90:346`, `meds_fast_be_stage.f90:187,200,203`, `meds_fast_rk45.f90:823` | S code + a verification run |
| 3 | **#306**: canopy-air depth is 20 m on day 1 of every run; `refresh_canopy_depth` is reachable only from the slow step. Exposed by #305 (the forcing is now moved to `can_depth`). One call at the end of `polygon_prepare`. | `meds_slow_dynamics.f90:61,91,108`; `meds_polygon.f90:99-131` | S |
| 4 | Region I/O phase on a failed polygon **error-stops instead of writing what closed** (§3 O1). | `meds_output_manager.f90:70-74`; `meds_region.f90:228-236` | S |
| 5 | **CHANGELOG**: the #294/#296 and #295 `### Fixed` entries are filed under `[0.2.1]`, `[Unreleased]` has no `Fixed`; two `[Unreleased]` sentences say `prep_era5land_forcing.py` "stays" / "now reads box files" while the same section removes it. Needs a `## [0.3.0] — date` heading. | `CHANGELOG.md:394-421, 115, 302, 307-308` | S |
| 6 | **Docs that describe the previous shape** (§5): `src/README.md` `main/` row and counts; `forcing.md` API text (`met_advance(drv, now)`, `ensure_month`); `configuration.md` `[site]` sentence; `.claude/rules/output.md` and `state-demography.md` say the tick runs after fuse-fission; `diagnostics.md` inventory counts; a dozen source comments narrating history (rule 1). | see §5 | S–M |
| 7 | Small code hygiene found on the way: `met_open` leaks `src%ncid` on one rejection path; `met_instant` carries a duplicated LW/ρ block in its CONST early return; `polygon_prepare` does not reset `restructure_pending`; `detail_polygons` validated after building every polygon; `make_forcing_file.py` re-floors wind and drops `u10/v10`; longitude-axis docstrings say −180…179.9 for an axis that is −179.9…180. | §1 R1, R8; §2 F2, F5; §4 P1, P4 | S each |
| 8 | Two open bugs whose cheapest half fits here: **#299** (slow-only runs weight patch rows with `w = 0` → one line) and **#275** (drop/relabel `ground_temp_site`, relabel the four variances — an output-API change that "wants a release boundary"). | §7 B4, B5 | S |
| 9 | **FAST-tier stamp**: fast records are stamped at the sample instant (`forcing_sample_frac = 0.5` → 07:30 into a 15-min window) under a `long_name` of "period start", while D/M/Y now stamp by window start (#294). Either stage `step_start + (isub−1)·dt_fast` or fix the label. | `meds_fast_dynamics.f90:490-498`; `meds_output_stream.f90:193-203` | S |
| 10 | **gfortran `region` segfault** (new code, mechanism verified by experiment, §6.2): `reg%poly(:)%out_bufs` — a strided component section of a type with nested allocatable components — is passed `intent(inout)` to `output_serialize_region` and `output_region_close`; gfortran's array-temporary copy-out leaves the buffers dangling. Give the region a contiguous `bufs(:)`, and add the construct to CLAUDE.md's trap list. | `meds_region.f90:254, 310`; `meds_output_manager.f90:64-67, 83` | S/M |
| 11 | **gfortran cannot build at all**: six lines over column 132 (pre-existing since v0.2.2, a hard error in gfortran 11.5 and 15.2). Wrap them; `main` has the same defect. | §6.1 | S |

**Time-box:** #298 (restart not exact) to persisting `patch%adapt_dt_last` plus one `ckpt_2yr`
comparison; if that does not close it, record what remains and defer to R4.

**Defer to issues** (the consolidation work, §1–§4 `[AFTER-MERGE]`): the one-container / one-loop
shape for site and region; the region serializer that duplicates the site serializer; the forcing
echo enumerated three times; `meds_met_driver.f90` as three sources in one module; the optional-argument
fan-out in `advance_one_step`. Most of these are best done with R3, which rewrites the loop anyway.

### Cross-cutting themes (criterion 1)

- **"A site run is a region of one polygon" is met only at `polygon_step`.** Above it there are
  two containers, two program loops, three copies of the output-file-set sequence, two serializers,
  and (in the C API) a third copy of the calendar-boundary sequence. §1 R2–R3, §3 O7, O11.
- **The same list of forcing quantities is spelled out three times** (`fast_sample_t` fields,
  `PY_*` block, registry rows + `extract` cases): adding one variable is ~9 edits in 5 files. §3 O6.
- **Backend branches instead of a source interface**: ten `backend ==` sites in the met driver. §2 F6.
- **Literals and tables in more than one place**: the hourly cadence `3600` in five places; the
  ERA5-Land variable table, unit factors and de-accumulation rule once in `era5land_common.py`
  and again in `make_forcing_file.py`; the wind floor in Python and Fortran. §2 F7, §4 P3, P4.
- **Comments that narrate what the code used to do** recur in every slice (rule 1). §1 R7, §2 F4, §3 O3, O5.
- **Portability is asserted only for ifx.** gfortran has not built either branch since the six
  over-length lines went in, and the one gfortran-only crash is in the new region code (§6.1–6.2).

---

## 1. Runtime — `meds_polygon`, `meds_region`, `meds_driver`, `meds_stepper`, C API

**As built.** `meds_polygon_t` holds the site, fast context, forcing cursor, output buffers, budgets,
ledger, seam, status and the pending-restructure flags; `polygon_prepare` / `polygon_step` /
`polygon_report`. `polygon_step` is one slow step with no file work: the restructuring owed from the
last boundary, acclimation, `advance_one_step`, the fast-tier replay and the slower-tier tick into
record queues, then the NaN and soil-carbon guards. `meds_driver` (`meds_run_t`) is the site run;
`driver_step` = `met_prefetch` + `polygon_step` + `driver_io_phase` at month ends. `meds_region`
(`meds_region_t`) selects the box, builds one `met_source_t` for every cell and N bare-ground
polygons, and runs `region_step_month` = one prefetch, polygon-outer / step-inner compute, then the
I/O phase. `meds_stepper` lost the cadence flags: `advance_one_step` and a new `advance_boundary`.
`meds_main` peeks `run.mode` and dispatches to one of two copies of open / loop / finalize.

Config keys `run.mode`, `region.box_nwse`, `region.land_fraction_min`, `region.detail_polygons` are
parsed (`meds_config_io.f90:707-734`) = documented (`configuration.md:195-230`) = present in
`meds_config_main.toml:493-502`, with matching defaults. No netCDF call inside the polygon loop.

**R1 `[BEFORE-MERGE]` A reused `meds_run_t` carries the previous run's pending restructuring — S.**
`polygon_step` leaves `poly%restructure_pending = is_new_month` (`meds_polygon.f90:205`); only the
restart branch resets it (`meds_io.f90:415`, via `meds_driver.f90:129`); `polygon_prepare` resets
budgets, ledger, seam and counters but not these flags (`:107-112`). The C API hands freed registry
slots out again (`meds_c_api_run.f90:84-90,123-128`) and `test_region` reuses one `run` for three
runs (`test/test_region.f90:83-91`); a run whose `end_time` falls on the 1st therefore leaves a
stale `.true.` for the next bare-ground or census run in the process, whose first step then calls
`advance_boundary` on a stand that owes nothing. Not exercised today (both tests end mid-month).
Fix: reset both flags in the counter-reset block at `meds_driver.f90:114`, before the `select case`,
and add a month-boundary `end_time` to the slot-reuse test.

**R2 `[BEFORE-MERGE]` "Build an output file set" is written three times; the region imports helpers from its peer — S/M.**
`manager_setup → manager_set_soil_params → apply_io_overrides → manager_finalize →
manager_alloc_buffers → activate_site_diag` at `meds_driver.f90:214-226`, `meds_region.f90:139-153`
(plus `manager_restrict_region`) and `:154-164` (detail polygons). To reach `apply_io_overrides` and
`ensure_output_dir` the region driver `use`s the site driver (`meds_region.f90:44`), so the two are
no longer peers over `meds_polygon`. One `open_output_files(files, cfg, prefix, soil,
restrict_region, verbose)` beside `polygon_prepare` (or a small `meds_run_setup` module) makes the
three sites one call each.

**R3 `[AFTER-MERGE]` Two parallel run containers and two program loops — M/L, with R3.**
`meds_run_t` (`meds_driver.f90:65-75`) and `meds_region_t` (`meds_region.f90:50-59`) hold the same
fields with `poly` scalar against `poly(:)`; `step_days` derived twice (`:194` / `:126`), the area
verdict twice (`:322` / `:293`), `met_open` and its message twice, and `meds_main.f90:62-95` holds two
copies of open → loop → status → finalize. Suggested shape: one container with `poly(:)`;
`[run].mode = "site"` builds one polygon from `[site]` (the restart/census branch becomes a
`polygon_init`); `driver_step` stays as the step-granular loop the C API needs over `poly(1)`,
`region_step_month` as the polygon-outer loop R3 threads; `meds_main` has one loop.

**R4 `[AFTER-MERGE]` A failed month is not consumed — M, with R5.** On a non-OK status
`region_step_month` returns before `reg%now = clk` (`meds_region.f90:228-237`), so a caller that
continues re-steps polygons `1..p-1` through the same month (double-integrated tiers, restructuring
fired twice). Latent: `meds_main` and `test_region` stop on any error. R5's failure isolation should
flag the polygon, skip it, and always advance `reg%now`.

**R5 `[AFTER-MERGE]` The month's steps are enumerated twice — S/M.** `meds_region.f90:205-219` walks
the month to prefetch and assert `n_loads`/`carry_rec` do not move, then `:222-227` walks it again.
The loaded month can change inside a model month only at the recycle wrap; with `dt_slow = 1 d` from
midnight already enforced for the archive (`meds_config.f90:783-784`), "in region mode the recycle
window starts on the 1st" is a `validate_config` rule, after which one enumeration suffices.

**R6 `[AFTER-MERGE]` Optional-argument fan-out in `advance_one_step`; one branch drops the latitude — S/M.**
`meds_stepper.f90:65-77, 86-109` spell the same two calls 2 + 4 times to cope with absent optionals,
and `polygon_step` adds a seventh (`meds_polygon.f90:249-267`) although `poly%met_cur` and
`poly%fast_ctx` always exist. The branch at `:96-97` omits `latitude_deg` while its three siblings
pass it — reachable only from tests today, but exactly the error this shape invites. Make
`step_start` required, hoist `rho_air`, pass the remaining optionals straight through.

**R7 `[BEFORE-MERGE]` Documentation describes the previous shape — S.**
- `src/README.md:3,16,70`: `main/` lists `meds_stepper`, `meds_driver`, `meds_main`; `meds_polygon`,
  `meds_region`, `meds_region_opts` absent; "86 modules, 19 CMake libraries" stale (the file is not in
  the diff at all).
- `CMakeLists.txt:352-357`: "This is meds_main's former body … owns … the output manager" for a target
  that is now three files and two drivers.
- `meds_driver.f90:3-11` (the 400-line-program history), `:88-90`, `:212-213` ("the step drains them" —
  the I/O phase does), `:238-240` ("drain the FAST tier … state checkpoint" — both moved out).
- `meds_polygon.f90:43-44`: "the conditions the model used to `error stop` on" (rule 1).
- `docs/science/forcing.md:600` names `ensure_month`, which no longer exists.
- `meds_c_api_run.f90:93` lists statuses 0/1/2 but not 4.

**R8 `[BEFORE-MERGE]` `region_open` validates `detail_polygons` after building every polygon, and leaves the source open — S.**
The ids are known at `meds_region.f90:115`; the check sits at `:129-135`, after N `init_bare_ground`
+ `polygon_prepare` calls, and returns `ok = .false.` without `met_close`. Move it to right after the
id loop.

**R9 `[AFTER-MERGE]` The polygon's location is not on the polygon — S.** It sits in `met_cur` (set
only when forcing is on, `meds_polygon.f90:116-119`); `polygon_prepare` takes four location scalars;
`region_finalize` reads lat/lon back through `poly%cell` (`meds_region.f90:298`). A small location
component on `meds_polygon_t` lets everything read one place.

**R10 `[AFTER-MERGE]` Polygon 1 stands in for the region; N identical fast contexts — M.**
`manager_set_soil_params(reg%out_files, reg%poly(1)%fast_ctx%col_config%soil)` (`meds_region.f90:142`)
because the column parameters depend on `cfg` alone; `build_fast_context` runs once per polygon
(`meds_polygon.f90:115`) although only the acclimated leaf table is mutated per polygon (`:360`).
The plan's §10.3.1 lists this as the memory item to measure before R3.

**R11 `[AFTER-MERGE]` `activate_site_diag` is not additive — S.** It recomputes the need flags from one
file set (`meds_output_registry.f90:1271-1293`); `meds_region.f90:164-166` calls it for the detail set
*or* the region set, correct only because `manager_restrict_region` only disables variables.

**R12 `[AFTER-MERGE]` Dead or write-only state — S.** `meds_run_t%is_open` written, never read;
`meds_region_t%istep`, `%prev` only written (`meds_region.f90:237`).

Nits: `DRIVER_*` status names now belong to `meds_polygon`; `polygon_step` prints the soil-carbon
error from inside the compute phase (`meds_polygon.f90:229-236`, against B7); `is_new_day` (`:195`)
is always true with `step_days >= 1`; `test_region` hardcodes `region-D-202101.nc` (`:96`) and does
not cover `[output].enabled = false` in region mode; `meds_main` parses the TOML twice;
`region_finalize` on a region that never opened indexes an unallocated `reg%poly`.

---

## 2. Forcing — `meds_met_driver`, `meds_era5land_reader`, `meds_co2_series`, `meds_lapse_rate`

`D` = `src/forcing/meds_met_driver.f90`, `R` = `meds_era5land_reader.f90`, `L` = `meds_lapse_rate.f90`,
`F` = `src/fast_dynamics/driver/meds_fast_dynamics.f90`.

**As built.** `met_open` reads the CO₂ series, then per backend: the netCDF file is validated and only
the run's record range read into `src%series` (D:216-229); the archive resolves its cell(s) from the
static file, lays the months of the recycle window or run period end to end as one hourly axis, checks
every month file exists (D:725-838) and loads axis month 1. `met_cursor_init` (D:267) binds one polygon
and loads its first bracket through `read_record` (D:954), the one ingest pipeline for both backends
(raw values per backend, NaN assertions, terrain lapse, shortwave partition). Per slow step
`met_prefetch` (D:292, archive only) loads the step's month plus the record before it into a
**one-month buffer** (`month_loaded`, D:881-885); per sub-step `met_advance` (D:404) slides the bracket
on the recycle-mapped axis and `met_instant` (D:488) interpolates, splits phases, reconstructs SW,
synthesizes LW, looks CO₂ up on **model** time. `fast_dynamics` hoists the `nsub` samples out of the
patch loop (F:487-494); per patch `met_to_cas_top` (L:148) moves wind and temperature to the
canopy-air top with `canopy_roughness` shared with the aerodynamics. During stepping `src` is
`intent(in)`; every mutable per-polygon datum lives in `met_cursor_t`.

Verified clean: lapse units and bounds; RH clipping; liquid-curve humidity everywhere; month seam and
01:00-anchor year wrap; leap handling; CO₂ on model time under recycled met; archive NaN check before
use; every netCDF id closed on every path across months; the doc's §8 worked example reproduces from
L:113-143. The corrections, partition, rain/snow split, lapse and CO₂ are applied in one place for
both backends (D:996-1050, 503-508, 546-600) — the design goal is met there.

**F1 `[BLOCKER]` Archive + recycling crashes at the year wrap for any anchor not at 01:00 — S.**
`met_prefetch` detects the wrap only when the **step start** already lies at or past the window's
last record (D:300-303) and then carries that record (D:310-316). For any other anchor the seam
interval `[irec_cycle_last, +dt)` falls inside a step: `met_advance` takes the seam branch
(D:434-439) → `load_wrap_bracket` (D:710-717) → `read_record(irec_cycle_first)` → `locate_record`
finds that record's month neither loaded (the buffer holds one month) nor carried →
`error stop 'a forcing record was not prefetched … (internal error)'` (D:912-918). Traced for
`recycle_start = 2021-01-01 06:00` (fires at 05:07 on the first model day) and for
`2021-01-01 00:00` — the most natural anchor — (fires at 23:07 on the last day of the first cycle
year). Every config in the tree and in `~/claude_workspace/meds_runs` uses `01:00`
(`examples/example_biophysics/meds_config_july.toml:201`, `test_met_era5land.f90:286`), which is why
it never fired; `validate_config` does not restrict the hour (`meds_config.f90:769-787`) and
`forcing.md` §9 promises "a cycle may begin anywhere in the calendar".
Fix: the window's first record always lies in axis month 1 (D:782, 361-365), so stash it in a second
carry slot at open (`src%seam_first(cell, var)`), let `locate_record` return a sentinel for
`irec == irec_cycle_first`, and `archive_value` read it; extend the synthetic archive by one month and
test 00:00 and 06:00 anchors across the wrap. Minimal alternative: reject non-01:00 anchors for
`format = "era5land"` in `validate_config` and say so in `configuration.md`.

**F2 `[BEFORE-MERGE]` `met_instant` duplicates the LW-synthesis + ρ block inside the CONST early return — S.**
D:510-529 is a mis-indented copy of D:584-600 (inherited from v0.2.2). On CONST it synthesizes LW
from the type defaults with `kt_last_day` never refreshed (`met_advance` returns at D:411). Wrap only
the file-dependent interpolation (D:531-582) in `if (backend /= CONST)` and let CO₂, LW and ρ run
once at the single exit.

**F3 `[BEFORE-MERGE]` Science/config docs describe the pre-split API — S.**
- `forcing.md:114-116` — `met_advance(drv, now)`, `met_instant(drv, now)`, "`met_open` … loads records
  #1–#2": the signatures are `(src, cur, now)` and the bracket is loaded by `met_cursor_init` (D:281).
- `forcing.md:600` lists `ensure_month` (does not exist; `met_prefetch`/`load_axis_month`) and omits
  `met_cursor_init`/`met_prefetch`.
- `forcing.md` §9 "no record is skipped or double-counted": `recycle_model_to_file` is one-directional
  (D:667), so a leap **file** year under a non-leap model year never reads its Feb-29 records.
- `configuration.md:224-225` "the rest of `[site]` (reference heights, profile and lapse switches)
  still applies" — heights and profile now live in `[forcing]`, and the old keys are rejected
  (`meds_config_io.f90:640-646`).
- `src/forcing/README.md:7-8` "links the shared foundation and the netCDF C bindings only" vs
  `CMakeLists.txt:273` linking `meds_config`; `meds_forcing_config.f90:6` "Placed in src/shared".

**F4 `[BEFORE-MERGE]` Comments narrating history (rule 1), one now wrong — S.** `file_lookup_sec`
header (D:623-628) describes three regimes including a "LEGACY absolute-seconds span-wrap" while the
body has two (D:633-637); `met_advance` header "EOF (recycle by whole file spans)" (D:402); D:39,
186-189, 325-330, 416-418, 650-653; `meds_forcing_config.f90:48-52` (deleted `SWPART_SIB`), `:126-139`
("the old behaviour"). The history is in the CHANGELOG.

**F5 `[BEFORE-MERGE]` One rejection path in `met_open` leaks `src%ncid` — S.** The
`validate_file_against_config` failure with `stat` present closes the raw id but leaves
`src%ncid ≥ 0` (D:190-196); the other three paths call `met_close(src)` (D:122, 137, 207). A later
`met_close` then `nc_close`s a closed id → `nc_check` error-stops. Use `met_close(src)` there too.

**F6 `[AFTER-MERGE]` `meds_met_driver.f90` carries three sources in one module — M.** 1256 lines (770 in
v0.2.2), 6 public entries, ten `backend ==` sites (D:127, 132, 276, 277, 298, 411, 510, 610, 962,
1035). The netCDF-file source (D:146-229, 685-695, 934-951, 1054-1132, 1147-1155, 1173-1254) and the
archive source with its axis (D:292-318, 725-930) are separable from the backend-independent
cursor/recycle/ingest core. Each source exposes `open`, `raw_record(src, cur, irec, raw)`, `close`;
`read_record` keeps the shared pipeline over `raw`, collapsing its two branch points to one call.

**F7 `[AFTER-MERGE]` The hourly cadence is a literal in five places — S.** `3600` at D:806, 849, 851,
R:380 and `meds_config.f90:773`; `era5land_month_hours` (R:104-107) encodes it a sixth way. One
`ERA_DT_SEC` in the reader. Same family: two hand-rolled `seconds since` parsers (D:1147-1155,
R:373-374) → one function in `meds_time`.

**F8 `[AFTER-MERGE]` `read_scalar`/`read_scalar_default` and the MEDS-file name list — S.** D:1091-1104
and D:1118-1132 differ by an optional `default`; the 13 field names appear as `FIELDS` (D:1058-1060)
and again as literals in `read_record` (D:975-1029), with `series_field` string-matching per field per
record (D:1082-1089). Resolve field indices once at open.

**F9 `[AFTER-MERGE]` Dead code and a dead config field — S.** `wind_log_profile` (L:38-48) has no
production caller since #305; `forcing_config_t%rad_sw_ground_const` (`meds_forcing_config.f90:125`)
is parsed nowhere and read nowhere.

**F10 `[AFTER-MERGE]` `met_instant` could be `pure` — S.** Every callee is pure and both state
arguments are `intent(in)` (D:488-602); `pure` makes the "safe to sample ahead of the threaded patch
loop" claim compiler-enforced.

**F11 `[AFTER-MERGE]` Archive fill values are trusted to be NaN — S.** The reader's only missing-value
test is `ieee_is_nan` (R:342-347); it never reads `_FillValue`. Let `check_month_file` assert
`_FillValue` is absent or NaN.

**F12 `[AFTER-MERGE]` A run can never start in the archive's first month — S.** `start_time` must be
00:00 (`meds_config.f90:787`) and `archive_month_of` steps back 1 s (D:851), so the axis always
includes the previous month and `met_open` demands `…_<YYYY−1>12.nc` (D:810-825). Document it, or
let the axis start at the run's month under `start_clamp = "hold"`.

**Config keys.** Parsed ∩ `meds_config_main.toml:403-489` ∩ `configuration.md` agree on the
user-facing set. Parsed but in neither the TOML nor `configuration.md`: `forcing.lw_clear_form`,
`forcing.lw_cloud_a` (only in `forcing.md` §11). In the TOML but not `configuration.md` (which
defers): `timestep`, `avg_convention`, `sw_partition`, `lwdown_source`, `start_clamp`,
`file_template`, `static_file`, the `[site]` location keys.

Nits: `assert_finite` reports `src%grid_index` (always 1 for the archive) rather than the cursor's
cell (D:996-1002); `meds_forcing_config.f90:99` documents `backend` as `"netcdf" | "const"` and `:33`
still says `"legacy_file"`; `meds_forcing_types.f90:10-11` claims to link `meds_shared` only;
`forcing.md` §8 quotes lapse rates in K/km beside a key in K/m; `MAX_RECYCLE_YEARS` lives in the
config root and pulls `use meds_config` into the driver (D:39). The CO₂ `const` branch (D:504-508,
242) is simpler than a synthetic one-row series — leave it.

---

## 3. Output layer and the slow-dynamics changes

**As built.** `output_files_t` (registry, `diag_params_t`, per-tier `stream_file_t`, polygon
coordinates) is shared read-only during a step; `output_buffers_t` (per-(var, tier) `integ_buffer_t`,
`t_open`, scratch record, closed-record queue, fast staging) is one per polygon. Only the serializer
writes netCDF ids. The tick (`output_integrate`, `meds_output_integrate.f90:794`) folds the step
`[prev, now)` into the window holding `prev` (#294), then closes the tiers `now` has left.
Restructuring moved out of the step (#297): `polygon_step` ticks, zeroes the diag blocks, sets
`restructure_pending`; the next step runs `advance_boundary → restructure_stand` first — #294's
slot-bound deferral and writer workaround are gone, so one code path serves all tiers. The `forcing`
group rides the ordinary registry: coarse tiers read a fixed-size `polygon_diag_block` (`PY_*`,
`FLD_PY_DIAG0 = 4100`), the fast tier reads `fast_sample_t` fields, three per-patch rows ride the
patch block. Region runs serialize the same queues through `output_serialize_region →
region_write_record`.

Verified clean: `meds_io_config.toml` = registry (253/253); `order_of_processes.md` §2/§5/§6 vs
`polygon_step` (162→173→177→196→204→209) and `restructure_stand`; the only cohort/patch-count
operators are called from `restructure_stand` (plus the C API), so no D/M record can mix two stands;
disturbance-row units follow the block's "rate × seconds" contract; the new patch forcing rows are
states, so no per-second/per-year slip; polygon plumbing is `bufs(:)` plus a handle, not an integer
threaded through signatures; first window after a restart `t_open = prev`.

**O1 `[BEFORE-MERGE]` Region I/O on a failed polygon error-stops instead of writing what closed — S.**
`meds_region.f90:228-236` calls `io_phase` on failure "so what closed before the failure is written",
but `output_serialize_region` (`meds_output_manager.f90:70-74`) demands every polygon's `queue(t)%n`
equal polygon 1's: polygon 1 has the month's records so far, the failed polygon its partial count
(the failing step's tick runs before the guards, `meds_polygon.f90:196` vs `:209`), later polygons
zero → `error stop`, and the partial month is lost for every polygon. Write
`nmin = minval(bufs(:)%queue(t)%n)` records on the exit path and keep the equality check for the
normal path.

**O2 `[BEFORE-MERGE]` FAST tier does not stamp by window start — S.** `meds_fast_dynamics.f90:490-498`
stage `fast_time(isub) = step_start + (isub−1+forcing_sample_frac)·dt_fast` (default 0.5,
`meds_config.f90:255`); `meds_output_integrate.f90:980` uses it as `t_open(1)`;
`meds_output_stream.f90:193-203, 414-420` write it as `time/hour/minute/second` with `long_name`
"(period start)". D/M/Y now stamp the true window start (#294) and `diagnostics.md:98` promises
"period start"; the tower comparison the forcing group exists for is offset by `frac·dt_fast`. Stage
`step_start + (isub−1)·dt_fast` (keep `t_sample` for the met lookup), and fix the
`meds_output_types.f90:336` comment "sub-step midpoint stamps" — or, at minimum, fix the label.

**O3 `[BEFORE-MERGE]` Rules and source headers still say the tick runs after fuse-fission — S.**
`.claude/rules/output.md:48`, `.claude/rules/state-demography.md:76-77`,
`meds_site_diag_types.f90:24-27, 119-125`, `meds_output_integrate.f90:134-136` state the pre-#297
order. The invariant "never read the tendency bundle" still holds (the daily `sort_cohorts` reorders;
`deriv` is not reordered) but its stated reason is now false, and rules files steer future edits.
`meds_output_integrate.f90:500-503` is an orphaned comment describing a `site%deriv` read that no
longer exists — it documents the forbidden pattern.

**O4 `[BEFORE-MERGE]` `diagnostics.md` inventory stale after #305 — S.** `:280, 286-288`: "248
registered", energy 47, patch 24, forcing 25; the registry has 253 / 49 / 29 / 28 (the five `*_patch`
rows). §5 (`:212`) lists the forcing site scalars and `*_fast` twins but not the three patch rows.
Regenerate from `--dump-io-config`.

**O5 `[BEFORE-MERGE]` Gap-inheritance rationale is stale for the driver path; the test seeds a case production no longer produces — S.**
`meds_demography_patch_fusefiss.f90:603-607`, `meds_site_diag_types.f90:511-518`,
`test/test_disturbance.f90` (seeds `PD_LE` with `w = dt_slow`). Since #297 the disturbance runs after
`reset_step_diagnostics`, so at `patch_diag_inherit` every donor has `w = 0` and `v` holds only this
boundary's events. The mechanism is still required and correct (Σ area·v conserved; all patches later
share `w = dt_slow`, so `patch_diag_blend` not blending `w` is harmless) but "its ground was donor
ground for every fast sub-step" now describes only the C API path. Reword; seed `PD_MORT_C_CULL` with
`w = 0` in the test. Structurally `patch_diag_inherit` is the N-donor form of `patch_diag_blend`.
The patch block declares one fusion kind for all rows (`diag_types.f90:497-499`) — "declared once" holds.

**O6 `[AFTER-MERGE]` Forcing echo is enumerated twice (three times for three fields) — M.**
`meds_output_types.f90:203-217` (12 `fast_sample_t` fields) + `meds_fast_dynamics.f90:503-508` (12
copies) + `meds_output_integrate.f90:214-226, 950-961` (12 ids, 12 cases) + `meds_output_registry.f90:847-870`
(12 rows), parallel to `PY_*` (`diag_types.f90:223-238`) + `accumulate_polygon_diag`
(`fast_dynamics.f90:1265-1285`) + 13 rows (`registry.f90:812-837`). `air_temp`, `sw_in`, `atm_co2`
additionally travel a third route, area-summed per patch through `red_fast` although site-uniform.
Group membership is split: `sw_in_site` (ENERGY), `precip_site` (WATER), `air_temp_fast`/`sw_in_fast`
(ENERGY), `atm_co2_fast` (CARBON) — `forcing = false` leaves five echoes on. Simpler shape: give
`fast_sample_t` a `forcing(N_PYDIAG)` array filled by the same kernel as the PY block and resolve
`SRC_F_PY0 + PY_*` by index in `extract_fast_scalar`, as `extract_scalar_source:771` already does for
`FLD_PY_DIAG0`; move the five strays to `GRP_FORCING`.

**O7 `[AFTER-MERGE]` Region serializer is a copy of the site serializer — M.** `meds_output_stream.f90`:
time companions 193-196 ≡ 536-539; coordinate definitions 212-238 ≡ 558-583; globals 248-250 ≡
621-623; coordinate writes 254-256 ≡ 632-634; `def_registry_var` 337-385 re-implemented inline
590-616; `write_one_record` 396-456 vs `region_write_one` 646-716 — ~130 duplicated lines,
`region_open_file` 148 lines. Shared `define_axes`/`write_axes`; `def_registry_var` with an optional
leading polygon dim; one `put_record_var` for both writers.

**O8 `[AFTER-MERGE]` Site writer unpacks the handle it is handed — S.** `meds_output_manager.f90:33-36`
passes nine fields of `files` into `stream_write_record(...)` (`stream.f90:36-44`) while
`region_write_record(files, bufs, t, i)` takes the handle.

**O9 `[AFTER-MERGE]` `strict_caps` is a dead knob; `cohort_max` is enforced only at serialization — S (pre-existing).**
`meds_output_config.f90:102` documents it; the only other reference is the parse
(`meds_config_io.f90:753`). `cohort_source_field` (`integrate.f90:449-457`) writes `x(1:site%cohort%n)`
into a `cohort_max` buffer and `fast_dynamics.f90:794-798` writes `fast_coh_*(i, isub)` at global slot
`i` with no cap check; the first guard is `stream.f90:169`, up to a month later. Same class as #247.

**O10 `[AFTER-MERGE]` `fast_interval_steps` is not validated against `n_fast_per_slow` — S.** Nothing checks
`mod(n_fast_per_slow, fast_interval_steps) == 0`; `meds_polygon.f90:280-281` closes the fast tier on
a running `fast_step_total`, so 96/5 makes a fast window straddle two slow steps and, across a month
boundary, the restructured slot set (`stream.f90:54-57` assumes that cannot happen within a day).
Reject in `validate_config`. (Unverified at runtime.)

**O11 `[AFTER-MERGE]` Three copies of the calendar-boundary sequence — S.** `restructure_stand`
(`meds_vegetation_dynamics.f90:302`) is the driver's; `meds_c_api_demography.f90:205-215` re-lists the
same operators inline (no ledger, no canopy-depth refresh); `meds_advance_slow` (`:116`) calls it
after the step with a comment (`:107-108`) claiming "as the driver does it", which is no longer the
driver's order.

**O12 `[AFTER-MERGE]` `apply_patch_disturbance` is 199 lines — M (pre-existing; #295 added to it).**
`meds_demography_patch_fusefiss.f90:417-616`. `patch_diag_grow` at `:608` is redundant
(`patch_ensure_capacity` at `:489` already grows it).

Nits: `meds_output_registry.f90:26-27` two `use meds_output_types` lines, `:645-646` mis-indented row;
`meds_output_config.f90:105-112` says eight groups / six on (nine / seven);
`meds_output_manager.f90:3-4` "the ONLY flush" (two now); `cas_depth_patch` and `cas_depth_site` are
one quantity under two source ids and groups; two "off" representations in the registry (`enabled`
vs `FREQ_NONE`); `MEDS_IO_DESIGN.md:607-608, 652, 676, 1122` still say "flush BEFORE fiss/fuse"
(acceptable for a never-renumbered record; a pointer at `:607` to the `:682-687` amendment would
help); `order_of_processes.md` §5 omits that `advance_boundary` is gated on `slow_on`
(`meds_stepper.f90:124`).

---

## 4. Python tools and examples

**As built.** `download_era5land_gdex.py` / `download_era5land_cds.py` fetch raw ERA5-Land into a
resumable pool with a JSON-lines log; `build_era5land_static.py` writes the static file;
`build_era5land_archive.py` writes one global float32 `(time, lat, lon)` file per variable-month,
de-accumulating `tp/ssrd/strd`, with a manifest and OD2 raw deletion — GDEX and CDS through a
`read_band` adapter; `postprocess_era5land.py` cuts either source into box files;
`make_forcing_file.py` writes the single `(time, grid)` file from the archive or box files;
`make_co2_file.py` builds the shipped CO₂ series. `py_compile` 12/12; SPDX header and `__main__`
guard on every file; no import-time side effects; 0 unused imports, bare `except`, mutable defaults.
`pyflakes` is not installed in any env here.

**P1 `[BEFORE-MERGE]` Archive longitude axis is −179.9 → 180.0, not −180 → 179.9 as documented — S.**
`build_era5land_static.py:69-70` and `build_era5land_archive.py:90` map `lon > 180 → lon − 360`, so
the native 0.0…359.9 axis keeps 180.0 and never produces −180.0. Every description says otherwise
(`build_era5land_archive.py:9`, `build_era5land_static.py:6`, design §14.2). Not a runtime bug — the
reader wraps into `[lon(1), lon(1)+360)` (`meds_era5land_reader.f90:228`) and `make_forcing_file.py:224`
uses modulo — but the archive is already built for 2002–2026, so **fix the docstrings and §14.2, not
the `>`** (that would orphan every existing file against a rebuilt static).

**P2 `[BEFORE-MERGE]` Example README numbers and figures predate the configs' physics change — M.**
Both example configs now set `apply_elevation_lapse = true`, `grid_elevation = 367.5` and the new
`[forcing]` height keys (`meds_config_july.toml:168-185`, `meds_config_spinup.toml:151-168`), which
`CHANGELOG.md:161-190` measures as +0.31 K, +2.2 W m⁻², friction velocity +33 %, GPP +0.6 %. Yet
`examples/example_biophysics/README.md:193-195` still quotes the stage-1 end state from `main`
(14 cohorts / 2 patches, LAI 4.16, AGB 9.6, 16.2 kgC m⁻²), none of the four PNGs changed, and the
"Notes on the configuration" section (`:284-330`) never mentions the lapse or the heights. Rerun
`run_example.py`, refresh numbers and figures, note that 367.5 m is the orography
`make_forcing_file.py` prints (`:238-239`).

**P3 `[BEFORE-MERGE]` De-accumulation, unit factors, clip flags and the variable/units tables exist twice — M.**
The rule: `build_era5land_archive.py:288-296` (vectorised) vs `make_forcing_file.py:95-108` (per-cell
loop). Factors: `RHO_W/SEC_PER_HOUR` at `make_forcing_file.py:196-198` vs
`ARCHIVE_VARIABLES[..]["factor"]` (`era5land_common.py:62-69`). Clip flags hard-coded at
`make_forcing_file.py:190-192` vs `clip_negative`. Variable lists at `:251, :353` and the `meta`
table `:382-390`, whose long names differ from the archive's. One `common.deaccumulate(...)`;
`make_forcing_file.py` takes factor, clip, units and long_name from `ARCHIVE_VARIABLES`.

**P4 `[BEFORE-MERGE]` `make_forcing_file.py` duplicates the Fortran's wind floor and drops the wind vector — S.**
`Wind = max(√(u²+v²), 0.1)` (`:61, 195, 275`); the reader floors every backend itself at `U_MIN = 0.1`
(`meds_met_driver.f90:66, 551`), and `forcing.md:66` documents an unfloored `Wind`. The single-file
format takes optional `u10`/`v10` (`forcing.md:34, 58-61`; D:182-183, 978-981) and the archive carries
the vector, but the file cut from the archive writes only `Wind`, so `has_wind_vector` is false. Write
`u10`, `v10` (keep `Wind`), delete `U_MIN`. The Fortran is canonical.

**P5 `[BEFORE-MERGE]` Stale statements inside the same release — S.** `CHANGELOG.md:115`
"`prep_era5land_forcing.py` stays until the forcing reader upgrade lands" and `:307-308` "now reads
their box files" are under `[Unreleased]` while `:302` removes it. `download_era5land_cds.py:4`
"Adapted from scripts/download_era5land.py" cites the removed file (rule 1).

**P6 `[AFTER-MERGE]` Source-specific naming and grid decoding are copied rather than shared — M.** The
CDS raw-file pattern is a literal in three scripts (`download_era5land_cds.py:181`,
`build_era5land_archive.py:199`, `postprocess_era5land.py:262-263`); GRIB grid keys decoded twice;
the 0…360 → ±180 reorder twice; `EPOCH` twice; "stamp − 1 h → month" twice; `postprocess:156-159`
re-derives `common.month_interval`; `make_forcing_file.py:248-249` re-derives `common.interval_stamps`;
`build:400-404` re-implements `common.parse_variables`; `PROCESSING_VERSION` in two files.

**P7 `[AFTER-MERGE]` Downloader asymmetry — S–M.** GDEX retries transient failures
(`era5land_common.py:268-324`); CDS has none and an exception inside `fetch` aborts the whole plan
(`download_era5land_cds.py:206-209`). GDEX keeps every netCDF/HDF5 call on the main thread
(`download_era5land_gdex.py:79-80`); CDS runs `raw_complete` inside worker threads (`:182, 196`).

**P8 `[AFTER-MERGE]` Long multi-responsibility functions — M.** `build_one` 118 lines
(`build_era5land_archive.py:237-354`), `main` 78 (`:381-458`), `make_forcing_file.main` 85
(`:322-406`), `build_era5land_static.main` 90 (`:39-128`).

**P9 `[AFTER-MERGE]` Silent site defaults — S.** `make_forcing_file.py:332-341` defaults `--lat/--lon`
to Ithaca and `--elevation` to 320 m; omitting `--lat/--lon` silently writes an Ithaca file. FD12
forbids a default period; a default location is the same trap. `elevation(grid)` means the cell's
orography from the archive but the user's site elevation from box files.

**P10 `[AFTER-MERGE]` Documentation gaps — S.** `build_era5land_static.py` needs a GDEX pool, so a
CDS-only archive still needs one GDEX download (`src/forcing/README.md:103-105` implies either
suffices); `environment.yml:2` says "download scripts" but serves the builders too; `forcing.md:30-42`
CDL omits `_FillValue = 1e20` and `wind_meas_height_m`; `make_forcing_file.py:70` cites
`meds_thermo%sat_vapor_pressure` (the module is `meds_therm_lib`). `--help` texts, CHANGELOG flags,
`meds_config_main.toml:392-401`, `run_example.py:178-182` and the example README commands all match
the argparse definitions.

Nits: `common.VARIABLES` fourth field unused; `json.load(open(path))` without closing
(`era5land_common.py:248, 259`); `gdex_reader` checks raw `longitude` but never `latitude`
(`build_era5land_archive.py:105-115`); `build_era5land_static.py:60` takes "the first non-coordinate
variable"; `fcntl` makes the tools POSIX-only. The CO₂ file header matches `FORMAT_BLOCK` and
`forcing.md` §12 exactly; 1023 rows, 278.007 at the join.

---

## 5. Documentation consistency

**Headline:** the docs `beta` wrote (the regional section of `configuration.md`, `forcing.md` §8/§12,
`order_of_processes.md`, the forcing README, the example README, the dev_plans tombstones) are
substantially right — ~100 constants and formulas in `forcing.md` match the code, all 55 `§` cites
from `src/` resolve, `order_of_processes.md` matches the call order literally, `meds_io_config.toml`
matches the registry exactly. The drift is in the docs `beta` did **not** touch (`README.md`,
`ed2_comparison.md`, `src/README.md`, the rules files, `CLAUDE.md` counts) plus a few
self-contradictions inside `beta`'s own CHANGELOG and `configuration.md`.

### 5.1 Config surface (parser vs `configuration.md` vs shipped TOML vs examples)

All `[run]/[forcing]/[site]/[region]/[output]` keys appear in all three, with matching defaults, except:

| key | parser | `configuration.md` | `meds_config_main.toml` | note |
|---|---|---|---|---|
| `output.fast_interval_steps` | `meds_config_io.f90:754` | `diagnostics.md:192` spells it `[output.fast].interval_steps` | **`:322` `[output.fast] interval_steps = 4` — never read**; the default 4 masks it | `[A]` pre-existing; fix toml + doc |
| `forcing.lw_clear_form`, `lw_cloud_a` | `:585-592` optional | absent | absent | `[A]` only `forcing.md:459-467` |
| `run.slow_on` | `:994` | absent | absent | `[A]` |
| `forcing.format = "const"` | `:470` | lists era5land/netcdf only | comment `:422` | `[A]` |
| `fast.integration_scheme`, `state.write_output`, `state.output_interval_years`, `state.cohort_max`, `state.patch_max` | **not read** | — | in `example_demography/example_config_main.toml:32, 74, 77-79` | `[A]` dead keys, pre-existing |

- `[B]` `configuration.md:224` "the rest of `[site]` (reference heights, profile and lapse switches)
  still applies" — `meds_config_io.f90:640-643` rejects `[site].reference_height`, `wind_meas_height`,
  `apply_wind_profile`, `wind_roughness_z0` (#305); the same doc says so at `:178`.
- `[B]` `configuration.md:121` "Around 208 variables across 8 groups", `meds_config_main.toml:270`
  "~210 … 8 groups", `README.md:56` "~208 variables" → 253 variables, 9 groups.
- `[A]` `configuration.md:137` "`[io]` … still loads in v0.2.x" — the shim is still in the code
  (`meds_config_io.f90:1151-1154`) and the release being cut is 0.3.0; see 5.5.

### 5.2 Output variables

Registry 253 rows / 9 groups / 7 axes = `meds_io_config.toml` exactly (name, units, default streams,
group) — clean. `diagnostics.md` names 52 variables, all registered with consistent units and
aggregation; it delegates the inventory to `--dump-io-config` by design. Defect: `[B]`
`diagnostics.md:280-288` "248 registered", energy 47, forcing 25, patch 24 → 253 / 49 / 28 / 29
(the five #305 patch rows); `CHANGELOG.md:74` "recounted to 248" is stale for the same reason.

### 5.3 Science docs vs code

`forcing.md`:
- `[B]` `:600` code map lists `ensure_month` (does not exist); the archive path is `met_prefetch`
  (`meds_met_driver.f90:292`), `load_axis_month` (`:888`), `locate_record` (`:904`);
  `met_cursor_init` (`:267`) is missing from the row.
- `[B]` `:67-68` "every value is read as a `[1,1]` hyperslab … never loads a whole series" —
  `read_series` (`:212-229, 1054-1079`) reads the whole window at `met_open` (R1); the doc's own §2
  (`:124-126`) says so.
- `[B]` The `format = "era5land"` restriction (`dt_slow = "1d"`, `start_time` 00:00,
  `meds_config.f90:783-787`) is stated in `configuration.md` and the TOML comment but nowhere in
  `forcing.md` or `src/forcing/README.md`.
- `[A]` `:41` "global attributes never read" (the reader checks `avg_convention`/`sw_input_kind`,
  `:1206-1230`); `:115` pre-split API (see §2 F3); `:76` omits the 0.1 m/s wind floor; `:13` "links
  the config leaf and netcdf_c and nothing else" (links `meds_config`, `CMakeLists.txt:273`); `:287`
  cites `snow_biophysics.md §1` (no numbered sections); `:94-95/:608` recipe omits
  `build_era5land_static.py`; equation (11) is used twice (`:307, :456`), (16) precedes (14)/(15).

`order_of_processes.md`: the per-step order, the calendar-boundary order (tick → reset → flags →
checkpoint → `advance_boundary` at the next step's start), the fast sub-step sequence, the
slow-dynamics list, the region month loop and `restructure_pending` none|month|year all match
literally; 25/26 code-map routines and 15/15 modules exist.
- `[B]` `:20, :220` `driver_init` — the routine is `driver_open` (`meds_driver.f90:54, 92`).
- `[A]` `:81-84` phenology air-T sum listed as reset before the sub-steps (zeroed in the fold after,
  `meds_fast_dynamics.f90:875-891`); `:96-98` per-patch setup order reversed (freeze `:604-608`
  precedes z_top/roughness `:635-636`); `:208-210` "I/O phase writes the checkpoint" — site runs only
  (`meds_region.f90:250-259` writes none); `:56/:169` the boundary gate is `slow_on`
  (`meds_stepper.f90:124`), not `demography_on`; `:6` "default 900 s" — `dt_fast` is required
  (`meds_config_io.f90:1005`); `:118-122` omits the GPP/respiration/ψ handover (`:830-838`); `:178`
  omits `sort_patches` (`meds_vegetation_dynamics.f90:359`).

### 5.4 Navigation / meta docs

1. `[B]` **Test count.** `CLAUDE.md:46` "45 tests" and "6 s"; `docs/ed2_comparison.md:81, 365` "49
   CTest targets". `ctest -N` on `beta` = **53** (§6), ~22 s. `docs/building.md` gives no number (clean).
2. `[B]` **`docs/ed2_comparison.md`** (beta changed only `:319`): `:63` "no gridded/regional runs.
   MEDS v0.2 is a site model", `:101` "supports one site", `:367` "Scale | One site", `:378` — region
   runs shipped; **`:321` "CO₂ … or from the forcing file; echoed into the output"** — a met file with
   `CO2air` is now *rejected* (`meds_met_driver.f90:1232-1243`); CO₂ is `co2_const` or `co2_file`.
   Title and `:7` pin "v0.2.0".
3. `[B]` **`README.md`**: `:61` "Status: v0.2.0", `:121` license "v0.1.0 through v0.2.1", `:56`
   "~208 variables"; the feature table has no row for regional runs / the ED_ERA5land archive /
   prescribed CO₂, and the forcing row omits the canopy-air-top move.
4. `[B]` **`.claude/rules/output.md:48`, `.claude/rules/state-demography.md:76`** "The output tick
   runs after the monthly fuse-fission" — inverted by #297 (§3 O3). `state-demography.md:69` "Order of
   operations" still gives the in-step order. `output.md:58` "`[io]` is the restart stream only" — the
   block is `[state]` (pre-existing).
5. `[A]` **`src/README.md`** `:16` "30.3 k lines · 86 modules · 19 CMake libraries" → 35.5 k / 92 /
   21 (new `meds_testarchive`); `:63` forcing/ "1 067" → 2 747 lines incl. `meds_era5land_reader`,
   `meds_co2_series`, `meds_lapse_rate`; `:70` main/ "meds_stepper, meds_driver, meds_main, 740" →
   1 372 incl. `meds_polygon`, `meds_region`; io/ 4 563 → 5 372; config/ 2 920 → 3 702
   (`meds_region_opts`). None of the six new modules is named. The library graph is still correct.
   (Tagged `[A]` by the audit; I would do it before the release — it is the map new readers open first.)
6. `[A]` **`src/forcing/README.md:72-73`** "`meds_main` opens the reader … refreshes a local context
   overlay per sub-step" — `met_open` is called from `meds_driver.f90:176` / `meds_region.f90:104`,
   and the overlay (`apply_met_to_ctx`) is gone (#292). `src/fast_dynamics/README.md`: clean.
7. `[B]` **`docs/dev_plans/README.md`**: the tables match the files; "Thirty-four" = 34; the
   `MEDS_FORCING_DESIGN.md` tombstone lists what shipped and defers #302/#303 — its
   last-open-item claim holds; the polygon-plan status matches the plan header. Gap:
   **`MEDS_V02_RELEASE_PLAN.md`** sits in the live directory and is in no table (pre-existing); by the
   directory's own rule it needs a tombstone and a move to `archive/` — a release is the moment.
8. Clean: `docs/README.md`, CLAUDE.md paths, `building.md` traps = CLAUDE.md traps,
   `examples/example_biophysics/README.md` recipe (every CLI flag exists; archive link updated).

### 5.5 CHANGELOG / ROADMAP / release readiness

1. `[B]` **Cutting 0.3.0:** no `## [0.3.0] — date` heading (`:15` is `[Unreleased]`); `:1858` link is
   `v0.2.2...beta`. Version strings to bump: `CMakeLists.txt:17` `VERSION 0.2.2`,
   `python/pyproject.toml:16`, `python/meds/__init__.py:19`, `README.md:61, 121`. No `CITATION` file;
   `AUTHORS`/`NOTICE` present.
2. `[B]` **Contradictions inside `[Unreleased]`:** `:115` "`scripts/prep_era5land_forcing.py` stays
   until the forcing reader upgrade lands" vs Removed `:302`; `:307` says the removed script "now reads
   their box files"; `:153` "Box selection … the config does not expose it yet" vs the Region-runs
   entry (`[region].box_nwse`); `:74` "recounted to 248" (now 253).
3. `[B]` **Misfiled and uncited entries.** The `### Fixed` entries for #294/#296 (`:396`) and #295
   (`:421`) sit under `## [0.2.1]` (`:353`) — `[Unreleased]` has no `Fixed` section at all (§7 B3).
   The PR-citation rule (`:11`) is not met by the entries for #281–#289, #291–#293 (no PR number).
   #304/#307 are docs-only, acceptable.
4. Every routine, attribute and option `[Unreleased]` names exists in the code (`forcing_qair`,
   `restructure_pending`, `has_wind_vector`, `met_cursor_init`, `output_files_t`/`output_buffers_t`,
   `reset_step_diagnostics`, `advance_boundary`, `restructure_stand`, `polygon_step`,
   `fill_forcing`/`fill_aenv`, `MET_PATH_LEN = 1024`, every script flag); numbers-moved entries give
   before/after. Clean.
5. `[B]` **`docs/ROADMAP.md`:** `:86` (#104), `:109` (#167), `:239` (#190/#146) "Deferred to v0.3.0"
   and `:158` "Open, v0.3.0 question" — none shipped, so the labels become false at release; retarget.
   `:143` "Remove the `[io]` deprecation shim — scheduled, post-v0.2.x" is now due and is the one
   deferred item **without an issue number** (also `configuration.md:137`, `output.md:58`): drop the
   shim in 0.3.0 or file an issue. Other deferred items with no issue: `:158` storage-maintenance
   default, `:165` per-layer root biomass, `:169` hydraulic-trait calibration.
6. ROADMAP §8 vs beta: #183 rewritten for R2 (R3–R6 listed), #184 removed, #182 → #257, #302/#303
   added with numbers. Clean.

### 5.6 Retired names (live references only)

| name | live references | verdict |
|---|---|---|
| `apply_met_to_ctx` | `MEDS_SNOW_DESIGN.md:315` (also cites `meds_fast_loop.f90:482`, a file that no longer exists); `MEDS_PRODUCTION_INTEGRATOR_PLAN.md:561, 2077, 2098` | `[A]` live plans; add a note |
| `prep_era5land_forcing.py`, `scripts/download_era5land.py` | CHANGELOG (5.5.2) and the provenance line `download_era5land_cds.py:4` | fix those two |
| `output_manager_t` | `MEDS_IO_DESIGN.md:164-201` (carries the R2 note at `:168`), `MEDS_IO_V01_PLAN.md:112, 416` | acceptable (reference records) |
| `met_driver_t` | polygon plan, describing the split | clean |
| bare `dev_plans/MEDS_FORCING_DESIGN.md` links | none remain | clean |
| `[site].reference_height` & co. | only the rejection code and docs describing the removal; `make_forcing_file.py:401` still writes an unread `wind_meas_height_m` attribute | harmless |

---

## 6. Build and test on a clean export of `beta`

Export of `e33cb56` in the job's tmp dir; the 1-CPU login node; ifx 2026.1.0, gfortran 11.5.0 and
gcc/15.2.0, CMake 3.31.8, netCDF from the `meds` env. **nvfortran is not installed here**
(`CLAUDE.local.md`), so the NVHPC back end CLAUDE.md asks for on new modules was not exercised.

| Build | Result | Tests |
|---|---|---|
| ifx Release (`-O2`) | OK, 132 s | **53/53** (21.7 s) |
| ifx Debug (`-stand f18 -check all -fpe0`) | OK, 109 s | **50/53**: `soil_column_config`, `soil_biogeochem` abort (FPE); `column_ark` FAIL — **all three identical on `main` v0.2.2** (rebuilt and rerun there) |
| gfortran 11.5 Release | **build fails**: line > 132 is a hard error | — |
| gfortran 15.2.0 Release | identical failure | — |
| gfortran 11.5 Release + `-ffree-line-length-none` (diagnostic) | OK, 129 s | **52/53**: `region` SEGFAULT |
| gfortran 11.5 Debug + line length (4 tests) | OK | `region` SEGFAULT; the same two FPE aborts; `column_ark` passes |
| `libmeds.so` (`-DMEDS_BUILD_PYLIB=ON`, ifx) | OK, 70 s | `python -m meds.plant` self-tests pass; 14/14 Python tests via a `pytest.skip` stand-in (`pytest` is not installed) |
| Region smoke (ifx Release, your 4-cell July-2024 config, outputs redirected) | exit 0, **5.05 s** | `era5-D-202407.nc`: `time = 31`, `polygon = 4`, `soil = 20`; "area conserved, no NaNs" |

Test count: `ctest -N` = **53** (v0.2.2's 51 + `region` + `met_era5land`; `meds_test_era5land_archive`
is a library). `CLAUDE.md:46` "45 tests" and "6 s" are stale — the suite is ~22 s (`region` 11.3 s,
`met_era5land` 3.0 s).

**6.1 `[BEFORE-MERGE]` gfortran cannot build: six lines over column 132 — S (pre-existing).**
The same six lines `CLAUDE.local.md` recorded at v0.2.2 (`meds_fast_be_stage.f90:430`,
`meds_fast_frozen.f90:426`, `test_column_ark.f90:519, 549`, `test_column_dynamics.f90:525`,
`test_column_rk45.f90:738`; 139–151 characters): ifx warning `#5268`, gfortran hard error
(`-Werror=line-truncation` is the free-form default and the project's `-std=f2018 -Wall` does not
relax it), 11.5 and 15.2 alike. `main` does not build with gfortran either. No `beta`-new file has a
line over 132. Wrap the six lines; consider adding the `awk 'length > 132'` check to CI.

**6.2 `[BEFORE-MERGE]` gfortran `region` segfault: a strided component section passed `intent(inout)` — S/M (new code, mechanism verified).**
Backtrace (gfortran Debug, `-fcheck=all` reports no bounds violation): `reset_buffer`
(`meds_output_integrate.f90:282`, `buf%slab(:) = buf%seed` inside `if (allocated(buf%slab))`) ←
`close_tier` (`:857`) ← `output_integrate` (`:826`, the monthly close) ← `tick_output`
(`meds_polygon.f90:286`) ← `polygon_step` (`:196`) ← `region_step_month` (`meds_region.f90:230`),
at the first step of month 2, right after January's I/O phase. That I/O phase is
`call output_serialize_region(reg%out_files, reg%poly(:)%out_bufs)` (`meds_region.f90:254`), for
which gfortran prints `Fortran runtime warning: An array temporary was created`: the actual is a
strided section of a derived-type component whose type (`output_buffers_t`) holds nested allocatable
components (`buf(:,:)` of `integ_buffer_t`, each with `slab(:)` …), and the dummy is
`intent(inout) :: bufs(:)` (`meds_output_manager.f90:64-67`). gfortran's copy-out of that temporary
leaves the originals' nested allocatables pointing at freed memory; ifx passes a descriptor and is
unaffected.
*Experiment (export copy only, restored afterwards):* with the dummy `intent(in)` and the queue
reset moved to the caller, the month-boundary crash **disappears** and every step assertion passes
("every polygon conserves area"); the process then aborts in `region_free` right after the other
`intent(inout)` site, `output_region_close(reg%out_files, reg%poly(:)%out_bufs)`
(`meds_region.f90:310`) — same mechanism, second site. A gfortran defect, but the construct is the
same family as CLAUDE.md's nvfortran trap ("never pass an array-valued function result straight
into a call"). Fix shapes, in order of preference: (a) the region owns a contiguous
`type(output_buffers_t), allocatable :: bufs(:)` and each polygon steps with `bufs(p)` (this is also
the shape R3's threaded loop wants); (b) the two manager routines take `intent(in)` arrays and the
per-polygon mutations (`queue%n = 0`, the final `close_tier`s) run through scalar calls from the
caller. Add "never pass a derived-type component section `a(:)%c` whose type has allocatable
components as an actual argument" to CLAUDE.md's trap list, and run the gfortran suite before the
release: it is the only second compiler available on this machine.

**6.3 `[AFTER-MERGE]` ifx Debug: three pre-existing failures — file issues.**
- `soil_column_config`: `forrtl: error (65): floating invalid` in `pow` ← `soil_theta_from_psi`
  (`meds_hydr_lib.f90:332`, the Campbell branch `(psi/par_a)**(-1/par_n)`) ← `build_soil_hydr_params`
  (`meds_column_params.f90:172`) ← `test_soil_column_config.f90:89/68`. The test's step 5 sets
  `retention = SOIL_RETENTION_CAMPBELL` (`:67`) while `curve_par_a/n` keep the van Genuchten defaults
  (`alpha > 0`), so a negative base under a fractional exponent is NaN. Release silently stores NaN
  `theta_fc`/`theta_wp` — Campbell retention with van Genuchten parameters is accepted by validation;
  that config gap deserves its own issue.
- `soil_biogeochem`: `floating invalid` in `soil_carbon_bad_pool` (`meds_soil_biogeochem.f90:598`) ←
  `test_pool_plausibility` (`test_soil_biogeochem.f90:670`): the test injects `ieee_quiet_nan` to
  exercise the NaN guard, but the guard's ordered comparisons raise *invalid* under `-fpe0`. Test
  `ieee_is_nan` first.
- `column_ark`: `FAIL : ARK PROG-WOOD: cap->0 recovers …` (`test_column_ark.f90:241-242`,
  `abs(tw_tiny - tw_diag) < 0.5`): `tw_tiny` and `tw_diag` are declared (`:208`) and **never
  assigned**; the assertion is vacuous and passes elsewhere by stack luck (gfortran
  `-Wuninitialized` says so).
- Consequence: "the Debug suite is green" is not true on either branch; CLAUDE.md's Debug recipe
  should say so until the three are fixed.

**6.4 Warnings in `beta`-new files.** ifx: none (the six `#5268` are the pre-existing lines).
gfortran: `meds_region.f90:239` `'new_year' may be used uninitialized` — assigned inside the
`do p`/`do s` loops and read after them, undefined when `nstep = 0` or there are no polygons
(`[AFTER-MERGE]`, S); `test_region.f90:40` constant-folded integer division (benign). Pre-existing
elsewhere: 21 `-Wcharacter-truncation` (20 in `meds_config_io.f90`, 64/1024 and 256/1024),
~25 unused variables, `-Wmaybe-uninitialized` at `meds_fast_dynamics.f90:314-316`,
`meds_met_driver.f90:928` (`h`), `meds_init.f90:144`.

**6.5 Smoke-run observations.** The summary reports `energy/water fails = 17856 0` for every
polygon = 6 patches × 2976 fast steps, i.e. the energy budget fails on every patch-step of this
bare-ground run — #290 (§7 B1) under `bottom_bc = "dirichlet"`, identical in your own
`r2region/smoke/run.log`. No `-F-` files despite `[output.fast] enabled = true`: by design —
`manager_restrict_region` strips `FREQ_FAST` and the cohort/patch axes (`meds_output_registry.f90:1131-1137`)
and `configuration.md` says the sub-daily tier is written only for `detail_polygons`. Your earlier
binary also wrote `era5-D-202408.nc`; the current `beta` does not — the #294 fix (the Jul 31 → Aug 1
step is now dated Jul 31). Expected.

---

## 7. Issue triage

`gh issue list --state open` returns 41 issues. `main` is an ancestor of `beta`, so every fix on
`main` is on `beta`. Issue bodies and comments are cached under the job's `tmp/issues/`.

### A. Close with the release

| Issue | Evidence | Suggested closing comment |
|---|---|---|
| **#184** constant 420 ppm CO₂ | PR #301: `meds_forcing_config.f90:122-124` (`co2_source`, `co2_const`, `co2_file`); `meds_co2_series.f90`; `data/co2/…`; `make_co2_file.py`; `configuration.md:182-183`. | "Shipped in v0.3.0 (#301): `[forcing].co2_source = \"const\" \| \"file\"`, a MEDS CO₂ file interpolated on model time (`src/forcing/meds_co2_series.f90`); the shipped CMIP7 series and `make_co2_file.py` build it. Format: `docs/science/forcing.md` §12." |
| **#197** single-precision experiment | Decided by the owner's 2026-09-13 comment (stay `real64`); ROADMAP §11 still lists it as *Candidate*. | "Decided 2026-09-13: MEDS stays `real64` throughout. Closing as won't-do; ROADMAP §11 records the decision." (edit the ROADMAP line in the same change) |

**#183 (multi-polygon runtime) — partly delivered; do not close as-is.** R2 shipped (PR #289), but
the plan header says R3–R6 open (`MEDS_POLYGON_RUNTIME_PLAN.md:3`), ROADMAP §8 uses #183 as the
number for R3–R6, and the issue body asks for MPI, which the plan dropped (§8: job arrays over
tiles). Either keep it open with a comment ("R2 shipped in v0.3.0 (#289); MPI out of scope; R3–R6
per plan §10") or close it and open one issue per phase, updating ROADMAP §8 in the same change.

Not closable: #290 (leak unchanged, see B1); #302/#303 are deferral tickets ROADMAP §8 cites. #182 and
#185, referenced from the forcing work, were already closed 2026-09-14.

### B. Fix before merge (ranked)

| # | Issue | Why now | Premise verified | Minimal fix | Effort |
|---|---|---|---|---|---|
| 1 | **#290** energy ledger leaks −0.3…−4 W/m² (`high priority`) | Every fast run prints `WARNING: N whole-column budget checks breached` (`meds_polygon.f90:312-313`; per polygon `meds_region.f90:296-301`), so the headline feature ships with a conservation warning at every step. | **Yes, by reading.** `bottom_heat_face` applies `−g_deep·(T_bot − deep_temp)` under Dirichlet (`meds_soil_energy.f90:119-129`, at `:182` implicit and `:256` explicit), but every ledger books the bottom face as `frozen%hydrology%geothermal`, hard-set to 0 (`meds_fast_frozen.f90:346`): `bf%soil_enth_in` (`meds_fast_be_stage.f90:187`), `whole_enth_in/out` (`:200, :203`, no bottom term at all), the explicit path (`meds_fast_time_derivs.f90:366`), RK45's `e_in/e_out` (`meds_fast_rk45.f90:823-828`). No other term books the bottom face. The examples switched to `bottom_bc = "dirichlet"`, `deep_temp = 284.75` in v0.2.1 (#267, `CHANGELOG.md:371`), one day after #244's "0 fails" measurement. Magnitude: `g_deep ≈ 1.2/(3.12 − 1.73) ≈ 0.86 W/m²/K`, July base ≈ 289.5 K → ≈ 4 W/m² out of the column, ≈ 0.3 W/m² annually; negative; water unaffected; forcing-independent — all six rows of the issue's table. Not a beta regression. The soil sub-ledger fails the same way (`meds_fast_ark.f90:560`) but its `n_fail` is never reported. | Ledger only, no physics moves: expose the committed bottom face (`energy_flux_t%bottom_heat` exists, `meds_soil_types.f90:154`, set only on the split path `meds_soil_energy.f90:204`) from the BE stage and the time-derivative path, use it instead of `geothermal` in `bf%soil_enth_in`, add it to `bf%whole_enth_in` and to RK45's `e_in`. Outputs unchanged except `resid_energy_site` (state before/after). Re-run the July case: expect 0 fails. | S code, M with the run |
| 2 | **#306** canopy-air depth 20 m on day 1 | Exposed by beta: #305 uses `can_depth` as the height the forcing is moved to (`meds_fast_dynamics.f90:635`), and every region polygon starts from bare ground (`meds_region.f90:118-119`), so day 1 runs at 20 m against the 5 m floor. Moves day-1 numbers in every fast run → belongs in the v0.3.0 baseline. | **Yes.** Default 20 m (`meds_column_state_types.f90:77`); the only writer is `refresh_canopy_depth` (`meds_slow_dynamics.f90:108`, private), called from the slow step (`:61`) and the boundary (`:91`); `polygon_prepare` and `init_fast_reservoirs` never set it; the checkpoint carries CAS enthalpy/shv/CO₂/temperature only (`meds_io.f90:179-182`). The issue's "#298 contributor" claim is doubtful: a boundary checkpoint carries `restructure_pending`, and the resumed run's first step runs `advance_boundary → refresh_canopy_depth` before any fast loop. | Make `refresh_canopy_depth` public; call it at the end of `polygon_prepare` after the stand and fast reservoirs exist (no ledger → plain geometry update). Regenerate regression references; CHANGELOG before/after on day 1. | S |
| 3 | CHANGELOG placement (not an issue) | Release notes are cut from `[Unreleased]`. | **Yes.** `CHANGELOG.md:353` heads 0.2.1; its `### Fixed` at `:394` opens with the #294/#296 entry (`:396`) and the #295 entry (`:421`); `[Unreleased]` (`:15`) has `Added`/`Changed`/`Removed` and no `Fixed`. | Move both entries into a new `### Fixed` under `[Unreleased]`. | minutes |
| 4 | **#299** slow-only runs report 0 for patch-sourced rates (`bug`) | Wrong numbers, silently. Region mode requires the fast loop (`meds_config.f90:709-710`), but `example_demography` is slow-only by design (#260), so the shipped example writes zeros for disturbance, litter, recruitment and mortality carbon. | **Yes.** `w` is written only by the fast loop (`meds_fast_dynamics.f90:814, 1260, 1284`); the reader returns 0 when `w = 0` (`meds_site_diag_types.f90:546-548`); the slow rows are correctly rate × `dt_slow`. | In the slow-row loop (`meds_vegetation_dynamics.f90:275-291`): `if (.not. cfg%fast_biophysics_on) w(ip) = w(ip) + cfg%dt_slow`. The "fast-only rows read `_FillValue`" half needs a per-variable registry flag — time-box or defer. | S (weight) / M (fill) |
| 5 | **#275** `AGG_VARIANCE` reads the tick; `ground_temp_site` is layer 1 | The issue says the rename "wants a release boundary" — this is it. | **Yes.** `soil_temp_top_site` and `ground_temp_site` both source `PD_GROUND_TEMP` (`meds_output_registry.f90:411-412, 725-726`); the three variances (`:429-433`) and `leaf_temp_var_site` (`:441`) read the tick. | Drop `ground_temp_site` (or relabel "top soil layer"); relabel the four variances "variance of end-of-step samples". Skin-temperature plumbing and `PD_*_SQ` rows: defer. | S |
| 6 | **#298** restart not bitwise (`bug`) | Region restarts do not exist yet (R4), so nothing on beta depends on it. | **Partly.** `patch%adapt_dt_last` is per-patch controller state (`meds_site_state_types.f90:324`), the warm start of both integrators (`meds_fast_ark.f90:361`, `meds_fast_rk45.f90:567`), cold-started for a fresh gap, and absent from the checkpoint (`meds_io.f90:97-208`). Matches the signature. Whether it is the only missing state is unverified (`growth_hist`, `meds_site_state_types.f90:155-157`, is another candidate). | Persist `adapt_dt_last` per patch (optional on read, 0 when absent, like the `cas_can_*` group at `meds_io.f90:521`); re-run `ckpt_2yr`. | S code, M verification |
| 7 | **#270** FAST tier has no patch axis | Diagnostic gap, not wrong numbers. | **Yes.** `output_integrate_fast` handles `DIM_SCALAR/DIM_SOIL/DIM_COHORT` only (`meds_output_integrate.f90:985-993`). | New staging layout + fast cases + registry twins. | M–L → defer |

**Recommended pre-merge set:** rows 1–5, plus the code items in §0. Time-box row 6. Rows 1, 2 and 4
each need a CHANGELOG entry with before/after numbers and regenerated regression references (only
#306 moves model state; #290 and #299 move diagnostics only).

### C. Defer

**Science.** #1 equal-height shading (superseded by the two-stream when it replaces the test light
sweep); #74 condensate destination (only with `canopy_water_on`); #96 Kelvin vapour pressure; #104
`rwc_floor` decision; #154 nitrogen twin; #155 layered soil carbon; #156 CWD; #157 fire; #165 canopy
film thermal state; #167 free-convection slope (re-measured to 1.02–1.07); #180 per-layer root nodes;
#181 hydraulic redistribution; #186/#187 snow P1/P2; #254 deep thermal layers; #255 spectral radiation
record; #256 per-cohort acclimation; #257 LWdown cloud term; #258 `root_phen_factor`; #265 sub-canopy
conductance port (real, matters in gaps); #268 litter layer; #269 layered canopy air; #302/#303
forcing products and pre-2002 archive years (decided deferrals).

**Performance.** #195 allocator traffic (R5 item B8); #188 params-through-config (as #195's first
step); #196 cohort-axis threading; #158 adaptive freeze cadence; #159 soil water in the tableau.

**Structure.** #146 packed state vector and #190 typed slices (paired); #189 residue; #270 fast-tier
patch axis; #183 R3–R6.

**ROADMAP omissions to fix in the release PR** (rule 3 requires every deferred item there with an
issue number): **#265, #268, #269, #270, #275** appear nowhere in `docs/ROADMAP.md`, nor do the bugs
#290/#298/#299/#306 — whichever are not fixed pre-merge belong in §1 "Known defects". Stale on
release day: the "Deferred to v0.3.0" tags at `ROADMAP.md:86` (#104), `:109` (#167), `:239`, `:248`
(#190/#146), the "v0.3.0 question" at `:158`, "post-v0.2.x" at `:143`; §11's #197 line once closed.

### D. Closed-issue check

All ten recent closures (#294, #267, #266, #264, #260, #247, #246, #245, #241, #239) have their fix on
`beta`, each spot-checked in the tree. #294 was closed on 2026-09-27 while only on `beta`, against
CLAUDE.md's "close when the release lands" — harmless since v0.3.0 carries it. #267 is the config
change that exposed #290 (B1); #260's successor is #299 (B4). None was closed prematurely.

---

## 8. Release checklist for v0.3.0

1. Land the pre-merge set (§0 rows 1–11; §7 B1–B5) as PRs into `beta` — phased in §9 — each with its
   CHANGELOG entry and before/after numbers; regenerate the r1 regression references after #306. Then run the gfortran
   suite (§6.1–6.2) as well as ifx — it is the only second compiler on this machine — and add the
   component-section trap to CLAUDE.md.
2. `CHANGELOG.md`: move the misfiled `Fixed` entries; reword the `prep_era5land_forcing.py`
   sentences; add `## [0.3.0] — <date>`; bump version strings (§5 lists them).
3. `docs/ROADMAP.md`: add the deferred items that lack an entry (§7 C); re-date the "v0.3.0" tags;
   record the #197 decision.
4. Docs in §5 marked `[BEFORE-MERGE]`; the rules files (O3); `src/README.md` (R7).
5. Merge `beta` → `main`, tag `v0.3.0`, close #184 and #197 with the comments above, comment on #183.
6. File the `[AFTER-MERGE]` items as issues (§1 R3–R6, R9–R12; §2 F6–F12; §3 O6–O12; §4 P6–P10;
   §6.3 the three Debug-suite failures and the Campbell-with-vG validation gap; §6.4 `new_year`),
   grouped: "runtime consolidation for R3", "forcing source interface", "output serializer dedup",
   "prepare_era5 dedup", "Debug suite green".

---

## 9. Phased plan for the pre-merge fixes

Written 2026-09-27, after the review was read. It sequences §0 rows 1–11 and §7 B1–B5 into five
phases plus the release cut, ordered so that (i) every phase is independently green, (ii) the
bitwise-safe changes land before anything that moves a number, so the regression references are
regenerated exactly once, and (iii) the second compiler is in the loop from the first phase that
touches the region code.

### 9.0 Ground rules

- **Branch and PRs.** Work on `cleanup/pre-v0.3.0` (from `beta` at `e33cb56`). One commit per item,
  its message naming the review id (`fix(forcing): carry the window's first record across the year
  wrap (review F1)`), so the r1 comparison can bisect. One pull request into `beta` per phase, opened
  when the phase's gate is green; a phase is rebased onto `beta` once the previous phase has merged.
- **Gates, in every phase.** (1) ifx Release: 53/53 plus the phase's new tests. (2) ifx Debug: 50/53
  until §6.3's three pre-existing failures are fixed — record that number in the PR so a fourth
  failure is visible. (3) gfortran Release and Debug, from Phase 0 on (CMake already gives gfortran
  Debug `-fcheck=all -ffpe-trap=invalid,zero,overflow`, `CMakeLists.txt:102-107`). (4) The r1
  harness: `run_cases.sh BIN TAG` then `compare_runs.py REF NEW`, which reports every variable,
  attribute and budget line that differs — "no differences" for a bitwise-safe phase, an enumerated
  list for the others, copied into the PR. (5) The 4-cell regional smoke (`r2region/smoke`). (6) A
  CHANGELOG entry per item under `[Unreleased]` with the PR number and before/after numbers for
  anything that moves one; rule-1 comments; ROADMAP for anything deferred out of the phase.
- **Commands** (`CLAUDE.local.md` has the toolchain activation):

  ```bash
  cmake -S . -B build-ifx -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
  cmake -S . -B build-debug -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Debug   -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
  cmake -S . -B build-gfortran -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
  cmake -S . -B build-gfortran-dbg -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Debug -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
  for b in build-ifx build-debug build-gfortran build-gfortran-dbg; do cmake --build $b && ctest --test-dir $b --output-on-failure; done
  R1=~/claude_workspace/meds_runs/r1
  $R1/run_cases.sh build-ifx/meds_main <TAG>                       # six cases: bare_july_archive, bare_july_single,
  ~/miniconda3/envs/meds/bin/python $R1/compare_runs.py $R1/<REF> $R1/<TAG>   #   ckpt_2yr, demography_30yr, est_july, est_year
  build-ifx/meds_main ~/claude_workspace/meds_runs/r2region/smoke/region_smoke.toml   # outputs go where its [output].dir says
  ```

| Phase | Content | Moves numbers? | Effort |
|---|---|---|---|
| 0 | Baseline references, gfortran back in the build, the counts | no | ½ day |
| 1 | Crashes and correctness: year wrap (F1), gfortran `bufs(:)` (6.2), region I/O exit path (O1), F5, R1, R8, (R2) | no — bitwise gate | 1 day |
| 2 | Ledger and diagnostics: #290, #299, #275, FAST stamp (O2), F2, O5 | diagnostics only | 1 day |
| 3 | Model state and checkpoint: #306, #298 time-box | yes, once | ½–1 day |
| 4 | Documentation, tools, examples: §5 `[B]`, R7, F3–F4, O3–O4, P1–P5 | no | 1 day |
| 5 | Release cut | no | ½ day |

### 9.1 Phase 0 — baseline and a second compiler

1. **Wrap the six lines over column 132** (§6.1): `meds_fast_be_stage.f90:430`, `meds_fast_frozen.f90:426`,
   `test_column_ark.f90:519, 549`, `test_column_dynamics.f90:525`, `test_column_rk45.f90:738`. Add
   `find src test -name '*.f90' | xargs awk 'length > 132 {print FILENAME": "FNR}'` to the PR
   description as the check (and to CI when there is one). ifx output is bitwise unaffected.
2. **Baseline references.** Build `e33cb56` (plus item 1) with ifx Release and run
   `run_cases.sh build-ifx/meds_main ref_e33cb56`. The existing `ref/` is `beta` at `2bb8bbd`, from
   before #305 changed the forcing height, so it no longer serves; keep it, do not overwrite it.
   Record for the CHANGELOG entries to come: `energy fails` and the `resid_energy_site` range in
   `est_july`, `bare_july_*` and the smoke run (#290, expected ≈ −0.3…−4 W/m²); day-1
   `cas_depth_patch` (20 m) and `wind_cas_top_patch` in `bare_july_single` (#306); the
   `demography_30yr` disturbance/mortality/litter rows (all 0, #299).
3. **gfortran state.** Configure and run all four builds above. Expected after item 1: gfortran
   Release 52/53 and Debug with the `region` segfault (fixed in Phase 1), ifx Debug 50/53. Open the
   issue for §6.3's three Debug failures now, so the gate numbers have a home.
4. **Counts that are measured, not designed:** `CLAUDE.md:46` → 53 tests, ~22 s; the three Debug
   failures noted beside the Debug recipe until the issue closes.

Gate: ifx Release 53/53; gfortran Release builds; `ref_e33cb56` saved; the baseline numbers written down.

### 9.2 Phase 1 — crashes and correctness that move no numbers

Every item here must leave the six r1 cases and the smoke run **bitwise identical** to `ref_e33cb56`.

1. **F1, the year wrap for anchors other than 01:00** (`meds_met_driver.f90:292-318, 404-450, 700-717, 888-918`).
   Design: the window's first record always lies in axis month 1 (`:782`), so after `load_axis_month(src, 1)`
   at open copy it into a second carry slot, `src%seam_first(cell, var)`, beside `carry`; let
   `locate_record` return a second sentinel for `irec == irec_cycle_first` when its month is not
   loaded, and `archive_value` read it. `met_prefetch`'s wrap branch stays (it carries
   `irec_cycle_last` when month 1 is the loaded month); the mid-step `error stop` remains the
   programming-error assertion it was meant to be. Delete the "01:00 only" assumption nowhere
   stated — and state the real rule in `configuration.md` ("any record stamp; the seam falls inside
   a step").
   Test: `meds_test_era5land_archive.f90` writes every month of 2021 with the year hard-coded in
   `write_month` (`:149-151`); parametrize the year and add **January 2022**. Then, in
   `test_met_era5land.f90` test 6 (`:283-312`), repeat the wrap walk for `recycle_start =
   2021-01-02 00:00` (window to 2022-01-02 00:00) and `2021-01-01 06:00` (to 2022-01-01 06:00): the
   seam half-hour must read `0.5·(last + first)` as the 01:00 case does at `:294-297`, with one
   month load, not an `error stop`. Region mode: `region_step_month`'s assertion that no month loads
   inside a model month (`meds_region.f90:205-219`) must still hold — the seam record now comes from
   the slot, not a load.
   Alternative if time is short: reject non-01:00 anchors for `format = "era5land"` in
   `validate_config` with a message, and say so in `configuration.md`. Do not ship the silent crash.
2. **6.2, the gfortran `region` segfault.** Design (a) of §6.2: `type(output_buffers_t), allocatable :: bufs(:)`
   owned by `meds_region_t`, `type(output_buffers_t) :: bufs` owned by `meds_run_t`, and `out_bufs`
   removed from `meds_polygon_t`; `polygon_prepare`, `polygon_step` (its `stepper` and
   `tick_output`, `meds_polygon.f90:177-197`) and `polygon_report` take `bufs` as an argument; the
   detail set stays on the polygon (`detail_files`/`detail_bufs` are scalar allocatable
   components, never a strided section). `io_phase` becomes
   `output_serialize_region(reg%out_files, reg%bufs)` and `region_finalize`
   `output_region_close(reg%out_files, reg%bufs)`: contiguous, no temporary. This is also the shape
   R3's threaded loop indexes (`bufs(p)`). Fold in `new_year`'s initialization before the loops
   (`meds_region.f90:239`, §6.4). Add to CLAUDE.md's trap list: "never pass a derived-type component
   section `a(:)%c` whose type has allocatable components as an actual argument; gfortran builds a
   temporary and its copy-out dangles the components (review 2026-09-27 §6.2)."
   Test: the existing `region` test under gfortran Release and Debug (`-fcheck=all`) is the test;
   ifx output bitwise unchanged.
3. **O1, the region I/O exit path** (`meds_output_manager.f90:64-80`, `meds_output_stream.f90:466-482`).
   Design: drop the equal-count assertion and the special case with it — write records
   `i = 1 .. maxval(bufs(:)%queue(t)%n)`; a polygon without record `i` contributes the fill value
   (`region_write_one` packs per polygon, so this is one branch per variable); the calendar of
   record `i` comes from the first polygon that holds it; `same_time` still checks the polygons that
   do. On the normal path every polygon holds every record and nothing changes; on the failure path
   what closed is written and the failed polygon's slice is fill from its last record on — which is
   exactly R5's steady state for an isolated failed polygon, so R5 inherits the writer unchanged.
   Test: in `test_region`, after a healthy month, poison one polygon (set a cohort's `bleaf` to
   `ieee_value(…, ieee_quiet_nan)` in `reg%poly(2)%site` before a mid-month step) so the NaN guard
   trips; assert `region_step_month` returns the failing status, the daily file holds the days the
   other polygons closed, and polygon 2's slice is fill after its last good day.
4. **F5, `met_open` leaks `src%ncid` on the `validate_file_against_config` rejection** (`meds_met_driver.f90:190-196`):
   use `met_close(src)` as the other three rejection paths do. Test: `test_met_driver.f90:351-368`
   calls `met_close` after the rejection and checks it does not error-stop.
5. **R1, `polygon_prepare` does not reset the restructure flags** (`meds_polygon.f90:107-112`,
   `meds_driver.f90:114`): reset `restructure_pending`/`restructure_new_year` in the counter-reset
   block before the `select case`, so the restart branch can still set them. Test: give the C-API
   slot-reuse test (`test_c_api_run`) an `end_time` on the 1st, then a second bare-ground run in the
   same slot; its first step must not call `advance_boundary` (assert on the ledger's phase counter
   or on `restructure_pending` after `driver_open`).
6. **R8, `detail_polygons` validated after building every polygon** (`meds_region.f90:115-135`):
   move the membership check to right after the id loop, and `met_close` on that return.
7. **R2, the output-file-set sequence written three times** (`meds_driver.f90:214-226`,
   `meds_region.f90:139-164`) — optional in this phase, but item 2 rewrites `region_open` anyway:
   one `open_output_files(files, cfg, prefix, soil, restrict_region, verbose)` beside
   `polygon_prepare`, with `apply_io_overrides` and `ensure_output_dir` moved out of the site driver
   so `meds_region` no longer `use`s `meds_driver`.

Gate: `compare_runs.py ref_e33cb56 phase1` reports no differences; smoke output identical; ifx
53/53 + the new tests; **gfortran Release and Debug 53/53**.

### 9.3 Phase 2 — ledger and diagnostics

Numbers move here only in diagnostics: the r1 comparison must list exactly the variables named
below and nothing in the state or flux records.

1. **#290, the bottom face in the ledgers** (§7 B1). Design: the soil solve already evaluates the
   committed bottom face, `bottom_heat_face` at `meds_soil_energy.f90:182` (implicit) and `:256`
   (explicit), positive up; `energy_flux_t%bottom_heat = −hf(n)` (`:204`) exists but is set only on
   the split path. Return the committed face from both solves and book it in the three ledger sites
   instead of `frozen%hydrology%geothermal` (`meds_fast_frozen.f90:346`, always 0):
   `bf%soil_enth_in` (`meds_fast_be_stage.f90:187`), `bf%whole_enth_in` (`:200`; a signed term on
   the "in" side is enough, the ledger sums in − out), and RK45's `e_in` (`meds_fast_rk45.f90:823`).
   The ARK soil sub-ledger (`meds_fast_ark.f90:560`) reads the same `bf`, so it closes with it. No
   physics moves. Measure first, per CLAUDE.md: on `est_july`, print one step's `resid` beside
   `−g_deep·(T_bot − deep_temp)`; they should agree to the ledger tolerance. If after the fix
   `n_fail` is not 0, a second leak exists — stop, measure the residual's dependence on the same
   quantity, and file it rather than widening the tolerance.
   Acceptance: `est_july`, `bare_july_*`, `est_year` and the smoke run report `energy fails = 0`;
   `resid_energy_site` before/after in the CHANGELOG (§6.5's 17856/17856 → 0 for the smoke run);
   everything else bitwise. Issue #290 closes with the release.
2. **#299, slow-only runs weight patch rows with `w = 0`** (§7 B4): in the slow-row loop
   (`meds_vegetation_dynamics.f90:275-291`), `if (.not. cfg%fast_biophysics_on) w(ip) = w(ip) + cfg%dt_slow`.
   Acceptance: `demography_30yr`'s disturbance, litter, recruitment and mortality rows are non-zero
   (before: 0; after: the values, in the CHANGELOG). The "fast-only rows read `_FillValue` in a
   slow-only run" half stays open on #299 with its scope narrowed in a comment.
3. **#275, the release-boundary renames** (§7 B5): drop `ground_temp_site` (it duplicates
   `soil_temp_top_site`, `meds_output_registry.f90:411-412, 725-726`; the example README's `:344`
   already points users at `soil_temp_top_site`, and no script reads the dropped name); relabel the
   four variances (`:429-433, :441`) "variance of end-of-step samples". Regenerate
   `meds_io_config.toml` with `meds_main --dump-io-config` (ctest `io_config_example` enforces the
   match) and `diagnostics.md`'s counts (§5.2). Skin temperature and `PD_*_SQ` stay on #275.
4. **O2, the FAST tier's stamp** (`meds_fast_dynamics.f90:490-498`, `meds_output_stream.f90:193-203`).
   Decide once: the recommendation is the period start, consistent with #294 — stage
   `fast_time(isub) = step_start + (isub−1)·dt_fast` and keep `t_sample` for the met lookup and the
   probe; fix the "sub-step midpoint stamps" comment (`meds_output_types.f90:336`). The `time`
   coordinate of every `-F-` file moves by `−forcing_sample_frac·dt_fast` (7.5 min at the defaults);
   the three `examples/example_biophysics/plot_*.py` read `-F-` files, so re-run them in Phase 4. If
   the output change is unwanted now, fix the `long_name` instead and file the stamp as an issue —
   but do not ship a label that contradicts the value.
5. **F2, the duplicated LW-synthesis + ρ block in `met_instant`'s CONST return** (`meds_met_driver.f90:510-529`
   vs `:584-600`): wrap only the file-dependent interpolation in `if (backend /= CONST)`, one exit.
   `test_fast_loop` block 8 (`test/test_fast_loop.f90:463`) exercises the CONST path: it, and any r1
   case, must stay bitwise; if the dedup changes a CONST number, the duplicate was load-bearing —
   find out why before deleting it.
6. **O5, the disturbance test's seed** (`test/test_disturbance.f90`): seed `PD_MORT_C_CULL` with
   `w = 0`, the case production produces since #297; reword the two comments (`meds_demography_patch_fusefiss.f90:603-607`,
   `meds_site_diag_types.f90:511-518`).

Gate: the r1 difference list is exactly {`resid_energy_site` (fast cases), `ground_temp_site`
absent, four `long_name` attributes, `-F-` `time` coordinates (if O2's stamp was chosen),
`demography_30yr`'s patch rows}; ifx and gfortran suites green; smoke `energy fails = 0`.

### 9.4 Phase 3 — model state and checkpoint completeness

The one phase that moves the model. Do it alone, so its difference list is its own.

1. **#306, the canopy-air depth on day 1** (§7 B2). `refresh_canopy_depth(site, cfg, ledger)`
   (`meds_slow_dynamics.f90:108`) takes the ledger as an optional argument and without it is a plain
   geometry update through `cas_set_depth(cas, depth_new)` (`:139-142`), intensive state invariant.
   Make it public and call `refresh_canopy_depth(poly%site, cfg)` at the end of `polygon_prepare`'s
   fast block (`meds_polygon.f90:114-126`), after `init_fast_reservoirs`/`seed_snow` and for every
   init mode: a bare-ground stand gets the 5 m floor instead of the 20 m type default, an
   established one its `h_top + freeboard`. Test: after `polygon_prepare` on bare ground,
   `patch%cas(1)%depth == cfg%aero%min_canopy_depth`; on a census stand, `h_top + freeboard`.
   Acceptance: `bare_july_single` day 1 — `cas_depth_patch` 20 → 5 m, `wind_cas_top_patch` and the
   day-1 fluxes in the CHANGELOG; the July monthly means and `est_july` (whose depth is refreshed by
   the first slow step anyway) should move little — record the largest relative change.
2. **#298, restart exactness, time-boxed to one day** (§7 B6). Two per-patch quantities are missing
   from the state file: `adapt_dt_last` (the integrators' warm start, `meds_site_state_types.f90:324`)
   and `can_depth` itself (item 1 shows the resumed run re-derives it, but only at the next boundary
   unless `polygon_prepare` refreshes it — after item 1 it does, so persisting the depth is about
   exactness under a non-boundary checkpoint). Persist both, optional on read with the
   `cas_can_*` pattern (`meds_io.f90:521`, absent → 0 / re-derived), so old state files still load.
   Run `ckpt_2yr`: compare the resumed run against the continuous one with `compare_runs.py`. If
   they now match, #298 closes with the release; if a residual remains, list the differing
   variables in ROADMAP §1 under #298 and stop — `growth_hist` (`meds_site_state_types.f90:155-157`)
   is the next candidate, for R4.
3. **Regenerate the references**: `run_cases.sh build-ifx/meds_main ref_phase3` becomes the
   reference for Phase 4 and for the release; keep `ref_e33cb56` beside it.

Gate: the difference list against `ref_e33cb56` is explained by day 1 (item 1) and the checkpoint
variables (item 2), with the numbers in the CHANGELOG; `ckpt_2yr` resumed == continuous, or the
residual is documented; both suites green; smoke completes.

### 9.5 Phase 4 — documentation, tools and examples

No Fortran behaviour changes. Grouped by file so each is one commit; the review's ids in brackets.

1. **Source comments that narrate history (rule 1):** `meds_driver.f90:3-11, 88-90, 212-213, 238-240`;
   `meds_polygon.f90:43-44`; `meds_c_api_run.f90:93` (status 4); `CMakeLists.txt:352-357` [R7];
   `meds_met_driver.f90:39, 186-189, 325-330, 402, 416-418, 623-628, 650-653`,
   `meds_forcing_config.f90:6, 33, 48-52, 99, 126-139` [F4, nits]; `meds_output_integrate.f90:134-136, 500-503`,
   `meds_site_diag_types.f90:24-27, 119-125` [O3]; `meds_output_config.f90:105-112`,
   `meds_output_manager.f90:3-4` [nits]; `download_era5land_cds.py:4` [P5].
2. **Rules and orientation:** `.claude/rules/output.md:48, 58` and `state-demography.md:69, 76-77`
   (the tick runs *before* the restructuring, which is between steps; the tendency-bundle rule's
   reason is the daily `sort_cohorts`) [O3]; `src/README.md` (`main/`, `forcing/`, `io/`, `config/`
   rows and the header counts — name the six new modules) [R7]; `src/forcing/README.md:7-8, 72-73, 103-105`
   [F3, P10]; `CLAUDE.md` (the trap from Phase 1, if not already there; the Debug note from Phase 0).
3. **Science and configuration docs:** `forcing.md:67-68, 114-116, 600` and the `era5land` rule
   (`dt_slow = 1d`, midnight start) [F3, §5.3]; `configuration.md:121, 224-225` [§5.1];
   `order_of_processes.md:20, 220` (`driver_open`) [§5.3]; `diagnostics.md:212, 280-288` (regenerate
   from `--dump-io-config` after Phase 2) [O4]; `ed2_comparison.md:63, 101, 321, 367, 378` and its
   version pin [§5.4]; `README.md:56, 61, 121` and a feature row for regional runs, the ED_ERA5land
   archive and prescribed CO₂ [§5.4].
4. **`docs/dev_plans/`:** the archive audit of §10 — six moves now (biogeochemistry, snow, GPU
   evaluation, veg-energy plan, the 2026-09-13 docs review, the v0.2 release plan), the integrator
   plan after its two issue actions, the structure design to the Reference table, four header
   refreshes, the README tables and count, and the bare-path links §10.4 lists [§5.4.7, §10].
5. **CHANGELOG hygiene** (the `[0.3.0]` heading waits for Phase 5): move the #294/#296 and #295
   `### Fixed` entries from `[0.2.1]` (`:394-421`) into a new `### Fixed` under `[Unreleased]`;
   reword `:115, :153, :307-308`; `:74` → 253; add the missing PR numbers to the #281–#293 entries [§5.5].
6. **`docs/ROADMAP.md`:** entries for #265, #268, #269, #270, #275 (narrowed), #298 (if residual),
   #299 (narrowed); retarget `:86, :109, :158, :239, :248`; decide the `[io]` shim (`:143`,
   `meds_config_io.f90:1151-1154`): remove it in 0.3.0 (the shim was "post-v0.2.x") or file an
   issue; record #197's decision in §11 [§5.5, §7 C].
7. **Python tools:** `make_forcing_file.py` writes `u10`/`v10` beside `Wind` and drops `U_MIN` (the
   Fortran floors) [P4]; the longitude docstrings (`build_era5land_archive.py:9`,
   `build_era5land_static.py:6`, design §14.2 — the axis is −179.9 → 180.0; do **not** change the
   `>`) [P1]; `forcing.md:30-42, 66, 76` CDL and the wind floor [P4, P10]. P3 (the de-accumulation
   and unit tables in one place) only if the phase has time; otherwise it is the "prepare_era5
   dedup" issue.
8. **The example, last:** rerun `examples/example_biophysics/run_example.py` (both stages; the
   50-year spin-up is ~10 min on a compute node) with the Phase 3 binary, refresh the README's
   numbers (`:193-195`) and the four PNGs, add the lapse/height note to "Notes on the configuration"
   (`:284-330`) [P2].

Gate: the audit's checks come back clean — `grep -rn "ensure_month\|driver_init\|apply_met_to_ctx\|prep_era5land_forcing" docs src .claude CLAUDE.md README.md`
(historical hits in CHANGELOG and `archive/` excepted); `meds_main --dump-io-config` equals
`meds_io_config.toml` (ctest `io_config_example`); `diagnostics.md`'s counts equal the registry's;
every `[Unreleased]` entry cites a PR; ROADMAP has an issue number on every deferred line it adds.

### 9.6 Phase 5 — the release cut

1. `CHANGELOG.md`: `## [0.3.0] — <date>` over the `[Unreleased]` block, the compare link `:1858`.
2. Version strings: `CMakeLists.txt:17`, `python/pyproject.toml:16`, `python/meds/__init__.py:19`;
   `README.md:61, 121`.
3. Final gates on the release candidate: all four builds and suites; `compare_runs.py ref_phase3 rc`
   with no differences; the smoke run; `python -m meds.plant` round trip.
4. Pull request `beta` → `main`; tag `v0.3.0`; close #184, #197 (§7 A), #290, #306, #299 (narrowed)
   and #275 (narrowed) with one-line comments naming the PR; comment on #183; #298 per Phase 3's
   outcome.
5. File the after-merge issues (§8 item 6) and update the untracked `CLAUDE.local.md` baseline.

### 9.7 Deliberately outside the plan

The `[AFTER-MERGE]` consolidations (§1 R3–R6, R9–R12; §2 F6–F12; §3 O6–O12; §4 P6–P10) — they
are R3's shape and belong with it; the three Debug-suite failures (§6.3, an issue from Phase 0);
#270; the fill-value half of #299; the skin-temperature half of #275. If a phase runs long, the
items to drop first are, in order: R2 (Phase 1.7), P3 (Phase 4.7), the O2 stamp change (keep the
label fix), and the #298 time-box.

---

## 10. `docs/dev_plans/`: what can move to `archive/`

Written 2026-09-27. The directory's rule (`docs/dev_plans/README.md`): a document stays only while it
has open items still intended, or is the as-built description that source comments cite by section
("Reference"); everything else goes to `archive/` with a tombstone. Section numbers are never
renumbered and most citations are by bare name, so a move breaks nothing but path links.

Method: two readers over the thirteen documents in the live directory, every open item the document
itself still presents traced to an open issue, a `ROADMAP.md` line, another file, or the PR that
shipped it; citations counted from `src`, `test`, `CMakeLists.txt`, `.claude`, `python`, `scripts`;
the load-bearing claims re-checked by hand (issue states with `gh`, the DAMM kernel, every
bare-path link). Line numbers are the `cleanup/pre-v0.3.0` tree.

### 10.1 Verdicts

| Document | README table | Open items → where they live now | Cited from code | Verdict |
|---|---|---|---|---|
| `MEDS_BIOGEOCHEMISTRY_DESIGN.md` | Live | P1 nitrogen #154, P2 vertical pools #155, CWD #156, fire #157 (ROADMAP §3). **P1 DAMM: the kernel was deleted as unreachable (#153, closed 2026-09-14; branch `archive/damm-hr`) — the header still lists it.** | 1 by § (`test_soil_biogeochem.f90:4`, the §8 test plan), 2 bare | **ARCHIVE** — live description: `docs/science/soil_carbon.md`, `src/slow_dynamics/soil/README.md` |
| `MEDS_SNOW_DESIGN.md` | Live | P1 multi-layer #186, P2 interception #187 (ROADMAP §9). Unfiled: §15 Q1–Q2 (`rho_snow`, `snowfac` calibration) — file or drop. | 3 bare "P0" (`meds_biophysics_opts.f90:119`, `meds_ground_biophysics.f90:11`, `meds_soil_types.f90:168`) | **ARCHIVE** — live: `docs/science/snow_biophysics.md`, `snow_params_t`; §13's keys differ from as-built |
| `MEDS_GPU_EVALUATION.md` | Live | §12: three done (NUMERICS §7 refuted, #194, `building.md:101-106`), #195, #196, #104, #197 (decided 2026-09-13, still open — §7 A), #183. README's "five of seven" overstates: two done, two decided. | 2 bare paths (`CMakeLists.txt:14`, `src/README.md:88`) | **ARCHIVE** — the verdict is restated in `docs/building.md` and the CMake header |
| `MEDS_VEG_ENERGY_INTEGRATION_PLAN.md` | Reference | #165, #167 (ROADMAP §5). Header still lists "retire `veg_energy_step_implicit`" (#166 closed, PR #120) and "wood sizing" (#168 closed, PR #125). Unfiled: §3's three allocation/deep-copy refactors (fold into #195) and §6 sunlit/shaded (drop). | 1 (`meds_fast_types.f90:430`, rationale) | **ARCHIVE**, ⚰️ for §9–§11, §14 — live: `docs/science/vegetation_energy_dynamics.md`, `numerical_scheme.md` §5a′ |
| `MEDS_DOCS_REVIEW_2026-09-13.md` | Reviews | Executed (29 moves, README/`src/README.md`/`CLAUDE.md` split, ROADMAP and CHANGELOG created) except its §5.3 comment sweep — 156 history-narrating comment lines remain in `src`+`test`, the same rule-1 residue §9.5 item 1 lists — and the §5.4 page items. | 0 | **ARCHIVE** — its outputs are the live files |
| `MEDS_V02_RELEASE_PLAN.md` | *none* | v0.2.0 shipped 2026-09-14 (PRs #203–#259). Never decided: milestone/labels. Never filed: the RK45 frozen-record oracle (§2.2), a fast-loop golden (§5.1; the six r1 cases may cover it). | 0 | **ARCHIVE** — CHANGELOG `[0.2.0]`, ROADMAP |
| `MEDS_PRODUCTION_INTEGRATOR_PLAN.md` | Live | #158, #159, #104; RK45 production warning shipped (#160, `meds_config.f90:684-685`) — header stale. **ψ_leaf/N2b: #162 was closed 2026-09-14 without a comment while ROADMAP:90 ("open question"), CHANGELOG:520 ("#162, open") and `numerical_scheme.md:292, 464` ("tracked as N2b in this plan") say open.** E5 (RK45 rescue snapshot, `:1999`): tracked nowhere. | 3 (`meds_soil_water.f90:665`, `meds_fast_types.f90:369`, and a path inside an error string at `meds_fast_frozen.f90:442`), all history | **ARCHIVE once #162 and E5 are re-homed** (reopen #162 or record why it closed; file or drop E5) — live: `docs/science/numerical_scheme.md` |
| `MEDS_CODE_STRUCTURE_DESIGN.md` | Live | §15 closed: Phase 1 #141, Phase 2 #145 → #254, Phase 3 #201, Phase 4 #147 → #146, Phase 5 → #188/#189/#190, Phase 6 #191 (`meds_test_assert.f90`, no local `check` left); decision #13 reversed (#193, `test/` is flat). | ~20 by decision/rule number: `CMakeLists.txt:132, 146, 171` (decisions #9, #2, #4), `python/pyproject.toml:7` and the FFI shims (#1), `meds_canopy_types.f90:12`, `meds_leaf_opts.f90:9`, `meds_allometry.f90:47` (#8, #10, #12), `meds_fast_config.f90:9`, `meds_demography_rates.f90:8` (rules 5, 8) | **→ REFERENCE** — placement rules 2, 5, 6, 7 and decisions #8–#12 exist only here (`src/README.md`/`CLAUDE.md` carry rules 1, 3, 4, 8 and the graph); §7.2 and #13 have drifted |
| `MEDS_FROZEN_SEAM_CONTRACT.md` | Reference | The note *is* #201's deliverable (PR #208). §4's "make debit-before-credit an assertion" never filed. `CODE_STRUCTURE` §15.4 still says "write it down". | **0** — the README's "cited by section from the code" does not hold | **REFERENCE for now** — the only statement of the rate-seam and arbitration rules; fold §2–§5 into `soil_carbon.md` §7 / `order_of_processes.md`, then **ARCHIVE** |
| `MEDS_NUMERICS_SCOPING.md` | Reference | Header lists MB2 (#163 closed, knobs deleted) and §11.3 bare arrays (#164 closed premise-false) as open; §8a/§8d–§8g, §9 #5, §12.2 all shipped or rejected. Unfiled: the §8b L2 "enforce conservation everywhere" sweep (nearest: #189, #290). | 20 files by § (§5.1, §5.3, §11, §8e, BB1, QW2/QW4) — as-built mechanisms | **REFERENCE (earned)** — refresh the header |
| `MEDS_IO_DESIGN.md` | Reference | Header "variance deferred" — #174 shipped (residue #275); async writer / `NC_FLOAT` / `AGG_WMEAN` dropped by decision; flush-order lines `:607-608, 652, 676, 1122` amended only at `:681-688` (#294). | 10 by § (registry, integrators, serializer; `.claude/rules/output.md:12`) | **REFERENCE (earned)** — refresh the header, pointer at `:607` |
| `MEDS_IO_V01_PLAN.md` | Reference | Its deferred list shipped whole in v0.2.0 (#169 #170 #171 #172 #173 #174 #175); the `[io]` shim removal (ROADMAP:143) has no issue; unfiled and low value: `_fast` suffix retirement, `add_variable_family`, the FvCB rate rows. | 13 by § (stages, reducer, diag blocks, 2-D axis, DBH binning) | **REFERENCE (earned)** — rewrite the header |
| `MEDS_POLYGON_RUNTIME_PLAN.md` | Live | R3–R6 → #183 (+ #290, #195); in-plan only: B9 nvfortran check, the memory split "before R3" (§10.3.1), OR2. #183's body asks for MPI, which the plan dropped. | 24 by § from `src`/`test`/CMake | **LIVE** — the as-built description of region mode and the plan in execution |
| `MEDS_CODE_REVIEW_2026-09-27.md` | Reviews | this document | 0 | LIVE until §9 is executed |

Net: **six moves now** (biogeochemistry, snow, GPU evaluation, veg-energy plan, the 2026-09-13 docs
review, the v0.2 release plan), **one more after two issue actions** (integrator plan), **one
re-classification** (structure design → Reference), **one fold-then-move** (seam contract), **four
header refreshes** (numerics scoping, IO design, IO v0.1 plan, structure design), and the
archive count goes 34 → 40 (41, 42). The live table is then the polygon plan alone, which is what
"live" should mean.

### 10.2 Tombstone drafts

House style: `> # 🗃️ ARCHIVED — <date>. <one line>` then *shipped* (PRs), *deferred* (issues),
*live description*, *do not trust*. ⚰️ where a section is actively wrong.

- **Biogeochemistry.** 🗃️ P0 and P3 shipped; the rest is issues. Shipped: pools, CENTURY matrix,
  EXPM, SASU (#35); per-patch state, `[soil_carbon]`, restart, the litter → daily step → fast-Rh
  seam (#64); `soil_carbon_on` default (#143); Λ/lignin audits (#141). Deleted: the DAMM kernel
  (#153). Deferred: #154 nitrogen, #155 vertical pools, #156 CWD, #157 fire (ROADMAP §3). Live:
  `docs/science/soil_carbon.md`, `src/slow_dynamics/soil/README.md`. Do not trust: §5.1/§10 paths,
  §9 "Ra is still 0", §7 P1 "DAMM (already in `meds_column_co2`)".
- **Snow.** 🗃️ P0 shipped; P1/P2 are issues. Shipped: the single-layer store, `snowfac`-ramped
  optics/BC/latent split, paired melt transfer, closed ledgers (#42); the shared snow stage on ARK
  and RK45, `[fast].snow_on` deleted (#77, #80); ice-curve sublimation (CHANGELOG:1211). Deferred:
  #186, #187 (ROADMAP §9); §15 Q1–Q2 never filed. Live: `docs/science/snow_biophysics.md`,
  `snow_params_t`. Do not trust: the 2026-07 status block, §6 "ARK stays snow-free", §2.2 module
  names, §5's `file:line`, §13's keys.
- **GPU evaluation.** 🗃️ The measurement stands; every recommendation is done, decided or an
  issue. Measured (#110): offload 1.4× slower, one kernel at 0.4 % occupancy, GPU-as-20-cores 26×
  slower. Done: BB2/BB3 refuted in NUMERICS §7; build/CLAUDE overselling fixed (#194);
  `MEDS_GPU=gpu` demoted. Decided: `wp` stays `real64` (#197). Deferred: #195, #196, #104; the
  regional axis is `MEDS_POLYGON_RUNTIME_PLAN.md` (#183). Do not trust: §2 counts, §13 paths.
- **Veg-energy integration plan.** 🗃️ §1–§8, §12–§13 shipped (store on 2026-07-31; selectors
  deleted, PR #88 era). ⚰️ §9–§11 and §14 were overturned the same day by PR #90 — one
  coefficient, refreshed per stage; default 900 s. Open: #165 film store, #167 slope (ROADMAP §5);
  #166, #168 closed. Live: `docs/science/vegetation_energy_dynamics.md`; stability by stand height
  in `numerical_scheme.md` §5a′.
- **Docs review 2026-09-13.** 🗃️ Executed 2026-09-13/14 (commit `bf85b45`; PRs #203–#208, #235,
  #253): 29 plans archived, README / `src/README.md` / `CLAUDE.md` split, ROADMAP and CHANGELOG
  created, three science pages written. Not executed: the §5.3 comment sweep and the §5.4 page
  items — now §9.5 item 1 of the 2026-09-27 review. Line numbers are the `b9596c5` tree.
- **v0.2 release plan.** 🗃️ v0.2.0 shipped 2026-09-14 (tag `v0.2.0`, PR #259; phases 0–6 in PRs
  #203–#253; CHANGELOG `[0.2.0]`). Deferred to v0.3+ with issues (ROADMAP); #164 later closed
  premise-false. Never decided: milestone/labels. Never filed: the RK45 frozen-record oracle
  (§2.2), the fast-loop golden (§5.1).
- **Production integrator plan** (after the #162/E5 actions). 🗃️ The stability question is
  closed; the remainder is issues. Shipped: N2a per-stage conductance refresh + 900 s (#90), N2b
  corrector (#91), E3/E1b thrash detector + E4 (#105), §7 C1–C5 patch threading (#107, #109),
  RK45 > 300 s warning (#160), the stomatal feedback that dissolved N2d/N2f (#95, #98). Refuted by
  measurement (read as findings): §1i.1, N1, N2e, N6, E1 as written, E2. Deferred: N5 → #158;
  tableau → #159; ψ clamp → #104; ψ_leaf → #162 (re-homed); E5 → issue or dropped. Live:
  `docs/science/numerical_scheme.md`; kernel convention NUMERICS §11. Do not trust: §9 paths, the
  §5 soil-T row (`:895-900`), anything marked SUPERSEDED (§1h, §2a).
- **Structure design** (header, not a tombstone). 📚 REFERENCE — §15 is closed (Phase 1 #141,
  Phase 2 #143/#145 → #254, Phase 3 #201, Phase 4 #147 → #146, Phase 5 → #188/#189/#190, Phase 6
  #191); decision #13 reversed (#193). What stays live: §1 decisions, §5 straddlers, §6 placement
  rules, §7.6 #3–#4 — cited by number from ~20 files. Do not trust §7.2 (the package is
  `meds.plant/.demography/.model`) or decision #13. Live tree: `src/README.md`.

### 10.3 Headers to refresh (Reference documents that overstate what is open)

- `MEDS_NUMERICS_SCOPING.md:20-21`: drop MB2 (#163 closed) and §11.3 (#164 closed); file or drop
  the §8b L2 sweep.
- `MEDS_IO_DESIGN.md:15`: variance shipped (#174; residue #275); add "see `:681-688`" at `:607`.
- `MEDS_IO_V01_PLAN.md:8-17, 45-50`: nothing deferred remains; point at #255, #275 and the shim line.
- `MEDS_FROZEN_SEAM_CONTRACT.md`: say it is #201's deliverable; mark `CODE_STRUCTURE` §15.4
  "delivered (#201)".

### 10.4 What the moves touch

Bare-path links that break on a move (all verified; historical hits in `CHANGELOG.md` and
`archive/` excepted, but the two CHANGELOG links are one word each to fix):

| Link | Document |
|---|---|
| `CMakeLists.txt:14`, `src/README.md:88`, `docs/building.md:106`, `docs/ROADMAP.md:257` | GPU evaluation |
| `src/slow_dynamics/soil/README.md:48`, `docs/ROADMAP.md:57` | biogeochemistry |
| `docs/science/snow_biophysics.md:16`, `docs/ROADMAP.md:220` | snow |
| `docs/ROADMAP.md:104`, `examples/example_biophysics/README.md:229` (cites the *overturned* §10 — add the note) | veg-energy plan |
| `src/fast_dynamics/numerics/meds_fast_frozen.f90:442` (inside an error string), `docs/ROADMAP.md:74`, `docs/science/numerical_scheme.md:292, 464`, `docs/science/canopy_aerodynamics.md:124`, `MEDS_NUMERICS_SCOPING.md:42`, eight `archive/` banners | integrator plan |
| `CHANGELOG.md:1479` | seam contract (later) |
| `CHANGELOG.md:1501` | v0.2 release plan |

`docs/dev_plans/README.md`: delete the Live rows for the four moved plans; move the structure
design to the Reference table; delete the Reference row for the veg-energy plan and the Reviews row
for the docs review; add a Reviews row for this document; set the `archive/` count (34 → 40) and add
the new tombstone grades. Gate for the phase: every `.md` in `dev_plans/` appears in exactly one
README table, `ls archive | wc -l` equals the count, and `grep -rn "dev_plans/MEDS_" --include=*.md --include=*.f90 --include=*.py --include=*.txt .`
resolves for every non-historical hit.

### 10.5 Stale cross-references found on the way (outside `dev_plans/`)

- `docs/ROADMAP.md:90` calls #162 an open question and `:96` lists #164 as open — both closed
  2026-09-14/15. `CHANGELOG.md:520` "(#162, open)"; `docs/science/numerical_scheme.md:292, 464`
  "tracked as N2b". Decision needed: reopen #162 (the ψ_leaf non-convergence is still a fact) or
  record why it was closed, then fix the four references.
- `docs/science/soil_carbon.md:250-253` describes `heterotrophic_respiration_damm` and its
  dispatcher as existing; the kernel was deleted (#153) — only an enum mention remains in
  `meds_biogeochem_types.f90`.
- `MEDS_CODE_STRUCTURE_DESIGN.md` §15.4 still says "write it down" for the seam contract that #201
  delivered.
- #197 is open with the decision made (§7 A); ROADMAP:263 still says *Candidate*.
- The `[io]` shim removal (ROADMAP:143) remains the one deferred item without an issue (§5.5).

All of this fits Phase 4 (§9.5 item 4): the moves and header refreshes are mechanical (½ day);
folding the seam contract into a science page and the #162/E5 decisions are the only parts that
need thought.
