# ERA5 forcing preparation — the `scripts/prepare_era5/` upgrade plan

> # 📐 DESIGN — written 2026-09-25. No code yet.
>
> Plans the reorganization of the ERA5-Land preparation scripts into `scripts/prepare_era5/`, and
> the new **continental** workflow: download a whole continent (North America first) once,
> convert it into a MEDS forcing store, and cut a single site out of that store in seconds.
>
> Extends `MEDS_FORCING_DESIGN.md` §7, which stays the reference for the forcing NetCDF format
> (§7.1) and the ERA5-Land de-accumulation recipe (§7.3). Nothing here changes the Fortran reader.
> Reader improvements that would help are listed separately, in §15.
>
> Baseline: `main` = `beta` at `92cad44` (v0.2.2). The `file:line` citations are to that commit.

---

## 0. Summary

**What exists today.** There are two scripts, both written for one grid cell or a small box:

- `scripts/download_era5land.py` fetches a small box around one point from the Copernicus Climate
  Data Store (CDS), one month per request.
- `scripts/prep_era5land_forcing.py` de-accumulates the fluxes, converts units, and writes the
  MEDS `(time, grid)` forcing NetCDF.

Neither scales past a handful of cells. Two further facts limit what they can be used for:

- **The reader handles one file, one grid cell per run.** It reads every forcing value as a single
  1×1 read at `[record, grid]`.
- **Performance is set by the file's chunking**, which the current writer never controls (§2.3).

**What this plan builds:**

1. **`scripts/prepare_era5/`**, a single home for all ERA5 preparation. It holds a small shared
   library, so de-accumulation, humidity, the writer and the reader-contract checks exist exactly
   once.
2. **A clearly labelled single-grid workflow.** These are the two existing scripts, renamed
   `*_single_grid.py`, with the file-layout defects fixed. They gain a new backend for the CDS
   *point time-series* dataset, which returns decades of one cell in a single cheap request.
3. **A continental workflow**, in three scripts:
   - **download**: GRIB tiles, one month per request, resumable across weeks of CDS queueing;
   - **build**: a land-only, time-series-chunked **MEDS forcing store**, one standard MEDS forcing
     file per year plus a static grid file;
   - **extract**: nearest land cell to a small single-site file, with a printed config snippet.

   The store is organised by named domain, so other continents are one registry entry away.

**Why extraction, not direct reads.** A single-site run against the continental data should use
extraction. The reader cannot span multiple files (§2.2), and a multi-decade continental single
file would be terabytes. Each yearly store file is still a valid forcing file in its own right,
so runs of up to one year can read it directly (§8.4).

## 1. Goals and non-goals

**Goals** (as requested, 2026-09-25):

- **G1.** Create the `scripts/prepare_era5/` subfolder.
- **G2.** Move the current single-grid scripts into it, and make their single-grid purpose
  unmistakable in the name, the docstring, the `--help` text and the README.
- **G3.** Add scripts that download every relevant forcing field for all of North America, and
  that process it into a MEDS-compatible form a single-site run can use efficiently. The design
  must extend to other continents without code changes.

**Non-goals for this plan:**

- **No change to the Fortran reader** (`src/forcing/`). Improvements it would benefit from are
  proposed in §15 for `docs/ROADMAP.md`.
- **No gap-filling.** MEDS never gap-fills (`MEDS_FORCING_DESIGN.md` §5.5), and neither do these
  scripts. A missing value is a hard error, and cells without data are excluded from the grid,
  never filled.
- **No transient CO₂.** That is ROADMAP #184. The stores omit `CO2air`, so the reader uses
  `forcing.co2_const`.
- **No multi-polygon runtime** (ROADMAP #183). A multi-grid extract (§8.3) is the input that
  runtime will eventually want.

## 2. Where things stand

### 2.1 Inventory

| File | Does | Limits |
|---|---|---|
| `scripts/download_era5land.py` | Fetches a CDS `reanalysis-era5-land` NetCDF over a `±half-box` around one point, one request per month. Handles the mixed-stepType zip fallback. Pads one trailing day for de-accumulation. | Imports `xarray`, which is not in `environment.yml`, so the script does not run in the `meds` env. NetCDF requests cost double (§6.1). No resume. |
| `scripts/prep_era5land_forcing.py` | De-accumulates `tp`/`ssrd`/`strd` (the 00Z trap handled). Converts units. Computes `Qair` from dewpoint (Bolton, matching `meds_thermo`) and wind magnitude. Writes a `(time, grid)` file. | Loads every variable fully into memory. Loops per cell and per hour in Python. Raises on any NaN, so `--all-cells` over a box containing ocean fails. Writes one constant `--elevation` to every cell. No chunking control (§2.3). |

These two paths are cited from about ten places (§12), all of which must move with the scripts.

### 2.2 The reader contract (what any file we write must satisfy)

Established by reading `src/forcing/meds_met_driver.f90` (cited as MD below) and its config leaf:

| Contract point | Consequence for the writers | Where |
|---|---|---|
| Dimensions must be named exactly `grid` and `time`. Variables are `(time, grid)`, time slowest. The order is **not checked**. | Always write `(time, grid)`. A transposed file is misread silently. | MD:101, 113, 638–639 |
| Always required: `Tair`, `Qair`, `PSurf`, `Wind`, `Rainf`. Also `LWdown` (unless it is synthesized) and `SWdown` (the total) when `sw_partition` is not `passthrough`. `CO2air` is optional and falls back to `co2_const`. | Write these seven. Omit `CO2air`. | MD:556–625, 646–659 |
| `time:units` is parsed after the word `since`. **The word before `since` is ignored, so values are always read as seconds.** | Write `"seconds since YYYY-MM-DD HH:MM:SS"`, with no `Z`, offset or fraction. The xarray default `hours since` would be silently wrong. | MD:674–682; `meds_time.f90:325–358` |
| Record spacing is checked against `forcing.timestep` (±0.5 s), and the axis must be uniform. | The time axis must be exactly hourly and gap-free. | MD:714–729 |
| Values are read with 1×1 point reads per variable per record. Both bracket records are re-read whenever the bracket moves. | **Disk cost is governed entirely by the HDF5 chunk shape.** | MD:515–522, 628–643 |
| No multi-file support: one `forcing.path` of at most 256 characters, and one `nc_open`. Past the last record, the reader silently clamps. | A multi-year run needs one file spanning the period. | `meds_forcing_config.f90:78`; MD:97, 279–282 |
| `_FillValue` / `missing_value` are **not read**. Only NaN is detected. `scale_factor` / `add_offset` are **not applied**. | Declare no `_FillValue`, and never pack to int16. Store float32. | MD:577–583, 662–671 |
| `latitude(grid)` / `longitude(grid)` are read only for `grid_match = "nearest"`. The distance is great-circle, and a tie keeps the lowest index. `elevation(grid)` is never read; the lapse correction uses `[site].grid_elevation`. | The extractor must mirror the nearest-match rule, and must *print* the matched cell's elevation for `grid_elevation`. | MD:536–553, 696–698; `meds_forcing_kernels.f90:382–391` |
| Global attributes checked when present: `avg_convention` must match the config, and `sw_input_kind` must be `"total"` or `"components"`. `MEDS_FORCING_DESIGN.md` §7.1 says `"fourstream"`, which the reader does not recognise. | Write `avg_convention = "end"` and `sw_input_kind = "total"`. | MD:733–757 |
| `recycle_start` must land exactly on a record stamp. ERA5-Land "end" files are stamped at 01:00. | Keep the 01:00 first-stamp convention (§7.3). | MD:199–213; `meds_config_main.toml:413–435` |

### 2.3 A defect that already costs single-site users

The current writer creates `time` as UNLIMITED and sets no chunk sizes. With netCDF4-python 1.7.4
and libnetcdf 4.10.1, that yields **`(1, 1)` chunks**. Measured on a two-year, one-cell,
seven-variable file:

| Layout | Size |
|---|---|
| unlimited `time` (current) | **6.27 MB** |
| fixed `time`, contiguous | **0.50 MB** |

That is 12.5× larger, and all of the difference is HDF5 chunk-index overhead. For a continental
file, the same default would give `(1, ngrid)` chunks: every 1×1 read would pull about 1.2 MB, and
a ten-year single-site run would move hundreds of gigabytes. The fix is part of P1 (§5.3).

### 2.4 What the CDS allows (researched 2026-09-25)

**Cost and limits** (checked against the live `/costing` endpoint):

- A request may cost at most **12,000**. Cost = variables × days × hours, and **netCDF counts
  double**.
  - 8 variables × 1 month costs 5,952 in GRIB and 11,904 in NetCDF.
  - 1 variable × 12 months costs 8,928 in GRIB (accepted) and 17,856 in NetCDF (rejected).
- **Area subsetting does not reduce the cost**, so a continental request costs the same as a
  global one.
- In practice about **one running request per user**, with the rest queued. Expect waits of hours
  to days.
- ECMWF recommends GRIB for volume. Its NetCDF output is labelled "Experimental".

**The point time-series dataset** (`reanalysis-era5-land-timeseries`):

- Contains all eight variables, with `tp`/`ssrd`/`strd` **already de-accumulated** to hourly
  amounts.
- Runs from 1950-01-02 to about five days before today.
- Returns the nearest 0.1° cell.
- A full-record, eight-variable NetCDF request costs **8 out of 500**.
- It is a beta service.

**ECMWF's analysis-ready Zarr copies ("ARCO")** of ERA5-Land:

- A time-series-chunked copy exists, with chunks of (33792 time, 4 lat, 8 lon), already
  de-accumulated.
- Access uses the CDS API key as a Bearer token.
- It is beta, with undisclosed rate limits, and ECMWF steers heavy bulk processing to traditional
  downloads.
- Earth Data Hub's copy is *not* de-accumulated, is quota-limited, and has balanced chunks. It is
  not recommended here.

**Grid and fields:**

- ERA5-Land is 1801×3600 at 0.1°, latitude 90→−90, longitude 0→359.9. Normalise longitudes in code.
- `land_sea_mask` and `geopotential` are available as invariant fields, the latter giving
  orography.

**ERA5 (not Land) single levels** accumulate over the preceding hour with **no** daily reset. A
future ERA5 adapter therefore needs a different de-accumulation rule (§9.2).

Sources are in Appendix B.

## 3. Target layout

```
scripts/prepare_era5/
├── README.md                             # which workflow to use; quick starts; data-root layout
├── environment.yml                       # meds-era5 env: cdsapi, eccodes, netcdf4, numpy, pytest
├── domains.toml                          # named domains (§9.1)
├── era5lib/                              # shared library — importable because the script dir is on sys.path
│   ├── __init__.py
│   ├── variables.py                      # variable catalogue per source (§4)
│   ├── physics.py                        # de-accumulation, humidity, wind, unit conversion
│   ├── timeaxis.py                       # interval-year stamps, month tiles, leap years
│   ├── grid.py                           # longitude normalisation, 0.1° alignment, great-circle nearest
│   ├── meds_writer.py                    # THE writer: the reader contract of §2.2, chunking, quantization
│   ├── contract.py                       # checks a file against §2.2 (also a CLI, below)
│   ├── cds.py                            # request building, retries, zip handling, costing, manifest
│   └── domains.py                        # loads domains.toml; antimeridian split
│
├── download_era5land_single_grid.py      # was scripts/download_era5land.py        (§5)
├── prep_era5land_single_grid.py          # was scripts/prep_era5land_forcing.py    (§5)
│
├── download_era5land_continental.py      # §6
├── build_era5land_continental_store.py   # §7
├── extract_era5land_site.py              # §8
│
├── verify_meds_forcing.py                # any forcing file vs the reader contract (§11.1)
└── tests/                                # pytest, offline, synthetic data (§11)
```

**Naming.** Every script name says `single_grid` or `continental`, so a listing alone tells you
which workflow a script belongs to. The prefix is `era5land_`, not `era5_`, because the data set is
ERA5-Land. `era5lib` stays source-generic so an ERA5 single-levels adapter (§9.2) can drop in.

**Layout choice.** The directory is flat, with one shared package, rather than `single_grid/` and
`continental/` subfolders. Subfolders would need `sys.path` manipulation or `python -m` invocations
to reach the shared package, because `scripts/` is not a package.

## 4. The shared library `era5lib/`

Each piece of logic lives here exactly once. Every script is a thin command-line layer over it.

- **`variables.py`** — one catalogue row per MEDS input. Each row gives the CDS request name, the
  GRIB `shortName`/`paramId`, the NetCDF name, the step type (`instant` or `accum`), and the
  **accumulation semantics of each source**:

  | Source | Accumulation semantics |
  |---|---|
  | ERA5-Land gridded | `since_00utc` |
  | ERA5-Land time series | `hourly_amount` (already de-accumulated) |
  | ERA5 single levels | `hourly_amount` (§9.2) |

  Downstream code never branches on a source name; it branches on these semantics.
- **`physics.py`** — the current prep script's formulas, vectorized over a `(time, cells)` block:
  - **`deaccumulate_since_00utc(raw, valid_hour, lookahead_row)`.** At 01 UTC the value is taken
    as-is; otherwise it is `raw(H) - raw(H-1)`. The 00Z stamp therefore yields the previous day's
    `[23Z, 00Z]` hour. Negative packing noise is clipped for `tp` and `ssrd` and left alone for
    `strd`, all as today. The `lookahead_row` argument carries the first hour of the *next* tile,
    so that month and year boundaries are exact (§7.3).
  - **`dewpoint_to_specific_humidity`.** Bolton 1980, **identical to `meds_thermo`**; a test pins
    the constants (§11.1).
  - **Wind.** Wind magnitude with the `U_MIN = 0.1 m/s` floor.
  - **Units.** `tp` goes from m/h to kg m⁻² s⁻¹ (× 1000 / 3600). `ssrd` and `strd` go from J m⁻²
    per hour to W m⁻² (÷ 3600).
- **`timeaxis.py`** — end-stamped hourly axes, the interval-year convention (§7.3), and month-tile
  enumeration, including the trailing tile.
- **`grid.py`** — the native ERA5-Land grid geometry (0.1° spacing, row and column indices) and
  longitude normalisation to −180..180. It aligns box edges to 0.1° multiples, because unaligned
  boxes come back on a shifted grid. Its great-circle nearest search is a line-for-line mirror of
  `nearest_grid_index`, including the lowest-index tie rule.
- **`meds_writer.py`** — `write_meds_forcing(path, time_stamps, latitude, longitude, elevation,
  fields, attrs, *, chunk_grid=None, quantize=None, compress=None)`. It enforces everything in §2.2
  by construction:
  - a **fixed** `time` dimension and `"seconds since …"` units;
  - float32, with **no `_FillValue`**;
  - a NaN or out-of-range check before anything is written (§11.2);
  - `avg_convention = "end"` and `sw_input_kind = "total"`, plus the provenance attributes;
  - a contiguous layout when `chunk_grid` is None (single-site files), otherwise `(ntime,
    chunk_grid)` chunks;
  - optional zlib with shuffle, and optional netCDF quantization (§7.2).
- **`contract.py`** — `check_meds_forcing(path) -> list[Finding]`, covering:
  - dimension names and order;
  - required variables and the `seconds since` units;
  - uniform spacing;
  - no NaN and no `|x| > 1e15` (which is what a stray `_FillValue` looks like to the reader);
  - the global-attribute vocabulary;
  - a chunk-shape warning when chunks are `(1, *)` or larger than about 4 MB.
- **`cds.py`** — builds requests and retrieves with retry and exponential backoff. It detects and
  merges zips, POSTs to the `/costing` endpoint for `--dry-run`, and keeps an atomic JSON manifest
  (§6.3).
- **`domains.py`** — parses `domains.toml` and splits any box that crosses the antimeridian.

## 5. Single-grid workflow (G2)

### 5.1 Move and rename

- Use `git mv`, so the history follows the files:
  - `scripts/download_era5land.py` becomes `scripts/prepare_era5/download_era5land_single_grid.py`;
  - `scripts/prep_era5land_forcing.py` becomes `scripts/prepare_era5/prep_era5land_single_grid.py`.
- **Existing flags stay**, so current invocations only change path: `--lat`, `--lon`, `--start`,
  `--end`, `--half-box`, `--split-requests`, `--cells` and `--all-cells`.
- **No shim** is left at the old paths. All references move in the same PR (§12). This is decision
  D5 in §14.

### 5.2 Make the single-grid scope unmistakable

- **Docstring.** The first line becomes *"SINGLE-GRID workflow: forcing for one site (or a small
  box of cells) …"*. A second line reads: *"For many sites, or regional or continental coverage,
  use the continental workflow: `download_era5land_continental.py` →
  `build_era5land_continental_store.py` → `extract_era5land_site.py`."*
- **Help text.** `--help` gets the same epilogue.
- **Guard.** A box larger than about 2°×2° prints a warning pointing to the continental workflow.
  The cost does not grow with area (§2.4), but the prep script's memory and per-cell loop do.

### 5.3 Fix the file layout (behaviour-preserving for the values)

The prep script is reimplemented on `era5lib`. Changes to its output:

- **Time dimension.** `time` becomes fixed and contiguous, removing the 12.5× overhead of §2.3.
- **Fill value.** No `_FillValue` is declared. Today's `fill_value = 1e20` is never *used*, because
  a NaN aborts first, but a partially written file would present `1e20` to the reader as data.
- **CO₂.** `CO2air` is **omitted by default**, so the reader uses `forcing.co2_const`, which is a
  required key anyway. `--co2 VALUE` writes the constant field when wanted. This change goes in the
  CHANGELOG.
- **Elevation.** `elevation(grid)` comes from ERA5-Land orography (geopotential ÷ 9.80665) when a
  static file is supplied (`--static`). The constant `--elevation` stays as a fallback. The script
  prints a config snippet with `grid_elevation`, because the reader never reads `elevation(grid)`.
- **Dependencies.** The downloader drops `xarray`: the zip-member merge and month concatenation use
  `netCDF4`, which is a small change. The whole workflow then runs in the `meds-era5` env with no
  extra packages.

**Acceptance:** on the same raw input, the new prep script's data variables are **bit-identical**
to the old script's (ignoring `CO2air`), and `verify_meds_forcing.py` passes.

### 5.4 New backend: the CDS point time series (`--source timeseries`)

`download_era5land_single_grid.py --source {box,timeseries}`:

- **`box`**: today's behaviour, using the gridded dataset over a small area.
- **`timeseries`**: requests `reanalysis-era5-land-timeseries` with
  `location = {latitude, longitude}`, a `date` range, and NetCDF output. It needs one request for
  the whole period: decades cost 8 out of 500, against 12 or more queued requests per year for
  `box`. It fetches the nearest 0.1° cell only.
- The raw file carries a global attribute `era5prep_accumulation = "hourly_amount"`. The prep
  script reads it and **skips de-accumulation** for those files. It still converts units: m and
  J m⁻² per hour into rates.

**Open points P0 must pin down** (§13) before `timeseries` becomes the default:

1. The time-stamp convention: is the hourly amount stamped at the *end* of its hour, as in the
   gridded data?
2. The units.
3. The start at 1950-01-02, not 01-01.
4. How the point service picks a cell next to the coast.

**The pin is an equivalence test.** Ithaca, June–August 2020, via `box` and via `timeseries`, must
give identical MEDS files within float32 round-off. If it passes, `timeseries` becomes the default
(decision D7).

## 6. Continental download (G3, part 1)

`download_era5land_continental.py --domain north_america --start 1981-01 --end 2025-12 --root /data/era5land`

### 6.1 Request design

- **GRIB, all eight variables, one calendar month per request**, over the domain's box. This costs
  5,952 out of 12,000 and gives 12 requests per year.
- GRIB is chosen for three reasons:
  - NetCDF costs twice as much.
  - GRIB returns mixed step types in one file, with no zip.
  - ECMWF labels its NetCDF conversion "Experimental".
- Two months per request (11,904) would also fit. That is an option, `--months-per-request 2`, only
  once P0 has confirmed the API accepts it.
- **Trailing tile.** If the range ends in month M, the tool also fetches day 1 of month M+1 (a small
  request). It supplies the 00Z stamp that closes the last interval (§7.3). This is the
  continental form of the existing "pad one trailing day" rule.
- **Static tile, once per domain.** Invariant `land_sea_mask` and `geopotential` over the same box.
  As a fallback, cut them from the ECMWF-hosted 0.1° files listed in Appendix B.
- **Multi-box domains issue one request per box per month.** Area does not reduce cost, so each box
  costs a full request. Keep domains to one box where possible (§9.1).

### 6.2 Throughput reality

The CDS allows about one running request per user. Throughput is therefore set by queue plus run
time per request, *not* by bandwidth:

- **At 2–6 h per request,** one data-year (12 requests) takes 1–3 days.
- **1981–2025 (45 years) is about 1.5–4.5 months** of wall time.

The tool is built for that:

- **Keep a submission window full.** It keeps up to `--max-queued` requests submitted (default 3),
  so the next request is always waiting, while only about one runs. P0 must confirm the submit and
  poll calls of the installed `cdsapi` / ECMWF data-stores client.
- **Survive interruption.** A long-lived process (tmux, `nohup`, or `sbatch` with requeue) resumes
  from the manifest after any interruption.
- **Order the work usefully.** The default order is newest year first, so recent years are usable
  early. `--order oldest` is available.
- **Offer a faster route to test.** The ARCO Zarr route (§6.5) exists largely because of this
  wall-time figure.

### 6.3 Manifest, resume, integrity

- **The manifest.** One JSON manifest per domain, at `<root>/<domain>/raw/manifest.json`. It holds
  one record per tile: key, request body, state, CDS request id, submit and finish times, bytes,
  sha256, and the verification result. The states are `planned`, `submitted`, `downloaded`,
  `verified` and `failed`.
- **Atomic writes.** Downloads go to `*.part`, and a rename is the commit. The manifest is written
  to a temporary file and renamed.
- **One instance at a time.** A lockfile stops two instances sharing a manifest.
- **Verification.** A tile is verified when all of the following hold:
  - the GRIB message count equals 8 × hours in the month;
  - every validity time is present exactly once;
  - the grid shape and corners match the domain;
  - the land points are non-missing.

  A tile that fails is re-queued; after `--max-attempts` it is marked `failed` and reported.
- **Reporting commands:**
  - `--status` prints the counts in each state and the estimated remaining wall time;
  - `--dry-run` prints the request plan and each request's cost from `/costing`, with no login
    needed.

### 6.4 Raw data layout (outside the repository)

```
<root>/era5land/<domain>/
  raw/manifest.json
  raw/static/era5land_<domain>_static.grib           # lsm + z
  raw/grib/<YYYY>/era5land_<domain>_<YYYYMM>.grib    # 8 vars × all hours of the month
  raw/grib/<YYYY+1>/era5land_<domain>_<YYYY+1>0101_trail.grib
  work/…                                             # §7.4 intermediates, deletable
  store/…                                            # §7.1 the product
```

### 6.5 Alternative source to evaluate in P0: ECMWF ARCO Zarr

**What it offers.** The time-series-chunked Zarr copy, with chunks of (33792 h ≈ 3.85 years × 4
lat × 8 lon), is already arranged the way the store needs it, and already de-accumulated. Reading
it would skip the CDS queue *and* the reordering pass of §7.4.

**What counts against it:**

- it is beta and rate-limited;
- ECMWF steers heavy bulk use to the CDS;
- partial periods waste whole 3.85-year chunks;
- it adds `zarr` and `fsspec` dependencies.

**The spike.** Pull one year for a 5°×5° North American box. Measure throughput, error rate and
rate limiting, and confirm the stamp convention against the gridded path. Decide in D4.

**If adopted, the scope is narrow.** It becomes a `--source arco` path in the downloader that writes
the same intermediate files as §7.4 pass A. Everything downstream is unchanged.

## 7. The continental forcing store (G3, part 2)

### 7.1 Store format

```
<root>/era5land/<domain>/store/
  era5land_<domain>_static.nc        # the grid definition (identical for every year)
  era5land_<domain>_<YYYY>.nc        # one standard MEDS forcing file per calendar year
  store_manifest.json                # version, domain, years, n_grid, checksums, script version
```

**Grid.**

- **Land cells only**, as an unstructured `grid` dimension. A cell is in the grid when ERA5-Land
  supplies data for it: non-missing in the GRIB bitmap, with an lsm threshold as an optional
  filter (decision D2).
- The grid is **fixed once, from the static tile and the first processed year**, and is identical
  in every year file. A given `grid_index` therefore means the same cell across the whole store.
- If a later year has a missing value in an included cell, the builder **hard-errors** with a
  report, following the no-gap-fill rule. It never fills and never drops the cell silently.
- Order is row-major, north to south and then west to east, which is deterministic and needs no
  extra dependency. Space-filling-curve ordering would improve locality for regional multi-site
  extracts; that is an option for later (§14).

**`era5land_<domain>_static.nc`** holds:

- `latitude(grid)`, `longitude(grid)` (−180..180) and `elevation(grid)` in metres (geopotential ÷
  9.80665);
- `land_fraction(grid)`, `native_row(grid)` and `native_col(grid)`;
- the native axes `lat(nlat)` and `lon(nlon)`;
- `grid_map(nlat, nlon)` as int32: the **1-based** `grid_index` of each native cell, or 0 where the
  cell is not in the store. This makes the "which grid index is at this location?" lookup O(1);
- provenance attributes.

**Year files `era5land_<domain>_<YYYY>.nc`**, written by `meds_writer`:

- dims `(time, grid)`, with `time` fixed;
- `Tair`, `Qair`, `PSurf`, `Wind`, `Rainf`, `SWdown` and `LWdown` as float32;
- `latitude`, `longitude` and `elevation(grid)` duplicated from the static file (about 8 MB), so
  **every year file is a valid standalone MEDS forcing file** that nearest-match works on;
- no `CO2air` and no `_FillValue`;
- the global attributes of §2.2, plus `store_domain`, `store_year`, `store_version`, and the
  static file's sha256.

### 7.2 Chunking, compression, precision

- **Chunks: `(ntime_year, G)`** — one whole year by `G` cells. The default is **`G = 32`**:
  8784 × 32 × 4 B ≈ 1.1 MB uncompressed. P0 settles the final value (§13), by benchmarking two
  access patterns:
  - **extraction**, where one chunk per variable per year is fastest, and larger `G` is tolerable;
  - **the reader's direct 1×1 reads**, where the chunk must fit the libnetcdf 4.10.1 per-variable
    chunk cache, or every point read re-inflates the chunk.
- **Compression.** zlib level 1–4 with shuffle.
- **Precision: netCDF quantization** (`significant_digits`, BitGroom or GranularBR). This needs
  libnetcdf ≥ 4.9 and netCDF4-python ≥ 1.6; the installed 4.10.1 and 1.7.4 qualify.
  - Values stay float32, so it is **transparent to the reader**, which cannot unpack int16 (§2.2).
  - Unlike int16 packing, it is relative precision, so light drizzle keeps its digits.
  - Proposed digits: `Tair` 5, `PSurf` 6, and 4 each for `Qair`, `Wind`, `Rainf`, `SWdown` and
    `LWdown`. P0 measures the error, and "off" is always available.

### 7.3 The time convention: interval years

Year file `Y` holds the **8760 or 8784 end-stamped hourly intervals that make up calendar year Y**:

- stamps run from **01:00 Jan 1 Y** to **00:00 Jan 1 Y+1**;
- units are `seconds since Y-01-01 00:00:00`, so values run 3600 … 3600·n.

For the flux variables, the stamp at hour H means the mean over [H−1, H]. This matches
`avg_convention = "end"` and the existing single-grid files' 01:00 first stamp, which the recycle
anchoring expects.

**De-accumulation needs exactly one hour of lookahead, never lookbehind:**

- **01 UTC** takes the raw value as-is.
- **00 UTC on day d+1** is `raw(00Z d+1) − raw(23Z d)`.
- **The last interval of month M** therefore needs the 00Z stamp of day 1 of month M+1, taken from
  the next tile or from the trailing tile (§6.1).
- **Year files are independent** given that one lookahead hour, so the builder parallelises across
  years.

### 7.4 `build_era5land_continental_store.py` — algorithm

```
build_era5land_continental_store.py --domain north_america --root /data/era5land \
    --years 1981-2025 [--workers 8] [--band-cells 16384] [--chunk-grid 32] [--quantize default]
```

Raw tiles are map-major (all cells at one time step), while the store is time-series-major. That
reordering is the whole cost, so it runs in two passes with bounded memory:

**Pass 0 — the grid, once.**

1. Decode the static tile and the first year's January.
2. Fix the land-cell list and write `era5land_<domain>_static.nc`.

**Pass A — land compaction, per month tile, parallel over months.**

1. Decode the GRIB with `eccodes`, streaming message by message, and order the messages by
   `validityDate` / `validityTime`.
2. Gather the land cells through the static index.
3. Write `work/land/<YYYY>/<YYYYMM>.nc` holding the **raw** values: `(hours, n_grid)` for the eight
   raw variables, chunked `(hours, 4096)`, with zlib.

About 8 GB per month uncompressed; deletable after pass B.

**Pass B — assembly, per year × grid band, parallel over years.**

For each band of `--band-cells` cells:

1. Read the band from the 12 monthly land files, plus the first hour of the following January (the
   lookahead).
2. Compute, vectorized over `(time, band)`:
   - `Qair` from `d2m` and `sp`;
   - `Wind` from `u10` and `v10`;
   - the de-accumulation of `tp`, `ssrd` and `strd` by validity hour;
   - the unit conversions.
3. Range-check the results (§11.2).
4. Write `[:, band]` into the year file.

With a 16k-cell band, peak memory is about 8785 × 16384 × 4 B × 10 arrays ≈ 5.8 GB per worker.

**Finish:**

1. Run `verify_meds_forcing.py` on each year file.
2. Write `store_manifest.json` with the checksums.
3. With `--delete-work`, remove the intermediates.

**Idempotence.** Each output is written to `*.part` and renamed. A rerun skips year files whose
manifest checksum matches, and `--force-year` rebuilds one year. Appending a new year never touches
the old ones.

**Implementation note.** The implementation is pure `numpy` + `netCDF4` + `eccodes`, deliberately
without `dask`/`xarray`/`rechunker`: the reordering is simple enough to write by hand, and every
extra package is a burden on HPC installs. Revisit only if P0 shows pass B is I/O-bound in a way a
library would fix.

## 8. Site extraction — efficient single-site runs (G3, part 3)

### 8.1 `extract_era5land_site.py`

```
extract_era5land_site.py --store /data/era5land/north_america/store \
    --lat 42.44 --lon -76.50 --start 2000-01-01 --end 2020-12-31 --out ithaca_2000_2020.nc \
    [--max-distance-km 15] [--name ithaca]
```

1. **Load the static file.** Go from the target's native row and column to `grid_map`, and search a
   ±k-cell window for in-store cells.
2. **Pick the nearest cell** by the **same great-circle rule and tie-break as the reader**
   (`era5lib.grid`, which mirrors `nearest_grid_index`). Report the distance. Hard-error beyond
   `--max-distance-km`, for example a coastal site whose covering cell ERA5-Land does not supply;
   this answers `MEDS_FORCING_DESIGN.md` §10 Q7 for the scripts.
3. **Read the data.** For each year touched, read `[time slice, g]` for the seven variables. That
   is one chunk per variable per year.
4. **Assemble.** Concatenate the years and trim to the interval `[start 01:00, end+1 00:00]`.
   Refuse a range the store does not cover; never clamp.
5. **Write** a **single-site MEDS file** (`grid = 1`, fixed `time`, contiguous) through
   `meds_writer`, carrying provenance (store version, domain, source cell, distance), then check it
   with `contract.py`.
6. **Print a ready-to-paste config fragment:**

   ```toml
   [forcing]
   path = "ithaca_2000_2020.nc"
   grid_match = "explicit"
   grid_index = 1
   [site]
   latitude = 42.44
   longitude = -76.50
   grid_elevation = 382.7      # illustrative value: the ERA5-Land orography of the matched cell
   ```

**Performance target:** at most 30 s for 45 years of one site on local disk. P0 sets the real
number. The first-order cost is 45 years × 7 chunks, about 1.1 MB each uncompressed.

### 8.2 Why extraction is the primary path

- **It needs no reader change.**
- **It produces the same kind of tiny file** today's single-site users already run on, so their
  existing configs keep working.
- **It isolates runs from the store.** Runs stay reproducible even if the store is rebuilt, because
  the extract records its provenance, and a run never holds a multi-terabyte store open.

### 8.3 Multi-site extraction

`--sites sites.csv` takes a CSV with columns `name,lat,lon[,start,end]`. `--layout` chooses the
output:

- **`per-site`** (default) writes one file per site.
- **`multi-grid`** writes one `(time, N)` file. Each run then selects `grid_index`, and the future
  multi-polygon runtime (#183) gets its input.

Sites are grouped by grid chunk, so each chunk is read once.

### 8.4 Direct reads of a year file (secondary)

A year file is a valid forcing file, so a run of up to one year can point `forcing.path` at it with
`grid_match = "nearest"`. The reader then scans all `n_grid` latitudes and longitudes (about 5 MB)
and performs its 1×1 reads. Whether that is fast depends on `G` against the chunk cache, which is
the P0 benchmark in §7.2. Multi-year direct reads need reader multi-file support (§15 F1).

## 9. Extending to other continents and sources

### 9.1 Domain registry (`domains.toml`)

```toml
[north_america]
description = "Canada, USA incl. Alaska, Mexico, Central America, Caribbean"
boxes = [[84.0, -170.0, 5.0, -50.0]]      # [N, W, S, E], edges on 0.1° multiples

[north_america_greenland]                  # optional add-on (Greenland reaches ~-11°E)
boxes = [[84.0, -75.0, 59.0, -10.0]]

[aleutians_west]                           # west of the antimeridian; W > E ⇒ split in two requests
boxes = [[56.0, 172.0, 51.0, -170.0]]

[south_america]
boxes = [[13.0, -82.0, -56.0, -34.0]]
[europe]
boxes = [[72.0, -25.0, 34.0, 45.0]]
[africa]
boxes = [[38.0, -18.0, -35.0, 52.0]]
# asia, oceania … and "global" (no area key) — defined, exercised by tests for parsing only
```

- **Nothing else knows about continents.** The downloader, builder and extractor read the domain
  name, and a new continent is one registry entry.
- **Cost.** Every box is a full-cost request each month (§2.4). If more than two continents are
  planned, it is cheaper in requests to download **global** tiles once and build each continental
  store from the same raw tiles. The builder takes `--raw-domain global --domain europe` for this.
  The trade is about 3–4× the North American transfer volume per tile for zero extra requests.
- **Antimeridian.** A box with W > E is split into two requests, and the builder stitches the
  halves into one longitude-normalised grid.

### 9.2 Other sources

`variables.py` and the `accumulation` semantic keep this additive:

- **ERA5 single levels** (`reanalysis-era5-single-levels`, 0.25°, 1940–):
  - hourly accumulations with no daily reset, so the `hourly_amount` path with no
    de-accumulation;
  - ocean coverage, so no land-only mask is needed for coastal sites;
  - different variable names in a few places.

  This is useful for sites ERA5-Land misses, or for years before 1950. It is a P6 option, not part
  of the North America deliverable.
- **Other reanalyses** would each get a catalogue row set and writer reuse. They are out of scope
  here.

## 10. Environment and dependencies

- **A separate environment.** `scripts/prepare_era5/environment.yml` defines env `meds-era5`:
  `python >=3.10`, `numpy`, `netcdf4`, `cdsapi >=0.7.7`, `python-eccodes`, and `pytest`. The
  `zarr` and `fsspec` packages are added only if D4 adopts ARCO.
  - It is separate because the main `environment.yml` serves the Fortran build and post-processing.
    Download tooling should not weigh on it.
  - `tomllib` (stdlib, Python ≥ 3.11) reads `domains.toml`, with a `tomli` fallback for 3.10.
- **Credentials.** `~/.cdsapirc` holds exactly two lines: `url: https://cds.climate.copernicus.eu/api`
  and `key: <token>`. Both ERA5-Land licences (gridded and time series) must be accepted once on the
  CDS website. The README says so first.
- **HPC.** The README gives an `sbatch` template for pass A and pass B (array job over months or
  years), and a tmux or `nohup` recipe for the long-lived downloader.

## 11. Validation and testing

### 11.1 Offline unit tests (`scripts/prepare_era5/tests/`, pytest, synthetic data)

- **`physics`:**
  - the 00Z trap and the 01Z rule;
  - month-boundary and year-boundary lookahead;
  - leap day, and the first sample of a record;
  - the clip policy for each variable.
- **Humidity parity with `meds_thermo`.** The constants are asserted against the Fortran source, and
  a table of values is checked.
- **`timeaxis`:** 8760 or 8784 stamps per interval year, with the first at 01:00 and the last at
  00:00 of Jan 1 Y+1.
- **`meds_writer` + `contract`:**
  - every writer output passes;
  - hand-made bad files each trigger their finding: `hours since`, a transposed variable, an
    unlimited `(1, 1)` layout, a declared `_FillValue`, a `1e20` value, NaN, and non-uniform
    spacing.
- **`grid`:** nearest-match parity with a Python port of `nearest_grid_index`, including ties;
  longitude normalisation; 0.1° alignment.
- **`domains`:** parsing, and the antimeridian split.
- **`cds`:** manifest state transitions, resume after a crash at each state, zip handling, and
  retry and backoff, all against a mock client with no network.
- **Store round trip.** Synthetic map-major tiles go through pass A, pass B and extraction, and the
  result must equal a direct computation on the synthetic truth, within the quantization tolerance.

**Where the tests run** is decision D6: standalone `pytest`, or registered as a CTest target when
Python and pytest are found.

### 11.2 Plausibility gates (the builder and both writers)

MEDS's lesson is that "a closed budget proves bookkeeping, not plausibility", so the writers check
values before writing them. Hard bounds (error), for example:

| Variable | Bound |
|---|---|
| `Tair` | 170–340 K |
| `Qair` | 0–0.05 |
| `PSurf` | 30–110 kPa |
| `Wind` | 0.1–75 m s⁻¹ |
| `Rainf` | 0–0.2 kg m⁻² s⁻¹ |
| `SWdown` | 0–1500 W m⁻² |
| `LWdown` | 30–650 W m⁻² |

Soft checks (warning): `SWdown > 0` at night by solar zenith; daily `SWdown` above top-of-atmosphere.
Each year file also gets a store summary with per-variable global means and extremes, and each is
compared to the previous year's to catch unit or stamp slips.

### 11.3 Integration checks (manual, need CDS data; recorded in this plan when run)

1. **Three-way equivalence.** Ithaca, June–August 2020, produced three ways: single-grid `box`,
   single-grid `timeseries`, and continental extraction from a North American store built for
   2020. The three must agree within float32 round-off, or the quantization tolerance for the
   store.
2. **Model check.** `examples/example_biophysics` run on the extracted file must reproduce the run
   on the single-grid file: identical diagnostics when the inputs are identical, otherwise within
   tolerance.
3. **Recorded benchmarks** (P0 and P5): extraction seconds per site-year, and direct-read seconds
   per model-year, for each candidate `G`.

## 12. Documentation and reference updates

**Moved in P1, each changed from `scripts/download_era5land.py` or `scripts/prep_era5land_forcing.py`
to the new path:**

| File | Lines |
|---|---|
| `meds_config_main.toml` | 385–386 |
| `src/forcing/README.md` | 69–70 |
| `docs/science/forcing.md` | 58–60 (paragraph), 396 (table row) |
| `docs/ed2_comparison.md` | 319 |
| `examples/example_biophysics/README.md` | 243–244 |
| `examples/example_biophysics/run_example.py` | 181–182 |
| `MEDS_FORCING_DESIGN.md` | 58–59, 318–319, 646, 727, 735, 1001–1006, 1043, 1065, 1142 |

- **`MEDS_FORCING_DESIGN.md`** gets its paths updated and a one-line pointer to this plan at the
  top of §7. **Its section numbers stay as they are.**
- **Archived documents stay as they are**, since archived line references are pre-reorganization by
  policy.

**New documentation:**

- **`scripts/prepare_era5/README.md`**: a decision table (one site versus many sites or a region),
  quick starts for both workflows, the data-root layout, the credentials and licences, and the HPC
  recipes.
- **`docs/science/forcing.md`** gets a short "Continental stores" paragraph pointing to the README.

**Every PR** carries a `CHANGELOG.md` entry. **The last PR** moves this plan to `archive/` with a
tombstone and files the §15 items in `docs/ROADMAP.md` (CLAUDE.md documentation rules 2–4).

## 13. Phasing — pull requests against `beta`

| Phase | PR content | Acceptance |
|---|---|---|
| **P0** measure | No production code. Throwaway notebooks or scripts under `$SCRATCH`. The results are recorded in Appendix A of this plan in one docs PR. **(a)** North American land-cell count from the lsm. **(b)** One North American month by CDS GRIB: queue plus run time, size, eccodes decode time. **(c)** `timeseries` stamp and units vs `box` (Ithaca, one month). **(d)** ARCO spike (§6.5). **(e)** Chunk benchmark over `G` ∈ {8, 16, 32, 64, 128}: extraction, and Fortran direct reads through a `test_met_driver`-style harness. **(f)** Quantization error for each variable. **(g)** The submit and poll calls of `cdsapi`. | Appendix A filled in, and decisions D1–D4 and D7 taken. |
| **P1** reorganize (G1, G2) | `scripts/prepare_era5/` with `era5lib/` (physics, timeaxis, meds_writer, contract, variables), the two renamed single-grid scripts on top of it, the §5.3 layout fixes, `verify_meds_forcing.py`, `environment.yml`, the README, the §11.1 unit tests for what exists, and the §12 reference moves. | Data **bit-identical** to the old script on the same raw input. The file is about 12× smaller. The contract check is clean. All references move and no old path remains. CHANGELOG. |
| **P2** single-grid time series | `--source timeseries`, the accumulation-semantics switch in prep, and the equivalence test (§11.3 #1, first two legs). | `box` and `timeseries` agree. The default switches if D7 says so. |
| **P3** continental download | `domains.toml` and `domains.py`, `cds.py` (manifest, costing, retries), `download_era5land_continental.py`, the static tile. | `--dry-run` cost plan is correct. Resume is proven by killing the process in each state (mock plus one real month). A tile verifies. |
| **P4** store builder | Passes 0, A and B, `store_manifest`, the plausibility gates, and idempotent reruns. It can be developed on synthetic tiles while the P3 download runs. | The synthetic round trip passes. The 2020 North American store builds, and each year file passes the contract check. The peak memory of the stated band size is measured. |
| **P5** extraction and end-to-end | `extract_era5land_site.py`, including multi-site and the config snippet, and §11.3 #1–#3. | The three-way equivalence holds. The MEDS example reproduces. The extraction performance target is met and recorded. |
| **P6** extend and close | Test one non-North-American domain on one month. Optionally, an ERA5 single-levels adapter and ARCO (per D4). File the §15 items in ROADMAP. Archive this plan. | Per-item. The plan is tombstoned. |

**Dependencies between phases.**

- P1 blocks everything else, because of the shared library.
- **Start the P3 real download as soon as P3 merges.** The download is the long pole at weeks to
  months (§6.2); P4 and P5 proceed on synthetic data and the first finished years in parallel.

## 14. Open decisions

| # | Decision | Recommendation |
|---|---|---|
| D1 | **Period** for the North American store: 1981–present, 1950–present, or a shorter target window. It drives download wall time and storage (Appendix A). | Start from the period your runs need. The downloader works newest year first, so a longer record can be added later without rebuilding. |
| D2 | **Domain extent and grid inclusion.** Include Greenland, the Caribbean and Central America, and the Aleutians west of 180°? Include coastal cells with small land fractions, or set an lsm threshold? | The `north_america` box as in §9.1, add-ons optional. Keep every cell ERA5-Land supplies, and record `land_fraction` so users can filter. |
| D3 | **Storage and compute location**: the data root (multi-TB) and the HPC scheduler, for the sbatch templates. | Your call. It is only a path and template in the README. |
| D4 | **Bulk source**: CDS GRIB (sanctioned, slow queue) or ARCO Zarr (no queue, beta, ECMWF discourages bulk use). | CDS GRIB as the baseline. Adopt ARCO only if the P0 spike shows it is reliable and the terms allow it. |
| D5 | **Shims** at the old script paths? | None. Move every reference in one PR. The old names are cited only in docs and one example. |
| D6 | **Test hosting**: standalone pytest, or a CTest target when Python and pytest are present. | A CTest target, so `ctest` covers the scripts too, skipped cleanly when pytest is absent. |
| D7 | **Default single-grid backend** after P2. | `timeseries`, if the equivalence holds: one cheap request per site and no queue. |

## 15. Proposed follow-ons outside `scripts/` (for `docs/ROADMAP.md`, not this plan)

These came out of reading the reader. None blocks this plan. Each becomes a ROADMAP entry with an
issue number when accepted.

- **F1 — multi-file forcing.** A file list or a `%Y` path template, with rollover at file edges.
  This would let multi-year runs read the yearly store directly, making extraction optional.
- **F2 — detect `_FillValue` / `missing_value`.** Today a `1e20` fill is used as data (MD:662–671).
- **F3 — reader I/O.**
  - Read both bracket records in one `count = [2, 1]` call, as `MEDS_FORCING_DESIGN.md` §4.1
    already describes.
  - Shift the bracket instead of re-reading both records (MD:515–522).
  - Optionally set a per-variable chunk cache (`nc_set_var_chunk_cache`) sized to the file's
    chunks.
- **F4 — validate the file shape.** Check the variables' dimension order, and reject a `time:units`
  whose unit word is not `seconds`.
- **F5 — vocabulary.** `MEDS_FORCING_DESIGN.md` §7.1 says `sw_input_kind = "fourstream"`, but the
  reader accepts `"components"`. Align the documentation.
- **F6 — optional use of `elevation(grid)`.** When `grid_match = "nearest"`, the reader could set
  `grid_elevation` from the file. The file's elevation is ignored deliberately today (MD:696–698),
  so this needs a decision, not just code.

---

## Appendix A — volume and throughput estimates (to be replaced by P0 measurements)

**Assumptions:** the `north_america` box [84, −170, 5, −50] has 791 × 1201 ≈ 0.95 M native cells.
**Land cells, estimated at 0.30–0.45 M; say 0.35 M** (P0 (a) measures this).

| Quantity | Estimate |
|---|---|
| Raw GRIB, 8 vars, 16-bit packing, ocean bitmapped | ≈ 8 × 8760 × 0.35 M × 2 B ≈ **50 GB / year** |
| Pass-A intermediates (float32, zlib) | ≈ 8 GB / month uncompressed; deleted after pass B |
| Store, uncompressed float32, 7 vars | ≈ 7 × 8760 × 0.35 M × 4 B ≈ **86 GB / year** |
| Store with quantization + zlib (3–5×, a guess) | ≈ **17–29 GB / year**, so 45 years ≈ **0.8–1.3 TB** |
| CDS requests | 12 per year (+1 trailing, +1 static per domain), so 45 years ≈ 540 |
| CDS wall time at 2–6 h per request, about 1 running | **1–3 days per data-year**, so 45 years ≈ **1.5–4.5 months** |
| Extraction read volume, 1 site × 45 years | 45 × 7 × 1.1 MB ≈ 350 MB uncompressed (a few seconds to tens of seconds) |
| Pass B memory per worker at 16k-cell bands | ≈ 5.8 GB |

## Appendix B — sources (retrieved 2026-09-25)

- CDS documentation and limits: https://confluence.ecmwf.int/display/CKB/Climate+Data+Store+(CDS)+documentation
- NetCDF cost change: https://forum.ecmwf.int/t/limitation-change-on-netcdf-era5-requests/12477
- "Request too large": https://forum.ecmwf.int/t/api-request-failing-with-your-request-is-too-large/14622
- Per-user running limit: https://forum.ecmwf.int/t/cds-requests-says-running-1-queued-1-but-no-job-is-running/15283
- Grid shift with unaligned areas: https://forum.ecmwf.int/t/lon-lat-coordinates-of-downloaded-era5-products-appear-to-depend-on-bounding-box/14602
- GRIB-to-NetCDF on the new CDS: https://confluence.ecmwf.int/display/CKB/GRIB+to+netCDF+conversion+on+new+CDS+and+ADS+systems
- ERA5-Land time series dataset: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land-timeseries
- ARCO access notebook: https://github.com/ecmwf-training/dss-notebooks/blob/main/datasets/reanalysis-era5-land/arco-access.ipynb
- Earth Data Hub ERA5-Land: https://earthdatahub.destine.eu/collections/era5/datasets/era5-land
- ERA5-Land documentation (grid, invariant files, accumulation convention): https://confluence.ecmwf.int/display/CKB/ERA5-Land:+data+documentation
- ERA5 accumulation convention: https://confluence.ecmwf.int/pages/viewpage.action?pageId=197702790
