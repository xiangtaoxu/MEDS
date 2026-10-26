# Several sites per polygon — plan

**Status:** revision 6, 2026-10-19. S0–S4 merged (#449–#453); S5 implemented in PR #454 (§10). Next: lateral water and seed dispersal between sites, a site map for a region.
**Base:** `beta` at `835458c` (the state hierarchy, PR #448).
**Decided (owner, 2026-10-18):**
- D1 — all patches of a polygon in one patch block; a site block records which patches are its own,
  the way the patch block records its cohorts.
- D7 — a patch's area is a fraction of its site, and a site's area a fraction of the polygon.
- D4 — each site has its own soil-layer grid: soil maps are finer than the forcing grid, so the sites
  of one polygon can have different soil depths and layers.
- D11 — soil-layer state is output as `(time, site, soil_layer)`; the layer geometry and each site's
  fixed properties carry no time axis.
- D2 — the polygon's state is `polygon_t` (argument `poly`); its run machinery is `polygon_run_t`.
- D12 — one run path: a run is a region of one polygon or many (S0, §3.6).
- D13 — output always writes the `polygon` dimension, a one-polygon run included.
- D14 — two forcing formats, each with its own role: the ERA5-Land archive, for a site or a region,
  and the site forcing file (`ED_default`), for a site only, read at one grid (`grid_index`). A tower's
  measurements stay in the site file, with their own measurement heights.

---

## 1. Summary

A polygon is one forcing cell. Today it holds one site, so every location in it shares one soil column
and one terrain. The owner wants several sites per polygon: each site has its own soil and surface
(soil depth, layers and texture, elevation, slope, aspect), and its area is a fraction of the polygon,
set from topography and soil maps when the user defines sites. One site with area 1 stays the default
and reproduces today's results bit for bit.

- **Shape.** The polygon's ecological state holds three blocks: the cohort block, the patch block and a
  new site block. Patches are grouped by site, and the site block maps each site to its patches with
  `patch_offset` and `patch_count`, as the patch block maps each patch to its cohorts. A patch records
  its site in `owner_site`, as a cohort records its patch in `owner_patch`.
- **One patch axis.** The fast loop stays one OpenMP loop over the patch block, now covering every
  site. A site's patches are one slice of the block.
- **Patch area is a fraction of its site, site area a fraction of the polygon.** The patch thresholds
  (`min_patch_area`, the disturbance gap) then act within a site, so a small site keeps its patches. A
  patch's share of the polygon is its site's area times its own, given by one function used by every
  polygon total; with one site of area 1 it equals the patch area exactly.
- **Each site has its own soil-layer grid.** Patches of one site share a grid, and nothing mixes two
  sites' patches, so the fast loop and the blends need no change; the soil-layer output becomes per
  site.
- **Site boundaries by slices, not checks.** Every slow operation that changes patches works on one
  site's slice at a time, as the cohort operations work on one patch's slice, so no operation can mix
  two sites' soil.
- **No interaction between sites at first.** Lateral water flow and seed dispersal between sites are
  later plans.
- **One run path (S0).** A site run is a region of one polygon: one driver, one restart and output
  layout, the `polygon` dimension always written.
- **Six steps (§8).** S0 merges the two run paths; S1 and S2 are bit-identical restructurings of the
  state with one site; S3 makes several sites run; S4 lets the user define them; S5 corrects each
  site's forcing for its terrain.

---

## 2. What the code assumes today

From a survey of `beta` on 2026-10-18:

- **`site_t` is the container.** It holds the cohort block, the patch block and the site-level
  numbers (counters, the growth clock, the growth temperature, accumulators, the slow-step scratch,
  the polygon diagnostics). 120 routines take `type(site_t)`; there are 973 `site%cohort` and 727
  `site%patch` references, about two thirds of them in tests. The kernels never see it (lint
  `kernels_no_site_state`).
- **One soil column per polygon.** `soil_column_params(cfg%soil_column, soil)` builds the soil
  parameters at 7 call sites; the fast loop builds them once into `fast_ctx%col_config`, shared by
  every patch. A region sets a polygon's soil from the maps in that polygon's config copy.
- **Terrain.** Elevation exists only as the forcing's lapse correction at ingest, zero in a region.
  There is no slope or aspect.
- **The fast loop** runs one OpenMP loop over the patches, then folds the per-patch sums in patch order,
  weighted by patch area.
- **The patch operations act on the whole block:** `sort_patches` (by age), `patch_fuse_pass` (checks
  only the disturbance type), `terminate_patches` (renormalizes areas to sum to 1), `patch_compact`,
  `apply_patch_disturbance` (one gap seeded from every patch). `apply_recruitment` uses one soil.
- **Restart** stores the site's metadata (`site_area`, counters, the growth clock and temperature); in
  a region every variable has a leading `polygon` axis.
- **Output.** 120 scalar variables named `*_site` are patch-area-weighted sums over the block, so they
  are polygon totals. `reduce_patch_column_to_site` averages soil layers by index.
- **Census files** carry a `site_id` column already; `init_from_census` keeps one site's rows.

---

## 3. The shape

### 3.1 Three blocks

| Block | One entry per | Map | Who records the owner |
|---|---|---|---|
| `site_block` | site | `patch_offset`, `patch_count` into the patch block | — |
| `patch_block` | patch | `cohort_offset`, `cohort_count` into the cohort block | `owner_site` |
| `cohort_block` | cohort | — | `owner_patch` |

The patches are sorted by site, then by age within a site, so a site's patches are one contiguous
slice, `patch_offset(is) .. patch_offset(is) + patch_count(is) - 1`. `rebuild_patch_map` derives the map
from `owner_site`, as `rebuild_cohort_map` derives the cohort map from `owner_patch` (§7.2).

### 3.2 Why this shape

- **The fast loop keeps its one patch axis.** The threads share every patch of every site, with one
  set of per-thread scratch pools and one fold.
- **A site's sums are a slice.** Per-site evapotranspiration, carbon or output means are sums over a
  contiguous range of patches, weighted by patch area.
- **A polygon total is the same fold with one more factor:** each patch weighted by its share of the
  polygon, `site_area(owner_site(ip)) * area(ip)` (§3.4).

The cost is that the slow operations that create, merge or remove patches must keep to one site. §7.2
makes that structural: each runs over a site's slice. The slow step's structural operations then run
site by site in one thread; they are a small part of the run's time.

### 3.3 The polygon and its run machinery (D2)

`site_t` no longer describes what it holds: it would hold every site of the polygon.

- **`polygon_t`**, in a new `meds_polygon_state_vars.f90`, is the polygon: the three blocks and the
  polygon-level numbers. Every routine below `src/main` takes it as `poly` (`site%patch` becomes
  `poly%patch`); the rename is mechanical and changes no result.
- **`site_block`** goes in `meds_site_state_vars.f90`, next to the per-site routines.
- **`polygon_run_t`** (src/main, `meds_polygon.f90`, replacing `meds_polygon_t`) is what the run
  needs to step and check one polygon. It cannot be part of `polygon_t`: the fast-loop context, the
  slow ledger and the soil-carbon seam are defined in the fast and slow libraries, which use the state
  library, so a state type holding them would need those libraries compiled first, and they it. The
  two sit side by side in the region (§3.6): neither holds the other, and `polygon_run_t` does not
  extend `polygon_t` (owner, 2026-10-18).

  | Role | Fields |
  |---|---|
  | forcing | `met_cur` (the cursor: cell, forcing-cell elevation, the two bracketing records, the recycle seam, the last daytime clearness index), `cell`, `label` |
  | fast-loop set-up | `fast_ctx` (column configuration, acclimated leaf table, starting soil state, canopy optics, ground albedo and emissivity, patch threads); `cfg`, the polygon's config copy, until S2 |
  | conservation checks | `energy_budget`, `water_budget`, `face_budget`, `slow_ledger`, `seam` with the worst seam's date and patch count |
  | run report | `most_patches`, `most_cohorts` and when, `fast_step_total`, `area_start`, `status` |

  To `polygon_t` instead, as the polygon's own: `id`, latitude and longitude (the cursor reads them
  rather than keeping copies), and `restructure_pending` / `restructure_new_year`, which a checkpoint
  stores. `elevation_m` stays with them until S2 makes it per site. The owner will revisit
  `fast_ctx` later in the refactor; S2 already moves its soil parameters to the sites and `cfg` to the
  region.

This revisits the state plan's §3.6 ("no polygon file"): that file was declined because the run's
polygon needs the forcing and the fast and slow libraries. `polygon_t` holds state only and lives in
the state library, below all of them.

### 3.4 Areas (D7)

- **A patch's area is a fraction of its site**, summing to 1 over the site's patches, as it sums to 1
  over the polygon's today. Every within-site operation reads it unchanged: the patch thresholds
  (`min_patch_area`), fusion, termination's renormalization and the disturbance gap's area. A small
  site keeps its patches; were patch area a fraction of the polygon, every patch of a site covering 1 %
  of the polygon would sit below `min_patch_area` and be merged away.
- **A site's area is a fraction of the polygon**, summing to 1 over the sites. It changes only when the
  user's site definition does (and later, perhaps, by land-use change).
- **A patch's share of the polygon** is `site%area(owner_site(ip)) * patch%area(ip)`. It is computed
  where it is needed, by one function (`patch_polygon_area`), never stored: it is cheap, read once per
  output tick or slow step, and a stored copy would have to be refreshed after every area change.
- **Every polygon total uses it:** about 45 lines today read `patch%area` to form a polygon total (the
  output's reductions, the fast loop's fold, the slow ledger, the budgets' merge, the C API's totals).
  With one site of area 1, the product equals the patch area exactly, so switching them is
  bit-identical (S1).

---

### 3.5 How the levels fit together

```
meds_region_t                       src/main: the run, one polygon or many
 ├─ cfg, met_src, out_files, now, ...         shared by every polygon
 ├─ poly(:) : polygon_t             meds_polygon_state_vars    each polygon's state  (below src/main: poly)
 │   ├─ site   : site_block         meds_site_state_vars       one entry per site
 │   │    area(:)                     fraction of the polygon
 │   │    patch_offset(:), patch_count(:)   its patches: one slice of the patch block
 │   │    surface(:)                  soil column and its parameters, elevation, slope, aspect (S2)
 │   ├─ patch  : patch_block        meds_patch_state_vars      one entry per patch, sorted by (site, age)
 │   │    owner_site(:)               its site
 │   │    area(:)                     fraction of its site
 │   │    cohort_offset(:), cohort_count(:) its cohorts: one slice of the cohort block
 │   │    cas, snow, soil_water, soil_energy, soil_carbon, ...   its stores, as today
 │   ├─ cohort : cohort_block       meds_cohort_state_vars     one entry per cohort, sorted by (patch, height)
 │   │    owner_patch(:)              its patch
 │   │    nplant(:)                   per m2 of its patch, as today
 │   └─ the polygon's own numbers: id, latitude, longitude, n_pft, the id counters, growth_hist_pos,
 │      t_growth_avg, pheno_tair_sum/n, the pending restructuring, deriv, the step accumulators, diag (PY_*)
 ├─ run(:)  : polygon_run_t         each polygon's machinery: forcing cursor, fast context, budgets,
 │                                    ledgers, status, counters
 └─ out_bufs(:)                     each polygon's share of the output files, as today
```

- **Two maps, one pattern.** Each level records its owner (`cohort%owner_patch`, `patch%owner_site`),
  and the level above holds the map back to it as an offset and a count. `rebuild_cohort_map` (today's
  `rebuild_csr`) builds the patch → cohort map; `rebuild_patch_map` builds the site → patch map.
- **Three area levels.** A cohort's `nplant` is per m2 of its patch, a patch's `area` a fraction of its
  site, a site's `area` a fraction of the polygon; `patch_polygon_area` multiplies the last two.
- **Loops.**

  ```fortran
  do ip = 1, poly%patch%n                                   ! the fast loop: one axis, every site
     is = poly%patch%owner_site(ip)                         ! the patch's soil is poly%site%surface(is)
  end do

  do is = 1, poly%site%n                                    ! a slow structural operation
     do ip = poly%site%patch_offset(is), poly%site%patch_offset(is) + poly%site%patch_count(is) - 1
        do i = poly%patch%cohort_offset(ip), poly%patch%cohort_offset(ip) + poly%patch%cohort_count(ip) - 1
        end do
     end do
  end do
  ```

- **Names.** `site_t` is gone: `polygon_t` takes its blocks and numbers, `site_block` its area.
  `site_alloc` and `site_free` become `alloc_polygon` and `free_polygon`; the site block has `alloc_`,
  `expand_`, `reorder_` and `free_site_block`, as the other blocks do.
- **In `src/main`** a polygon is `region%poly(ip)` and its machinery `region%run(ip)`; below `src/main`
  the polygon is the argument `poly`, so a field reads `poly%patch%area(ip)` everywhere.
- **The kernel lint** adds `polygon_t`, `site_block` and `meds_polygon_state_vars` to what a kernel
  may not name.

### 3.6 One run path (D12, D13, D14)

Today a site run (`meds_run_t`, `meds_driver`) and a region run (`meds_region_t`, `meds_region`) are
two drivers with the same parts: the config, the forcing reader, the output files, the clock, and one
polygon or an array of them (`meds_run_t` already holds its output buffers as "an array of one, as a
region's"). `meds_main` picks one by `[run].mode`, and the restart and output code branch on which.

- **A run is a region of one polygon or many.** `meds_region_t` is the run; `meds_driver` and
  `meds_run_t` go, and so do `meds_main`'s mode switch and the site-or-region branches of the restart
  and output code (about 25 lines).
- **Stepping.** The region steps a month at a time, its polygons in parallel; the C API and Python
  step one slow step at a time. The region gains a one-step entry beside the monthly one.
- **Forcing (D14).** Two formats, as today. The ERA5-Land archive serves a site (the cell nearest the
  `[site]`) or a region (the cells of a box). The site forcing file (`ED_default`) serves a site only,
  at one grid (`grid_index`, or the grid nearest `[site]`): a tower's measurements, or a reanalysis
  extracted at a point. Its heights (`[forcing].tq_height`, `wind_height`, `height_above`) are the
  run's, which suits one site's instruments; a region over a file's grids would need them per grid, so
  a region stays on the archive (owner, 2026-10-18).
- **Output and restart (D13)** always carry the `polygon` dimension, a one-polygon run included, as the
  region files do today: one layout for every run, in the order ED2's output also follows (polygon,
  site, patch, cohort). A one-polygon run's files change layout once; their values do not change.
- **The polygon's own config copy.** Each polygon carries a copy of the config today, because the soil
  maps had nowhere else to go. Once each site carries its own soil (S2), the region keeps one config.
- **Why a record per polygon, not one array per field.** The region steps its polygons on parallel
  threads, and each thread uses one polygon's machinery whole. Separate arrays would put neighbouring
  polygons' counters and sums on cache lines that other threads write, and no code reads one of these
  fields across polygons in a hot loop.

## 4. What lives where (D3)

| Field | Today | After | Why |
|---|---|---|---|
| cohort and patch blocks | `site_t` | `polygon_t` | one block each per polygon |
| `owner_site` | — | patch block | the patch's site |
| `patch_offset`, `patch_count` | — | site block | the site → patch map |
| site area | `site_t%site_area`, unused | site block, a fraction of the polygon | patch area is a fraction of its site (§3.4) |
| surface: soil column, elevation, slope, aspect | config (soil); polygon (elevation) | site block, a `site_surface_t` per site (§5) | what defines a site |
| `next_cohort_id`, `next_patch_id` | `site_t` | `polygon_t` | ids unique within the polygon, so a tracker needs only the id |
| `growth_hist_pos` | `site_t` | `polygon_t` | every cohort's shared sample clock |
| `deriv` | `site_t` | `polygon_t` | indexed by cohort |
| `t_growth_avg`, `pheno_tair_sum/n` | `site_t` | `polygon_t` through S4; site block in S5 | they follow the air temperature, which differs between sites only once S5 corrects it |
| `et_accum`, `runoff_surf_accum`, `drainage_accum`, `work_*` | `site_t` | `polygon_t` | polygon sums over every patch, unchanged |
| forcing echo (`diag`, `PY_*`) | `site_t` | `polygon_t` | one forcing per polygon; `SI_*` rows on the site block from S5 |
| `n_pft` | `site_t` | `polygon_t` | one PFT set per polygon |

The output-only accumulators move into the diagnostic blocks as their own change (the state plan's §9
follow-up), because that changes the output arithmetic.

---

## 5. A site's surface

A `site_surface_t` per site:

- **The soil column:** the `soil_column_config_t` the config holds today (depth, layer bottoms,
  per-layer texture, retention choice), the bottom boundary and drain factor, the deep temperature and
  its depth. Its parameters (`soil_params_t`, `soil_thermal_params_t`) are derived once, when the site
  is made, and kept with it.
- **The terrain:** elevation [m], slope [degrees], aspect [degrees from north]. Unused until S5.
- **To decide (D10):** the ground albedo and emissivity (today `[soil]`).

Consequences:

- **The 7 calls of `soil_column_params(cfg%soil_column, ...)`** read the patch's site's parameters.
- **The fast-loop context holds one soil column per site**, and each patch's work picks its site's
  (`owner_site`). The kernels are unchanged: they take soil parameters as arguments already.
- **Each site has its own soil-layer grid (D4).** A soil map is finer than the forcing grid, so the
  sites of one polygon can differ in soil depth and layering. Within a site every patch shares the
  grid, and no operation mixes two sites' patches, so the fast loop, the blends and the patch stores
  (already sized `n_soil_layer_max`) need no change. What changes:
  - the slow ledger's soil geometry is the patch's site's, not one per polygon;
  - the soil-layer output (7 soil and 7 soil-by-patch variables) averages patches by layer index,
    which is meaningless across grids, so it is reported per site (§7.4);
  - the region's map echo (`polygon_soil_z`, `_face`, `_depth`) becomes per site.
- **The surface is configuration held on the site.** It is never integrated, blended or reordered; the
  restart stores it so a resumed run can check it against its config.

---

## 6. How a user defines sites (D5)

- **A site run.** A `[[site]]` list in the run's TOML, one table per site: `area_fraction`, the soil
  column keys `[soil]` takes today, and `elevation`, `slope`, `aspect`. Without the list, the run has
  one site with area 1 and today's `[soil]`, so every existing config keeps working.
- **A region.** A site map gives each cell its sites and their fractions, and each site its soil and
  terrain, read by an extension of `apply_map_inputs`. Which product defines sites is the owner's
  choice: soil map units, height above the nearest drainage, topographic wetness or slope classes, or
  their combination. S4 builds the interface with one site per cell, so a region run is unchanged until
  a real site map exists.
- **A census.** The `site_id` column names each row's site.

---

## 7. Through the code

### 7.1 Fast loop

- One OpenMP loop over the patch block, as today. Each patch's work reads its site's soil column
  through `owner_site`.
- The per-patch sums fold in patch order, as today, weighted by each patch's share of the polygon
  (§3.4); a per-site sum is the same fold over the site's slice, weighted by patch area.
- From S5, the forcing is corrected per site after `met_instant` (elevation lapse; beam radiation on a
  sloped surface), and each patch takes its site's.

### 7.2 Slow loop: every structural operation works on one site's slice

| Operation | Today | After |
|---|---|---|
| `sort_patches` | by age, whole block | by (site, age); then `rebuild_patch_map` and `rebuild_cohort_map` |
| `patch_fuse_pass` / `fuse_2_patches` | every pair in the block | every pair within one site's slice |
| `terminate_patches` | merges into the largest survivor; areas renormalized to sum to 1 | per site: into the site's largest survivor; areas renormalized to the site's area |
| `apply_patch_disturbance` | one gap from every patch | one gap per site, from that site's patches |
| `patch_compact` | one keep mask | the same, then `rebuild_patch_map` |
| `apply_recruitment` | one soil | each patch's site's soil |
| cohort operations, canopy depth, soil carbon | one patch's slice | unchanged |

The pattern is the one the cohort operations already follow per patch: an outer loop over sites, an
inner loop over the site's slice. A routine that never sees more than one site's slice cannot blend two
soil columns.

### 7.3 Restart

- A `site` dimension (padded behind `polygon` in a region) with each site's area and surface, and
  `owner_site` per patch. A file without them reads as one site with area 1.

### 7.4 Output (D6)

- **The `*_site` variables** are polygon totals today and stay so, weighted by each patch's share of
  the polygon (§3.4). The recommendation is to keep their names and say "polygon" in their
  descriptions.
- **Per-site output** is an optional `site` dimension for the stand-level variables, each a fold over a
  site's slice weighted by patch area.
- **Soil-layer output (D11, decided):** the state that changes is `(time, site, soil_layer)` (soil
  water, potential, temperature, liquid fraction), a one-site run included. The layer geometry does not
  change, so it has no time axis, as today: `soil_z` and `soil_z_bnds` become `(site, soil_layer)`,
  `(polygon, site, soil_layer)` in a region.
- **A site's fixed properties** are written once, with no time axis: its area fraction, elevation,
  slope, aspect, soil depth and per-layer texture, as a region echoes each polygon's map inputs today.
- **Patch output** gains `patch_site`; **cohort output** keeps `owner_patch`, which leads to the site.

### 7.5 C API and Python

- Entry points that read the stand gain an optional site index (0 for the polygon). Python's `Run` gains
  `n_site` and per-site views.

---

## 8. Steps (each its own pull request)

| Step | What | Results |
|---|---|---|
| **S0** | One run path (§3.6): `meds_region_t` runs one polygon or many; `meds_driver` and `meds_run_t` go; `meds_polygon_t` becomes `polygon_run_t` beside the state; a one-step entry for the C API; the `polygon` dimension in every output and restart file. | Values bit-identical; one-polygon files gain the `polygon` dimension |
| **S1** | `polygon_t` holds the three blocks and the polygon-level numbers; `site_block` with one site of area 1 covering every patch; `owner_site` = 1; `sort_patches` sorts by (site, age); `rebuild_patch_map`. The 120 routines take `polygon_t` (`site%` → `poly%`). Every polygon total weights by `patch_polygon_area` (§3.4). | Bit-identical |
| **S2** | Each site owns its surface (§5), its own layer grid included: the soil column moves from `cfg` to the site block, the fast context and the slow ledger hold one soil column per site, the 7 call sites read the patch's site's. | Bit-identical |
| **S3** | Several sites run: the structural operations by site slice (§7.2), initialization per site (bare ground or census by `site_id`), the restart's `site` axis, `patch_site`, per-site soil-layer output (D11) and optional per-site stand output. Tests: two identical sites with areas 0.3 and 0.7 reproduce one site's polygon output, to round-off; two soils give two different stands; thread-count invariance. | Bit-identical with one site |
| **S4** | Users define sites: `[[site]]` in a site run, a site map in a region (one site per cell until a real map exists). | Bit-identical with one site |
| **S5** | Terrain corrects each site's forcing: the elevation lapse per site, slope and aspect on the beam radiation (a new section of `docs/science/forcing.md`); the growth temperature and the phenology's air temperature move to the site block; `SI_*` rows. | New science |
| later | Lateral water flow from upslope to downslope sites; seed dispersal between sites. | Their own plans |

S0 changes the drivers and the I/O; S1 is the largest diff (the rename touches every routine that
takes a site) and the simplest to verify; S3 holds the new behaviour.

---

## 9. Decisions

| # | Decision | Recommendation |
|---|---|---|
| D1 | Shape | **Decided:** one patch block per polygon, a site block mapping sites to patches |
| D2 | The polygon | **Decided:** `polygon_t` (state, argument `poly`) beside `polygon_run_t` (machinery, src/main), both in the region; no type extension (§3.3) |
| D3 | Where each field lives | §4 |
| D4 | Soil-layer grid | **Decided:** each site its own |
| D5 | How sites are defined | **Taken (S4, on the recommendation):** `[[sites]]` in a site run (the name departs from §6, see §10); a site map in a region, its product the owner's choice (§6), one site per cell until then |
| D6 | Output names | **Taken (S3, on the recommendation):** keep `*_site` as the polygon totals; `patch_site`; the optional per-site stand output deferred |
| D7 | Areas | **Decided:** patch area a fraction of its site, site area a fraction of the polygon |
| D8 | Ids | **Taken (S3, on the recommendation):** unique within the polygon |
| D9 | Interaction between sites | **Taken (S3, on the recommendation):** none at first; lateral water and seed dispersal later |
| D10 | Ground albedo and emissivity | On the site's surface, or shared by the polygon |
| D11 | Soil-layer output layout | **Decided:** state `(time, site, soil_layer)`; geometry and site properties without time |
| D12 | Run paths | **Decided:** one; a run is a region of one polygon or many (S0) |
| D13 | `polygon` dimension | **Decided:** always written |
| D14 | Forcing formats | **Decided:** two: the archive (site or region) and the site file (`ED_default`, site only, one grid); a region stays on the archive |

---

## 10. Implementation log

**S0, PR #449 (2026-10-18).** Four commits on `refactor/one-run-path`, each checked against `beta` at
`835458c` on seven cases: a two-year census run, a one-year run with every output group on, example
04's evaluation month, a Palo Verde inventory site on the ERA5-Land archive with the soil maps, a
four-cell region, a site restart from an old-format state file and a region restart. Every variable
the old files hold matches bit for bit; the new files add only the polygon coordinates. A site
restarted from a new one-polygon state file runs as one restarted from the old format.

1. One run path: `meds_region` runs every run, `meds_driver` is deleted, `region_step` is the C API's
   one-step entry, and a one-polygon run writes the polygon axis.
2. The output and state writers lose their site-only branches; `manager_set_polygons` lays out a file
   set's polygon axis.
3. `polygon_run_t` beside the state: the region holds `poly(:)` (a `site_t` until S1) and `run(:)`.
4. The scripts, tests and documents read the polygon axis.

Departures from the plan:

- **`ED_default` sources in a region are dropped** (owner, 2026-10-18; D14). The site file stays a
  site's: one grid, read by `grid_index`, with the run's measurement heights. A region needs the
  ERA5-Land archive, as before.
- **`id`, the location and the restructure flags stay in `polygon_run_t`** until S1 creates
  `polygon_t` (§3.3).
- **The `meds_driver` CMake library keeps its name**; it now builds `meds_polygon` and `meds_region`.

**S1, PR #450 (2026-10-19).** Five commits on `refactor/polygon-state`, each checked against the merged
S0 build (`5e40909`) on S0's seven cases and a site restarted from its own one-polygon state file:
every variable of every output and state file matches bit for bit.

1. `polygon_t`, renamed from `site_t`, in `meds_polygon_state_vars`; its argument is `poly`;
   `alloc_polygon`, `free_polygon` and `rebuild_cohort_map` (was `rebuild_csr`). The kernel lint
   forbids the new names. Every changed code line is the rename, string literals byte-identical.
2. The site block (`meds_site_state_vars`): one site of area 1, its patch map; `owner_site` on the
   patch block; `sort_patches` by (site, age); `rebuild_patch_map` wherever the patch list changes.
   The unused `site_area` scalar is the site's area.
3. The polygon's id, location and owed restructuring on `polygon_t`; `io_read_state` sets the
   restructuring owed on the polygon it reads.
4. `patch_polygon_area` in every polygon total: 48 lines in the output reductions, the fast loop's
   staging and fold, the budgets' merge, the slow ledger and its declarations, the recruitment and
   termination totals, the disturbance's shed water and the C API's soil carbon. The output
   reductions take the shares as one contiguous array (`patch_polygon_areas`): with the function
   inside their loops, the compiler summed some patch reductions in a different order, changing up to
   38 daily `*_site` variables by up to 2e-13 relative; with the array the loops are what they were.
5. The documents.

Notes:

- **Names.** The map builders are §3.5's: `rebuild_patch_map` and `rebuild_cohort_map` (§3.1 and §7.2
  said `rebuild_site_map`; corrected). The site block has `alloc_` and `free_site_block`; `expand_` and
  `reorder_site_block` come with their first caller (S3, S4).
- **`owner_site` defaults to 1** when the patch block is allocated, so a patch built by hand (tests)
  is in the one site; the disturbance gap sets it explicitly, because a reused slot can hold a moved
  patch's value. A state file has no site axis until S3: every patch reads into site 1, and
  `meta_real(1)` carries that site's area.
- **The forcing cursor keeps its own copy of the location** (it lives in the forcing library); it is
  set from `polygon_t` when the cursor is made.
- **Found for S3 while switching the totals.** These act on the whole polygon today and must act
  within a site:
  - `patch_fuse_pass` pairs patches by disturbance type alone; it must also match `owner_site`.
  - `terminate_patches` merges a sliver into the polygon's largest survivor and renormalizes the
    areas over the polygon; per site, both.
  - `apply_patch_disturbance` cuts one gap from every patch; one gap per site.
  - The slow ledger takes one soil column for every patch; per site (S2).
  - The output `area_patch` and the C API's patch-area export give a patch's own area, a fraction of
    its site; Python's `Run` uses it as a polygon weight, so it needs `owner_site` and the site areas
    beside it (S3, with `patch_site`).

**S2, PR #451 (2026-10-19).** Three commits on `refactor/site-surface`, each bit-identical to the
merged S1 build (`6fa01cc`) on the eight cases, the region of cells with different soils included.

1. `site_surface_t` on the site block (`meds_site_state_vars`): the soil column with its derived
   parameters and each layer's starting water, the column's bottom (water boundary, drainage factor,
   deep temperature and depth), and the terrain (elevation, slope, aspect). `make_site_surface` makes
   it from a configuration wherever a polygon is made from one (bare ground, a census, a restart); a
   region sets each site's elevation to its cell's.
2. Every reader takes the patch's site's surface: the fast loop's column configuration
   (`site_column_configs`, made at each fast step from the polygon's, so each copy carries the
   acclimated leaf table), the fast reservoirs (`theta_init` leaves the fast context), the slow
   ledger's stores, a new cohort's birth state (`patch_physical_state`), recruitment and the
   restart's tissue fallback. The output echo reads site 1's surface.
3. The region keeps one configuration: `polygon_run_t%cfg` goes, and the polygon routines take the
   run's. A cell's own inputs reach only the initial community (its sites' surfaces) and, through
   the cohort block's allometry, the height allometry's site index.

Notes:

- **The thermal bottom boundary stays run-wide** (`[soil].thermal_bottom_bc`; the soil maps do not set
  it). **The ground albedo and emissivity stay run-wide too** (D10): `[soil]` sets them for every site.
- **A run without the fast loop** now writes its real soil layers in `soil_z` (it wrote the unbuilt
  column of a fast context it never made).
- **The restart does not store the surface yet:** it is made again from the configuration and maps
  (S3 adds the site axis).

**S3, PR #452 (2026-10-19).** Five commits on `feat/several-sites`, each checked against the merged
S2 build (`f8fd973`) on the eight cases: every variable of every output and state file matches bit
for bit, the soil columns behind a site axis of length 1 (the comparison drops a length-1 `site`
axis as it drops the `polygon` one).

1. Several sites from an initial community, and the structural operations by site slice (§7.2):
   `alloc_polygon` takes a site count, `init_bare_sites` builds bare sites with their own
   configurations and areas, and `init_from_census` reads one site per listed `site_id`; fusion pairs
   patches of one site, termination merges into the site's largest patch and renormalizes each
   site, a disturbance cuts one gap per site, and `sort_patches` puts the gaps in their sites.
   `test_several_sites`: two identical sites of areas 0.3 and 0.7 each evolve as a one-site polygon
   over 60 days with two disturbances and four restructurings, to 1e-12; two soils give two stands,
   each its one-site polygon's; four patch threads give one thread's polygon bit for bit.
2. The restart's site axis (§7.3): `n_site`, each site's area and surface echo (soil layers and
   depth, elevation, slope, aspect), `owner_site` per patch. `io_read_state` makes each site from
   `site_cfg` (the run's configuration when absent) and stops on a site count or a site's layer
   count that differs from the configuration's. A file without the axis is one site.
3. Output (§7.4, D11): a `site` axis in every file; the soil columns `(time, polygon, site, soil)`,
   each site's on its own layers, its patches weighted by their areas in the site; `soil_z`,
   `soil_z_bnds` and the soil echoes per site; `site_area`, `site_elevation`, `site_slope`,
   `site_aspect`, `n_site`; `patch_site`. The soil potential reads each patch's site's curve from its
   surface, so the output's copy of one curve goes.
4. Comment boxes.
5. C API and Python (§7.5): `n_site`, the sites' areas, each patch's site and share of the polygon.

Decisions taken on the plan's recommendations (the owner asked for S3–S5 without stopping):

- **D6:** the `*_site` scalars keep their names and stay the polygon's (`docs/science/diagnostics.md`
  §2 says so); the soil columns are per site, as D11 decided.
- **D8:** ids stay unique within the polygon (`next_patch_id`, `next_cohort_id` are the polygon's).
- **D9:** no interaction between sites.

Notes:

- **The fast tier's soil columns** weight a patch by its area in its site (the pre-S1 weight, bit
  for bit); the polygon rows keep `patch_polygon_area`.
- **Deferred:** the optional per-site stand output (§7.4: a `site` axis for the stand variables) has
  no user until sites can be defined; the soil texture echo per layer (§7.4) waits for a site map
  that varies it. Scripts that weight cohorts by a state file's `patch_area` treat it as the
  polygon's share, right for one site.

**S4, PR #453 (2026-10-19).** One commit on `feat/site-definitions`, checked against the merged S3
build (`0fb0ce8`) on the eight cases (none defines sites): every variable matches bit for bit.

- **`[[sites]]`, not `[[site]]`** (departure from §6): `[site]` is already the polygon's place, and a
  TOML file cannot hold a table and an array of tables of one name. Each table: `area_fraction`
  (required, the shares summing to 1 within 1e-6, made exact when the polygon is built), `elevation`
  (default the polygon's), `slope`, `aspect` (default flat), `census_site_id` (required in a census
  run) and any of `[soil]`'s column keys under their own names (`SITE_GROUND_KEYS`: the column, its
  water bottom, `deep_temp`, `deep_depth`). A table with any of them has its own ground,
  `[soil]`'s values for the rest; one without takes the polygon's (so the soil maps' cell soil in a
  mapped site run, where those keys are refused).
- **Reading.** `meds_toml` numbers the tables of `[[name]]` (`sites.2.depth`); `check_config_keys`
  checks `sites.2.depth` as `sites.depth` or as a ground key; the CMake key lister reads `[[name]]`
  headers. Each table is read through a view under `[soil]`'s names (`site_view`), so the `[soil]`
  readers serve it and the parameter record lists its values with `<file> [[sites]] <i>` as the
  source.
- **Building.** `site_config(cfg, is)` is the polygon's configuration with the site's ground; the
  region passes the sites to `init_bare_sites`, `init_from_census` (`site_ids`) and `io_read_state`
  (`site_cfg`), and `set_site_terrain` sets each site's elevation, slope and aspect.
  `validate_config` checks the shares, terrain and census ids, and validates each site's own ground
  as the run with that ground.
- **A region refuses `[[sites]]`** (D5): the site-map interface is a cell's one site, its config the
  cell's; a real map adds sites in `apply_map_inputs` once the owner picks its product.
- **Python:** `meds.config.dumps` writes a list of tables as `[[name]]` tables.
- **Tests:** `test_site_definitions` (two sites from example 04's config: the config, the polygon,
  its output's site axis, its state file and a restart from it, a census by `census_site_id`, and
  four refusals by `meds_main`); `test_region` refuses `[[sites]]` in a region and a site's soil key
  beside the maps; `test_config.py` round-trips `[[sites]]`.

**S5, PR #454 (2026-10-19).** One commit on `feat/site-terrain`, checked against the merged S4 build
(`0a96c77`) on the eight cases: every variable of every output and state file matches bit for bit
(`t_growth_avg` moves from the meta axis to the site axis, its one value unchanged), the new
`*_by_site` variables aside.

- **The forcing per site** (`docs/science/forcing.md` §8, Each site of a polygon): the fast loop makes
  each site's forcing per sub-step from the polygon's sample (`met_at_site`): the terrain lapse from
  the polygon's elevation to the site's, and the beam on the slope (`slope_incidence_cosz`,
  `slope_beam_factor`; the beam normal to the sun held at the solar constant). A flat site at the
  polygon's elevation copies the sample, so a run without sites is unchanged bit for bit. The
  diffuse light, wind and precipitation are the polygon's, and the canopy keeps the sun's zenith.
- **The site's climate**: `t_growth_avg`, `pheno_tair_sum` and `pheno_tair_n` move from `polygon_t`
  to the site block; the phenology reads its cohort's site's daily mean; each site's leaf table is
  acclimated to its own growth temperature (`fast_context_t%site_leaf_photo`, which
  `site_column_configs` hands each site's column).
- **SI_* rows**: a per-site diagnostic block (`site_diag_block`) holds each site's forcing; the
  output's new `DIM_SITE` axis kind writes it and the growth temperature as `*_by_site` variables.
- **Not done**: the diffuse light's sky view and the light reflected by the surrounding terrain (a
  slope sees less sky and some ground); lateral water and seed dispersal between sites (§8, later).

