# `src/io/input/maps/` — static per-cell maps

The one place MEDS opens a map: a static, per-cell input on the forcing grid, as opposed to the
time-varying drivers in [`../forcing/`](../forcing/README.md). Library `meds_maps`; it links
`meds_forcing` for the archive grid it checks maps against and the cell lists it reads them at.
The design is `docs/dev_plans/MEDS_POLYGON_RUNTIME_PLAN.md` §14 (RP2, RP3).

| Module | What it does |
|---|---|
| `meds_surface_maps` | Reads a map variable at a list of cells: 2-D `(lat, lon)` or layered `(layer, lat, lon)`, one hyperslab over the rows the cells span. Applies a region's cell mask. Reads the `[surface]` maps at a run's cells (`read_surface_inputs`) and turns one cell's values into its polygon's soil column, bottom boundaries and Chave E (`apply_surface_inputs`; `site_surface_inputs` for a site). |

## Rules

- **A map is on the forcing grid.** The file's `lat` and `lon` axes must be the ED_ERA5land static
  file's, so a cell's value sits at the cell's own `(row, col)`. MEDS never regrids: a map on
  another grid is refused (`MAP_ERR_GRID`).
- **An empty cell reads as NaN.** The variable's `_FillValue` or `missing_value`, or a NaN, marks a
  cell the map does not hold. Filling gaps is the product's job, in its own pipeline, not MEDS's; the
  caller decides whether a missing value is an error. Packed values (`scale_factor`, `add_offset`)
  are unpacked.
- **Nothing else opens a map.** A new static input is a new routine here, read the same way.

## Used by

- **The cell mask** (`[region].cell_mask`, `cell_mask_var`, `cell_mask_min`): `meds_region` keeps the
  cells of its box whose mask value is at least the threshold, through `apply_cell_mask`, which
  regroups the kept cells by archive chunk (`era5land_keep_cells`).
- **The surface inputs** (`[surface]`): `meds_region` reads them for every polygon once and gives each
  polygon its cell's; `meds_driver` reads them at a site's forcing cell. A run cell a map leaves empty
  is an error that names the map, the variable and the cell (OR17). The products' variable names are
  fixed here: ED-SoilMaps v2 (`soil_depth`, `bc_class`, `deep_temp`; `sand`, `silt`, `clay`,
  `bulk_density`, `coarse_fragments`; `organic_carbon`, `cec`, `ph_h2o`, on `depth` with `depth_bnds`)
  and Chave et al. (2014)'s `chave_e`.
