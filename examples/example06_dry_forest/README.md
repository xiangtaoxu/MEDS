# Example 06 (draft): dry-forest leaf habits at Palo Verde

Work toward the Central America regional example (plan:
[`docs/dev_plans/MEDS_EXAMPLE06_DRY_FOREST_PHENOLOGY_PLAN.md`](../../docs/dev_plans/MEDS_EXAMPLE06_DRY_FOREST_PHENOLOGY_PLAN.md)).
This folder holds the Palo Verde design run, which sets the four PFTs' leaf habits before the regional run.

## What it runs

- **Stand:** the 2009 inventory of the Palo Verde dry forest, Costa Rica (`data/pv_2009_census.csv`:
  6 plots, 12 patches, 336 cohorts). The structure is frozen: trees grow and die, none recruit.
- **Weather:** ERA5-Land at the Palo Verde cell, July 2008 - December 2010.
- **PFTs** (`pft_parameters.toml`): Xu et al. (2016)'s Ds, D, BD and E, with wood densities 0.35, 0.52,
  0.74 and 0.87 g/cm3. Leaves follow the water cue alone, on each cohort's predawn leaf potential.
- **Soil:** 2 m of clay loam.

## Run

From this folder, with the ED_ERA5land archive in place:

```
../../build/meds_main meds_config_pv.toml
```

The output goes to `output/`, with daily per-cohort records and hourly leaf potential (`io_config.toml`).

## Current results

Lowest leaf cover of each PFT in the 2009 and 2010 dry seasons (leaf carbon against the preceding
September-November), and its daily minimum leaf potential in April (2009-2010 mean):

| PFT | Lowest cover, 2009 / 2010 | Target | Minimum leaf potential, April |
|---|---|---|---|
| Ds | 0.07 / 0.08 | about bare | -1.70 MPa |
| D | 0.07 / 0.10 | about bare | -2.28 MPa |
| BD | 0.46 / 0.63 | 0.25-0.50 | -2.54 MPa |
| E | 0.91 / 0.96 | 0.50-1.00 | -1.90 MPa |

- The stand's lowest monthly LAI is 3.3, against MODIS's 1.5 for the cell; its October LAI is 5.2
  (MODIS 5.3).
- The canopy thins from February-March; MODIS shows it thinning from January.
- Still open: the height allometry (Chave et al. 2014) and the regional run.
