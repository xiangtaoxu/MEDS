# Examples

Self-contained example runs of MEDS, each in its own folder with the config, outputs, figures, and a
README that shows how to reproduce it (run from the repository root):

- **[`example01_leaf_gas_exchange/`](example01_leaf_gas_exchange/)** — the leaf gas-exchange module
  on its own, driven from Python: C3 photosynthesis, Medlyn stomata and the coupled Cᵢ solve
  reproduce the A–Cᵢ curve and the leaf-temperature responses of four tropical tree species in
  Slot & Winter (2017).
- **[`example02_canopy_phenology/`](example02_canopy_phenology/)** — the leaf-phenology module on
  its own, driven from Python: one kernel, four forests, four leaf habits, every site's light from
  ERA5-Land. A deciduous broadleaf on warmth and day length (Harvard Forest: MODIS LAI, leaf fall,
  litter baskets), an evergreen pine on warmth and the bright hours (Hyytiälä: needle litter), a
  tropical tree exchanging its leaves in the bright dry season (Barro Colorado Island: litter traps)
  and a drought-deciduous dry forest on water and day length (Palo Verde: leaf litter, MODIS LAI).
  They differ only in parameter values; evergreen is a leaf-cover floor.
- **[`example_demography/`](example_demography/)** — a 250-year demographic spin-up from near-bare
  ground (cohort/patch dynamics, succession), with the site-timeseries, per-PFT AGB, and animated
  stand-structure figures.
- **[`example04_column_biophysics/`](example04_column_biophysics/)** — the coupled column at a
  flux tower: canopy radiation, leaf gas exchange, the energy balances of leaves, canopy air and
  soil, plant hydraulics and soil water, every 15 minutes, at Barro Colorado Island, Panama. The
  forcing is built from the tower's own meteorology, the stand is the 2010 census of the 50-ha plot
  (no spin-up), the fast parameters are calibrated against the tower, and both runs are compared
  with its carbon, water and energy fluxes over five years. Ten days restarted at half-hourly output
  show the states the column solves: leaf, canopy-air and soil temperatures, canopy-air CO₂, leaf
  water potential and soil moisture. The data are downloaded or read in place, never committed.
