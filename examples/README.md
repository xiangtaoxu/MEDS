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
- **[`example03_demography/`](example03_demography/)** — the demography module on its own, driven
  from Python: growth, mortality (Camac et al. 2018) and recruitment are fitted to five intervals of
  the Barro Colorado Island 50-ha plot census, and the cohort-and-patch engine runs them from the
  1985 census (it tracks the next 25 years) and from near-bare ground.
- **[`example04_column_biophysics/`](example04_column_biophysics/)** — the coupled column at the
  Barro Colorado Island flux tower, Panama: canopy radiation, leaf gas exchange, energy balances,
  plant hydraulics and soil water every 15 minutes. It starts from the plot's 2010 census, calibrates
  the fast parameters against the tower, compares five years with its fluxes, and shows ten days
  at half-hourly output.
- **[`example05_forest_regeneration/`](example05_forest_regeneration/)** — the whole model at once,
  driven from Python: a forest regrowing from bare ground at Barro Colorado Island from 1600 to 2020,
  with the column of example 04 every 15 minutes and the demography of example 03 every day, three
  PFTs (early, mid and late successional) with Panama's sun leaves and a sink limit on diameter
  growth, checked against the census's growth and the forest's carbon budget, on ERA5-Land weather
  repeated and the CO₂ of its history.
