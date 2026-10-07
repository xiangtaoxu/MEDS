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
- **[`example_biophysics/`](example_biophysics/)** — the fast (sub-daily) loop at hourly resolution:
  a 50-year spin-up at Ithaca NY, then one July restarted for hourly output, plotting air, canopy-air,
  tallest-cohort leaf, and soil-surface temperature. Shows the coupled canopy energy balance producing
  leaf and canopy-air temperatures that the meteorological forcing never contained. This is the one
  example that drives the **full coupled model from Python** (`meds.model.Run`), with the time loop
  on the Python side — so it also plots a spin-up trajectory sampled from the running model.
- **[`example_flux_tower_bci/`](example_flux_tower_bci/)** — forcing built from a flux tower's own
  meteorology at Barro Colorado Island, Panama: a site TOML declares the data, the tool checks it
  against the sun and the data, fills gaps with a flag on every value, and states the tower's 41 m
  heights so MEDS moves each sample to every patch's canopy-air top. Also scores its longwave gap
  fill on hidden observations. The run starts from the 2010 census of the BCI 50-ha plot, with no
  spin-up, and is compared with the tower's carbon, water and energy fluxes. The data are
  downloaded or read in place, never committed.
