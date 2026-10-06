# Examples

Self-contained example runs of MEDS, each in its own folder with the config, outputs, figures, and a
README that shows how to reproduce it (run from the repository root):

- **[`example01_leaf_gas_exchange/`](example01_leaf_gas_exchange/)** — the leaf gas-exchange module
  on its own, driven from Python: C3 photosynthesis, Medlyn stomata and the coupled Cᵢ solve
  reproduce the A–Cᵢ curve and the leaf-temperature responses of four tropical tree species in
  Slot & Winter (2017).
- **[`example02_canopy_phenology/`](example02_canopy_phenology/)** — the leaf-phenology module on
  its own, driven from Python: one kernel with per-PFT cue masks gives temperate deciduous,
  temperate evergreen, tropical drought-deciduous and light-driven leaf-exchanging canopies over
  synthetic climates.
- **[`example_demography/`](example_demography/)** — a 250-year demographic spin-up from near-bare
  ground (cohort/patch dynamics, succession), with the site-timeseries, per-PFT AGB, and animated
  stand-structure figures.
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
