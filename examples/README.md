# Examples

Self-contained example runs of MEDS, each in its own folder with the config, outputs, figures, and a
README that shows how to reproduce it (run from the repository root):

- **[`example_demography/`](example_demography/)** — a 250-year demographic spin-up from near-bare
  ground (cohort/patch dynamics, succession), with the site-timeseries, per-PFT AGB, and animated
  stand-structure figures.
- **[`example_leaf_gas_exchange/`](example_leaf_gas_exchange/)** — the standalone leaf-level
  photosynthesis + stomatal-conductance module: A–Ci, A–PAR, A–temperature and gs–VPD response curves
  (FvCB C3 / Collatz C4; Leuning / Medlyn / Katul stomata). Not yet coupled to the demographic spin-up.
- **[`example_biophysics/`](example_biophysics/)** — the fast (sub-daily) loop at hourly resolution:
  a 50-year spin-up at Ithaca NY, then one July restarted for hourly output, plotting air, canopy-air,
  tallest-cohort leaf, and soil-surface temperature. Shows the coupled canopy energy balance producing
  leaf and canopy-air temperatures that the meteorological forcing never contained. (See its README for
  a known bug in the soil-surface trace during rain.)
- **[`example_phenology/`](example_phenology/)** — the leaf-phenology kernel driven over four synthetic
  climates, reproducing the four strategies (temperate deciduous / evergreen, tropical drought-deciduous
  / light-driven leaf-exchanging): relative LAI + the flush and shed rate tendencies for each.
