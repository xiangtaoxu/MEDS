# Examples

Reference figures from a default 250-year MEDS spin-up
([`meds_config.toml`](../meds_config.toml), daily time step), produced by the netCDF output +
post-processing pipeline:

```bash
LD_LIBRARY_PATH=$CONDA_PREFIX/lib ./build-io/meds_io_demo meds_config.toml meds_state.nc
python post_proc/plot_site_timeseries.py meds_state.nc -o examples/meds_timeseries.png
```

- **`meds_timeseries.png`** — site-level totals over time: plant number, leaf area index,
  aboveground biomass, basal area, mean DBH, and cohort/patch structure counts.
- **`meds_timeseries_pft.png`** — per-PFT aboveground-biomass stack showing the successional
  composition (pioneer → mid → climax, the wood-density growth–mortality trade-off in action).
