# MEDS — Modular Ecosystem Dynamics Simulator

*(equally: **Modern Ecosystem Demography Simulator**, a nod to its descent from ED / ED2.)*

MEDS is a ground-up reimplementation in **Fortran 2018** of the Ecosystem Demography model
[ED2](https://github.com/EDmodel/ED2). It keeps ED2's size- and age-structured representation of a
forest — cohorts of plants, patches of land, fusion and fission to keep both bounded — and puts a
**solved land-surface physics loop** underneath it, at a sub-daily timestep.

The distinction that matters to an ecologist: MEDS does not read leaf temperature, canopy humidity
or soil moisture from a driver file and does not prescribe growth. It solves a coupled surface
energy, water and CO₂ balance every sub-daily step, and the demography runs on the carbon that
balance produces.

<p align="center">
  <img src="examples/example_demography/example_output_forest.gif" height="230" alt="Canopy-layer stand profile with vertical LAI">
  &nbsp;&nbsp;
  <img src="examples/example_demography/forest3d_growth.gif" height="230" alt="3D landscape, trees growing in place">
</p>

*A 250-year spin-up from near-bare ground. Left: the vertical leaf-area profile beside a stand
cross-section, each bar a cohort's canopy disk seen edge-on. Right: the same run as a synthetic
landscape, with every tree tracked by a persistent id so it grows in place. Colour is plant
functional type — green pioneer, blue mid-successional, magenta climax — and the sequence is the
textbook ED succession: bare ground, a pioneer flush, canopy closure. Reproduce it from
[`examples/example_demography/`](examples/example_demography/).*

[![Hourly canopy energy balance for one July](examples/example_biophysics/biophysics_july.png)](examples/example_biophysics/)

*One July at hourly resolution, year 50 of a run at Ithaca NY. **Only the grey curve is an input** —
above-canopy air temperature from the ERA5-Land forcing. The canopy air space, the sunlit leaf and
the soil surface are all solved from it and the incoming shortwave. They separate from the forcing
in different directions and with different phase, which is what a real surface does and what a
meteorological file cannot tell you. The same run's carbon and soil-water figures come from the
same output stream, not three separate studies:
[`examples/example_biophysics/`](examples/example_biophysics/).*

## What is in the model

| | Scheme | Documented in |
|---|---|---|
| Canopy radiation | ED2 two-stream (`icanrad = 2`), multi-band, Beta leaf-angle distribution | [canopy_radiation_transfer](docs/science/canopy_radiation_transfer.md) |
| Aerodynamics | CLM5 Monin-Obukhov, ED2 Nusselt leaf and wood boundary layers | [canopy_aerodynamics](docs/science/canopy_aerodynamics.md) |
| Canopy air space | three prognostic twins: enthalpy, humidity, CO₂ | [canopy_air_space](docs/science/canopy_air_space_biophysics.md) |
| Photosynthesis | FvCB C3 and Collatz C4, three stomatal models, two-limb water stress | [leaf_gas_exchange](docs/science/leaf_gas_exchange.md) |
| Plant hydraulics | matrix-exponential network solve, multi-layer root uptake | [plant_hydraulics](docs/science/plant_hydraulics.md) |
| Soil water | implicit Richards, van Genuchten or Campbell retention | [soil_biophysics](docs/science/soil_biophysics.md) |
| Energy and snow | **internal energy, not temperature** — freeze/thaw is a read-off, not a branch | [column_biophysics](docs/science/column_biophysics.md) |
| Meteorological forcing | interval-mean-conserving disaggregation; never gap-fills | [forcing](docs/science/forcing.md) |
| Phenology | rate-based signal kernel: cues → flush and shed tendencies | [plant_phenology](docs/science/plant_phenology.md) |
| Carbon allocation | daily FATES PARTEH-H1 priority ladder | [plant_carbon_allocation](docs/science/plant_carbon_allocation.md) |
| Respiration | leaf, stem and fine-root maintenance; growth charged on realized growth | [plant_respiration](docs/science/plant_respiration.md) |
| Soil carbon | ED2's CENTURY network as a carbon matrix ODE | [soil_carbon](docs/science/soil_carbon.md) |
| Demography | cohort and patch fusion/fission, treefall disturbance, wood-density PFT axis | [ed2_comparison](docs/ed2_comparison.md) |
| Time integration | ESDIRK2 by default, adaptive Cash-Karp RK45 as the accuracy baseline | [numerical_scheme](docs/science/numerical_scheme.md) |
| Output | ~208 variables, 7 axes, individually switchable per timescale | [diagnostics](docs/science/diagnostics.md) |

Every process closes a conservation budget each step, asserted in the test suite. The patch axis is
threaded and the output is byte-identical at any thread count.

**Status: v0.2.0, and unbenchmarked.** No EDTS-equivalent regression suite has been run, no site has
been compared flux-for-flux, and no output has been scored against observations. What is verified is
internal: the test suite on two compilers, per-step conservation ledgers, and thread invariance.
Treat the numbers a MEDS run produces as a working model's numbers.
[`CHANGELOG.md`](CHANGELOG.md) states the release's known limitations, and
[`docs/ed2_comparison.md` §0](docs/ed2_comparison.md) says which v0.1.0 numbers are no longer
comparable.

## Quick start

```bash
./scripts/install_netcdf.sh                    # if you do not have netCDF; prints the prefix

cmake -S . -B build -DCMAKE_Fortran_COMPILER=ifx -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
cmake --build build -j
ctest --test-dir build --output-on-failure

./build/meds_main examples/example_demography/example_config_main.toml
python post_proc/plot_site_timeseries.py \
       examples/example_demography/example_output/example_output-D-output.nc -o timeseries.png
```

netCDF is a hard dependency and gfortran, ifx and nvfortran are all supported. Details, including
the parallel builds and three compiler traps that cost real time, are in
[`docs/building.md`](docs/building.md).

## Where to go next

| If you want to | Read |
|---|---|
| run MEDS, coming from ED2 | [`docs/ed2_comparison.md`](docs/ed2_comparison.md) — process by process, keyed to `ED2IN` |
| set up a run | [`docs/configuration.md`](docs/configuration.md) — the two TOML files and the three ways to start |
| see worked examples | [`examples/`](examples/) — demography, biophysics, leaf gas exchange, phenology |
| understand the equations | [`docs/science/`](docs/science/) |
| change the code | [`src/README.md`](src/README.md) — the layout, the placement rules, the library graph |
| drive it from Python | [`python/README.md`](python/README.md) |
| plot the output | [`post_proc/README.md`](post_proc/README.md) |
| know what is deferred | [`docs/ROADMAP.md`](docs/ROADMAP.md) |
| know what changed | [`CHANGELOG.md`](CHANGELOG.md) |

## Design goals

- **Modern Fortran 2018** — modules, derived-type encapsulation, explicit interfaces, `allocatable`
  ownership, parameterized kinds, `pure` and `elemental` kernels, OpenMP `target` array kernels.
- **Modular** — one responsibility per module, no hidden global mutable state, parameters and rates
  passed explicitly as data.
- **Testable** — the kernel libraries build and unit-test standalone, without the demographic state.
- **Reproducible** — no hard-coded model parameters; every run writes back the per-PFT table it
  actually used.

## Scientific reference

- Moorcroft, Hurtt & Pacala (2001), *Ecological Monographs* — the original ED formulation.
- Medvigy et al. (2009), *JGR Biogeosciences* — ED2.
- Longo et al. (2019), *Geoscientific Model Development* 12:4309 — the ED-2.2 technical description.

## License

MEDS is released under the [Apache License 2.0](LICENSE), and the license also covers every earlier
release, v0.1.0 through v0.2.1. [`NOTICE`](NOTICE) carries the copyright notice and the attribution
to ED2, from which parts of MEDS are adapted, and [`AUTHORS`](AUTHORS) lists the copyright holders.
Contributions are welcome under the same license; see [`CONTRIBUTING.md`](CONTRIBUTING.md).
