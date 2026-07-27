# example_biophysics — the fast loop, at hourly resolution

What MEDS does with a meteorological forcing file: solve a coupled canopy energy balance every
30 minutes and hand back leaf, canopy-air and soil temperatures that a met file never contained.

![Hourly canopy energy balance for July of year 50](biophysics_july.png)

Four temperatures over one July at 1 h resolution. Only the first is an input:

| | | |
|---|---|---|
| **Air** | above-canopy air temperature, straight from the ERA5-Land forcing | *boundary condition* |
| **Canopy air space** | the prognostic CAS temperature | *solved* |
| **Leaf, tallest cohort** | leaf temperature of the tallest cohort — the sunlit upper canopy | *solved* |
| **Soil surface** | top soil-layer temperature | *solved* |

The three solved curves separate from the forcing in different directions and with different
phase. Sunlit leaves run above air by day and below it at night — shortwave absorption and
longwave loss against a finite boundary-layer conductance, offset by transpirational cooling. The
canopy air space sits between leaf and soil, ventilated toward the free atmosphere at a rate the
aerodynamic scheme sets. The soil surface is damped and lagged by its heat capacity. Reproducing
that structure from nothing but a met file is the fast loop's whole job, and the third panel —
each store's departure from the driving air temperature — is where it is easiest to read.

## Running it

```bash
cd examples/example_biophysics
./run_example.sh
```

Two stages, both driven by the same recycled year of ERA5-Land forcing for Ithaca NY (42.44 °N,
76.50 °W):

1. **`meds_config_spinup.toml`** — 50 years from bare ground, 2024-07-01 → 2074-07-01. Writes no
   diagnostics at all; its only product is the restart checkpoint `spinup-S-20740701000000.nc`.
   Roughly 15 minutes.
2. **`meds_config_july.toml`** — restarts from that checkpoint and runs July 2074 alone, writing
   the FAST output tier hourly. Seconds.

Then `plot_biophysics.py` builds the figure. `./run_example.sh --replot` skips the model entirely
and rebuilds it from existing output; stage 1 is also skipped automatically whenever its state
file is already present, so iterating on the figure costs seconds rather than the full spin-up.

### Requirements

- A built `meds_main` (`../../build-ifx/meds_main` by default; override with `MEDS_BIN=...`).
- The forcing file `../../data/forcing/ithaca_forcing.nc`. NetCDF files are git-ignored, so it is
  not in the repo — build it with `scripts/download_era5land.py` (needs a CDS API key) followed by
  `scripts/prep_era5land_forcing.py`.
- Python with `numpy`, `netCDF4`, `matplotlib`.

## Why it is split into two stages

Fifty years of hourly output would be ~440 000 records for a figure that needs 744. Stage 1
therefore writes only its final state, and stage 2 restarts into exactly the month being plotted —
which is also why stage 1 *ends* on July 1 rather than January 1.

The restart is exact rather than approximate: the state file carries the fast reservoirs (canopy
air space, soil column, snow, and per-cohort leaf/wood water and temperature), so July 1 continues
the spin-up instead of re-seeding a cold surface and spending the first days of the month relaxing
out of an artificial transient. If you plot the first 48 hours and see no start-up kink, that is
what you are looking at.

## Notes on the configuration

**Forcing recycling.** One calendar year of ERA5-Land drives all 50 years. The recycle window is
*declared*, never inferred:

```toml
recycle       = true
recycle_start = "2024-01-01 01:00:00"    # an EXACT record stamp
recycle_end   = "2025-01-01 01:00:00"    # exclusive; whole calendar years
```

`recycle_start` is `01:00:00`, not `00:00:00`, because `avg_convention = "end"` means each record
is stamped at the *end* of the hour it averages — so the first record of calendar year 2024 is
01:00. MEDS validates this against the file rather than guessing, and rejects a mismatch outright.
The window must also span a whole number of calendar years, so that hour-of-day and day-of-year
survive every wrap. (They do not survive a wrap on any other span, and the failure is quiet: the
daily *mean* shortwave stays correct while the sub-daily phase drifts. See
`docs/dev_plans/MEDS_FORCING_DESIGN.md` §P3.) Model year 2074 is 50 wraps past the file year and
reads the correct hour of the correct day.

**`energy_fluxes = true`** in `[output]` is required — every temperature plotted here belongs to
the `GRP_ENERGY` output group and is silently absent without it.

**Surface temperature** is `soil_temp_top_site`, the top soil layer. MEDS does not currently
expose a separate ground-skin temperature diagnostic, so that is the surface temperature available
and the figure labels it accordingly.

**Tallest cohort.** Cohort composition changes as the stand develops and cohorts fuse and split,
so cohort index 1 is not a stable identity. `plot_biophysics.py` resolves the tallest cohort *per
record* from `height_cohort_fast`, masking the unused slots beyond `n_cohort`.

## Known issue — the soil-surface curve during rain

**The soil-surface trace is not trustworthy during infiltration events, and the spikes in it are a
model bug rather than physics.** Building this example surfaced it; it is not a plotting artifact.

Across July 2074, the 47 hours with appreciable infiltration (Δθ₁ > 0.004) warm the top soil layer
by **+1.23 K/h** on average, against −0.075 K/h in every other hour (correlation of Δθ₁ with
ΔT₁ = +0.55). That sign is wrong twice over: the soil is already ~1.6 K *warmer* than the air
during those events, so rain arriving at roughly air temperature should cool it, and the largest
excursions land at 23:00–01:00 with zero shortwave — one of them after a cloudy day whose peak
shortwave never exceeded 220 W/m².

The layer-wise pattern points at the mechanism. During the five largest events, layer 1 warms while
layers 2–3 cool:

```
dtheta1=+0.048  ->   +3.54  -1.75  -0.37  +0.06  +0.04  +0.00   (K/h, layers 1-6)
dtheta1=+0.040  ->   +4.65  -1.95  -0.79  -0.03  +0.00  +0.02
dtheta1=+0.034  ->   +8.69  -2.05  -3.84  -0.60  -0.01  +0.00
```

Enthalpy is moving *upward* while water percolates *downward*. Two candidates, both in
`soil_energy_step_implicit` (`src/biophysics/meds_soil_energy.f90`):

- the upwind advective-enthalpy convention at lines 64–80 (its comment records that the sign was
  checked once and a suspected bug refuted — worth rechecking against this evidence), and
- `qwf(0) = 0.0_wp` at line 66, which gives the **top** face no water-enthalpy flux at all, so
  infiltrating rain enters layer 1 as mass carrying none of its own enthalpy.

The air, canopy-air-space and leaf curves are unaffected — they are solved by different kernels and
show the expected radiative and transpirational structure throughout.
