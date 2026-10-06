# Example 02 — Canopy phenology

The leaf-phenology module of MEDS on its own, at a deciduous and an evergreen forest. Both run the
same kernel ([`meds_phenology.f90`](../../src/slow_dynamics/plant/meds_phenology.f90)) on the same two
cues, daily air temperature and day length, with no water stress; they differ only in parameter
values. A Python script drives the compiled kernel through
[`meds.plant.pheno`](../../python/meds/plant/pheno.py), and `Phenology.step` applies the coupled
model's leaf rule (`leaf_turnover_step`) with carbon never limiting the flush. The equations are in
[`docs/science/plant_phenology.md`](../../docs/science/plant_phenology.md).

## The model in brief

Each day the kernel advances two tendencies in [0, 1]:

| | Flush tendency | Senescence (shed) tendency |
|---|---|---|
| Temperature | warmth sum above 5 °C since midwinter, past a centre | cold sum below a base since midsummer, past a centre |
| Day length | longer than a threshold | shorter than a threshold |
| Combined | product: warm enough *and* long enough days | product: cold enough *and* short enough days |

Every switch is a logistic σ(s (x − x*)) with a centre x* and a sharpness s, and each signal is
smoothed over 5 days. The leaf rule then loses leaves to senescence at `shed_rate_max` × the shed
tendency, down to `min_leaf_cover` of the full canopy, plus a background turnover while the canopy
flushes, and grows leaves at up to `flush_rate_max` × the flush tendency. **The pine is evergreen
because its senescence stops at `min_leaf_cover`**, not because of a flag.

## Sites and data

| | Harvard Forest, Massachusetts | Hyytiälä, Finland |
|---|---|---|
| Forest | deciduous broadleaf (red oak, red maple), 42.5° N | Scots pine, 61.8° N |
| Driver | daily air temperature, HF001, 2002–2023 (2002 spins up) | daily air temperature, ICOS FI-Hyy, 2018 – July 2024 |
| Observations | MODIS LAI (MCD15A3H, 3 × 3 pixels at the tower); leaf fall on tagged trees (HF003); broadleaf litter baskets (HF069) | needle litter traps, about monthly (ICOS ancillary data) |
| Fit targets | MODIS LAI scaled each year between its winter and summer levels; fraction of leaves fallen; cumulative basket fraction within each leaf year; **one canopy of litter a year** | needle-fall rate of each collection over its mean; **0.30 canopy of needles a year** |

The annual amounts are not in the timing data. A deciduous canopy is built once a year and every leaf
falls. Scots pine in southern Finland keeps 3.4–4.2 needle cohorts and needles live about three years
(Pensa & Jalkanen 1999, *Silva Fennica* 33:654), so about 0.3 of the canopy falls each year.

Ten parameters are fitted at Hyytiälä and nine at Harvard, which has no `min_leaf_cover`. The fit
minimises the summed squared errors by differential evolution (best of three seeds). The flush
day-length threshold is kept at least an hour below the longest day, so the gate opens fully in
summer.

## Run

Build the Python library once (as in [example 01](../example01_leaf_gas_exchange/README.md#run)),
then from the repository root:

```bash
python examples/example02_canopy_phenology/fetch_phenology_data.py          # once: data/, < 1 min
PYTHONPATH=python python examples/example02_canopy_phenology/run_phenology.py
```

The run takes a few seconds: it reads [`fitted_parameters.json`](fitted_parameters.json), prints the
parameters and scores, and writes `harvard_forest.png` and `hyytiala.png`. `--fit` refits both sites
first and rewrites the JSON (about 4 minutes on 40 cores; `--workers` sets the processes).

## One model, two parameter sets

Both sites run the same equations on the same cues. Every parameter is listed below; ★ marks the
fitted ones (the values are in [`fitted_parameters.json`](fitted_parameters.json)), and the rest
are shared.

| Process | Parameter | Harvard Forest | Hyytiälä | What the difference does |
|---|---|---|---|---|
| Cues | `flush_cue_mask`, `shed_cue_mask` | 3, 3 | 3, 3 | — (TEMP + LIGHT on both sides) |
| | `light_variable` | 1 | 1 | — (day length) |
| | `flush_cue_timescale`, `shed_cue_timescale` [day] | 5, 5 | 5, 5 | — |
| Flush: warmth | `flush_base_temp` | 5 °C | 5 °C | — |
| | `flush_degree_days` ★ [K day] | 94 | 138 | the pine needs more warmth, met in late May rather than late April |
| | `flush_temp_sharpness` [(K day)⁻¹] | 0.04 | 0.04 | — |
| Flush: day length | `flush_light_threshold` ★ [h] | 13.75 | 13.08 | gate open 28 Apr – 14 Aug at Harvard, 2 Apr – 9 Sep at Hyytiälä |
| | `flush_light_sharpness` ★ [h⁻¹] | 8.0 | 6.6 | — (both near a step: the gate, not warmth, ends the flush in late summer) |
| Senescence: cold | `shed_base_temp` ★ | 17.0 °C | 6.8 °C | every day below 17 °C counts at Harvard, so the clock starts in late summer; only days below 7 °C count at Hyytiälä |
| | `shed_degree_days` ★ [K day] | 46 | 37 | — |
| | `shed_temp_sharpness` [(K day)⁻¹] | 0.1 | 0.1 | — |
| Senescence: day length | `shed_light_threshold` ★ [h] | 9.13 | 14.49 | Harvard's days stay above it until late November, so it only scales the autumn rate (switch 0.07–0.24 through October); at Hyytiälä it switches on when days fall below 14.5 h, on 27 Aug |
| | `shed_light_sharpness` [h⁻¹] | −1 | −1 | — |
| Rates | `flush_rate_max` ★ [day⁻¹] | 0.034 | 0.024 | a bare broadleaf canopy fills in about a month; the pine only regrows the needles it lost |
| | `shed_rate_max` ★ [day⁻¹] | 0.333 | 0.159 | the broadleaf canopy drops within weeks even at a switch of 0.15 |
| | `leaf_turnover_rate` ★ [yr⁻¹] | 0.0003 | 0.24 | none in the broadleaf season; the pine sheds needles all summer |
| Leaf cover | `min_leaf_cover` ★ | 0 | 0.88 | **the leaf habit**: the broadleaf senesces to bare, the pine stops at 88 % |
| | `bare_leaf_cover` | 0.02 | 0.02 | — |

What the two parameter sets produce (median year):

| | Harvard Forest | Hyytiälä |
|---|---|---|
| Warmth requirement met | 25 Apr (25 Mar – 4 May) | 23 May (19 May – 3 Jun) |
| What sets the spring flush | day length in 15 of 21 years | warmth in every year |
| Flushing (tendency above 0.5) | 4 May – 19 Aug | 29 May – 13 Sep |
| Senescence | cold sum met 21 Sep; leaves half fallen 17 Oct, all by 20 Nov | from late August; needle fall over by 10 Oct |
| Leaf cover over the year | 0 to 1 | 0.88 to 1 |
| Leaf loss per year | 1.01 canopy, in autumn | 0.30 canopy: 0.23 senescence, 0.07 summer turnover |
| Leaf lifespan (mean cover / annual loss) | 0.41 yr | 3.0 yr |

The leaf lifespan is an output: 5 months against 3 years, from the same model with ten numbers
changed.

## Results

![Harvard Forest](harvard_forest.png)

**Harvard Forest.** (a) Model leaf cover against MODIS LAI, with the two tendencies, 2014–2018.
(b) The day the canopy is half green (MODIS) and half fallen (HF003), year by year. (c) The fraction
of leaves fallen in autumn, each year of the model against the tagged trees.

| Score (2003–2023) | Value |
|---|---|
| RMSE: MODIS LAI / leaf fall / baskets | 0.184 / 0.063 / 0.076 |
| Half-green day: r, RMSE | 0.34, 5.5 days |
| Half-fallen day: r, RMSE | 0.44, 2.9 days |
| Leaf litter | 1.01 canopies a year |

![Hyytiälä](hyytiala.png)

**Hyytiälä.** (a) Model leaf cover and the two tendencies. (b) Needle-fall rate of each collection
over its mean, traps against model. (c) The same, averaged by the month a collection ended.

| Score (collections from July 2018) | Value |
|---|---|
| RMSE of the needle-fall rate over its mean | 1.02 |
| Correlation, model against traps | 0.85 |
| Needle fall | 0.30 canopies a year |
| Lowest leaf cover | 0.88 |

- **The flush gate is a photoperiod step at both sites** (13.75 and 13.1 h, sharpness 8.0 and 6.6 per
  hour). A gradual gate is still partly open in October, when the warmth sum is high: the canopy then
  refills while it senesces and drops about 1.9 canopies of leaves a year, which no timing
  observation shows. Holding the annual litter to one canopy is what pins the gate.
- **Autumn is timed to within 3 days, spring to within 6**, but neither follows the year-to-year
  variation closely (r 0.44 and 0.34). At Harvard the 13.75 h gate sets the spring date in most
  years, and the 17 °C senescence base makes the cold sum a clock that starts in late summer.
- **The Hyytiälä needle fall peaks in September–October as observed**; the model puts too much in
  July–August and misses the small winter fall. The largest residuals are timing errors of about two
  weeks: in 2022 the model's needle fall peaked in the second half of September, the traps' in early
  October.
- **The pine's flush and needle fall overlap** from late August to mid-September, so about 0.12
  canopy is shed and regrown then. That is how the fit reconciles 0.30 canopy of needle fall a year
  with a canopy that only dips to 0.88; the floor is not otherwise constrained by the traps.

## Files

- [`fetch_phenology_data.py`](fetch_phenology_data.py) — downloads the site data into `data/`
  (not tracked) and cites the sources.
- [`run_phenology.py`](run_phenology.py) — the site drivers and observations, the fit and the figures.
- [`fitted_parameters.json`](fitted_parameters.json) — the fitted parameters of both sites.
- `harvard_forest.png`, `hyytiala.png` — the figures.
