# Example 02 — Canopy phenology

The leaf-phenology module of MEDS on its own, at three forests and four leaf habits: a temperate
deciduous broadleaf (Harvard Forest), an evergreen conifer (Hyytiälä), and two tropical species
under one climate (Barro Colorado Island), one drought-deciduous and one exchanging its leaves in
the bright dry season. All four run the same kernel
([`meds_phenology.f90`](../../src/slow_dynamics/plant/meds_phenology.f90)); they differ only in
parameter values. A Python script drives the compiled kernel through
[`meds.plant.pheno`](../../python/meds/plant/pheno.py), and `Phenology.step` applies the coupled
model's leaf rule (`leaf_turnover_step`) with carbon never limiting the flush. The equations are in
[`docs/science/plant_phenology.md`](../../docs/science/plant_phenology.md).

## The model in brief

Each day the kernel advances two tendencies in [0, 1], a flush tendency and a senescence (shed)
tendency, from up to four cues, chosen for each side by a bit mask:

| Cue (bit) | Driver | Flush switch | Senescence switch |
|---|---|---|---|
| TEMP (1) | daily air temperature | warmth sum above 5 °C since midwinter, past a centre | cold sum below a base since midsummer, past a centre |
| DAYLENGTH (2) | day length (a calendar) | days longer than a threshold | days shorter than a threshold |
| WATER (4) | predawn leaf water potential | wet sum above a threshold potential | dry sum below it |
| PAR (8) | running-mean PAR at the cohort's top | bright enough | bright (or dim) enough |

Every switch is a logistic σ(s (x − x*)) with a centre x* and a sharpness s. The flush signal is the
product of the flush switches; the senescence signal is the larger of the seasonal trigger (the
product of its TEMP, DAYLENGTH and PAR switches) and the water trigger. Each signal is smoothed over
a few days. The leaf rule then loses leaves to senescence at `shed_rate_max` × the shed tendency,
down to `min_leaf_cover` of the full canopy, plus a background turnover while the canopy flushes,
and grows leaves at up to `flush_rate_max` × the flush tendency. **An evergreen is a PFT whose
senescence stops at `min_leaf_cover`**, not a flag.

## Sites and data

| | Harvard Forest, Massachusetts | Hyytiälä, Finland | Barro Colorado Island, Panama |
|---|---|---|---|
| Forest | deciduous broadleaf (red oak, red maple), 42.5° N | Scots pine, 61.8° N | moist tropical forest, 9.2° N |
| Drivers | air temperature (HF001), day length; 2002–2023 | air temperature and PAR (ICOS FI-Hyy); 2018 – July 2024 | PAR and soil water content (BCI tower, Detto); July 2012 – August 2017 |
| Observations | MODIS LAI; leaf fall on tagged trees (HF003); broadleaf litter baskets (HF069) | needle litter traps, about monthly (ICOS) | fine litter traps, about monthly (GLiMP control plots, Gigante, next to BCI) |
| Fit targets | MODIS LAI scaled each year between its winter and summer levels; fraction of leaves fallen; cumulative basket fraction within each leaf year; **one canopy of litter a year** | needle-fall rate of each collection over its mean; **0.30 canopy of needles a year** | litter rate of each collection over its mean; **one canopy of litter a year** |

The annual amounts are not in the timing data. A deciduous canopy is built once a year and every leaf
falls. Scots pine in southern Finland keeps 3.4–4.2 needle cohorts and needles live about three years
(Pensa & Jalkanen 1999, *Silva Fennica* 33:654), so about 0.3 of the canopy falls each year. The
leaves falling at BCI each year have about the canopy's area, 7.3 m² per m² of ground (Leigh 1999,
ORNL NPP data set BRR).

**The BCI water driver is a surrogate.** Predawn leaf water potential is not measured, so the
tower's soil water content stands in for it, through a Campbell retention curve
ψ = −e<sup>−3.74</sup> θ<sub>g</sub><sup>−2.58</sup> MPa fitted to the 1,020 paired samples of
soil water content and potential that Kupers et al. (2019) took in the BCI 50-ha plot. The tower's
volumetric water content becomes gravimetric by its ratio to the plot's on their four sampling
dates (0.83, a bulk density). The surrogate runs from −0.05 MPa in the wet season to −0.4 to −0.9 MPa
at the end of a dry season, a root-zone soil potential, less negative than a leaf's turgor-loss point;
the species' water thresholds are set on its scale. BCI's PAR is the tower shortwave times its
PAR/shortwave ratio (2.11).

The fits minimise the summed squared errors by differential evolution (best of three seeds). At BCI
only the light exchanger is fitted; the drought-deciduous species is set by hand (below), because
no observation separates its leaves from the community's.

## Run

Build the Python library once (as in [example 01](../example01_leaf_gas_exchange/README.md#run)),
then from the repository root:

```bash
python examples/example02_canopy_phenology/fetch_phenology_data.py          # once: data/, ~1 min
PYTHONPATH=python python examples/example02_canopy_phenology/run_phenology.py
```

The run takes a few seconds: it reads [`fitted_parameters.json`](fitted_parameters.json), prints the
parameters and scores, and writes `harvard_forest.png`, `hyytiala.png` and `bci.png`. `--fit` refits
the fitted parameters first and rewrites the JSON (about 4 minutes on 40 cores; `--workers` sets the
processes).

## Which light cue: day length or PAR

Day length is a calendar: the same every year and for every cohort. PAR is the light a cohort
actually receives, with its clouds, its year-to-year swings and, in the coupled model, its place in
the canopy. Each site was fitted with each cue on each side (same loss, same data; the PAR fits
carry one more parameter, the averaging window):

| Light cue, flush / senescence | Harvard Forest loss | Hyytiälä loss | BCI exchanger loss |
|---|---|---|---|
| day length / day length | **0.043** (spring RMSE 5.5 d) | 1.03 (needle-fall r 0.85) | 0.101 |
| PAR / PAR | 0.052 (spring RMSE 9.6 d) | **0.70** (r 0.91) | **0.065** |
| day length / PAR | 0.049 | 0.75 | — |
| PAR / day length | 0.046 | 0.81 | — |
| none / day length, none / PAR | — | — | 0.138 / 0.067 |

At BCI this comparison used the light exchanger alone on the 17 years of GLiMP collections
(2003–2019), with the Lutz tower's shortwave (STRI) as the PAR source and its canopy kept above 0.8.

No single light variable serves all three: day length times the Harvard spring (a PAR gate nearly
doubles its error), PAR the pine's needle fall and the tropical leaf exchange, where day length can
only repeat the same calendar each year (its year-to-year correlation of dry-season litter is 0.08
against 0.59 for PAR). A mixed pair never beat the best pure one. So MEDS keeps day length and PAR as
two separate cues, and each habit below uses the one its site supports.

## Four habits, one model

★ marks the fitted values (in [`fitted_parameters.json`](fitted_parameters.json)); the rest are set
or left at their defaults (5 °C warmth base, sharpness 0.04 per K day for warmth and 0.1 for cold,
−1 per hour for the senescence day length, 5-day smoothing). "—" means the cue is off.

| Parameter | Harvard Forest | Hyytiälä | BCI drought-deciduous | BCI light exchanger |
|---|---|---|---|---|
| flush cues, senescence cues | TEMP + DAYLENGTH, TEMP + DAYLENGTH | TEMP + PAR, TEMP + PAR | WATER, PAR + WATER | WATER, PAR + WATER |
| `flush_degree_days` [K day] | 94 ★ | 88 ★ | — | — |
| `shed_base_temp`, `shed_degree_days` | 17.0 °C, 46 ★ | 26.5 °C, 456 ★ | — | — |
| `flush_daylength_threshold`, sharpness | 13.75 h, 8.0 h⁻¹ ★ | — | — | — |
| `shed_daylength_threshold` | 9.13 h ★ | — | — | — |
| `flush_par_threshold`, sharpness | — | 307, 0.37 ★ | — | — |
| `shed_par_threshold`, sharpness [µmol m⁻² s⁻¹] | — | 124, −0.013 ★ (dim light) | 5000 (never) | 445, +0.50 ★ (bright light) |
| `par_window` [day] | — | 16 ★ | — | 3.7 ★ |
| water threshold `leaf_psi_tlp` [MPa] | — | — | −0.33 | −1.5 (never reached) |
| water sums: senescence / flush [MPa day] | — | — | 1 / 3 | 1 / 3 |
| `flush_rate_max`, `shed_rate_max` [day⁻¹] | 0.034, 0.33 ★ | 0.29, 0.031 ★ | 0.067, 0.10 | 0.17, 0.0029 ★ |
| `leaf_turnover_rate` [yr⁻¹] | 0.0003 ★ | 0.27 ★ | 0 | 0.61 ★ |
| `min_leaf_cover` | 0 | 0.90 ★ | 0 | 0.83 ★ |

What the four parameter sets produce (median year):

| | Harvard Forest | Hyytiälä | BCI drought-deciduous | BCI light exchanger |
|---|---|---|---|---|
| What sets the flush | day length in 15 of 21 years, warmth in the rest | warmth (met 16 May), then PAR ends it | the first weeks of rain | always on |
| Flushing (tendency above 0.5) | 4 May – 19 Aug | 19 May – 10 Sep | after the rains (May–June) | all year |
| Senescence | cold sum met 21 Sep; half fallen 17 Oct, all by 20 Nov | from early September as PAR dims; over by 9 Oct | as the soil dries: leaf drop begins between 21 Feb (2013) and 25 Apr (2017) | whenever the 4-day PAR mean passes 445: about 25 days a month in January–April, under 8 from May to November |
| Leaf cover | 0 to 1 | 0.90 to 1 | 0 to 1; below half for 29–100 days each dry season | 1 |
| Leaf loss per year | 1.01 canopy | 0.30 (0.22 senescence, 0.08 turnover) | 1.16 | 1.00, half of it in January–April |
| Leaf lifespan (mean cover / loss) | 0.41 yr | 3.1 yr | 0.68 yr | 1.0 yr |

At BCI the two species see the same soil and the same light. The deciduous one sheds at −0.33 MPa,
which the soil passes in every dry season, and ignores light; the exchanger's −1.5 MPa is never
reached, but bright light triggers its senescence while its fast flush keeps the canopy full, so it
renews its leaves in the dry season instead of losing them.

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
| RMSE of the needle-fall rate over its mean | 0.85 |
| Correlation, model against traps | 0.90 |
| Needle fall | 0.30 canopies a year |
| Lowest leaf cover | 0.90 |

![Barro Colorado Island](bci.png)

**BCI.** (a) The drivers: the soil water potential surrogate (the deciduous species' −0.33 MPa
threshold dotted) and the light exchanger's running-mean PAR (its threshold dotted). (b) Leaf cover
and senescence tendency of both species. (c) The light exchanger's litter rate per collection over
its mean, traps against model. (d) The same, by month.

| Score, light exchanger (collections 2013 – August 2017) | Value |
|---|---|
| RMSE of the litter rate over its mean | 0.19 |
| Correlation, model against traps | 0.91 |
| Leaf litter | 1.00 canopy a year |

- **The flush gate is a photoperiod step at Harvard** (13.75 h, 8 per hour). A gradual gate is still
  partly open in October, when the warmth sum is high: the canopy then refills while it senesces and
  drops about 1.9 canopies of leaves a year, which no timing observation shows. Holding the annual
  litter to one canopy is what pins the gate.
- **Autumn is timed to within 3 days, spring to within 6**, but neither follows the year-to-year
  variation closely (r 0.44 and 0.34). At Harvard the 13.75 h gate sets the spring date in most
  years, and the 17 °C senescence base makes the cold sum a clock that starts in late summer.
- **At Hyytiälä the cold sum is only a clock.** Its fitted base, 26.5 °C, is above almost every daily
  mean, so it counts days from midsummer; dimming PAR times the needle fall. The model puts too much
  needle fall in July–August and misses the small winter fall, and its flush (to mid-September)
  overlaps the start of needle fall, so some needles are shed and regrown then.
- **The BCI litter traps catch all fine litter**, leaves of every species plus twigs, flowers and
  fruit, so the light exchanger stands for the community. It reproduces the dry-season high and the
  wet-season low; it misses the February peak and puts too much in May.
- **The drought-deciduous species is a demonstration, not a fit.** Its threshold and sums give one
  leafless period per dry season, longest in the dry 2013 (100 days) and shortest in the mild 2017
  (29 days). A fast reflush needs sharp water switches: a sharpness of 5 per MPa day leaves a 0.7 %
  senescence floor that, refilled by the flush, sheds a quarter canopy a year in the wet season.

## Files

- [`fetch_phenology_data.py`](fetch_phenology_data.py) — downloads the site data into `data/`
  (not tracked) and cites the sources.
- [`run_phenology.py`](run_phenology.py) — the site drivers and observations, the fit and the figures.
- [`fitted_parameters.json`](fitted_parameters.json) — the parameters of the four habits.
- `harvard_forest.png`, `hyytiala.png`, `bci.png` — the figures.
