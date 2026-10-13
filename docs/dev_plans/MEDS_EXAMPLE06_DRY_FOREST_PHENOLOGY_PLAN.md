# Example 06: dry-forest leaf habits across Central America — PFT and phenology design

Written 2026-10-10. Design only: nothing here is implemented. The owner asked for an example06
that reproduces the regional simulation of Xu et al. (2016, *New Phytol.* 212:80–95) over Central
America. Before any regional run, the PFTs and the phenology have to be designed so the model can
follow how leaves respond to rainfall seasonality. This plan covers that first step in full
(§6–§7). It only outlines the regional run (§8), which gets its own revision once the grid-cell
phase has settled. Decisions D1–D16 are in §9; revisions 2–5 settled D1, D5 and D12–D15. §11 holds the PV tests.

**Revision 1 (owner, 2026-10-10).**
- Initial testing uses ERA5-Land and gridded observations only. No Palo Verde site data are used:
  no station met, no census start, no litter traps. Palo Verde is tested as its ERA5-Land grid
  cell, scored against MODIS (§6.3, §6.4, P2).
- Regional runs drop every non-vegetated cell first, at 0.1° or 0.5°. The regions are then
  designed around what is left, or region mode masks cells itself (§8.1). From bare ground, every
  cell in a box grows a forest — pasture, cropland and town included — so each cell dropped saves a
  full stand's cost.

**Revision 2 (owner, 2026-10-10): the PFT design.**
1. WD and SLA come from Xu et al. (2016).
2. TLP and the plant hydraulics come from Xu's relations as coded in ED2.
3. Height follows Chave et al. (2014)'s pantropical curve, with its climate term.
4. Every other trait is inherited from example05.
5. Phenology: Ds responds only to photoperiod; D, BD and E respond only to leaf water potential.
   Revision 3 replaced this: every PFT responds only to water.

The trait tables are in §6.2; the questions they raise are D12–D15.

**Revision 3 (owner, 2026-10-10).**
- D12: Ds takes WD 0.35 everywhere.
- D14: the first try uses only the water cue, for all four PFTs; the light cue is off.
- D13: Chave's E is computed from ERA5-Land.
- D15: rooting depth deepens with WD (§6.2).

**Revision 4 (owner, 2026-10-10): the PV test.**
- **Starting stand:** the PV test starts from the 2009 inventory of Xu's six plots, the owner's ED2
  files merged into one MEDS census file (§6.4). This relaxes revision 1 for the PV cell only:
  forcing stays ERA5-Land, and scoring stays MODIS.
- **PFT mapping:** ED2 PFT codes map as 29 → Ds, 27 → D, 26 → BD, 24 → E.
- **Soil:** the soil at PV is about 2 m deep.

**Revision 5 (owner, 2026-10-10).**
- **Roots and soil:** the PV soil column is 1.5 m and `root_beta` is 0.001, with the same value for every PFT.
- **Vcmax25 from Xu's A_area–WD relation** (§6.2). Jmax, TPU and Rd follow as ratios, and the
  shade slopes are recomputed.

**Revision 6 (owner, 2026-10-11): settled in beta (#386–#388) and in the example.**
- **Height:** Chave et al. (2014) eq. 6a is a MEDS option (#386). E comes from Chave's own E raster on
  the ERA5-Land grid (`/ibstorage/SharedData/ED_ERA5land_Chave2014_E.nc`), not from ERA5-Land's climate
  (D13): 0.134 at the Palo Verde forcing cell. The owner is developing a height allometry to replace it.
- **Leaf and root turnover (#387):** leaf lifespan 1 yr and fine roots 1.5 per yr for every PFT;
  background turnover runs only while no leaves drop, fine roots age on the same gate and are not shed
  with the leaves (root length sets the root–soil conductance); 30 % of all leaf and root loss is resorbed.
- **Reproduction (#388):** from 0.33 of the local maximum height (11.7 m, about 12 cm dbh here).
- **Soil depth:** 2 m (the design runs settled on 2 m after revision 5's 1.5 m).
- **Open:** in April the nearly bare Ds and D cohorts' leaf potential collapses to -25 to -100 MPa
  (their water stores hit the floor while the soil throttles uptake but not sapflow, #104). This goes to
  the coupled soil–plant water solve (branch `feat/spac-water-solver`), not to a fix here.

**Revision 7 (owner, 2026-10-12): calibrated on the contemporary stand, and the soil map.**
- **Soil:** the ED soil map v2 at the forcing cell: a 1.57 m Cambisol, ~40 % clay, with MEDS's
  pedotransfer function (HT02 + Ottoni) per SoilGrids layer fitted to one column curve (122 mm of
  plant-available water, against 239 mm for the 2 m clay loam); the layers are MEDS's SoilGrids default.
- **Targets** for the 2009 census stand (2009-2010, dry mass): leaf litter ~0.5 kg m⁻² yr⁻¹ (PV traps
  0.54), AGB growth of trees ≥ 10 cm 0.3-0.5 kg m⁻² yr⁻¹, relative AGB growth Ds ~5, D and BD 2-3, E
  1-2 % yr⁻¹; Ds and D fully deciduous, BD dropping 50-75 %, E 0-50 %. The December-January drop is not
  a target.
- **Traits (§6.2):** leaf lifespan 0.5, 0.5, 0.8, 1.5 yr; growth cap 0.25, 0.20, 0.15, 0.10 cm yr⁻¹ at
  1 cm (falling with WD); shed rate 0.15, 0.15, 0.025, 0.01 day⁻¹; `root_beta` 0.001, 0.002, 0.005,
  0.01 (deeper with WD); mortality γ, α, β along WD (Camac's power-law form fitted to example03's
  three BCI classes, so BD and E differ); E keeps Vcmax25 30.2.
- **Mortality window:** `growth_memory_days` 365 (was 90) in every example: with 90 days the growth
  average fell to ~0 after each dry season and the Camac hazard rose from 0.6 to ~10 % yr⁻¹ for E and BD
  in May-June.
- **Litter outputs:** `litter_leaf_site` is leaf litterfall after the #399 fix (it was the labile
  above-ground stream, 0.49 of it).

## 1. What example06 has to show

One set of PFTs and one phenology parameter set, run unchanged in every grid cell, should give:

- deciduous canopies where the dry season is long, evergreen ones where it is short or absent;
- a site's year-to-year swing in deciduousness (dry against wet years);
- leaf habit as something that emerges from hydraulic traits, not a calendar set per site or per PFT.
  In the first try every PFT follows leaf water potential (rev. 3). A photoperiod habit for the stem
  succulents (Borchert & Rivera 2001) is a later option.

That is the claim of Xu et al. (2016): diverse hydraulic traits, together with a phenology driven
by leaf water potential, explain the spatial pattern of deciduousness. MEDS already follows the same
rule: "there is no evergreen flag" (`docs/science/plant_phenology.md` §6).

## 2. The target: Xu et al. (2016)

The main text has been read in full; a copy is at `~/claude_workspace/meds_pheno_proto/xu2016/`. The
Supporting Information was **not** obtained, because Wiley refuses scripted downloads. It holds the
equations (Notes S1), the per-PFT values (Table S3), the species-to-PFT list (Table S2) and the code
(Notes S2).

| | Xu et al. (2016) |
|---|---|
| Model | ED2 plus a new plant-hydraulics module, run per cohort every 10 min. It carries leaf and stem water potential, capacitance, Darcy sap flow, and an S-shaped xylem vulnerability curve with full recovery. Stomata follow Katul-type optimisation: λ rises and Vcmax/Jmax fall as leaf ψ drops. A leafless cohort takes leaf ψ = stem ψ minus gravity. |
| Site | Palo Verde (PV), 10.35° N, 85.35° W; 1,440 mm/yr, dry December–April. Six 0.1-ha plots, started from the 2008 census, spun up to 2009, scored 2009–2013 (2009 and 2012 were dry years). Station forcing; soil 3 m deep, 33 % clay, 39 % sand. |
| PFTs | Four leaf-habit groups (after Borchert et al. 2002):<br>• **Ds**, deciduous stem-succulent: very light wood, sheds early;<br>• **D**, deciduous: leafless about 4 months;<br>• **BD**, brevi-deciduous: leafless 1–2 months;<br>• **E**, evergreen.<br>SLA and wood density (WD) are the PV abundance-weighted species means (Fig. 2b, Table S3). |
| Traits from WD and SLA (Table 1, neotropical fits) | TLP = −4.08 + 0.571 ln SLA − 1.191 ln WD (R² 0.34, RMSE 0.63 MPa)<br>ln Ks,sat = −2.46 WD + 2.35 (R² 0.21)<br>P50 = −3.00 WD − 0.60 (R² 0.25)<br>ln A_area = −1.47 WD + 3.01 (R² 0.21)<br>Also: Vcmax falls with WD; λ rises from Ds to E; roots go deeper with WD at a given size; softwood stores more stem water.<br>Units: SLA in m² kg⁻¹ dry mass, WD in g cm⁻³. Only the growth-respiration fraction was tuned. |
| Phenology | One rule for all PFTs, with no light cue (this is ED2 `phenology = 5`):<br>• shed when the predawn leaf ψ stays below TLP for 10 days in a row;<br>• reflush when it stays above ½ TLP for 10 days;<br>• leaves go 1/20 of the canopy a day in either direction;<br>• `storage_reflush_times` is 2. |
| Site results (new PFTs + hydraulics) | • Leaf litter: 2.9 of the observed 3.1 kg m⁻² over 5 yr (r 0.80).<br>• LAI: 4.14 in the wet season, 0.6 in the dry; the MODIS seasonal range is about 4.<br>• Ds and D shed from December and are nearly leafless by mid dry season.<br>• BD starts in January and keeps 15 %.<br>• E keeps about 50 %, and grows 50 % more in dry years. |
| Regional run | • Domain: 5–25° N, 80–100° W at 0.5°.<br>• Seedlings of every PFT, grown for 150 yr (1860–2010).<br>• Sheffield 0.5° forcing, recycling 1970–2010; each cell's rain scaled to TRMM's 2001–2010 mean.<br>• Historical CO₂; ageing mortality by WD; soil 3 m.<br>• Simulated LAI multiplied by the Hurtt 2005 forest fraction; cells under 25 % forest masked. |
| Regional score (2001–2010, against MODIS MOD15A2 monthly maxima) | • Mean LAI 4.25 against MODIS 3.37 (spatial r 0.61; the original ED2 gave 2.25, r 0.49).<br>• **Minimum leaf cover** (lowest monthly LAI over highest) 0.72 against 0.68 (original ED2 0.86).<br>• Still too evergreen on Mexico's Caribbean coast and in central Nicaragua and Costa Rica. |

The authors list the missing light and photoperiod cues among the model's limits.

## 3. What MEDS has today (beta c20bbd8)

**The kernel already contains Xu's rule, in a smoother form.** The WATER cue sums the predawn ψ
against the PFT's TLP: M_dry += max(0, TLP − ψ_pd) and M_wet += max(0, ψ_pd − TLP). Each sum restarts
when the other side's switch crosses 0.5. TLP is not a phenology key: it is derived from the leaf
pressure–volume curve (`leaf_pi0`, `leaf_elastic_mod`). The same TLP also sets:

- where stomatal stress begins: ½ TLP by default (`pft_stomata_psi_onset`, `meds_config.f90`);
- the low-ψ limb of conductance, which falls from TLP to 2·TLP.

So one trait already drives both the stomata and the leaf habit, as in Xu.

**Example02's Palo Verde fit** (`fitted_parameters.json`, entry `palo_verde`):

- WATER and LIGHT cues on both sides, with `par_min` 5, so the light cue is the photoperiod;
- a sharp flush gate at 12.38 h and short-day shedding below 11.44 h;
- a fitted TLP of −2.09 MPa and a leaf-cover floor (`min_leaf_cover`) of 0.27.

Its water driver is a hypothetical canopy predawn ψ taken from a MEDS run of an *evergreen* stand,
so shedding cannot feed back on ψ. Losses: water alone 0.064, photoperiod alone 0.025, both 0.023.
The reason: litter peaks in December–January, but that ψ only falls below −3 MPa in February–March.

**Not yet done** (`MEDS_PHENOLOGY_SENESCENCE_PLAN.md` §10 and §11 #5): the water cue running on a
model run's own predawn ψ, with the shed → water feedback, and any scored leaf cycle from a coupled
run. That plan names "the Guanacaste example" as the place to do this; example06 is that example.

**Hydraulic traits.**
- Every hydraulic key is an optional per-PFT array: `leaf_pi0`, `leaf_elastic_mod`,
  `leaf_water_sat`, the `wood_*` set, `wood_psi50`, `wood_kexp`, `root_beta`,
  `root_depth_b1`/`b2`, and so on.
- No shipped PFT file sets any of them. In every example, the PFTs differ only in `wood_kmax`, which
  MEDS derives from WD with Xu's neotropical Ks fit (`meds_hydr_lib.f90`).
- `wood_psi50` is not derived; the ROADMAP lists it as "next".
- There is no mortality from hydraulic failure or carbon starvation.

**Leafless cohorts — to verify first.**
- The leaf ψ is worked out from the leaf's water and its leaf carbon (`psi_from_water_content`).
- At zero leaf carbon the relative water content floors at 1e-4, so ψ = π0/1e-4 ≈ −15,000 MPa
  (`meds_water_retention.f90`).
- The snap-to-bare rule can produce exactly zero leaf carbon.
- Xu's model gave a leafless cohort its stem ψ minus gravity.
- This path has not been tested (P0.1).

**Region mode** (`MEDS_POLYGON_RUNTIME_PLAN.md` R2, #289):
- What it does: runs one lat/lon box of `ED_ERA5land` cells, one polygon per cell, with OpenMP
  threads over polygons, and writes NetCDF with a polygon axis.
- What it lacks: per-polygon soil or parameters (one set serves every polygon), any start other
  than bare ground, restarts (R4 is open), a Python entry point (R6), and a mask or stride option.
- Xu's box holds 9,461 land cells at 0.1°, or 392 if only every fifth row and column is kept (0.5°).

## 4. Can one parameter set give different leaf habits from TLP alone? (offline test, 2026-10-10)

This is the owner's question.

**Set-up:**
- The compiled kernel (`meds.plant.pheno`) on example02's Palo Verde drivers, 2008–2013.
- Every parameter shared; only `leaf_psi_tlp` varied, from −0.5 to −3.5 MPa.
- Each dry season (October–September) scored from 2009 to 2013.
- Script, figure and tables are in `~/claude_workspace/meds_pheno_proto/tlp_sweep/`.

| Variant | What changes | Result across TLP −3.5 → −0.5 MPa |
|---|---|---|
| A | example02's PV fit as is | **Shedding starts on the same days (21–25 Nov) at every TLP**: that is the short-day trigger. TLP moves only the refill date (late April at −3.5 MPa, late May to mid-June at −0.5) and the time spent below half cover (74 → 133 days). The lowest cover is 0.38–0.27, held up by the floor: no PFT is evergreen and none goes bare. |
| B | A with the floor at 0 | Same start dates; the lowest cover now runs 0.38 → 0.06. Still no evergreen. |
| C | water cue only on both sides, floor 0, A's water numbers | A spectrum. TLP ≤ −3 stays evergreen (lowest cover 0.99–1.00); −2.5 drops to 0.78; −2.09 to 0.58; −1.5 to 0.40; −0.5 to 0.20. Shedding starts earlier as TLP rises, from April to early March. Nothing goes bare: `shed_rate_max` 0.019 d⁻¹ is too slow for a drought of about 3 months. |
| D | water only, floor 0, Xu-like rates: shed and flush 0.05 d⁻¹, dry sum 2 MPa·d, wet sum 5 MPa·d, smoothing 2 d | Xu's spectrum appears:<br>• −3.5 MPa: evergreen (lowest cover 0.92);<br>• −3.0: brevi-deciduous (lowest 0.21, 34 days below half);<br>• −2.5: deciduous (65 days below half, 22 below 10 %);<br>• −1.0 to −0.5: long-deciduous (102–124 days below half, 63–85 below 10 %).<br>Shedding starts in late March at −3.0 MPa and between 16 January and 14 February at −0.5. Between years, at −3.0 the lowest cover ranges 0.04–0.38. |

**Answer.** With the water cue alone, yes. On a single driver, TLP alone spans the range from evergreen
to four months leafless, with shedding earlier and refilling later as TLP rises. That is how Xu's PFTs
behave. With example02's PV parameter set, no, for three reasons:

- the photoperiod branch of senescence, which takes a maximum with the water branch, fires on the same
  date for every TLP;
- the flush gate holds every PFT to the same lengthening days;
- `min_leaf_cover` sets the depth of shedding per PFT, whatever the traits.

**Caveats** (this is a screen, not a result):

1. **One ψ trajectory serves every TLP.** In a coupled run each PFT has its own predawn ψ, set by its
   rooting depth, conductance, capacitance and stomatal closure, and a bare canopy slows the soil's
   drying. Deep-rooted, dense-wooded E trees should see a higher ψ than this driver, which adds to the
   effect of their low TLP.
2. **The driver dries late.** Monthly mean ψ is −0.62 MPa in January, −2.16 in February and −3.11 in
   March, so even TLP −0.5 starts shedding between mid-January and mid-February, one to two months after the litter peak. ERA5-Land rain at PV
   is 1,430–2,800 mm/yr against the station's ~1,440: a wet bias that delays the dry-down. MODIS at
   the PV cell starts dropping in January (P2). Initial tests take ERA5-Land as it is (rev. 1), so a
   late drop at one cell may be the forcing rather than the scheme. The test of the scheme is the
   pattern across the gradient cells (P3); per-cell rain scaling is a later regional question (D6).
3. **The wet sum counts ψ − TLP in absolute MPa**, so a low-TLP PFT earns wet credit quickly. In
   variant D, the canopy at TLP −3.0 reflushed on the March 2009 rain pulse (ψ back to −1.7 MPa) and
   then shed again. Xu et al. blamed exactly this single-rain reflush for the original ED2's carbon
   starvation. Xu's rule asks for ψ above ½ TLP for 10 days. D4 asks whether MEDS's sums should scale
   with TLP.
4. **The realistic TLP range is narrower than the sweep.** Xu's fit gives −1.1 MPa for WD 0.35 g cm⁻³
   and SLA 20 m² kg⁻¹, and −2.5 MPa for WD 0.8 and SLA 10 (illustrative traits; PV's are in Table S3),
   with an RMSE of 0.63 MPa. On this driver variant D puts the whole −1.1 to −2.5 range in the
   deciduous group, so separating BD and E has to come from each PFT's own ψ (caveat 1).

## 5. Why example02's light cue cannot go regional

The PV fit's light switches are thresholds on day length. Across Xu's domain the astronomical day
length is as follows (`par_min` 5 counts about the same hours):

| Latitude | Shortest / longest day [h] | Flush gate open (> 12.38 h) | Short-day shedding (< 11.44 h) |
|---|---|---|---|
| 5° N | 11.71 / 12.29 | never | never |
| 7° N | 11.59 / 12.41 | 40 d (2 Jun–11 Jul) | never |
| 9° N | 11.47 / 12.53 | 85 d | never |
| 10.35° N (PV) | 11.39 / 12.61 | 100 d (3 May–10 Aug) | 43 d around the December solstice |
| 15° N | 11.11 / 12.89 | 128 d | 100 d |
| 25° N | 10.44 / 13.56 | 152 d | 137 d |

South of about 6.5° N a canopy that sheds could never reflush, though no land cell of Xu's box lies
that far south. The short-day trigger is the real limit: it cannot fire south of 9.6° N, which holds
424 of the box's 9,461 land cells (Panama and southern Costa Rica). North of PV the shedding season
grows with latitude whatever the rain does. For a whole forest the light cue would draw a map of
latitude, not of rainfall. Revision 3 leaves it off for every PFT in the first try (D14). At
PV it also stood in for a water driver that dries too late (§4, caveat 2), and the coupled model may
not need that. A light cue could come back later, either as a flush gate keyed to the *change* in day
length rather than an absolute length (Borchert & Rivera 2001), or with a latitude-relative threshold.
That should happen only if water-only coupled runs show a timing error that the forcing does not
explain.

## 6. Proposed design

### 6.1 Phenology: one parameter set, water only (rev. 3)

- **Cues (rev. 3):** `flush_cue_mask = shed_cue_mask = WATER` for all four PFTs, with every phenology
  number shared. The leaf habit then comes only from each PFT's TLP and its own ψ, as in Xu. A
  photoperiod cue for Ds (Borchert & Rivera 2001) is left for a later try (D14).
- **Rates:** start from variant D (§4): flush and shed at 0.05 d⁻¹ at full tendency, 2-day smoothing,
  and sharp water switches of at least 10 per MPa·d (example02 found softer switches leak and churn
  leaves). MEDS sheds a fraction of the *current* leaves each day, so at 0.05 d⁻¹ the canopy halves in
  about 14 days and reaches 10 % in about 46. ED2's linear 1/20 a day empties it in 20 days.
- **Floor:** `min_leaf_cover` = 0 for every PFT, so the depth of shedding emerges, but only once a
  leafless cohort has a defined ψ (P0.1). Until then, use a tiny shared floor (0.01) with
  `bare_leaf_cover` below it, so snap-to-bare never fires.
- **Background turnover:** decided 2026-10-11: `leaf_lifespan_toc` = 1.0 yr for every PFT, and
  background turnover runs only while no leaves drop (the turnover gate, feat/leaf-root-turnover), so
  the drought drop is the phenology's alone and is not also counted as ageing. 30 % of all leaf and
  fine-root loss is resorbed (`retained_carbon_fraction`, the new default).
- **Fine roots:** decided 2026-10-11: fine roots keep their carbon through a leaf drop, because root
  length sets the root–soil (Gardner) conductance and a drop would cut a leafless tree's uptake. They
  turn over at 1.5 / leaf lifespan = 1.5 yr⁻¹ on the leaves' turnover gate and regrow with the flush.
- **Storage:** carbon storage must pay for a reflush; ED2 set `storage_reflush_times` to 2 for this
  scheme. P2 checks that a D cohort ends the dry season with enough storage to reflush, and that a
  failed reflush shows up in the output (Xu's carbon starvation).

### 6.2 PFTs: Xu's four (rev. 2)

All values come from `~/claude_workspace/meds_pheno_proto/tlp_sweep/pft_traits.py`.

**Identity: WD and SLA.** Digitised from Xu et al. (2016) Fig. 2b (to about ±0.2 m² kg⁻¹ and ±0.02
g cm⁻³); Table S3 holds the exact values.

| | Ds | D | BD | E | MEDS key |
|---|---|---|---|---|---|
| WD [g cm⁻³] | 0.35 (Xu: 0.25; rev. 3) | 0.52 | 0.74 | 0.87 | `wood_density` |
| SLA [m² kg⁻¹ dry mass] | 10.4 | 10.3 | 9.3 | 7.9 | — |
| SLA [m² kgC⁻¹] | 20.8 | 20.6 | 18.6 | 15.8 | `sla` |
| LMA [g m⁻²] | 96 | 97 | 108 | 127 | — |

**Hydraulics: Xu's relations as ED2 codes them** (`ed_params.f90`, `plant_hydro_scheme = 2`).
- ED2 clamps WD to 0.35–0.95 before using it, and so does MEDS for `wood_kmax`. Ds takes WD 0.35
  everywhere: hydraulics, biomass and mortality alike (rev. 3, D12).
- MEDS describes tissue water with pressure–volume curves, not ED2's constant capacitances. The
  curve parameters therefore come from ED2's base relations, which also depend on WD. The leaf π0 is
  solved so that MEDS's TLP (π0 ε / (π0 + ε)) equals the X16 TLP exactly (D2).

| | Ds | D | BD | E | MEDS key |
|---|---|---|---|---|---|
| TLP [MPa] = −4.59 + 0.62 ln SLA_C − 1.15 ln WD | −1.50 | −1.96 | −2.43 | −2.72 | via `leaf_pi0`, `leaf_elastic_mod` |
| leaf ε [MPa] = 2.5 + 37.5 / (1 + e^(5.7 − 8 WD)) | 4.46 | 9.12 | 23.3 | 31.7 | `leaf_elastic_mod` |
| leaf π0 [MPa], solved for the TLP above | −1.12 | −1.62 | −2.20 | −2.50 | `leaf_pi0` |
| stomatal onset = ½ TLP [MPa] | −0.75 | −0.98 | −1.22 | −1.36 | derived (default) |
| P50 [MPa] = −3.0 WD − 0.599 | −1.65 | −2.16 | −2.82 | −3.21 | `wood_psi50` |
| safety margin TLP − P50 [MPa] | 0.15 | 0.20 | 0.39 | 0.49 | — |
| Ks [kg m⁻¹ s⁻¹ MPa⁻¹] = e^(2.348 − 2.455 WD + 0.309) | 6.04 | 3.98 | 2.32 | 1.68 | `wood_kmax` (MEDS derives it already) |
| vulnerability-curve exponent | 4 | 4 | 4 | 4 | `wood_kexp` (default 2) |
| leaf water at saturation [kg H₂O kgC⁻¹] (X16: 1.85 kg kg⁻¹ dry) | 3.70 | 3.70 | 3.70 | 3.70 | `leaf_water_sat` |
| wood π0 [MPa] = 0.52 − 4.16 WD | −0.94 | −1.64 | −2.56 | −3.10 | `wood_pi0` |
| wood ε [MPa] = √(1.02 e^(8.5 WD) − 2.89) | 4.13 | 9.05 | 23.4 | 40.7 | `wood_elastic_mod` |
| wood water at saturation [kg H₂O kgC⁻¹] = 2 (1 − WD/1.53) / WD | 4.41 | 2.54 | 1.40 | 0.99 | `wood_water_sat` |
| resulting wood capacitance in MEDS near saturation ≈ (1 − 0.2) W_sat / (ε + \|π0\|) [kg kgC⁻¹ MPa⁻¹] | 0.70 | 0.19 | 0.043 | 0.018 | derived |
| X16's own wood capacitance (Scholz 2011, ÷ 3), for comparison | 0.70 | 0.32 | 0.083 | 0.038 | — |

**Phenology (rev. 3).** All four PFTs share one set:

| Cue | Ds | D | BD | E |
|---|---|---|---|---|
| flush; senescence | WATER; WATER | WATER; WATER | WATER; WATER | WATER; WATER |

- Dry sum 2 MPa·d and wet sum 5 MPa·d, both with sharpness 10 (variant D, §4).
- Rates 0.05 d⁻¹ in both directions, with 2-day smoothing.
- Floor 0, or 0.01 until P0.1 is done.

The habits differ only through TLP (−1.50 / −1.96 / −2.43 / −2.72 MPa) and each PFT's own ψ.

**Height: Chave et al. (2014) eq. 6a**, one curve for all four PFTs.
- ln H = 0.893 − E + 0.760 ln D − 0.0340 (ln D)², with D in cm and H in m. BIOMASS's `retrieveH`
  multiplies by e^(0.5·0.243²) when back-transforming, and the table below does too.
- The climate stress index is E = 10⁻³ (0.178 TS − 0.938 CWD − 6.61 PS): temperature seasonality
  (WorldClim bio4), climatic water deficit (mm yr⁻¹, negative) and precipitation seasonality (bio15).
- E is computed from ERA5-Land's monthly climate, cell by cell, over 2003–2022 (rev. 3, D13; P0.5):
  - TS = 100 × the standard deviation of the 12 monthly mean temperatures [°C];
  - PS = the coefficient of variation of the 12 monthly rainfall totals [%];
  - CWD = Σ min(0, P_month − ET_month) over the year [mm yr⁻¹, ≤ 0]. ET is a fixed 100 mm a month
    (the Amazon-drought convention), which is an assumption until it is checked against Chave's
    Methods.
- The height cap `hgt_max` stays at 46 m (example05).
- MEDS has a `power` and a `gmm` height form; this is a third (P0.5). It pairs with the Chave (2014)
  biomass equation example05 already uses.

| E | H at D = 10 cm | 30 cm | 60 cm | 100 cm |
|---|---|---|---|---|
| 0.0 | 12.1 m | 22.5 m | 32.0 m | 40.5 m |
| 0.2 | 9.9 m | 18.4 m | 26.2 m | 33.2 m |
| 0.4 | 8.1 m | 15.1 m | 21.4 m | 27.2 m |

**Rooting depth deepens with WD (rev. 3, D15).** MEDS gives each cohort a rooting depth
D = b1 · h^b2, capped at the soil column, with `root_depth_b1` 1.114 and `root_depth_b2` 0.4223 by
default (ED2 IALLOM 1, Christoffersen 2013; 5 m at 35 m height). Within D, the fine roots fall off as
Y(d) = (1 − β^(d/D)) / (1 − β), with β = `root_beta` 0.1, so the whole profile stretches with D.
Neither ED2 version (local, or the owner's fork) makes depth depend on WD. Xu's per-PFT depths are in
Table S3. The proposed form is:

  D(h, WD) = 1.114 · (WD / 0.6) · h^0.4223, capped at the soil column.

- Depth is proportional to WD at a given height. The proportionality is an assumption, not a fit;
  evergreen, dense-wooded dry-forest trees draw deeper water in isotope studies.
- A WD-0.6 tree keeps MEDS's default curve; 0.6 is example04's WD and the Camac reference.
- It needs no code: the per-PFT `root_depth_b1` key exists.

| | Ds | D | BD | E | MEDS key |
|---|---|---|---|---|---|
| b1 [m] | 0.650 | 0.966 | 1.374 | 1.615 | `root_depth_b1` |
| b2 | 0.4223 | 0.4223 | 0.4223 | 0.4223 | `root_depth_b2` |
| D at h = 2 m | 0.87 m | 1.29 m | 1.84 m | 2.16 m | |
| D at h = 5 m | 1.28 m | 1.91 m | 2.71 m | 3.19 m | |
| D at h = 10 m | 1.72 m | 2.55 m | 3.63 m | 4.27 m | |
| D at h = 20 m | 2.30 m | 3.42 m | 4.87 m | 5.72 m | |
| D at h = 30 m | 2.73 m | 4.06 m | 5.78 m | 6.79 m | |
| height reaching 3 m / 6 m | 37 m / never | 15 m / 76 m | 6 m / 33 m | 4 m / 22 m | |

On Chave's curve with E = 0.2, trees of 10, 30 and 60 cm DBH stand about 10, 18 and 26 m tall.
- **In Xu's 3 m soil,** every D, BD and E canopy tree hits the floor of the column. Their profiles
  are then identical, and the WD effect vanishes above about 15 m.
- **In a 6 m column** (example04's), it survives to canopy height. Hence D11's recommendation.

**Photosynthesis (rev. 5): Vcmax from Xu's A_area.**
- **A_area:** exp(3.01 − 1.47 WD) × exp(0.5 · 0.38²), from Table 1 (neotropical), back-transformed
  with the log-regression correction.
- **Vcmax25:** the value at which MEDS's own leaf kernel reproduces A_area for a sunlit,
  well-watered leaf, using example05's leaf (Medlyn g1 2.82, peaked temperature response).
  Conditions: PAR 2,000 µmol m⁻² s⁻¹, 25 °C, 1 kPa VPD, CO₂ 385 ppm. A 30 °C leaf changes it by
  less than 1 %.
- **Jmax, TPU and Rd** keep example05's ratios.
- **The Vcmax and Rd shade slopes** come from Ma et al. (2025): −(0.00242 V + 0.06212) and
  −(0.18974 · 0.015 V + 0.11744).

| | Ds | D | BD | E | MEDS key |
|---|---|---|---|---|---|
| A_area [µmol m⁻² s⁻¹] | 13.04 | 10.15 | 7.35 | 6.07 | — |
| Vcmax25 [µmol m⁻² s⁻¹] | 65.64 | 50.85 | 36.62 | 30.18 | `vcmax25` |
| Jmax25 = 1.6983 Vcmax25 | 111.5 | 86.4 | 62.2 | 51.3 | `jmax_vcmax_ratio` |
| TPU25 = 0.167 Vcmax25 | 10.96 | 8.49 | 6.12 | 5.04 | `tpu_vcmax_ratio` |
| Rd25 = 0.015 Vcmax25 | 0.985 | 0.763 | 0.549 | 0.453 | `rd_vcmax_ratio` |
| shade slope of Vcmax [per unit LAI above] | −0.2210 | −0.1852 | −0.1507 | −0.1352 | `[derived].kplastic_vcmax` |
| shade slope of Rd | −0.3042 | −0.2622 | −0.2217 | −0.2033 | `[derived].kplastic_rd` |

The script is `ex06_pv_test/vcmax_from_aarea.py`. Example05's 31.33 maps back to its slopes
(−0.1379, −0.2066).

**Everything else: example05's values.** Example05 varies a few traits by its three WD classes
(early 0.36, mid 0.50, late 0.68). The new PFTs take the nearest class: Ds → early, D → mid, BD and
E → late.

| Inherited by nearest WD class | Ds (early) | D (mid) | BD (late) | E (late) |
|---|---|---|---|---|
| `leaf_lifespan_toc` [yr] (rev. 7; 1.0 for all on 2026-10-11) | 0.5 | 0.5 | 0.8 | 1.5 |
| `fineroot_turnover_rate` [yr⁻¹] | 1.5 | 1.5 | 1.5 | 1.5 |
| `max_absolute_growth_rate` [cm yr⁻¹ at 1 cm] (rev. 7; was 0.422, 0.306, 0.118, 0.118) | 0.25 | 0.20 | 0.15 | 0.10 |
| `shed_rate_max` [day⁻¹] (rev. 7; was 0.15, 0.15, 0.05, 0.02) | 0.15 | 0.15 | 0.025 | 0.01 |
| `root_beta` (rev. 7; 0.001 for all on 2026-10-10) | 0.001 | 0.002 | 0.005 | 0.01 |
| `mort_gamma` [yr⁻¹] (census γ − 0.014 treefall) | 0.0108 | 0.0074 | 0.0047 | 0.0037 |
| `mort_alpha` [yr⁻¹] | 0.196 | 0.143 | 0.108 | 0.095 |
| `mort_beta` [yr cm⁻¹] | 12.3 | 28.5 | 60.2 | 84.9 |

The mortality rows are rev. 7: each is Camac et al. (2018)'s form p₀ (WD/0.6)^b fitted (ln p on ln WD)
to example03's BCI fits for its three classes (WD 0.36, 0.50, 0.68) and read at these WDs: census γ
0.02027 (WD/0.6)^−0.372, α 0.1276 (WD/0.6)^−0.793, β 38.58 (WD/0.6)^2.123. They replace the nearest-class
pins (early, mid, late, late), which gave BD and E the same mortality.

The same for all four, as in example05:
- Medlyn g1 2.82 with its g0, and the Jmax, TPU and Rd ratios to Vcmax;
- root:leaf ratio 1.0;
- sapwood stem respiration 31.25, root respiration 2.7, reproduction 0.10;
- the growth-cap exponent 0.5, the SLA and leaf-lifespan shade slopes, optics, allocation, and Chave's biomass
  equation (eq. 4).

**Where this departs from Xu (2016)**, kept visible:
- Vcmax falls with WD as in Xu, but A_area is turned into Vcmax with MEDS's own leaf model rather
  than ED2's (rev. 5).
- Rooting depth deepens with WD as in Xu, but in an assumed form (proportional to WD); Xu's values
  are in Table S3.
- Ds has WD 0.35 rather than Xu's 0.25 (rev. 3).
- Stomata are Medlyn with a TLP-derived onset, not Katul with λ rising from Ds to E.

### 6.3 Forcing and soil

- **Forcing:** ERA5-Land from the archive for every test, site cells and region alike (rev. 1). It
  is the only regional source MEDS reads, so the grid-cell tests run on exactly the forcing the
  region will use. Runs cover the archive's years from July 2002 on, recycled for spin-up.
- **Rain bias (D6), later:** Xu scaled each cell's rain to TRMM's 2001–2010 mean. For the regional
  stage the choice is ERA5-Land as is, or a per-cell annual rain factor from a gauge–satellite
  climatology (CHIRPS or IMERG). The second needs either a forcing preprocessing step or a
  per-polygon rain multiplier; MEDS has neither yet. Not part of the initial tests.
- **Soil at the PV cell (rev. 4):** 2 m deep. In 2 m of soil the WD-dependent rooting depths of §6.2
  reach the bottom early: E at 1.7 m height, BD at 2.4 m, D at 5.6 m and Ds at 14 m. So at PV they
  separate only small trees, and BD and E canopy trees have to differ through TLP, P50 and
  capacitance.
- **Soil for the region (D11):** one column shared by every cell (region mode has one soil). The
  candidates are PV's 2 m, Xu's 3 m, or example04's 6 m column, which keeps the rooting depths apart
  for canopy trees. Its texture is a choice, not a site measurement.

### 6.4 Starting state

Region mode starts from bare ground, as Xu did. Site tests that iterate on the scheme need short runs,
so there are two kinds of run (D7):

- **Design runs at the PV cell (P2; rev. 4):** start from the PV 2009 inventory, with demography off
  (`demography.demography_on = false`), a few minutes each.
  - **The file:** `~/claude_workspace/PV_demography/pv_2009_census.csv`, written by
    `merge_ed2_pss_css.py` from the owner's ED2 files `PV_2009_2P_site<k>.lat10.4lon-85.4.{pss,css}`.
  - **Layout:** one site holding the six 0.1-ha plots, two equal patches per plot, so 12 patches. The
    file has 336 cohorts, carries ED2's `n` as `nplant`, and maps the ED2 PFT codes as
    29 → Ds, 27 → D, 26 → BD, 24 → E.
  - **Not carried:** the plots' stand ages (7–60 yr) and ED2 soil pools have no census column, so
    they are recorded as comments in the file. MEDS census patches start at age 0, and soil carbon
    comes from the `[init]` keys.
  - **The stand** (stems from 1 cm dbh, per ha):

    | | Ds | D | BD | E | total |
    |---|---|---|---|---|---|
    | stems | 495 | 1,928 | 3,967 | 298 | 6,688 |
    | stems ≥ 10 cm | 12 | 135 | 307 | 7 | 461 |
    | basal area [m²] | 0.99 | 5.69 | 13.30 | 1.80 | 21.8 |

- **Design runs at the gradient cells (P3):** the same fixed-stand approach, the same at every cell so
  that only the climate differs between cells. The stand is the PV inventory or a synthetic stand
  (D7).
- **Emergence runs (P4):** bare ground, 150 years of recycled forcing, all four PFTs seeding, to see
  which habit wins along the gradient.

## 7. Phases

**P0 — prerequisites (small code changes).**
1. **Leafless ψ:** write a unit test for a cohort snapped to zero leaf carbon. If its ψ is undefined,
   give a leafless cohort its wood ψ minus its gravity head (Xu's rule), defined once in the fast
   prepass rather than as a special case in phenology.
2. **One coupled drought-deciduous PFT, end to end:** a PFT file with a `[phenology]` section using
   WATER masks, and a 2-year run at the PV grid cell that writes per-cohort `dmax_psi_leaf_cohort`,
   `leaf_flush_tendency_cohort`, `leaf_shed_tendency_cohort`, leaf carbon and storage. Confirm the
   shed → transpiration → soil-water feedback and the one-day lag.
3. **(If D4 = relative)** change the kernel so the water sums scale with TLP, and add
   `test_plant_phenology` cases.
4. **(If D3 = derive)** derive `wood_psi50` from WD.
5. **Chave (2014) height (rev. 2):** a third `[allometry].height_allometry` form, `"chave2014"`, with
   the stress index as one config number (`chave_e`) for site and grid-cell runs. Its inverse (height
   → dbh) is a quadratic in ln D. Unit-test it against BIOMASS's `retrieveH`. A script computes E
   on the forcing grid from the ERA5-Land archive's monthly temperature and rainfall (§6.2). For
   region runs, the per-cell E travels with §8.1's mask file.

**P1 — the PFT file.**
- `pft_parameters.toml` for Ds, D, BD and E, from §6.2's tables.
- A table of each PFT's resulting TLP, P50, Ks, Vcmax and rooting depth.
- A one-day check of each PFT's stomatal response to ψ with the leaf kernel (example01's tool).

**P2 — the Palo Verde grid cell, coupled (design runs).**
- **Runs:** the ERA5-Land cell at 10.35° N, 85.35° W, from the PV 2009 inventory on 2 m of soil.
  The run starts in January 2008, one year for soil water (with the 2009 stand), and is scored
  2009–2013 (rev. 4).
- **Target — MODIS LAI at the cell** (2009–2013 monthly maxima, as in example02's `drivers/`):
  - LAI falls from 4.7 in December to 2.8 in January;
  - it bottoms at 1.5 in March–April;
  - it refills to 4.5 in May and peaks at 5.3 in October;
  - mean LAI 3.8; minimum leaf cover 0.23–0.30 a year (0.28 in the climatology).
- **Scored against:**
  - that cycle: the minimum leaf cover, the months of the drop, the minimum and the refill;
  - dry against wet years;
  - as a qualitative check only, each PFT's dry-season leaf cover against Xu's published behaviour:
    Ds and D bare by mid dry season, BD about 15 %, E about 50 %.
- **Questions answered:**
  - Does water alone, on each PFT's own ψ, follow the MODIS drop in January and the refill in May?
  - Does the canopy floor come out near 0.28 without being set?
  - Do the shared parameters need refitting? Only the common ones (sums and rates) are candidates.

**P3 — the rainfall-seasonality gradient (design runs).** The same PFT file, parameters and synthetic
stand at 5–6 ERA5-Land cells, chosen to spread across dry-season length and latitude (D8). Each cell
is named after a well-known site only to locate it; no site data are used. The rainfall figures
below are approximate and are to be read from ERA5-Land itself:

| Site | Lat, lon | Rain, dry season | Role |
|---|---|---|---|
| La Selva, Costa Rica | 10.43, −84.00 | ~4,000 mm, none | should stay evergreen |
| Palo Verde, Costa Rica | 10.35, −85.35 | ~1,440 mm, 5 months | the P2 site |
| Santa Rosa, Costa Rica | 10.84, −85.62 | ~1,700 mm, 5–6 months | a near-duplicate check |
| BCI, Panama | 9.15, −79.85 (0.15° east of Xu's box) | ~2,600 mm, 4 months | mostly evergreen; light-driven leaf exchange is outside this scheme |
| Kaxil Kiuic, Yucatán | 20.08, −89.55 | ~1,100 mm, 6 months | similar seasonality, 10° further north |
| Mexican Caribbean coast | ~19.6, −88.0 | ~1,300 mm | where Xu's run was too evergreen |

- **Scoring:** each site against MODIS LAI (monthly maxima from ORNL point subsets) on four measures —
  mean LAI, minimum leaf cover (Xu's measure), the month of the LAI minimum, and the month of green-up.
- **Pass:** the ranking of minimum leaf cover follows dry-season length, with no site-specific
  parameter.
- **Guard:** if any shared number is refitted, refit it leaving one site out at a time, so nothing is
  tuned to PV.

**P4 — emergence at the P3 cells (bare ground, 150 years).**
- Which PFTs dominate where.
- The stand's minimum leaf cover against MODIS.
- Whether the deciduous PFTs starve of carbon in long dry seasons.
- How many years the minimum leaf cover takes to settle (Xu used 150; fewer may do).
- The cost per cell, which sets the regional design (§8). This phase is the regional run in miniature.

**P5 — write-up into example06.** A README ordered by the pipeline, and one figure. The regional half
follows §8.

## 8. The regional run (outline; its own revision later)

**Domain and cost.**
- Xu's domain is 5–25° N, 80–100° W: 9,461 land cells at 0.1°, or 392 at 0.5° (one cell in 25).
- A bare-ground cell costs about 8 s per simulated year on one thread (the R2 demo at Ithaca).
- A closed tropical stand costs far more. Example05's BCI run (420 yr, up to 25 patches) took 121 min
  on 20 threads: about 6 min per simulated year on one thread.
- At that cost, 392 cells take about 1 hour per simulated year on 40 threads (about 6 days for 150
  years), and the full 0.1° box is out of reach.
- Two things cut this: dropping the non-vegetated cells (§8.1), and fewer patches per cell (for
  example a maximum of 6–10). P4 measures the real per-cell cost before the grid is chosen (D9).

### 8.1 Drop the non-vegetated cells first (rev. 1)

From bare ground MEDS grows a forest in every cell it runs: it knows nothing of land use, so
pasture, cropland, town and bare rock all become forest stands, each as costly as a real one. Such
cells are dropped before anything is designed, at whichever resolution is chosen.

1. **A vegetation mask on the forcing grid.** Build a 0/1 raster on the ED_ERA5land 0.1° grid: 1
   where the woody (forest) cover passes a threshold. One script, run once, kept with the example.
   Candidate sources (D9):
   - ERA5-Land's own high-vegetation cover `cvh`, a time-invariant ERA5-Land field. It keeps every
     input in ERA5, but it is not on disk and needs a one-time fetch.
   - MODIS land cover (MCD12Q1) or tree cover (MOD44B).
   - LUH2 forest fraction, the successor of the Hurtt et al. data Xu used (25 % threshold).
2. **At 0.5°:** keep one 0.1° cell per 0.5° block, for example the vegetated cell nearest the
   block's centre, and only where the block's vegetated share passes the threshold. The same 0/1
   raster then describes both grids, so no separate stride option is needed.
3. **Then run only the masked cells**, in one of two ways (D9):
   - **(a) Design regions around the mask:** cut boxes and tiles so that they hold as few
     non-vegetated cells as possible, and run each as a region. This works with today's region mode
     but still pays for the non-vegetated cells inside each box.
   - **(b) Grid-level masking in region mode (recommended).** Add a `[region].mask_file` key: a
     raster on the forcing grid, where only cells marked 1 become polygons. It is one more selection
     rule beside `land_fraction_min`. `MEDS_POLYGON_RUNTIME_PLAN.md` §5 already plans "a mask or ID
     raster later", the polygon id stays the cell's global index (OR6), and the output keeps its
     `polygon` axis. It also serves tiling: a job array can split one mask into tiles.

The same mask scores the run against MODIS, so the model and the observations cover the same
cells. The cost of the run scales with the number of masked cells, to be counted once the source
and threshold are chosen.

### 8.2 Restarts, observations and output

**Restarts.** Region mode has no restarts (R4). A 150-year spin-up must finish in one job; otherwise
R4 has to come first, or the masked cells are split into tiles run as a job array.

**Observations.** MODIS LAI (C6.1) for 2003–2012, since the ERA5-Land archive starts in June 2002,
taken as monthly maxima. A regional extract needs NASA Earthdata access (D9).

**Output.** Monthly `lai_site`, `lai_pft`, `litter_leaf_site`, `gpp_site` and `et_site`.
`psi_leaf_site` is a time mean, not predawn; a predawn site variable would help diagnosis.

## 9. Decisions for the owner

- **D1 PFTs. Settled (rev. 2):** Xu's four (Ds, D, BD, E), with the WD and SLA of Xu's Fig. 2b.
  Replace the digitised values with Table S3's if the owner has them.
- **D2 TLP in MEDS keys. Proposed (rev. 2):** the leaf ε from ED2's WD relation, with π0 solved so
  MEDS's TLP equals the X16 TLP (§6.2). The alternative is one shared ε (12 MPa), which gives
  π0 −1.33 / −1.69 / −2.02 / −2.22. ε also sets leaf capacitance.
- **D3 P50 and wood water storage.** Recommended: the values the ED2 relations give (§6.2), written
  into example06's file, which leaves the other examples bit-identical. Deriving them from WD in
  code, like `wood_kmax`, is the cleaner long-run design. But it changes example04 and example05
  results, so it belongs in its own pull request.
- **D4 Water sums.** Keep absolute MPa·d sums (current), or scale them by |TLP| so a threshold means
  the same strain for every PFT, as Xu's ½ TLP does. Recommended: prototype the scaled version
  offline in P0.3, and adopt it only if the low-TLP churn of §4 caveat 3 also appears in P2.
- **D5 Stomata and Vcmax. Settled (rev. 2; Vcmax rev. 5):** Medlyn stomata with the TLP-derived
  onset, inherited from example05. Vcmax25 comes from Xu's A_area–WD fit, 65.6 / 50.9 / 36.6 / 30.2
  (§6.2). The traits example05 varies by WD class map by the nearest class.
- **D6 Rain (regional stage).** ERA5-Land as is for every initial test (rev. 1). For the regional
  run: keep it as is, or scale each cell's rain to a gauge–satellite climatology (as Xu did with
  TRMM). Decide after P3, if the gradient cells show a timing error that follows ERA5-Land's rain.
- **D7 Starting state for the design runs.** The PV cell starts from the PV 2009 inventory (rev. 4).
  The gradient cells (P3) are open: the PV inventory at every cell, or a synthetic four-PFT stand.
- **D8 Gradient cells.** The P3 list, or the owner's own list.
- **D9 Regional grid and vegetation mask (§8.1).**
  - The mask's source and threshold: ERA5-Land `cvh`, MODIS land or tree cover, or LUH2.
  - The grid: 0.5° (one cell per block, Xu's resolution) or every vegetated 0.1° cell.
  - Regions designed around the mask, or `[region].mask_file` in region mode (recommended).
  - Restarts (R4) first or not, and the MODIS regional extract.
- **D10 Scope.** Recommended: one example06 in two parts, grid-cell design then regional run. The
  alternative is a split: example06 the gradient cells, example07 the region.
- **D11 Soil column.** PV: 2 m (rev. 4). For the gradient cells and the region, one column for every
  cell. A deep column (example04's 6 m) keeps WD-dependent rooting depths apart for canopy trees;
  PV's 2 m and Xu's 3 m collapse them (§6.3). Texture: example04's BCI clay or a generic medium loam.
  Spatially varying soil is a later plan (runtime plan OR2).
- **D12 Ds's wood density. Settled (rev. 3):** 0.35 everywhere (ED2's floor). Xu's 0.25 is not used.
- **D13 Chave's E. Settled (rev. 6):** Chave's own E raster on the ERA5-Land grid. Computing E from
  ERA5-Land (rev. 3) gave -0.16 at Palo Verde against the raster's 0.137: Chave's TS is WorldClim bio4 as
  stored (x1000) and his CWD does not use a fixed 100 mm/month ET. For site and
  grid-cell runs E is one config number. Region runs need E per cell, the first spatially varying
  parameter, carried as a second field of §8.1's mask raster.
- **D14 Ds's photoperiod. Settled (rev. 3):** off for the first try; every PFT follows water only.
  If it comes back later: example02's short-day trigger (11.44 h) cannot fire south of 9.6° N
  (424 of 9,461 land cells), so a latitude-relative day length would be needed.
- **D15 Rooting depth. Settled (rev. 3):** it deepens with WD, as D = 1.114 (WD/0.6) h^0.4223, through
  the per-PFT `root_depth_b1` (§6.2). The proportionality is assumed; Table S3's values would replace
  it.

- **D16 The PV drought is too weak (§11).** What should bring the trees' water stress into
  December–January:
  - rooting depth or root profile (`root_depth_b1`, `root_beta`);
  - soil storage (depth, texture);
  - the water cue's signal (TLP vs ½ TLP, predawn vs daily mean ψ);
  - photoperiod for Ds.

## 10. Data the plan needs

- ERA5-Land: the archive at `/ibstorage/SharedData/ED_ERA5land` (July 2002 on), the only forcing.
- MODIS LAI: monthly maxima at the PV cell (already in example02's `drivers/`), point subsets for the
  P3 cells (ORNL API), and later a regional extract.
- A vegetation mask source (D9). ERA5-Land `cvh` is not on disk; MODIS land cover and LUH2 need
  downloads.
- PFT traits: WD and SLA digitised from Xu et al. (2016) Fig. 2b. Table S3 (exact values, rooting
  depths) would replace them if available; Wiley blocks scripted downloads.
- ERA5-Land monthly temperature and rainfall for Chave's E (D13): from the archive, no download.
- The PV 2009 inventory (rev. 4): the owner's ED2 .pss/.css files, merged in
  `~/claude_workspace/PV_demography/`.
- Not used in the initial tests: PV station met and litter traps.

## 11. First PV test (2026-10-10)

**Set-up.**
- **Code fix first (P0.1):** a leafless cohort's leaf node now reads the wood potential less the
  gravity head (`leaf_node_psi`, in `meds_plant_hydraulics`). Before the fix, an empty leaf store was
  read as π0 × 10⁴. The fix has a unit test in `test_plant_hydraulics`; ctest passes 67/67.
- **Run:** the PV 2009 inventory at the ERA5-Land cell (7.8 km away), July 2008 – December 2010.
  - Four PFTs from §6.2, all on the water cue (variant D).
  - Example05's height curve (Chave 2014 was not yet a MEDS option; it is now, #386).
  - 2 m of clay loam (Carsel & Parrish 1988), starting wet.
  - Demography restructuring off, so growth and mortality only.
- **Cost:** 68 s on 12 threads. Energy and water closed with 0 failures in about 1.05 million
  checks. Files: `~/claude_workspace/meds_runs/ex06_pv_test/`.

**Survival: yes.** Every one of the 336 cohorts is alive at the end. Plants fell 4.5–14 % over 2.5
years (Camac hazard 1.9–6.2 % yr⁻¹, Ds the highest), and stems of 10 cm and more rose for Ds and D.
Storage per plant fell 5–25 %, most of it in April 2009, when stomata closed (β 0.44–0.93) while the
leaves kept respiring.

**Seasonality: no.** Leaf cover stayed at 0.99–1.15 for every PFT in both dry seasons. Only Ds dipped,
to 0.89, in May 2009. MODIS at the cell falls to 0.24–0.26 of its peak, with the drop in January.

**Why: the 2 m soil carries the stand through the dry season.**
- Only the top 5 cm dries hard (−5.6 MPa in March 2009).
- Below 10 cm the soil stays at −0.2 to −1.3 MPa until the May rains. Predawn leaf ψ follows the
  wettest rooted layer, so its monthly mean reaches only −1.3 MPa (April 2009) and −0.9 MPa
  (April 2010).
- No PFT spends long below its TLP: Ds 21 days, D 8, BD 3, E 0 in 2009.
- ERA5-Land is not too wet here: 49 mm fell from December 2008 to April 2009, and 1,593 mm in all of
  2009.
- Dry-season ET stays at 80–97 mm a month from a full canopy.
- A run with example05's BCI clay (saturation 0.60) is wetter still: predawn ψ −0.8 MPa in April.

The water cue works as designed, but the simulated drought reaches the trees about three months
after MODIS shows the canopy thinning. Open choices for the owner (§9 D16):
- shallower effective rooting, or more roots near the surface;
- a soil that stores less water;
- the cue on a less conservative water signal;
- photoperiod for Ds after all.

**Fine roots and soil depth (2026-10-10).**
- **Stand fine-root profile** (the 2009 inventory, MEDS's own rules; `ex06_pv_test/root_profile.py`):
  - 2.46 MgC ha⁻¹ of fine roots;
  - 12 % of fine roots in the top 10 cm, 33 % in 30 cm, 49 % in 50 cm, 76 % in 1 m; half are above
    51 cm and 95 % above 167 cm.
  - The cause: 94 % of fine-root carbon is in cohorts rooted to the 2 m floor, and `root_beta` 0.1
    spreads each profile over its whole depth.
- **Jackson et al. (1996)**, tropical deciduous forest (β 0.961): 33 / 70 / 86 / 98 % at the same
  depths, half above 17 cm.
- **The model roots are about three times too deep**, which is why the predawn ψ follows the wet deep
  layers. With the same rooting depths, `root_beta` 3.8 × 10⁻⁴ reproduces Jackson's curve: half above
  17 cm, 95 % above 75 cm, 70 % in the top 30 cm.
- **Soil depth at PV: no measurement exists for the Powers plots.** A literature and gridded-data
  search (sources in `ex06_pv_test/pvsoil_refs/`) found:
  - **Xu et al. (2016):** "depth to bedrock can reach several meters", with no source given; their
    3 m was a "best estimate".
  - **Terrain:** the plots sit on flat colluvium or alluvium below karstic limestone hills (Becknell
    et al. 2021; Powers et al. 2009). The nearby hillslopes are 30–60 cm deep, about half under 30 cm
    (Rojas-Zeledón 2018).
  - **SoilGrids 2017 at the cell:** depth to an R horizon 1.4 m, a 60 % chance of rock within 2 m;
    depth to bedrock including weathered rock 8 m.
  - **SoilGrids v2.0:** texture clay to clay loam (38–43 % clay) to 2 m, with stones.
  - **Analogue:** at Horizontes, 45 km away on deep soils, deciduous trees root to about 1.4 m and
    evergreen trees to about 2.2 m (Smith-Martin et al. 2020).
  - **Estimate:** 1–3 m, about 2 m most likely; 1 and 3 m are the sensitivity bounds.

**Rerun on 1.5 m of soil with `root_beta` 0.001 (owner, 2026-10-10).** Half of the fine roots now lie
above 15 cm and 95 % above 65 cm, close to Jackson's tropical deciduous profile.

The first rerun exposed a second bug, which does not depend on the leafless fix.
- **What happens:** a big E tree overdraws its leaf water in one fast step, down to the apoplastic
  share. Its leaf then reads π0 × 10⁴ MPa and never refills, although the wood is at −1.3 MPa.
- **Cause:** for a vulnerability exponent other than 1 or 2 (X16 uses 4), the Kirchhoff table
  extrapolated its last slope past r = 8. The one-rule quadrature also missed the integral at large r.
  Together they turned the potential's sign, so the edge conductance went negative and was floored.
- **Fix** (`meds_hydr_lib`):
  - one 7-point panel per table interval, so the table rises monotonically;
  - panels in `kirchhoff_integral`;
  - an analytic tail (`kirchhoff_tail`) past r = 8.
  Runs with exponent 2 are unchanged.
- **Tests:** a new far-tail test; ctest 67/67. The desiccated leaf now refills in one step.
- **Left open:** the overdraw itself (no limit on transpiration against the stored leaf water within
  a step) remains, but no longer traps a cohort.

**Results with the fix** (40 s run):
- **Survival:** every cohort is alive. Plants fell 21.6 % for Ds (hazard 9.7 % yr⁻¹), 9.0 % for D,
  8.5 % for BD and 3.3 % for E.
- **Storage** falls 65–90 % each March–May and recovers only partly. Ds runs out in 2010 and regrows
  only 56 % of its canopy.
- **Lowest leaf cover**, 2009 / 2010:

  | PFT | Ds | D | BD | E |
  |---|---|---|---|---|
  | lowest leaf cover | 0.34 / 0.21 | 0.48 / 0.49 | 0.44 / 0.47 | 0.71 / 0.41 |

- **No PFT goes below 10 %.** Shedding starts in March and bottoms out in mid-April to early May.
- **Stand LAI** dips to about 2.3, against MODIS's 1.3–1.9 from February to April, with the drop
  beginning in January.
- **Water status:** predawn ψ reaches −2.0 to −2.5 MPa in March and April, and at least −7.4 MPa in
  one cohort. Only the top layer dries by January–February.

So the habits now order roughly as in Xu (Ds the most deciduous, E the least in 2009), but the water
cue still acts two months after MODIS. The January drop precedes any soil drying in ERA5-Land, as
example02 found.

**Vcmax from Xu's A_area (owner, 2026-10-10).**
- **Method:** A_area = exp(3.01 − 1.47 WD) × exp(0.5 · 0.38²) gives 13.0 / 10.2 / 7.35 / 6.07
  µmol m⁻² s⁻¹. Vcmax25 is solved with MEDS's own leaf kernel (Medlyn, peaked temperature
  response; PAR 2,000, 25 °C, 1 kPa VPD, CO₂ 385 ppm).
- **Values:** Vcmax25 65.6 / 50.9 / 36.6 / 30.2 µmol m⁻² s⁻¹. Jmax, TPU and Rd follow their ratios
  (1.6983, 0.167, 0.015), and the Ma et al. shade slopes are recomputed per PFT.
- **Run** (`ex06_pv_test/output_fast`, report `report_vcxu.py`):
  - Stand monthly LAI falls from 5.0 in February to 2.8 (May 2009) and 2.1 (April 2010), relative
    0.57 / 0.42.
  - Monthly minima by PFT: Ds 0.36 / 0.11, D 0.48 / 0.35, BD 0.61 / 0.47, E 0.70 / 0.54.
  - Timing is still March–May. Daily-minimum ψ stays at −0.4 to −1.1 MPa through January and
    reaches −1.6 to −2.5 MPa in February.
  - Survival: Ds 0.76, D 0.90, BD 0.91, E 0.97. Stand NPP 13.8 MgC ha⁻¹ yr⁻¹ (2009).
  - The one-step leaf-water overdraw still drops 4–90 cohort-hours a month to the floor in BD and E
    (March–April).

**Water-stress consolidation (2026-10-10, merged into beta as #384, 28f550a; CHANGELOG).**
- **The floor events were a rain-pulse runaway**, not a numerical overdraw. The stomata saw only the
  previous day's predawn potential, and the soil-side wilting ramp and −1.5 MPa cap made the wood pay
  for root uptake the soil refused.
- **The change:**
  - The capacity limb on the current leaf potential (TLP → 2·TLP) is the always-on fast control.
  - The predawn linear gs decline is removed.
  - The soil has no wilting point: the root sink is the plant's request, limited only by the water a
    layer holds.
- **Run** (`ex06_pv_test/output_new`; `output_soilonly` isolates the soil part):
  - Canopy-tree floor events are gone: 249 cohort-hours before; 25 hours now, all in near-leafless Ds
    saplings.
  - Transpiration 592 → 517 mm/yr on the old leaf-kernel diagnostic. #384 also makes `transp_cohort`
    report the water the leaves lose: 401 mm/yr for this run (1003 mm over 914 days).
  - The roots now dry the shallow root zone, so predawn ψ crosses TLP sooner and the water cue sheds
    harder. Stand LAI minimum 2.83 → 1.26 (2009) and 2.10 → 1.14 (2010). Wet-season LAI recovers only
    to ~3.5, with E at 0.3 of its leaves.
  - This comes almost entirely from the soil change. The phenology (`shed_water_sum` 2 MPa·day) and
    rooting (`root_beta` 0.001, 1.5 m soil) were tuned while the cap held the root zone at −1.5 MPa,
    and need revisiting.

**The adopted design-run set (2026-10-10; `examples/example06_regional_run/`).**
- **Carbon reserve:** `storage_cushion` counts leaf and fine root and is 3 (#385). At 1 x leaf, the large
  trees of every PFT could not pay for their next flush and stayed leafless through the wet season.
- **Soil:** 2 m. The four habits separate only around 2-2.5 m.
  - At 1.5 m every PFT reaches -3.3 to -3.8 MPa and sheds 79-97 % of its leaves.
  - At 3 m only Ds sheds.
- **Water cue:** `shed_rate_max` 0.15 / 0.15 / 0.05 / 0.02 per day. At one shared 0.05, D never spends
  enough days below its TLP to go bare.
- **Stomata:** g1 3.49 / 3.04 / 2.44 / 2.09 kPa^0.5 (Lin et al. 2015's wood-density slope, anchored to
  example 04's 2.82 at 0.60); `wstress_sref_stomata` 0.4 per MPa (Sabot et al. 2022's fits for
  tropical trees); `stomata_psi_onset` 0. g0 stays example 04's 0.00046 mol m-2 s-1.
  - Raising g0 to the measured g_min (0.005, Duursma et al. 2019) drove BD and E to -9 and -7 MPa in
    April.
  - The leaves dry out in days, but the cue sheds at most 5 % (BD) or 2 % (E) a day, and nothing drops a
    leaf that has failed hydraulically. A realistic g0 needs leaf shedding that responds to leaf water.
- **Result:** lowest cover Ds 0.07 / 0.08, D 0.07 / 0.10, BD 0.46 / 0.63, E 0.91 / 0.96 (2009 / 2010);
  stand minimum LAI 3.3 against MODIS 1.5.

## References

- Xu, X., Medvigy, D., Powers, J.S., Becknell, J.M. & Guan, K. (2016) Diversity in plant hydraulic
  traits explains seasonal and inter-annual variations of vegetation dynamics in seasonally dry
  tropical forests. *New Phytol.* 212: 80–95.
- Borchert, R. & Rivera, G. (2001) *Tree Physiol.* 21: 213–221.
- Bartlett, M.K., Scoffoni, C. & Sack, L. (2012) *Ecol. Lett.* 15: 393–405.
- Markesteijn, L. & Poorter, L. (2009) *J. Ecol.* 97: 311–325.
- Chave, J. et al. (2014) Improved allometric models to estimate the aboveground biomass of tropical
  trees. *Glob. Change Biol.* 20: 3177–3190 (eq. 6a, height; eq. 4, biomass).
- Scholz, F.G. et al. (2011) in *Size- and Age-Related Changes in Tree Structure and Function*
  (stem capacitance against WD), as cited in ED2's `ed_params.f90`.
- ED2 `ED/src/dynamics/phenology_driv.f90` (`phenology = 5`) and `ED/src/init/ed_params.f90`
  (`plant_hydro_scheme = 2`, the X16 trait relations).
- `MEDS_PHENOLOGY_SENESCENCE_PLAN.md` §10–§13; `MEDS_POLYGON_RUNTIME_PLAN.md`;
  `examples/example02_canopy_phenology`.
