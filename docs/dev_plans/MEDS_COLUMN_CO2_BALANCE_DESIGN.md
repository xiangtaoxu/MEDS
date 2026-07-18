# MEDS Column CO2 Balance — Design Document

*Definitive implementation plan for the canopy-air-space CO2 balance: a prognostic canopy-air CO2
mixing ratio (the third CAS twin), advanced from cohort carbon fluxes (GPP − autotrophic
respiration) and soil heterotrophic respiration, with a closed CO2 budget and a turbulent flux to
the free atmosphere.*

Design-only (no source committed by this document). Synthesizes three design angles — ED2-faithful
fast-loop port, minimalist CAS third-twin, and biogeochemistry column-carbon — against research on
ED2 (`rk4_derivs`, `canopy_struct_dynamics`, `soil_respiration`), CLM5/CTSM-FATES, and CliMA Land,
verified against the live MEDS tree and hardened by a four-lens adversarial review. Companion to
`MEDS_ENERGY_BALANCE_DESIGN.md` (the CAS energy/water twins this mirrors) and
`MEDS_COLUMN_HYDROLOGY_DESIGN.md`.

---

## 0. Executive summary

**What / where / why (one paragraph).** This module adds a **prognostic canopy-air-space CO2 mixing
ratio** — `can_co2` in `[umol/mol]` (= ppm), the **third CAS twin** beside `can_enthalpy` and
`can_shv` — advanced each fast sub-daily substep from (i) the net biotic source `Reco − GPP` (cohort
autotrophic respiration + soil heterotrophic respiration − gross primary production, all reduced to
`[umol CO2 / m2 ground / s]`) and (ii) turbulent exchange with the free atmosphere, closing a
machine-precision CO2 budget with `NEE` / `NEP` / `loss2atm` diagnostics. The **module** — all
stateless `pure` kernels (the CAS CO2 twin update, the cohort→ground flux aggregation, the MVP
heterotrophic-respiration flux, the column assembler), the slow soil-carbon pool, the `co2_opts_t`
config, and the unit tests — lives in **`src/biogeochemistry/`** (`libmeds_biogeochemistry`, links
`src/shared` only), honoring the stated placement and owning the whole ground-to-atmosphere carbon
path. The **only** biophysics change is two one-line field additions (`can_co2` to `cas_state_t`,
`co2_atm` to `cas_atm_forcing_t`), because `can_co2` is *physically* the third CAS twin and must ride
the P3 patch-state lockstep beside its enthalpy/humidity siblings. **Why build it at all:** ED2 — the
model MEDS reimplements — carries a prognostic CAS CO2. It buys a closed CO2 budget (a `resid ~ 0`
like every other biophysics kernel), qualitative nocturnal sub-canopy CO2 accumulation, and the
`ca = can_co2` feedback into photosynthesis that CLM5 / FATES / CliMA deliberately omit by pinning
leaf CO2 to the atmospheric forcing.

**The four decisive calls:**

1. **Placement — module in `biogeochemistry`, state field in `biophysics`.** Column carbon belongs in
   biogeochemistry; the CO2 twin is *physically* a CAS object. Because the kernels are **stateless
   bare-scalar** (the shipped `canopy_air_update` idiom), the *code* home is decoupled from the
   *state* home: kernels go to biogeochemistry (need only `meds_shared`); the prognostic `can_co2`
   field rides in `cas_state_t`. No `biogeochem → biophysics` library edge; the driver (`meds_aux`,
   which links both) passes `cas%can_co2` by reference. **This requires formally widening the
   biogeochemistry charter** to own *both* fast carbon exchange and slow carbon pools (§2.4).
2. **Units — dry-air molar mixing ratio `[umol/mol]` with a MOLAR capacity.** Faithful to ED2; a
   mixing ratio is a conserved *intensive* quantity under T/p change, exactly parallel to the
   mass-specific enthalpy `[J/kg]` and humidity `[kg/kg]` twins. Rejected: CLM's partial-pressure
   `[Pa]` (not conserved under T/p change).
3. **Heterotrophic respiration — MVP Q10 (or ED2 capped-exponential) × moisture modifier on a
   possibly-stubbed soil-C pool.** Exploits ED2's own fast/slow seam: within a day the soil-carbon
   pools are frozen, so the fast-loop variability is *purely* the temperature/moisture modifier. Full
   CENTURY multi-pool dynamics defer to a slow biogeochem phase.
4. **Coupling — implicit-in-atmosphere update now, `ca = can_co2` feedback + variable-density
   `denseffect` at P3.** The fast source (`Reco − GPP`) is explicit/lagged; only the
   atmospheric-exchange term is implicit (L-stable), byte-for-byte mirroring `canopy_air_update`'s
   enthalpy branch.

**Two forward-designs baked in (added after the first review — §3.5 and §4c′).** (i) **Multi-layer
canopy air space:** the scalar CAS box is the `n_can_layer = 1` degenerate case of a 1-D vertical
K-theory transport column that reuses the soil column's Thomas machinery (`meds_soil_solver` promoted to
`meds_shared`) — proven bit-identical at `n = 1`, with the shared canopy grid + eddy-diffusivity `K(z)`
intended for all three CAS twins (CO2 pioneers it). (ii) **DAMM heterotrophic respiration** as a
selectable `hr_model ∈ {HR_Q10, HR_EXP_ED2, HR_DAMM}` — Davidson et al. 2012's dual-Arrhenius/
Michaelis-Menten, whose mechanistic moisture response *removes* the empirical `theta_dry` floor problem,
reframed as one instance of a generalizable soil-gas-flux seam so future CH4/N2O stay additive. Both are
**P2 targets**; P0/P1 ship the scalar box + empirical Q10 unchanged.

**Honesty scope of the MVP (do not over-claim).** The MVP is a *structurally complete, dimensionally
closed* CO2 twin — not yet a validated NEE product. Three known biases keep it out of observation
comparison until P3: it omits committed **growth and storage** respiration (§4b), it drops ED2's
Monin-Obukhov stability coefficient `c3` from the exchange (unity for all three twins — damps
nocturnal accumulation vs ED2; §1.5), and its moisture modifier floors on the van Genuchten residual
`theta_res` rather than ED2's air-dry `soilcp` (§4c). Each is documented, tested where testable, and
carries a P2/P3 upgrade slot.

---

## 1. Scientific background

### 1.1 The well-mixed canopy-air-space CO2 budget

MEDS already treats the canopy air space (CAS) as a single well-mixed box of depth `can_depth [m]` per
unit ground area, advancing two prognostic *specific* (per-kg-air) twins — enthalpy and humidity — in
`canopy_air_update` (`meds_column_energy.f90:253`). CO2 is the natural third twin. Store it as a
**dry-air molar mixing ratio** `can_co2 [umol/mol]` (identical to ppm by volume). Its budget over the
box is

```
d(can_co2)/dt = ( f_bio + gatm_co2 * (co2_atm − can_co2) ) / ccapcan          [umol/mol/s]
```

with

- **molar CAS capacity** `ccapcan = can_dmol * can_depth  [mol_dryair/m2]` — the CO2 analog of the
  energy/water twins' *mass* capacity `wcapcan = rho_air*can_depth [kg/m2]`;
- **dry-air molar density** `can_dmol = rho_air*(1 − can_shv)/mmdry  [mol_dryair/m3]` (see §1.4);
- **net biotic source (source-positive)** `f_bio = Reco − GPP  [umol/m2/s]` (respiration adds CO2,
  GPP removes it);
- **molar atmosphere↔CAS conductance** `gatm_co2 = can_dmol*ustar  [mol_dryair/m2/s]` — the CO2 analog
  of the mass conductance `gatm = rho_air*ustar`, using the **same** `ustar` so all three twins share
  one turbulent basis.

This is exactly ED2's fast-loop CO2 ODE `d(can_co2)/dt = (nee_tot + cflxac)*ccapcani`
(`rk4_derivs.f90:2163`), with `cflxac = gatm_co2*(co2_atm − can_co2)` the atmosphere→CAS eddy flux and
`nee_tot` the net biological exchange. ED2's molar capacity is `ccapcan = can_dmol*can_depth`
(`canopy_struct_dynamics.f90:2181`), `can_dmol = can_rhos*mmdryi` — dry-air molar density — confirming
the molar (not mass) capacity is the reference structure.

### 1.2 NEE / NEP definitions (sign discipline)

- **NEE** (net ecosystem exchange, atmospheric sign) `= f_bio = Reco − GPP`. **Positive = net source
  to the atmosphere.**
- **NEP** (net ecosystem production, biosphere sign) `= −NEE = GPP − Reco`. **Positive = net ecosystem
  uptake.** (Kept as a derived diagnostic; see the naming note in §2.5.)
- **loss2atm** `= gatm_co2*(can_co2_new − co2_atm)  [umol/m2/s]`. **Positive = the CAS vents CO2
  upward** to the free atmosphere (ED2's `−cflxac`, `rk4_derivs.f90:2236`). At night with respiration
  only and weak `ustar`, `can_co2 > co2_atm` and `loss2atm > 0`.
- **Storage** `= ccapcan*can_co2  [umol/m2]` — the physical CAS CO2 inventory.

Carbon-mass outputs use `umol_2_kgC = 1.20107e-8 kgC/umol` (molar mass of C, 12.0107 g/mol) — applied
**only at the I/O seam**, never inside a kernel.

### 1.3 ED2 reference vs CLM5 / FATES / CliMA — is the prognostic store worth it?

| Model | Canopy CO2 for photosynthesis | Prognostic CAS-CO2 store? | Soil / hetero respiration | NEE |
|---|---|---|---|---|
| **ED2** (MEDS reference) | `ca = can_co2` (the CAS value) | **Yes** — `can_co2 [umol/mol]`, molar `ccapcan`, atm exchange via `ustar`, closed `co2budget` | in-model, fast modifier × slow pools | closed inside the CAS: `nee_tot = Reco − GPP`, `loss2atm = −cflxac` |
| **CLM5** | `co2(p) = forc_pco2(g)` (pinned to forcing) | **No** — leaf CO2 drops from `ca` only via the `1.4/gb + 1.6/gs` diffusion cascade | host soil BGC (CENTURY / MIMICS) | diagnostic sum `Ra + Rh − GPP` |
| **CTSM-FATES** | single well-mixed atmospheric CO2 held constant through the canopy | **No** | handed to the host LSM | assembled at host level |
| **CliMA Land** | prescribed `c_co2` as leaf `ca`; a *separate* soil-CO2 diffusion column vents to the atmosphere, **decoupled** from the canopy | **No** (canopy); soil column is one-way to atm | DAMM in the soil column | no NEE closure through a canopy store |

**Verdict.** A prognostic CAS-CO2 store is an **ED2-lineage feature** the reference land models omit.
It is *optional* for GPP/photosynthesis parity (for CLM/FATES-like behavior, pin leaf `ca` to the
forcing and skip this module). It is *required* for MEDS's stated goals: ED2 fidelity, a closed CO2
budget with a machine-precision residual, nocturnal sub-canopy CO2 accumulation, and the
`ca = can_co2` feedback on leaf `cs`/`ci`. **Build it.** From CliMA, adopt only the soil-CO2 *source*
structure (DAMM/Q10 heterotrophic respiration) — but redirect its top boundary into the CAS store,
not the atmosphere.

### 1.4 Why dry-air molar density (and why it is exact here)

`can_co2` is defined per **mole of dry air**, so its capacity must be the dry-air molar column. Given
the CAS *moist* air density `rho_air [kg/m3]` and specific humidity `can_shv [kg vapor/kg moist air]`,
the dry-air mass fraction is `(1 − can_shv)` **by definition of specific humidity**, so

```
rho_dry  = rho_air*(1 − can_shv)                          [kg_dryair/m3]     (exact)
can_dmol = rho_dry / mmdry = rho_air*(1 − can_shv)/mmdry  [mol_dryair/m3]
```

This is **exact** and needs no CAS pressure field — it reuses `rho_air` (already in
`cas_atm_forcing_t`) and `can_shv` (the second twin). It reproduces ED2's `can_dmol = pdry/(rmol*T)`
without carrying `pdry` separately. Using *moist* `rho_air` (dropping the `(1 − can_shv)` factor) is a
silent ~1–2 % mixing-ratio bias — one of the biggest risks (§9). At MVP with the twins uncoupled,
`can_shv` may be ~0, giving `can_dmol = rho_air/mmdry` (dry = moist), which is safe.

**Reference-level vs CAS-level density (documented divergence).** `rho_air` in `cas_atm_forcing_t` is
the reference-level (near-surface atmospheric) density; ED2 recomputes `can_dmol` from the CAS state
itself (`can_prss`/`can_temp`/`can_shv`), a few percent different. MEDS uses the reference-level
`rho_air` for the CO2 twin **exactly as the energy twin already does** (`meds_column_energy.f90:267`
uses `rho_air` in `gatm`), so all three twins stay consistent. When a genuine CAS air density becomes
available at P3, revisit **all three twins together**, not CO2 alone.

### 1.5 The dropped stability coefficient `c3` (a documented MVP simplification)

ED2's atmosphere↔CAS carbon flux is `cflxac = can_dmol*ustar*cstar` with `cstar = c3*(atm − can)`,
where `c3` is a dimensionless Monin-Obukhov scalar-transfer coefficient (< 1 under stable nocturnal
stratification). MEDS folds `c3 → 1` into `gatm_co2 = can_dmol*ustar`, **matching what the energy twin
already does** (`meds_column_energy.f90:267` has no `c3`). Consequence: under stable nocturnal
conditions the model under-traps CO2 relative to ED2, so **nocturnal sub-canopy accumulation is
qualitatively present but quantitatively damped** until a *shared* `c3`-corrected conductance lands at
P2 for all three twins at once. The executive-summary "nocturnal buildup" claim is scoped
accordingly.

---

## 2. Placement, libraries, files, naming

### 2.1 Decision

**The column CO2 balance MODULE lives entirely in `src/biogeochemistry/`; only the prognostic
`can_co2` field and the `co2_atm` forcing are added to biophysics.** This honors the stated placement,
keeps every kernel stateless / `pure` / GPU-eligible, and — with the charter amendment in §2.4 —
keeps the fast/slow domain boundary clean and honest.

| File | Contents | Cadence | Library |
|---|---|---|---|
| `src/biogeochemistry/meds_biogeochem_types.f90` | `co2_opts_t`, `soil_carbon_t`, `column_co2_budget_t`, `cohort_co2_flux_t`; `HR_*` selector params | types | `meds_biogeochemistry` |
| `src/biogeochemistry/meds_column_co2.f90` | **FAST pure kernels**: `canopy_air_co2_update`, `aggregate_cohort_co2_fluxes`, `heterotrophic_respiration_flux`, `column_co2_step` | fast (every substep) | `meds_biogeochemistry` |
| `src/biogeochemistry/meds_soil_carbon.f90` | **SLOW** `soil_carbon_step` (litter in − decomposition out); CENTURY seam — **P2 stub** | slow (daily) | `meds_biogeochemistry` |
| `src/biophysics/meds_biophysics_types.f90` | **one line each**: `can_co2` on `cas_state_t`; `co2_atm` on `cas_atm_forcing_t` | (edit) | `meds_biophysics` |
| `src/shared/meds_constants.f90` | add `mmdry`, `umol_2_kgC`, `kgC_2_umol`, `kgCday_2_umols`, `umols_per_kgCyr` | (edit) | `meds_shared` |
| `test/test_column_co2.f90` | unit tests (§7) | — | links `meds_biogeochemistry` |

New library (mirror the biophysics block at `CMakeLists.txt:96–100`):

```cmake
#----- Self-contained biogeochemistry kernels (links shared only, like biophysics). ---------#
file(GLOB BIOGEOCHEM_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/biogeochemistry/*.f90)
add_library(meds_biogeochemistry STATIC ${BIOGEOCHEM_SOURCES})
target_link_libraries(meds_biogeochemistry PUBLIC meds_shared)
meds_fortran_flags(meds_biogeochemistry)
```

Test wiring (mirror the energy-test foreach at `CMakeLists.txt:227–230`):

```cmake
add_executable(test_column_co2 test/test_column_co2.f90)
target_link_libraries(test_column_co2 PRIVATE meds_biogeochemistry)
meds_fortran_flags(test_column_co2)
add_test(NAME column_co2 COMMAND test_column_co2)
```

At P3, add `meds_biogeochemistry` (and `meds_biophysics`) to the `meds_aux` link line so the driver
can call the kernels.

### 2.2 Why the code home is decoupled from the state home (no forbidden edge)

- **Library wall holds.** Both `biophysics/` and `biogeochemistry/` link `src/shared` **only**.
  `cas_state_t` / `cas_atm_forcing_t` live in `meds_biophysics_types` (biophysics), so a
  biogeochemistry module cannot `use` them. Resolution: the kernels take **bare inout scalars** — the
  *shipped* `canopy_air_update` idiom (its actual signature is bare `cas_enthalpy, cas_shv, cas_temp`,
  **not** `type(cas_state_t)`; `meds_column_energy.f90:253`). The CO2 kernels need only `meds_kinds`
  and `meds_constants` (both shared). **No forbidden edge.**
- **State rides with its physical siblings.** `can_co2` is added to `cas_state_t` (one line) so that at
  P3 it threads through the per-patch SoA and the fusion/disturbance lockstep **beside**
  `can_enthalpy` / `can_shv` — they share `can_depth`, `ustar`, `rho_air`, and eventually
  `ca = can_co2`. The driver (`meds_aux`, links both libraries) owns the `cas_state_t` and passes
  `cas%can_co2` by reference into the biogeochem kernel.

### 2.3 The rejected alternative (defensible)

**Co-locate the twin kernel in `meds_column_energy.f90`** as the next `pure subroutine` after
`canopy_air_update` (design Angle B). This maximizes consistency with the energy code and keeps the
three twins reviewed as one unit. **Why rejected:** it splits the "column CO2 balance module" across
two libraries (twin in biophysics, soil-C source in biogeochemistry), separates the twin from the
soil-carbon source it is budgeted against, and forces the `column_co2_step` assembler up into the
driver layer (it may not call a biophysics routine from biogeochemistry). Because the kernel is
stateless, we lose *nothing* physical by housing the logic in biogeochemistry — the shared
`can_depth` / `ustar` / `rho_air` arrive as forced scalar inputs exactly as the sibling temperatures
already do for the energy kernels ("the coupled fixed point is deferred to P3"). This is the one
genuine placement fork; it is settled here in favor of the stated placement, with the alternative
recorded for the record.

### 2.4 Charter amendment — biogeochemistry owns fast *and* slow carbon

The current `src/biogeochemistry/README.md` scopes the library to *"slow, stateful pools ... kept
separate from `biophysics/` so the fast/slow timescale split stays clean."* Placing a fast CAS-CO2
kernel here would violate that charter as written. **Resolve it honestly by widening the charter,
not by smuggling a fast kernel into a slow-only library:** carbon has *both* a fast component (CO2
exchange between the surfaces, the CAS, and the atmosphere) and a slow component (litter and soil
organic-matter decomposition). The clean boundary is therefore **by domain, not by timescale**:

- `biophysics/` — the fast **energy / water / momentum** physics of the column (radiation, CAS heat &
  moisture, soil thermal & hydrology).
- `biogeochemistry/` — the **carbon / nutrient** cycle of the column, spanning the fast CAS-CO2
  exchange *and* the slow soil-C / litter pools.

Update `src/biogeochemistry/README.md` and the `CLAUDE.md` source-layout section accordingly (the
`CLAUDE.md` library-DAG line currently omits `biogeochemistry` entirely — add it as a shared-only
sibling of `biophysics`). The stateless / `pure` / GPU-eligible discipline is unchanged; the fast/slow
distinction survives *within* the library as the file split (`meds_column_co2.f90` fast,
`meds_soil_carbon.f90` slow) and the `intent(in)` read-only-pool contract (§9).

### 2.5 Naming

Per the MEDS convention (spell out unconventional ED2 acronyms; keep established domain tokens):
spell out `heterotrophic_respiration`, `gross_primary_prod`, `growth_respiration`,
`storage_respiration`, `leaf_respiration`, `stem_respiration`, `root_respiration`; keep the accepted
tokens `co2`, `nplant`, `pft`, `dbh`, and the diagnostic tokens `nee` / `gpp`. **`nep`** is not on the
sanctioned token list — keep it as a *derived* diagnostic (`nep = −nee`), computed at the budget/I/O
seam and documented as such, rather than introducing a new prognostic identifier. `biogeochem` is the
accepted library/domain token (directory `biogeochemistry`, parallel to `biophysics`).

---

## 3. Prognostic state + derived types

### 3.1 Biophysics edits (two one-line additions)

```fortran
!----- cas_state_t (meds_biophysics_types.f90:223-228) — ADD the third twin beside its siblings. -----!
type :: cas_state_t
   real(wp) :: can_enthalpy = 0.0_wp     !< [J/kg]     specific enthalpy (PROGNOSTIC)
   real(wp) :: can_shv      = 0.0_wp     !< [kg/kg]    specific humidity (PROGNOSTIC twin)
   real(wp) :: can_co2      = 400.0_wp   !< [umol/mol] dry-air CO2 mixing ratio (PROGNOSTIC third twin)
   real(wp) :: can_temp     = 0.0_wp     !< [K]        diagnosed
   real(wp) :: can_depth    = 20.0_wp    !< [m]        CAS depth (from canopy height; forcing)
end type cas_state_t

!----- cas_atm_forcing_t (meds_biophysics_types.f90:259-264) — ADD the free-atmosphere reference. ----!
   real(wp) :: co2_atm = 400.0_wp        !< [umol/mol] reference-level (free-atmosphere) CO2
```

`can_co2` is directly prognostic with **no diagnosed sibling** (unlike `can_temp`, which is inverted
from `can_enthalpy`+`can_shv`); it behaves like `can_shv`. No `can_dmol` field is added — it is
computed inline from `rho_air` + `can_shv` (§1.4), so the forcing type gains only `co2_atm`. Both
fields are **defaulted components**, so existing energy/hydrology tests that build a `cas_state_t` /
`cas_atm_forcing_t` continue to compile and run unchanged.

While here, annotate the `leaf_flux_t` comments (`meds_plant_types.f90:63-69`) so the aggregation's
leaf-area basis is self-documenting: `a_net` / `a_gross` / `rd` are `[umol CO2/m2 leaf/s]` (the word
"leaf" is currently implicit, unlike `stem_resp` / `root_resp` which already say "per plant").

### 3.2 Biogeochemistry types (`meds_biogeochem_types.f90`)

```fortran
module meds_biogeochem_types
   use meds_kinds, only : wp, ik
   implicit none
   private
   public :: HR_Q10, HR_EXP_ED2
   public :: co2_opts_t, soil_carbon_t, column_co2_budget_t, cohort_co2_flux_t

   !----- Heterotrophic-respiration MODEL selector codes (the `hr_model` field; see §4c′). -------!
   integer(ik), parameter :: HR_Q10     = 1_ik   !< Collatz/K13 Q10  q10**((T-T_ref)/10), ref 15 C
   integer(ik), parameter :: HR_EXP_ED2 = 2_ik   !< ED2 scheme-0 capped exponential min(1, exp(a*(T-T_sat)))
   integer(ik), parameter :: HR_DAMM    = 3_ik   !< Davidson 2012 dual-Arrhenius/Michaelis-Menten (§4c′; mechanistic moisture)

   !----- Slow, stateful per-patch soil-carbon pool (written DAILY, read-only in the fast loop). ----!
   !      MVP: one lumped decomposable pool. The type is shaped so the CENTURY expansion is additive.  !
   type :: soil_carbon_t
      real(wp) :: fast_soil_carbon = 0.0_wp   !< [kgC/m2] lumped decomposable pool (MVP)
      ! P2/P3 CENTURY seam: structural_soil_carbon, slow_soil_carbon, (+microbial/passive, lignin) ...
   end type soil_carbon_t

   !----- Pre-summed patch-ground cohort CO2 fluxes (filled by aggregate_cohort_co2_fluxes). -------!
   type :: cohort_co2_flux_t
      real(wp) :: gross_primary_prod  = 0.0_wp   !< [umol/m2/s] SUM a_gross*leaf_area*nplant (uptake, >=0)
      real(wp) :: leaf_respiration    = 0.0_wp   !< [umol/m2/s] SUM rd*leaf_area*nplant       (source, >=0)
      real(wp) :: stem_respiration    = 0.0_wp   !< [umol/m2/s] SUM stem_resp*nplant          (source, >=0)
      real(wp) :: root_respiration    = 0.0_wp   !< [umol/m2/s] SUM root_resp*nplant          (source, >=0)
      real(wp) :: growth_respiration  = 0.0_wp   !< [umol/m2/s] committed daily growth resp, amortized; MVP=0
      real(wp) :: storage_respiration = 0.0_wp   !< [umol/m2/s] committed daily storage resp, amortized; MVP=0
   end type cohort_co2_flux_t

   !----- Column CO2 budget + diagnostics (mirrors energy_flux_t / chydro_flux_t). ----------------!
   type :: column_co2_budget_t
      real(wp) :: nee      = 0.0_wp   !< [umol/m2/s] Reco - GPP  (atmospheric sign: >0 = source to atm)
      real(wp) :: nep      = 0.0_wp   !< [umol/m2/s] GPP - Reco  (= -nee; >0 = ecosystem uptake; DERIVED)
      real(wp) :: loss2atm = 0.0_wp   !< [umol/m2/s] CAS->free-atm venting = gatm_co2*(can_co2 - co2_atm)
      real(wp) :: storage  = 0.0_wp   !< [umol/m2]   ccapcan*can_co2 (physical CAS CO2 inventory)
      real(wp) :: resid    = 0.0_wp   !< [umol/m2]   closed-budget residual (~0 by construction)
      ! P2: real(wp) :: denseffect = 0.0_wp  !< [umol/m2] density-effect correction (variable can_dmol/can_depth)
   end type column_co2_budget_t

   !----- Pre-extracted CO2/decomposition selectors + parameters (NOT the whole config). ---------!
   !      In-type defaults so the standalone kernels/tests compile pre-P3 (like soil_opts_t).       !
   type :: co2_opts_t
      integer(ik) :: hr_model              = HR_Q10        !< {HR_Q10, HR_EXP_ED2, HR_DAMM}; §4c′ adds a `damm` sub-block
      real(wp)    :: rh_k_base             = 0.0_wp        !< [1/day]  effective decomposition rate (x pool)
      real(wp)    :: rh_q10                = 1.5_wp        !< [-]      Collatz/K13 Q10                (HR_Q10)
      real(wp)    :: rh_t_ref              = 288.15_wp     !< [K]      15 C Q10 reference             (HR_Q10)
      real(wp)    :: resp_temp_increase    = 0.0757_wp     !< [1/K]    ED2 scheme-0 slope           (HR_EXP_ED2)
      real(wp)    :: resp_temp_ref         = 318.15_wp     !< [K]      45 C saturation              (HR_EXP_ED2)
      real(wp)    :: resp_opt_water        = 0.8938_wp     !< [-]      moisture optimum (relative)
      real(wp)    :: resp_water_below_opt  = 5.0786_wp     !< [-]      dry-side exponential slope
      real(wp)    :: resp_water_above_opt  = 4.5139_wp     !< [-]      wet-side (anoxia) exponential slope
      real(wp)    :: co2_atm_ref           = 400.0_wp      !< [umol/mol] fixed atm CO2 when not met-forced
      real(wp)    :: rtol = 1.0e-8_wp                      !< [-]       relative closure tolerance
      real(wp)    :: atol = 1.0e-3_wp                      !< [umol/m2] absolute closure floor
      logical     :: debug_error = .false.
   end type co2_opts_t
end module meds_biogeochem_types
```

### 3.3 Constants to add (`meds_constants.f90`; `day_sec` line 18, `yr_sec` line 19 already present)

```fortran
!----- Canopy-air CO2 balance + soil respiration (carbon <-> mole conversions). --------!
real(wp), parameter :: mmdry          = 0.0289655_wp             !< [kg/mol] dry-air molar mass (~28.97 g/mol)
real(wp), parameter :: umol_2_kgC     = 1.20107e-8_wp            !< [kgC/umol] carbon mass per umol CO2 (12.0107 g/mol)
real(wp), parameter :: kgC_2_umol     = 1.0_wp / umol_2_kgC      !< [umol/kgC]  (~8.3259e7)
real(wp), parameter :: kgCday_2_umols = kgC_2_umol / day_sec     !< [ (umol/m2/s) per (kgC/m2/day) ]  (~963.6)
real(wp), parameter :: umols_per_kgCyr = kgC_2_umol / yr_sec     !< [ (umol/m2/s) per (kgC/m2/yr)  ]  (~2.638)
```

`mmco2 ≈ 0.0440095 kg/mol` is *not* required — nothing on the fast path converts µmol↔kg-CO2; add it
only if a mass-of-CO2 output is ever requested. `umols_per_kgCyr` is the missing conversion for the
`gpp_ref [kgC/m2 leaf/yr]` stub (§4b/P1).

---

## 3.5 Multi-layer canopy-air-space forward-design — the CO2 twin pioneers the shared CAS column

*The single well-mixed box of §1.1 / §4a is not a different model from a resolved canopy — it is the
`n_can_layer = 1` degenerate case of a 1-D vertical CAS transport column that mirrors the soil column
**exactly** (`soil_energy_flux` + the bare-array `soil_heat_be_step` over `soil_params_t`, advanced by
`meds_soil_solver%thomas_solve`; `meds_column_energy.f90:44,123`). This subsection forward-designs that
column so `can_co2` becomes layer-capable with **zero change to the shipped scalar behavior**, and so the
canopy-air grid + eddy-diffusivity `K(z)` it introduces is the **shared** transport abstraction the
enthalpy and humidity twins adopt later. CO2 pioneers it because it is the newest twin and has no
diagnosed sibling to invert (§3.1).*

This is a **forward-design (P2 target), not an MVP deliverable.** P0/P1 ship the `n = 1` scalar box of
§4a unchanged; everything here is the reserved multi-layer upgrade slot, held to the same MVP-honesty
discipline (each deferred piece is named, each closes `resid ~ 0`, nothing is claimed working before it is).

### 3.5.1 Governing equation and the closure choice (first-order K-theory MVP)

Every multilayer canopy model (CLM-ml v0/v1, CANOAK/CANVEG) advances a 1-D vertical scalar-conservation
equation for each transported CAS quantity `C(z)` — the same PDE for heat, vapor, and CO2:

```
rho_dry * dC/dt = d/dz( rho_dry * K(z) * dC/dz ) + S(z)          [umol/mol/s for CO2]
```

- `C` = dry-air CO2 mole fraction `[umol/mol]` (the intensive `can_co2`; no conversion),
- `rho_dry` = dry-air **molar** density `can_dmol = rho_air*(1 − can_shv)/mmdry [mol/m3]` (§1.4),
- `K(z)` = turbulent eddy diffusivity `[m2/s]` (§3.5.5),
- `S(z)` = leaf-area-density-weighted net biotic source `[umol/m3/s]` = `a(z)*(R_leaf − A)` (§3.5.9),
  with the **soil CO2 efflux entering as the bottom-face flux BC**, not an interior source.

**Closure = first-order (down-gradient) K-theory** — the flux is `F_c = −rho_dry*K(z)*dC/dz`. This is
the only closure that discretizes to the **fixed-size, implicit backward-Euler, harmonic-face,
tridiagonal (Thomas)** solve MEDS already ships for soil, so it reuses `meds_soil_solver%thomas_solve`
and a canopy grid built exactly like `soil_params_t`. The single well-mixed box is its `n = 1` limit
(§3.5.8, proven bit-identical).

**Documented K-theory limitation + the deferred alternative (MVP honesty).** Down-gradient K-theory is
known to *fail inside canopies*: transport is dominated by canopy-scale coherent eddies larger than the
scale over which `dC/dz` varies, producing **counter-gradient** daytime CO2 transport that a positive
`K(z)` cannot represent. The rigorous fix is Raupach's (1989) **localized-near-field / Lagrangian**
theory (CANOAK/CANVEG), a precomputed dense **dispersion matrix** `C_i − C_ref = Σ_j S_j·dz_j·D_ij` —
higher fidelity but **not tridiagonal** (dense solve, stability-binned lookup tables, poor Thomas/GPU
fit). **Verdict:** K-theory is the tractable first cut; the Lagrangian dispersion-matrix is the reserved
fidelity slot, carrying a standing note that *daytime sub-canopy CO2 profiles are a known K-theory bias*
(the CAS analog of §1.5's `c3` and §4b's committed-respiration honesty caveats).

### 3.5.2 Structural prerequisite — promote the Thomas sweep + CAS grid to `meds_shared`

**The one genuinely non-additive edit of this section (cf. §2.2), resolved honestly, not smuggled.** The
`n = 1` scalar twin of §4a is a bare-scalar kernel needing no solver, so "biogeochemistry links `shared`
only" holds today. The **multi-layer** extension needs `thomas_solve`, which currently lives in
`meds_soil_solver` (**biophysics**) and is dimensioned to `n_soil_layer_max` (from
`meds_biophysics_types`, also biophysics). A biogeochemistry module (`meds_column_co2`) **cannot `use`
either**. Since `thomas_solve`'s only callers are `soil_energy_flux` and `column_hydrology_flux` (both
biophysics), promotion is safe. This also *realizes* the maintainer's mandate that the grid + `K(z)` be
**one shared abstraction for all three CAS twins**. Promote to `meds_shared`:

| Promote to `meds_shared` | New module | Rationale |
|---|---|---|
| `thomas_solve`, rewritten **assumed-size** (`a(*)…`) with fixed locals sized `n_column_layer_max`; zero **only** `x(1:n)` (the current whole-array `x = 0.0_wp` at `meds_soil_solver.f90:27` would overrun a shorter actual) | `meds_column_solver` | any column (soil OR canopy) calls it; device-safe (assumed-size = by-address, fixed-size locals stay OpenMP-target eligible) |
| `n_can_layer_max = 20`, master ceiling `n_column_layer_max = 20` (must satisfy `n_column_layer_max ≥ n_soil_layer_max` and `≥ n_can_layer_max`, asserted at init) | `meds_column_solver` | solver work-array ceiling shared by soil + canopy |
| `canopy_air_grid_t` + `init_canopy_air_grid`, the shared bare-array transport operator `canopy_air_be_step`, and `eddy_diffusivity_profile` | `meds_cas_transport` | the shared grid + `K(z)` + BE-step all three twins consume |

**Soil callers are unaffected** — their `n_soil_layer_max` arrays (≤ `n_column_layer_max`) are valid
actuals for the assumed-size dummies. Biophysics may keep its own `n_soil_layer_max` alias
(`= n_column_layer_max` via `use`, or an independent value `≤` it). This is the one structural edit that
lets a shared-only module run a multi-layer solve at all; the rest of §3.5 is additive. **It is deferred
to P2** — P0/P1 never build the solver-shared path (the scalar box of §4a needs no solver).

### 3.5.3 The canopy-air grid type `canopy_air_grid_t` (analog of `soil_params_t`)

Positive-z upward (canopy top at `k = 1`, ground at `k = n+1`) — the sign mirror of soil's negative-z,
but `dz`/`dz_node` are positive magnitudes in both, so the solver is sign-agnostic. Indexing keeps
`k = 1` atmosphere-adjacent so the **top face is the atmosphere-exchange face** and the **bottom face is
the ground face**, making the `r(1)`/`r(n)` BC placement bit-identical to soil.

```fortran
!----- Fixed-size 1-D canopy-air column geometry + source weights (meds_cas_transport). --------!
!      Analog of soil_params_t (meds_biophysics_types.f90:137); positive-z (k=1 canopy top).      !
integer(ik), parameter :: n_can_layer_max = 20_ik      !< compile-time CAS column ceiling (mirror n_soil_layer_max)

type :: canopy_air_grid_t
   integer(ik) :: n_active = 1_ik                          !< active layers; 1 = well-mixed box (DEFAULT)
   real(wp) :: can_layer_z(n_can_layer_max+1) = 0.0_wp     !< [m] interface heights >= 0 (k=1 top .. k=n+1 ground=0)
   real(wp) :: z_can(n_can_layer_max)   = 0.0_wp           !< [m] node (mid) heights (> 0)
   real(wp) :: dz_can(n_can_layer_max)  = 0.0_wp           !< [m] layer thickness (> 0); sum(dz_can) = can_depth
   real(wp) :: dz_node(n_can_layer_max) = 0.0_wp           !< [m] internode spacing (> 0); dz_node(n)=dz_can(n) filler
   real(wp) :: lad(n_can_layer_max)     = 1.0_wp           !< [-] normalized leaf-area-density source weights, sum = 1
end type canopy_air_grid_t
```

Builder `init_canopy_air_grid` mirrors `init_soil_params` (`meds_soil_parameters.f90:144-166`)
term-for-term with `+can_depth` replacing `−soil_depth` and an `lad` normalization (`Σ lad = 1`)
replacing the texture broadcast. `lad(k)` is the source-distribution analog of soil's
`root_frac(k)` (`meds_biophysics_types.f90:153`); at `n = 1`, `lad(1) = 1`. **Dependency note:** the
per-layer `lad(z)` profile needs a multi-layer crown/leaf-area vertical structure — the *same* one a
multi-layer canopy RT will need — so building `lad` is coupled to that P2 canopy-vertical-grid work
(§3.5.9), not a standalone add.

### 3.5.4 Layer-capable `can_co2` + the twin-symmetry plan

At P0/P1 `can_co2` stays the **scalar** of §3.1 (bit-identical to the shipped box). The multi-layer
migration is a **P2** step: `can_co2` becomes an array with `n_can_layer = 1` the default, so every P0/P1
caller and test remains bit-identical (the array shape is inert until the grid is configured), and —
ideally in lockstep with the enthalpy/humidity twins — the three twins go multi-layer on the same grid:

```fortran
!----- cas_state_t (meds_biophysics_types.f90:223) — the LAYER-CAPABLE third twin (P2 form). ----!
real(wp) :: can_co2(n_can_layer_max) = 400.0_wp   !< [umol/mol] dry-air CO2 mixing ratio (PROGNOSTIC; layer 1 = box)
```

**Twin symmetry (the shared-abstraction payoff):** the *same* `canopy_air_grid_t` and the *same*
`K(z)`/`kcond` face conductances transport all three twins — CLM-ml solves temperature and vapor on one
shared aerodynamic-conductance set, and CO2 is the identical third scalar. When enthalpy/shv go
multi-layer, they become arrays on the **same grid**, call the **same** `canopy_air_be_step` (§3.5.6)
with their own capacity density (`rho_air` for mass twins vs `can_dmol` for the mole twin) and their own
source/`s_atm`, and reduce to today's scalar `canopy_air_update` (`meds_column_energy.f90:253`) at
`n = 1`. This mirrors exactly how soil heat and soil water already share `meds_soil_solver` + the
negative-z geometry — CO2 just gets there first.

### 3.5.5 The eddy-diffusivity profile `K(z)` (shared by all three twins)

A separate shared physics helper `eddy_diffusivity_profile(grid, ustar, cum_lai, k_opts, kdiff)` fills
`kdiff(k) [m2/s]` at nodes — decoupled from the transport operator (which never inspects it) and
**never referenced at `n = 1`** (a well-mixed box has no interior gradient). MVP recommendation:
Harman-Finnigan / CLM-ml **exponential-decay** K-theory, the canonical multilayer closure:

```
K(z)   = (Sc/Pr) * ustar * lm * exp( (z − h)/lambda ),   z <= h       [m2/s]
beta   = ustar / u(h)                     [-]   neutral ~ 0.30 (HF cubic in Lc/L)
Lc     = h / (cd*(LAI + SAI)) = 1/(cd*a)  [m]   canopy drag length
lm     = 2*beta^3 * Lc                    [m]   within-canopy mixing length
lambda = lm/beta = 2*beta^2 * Lc          [m]   K(z) e-folding depth
K(h)  ~ ustar * lm ~ k_von*ustar*(h − d)  [m2/s] canopy-top value (Kaimal & Finnigan 1994)
```

**Ship a fixed neutral `beta` as the day-1 stub**; add the full HF `beta(Lc/L)` stability solver at a
later phase — and that solver is *also* where the deferred Monin-Obukhov `c3` correction (§1.5) naturally
lands for **all three twins at once**, since `beta` carries the stability dependence. All `K`-profile
parameters (`beta`, `Sc/Pr`, `cd`) are TOML-exposed under `[cas_turbulence]` (the no-hard-coded-params
rule, §6).

### 3.5.6 The shared bare-array transport operator `canopy_air_be_step`

Line-for-line the `soil_heat_be_step` idiom (`meds_column_energy.f90:123-145`) with **one structural
change**: the soil top BC is pure Neumann (`r(1) += g_top`); the CAS top BC is an **implicit Robin
atmosphere exchange** (flux into layer 1 `= g_atm*(s_atm − s(1))` depends on the unknown `s(1)`), so the
`−g_atm*s(1)` term moves to the diagonal. Scalar-agnostic (transports any CAS quantity `s`), `pure`,
device-eligible, fixed-size:

```fortran
!---------------------------------------------------------------------------------------!
! Shared 1-D implicit-BE CAS transport (meds_cas_transport). Harmonic-face eddy conductance +  !
! explicit volumetric source q_src, IMPLICIT Robin atmosphere top (g_atm, s_atm), Neumann ground !
! bottom (q_ground). Bit-identical to soil_heat_be_step except the two top-boundary lines. Used   !
! by the CO2 twin now; enthalpy/shv adopt it at their multi-layer phase. Bare arrays -> GPU-eligible.!
!---------------------------------------------------------------------------------------!
pure subroutine canopy_air_be_step(s_n, dz, dz_node, kcond, dens, q_src,                   &
                                   g_atm, s_atm, q_ground, dt, n, s_new)
   use meds_kinds,          only : wp, ik
   use meds_column_solver,  only : thomas_solve, n_column_layer_max
   integer(ik), intent(in)  :: n
   real(wp),    intent(in)  :: s_n(n_column_layer_max), dz(n_column_layer_max)         !< [C], [m]
   real(wp),    intent(in)  :: dz_node(n_column_layer_max), kcond(n_column_layer_max)  !< [m], [mol/(m*s)]
   real(wp),    intent(in)  :: dens(n_column_layer_max), q_src(n_column_layer_max)     !< [mol/m3], [C/m3/s]
   real(wp),    intent(in)  :: g_atm, s_atm, q_ground, dt   !< [mol/m2/s], [C], [umol/m2/s], [s]
   real(wp),    intent(out) :: s_new(n_column_layer_max)
   real(wp)    :: kf(n_column_layer_max), a(n_column_layer_max), b(n_column_layer_max)
   real(wp)    :: c(n_column_layer_max), r(n_column_layer_max)
   integer(ik) :: k
   do k = 1_ik, n - 1_ik
      kf(k) = (dz(k) + dz(k+1)) / (dz(k)/kcond(k) + dz(k+1)/kcond(k+1))   ! series-resistor (harmonic) face
   end do
   do k = 1_ik, n
      if (k >= 2_ik)     then ; a(k) = -kf(k-1) / dz_node(k-1) ; else ; a(k) = 0.0_wp ; end if
      if (k <= n - 1_ik) then ; c(k) = -kf(k)   / dz_node(k)   ; else ; c(k) = 0.0_wp ; end if
      b(k) = dens(k) * dz(k) / dt - a(k) - c(k)
      r(k) = dens(k) * dz(k) / dt * s_n(k) + dz(k) * q_src(k)
   end do
   b(1) = b(1) + g_atm            ! IMPLICIT Robin atm exchange (diagonal) -- THE only change vs soil_heat_be_step
   r(1) = r(1) + g_atm * s_atm    ! atmosphere reservoir source (RHS)
   r(n) = r(n) + q_ground         ! ground CO2 efflux (bottom Neumann) -- the geothermal slot of soil_heat_be_step
   call thomas_solve(a, b, c, r, s_new, n)
end subroutine canopy_air_be_step
```

**Encoding the eddy transport as soil's `kappa`.** The molar eddy flux of a mixing ratio is
`F = −(can_dmol·K)·dC/dz`, so the node conductance that plays soil `kappa`'s role is
`kcond(k) = can_dmol*kdiff(k) [mol/(m·s)]` and the per-layer capacity density is `dens(k) = can_dmol
[mol/m3]` → `ccapcan(k) = can_dmol·dz_can(k) [mol/m2]`. The harmonic face `kf(k)` then gives the **face
molar conductance** `g_face(k) = kf(k)/dz_node(k) [mol/m2/s]`. With a uniform (well-mixed) `can_dmol`
this equals `g_face(k) = can_dmol_face·K_face(k)/dz_node(k)` exactly, `can_dmol` factoring out of the
harmonic mean; the variable-density case keeps it inside `kcond` (the more correct form, and identical to
soil). `a(k) = −g_face(k−1)`, `c(k) = −g_face(k)`.

### 3.5.7 The CO2 wrapper `canopy_air_co2_column_update` + closed budget (analog of `soil_energy_flux`)

The biogeochem wrapper (in `meds_column_co2`) assembles the per-layer `dens`/`kcond`/`q_src`, calls the
shared operator, then closes the mole budget by the same telescoping argument as `soil_energy_flux`
(`meds_column_energy.f90:91-110`) — interior faces cancel in the implicit solve, leaving only the top
Robin flux, the bottom Neumann ground efflux, and the distributed biotic sources:

```fortran
subroutine canopy_air_co2_column_update(can_co2, grid, kdiff, can_shv, gpp_layer,          &
                                        plant_resp_layer, heterotrophic_respiration,        &
                                        ustar, co2_atm, rho_air, dt, opts, budget)
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mmdry
   use meds_cas_transport,    only : canopy_air_grid_t, n_can_layer_max, canopy_air_be_step
   use meds_biogeochem_types, only : co2_opts_t, column_co2_budget_t
   real(wp), intent(inout) :: can_co2(n_can_layer_max)            !< [umol/mol] per-layer prognostic twin
   type(canopy_air_grid_t), intent(in) :: grid
   real(wp), intent(in) :: kdiff(n_can_layer_max)                !< [m2/s] eddy diffusivity at nodes (unused at n=1)
   real(wp), intent(in) :: can_shv, ustar, co2_atm, rho_air, dt
   real(wp), intent(in) :: gpp_layer(n_can_layer_max)            !< [umol/m2/s] per-layer GPP (uptake, >= 0)
   real(wp), intent(in) :: plant_resp_layer(n_can_layer_max)     !< [umol/m2/s] per-layer autotrophic respiration
   real(wp), intent(in) :: heterotrophic_respiration             !< [umol/m2/s] soil Rh (bottom Neumann source)
   type(co2_opts_t),          intent(in)  :: opts
   type(column_co2_budget_t), intent(out) :: budget
   real(wp)    :: dens(n_can_layer_max), kcond(n_can_layer_max), q_src(n_can_layer_max), co2_new(n_can_layer_max)
   real(wp)    :: can_dmol, gatm_co2, s0, s1, ftop, f_bio_tot, scale
   integer(ik) :: n, k
   n        = grid%n_active
   can_dmol = rho_air * (1.0_wp - can_shv) / mmdry               ! [mol/m3] DRY-air molar density (exact; §1.4)
   gatm_co2 = can_dmol * ustar                                   ! [mol/m2/s] atm<->CAS exchange (c3 = 1; §1.5)
   do k = 1_ik, n
      dens(k)  = can_dmol                                        ! [mol/m3]    per-layer molar capacity density
      kcond(k) = can_dmol * kdiff(k)                             ! [mol/(m*s)] molar eddy conductivity (kappa analog)
      q_src(k) = (plant_resp_layer(k) - gpp_layer(k)) / grid%dz_can(k)  ! [umol/m3/s] volumetric biotic source
   end do
   !----- Shared 1-D implicit transport: Robin atm top, ground-Rh bottom; n=1 -> scalar box. -----!
   call canopy_air_be_step(can_co2, grid%dz_can, grid%dz_node, kcond, dens, q_src,          &
                           gatm_co2, co2_atm, heterotrophic_respiration, dt, n, co2_new)
   !----- Conservative molar storage + closed CO2 budget (interior faces telescope). -------------!
   s0 = 0.0_wp ; s1 = 0.0_wp ; f_bio_tot = heterotrophic_respiration
   do k = 1_ik, n
      s0 = s0 + can_dmol * can_co2(k) * grid%dz_can(k)
      s1 = s1 + can_dmol * co2_new(k) * grid%dz_can(k)
      f_bio_tot = f_bio_tot + (plant_resp_layer(k) - gpp_layer(k))       ! = Reco - GPP over the column
   end do
   ftop            = gatm_co2 * (co2_atm - co2_new(1))                    ! atm -> CAS flux at the TOP layer
   budget%nee      = f_bio_tot
   budget%nep      = -f_bio_tot
   budget%loss2atm = -ftop                                               ! = gatm_co2*(co2_new(1) - co2_atm)
   budget%storage  = s1
   budget%resid    = (s1 - s0) - dt * (f_bio_tot + ftop)                 ! = 0 by construction
   can_co2(1:n)    = co2_new(1:n)
   scale = max(abs(budget%storage), 1.0_wp)                              ! mixed rtol/atol form (cf. §4d)
   if (opts%debug_error .and. abs(budget%resid) > opts%rtol*scale + opts%atol) then
      error stop 'canopy_air_co2_column_update: CO2 budget did not close'
   end if
end subroutine canopy_air_co2_column_update
```

**Why the budget closes.** The BE update per layer is discrete conservation
`can_dmol·dz_can(k)·(co2_new(k) − can_co2(k))/dt = [face flux in] − [face flux out] + dz_can(k)·q_src(k)
+ BC`. Summing over `k`, interior faces telescope (out-of-`k` = into-`k+1`), leaving only the top Robin
flux `gatm_co2·(co2_atm − co2_new(1))`, the bottom Neumann `heterotrophic_respiration`, and
`Σ(plant_resp_layer − gpp_layer)`. Hence `(s1 − s0)/dt = ftop + f_bio_tot` ⟹ `resid ≡ 0` — the same
argument that certifies `soil_energy_flux`, now keyed on the **top layer's** mixing ratio (not the box
mean) for `loss2atm`. **The two box-era guards survive:** dry-air `can_dmol` via `(1 − can_shv)` (§1.4),
and the `denseffect` multi-step term (§9, risk 2) when `can_dmol`/`dz_can` vary between substeps.

### 3.5.8 The `n_can_layer = 1` reduction is bit-identical to §4a (proof)

With `n = 1` there are **no interior faces**, so `a(1) = c(1) = 0` and `kcond`/`kdiff` are never
referenced. `dz_can(1) = can_depth`, `lad(1) = 1`. The tridiagonal collapses to the single row
`b(1)·x(1) = r(1)`:

```
b(1) = can_dmol*dz_can(1)/dt + gatm_co2  = ccapcan/dt + gatm_co2          (Robin top added to diagonal)
r(1) = ccapcan/dt*can_co2 + (plant_resp - gpp) + heterotrophic_respiration + gatm_co2*co2_atm
     = ccapcan/dt*can_co2 + f_bio + gatm_co2*co2_atm                       (f_bio = Reco - GPP; Rh = bottom Neumann)
```

so, dividing numerator and denominator by `ccapcan/dt` (with `cci = 1/ccapcan`):

```
co2_new = x(1) = r(1)/b(1)
        = ( can_co2 + dt*cci*(f_bio + gatm_co2*co2_atm) ) / ( 1 + dt*cci*gatm_co2 )
```

which is **line-for-line** the §4a `canopy_air_co2_update` update (and the `enth_new` branch of
`canopy_air_update`, `meds_column_energy.f90:269`, under the `÷mmdry` dry-air substitution).
`budget%resid`, `loss2atm`, `storage`, `nee`/`nep` all match term-for-term. **So the shipped scalar box
is the `n = 1` special case with no separate code path** — the multi-layer wrapper at `grid%n_active = 1`
reproduces §4a to the last bit. Test `co2_column_reduces_to_box` (§7, new #12) asserts this to machine
precision — the regression guard that keeps the box the exercised default. *(Verified by hand during
integration: the `b(1)`/`r(1)` collapse yields the §4a formula exactly.)*

### 3.5.9 Source distribution `S(z)` — generalizing `aggregate_cohort_co2_fluxes` to a per-layer bin

Today `aggregate_cohort_co2_fluxes` (§4b) reduces all cohort fluxes to **one** patch total. For `N`
layers, the *same* per-cohort fluxes are binned into canopy-air layer `k` by the fraction of each
cohort's crown (its leaf-area-density profile) in layer `k`'s height band:

```
q_src(k)*dz_can(k) = plant_resp_layer(k) − gpp_layer(k)
                   = Σ_cohorts [ (rd − a_gross)*leaf_area*nplant*frac_LAI_in_k ]   [umol/m2/s]
                   + stem/root respiration placed at the cohort's stem height
```

GPP/`Rd` stay leaf-area-weighted (per `dLAI(k)`); stem/root respiration are `nplant`-weighted at height.
The **soil efflux** (`heterotrophic_respiration` + surface root respiration) is **not** an interior
`q_src` — it is the bottom Neumann `q_ground` into layer `n`. The reduction stays an **exact linear map**
(no residual of its own, §4b budget note), and at `n = 1` all `lad(k)` collapse to `lad(1) = 1`, so the
per-layer sum equals the single `f_bio` — preserving the §3.5.8 identity. The `lad(z)` profile is the
*same* leaf-area-density MEDS already needs for multilayer canopy RT (`meds_canopy_radiation`), so the
crown-binning geometry is shared with the two-stream grid — not new machinery (but it does depend on that
multi-layer canopy vertical grid landing, §3.5.3).

### 3.5.10 Phasing (slots into §8)

- **P0/P1 — the `n = 1` scalar box only** (as §8 already states). `can_co2` stays the **scalar** of §3.1;
  `canopy_air_grid_t`/`K(z)`/`canopy_air_be_step` are **not built**. Zero multi-layer code ships; §4a is
  the whole CAS-CO2 story.
- **P2 — the shared CAS transport column (this subsection), CO2 first.** Promote `thomas_solve` +
  ceilings to `meds_column_solver`; add `meds_cas_transport` (`canopy_air_grid_t`, `init_canopy_air_grid`,
  `eddy_diffusivity_profile`, `canopy_air_be_step`); **migrate `can_co2` to the array form** (§3.5.4); add
  `canopy_air_co2_column_update` to `meds_column_co2`; generalize `aggregate_cohort_co2_fluxes` to
  `lad`-binned per-layer sources; add `[cas_turbulence]` TOML. Ship a **fixed neutral `beta`** `K(z)`;
  `n_can_layer_max = 20` but `n_active = 1` remains the exercised default (test #12). **This is a
  shared-with-energy/water upgrade** — the grid + `K(z)` + operator are `meds_shared`, so the
  enthalpy/humidity twins adopt them next with no new transport code.
- **P3 — state lockstep + coupling.** `can_co2(:)` (and the grid's `n_active`) thread the per-patch CAS
  SoA through `sort_patches`/`patch_compact` beside `can_enthalpy`/`can_shv` (§8-P3); the enthalpy/vapor
  twins go multi-layer on the *same* `canopy_air_grid_t`; the HF `beta(Lc/L)` stability solver (with the
  shared `c3`, §1.5) lands for all three at once. The Lagrangian dispersion-matrix closure (§3.5.1)
  remains the reserved fidelity slot.

**Test additions (§7):** #12 `co2_column_reduces_to_box` (`n_active = 1` reproduces §4a to machine
precision); #13 `co2_column_resid_zero` (randomized `n > 1`, `|resid| < 1e-10·max(storage,1)`); #14
`co2_two_layer_gradient` (Rh at the ground + a GPP sink aloft drives a steady `co2(bottom) > co2(top)`
gradient — the K-theory transport sanity check); all built under **both ifx and nvfortran multicore**
(the bare-array BE-step is exactly the array-temporary miscompile class of CLAUDE.md issue #7).

---

## 4. The kernels

All kernels are `pure`, take bare scalars/arrays, and link `src/shared` only — GPU-eligible, matching
the shipped `canopy_air_update`. Each budget-bearing kernel returns a `resid` that is **zero by
algebraic construction**; the non-`pure` assembler `error stop`s in Debug when `|resid|` exceeds
tolerance (the uniform MEDS discipline — cf. `soil_energy_flux`, `meds_column_energy.f90:113`).

### 4a. `canopy_air_co2_update` — CAS storage + atmospheric flux (implicit-atm, closed residual)

This is the twin update **and** the budget closure in one routine, exactly as `canopy_air_update` both
sums `f_sens` and closes its residual.

```fortran
!---------------------------------------------------------------------------------------!
! Canopy-air-space CO2 update -- the THIRD prognostic twin. Advances the dry-air CO2      !
! mixing ratio one step from the net biotic source (respiration - GPP) and the atmospheric !
! exchange (IMPLICIT in the atm term, like canopy_air_update's enthalpy branch). Molar      !
! capacity because can_co2 is a mole fraction. Fills the closed CO2 budget (resid ~ 0).      !
!---------------------------------------------------------------------------------------!
pure subroutine canopy_air_co2_update(can_co2, can_depth, can_shv,                        &
                                      gross_primary_prod, plant_respiration,              &
                                      heterotrophic_respiration,                          &
                                      ustar, co2_atm, rho_air, dt, budget)
   use meds_kinds,            only : wp
   use meds_constants,        only : mmdry, tiny_num
   use meds_biogeochem_types, only : column_co2_budget_t
   real(wp), intent(inout) :: can_co2                    !< [umol/mol]  prognostic third twin
   real(wp), intent(in)    :: can_depth                  !< [m]         CAS depth (fixed within a step at MVP)
   real(wp), intent(in)    :: can_shv                    !< [kg/kg]     CAS humidity (-> dry-air molar density)
   real(wp), intent(in)    :: gross_primary_prod         !< [umol/m2/s] patch GPP (uptake, >=0)
   real(wp), intent(in)    :: plant_respiration          !< [umol/m2/s] leaf+stem+root+growth+storage (source, >=0)
   real(wp), intent(in)    :: heterotrophic_respiration  !< [umol/m2/s] soil Rh (source, >=0)
   real(wp), intent(in)    :: ustar                      !< [m/s]       friction velocity (shared twin conductance)
   real(wp), intent(in)    :: co2_atm                    !< [umol/mol]  free-atmosphere reference
   real(wp), intent(in)    :: rho_air                    !< [kg/m3]     moist CAS air density
   real(wp), intent(in)    :: dt                         !< [s]
   type(column_co2_budget_t), intent(out) :: budget      !< nee/nep/loss2atm/storage/resid
   real(wp) :: can_dmol, ccapcan, cci, f_bio, gatm_co2, co2_new
   can_dmol = rho_air * (1.0_wp - can_shv) / mmdry        ! [mol_dryair/m3] DRY-air molar density (exact; §1.4)
   ccapcan  = can_dmol * can_depth                        ! [mol_air/m2]    MOLAR CAS capacity (cf wcapcan)
   cci      = 1.0_wp / max(ccapcan, tiny_num)             ! [m2/mol]
   f_bio    = heterotrophic_respiration + plant_respiration - gross_primary_prod  ! [umol/m2/s] Reco - GPP
   gatm_co2 = can_dmol * ustar                            ! [mol_air/m2/s]  atm<->CAS molar exchange (c3 = 1; §1.5)
   !----- Implicit in the atmospheric-exchange term (L-stable); biotic source explicit. --------!
   co2_new  = (can_co2 + dt*cci*(f_bio + gatm_co2*co2_atm)) / (1.0_wp + dt*cci*gatm_co2)
   !----- Diagnostics + closed budget (resid = 0 by substitution). -----------------------------!
   budget%resid    = ccapcan*(co2_new - can_co2) - dt*(f_bio + gatm_co2*(co2_atm - co2_new))  ! = 0
   budget%nee      = f_bio
   budget%nep      = -f_bio
   budget%loss2atm = gatm_co2*(co2_new - co2_atm)
   budget%storage  = ccapcan*co2_new
   can_co2  = co2_new
end subroutine canopy_air_co2_update
```

**Line-by-line correspondence to `canopy_air_update`:** `ccapcan ↔ wcapcan` (×`1/mmdry` and
dry-air), `cci ↔ wci`, `f_bio ↔ f_sens`, `gatm_co2 ↔ gatm` (÷`mmdry`, dry-air), `co2_new ↔ enth_new`
(identical algebra), `budget%resid ↔ resid` (identical algebra, `meds_column_energy.f90:272-273`).

**How the budget closes.** Substituting `co2_new` from the update into `budget%resid`:
```
(1 + dt*cci*gatm_co2)*co2_new = can_co2 + dt*cci*(f_bio + gatm_co2*co2_atm)
  ⟹ ccapcan*(co2_new − can_co2) = dt*(f_bio + gatm_co2*(co2_atm − co2_new))   ⟹ resid ≡ 0.
```
Equivalently `Δstorage = dt*(NEE − loss2atm)` — ED2's `co2budget` closure in one line, with the
`denseffect` / `zcaneffect` terms identically zero at MVP because `can_dmol` and `can_depth` are held
fixed within a step (§9, risk 2).

### 4b. `aggregate_cohort_co2_fluxes` — per-cohort → µmol CO2 m⁻² ground s⁻¹

The per-cohort reduction runs **host-side over the patch's CSR slice** (`cohort_offset` /
`cohort_count`), so the twin kernel stays scalar / `pure`. The four fast fluxes are **natively µmol**
(verified against `meds_plant_types`): leaf photosynthesis/respiration are per leaf area, stem/root
respiration are per plant. **Growth and storage respiration are the only `kgC` terms** — they are
ED2's *committed* daily respiration (`commit_growth_resp` / `commit_storage_resp`), amortized to a
constant `[umol/m2/s]` across the day; both are **0 at MVP** and turn on at P3 when carbon allocation
couples to the fast loop.

| Flux | Source field | Native unit | ×factor → `[umol/m2 ground/s]` |
|---|---|---|---|
| GPP | `leaf_flux_t%a_gross` | `[umol CO2/m2 leaf/s]` | `× leaf_area × nplant` (= LAI weighting) |
| Leaf resp `Rd` | `leaf_flux_t%rd` | `[umol CO2/m2 leaf/s]` | `× leaf_area × nplant` |
| Stem resp | `wood_flux_t%stem_resp` | `[umol CO2/plant/s]` | `× nplant` |
| Root resp | `root_flux_t%root_resp` | `[umol CO2/plant/s]` | `× nplant` |
| Growth resp | allocation `rg` (committed) | `[kgC/plant/day]` | `× nplant × kgC_2_umol / day_sec` (MVP = 0) |
| Storage resp | committed storage resp | `[kgC/plant/day]` | `× nplant × kgC_2_umol / day_sec` (MVP = 0) |

```fortran
pure subroutine aggregate_cohort_co2_fluxes(n, nplant, leaf_area, a_gross, rd,          &
                                            stem_resp, root_resp, out)
   use meds_kinds,            only : wp, ik
   use meds_biogeochem_types, only : cohort_co2_flux_t
   integer(ik), intent(in) :: n
   real(wp), intent(in)  :: nplant(n), leaf_area(n)        !< [plant/m2], [m2 leaf/plant]
   real(wp), intent(in)  :: a_gross(n), rd(n)              !< [umol CO2/m2 leaf/s]
   real(wp), intent(in)  :: stem_resp(n), root_resp(n)     !< [umol CO2/plant/s]
   type(cohort_co2_flux_t), intent(out) :: out
   integer(ik) :: i
   out%gross_primary_prod = 0.0_wp ; out%leaf_respiration = 0.0_wp
   out%stem_respiration   = 0.0_wp ; out%root_respiration = 0.0_wp
   out%growth_respiration = 0.0_wp ; out%storage_respiration = 0.0_wp
   do i = 1_ik, n
      out%gross_primary_prod = out%gross_primary_prod + a_gross(i)  *leaf_area(i)*nplant(i)
      out%leaf_respiration   = out%leaf_respiration   + rd(i)       *leaf_area(i)*nplant(i)
      out%stem_respiration   = out%stem_respiration   + stem_resp(i)             *nplant(i)
      out%root_respiration   = out%root_respiration   + root_resp(i)             *nplant(i)
   end do
end subroutine aggregate_cohort_co2_fluxes
```

**Budget note.** This kernel has no residual of its own — it is an exact linear reduction; its
correctness is guarded by test 7 (a hand-checked multi-cohort patch). Explicit `do`-loop (not
`sum(a_gross*leaf_area*nplant)`) for clarity and to sidestep the nvfortran array-temporary trap
(CLAUDE.md issue #7). `growth_respiration` / `storage_respiration` are set by the caller (§4d), not in
this kernel — they arrive already reduced to `[umol/m2/s]`.

### 4c. `heterotrophic_respiration_flux` — MVP Q10 + moisture modifier (the `HR_Q10`/`HR_EXP_ED2` branch)

> **Superseded framing — see §4c′.** This subsection describes the two *empirical* branches
> (`HR_Q10`, `HR_EXP_ED2`). §4c′ reframes `heterotrophic_respiration_flux` as a pluggable dispatcher over
> `hr_model ∈ {HR_Q10, HR_EXP_ED2, HR_DAMM}` and adds the mechanistic **DAMM** kernel (Davidson et al.
> 2012). The selector field named `hr_temp_response` in §3.2 is renamed `hr_model` there; the kernel body
> below is the empirical tail, unchanged in meaning.

Exploits ED2's fast/slow seam: the soil-carbon pool is **frozen within the fast loop** (`intent(in)`),
so sub-daily variability is *purely* the temperature/moisture modifier. Reads soil state collapsed to
one representative (root-weighted) temperature/moisture at MVP — the per-layer A/B above/below-ground
split defers to P2.

```fortran
!---------------------------------------------------------------------------------------!
! MVP heterotrophic respiration [umol/m2/s]: a frozen decomposable pool x an intrinsic     !
! decay rate x a dimensionless temperature modifier x a dimensionless moisture modifier.    !
! ED2 default scheme-0 (capped exponential) OR Collatz Q10, selectable. Pool state is SLOW   !
! (read-only here); this is the fast coupling flux from the slow pool to the CAS twin.        !
!---------------------------------------------------------------------------------------!
pure function heterotrophic_respiration_flux(fast_soil_carbon, soil_temp, theta,        &
                                             theta_dry, theta_sat, opts) result(rh)
   use meds_kinds,            only : wp
   use meds_constants,        only : tiny_num, kgCday_2_umols
   use meds_biogeochem_types, only : co2_opts_t, HR_EXP_ED2
   real(wp), intent(in) :: fast_soil_carbon      !< [kgC/m2]  decomposable pool (frozen in the fast loop)
   real(wp), intent(in) :: soil_temp             !< [K]       representative (root-weighted) soil temperature
   real(wp), intent(in) :: theta, theta_dry, theta_sat  !< [m3/m3] moisture, air-dry floor, porosity
   type(co2_opts_t), intent(in) :: opts
   real(wp) :: rh                                !< [umol/m2/s]
   real(wp) :: f_temp, f_water, rel
   select case (opts%hr_model)                   ! (renamed from hr_temp_response; §4c′ adds HR_DAMM)
   case (HR_EXP_ED2)                             ! ED2 scheme-0: min(1, exp(a*(T - T_sat)))
      f_temp = min(1.0_wp, exp(opts%resp_temp_increase*(soil_temp - opts%resp_temp_ref)))
   case default                                  ! HR_Q10 (Collatz/K13): q10**((T - T_ref)/10)
      f_temp = opts%rh_q10 ** ((soil_temp - opts%rh_t_ref)*0.1_wp)
   end select
   rel = min(1.0_wp, max(0.0_wp, (theta - theta_dry)/max(theta_sat - theta_dry, tiny_num)))
   if (rel <= opts%resp_opt_water) then          ! one-sided exponential about the optimum
      f_water = exp((rel - opts%resp_opt_water)*opts%resp_water_below_opt)
   else
      f_water = exp((opts%resp_opt_water - rel)*opts%resp_water_above_opt)
   end if
   rh = fast_soil_carbon * opts%rh_k_base * f_temp * f_water * kgCday_2_umols   ! kgC/m2/day -> umol/m2/s
end function heterotrophic_respiration_flux
```

**Budget note.** Rh is a *source term*, not a stored quantity — it has no residual of its own; its
magnitude enters `f_bio` in 4a, whose residual closes. **The `kgCday_2_umols ≈ 963.6` factor is
mandatory** — omitting it is a silent ~1000× error, guarded by test 8.

**The moisture floor — a corrected, documented choice (Lens 1 + Lens 2).** The relative moisture
normalizes against an **air-dry floor `theta_dry`**, *not* the van Genuchten residual `theta_res`.
ED2's `het_resp_weight` uses `soilcp` (air-dry water content, ~ ψ ≈ −3.1 MPa), which for most textures
is a **different (usually larger) floor** than the vG retention residual `theta_res`
(`meds_soil_parameters.f90:56`, `meds_biophysics_types.f90:145`). Passing `theta_res` where an air-dry
floor is meant rescales `rel` and silently mis-tunes the whole moisture response. **MVP policy:** the
caller supplies `theta_dry` as an explicit argument. Two acceptable sources, decided at P1 wiring
time: (a) if `soil_params_t` exposes an air-dry / `soilcp` equivalent, pass it; (b) otherwise derive
`theta_dry = theta_of_psi(ψ ≈ −3.1 MPa)` from the retention curve in `meds_soil_parameters`. Do **not**
default it to `theta_res` and call it `soilcp`. This is tracked as a **correctness dependency** (§9),
verified by test 8's moisture branch.

**The poolless stub** (the leanest day-1 target) is a special case: set `fast_soil_carbon = 1` and
`rh_k_base = R_ref [kgC/m2/day]` — no soil-C state needed. The temperature form is inlined (not routed
through `meds_temp_response`) because ED2's capped exponential is not one of that module's 25 °C-
referenced photosynthetic forms; if a second consumer needs a Q10, promote a
`q10_scale(q10, t, t_ref)` to `meds_temp_response` (shared) at P2. **This kernel uses only
`meds_kinds`, `meds_constants`, and `meds_biogeochem_types` — no `meds_thermo` dependency** (correcting
an over-stated dependency from the draft).

### 4c′. DAMM option + the pluggable soil-gas-flux family (`hr_model`)

*Reframes `heterotrophic_respiration_flux` as one instance of a pluggable soil-gas-flux family (selector
`hr_model ∈ {HR_Q10, HR_EXP_ED2, HR_DAMM}`), adds the mechanistic **DAMM** kernel (Davidson, Samanta,
Caramori & Savage 2012, GCB 18:371-384), and shows the seam stays additive for future CH4/N2O. Extends
§3.2 (co2_opts_t), §6 (config), §7 (tests), §8 (phasing).*

`heterotrophic_respiration_flux` becomes a **thin dispatcher** over `hr_model`. The two existing branches
are **unchanged in meaning** — the §3.2 `hr_temp_response` field is renamed `hr_model` and its codes
`HR_Q10 = 1` / `HR_EXP_ED2 = 2` keep their exact semantics; a third code `HR_DAMM = 3` is added.
`HR_Q10`/`HR_EXP_ED2` share the empirical `pool × rate_base × f_temp(T) × f_water(θ)` form; `HR_DAMM`
computes the **absolute** rate from diffusion physics and takes a *different* argument set, so it branches
cleanly (it must **not** double-apply an empirical `f_water` — DAMM supplies the moisture response
mechanistically).

#### 4c′.1 Selector, extended `co2_opts_t`, and the dispatcher

```fortran
!----- Heterotrophic-respiration MODEL selector (meds_biogeochem_types; renames hr_temp_response). --!
integer(ik), parameter :: HR_Q10     = 1_ik   !< Collatz/K13 Q10  q10**((T-T_ref)/10) x f_water x pool  (unchanged)
integer(ik), parameter :: HR_EXP_ED2 = 2_ik   !< ED2 scheme-0 capped exp min(1,exp(a*(T-T_sat))) x f_water x pool (unchanged)
integer(ik), parameter :: HR_DAMM    = 3_ik   !< Davidson 2012 dual-Arrhenius/Michaelis-Menten (mechanistic moisture)
```

`co2_opts_t` (§3.2) gains the selector rename (`hr_temp_response` → `hr_model`) plus a nested DAMM
parameter block (defaults §4c′.3); the existing empirical fields are untouched:

```fortran
type :: co2_opts_t
   integer(ik) :: hr_model = HR_Q10             !< {HR_Q10, HR_EXP_ED2, HR_DAMM}   (was hr_temp_response)
   real(wp)    :: rh_k_base            = 0.0_wp     !< [1/day] decomposition rate  (HR_Q10 / HR_EXP_ED2 only)
   real(wp)    :: rh_q10               = 1.5_wp     !< [-]     Collatz/K13 Q10      (HR_Q10)
   real(wp)    :: rh_t_ref             = 288.15_wp  !< [K]     15 C Q10 reference   (HR_Q10)
   real(wp)    :: resp_temp_increase   = 0.0757_wp  !< [1/K]   ED2 scheme-0 slope   (HR_EXP_ED2)
   real(wp)    :: resp_temp_ref        = 318.15_wp  !< [K]     45 C saturation      (HR_EXP_ED2)
   real(wp)    :: resp_opt_water       = 0.8938_wp  !< [-]     moisture optimum     (HR_Q10 / HR_EXP_ED2)
   real(wp)    :: resp_water_below_opt = 5.0786_wp  !< [-]     dry-side slope        (HR_Q10 / HR_EXP_ED2)
   real(wp)    :: resp_water_above_opt = 4.5139_wp  !< [-]     wet-side slope        (HR_Q10 / HR_EXP_ED2)
   real(wp)    :: co2_atm_ref          = 400.0_wp   !< [umol/mol]
   real(wp)    :: rtol = 1.0e-8_wp
   real(wp)    :: atol = 1.0e-3_wp
   logical     :: debug_error = .false.
   type(damm_params_t) :: damm                  !< DAMM diffusion parameters (HR_DAMM only; §4c′.3)
end type co2_opts_t
```

The dispatcher keeps the `HR_Q10`/`HR_EXP_ED2` tail **byte-for-byte** as the §4c code; only the `HR_DAMM`
case is new. DAMM ignores `rh_k_base`, `resp_opt_water`, `resp_water_below_opt`, `resp_water_above_opt`,
and `theta_dry` (its diffusion physics supplies the moisture response):

```fortran
pure function heterotrophic_respiration_flux(fast_soil_carbon, soil_temp, theta,           &
                                             theta_dry, theta_sat, opts) result(rh)
   use meds_kinds,            only : wp
   use meds_constants,        only : tiny_num, kgCday_2_umols
   use meds_biogeochem_types, only : co2_opts_t, HR_EXP_ED2, HR_DAMM
   real(wp), intent(in) :: fast_soil_carbon      !< [kgC/m2] decomposable pool (frozen in the fast loop)
   real(wp), intent(in) :: soil_temp             !< [K]      representative (root-weighted) soil temperature
   real(wp), intent(in) :: theta, theta_dry, theta_sat  !< [m3/m3] moisture, air-dry floor, porosity
   type(co2_opts_t), intent(in) :: opts
   real(wp) :: rh                                !< [umol/m2/s]
   real(wp) :: f_temp, f_water, rel
   select case (opts%hr_model)
   case (HR_DAMM)                                ! ---- mechanistic (moisture emergent; theta_dry UNUSED) ----
      rh = heterotrophic_respiration_damm(fast_soil_carbon, soil_temp, theta, theta_sat, opts%damm)
   case (HR_EXP_ED2)                             ! ---- empirical, ED2 scheme-0 capped exp (UNCHANGED) ----
      f_temp  = min(1.0_wp, exp(opts%resp_temp_increase*(soil_temp - opts%resp_temp_ref)))
      rel     = min(1.0_wp, max(0.0_wp, (theta - theta_dry)/max(theta_sat - theta_dry, tiny_num)))
      f_water = water_modifier(rel, opts)        ! one-sided exponential about resp_opt_water (§4c code)
      rh      = fast_soil_carbon * opts%rh_k_base * f_temp * f_water * kgCday_2_umols
   case default                                  ! ---- HR_Q10 (Collatz/K13), UNCHANGED ----
      f_temp  = opts%rh_q10 ** ((soil_temp - opts%rh_t_ref)*0.1_wp)
      rel     = min(1.0_wp, max(0.0_wp, (theta - theta_dry)/max(theta_sat - theta_dry, tiny_num)))
      f_water = water_modifier(rel, opts)
      rh      = fast_soil_carbon * opts%rh_k_base * f_temp * f_water * kgCday_2_umols
   end select
end function heterotrophic_respiration_flux
```

(`water_modifier` is the two-branch `exp` about `resp_opt_water` factored from the §4c body — a cosmetic
extraction so both empirical cases share it verbatim; no behavior change.)

#### 4c′.2 The DAMM kernel `heterotrophic_respiration_damm`

DAMM computes an **aerobic** decomposition rate as an Arrhenius maximum velocity times **two**
Michaelis-Menten saturation factors — one for soluble-C substrate (rises with moisture via a liquid-film
`θ³` diffusion term), one for O2 (falls with moisture via an air-filled-porosity `(θ_sat − θ)^{4/3}` gas
term). The competing moisture dependencies make the response **unimodal with no empirical `f_water`** (the
paper's headline). Seven equations, extracted from the primary PDF and cross-checked against the
`colinaverill/DAMM-model` R reference:

```
Eqn 1  R_Sx  = Vmax_Sx * Sx/(kM_Sx + Sx) * O2/(kM_O2 + O2)         [mgC cm-3 soil h-1]
Eqn 2  Vmax  = alpha_Sx * exp( -Ea_Sx / (R * T) )                  [mgC cm-3 h-1]   (T in K; Ea,R in kJ/mol)
Eqn 4  [Sx_soluble] = p * [Sx_total]
Eqn 5  Sx    = p * [Sx_total] * D_liq * theta^3                    [gC cm-3]        (theta CUBED)
Eqn 6  O2    = D_gas * 0.209 * a^(4/3)                             [cm3 O2 cm-3 air]
Eqn 7  a     = (1 - BD/PD) - theta = theta_sat - theta             [m3/m3]  (air-filled porosity; CLAMP >= 0)
MEDS   Rh    = R_Sx * depth_cm * 231.269                           [umol CO2 m-2 s-1]
```

The MEDS unit chain (R1-research, **independently re-verified during integration**):
`231.269 = (1e4 cm²/m²)/(3600 s/h)·(1000/12.011 µmol/mgC)`. The **respiring depth appears twice** and does
**not** cancel: once mapping the column pool `fast_soil_carbon [kgC/m2]` to a volumetric concentration
`[Sx_total] = fast_soil_carbon·0.1/depth_cm [gC cm⁻³]`, once integrating the volumetric rate over the
depth — physically correct (SOC spread over the depth, flux integrated over it) and non-cancelling because
`Sx` passes through the nonlinear Michaelis term.

```fortran
!---------------------------------------------------------------------------------------!
! DAMM heterotrophic respiration [umol CO2 m-2 s-1] (Davidson et al. 2012 GCB 18:371-384).   !
! Arrhenius Vmax x MM(soluble-C, liquid-diffusion theta^3) x MM(O2, gas-diffusion (theta_sat-  !
! theta)^{4/3}). The unimodal moisture response is EMERGENT -- no empirical f_water. Needs only  !
! soil_temp[K], theta, theta_sat, and the frozen soil-C pool: all MEDS already has. Pure/GPU-ok.  !
!---------------------------------------------------------------------------------------!
pure function heterotrophic_respiration_damm(fast_soil_carbon, soil_temp, theta,           &
                                             theta_sat, damm) result(rh)
   use meds_kinds,            only : wp
   use meds_constants,        only : r_gas_kj, o2_air_frac, damm_flux_factor
   use meds_biogeochem_types, only : damm_params_t
   real(wp), intent(in) :: fast_soil_carbon      !< [kgC/m2] frozen decomposable pool (= DAMM total soil C)
   real(wp), intent(in) :: soil_temp             !< [K]      soil temperature (NATIVE Kelvin -> Arrhenius directly)
   real(wp), intent(in) :: theta                 !< [m3/m3]  volumetric soil moisture (native; no conversion)
   real(wp), intent(in) :: theta_sat             !< [m3/m3]  porosity = operative air-filled-porosity ceiling
   type(damm_params_t), intent(in) :: damm
   real(wp) :: rh                                !< [umol/m2/s]
   real(wp) :: vmax, sx_total, sx, a_air, o2, mm_sx, mm_o2, r_sx
   !----- Arrhenius max velocity (Ea & R BOTH in kJ/mol; r_gas_kj = r_gas*1e-3). ---------------!
   vmax     = damm%alpha_sx * exp( -damm%ea_sx / (r_gas_kj * soil_temp) )       ! [mgC cm-3 h-1]
   !----- Soluble-C substrate: SOC -> volumetric conc, soluble fraction, liquid diffusion (^3). !
   sx_total = fast_soil_carbon * 0.1_wp / damm%depth_cm                         ! [gC cm-3] (kgC/m2 over depth)
   sx       = damm%p_soluble * sx_total * damm%d_liq * theta**3                 ! [gC cm-3] Eqn 4-5
   !----- Oxygen: air-filled porosity CLAMPED >= 0 before the 4/3 power (else NaN when theta>sat).!
   a_air    = max(theta_sat - theta, 0.0_wp)                                    ! [m3/m3] Eqn 7
   o2       = damm%d_gas * o2_air_frac * a_air**(4.0_wp/3.0_wp)                 ! [cm3 O2 cm-3 air] Eqn 6
   !----- Dual Michaelis-Menten (both self-bounded in [0,1); kM>0 => safe at conc=0). ----------!
   mm_sx    = sx / (damm%km_sx + sx)
   mm_o2    = o2 / (damm%km_o2 + o2)
   r_sx     = vmax * mm_sx * mm_o2                                              ! [mgC cm-3 h-1] Eqn 1
   !----- Depth-integrate + convert to MEDS units (231.269 folds cm,h,mgC->umol). --------------!
   rh       = r_sx * damm%depth_cm * damm_flux_factor                          ! [umol CO2 m-2 s-1]
end function heterotrophic_respiration_damm
```

**New constants (`meds_constants.f90`, beside the §3.3 block):**

```fortran
real(wp), parameter :: r_gas_kj        = r_gas * 1.0e-3_wp   !< [kJ/mol/K] gas const in kJ (DAMM Arrhenius; = 8.3145e-3)
real(wp), parameter :: o2_air_frac     = 0.209_wp            !< [-] volume fraction of O2 in air (DAMM Eqn 6)
real(wp), parameter :: damm_flux_factor = 1.0e4_wp/3600.0_wp * 1000.0_wp/12.011_wp
                                                             !< [(umol/m2/s) per (mgC cm-3 h-1 . cm)] ~ 231.27
```

Reusing `r_gas = 8.314462618` (`meds_constants.f90:24`) for `r_gas_kj` keeps provenance single-sourced.
(The `12.011` g/mol here vs `12.0107` in `umol_2_kgC` is Davidson's value; identical to 5 sig figs — use
`12.0107` if strict single-sourcing is preferred, the difference is < 3e-5 relative.) **Numerical guards
(GPU-safe, branchless):** the only fragility is `a_air^{4/3}` of a negative when `θ > θ_sat` from rounding
— killed by `max(θ_sat − θ, 0)`; the MM denominators are `kM + conc` with `kM > 0`, so they never divide
by zero even at `conc = 0`. DAMM adds only intrinsics (`exp`, two divides, two powers), so the kernel
stays `pure`/bare-scalar/OpenMP-target-eligible.

**Budget note.** Like the Q10/EXP_ED2 Rh, DAMM is a **source term** with no residual of its own — its
value enters `f_bio` / the bottom Neumann `q_ground`, whose column residual (§4a / §3.5.7) closes.

**Hand-verified value (integration check).** At the Harvard-Forest point (`T = 288.15 K`, `θ = 0.229`,
`θ_sat = 0.6825`, `fast_soil_carbon = 4.8 kgC/m²`, `depth = 10 cm`): `Vmax = 4.29e-3`, `Sx = 7.57e-7`,
`mm_sx = 0.432`, `a = 0.4535`, `O2 = 0.1216`, `mm_o2 = 0.501`, `R_Sx = 9.28e-4 mgC cm⁻³ h⁻¹` →
`Rh = 9.28e-4 × 10 × 231.269 = 2.146 µmol/m²/s`, matching the paper's fitted Q10 `Rref`. This is test #15.

#### 4c′.3 DAMM parameters `damm_params_t` (Davidson-2012 defaults) + TOML

Provenance tiers (per the no-hard-coded-params rule, §6): **physically fixed** (`d_liq`, `d_gas`,
`o2_air_frac`, exponents 3 & 4/3, `bd`/`pd` → porosity), **calibrated** (`alpha_sx`, `ea_sx` — strongly
correlated, log-fit), **weakly-identifiable fixed priors** (`km_sx`, `km_o2`). `bd`/`pd` are retained for
provenance only — they enter the runtime **solely** through `theta_sat` (the operative porosity, from
`soil_params_t`), so the kernel signature needs only `soil_temp`/`theta`/`theta_sat`/pool. Where a site's
`theta_sat ≠ 1 − bd/pd = 0.6825`, `d_gas`'s dry-soil calibration shifts slightly (documented, P2 recal slot).

```fortran
!----- DAMM diffusion parameters (Davidson et al. 2012 Table 2, Harvard Forest). ---------------!
type :: damm_params_t
   real(wp) :: alpha_sx  = 5.38e10_wp   !< [mgC cm-3 soil h-1] Arrhenius pre-exponential (CALIBRATED; absorbs depth)
   real(wp) :: ea_sx     = 72.26_wp     !< [kJ/mol]            activation energy          (CALIBRATED)
   real(wp) :: km_sx     = 9.95e-7_wp   !< [gC cm-3 soil]      soluble-C half-saturation  (weak prior)
   real(wp) :: km_o2     = 0.121_wp     !< [cm3 O2 cm-3 air]   O2 half-saturation         (weak prior)
   real(wp) :: p_soluble = 4.14e-4_wp   !< [-] soluble fraction of total soil C  (CONFIRMED 4.14e-4, NOT 0.0414)
   real(wp) :: d_liq     = 3.17_wp      !< [-] liquid-diffusion coefficient      (physically fixed)
   real(wp) :: d_gas     = 1.67_wp      !< [-] gas-diffusion coefficient         (physically fixed)
   real(wp) :: depth_cm  = 10.0_wp      !< [cm] effective respiring depth (SOC->conc AND flux depth-integration)
   real(wp) :: bd        = 0.80_wp      !< [g cm-3] bulk density   (provenance; enters only via theta_sat)
   real(wp) :: pd        = 2.52_wp      !< [g cm-3] particle density(provenance; 1 - bd/pd = 0.6825 = porosity)
end type damm_params_t
```

```toml
[co2]
hr_model = "damm"                # -> req_hr_model mapper: "q10" | "exp_ed2" | "damm"; case default = hard error

[co2.damm]                       # required only when hr_model = "damm" (presence-mapped)
alpha_sx  = 5.38e10              # [mgC cm-3 soil h-1]
ea_sx     = 72.26                # [kJ/mol]
km_sx     = 9.95e-7              # [gC cm-3 soil]
km_o2     = 0.121                # [cm3 O2 cm-3 air]
p_soluble = 4.14e-4              # [-]
d_liq     = 3.17                 # [-]
d_gas     = 1.67                 # [-]
depth_cm  = 10.0                 # [cm]
bd        = 0.80                 # [g cm-3]
pd        = 2.52                 # [g cm-3]
```

`validate_config` guards (DAMM branch): `alpha_sx > 0`, `ea_sx > 0`, `km_sx > 0`, `km_o2 > 0`,
`0 < p_soluble < 1`, `d_liq > 0`, `d_gas > 0`, `depth_cm > 0`, `0 < bd < pd`, and a **soft warn** if
`|theta_sat − (1 − bd/pd)| > 0.05` (the `d_gas`-recal flag). The `req_hr_model` string→enum mapper is
cloned from `req_temp_response` (§6), `case default → note_missing` so an unknown string is a hard error.
When `hr_model = "damm"`, `rh_k_base`/`resp_*` are unused; `validate_config` does **not** require
`rh_k_base > 0` in that mode (avoiding a false constraint on a would-be-ignored field).

#### 4c′.4 DAMM removes the empirical moisture-floor problem (a bonus over the Q10 path)

The §4c empirical branches normalize moisture against an **air-dry floor `theta_dry`** that is a standing
**correctness dependency** (§9): passing the van Genuchten residual `theta_res` where ED2's air-dry
`soilcp` (ψ ≈ −3.1 MPa) is meant silently mis-tunes the whole moisture response, and MVP must resolve it
to a real `soilcp`/retention-curve value at P1. **`HR_DAMM` deletes that dependency entirely.** The
moisture response is emergent from two diffusion terms:

- as `θ → 0`: `Sx = p·SOC·D_liq·θ³ → 0`, `mm_sx → 0` (substrate-diffusion limited),
- as `θ → θ_sat`: `a = θ_sat − θ → 0`, `O2 → 0`, `mm_o2 → 0` (O2/anoxia limited),

so `mm_sx(θ)·mm_o2(θ)` peaks at intermediate moisture **with no `resp_opt_water` / `resp_water_below_opt`
/ `resp_water_above_opt` / `theta_dry`**. DAMM needs only `soil_temp` (native Kelvin → Arrhenius
directly), `theta` and `theta_sat` (native `m³/m³`, no conversion; `theta_sat` from `soil_params_t`), and
`fast_soil_carbon` — **exactly the inputs the Q10 path already receives**, minus the problematic floor.
Documented limitation (not adopted blindly): the `θ³` / `(θ_sat−θ)^{4/3}` exponents are approximations and
the optimum's location is sensitive to `km_o2`; both exponents/diffusion coefs are TOML-exposed and the
hump shape is regression-guarded (§4c′.8 test #16).

#### 4c′.5 The generalizable soil-gas-flux seam (CH4/N2O stay additive)

DAMM's authors frame it as *"a module of a larger model"* that *"could be modified to include anaerobic
processes by substituting different electron acceptors in place of O2."* Every published generalization
(DAMM-MCNiP for C+N; M3D-DAMM for CH4; denitrification enzyme kinetics for N2O) keeps the identical
skeleton — **Arrhenius (or Q10) Vmax × product of Michaelis-Menten factors × diffusion-limited
concentrations** — and only changes *which pools feed the factors* and *how many pathways sum into a gas*.
So the MVP ships the concrete `heterotrophic_respiration_damm` (§4c′.2), and **reserves** (does not build)
a descriptor-driven `soil_gas_flux` kernel of which DAMM-CO2 is the first instance; CH4/N2O then add only
parameter rows + pool ids, with **no kernel change** (the MVP-honesty split: ship the concrete function
now, reserve the general seam, write no speculative CH4/N2O code):

```fortran
!----- Additive soil-gas-flux descriptor (RESERVED seam, P3+; fixed-size, no allocatables). ----!
integer(ik), parameter :: PHASE_AQUEOUS = 1_ik, PHASE_GASEOUS  = 2_ik    !< diffusion phase
integer(ik), parameter :: FACTOR_SUBSTRATE = 1_ik, FACTOR_INHIBITOR = 2_ik  !< MM (increasing) vs kI/(kI+C)
integer(ik), parameter :: max_gas_factor = 4_ik, max_gas_pathway = 3_ik
integer(ik), parameter :: GAS_CO2 = 1_ik, GAS_CH4 = 2_ik, GAS_N2O = 3_ik  !< planned species enum

type :: gas_factor_t
   integer(ik) :: phase = PHASE_AQUEOUS      !< aqueous (pool*D_liq*theta^expo) or gaseous (atm_frac*D_gas*a^expo)
   integer(ik) :: mode  = FACTOR_SUBSTRATE   !< increasing MM  or  decreasing inhibitor
   real(wp)    :: k_half = 0.0_wp            !< Michaelis half-saturation
   real(wp)    :: d_coef = 1.0_wp            !< diffusion coefficient (D_liq / D_gas)
   real(wp)    :: expo   = 3.0_wp            !< moisture exponent (3 aqueous, 4/3 gaseous)
   real(wp)    :: atm_frac = 0.0_wp          !< atmospheric volume fraction (gaseous; O2 = 0.209)
   integer(ik) :: pool_id  = 0_ik            !< index into the frozen substrate/pool array (aqueous)
end type gas_factor_t

type :: gas_pathway_t
   integer(ik) :: n_factor = 0_ik
   type(gas_factor_t) :: factor(max_gas_factor)
   real(wp)    :: alpha = 0.0_wp, ea = 0.0_wp, q10 = 0.0_wp, t_ref = 288.15_wp  !< Arrhenius OR Q10 Vmax
   real(wp)    :: yield = 1.0_wp             !< signed stoichiometry: +production, -consumption
end type gas_pathway_t

type :: soil_gas_kinetics_t
   integer(ik) :: n_pathway = 0_ik
   type(gas_pathway_t) :: pathway(max_gas_pathway)
   real(wp)    :: depth_eff = 0.10_wp        !< [m] effective respiring depth (or per-layer dz(k), §4c′.6)
end type soil_gas_kinetics_t
```

Two shared helpers hold the max-clamp and the Vmax form in one place each:
`diffusion_concentration(phase, bulk_or_atmfrac, d_coef, expo, theta, theta_sat)` (aqueous
`bulk*d_coef*theta**expo`; gaseous `atmfrac*d_coef*max(theta_sat-theta,0)**expo` — the clamp lives
**here, once**) and `vmax_arrhenius_or_q10(alpha, ea, q10, t_ref, temp)` (Arrhenius if `q10 <= 0`, else
`alpha*q10**((temp-t_ref)*0.1)`). The reserved kernel is one pure loop (pseudocode):

```
soil_gas_flux(kin, temp, theta, theta_sat, pool):
   flux = 0
   for p in 1..kin%n_pathway:
      rate = vmax_arrhenius_or_q10(pathway(p)%alpha, %ea, %q10, %t_ref, temp)
      for k in 1..pathway(p)%n_factor:
         conc = diffusion_concentration(f%phase, aqueous? pool(f%pool_id):0, f%d_coef, f%expo, theta, theta_sat, f%atm_frac)
         rate = rate * ( f%mode==SUBSTRATE ? conc/(f%k_half+conc) : f%k_half/(f%k_half+conc) )
      flux += pathway(p)%yield * rate
   return flux * kin%depth_eff        ! [umol m-2 s-1] (species-specific molar conversion at ingest)
```

**Instance mapping (all additive, ZERO kernel change):**

- **HR_DAMM CO2** = 1 pathway, 2 factors: `{AQUEOUS soluble-C substrate, pool_id=SOLUBLE_C, kM_Sx,
  D_liq=3.17, expo=3}` and `{GASEOUS O2 substrate, atm_frac=0.209, kM_O2=0.121, D_gas=1.67, expo=4/3}`;
  `yield=+1` — the first shipped instance (equivalent to §4c′.2's direct function).
- **CH4 (later)** = 3 pathways: aceticlastic `{AQUEOUS acetate + GASEOUS O2 INHIBITOR}`, hydrogenotrophic
  `{GASEOUS H2 + GASEOUS CO2 + O2 INHIBITOR}`, methanotrophy `{GASEOUS CH4 + GASEOUS O2, yield=−1}`; Q10
  form + optional pH factor.
- **N2O (later)** = nitrification `{AQUEOUS NH4 + GASEOUS O2}` + denitrification cascade
  `{AQUEOUS DOC donor + AQUEOUS N-oxide acceptor + O2 INHIBITOR}`; net N2O = production − reduction.

The common input vector is `{a substrate/acceptor concentration per MM factor, soil_temp, theta,
theta_sat}` — adding a gas adds `pool_id`s + pathway rows, not a signature. Pools stay
`intent(in)`/**frozen in the fast loop** (the fast/slow seam, §2.4/§9) — a future prognostic microbial/DOC
(Millennial/CENTURY) network updates them **daily** in `meds_soil_carbon.f90`; the fast `soil_gas_flux`
reads them, needs **no sub-stepping**, and stays GPU-eligible.

#### 4c′.6 Multi-layer DAMM — the P2 upgrade that parallels §3.5

DAMM is per unit soil **volume**, so the MVP's single `depth_cm` multiplier is the CAS box's `n = 1`
analog. The natural P2 upgrade — **paralleling §3.5's multi-layer canopy** — is to evaluate
`heterotrophic_respiration_damm` **per soil layer** over the *existing* negative-z `soil_params_t` grid,
each layer using its own `soil_temp(k)` (from `soil_energy_flux`) and `theta(k)` (from
`column_hydrology_flux`) and its own SOC share, summing production over the column:

```
heterotrophic_respiration = Σ_k  R_Sx( soil_temp(k), theta(k), theta_sat(k), SOC(k) ) * dz_cm(k) * damm_flux_factor
```

with `dz_cm(k) = dz(k)·100` replacing the scalar `depth_cm`. This reuses the soil-column machinery (same
`dz(k)`, same per-layer `theta`/`T` the energy/hydrology columns already produce) exactly as
`soil_energy_flux`/`column_hydrology_flux` do — the biogeochem sibling of §3.5's multi-layer CAS. It also
**restores** the per-layer root-respiration integral and the near-surface-vs-deep temperature split that
§9 flags as MVP biases (deeper layers wetter/colder/more O2-limited), replacing the single root-weighted
`T_soil` collapse. Structurally identical to §3.5: MVP = 1 effective depth; P2 = per-layer over the
fixed-size column; both close the same source→`f_bio`→twin residual.

#### 4c′.7 Config (§6 revision) + provenance

Add to §6: the `hr_model` mapper (`"q10"|"exp_ed2"|"damm"`), the `[co2.damm]` presence-mapped block
(required only when `hr_model = "damm"`), and the DAMM `validate_config` guards (§4c′.3). The
provenance-tier documentation (physically-fixed / calibrated / weakly-identifiable) goes in the
`[co2.damm]` comment header, and the three new `meds_constants` (`r_gas_kj`, `o2_air_frac`,
`damm_flux_factor`) join the §3.3 carbon-conversion block — `damm_flux_factor` **written as the
expression** `1.0e4/3600 * 1000/12.011` so the molar-mass choice is auditable.

#### 4c′.8 Test additions (§7)

Extend `test_column_co2` (built under **both ifx and nvfortran multicore**, CLAUDE.md issue #7):

15. **`damm_hand_value`** — at `(soil_temp = 288.15 K, theta = 0.229, theta_sat = 0.6825,
    fast_soil_carbon = 4.8 kgC/m2, depth_cm = 10)` assert `Rh ≈ 2.15 µmol/m²/s` (tol ~1%); intermediates
    per §4c′.2. Guards the whole unit chain and the `231.269` factor.
16. **`damm_moisture_unimodality`** — sweep `theta` from `0` to `theta_sat` at fixed `T`; assert `Rh`
    **rises then falls** with a single interior maximum, `Rh(0) ≈ 0` (`θ³→0`) and `Rh(theta_sat) ≈ 0`
    (`(θ_sat−θ)^{4/3}→0`), and **stays finite** (no NaN) at `theta = theta_sat`. The mechanistic-moisture
    regression guard that *replaces* the empirical-`f_water` branch test for the DAMM path.
17. **`damm_arrhenius`** — hold `theta` fixed, raise `soil_temp`; assert the `Rh` ratio equals the analytic
    `exp((ea_sx/r_gas_kj)·(1/T1 − 1/T2))` (guards the kJ/mol `Ea` ↔ `r_gas_kj` pairing — a units slip here
    is a silent large error).
18. **`damm_anoxia_limit`** — as `theta → theta_sat`, `mm_o2 → 0` so `Rh → 0` (O2 limitation), and the
    clamp holds `Rh = 0` (not NaN) for `theta` slightly above `theta_sat` from rounding.

**Phasing (§8):** `HR_DAMM` is a **P1** addition (selectable alongside the P0 Q10/EXP_ED2, still a bare
`pure` function on a frozen pool — no new state); multi-layer DAMM (§4c′.6) is **P2** beside the
multi-layer CAS (§3.5.10); the descriptor-driven `soil_gas_flux` seam + CH4/N2O parameter sets are
**P3+** (reserved, additive, no kernel change).

---

### 4d. `column_co2_step` — the assembler + NEE/NEP + closure check

The thin (non-`pure`) host-side seam that runs the whole per-patch column CO2 step: aggregate cohorts
(4b), attach the committed growth/storage respiration, compute Rh (4c), sum `plant_respiration`,
advance the twin (4a), and enforce the closure in Debug. This is the fast-loop coupling point at P3.

```fortran
subroutine column_co2_step(cas_can_co2, can_depth, can_shv, ustar, co2_atm, rho_air, dt, &
                           n, nplant, leaf_area, a_gross, rd, stem_resp, root_resp,       &
                           growth_resp_committed, storage_resp_committed,                 &
                           fast_soil_carbon, soil_temp, theta, theta_dry, theta_sat,      &
                           opts, budget)
   use meds_kinds,            only : wp, ik
   use meds_biogeochem_types, only : co2_opts_t, column_co2_budget_t, cohort_co2_flux_t
   real(wp),    intent(inout) :: cas_can_co2            ! [umol/mol] = cas%can_co2 (passed by reference)
   real(wp),    intent(in)    :: can_depth, can_shv, ustar, co2_atm, rho_air, dt
   integer(ik), intent(in)    :: n
   real(wp),    intent(in)    :: nplant(n), leaf_area(n), a_gross(n), rd(n), stem_resp(n), root_resp(n)
   real(wp),    intent(in)    :: growth_resp_committed, storage_resp_committed  ! [umol/m2/s] MVP = 0
   real(wp),    intent(in)    :: fast_soil_carbon, soil_temp, theta, theta_dry, theta_sat
   type(co2_opts_t),          intent(in)  :: opts
   type(column_co2_budget_t), intent(out) :: budget
   type(cohort_co2_flux_t) :: coh
   real(wp) :: hetero, plant_resp, scale
   call aggregate_cohort_co2_fluxes(n, nplant, leaf_area, a_gross, rd, stem_resp, root_resp, coh)
   coh%growth_respiration  = growth_resp_committed
   coh%storage_respiration = storage_resp_committed
   hetero     = heterotrophic_respiration_flux(fast_soil_carbon, soil_temp, theta, theta_dry, theta_sat, opts)
   plant_resp = coh%leaf_respiration + coh%stem_respiration + coh%root_respiration        &
              + coh%growth_respiration + coh%storage_respiration
   call canopy_air_co2_update(cas_can_co2, can_depth, can_shv, coh%gross_primary_prod, plant_resp, &
                              hetero, ustar, co2_atm, rho_air, dt, budget)
   !----- Closed-budget guard (the uniform biophysics discipline; mixed rtol/atol form). -------!
   scale = max(abs(budget%storage), 1.0_wp)
   if (opts%debug_error .and. abs(budget%resid) > opts%rtol*scale + opts%atol) then
      error stop 'column_co2_step: CO2 budget did not close'
   end if
end subroutine column_co2_step
```

**Tolerance convention (Lens 1 fix).** The guard uses the **standard mixed form**
`rtol*scale + atol` (`rtol` dimensionless, `atol [umol/m2]`), mirroring the relative-to-capacity
convention of `soil_energy_flux` (`atol * sum(c_eff*dz)`, `meds_column_energy.f90:113`). Both
tolerances are now used and unambiguously typed.

**Closure.** `column_co2_step` inherits 4a's `resid ≡ 0`; the `error stop` fires only on a *numerical*
fault (the algebra guarantees the physics). NEE/NEP/loss2atm are in `budget`; convert to `kgC` at the
I/O seam with `umol_2_kgC`.

---

## 5. Units + sign conventions

Every flux is per **ground** area; every store is per ground area; every kernel's `resid ~ 0`.

| Quantity | Symbol / field | Unit | Sign convention |
|---|---|---|---|
| CAS CO2 (state) | `can_co2` | `[umol/mol]` dry-air mixing ratio (= ppm) | — (prognostic) |
| Free-atm CO2 (forcing) | `co2_atm` | `[umol/mol]` | — |
| Molar CAS capacity | `ccapcan = can_dmol·can_depth` | `[mol_dryair/m2]` | > 0 |
| Dry-air molar density | `can_dmol = rho_air·(1−can_shv)/mmdry` | `[mol_dryair/m3]` | > 0 |
| Atm↔CAS molar conductance | `gatm_co2 = can_dmol·ustar` (c3 = 1) | `[mol_dryair/m2/s]` | > 0 |
| Gross primary production | `gross_primary_prod` | `[umol/m2/s]` | ≥ 0; enters `f_bio` with **−** (uptake removes CO2) |
| Autotrophic respiration | `plant_respiration` (leaf+stem+root+growth+storage) | `[umol/m2/s]` | ≥ 0; enters `f_bio` with **+** |
| Heterotrophic respiration | `heterotrophic_respiration` | `[umol/m2/s]` | ≥ 0; enters `f_bio` with **+** |
| Net biotic source | `f_bio = Reco − GPP` | `[umol/m2/s]` | source-positive |
| **NEE** | `budget%nee = f_bio` | `[umol/m2/s]` | **> 0 = net source to atmosphere** |
| **NEP** (derived) | `budget%nep = −f_bio` | `[umol/m2/s]` | **> 0 = net ecosystem uptake** |
| Flux to free atmosphere | `budget%loss2atm = gatm_co2·(can_co2−co2_atm)` | `[umol/m2/s]` | **> 0 = CAS vents upward** |
| CAS CO2 storage | `budget%storage = ccapcan·can_co2` | `[umol/m2]` | > 0 |
| **Closure residual** | `budget%resid` | `[umol/m2]` | **≡ 0** by construction |
| Soil-carbon pool (slow) | `fast_soil_carbon` | `[kgC/m2]` | ≥ 0; `intent(in)` in the fast loop |
| Decomposition rate | `rh_k_base` | `[1/day]` | > 0 |
| Temperature/moisture modifiers | `f_temp`, `f_water` | `[-]` | — (`f_water ∈ [0,1]`; `f_temp ≤ 1` for `HR_EXP_ED2`) |
| Carbon-mass output | `umol_2_kgC` | `[kgC/umol]` | applied at I/O only |

---

## 6. Config: `[co2]` TOML block + `co2_opts_t`, and the provenance rule

**Two-part idiom, matching `soil_opts_t` / `energy_opts_t`:** the standalone `co2_opts_t` (§3.2)
carries in-type defaults *only* so the kernels and tests compile pre-P3; the run values are
**REQUIRED from TOML** and presence-mapped at P3 (no hard-coded model parameters in `meds_config_t` —
the MEDS provenance rule).

**P3 wiring into `meds_config_io%load_meds_config`** (the pattern exists — `req_r` / `req_i` / `req_l`,
`note_missing`, and the string→enum mappers `req_growth_source` / `req_temp_response`):

```toml
[co2]
hr_model              = "q10"        # -> req_hr_model mapper: "q10"|"exp_ed2"|"damm"; case default = hard error
                                     #    ("damm" additionally requires the [co2.damm] block, §4c′.3)
rh_k_base             = 0.0301       # [1/day] effective decomposition rate (e.g. 11/yr for an ED2 fast pool)
rh_q10                = 1.5          # [-]  Collatz/K13
rh_t_ref              = 288.15       # [K]  15 C
resp_temp_increase    = 0.0757       # [1/K] ED2 scheme-0
resp_temp_ref         = 318.15       # [K]  45 C saturation
resp_opt_water        = 0.8938
resp_water_below_opt  = 5.0786
resp_water_above_opt  = 4.5139
co2_atm_ref           = 400.0        # [umol/mol] fixed atm CO2 when not met-forced
rtol                  = 1.0e-8
atol                  = 1.0e-3       # [umol/m2] absolute closure floor
```

- Scalars via `call req_r(tm,'co2.rh_k_base',cfg%co2%rh_k_base,miss)` etc.; the selector via a new
  `req_hr_model` mapper **cloned from `req_temp_response`** — default value, `note_missing` if
  absent, `select case` with `case default → note_missing` so an unknown string is a hard error. When
  `hr_model = "damm"`, presence-map the `[co2.damm]` block (§4c′.3) and skip the empirical `rh_*`/`resp_*`
  requirements.
- **GPP provenance** is already handled: `cfg%gpp_ref [kgC/m2 leaf/yr]` (carbon mode) feeds the stub
  GPP until canopy RT supplies absorbed-PAR-driven `a_gross`; the P1 stub applies `umols_per_kgCyr`
  (§3.3) for the `[kgC/m2 leaf/yr] → [umol/m2 leaf/s]` conversion.
- **Per-PFT decomposition/respiration traits** (when the CENTURY pool lands) go in the PFT file via
  `req_pa('pft.<trait>', …, npft, miss)`.
- **`validate_config` guards:** `rh_q10 > 0`, `rh_k_base ≥ 0` (see the stub-zeroing caveat below),
  `co2_atm_ref > 0`, `0 < resp_opt_water < 1`, `theta_sat > theta_dry ≥ 0`, `rtol ≥ 0`, `atol > 0`.
  Any derived quantity goes in `derive_config`.
- **Stub-zeroing caveat (Lens 1).** `rh_k_base = 0` (the type default) makes Rh ≡ 0 silently. Keep the
  loose `rh_k_base ≥ 0` guard (a legitimately GPP-only test may want it), but defend against a
  mis-wire two ways: (1) test 8 asserts that a **nonzero** `(fast_soil_carbon, rh_k_base)` yields a
  nonzero Rh; (2) when a soil-carbon source is enabled in config at P3, `validate_config` requires
  `rh_k_base > 0`.

---

## 7. Test plan (`test/test_column_co2.f90`)

Program `test_column_co2` in the house style (`check` / `check_true`, `nfail`, `error stop 1`),
linking `meds_biogeochemistry`. **Build under both ifx *and* nvfortran multicore** — a green ifx run
is not sufficient (CLAUDE.md issue #7; the array-temporary miscompile is silent at `-O2`).

1. **`co2_resid_zero`** — randomized `can_co2`, fluxes, `ustar`, `rho_air`, `dt`: assert
   `|budget%resid| < 1e-10·max(budget%storage,1)`. The core conservation guarantee.
2. **`co2_steady_state`** — `f_bio = 0`, `can_co2 = co2_atm`: `can_co2` unchanged, `loss2atm = 0`.
3. **`co2_atm_relaxation`** — `f_bio = 0`, `can_co2 ≠ co2_atm`: monotone relaxation toward `co2_atm`;
   **L-stable for large `dt·cci·gatm_co2`** (no overshoot, never negative — the payoff of the implicit
   form).
4. **`co2_steady_ca`** — constant `f_bio`, iterate to fixed point:
   `can_co2 → co2_atm + f_bio/gatm_co2` (the analytic steady canopy CO2, ED2's `cstar` balance).
5. **`co2_sign_discipline`** — GPP-only pulls `can_co2` below `co2_atm` (`nep > 0`); respiration-only
   pushes it above (`nee > 0`, `loss2atm > 0` = nocturnal venting).
6. **`co2_conservation_identity`** — recompute storage independently as `ccapcan·can_co2` and assert
   `Δstorage = dt·(nee − loss2atm)` to machine precision, using **independently assembled** `nee` /
   `loss2atm` (catches a mis-assembled `f_bio` that the kernel's internal algebra would mask — §9).
7. **`aggregate_units`** — a hand-built 3-cohort patch: GPP/`Rd` leaf-area-weighted, stem/root
   `nplant`-weighted; matches the analytic µmol/m²/s sum (guards the LAI weighting).
8. **`heterotrophic_response`** — Rh ≈ doubles per 10 K at `rh_q10 = 2`; `f_temp ≤ 1` and saturates at
   45 °C for `HR_EXP_ED2`; `f_water` peaks at `resp_opt_water` and falls both sides; `Rh = 0` at
   `fast_soil_carbon = 0`; **`Rh > 0` for nonzero `(C, k)`** (the stub-zeroing guard); a known
   `(C,k,T,θ,theta_dry,theta_sat)` reproduces the hand-computed µmol/m²/s (guards the ~963.6 factor —
   the silent ~1000× trap — and the air-dry floor).
9. **`nep_identity`** — `nep = −nee = GPP − Reco` exactly.
10. **`co2_multistep_density` (expected-fail scaffold)** — vary `can_shv` between substeps and assert
    the *multi-step* budget does **not** reconcile without a `denseffect` term (documents risk 2 and
    prevents the P2 density correction from being silently skipped). Flipped to a passing test when the
    P2 `denseffect` / `zcaneffect` terms land.
11. **`nvfortran_multicore_build`** — CI target: the new modules compile and tests pass under
    `nvfortran -mp` (not just ifx).

---

## 8. Phasing P0 → P3 and the fast-loop coupling seam

Consistent with the biophysics P0/P1/P2/P3 phasing (RT, hydrology, energy all followed this):

- **P0 — stateless kernels + tests (the bulk of the value, ~standalone).** `meds_biogeochem_types`
  (types + `HR_*`), `meds_column_co2` (all four kernels), the five `meds_constants` additions, the two
  one-line `cas_state_t` / `cas_atm_forcing_t` edits, the `leaf_flux_t` comment annotation, the CMake
  `meds_biogeochemistry` library + `test_column_co2` wiring, the README/`CLAUDE.md` charter amendment
  (§2.4), tests 1–11. Soil-C pool is a **prescribed scalar input** (`fast_soil_carbon`, no state
  update). Everything bare-scalar, `pure`, GPU-eligible. No per-patch state, no config wiring — exactly
  how soil hydrology/energy landed.
- **P1 — aggregation seam + diagnostics.** Wire the host-side cohort reduction over the CSR slice
  (`cohort_offset` / `cohort_count`) using the stub GPP (`gpp_ref × umols_per_kgCyr`) +
  `meds_plant_respiration` sums; resolve `theta_dry` (air-dry floor, §4c); expose
  `nee` / `nep` / `loss2atm` / `storage` through `meds_io` (convert to `kgC` with `umol_2_kgC`). Still
  stateless (bare arrays in, budget out).
- **P2 — fidelity upgrades, still stateless.** Slow `meds_soil_carbon%soil_carbon_step` (daily: litter
  in from demography turnover/deficit − decomposition out) — the fast Rh now reads the slowly-updated
  `fast_soil_carbon`; the CENTURY multi-pool (fast/structural/slow, lignin, N-immobilization, respired
  fractions, A/B above/below split) begins here. Per-layer root-weighted soil temperature/moisture for
  Rh, and the per-layer root-respiration integral (replacing the single-`T_soil` collapse). The
  **`denseffect` / `zcaneffect`** residual terms (`can_depth·(dmol_new − dmol_old)·can_co2`) so the
  *multi-step* budget stays closed once `rho_air` / `can_shv` / `can_depth` vary between steps — gated
  in by flipping test 10 to passing. A **shared `c3`-corrected conductance** for all three twins
  (§1.5). Optional `q10_scale` promotion to `meds_temp_response`.
- **P3 — state + config + coupling (the deferred integration all biophysics stores share).**
  - **Per-patch state lockstep:** `can_co2` joins the patch CAS SoA; the slow `soil_carbon_t` joins the
    patch state. Thread both through the patch permute/pack sites `sort_patches` **and**
    `patch_compact` (patches have no single reorder routine — CLAUDE.md), `patch_ensure_capacity`, the
    `fuse_2_patches` **area-fraction blend**, and the disturbance **donor-copy**, in lockstep beside
    `can_enthalpy` / `can_shv` (the ED2 "forgot to reallocate" trap).
  - **Config:** the `[co2]` block into `meds_config_t` / `meds_config_io` / `validate_config` /
    `derive_config`; add `meds_biogeochemistry` to the `meds_aux` link line.
  - **Fast-loop master step (the seam).** In the future fast loop (`meds_stepper`, the ED2 `DTLSM`
    analogue), per patch per substep: canopy RT → `meds_leaf_physiology%leaf_gas_exchange` →
    `canopy_air_update` (enthalpy/shv) → `column_hydrology_flux` → **`column_co2_step`** (the CO2
    twin), all sharing one `ustar` / `rho_air` / `can_depth`, and one `can_shv` **taken after the
    energy twin's within-substep update** (so `can_dmol` sees the freshest humidity; irrelevant to the
    single-step residual, which is internally consistent, but it fixes the ordering for reproducibility
    — §9). **The feedback that makes the twin matter:** set `leaf_env_t%ca = cas%can_co2` (currently
    the atm reference) for the next photosynthesis solve — the single substitution that turns a
    CLM/CliMA-style forced-CO2 leaf model into an ED2-style coupled one (nocturnal sub-canopy buildup
    feeds back on `cs` / `ci`). The slow soil-C pool updates once per day from demography litter,
    closing the ground↔CAS↔leaf carbon loop.

---

## 9. Open questions, risks, deferred work

**The two biggest risks — both *outside* the kernel's self-certifying algebra:**

1. **`resid ≡ 0` gives false confidence at the fast/slow seam.** The kernel's residual certifies only
   the twin's *internal* algebra — it is zero **regardless of whether `f_bio` was assembled
   correctly**. The real failure surface is where `f_bio` is built: a forgotten `kgCday_2_umols` /
   `umols_per_kgCyr` / `nplant` / `leaf_area` factor (silent ~1000× or LAI error), a sign slip (GPP
   must be **−**; a flip inverts NEE), or **wet `rho_air`** where the molar capacity needs dry-air
   density (§1.4 — a silent ~1–2 % mixing-ratio drift per strong diurnal cycle). Each produces a
   physically wrong `can_co2` while `resid` reads a contented ~0. **Mitigation:** test 6 recomputes
   storage *independently* and reconciles `Δstorage` against `dt·(nee − loss2atm)` assembled from
   independent inputs — a check the kernel's algebra cannot mask; tests 7–8 pin the unit factors; the
   dry-air `(1 − can_shv)` factor is in the kernel from day one.

2. **Variable-density budget drift at P3.** Once `rho_air` / `can_shv` / `can_depth` vary between
   substeps, storage changes with no flux and the *multi-step* budget will not reconcile unless the
   `denseffect` / `zcaneffect` terms are added. **Mitigation:** hold `can_dmol` / `can_depth` fixed
   *within* a step at MVP (making both corrections identically zero and the single-step residual
   exact); test 10 is a standing expected-fail that *forces* the P2 density correction before variable
   density is enabled in the coupled loop.

**Cadence risk (the fast/slow wall around the pool).** Correctness of the whole budget hinges on
`fast_soil_carbon` being **read-only within the fast loop and written exactly once per day** (ED2's
exploited seam). If a contributor lets `heterotrophic_respiration_flux` mutate the pool, or
double-applies daily decomposition (once slowly *and* implicitly through fast Rh), the per-substep CO2
`resid` still reads ~0 while carbon-mass conservation across the fast/slow boundary silently drifts.
**Mitigation:** `fast_soil_carbon` is `intent(in)` in every fast kernel (enforced by the `pure`
signatures); add a **daily carbon-mass audit** (`ΔC_pool = litter_in − ∫Rh·dt`) as a separate P2 test,
distinct from the per-substep residual.

**Known MVP biases (documented, not bugs — do not compare to observed NEE until P3):**
- **Committed growth + storage respiration = 0** (§4b). ED2's `nee_tot` includes both
  (`rk4_derivs.f90:1394–1408`); the MVP omits them, biasing NEE toward uptake. Slots reserved in
  `cohort_co2_flux_t`; turned on when carbon allocation couples to the fast loop (P3).
- **Stability coefficient `c3` = 1** (§1.5). Nocturnal CO2 accumulation is qualitatively present but
  quantitatively damped vs ED2; a shared `c3` for all three twins is a P2 item.
- **Air-dry moisture floor `theta_dry`** (§4c) is a **correctness dependency**, not a free choice:
  resolve it at P1 to a real `soilcp` / retention-curve air-dry value, never silently to `theta_res`.
- **Root respiration single-layer temperature collapse** (§4c): MVP uses one root-weighted `T_soil`;
  a near-surface-weighted temperature inflates the diurnal amplitude of root respiration vs ED2's
  per-layer integral (`soil_respiration.f90:190`). Restored at P2.

**Deferred work.** Full CENTURY biogeochemistry (multi-pool dynamics, N cycle, lignin, inter-pool
transfers, litter routing from the demography turnover/deficit streams that are currently dropped) is
a whole slow-loop phase, explicitly out of scope for this fast CO2 module. GPU-eligibility of the CO2
kernels is asserted by construction (bare-scalar, `pure`, intrinsics only) but is only *tested* under
`nvfortran -mp` at P0 (test 11); actual `!$omp target` offload of the fast loop is P3.

---

## 10. File manifest (all absolute)

**New (P0/P1 — the MVP scalar box + DAMM option).**
- `/home/xiangtao/projects/MEDS/src/biogeochemistry/meds_biogeochem_types.f90` (incl. `damm_params_t`, `HR_DAMM`)
- `/home/xiangtao/projects/MEDS/src/biogeochemistry/meds_column_co2.f90` (incl. `heterotrophic_respiration_damm`)
- `/home/xiangtao/projects/MEDS/src/biogeochemistry/meds_soil_carbon.f90` (P2 stub)
- `/home/xiangtao/projects/MEDS/test/test_column_co2.f90`

**New (P2 — multi-layer CAS forward-design, §3.5; shared with energy/water).**
- `/home/xiangtao/projects/MEDS/src/shared/meds_column_solver.f90` — `thomas_solve` (promoted from
  `meds_soil_solver`, rewritten assumed-size) + `n_column_layer_max`/`n_can_layer_max`.
- `/home/xiangtao/projects/MEDS/src/shared/meds_cas_transport.f90` — `canopy_air_grid_t`,
  `init_canopy_air_grid`, `eddy_diffusivity_profile`, the shared `canopy_air_be_step`.
- (edit) `meds_column_co2.f90` gains `canopy_air_co2_column_update`; `meds_column_energy.f90` /
  `meds_column_hydrology.f90` switch their `use meds_soil_solver` → `meds_column_solver`.

**Edit (one line each).**
- `/home/xiangtao/projects/MEDS/src/biophysics/meds_biophysics_types.f90` — `can_co2` on
  `cas_state_t:223` (scalar at P0/P1; array `can_co2(n_can_layer_max)` at P2, §3.5.4), `co2_atm` on
  `cas_atm_forcing_t:259`.
- `/home/xiangtao/projects/MEDS/src/plant/meds_plant_types.f90` — annotate `leaf_flux_t:63-69`
  comments as `[umol CO2/m2 leaf/s]`.

**Edit (constants / build / docs).**
- `/home/xiangtao/projects/MEDS/src/shared/meds_constants.f90` — add `mmdry`, `umol_2_kgC`,
  `kgC_2_umol`, `kgCday_2_umols`, `umols_per_kgCyr` after line 30; and (DAMM, §4c′.2) `r_gas_kj`,
  `o2_air_frac`, `damm_flux_factor`.
- `/home/xiangtao/projects/MEDS/CMakeLists.txt` — new `meds_biogeochemistry` library after the
  biophysics block (~line 100); `test_column_co2` after the energy tests (~line 230).
- `/home/xiangtao/projects/MEDS/src/biogeochemistry/README.md` and
  `/home/xiangtao/projects/MEDS/CLAUDE.md` — charter amendment (§2.4): biogeochemistry owns fast CO2
  exchange *and* slow carbon pools; add it to the library-DAG line as a shared-only sibling of
  biophysics.

**Reference kernel to mirror.**
- `/home/xiangtao/projects/MEDS/src/biophysics/meds_column_energy.f90:253-277` (`canopy_air_update`).

**ED2 references.**
- `/home/xiangtao/projects/ED2/ED/src/dynamics/rk4_derivs.f90` (`can_co2` ODE :2163, `nee_tot`
  :1394–1408, `co2budget` :2236).
- `/home/xiangtao/projects/ED2/ED/src/dynamics/canopy_struct_dynamics.f90` (`ccapcan = can_dmol·can_depth`
  :2181).
- `/home/xiangtao/projects/ED2/ED/src/dynamics/soil_respiration.f90` (`het_resp_weight` :1486).
