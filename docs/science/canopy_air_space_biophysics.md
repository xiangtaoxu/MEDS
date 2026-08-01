# Canopy-air-space biophysics

The **canopy air space (CAS)** is the well-mixed sub-canopy air box that acts as the shared coupling
reservoir of the fast loop: every surface store (leaf, wood, ground, snow) exchanges heat, water vapour,
and CO₂ with it, and it in turn vents to the free atmosphere. It carries three **implicit prognostic
twins** — specific enthalpy, specific humidity, and CO₂ mixing ratio — each advanced from its summed
surface sources and an implicit atmospheric-exchange term that shares the
$`u_*\cdot\mathrm{temp1/temp2}`$ conductance from the aerodynamics kernel
([canopy_aerodynamics](canopy_aerodynamics.md)).

This page is one of the per-store pages under [column_biophysics](column_biophysics.md), which describes
the two cross-cutting principles (prognostic internal energy; closed-budget discipline) and how the
stores are woven together each `dt_fast`.

---

## The three CAS twins

The CAS is a well-mixed box of air of depth $D_{can}$ and mass per ground area
$`W_{cap}=\rho\,D_{can}`$.

**$D_{can}$ is the tallest cohort's height plus a freeboard** (`[aerodynamics].canopy_freeboard`,
default 5 m, floored by `min_canopy_depth`): the canopy air space is the well-mixed layer the canopy
exchanges with, which extends *above* the crowns — it is not the canopy volume. It is **per-patch
state owned by the slow loop** (`refresh_canopy_depth`), recomputed after growth, mortality,
recruitment, fusion and disturbance have settled, so the fast loop sees a constant box across every
sub-step of a day.

Resizing that box is an **open-control-volume** operation, and the invariant is worth stating: growing
the canopy does not create air, it entrains air from just above at essentially the canopy-air state, so
the **intensive** state (specific enthalpy, specific humidity, CO₂ mixing ratio) is what carries over
unchanged; the extensive content changes and `cas_set_depth` reports that exchange. Conserving *total*
energy instead — rescaling `can_enthalpy` by the mass ratio — would cool the canopy air simply because
the trees grew.

> $`W_{cap}`$ is not a bookkeeping detail: it is the canopy air's heat capacity, so it sets how fast
> the CAS responds. It was a **hardcoded 20 m with no writer** until 2026-07-31 — the aerodynamics
> computed the right value and it was discarded — which meant a 1 m regenerating gap and a 35 m
> tropical canopy were given the same box.
>
> $`W_{cap}`$ used to set a hard *stability* bound on `dt_fast`, because it is the damping term
> opposing a lagged ventilation feedback. That feedback is no longer lagged (the surface layer is
> re-solved per stage, [canopy_aerodynamics](canopy_aerodynamics.md) §2), so the bound is gone.
> Measured across stand heights, the worst case was in any case a *mid-height* canopy rather than a
> short one — leaf area falls with height faster than canopy-air volume does at the low end — so the
> intuition that a regenerating gap is always the strictest case is wrong. See
> [numerical_scheme](numerical_scheme.md) §5a′.

All three twins (enthalpy, humidity, CO₂) are advanced by the shared
two-form box kernel `cas_column_step_implicit` / `cas_column_time_deriv` (`meds_cas_biophysics`),
implicit in the atmosphere exchange. Each twin obeys

```math
\mathcal{C}\,\frac{dX}{dt} = F_{surf} + g_{atm}\,(X_{atm}-X) \qquad(1)
```

with capacity $\mathcal{C}$, summed surface source $F_{surf}$, and the atm↔CAS conductance $g_{atm}$.
The atmospheric term is taken **implicit** for L-stability (it is the stiff coupling); the surface
sources are explicit. For the enthalpy twin, with $`F_{sens}`$ the summed cohort+ground sensible+latent
flux and $`g_{atm}=\rho\,u_*\,\mathrm{temp1}`$:

```math
H^{n+1} = \frac{W_{cap}H^n + \Delta t\,(F_{sens} + g_{atm}\,H_{atm})}{W_{cap} + \Delta t\,g_{atm}},
\qquad
\text{resid} = W_{cap}(H^{n+1}-H^n) - \Delta t\big(F_{sens}+g_{atm}(H_{atm}-H^{n+1})\big) = 0 \qquad(2)
```

The **vapour** twin is identical with $`g_{aw}=\rho u_*\,\mathrm{temp2}`$ and a non-negativity clamp.
The **CO₂** twin differs only in units: CO₂ is a molar mixing ratio, so the capacity is the dry-air
**molar** column $`C_{cap}=\rho\,(1-q)/M_{d}\cdot D_{can}`$ [mol m⁻²] and the conductance is
$`g_{ac}=\rho_{mol}\,u_*\,\mathrm{temp2}`$. Its biotic source is
$`F_{bio}=R_a+R_h-\mathrm{GPP}`$ (net ecosystem exchange), assembled by the driver (`column_fast_step`)
from the aggregated cohort GPP/leaf-respiration ($`\times`$ LAI), stem/root maintenance respiration
($`\times n_{plant}`$), and the heterotrophic soil flux (see [soil_biophysics](soil_biophysics.md#soil-heterotrophic-respiration)).
CAS temperature is re-diagnosed from $(H^{n+1},q^{n+1})$ each step. Sharing $`u_*`$ **and** the profile
factor across all three twins keeps them on one turbulence basis — the fix for the nocturnal-CO₂
over-coupling bug (see [canopy_aerodynamics](canopy_aerodynamics.md)).

The two forms of the box kernel are the two-integrator seam: `cas_column_step_implicit` is the
backward-Euler advance; `cas_column_time_deriv` is the explicit RHS the ESDIRK2 (`ark`) stepper calls. Both take the pre-assembled sources (`cas_source_t`) and the
capacities/conductances/atmospheric BCs (`cas_column_t`) — any scheme-specific adjustment (the ARK
condensation sink, snow sublimation) is folded into the source by the driver, so the box math
itself is scheme-shared.

---

## Prognostic state

| Store | Type | Prognostic variable(s) | Diagnosed |
|---|---|---|---|
| Canopy air space | `cas_state_t` | `can_enthalpy`, `can_shv`, `can_co2` | `can_temp` |

## Code map

| Concept | Routine |
|---|---|
| CAS enthalpy + vapour + CO₂ twins (two-form box) | `meds_cas_biophysics`: `cas_column_step_implicit`, `cas_column_time_deriv` |
| source / capacity / conductance bundles | `meds_cas_biophysics`: `cas_source_t`, `cas_column_t` |
| $`u_*`$ conductance (temp1/temp2) | `meds_canopy_aerodynamics` (see [canopy_aerodynamics](canopy_aerodynamics.md)) |

## References
- ED2 `../ED2/ED/src/dynamics/canopy_struct_dynamics.f90` — the canopy-air budget in the reference model.
- Design doc: `MEDS_COLUMN_CO2_BALANCE_DESIGN.md` (the CO₂ twin), `MEDS_COLUMN_DYNAMICS_DESIGN.md`.
