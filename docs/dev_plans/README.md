# archive

Design records and planning documents kept for provenance (tracked in git). These are
point-in-time design rationales, not living documentation — the authoritative description of
the code lives in `CLAUDE.md`, the per-directory `README.md` files, and the source itself.

## `archive/` — retired, not merely dated

Everything in this directory is a point-in-time record, so being *superseded* is normal and does not
warrant moving a file. `archive/` is narrower: it holds documents whose **conclusions are actively
misleading**, not just old — typically because they were measured against a reference that turned out
to be invalid. Each carries a tombstone at the top saying what survives and what does not.

- `archive/MEDS_INTEGRATOR_PARITY.md` — **RETIRED 2026-07-31.** Its central comparison scored `ark`
  and `rk45` against `split`, and split converged to a different, never-attributed limit (~0.45 K in
  canopy-air temperature) before being retired outright. Independently, every measurement in it was
  taken at `dt_fast` 900–1800 s, inside the freeze-cadence instability found later. Its *structural*
  content survives and is still cited from source comments (the Class 1/2/3 inventory and its row
  numbering, Phase A telemetry, the RK45-is-a-hybrid finding, the sapflow advected-enthalpy fix);
  none of its cross-scheme accuracy numbers do.

---

- `radiative_transfer_design.md` — design & implementation plan for the `src/biophys/` canopy
  radiative-transfer scheme: a faithful reimplementation of ED2's two-stream (`icanrad=2`) RT,
  modernized (unified multi-band solver, SCOPE/4SAIL leaf-angle scattering with a Beta leaf-angle
  distribution, in-house block-tridiagonal solver). Includes the ED2 bugs found during the port.
