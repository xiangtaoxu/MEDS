# archive

Design records and planning documents kept for provenance (tracked in git). These are
point-in-time design rationales, not living documentation — the authoritative description of
the code lives in `CLAUDE.md`, the per-directory `README.md` files, and the source itself.

- `radiative_transfer_design.md` — design & implementation plan for the `src/biophys/` canopy
  radiative-transfer scheme: a faithful reimplementation of ED2's two-stream (`icanrad=2`) RT,
  modernized (unified multi-band solver, SCOPE/4SAIL leaf-angle scattering with a Beta leaf-angle
  distribution, in-house block-tridiagonal solver). Includes the ED2 bugs found during the port.
