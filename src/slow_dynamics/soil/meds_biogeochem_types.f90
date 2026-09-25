! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_biogeochem_types -- shared derived types + selector codes for the biogeochemistry     !
! domain: the SLOW soil-carbon / nutrient cycle of the ecosystem column. Pure DATA +           !
! parameters, no methods -- the biogeochem analogue of meds_canopy_types/meds_soil_types. Links              !
! src/shared (meds_kinds, meds_column_state_types) ONLY.                                          !
!                                                                                          !
! (The FAST canopy-air-space CO2 exchange -- a sub-daily biophysical diffusion/venting process --  !
! moved to meds_canopy_types/meds_soil_types + meds_cas_biophysics under src/biophysics; this module now holds      !
! only the SLOW soil-carbon pools.)                                                                 !
!                                                                                          !
! SLOW soil carbon (design MEDS_BIOGEOCHEMISTRY_DESIGN.md): the CENTURY-family multi-pool soil-    !
! carbon state advanced daily by the carbon matrix ODE dX/dt = B*I + A*xi*K*X. `soil_carbon_t` (the   !
! 7-pool vector + lignin sub-state + optional N) is DEFINED in meds_column_state_types (state/column, !
! docs/dev_plans/archive/MEDS_SLOW_DYNAMICS_DESIGN.md Part I §8.1) and RE-EXPORTED here; field                !
! `fast_soil_carbon` KEEPS its name/index (2) so meds_cas_biophysics -- which reads it as a BARE       !
! SCALAR in the fast loop -- compiles unchanged. `decomp_opts_t` + its selector codes are similarly    !
! DEFINED in meds_biogeochem_opts (shared/config, MEDS_SLOW_DYNAMICS_DESIGN.md Part II B0) so           !
! meds_config can carry it with no shared->biogeochemistry edge, and RE-EXPORTED here. Slow types      !
! owned HERE: litter_input_t / soilc_audit_t / soilc_diag_t, and the fixed pool count `n_soil_pool` +   !
! `IP_*`.                                                                                          !
!==========================================================================================!
module meds_biogeochem_types
   use meds_kinds, only : wp, ik
   use meds_column_state_types, only : litter_input_t
   implicit none
   private

   !----- Slow soil-carbon matrix additions (P0). ------------------------------------------------!
   public :: n_soil_pool
   public :: IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL, IP_MICR, IP_SLOW, IP_PASSIVE
   public :: litter_input_t, soilc_audit_t, soilc_seam_t, soilc_diag_t
   !----- Fast heterotrophic-respiration selectors + params (kernels in meds_soil_biogeochem). ----!

   !----- There is exactly ONE production heterotrophic-respiration authority: the CENTURY matrix   !
   !      (heterotrophic_respiration_matrix), reached through patch_heterotrophic_respiration.       !
   !                                                                                          !
   !      A Q10 form, an ED2 capped-exponential form, a DAMM mechanistic form, their HR_* selector   !
   !      codes, co2_opts_t and damm_params_t all used to live here. None was reachable: there was    !
   !      no `hr_model` TOML key anywhere and co2_opts_t was never carried by meds_config, so the      !
   !      selector could not be set and the kernels could not be called. They were exercised only by   !
   !      their own unit tests -- the worst of both worlds, carrying maintenance and review weight     !
   !      while contributing nothing, and reading to a newcomer as available options (#153).           !
   !                                                                                          !
   !      Deleted 2026-09-13. The implementations are preserved on branch `archive/damm-hr` (at        !
   !      2fcb647) for whoever wants DAMM later; recovering them means reinstating the selector and    !
   !      the fast/slow rh_seam_gap contract for all three forms, which is the work that was never     !
   !      done. -------------------------------------------------------------------------------------!

   !==========================================================================================!
   !  SLOW soil-carbon matrix: fixed pool count + index parameters (the matrix state vector X).   !
   !  Ordered litter -> SOM -> passive so the default-scheme ACTIVE block of A is lower-triangular  !
   !  (back-substitution / GE solvable). The 7-pool ceiling matches ED2's CENTURY set exactly.      !
   !==========================================================================================!
   integer(ik), parameter :: n_soil_pool = 7_ik   !< length of the carbon matrix state vector X

   integer(ik), parameter :: IP_FAST_GRND   = 1_ik   !< metabolic litter, above-ground (flammable)
   integer(ik), parameter :: IP_FAST_SOIL   = 2_ik   !< metabolic litter, below-ground  (the P0 lumped pool)
   integer(ik), parameter :: IP_STRUCT_GRND = 3_ik   !< structural litter + CWD, above-ground
   integer(ik), parameter :: IP_STRUCT_SOIL = 4_ik   !< structural litter + CWD, below-ground
   integer(ik), parameter :: IP_MICR        = 5_ik   !< microbial SOM        (active only in scheme 5)
   integer(ik), parameter :: IP_SLOW        = 6_ik   !< slow / humified SOM
   integer(ik), parameter :: IP_PASSIVE     = 7_ik   !< passive SOM          (active only in scheme 5)

   !==========================================================================================!
   !  Slow, stateful per-patch soil-carbon pools (written DAILY, read-only in the fast loop).       !
   !  DEFINED in meds_column_state_types (src/state/column) -- re-exported here (its conceptual      !
   !  home) so meds_soil_biogeochem's kernels + meds_cas_biophysics compile unchanged. It lives in    !
   !  state/column, not here, so `patch_block` can carry a per-patch                                  !
   !  soil_carbon_t with no core->biogeochemistry library edge (the same reason cas_state_t/           !
   !  soil_column_t/soil_energy_column_t/snow_column_t live there). The ACTIVE pool count is set by     !
   !  decomp_scheme (default 3-active: idx 5,7 inert).                                                  !
   !==========================================================================================!

   !----- litter_input_t now lives in meds_column_state_types, beside the soil_carbon_t pools it
   !      feeds. It moved because it became a PER-PATCH FIELD (patch%litter_in) and the site state
   !      module cannot depend on this one without inverting the library DAG. Re-exported here so
   !      every existing `use meds_biogeochem_types, only : litter_input_t` keeps working.

   !----- Daily carbon-mass conservation guard (the fast/slow contract). --------------------------!
   type :: soilc_audit_t
      !----- `litter_in_matrix`, not `litter_in`: this is sum(u), the litter that enters through the  !
      !      MATRIX source term -- turnover plus continuous background mortality. Cull-termination    !
      !      and disturbance-kill necromass are added DIRECTLY onto patch%soil_carbon by the          !
      !      demography operators and never pass through `u`, so this is not the patch's litter in.   !
      !      The name says which of the two it is; the old one did not.  ---------------------------!
      real(wp) :: litter_in_matrix = 0.0_wp !< [kgC/m2/day] sum(u) -- the matrix source term ONLY
      real(wp) :: rh_out       = 0.0_wp   !< [kgC/m2/day] Rh reported by soil_carbon_step (= litter_in_matrix - dC_pool)
      real(wp) :: rh_fast_accum= 0.0_wp   !< [kgC/m2/day] fast loop's accumulated today_rh (the CAS-fed flux)
      real(wp) :: dC_pool      = 0.0_wp   !< [kgC/m2/day] net pool change
      !----- rh_out - rh_fast_accum: the fast/slow reconciliation check. Zero BY CONSTRUCTION on an  !
      !      ordinary day (measured 1.1e-14 kgC/m2 over a 31-day July), because both ends read the    !
      !      same frozen pool and the same per-pool environmental scalar.                              !
      !                                                                                          !
      !      IT IS NOT ZERO ON A DAY WHEN PATCH STRUCTURE CHANGES, and that is expected rather than    !
      !      a defect (#192). Measured 8.370e-4 kgC/m2 at a year rollover over a 50-year spin-up,      !
      !      when the annual patch cadence fires: the fast loop accumulates rh_fast_accum against one  !
      !      patch composition and the daily step debits a different one, with blend_xi_accum and      !
      !      blend_soil_carbon area-weighting the two ends separately through a matrix that is         !
      !      NONLINEAR in the lignin fraction. The contract is "equal by construction given the same   !
      !      patch composition", and a structural-change day does not supply that premise.             !
      !                                                                                          !
      !      Do NOT widen a tolerance to absorb it -- that would blind the check on every ordinary     !
      !      day, which is where it earns its keep. Read a nonzero gap on a non-structural day as a    !
      !      broken contract; on a structural-change day, read it as the composition blend.            !
      real(wp) :: rh_seam_gap  = 0.0_wp   !< [kgC/m2/day] rh_out - rh_fast_accum (see above)
      real(wp) :: lignin_resid = 0.0_wp   !< [kgC/m2/day] max_s |dL_s - (lignin_in_s - d_s*L_s)|; passive-tracer check
      !----- THE FREEZE NUMBER. dvec_j = xi_int_j * K_j is the FRACTION of pool j this slow step       !
      !      withdraws, and it is the number that says whether freezing the pool across the step is    !
      !      sound at all. It is dimensionless and scale-free: because the flux is linear in the pool  !
      !      (F = k.S), lambda does NOT grow as the pool empties, which is why bare-ground spin-up     !
      !      works. MEASURED 3.6e-3/day on the mature Ithaca stand (fast_grnd, July), i.e. a daily      !
      !      step withdraws ~0.4% of the fastest pool; approaching 1 means dt_slow is too long          !
      !      for that pool, and past 1 a forward-Euler step drives it negative.                        !
      !      Distinct from rh_seam_gap: the seam gap catches a BROKEN CONTRACT, lambda catches the     !
      !      APPROXIMATION DEGRADING, before anything breaks.  ---------------------------------------!
      real(wp)    :: lambda_max  = 0.0_wp   !< [-] max_j dvec_j = max_j (xi_int_j * K_j)
      integer(ik) :: lambda_pool = 0_ik     !< which pool attained it (IP_* index)
   end type soilc_audit_t

   !----- Per-RUN worsts for the three soil-carbon seam diagnostics, so the drivers thread ONE object !
   !      rather than an optional real per number (which is how worst_rh_seam_gap spent its life      !
   !      unreported: adding a second was going to mean a second argument through three layers).      !
   type :: soilc_seam_t
      real(wp)    :: worst_rh_gap  = 0.0_wp  !< [kgC/m2/day] max |rh_out - rh_fast_accum|
      real(wp)    :: worst_lignin  = 0.0_wp  !< [kgC/m2/day] max |lignin passive-tracer residual|
      real(wp)    :: worst_lambda  = 0.0_wp  !< [-] max fraction of any pool withdrawn in one step
      integer(ik) :: lambda_pool   = 0_ik    !< the pool that attained worst_lambda
   end type soilc_seam_t

   !----- Traceability diagnostics (pure post-processing off the assembled matrices). -------------!
   type :: soilc_diag_t
      real(wp) :: x_c(n_soil_pool) = 0.0_wp   !< [kgC/m2] storage CAPACITY (equilibrium under current climate+input)
      real(wp) :: x_p(n_soil_pool) = 0.0_wp   !< [kgC/m2] storage POTENTIAL = x_c - X (remaining sink; <0 = source)
      real(wp) :: tau_e            = 0.0_wp   !< [yr]     ecosystem soil-carbon residence time = sum(X_c)/sum(u)
      real(wp) :: rh               = 0.0_wp   !< [kgC/m2/day] instantaneous Rh = -1^T*A*xi*K*X
   end type soilc_diag_t

end module meds_biogeochem_types
