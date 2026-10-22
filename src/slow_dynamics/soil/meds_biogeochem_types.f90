! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_biogeochem_types -- the audit and diagnostic records of the SLOW soil-carbon cycle: the    !
! CENTURY-family multi-pool soil carbon advanced daily by the carbon matrix ODE                    !
! dX/dt = B*I + A*xi*K*X (docs/science/soil_carbon.md eq. (1)). Pure DATA, no methods -- the       !
! biogeochem analogue of meds_canopy_types/meds_soil_types. Links meds_kinds,                     !
! meds_patch_state_vars and meds_biogeochem_params only.                                           !
!                                                                                          !
! Owned HERE: soilc_audit_t, soilc_seam_t and soilc_diag_t. The pools themselves (soil_carbon_t)  !
! and litter_input_t are per-patch state in meds_patch_state_vars (litter_input_t is re-exported   !
! here; see the notes below), and the pool layout (N_SOIL_POOL, IP_*) and soil_carbon_params_t +   !
! its selector codes are config, in meds_biogeochem_params.                                        !
!==========================================================================================!
module meds_biogeochem_types
   use meds_kinds, only : wp, ik
   use meds_patch_state_vars,   only : litter_input_t
   use meds_biogeochem_params,  only : N_SOIL_POOL
   implicit none
   private

   !----- Slow soil-carbon matrix: its audit and diagnostic records. ----------------------------!
   public :: litter_input_t, soilc_audit_t, soilc_seam_t, soilc_diag_t
   !----- There is exactly ONE heterotrophic-respiration form: the CENTURY matrix                    !
   !      (heterotrophic_respiration_matrix), reached through patch_heterotrophic_respiration.       !
   !      No Q10, capped-exponential or DAMM kernel is available (docs/science/soil_carbon.md        !
   !      §7, §10). ---------------------------------------------------------------------------------!

   !==========================================================================================!
   !  The slow, stateful per-patch soil-carbon pools themselves (soil_carbon_t: written DAILY,       !
   !  read-only in the fast loop) are defined in meds_patch_state_vars, with the patch's other       !
   !  stores: they are what a patch holds (src/README.md rule 3), and the kernels take them from     !
   !  there. The ACTIVE pool count is set by decomp_scheme (default 3-active: idx 5, 7 inert).       !
   !==========================================================================================!

   !----- litter_input_t is defined in meds_patch_state_vars, beside the soil_carbon_t pools it
   !      feeds: it is a PER-PATCH FIELD (patch%litter_in), and the state layer cannot depend on
   !      this kernel module without inverting the library graph. Re-exported here.

   !----- Daily carbon-mass conservation guard (the fast/slow contract). --------------------------!
   type :: soilc_audit_t
      !----- `litter_in_matrix`, not `litter_in`: this is sum(u), the litter that enters through the  !
      !      MATRIX source term -- turnover plus continuous background mortality. Cull-termination    !
      !      and disturbance-kill necromass are added DIRECTLY onto patch%soil_carbon by the          !
      !      demography operators and never pass through `u`, so this is not the patch's litter in. -!
      real(wp) :: litter_in_matrix = 0.0_wp !< [kgC/m2/day] sum(u) -- the matrix source term ONLY
      real(wp) :: rh_out       = 0.0_wp   !< [kgC/m2/day] Rh reported by soil_carbon_step (= litter_in_matrix - dC_pool)
      real(wp) :: rh_fast_accum= 0.0_wp   !< [kgC/m2/day] fast loop's accumulated matrix Rh (the CAS-fed flux)
      real(wp) :: dC_pool      = 0.0_wp   !< [kgC/m2/day] net pool change
      !----- rh_out - rh_fast_accum: the fast/slow reconciliation check (docs/science/soil_carbon.md   !
      !      §7). Zero BY CONSTRUCTION on an ordinary day, because both ends read the same frozen      !
      !      pool and the same per-pool environmental scalar. NOT zero on a day the patch structure    !
      !      changes, and that is expected: the fast loop respired one patch composition and the       !
      !      daily step debits the area-blended one. Do NOT widen a tolerance to absorb it -- that     !
      !      would blind the check on every ordinary day, which is where it earns its keep.            !
      real(wp) :: rh_seam_gap  = 0.0_wp   !< [kgC/m2/day] rh_out - rh_fast_accum (see above)
      real(wp) :: lignin_resid = 0.0_wp   !< [kgC/m2/day] max_s |dL_s - (lignin_in_s - (1-surv_s)*L_s)|; tracer check
      !----- THE FREEZE NUMBER (docs/science/soil_carbon.md §6). dvec_j = xi_int_j * K_j is the        !
      !      FRACTION of pool j this slow step withdraws, the number that says whether freezing the    !
      !      pool across the step is sound. Scale-free (the flux is linear in the pool), so it does    !
      !      not grow as a pool empties; approaching 1 means dt_slow is too long for that pool, and    !
      !      past 1 a forward-Euler step drives it negative. Distinct from rh_seam_gap: the seam gap   !
      !      catches a BROKEN CONTRACT, lambda the APPROXIMATION DEGRADING, before anything breaks. ---!
      real(wp)    :: lambda_max  = 0.0_wp   !< [-] max_j dvec_j = max_j (xi_int_j * K_j)
      integer(ik) :: lambda_pool = 0_ik     !< which pool attained it (IP_* index)
   end type soilc_audit_t

   !----- Per-RUN worsts for the three soil-carbon seam diagnostics, so the drivers thread ONE object !
   !      through three layers rather than an optional real per number. ------------------------------!
   type :: soilc_seam_t
      real(wp)    :: worst_rh_gap  = 0.0_wp  !< [kgC/m2/day] max |rh_out - rh_fast_accum|
      real(wp)    :: worst_lignin  = 0.0_wp  !< [kgC/m2/day] max |lignin passive-tracer residual|
      real(wp)    :: worst_lambda  = 0.0_wp  !< [-] max fraction of any pool withdrawn in one step
      integer(ik) :: lambda_pool   = 0_ik    !< the pool that attained worst_lambda
   end type soilc_seam_t

   !----- Traceability diagnostics (pure post-processing off the assembled matrices). -------------!
   type :: soilc_diag_t
      real(wp) :: x_c(N_SOIL_POOL) = 0.0_wp   !< [kgC/m2] storage CAPACITY (equilibrium under current climate+input)
      real(wp) :: x_p(N_SOIL_POOL) = 0.0_wp   !< [kgC/m2] storage POTENTIAL = x_c - X (remaining sink; <0 = source)
      real(wp) :: tau_e            = 0.0_wp   !< [yr]     ecosystem soil-carbon residence time = sum(X_c)/sum(u)
      real(wp) :: rh               = 0.0_wp   !< [kgC/m2/day] instantaneous Rh = -1^T*A*xi*K*X
   end type soilc_diag_t

end module meds_biogeochem_types
