! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_biogeochem_types -- the audit and diagnostic records of the daily soil-carbon step  !
! (docs/science/soil_carbon.md §6-§9): soilc_audit_t (one step's checks), soilc_seam_t     !
! (the run's worst values of them) and soilc_diag_t (the traceability diagnostics). Data   !
! only. The pools (soil_carbon_t) and the litter input (litter_input_t, re-exported here)  !
! are per-patch state in meds_patch_state_vars; the pool layout (N_SOIL_POOL, IP_*) and    !
! the parameters are in meds_biogeochem_params.                                            !
!==========================================================================================!
module meds_biogeochem_types
   use meds_kinds, only : wp, ik
   use meds_patch_state_vars,   only : litter_input_t
   use meds_biogeochem_params,  only : N_SOIL_POOL
   implicit none
   private

   public :: litter_input_t, soilc_audit_t, soilc_seam_t, soilc_diag_t

   !----- One daily step's carbon-mass checks (soil_carbon.md §7, §9). --------------------!
   type :: soilc_audit_t
      !----- litter_in_matrix is sum(u) only: cull and disturbance kills add their         !
      !      necromass straight onto patch%soil_carbon, outside u (§5). -------------------!
      real(wp) :: litter_in_matrix = 0.0_wp !< [kgC/m2/day] sum(u) -- the matrix source term only
      real(wp) :: rh_out       = 0.0_wp   !< [kgC/m2/day] Rh reported by soil_carbon_step (= litter_in_matrix - dC_pool)
      real(wp) :: rh_fast_accum= 0.0_wp   !< [kgC/m2/day] fast loop's accumulated matrix Rh (the flux fed to the canopy air)
      real(wp) :: dC_pool      = 0.0_wp   !< [kgC/m2/day] net pool change
      !----- rh_seam_gap: zero by construction on an ordinary day, nonzero as expected on a !
      !      day the patch structure changes; do not widen a tolerance to absorb it (§7). ---!
      real(wp) :: rh_seam_gap  = 0.0_wp   !< [kgC/m2/day] rh_out - rh_fast_accum (see above)
      real(wp) :: lignin_resid = 0.0_wp   !< [kgC/m2/day] max_s |dL_s - (lignin_in_s - (1-surv_s)*L_s)|; tracer check
      !----- lambda, the largest fraction of any pool the step withdraws, says whether     !
      !      freezing the pool is sound; past 1 forward Euler goes negative (§6). ---------!
      real(wp)    :: lambda_max  = 0.0_wp   !< [-] max_j dvec_j = max_j (xi_int_j * K_j)
      integer(ik) :: lambda_pool = 0_ik     !< which pool attained it (IP_* index)
   end type soilc_audit_t

   !----- The run's worst values of the three soil-carbon checks, carried as one object. ---!
   type :: soilc_seam_t
      real(wp)    :: worst_rh_gap  = 0.0_wp  !< [kgC/m2/day] max |rh_out - rh_fast_accum|
      real(wp)    :: worst_lignin  = 0.0_wp  !< [kgC/m2/day] max |lignin passive-tracer residual|
      real(wp)    :: worst_lambda  = 0.0_wp  !< [-] max fraction of any pool withdrawn in one step
      integer(ik) :: lambda_pool   = 0_ik    !< the pool that attained worst_lambda
   end type soilc_seam_t

   !----- Traceability diagnostics (pure post-processing of the assembled matrices; §8). ---!
   type :: soilc_diag_t
      real(wp) :: x_c(N_SOIL_POOL) = 0.0_wp   !< [kgC/m2] storage capacity (equilibrium under current climate+input)
      real(wp) :: x_p(N_SOIL_POOL) = 0.0_wp   !< [kgC/m2] storage potential = x_c - X (remaining sink; <0 = source)
      real(wp) :: tau_e            = 0.0_wp   !< [yr]     ecosystem soil-carbon residence time = sum(X_c)/sum(u)
      real(wp) :: rh               = 0.0_wp   !< [kgC/m2/day] instantaneous Rh = -1^T*A*xi*K*X
   end type soilc_diag_t

end module meds_biogeochem_types
