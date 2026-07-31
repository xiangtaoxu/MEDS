!==========================================================================================!
! meds_fast_step -- the fast-loop TIME-INTEGRATOR DISPATCH. One thin routine, `column_fast_step`, !
! that hands a dt_fast sub-step to the chosen scheme and owns the RK45 -> ARK stiff rescue.        !
!                                                                                          !
! WHAT USED TO BE HERE. This file replaces `meds_fast_split.f90`, which held BOTH the dispatcher   !
! and a third integrator -- the operator-split + Picard scheme that was the historical default.    !
! That scheme is RETIRED (2026-07-31). The reasons, in the order they matter:                      !
!                                                                                          !
!   1. It converges to a DIFFERENT limit than ARK/RK45 (~0.45 K on the canopy air), which was      !
!      established earlier and means every "ARK agrees with split" assertion was anchored on       !
!      something that is not the answer.                                                           !
!   2. It cannot carry the prognostic tissue heat store, and not as a fixable bug: with correctly  !
!      sized wood, cap_wood/cap_cas ~ 0.5 while tau_wood << dt_fast, so tissue and canopy air are  !
!      a COUPLED stiff pair that needs a joint implicit solve. A single explicit pass diverged     !
!      (canopy air alternating 291 K <-> 311 K per step, soil DRYING under a prescribed water      !
!      input); forcing the Picard iterate diverged too; Schur-preconditioning restored contraction !
!      but left the budgets open. ARK closes every budget on the same physics because              !
!      newton_surface_solve already couples the two. See MEDS_VEG_ENERGY_INTEGRATION_PLAN.md §8.   !
!   3. It made every physics addition a three-way wiring exercise, and it is where parity broke.   !
!                                                                                          !
! ED2 PROVENANCE, stated accurately: ED2's default is INTEGRATION_SCHEME = 1, the 4th-order        !
! Runge-Kutta (ED2IN:803, which explicitly recommends it). ED2 DOES have an operator split -- its  !
! INTEGRATION_SCHEME = 3 "hybrid" runs the canopy implicitly via bdf2_solver and the rest through  !
! the RK path -- but it is an alternative there, not the default. MEDS's RK45 covers the ED2       !
! default lineage; nothing in ED2 requires us to keep our own split.                               !
!==========================================================================================!
module meds_fast_step
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : cp_air, latent_heat_vap
   use meds_config,           only : meds_config_t, INTEG_ARK, INTEG_RK4
   use meds_biophysics_types, only : aero_env_t, aero_geom_t, aero_out_t, patch_biophys_t
   use meds_fast_types,       only : column_config_t, column_cohort_t, column_forcing_t,          &
                                     column_budget_t
   use meds_fast_ark,         only : column_fast_step_ark
   use meds_fast_rk45,        only : column_fast_step_rk45, rk45_state_railed
   implicit none
   private

   public :: column_fast_step

contains

   !---------------------------------------------------------------------------------------!
   ! Advance ONE dt_fast sub-step with the configured scheme.                                   !
   !                                                                                          !
   !   time_integrator = "ark"  (DEFAULT) -- the coupled implicit scheme. Despite the historical  !
   !                                        name it is NOT an IMEX method: the biotic CO2 source   !
   !                                        is folded implicit, so f_E == 0 and the tableau is a    !
   !                                        clean 2-solve ESDIRK2 with gamma = 1 - 1/sqrt(2)        !
   !                                        (the ARS(2,2,2) value). The config value stays "ark"    !
   !                                        for compatibility; the SCHEME is ESDIRK2.               !
   !   time_integrator = "rk45"           -- the fully explicit adaptive Cash-Karp march, kept as    !
   !                                        the ACCURACY BASELINE. Deliberately not optimised.       !
   !---------------------------------------------------------------------------------------!
   subroutine column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh, &
                               leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters,        &
                               le_flux, h_flux)
      real(wp),                intent(in)    :: dt_fast
      type(meds_config_t),     intent(in)    :: cfg          !< PFT traits for leaf gas exchange
      type(column_config_t),   intent(in)    :: ccfg
      type(aero_env_t),        intent(inout) :: aenv         !< can_* fields refreshed from CAS state
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: coh
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(inout) :: bio
      type(aero_out_t),        intent(inout) :: aero         !< preallocated (alloc_aero_out)
      type(column_budget_t),   intent(inout) :: budg
      real(wp), optional,      intent(out)   :: gpp_coh(:)   !< [umol CO2/plant/s] per-cohort GROSS GPP (fast->slow)
      real(wp), optional,      intent(out)   :: leaf_resp_coh(:) !< [umol CO2/plant/s] leaf dark respiration
      real(wp), optional,      intent(out)   :: stem_resp_coh(:) !< [umol CO2/plant/s] stem maintenance resp
      real(wp), optional,      intent(out)   :: root_resp_coh(:) !< [umol CO2/plant/s] fine-root maint. resp
      logical,     optional,   intent(out)   :: converged    !< scheme's inner solve converged this sub-step
      integer(ik), optional,   intent(out)   :: iters        !< outer-iteration count taken
      real(wp),    optional,   intent(out)   :: le_flux      !< [W/m2] CAS->atm latent-heat (ET) flux
      real(wp),    optional,   intent(out)   :: h_flux       !< [W/m2] CAS->atm sensible-heat flux

      budg%rk45_rescue = 0_ik

      if (cfg%time_integrator == INTEG_RK4) then
         !----- RK45 is FULLY EXPLICIT over the whole column (no implicit canopy-air box). At high LAI  !
         !      plus cold, the coupled leaf<->CAS exchange is stiff enough that the explicit stages      !
         !      cannot resolve it even at the sub-step floor: the state rails to the clamp bounds        !
         !      (CAS/soil pinned at 180 or 350 K) and -- worse -- the 5th- and 4th-order embedded        !
         !      solutions rail TOGETHER, so the controller sees err~0 and commits the railed garbage,    !
         !      locking into a cold-dead attractor (GPP -> 0) that never recovers.                       !
         !                                                                                          !
         !      HYBRID RESCUE: snapshot state^n, take the explicit step, and roll back + REDO this       !
         !      dt_fast on ARK when either trigger fires:                                                !
         !        * stiff_bail        -- the march burned its whole work budget without resolving.        !
         !        * rk45_state_railed -- it finished, but committed a clamp-pinned CAS/soil state.        !
         !                                                                                          !
         !      The rescue target USED TO BE the operator-split path. It is ARK now, and ARK is the      !
         !      better target on its own merits: the thing being rescued from is a COUPLED STIFF         !
         !      canopy-air problem, and ARK's newton_surface_solve solves exactly that pair implicitly,  !
         !      which is also why ARK carries the tissue heat store while split could not. RK45 keeps    !
         !      its fast explicit path everywhere it is stable; only genuinely stiff steps pay. ---------!
         block
            type(patch_biophys_t) :: bio_save
            type(column_budget_t) :: budg_save
            logical               :: rk45_stiff
            bio_save = bio ; budg_save = budg
            call column_fast_step_rk45(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, &
                                       gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh,        &
                                       converged, iters, stiff_bail=rk45_stiff)
            if (.not. rk45_stiff .and. .not. rk45_state_railed(bio, ccfg%soil%n_active)) then
               call atm_fluxes(aenv, aero, bio, forc, le_flux, h_flux)
               return
            end if
            bio = bio_save ; budg = budg_save         ! discard the railed/bailed RK45 step (rollback)
            budg%integ_nrej  = budg%integ_nrej  + 1_ik   ! count the rescue as a rejected integrator step
            budg%rk45_rescue = budg%rk45_rescue + 1_ik   ! ...and as an RK45->ARK rescue (diagnostic)
         end block
      end if

      !----- ARK (ESDIRK2): the default, and the RK45 rescue target. ----------------------------!
      call column_fast_step_ark(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg,      &
                                gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters)
      call atm_fluxes(aenv, aero, bio, forc, le_flux, h_flux)
   end subroutine column_fast_step

   !----- CAS->atmosphere turbulent fluxes, reported identically whichever scheme ran. --------!
   pure subroutine atm_fluxes(aenv, aero, bio, forc, le_flux, h_flux)
      type(aero_env_t),       intent(in)  :: aenv
      type(aero_out_t),       intent(in)  :: aero
      type(patch_biophys_t),  intent(in)  :: bio
      type(column_forcing_t), intent(in)  :: forc
      real(wp), optional,     intent(out) :: le_flux, h_flux
      if (present(le_flux)) le_flux = aenv%rho_air * aero%ustar * aero%temp2                      &
                                      * (bio%cas%can_shv - forc%shv_atm) * latent_heat_vap
      if (present(h_flux))  h_flux  = aenv%rho_air * aero%ustar * aero%temp2                      &
                                      * (bio%cas%can_temp - aenv%theta_atm) * cp_air
   end subroutine atm_fluxes

end module meds_fast_step
