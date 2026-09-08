!==========================================================================================!
! meds_fast_step -- the fast-loop TIME-INTEGRATOR DISPATCH. One thin routine, `column_fast_step`, !
! hands a dt_fast sub-step to the chosen scheme (INTEG_ARK, the default; INTEG_RK4 = the adaptive  !
! Cash-Karp RK45 in meds_fast_rk45) and owns the RK45 -> ARK stiff rescue: when the explicit march !
! bails or commits a clamp-railed CAS/soil state, the step is rolled back and redone on ARK.       !
! It also reports the per-cohort psi_leaf for the daily-max accumulator (both schemes) and the      !
! CAS -> atmosphere LE/H diagnostics. The operator-split integrator that used to live here was      !
! retired 2026-07-31 (docs/science/numerical_scheme.md records why and the ED2 provenance).         !
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
   use meds_hydr_lib,        only : psi_from_water_content
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
                               le_flux, h_flux, psi_leaf_coh, cdiag)
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
      real(wp), optional,      intent(out)   :: psi_leaf_coh(:) !< [MPa] this step's per-cohort psi_leaf (daily-max accumulator)
      real(wp), optional,      intent(out)   :: leaf_resp_coh(:) !< [umol CO2/plant/s] leaf dark respiration
      real(wp), optional,      intent(out)   :: stem_resp_coh(:) !< [umol CO2/plant/s] stem maintenance resp
      real(wp), optional,      intent(out)   :: root_resp_coh(:) !< [umol CO2/plant/s] fine-root maint. resp
      logical,     optional,   intent(out)   :: converged    !< scheme's inner solve converged this sub-step
      integer(ik), optional,   intent(out)   :: iters        !< outer-iteration count taken
      real(wp),    optional,   intent(out)   :: le_flux      !< [W/m2] CAS->atm latent-heat (ET) flux
      real(wp),    optional,   intent(out)   :: h_flux       !< [W/m2] CAS->atm sensible-heat flux
      !----- OPTIONAL per-cohort diagnostic capture (MEDS_IO_V01_PLAN.md section 3.4). Filled by the  !
      !      Act-1 pre-pass with the leaf gas-exchange + hydraulics quantities that were previously    !
      !      recomputed every dt_fast and discarded. Absent => nothing extra is computed at all.  ----!
      real(wp),    optional,   intent(inout) :: cdiag(:,:)
      integer(ik) :: jcoh

      budg%rk45_rescue = 0_ik

      !----- Report this step's per-cohort psi_leaf for the daily-MAX accumulator that drives      !
      !      beta_stomata on the NEXT day (issue #95). Diagnosed from leaf_water_mass at state^n,     !
      !      which is exactly the value the leaf kernel was handed this step (column_prepass freezes  !
      !      psi_leaf once per dt_fast), so the accumulator and the kernel can never disagree.        !
      !      It runs BEFORE the scheme dispatch so that BOTH schemes fill it: it used to sit after    !
      !      the RK45 block, whose success path returns early, leaving the caller's per-thread        !
      !      buffer with the previous patch's (or uninitialised) values under time_integrator=rk45. -!
      if (present(psi_leaf_coh)) then
         do jcoh = 1_ik, coh%n
            psi_leaf_coh(jcoh) = psi_from_water_content(bio%leaf_water_mass(jcoh),                  &
                 ccfg%hydro_p%leaf_pi0, ccfg%hydro_p%leaf_elastic_mod,                              &
                 ccfg%hydro_p%leaf_apoplast_frac, ccfg%hydro_p%leaf_water_sat, coh%bleaf(jcoh))
         end do
      end if

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
                                       converged, iters, stiff_bail=rk45_stiff, cdiag=cdiag)
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
                                gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged,   &
                                iters, cdiag)
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
