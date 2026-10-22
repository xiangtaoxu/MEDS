! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_step -- the fast-loop TIME-INTEGRATOR DISPATCH. One thin routine, `column_fast_step`, !
! hands a dt_fast sub-step to the chosen scheme (INTEG_ARK, the default; INTEG_RK45 = the adaptive  !
! Cash-Karp RK45 in meds_fast_rk45) and owns the RK45 -> ARK stiff rescue: when the explicit march !
! bails or commits a clamp-railed CAS/soil state, the step is rolled back and redone on ARK.       !
! It also reports the per-cohort psi_leaf for the daily-max accumulator (both schemes) and the      !
! CAS -> atmosphere LE/H diagnostics. There are TWO schemes and no third: see                       !
! docs/science/numerical_scheme.md §3 for why there is no operator-split path.                      !
!==========================================================================================!
module meds_fast_step
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : latent_heat_vap
   use meds_config,           only : meds_config_t
   use meds_fast_params,      only : INTEG_RK45
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t
   use meds_fast_types, only : patch_biophys_t, column_config_t, column_cohort_t, column_forcing_t,               &
                               column_budget_t, fast_workspace_t
   use meds_fast_ark,         only : column_fast_step_ark
   use meds_fast_rk45,        only : column_fast_step_rk45, rk45_state_railed, rk45_ledgers_stop
   use meds_plant_hydraulics, only : leaf_node_psi
   implicit none
   private

   public :: column_fast_step

contains

   !---------------------------------------------------------------------------------------!
   ! Advance ONE dt_fast sub-step with the configured scheme.                                       !
   !                                                                                                !
   !   time_integrator = "ark"  (DEFAULT) -- the coupled implicit scheme. Despite the `ark` label   !
   !                                        it is NOT an IMEX method: the biotic CO2 source is      !
   !                                        folded implicit, so f_E == 0 and the tableau is a       !
   !                                        clean 2-solve ESDIRK2 with gamma = 1 - 1/sqrt(2)        !
   !                                        (the ARS(2,2,2) value; docs/science/numerical_scheme.md !
   !                                        eq. (1)). The config value stays "ark" for              !
   !                                        compatibility; the SCHEME is ESDIRK2.                   !
   !   time_integrator = "rk45"           -- the fully explicit adaptive Cash-Karp march, kept as   !
   !                                        the ACCURACY BASELINE. Deliberately not optimised.      !
   !---------------------------------------------------------------------------------------!
   subroutine column_fast_step(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, &
                               gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters, &
                               le_flux, h_flux, psi_leaf_coh, cdiag)
      real(wp),                intent(in)    :: dt_fast
      type(meds_config_t),     intent(in)    :: cfg          !< PFT traits for leaf gas exchange
      type(column_config_t),   intent(in)    :: col_config
      type(aero_env_t),        intent(inout) :: aenv         !< can_* fields refreshed from CAS state
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: col_cohort
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(inout) :: biophys
      type(aero_out_t),        intent(inout) :: aero         !< preallocated (alloc_aero_out)
      type(column_budget_t),   intent(inout) :: budget
      type(fast_workspace_t),  intent(inout) :: ws           !< this thread's step storage, reused step to step
      real(wp), optional,      intent(out)   :: gpp_coh(:)   !< [umol CO2/plant/s] per-cohort GROSS GPP (fast->slow)
      real(wp), optional,      intent(out)   :: psi_leaf_coh(:) !< [MPa] this step's per-cohort psi_leaf (daily-max accumulator)
      real(wp), optional,      intent(out)   :: leaf_resp_coh(:) !< [umol CO2/plant/s] leaf dark respiration
      real(wp), optional,      intent(out)   :: stem_resp_coh(:) !< [umol CO2/plant/s] stem maintenance resp
      real(wp), optional,      intent(out)   :: root_resp_coh(:) !< [umol CO2/plant/s] fine-root maint. resp
      logical,     optional,   intent(out)   :: converged    !< scheme's inner solve converged this sub-step
      integer(ik), optional,   intent(out)   :: iters        !< outer-iteration count taken
      real(wp),    optional,   intent(out)   :: le_flux      !< [W/m2] CAS->atm latent-heat (ET) flux
      real(wp),    optional,   intent(out)   :: h_flux       !< [W/m2] CAS->atm sensible-heat flux
      !----- OPTIONAL per-cohort diagnostic capture (docs/science/diagnostics.md §1, stage [2]).   !
      !      Filled by the pre-pass with the leaf gas-exchange + hydraulics quantities it computes !
      !      every dt_fast anyway. Absent => nothing extra is computed at all. --------------------!
      real(wp),    optional,   intent(inout) :: cdiag(:,:)
      integer(ik) :: jcoh

      budget%rk45_rescue = 0_ik

      !----- Report this step's per-cohort psi_leaf for the daily-MAX accumulator that drives        !
      !      beta_stomata on the NEXT day (docs/science/leaf_gas_exchange.md §4.2). Diagnosed from   !
      !      leaf_water_mass at state^n, which is exactly the value the leaf kernel was handed this  !
      !      step (column_prepass freezes psi_leaf once per dt_fast), so the accumulator and the     !
      !      kernel can never disagree. It runs BEFORE the scheme dispatch so that BOTH schemes fill !
      !      it: the RK45 success path returns early, which would leave the caller's per-thread      !
      !      buffer with the previous patch's (or uninitialised) values. ----------------------------!
      if (present(psi_leaf_coh)) then
         do jcoh = 1_ik, col_cohort%n
            psi_leaf_coh(jcoh) = leaf_node_psi(biophys%leaf_water_mass(jcoh), col_cohort%bleaf(jcoh),       &
                  biophys%wood_water_mass(jcoh), col_cohort%bsap(jcoh) + col_cohort%broot(jcoh),          &
                  col_cohort%height(jcoh), col_config%hydraulics_table%pft(col_cohort%pft(jcoh)),         &
                  col_config%hydraulics_opts)
         end do
      end if

      if (col_config%integrator%scheme == INTEG_RK45) then
         !----- RK45 is FULLY EXPLICIT over the whole column (no implicit canopy-air box). At high LAI  !
         !      plus cold, the coupled leaf<->CAS exchange is stiff enough that the explicit stages     !
         !      cannot resolve it even at the sub-step floor: the state rails to the clamp bounds       !
         !      (CAS/soil pinned at 180 or 350 K) and -- worse -- the 5th- and 4th-order embedded       !
         !      solutions rail TOGETHER, so the controller sees err~0 and commits the railed garbage,   !
         !      locking into a cold-dead attractor (GPP -> 0) that never recovers.                      !
         !                                                                                              !
         !      HYBRID RESCUE (docs/science/numerical_scheme.md §3, `rk45`): snapshot state^n, take the !
         !      explicit step, and roll back + REDO this dt_fast on ARK when either trigger fires:      !
         !        * stiff_bail        -- the march burned its whole work budget without resolving.      !
         !        * rk45_state_railed -- it finished, but committed a clamp-pinned CAS/soil state.      !
         !                                                                                              !
         !      ARK is the rescue target because the thing being rescued from is a COUPLED STIFF        !
         !      canopy-air problem, and ARK's newton_surface_solve solves exactly that pair implicitly. !
         !      RK45 keeps its fast explicit path everywhere it is stable; only genuinely stiff steps   !
         !      pay. -----------------------------------------------------------------------------------!
         block
            type(patch_biophys_t) :: bio_save
            type(column_budget_t) :: budg_save
            logical               :: rk45_stiff
            bio_save = biophys
            budg_save = budget
            call column_fast_step_rk45(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, &
                                       ws, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh,        &
                                       converged, iters, stiff_bail=rk45_stiff, cdiag=cdiag)
            if (.not. rk45_stiff .and. .not. rk45_state_railed(biophys, col_config%soil%n_active)) then
               call rk45_ledgers_stop(budget, dt_fast, col_config)   ! kept, so a breach may now stop the run
               call atm_fluxes(budget, dt_fast, le_flux, h_flux)
               return
            end if
            !----- Whole-dt_fast rollback, not a resume from the point of failure. A resume would        !
            !      save redoing the sub-steps RK45 had already accepted, but there is no work to save on !
            !      the reference workload: over a full simulated year at Ithaca on rk45 at dt_fast =     !
            !      900 s the rescue fires ZERO times (mean_work_rk45_rescue_py = 0, against 70 577          !
            !      integrator sub-steps on the same run, so the counter is live and the zero is real).   !
            !                                                                                            !
            !      A resume is also not free. Retaining part of RK45's boundary-flux accumulation while  !
            !      ARK finishes the interval means one dt_fast's ledger summing two schemes'             !
            !      contributions -- the "borrow one solve's numbers, commit another's trajectory" class  !
            !      (docs/science/numerical_scheme.md §5b). And it is only well-defined for the           !
            !      stiff_bail trigger: rk45_state_railed tests the FINAL state, so on that path there    !
            !      is no identified sub-step to resume from. It would pay only on a forcing where        !
            !      rescues are common. ------------------------------------------------------------------!
            biophys = bio_save
            budget = budg_save         ! discard the railed/bailed RK45 step (rollback)
            budget%integ_nrej  = budget%integ_nrej  + 1_ik   ! count the rescue as a rejected integrator step
            budget%rk45_rescue = budget%rk45_rescue + 1_ik   ! ...and as an RK45->ARK rescue (diagnostic)
         end block
      end if

      !----- ARK (ESDIRK2): the default, and the RK45 rescue target. ----------------------------!
      call column_fast_step_ark(dt_fast, cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws,  &
                                gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged,   &
                                iters, cdiag)
      call atm_fluxes(budget, dt_fast, le_flux, h_flux)
   end subroutine column_fast_step

   !---------------------------------------------------------------------------------------!
   ! atm_fluxes -- the CAS -> atmosphere turbulent fluxes the run REPORTS, from the conservation        !
   ! ledger's own boundary export over the step: the flux the integrator charged the canopy air         !
   ! with, b-weighted over the stages at the stage states and (on a live-MO stage) the stage            !
   ! conductances, so the reported ET IS the conserved vapour export. Recomputing it from the           !
   ! END-of-step CAS state with the pre-pass conductance (rho*ustar*temp2) would give a different       !
   ! number.                                                                                            !
   !                                                                                                    !
   !   LE = L_v * (net vapour export) / dt   -- the vapour the CAS actually shed, times L_v;            !
   !   H  = (net sensible export) / dt        -- g_atm_heat*cp_air*(T_cas - theta_atm), b-weighted      !
   !                                             the same way (stage_bflux_t%atm_heat_out).             !
   !                                                                                                    !
   ! H is NOT "enthalpy export minus LE": the CAS enthalpy values vapour at cp_vap*(T - tsupercool_vap) !
   ! ~ 3.4 MJ/kg, so that difference would carry the water's liquid-datum enthalpy (~40% of LE) inside  !
   ! the sensible flux.                                                                                 !
   !---------------------------------------------------------------------------------------!
   pure subroutine atm_fluxes(budget, dt_fast, le_flux, h_flux)
      type(column_budget_t), intent(in)  :: budget
      real(wp),              intent(in)  :: dt_fast
      real(wp), optional,    intent(out) :: le_flux, h_flux   !< [W/m2] step-mean CAS -> atmosphere
      if (present(le_flux)) le_flux = budget%atm_vap_export  / dt_fast * latent_heat_vap
      if (present(h_flux))  h_flux  = budget%atm_heat_export / dt_fast
   end subroutine atm_fluxes

end module meds_fast_step
