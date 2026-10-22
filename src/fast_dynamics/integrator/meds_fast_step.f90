! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_step -- the fast-loop integrator dispatch. column_fast_step hands a dt_fast step !
! to the scheme the run chose: ARK (meds_fast_ark, the default) or RK45 (meds_fast_rk45,   !
! the explicit Cash-Karp integrator, also the explicit reference the tests compare ARK     !
! against). It owns the RK45 -> ARK rescue, reports each cohort's psi_leaf for the daily-  !
! maximum accumulator, and reports H and LE (atm_fluxes). There is no operator-split       !
! scheme (docs/science/numerical_scheme.md §3).                                            !
!==========================================================================================!
module meds_fast_step
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : latent_heat_vap
   use meds_fast_params,      only : INTEG_RK45
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t
   use meds_fast_types, only : patch_biophys_t, column_config_t, column_cohort_t, column_forcing_t,               &
                               column_budget_t, fast_workspace_t, column_site_t
   use meds_fast_ark,         only : column_fast_step_ark
   use meds_fast_rk45,        only : column_fast_step_rk45, rk45_state_railed, rk45_ledgers_stop
   use meds_plant_hydraulics, only : leaf_node_psi
   implicit none
   private

   public :: column_fast_step

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one dt_fast with the configured scheme ([fast].time_integrator):                 !
   !   "ark" (default) -- the 2-solve ESDIRK2, gamma = 1 - 1/sqrt(2); not IMEX, since the     !
   !                      biotic CO2 source is implicit and the explicit tableau is empty     !
   !                      (docs/science/numerical_scheme.md eq. (1)).                         !
   !   "rk45"          -- the adaptive explicit Cash-Karp march, a full alternative a run can  !
   !                      choose and the explicit reference the tests compare ARK against.    !
   !---------------------------------------------------------------------------------------!
   subroutine column_fast_step(dt_fast, col_config, col_site, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws, &
                               gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters, &
                               le_flux, h_flux, psi_leaf_coh, cdiag)
      real(wp),                intent(in)    :: dt_fast
      type(column_config_t),   intent(in)    :: col_config
      type(column_site_t),     intent(in)    :: col_site
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
      !----- Optional per-cohort diagnostics (docs/science/diagnostics.md §1, stage [2]), filled by  !
      !      the pre-pass from what it computes anyway; absent, nothing extra is computed. -------!
      real(wp),    optional,   intent(inout) :: cdiag(:,:)
      integer(ik) :: jcoh

      budget%rk45_rescue = 0_ik

      !----- Each cohort's psi_leaf for the daily-maximum accumulator behind the next day's stomatal  !
      !      control (docs/science/leaf_gas_exchange.md §4.2), at state^n, the value the leaf kernel   !
      !      was handed. Before the dispatch, so both schemes fill it (RK45's success path returns). -!
      if (present(psi_leaf_coh)) then
         do jcoh = 1_ik, col_cohort%n
            psi_leaf_coh(jcoh) = leaf_node_psi(biophys%leaf_water_mass(jcoh), col_cohort%bleaf(jcoh),       &
                  biophys%wood_water_mass(jcoh), col_cohort%bsap(jcoh) + col_cohort%broot(jcoh),          &
                  col_cohort%height(jcoh), col_config%hydraulics_table%pft(col_cohort%pft(jcoh)))
         end do
      end if

      if (col_config%integrator%scheme == INTEG_RK45) then
         !----- RK45 is explicit over the whole column. On a dense, cold canopy the leaf <-> canopy-air !
         !      coupling can be too stiff even at the sub-step floor: the state rails, and its 5th- and  !
         !      4th-order solutions rail together, so the controller would accept it. The rescue        !
         !      (docs/science/numerical_scheme.md §3, `rk45`): keep state^n, take the explicit step,    !
         !      and redo the dt_fast on ARK, whose Newton solves that coupled pair, when the march      !
         !      gave up (stiff_bail) or committed a railed state (rk45_state_railed). -----------------!
         block
            type(patch_biophys_t) :: bio_save
            type(column_budget_t) :: budg_save
            logical               :: rk45_stiff
            bio_save = biophys
            budg_save = budget
            call column_fast_step_rk45(dt_fast, col_config, col_site, aenv, ageom, col_cohort, forc, biophys, aero, budget, &
                                       ws, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh,        &
                                       converged, iters, stiff_bail=rk45_stiff, cdiag=cdiag)
            if (.not. rk45_stiff .and. .not. rk45_state_railed(biophys, col_site%soil_profile%n_active)) then
               call rk45_ledgers_stop(budget, dt_fast, col_config, col_site)   ! kept, so a breach may now stop the run
               call atm_fluxes(budget, dt_fast, le_flux, h_flux)
               return
            end if
            !----- Roll back the whole dt_fast rather than resume from the failure: rescues are rare, a  !
            !      resume would sum two schemes' boundary fluxes in one step's ledger                   !
            !      (docs/science/numerical_scheme.md §5b), and a railed final state has no sub-step to  !
            !      resume from. -------------------------------------------------------------------------!
            biophys = bio_save
            budget = budg_save         ! discard the railed/bailed RK45 step (rollback)
            budget%integ_nrej  = budget%integ_nrej  + 1_ik   ! count the rescue as a rejected integrator step
            budget%rk45_rescue = budget%rk45_rescue + 1_ik   ! ...and as an RK45->ARK rescue (diagnostic)
         end block
      end if

      !----- ARK (ESDIRK2): the default, and the RK45 rescue target. ----------------------------!
      call column_fast_step_ark(dt_fast, col_config, col_site, aenv, ageom, col_cohort, forc, biophys, aero, budget, ws,  &
                                gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged,   &
                                iters, cdiag)
      call atm_fluxes(budget, dt_fast, le_flux, h_flux)
   end subroutine column_fast_step

   !---------------------------------------------------------------------------------------!
   ! atm_fluxes -- the canopy-air -> atmosphere fluxes the run reports, from the ledger's own  !
   ! boundary export over the step, so the reported ET is the conserved vapour export:         !
   ! LE = L_v * vapour export / dt and H = sensible export / dt, the b-weighted              !
   ! g_atm_heat*cp_air*(T_cas - theta_atm), not the enthalpy export minus LE                   !
   ! (docs/science/numerical_scheme.md §5b).                                                    !
   !---------------------------------------------------------------------------------------!
   pure subroutine atm_fluxes(budget, dt_fast, le_flux, h_flux)
      type(column_budget_t), intent(in)  :: budget
      real(wp),              intent(in)  :: dt_fast
      real(wp), optional,    intent(out) :: le_flux, h_flux   !< [W/m2] step-mean CAS -> atmosphere
      if (present(le_flux)) le_flux = budget%atm_vap_export  / dt_fast * latent_heat_vap
      if (present(h_flux))  h_flux  = budget%atm_heat_export / dt_fast
   end subroutine atm_fluxes

end module meds_fast_step
