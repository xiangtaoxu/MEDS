! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_step -- one dt_fast of one patch. fast_step takes the run's configuration, the !
! patch's site, its static state (what is fixed for the slow step), the forcing, the       !
! dynamic state it advances, the output it fills and this thread's workspace. It reports   !
! each cohort's psi_leaf for the daily-maximum accumulator, hands the step to the scheme   !
! the run chose (ARK, meds_fast_ark, the default; or RK45, meds_fast_rk45, the explicit    !
! Cash-Karp integrator and the tests' explicit reference), owns the RK45 -> ARK rescue,    !
! and reports H and LE (atm_fluxes). There is no operator-split scheme                     !
! (docs/science/numerical_scheme.md §3).                                                   !
!==========================================================================================!
module meds_fast_step
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : latent_heat_vap
   use meds_config,           only : meds_config_t
   use meds_fast_params,      only : INTEG_RK45
   use meds_site_state_vars,  only : site_surface_t
   use meds_fast_types,       only : fast_static_state_t, fast_dynamic_state_t, fast_forcing_t, fast_output_t,     &
                                     fast_workspace_t, step_report_t, step_exports_t, WORK_INTEG_REJECT,         &
                                     WORK_RK45_RESCUE
   use meds_fast_ark,         only : column_fast_step_ark
   use meds_fast_rk45,        only : column_fast_step_rk45, rk45_state_railed, rk45_ledgers_stop
   use meds_plant_hydraulics, only : leaf_node_psi
   implicit none
   private

   public :: fast_step

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one dt_fast with the configured scheme ([fast].time_integrator):                 !
   !   "ark" (default) -- the 2-solve ESDIRK2, gamma = 1 - 1/sqrt(2); not IMEX, since the     !
   !                      biotic CO2 source is implicit and the explicit tableau is empty     !
   !                      (docs/science/numerical_scheme.md eq. (1)).                         !
   !   "rk45"          -- the adaptive explicit Cash-Karp march, a full alternative a run can  !
   !                      choose and the explicit reference the tests compare ARK against.    !
   !---------------------------------------------------------------------------------------!
   subroutine fast_step(dt_fast, cfg, site, static, forcing, dynamic, out, work)
      real(wp),                   intent(in)    :: dt_fast
      type(meds_config_t),        intent(in)    :: cfg
      type(site_surface_t),       intent(in)    :: site      !< the patch's site: its soil and its ground
      type(fast_static_state_t),  intent(in)    :: static
      type(fast_forcing_t),       intent(in)    :: forcing
      type(fast_dynamic_state_t), intent(inout) :: dynamic
      type(fast_output_t),        intent(inout) :: out       !< sized by the caller (ensure_output_capacity)
      type(fast_workspace_t),     intent(inout) :: work      !< this thread's step storage, reused step to step
      integer(ik) :: jcoh

      out%report%work(WORK_RK45_RESCUE) = 0_ik
      if (out%capture) out%cdiag(:, 1:static%n) = 0.0_wp

      !----- Each cohort's psi_leaf for the daily-maximum accumulator behind the next day's stomatal  !
      !      control (docs/science/leaf_gas_exchange.md §4.2), at state^n, the value the leaf kernel   !
      !      was handed. Before the dispatch, so both schemes fill it (RK45's success path returns). -!
      do jcoh = 1_ik, static%n
         out%psi_leaf_coh(jcoh) = leaf_node_psi(dynamic%leaf_water_mass(jcoh), static%bleaf(jcoh),              &
               dynamic%wood_water_mass(jcoh), static%bwood_water(jcoh),                                            &
               static%height(jcoh), cfg%pft%hydraulics%pft(static%pft(jcoh)))
      end do

      if (cfg%time_integrator == INTEG_RK45) then
         !----- RK45 is explicit over the whole column. On a dense, cold canopy the leaf <-> canopy-air !
         !      coupling can be too stiff even at the sub-step floor: the state rails, and its 5th- and  !
         !      4th-order solutions rail together, so the controller would accept it. The rescue        !
         !      (docs/science/numerical_scheme.md §3, `rk45`): keep state^n, take the explicit step,    !
         !      and redo the dt_fast on ARK, whose Newton solves that coupled pair, when the march      !
         !      gave up (stiff_bail) or committed a railed state (rk45_state_railed). -----------------!
         block
            type(fast_dynamic_state_t) :: dynamic_save
            type(step_report_t)        :: report_save
            logical                    :: rk45_stiff
            dynamic_save = dynamic
            report_save = out%report
            call column_fast_step_rk45(dt_fast, cfg, site, static, forcing, dynamic, out, work, rk45_stiff)
            if (.not. rk45_stiff .and. .not. rk45_state_railed(dynamic, site%soil_profile%n_active)) then
               call rk45_ledgers_stop(out%report, dt_fast, cfg%debug_error)   ! kept, so a breach may now stop the run
               call atm_fluxes(out%report%exports, dt_fast, out%le_flux, out%h_flux)
               return
            end if
            !----- Roll back the whole dt_fast rather than resume from the failure: rescues are rare, a  !
            !      resume would sum two schemes' boundary fluxes in one step's ledger                   !
            !      (docs/science/numerical_scheme.md §5b), and a railed final state has no sub-step to  !
            !      resume from. -------------------------------------------------------------------------!
            dynamic = dynamic_save
            out%report = report_save   ! discard the railed/bailed RK45 step (rollback)
            out%report%work(WORK_INTEG_REJECT) = out%report%work(WORK_INTEG_REJECT) + 1_ik   ! a rejected step
            out%report%work(WORK_RK45_RESCUE)  = out%report%work(WORK_RK45_RESCUE) + 1_ik    ! ...and a rescue
         end block
      end if

      !----- ARK (ESDIRK2): the default, and the RK45 rescue target. ----------------------------!
      call column_fast_step_ark(dt_fast, cfg, site, static, forcing, dynamic, out, work)
      call atm_fluxes(out%report%exports, dt_fast, out%le_flux, out%h_flux)
   end subroutine fast_step

   !---------------------------------------------------------------------------------------!
   ! atm_fluxes -- the canopy-air -> atmosphere fluxes the run reports, from the ledger's own  !
   ! boundary export over the step, so the reported ET is the conserved vapour export:         !
   ! LE = L_v * vapour export / dt and H = sensible export / dt, the b-weighted              !
   ! g_atm_heat*cp_air*(T_cas - theta_atm), not the enthalpy export minus LE                   !
   ! (docs/science/numerical_scheme.md §5b).                                                    !
   !---------------------------------------------------------------------------------------!
   pure subroutine atm_fluxes(exports, dt_fast, le_flux, h_flux)
      type(step_exports_t), intent(in)  :: exports
      real(wp),             intent(in)  :: dt_fast
      real(wp),             intent(out) :: le_flux, h_flux   !< [W/m2] step-mean CAS -> atmosphere
      le_flux = exports%atm_vap_export  / dt_fast * latent_heat_vap
      h_flux  = exports%atm_heat_export / dt_fast
   end subroutine atm_fluxes

end module meds_fast_step
