! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_biogeochem_dynamics -- the slow loop's soil biogeochemistry, run after the          !
! vegetation driver (docs/science/order_of_processes.md §4). advance_biogeochem_dynamics   !
! advances each patch's CENTURY pools one daily step: the step's litter                    !
! (poly%patch%litter_in) and the fast loop's decomposition integral (poly%patch%decomp_xi) !
! drive soil_carbon_step (docs/science/soil_carbon.md §5-§6). The fast loop respires a     !
! frozen daily copy and never writes the pools, so the fast/slow seam gap is an assertion  !
! guard (§7).                                                                              !
!==========================================================================================!
module meds_biogeochem_dynamics
   use meds_kinds,            only : wp, ik
   use meds_config,           only : meds_config_t
   use meds_polygon_state_vars, only : polygon_t, patch_polygon_area
   use meds_biogeochem_types, only : litter_input_t, soilc_audit_t, soilc_seam_t
   use meds_biogeochem_params, only : N_SOIL_POOL
   use meds_soil_biogeochem,  only : build_litter_input, soil_carbon_step
   use meds_slow_ledger,      only : slow_ledger_t, slow_ledger_declare
   implicit none
   private

   public :: advance_biogeochem_dynamics

contains

   !---------------------------------------------------------------------------------------!
   ! Advance each patch's soil-carbon pools one daily step on the step's litter,           !
   ! poly%patch%litter_in (turnover and mortality; the cull and the disturbance put theirs !
   ! on the pools directly). `seam` keeps the run's worst Rh gap, lignin residual and      !
   ! lambda (soil_carbon.md §6-§7, §9).                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine advance_biogeochem_dynamics(poly, cfg, seam, ledger)
      type(polygon_t),       intent(inout) :: poly
      type(meds_config_t),   intent(in)    :: cfg
      type(soilc_seam_t), optional, intent(inout) :: seam   !< the run's worst fast/slow seam checks (soil_carbon.md §7)
      !----- The slow ledger: the litter declared leaving the live pools is declared       !
      !      arriving here, and the day's Rh is declared leaving (conservation.md §3). ----!
      type(slow_ledger_t), optional, intent(inout) :: ledger
      real(wp)             :: u(N_SOIL_POOL), lignin_in(2), rh_today
      type(soilc_audit_t)  :: audit
      integer(ik)          :: ip

      if (present(ledger)) then
         do ip = 1_ik, poly%patch%n
            associate (lit => poly%patch%litter_in(ip))
               call slow_ledger_declare(ledger, carbon_in = patch_polygon_area(poly, ip)           &
                        * (lit%labile_grnd + lit%labile_soil + lit%struct_grnd + lit%struct_soil))
            end associate
         end do
      end if
      do ip = 1_ik, poly%patch%n
         call build_litter_input(poly%patch%litter_in(ip), u, lignin_in)
         associate (xa => poly%patch%decomp_xi(ip))
            call soil_carbon_step(poly%patch%soil_carbon(ip), u, lignin_in, xa%xi_days, cfg%soil_carbon, &
                                  rh_today, audit)
            !----- The assertion guard (soil_carbon.md §7): rh_today, the matrix's own     !
            !      debit, should equal the fast loop's accumulated Rh. --------------------!
            audit%rh_fast_accum = xa%rh_fast
            audit%rh_seam_gap   = rh_today - xa%rh_fast
            if (present(seam)) then
               seam%worst_rh_gap = max(seam%worst_rh_gap, abs(audit%rh_seam_gap))
               seam%worst_lignin = max(seam%worst_lignin, abs(audit%lignin_resid))
               if (audit%lambda_max > seam%worst_lambda) then
                  seam%worst_lambda = audit%lambda_max
                  seam%lambda_pool  = audit%lambda_pool
               end if
            end if
            !----- Declare the matrix's own Rh debit (conservation.md §3): what is left in !
            !      the residual is the fast/slow seam gap. --------------------------------!
            if (present(ledger)) call slow_ledger_declare(ledger,                                  &
                                          carbon_out = patch_polygon_area(poly, ip) * rh_today)
         end associate
      end do
   end subroutine advance_biogeochem_dynamics

end module meds_biogeochem_dynamics
