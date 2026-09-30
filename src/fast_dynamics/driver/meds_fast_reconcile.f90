! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_reconcile -- make state CONSISTENT before the fast loop marches on it.            !
!                                                                                          !
! One reconciliation lives here today: each cohort's stored tissue water against the CAPACITY  !
! its current biomass allows. The module is named for the act rather than the quantity because  !
! this is where any future "repair the state before the march" step belongs -- and named        !
! `reconcile`, not `check`, because everything here WRITES. A name promising a check invites the !
! next reader to skip or reorder it, and these are mass edits.                                   !
!                                                                                          !
! WHY THIS IS NOT IN THE FAST LOOP. Tissue water capacity is a function of leaf, sapwood and  !
! fine-root carbon, all of which only change in the SLOW loop. The fast loop nevertheless ran  !
! this test inside its per-patch gather, on every cohort, on every sub-step -- and it is not a  !
! read: both branches WRITE the cohort's water back. So an unbooked mass edit lived in the      !
! inner loop of the integrator, where the fast loop's own whole-column ledger cannot see it      !
! (that ledger spans one dt_fast and opens after the gather).                                    !
!                                                                                          !
! Running it once per slow step is equivalent, not an approximation: within a slow step the      !
! capacity is constant, so the clamp is idempotent after its first application and the seed can   !
! fire at most once. What changes is that the edit happens in ONE place, at a declared point in    !
! the cadence, and reports the mass it moved.                                                      !
!                                                                                          !
! WHY IT SITS IN THE FAST DRIVER AND NOT THE SLOW ONE. The structure plan files it under the slow  !
! loop (§10.1), which is right about ownership -- the EVENT is a biomass change, and the slow-loop  !
! water ledger of §10.2 is where the seeded and discarded mass should eventually be booked. But     !
! the CONSUMER is the fast loop, and a caller that runs the fast loop without a preceding slow one  !
! (a test, a probe, a future C-API fast verb) reads zero tissue water with no error at all -- which  !
! is what happened to test_fast_loop the moment this left the gather. So it is called at the top of !
! `fast_dynamics`, once per call, where it cannot be skipped. When the slow ledger lands, the call   !
! moves and these two outputs become ledger terms.                                                   !
!                                                                                          !
! LEAF AND WOOD ARE INDEPENDENT (2026-09 review item 1B #3). One shared `leaf_water_mass <= 0`    !
! test used to re-seed BOTH stores, so a dormant deciduous cohort -- whose leaf water is legitimately!
! zero after the snap-to-bare shed -- had yesterday's integrated WOOD water overwritten by the      !
! seed every day of dormancy, and never carried a water deficit through winter.                     !
!==========================================================================================!
module meds_fast_reconcile
   use meds_kinds,            only : wp, ik
   use meds_site_state_types, only : site_t
   use meds_plant_types,      only : hydro_params_table_t
   use meds_column_params, only : PSI_INIT
   use meds_water_retention,  only : water_content, clamp_water_to_capacity
   implicit none
   private

   public :: reconcile_tissue_water_capacity

contains

   !---------------------------------------------------------------------------------------!
   ! Seed an empty tissue store to its PSI_INIT-equivalent mass, and clamp a full one to the  !
   ! capacity today's biomass allows. Returns the mass each branch moved [kg/m2 ground], which !
   ! is a SOURCE and a SINK respectively -- neither is currently booked anywhere.              !
   !                                                                                          !
   ! The clamp is needed because mass, not potential, is the slow/fast seam-continuous          !
   ! quantity: yesterday's water carries forward unchanged into today's grown biomass, which     !
   ! simply reads as a slightly lower relative water content. The only unreachable state is       !
   ! water ABOVE saturation, which a discontinuous biomass SHRINK can produce -- the phenology     !
   ! dormant-canopy snap-to-bare drops leaf carbon far enough in one step to do it.                !
   !---------------------------------------------------------------------------------------!
   subroutine reconcile_tissue_water_capacity(site, hydraulics_table, seeded, discarded)
      type(site_t),               intent(inout) :: site
      type(hydro_params_table_t), intent(in)    :: hydraulics_table   !< each PFT's tissue curves
      real(wp), optional,  intent(out)   :: seeded      !< [kg/m2] water CREATED by the lazy seed
      real(wp), optional,  intent(out)   :: discarded   !< [kg/m2] water DESTROYED by the clamp
      real(wp)    :: w_new, w_seed, w_lost, wood_c
      integer(ik) :: i
      w_seed = 0.0_wp ; w_lost = 0.0_wp
      associate (c => site%cohort)
         do i = 1_ik, c%n
            associate (leaf_curve => hydraulics_table%pft(c%pft(i))%leaf_curve,                    &
                       wood_curve => hydraulics_table%pft(c%pft(i))%wood_curve)
               !----- WOOD: the sapwood ring plus the fine roots share one store. ---------------!
               wood_c = c%sapwood_carbon(i) + c%fineroot_carbon(i)
               if (c%wood_water_mass(i) <= 0.0_wp) then
                  w_new = water_content(PSI_INIT, wood_curve, wood_c)
                  w_seed = w_seed + c%nplant(i) * w_new
               else
                  w_new = clamp_water_to_capacity(c%wood_water_mass(i), wood_curve, wood_c)
                  w_lost = w_lost + c%nplant(i) * (c%wood_water_mass(i) - w_new)
               end if
               c%wood_water_mass(i) = w_new
               !----- LEAF: tested independently of wood, see the module header. ----------------!
               if (c%leaf_water_mass(i) <= 0.0_wp) then
                  w_new = water_content(PSI_INIT, leaf_curve, c%leaf_carbon(i))
                  w_seed = w_seed + c%nplant(i) * w_new
               else
                  w_new = clamp_water_to_capacity(c%leaf_water_mass(i), leaf_curve, c%leaf_carbon(i))
                  w_lost = w_lost + c%nplant(i) * (c%leaf_water_mass(i) - w_new)
               end if
               c%leaf_water_mass(i) = w_new
            end associate
         end do
      end associate
      if (present(seeded))    seeded    = w_seed
      if (present(discarded)) discarded = w_lost
   end subroutine reconcile_tissue_water_capacity

end module meds_fast_reconcile
