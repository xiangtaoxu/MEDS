! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_diagnostic_kernels -- stage [1] of the diagnostic output: quantities stored nowhere   !
! in the state but a closed-form function of it (cohort_lai, cohort_npp_per_plant,          !
! soil_wetness, dbh_class_index). Plain arguments only, so each routine is unit-testable.   !
! It holds only what no physics library owns, and sits above the physics, which cannot call !
! it (docs/science/diagnostics.md §1, Where the stages live).                               !
!==========================================================================================!
module meds_diagnostic_kernels
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : tiny_num
   implicit none
   private

   public :: cohort_lai, cohort_npp_per_plant, soil_wetness, dbh_class_index

contains

   !=======================================================================================!
   !  Structural / carbon                                                                   !
   !=======================================================================================!

   !----- Per-cohort leaf area index [m2 leaf / m2 ground of its own patch]. ----------------!
   elemental real(wp) function cohort_lai(nplant, leaf_area) result(lai)
      real(wp), intent(in) :: nplant     !< [plant/m2]
      real(wp), intent(in) :: leaf_area  !< [m2/plant]
      lai = nplant * leaf_area
   end function cohort_lai

   !----- Per-plant NPP [kgC/plant] over the accumulators' period: gross GPP minus the three   !
   !      maintenance respirations. Growth respiration is charged in the slow step's          !
   !      allocation (growth_resp_site), so subtracting it here would count it twice. --------!
   elemental real(wp) function cohort_npp_per_plant(gpp, leaf_resp, stem_resp, root_resp) &
                                    result(npp)
      real(wp), intent(in) :: gpp, leaf_resp, stem_resp, root_resp   !< [kgC/plant]
      npp = gpp - leaf_resp - stem_resp - root_resp
   end function cohort_npp_per_plant

   !=======================================================================================!
   !  Soil                                                                                  !
   !=======================================================================================!

   !----- Relative saturation [-]: where the layer sits between residual and saturation.     !
   !      Clamped to [0,1] so a solver excursion outside the constitutive domain reports as    !
   !      exactly dry or saturated rather than as an out-of-range number. --------------------!
   elemental real(wp) function soil_wetness(theta, theta_res, theta_sat) result(w)
      real(wp), intent(in) :: theta, theta_res, theta_sat
      w = (theta - theta_res) / max(theta_sat - theta_res, tiny_num)
      w = min(max(w, 0.0_wp), 1.0_wp)
   end function soil_wetness

   !=======================================================================================!
   !  Size-class binning (DIM_SIZE)                                                          !
   !=======================================================================================!

   !----- 1-based DBH class of a dbh (docs/science/diagnostics.md §3). Over n_class+1 ascending !
   !      edges, bin k is [edges(k), edges(k+1)); the last bin also takes every dbh above its  !
   !      top edge, and a dbh below edges(1) falls in bin 1, so no plant is dropped and        !
   !      sum(nplant_size) == nplant_site holds. ----------------------------------------------!
   pure integer(ik) function dbh_class_index(dbh, edges, n_class) result(k)
      real(wp),    intent(in) :: dbh
      real(wp),    intent(in) :: edges(:)   !< (n_class+1) ascending
      integer(ik), intent(in) :: n_class
      integer(ik) :: j
      k = 1_ik
      if (n_class <= 0_ik) return
      do j = n_class, 1_ik, -1_ik
         if (dbh >= edges(j)) then
            k = j
            return
         end if
      end do
   end function dbh_class_index

end module meds_diagnostic_kernels
