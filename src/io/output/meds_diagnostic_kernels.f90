! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_diagnostic_kernels -- pure/elemental DERIVED diagnostic quantities.                  !
!                                                                                          !
! Stage [1] of the diagnostic wall (docs/science/diagnostics.md §1): quantities that are not        !
! stored anywhere in state but are a closed-form function of things that are. Plain arguments      !
! only -- no polygon_t, no config aggregator -- so every routine here is directly unit-testable.   !
!                                                                                          !
! CONTENT RULE. Anything a physics library already owns is CALLED, never re-implemented:   !
! the soil's water potential lives in meds_soil_hydr_lib, tissue_water_to_psi in           !
! meds_plant_hydr_lib, internal_energy_to_temp and the moist-air thermodynamics in         !
! meds_therm_lib, the allometry in meds_allometry_lib. This module holds only diagnostics  !
! that no physics kernel owns. Where a thin re-export is more readable than making every   !
! caller `use` two modules, the wrapper is marked as such and simply forwards.             !
!                                                                                          !
! PLACEMENT. src/io/ (the meds_io_prep target, which links meds_demography only -- no netCDF), NOT  !
! src/share/functions/. Accepted consequence: nothing BELOW io in the library DAG can call these.    !
! That is correct rather than merely tolerated -- a diagnostic is by definition something no        !
! physics kernel needs, and putting it below the wall would be the first step toward physics        !
! depending on its own reporting layer. If a quantity here ever turns out to be needed by physics,  !
! that is the signal it was never a diagnostic and belongs in the owning physics library instead.   !
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

   !----- Per-plant NPP [kgC/plant] over whatever period the accumulators span: GROSS GPP    !
   !      minus the three autotrophic MAINTENANCE respiration terms. Growth respiration is     !
   !      NOT included here -- it is charged inside the slow-loop allocation step, and adding   !
   !      it would double-count against npp_*_site. Mirrors total_npp's definition exactly.     !
   elemental real(wp) function cohort_npp_per_plant(gpp, leaf_resp, stem_resp, root_resp) &
                                    result(npp)
      real(wp), intent(in) :: gpp, leaf_resp, stem_resp, root_resp   !< [kgC/plant]
      npp = gpp - leaf_resp - stem_resp - root_resp
   end function cohort_npp_per_plant

   !=======================================================================================!
   !  Soil                                                                                  !
   !=======================================================================================!

   !----- Relative saturation (ED2 `soil_wetness`) [-]: where the layer sits between residual !
   !      and saturation. Clamped to [0,1] so a solver excursion outside the constitutive       !
   !      domain reports as exactly dry/saturated rather than as an out-of-range number.        !
   elemental real(wp) function soil_wetness(theta, theta_res, theta_sat) result(w)
      real(wp), intent(in) :: theta, theta_res, theta_sat
      w = (theta - theta_res) / max(theta_sat - theta_res, tiny_num)
      w = min(max(w, 0.0_wp), 1.0_wp)
   end function soil_wetness

   !=======================================================================================!
   !  Size-class binning (DIM_SIZE)                                                          !
   !=======================================================================================!

   !----- 1-based DBH class of one plant, ED2-style (docs/science/diagnostics.md §3).        !
   !                                                                                          !
   !      `edges` holds n_class+1 ASCENDING boundaries. Bin k is the half-open interval          !
   !      [edges(k), edges(k+1)), EXCEPT the last, which is closed at the top so the largest      !
   !      tree in the stand is never dropped. A dbh below edges(1) is clamped into bin 1 rather   !
   !      than discarded -- silently losing the smallest cohorts would break the closure           !
   !      sum(nplant_size) == nplant_site that makes this axis trustworthy.                       !
   !                                                                                          !
   !      A cohort is NEVER split across bins: it lands whole in the bin of its mean dbh. That     !
   !      is the standard ED2 treatment, and it is what makes nplant_size a genuine stem-density   !
   !      distribution comparable to a forest inventory.                                          !
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
