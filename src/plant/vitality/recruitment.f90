!==========================================================================================!
! meds_recruitment -- RECRUITMENT kernels of the phenomenological vital rates (incl. the       !
! reproduction flux; a dedicated reproduction module can split out later).                     !
!                                                                                          !
! Pure, stateless functions of (structure, traits) via the shared allometry. Two ingredients:  !
!   * the carbon of one minimum-size cohort (the recruit "unit"), per PFT; and                  !
!   * the reproduction recruit flux from one reproductive cohort -- the carbon diverted from     !
!     growth to reproduction over one year, in allometric AGB terms, converted to min-size       !
!     recruits. The baseline external seed rain is a plain per-PFT trait handled by the caller.  !
!==========================================================================================!
module meds_recruitment
   use meds_kinds,     only : wp
   use meds_allometry, only : dbh_to_height, dbh_to_agb, height_to_dbh
   implicit none
   private

   public :: min_cohort_carbon, reproduction_recruits

contains

   !----- Carbon of one minimum-size cohort [kgC/plant] (the recruit unit), for a given PFT.  !
   elemental pure function min_cohort_carbon(min_cohort_height, wood_density) result(carbon_min)
      real(wp), intent(in) :: min_cohort_height, wood_density
      real(wp)             :: carbon_min
      carbon_min = dbh_to_agb(height_to_dbh(min_cohort_height), min_cohort_height, wood_density)
   end function min_cohort_carbon

   !----- Reproduction recruit flux from one cohort [plant/m2/yr]: grow the diameter by the    !
   !       repro-allocated increment (repro_dbh, one year's worth), take the resulting AGB      !
   !       gain, and convert that reproductive carbon (times the establishment efficiency) into !
   !       minimum-size recruits. Caller has already gated this on the maturity height.  -------!
   elemental pure function reproduction_recruits(dbh, hgt_max, agb, nplant, wood_density,        &
                                repro_dbh, repro_carbon_efficiency, carbon_min) result(rec)
      real(wp), intent(in) :: dbh, hgt_max, agb, nplant, wood_density
      real(wp), intent(in) :: repro_dbh, repro_carbon_efficiency, carbon_min
      real(wp)             :: rec, new_dbh, dagb
      new_dbh = dbh + repro_dbh
      dagb    = dbh_to_agb(new_dbh, dbh_to_height(new_dbh, hgt_max), wood_density) - agb
      rec     = nplant * dagb * repro_carbon_efficiency / carbon_min
   end function reproduction_recruits

end module meds_recruitment
