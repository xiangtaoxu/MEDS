!==========================================================================================!
! meds_competition -- the law-free STRUCTURAL competition sweep. It fills the stored per-cohort !
! `overtopping_lai` diagnostic: for each patch, a single height-sorted top-down pass accumulates !
! the overtopping leaf-area index (the cumulative LAI of all TALLER cohorts), with equal-height  !
! cohorts sharing one layer (co-dominance). This is a pure GEOMETRIC reduction over the cohort   !
! SoA + patch CSR -- NO rate law, NO plant physiology -- so it lives in the demography engine     !
! (⊥ plant). It was lifted verbatim out of the former empirical_vital_rates sweep; the driver    !
! and the Python empirical example read `overtopping_lai` as the Beer competition context.        !
!                                                                                          !
! ASSUMES the cohorts are height-sorted (descending) within each patch -- the fuse/fiss           !
! sort_cohorts invariant. Call it after any structural change (and once per slow step).           !
!==========================================================================================!
module meds_competition
   use meds_kinds,           only : wp, ik
   use meds_ecosystem_state, only : site_t
   implicit none
   private

   public :: compute_overtopping_lai

contains

   !---------------------------------------------------------------------------------------!
   ! Fill cohort%overtopping_lai for every cohort from the current (height-sorted) stand.    !
   !---------------------------------------------------------------------------------------!
   subroutine compute_overtopping_lai(site)
      type(site_t), intent(inout) :: site
      integer(ik) :: ip, i0, i1, i, j, k
      real(wp)    :: cum, over_lai, layer_lai

      associate (cohort => site%cohort, patch => site%patch)
         do ip = 1_ik, patch%n
            i0  = patch%cohort_offset(ip)
            i1  = i0 + patch%cohort_count(ip) - 1_ik
            cum = 0.0_wp                                   ! overtopping LAI accumulator
            i   = i0
            do while (i <= i1)
               !----- Equal-height cohorts are co-dominant (share the overtopping LAI). -----!
               k = i
               do while (k < i1)
                  if (cohort%height(k + 1_ik) /= cohort%height(i)) exit
                  k = k + 1_ik
               end do
               over_lai  = cum
               layer_lai = 0.0_wp
               do j = i, k
                  cohort%overtopping_lai(j) = over_lai
                  layer_lai = layer_lai + cohort%nplant(j) * cohort%leaf_area(j)
               end do
               cum = cum + layer_lai                        ! whole layer shades those below
               i   = k + 1_ik
            end do
         end do
      end associate
   end subroutine compute_overtopping_lai

end module meds_competition
