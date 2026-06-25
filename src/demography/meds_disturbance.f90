!==========================================================================================!
! meds_disturbance -- treefall patch disturbance (ED2 analogue).                           !
!                                                                                          !
! Once per structural (annual) step a fraction f = 1 - exp(-patch_disturbance_rate*dt) of    !
! EVERY patch's area is disturbed and aggregated into ONE new age-0 gap patch (dist_type      !
! DIST_TREEFALL). Survivorship follows ED2 treefall: tall canopy cohorts                       !
! (height >= disturbance_survive_height) are killed in the gap, while the short understory    !
! (height < threshold) survives into it. Each donor keeps area (1-f)*area with ALL its         !
! cohorts intact; the new patch receives the survivors at an area-weighted density.            !
!                                                                                          !
! Conservation: site area is exactly conserved (donors shed Sum f*area; the gap gains the     !
! same). Plant number is conserved for survivors and reduced by the disturbed fraction for     !
! the killed canopy -- that loss IS the disturbance. The gap is later consolidated by patch    !
! and cohort fusion (called after this routine in update_demography).                          !
!==========================================================================================!
module meds_disturbance
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   use meds_config,    only : meds_config_t, DIST_TREEFALL
   use meds_demography_types, only : site, patch_ensure_capacity, cohort_ensure_capacity,     &
                                     copy_cohort_slot, rebuild_csr
   use meds_sort,      only : sort_cohorts
   implicit none
   private

   public :: apply_patch_disturbance

contains

   subroutine apply_patch_disturbance(comm, cfg, dt_yr)
      type(site),          intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: dt_yr
      integer(ik) :: np0, newp, d, i, i0, i1, m, nsurv
      real(wp)    :: frac, new_area, wsum, tavg, tmin, w

      np0 = comm%pat%n
      if (np0 < 1_ik .or. cfg%patch_disturbance_rate <= 0.0_wp) return

      !----- Disturbed area fraction this step, and the aggregate new-patch area. ----------!
      frac = 1.0_wp - exp(-cfg%patch_disturbance_rate * dt_yr)
      if (frac <= tiny_num) return
      new_area = frac * sum(comm%pat%area(1:np0))
      if (new_area <= tiny_num) return

      !----- Count understorey survivors (all current cohorts live in donor patches). ------!
      nsurv = 0_ik
      do i = 1_ik, comm%coh%n
         if (comm%coh%height(i) < cfg%disturbance_survive_height) nsurv = nsurv + 1_ik
      end do

      call patch_ensure_capacity(comm%pat, np0 + 1_ik, comm%n_pft)
      call cohort_ensure_capacity(comm%coh, comm%coh%n + nsurv)
      newp = np0 + 1_ik

      associate (c => comm%coh, p => comm%pat)
         !----- New age-0 gap patch: area-weighted temperatures over the donors. -----------!
         wsum = 0.0_wp ; tavg = 0.0_wp ; tmin = huge(1.0_wp)
         do d = 1_ik, np0
            w    = frac * p%area(d)
            wsum = wsum + w
            tavg = tavg + w * p%avg_daily_temp(d)
            tmin = min(tmin, p%min_month_temp(d))
         end do
         tavg = tavg / max(wsum, tiny_num)
         p%area(newp)           = new_area
         p%age(newp)            = 0.0_wp
         p%dist_type(newp)      = DIST_TREEFALL
         p%avg_daily_temp(newp) = tavg
         p%min_month_temp(newp) = tmin
         p%recruit_pool(:,newp) = 0.0_wp
         p%n = newp

         !----- Move understorey survivors into the gap at area-weighted density. ----------!
         m = c%n
         do d = 1_ik, np0
            i0 = p%cohort_offset(d) ; i1 = i0 + p%cohort_count(d) - 1_ik
            do i = i0, i1
               if (c%height(i) >= cfg%disturbance_survive_height) cycle   ! canopy dies in gap
               m = m + 1_ik
               call copy_cohort_slot(c, m, i)
               c%nplant(m)      = c%nplant(i) * (frac * p%area(d) / new_area)
               c%owner_patch(m) = newp
            end do
         end do
         c%n = m

         !----- Donors keep the undisturbed remainder (all their cohorts intact). ----------!
         do d = 1_ik, np0
            p%area(d) = (1.0_wp - frac) * p%area(d)
         end do
      end associate

      call rebuild_csr(comm)
      call sort_cohorts(comm)
   end subroutine apply_patch_disturbance

end module meds_disturbance
