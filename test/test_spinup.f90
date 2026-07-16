!----- Full spin-up (1-day, 7-day, 28-day slow steps): area conserved, no NaN, bounded, BA>0. -!
program test_spinup
   use meds_kinds,                only : wp, ik
   use meds_constants,            only : yr_day, day_sec
   use meds_config,               only : meds_config_t
   use meds_demography_interface, only : site_t
   use meds_ecosystem_state,     only : site_free
   use meds_init,                 only : init_bare_ground
   use meds_stepper,              only : advance_one_step
   use meds_demography_fusefiss, only : max_cohort_count
   use meds_demography_diagnostics, only : total_area, total_basal_area, total_agb, total_lai, &
                                         total_nplant, has_nan
   use meds_test_support, only : build_test_config, check, banner
   implicit none

   call banner('full spin-up conservation (dt_slow = 1 d, 7 d, 28 d)')
   call run_and_check(build_test_config(dt_slow =  1.0_wp * day_sec), 40_ik)
   call run_and_check(build_test_config(dt_slow =  7.0_wp * day_sec), 40_ik)
   call run_and_check(build_test_config(dt_slow = 28.0_wp * day_sec), 40_ik)
   write(*,'(a)') '   PASS'

contains

   subroutine run_and_check(cfg, nyears)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: nyears
      type(site_t)  :: site
      integer(ik) :: yday, nday, istep, nsteps, step_days, month, prev_month
      real(wp)    :: ba
      logical     :: nm, ny

      call init_bare_ground(site, cfg, 4_ik)
      nday       = nint(yr_day, ik)
      step_days  = max(1_ik, nint(cfg%dt_slow / day_sec, ik))
      nsteps     = (nyears * nday) / step_days
      prev_month = 0_ik
      yday       = 0_ik

      do istep = 1_ik, nsteps
         !----- Fold the day counter into the month/year fuse-fiss cadence flags. -----------!
         yday = yday + step_days
         ny   = .false.
         if (yday > nday) then
            yday = yday - nday
            ny   = .true.
         end if
         month = min(12_ik, (yday - 1_ik) * 12_ik / nday + 1_ik)
         nm    = (month /= prev_month) .or. ny .or. (istep == 1_ik)
         if (istep == 1_ik) ny = .true.
         prev_month = month

         call advance_one_step(site, cfg, nm, ny)

         if (ny) then
            call check(abs(total_area(site) - 1.0_wp) < 1.0e-5_wp, 'site_t area drifted from 1')
            call check(.not. has_nan(site), 'NaN in state during spin-up')
            call check(max_cohort_count(site) <= abs(cfg%max_cohort), 'per-patch cohort count exceeded cap')
            ba = total_basal_area(site)
            call check(ba == ba .and. ba >= 0.0_wp, 'basal area non-finite or negative')
            call check(total_nplant(site) >= 0.0_wp, 'negative plant number')
         end if
      end do

      !----- After a long spin-up the stand should be populated (leaf area + carbon). -----!
      call check(total_lai(site) > 0.1_wp, 'stand failed to develop leaf area')
      call check(total_agb(site) > 0.0_wp, 'stand failed to accumulate biomass')
      call site_free(site)
   end subroutine run_and_check

end program test_spinup
