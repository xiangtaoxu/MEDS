!----- Full spin-up (daily, weekly, monthly): area conserved, no NaN, bounded counts, BA>0. --!
program test_spinup
   use meds_kinds,                only : wp, ik
   use meds_constants,            only : yr_day
   use meds_config,               only : meds_config_t, build_config, TS_DAILY, TS_WEEKLY, TS_MONTHLY
   use meds_demography_interface, only : site
   use meds_demography_types,     only : site_free
   use meds_setup,                only : init_bare_ground
   use meds_step,                 only : advance_one_step
   use meds_cohort_dynamics,      only : max_cohort_count
   use meds_diagnostics,          only : total_area, total_basal_area, total_nplant, has_nan
   use meds_test_support,         only : check, banner
   implicit none

   call banner('full spin-up conservation (daily + weekly + monthly)')
   call run_and_check(build_config(ts_mode = TS_DAILY),   40_ik)
   call run_and_check(build_config(ts_mode = TS_WEEKLY),  40_ik)
   call run_and_check(build_config(ts_mode = TS_MONTHLY), 40_ik)
   write(*,'(a)') '   PASS'

contains

   subroutine run_and_check(cfg, nyears)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: nyears
      type(site)  :: comm
      integer(ik) :: yday, nday, istep, nsteps, step_days, month, prev_month
      real(wp)    :: ba
      logical     :: nm, ny

      call init_bare_ground(comm, cfg, 4_ik, 298.15_wp, 295.15_wp)
      nday       = nint(yr_day, ik)
      prev_month = 0_ik
      yday       = 0_ik

      select case (cfg%ts_mode)
      case (TS_MONTHLY)
         nsteps = nyears * 12_ik
      case (TS_WEEKLY)
         step_days = 7_ik ; nsteps = (nyears * nday) / step_days
      case default
         step_days = 1_ik ; nsteps = nyears * nday
      end select

      do istep = 1_ik, nsteps
         if (cfg%ts_mode == TS_MONTHLY) then
            month = mod(istep - 1_ik, 12_ik) + 1_ik
            nm    = .true.
            ny    = (month == 1_ik)
         else
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
         end if

         call advance_one_step(comm, cfg, nm, ny)

         if (ny) then
            call check(abs(total_area(comm) - 1.0_wp) < 1.0e-5_wp, 'site area drifted from 1')
            call check(.not. has_nan(comm), 'NaN in state during spin-up')
            call check(max_cohort_count(comm) <= abs(cfg%max_cohort), 'per-patch cohort count exceeded cap')
            ba = total_basal_area(comm)
            call check(ba == ba .and. ba >= 0.0_wp, 'basal area non-finite or negative')
            call check(total_nplant(comm) >= 0.0_wp, 'negative plant number')
         end if
      end do

      !----- After a long spin-up the stand should be populated. -------------------------!
      call check(total_basal_area(comm) > 0.1_wp, 'stand failed to develop basal area')
      call site_free(comm)
   end subroutine run_and_check

end program test_spinup
