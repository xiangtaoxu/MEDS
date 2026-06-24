!----- Full spin-up (daily AND monthly): area conserved, no NaN, bounded counts, BA>0. ------!
program test_spinup
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : yr_day
   use meds_config,          only : meds_config_t, build_config, TS_DAILY, TS_MONTHLY
   use meds_demography_types,           only : community, community_free
   use meds_rates_empirical, only : empirical_provider_t, make_empirical_provider
   use meds_setup,           only : init_bare_ground
   use meds_step,            only : advance_one_step
   use meds_cohort_dynamics, only : max_ccount
   use meds_diagnostics,     only : total_area, total_basal_area, total_nplant, has_nan
   use meds_test_support,    only : check, banner
   implicit none

   call banner('full spin-up conservation (daily + monthly)')
   call run_and_check(build_config(ts_mode = TS_DAILY),   40_ik)
   call run_and_check(build_config(ts_mode = TS_MONTHLY), 40_ik)
   write(*,'(a)') '   PASS'

contains

   subroutine run_and_check(cfg, nyears)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: nyears
      type(community)            :: comm
      type(empirical_provider_t) :: prov
      integer(ik) :: spm, istep, month, prev_month, diy
      real(wp)    :: ba
      logical     :: nm, ny

      prov = make_empirical_provider(cfg%pft)
      call init_bare_ground(comm, cfg, 4_ik, 298.15_wp, 295.15_wp)
      if (cfg%ts_mode == TS_MONTHLY) then
         spm = 12_ik
      else
         spm = nint(yr_day, ik)
      end if
      prev_month = 1_ik ; diy = 0_ik

      do istep = 1_ik, nyears * spm
         if (cfg%ts_mode == TS_MONTHLY) then
            month = mod(istep - 1_ik, 12_ik) + 1_ik
            nm    = .true.
            ny    = (month == 1_ik)
         else
            diy   = diy + 1_ik
            month = min(12_ik, (diy - 1_ik) * 12_ik / spm + 1_ik)
            nm    = (month /= prev_month) .or. (istep == 1_ik)
            ny    = (diy == 1_ik)
            if (diy >= spm) diy = 0_ik
            prev_month = month
         end if

         call advance_one_step(comm, prov, cfg, nm, ny)

         if (ny) then
            call check(abs(total_area(comm) - 1.0_wp) < 1.0e-5_wp, 'site area drifted from 1')
            call check(.not. has_nan(comm), 'NaN in state during spin-up')
            call check(max_ccount(comm) <= abs(cfg%maxcohort), 'per-patch cohort count exceeded cap')
            ba = total_basal_area(comm)
            call check(ba == ba .and. ba >= 0.0_wp, 'basal area non-finite or negative')
            call check(total_nplant(comm) >= 0.0_wp, 'negative plant number')
         end if
      end do

      !----- After a long spin-up the stand should be populated. -------------------------!
      call check(total_basal_area(comm) > 0.1_wp, 'stand failed to develop basal area')
      call community_free(comm)
   end subroutine run_and_check

end program test_spinup
