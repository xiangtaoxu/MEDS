!==========================================================================================!
! capture_golden -- P0 of the reorg: snapshot the CURRENT empirical demography spin-up as a  !
! reproducible golden trajectory. Mirrors test_spinup's deterministic run (build_test_config, !
! no forcing, GS_EMPIRICAL) at dt_slow = 1 day for 40 years and writes per-year site totals   !
! plus the final per-cohort state to CSV. The reorganized Fortran can no longer run empirical  !
! (it becomes a Python example), so this golden is the anchor the Python `example_demography`  !
! must reproduce to tolerance (reorg design doc S8).                                           !
!==========================================================================================!
program capture_golden
   use meds_kinds,                  only : wp, ik
   use meds_constants,              only : yr_day, day_sec
   use meds_config,                 only : meds_config_t
   use meds_demography_interface,   only : site_t
   use meds_demography_types,       only : site_free
   use meds_init,                   only : init_bare_ground
   use meds_stepper,                only : advance_one_step
   use meds_demography_diagnostics, only : total_basal_area, total_agb, total_lai, total_nplant, count_cohorts
   use meds_test_support,           only : build_test_config
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)        :: site
   integer(ik) :: yday, nday, istep, nsteps, step_days, month, prev_month, year, u, i
   logical     :: nm, ny

   cfg       = build_test_config(dt_slow = 1.0_wp * day_sec)
   call init_bare_ground(site, cfg, 4_ik)
   nday      = nint(yr_day, ik)
   step_days = max(1_ik, nint(cfg%dt_slow / day_sec, ik))
   nsteps    = (40_ik * nday) / step_days
   prev_month = 0_ik ; yday = 0_ik ; year = 0_ik

   open(newunit=u, file='test/golden/empirical_spinup_golden.csv', status='replace', action='write')
   write(u,'(a)') 'year,n_cohort,total_agb,total_lai,total_nplant,total_basal_area'

   do istep = 1_ik, nsteps
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
         year = year + 1_ik
         write(u,'(i0,",",i0,4(",",es24.16e3))') year, count_cohorts(site), &
              total_agb(site), total_lai(site), total_nplant(site), total_basal_area(site)
      end if
   end do
   close(u)

   !----- Final per-cohort state (a richer anchor: pft, dbh, height, nplant, agb). ---------!
   open(newunit=u, file='test/golden/empirical_spinup_final_cohorts.csv', status='replace', action='write')
   write(u,'(a)') 'cohort,pft,owner_patch,dbh,height,nplant,agb'
   do i = 1_ik, site%cohort%n
      write(u,'(i0,2(",",i0),4(",",es24.16e3))') i, site%cohort%pft(i), site%cohort%owner_patch(i), &
           site%cohort%dbh(i), site%cohort%height(i), site%cohort%nplant(i), site%cohort%agb(i)
   end do
   close(u)

   write(*,'(a,i0,a)') 'capture_golden: wrote ', year, ' years + final cohorts to test/golden/'
   call site_free(site)
end program capture_golden
