!==========================================================================================!
! meds_demo -- standalone demographic spin-up driver (a minimal master stepper).           !
!                                                                                          !
! Starts from near-bare ground and runs the demographic engine for N years, printing the     !
! evolving stand structure annually. Demonstrates the full pipeline: recruitment -> growth   !
! (light competition via overtopping LAI) -> shade-driven mortality -> cohort & patch         !
! fusion/fission -> treefall disturbance. The per-step rates come from the empirical          !
! evaluator (meds_rates_empirical); demography is driven only through meds_stepper -> the     !
! meds_demography_interface seam.                                                            !
!                                                                                          !
! An explicit calendar lets daily / weekly / monthly modes share the exact same engine call. !
! Usage:  meds_demo [n_years] [daily|weekly|monthly]                                          !
!==========================================================================================!
program meds_demo
   use meds_kinds,                only : wp, ik
   use meds_constants,            only : yr_day
   use meds_config,               only : meds_config_t, build_config, validate_config,      &
                                         TS_DAILY, TS_WEEKLY, TS_MONTHLY
   use meds_demography_interface, only : site
   use meds_demography_types,     only : site_free
   use meds_setup,                only : init_bare_ground
   use meds_stepper,              only : advance_one_step
   use meds_diagnostics,          only : print_summary, total_area, has_nan
   implicit none

   type(meds_config_t) :: cfg
   type(site)          :: comm
   integer(ik)         :: n_years, n_patch, ts_mode, steps_per_year, nday, step_days
   integer(ik)         :: istep, nsteps, iyear, yday, month, prev_month
   logical             :: is_new_month, is_new_year
   real(wp)            :: a0, a1
   character(len=32)   :: arg

   n_years = 60_ik
   n_patch = 6_ik
   ts_mode = TS_DAILY

   !----- Optional CLI: arg1 = number of years; arg2 = daily|weekly|monthly. --------------!
   if (command_argument_count() >= 1_ik) then
      call get_command_argument(1, arg) ; read(arg, *) n_years
   end if
   if (command_argument_count() >= 2_ik) then
      call get_command_argument(2, arg)
      select case (trim(arg))
      case ('weekly') ; ts_mode = TS_WEEKLY
      case ('monthly'); ts_mode = TS_MONTHLY
      case default    ; ts_mode = TS_DAILY
      end select
   end if

   cfg = build_config(ts_mode = ts_mode)
   call validate_config(cfg)

   nday = nint(yr_day, ik)
   select case (cfg%ts_mode)
   case (TS_MONTHLY) ; step_days = 0_ik ; steps_per_year = 12_ik
   case (TS_WEEKLY)  ; step_days = 7_ik ; steps_per_year = nday / step_days
   case default      ; step_days = 1_ik ; steps_per_year = nday
   end select
   nsteps = n_years * steps_per_year

   !----- Warm tropical site (no frost; recruitment enabled). -----------------------------!
   call init_bare_ground(comm, cfg, n_patch, avg_temp = 298.15_wp, min_temp = 295.15_wp)
   a0 = total_area(comm)

   write(*,'(a)') '==================== MEDS standalone demographic spin-up ===================='
   write(*,'(a,i0,a,i0,a,i0)') ' years=', n_years, '  patches(init)=', n_patch,             &
                               '  steps/yr=', steps_per_year
   write(*,'(a)') '-----------------------------------------------------------------------------'
   call print_summary(comm, 'year 0')

   prev_month = 0_ik
   yday       = 0_ik
   iyear      = 0_ik

   do istep = 1_ik, nsteps
      if (cfg%ts_mode == TS_MONTHLY) then
         month        = mod(istep - 1_ik, 12_ik) + 1_ik
         is_new_month = .true.
         is_new_year  = (month == 1_ik)
      else
         yday        = yday + step_days
         is_new_year = .false.
         if (yday > nday) then
            yday        = yday - nday
            is_new_year = .true.
         end if
         month        = min(12_ik, (yday - 1_ik) * 12_ik / nday + 1_ik)
         is_new_month = (month /= prev_month) .or. is_new_year .or. (istep == 1_ik)
         if (istep == 1_ik) is_new_year = .true.
         prev_month   = month
      end if

      call advance_one_step(comm, cfg, is_new_month, is_new_year)

      if (is_new_year) then
         iyear = iyear + 1_ik
         if (mod(iyear, 5_ik) == 0_ik .or. iyear == 1_ik) call print_summary(comm, 'year '//itoa(iyear))
         if (has_nan(comm)) error stop 'meds_demo: NaN detected in state'
      end if
   end do

   call print_summary(comm, 'final')
   a1 = total_area(comm)
   write(*,'(a)') '-----------------------------------------------------------------------------'
   write(*,'(a,f12.9,a,f12.9)') ' site area start=', a0, '  end=', a1
   if (abs(a1 - 1.0_wp) > 1.0e-5_wp) error stop 'meds_demo: site area not conserved'
   write(*,'(a)') ' OK: spin-up completed, area conserved, no NaNs.'

   call site_free(comm)

contains

   function itoa(i) result(s)
      integer(ik), intent(in)       :: i
      character(len=:), allocatable :: s
      character(len=12)             :: buf
      write(buf,'(i0)') i
      s = trim(buf)
   end function itoa

end program meds_demo
