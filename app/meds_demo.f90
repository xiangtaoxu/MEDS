!==========================================================================================!
! meds_demo -- standalone demographic spin-up driver.                                      !
!                                                                                          !
! Starts from near-bare ground and runs the demographic engine for N years with the test    !
! empirical rate provider, printing the evolving stand structure annually. Demonstrates the  !
! full pipeline: recruitment -> growth (size+competition) -> growth-dependent mortality ->   !
! cohort & patch fusion/fission. The provider is passed polymorphically, so swapping in a    !
! mechanistic rate model later requires no change here.                                     !
!                                                                                          !
! Time stepping is driven by an explicit calendar so daily and monthly modes share the       !
! exact same engine call.                                                                   !
!==========================================================================================!
program meds_demo
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : yr_day
   use meds_config,        only : meds_config_t, build_config, validate_config,            &
                                  TS_DAILY, TS_MONTHLY
   use meds_pft_params,    only : pft_table_t
   use meds_types,         only : community, community_free
   use meds_rates_empirical, only : empirical_provider_t, make_empirical_provider
   use meds_setup,         only : init_bare_ground
   use meds_step,          only : advance_one_step
   use meds_diagnostics,   only : print_summary, total_area, has_nan
   implicit none

   type(meds_config_t)        :: cfg
   type(community)            :: comm
   type(empirical_provider_t) :: prov
   integer(ik)                :: n_years, n_patch, steps_per_year, istep, total_steps
   integer(ik)                :: iyear, day_in_year, month, prev_month
   logical                    :: is_new_month, is_new_year
   real(wp)                   :: a0, a1
   character(len=32)          :: arg

   n_years = 60_ik
   n_patch = 6_ik

   !----- Optional first CLI argument overrides the number of years (for parity runs). ----!
   if (command_argument_count() >= 1_ik) then
      call get_command_argument(1, arg)
      read(arg, *) n_years
   end if

   cfg = build_config(ts_mode = TS_DAILY)
   call validate_config(cfg)
   prov = make_empirical_provider(cfg%pft)

   if (cfg%ts_mode == TS_MONTHLY) then
      steps_per_year = 12_ik
   else
      steps_per_year = nint(yr_day, ik)
   end if
   total_steps = n_years * steps_per_year

   !----- Warm tropical site (no frost; recruitment enabled). -----------------------------!
   call init_bare_ground(comm, cfg, n_patch, avg_temp = 298.15_wp, min_temp = 295.15_wp)
   a0 = total_area(comm)

   write(*,'(a)') '==================== MEDS standalone demographic spin-up ===================='
   write(*,'(a,i0,a,i0,a,i0,a)') ' years=', n_years, '  patches(init)=', n_patch,           &
                                 '  steps/yr=', steps_per_year, '  (provider: empirical)'
   write(*,'(a)') '-----------------------------------------------------------------------------'
   call print_summary(comm, 'year 0')

   prev_month  = 1_ik
   day_in_year = 0_ik
   iyear       = 0_ik

   do istep = 1_ik, total_steps
      if (cfg%ts_mode == TS_MONTHLY) then
         month        = mod(istep - 1_ik, 12_ik) + 1_ik
         is_new_month = .true.
         is_new_year  = (month == 1_ik)
      else
         day_in_year  = day_in_year + 1_ik
         month        = min(12_ik, (day_in_year - 1_ik) * 12_ik / steps_per_year + 1_ik)
         is_new_month = (month /= prev_month) .or. (istep == 1_ik)
         is_new_year  = (day_in_year == 1_ik)
         if (day_in_year >= steps_per_year) day_in_year = 0_ik
         prev_month   = month
      end if

      call advance_one_step(comm, prov, cfg, is_new_month, is_new_year)

      if (is_new_year) then
         iyear = iyear + 1_ik
         if (mod(iyear, 5_ik) == 0_ik .or. iyear == 1_ik) then
            call print_summary(comm, 'year '//itoa(iyear))
         end if
         if (has_nan(comm)) error stop 'meds_demo: NaN detected in state'
      end if
   end do

   call print_summary(comm, 'final')
   a1 = total_area(comm)
   write(*,'(a)') '-----------------------------------------------------------------------------'
   write(*,'(a,f12.9,a,f12.9)') ' site area start=', a0, '  end=', a1
   if (abs(a1 - 1.0_wp) > 1.0e-5_wp) error stop 'meds_demo: site area not conserved'
   write(*,'(a)') ' OK: spin-up completed, area conserved, no NaNs.'

   call community_free(comm)

contains

   function itoa(i) result(s)
      integer(ik), intent(in)       :: i
      character(len=:), allocatable :: s
      character(len=12)             :: buf
      write(buf,'(i0)') i
      s = trim(buf)
   end function itoa

end program meds_demo
