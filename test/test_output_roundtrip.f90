!==========================================================================================!
! test_output_roundtrip -- the netCDF serializer end to end: build a manager, walk a 3-day        !
! calendar folding site state, flush per-stream files, then RE-OPEN them and assert record counts,  !
! averaged values, calendar, and CF layout. MEDS_IO_DESIGN.md tests 5 & 6 (P0 subset).             !
!==========================================================================================!
module test_ro_support
   use meds_kinds,            only : wp, ik
   use meds_core_state_types, only : site_t
   implicit none
contains
   !----- Opaque per-cohort agb mutation (production mutates through veg dynamics; see the        !
   !      test_output_integrate header note on the nvfortran -O2 aliasing hazard).                !
   subroutine set_site_agb(site, a)
      type(site_t), intent(inout) :: site
      real(wp),     intent(in)    :: a
      site%cohort%agb(1) = a
   end subroutine set_site_agb

   !----- Opaque setters for the P1 flux/state sources (fast-loop accumulators + soil column). ----!
   subroutine set_site_gpp(site, g)
      type(site_t), intent(inout) :: site
      real(wp),     intent(in)    :: g
      site%cohort%gpp_accum(1) = g
   end subroutine set_site_gpp

   subroutine set_site_soil_temp(site, t)
      type(site_t), intent(inout) :: site
      real(wp),     intent(in)    :: t
      site%patch%soil_e(1)%soil_temp(:) = t
   end subroutine set_site_soil_temp
end module test_ro_support

program test_output_roundtrip
   use iso_c_binding,         only : c_int, c_size_t, c_double
   use meds_kinds,            only : wp, ik
   use meds_config,           only : meds_config_t
   use meds_time,             only : meds_time_t, time_advance_days
   use meds_core_state_types, only : site_t, site_alloc, site_free
   use meds_output_config,    only : FC_DAY, FC_MONTH, FC_YEAR, FC_RUN, SYNC_FLUSH
   use meds_output_types,     only : output_manager_t
   use meds_output_registry,  only : manager_alloc
   use meds_output_integrate, only : output_integrate
   use meds_output_manager,   only : output_serialize_pending, output_manager_close
   use meds_netcdf_c
   use meds_test_support,     only : check, check_close, banner
   use test_ro_support,       only : set_site_agb, set_site_gpp, set_site_soil_temp
   implicit none

   type(meds_config_t)    :: cfg
   type(site_t)           :: site
   type(output_manager_t) :: mgr
   type(meds_time_t)      :: now
   real(wp)               :: dt

   call banner('output_roundtrip')

   !----- A fixed 1-cohort / 1-patch site (constant slot set). -----!
   call site_alloc(site, 1_ik, 16_ik, 4_ik, 4_ik)
   site%cohort%n = 1_ik ; site%patch%n = 1_ik
   site%cohort%nplant(1) = 1.0_wp ; site%cohort%dbh(1) = 5.0_wp
   site%cohort%height(1) = 4.0_wp ; site%cohort%basal_area(1) = 19.6_wp
   site%cohort%leaf_area(1) = 2.0_wp ; site%cohort%growth_avg(1) = 0.1_wp
   site%cohort%pft(1) = 1_ik ; site%cohort%owner_patch(1) = 1_ik ; site%cohort%global_id(1) = 1_ik
   site%patch%area(1) = 1.0_wp ; site%patch%cohort_offset(1) = 1_ik ; site%patch%cohort_count(1) = 1_ik
   site%patch%age(1) = 0.0_wp ; site%patch%dist_type(1) = 1_ik ; site%patch%global_id(1) = 1_ik

   !----- Constant carbon + soil sources (fast loop is off in this unit test, so set them here    !
   !      via opaque setters, as the fast loop would). gpp_accum=2 -> daily AGG_SUM gpp_site=2;     !
   !      soil_temp=290 -> daily AGG_TMEAN soil_temp_site=290 per layer.                            !
   call set_site_gpp(site, 2.0_wp)
   call set_site_soil_temp(site, 290.0_wp)

   !----- Output config: daily + annual streams; enable the energy group (soil). -----!
   cfg = build_cfg()
   dt  = 86400.0_wp
   call manager_alloc(mgr, cfg)

   !----- Walk 3 days of 2000-01; agb 10 / 20 / 30. is_new_day closes the previous day. -----!
   now = meds_time_t(year=2000_ik, month=1_ik, day=1_ik)
   call set_site_agb(site, 10.0_wp)
   call output_integrate(mgr, site, now, dt, .false., .false., .false.)
   call output_serialize_pending(mgr)

   now = time_advance_days(now, 1_ik)                  ! 2000-01-02
   call set_site_agb(site, 20.0_wp)
   call output_integrate(mgr, site, now, dt, .true., .false., .false.)
   call output_serialize_pending(mgr)

   now = time_advance_days(now, 1_ik)                  ! 2000-01-03
   call set_site_agb(site, 30.0_wp)
   call output_integrate(mgr, site, now, dt, .true., .false., .false.)
   call output_serialize_pending(mgr)

   call output_manager_close(mgr, .true.)              ! flush the day-3 daily + the 2000 annual partial

   !----- Re-open and assert. -----!
   call check_daily('test_ro-D-200001.nc')
   call check_annual('test_ro-Y.nc')

   call site_free(site)
   write(*,'(a)') 'test_output_roundtrip: ALL PASSED'

contains

   function build_cfg() result(c)
      type(meds_config_t) :: c
      c%output%enabled    = .true.
      c%output%dir        = '.'
      c%output%prefix     = 'test_ro'
      c%output%cohort_max = 16_ik
      c%output%patch_max  = 4_ik
      c%output%sync_every = SYNC_FLUSH
      c%output%freq_on    = [.false., .true., .false., .true.]     ! daily + annual
      c%output%file_chunk = [FC_DAY, FC_MONTH, FC_YEAR, FC_RUN]
      c%output%grp_on     = [.true., .true., .false., .true., .false.]      ! structure + carbon + ENERGY (soil)
   end function build_cfg

   !----- Daily file: 3 records, agb_site = [10,20,30], n_cohort=1, agb_cohort slab = [10,20,30]. !
   subroutine check_daily(path)
      character(len=*), intent(in) :: path
      integer(c_int)    :: ncid, vid
      integer(c_size_t) :: nt
      real(c_double)    :: agbs(3), agbc(3), gpp(3), soilt(3)
      integer(c_int)    :: ncoh(3)
      call nc_check(nc_open_f(trim(path), NC_NOWRITE, ncid), 'open daily')
      call nc_check(nc_inq_dimlen_f(ncid, 'time', nt), 'daily time len')
      call check(int(nt, ik) == 3_ik, 'daily has 3 records')
      call nc_check(nc_inq_varid_f(ncid, 'agb_site', vid), 'daily agb_site id')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t], [3_c_size_t], agbs), 'get agb_site')
      call check_close(real(agbs(1), wp), 10.0_wp, 1.0e-9_wp, 'daily agb_site day1')
      call check_close(real(agbs(2), wp), 20.0_wp, 1.0e-9_wp, 'daily agb_site day2')
      call check_close(real(agbs(3), wp), 30.0_wp, 1.0e-9_wp, 'daily agb_site day3')
      call nc_check(nc_inq_varid_f(ncid, 'n_cohort', vid), 'daily n_cohort id')
      call nc_check(nc_get_vara_int(ncid, vid, [0_c_size_t], [3_c_size_t], ncoh), 'get n_cohort')
      call check(all(ncoh == 1_c_int), 'daily n_cohort == 1')
      call nc_check(nc_inq_varid_f(ncid, 'agb_cohort', vid), 'daily agb_cohort id')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t], [3_c_size_t, 1_c_size_t], agbc), 'get agb_cohort')
      call check_close(real(agbc(1), wp), 10.0_wp, 1.0e-9_wp, 'daily agb_cohort day1 slot1')
      call check_close(real(agbc(3), wp), 30.0_wp, 1.0e-9_wp, 'daily agb_cohort day3 slot1')
      !----- P1 carbon: gpp_site AGG_SUM over 1 step/day = 2.0 each day. -----!
      call nc_check(nc_inq_varid_f(ncid, 'gpp_site', vid), 'daily gpp_site id')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t], [3_c_size_t], gpp), 'get gpp_site')
      call check_close(real(gpp(2), wp), 2.0_wp, 1.0e-9_wp, 'daily gpp_site (AGG_SUM) = 2')
      !----- P1 soil (DIM_SOIL): soil_temp_site(time, soil) AGG_TMEAN = 290 K per layer. -----!
      call nc_check(nc_inq_varid_f(ncid, 'soil_temp_site', vid), 'daily soil_temp_site id')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t], [3_c_size_t, 1_c_size_t], soilt), 'get soil_temp')
      call check_close(real(soilt(1), wp), 290.0_wp, 1.0e-9_wp, 'daily soil_temp_site layer1 = 290')
      call nc_check(nc_close(ncid), 'close daily')
   end subroutine check_daily

   !----- Annual file: 1 record, site-level ONLY (no cohort dim), agb_site = mean(10,20,30)=20. --!
   subroutine check_annual(path)
      character(len=*), intent(in) :: path
      integer(c_int)    :: ncid, vid, st
      integer(c_size_t) :: nt
      real(c_double)    :: agbs(1)
      call nc_check(nc_open_f(trim(path), NC_NOWRITE, ncid), 'open annual')
      call nc_check(nc_inq_dimlen_f(ncid, 'time', nt), 'annual time len')
      call check(int(nt, ik) == 1_ik, 'annual has 1 record')
      call nc_check(nc_inq_varid_f(ncid, 'agb_site', vid), 'annual agb_site id')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t], [1_c_size_t], agbs), 'get annual agb_site')
      call check_close(real(agbs(1), wp), 20.0_wp, 1.0e-9_wp, 'annual agb_site = mean(10,20,30)')
      !----- annual is site-level only: agb_cohort must NOT exist. -----!
      st = nc_inq_varid_f(ncid, 'agb_cohort', vid)
      call check(st /= NC_NOERR, 'annual has NO agb_cohort (site-level only)')
      call nc_check(nc_close(ncid), 'close annual')
   end subroutine check_annual

end program test_output_roundtrip
