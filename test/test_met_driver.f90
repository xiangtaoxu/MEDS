!==========================================================================================!
! test_met_driver -- unit tests for the P0 meteorological-forcing library (meds_forcing):       !
! the pure disaggregation kernels, the CONST backend, and a NetCDF round-trip (write a small     !
! (time,grid) MEDS forcing file, read it back through the reader). Mirrors the design test plan   !
! (MEDS_FORCING_DESIGN.md section 9). House style: check/check_true/nfail/error stop 1.            !
!==========================================================================================!
program test_met_driver
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : t_3ple
   use meds_time,            only : meds_time_t, seconds_between, seconds_into_day,             &
                                    time_advance_seconds
   use meds_therm_lib,          only : sat_vapor_pressure
   use meds_forcing_config,  only : forcing_config_t, MET_BACKEND_CONST, MET_BACKEND_NETCDF,    &
                                    SWPART_CLEARIDX, SWPART_WEISS_NORMAN, INTERP_LINEAR,        &
                                    INTERP_STEP, METAVG_END, CLAMP_HOLD, CLAMP_ERROR,           &
                                    GRIDMATCH_EXPLICIT, GRIDMATCH_NEAREST
   use meds_forcing_types,   only : met_forcing_t, met_driver_t
   use meds_forcing_kernels, only : interpolate_forcing, dewpoint_to_specific_humidity,         &
                                    rh_to_specific_humidity, precip_phase, partition_shortwave, &
                                    met_solar_cosz, cosz_reconstruct_factor, disaggregate_sw,   &
                                    great_circle_distance, nearest_grid_index, wind_log_profile, &
                                    lapse_air_temperature, lapse_pressure
   use meds_met_driver,      only : met_open, met_advance, met_instant, met_close,             &
                                   MET_OK, MET_ERR_WINDOW_NOT_WHOLE_YEARS,                     &
                                   MET_ERR_START_NOT_A_RECORD, MET_ERR_WINDOW_NOT_COVERED
   use meds_netcdf_c
   use iso_c_binding,        only : c_int, c_size_t, c_double
   implicit none
   integer(ik) :: nfail
   character(len=*), parameter :: NCFILE = 'test_met_driver_tmp.nc'
   nfail = 0_ik

   call test_interpolation()
   call test_humidity()
   call test_precip_phase()
   call test_sw_partition()
   call test_sw_partition_weissnorman()
   call test_cosz_reconstruction()
   call test_const_backend()
   call test_netcdf_roundtrip()
   call test_nearest_grid()
   call test_wind_lapse()
   call test_multiyear_cycling()
   call test_recycle_anchor_phase()

   if (nfail == 0_ik) then
      print '(a)', 'test_met_driver: ALL PASSED'
   else
      print '(a,i0,a)', 'test_met_driver: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   subroutine check(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) <= atol) then
         print '(a,a,a,es13.5,a,es13.5)', '  ok   : ', name, '  (', got, ' ~ ', expect, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,es13.5,a,es13.5)', '  FAIL : ', name, '  got ', got, ' expected ', expect
      end if
   end subroutine check

   subroutine check_true(name, cond, val)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      real(wp),         intent(in) :: val
      if (cond) then
         print '(a,a,a,es13.5,a)', '  ok   : ', name, '  (', val, ')'
      else
         nfail = nfail + 1_ik ; print '(a,a,a,es13.5,a)', '  FAIL : ', name, '  (', val, ')'
      end if
   end subroutine check_true

   !----- 1. temporal interpolation policies. ---------------------------------------------!
   subroutine test_interpolation()
      print '(a)', '-- test 1: interpolation --'
      call check('linear midpoint', interpolate_forcing(INTERP_LINEAR, 10.0_wp, 20.0_wp, 0.5_wp), 15.0_wp, 1.0e-12_wp)
      call check('linear at prev',  interpolate_forcing(INTERP_LINEAR, 10.0_wp, 20.0_wp, 0.0_wp), 10.0_wp, 1.0e-12_wp)
      call check('step holds prev', interpolate_forcing(INTERP_STEP,   10.0_wp, 20.0_wp, 0.9_wp), 10.0_wp, 1.0e-12_wp)
   end subroutine test_interpolation

   !----- 2. humidity from dewpoint (matches independent Bolton). --------------------------!
   subroutine test_humidity()
      real(wp) :: q, e, td, p
      print '(a)', '-- test 2: humidity --'
      td = 285.0_wp ; p = 98000.0_wp
      q  = dewpoint_to_specific_humidity(td, p)
      e  = sat_vapor_pressure(td)
      call check('q from dewpoint == 0.622 e/(p-0.378e)', q, 0.622_wp*e/(p-0.378_wp*e), 1.0e-14_wp)
      call check_true('q in physical range', q > 0.001_wp .and. q < 0.02_wp, q)
      ! RH=100% at T=Td gives the same q as dewpoint=T
      call check('rh=1 at T=Td matches dewpoint q', rh_to_specific_humidity(1.0_wp, td, p), q, 1.0e-14_wp)
   end subroutine test_humidity

   !----- 3. precip phase split (mass-conserving). -----------------------------------------!
   subroutine test_precip_phase()
      real(wp) :: rain, snow, tot
      print '(a)', '-- test 3: precip phase --'
      tot = 1.0e-4_wp
      call precip_phase(tot, t_3ple + 5.0_wp, rain, snow)
      call check('warm -> all rain', rain, tot, 1.0e-16_wp)
      call check('warm -> no snow',  snow, 0.0_wp, 1.0e-16_wp)
      call precip_phase(tot, t_3ple - 5.0_wp, rain, snow)
      call check('cold -> all snow', snow, tot, 1.0e-16_wp)
      call precip_phase(tot, t_3ple, rain, snow)
      call check('at t_3ple mass conserved', rain + snow, tot, 1.0e-16_wp)
      call check_true('partial snow at t_3ple', snow > 0.0_wp .and. snow < tot, snow)
   end subroutine test_precip_phase

   !----- 4. shortwave partition (Erbs clearness-index). -----------------------------------!
   subroutine test_sw_partition()
      real(wp) :: pb, pd, nb, nd, cosz
      real(wp), parameter :: PS = 98000.0_wp
      print '(a)', '-- test 4: SW partition (Erbs) --'
      cosz = 0.8_wp
      call partition_shortwave(800.0_wp, cosz, PS, SWPART_CLEARIDX, pb, pd, nb, nd)
      call check('4 streams sum to total SW', pb+pd+nb+nd, 800.0_wp, 1.0e-9_wp)
      call check_true('all streams >= 0', pb>=0 .and. pd>=0 .and. nb>=0 .and. nd>=0, min(pb,pd,nb,nd))
      call check_true('clear sky -> mostly beam', (pb+nb) > (pd+nd), (pb+nb)-(pd+nd))
      ! overcast (low kt) -> mostly diffuse
      call partition_shortwave(150.0_wp, cosz, PS, SWPART_CLEARIDX, pb, pd, nb, nd)
      call check_true('overcast -> mostly diffuse', (pd+nd) > (pb+nb), (pd+nd)-(pb+nb))
      ! night -> all zero
      call partition_shortwave(0.0_wp, 0.0_wp, PS, SWPART_CLEARIDX, pb, pd, nb, nd)
      call check('night -> zero SW', pb+pd+nb+nd, 0.0_wp, 1.0e-30_wp)
   end subroutine test_sw_partition

   !----- Weiss-Norman band-specific partition: exact energy conservation, clear/overcast beam    !
   !      trend, dawn all-diffuse, and pressure sensitivity of the beam fraction.                   !
   subroutine test_sw_partition_weissnorman()
      real(wp) :: pb, pd, nb, nd, cosz, beam_hi, beam_lo
      print '(a)', '-- test 4b: SW partition (Weiss-Norman) --'
      cosz = 0.8_wp
      call partition_shortwave(800.0_wp, cosz, 98000.0_wp, SWPART_WEISS_NORMAN, pb, pd, nb, nd)
      call check('WN: 4 streams sum to total SW', pb+pd+nb+nd, 800.0_wp, 1.0e-8_wp)
      call check_true('WN: all streams >= 0', pb>=0 .and. pd>=0 .and. nb>=0 .and. nd>=0, min(pb,pd,nb,nd))
      call check_true('WN: clear sky -> mostly beam', (pb+nb) > (pd+nd), (pb+nb)-(pd+nd))
      call check_true('WN: PAR beam > NIR-diffuse in clear sky', pb > nd, pb-nd)
      beam_hi = pb + nb
      ! overcast (low observed vs potential) -> beam collapses, diffuse dominates
      call partition_shortwave(120.0_wp, cosz, 98000.0_wp, SWPART_WEISS_NORMAN, pb, pd, nb, nd)
      call check('WN overcast: streams sum to total', pb+pd+nb+nd, 120.0_wp, 1.0e-8_wp)
      call check_true('WN overcast -> mostly diffuse', (pd+nd) > (pb+nb), (pd+nd)-(pb+nb))
      ! grazing sun (cosz below WN_COSZ_MIN) -> all diffuse, still conserving
      call partition_shortwave(50.0_wp, 0.01_wp, 98000.0_wp, SWPART_WEISS_NORMAN, pb, pd, nb, nd)
      call check('WN dawn: streams sum to total', pb+pd+nb+nd, 50.0_wp, 1.0e-9_wp)
      call check('WN dawn: no beam', pb+nb, 0.0_wp, 1.0e-30_wp)
      ! pressure sensitivity: psurf_pa enters the WN optical depth, so the split MUST change with it
      call partition_shortwave(800.0_wp, cosz, 101325.0_wp, SWPART_WEISS_NORMAN, pb, pd, nb, nd)
      beam_lo = pb + nb
      call check_true('WN: surface pressure changes the beam split', abs(beam_hi-beam_lo) > 1.0e-6_wp, &
                      abs(beam_hi-beam_lo))
   end subroutine test_sw_partition_weissnorman

   !----- 5. cosz reconstruction conserves the interval mean (the load-bearing method). -----!
   subroutine test_cosz_reconstruction()
      type(meds_time_t) :: t
      real(wp) :: factor, f_avg, win_start, dt_win, dt_sub, sec, cosz, fsum, secz_sum, secz, mean_secz
      integer(ik) :: i, n
      real(wp), parameter :: LAT = 42.44_wp, LON = -76.50_wp
      print '(a)', '-- test 5: cosz reconstruction (interval-mean conserving) --'
      t = meds_time_t(year=2020_ik, month=7_ik, day=1_ik)
      ! A SUNRISE window (10-11 UTC ~ 05-06 local at Ithaca), where cosz rises steeply from ~0.07 --
      ! this is where ED2's <sec z> mistake is worst (design §9). Solar noon is ~17 UTC here.
      f_avg = 500.0_wp ; win_start = 10.0_wp*3600.0_wp ; dt_win = 3600.0_wp
      n = 12_ik ; dt_sub = dt_win / real(n, wp)
      factor = cosz_reconstruct_factor(t, win_start, dt_sub, dt_win, LAT, LON, 0.0_wp, .true.)  ! = 1/<cosz>
      ! (a) the CORRECT 1/<cosz> factor conserves the interval mean EXACTLY over the same subsamples:
      fsum = 0.0_wp ; secz_sum = 0.0_wp
      do i = 1_ik, n
         sec  = win_start + (real(i,wp)-0.5_wp)*dt_sub
         cosz = met_solar_cosz(t, sec, LAT, LON, 0.0_wp, .true.)
         fsum = fsum + disaggregate_sw(f_avg, cosz, factor)
         secz_sum = secz_sum + 1.0_wp/max(cosz, 0.03_wp)          ! <sec z> = <1/cosz> (ED2 mean_daysecz, clamped)
      end do
      call check('1/<cosz> conserves interval mean (sunrise window)', fsum/real(n,wp), f_avg, 1.0e-8_wp)
      ! (b) the GENUINE ED2 error is F = F_avg*cosz*<sec z> (multiply by the MEAN SECANT, ed_met_driver
      !     fperp=nbdsf*secz; flux=fperp*cosz). By Cauchy-Schwarz <cosz><sec z> >= 1, so it is BIASED HIGH:
      fsum = 0.0_wp
      do i = 1_ik, n
         sec  = win_start + (real(i,wp)-0.5_wp)*dt_sub
         cosz = met_solar_cosz(t, sec, LAT, LON, 0.0_wp, .true.)
         fsum = fsum + f_avg*cosz*(secz_sum/real(n,wp))
      end do
      mean_secz = fsum/real(n,wp)
      call check_true('<sec z> form is BIASED HIGH (does NOT conserve)', mean_secz > f_avg + 1.0_wp, mean_secz - f_avg)
      ! (c) a fully-night window -> factor 0 -> all SW routes to 0
      factor = cosz_reconstruct_factor(t, 4.0_wp*3600.0_wp, dt_sub, dt_win, LAT, LON, 0.0_wp, .true.) ! 04-05 UTC = night
      call check('night window -> factor 0', factor, 0.0_wp, 1.0e-30_wp)
   end subroutine test_cosz_reconstruction

   !----- 6. CONST backend returns the reference climate (SW = 400, held flat), cosz derived. -!
   subroutine test_const_backend()
      type(met_driver_t)     :: drv
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met
      type(meds_time_t)      :: now
      print '(a)', '-- test 6: CONST backend --'
      fc%backend = MET_BACKEND_CONST ; fc%latitude_deg = 42.44_wp ; fc%longitude_deg = -76.50_wp
      call met_open(drv, fc)
      now = meds_time_t(year=2020_ik, month=7_ik, day=1_ik, hour=17_ik)   ! ~solar noon at Ithaca
      call met_advance(drv, now)
      met = met_instant(drv, now)
      call check('CONST swdown = 400', met%swdown(), 400.0_wp, 1.0e-9_wp)
      call check_true('CONST cosz recomputed > 0 at noon', met%cosz > 0.5_wp, met%cosz)
      call check_true('CONST rho_air ~ 1.2', abs(met%rho_air - 1.2_wp) < 0.2_wp, met%rho_air)
      call met_close(drv)
   end subroutine test_const_backend

   !----- 7. NetCDF round-trip: write a (time=25, grid=2) file, read it back. ----------------!
   subroutine test_netcdf_roundtrip()
      type(met_driver_t)     :: drv
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met_noon, met_night, met_g2
      type(meds_time_t)      :: base, now
      real(wp) :: sw_g1
      print '(a)', '-- test 7: NetCDF round-trip (multi-grid) --'
      base = meds_time_t(year=2020_ik, month=7_ik, day=1_ik)
      call write_synthetic_forcing(NCFILE, base)

      fc%backend = MET_BACKEND_NETCDF ; fc%path = NCFILE ; fc%grid_index = 1_ik
      fc%dt_forcing = 3600.0_wp ; fc%avg_convention = METAVG_END ; fc%sw_partition = SWPART_CLEARIDX
      fc%latitude_deg = 42.44_wp ; fc%longitude_deg = -76.50_wp ; fc%utc_offset_h = 0.0_wp
      fc%apply_solar_longitude = .true. ; fc%recycle = .false.
      call met_open(drv, fc)
      call check_true('opened file: nrec=25', drv%nrec == 25_ik, real(drv%nrec, wp))
      call check_true('opened file: ngrid=2', drv%ngrid == 2_ik, real(drv%ngrid, wp))

      ! noon (17 UTC): SW positive, Tair follows the diurnal input
      now = time_advance_seconds(base, 17.0_wp*3600.0_wp)
      call met_advance(drv, now) ; met_noon = met_instant(drv, now)
      call check_true('noon SWdown > 0', met_noon%swdown() > 100.0_wp, met_noon%swdown())
      call check_true('noon cosz > 0.5', met_noon%cosz > 0.5_wp, met_noon%cosz)
      call check_true('Tair in [285,305]', met_noon%tair_k > 285.0_wp .and. met_noon%tair_k < 305.0_wp, met_noon%tair_k)
      sw_g1 = met_noon%swdown()

      ! night (04 UTC): SW exactly 0
      now = time_advance_seconds(base, 4.0_wp*3600.0_wp)
      call met_advance(drv, now) ; met_night = met_instant(drv, now)
      call check('night SWdown = 0', met_night%swdown(), 0.0_wp, 1.0e-12_wp)
      call met_close(drv)

      ! multi-grid: grid_index=2 carries a distinct (scaled) SW series
      fc%grid_index = 2_ik
      call met_open(drv, fc)
      now = time_advance_seconds(base, 17.0_wp*3600.0_wp)
      call met_advance(drv, now) ; met_g2 = met_instant(drv, now)
      call check_true('grid 2 SW differs from grid 1', abs(met_g2%swdown() - sw_g1) > 1.0_wp, met_g2%swdown()-sw_g1)
      call met_close(drv)

      call test_edge_paths(base)
   end subroutine test_netcdf_roundtrip

   !----- 8. Edge paths the review flagged (regression tests for the CLAMP_HOLD + recycle fixes). !
   subroutine test_edge_paths(base)
      type(meds_time_t), intent(in) :: base
      type(met_driver_t)     :: drv
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met
      type(meds_time_t)      :: now
      integer(ik) :: st
      real(wp) :: tair0, tair1, expect_hold
      print '(a)', '-- test 8: edge paths (start-hold + recycle) --'
      fc%backend = MET_BACKEND_NETCDF ; fc%path = NCFILE ; fc%grid_index = 1_ik
      fc%dt_forcing = 3600.0_wp ; fc%avg_convention = METAVG_END ; fc%sw_partition = SWPART_CLEARIDX
      fc%latitude_deg = 42.44_wp ; fc%longitude_deg = -76.50_wp ; fc%apply_solar_longitude = .true.
      ! synthetic Tair(hh) = 288 + 8 sin(2 pi (hh-15)/24)  -- reproduce the file's values for the checks:
      tair0 = 288.0_wp + 8.0_wp*sin(2.0_wp*3.14159265_wp*(0.0_wp-15.0_wp)/24.0_wp)
      tair1 = 288.0_wp + 8.0_wp*sin(2.0_wp*3.14159265_wp*(1.0_wp-15.0_wp)/24.0_wp)

      ! (a) CLAMP_HOLD: start before the first record holds rec1; then MARCHING into [t0,t1) must
      !     interpolate (the bug held rec1 for the whole first interval).
      fc%start_clamp = CLAMP_HOLD ; fc%recycle = .false.
      call met_open(drv, fc)
      now = time_advance_seconds(base, -12.0_wp*3600.0_wp)          ! 12 h before the first record
      call met_advance(drv, now) ; met = met_instant(drv, now)
      call check('hold before t0 -> record 1', met%tair_k, tair0, 1.0e-4_wp)
      now = time_advance_seconds(base, 0.5_wp*3600.0_wp)            ! 00:30 -> midpoint of [t0,t1)
      call met_advance(drv, now) ; met = met_instant(drv, now)
      expect_hold = 0.5_wp*(tair0 + tair1)
      call check('march into [t0,t1) interpolates (not stuck at rec1)', met%tair_k, expect_hold, 1.0e-4_wp)
      call met_close(drv)

      ! (b) A 25-record (24 h span) file CANNOT be recycled: a sub-year window drifts both
      !     hour-of-day and day-of-year on every wrap. This used to fall through to an
      !     absolute-seconds span-wrap SILENTLY, which is exactly how a 29-yr production run came
      !     to be driven by ~10 h-anti-phased sub-daily shortwave. It must now be REJECTED.
      fc%start_clamp = CLAMP_ERROR ; fc%recycle = .true.
      fc%recycle_start = base
      fc%recycle_end   = time_advance_seconds(base, 24.0_wp*3600.0_wp)   ! 1 day, not a whole year
      call met_open(drv, fc, stat=st)
      call check_true('sub-year recycle window rejected (not silently span-wrapped)',            &
                      st == MET_ERR_WINDOW_NOT_WHOLE_YEARS, real(st, wp))

      ! (c) A whole-year window whose start does NOT match a record stamp is rejected too --
      !     MEDS never guesses the start of the day (e.g. a config saying 00:00 against an
      !     end-of-interval ERA5-Land file whose records are stamped 01:00).
      fc%recycle_start = time_advance_seconds(base, 1800.0_wp)           ! 00:30, between records
      fc%recycle_end   = meds_time_t(2021_ik, 7_ik, 1_ik, 0_ik, 30_ik)
      call met_open(drv, fc, stat=st)
      call check_true('recycle_start off the record grid rejected',                              &
                      st == MET_ERR_START_NOT_A_RECORD, real(st, wp))

      ! (d) A whole-year window ON the record grid, but the file only holds 24 h of it.
      fc%recycle_start = base
      fc%recycle_end   = meds_time_t(2021_ik, 7_ik, 1_ik)
      call met_open(drv, fc, stat=st)
      call check_true('window not covered by the file rejected',                                 &
                      st == MET_ERR_WINDOW_NOT_COVERED, real(st, wp))
      call met_close(drv)
   end subroutine test_edge_paths

   !----- The recycle mapping must preserve hour-of-day and day-of-year EXACTLY, for an anchor    !
   !      anywhere in the calendar -- not only Jan-1. This is the regression for the met-recycle    !
   !      phase bug: the old classifier accepted only Jan-1 00:00 files and silently span-wrapped    !
   !      everything else, which shifts hour-of-day whenever the span is not a whole number of days. !
   subroutine test_recycle_anchor_phase()
      character(len=*), parameter :: YF = 'test_met_anchorfile_tmp.nc'
      type(met_driver_t)     :: drv
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: m_ref, m_map
      integer(ik) :: st
      print '(a)', '-- test: recycle anchor + sub-daily phase --'
      call write_yearfile(YF, 2021_ik)                       ! 365 daily records from 2021-01-01 00:00
      fc%backend = MET_BACKEND_NETCDF ; fc%path = YF ; fc%grid_index = 1_ik
      fc%dt_forcing = 86400.0_wp ; fc%avg_convention = METAVG_END ; fc%sw_partition = SWPART_CLEARIDX
      fc%recycle = .true. ; fc%apply_solar_longitude = .false.

      !----- A MID-YEAR window: 2021-03-01 .. 2022-03-01. Plain year substitution would map a
      !      January instant to 2021-01, which is OUTSIDE this window; the anchored form sends it
      !      to 2022-01 instead. The file only reaches 2021-12-31, so this window is (correctly)
      !      rejected as not covered -- which is itself the check that the coverage test bites.
      fc%recycle_start = meds_time_t(2021_ik, 3_ik, 1_ik)
      fc%recycle_end   = meds_time_t(2022_ik, 3_ik, 1_ik)
      call met_open(drv, fc, stat=st)
      call check_true('mid-year window beyond the file is rejected',                             &
                      st == MET_ERR_WINDOW_NOT_COVERED, real(st, wp))

      !----- Whole-year window on the file: every model year must read the SAME calendar day, for
      !      model years both before and after the window (the modulo must not go negative).
      fc%recycle_start = meds_time_t(2021_ik, 1_ik, 1_ik)
      fc%recycle_end   = meds_time_t(2022_ik, 1_ik, 1_ik)
      call met_open(drv, fc, stat=st)
      call check_true('whole-year window accepted', st == MET_OK, real(st, wp))
      call met_advance(drv, meds_time_t(2021_ik,9_ik,17_ik)) ; m_ref = met_instant(drv, meds_time_t(2021_ik,9_ik,17_ik))
      call met_advance(drv, meds_time_t(2049_ik,9_ik,17_ik)) ; m_map = met_instant(drv, meds_time_t(2049_ik,9_ik,17_ik))
      call check('29 yr later reads the same calendar day', m_map%tair_k, m_ref%tair_k, 1.0e-9_wp)
      call met_advance(drv, meds_time_t(2013_ik,9_ik,17_ik)) ; m_map = met_instant(drv, meds_time_t(2013_ik,9_ik,17_ik))
      call check('a model year BEFORE the window maps in too', m_map%tair_k, m_ref%tair_k, 1.0e-9_wp)
      call met_close(drv)
   end subroutine test_recycle_anchor_phase

   !----- Write a synthetic (time=25 hourly, grid=2) MEDS forcing NetCDF via meds_netcdf_c. ----!
   subroutine write_synthetic_forcing(path, base)
      character(len=*),  intent(in) :: path
      type(meds_time_t), intent(in) :: base
      integer, parameter :: NT = 25, NG = 2
      integer(c_int) :: st, ncid, td, gd, vt, vla, vlo, vv(8)
      integer(c_int) :: dims2(2), dims1(1)
      real(c_double) :: tsec(NT), lat(NG), lon(NG)
      real(c_double) :: dat(NG, NT)        ! (grid, time) column-major == C [time][grid]
      integer :: it, ig, k
      real(wp) :: hh, sw_mean
      character(len=8), parameter :: vnames(8) = ['Tair    ','Qair    ','PSurf   ','Wind    ', &
                                                  'Rainf   ','SWdown  ','LWdown  ','CO2air  ']
      integer(c_size_t) :: start2(2), count2(2), start1(1), count1(1)

      st = nc_create_f(path, NC_NETCDF4, ncid) ; call nc_check(st, 'write: create')
      st = nc_def_dim_f(ncid, 'time', int(NT, c_size_t), td) ; call nc_check(st, 'write: time dim')
      st = nc_def_dim_f(ncid, 'grid', int(NG, c_size_t), gd) ; call nc_check(st, 'write: grid dim')
      dims1(1) = td
      st = nc_def_var_f(ncid, 'time', NC_DOUBLE, 1, dims1, vt) ; call nc_check(st, 'write: time var')
      st = nc_put_att_text_f(ncid, vt, 'units', int(len_trim('seconds since 2020-07-01 00:00:00'), c_size_t), &
                             'seconds since 2020-07-01 00:00:00') ; call nc_check(st, 'write: units')
      dims1(1) = gd
      st = nc_def_var_f(ncid, 'latitude',  NC_DOUBLE, 1, dims1, vla) ; call nc_check(st, 'write: lat var')
      st = nc_def_var_f(ncid, 'longitude', NC_DOUBLE, 1, dims1, vlo) ; call nc_check(st, 'write: lon var')
      dims2(1) = td ; dims2(2) = gd        ! [time, grid] (slowest first)
      do k = 1, 8
         st = nc_def_var_f(ncid, trim(vnames(k)), NC_DOUBLE, 2, dims2, vv(k))
         call nc_check(st, 'write: var '//trim(vnames(k)))
      end do
      st = nc_enddef(ncid) ; call nc_check(st, 'write: enddef')

      do it = 1, NT
         tsec(it) = real(it - 1, c_double) * 3600.0_c_double
      end do
      lat = [42.5_c_double, 42.4_c_double] ; lon = [-76.6_c_double, -76.5_c_double]
      start1(1) = 0_c_size_t ; count1(1) = int(NT, c_size_t)
      st = nc_put_vara_double(ncid, vt, start1, count1, tsec) ; call nc_check(st, 'write: time vals')
      count1(1) = int(NG, c_size_t)
      st = nc_put_vara_double(ncid, vla, start1, count1, lat) ; call nc_check(st, 'write: lat vals')
      st = nc_put_vara_double(ncid, vlo, start1, count1, lon) ; call nc_check(st, 'write: lon vals')

      start2 = [0_c_size_t, 0_c_size_t] ; count2 = [int(NT, c_size_t), int(NG, c_size_t)]
      do k = 1, 8
         do it = 1, NT
            hh = real(it - 1, wp)                      ! UTC hour 0..24
            do ig = 1, NG
               select case (k)
               case (1) ; dat(ig,it) = 288.0_wp + 8.0_wp*sin(2.0_wp*3.14159265_wp*(hh-15.0_wp)/24.0_wp) ! Tair
               case (2) ; dat(ig,it) = 0.008_wp        ! Qair
               case (3) ; dat(ig,it) = 98000.0_wp      ! PSurf
               case (4) ; dat(ig,it) = 3.0_wp          ! Wind
               case (5) ; dat(ig,it) = 0.0_wp          ! Rainf
               case (6)                                 ! SWdown: daytime bump, scaled per grid cell
                  sw_mean = max(0.0_wp, 900.0_wp*sin(3.14159265_wp*(hh-11.0_wp)/12.0_wp))
                  if (hh < 11.0_wp .or. hh > 23.0_wp) sw_mean = 0.0_wp
                  dat(ig,it) = sw_mean * (1.0_wp + 0.3_wp*real(ig-1, wp))
               case (7) ; dat(ig,it) = 320.0_wp        ! LWdown
               case (8) ; dat(ig,it) = 415.0_wp        ! CO2air
               end select
            end do
         end do
         st = nc_put_vara_double(ncid, vv(k), start2, count2, dat)
         call nc_check(st, 'write: vals '//trim(vnames(k)))
      end do
      st = nc_close(ncid) ; call nc_check(st, 'write: close')
   end subroutine write_synthetic_forcing

   !----- Nearest-grid match: pure kernel (argmin + great-circle) and the reader override. --------!
   subroutine test_nearest_grid()
      type(met_driver_t)     :: drv
      type(forcing_config_t) :: fc
      real(wp)    :: lon(3), lat(3), d
      integer(ik) :: idx
      print '(a)', '-- test: nearest-grid match --'
      lon = [-80.0_wp, -76.5_wp, -70.0_wp] ; lat = [40.0_wp, 42.44_wp, 45.0_wp]
      call check('argmin picks the coincident cell', real(nearest_grid_index(-76.5_wp,42.44_wp,lon,lat),wp), 2.0_wp, 0.5_wp)
      call check('argmin picks the SW cell',         real(nearest_grid_index(-79.0_wp,40.5_wp,lon,lat),wp), 1.0_wp, 0.5_wp)
      call check('great-circle self-distance = 0', great_circle_distance(-76.5_wp,42.44_wp,-76.5_wp,42.44_wp), 0.0_wp, 1.0e-6_wp)
      call check('great-circle 1 deg lat ~ 111 km', great_circle_distance(0.0_wp,0.0_wp,0.0_wp,1.0_wp), 111195.0_wp, 500.0_wp)
      !----- reader override: the 2-grid file has lon/lat cells (-76.6,42.5) and (-76.5,42.4). ---!
      call write_synthetic_forcing(NCFILE, meds_time_t(2020_ik,7_ik,1_ik))
      fc%backend = MET_BACKEND_NETCDF ; fc%path = NCFILE ; fc%dt_forcing = 3600.0_wp
      fc%grid_match = GRIDMATCH_NEAREST ; fc%grid_index = 1_ik      ! grid_index deliberately WRONG for cell 2
      fc%latitude_deg = 42.41_wp ; fc%longitude_deg = -76.49_wp     ! nearest to cell 2
      call met_open(drv, fc)
      call check('reader resolves nearest -> grid_index 2', real(drv%grid_index,wp), 2.0_wp, 0.5_wp)
      call met_close(drv)
      fc%latitude_deg = 42.51_wp ; fc%longitude_deg = -76.61_wp     ! nearest to cell 1
      call met_open(drv, fc)
      call check('reader resolves nearest -> grid_index 1', real(drv%grid_index,wp), 1.0_wp, 0.5_wp)
      call met_close(drv)
   end subroutine test_nearest_grid

   !----- Wind log-profile + elevation lapse (pure kernels). --------------------------------------!
   subroutine test_wind_lapse()
      real(wp) :: u, t, p
      real(wp), parameter :: DZ = 500.0_wp, G = 0.0065_wp
      print '(a)', '-- test: wind-height + elevation lapse --'
      u = wind_log_profile(3.0_wp, 10.0_wp, 40.0_wp, 0.1_wp)
      call check('wind log-profile value', u, 3.0_wp*log(40.0_wp/0.1_wp)/log(10.0_wp/0.1_wp), 1.0e-9_wp)
      call check_true('wind lifts to a higher reference height', u > 3.0_wp, u-3.0_wp)
      call check('degenerate z0 -> wind unchanged', wind_log_profile(3.0_wp,10.0_wp,40.0_wp,0.0_wp), 3.0_wp, 1.0e-30_wp)
      t = lapse_air_temperature(290.0_wp, DZ, G)
      call check('lapse cools a higher site', t, 290.0_wp - G*DZ, 1.0e-9_wp)
      call check('lapse dz=0 -> T unchanged', lapse_air_temperature(290.0_wp,0.0_wp,G), 290.0_wp, 1.0e-30_wp)
      p = lapse_pressure(101325.0_wp, 290.0_wp, DZ, G)
      call check_true('pressure drops at a higher site', p < 101325.0_wp, 101325.0_wp - p)
      call check('lapse dz=0 -> P unchanged', lapse_pressure(101325.0_wp,290.0_wp,0.0_wp,G), 101325.0_wp, 1.0e-6_wp)
      call check('isothermal (gamma=0) = barometric', lapse_pressure(101325.0_wp,290.0_wp,DZ,0.0_wp), &
                 101325.0_wp*exp(-9.80665_wp*DZ/(287.04_wp*290.0_wp)), 1.0e-2_wp)
   end subroutine test_wind_lapse

   !----- Multi-year CALENDAR recycling + Feb-29 reconciliation (whole-year daily file). ----------!
   subroutine test_multiyear_cycling()
      character(len=*), parameter :: YF = 'test_met_yearfile_tmp.nc'
      type(met_driver_t)     :: drv
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: m_ref, m_map
      print '(a)', '-- test: multi-year cycling + Feb-29 --'
      call write_yearfile(YF, 2021_ik)                             ! 365 daily records, Jan-1 2021 aligned
      fc%backend = MET_BACKEND_NETCDF ; fc%path = YF ; fc%grid_index = 1_ik
      fc%dt_forcing = 86400.0_wp ; fc%avg_convention = METAVG_END ; fc%sw_partition = SWPART_CLEARIDX
      fc%recycle = .true. ; fc%apply_solar_longitude = .false.
      !----- The cycle is DECLARED, never inferred from the file. --------------------------------!
      fc%recycle_start = meds_time_t(2021_ik, 1_ik, 1_ik)
      fc%recycle_end   = meds_time_t(2022_ik, 1_ik, 1_ik)
      call met_open(drv, fc)
      call check('n_cycle_years = 1', real(drv%n_cycle_years,wp), 1.0_wp, 0.5_wp)
      call check('cycle anchor year = 2021', real(drv%cycle_anchor%year,wp), 2021.0_wp, 0.5_wp)
      call check('cycle first record = #1', real(drv%irec_cycle_first,wp), 1.0_wp, 0.5_wp)
      !----- recycle identity: model 2023-07-01 reads the 2021-07-01 record (day-of-year exact). --!
      call met_advance(drv, meds_time_t(2021_ik,7_ik,1_ik)) ; m_ref = met_instant(drv, meds_time_t(2021_ik,7_ik,1_ik))
      call met_advance(drv, meds_time_t(2023_ik,7_ik,1_ik)) ; m_map = met_instant(drv, meds_time_t(2023_ik,7_ik,1_ik))
      call check('recycle: 2023-07-01 reads the 2021-07-01 record', m_map%tair_k, m_ref%tair_k, 1.0e-9_wp)
      !----- Feb-29 reconciliation: model 2024-02-29 (leap) maps to file 2021-02-28. --------------!
      call met_advance(drv, meds_time_t(2021_ik,2_ik,28_ik)) ; m_ref = met_instant(drv, meds_time_t(2021_ik,2_ik,28_ik))
      call met_advance(drv, meds_time_t(2024_ik,2_ik,29_ik)) ; m_map = met_instant(drv, meds_time_t(2024_ik,2_ik,29_ik))
      call check('Feb-29 (leap model) maps to file Feb-28', m_map%tair_k, m_ref%tair_k, 1.0e-9_wp)
      call met_close(drv)
   end subroutine test_multiyear_cycling

   !----- Write a whole-year DAILY (365 records, grid=1) forcing file; Tair encodes the day index. !
   subroutine write_yearfile(path, year)
      character(len=*), intent(in) :: path
      integer(ik),      intent(in) :: year
      integer, parameter :: NT = 365, NG = 1
      integer(c_int)    :: st, ncid, td, gd, vt, vla, vlo, vv(7), dims2(2), dims1(1)
      integer(c_size_t) :: start2(2), count2(2), start1(1), count1(1)
      real(c_double)    :: tsec(NT), la(NG), lo(NG), dat(NG, NT)
      integer :: it, k
      character(len=8), parameter :: vnames(7) = ['Tair    ','Qair    ','PSurf   ','Wind    ', &
                                                  'Rainf   ','SWdown  ','LWdown  ']
      character(len=40) :: units
      write(units,'(a,i0,a)') 'seconds since ', year, '-01-01 00:00:00'
      st = nc_create_f(path, NC_NETCDF4, ncid) ; call nc_check(st, 'yf: create')
      st = nc_def_dim_f(ncid, 'time', int(NT,c_size_t), td) ; call nc_check(st, 'yf: time dim')
      st = nc_def_dim_f(ncid, 'grid', int(NG,c_size_t), gd) ; call nc_check(st, 'yf: grid dim')
      dims1(1) = td ; st = nc_def_var_f(ncid, 'time', NC_DOUBLE, 1, dims1, vt) ; call nc_check(st, 'yf: time var')
      st = nc_put_att_text_f(ncid, vt, 'units', int(len_trim(units),c_size_t), trim(units)) ; call nc_check(st, 'yf: units')
      dims1(1) = gd
      st = nc_def_var_f(ncid, 'latitude',  NC_DOUBLE, 1, dims1, vla) ; call nc_check(st, 'yf: lat var')
      st = nc_def_var_f(ncid, 'longitude', NC_DOUBLE, 1, dims1, vlo) ; call nc_check(st, 'yf: lon var')
      dims2(1) = td ; dims2(2) = gd
      do k = 1, 7 ; st = nc_def_var_f(ncid, trim(vnames(k)), NC_DOUBLE, 2, dims2, vv(k)) ; call nc_check(st, 'yf: var') ; end do
      st = nc_enddef(ncid) ; call nc_check(st, 'yf: enddef')
      do it = 1, NT ; tsec(it) = real(it-1, c_double) * 86400.0_c_double ; end do
      la = 42.44_c_double ; lo = -76.50_c_double
      start1(1) = 0_c_size_t ; count1(1) = int(NT,c_size_t)
      st = nc_put_vara_double(ncid, vt, start1, count1, tsec) ; call nc_check(st, 'yf: time vals')
      count1(1) = int(NG,c_size_t)
      st = nc_put_vara_double(ncid, vla, start1, count1, la) ; call nc_check(st, 'yf: lat vals')
      st = nc_put_vara_double(ncid, vlo, start1, count1, lo) ; call nc_check(st, 'yf: lon vals')
      start2 = [0_c_size_t, 0_c_size_t] ; count2 = [int(NT,c_size_t), int(NG,c_size_t)]
      do k = 1, 7
         do it = 1, NT
            select case (k)
            case (1) ; dat(1,it) = 280.0_wp + real(it-1, wp)   ! Tair encodes the day index (unique per day)
            case (2) ; dat(1,it) = 0.006_wp                    ! Qair
            case (3) ; dat(1,it) = 99000.0_wp                  ! PSurf
            case (4) ; dat(1,it) = 2.0_wp                      ! Wind
            case (5) ; dat(1,it) = 0.0_wp                      ! Rainf
            case (6) ; dat(1,it) = 200.0_wp                    ! SWdown (daily mean)
            case (7) ; dat(1,it) = 300.0_wp                    ! LWdown
            end select
         end do
         st = nc_put_vara_double(ncid, vv(k), start2, count2, dat) ; call nc_check(st, 'yf: vals')
      end do
      st = nc_close(ncid) ; call nc_check(st, 'yf: close')
   end subroutine write_yearfile

end program test_met_driver
