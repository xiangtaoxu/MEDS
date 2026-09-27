! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_met_driver -- unit tests for the P0 meteorological-forcing library (meds_forcing):       !
! the pure disaggregation kernels, the CONST backend, and a NetCDF round-trip (write a small     !
! (time,grid) MEDS forcing file, read it back through the reader). Mirrors the design test plan   !
! (MEDS_FORCING_DESIGN.md section 9). House style: check/check_true/nfail/error stop 1.            !
! Argument 1 (from CMake): the shipped CO2 series, data/co2/, checked by test_co2_series.          !
!==========================================================================================!
program test_met_driver
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : t_3ple
   use meds_time,            only : meds_time_t, seconds_between, seconds_into_day,             &
                                    time_advance_seconds
   use meds_therm_lib,          only : sat_vapor_pressure, air_density
   use meds_forcing_config,  only : LW_CLEAR_BRUTSAERT, LW_CLEAR_IDSO
   use meds_forcing_config,  only : forcing_config_t, MET_BACKEND_CONST, MET_BACKEND_NETCDF,    &
                                    SWPART_CLEARIDX, SWPART_WEISS_NORMAN, INTERP_LINEAR,        &
                                    INTERP_STEP, METAVG_END, METAVG_BEGIN, SWPART_PASSTHROUGH,  &
                                    CLAMP_HOLD, CLAMP_ERROR,                                   &
                                    GRIDMATCH_EXPLICIT, GRIDMATCH_NEAREST, CO2_SOURCE_FILE,      &
                                    HEIGHT_ABOVE_ZERO_PLANE, HEIGHT_ABOVE_GROUND,               &
                                    WIND_EXPOSURE_OPEN_TERRAIN, WIND_EXPOSURE_LOCAL
   use meds_forcing_types,   only : met_forcing_t, met_source_t, met_cursor_t
   use meds_forcing_kernels, only : interpolate_forcing, dewpoint_to_specific_humidity,         &
                                    rh_to_specific_humidity, precip_phase, partition_shortwave, &
                                    met_solar_cosz, cosz_reconstruct_factor, disaggregate_sw,   &
                                    great_circle_distance, nearest_grid_index,                   &
                                    clearness_index, clear_sky_emissivity, synthesize_lwdown
   use meds_lapse_rate,      only : wind_log_profile, lapse_air_temperature, lapse_pressure,        &
                                    monthly_lapse_rate, lapse_specific_humidity, lapse_longwave,    &
                                    cas_top_wind_factor, cas_top_air_temperature, met_to_cas_top
   use meds_constants,       only : grav, cp_air
   use meds_met_driver,      only : met_open, met_cursor_init, met_advance, met_instant, met_close, &
                                   MET_OK, MET_ERR_WINDOW_NOT_WHOLE_YEARS,                     &
                                   MET_ERR_START_NOT_A_RECORD, MET_ERR_WINDOW_NOT_COVERED,       &
                                   MET_ERR_DT_MISMATCH, MET_ERR_ATTR_MISMATCH, MET_ERR_CO2_FILE,  &
                                   MET_ERR_CO2_NOT_COVERED, MET_ERR_CO2_IN_MET_FILE
   use meds_netcdf_c
   use iso_c_binding,        only : c_int, c_size_t, c_double
   implicit none
   character(len=*), parameter :: NCFILE = 'test_met_driver_tmp.nc'
   character(len=1024) :: co2_default_file = ''

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
   call test_terrain_lapse()
   call test_cas_top()
   call test_multiyear_cycling()
   call test_lwdown_synthesis()
   call test_recycle_anchor_phase()
   call test_file_config_agreement()
   if (command_argument_count() >= 1) call get_command_argument(1, co2_default_file)
   call test_co2_series(trim(co2_default_file))

   call test_report('test_met_driver')

contains

   !----- The site's cursor into an opened source: record-by-record access goes through a cursor   !
   !      (MEDS_POLYGON_RUNTIME_PLAN.md R2); a site run has one, at cell 1 and the [site] location. -!
   subroutine site_cursor(src, cur, fc)
      type(met_source_t),     intent(in)  :: src
      type(met_cursor_t),     intent(out) :: cur
      type(forcing_config_t), intent(in)  :: fc
      call met_cursor_init(src, cur, 1_ik, fc%latitude_deg, fc%longitude_deg, fc%utc_offset_h,  &
                           fc%elevation_m)
   end subroutine site_cursor

   !----- LONGWAVE SYNTHESIS (#182). The clear-sky emissivities are closed-form, so these are     !
   !      known-answer checks against Brutsaert (1975) and Idso & Jackson (1969) evaluated by      !
   !      hand -- not against the code's own output. What they pin beyond arithmetic: the cloud    !
   !      term only ever ADDS, the night sentinel is distinguishable from an overcast sky, and     !
   !      both emissivity forms stay inside physical bounds.                                       !
   subroutine test_lwdown_synthesis()
      real(wp) :: eps, lw_clear, lw_cloudy, kt
      real(wp), parameter :: SIG = 5.670374419e-8_wp
      print '(a)', '-- 12. LWdown synthesis (Brutsaert / Idso + cloud term) --'

      !----- Brutsaert at 20 C, q = 0.010, 1013.25 hPa: eps = 1.24 (e/T)^(1/7) = 0.81985. -------!
      eps = clear_sky_emissivity(LW_CLEAR_BRUTSAERT, 293.15_wp, 0.010_wp, 101325.0_wp)
      call check('Brutsaert eps at 20 C, q = 0.010', eps, 0.81985_wp, 1.0e-4_wp)
      !----- Idso-Jackson is temperature only: 1 - 0.261 exp(-7.77e-4 dT^2) = 0.80866. ----------!
      eps = clear_sky_emissivity(LW_CLEAR_IDSO, 293.15_wp, 0.010_wp, 101325.0_wp)
      call check('Idso eps at 20 C (humidity-independent)', eps, 0.80866_wp, 1.0e-4_wp)
      call check('Idso ignores humidity entirely',                                          &
                       clear_sky_emissivity(LW_CLEAR_IDSO, 293.15_wp, 0.001_wp, 101325.0_wp),     &
                       clear_sky_emissivity(LW_CLEAR_IDSO, 293.15_wp, 0.030_wp, 101325.0_wp), 1.0e-12_wp)
      !----- Brutsaert does NOT: a moister sky emits more. --------------------------------------!
      call check_true('Brutsaert: a moister sky is more emissive',                                &
             clear_sky_emissivity(LW_CLEAR_BRUTSAERT, 293.15_wp, 0.030_wp, 101325.0_wp) >         &
             clear_sky_emissivity(LW_CLEAR_BRUTSAERT, 293.15_wp, 0.001_wp, 101325.0_wp))

      !----- The full synthesis, clear sky (kt = 1): LW = eps sigma T^4 = 343.33 W/m2. ----------!
      lw_clear = synthesize_lwdown(LW_CLEAR_BRUTSAERT, 293.15_wp, 0.010_wp, 101325.0_wp,          &
                                   1.0_wp, 0.22_wp)
      call check('clear-sky LWdown at 20 C', lw_clear, 0.81985_wp * SIG * 293.15_wp**4, 1.0e-2_wp)
      !----- Overcast (kt = 0.3) with a = 0.22: x (1 + 0.22*0.7) = 396.20 W/m2. -----------------!
      lw_cloudy = synthesize_lwdown(LW_CLEAR_BRUTSAERT, 293.15_wp, 0.010_wp, 101325.0_wp,         &
                                    0.3_wp, 0.22_wp)
      call check('cloudy-sky LWdown (kt = 0.3, a = 0.22)', lw_cloudy, 396.199_wp, 1.0e-2_wp)
      call check_true('the cloud term only ADDS to the clear-sky flux', lw_cloudy > lw_clear)
      !----- a = 0 recovers the clear-sky value for ANY kt: the term is genuinely optional. -----!
      call check('cloud coefficient 0 => clear sky at any kt',                              &
                       synthesize_lwdown(LW_CLEAR_BRUTSAERT, 293.15_wp, 0.010_wp, 101325.0_wp,    &
                                         0.0_wp, 0.0_wp), lw_clear, 1.0e-9_wp)

      !----- The clearness index: a real value by day, a NEGATIVE sentinel at night. A zero here  !
      !      would read as fully overcast and give every night the maximum cloud correction. -----!
      kt = clearness_index(500.0_wp, 0.5_wp)
      call check('clearness index by day = SW/(S0 cosz)', kt, 500.0_wp / (1361.0_wp*0.5_wp), 1.0e-9_wp)
      call check_true('clearness index is a NEGATIVE sentinel at night, not 0',                    &
                      clearness_index(0.0_wp, 0.0_wp) < 0.0_wp, clearness_index(0.0_wp, 0.0_wp))
      call check_true('clearness index is capped at 1', clearness_index(5000.0_wp, 0.5_wp) <= 1.0_wp)
   end subroutine test_lwdown_synthesis



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

   !----- 3. rainfall phase split (mass-conserving). -----------------------------------------!
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
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met
      type(meds_time_t)      :: now
      print '(a)', '-- test 6: CONST backend --'
      fc%backend = MET_BACKEND_CONST ; fc%latitude_deg = 42.44_wp ; fc%longitude_deg = -76.50_wp
      call met_open(src, fc) ; call site_cursor(src, cur, fc)
      now = meds_time_t(year=2020_ik, month=7_ik, day=1_ik, hour=17_ik)   ! ~solar noon at Ithaca
      call met_advance(src, cur, now)
      met = met_instant(src, cur, now)
      call check('CONST swdown = 400', met%swdown(), 400.0_wp, 1.0e-9_wp)
      call check_true('CONST cosz recomputed > 0 at noon', met%cosz > 0.5_wp, met%cosz)
      call check_true('CONST rho_air ~ 1.2', abs(met%rho_air - 1.2_wp) < 0.2_wp, met%rho_air)
      call met_close(src)
   end subroutine test_const_backend

   !----- 7. NetCDF round-trip: write a (time=25, grid=2) file, read it back. ----------------!
   subroutine test_netcdf_roundtrip()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
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
      call met_open(src, fc) ; call site_cursor(src, cur, fc)
      call check_true('opened file: nrec=25', src%nrec == 25_ik, real(src%nrec, wp))
      call check_true('opened file: ngrid=2', src%ngrid == 2_ik, real(src%ngrid, wp))
      call check_true('a file with Wind only supplies no wind vector', .not. src%has_wind_vector)

      ! noon (17 UTC): SW positive, Tair follows the diurnal input
      now = time_advance_seconds(base, 17.0_wp*3600.0_wp)
      call met_advance(src, cur, now) ; met_noon = met_instant(src, cur, now)
      call check_true('noon SWdown > 0', met_noon%swdown() > 100.0_wp, met_noon%swdown())
      call check_true('noon cosz > 0.5', met_noon%cosz > 0.5_wp, met_noon%cosz)
      call check_true('Tair in [285,305]', met_noon%tair_k > 285.0_wp .and. met_noon%tair_k < 305.0_wp, met_noon%tair_k)
      sw_g1 = met_noon%swdown()

      ! night (04 UTC): SW exactly 0
      now = time_advance_seconds(base, 4.0_wp*3600.0_wp)
      call met_advance(src, cur, now) ; met_night = met_instant(src, cur, now)
      call check('night SWdown = 0', met_night%swdown(), 0.0_wp, 1.0e-12_wp)
      call check_true('an instant from a Wind-only file has no vector', .not. met_night%has_wind_vector)
      call met_close(src)

      !----- A file carrying u10 and v10 supplies the vector, and the speed comes from it (§7.1). -!
      call write_synthetic_forcing(NCFILE, base, with_components=.true.)
      fc%grid_index = 1_ik
      call met_open(src, fc) ; call site_cursor(src, cur, fc)
      call check_true('a file with u10 and v10 supplies the wind vector', src%has_wind_vector)
      now = time_advance_seconds(base, 6.0_wp*3600.0_wp)
      call met_advance(src, cur, now) ; met_night = met_instant(src, cur, now)
      call check('wind_u = u10', met_night%wind_u, 3.0_wp, 1.0e-12_wp)
      call check('wind_v = v10', met_night%wind_v, -4.0_wp, 1.0e-12_wp)
      call check('speed from the components, not Wind', met_night%wind, 5.0_wp, 1.0e-12_wp)
      call met_close(src)
      call write_synthetic_forcing(NCFILE, base)

      ! multi-grid: grid_index=2 carries a distinct (scaled) SW series
      fc%grid_index = 2_ik
      call met_open(src, fc) ; call site_cursor(src, cur, fc)
      now = time_advance_seconds(base, 17.0_wp*3600.0_wp)
      call met_advance(src, cur, now) ; met_g2 = met_instant(src, cur, now)
      call check_true('grid 2 SW differs from grid 1', abs(met_g2%swdown() - sw_g1) > 1.0_wp, met_g2%swdown()-sw_g1)
      call met_close(src)

      call test_edge_paths(base)
   end subroutine test_netcdf_roundtrip

   !----- 13. THE FILE AND THE CONFIG MUST AGREE (#185). Three things the reader used to take on   !
   !      trust: the record spacing (dt_forcing was read straight from the config and never checked !
   !      against the file at all), and the two self-describing global attributes the prep script    !
   !      writes and nothing read. A file that disagreed with its config was silently mis-timed or   !
   !      mis-partitioned -- no later check can see either, because both produce a plausible run.    !
   subroutine test_file_config_agreement()
      type(met_source_t)     :: src
      type(forcing_config_t) :: fc
      type(meds_time_t)      :: base
      integer(ik)            :: st
      print '(a)', '-- test 13: the forcing file must agree with [forcing] --'
      base = meds_time_t(year=2020_ik, month=7_ik, day=1_ik)
      call write_synthetic_forcing(NCFILE, base)

      fc%backend = MET_BACKEND_NETCDF ; fc%path = NCFILE ; fc%grid_index = 1_ik
      fc%dt_forcing = 3600.0_wp ; fc%avg_convention = METAVG_END ; fc%sw_partition = SWPART_CLEARIDX
      fc%latitude_deg = 42.44_wp ; fc%longitude_deg = -76.50_wp ; fc%utc_offset_h = 0.0_wp
      fc%apply_solar_longitude = .true. ; fc%recycle = .false.

      !----- the honest config opens cleanly. A validator that rejects the good case is useless. --!
      call met_open(src, fc, stat=st)
      call check_true('a matching config opens (MET_OK)', st == MET_OK, real(st, wp))
      call met_close(src)

      !----- dt_forcing lying about the cadence. THIS is the one nothing checked: dt_forcing places !
      !      interval midpoints, disaggregates shortwave and brackets the recycle seam. -------------!
      fc%dt_forcing = 1800.0_wp
      call met_open(src, fc, stat=st)
      call check_true('a half-hourly config against an hourly file is rejected',                   &
                      st == MET_ERR_DT_MISMATCH, real(st, wp))
      fc%dt_forcing = 3600.0_wp

      !----- the file says its flux means end at the stamp; the config says they begin there. ------!
      fc%avg_convention = METAVG_BEGIN
      call met_open(src, fc, stat=st)
      call check_true('avg_convention contradicting the file attribute is rejected',               &
                      st == MET_ERR_ATTR_MISMATCH, real(st, wp))
      fc%avg_convention = METAVG_END

      !----- the file carries TOTAL shortwave; passthrough would read component fields it lacks. ---!
      fc%sw_partition = SWPART_PASSTHROUGH
      call met_open(src, fc, stat=st)
      call check_true('passthrough against a total-shortwave file is rejected',                    &
                      st == MET_ERR_ATTR_MISMATCH, real(st, wp))
      fc%sw_partition = SWPART_CLEARIDX
   end subroutine test_file_config_agreement

   !----- 8. Edge paths the review flagged (regression tests for the CLAMP_HOLD + recycle fixes). !
   subroutine test_edge_paths(base)
      type(meds_time_t), intent(in) :: base
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
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
      call met_open(src, fc) ; call site_cursor(src, cur, fc)
      now = time_advance_seconds(base, -12.0_wp*3600.0_wp)          ! 12 h before the first record
      call met_advance(src, cur, now) ; met = met_instant(src, cur, now)
      call check('hold before t0 -> record 1', met%tair_k, tair0, 1.0e-4_wp)
      now = time_advance_seconds(base, 0.5_wp*3600.0_wp)            ! 00:30 -> midpoint of [t0,t1)
      call met_advance(src, cur, now) ; met = met_instant(src, cur, now)
      expect_hold = 0.5_wp*(tair0 + tair1)
      call check('march into [t0,t1) interpolates (not stuck at rec1)', met%tair_k, expect_hold, 1.0e-4_wp)
      call met_close(src)

      ! (b) A 25-record (24 h span) file CANNOT be recycled: a sub-year window drifts both
      !     hour-of-day and day-of-year on every wrap. This used to fall through to an
      !     absolute-seconds span-wrap SILENTLY, which is exactly how a 29-yr production run came
      !     to be driven by ~10 h-anti-phased sub-daily shortwave. It must now be REJECTED.
      fc%start_clamp = CLAMP_ERROR ; fc%recycle = .true.
      fc%recycle_start = base
      fc%recycle_end   = time_advance_seconds(base, 24.0_wp*3600.0_wp)   ! 1 day, not a whole year
      call met_open(src, fc, stat=st)
      call check_true('sub-year recycle window rejected (not silently span-wrapped)',            &
                      st == MET_ERR_WINDOW_NOT_WHOLE_YEARS, real(st, wp))

      ! (c) A whole-year window whose start does NOT match a record stamp is rejected too --
      !     MEDS never guesses the start of the day (e.g. a config saying 00:00 against an
      !     end-of-interval ERA5-Land file whose records are stamped 01:00).
      fc%recycle_start = time_advance_seconds(base, 1800.0_wp)           ! 00:30, between records
      fc%recycle_end   = meds_time_t(2021_ik, 7_ik, 1_ik, 0_ik, 30_ik)
      call met_open(src, fc, stat=st)
      call check_true('recycle_start off the record grid rejected',                              &
                      st == MET_ERR_START_NOT_A_RECORD, real(st, wp))

      ! (d) A whole-year window ON the record grid, but the file only holds 24 h of it.
      fc%recycle_start = base
      fc%recycle_end   = meds_time_t(2021_ik, 7_ik, 1_ik)
      call met_open(src, fc, stat=st)
      call check_true('window not covered by the file rejected',                                 &
                      st == MET_ERR_WINDOW_NOT_COVERED, real(st, wp))
      call met_close(src)
   end subroutine test_edge_paths

   !----- The recycle mapping must preserve hour-of-day and day-of-year EXACTLY, for an anchor    !
   !      anywhere in the calendar -- not only Jan-1. This is the regression for the met-recycle    !
   !      phase bug: the old classifier accepted only Jan-1 00:00 files and silently span-wrapped    !
   !      everything else, which shifts hour-of-day whenever the span is not a whole number of days. !
   subroutine test_recycle_anchor_phase()
      character(len=*), parameter :: YF = 'test_met_anchorfile_tmp.nc'
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
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
      call met_open(src, fc, stat=st)
      call check_true('mid-year window beyond the file is rejected',                             &
                      st == MET_ERR_WINDOW_NOT_COVERED, real(st, wp))

      !----- Whole-year window on the file: every model year must read the SAME calendar day, for
      !      model years both before and after the window (the modulo must not go negative).
      fc%recycle_start = meds_time_t(2021_ik, 1_ik, 1_ik)
      fc%recycle_end   = meds_time_t(2022_ik, 1_ik, 1_ik)
      call met_open(src, fc, stat=st)
      call check_true('whole-year window accepted', st == MET_OK, real(st, wp))
      call site_cursor(src, cur, fc)
      call met_advance(src, cur, meds_time_t(2021_ik,9_ik,17_ik)) ; m_ref = met_instant(src, cur, meds_time_t(2021_ik,9_ik,17_ik))
      call met_advance(src, cur, meds_time_t(2049_ik,9_ik,17_ik)) ; m_map = met_instant(src, cur, meds_time_t(2049_ik,9_ik,17_ik))
      call check('29 yr later reads the same calendar day', m_map%tair_k, m_ref%tair_k, 1.0e-9_wp)
      call met_advance(src, cur, meds_time_t(2013_ik,9_ik,17_ik)) ; m_map = met_instant(src, cur, meds_time_t(2013_ik,9_ik,17_ik))
      call check('a model year BEFORE the window maps in too', m_map%tair_k, m_ref%tair_k, 1.0e-9_wp)
      call met_close(src)
   end subroutine test_recycle_anchor_phase

   !----- Write a synthetic (time=25 hourly, grid=2) MEDS forcing NetCDF via meds_netcdf_c. ----!
   subroutine write_synthetic_forcing(path, base, with_components, with_co2air)
      character(len=*),  intent(in) :: path
      type(meds_time_t), intent(in) :: base
      logical, optional, intent(in) :: with_components   !< also write u10 = 3, v10 = -4 (speed 5, not Wind's 3)
      logical, optional, intent(in) :: with_co2air       !< also write CO2air, which the reader rejects (#184)
      integer, parameter :: NT = 25, NG = 2
      integer(c_int) :: st, ncid, td, gd, vt, vla, vlo, vv(7), vu, vw, vc
      integer(c_int) :: dims2(2), dims1(1)
      real(c_double) :: tsec(NT), lat(NG), lon(NG)
      real(c_double) :: dat(NG, NT)        ! (grid, time) column-major == C [time][grid]
      integer :: it, ig, k
      real(wp) :: hh, sw_mean
      character(len=8), parameter :: vnames(7) = ['Tair    ','Qair    ','PSurf   ','Wind    ', &
                                                  'Rainf   ','SWdown  ','LWdown  ']
      integer(c_size_t) :: start2(2), count2(2), start1(1), count1(1)
      logical :: co2air

      co2air = .false.
      if (present(with_co2air)) co2air = with_co2air
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
      do k = 1, 7
         st = nc_def_var_f(ncid, trim(vnames(k)), NC_DOUBLE, 2, dims2, vv(k))
         call nc_check(st, 'write: var '//trim(vnames(k)))
      end do
      if (co2air) then
         st = nc_def_var_f(ncid, 'CO2air', NC_DOUBLE, 2, dims2, vc) ; call nc_check(st, 'write: var CO2air')
      end if
      if (present(with_components)) then
         if (with_components) then
            st = nc_def_var_f(ncid, 'u10', NC_DOUBLE, 2, dims2, vu) ; call nc_check(st, 'write: var u10')
            st = nc_def_var_f(ncid, 'v10', NC_DOUBLE, 2, dims2, vw) ; call nc_check(st, 'write: var v10')
         end if
      end if
      !----- The two SELF-DESCRIBING global attributes the prep script writes. The fixture has to    !
      !      carry them or it is not the file the reader validates against, and the #185 checks      !
      !      would be skipped here while firing in production -- the fixture-must-mirror-the-driver  !
      !      trap. This file is hourly, ERA5-Land-style, with total shortwave.                        !
      st = nc_put_att_text_f(ncid, NC_GLOBAL, 'avg_convention', int(len_trim('end'), c_size_t), 'end')
      call nc_check(st, 'write: avg_convention')
      st = nc_put_att_text_f(ncid, NC_GLOBAL, 'sw_input_kind', int(len_trim('total'), c_size_t), 'total')
      call nc_check(st, 'write: sw_input_kind')
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
      do k = 1, 7
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
               end select
            end do
         end do
         st = nc_put_vara_double(ncid, vv(k), start2, count2, dat)
         call nc_check(st, 'write: vals '//trim(vnames(k)))
      end do
      if (present(with_components)) then
         if (with_components) then
            dat = 3.0_wp
            st = nc_put_vara_double(ncid, vu, start2, count2, dat) ; call nc_check(st, 'write: vals u10')
            dat = -4.0_wp
            st = nc_put_vara_double(ncid, vw, start2, count2, dat) ; call nc_check(st, 'write: vals v10')
         end if
      end if
      if (co2air) then
         dat = 415.0_wp
         st = nc_put_vara_double(ncid, vc, start2, count2, dat) ; call nc_check(st, 'write: vals CO2air')
      end if
      st = nc_close(ncid) ; call nc_check(st, 'write: close')
   end subroutine write_synthetic_forcing

   !----- Nearest-grid match: pure kernel (argmin + great-circle) and the reader override. --------!
   subroutine test_nearest_grid()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
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
      call met_open(src, fc) ; call site_cursor(src, cur, fc)
      call check('reader resolves nearest -> grid_index 2', real(src%grid_index,wp), 2.0_wp, 0.5_wp)
      call met_close(src)
      fc%latitude_deg = 42.51_wp ; fc%longitude_deg = -76.61_wp     ! nearest to cell 1
      call met_open(src, fc) ; call site_cursor(src, cur, fc)
      call check('reader resolves nearest -> grid_index 1', real(src%grid_index,wp), 1.0_wp, 0.5_wp)
      call met_close(src)
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

   !----- The terrain lapse's humidity, longwave and monthly rate (docs/science/forcing.md §8). ---!
   subroutine test_terrain_lapse()
      real(wp), parameter :: T0 = 285.0_wp, P0 = 95000.0_wp, Q0 = 0.007_wp, DZ = -500.0_wp, G = 0.0065_wp
      real(wp) :: t1, p1, q1, e0, e1, lw1, gm(12)
      integer  :: k
      print '(a)', '-- test: terrain lapse (RH-held humidity, longwave, monthly rate) --'
      t1 = lapse_air_temperature(T0, DZ, G) ; p1 = lapse_pressure(P0, T0, DZ, G)
      q1 = lapse_specific_humidity(Q0, T0, P0, t1, p1)
      e0 = Q0 * P0 / (0.622_wp + 0.378_wp * Q0) ; e1 = q1 * p1 / (0.622_wp + 0.378_wp * q1)
      call check('relative humidity is held', e1 / sat_vapor_pressure(t1), e0 / sat_vapor_pressure(T0), 1.0e-12_wp)
      call check_true('a site 500 m below its cell is warmer and holds more vapour', q1 > Q0, q1 - Q0)
      call check('dz = 0 leaves q unchanged', lapse_specific_humidity(Q0, T0, P0, T0, P0), Q0, 1.0e-15_wp)
      lw1 = lapse_longwave(300.0_wp, LW_CLEAR_BRUTSAERT, T0, Q0, P0, t1, q1, p1)
      call check('file longwave scales by eps*T^4', lw1,                                           &
                 300.0_wp * clear_sky_emissivity(LW_CLEAR_BRUTSAERT, t1, q1, p1) * t1**4           &
                          / (clear_sky_emissivity(LW_CLEAR_BRUTSAERT, T0, Q0, P0) * T0**4), 1.0e-9_wp)
      call check_true('the warmer, moister site gets more longwave', lw1 > 300.0_wp, lw1 - 300.0_wp)
      call check('dz = 0 leaves longwave unchanged',                                               &
                 lapse_longwave(300.0_wp, LW_CLEAR_IDSO, T0, Q0, P0, T0, Q0, P0), 300.0_wp, 1.0e-9_wp)
      gm = [(0.001_wp * real(k, wp), k = 1, 12)]
      call check('monthly rate: July', monthly_lapse_rate(gm, 7_ik), 0.007_wp, 1.0e-15_wp)
      call check('monthly rate: December', monthly_lapse_rate(gm, 12_ik), 0.012_wp, 1.0e-15_wp)
   end subroutine test_terrain_lapse

   !----- The move to a patch's canopy-air top (docs/science/forcing.md §8). A 25 m canopy:         !
   !      d = 15.75 m, z0 = 3.25 m, top z_c = 30 m.                                                !
   subroutine test_cas_top()
      real(wp), parameter :: ZC = 30.0_wp, D = 15.75_wp, Z0 = 3.25_wp, GC = grav / cp_air
      type(forcing_config_t) :: f
      type(met_forcing_t)    :: m, top
      real(wp) :: fac
      print '(a)', '-- test: forcing to the canopy-air top --'
      !----- ERA5-Land: open-terrain 10 m wind, heights above the zero plane. -------------------!
      f%tq_height = 2.0_wp ; f%wind_height = 10.0_wp ; f%height_above = HEIGHT_ABOVE_ZERO_PLANE
      f%wind_exposure = WIND_EXPOSURE_OPEN_TERRAIN ; f%wind_exposure_z0 = 0.03_wp ; f%wind_blending_height = 40.0_wp
      fac = cas_top_wind_factor(f, ZC, D, Z0)
      call check('ERA5 wind: open-terrain step, then the patch profile', fac,                     &
                 log(40.0_wp/0.03_wp) / log(10.0_wp/0.03_wp) * log((ZC - D)/Z0) / log(40.0_wp/Z0), 1.0e-12_wp)
      call check('ERA5 wind: ~0.73 u10 at the top of a 25 m canopy', fac, 0.729_wp, 1.0e-3_wp)
      call check('temperature conserves theta from d + 2 m', cas_top_air_temperature(f, 290.0_wp, ZC, D), &
                 290.0_wp - GC * (ZC - (D + 2.0_wp)), 1.0e-12_wp)
      call check('theta at the top equals theta at the forcing height',                            &
                 cas_top_air_temperature(f, 290.0_wp, ZC, D) + GC * ZC, 290.0_wp + GC * (D + 2.0_wp), 1.0e-9_wp)
      !----- A flux tower: local wind at 40 m above the ground; the 2 z0 floor on a tall stand. --!
      f%height_above = HEIGHT_ABOVE_GROUND ; f%wind_exposure = WIND_EXPOSURE_LOCAL
      f%tq_height = 40.0_wp ; f%wind_height = 40.0_wp
      call check('tower wind: the patch profile from 40 m above the ground', cas_top_wind_factor(f, ZC, D, Z0), &
                 log((ZC - D)/Z0) / log((40.0_wp - D)/Z0), 1.0e-12_wp)
      call check('tower temperature: along the dry adiabat from 40 m', cas_top_air_temperature(f, 290.0_wp, ZC, D), &
                 290.0_wp + GC * 10.0_wp, 1.0e-12_wp)
      call check('tower at the top: no change', cas_top_wind_factor(f, 40.0_wp, D, Z0), 1.0_wp, 1.0e-15_wp)
      call check('heights floored at 2 z0 above d', cas_top_wind_factor(f, D + 1.0_wp, D, Z0),       &
                 log(2.0_wp) / log((40.0_wp - D)/Z0), 1.0e-12_wp)
      !----- The whole record: wind and vector scaled, temperature moved, rho re-derived, the rest kept. -!
      f%height_above = HEIGHT_ABOVE_ZERO_PLANE ; f%wind_exposure = WIND_EXPOSURE_OPEN_TERRAIN
      f%tq_height = 2.0_wp ; f%wind_height = 10.0_wp
      m%wind = 5.0_wp ; m%wind_u = 3.0_wp ; m%wind_v = -4.0_wp ; m%has_wind_vector = .true.
      m%tair_k = 290.0_wp ; m%qair = 0.009_wp ; m%psurf_pa = 98000.0_wp
      top = met_to_cas_top(m, f, ZC, D, Z0)
      fac = cas_top_wind_factor(f, ZC, D, Z0)
      call check('record: wind speed scaled', top%wind, 5.0_wp * fac, 1.0e-12_wp)
      call check('record: vector scaled (u)', top%wind_u, 3.0_wp * fac, 1.0e-12_wp)
      call check('record: vector scaled (v)', top%wind_v, -4.0_wp * fac, 1.0e-12_wp)
      call check('record: temperature moved', top%tair_k, cas_top_air_temperature(f, 290.0_wp, ZC, D), 1.0e-12_wp)
      call check('record: humidity conserved', top%qair, 0.009_wp, 1.0e-15_wp)
      call check('record: pressure stays at the ground', top%psurf_pa, 98000.0_wp, 1.0e-12_wp)
      call check('record: shortwave unchanged', top%swdown(), m%swdown(), 1.0e-12_wp)
      call check('record: air density re-derived', top%rho_air, air_density(top%tair_k, 98000.0_wp, 0.009_wp), 1.0e-12_wp)
   end subroutine test_cas_top

   !----- Multi-year CALENDAR recycling + Feb-29 reconciliation (whole-year daily file). ----------!
   subroutine test_multiyear_cycling()
      character(len=*), parameter :: YF = 'test_met_yearfile_tmp.nc'
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
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
      call met_open(src, fc) ; call site_cursor(src, cur, fc)
      call check('n_cycle_years = 1', real(src%n_cycle_years,wp), 1.0_wp, 0.5_wp)
      call check('cycle anchor year = 2021', real(src%cycle_anchor%year,wp), 2021.0_wp, 0.5_wp)
      call check('cycle first record = #1', real(src%irec_cycle_first,wp), 1.0_wp, 0.5_wp)
      !----- recycle identity: model 2023-07-01 reads the 2021-07-01 record (day-of-year exact). --!
      call met_advance(src, cur, meds_time_t(2021_ik,7_ik,1_ik)) ; m_ref = met_instant(src, cur, meds_time_t(2021_ik,7_ik,1_ik))
      call met_advance(src, cur, meds_time_t(2023_ik,7_ik,1_ik)) ; m_map = met_instant(src, cur, meds_time_t(2023_ik,7_ik,1_ik))
      call check('recycle: 2023-07-01 reads the 2021-07-01 record', m_map%tair_k, m_ref%tair_k, 1.0e-9_wp)
      !----- Feb-29 reconciliation: model 2024-02-29 (leap) maps to file 2021-02-28. --------------!
      call met_advance(src, cur, meds_time_t(2021_ik,2_ik,28_ik)) ; m_ref = met_instant(src, cur, meds_time_t(2021_ik,2_ik,28_ik))
      call met_advance(src, cur, meds_time_t(2024_ik,2_ik,29_ik)) ; m_map = met_instant(src, cur, meds_time_t(2024_ik,2_ik,29_ik))
      call check('Feb-29 (leap model) maps to file Feb-28', m_map%tair_k, m_ref%tair_k, 1.0e-9_wp)
      call met_close(src)
   end subroutine test_multiyear_cycling

   !----- 15. Prescribed CO2 (#184): co2_const on every backend, or a MEDS CO2 file read on MODEL  !
   !      time -- each value at its period's middle, linear between middles, the end values held  !
   !      over the outer half-periods, a run the file does not cover rejected at open.             !
   subroutine test_co2_series(default_file)
      character(len=*), intent(in) :: default_file      !< the shipped series ('' -> its check is skipped)
      character(len=*), parameter :: CF = 'test_met_co2_tmp.txt', YF = 'test_met_co2_year_tmp.nc'
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met, m_ref, m_map
      integer(ik)            :: st
      print '(a)', '-- test 15: prescribed CO2 (#184) --'

      !----- "const": co2_const, the CONST backend included (it used to return the type's 420). -!
      fc%backend = MET_BACKEND_CONST ; fc%co2_const = 284.3_wp
      call met_open(src, fc) ; call site_cursor(src, cur, fc)
      met = met_instant(src, cur, meds_time_t(2020_ik, 7_ik, 1_ik))
      call check('const: CO2 = co2_const on the CONST backend', met%co2, 284.3_wp, 1.0e-12_wp)
      call met_close(src)

      !----- Annual means, with comments and a blank line. 2000 is a leap year, so its middle is   !
      !      Jul 2 00:00; 2001's is Jul 2 12:00. On 2001-01-01 the weight is 183 d / 365.5 d. -------!
      call write_text(CF, [character(len=48) :: '# a comment', 'timestep  1 year', 'units umol/mol', &
                           '2000  300.0   # a trailing comment', '2001  310.0', '', '2002  330.0'])
      call open_file_co2(CF, meds_time_t(2000_ik,1_ik,1_ik), meds_time_t(2003_ik,1_ik,1_ik), src, cur, st)
      call check_true('year: the file opens and covers 2000-01-01 .. 2003-01-01', st == MET_OK, real(st, wp))
      call check('year: at the middle of 2001, the 2001 value', co2_now(src, cur,                    &
                 meds_time_t(2001_ik,7_ik,2_ik,12_ik)), 310.0_wp, 1.0e-9_wp)
      call check('year: linear between middles', co2_now(src, cur, meds_time_t(2001_ik,1_ik,1_ik)), &
                 300.0_wp + 10.0_wp * 183.0_wp / 365.5_wp, 1.0e-9_wp)
      call check('year: held before the first middle', co2_now(src, cur,                           &
                 meds_time_t(2000_ik,3_ik,1_ik)), 300.0_wp, 1.0e-12_wp)
      call check('year: held after the last middle', co2_now(src, cur,                             &
                 meds_time_t(2002_ik,12_ik,31_ik,23_ik)), 330.0_wp, 1.0e-12_wp)
      call met_close(src)
      call open_file_co2(CF, meds_time_t(1999_ik,12_ik,31_ik), meds_time_t(2003_ik,1_ik,1_ik), src, cur, st)
      call check_true('year: a run starting before the file is rejected', st == MET_ERR_CO2_NOT_COVERED, real(st, wp))
      call open_file_co2(CF, meds_time_t(2000_ik,1_ik,1_ik), meds_time_t(2003_ik,1_ik,1_ik,1_ik), src, cur, st)
      call check_true('year: a run ending after the file is rejected', st == MET_ERR_CO2_NOT_COVERED, real(st, wp))

      !----- Monthly, across a leap February: Jan's middle is the 16th 12:00 (15.5 d), Feb's the   !
      !      15th 12:00 (14.5 d), so Feb 1 00:00 is 15.5/30 of the way. -----------------------------!
      call write_text(CF, [character(len=48) :: 'timestep 1 month', 'units umol/mol', '2024-01 400', &
                           '2024-02 410', '2024-03 420'])
      call open_file_co2(CF, meds_time_t(2024_ik,1_ik,1_ik), meds_time_t(2024_ik,4_ik,1_ik), src, cur, st)
      call check_true('month: opens', st == MET_OK, real(st, wp))
      call check('month: at the middle of a leap February', co2_now(src, cur,                       &
                 meds_time_t(2024_ik,2_ik,15_ik,12_ik)), 410.0_wp, 1.0e-9_wp)
      call check('month: Feb 1, linear by calendar days', co2_now(src, cur,                         &
                 meds_time_t(2024_ik,2_ik,1_ik)), 400.0_wp + 10.0_wp * 15.5_wp / 30.0_wp, 1.0e-9_wp)
      call met_close(src)

      !----- 5-year means from the year 1000: both spans are 1826 d, so 1005-01-01 is halfway. -----!
      call write_text(CF, [character(len=48) :: 'timestep 5 year', 'units umol/mol', '1000 280', '1005 290'])
      call open_file_co2(CF, meds_time_t(1000_ik,1_ik,1_ik), meds_time_t(1010_ik,1_ik,1_ik), src, cur, st)
      call check_true('5 year: opens and covers 1000 .. 1010', st == MET_OK, real(st, wp))
      call check('5 year: halfway between the two middles', co2_now(src, cur,                       &
                 meds_time_t(1005_ik,1_ik,1_ik)), 285.0_wp, 1.0e-9_wp)
      call met_close(src)

      !----- Daily, and 30-minute (a flux tower's cadence). ------------------------------------!
      call write_text(CF, [character(len=48) :: 'timestep 1 day', 'units umol/mol', '2024-02-28 400', &
                           '2024-02-29 402'])
      call open_file_co2(CF, meds_time_t(2024_ik,2_ik,28_ik), meds_time_t(2024_ik,3_ik,1_ik), src, cur, st)
      call check_true('day: opens (a leap day row)', st == MET_OK, real(st, wp))
      call check('day: midnight between two daily means', co2_now(src, cur,                         &
                 meds_time_t(2024_ik,2_ik,29_ik)), 401.0_wp, 1.0e-9_wp)
      call met_close(src)
      call write_text(CF, [character(len=48) :: 'timestep 30 minute', 'units umol/mol',             &
                           '2024-07-01T00:00 400', '2024-07-01T00:30 402', '2024-07-01T01:00 404'])
      call open_file_co2(CF, meds_time_t(2024_ik,7_ik,1_ik), meds_time_t(2024_ik,7_ik,1_ik,1_ik,30_ik), src, cur, st)
      call check_true('30 minute: opens and covers 00:00 .. 01:30', st == MET_OK, real(st, wp))
      call check('30 minute: between the 00:15 and 00:45 middles', co2_now(src, cur,                &
                 meds_time_t(2024_ik,7_ik,1_ik,0_ik,30_ik)), 401.0_wp, 1.0e-9_wp)
      call check('30 minute: held over the last half-period', co2_now(src, cur,                     &
                 meds_time_t(2024_ik,7_ik,1_ik,1_ik,25_ik)), 404.0_wp, 1.0e-12_wp)
      call met_close(src)
      call open_file_co2(CF, meds_time_t(2024_ik,7_ik,1_ik), meds_time_t(2024_ik,7_ik,1_ik,1_ik,31_ik), src, cur, st)
      call check_true('30 minute: a run past 01:30 is rejected', st == MET_ERR_CO2_NOT_COVERED, real(st, wp))

      !----- Files that break format 1 are rejected, never guessed at. ------------------------!
      call expect_bad('no timestep line', [character(len=48) :: 'units umol/mol', '2000 300', '2001 310'])
      call expect_bad('a second timestep line', [character(len=48) :: 'timestep 1 year', 'timestep 1 year', &
                      'units umol/mol', '2000 300', '2001 310'])
      call expect_bad('an unknown unit', [character(len=48) :: 'timestep 1 week', 'units umol/mol', &
                      '2000 300', '2001 310'])
      call expect_bad('units other than umol/mol', [character(len=48) :: 'timestep 1 year', 'units ppm', &
                      '2000 300', '2001 310'])
      call expect_bad('a keyword after the data', [character(len=48) :: 'timestep 1 year', 'units umol/mol', &
                      '2000 300', '2001 310', 'units umol/mol'])
      call expect_bad('a start at the wrong precision', [character(len=48) :: 'timestep 1 year',       &
                      'units umol/mol', '2000-01 300', '2001-01 310'])
      call expect_bad('an empty date part', [character(len=48) :: 'timestep 1 month', 'units umol/mol', &
                      '2000-01 300', '2000- 310'])
      call expect_bad('a gap', [character(len=48) :: 'timestep 1 year', 'units umol/mol', '2000 300', '2002 310'])
      call expect_bad('one row', [character(len=48) :: 'timestep 1 year', 'units umol/mol', '2000 300'])
      call expect_bad('a value that is not a number', [character(len=48) :: 'timestep 1 year',         &
                      'units umol/mol', '2000 300', '2001 abc'])
      call expect_bad('a value that is not positive', [character(len=48) :: 'timestep 1 year',         &
                      'units umol/mol', '2000 300', '2001 -5'])
      call expect_bad('a third column', [character(len=48) :: 'timestep 1 year', 'units umol/mol',     &
                      '2000 300 1', '2001 310 1'])
      call open_file_co2('no_such_co2_file.txt', meds_time_t(2000_ik,1_ik,1_ik),                      &
                         meds_time_t(2001_ik,1_ik,1_ik), src, cur, st)
      call check_true('bad file: a missing file is rejected', st == MET_ERR_CO2_FILE, real(st, wp))

      !----- THE regression this exists for: recycled met repeats, the CO2 does not. Model       !
      !      2023-07-02 12:00 reads the 2021 met record but the 2023 CO2. ---------------------------!
      call write_yearfile(YF, 2021_ik)
      call write_text(CF, [character(len=48) :: 'timestep 1 year', 'units umol/mol', '2021 400',    &
                           '2022 410', '2023 420', '2024 430'])
      fc = forcing_config_t()
      fc%backend = MET_BACKEND_NETCDF ; fc%path = YF ; fc%grid_index = 1_ik
      fc%dt_forcing = 86400.0_wp ; fc%avg_convention = METAVG_END ; fc%sw_partition = SWPART_CLEARIDX
      fc%recycle = .true. ; fc%apply_solar_longitude = .false.
      fc%recycle_start = meds_time_t(2021_ik, 1_ik, 1_ik) ; fc%recycle_end = meds_time_t(2022_ik, 1_ik, 1_ik)
      fc%co2_source = CO2_SOURCE_FILE ; fc%co2_file = CF
      call met_open(src, fc, stat=st, run_start=meds_time_t(2021_ik,1_ik,1_ik),                      &
                    run_end=meds_time_t(2025_ik,1_ik,1_ik))
      call check_true('recycle: opens with a CO2 file', st == MET_OK, real(st, wp))
      call site_cursor(src, cur, fc)
      call met_advance(src, cur, meds_time_t(2021_ik,7_ik,2_ik,12_ik))
      m_ref = met_instant(src, cur, meds_time_t(2021_ik,7_ik,2_ik,12_ik))
      call met_advance(src, cur, meds_time_t(2023_ik,7_ik,2_ik,12_ik))
      m_map = met_instant(src, cur, meds_time_t(2023_ik,7_ik,2_ik,12_ik))
      call check('recycle: the met repeats (2023 reads the 2021 record)', m_map%tair_k, m_ref%tair_k, 1.0e-9_wp)
      call check('recycle: 2021 CO2', m_ref%co2, 400.0_wp, 1.0e-9_wp)
      call check('recycle: the CO2 does not repeat (2023 CO2)', m_map%co2, 420.0_wp, 1.0e-9_wp)
      call met_close(src)

      !----- A met file that carries CO2air is rejected: nothing would read it. ---------------!
      call write_synthetic_forcing(NCFILE, meds_time_t(year=2020_ik, month=7_ik, day=1_ik), with_co2air=.true.)
      fc = forcing_config_t()
      fc%backend = MET_BACKEND_NETCDF ; fc%path = NCFILE ; fc%grid_index = 1_ik
      fc%dt_forcing = 3600.0_wp ; fc%avg_convention = METAVG_END ; fc%sw_partition = SWPART_CLEARIDX
      fc%recycle = .false.
      call met_open(src, fc, stat=st)
      call check_true('CO2air: a met file carrying it is rejected', st == MET_ERR_CO2_IN_MET_FILE, real(st, wp))
      call write_synthetic_forcing(NCFILE, meds_time_t(year=2020_ik, month=7_ik, day=1_ik))

      !----- The shipped default series (data/co2/): it loads, covers 1000 .. 2023, and reads   !
      !      back the CMIP7 values at the years' middles (1700 and 1850 are not leap years). --------!
      if (len_trim(default_file) == 0) then
         print '(a)', '   NOTE: no path to the shipped CO2 file was given; its check is skipped'
      else
         call open_file_co2(default_file, meds_time_t(1000_ik,1_ik,1_ik), meds_time_t(2023_ik,1_ik,1_ik), &
                            src, cur, st)
         call check_true('shipped: loads and covers 1000-01-01 .. 2023-01-01', st == MET_OK, real(st, wp))
         if (st == MET_OK) then
            call check('shipped: 1000 (held)', co2_now(src, cur, meds_time_t(1000_ik,3_ik,1_ik)), 282.437_wp, 1.0e-9_wp)
            call check('shipped: 1700', co2_now(src, cur, meds_time_t(1700_ik,7_ik,2_ik,12_ik)), 277.537_wp, 1.0e-9_wp)
            call check('shipped: 1850', co2_now(src, cur, meds_time_t(1850_ik,7_ik,2_ik,12_ik)), 284.297_wp, 1.0e-9_wp)
            call check('shipped: 2022 (held)', co2_now(src, cur, meds_time_t(2022_ik,12_ik,31_ik)), 417.320_wp, 1.0e-9_wp)
            call met_close(src)
         end if
         call open_file_co2(default_file, meds_time_t(1000_ik,1_ik,1_ik), meds_time_t(2023_ik,1_ik,2_ik), &
                            src, cur, st)
         call check_true('shipped: a run past 2022 is rejected', st == MET_ERR_CO2_NOT_COVERED, real(st, wp))
      end if
   end subroutine test_co2_series

   !----- Open the CONST backend with co2_file = path over the run [t0, t1]; `st` is met_open's.  !
   subroutine open_file_co2(path, t0, t1, src, cur, st)
      character(len=*),   intent(in)    :: path
      type(meds_time_t),  intent(in)    :: t0, t1
      type(met_source_t), intent(inout) :: src
      type(met_cursor_t), intent(out)   :: cur
      integer(ik),        intent(out)   :: st
      type(forcing_config_t) :: fc
      fc%backend = MET_BACKEND_CONST ; fc%co2_source = CO2_SOURCE_FILE ; fc%co2_file = path
      call met_open(src, fc, stat=st, run_start=t0, run_end=t1)
      if (st == MET_OK) call site_cursor(src, cur, fc)
   end subroutine open_file_co2

   real(wp) function co2_now(src, cur, t) result(co2)
      type(met_source_t), intent(in) :: src
      type(met_cursor_t), intent(in) :: cur
      type(meds_time_t),  intent(in) :: t
      type(met_forcing_t) :: met
      met = met_instant(src, cur, t)
      co2 = met%co2
   end function co2_now

   !----- A CO2 file with `lines` must be rejected at open (MET_ERR_CO2_FILE). ----------------!
   subroutine expect_bad(what, lines)
      character(len=*), intent(in) :: what, lines(:)
      character(len=*), parameter  :: BF = 'test_met_co2_bad_tmp.txt'
      type(met_source_t) :: src
      type(met_cursor_t) :: cur
      integer(ik)        :: st
      call write_text(BF, lines)
      call open_file_co2(BF, meds_time_t(2000_ik,1_ik,1_ik), meds_time_t(2001_ik,1_ik,1_ik), src, cur, st)
      call check_true('bad file: '//what, st == MET_ERR_CO2_FILE, real(st, wp))
      if (st == MET_OK) call met_close(src)
   end subroutine expect_bad

   subroutine write_text(path, lines)
      character(len=*), intent(in) :: path, lines(:)
      integer :: u, j
      open(newunit=u, file=path, status='replace', action='write')
      do j = 1, size(lines)
         write(u, '(a)') trim(lines(j))
      end do
      close(u)
   end subroutine write_text

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
