! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_met_tower -- the contract an ED_default forcing file built from flux-tower data relies on !
! (MEDS_FLUX_TOWER_FORCING_PLAN.md §5 and §9 V6):                                            !
!                                                                                          !
!   1. the three humidity forms, each converted with the model's own curve, and the files that   !
!      carry none, two, or a percentage refused;                                                  !
!   2. the UTC requirement;                                                                       !
!   3. the heights a file states, checked against [forcing];                                      !
!   4. rain, shortwave and the clearness the longwave synthesis remembers, read from the interval   !
!      that CONTAINS the instant, on an end-stamped and a begin-stamped file alike;                !
!   5. the round trip of a tower file: the model's relative humidity at the forcing temperature is   !
!      the file's RHair, a saturated record reads back as VPD = 0, and the move from the tower      !
!      height to a canopy-air top follows the dry adiabat and the patch's log wind profile.        !
!==========================================================================================!
program test_met_tower
   use meds_test_assert,     only : check, check_true, test_report
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : grav, cp_air
   use meds_time,            only : meds_time_t, time_advance_seconds
   use meds_therm_lib,       only : sat_vapor_pressure, specific_humidity_to_vpd
   use meds_forcing_config,  only : forcing_config_t, MET_BACKEND_ED_DEFAULT, SWPART_CLEARIDX,  &
                                    METAVG_END, METAVG_BEGIN, HEIGHT_ABOVE_GROUND,              &
                                    HEIGHT_ABOVE_ZERO_PLANE, WIND_EXPOSURE_LOCAL, LW_SYNTHESIZE
   use meds_forcing_types,   only : met_forcing_t, met_source_t, met_cursor_t
   use meds_forcing_kernels, only : rh_to_specific_humidity, dewpoint_to_specific_humidity,        &
                                    clearness_index, met_solar_cosz
   use meds_lapse_rate,      only : met_to_cas_top
   use meds_met_driver,      only : met_open, met_cursor_init, met_advance, met_instant, met_close, &
                                    MET_OK, MET_ERR_NOT_UTC, MET_ERR_HUMIDITY, MET_ERR_ATTR_MISMATCH
   use meds_netcdf_c
   use iso_c_binding,        only : c_int, c_size_t, c_double
   implicit none
   character(len=*), parameter :: TOWER_FILE = 'test_met_tower_tmp.nc'
   real(wp),    parameter :: DT           = 1800.0_wp    !< [s] half-hourly, as a tower reports
   integer,     parameter :: NT           = 96           !< two days of records
   integer,     parameter :: RAIN_REC     = 30           !< the one record with rain
   real(wp),    parameter :: RAIN_RATE    = 2.0e-3_wp    !< [kg/m2/s] its rate
   integer,     parameter :: SATURATED_REC = 10          !< the record whose RHair is exactly 1
   real(wp),    parameter :: TOWER_HEIGHT = 41.0_wp      !< [m] above the ground (BCI)
   real(wp),    parameter :: PSURF        = 98800.0_wp   !< [Pa]
   type(meds_time_t)      :: base

   base = meds_time_t(year=2020_ik, month=7_ik, day=1_ik)
   call test_humidity_forms()
   call test_utc_required()
   call test_stated_heights()
   call test_flux_interval(METAVG_END)
   call test_flux_interval(METAVG_BEGIN)
   call test_tower_round_trip()
   call test_report('test_met_tower')

contains

   !----- The [forcing] block a tower file is read with: BCI's location, 41 m heights above the   !
   !      ground, a wind measured over this canopy. --------------------------------------------!
   function tower_config(convention) result(fc)
      integer(ik), intent(in) :: convention
      type(forcing_config_t)  :: fc
      fc%backend = MET_BACKEND_ED_DEFAULT ; fc%path = TOWER_FILE ; fc%grid_index = 1_ik
      fc%dt_forcing = DT ; fc%avg_convention = convention ; fc%sw_partition = SWPART_CLEARIDX
      fc%latitude_deg = 9.1568_wp ; fc%longitude_deg = -79.8486_wp ; fc%elevation_m = 150.0_wp
      fc%recycle = .false.
      fc%tq_height = TOWER_HEIGHT ; fc%wind_height = TOWER_HEIGHT
      fc%height_above = HEIGHT_ABOVE_GROUND ; fc%wind_exposure = WIND_EXPOSURE_LOCAL
   end function tower_config

   !----- The fixture's values at record k (1-based), so a test can name what it expects. ---------!
   pure real(wp) function tair_at(k) result(v)
      integer, intent(in) :: k
      v = 299.0_wp + 3.0_wp * sin(2.0_wp * 3.14159265_wp * (0.5_wp * real(k - 1, wp) - 13.0_wp) / 24.0_wp)
   end function tair_at

   pure real(wp) function rh_at(k) result(v)
      integer, intent(in) :: k
      v = 0.80_wp - 0.15_wp * sin(2.0_wp * 3.14159265_wp * (0.5_wp * real(k - 1, wp) - 13.0_wp) / 24.0_wp)
      if (k == SATURATED_REC) v = 1.0_wp
   end function rh_at

   !----- Distinct daytime means, rising record by record, so a read from the wrong interval is off  !
   !      by 10 W/m2 rather than by a rounding error. The daylit records at BCI are 11:00-23:30 UTC. -!
   pure real(wp) function sw_at(k) result(v)
      integer, intent(in) :: k
      real(wp) :: hh
      hh = modulo(0.5_wp * real(k - 1, wp), 24.0_wp)
      v = 0.0_wp
      if (hh >= 13.0_wp .and. hh <= 21.0_wp) v = 100.0_wp + 10.0_wp * real(k, wp)
   end function sw_at

   pure real(wp) function stamp_seconds(k) result(s)
      integer, intent(in) :: k
      s = real(k - 1, wp) * DT
   end function stamp_seconds

   !----- A one-point ED_default file, half-hourly from 2020-07-01 00:00 UTC. `humidity` names the  !
   !      variable(s) it carries ('RHair', 'Tdew', 'Qair', 'RHair+Qair' or 'none'); an empty         !
   !      `time_zone` or `height_above` leaves that attribute out, and so does stated_height <= 0.    !
   subroutine write_tower_file(convention, humidity, time_zone, rh_scale, stated_height, height_above)
      character(len=*), intent(in) :: convention, humidity, time_zone, height_above
      real(wp),         intent(in) :: rh_scale, stated_height
      character(len=8)  :: names(9)
      integer(c_int)    :: st, ncid, td, gd, vt, vla, vlo, vid(9), dims1(1), dims2(2)
      integer(c_size_t) :: s1(1), c1(1), s2(2), c2(2)
      real(c_double)    :: tsec(NT), dat(1, NT), one(1)
      integer :: nvar, j, k
      names(1:6) = [character(len=8) :: 'Tair', 'PSurf', 'Wind', 'Rainf', 'SWdown', 'LWdown']
      nvar = 6
      select case (humidity)
      case ('RHair')      ; nvar = 7 ; names(7) = 'RHair'
      case ('Tdew')       ; nvar = 7 ; names(7) = 'Tdew'
      case ('Qair')       ; nvar = 7 ; names(7) = 'Qair'
      case ('RHair+Qair') ; nvar = 8 ; names(7) = 'RHair' ; names(8) = 'Qair'
      end select
      st = nc_create_f(TOWER_FILE, NC_NETCDF4, ncid) ; call nc_check(st, 'tower: create')
      st = nc_def_dim_f(ncid, 'time', int(NT, c_size_t), td) ; call nc_check(st, 'tower: time dim')
      st = nc_def_dim_f(ncid, 'grid', 1_c_size_t, gd) ; call nc_check(st, 'tower: grid dim')
      dims1(1) = td
      st = nc_def_var_f(ncid, 'time', NC_DOUBLE, 1, dims1, vt) ; call nc_check(st, 'tower: time var')
      st = nc_put_att_text_f(ncid, vt, 'units', 33_c_size_t, 'seconds since 2020-07-01 00:00:00')
      call nc_check(st, 'tower: units')
      dims1(1) = gd
      st = nc_def_var_f(ncid, 'latitude',  NC_DOUBLE, 1, dims1, vla) ; call nc_check(st, 'tower: lat')
      st = nc_def_var_f(ncid, 'longitude', NC_DOUBLE, 1, dims1, vlo) ; call nc_check(st, 'tower: lon')
      dims2(1) = td ; dims2(2) = gd
      do j = 1, nvar
         st = nc_def_var_f(ncid, trim(names(j)), NC_DOUBLE, 2, dims2, vid(j)) ; call nc_check(st, 'tower: var')
      end do
      !----- The global attributes make_tower_forcing.py writes; the fixture must carry them or the  !
      !      checks under test would be skipped here and fire in production. ---------------------!
      st = nc_put_att_text_f(ncid, NC_GLOBAL, 'avg_convention', int(len(convention), c_size_t), convention)
      call nc_check(st, 'tower: avg_convention')
      st = nc_put_att_text_f(ncid, NC_GLOBAL, 'sw_input_kind', 5_c_size_t, 'total')
      call nc_check(st, 'tower: sw_input_kind')
      if (len(time_zone) > 0) then
         st = nc_put_att_text_f(ncid, NC_GLOBAL, 'time_zone', int(len(time_zone), c_size_t), time_zone)
         call nc_check(st, 'tower: time_zone')
      end if
      if (stated_height > 0.0_wp) then
         st = nc_put_att_double_f(ncid, NC_GLOBAL, 'tq_height_m', NC_DOUBLE, real(stated_height, c_double))
         call nc_check(st, 'tower: tq_height_m')
         st = nc_put_att_double_f(ncid, NC_GLOBAL, 'wind_height_m', NC_DOUBLE, real(stated_height, c_double))
         call nc_check(st, 'tower: wind_height_m')
      end if
      if (len(height_above) > 0) then
         st = nc_put_att_text_f(ncid, NC_GLOBAL, 'height_above', int(len(height_above), c_size_t), height_above)
         call nc_check(st, 'tower: height_above')
      end if
      st = nc_enddef(ncid) ; call nc_check(st, 'tower: enddef')

      do k = 1, NT ; tsec(k) = real(stamp_seconds(k), c_double) ; end do
      s1(1) = 0_c_size_t ; c1(1) = int(NT, c_size_t)
      st = nc_put_vara_double(ncid, vt, s1, c1, tsec) ; call nc_check(st, 'tower: time vals')
      c1(1) = 1_c_size_t
      one = 9.1568_c_double  ; st = nc_put_vara_double(ncid, vla, s1, c1, one) ; call nc_check(st, 'tower: lat vals')
      one = -79.8486_c_double ; st = nc_put_vara_double(ncid, vlo, s1, c1, one) ; call nc_check(st, 'tower: lon vals')
      s2 = 0_c_size_t ; c2 = [int(NT, c_size_t), 1_c_size_t]
      do j = 1, nvar
         do k = 1, NT
            select case (trim(names(j)))
            case ('Tair')   ; dat(1, k) = tair_at(k)
            case ('PSurf')  ; dat(1, k) = PSURF
            case ('Wind')   ; dat(1, k) = 2.5_wp
            case ('Rainf')  ; dat(1, k) = merge(RAIN_RATE, 0.0_wp, k == RAIN_REC)
            case ('SWdown') ; dat(1, k) = sw_at(k)
            case ('LWdown') ; dat(1, k) = 450.0_wp
            case ('RHair')  ; dat(1, k) = rh_scale * rh_at(k)
            case ('Tdew')   ; dat(1, k) = tair_at(k) - 3.0_wp
            case ('Qair')   ; dat(1, k) = 0.018_wp
            end select
         end do
         st = nc_put_vara_double(ncid, vid(j), s2, c2, dat) ; call nc_check(st, 'tower: vals')
      end do
      st = nc_close(ncid) ; call nc_check(st, 'tower: close')
   end subroutine write_tower_file

   !----- The forcing at an instant `sec` seconds after the file's first stamp. --------------------!
   function sample(src, cur, sec) result(met)
      type(met_source_t), intent(in)    :: src
      type(met_cursor_t), intent(inout) :: cur
      real(wp),           intent(in)    :: sec
      type(met_forcing_t) :: met
      type(meds_time_t)   :: now
      now = time_advance_seconds(base, sec)
      call met_advance(src, cur, now)
      met = met_instant(src, cur, now)
   end function sample

   !----- 1. Each humidity form becomes q through the model's own conversion, at the record's own   !
   !      temperature and pressure; a file with none, with two, or with RHair in percent is refused. -!
   subroutine test_humidity_forms()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met
      integer(ik) :: st
      integer, parameter :: K = 40
      print '(a)', '-- test 1: humidity forms --'
      fc = tower_config(METAVG_BEGIN)

      call write_tower_file('begin', 'RHair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('an RHair file opens', st == MET_OK, real(st, wp))
      call met_cursor_init(src, cur, 1_ik, fc%latitude_deg, fc%longitude_deg, fc%elevation_m)
      met = sample(src, cur, stamp_seconds(K))
      call check('RHair -> q by rh_to_specific_humidity at the record''s T and P', met%qair,           &
                 rh_to_specific_humidity(rh_at(K), tair_at(K), PSURF), 1.0e-15_wp)
      call met_close(src)

      call write_tower_file('begin', 'Tdew', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a Tdew file opens', st == MET_OK, real(st, wp))
      call met_cursor_init(src, cur, 1_ik, fc%latitude_deg, fc%longitude_deg, fc%elevation_m)
      met = sample(src, cur, stamp_seconds(K))
      call check('Tdew -> q by dewpoint_to_specific_humidity', met%qair,                              &
                 dewpoint_to_specific_humidity(tair_at(K) - 3.0_wp, PSURF), 1.0e-15_wp)
      call met_close(src)

      call write_tower_file('begin', 'Qair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a Qair file opens', st == MET_OK, real(st, wp))
      call met_cursor_init(src, cur, 1_ik, fc%latitude_deg, fc%longitude_deg, fc%elevation_m)
      met = sample(src, cur, stamp_seconds(K))
      call check('Qair is used as it is', met%qair, 0.018_wp, 1.0e-15_wp)
      call met_close(src)

      call write_tower_file('begin', 'RHair+Qair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a file with two humidity variables is refused', st == MET_ERR_HUMIDITY, real(st, wp))
      call met_close(src)

      call write_tower_file('begin', 'none', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a file with no humidity variable is refused', st == MET_ERR_HUMIDITY, real(st, wp))
      call met_close(src)

      call write_tower_file('begin', 'RHair', 'UTC', 100.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('RHair written in percent is refused', st == MET_ERR_HUMIDITY, real(st, wp))
      call met_close(src)
   end subroutine test_humidity_forms

   !----- 2. MEDS runs in UTC: a file that does not say time_zone = "UTC" is refused. -----------!
   subroutine test_utc_required()
      type(met_source_t)     :: src
      type(forcing_config_t) :: fc
      integer(ik) :: st
      print '(a)', '-- test 2: the UTC requirement --'
      fc = tower_config(METAVG_BEGIN)
      call write_tower_file('begin', 'RHair', '', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a file without time_zone is refused', st == MET_ERR_NOT_UTC, real(st, wp))
      call met_close(src)
      call write_tower_file('begin', 'RHair', 'UTC-05:00', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a local-time file is refused', st == MET_ERR_NOT_UTC, real(st, wp))
      call met_close(src)
      call write_tower_file('begin', 'RHair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      call met_open(src, fc, stat=st)
      call check_true('a UTC file opens', st == MET_OK, real(st, wp))
      call met_close(src)
   end subroutine test_utc_required

   !----- 3. The heights a file states must be the ones [forcing] moves the forcing from. ---------!
   subroutine test_stated_heights()
      type(met_source_t)     :: src
      type(forcing_config_t) :: fc
      integer(ik) :: st
      print '(a)', '-- test 3: the stated heights --'
      call write_tower_file('begin', 'RHair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      fc = tower_config(METAVG_BEGIN)
      fc%tq_height = 2.0_wp
      call met_open(src, fc, stat=st)
      call check_true('tq_height against a file stating 41 m is refused', st == MET_ERR_ATTR_MISMATCH, real(st, wp))
      call met_close(src)
      fc = tower_config(METAVG_BEGIN)
      fc%wind_height = 10.0_wp
      call met_open(src, fc, stat=st)
      call check_true('wind_height against a file stating 41 m is refused', st == MET_ERR_ATTR_MISMATCH, real(st, wp))
      call met_close(src)
      fc = tower_config(METAVG_BEGIN)
      fc%height_above = HEIGHT_ABOVE_ZERO_PLANE
      call met_open(src, fc, stat=st)
      call check_true('height_above against a file stating "ground" is refused',                     &
                      st == MET_ERR_ATTR_MISMATCH, real(st, wp))
      call met_close(src)
      call write_tower_file('begin', 'RHair', 'UTC', 1.0_wp, 0.0_wp, '')
      fc = tower_config(METAVG_BEGIN)
      fc%tq_height = 2.0_wp
      call met_open(src, fc, stat=st)
      call check_true('a file that states no heights is not checked', st == MET_OK, real(st, wp))
      call met_close(src)
   end subroutine test_stated_heights

   !----- 4. Rain and shortwave are interval means, read from the record whose interval CONTAINS   !
   !      the instant: record k covers [t_k, t_k+dt) on a begin-stamped file and (t_k-dt, t_k] on an !
   !      end-stamped one. A read one record off puts every storm in the wrong half-hour; the       !
   !      shortwave check averages the instantaneous flux over the reconstruction's own ten         !
   !      midpoints, which returns the interval mean exactly when the right record is used.         !
   subroutine test_flux_interval(convention)
      integer(ik), intent(in) :: convention
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met
      real(wp) :: t0, sw_mean
      integer  :: i, k_sw
      integer(ik) :: st
      if (convention == METAVG_BEGIN) then
         print '(a)', '-- test 4b: rain and shortwave intervals, begin-stamped --'
         call write_tower_file('begin', 'RHair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
         t0 = stamp_seconds(RAIN_REC)                        ! the rain interval starts at its stamp
      else
         print '(a)', '-- test 4a: rain and shortwave intervals, end-stamped --'
         call write_tower_file('end', 'RHair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
         t0 = stamp_seconds(RAIN_REC) - DT                   ! ... ends at its stamp
      end if
      fc = tower_config(convention)
      call met_open(src, fc, stat=st)
      call check_true('the file opens', st == MET_OK, real(st, wp))
      call met_cursor_init(src, cur, 1_ik, fc%latitude_deg, fc%longitude_deg, fc%elevation_m)
      met = sample(src, cur, t0 - 0.5_wp * DT)
      call check('no rain in the interval before', met%rainf + met%snowfall, 0.0_wp, 1.0e-15_wp)
      met = sample(src, cur, t0 + 0.25_wp * DT)
      call check('the rain falls in its own interval (early)', met%rainf + met%snowfall, RAIN_RATE, 1.0e-15_wp)
      met = sample(src, cur, t0 + 0.75_wp * DT)
      call check('the rain falls in its own interval (late)', met%rainf + met%snowfall, RAIN_RATE, 1.0e-15_wp)
      met = sample(src, cur, t0 + 1.5_wp * DT)
      call check('no rain in the interval after', met%rainf + met%snowfall, 0.0_wp, 1.0e-15_wp)

      !----- A midday interval, 17:00-17:30 UTC (about solar noon at BCI). -------------------------!
      t0 = 17.0_wp * 3600.0_wp
      k_sw = nint(t0 / DT) + 1                               ! begin: the record stamped 17:00
      if (convention == METAVG_END) k_sw = k_sw + 1          ! end:   the record stamped 17:30
      sw_mean = 0.0_wp
      do i = 1, 10
         met = sample(src, cur, t0 + (real(i, wp) - 0.5_wp) * DT / 10.0_wp)
         sw_mean = sw_mean + (met%par_beam + met%par_diffuse + met%nir_beam + met%nir_diffuse) / 10.0_wp
      end do
      call check('the shortwave of the interval is its own record''s mean', sw_mean, sw_at(k_sw), 1.0e-9_wp)
      call met_close(src)

      !----- The clearness index the longwave synthesis holds through the night is the containing   !
      !      interval's too: on a begin file, the next record's would be the dark one after sunset.  !
      fc%lwdown_source = LW_SYNTHESIZE
      call met_open(src, fc, stat=st)
      call met_cursor_init(src, cur, 1_ik, fc%latitude_deg, fc%longitude_deg, fc%elevation_m)
      met = sample(src, cur, t0 + 0.4_wp * DT)
      call check('the remembered clearness is the containing interval''s', cur%kt_last_day,              &
                 clearness_index(sw_at(k_sw), met_solar_cosz(time_advance_seconds(base, t0 + 0.4_wp * DT),  &
                                 t0 + 0.4_wp * DT, fc%latitude_deg, fc%longitude_deg)), 1.0e-12_wp)
      call met_close(src)
   end subroutine test_flux_interval

   !----- 5. The round trip of a tower file (V6). ------------------------------------------------!
   subroutine test_tower_round_trip()
      type(met_source_t)     :: src
      type(met_cursor_t)     :: cur
      type(forcing_config_t) :: fc
      type(met_forcing_t)    :: met, top
      real(wp) :: vpd, rh_model, worst, z_top, displace, rough, factor
      integer  :: k
      integer(ik) :: st
      print '(a)', '-- test 5: the round trip of a tower file --'
      call write_tower_file('begin', 'RHair', 'UTC', 1.0_wp, TOWER_HEIGHT, 'ground')
      fc = tower_config(METAVG_BEGIN)
      call met_open(src, fc, stat=st)
      call check_true('the tower file opens', st == MET_OK, real(st, wp))
      call met_cursor_init(src, cur, 1_ik, fc%latitude_deg, fc%longitude_deg, fc%elevation_m)

      !----- At every stamp the model's own inverse returns the tower's relative humidity. This is  !
      !      what storing RHair buys: the provider's saturation curve never enters. ----------------!
      worst = 0.0_wp
      do k = 1, NT - 1
         if (k == SATURATED_REC) cycle
         met = sample(src, cur, stamp_seconds(k))
         vpd = specific_humidity_to_vpd(met%tair_k, met%qair, met%psurf_pa)
         rh_model = 1.0_wp - vpd / sat_vapor_pressure(met%tair_k)
         worst = max(worst, abs(rh_model - rh_at(k)))
      end do
      call check('the model''s RH at the forcing temperature is RHair at every stamp', worst, 0.0_wp, 1.0e-12_wp)
      met = sample(src, cur, stamp_seconds(SATURATED_REC))
      vpd = specific_humidity_to_vpd(met%tair_k, met%qair, met%psurf_pa)
      call check('a saturated tower record reads back as VPD = 0', vpd, 0.0_wp, 1.0e-9_wp)

      !----- The move to the top of a 25 m stand's canopy air (30 m) from 41 m above the ground. ----!
      met = sample(src, cur, stamp_seconds(40))
      z_top = 30.0_wp ; displace = 0.63_wp * 25.0_wp ; rough = 0.13_wp * 25.0_wp
      top = met_to_cas_top(met, fc, z_top, displace, rough)
      call check('the canopy-air top is warmer along the dry adiabat', top%tair_k - met%tair_k,         &
                 -(grav / cp_air) * (z_top - TOWER_HEIGHT), 1.0e-12_wp)
      call check('specific humidity is carried unchanged', top%qair, met%qair, 1.0e-15_wp)
      factor = log((z_top - displace) / rough) / log((TOWER_HEIGHT - displace) / rough)
      call check('the wind follows the patch log profile from the tower height', top%wind,           &
                 factor * met%wind, 1.0e-12_wp)
      call met_close(src)
   end subroutine test_tower_round_trip

end program test_met_tower
