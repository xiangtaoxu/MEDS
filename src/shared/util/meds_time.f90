!==========================================================================================!
! meds_time -- calendar date tracking with real (leap-year-aware) Gregorian dates.         !
!                                                                                          !
! MEDS advances on a real calendar: the run is bounded by a start and end DATE (not a count !
! of "years"), and both output streams stamp the simulated calendar date. This module owns   !
! the date arithmetic so the driver/IO never re-derive leap years by hand.                    !
!                                                                                          !
! The atom is meds_time_t (year/month/day + hour/minute/second). All date stepping goes        !
! through the proleptic-Gregorian Julian Day Number (to_jdn/from_jdn), which makes leap years  !
! exact and month/year roll-over trivial; the JDN arithmetic uses 64-bit (lk) intermediates    !
! because the inverse formula multiplies by ~4000 (overflows int32 for JDN ~ 2.45e6). Helpers: !
! parse/format ISO-ish strings, a 14-char YYYYMMDDHHMMSS stamp for filenames, decimal calendar !
! year for plotting, and lexicographic comparisons for the run-loop bound.                      !
!==========================================================================================!
module meds_time
   use meds_kinds,     only : wp, ik, lk
   use meds_constants, only : pi
   implicit none
   private

   public :: meds_time_t
   public :: is_leap_year, days_in_month, days_in_year, day_of_year
   public :: time_advance_days, time_advance_months, days_between, years_between
   public :: seconds_into_day, seconds_between, time_advance_seconds
   public :: time_lt, time_le, time_eq, time_valid
   public :: time_from_string, time_to_string, time_to_stamp, time_to_decimal_year
   public :: solar_cosz, daylength, doy_effective

   !----- A calendar instant. Defaults give a valid date so meds_time_t() is usable. ------!
   type :: meds_time_t
      integer(ik) :: year   = 2000_ik
      integer(ik) :: month  = 1_ik
      integer(ik) :: day    = 1_ik
      integer(ik) :: hour   = 0_ik
      integer(ik) :: minute = 0_ik
      integer(ik) :: second = 0_ik
   end type meds_time_t

   integer(ik), parameter :: month_len(12) = [31_ik, 28_ik, 31_ik, 30_ik, 31_ik, 30_ik,       &
                                               31_ik, 31_ik, 30_ik, 31_ik, 30_ik, 31_ik]

contains

   !---------------------------------------------------------------------------------------!
   ! Proleptic-Gregorian leap-year rule.                                                     !
   !---------------------------------------------------------------------------------------!
   elemental logical function is_leap_year(year) result(leap)
      integer(ik), intent(in) :: year
      leap = (mod(year, 4_ik) == 0_ik .and. mod(year, 100_ik) /= 0_ik) .or.                   &
             mod(year, 400_ik) == 0_ik
   end function is_leap_year

   elemental integer(ik) function days_in_month(year, month) result(d)
      integer(ik), intent(in) :: year, month
      d = month_len(month)
      if (month == 2_ik .and. is_leap_year(year)) d = 29_ik
   end function days_in_month

   elemental integer(ik) function days_in_year(year) result(d)
      integer(ik), intent(in) :: year
      d = 365_ik
      if (is_leap_year(year)) d = 366_ik
   end function days_in_year

   !----- 1-based day of the year (1 = Jan 1; 365 or 366 = Dec 31). ------------------------!
   pure integer(ik) function day_of_year(t) result(doy)
      type(meds_time_t), intent(in) :: t
      type(meds_time_t) :: jan1
      jan1 = meds_time_t(year = t%year)        ! Jan 1, 00:00 of the same year
      doy  = int(to_jdn(t) - to_jdn(jan1), ik) + 1_ik
   end function day_of_year

   !---------------------------------------------------------------------------------------!
   ! Cosine of the solar zenith angle from the TIME DIMENSION (date -> declination, sub-daily !
   ! cursor -> hour angle) + site latitude. cosz is a DERIVED quantity, not a met forcing:     !
   ! the fast loop calls this each dt_fast with `t_sec` = local-solar seconds into the day.     !
   ! Cooper (1969) declination; local apparent solar time (no longitude / equation-of-time --  !
   ! adequate for a single-site diurnal cycle). Floored at 0 (night).                           !
   !---------------------------------------------------------------------------------------!
   pure real(wp) function solar_cosz(t, t_sec, latitude_deg) result(cosz)
      type(meds_time_t), intent(in) :: t
      real(wp),          intent(in) :: t_sec          !< [s] local solar seconds into the day (0 = midnight)
      real(wp),          intent(in) :: latitude_deg   !< [deg] site latitude (+N)
      real(wp) :: deg2rad, lat, decl, hourangle, frac_day
      deg2rad   = pi / 180.0_wp
      lat       = latitude_deg * deg2rad
      decl      = 23.45_wp * deg2rad * sin(2.0_wp * pi * (284.0_wp + real(day_of_year(t), wp)) / 365.0_wp)
      frac_day  = t_sec / 86400.0_wp                  ! 0 at midnight, 0.5 at solar noon
      hourangle = 2.0_wp * pi * (frac_day - 0.5_wp)   ! -pi (midnight) .. 0 (noon) .. +pi
      cosz      = sin(lat) * sin(decl) + cos(lat) * cos(decl) * cos(hourangle)
      cosz      = max(cosz, 0.0_wp)
   end function solar_cosz

   !---------------------------------------------------------------------------------------!
   ! Daylength [h] from latitude [deg] + day-of-year (White et al. 1997 form). Relocated from  !
   ! meds_phenology so any seasonal process can call it without a plant-library edge. The ED2   !
   ! polar branch is FIXED: polar DAY is arg <= -1 (ED2 wrote arg <= 1, a bug). NOTE: this uses  !
   ! the White-1997 declination -23.44*cos(2pi(doy+9)/365), DISTINCT from solar_cosz's Cooper    !
   ! 23.45*sin(2pi(284+doy)/365) form -- the two standard approximations are kept separate on     !
   ! purpose (unifying would re-baseline the daylength golden values; see the phenology design).  !
   !---------------------------------------------------------------------------------------!
   elemental real(wp) function daylength(lat_deg, doy) result(dl)
      real(wp),    intent(in) :: lat_deg
      integer(ik), intent(in) :: doy
      real(wp) :: latr, decl, arg
      latr = lat_deg * pi / 180.0_wp
      decl = -23.44_wp * pi / 180.0_wp * cos(2.0_wp * pi / 365.0_wp * (real(doy, wp) + 9.0_wp))
      arg  = -tan(latr) * tan(decl)
      if (arg >= 1.0_wp) then
         dl = 0.0_wp                 ! polar night
      else if (arg <= -1.0_wp) then
         dl = 24.0_wp                ! polar day (FIX: ED2's branch read arg <= 1)
      else
         dl = 24.0_wp / pi * acos(arg)
      end if
   end function daylength

   !----- Effective (northern-hemisphere-equivalent) day-of-year: shift SH by half a year, so   !
   !      season-gated accumulators (GDD/chilling) can share one hemisphere-agnostic calendar. --!
   elemental integer(ik) function doy_effective(doy, hemis_north) result(de)
      integer(ik), intent(in) :: doy
      logical,     intent(in) :: hemis_north
      if (hemis_north) then
         de = doy
      else
         de = modulo(doy - 1_ik + 182_ik, 365_ik) + 1_ik
      end if
   end function doy_effective

   !=======================================================================================!
   ! Julian Day Number conversions (Fliegel & Van Flandern; 64-bit intermediates).          !
   !=======================================================================================!
   pure integer(lk) function to_jdn(t) result(jdn)
      type(meds_time_t), intent(in) :: t
      integer(lk) :: y, m, d
      y = int(t%year, lk) ; m = int(t%month, lk) ; d = int(t%day, lk)
      jdn = (1461_lk * (y + 4800_lk + (m - 14_lk) / 12_lk)) / 4_lk                            &
          + (367_lk * (m - 2_lk - 12_lk * ((m - 14_lk) / 12_lk))) / 12_lk                     &
          - (3_lk * ((y + 4900_lk + (m - 14_lk) / 12_lk) / 100_lk)) / 4_lk                    &
          + d - 32075_lk
   end function to_jdn

   pure subroutine from_jdn(jdn, year, month, day)
      integer(lk), intent(in)  :: jdn
      integer(ik), intent(out) :: year, month, day
      integer(lk) :: l, n, i, j, k
      l = jdn + 68569_lk
      n = (4_lk * l) / 146097_lk
      l = l - (146097_lk * n + 3_lk) / 4_lk
      i = (4000_lk * (l + 1_lk)) / 1461001_lk
      l = l - (1461_lk * i) / 4_lk + 31_lk
      j = (80_lk * l) / 2447_lk
      k = l - (2447_lk * j) / 80_lk
      l = j / 11_lk
      j = j + 2_lk - 12_lk * l
      i = 100_lk * (n - 49_lk) + i + l
      year = int(i, ik) ; month = int(j, ik) ; day = int(k, ik)
   end subroutine from_jdn

   !=======================================================================================!
   ! Date stepping.                                                                         !
   !=======================================================================================!
   !----- Advance by a whole number of days (leap years exact); time-of-day unchanged. -----!
   pure function time_advance_days(t, ndays) result(t2)
      type(meds_time_t), intent(in) :: t
      integer(ik),       intent(in) :: ndays
      type(meds_time_t) :: t2
      t2 = t
      call from_jdn(to_jdn(t) + int(ndays, lk), t2%year, t2%month, t2%day)
   end function time_advance_days

   !----- Advance by whole months; the day is clamped into the target month (e.g. Jan 31 ---!
   !      + 1 month -> Feb 28/29), as in ED2's monthly stepping.                             !
   pure function time_advance_months(t, nmonths) result(t2)
      type(meds_time_t), intent(in) :: t
      integer(ik),       intent(in) :: nmonths
      type(meds_time_t) :: t2
      integer(ik) :: total
      t2 = t
      total    = t%year * 12_ik + (t%month - 1_ik) + nmonths
      t2%year  = total / 12_ik
      t2%month = mod(total, 12_ik) + 1_ik
      t2%day   = min(t%day, days_in_month(t2%year, t2%month))
   end function time_advance_months

   !----- Whole days from a to b (positive if b is later). ---------------------------------!
   pure integer(ik) function days_between(a, b) result(d)
      type(meds_time_t), intent(in) :: a, b
      d = int(to_jdn(b) - to_jdn(a), ik)
   end function days_between

   !----- Difference in decimal calendar years (b - a); for progress reporting/labels. -----!
   pure real(wp) function years_between(a, b) result(y)
      type(meds_time_t), intent(in) :: a, b
      y = time_to_decimal_year(b) - time_to_decimal_year(a)
   end function years_between

   !=======================================================================================!
   ! Comparisons (lexicographic on the full instant).                                       !
   !=======================================================================================!
   pure integer(lk) function sec_of_day(t) result(s)
      type(meds_time_t), intent(in) :: t
      s = int(t%hour, lk) * 3600_lk + int(t%minute, lk) * 60_lk + int(t%second, lk)
   end function sec_of_day

   !=======================================================================================!
   ! Second-level helpers (the met-forcing interpolation / window slide needs sub-daily,      !
   ! second-accurate arithmetic across day boundaries -- meds_forcing / meds_met_driver).      !
   !=======================================================================================!
   !----- Seconds elapsed since midnight of the instant's own day (0 .. 86399). -------------!
   pure real(wp) function seconds_into_day(t) result(s)
      type(meds_time_t), intent(in) :: t
      s = real(sec_of_day(t), wp)
   end function seconds_into_day

   !----- Total seconds from a to b (positive if b is later); exact across days via the JDN. -!
   pure real(wp) function seconds_between(a, b) result(s)
      type(meds_time_t), intent(in) :: a, b
      s = real((to_jdn(b) - to_jdn(a)) * 86400_lk + sec_of_day(b) - sec_of_day(a), wp)
   end function seconds_between

   !----- Advance by dsec seconds (may be negative or many days); rolls date + time-of-day. --!
   pure function time_advance_seconds(t, dsec) result(t2)
      type(meds_time_t), intent(in) :: t
      real(wp),          intent(in) :: dsec
      type(meds_time_t) :: t2
      integer(lk) :: total, ndays, rem
      total = sec_of_day(t) + nint(dsec, lk)             ! seconds from midnight of t's day
      ndays = int(floor(real(total, wp) / 86400.0_wp), lk)
      rem   = total - ndays * 86400_lk                   ! guaranteed 0 .. 86399
      t2 = time_advance_days(t, int(ndays, ik))
      t2%hour   = int(rem / 3600_lk, ik)
      t2%minute = int(mod(rem, 3600_lk) / 60_lk, ik)
      t2%second = int(mod(rem, 60_lk), ik)
   end function time_advance_seconds

   pure logical function time_lt(a, b) result(yes)
      type(meds_time_t), intent(in) :: a, b
      integer(lk) :: ja, jb
      ja = to_jdn(a) ; jb = to_jdn(b)
      if (ja /= jb) then
         yes = ja < jb
      else
         yes = sec_of_day(a) < sec_of_day(b)
      end if
   end function time_lt

   pure logical function time_eq(a, b) result(yes)
      type(meds_time_t), intent(in) :: a, b
      yes = to_jdn(a) == to_jdn(b) .and. sec_of_day(a) == sec_of_day(b)
   end function time_eq

   pure logical function time_le(a, b) result(yes)
      type(meds_time_t), intent(in) :: a, b
      yes = time_lt(a, b) .or. time_eq(a, b)
   end function time_le

   !----- Range check on the fields (used to validate parsed configuration). ---------------!
   pure logical function time_valid(t) result(ok)
      type(meds_time_t), intent(in) :: t
      ok = t%month >= 1_ik  .and. t%month <= 12_ik                                           &
           .and. t%day  >= 1_ik  .and. t%day   <= days_in_month(t%year, max(min(t%month,12_ik),1_ik)) &
           .and. t%hour   >= 0_ik .and. t%hour   <= 23_ik                                    &
           .and. t%minute >= 0_ik .and. t%minute <= 59_ik                                    &
           .and. t%second >= 0_ik .and. t%second <= 59_ik
   end function time_valid

   !=======================================================================================!
   ! String <-> time.                                                                       !
   !=======================================================================================!
   !----- Parse "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS" (also accepts '/', 'T' separators).    !
   !      Missing trailing fields default (month/day -> 1, time -> 0). ok reports a clean,    !
   !      in-range parse so the caller can fall back to the default date.                    !
   subroutine time_from_string(s, t, ok)
      character(len=*),  intent(in)  :: s
      type(meds_time_t), intent(out) :: t
      logical,           intent(out) :: ok
      character(len=len(s)) :: buf
      integer(ik) :: i, ntok, ios
      integer(ik) :: y, mo, d, h, mi, se

      t = meds_time_t()
      y = 2000_ik ; mo = 1_ik ; d = 1_ik ; h = 0_ik ; mi = 0_ik ; se = 0_ik

      !----- Turn every separator into a blank, then read with list-directed input. --------!
      buf = s
      do i = 1_ik, int(len_trim(buf), ik)
         select case (buf(i:i))
         case ('-', '/', ':', 'T', 't', ',') ; buf(i:i) = ' '
         end select
      end do

      ntok = count_tokens(buf)
      ios  = 0_ik
      select case (ntok)
      case (6_ik:) ; read(buf, *, iostat=ios) y, mo, d, h, mi, se
      case (5_ik)  ; read(buf, *, iostat=ios) y, mo, d, h, mi
      case (4_ik)  ; read(buf, *, iostat=ios) y, mo, d, h
      case (3_ik)  ; read(buf, *, iostat=ios) y, mo, d
      case (2_ik)  ; read(buf, *, iostat=ios) y, mo
      case (1_ik)  ; read(buf, *, iostat=ios) y
      case default ; ios = -1_ik
      end select

      t = meds_time_t(year=y, month=mo, day=d, hour=h, minute=mi, second=se)
      ok = (ios == 0_ik) .and. time_valid(t)
   end subroutine time_from_string

   !----- "YYYY-MM-DD HH:MM:SS" (human-readable, e.g. for attributes/log lines). -----------!
   pure function time_to_string(t) result(s)
      type(meds_time_t), intent(in) :: t
      character(len=19) :: s
      write(s, '(i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)')                         &
            t%year, t%month, t%day, t%hour, t%minute, t%second
   end function time_to_string

   !----- "YYYYMMDDHHMMSS" (14 chars; the state-file timestamp <prefix>-S-<stamp>.nc). ------!
   pure function time_to_stamp(t) result(s)
      type(meds_time_t), intent(in) :: t
      character(len=14) :: s
      write(s, '(i4.4,5i2.2)') t%year, t%month, t%day, t%hour, t%minute, t%second
   end function time_to_stamp

   !----- Decimal calendar year (year + elapsed fraction); the continuous plotting axis. ---!
   pure real(wp) function time_to_decimal_year(t) result(yr)
      type(meds_time_t), intent(in) :: t
      real(wp) :: frac
      frac = (real(day_of_year(t) - 1_ik, wp)                                                &
              + real(sec_of_day(t), wp) / 86400.0_wp) / real(days_in_year(t%year), wp)
      yr = real(t%year, wp) + frac
   end function time_to_decimal_year

   !----- Count whitespace-separated tokens. -----------------------------------------------!
   pure integer(ik) function count_tokens(s) result(nt)
      character(len=*), intent(in) :: s
      integer(ik) :: i
      logical     :: in_tok
      nt = 0_ik ; in_tok = .false.
      do i = 1_ik, int(len_trim(s), ik)
         if (s(i:i) /= ' ') then
            if (.not. in_tok) nt = nt + 1_ik
            in_tok = .true.
         else
            in_tok = .false.
         end if
      end do
   end function count_tokens

end module meds_time
