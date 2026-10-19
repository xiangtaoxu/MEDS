! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_time -- the calendar. A run goes from a start date to an end date on the proleptic  !
! Gregorian calendar, leap years included, and every output record carries its date. This  !
! module does the date arithmetic, through the Julian Day Number, and the solar geometry   !
! that follows from the date: the declination, the zenith cosine and the daylength         !
! (docs/science/forcing.md §4).                                                            !
!==========================================================================================!
module meds_time
   use meds_kinds,     only : wp, ik, lk
   use meds_constants, only : pi
   implicit none
   private

   public :: meds_time_t
   public :: is_leap_year, days_in_month, days_in_year, day_of_year
   public :: time_advance_days, time_advance_months, days_between, years_between
   public :: time_advance_years, whole_years_between
   public :: seconds_into_day, seconds_between, time_advance_seconds
   public :: time_lt, time_le, time_eq, time_valid
   public :: time_from_string, time_units_base, time_to_string, time_to_stamp, time_to_decimal_year
   public :: solar_cosz, daylength, doy_effective, solar_declination

   !----- A calendar instant; the defaults make meds_time_t() a valid date. ---------------!
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

   !----- The proleptic-Gregorian leap-year rule. -----------------------------------------!
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

   !----- The day of the year, 1 on 1 January. --------------------------------------------!
   pure integer(ik) function day_of_year(t) result(doy)
      type(meds_time_t), intent(in) :: t
      type(meds_time_t) :: jan1
      jan1 = meds_time_t(year = t%year)        ! 1 January, 00:00, of the same year
      doy  = int(to_jdn(t) - to_jdn(jan1), ik) + 1_ik
   end function day_of_year

   !----- The solar declination [rad], Cooper (1969): 23.45 deg sin(2 pi (284 + doy)/365)  !
   !      (docs/science/forcing.md eq. (4)). It is the model's only declination:           !
   !      solar_cosz and daylength both use it. -------------------------------------------!
   elemental real(wp) function solar_declination(doy) result(decl)
      integer(ik), intent(in) :: doy          !< day of year (1 = 1 January)
      decl = 23.45_wp * pi / 180.0_wp * sin(2.0_wp * pi * (284.0_wp + real(doy, wp)) / 365.0_wp)
   end function solar_declination

   !----- The cosine of the solar zenith angle, 0 at night (docs/science/forcing.md        !
   !      eq. (4)). t_sec is local apparent solar time, which the forcing makes from UTC,  !
   !      the longitude and the equation of time (eqs. (2)-(3)). --------------------------!
   pure real(wp) function solar_cosz(t, t_sec, latitude_deg) result(cosz)
      type(meds_time_t), intent(in) :: t
      real(wp),          intent(in) :: t_sec          !< [s] local solar seconds into the day (0 = midnight)
      real(wp),          intent(in) :: latitude_deg   !< [deg] site latitude (+N)
      real(wp) :: deg2rad, lat, decl, hourangle, frac_day
      deg2rad   = pi / 180.0_wp
      lat       = latitude_deg * deg2rad
      decl      = solar_declination(day_of_year(t))
      frac_day  = t_sec / 86400.0_wp                  ! 0 at midnight, 0.5 at solar noon
      hourangle = 2.0_wp * pi * (frac_day - 0.5_wp)   ! -pi (midnight) .. 0 (noon) .. +pi
      cosz      = sin(lat) * sin(decl) + cos(lat) * cos(decl) * cos(hourangle)
      cosz      = max(cosz, 0.0_wp)
   end function solar_cosz

   !----- Daylength [h], 24/pi acos(-tan(lat) tan(decl)): 0 in the polar night, 24 in the  !
   !      polar day. ----------------------------------------------------------------------!
   elemental real(wp) function daylength(lat_deg, doy) result(dl)
      real(wp),    intent(in) :: lat_deg
      integer(ik), intent(in) :: doy
      real(wp) :: latr, decl, arg
      latr = lat_deg * pi / 180.0_wp
      decl = solar_declination(doy)
      arg  = -tan(latr) * tan(decl)
      if (arg >= 1.0_wp) then
         dl = 0.0_wp                 ! polar night
      else if (arg <= -1.0_wp) then
         dl = 24.0_wp                ! polar day
      else
         dl = 24.0_wp / pi * acos(arg)
      end if
   end function daylength

   !----- The day of year, shifted by half a year in the southern hemisphere, so a         !
   !      seasonal cue written for the north serves both hemispheres. ---------------------!
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
   ! Julian Day Number conversions (Fliegel & Van Flandern 1968). Day numbers are          !
   ! 64-bit, so a span in seconds stays exact over centuries: 32-bit seconds overflow      !
   ! after 68 years.                                                                       !
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
   ! Date stepping.                                                                        !
   !=======================================================================================!
   !----- Advance by ndays whole days; the time of day is kept. ---------------------------!
   pure function time_advance_days(t, ndays) result(t2)
      type(meds_time_t), intent(in) :: t
      integer(ik),       intent(in) :: ndays
      type(meds_time_t) :: t2
      t2 = t
      call from_jdn(to_jdn(t) + int(ndays, lk), t2%year, t2%month, t2%day)
   end function time_advance_days

   !----- Advance by whole months; a day past the end of the target month becomes its last !
   !      day (31 January + 1 month is 28 or 29 February). --------------------------------!
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

   !----- Whole days from a to b (positive if b is later). --------------------------------!
   pure integer(ik) function days_between(a, b) result(d)
      type(meds_time_t), intent(in) :: a, b
      d = int(to_jdn(b) - to_jdn(a), ik)
   end function days_between

   !----- b - a in decimal calendar years, for progress messages. -------------------------!
   pure real(wp) function years_between(a, b) result(y)
      type(meds_time_t), intent(in) :: a, b
      y = time_to_decimal_year(b) - time_to_decimal_year(a)
   end function years_between

   !=======================================================================================!
   ! Arithmetic to the second, across day boundaries (the forcing's interpolation window). !
   !=======================================================================================!
   pure integer(lk) function sec_of_day(t) result(s)
      type(meds_time_t), intent(in) :: t
      s = int(t%hour, lk) * 3600_lk + int(t%minute, lk) * 60_lk + int(t%second, lk)
   end function sec_of_day

   !----- Seconds since midnight of the instant's own day (0 .. 86399). -------------------!
   pure real(wp) function seconds_into_day(t) result(s)
      type(meds_time_t), intent(in) :: t
      s = real(sec_of_day(t), wp)
   end function seconds_into_day

   !----- Seconds from a to b (positive if b is later). -----------------------------------!
   pure real(wp) function seconds_between(a, b) result(s)
      type(meds_time_t), intent(in) :: a, b
      s = real((to_jdn(b) - to_jdn(a)) * 86400_lk + sec_of_day(b) - sec_of_day(a), wp)
   end function seconds_between

   !----- Advance by dsec seconds, rounded to the nearest second; dsec may be negative or  !
   !      span many days. -----------------------------------------------------------------!
   pure function time_advance_seconds(t, dsec) result(t2)
      type(meds_time_t), intent(in) :: t
      real(wp),          intent(in) :: dsec
      type(meds_time_t) :: t2
      integer(lk) :: total, ndays, rem
      total = sec_of_day(t) + nint(dsec, lk)             ! seconds from midnight of t's day
      ndays = int(floor(real(total, wp) / 86400.0_wp), lk)
      rem   = total - ndays * 86400_lk                   ! 0 .. 86399, also when total < 0
      t2 = time_advance_days(t, int(ndays, ik))
      t2%hour   = int(rem / 3600_lk, ik)
      t2%minute = int(mod(rem, 3600_lk) / 60_lk, ik)
      t2%second = int(mod(rem, 60_lk), ik)
   end function time_advance_seconds

   !=======================================================================================!
   ! Comparisons, on the full instant.                                                     !
   !=======================================================================================!
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

   !----- Advance by n calendar years, keeping the month, day and time of day; 29 February !
   !      in a non-leap target year becomes 28 February, as in the forcing's recycling     !
   !      (docs/science/forcing.md §9). ---------------------------------------------------!
   pure function time_advance_years(t, n) result(out)
      type(meds_time_t), intent(in) :: t
      integer(ik),       intent(in) :: n
      type(meds_time_t) :: out
      out = t
      out%year = t%year + n
      if (t%month == 2_ik .and. t%day == 29_ik .and. .not. is_leap_year(out%year)) out%day = 28_ik
   end function time_advance_years

   !----- n when b is exactly n calendar years after a (1 <= n <= nmax), else 0. The       !
   !      forcing's recycle window must span whole years, or every wrap shifts the hour of !
   !      day and the season (docs/science/forcing.md §9). --------------------------------!
   pure function whole_years_between(a, b, nmax) result(n)
      type(meds_time_t), intent(in) :: a, b
      integer(ik),       intent(in) :: nmax
      integer(ik) :: n, k
      n = 0_ik
      do k = 1_ik, nmax
         if (time_eq(time_advance_years(a, k), b)) then ; n = k ; return ; end if
         if (time_lt(b, time_advance_years(a, k))) return       ! past b: not a whole number of years
      end do
   end function whole_years_between

   !----- True when the month, day, hour, minute and second are in range (a parsed config  !
   !      date); any year passes. ---------------------------------------------------------!
   pure logical function time_valid(t) result(ok)
      type(meds_time_t), intent(in) :: t
      ok = t%month >= 1_ik  .and. t%month <= 12_ik                                           &
           .and. t%day  >= 1_ik  .and. t%day   <= days_in_month(t%year, max(min(t%month,12_ik),1_ik)) &
           .and. t%hour   >= 0_ik .and. t%hour   <= 23_ik                                    &
           .and. t%minute >= 0_ik .and. t%minute <= 59_ik                                    &
           .and. t%second >= 0_ik .and. t%second <= 59_ik
   end function time_valid

   !=======================================================================================!
   ! String <-> time.                                                                      !
   !=======================================================================================!
   !----- Parse "YYYY-MM-DD" or "YYYY-MM-DD HH:MM:SS" ('/', 'T' and ',' also separate). A  !
   !      missing trailing field defaults (month and day to 1, time to 0). ok is false on  !
   !      an unreadable or out-of-range date, so the caller can use its default. ----------!
   subroutine time_from_string(s, t, ok)
      character(len=*),  intent(in)  :: s
      type(meds_time_t), intent(out) :: t
      logical,           intent(out) :: ok
      character(len=len(s)) :: buf
      integer(ik) :: i, ntok, ios
      integer(ik) :: y, mo, d, h, mi, se

      t = meds_time_t()
      y = 2000_ik ; mo = 1_ik ; d = 1_ik ; h = 0_ik ; mi = 0_ik ; se = 0_ik

      !----- Separators become blanks; a list-directed read then takes the numbers. -------!
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

   !----- The base date of a CF time-units string, "<unit> since YYYY-MM-DD[ HH:MM:SS]";   !
   !      the caller checks the unit. ok is false without "since" or with an unreadable    !
   !      date. ---------------------------------------------------------------------------!
   subroutine time_units_base(units, base, ok)
      character(len=*),  intent(in)  :: units
      type(meds_time_t), intent(out) :: base
      logical,           intent(out) :: ok
      integer :: k
      k = index(units, 'since')
      if (k == 0) then ; ok = .false. ; return ; end if
      call time_from_string(adjustl(units(k + 5:)), base, ok)
   end subroutine time_units_base

   !----- "YYYY-MM-DD HH:MM:SS", for attributes and log lines. ----------------------------!
   pure function time_to_string(t) result(s)
      type(meds_time_t), intent(in) :: t
      character(len=19) :: s
      write(s, '(i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2)')                         &
            t%year, t%month, t%day, t%hour, t%minute, t%second
   end function time_to_string

   !----- "YYYYMMDDHHMMSS", the stamp in a state file's name (<prefix>-S-<stamp>.nc). -----!
   pure function time_to_stamp(t) result(s)
      type(meds_time_t), intent(in) :: t
      character(len=14) :: s
      write(s, '(i4.4,5i2.2)') t%year, t%month, t%day, t%hour, t%minute, t%second
   end function time_to_stamp

   !----- The year plus the elapsed fraction of it: a continuous time axis for plots. -----!
   pure real(wp) function time_to_decimal_year(t) result(yr)
      type(meds_time_t), intent(in) :: t
      real(wp) :: frac
      frac = (real(day_of_year(t) - 1_ik, wp)                                                &
              + real(sec_of_day(t), wp) / 86400.0_wp) / real(days_in_year(t%year), wp)
      yr = real(t%year, wp) + frac
   end function time_to_decimal_year

   !----- The number of blank-separated tokens in s. --------------------------------------!
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
