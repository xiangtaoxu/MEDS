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
   use meds_kinds, only : wp, ik, lk
   implicit none
   private

   public :: meds_time_t
   public :: is_leap_year, days_in_month, days_in_year, day_of_year
   public :: time_advance_days, time_advance_months, days_between, years_between
   public :: time_lt, time_le, time_eq, time_valid
   public :: time_from_string, time_to_string, time_to_stamp, time_to_decimal_year

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
