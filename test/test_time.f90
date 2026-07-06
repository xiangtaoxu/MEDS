!----- Calendar arithmetic: leap years, month/day roll-over, day-of-year, parse/format. ------!
program test_time
   use meds_kinds,        only : wp, ik
   use meds_time,         only : meds_time_t, is_leap_year, days_in_month, days_in_year,       &
                                 day_of_year, time_advance_days, time_advance_months,          &
                                 days_between, time_lt, time_le, time_eq,                       &
                                 time_from_string, time_to_string, time_to_stamp,              &
                                 time_to_decimal_year, solar_cosz
   use meds_test_support, only : check, check_close, banner
   implicit none

   type(meds_time_t) :: a, b
   logical           :: ok
   real(wp)          :: cz_noon, cz_mid, cz_summer, cz_winter

   call banner('calendar / time tracking')

   !=== Leap-year rule (incl. the century exceptions). ====================================!
   call check(      is_leap_year(2000_ik), '2000 is a leap year (div by 400)')
   call check(.not. is_leap_year(1900_ik), '1900 is NOT a leap year (div by 100, not 400)')
   call check(      is_leap_year(2004_ik), '2004 is a leap year')
   call check(.not. is_leap_year(2001_ik), '2001 is not a leap year')

   !=== Days in month / year. ============================================================!
   call check(days_in_month(2000_ik, 2_ik) == 29_ik, 'Feb 2000 has 29 days')
   call check(days_in_month(1900_ik, 2_ik) == 28_ik, 'Feb 1900 has 28 days')
   call check(days_in_month(2001_ik, 4_ik) == 30_ik, 'April has 30 days')
   call check(days_in_month(2001_ik, 1_ik) == 31_ik, 'January has 31 days')
   call check(days_in_year(2000_ik) == 366_ik, '2000 has 366 days')
   call check(days_in_year(2001_ik) == 365_ik, '2001 has 365 days')

   !=== Day stepping across leap / non-leap February and the year boundary. ==============!
   a = time_advance_days(meds_time_t(2000_ik, 2_ik, 28_ik), 1_ik)
   call check(a%month == 2_ik .and. a%day == 29_ik, '2000-02-28 + 1 day -> 02-29 (leap)')
   a = time_advance_days(meds_time_t(1900_ik, 2_ik, 28_ik), 1_ik)
   call check(a%month == 3_ik .and. a%day == 1_ik,  '1900-02-28 + 1 day -> 03-01 (non-leap)')
   a = time_advance_days(meds_time_t(2000_ik, 12_ik, 31_ik), 1_ik)
   call check(a%year == 2001_ik .and. a%month == 1_ik .and. a%day == 1_ik, 'year roll-over on +1 day')
   a = time_advance_days(meds_time_t(2000_ik, 1_ik, 1_ik), days_in_year(2000_ik))
   call check(a%year == 2001_ik .and. a%month == 1_ik .and. a%day == 1_ik, '+366 days from a leap Jan 1')

   !=== Month stepping clamps the day into the shorter target month. =====================!
   a = time_advance_months(meds_time_t(2000_ik, 1_ik, 31_ik), 1_ik)
   call check(a%month == 2_ik .and. a%day == 29_ik, 'Jan 31 + 1 month -> Feb 29 (leap clamp)')
   a = time_advance_months(meds_time_t(2001_ik, 1_ik, 31_ik), 1_ik)
   call check(a%month == 2_ik .and. a%day == 28_ik, 'Jan 31 + 1 month -> Feb 28 (non-leap clamp)')
   a = time_advance_months(meds_time_t(2000_ik, 12_ik, 15_ik), 1_ik)
   call check(a%year == 2001_ik .and. a%month == 1_ik .and. a%day == 15_ik, 'Dec + 1 month rolls the year')
   a = time_advance_months(meds_time_t(2000_ik, 6_ik, 10_ik), 12_ik)
   call check(a%year == 2001_ik .and. a%month == 6_ik .and. a%day == 10_ik, '+12 months -> same day next year')

   !=== Day-of-year. =====================================================================!
   call check(day_of_year(meds_time_t(2000_ik, 1_ik,  1_ik))  == 1_ik,   'Jan 1 is DOY 1')
   call check(day_of_year(meds_time_t(2000_ik, 3_ik,  1_ik))  == 61_ik,  'Mar 1 (leap) is DOY 61')
   call check(day_of_year(meds_time_t(2000_ik, 12_ik, 31_ik)) == 366_ik, 'Dec 31 (leap) is DOY 366')
   call check(day_of_year(meds_time_t(2001_ik, 12_ik, 31_ik)) == 365_ik, 'Dec 31 (non-leap) is DOY 365')

   !=== Solar zenith cosine (derived from the time dimension). ============================!
   cz_noon   = solar_cosz(meds_time_t(2001_ik, 3_ik, 21_ik),  43200.0_wp,  0.0_wp)   ! equinox, equator, noon
   cz_mid    = solar_cosz(meds_time_t(2001_ik, 3_ik, 21_ik),      0.0_wp,  0.0_wp)   ! equinox, equator, midnight
   cz_summer = solar_cosz(meds_time_t(2001_ik, 6_ik, 21_ik),  43200.0_wp, 45.0_wp)   ! N-summer solstice noon, 45N
   cz_winter = solar_cosz(meds_time_t(2001_ik, 12_ik, 21_ik), 43200.0_wp, 45.0_wp)   ! N-winter solstice noon, 45N
   call check(abs(cz_noon - 1.0_wp) < 0.02_wp, 'equinox equator noon: cosz ~ 1')
   call check(cz_mid == 0.0_wp,                'midnight: cosz floored to 0')
   call check(cz_summer > cz_winter .and. cz_winter > 0.0_wp, '45N summer noon sun higher than winter')

   !=== Days between counts leap days. ===================================================!
   call check(days_between(meds_time_t(2000_ik,1_ik,1_ik), meds_time_t(2001_ik,1_ik,1_ik)) == 366_ik, &
              '2000 spans 366 days')
   call check(days_between(meds_time_t(2001_ik,1_ik,1_ik), meds_time_t(2002_ik,1_ik,1_ik)) == 365_ik, &
              '2001 spans 365 days')

   !=== Comparisons. =====================================================================!
   a = meds_time_t(2000_ik, 1_ik, 1_ik) ; b = meds_time_t(2000_ik, 1_ik, 2_ik)
   call check(time_lt(a, b),       'Jan 1 < Jan 2')
   call check(.not. time_lt(b, a), 'Jan 2 not < Jan 1')
   call check(time_le(a, a),       'a <= a')
   call check(time_eq(a, meds_time_t(2000_ik,1_ik,1_ik)), 'equal dates compare equal')
   a = meds_time_t(2000_ik,1_ik,1_ik,10_ik,0_ik,0_ik) ; b = meds_time_t(2000_ik,1_ik,1_ik,11_ik,0_ik,0_ik)
   call check(time_lt(a, b), 'same day, earlier hour is less')

   !=== Decimal calendar year (mid-leap-year is exactly x.5). ============================!
   call check_close(time_to_decimal_year(meds_time_t(2000_ik,1_ik,1_ik)), 2000.0_wp, 1.0e-12_wp, &
                    'Jan 1 -> integer year')
   call check_close(time_to_decimal_year(meds_time_t(2000_ik,7_ik,2_ik)), 2000.5_wp, 1.0e-12_wp, &
                    'DOY 184 of 366 -> year + 0.5')

   !=== String round-trips and validation. ===============================================!
   call time_from_string('2000-01-01', a, ok)
   call check(ok .and. a%year==2000_ik .and. a%month==1_ik .and. a%day==1_ik, 'parse YYYY-MM-DD')
   call time_from_string('1999-12-31 23:59:59', a, ok)
   call check(ok .and. a%hour==23_ik .and. a%minute==59_ik .and. a%second==59_ik, 'parse date + time')
   call time_from_string('2000/06/15', a, ok)
   call check(ok .and. a%month==6_ik .and. a%day==15_ik, 'parse with / separators')
   call time_from_string('2000-02-29', a, ok)
   call check(ok, '2000-02-29 is a valid (leap) date')
   call time_from_string('2001-02-29', a, ok)
   call check(.not. ok, '2001-02-29 is rejected (non-leap)')
   call time_from_string('2000-13-01', a, ok)
   call check(.not. ok, 'month 13 is rejected')
   call time_from_string('not-a-date', a, ok)
   call check(.not. ok, 'garbage is rejected')

   call check(time_to_stamp(meds_time_t(2050_ik,1_ik,1_ik)) == '20500101000000', '14-char stamp')
   call check(time_to_stamp(meds_time_t(1999_ik,12_ik,31_ik,23_ik,59_ik,59_ik)) == '19991231235959', &
              'stamp carries the time of day')
   call check(time_to_string(meds_time_t(2000_ik,1_ik,1_ik)) == '2000-01-01 00:00:00', 'human string')

   write(*,'(a)') '[test] calendar / time tracking: PASS'
end program test_time
