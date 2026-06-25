!----- Pan-tropical allometry: dbh<->height round-trip, AGB monotonicity, agb_to_dbh inverse. !
program test_allometry
   use meds_kinds,       only : wp, ik
   use meds_allometry,   only : dbh_to_height, height_to_dbh, dbh_to_agb, agb_to_dbh,         &
                                dbh_to_leaf_area, hgt_max
   use meds_test_support, only : check, check_close, banner
   implicit none

   real(wp), parameter :: rho = 0.60_wp
   real(wp) :: dbh, h, agb, dbh_back
   integer(ik) :: i
   real(wp), dimension(3) :: dvec = [ 5.0_wp, 30.0_wp, 80.0_wp ]

   call banner('pan-tropical allometry')

   !=== Height <-> diameter round-trip below the height cap. ==============================!
   do i = 1_ik, 3_ik
      dbh = dvec(i)
      h   = dbh_to_height(dbh)
      call check(h < hgt_max, 'test sizes must stay below the height cap')
      call check_close(height_to_dbh(h), dbh, 1.0e-9_wp, 'dbh<->height round-trip failed')
   end do

   !=== Height saturates at the cap for very large diameters. ============================!
   call check_close(dbh_to_height(500.0_wp), hgt_max, 1.0e-12_wp, 'height did not cap at hgt_max')

   !=== AGB and leaf area increase monotonically with diameter. ==========================!
   call check(dbh_to_agb(20.0_wp, dbh_to_height(20.0_wp), rho) >                             &
              dbh_to_agb(10.0_wp, dbh_to_height(10.0_wp), rho), 'AGB not increasing with dbh')
   call check(dbh_to_leaf_area(20.0_wp, dbh_to_height(20.0_wp)) >                            &
              dbh_to_leaf_area(10.0_wp, dbh_to_height(10.0_wp)), 'leaf area not increasing with dbh')
   call check(dbh_to_agb(30.0_wp, dbh_to_height(30.0_wp), rho) > 0.0_wp, 'AGB must be positive')

   !=== AGB heavier for denser wood at the same size. ====================================!
   call check(dbh_to_agb(30.0_wp, dbh_to_height(30.0_wp), 0.90_wp) >                         &
              dbh_to_agb(30.0_wp, dbh_to_height(30.0_wp), 0.40_wp), 'denser wood should weigh more')

   !=== agb_to_dbh inverts dbh_to_agb in BOTH regimes. ===================================!
   !----- uncapped (small): height(30) < hgt_max. ------------------------------------------!
   dbh      = 30.0_wp
   agb      = dbh_to_agb(dbh, dbh_to_height(dbh), rho)
   dbh_back = agb_to_dbh(agb, rho)
   call check_close(dbh_back, dbh, 1.0e-8_wp, 'agb_to_dbh inverse failed (uncapped)')
   !----- capped (large): height(130) saturates at hgt_max. --------------------------------!
   dbh      = 130.0_wp
   call check_close(dbh_to_height(dbh), hgt_max, 1.0e-12_wp, 'large stem should be at the cap')
   agb      = dbh_to_agb(dbh, dbh_to_height(dbh), rho)
   dbh_back = agb_to_dbh(agb, rho)
   call check_close(dbh_back, dbh, 1.0e-8_wp, 'agb_to_dbh inverse failed (capped)')

   write(*,'(a)') '   PASS'
end program test_allometry
