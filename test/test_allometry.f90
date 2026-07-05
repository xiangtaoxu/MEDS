!----- Pan-tropical allometry: dbh<->height round-trip, AGB monotonicity, agb_to_dbh inverse. !
program test_allometry
   use meds_kinds,       only : wp, ik
   use meds_allometry,   only : dbh_to_height, height_to_dbh, dbh_to_agb, agb_to_dbh,         &
                                dbh_to_leaf_area, set_allometry,                              &
                                size2leaf_carbon, size2wood_carbon, wood_to_dbh
   use meds_test_support, only : check, check_close, banner
   implicit none

   real(wp), parameter :: rho  = 0.60_wp
   real(wp), parameter :: hcap = 46.0_wp        ! per-PFT asymptotic height cap (an argument now)
   real(wp), parameter :: sla  = 13.0_wp        ! [m2/kgC] specific leaf area
   real(wp), parameter :: agf  = 0.70_wp        ! aboveground fraction of woody carbon
   real(wp) :: dbh, h, agb, dbh_back, wood_c
   integer(ik) :: i
   real(wp), dimension(3) :: dvec = [ 5.0_wp, 30.0_wp, 80.0_wp ]

   call banner('pan-tropical allometry')
   !----- Install the global allometry coefficients (the model no longer hard-codes them). --!
   call set_allometry(1.139963_wp, 0.564899_wp, 0.06080334_wp, 1.0044785_wp,         &
                      0.370_wp, 0.464_wp, 0.46769540_wp, 0.6410495_wp, 0.5_wp)

   !=== Height <-> diameter round-trip below the height cap. ==============================!
   do i = 1_ik, 3_ik
      dbh = dvec(i)
      h   = dbh_to_height(dbh, hcap)
      call check(h < hcap, 'test sizes must stay below the height cap')
      call check_close(height_to_dbh(h), dbh, 1.0e-9_wp, 'dbh<->height round-trip failed')
   end do

   !=== Height saturates at the cap for very large diameters. ============================!
   call check_close(dbh_to_height(500.0_wp, hcap), hcap, 1.0e-12_wp, 'height did not cap at hgt_max')

   !=== AGB and leaf area increase monotonically with diameter. ==========================!
   call check(dbh_to_agb(20.0_wp, dbh_to_height(20.0_wp, hcap), rho) >                        &
              dbh_to_agb(10.0_wp, dbh_to_height(10.0_wp, hcap), rho), 'AGB not increasing with dbh')
   call check(dbh_to_leaf_area(20.0_wp, dbh_to_height(20.0_wp, hcap)) >                       &
              dbh_to_leaf_area(10.0_wp, dbh_to_height(10.0_wp, hcap)), 'leaf area not increasing with dbh')
   call check(dbh_to_agb(30.0_wp, dbh_to_height(30.0_wp, hcap), rho) > 0.0_wp, 'AGB must be positive')

   !=== AGB heavier for denser wood at the same size. ====================================!
   call check(dbh_to_agb(30.0_wp, dbh_to_height(30.0_wp, hcap), 0.90_wp) >                    &
              dbh_to_agb(30.0_wp, dbh_to_height(30.0_wp, hcap), 0.40_wp), 'denser wood should weigh more')

   !=== agb_to_dbh inverts dbh_to_agb in BOTH regimes. ===================================!
   !----- uncapped (small): height(30) < hcap. ------------------------------------------------!
   dbh      = 30.0_wp
   agb      = dbh_to_agb(dbh, dbh_to_height(dbh, hcap), rho)
   dbh_back = agb_to_dbh(agb, rho, hcap)
   call check_close(dbh_back, dbh, 1.0e-8_wp, 'agb_to_dbh inverse failed (uncapped)')
   !----- capped (large): height(130) saturates at hcap. --------------------------------------!
   dbh      = 130.0_wp
   call check_close(dbh_to_height(dbh, hcap), hcap, 1.0e-12_wp, 'large stem should be at the cap')
   agb      = dbh_to_agb(dbh, dbh_to_height(dbh, hcap), rho)
   dbh_back = agb_to_dbh(agb, rho, hcap)
   call check_close(dbh_back, dbh, 1.0e-8_wp, 'agb_to_dbh inverse failed (capped)')

   !=== Carbon-pool targets: EXACT re-expressions of the size allometry (behavior-preserving). =!
   do i = 1_ik, 3_ik
      dbh = dvec(i)
      h   = dbh_to_height(dbh, hcap)
      !----- leaf_carbon * SLA reproduces the leaf area exactly (LAI unchanged). ---------------!
      call check_close(size2leaf_carbon(dbh, h, sla) * sla, dbh_to_leaf_area(dbh, h), 1.0e-12_wp, &
                       'size2leaf_carbon * sla /= leaf_area')
      !----- aboveground_frac * wood_carbon reproduces the Chave AGB exactly. ------------------!
      call check_close(size2wood_carbon(dbh, h, rho, agf) * agf, dbh_to_agb(dbh, h, rho), 1.0e-12_wp, &
                       'aboveground_frac * wood_carbon /= agb')
   end do

   !=== Carbon targets increase with diameter. ==========================================!
   call check(size2leaf_carbon(20.0_wp, dbh_to_height(20.0_wp, hcap), sla) >                   &
              size2leaf_carbon(10.0_wp, dbh_to_height(10.0_wp, hcap), sla), 'leaf carbon not increasing')
   call check(size2wood_carbon(20.0_wp, dbh_to_height(20.0_wp, hcap), rho, agf) >              &
              size2wood_carbon(10.0_wp, dbh_to_height(10.0_wp, hcap), rho, agf), 'wood carbon not increasing')

   !=== wood_to_dbh inverts size2wood_carbon in BOTH regimes (the carbon-prognostic anchor). ==!
   dbh    = 30.0_wp                                        ! uncapped
   wood_c = size2wood_carbon(dbh, dbh_to_height(dbh, hcap), rho, agf)
   call check_close(wood_to_dbh(wood_c, rho, hcap, agf), dbh, 1.0e-8_wp, 'wood_to_dbh inverse failed (uncapped)')
   dbh    = 130.0_wp                                       ! capped
   wood_c = size2wood_carbon(dbh, dbh_to_height(dbh, hcap), rho, agf)
   call check_close(wood_to_dbh(wood_c, rho, hcap, agf), dbh, 1.0e-8_wp, 'wood_to_dbh inverse failed (capped)')

   write(*,'(a)') '   PASS'
end program test_allometry
