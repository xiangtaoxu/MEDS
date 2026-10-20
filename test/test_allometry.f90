! SPDX-License-Identifier: Apache-2.0
!----- Pan-tropical allometry: dbh<->height round-trip, AGB monotonicity, agb_to_dbh inverse. !
program test_allometry
   use meds_kinds,       only : wp, ik
   use meds_allometry_lib, only : dbh_to_height, height_to_dbh, dbh_to_agb, agb_to_dbh, dbh_to_leaf_area,  &
                                  allometry_t, size_to_leaf_carbon, size_to_wood_carbon, wood_to_dbh,      &
                                  HEIGHT_CHAVE
   use meds_test_support, only : banner, check, check_close
   implicit none

   real(wp), parameter :: rho  = 0.60_wp
   real(wp), parameter :: hcap = 46.0_wp        ! the record's height cap (its initializer)
   real(wp), parameter :: sla  = 13.0_wp        ! [m2/kgC] specific leaf area
   real(wp), parameter :: agf  = 0.70_wp        ! aboveground fraction of woody carbon
   real(wp) :: dbh, h, agb, dbh_back, wood_c
   integer(ik) :: i
   type(allometry_t) :: allom, two(2)
   real(wp)          :: hh(2)
   real(wp), dimension(3) :: dvec = [ 5.0_wp, 30.0_wp, 80.0_wp ]

   call banner('pan-tropical allometry')
   !----- The shipped pan-tropical coefficients: the record's initializers. ----------------!
   allom = allometry_t()
   call check(allom%hgt_max == hcap .and. allom%aboveground_frac == agf, 'the record''s cap and aboveground share')

   !=== Height <-> diameter round-trip below the height cap. ==============================!
   do i = 1_ik, 3_ik
      dbh = dvec(i)
      h   = dbh_to_height(allom, dbh)
      call check(h < hcap, 'test sizes must stay below the height cap')
      call check_close(height_to_dbh(allom, h), dbh, 1.0e-9_wp, 'dbh<->height round-trip failed')
   end do

   !=== Height saturates at the cap for very large diameters. ============================!
   call check_close(dbh_to_height(allom, 500.0_wp), hcap, 1.0e-12_wp, 'height did not cap at hgt_max')

   !=== AGB and leaf area increase monotonically with diameter. ==========================!
   call check(dbh_to_agb(allom, 20.0_wp, dbh_to_height(allom, 20.0_wp), rho) >                              &
              dbh_to_agb(allom, 10.0_wp, dbh_to_height(allom, 10.0_wp), rho), 'AGB not increasing with dbh')
   call check(dbh_to_leaf_area(allom, 20.0_wp, dbh_to_height(allom, 20.0_wp)) >                             &
              dbh_to_leaf_area(allom, 10.0_wp, dbh_to_height(allom, 10.0_wp)), 'leaf area not increasing with dbh')
   call check(dbh_to_agb(allom, 30.0_wp, dbh_to_height(allom, 30.0_wp), rho) > 0.0_wp, 'AGB must be positive')

   !=== AGB heavier for denser wood at the same size. ====================================!
   call check(dbh_to_agb(allom, 30.0_wp, dbh_to_height(allom, 30.0_wp), 0.90_wp) >                          &
              dbh_to_agb(allom, 30.0_wp, dbh_to_height(allom, 30.0_wp), 0.40_wp), 'denser wood should weigh more')

   !=== agb_to_dbh inverts dbh_to_agb in BOTH regimes. ===================================!
   !----- uncapped (small): height(30) < hcap. ------------------------------------------------!
   dbh      = 30.0_wp
   agb      = dbh_to_agb(allom, dbh, dbh_to_height(allom, dbh), rho)
   dbh_back = agb_to_dbh(allom, agb, rho)
   call check_close(dbh_back, dbh, 1.0e-8_wp, 'agb_to_dbh inverse failed (uncapped)')
   !----- capped (large): height(130) saturates at hcap. --------------------------------------!
   dbh      = 130.0_wp
   call check_close(dbh_to_height(allom, dbh), hcap, 1.0e-12_wp, 'large stem should be at the cap')
   agb      = dbh_to_agb(allom, dbh, dbh_to_height(allom, dbh), rho)
   dbh_back = agb_to_dbh(allom, agb, rho)
   call check_close(dbh_back, dbh, 1.0e-8_wp, 'agb_to_dbh inverse failed (capped)')

   !=== Carbon-pool targets: EXACT re-expressions of the size allometry. =====================!
   do i = 1_ik, 3_ik
      dbh = dvec(i)
      h   = dbh_to_height(allom, dbh)
      !----- leaf_carbon * SLA reproduces the leaf area exactly. -------------------------------!
      call check_close(size_to_leaf_carbon(allom, dbh, h, sla) * sla, dbh_to_leaf_area(allom, dbh, h), 1.0e-12_wp, &
                       'size_to_leaf_carbon * sla /= leaf_area')
      !----- aboveground_frac * wood_carbon reproduces the Chave AGB exactly. ------------------!
      call check_close(size_to_wood_carbon(allom, dbh, h, rho) * agf, dbh_to_agb(allom, dbh, h, rho), 1.0e-12_wp, &
                       'aboveground_frac * wood_carbon /= agb')
   end do

   !=== Carbon targets increase with diameter. ==========================================!
   call check(size_to_leaf_carbon(allom, 20.0_wp, dbh_to_height(allom, 20.0_wp), sla) >                         &
              size_to_leaf_carbon(allom, 10.0_wp, dbh_to_height(allom, 10.0_wp), sla), 'leaf carbon not increasing')
   call check(size_to_wood_carbon(allom, 20.0_wp, dbh_to_height(allom, 20.0_wp), rho) >                         &
              size_to_wood_carbon(allom, 10.0_wp, dbh_to_height(allom, 10.0_wp), rho), 'wood carbon not increasing')

   !=== wood_to_dbh inverts size_to_wood_carbon in BOTH regimes (the carbon anchor). =============!
   dbh    = 30.0_wp                                        ! uncapped
   wood_c = size_to_wood_carbon(allom, dbh, dbh_to_height(allom, dbh), rho)
   call check_close(wood_to_dbh(allom, wood_c, rho), dbh, 1.0e-8_wp, 'wood_to_dbh inverse failed (uncapped)')
   dbh    = 130.0_wp                                       ! capped
   wood_c = size_to_wood_carbon(allom, dbh, dbh_to_height(allom, dbh), rho)
   call check_close(wood_to_dbh(allom, wood_c, rho), dbh, 1.0e-8_wp, 'wood_to_dbh inverse failed (capped)')

   !=== Chave et al. (2014) eq. 6a at E = 0.2 and its two inverses. The expected heights are     !
   !    BIOMASS's retrieveH: exp(0.893 - E + 0.760 ln D - 0.0340 (ln D)^2 + 0.243^2/2). =========!
   allom = allometry_t(height_form = HEIGHT_CHAVE, chave_e = 0.2_wp, hgt_max = 100.0_wp)
   call check_close(dbh_to_height(allom, 1.0_wp), 2.0596261794_wp, 1.0e-9_wp, 'Chave height at 1 cm')
   call check_close(dbh_to_height(allom, 10.0_wp), 9.8969302660_wp, 1.0e-9_wp, 'Chave height at 10 cm')
   call check_close(dbh_to_height(allom, 30.0_wp), 18.4325408427_wp, 1.0e-9_wp, 'Chave height at 30 cm')
   call check_close(dbh_to_height(allom, 60.0_wp), 26.1607976122_wp, 1.0e-9_wp, 'Chave height at 60 cm')
   call check_close(dbh_to_height(capped(allom, 30.0_wp), 100.0_wp), 30.0_wp, 1.0e-12_wp, 'Chave height caps at hgt_max')
   do i = 1_ik, 3_ik
      dbh = dvec(i)
      call check_close(height_to_dbh(allom, dbh_to_height(allom, dbh)), dbh, 1.0e-9_wp, 'Chave dbh<->height round trip')
      agb = dbh_to_agb(allom, dbh, dbh_to_height(allom, dbh), rho)
      call check_close(agb_to_dbh(allom, agb, rho), dbh, 1.0e-9_wp, 'Chave agb_to_dbh inverse (uncapped)')
   end do
   dbh = 120.0_wp
   call check_close(dbh_to_height(capped(allom, 30.0_wp), dbh), 30.0_wp, 1.0e-12_wp, 'a 120 cm stem is at a 30 m cap')
   agb = dbh_to_agb(allom, dbh, dbh_to_height(capped(allom, 30.0_wp), dbh), rho)
   call check_close(agb_to_dbh(capped(allom, 30.0_wp), agb, rho), dbh, 1.0e-9_wp, 'Chave agb_to_dbh inverse (capped)')
   wood_c = size_to_wood_carbon(allom, 80.0_wp, dbh_to_height(allom, 80.0_wp), rho)
   call check_close(wood_to_dbh(allom, wood_c, rho), 80.0_wp, 1.0e-9_wp, 'Chave wood_to_dbh inverse')
   !----- A PFT's own log-quadratic coefficients replace Chave's: a flatter curve is shorter. --!
   call check(dbh_to_height(allometry_t(height_form = HEIGHT_CHAVE, chave_b = 0.6_wp), 30.0_wp) <              &
              dbh_to_height(allometry_t(height_form = HEIGHT_CHAVE), 30.0_wp), 'a smaller chave_b gives shorter trees')

   !=== A local log-quadratic (Barro Colorado's: chave_a, chave_b, chave_c, no E, no RSE term)   !
   !    follows its own coefficients, and both inverses hold on it. =======================!
   allom = allometry_t(height_form = HEIGHT_CHAVE, chave_a = 0.8898_wp, chave_b = 0.7908_wp,          &
                       chave_c = -0.04923_wp, chave_rse = 0.0_wp, chave_e = 0.0_wp, hgt_max = 100.0_wp)
   call check_close(dbh_to_height(allom, 1.0_wp), exp(0.8898_wp), 1.0e-12_wp, 'local fit: height at 1 cm is exp(chave_a)')
   call check_close(dbh_to_height(allom, 100.0_wp),                                                          &
                    exp(0.8898_wp + 0.7908_wp * log(100.0_wp) - 0.04923_wp * log(100.0_wp)**2), 1.0e-12_wp, &
                    'local fit: height at 100 cm follows its coefficients')
   do i = 1_ik, 3_ik
      dbh = dvec(i)
      call check_close(height_to_dbh(allom, dbh_to_height(allom, dbh)), dbh, 1.0e-9_wp, 'local fit: dbh<->height')
      agb = dbh_to_agb(allom, dbh, dbh_to_height(allom, dbh), rho)
      call check_close(agb_to_dbh(allom, agb, rho), dbh, 1.0e-9_wp, 'local fit: agb_to_dbh inverse')
   end do

   !=== Each PFT has its own record: a short log-quadratic understorey PFT and a power-law PFT     !
   !    keep their own heights and caps side by side (elemental over the records). ============!
   two(1) = allometry_t(height_form = HEIGHT_CHAVE, chave_a = 0.6_wp, chave_b = 0.6_wp, chave_c = -0.03_wp,  &
                        chave_rse = 0.0_wp, hgt_max = 12.0_wp)
   two(2) = allometry_t()
   hh = dbh_to_height(two, 40.0_wp)
   call check_close(hh(1), min(exp(0.6_wp + 0.6_wp * log(40.0_wp) - 0.03_wp * log(40.0_wp)**2), 12.0_wp), 1.0e-12_wp, &
                    'PFT 1 follows its own log-quadratic and cap')
   call check_close(hh(2), dbh_to_height(allometry_t(), 40.0_wp), 0.0_wp, 'PFT 2 follows the power law')
   call check(hh(1) < hh(2), 'the understorey PFT is shorter at the same diameter')

   write(*,'(a)') '   PASS'

contains

   !----- A copy of a record with another height cap. ------------------------------------!
   pure function capped(a, hgt_max) result(b)
      type(allometry_t), intent(in) :: a
      real(wp),          intent(in) :: hgt_max
      type(allometry_t)             :: b
      b = a
      b%hgt_max = hgt_max
   end function capped

end program test_allometry
