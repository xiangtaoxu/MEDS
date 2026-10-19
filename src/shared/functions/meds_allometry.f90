! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_allometry -- a stem's size from its diameter: height, aboveground biomass (AGB),    !
! leaf area, wood area and sapwood; the way back from carbon to diameter; and the carbon-  !
! pool targets the allocator fills. The curves, their shipped coefficients and the         !
! inversions are docs/science/plant_carbon_allocation.md §2 (The size allometry).          !
!                                                                                          !
! The coefficients are one record, allometry_t, the first argument of every routine, and   !
! are the same for every PFT. A run's record is cfg%allom; each site carries a copy        !
! (cohort%allom), so the polygons of a region can each have their own                      !
! (docs/configuration.md §Regional runs). Biomass and leaf area are in carbon.             !
!==========================================================================================!
module meds_allometry
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num, pio4
   implicit none
   private

   public :: allometry_t
   public :: dbh_to_height, height_to_dbh, dbh_to_agb, agb_to_dbh, dbh_to_leaf_area, local_hmax
   public :: size2leaf_carbon, size2wood_carbon, wood_to_dbh, carbon_to_structure, min_cohort_carbon
   public :: storage_carbon_target
   public :: dbh_to_wai, sapwood_fraction
   public :: HEIGHT_POWER, HEIGHT_GMM, HEIGHT_CHAVE

   !----- Height curves: the pan-tropical power law, a saturating gMM, or Chave (2014). ---!
   integer(ik), parameter :: HEIGHT_POWER = 1_ik, HEIGHT_GMM = 2_ik, HEIGHT_CHAVE = 3_ik
   !----- Chave et al. (2014) eq. 6a: ln H = CHAVE_A - E + CHAVE_B x + CHAVE_C x^2, with   !
   !      x = ln(dbh), and the residual standard error that turns its median into the mean !
   !      (RSE^2 / 2). --------------------------------------------------------------------!
   real(wp),    parameter :: CHAVE_A   = 0.893_wp
   real(wp),    parameter :: CHAVE_B   = 0.760_wp
   real(wp),    parameter :: CHAVE_C   = -0.0340_wp
   real(wp),    parameter :: CHAVE_RSE = 0.243_wp
   !----- The local maximum height is the height of a tree this wide on the local curve:   !
   !      one definition for every curve, saturating or not (Chave's has no asymptote). ---!
   real(wp),    parameter :: HMAX_DBH  = 100.0_wp     !< [cm]
   !----- Newton's method for agb -> dbh under the gMM height (agb_to_dbh_gmm). -----------!
   integer,     parameter :: NEWTON_MAX_ITER = 50
   real(wp),    parameter :: NEWTON_TOL      = 1.0e-12_wp    !< [ln cm] change in ln(dbh)

   !----- The allometry coefficients, from the PFT file's [allometry] block. The defaults  !
   !      are the shipped pan-tropical values, so a record made without a config (a unit   !
   !      test, a Python call) is usable. -------------------------------------------------!
   type :: allometry_t
      real(wp)    :: b1Ht        = 1.139963_wp    !< [--]  power-law height: intercept of ln H on ln dbh
      real(wp)    :: b2Ht        = 0.564899_wp    !< [--]  power-law height: slope of ln H on ln dbh
      real(wp)    :: agb_c1      = 0.03365_wp     !< [kgC] AGB scale: Chave et al. (2014) eq. 4, halved to carbon
      real(wp)    :: agb_c2      = 0.976_wp       !< [--]  AGB exponent on rho and on dbh^2*h (Chave 2014)
      real(wp)    :: lai_b1      = 0.23384770_wp  !< [--]  per-stem leaf-area scale
      real(wp)    :: lai_b2      = 0.6410495_wp   !< [--]  per-stem leaf-area exponent
      real(wp)    :: light_ext   = 0.5_wp         !< [--]  Beer-Lambert extinction through overtopping LAI
      integer(ik) :: height_form = HEIGHT_POWER   !< the height curve: power law, gMM or Chave (2014)
      !----- The gMM height (Cano et al. 2019 eq. 7), read only under HEIGHT_GMM. ---------!
      real(wp)    :: gmm_a       = 58.0_wp        !< [m]   asymptotic height
      real(wp)    :: gmm_b       = 0.73_wp        !< [--]  exponent on dbh
      real(wp)    :: gmm_k       = 21.8_wp        !< [cm^gmm_b] half-saturation
      !----- Chave's stress index E, read only when height_form = HEIGHT_CHAVE. -----------!
      real(wp)    :: chave_e     = 0.0_wp         !< [--]  site stress index: the larger, the shorter the trees
   end type allometry_t

contains

   !----- Diameter [cm] -> height [m] on the configured curve, capped at hgt_max. ---------!
   elemental pure function dbh_to_height(allom, dbh, hgt_max) result(h)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, hgt_max
      real(wp)                      :: h, db, x
      select case (allom%height_form)
      case (HEIGHT_GMM)
         db = max(dbh, tiny_num) ** allom%gmm_b
         h  = allom%gmm_a * db / (allom%gmm_k + db)
      case (HEIGHT_CHAVE)
         x  = log(max(dbh, tiny_num))
         h  = exp(chave_intercept(allom) + CHAVE_B * x + CHAVE_C * x * x)
      case default
         h  = exp(allom%b1Ht + allom%b2Ht * log(max(dbh, tiny_num)))
      end select
      h = min(h, hgt_max)
   end function dbh_to_height

   !----- The local maximum height [m]: the height of a HMAX_DBH-wide tree on the local    !
   !      curve, capped at hgt_max (docs/science/plant_carbon_allocation.md §2). ----------!
   elemental pure function local_hmax(allom, hgt_max) result(h)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: hgt_max
      real(wp)                      :: h
      h = dbh_to_height(allom, HMAX_DBH, hgt_max)
   end function local_hmax

   !----- Height [m] -> diameter [cm], the inverse of the uncapped curve. -----------------!
   elemental pure function height_to_dbh(allom, h) result(dbh)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: h
      real(wp)                      :: dbh
      select case (allom%height_form)
      case (HEIGHT_GMM)
         dbh = (allom%gmm_k * max(h, tiny_num) / max(allom%gmm_a - h, tiny_num)) ** (1.0_wp / allom%gmm_b)
      case (HEIGHT_CHAVE)
         dbh = exp(chave_small_root(CHAVE_B, chave_intercept(allom) - log(max(h, tiny_num))))
      case default
         dbh = exp((log(max(h, tiny_num)) - allom%b1Ht) / allom%b2Ht)
      end select
   end function height_to_dbh

   !----- The intercept a of Chave's curve, ln H = a + CHAVE_B x + CHAVE_C x^2 with        !
   !      x = ln(dbh). --------------------------------------------------------------------!
   pure real(wp) function chave_intercept(allom)
      type(allometry_t), intent(in) :: allom
      chave_intercept = CHAVE_A - allom%chave_e + 0.5_wp * CHAVE_RSE * CHAVE_RSE
   end function chave_intercept

   !----- The smaller root x of CHAVE_C x^2 + b x + c = 0: the curve's rising branch,      !
   !      which holds up to dbh = exp(-CHAVE_B / (2 CHAVE_C)), about 70,000 cm. Written as !
   !      2(-c)/(b + sqrt(.)) so it keeps its digits when c is small. ---------------------!
   pure real(wp) function chave_small_root(b, c) result(x)
      real(wp), intent(in) :: b, c
      x = -2.0_wp * c / (b + sqrt(max(b * b - 4.0_wp * CHAVE_C * c, 0.0_wp)))
   end function chave_small_root

   !----- Aboveground biomass [kgC/plant], Chave et al. (2014) eq. 4 in carbon. -----------!
   elemental pure function dbh_to_agb(allom, dbh, h, rho) result(agb)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h, rho
      real(wp)                      :: agb
      agb = allom%agb_c1 * rho ** allom%agb_c2 * (dbh * dbh * h) ** allom%agb_c2
   end function dbh_to_agb

   !----- The AGB [kgC/plant] of a recruit at min_cohort_height: the carbon per recruit    !
   !      that turns reproduction carbon into recruit numbers. ----------------------------!
   elemental pure function min_cohort_carbon(allom, min_cohort_height, rho) result(carbon_min)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: min_cohort_height, rho
      real(wp)                      :: carbon_min
      carbon_min = dbh_to_agb(allom, height_to_dbh(allom, min_cohort_height), min_cohort_height, rho)
   end function min_cohort_carbon

   !----- One-sided leaf area of one stem [m2]; the cohort's LAI is nplant times it. ------!
   elemental pure function dbh_to_leaf_area(allom, dbh, h) result(leaf_area)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h
      real(wp)                      :: leaf_area
      leaf_area = allom%lai_b1 * (dbh * dbh * h) ** allom%lai_b2
   end function dbh_to_leaf_area

   !----- AGB -> diameter [cm], which recovers a diameter from the conserved carbon. Under !
   !      the power law agb = K dbh^P on both sides of the height cap: solve the uncapped  !
   !      form, and redo it in the capped form if that stem is taller than hgt_max. -------!
   elemental pure function agb_to_dbh(allom, agb, rho, hgt_max) result(dbh)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: agb, rho, hgt_max
      real(wp)                      :: dbh, k_un, p_un, k_cap, p_cap, h
      if (allom%height_form == HEIGHT_GMM) then
         dbh = agb_to_dbh_gmm(allom, agb, rho, hgt_max)
         return
      end if
      if (allom%height_form == HEIGHT_CHAVE) then
         dbh = agb_to_dbh_chave(allom, agb, rho, hgt_max)
         return
      end if
      k_un = allom%agb_c1 * rho ** allom%agb_c2 * exp(allom%agb_c2 * allom%b1Ht)
      p_un = allom%agb_c2 * (2.0_wp + allom%b2Ht)
      dbh  = (max(agb, tiny_num) / k_un) ** (1.0_wp / p_un)
      h    = exp(allom%b1Ht + allom%b2Ht * log(max(dbh, tiny_num)))
      if (h > hgt_max) then
         k_cap = allom%agb_c1 * rho ** allom%agb_c2 * hgt_max ** allom%agb_c2
         p_cap = 2.0_wp * allom%agb_c2
         dbh   = (max(agb, tiny_num) / k_cap) ** (1.0_wp / p_cap)
      end if
   end function agb_to_dbh

   !----- AGB -> diameter under the gMM height, which has no closed-form inverse: Newton's !
   !      method on x = ln(dbh). The biomass law, 2x + ln h(x) = target, is rising and     !
   !      concave, so after one step the iterates climb to the root without overshoot. ----!
   elemental pure function agb_to_dbh_gmm(allom, agb, rho, hgt_max) result(dbh)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: agb, rho, hgt_max
      real(wp)                      :: dbh, target, x, db, h, f, slope, dx
      integer                       :: it
      target = log(max(agb, tiny_num) / (allom%agb_c1 * rho ** allom%agb_c2)) / allom%agb_c2
      x = target / (2.0_wp + allom%gmm_b)
      do it = 1, NEWTON_MAX_ITER
         db = exp(allom%gmm_b * x)
         h  = allom%gmm_a * db / (allom%gmm_k + db)
         if (h < hgt_max) then
            f     = 2.0_wp * x + log(h) - target
            slope = 2.0_wp + allom%gmm_b * allom%gmm_k / (allom%gmm_k + db)
         else
            f     = 2.0_wp * x + log(hgt_max) - target
            slope = 2.0_wp
         end if
         dx = f / slope
         x  = x - dx
         if (abs(dx) < NEWTON_TOL) exit
      end do
      dbh = exp(x)
   end function agb_to_dbh_gmm

   !----- AGB -> diameter under Chave's height, in closed form: below the cap, the biomass !
   !      law 2x + ln h(x) = target is quadratic in x = ln(dbh); x is its smaller root. A  !
   !      stem that root makes taller than hgt_max has 2x + ln(hgt_max) = target. ---------!
   elemental pure function agb_to_dbh_chave(allom, agb, rho, hgt_max) result(dbh)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: agb, rho, hgt_max
      real(wp)                      :: dbh, target, x
      target = log(max(agb, tiny_num) / (allom%agb_c1 * rho ** allom%agb_c2)) / allom%agb_c2
      x = chave_small_root(2.0_wp + CHAVE_B, chave_intercept(allom) - target)
      if (chave_intercept(allom) + CHAVE_B * x + CHAVE_C * x * x > log(hgt_max))                  &
         x = 0.5_wp * (target - log(hgt_max))
      dbh = exp(x)
   end function agb_to_dbh_chave

   !----- The leaf carbon target [kgC/plant]: the allometric leaf area over SLA, so a      !
   !      cohort at its target shows exactly that leaf area                                !
   !      (docs/science/plant_carbon_allocation.md §2). -----------------------------------!
   elemental pure function size2leaf_carbon(allom, dbh, h, sla) result(leaf_carbon)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h, sla        !< [cm],[m],[m2/kgC]
      real(wp)                      :: leaf_carbon         !< [kgC/plant]
      leaf_carbon = dbh_to_leaf_area(allom, dbh, h) / max(sla, tiny_num)
   end function size2leaf_carbon

   !----- A cohort's storage target: storage_cushion times the carbon one full flush       !
   !      rebuilds, its leaves plus its fine roots. It is also the starting storage of a   !
   !      cohort set on allometry. --------------------------------------------------------!
   elemental pure function storage_carbon_target(storage_cushion, leaf_carbon, fineroot_carbon) result(storage)
      real(wp), intent(in) :: storage_cushion              !< [--]
      real(wp), intent(in) :: leaf_carbon, fineroot_carbon !< [kgC/plant]
      real(wp)             :: storage                      !< [kgC/plant]
      storage = storage_cushion * (leaf_carbon + fineroot_carbon)
   end function storage_carbon_target

   !----- A stem's total wood carbon [kgC/plant], coarse roots included: its AGB over the  !
   !      aboveground fraction, so the aboveground share is exactly the AGB. --------------!
   elemental pure function size2wood_carbon(allom, dbh, h, rho, aboveground_frac) result(wood_carbon)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h, rho, aboveground_frac
      real(wp)                      :: wood_carbon         !< [kgC/plant] total woody (above + belowground)
      wood_carbon = dbh_to_agb(allom, dbh, h, rho) / max(aboveground_frac, tiny_num)
   end function size2wood_carbon

   !----- The inverse of size2wood_carbon: wood carbon is the prognostic size, and the     !
   !      diameter follows from it. -------------------------------------------------------!
   elemental pure function wood_to_dbh(allom, wood_carbon, rho, hgt_max, aboveground_frac) result(dbh)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: wood_carbon, rho, hgt_max, aboveground_frac
      real(wp)                      :: dbh                 !< [cm]
      dbh = agb_to_dbh(allom, wood_carbon * aboveground_frac, rho, hgt_max)
   end function wood_to_dbh

   !----- A cohort's geometry from its carbon pools: the diameter from wood carbon, then   !
   !      the height, basal area and AGB, and the leaf area from leaf carbon. It is the    !
   !      body of set_cohort_size_from_carbon, the one writer of a cohort's geometry. -----!
   pure subroutine carbon_to_structure(allom, wood_carbon, leaf_carbon, rho, hgt_max, aboveground_frac, &
                                       sla, dbh, height, basal_area, agb, leaf_area)
      type(allometry_t), intent(in)  :: allom
      real(wp),          intent(in)  :: wood_carbon, leaf_carbon    !< [kgC/plant] prognostic pools
      real(wp),          intent(in)  :: rho, hgt_max, aboveground_frac, sla    !< PFT traits
      real(wp),          intent(out) :: dbh, height, basal_area, agb, leaf_area
      dbh        = wood_to_dbh(allom, wood_carbon, rho, hgt_max, aboveground_frac)
      height     = dbh_to_height(allom, dbh, hgt_max)
      basal_area = pio4 * dbh * dbh
      agb        = aboveground_frac * wood_carbon
      leaf_area  = leaf_carbon * sla
   end subroutine carbon_to_structure

   !----- Wood area index [m2/m2], nplant b1wai dbh^b2wai. It follows the stem, not the    !
   !      leaf area, because the wood's heat store and exchanges scale with the stem       !
   !      (docs/science/vegetation_energy_dynamics.md, below eq. (8)). --------------------!
   elemental pure function dbh_to_wai(dbh, nplant, b1wai, b2wai) result(wai)
      real(wp), intent(in) :: dbh      !< [cm]    diameter at breast height
      real(wp), intent(in) :: nplant   !< [pl/m2] plant density
      real(wp), intent(in) :: b1wai    !< [--]    per-PFT WAI intercept
      real(wp), intent(in) :: b2wai    !< [--]    per-PFT WAI exponent
      real(wp)             :: wai      !< [m2/m2] wood area index
      wai = nplant * b1wai * max(dbh, tiny_num) ** b2wai
   end function dbh_to_wai

   !----- Sapwood share of basal area, b1sa dbh^b2sa/(pi/4 dbh^2) (Xu 2018), capped at one !
   !      since a small stem is all sapwood                                                !
   !      (docs/science/vegetation_energy_dynamics.md, below eq. (8)). --------------------!
   elemental pure function sapwood_fraction(dbh, b1sa, b2sa) result(f_sap)
      real(wp), intent(in) :: dbh      !< [cm] diameter at breast height
      real(wp), intent(in) :: b1sa     !< [--] per-PFT sapwood-area intercept
      real(wp), intent(in) :: b2sa     !< [--] per-PFT sapwood-area exponent
      real(wp)             :: f_sap    !< [--] sapwood area / basal area, in (0,1]
      real(wp) :: d, sapw_area, basal_area
      d          = max(dbh, tiny_num)
      sapw_area  = b1sa * d ** b2sa
      basal_area = pio4 * d * d
      f_sap      = min(1.0_wp, sapw_area / max(basal_area, tiny_num))
   end function sapwood_fraction

end module meds_allometry
