! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_allometry -- pan-tropical (ED2 iallom==3) size allometry, shared across the engine.  !
!                                                                                          !
! The coefficients are one record, allometry_t, that every routine here takes as its first  !
! argument. A run's record is cfg%allom; each site carries a copy on its cohort block        !
! (cohort%allom), so the polygons of a region can each have their own (polygon runtime plan  !
! §14, RP0). The coefficients are the same for every PFT, exactly as ED2's iallom==3; the    !
! only per-stem inputs are diameter, the derived height, and -- for biomass -- wood density  !
! (rho). The relationships (source: ED2/ED/src/init/ed_params.f90 and utils/allometry.f90):  !
!                                                                                          !
!   height(dbh)        = exp(b1Ht + b2Ht*ln(dbh)),  capped at hgt_max (per-PFT arg) [m]           !
!   dbh(height)        = exp((ln(height) - b1Ht)/b2Ht)                          [cm]          !
!   ... or, with height_form = HEIGHT_GMM, a saturating generalized Michaelis-Menten height:   !
!   height(dbh)        = gmm_a*dbh^gmm_b / (gmm_k + dbh^gmm_b), capped at hgt_max     [m]       !
!   dbh(height)        = (gmm_k*height / (gmm_a - height))^(1/gmm_b)                   [cm]      !
!   (Cano et al. 2019, Biogeosciences 16:847, fitted it to the Barro Colorado Nature Monument.)  !
!   ... or, with height_form = HEIGHT_CHAVE, Chave et al. (2014) eq. 6a, the pantropical curve with !
!   its site stress index E (chave_e: dry, seasonal or hot sites have shorter trees):            !
!   ln height(dbh)     = 0.893 - E + 0.760 ln(dbh) - 0.0340 ln(dbh)^2 + 0.243^2/2, capped    [m]  !
!   The 0.243^2/2 turns the fit's log-scale median into the mean, as BIOMASS's retrieveH does.   !
!   dbh(height) is the smaller root of that quadratic in ln(dbh), the rising branch.            !
!   local_hmax         = height(HMAX_DBH = 100 cm), capped: the local maximum height    [m]       !
!   crown_area(dbh,h)  = ca_b1 * (dbh^2*h)^ca_b2                                [m2]          !
!   agb(dbh,h,rho)     = agb_c1 * rho^agb_c2 * (dbh^2*h)^agb_c2                  [kgC/plant]   !
!   leaf_area(dbh,h)   = lai_b1 * (dbh^2*h)^lai_b2                              [m2/plant]    !
!                                                                                          !
! `leaf_area` is the per-stem one-sided leaf area; a cohort's LAI contribution is            !
! nplant*leaf_area. `agb` is the conserved carbon currency. `agb_to_dbh` inverts agb -> dbh    !
! and is the only place fusion/fission needs.                                                !
!                                                                                          !
! The default coefficients, all in carbon at the carbon fraction C2B = 2:                    !
!   agb    Chave et al. (2014, Glob. Change Biol.) eq. 4, AGB = 0.0673*(rho*dbh^2*h)^0.976   !
!          in kg dry mass, so agb_c1 = 0.0673/2 = 0.03365 and agb_c2 = 0.976. The MEDS form    !
!          agb_c1*rho^agb_c2*(dbh^2*h)^agb_c2 is Chave's exactly.                             !
!   leaf   ED2's pan-tropical BAAD fit as ED2 applies it: size2bl uses                        !
!          c14f15_bl_xx(1)/(SLA*C2B), so lai_b1 = 0.46769540/2 = 0.23384770 (= ED2 SLA*bleaf; !
!          the specific leaf area cancels) and lai_b2 = 0.6410495.                            !
! Taken without the division, AGB is dry mass labelled as carbon and leaf area twice ED2's.  !
!==========================================================================================!
module meds_allometry
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num, pio4
   implicit none
   private

   public :: allometry_t
   public :: dbh_to_height, height_to_dbh, dbh_to_crown_area, dbh_to_agb, agb_to_dbh,         &
             dbh_to_leaf_area, local_hmax
   public :: size2leaf_carbon, size2wood_carbon, wood_to_dbh, carbon_to_structure, min_cohort_carbon
   public :: storage_carbon_target
   public :: dbh_to_wai, sapwood_fraction
   public :: HEIGHT_POWER, HEIGHT_GMM, HEIGHT_CHAVE

   !----- The height curve: the pan-tropical power law, a saturating gMM, or Chave (2014). ---!
   integer(ik), parameter :: HEIGHT_POWER = 1_ik, HEIGHT_GMM = 2_ik, HEIGHT_CHAVE = 3_ik
   !----- Chave et al. (2014) eq. 6a: ln H = CHAVE_A - E + CHAVE_B ln D + CHAVE_C (ln D)^2, and the  !
   !       residual standard error that turns its median into the mean (RSE^2 / 2). --------------!
   real(wp),    parameter :: CHAVE_A   = 0.893_wp
   real(wp),    parameter :: CHAVE_B   = 0.760_wp
   real(wp),    parameter :: CHAVE_C   = -0.0340_wp
   real(wp),    parameter :: CHAVE_RSE = 0.243_wp
   !----- The local maximum height is the height of a tree this wide on the local curve: one     !
   !       definition for every curve, saturating or not (Chave's has no asymptote). -----------!
   real(wp),    parameter :: HMAX_DBH  = 100.0_wp     !< [cm]
   !----- Newton's method for agb -> dbh under the gMM height (agb_to_dbh_gmm). -------------!
   integer,     parameter :: NEWTON_MAX_ITER = 50
   real(wp),    parameter :: NEWTON_TOL      = 1.0e-12_wp    !< [ln cm] change in ln(dbh)

   !----- The allometry coefficients. A run's record is cfg%allom, filled from the PFT file's     !
   !       [allometry] block; a site carries a copy (cohort%allom). The initializers are the shipped !
   !       pan-tropical defaults, the numbers meds_config_pft.toml ships, so a record made without a  !
   !       config (a unit test, a Python call) is correct rather than zero. -----------------------!
   type :: allometry_t
      real(wp)    :: b1Ht        = 1.139963_wp    !< [--]  height <-> diameter intercept
      real(wp)    :: b2Ht        = 0.564899_wp    !< [--]  height <-> diameter slope
      real(wp)    :: agb_c1      = 0.03365_wp     !< [kgC] AGB scale (Chave 2014 eq. 4 / C2B)
      real(wp)    :: agb_c2      = 0.976_wp       !< [--]  AGB exponent on rho and on dbh^2*h (Chave 2014)
      real(wp)    :: ca_b1       = 0.370_wp       !< [--]  crown-area scale
      real(wp)    :: ca_b2       = 0.464_wp       !< [--]  crown-area exponent
      real(wp)    :: lai_b1      = 0.23384770_wp  !< [--]  per-stem leaf-area scale (= ED2 SLA*bleaf)
      real(wp)    :: lai_b2      = 0.6410495_wp   !< [--]  per-stem leaf-area exponent
      real(wp)    :: light_ext   = 0.5_wp         !< [--]  Beer-Lambert extinction through overtopping LAI
      integer(ik) :: height_form = HEIGHT_POWER   !< the height curve: power law, gMM or Chave (2014)
      !----- gMM coefficients, read only when height_form = HEIGHT_GMM: Cano et al. (2019) eq. 7. !
      real(wp)    :: gmm_a       = 58.0_wp        !< [m]   asymptotic height
      real(wp)    :: gmm_b       = 0.73_wp        !< [--]  exponent on dbh
      real(wp)    :: gmm_k       = 21.8_wp        !< [cm^gmm_b] half-saturation
      !----- Chave's stress index E, read only when height_form = HEIGHT_CHAVE. ----------------!
      real(wp)    :: chave_e     = 0.0_wp         !< [--]  site stress index: the larger, the shorter the trees
   end type allometry_t

contains

   !----- Diameter -> height [m], capped at the per-PFT asymptote hgt_max. -----------------!
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

   !----- The local maximum height [m]: a HMAX_DBH-wide tree's height on the local curve, capped  !
   !       at hgt_max. Chave's E scales every height by the same factor, so a fixed share of it is    !
   !       reached at the same diameter at every site, at a height that follows the local stature. !
   elemental pure function local_hmax(allom, hgt_max) result(h)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: hgt_max
      real(wp)                      :: h
      h = dbh_to_height(allom, HMAX_DBH, hgt_max)
   end function local_hmax

   !----- Height -> diameter [cm] (inverse of the uncapped curve). -----------------------!
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

   !----- Chave's curve: ln H = chave_intercept + CHAVE_B x + CHAVE_C x^2 with x = ln(dbh). ---!
   pure real(wp) function chave_intercept(allom)
      type(allometry_t), intent(in) :: allom
      chave_intercept = CHAVE_A - allom%chave_e + 0.5_wp * CHAVE_RSE * CHAVE_RSE
   end function chave_intercept

   !----- The smaller root x of CHAVE_C x^2 + b x + c = 0: the curve's rising branch, which holds  !
   !       up to dbh = exp(-CHAVE_B / (2 CHAVE_C)), about 70,000 cm. Written as 2(-c)/(b + sqrt(.)) !
   !       so it does not lose digits when c is small. -----------------------------------------!
   pure real(wp) function chave_small_root(b, c) result(x)
      real(wp), intent(in) :: b, c
      x = -2.0_wp * c / (b + sqrt(max(b * b - 4.0_wp * CHAVE_C * c, 0.0_wp)))
   end function chave_small_root

   !----- Crown area [m2]. ----------------------------------------------------------------!
   elemental pure function dbh_to_crown_area(allom, dbh, h) result(ca)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h
      real(wp)                      :: ca
      ca = allom%ca_b1 * (dbh * dbh * h) ** allom%ca_b2
   end function dbh_to_crown_area

   !----- Aboveground biomass [kgC/plant]. ------------------------------------------------!
   elemental pure function dbh_to_agb(allom, dbh, h, rho) result(agb)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h, rho
      real(wp)                      :: agb
      agb = allom%agb_c1 * rho ** allom%agb_c2 * (dbh * dbh * h) ** allom%agb_c2
   end function dbh_to_agb

   !----- AGB carbon [kgC/plant] of a minimum-size (recruit) plant at min_cohort_height: the    !
   !       recruit "unit" carbon used to convert a reproduction-carbon flux into recruit numbers. !
   elemental pure function min_cohort_carbon(allom, min_cohort_height, rho) result(carbon_min)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: min_cohort_height, rho
      real(wp)                      :: carbon_min
      carbon_min = dbh_to_agb(allom, height_to_dbh(allom, min_cohort_height), min_cohort_height, rho)
   end function min_cohort_carbon

   !----- Per-stem one-sided leaf area [m2]; cohort LAI = nplant*this. ---------------------!
   elemental pure function dbh_to_leaf_area(allom, dbh, h) result(leaf_area)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h
      real(wp)                      :: leaf_area
      leaf_area = allom%lai_b1 * (dbh * dbh * h) ** allom%lai_b2
   end function dbh_to_leaf_area

   !---------------------------------------------------------------------------------------!
   ! Invert AGB -> DBH, the operation fusion/fission use to recover a diameter from the     !
   ! conserved carbon. Two regimes: below the height cap height(dbh) follows the log-linear  !
   ! law so agb = K*dbh^P (K, P below); once height saturates at hgt_max, agb scales as          !
   ! dbh^(2*agb_c2). We solve the uncapped form, then redo it in the capped form if the       !
   ! resulting stem would be taller than hgt_max. The two branches agree at the transition.      !
   !---------------------------------------------------------------------------------------!
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

   !---------------------------------------------------------------------------------------!
   ! agb -> dbh under the gMM height, which has no closed-form inverse. With x = ln(dbh), the  !
   ! biomass law reads 2x + ln h(x) = ln(agb / (agb_c1 rho^agb_c2)) / agb_c2. The left side     !
   ! rises with x at a slope between 2 and 2 + gmm_b and is concave (the hgt_max cap only        !
   ! flattens it), so Newton's method lands at or below the root after its first step and then  !
   ! climbs to it without overshooting.                                                      !
   !---------------------------------------------------------------------------------------!
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

   !---------------------------------------------------------------------------------------!
   ! agb -> dbh under Chave's height, in closed form. With x = ln(dbh) the biomass law reads    !
   ! 2x + ln h(x) = T, T = ln(agb / (agb_c1 rho^agb_c2)) / agb_c2; below the cap ln h(x) is       !
   ! quadratic in x, so x is the smaller root of CHAVE_C x^2 + (2 + CHAVE_B) x + a - T = 0. A stem !
   ! that root would make taller than hgt_max has 2x + ln(hgt_max) = T instead.                 !
   !---------------------------------------------------------------------------------------!
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

   !---------------------------------------------------------------------------------------!
   ! Carbon-pool size targets for the daily allocator (meds_plant_carbon_allocation).        !
   ! All are CARBON [kgC/plant] and are thin, EXACT re-expressions of the size allometry      !
   ! above, so the two descriptions of a plant cannot disagree.                                !
   ! All-carbon: SLA / density conversions are the caller's traits, done                       !
   ! once at parameter init.                                                                  !
   !                                                                                          !
   !   * size2leaf_carbon = leaf_area / SLA  -- UN-FOLDS the SLA that lai_b1 folds in, so       !
   !     leaf_area = leaf_carbon*SLA exactly and LAI = nplant*leaf_carbon*SLA is unchanged.     !
   !   * size2wood_carbon = agb / aboveground_frac -- TOTAL woody carbon (incl. the belowground  !
   !     coarse fraction); its aboveground share (aboveground_frac * wood_carbon) is EXACTLY     !
   !     the Chave dbh_to_agb, so the AGB currency is reproduced with no calibration.            !
   !   * wood_to_dbh = agb_to_dbh(wood_carbon * aboveground_frac) -- the analytic inverse of     !
   !     size2wood_carbon (reuses the capped-height inversion). wood_carbon is the carbon-        !
   !     prognostic SIZE ANCHOR: dbh is DERIVED from it (the mirror of dbh -> agb today).         !
   !---------------------------------------------------------------------------------------!
   elemental pure function size2leaf_carbon(allom, dbh, h, sla) result(leaf_carbon)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h, sla        !< [cm],[m],[m2/kgC]
      real(wp)                      :: leaf_carbon         !< [kgC/plant]
      leaf_carbon = dbh_to_leaf_area(allom, dbh, h) / max(sla, tiny_num)
   end function size2leaf_carbon

   !----- A cohort's storage target: storage_cushion times the carbon one full flush rebuilds, its    !
   !      leaves plus its fine roots. It is also the starting storage of a cohort set on allometry. --!
   elemental pure function storage_carbon_target(storage_cushion, leaf_carbon, fineroot_carbon) result(storage)
      real(wp), intent(in) :: storage_cushion              !< [--]
      real(wp), intent(in) :: leaf_carbon, fineroot_carbon !< [kgC/plant]
      real(wp)             :: storage                      !< [kgC/plant]
      storage = storage_cushion * (leaf_carbon + fineroot_carbon)
   end function storage_carbon_target

   elemental pure function size2wood_carbon(allom, dbh, h, rho, aboveground_frac) result(wood_carbon)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: dbh, h, rho, aboveground_frac
      real(wp)                      :: wood_carbon         !< [kgC/plant] total woody (above + belowground)
      wood_carbon = dbh_to_agb(allom, dbh, h, rho) / max(aboveground_frac, tiny_num)
   end function size2wood_carbon

   elemental pure function wood_to_dbh(allom, wood_carbon, rho, hgt_max, aboveground_frac) result(dbh)
      type(allometry_t), intent(in) :: allom
      real(wp),          intent(in) :: wood_carbon, rho, hgt_max, aboveground_frac
      real(wp)                      :: dbh                 !< [cm] derived from the total woody carbon
      dbh = agb_to_dbh(allom, wood_carbon * aboveground_frac, rho, hgt_max)
   end function wood_to_dbh

   !---------------------------------------------------------------------------------------!
   ! CARBON -> full cached geometry (the carbon-prognostic flip): wood_carbon is the size    !
   ! anchor, so dbh = wood_to_dbh(wood_carbon), then height/basal_area/agb follow, and        !
   ! leaf_area comes straight from the prognostic leaf_carbon. Consolidates the carbon        !
   ! allometry (wood_to_dbh + size2*carbon above) into one composite; it IS the body of the   !
   ! former set_cohort_size_from_carbon, lifted here so the engine + plant share one path.    !
   !---------------------------------------------------------------------------------------!
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

   !---------------------------------------------------------------------------------------!
   ! WOOD AREA INDEX from stem size (ED2 allometry.f90:1475, `cpatch%wai`).                  !
   !                                                                                        !
   !     wai = nplant * b1WAI * dbh ** b2WAI                                                 !
   !                                                                                        !
   ! WAI is a STEM-SIZE relation, which is the whole point: the fast loop used to set        !
   ! wai = 0.20 * lai, tying wood area to LEAF area. That is what made the modelled wood     !
   ! thermal timescale ~6-20x too short, because tau_wood ~ (wood mass)/(wood AREA) and a    !
   ! leaf-tied area makes the ratio nearly size-independent. It also mis-scales the wood     !
   ! boundary layer, the wood longwave emission area and the wood sensible-heat coefficient, !
   ! all of which take WAI directly.                                                          !
   !---------------------------------------------------------------------------------------!
   elemental pure function dbh_to_wai(dbh, nplant, b1wai, b2wai) result(wai)
      real(wp), intent(in) :: dbh      !< [cm]    diameter at breast height
      real(wp), intent(in) :: nplant   !< [pl/m2] plant density
      real(wp), intent(in) :: b1wai    !< [--]    per-PFT WAI intercept
      real(wp), intent(in) :: b2wai    !< [--]    per-PFT WAI exponent
      real(wp)             :: wai      !< [m2/m2] wood area index
      wai = nplant * b1wai * max(dbh, tiny_num) ** b2wai
   end function dbh_to_wai

   !---------------------------------------------------------------------------------------!
   ! SAPWOOD FRACTION of basal area (ED2 allometry.f90:202 `dbh2sf`, Xu 2018).               !
   !                                                                                        !
   !     f_sap = min(1, b1SA * dbh**b2SA / (pi/4 * dbh^2))                                   !
   !                                                                                        !
   ! Capped at 1 because a small stem is sapwood all the way through -- which is also why    !
   ! this is a defensible proxy for the THERMALLY ACTIVE wood mass, not just the hydraulic   !
   ! one. The thermal store wants the wood within a diurnal damping depth of the surface     !
   ! (~4.5 cm in wet wood); the sapwood ring is ~2-5 cm on a mature bole and becomes the      !
   ! whole stem on a small one, so the two quantities are comparable at large size and       !
   ! converge exactly at small size. Documented as a proxy rather than an identity -- a      !
   ! separate thermally-active mass would additionally need bole/branch partitioning, which  !
   ! MEDS does not carry.                                                                     !
   !---------------------------------------------------------------------------------------!
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
