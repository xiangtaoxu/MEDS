!==========================================================================================!
! meds_allometry -- pan-tropical (ED2 iallom==3) size allometry, shared across the engine.  !
!                                                                                          !
! All coefficients are GLOBAL (the same for every tropical PFT, exactly as ED2's iallom==3  !
! parameterization); the only per-stem inputs are diameter, the derived height, and -- for   !
! biomass -- wood density (rho). The relationships used here (source:                        !
! ED2/ED/src/init/ed_params.f90 and utils/allometry.f90):                                   !
!                                                                                          !
!   height(dbh)        = exp(b1Ht + b2Ht*ln(dbh)),  capped at hgt_max          [m]           !
!   dbh(height)        = exp((ln(height) - b1Ht)/b2Ht)                          [cm]          !
!   crown_area(dbh,h)  = ca_b1 * (dbh^2*h)^ca_b2                                [m2]          !
!   agb(dbh,h,rho)     = agb_c1 * rho^agb_c2 * (dbh^2*h)^agb_c2                  [kgC/plant]   !
!   leaf_area(dbh,h)   = lai_b1 * (dbh^2*h)^lai_b2                              [m2/plant]    !
!                                                                                          !
! `leaf_area` is the per-stem one-sided leaf area; a cohort's LAI contribution is            !
! nplant*leaf_area. lai_b1 (0.46769540) IS the ED2 pan-tropical SLA*bleaf coefficient (the    !
! specific leaf area cancels), so MEDS gets a faithful LAI WITHOUT carrying a leaf-biomass     !
! pool. `agb` is the Chave-2014 structural-biomass form used here as the conserved carbon     !
! currency. `agb_to_dbh` inverts agb -> dbh and is the only place fusion/fission needs.       !
!==========================================================================================!
module meds_allometry
   use meds_kinds,     only : wp
   use meds_constants, only : tiny_num
   implicit none
   private

   public :: dbh_to_height, height_to_dbh, dbh_to_crown_area, dbh_to_agb, agb_to_dbh,         &
             dbh_to_leaf_area
   public :: b1Ht, b2Ht, hgt_max, agb_c1, agb_c2, ca_b1, ca_b2, lai_b1, lai_b2, light_ext

   !----- Height <-> diameter (log-linear). -----------------------------------------------!
   real(wp), parameter :: b1Ht    = 1.139963_wp     !< [--] height intercept
   real(wp), parameter :: b2Ht    = 0.564899_wp     !< [--] height slope
   real(wp), parameter :: hgt_max = 46.0_wp         !< [m]  asymptotic tropical height cap
   !----- Aboveground (structural) biomass: Chave-2014, kgC/plant. ------------------------!
   real(wp), parameter :: agb_c1  = 0.06080334_wp   !< [kgC] scale
   real(wp), parameter :: agb_c2  = 1.0044785_wp    !< [--]  exponent on rho and on dbh^2*h
   !----- Crown area, m2. -----------------------------------------------------------------!
   real(wp), parameter :: ca_b1   = 0.370_wp
   real(wp), parameter :: ca_b2   = 0.464_wp
   !----- Per-stem leaf area, m2 (= ED2 SLA*bleaf, SLA cancelled). -------------------------!
   real(wp), parameter :: lai_b1  = 0.46769540_wp
   real(wp), parameter :: lai_b2  = 0.6410495_wp
   !----- Beer-Lambert light extinction through overtopping LAI. --------------------------!
   real(wp), parameter :: light_ext = 0.5_wp

contains

   !----- Diameter -> height [m], capped at the tropical asymptote. -----------------------!
   elemental pure function dbh_to_height(dbh) result(h)
      real(wp), intent(in) :: dbh
      real(wp)             :: h
      h = min(exp(b1Ht + b2Ht * log(max(dbh, tiny_num))), hgt_max)
   end function dbh_to_height

   !----- Height -> diameter [cm] (inverse of the uncapped branch). -----------------------!
   elemental pure function height_to_dbh(h) result(dbh)
      real(wp), intent(in) :: h
      real(wp)             :: dbh
      dbh = exp((log(max(h, tiny_num)) - b1Ht) / b2Ht)
   end function height_to_dbh

   !----- Crown area [m2]. ----------------------------------------------------------------!
   elemental pure function dbh_to_crown_area(dbh, h) result(ca)
      real(wp), intent(in) :: dbh, h
      real(wp)             :: ca
      ca = ca_b1 * (dbh * dbh * h) ** ca_b2
   end function dbh_to_crown_area

   !----- Aboveground biomass [kgC/plant]. ------------------------------------------------!
   elemental pure function dbh_to_agb(dbh, h, rho) result(agb)
      real(wp), intent(in) :: dbh, h, rho
      real(wp)             :: agb
      agb = agb_c1 * rho ** agb_c2 * (dbh * dbh * h) ** agb_c2
   end function dbh_to_agb

   !----- Per-stem one-sided leaf area [m2]; cohort LAI = nplant*this. ---------------------!
   elemental pure function dbh_to_leaf_area(dbh, h) result(leaf_area)
      real(wp), intent(in) :: dbh, h
      real(wp)             :: leaf_area
      leaf_area = lai_b1 * (dbh * dbh * h) ** lai_b2
   end function dbh_to_leaf_area

   !---------------------------------------------------------------------------------------!
   ! Invert AGB -> DBH, the operation fusion/fission use to recover a diameter from the     !
   ! conserved carbon. Two regimes: below the height cap height(dbh) follows the log-linear  !
   ! law so agb = K*dbh^P (K, P below); once height saturates at hgt_max, agb scales as       !
   ! dbh^(2*agb_c2). We solve the uncapped form, then redo it in the capped form if the       !
   ! resulting stem would be taller than hgt_max. The two branches agree at the transition.   !
   !---------------------------------------------------------------------------------------!
   elemental pure function agb_to_dbh(agb, rho) result(dbh)
      real(wp), intent(in) :: agb, rho
      real(wp)             :: dbh, k_un, p_un, k_cap, p_cap, h
      k_un = agb_c1 * rho ** agb_c2 * exp(agb_c2 * b1Ht)
      p_un = agb_c2 * (2.0_wp + b2Ht)
      dbh  = (max(agb, tiny_num) / k_un) ** (1.0_wp / p_un)
      h    = exp(b1Ht + b2Ht * log(max(dbh, tiny_num)))
      if (h > hgt_max) then
         k_cap = agb_c1 * rho ** agb_c2 * hgt_max ** agb_c2
         p_cap = 2.0_wp * agb_c2
         dbh   = (max(agb, tiny_num) / k_cap) ** (1.0_wp / p_cap)
      end if
   end function agb_to_dbh

end module meds_allometry
