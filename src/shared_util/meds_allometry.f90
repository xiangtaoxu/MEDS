!==========================================================================================!
! meds_allometry -- pan-tropical (ED2 iallom==3) size allometry, shared across the engine.  !
!                                                                                          !
! All coefficients are GLOBAL (the same for every tropical PFT, exactly as ED2's iallom==3  !
! parameterization); the only per-stem inputs are diameter, the derived height, and -- for   !
! biomass -- wood density (rho). The relationships used here (source:                        !
! ED2/ED/src/init/ed_params.f90 and utils/allometry.f90):                                   !
!                                                                                          !
!   height(dbh)        = exp(b1Ht + b2Ht*ln(dbh)),  capped at height_max          [m]           !
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
   public :: b1Ht, b2Ht, height_max, agb_c1, agb_c2, ca_b1, ca_b2, lai_b1, lai_b2, light_ext
   public :: set_allometry

   !----- Allometry coefficients are RUNTIME CONFIGURATION, not hard-coded: they are set once  !
   !       at config load by set_allometry (from the PFT config file) and read thereafter.      !
   !       `protected` => read-only outside this module; the pan-tropical (ED2 iallom==3)        !
   !       values are the canonical defaults shipped in the config, not baked into the source.  !
   !       The host allometry functions read these directly; the offloaded growth kernel takes  !
   !       them as scalar arguments instead (it cannot read host module state on the device).   !
   real(wp), protected :: b1Ht       !< [--]  height <-> diameter intercept
   real(wp), protected :: b2Ht       !< [--]  height <-> diameter slope
   real(wp), protected :: height_max !< [m]   asymptotic tropical height cap
   real(wp), protected :: agb_c1     !< [kgC] AGB scale (Chave-2014)
   real(wp), protected :: agb_c2     !< [--]  AGB exponent on rho and on dbh^2*h
   real(wp), protected :: ca_b1      !< [--]  crown-area scale
   real(wp), protected :: ca_b2      !< [--]  crown-area exponent
   real(wp), protected :: lai_b1     !< [--]  per-stem leaf-area scale (= ED2 SLA*bleaf)
   real(wp), protected :: lai_b2     !< [--]  per-stem leaf-area exponent
   real(wp), protected :: light_ext  !< [--]  Beer-Lambert extinction through overtopping LAI

contains

   !----- Install the allometry coefficients (called once at config load). ----------------!
   subroutine set_allometry(b1Ht_in, b2Ht_in, height_max_in, agb_c1_in, agb_c2_in,           &
                            ca_b1_in, ca_b2_in, lai_b1_in, lai_b2_in, light_ext_in)
      real(wp), intent(in) :: b1Ht_in, b2Ht_in, height_max_in, agb_c1_in, agb_c2_in
      real(wp), intent(in) :: ca_b1_in, ca_b2_in, lai_b1_in, lai_b2_in, light_ext_in
      b1Ht = b1Ht_in ; b2Ht = b2Ht_in ; height_max = height_max_in
      agb_c1 = agb_c1_in ; agb_c2 = agb_c2_in
      ca_b1 = ca_b1_in ; ca_b2 = ca_b2_in
      lai_b1 = lai_b1_in ; lai_b2 = lai_b2_in ; light_ext = light_ext_in
   end subroutine set_allometry

   !----- Diameter -> height [m], capped at the tropical asymptote. -----------------------!
   elemental pure function dbh_to_height(dbh) result(h)
      real(wp), intent(in) :: dbh
      real(wp)             :: h
      h = min(exp(b1Ht + b2Ht * log(max(dbh, tiny_num))), height_max)
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
   ! law so agb = K*dbh^P (K, P below); once height saturates at height_max, agb scales as       !
   ! dbh^(2*agb_c2). We solve the uncapped form, then redo it in the capped form if the       !
   ! resulting stem would be taller than height_max. The two branches agree at the transition.   !
   !---------------------------------------------------------------------------------------!
   elemental pure function agb_to_dbh(agb, rho) result(dbh)
      real(wp), intent(in) :: agb, rho
      real(wp)             :: dbh, k_un, p_un, k_cap, p_cap, h
      k_un = agb_c1 * rho ** agb_c2 * exp(agb_c2 * b1Ht)
      p_un = agb_c2 * (2.0_wp + b2Ht)
      dbh  = (max(agb, tiny_num) / k_un) ** (1.0_wp / p_un)
      h    = exp(b1Ht + b2Ht * log(max(dbh, tiny_num)))
      if (h > height_max) then
         k_cap = agb_c1 * rho ** agb_c2 * height_max ** agb_c2
         p_cap = 2.0_wp * agb_c2
         dbh   = (max(agb, tiny_num) / k_cap) ** (1.0_wp / p_cap)
      end if
   end function agb_to_dbh

end module meds_allometry
