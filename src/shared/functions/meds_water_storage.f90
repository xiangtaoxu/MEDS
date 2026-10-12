! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_water_storage -- what each kind of water store holds, written in the unknown the soil-  !
! plant water network solves for (MEDS_SPAC_WATER_NETWORK_PLAN.md §12.1 D1).                  !
!                                                                                          !
! Newton converges well only when its unknown moves smoothly with both the water a store holds !
! and its water potential. No one variable does that for every store, so each kind picks its   !
! own:                                                                                      !
!   * a soil layer uses its effective saturation Se on the dry side, where the potential runs  !
!     off towards minus infinity while Se barely moves, and its potential on the wet side, where !
!     the retention curve is flat. The two are joined smoothly at the curve's point of largest   !
!     capacity. Above saturation the layer stores water through a small specific storage, so a  !
!     layer pushed past saturation raises its potential instead of being clipped.               !
!   * a plant tissue uses its symplastic relative water content R: the pressure-volume curve     !
!     gives a finite, smooth potential on every branch, above saturation too, and R means the    !
!     same thing for a leafless cohort whose leaf store holds nothing.                           !
!   * the pond uses its head of water.                                                         !
!                                                                                          !
! For one value of the unknown each routine returns the water potential [MPa], the water held   !
! [kg per m2 of ground], and their derivatives with respect to the unknown; a soil layer also   !
! returns its hydraulic conductivity, which follows the same curve. Nothing is floored: the      !
! solver keeps each unknown inside its range. The curves are the ones meds_water_retention uses  !
! for the diagnostic read-offs, without those read-offs' floors.                                  !
!==========================================================================================!
module meds_water_storage
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : grav_head, rho_h2o, tiny_num
   use meds_water_retention, only : SOIL_RETENTION_CAMPBELL, water_curve_t
   implicit none
   private

   public :: SOIL_SPECIFIC_STORAGE, POND_DEPRESSION_DEPTH
   public :: soil_curve_t, soil_curve_setup, soil_store, soil_unknown_from_theta
   public :: tissue_store, tissue_unknown_from_rwc, tissue_rwc_from_water
   public :: pond_store, pond_unknown_from_water

   !----- Water a saturated soil layer gains per metre of pressure head [1/m]: the compressibility !
   !      of the soil matrix and of water, a typical value for clays and loams (Freeze and Cherry   !
   !      1979, table 2.5). It matters only above saturation, where it lets a layer hold a little   !
   !      more water at a positive potential. -------------------------------------------------!
   real(wp), parameter :: SOIL_SPECIFIC_STORAGE = 1.0e-4_wp
   !----- [m] Depth over which the pond's water fills the surface's small hollows: the pond holds  !
   !      h^2/(2 delta) at a head h below it, h - delta/2 above, and nothing at a head of zero or  !
   !      less, which keeps its storage smooth when it fills or empties. -----------------------!
   real(wp), parameter :: POND_DEPRESSION_DEPTH = 1.0e-3_wp
   !----- The largest logarithm the dry-soil potential is evaluated at: e^700 m of head is far    !
   !      beyond any physical soil and still finite in double precision. -----------------------!
   real(wp), parameter :: LOG_HEAD_MAX = 700.0_wp

   !----- One soil layer's curve, set up once per layer: the retention parameters and the point  !
   !      where the unknown changes from saturation to potential. ----------------------------!
   type :: soil_curve_t
      integer(ik) :: retention = 1_ik   !< SOIL_RETENTION_VG or SOIL_RETENTION_CAMPBELL
      real(wp) :: theta_sat = 0.0_wp    !< [m3/m3]
      real(wp) :: theta_res = 0.0_wp    !< [m3/m3]
      real(wp) :: par_a     = 0.0_wp    !< vG alpha [1/m], or Campbell's air-entry potential [m, < 0]
      real(wp) :: par_n     = 0.0_wp    !< vG n [-], or Campbell's b [-]
      real(wp) :: ksat      = 0.0_wp    !< [m/s] saturated conductivity
      real(wp) :: dz        = 0.0_wp    !< [m] layer thickness
      real(wp) :: se_switch  = 1.0_wp   !< [-] saturation where the unknown changes branch
      real(wp) :: psi_switch = 0.0_wp   !< [m] the potential there
      real(wp) :: slope_switch = 1.0_wp !< [1/m] dSe/dpsi there, the unknown's scale on the wet side
   end type soil_curve_t

contains

   !=======================================================================================!
   !     Soil layers.                                                                       !
   !=======================================================================================!

   !----- Set up one layer's curve. The branch point is where the retention curve holds the   !
   !      most water per unit potential: for van Genuchten (alpha |psi|)^n = m, Se = (1+m)^-m;   !
   !      for Campbell the air-entry point, where Se reaches one. On either side of it one of    !
   !      the two unknowns is the well-scaled one. -------------------------------------------!
   pure subroutine soil_curve_setup(retention, theta_sat, theta_res, par_a, par_n, ksat, dz, c)
      integer(ik),        intent(in)  :: retention
      real(wp),           intent(in)  :: theta_sat, theta_res, par_a, par_n, ksat, dz
      type(soil_curve_t), intent(out) :: c
      real(wp) :: m, x
      c%retention = retention ; c%theta_sat = theta_sat ; c%theta_res = theta_res
      c%par_a = par_a ; c%par_n = par_n ; c%ksat = ksat ; c%dz = dz
      if (retention == SOIL_RETENTION_CAMPBELL) then
         c%se_switch    = 1.0_wp
         c%psi_switch   = par_a
         c%slope_switch = -1.0_wp / (par_n * par_a)
      else
         m = 1.0_wp - 1.0_wp / par_n
         x = m ** (1.0_wp / par_n)
         c%se_switch    = (1.0_wp + m) ** (-m)
         c%psi_switch   = -x / par_a
         c%slope_switch = par_a * m * par_n * c%se_switch * (m / (1.0_wp + m)) / x
      end if
   end subroutine soil_curve_setup

   !----- A soil layer at unknown tau: potential [MPa], water [kg/m2], conductivity [m/s], and   !
   !      their derivatives with respect to tau. Below the branch point tau is Se; above it tau   !
   !      runs linearly with the potential, at the slope that makes the two branches meet with    !
   !      the same derivative. --------------------------------------------------------------!
   elemental subroutine soil_store(c, tau, psi, dpsi, w, dw, kcond, dkcond)
      type(soil_curve_t), intent(in)  :: c
      real(wp),           intent(in)  :: tau
      real(wp),           intent(out) :: psi, dpsi    !< [MPa], [MPa per unit tau]
      real(wp),           intent(out) :: w, dw        !< [kg/m2], [kg/m2 per unit tau]
      real(wp),           intent(out) :: kcond, dkcond !< [m/s], [m/s per unit tau]
      real(wp) :: head, dhead, se, dse_dhead, dk_dse, dk_dhead, dtheta
      dtheta = c%theta_sat - c%theta_res
      if (tau <= c%se_switch) then
         !----- Dry branch: tau is the saturation itself. -----------------------------------!
         se = tau
         call head_of_saturation(c, se, head, dse_dhead)
         call conductivity_of_saturation(c, se, kcond, dk_dse)
         psi    = grav_head * head
         dpsi   = grav_head / max(dse_dhead, tiny_num)   ! clamp-ok: a slope underflowed in extreme dry soil
         w      = rho_h2o * c%dz * (c%theta_res + dtheta * se)
         dw     = rho_h2o * c%dz * dtheta
         dkcond = dk_dse
      else
         !----- Wet branch: tau measures the potential. -------------------------------------!
         head  = c%psi_switch + (tau - c%se_switch) / c%slope_switch
         dhead = 1.0_wp / c%slope_switch
         psi   = grav_head * head
         dpsi  = grav_head * dhead
         if (head < 0.0_wp) then
            call saturation_of_head(c, head, se, dse_dhead, kcond, dk_dhead)
            w      = rho_h2o * c%dz * (c%theta_res + dtheta * se)
            dw     = rho_h2o * c%dz * dtheta * dse_dhead * dhead
            dkcond = dk_dhead * dhead
         else
            w      = rho_h2o * c%dz * (c%theta_sat + SOIL_SPECIFIC_STORAGE * head)
            dw     = rho_h2o * c%dz * SOIL_SPECIFIC_STORAGE * dhead
            kcond  = c%ksat
            dkcond = 0.0_wp
         end if
      end if
   end subroutine soil_store

   !----- The unknown of a layer holding theta [m3/m3]: the inverse of soil_store's water. A     !
   !      layer at or below residual -- which an earlier model's floor could leave -- starts just  !
   !      inside the curve, at SE_START. A Campbell layer at exactly saturation starts at its      !
   !      air-entry point, the one potential its water content can name. -----------------------!
   elemental function soil_unknown_from_theta(c, theta) result(tau)
      type(soil_curve_t), intent(in) :: c
      real(wp),           intent(in) :: theta
      real(wp) :: tau
      real(wp), parameter :: SE_START = 1.0e-9_wp   ! the read-off's own lowest saturation (SE_MIN)
      real(wp) :: se, head, x, m
      se = (theta - c%theta_res) / max(c%theta_sat - c%theta_res, tiny_num)
      if (se <= c%se_switch) then
         tau = max(se, SE_START)   ! clamp-ok: a state at residual starts just inside the curve
         return
      end if
      if (theta > c%theta_sat) then
         head = (theta - c%theta_sat) / SOIL_SPECIFIC_STORAGE
      else if (c%retention == SOIL_RETENTION_CAMPBELL) then
         head = c%par_a * se ** (-c%par_n)
      else if (se >= 1.0_wp) then
         head = 0.0_wp
      else
         m = 1.0_wp - 1.0_wp / c%par_n
         x = (se ** (-1.0_wp / m) - 1.0_wp) ** (1.0_wp / c%par_n)
         head = -x / c%par_a
      end if
      tau = c%se_switch + c%slope_switch * (head - c%psi_switch)
   end function soil_unknown_from_theta

   !----- Matric head [m] at saturation se in (0, 1), and dSe/dhead [1/m]. For van Genuchten the   !
   !      head is evaluated through its logarithm, so a very dry layer gives a large but finite    !
   !      number; Se (1 - Se^(1/m)) / x is the derivative written without a cancellation. ------!
   elemental subroutine head_of_saturation(c, se, head, dse_dhead)
      type(soil_curve_t), intent(in)  :: c
      real(wp),           intent(in)  :: se
      real(wp),           intent(out) :: head, dse_dhead
      real(wp) :: m, u, logx, x
      if (c%retention == SOIL_RETENTION_CAMPBELL) then
         head      = c%par_a * se ** (-c%par_n)
         dse_dhead = -se ** (1.0_wp + c%par_n) / (c%par_n * c%par_a)
      else
         m    = 1.0_wp - 1.0_wp / c%par_n
         u    = 1.0_wp - se ** (1.0_wp / m)
         logx = (-log(se) / m + log(u)) / c%par_n
         x    = exp(min(logx, LOG_HEAD_MAX))   ! clamp-ok: keeps an extreme dry head finite
         head = -x / c%par_a
         dse_dhead = c%par_a * m * c%par_n * se * u / x
      end if
   end subroutine head_of_saturation

   !----- Saturation, dSe/dhead, conductivity and dK/dhead at a head below zero (wet branch).     !
   !      Written in x = alpha |head| so that 1 - Se^(1/m) = x^n / (1 + x^n) is exact near        !
   !      saturation, where the van Genuchten-Mualem conductivity's slope grows without bound for  !
   !      n < 2. A Campbell layer between its air-entry point and zero is saturated. -------------!
   elemental subroutine saturation_of_head(c, head, se, dse_dhead, kcond, dk_dhead)
      type(soil_curve_t), intent(in)  :: c
      real(wp),           intent(in)  :: head
      real(wp),           intent(out) :: se, dse_dhead, kcond, dk_dhead
      real(wp) :: m, x, xn, u, t, dk_dse
      if (c%retention == SOIL_RETENTION_CAMPBELL) then
         if (head < c%par_a) then
            se        = (head / c%par_a) ** (-1.0_wp / c%par_n)
            dse_dhead = -se ** (1.0_wp + c%par_n) / (c%par_n * c%par_a)
            call conductivity_of_saturation(c, se, kcond, dk_dse)
            dk_dhead  = dk_dse * dse_dhead
         else
            se = 1.0_wp ; dse_dhead = 0.0_wp ; kcond = c%ksat ; dk_dhead = 0.0_wp
         end if
         return
      end if
      m  = 1.0_wp - 1.0_wp / c%par_n
      x  = c%par_a * abs(head)
      xn = x ** c%par_n
      se = (1.0_wp + xn) ** (-m)
      u  = xn / (1.0_wp + xn)                      ! 1 - Se^(1/m)
      dse_dhead = c%par_a * m * c%par_n * se * u / max(x, tiny_num)
      t     = 1.0_wp - u ** m
      kcond = c%ksat * sqrt(se) * t * t
      !----- dK/dhead = dK/dSe dSe/dhead, the Mualem term rewritten in u and x. -----------------!
      dk_dhead = c%ksat * (0.5_wp * t * t / sqrt(se) * dse_dhead                                 &
                 + 2.0_wp * c%par_a * m * c%par_n * t * se ** (0.5_wp + 1.0_wp / m) * u ** m      &
                   / max(x, tiny_num))
   end subroutine saturation_of_head

   !----- Conductivity [m/s] and dK/dSe at saturation se in (0, 1]: Mualem-van Genuchten (pore    !
   !      connectivity 0.5) or Campbell's Se^(2b+3). No floor: it goes to zero with Se. ---------!
   elemental subroutine conductivity_of_saturation(c, se, kcond, dk_dse)
      type(soil_curve_t), intent(in)  :: c
      real(wp),           intent(in)  :: se
      real(wp),           intent(out) :: kcond, dk_dse
      real(wp) :: m, u, t
      if (c%retention == SOIL_RETENTION_CAMPBELL) then
         kcond  = c%ksat * se ** (2.0_wp * c%par_n + 3.0_wp)
         dk_dse = (2.0_wp * c%par_n + 3.0_wp) * kcond / se
      else
         m  = 1.0_wp - 1.0_wp / c%par_n
         u  = 1.0_wp - se ** (1.0_wp / m)
         t  = 1.0_wp - u ** m
         kcond  = c%ksat * sqrt(se) * t * t
         dk_dse = c%ksat * (0.5_wp * t * t / sqrt(se)                                             &
                  + 2.0_wp * sqrt(se) * t * u ** (m - 1.0_wp) * se ** (1.0_wp / m - 1.0_wp))
      end if
   end subroutine conductivity_of_saturation

   !=======================================================================================!
   !     Plant tissues.                                                                     !
   !=======================================================================================!

   !----- A tissue at symplastic relative water content r: potential [MPa], water [kg/m2 of      !
   !      ground] for nplant plants of `biomass` each, and their derivatives with respect to r.     !
   !      The pressure-volume curve without meds_water_retention's RWC floor: the solver keeps r    !
   !      above zero. Turgid above the turgor-loss point, psi = E (r - r_tlp) + pi0/r; flaccid     !
   !      below, pi0/r. -----------------------------------------------------------------------!
   elemental subroutine tissue_store(curve, biomass, nplant, r, psi, dpsi, w, dw)
      type(water_curve_t), intent(in)  :: curve
      real(wp),            intent(in)  :: biomass   !< [kgC per plant]
      real(wp),            intent(in)  :: nplant    !< [plants/m2]
      real(wp),            intent(in)  :: r         !< [-] symplastic relative water content
      real(wp),            intent(out) :: psi, dpsi, w, dw
      real(wp) :: r_tlp, w_sat
      associate (pi0 => curve%pi0, emod => curve%elastic_mod, apo => curve%apoplast_frac)
         r_tlp = (pi0 + emod) / emod
         if (r >= r_tlp) then
            psi  = emod * (r - r_tlp) + pi0 / r
            dpsi = emod - pi0 / (r * r)
         else
            psi  = pi0 / r
            dpsi = -pi0 / (r * r)
         end if
         w_sat = nplant * curve%water_sat * biomass
         w     = w_sat * ((1.0_wp - apo) * r + apo)
         dw    = w_sat * (1.0_wp - apo)
      end associate
   end subroutine tissue_store

   !----- Symplastic relative water content of a tissue holding w [kg per plant]: the inverse of  !
   !      tissue_store's water. Undefined for a tissue with no biomass, whose caller names its     !
   !      potential instead (tissue_unknown_from_rwc of rwc_from_psi). --------------------------!
   elemental function tissue_rwc_from_water(curve, biomass, w) result(r)
      type(water_curve_t), intent(in) :: curve
      real(wp),            intent(in) :: biomass, w
      real(wp) :: r, w_sat
      w_sat = curve%water_sat * biomass
      r = (w / w_sat - curve%apoplast_frac) / (1.0_wp - curve%apoplast_frac)
   end function tissue_rwc_from_water

   !----- The unknown is the relative water content itself; named for symmetry with the other   !
   !      stores, and kept strictly positive (an RWC of zero is an infinitely dry tissue). ------!
   elemental function tissue_unknown_from_rwc(r) result(tau)
      real(wp), intent(in) :: r
      real(wp) :: tau
      real(wp), parameter :: R_START = 1.0e-6_wp
      tau = max(r, R_START)   ! clamp-ok: an emptied store starts just inside the curve
   end function tissue_unknown_from_rwc

   !=======================================================================================!
   !     The pond.                                                                          !
   !=======================================================================================!

   !----- The pond at head h [m]: potential [MPa], water [kg/m2] and derivatives. Below zero head  !
   !      it holds nothing and its potential follows the soil it drains into. -----------------!
   elemental subroutine pond_store(h, psi, dpsi, w, dw)
      real(wp), intent(in)  :: h
      real(wp), intent(out) :: psi, dpsi, w, dw
      psi  = grav_head * h
      dpsi = grav_head
      if (h <= 0.0_wp) then
         w = 0.0_wp ; dw = 0.0_wp
      else if (h < POND_DEPRESSION_DEPTH) then
         w  = rho_h2o * h * h / (2.0_wp * POND_DEPRESSION_DEPTH)
         dw = rho_h2o * h / POND_DEPRESSION_DEPTH
      else
         w  = rho_h2o * (h - 0.5_wp * POND_DEPRESSION_DEPTH)
         dw = rho_h2o
      end if
   end subroutine pond_store

   !----- The head of a pond holding w [kg/m2]. An empty pond has no head to name; it returns     !
   !      zero, and the solver finds the head from the soil below. -----------------------------!
   elemental function pond_unknown_from_water(w) result(h)
      real(wp), intent(in) :: w
      real(wp) :: h, depth
      depth = w / rho_h2o
      if (depth <= 0.0_wp) then
         h = 0.0_wp
      else if (depth < 0.5_wp * POND_DEPRESSION_DEPTH) then
         h = sqrt(2.0_wp * POND_DEPRESSION_DEPTH * depth)
      else
         h = depth + 0.5_wp * POND_DEPRESSION_DEPTH
      end if
   end function pond_unknown_from_water

end module meds_water_storage
