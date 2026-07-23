!==========================================================================================!
! meds_hydr_lib -- hydraulics CONSTITUTIVE relations (material properties of tissue AND soil): !
! nonlinear pressure-volume curves (Bartlett/Tyree-Hammel), the xylem vulnerability curve and    !
! its Kirchhoff (matric flux) potential, and a fixed-grid lookup table + linear interpolant for   !
! the general-exponent Kirchhoff integral. All are stateless, elemental/pure, scalar-in kernels   !
! -- the hydraulics analogue of meds_allometry -- so they live in shared/functions and can be      !
! evaluated at config-load (table build) and on the GPU hot path alike. The stateful NETWORK        !
! SOLVER (solve_plant_water / plant_water_tendency) that assembles these into an ODE stays in        !
! src/plant/meds_plant_hydraulics; it `use`s this module.                                            !
!==========================================================================================!
module meds_hydr_lib
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   use meds_numerics,  only : gauss_legendre_7, bisect_root
   implicit none
   private

   !----- Vulnerability / Kirchhoff-conductance family. ------------------------------------!
   public :: plc_retained, dplc_dpsi, flux_potential, phi_inverse, kirchhoff_edge
   !----- Pressure-volume family. ----------------------------------------------------------!
   public :: pv_psi_tlp, pv_rwc_tlp, rwc_from_psi, psi_from_rwc
   public :: water_content, capacitance, pv_water_cap_from_traits, psi_from_water_content
   !----- Precomputed lookup table (dormant until the solver adopts it for general kexp). ---!
   public :: hydro_table_t, build_hydro_table, flux_potential_lin, kirchhoff_edge_tab
   public :: HYDRO_TABLE_NTAB, HYDRO_TABLE_RMAX
   !----- SOIL retention family (van Genuchten / Campbell): the soil-water analogue of the     !
   !      tissue PV curve -- theta(psi)/psi(theta) retention + K(theta) + C(psi) capacity, all   !
   !      elemental scalar-in kernels, so they belong beside the plant curves in this module.    !
   public :: SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   public :: soil_theta_from_psi, soil_psi_from_theta, soil_hydr_cond_from_theta, soil_moist_cap_from_psi

   !----- Retention-curve family selector codes (the soil_params_t%retention field). -------!
   integer(ik), parameter :: SOIL_RETENTION_VG       = 1_ik  !< van Genuchten-Mualem (default)
   integer(ik), parameter :: SOIL_RETENTION_CAMPBELL = 2_ik  !< Campbell / Clapp-Hornberger (option)

   !----- Numerical floors (a property of the curve regularization, not of uptake). --------!
   real(wp), parameter :: K_MIN  = 1.16e-13_wp    !< [m/s] ED2 hydcond_min (no zero-conductance lock)
   real(wp), parameter :: C_MIN  = 1.0e-9_wp      !< [1/m] specific-capacity floor (vG C -> 0 at saturation)
   real(wp), parameter :: SE_MIN = 1.0e-9_wp      !< effective-saturation floor (residual/air-dry)

   real(wp), parameter :: dpsi_eps = 1.0e-6_wp   !< |up-down| below which K_eff -> pointwise limit
   real(wp), parameter :: rwc_floor = 1.0e-4_wp   !< keep R strictly positive in the flaccid tail

   !----- Fixed-grid Kirchhoff lookup table (POD value type: trivially copyable + GPU-mappable). -!
   integer,  parameter :: HYDRO_TABLE_NTAB = 512        !< uniform intervals over [0, R_MAX]
   real(wp), parameter :: HYDRO_TABLE_RMAX = 8.0_wp     !< table domain in r = psi/psi50
   type :: hydro_table_t
      real(wp)    :: kexp  = 0.0_wp                       !< shape the table was built for (staleness guard)
      real(wp)    :: r_max = HYDRO_TABLE_RMAX             !< table domain [0, r_max]
      real(wp)    :: inv_h = 0.0_wp                       !< NTAB / r_max (uniform-grid index scale)
      integer(ik) :: n     = int(HYDRO_TABLE_NTAB, ik)    !< number of intervals
      real(wp)    :: g(0:HYDRO_TABLE_NTAB) = 0.0_wp       !< g(j) = integral_0^{j*h} du/(1+u^kexp)
   end type hydro_table_t

contains

   !----- Fraction of conductance retained (1 - PLC). psi<0, psi50<0 => r>0; clamp for psi>0. -!
   elemental real(wp) function plc_retained(psi, psi50, kexp) result(f)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r = max(psi/psi50, 0.0_wp)
      f = 1.0_wp / (1.0_wp + r**kexp)
   end function plc_retained

   !----- d(plc_retained)/d(psi); finite at psi=0 for kexp>1 (avoids the -(a/psi)f(1-f) NaN). --!
   elemental real(wp) function dplc_dpsi(psi, psi50, kexp) result(df)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r  = max(psi/psi50, 0.0_wp)
      df = -kexp * r**(kexp - 1.0_wp) / (psi50 * (1.0_wp + r**kexp)**2)
   end function dplc_dpsi

   !----- Kirchhoff (matric flux) potential Phi(psi) = integral_0^psi plc ds [MPa]. ----------!
   !      Closed form for kexp in {1,2}; fixed Gauss-Legendre on [0,r] otherwise. Phi(0)=0,     !
   !      Phi(psi<0)<0, strictly increasing in psi.                                             !
   pure real(wp) function flux_potential(psi, psi50, kexp) result(phi)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r = max(psi/psi50, 0.0_wp)                 ! r >= 0
      if (abs(kexp - 1.0_wp) < 1.0e-9_wp) then
         phi = psi50 * log(1.0_wp + r)
      else if (abs(kexp - 2.0_wp) < 1.0e-9_wp) then
         phi = psi50 * atan(r)
      else
         !----- integral_0^r du/(1+u^kexp) via the shared 7-pt Gauss-Legendre quadrature; the    !
         !       integrand is a pure internal fn carrying kexp by host association. --------------!
         phi = psi50 * gauss_legendre_7(kirchhoff_integrand, 0.0_wp, r)
      end if
   contains
      pure real(wp) function kirchhoff_integrand(u) result(y)
         real(wp), intent(in) :: u
         y = 1.0_wp / (1.0_wp + u**kexp)          ! kexp host-associated
      end function kirchhoff_integrand
   end function flux_potential

   !----- Inverse of Phi: find psi in [psi_lo, 0] with flux_potential(psi)=phi_target. --------!
   !      Monotone => robust bisection. Used by the vertical-profile diagnostic (3-node).       !
   pure real(wp) function phi_inverse(phi_target, psi50, kexp, psi_lo) result(psi)
      real(wp), intent(in) :: phi_target, psi50, kexp, psi_lo
      logical :: converged
      !----- Phi is monotone in psi, so a bracketed bisection on [psi_lo, 0] is robust. The      !
      !       residual is a pure internal fn carrying psi50/kexp/phi_target by host association. --!
      call bisect_root(phi_residual, psi_lo, 0.0_wp, 1.0e-9_wp, 60_ik, psi, converged)
   contains
      pure real(wp) function phi_residual(x) result(y)
         real(wp), intent(in) :: x
         y = flux_potential(x, psi50, kexp) - phi_target
      end function phi_residual
   end function phi_inverse

   !----- Kirchhoff edge conductance K_eff = k_cond * <plc> [kg/s/MPa]. k_cond is the maximum    !
   !      (plc=1) whole-plant/segment conductance already scaled to per-plant [kg/s/MPa]. The     !
   !      |dpsi|->0 limit is the pointwise value (L'Hopital of DeltaPhi/Deltapsi).                 !
   pure real(wp) function kirchhoff_edge(psi_up, psi_down, k_cond, psi50, kexp) result(keff)
      real(wp), intent(in) :: psi_up, psi_down, k_cond, psi50, kexp
      real(wp) :: dpsi
      dpsi = psi_up - psi_down
      if (abs(dpsi) > dpsi_eps) then
         keff = k_cond * ( flux_potential(psi_up,   psi50, kexp)                              &
                         - flux_potential(psi_down, psi50, kexp) ) / dpsi
      else
         keff = k_cond * plc_retained(0.5_wp*(psi_up + psi_down), psi50, kexp)
      end if
   end function kirchhoff_edge

   !=======================================================================================!
   !     Fixed-grid lookup table for the general-exponent Kirchhoff integral.               !
   !=======================================================================================!

   !----- Build the table G(r) = integral_0^r du/(1+u^kexp) on a uniform grid over [0, R_MAX]. --!
   !      Uses the exact 7-pt quadrature at each node (one-time cost; psi50=1 => flux_potential     !
   !      returns G(r) directly). A subroutine (not an array-returning fn) to avoid the nvfortran    !
   !      array-temp trap.                                                                            !
   pure subroutine build_hydro_table(tab, kexp)
      type(hydro_table_t), intent(out) :: tab
      real(wp),            intent(in)  :: kexp
      real(wp) :: h
      integer  :: j
      tab%kexp  = kexp
      tab%r_max = HYDRO_TABLE_RMAX
      tab%n     = int(HYDRO_TABLE_NTAB, ik)
      h         = HYDRO_TABLE_RMAX / real(HYDRO_TABLE_NTAB, wp)
      tab%inv_h = 1.0_wp / h
      do j = 0, HYDRO_TABLE_NTAB
         tab%g(j) = flux_potential(real(j, wp)*h, 1.0_wp, kexp)   ! psi50=1 => G(r) itself
      end do
   end subroutine build_hydro_table

   !----- Kirchhoff potential Phi(psi) via linear interpolation of the table [MPa]. ------------!
   !      G depends only on kexp; psi50 enters here as a runtime multiply (so a psi50 change needs   !
   !      no rebuild). r beyond r_max clamps onto the final interval (G asymptotes for kexp>1).       !
   pure real(wp) function flux_potential_lin(psi, psi50, tab) result(phi)
      real(wp),            intent(in) :: psi, psi50
      type(hydro_table_t), intent(in) :: tab
      real(wp) :: r, x, t
      integer  :: idx
      r   = max(psi/psi50, 0.0_wp)
      x   = r * tab%inv_h
      idx = int(x)
      if (idx >= HYDRO_TABLE_NTAB) idx = HYDRO_TABLE_NTAB - 1
      t   = x - real(idx, wp)
      phi = psi50 * ( tab%g(idx)*(1.0_wp - t) + tab%g(idx+1)*t )
   end function flux_potential_lin

   !----- Table-backed edge conductance: identical to kirchhoff_edge but consulting the lookup    !
   !      table (linear interp) instead of the 7-pt quadrature. The dpsi->0 pointwise limit uses     !
   !      the exact plc at the stored kexp.                                                           !
   pure real(wp) function kirchhoff_edge_tab(psi_up, psi_down, k_cond, psi50, tab) result(keff)
      real(wp),            intent(in) :: psi_up, psi_down, k_cond, psi50
      type(hydro_table_t), intent(in) :: tab
      real(wp) :: dpsi
      dpsi = psi_up - psi_down
      if (abs(dpsi) > dpsi_eps) then
         keff = k_cond * ( flux_potential_lin(psi_up,   psi50, tab)                            &
                         - flux_potential_lin(psi_down, psi50, tab) ) / dpsi
      else
         keff = k_cond * plc_retained(0.5_wp*(psi_up + psi_down), psi50, tab%kexp)
      end if
   end function kirchhoff_edge_tab

   !=======================================================================================!
   !     Pressure-volume (Bartlett / Tyree-Hammel) constitutive curves.                     !
   !=======================================================================================!

   !----- Turgor loss point [MPa] (Bartlett eqn 1). ---------------------------------------!
   elemental real(wp) function pv_psi_tlp(pi0, elastic_mod) result(psi_tlp)
      real(wp), intent(in) :: pi0, elastic_mod
      psi_tlp = pi0*elastic_mod / (pi0 + elastic_mod)
   end function pv_psi_tlp

   !----- Symplastic relative water content at turgor loss (Bartlett eqn 2). --------------!
   elemental real(wp) function pv_rwc_tlp(pi0, elastic_mod) result(rwc_tlp)
      real(wp), intent(in) :: pi0, elastic_mod
      rwc_tlp = (pi0 + elastic_mod) / elastic_mod
   end function pv_rwc_tlp

   !----- Water potential [MPa] from symplastic RWC. --------------------------------------!
   elemental real(wp) function psi_from_rwc(rwc, pi0, elastic_mod) result(psi)
      real(wp), intent(in) :: rwc, pi0, elastic_mod
      real(wp) :: rwc_tlp, r
      rwc_tlp = (pi0 + elastic_mod) / elastic_mod
      r       = max(rwc, rwc_floor)
      if (r >= rwc_tlp) then
         psi = elastic_mod*(r - rwc_tlp) + pi0/r      ! turgor + osmotic
      else
         psi = pi0/r                          ! turgor lost
      end if
   end function psi_from_rwc

   !----- Symplastic RWC from water potential (closed-form turgid inverse). ----------------!
   elemental real(wp) function rwc_from_psi(psi, pi0, elastic_mod) result(rwc)
      real(wp), intent(in) :: psi, pi0, elastic_mod
      real(wp) :: psi_tlp, b, disc
      psi_tlp = pi0*elastic_mod / (pi0 + elastic_mod)
      if (psi >= psi_tlp) then                ! turgid (psi less negative than the TLP)
         !----- elastic_mod*R^2 - (psi+elastic_mod+pi0)*R + pi0 = 0; take the + root in [rwc_tlp, 1]. ------!
         b    = psi + elastic_mod + pi0
         disc = max(b*b - 4.0_wp*elastic_mod*pi0, 0.0_wp)
         rwc  = (b + sqrt(disc)) / (2.0_wp*elastic_mod)
      else                                    ! flaccid
         rwc  = pi0/psi
      end if
      rwc = max(rwc, rwc_floor)
   end function rwc_from_psi

   !----- Total tissue water [kg, per plant] at potential psi. -----------------------------!
   !      W = W_sat_sym*R + W_apoplast, with W_sat = water_sat*biomass, W_sat_sym =           !
   !      (1-apoplast_frac)*W_sat and W_apoplast = apoplast_frac*W_sat (a constant reservoir).                       !
   elemental real(wp) function water_content(psi, pi0, elastic_mod, apoplast_frac, water_sat, biomass) result(w)
      real(wp), intent(in) :: psi, pi0, elastic_mod, apoplast_frac, water_sat, biomass
      real(wp) :: w_sat, r
      w_sat = water_sat*biomass
      r     = rwc_from_psi(psi, pi0, elastic_mod)
      w     = (1.0_wp - apoplast_frac)*w_sat*r + apoplast_frac*w_sat
   end function water_content

   !----- Capacitance C = dW/dpsi [kg/MPa, per plant] at potential psi. ---------------------!
   elemental real(wp) function capacitance(psi, pi0, elastic_mod, apoplast_frac, water_sat, biomass) result(c)
      real(wp), intent(in) :: psi, pi0, elastic_mod, apoplast_frac, water_sat, biomass
      real(wp) :: psi_tlp, w_sat_sym, r, cr
      psi_tlp   = pi0*elastic_mod / (pi0 + elastic_mod)
      w_sat_sym = (1.0_wp - apoplast_frac)*water_sat*biomass
      r         = rwc_from_psi(psi, pi0, elastic_mod)
      if (psi >= psi_tlp) then
         cr = 1.0_wp / (elastic_mod - pi0/(r*r))      ! dR/dpsi, turgid  (-> 1/(elastic_mod+|pi0|) at R=1)
      else
         cr = -(r*r)/pi0                       ! dR/dpsi, flaccid (= |pi0|/psi^2)
      end if
      c = w_sat_sym*cr
   end function capacitance

   !----- Legacy-linear mapping: the constant capacitance per biomass [kg/kgC/MPa] that a       !
   !      pure-elastic (X16) reservoir would need to match this curve's full-turgor capacitance. !
   elemental real(wp) function pv_water_cap_from_traits(pi0, elastic_mod, apoplast_frac, water_sat) result(water_cap)
      real(wp), intent(in) :: pi0, elastic_mod, apoplast_frac, water_sat
      water_cap = (1.0_wp - apoplast_frac)*water_sat / (elastic_mod + abs(pi0))
   end function pv_water_cap_from_traits

   !----- Exact inverse of water_content: water potential [MPa] from total tissue water [kg,      !
   !      per plant] (MEDS_ED2_RK45_DESIGN.md sec 4 -- diagnoses psi from the prognostic mass      !
   !      state). Recovers the symplastic RWC by removing the constant apoplastic reservoir, then   !
   !      composes the existing exact psi_from_rwc inverse -- every constituent already exists,     !
   !      this only inverts water_content's own forward map (w = (1-apo)*w_sat*rwc + apo*w_sat).    !
   elemental real(wp) function psi_from_water_content(w, pi0, elastic_mod, apoplast_frac, water_sat, biomass) &
                              result(psi)
      real(wp), intent(in) :: w, pi0, elastic_mod, apoplast_frac, water_sat, biomass
      real(wp) :: w_sat, rwc
      w_sat = water_sat*biomass
      rwc   = (w - apoplast_frac*w_sat) / max((1.0_wp - apoplast_frac)*w_sat, tiny_num)
      psi   = psi_from_rwc(rwc, pi0, elastic_mod)
   end function psi_from_water_content

   !=======================================================================================!
   !  SOIL retention curves (van Genuchten-Mualem default / Campbell-Clapp-Hornberger option). !
   !  Closed-form theta(psi)/psi(theta) inverses + K(theta) + C(psi), all elemental, branch-    !
   !  light and FPE-safe under -fpe0 / -Ktrap=fp (Se clamped, K floored, C guarded). Coupling to !
   !  the plant hydraulics kernel is through the POTENTIAL psi (curve-independent), so the curve  !
   !  is a free run-time choice. par_a/par_n carry {alpha,n} for vG or {psi_sat,b} for Campbell.  !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! theta(psi): volumetric water content from matric potential [m, <= 0].                  !
   !---------------------------------------------------------------------------------------!
   elemental function soil_theta_from_psi(retention, psi, theta_sat, theta_res, par_a, par_n)  &
                      result(theta)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: psi, theta_sat, theta_res, par_a, par_n
      real(wp)                :: theta, se, m, ah
      if (psi >= 0.0_wp) then
         theta = theta_sat
         return
      end if
      if (retention == SOIL_RETENTION_CAMPBELL) then
         !----- par_a = psi_sat (< 0), par_n = b. --------------------------------------!
         se = min(1.0_wp, (psi / par_a) ** (-1.0_wp / par_n))
      else
         !----- van Genuchten: par_a = alpha [1/m], par_n = n. -------------------------!
         m  = 1.0_wp - 1.0_wp / par_n
         ah = par_a * abs(psi)
         se = (1.0_wp + ah ** par_n) ** (-m)
      end if
      se    = min(max(se, SE_MIN), 1.0_wp)
      theta = theta_res + (theta_sat - theta_res) * se
   end function soil_theta_from_psi

   !---------------------------------------------------------------------------------------!
   ! psi(theta): matric potential [m, <= 0] from water content -- closed-form inverse.      !
   !---------------------------------------------------------------------------------------!
   elemental function soil_psi_from_theta(retention, theta, theta_sat, theta_res, par_a, par_n) &
                      result(psi)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: theta, theta_sat, theta_res, par_a, par_n
      real(wp)                :: psi, se, m
      se = (theta - theta_res) / max(theta_sat - theta_res, tiny_num)
      se = min(max(se, SE_MIN), 1.0_wp)
      if (se >= 1.0_wp) then
         psi = 0.0_wp
         return
      end if
      if (retention == SOIL_RETENTION_CAMPBELL) then
         psi = par_a * se ** (-par_n)                       ! par_a = psi_sat, par_n = b
      else
         m   = 1.0_wp - 1.0_wp / par_n
         psi = -(1.0_wp / par_a) * (se ** (-1.0_wp / m) - 1.0_wp) ** (1.0_wp / par_n)
      end if
   end function soil_psi_from_theta

   !---------------------------------------------------------------------------------------!
   ! K(theta): unsaturated hydraulic conductivity [m/s], floored at K_MIN.                  !
   !---------------------------------------------------------------------------------------!
   elemental function soil_hydr_cond_from_theta(retention, theta, theta_sat, theta_res, par_a, par_n,    &
                      ksat) result(kcond)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: theta, theta_sat, theta_res, par_a, par_n, ksat
      real(wp)                :: kcond, se, m, tmp
      se = (theta - theta_res) / max(theta_sat - theta_res, tiny_num)
      se = min(max(se, SE_MIN), 1.0_wp)
      if (retention == SOIL_RETENTION_CAMPBELL) then
         kcond = ksat * se ** (2.0_wp * par_n + 3.0_wp)     ! par_n = b
      else
         m     = 1.0_wp - 1.0_wp / par_n
         tmp   = 1.0_wp - (1.0_wp - se ** (1.0_wp / m)) ** m
         kcond = ksat * sqrt(se) * tmp * tmp                ! l = 0.5 Mualem pore-connectivity
      end if
      kcond = max(kcond, K_MIN)
   end function soil_hydr_cond_from_theta

   !---------------------------------------------------------------------------------------!
   ! C(psi) = dtheta/dpsi: specific moisture capacity [1/m], >= 0, floored at C_MIN.        !
   !---------------------------------------------------------------------------------------!
   elemental function soil_moist_cap_from_psi(retention, psi, theta_sat, theta_res, par_a, par_n)      &
                      result(cap)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: psi, theta_sat, theta_res, par_a, par_n
      real(wp)                :: cap, m, ah, dtr
      dtr = theta_sat - theta_res
      if (psi >= 0.0_wp) then
         cap = C_MIN
         return
      end if
      if (retention == SOIL_RETENTION_CAMPBELL) then
         !----- C = -(theta_sat-theta_res)/(b*psi_sat) * (psi/psi_sat)^(-1/b - 1). ------!
         cap = -dtr / (par_n * par_a) * (psi / par_a) ** (-1.0_wp / par_n - 1.0_wp)
      else
         m   = 1.0_wp - 1.0_wp / par_n
         ah  = par_a * abs(psi)
         cap = par_a * m * par_n * dtr * ah ** (par_n - 1.0_wp)                            &
               * (1.0_wp + ah ** par_n) ** (-m - 1.0_wp)
      end if
      cap = max(cap, C_MIN)
   end function soil_moist_cap_from_psi

end module meds_hydr_lib
