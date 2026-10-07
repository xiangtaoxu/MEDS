! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_hydr_lib -- how water MOVES through the plant: the xylem vulnerability curve and its      !
! Kirchhoff (matric flux) potential, a fixed-grid lookup table + linear interpolant for the        !
! general-exponent Kirchhoff integral, the sapwood conductivity of a wood density, and the root    !
! profile. How much water a tissue or soil layer holds is in meds_water_retention. All are         !
! stateless, elemental/pure, scalar-in kernels                                                     !
! -- the hydraulics analogue of meds_allometry -- so they live in src/shared/functions and can be  !
! evaluated at config-load (table build) and on the GPU hot path alike. The stateful NETWORK        !
! SOLVER (solve_plant_water / plant_water_tendency) that assembles these into an ODE stays in        !
! src/fast_dynamics/plant/meds_plant_hydraulics; it `use`s this module.                              !
!==========================================================================================!
module meds_hydr_lib
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   implicit none
   private

   !----- Vulnerability / Kirchhoff-conductance family. ------------------------------------!
   public :: plc_retained, dplc_dpsi, flux_potential, kirchhoff_edge
   public :: root_fraction_profile, wood_kmax_from_density
   !----- Precomputed lookup table (dormant until the solver adopts it for general kexp). ---!
   public :: hydro_table_t, build_hydro_table, flux_potential_lin, kirchhoff_edge_tab
   public :: HYDRO_TABLE_NTAB, HYDRO_TABLE_RMAX
   real(wp), parameter :: dpsi_eps = 1.0e-6_wp   !< |up-down| below which K_eff -> pointwise limit

   !----- Sapwood specific conductivity against wood density (Xu et al. 2016, New Phytol 212:80;   !
   !      ED2 plant_hydro_scheme 2): ln Ks = a + b rho, residual variance s2, fitted to 72 species  !
   !      of neotropical seasonally dry forests (R2 0.21). Panama canopy branches give 2-4 at 0.6. -!
   real(wp), parameter :: kmax_ln_a     =  2.348_wp   !< [ln(kg/m/s/MPa)] intercept
   real(wp), parameter :: kmax_ln_b     = -2.455_wp   !< [ln(kg/m/s/MPa) per g/cm3] slope
   real(wp), parameter :: kmax_ln_var   =  0.6186_wp  !< [-] residual variance of ln Ks
   real(wp), parameter :: kmax_rho_min  =  0.35_wp    !< [g/cm3] wood densities outside the fitted
   real(wp), parameter :: kmax_rho_max  =  0.95_wp    !<   range take its edge (ED2's bounds)

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

   !----- Sapwood specific conductivity [kg/m/s/MPa] of a wood density rho [g/cm3]: the mean of   !
   !      the lognormal fit, exp(a + b rho + s2/2) -- 6.0 at rho 0.35, 3.3 at 0.60, 1.4 at 0.95. -!
   elemental real(wp) function wood_kmax_from_density(rho) result(kmax)
      real(wp), intent(in) :: rho
      kmax = exp(kmax_ln_a + kmax_ln_b*min(max(rho, kmax_rho_min), kmax_rho_max) + 0.5_wp*kmax_ln_var)
   end function wood_kmax_from_density

   !----- Fraction of conductance retained (1 - PLC). psi<0, psi50<0 => r>0; clamp for psi>0.    !
   !      r=0 is guarded explicitly (0**kexp=0 for any kexp>0): nvfortran's real**real codegen      !
   !      for a general (non-integer-recognized) exponent routes through log(u), and log(0) trips   !
   !      its strict -Ktrap=fp even though the OVERALL mathematical result is well-defined. The      !
   !      Kirchhoff integrand below guards u=0 for the same reason. -------------------------------!
   elemental real(wp) function plc_retained(psi, psi50, kexp) result(f)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r = max(psi/psi50, 0.0_wp)
      if (r <= 0.0_wp) then
         f = 1.0_wp
      else
         f = 1.0_wp / (1.0_wp + r**kexp)
      end if
   end function plc_retained

   !----- d(plc_retained)/d(psi); finite at psi=0 for kexp>1 (avoids the -(a/psi)f(1-f) NaN). --!
   !      NOTE: shares plc_retained's r=0 real**real hazard (r**(kexp-1.0) is a further 0**0 case  !
   !      when kexp=1) but is currently unused by any caller, so it is deliberately left as-is --   !
   !      guard it the same way before giving it a real call site. --------------------------------!
   elemental real(wp) function dplc_dpsi(psi, psi50, kexp) result(df)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r  = max(psi/psi50, 0.0_wp)
      df = -kexp * r**(kexp - 1.0_wp) / (psi50 * (1.0_wp + r**kexp)**2)
   end function dplc_dpsi

   !----- Kirchhoff (matric flux) potential Phi(psi) = integral_0^psi plc ds [MPa]. ----------!
   !      Closed form for kexp in {1,2}; a fixed quadrature otherwise. Phi(0)=0, Phi(psi<0)<0,  !
   !      strictly increasing in psi.                                                           !
   pure real(wp) function flux_potential(psi, psi50, kexp) result(phi)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r = max(psi/psi50, 0.0_wp)                 ! r >= 0
      if (abs(kexp - 1.0_wp) < 1.0e-9_wp) then
         phi = psi50 * log(1.0_wp + r)
      else if (abs(kexp - 2.0_wp) < 1.0e-9_wp) then
         phi = psi50 * atan(r)
      else
         phi = psi50 * kirchhoff_integral(r, kexp)
      end if
   end function flux_potential

   !----- integral_0^r du/(1+u^kexp) by 7-point Gauss-Legendre quadrature: exact for polynomials !
   !      up to degree 13, and far below any model tolerance for this smooth integrand. The sum is !
   !      written out rather than handed a function, because ifx allocates a lock-guarded record   !
   !      on every call that passes a contained function as an argument, which made more than four !
   !      threads slower than four (#325). u=0 is guarded as in plc_retained. -------------------!
   pure real(wp) function kirchhoff_integral(r, kexp) result(g)
      real(wp), intent(in) :: r, kexp
      real(wp), parameter :: node(7)   = [ -0.9491079123427585_wp, -0.7415311855993945_wp,     &
                                           -0.4058451513773972_wp,  0.0000000000000000_wp,     &
                                            0.4058451513773972_wp,  0.7415311855993945_wp,     &
                                            0.9491079123427585_wp ]
      real(wp), parameter :: weight(7) = [  0.1294849661688697_wp,  0.2797053914892766_wp,     &
                                            0.3818300505051189_wp,  0.4179591836734694_wp,     &
                                            0.3818300505051189_wp,  0.2797053914892766_wp,     &
                                            0.1294849661688697_wp ]
      real(wp)    :: mid, half, u, acc
      integer(ik) :: k
      mid  = 0.5_wp * r                          ! the nodes on [-1,1], mapped onto [0,r]
      half = 0.5_wp * r
      acc  = 0.0_wp
      do k = 1_ik, 7_ik
         u = mid + half * node(k)
         if (u <= 0.0_wp) then
            acc = acc + weight(k)
         else
            acc = acc + weight(k) / (1.0_wp + u**kexp)
         end if
      end do
      g = half * acc
   end function kirchhoff_integral

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

   !----- ED2 cumulative-exponential root fraction in a soil layer spanning depths [z_top, z_bot]     !
   !      (both >= 0, below surface, z_bot > z_top): frac = beta^(z_top/D) - beta^(z_bot/D), with      !
   !      D = root_depth and beta in (0,1). Shallower layers get more roots; depths clamp to [0, D] so  !
   !      layers below the rooting depth contribute 0. Summed over [0, D] the profile telescopes to     !
   !      1 - beta (ED2 convention; only the RELATIVE distribution enters the conductance weights).     !
   elemental real(wp) function root_fraction_profile(root_beta, root_depth, z_top, z_bot) result(frac)
      real(wp), intent(in) :: root_beta   !< [-]  root-profile decay (0,1); smaller => shallower
      real(wp), intent(in) :: root_depth  !< [m]  maximum rooting depth (> 0)
      real(wp), intent(in) :: z_top       !< [m]  depth of the layer top    (>= 0)
      real(wp), intent(in) :: z_bot       !< [m]  depth of the layer bottom (> z_top)
      real(wp) :: inv_d, a, b
      inv_d = 1.0_wp / max(root_depth, tiny_num)
      a = min(max(z_top, 0.0_wp), root_depth) * inv_d
      b = min(max(z_bot, 0.0_wp), root_depth) * inv_d
      frac = root_beta**a - root_beta**b
   end function root_fraction_profile

end module meds_hydr_lib
