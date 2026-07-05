!==========================================================================================!
! meds_soil_parameters -- soil hydraulic constitutive curves and the per-column parameter    !
! assembly for the vertical-hydrology kernel (design MEDS_VERTICAL_HYDROLOGY_DESIGN.md, 3b/7). !
!                                                                                          !
! Two retention closures, selected per run by SOIL_RETENTION_VG (default) / _CAMPBELL:        !
!   * van Genuchten-Mualem (van Genuchten 1980; Mualem 1976) -- the default, the direction     !
!     the ED2 community is moving; smooth through air-entry, better fit to data.                !
!   * Campbell / Clapp-Hornberger (Campbell 1974; Clapp & Hornberger 1978) -- the ED2/CLM        !
!     option, for BC64 reproducibility.                                                          !
! Both have closed-form theta(psi) and psi(theta) inverses, so all four curve functions are      !
! `elemental`, branch-light and FPE-safe under -fpe0 / -Ktrap=fp (Se clamped, K floored, C        !
! guarded). Coupling to the plant hydraulics kernel is through the POTENTIAL psi (curve-           !
! independent at the interface), so the curve is a free choice.                                    !
!==========================================================================================!
module meds_soil_parameters
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num
   use meds_biophysics_types, only : soil_params_t, n_soil_layer_max,                         &
                                     SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   implicit none
   private

   public :: soil_theta_of_psi, soil_psi_of_theta, soil_hydr_cond, soil_moist_cap
   public :: build_soil_params

   !----- Numerical floors (all a property of the curve regularization, not of uptake). ----!
   real(wp), parameter :: K_MIN  = 1.16e-13_wp    !< [m/s] ED2 hydcond_min (no zero-conductance lock)
   real(wp), parameter :: C_MIN  = 1.0e-9_wp      !< [1/m] specific-capacity floor (vG C -> 0 at saturation)
   real(wp), parameter :: SE_MIN = 1.0e-9_wp      !< effective-saturation floor (residual/air-dry)
   real(wp), parameter :: PSI_WP = -152.96_wp     !< [m] -1.5 MPa in head (wilting-point definition)

contains

   !---------------------------------------------------------------------------------------!
   ! theta(psi): volumetric water content from matric potential [m, <= 0].                  !
   !---------------------------------------------------------------------------------------!
   elemental function soil_theta_of_psi(retention, psi, theta_sat, theta_res, par_a, par_n)  &
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
   end function soil_theta_of_psi

   !---------------------------------------------------------------------------------------!
   ! psi(theta): matric potential [m, <= 0] from water content -- closed-form inverse.      !
   !---------------------------------------------------------------------------------------!
   elemental function soil_psi_of_theta(retention, theta, theta_sat, theta_res, par_a, par_n) &
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
   end function soil_psi_of_theta

   !---------------------------------------------------------------------------------------!
   ! K(theta): unsaturated hydraulic conductivity [m/s], floored at K_MIN.                  !
   !---------------------------------------------------------------------------------------!
   elemental function soil_hydr_cond(retention, theta, theta_sat, theta_res, par_a, par_n,    &
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
   end function soil_hydr_cond

   !---------------------------------------------------------------------------------------!
   ! C(psi) = dtheta/dpsi: specific moisture capacity [1/m], >= 0, floored at C_MIN.        !
   !---------------------------------------------------------------------------------------!
   elemental function soil_moist_cap(retention, psi, theta_sat, theta_res, par_a, par_n)      &
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
   end function soil_moist_cap

   !---------------------------------------------------------------------------------------!
   ! Assemble a per-column soil_params_t from plain inputs: the ED2 negative-z geometry      !
   ! (exponential generator OR an explicit soil_layer_z interface array), uniform texture     !
   ! broadcast to every layer, an exponential root profile, and the DERIVED thresholds        !
   ! theta_fc/theta_wp. (Config wiring -- per-layer texture, TOML -- is deferred to P3; this   !
   ! standalone builder keeps the kernel testable now.)                                        !
   !---------------------------------------------------------------------------------------!
   subroutine build_soil_params(n_active, retention, soil_depth, grid_growth, theta_sat,      &
                                theta_res, ksat, par_a, par_n, root_beta, psi_fc_m, params,    &
                                soil_layer_z_in)
      integer(ik),         intent(in)  :: n_active, retention
      real(wp),            intent(in)  :: soil_depth, grid_growth, theta_sat, theta_res
      real(wp),            intent(in)  :: ksat, par_a, par_n, root_beta, psi_fc_m
      type(soil_params_t), intent(out) :: params
      real(wp), optional,  intent(in)  :: soil_layer_z_in(:)
      integer(ik) :: k
      real(wp)    :: denom, rsum

      params%n_active  = n_active
      params%retention = retention

      !----- Interface elevations soil_layer_z(k) <= 0 (k=1 surface, deeper = more negative). -!
      if (present(soil_layer_z_in)) then
         params%soil_layer_z(1:n_active+1) = soil_layer_z_in(1:n_active+1)
      else
         denom = exp(grid_growth) - 1.0_wp
         do k = 1_ik, n_active + 1_ik
            params%soil_layer_z(k) = -soil_depth                                            &
               * (exp(grid_growth * real(k - 1_ik, wp) / real(n_active, wp)) - 1.0_wp) / denom
         end do
      end if

      !----- Thicknesses, node elevations, internode spacings (dz, dz_node > 0). ----------!
      do k = 1_ik, n_active
         params%dz(k)     = params%soil_layer_z(k) - params%soil_layer_z(k+1)
         params%z_node(k) = 0.5_wp * (params%soil_layer_z(k) + params%soil_layer_z(k+1))
      end do
      do k = 1_ik, n_active - 1_ik
         params%dz_node(k) = params%z_node(k) - params%z_node(k+1)
      end do
      params%dz_node(n_active) = params%dz(n_active)      ! unused at the bottom face; kept > 0

      !----- Uniform texture broadcast to every active layer. -----------------------------!
      params%theta_sat(1:n_active) = theta_sat
      params%theta_res(1:n_active) = theta_res
      params%ksat(1:n_active)      = ksat
      if (retention == SOIL_RETENTION_CAMPBELL) then
         params%psi_sat(1:n_active) = par_a
         params%b_camp(1:n_active)  = par_n
      else
         params%vg_alpha(1:n_active) = par_a
         params%vg_n(1:n_active)     = par_n
      end if

      !----- Derived thresholds (from the curve, per design 7). ---------------------------!
      params%theta_fc(1:n_active) =                                                          &
         soil_theta_of_psi(retention, psi_fc_m, theta_sat, theta_res, par_a, par_n)
      params%theta_wp(1:n_active) =                                                          &
         soil_theta_of_psi(retention, PSI_WP, theta_sat, theta_res, par_a, par_n)

      !----- Exponential root profile (z_node <= 0 => decays with depth), normalized. -----!
      rsum = 0.0_wp
      do k = 1_ik, n_active
         params%root_frac(k) = exp(root_beta * params%z_node(k)) * params%dz(k)
         rsum = rsum + params%root_frac(k)
      end do
      if (rsum > tiny_num) params%root_frac(1:n_active) = params%root_frac(1:n_active) / rsum
   end subroutine build_soil_params

end module meds_soil_parameters
