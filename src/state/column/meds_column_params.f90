!==========================================================================================!
! meds_column_params -- the static per-column PARAMETER bundles (geometry, texture, thermal    !
! properties) and their `pure` assemblers. These describe the same stores meds_column_reservoirs !
! holds, but they are derived ONCE per column and never integrated, so they are a separate       !
! module from the reservoirs (the plan's placement rule 3: parameters are not state).            !
!                                                                                          !
! The builders are state-free constructors -- they read plain scalar texture/geometry inputs   !
! and the shared retention curve, never soil state -- so this module stays below the kernels.  !
!==========================================================================================!
module meds_column_params
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num
   use meds_column_constants, only : n_soil_layer_max
   use meds_hydr_lib,         only : soil_theta_from_psi, SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL
   implicit none
   private

   public :: soil_params_t, soil_thermal_params_t
   public :: build_soil_hydr_params, build_soil_therm_params
   public :: curve_a, curve_n

   real(wp),    parameter :: PSI_WP = -152.96_wp           !< [m] -1.5 MPa head (wilting-point derivation)

   !----- Per-column geometry + texture (assembled once per site; ED2 negative-z convention:   !
   !      interface elevations z <= 0 below ground; dz, dz_node are positive magnitudes). Fixed- !
   !      size (n_soil_layer_max) so the hydrology kernel stays allocatable-free and GPU-eligible.!
   type :: soil_params_t
      integer(ik) :: n_active  = n_soil_layer_max            !< active layer count (<= n_soil_layer_max)
      integer(ik) :: retention = SOIL_RETENTION_VG           !< curve family
      real(wp) :: soil_layer_z(n_soil_layer_max+1) = 0.0_wp  !< [m] interface elevations (<= 0, ED2 slz)
      real(wp) :: z_node(n_soil_layer_max)  = 0.0_wp         !< [m] node (mid) elevations (<= 0)
      real(wp) :: dz(n_soil_layer_max)      = 0.0_wp         !< [m] layer thickness (> 0)
      real(wp) :: dz_node(n_soil_layer_max) = 0.0_wp         !< [m] internode spacing (> 0)
      real(wp) :: theta_sat(n_soil_layer_max) = 0.0_wp       !< [m3/m3] porosity
      real(wp) :: theta_res(n_soil_layer_max) = 0.0_wp       !< [m3/m3] residual water content
      real(wp) :: ksat(n_soil_layer_max)      = 0.0_wp       !< [m/s] saturated conductivity
      real(wp) :: vg_alpha(n_soil_layer_max)  = 0.0_wp       !< [1/m] van Genuchten inverse air-entry
      real(wp) :: vg_n(n_soil_layer_max)      = 0.0_wp       !< [-] van Genuchten pore-size index (> 1)
      real(wp) :: psi_sat(n_soil_layer_max)   = 0.0_wp       !< [m] Campbell air-entry potential (option)
      real(wp) :: b_camp(n_soil_layer_max)    = 0.0_wp       !< [-] Campbell exponent (option)
      real(wp) :: theta_fc(n_soil_layer_max)  = 0.0_wp       !< [m3/m3] field capacity (DERIVED)
      real(wp) :: theta_wp(n_soil_layer_max)  = 0.0_wp       !< [m3/m3] wilting point (DERIVED)
      real(wp) :: root_frac(n_soil_layer_max) = 0.0_wp       !< [-] normalized root fraction (sum = 1)
   end type soil_params_t

   !----- Per-column soil THERMAL texture (geometry + porosity arrive via soil_params_t). ----!
   type :: soil_thermal_params_t
      integer(ik) :: nzg_active = n_soil_layer_max
      real(wp) :: soil_solid_conductivity(n_soil_layer_max) = 0.0_wp   !< [W/m/K] kappa_solid
      real(wp) :: soil_dry_conductivity(n_soil_layer_max)   = 0.0_wp   !< [W/m/K] kappa_dry
      real(wp) :: soil_dry_heat_capacity(n_soil_layer_max)  = 0.0_wp   !< [J/m3/K] dry-matrix vol. heat cap
   end type soil_thermal_params_t

contains

   !=======================================================================================!
   !  Per-column soil PARAMETER assemblers. Both are `pure` -- they take plain scalar texture/    !
   !  geometry inputs and fill a parameter struct; they depend on NO soil STATE (theta/energy),  !
   !  so they are state-free constructors, not part of the fast loop. Per-layer texture + TOML    !
   !  wiring land at P3; these standalone builders keep the kernels testable now.                  !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Assemble a per-column soil_params_t: the ED2 negative-z geometry (exponential generator !
   ! OR an explicit soil_layer_z interface array), uniform texture broadcast, an exponential  !
   ! root profile, and the DERIVED thresholds theta_fc/theta_wp (from the retention curve).   !
   !---------------------------------------------------------------------------------------!
   pure subroutine build_soil_hydr_params(n_active, retention, soil_depth, grid_growth, theta_sat,  &
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
         if (abs(grid_growth) < tiny_num) then          ! uniform grid: the exact 0/0 limit of below
            do k = 1_ik, n_active + 1_ik
               params%soil_layer_z(k) = -soil_depth * real(k - 1_ik, wp) / real(n_active, wp)
            end do
         else
            denom = exp(grid_growth) - 1.0_wp
            do k = 1_ik, n_active + 1_ik
               params%soil_layer_z(k) = -soil_depth                                         &
                  * (exp(grid_growth * real(k - 1_ik, wp) / real(n_active, wp)) - 1.0_wp) / denom
            end do
         end if
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

      !----- Derived thresholds (from the retention curve). -------------------------------!
      params%theta_fc(1:n_active) =                                                          &
         soil_theta_from_psi(retention, psi_fc_m, theta_sat, theta_res, par_a, par_n)
      params%theta_wp(1:n_active) =                                                          &
         soil_theta_from_psi(retention, PSI_WP, theta_sat, theta_res, par_a, par_n)

      !----- Exponential root profile (z_node <= 0 => decays with depth), normalized. -----!
      rsum = 0.0_wp
      do k = 1_ik, n_active
         params%root_frac(k) = exp(root_beta * params%z_node(k)) * params%dz(k)
         rsum = rsum + params%root_frac(k)
      end do
      if (rsum > tiny_num) params%root_frac(1:n_active) = params%root_frac(1:n_active) / rsum
   end subroutine build_soil_hydr_params

   !---------------------------------------------------------------------------------------!
   ! Assemble a soil_thermal_params_t from plain (uniform) inputs -- the thermal twin of      !
   ! build_soil_hydr_params (per-layer texture + TOML land at P3, with hydrology).            !
   !---------------------------------------------------------------------------------------!
   pure subroutine build_soil_therm_params(n_active, k_solid, k_dry, dry_cvol, therm)
      integer(ik),                 intent(in)  :: n_active
      real(wp),                    intent(in)  :: k_solid, k_dry, dry_cvol
      type(soil_thermal_params_t), intent(out) :: therm
      therm%nzg_active = n_active
      therm%soil_solid_conductivity(1:n_active) = k_solid
      therm%soil_dry_conductivity(1:n_active)   = k_dry
      therm%soil_dry_heat_capacity(1:n_active)  = dry_cvol
   end subroutine build_soil_therm_params

   !=======================================================================================!
   !  Retention-curve parameter accessors: (alpha, n) for van Genuchten, (psi_sat, b) for       !
   !  Campbell. They live HERE, beside soil_params_t, because they are a property of that type   !
   !  -- which family its two generic curve parameters mean. They were previously private        !
   !  duplicates inside meds_soil_water; every consumer (the Richards solver, the diagnostic      !
   !  psi read-off) now shares this ONE mapping, so a family added later cannot be taught to       !
   !  one caller and not the other.                                                                !
   !=======================================================================================!
   pure function curve_a(params, k) result(a)
      type(soil_params_t), intent(in) :: params
      integer(ik),         intent(in) :: k
      real(wp)                        :: a
      if (params%retention == SOIL_RETENTION_CAMPBELL) then
         a = params%psi_sat(k)
      else
         a = params%vg_alpha(k)
      end if
   end function curve_a

   pure function curve_n(params, k) result(nn)
      type(soil_params_t), intent(in) :: params
      integer(ik),         intent(in) :: k
      real(wp)                        :: nn
      if (params%retention == SOIL_RETENTION_CAMPBELL) then
         nn = params%b_camp(k)
      else
         nn = params%vg_n(k)
      end if
   end function curve_n

end module meds_column_params
