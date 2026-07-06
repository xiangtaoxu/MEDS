!==========================================================================================!
! meds_soil_thermal -- soil thermal properties (conductivity + volumetric heat capacity) and  !
! the per-column thermal-parameter assembly for the soil-heat kernel (design                   !
! MEDS_ENERGY_BALANCE_DESIGN.md, 4c/8). The thermal twin of meds_soil_parameters.               !
!                                                                                          !
!   * soil_thermal_cond -- Johansen/Farouki interpolation between dry and saturated with an     !
!     ice-aware saturated geometric mean and a Kersten number clamped to [0,1] (FPE-safe at the  !
!     dry end where log10 -> -inf).                                                              !
!   * soil_heat_cap_vol -- dry matrix + phase-appropriate soil-water heat capacity.              !
! Both `elemental`; the Bolton-style empirical coefficients (0.05 saturation floor) are          !
! numerical regularization, not tunable model parameters.                                        !
!==========================================================================================!
module meds_soil_thermal
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, cp_liq, cp_ice, k_water, k_ice, tiny_num
   use meds_biophysics_types, only : soil_thermal_params_t, n_soil_layer_max
   implicit none
   private

   public :: soil_thermal_cond, soil_heat_cap_vol, build_soil_thermal

   real(wp), parameter :: SR_FLOOR = 0.05_wp     !< Kersten log10 floor (Johansen validity limit)

contains

   !---------------------------------------------------------------------------------------!
   ! Soil thermal conductivity [W/m/K] (Johansen). fliq = liquid fraction of the soil water. !
   !---------------------------------------------------------------------------------------!
   elemental function soil_thermal_cond(theta, fliq, theta_sat, k_solid, k_dry) result(kappa)
      real(wp), intent(in) :: theta, fliq, theta_sat, k_solid, k_dry
      real(wp)             :: kappa, s_r, k_e, kappa_sat
      s_r = min(max(theta / max(theta_sat, tiny_num), SR_FLOOR), 1.0_wp)      ! relative saturation
      k_e = min(max(log10(s_r) + 1.0_wp, 0.0_wp), 1.0_wp)                     ! Kersten number, clamped
      !----- Saturated geometric mean, ice-aware (fixes ED2's liquid-only kappa_sat). -----!
      kappa_sat = k_solid ** (1.0_wp - theta_sat)                                            &
                  * k_water ** (theta_sat * fliq) * k_ice ** (theta_sat * (1.0_wp - fliq))
      kappa = k_e * kappa_sat + (1.0_wp - k_e) * k_dry
   end function soil_thermal_cond

   !---------------------------------------------------------------------------------------!
   ! Effective volumetric heat capacity [J/m3/K] = dry matrix + phase-appropriate water.    !
   !---------------------------------------------------------------------------------------!
   elemental function soil_heat_cap_vol(theta, fliq, dry_cvol) result(c_eff)
      real(wp), intent(in) :: theta, fliq, dry_cvol
      real(wp)             :: c_eff
      c_eff = dry_cvol + theta * rho_h2o * (fliq * cp_liq + (1.0_wp - fliq) * cp_ice)
   end function soil_heat_cap_vol

   !---------------------------------------------------------------------------------------!
   ! Assemble a soil_thermal_params_t from plain (uniform) inputs -- the standalone builder  !
   ! that keeps the kernel testable now (per-layer texture + TOML land at P3, with hydrology).!
   !---------------------------------------------------------------------------------------!
   subroutine build_soil_thermal(n_active, k_solid, k_dry, dry_cvol, therm)
      integer(ik),                intent(in)  :: n_active
      real(wp),                   intent(in)  :: k_solid, k_dry, dry_cvol
      type(soil_thermal_params_t), intent(out) :: therm
      therm%nzg_active = n_active
      therm%soil_solid_conductivity(1:n_active) = k_solid
      therm%soil_dry_conductivity(1:n_active)   = k_dry
      therm%soil_dry_heat_capacity(1:n_active)  = dry_cvol
   end subroutine build_soil_thermal

end module meds_soil_thermal
