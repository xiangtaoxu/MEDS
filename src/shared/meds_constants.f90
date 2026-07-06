!==========================================================================================!
! meds_constants -- compile-time constants and a few pure scalar helpers.                  !
!                                                                                          !
! No mutable state lives here. Calendar constants follow ED2 (proleptic Gregorian mean     !
! year). Tolerances default to the values used by the fusion/fission conservation checks.   !
!==========================================================================================!
module meds_constants
   use meds_kinds, only : wp
   implicit none
   public

   !----- Geometry. -----------------------------------------------------------------------!
   real(wp), parameter :: pi   = 3.14159265358979323846_wp
   real(wp), parameter :: pio4 = 0.25_wp * pi                 !< basal area = pio4 * dbh^2

   !----- Calendar (mean Gregorian year). -------------------------------------------------!
   real(wp), parameter :: yr_day      = 365.2425_wp           !< days per year
   real(wp), parameter :: day_sec     = 86400.0_wp            !< seconds per day
   real(wp), parameter :: yr_sec      = yr_day * day_sec      !< seconds per year
   real(wp), parameter :: mon_per_yr  = 12.0_wp
   real(wp), parameter :: day_per_mon = yr_day / mon_per_yr   !< ~30.4 days per month

   !----- Physical constants (leaf biophysics: photosynthesis / stomatal conductance). ----!
   real(wp), parameter :: r_gas       = 8.314462618_wp       !< [J/mol/K] universal gas constant
   real(wp), parameter :: t_kelvin    = 273.15_wp            !< [K] 0 degC in Kelvin
   real(wp), parameter :: p_std       = 101325.0_wp          !< [Pa] standard atmospheric pressure
   real(wp), parameter :: t_ref_photo = 298.15_wp            !< [K] 25 degC reference for photosynthesis
   real(wp), parameter :: stefan      = 5.670374419e-8_wp    !< [W/m2/K4] Stefan-Boltzmann (canopy RT longwave)
   real(wp), parameter :: halfpi      = 0.5_wp * pi          !< [rad] pi/2 (leaf-inclination domain upper bound)
   real(wp), parameter :: grav_head   = 9.804e-3_wp          !< [MPa/m] hydrostatic head (rho_w*g); plant hydraulics

   !----- Canopy-air CO2 balance + soil respiration (carbon <-> mole conversions). --------!
   real(wp), parameter :: mmdry          = 0.0289655_wp        !< [kg/mol] dry-air molar mass (~28.97 g/mol)
   real(wp), parameter :: umol_2_kgC     = 1.20107e-8_wp       !< [kgC/umol] carbon mass per umol CO2 (12.0107 g/mol)
   real(wp), parameter :: kgC_2_umol     = 1.0_wp / umol_2_kgC !< [umol/kgC]  (~8.3259e7)
   real(wp), parameter :: kgCday_2_umols = kgC_2_umol / day_sec  !< [(umol/m2/s) per (kgC/m2/day)]  (~963.6)
   real(wp), parameter :: umols_per_kgCyr = kgC_2_umol / yr_sec  !< [(umol/m2/s) per (kgC/m2/yr)]   (~2.638)
   !----- DAMM heterotrophic respiration (Davidson et al. 2012 GCB 18:371). ---------------!
   real(wp), parameter :: r_gas_kj        = r_gas * 1.0e-3_wp  !< [kJ/mol/K] gas const in kJ (DAMM Arrhenius; ~8.3145e-3)
   real(wp), parameter :: o2_air_frac     = 0.209_wp           !< [-] volume fraction of O2 in air (DAMM Eqn 6)
   real(wp), parameter :: damm_flux_factor = 1.0e4_wp / 3600.0_wp * 1000.0_wp / 12.011_wp
                                     !< [(umol/m2/s) per (mgC cm-3 h-1 . cm)] ~ 231.27: cm2/m2, /h->/s, mgC->umol

   !----- Soil hydrology (vertical water column). -----------------------------------------!
   real(wp), parameter :: rho_h2o     = 1000.0_wp            !< [kg/m3] liquid-water density (kg/m2 <-> m <-> m3/m3)
   real(wp), parameter :: grav        = 9.80665_wp           !< [m/s2] gravitational acceleration (head + alpha_soil)
   real(wp), parameter :: r_wv        = 461.5_wp             !< [J/kg/K] water-vapour gas constant (Philip alpha_soil)
   real(wp), parameter :: latent_heat_vap = 2.501e6_wp       !< [J/kg] latent heat of vaporization

   !----- Energy balance (thermal). specific heats [J/kg/K], conductivities [W/m/K]. ------!
   real(wp), parameter :: cp_air = 1004.6_wp                 !< dry-air cp (sensible fluxes, CAS enthalpy)
   real(wp), parameter :: cp_vap = 1859.0_wp                 !< water-VAPOUR cp (CAS + vapour-enthalpy twins)
   real(wp), parameter :: cp_liq = 4186.0_wp                 !< liquid-water cp (store water, advected enthalpy)
   real(wp), parameter :: cp_ice = 2093.0_wp                 !< ice cp (frozen-store heat capacity)
   real(wp), parameter :: latent_heat_fusion = 3.34e5_wp     !< [J/kg] latent heat of fusion (freeze/thaw plateau)
   real(wp), parameter :: k_water = 0.57_wp                  !< liquid-water thermal conductivity
   real(wp), parameter :: k_ice   = 2.29_wp                  !< ice thermal conductivity
   real(wp), parameter :: k_air   = 0.025_wp                 !< air thermal conductivity
   real(wp), parameter :: t_3ple  = 273.16_wp                !< [K] triple point (distinct from t_kelvin = 273.15)
   !----- Internal-energy zero references (ED2 cmtl2uext-consistent; see meds_thermo, design 3.1). --!
   real(wp), parameter :: tsupercool_liq = t_3ple - (cp_ice * t_3ple + latent_heat_fusion) / cp_liq
   real(wp), parameter :: tsupercool_vap = t_3ple                                                     &
                          - (cp_ice * t_3ple + latent_heat_vap + latent_heat_fusion) / cp_vap

   !----- Numerical safety. ---------------------------------------------------------------!
   real(wp), parameter :: tiny_num   = 1.0e-30_wp
   real(wp), parameter :: almost_one = 1.0_wp - 1.0e-7_wp
   real(wp), parameter :: lnexp_min  = -38.0_wp               !< clamp for exp() arguments
   real(wp), parameter :: lnexp_max  =  38.0_wp

   !----- Default conservation tolerances. ------------------------------------------------!
   real(wp), parameter :: size_tol = 0.01_wp                  !< 1% basal-area conservation
   real(wp), parameter :: area_tol = 1.0e-5_wp                !< patch-area sum tolerance

contains

   !---------------------------------------------------------------------------------------!
   ! Numerically safe exponential: clamps the argument to avoid under/overflow.            !
   !---------------------------------------------------------------------------------------!
   elemental pure function safe_exp(x) result(y)
      real(wp), intent(in) :: x
      real(wp)             :: y
      y = exp(min(max(x, lnexp_min), lnexp_max))
   end function safe_exp

end module meds_constants
