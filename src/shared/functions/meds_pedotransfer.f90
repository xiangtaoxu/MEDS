! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_pedotransfer -- a soil layer's hydraulic and thermal properties from what soil maps !
! report: texture, organic carbon, CEC, pH, bulk density and coarse fragments. The fine    !
! earth's van Genuchten curve is Hodnett & Tomasella (2002) and its saturated conductivity !
! Ottoni et al. (2019); coarse fragments hold no water and slow the flow (Saxton & Rawls   !
! 2006). The thermal properties are CLM5's mineral-soil formulas. The equations, and the   !
! tests behind the choice, are docs/science/soil_biophysics.md §Hydraulic parameters from  !
! soil texture and §The column from the soil maps, eqs. (6)-(9).                           !
!==========================================================================================!
module meds_pedotransfer
   use meds_kinds,           only : wp
   use meds_constants,       only : grav_head, day_sec, tiny_num
   use meds_water_retention, only : SOIL_RETENTION_VG, soil_theta_from_psi
   implicit none
   private

   public :: soil_texture_t, soil_hydraulics_from_texture, soil_thermal_from_texture

   !----- One soil layer as the soil maps describe it (SoilGrids 2.0's quantities). -------!
   type :: soil_texture_t
      real(wp) :: sand         = 0.0_wp   !< [%]           sand share of the fine earth
      real(wp) :: silt         = 0.0_wp   !< [%]           silt share of the fine earth
      real(wp) :: clay         = 0.0_wp   !< [%]           clay share of the fine earth
      real(wp) :: org_carbon   = 0.0_wp   !< [g/kg]        organic carbon of the fine earth
      real(wp) :: cec          = 0.0_wp   !< [cmol(c)/kg]  cation exchange capacity
      real(wp) :: ph           = 0.0_wp   !< [-]           pH in water
      real(wp) :: bulk_density = 0.0_wp   !< [g/cm3]       bulk density of the fine earth
      real(wp) :: coarse_frag  = 0.0_wp   !< [% volume]    coarse fragments (> 2 mm) in the layer
   end type soil_texture_t

   real(wp), parameter :: kpa_per_mpa      = 1000.0_wp                !< [kPa/MPa]

   !----- Ottoni et al. (2019): Ks = a (theta_sat - theta(-33 kPa))^b, Ks in m/day. -------!
   real(wp), parameter :: psi_33kpa        = -0.033_wp / grav_head   !< [m] -33 kPa as water head
   real(wp), parameter :: ottoni_ks_coef   = 19.31_wp                !< [m/day]
   real(wp), parameter :: ottoni_ks_exp    = 1.948_wp                !< [-]
   !----- Density of the mineral grains and of the fragments (Saxton & Rawls 2006). -------!
   real(wp), parameter :: particle_density = 2.65_wp                 !< [g/cm3]

   !----- CLM5 mineral soil (Lawrence & Slater 2008; CLM5 technote 2.7.2): the solids'     !
   !      conductivity and heat capacity per sand and clay, and the dry soil's             !
   !      conductivity from its bulk density (Johansen 1975); grains are 2700 kg/m3. ------!
   real(wp), parameter :: k_quartz         = 8.80_wp     !< [W/m/K]    sand grains
   real(wp), parameter :: k_clay_mineral   = 2.92_wp     !< [W/m/K]    clay grains
   real(wp), parameter :: cv_sand          = 2.128e6_wp  !< [J/m3/K]   sand grains
   real(wp), parameter :: cv_clay          = 2.385e6_wp  !< [J/m3/K]   clay grains
   real(wp), parameter :: clm_grain_density = 2700.0_wp  !< [kg/m3]

contains

   !----- The van Genuchten curve and Ks of a soil layer (docs/science/soil_biophysics.md  !
   !      eqs. (6)-(8)). Alpha comes out per kPa and is turned into per metre of water     !
   !      head, the unit of the model's soil potential. -----------------------------------!
   elemental subroutine soil_hydraulics_from_texture(texture, theta_sat, theta_res, vg_alpha, vg_n, ksat)
      type(soil_texture_t), intent(in)  :: texture
      real(wp),             intent(out) :: theta_sat   !< [m3/m3] water content at saturation, per layer volume
      real(wp),             intent(out) :: theta_res   !< [m3/m3] residual water content, per layer volume
      real(wp),             intent(out) :: vg_alpha    !< [1/m]   van Genuchten alpha
      real(wp),             intent(out) :: vg_n        !< [-]     van Genuchten n
      real(wp),             intent(out) :: ksat        !< [m/s]   saturated hydraulic conductivity of the layer
      real(wp) :: oc_percent, theta_33, drainable, rock_volume, rock_mass, matrix_density_ratio

      associate(sand => texture%sand, silt => texture%silt, clay => texture%clay, cec => texture%cec,  &
                ph => texture%ph, bulk_density => texture%bulk_density)
         oc_percent = texture%org_carbon / 10.0_wp

         !----- The fine earth's curve (Hodnett & Tomasella 2002, Table 8). ---------------!
         vg_alpha  = exp((-2.294_wp - 3.526_wp*silt + 2.440_wp*oc_percent - 0.076_wp*cec - 11.331_wp*ph  &
                          + 0.019_wp*silt**2) / 100.0_wp) * grav_head * kpa_per_mpa
         vg_n      = exp((62.986_wp - 0.833_wp*clay - 0.529_wp*oc_percent + 0.593_wp*ph                 &
                          + 0.007_wp*clay**2 - 0.014_wp*sand*silt) / 100.0_wp)
         theta_sat = (81.799_wp + 0.099_wp*clay - 31.42_wp*bulk_density + 0.018_wp*cec + 0.451_wp*ph     &
                      - 0.0005_wp*sand*clay) / 100.0_wp
         !----- The regression dips below zero on sandy, alkaline, low-CEC soils, which    !
         !      hold no residual water. ---------------------------------------------------!
         theta_res = max((22.733_wp - 0.164_wp*sand + 0.235_wp*cec - 0.831_wp*ph + 0.0018_wp*clay**2  &
                          + 0.0026_wp*sand*clay) / 100.0_wp, 0.0_wp)

         !----- Ks of the fine earth from the water that drains between saturation, at     !
         !      psi = 0 on this curve, and -33 kPa (Ottoni et al. 2019). ------------------!
         theta_33  = soil_theta_from_psi(SOIL_RETENTION_VG, psi_33kpa, theta_sat, theta_res, vg_alpha, vg_n)
         drainable = max(theta_sat - theta_33, tiny_num)
         ksat      = ottoni_ks_coef / day_sec * drainable**ottoni_ks_exp

         !----- Coarse fragments hold no water and slow the flow: Saxton & Rawls (2006)    !
         !      eq. 22, on the fragments' share of the mass, which eq. 19 gives from their !
         !      share of the volume. ------------------------------------------------------!
         rock_volume          = texture%coarse_frag / 100.0_wp
         matrix_density_ratio = bulk_density / particle_density
         rock_mass = rock_volume / max(matrix_density_ratio*(1.0_wp - rock_volume) + rock_volume, tiny_num)
         theta_sat = theta_sat * (1.0_wp - rock_volume)
         theta_res = theta_res * (1.0_wp - rock_volume)
         ksat      = ksat * (1.0_wp - rock_mass) / (1.0_wp - rock_mass*(1.0_wp - 1.5_wp*matrix_density_ratio))
      end associate
   end subroutine soil_hydraulics_from_texture

   !----- The thermal properties of one soil layer by CLM5's mineral-soil formulas         !
   !      (docs/science/soil_biophysics.md eq. (9)). The heat capacity holds the solids    !
   !      only, since the model counts the water and ice itself. --------------------------!
   elemental subroutine soil_thermal_from_texture(texture, theta_sat, k_solid, k_dry, dry_cvol)
      type(soil_texture_t), intent(in)  :: texture
      real(wp),             intent(in)  :: theta_sat  !< [m3/m3] porosity of the layer
      real(wp),             intent(out) :: k_solid    !< [W/m/K]  conductivity of the soil solids
      real(wp),             intent(out) :: k_dry      !< [W/m/K]  conductivity of the dry soil
      real(wp),             intent(out) :: dry_cvol   !< [J/m3/K] heat capacity of the solids per layer volume
      real(wp) :: bulk_density

      associate(sand => texture%sand, clay => texture%clay)
         k_solid      = (k_quartz*sand + k_clay_mineral*clay) / (sand + clay)
         bulk_density = clm_grain_density * (1.0_wp - theta_sat)
         k_dry        = (0.135_wp*bulk_density + 64.7_wp) / (clm_grain_density - 0.947_wp*bulk_density)
         dry_cvol     = (cv_sand*sand + cv_clay*clay) / (sand + clay) * (1.0_wp - theta_sat)
      end associate
   end subroutine soil_thermal_from_texture

end module meds_pedotransfer
