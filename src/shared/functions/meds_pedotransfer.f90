! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_pedotransfer -- a soil layer's water-retention curve and saturated conductivity     !
! from what soil maps report: texture, organic carbon, CEC, pH, bulk density and coarse    !
! fragments. One pedotransfer function, closed-form and elemental:                         !
!                                                                                          !
!   * the retention curve of the fine earth (< 2 mm) is Hodnett & Tomasella (2002,         !
!     Geoderma 108:155) Table 8: van Genuchten theta_sat, theta_res, alpha and n, fitted   !
!     to tropical soils -- ED2's soil_hydro_scheme 2;                                      !
!   * the saturated conductivity of the fine earth is Ottoni et al. (2019, J. Hydrol.      !
!     575:1345): a power of the water that drains between saturation and -33 kPa, as ED2   !
!     uses it with that curve;                                                             !
!   * coarse fragments (> 2 mm) hold no water, so the layer's curve is the fine earth's    !
!     scaled to the fine earth's share of the volume; they slow the flow by Saxton & Rawls !
!     (2006, SSSA J. 70:1569) eq. 22.                                                      !
!                                                                                          !
! The result is the van Genuchten-Mualem curve meds_water_retention evaluates. The science !
! page (docs/science/soil_biophysics.md) has the equations and the tests behind the choice.!
!==========================================================================================!
module meds_pedotransfer
   use meds_kinds,           only : wp
   use meds_constants,       only : grav_head, day_sec, tiny_num
   use meds_water_retention, only : SOIL_RETENTION_VG, soil_theta_from_psi
   implicit none
   private

   public :: soil_texture_t, soil_hydraulics_from_texture

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

   real(wp), parameter :: kpa_per_mpa      = 1000.0_wp

   !----- Ottoni et al. (2019): Ks = a (theta_sat - theta(-33 kPa))^b, Ks in m/day. -------!
   real(wp), parameter :: psi_33kpa        = -0.033_wp / grav_head   !< [m] -33 kPa as water head
   real(wp), parameter :: ottoni_ks_coef   = 19.31_wp                !< [m/day]
   real(wp), parameter :: ottoni_ks_exp    = 1.948_wp                !< [-]
   !----- Density of the mineral grains and of the fragments (Saxton & Rawls 2006). -------!
   real(wp), parameter :: particle_density = 2.65_wp                 !< [g/cm3]

contains

   !---------------------------------------------------------------------------------------!
   ! The van Genuchten curve and Ks of one soil layer. Hodnett & Tomasella fit every       !
   ! parameter as a polynomial in percent, with sand, silt, clay and organic carbon in %,  !
   ! CEC in cmol(c)/kg and bulk density in g/cm3. Alpha comes out per kPa and is turned    !
   ! into per metre of water head, the unit of MEDS's soil potential.                      !
   !---------------------------------------------------------------------------------------!
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
         !----- The regression dips below zero on sandy, alkaline soils of little CEC,     !
         !      which hold no residual water; ED2 floors it the same way. -----------------!
         theta_res = max((22.733_wp - 0.164_wp*sand + 0.235_wp*cec - 0.831_wp*ph + 0.0018_wp*clay**2  &
                          + 0.0026_wp*sand*clay) / 100.0_wp, 0.0_wp)

         !----- Ks of the fine earth from the water that drains between saturation and     !
         !      -33 kPa (Ottoni et al. 2019). MEDS's curve saturates at psi = 0, so        !
         !      that is where the drainable water starts. ---------------------------------!
         theta_33  = soil_theta_from_psi(SOIL_RETENTION_VG, psi_33kpa, theta_sat, theta_res, vg_alpha, vg_n)
         drainable = max(theta_sat - theta_33, tiny_num)
         ksat      = ottoni_ks_coef / day_sec * drainable**ottoni_ks_exp

         !----- Coarse fragments: the layer holds the fine earth's water in the fine       !
         !      earth's share of its volume. Ks follows Saxton & Rawls (2006) eq. 22,      !
         !      written for the fragments' share of the mass, which eq. 19 gives from      !
         !      their share of the volume. ------------------------------------------------!
         rock_volume          = texture%coarse_frag / 100.0_wp
         matrix_density_ratio = bulk_density / particle_density
         rock_mass = rock_volume / max(matrix_density_ratio*(1.0_wp - rock_volume) + rock_volume, tiny_num)
         theta_sat = theta_sat * (1.0_wp - rock_volume)
         theta_res = theta_res * (1.0_wp - rock_volume)
         ksat      = ksat * (1.0_wp - rock_mass) / (1.0_wp - rock_mass*(1.0_wp - 1.5_wp*matrix_density_ratio))
      end associate
   end subroutine soil_hydraulics_from_texture

end module meds_pedotransfer
