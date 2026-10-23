! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_respiration -- stateless non-leaf maintenance respiration, per plant          !
! [umol CO2/plant/s], at rates referenced to 25 degC (docs/science/plant_respiration.md):  !
!   * stem_maintenance_respiration      -- stem surface and sapwood terms, eq. (3);        !
!   * fine_root_maintenance_respiration -- per unit fine-root carbon, eq. (4), with the    !
!                                          layer-summed scale root_zone_temp_scale.        !
! Leaf Rd is in the leaf solver; growth respiration is in meds_plant_carbon_allocation.    !
!==========================================================================================!
module meds_plant_respiration
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : pi, tiny_num
   use meds_therm_lib,     only : peaked_arrhenius_scale
   implicit none
   private

   public :: stem_maintenance_respiration, fine_root_maintenance_respiration
   public :: root_zone_temp_scale

contains

   !---------------------------------------------------------------------------------------!
   ! Stem maintenance respiration [umol CO2/plant/s], docs/science/plant_respiration.md    !
   ! eq. (3). Elemental: a scalar call does one cohort, an array call a patch, with params !
   ! broadcast. A grass (is_woody = .false.) respires nothing here.                        !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine stem_maintenance_respiration(wood_temp, dbh, height, wai, nplant,          &
                                                         aboveground_frac, is_woody, resp_factor25,     &
                                                         sap_volume, resp_sapwood25, ea, hd, ds, stem_resp)
      real(wp),            intent(in)  :: wood_temp   !< [K]  woody-tissue temperature
      real(wp),            intent(in)  :: dbh         !< [cm] stem diameter at breast height
      real(wp),            intent(in)  :: height      !< [m]  cohort height
      real(wp),            intent(in)  :: wai         !< [m2/m2 ground] wood area index
      real(wp),            intent(in)  :: nplant      !< [plant/m2] stem density
      real(wp),            intent(in)  :: aboveground_frac  !< [--] cohort PFT's aboveground fraction
      logical,             intent(in)  :: is_woody          !< cohort PFT is woody (grass => 0)
      real(wp),            intent(in)  :: resp_factor25     !< [umol CO2/m2 stem/s @25C] cohort PFT's surface rate
      !----- The sapwood term: the living wood respires by its volume, coarse roots included (§2). ---!
      real(wp),            intent(in)  :: sap_volume        !< [m3/plant] sapwood volume
      real(wp),            intent(in)  :: resp_sapwood25    !< [umol CO2/m3 sapwood/s @25C] cohort PFT's sapwood rate
      !----- Leaf Rd's temperature parameters, always in the peaked form (plant_respiration.md §4). -!
      real(wp),            intent(in)  :: ea          !< [J/mol]   activation energy
      real(wp),            intent(in)  :: hd          !< [J/mol]   deactivation energy
      real(wp),            intent(in)  :: ds          !< [J/mol/K] entropy term
      real(wp),            intent(out) :: stem_resp   !< [umol CO2 / plant / s]
      real(wp) :: tscale, stem_area

      if (.not. is_woody) then
         stem_resp = 0.0_wp
         return
      end if

      !----- Peaked temperature response (= 1 at 25 degC) on leaf Rd's parameters. -----------!
      tscale = peaked_arrhenius_scale(1.0_wp, ea, hd, ds, wood_temp)
      !----- Per-plant stem surface: the cylinder plus the WAI branch term (WAI is per ground, !
      !      hence /nplant), over the aboveground fraction to cover the whole woody plant.      !
      stem_area = ( pi * (dbh * 1.0e-2_wp) * height                                             &
                  + pi * wai / max(nplant, tiny_num) ) / max(aboveground_frac, tiny_num)
      !----- The surface and sapwood terms share the temperature response (eq. (3)); a PFT with   !
      !      no sapwood rate respires by its surface alone. --------------------------------------!
      stem_resp = resp_factor25 * tscale * stem_area + resp_sapwood25 * tscale * sap_volume
   end subroutine stem_maintenance_respiration

   !---------------------------------------------------------------------------------------!
   ! Fine-root maintenance respiration [umol CO2/plant/s], docs/science/plant_respiration.md !
   ! eq. (4). Elemental; the patch-uniform tscale_root broadcasts over the cohorts.        !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine fine_root_maintenance_respiration(tscale_root, broot, resp_factor25,      &
                                                              root_resp)
      real(wp),            intent(in)  :: tscale_root !< [-] root-weighted temperature scale (root_zone_temp_scale)
      real(wp),            intent(in)  :: broot       !< [kgC/plant] fine-root biomass
      real(wp),            intent(in)  :: resp_factor25 !< [umol CO2/kgC root/s @25C] cohort PFT's baseline
      real(wp),            intent(out) :: root_resp   !< [umol CO2 / plant / s]
      root_resp = resp_factor25 * tscale_root * broot
   end subroutine fine_root_maintenance_respiration

   !---------------------------------------------------------------------------------------!
   ! The root-weighted temperature scale sum_k root_frac_k * f(T_k): summed over layers,   !
   ! not f at the mean temperature, whose error flips sign with the season                 !
   ! (docs/science/plant_respiration.md §3). Patch-uniform: evaluated once per patch.      !
   !---------------------------------------------------------------------------------------!
   pure function root_zone_temp_scale(soil_temp, root_frac, n, ea, hd, ds) result(tscale)
      real(wp),            intent(in) :: soil_temp(:)  !< [K] per-layer soil temperature
      real(wp),            intent(in) :: root_frac(:)  !< [-] per-layer root fraction (sums to 1)
      integer(ik),         intent(in) :: n             !< active layers
      real(wp),            intent(in) :: ea, hd, ds    !< [J/mol], [J/mol], [J/mol/K] leaf Rd's response
      real(wp)    :: tscale
      integer(ik) :: k
      tscale = 0.0_wp
      do k = 1_ik, n
         tscale = tscale + root_frac(k)                                                          &
                * peaked_arrhenius_scale(1.0_wp, ea, hd, ds, soil_temp(k))
      end do
   end function root_zone_temp_scale

end module meds_plant_respiration
