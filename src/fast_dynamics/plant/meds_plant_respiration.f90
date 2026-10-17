! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_respiration -- non-leaf autotrophic-respiration COMPUTE kernels (an ED2 port; !
! docs/science/plant_respiration.md).                                                      !
!                                                                                          !
!   * stem_maintenance_respiration -- woody-tissue maintenance respiration, eq. (3): a     !
!       peaked temperature response times two terms, each with its own 25 degC rate -- the !
!       ED2 Chambers (2004) surface-area form on the stem                                  !
!       surface (cylinder + WAI branch term), and a sapwood term on the sapwood volume     !
!       (the living wood; Ryan 1990). A PFT sets either rate to 0 to use the other basis   !
!       alone.                                                                             !
!   * fine_root_maintenance_respiration -- fine-root maintenance respiration, eq. (4), ED2 !
!       per-broot form: a 25 degC per-kgC baseline times the layer-summed temperature      !
!       scale (root_zone_temp_scale), times fine-root biomass.                             !
! (Growth/construction respiration lives with the growth it charges, in                    !
! meds_plant_carbon_allocation.)                                                           !
!                                                                                          !
! The maintenance factors are 25 degC-referenced (MEDS's single model-wide reference). ED2 !
! references stem/root respiration at 15 degC, so seeding a MEDS default from an           !
! ED2/Chambers number is a one-time conversion done at PARAMETER INIT (where PFT params    !
! are chosen), NOT here -- the kernels just consume the 25 degC factor. All fluxes are per !
! plant [umol CO2 / plant / s]; x nplant -> per ground. Stateless. Leaf Rd is in the leaf  !
! solver, not here. The fast-loop pre-pass (meds_fast_prepass) is the caller.              !
!==========================================================================================!
module meds_plant_respiration
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : pi, tiny_num
   use meds_plant_types, only : wood_params_t, root_params_t
   use meds_temp_response, only : peaked_arrhenius_scale
   implicit none
   private

   public :: stem_maintenance_respiration, fine_root_maintenance_respiration
   public :: root_zone_temp_scale

contains

   !---------------------------------------------------------------------------------------!
   ! Stem maintenance respiration [umol CO2 / plant / s],                                  !
   ! docs/science/plant_respiration.md eq. (3). `elemental pure` over the per-plant SCALAR !
   ! inputs (MEDS_NUMERICS_SCOPING.md §11): a scalar call does one cohort; an array call   !
   ! (`wood_temp(:)`, `dbh(:)`, ...) does a whole patch, with the uniform `params` POD     !
   ! broadcast -- so no separate batch wrapper is needed (the elemental broadcast IS the   !
   ! batch), and a Python/ctypes wrapper vectorises over numpy arrays. Grasses             !
   ! (is_woody=.false.) have no stem => 0.                                                 !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine stem_maintenance_respiration(wood_temp, dbh, height, wai, nplant,          &
                                                         aboveground_frac, is_woody, resp_factor25,     &
                                                         sap_volume, resp_sapwood25, params, stem_resp)
      real(wp),            intent(in)  :: wood_temp   !< [K]  woody-tissue temperature
      real(wp),            intent(in)  :: dbh         !< [cm] stem diameter at breast height
      real(wp),            intent(in)  :: height      !< [m]  cohort height
      real(wp),            intent(in)  :: wai         !< [m2/m2 ground] wood area index
      real(wp),            intent(in)  :: nplant      !< [plant/m2] stem density
      !----- PER-COHORT, not a run constant: this is the PFT's aboveground fraction of woody     !
      !      carbon, the same trait demography and cohort fusion read, so a run whose PFTs       !
      !      differ in allocation uses their values here too. -----------------------------------!
      real(wp),            intent(in)  :: aboveground_frac  !< [--] cohort PFT's aboveground fraction
      logical,             intent(in)  :: is_woody          !< cohort PFT is woody (grass => 0)
      real(wp),            intent(in)  :: resp_factor25     !< [umol CO2/m2 stem/s @25C] cohort PFT's surface rate
      !----- The sapwood term: the living wood respires by its volume, so a small stem, all of it  !
      !      sapwood but little of it per m2 of surface, pays less than the surface form charges   !
      !      it, and a large one more. The volume includes the coarse roots, as the surface does. --!
      real(wp),            intent(in)  :: sap_volume        !< [m3/plant] sapwood volume
      real(wp),            intent(in)  :: resp_sapwood25    !< [umol CO2/m3 sapwood/s @25C] cohort PFT's sapwood rate
      type(wood_params_t), intent(in)  :: params      !< run-uniform trait POD (broadcast)
      real(wp),            intent(out) :: stem_resp   !< [umol CO2 / plant / s]
      real(wp) :: tscale, stem_area

      if (.not. is_woody) then
         stem_resp = 0.0_wp
         return
      end if

      !----- Peaked temperature response (= 1 at 25 degC), shared with leaf Rd. ---------------!
      tscale = peaked_arrhenius_scale(1.0_wp, params%ea, params%hd, params%ds, wood_temp)
      !----- Per-plant stem surface area: cylinder lateral area + the WAI branch term, scaled !
      !      by the aboveground structural fraction (ED2). WAI is per-ground => /nplant.        !
      stem_area = ( pi * (dbh * 1.0e-2_wp) * height                                             &
                  + pi * wai / max(nplant, tiny_num) ) / max(aboveground_frac, tiny_num)
      !----- The surface and sapwood terms share the temperature response (eq. (3)); a PFT with   !
      !      no sapwood rate respires by its surface alone. --------------------------------------!
      stem_resp = resp_factor25 * tscale * stem_area + resp_sapwood25 * tscale * sap_volume
   end subroutine stem_maintenance_respiration

   !---------------------------------------------------------------------------------------!
   ! Fine-root maintenance respiration [umol CO2 / plant / s],                             !
   ! docs/science/plant_respiration.md eq. (4). `elemental pure` over the per-plant SCALAR !
   ! inputs (MEDS_NUMERICS_SCOPING.md §11): scalar => one cohort, array => a patch (the    !
   ! root-zone temperature scale is patch-uniform, so a scalar `tscale_root` broadcasts    !
   ! over the `broot(:)` array). broot=0 => 0.                                             !
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
   ! root_zone_temp_scale -- the root-weighted temperature scale for fine-root maintenance !
   ! respiration, summed OVER LAYERS rather than evaluated at a mean temperature:          !
   !                                                                                       !
   !     sum_k root_frac_k * f(T_k)     NOT     f( sum_k root_frac_k * T_k )               !
   !                                                                                       !
   ! f is nonlinear, and the peaked form is convex below its optimum and concave near it,  !
   ! so a mean temperature would make a SEASONAL error that flips sign                     !
   ! (docs/science/plant_respiration.md §3, eq. (4)). Patch-uniform, so the caller         !
   ! evaluates it ONCE and broadcasts it over the cohort array.                            !
   !---------------------------------------------------------------------------------------!
   pure function root_zone_temp_scale(soil_temp, root_frac, n, params) result(tscale)
      real(wp),            intent(in) :: soil_temp(:)  !< [K] per-layer soil temperature
      real(wp),            intent(in) :: root_frac(:)  !< [-] per-layer root fraction (sums to 1)
      integer(ik),         intent(in) :: n             !< active layers
      type(root_params_t), intent(in) :: params
      real(wp)    :: tscale
      integer(ik) :: k
      tscale = 0.0_wp
      do k = 1_ik, n
         tscale = tscale + root_frac(k)                                                          &
                * peaked_arrhenius_scale(1.0_wp, params%ea, params%hd, params%ds, soil_temp(k))
      end do
   end function root_zone_temp_scale

end module meds_plant_respiration
