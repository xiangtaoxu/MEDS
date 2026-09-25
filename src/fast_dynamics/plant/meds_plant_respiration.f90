! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_respiration -- non-leaf autotrophic-respiration COMPUTE kernels (an ED2 port).  !
!                                                                                          !
!   * stem_maintenance_respiration      -- woody-tissue maintenance respiration, ED2 Chambers  !
!       (2004) surface-area form: a 25 degC baseline (optionally DBH-size-scaled) times a peaked  !
!       temperature response, times the per-plant stem surface area (cylinder + WAI branch term). !
!   * fine_root_maintenance_respiration -- fine-root maintenance respiration, ED2 per-broot form: !
!       a 25 degC per-kgC baseline times the peaked temperature response, times fine-root biomass.  !
! (Growth/construction respiration lives with the growth it charges, in meds_plant_carbon_allocation.)!
!                                                                                          !
! The maintenance factors are 25 degC-referenced (MEDS's single model-wide reference). ED2 references    !
! stem/root respiration at 15 degC, so seeding a MEDS default from an ED2/Chambers number is a one-time   !
! conversion done at PARAMETER INIT (where PFT params are chosen), NOT here -- the kernels just consume    !
! the 25 degC factor. All fluxes are per plant [umol CO2 / plant / s]; x nplant -> per ground. Stateless   !
! (the reserved t_acclim env fields are unused in v1). Leaf Rd is in the leaf solver, not here.            !
! The public seams are re-exported through meds_fast_config.                                        !
!==========================================================================================!
module meds_plant_respiration
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : pi
   use meds_plant_types, only : wood_params_t, root_params_t
   use meds_temp_response, only : peaked_arrhenius_scale
   implicit none
   private

   public :: stem_maintenance_respiration, fine_root_maintenance_respiration
   public :: root_zone_temp_scale

contains

   !---------------------------------------------------------------------------------------!
   ! Stem maintenance respiration [umol CO2 / plant / s]. `elemental pure` over the per-plant      !
   ! SCALAR inputs (MEDS_NUMERICS_SCOPING.md §11): a scalar call does one cohort; an array call     !
   ! (`wood_temp(:)`, `dbh(:)`, ...) does a whole patch, with the uniform `params` POD broadcast --  !
   ! so no separate batch wrapper is needed (the elemental broadcast IS the batch), and a Python/    !
   ! ctypes wrapper vectorises over numpy arrays. Bit-identical to the former derived-type kernel:   !
   ! same arithmetic, the wood_env_t fields are now bare scalar dummies. Grasses (is_woody=.false.)  !
   ! have no stem => 0. The reserved acclimation temperature (unused v1) is dropped from the args.    !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine stem_maintenance_respiration(wood_temp, dbh, height, wai, nplant,          &
                                                         aboveground_frac, is_woody, resp_factor25,     &
                                                         params, stem_resp)
      real(wp),            intent(in)  :: wood_temp   !< [K]  woody-tissue temperature
      real(wp),            intent(in)  :: dbh         !< [cm] stem diameter at breast height
      real(wp),            intent(in)  :: height      !< [m]  cohort height
      real(wp),            intent(in)  :: wai         !< [m2/m2 ground] wood area index
      real(wp),            intent(in)  :: nplant      !< [plant/m2] stem density
      !----- PER-COHORT, not a run constant: this is the PFT's aboveground fraction of woody      !
      !      carbon, the same trait demography and cohort fusion read. It used to be a field on    !
      !      wood_params_t set from a hard-coded 0.7, so a run whose PFTs differed in allocation   !
      !      used their values everywhere EXCEPT here (issue #128). ------------------------------!
      real(wp),            intent(in)  :: aboveground_frac  !< [--] cohort PFT's aboveground fraction
      logical,             intent(in)  :: is_woody          !< cohort PFT is woody (grass => 0)
      real(wp),            intent(in)  :: resp_factor25     !< [umol CO2/m2 stem/s @25C] cohort PFT's baseline
      type(wood_params_t), intent(in)  :: params      !< run-uniform trait POD (broadcast)
      real(wp),            intent(out) :: stem_resp   !< [umol CO2 / plant / s]
      real(wp) :: srf25, tscale, stem_area

      if (.not. is_woody) then
         stem_resp = 0.0_wp
         return
      end if

      !----- Size-dependent baseline at 25 degC (scaler = 0 => flat), Chambers et al. 2004. --!
      srf25  = resp_factor25 * 10.0_wp ** (params%stem_resp_size_scaler * dbh)
      !----- Peaked temperature response (= 1 at 25 degC), shared with leaf Rd. ---------------!
      tscale = peaked_arrhenius_scale(1.0_wp, params%ea, params%hd, params%ds, wood_temp)
      !----- Per-plant stem surface area: cylinder lateral area + the WAI branch term, scaled !
      !      by the aboveground structural fraction (ED2). WAI is per-ground => /nplant.        !
      stem_area = ( pi * (dbh * 1.0e-2_wp) * height                                             &
                  + pi * wai / max(nplant, tiny(1.0_wp)) ) / max(aboveground_frac, tiny(1.0_wp))
      stem_resp = srf25 * tscale * stem_area
   end subroutine stem_maintenance_respiration

   !---------------------------------------------------------------------------------------!
   ! Fine-root maintenance respiration [umol CO2 / plant / s]. `elemental pure` over the per-plant   !
   ! SCALAR inputs (§11): scalar => one cohort, array => a patch (the root-weighted mean soil_temp    !
   ! is patch-uniform, so a scalar `soil_temp` broadcasts over the `broot(:)` array). broot=0 => 0.   !
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
   ! root_zone_temp_scale -- the root-weighted temperature scale for fine-root maintenance      !
   ! respiration (#178), summed OVER LAYERS rather than evaluated at a mean temperature:         !
   !                                                                                          !
   !     sum_k root_frac_k * f(T_k)     NOT     f( sum_k root_frac_k * T_k )                    !
   !                                                                                          !
   ! The model already resolves a soil temperature per layer; collapsing it to a mean before     !
   ! the response throws that away, and the two differ because f is not linear. It is the same   !
   ! Jensen argument #145 makes for Rh, with one extra turn: the PEAKED form is convex below its !
   ! optimum and CONCAVE near it, so the sign of the error flips with season rather than biasing !
   ! one way. Measured at Ithaca on a 2 m column with beta = 2 rooting: +4.8 % in December,      !
   ! -8.3 % in June, and only -0.17 % in the annual mean -- so this is a SEASONAL correction to  !
   ! root respiration, not an annual-budget one, and reporting it as the latter would understate !
   ! it by a factor of thirty.                                                                   !
   !                                                                                          !
   ! Patch-uniform, so the caller evaluates it ONCE and broadcasts it over the cohort array --   !
   ! which is also why the respiration kernel above now takes a scale rather than a temperature. !
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
