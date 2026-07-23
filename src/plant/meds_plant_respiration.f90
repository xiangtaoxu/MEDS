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
! The public seams are re-exported through meds_plant_interface.                                        !
!==========================================================================================!
module meds_plant_respiration
   use meds_kinds,         only : wp
   use meds_constants,     only : pi
   use meds_plant_types,   only : wood_params_t, root_params_t
   use meds_temp_response, only : peaked_arrhenius_scale
   implicit none
   private

   public :: stem_maintenance_respiration, fine_root_maintenance_respiration

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
   elemental pure subroutine stem_maintenance_respiration(wood_temp, dbh, height, wai, nplant, params, stem_resp)
      real(wp),            intent(in)  :: wood_temp   !< [K]  woody-tissue temperature
      real(wp),            intent(in)  :: dbh         !< [cm] stem diameter at breast height
      real(wp),            intent(in)  :: height      !< [m]  cohort height
      real(wp),            intent(in)  :: wai         !< [m2/m2 ground] wood area index
      real(wp),            intent(in)  :: nplant      !< [plant/m2] stem density
      type(wood_params_t), intent(in)  :: params      !< run-uniform trait POD (broadcast)
      real(wp),            intent(out) :: stem_resp   !< [umol CO2 / plant / s]
      real(wp) :: srf25, tscale, stem_area

      if (.not. params%is_woody) then
         stem_resp = 0.0_wp
         return
      end if

      !----- Size-dependent baseline at 25 degC (scaler = 0 => flat), Chambers et al. 2004. --!
      srf25  = params%stem_resp_factor25 * 10.0_wp ** (params%stem_resp_size_scaler * dbh)
      !----- Peaked temperature response (= 1 at 25 degC), shared with leaf Rd. ---------------!
      tscale = peaked_arrhenius_scale(1.0_wp, params%ea, params%hd, params%ds, wood_temp)
      !----- Per-plant stem surface area: cylinder lateral area + the WAI branch term, scaled !
      !      by the aboveground structural fraction (ED2). WAI is per-ground => /nplant.        !
      stem_area = ( pi * (dbh * 1.0e-2_wp) * height                                             &
                  + pi * wai / max(nplant, tiny(1.0_wp)) ) / params%agf_bs
      stem_resp = srf25 * tscale * stem_area
   end subroutine stem_maintenance_respiration

   !---------------------------------------------------------------------------------------!
   ! Fine-root maintenance respiration [umol CO2 / plant / s]. `elemental pure` over the per-plant   !
   ! SCALAR inputs (§11): scalar => one cohort, array => a patch (the root-weighted mean soil_temp    !
   ! is patch-uniform, so a scalar `soil_temp` broadcasts over the `broot(:)` array). broot=0 => 0.   !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine fine_root_maintenance_respiration(soil_temp, broot, params, root_resp)
      real(wp),            intent(in)  :: soil_temp   !< [K] effective (root-weighted mean) soil temperature
      real(wp),            intent(in)  :: broot       !< [kgC/plant] fine-root biomass
      type(root_params_t), intent(in)  :: params      !< run-uniform trait POD (broadcast)
      real(wp),            intent(out) :: root_resp   !< [umol CO2 / plant / s]
      real(wp) :: tscale
      tscale = peaked_arrhenius_scale(1.0_wp, params%ea, params%hd, params%ds, soil_temp)
      root_resp = params%root_resp_factor25 * tscale * broot
   end subroutine fine_root_maintenance_respiration

end module meds_plant_respiration
