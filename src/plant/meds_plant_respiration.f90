!==========================================================================================!
! meds_plant_respiration -- non-leaf autotrophic-respiration COMPUTE kernels (an ED2 port).  !
!                                                                                          !
!   * stem_maintenance_respiration      -- woody-tissue maintenance respiration, ED2 Chambers  !
!       (2004) surface-area form: a 25 degC baseline (optionally DBH-size-scaled) times a peaked  !
!       temperature response, times the per-plant stem surface area (cylinder + WAI branch term). !
!   * fine_root_maintenance_respiration -- fine-root maintenance respiration, ED2 per-broot form: !
!       a 25 degC per-kgC baseline times the peaked temperature response, times fine-root biomass.  !
!   * growth_respiration                -- construction respiration as a fraction of the input       !
!       after-maintenance carbon (ED2 form); the caller supplies npp_in, owns no carbon pool here.     !
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
   use meds_plant_types,   only : wood_env_t, wood_params_t, wood_flux_t,                       &
                                  root_env_t, root_params_t, root_flux_t
   use meds_temp_response, only : peaked_arrhenius_scale
   implicit none
   private

   public :: stem_maintenance_respiration, fine_root_maintenance_respiration
   public :: growth_respiration

contains

   !---------------------------------------------------------------------------------------!
   ! Stem maintenance respiration for one cohort/individual [umol CO2 / plant / s].          !
   ! Grasses (is_woody = .false.) have no stem, so the flux is identically zero.              !
   !---------------------------------------------------------------------------------------!
   subroutine stem_maintenance_respiration(env, params, out)
      type(wood_env_t),    intent(in)  :: env
      type(wood_params_t), intent(in)  :: params
      type(wood_flux_t),   intent(out) :: out
      real(wp) :: srf25, tscale, stem_area

      if (.not. params%is_woody) then
         out%stem_resp = 0.0_wp
         return
      end if

      !----- Size-dependent baseline at 25 degC (scaler = 0 => flat), Chambers et al. 2004. --!
      srf25  = params%stem_resp_factor25 * 10.0_wp ** (params%stem_resp_size_scaler * env%dbh)
      !----- Peaked temperature response (= 1 at 25 degC), shared with leaf Rd. ---------------!
      tscale = peaked_arrhenius_scale(1.0_wp, params%ea, params%hd, params%ds, env%wood_temp)
      !----- Per-plant stem surface area: cylinder lateral area + the WAI branch term, scaled !
      !      by the aboveground structural fraction (ED2). WAI is per-ground => /nplant.        !
      stem_area = ( pi * (env%dbh * 1.0e-2_wp) * env%height                                     &
                  + pi * env%wai / max(env%nplant, tiny(1.0_wp)) ) / params%agf_bs
      out%stem_resp = srf25 * tscale * stem_area
   end subroutine stem_maintenance_respiration

   !---------------------------------------------------------------------------------------!
   ! Fine-root maintenance respiration for one cohort/individual [umol CO2 / plant / s].      !
   ! Uses ONE effective (root-biomass-weighted mean) soil temperature; the layered vertical   !
   ! integral is a caller-side policy (kept outside the kernel). broot = 0 => zero flux.       !
   !---------------------------------------------------------------------------------------!
   subroutine fine_root_maintenance_respiration(env, params, out)
      type(root_env_t),    intent(in)  :: env
      type(root_params_t), intent(in)  :: params
      type(root_flux_t),   intent(out) :: out
      real(wp) :: tscale
      tscale = peaked_arrhenius_scale(1.0_wp, params%ea, params%hd, params%ds, env%soil_temp)
      out%root_resp = params%root_resp_factor25 * tscale * env%broot
   end subroutine fine_root_maintenance_respiration

   !---------------------------------------------------------------------------------------!
   ! Growth (construction) respiration as a fraction of the input after-maintenance carbon    !
   ! (ED2 form). npp_in = GPP - Rleaf - Rm_wood - Rm_root, supplied by the caller; rg carries  !
   ! npp_in's units. NO carbon pool is owned here. True NPP = npp_in - rg.                     !
   !---------------------------------------------------------------------------------------!
   elemental pure function growth_respiration(npp_in, growth_resp_factor) result(rg)
      real(wp), intent(in) :: npp_in              !< [carbon flux / plant] after-maintenance carbon (A)
      real(wp), intent(in) :: growth_resp_factor  !< [--] per-PFT construction-cost fraction
      real(wp)             :: rg                   !< [same units as npp_in] growth respiration
      rg = growth_resp_factor * max(0.0_wp, npp_in)
   end function growth_respiration

end module meds_plant_respiration
