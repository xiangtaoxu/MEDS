!==========================================================================================!
! meds_surface_optics -- ground (soil / litter / water / snow) reflectance and thermal        !
! emission, the lower boundary condition of the canopy two-stream.                            !
!                                                                                          !
! PLACEHOLDER: this first cut models BARE SOIL only -- a configured per-band soil albedo for    !
! shortwave bands, and reflectance = 1 - emissivity with a blackbody source for the thermal      !
! band. The `surface_state_t` reserves the fields a full implementation will use (soil moisture, !
! litter, standing water, snow), but they are not consulted yet. When soil state exists this is   !
! the one place to grow -- and the place to implement ED2's soil-moisture / colour / snow albedo    !
! CORRECTLY (ED2's radiate_driver.f90:639 has an `albedo_damp_nir = albedo_damp_nir` self-          !
! assignment bug in its bedrock branch; do not reproduce it).                                       !
!==========================================================================================!
module meds_surface_optics
   use meds_kinds,     only : wp, ik
   use meds_constants, only : stefan
   implicit none
   private

   public :: surface_state_t, ground_optics

   !---------------------------------------------------------------------------------------!
   ! Surface state. Only the bare-soil fields are used now; the rest are reserved.          !
   !---------------------------------------------------------------------------------------!
   type :: surface_state_t
      integer(ik)           :: n_band = 0_ik
      real(wp), allocatable :: soil_albedo(:)    !< (band) shortwave soil albedo; unused for emission bands
      real(wp)              :: soil_emiss = 0.96_wp   !< thermal emissivity of the ground
      real(wp)              :: soil_temp  = 298.0_wp  !< [K] ground (skin) temperature
      !----- Reserved for the full implementation (not consulted yet). -------------------!
      real(wp)              :: soil_moisture = 0.0_wp !< [m3/m3] top-layer volumetric water
      real(wp)              :: snow_frac     = 0.0_wp !< [--] snow cover fraction
      real(wp)              :: water_frac    = 0.0_wp !< [--] standing-water fraction
   end type surface_state_t

contains

   !---------------------------------------------------------------------------------------!
   ! Fill the per-band ground reflectance and thermal emission. For a shortwave band the      !
   ! reflectance is the soil albedo and the emission is zero; for a thermal (emission) band    !
   ! the reflectance is 1 - emissivity and the emission is emissivity * sigma * T_soil^4.       !
   !---------------------------------------------------------------------------------------!
   subroutine ground_optics(surf, n_band, has_emission, grnd_refl, grnd_emiss)
      type(surface_state_t), intent(in)  :: surf
      integer(ik),           intent(in)  :: n_band
      logical,               intent(in)  :: has_emission(n_band)
      real(wp),              intent(out) :: grnd_refl(n_band)
      real(wp),              intent(out) :: grnd_emiss(n_band)
      integer(ik) :: b
      do b = 1_ik, n_band
         if (has_emission(b)) then
            grnd_refl(b)  = 1.0_wp - surf%soil_emiss
            grnd_emiss(b) = surf%soil_emiss * stefan * surf%soil_temp ** 4
         else
            grnd_refl(b)  = surf%soil_albedo(b)
            grnd_emiss(b) = 0.0_wp
         end if
      end do
   end subroutine ground_optics

end module meds_surface_optics
