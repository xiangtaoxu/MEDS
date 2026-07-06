!==========================================================================================!
! meds_column_state_types -- the PROGNOSTIC per-store column-state derived types of the fast   !
! (sub-daily) biophysics loop: the canopy-air-space thermal twins and the two soil columns      !
! (water + thermal internal energy). Pure DATA (no methods, no hidden state); they depend only   !
! on meds_kinds, so they live in `shared` -- the ONE place reachable by BOTH the biophysics       !
! kernels (which mutate them) AND the demographic state hub (which will OWN them, per-patch, and   !
! thread them through the lockstep reorder). meds_biophysics_types re-exports these names so the   !
! fast kernels compile unchanged; the state hub `use`s this module directly (no biophysics edge).  !
!                                                                                          !
! Fixed-size (n_soil_layer_max) so the soil columns stay allocatable-free and GPU-eligible.     !
!==========================================================================================!
module meds_column_state_types
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: n_soil_layer_max
   public :: cas_state_t, soil_column_t, soil_energy_column_t

   integer(ik), parameter :: n_soil_layer_max = 20_ik      !< compile-time soil-column-depth ceiling

   !----- Prognostic per-patch soil WATER column (the value the hydrology kernel updates). --!
   type :: soil_column_t
      real(wp) :: theta(n_soil_layer_max) = 0.0_wp   !< [m3/m3] volumetric soil moisture (PROGNOSTIC)
      real(wp) :: w_surface = 0.0_wp                 !< [kg/m2] ponded surface water
      real(wp) :: w_aquifer = 0.0_wp                 !< [kg/m2] lumped aquifer store (SOIL_BC_AQUIFER, P2)
      real(wp) :: z_wt      = 0.0_wp                 !< [m] water-table elevation (<= 0; P2)
   end type soil_column_t

   !----- Prognostic per-patch soil THERMAL column (internal energy; temp/fliq diagnosed). --!
   type :: soil_energy_column_t
      real(wp) :: soil_energy(n_soil_layer_max) = 0.0_wp    !< [J/m3] volumetric internal energy (PROGNOSTIC)
      real(wp) :: soil_temp(n_soil_layer_max)   = 0.0_wp    !< [K]    diagnosed each step
      real(wp) :: soil_fliq(n_soil_layer_max)   = 1.0_wp    !< [-]    diagnosed liquid fraction
   end type soil_energy_column_t

   !----- Prognostic per-patch canopy-air-space thermal state (three implicit twins). -------!
   type :: cas_state_t
      real(wp) :: can_enthalpy = 0.0_wp                     !< [J/kg] specific enthalpy (PROGNOSTIC)
      real(wp) :: can_shv      = 0.0_wp                     !< [kg/kg] specific humidity (PROGNOSTIC twin)
      real(wp) :: can_co2      = 400.0_wp                   !< [umol/mol] dry-air CO2 mixing ratio (PROGNOSTIC third twin)
      real(wp) :: can_temp     = 0.0_wp                     !< [K]    diagnosed
      real(wp) :: can_depth    = 20.0_wp                    !< [m]    CAS depth (from canopy height; forcing)
   end type cas_state_t

end module meds_column_state_types
