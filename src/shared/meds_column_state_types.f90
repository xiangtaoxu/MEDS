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

   public :: n_soil_layer_max, n_snow_layer_max, N_HYDRO_NODE, LEAF_TEMP_INIT, PSI_INIT
   public :: cas_state_t, soil_column_t, soil_energy_column_t, snow_column_t
   public :: blend_cas, blend_soil_w, blend_soil_e, blend_snow  !< area-weighted mix (patch fusion / disturbance seed)

   integer(ik), parameter :: n_soil_layer_max = 20_ik      !< compile-time soil-column-depth ceiling

   !----- Per-cohort fast state carried on cohort_block (rides the cohort lockstep). N_HYDRO_NODE !
   !      MUST equal meds_plant_types%N_HYDRO (the psi node count); a fresh cohort starts at a     !
   !      mild tension + a neutral leaf temperature (both relax within one fast step).            !
   integer(ik), parameter :: N_HYDRO_NODE   = 3_ik          !< == N_HYDRO (leaf/wood/root psi nodes)
   real(wp),    parameter :: LEAF_TEMP_INIT = 288.15_wp     !< [K]   fresh-cohort leaf temperature
   real(wp),    parameter :: PSI_INIT       = -0.1_wp       !< [MPa] fresh-cohort node water potential

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

   integer(ik), parameter :: n_snow_layer_max = 1_ik        !< MVP single bulk layer; P1 raises to ED2 nzs~8

   !----- Temporary-surface-water / SNOW store: a stacked mass+energy reservoir between the CAS      !
   !      and soil_energy(1). PROGNOSTIC = water-equivalent mass (swe) + EXTENSIVE internal energy    !
   !      (J/m2, unlike soil_energy's J/m3); temperature + liquid fraction are read-offs of           !
   !      uext_to_temp (dry_hcap=0). nlayer=0 is the "no snow" state (store present but empty).        !
   type :: snow_column_t
      real(wp) :: swe(n_snow_layer_max)         = 0.0_wp   !< [kg/m2] water-equivalent mass  (PROGNOSTIC)
      real(wp) :: snow_energy(n_snow_layer_max) = 0.0_wp   !< [J/m2]  extensive internal energy (PROGNOSTIC)
      real(wp) :: snow_depth(n_snow_layer_max)  = 0.0_wp   !< [m]     geometric depth (PROGNOSTIC, = swe/rho_snow in P0)
      real(wp) :: snow_temp(n_snow_layer_max)   = 273.16_wp !< [K]    diagnosed each step (read-off)
      real(wp) :: snow_fliq(n_snow_layer_max)   = 0.0_wp   !< [-]     diagnosed liquid fraction (read-off)
      integer(ik) :: nlayer                     = 0_ik     !< active layer count (0 = no snow)
   end type snow_column_t

   !----- Prognostic per-patch canopy-air-space thermal state (three implicit twins). -------!
   type :: cas_state_t
      real(wp) :: can_enthalpy = 0.0_wp                     !< [J/kg] specific enthalpy (PROGNOSTIC)
      real(wp) :: can_shv      = 0.0_wp                     !< [kg/kg] specific humidity (PROGNOSTIC twin)
      real(wp) :: can_co2      = 400.0_wp                   !< [umol/mol] dry-air CO2 mixing ratio (PROGNOSTIC third twin)
      real(wp) :: can_temp     = 0.0_wp                     !< [K]    diagnosed
      real(wp) :: can_depth    = 20.0_wp                    !< [m]    CAS depth (from canopy height; forcing)
   end type cas_state_t

contains

   !=======================================================================================!
   !  Area-weighted linear mixes: result = w1*a + w2*b. The caller passes NORMALIZED weights !
   !  (w1 = a1/(a1+a2), w2 = a2/(a1+a2)) so an intensive quantity (theta, enthalpy, [J/m3])   !
   !  is conserved on an area basis when two patches fuse or a disturbance gap is carved from  !
   !  its donors. Diagnosed fields (temp/fliq) mix too and are re-diagnosed next fast step.    !
   !=======================================================================================!
   pure function blend_cas(w1, a, w2, b) result(c)
      real(wp),          intent(in) :: w1, w2
      type(cas_state_t), intent(in) :: a, b
      type(cas_state_t)             :: c
      c%can_enthalpy = w1 * a%can_enthalpy + w2 * b%can_enthalpy
      c%can_shv      = w1 * a%can_shv      + w2 * b%can_shv
      c%can_co2      = w1 * a%can_co2      + w2 * b%can_co2
      c%can_temp     = w1 * a%can_temp     + w2 * b%can_temp
      c%can_depth    = w1 * a%can_depth    + w2 * b%can_depth
   end function blend_cas

   pure function blend_soil_w(w1, a, w2, b) result(c)
      real(wp),           intent(in) :: w1, w2
      type(soil_column_t), intent(in) :: a, b
      type(soil_column_t)             :: c
      c%theta     = w1 * a%theta     + w2 * b%theta
      c%w_surface = w1 * a%w_surface + w2 * b%w_surface
      c%w_aquifer = w1 * a%w_aquifer + w2 * b%w_aquifer
      c%z_wt      = w1 * a%z_wt      + w2 * b%z_wt
   end function blend_soil_w

   pure function blend_soil_e(w1, a, w2, b) result(c)
      real(wp),                  intent(in) :: w1, w2
      type(soil_energy_column_t), intent(in) :: a, b
      type(soil_energy_column_t)             :: c
      c%soil_energy = w1 * a%soil_energy + w2 * b%soil_energy
      c%soil_temp   = w1 * a%soil_temp   + w2 * b%soil_temp
      c%soil_fliq   = w1 * a%soil_fliq   + w2 * b%soil_fliq
   end function blend_soil_e

   !----- Area-weighted mix of two snow stores (patch fusion / disturbance seed). The extensive       !
   !      mass + energy blend additively-per-area (like w_surface); depth blends by area weight;       !
   !      temp/fliq are RE-DIAGNOSED by the caller from the blended (swe, snow_energy), never blended.  !
   !      nlayer of the mix = 1 if any blended swe remains, else 0 (caller re-diagnoses / cleans up).   !
   pure function blend_snow(w1, a, w2, b) result(c)
      real(wp),            intent(in) :: w1, w2
      type(snow_column_t), intent(in) :: a, b
      type(snow_column_t)             :: c
      c%swe         = w1 * a%swe         + w2 * b%swe
      c%snow_energy = w1 * a%snow_energy + w2 * b%snow_energy
      c%snow_depth  = w1 * a%snow_depth  + w2 * b%snow_depth
      c%snow_temp   = w1 * a%snow_temp   + w2 * b%snow_temp   ! provisional; caller re-diagnoses from (swe, energy)
      c%snow_fliq   = w1 * a%snow_fliq   + w2 * b%snow_fliq   ! provisional; caller re-diagnoses
      c%nlayer      = max(a%nlayer, b%nlayer)                 ! caller collapses to 0 if blended swe < tiny
   end function blend_snow

end module meds_column_state_types
