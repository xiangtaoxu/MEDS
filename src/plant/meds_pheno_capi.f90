!==========================================================================================!
! meds_pheno_capi -- a thin ISO_C_BINDING shim that exposes the leaf-phenology SIGNAL kernel  !
! to C / Python (ctypes). Like meds_plant_capi it is NOT part of libmeds_plant; it is compiled  !
! only into the optional shared library libmeds_plant_c (CMake -DMEDS_BUILD_PYLIB=ON), via the   !
! src/plant/*_capi.f90 glob, so the core model stays foreign-call-free.                          !
!                                                                                          !
! Four C-interoperable derived types mirror pheno_env_t / pheno_params_t / pheno_state_t /       !
! pheno_out_t (real64 -> c_double, int32 -> c_int; the two logicals hemis_north /                !
! water_use_potential -> c_int 0/1). One exported procedure advances the kernel ONE daily step:  !
!   * meds_phenology_step -- unpack the C structs, call phenology_kernel (which advances the two  !
!       governor drives in state IN PLACE), pack state + the two relative rates back out.         !
! The state struct is intent(inout): the caller keeps it across days (the phenological memory).   !
! The Python side (python/meds/pheno/_ffi.py) mirrors these structs field-for-field.              !
!==========================================================================================!
module meds_pheno_capi
   use iso_c_binding, only : c_double, c_int
   use meds_kinds,    only : wp, ik
   use meds_plant_types, only : pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t
   use meds_phenology,   only : phenology_kernel
   implicit none
   private

   public :: pheno_env_c, pheno_params_c, pheno_state_c, pheno_out_c
   public :: meds_phenology_step

   !----- C-interoperable mirror of pheno_env_t (6 doubles + 2 ints; hemis_north 0/1). -----!
   type, bind(c) :: pheno_env_c
      real(c_double) :: temp_day, soil_temp, avail_water, dmax_leaf_psi, rad, daylength
      integer(c_int) :: doy, hemis_north
   end type pheno_env_c

   !----- C-interoperable mirror of pheno_params_t (masks + rate scales + all cue params). --!
   type, bind(c) :: pheno_params_c
      integer(c_int) :: flush_cue_mask, shed_cue_mask
      real(c_double) :: cue_sharpness, k_flush_max, k_shed_max, tau_flush, tau_shed
      real(c_double) :: gdd_base_temp, chill_base_temp, phen_a, phen_b, phen_c
      real(c_double) :: cold_drop_daylength, cold_drop_soiltemp1, cold_drop_soiltemp2
      integer(c_int) :: water_use_potential
      real(c_double) :: water_off_threshold, water_on_threshold, water_window, water_width
      real(c_double) :: leaf_psi_tlp, low_psi_threshold, high_psi_threshold
      real(c_double) :: photo_crit, photo_slope
      real(c_double) :: light_on_threshold, light_width, light_window
   end type pheno_params_c

   !----- C-interoperable mirror of pheno_state_t (8 doubles; the prognostic memory). ------!
   type, bind(c) :: pheno_state_c
      real(c_double) :: flush_drive, shed_drive, gdd, chill, water_avg, low_psi_days, &
                        high_psi_days, light_avg
   end type pheno_state_c

   !----- C-interoperable mirror of pheno_out_t (2 doubles + 1 int). -----------------------!
   type, bind(c) :: pheno_out_c
      real(c_double) :: leaf_flush_rate, leaf_shed_rate
      integer(c_int) :: cue_limiting
   end type pheno_out_c

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one cohort's phenology ONE step: unpack the C structs, call phenology_kernel     !
   ! (state advanced in place), pack the advanced state + the two relative rates back out.    !
   !---------------------------------------------------------------------------------------!
   subroutine meds_phenology_step(env_c, p_c, dt, state_c, out_c) bind(c, name="meds_phenology_step")
      type(pheno_env_c),    intent(in)    :: env_c
      type(pheno_params_c), intent(in)    :: p_c
      real(c_double), value, intent(in)   :: dt
      type(pheno_state_c),  intent(inout) :: state_c
      type(pheno_out_c),    intent(out)   :: out_c
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: p
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out

      !----- Unpack the environment. -----------------------------------------------------!
      env%temp_day      = env_c%temp_day
      env%soil_temp     = env_c%soil_temp
      env%avail_water   = env_c%avail_water
      env%dmax_leaf_psi = env_c%dmax_leaf_psi
      env%rad           = env_c%rad
      env%daylength     = env_c%daylength
      env%doy           = int(env_c%doy, ik)
      env%hemis_north   = env_c%hemis_north /= 0_c_int

      !----- Unpack the parameters. ------------------------------------------------------!
      p%flush_cue_mask      = int(p_c%flush_cue_mask, ik)
      p%shed_cue_mask       = int(p_c%shed_cue_mask, ik)
      p%cue_sharpness       = p_c%cue_sharpness
      p%k_flush_max         = p_c%k_flush_max
      p%k_shed_max          = p_c%k_shed_max
      p%tau_flush           = p_c%tau_flush
      p%tau_shed            = p_c%tau_shed
      p%gdd_base_temp       = p_c%gdd_base_temp
      p%chill_base_temp     = p_c%chill_base_temp
      p%phen_a              = p_c%phen_a
      p%phen_b              = p_c%phen_b
      p%phen_c              = p_c%phen_c
      p%cold_drop_daylength = p_c%cold_drop_daylength
      p%cold_drop_soiltemp1 = p_c%cold_drop_soiltemp1
      p%cold_drop_soiltemp2 = p_c%cold_drop_soiltemp2
      p%water_use_potential = p_c%water_use_potential /= 0_c_int
      p%water_off_threshold = p_c%water_off_threshold
      p%water_on_threshold  = p_c%water_on_threshold
      p%water_window        = p_c%water_window
      p%water_width         = p_c%water_width
      p%leaf_psi_tlp        = p_c%leaf_psi_tlp
      p%low_psi_threshold   = p_c%low_psi_threshold
      p%high_psi_threshold  = p_c%high_psi_threshold
      p%photo_crit          = p_c%photo_crit
      p%photo_slope         = p_c%photo_slope
      p%light_on_threshold  = p_c%light_on_threshold
      p%light_width         = p_c%light_width
      p%light_window        = p_c%light_window

      !----- Unpack the prognostic state (the caller keeps it across days). ---------------!
      state%flush_drive   = state_c%flush_drive
      state%shed_drive    = state_c%shed_drive
      state%gdd           = state_c%gdd
      state%chill         = state_c%chill
      state%water_avg     = state_c%water_avg
      state%low_psi_days  = state_c%low_psi_days
      state%high_psi_days = state_c%high_psi_days
      state%light_avg     = state_c%light_avg

      call phenology_kernel(env, p, real(dt, wp), state, out)

      !----- Pack the advanced state back. -----------------------------------------------!
      state_c%flush_drive   = state%flush_drive
      state_c%shed_drive    = state%shed_drive
      state_c%gdd           = state%gdd
      state_c%chill         = state%chill
      state_c%water_avg     = state%water_avg
      state_c%low_psi_days  = state%low_psi_days
      state_c%high_psi_days = state%high_psi_days
      state_c%light_avg     = state%light_avg

      !----- Pack the outputs. -----------------------------------------------------------!
      out_c%leaf_flush_rate = out%leaf_flush_rate
      out_c%leaf_shed_rate  = out%leaf_shed_rate
      out_c%cue_limiting    = int(out%cue_limiting, c_int)
   end subroutine meds_phenology_step

end module meds_pheno_capi
