!==========================================================================================!
! meds_phenology_driver -- MEDS's slow-loop LEAF-PHENOLOGY driver (the analogue of ED2's        !
! phenology_driv). It lives in src/driver/ (compiled into meds_aux) and advances every cohort's  !
! prognostic phenology memory ONE daily step, writing back the directional status the carbon     !
! leaf-flush gate consumes.                                                                       !
!                                                                                          !
! It is the marshalling seam between the flat cohort state (src/state) and the stateless plant    !
! phenology kernel (src/plant): per cohort it FLATTENS the PFT cue params (cfg%pft%pheno_*) into a  !
! meds_plant pheno_params_t -- exactly as meds_plant_interface's leaf seam flattens the photo       !
! traits -- PACKS the cohort's (gdd, chill) into a pheno_state_t, builds the daily environment       !
! from the site daily-mean air temperature (the fast loop's accumulator) + the site latitude and     !
! day-of-year, calls update_phenology, and UNPACKS the advanced memory + status back to the cohort.   !
!                                                                                          !
! v1 wires the TEMP (GDD/chilling) and PHOTO (photoperiod) cues; validate_config rejects the WATER    !
! and HYDRO cue bits until their soil-water / leaf-psi drivers are threaded in here. The cold-drop      !
! trigger's soil temperature uses the daily air temperature as a documented v1 proxy.                   !
!==========================================================================================!
module meds_phenology_driver
   use meds_kinds,                only : wp, ik
   use meds_constants,            only : day_sec
   use meds_config,               only : meds_config_t
   use meds_core_interface, only : site_t
   use meds_plant_interface,      only : pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t,  &
                                         update_phenology, daylength, CUE_NONE
   implicit none
   private

   public :: leaf_phenology

contains

   !---------------------------------------------------------------------------------------!
   ! Advance the leaf phenology of every cohort over one slow step. `doy` is the day-of-year   !
   ! at the START of the step (used for the thermal-sum season gating + daylength). Reads the   !
   ! site daily-mean air temperature the fast loop accumulated; a no-temperature step (no fast   !
   ! sub-steps ran) is skipped so the memory is never advanced on a bogus 0/0 mean.             !
   !---------------------------------------------------------------------------------------!
   subroutine leaf_phenology(site, cfg, doy)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      integer(ik),         intent(in)    :: doy
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: params
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out
      integer(ik) :: i, pf
      real(wp)    :: dt_days, temp_day, dlen
      logical     :: north

      if (site%pheno_tair_n < 1_ik) return             ! no fast sub-steps this slow step -> no drivers
      dt_days  = cfg%dt_slow / day_sec
      temp_day = site%pheno_tair_sum / real(site%pheno_tair_n, wp)
      north    = cfg%forcing%latitude_deg >= 0.0_wp
      dlen     = daylength(cfg%forcing%latitude_deg, doy)

      do i = 1_ik, site%cohort%n
         pf = site%cohort%pft(i)
         call flatten_pheno_params(cfg, pf, params)
         if (params%cue_mask == CUE_NONE) cycle        ! evergreen PFT: leave status untouched (stays ON)

         !----- Daily environment (v1: air temp drives GDD/chilling AND stands in for the soil-temp  !
         !      cold-drop trigger; WATER/HYDRO drivers are unused because those cue bits are rejected). !
         env%temp_day    = temp_day
         env%soil_temp   = temp_day
         env%daylength   = dlen
         env%doy         = doy
         env%hemis_north = north
         env%avail_water = 0.0_wp
         env%psi_leaf    = 0.0_wp

         !----- Pack the cohort's cue memory, advance, unpack. (Reset the whole state each cohort so  !
         !      the unused WATER/HYDRO accumulators cannot leak across cohorts.)  -----------------!
         state       = pheno_state_t()
         state%gdd   = site%cohort%pheno_gdd(i)
         state%chill = site%cohort%pheno_chill(i)
         call update_phenology(env, params, dt_days, state, out)
         site%cohort%pheno_gdd(i)        = state%gdd
         site%cohort%pheno_chill(i)      = state%chill
         site%cohort%phenology_status(i) = out%phenology_status
      end do
   end subroutine leaf_phenology

   !---------------------------------------------------------------------------------------!
   ! Flatten the per-PFT phenology traits (cfg%pft%pheno_*) into the self-contained kernel     !
   ! param set. The WATER/HYDRO param fields keep their pheno_params_t defaults (their cues are  !
   ! rejected in v1). Mirrors meds_plant_interface's leaf-trait flattening.                     !
   !---------------------------------------------------------------------------------------!
   subroutine flatten_pheno_params(cfg, ipft, p)
      type(meds_config_t),  intent(in)  :: cfg
      integer(ik),          intent(in)  :: ipft
      type(pheno_params_t), intent(out) :: p
      associate (t => cfg%pft)
         p%cue_mask            = t%pheno_cue_mask(ipft)
         p%cue_sharpness       = t%pheno_cue_sharpness(ipft)
         p%phen_on_threshold   = t%pheno_on_threshold(ipft)
         p%phen_off_threshold  = t%pheno_off_threshold(ipft)
         p%gdd_base_temp       = t%pheno_gdd_base_temp(ipft)
         p%chill_base_temp     = t%pheno_chill_base_temp(ipft)
         p%phen_a              = t%pheno_phen_a(ipft)
         p%phen_b              = t%pheno_phen_b(ipft)
         p%phen_c              = t%pheno_phen_c(ipft)
         p%cold_drop_daylength = t%pheno_cold_drop_daylength(ipft)
         p%cold_drop_soiltemp1 = t%pheno_cold_drop_soiltemp1(ipft)
         p%cold_drop_soiltemp2 = t%pheno_cold_drop_soiltemp2(ipft)
         p%photo_crit          = t%pheno_photo_crit(ipft)
         p%photo_slope         = t%pheno_photo_slope(ipft)
      end associate
   end subroutine flatten_pheno_params

end module meds_phenology_driver
