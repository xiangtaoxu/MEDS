!==========================================================================================!
! test_phenology_driver -- integration test for the slow-loop leaf-phenology WIRING            !
! (meds_vegetation_dynamics.advance_leaf_phenology, the folded phenology driver), as opposed to  !
! the stateless kernel (test_plant_phenology). It drives a two-cohort site through one synthetic  !
! Ithaca year, feeding a daily air temperature via the site accumulator the fast loop fills, and  !
! checks the two GOVERNOR drives:                                                                 !
!                                                                                          !
!   1. TEMP-DECIDUOUS : a cohort with flush/shed masks = CUE_TEMP flushes (flush_drive high) in     !
!                       mid-summer and sheds (shed_drive high) in late autumn; GDD accumulates.     !
!   2. EVERGREEN      : a cohort with both masks CUE_NONE holds flush_drive~1, shed_drive~0 all year.!
!   3. NO-TEMPERATURE : a step with no accumulated air temperature (pheno_tair_n = 0) is skipped.    !
!   4. THE FOUR CUE DRIVERS (#150), each against a HAND-COMPUTED value, because nothing else       !
!      checks that the numbers the fast loop reduces are the numbers the kernel was written        !
!      against: the soil-temperature cold-drop trigger (no longer the air-temperature proxy), the  !
!      soil-water running mean, the consecutive-dry-day counter, and the radiation running mean.   !
!==========================================================================================!
program test_phenology_driver
   use meds_kinds,                only : wp, ik
   use meds_config,               only : meds_config_t
   use meds_site_state_types, only : site_t
   use meds_init,                 only : init_bare_ground, add_cohort
   use meds_vegetation_dynamics,  only : advance_leaf_phenology
   use meds_phenology_types, only : CUE_TEMP, CUE_NONE, CUE_WATER, CUE_HYDRO, CUE_LIGHT
   use meds_test_support, only : build_test_config, check_close, check_int, check_true, test_report
   implicit none

   real(wp), parameter :: twopi = 6.283185307179586_wp
   type(meds_config_t) :: cfg
   type(site_t)        :: site
   integer(ik) :: doy
   real(wp)    :: fl_temp_200, sh_temp_200, fl_temp_340, sh_temp_340
   real(wp)    :: fl_ever_200, sh_ever_200, fl_ever_340, sh_ever_340, gdd_summer


   !----- Config: one temperature-deciduous PFT (1), the rest evergreen; phenology is             !
   !       unconditional now (docs/dev_plans/archive/MEDS_SLOW_DYNAMICS_DESIGN.md Part I) -- this test     !
   !       calls advance_leaf_phenology directly, so no config flag is needed to enable it. --------!
   cfg = build_test_config()
   cfg%forcing%latitude_deg      = 42.44_wp            ! Ithaca NY (northern hemisphere)
   cfg%pft%pheno_flush_cue_mask  = CUE_NONE            ! default: permissive flush (evergreen)
   cfg%pft%pheno_shed_cue_mask   = CUE_NONE            ! default: no active shed
   cfg%pft%pheno_flush_cue_mask(1) = CUE_TEMP          ! PFT 1: cold-deciduous, flush on GDD
   cfg%pft%pheno_shed_cue_mask(1)  = CUE_TEMP          ! PFT 1: cold-deciduous, shed on cold-drop

   !----- A site with two cohorts: cohort 1 = PFT 1 (deciduous), cohort 2 = PFT 2 (evergreen). !
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
   call add_cohort(site, cfg, 1_ik, 2_ik, 0.5_wp, 10.0_wp)
   call check_int('two cohorts created', int(site%cohort%n, ik), 2_ik)
   call check_close('cohort 1 born flushing (flush_drive=1)', site%cohort%pheno_flush_drive(1), 1.0_wp, 1.0e-9_wp)
   call check_close('cohort 1 born with no active shed',      site%cohort%pheno_shed_drive(1),  0.0_wp, 1.0e-9_wp)

   !----- Drive one synthetic year (dt_slow = 1 day). Each day set the site daily-mean air temp !
   !      (the fast loop's accumulator) to an Ithaca-like sinusoid, then advance the phenology.  !
   fl_temp_200 = -99.0_wp ; sh_temp_200 = -99.0_wp ; fl_temp_340 = -99.0_wp ; sh_temp_340 = -99.0_wp
   fl_ever_200 = -99.0_wp ; sh_ever_200 = -99.0_wp ; fl_ever_340 = -99.0_wp ; sh_ever_340 = -99.0_wp
   gdd_summer  = 0.0_wp
   do doy = 1_ik, 365_ik
      call set_daily_drivers(site, daily_tair(doy), daily_tsoil(doy), 0.8_wp, 300.0_wp)
      call advance_leaf_phenology(site, cfg, doy)
      if (doy == 200_ik) then
         fl_temp_200 = site%cohort%pheno_flush_drive(1) ; sh_temp_200 = site%cohort%pheno_shed_drive(1)
         fl_ever_200 = site%cohort%pheno_flush_drive(2) ; sh_ever_200 = site%cohort%pheno_shed_drive(2)
         gdd_summer  = site%cohort%pheno_gdd(1)
      end if
      if (doy == 340_ik) then
         fl_temp_340 = site%cohort%pheno_flush_drive(1) ; sh_temp_340 = site%cohort%pheno_shed_drive(1)
         fl_ever_340 = site%cohort%pheno_flush_drive(2) ; sh_ever_340 = site%cohort%pheno_shed_drive(2)
      end if
   end do

   !----- 1. Temperature-deciduous: flushing in summer, shedding in autumn. -----------------!
   call check_true('deciduous flush_drive HIGH in mid-summer (doy 200)', fl_temp_200 > 0.5_wp)
   call check_true('deciduous shed_drive  LOW  in mid-summer',           sh_temp_200 < 0.2_wp)
   call check_true('deciduous shed_drive  HIGH in late autumn (doy 340)', sh_temp_340 > 0.5_wp)
   call check_true('deciduous flush_drive LOW  in late autumn',           fl_temp_340 < 0.5_wp)
   call check_true('deciduous GDD accumulated by summer',                 gdd_summer > 100.0_wp)

   !----- 2. Evergreen cohort: flush_drive ~1, shed_drive ~0 all year. ----------------------!
   call check_true('evergreen flush_drive ~1 in summer', fl_ever_200 > 0.9_wp)
   call check_true('evergreen shed_drive  ~0 in summer', sh_ever_200 < 0.1_wp)
   call check_true('evergreen flush_drive ~1 in autumn', fl_ever_340 > 0.9_wp)
   call check_true('evergreen shed_drive  ~0 in autumn', sh_ever_340 < 0.1_wp)

   !----- 3. A no-temperature step is skipped (drives + memory unchanged). ------------------!
   block
      real(wp) :: fl_before, sh_before, gdd_before
      fl_before  = site%cohort%pheno_flush_drive(1)
      sh_before  = site%cohort%pheno_shed_drive(1)
      gdd_before = site%cohort%pheno_gdd(1)
      call set_daily_drivers(site, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp)
      site%pheno_tair_n = 0_ik                                    ! no fast sub-steps ran
      call advance_leaf_phenology(site, cfg, 1_ik)
      call check_true('no-temperature step leaves flush_drive unchanged', &
                      abs(site%cohort%pheno_flush_drive(1) - fl_before) < tiny(1.0_wp))
      call check_true('no-temperature step leaves shed_drive unchanged', &
                      abs(site%cohort%pheno_shed_drive(1) - sh_before) < tiny(1.0_wp))
      call check_true('no-temperature step leaves GDD unchanged', &
                      abs(site%cohort%pheno_gdd(1) - gdd_before) < tiny(1.0_wp))
   end block

   !=== 4. THE CUE DRIVERS, each against a hand-computed value. ============================!
   !     These exist because #150's drivers are the one part of the phenology chain nothing had
   !     checked: the kernel is tested on cue values fed in directly, and the strategies are
   !     tested end to end for TEMP only. What was never asserted is that the value the fast
   !     loop reduces into site%pheno_*_sum is the value the kernel reads out of pheno_env_t.

   !----- 4a. The cold-drop trigger reads SOIL temperature, not air. A warm-air / cold-soil day  !
   !          must shed; if the driver still passed temp_day this cohort would stay flushed.     !
   block
      real(wp) :: shed_cold_soil, shed_warm_soil
      call reset_pheno_memory(site)
      do doy = 1_ik, 20_ik                       ! warm AIR, cold SOIL -> cold-drop must fire
         call set_daily_drivers(site, 295.0_wp, 270.0_wp, 0.8_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      shed_cold_soil = site%cohort%pheno_shed_drive(1)
      call reset_pheno_memory(site)
      do doy = 1_ik, 20_ik                       ! same air, warm soil -> no cold drop
         call set_daily_drivers(site, 295.0_wp, 295.0_wp, 0.8_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      shed_warm_soil = site%cohort%pheno_shed_drive(1)
      call check_true('cold SOIL under warm air sheds (driver reads soil, not air)',             &
                      shed_cold_soil > 0.8_wp, shed_cold_soil)
      call check_true('warm soil under the same air does not shed',                              &
                      shed_warm_soil < 0.1_wp, shed_warm_soil)
   end block

   !----- 4b. CUE_WATER: the running mean is x += w*(env - x) with w = dt/window. From x = 0     !
   !          with a constant input r and w = 0.1, after n days x = r*(1 - 0.9^n) EXACTLY.       !
   block
      real(wp) :: expect
      cfg%pft%pheno_shed_cue_mask(1) = CUE_WATER
      cfg%pft%pheno_water_window     = 10.0_wp
      call reset_pheno_memory(site)
      do doy = 1_ik, 5_ik
         call set_daily_drivers(site, 290.0_wp, 290.0_wp, 0.60_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      expect = 0.60_wp * (1.0_wp - 0.9_wp**5)
      call check_close('CUE_WATER running mean after 5 days at 0.60', site%cohort%pheno_water_avg(1), &
                       expect, 1.0e-12_wp)
      !----- And it PERSISTS: before #150 the accumulator was a local re-zeroed every day, so     !
      !      this would read 0.060 (one day's worth) instead of 0.246.  -------------------------!
      call check_true('the running mean persisted across days (not re-zeroed)',                  &
                      site%cohort%pheno_water_avg(1) > 0.2_wp, site%cohort%pheno_water_avg(1))
   end block

   !----- 4c. CUE_HYDRO: consecutive days with dmax_psi_leaf below the turgor-loss point. With    !
   !          dt = 1 day, low_psi_days after n dry days is exactly n -- and resets to 0 on one    !
   !          wet day, which is what "consecutive" means and what a daily re-zero could not show. !
   block
      cfg%pft%pheno_shed_cue_mask(1) = CUE_HYDRO
      call reset_pheno_memory(site)
      site%cohort%dmax_psi_leaf(1) = -4.0_wp        ! well below psi_tlp
      do doy = 1_ik, 7_ik
         call set_daily_drivers(site, 290.0_wp, 290.0_wp, 0.3_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      call check_close('CUE_HYDRO dry-day counter after 7 dry days', site%cohort%pheno_low_psi_days(1), &
                       7.0_wp, 1.0e-12_wp)
      site%cohort%dmax_psi_leaf(1) = -0.1_wp        ! one wet day
      call set_daily_drivers(site, 290.0_wp, 290.0_wp, 0.3_wp, 300.0_wp)
      call advance_leaf_phenology(site, cfg, 200_ik)
      call check_close('one wet day resets the CONSECUTIVE dry-day counter',                     &
                       site%cohort%pheno_low_psi_days(1), 0.0_wp, 1.0e-12_wp)
   end block

   !----- 4d. CUE_LIGHT: same exponential mean, on incident shortwave. ----------------------!
   block
      real(wp) :: expect
      cfg%pft%pheno_shed_cue_mask(1) = CUE_LIGHT
      cfg%pft%pheno_light_window     = 10.0_wp
      call reset_pheno_memory(site)
      do doy = 1_ik, 5_ik
         call set_daily_drivers(site, 290.0_wp, 290.0_wp, 0.8_wp, 400.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      expect = 400.0_wp * (1.0_wp - 0.9_wp**5)
      call check_close('CUE_LIGHT running mean after 5 days at 400 W/m2',                        &
                       site%cohort%pheno_light_avg(1), expect, 1.0e-9_wp)
   end block

   !=== 5. ACCEPTANCE (#150): the two strategies that could not be selected before now run, and  !
   !       reproduce the behaviour the design's patterns 3 and 4 describe. =====================!

   !----- Pattern 3, facultative drought-deciduous (flush and shed both CUE_HYDRO): full when     !
   !      watered, sheds under sustained drought, REFLUSHES on rewet. The reflush is the part     !
   !      that needs the persisted counters -- it depends on high_psi_days building back up.      !
   block
      real(wp) :: shed_wet, shed_dry, flush_dry, flush_rewet
      integer(ik) :: d
      cfg%pft%pheno_flush_cue_mask(1) = CUE_HYDRO
      cfg%pft%pheno_shed_cue_mask(1)  = CUE_HYDRO
      call reset_pheno_memory(site)
      site%cohort%dmax_psi_leaf(1) = -0.2_wp                    ! well watered
      do d = 1_ik, 30_ik
         call set_daily_drivers(site, 295.0_wp, 295.0_wp, 0.9_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      shed_wet = site%cohort%pheno_shed_drive(1)
      site%cohort%dmax_psi_leaf(1) = -4.0_wp                    ! sustained drought, past psi_tlp
      do d = 1_ik, 30_ik
         call set_daily_drivers(site, 295.0_wp, 295.0_wp, 0.1_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      shed_dry  = site%cohort%pheno_shed_drive(1)
      flush_dry = site%cohort%pheno_flush_drive(1)
      site%cohort%dmax_psi_leaf(1) = -0.2_wp                    ! rewet
      do d = 1_ik, 30_ik
         call set_daily_drivers(site, 295.0_wp, 295.0_wp, 0.9_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      flush_rewet = site%cohort%pheno_flush_drive(1)
      call check_true('pattern 3: no active shed when watered', shed_wet < 0.1_wp, shed_wet)
      call check_true('pattern 3: sheds under sustained drought', shed_dry > 0.8_wp, shed_dry)
      call check_true('pattern 3: flush suppressed while droughted', flush_dry < 0.2_wp, flush_dry)
      call check_true('pattern 3: REFLUSHES on rewet', flush_rewet > 0.8_wp, flush_rewet)
   end block

   !----- Pattern 4, light-driven leaf-exchanging: flush stays permissive (the fixed high        !
   !      k_flush_max the design decided on), shed RISES WITH LIGHT. ------------------------!
   block
      real(wp) :: shed_dim, shed_bright, flush_bright
      integer(ik) :: d
      cfg%pft%pheno_flush_cue_mask(1) = CUE_NONE                ! permissive flush (design sec 10.1)
      cfg%pft%pheno_shed_cue_mask(1)  = CUE_LIGHT
      cfg%pft%pheno_light_on_threshold = 200.0_wp
      call reset_pheno_memory(site)
      do d = 1_ik, 60_ik
         call set_daily_drivers(site, 295.0_wp, 295.0_wp, 0.8_wp, 60.0_wp)    ! dim
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      shed_dim = site%cohort%pheno_shed_drive(1)
      call reset_pheno_memory(site)
      do d = 1_ik, 60_ik
         call set_daily_drivers(site, 295.0_wp, 295.0_wp, 0.8_wp, 500.0_wp)   ! bright
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      shed_bright  = site%cohort%pheno_shed_drive(1)
      flush_bright = site%cohort%pheno_flush_drive(1)
      call check_true('pattern 4: little shed under dim light', shed_dim < 0.2_wp, shed_dim)
      call check_true('pattern 4: shed rises with light',        shed_bright > 0.8_wp, shed_bright)
      call check_true('pattern 4: flush stays permissive while exchanging leaves',                &
                      flush_bright > 0.9_wp, flush_bright)
   end block

   call test_report('test_phenology_driver')

contains

   !----- Fill EVERY site accumulator the fast loop fills, exactly as the fast loop fills it: the !
   !      air-temperature sum is not area-weighted (site-uniform forcing) while the other three   !
   !      are, so with one patch of area 1 all four are just the value itself. A fixture that     !
   !      sets only the air temperature leaves soil temperature at 0 K, which trips the           !
   !      unconditional cold-soil drop in midsummer -- which is exactly what it did.  ------------!
   subroutine set_daily_drivers(site, tair, tsoil, swater, rad)
      type(site_t), intent(inout) :: site
      real(wp),     intent(in)    :: tair, tsoil, swater, rad
      site%pheno_tair_sum   = tair ; site%pheno_tair_n = 1_ik
      site%pheno_soilt_sum  = tsoil
      site%pheno_swater_sum = swater
      site%pheno_rad_sum    = rad
   end subroutine set_daily_drivers

   !----- Clear every phenology memory so each cue block starts from a known state. ---------!
   subroutine reset_pheno_memory(site)
      type(site_t), intent(inout) :: site
      site%cohort%pheno_flush_drive(1:site%cohort%n)   = 1.0_wp
      site%cohort%pheno_shed_drive(1:site%cohort%n)    = 0.0_wp
      site%cohort%pheno_gdd(1:site%cohort%n)           = 0.0_wp
      site%cohort%pheno_chill(1:site%cohort%n)         = 0.0_wp
      site%cohort%pheno_water_avg(1:site%cohort%n)     = 0.0_wp
      site%cohort%pheno_low_psi_days(1:site%cohort%n)  = 0.0_wp
      site%cohort%pheno_high_psi_days(1:site%cohort%n) = 0.0_wp
      site%cohort%pheno_light_avg(1:site%cohort%n)     = 0.0_wp
   end subroutine reset_pheno_memory

   !----- Soil temperature: the air sinusoid damped and lagged, as a real column would be. ---!
   pure real(wp) function daily_tsoil(doy) result(t)
      integer(ik), intent(in) :: doy
      t = 283.15_wp + 9.0_wp * sin(twopi * (real(doy, wp) - 130.0_wp) / 365.0_wp)
   end function daily_tsoil


   !----- Ithaca-like daily-mean air temperature [K]: ~270 in winter, ~297 in summer. --------!
   pure real(wp) function daily_tair(doy) result(t)
      integer(ik), intent(in) :: doy
      t = 283.15_wp + 14.0_wp * sin(twopi * (real(doy, wp) - 110.0_wp) / 365.0_wp)
   end function daily_tair




end program test_phenology_driver
