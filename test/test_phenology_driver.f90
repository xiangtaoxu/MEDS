! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_phenology_driver -- integration test for the slow-loop leaf-phenology WIRING            !
! (meds_vegetation_dynamics.advance_leaf_phenology, the folded phenology driver), as opposed to  !
! the kernel itself (test_plant_phenology). It drives a two-cohort site through one synthetic    !
! Ithaca year, feeding the daily drivers through the site accumulators the fast loop fills, and   !
! checks the two leaf tendencies:                                                                 !
!                                                                                          !
!   1. TEMPERATE DECIDUOUS : warmth x day-length flush, cold x short-day shed -- flushing in      !
!                            mid-summer, senescing in late autumn.                               !
!   2. NO CUES             : a cohort with both masks CUE_NONE stays flushing, never senescing.  !
!   3. NO-TEMPERATURE      : a step with no accumulated air temperature (pheno_tair_n = 0) is    !
!                            skipped.                                                            !
!   4. THE CUE DRIVERS, each against a HAND-COMPUTED value, because nothing else checks that the !
!      numbers the fast loop reduces are the numbers the kernel reads: the warmth and cold sums,  !
!      the predawn-psi sums against the PFT's derived turgor-loss point, the shortwave mean.     !
!   5. ACCEPTANCE          : a drought-deciduous and a light-exchanging PFT run end to end.      !
!==========================================================================================!
program test_phenology_driver
   use meds_kinds,                only : wp, ik
   use meds_config,               only : meds_config_t, pft_leaf_psi_tlp
   use meds_site_state_types,     only : site_t
   use meds_init,                 only : init_bare_ground, add_cohort
   use meds_vegetation_dynamics,  only : advance_leaf_phenology
   use meds_phenology_types,      only : CUE_TEMP, CUE_NONE, CUE_WATER, CUE_LIGHT, LIGHT_RADIATION
   use meds_test_support, only : build_test_config, check_close, check_int, check_true, test_report
   implicit none

   real(wp), parameter :: twopi = 6.283185307179586_wp
   type(meds_config_t) :: cfg
   type(site_t)        :: site
   integer(ik) :: doy
   real(wp)    :: fl_decid_200, sh_decid_200, fl_decid_340, sh_decid_340
   real(wp)    :: fl_none_200, sh_none_200, fl_none_340, sh_none_340, gdd_summer

   !----- Config: PFT 1 temperate deciduous (warmth x day length, cold x short days), the    !
   !       with no cues. This test calls advance_leaf_phenology directly. -------------------!
   cfg = build_test_config()
   cfg%forcing%latitude_deg               = 42.44_wp           ! Ithaca NY (northern hemisphere)
   cfg%pft%pheno_flush_cue_mask(1)        = CUE_TEMP + CUE_LIGHT
   cfg%pft%pheno_shed_cue_mask(1)         = CUE_TEMP + CUE_LIGHT
   cfg%pft%pheno_flush_degree_days(1)     = 92.0_wp
   cfg%pft%pheno_shed_base_temp(1)        = 290.37_wp
   cfg%pft%pheno_shed_degree_days(1)      = 48.0_wp
   cfg%pft%pheno_flush_light_threshold(1) = 10.35_wp
   cfg%pft%pheno_shed_light_threshold(1)  = 9.83_wp

   !----- A site with two cohorts: cohort 1 = PFT 1 (deciduous), cohort 2 = PFT 2 (no cues). -!
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
   call add_cohort(site, cfg, 1_ik, 2_ik, 0.5_wp, 10.0_wp)
   call check_int('two cohorts created', int(site%cohort%n, ik), 2_ik)
   call check_close('cohort 1 born flushing',      site%cohort%leaf_flush_tendency(1), 1.0_wp, 1.0e-9_wp)
   call check_close('cohort 1 born not senescing', site%cohort%leaf_shed_tendency(1),  0.0_wp, 1.0e-9_wp)

   !----- Drive one synthetic year (dt_slow = 1 day). ------------------------------------!
   do doy = 1_ik, 365_ik
      call set_daily_drivers(site, daily_tair(doy), 300.0_wp)
      call advance_leaf_phenology(site, cfg, doy)
      if (doy == 200_ik) then
         fl_decid_200 = site%cohort%leaf_flush_tendency(1) ; sh_decid_200 = site%cohort%leaf_shed_tendency(1)
         fl_none_200  = site%cohort%leaf_flush_tendency(2) ; sh_none_200  = site%cohort%leaf_shed_tendency(2)
         gdd_summer   = site%cohort%growing_degree_days(1)
      end if
      if (doy == 340_ik) then
         fl_decid_340 = site%cohort%leaf_flush_tendency(1) ; sh_decid_340 = site%cohort%leaf_shed_tendency(1)
         fl_none_340  = site%cohort%leaf_flush_tendency(2) ; sh_none_340  = site%cohort%leaf_shed_tendency(2)
      end if
   end do

   !----- 1. Temperate deciduous: flushing in summer, senescing in autumn. -----------------!
   call check_true('deciduous flush tendency HIGH in mid-summer (doy 200)',  fl_decid_200 > 0.9_wp,  fl_decid_200)
   call check_true('deciduous shed tendency  LOW  in mid-summer',            sh_decid_200 < 0.05_wp, sh_decid_200)
   call check_true('deciduous shed tendency  HIGH in late autumn (doy 340)', sh_decid_340 > 0.5_wp,  sh_decid_340)
   call check_true('deciduous flush tendency LOW  in late autumn',           fl_decid_340 < 0.5_wp,  fl_decid_340)
   call check_true('deciduous warmth sum accumulated by summer',             gdd_summer > 92.0_wp,   gdd_summer)

   !----- 2. No cues: flushing, never senescing, all year. ---------------------------------!
   call check_true('no cues: flush tendency 1 in summer', fl_none_200 > 1.0_wp - 1.0e-12_wp)
   call check_true('no cues: shed tendency 0 in summer',  sh_none_200 < 1.0e-12_wp)
   call check_true('no cues: flush tendency 1 in autumn', fl_none_340 > 1.0_wp - 1.0e-12_wp)
   call check_true('no cues: shed tendency 0 in autumn',  sh_none_340 < 1.0e-12_wp)

   !----- 3. A no-temperature step is skipped (tendencies + memory unchanged). -------------!
   block
      real(wp) :: fl_before, sh_before, gdd_before
      fl_before  = site%cohort%leaf_flush_tendency(1)
      sh_before  = site%cohort%leaf_shed_tendency(1)
      gdd_before = site%cohort%growing_degree_days(1)
      call set_daily_drivers(site, 0.0_wp, 0.0_wp)
      site%pheno_tair_n = 0_ik                                    ! no fast sub-steps ran
      call advance_leaf_phenology(site, cfg, 1_ik)
      call check_true('no-temperature step leaves the flush tendency unchanged', &
                      abs(site%cohort%leaf_flush_tendency(1) - fl_before) < tiny(1.0_wp))
      call check_true('no-temperature step leaves the shed tendency unchanged', &
                      abs(site%cohort%leaf_shed_tendency(1) - sh_before) < tiny(1.0_wp))
      call check_true('no-temperature step leaves the warmth sum unchanged', &
                      abs(site%cohort%growing_degree_days(1) - gdd_before) < tiny(1.0_wp))
   end block

   !=== 4. THE CUE DRIVERS, each against a hand-computed value. ============================!

   !----- 4a. Temperature: 10 days at 288.15 K add 10 x (288.15 - 278.15) to the warmth sum;   !
   !          before midsummer the cold sum stays 0, after it 10 days at 280.15 K add         !
   !          10 x (290.37 - 280.15). -------------------------------------------------------!
   block
      call reset_pheno_memory(site)
      do doy = 100_ik, 109_ik
         call set_daily_drivers(site, 288.15_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, doy)
      end do
      call check_close('warmth sum after 10 days 10 K above base', site%cohort%growing_degree_days(1), &
                       100.0_wp, 1.0e-9_wp)
      call check_close('no cold sum before midsummer', site%cohort%cold_degree_days(1), 0.0_wp, 0.0_wp)
      call reset_pheno_memory(site)
      do doy = 250_ik, 259_ik
         call set_daily_drivers(site, 280.15_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, doy)
      end do
      call check_close('cold sum after 10 days below the shed base', site%cohort%cold_degree_days(1), &
                       10.0_wp * (290.37_wp - 280.15_wp), 1.0e-9_wp)
   end block

   !----- 4b. Water: the dry sum adds (psi_tlp - psi_pd) per day, against the PFT's DERIVED     !
   !          turgor-loss point; one wet day adds wet credit without wiping it. -------------!
   block
      real(wp) :: tlp, dry7
      cfg%pft%pheno_flush_cue_mask(1) = CUE_NONE
      cfg%pft%pheno_shed_cue_mask(1)  = CUE_WATER
      tlp = pft_leaf_psi_tlp(cfg, 1_ik)
      call reset_pheno_memory(site)
      site%cohort%dmax_psi_leaf(1) = tlp - 1.0_wp           ! 1 MPa below the TLP
      do doy = 1_ik, 7_ik
         call set_daily_drivers(site, 290.0_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      dry7 = site%cohort%dry_psi_sum(1)
      call check_close('dry sum after 7 days 1 MPa below the TLP', dry7, 7.0_wp, 1.0e-12_wp)
      site%cohort%dmax_psi_leaf(1) = tlp + 0.5_wp           ! one wet day
      call set_daily_drivers(site, 290.0_wp, 300.0_wp)
      call advance_leaf_phenology(site, cfg, 200_ik)
      call check_close('one wet day adds wet credit', site%cohort%wet_psi_sum(1), 0.5_wp, 1.0e-12_wp)
      call check_close('one wet day does not wipe the dry sum', site%cohort%dry_psi_sum(1), dry7, 0.0_wp)
   end block

   !----- 4c. Light on shortwave: x += w*(rad - x) with w = dt/window. From x = 0 with a       !
   !          constant input r and w = 0.1, after n days x = r*(1 - 0.9^n) EXACTLY. ----------!
   block
      real(wp) :: expect
      cfg%pft%pheno_shed_cue_mask(1)  = CUE_LIGHT
      cfg%pft%pheno_light_variable(1) = LIGHT_RADIATION
      cfg%pft%pheno_light_window(1)   = 10.0_wp
      call reset_pheno_memory(site)
      do doy = 1_ik, 5_ik
         call set_daily_drivers(site, 290.0_wp, 400.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      expect = 400.0_wp * (1.0_wp - 0.9_wp**5)
      call check_close('shortwave running mean after 5 days at 400 W/m2',                       &
                       site%cohort%shortwave_mean(1), expect, 1.0e-9_wp)
   end block

   !=== 5. ACCEPTANCE: the drought-deciduous and light-exchanging habits, end to end. =======!

   !----- Drought deciduous (flush and shed both on water): full when watered, senescing under  !
   !      sustained drought, flushing again on rewet. -----------------------------------------!
   block
      real(wp) :: shed_wet, shed_dry, flush_dry, flush_rewet, tlp
      integer(ik) :: d
      cfg%pft%pheno_flush_cue_mask(1) = CUE_WATER
      cfg%pft%pheno_shed_cue_mask(1)  = CUE_WATER
      tlp = pft_leaf_psi_tlp(cfg, 1_ik)
      call reset_pheno_memory(site)
      site%cohort%dmax_psi_leaf(1) = 0.5_wp * tlp               ! well watered
      do d = 1_ik, 60_ik
         call set_daily_drivers(site, 295.0_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      shed_wet = site%cohort%leaf_shed_tendency(1)
      site%cohort%dmax_psi_leaf(1) = tlp - 2.0_wp               ! sustained drought, past the TLP
      do d = 1_ik, 30_ik
         call set_daily_drivers(site, 295.0_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      shed_dry  = site%cohort%leaf_shed_tendency(1)
      flush_dry = site%cohort%leaf_flush_tendency(1)
      site%cohort%dmax_psi_leaf(1) = 0.5_wp * tlp               ! rewet
      do d = 1_ik, 30_ik
         call set_daily_drivers(site, 295.0_wp, 300.0_wp)
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      flush_rewet = site%cohort%leaf_flush_tendency(1)
      call check_true('drought-deciduous: no senescence when watered',  shed_wet < 0.05_wp,   shed_wet)
      call check_true('drought-deciduous: senesces under drought',      shed_dry > 0.8_wp,    shed_dry)
      call check_true('drought-deciduous: flush suppressed in drought', flush_dry < 0.2_wp,   flush_dry)
      call check_true('drought-deciduous: flushes again on rewet',      flush_rewet > 0.8_wp, flush_rewet)
   end block

   !----- Light-driven leaf exchange: flushing throughout, senescence rising with light. ----!
   block
      real(wp) :: shed_dim, shed_bright, flush_bright
      integer(ik) :: d
      cfg%pft%pheno_flush_cue_mask(1)       = CUE_NONE
      cfg%pft%pheno_shed_cue_mask(1)        = CUE_LIGHT
      cfg%pft%pheno_shed_light_threshold(1) = 200.0_wp
      cfg%pft%pheno_shed_light_sharpness(1) = 0.05_wp        ! > 0: bright light triggers
      call reset_pheno_memory(site)
      do d = 1_ik, 60_ik
         call set_daily_drivers(site, 295.0_wp, 60.0_wp)      ! dim
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      shed_dim = site%cohort%leaf_shed_tendency(1)
      call reset_pheno_memory(site)
      do d = 1_ik, 60_ik
         call set_daily_drivers(site, 295.0_wp, 500.0_wp)     ! bright
         call advance_leaf_phenology(site, cfg, 200_ik)
      end do
      shed_bright  = site%cohort%leaf_shed_tendency(1)
      flush_bright = site%cohort%leaf_flush_tendency(1)
      call check_true('leaf exchange: little senescence under dim light', shed_dim < 0.05_wp,     shed_dim)
      call check_true('leaf exchange: senescence rises with light',       shed_bright > 0.8_wp,   shed_bright)
      call check_true('leaf exchange: still flushing while exchanging',   flush_bright > 0.99_wp, flush_bright)
   end block

   call test_report('test_phenology_driver')

contains

   !----- Fill the site accumulators the fast loop fills: the air-temperature sum (not area-  !
   !      weighted, site-uniform forcing) and the area-weighted shortwave; with one patch of   !
   !      area 1 both are just the value itself. ------------------------------------------!
   subroutine set_daily_drivers(site, tair, rad)
      type(site_t), intent(inout) :: site
      real(wp),     intent(in)    :: tair, rad
      site%pheno_tair_sum = tair ; site%pheno_tair_n = 1_ik
      site%pheno_rad_sum  = rad
   end subroutine set_daily_drivers

   !----- Clear every phenology memory so each cue block starts from a known state. ---------!
   subroutine reset_pheno_memory(site)
      type(site_t), intent(inout) :: site
      associate (n => site%cohort%n)
         site%cohort%leaf_flush_tendency(1:n) = 1.0_wp
         site%cohort%leaf_shed_tendency(1:n)  = 0.0_wp
         site%cohort%growing_degree_days(1:n) = 0.0_wp
         site%cohort%cold_degree_days(1:n)    = 0.0_wp
         site%cohort%dry_psi_sum(1:n)         = 0.0_wp
         site%cohort%wet_psi_sum(1:n)         = 0.0_wp
         site%cohort%shortwave_mean(1:n)      = 0.0_wp
      end associate
   end subroutine reset_pheno_memory

   !----- Ithaca-like daily-mean air temperature [K]: ~270 in winter, ~297 in summer. --------!
   pure real(wp) function daily_tair(doy) result(t)
      integer(ik), intent(in) :: doy
      t = 283.15_wp + 14.0_wp * sin(twopi * (real(doy, wp) - 110.0_wp) / 365.0_wp)
   end function daily_tair

end program test_phenology_driver
