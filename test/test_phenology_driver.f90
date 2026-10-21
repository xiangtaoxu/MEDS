! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_phenology_driver -- integration test for the slow-loop leaf-phenology WIRING            !
! (meds_vegetation_dynamics.advance_leaf_phenology, the phenology driver), as opposed to         !
! the kernel itself (test_plant_phenology). It drives a two-cohort site through one synthetic    !
! Ithaca year, feeding the daily drivers through the site accumulators the fast loop fills, and   !
! checks the two leaf tendencies:                                                                 !
!                                                                                          !
!   1. TEMPERATE DECIDUOUS : warmth x light flush, cold x short-day shed (the hours of light are  !
!                            the day length here, as a low par_min counts) -- flushing in      !
!                            mid-summer, senescing in late autumn.                               !
!   2. NO CUES             : a cohort with both masks CUE_NONE stays flushing, never senescing.  !
!   3. NO-TEMPERATURE      : a step with no accumulated air temperature (pheno_tair_n = 0) is    !
!                            skipped.                                                            !
!   4. THE CUE DRIVERS (docs/science/plant_phenology.md eqs. (3)-(5)), each against a            !
!      HAND-COMPUTED value, because nothing else checks that the numbers the fast loop reduces   !
!      are the numbers the kernel reads: the warmth and cold sums, the predawn-psi sums against  !
!      the PFT's derived turgor-loss point, each cohort's own light mean, and the light memory a !
!      new cohort starts from.                                                                   !
!   5. ACCEPTANCE          : a drought-deciduous and a light-exchanging PFT run end to end.      !
!==========================================================================================!
program test_phenology_driver
   use meds_kinds,                only : wp, ik
   use meds_config,               only : meds_config_t, pft_leaf_psi_tlp
   use meds_polygon_state_vars,   only : polygon_t
   use meds_init_bare_ground, only : init_bare_ground, add_cohort
   use meds_vegetation_dynamics,  only : advance_leaf_phenology
   use meds_constants,            only : day_sec
   use meds_time,                 only : daylength
   use meds_phenology_types,      only : CUE_TEMP, CUE_NONE, CUE_WATER, CUE_LIGHT
   use meds_test_support, only : build_test_config, check_close, check_int, check_true, test_report
   implicit none

   real(wp), parameter :: TWOPI = 6.283185307179586_wp
   type(meds_config_t) :: cfg
   type(polygon_t)     :: poly
   integer(ik) :: doy
   real(wp)    :: fl_decid_200, sh_decid_200, fl_decid_340, sh_decid_340
   real(wp)    :: fl_none_200, sh_none_200, fl_none_340, sh_none_340, gdd_summer

   !----- Config: PFT 1 temperate deciduous (warmth x light, cold x short days), the rest   !
   !       with no cues. This test calls advance_leaf_phenology directly. ------------------!
   cfg = build_test_config()
   cfg%forcing%latitude               = 42.44_wp           ! Ithaca NY (northern hemisphere)
   cfg%pft%pheno(1)%flush_cue_mask        = CUE_TEMP + CUE_LIGHT
   cfg%pft%pheno(1)%shed_cue_mask         = CUE_TEMP + CUE_LIGHT
   cfg%pft%pheno(1)%flush_degree_days     = 92.0_wp
   cfg%pft%pheno(1)%shed_base_temp        = 290.37_wp
   cfg%pft%pheno(1)%shed_degree_days      = 48.0_wp
   cfg%pft%pheno(1)%flush_light_hours     = 10.35_wp
   cfg%pft%pheno(1)%shed_light_hours      = 9.83_wp

   !----- A site with two cohorts: cohort 1 = PFT 1 (deciduous), cohort 2 = PFT 2 (no cues). -!
   call init_bare_ground(poly, cfg, 1_ik)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
   call add_cohort(poly, cfg, 1_ik, 2_ik, 0.5_wp, 10.0_wp)
   call check_int('two cohorts created', int(poly%cohort%n, ik), 2_ik)
   call check_close('cohort 1 born flushing',      poly%cohort%leaf_flush_tendency(1), 1.0_wp, 1.0e-9_wp)
   call check_close('cohort 1 born not senescing', poly%cohort%leaf_shed_tendency(1),  0.0_wp, 1.0e-9_wp)

   !----- Drive one synthetic year (dt_slow = 1 day). ------------------------------------!
   do doy = 1_ik, 365_ik
      call set_daily_drivers(poly, daily_tair(doy), daylength(cfg%forcing%latitude, doy))
      call advance_leaf_phenology(poly, cfg, doy)
      if (doy == 200_ik) then
         fl_decid_200 = poly%cohort%leaf_flush_tendency(1) ; sh_decid_200 = poly%cohort%leaf_shed_tendency(1)
         fl_none_200  = poly%cohort%leaf_flush_tendency(2) ; sh_none_200  = poly%cohort%leaf_shed_tendency(2)
         gdd_summer   = poly%cohort%growing_degree_days(1)
      end if
      if (doy == 340_ik) then
         fl_decid_340 = poly%cohort%leaf_flush_tendency(1) ; sh_decid_340 = poly%cohort%leaf_shed_tendency(1)
         fl_none_340  = poly%cohort%leaf_flush_tendency(2) ; sh_none_340  = poly%cohort%leaf_shed_tendency(2)
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
      fl_before  = poly%cohort%leaf_flush_tendency(1)
      sh_before  = poly%cohort%leaf_shed_tendency(1)
      gdd_before = poly%cohort%growing_degree_days(1)
      call set_daily_drivers(poly, 0.0_wp, 0.0_wp)
      poly%site%pheno_tair_n(1) = 0_ik                            ! no fast sub-steps ran
      call advance_leaf_phenology(poly, cfg, 1_ik)
      call check_true('no-temperature step leaves the flush tendency unchanged', &
                      abs(poly%cohort%leaf_flush_tendency(1) - fl_before) < tiny(1.0_wp))
      call check_true('no-temperature step leaves the shed tendency unchanged', &
                      abs(poly%cohort%leaf_shed_tendency(1) - sh_before) < tiny(1.0_wp))
      call check_true('no-temperature step leaves the warmth sum unchanged', &
                      abs(poly%cohort%growing_degree_days(1) - gdd_before) < tiny(1.0_wp))
   end block

   !=== 4. THE CUE DRIVERS, each against a hand-computed value. ============================!

   !----- 4a. Temperature: 10 days at 288.15 K add 10 x (288.15 - 278.15) to the warmth sum;   !
   !          before midsummer the cold sum stays 0, after it 10 days at 280.15 K add         !
   !          10 x (290.37 - 280.15). -------------------------------------------------------!
   block
      call reset_pheno_memory(poly)
      do doy = 100_ik, 109_ik
         call set_daily_drivers(poly, 288.15_wp, 12.0_wp)
         call advance_leaf_phenology(poly, cfg, doy)
      end do
      call check_close('warmth sum after 10 days 10 K above base', poly%cohort%growing_degree_days(1), &
                       100.0_wp, 1.0e-9_wp)
      call check_close('no cold sum before midsummer', poly%cohort%cold_degree_days(1), 0.0_wp, 0.0_wp)
      call reset_pheno_memory(poly)
      do doy = 250_ik, 259_ik
         call set_daily_drivers(poly, 280.15_wp, 12.0_wp)
         call advance_leaf_phenology(poly, cfg, doy)
      end do
      call check_close('cold sum after 10 days below the shed base', poly%cohort%cold_degree_days(1), &
                       10.0_wp * (290.37_wp - 280.15_wp), 1.0e-9_wp)
   end block

   !----- 4b. Water: the dry sum adds (psi_tlp - psi_pd) per day, against the PFT's DERIVED     !
   !          turgor-loss point; one wet day adds wet credit without wiping it. -------------!
   block
      real(wp) :: tlp, dry7
      cfg%pft%pheno(1)%flush_cue_mask = CUE_NONE
      cfg%pft%pheno(1)%shed_cue_mask  = CUE_WATER
      tlp = pft_leaf_psi_tlp(cfg, 1_ik)
      call reset_pheno_memory(poly)
      poly%cohort%dmax_psi_leaf(1) = tlp - 1.0_wp           ! 1 MPa below the TLP
      do doy = 1_ik, 7_ik
         call set_daily_drivers(poly, 290.0_wp, 12.0_wp)
         call advance_leaf_phenology(poly, cfg, 200_ik)
      end do
      dry7 = poly%cohort%dry_psi_sum(1)
      call check_close('dry sum after 7 days 1 MPa below the TLP', dry7, 7.0_wp, 1.0e-12_wp)
      poly%cohort%dmax_psi_leaf(1) = tlp + 0.5_wp           ! one wet day
      call set_daily_drivers(poly, 290.0_wp, 12.0_wp)
      call advance_leaf_phenology(poly, cfg, 200_ik)
      call check_close('one wet day adds wet credit', poly%cohort%wet_psi_sum(1), 0.5_wp, 1.0e-12_wp)
      call check_close('one wet day does not wipe the dry sum', poly%cohort%dry_psi_sum(1), dry7, 0.0_wp)
   end block

   !----- 4c. Light: x += w*(h - x) with w = dt/window. From x = 0 with a constant input h     !
   !          and w = 0.1, after n days x = h*(1 - 0.9^n) EXACTLY. Each cohort reads its OWN   !
   !          hours of light (light_hours_accum, counted at its top): a shaded cohort's mean   !
   !          stays lower. A cohort with no light memory starts from its first day's hours. ----!
   block
      real(wp) :: expect
      cfg%pft%pheno(1:2)%shed_cue_mask = CUE_LIGHT
      cfg%pft%pheno(1:2)%light_window  = 10.0_wp
      call reset_pheno_memory(poly)
      call set_daily_drivers(poly, 290.0_wp, 9.0_wp)
      call advance_leaf_phenology(poly, cfg, 200_ik)
      call check_close('a cohort with no light memory starts from its first day', &
                       poly%cohort%light_hours_mean(1), 9.0_wp, 0.0_wp)
      poly%cohort%light_hours_mean(1:2) = 0.0_wp
      do doy = 1_ik, 5_ik
         call set_daily_drivers(poly, 290.0_wp, 12.0_wp)
         poly%cohort%light_hours_accum(2) = 3.0_wp * cfg%dt_slow / day_sec    ! cohort 2 in shade
         call advance_leaf_phenology(poly, cfg, 200_ik)
      end do
      expect = 12.0_wp * (1.0_wp - 0.9_wp**5)
      call check_close('light running mean after 5 days of 12 h',                                &
                       poly%cohort%light_hours_mean(1), expect, 1.0e-9_wp)
      call check_close('a shaded cohort reads its own hours of light (3 h)',                     &
                       poly%cohort%light_hours_mean(2), expect / 4.0_wp, 1.0e-9_wp)
      cfg%pft%pheno(2)%shed_cue_mask = CUE_NONE
   end block

   !=== 5. ACCEPTANCE: the drought-deciduous and light-exchanging habits, end to end. =======!

   !----- Drought deciduous (flush and shed both on water): full when watered, senescing under  !
   !      sustained drought, flushing again on rewet. -----------------------------------------!
   block
      real(wp) :: shed_wet, shed_dry, flush_dry, flush_rewet, tlp
      integer(ik) :: d
      cfg%pft%pheno(1)%flush_cue_mask = CUE_WATER
      cfg%pft%pheno(1)%shed_cue_mask  = CUE_WATER
      tlp = pft_leaf_psi_tlp(cfg, 1_ik)
      call reset_pheno_memory(poly)
      poly%cohort%dmax_psi_leaf(1) = 0.5_wp * tlp               ! well watered
      do d = 1_ik, 60_ik
         call set_daily_drivers(poly, 295.0_wp, 12.0_wp)
         call advance_leaf_phenology(poly, cfg, 200_ik)
      end do
      shed_wet = poly%cohort%leaf_shed_tendency(1)
      poly%cohort%dmax_psi_leaf(1) = tlp - 2.0_wp               ! sustained drought, past the TLP
      do d = 1_ik, 30_ik
         call set_daily_drivers(poly, 295.0_wp, 12.0_wp)
         call advance_leaf_phenology(poly, cfg, 200_ik)
      end do
      shed_dry  = poly%cohort%leaf_shed_tendency(1)
      flush_dry = poly%cohort%leaf_flush_tendency(1)
      poly%cohort%dmax_psi_leaf(1) = 0.5_wp * tlp               ! rewet
      do d = 1_ik, 30_ik
         call set_daily_drivers(poly, 295.0_wp, 12.0_wp)
         call advance_leaf_phenology(poly, cfg, 200_ik)
      end do
      flush_rewet = poly%cohort%leaf_flush_tendency(1)
      call check_true('drought-deciduous: no senescence when watered',  shed_wet < 0.05_wp,   shed_wet)
      call check_true('drought-deciduous: senesces under drought',      shed_dry > 0.8_wp,    shed_dry)
      call check_true('drought-deciduous: flush suppressed in drought', flush_dry < 0.2_wp,   flush_dry)
      call check_true('drought-deciduous: flushes again on rewet',      flush_rewet > 0.8_wp, flush_rewet)
   end block

   !----- Light-driven leaf exchange: flushing throughout, senescence rising with light. ----!
   block
      real(wp) :: shed_dim, shed_bright, flush_bright
      integer(ik) :: d
      cfg%pft%pheno(1)%flush_cue_mask       = CUE_NONE
      cfg%pft%pheno(1)%shed_cue_mask        = CUE_LIGHT
      cfg%pft%pheno(1)%shed_light_hours     = 6.0_wp
      cfg%pft%pheno(1)%shed_light_sharpness = 2.0_wp         ! > 0: many bright hours trigger
      call reset_pheno_memory(poly)
      do d = 1_ik, 60_ik
         call set_daily_drivers(poly, 295.0_wp, 2.0_wp)       ! few bright hours
         call advance_leaf_phenology(poly, cfg, 200_ik)
      end do
      shed_dim = poly%cohort%leaf_shed_tendency(1)
      call reset_pheno_memory(poly)
      do d = 1_ik, 60_ik
         call set_daily_drivers(poly, 295.0_wp, 10.0_wp)      ! many bright hours
         call advance_leaf_phenology(poly, cfg, 200_ik)
      end do
      shed_bright  = poly%cohort%leaf_shed_tendency(1)
      flush_bright = poly%cohort%leaf_flush_tendency(1)
      call check_true('leaf exchange: little senescence under dim light', shed_dim < 0.05_wp,     shed_dim)
      call check_true('leaf exchange: senescence rises with light',       shed_bright > 0.8_wp,   shed_bright)
      call check_true('leaf exchange: still flushing while exchanging',   flush_bright > 0.99_wp, flush_bright)
   end block

   call test_report('test_phenology_driver')

contains

   !----- Fill what the fast loop fills: the site air-temperature sum and, per cohort, the time !
   !      [h] its top's PAR exceeded par_min over the step -- here `hours` a day at every cohort, !
   !      so the driver's per-day value recovers it exactly. ------------------------------------!
   subroutine set_daily_drivers(poly, tair, hours)
      type(polygon_t), intent(inout) :: poly
      real(wp),     intent(in)    :: tair, hours
      poly%site%pheno_tair_sum(1) = tair
      poly%site%pheno_tair_n(1) = 1_ik
      poly%cohort%light_hours_accum(1:poly%cohort%n) = hours * cfg%dt_slow / day_sec
   end subroutine set_daily_drivers

   !----- Clear every phenology memory so each cue block starts from a known state. ---------!
   subroutine reset_pheno_memory(poly)
      type(polygon_t), intent(inout) :: poly
      associate (n => poly%cohort%n)
         poly%cohort%leaf_flush_tendency(1:n) = 1.0_wp
         poly%cohort%leaf_shed_tendency(1:n)  = 0.0_wp
         poly%cohort%growing_degree_days(1:n) = 0.0_wp
         poly%cohort%cold_degree_days(1:n)    = 0.0_wp
         poly%cohort%dry_psi_sum(1:n)         = 0.0_wp
         poly%cohort%wet_psi_sum(1:n)         = 0.0_wp
         poly%cohort%light_hours_mean(1:n)    = -1.0_wp      ! no light memory
      end associate
   end subroutine reset_pheno_memory

   !----- Ithaca-like daily-mean air temperature [K]: ~270 in winter, ~297 in summer. --------!
   pure real(wp) function daily_tair(doy) result(t)
      integer(ik), intent(in) :: doy
      t = 283.15_wp + 14.0_wp * sin(TWOPI * (real(doy, wp) - 110.0_wp) / 365.0_wp)
   end function daily_tair

end program test_phenology_driver
