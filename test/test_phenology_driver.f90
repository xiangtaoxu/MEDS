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
!==========================================================================================!
program test_phenology_driver
   use meds_kinds,                only : wp, ik
   use meds_config,               only : meds_config_t
   use meds_core_interface,       only : site_t
   use meds_init,                 only : init_bare_ground, add_cohort
   use meds_vegetation_dynamics,  only : advance_leaf_phenology
   use meds_plant_interface,      only : CUE_TEMP, CUE_NONE
   use meds_test_support,         only : build_test_config
   implicit none

   real(wp), parameter :: twopi = 6.283185307179586_wp
   type(meds_config_t) :: cfg
   type(site_t)        :: site
   integer(ik) :: nfail, doy
   real(wp)    :: fl_temp_200, sh_temp_200, fl_temp_340, sh_temp_340
   real(wp)    :: fl_ever_200, sh_ever_200, fl_ever_340, sh_ever_340, gdd_summer

   nfail = 0_ik

   !----- Config: one temperature-deciduous PFT (1), the rest evergreen; phenology is             !
   !       unconditional now (docs/dev_plans/MEDS_SLOW_DYNAMICS_DESIGN.md Part I) -- this test     !
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
      site%pheno_tair_sum = daily_tair(doy)
      site%pheno_tair_n   = 1_ik
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
      site%pheno_tair_sum = 0.0_wp ; site%pheno_tair_n = 0_ik      ! no fast sub-steps ran
      call advance_leaf_phenology(site, cfg, 1_ik)
      call check_true('no-temperature step leaves flush_drive unchanged', &
                      abs(site%cohort%pheno_flush_drive(1) - fl_before) < tiny(1.0_wp))
      call check_true('no-temperature step leaves shed_drive unchanged', &
                      abs(site%cohort%pheno_shed_drive(1) - sh_before) < tiny(1.0_wp))
      call check_true('no-temperature step leaves GDD unchanged', &
                      abs(site%cohort%pheno_gdd(1) - gdd_before) < tiny(1.0_wp))
   end block

   if (nfail == 0_ik) then
      print '(a)', 'test_phenology_driver: ALL PASSED'
   else
      print '(a,i0,a)', 'test_phenology_driver: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   !----- Ithaca-like daily-mean air temperature [K]: ~270 in winter, ~297 in summer. --------!
   pure real(wp) function daily_tair(doy) result(t)
      integer(ik), intent(in) :: doy
      t = 283.15_wp + 14.0_wp * sin(twopi * (real(doy, wp) - 110.0_wp) / 365.0_wp)
   end function daily_tair

   subroutine check_true(name, cond)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      if (cond) then
         print '(a,a)', '  ok   : ', name
      else
         nfail = nfail + 1_ik
         print '(a,a)', '  FAIL : ', name
      end if
   end subroutine check_true

   subroutine check_int(name, got, expect)
      character(len=*), intent(in) :: name
      integer(ik),      intent(in) :: got, expect
      if (got == expect) then
         print '(a,a)', '  ok   : ', name
      else
         nfail = nfail + 1_ik
         print '(a,a,i0,a,i0)', '  FAIL : ', name, got, ' expected ', expect
      end if
   end subroutine check_int

   subroutine check_close(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) <= atol) then
         print '(a,a)', '  ok   : ', name
      else
         nfail = nfail + 1_ik
         print '(a,a,es13.6,a,es13.6)', '  FAIL : ', name, got, ' expected ', expect
      end if
   end subroutine check_close

end program test_phenology_driver
