!==========================================================================================!
! test_phenology_driver -- integration test for the slow-loop leaf-phenology WIRING           !
! (meds_phenology_driver.leaf_phenology), as opposed to the stateless kernel (test_plant_        !
! phenology). It drives a two-cohort site through one synthetic Ithaca year, feeding a daily     !
! air temperature via the site accumulator the fast loop normally fills, and checks:             !
!                                                                                          !
!   1. TEMP-DECIDUOUS : a CUE_TEMP cohort flushes (PHEN_ON) in mid-summer and drops (PHEN_OFF)    !
!                       in late autumn; its GDD memory accumulates during the growing season.     !
!   2. EVERGREEN      : a CUE_NONE cohort keeps its born-leafed PHEN_ON status all year.           !
!   3. NO-TEMPERATURE : a step with no accumulated air temperature (pheno_tair_n = 0) is skipped   !
!                       (the memory is not advanced on a bogus 0/0 mean).                          !
!==========================================================================================!
program test_phenology_driver
   use meds_kinds,                only : wp, ik
   use meds_config,               only : meds_config_t
   use meds_demography_interface, only : site_t
   use meds_init,                 only : init_bare_ground, add_cohort
   use meds_phenology_driver,     only : leaf_phenology
   use meds_plant_interface,      only : CUE_TEMP, CUE_NONE, PHEN_ON, PHEN_OFF
   use meds_test_support,         only : build_test_config
   implicit none

   real(wp), parameter :: twopi = 6.283185307179586_wp
   type(meds_config_t) :: cfg
   type(site_t)        :: site
   integer(ik) :: nfail, doy, st_temp_200, st_temp_340, st_ever_200, st_ever_340
   real(wp)    :: gdd_summer

   nfail = 0_ik

   !----- Config: one temperature-deciduous PFT (1), the rest evergreen; phenology ON. --------!
   cfg = build_test_config()
   cfg%phenology_on            = .true.
   cfg%forcing%latitude_deg    = 42.44_wp            ! Ithaca NY (northern hemisphere)
   cfg%pft%pheno_cue_mask(1)   = CUE_TEMP            ! PFT 1: cold-deciduous
   cfg%pft%pheno_cue_mask(2:)  = CUE_NONE            ! PFT 2,3: evergreen (unchanged default)

   !----- A site with two cohorts: cohort 1 = PFT 1 (deciduous), cohort 2 = PFT 2 (evergreen). !
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.5_wp, 10.0_wp)
   call add_cohort(site, cfg, 1_ik, 2_ik, 0.5_wp, 10.0_wp)
   call check_int('two cohorts created', int(site%cohort%n, ik), 2_ik)
   call check_int('cohort 1 starts leafed (PHEN_ON)', site%cohort%phenology_status(1), PHEN_ON)

   !----- Drive one synthetic year (dt_slow = 1 day). Each day set the site daily-mean air temp !
   !      (the fast loop's accumulator) to an Ithaca-like sinusoid, then advance the phenology.  !
   st_temp_200 = -99_ik ; st_temp_340 = -99_ik ; st_ever_200 = -99_ik ; st_ever_340 = -99_ik
   gdd_summer = 0.0_wp
   do doy = 1_ik, 365_ik
      site%pheno_tair_sum = daily_tair(doy)
      site%pheno_tair_n   = 1_ik
      call leaf_phenology(site, cfg, doy)
      if (doy == 200_ik) then
         st_temp_200 = site%cohort%phenology_status(1)
         st_ever_200 = site%cohort%phenology_status(2)
         gdd_summer  = site%cohort%pheno_gdd(1)
      end if
      if (doy == 340_ik) then
         st_temp_340 = site%cohort%phenology_status(1)
         st_ever_340 = site%cohort%phenology_status(2)
      end if
   end do

   !----- 1. Temperature-deciduous seasonal cycle. -----------------------------------------!
   call check_int('deciduous cohort ON in mid-summer (doy 200)', st_temp_200, PHEN_ON)
   call check_int('deciduous cohort OFF in late autumn (doy 340)', st_temp_340, PHEN_OFF)
   call check_true('deciduous GDD accumulated by summer', gdd_summer > 100.0_wp)

   !----- 2. Evergreen cohort never leaves PHEN_ON. ----------------------------------------!
   call check_int('evergreen cohort ON in summer', st_ever_200, PHEN_ON)
   call check_int('evergreen cohort ON in autumn', st_ever_340, PHEN_ON)

   !----- 3. A no-temperature step is skipped (status + memory unchanged). ------------------!
   block
      integer(ik) :: st_before
      real(wp)    :: gdd_before
      st_before  = site%cohort%phenology_status(1)
      gdd_before = site%cohort%pheno_gdd(1)
      site%pheno_tair_sum = 0.0_wp ; site%pheno_tair_n = 0_ik      ! no fast sub-steps ran
      call leaf_phenology(site, cfg, 1_ik)
      call check_int('no-temperature step leaves status unchanged', &
                     site%cohort%phenology_status(1), st_before)
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

end program test_phenology_driver
