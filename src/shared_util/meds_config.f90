!==========================================================================================!
! meds_config -- immutable run configuration, threaded read-only through the engine.       !
!                                                                                          !
! Holds the time-stepping mode, the fusion/fission tunables (the diameter & size-           !
! distribution analogues of ED2's LAI/light tolerances), and the PFT trait table.          !
! `build_config` fills sensible defaults; `validate_config` error-stops on inconsistent     !
! settings (e.g. recruits born below the termination size, which would churn forever).      !
!==========================================================================================!
module meds_config
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : yr_day, size_tol
   use meds_pft_params, only : pft_table_t, init_default_pfts
   implicit none
   private

   public :: meds_config_t, build_config, validate_config
   public :: TS_DAILY, TS_WEEKLY, TS_MONTHLY, BK_SERIAL, BK_MULTICORE, BK_GPU

   !----- Time-step modes. ----------------------------------------------------------------!
   integer(ik), parameter :: TS_DAILY   = 1_ik
   integer(ik), parameter :: TS_WEEKLY  = 3_ik
   integer(ik), parameter :: TS_MONTHLY = 2_ik
   !----- Parallel backend labels (the actual backend is chosen at COMPILE time via the    !
   !      compiler's do-concurrent target; this is for reporting/reproducibility only).    !
   integer(ik), parameter :: BK_SERIAL    = 0_ik
   integer(ik), parameter :: BK_MULTICORE = 1_ik
   integer(ik), parameter :: BK_GPU       = 2_ik

   type :: meds_config_t
      !----- Time stepping. ---------------------------------------------------------------!
      integer(ik) :: ts_mode  = TS_DAILY
      real(wp)    :: dt_years = 1.0_wp / yr_day     !< set by build_config from ts_mode
      logical     :: vegetation_dynamics_on = .true. !< if .false. structure is frozen
      integer(ik) :: backend  = BK_SERIAL           !< reporting only

      !----- Fission/fusion master switches (passed by the stepper to update_demography). --!
      logical     :: do_cohort_fissfuse = .true.    !< run cohort fusion + split each month
      logical     :: do_patch_fissfuse  = .true.    !< run patch fusion/termination each year

      !----- Cohort fusion / termination. -------------------------------------------------!
      integer(ik) :: max_cohort       = 60_ik        !< >0 target, 0 disable, <0 force-merge
      integer(ik) :: n_cohort_fusion_iter    = 6_ik         !< tolerance-relaxation iterations
      real(wp)    :: cohort_size_tol_min = 0.02_wp     !< relative DBH/height tolerance, min
      real(wp)    :: cohort_size_tol_max = 0.10_wp     !< relative DBH/height tolerance, max
      real(wp)    :: cohort_size_tol_mult = 1.0_wp     !< geometric multiplier (derived)
      real(wp)    :: basal_area_bin_cap   = 2000.0_wp       !< [cm2/m2] single-cohort basal-area cap
                                                    !  (fusion will not merge beyond it; a cohort
                                                    !   above it splits). Set well above the per-
                                                    !   cohort mean so max_cohort governs the count.
      real(wp)    :: min_cohort_size = 1.0e-3_wp    !< [cm2/m2] cull below nplant*basal_area
      real(wp)    :: negligible_nplant = 1.0e-8_wp  !< [plant/m2] absolute density floor
      real(wp)    :: split_eps    = 1.0e-4_wp       !< symmetric DBH perturbation on split
      logical     :: enable_cohort_fission = .true.

      !----- Size-distribution profile bins (replace ED2 cumulative-LAI profile). ---------!
      integer(ik)           :: n_dbh_bins = 8_ik
      real(wp), allocatable :: dbh_edges(:)         !< ascending interior bin edges [cm]

      !----- Patch fusion / termination. --------------------------------------------------!
      integer(ik) :: max_patch       = 12_ik         !< >0 target, 0 disable, <0 force
      integer(ik) :: n_patch_fusion_iter   = 6_ik
      real(wp)    :: patch_profile_tol    = 0.20_wp      !< avg cumulative-BA profile distance
      real(wp)    :: patch_profile_maxdev_factor = 4.0_wp      !< max-deviation multiplier
      real(wp)    :: patch_diff_age_tol = 1.0_wp      !< [yr] same-age phase window
      real(wp)    :: min_patch_area  = 1.0e-4_wp    !< cull patches below this area fraction
      real(wp)    :: patch_min_area_remain = 0.99_wp  !< stop fusing once this area is kept
      logical     :: enable_patch_fission = .false. !< no clean ED2 analog; off by default

      !----- Recruitment. -----------------------------------------------------------------!
      real(wp) :: min_recruit_size = 1.0e-2_wp      !< [plant/m2] spawn threshold on the pool

      !----- Conservation check tolerance. ------------------------------------------------!
      real(wp) :: conservation_tol = size_tol               !< 1% basal-area / individuals tolerance

      !----- PFT traits. ------------------------------------------------------------------!
      type(pft_table_t) :: pft
   end type meds_config_t

contains

   !---------------------------------------------------------------------------------------!
   ! Build a default configuration for a given time-step mode.                             !
   !---------------------------------------------------------------------------------------!
   function build_config(ts_mode, backend) result(cfg)
      integer(ik), intent(in), optional :: ts_mode, backend
      type(meds_config_t)               :: cfg
      integer(ik)                       :: i

      if (present(ts_mode)) cfg%ts_mode = ts_mode
      if (present(backend)) cfg%backend = backend

      select case (cfg%ts_mode)
      case (TS_MONTHLY)
         cfg%dt_years = 1.0_wp / 12.0_wp
      case (TS_WEEKLY)
         cfg%dt_years = 7.0_wp / yr_day
      case default
         cfg%dt_years = 1.0_wp / yr_day
      end select

      !----- Geometric tolerance growth from min to max over niter iterations. ------------!
      if (cfg%n_cohort_fusion_iter > 1_ik) then
         cfg%cohort_size_tol_mult = (cfg%cohort_size_tol_max / cfg%cohort_size_tol_min)             &
                               ** (1.0_wp / real(cfg%n_cohort_fusion_iter - 1_ik, wp))
      else
         cfg%cohort_size_tol_mult = 1.0_wp
      end if

      !----- Default DBH bin edges: geometric-ish spacing up to ~150 cm. ------------------!
      allocate(cfg%dbh_edges(cfg%n_dbh_bins - 1_ik))
      do i = 1_ik, cfg%n_dbh_bins - 1_ik
         cfg%dbh_edges(i) = 2.0_wp * (2.0_wp ** real(i - 1_ik, wp))    ! 2,4,8,16,32,64,128
      end do

      call init_default_pfts(cfg%pft)
   end function build_config

   !---------------------------------------------------------------------------------------!
   ! Validate a configuration; halt on a setting that would corrupt the run.               !
   !---------------------------------------------------------------------------------------!
   subroutine validate_config(cfg)
      type(meds_config_t), intent(in) :: cfg
      character(len=*), parameter :: tag = 'meds_config: '

      if (cfg%pft%n < 1_ik)                          error stop tag//'empty PFT table'
      if (cfg%cohort_size_tol_min <= 0.0_wp)            error stop tag//'cohort_size_tol_min <= 0'
      if (cfg%cohort_size_tol_max < cfg%cohort_size_tol_min) error stop tag//'cohort_size_tol_max < min'
      if (cfg%n_cohort_fusion_iter < 1_ik)                   error stop tag//'n_cohort_fusion_iter < 1'
      if (cfg%n_patch_fusion_iter < 1_ik)                   error stop tag//'n_patch_fusion_iter < 1'
      if (cfg%n_dbh_bins < 2_ik)                        error stop tag//'n_dbh_bins < 2'
      if (cfg%min_patch_area <= 0.0_wp)              error stop tag//'min_patch_area <= 0'
      !----- A recruit must survive its own birth: pool threshold must exceed the cull. ---!
      if (cfg%min_recruit_size <= cfg%negligible_nplant)                                   &
         error stop tag//'min_recruit_size must exceed negligible_nplant'
   end subroutine validate_config

end module meds_config
