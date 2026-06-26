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
   use meds_allometry,  only : height_max
   use meds_pft_params, only : pft_table_t, init_default_pfts
   implicit none
   private

   public :: meds_config_t, build_config, validate_config
   public :: TS_DAILY, TS_WEEKLY, TS_MONTHLY, BK_SERIAL, BK_MULTICORE, BK_GPU
   public :: DIST_PRIMARY, DIST_TREEFALL

   !----- Time-step modes. ----------------------------------------------------------------!
   integer(ik), parameter :: TS_DAILY   = 1_ik
   integer(ik), parameter :: TS_WEEKLY  = 3_ik
   integer(ik), parameter :: TS_MONTHLY = 2_ik
   !----- Parallel backend labels (the actual backend is chosen at COMPILE time via the    !
   !      compiler's do-concurrent target; this is for reporting/reproducibility only).    !
   integer(ik), parameter :: BK_SERIAL    = 0_ik
   integer(ik), parameter :: BK_MULTICORE = 1_ik
   integer(ik), parameter :: BK_GPU       = 2_ik
   !----- Patch disturbance / land-use classes. -------------------------------------------!
   integer(ik), parameter :: DIST_PRIMARY  = 1_ik   !< undisturbed / primary stand
   integer(ik), parameter :: DIST_TREEFALL = 2_ik   !< treefall-gap (age-0) patch

   type :: meds_config_t
      !----- Time stepping. ---------------------------------------------------------------!
      integer(ik) :: ts_mode  = TS_DAILY
      real(wp)    :: dt_years = 1.0_wp / yr_day     !< set by build_config from ts_mode
      logical     :: vegetation_dynamics_on = .true. !< if .false. structure is frozen
      integer(ik) :: backend  = BK_SERIAL           !< reporting only

      !----- Structural master switches (passed by the stepper to update_demography). ------!
      logical     :: do_cohort_fissfuse  = .true.   !< run cohort fusion + split each month
      logical     :: do_patch_fissfuse   = .true.   !< run patch fusion/termination each year
      logical     :: do_patch_disturbance = .true.  !< open treefall gaps each year

      !----- Cohort fusion / termination. -------------------------------------------------!
      integer(ik) :: max_cohort       = 60_ik        !< >0 target, 0 disable, <0 force-merge
      integer(ik) :: n_cohort_fusion_iter    = 6_ik         !< tolerance-relaxation iterations
      real(wp)    :: cohort_size_tol_min = 0.02_wp     !< relative DBH/height tolerance, min
      real(wp)    :: cohort_size_tol_max = 0.10_wp     !< relative DBH/height tolerance, max
      real(wp)    :: cohort_size_tol_mult = 1.0_wp     !< geometric multiplier (derived)
      real(wp)    :: cohort_lai_cap = 5.0_wp       !< [m2/m2] single-cohort LAI cap (fusion will
                                                    !  not merge beyond it; a cohort above it
                                                    !  splits). Set well above the per-cohort mean
                                                    !  so max_cohort governs the working count.
      real(wp)    :: min_cohort_agb = 1.0e-6_wp    !< [kgC/m2] cull below nplant*agb
      real(wp)    :: negligible_nplant = 1.0e-8_wp  !< [plant/m2] absolute density floor
      real(wp)    :: split_eps    = 1.0e-4_wp       !< symmetric DBH perturbation on split
      logical     :: enable_cohort_fission = .true.

      !----- Vertical light profile for patch fusion (cumulative-LAI by height layer, ED2). !
      integer(ik)           :: n_height_layers = 16_ik
      real(wp), allocatable :: height_edges(:)      !< ascending interior layer edges [m]

      !----- Patch fusion / termination. --------------------------------------------------!
      integer(ik) :: max_patch       = 12_ik         !< >0 target, 0 disable, <0 force
      integer(ik) :: n_patch_fusion_iter   = 6_ik
      real(wp)    :: patch_light_tol    = 0.10_wp      !< avg light-profile distance to fuse
      real(wp)    :: patch_light_maxdev_factor = 1.5_wp  !< max single-layer deviation multiplier
      real(wp)    :: patch_diff_age_tol = 1.0_wp      !< [yr] same-age phase window
      real(wp)    :: min_patch_area  = 1.0e-4_wp    !< cull patches below this area fraction
      real(wp)    :: patch_min_area_remain = 0.99_wp  !< stop fusing once this area is kept
      logical     :: enable_patch_fission = .false. !< no clean ED2 analog; off by default

      !----- Patch disturbance (ED2 treefall): a new age-0 gap patch per structural step. --!
      real(wp) :: patch_disturbance_rate     = 0.014_wp  !< [1/yr] fraction of area disturbed
      real(wp) :: disturbance_survive_height = 10.0_wp   !< [m] tall cohorts (>=) die in the gap;
                                                         !  short understory survives (ED2 treefall)

      !----- Recruitment. -----------------------------------------------------------------!
      real(wp) :: min_recruit_size = 1.0e-2_wp      !< [plant/m2] spawn threshold on the pool

      !----- Conservation check tolerance. ------------------------------------------------!
      real(wp) :: conservation_tol = size_tol               !< 1% AGB / individuals tolerance

      !----- Initial conditions. ----------------------------------------------------------!
      character(len=256) :: init_census_file = ''   !< non-empty => start from this cohort census CSV
                                                     !  (else near-bare ground); see meds_init.

      !----- netCDF output (written by meds_main; effective only with MEDS_ENABLE_IO). ----!
      logical            :: io_write_output  = .true.        !< write netCDF state output at all
      character(len=256) :: io_output_dir    = '.'           !< directory for output files (created if missing)
      character(len=256) :: io_output_prefix = 'meds_output' !< output file is <io_output_dir>/<io_output_prefix>.nc
      integer(ik) :: io_output_interval_years = 1_ik  !< append a full snapshot every N years
      integer(ik) :: io_cohort_max = 2048_ik          !< fixed netCDF cohort dimension (cap+slack)
      integer(ik) :: io_patch_max  = 64_ik            !< fixed netCDF patch dimension (cap+slack)

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

      !----- Evenly spaced height-layer edges from 0 to the canopy-height cap. ------------!
      allocate(cfg%height_edges(cfg%n_height_layers - 1_ik))
      do i = 1_ik, cfg%n_height_layers - 1_ik
         cfg%height_edges(i) = real(i, wp) * height_max / real(cfg%n_height_layers, wp)
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
      if (cfg%n_height_layers < 2_ik)                error stop tag//'n_height_layers < 2'
      if (cfg%min_patch_area <= 0.0_wp)              error stop tag//'min_patch_area <= 0'
      if (cfg%cohort_lai_cap <= 0.0_wp)              error stop tag//'cohort_lai_cap <= 0'
      if (cfg%patch_disturbance_rate < 0.0_wp)       error stop tag//'patch_disturbance_rate < 0'
      if (cfg%disturbance_survive_height <= 0.0_wp)  error stop tag//'disturbance_survive_height <= 0'
      if (any(cfg%pft%wood_density <= 0.0_wp))       error stop tag//'wood_density <= 0'
      !----- A recruit must survive its own birth: pool threshold must exceed the cull. ---!
      if (cfg%min_recruit_size <= cfg%negligible_nplant)                                   &
         error stop tag//'min_recruit_size must exceed negligible_nplant'
   end subroutine validate_config

end module meds_config
