!==========================================================================================!
! meds_output_types -- the pure DATA descriptors of the diagnostic-aggregation subsystem:      !
! the per-variable registry descriptor (var_desc_t), the running per-(variable,tier) reduction   !
! buffer (integ_buffer_t), and the registry container (output_registry_t).                        !
!                                                                                          !
! netCDF-FREE by construction (links meds_kinds + meds_output_config ONLY, no meds_netcdf_c):      !
! so the stepper edge that references the integrate kernels over these types pulls no C           !
! dependency (the DAG-hygiene wall, §2). The io-only reduction/axis codes AGG_*/DIM_* live here;    !
! the config-shared FREQ_*/GRP_*/FC_* live in src/shared/meds_output_config. Because the C          !
! bindings are NOT visible here, xtype and the fill/missing sentinels are declared LOCALLY          !
! (XTYPE_*, MISSING_*) and mapped to NC_* only in the serializer. Design: MEDS_IO_DESIGN.md §3.     !
!==========================================================================================!
module meds_output_types
   use meds_kinds,         only : wp, ik
   use meds_time,          only : meds_time_t
   use meds_output_config, only : N_FREQ
   implicit none
   private

   public :: var_desc_t, integ_buffer_t, output_registry_t
   public :: pending_record_t, stream_file_t, output_manager_t, fast_sample_t
   public :: AGG_MEAN, AGG_SUM, AGG_MIN, AGG_MAX, AGG_LAST, AGG_MEANSQ, AGG_TMEAN, AGG_FLUXSUM
   public :: DIM_SCALAR, DIM_COHORT, DIM_PATCH, DIM_SOIL, DIM_PFT
   public :: XTYPE_DOUBLE, XTYPE_INT
   public :: MISSING_VALUE, MISSING_INT, MAX_OUTPUT_VARS
   public :: agg_is_slabwise

   !----- Temporal reduction operators (the `agg` of a variable). -----------------------------!
   integer(ik), parameter :: AGG_MEAN    = 1_ik  !< equal-weight arithmetic mean
   integer(ik), parameter :: AGG_SUM     = 2_ik  !< plain sum (integer count tallies ONLY)
   integer(ik), parameter :: AGG_MIN     = 3_ik  !< period minimum
   integer(ik), parameter :: AGG_MAX     = 4_ik  !< period maximum
   integer(ik), parameter :: AGG_LAST    = 5_ik  !< end-of-period snapshot (ids / CSR / instantaneous)
   integer(ik), parameter :: AGG_MEANSQ  = 6_ik  !< second moment -> variance companion (DEFERRED, §3.5)
   integer(ik), parameter :: AGG_TMEAN   = 7_ik  !< dt-weighted state mean (the physical-stock default)
   integer(ik), parameter :: AGG_FLUXSUM = 8_ik  !< dt-weighted integral of a rate (period total)

   !----- Trailing (per-record) axis of a variable. DIM_COHORT/DIM_PATCH are fixed WITHIN a       !
   !      window (§4.4), so all non-scalar dims fold by direct slot index -- no id-keying.         !
   integer(ik), parameter :: DIM_SCALAR = 0_ik
   integer(ik), parameter :: DIM_COHORT = 1_ik
   integer(ik), parameter :: DIM_PATCH  = 2_ik
   integer(ik), parameter :: DIM_SOIL   = 3_ik
   integer(ik), parameter :: DIM_PFT    = 4_ik

   !----- On-disk numeric type. Kept LOCAL (not NC_*) so this module stays netCDF-free; the        !
   !      serializer maps XTYPE_DOUBLE->NC_DOUBLE, XTYPE_INT->NC_INT.                               !
   integer(ik), parameter :: XTYPE_DOUBLE = 1_ik
   integer(ik), parameter :: XTYPE_INT    = 2_ik

   !----- Fill / missing sentinels. Values MATCH netCDF NC_FILL_DOUBLE / NC_FILL_INT so a reader     !
   !      recognizes them, but are declared here to keep the module netCDF-free (§4.3, §5.3).        !
   real(wp),    parameter :: MISSING_VALUE = 9.9692099683868690e+36_wp
   integer(ik), parameter :: MISSING_INT   = -2147483647_ik

   !----- Registry capacity (the source registration list is short; §3.4). --------------------!
   integer(ik), parameter :: MAX_OUTPUT_VARS = 128_ik

   !==========================================================================================!
   ! One diagnostic variable: a pure DATA descriptor, NO pointers into the SoA (§3.1, §3.3).      !
   !==========================================================================================!
   type :: var_desc_t
      character(len=32) :: name      = ''          !< netCDF variable name (scale-suffixed, e.g. 'agb_cohort')
      character(len=96) :: long_name = ''          !< CF long_name attribute
      character(len=24) :: units     = ''          !< CF units attribute
      integer(ik) :: dim             = DIM_SCALAR   !< DIM_* : the trailing axis
      integer(ik) :: agg             = AGG_MEAN     !< AGG_* : temporal reduction operator
      integer(ik) :: group           = 1_ik         !< GRP_* : the flux group a high-level toggle switches (§6)
      integer(ik) :: streams         = 0_ik         !< ior(FREQ_*) EFFECTIVE membership (after §6.1 resolution)
      integer(ik) :: streams_default = 0_ik         !< the registry default membership (restored by a `true` override)
      logical     :: enabled         = .true.       !< master per-variable on/off
      integer(ik) :: xtype           = XTYPE_DOUBLE !< XTYPE_DOUBLE | XTYPE_INT
      integer(ik) :: source_id       = 0_ik         !< which SoA field / reduction supplies it (§3.3)
   end type var_desc_t

   !==========================================================================================!
   ! One running reduction for one (variable, tier) pair (§4.2). Single addressing mode: every    !
   ! non-scalar dim folds by direct slot index (cohort/patch fixed within a window, §4.4).         !
   !==========================================================================================!
   type :: integ_buffer_t
      integer(ik) :: var_id = 0_ik          !< index into registry%var(:)
      integer(ik) :: freq   = 0_ik          !< the single FREQ_* bit this buffer serves
      integer(ik) :: agg    = AGG_MEAN
      integer(ik) :: dim    = DIM_SCALAR
      logical     :: active = .false.       !< allocated + participating (this var is on in this tier)
      !----- scalar accumulators (DIM_SCALAR). -------------------------------------------------!
      real(wp)    :: scal  = 0.0_wp         !< running reduction
      real(wp)    :: scal2 = 0.0_wp         !< second moment (AGG_MEANSQ)
      real(wp)    :: wsum  = 0.0_wp         !< Sum(dt) weight (AGG_TMEAN/FLUXSUM/MEANSQ)
      integer(ik) :: nsamp = 0_ik           !< equal-weight sample count (MEAN/SUM/MIN/MAX/LAST)
      real(wp)    :: seed  = 0.0_wp         !< re-seed value (MIN=+huge, MAX=-huge, else 0)
      !----- per-index accumulators (DIM_COHORT/PATCH/SOIL/PFT), sized to the relevant cap. ------!
      real(wp),    allocatable :: slab(:)
      real(wp),    allocatable :: slab2(:)
      real(wp),    allocatable :: wsum_slab(:)
      integer(ik), allocatable :: hits(:)
      integer(ik) :: n_slab = 0_ik          !< live index count this window (the [1:n] cohort/patch count)
   end type integ_buffer_t

   !==========================================================================================!
   ! One per-(patch,sub-step) LIVE sample of the site-scalar fast-loop diagnostics, assembled     !
   ! inside run_fast_biophysics and area-weighted-accumulated into the manager's fast(:) staging   !
   ! (§FAST tier). Pure scalars -- the DIM_SOIL / DIM_COHORT fast slabs live as 2-D arrays on the   !
   ! manager, not here, so this type stays allocatable-free and trivially default-constructs to 0.  !
   !==========================================================================================!
   type :: fast_sample_t
      real(wp) :: cas_temp      = 0.0_wp   !< [K]         area-weighted CAS temperature
      real(wp) :: soil_temp_top = 0.0_wp   !< [K]         area-weighted top-soil temperature
      real(wp) :: gpp_rate      = 0.0_wp   !< [umol/m2/s] instantaneous site GPP rate
      real(wp) :: le_flux       = 0.0_wp   !< [W/m2]      latent-heat (ET) flux
      real(wp) :: h_flux        = 0.0_wp   !< [W/m2]      sensible-heat flux
      real(wp) :: rnet          = 0.0_wp   !< [W/m2]      net all-wave radiation absorbed by the column
      real(wp) :: sw_in         = 0.0_wp   !< [W/m2]      incident shortwave at canopy top
      real(wp) :: ustar         = 0.0_wp   !< [m/s]       friction velocity
   end type fast_sample_t

   !==========================================================================================!
   ! The immutable registration list + the precomputed per-tier live-variable index (§3.2).       !
   !==========================================================================================!
   type :: output_registry_t
      type(var_desc_t), allocatable :: var(:)      !< the registration list [1:nvar]
      integer(ik)                   :: nvar = 0_ik
      integer(ik), allocatable      :: idx_freq(:,:) !< (MAX_OUTPUT_VARS, N_FREQ): var indices live in each tier
      integer(ik)                   :: nidx(N_FREQ) = 0_ik
   end type output_registry_t

   !==========================================================================================!
   ! A closed period staged for the serializer (netCDF-FREE plain data). Filled at a roll-over    !
   ! by the stepper-side tick (output_integrate); drained by main (output_serialize_pending). One  !
   ! per tier can be pending at a time (main drains every step). Payload is indexed by REGISTRY     !
   ! var index; only variables live in the tier are filled (§4.5, §2).                              !
   !==========================================================================================!
   type :: pending_record_t
      logical           :: used     = .false.
      integer(ik)       :: freq     = 0_ik      !< the FREQ_* bit of the tier this record closes
      type(meds_time_t) :: t_open              !< period-start stamp (calendar companions, §5.3)
      integer(ik)       :: n_cohort = 0_ik, n_patch = 0_ik  !< live slab lengths this record
      real(wp),    allocatable :: sval(:)       !< (nvar) normalized scalar value
      logical,     allocatable :: svalid(:)     !< (nvar)
      real(wp),    allocatable :: slab(:,:)     !< (max_slab, nvar) normalized slab values
      logical,     allocatable :: slabvalid(:,:)!< (max_slab, nvar)
      integer(ik), allocatable :: nslab(:)      !< (nvar) slab length (0 for scalar vars)
   end type pending_record_t

   !==========================================================================================!
   ! One open netCDF stream file (the serializer's per-tier handle). ncids/varids kept as plain    !
   ! integers so this type stays netCDF-agnostic; the serializer casts to c_int (§7.3).            !
   !==========================================================================================!
   type :: stream_file_t
      integer(ik) :: freq         = 0_ik
      integer(ik) :: ncid         = -1_ik   !< open handle (-1 = none)
      integer(ik) :: nrec         = 0_ik    !< records written to the current file
      integer(ik) :: chunk_bucket = -1_ik   !< current time-chunk bucket key (-1 = no file open)
      logical     :: has_cohort   = .false. !< this tier defines the cohort dim
      logical     :: has_patch    = .false. !< this tier defines the patch dim
      logical     :: has_soil     = .false. !< this tier defines the soil dim
      integer(ik) :: d_time = -1_ik, d_cohort = -1_ik, d_patch = -1_ik, d_soil = -1_ik
      integer(ik) :: cohort_dim = 0_ik, patch_dim = 0_ik   !< the file's ACTUAL trimmed cohort/patch axis length
      integer(ik) :: v_time = -1_ik, v_year = -1_ik, v_month = -1_ik, v_day = -1_ik
      integer(ik) :: v_hour = -1_ik, v_minute = -1_ik, v_second = -1_ik   !< FAST-tier sub-daily companions
      integer(ik) :: v_ncohort = -1_ik, v_npatch = -1_ik
      integer(ik), allocatable :: vid(:)    !< (nvar) netCDF varid per registry var (-1 if not in tier)
   end type stream_file_t

   !==========================================================================================!
   ! The output manager: netCDF-FREE plain-data glue (registry + integrators + pending stage +     !
   ! stream handles). main owns it; the stepper ticks it; only output_serialize_pending touches C.  !
   !==========================================================================================!
   type :: output_manager_t
      logical                 :: enabled = .false.
      type(output_registry_t) :: reg
      type(integ_buffer_t), allocatable :: buf(:,:)   !< (nvar, N_FREQ) running reductions
      logical           :: has_data(N_FREQ) = .false. !< tier's current window has >=1 sample
      type(meds_time_t) :: t_open(N_FREQ)             !< period-start of each tier's current window
      integer(ik)       :: cohort_max = 0_ik, patch_max = 0_ik, max_slab = 0_ik
      type(pending_record_t) :: pending(N_FREQ)
      type(stream_file_t)    :: stream(N_FREQ)
      character(len=256)     :: dir = '.', prefix = 'meds'
      integer(ik)           :: file_chunk(N_FREQ) = 0_ik
      integer(ik)           :: sync_every = 1_ik
      integer(ik)           :: fast_interval_steps = 4_ik   !< fast tier closes every N*dt_fast sub-steps
      !----- FAST (sub-daily) tier staging (netCDF-free): filled per (patch,sub-step) by the fast     !
      !      loop, replayed into buf(:,1) by main via output_integrate_fast. Site scalars in fast(:);  !
      !      the area-weighted soil column + per-cohort slabs in 2-D [slot, sub-step] arrays.  --------!
      type(fast_sample_t), allocatable :: fast(:)              !< (n_fast_sub) site-scalar samples
      type(meds_time_t),   allocatable :: fast_time(:)         !< (n_fast_sub) sub-step midpoint stamps
      real(wp),            allocatable :: fast_soil_temp(:,:)   !< (n_soil, n_fast_sub)  area-weighted [K]
      real(wp),            allocatable :: fast_soil_water(:,:)  !< (n_soil, n_fast_sub)  area-weighted [m3/m3]
      real(wp),            allocatable :: fast_coh_ltemp(:,:)   !< (cohort_max, n_fast_sub) per-cohort leaf temp [K]
      real(wp),            allocatable :: fast_coh_gpp(:,:)     !< (cohort_max, n_fast_sub) per-cohort GPP [umol/plant/s]
      integer(ik)          :: n_fast_sub   = 0_ik              !< sub-steps staged this slow step
      integer(ik)          :: fast_n_soil  = 0_ik              !< live soil layers in the fast slabs
      integer(ik)          :: fast_n_cohort = 0_ik             !< live site cohorts in the fast cohort slabs
      logical              :: fast_ready   = .false.           !< fast(:) filled + awaiting replay
   end type output_manager_t

contains

   !----- .true. if the operator averages/weights (needs seed/wsum handling), used by callers. --!
   pure logical function agg_is_slabwise(dim) result(yes)
      integer(ik), intent(in) :: dim
      yes = (dim /= DIM_SCALAR)
   end function agg_is_slabwise

end module meds_output_types
