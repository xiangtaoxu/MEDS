! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_output_types -- the pure DATA descriptors of the diagnostic-aggregation subsystem:      !
! the per-variable registry descriptor (var_desc_t), the running per-(variable,tier) reduction   !
! buffer (integ_buffer_t), and the registry container (output_registry_t).                        !
!                                                                                          !
! netCDF-FREE by construction (links meds_kinds + meds_output_config ONLY, no meds_netcdf_c):      !
! so the stepper edge that references the integrate kernels over these types pulls no C           !
! dependency (the DAG-hygiene wall, §2). The io-only reduction/axis codes AGG_*/DIM_* live here;    !
! the config-shared FREQ_*/GRP_*/FC_* live in src/config/meds_output_config. Because the C          !
! bindings are NOT visible here, xtype and the fill/missing sentinels are declared LOCALLY          !
! (XTYPE_*, MISSING_*) and mapped to NC_* only in the serializer. Design: MEDS_IO_DESIGN.md §3.     !
!==========================================================================================!
module meds_output_types
   use meds_kinds,         only : wp, ik
   use meds_time,          only : meds_time_t
   use meds_output_config,      only : N_FREQ
   use meds_column_params, only : n_soil_layer_max
   implicit none
   private

   public :: var_desc_t, integ_buffer_t, output_registry_t, diag_params_t
   public :: pending_record_t, record_queue_t, stream_file_t
   public :: output_files_t, output_buffers_t
   public :: AGG_MEAN, AGG_SUM, AGG_MIN, AGG_MAX, AGG_LAST, AGG_TMEAN, AGG_FLUXSUM
   public :: DIM_SCALAR, DIM_COHORT, DIM_PATCH, DIM_SOIL, DIM_PFT, DIM_SIZE, DIM_SOIL_PATCH
   public :: XTYPE_DOUBLE, XTYPE_INT
   public :: MISSING_VALUE, MISSING_INT, MAX_OUTPUT_VARS, MAX_DBH_CLASS
   public :: agg_is_slabwise, ragged_dim

   !----- Temporal reduction operators (the `agg` of a variable). -----------------------------!
   integer(ik), parameter :: AGG_MEAN    = 1_ik  !< equal-weight arithmetic mean
   integer(ik), parameter :: AGG_SUM     = 2_ik  !< plain sum (integer count tallies ONLY)
   integer(ik), parameter :: AGG_MIN     = 3_ik  !< period minimum
   integer(ik), parameter :: AGG_MAX     = 4_ik  !< period maximum
   integer(ik), parameter :: AGG_LAST    = 5_ik  !< end-of-period snapshot (ids / CSR / instantaneous)
   integer(ik), parameter :: AGG_TMEAN   = 6_ik  !< dt-weighted state mean (the physical-stock default)
   integer(ik), parameter :: AGG_FLUXSUM = 7_ik  !< dt-weighted integral of a rate (period total)

   !----- Trailing (per-record) axis of a variable. DIM_COHORT/DIM_PATCH are fixed WITHIN a       !
   !      window (§4.4), so all non-scalar dims fold by direct slot index -- no id-keying.         !
   integer(ik), parameter :: DIM_SCALAR = 0_ik
   integer(ik), parameter :: DIM_COHORT = 1_ik
   integer(ik), parameter :: DIM_PATCH  = 2_ik
   integer(ik), parameter :: DIM_SOIL   = 3_ik
   integer(ik), parameter :: DIM_PFT    = 4_ik
   integer(ik), parameter :: DIM_SIZE   = 5_ik   !< DBH size class (ED2-style; edges from TOML)
   !----- 2-D (soil layer x patch), FLATTENED into the same 1-D slab as every other axis at        !
   !      index = (ip-1)*n_soil_layer_max + k. The stride is the COMPILE-TIME layer ceiling, never   !
   !      a live count -- striding by a live count would silently re-map every profile the moment     !
   !      the active layer count differed from what a reader assumed, producing plausible shifted     !
   !      columns rather than an error. Because it is a plain 1-D slab here, integrate_slab /         !
   !      normalize_slab / reset_buffer serve it with NO new code; only the serializer's dim mapping   !
   !      differs (it writes rank-3 (time, patch, soil)). MEDS_IO_V01_PLAN.md section 3.6.1. ---------!
   integer(ik), parameter :: DIM_SOIL_PATCH = 6_ik

   !----- On-disk numeric type. Kept LOCAL (not NC_*) so this module stays netCDF-free; the        !
   !      serializer maps XTYPE_DOUBLE->NC_DOUBLE, XTYPE_INT->NC_INT.                               !
   integer(ik), parameter :: XTYPE_DOUBLE = 1_ik
   integer(ik), parameter :: XTYPE_INT    = 2_ik

   !----- Fill / missing sentinels. Values MATCH netCDF NC_FILL_DOUBLE / NC_FILL_INT so a reader     !
   !      recognizes them, but are declared here to keep the module netCDF-free (§4.3, §5.3).        !
   real(wp),    parameter :: MISSING_VALUE = 9.9692099683868690e+36_wp
   integer(ik), parameter :: MISSING_INT   = -2147483647_ik

   !----- Registry capacity. Raised 128 -> 512 for the v0.1 variable set (~230 registered, plus     !
   !      headroom); the list is still short enough that a linear find_var_index is fine.            !
   integer(ik), parameter :: MAX_OUTPUT_VARS = 512_ik
   !----- DBH size-class ceiling (DIM_SIZE). Edges come from [output].dbh_class_edges. ---------!
   integer(ik), parameter :: MAX_DBH_CLASS   = 32_ik

   !==========================================================================================!
   ! The run-dependent parameters the DERIVED diagnostics need, resolved once at manager setup   !
   ! and carried as plain data (netCDF-free, allocatable-free).                                   !
   !                                                                                          !
   ! WHY THIS EXISTS. extract_variable is handed `site` and a descriptor -- it has no config      !
   ! aggregator and no fast context. But a soil matric potential needs the retention curve's       !
   ! texture parameters, a PFT-axis reduction needs the run's PFT count, and a size-class          !
   ! reduction needs the class edges. Passing the whole meds_config_t down would drag the config    !
   ! aggregator into the per-step tick; passing this small bundle keeps the tick's dependencies     !
   ! exactly as narrow as they were.                                                                !
   !                                                                                          !
   ! The soil block is COPIED FROM THE FAST CONTEXT the physics actually ran (meds_main wires it),  !
   ! not re-derived from the TOML. That is deliberate: a psi reported here and a psi the roots saw  !
   ! are then the same curve by construction, rather than two independent derivations that agree     !
   ! until someone edits one of them.                                                                !
   !==========================================================================================!
   type :: diag_params_t
      integer(ik) :: n_pft       = 0_ik      !< run-time PFT count (the DIM_PFT axis length)
      integer(ik) :: n_soil      = 0_ik      !< active soil layers
      !----- Soil retention curve: the SOIL_RETENTION_* family code plus the van Genuchten /     !
      !      Campbell parameters, PER LAYER (soil_params_t stores them per layer, and texture can   !
      !      legitimately vary with depth -- collapsing them to a scalar here would silently        !
      !      report the wrong curve for every layer but one).  --------------------------------!
      integer(ik) :: retention   = 1_ik
      real(wp)    :: theta_sat(n_soil_layer_max) = 0.0_wp
      real(wp)    :: theta_res(n_soil_layer_max) = 0.0_wp
      real(wp)    :: par_a(n_soil_layer_max)     = 0.0_wp
      real(wp)    :: par_n(n_soil_layer_max)     = 0.0_wp
      !----- Layer NODE depths [m], negative downward -- the `soil` axis coordinate. Written to any  !
      !      file carrying the soil dim, for the same reason the pft and dbh_class axes carry theirs:  !
      !      the vertical grid is a run-time property (depth, layer count, geometric growth), so a     !
      !      file without it cannot be plotted against depth without re-deriving the grid by hand.     !
      real(wp)    :: soil_z(n_soil_layer_max)    = 0.0_wp
      logical     :: soil_ready  = .false.   !< .false. => psi/wetness diagnostics emit _FillValue
      !----- DBH size classes (DIM_SIZE): n_class+1 ascending edges [cm]. ----------------------!
      integer(ik) :: n_dbh_class = 0_ik
      real(wp)    :: dbh_edges(MAX_DBH_CLASS + 1_ik) = 0.0_wp
   end type diag_params_t

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
      integer(ik) :: source_id       = 0_ik         !< which SoA FIELD supplies it (§3.3)
      !----- HOW the field is aggregated to this variable's scale. These three are what let ONE     !
      !      per-cohort field emit its cohort / patch / site / PFT / size-class variants: the        !
      !      reduction is data on the descriptor, not a hand-written routine per variable.           !
      !                                                                                          !
      !      `weight` (W_* in meds_diagnostic_reduce) and `mean` together encode the EXTENSIVE vs     !
      !      INTENSIVE distinction, which is a physical statement and the classic place diagnostics   !
      !      go wrong: agb is extensive (weighted SUM over nplant -> per ground area), whereas dbh     !
      !      is intensive (basal-area-weighted MEAN) and leaf_temp is intensive (leaf-area-weighted    !
      !      MEAN). Declaring it here means it is stated once, at registration, beside the units.      !
      integer(ik) :: weight          = 0_ik         !< W_* per-cohort weight kind (0 = W_NONE)
      logical     :: mean            = .false.      !< .true. weighted MEAN (intensive); .false. SUM
      real(wp)    :: scale           = 1.0_wp       !< unit conversion applied inside the patch loop
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
      real(wp)    :: wsum  = 0.0_wp         !< Sum(dt) weight (AGG_TMEAN/FLUXSUM)
      integer(ik) :: nsamp = 0_ik           !< equal-weight sample count (MEAN/SUM/MIN/MAX/LAST)
      real(wp)    :: seed  = 0.0_wp         !< re-seed value (MIN=+huge, MAX=-huge, else 0)
      !----- per-index accumulators (DIM_COHORT/PATCH/SOIL/PFT), sized to the relevant cap. ------!
      real(wp),    allocatable :: slab(:)
      real(wp),    allocatable :: wsum_slab(:)
      integer(ik), allocatable :: hits(:)
      integer(ik) :: n_slab = 0_ik          !< live index count this window (the [1:n] cohort/patch count)
   end type integ_buffer_t

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
   ! by the stepper-side tick (output_integrate) into the tier's scratch record, then copied into   !
   ! the tier's queue; drained by main (output_serialize) in the I/O phase at month                  !
   ! boundaries (MEDS_POLYGON_RUNTIME_PLAN.md §4, B11). Payload is indexed by REGISTRY var index;   !
   ! only variables live in the tier are filled (§4.5, §2).                                         !
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
      !----- A queued record stores only its tier's slab variables: slab(:, col(k)) is registry     !
      !      variable k's slab. The scratch record leaves col unallocated and indexes slab(:, k). ---!
      integer(ik), allocatable :: col(:)        !< (nvar) slab column of each variable (0 = none)
   end type pending_record_t

   public :: slab_col

   !----- The closed records of one tier waiting for the I/O phase, in closing order. A queued     !
   !      record keeps only the slab rows the writer reads, not the scratch's max_slab rows, so a   !
   !      month of fast-tier records costs the live cohort count, not cohort_max. Elements are       !
   !      reused across months; `n` counts the live ones. --------------------------------------------!
   type :: record_queue_t
      integer(ik) :: n = 0_ik
      type(pending_record_t), allocatable :: rec(:)
   end type record_queue_t

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
      logical     :: has_pft      = .false. !< this tier defines the pft dim
      logical     :: has_size     = .false. !< this tier defines the dbh-class dim
      integer(ik) :: d_time = -1_ik, d_cohort = -1_ik, d_patch = -1_ik, d_soil = -1_ik
      integer(ik) :: d_pft = -1_ik, d_size = -1_ik
      integer(ik) :: d_polygon = -1_ik       !< a region file's polygon dimension
      integer(ik) :: v_pft = -1_ik, v_dbh_lower = -1_ik, v_dbh_upper = -1_ik  !< self-describing axis coords
      integer(ik) :: v_soil_z = -1_ik                                        !< soil layer node depths [m]
      integer(ik) :: cohort_dim = 0_ik, patch_dim = 0_ik   !< the file's ACTUAL trimmed cohort/patch axis length
      integer(ik) :: v_time = -1_ik, v_year = -1_ik, v_month = -1_ik, v_day = -1_ik
      integer(ik) :: v_hour = -1_ik, v_minute = -1_ik, v_second = -1_ik   !< FAST-tier sub-daily companions
      integer(ik) :: v_ncohort = -1_ik, v_npatch = -1_ik
      integer(ik), allocatable :: vid(:)    !< (nvar) netCDF varid per registry var (-1 if not in tier)
   end type stream_file_t

   !==========================================================================================!
   ! The output manager comes in two types (MEDS_POLYGON_RUNTIME_PLAN.md §10.3, R2). Both are        !
   ! netCDF-FREE plain data; only output_serialize touches C.                                        !
   !   output_files_t   -- one set of output files (a site's, a region's, a detail polygon's) and    !
   !                       what writing them needs: registry, diagnostic parameters, file settings,  !
   !                       stream handles. Read-only while a step runs.                              !
   !   output_buffers_t -- one polygon's buffers for one set of files: its running reductions, open  !
   !                       windows, closed records and fast-tier staging. The only output state a    !
   !                       step writes. Not a cache: the reductions exist nowhere else.              !
   ! A site run has one of each. main owns them; the stepper ticks the buffers.                      !
   !==========================================================================================!
   type :: output_files_t
      logical                 :: enabled = .false.
      !----- The run has a fast loop. Without one the fast-loop diagnostic blocks stay off, so their  !
      !      variables read as missing rather than as a 0 nothing computed (#299). ------------------!
      logical                 :: fast_loop_on = .true.
      type(output_registry_t) :: reg
      type(diag_params_t)     :: diag       !< run-dependent params the derived diagnostics need
      integer(ik)       :: cohort_max = 0_ik, patch_max = 0_ik, max_slab = 0_ik
      logical           :: cohort_axis = .false., patch_axis = .false.  !< a live variable has that axis
      type(stream_file_t)    :: stream(N_FREQ)
      character(len=256)     :: dir = '.', prefix = 'meds'
      !----- Forcing provenance written as a global attribute on every output file. With the ED_ERA5land  !
      !      archive the model computes qair itself, so a run's humidity depends on the model version and  !
      !      the formula is recorded (MEDS_FORCING_DESIGN.md §15.4, FD8); empty -> no attribute.          !
      character(len=512)     :: forcing_qair = ''
      integer(ik)           :: file_chunk(N_FREQ) = 0_ik
      integer(ik)           :: sync_every = 1_ik
      integer(ik)           :: fast_interval_steps = 4_ik   !< fast tier closes every N*dt_fast sub-steps
      !----- A REGION's polygon axis (MEDS_POLYGON_RUNTIME_PLAN.md §6): its files carry a `polygon`  !
      !      dimension, one entry per polygon's buffers, with these coordinates (0 for site files). --!
      integer(ik)              :: n_polygon = 0_ik
      integer(ik), allocatable :: polygon_id(:)             !< the cell's row-major index on the forcing grid
      integer(ik), allocatable :: polygon_row(:), polygon_col(:)   !< 0-based grid indices
      real(wp),    allocatable :: polygon_lat(:), polygon_lon(:)   !< [deg] cell centre
   end type output_files_t

   type :: output_buffers_t
      type(integ_buffer_t), allocatable :: buf(:,:)   !< (nvar, N_FREQ) running reductions
      logical           :: has_data(N_FREQ) = .false. !< tier's current window has >=1 sample
      type(meds_time_t) :: t_open(N_FREQ)             !< period-start of each tier's current window
      type(pending_record_t) :: pending(N_FREQ)          !< per-tier scratch that close_tier normalizes into
      type(record_queue_t)   :: queue(N_FREQ)            !< closed records awaiting the I/O phase
      !----- FAST (sub-daily) tier staging (netCDF-free): filled per (patch, sub-step) by the fast   !
      !      loop, replayed into buf(:,1) by main via output_integrate_fast. Each patch's sub-step    !
      !      is the patch block's own row (PD_*) and its soil column; the site means are their        !
      !      area-weighted sums, and the forcing is the polygon block's row (PY_*). The per-cohort   !
      !      slabs are indexed by the site's cohort slot. The buffers follow the live patch and       !
      !      cohort counts; the output layer checks those against output.patch_max and cohort_max.   !
      logical              :: fast_on = .false.                !< output on and the FAST tier has live variables
      type(meds_time_t),   allocatable :: fast_time(:)         !< (n_fast_sub) each sub-step's start
      real(wp),            allocatable :: fast_forcing(:,:)     !< (PY_*, n_fast_sub) the sub-step's forcing
      real(wp),            allocatable :: fast_patch(:,:,:)     !< (PD_*, patch, n_fast_sub) each patch's row
      real(wp),            allocatable :: fast_site(:,:)        !< (PD_*, n_fast_sub) their area-weighted sum
      !----- The soil columns, per patch in the soil-by-patch slab's layout, (ip-1)*n_soil_layer_max  !
      !      + k (DIM_SOIL_PATCH), and area-weighted over the patches. ------------------------------!
      real(wp),            allocatable :: fast_soil_temp_patch(:,:)   !< (layer x patch, n_fast_sub) [K]
      real(wp),            allocatable :: fast_soil_water_patch(:,:)  !< (layer x patch, n_fast_sub) [m3/m3]
      real(wp),            allocatable :: fast_soil_temp(:,:)   !< (n_soil, n_fast_sub)  area-weighted [K]
      real(wp),            allocatable :: fast_soil_water(:,:)  !< (n_soil, n_fast_sub)  area-weighted [m3/m3]
      real(wp),            allocatable :: fast_coh_ltemp(:,:)   !< (cohort, n_fast_sub) per-cohort leaf temp [K]
      real(wp),            allocatable :: fast_coh_gpp(:,:)     !< (cohort, n_fast_sub) per-cohort GPP [umol/plant/s]
      !< (cohort, n_fast_sub) per-cohort height [m] (tallest post-proc)
      real(wp),            allocatable :: fast_coh_height(:,:)
      integer(ik)          :: n_fast_sub   = 0_ik              !< sub-steps staged this slow step
      integer(ik)          :: fast_n_soil  = 0_ik              !< live soil layers in the fast slabs
      integer(ik)          :: fast_n_cohort = 0_ik             !< live site cohorts in the fast cohort slabs
      integer(ik)          :: fast_n_patch  = 0_ik             !< live patches in the fast patch rows
      logical              :: fast_ready   = .false.           !< the staging is filled and awaits replay
   end type output_buffers_t

contains

   !----- The slab column of registry variable k in a record (identity for the scratch record). --!
   pure integer(ik) function slab_col(pr, k) result(c)
      type(pending_record_t), intent(in) :: pr
      integer(ik),            intent(in) :: k
      if (allocated(pr%col)) then ; c = pr%col(k) ; else ; c = k ; end if
   end function slab_col

   !----- .true. if the operator averages/weights (needs seed/wsum handling), used by callers. --!
   pure logical function agg_is_slabwise(dim) result(yes)
      integer(ik), intent(in) :: dim
      yes = (dim /= DIM_SCALAR)
   end function agg_is_slabwise

   !----- The axes whose length and slots follow the stand: cohort, patch, and soil by patch. They  !
   !      are fixed only within a month (§4.4), so a file of them spans at most a month, the annual  !
   !      stream refuses them, and a region file, whose polygons differ, cannot hold them. --------!
   pure logical function ragged_dim(dim) result(yes)
      integer(ik), intent(in) :: dim
      yes = dim == DIM_COHORT .or. dim == DIM_PATCH .or. dim == DIM_SOIL_PATCH
   end function ragged_dim

end module meds_output_types
