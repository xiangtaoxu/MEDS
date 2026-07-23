!==========================================================================================!
! meds_output_config -- the [output] configuration block for the diagnostic-aggregation      !
! subsystem, plus the selector codes shared by the config, the reader, the registry, the       !
! integrators and the serializer.                                                               !
!                                                                                          !
! Placed in src/shared (NOT src/io) so meds_config -- the DAG ROOT -- can carry output_config_t  !
! as a plain-scalar component with NO backward `shared -> io` edge (the SAME rule                 !
! meds_forcing_config uses for forcing_config_t). Pure DATA + parameters; links meds_kinds only.  !
! Design: docs/dev_plans/MEDS_IO_DESIGN.md sections 3.1, 6.4.                                              !
!                                                                                          !
! The temporal-tier (FREQ_*), variable-group (GRP_*), file-chunk (FC_*) and sync (SYNC_*) codes   !
! live HERE because output_config_t references them; the io-only reduction/axis codes (AGG_*/      !
! DIM_*) live with the io types in meds_output_types.                                              !
!==========================================================================================!
module meds_output_config
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: output_config_t
   public :: FREQ_FAST, FREQ_DAILY, FREQ_MONTHLY, FREQ_ANNUAL, FREQ_NONE, N_FREQ
   public :: GRP_STRUCTURE, GRP_CARBON, GRP_WATER, GRP_ENERGY, GRP_NUMERICS, N_GRP
   public :: FC_DAY, FC_MONTH, FC_YEAR, FC_RUN
   public :: SYNC_NEVER, SYNC_FLUSH
   public :: freq_tier_index, freq_letter

   !----- Temporal tiers, as bit positions of a per-variable membership mask (ior/iand). ------!
   !      FAST = a user-set multiple of dt_fast (§4.1); the tiers chain fast->daily->monthly->annual. !
   integer(ik), parameter :: FREQ_NONE    = 0_ik
   integer(ik), parameter :: FREQ_FAST    = 1_ik
   integer(ik), parameter :: FREQ_DAILY   = 2_ik
   integer(ik), parameter :: FREQ_MONTHLY = 4_ik
   integer(ik), parameter :: FREQ_ANNUAL  = 8_ik
   integer(ik), parameter :: N_FREQ       = 4_ik   !< number of temporal tiers (FAST/DAILY/MONTHLY/ANNUAL)

   !----- Variable GROUPS the main-config high-level flux toggles switch on/off (§6). ----------!
   integer(ik), parameter :: GRP_STRUCTURE = 1_ik  !< demography/size: nplant, dbh, agb, patch geometry, ids
   integer(ik), parameter :: GRP_CARBON    = 2_ik  !< carbon fluxes: gpp, nee (+ future npp/Rh)
   integer(ik), parameter :: GRP_WATER     = 3_ik  !< water fluxes: transpiration, soil moisture (+ future evap)
   integer(ik), parameter :: GRP_ENERGY    = 4_ik  !< energy fluxes/state: soil/leaf temperature (+ sensible/latent)
   !----- NUMERICS instrumentation (MEDS_NUMERICS_SCOPING.md section 5.3 "work" metrics): integrator  !
   !      WORK per step -- accepted/rejected substeps, sub-solver substeps, non-convergence events.    !
   !      Not a physical group; it is the cost axis of the goal-(b) benchmark, without which cost per  !
   !      unit accuracy (and therefore goal (c)'s selection rule) cannot be computed. Default OFF.     !
   integer(ik), parameter :: GRP_NUMERICS  = 5_ik
   integer(ik), parameter :: N_GRP         = 5_ik

   !----- File-chunk buckets: how many records per file for a stream (the time-chunk, §5.1). ----!
   integer(ik), parameter :: FC_DAY   = 1_ik  !< one file per sim-day
   integer(ik), parameter :: FC_MONTH = 2_ik  !< one file per sim-month
   integer(ik), parameter :: FC_YEAR  = 3_ik  !< one file per sim-year
   integer(ik), parameter :: FC_RUN   = 4_ik  !< ONE file for the whole run (annual default; §5.1)

   !----- nc_sync durability policy (§5.5): never, or after every flush. ------------------------!
   integer(ik), parameter :: SYNC_NEVER = 0_ik  !< no nc_sync (fastest; crash forfeits the open chunk)
   integer(ik), parameter :: SYNC_FLUSH = 1_ik  !< nc_sync after each flush (durable; default)

   !==========================================================================================!
   !  The [output] block. Plain scalars + small fixed arrays (no allocatables), so meds_config    !
   !  carries it trivially. The per-variable overrides are NOT stored here -- they live in the      !
   !  optional meds_io_config.toml, loaded by build_output_registry (§6.4). enabled defaults        !
   !  .false. so a config with no [output] block runs the legacy [io] path unchanged (§6.1).        !
   !==========================================================================================!
   type :: output_config_t
      logical            :: enabled     = .false.               !< master switch (replaces io.write_output)
      character(len=256) :: dir         = '.'                   !< output directory
      character(len=256) :: prefix      = 'meds'                !< filename stem
      character(len=256) :: io_config   = ''                    !< path to meds_io_config.toml ('' -> group/freq defaults only)
      integer(ik)        :: cohort_max  = 4096_ik               !< fixed netCDF cohort dim (live n must not exceed it)
      integer(ik)        :: patch_max   = 256_ik                !< fixed netCDF patch dim
      logical            :: strict_caps = .false.               !< .false. warn+truncate on n>cap; .true. error stop
      integer(ik)        :: sync_every  = SYNC_FLUSH            !< nc_sync policy: SYNC_NEVER | SYNC_FLUSH (§5.5)
      integer(ik)        :: fast_interval_steps = 4_ik          !< fast tier flushes every N*dt_fast (§4.1)
      !----- High-level flux-group toggles: switch whole GRP_* groups (main config, §6). ------!
      logical            :: grp_on(N_GRP) = [.true., .true., .false., .false., .false.] ! STRUCT/CARBON/WATER/ENERGY/NUMERICS
      !----- Per-tier enable + records-per-file bucket (index order FAST/DAILY/MONTHLY/ANNUAL).-!
      logical            :: freq_on(N_FREQ)    = [.false., .true., .true., .true.]
      integer(ik)        :: file_chunk(N_FREQ) = [FC_DAY, FC_MONTH, FC_YEAR, FC_RUN]
   end type output_config_t

contains

   !----- Tier index (1..N_FREQ) of a single FREQ_* bit; 0 if not a single tier bit. ----------!
   pure integer(ik) function freq_tier_index(freq_bit) result(idx)
      integer(ik), intent(in) :: freq_bit
      select case (freq_bit)
      case (FREQ_FAST)   ; idx = 1_ik
      case (FREQ_DAILY)  ; idx = 2_ik
      case (FREQ_MONTHLY); idx = 3_ik
      case (FREQ_ANNUAL) ; idx = 4_ik
      case default       ; idx = 0_ik
      end select
   end function freq_tier_index

   !----- Single-letter stream tag F/D/M/Y for a FREQ_* bit (file naming, §5.1). --------------!
   pure character(len=1) function freq_letter(freq_bit) result(c)
      integer(ik), intent(in) :: freq_bit
      select case (freq_bit)
      case (FREQ_FAST)   ; c = 'F'
      case (FREQ_DAILY)  ; c = 'D'
      case (FREQ_MONTHLY); c = 'M'
      case (FREQ_ANNUAL) ; c = 'Y'
      case default       ; c = '?'
      end select
   end function freq_letter

end module meds_output_config
