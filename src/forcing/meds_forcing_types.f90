!==========================================================================================!
! meds_forcing_types -- the RUNTIME meteorological-forcing types (design MEDS_FORCING_DESIGN   !
! .md sections 3.1-3.2): the instantaneous per-site record met_forcing_t (a read-only boundary   !
! condition, the analogue of rad_forcing_t / chydro_forcing_t), one raw file record met_record_t, !
! and the MUTABLE per-polygon reader buffer met_driver_t (the ED2 cgrid%metinput analogue).        !
!                                                                                          !
! Lives in src/forcing (libmeds_forcing); links meds_shared only (meds_kinds, meds_time for the     !
! meds_time_t timestamp, meds_forcing_config for the [forcing]/[site] config). No physics library    !
! ever names these types; the reader (meds_met_driver) owns the only mutable forcing state.           !
!==========================================================================================!
module meds_forcing_types
   use meds_kinds,          only : wp, ik
   use meds_time,           only : meds_time_t
   use meds_forcing_config, only : forcing_config_t, MET_BACKEND_CONST
   implicit none
   private

   public :: met_forcing_t, met_record_t, met_driver_t

   !==========================================================================================!
   !  met_forcing_t -- the instantaneous per-SITE atmospheric state the fast loop consumes.       !
   !  Per-site scalars + four SW-band reals (no cohort-sized allocatables) -> trivially copyable,   !
   !  GPU-mappable. All defaults are a valid "reference climate" box (the four SW streams SUM to      !
   !  400 W/m2 = today's fast_context_t%rad_sw_top), so met_forcing_t() is usable as the CONST         !
   !  backend value. cosz and rho_air are DERIVED each substep (not read).                             !
   !==========================================================================================!
   type :: met_forcing_t
      real(wp) :: tair_k       = 288.0_wp     !< [K]        air temperature at reference height
      real(wp) :: qair         = 0.008_wp     !< [kg/kg]    specific humidity
      real(wp) :: psurf_pa     = 101325.0_wp  !< [Pa]       surface pressure
      real(wp) :: rainf        = 0.0_wp       !< [kg/m2/s]  liquid precipitation rate (post phase-split)
      real(wp) :: snowf        = 0.0_wp       !< [kg/m2/s]  frozen precip
      real(wp) :: wind         = 2.0_wp       !< [m/s]      wind speed at reference height
      real(wp) :: lwdown       = 380.0_wp     !< [W/m2]     downwelling longwave (positive down)
      real(wp) :: par_beam     = 180.0_wp     !< [W/m2]     direct-beam PAR at canopy top
      real(wp) :: par_diffuse  = 40.0_wp      !< [W/m2]     diffuse PAR
      real(wp) :: nir_beam     = 150.0_wp     !< [W/m2]     direct-beam NIR
      real(wp) :: nir_diffuse  = 30.0_wp      !< [W/m2]     diffuse NIR   (Sigma = 400 W/m2)
      real(wp) :: co2          = 420.0_wp     !< [umol/mol] free-atmosphere CO2
      real(wp) :: cosz         = 0.0_wp       !< [-]        cos(solar zenith); DERIVED each substep
      real(wp) :: rho_air      = 1.2_wp       !< [kg/m3]    DERIVED from tair/psurf/qair
   contains
      procedure :: swdown         => met_swdown          !< total downward SW = sum of the four streams
      procedure :: rshort_diffuse => met_rshort_diffuse  !< diffuse SW = par_diffuse + nir_diffuse
   end type met_forcing_t

   !----- One raw timestamped record as read (pre-interpolation). The four SW streams are      !
   !      already split at INGEST (partition_shortwave) from the file's total SWdown.           !
   type :: met_record_t
      type(meds_time_t) :: when
      real(wp) :: tair_k = 288.0_wp, qair = 0.008_wp, psurf_pa = 101325.0_wp
      real(wp) :: rainf = 0.0_wp, wind = 2.0_wp, lwdown = 380.0_wp, co2 = 420.0_wp
      real(wp) :: par_beam = 180.0_wp, par_diffuse = 40.0_wp, nir_beam = 150.0_wp, nir_diffuse = 30.0_wp
   end type met_record_t

   !==========================================================================================!
   !  met_driver_t -- the MUTABLE per-POLYGON reader state (the ONLY mutable forcing state).      !
   !  Holds the two records bracketing the model time at this polygon's grid_index, plus the       !
   !  cursor and the cached time axis. For a multi-polygon run there is one per polygon, each        !
   !  bound to its grid_index (all reading the same file, different `grid` slices).                  !
   !==========================================================================================!
   type :: met_driver_t
      type(forcing_config_t)  :: fcfg                       !< [forcing]/[site] config
      integer(ik) :: backend    = MET_BACKEND_CONST         !< NETCDF | CONST
      integer(ik) :: ncid       = -1_ik                     !< NetCDF handle (via meds_netcdf_c), -1 if const
      integer(ik) :: grid_index = 1_ik                      !< which location (1..ngrid) this polygon reads
      integer(ik) :: ngrid      = 1_ik                      !< total locations in the file
      integer(ik) :: nrec       = 0_ik                      !< total time records in the file
      integer(ik) :: irec_prev  = 0_ik                      !< cursor: index of the "previous" bracketing record
      real(wp)    :: dt_forcing  = 3600.0_wp                !< [s] native forcing interval
      type(met_record_t) :: rec_prev, rec_next              !< the two records bracketing the model time
      type(meds_time_t)  :: base_time                       !< time-coordinate anchor ("seconds since base_time")
      real(wp), allocatable :: time_sec(:)                  !< [s] cached time coordinate (seconds since base_time)
      !----- multi-year calendar recycling (derived at open; see derive_cycle_years). --------!
      integer(ik) :: file_year1      = 0_ik                 !< calendar year of record #1 (Jan-1-aligned whole-year file)
      integer(ik) :: n_cycle_years   = 0_ik                 !< number of whole calendar years the file spans
      logical     :: recycle_calendar = .false.             !< recycle in the CALENDAR domain (else legacy span-wrap)
      logical     :: at_wrap_seam    = .false.              !< current bracket is the cycle-boundary seam (rec nrec -> rec 1)
   end type met_driver_t

contains

   pure real(wp) function met_swdown(met) result(sw)
      class(met_forcing_t), intent(in) :: met
      sw = met%par_beam + met%par_diffuse + met%nir_beam + met%nir_diffuse
   end function met_swdown

   pure real(wp) function met_rshort_diffuse(met) result(diff)
      class(met_forcing_t), intent(in) :: met
      diff = met%par_diffuse + met%nir_diffuse
   end function met_rshort_diffuse

end module meds_forcing_types
