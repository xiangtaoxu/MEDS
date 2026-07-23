!==========================================================================================!
! meds_core_state_types -- the demographic state, as a flat site_t-wide Structure-of-Arrays.            !
!                                                                                          !
! ALL cohorts of the WHOLE site_t live in one contiguous set of 1-D arrays (`cohort_block`),  !
! so the dominant daily kernels are a single unit-stride sweep, ideal for SIMD and GPU.     !
! Patch membership is a CSR map (`cohort_offset`,`cohort_count`) over the flat cohort arrays. Cohorts of   !
! a patch occupy a contiguous slice, so per-patch operations (sort, fusion) work on slices.  !
!                                                                                          !
! Every structural change goes through ONE centralized routine, `cohort_reorder`, which     !
! permutes EVERY per-cohort array in lockstep -- this is the single place to update when a   !
! field is added, eliminating the ED2 "forgot to reallocate an array" class of bug.         !
!==========================================================================================!
module meds_core_state_types
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : pio4, tiny_num
   use meds_pft_params, only : pft_table_t
   use meds_allometry,  only : dbh_to_height, dbh_to_agb, dbh_to_leaf_area, wood_to_dbh, carbon_to_structure
   use meds_column_state_types, only : cas_state_t, soil_column_t, soil_energy_column_t,        &
                                       snow_column_t, soil_carbon_t, xi_accum_t, N_HYDRO_NODE, LEAF_TEMP_INIT, PSI_INIT
   implicit none
   private

   public :: cohort_block, patch_block, site_t
   public :: site_alloc, site_free
   public :: cohort_ensure_capacity, cohort_reorder, cohort_compact, gather_pft_params
   public :: patch_ensure_capacity, rebuild_csr, copy_cohort_slot, set_cohort_size, init_cohort
   public :: set_cohort_size_from_carbon, carbon_flux_block
   public :: cohort_deriv_block, cohort_deriv_alloc
   public :: assign_cohort_id, assign_patch_id
   public :: GROWTH_AVG_UNSET, PHENO_FLUSH_INIT, PHENO_SHED_INIT

   !----- Sentinel for a not-yet-sampled moving-average growth: a freshly created cohort holds  !
   !      this until its first growth step adds a sample (negative => "unset"; a real growth     !
   !      average is >= 0). Used by the growth kernel and the rate evaluator.                    !
   real(wp), parameter :: GROWTH_AVG_UNSET = -1.0_wp

   !----- Initial phenology GOVERNOR drives of a freshly created cohort: born leafed / flushing,  !
   !      no active shed (= the evergreen fixed point). Plain literals because `state` may not     !
   !      depend on the plant library (the DAG wall); the phenology driver advances them once real !
   !      drivers arrive. flush_drive=1 with no active shed keeps a fresh canopy building.          !
   real(wp), parameter :: PHENO_FLUSH_INIT = 1.0_wp
   real(wp), parameter :: PHENO_SHED_INIT  = 0.0_wp

   !---------------------------------------------------------------------------------------!
   ! All cohorts of the site_t (contiguous SoA). `dbh` is the prognostic size axis; `height`,   !
   ! `basal_area`, `agb` and `leaf_area` are derived but stored so the hot kernels stay         !
   ! arithmetic-only. `agb` (per-plant carbon) is the quantity conserved across fusion/        !
   ! fission; `leaf_area` (per-plant) gives the cohort's LAI contribution (nplant*leaf_area).   !
   !---------------------------------------------------------------------------------------!
   type :: cohort_block
      integer(ik) :: n   = 0_ik       !< cohorts in use
      integer(ik) :: cap = 0_ik       !< allocated capacity
      !----- Prognostic / primary state. --------------------------------------------------!
      integer(ik), allocatable :: pft(:)            !< PFT index
      real(wp),    allocatable :: nplant(:)         !< [plant/m2] PRIMARY demographic state
      real(wp),    allocatable :: dbh(:)            !< [cm]       prognostic size
      real(wp),    allocatable :: height(:)           !< [m]        = dbh_to_height(dbh), cached
      real(wp),    allocatable :: basal_area(:)        !< [cm2/plant]= pio4*dbh^2, cached
      real(wp),    allocatable :: agb(:)             !< [kgC/plant] conserved carbon, cached
      real(wp),    allocatable :: leaf_area(:)           !< [m2/plant]  leaf area, cached (LAI=nplant*leaf_area)
      real(wp),    allocatable :: overtopping_lai(:)     !< [m2/m2] cumulative LAI of all TALLER cohorts in the patch
                                                         !<         (Beer competition context); a RECOMPUTED diagnostic
                                                         !<         filled by meds_competition%update_overtopping_lai
      !----- Carbon pools [kgC/plant], on-allometry cached diagnostics (derived in set_cohort_size   !
      !      from agb & leaf_area). PR3 introduces the fields; PR4 promotes wood_carbon to the        !
      !      prognostic size anchor. leaf_carbon*sla = leaf_area; wood_carbon*aboveground_frac = agb.  !
      real(wp),    allocatable :: leaf_carbon(:)          !< [kgC/plant] = leaf_area / sla
      real(wp),    allocatable :: fineroot_carbon(:)      !< [kgC/plant] = root_to_leaf_ratio*leaf_carbon
      real(wp),    allocatable :: wood_carbon(:)          !< [kgC/plant] = agb / aboveground_frac
      real(wp),    allocatable :: nonstructural_carbon(:) !< [kgC/plant] = storage_cushion*leaf_carbon
      !----- Sliding simple-moving-average growth rate [cm/yr] driving growth-dependent          !
      !      mortality. growth_hist is the per-cohort ring buffer of the last growth_window steps !
      !      (column per cohort; the write slot is the site-global growth_hist_pos), growth_accum !
      !      its running sum and growth_count its fill level (<= window); growth_avg =            !
      !      accum/count is the SMA (GROWTH_AVG_UNSET until the first sample). A PROGNOSTIC        !
      !      history (not derivable from dbh).                                                    !
      real(wp),    allocatable :: growth_avg(:)
      real(wp),    allocatable :: growth_accum(:)   !< [cm/yr] running sum of the samples in the window
      integer(ik), allocatable :: growth_count(:)   !< samples currently in the window (<= growth_window)
      real(wp),    allocatable :: growth_hist(:,:)  !< (growth_window, cohort) ring buffer of growth samples
      !----- Per-cohort gathered PFT params (kernels never index the PFT table). ----------!
      real(wp),    allocatable :: p_dbh_critical(:)
      real(wp),    allocatable :: p_wood_density(:)
      real(wp),    allocatable :: p_hgt_max(:)        !< [m] gathered per-PFT asymptotic max height
      !----- Gathered carbon-allometry traits (for the on-allometry pool diagnostics above). --!
      real(wp),    allocatable :: p_aboveground_frac(:)   !< [--] aboveground fraction of woody carbon
      real(wp),    allocatable :: p_root_to_leaf_ratio(:) !< [--] fine-root:leaf target ratio
      real(wp),    allocatable :: p_storage_cushion(:)    !< [--] storage target as multiple of leaf target
      !----- DYNAMIC leaf traits (no p_ prefix => mutable): seeded from the PFT top-of-canopy values !
      !       at birth and acclimated to light by meds_plant_trait_dynamics; leaf-area-weighted on    !
      !       cohort fusion. sla enters the leaf carbon<->area map, llspan sets baseline leaf turnover,!
      !       vcmax25/rd25 feed leaf gas exchange. -------------------------------------------------!
      real(wp),    allocatable :: sla(:)            !< [m2/kgC]     specific leaf area
      real(wp),    allocatable :: vcmax25(:)        !< [umol/m2/s]  max carboxylation at 25 degC
      real(wp),    allocatable :: rd25(:)           !< [umol/m2/s]  leaf dark respiration at 25 degC
      real(wp),    allocatable :: llspan(:)         !< [yr]         leaf lifespan (baseline turnover = 1/llspan)
      !----- PROGNOSTIC fast-biophysics per-cohort state (owned here so it rides the cohort     !
      !      lockstep; mutated by the fast loop, leaf-area-weighted on cohort fusion). psi carries !
      !      genuine sub-slow-step hydraulic memory; leaf_temp is a warm-start for the VPD lag.   !
      real(wp),    allocatable :: leaf_temp(:)      !< [K]   cohort leaf temperature
      real(wp),    allocatable :: wood_temp(:)      !< [K]   cohort wood/branch temperature (own energy store)
      real(wp),    allocatable :: psi(:,:)          !< [MPa] (N_HYDRO_NODE, cohort) node water potentials
      real(wp),    allocatable :: gpp_accum(:)      !< [kgC/plant] GROSS GPP accumulated over the slow step
                                                    !<            (fast->slow carbon bridge; reset each slow step)
      !----- Autotrophic MAINTENANCE respiration accumulated over the slow step (fast->slow bridge, !
      !      reset each slow step; kept SEPARATE from gross gpp_accum, mirroring ED2 today_*_resp).  !
      real(wp),    allocatable :: leaf_resp_accum(:) !< [kgC/plant] leaf dark respiration  (ED2 today_leaf_resp)
      real(wp),    allocatable :: stem_resp_accum(:) !< [kgC/plant] stem maintenance resp  (ED2 today_stem_resp)
      real(wp),    allocatable :: root_resp_accum(:) !< [kgC/plant] fine-root maint. resp  (ED2 today_root_resp)
      !----- PROGNOSTIC leaf-phenology GOVERNOR drives + thermal cue memory (owned here so they ride !
      !      the cohort lockstep). Advanced daily by the slow-loop phenology advance from the cohort's !
      !      PFT cue params + drivers; the two rates (k*drive) are derived at the carbon consumer.     !
      !      Survivor-keeps on cohort fusion (like growth_avg -- the donor's memory is discarded).     !
      real(wp),    allocatable :: pheno_flush_drive(:)  !< [-] smoothed flush governor in [0,1] (0 => dormant)
      real(wp),    allocatable :: pheno_shed_drive(:)   !< [-] smoothed active-shed governor in [0,1] (0 => none)
      real(wp),    allocatable :: pheno_gdd(:)          !< [K day] growing-degree-day sum   (CUE_TEMP)
      real(wp),    allocatable :: pheno_chill(:)        !< [day]   chilling-day count        (CUE_TEMP)
      !----- Host-only back-index used to regroup the flat array by patch. ----------------!
      integer(ik), allocatable :: owner_patch(:)
      !----- Persistent identity: a global id stamped at creation and carried (in lockstep   !
      !      with every other field) through sorts/fusion/compaction, so an external tracker  !
      !      can follow one cohort across output records until it fuses away or is culled.     !
      integer(ik), allocatable :: global_id(:)
   end type cohort_block

   !---------------------------------------------------------------------------------------!
   ! Patches of the site_t (SoA) with a CSR map into the cohort block.                        !
   !---------------------------------------------------------------------------------------!
   type :: patch_block
      integer(ik) :: n   = 0_ik
      integer(ik) :: cap = 0_ik
      real(wp),    allocatable :: area(:)           !< fraction, sum = 1 (kept real64)
      real(wp),    allocatable :: age(:)            !< [yr]
      integer(ik), allocatable :: dist_type(:)      !< disturbance/land-use class
      integer(ik), allocatable :: cohort_offset(:)           !< CSR: first cohort of patch (1-based)
      integer(ik), allocatable :: cohort_count(:)         !< cohorts in patch
      real(wp),    allocatable :: recruit_pool(:,:)  !< (pft, patch) carry-forward density
      !----- Persistent identity (as for cohorts): stamped at creation, carried through       !
      !      patch sort/fusion/compaction so a tracker can follow one patch across records.    !
      integer(ik), allocatable :: global_id(:)
      !----- PROGNOSTIC fast-biophysics reservoirs (owned here so they ride the patch lockstep  !
      !      alongside area/age; mutated by the sub-daily fast loop, area-weighted on fusion).   !
      type(cas_state_t),          allocatable :: cas(:)      !< canopy-air-space twins per patch
      type(soil_energy_column_t), allocatable :: soil_e(:)   !< soil thermal column per patch
      type(soil_column_t),        allocatable :: soil_w(:)   !< soil water column per patch
      type(snow_column_t),        allocatable :: snow(:)     !< temporary-surface-water / snow store per patch
      !----- PROGNOSTIC slow biogeochemistry reservoir (owned here so it rides the patch lockstep    !
      !      alongside the fast reservoirs above; written DAILY by the slow soil-carbon step,          !
      !      area-weighted on fusion; docs/dev_plans/MEDS_SLOW_DYNAMICS_DESIGN.md Part II §8.1). ------!
      type(soil_carbon_t),        allocatable :: soil_carbon(:) !< soil-carbon pools per patch
      !----- Daily fast->slow accumulator for the soil-carbon matrix (B2): reset each slow step,   !
      !      accumulated once per (patch, fast sub-step) by column_prepass, consumed by the daily     !
      !      soil_carbon_step. Rides the SAME lockstep as soil_carbon (area-weighted on fusion). ------!
      type(xi_accum_t),           allocatable :: xi_accum(:)   !< [day]/[kgC/m2] per patch
   end type patch_block

   !----- TRANSIENT per-cohort TIME-DERIVATIVE bundle (all rates [unit/yr]). The slow-loop driver !
   !      COMPUTES these (carbon or empirical) each step; update_cohort_states then advances the   !
   !      cohort state by them. It is index-aligned with cohort_block but is NOT persistent state:  !
   !      filled the same step it is consumed (before any sort/fuse) and NOT lockstep-reordered.    !
   !      dln_nplant_dt is the LOG-space mortality rate (nplant *= exp(dln*dt)); every other field  !
   !      is an additive rate (state += rate*dt).                                                   !
   type :: cohort_deriv_block
      real(wp), allocatable :: d_dbh_dt(:), d_height_dt(:), d_basal_area_dt(:)
      real(wp), allocatable :: d_agb_dt(:), d_leaf_area_dt(:)
      real(wp), allocatable :: d_leaf_carbon_dt(:), d_fineroot_carbon_dt(:)
      real(wp), allocatable :: d_wood_carbon_dt(:), d_nonstructural_carbon_dt(:)
      real(wp), allocatable :: dln_nplant_dt(:)
   end type cohort_deriv_block

   type :: site_t
      type(cohort_block) :: cohort
      type(patch_block)  :: patch
      !----- TRANSIENT slow-loop scratch: the per-cohort tendency bundle the driver fills each     !
      !      step and update_cohort_states applies. Site-carried so it allocates once (resize-on-   !
      !      grow) instead of per step; refilled every step, so it is NOT lockstep-reordered and    !
      !      NOT serialized (a scratch buffer, not prognostic state).                               !
      type(cohort_deriv_block) :: deriv
      real(wp)           :: site_area = 1.0_wp
      integer(ik)        :: n_pft     = 0_ik
      !----- Monotonic counters for the persistent global ids (never reused within a run). ---!
      integer(ik)        :: next_cohort_id = 1_ik
      integer(ik)        :: next_patch_id  = 1_ik
      !----- Global ring-buffer write slot for the moving-average growth (same for all cohorts !
      !      since they all take one sample per step); advanced by update_demography.           !
      integer(ik)        :: growth_hist_pos = 0_ik
      !----- Site daily-mean air-temperature accumulator for phenology (the deferred "daily      !
      !      accumulator"): the fast loop sums the met air temperature over the slow step (reset  !
      !      each step), and the slow-loop phenology driver reads the mean (sum/n) then. Transient !
      !      (never restarted). Air temperature is site-uniform (single-site forcing).            !
      real(wp)           :: pheno_tair_sum = 0.0_wp
      integer(ik)        :: pheno_tair_n   = 0_ik
      !----- Site evapotranspiration accumulator [kg/m2 = mm] (site-uniform, single-site): the fast   !
      !      loop sums the area-weighted canopy-air -> atmosphere water-vapour flux * dt_fast over the  !
      !      slow step (reset each step, mirrors the gpp_accum/pheno lifecycle); read as a diagnostic   !
      !      (total_et) at the output tick. Transient (never restarted).                                !
      real(wp)           :: et_accum       = 0.0_wp
      !----- Integrator WORK accumulators (MEDS_NUMERICS_SCOPING.md section 5.3): summed over patches   !
      !      and fast sub-steps within a slow step, reset alongside et_accum. These are the COST axis   !
      !      of the scheme benchmark -- wall-clock is too coarse and too machine-dependent to rank      !
      !      integrators, whereas step/sub-solver counts are exact and reproducible. Transient          !
      !      (diagnostic only, never restarted). Real-valued so the output integrators can average      !
      !      them over a period like any other diagnostic. --------------------------------------------!
      real(wp)           :: work_integ_steps = 0.0_wp   !< accepted integrator sub-steps
      real(wp)           :: work_integ_rej   = 0.0_wp   !< rejected integrator steps
      real(wp)           :: work_soil_nsub   = 0.0_wp   !< soil-water Richards sub-steps
      real(wp)           :: work_hydro_nsub  = 0.0_wp   !< plant-hydraulics sub-steps (summed over cohorts)
      real(wp)           :: work_nonconv     = 0.0_wp   !< non-convergence events (hydraulics + Picard)
   end type site_t

   !----- A small SoA of the per-cohort carbon-NPP pool fluxes [kgC/plant/step], produced by     !
   !      the carbon rate provider and folded into the tendency bundle by the driver. -----------!
   type :: carbon_flux_block
      real(wp), allocatable :: leaf(:), fineroot(:), wood(:), nonstructural(:)
   end type carbon_flux_block

contains

   !----- Allocate (or resize-on-grow) the tendency bundle to hold n cohorts. Resize-on-grow so   !
   !      the site-carried scratch buffer allocates once and then reuses. The arrays are ALWAYS    !
   !      allocated/freed together, so a partially-allocated bundle (e.g. after a piecemeal edit)   !
   !      forces a full re-allocation rather than a stale reuse. -------------------------------!
   subroutine cohort_deriv_alloc(deriv, n)
      type(cohort_deriv_block), intent(inout) :: deriv
      integer(ik),              intent(in)    :: n
      integer(ik) :: m
      logical     :: ok
      m  = max(n, 1_ik)
      ok = allocated(deriv%d_dbh_dt)             .and. allocated(deriv%d_height_dt)            .and. &
           allocated(deriv%d_basal_area_dt)      .and. allocated(deriv%d_agb_dt)               .and. &
           allocated(deriv%d_leaf_area_dt)       .and. allocated(deriv%d_leaf_carbon_dt)       .and. &
           allocated(deriv%d_fineroot_carbon_dt) .and. allocated(deriv%d_wood_carbon_dt)       .and. &
           allocated(deriv%d_nonstructural_carbon_dt) .and. allocated(deriv%dln_nplant_dt)
      if (ok) then
         if (size(deriv%d_dbh_dt) >= m) return   ! big enough (all co-sized) -> reuse
         deallocate(deriv%d_dbh_dt, deriv%d_height_dt, deriv%d_basal_area_dt, deriv%d_agb_dt,        &
                    deriv%d_leaf_area_dt, deriv%d_leaf_carbon_dt, deriv%d_fineroot_carbon_dt,        &
                    deriv%d_wood_carbon_dt, deriv%d_nonstructural_carbon_dt, deriv%dln_nplant_dt)
      else                                       ! none or partial -> drop whatever is there
         if (allocated(deriv%d_dbh_dt))                  deallocate(deriv%d_dbh_dt)
         if (allocated(deriv%d_height_dt))               deallocate(deriv%d_height_dt)
         if (allocated(deriv%d_basal_area_dt))           deallocate(deriv%d_basal_area_dt)
         if (allocated(deriv%d_agb_dt))                  deallocate(deriv%d_agb_dt)
         if (allocated(deriv%d_leaf_area_dt))            deallocate(deriv%d_leaf_area_dt)
         if (allocated(deriv%d_leaf_carbon_dt))          deallocate(deriv%d_leaf_carbon_dt)
         if (allocated(deriv%d_fineroot_carbon_dt))      deallocate(deriv%d_fineroot_carbon_dt)
         if (allocated(deriv%d_wood_carbon_dt))          deallocate(deriv%d_wood_carbon_dt)
         if (allocated(deriv%d_nonstructural_carbon_dt)) deallocate(deriv%d_nonstructural_carbon_dt)
         if (allocated(deriv%dln_nplant_dt))             deallocate(deriv%dln_nplant_dt)
      end if
      allocate(deriv%d_dbh_dt(m), deriv%d_height_dt(m), deriv%d_basal_area_dt(m), deriv%d_agb_dt(m), &
               deriv%d_leaf_area_dt(m), deriv%d_leaf_carbon_dt(m), deriv%d_fineroot_carbon_dt(m),    &
               deriv%d_wood_carbon_dt(m), deriv%d_nonstructural_carbon_dt(m), deriv%dln_nplant_dt(m))
   end subroutine cohort_deriv_alloc

   !=======================================================================================!
   !  Allocation                                                                            !
   !=======================================================================================!
   subroutine site_alloc(site, n_pft, coh_cap, pat_cap, n_growth_window)
      type(site_t), intent(out) :: site
      integer(ik),     intent(in)  :: n_pft, coh_cap, pat_cap, n_growth_window
      site%n_pft = n_pft
      call cohort_alloc(site%cohort, max(coh_cap, 1_ik), max(n_growth_window, 1_ik))
      call patch_alloc(site%patch, max(pat_cap, 1_ik), n_pft)
   end subroutine site_alloc

   subroutine site_free(site)
      type(site_t), intent(inout) :: site
      site%cohort%n = 0_ik ; site%cohort%cap = 0_ik
      site%patch%n = 0_ik ; site%patch%cap = 0_ik
      if (allocated(site%cohort%nplant)) deallocate(site%cohort%pft, site%cohort%nplant, site%cohort%dbh, &
         site%cohort%height, site%cohort%basal_area, site%cohort%agb, site%cohort%leaf_area,                       &
         site%cohort%growth_avg, site%cohort%growth_accum, site%cohort%growth_count,                     &
         site%cohort%growth_hist, site%cohort%p_dbh_critical, site%cohort%p_wood_density,                &
         site%cohort%p_hgt_max, site%cohort%sla, site%cohort%vcmax25, site%cohort%rd25,               &
         site%cohort%llspan, site%cohort%p_aboveground_frac,                                            &
         site%cohort%p_root_to_leaf_ratio, site%cohort%p_storage_cushion,                                &
         site%cohort%leaf_carbon, site%cohort%fineroot_carbon, site%cohort%wood_carbon,                  &
         site%cohort%nonstructural_carbon, site%cohort%owner_patch, site%cohort%global_id,               &
         site%cohort%overtopping_lai,                                                             &
         site%cohort%leaf_temp, site%cohort%wood_temp, site%cohort%psi, site%cohort%gpp_accum,   &
         site%cohort%leaf_resp_accum, site%cohort%stem_resp_accum, site%cohort%root_resp_accum,  &
         site%cohort%pheno_flush_drive, site%cohort%pheno_shed_drive,                            &
         site%cohort%pheno_gdd, site%cohort%pheno_chill)
      if (allocated(site%patch%area)) deallocate(site%patch%area, site%patch%age, site%patch%dist_type, &
         site%patch%cohort_offset, site%patch%cohort_count, site%patch%recruit_pool, site%patch%global_id, &
         site%patch%cas, site%patch%soil_e, site%patch%soil_w, site%patch%snow, site%patch%soil_carbon, &
         site%patch%xi_accum)
   end subroutine site_free

   subroutine cohort_alloc(cohort, cap, nwin)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: cap, nwin
      cohort%cap = cap ; cohort%n = 0_ik
      allocate(cohort%pft(cap), cohort%owner_patch(cap), cohort%global_id(cap))
      allocate(cohort%nplant(cap), cohort%dbh(cap), cohort%height(cap), cohort%basal_area(cap),            &
               cohort%agb(cap), cohort%leaf_area(cap), cohort%overtopping_lai(cap), cohort%growth_avg(cap),&
               cohort%growth_accum(cap), cohort%growth_count(cap), cohort%growth_hist(nwin, cap))
      allocate(cohort%p_dbh_critical(cap), cohort%p_wood_density(cap), cohort%p_hgt_max(cap))
      allocate(cohort%leaf_carbon(cap), cohort%fineroot_carbon(cap), cohort%wood_carbon(cap),        &
               cohort%nonstructural_carbon(cap))
      allocate(cohort%sla(cap), cohort%p_aboveground_frac(cap), cohort%p_root_to_leaf_ratio(cap),  &
               cohort%p_storage_cushion(cap))
      allocate(cohort%vcmax25(cap), cohort%rd25(cap), cohort%llspan(cap))
      allocate(cohort%leaf_temp(cap), cohort%wood_temp(cap), cohort%psi(N_HYDRO_NODE, cap), cohort%gpp_accum(cap))
      allocate(cohort%leaf_resp_accum(cap), cohort%stem_resp_accum(cap), cohort%root_resp_accum(cap))
      allocate(cohort%pheno_flush_drive(cap), cohort%pheno_shed_drive(cap),                       &
               cohort%pheno_gdd(cap), cohort%pheno_chill(cap))
      cohort%leaf_temp = LEAF_TEMP_INIT ; cohort%wood_temp = LEAF_TEMP_INIT
      cohort%psi = PSI_INIT ; cohort%gpp_accum = 0.0_wp
      cohort%leaf_resp_accum = 0.0_wp ; cohort%stem_resp_accum = 0.0_wp ; cohort%root_resp_accum = 0.0_wp
      cohort%pheno_flush_drive = PHENO_FLUSH_INIT ; cohort%pheno_shed_drive = PHENO_SHED_INIT
      cohort%pheno_gdd = 0.0_wp ; cohort%pheno_chill = 0.0_wp
      cohort%pft = 0_ik ; cohort%owner_patch = 0_ik ; cohort%global_id = 0_ik
      cohort%nplant = 0.0_wp ; cohort%dbh = 0.0_wp ; cohort%height = 0.0_wp ; cohort%basal_area = 0.0_wp
      cohort%agb = 0.0_wp ; cohort%leaf_area = 0.0_wp ; cohort%overtopping_lai = 0.0_wp
      cohort%growth_avg = GROWTH_AVG_UNSET
      cohort%growth_accum = 0.0_wp ; cohort%growth_count = 0_ik ; cohort%growth_hist = 0.0_wp
      cohort%p_dbh_critical = 0.0_wp ; cohort%p_wood_density = 0.0_wp ; cohort%p_hgt_max = 0.0_wp
      cohort%leaf_carbon = 0.0_wp ; cohort%fineroot_carbon = 0.0_wp ; cohort%wood_carbon = 0.0_wp
      cohort%nonstructural_carbon = 0.0_wp
      cohort%sla = 0.0_wp ; cohort%p_aboveground_frac = 0.0_wp
      cohort%p_root_to_leaf_ratio = 0.0_wp ; cohort%p_storage_cushion = 0.0_wp
      cohort%vcmax25 = 0.0_wp ; cohort%rd25 = 0.0_wp ; cohort%llspan = 0.0_wp
   end subroutine cohort_alloc

   subroutine patch_alloc(patch, cap, n_pft)
      type(patch_block), intent(inout) :: patch
      integer(ik),       intent(in)    :: cap, n_pft
      patch%cap = cap ; patch%n = 0_ik
      allocate(patch%area(cap), patch%age(cap), patch%dist_type(cap), patch%global_id(cap))
      allocate(patch%cohort_offset(cap), patch%cohort_count(cap), patch%recruit_pool(n_pft, cap))
      allocate(patch%cas(cap), patch%soil_e(cap), patch%soil_w(cap), patch%snow(cap))   !< default-initialised reservoirs
      allocate(patch%soil_carbon(cap))                                                 !< default-initialised (0)
      allocate(patch%xi_accum(cap))                                                    !< default-initialised (0)
      patch%area = 0.0_wp ; patch%age = 0.0_wp ; patch%dist_type = 1_ik ; patch%global_id = 0_ik
      patch%cohort_offset = 0_ik ; patch%cohort_count = 0_ik ; patch%recruit_pool = 0.0_wp
   end subroutine patch_alloc

   !=======================================================================================!
   !  cohort_ensure_capacity -- ensure the cohort SoA can hold at least `need` cohorts,    !
   !  growing it only if it must. This is the growable-array primitive (like C++           !
   !  vector::reserve) for the flat cohort block, where `n` is the LIVE count and          !
   !  `cap` the ALLOCATED count (slots n+1..cap are dead reserve). Callers about to        !
   !  APPEND cohorts call it first; it NO-OPs when room already exists (need <= cap),      !
   !  so it is safe to call unconditionally; otherwise it grows GEOMETRICALLY (1.5x)       !
   !  for amortized-O(1) appends, copying the live 1:n prefix and move_alloc-ing in.       !
   !  LOCKSTEP SITE: a new per-cohort field must be added to the 1:n copy list AND to      !
   !  move_alloc_block (in both = kept on grow; move_alloc_block only = reset).            !
   !=======================================================================================!
   subroutine cohort_ensure_capacity(cohort, need)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: need
      type(cohort_block)                :: tmp
      integer(ik)                       :: newcap, m
      if (need <= cohort%cap) return
      newcap = max(need, (cohort%cap * 3_ik) / 2_ik + 1_ik)
      call cohort_alloc(tmp, newcap, int(size(cohort%growth_hist, 1), ik))  ! keep the window size
      m = cohort%n
      tmp%n = m
      tmp%pft(1:m)            = cohort%pft(1:m)
      tmp%owner_patch(1:m)    = cohort%owner_patch(1:m)
      tmp%nplant(1:m)         = cohort%nplant(1:m)
      tmp%dbh(1:m)            = cohort%dbh(1:m)
      tmp%height(1:m)           = cohort%height(1:m)
      tmp%basal_area(1:m)        = cohort%basal_area(1:m)
      tmp%agb(1:m)            = cohort%agb(1:m)
      tmp%leaf_area(1:m)          = cohort%leaf_area(1:m)
      tmp%growth_avg(1:m)     = cohort%growth_avg(1:m)
      tmp%growth_accum(1:m)   = cohort%growth_accum(1:m)
      tmp%growth_count(1:m)   = cohort%growth_count(1:m)
      tmp%growth_hist(:,1:m)  = cohort%growth_hist(:,1:m)
      tmp%p_dbh_critical(1:m)     = cohort%p_dbh_critical(1:m)
      tmp%p_wood_density(1:m) = cohort%p_wood_density(1:m)
      tmp%p_hgt_max(1:m)      = cohort%p_hgt_max(1:m)
      tmp%leaf_carbon(1:m)           = cohort%leaf_carbon(1:m)
      tmp%fineroot_carbon(1:m)       = cohort%fineroot_carbon(1:m)
      tmp%wood_carbon(1:m)           = cohort%wood_carbon(1:m)
      tmp%nonstructural_carbon(1:m)  = cohort%nonstructural_carbon(1:m)
      tmp%sla(1:m)                 = cohort%sla(1:m)
      tmp%vcmax25(1:m)               = cohort%vcmax25(1:m)
      tmp%rd25(1:m)                  = cohort%rd25(1:m)
      tmp%llspan(1:m)                = cohort%llspan(1:m)
      tmp%p_aboveground_frac(1:m)    = cohort%p_aboveground_frac(1:m)
      tmp%p_root_to_leaf_ratio(1:m)  = cohort%p_root_to_leaf_ratio(1:m)
      tmp%p_storage_cushion(1:m)     = cohort%p_storage_cushion(1:m)
      tmp%global_id(1:m)      = cohort%global_id(1:m)
      tmp%leaf_temp(1:m)      = cohort%leaf_temp(1:m)
      tmp%wood_temp(1:m)      = cohort%wood_temp(1:m)
      tmp%psi(:,1:m)          = cohort%psi(:,1:m)
      tmp%gpp_accum(1:m)      = cohort%gpp_accum(1:m)
      tmp%leaf_resp_accum(1:m) = cohort%leaf_resp_accum(1:m)
      tmp%stem_resp_accum(1:m) = cohort%stem_resp_accum(1:m)
      tmp%root_resp_accum(1:m) = cohort%root_resp_accum(1:m)
      tmp%pheno_flush_drive(1:m) = cohort%pheno_flush_drive(1:m)
      tmp%pheno_shed_drive(1:m)  = cohort%pheno_shed_drive(1:m)
      tmp%pheno_gdd(1:m)        = cohort%pheno_gdd(1:m)
      tmp%pheno_chill(1:m)      = cohort%pheno_chill(1:m)
      call move_alloc_block(tmp, cohort)
   end subroutine cohort_ensure_capacity

   !----- Move every allocatable from src into dst (src is reset). ------------------------!
   subroutine move_alloc_block(src, dst)
      type(cohort_block), intent(inout) :: src
      type(cohort_block), intent(inout) :: dst
      dst%n = src%n ; dst%cap = src%cap
      call move_alloc(src%pft, dst%pft)
      call move_alloc(src%owner_patch, dst%owner_patch)
      call move_alloc(src%nplant, dst%nplant)
      call move_alloc(src%dbh, dst%dbh)
      call move_alloc(src%height, dst%height)
      call move_alloc(src%basal_area, dst%basal_area)
      call move_alloc(src%agb, dst%agb)
      call move_alloc(src%leaf_area, dst%leaf_area)
      call move_alloc(src%overtopping_lai, dst%overtopping_lai)
      call move_alloc(src%growth_avg, dst%growth_avg)
      call move_alloc(src%growth_accum, dst%growth_accum)
      call move_alloc(src%growth_count, dst%growth_count)
      call move_alloc(src%growth_hist, dst%growth_hist)
      call move_alloc(src%p_dbh_critical, dst%p_dbh_critical)
      call move_alloc(src%p_wood_density, dst%p_wood_density)
      call move_alloc(src%p_hgt_max, dst%p_hgt_max)
      call move_alloc(src%leaf_carbon, dst%leaf_carbon)
      call move_alloc(src%fineroot_carbon, dst%fineroot_carbon)
      call move_alloc(src%wood_carbon, dst%wood_carbon)
      call move_alloc(src%nonstructural_carbon, dst%nonstructural_carbon)
      call move_alloc(src%sla, dst%sla)
      call move_alloc(src%vcmax25, dst%vcmax25)
      call move_alloc(src%rd25, dst%rd25)
      call move_alloc(src%llspan, dst%llspan)
      call move_alloc(src%p_aboveground_frac, dst%p_aboveground_frac)
      call move_alloc(src%p_root_to_leaf_ratio, dst%p_root_to_leaf_ratio)
      call move_alloc(src%p_storage_cushion, dst%p_storage_cushion)
      call move_alloc(src%global_id, dst%global_id)
      call move_alloc(src%leaf_temp, dst%leaf_temp)
      call move_alloc(src%wood_temp, dst%wood_temp)
      call move_alloc(src%psi, dst%psi)
      call move_alloc(src%gpp_accum, dst%gpp_accum)
      call move_alloc(src%leaf_resp_accum, dst%leaf_resp_accum)
      call move_alloc(src%stem_resp_accum, dst%stem_resp_accum)
      call move_alloc(src%root_resp_accum, dst%root_resp_accum)
      call move_alloc(src%pheno_flush_drive, dst%pheno_flush_drive)
      call move_alloc(src%pheno_shed_drive, dst%pheno_shed_drive)
      call move_alloc(src%pheno_gdd, dst%pheno_gdd)
      call move_alloc(src%pheno_chill, dst%pheno_chill)
   end subroutine move_alloc_block

   !=======================================================================================!
   !  patch_ensure_capacity -- the patch-SoA twin of cohort_ensure_capacity: ensure the    !
   !  patch SoA can hold at least `need` patches, growing 1.5x if it must. Same n/cap      !
   !  reserve model, no-op fast path, geometric growth, and move_alloc as in the           !
   !  cohort routine. LOCKSTEP SITE: a new per-patch field must be copied here (1:n)       !
   !  and move_alloc'd, then registered at sort_patches / patch_compact -- the patch       !
   !  lockstep sites, since patches have no single reorder routine like cohorts do.        !
   !=======================================================================================!
   subroutine patch_ensure_capacity(patch, need, n_pft)
      type(patch_block), intent(inout) :: patch
      integer(ik),       intent(in)    :: need, n_pft
      type(patch_block)                :: tmp
      integer(ik)                       :: newcap, m
      if (need <= patch%cap) return
      newcap = max(need, (patch%cap * 3_ik) / 2_ik + 1_ik)
      call patch_alloc(tmp, newcap, n_pft)
      m = patch%n ; tmp%n = m
      tmp%area(1:m)           = patch%area(1:m)
      tmp%age(1:m)            = patch%age(1:m)
      tmp%dist_type(1:m)      = patch%dist_type(1:m)
      tmp%cohort_offset(1:m)           = patch%cohort_offset(1:m)
      tmp%cohort_count(1:m)         = patch%cohort_count(1:m)
      tmp%recruit_pool(:,1:m) = patch%recruit_pool(:,1:m)
      tmp%global_id(1:m)      = patch%global_id(1:m)
      tmp%cas(1:m)            = patch%cas(1:m)
      tmp%soil_e(1:m)         = patch%soil_e(1:m)
      tmp%soil_w(1:m)         = patch%soil_w(1:m)
      tmp%snow(1:m)           = patch%snow(1:m)
      tmp%soil_carbon(1:m)    = patch%soil_carbon(1:m)
      tmp%xi_accum(1:m)       = patch%xi_accum(1:m)
      patch%n = tmp%n ; patch%cap = tmp%cap
      call move_alloc(tmp%area, patch%area)             ; call move_alloc(tmp%age, patch%age)
      call move_alloc(tmp%dist_type, patch%dist_type)
      call move_alloc(tmp%cohort_offset, patch%cohort_offset)             ; call move_alloc(tmp%cohort_count, patch%cohort_count)
      call move_alloc(tmp%recruit_pool, patch%recruit_pool)
      call move_alloc(tmp%global_id, patch%global_id)
      call move_alloc(tmp%cas, patch%cas) ; call move_alloc(tmp%soil_e, patch%soil_e)
      call move_alloc(tmp%soil_w, patch%soil_w) ; call move_alloc(tmp%snow, patch%snow)
      call move_alloc(tmp%soil_carbon, patch%soil_carbon)
      call move_alloc(tmp%xi_accum, patch%xi_accum)
   end subroutine patch_ensure_capacity

   !=======================================================================================!
   !  Centralized lockstep reorder: place old index perm(k) at new position k, for k=1..m.  !
   !  Array assignment with a vector subscript on the RHS is evaluated to a temporary       !
   !  first, so self-permutation is safe. THIS is the one place that knows every field.     !
   !=======================================================================================!
   subroutine cohort_reorder(cohort, perm, m)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: perm(:)
      integer(ik),        intent(in)    :: m
      cohort%pft(1:m)            = cohort%pft(perm(1:m))
      cohort%owner_patch(1:m)    = cohort%owner_patch(perm(1:m))
      cohort%nplant(1:m)         = cohort%nplant(perm(1:m))
      cohort%dbh(1:m)            = cohort%dbh(perm(1:m))
      cohort%height(1:m)           = cohort%height(perm(1:m))
      cohort%basal_area(1:m)        = cohort%basal_area(perm(1:m))
      cohort%agb(1:m)            = cohort%agb(perm(1:m))
      cohort%leaf_area(1:m)          = cohort%leaf_area(perm(1:m))
      cohort%overtopping_lai(1:m)    = cohort%overtopping_lai(perm(1:m))
      cohort%growth_avg(1:m)     = cohort%growth_avg(perm(1:m))
      cohort%growth_accum(1:m)   = cohort%growth_accum(perm(1:m))
      cohort%growth_count(1:m)   = cohort%growth_count(perm(1:m))
      cohort%growth_hist(:,1:m)  = cohort%growth_hist(:,perm(1:m))
      cohort%p_dbh_critical(1:m)     = cohort%p_dbh_critical(perm(1:m))
      cohort%p_wood_density(1:m) = cohort%p_wood_density(perm(1:m))
      cohort%p_hgt_max(1:m)      = cohort%p_hgt_max(perm(1:m))
      cohort%leaf_carbon(1:m)           = cohort%leaf_carbon(perm(1:m))
      cohort%fineroot_carbon(1:m)       = cohort%fineroot_carbon(perm(1:m))
      cohort%wood_carbon(1:m)           = cohort%wood_carbon(perm(1:m))
      cohort%nonstructural_carbon(1:m)  = cohort%nonstructural_carbon(perm(1:m))
      cohort%sla(1:m)                 = cohort%sla(perm(1:m))
      cohort%vcmax25(1:m)               = cohort%vcmax25(perm(1:m))
      cohort%rd25(1:m)                  = cohort%rd25(perm(1:m))
      cohort%llspan(1:m)                = cohort%llspan(perm(1:m))
      cohort%p_aboveground_frac(1:m)    = cohort%p_aboveground_frac(perm(1:m))
      cohort%p_root_to_leaf_ratio(1:m)  = cohort%p_root_to_leaf_ratio(perm(1:m))
      cohort%p_storage_cushion(1:m)     = cohort%p_storage_cushion(perm(1:m))
      cohort%global_id(1:m)      = cohort%global_id(perm(1:m))
      cohort%leaf_temp(1:m)      = cohort%leaf_temp(perm(1:m))
      cohort%wood_temp(1:m)      = cohort%wood_temp(perm(1:m))
      cohort%psi(:,1:m)          = cohort%psi(:,perm(1:m))
      cohort%gpp_accum(1:m)      = cohort%gpp_accum(perm(1:m))
      cohort%leaf_resp_accum(1:m) = cohort%leaf_resp_accum(perm(1:m))
      cohort%stem_resp_accum(1:m) = cohort%stem_resp_accum(perm(1:m))
      cohort%root_resp_accum(1:m) = cohort%root_resp_accum(perm(1:m))
      cohort%pheno_flush_drive(1:m) = cohort%pheno_flush_drive(perm(1:m))
      cohort%pheno_shed_drive(1:m)  = cohort%pheno_shed_drive(perm(1:m))
      cohort%pheno_gdd(1:m)        = cohort%pheno_gdd(perm(1:m))
      cohort%pheno_chill(1:m)      = cohort%pheno_chill(perm(1:m))
      cohort%n = m
   end subroutine cohort_reorder

   !----- Compact by keep-mask (length n): drop entries where keep is .false. -------------!
   subroutine cohort_compact(cohort, keep)
      type(cohort_block), intent(inout) :: cohort
      logical,            intent(in)    :: keep(:)
      integer(ik), allocatable :: perm(:)
      integer(ik)              :: i, m
      allocate(perm(cohort%n)) ; m = 0_ik
      do i = 1_ik, cohort%n
         if (keep(i)) then
            m = m + 1_ik
            perm(m) = i
         end if
      end do
      call cohort_reorder(cohort, perm, m)
   end subroutine cohort_compact

   !----- Copy every per-cohort field from slot src to slot dst (centralized). ------------!
   subroutine copy_cohort_slot(cohort, dst, src)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: dst, src
      cohort%pft(dst)            = cohort%pft(src)
      cohort%owner_patch(dst)    = cohort%owner_patch(src)
      cohort%nplant(dst)         = cohort%nplant(src)
      cohort%dbh(dst)            = cohort%dbh(src)
      cohort%height(dst)           = cohort%height(src)
      cohort%basal_area(dst)        = cohort%basal_area(src)
      cohort%agb(dst)            = cohort%agb(src)
      cohort%leaf_area(dst)          = cohort%leaf_area(src)
      cohort%overtopping_lai(dst)    = cohort%overtopping_lai(src)
      cohort%growth_avg(dst)     = cohort%growth_avg(src)
      cohort%growth_accum(dst)   = cohort%growth_accum(src)
      cohort%growth_count(dst)   = cohort%growth_count(src)
      cohort%growth_hist(:,dst)  = cohort%growth_hist(:,src)
      cohort%p_dbh_critical(dst)     = cohort%p_dbh_critical(src)
      cohort%p_wood_density(dst) = cohort%p_wood_density(src)
      cohort%p_hgt_max(dst)      = cohort%p_hgt_max(src)
      cohort%leaf_carbon(dst)           = cohort%leaf_carbon(src)
      cohort%fineroot_carbon(dst)       = cohort%fineroot_carbon(src)
      cohort%wood_carbon(dst)           = cohort%wood_carbon(src)
      cohort%nonstructural_carbon(dst)  = cohort%nonstructural_carbon(src)
      cohort%sla(dst)                 = cohort%sla(src)
      cohort%vcmax25(dst)               = cohort%vcmax25(src)
      cohort%rd25(dst)                  = cohort%rd25(src)
      cohort%llspan(dst)                = cohort%llspan(src)
      cohort%p_aboveground_frac(dst)    = cohort%p_aboveground_frac(src)
      cohort%p_root_to_leaf_ratio(dst)  = cohort%p_root_to_leaf_ratio(src)
      cohort%p_storage_cushion(dst)     = cohort%p_storage_cushion(src)
      cohort%global_id(dst)      = cohort%global_id(src)
      cohort%leaf_temp(dst)      = cohort%leaf_temp(src)
      cohort%wood_temp(dst)      = cohort%wood_temp(src)
      cohort%psi(:,dst)          = cohort%psi(:,src)
      cohort%gpp_accum(dst)      = cohort%gpp_accum(src)
      cohort%leaf_resp_accum(dst) = cohort%leaf_resp_accum(src)
      cohort%stem_resp_accum(dst) = cohort%stem_resp_accum(src)
      cohort%root_resp_accum(dst) = cohort%root_resp_accum(src)
      cohort%pheno_flush_drive(dst) = cohort%pheno_flush_drive(src)
      cohort%pheno_shed_drive(dst)  = cohort%pheno_shed_drive(src)
      cohort%pheno_gdd(dst)        = cohort%pheno_gdd(src)
      cohort%pheno_chill(dst)      = cohort%pheno_chill(src)
   end subroutine copy_cohort_slot

   !----- Fill the gathered per-cohort PFT params from the trait table. -------------------!
   subroutine gather_pft_params(cohort, pft)
      type(cohort_block), intent(inout) :: cohort
      type(pft_table_t),  intent(in)    :: pft
      integer(ik) :: i, p
      do i = 1_ik, cohort%n
         p = cohort%pft(i)
         cohort%p_dbh_critical(i)     = pft%dbh_critical(p)
         cohort%p_wood_density(i) = pft%wood_density(p)
         cohort%p_hgt_max(i)      = pft%hgt_max(p)
         cohort%p_aboveground_frac(i)   = pft%aboveground_frac(p)
         cohort%p_root_to_leaf_ratio(i) = pft%root_to_leaf_ratio(p)
         cohort%p_storage_cushion(i)    = pft%storage_cushion(p)
      end do
   end subroutine gather_pft_params

   !---------------------------------------------------------------------------------------!
   ! Re-derive the cached geometry (height, basal area, AGB, leaf area) of ONE cohort slot   !
   ! from its prognostic diameter and gathered wood density. Single-slot host analogue of    !
   ! the array math in growth_step; used by recruitment, setup, fusion and fission so the     !
   ! allometry lives in exactly one (shared) place.                                          !
   !---------------------------------------------------------------------------------------!
   subroutine set_cohort_size(cohort, i)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: i
      cohort%height(i)      = dbh_to_height(cohort%dbh(i), cohort%p_hgt_max(i))
      cohort%basal_area(i)  = pio4 * cohort%dbh(i) * cohort%dbh(i)
      cohort%agb(i)         = dbh_to_agb(cohort%dbh(i), cohort%height(i), cohort%p_wood_density(i))
      cohort%leaf_area(i)       = dbh_to_leaf_area(cohort%dbh(i), cohort%height(i))
      !----- On-allometry carbon pools (PR3 cached diagnostics; reuse the just-computed agb &      !
      !      leaf_area, so leaf_carbon*sla = leaf_area and wood_carbon*aboveground_frac = agb). -----!
      cohort%leaf_carbon(i)          = cohort%leaf_area(i) / max(cohort%sla(i), tiny_num)
      cohort%wood_carbon(i)          = cohort%agb(i) / max(cohort%p_aboveground_frac(i), tiny_num)
      cohort%fineroot_carbon(i)      = cohort%p_root_to_leaf_ratio(i) * cohort%leaf_carbon(i)
      cohort%nonstructural_carbon(i) = cohort%p_storage_cushion(i) * cohort%leaf_carbon(i)
   end subroutine set_cohort_size

   !---------------------------------------------------------------------------------------!
   ! Initialize ONE freshly-created cohort slot at BIRTH from (pft, patch, dbh): set the      !
   ! per-slot identity + prognostic + gathered-PFT fields to their birth values, then derive  !
   ! the on-allometry geometry via set_cohort_size. The SINGLE canonical birth-field policy    !
   ! shared by add_cohort (setup / census) and apply_recruitment -- add any new per-cohort     !
   ! field that needs a birth value HERE. overtopping_lai is zeroed (a reused, culled slot may  !
   ! hold a stale value; it is recomputed by update_overtopping_lai every slow step anyway).    !
   ! The caller ensures capacity, stamps the global id, and bumps cohort%n.                    !
   !---------------------------------------------------------------------------------------!
   subroutine init_cohort(cohort, m, pft, ipft, owner_patch, nplant, dbh)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: m, ipft, owner_patch
      type(pft_table_t),  intent(in)    :: pft
      real(wp),           intent(in)    :: nplant, dbh
      cohort%pft(m)              = ipft
      cohort%owner_patch(m)      = owner_patch
      cohort%nplant(m)           = nplant
      cohort%dbh(m)              = dbh
      cohort%growth_avg(m)       = GROWTH_AVG_UNSET   ! set on its first growth step
      cohort%growth_accum(m)     = 0.0_wp
      cohort%growth_count(m)     = 0_ik
      cohort%overtopping_lai(m)  = 0.0_wp             ! fresh competition context (recomputed each slow step)
      cohort%pheno_flush_drive(m) = PHENO_FLUSH_INIT  ! born flushing (evergreen fixed point)
      cohort%pheno_shed_drive(m)  = PHENO_SHED_INIT   ! no active shed at birth
      cohort%pheno_gdd(m)        = 0.0_wp             ! fresh phenology memory
      cohort%pheno_chill(m)      = 0.0_wp
      cohort%p_dbh_critical(m)       = pft%dbh_critical(ipft)
      cohort%p_wood_density(m)       = pft%wood_density(ipft)
      cohort%p_hgt_max(m)            = pft%hgt_max(ipft)
      cohort%sla(m)                = pft%sla(ipft)
      cohort%vcmax25(m)              = pft%vcmax25(ipft)          ! dynamic leaf traits born at
      cohort%rd25(m)                 = pft%rd25(ipft)             ! top-of-canopy; plasticity relaxes
      cohort%llspan(m)               = pft%leaf_lifespan_toc(ipft)! them toward the shaded target
      cohort%p_aboveground_frac(m)   = pft%aboveground_frac(ipft)
      cohort%p_root_to_leaf_ratio(m) = pft%root_to_leaf_ratio(ipft)
      cohort%p_storage_cushion(m)    = pft%storage_cushion(ipft)
      cohort%leaf_temp(m)        = LEAF_TEMP_INIT     ! fresh fast state (slot may be a reused, stale cull)
      cohort%wood_temp(m)        = LEAF_TEMP_INIT     ! ditto -- reset like cohort_alloc, else a reused slot keeps a dead cohort's wood_temp
      cohort%psi(:,m)            = PSI_INIT
      call set_cohort_size(cohort, m)                 ! height/basal_area/agb/leaf_area + carbon pools from dbh
   end subroutine init_cohort

   !---------------------------------------------------------------------------------------!
   ! CARBON-MODE geometry (the inverse of set_cohort_size): with wood_carbon the prognostic    !
   ! size anchor, derive dbh from it (wood_to_dbh) and the rest of the cached geometry, and take !
   ! leaf_area straight from the prognostic leaf_carbon. Used by apply_growth. The carbon     !
   ! pools are the INPUTS here, so they are NOT re-derived (unlike set_cohort_size).             !
   !---------------------------------------------------------------------------------------!
   subroutine set_cohort_size_from_carbon(cohort, i)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: i
      !----- The carbon->geometry flip is now the single shared kernel carbon_to_structure       !
      !       (meds_allometry); this per-slot wrapper just packs/unpacks the SoA fields.           !
      call carbon_to_structure(cohort%wood_carbon(i), cohort%leaf_carbon(i),                      &
                               cohort%p_wood_density(i), cohort%p_hgt_max(i),                     &
                               cohort%p_aboveground_frac(i), cohort%sla(i),                     &
                               cohort%dbh(i), cohort%height(i), cohort%basal_area(i),             &
                               cohort%agb(i), cohort%leaf_area(i))
   end subroutine set_cohort_size_from_carbon

   !=======================================================================================!
   !  Rebuild the CSR patch map by stably regrouping cohorts by owner_patch (counting sort).!
   !=======================================================================================!
   subroutine rebuild_csr(site)
      type(site_t), intent(inout) :: site
      integer(ik), allocatable :: perm(:), pos(:), slot(:)
      integer(ik)              :: i, ip, np, nc
      np = site%patch%n ; nc = site%cohort%n
      if (np < 1_ik) then                                  ! no patches: nothing to group
         if (nc > 0_ik) error stop 'rebuild_csr: cohorts present with no patches'
         return
      end if
      site%patch%cohort_count(1:np) = 0_ik
      do i = 1_ik, nc
         ip = site%cohort%owner_patch(i)
         if (ip < 1_ik .or. ip > np) error stop 'rebuild_csr: owner_patch out of range'
         site%patch%cohort_count(ip) = site%patch%cohort_count(ip) + 1_ik
      end do
      site%patch%cohort_offset(1) = 1_ik
      do ip = 2_ik, np
         site%patch%cohort_offset(ip) = site%patch%cohort_offset(ip-1) + site%patch%cohort_count(ip-1)
      end do
      !----- Stable placement: preserve within-patch order. -------------------------------!
      allocate(pos(np), slot(nc), perm(nc))
      pos(1:np) = site%patch%cohort_offset(1:np)
      do i = 1_ik, nc
         ip = site%cohort%owner_patch(i)
         slot(i) = pos(ip)
         pos(ip) = pos(ip) + 1_ik
      end do
      do i = 1_ik, nc
         perm(slot(i)) = i
      end do
      call cohort_reorder(site%cohort, perm, nc)
   end subroutine rebuild_csr

   !=======================================================================================!
   !  Persistent identity: stamp a slot with the next unused global id and bump the         !
   !  monotonic counter. Ids are never reused within a run, so a terminated/fused-away id    !
   !  simply never reappears -- exactly what an external tracker needs.                      !
   !=======================================================================================!
   subroutine assign_cohort_id(site, slot)
      type(site_t), intent(inout) :: site
      integer(ik),  intent(in)    :: slot
      site%cohort%global_id(slot) = site%next_cohort_id
      site%next_cohort_id = site%next_cohort_id + 1_ik
   end subroutine assign_cohort_id

   subroutine assign_patch_id(site, slot)
      type(site_t), intent(inout) :: site
      integer(ik),  intent(in)    :: slot
      site%patch%global_id(slot) = site%next_patch_id
      site%next_patch_id = site%next_patch_id + 1_ik
   end subroutine assign_patch_id

end module meds_core_state_types
