! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_cohort_state_vars -- the state of every cohort of a site, as one flat               !
! Structure-of-Arrays (`cohort_block`), and every routine that lists the cohort fields.    !
!                                                                                          !
! ALL cohorts of the site live in one contiguous set of 1-D arrays, so the dominant daily  !
! kernels are a single unit-stride sweep, ideal for SIMD and GPU. The patch block's CSR    !
! map (meds_patch_state_vars) groups them: the cohorts of a patch occupy one slice.        !
!                                                                                          !
! Every structural change goes through ONE centralized routine, `cohort_reorder`, which    !
! permutes EVERY per-cohort array in lockstep -- this is the single place to update when   !
! a field is added, eliminating the ED2 "forgot to reallocate an array" class of bug.      !
!                                                                                          !
! The name follows ED2's ed_state_vars, which holds the same kind of types with their      !
! allocate and copy routines. Unlike ED2 there is no global state (no edgrid_g): this      !
! module defines the types and the routines on them, and the state is passed as an         !
! argument.                                                                                !
!==========================================================================================!
module meds_cohort_state_vars
   use meds_diag_vars,  only : cohort_diag_block, cohort_diag_grow, cohort_diag_reorder,         &
                               cohort_diag_copy_slot, cohort_diag_clear_slot
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : cp_liq, C2B_WOOD, TISSUE_C_LEAF, TISSUE_C_WOOD_DRY
   use meds_water_retention, only : wood_thermal_water, water_curve_t, water_content
   use meds_numerics,   only : pair_mean
   use meds_pft_params, only : pft_table_t
   use meds_allometry,  only : allometry_t, dbh_to_height, size2leaf_carbon, size2wood_carbon,           &
                               carbon_to_structure, dbh_to_wai, sapwood_fraction, storage_carbon_target
   use meds_config,     only : hydraulics_config_t, pft_tissue_curves
   implicit none
   private

   public :: cohort_block, cohort_alloc, cohort_free, cohort_ensure_capacity, cohort_reorder, cohort_compact
   public :: copy_cohort_slot, set_cohort_size, set_cohort_size_from_carbon, init_cohort
   public :: scale_cohort_ground_fields, fuse_cohort_fast_state
   public :: cohort_deriv_block, cohort_deriv_alloc
   public :: cohort_tissue_heat_capacity, cohort_tissue_water, cohort_wood_store_carbon, set_cohort_tissue_water
   public :: cohort_bad_field, cohort_field_name

   public :: GROWTH_AVG_UNSET, DMAX_PSI_LEAF_UNSET, DMAX_PSI_LEAF_ACCUM_RESET, PHENO_FLUSH_INIT, PHENO_SHED_INIT

   !----- Sentinel for a not-yet-sampled moving-average growth: a freshly created cohort holds  !
   !      this until its first growth step adds a sample (negative => "unset"; a real growth     !
   !      average is >= 0). Used by the growth kernel and the rate evaluator.                    !
   real(wp), parameter :: GROWTH_AVG_UNSET = -1.0_wp
   !----- Sentinel: 'no predawn history yet'. Positive, so it is unmistakable -- a real leaf water  !
   !      potential is always <= 0. The fast loop replaces it with the surface-layer soil potential. !
   real(wp), parameter :: DMAX_PSI_LEAF_UNSET = 1.0_wp
   !----- Reset value for the running daily MAX. Must be VERY NEGATIVE, not the UNSET sentinel:      !
   !      the accumulator is a max() against psi_leaf <= 0, so seeding it with the POSITIVE sentinel  !
   !      would pin it there forever and the daily max would never update.                            !
   real(wp), parameter :: DMAX_PSI_LEAF_ACCUM_RESET = -1.0e30_wp

   !----- Initial phenology tendencies of a freshly created cohort: born flushing, not senescing   !
   !      (the always-flushing fixed point), so a fresh canopy keeps building. Plain literals       !
   !      because `state` may not depend on the plant library (the DAG wall); the phenology driver   !
   !      advances them once real drivers arrive.                                                    !
   real(wp), parameter :: PHENO_FLUSH_INIT = 1.0_wp
   real(wp), parameter :: PHENO_SHED_INIT  = 0.0_wp

   !---------------------------------------------------------------------------------------!
   ! All cohorts of the site_t (contiguous SoA). The demographic state is nplant and the four   !
   ! carbon pools; every slow operator moves only those. The geometry (dbh, height, basal_area, !
   ! agb, leaf_area and the wood caches) is a function of the pools, stored so the hot kernels  !
   ! stay arithmetic-only, and written only by set_cohort_size_from_carbon -- after every       !
   ! operator that moves carbon, never integrated or read back on its own.                      !
   !---------------------------------------------------------------------------------------!
   type :: cohort_block
      !----- The site's allometry coefficients (cfg%allom of the config that built the site, set by  !
      !      site_alloc). One record for the whole block: every cohort of a site shares it, and the   !
      !      polygons of a region each hold their own. Parameters, not state: never in a state file.  !
      type(allometry_t) :: allom
      !----- The site's PFT trait table (cfg%pft of the config that built the site, set by site_alloc  !
      !      beside allom). A cohort's PFT-level traits are read from it by the cohort's PFT index,    !
      !      pft_table%x(pft(i)), so no trait is copied onto every cohort. Parameters, not state: never !
      !      in a state file.                                                                           !
      type(pft_table_t) :: pft_table
      integer(ik) :: n   = 0_ik       !< cohorts in use
      integer(ik) :: cap = 0_ik       !< allocated capacity
      !----- Identity, plant density, and the geometry derived from the carbon pools. --------!
      integer(ik), allocatable :: pft(:)            !< PFT index
      real(wp),    allocatable :: nplant(:)         !< [plant/m2] PRIMARY demographic state
      real(wp),    allocatable :: dbh(:)            !< [cm]       = wood_to_dbh(wood_carbon), derived
      real(wp),    allocatable :: height(:)           !< [m]        = dbh_to_height(dbh), derived
      real(wp),    allocatable :: basal_area(:)        !< [cm2/plant]= pio4*dbh^2, derived
      real(wp),    allocatable :: agb(:)             !< [kgC/plant] = aboveground_frac*wood_carbon, derived
      real(wp),    allocatable :: leaf_area(:)           !< [m2/plant]  = sla*leaf_carbon, derived (LAI=nplant*leaf_area)
      !----- The wood twins of leaf_area, derived the same way and for the same reason: the fast   !
      !      loop needs them every dt_fast and they are pure functions of dbh and the PFT's        !
      !      traits. Kept PER PLANT, exactly like leaf_area, so nothing here goes stale when        !
      !      mortality changes nplant -- the per-ground index is nplant*area at the point of use     !
      !      (WAI = nplant*wood_area, mirroring LAI = nplant*leaf_area).                             !
      real(wp),    allocatable :: wood_area(:)            !< [m2/plant]  wood area, cached (WAI=nplant*wood_area)
      real(wp),    allocatable :: sapwood_carbon(:)       !< [kgC/plant] sapwood ring (hydraulic capacitance)
      real(wp),    allocatable :: sapwood_area(:)         !< [cm2/plant] sapwood cross-section (sap-flow path)
      real(wp),    allocatable :: overtopping_lai(:)     !< [m2/m2] cumulative LAI of all TALLER cohorts in the patch
                                                         !<         (Beer competition context); a RECOMPUTED diagnostic
                                                         !<         filled by meds_competition%update_overtopping_lai
      !----- Carbon pools [kgC/plant]: the prognostic state. wood_carbon is the size anchor; a     !
      !      cohort placed on the allometry (set_cohort_size) starts with the pools of its dbh.     !
      real(wp),    allocatable :: leaf_carbon(:)          !< [kgC/plant] prognostic; on allometry = leaf_area(dbh)/sla
      real(wp),    allocatable :: fineroot_carbon(:)      !< [kgC/plant] prognostic; on allometry = root_to_leaf_ratio*leaf_carbon
      real(wp),    allocatable :: wood_carbon(:)          !< [kgC/plant] prognostic size anchor
      real(wp),    allocatable :: nonstructural_carbon(:) !< [kgC/plant] prognostic; starts at storage_cushion*(leaf + fine root)
      !----- Sliding simple-moving-average growth rate [cm/yr] driving growth-dependent          !
      !      mortality. growth_hist is the per-cohort ring buffer of the last growth_window steps !
      !      (column per cohort; the write slot is the site-global growth_hist_pos), growth_accum !
      !      its sum, taken from the buffer each step, and growth_count its fill level (<= window); !
      !      growth_avg = accum/count is the SMA (GROWTH_AVG_UNSET until the first sample). A      !
      !      PROGNOSTIC history (not derivable from dbh).                                         !
      real(wp),    allocatable :: growth_avg(:)
      !----- PLANT WATER STATUS driving the stomatal stress limb                                      !
      !      (docs/science/leaf_gas_exchange.md §4.2, eq. (10)). dmax_psi_leaf is YESTERDAY's         !
      !      daily-MAXIMUM leaf water potential. The daily max occurs near dawn, after                !
      !      transpiration has stopped and the plant has relaxed toward the soil, so it is the        !
      !      model's analogue of measured PREDAWN leaf water potential -- what the fast loop          !
      !      passes the leaf kernel as leaf_env_t%psi, beta_stomata's driver.                         !
      !                                                                                               !
      !      PER-COHORT, not per-patch. A per-patch value would rest on the argument that             !
      !      overnight every cohort relaxes to the same soil -- which is TRUE ONLY IN WET SOIL,       !
      !      and therefore only where this feature does not matter. The wood<->soil relaxation        !
      !      time is tau_w = C_wood/rhizo ~ 9 s at theta 0.25 but ~4.8 DAYS at theta 0.10, so         !
      !      under drought cohorts do not equilibrate overnight at all. Add the gravity head          !
      !      (~0.3 MPa between a 30 m tree and a sapling), differing rooting depth and per-PFT        !
      !      vulnerability, and a patch mean would let a stressed canopy tree hide behind an          !
      !      unstressed sapling. ---------------------------------------------------------------------!
      !      WHY TWO FIELDS AND NOT ONE. They are a DOUBLE BUFFER, and within a day their roles are   !
      !      exact opposites: `dmax_psi_leaf` is READ every fast step by the leaf kernel and written   !
      !      once, at the daily rollover; `dmax_psi_leaf_accum` is WRITTEN every fast step and read    !
      !      once, at that same rollover. A single field cannot do both, because the reset destroys    !
      !      the value the kernel needs:                                                               !
      !                                                                                          !
      !        * reset it at day start  -> nothing left to read for the next 24 h;                     !
      !        * never reset it         -> max() is monotone FOREVER, so it ratchets to the least-      !
      !          negative psi the run has ever seen and never reports stress again. Under drought it    !
      !          would pin at the wettest value in the run's history -- a silent, permanent disabling   !
      !          of the whole closure.                                                                  !
      !                                                                                          !
      !      A merged field also changes the physics: the kernel would read TODAY's partial max, which  !
      !      early in the morning is barely distinguishable from instantaneous psi_leaf -- the tight    !
      !      sub-step feedback loop the design avoids, arriving by accident.                            !
      !                                                                                          !
      !      A genuine ONE-field design exists, but it is a different algorithm, not a merge: SAMPLE    !
      !      psi_leaf once a day at true predawn instead of accumulating a max. That trades this 8      !
      !      bytes/cohort for a clock dependency in the fast loop plus sensitivity to which sub-step    !
      !      is caught and whether the plant has equilibrated. The daily max is clock-free and robust   !
      !      to both, which is why it is what ships. -------------------------------------------------!
      !< [MPa] YESTERDAY's completed daily-max psi_leaf (<= 0); the PUBLISHED value the leaf kernel reads
      real(wp),    allocatable :: dmax_psi_leaf(:)
      real(wp),    allocatable :: dmax_psi_leaf_accum(:) !< [MPa] running max over the CURRENT day; write-only until the rollover
      real(wp),    allocatable :: growth_accum(:)   !< [cm/yr] sum of the samples in the window
      integer(ik), allocatable :: growth_count(:)   !< samples currently in the window (<= growth_window)
      real(wp),    allocatable :: growth_hist(:,:)  !< (growth_window, cohort) ring buffer of growth samples
      !----- DYNAMIC leaf traits, one value per cohort: seeded from the PFT top-of-canopy values     !
      !       at birth and acclimated to light by meds_plant_trait_dynamics; blended on cohort fusion  !
      !       on the quantity each is per (meds_demography_cohort_fusefiss). sla enters the leaf       !
      !       carbon<->area map, leaf_turnover_rate sets baseline leaf turnover, vcmax25/rd25 feed     !
      !       leaf gas exchange. ------------------------------------------------------------------!
      real(wp),    allocatable :: sla(:)            !< [m2/kgC]     specific leaf area
      real(wp),    allocatable :: vcmax25(:)        !< [umol/m2/s]  max carboxylation at 25 degC
      real(wp),    allocatable :: rd25(:)           !< [umol/m2/s]  leaf dark respiration at 25 degC
      real(wp),    allocatable :: leaf_turnover_rate(:) !< [1/yr]   baseline leaf turnover
      !----- PROGNOSTIC fast-biophysics per-cohort state (owned here so it rides the cohort     !
      !      lockstep; mutated by the fast loop, leaf-area-weighted on cohort fusion for the       !
      !      intensive quantities). leaf_temp is a warm-start for the VPD lag. -----------------!
      real(wp),    allocatable :: leaf_temp(:)      !< [K]   cohort leaf temperature
      real(wp),    allocatable :: wood_temp(:)      !< [K]   cohort wood/branch temperature (own energy store)
      !----- Leaf/wood INTERNAL (xylem/symplast) water mass [kg/plant] -- the prognostic             !
      !      hydraulic state; psi is diagnosed from mass via psi_from_water_content wherever         !
      !      needed (docs/science/plant_hydraulics.md §2). EXTENSIVE per plant (like AGB):           !
      !      fusion conserves the nplant-weighted TOTAL, not a leaf-area-weighted mean               !
      !      (fuse_cohort_fast_state). Genuine sub-slow-step hydraulic memory. A cohort is born      !
      !      holding the water its tissues hold at its patch's top-soil potential (init_cohort,      !
      !      patch_birth_state), so 0 is a real amount: a leafless cohort's leaves hold none. --------!
      real(wp),    allocatable :: leaf_water_mass(:) !< [kg/plant] internal leaf water
      real(wp),    allocatable :: wood_water_mass(:) !< [kg/plant] internal wood (sapwood + fine-root) water
      !----- Leaf/wood SURFACE (interception film) water [kg/m2 GROUND] -- a store DISTINCT from     !
      !      the internal water above (docs/science/soil_biophysics.md, "Canopy interception"):      !
      !      rain intercepted onto the canopy, evaporating through the boundary layer (no            !
      !      stomata) rather than through the xylem. Already ground-area-referenced (matches how     !
      !      intercept_canopy_layer's own leaf_water/rain_above/rainfall are all per m2 of the       !
      !      WHOLE patch, not per plant), so fusion SUMS it (like two ground-area contributions      !
      !      adding up), NOT nplant-weights it -- the opposite convention from                       !
      !      leaf_water_mass/wood_water_mass just above. A new cohort's films are dry. ---------------!
      real(wp),    allocatable :: leaf_surf_water(:) !< [kg/m2 ground] leaf interception film
      real(wp),    allocatable :: wood_surf_water(:) !< [kg/m2 ground] wood interception film
      real(wp),    allocatable :: gpp_accum(:)      !< [kgC/plant] GROSS GPP accumulated over the slow step
                                                    !<            (fast->slow carbon bridge; reset each slow step)
      !----- Autotrophic MAINTENANCE respiration accumulated over the slow step (fast->slow bridge, !
      !      reset each slow step; kept SEPARATE from gross gpp_accum, mirroring ED2 today_*_resp).  !
      real(wp),    allocatable :: leaf_resp_accum(:) !< [kgC/plant] leaf dark respiration  (ED2 today_leaf_resp)
      real(wp),    allocatable :: stem_resp_accum(:) !< [kgC/plant] stem maintenance resp  (ED2 today_stem_resp)
      real(wp),    allocatable :: root_resp_accum(:) !< [kgC/plant] fine-root maint. resp  (ED2 today_root_resp)
      !----- The PAR reaching the cohort's top, integrated over the slow step (reset each slow step, !
      !      like gpp_accum): the phenology PAR cue's driver, so each cohort sees its own place in    !
      !      the canopy's light gradient. Consumed by the phenology advance before any fusion. ------!
      real(wp),    allocatable :: light_hours_accum(:) !< [h] time the PAR at the cohort's top exceeded par_min
      !----- PROGNOSTIC leaf phenology (owned here so it rides the cohort lockstep): the two       !
      !      smoothed tendencies and the cue memory, advanced daily by the slow-loop phenology        !
      !      advance; the carbon layer turns the tendencies into leaf growth and loss. The memory     !
      !      must persist across steps -- a running mean or a seasonal sum reset every step would be  !
      !      its own instantaneous input. FUSION POLICY: survivor-keeps (see fuse_cohort_fast_state). !
      real(wp),    allocatable :: leaf_flush_tendency(:)  !< [-]       smoothed flush signal in [0,1]
      real(wp),    allocatable :: leaf_shed_tendency(:)   !< [-]       smoothed senescence signal in [0,1]
      real(wp),    allocatable :: growing_degree_days(:)  !< [K day]   warmth sum since midwinter    (TEMP)
      real(wp),    allocatable :: cold_degree_days(:)     !< [K day]   cold sum since midsummer      (TEMP)
      real(wp),    allocatable :: dry_psi_sum(:)          !< [MPa day] predawn psi below the TLP     (WATER)
      real(wp),    allocatable :: wet_psi_sum(:)          !< [MPa day] predawn psi above the TLP     (WATER)
      real(wp),    allocatable :: light_hours_mean(:)     !< [h/day] running-mean hours of light (LIGHT); < 0: none yet
      !----- Host-only back-index used to regroup the flat array by patch. ----------------!
      integer(ik), allocatable :: owner_patch(:)
      !----- Persistent identity: a global id stamped at creation and carried (in lockstep   !
      !      with every other field) through sorts/fusion/compaction, so an external tracker  !
      !      can follow one cohort across output records until it fuses away or is culled.     !
      integer(ik), allocatable :: global_id(:)
      !----- TRANSIENT fast-loop DIAGNOSTIC accumulators (meds_diag_vars). Not prognostic          !
      !      state and never restarted, but they must ride this lockstep because restructuring       !
      !      happens INSIDE the slow step, between the fast loop filling them and the output tick     !
      !      reading them. Stored as one 2-D array so every permutation below is a single statement   !
      !      that cannot omit a field. Inactive (and unallocated) unless the run reports per-cohort    !
      !      ecophysiology.  --------------------------------------------------------------------!
      type(cohort_diag_block) :: diag
      !----- SLOW-loop tendency diagnostics (growth / mortality / NPP allocation), same block type   !
      !      with its own field set. Separate from `diag` because it takes ONE sample per slow step   !
      !      where that takes ~48, so one shared dt weight could not normalize both.  ---------------!
      type(cohort_diag_block) :: sdiag
   end type cohort_block

   !----- TRANSIENT per-cohort TIME-DERIVATIVE bundle (all rates [unit/yr]). The slow-loop driver !
   !      COMPUTES these (carbon or empirical) each step; update_cohort_states then advances the   !
   !      cohort state by them. It is index-aligned with cohort_block but is NOT persistent state:  !
   !      filled the same step it is consumed (before any sort/fuse) and NOT lockstep-reordered.    !
   !      dln_nplant_dt is the LOG-space mortality rate (nplant *= exp(dln*dt)); the carbon pools'  !
   !      are additive rates (state += rate*dt). Only PROGNOSTIC state has a tendency: the geometry  !
   !      is derived from the carbon after the update (update_cohort_states), never integrated.     !
   type :: cohort_deriv_block
      real(wp), allocatable :: d_leaf_carbon_dt(:), d_fineroot_carbon_dt(:)
      real(wp), allocatable :: d_wood_carbon_dt(:), d_nonstructural_carbon_dt(:)
      real(wp), allocatable :: dln_nplant_dt(:)
   end type cohort_deriv_block

   !----- The bounded cohort fields cohort_bad_field checks, in its order. ------------------!
   integer(ik), parameter :: N_COHORT_BOUNDED = 10_ik
   character(len=20), parameter :: COHORT_BOUNDED_NAME(N_COHORT_BOUNDED) = [character(len=20) ::   &
      'nplant', 'leaf_carbon', 'fineroot_carbon', 'wood_carbon', 'nonstructural_carbon', 'sla',     &
      'leaf_water_mass', 'wood_water_mass', 'leaf_surf_water', 'wood_surf_water']

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
      ok = allocated(deriv%d_leaf_carbon_dt)       .and. allocated(deriv%d_fineroot_carbon_dt) .and. &
           allocated(deriv%d_wood_carbon_dt)       .and. allocated(deriv%d_nonstructural_carbon_dt) .and. &
           allocated(deriv%dln_nplant_dt)
      if (ok) then
         if (size(deriv%d_leaf_carbon_dt) >= m) return   ! big enough (all co-sized) -> reuse
         deallocate(deriv%d_leaf_carbon_dt, deriv%d_fineroot_carbon_dt, deriv%d_wood_carbon_dt,      &
                    deriv%d_nonstructural_carbon_dt, deriv%dln_nplant_dt)
      else                                       ! none or partial -> drop whatever is there
         if (allocated(deriv%d_leaf_carbon_dt))          deallocate(deriv%d_leaf_carbon_dt)
         if (allocated(deriv%d_fineroot_carbon_dt))      deallocate(deriv%d_fineroot_carbon_dt)
         if (allocated(deriv%d_wood_carbon_dt))          deallocate(deriv%d_wood_carbon_dt)
         if (allocated(deriv%d_nonstructural_carbon_dt)) deallocate(deriv%d_nonstructural_carbon_dt)
         if (allocated(deriv%dln_nplant_dt))             deallocate(deriv%dln_nplant_dt)
      end if
      allocate(deriv%d_leaf_carbon_dt(m), deriv%d_fineroot_carbon_dt(m), deriv%d_wood_carbon_dt(m),  &
               deriv%d_nonstructural_carbon_dt(m), deriv%dln_nplant_dt(m))
   end subroutine cohort_deriv_alloc

   subroutine cohort_alloc(cohort, cap, nwin)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: cap, nwin
      cohort%cap = cap ; cohort%n = 0_ik
      allocate(cohort%pft(cap), cohort%owner_patch(cap), cohort%global_id(cap))
      allocate(cohort%nplant(cap), cohort%dbh(cap), cohort%height(cap), cohort%basal_area(cap),            &
               cohort%agb(cap), cohort%leaf_area(cap), cohort%overtopping_lai(cap), cohort%growth_avg(cap),&
               cohort%dmax_psi_leaf(cap), cohort%dmax_psi_leaf_accum(cap),                        &
               cohort%growth_accum(cap), cohort%growth_count(cap), cohort%growth_hist(nwin, cap))
      allocate(cohort%leaf_carbon(cap), cohort%fineroot_carbon(cap), cohort%wood_carbon(cap),        &
               cohort%nonstructural_carbon(cap))
      allocate(cohort%sla(cap))
      allocate(cohort%wood_area(cap), cohort%sapwood_carbon(cap), cohort%sapwood_area(cap))
      allocate(cohort%vcmax25(cap), cohort%rd25(cap), cohort%leaf_turnover_rate(cap))
      allocate(cohort%leaf_temp(cap), cohort%wood_temp(cap), cohort%gpp_accum(cap), cohort%light_hours_accum(cap))
      allocate(cohort%leaf_water_mass(cap), cohort%wood_water_mass(cap))
      allocate(cohort%leaf_surf_water(cap), cohort%wood_surf_water(cap))
      allocate(cohort%leaf_resp_accum(cap), cohort%stem_resp_accum(cap), cohort%root_resp_accum(cap))
      allocate(cohort%leaf_flush_tendency(cap), cohort%leaf_shed_tendency(cap),                       &
               cohort%growing_degree_days(cap), cohort%cold_degree_days(cap),                                   &
               cohort%dry_psi_sum(cap),                                                   &
               cohort%wet_psi_sum(cap), cohort%light_hours_mean(cap))
      cohort%leaf_temp = 0.0_wp ; cohort%wood_temp = 0.0_wp      ! a slot's tissue state is set at birth
      cohort%leaf_water_mass = 0.0_wp ; cohort%wood_water_mass = 0.0_wp ; cohort%gpp_accum = 0.0_wp
      cohort%light_hours_accum = 0.0_wp
      cohort%leaf_surf_water = 0.0_wp ; cohort%wood_surf_water = 0.0_wp
      cohort%leaf_resp_accum = 0.0_wp ; cohort%stem_resp_accum = 0.0_wp ; cohort%root_resp_accum = 0.0_wp
      cohort%leaf_flush_tendency = PHENO_FLUSH_INIT ; cohort%leaf_shed_tendency = PHENO_SHED_INIT
      cohort%growing_degree_days = 0.0_wp ; cohort%cold_degree_days = 0.0_wp
      cohort%dry_psi_sum = 0.0_wp ; cohort%wet_psi_sum = 0.0_wp ; cohort%light_hours_mean = -1.0_wp
      cohort%pft = 0_ik ; cohort%owner_patch = 0_ik ; cohort%global_id = 0_ik
      cohort%nplant = 0.0_wp ; cohort%dbh = 0.0_wp ; cohort%height = 0.0_wp ; cohort%basal_area = 0.0_wp
      cohort%agb = 0.0_wp ; cohort%leaf_area = 0.0_wp ; cohort%overtopping_lai = 0.0_wp
      cohort%growth_avg = GROWTH_AVG_UNSET
      cohort%dmax_psi_leaf = DMAX_PSI_LEAF_UNSET ; cohort%dmax_psi_leaf_accum = DMAX_PSI_LEAF_ACCUM_RESET
      cohort%growth_accum = 0.0_wp ; cohort%growth_count = 0_ik ; cohort%growth_hist = 0.0_wp
      cohort%leaf_carbon = 0.0_wp ; cohort%fineroot_carbon = 0.0_wp ; cohort%wood_carbon = 0.0_wp
      cohort%nonstructural_carbon = 0.0_wp
      cohort%sla = 0.0_wp
      cohort%wood_area = 0.0_wp ; cohort%sapwood_carbon = 0.0_wp ; cohort%sapwood_area = 0.0_wp
      cohort%vcmax25 = 0.0_wp ; cohort%rd25 = 0.0_wp ; cohort%leaf_turnover_rate = 0.0_wp
   end subroutine cohort_alloc

   !----- Free every per-cohort array. ------------------------------------------------------!
   subroutine cohort_free(cohort)
      type(cohort_block), intent(inout) :: cohort
      cohort%n = 0_ik ; cohort%cap = 0_ik
      if (allocated(cohort%nplant)) deallocate(cohort%pft, cohort%nplant, cohort%dbh, cohort%height,     &
         cohort%basal_area, cohort%agb, cohort%leaf_area, cohort%growth_avg, cohort%growth_accum,        &
         cohort%growth_count, cohort%dmax_psi_leaf, cohort%dmax_psi_leaf_accum, cohort%growth_hist,      &
         cohort%sla, cohort%vcmax25, cohort%rd25, cohort%leaf_turnover_rate, cohort%wood_area,           &
         cohort%sapwood_carbon, cohort%sapwood_area, cohort%leaf_carbon, cohort%fineroot_carbon,         &
         cohort%wood_carbon, cohort%nonstructural_carbon, cohort%owner_patch, cohort%global_id,          &
         cohort%overtopping_lai, cohort%leaf_temp, cohort%wood_temp, cohort%leaf_water_mass,             &
         cohort%wood_water_mass, cohort%leaf_surf_water, cohort%wood_surf_water, cohort%gpp_accum,       &
         cohort%light_hours_accum, cohort%leaf_resp_accum, cohort%stem_resp_accum,                       &
         cohort%root_resp_accum, cohort%leaf_flush_tendency, cohort%leaf_shed_tendency,                  &
         cohort%growing_degree_days, cohort%cold_degree_days, cohort%dry_psi_sum, cohort%wet_psi_sum,    &
         cohort%light_hours_mean)
   end subroutine cohort_free

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
      !----- Grow the diagnostic block IN PLACE first, so move_alloc_block below can hand its      !
      !      (already correctly sized, prefix-preserving) arrays straight across.  ---------------!
      call cohort_diag_grow(cohort%diag,  newcap)
      call cohort_diag_grow(cohort%sdiag, newcap)
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
      tmp%dmax_psi_leaf(1:m)   = cohort%dmax_psi_leaf(1:m)
      tmp%dmax_psi_leaf_accum(1:m) = cohort%dmax_psi_leaf_accum(1:m)
      tmp%growth_accum(1:m)   = cohort%growth_accum(1:m)
      tmp%growth_count(1:m)   = cohort%growth_count(1:m)
      tmp%growth_hist(:,1:m)  = cohort%growth_hist(:,1:m)
      tmp%leaf_carbon(1:m)           = cohort%leaf_carbon(1:m)
      tmp%fineroot_carbon(1:m)       = cohort%fineroot_carbon(1:m)
      tmp%wood_carbon(1:m)           = cohort%wood_carbon(1:m)
      tmp%nonstructural_carbon(1:m)  = cohort%nonstructural_carbon(1:m)
      tmp%sla(1:m)                 = cohort%sla(1:m)
      tmp%vcmax25(1:m)               = cohort%vcmax25(1:m)
      tmp%rd25(1:m)                  = cohort%rd25(1:m)
      tmp%leaf_turnover_rate(1:m)    = cohort%leaf_turnover_rate(1:m)
      tmp%wood_area(1:m) = cohort%wood_area(1:m)
      tmp%sapwood_carbon(1:m) = cohort%sapwood_carbon(1:m)
      tmp%sapwood_area(1:m) = cohort%sapwood_area(1:m)
      tmp%global_id(1:m)      = cohort%global_id(1:m)
      tmp%leaf_temp(1:m)      = cohort%leaf_temp(1:m)
      tmp%wood_temp(1:m)      = cohort%wood_temp(1:m)
      tmp%leaf_water_mass(1:m) = cohort%leaf_water_mass(1:m)
      tmp%wood_water_mass(1:m) = cohort%wood_water_mass(1:m)
      tmp%leaf_surf_water(1:m) = cohort%leaf_surf_water(1:m)
      tmp%wood_surf_water(1:m) = cohort%wood_surf_water(1:m)
      tmp%gpp_accum(1:m)      = cohort%gpp_accum(1:m)
      tmp%light_hours_accum(1:m)      = cohort%light_hours_accum(1:m)
      tmp%leaf_resp_accum(1:m) = cohort%leaf_resp_accum(1:m)
      tmp%stem_resp_accum(1:m) = cohort%stem_resp_accum(1:m)
      tmp%root_resp_accum(1:m) = cohort%root_resp_accum(1:m)
      tmp%leaf_flush_tendency(1:m) = cohort%leaf_flush_tendency(1:m)
      tmp%leaf_shed_tendency(1:m)  = cohort%leaf_shed_tendency(1:m)
      tmp%growing_degree_days(1:m) = cohort%growing_degree_days(1:m)
      tmp%cold_degree_days(1:m)    = cohort%cold_degree_days(1:m)
      tmp%dry_psi_sum(1:m)         = cohort%dry_psi_sum(1:m)
      tmp%wet_psi_sum(1:m)         = cohort%wet_psi_sum(1:m)
      tmp%light_hours_mean(1:m)      = cohort%light_hours_mean(1:m)
      !----- Carry the (already grown) diagnostic block over the fresh tmp, whose own diag is      !
      !      inactive -- move_alloc_block copies tmp INTO cohort, so it must be the one holding it. !
      tmp%diag = cohort%diag ; tmp%sdiag = cohort%sdiag
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
      call move_alloc(src%dmax_psi_leaf, dst%dmax_psi_leaf)
      call move_alloc(src%dmax_psi_leaf_accum, dst%dmax_psi_leaf_accum)
      call move_alloc(src%growth_accum, dst%growth_accum)
      call move_alloc(src%growth_count, dst%growth_count)
      call move_alloc(src%growth_hist, dst%growth_hist)
      call move_alloc(src%leaf_carbon, dst%leaf_carbon)
      call move_alloc(src%fineroot_carbon, dst%fineroot_carbon)
      call move_alloc(src%wood_carbon, dst%wood_carbon)
      call move_alloc(src%nonstructural_carbon, dst%nonstructural_carbon)
      call move_alloc(src%sla, dst%sla)
      call move_alloc(src%vcmax25, dst%vcmax25)
      call move_alloc(src%rd25, dst%rd25)
      call move_alloc(src%leaf_turnover_rate, dst%leaf_turnover_rate)
      call move_alloc(src%wood_area, dst%wood_area)
      call move_alloc(src%sapwood_carbon, dst%sapwood_carbon)
      call move_alloc(src%sapwood_area, dst%sapwood_area)
      call move_alloc(src%global_id, dst%global_id)
      call move_alloc(src%leaf_temp, dst%leaf_temp)
      call move_alloc(src%wood_temp, dst%wood_temp)
      call move_alloc(src%leaf_water_mass, dst%leaf_water_mass)
      call move_alloc(src%wood_water_mass, dst%wood_water_mass)
      call move_alloc(src%leaf_surf_water, dst%leaf_surf_water)
      call move_alloc(src%wood_surf_water, dst%wood_surf_water)
      call move_alloc(src%gpp_accum, dst%gpp_accum)
      call move_alloc(src%light_hours_accum, dst%light_hours_accum)
      call move_alloc(src%leaf_resp_accum, dst%leaf_resp_accum)
      call move_alloc(src%stem_resp_accum, dst%stem_resp_accum)
      call move_alloc(src%root_resp_accum, dst%root_resp_accum)
      call move_alloc(src%leaf_flush_tendency, dst%leaf_flush_tendency)
      call move_alloc(src%leaf_shed_tendency, dst%leaf_shed_tendency)
      call move_alloc(src%growing_degree_days, dst%growing_degree_days)
      call move_alloc(src%cold_degree_days, dst%cold_degree_days)
      call move_alloc(src%dry_psi_sum, dst%dry_psi_sum)
      call move_alloc(src%wet_psi_sum, dst%wet_psi_sum)
      call move_alloc(src%light_hours_mean, dst%light_hours_mean)
   end subroutine move_alloc_block

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
      cohort%dmax_psi_leaf(1:m)   = cohort%dmax_psi_leaf(perm(1:m))
      cohort%dmax_psi_leaf_accum(1:m) = cohort%dmax_psi_leaf_accum(perm(1:m))
      cohort%growth_accum(1:m)   = cohort%growth_accum(perm(1:m))
      cohort%growth_count(1:m)   = cohort%growth_count(perm(1:m))
      cohort%growth_hist(:,1:m)  = cohort%growth_hist(:,perm(1:m))
      cohort%leaf_carbon(1:m)           = cohort%leaf_carbon(perm(1:m))
      cohort%fineroot_carbon(1:m)       = cohort%fineroot_carbon(perm(1:m))
      cohort%wood_carbon(1:m)           = cohort%wood_carbon(perm(1:m))
      cohort%nonstructural_carbon(1:m)  = cohort%nonstructural_carbon(perm(1:m))
      cohort%sla(1:m)                 = cohort%sla(perm(1:m))
      cohort%vcmax25(1:m)               = cohort%vcmax25(perm(1:m))
      cohort%rd25(1:m)                  = cohort%rd25(perm(1:m))
      cohort%leaf_turnover_rate(1:m)    = cohort%leaf_turnover_rate(perm(1:m))
      cohort%wood_area(1:m) = cohort%wood_area(perm(1:m))
      cohort%sapwood_carbon(1:m) = cohort%sapwood_carbon(perm(1:m))
      cohort%sapwood_area(1:m) = cohort%sapwood_area(perm(1:m))
      cohort%global_id(1:m)      = cohort%global_id(perm(1:m))
      cohort%leaf_temp(1:m)      = cohort%leaf_temp(perm(1:m))
      cohort%wood_temp(1:m)      = cohort%wood_temp(perm(1:m))
      cohort%leaf_water_mass(1:m) = cohort%leaf_water_mass(perm(1:m))
      cohort%wood_water_mass(1:m) = cohort%wood_water_mass(perm(1:m))
      cohort%leaf_surf_water(1:m) = cohort%leaf_surf_water(perm(1:m))
      cohort%wood_surf_water(1:m) = cohort%wood_surf_water(perm(1:m))
      cohort%gpp_accum(1:m)      = cohort%gpp_accum(perm(1:m))
      cohort%light_hours_accum(1:m)      = cohort%light_hours_accum(perm(1:m))
      cohort%leaf_resp_accum(1:m) = cohort%leaf_resp_accum(perm(1:m))
      cohort%stem_resp_accum(1:m) = cohort%stem_resp_accum(perm(1:m))
      cohort%root_resp_accum(1:m) = cohort%root_resp_accum(perm(1:m))
      cohort%leaf_flush_tendency(1:m) = cohort%leaf_flush_tendency(perm(1:m))
      cohort%leaf_shed_tendency(1:m)  = cohort%leaf_shed_tendency(perm(1:m))
      cohort%growing_degree_days(1:m) = cohort%growing_degree_days(perm(1:m))
      cohort%cold_degree_days(1:m)    = cohort%cold_degree_days(perm(1:m))
      cohort%dry_psi_sum(1:m)         = cohort%dry_psi_sum(perm(1:m))
      cohort%wet_psi_sum(1:m)         = cohort%wet_psi_sum(perm(1:m))
      cohort%light_hours_mean(1:m)      = cohort%light_hours_mean(perm(1:m))
      !----- The diagnostic accumulators ride the SAME permutation. One call, and it cannot omit  !
      !      a field: they are rows of one 2-D array (meds_diag_vars, decision 4).  --------      !
      call cohort_diag_reorder(cohort%diag,  perm, m)
      call cohort_diag_reorder(cohort%sdiag, perm, m)
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
      cohort%dmax_psi_leaf(dst)   = cohort%dmax_psi_leaf(src)
      cohort%dmax_psi_leaf_accum(dst) = cohort%dmax_psi_leaf_accum(src)
      cohort%growth_accum(dst)   = cohort%growth_accum(src)
      cohort%growth_count(dst)   = cohort%growth_count(src)
      cohort%growth_hist(:,dst)  = cohort%growth_hist(:,src)
      cohort%leaf_carbon(dst)           = cohort%leaf_carbon(src)
      cohort%fineroot_carbon(dst)       = cohort%fineroot_carbon(src)
      cohort%wood_carbon(dst)           = cohort%wood_carbon(src)
      cohort%nonstructural_carbon(dst)  = cohort%nonstructural_carbon(src)
      cohort%sla(dst)                 = cohort%sla(src)
      cohort%vcmax25(dst)               = cohort%vcmax25(src)
      cohort%rd25(dst)                  = cohort%rd25(src)
      cohort%leaf_turnover_rate(dst)    = cohort%leaf_turnover_rate(src)
      cohort%wood_area(dst) = cohort%wood_area(src)
      cohort%sapwood_carbon(dst) = cohort%sapwood_carbon(src)
      cohort%sapwood_area(dst) = cohort%sapwood_area(src)
      cohort%global_id(dst)      = cohort%global_id(src)
      call cohort_diag_copy_slot(cohort%diag,  dst, src)
      call cohort_diag_copy_slot(cohort%sdiag, dst, src)
      cohort%leaf_temp(dst)      = cohort%leaf_temp(src)
      cohort%wood_temp(dst)      = cohort%wood_temp(src)
      cohort%leaf_water_mass(dst) = cohort%leaf_water_mass(src)
      cohort%wood_water_mass(dst) = cohort%wood_water_mass(src)
      cohort%leaf_surf_water(dst) = cohort%leaf_surf_water(src)
      cohort%wood_surf_water(dst) = cohort%wood_surf_water(src)
      cohort%gpp_accum(dst)      = cohort%gpp_accum(src)
      cohort%light_hours_accum(dst)      = cohort%light_hours_accum(src)
      cohort%leaf_resp_accum(dst) = cohort%leaf_resp_accum(src)
      cohort%stem_resp_accum(dst) = cohort%stem_resp_accum(src)
      cohort%root_resp_accum(dst) = cohort%root_resp_accum(src)
      cohort%leaf_flush_tendency(dst) = cohort%leaf_flush_tendency(src)
      cohort%leaf_shed_tendency(dst)  = cohort%leaf_shed_tendency(src)
      cohort%growing_degree_days(dst) = cohort%growing_degree_days(src)
      cohort%cold_degree_days(dst)    = cohort%cold_degree_days(src)
      cohort%dry_psi_sum(dst)         = cohort%dry_psi_sum(src)
      cohort%wet_psi_sum(dst)         = cohort%wet_psi_sum(src)
      cohort%light_hours_mean(dst)      = cohort%light_hours_mean(src)
   end subroutine copy_cohort_slot

   !---------------------------------------------------------------------------------------!
   ! Rescale the per-cohort fields that are referenced to the PATCH GROUND AREA [kg/m2     !
   ! ground] rather than to the plant -- the canopy interception films. Every operation    !
   ! that changes the ground area a cohort's density refers to must apply the same factor  !
   ! to these fields as it applies to nplant -- cohort fission (0.5 each daughter), patch  !
   ! fusion (area weights), and the survivor copy into a disturbance gap. Copying them     !
   ! verbatim after halving nplant would double the site's film water on every split.      !
   !---------------------------------------------------------------------------------------!
   pure subroutine scale_cohort_ground_fields(cohort, i, factor)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: i
      real(wp),           intent(in)    :: factor
      !----- GROUND-referenced fields only -- see fuse_cohort_fast_state for the one place the      !
      !      per-field policy is declared. A field that is per PLANT must NOT appear here: it       !
      !      carries no area normalization to rescale.  ------------------------------------------!
      cohort%leaf_surf_water(i) = cohort%leaf_surf_water(i) * factor
      cohort%wood_surf_water(i) = cohort%wood_surf_water(i) * factor
   end subroutine scale_cohort_ground_fields

   !=========================================================================================!
   !  THE PER-FIELD FUSION POLICY OF THE FAST-LOOP-OWNED COHORT STATE -- declared once,      !
   !  here.                                                                                  !
   !                                                                                         !
   !  These fields are written by the fast loop and read back by it; the slow loop only has  !
   !  to carry them correctly through a fusion. Each is one of three kinds, and getting the  !
   !  kind wrong is invisible: the AGB assert still passes, the conservation ledgers still   !
   !  close, and the answer is quietly wrong.                                                !
   !                                                                                         !
   !    FK_INTENSIVE  a per-mass property (a temperature). Weighted by the HEAT CAPACITY     !
   !                  that carries it, which is what conserves the tissue's energy through   !
   !                  the merge.                                                             !
   !    FK_EXTENSIVE  a per-PLANT amount (tissue water, an accumulated per-plant flux).      !
   !                  nplant-weighted mean, so the site total (nplant * value, summed) is    !
   !                  conserved.                                                             !
   !    FK_GROUND     ALREADY per m2 of patch ground (the interception film). Two cohorts    !
   !                  over the same ground simply ADD; weighting would double-count the      !
   !                  normalization.                                                         !
   !                                                                                         !
   !  The diagnostic twins are fused from their own declared table (CDIAG_FUSE,              !
   !  meds_diag_vars) on weights the caller computes: the tissue heat capacities for         !
   !  the leaf and wood temperatures, as here, and leaf area for a per-LEAF-AREA diagnostic  !
   !  (transpiration per m2 of leaf, say), which is not a temperature. Adding a fast-loop    !
   !  field is one line here; the compiler will not remind you, so the line is the reminder. !
   !=========================================================================================!
   !----- THE per-field policy for the fast-loop cohort state, declared once. Three kinds, and    !
   !      the headings below are the declaration: INTENSIVE (weighted by the capacity that carries  !
   !      the quantity), EXTENSIVE (nplant-weighted), GROUND (added, never weighted).                !
   !                                                                                          !
   !      It is centralised but not enforced: a per-cohort fast field left out of this routine       !
   !      compiles clean and is silently wrong at the next fusion. A new field needs a line here, a   !
   !      line in scale_cohort_ground_fields if it is per m2 of ground, and a witness in             !
   !      test_fusion_cohort. A table the blend iterated would be no safer (the compiler cannot      !
   !      check that a table lists every field either) and would be longer than these lines. -------!
   pure subroutine fuse_cohort_fast_state(cohort, recc, donc, np_r, np_d)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: recc, donc   !< survivor, donor
      real(wp),           intent(in)    :: np_r, np_d   !< nplant weights
      real(wp) :: clr, cwr, cld, cwd
      !----- INTENSIVE: tissue temperatures, weighted by the HEAT CAPACITY that carries them.    !
      !                                                                                          !
      !      Leaf-area weighting would be roughly right for leaves, since leaf heat capacity         !
      !      tracks leaf carbon and leaf area tracks it through sla. For WOOD it is not even         !
      !      approximately right: wood heat capacity follows wood carbon and the sapwood ring,       !
      !      which have no reason to scale with leaf area -- so fusing two cohorts of different      !
      !      wood:leaf ratio on leaf area would move the tissue energy by ~1e4 J an event.           !
      !                                                                                          !
      !      Capacity weighting conserves `cap x T` EXACTLY for the additive part, because the    !
      !      merged capacity IS the sum: leaf_carbon, wood_carbon and the tissue waters are all   !
      !      nplant-weighted below, so nplant x pool adds across the merge. What does not add is  !
      !      the sapwood ring, which set_cohort_size_from_carbon re-derives from the merged dbh;  !
      !      that leaves a small thermal-mass change, of the same kind growth makes and declared   !
      !      the same way by the caller.  --------------------------------------------------!
      call cohort_tissue_heat_capacity(cohort, recc, clr, cwr)
      call cohort_tissue_heat_capacity(cohort, donc, cld, cwd)
      call blend(cohort%leaf_temp,        clr, cld)
      call blend(cohort%wood_temp,        cwr, cwd)
      !----- EXTENSIVE: internal tissue water [kg/plant] and the per-plant flux accumulators. --!
      call blend(cohort%leaf_water_mass,  np_r, np_d)
      call blend(cohort%wood_water_mass,  np_r, np_d)
      call blend(cohort%gpp_accum,        np_r, np_d)
      call blend(cohort%leaf_resp_accum,  np_r, np_d)
      call blend(cohort%stem_resp_accum,  np_r, np_d)
      call blend(cohort%root_resp_accum,  np_r, np_d)
      !----- GROUND: the interception film [kg/m2 ground] -- add, never weight. ----------------!
      cohort%leaf_surf_water(recc) = cohort%leaf_surf_water(recc) + cohort%leaf_surf_water(donc)
      cohort%wood_surf_water(recc) = cohort%wood_surf_water(recc) + cohort%wood_surf_water(donc)
      !----- SURVIVOR-KEEPS, declared rather than left implicit: the phenology memories            !
      !      (the two tendencies, the degree-day and psi sums, the light mean), light_hours_accum    !
      !      (already consumed this step) and dmax_psi_leaf are absent from this routine ON PURPOSE,  !
      !      so the survivor keeps its own. Two cohorts only fuse when they already match in size,   !
      !      height and PFT, and every one of these is driven by a signal both of them saw (the PAR at a !
      !      shared height included), so the two memories are near-identical at fusion and a blend would !
      !      move nothing. Blending is also not obviously right for a sum that resets on a trigger:  !
      !      an nplant-weighted average of a sum that has just reset and one that has not is a       !
      !      history neither cohort had. If a future cue is driven by something genuinely           !
      !      per-cohort and slow, revisit this -- and add it here rather than leaving it implicit.   !
   contains
      pure subroutine blend(a, wr, wd)
         real(wp), intent(inout) :: a(:)
         real(wp), intent(in)    :: wr, wd
         a(recc) = pair_mean(a(recc), wr, a(donc), wd)
      end subroutine blend
   end subroutine fuse_cohort_fast_state

   !---------------------------------------------------------------------------------------!
   ! The cached WOOD geometry of one slot, from the just-set dbh / basal_area / wood_carbon.  !
   ! Part of set_cohort_size_from_carbon, so the leaf and wood caches can never disagree about !
   ! which dbh they were derived from. Everything here is PER PLANT (see the declarations):   !
   ! `wood_area` is the per-plant twin of `leaf_area`, so WAI = nplant*wood_area at the point  !
   ! of use and no cached field carries a stale plant density.                                 !
   !---------------------------------------------------------------------------------------!
   pure subroutine set_cohort_wood_geometry(cohort, i)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: i
      real(wp)    :: f_sap
      integer(ik) :: pf
      pf = cohort%pft(i)
      associate (pft => cohort%pft_table)
         !----- dbh_to_wai is exactly linear in nplant, so passing 1 gives the per-plant area. ---!
         cohort%wood_area(i)      = dbh_to_wai(cohort%dbh(i), 1.0_wp, pft%wai_b1(pf), pft%wai_b2(pf))
         f_sap                    = sapwood_fraction(cohort%dbh(i), pft%sapwood_area_b1(pf),          &
                                                     pft%sapwood_area_b2(pf))
         cohort%sapwood_carbon(i) = f_sap * cohort%wood_carbon(i)
         cohort%sapwood_area(i)   = f_sap * cohort%basal_area(i)
      end associate
   end subroutine set_cohort_wood_geometry

   !---------------------------------------------------------------------------------------!
   ! Place ONE cohort slot on the allometry from its dbh (setup, census, recruitment, a state  !
   ! file without pools): the four carbon pools a plant of that dbh holds, then the geometry   !
   ! derived from those pools by set_cohort_size_from_carbon, exactly as after every other     !
   ! carbon change. Geometry has that one writer, so a cohort's geometry is always the same    !
   ! function of its carbon and a restart that re-derives it gets the same bits.               !
   !---------------------------------------------------------------------------------------!
   subroutine set_cohort_size(cohort, i)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: i
      real(wp)    :: height
      integer(ik) :: pf
      pf = cohort%pft(i)
      associate (pft => cohort%pft_table)
         height = dbh_to_height(cohort%allom, cohort%dbh(i), pft%hgt_max(pf))
         cohort%leaf_carbon(i)          = size2leaf_carbon(cohort%allom, cohort%dbh(i), height, cohort%sla(i))
         cohort%wood_carbon(i)          = size2wood_carbon(cohort%allom, cohort%dbh(i), height,              &
                                                           pft%wood_density(pf), pft%aboveground_frac(pf))
         cohort%fineroot_carbon(i)      = pft%root_to_leaf_ratio(pf) * cohort%leaf_carbon(i)
         cohort%nonstructural_carbon(i) = storage_carbon_target(pft%storage_cushion(pf), cohort%leaf_carbon(i), &
                                                                cohort%fineroot_carbon(i))
      end associate
      call set_cohort_size_from_carbon(cohort, i)
   end subroutine set_cohort_size

   !=========================================================================================!
   ! cohort_tissue_heat_capacity -- one cohort's leaf and wood heat capacity [J/m2           !
   ! ground/K].                                                                              !
   !                                                                                         !
   ! WHY IT LIVES HERE. Three consumers need the same number and none of them may link the   !
   ! others: the slow ledger values the tissue energy store with it, the demography          !
   ! operators need it to say how much heat a dying cohort takes with it, and the fast       !
   ! loop's frozen record builds the same quantity to relax the tissue against. This is the  !
   ! state layer, which all three already import, so it is the one place all three can       !
   ! share.                                                                                  !
   !                                                                                         !
   ! It reproduces meds_fast_frozen's construction exactly                                   !
   ! (docs/science/vegetation_energy_dynamics.md eq. (8)), including two choices worth       !
   ! naming. The WOOD's dry tissue is all the wood, at the dry specific heat, and its water  !
   ! is the sapwood share of the hydraulic store `wood_water_mass` (wood_thermal_water, one  !
   ! function for both): the fine roots' share sits in the soil and heartwood holds none.    !
   ! Every term scales with plant number, so splitting a cohort in two changes neither the   !
   ! sum nor a fusion weight. The fast loop builds its own copy from the column view; this   !
   ! function is the half of that pair that three libraries can reach.                      !
   !=========================================================================================!
   pure subroutine cohort_tissue_heat_capacity(cohort, i, cap_leaf, cap_wood)
      type(cohort_block), intent(in)  :: cohort
      integer(ik),        intent(in)  :: i
      real(wp),           intent(out) :: cap_leaf, cap_wood             !< [J/m2 ground/K]
      real(wp) :: n
      n = cohort%nplant(i)
      cap_leaf = cohort%leaf_carbon(i) * n * C2B_WOOD * TISSUE_C_LEAF                                &
               + max(cohort%leaf_water_mass(i), 0.0_wp) * n * cp_liq
      cap_wood = cohort%wood_carbon(i) * n * C2B_WOOD * TISSUE_C_WOOD_DRY                            &
               + n * wood_thermal_water(cohort%wood_water_mass(i), cohort%sapwood_carbon(i),              &
                                        cohort%fineroot_carbon(i)) * cp_liq
   end subroutine cohort_tissue_heat_capacity

   !----- One cohort's INTERNAL (hydraulic) tissue water [kg/m2 ground]: the stores the plant     !
   !       hydraulics carries, of which the heat capacity above counts the leaf's and the wood     !
   !       store's sapwood share. -----------------------------------------------------------------!
   pure function cohort_tissue_water(cohort, i) result(w)
      type(cohort_block), intent(in) :: cohort
      integer(ik),        intent(in) :: i
      real(wp)                       :: w
      w = cohort%nplant(i) * (max(cohort%leaf_water_mass(i), 0.0_wp)                               &
                            + max(cohort%wood_water_mass(i), 0.0_wp))
   end function cohort_tissue_water

   !---------------------------------------------------------------------------------------!
   ! The bounds of a cohort's state: its plant density, carbon pools, SLA, and the water it    !
   ! holds (internal stores and the films on its leaves and wood) are amounts, never below     !
   ! zero and never infinite. Every operator keeps them in bounds by construction -- the        !
   ! applier zeroes only a pool's round-off, mortality never takes nplant below `negligible`,   !
   ! the fast loop clamps the films into [0, capacity] -- so a violation is an operator bug.    !
   ! The driver checks this every step and names the field, so the bug is caught the step it    !
   ! happens instead of through a quantity derived from it, somewhere else, years later. The    !
   ! geometry is not checked: it is a function of these. Returns the first bad field's index,   !
   ! 0 if none.                                                                                 !
   !---------------------------------------------------------------------------------------!
   pure function cohort_bad_field(cohort, i) result(k)
      type(cohort_block), intent(in) :: cohort
      integer(ik),        intent(in) :: i
      integer(ik)                    :: k
      real(wp)    :: x(N_COHORT_BOUNDED)
      integer(ik) :: j
      x = [cohort%nplant(i), cohort%leaf_carbon(i), cohort%fineroot_carbon(i),                      &
           cohort%wood_carbon(i), cohort%nonstructural_carbon(i), cohort%sla(i),                    &
           cohort%leaf_water_mass(i), cohort%wood_water_mass(i),                                    &
           cohort%leaf_surf_water(i), cohort%wood_surf_water(i)]
      k = 0_ik
      do j = 1_ik, N_COHORT_BOUNDED
         !----- NaN first, by classification, as soil_carbon_bad_pool does: a comparison with a NaN !
         !      raises IEEE invalid, which a trapping build stops on. ---------------------------!
         if (ieee_is_nan(x(j))) then
            k = j
            return
         end if
         if (x(j) < 0.0_wp .or. x(j) > huge(x(j))) then
            k = j
            return
         end if
      end do
   end function cohort_bad_field

   !----- Name for the index cohort_bad_field returns, so a failure message is readable. -------!
   pure function cohort_field_name(k) result(nm)
      integer(ik), intent(in) :: k
      character(len=20)       :: nm
      if (k >= 1_ik .and. k <= N_COHORT_BOUNDED) then
         nm = COHORT_BOUNDED_NAME(k)
      else
         nm = '(none)'
      end if
   end function cohort_field_name

   !---------------------------------------------------------------------------------------!
   ! Initialize ONE freshly-created cohort slot at BIRTH from (pft, patch, dbh): set the      !
   ! per-slot identity + prognostic fields to their birth values, then derive                 !
   ! the on-allometry geometry via set_cohort_size. The SINGLE canonical birth-field policy    !
   ! shared by add_cohort (setup / census) and apply_recruitment -- add any new per-cohort     !
   ! field that needs a birth value HERE. overtopping_lai is zeroed (a reused, culled slot may  !
   ! hold a stale value; it is recomputed by update_overtopping_lai every slow step anyway).    !
   ! The caller ensures capacity, stamps the global id, and bumps cohort%n.                    !
   !---------------------------------------------------------------------------------------!
   subroutine init_cohort(cohort, m, pft, ipft, owner_patch, nplant, dbh, hydraulics, birth_temp, birth_psi)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: m, ipft, owner_patch
      type(pft_table_t),  intent(in)    :: pft
      real(wp),           intent(in)    :: nplant, dbh
      !----- The shared [hydraulics] traits, which with the PFT's own give its tissue curves. -------!
      type(hydraulics_config_t), intent(in) :: hydraulics
      !----- The state the cohort is born in, its patch's (patch_birth_state): its leaves and wood  !
      !      start at birth_temp, holding the water their tissues hold at birth_psi. ---------------!
      real(wp),           intent(in)    :: birth_temp   !< [K]
      real(wp),           intent(in)    :: birth_psi    !< [MPa]
      type(water_curve_t) :: leaf_curve, wood_curve
      cohort%pft(m)              = ipft
      cohort%owner_patch(m)      = owner_patch
      cohort%nplant(m)           = nplant
      cohort%dbh(m)              = dbh
      cohort%growth_avg(m)       = GROWTH_AVG_UNSET   ! set on its first growth step
      !----- UNSET: the fast loop seeds it from the surface-layer soil potential on the cohort's  !
      !      first step, so a recruit starts at its patch's actual water status, not at 0 (turgid). !
      cohort%dmax_psi_leaf(m)   = DMAX_PSI_LEAF_UNSET
      cohort%dmax_psi_leaf_accum(m) = DMAX_PSI_LEAF_ACCUM_RESET
      cohort%growth_accum(m)     = 0.0_wp
      cohort%growth_count(m)     = 0_ik
      cohort%growth_hist(:, m)   = 0.0_wp             ! a reused slot keeps a dead cohort's samples otherwise
      cohort%overtopping_lai(m)  = 0.0_wp             ! fresh competition context (recomputed each slow step)
      cohort%light_hours_accum(m)        = 0.0_wp             ! no light seen yet this step
      call cohort_diag_clear_slot(cohort%diag,  m)    ! fresh diagnostics (slot may be a reused, stale cull)
      call cohort_diag_clear_slot(cohort%sdiag, m)
      cohort%leaf_flush_tendency(m) = PHENO_FLUSH_INIT  ! born flushing
      cohort%leaf_shed_tendency(m)  = PHENO_SHED_INIT   ! not senescing at birth
      cohort%growing_degree_days(m) = 0.0_wp            ! fresh phenology memory
      cohort%cold_degree_days(m)    = 0.0_wp
      cohort%dry_psi_sum(m)         = 0.0_wp
      cohort%wet_psi_sum(m)         = 0.0_wp
      cohort%light_hours_mean(m) = -1.0_wp          ! no light memory: its first day sets it
      cohort%sla(m)                = pft%sla(ipft)
      cohort%vcmax25(m)              = pft%vcmax25(ipft)          ! dynamic leaf traits born at
      cohort%rd25(m)                 = pft%rd25(ipft)             ! top-of-canopy; plasticity relaxes
      cohort%leaf_turnover_rate(m)   = pft%leaf_turnover_rate_toc(ipft) ! them toward the shaded target
      cohort%leaf_temp(m)        = birth_temp         ! fresh fast state (slot may be a reused, stale cull)
      cohort%wood_temp(m)        = birth_temp
      cohort%leaf_surf_water(m)  = 0.0_wp             ! a new cohort's canopy starts dry
      cohort%wood_surf_water(m)  = 0.0_wp
      call set_cohort_size(cohort, m)                 ! carbon pools from dbh, then the geometry from them
      call pft_tissue_curves(hydraulics, pft, ipft, leaf_curve, wood_curve)
      call set_cohort_tissue_water(cohort, m, birth_psi, leaf_curve, wood_curve)
   end subroutine init_cohort

   !----- The carbon whose water a cohort's wood store holds [kgC/plant]: the sapwood ring and the   !
   !      fine roots, one store in the plant-water network (docs/science/soil_plant_water.md §1). ---!
   elemental pure real(wp) function cohort_wood_store_carbon(sapwood_carbon, fineroot_carbon) result(c)
      real(wp), intent(in) :: sapwood_carbon, fineroot_carbon   !< [kgC/plant]
      c = sapwood_carbon + fineroot_carbon
   end function cohort_wood_store_carbon

   !----- Fill a cohort's leaf and wood stores with the water they hold at potential `psi` on their  !
   !      own curves (docs/science/plant_hydraulics.md eq. (5)): a cohort's water at birth. ---------!
   pure subroutine set_cohort_tissue_water(cohort, i, psi, leaf_curve, wood_curve)
      type(cohort_block),  intent(inout) :: cohort
      integer(ik),         intent(in)    :: i
      real(wp),            intent(in)    :: psi                    !< [MPa]
      type(water_curve_t), intent(in)    :: leaf_curve, wood_curve
      cohort%leaf_water_mass(i) = water_content(psi, leaf_curve, cohort%leaf_carbon(i))
      cohort%wood_water_mass(i) = water_content(psi, wood_curve,                                     &
                                     cohort_wood_store_carbon(cohort%sapwood_carbon(i), cohort%fineroot_carbon(i)))
   end subroutine set_cohort_tissue_water

   !---------------------------------------------------------------------------------------!
   ! The ONE writer of a cohort's geometry: with wood_carbon the size anchor, derive dbh from   !
   ! it (wood_to_dbh) and the rest of the geometry, and take leaf_area straight from leaf_carbon. !
   ! Called after every operator that moves carbon or SLA (birth, the daily update, trait       !
   ! acclimation, fusion, fission, a restart read). The pools are the inputs and are not        !
   ! changed.                                                                                     !
   !---------------------------------------------------------------------------------------!
   subroutine set_cohort_size_from_carbon(cohort, i)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: i
      !----- The carbon->geometry flip is the single shared kernel carbon_to_structure             !
      !       (meds_allometry); this per-slot wrapper just packs/unpacks the SoA fields.           !
      integer(ik) :: pf
      pf = cohort%pft(i)
      call carbon_to_structure(cohort%allom, cohort%wood_carbon(i), cohort%leaf_carbon(i),        &
                               cohort%pft_table%wood_density(pf), cohort%pft_table%hgt_max(pf),  &
                               cohort%pft_table%aboveground_frac(pf), cohort%sla(i),             &
                               cohort%dbh(i), cohort%height(i), cohort%basal_area(i),             &
                               cohort%agb(i), cohort%leaf_area(i))
      call set_cohort_wood_geometry(cohort, i)
   end subroutine set_cohort_size_from_carbon

end module meds_cohort_state_vars
