! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_patch_state_vars -- what one patch holds, and the block of all the site's patches.  !
!                                                                                          !
! One patch's stores, each its own type: the canopy-air-space twins, the soil water and    !
! soil energy columns, the snow store, the slow soil-carbon pools, the daily fast->slow    !
! decomposition accumulator and the litter input -- plus the area-weighted blends that     !
! mix two of each when patches fuse or a disturbance gap is carved. The kernels take       !
! these store types from here, one patch at a time, and nothing else from the state        !
! layer. The parameters that describe the stores (soil_params_t) are derived from the      !
! config and live with it, in meds_soil_column_config -- a parameter is not state          !
! (src/README.md, rule 3).                                                                 !
!                                                                                          !
! `patch_block` holds every patch of the site as one Structure-of-Arrays, with a CSR map   !
! into the cohort block (meds_cohort_state_vars).                                          !
!==========================================================================================!
module meds_patch_state_vars
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : n_soil_layer_max, n_snow_layer_max, tiny_num, t_3ple
   use meds_therm_lib,        only : cas_molar_density
   use meds_diag_vars,        only : patch_diag_block, patch_diag_grow, patch_diag_reorder
   use meds_state_defaults,   only : CAN_CO2_UNSEEDED, CAN_DEPTH_UNSET
   use meds_config,           only : meds_config_t, DIST_PRIMARY
   use meds_soil_column_config, only : soil_params_t, soil_column_start_theta, top_soil_psi
   implicit none
   private

   public :: cas_state_t, cas_set_depth
   public :: snow_column_t, soil_water_column_t, soil_energy_column_t, soil_carbon_t
   public :: xi_accum_t   !< daily fast->slow accumulator for the soil-carbon matrix
   public :: litter_input_t  !< per-patch, per-day litter arriving at the CENTURY pools
   public :: blend_cas, blend_snow, blend_soil_water, blend_soil_energy, blend_soil_carbon, blend_xi_accum
   public :: blend_litter_input
                                                                    !< area-weighted mix (patch fusion / disturbance seed)
   public :: patch_block, alloc_patch_block, expand_patch_block, reorder_patch_block, free_patch_block, patch_physical_state

   !----- Prognostic per-patch canopy-air-space thermal state (three implicit twins). -------!
   type :: cas_state_t
      real(wp) :: can_enthalpy = 0.0_wp                     !< [J/kg] specific enthalpy (PROGNOSTIC)
      real(wp) :: can_shv      = 0.0_wp                     !< [kg/kg] specific humidity (PROGNOSTIC twin)
      real(wp) :: can_co2      = CAN_CO2_UNSEEDED           !< [umol/mol] dry-air CO2 mixing ratio (PROGNOSTIC third twin)
      real(wp) :: can_temp     = 0.0_wp                     !< [K]    diagnosed
      !----- CAS DEPTH is PROGNOSTIC-ish: a per-patch geometry state the SLOW loop owns, set from     !
      !      the tallest cohort plus a freeboard (aero_cfg_t%canopy_freeboard, floored at             !
      !      min_canopy_depth). It sizes cas_mass_capacity and cas_molar_capacity: one fixed          !
      !      depth would hold 4x too much air over a 1 m regenerating gap and ~1.8x too little        !
      !      over a 35 m tropical canopy. The 20 m here is only the value before the slow loop        !
      !      first sets it. --------------------------------------------------------------------------!
      real(wp) :: can_depth    = CAN_DEPTH_UNSET            !< [m]    CAS depth (slow loop owns it)
   end type cas_state_t

   !----- Temporary-surface-water / SNOW store: a stacked mass+energy reservoir between the       !
   !      CAS and soil_energy(1) (docs/science/snow_biophysics.md). PROGNOSTIC = water-           !
   !      equivalent mass (swe) + EXTENSIVE internal energy (J/m2, unlike soil_energy's           !
   !      J/m3); temperature + liquid fraction are read-offs of internal_energy_to_temp           !
   !      (dry_hcap=0). nlayer=0 is the "no snow" state (store present but empty). ---------------!
   type :: snow_column_t
      real(wp) :: swe(n_snow_layer_max)         = 0.0_wp   !< [kg/m2] water-equivalent mass  (PROGNOSTIC)
      real(wp) :: snow_energy(n_snow_layer_max) = 0.0_wp   !< [J/m2]  extensive internal energy (PROGNOSTIC)
      real(wp) :: snow_depth(n_snow_layer_max)  = 0.0_wp   !< [m]     geometric depth (= swe/rho_snow, the fixed density)
      real(wp) :: snow_temp(n_snow_layer_max)   = t_3ple    !< [K]    diagnosed each step (read-off)
      real(wp) :: snow_fliq(n_snow_layer_max)   = 0.0_wp   !< [-]     diagnosed liquid fraction (read-off)
      integer(ik) :: nlayer                     = 0_ik     !< active layer count (0 = no snow)
   end type snow_column_t

   !----- Prognostic per-patch soil WATER column (the value the hydrology kernel updates). --!
   type :: soil_water_column_t
      real(wp) :: theta(n_soil_layer_max) = 0.0_wp   !< [m3/m3] volumetric soil moisture (PROGNOSTIC)
      real(wp) :: w_surface = 0.0_wp                 !< [kg/m2] ponded surface water
      !----- ENTHALPY of the ponded water, PROGNOSTIC so that every pond seam is a PAIRED (mass,       !
      !      enthalpy) transfer, the same discipline the snow pack follows                             !
      !      (docs/science/numerical_scheme.md §5b, "The pond holds heat"). A mass-only pond lets      !
      !      water cross a seam and shed its enthalpy at the source layer's temperature while the      !
      !      water itself still holds that heat on the surface: both sides of the books agree and      !
      !      the column loses energy that has not gone anywhere. EXTENSIVE [J/m2], not volumetric,     !
      !      so the transfers are plain additions; temperature is a read-off of                        !
      !      internal_energy_to_temp with dry_hcap = 0, exactly as for snow -- which means             !
      !      freeze/thaw of ponded water comes for free. ----------------------------------------------!
      real(wp) :: w_surface_energy = 0.0_wp            !< [J/m2] ponded-water internal energy (PROGNOSTIC)
   end type soil_water_column_t

   !----- Prognostic per-patch soil THERMAL column (internal energy; temp/fliq diagnosed). --!
   type :: soil_energy_column_t
      real(wp) :: soil_energy(n_soil_layer_max) = 0.0_wp    !< [J/m3] volumetric internal energy (PROGNOSTIC)
      real(wp) :: soil_temp(n_soil_layer_max)   = 0.0_wp    !< [K]    diagnosed each step
      real(wp) :: soil_fliq(n_soil_layer_max)   = 1.0_wp    !< [-]    diagnosed liquid fraction
   end type soil_energy_column_t

   !==========================================================================================!
   !  Slow, stateful per-patch soil-carbon pools (docs/science/soil_carbon.md §2), written    !
   !  DAILY by meds_soil_biogeochem%soil_carbon_step and READ-ONLY, frozen across the day, by !
   !  the fast loop's heterotrophic respiration. Its conceptual home is                       !
   !  meds_biogeochem_types, but the pools are what a patch holds, so they live here with the !
   !  patch's other stores and the soil kernels take them from here. The matrix state vector  !
   !  X, ordered litter -> SOM -> passive (n_soil_pool=7 in meds_biogeochem_types);           !
   !  pack_pool_vector / unpack_pool_vector marshal the named fields to and from the array    !
   !  the kernels use. Lignin is a passive tracer of the structural pools (0 <= L <= C),      !
   !  outside the carbon-mass vector. There is no nitrogen (§10).                             !
   !==========================================================================================!
   type :: soil_carbon_t
      ! carbon pools [kgC/m2] -- the matrix state vector X, ordered litter -> SOM -> passive
      real(wp) :: fast_grnd_carbon    = 0.0_wp   !< X(1) metabolic litter, above-ground (flammable)
      real(wp) :: fast_soil_carbon    = 0.0_wp   !< X(2) metabolic litter, below-ground
      real(wp) :: struct_grnd_carbon  = 0.0_wp   !< X(3) structural litter + CWD, above
      real(wp) :: struct_soil_carbon  = 0.0_wp   !< X(4) structural litter + CWD, below
      real(wp) :: microbial_carbon    = 0.0_wp   !< X(5) microbial SOM        (scheme-5 only)
      real(wp) :: slow_carbon         = 0.0_wp   !< X(6) slow / humified SOM
      real(wp) :: passive_carbon      = 0.0_wp   !< X(7) passive SOM          (scheme-5 only)
      ! lignin sub-state of the structural pools [kgC/m2] (fraction f_lignin = L/C brakes decomposition)
      real(wp) :: struct_grnd_lignin  = 0.0_wp
      real(wp) :: struct_soil_lignin  = 0.0_wp
   end type soil_carbon_t

   !==========================================================================================!
   !  Daily fast->slow accumulator for the soil-carbon matrix (docs/science/soil_carbon.md    !
   !  §7, eq. (5)): the per-pool day-INTEGRAL of the environmental decomposition scalar,      !
   !  `xi_int_j = INT_day xi_j dt` [day], accumulated once per (patch, fast sub-step) from    !
   !  the column pre-pass (patch_heterotrophic_respiration) over the SAME frozen pool the     !
   !  fast loop's heterotrophic_respiration_matrix respires against. Reset to 0 at the start  !
   !  of each slow step (mirrors cohort%gpp_accum's reset in fast_dynamics), consumed by the  !
   !  daily soil_carbon_step, which decrements each donor pool by dvec_j = xi_int_j*K_j.      !
   !  `rh_fast_accum` is the day's ACCUMULATED fast-loop Rh (audit-only cross-check against   !
   !  soil_carbon_step's own daily Rh: the rh_seam_gap of §7). Named fields (not an           !
   !  n_soil_pool-sized array) for the SAME reason soil_carbon_t uses named fields: the pool  !
   !  count lives in meds_biogeochem_types, a soil-kernel module the state layer must not     !
   !  use; meds_soil_biogeochem's pack/unpack_pool_vector marshal to/from the array form      !
   !  kernels use.                                                                            !
   !==========================================================================================!
   type :: xi_accum_t
      real(wp) :: fast_grnd   = 0.0_wp   !< INT_day xi_1 dt [day]
      real(wp) :: fast_soil   = 0.0_wp   !< INT_day xi_2 dt [day]
      real(wp) :: struct_grnd = 0.0_wp   !< INT_day xi_3 dt [day]
      real(wp) :: struct_soil = 0.0_wp   !< INT_day xi_4 dt [day]
      real(wp) :: microbial   = 0.0_wp   !< INT_day xi_5 dt [day]
      real(wp) :: slow        = 0.0_wp   !< INT_day xi_6 dt [day]
      real(wp) :: passive     = 0.0_wp   !< INT_day xi_7 dt [day]
      real(wp) :: rh_fast_accum = 0.0_wp !< [kgC/m2] today's accumulated fast-loop Rh (audit-only)
   end type xi_accum_t

   !----- Per-patch, per-day litter input -- ALREADY PARTITIONED to pool destinations               !
   !      (docs/science/soil_carbon.md §5). The driver sums each cohort's                           !
   !      leaf/fineroot/wood/storage necromass into these bins USING that cohort's PFT              !
   !      f_labile_leaf/f_labile_stem and aboveground_frac (necromass_to_litter); the kernel        !
   !      maps them straight onto u. Units [kgC/m2/day]. -------------------------------------------!
   type :: litter_input_t
      real(wp) :: labile_grnd  = 0.0_wp      !< -> X(1) fast_grnd  (labile leaf, fruit; above-ground storage + wood)
      real(wp) :: labile_soil  = 0.0_wp      !< -> X(2) fast_soil  (labile fine root; below-ground storage + wood)
      real(wp) :: struct_grnd  = 0.0_wp      !< -> X(3) struct_grnd (structural leaf + CWD, above)
      real(wp) :: struct_soil  = 0.0_wp      !< -> X(4) struct_soil (structural fineroot + belowground CWD)
      real(wp) :: lignin_grnd  = 0.0_wp      !< lignin flux to struct_grnd (sets f_lignin of incoming litter)
      real(wp) :: lignin_soil  = 0.0_wp      !< lignin flux to struct_soil
      ! the same carbon by the tissue it came from (diagnostic; the soil reads the pool streams
      ! above). Leaf litter here is what a litter trap collects: after resorption. Summed, these
      ! equal labile_grnd + labile_soil + struct_grnd + struct_soil.
      real(wp) :: leaf_c     = 0.0_wp      !< leaves (turnover, drop, dead plants)
      real(wp) :: fineroot_c = 0.0_wp      !< fine roots
      real(wp) :: wood_c     = 0.0_wp      !< wood of dead plants (stem + coarse roots)
      real(wp) :: storage_c  = 0.0_wp      !< non-structural carbon of dead plants
      real(wp) :: repro_c    = 0.0_wp      !< reproduction carbon that made no recruit (fruit, seed)
      real(wp) :: exudate_c  = 0.0_wp      !< root exudate (carbon the growth sink could not use)
   end type litter_input_t
   !---------------------------------------------------------------------------------------!
   ! Patches of the site_t (SoA) with a CSR map into the cohort block.                        !
   !---------------------------------------------------------------------------------------!
   type :: patch_block
      !==========================================================================================!
      ! 0. BLOCK RECORD                                                                          !
      !==========================================================================================!
      integer(ik) :: n   = 0_ik       !< patches in use
      integer(ik) :: cap = 0_ik       !< allocated capacity
      !==========================================================================================!
      ! 1. IDENTITY AND AREA                                                                     !
      !==========================================================================================!
      !----- Persistent identity (as for cohorts): stamped at creation, carried through       !
      !      patch sort/fusion/compaction so a tracker can follow one patch across records.    !
      integer(ik), allocatable :: global_id(:)
      integer(ik), allocatable :: dist_type(:)      !< disturbance/land-use class
      real(wp),    allocatable :: area(:)           !< fraction, sum = 1 (kept real64)
      real(wp),    allocatable :: age(:)            !< [yr]
      !==========================================================================================!
      ! 2. COHORT MAP: a CSR map into the cohort block (rebuild_csr derives it)                  !
      !==========================================================================================!
      integer(ik), allocatable :: cohort_offset(:)  !< first cohort of the patch (1-based)
      integer(ik), allocatable :: cohort_count(:)   !< cohorts in the patch
      !==========================================================================================!
      ! 3. RECRUIT BANK                                                                          !
      !==========================================================================================!
      real(wp),    allocatable :: recruit_pool(:,:)  !< (pft, patch) carry-forward density
      !==========================================================================================!
      ! 4. CANOPY AIR                                                                            !
      !==========================================================================================!
      !----- PROGNOSTIC fast-biophysics reservoirs, here and in sections 5-6 (owned here so they ride
      !      the patch lockstep alongside area/age; mutated by the sub-daily fast loop, area-weighted
      !      on fusion).
      type(cas_state_t),          allocatable :: cas(:)         !< canopy-air-space twins per patch
      !==========================================================================================!
      ! 5. SNOW                                                                                  !
      !==========================================================================================!
      type(snow_column_t),        allocatable :: snow(:)        !< temporary-surface-water / snow store per patch
      !==========================================================================================!
      ! 6. SOIL WATER AND HEAT                                                                   !
      !==========================================================================================!
      type(soil_water_column_t),  allocatable :: soil_water(:)      !< soil water column per patch
      type(soil_energy_column_t), allocatable :: soil_energy(:)      !< soil thermal column per patch
      !==========================================================================================!
      ! 7. SOIL CARBON                                                                           !
      !==========================================================================================!
      !----- PROGNOSTIC slow biogeochemistry reservoir (owned here so it rides the patch             !
      !      lockstep alongside the fast reservoirs above; written DAILY by the slow soil-carbon     !
      !      step, area-weighted on fusion; docs/science/soil_carbon.md). ---------------------------!
      type(soil_carbon_t),        allocatable :: soil_carbon(:) !< soil-carbon pools per patch
      !==========================================================================================!
      ! 8. DAILY SLOW-TO-FAST RATES (set once per slow step, held fixed all day)                 !
      !==========================================================================================!
      !----- Daily slow->fast bridge for leaf/root-turnover shed water: SET (not accumulated)        !
      !      once per slow step by meds_vegetation_dynamics from this step's net leaf/fineroot       !
      !      carbon LOSS (proportional water shed, keeping the remaining tissue's rwc                !
      !      unchanged); consumed once per day by the fast loop's gather as a FROZEN additive        !
      !      input to the ground-water pathway (like rainfall/throughfall), held constant across     !
      !      every dt_fast sub-step of that day -- same "frozen daily quantity" convention as        !
      !      soil_carbon/xi_accum, but area-weighted like `age` on fusion (not blended from          !
      !      donors on a disturbance gap: a brand-new gap patch has no cohorts of its own that       !
      !      contributed to today's rate). ----------------------------------------------------------!
      real(wp),    allocatable :: shed_water_rate(:)   !< [kg/m2 ground/s] today's shed tissue water
      !----- The CO2 the SLOW tier owes the canopy air, as a frozen daily RATE the fast loop       !
      !      adds to nee_biotic (docs/science/plant_respiration.md eq. (7)) -- the carbon          !
      !      analogue of shed_water_rate, and set the same way. Two terms, both computed once      !
      !      per slow step by the allocator, that would otherwise go nowhere: GROWTH               !
      !      RESPIRATION (debited from the plant, never exhaled) and, with the opposite sign,      !
      !      the STARVATION DEFICIT (exhaled by the fast loop, debited from no pool because        !
      !      there is nothing left to debit). Signed, so it can be negative on a day the           !
      !      correction outweighs the growth. -----------------------------------------------------!
      real(wp),    allocatable :: slow_co2_rate(:)     !< [umol/m2 ground/s] growth resp - starvation over-report
      !==========================================================================================!
      ! 9. WITHIN-STEP ACCUMULATORS (reset each step, never restarted)                           !
      !==========================================================================================!
      !----- Daily fast->slow accumulator for the soil-carbon matrix: reset each slow step,           !
      !      accumulated once per (patch, fast sub-step) by column_prepass, consumed by the daily     !
      !      soil_carbon_step. Rides the SAME lockstep as soil_carbon (area-weighted on fusion). ------!
      type(xi_accum_t),           allocatable :: xi_accum(:)   !< [day]/[kgC/m2] per patch
      !----- Per-patch litter arriving at the CENTURY pools THIS slow step [kgC/m2/day]. Zeroed at !
      !      the top of the step by the vegetation driver, accumulated by turnover and mortality,   !
      !      consumed by advance_biogeochem_dynamics as the matrix source term u.                    !
      !                                                                                          !
      !      It lives HERE, on the patch block, because patch DISTURBANCE creates a gap partway      !
      !      through the step: a driver-local array sized to the patch count at the START of the     !
      !      step would be read past its end -- undefined memory into the litter input -- and on     !
      !      a FUSION would misattribute one patch's litter to another. As a patch field it          !
      !      rides the same lockstep as soil_carbon -- reordered on sort, packed on terminate,       !
      !      area-weighted on fusion and inherited by a gap. A slow ledger that declares the         !
      !      read as a boundary term would balance against the garbage, so the placement, not        !
      !      the ledger, is what keeps it right. ----------------------------------------------------!
      type(litter_input_t),       allocatable :: litter_in(:)  !< [kgC/m2/day] this step's litter
      !==========================================================================================!
      ! 10. SOLVER MEMORY                                                                        !
      !==========================================================================================!
      !----- Adaptive-controller warm start, PER PATCH. It is controller state belonging to a           !
      !      column, so it lives beside that column's other per-patch reservoirs and rides the          !
      !      patch lockstep (sort/fuse/compact/new-patch). Holding it on the shared per-patch scratch   !
      !      instead makes it loop-carried -- patch 1 cold-starts and the rest inherit whatever step    !
      !      their predecessor left -- which makes the answer depend on patch ORDER and, once the patch !
      !      axis is threaded, on THREAD SCHEDULING. --------------------------------------------------!
      real(wp),    allocatable :: adapt_dt_last(:)     !< [s] last accepted fast-loop controller step
      !==========================================================================================!
      ! 11. OUTPUT ACCUMULATORS (meds_diag_vars; never restarted)                                !
      !==========================================================================================!
      !----- TRANSIENT fast-loop DIAGNOSTIC accumulators (the patch twin of cohort%diag). -------!
      type(patch_diag_block) :: diag
   end type patch_block

contains

   !=======================================================================================!
   !  BLENDS: mixing two patches' stores (a fusion, a disturbance gap), top-down           !
   !=======================================================================================!


   !=======================================================================================!
   !  Area-weighted linear mixes: result = w1*a + w2*b. The caller passes NORMALIZED weights !
   !  (w1 = a1/(a1+a2), w2 = a2/(a1+a2)) so an intensive quantity (theta, enthalpy, [J/m3])   !
   !  is conserved on an area basis when two patches fuse or a disturbance gap is carved from  !
   !  its donors. Diagnosed fields (temp/fliq) mix too and are re-diagnosed next fast step.    !
   !=======================================================================================!
   !=========================================================================================!
   ! blend_cas -- mix two canopy-air control volumes on their AIR MASS, not their ground area.   !
   !                                                                                          !
   ! Enthalpy, specific humidity and CO2 mixing ratio are INTENSIVE -- per kg of air. What the   !
   ! merged volume must conserve is the EXTENSIVE content, area x rho x depth x value, and the   !
   ! merged depth is the area-weighted mean (the volumes add over the summed ground). Putting     !
   ! those together, the weight that conserves is area x depth:                                   !
   !                                                                                          !
   !     (a1+a2) . d_new . v_new = a1.d1.v1 + a2.d2.v2 ,   (a1+a2).d_new = a1.d1 + a2.d2          !
   !     =>  v_new = (a1.d1.v1 + a2.d2.v2) / (a1.d1 + a2.d2)                                      !
   !                                                                                          !
   ! Weighting on AREA alone drops the covariance term w1.w2.(d1-d2).(v1-v2), so fusing a    !
   ! tall patch with a gap would corrupt the canopy air's energy, humidity AND CO2 together, !
   ! all three by the same relative error -- of order 1e4 J per patch-fusion event. It is    !
   ! exact only when the two depths match, which is precisely the case that does not need a  !
   ! blend.                                                                                  !
   !                                                                                         !
   ! rho cancels: it is site-uniform, so it never has to be passed in. `can_depth` stays AREA-      !
   ! weighted -- it is the merged volume over the merged ground, which is what a depth is -- and    !
   ! that value equals the total mass weight, so the two are computed once.                         !
   !                                                                                          !
   ! `can_temp` is re-diagnosed from the blended enthalpy by the caller; it is mixed on the same    !
   ! weights only so an un-refreshed read is not wildly wrong.                                      !
   !=========================================================================================!
   pure function blend_cas(w1, a, w2, b) result(c)
      real(wp),          intent(in) :: w1, w2
      type(cas_state_t), intent(in) :: a, b
      type(cas_state_t)             :: c
      real(wp) :: m1, m2, mt
      m1 = w1 * a%can_depth
      m2 = w2 * b%can_depth
      mt = m1 + m2
      c%can_depth = mt                                   ! = w1*d1 + w2*d2, the area-weighted depth
      if (mt > tiny_num) then        ! below it, a blend has no air mass to conserve
         c%can_enthalpy = (m1 * a%can_enthalpy + m2 * b%can_enthalpy) / mt
         c%can_shv      = (m1 * a%can_shv      + m2 * b%can_shv     ) / mt
         c%can_co2      = (m1 * a%can_co2      + m2 * b%can_co2     ) / mt
         c%can_temp     = (m1 * a%can_temp     + m2 * b%can_temp    ) / mt
      else
         !----- Both volumes are empty: nothing to conserve, so fall back to the area weights     !
         !      rather than divide by zero. The result is carried only until the next refresh.    !
         c%can_enthalpy = w1 * a%can_enthalpy + w2 * b%can_enthalpy
         c%can_shv      = w1 * a%can_shv      + w2 * b%can_shv
         c%can_co2      = w1 * a%can_co2      + w2 * b%can_co2
         c%can_temp     = w1 * a%can_temp     + w2 * b%can_temp
      end if
   end function blend_cas

   !=========================================================================================!
   ! cas_set_depth -- resize the canopy-air control volume when the canopy height changes.       !
   !                                                                                          !
   ! The CAS is an OPEN control volume. Growing the canopy does not create air: it enlarges the  !
   ! well-mixed layer by entraining air from just above, and shrinking it detrains air back. The !
   ! entrained/detrained air is at essentially the canopy-air state (it is the air the CAS was    !
   ! already exchanging with, one aerodynamic timescale away), so the INTENSIVE state -- specific !
   ! enthalpy, specific humidity, CO2 mixing ratio -- is the invariant, and it is left untouched. !
   !                                                                                          !
   ! The EXTENSIVE content therefore changes: mass by rho*d(depth), energy by rho*d(depth)*enth. !
   ! That is a real exchange with the atmosphere, not a leak, and the three *_open outputs report !
   ! it -- one per twin, because all three intensive quantities are invariant across the resize    !
   ! and all three extensive contents therefore move together -- so a caller can book it rather    !
   ! than discover it as an unexplained jump. Two reasons it is done HERE and                      !
   ! on the slow step rather than inside the fast loop:                                          !
   !                                                                                          !
   !   * canopy height only changes on a slow step, so the fast ledger never has to carry a      !
   !     moving control volume -- its cas_mass_capacity is constant across every sub-step of a day;           !
   !   * the jumps that matter are not growth (a 0.003 m/day increment is ~7 J/m2, negligible)    !
   !     but DISTURBANCE and FUSION, where a 20 m canopy can become a 1 m gap in one step. Those  !
   !     happen in the slow loop, so this is where the term belongs.                             !
   !                                                                                          !
   ! Conserving TOTAL energy instead (rescaling can_enthalpy by the mass ratio) would be wrong:   !
   ! it would cool the canopy air simply because the trees grew.                                 !
   !=========================================================================================!
   pure subroutine cas_set_depth(cas, depth_new, rho_air, de_open, dw_open, dc_open)
      type(cas_state_t),  intent(inout) :: cas
      real(wp),           intent(in)    :: depth_new  !< [m] new CAS depth
      !----- All optional, and the *_open outputs need `rho_air` to mean anything. A caller that   !
      !      just wants the geometry updated omits them; one booking the exchange asks for the      !
      !      twins it tracks. ------------------------------------------------------------------!
      real(wp), optional, intent(in)    :: rho_air    !< [kg/m3]
      real(wp), optional, intent(out)   :: de_open    !< [J/m2]     energy    entrained (+) / detrained (-)
      real(wp), optional, intent(out)   :: dw_open    !< [kg/m2]    vapour    entrained (+) / detrained (-)
      real(wp), optional, intent(out)   :: dc_open    !< [umol/m2]  CO2       entrained (+) / detrained (-)
      real(wp) :: dmass, dmol
      dmass = 0.0_wp
      dmol = 0.0_wp
      if (present(rho_air)) then
         dmass = rho_air * (depth_new - cas%can_depth)
         dmol  = cas_molar_density(rho_air, cas%can_shv) * (depth_new - cas%can_depth)
      end if
      if (present(de_open)) de_open = dmass * cas%can_enthalpy
      if (present(dw_open)) dw_open = dmass * cas%can_shv
      if (present(dc_open)) dc_open = dmol  * cas%can_co2
      cas%can_depth = depth_new
   end subroutine cas_set_depth

   !----- Area-weighted mix of two snow stores (patch fusion / disturbance seed). The extensive       !
   !      mass + energy blend additively-per-area (like w_surface); depth blends by area weight;       !
   !      temp/fliq are RE-DIAGNOSED by the caller from the blended (swe, snow_energy), never blended.  !
   !      nlayer of the mix = 1 if any blended swe remains, else 0 (caller re-diagnoses / cleans up).   !
   pure function blend_snow(w1, a, w2, b) result(c)
      real(wp),            intent(in) :: w1, w2
      type(snow_column_t), intent(in) :: a, b
      type(snow_column_t)             :: c
      c%swe         = w1 * a%swe         + w2 * b%swe
      c%snow_energy = w1 * a%snow_energy + w2 * b%snow_energy
      c%snow_depth  = w1 * a%snow_depth  + w2 * b%snow_depth
      c%snow_temp   = w1 * a%snow_temp   + w2 * b%snow_temp   ! provisional; caller re-diagnoses from (swe, energy)
      c%snow_fliq   = w1 * a%snow_fliq   + w2 * b%snow_fliq   ! provisional; caller re-diagnoses
      c%nlayer      = max(a%nlayer, b%nlayer)                 ! caller collapses to 0 if blended swe < tiny
   end function blend_snow

   pure function blend_soil_water(w1, a, w2, b) result(c)
      real(wp),                  intent(in) :: w1, w2
      type(soil_water_column_t), intent(in) :: a, b
      type(soil_water_column_t)             :: c
      c%theta     = w1 * a%theta     + w2 * b%theta
      c%w_surface = w1 * a%w_surface + w2 * b%w_surface
      !----- EXTENSIVE, so it blends additively per area exactly like the mass it belongs to. ---------!
      c%w_surface_energy = w1 * a%w_surface_energy + w2 * b%w_surface_energy
   end function blend_soil_water

   pure function blend_soil_energy(w1, a, w2, b) result(c)
      real(wp),                  intent(in) :: w1, w2
      type(soil_energy_column_t), intent(in) :: a, b
      type(soil_energy_column_t)             :: c
      c%soil_energy = w1 * a%soil_energy + w2 * b%soil_energy
      c%soil_temp   = w1 * a%soil_temp   + w2 * b%soil_temp
      c%soil_fliq   = w1 * a%soil_fliq   + w2 * b%soil_fliq
   end function blend_soil_energy

   !----- Area-weighted mix of two soil-carbon columns (patch fusion / disturbance seed). Every    !
   !      field is a plain per-area density [kgC/m2] (no diagnosed/re-derived fields, unlike        !
   !      temp/fliq elsewhere), so a straight area-weighted average of all nine fields conserves    !
   !      total site-wide soil carbon exactly, mirroring blend_soil_water/blend_soil_energy.                 !
   pure function blend_soil_carbon(w1, a, w2, b) result(c)
      real(wp),             intent(in) :: w1, w2
      type(soil_carbon_t),  intent(in) :: a, b
      type(soil_carbon_t)              :: c
      c%fast_grnd_carbon   = w1 * a%fast_grnd_carbon   + w2 * b%fast_grnd_carbon
      c%fast_soil_carbon   = w1 * a%fast_soil_carbon   + w2 * b%fast_soil_carbon
      c%struct_grnd_carbon = w1 * a%struct_grnd_carbon + w2 * b%struct_grnd_carbon
      c%struct_soil_carbon = w1 * a%struct_soil_carbon + w2 * b%struct_soil_carbon
      c%microbial_carbon   = w1 * a%microbial_carbon   + w2 * b%microbial_carbon
      c%slow_carbon        = w1 * a%slow_carbon        + w2 * b%slow_carbon
      c%passive_carbon     = w1 * a%passive_carbon     + w2 * b%passive_carbon
      c%struct_grnd_lignin = w1 * a%struct_grnd_lignin + w2 * b%struct_grnd_lignin
      c%struct_soil_lignin = w1 * a%struct_soil_lignin + w2 * b%struct_soil_lignin
   end function blend_soil_carbon

   !----- Area-weighted mix of two daily xi accumulators (patch fusion / disturbance seed) --   !
   !      same rationale as blend_cas/blend_soil_water: an intra-day fusion should blend the         !
   !      partial-day accumulation exactly like the other per-patch fast reservoirs. -----------!
   pure function blend_xi_accum(w1, a, w2, b) result(c)
      real(wp),         intent(in) :: w1, w2
      type(xi_accum_t), intent(in) :: a, b
      type(xi_accum_t)             :: c
      c%fast_grnd     = w1 * a%fast_grnd     + w2 * b%fast_grnd
      c%fast_soil     = w1 * a%fast_soil     + w2 * b%fast_soil
      c%struct_grnd   = w1 * a%struct_grnd   + w2 * b%struct_grnd
      c%struct_soil   = w1 * a%struct_soil   + w2 * b%struct_soil
      c%microbial     = w1 * a%microbial     + w2 * b%microbial
      c%slow          = w1 * a%slow          + w2 * b%slow
      c%passive       = w1 * a%passive       + w2 * b%passive
      c%rh_fast_accum = w1 * a%rh_fast_accum + w2 * b%rh_fast_accum
   end function blend_xi_accum

   !----- Area-weighted mix of two per-patch litter accumulators (patch fusion / disturbance seed).!
   !                                                                                               !
   !      The litter accumulator is per-patch state the patch block owns, so like every other      !
   !      per-patch quantity it survives a patch operator through the reorder/pack/blend           !
   !      paths. A driver-local array sized to the patch count at the START of the slow step       !
   !      would not: patch disturbance CREATES a patch later in the same step, and the             !
   !      consumer would read past the array's end, feeding undefined memory into the CENTURY      !
   !      litter input.                                                                            !
   !                                                                                               !
   !      Every field is a per-area rate [kgC/m2/day], so an area-weighted mean is what conserves  !
   !      the site-wide litter flux -- the same rule as blend_soil_carbon, which is where this     !
   !      litter is headed.                                                                        !
   pure function blend_litter_input(w1, a, w2, b) result(c)
      real(wp),             intent(in) :: w1, w2
      type(litter_input_t), intent(in) :: a, b
      type(litter_input_t)             :: c
      c%labile_grnd   = w1 * a%labile_grnd   + w2 * b%labile_grnd
      c%labile_soil   = w1 * a%labile_soil   + w2 * b%labile_soil
      c%struct_grnd   = w1 * a%struct_grnd   + w2 * b%struct_grnd
      c%struct_soil   = w1 * a%struct_soil   + w2 * b%struct_soil
      c%lignin_grnd   = w1 * a%lignin_grnd   + w2 * b%lignin_grnd
      c%lignin_soil   = w1 * a%lignin_soil   + w2 * b%lignin_soil
      c%leaf_c        = w1 * a%leaf_c        + w2 * b%leaf_c
      c%fineroot_c    = w1 * a%fineroot_c    + w2 * b%fineroot_c
      c%wood_c        = w1 * a%wood_c        + w2 * b%wood_c
      c%storage_c     = w1 * a%storage_c     + w2 * b%storage_c
      c%repro_c       = w1 * a%repro_c       + w2 * b%repro_c
      c%exudate_c     = w1 * a%exudate_c     + w2 * b%exudate_c
   end function blend_litter_input

   !=======================================================================================!
   !  THE PATCH LOCKSTEP: the routines that list every patch field, in section order       !
   !=======================================================================================!

   subroutine alloc_patch_block(patch, cap, n_pft)
      type(patch_block), intent(inout) :: patch
      integer(ik),       intent(in)    :: cap, n_pft
      patch%cap = cap
      patch%n = 0_ik
      !----- 1. Identity and area ---------------------------------------------------------------!
      allocate(patch%global_id(cap), patch%dist_type(cap), patch%area(cap), patch%age(cap))
      patch%global_id = 0_ik
      patch%dist_type = DIST_PRIMARY
      patch%area = 0.0_wp
      patch%age = 0.0_wp
      !----- 2. Cohort map ----------------------------------------------------------------------!
      allocate(patch%cohort_offset(cap), patch%cohort_count(cap))
      patch%cohort_offset = 0_ik
      patch%cohort_count = 0_ik
      !----- 3. Recruit bank --------------------------------------------------------------------!
      allocate(patch%recruit_pool(n_pft, cap))
      patch%recruit_pool = 0.0_wp
      !----- 4. Canopy air ----------------------------------------------------------------------!
      allocate(patch%cas(cap))                              ! default-initialised by its type
      !----- 5. Snow ----------------------------------------------------------------------------!
      allocate(patch%snow(cap))                             ! default-initialised by its type
      !----- 6. Soil water and heat -------------------------------------------------------------!
      allocate(patch%soil_water(cap), patch%soil_energy(cap))        ! default-initialised by its type
      !----- 7. Soil carbon ---------------------------------------------------------------------!
      allocate(patch%soil_carbon(cap))                      ! default-initialised by its type
      !----- 8. Daily slow-to-fast rates --------------------------------------------------------!
      allocate(patch%shed_water_rate(cap), patch%slow_co2_rate(cap))
      patch%shed_water_rate = 0.0_wp
      patch%slow_co2_rate = 0.0_wp
      !----- 9. Within-step accumulators --------------------------------------------------------!
      allocate(patch%xi_accum(cap), patch%litter_in(cap))   ! default-initialised by its type
      !----- 10. Solver memory ------------------------------------------------------------------!
      allocate(patch%adapt_dt_last(cap))
      patch%adapt_dt_last = 0.0_wp
   end subroutine alloc_patch_block

   !=======================================================================================!
   !  expand_patch_block -- the patch-SoA twin of expand_cohort_block: ensure the          !
   !  patch SoA can hold at least `need` patches, growing 1.5x if it must. Same n/cap      !
   !  reserve model, no-op fast path, geometric growth, and move_alloc as in the           !
   !  cohort routine. LOCKSTEP SITE: a new per-patch field must be copied here (1:n)       !
   !  and move_alloc'd, and listed in alloc_patch_block, reorder_patch_block and free_patch_block.           !
   !=======================================================================================!
   subroutine expand_patch_block(patch, need, n_pft)
      type(patch_block), intent(inout) :: patch
      integer(ik),       intent(in)    :: need, n_pft
      type(patch_block)                :: tmp
      integer(ik)                       :: newcap, m
      if (need <= patch%cap) return
      newcap = max(need, (patch%cap * 3_ik) / 2_ik + 1_ik)
      !----- Grow the diagnostic block in place; the move_alloc list below leaves it alone.  -----!
      call patch_diag_grow(patch%diag, newcap)
      call alloc_patch_block(tmp, newcap, n_pft)
      m = patch%n
      tmp%n = m
      !----- 1. Identity and area ---------------------------------------------------------------!
      tmp%global_id(1:m) = patch%global_id(1:m)
      tmp%dist_type(1:m) = patch%dist_type(1:m)
      tmp%area(1:m)      = patch%area(1:m)
      tmp%age(1:m)       = patch%age(1:m)
      !----- 2. Cohort map ----------------------------------------------------------------------!
      tmp%cohort_offset(1:m) = patch%cohort_offset(1:m)
      tmp%cohort_count(1:m)  = patch%cohort_count(1:m)
      !----- 3. Recruit bank --------------------------------------------------------------------!
      tmp%recruit_pool(:,1:m) = patch%recruit_pool(:,1:m)
      !----- 4. Canopy air ----------------------------------------------------------------------!
      tmp%cas(1:m) = patch%cas(1:m)
      !----- 5. Snow ----------------------------------------------------------------------------!
      tmp%snow(1:m) = patch%snow(1:m)
      !----- 6. Soil water and heat -------------------------------------------------------------!
      tmp%soil_water(1:m) = patch%soil_water(1:m)
      tmp%soil_energy(1:m) = patch%soil_energy(1:m)
      !----- 7. Soil carbon ---------------------------------------------------------------------!
      tmp%soil_carbon(1:m) = patch%soil_carbon(1:m)
      !----- 8. Daily slow-to-fast rates --------------------------------------------------------!
      tmp%shed_water_rate(1:m) = patch%shed_water_rate(1:m)
      tmp%slow_co2_rate(1:m)   = patch%slow_co2_rate(1:m)
      !----- 9. Within-step accumulators --------------------------------------------------------!
      tmp%xi_accum(1:m)  = patch%xi_accum(1:m)
      tmp%litter_in(1:m) = patch%litter_in(1:m)
      !----- 10. Solver memory ------------------------------------------------------------------!
      tmp%adapt_dt_last(1:m) = patch%adapt_dt_last(1:m)
      patch%n = tmp%n
      patch%cap = tmp%cap
      !----- 1. Identity and area ---------------------------------------------------------------!
      call move_alloc(tmp%global_id, patch%global_id)
      call move_alloc(tmp%dist_type, patch%dist_type)
      call move_alloc(tmp%area, patch%area)
      call move_alloc(tmp%age, patch%age)
      !----- 2. Cohort map ----------------------------------------------------------------------!
      call move_alloc(tmp%cohort_offset, patch%cohort_offset)
      call move_alloc(tmp%cohort_count, patch%cohort_count)
      !----- 3. Recruit bank --------------------------------------------------------------------!
      call move_alloc(tmp%recruit_pool, patch%recruit_pool)
      !----- 4. Canopy air ----------------------------------------------------------------------!
      call move_alloc(tmp%cas, patch%cas)
      !----- 5. Snow ----------------------------------------------------------------------------!
      call move_alloc(tmp%snow, patch%snow)
      !----- 6. Soil water and heat -------------------------------------------------------------!
      call move_alloc(tmp%soil_water, patch%soil_water)
      call move_alloc(tmp%soil_energy, patch%soil_energy)
      !----- 7. Soil carbon ---------------------------------------------------------------------!
      call move_alloc(tmp%soil_carbon, patch%soil_carbon)
      !----- 8. Daily slow-to-fast rates --------------------------------------------------------!
      call move_alloc(tmp%shed_water_rate, patch%shed_water_rate)
      call move_alloc(tmp%slow_co2_rate, patch%slow_co2_rate)
      !----- 9. Within-step accumulators --------------------------------------------------------!
      call move_alloc(tmp%xi_accum, patch%xi_accum)
      call move_alloc(tmp%litter_in, patch%litter_in)
      !----- 10. Solver memory ------------------------------------------------------------------!
      call move_alloc(tmp%adapt_dt_last, patch%adapt_dt_last)
   end subroutine expand_patch_block

   !=======================================================================================!
   !  The patch twin of reorder_cohort_block: place old index perm(k) at new position k, for !
   !  k=1..m, in every per-patch array and the diagnostic block. Array assignment with a     !
   !  vector subscript on the RHS is evaluated to a temporary first, so self-permutation is  !
   !  safe. The CSR map is not permuted: rebuild_csr derives it from the cohorts afterwards. !
   !=======================================================================================!
   subroutine reorder_patch_block(patch, perm, m)
      type(patch_block), intent(inout) :: patch
      integer(ik),       intent(in)    :: perm(:)
      integer(ik),       intent(in)    :: m
      !----- 1. Identity and area ---------------------------------------------------------------!
      patch%global_id(1:m) = patch%global_id(perm(1:m))
      patch%dist_type(1:m) = patch%dist_type(perm(1:m))
      patch%area(1:m)      = patch%area(perm(1:m))
      patch%age(1:m)       = patch%age(perm(1:m))
      !----- 2. Cohort map ----------------------------------------------------------------------!
      !      (not permuted: rebuild_csr derives the map from the cohorts afterwards)
      !----- 3. Recruit bank --------------------------------------------------------------------!
      patch%recruit_pool(:,1:m) = patch%recruit_pool(:,perm(1:m))
      !----- 4. Canopy air ----------------------------------------------------------------------!
      patch%cas(1:m) = patch%cas(perm(1:m))
      !----- 5. Snow ----------------------------------------------------------------------------!
      patch%snow(1:m) = patch%snow(perm(1:m))
      !----- 6. Soil water and heat -------------------------------------------------------------!
      patch%soil_water(1:m) = patch%soil_water(perm(1:m))
      patch%soil_energy(1:m) = patch%soil_energy(perm(1:m))
      !----- 7. Soil carbon ---------------------------------------------------------------------!
      patch%soil_carbon(1:m) = patch%soil_carbon(perm(1:m))
      !----- 8. Daily slow-to-fast rates --------------------------------------------------------!
      patch%shed_water_rate(1:m) = patch%shed_water_rate(perm(1:m))
      patch%slow_co2_rate(1:m)   = patch%slow_co2_rate(perm(1:m))
      !----- 9. Within-step accumulators --------------------------------------------------------!
      patch%xi_accum(1:m)  = patch%xi_accum(perm(1:m))
      patch%litter_in(1:m) = patch%litter_in(perm(1:m))
      !----- 10. Solver memory ------------------------------------------------------------------!
      patch%adapt_dt_last(1:m) = patch%adapt_dt_last(perm(1:m))
      !----- 11. Output accumulators ------------------------------------------------------------!
      call patch_diag_reorder(patch%diag, perm, m)
      patch%n = m
   end subroutine reorder_patch_block

   !----- Free every per-patch array. -------------------------------------------------------!
   subroutine free_patch_block(patch)
      type(patch_block), intent(inout) :: patch
      patch%n = 0_ik
      patch%cap = 0_ik
      if (allocated(patch%area)) deallocate(patch%global_id, patch%dist_type, patch%area, patch%age,           &
                                            patch%cohort_offset, patch%cohort_count, patch%recruit_pool,       &
                                            patch%cas, patch%snow, patch%soil_water, patch%soil_energy,                 &
                                            patch%soil_carbon, patch%shed_water_rate, patch%slow_co2_rate,     &
                                            patch%xi_accum, patch%litter_in, patch%adapt_dt_last)
   end subroutine free_patch_block

   !=======================================================================================!
   !  BIRTH: the state a new cohort starts in                                              !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! The state a new cohort is born in, from its patch's canopy air `cas` and soil water      !
   ! `soil_water`: the canopy-air temperature [K] and the water potential of the top soil layer   !
   ! [MPa] (top_soil_psi, the value the stomatal predawn potential is seeded with), so the    !
   ! cohort starts in equilibrium with where it stands. A patch whose fast state has not been !
   ! seeded yet -- at setup, or in a run without the fast loop -- has neither, and gives the  !
   ! run's starting state from [init] instead: soil_temp, and the top layer's potential at    !
   ! the starting soil water. init_fast_reservoirs re-seeds the setup cohorts from their      !
   ! patch once it has seeded it. `soil` is the column's parameters (soil_column_params).     !
   !---------------------------------------------------------------------------------------!
   pure subroutine patch_physical_state(cas, soil_water, cfg, soil, temp, psi)
      type(cas_state_t),         intent(in)  :: cas      !< the patch's canopy air
      type(soil_water_column_t), intent(in)  :: soil_water   !< the patch's soil water
      type(meds_config_t),       intent(in)  :: cfg
      type(soil_params_t),       intent(in)  :: soil
      real(wp),                  intent(out) :: temp, psi
      real(wp) :: theta(n_soil_layer_max)
      !----- A seeded canopy air has a temperature: the type's default is 0 K. ----------------!
      if (cas%can_temp > 0.0_wp) then
         temp = cas%can_temp
         psi  = top_soil_psi(soil, soil_water%theta(1))
      else
         temp = cfg%init_soil_temp
         call soil_column_start_theta(cfg%soil_column, soil, cfg%init_soil_theta, theta)
         psi  = top_soil_psi(soil, theta(1))
      end if
   end subroutine patch_physical_state


end module meds_patch_state_vars
