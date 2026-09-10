!==========================================================================================!
! meds_slow_ledger -- the SITE conservation ledger across a slow step (plan §10.2).            !
!                                                                                          !
! WHY A SECOND LEDGER EXISTS AT ALL. `budget_t` already gives every fast store a closure test, !
! and fast_dynamics merges those into the run-level energy/water ledgers. But those accumulate  !
! per-fast-step FLUX RESIDUALS: each check asks whether the change over one dt_fast matches the !
! fluxes that crossed the store's boundary during it. Nothing anywhere compares a STORE across   !
! the SLOW step. A discontinuous jump between the end of one fast window and the start of the    !
! next -- a cull discarding tissue water, a fusion mixing canopy air on the wrong weights, a      !
! recruit endowed with carbon nobody paid for -- is therefore invisible BY CONSTRUCTION, not by   !
! oversight. That is why a conservation-focused review and two conservation PRs all passed over   !
! the items §10.2 lists.                                                                          !
!                                                                                          !
! So this ledger is the fast one turned inside out: a SNAPSHOT ledger. Store before, declared     !
! boundary terms, store after; whatever is left is the imbalance. It is SITE-LEVEL and AREA-       !
! WEIGHTED because patch identity does not survive the step -- apply_patch_disturbance creates a    !
! patch, fuse_2_patches destroys one, terminate_patches renormalizes every area -- so a per-patch   !
! ledger would have nothing to compare against.                                                     !
!                                                                                          !
! WHAT THIS MODULE DELIBERATELY DOES NOT DO: fix anything. Every term §10.2 identifies is still     !
! wrong after this lands. The point is to make each one MEASURABLE, attributed to the phase that    !
! produced it, before deciding what it deserves -- measurement inverted the expected ranking twice  !
! in #128/#129, and the prediction that growth respiration dominates the carbon bucket is recorded  !
! in §10.2 precisely so the numbers can overturn it. Terms that are already honest are DECLARED     !
! (see slow_ledger_declare); everything else lands in the phase's residual, which is the model's    !
! own admission of what it cannot account for.                                                      !
!                                                                                          !
! PHASE ATTRIBUTION is the whole value. A single number for the step would say "carbon does not     !
! close" and leave the search to the reader. Marking at each operator boundary says WHICH operator, !
! in WHICH currency, and by how much -- which is a bug report rather than a symptom.                 !
!==========================================================================================!
module meds_slow_ledger
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : rho_h2o
   use meds_config,             only : meds_config_t
   use meds_site_state_types,   only : site_t, cohort_tissue_heat_capacity,                    &
                                       TISSUE_C_LEAF, TISSUE_C_SAPW, TISSUE_HCAP_MIN
   use meds_column_params,      only : soil_params_t, build_soil_hydr_params
   use meds_therm_lib,          only : cas_molar_density
   use meds_allometry,          only : min_cohort_carbon
   use meds_budget_check,       only : budget_rtol_flux
   implicit none
   private

   public :: slow_store_t, slow_ledger_t, slow_site_store, slow_fast_carbon_handover
   public :: slow_tissue_heat
   public :: slow_ledger_open, slow_ledger_declare, slow_ledger_mark, slow_ledger_report
   public :: N_SLOW_PHASE, slow_phase_name, KGC_PER_UMOL_C
   public :: SLOW_PHASE_ALLOCATE, SLOW_PHASE_GROW, SLOW_PHASE_RECRUIT, SLOW_PHASE_COHORT,          &
             SLOW_PHASE_DISTURB, SLOW_PHASE_PATCH, SLOW_PHASE_CANOPY, SLOW_PHASE_SOILC

   !----- The phases a slow step is divided into. Each is marked at the boundary of the operator  !
   !      group that owns it, so a residual names the operator that produced it. -----------------!
   integer(ik), parameter :: SLOW_PHASE_ALLOCATE = 1_ik  !< phenology, traits, carbon allocation, turnover shed
   integer(ik), parameter :: SLOW_PHASE_GROW     = 2_ik  !< tendency build + state commit (growth AND mortality)
   integer(ik), parameter :: SLOW_PHASE_RECRUIT  = 3_ik  !< apply_recruitment
   integer(ik), parameter :: SLOW_PHASE_COHORT   = 4_ik  !< cohort fuse / terminate / split
   integer(ik), parameter :: SLOW_PHASE_DISTURB  = 5_ik  !< apply_patch_disturbance
   integer(ik), parameter :: SLOW_PHASE_PATCH    = 6_ik  !< patch sort / fuse / terminate (+ trailing cohort ops)
   integer(ik), parameter :: SLOW_PHASE_CANOPY   = 7_ik  !< refresh_canopy_depth (the open CAS control volume)
   integer(ik), parameter :: SLOW_PHASE_SOILC    = 8_ik  !< the CENTURY daily step
   integer(ik), parameter :: N_SLOW_PHASE        = 8_ik

   !----- Carbon of one mole of CO2 [kgC/umol]: 1e-6 mol * 0.012 kgC/mol. Turns the CAS mixing     !
   !      ratio into the same currency as the plant and soil pools, so all four carbon stores add. !
   real(wp), parameter :: KGC_PER_UMOL_C = 1.2e-8_wp
   !----- Absolute floors PER MARK, so the tolerance grows with the number of checks rather than  !
   !      demanding that a thousand round-off errors cancel. Judged alongside a relative test on  !
   !      the declared flux -- the same policy budget_check uses, and for the same reason: a      !
   !      residual means nothing measured against a store, only against what crossed the boundary.!
   real(wp), parameter :: SLOW_ATOL_CARBON = 1.0e-12_wp   !< [kgC/m2] per mark
   real(wp), parameter :: SLOW_ATOL_WATER  = 1.0e-9_wp    !< [kg/m2]  per mark
   real(wp), parameter :: SLOW_ATOL_ENERGY = 1.0e-3_wp    !< [J/m2]   per mark
   !----- ...plus a ROUND-OFF allowance proportional to the store. Those floors alone are wrong   !
   !      for a phase that DECLARES nothing but rearranges a large store: the structural          !
   !      operators permute, merge and renormalise ~25 kgC/m2 of live carbon, and the arithmetic  !
   !      cost of that is ~1e-11 kgC/m2 -- above a 3e-12 floor and utterly below anything the     !
   !      ledger exists to find (the smallest REAL term it caught was 2.96e-6). This is NOT the   !
   !      store-relative tolerance budget_check's header warns against: 1e-11 is a round-off      !
   !      bound, five orders tighter than the 1e-6 that let a sustained leak hide.  -------------!
   real(wp), parameter :: SLOW_ROUNDOFF_FRAC = 1.0e-11_wp !< [-] of the store, per mark

   !----- One site total per currency, per m2 of SITE (already area-weighted over patches). -------!
   type :: slow_store_t
      real(wp) :: carbon = 0.0_wp   !< [kgC/m2] live pools + necromass + canopy air + recruit transit
      real(wp) :: water  = 0.0_wp   !< [kg/m2]  soil + pond + snow + canopy vapour + tissue + film + shed transit
      real(wp) :: energy = 0.0_wp   !< [J/m2]   soil + pond + snow + canopy air + tissue heat
   end type slow_store_t

   !----- The accumulator. Lives for the run and is threaded from meds_main exactly like the fast   !
   !      loop's run_energy_budget / run_water_budget, so its lifetime is visible at the top level  !
   !      rather than hidden in site state.  ------------------------------------------------------!
   type :: slow_ledger_t
      logical             :: active = .false.       !< off => every entry point is a no-op
      !----- Soil layer geometry, built once from [soil_column]; only `dz` is used (store depth).  !
      type(soil_params_t) :: soil
      logical             :: soil_ready = .false.
      !----- The air density the canopy-air store is valued at. Supplied by the caller from the    !
      !      fast context (site-uniform, from the forcing). Held CONSTANT across the step, which   !
      !      is exact: no fast sub-step runs during a slow step, so rho cannot change inside the   !
      !      window, and the before/after difference is the true change in canopy-air content.     !
      real(wp)            :: rho_air = 0.0_wp
      !----- Rolling state within the current step. --------------------------------------------!
      type(slow_store_t)  :: at_mark                !< store when the current phase opened
      type(slow_store_t)  :: decl_in, decl_out      !< declared boundary terms for the CURRENT phase
      !----- Per-phase attribution over the run. A SIGNED sum as well as an absolute one, because  !
      !      a one-signed bias below any per-step threshold is invisible to a running maximum and  !
      !      only shows up as a drift over a season (the lesson budget_t's header records).        !
      type(slow_store_t)  :: resid_sum(N_SLOW_PHASE)
      type(slow_store_t)  :: resid_abs(N_SLOW_PHASE)
      type(slow_store_t)  :: resid_worst(N_SLOW_PHASE)
      type(slow_store_t)  :: decl_gross(N_SLOW_PHASE)
      integer(ik)         :: n_mark(N_SLOW_PHASE) = 0_ik
      integer(ik)         :: n_step = 0_ik
      !----- Whole-run endpoints, for the drift line of the report. -----------------------------!
      type(slow_store_t)  :: store_first, store_last
      logical             :: has_first = .false.
   end type slow_ledger_t

contains

   !---------------------------------------------------------------------------------------!
   ! Open a slow step: snapshot the store as the first phase's baseline and clear the declared  !
   ! terms. `rho_air` values the canopy-air store (see the type's comment on why holding it       !
   ! constant across the step is exact rather than an approximation).                             !
   !---------------------------------------------------------------------------------------!
   subroutine slow_ledger_open(ledger, site, cfg, rho_air)
      type(slow_ledger_t), intent(inout) :: ledger
      type(site_t),        intent(in)    :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp), optional,  intent(in)    :: rho_air
      if (.not. ledger%active) return
      if (present(rho_air)) ledger%rho_air = rho_air
      call ensure_soil_geometry(ledger, cfg)
      ledger%at_mark  = slow_site_store(site, cfg, ledger%soil, ledger%rho_air)
      ledger%decl_in  = slow_store_t()
      ledger%decl_out = slow_store_t()
      ledger%n_step   = ledger%n_step + 1_ik
      if (.not. ledger%has_first) then
         ledger%store_first = ledger%at_mark
         ledger%has_first   = .true.
      end if
   end subroutine slow_ledger_open

   !---------------------------------------------------------------------------------------!
   ! Declare a boundary term for the phase now open. A term belongs here ONLY if the quantity   !
   ! genuinely crosses the site boundary and the model means it to: external seed rain arriving,  !
   ! necromass leaving a run that models no soil carbon, canopy air entrained as the control      !
   ! volume grows. Anything declared here is subtracted from the phase's residual, so declaring    !
   ! a defect would hide it -- which is exactly what this module exists to prevent.                !
   !---------------------------------------------------------------------------------------!
   pure subroutine slow_ledger_declare(ledger, carbon_in, carbon_out, water_in, water_out,         &
                                       energy_in, energy_out)
      type(slow_ledger_t), intent(inout) :: ledger
      real(wp), optional,  intent(in)    :: carbon_in, carbon_out, water_in, water_out
      real(wp), optional,  intent(in)    :: energy_in, energy_out
      if (.not. ledger%active) return
      if (present(carbon_in))  ledger%decl_in%carbon  = ledger%decl_in%carbon  + carbon_in
      if (present(carbon_out)) ledger%decl_out%carbon = ledger%decl_out%carbon + carbon_out
      if (present(water_in))   ledger%decl_in%water   = ledger%decl_in%water   + water_in
      if (present(water_out))  ledger%decl_out%water  = ledger%decl_out%water  + water_out
      if (present(energy_in))  ledger%decl_in%energy  = ledger%decl_in%energy  + energy_in
      if (present(energy_out)) ledger%decl_out%energy = ledger%decl_out%energy + energy_out
   end subroutine slow_ledger_declare

   !---------------------------------------------------------------------------------------!
   ! Close the phase now open and attribute its imbalance:                                     !
   !                                                                                          !
   !     resid = (store_now - store_at_mark) - (declared_in - declared_out)                     !
   !                                                                                          !
   ! then reopen at the new store with the declared terms cleared, so the next phase starts from !
   ! what this one actually left behind. Marking with no intervening operator is harmless (a zero !
   ! residual), which is what lets a caller mark defensively around a conditional block.          !
   !---------------------------------------------------------------------------------------!
   subroutine slow_ledger_mark(ledger, site, cfg, phase)
      type(slow_ledger_t), intent(inout) :: ledger
      type(site_t),        intent(in)    :: site
      type(meds_config_t), intent(in)    :: cfg
      integer(ik),         intent(in)    :: phase
      type(slow_store_t) :: now, r
      if (.not. ledger%active) return
      if (phase < 1_ik .or. phase > N_SLOW_PHASE) return
      now = slow_site_store(site, cfg, ledger%soil, ledger%rho_air)
      r%carbon = (now%carbon - ledger%at_mark%carbon) - (ledger%decl_in%carbon - ledger%decl_out%carbon)
      r%water  = (now%water  - ledger%at_mark%water ) - (ledger%decl_in%water  - ledger%decl_out%water )
      r%energy = (now%energy - ledger%at_mark%energy) - (ledger%decl_in%energy - ledger%decl_out%energy)

      associate (s => ledger%resid_sum(phase), a => ledger%resid_abs(phase),                       &
                 w => ledger%resid_worst(phase), g => ledger%decl_gross(phase))
         s%carbon = s%carbon + r%carbon ; s%water = s%water + r%water ; s%energy = s%energy + r%energy
         a%carbon = a%carbon + abs(r%carbon)
         a%water  = a%water  + abs(r%water)
         a%energy = a%energy + abs(r%energy)
         w%carbon = max(w%carbon, abs(r%carbon))
         w%water  = max(w%water,  abs(r%water))
         w%energy = max(w%energy, abs(r%energy))
         g%carbon = g%carbon + abs(ledger%decl_in%carbon) + abs(ledger%decl_out%carbon)
         g%water  = g%water  + abs(ledger%decl_in%water)  + abs(ledger%decl_out%water)
         g%energy = g%energy + abs(ledger%decl_in%energy) + abs(ledger%decl_out%energy)
      end associate
      ledger%n_mark(phase) = ledger%n_mark(phase) + 1_ik

      ledger%at_mark    = now
      ledger%store_last = now
      ledger%decl_in    = slow_store_t()
      ledger%decl_out   = slow_store_t()
   end subroutine slow_ledger_mark

   !=======================================================================================!
   !  THE STORE. Every conserved quantity the site holds, area-weighted to one m2 of site.     !
   !                                                                                          !
   !  The water and energy store lists mirror meds_fast_ark's whole_water / whole_energy check  !
   !  term for term, so the two tiers cannot disagree about what exists; carbon is added here,   !
   !  since the fast whole-column ledger does not track it. One entry is a store only this   !
   !  tier can see, and omitting it would report a real carry-forward as a leak:                !
   !                                                                                          !
   !    * `recruit_pool * carbon_min` -- reproduction carbon converted to recruit density and     !
   !      carried forward until a pool crosses the spawn threshold. Valued at carbon_min because  !
   !      that is the conversion the model itself used; the gap between that and what init_cohort !
   !      actually endows is §10.2.2 item 3, and valuing it this way is what makes that gap land  !
   !      in the RECRUIT phase's residual instead of being absorbed silently.                     !
   !                                                                                          !
   !  `shed_water_rate` is deliberately NOT a store, though it looks like one. It is a HANDOFF   !
   !  to the fast tier, live across the FAST window rather than across the slow step: by the      !
   !  time the next slow step opens, that water is already in the soil and the rate variable      !
   !  still holds it. Carrying it as a store double-counts it at every open, which reads as a     !
   !  steady one-signed water leak in the allocate phase -- measured at ~8e-7 kg/step before this !
   !  was corrected. It is DECLARED instead, the exact mirror of the fast->slow carbon handover.  !
   !=======================================================================================!
   pure function slow_site_store(site, cfg, soil, rho_air) result(store)
      type(site_t),        intent(in) :: site
      type(meds_config_t), intent(in) :: cfg
      type(soil_params_t), intent(in) :: soil
      real(wp),            intent(in) :: rho_air
      type(slow_store_t)              :: store
      real(wp)    :: a, dmol, cap_mass, cap_mol, e_soil, w_soil
      real(wp)    :: cap_leaf, cap_wood, c_min(cfg%pft%n)
      integer(ik) :: ip, i, i0, i1, k, pf

      !----- Recruit-pool valuation, one per PFT (the same conversion compute_vital_rates used). -!
      do pf = 1_ik, cfg%pft%n
         c_min(pf) = min_cohort_carbon(cfg%pft%min_cohort_height, cfg%pft%wood_density(pf))
      end do

      do ip = 1_ik, site%patch%n
         a = site%patch%area(ip)
         if (a <= 0.0_wp) cycle

         !----- Soil column: theta is volumetric and soil_energy is per m3, so both integrate     !
         !      over the layer thicknesses (identical to soil_water_store / soil_energy_store).   !
         w_soil = 0.0_wp ; e_soil = 0.0_wp
         do k = 1_ik, soil%n_active
            w_soil = w_soil + site%patch%soil_w(ip)%theta(k)       * soil%dz(k) * rho_h2o
            e_soil = e_soil + site%patch%soil_e(ip)%soil_energy(k) * soil%dz(k)
         end do
         store%water  = store%water  + a * w_soil
         store%energy = store%energy + a * e_soil

         !----- Pond and snow: both already extensive per m2 of patch ground. --------------------!
         store%water  = store%water  + a * site%patch%soil_w(ip)%w_surface
         store%energy = store%energy + a * site%patch%soil_w(ip)%w_surface_enth
         store%water  = store%water  + a * sum(site%patch%snow(ip)%swe)
         store%energy = store%energy + a * sum(site%patch%snow(ip)%snow_energy)

         !----- Canopy air. The capacities are rho*depth and dmol*depth, exactly as the fast       !
         !      prepass forms them, so a depth change moves the store by the amount the open        !
         !      control volume entrained or detrained (plan §10.2.4).  ---------------------------!
         associate (cas => site%patch%cas(ip))
            dmol     = cas_molar_density(rho_air, cas%can_shv)
            cap_mass = rho_air * cas%can_depth
            cap_mol  = dmol    * cas%can_depth
            store%water  = store%water  + a * cap_mass * cas%can_shv
            store%energy = store%energy + a * cap_mass * cas%can_enthalpy
            store%carbon = store%carbon + a * cap_mol  * cas%can_co2 * KGC_PER_UMOL_C
         end associate

         !----- Necromass: the seven CENTURY pools (lignin is a sub-state of the structural pools, !
         !      not carbon of its own, so it is deliberately NOT summed here). -------------------!
         associate (sc => site%patch%soil_carbon(ip))
            store%carbon = store%carbon + a * (sc%fast_grnd_carbon + sc%fast_soil_carbon           &
                                             + sc%struct_grnd_carbon + sc%struct_soil_carbon       &
                                             + sc%microbial_carbon + sc%slow_carbon                &
                                             + sc%passive_carbon)
         end associate

         !----- Carbon in transit: the recruit carry-forward pool. -------------------------------!
         do pf = 1_ik, cfg%pft%n
            store%carbon = store%carbon + a * site%patch%recruit_pool(pf, ip) * c_min(pf)
         end do

         !----- Cohorts of this patch. nplant is per m2 of PATCH ground, so the patch area weight  !
         !      converts to site; the interception films are ALREADY ground-referenced and must    !
         !      not take the nplant factor (the distinction PR #119 had to fix by hand).           !
         i0 = site%patch%cohort_offset(ip)
         i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         associate (c => site%cohort)
            do i = i0, i1
               store%carbon = store%carbon + a * c%nplant(i) * (c%leaf_carbon(i)                   &
                            + c%fineroot_carbon(i) + c%wood_carbon(i) + c%nonstructural_carbon(i))
               store%water  = store%water  + a * c%nplant(i) * (c%leaf_water_mass(i) + c%wood_water_mass(i))
               store%water  = store%water  + a * (c%leaf_surf_water(i) + c%wood_surf_water(i))
               !----- Tissue heat, through the SHARED capacity (state/site), so the ledger, the    !
               !      demography operators and the slow driver cannot drift apart on what a        !
               !      cohort's thermal mass is. This module carried its own copy until the shared  !
               !      one existed; keeping two was always going to end one way.  -----------------!
               call cohort_tissue_heat_capacity(c, i, TISSUE_C_LEAF, TISSUE_C_SAPW,               &
                                                TISSUE_HCAP_MIN, cap_leaf, cap_wood)
               store%energy = store%energy + a * (cap_leaf * c%leaf_temp(i) + cap_wood * c%wood_temp(i))
            end do
         end associate
      end do

   end function slow_site_store

   !---------------------------------------------------------------------------------------!
   ! The site's TISSUE HEAT alone [J/m2 site] -- the same slice slow_site_store folds into the   !
   ! energy total, exposed so the driver can measure how much of it moved when biomass changed.   !
   !---------------------------------------------------------------------------------------!
   pure function slow_tissue_heat(site) result(e)
      type(site_t), intent(in) :: site
      real(wp)    :: e, cap_leaf, cap_wood
      integer(ik) :: ip, i, i0, i1
      e = 0.0_wp
      do ip = 1_ik, site%patch%n
         i0 = site%patch%cohort_offset(ip)
         i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         do i = i0, i1
            call cohort_tissue_heat_capacity(site%cohort, i, TISSUE_C_LEAF, TISSUE_C_SAPW,        &
                                             TISSUE_HCAP_MIN, cap_leaf, cap_wood)
            e = e + site%patch%area(ip) * (cap_leaf * site%cohort%leaf_temp(i)                    &
                                         + cap_wood * site%cohort%wood_temp(i))
         end do
      end do
   end function slow_tissue_heat

   !---------------------------------------------------------------------------------------!
   ! The carbon the FAST tier handed to the slow tier this step [kgC/m2 site]: gross GPP minus   !
   ! the maintenance respiration it has already paid back into the canopy air. That net sits in   !
   ! the per-cohort accumulators between the two tiers -- not in any store this ledger snapshots  !
   ! -- so without declaring it the day's photosynthate reads as carbon appearing from nowhere.   !
   !                                                                                          !
   ! It mirrors cohort_carbon_demand's own branch exactly, including the fast_biophysics_on=off   !
   ! stub, because the ledger has to declare what the allocator actually consumed rather than     !
   ! what the model would consume if the fast loop were running.                                  !
   !---------------------------------------------------------------------------------------!
   pure function slow_fast_carbon_handover(site, cfg) result(net_c)
      type(site_t),        intent(in) :: site
      type(meds_config_t), intent(in) :: cfg
      real(wp)    :: net_c, per_plant, a
      integer(ik) :: ip, i, i0, i1
      net_c = 0.0_wp
      do ip = 1_ik, site%patch%n
         a  = site%patch%area(ip)
         i0 = site%patch%cohort_offset(ip)
         i1 = i0 + site%patch%cohort_count(ip) - 1_ik
         associate (ch => site%cohort)
            do i = i0, i1
               if (cfg%fast_biophysics_on) then
                  per_plant = ch%gpp_accum(i) - (ch%leaf_resp_accum(i) + ch%stem_resp_accum(i)     &
                                                 + ch%root_resp_accum(i))
               else
                  per_plant = cfg%gpp_ref * ch%leaf_area(i) * cfg%dt_years
               end if
               net_c = net_c + a * ch%nplant(i) * per_plant
            end do
         end associate
      end do
   end function slow_fast_carbon_handover

   !----- Build the soil layer geometry once, from the same [soil_column] block the fast context  !
   !      uses, so the two tiers integrate over identical layer thicknesses. --------------------!
   subroutine ensure_soil_geometry(ledger, cfg)
      type(slow_ledger_t), intent(inout) :: ledger
      type(meds_config_t), intent(in)    :: cfg
      if (ledger%soil_ready) return
      associate (sc => cfg%soil_column)
         call build_soil_hydr_params(sc%n_layer, sc%retention, sc%depth, sc%grid_growth,           &
                                     sc%theta_sat, sc%theta_res, sc%ksat, sc%curve_par_a,          &
                                     sc%curve_par_n, sc%root_beta, sc%psi_fc, ledger%soil)
      end associate
      ledger%soil_ready = .true.
   end subroutine ensure_soil_geometry

   !----- Phase labels for the report. ---------------------------------------------------------!
   pure function slow_phase_name(phase) result(nm)
      integer(ik), intent(in) :: phase
      character(len=16)       :: nm
      select case (phase)
      case (SLOW_PHASE_ALLOCATE) ; nm = 'allocate'
      case (SLOW_PHASE_GROW)     ; nm = 'grow+mortality'
      case (SLOW_PHASE_RECRUIT)  ; nm = 'recruit'
      case (SLOW_PHASE_COHORT)   ; nm = 'cohort fuse/fiss'
      case (SLOW_PHASE_DISTURB)  ; nm = 'disturbance'
      case (SLOW_PHASE_PATCH)    ; nm = 'patch fuse/term'
      case (SLOW_PHASE_CANOPY)   ; nm = 'canopy depth'
      case (SLOW_PHASE_SOILC)    ; nm = 'soil carbon'
      case default               ; nm = '?'
      end select
   end function slow_phase_name

   !---------------------------------------------------------------------------------------!
   ! End-of-run report: one row per phase per currency that actually moved. Reports the SIGNED  !
   ! cumulative residual first, because that is the quantity a leak shows up in -- an unsigned    !
   ! maximum cannot distinguish round-off that cancels from a bias that accumulates.              !
   !---------------------------------------------------------------------------------------!
   subroutine slow_ledger_report(ledger)
      type(slow_ledger_t), intent(in) :: ledger
      integer(ik) :: p, nbad
      if (.not. ledger%active .or. ledger%n_step == 0_ik) return
      print '(a)', ''
      print '(a)', '=== slow-loop conservation ledger (site totals, per m2) ==================='
      print '(a,i0,a)', 'slow steps: ', ledger%n_step,                                             &
                        '   [residual = store change - declared boundary terms]'
      !----- `declared` is the gross boundary flux the phase DID account for. It is the scale the   !
      !      residual has to be read against: -3e-3 beside a declared 0.5 is a 0.6% loss, and the   !
      !      same -3e-3 beside a declared 0 is carbon with no provenance at all. -------------------!
      print '(a)', 'phase              currency    signed sum      sum |r|       worst     declared   marks'
      do p = 1_ik, N_SLOW_PHASE
         if (ledger%n_mark(p) == 0_ik) cycle
         call row(p, 'carbon [kgC]', ledger%resid_sum(p)%carbon, ledger%resid_abs(p)%carbon,       &
                  ledger%resid_worst(p)%carbon, ledger%decl_gross(p)%carbon)
         call row(p, 'water   [kg]', ledger%resid_sum(p)%water,  ledger%resid_abs(p)%water,        &
                  ledger%resid_worst(p)%water,  ledger%decl_gross(p)%water)
         call row(p, 'energy   [J]', ledger%resid_sum(p)%energy, ledger%resid_abs(p)%energy,       &
                  ledger%resid_worst(p)%energy, ledger%decl_gross(p)%energy)
      end do
      print '(a)', '--------------------------------------------------------------------------'
      print '(a,es12.4,a,es12.4,a,es12.4)', 'store drift   carbon ',                               &
            ledger%store_last%carbon - ledger%store_first%carbon, '   water ',                     &
            ledger%store_last%water  - ledger%store_first%water,  '   energy ',                    &
            ledger%store_last%energy - ledger%store_first%energy
      !----- The verdict. Every phase/currency is judged against the flux it declared, plus a     !
      !      per-mark absolute floor -- a residual measured against a STORE means nothing, which  !
      !      is the policy budget_check's header argues for at length. Reported rather than       !
      !      fatal, matching the fast loop's whole-column ledgers: a run that breaches this is    !
      !      telling you something, and stopping it mid-way tells you less than finishing it. ----!
      nbad = 0_ik
      do p = 1_ik, N_SLOW_PHASE
         if (ledger%n_mark(p) == 0_ik) cycle
         nbad = nbad + fails(ledger%resid_sum(p)%carbon, ledger%decl_gross(p)%carbon,             &
                             SLOW_ATOL_CARBON, ledger%store_last%carbon, ledger%n_mark(p))
         nbad = nbad + fails(ledger%resid_sum(p)%water,  ledger%decl_gross(p)%water,              &
                             SLOW_ATOL_WATER,  ledger%store_last%water,  ledger%n_mark(p))
         nbad = nbad + fails(ledger%resid_sum(p)%energy, ledger%decl_gross(p)%energy,             &
                             SLOW_ATOL_ENERGY, ledger%store_last%energy, ledger%n_mark(p))
      end do
      if (nbad == 0_ik) then
         print '(a)', 'VERDICT: every phase closes within tolerance on all three currencies.'
      else
         print '(a,i0,a)', 'WARNING: ', nbad, ' phase/currency residual(s) exceed tolerance --'
         print '(a)', '         a slow-step operator is moving something it does not declare.'
      end if
      print '(a)', '=========================================================================='
   contains
      !----- 1 if this residual breaches the mixed relative/absolute tolerance, else 0. ---------!
      pure integer(ik) function fails(resid, gross, atol_per_mark, store, nmark)
         real(wp),    intent(in) :: resid, gross, atol_per_mark, store
         integer(ik), intent(in) :: nmark
         real(wp) :: tol
         tol = budget_rtol_flux * abs(gross)                                                       &
             + real(nmark, wp) * (atol_per_mark + SLOW_ROUNDOFF_FRAC * abs(store))
         fails = 0_ik
         if (abs(resid) > tol) fails = 1_ik
      end function fails

      subroutine row(pp, label, s, aa, w, g)
         integer(ik),      intent(in) :: pp
         character(len=*), intent(in) :: label
         real(wp),         intent(in) :: s, aa, w, g
         if (aa <= 0.0_wp .and. g <= 0.0_wp) return
         print '(a18,1x,a12,4(1x,es12.4),1x,i7)', slow_phase_name(pp), label, s, aa, w, g,         &
               ledger%n_mark(pp)
      end subroutine row
   end subroutine slow_ledger_report

end module meds_slow_ledger
