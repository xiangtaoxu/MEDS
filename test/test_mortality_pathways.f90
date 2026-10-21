! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_mortality_pathways -- mortality carbon SPLIT BY PATHWAY.                                 !
!                                                                                          !
! A MEDS plant can die three ways, and without the split a stand thinning continuously and a      !
! stand being knocked over look alike in the output.                                               !
! The three pathways are mutually exclusive and jointly exhaustive:                                !
!                                                                                          !
!   BACKGROUND  the continuous hazard, applied every slow step to every cohort's density           !
!   CULL        a cohort dropping below the tracked size floor, on the monthly cohort restructure  !
!   DISTURB     the canopy killed when treefall opens a gap, on the annual patch restructure       !
!                                                                                          !
! Each is asserted against an oracle computed directly from the state it consumed -- the density   !
! the applier actually removed, the culled cohort's own pools, the disturbed fraction of the       !
! canopy's own pools -- so no physics number is baked in and the checks cannot go stale with the   !
! parameters.                                                                                     !
!                                                                                          !
! Each part also asserts that the OTHER TWO slots stay untouched. That is the real failure mode of !
! a pathway split: not a wrong total, but the same death counted twice, which no sum-to-the-total  !
! check can see and which is exactly what a user would be reading the split to rule out.           !
!                                                                                          !
! Units follow the block's contract: the slots hold (rate x weight in seconds) and are read        !
! through patch_diag_value, the same reader the output layer uses, never off the raw accumulator.  !
!==========================================================================================!
program test_mortality_pathways
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : yr_sec
   use meds_config,              only : meds_config_t
   use meds_cohort_state_vars,   only : set_cohort_size_from_carbon
   use meds_polygon_state_vars,  only : polygon_t
   use meds_diag_vars,           only : alloc_patch_diag, patch_diag_value, PA_MORT_C_BACKGROUND, PA_MORT_C_CULL, &
                                        PA_MORT_C_DISTURB
   use meds_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use meds_demography_cohort_fusefiss, only : terminate_cohorts
   use meds_demography_patch_fusefiss,  only : apply_patch_disturbance
   use meds_slow_dynamics,       only : advance_slow_dynamics, advance_boundary_dynamics
   use meds_test_support,        only : banner, build_test_config, check, check_close, test_report
   implicit none

   call banner('mortality carbon by pathway')
   call check_cull()
   call check_disturbance()
   call check_background()
   call check_exhausted()
   call test_report('test_mortality_pathways')

contains

   !----- Arm the per-patch diagnostic block the way a real run does, minus the registry: active, !
   !      sized, `n` set to the patch count, and carrying the weight ONE SLOW STEP of fast         !
   !      sub-steps accumulates. The weight is what patch_diag_value divides out, so a fixture     !
   !      that leaves it at zero reads every slot as 0 and passes nothing.  ---------------------!
   subroutine arm_patch_diag(poly, cfg)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      call alloc_patch_diag(poly%patch%diag, poly%patch%n + 2_ik, .true.)
      poly%patch%diag%n = poly%patch%n
      poly%patch%diag%w(1:poly%patch%diag%cap) = cfg%dt_slow
   end subroutine arm_patch_diag

   !----- Total carbon one individual of cohort `i` carries. A whole-individual death takes every !
   !      pool, which is what separates it from a turnover shed.  -------------------------------!
   pure real(wp) function plant_carbon(poly, i) result(c)
      type(polygon_t), intent(in) :: poly
      integer(ik),  intent(in) :: i
      c = poly%cohort%leaf_carbon(i) + poly%cohort%fineroot_carbon(i)                              &
        + poly%cohort%wood_carbon(i) + poly%cohort%nonstructural_carbon(i)
   end function plant_carbon

   !----- Site total of a patch diagnostic: Sum(area * value), the aggregation the output layer   !
   !      performs, times the elapsed years to turn the declared [kgC/m2/yr] rate back into the    !
   !      amount that flowed.  ------------------------------------------------------------------!
   real(wp) function site_amount(poly, cfg, slot) result(a)
      type(polygon_t),     intent(in) :: poly
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: slot
      real(wp)    :: x(16)
      integer(ik) :: n, p
      call patch_diag_value(poly%patch%diag, slot, x, n)
      a = 0.0_wp
      do p = 1_ik, n
         a = a + poly%patch%area(p) * x(p)
      end do
      a = a * (cfg%dt_slow / yr_sec)
   end function site_amount

   !=== CULL: a cohort below the tracked size floor. ==========================================!
   subroutine check_cull()
      type(meds_config_t) :: cfg
      type(polygon_t)     :: poly
      real(wp)            :: expect

      cfg = build_test_config()
      !----- soil_carbon_on OFF deliberately. Both this pathway's litter and the disturbance        !
      !      pathway's sit behind that switch, and the diagnostic must not: how much biomass died   !
      !      is demography. With the default (.true.) this part would pass even with the write on   !
      !      the wrong side of the guard, which is the one placement error it exists to catch.      !
      cfg%soil_carbon%soil_carbon_on    = .false.
      !----- Raise the density floor so a NORMAL cohort is culled: the alternative is a cohort so  !
      !      small its carbon is at the rounding floor, which would make the assertion vacuous.    !
      cfg%negligible_nplant = 0.30_wp

      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.20_wp, 20.0_wp)   ! below the floor -> culled
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.50_wp, 20.0_wp)   ! above it        -> kept
      call finalize_init(poly)
      call arm_patch_diag(poly, cfg)

      expect = 0.20_wp * plant_carbon(poly, 1_ik)
      call check(expect > 0.0_wp, 'cull fixture: the doomed cohort carries carbon')

      call terminate_cohorts(poly, cfg)
      call check(poly%cohort%n == 1_ik, 'cull fixture: exactly one cohort was culled')
      call check_close(site_amount(poly, cfg, PA_MORT_C_CULL), expect, 1.0e-12_wp,                 &
                       'PA_MORT_C_CULL == the culled cohort''s every pool x its density')
      call check_close(site_amount(poly, cfg, PA_MORT_C_BACKGROUND), 0.0_wp, 1.0e-12_wp,           &
                       'a cull is not also counted as background mortality')
      call check_close(site_amount(poly, cfg, PA_MORT_C_DISTURB), 0.0_wp, 1.0e-12_wp,              &
                       'a cull is not also counted as disturbance mortality')
   end subroutine check_cull

   !=== DISTURBANCE: the canopy killed when a gap opens. ======================================!
   subroutine check_disturbance()
      type(meds_config_t) :: cfg
      type(polygon_t)     :: poly
      real(wp)            :: frac, expect

      cfg = build_test_config()
      cfg%soil_carbon%soil_carbon_on = .false.          ! the split is demography -- see check_cull
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.20_wp, 40.0_wp)   ! tall  -> dies in the gap
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.50_wp,  3.0_wp)   ! short -> survives into it
      call finalize_init(poly)
      call check(poly%cohort%height(1) >= cfg%disturbance_survive_height,                          &
                 'disturbance fixture: cohort 1 is canopy')
      call check(poly%cohort%height(2) <  cfg%disturbance_survive_height,                          &
                 'disturbance fixture: cohort 2 is understorey')
      call arm_patch_diag(poly, cfg)

      !----- The oracle, on the state BEFORE the call, since the call removes the canopy cohort.   !
      !      One patch of unit area, so the site amount is just frac x its density x its pools.    !
      frac   = 1.0_wp - exp(-cfg%patch_disturbance_rate * 1.0_wp)
      expect = poly%patch%area(1) * frac * poly%cohort%nplant(1) * plant_carbon(poly, 1_ik)
      call check(expect > 0.0_wp, 'disturbance fixture: the canopy carries carbon')

      call apply_patch_disturbance(poly, cfg, 1.0_wp)

      !----- Asserted on the SITE AGGREGATE, not on the donor slot, because that is where the      !
      !      donor's area shrink lands: the carbon died on the frac of the donor that became the   !
      !      gap, and only Sum(area * value) after the shrink can show whether it was valued on    !
      !      the right ground.  --------------------------------------------------------------!
      call check_close(site_amount(poly, cfg, PA_MORT_C_DISTURB), expect, 1.0e-12_wp,              &
                       'PA_MORT_C_DISTURB site total == the canopy carbon the gap consumed')
      call check_close(site_amount(poly, cfg, PA_MORT_C_BACKGROUND), 0.0_wp, 1.0e-12_wp,           &
                       'disturbance is not also counted as background mortality')
      call check_close(site_amount(poly, cfg, PA_MORT_C_CULL), 0.0_wp, 1.0e-12_wp,                 &
                       'disturbance is not also counted as a cull')
   end subroutine check_disturbance

   !=== BACKGROUND: the continuous hazard, through the real slow driver. =======================!
   subroutine check_background()
      type(meds_config_t) :: cfg
      type(polygon_t)     :: poly
      real(wp)            :: n0, expect

      cfg = build_test_config()
      cfg%fast_biophysics_on = .false.       ! the hazard is a state question, not a GPP one
      cfg%soil_carbon%soil_carbon_on     = .false.       ! and the split is emitted whatever the soil pools do

      !----- ONE cohort, so the per-step sort_cohorts cannot permute the axis between the density  !
      !      captured before the step and the pools read after it.  -------------------------!
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.50_wp, 20.0_wp)
      call finalize_init(poly)
      call arm_patch_diag(poly, cfg)
      !----- With the fast loop off, the slow step weights the block itself, so the fixture's         !
      !      stand-in for the fast sub-steps' weight would count the step twice. --------------------!
      poly%patch%diag%w = 0.0_wp

      n0 = poly%cohort%nplant(1)
      call advance_slow_dynamics(poly, cfg)                      ! no month/year boundary
      call check(poly%cohort%n == 1_ik, 'background fixture: the cohort survived the step')
      call check(n0 - poly%cohort%nplant(1) > 0.0_wp, 'the step actually killed some density')

      !----- Valued on the density the applier REMOVED and the pools it LEFT: growth happens at the !
      !      old density and mortality at the new pools, which is the decomposition                 !
      !      accumulate_mortality_litter's header sets out. Reading the pools after the step is     !
      !      therefore not a convenience -- it is the definition.  ------------------------------!
      expect = (n0 - poly%cohort%nplant(1)) * plant_carbon(poly, 1_ik)
      call check_close(site_amount(poly, cfg, PA_MORT_C_BACKGROUND), expect, 1.0e-12_wp,           &
                       'PA_MORT_C_BACKGROUND == the removed density x the pools it left')
      call check_close(site_amount(poly, cfg, PA_MORT_C_CULL), 0.0_wp, 1.0e-12_wp,                 &
                       'background mortality is not also counted as a cull')
      call check_close(site_amount(poly, cfg, PA_MORT_C_DISTURB), 0.0_wp, 1.0e-12_wp,              &
                       'background mortality is not also counted as disturbance')
   end subroutine check_background

   !=== A cohort with no leaves and no storage dies; one living on its storage does not. =========!
   !    Both are leafless and dormant (no flush, so nothing grows) with a day's stem respiration    !
   !    and no photosynthesis, as in a dry season; only their storage differs. The first cannot pay  !
   !    its upkeep, so its density goes to the floor on its first step, through the background     !
   !    pathway, and the next month boundary culls it. The second pays its upkeep from storage      !
   !    through a 150-day dry season and loses density only at the growth-driven hazard.           !
   subroutine check_exhausted()
      type(meds_config_t) :: cfg
      type(polygon_t)     :: poly
      real(wp)            :: n0, expect
      integer(ik)         :: day

      cfg = build_test_config()
      cfg%fast_biophysics_on = .true.        ! the allocator reads the accumulators set below
      cfg%soil_carbon%soil_carbon_on     = .false.

      !----- No storage: dead on the first step, at the floor, its carbon down the background path. -!
      call leafless_stand(poly, cfg, 0.0_wp)
      call arm_patch_diag(poly, cfg)        ! fast loop on: the fixture's weight stands in for it
      n0     = poly%cohort%nplant(1)
      expect = (n0 - cfg%negligible_nplant) * plant_carbon(poly, 1_ik)
      call advance_slow_dynamics(poly, cfg)
      call check(poly%cohort%nplant(1) == cfg%negligible_nplant,                                   &
                 'no leaves, no storage: the density reaches the floor on the first step')
      call check_close(site_amount(poly, cfg, PA_MORT_C_BACKGROUND), expect, 1.0e-9_wp,            &
                       'no leaves, no storage: its carbon dies down the background pathway')
      do day = 2_ik, 30_ik
         call advance_slow_dynamics(poly, cfg)
      end do
      call advance_boundary_dynamics(poly, cfg, .true., .false.)
      call check(poly%cohort%n == 0_ik, 'no leaves, no storage: culled at the month boundary')

      !----- Storage for the season: the cohort lives on it, leafless, for 150 days. --------------!
      call leafless_stand(poly, cfg, 0.5_wp)
      n0 = poly%cohort%nplant(1)
      do day = 1_ik, 150_ik
         call advance_slow_dynamics(poly, cfg)
         if (mod(day, 30_ik) == 0_ik) call advance_boundary_dynamics(poly, cfg, .true., .false.)
      end do
      call check(poly%cohort%n == 1_ik, 'leafless with storage: survives a 150-day dry season')
      call check(poly%cohort%nonstructural_carbon(1) > 0.0_wp .and. poly%cohort%leaf_carbon(1) == 0.0_wp, &
                 'leafless with storage: it paid its upkeep from storage and stayed leafless')
      call check(poly%cohort%nplant(1) > 0.9_wp * n0,                                              &
                 'leafless with storage: only the growth-driven hazard thinned it')
   end subroutine check_exhausted

   !----- One leafless, dormant cohort holding `storage` [kgC/plant], with a day's stem respiration !
   !      and no photosynthesis on every step (there is no fast loop to overwrite the accumulators). -!
   subroutine leafless_stand(poly, cfg, storage)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: storage
      real(wp), parameter :: STEM_RESP_DAY = 1.0e-3_wp   !< [kgC/plant] one day's stem respiration
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.05_wp, 20.0_wp)
      call finalize_init(poly)
      poly%cohort%leaf_carbon(1)          = 0.0_wp
      poly%cohort%nonstructural_carbon(1) = storage
      call set_cohort_size_from_carbon(poly%cohort, 1_ik)
      poly%cohort%leaf_flush_tendency(1) = 0.0_wp      ! dormant: no flush, no turnover
      poly%cohort%gpp_accum(1)           = 0.0_wp
      poly%cohort%leaf_resp_accum(1)     = 0.0_wp
      poly%cohort%stem_resp_accum(1)     = STEM_RESP_DAY
      poly%cohort%root_resp_accum(1)     = 0.0_wp
   end subroutine leafless_stand

end program test_mortality_pathways
