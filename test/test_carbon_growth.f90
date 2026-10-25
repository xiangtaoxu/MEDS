! SPDX-License-Identifier: Apache-2.0
!----- Carbon-driven growth (driver level): the growth-only allocation closure and a carbon-mode --!
!----- step that grows wood_carbon -> dbh. (Leaf loss: leaf_turnover_step, test_plant_phenology.) --!
program test_carbon_growth
   use meds_kinds,                  only : wp, ik
   use meds_therm_lib,              only : internal_energy_liquid
   use meds_allometry_lib, only : dbh_to_wai, sapwood_fraction, carbon_to_structure
   use meds_config,                 only : meds_config_t, growth_window_steps
   use meds_cohort_state_vars, only : cohort_deriv_alloc, cohort_bad_field, cohort_field_finite, cohort_field_name
   use meds_patch_state_vars,  only : canopy_air_bad_field
   use meds_polygon_state_vars, only : polygon_t
   use meds_slow_update, only : update_cohort_states, terminate_cohorts, update_cohort_derivatives, commit_trait_carbon
   use meds_plant_carbon_allocation, only : plant_carbon_allocation
   use meds_slow_state, only : advance_plant_traits, shed_turnover_water, reacclimate_plant_traits,      &
                               slow_step_data_t
   use meds_init_bare_ground, only : init_bare_ground, add_cohort, finalize_init
   use meds_stepper,                only : advance_one_step
   use meds_test_support, only : banner, build_test_config, check, check_close
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use meds_slow_ledger,  only : slow_ledger_t
   implicit none
   real(wp), allocatable :: leaf_new(:), storage_new(:)   !< the traits' leaf and storage carbon
   type(slow_ledger_t) :: no_ledger    !< inactive: a test step books nothing
   real(wp)            :: shed_water   !< an operator's shed water, unread here

   type(meds_config_t) :: cfg
   type(polygon_t)     :: poly
   real(wp)    :: gl, gf, gw, gs, gr, gresp, def, ex, wood0, dbh0
   integer(ik) :: istep

   call banner('carbon-driven growth')
   cfg = build_test_config()

   !=== 1. Growth-only allocation closes: (growth pools + npp_store) - deficit = net - growth_resp.=!
   call plant_carbon_allocation(1.0_wp, 0.0_wp, 0.3_wp, 0.0_wp, 0.05_wp, 0.05_wp, 0.02_wp, 0.3_wp, &
        gl, gf, gw, gs, gr, gresp, def, huge(1.0_wp), ex)
   call check_close((gl+gf+gw+gr+gs) - def, 1.0_wp - gresp, 1.0e-9_wp,                              &
                    'growth allocation carbon closure')
   call check(gr > 0.0_wp, 'mature cohort allocates to reproduction')
   call check(gw > 0.0_wp, 'surplus carbon grows wood')
   call check(gresp > 0.0_wp, 'growth respiration charged on realized growth')

   !=== 2. A carbon-mode step grows wood_carbon -> dbh, and leaf_area stays leaf_carbon*sla. ==!
   cfg%gpp_ref = 0.5_wp
   call init_bare_ground(poly, cfg, 1_ik)
   call add_cohort(poly, cfg, 1_ik, 3_ik, 0.1_wp, 20.0_wp)       ! climax cohort
   call finalize_init(poly)
   wood0 = poly%cohort%wood_carbon(1) ; dbh0 = poly%cohort%dbh(1)
   do istep = 1_ik, 30_ik
      call advance_one_step(poly, cfg)                           ! growth + mortality only
   end do
   call check(poly%cohort%wood_carbon(1) > wood0, 'carbon growth did not add wood_carbon')
   call check(poly%cohort%dbh(1) > dbh0,          'dbh did not increase with wood_carbon (the flip)')
   call check_close(poly%cohort%leaf_area(1), poly%cohort%leaf_carbon(1) * cfg%pft%sla(3_ik),   &
                    1.0e-9_wp, 'carbon-mode leaf_area /= leaf_carbon*sla')
   call check(cohort_bad_field(poly%cohort, 1_ik) == 0_ik, 'carbon-mode step left the cohort out of bounds')

   !----- The driver's checks see the fast loop's states too: a non-finite canopy air or tissue     !
   !      energy is what the integrators commit when even their smallest step fails. ---------------!
   block
      real(wp)    :: keep
      integer(ik) :: k
      if (allocated(poly%patch%cas)) then
         keep = poly%patch%cas(1)%can_shv
         poly%patch%cas(1)%can_shv = ieee_value(keep, ieee_quiet_nan)
         call check(canopy_air_bad_field(poly%patch%cas(1)) == 2_ik, 'a NaN canopy-air humidity is caught')
         poly%patch%cas(1)%can_shv = keep
         call check(canopy_air_bad_field(poly%patch%cas(1)) == 0_ik, 'the restored canopy air passes')
      end if
      keep = poly%cohort%leaf_energy(1)
      poly%cohort%leaf_energy(1) = ieee_value(keep, ieee_positive_inf)
      k = cohort_bad_field(poly%cohort, 1_ik)
      call check(trim(cohort_field_name(k)) == 'leaf_energy' .and. .not. cohort_field_finite(poly%cohort, 1_ik, k), &
                 'an infinite leaf energy is caught, as not finite')
      poly%cohort%leaf_energy(1) = keep
      call check(cohort_bad_field(poly%cohort, 1_ik) == 0_ik, 'the restored cohort passes')
   end block

   !=== 2b. The cached WOOD geometry must describe TODAY's tree, not the one 30 steps ago. ======!
   !        wood_area / sapwood_carbon / sapwood_area are derived once per size change rather    !
   !        than recomputed every dt_fast, so the fast loop reads them in place. That is only     !
   !        safe while every path that moves dbh, basal_area or wood_carbon re-derives them,      !
   !        growth included. Assert it directly: a stale cache is invisible in a conservation     !
   !        ledger and shows up only as a slowly wrong tree.                                      !
   block
      real(wp) :: wa_fresh, fs_fresh
      associate (pft => poly%cohort%pft_table, pf => poly%cohort%pft(1))
         wa_fresh = dbh_to_wai(pft%allom(pf), poly%cohort%dbh(1), 1.0_wp)
         fs_fresh = sapwood_fraction(pft%allom(pf), poly%cohort%dbh(1))
      end associate
      call check(wa_fresh > 0.0_wp, 'wood-geometry fixture is degenerate (zero wood area)')
      call check_close(poly%cohort%wood_area(1), wa_fresh, 1.0e-13_wp,                           &
                       'cached wood_area is stale after growth')
      call check_close(poly%cohort%sapwood_carbon(1), fs_fresh * poly%cohort%wood_carbon(1),     &
                       1.0e-13_wp, 'cached sapwood_carbon is stale after growth')
      call check_close(poly%cohort%sapwood_area(1), fs_fresh * poly%cohort%basal_area(1),        &
                       1.0e-13_wp, 'cached sapwood_area is stale after growth')
   end block

   !=== 3. Trait plasticity: a SHADED cohort acclimates -- SLA up, Vcmax down, turnover down. =====!
   block
      real(wp) :: sla0, vc0, ll0, ns0
      cfg%trait_plasticity_on = .true.
      poly%cohort%overtopping_lai(1) = 4.0_wp                  ! deep shade
      sla0 = poly%cohort%sla(1) ; vc0 = poly%cohort%vcmax25(1) ; ll0 = poly%cohort%leaf_turnover_rate(1)
      ns0  = poly%cohort%nonstructural_carbon(1)
      do istep = 1_ik, 200_ik                                  ! many steps so the slow relaxation shows
         call advance_plant_traits(poly, 1.0_wp, leaf_new, storage_new)
         call commit_trait_carbon(poly, leaf_new, storage_new)
      end do
      call check(poly%cohort%sla(1)     > sla0, 'shaded cohort: SLA acclimates upward')
      call check(poly%cohort%vcmax25(1) < vc0,  'shaded cohort: Vcmax acclimates downward')
      ! PFT 3 (climax) turns its leaves over at 0.33 /yr, below the ~0.39 /yr ED2 crossover, so
      ! kplastic_leaf_turnover is capped at 0 (docs/science/plant_traits.md eq. (4)) => a
      ! slow-turnover PFT keeps ONE rate through the canopy (never faster in shade).
      call check_close(poly%cohort%leaf_turnover_rate(1), ll0, 1.0e-12_wp,                       &
                       'slow-turnover PFT: leaf turnover rate is canopy-invariant (capped)')
      ! SLA rose => the leaf pool overshot allometry; the excess was resorbed to storage.
      call check(poly%cohort%nonstructural_carbon(1) > ns0, 'SLA overshoot resorbed to storage')
      call check_close(poly%cohort%leaf_area(1), poly%cohort%leaf_carbon(1) * poly%cohort%sla(1),   &
                       1.0e-9_wp, 'leaf_area stays leaf_carbon*sla after resorption')
   end block

   !=== 4. shed_turnover_water: a NET leaf/root LOSS sheds water in the SAME proportion as the    !
   !    carbon loss (remaining tissue's rwc unchanged), credited to the patch's shed_water_rate;   !
   !    a NET GAIN touches nothing (the tissue keeps its water mass; its capacity grows). ========!
   block
      type(slow_step_data_t) :: step
      real(wp) :: leaf_water0, wood_water0, leaf_carbon0, fineroot_carbon0, nplant0
      real(wp) :: expect_leaf_shed, expect_wood_shed

      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
      call finalize_init(poly)
      poly%cohort%leaf_water_mass(1) = 4.0_wp
      poly%cohort%wood_water_mass(1) = 9.0_wp
      leaf_carbon0     = poly%cohort%leaf_carbon(1)
      fineroot_carbon0 = poly%cohort%fineroot_carbon(1)
      nplant0          = poly%cohort%nplant(1)
      leaf_water0      = poly%cohort%leaf_water_mass(1)
      wood_water0      = poly%cohort%wood_water_mass(1)

      allocate(step%d_leaf(1), step%d_fineroot(1), step%d_wood(1), step%d_nonstructural(1))
      step%d_wood = 0.0_wp ; step%d_nonstructural = 0.0_wp

      !----- (a) Net LOSS on both leaf and fine root: shed proportional to the carbon fraction lost. --!
      step%d_leaf(1)     = -0.25_wp * leaf_carbon0        ! sheds 25% of the leaf pool this step
      step%d_fineroot(1) = -0.10_wp * fineroot_carbon0    ! sheds 10% of the fine-root pool
      call shed_turnover_water(poly, cfg, step, no_ledger)
      expect_leaf_shed = leaf_water0 * 0.25_wp
      !----- The fine roots hold their carbon's share of the store they share with the sapwood, and   !
      !      shed a tenth of that: the store's water per unit of its carbon is unchanged once the     !
      !      fine roots have lost their tenth. ------------------------------------------------------!
      expect_wood_shed = wood_water0 * 0.10_wp * fineroot_carbon0                                    &
                       / (poly%cohort%sapwood_carbon(1) + fineroot_carbon0)
      call check_close(poly%cohort%leaf_water_mass(1), leaf_water0 - expect_leaf_shed, 1.0e-9_wp, &
                       'shed_turnover_water: leaf water sheds proportional to the carbon loss')
      call check_close(poly%cohort%wood_water_mass(1), wood_water0 - expect_wood_shed, 1.0e-9_wp, &
                       'shed_turnover_water: the fine roots shed their own share of the wood store')
      call check_close(poly%cohort%wood_water_mass(1)                                               &
                       / (poly%cohort%sapwood_carbon(1) + 0.9_wp * fineroot_carbon0),               &
                       wood_water0 / (poly%cohort%sapwood_carbon(1) + fineroot_carbon0), 1.0e-12_wp,   &
                       'shed_turnover_water: the wood store keeps its water per unit carbon')
      call check_close(poly%patch%shed_water_rate(1),                                              &
                       nplant0 * (expect_leaf_shed + expect_wood_shed) / cfg%dt_slow, 1.0e-9_wp,    &
                       'shed_turnover_water: patch rate = nplant-weighted shed amount / dt_slow')
      call check_close(poly%patch%shed_water_enth(1),                                              &
                       nplant0 * (expect_leaf_shed * internal_energy_liquid(poly%cohort%leaf_temp(1))  &
                                + expect_wood_shed * internal_energy_liquid(poly%cohort%wood_temp(1))) &
                       / cfg%dt_slow, 1.0e-9_wp,                                                   &
                       'shed_turnover_water: the water leaves with its enthalpy at its tissue''s temperature')

      !----- (b) Net GAIN: touches nothing (mass conserved; the larger capacity reads a lower rwc). !
      poly%cohort%leaf_water_mass(1) = leaf_water0    ! reset
      poly%cohort%wood_water_mass(1) = wood_water0
      step%d_leaf(1)     = 0.01_wp * leaf_carbon0
      step%d_fineroot(1) = 0.01_wp * fineroot_carbon0
      call shed_turnover_water(poly, cfg, step, no_ledger)
      call check_close(poly%cohort%leaf_water_mass(1), leaf_water0, 1.0e-12_wp,                    &
                       'shed_turnover_water: net GAIN leaves leaf water mass untouched')
      call check_close(poly%cohort%wood_water_mass(1), wood_water0, 1.0e-12_wp,                    &
                       'shed_turnover_water: net GAIN leaves wood water mass untouched')
      call check_close(poly%patch%shed_water_rate(1), 0.0_wp, 1.0e-12_wp,                           &
                       'shed_turnover_water: net GAIN sets patch rate back to 0 (SET, not accumulated)')

      !----- (c) Full snap-to-bare (net loss == the entire pool): sheds to EXACTLY 0. ---------------!
      poly%cohort%leaf_water_mass(1) = leaf_water0
      step%d_leaf(1)     = -leaf_carbon0
      step%d_fineroot(1) = 0.0_wp
      call shed_turnover_water(poly, cfg, step, no_ledger)
      call check_close(poly%cohort%leaf_water_mass(1), 0.0_wp, 1.0e-9_wp,                          &
                       'shed_turnover_water: full leaf-carbon loss sheds ALL leaf water')
   end block

   !=== The applier advances only the carbon and derives the geometry from it: a cohort         !
   !    that sheds every leaf ends with leaf area exactly 0, never a round-off negative, and its     !
   !    dbh and AGB are the allometry of its committed wood carbon.  ===============================!
   block
      real(wp) :: dt, dbh_c, h_c, ba_c, agb_c, la_c
      dt = 1.0_wp / 365.0_wp
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.5_wp, 3.0_wp)
      call finalize_init(poly)
      call cohort_deriv_alloc(poly%deriv, poly%cohort%n)
      !----- Shed a hair more leaf carbon than there is, as a rounded loss rate can. ------------!
      poly%deriv%d_leaf_carbon_dt(1)          = -poly%cohort%leaf_carbon(1) / dt * (1.0_wp + 1.0e-15_wp)
      poly%deriv%d_fineroot_carbon_dt(1)      = 0.0_wp
      poly%deriv%d_wood_carbon_dt(1)          = 0.01_wp
      poly%deriv%d_nonstructural_carbon_dt(1) = 0.0_wp
      poly%deriv%dln_nplant_dt(1)             = 0.0_wp
      call update_cohort_states(poly%cohort, poly%deriv, dt, cfg%negligible_nplant)
      call check(poly%cohort%leaf_carbon(1) == 0.0_wp, 'applier: a full shed leaves no leaf carbon')
      call check(poly%cohort%leaf_area(1) == 0.0_wp, 'applier: and exactly no leaf area')
      call carbon_to_structure(poly%cohort%pft_table%allom(poly%cohort%pft(1)), poly%cohort%wood_carbon(1),        &
                               poly%cohort%leaf_carbon(1), cfg%pft%wood_density(poly%cohort%pft(1)),     &
                               poly%cohort%sla(1), dbh_c, h_c, ba_c, agb_c, la_c)
      call check(poly%cohort%dbh(1) == dbh_c .and. poly%cohort%height(1) == h_c                     &
                 .and. poly%cohort%agb(1) == agb_c, 'applier: geometry is the allometry of the carbon')

      !----- An overdraw beyond round-off is an operator bug: the applier does not floor it (a   !
      !      floor would create the carbon it took), and the bounds guard names the pool. --------!
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.5_wp, 3.0_wp)
      call finalize_init(poly)
      call cohort_deriv_alloc(poly%deriv, poly%cohort%n)
      poly%deriv%d_leaf_carbon_dt(1)          = 0.0_wp
      poly%deriv%d_fineroot_carbon_dt(1)      = 0.0_wp
      poly%deriv%d_wood_carbon_dt(1)          = 0.0_wp
      poly%deriv%d_nonstructural_carbon_dt(1) = -poly%cohort%nonstructural_carbon(1) / dt * 1.01_wp
      poly%deriv%dln_nplant_dt(1)             = 0.0_wp
      call check(cohort_bad_field(poly%cohort, 1_ik) == 0_ik, 'applier: the cohort starts in bounds')
      call update_cohort_states(poly%cohort, poly%deriv, dt, cfg%negligible_nplant)
      call check(poly%cohort%nonstructural_carbon(1) < 0.0_wp, 'applier: a 1 % overdraw is not floored')
      call check(trim(cohort_field_name(cohort_bad_field(poly%cohort, 1_ik))) == 'nonstructural_carbon', &
                 'applier: the bounds guard names the overdrawn pool')
   end block

   !=== The driver's tendency follows the applier's rule: a step that empties a pool lands it on  !
   !    0, and an overdraw reaches the applier as an overdraw, which the bounds guard names. A     !
   !    floor on the driver's tentative pool would create the carbon before the guard could see it. =!
   block
      type(slow_step_data_t) :: step
      real(wp) :: dt
      dt = 1.0_wp / 365.0_wp
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.5_wp, 3.0_wp)
      call finalize_init(poly)
      allocate(step%d_leaf(1), step%d_fineroot(1), step%d_wood(1), step%d_nonstructural(1))
      step%d_leaf(1)          = -poly%cohort%leaf_carbon(1) * (1.0_wp + 1.0e-15_wp)   ! a rounded full shed
      step%d_fineroot(1)      = 0.0_wp
      step%d_wood(1)          = 0.0_wp
      step%d_nonstructural(1) = -poly%cohort%nonstructural_carbon(1) * 1.01_wp        ! a 1 % overdraw
      call cohort_deriv_alloc(poly%deriv, poly%cohort%n)
      step%mortality = [0.0_wp]
      call update_cohort_derivatives(poly, step, dt, growth_window_steps(cfg), 1_ik)
      call update_cohort_states(poly%cohort, poly%deriv, dt, cfg%negligible_nplant)
      call check(poly%cohort%leaf_carbon(1) == 0.0_wp, 'driver + applier: a rounded full shed lands on 0')
      call check(poly%cohort%nonstructural_carbon(1) < 0.0_wp,                                     &
                 'driver + applier: a 1 % overdraw reaches the applier, not floored on the way')
      call check(trim(cohort_field_name(cohort_bad_field(poly%cohort, 1_ik))) == 'nonstructural_carbon', &
                 'driver + applier: the bounds guard names the overdrawn pool')
   end block

   !=== Re-acclimation to a lower SLA asks storage for more leaf carbon: it takes only what     !
   !    storage holds, so the cohort's carbon is unchanged and its leaf area falls short. =======!
   block
      real(wp) :: carbon0, leaf0, store0
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 0.5_wp, 3.0_wp)
      call finalize_init(poly)
      poly%cohort%sla(1) = 2.0_wp * cfg%pft%sla(2_ik)       ! the file's SLA is half the cohort's
      leaf0  = poly%cohort%leaf_carbon(1)
      store0 = 0.25_wp * leaf0                              ! storage covers a quarter of the doubling
      poly%cohort%nonstructural_carbon(1) = store0
      carbon0 = leaf0 + store0
      call reacclimate_plant_traits(poly, cfg, leaf_new, storage_new)
      call commit_trait_carbon(poly, leaf_new, storage_new)
      call check(poly%cohort%sla(1) == cfg%pft%sla(2_ik), 'reacclimation takes the file''s SLA')
      call check(poly%cohort%nonstructural_carbon(1) == 0.0_wp, 'reacclimation empties storage, no further')
      call check_close(poly%cohort%leaf_carbon(1) + poly%cohort%nonstructural_carbon(1), carbon0,      &
                       1.0e-12_wp, 'reacclimation conserves the cohort''s leaf and storage carbon')
      call check(poly%cohort%leaf_carbon(1) < 2.0_wp * leaf0, 'the leaves take only what storage has')
   end block

   !=== A density that declines to the floor ends ON it, and the cull removes a cohort there. ====!
   !    nplant * exp(log(floor/nplant)) is the floor only to an epsilon either side; a cohort an   !
   !    epsilon above the floor would never be culled, and would decline no further. ===========!
   block
      real(wp)    :: dt, n0
      integer(ik) :: k, n_on, n_cut
      dt = 1.0_wp / 365.0_wp
      n_on = 0_ik
      do k = 1_ik, 200_ik
         !----- 200 starting densities across six decades, each aimed exactly at the floor. -------!
         n0 = cfg%negligible_nplant * 10.0_wp ** (6.0_wp * real(k, wp) / 200.0_wp)
         call init_bare_ground(poly, cfg, 1_ik)
         call add_cohort(poly, cfg, 1_ik, 2_ik, n0, 3.0_wp)
         call finalize_init(poly)
         call cohort_deriv_alloc(poly%deriv, poly%cohort%n)
         poly%deriv%d_leaf_carbon_dt(1)          = 0.0_wp
         poly%deriv%d_fineroot_carbon_dt(1)      = 0.0_wp
         poly%deriv%d_wood_carbon_dt(1)          = 0.0_wp
         poly%deriv%d_nonstructural_carbon_dt(1) = 0.0_wp
         poly%deriv%dln_nplant_dt(1)             = -(log(n0 / cfg%negligible_nplant) / dt)
         call update_cohort_states(poly%cohort, poly%deriv, dt, cfg%negligible_nplant)
         if (poly%cohort%nplant(1) == cfg%negligible_nplant) n_on = n_on + 1_ik
      end do
      call check(n_on == 200_ik, 'applier: a decline aimed at the floor ends exactly on it')

      !----- Two cohorts of a tracked size: one at the floor (culled), one a hair above (kept). ---!
      call init_bare_ground(poly, cfg, 1_ik)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 1.0_wp, 60.0_wp)
      call add_cohort(poly, cfg, 1_ik, 2_ik, 1.0_wp, 61.0_wp)
      call finalize_init(poly)
      poly%cohort%nplant(1) = cfg%negligible_nplant
      poly%cohort%nplant(2) = nearest(cfg%negligible_nplant, 1.0_wp)
      call check(all(poly%cohort%nplant(1:2) * poly%cohort%agb(1:2) >= cfg%min_cohort_agb),        &
                 'cull fixture: both cohorts clear the AGB floor, so only the density floor acts')
      n_cut = poly%cohort%n
      call terminate_cohorts(poly, cfg, shed_water)
      call check(poly%cohort%n == n_cut - 1_ik, 'cull: exactly one cohort removed')
      call check(poly%cohort%nplant(1) > cfg%negligible_nplant, 'cull: the cohort AT the floor is the one removed')
   end block

   write(*,'(a)') '   PASS'
end program test_carbon_growth
