! SPDX-License-Identifier: Apache-2.0
!----- Carbon-driven growth (driver level): the growth-only allocation closure and a carbon-mode --!
!----- step that grows wood_carbon -> dbh. (Leaf loss: leaf_turnover_step, test_plant_phenology.) --!
program test_carbon_growth
   use meds_kinds,                  only : wp, ik
   use meds_allometry, only : dbh_to_wai, sapwood_fraction, carbon_to_structure
   use meds_config,                 only : meds_config_t
   use meds_site_state_types, only : site_t, carbon_flux_block, cohort_deriv_alloc, cohort_bad_field,    &
                                     cohort_field_name
   use meds_demography_update, only : update_cohort_states
   use meds_plant_carbon_allocation, only : plant_carbon_allocation
   use meds_vegetation_dynamics,    only : advance_plant_traits, shed_turnover_water, reacclimate_plant_traits
   use meds_init,                   only : init_bare_ground, add_cohort, finalize_init
   use meds_stepper,                only : advance_one_step
   use meds_diagnostic_reduce,     only : has_nan
   use meds_test_support, only : banner, build_test_config, check, check_close
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)        :: site
   real(wp)    :: gl, gf, gw, gs, gr, gresp, def, wood0, dbh0
   logical     :: starv
   integer(ik) :: istep

   call banner('carbon-driven growth')
   cfg = build_test_config()

   !=== 1. Growth-only allocation closes: (growth pools + npp_store) - deficit = net - growth_resp.=!
   call plant_carbon_allocation(1.0_wp, 0.0_wp, 0.3_wp, 0.0_wp, 0.05_wp, 0.05_wp, 0.02_wp, 0.3_wp, &
        gl, gf, gw, gs, gr, gresp, def, starv)
   call check_close((gl+gf+gw+gr+gs) - def, 1.0_wp - gresp, 1.0e-9_wp,                              &
                    'growth allocation carbon closure')
   call check(gr > 0.0_wp, 'mature cohort allocates to reproduction')
   call check(gw > 0.0_wp, 'surplus carbon grows wood')
   call check(gresp > 0.0_wp, 'growth respiration charged on realized growth')

   !=== 2. A carbon-mode step grows wood_carbon -> dbh, and leaf_area stays leaf_carbon*sla. ==!
   cfg%gpp_ref = 0.5_wp
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 3_ik, 0.1_wp, 20.0_wp)       ! climax cohort
   call finalize_init(site)
   wood0 = site%cohort%wood_carbon(1) ; dbh0 = site%cohort%dbh(1)
   do istep = 1_ik, 30_ik
      call advance_one_step(site, cfg)                           ! growth + mortality only
   end do
   call check(site%cohort%wood_carbon(1) > wood0, 'carbon growth did not add wood_carbon')
   call check(site%cohort%dbh(1) > dbh0,          'dbh did not increase with wood_carbon (the flip)')
   call check_close(site%cohort%leaf_area(1), site%cohort%leaf_carbon(1) * cfg%pft%sla(3_ik),   &
                    1.0e-9_wp, 'carbon-mode leaf_area /= leaf_carbon*sla')
   call check(.not. has_nan(site), 'carbon-mode step produced NaN')

   !----- has_nan sees the fast loop's states too: a non-finite canopy air or leaf temperature is  !
   !      what the integrators commit when even their smallest step fails. ------------------------!
   block
      real(wp) :: keep
      if (allocated(site%patch%cas)) then
         keep = site%patch%cas(1)%can_shv
         site%patch%cas(1)%can_shv = ieee_value(keep, ieee_quiet_nan)
         call check(has_nan(site), 'has_nan missed a NaN canopy-air humidity')
         site%patch%cas(1)%can_shv = keep
      end if
      if (allocated(site%cohort%leaf_temp)) then
         keep = site%cohort%leaf_temp(1)
         site%cohort%leaf_temp(1) = ieee_value(keep, ieee_positive_inf)
         call check(has_nan(site), 'has_nan missed an infinite leaf temperature')
         site%cohort%leaf_temp(1) = keep
      end if
      call check(.not. has_nan(site), 'has_nan still true after the state was restored')
   end block

   !=== 2b. The cached WOOD geometry must describe TODAY's tree, not the one 30 steps ago. ======!
   !        wood_area / sapwood_carbon / sapwood_area are derived once per size change rather    !
   !        than recomputed every dt_fast, so the fast loop reads them in place. That is only     !
   !        safe while every path that moves dbh, basal_area or wood_carbon re-derives them,      !
   !        growth included. Assert it directly: a stale cache is invisible in a conservation     !
   !        ledger and shows up only as a slowly wrong tree.                                      !
   block
      real(wp) :: wa_fresh, fs_fresh
      wa_fresh = dbh_to_wai(site%cohort%dbh(1), 1.0_wp, site%cohort%p_wai_b1(1), site%cohort%p_wai_b2(1))
      fs_fresh = sapwood_fraction(site%cohort%dbh(1), site%cohort%p_sapwood_area_b1(1),          &
                                  site%cohort%p_sapwood_area_b2(1))
      call check(wa_fresh > 0.0_wp, 'wood-geometry fixture is degenerate (zero wood area)')
      call check_close(site%cohort%wood_area(1), wa_fresh, 1.0e-13_wp,                           &
                       'cached wood_area is stale after growth')
      call check_close(site%cohort%sapwood_carbon(1), fs_fresh * site%cohort%wood_carbon(1),     &
                       1.0e-13_wp, 'cached sapwood_carbon is stale after growth')
      call check_close(site%cohort%sapwood_area(1), fs_fresh * site%cohort%basal_area(1),        &
                       1.0e-13_wp, 'cached sapwood_area is stale after growth')
   end block

   !=== 3. Trait plasticity: a SHADED cohort acclimates -- SLA up, Vcmax down, leaf lifespan up. ==!
   block
      real(wp) :: sla0, vc0, ll0, ns0
      cfg%trait_plasticity_on = .true.
      site%cohort%overtopping_lai(1) = 4.0_wp                  ! deep shade
      sla0 = site%cohort%sla(1) ; vc0 = site%cohort%vcmax25(1) ; ll0 = site%cohort%llspan(1)
      ns0  = site%cohort%nonstructural_carbon(1)
      do istep = 1_ik, 200_ik                                  ! many steps so the slow relaxation shows
         call advance_plant_traits(site, cfg, 1.0_wp)
      end do
      call check(site%cohort%sla(1)     > sla0, 'shaded cohort: SLA acclimates upward')
      call check(site%cohort%vcmax25(1) < vc0,  'shaded cohort: Vcmax acclimates downward')
      ! PFT 3 (climax) has leaf lifespan ~3 yr > the ~2.6-yr ED2 crossover, so kplastic_llspan is
      ! capped at 0 (docs/science/plant_traits.md eq. (4)) => a long-lived PFT keeps ONE lifespan
      ! through the canopy (no shortening).
      call check_close(site%cohort%llspan(1), ll0, 1.0e-12_wp, 'long-lived PFT: leaf lifespan is canopy-invariant (capped)')
      ! SLA rose => the leaf pool overshot allometry; the excess was resorbed to storage.
      call check(site%cohort%nonstructural_carbon(1) > ns0, 'SLA overshoot resorbed to storage')
      call check_close(site%cohort%leaf_area(1), site%cohort%leaf_carbon(1) * site%cohort%sla(1),   &
                       1.0e-9_wp, 'leaf_area stays leaf_carbon*sla after resorption')
   end block

   !=== 4. shed_turnover_water: a NET leaf/root LOSS sheds water in the SAME proportion as the    !
   !    carbon loss (remaining tissue's rwc unchanged), credited to the patch's shed_water_rate;   !
   !    a NET GAIN touches nothing (the tissue keeps its water mass; its capacity grows). ========!
   block
      type(carbon_flux_block) :: npp
      real(wp) :: leaf_water0, wood_water0, leaf_carbon0, fineroot_carbon0, nplant0
      real(wp) :: expect_leaf_shed, expect_wood_shed

      call init_bare_ground(site, cfg, 1_ik)
      call add_cohort(site, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
      call finalize_init(site)
      site%cohort%leaf_water_mass(1) = 4.0_wp
      site%cohort%wood_water_mass(1) = 9.0_wp
      leaf_carbon0     = site%cohort%leaf_carbon(1)
      fineroot_carbon0 = site%cohort%fineroot_carbon(1)
      nplant0          = site%cohort%nplant(1)
      leaf_water0      = site%cohort%leaf_water_mass(1)
      wood_water0      = site%cohort%wood_water_mass(1)

      allocate(npp%leaf(1), npp%fineroot(1), npp%wood(1), npp%nonstructural(1))
      npp%wood = 0.0_wp ; npp%nonstructural = 0.0_wp

      !----- (a) Net LOSS on both leaf and fine root: shed proportional to the carbon fraction lost. --!
      npp%leaf(1)     = -0.25_wp * leaf_carbon0        ! sheds 25% of the leaf pool this step
      npp%fineroot(1) = -0.10_wp * fineroot_carbon0    ! sheds 10% of the fine-root pool
      call shed_turnover_water(site, cfg, npp)
      expect_leaf_shed = leaf_water0 * 0.25_wp
      expect_wood_shed = wood_water0 * 0.10_wp
      call check_close(site%cohort%leaf_water_mass(1), leaf_water0 - expect_leaf_shed, 1.0e-9_wp, &
                       'shed_turnover_water: leaf water sheds proportional to the carbon loss')
      call check_close(site%cohort%wood_water_mass(1), wood_water0 - expect_wood_shed, 1.0e-9_wp, &
                       'shed_turnover_water: wood/root water sheds proportional to the carbon loss')
      call check_close(site%patch%shed_water_rate(1),                                              &
                       nplant0 * (expect_leaf_shed + expect_wood_shed) / cfg%dt_slow, 1.0e-9_wp,    &
                       'shed_turnover_water: patch rate = nplant-weighted shed amount / dt_slow')

      !----- (b) Net GAIN: touches nothing (mass conserved; the larger capacity reads a lower rwc). !
      site%cohort%leaf_water_mass(1) = leaf_water0    ! reset
      site%cohort%wood_water_mass(1) = wood_water0
      npp%leaf(1)     = 0.01_wp * leaf_carbon0
      npp%fineroot(1) = 0.01_wp * fineroot_carbon0
      call shed_turnover_water(site, cfg, npp)
      call check_close(site%cohort%leaf_water_mass(1), leaf_water0, 1.0e-12_wp,                    &
                       'shed_turnover_water: net GAIN leaves leaf water mass untouched')
      call check_close(site%cohort%wood_water_mass(1), wood_water0, 1.0e-12_wp,                    &
                       'shed_turnover_water: net GAIN leaves wood water mass untouched')
      call check_close(site%patch%shed_water_rate(1), 0.0_wp, 1.0e-12_wp,                           &
                       'shed_turnover_water: net GAIN sets patch rate back to 0 (SET, not accumulated)')

      !----- (c) Full snap-to-bare (net loss == the entire pool): sheds to EXACTLY 0. ---------------!
      site%cohort%leaf_water_mass(1) = leaf_water0
      npp%leaf(1)     = -leaf_carbon0
      npp%fineroot(1) = 0.0_wp
      call shed_turnover_water(site, cfg, npp)
      call check_close(site%cohort%leaf_water_mass(1), 0.0_wp, 1.0e-9_wp,                          &
                       'shed_turnover_water: full leaf-carbon loss sheds ALL leaf water')
   end block

   !=== The applier advances only the carbon and derives the geometry from it: a cohort         !
   !    that sheds every leaf ends with leaf area exactly 0, never a round-off negative, and its     !
   !    dbh and AGB are the allometry of its committed wood carbon.  ===============================!
   block
      real(wp) :: dt, dbh_c, h_c, ba_c, agb_c, la_c
      dt = 1.0_wp / 365.0_wp
      call init_bare_ground(site, cfg, 1_ik)
      call add_cohort(site, cfg, 1_ik, 2_ik, 0.5_wp, 3.0_wp)
      call finalize_init(site)
      call cohort_deriv_alloc(site%deriv, site%cohort%n)
      !----- Shed a hair more leaf carbon than there is, as a rounded loss rate can. ------------!
      site%deriv%d_leaf_carbon_dt(1)          = -site%cohort%leaf_carbon(1) / dt * (1.0_wp + 1.0e-15_wp)
      site%deriv%d_fineroot_carbon_dt(1)      = 0.0_wp
      site%deriv%d_wood_carbon_dt(1)          = 0.01_wp
      site%deriv%d_nonstructural_carbon_dt(1) = 0.0_wp
      site%deriv%dln_nplant_dt(1)             = 0.0_wp
      call update_cohort_states(site%cohort, site%deriv, dt, cfg%negligible_nplant)
      call check(site%cohort%leaf_carbon(1) == 0.0_wp, 'applier: a full shed leaves no leaf carbon')
      call check(site%cohort%leaf_area(1) == 0.0_wp, 'applier: and exactly no leaf area')
      call carbon_to_structure(site%cohort%allom, site%cohort%wood_carbon(1), site%cohort%leaf_carbon(1),              &
                               site%cohort%p_wood_density(1), site%cohort%p_hgt_max(1),             &
                               site%cohort%p_aboveground_frac(1), site%cohort%sla(1),               &
                               dbh_c, h_c, ba_c, agb_c, la_c)
      call check(site%cohort%dbh(1) == dbh_c .and. site%cohort%height(1) == h_c                     &
                 .and. site%cohort%agb(1) == agb_c, 'applier: geometry is the allometry of the carbon')

      !----- An overdraw beyond round-off is an operator bug: the applier does not floor it (a   !
      !      floor would create the carbon it took), and the bounds guard names the pool. --------!
      call init_bare_ground(site, cfg, 1_ik)
      call add_cohort(site, cfg, 1_ik, 2_ik, 0.5_wp, 3.0_wp)
      call finalize_init(site)
      call cohort_deriv_alloc(site%deriv, site%cohort%n)
      site%deriv%d_leaf_carbon_dt(1)          = 0.0_wp
      site%deriv%d_fineroot_carbon_dt(1)      = 0.0_wp
      site%deriv%d_wood_carbon_dt(1)          = 0.0_wp
      site%deriv%d_nonstructural_carbon_dt(1) = -site%cohort%nonstructural_carbon(1) / dt * 1.01_wp
      site%deriv%dln_nplant_dt(1)             = 0.0_wp
      call check(cohort_bad_field(site%cohort, 1_ik) == 0_ik, 'applier: the cohort starts in bounds')
      call update_cohort_states(site%cohort, site%deriv, dt, cfg%negligible_nplant)
      call check(site%cohort%nonstructural_carbon(1) < 0.0_wp, 'applier: a 1 % overdraw is not floored')
      call check(trim(cohort_field_name(cohort_bad_field(site%cohort, 1_ik))) == 'nonstructural_carbon', &
                 'applier: the bounds guard names the overdrawn pool')
   end block

   !=== Re-acclimation to a lower SLA asks storage for more leaf carbon: it takes only what     !
   !    storage holds, so the cohort's carbon is unchanged and its leaf area falls short. =======!
   block
      real(wp) :: carbon0, leaf0, store0
      call init_bare_ground(site, cfg, 1_ik)
      call add_cohort(site, cfg, 1_ik, 2_ik, 0.5_wp, 3.0_wp)
      call finalize_init(site)
      site%cohort%sla(1) = 2.0_wp * cfg%pft%sla(2_ik)       ! the file's SLA is half the cohort's
      leaf0  = site%cohort%leaf_carbon(1)
      store0 = 0.25_wp * leaf0                              ! storage covers a quarter of the doubling
      site%cohort%nonstructural_carbon(1) = store0
      carbon0 = leaf0 + store0
      call reacclimate_plant_traits(site, cfg)
      call check(site%cohort%sla(1) == cfg%pft%sla(2_ik), 'reacclimation takes the file''s SLA')
      call check(site%cohort%nonstructural_carbon(1) == 0.0_wp, 'reacclimation empties storage, no further')
      call check_close(site%cohort%leaf_carbon(1) + site%cohort%nonstructural_carbon(1), carbon0,      &
                       1.0e-12_wp, 'reacclimation conserves the cohort''s leaf and storage carbon')
      call check(site%cohort%leaf_carbon(1) < 2.0_wp * leaf0, 'the leaves take only what storage has')
   end block

   write(*,'(a)') '   PASS'
end program test_carbon_growth
