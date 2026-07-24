!----- Carbon-driven growth (driver level): the growth-only allocation closure, the driver's -----!
!----- update_biomass_turnover shed/snap, and a carbon-mode step that grows wood_carbon -> dbh. ---!
program test_carbon_growth
   use meds_kinds,                  only : wp, ik
   use meds_config,                 only : meds_config_t
   use meds_core_interface,         only : site_t, carbon_flux_block
   use meds_plant_interface,        only : plant_carbon_allocation
   use meds_vegetation_dynamics,    only : update_biomass_turnover, advance_plant_traits,          &
                                            shed_turnover_water
   use meds_init,                   only : init_bare_ground, add_cohort, finalize_init
   use meds_stepper,                only : advance_one_step
   use meds_output_diagnostics,     only : has_nan
   use meds_test_support,           only : build_test_config, check, check_close, banner
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)        :: site
   real(wp)    :: gl, gf, gw, gs, gr, gresp, def, wood0, dbh0
   real(wp)    :: leaf_shed, root_shed
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

   !=== 2. update_biomass_turnover: current-pool decay, clamp, dormancy snap-to-bare. ==========!
   ! not dormant (flush high): leaf 0.1/d over pool 1 => 0.1; fine root 0.05/d over pool 0.5 => 0.025.
   call update_biomass_turnover(0.1_wp, 0.05_wp, 1.0e6_wp, 1.0_wp, 0.5_wp, 1.0_wp, 0.02_wp, 1.0_wp, &
        leaf_shed, root_shed)
   call check_close(leaf_shed, 0.1_wp,   1.0e-12_wp, 'turnover: leaf shed = rate*pool*dt')
   call check_close(root_shed, 0.025_wp, 1.0e-12_wp, 'turnover: fine-root shed = rate*pool*dt')
   ! clamp: rate*dt > 1 removes at most the whole pool.
   call update_biomass_turnover(10.0_wp, 0.0_wp, 1.0e6_wp, 0.4_wp, 0.0_wp, 1.0_wp, 0.02_wp, 1.0_wp, &
        leaf_shed, root_shed)
   call check_close(leaf_shed, 0.4_wp, 1.0e-12_wp, 'turnover: leaf shed clamped to pool')
   ! dormant (flush ~ 0) + residual below snap floor => snap the whole pool to bare.
   call update_biomass_turnover(0.1_wp, 0.0_wp, 0.0_wp, 0.021_wp, 0.0_wp, 1.0_wp, 0.02_wp, 1.0_wp,  &
        leaf_shed, root_shed)
   call check_close(leaf_shed, 0.021_wp, 1.0e-12_wp, 'turnover: dormant canopy snaps to bare')
   ! same pool but FLUSHING => no snap, just the decrement.
   call update_biomass_turnover(0.1_wp, 0.0_wp, 1.0_wp, 0.021_wp, 0.0_wp, 1.0_wp, 0.02_wp, 1.0_wp,  &
        leaf_shed, root_shed)
   call check_close(leaf_shed, 0.1_wp*0.021_wp, 1.0e-12_wp, 'turnover: no snap while flushing')

   !=== 3. A carbon-mode step grows wood_carbon -> dbh, and leaf_area stays leaf_carbon*sla. ==!
   cfg%gpp_ref = 0.5_wp
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 3_ik, 0.1_wp, 20.0_wp)       ! climax cohort
   call finalize_init(site)
   wood0 = site%cohort%wood_carbon(1) ; dbh0 = site%cohort%dbh(1)
   do istep = 1_ik, 30_ik
      call advance_one_step(site, cfg, .false., .false.)         ! growth + mortality only
   end do
   call check(site%cohort%wood_carbon(1) > wood0, 'carbon growth did not add wood_carbon')
   call check(site%cohort%dbh(1) > dbh0,          'dbh did not increase with wood_carbon (the flip)')
   call check_close(site%cohort%leaf_area(1), site%cohort%leaf_carbon(1) * cfg%pft%sla(3_ik),   &
                    1.0e-9_wp, 'carbon-mode leaf_area /= leaf_carbon*sla')
   call check(.not. has_nan(site), 'carbon-mode step produced NaN')

   !=== 4. Trait plasticity: a SHADED cohort acclimates -- SLA up, Vcmax down, leaf lifespan up. ==!
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
      ! capped at 0 => a long-lived PFT keeps ONE lifespan through the canopy (no shortening).
      call check_close(site%cohort%llspan(1), ll0, 1.0e-12_wp, 'long-lived PFT: leaf lifespan is canopy-invariant (capped)')
      ! SLA rose => the leaf pool overshot allometry; the excess was resorbed to storage.
      call check(site%cohort%nonstructural_carbon(1) > ns0, 'SLA overshoot resorbed to storage')
      call check_close(site%cohort%leaf_area(1), site%cohort%leaf_carbon(1) * site%cohort%sla(1),   &
                       1.0e-9_wp, 'leaf_area stays leaf_carbon*sla after resorption')
   end block

   !=== 5. shed_turnover_water (P4): a NET leaf/root LOSS sheds water in the SAME proportion as   !
   !    the carbon loss (remaining tissue's rwc unchanged), credited to the patch's shed_water_rate; !
   !    a NET GAIN touches nothing (P3's own mass-conserving seam handles that direction instead). ===!
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

      !----- (b) Net GAIN: touches nothing (mass conserved, P3's own seam handles the psi response). !
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

   write(*,'(a)') '   PASS'
end program test_carbon_growth
