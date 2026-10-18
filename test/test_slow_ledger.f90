! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_slow_ledger -- the site conservation ledger's ARITHMETIC and its STORE COVERAGE.       !
!                                                                                          !
! The ledger reports numbers nobody can check by inspection, so what it must never do is be   !
! quietly BLIND. Two ways it can be, and one test for each:                                   !
!                                                                                          !
!   1. THE STORE MISSES A RESERVOIR. If slow_site_store forgets a field, every operator that   !
!      moves that field reads as perfectly conservative -- the ledger prints zeros and the      !
!      defect it exists to find is the one thing it cannot see. Fortran cannot enumerate a      !
!      derived type's components, so this is done the way test_state_combinators does it:      !
!      perturb each reservoir in turn and assert the total MOVES, by the right amount and in    !
!      the right currency. Adding a state field means adding a row here.                        !
!                                                                                          !
!   2. THE DECLARATION ARITHMETIC IS WRONG. A declared term is subtracted from the residual, so  !
!      a sign error would let a real leak cancel against a real source. Asserted directly.       !
!                                                                                          !
! The third test binds the canopy-air open-volume terms to the store they are supposed to        !
! explain: resize the control volume, and the three declared terms must equal the three store    !
! changes EXACTLY. That is the one declaration in the skeleton that does real arithmetic rather   !
! than passing a quantity through, and a production run closes it to 1e-17.                      !
!==========================================================================================!
program test_slow_ledger
   use meds_kinds,              only : wp, ik
   use meds_config,             only : meds_config_t
   use meds_site_state_vars,    only : site_t
   use meds_patch_state_vars,   only : cas_set_depth, patch_physical_state
   use meds_state_defaults,     only : UNSET_POSITIVE
   use meds_soil_column_config, only : soil_params_t, build_soil_hydr_params, soil_column_layers,     &
                                       soil_column_params, soil_column_start_theta, top_soil_psi
   use meds_init,               only : init_bare_ground, add_cohort, finalize_init
   use meds_cohort_state_vars,  only : init_cohort, set_cohort_size_from_carbon
   use meds_demography_cohort_fusefiss, only : terminate_cohorts
   use meds_demography_patch_fusefiss,  only : fuse_2_patches
   use meds_vegetation_dynamics,        only : vegetation_dynamics, restructure_stand, shed_turnover_water, pool_change_block
   use meds_biogeochem_types,           only : litter_input_t
   use meds_slow_ledger,        only : slow_store_t, slow_ledger_t, slow_site_store,               &
                                       slow_ledger_open, slow_ledger_declare, slow_ledger_mark,    &
                                       slow_fast_carbon_handover, slow_tissue_heat, SLOW_PHASE_GROW, &
                                       SLOW_PHASE_ALLOCATE, SLOW_PHASE_RECRUIT, SLOW_PHASE_TISSUE_WATER
   use meds_slow_dynamics,      only : settle_tissue_water
   use meds_config,             only : pft_tissue_curves
   use meds_water_retention,    only : water_curve_t, psi_from_water_content
   use meds_constants,          only : grav_head, n_soil_layer_max
   use meds_test_support, only : banner, build_test_config, check, check_close
   implicit none
   real(wp)    :: soil_bottoms(n_soil_layer_max)   ! [m] the column's layer bottoms, as the ledger builds them
   integer(ik) :: n_soil_layer

   real(wp), parameter :: RHO = 1.2_wp        !< the density the canopy-air store is valued at
   type(meds_config_t) :: cfg
   type(site_t)        :: site
   type(soil_params_t) :: soil
   type(slow_store_t)  :: s0, s1
   type(slow_ledger_t) :: led
   real(wp)            :: de, dw, dc, expect, hand

   call banner('slow ledger: store coverage and declaration arithmetic')

   cfg = build_test_config()
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 2_ik, 0.3_wp, 18.0_wp)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.8_wp,  6.0_wp)
   call finalize_init(site)
   call check(site%cohort%n == 2_ik .and. site%patch%n == 1_ik, 'fixture: 2 cohorts in 1 patch')

   associate (sc => cfg%soil_column)
      call soil_column_layers(sc, soil_bottoms, n_soil_layer)
      call build_soil_hydr_params(n_soil_layer, sc%retention, soil_bottoms, sc%theta_sat, sc%theta_res, &
                                  sc%ksat, sc%curve_par_a, sc%curve_par_n, sc%psi_fc, soil)
   end associate

   !----- Give every reservoir a non-zero value, so a perturbation below is a change to a store !
   !      that already exists rather than the creation of one. --------------------------------!
   site%patch%soil_water(1)%theta(1:soil%n_active)      = 0.30_wp
   site%patch%soil_energy(1)%soil_energy(1:soil%n_active) = 5.0e7_wp
   site%patch%soil_water(1)%w_surface      = 2.0_wp
   site%patch%soil_water(1)%w_surface_energy = 8.0e5_wp
   site%patch%snow(1)%swe(1)           = 5.0_wp
   site%patch%snow(1)%snow_energy(1)   = -1.0e6_wp
   site%patch%cas(1)%can_shv      = 0.008_wp
   site%patch%cas(1)%can_enthalpy = 2.9e5_wp
   site%patch%cas(1)%can_co2      = 400.0_wp
   site%patch%cas(1)%can_depth    = 20.0_wp
   site%patch%soil_carbon(1)%slow_carbon = 3.0_wp
   site%patch%recruit_pool(:,1)          = 0.0_wp
   site%cohort%leaf_water_mass(1:2) = 0.02_wp
   site%cohort%wood_water_mass(1:2) = 0.05_wp
   site%cohort%leaf_surf_water(1:2) = 0.01_wp
   site%cohort%wood_surf_water(1:2) = 0.03_wp

   !=== 1. STORE COVERAGE: every reservoir must reach the total, in its own currency. =========!
   call moves_water ('soil theta',        perturb_theta,        0.01_wp * soil%dz(1) * 1000.0_wp)
   call moves_water ('pond w_surface',    perturb_pond,         1.0_wp)
   call moves_water ('snow swe',          perturb_swe,          1.0_wp)
   call moves_water ('canopy vapour',     perturb_shv,          RHO * 20.0_wp * 0.001_wp)
   call moves_water ('tissue water',      perturb_tissue_water, site%cohort%nplant(1) * 0.01_wp)
   call moves_water ('interception film', perturb_film,         0.01_wp)
   call moves_energy('soil energy',       perturb_soil_energy,  1.0e6_wp * soil%dz(1))
   call moves_energy('pond energy',       perturb_pond_energy,  1.0e5_wp)
   call moves_energy('snow energy',       perturb_snow_energy,  1.0e5_wp)
   call moves_energy('canopy enthalpy',   perturb_cas_enth,     RHO * 20.0_wp * 1.0e4_wp)
   call moves_carbon('live pools',        perturb_live_c,       site%cohort%nplant(1) * 0.5_wp)
   call moves_carbon('CENTURY pools',     perturb_necromass,    1.0_wp)
   !----- Tissue HEAT and canopy CO2 have no closed-form expectation worth restating here (they  !
   !      are the store's own formulae); assert only that they reach the total, which is the     !
   !      blindness this test is about. ----------------------------------------------------!
   call moves('tissue heat',   perturb_tissue_temp, 'energy')
   call moves('canopy CO2',    perturb_cas_co2,     'carbon')
   call moves('recruit pool',  perturb_recruit,     'carbon')

   !=== 2. DECLARATION ARITHMETIC: a declared term is subtracted from the residual. ===========!
   led%active = .true.
   call slow_ledger_open(led, site, cfg, RHO)
   call slow_ledger_declare(led, carbon_in = 1.0_wp, water_out = 2.0_wp, energy_in = 3.0_wp)
   call slow_ledger_mark(led, site, cfg, SLOW_PHASE_GROW)   ! nothing moved: residual = -(in - out)
   call check_close(led%resid_sum(SLOW_PHASE_GROW)%carbon, -1.0_wp, 1.0e-12_wp, 'declared carbon_in is subtracted')
   call check_close(led%resid_sum(SLOW_PHASE_GROW)%water,   2.0_wp, 1.0e-12_wp, 'declared water_out is added back')
   call check_close(led%resid_sum(SLOW_PHASE_GROW)%energy, -3.0_wp, 1.0e-12_wp, 'declared energy_in is subtracted')

   !----- An UNDECLARED, UNCHANGED step must be exactly zero, or every phase carries a bias. ----!
   led = slow_ledger_t()
   led%active = .true.
   call slow_ledger_open(led, site, cfg, RHO)
   call slow_ledger_mark(led, site, cfg, SLOW_PHASE_GROW)
   call check(led%resid_sum(SLOW_PHASE_GROW)%carbon == 0.0_wp .and.                                &
              led%resid_sum(SLOW_PHASE_GROW)%water  == 0.0_wp .and.                                &
              led%resid_sum(SLOW_PHASE_GROW)%energy == 0.0_wp, 'a step that changes nothing closes at exactly zero')

   !=== 3. The canopy-air OPEN-VOLUME terms must equal the store change they explain. ==========!
   s0 = slow_site_store(site, cfg, soil, RHO)
   call cas_set_depth(site%patch%cas(1), 8.0_wp, rho_air=RHO, de_open=de, dw_open=dw, dc_open=dc)
   s1 = slow_site_store(site, cfg, soil, RHO)
   call check_close(s1%energy - s0%energy, de,             1.0e-6_wp, 'de_open explains the canopy energy change')
   call check_close(s1%water  - s0%water,  dw,             1.0e-12_wp, 'dw_open explains the canopy vapour change')
   call check_close(s1%carbon - s0%carbon, dc * 1.2e-8_wp, 1.0e-14_wp, 'dc_open explains the canopy CO2 change')
   call check(de < 0.0_wp, 'shrinking the canopy DETRAINS (negative)')

   !=== 4. The fast->slow handover is nplant- and area-weighted. ==============================!
   cfg%fast_biophysics_on = .true.
   site%cohort%gpp_accum(1:2)       = 1.0e-3_wp
   site%cohort%leaf_resp_accum(1:2) = 2.0e-4_wp
   site%cohort%stem_resp_accum(1:2) = 1.0e-4_wp
   site%cohort%root_resp_accum(1:2) = 1.0e-4_wp
   expect = site%patch%area(1) * (site%cohort%nplant(1) + site%cohort%nplant(2)) * 6.0e-4_wp
   hand   = slow_fast_carbon_handover(site, cfg)
   call check_close(hand, expect, 1.0e-14_wp, 'handover = area * sum(nplant * (gpp - maintenance resp))')

   !=== 5. A CULL hands its water to the ground channel, and the store sees its heat go. =======!
   !       The assertion is the TRANSFER, not the disappearance: what leaves the cohort must     !
   !       arrive in the patch's shed channel, and the tissue heat store must fall by what the   !
   !       cohort carried.                                                                       !
   call check_cull()

   !=== 6. A recruit is born at its PATCH's temperature, not at the global constant. ===========!
   call check_birth_temp()

   !=== 7. The allocator's outputs actually LEAVE the slow tier. ==============================!
   !       Neither ledger can catch a break here: the slow ledger declares the HANDOFF, so it      !
   !       closes whether or not the fast loop ever picks the rate up, and the fast CAS ledger      !
   !       closes around whatever nee_biotic it is given. Only a test binds the two ends, so these  !
   !       assert that the trait CHANGES the routed quantity -- the same differential shape the      !
   !       PFT-trait tests use, because a silently-zero channel is the realistic failure.            !
   call check_routing()

   !=== 8. Tissue THERMAL MASS tracks biomass, and is the whole of the store's tissue term. ===!
   call check_tissue_thermal_mass()

   !=== 9. Mortality's water LEAVES the tissue and ARRIVES in the shed channel, exactly. ======!
   !       It guards two errors every other assertion in this file passes through: not routing the !
   !       water at all, and scaling the loss with the SURVIVORS' density instead of the density   !
   !       DROP. It is the one property that pins the basis.                                        !
   call check_mortality_water()
   call check_litter_patch_lockstep()
   call check_litter_without_soil_carbon()

   !=== 10. Reproduction carbon reaches the recruit pool, every step and in proportion. =======!
   call check_recruit_pool()

   !=== 11. The turnover shed takes its water's HEAT with it, and declares it. =================!
   !        The leaf's heat capacity counts its water and the wood's counts the sapwood share of   !
   !        its store, so shedding either moves the tissue heat store; the declared energy must    !
   !        explain that move exactly.                                                              !
   call check_shed_heat()
   !=== 12. Tissue water is booked when it is made, shed or moved: a recruit's birth water, a    !
   !        canopy's snap to bare, and water above a store's capacity. ======================!
   call check_tissue_water_ledger()

   print '(a)', 'test_slow_ledger: ALL PASSED'

contains

   !----- Every water amount a slow step moves into or out of the tissue is declared, so the       !
   !       ledger's water residual is zero across each: (a) a recruitment, whose recruits are born   !
   !       holding their patch's top-soil potential, drawn from outside like their carbon; (b) a     !
   !       canopy that snaps to bare, whose leaf water goes to the ground with the leaves; (c) a      !
   !       store left above its capacity, whose excess the settle before the fast loop sheds to the   !
   !       ground.                                                                                   !
   subroutine check_tissue_water_ledger()
      type(meds_config_t) :: c
      type(site_t)        :: st
      type(slow_ledger_t) :: lg
      type(water_curve_t) :: leaf_curve, wood_curve
      real(wp) :: shed0, excess, n_new, tol, theta_start(n_soil_layer_max)
      integer(ik) :: p, ph(3)
      type(soil_params_t) :: sp

      c = build_test_config()
      tol = 1.0e-12_wp

      !----- (a) A recruitment: water in, declared, and the recruit born in its patch's state -- here  !
      !      a patch the fast loop has not seeded, so the run's starting state (check_birth_temp      !
      !      covers a seeded one). -------------------------------------------------------------------!
      call init_bare_ground(st, c, 1_ik)
      st%patch%recruit_pool(:, 1) = 0.0_wp
      st%patch%recruit_pool(2, 1) = 5.0_wp * c%min_recruit_size
      lg = slow_ledger_t() ; lg%active = .true.
      call slow_ledger_open(lg, st, c, RHO)
      call restructure_stand(st, c, .true., .false., lg)
      call check(st%cohort%n == 1_ik, 'tissue water (a): one recruit is born')
      call check(lg%decl_gross(SLOW_PHASE_RECRUIT)%water > 0.0_wp,                                 &
                 'tissue water (a): the recruit''s birth water is declared as drawn in')
      call check(abs(lg%resid_sum(SLOW_PHASE_RECRUIT)%water) <= tol * lg%decl_gross(SLOW_PHASE_RECRUIT)%water, &
                 'tissue water (a): the recruit phase''s water residual is zero')
      call pft_tissue_curves(c%hydraulics, c%pft, 2_ik, leaf_curve, wood_curve)
      call soil_column_params(c%soil_column, sp)
      call soil_column_start_theta(c%soil_column, sp, c%init_soil_theta, theta_start)
      call check_close(psi_from_water_content(st%cohort%leaf_water_mass(1), leaf_curve,            &
                       st%cohort%leaf_carbon(1)), top_soil_psi(sp, theta_start(1)), 1.0e-9_wp,     &
                       'tissue water (a): a recruit''s leaves are born at the starting top-soil potential')
      call check_close(st%cohort%leaf_temp(1), c%init_soil_temp, 1.0e-12_wp,                       &
                       'tissue water (a): a recruit is born at the starting temperature')

      !----- (b) A dormant canopy sheds every leaf in one step, its water with it. ------------------!
      c%fast_biophysics_on = .true.
      c%demography_on      = .false.
      c%pft%pheno_shed_rate_max(:) = 1.0_wp          ! [1/day]: the whole canopy in one day
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 20.0_wp)
      call finalize_init(st)
      st%cohort%leaf_flush_tendency(1) = 0.0_wp
      st%cohort%leaf_shed_tendency(1)  = 1.0_wp
      st%cohort%gpp_accum(1)       = 0.0_wp
      st%cohort%leaf_resp_accum(1) = 0.0_wp
      st%cohort%stem_resp_accum(1) = 0.0_wp
      st%cohort%root_resp_accum(1) = 0.0_wp
      lg = slow_ledger_t() ; lg%active = .true.
      call slow_ledger_open(lg, st, c, RHO)
      call vegetation_dynamics(st, c, ledger=lg)
      call settle_tissue_water(st, c, lg)
      call check(st%cohort%leaf_carbon(1) == 0.0_wp .and. st%cohort%leaf_water_mass(1) == 0.0_wp, &
                 'tissue water (b): the canopy went bare, and its leaves hold no water')
      ph = [SLOW_PHASE_ALLOCATE, SLOW_PHASE_GROW, SLOW_PHASE_TISSUE_WATER]
      do p = 1_ik, 3_ik
         call check(abs(lg%resid_sum(ph(p))%water) <= tol * max(lg%decl_gross(ph(p))%water, 1.0_wp), &
                    'tissue water (b): the snap to bare closes the water ledger, phase by phase')
      end do
      call check(lg%decl_gross(SLOW_PHASE_ALLOCATE)%water > 0.0_wp,                                &
                 'tissue water (b): the shed leaves'' water is declared as leaving')

      !----- (c) Leaf carbon taken without its water: the settle sheds the excess, declared. -----!
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 20.0_wp)
      call finalize_init(st)
      call pft_tissue_curves(c%hydraulics, c%pft, 1_ik, leaf_curve, wood_curve)
      st%cohort%leaf_carbon(1) = 0.5_wp * st%cohort%leaf_carbon(1)
      call set_cohort_size_from_carbon(st%cohort, 1_ik)
      excess = st%cohort%leaf_water_mass(1) - leaf_curve%water_sat * st%cohort%leaf_carbon(1)
      call check(excess > 0.0_wp, 'tissue water (c): the fixture leaves the leaves above their capacity')
      shed0 = shed_total(st, c)
      n_new = st%cohort%nplant(1)
      lg = slow_ledger_t() ; lg%active = .true.            ! no window yet: the settle takes its reference
      call settle_tissue_water(st, c, lg)
      call check_close(lg%decl_gross(SLOW_PHASE_TISSUE_WATER)%water, st%patch%area(1) * n_new * excess, &
                       1.0e-12_wp, 'tissue water (c): the excess is declared as leaving')
      call check_close(shed_total(st, c) - shed0, st%patch%area(1) * n_new * excess, 1.0e-12_wp,    &
                       'tissue water (c): and it arrives in the patch shed channel')
      call check(abs(lg%resid_sum(SLOW_PHASE_TISSUE_WATER)%water) <= tol * lg%decl_gross(SLOW_PHASE_TISSUE_WATER)%water, &
                 'tissue water (c): the settle phase''s water residual is zero')
      call check(abs(lg%resid_sum(SLOW_PHASE_TISSUE_WATER)%energy) <= 1.0e-12_wp                   &
                 * lg%decl_gross(SLOW_PHASE_TISSUE_WATER)%energy, 'tissue water (c): and its energy residual')

      !----- (d) A bare cohort re-flushes from storage: its new leaves take their water from the  !
      !      wood, at the potential the two share at rest, so the leaves start where the stem is   !
      !      (less the crown's gravity head) rather than at the dry end of their curve, and the    !
      !      plant's water is unchanged. --------------------------------------------------------!
      call check_reflush_water()
   end subroutine check_tissue_water_ledger

   subroutine check_reflush_water()
      type(meds_config_t) :: c
      type(site_t)        :: st
      type(slow_ledger_t) :: lg
      type(water_curve_t) :: leaf_curve, wood_curve
      real(wp) :: w_plant0, w_plant1, psi_leaf, psi_wood
      c = build_test_config()
      c%fast_biophysics_on = .true.
      c%demography_on      = .false.
      c%pft%fineroot_turnover_rate(:) = 0.0_wp      ! no fine-root shed: only the fill moves water
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 20.0_wp)
      call finalize_init(st)
      call pft_tissue_curves(c%hydraulics, c%pft, 1_ik, leaf_curve, wood_curve)
      st%cohort%leaf_carbon(1)     = 0.0_wp           ! bare, with its storage, at a drier stem
      st%cohort%leaf_water_mass(1) = 0.0_wp
      call set_cohort_size_from_carbon(st%cohort, 1_ik)
      st%cohort%wood_water_mass(1) = 0.9_wp * st%cohort%wood_water_mass(1)
      st%cohort%leaf_flush_tendency(1) = 1.0_wp
      st%cohort%leaf_shed_tendency(1)  = 0.0_wp
      st%cohort%gpp_accum(1)       = 0.0_wp
      st%cohort%leaf_resp_accum(1) = 0.0_wp
      st%cohort%stem_resp_accum(1) = 0.0_wp
      st%cohort%root_resp_accum(1) = 0.0_wp
      w_plant0 = st%cohort%leaf_water_mass(1) + st%cohort%wood_water_mass(1)
      lg = slow_ledger_t() ; lg%active = .true.
      call slow_ledger_open(lg, st, c, RHO)
      call vegetation_dynamics(st, c, ledger=lg)
      call check(st%cohort%leaf_carbon(1) > 0.0_wp, 'reflush: the bare cohort grew leaves')
      w_plant1 = st%cohort%leaf_water_mass(1) + st%cohort%wood_water_mass(1)
      call check_close(w_plant1, w_plant0, 1.0e-14_wp, 'reflush: the plant''s water is unchanged')
      psi_leaf = psi_from_water_content(st%cohort%leaf_water_mass(1), leaf_curve, st%cohort%leaf_carbon(1))
      psi_wood = psi_from_water_content(st%cohort%wood_water_mass(1), wood_curve,                   &
                                        st%cohort%sapwood_carbon(1) + st%cohort%fineroot_carbon(1))
      call check(psi_wood < -0.1_wp, 'reflush: the fixture''s stem is drier than the birth potential')
      call check(abs(psi_leaf - (psi_wood - grav_head * st%cohort%height(1))) <= 1.0e-9_wp,        &
                 'reflush: the new leaves start at the wood''s potential less the crown head')
      call check(abs(lg%resid_sum(SLOW_PHASE_GROW)%water) <= 1.0e-12_wp,                           &
                 'reflush: the fill moves water inside the plant; the ledger sees no residual')
   end subroutine check_reflush_water

   !----- The recruit pool is credited EVERY step, and what it holds is the reproduction carbon   !
   !       the parents were debited for. Two assertions, one per half of the contract: the       !
   !       CADENCE (credited every step, not sampled on month boundaries from one day's rate      !
   !       scaled up to stand for the month) and the CARBON LINK (the credit must follow          !
   !       recruit_carbon_efficiency, which is what makes pool x carbon_min the debited carbon). ---!
   !       A third: under recruit_shade_decay the credit falls as exp(-k * LAI) of the patch.      !
   subroutine check_recruit_pool()
      real(wp) :: p1, p2, p_lo, p_hi, p_one, p_shade, lai
      call pool_after_steps(3_ik, 1.0e-3_wp, p1)
      call pool_after_steps(6_ik, 1.0e-3_wp, p2)
      call check(p1 > 0.0_wp, 'recruit pool: credited on an ordinary step, not only at month end')
      !----- Six steps credit about twice what three do. NOT exactly twice: the stand grows       !
      !      between steps, so later steps produce a little more reproduction carbon. The point   !
      !      is that every step contributes -- a monthly sampling would leave p1 and p2 both      !
      !      identically zero, since neither run crosses a month boundary.  --------------------!
      call check(abs(p2 - 2.0_wp * p1) <= 0.01_wp * p2,                                            &
                 'recruit pool: six steps credit ~twice what three do (it integrates, not samples)')
      call pool_after_steps(3_ik, 1.0e-3_wp, p_lo)
      call pool_after_steps(3_ik, 2.0e-3_wp, p_hi)
      call check_close(p_hi, 2.0_wp * p_lo, 1.0e-9_wp,                                             &
                       'recruit pool: the credit follows recruit_carbon_efficiency')
      !----- One step, so both runs see the same stand and its LAI exactly. ----------------------!
      call pool_after_steps(1_ik, 1.0e-3_wp, p_one, lai=lai)
      call pool_after_steps(1_ik, 1.0e-3_wp, p_shade, shade_decay=0.5_wp)
      call check(lai > 0.1_wp, 'recruit pool: the shade test stand has leaves')
      call check_close(p_shade, p_one * exp(-0.5_wp * lai), 1.0e-9_wp,                             &
                       'recruit pool: the credit falls as exp(-recruit_shade_decay * patch LAI)')
   end subroutine check_recruit_pool

   !----- Run `nstep` ordinary (non-month-boundary) slow steps on a fresh stand and return the    !
   !       accumulated recruit pool. Seed rain is off, so what accumulates is reproduction alone. !
   subroutine pool_after_steps(nstep, repro_eff, pool, shade_decay, lai)
      integer(ik), intent(in)  :: nstep
      real(wp),    intent(in)  :: repro_eff
      real(wp),    intent(out) :: pool
      real(wp), optional, intent(in)  :: shade_decay   !< pft.recruit_shade_decay (default 0)
      real(wp), optional, intent(out) :: lai           !< the stand's LAI before the first step
      type(meds_config_t) :: c
      type(site_t)        :: st
      integer(ik)         :: k
      c = build_test_config()
      c%fast_biophysics_on = .true.
      c%demography_on      = .false.
      c%pft%recruit_carbon_efficiency(:) = repro_eff
      if (present(shade_decay)) c%pft%recruit_shade_decay(:) = shade_decay
      c%pft%seed_rain_recruits(:)      = 0.0_wp
      c%pft%leaf_turnover_rate_toc(:)  = 0.01_wp    ! [1/yr] keep turnover from eating the supply
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 40.0_wp)   ! above the reproduction height
      call finalize_init(st)
      if (present(lai)) lai = sum(st%cohort%nplant(1:st%cohort%n) * st%cohort%leaf_area(1:st%cohort%n))
      do k = 1_ik, nstep
         st%cohort%gpp_accum(1:st%cohort%n)       = 1.0_wp
         st%cohort%leaf_resp_accum(1:st%cohort%n) = 0.0_wp
         st%cohort%stem_resp_accum(1:st%cohort%n) = 0.0_wp
         st%cohort%root_resp_accum(1:st%cohort%n) = 0.0_wp
         call vegetation_dynamics(st, c)
      end do
      pool = st%patch%recruit_pool(1, 1)
   end subroutine pool_after_steps


   !---------------------------------------------------------------------------------------!
   ! THE LITTER ACCUMULATOR MUST RIDE THE PATCH ARRAY.                                           !
   !                                                                                          !
   ! restructure_stand (apply_patch_disturbance) CREATES patches after vegetation_dynamics has   !
   ! filled the step's litter, and the consumer (advance_biogeochem_dynamics) loops to the NEW   !
   ! patch count. An accumulator sized to the patch count when it was filled is read past its    !
   ! end, feeding undefined memory into the CENTURY litter input (1.8e9 kgC/m2 on a gap patch    !
   ! born that step, over a 50-year spin-up).                                                    !
   !                                                                                          !
   ! The slow ledger cannot catch that: it DECLARES the same out-of-bounds read as a boundary    !
   ! term and so balances perfectly against the garbage. Only `-check bounds` sees it.           !
   !                                                                                          !
   ! So the assertion is the INVARIANT -- the accumulator is patch state, at least as long as    !
   ! the patch array after the operators that change the patch count have run -- plus the        !
   ! quieter half: a FUSION must carry the litter across area-weighted rather than dropping or   !
   ! duplicating it (a local array misattributes one patch's litter to another).                 !
   !---------------------------------------------------------------------------------------!
   subroutine check_litter_patch_lockstep()
      type(meds_config_t) :: c
      type(site_t)        :: st
      real(wp)    :: tot0, tot1
      integer(ik) :: np0

      c = build_test_config()
      c%fast_biophysics_on   = .true.
      c%demography_on        = .true.        ! the structural operators are the point here
      c%do_patch_disturbance = .true.
      c%soil_carbon_on       = .true.        ! the only consumer of the accumulator
      call init_bare_ground(st, c, 2_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 40.0_wp)
      call add_cohort(st, c, 2_ik, 1_ik, 0.3_wp, 35.0_wp)
      call finalize_init(st)
      st%cohort%gpp_accum(1:st%cohort%n)       = 1.0_wp
      st%cohort%leaf_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%stem_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%root_resp_accum(1:st%cohort%n) = 0.0_wp
      np0 = st%patch%n

      call vegetation_dynamics(st, c)
      call restructure_stand(st, c, .true., .true.)   ! monthly + annual: every operator fires

      !----- THE invariant. `>=` not `==`: the patch block is capacity-allocated, so the array is  !
      !      allowed to be longer than the live patch count -- it must never be SHORTER.  ---------!
      call check(size(st%patch%litter_in) >= st%patch%n,                                           &
                 'the litter accumulator is at least as long as the patch array')
      call check(allocated(st%patch%litter_in), 'the litter accumulator is patch state, not a local')

      !----- A fusion must MOVE the litter, not lose it. Fusing two patches area-weighted leaves   !
      !      the site-wide litter flux unchanged, which is the same rule blend_soil_carbon obeys    !
      !      -- and has to, since this litter is on its way into those pools.  ---------------------!
      if (st%patch%n >= 2_ik) then
         tot0 = litter_site_total(st)
         call fuse_2_patches(st, 1_ik, 2_ik)
         tot1 = litter_site_total(st)
         call check_close(tot1, tot0, 1.0e-12_wp * max(abs(tot0), 1.0_wp),                         &
                          'patch fusion conserves the site-wide litter flux')
      end if
   end subroutine check_litter_patch_lockstep

   !----- With the soil-carbon pools off, the litter is still recorded and still declared leaving at !
   !      the growth phase: the same litter, and the same growth-phase residual, as with them on.     !
   !      Undeclared, the phase would read it as a carbon leak. -----------------------------------!
   subroutine check_litter_without_soil_carbon()
      real(wp) :: lit_on, lit_off, res_on, res_off
      call grow_one_step(.true.,  lit_on,  res_on)
      call grow_one_step(.false., lit_off, res_off)
      call check(lit_off > 0.0_wp, 'soil carbon off: the litter is still recorded')
      call check_close(lit_off, lit_on, 1.0e-12_wp * lit_on, 'soil carbon off: the same litter as with it on')
      call check_close(res_off, res_on, 1.0e-12_wp * max(lit_on, 1.0e-12_wp),                        &
                       'soil carbon off: the growth phase closes as it does with it on')
   end subroutine check_litter_without_soil_carbon

   subroutine grow_one_step(soil_carbon_on, litter, grow_residual)
      logical,  intent(in)  :: soil_carbon_on
      real(wp), intent(out) :: litter, grow_residual
      type(meds_config_t) :: c
      type(site_t)        :: st
      type(slow_ledger_t) :: lg
      c = build_test_config()
      c%fast_biophysics_on = .true.
      c%demography_on      = .false.        ! one step of growth, turnover and mortality, no restructuring
      c%soil_carbon_on     = soil_carbon_on
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 40.0_wp)
      call finalize_init(st)
      st%cohort%gpp_accum(1:st%cohort%n)       = 1.0_wp
      st%cohort%leaf_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%stem_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%root_resp_accum(1:st%cohort%n) = 0.0_wp
      lg%active = .true.
      call slow_ledger_open(lg, st, c, RHO)
      call vegetation_dynamics(st, c, ledger=lg)
      litter        = litter_site_total(st)
      grow_residual = lg%resid_sum(SLOW_PHASE_GROW)%carbon
   end subroutine grow_one_step

   !----- Area-weighted site total of today's litter [kgC/m2/day]. ---------------------------!
   pure function litter_site_total(st) result(t)
      type(site_t), intent(in) :: st
      real(wp)    :: t
      integer(ik) :: ip
      t = 0.0_wp
      do ip = 1_ik, st%patch%n
         t = t + st%patch%area(ip) * (st%patch%litter_in(ip)%labile_grnd                           &
                                    + st%patch%litter_in(ip)%labile_soil                           &
                                    + st%patch%litter_in(ip)%struct_grnd                           &
                                    + st%patch%litter_in(ip)%struct_soil)
      end do
   end function litter_site_total

   subroutine check_mortality_water()
      type(meds_config_t) :: c
      type(site_t)        :: st
      real(wp) :: tis0, tis1, shed0, shed1
      c = build_test_config()
      c%fast_biophysics_on = .true.
      c%demography_on      = .false.        ! isolate the commit from the structural operators
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 40.0_wp)
      call finalize_init(st)
      st%cohort%gpp_accum(1:st%cohort%n)       = 1.0_wp     ! grow, so p1 /= p0 and the basis matters
      st%cohort%leaf_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%stem_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%root_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%leaf_water_mass(1:st%cohort%n) = 0.7_wp
      st%cohort%wood_water_mass(1:st%cohort%n) = 2.3_wp
      tis0  = tissue_water_total(st)
      shed0 = shed_total(st, c)
      call vegetation_dynamics(st, c)
      tis1  = tissue_water_total(st)
      shed1 = shed_total(st, c)
      call check(tis1 < tis0, 'mortality water: the tissue store actually falls')
      call check_close(shed1 - shed0, tis0 - tis1, 1.0e-12_wp,                                     &
                       'every kg that leaves tissue must arrive in the patch shed channel')
   end subroutine check_mortality_water

   subroutine check_shed_heat()
      type(pool_change_block) :: dpool
      type(slow_ledger_t)     :: lg
      integer(ik) :: nc
      nc = site%cohort%n
      allocate(dpool%leaf(nc), dpool%fineroot(nc), dpool%wood(nc), dpool%nonstructural(nc))
      dpool%wood = 0.0_wp ; dpool%nonstructural = 0.0_wp
      dpool%leaf     = -0.25_wp * site%cohort%leaf_carbon(1:nc)       ! a quarter of the leaves shed
      dpool%fineroot = -0.10_wp * site%cohort%fineroot_carbon(1:nc)   ! a tenth of the fine roots
      site%cohort%leaf_water_mass(1:nc) = 0.7_wp
      site%cohort%wood_water_mass(1:nc) = 2.3_wp
      lg%active = .true.
      call slow_ledger_open(lg, site, cfg, RHO)
      call shed_turnover_water(site, cfg, dpool, lg)
      call slow_ledger_mark(lg, site, cfg, SLOW_PHASE_ALLOCATE)
      call check(lg%decl_gross(SLOW_PHASE_ALLOCATE)%energy > 0.0_wp, 'shed heat: the shed declares energy')
      call check(abs(lg%resid_sum(SLOW_PHASE_ALLOCATE)%energy)                                      &
                 <= 1.0e-12_wp * lg%decl_gross(SLOW_PHASE_ALLOCATE)%energy,                            &
                 'shed heat: the declared energy is the tissue heat store''s change')
      call check(abs(lg%resid_sum(SLOW_PHASE_ALLOCATE)%water)                                       &
                 <= 1.0e-12_wp * lg%decl_gross(SLOW_PHASE_ALLOCATE)%water,                             &
                 'shed heat: the declared water is the tissue water store''s change')
   end subroutine check_shed_heat

   pure function tissue_water_total(st) result(w)
      type(site_t), intent(in) :: st
      real(wp)    :: w
      integer(ik) :: ip, i, i0, i1
      w = 0.0_wp
      do ip = 1_ik, st%patch%n
         i0 = st%patch%cohort_offset(ip) ; i1 = i0 + st%patch%cohort_count(ip) - 1_ik
         do i = i0, i1
            w = w + st%patch%area(ip) * st%cohort%nplant(i)                                        &
                  * (st%cohort%leaf_water_mass(i) + st%cohort%wood_water_mass(i))
         end do
      end do
   end function tissue_water_total

   pure function shed_total(st, c) result(w)
      type(site_t),        intent(in) :: st
      type(meds_config_t), intent(in) :: c
      real(wp)    :: w
      integer(ik) :: ip
      w = 0.0_wp
      do ip = 1_ik, st%patch%n
         w = w + st%patch%area(ip) * st%patch%shed_water_rate(ip) * c%dt_slow
      end do
   end function shed_total


   !----- A cohort's heat capacity is a function of its biomass, so changing biomass alone moves   !
   !       the tissue store with no flux. That is the exchange the growth commit declares; these   !
   !       assert the two properties the declaration rests on -- that the term is real (biomass    !
   !       moves it) and that slow_tissue_heat is exactly the store's tissue slice (so declaring   !
   !       the one cancels the other).  ----------------------------------------------------------!
   subroutine check_tissue_thermal_mass()
      type(slow_store_t) :: a, b
      real(wp) :: e0, e1, de_store, de_tissue
      e0 = slow_tissue_heat(site)
      a  = slow_site_store(site, cfg, soil, RHO)
      !----- Grow one cohort's wood: temperatures untouched, so ONLY the thermal mass changes. ----!
      site%cohort%wood_carbon(1) = site%cohort%wood_carbon(1) * 1.5_wp
      e1 = slow_tissue_heat(site)
      b  = slow_site_store(site, cfg, soil, RHO)
      de_tissue = e1 - e0
      de_store  = b%energy - a%energy
      call check(de_tissue > 0.0_wp, 'thermal mass: growing biomass raises the tissue heat store')
      call check_close(de_store, de_tissue, 1.0e-9_wp,                                             &
                       'thermal mass: slow_tissue_heat is exactly the store''s tissue slice')
      !----- and it is the BIOMASS doing it, not the temperature: cool the cohort and the store    !
      !      must fall, so a declaration written against biomass alone would not cancel that.  ----!
      site%cohort%leaf_temp(1) = site%cohort%leaf_temp(1) - 10.0_wp
      call check(slow_tissue_heat(site) < e1, 'thermal mass: temperature moves it too (so T is not frozen in)')
   end subroutine check_tissue_thermal_mass

   subroutine check_routing()
      real(wp) :: co2_on, co2_off, lit_lossy, lit_perfect, dummy
      real(wp) :: lab0, lab1, wood0, wood1, n0, dd, dty
      !----- Each variant runs on a FRESH stand. vegetation_dynamics COMMITS growth, so calling it !
      !      twice on one site compares two different forests and the second comparison passes for !
      !      the wrong reason (a mutation that should fail it does not).  -------------------------!
      call run_variant(0.3_wp, 1.0_wp, co2_on,  dummy)
      call run_variant(0.0_wp, 1.0_wp, co2_off, dummy)
      call check(co2_on > 0.0_wp, 'routing: growth respiration reaches patch%slow_co2_rate')
      call check(abs(co2_off) < abs(co2_on),                                                       &
                 'routing: zero construction cost => nothing owed (the channel is not a constant)')

      call run_variant(0.3_wp, 0.0_wp, dummy, lit_lossy)     ! every seed dies -> all necromass
      call run_variant(0.3_wp, 1.0_wp, dummy, lit_perfect)   ! every seed establishes -> none
      call check(lit_lossy > lit_perfect,                                                          &
                 'routing: the unestablished seed fraction becomes litter, not nothing')

      !----- STARVATION rides the same channel with the OPPOSITE sign. The fast loop has already   !
      !      exhaled the full maintenance respiration; a stand that could not fund it means the    !
      !      loop OVER-reported, so the correction owed to the canopy air is NEGATIVE. A sign      !
      !      error here would quietly turn a starving forest into a CO2 source.  -----------------!
      call run_variant(0.3_wp, 1.0_wp, co2_on, dummy, starve=.true.)
      call check(co2_on < 0.0_wp, 'routing: a starving stand OWES the canopy air a negative flux')

      !----- SINK LIMIT: the wood a binding limit leaves unbuilt reaches the soil as exudate, one for !
      !      one, and since exudate costs (1+g) like wood the CO2 owed does not change. Mortality is   !
      !      off in both runs: it rises as growth falls, so the limited stand would also lose more     !
      !      trees to litter. -------------------------------------------------------------------!
      call run_variant(0.3_wp, 1.0_wp, co2_on, dummy, lab_soil=lab0, wood=wood0, nplant0=n0,       &
                       immortal=.true.)
      call run_variant(0.3_wp, 1.0_wp, co2_off, dummy, max_rgr=1.0e-3_wp, lab_soil=lab1, wood=wood1, &
                       immortal=.true.)
      call check(wood1 < wood0, 'sink limit: a binding limit builds less wood')
      call check_close(lab1 - lab0, (wood0 - wood1) * n0, 1.0e-9_wp,                               &
                       'sink limit: the unbuilt wood arrives as below-ground labile litter')
      call check_close(co2_off, co2_on, 1.0e-9_wp,                                                 &
                       'sink limit: the exudate pays the growth respiration the wood would have')
      !----- An ABSOLUTE limit that binds grows the stem by exactly the limit: 0.2 cm/yr for a step. -!
      call run_variant(0.3_wp, 1.0_wp, dummy, dummy, max_agr=0.2_wp, ddbh=dd, dt_yr=dty, immortal=.true.)
      call check_close(dd, 0.2_wp * dty, 1.0e-9_wp, 'sink limit: a binding absolute limit grows dbh by rate x dt')
   end subroutine check_routing

   !----- One slow step on a FRESH stand, returning the two routed quantities (CO2 rate, litter). --!
   subroutine run_variant(g_resp, repro_eff, co2_rate, litter_total, starve, max_rgr,              &
                          lab_soil, wood, nplant0, immortal, max_agr, ddbh, dt_yr)
      real(wp), intent(in)  :: g_resp, repro_eff
      real(wp), intent(out) :: co2_rate, litter_total
      logical, optional, intent(in) :: starve
      real(wp), optional, intent(in)  :: max_rgr    !< pft.max_relative_growth_rate (default 0, no limit)
      real(wp), optional, intent(out) :: lab_soil   !< patch 1's below-ground labile litter
      real(wp), optional, intent(out) :: wood       !< the cohort's wood carbon after the step
      real(wp), optional, intent(out) :: nplant0    !< the cohort's density before the step
      logical,  optional, intent(in)  :: immortal   !< .true. => no background or low-growth mortality
      real(wp), optional, intent(in)  :: max_agr    !< pft.max_absolute_growth_rate (default 0, no limit)
      real(wp), optional, intent(out) :: ddbh       !< the cohort's dbh growth over the step
      real(wp), optional, intent(out) :: dt_yr      !< the step's length [yr]
      real(wp) :: dbh_start
      type(meds_config_t) :: c
      type(site_t)        :: st
      c = build_test_config()
      c%fast_biophysics_on = .true.
      c%demography_on      = .false.        ! isolate allocation from the structural operators
      c%soil_carbon_on     = .true.
      c%pft%growth_resp_factor(:)      = g_resp
      c%pft%recruit_carbon_efficiency(:) = repro_eff
      !----- A SLOW leaf turnover, so turnover does not eat the whole supply before reproduction !
      !      is reached. The allocator's priority order is leaf/root growth, then storage, then    !
      !      reproduction, then wood -- with the default ~1 /yr turnover this cohort's daily leaf  !
      !      replacement is 14x the photosynthate below and `avail` reaches zero before any seed   !
      !      is made, so the seed-loss assertion would pass while asserting nothing.               !
      c%pft%leaf_turnover_rate_toc(:) = 0.01_wp
      if (present(max_rgr)) c%pft%max_relative_growth_rate(:) = max_rgr
      if (present(max_agr)) c%pft%max_absolute_growth_rate(:) = max_agr
      if (present(dt_yr))   dt_yr = c%dt_years
      if (present(immortal)) then
         if (immortal) then
            c%pft%mort_gamma(:) = 0.0_wp
            c%pft%mort_alpha(:) = 0.0_wp
         end if
      end if
      call init_bare_ground(st, c, 1_ik)
      !----- 40 cm => ~25 m, clear of the reproduction height (20 m): a cohort below it allocates  !
      !      NOTHING to reproduction, and the seed-loss assertion would then be vacuously true.  --!
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 40.0_wp)
      call finalize_init(st)
      st%cohort%gpp_accum(1:st%cohort%n)       = 1.0_wp      ! generous: a clear surplus to allocate
      st%cohort%leaf_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%stem_resp_accum(1:st%cohort%n) = 0.0_wp
      st%cohort%root_resp_accum(1:st%cohort%n) = 0.0_wp
      if (present(starve)) then
         if (starve) then
            !----- No photosynthate, a real maintenance bill, and no reserves to pay it with. ----!
            st%cohort%gpp_accum(1:st%cohort%n)            = 0.0_wp
            st%cohort%leaf_resp_accum(1:st%cohort%n)      = 1.0_wp
            st%cohort%nonstructural_carbon(1:st%cohort%n) = 0.0_wp
         end if
      end if
      if (present(nplant0)) nplant0 = st%cohort%nplant(1)
      dbh_start = st%cohort%dbh(1)
      call vegetation_dynamics(st, c)
      if (present(ddbh)) ddbh = st%cohort%dbh(1) - dbh_start
      if (present(lab_soil)) lab_soil = st%patch%litter_in(1)%labile_soil
      if (present(wood))     wood     = st%cohort%wood_carbon(1)
      co2_rate     = st%patch%slow_co2_rate(1)
      !----- The litter accumulator is patch state, so it is read off the site rather than
      !      returned: st%patch%litter_in(1), in lockstep with patch 1.
      litter_total = st%patch%litter_in(1)%labile_grnd + st%patch%litter_in(1)%labile_soil          &
                   + st%patch%litter_in(1)%struct_grnd + st%patch%litter_in(1)%struct_soil
   end subroutine run_variant

   subroutine check_cull()
      real(wp) :: w_before, shed_before, w_shed, expect_w, heat_before, heat_after
      integer(ik) :: n0
      !----- Drive cohort 2 below the density floor so the cull takes it, and nothing else. ------!
      site%cohort%nplant(2)   = 0.5_wp * cfg%negligible_nplant
      site%cohort%leaf_temp(2) = 291.0_wp ; site%cohort%wood_temp(2) = 289.0_wp
      n0          = site%cohort%n
      shed_before = site%patch%shed_water_rate(1)
      heat_before = slow_tissue_heat(site)
      w_before    = site%cohort%nplant(2) * (site%cohort%leaf_water_mass(2) + site%cohort%wood_water_mass(2)) &
                  + site%cohort%leaf_surf_water(2) + site%cohort%wood_surf_water(2)
      expect_w    = site%patch%area(1) * w_before

      call terminate_cohorts(site, cfg, w_shed)

      call check(site%cohort%n == n0 - 1_ik, 'cull: the sub-floor cohort is removed')
      call check_close(w_shed, expect_w, 1.0e-12_wp, 'cull: reports the tissue + film water it carried')
      call check_close(site%patch%shed_water_rate(1) - shed_before, w_before / cfg%dt_slow,       &
                       1.0e-12_wp, 'cull: the water reaches the patch shed channel, not the void')
      !----- The cull's tissue HEAT is not reported by the operator: it is one part of the        !
      !      phase's thermal-mass change, which the driver brackets with slow_tissue_heat. Assert !
      !      that the store SEES the loss, which is what makes that bracket able to declare it.   !
      heat_after = slow_tissue_heat(site)
      call check(heat_after < heat_before, 'cull: the tissue heat store falls when a cohort goes')
   end subroutine check_cull

   !----- A cohort is born in its patch's state (patch_physical_state): at the canopy-air temperature, !
   !      its tissues at the top soil layer's potential. A patch not seeded yet gives the run's       !
   !      starting state from [init]. ------------------------------------------------------------------!
   subroutine check_birth_temp()
      integer(ik) :: m
      real(wp)    :: t_birth, psi_birth, t_save, theta_start(n_soil_layer_max)
      type(water_curve_t) :: leaf_curve, wood_curve
      m = site%cohort%n + 1_ik
      t_save = site%patch%cas(1)%can_temp
      site%patch%cas(1)%can_temp = 268.0_wp        ! a cold, seeded patch
      call patch_physical_state(site%patch%cas(1), site%patch%soil_water(1), cfg, soil, t_birth, psi_birth)
      call check_close(t_birth, 268.0_wp, 1.0e-12_wp, 'birth: a seeded patch gives its canopy-air temperature')
      call check_close(psi_birth, top_soil_psi(soil, site%patch%soil_water(1)%theta(1)), 1.0e-12_wp,   &
                       'birth: a seeded patch gives its top soil layer''s potential')
      call init_cohort(site%cohort, m, cfg%pft, 1_ik, 1_ik, 0.01_wp, 0.45_wp, cfg%hydraulics, t_birth, psi_birth)
      call check_close(site%cohort%leaf_temp(m), 268.0_wp, 1.0e-12_wp, 'birth: leaf starts at the patch temperature')
      call check_close(site%cohort%wood_temp(m), 268.0_wp, 1.0e-12_wp, 'birth: wood starts at the patch temperature')
      call pft_tissue_curves(cfg%hydraulics, cfg%pft, 1_ik, leaf_curve, wood_curve)
      call check_close(psi_from_water_content(site%cohort%leaf_water_mass(m), leaf_curve,                &
                       site%cohort%leaf_carbon(m)), psi_birth, 1.0e-9_wp,                               &
                       'birth: the leaves hold the water of that potential')
      !----- A patch not seeded yet (its canopy air at the type's UNSET_POSITIVE) gives the starting state. !
      site%patch%cas(1)%can_temp = UNSET_POSITIVE
      call patch_physical_state(site%patch%cas(1), site%patch%soil_water(1), cfg, soil, t_birth, psi_birth)
      call soil_column_start_theta(cfg%soil_column, soil, cfg%init_soil_theta, theta_start)
      call check_close(t_birth, cfg%init_soil_temp, 1.0e-12_wp, 'birth: unseeded, the run''s starting temperature')
      call check_close(psi_birth, top_soil_psi(soil, theta_start(1)), 1.0e-12_wp,                       &
                       'birth: unseeded, the starting soil water''s potential')
      site%patch%cas(1)%can_temp = t_save
   end subroutine check_birth_temp

   !----- Assert a perturbation moves the named currency by `amount` (and moves it at all). ----!
   subroutine moves_water(name, perturb, amount)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: amount
      interface
         subroutine perturb()
         end subroutine perturb
      end interface
      type(slow_store_t) :: a, b
      a = slow_site_store(site, cfg, soil, RHO)
      call perturb()
      b = slow_site_store(site, cfg, soil, RHO)
      call check_close(b%water - a%water, amount, 1.0e-9_wp * max(1.0_wp, abs(amount)),            &
                       'store: '//name//' reaches the water total')
   end subroutine moves_water

   subroutine moves_energy(name, perturb, amount)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: amount
      interface
         subroutine perturb()
         end subroutine perturb
      end interface
      type(slow_store_t) :: a, b
      a = slow_site_store(site, cfg, soil, RHO)
      call perturb()
      b = slow_site_store(site, cfg, soil, RHO)
      call check_close(b%energy - a%energy, amount, 1.0e-6_wp * max(1.0_wp, abs(amount)),          &
                       'store: '//name//' reaches the energy total')
   end subroutine moves_energy

   subroutine moves_carbon(name, perturb, amount)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: amount
      interface
         subroutine perturb()
         end subroutine perturb
      end interface
      type(slow_store_t) :: a, b
      a = slow_site_store(site, cfg, soil, RHO)
      call perturb()
      b = slow_site_store(site, cfg, soil, RHO)
      call check_close(b%carbon - a%carbon, amount, 1.0e-9_wp * max(1.0_wp, abs(amount)),          &
                       'store: '//name//' reaches the carbon total')
   end subroutine moves_carbon

   subroutine moves(name, perturb, currency)
      character(len=*), intent(in) :: name, currency
      interface
         subroutine perturb()
         end subroutine perturb
      end interface
      type(slow_store_t) :: a, b
      real(wp) :: d
      a = slow_site_store(site, cfg, soil, RHO)
      call perturb()
      b = slow_site_store(site, cfg, soil, RHO)
      select case (currency)
      case ('carbon') ; d = b%carbon - a%carbon
      case ('water')  ; d = b%water  - a%water
      case default    ; d = b%energy - a%energy
      end select
      call check(abs(d) > 0.0_wp, 'store: '//name//' reaches the '//currency//' total')
   end subroutine moves

   subroutine perturb_theta()        ; site%patch%soil_water(1)%theta(1) = site%patch%soil_water(1)%theta(1) + 0.01_wp
   end subroutine perturb_theta
   subroutine perturb_pond()         ; site%patch%soil_water(1)%w_surface = site%patch%soil_water(1)%w_surface + 1.0_wp
   end subroutine perturb_pond
   subroutine perturb_swe()          ; site%patch%snow(1)%swe(1) = site%patch%snow(1)%swe(1) + 1.0_wp
   end subroutine perturb_swe
   subroutine perturb_shv()          ; site%patch%cas(1)%can_shv = site%patch%cas(1)%can_shv + 0.001_wp
   end subroutine perturb_shv
   subroutine perturb_tissue_water() ; site%cohort%leaf_water_mass(1) = site%cohort%leaf_water_mass(1) + 0.01_wp
   end subroutine perturb_tissue_water
   subroutine perturb_film()         ; site%cohort%leaf_surf_water(1) = site%cohort%leaf_surf_water(1) + 0.01_wp
   end subroutine perturb_film
   subroutine perturb_soil_energy()
      site%patch%soil_energy(1)%soil_energy(1) = site%patch%soil_energy(1)%soil_energy(1) + 1.0e6_wp
   end subroutine perturb_soil_energy
   subroutine perturb_pond_energy()
      site%patch%soil_water(1)%w_surface_energy = site%patch%soil_water(1)%w_surface_energy + 1.0e5_wp
   end subroutine perturb_pond_energy
   subroutine perturb_snow_energy()
      site%patch%snow(1)%snow_energy(1) = site%patch%snow(1)%snow_energy(1) + 1.0e5_wp
   end subroutine perturb_snow_energy
   subroutine perturb_cas_enth()     ; site%patch%cas(1)%can_enthalpy = site%patch%cas(1)%can_enthalpy + 1.0e4_wp
   end subroutine perturb_cas_enth
   subroutine perturb_cas_co2()      ; site%patch%cas(1)%can_co2 = site%patch%cas(1)%can_co2 + 50.0_wp
   end subroutine perturb_cas_co2
   subroutine perturb_live_c()       ; site%cohort%wood_carbon(1) = site%cohort%wood_carbon(1) + 0.5_wp
   end subroutine perturb_live_c
   subroutine perturb_necromass()    ; site%patch%soil_carbon(1)%slow_carbon = site%patch%soil_carbon(1)%slow_carbon + 1.0_wp
   end subroutine perturb_necromass
   subroutine perturb_tissue_temp()  ; site%cohort%leaf_temp(1) = site%cohort%leaf_temp(1) + 5.0_wp
   end subroutine perturb_tissue_temp
   subroutine perturb_recruit()      ; site%patch%recruit_pool(1,1) = site%patch%recruit_pool(1,1) + 0.5_wp
   end subroutine perturb_recruit

end program test_slow_ledger
