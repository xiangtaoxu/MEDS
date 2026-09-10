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
! than passing a quantity through, and the run it was written against closes it to 1e-17.         !
!==========================================================================================!
program test_slow_ledger
   use meds_kinds,              only : wp, ik
   use meds_config,             only : meds_config_t
   use meds_site_state_types,   only : site_t
   use meds_column_state_types, only : cas_set_depth
   use meds_column_params,      only : soil_params_t, build_soil_hydr_params
   use meds_init,               only : init_bare_ground, add_cohort, finalize_init
   use meds_site_state_types,   only : init_cohort, cohort_tissue_heat_capacity,                  &
                                       TISSUE_C_LEAF, TISSUE_C_SAPW, TISSUE_HCAP_MIN
   use meds_demography_cohort_fusefiss, only : terminate_cohorts
   use meds_demography_patch_fusefiss,  only : fuse_2_patches
   use meds_vegetation_dynamics,        only : vegetation_dynamics
   use meds_biogeochem_types,           only : litter_input_t
   use meds_slow_ledger,        only : slow_store_t, slow_ledger_t, slow_site_store,               &
                                       slow_ledger_open, slow_ledger_declare, slow_ledger_mark,    &
                                       slow_fast_carbon_handover, slow_tissue_heat, SLOW_PHASE_GROW
   use meds_test_support,       only : build_test_config, check, check_close, banner
   implicit none

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
      call build_soil_hydr_params(sc%n_layer, sc%retention, sc%depth, sc%grid_growth, sc%theta_sat, &
                                  sc%theta_res, sc%ksat, sc%curve_par_a, sc%curve_par_n,            &
                                  sc%root_beta, sc%psi_fc, soil)
   end associate

   !----- Give every reservoir a non-zero value, so a perturbation below is a change to a store !
   !      that already exists rather than the creation of one. --------------------------------!
   site%patch%soil_w(1)%theta(1:soil%n_active)      = 0.30_wp
   site%patch%soil_e(1)%soil_energy(1:soil%n_active) = 5.0e7_wp
   site%patch%soil_w(1)%w_surface      = 2.0_wp
   site%patch%soil_w(1)%w_surface_enth = 8.0e5_wp
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
   call moves_energy('soil energy',       perturb_soil_e,       1.0e6_wp * soil%dz(1))
   call moves_energy('pond enthalpy',     perturb_pond_enth,    1.0e5_wp)
   call moves_energy('snow enthalpy',     perturb_snow_enth,    1.0e5_wp)
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

   !=== 5. A CULL hands its water to the ground channel and reports its heat. =================!
   !       Before this, cohort_compact simply dropped both (review item 1B #7). The assertion is  !
   !       the TRANSFER, not the disappearance: what leaves the cohort must arrive in the patch's !
   !       shed channel, and the reported heat must be the capacity the cohort actually had.      !
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
   !       This identity broke twice while it was being written -- once by not routing the water  !
   !       at all, once by scaling the loss with the SURVIVORS' density instead of the density     !
   !       DROP -- and in both cases every other assertion in this file still passed. It is the    !
   !       one property that pins the basis.                                                        !
   call check_mortality_water()
   call check_litter_patch_lockstep()

   !=== 10. Reproduction carbon reaches the recruit pool, every step and in proportion. =======!
   call check_recruit_pool()

   print '(a)', 'test_slow_ledger: ALL PASSED'

contains

   !----- The recruit pool is credited EVERY step, and what it holds is the reproduction carbon   !
   !       the parents were debited for. Two assertions, because the fix has two halves: the      !
   !       CADENCE (it used to be credited only on month boundaries, from that one day's rate     !
   !       scaled up to stand for the month) and the CARBON LINK (the credit must follow          !
   !       repro_carbon_efficiency, which is what makes pool x carbon_min the debited carbon). ---!
   subroutine check_recruit_pool()
      real(wp) :: p1, p2, p_lo, p_hi
      call pool_after_steps(3_ik, 1.0e-3_wp, p1)
      call pool_after_steps(6_ik, 1.0e-3_wp, p2)
      call check(p1 > 0.0_wp, 'recruit pool: credited on an ordinary step, not only at month end')
      !----- Six steps credit about twice what three do. NOT exactly twice: the stand grows       !
      !      between steps, so later steps produce a little more reproduction carbon. The point   !
      !      is that every step contributes -- under the old monthly sampling p1 and p2 were both !
      !      identically zero, since neither run crosses a month boundary.  --------------------!
      call check(abs(p2 - 2.0_wp * p1) <= 0.01_wp * p2,                                            &
                 'recruit pool: six steps credit ~twice what three do (it integrates, not samples)')
      call pool_after_steps(3_ik, 1.0e-3_wp, p_lo)
      call pool_after_steps(3_ik, 2.0e-3_wp, p_hi)
      call check_close(p_hi, 2.0_wp * p_lo, 1.0e-9_wp,                                             &
                       'recruit pool: the credit follows repro_carbon_efficiency')
   end subroutine check_recruit_pool

   !----- Run `nstep` ordinary (non-month-boundary) slow steps on a fresh stand and return the    !
   !       accumulated recruit pool. Seed rain is off, so what accumulates is reproduction alone. !
   subroutine pool_after_steps(nstep, repro_eff, pool)
      integer(ik), intent(in)  :: nstep
      real(wp),    intent(in)  :: repro_eff
      real(wp),    intent(out) :: pool
      type(meds_config_t) :: c
      type(site_t)        :: st
      integer(ik)         :: k
      c = build_test_config()
      c%fast_biophysics_on = .true.
      c%demography_on      = .false.
      c%pft%repro_carbon_efficiency(:) = repro_eff
      c%pft%seed_rain_recruits(:)      = 0.0_wp
      c%pft%leaf_lifespan_toc(:)       = 100.0_wp   ! keep turnover from eating the supply
      call init_bare_ground(st, c, 1_ik)
      call add_cohort(st, c, 1_ik, 1_ik, 0.3_wp, 40.0_wp)   ! above min_reproduction_height
      call finalize_init(st)
      do k = 1_ik, nstep
         st%cohort%gpp_accum(1:st%cohort%n)       = 1.0_wp
         st%cohort%leaf_resp_accum(1:st%cohort%n) = 0.0_wp
         st%cohort%stem_resp_accum(1:st%cohort%n) = 0.0_wp
         st%cohort%root_resp_accum(1:st%cohort%n) = 0.0_wp
         call vegetation_dynamics(st, c, .false., .false.)
      end do
      pool = st%patch%recruit_pool(1, 1)
   end subroutine pool_after_steps


   !---------------------------------------------------------------------------------------!
   ! THE LITTER ACCUMULATOR MUST RIDE THE PATCH ARRAY.                                           !
   !                                                                                          !
   ! `lit` used to be a LOCAL array in vegetation_dynamics, sized to the patch count on entry --  !
   ! and apply_patch_disturbance CREATES a patch partway through that same routine. The consumer  !
   ! (advance_biogeochem_dynamics) then looped to the NEW patch count and read past the end,       !
   ! feeding undefined memory into the CENTURY litter input. A 50-year spin-up reached 1.8e9        !
   ! kgC/m2 on a gap patch born that step.                                                          !
   !                                                                                          !
   ! Nothing caught it. The whole suite was green; the slow ledger closed on every phase and         !
   ! currency, because the ledger DECLARED the same out-of-bounds read as a boundary term and so     !
   ! balanced perfectly against the garbage. Only `-check bounds` saw it, and only a config with     !
   ! soil carbon on ever reached the code.                                                            !
   !                                                                                          !
   ! So the assertion is the INVARIANT that was violated -- the accumulator is at least as long as   !
   ! the patch array, after the operators that change the patch count have run -- plus the quieter   !
   ! half: a FUSION must carry the litter across area-weighted rather than dropping or duplicating   !
   ! it (the old local array silently misattributed one patch's litter to another).                   !
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

      call vegetation_dynamics(st, c, .true., .true.)   ! monthly + annual: every operator fires

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
      call vegetation_dynamics(st, c, .false., .false.)
      tis1  = tissue_water_total(st)
      shed1 = shed_total(st, c)
      call check(tis1 < tis0, 'mortality water: the tissue store actually falls')
      call check_close(shed1 - shed0, tis0 - tis1, 1.0e-12_wp,                                     &
                       'every kg that leaves tissue must arrive in the patch shed channel')
   end subroutine check_mortality_water

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
      !----- Each variant runs on a FRESH stand. vegetation_dynamics COMMITS growth, so calling it !
      !      twice on one site compares two different forests and the second comparison passes for !
      !      the wrong reason -- which is exactly what happened before this was rebuilt, and the    !
      !      mutation that should have failed it did not.  ---------------------------------------!
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
   end subroutine check_routing

   !----- One slow step on a FRESH stand, returning the two quantities this PR routes. ------------!
   subroutine run_variant(g_resp, repro_eff, co2_rate, litter_total, starve)
      real(wp), intent(in)  :: g_resp, repro_eff
      real(wp), intent(out) :: co2_rate, litter_total
      logical, optional, intent(in) :: starve
      type(meds_config_t) :: c
      type(site_t)        :: st
      c = build_test_config()
      c%fast_biophysics_on = .true.
      c%demography_on      = .false.        ! isolate allocation from the structural operators
      c%soil_carbon_on     = .true.
      c%pft%growth_resp_factor(:)      = g_resp
      c%pft%repro_carbon_efficiency(:) = repro_eff
      !----- A LONG leaf lifespan, so turnover does not eat the whole supply before reproduction  !
      !      is reached. The allocator's priority order is leaf/root growth, then storage, then    !
      !      reproduction, then wood -- with the default ~1 yr lifespan this cohort's daily leaf   !
      !      replacement is 14x the photosynthate below and `avail` reaches zero before any seed   !
      !      is made, which is how the first version of this test passed while asserting nothing.  !
      c%pft%leaf_lifespan_toc(:) = 100.0_wp
      call init_bare_ground(st, c, 1_ik)
      !----- 40 cm => ~25 m, clear of min_reproduction_height (20 m): a cohort below it allocates  !
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
      call vegetation_dynamics(st, c, .false., .false.)
      co2_rate     = st%patch%slow_co2_rate(1)
      !----- The litter accumulator is patch state now, so it is read off the site rather than
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
      !----- The cull's tissue HEAT is no longer reported by the operator: it is one part of the  !
      !      phase's thermal-mass change, which the driver brackets with slow_tissue_heat. Assert !
      !      that the store SEES the loss, which is what makes that bracket able to declare it.   !
      heat_after = slow_tissue_heat(site)
      call check(heat_after < heat_before, 'cull: the tissue heat store falls when a cohort goes')
   end subroutine check_cull

   subroutine check_birth_temp()
      integer(ik) :: m
      m = site%cohort%n + 1_ik
      site%patch%cas(1)%can_temp = 268.0_wp        ! a cold patch: LEAF_TEMP_INIT would be 20 K out
      call init_cohort(site%cohort, m, cfg%pft, 1_ik, 1_ik, 0.01_wp, 0.45_wp, birth_temp=268.0_wp)
      call check_close(site%cohort%leaf_temp(m), 268.0_wp, 1.0e-12_wp, 'birth: leaf starts at the patch temperature')
      call check_close(site%cohort%wood_temp(m), 268.0_wp, 1.0e-12_wp, 'birth: wood starts at the patch temperature')
      !----- and the default is still the constant, for a setup call with no environment yet. ----!
      call init_cohort(site%cohort, m, cfg%pft, 1_ik, 1_ik, 0.01_wp, 0.45_wp)
      call check(site%cohort%leaf_temp(m) > 280.0_wp, 'birth: absent a patch temperature, the constant stands')
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

   subroutine perturb_theta()        ; site%patch%soil_w(1)%theta(1) = site%patch%soil_w(1)%theta(1) + 0.01_wp
   end subroutine perturb_theta
   subroutine perturb_pond()         ; site%patch%soil_w(1)%w_surface = site%patch%soil_w(1)%w_surface + 1.0_wp
   end subroutine perturb_pond
   subroutine perturb_swe()          ; site%patch%snow(1)%swe(1) = site%patch%snow(1)%swe(1) + 1.0_wp
   end subroutine perturb_swe
   subroutine perturb_shv()          ; site%patch%cas(1)%can_shv = site%patch%cas(1)%can_shv + 0.001_wp
   end subroutine perturb_shv
   subroutine perturb_tissue_water() ; site%cohort%leaf_water_mass(1) = site%cohort%leaf_water_mass(1) + 0.01_wp
   end subroutine perturb_tissue_water
   subroutine perturb_film()         ; site%cohort%leaf_surf_water(1) = site%cohort%leaf_surf_water(1) + 0.01_wp
   end subroutine perturb_film
   subroutine perturb_soil_e()       ; site%patch%soil_e(1)%soil_energy(1) = site%patch%soil_e(1)%soil_energy(1) + 1.0e6_wp
   end subroutine perturb_soil_e
   subroutine perturb_pond_enth()    ; site%patch%soil_w(1)%w_surface_enth = site%patch%soil_w(1)%w_surface_enth + 1.0e5_wp
   end subroutine perturb_pond_enth
   subroutine perturb_snow_enth()    ; site%patch%snow(1)%snow_energy(1) = site%patch%snow(1)%snow_energy(1) + 1.0e5_wp
   end subroutine perturb_snow_enth
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
