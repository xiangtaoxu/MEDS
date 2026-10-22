! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_spac_solver -- the soil-plant water solver, meds_spac_solver                          !
! (docs/science/soil_plant_water.md).                                                      !
!                                                                                          !
!   1. ARROWHEAD SOLVE : the Schur solve equals dense elimination on a random arrowhead system. !
!   2. JACOBIAN        : the assembled Jacobian equals the residual's derivative (wet, dry,      !
!                        ponded, leafless).                                                       !
!   3. CONSERVATION    : two days of rain, transpiration, evaporation and drainage close to       !
!                        round-off, step by step.                                                 !
!   4. PLANT ONLY      : soil held at its potential (prescribed core nodes): the steady state     !
!                        carries the transpiration through every path.                              !
!   5. SOIL ONLY       : no cohorts, a sealed bottom: the column relaxes to hydrostatic equilibrium. !
!   6. STIFF EXCHANGE  : a wet, densely rooted stand whose wood-soil exchange takes a fraction of   !
!                        a second costs a handful of sub-steps, not thousands.                     !
!   7. INFEASIBLE DEMAND: a leaf asked for many times the water it and its xylem can supply stays  !
!                        above empty and transpires less than asked.                               !
!   8. LEAFLESS        : a cohort with no leaves: its leaf sits at the wood less the crown head.    !
!   9. DRY-LAYER STORM : a thin top layer near residual, a stand refilling its wood, then a       !
!                        storm: converges, no store leaves its range.                              !
!  10. SATURATED       : rain on a saturated column ponds, spills, and conserves.                   !
!  11. ACCURACY        : the default tolerances track a tight-tolerance run.                       !
!  12. VANISHING STORES: a near-leafless leaf and a trace pond stay in range, and water closes.   !
!  13. CAMPBELL        : on a Campbell clay the soil potentials are the Campbell curve's, finite  !
!                        and below zero; the stand takes up water, and every step closes.         !
!==========================================================================================!
program test_spac_solver
   use meds_test_assert, only : check_true, check_close, check_int, test_report
   use, intrinsic :: iso_fortran_env, only : int64
   use meds_kinds,       only : wp, ik
   use meds_constants,   only : grav_head, rho_h2o, pi, tiny_num, n_soil_layer_max
   use meds_linear_algebra,  only : arrowhead_system_t, arrowhead_allocate, arrowhead_factor, arrowhead_solve
   use meds_spac_solver, only : spac_solver_t, spac_allocate, spac_set_state, spac_advance,      &
                                spac_jacobian_error, POND, WOOD, LEAF
   use meds_plant_hydr_lib, only : xylem_max_conductance, rhizosphere_conductance, build_hydro_table,  &
                                   tissue_water_relations, water_curve_t, tissue_psi_to_rwc, tissue_psi_tlp
   use meds_soil_hydr_lib,   only : soil_water_relations, SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL,  &
                                    soil_theta_to_head
   use meds_pft_params,    only : hydraulics_table_t
   use meds_plant_hydraulics, only : cohort_root_profile
   use meds_soil_params, only : soil_profile_t, soil_params_t, SOIL_BC_DRAIN
   use meds_soil_derive,        only : build_profile_hydraulics
   use meds_test_soil_grid, only : geometric_layer_bottoms
   use meds_soil_water,    only : ground_evap_env_t
   implicit none

   type(hydraulics_table_t)   :: hyd
   type(soil_profile_t)       :: soil
   type(soil_params_t)        :: sopts
   type(ground_evap_env_t)     :: forc
   !----- The water network's tolerances: [fast].rtol and water_atol at their defaults. ---------!
   real(wp), parameter :: WATER_RTOL = 1.0e-3_wp, WATER_ATOL = 1.0e-4_wp

   call make_pft(hyd)
   call make_soil(soil, 10_ik, 1.44_wp, 3.0_wp)
   forc%t_ground = 299.0_wp ; forc%q_air = 0.016_wp ; forc%rho_air = 1.15_wp ; forc%r_aero = 80.0_wp

   call test_arrowhead_solve()
   call test_jacobian()
   call test_conservation()
   call test_plant_only()
   call test_soil_only()
   call test_stiff_exchange()
   call test_infeasible_demand()
   call test_leafless()
   call test_dry_layer_storm()
   call test_saturated()
   call test_accuracy()
   call test_vanishing_stores()
   call test_campbell()
   call test_report('test_spac_solver')

contains

   !=======================================================================================!
   !     Fixtures.                                                                          !
   !=======================================================================================!

   !----- One PFT with the shipped [hydraulics] values and a tropical wood conductivity. -----!
   subroutine make_pft(h)
      type(hydraulics_table_t),   intent(out) :: h
      h%n_pft = 1_ik
      allocate(h%pft(1))
      h%pft(1)%leaf_curve = water_curve_t(pi0=-1.5_wp, elastic_mod=12.0_wp, apoplast_frac=0.30_wp, water_sat=2.0_wp)
      h%pft(1)%wood_curve = water_curve_t(pi0=-1.0_wp, elastic_mod=8.0_wp,  apoplast_frac=0.20_wp, water_sat=1.0_wp)
      h%pft(1)%wood_psi50  = -2.0_wp ; h%pft(1)%wood_kexp = 2.0_wp
      h%pft(1)%wood_kmax   = 3.3_wp  ; h%pft(1)%vessel_curl = 1.5_wp
      h%pft(1)%root_beta = 0.1_wp
      h%pft(1)%root_depth_b1 = 1.1140580_wp
      h%pft(1)%root_depth_b2 = 0.4223014_wp
      h%pft(1)%specific_root_length = 2.0e4_wp
      h%pft(1)%fine_root_radius = 2.5e-4_wp
      call build_hydro_table(h%pft(1)%vuln_table, h%pft(1)%wood_kexp)
   end subroutine make_pft

   !----- BCI's clay (example 04), on a grid whose top layer is 2.6 cm. ----------------------!
   subroutine make_soil(s, nsl, depth, growth)
      type(soil_profile_t), intent(out) :: s
      integer(ik),          intent(in)  :: nsl
      real(wp),             intent(in)  :: depth, growth
      real(wp) :: bottoms(n_soil_layer_max)
      call geometric_layer_bottoms(depth, nsl, growth, bottoms)
      call build_profile_hydraulics(nsl, SOIL_RETENTION_VG, bottoms, 0.4314_wp, 0.0981_wp, 4.54e-6_wp,        &
                                  0.7607_wp, 1.3554_wp, -3.37_wp, s)
   end subroutine make_soil

   !----- A stand of n cohorts, a canopy tree first and smaller ones after, built as             !
   !      build_column_frozen builds them: rooting depth and profile from height, root length    !
   !      from fine-root carbon, root spacing from the patch's root-length density. --------------!
   subroutine make_stand(spac, n, theta, w_pond, leaf_rwc, wood_rwc, e_total, leafless)
      type(spac_solver_t), intent(inout) :: spac
      integer(ik),           intent(in)    :: n
      real(wp),              intent(in)    :: theta(:), w_pond, leaf_rwc, wood_rwc, e_total
      logical, optional,     intent(in)    :: leafless
      real(wp) :: height(n), nplant(n), bleaf(n), bsap(n), broot(n), sap_area(n)
      real(wp) :: frac(soil%n_active), root_len(soil%n_active, n), half(soil%n_active)
      real(wp) :: leaf_w(n), wood_w(n), lai_share(n)
      integer(ik) :: i, k, nsl
      nsl = soil%n_active
      call spac_allocate(spac, nsl, n)
      do i = 1_ik, n
         height(i)   = 28.0_wp / real(i, wp)
         nplant(i)   = 0.02_wp * real(i, wp) ** 1.5_wp
         bleaf(i)    = 4.0_wp / real(i, wp) ** 2
         bsap(i)     = 30.0_wp / real(i, wp) ** 2.5_wp
         broot(i)    = 2.0_wp / real(i, wp) ** 2
         sap_area(i) = 0.05_wp / real(i, wp) ** 2
      end do
      if (present(leafless)) then
         if (leafless) bleaf(n) = 0.0_wp
      end if
      lai_share = nplant * bleaf / max(sum(nplant * bleaf), tiny_num)
      do i = 1_ik, n
         spac%pft(i) = 1_ik ; spac%nplant(i) = nplant(i) ; spac%bleaf(i) = bleaf(i)
         spac%bwood(i) = bsap(i) + broot(i)
         spac%grav(i)  = grav_head * height(i)
         spac%k_xylem(i) = nplant(i) * xylem_max_conductance(hyd%pft(1)%wood_kmax, sap_area(i), height(i),     &
                                                            hyd%pft(1)%vessel_curl)
         spac%e_demand(i) = e_total * lai_share(i)
         spac%r_close(i)  = tissue_psi_to_rwc(2.0_wp * tissue_psi_tlp(hyd%pft(1)%leaf_curve%pi0,                        &
                                        hyd%pft(1)%leaf_curve%elastic_mod), hyd%pft(1)%leaf_curve%pi0,         &
                                        hyd%pft(1)%leaf_curve%elastic_mod)
         call cohort_root_profile(height(i), hyd%pft(1), -soil%soil_layer_z(1:nsl), -soil%soil_layer_z(2:nsl+1), frac)
         root_len(:, i) = broot(i) * frac * hyd%pft(1)%specific_root_length
         leaf_w(i) = bleaf(i) * hyd%pft(1)%leaf_curve%water_sat                                              &
                     * ((1.0_wp - hyd%pft(1)%leaf_curve%apoplast_frac) * leaf_rwc + hyd%pft(1)%leaf_curve%apoplast_frac)
         wood_w(i) = spac%bwood(i) * hyd%pft(1)%wood_curve%water_sat                                          &
                     * ((1.0_wp - hyd%pft(1)%wood_curve%apoplast_frac) * wood_rwc + hyd%pft(1)%wood_curve%apoplast_frac)
      end do
      do k = 1_ik, nsl
         half(k) = 1.0_wp / sqrt(pi * max(sum(nplant * root_len(k, :)) / soil%dz(k), tiny_num))
         do i = 1_ik, n
            spac%g_rhizo(k, i) = rhizosphere_conductance(nplant(i), root_len(k, i), half(k),                  &
                                                         hyd%pft(1)%fine_root_radius)
         end do
      end do
      spac%rtol = WATER_RTOL ; spac%atol = WATER_ATOL ; spac%w_pond_max = sopts%w_pond_max
      spac%t_soil = 299.0_wp ; spac%t_rain = 297.0_wp
      call spac_set_state(spac, soil, hyd, theta, w_pond, 0.0_wp, leaf_w, wood_w)
   end subroutine make_stand

   !----- All the water the network holds [kg/m2]. ---------------------------------------!
   real(wp) function total_water(spac)
      type(spac_solver_t), intent(in) :: spac
      total_water = sum(spac%water_soil) + sum(spac%water_coh)
   end function total_water

   !=======================================================================================!
   subroutine test_arrowhead_solve()
      type(arrowhead_system_t) :: s
      integer(ik), parameter :: NC = 7, NG = 5, MG = 2, NT = NC + NG*MG
      real(wp) :: full(NT, NT), rhs(NT), x_ref(NT), r_c(NC), r_g(MG, NG), x_c(NC), x_g(MG, NG), v
      real(wp) :: add_c(NC), add_g(MG, NG), a_kept(NC, NC)
      integer(ik) :: i, j, g, l, seed
      logical :: ok
      call arrowhead_allocate(s, NC, NG, MG)
      seed = 12345_ik ; full = 0.0_wp
      do j = 1_ik, NC
         do i = 1_ik, NC
            v = pseudo(seed) ; s%a(i, j) = v ; full(i, j) = v
         end do
         s%a(j, j) = s%a(j, j) + 8.0_wp ; full(j, j) = full(j, j) + 8.0_wp
      end do
      do g = 1_ik, NG
         do i = 1_ik, NC            ! each group meets the core through its first unknown
            v = pseudo(seed) ; s%b(i, g) = v ; full(i, NC + (g-1)*MG + 1) = v
            v = pseudo(seed) ; s%c(i, g) = v ; full(NC + (g-1)*MG + 1, i) = v
         end do
         do l = 1_ik, MG
            do j = 1_ik, MG
               v = pseudo(seed) + merge(5.0_wp, 0.0_wp, j == l)
               s%d(l, j, g) = v ; full(NC + (g-1)*MG + l, NC + (g-1)*MG + j) = v
            end do
         end do
      end do
      !----- Additions to the blocks' diagonals, as a Newton matrix's storage terms are added. --------!
      do i = 1_ik, NC
         add_c(i) = 1.5_wp + 0.1_wp * real(i, wp) ; full(i, i) = full(i, i) + add_c(i)
      end do
      do g = 1_ik, NG
         do l = 1_ik, MG
            add_g(l, g) = 0.7_wp ; full(NC + (g-1)*MG + l, NC + (g-1)*MG + l) = full(NC + (g-1)*MG + l, NC + (g-1)*MG + l) + 0.7_wp
         end do
      end do
      do i = 1_ik, NT
         rhs(i) = pseudo(seed)
      end do
      r_c = rhs(1:NC) ; r_g = reshape(rhs(NC+1:NT), [MG, NG])
      call dense_solve(full, rhs, x_ref, NT)
      a_kept = s%a
      call arrowhead_factor(s, add_c, add_g, ok)
      call check_true('arrowhead solve: factorisation succeeds', ok)
      call check_true('arrowhead solve: factoring leaves the blocks as they were', all(s%a == a_kept))
      call arrowhead_solve(s, r_c, r_g, x_c, x_g)
      call check_true('arrowhead solve: equals dense elimination',                                             &
                      maxval(abs([x_c, reshape(x_g, [NG*MG])] - x_ref)) < 1.0e-12_wp * maxval(abs(x_ref)))
   end subroutine test_arrowhead_solve

   !=======================================================================================!
   subroutine test_jacobian()
      type(spac_solver_t) :: spac
      real(wp) :: theta(soil%n_active), err
      theta = 0.36_wp
      call make_stand(spac, 3_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 4.0e-5_wp)
      err = spac_jacobian_error(spac, soil, sopts, forc, hyd, 300.0_wp)
      call check_true('jacobian: wet soil, transpiring stand', err < 1.0e-5_wp, err)
      theta = 0.12_wp ; theta(1) = 0.1 + 0.003_wp
      call make_stand(spac, 3_ik, theta, 0.0_wp, 0.4_wp, 0.35_wp, 2.0e-5_wp)
      err = spac_jacobian_error(spac, soil, sopts, forc, hyd, 300.0_wp)
      call check_true('jacobian: dry soil, stressed stand', err < 1.0e-5_wp, err)
      theta = 0.4314_wp - 1.0e-3_wp ; theta(1) = 0.4314_wp + 1.0e-5_wp
      call make_stand(spac, 2_ik, theta, 12.0_wp, 0.98_wp, 0.97_wp, 0.0_wp)
      err = spac_jacobian_error(spac, soil, sopts, forc, hyd, 300.0_wp)
      call check_true('jacobian: ponded, supersaturated top layer', err < 1.0e-5_wp, err)
      theta = 0.3_wp
      call make_stand(spac, 3_ik, theta, 0.0_wp, 0.9_wp, 0.85_wp, 1.0e-5_wp, leafless=.true.)
      err = spac_jacobian_error(spac, soil, sopts, forc, hyd, 300.0_wp)
      call check_true('jacobian: a leafless cohort', err < 1.0e-5_wp, err)
   end subroutine test_jacobian

   !=======================================================================================!
   subroutine test_conservation()
      type(spac_solver_t) :: spac
      real(wp) :: theta(soil%n_active), w0, w1, dt, hour, resid, worst, flows, worst_coh
      real(wp), allocatable :: wc0(:), wc1(:)
      integer(ik) :: step, nfail
      theta = 0.33_wp
      call make_stand(spac, 4_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp)
      dt = 900.0_wp ; worst = 0.0_wp ; worst_coh = 0.0_wp ; nfail = 0_ik
      do step = 1_ik, 192
         hour = mod(real(step, wp) * dt / 3600.0_wp, 24.0_wp)
         spac%e_demand = 5.0e-5_wp * max(0.0_wp, sin(pi * (hour - 6.0_wp) / 12.0_wp))                        &
                        * spac%nplant * spac%bleaf / sum(spac%nplant * spac%bleaf)
         spac%rain = 0.0_wp
         if (step >= 60 .and. step < 68) spac%rain = 2.0e-2_wp      ! 72 mm/h for two hours
         if (step >= 150 .and. step < 160) spac%rain = 1.0e-3_wp    ! a long drizzle
         w0 = total_water(spac) ; wc0 = sum(spac%water_coh, dim=1)
         call spac_advance(spac, soil, sopts, forc, hyd, dt)
         w1 = total_water(spac) ; wc1 = sum(spac%water_coh, dim=1)
         !----- Each cohort's leaf + wood water changes by its own uptake less its transpiration. ------!
         worst_coh = max(worst_coh, maxval(abs((wc1 - wc0) - (sum(spac%uptake, dim=1) - spac%transp)))    &
                                    / max(maxval(wc0), 1.0e-6_wp))
         flows = spac%rain * dt - spac%drain - spac%evap - spac%spill - sum(spac%transp)
         resid = abs((w1 - w0) - flows)
         worst = max(worst, resid / max(w0, 1.0_wp))
         nfail = nfail + spac%nfail
      end do
      call check_true('conservation: every step closes to round-off', worst < 1.0e-13_wp, worst)
      call check_true('conservation: every cohort''s tissue water closes to round-off', worst_coh < 1.0e-12_wp, &
                      worst_coh)
      call check_int('conservation: no sub-step forced through', nfail, 0_ik)
   end subroutine test_conservation

   !=======================================================================================!
   !----- The soil held at its potential (prescribed core nodes: a plant-only network). Under a      !
   !      steady demand the plant reaches its steady state, where every path carries the            !
   !      transpiration: root uptake = sapflow = transpiration, and the stores stop changing. --------!
   subroutine test_plant_only()
      type(spac_solver_t) :: spac
      real(wp) :: theta(soil%n_active), w_plant0, w_plant1, e, u, sap
      integer(ik) :: step
      theta = 0.28_wp
      call make_stand(spac, 1_ik, theta, 0.0_wp, 0.95_wp, 0.92_wp, 3.0e-5_wp)
      spac%soil_fixed = .true.
      spac%rtol = 1.0e-8_wp ; spac%atol = 1.0e-9_wp
      do step = 1_ik, 96_ik
         w_plant0 = sum(spac%water_coh)
         call spac_advance(spac, soil, sopts, forc, hyd, 900.0_wp)
         w_plant1 = sum(spac%water_coh)
      end do
      e = spac%transp(1) ; u = sum(spac%uptake(:, 1)) ; sap = spac%sapflow(1)
      call check_close('plant only: at steady state root uptake = transpiration', u, e, 1.0e-6_wp * e)
      call check_close('plant only: at steady state sapflow = transpiration', sap, e, 1.0e-6_wp * e)
      call check_close('plant only: the plant water stops changing', w_plant1 - w_plant0, 0.0_wp, 1.0e-8_wp * w_plant1)
      call check_true('plant only: the leaf sits below the wood by more than the crown head',                   &
                      spac%psi_coh(LEAF, 1) < spac%psi_coh(WOOD, 1) - spac%grav(1),                &
                      spac%psi_coh(LEAF, 1))
   end subroutine test_plant_only

   !=======================================================================================!
   !----- No cohorts, a sealed bottom, no rain and no evaporation: the column relaxes to hydrostatic !
   !      equilibrium, the total head (matric head plus elevation) the same in every layer. ---------!
   subroutine test_soil_only()
      type(spac_solver_t)    :: spac
      type(soil_params_t)    :: sealed
      type(ground_evap_env_t) :: still
      real(wp) :: theta(soil%n_active), head(soil%n_active), w0
      integer(ik) :: step, nsl
      nsl = soil%n_active
      theta = 0.40_wp ; theta(1:3) = 0.425_wp          ! wet: BCI's clay settles in days, not years
      sealed = sopts ; sealed%water_bottom_bc = SOIL_BC_DRAIN ; sealed%drain_factor = 0.0_wp
      still = forc ; still%q_air = 0.05_wp                 ! air wetter than any surface: no evaporation
      call make_stand(spac, 0_ik, theta, 0.0_wp, 0.9_wp, 0.9_wp, 0.0_wp)
      w0 = total_water(spac)
      do step = 1_ik, 96_ik * 30_ik
         call spac_advance(spac, soil, sealed, still, hyd, 900.0_wp)
      end do
      head = spac%psi_soil(2:nsl+1) / grav_head + soil%z_node(1:nsl)
      call check_true('soil only: hydrostatic equilibrium (total head uniform to 1 cm)',                        &
                      maxval(head) - minval(head) < 1.0e-2_wp, maxval(head) - minval(head))
      call check_close('soil only: sealed column keeps its water', total_water(spac), w0, 1.0e-11_wp * w0)
   end subroutine test_soil_only

   !=======================================================================================!
   subroutine test_stiff_exchange()
      type(spac_solver_t) :: spac
      real(wp) :: theta(soil%n_active)
      integer(ik) :: step, nsub
      theta = 0.42_wp
      call make_stand(spac, 4_ik, theta, 0.0_wp, 0.97_wp, 0.95_wp, 0.0_wp)
      spac%g_rhizo = 50.0_wp * spac%g_rhizo                ! a fifty-fold root density: a sub-second exchange
      nsub = 0_ik
      do step = 1_ik, 8_ik
         spac%e_demand = merge(6.0e-5_wp, 0.0_wp, mod(step, 2_ik) == 0_ik) * spac%nplant / sum(spac%nplant)
         call spac_advance(spac, soil, sopts, forc, hyd, 900.0_wp)
         nsub = nsub + spac%nsub + spac%nreject
      end do
      call check_true('stiff exchange: demand jumps cost at most 6 sub-steps a step', nsub <= 48_ik, real(nsub, wp))
   end subroutine test_stiff_exchange

   !=======================================================================================!
   subroutine test_infeasible_demand()
      type(spac_solver_t) :: spac
      real(wp) :: theta(soil%n_active), demand
      theta = 0.15_wp
      call make_stand(spac, 1_ik, theta, 0.0_wp, 0.6_wp, 0.5_wp, 0.0_wp)
      spac%bleaf(1) = 1.0e-3_wp                          ! a few grams of leaf
      spac%k_xylem(1) = 1.0e-3_wp * spac%k_xylem(1)        ! and a xylem that can hardly feed it
      call spac_set_state(spac, soil, hyd, theta, 0.0_wp, 0.0_wp,                                      &
                               [1.0e-3_wp * 2.0_wp * 0.75_wp], [spac%water_coh(WOOD, 1) / spac%nplant(1)])
      demand = 1.0e-4_wp
      spac%e_demand(1) = demand
      call spac_advance(spac, soil, sopts, forc, hyd, 900.0_wp)
      call check_int('infeasible demand: no sub-step forced through', spac%nfail, 0_ik)
      call check_true('infeasible demand: the leaf stays above empty', spac%tau_coh(LEAF, 1) > 0.0_wp,     &
                      spac%tau_coh(LEAF, 1))
      call check_true('infeasible demand: transpiration falls short of the demand',                          &
                      spac%transp(1) < 0.5_wp * demand * 900.0_wp, spac%transp(1) / (demand * 900.0_wp))
      call check_true('infeasible demand: leaf potential is finite', abs(spac%psi_coh(LEAF, 1)) < 1.0e6_wp,  &
                      spac%psi_coh(LEAF, 1))
   end subroutine test_infeasible_demand

   !=======================================================================================!
   subroutine test_leafless()
      type(spac_solver_t) :: spac
      real(wp) :: theta(soil%n_active)
      theta = 0.30_wp
      call make_stand(spac, 3_ik, theta, 0.0_wp, 0.9_wp, 0.8_wp, 2.0e-5_wp, leafless=.true.)
      spac%e_demand(3) = 0.0_wp
      call spac_advance(spac, soil, sopts, forc, hyd, 900.0_wp)
      call check_close('leafless: leaf at the wood less the crown head', spac%psi_coh(LEAF, 3),            &
                       spac%psi_coh(WOOD, 3) - spac%grav(3), 1.0e-6_wp)
      call check_close('leafless: its empty leaf store stays empty', spac%water_coh(LEAF, 3), 0.0_wp, 1.0e-12_wp)
   end subroutine test_leafless

   !=======================================================================================!
   !----- A dry-layer storm, from a BCI regrowth stand: the top 2.6 cm a hair above residual, a !
   !      stand whose wood is half empty and refilling, then 30 mm/h of rain for two steps. ----!
   subroutine test_dry_layer_storm()
      type(spac_solver_t) :: spac
      real(wp) :: theta(soil%n_active), w0, flows, theta_min
      integer(ik) :: step, nsl
      logical :: finite
      nsl = soil%n_active
      theta = 0.16_wp ; theta(1) = 0.0981_wp + 2.0e-3_wp
      call make_stand(spac, 6_ik, theta, 0.0_wp, 0.45_wp, 0.40_wp, 0.0_wp)
      w0 = total_water(spac) ; flows = 0.0_wp
      finite = .true. ; theta_min = 1.0_wp
      do step = 1_ik, 6_ik
         spac%rain = merge(30.0_wp / 3600.0_wp, 0.0_wp, step == 3 .or. step == 4)
         call spac_advance(spac, soil, sopts, forc, hyd, 900.0_wp)
         flows = flows + spac%rain * 900.0_wp - spac%drain - spac%evap - spac%spill - sum(spac%transp)
         finite = finite .and. all(spac%water_soil == spac%water_soil) .and. all(spac%psi_coh == spac%psi_coh)
         theta_min = min(theta_min, minval(spac%water_soil(2:nsl+1) / (rho_h2o * soil%dz(1:nsl)) - soil%theta_res(1:nsl)))
         call check_int('dry-layer storm: no sub-step forced through', spac%nfail, 0_ik)
      end do
      call check_true('dry-layer storm: every state finite', finite)
      call check_true('dry-layer storm: no layer below residual', theta_min > 0.0_wp, theta_min)
      call check_close('dry-layer storm: water conserved', total_water(spac) - w0, flows, 1.0e-10_wp * w0)
      call check_true('dry-layer storm: the wood refilled', spac%tau_coh(WOOD, 1) > 0.40_wp,                &
                      spac%tau_coh(WOOD, 1))
   end subroutine test_dry_layer_storm

   !=======================================================================================!
   subroutine test_saturated()
      type(spac_solver_t) :: spac
      real(wp) :: theta(soil%n_active), w0, flows, spill, over
      integer(ik) :: step, nsl
      nsl = soil%n_active
      theta = 0.4314_wp - 1.0e-3_wp
      call make_stand(spac, 2_ik, theta, 0.0_wp, 0.97_wp, 0.95_wp, 1.0e-5_wp)
      w0 = total_water(spac) ; flows = 0.0_wp ; spill = 0.0_wp ; over = 0.0_wp
      do step = 1_ik, 8_ik
         spac%rain = merge(50.0_wp / 3600.0_wp, 0.0_wp, step <= 4)
         call spac_advance(spac, soil, sopts, forc, hyd, 900.0_wp)
         flows = flows + spac%rain * 900.0_wp - spac%drain - spac%evap - spac%spill - sum(spac%transp)
         spill = spill + spac%spill
         over = max(over, maxval(spac%water_soil(2:nsl+1) / (rho_h2o * soil%dz(1:nsl)) - soil%theta_sat(1:nsl)))
         call check_int('saturated: no sub-step forced through', spac%nfail, 0_ik)
      end do
      call check_true('saturated: the pond spilled', spill > 0.0_wp, spill)
      call check_true('saturated: no layer more than 1e-3 above saturation', over < 1.0e-3_wp, over)
      call check_true('saturated: pond at or below its capacity', spac%water_soil(POND) <= sopts%w_pond_max,  &
                      spac%water_soil(POND))
      call check_close('saturated: water conserved', total_water(spac) - w0, flows, 1.0e-10_wp * w0)
   end subroutine test_saturated

   !=======================================================================================!
   !----- A transpiring day with a midday storm at the default tolerances, against the same day at  !
   !      tolerances a thousand times tighter. ------------------------------------------------!
   subroutine test_accuracy()
      type(spac_solver_t) :: spac, ref
      real(wp) :: theta(soil%n_active), hour, worst_theta, worst_rwc
      integer(ik) :: step, nsl, isub
      nsl = soil%n_active
      theta = 0.30_wp
      call make_stand(spac, 3_ik, theta, 0.0_wp, 0.9_wp, 0.85_wp, 0.0_wp)
      call make_stand(ref, 3_ik, theta, 0.0_wp, 0.9_wp, 0.85_wp, 0.0_wp)
      ref%rtol = 1.0e-3_wp * spac%rtol ; ref%atol = 1.0e-3_wp * spac%atol
      isub = 0_ik
      do step = 1_ik, 96_ik
         hour = real(step, wp) * 0.25_wp
         spac%e_demand = 5.0e-5_wp * max(0.0_wp, sin(pi * (hour - 6.0_wp) / 12.0_wp))                        &
                        * spac%nplant * spac%bleaf / sum(spac%nplant * spac%bleaf)
         ref%e_demand = spac%e_demand
         spac%rain = merge(1.0e-2_wp, 0.0_wp, step >= 50 .and. step < 54) ; ref%rain = spac%rain
         call spac_advance(spac, soil, sopts, forc, hyd, 900.0_wp)
         call spac_advance(ref, soil, sopts, forc, hyd, 900.0_wp)
         isub = isub + spac%nsub
      end do
      worst_theta = maxval(abs(spac%water_soil(2:nsl+1) - ref%water_soil(2:nsl+1)) / (rho_h2o * soil%dz(1:nsl)))
      worst_rwc   = maxval(abs(spac%tau_coh - ref%tau_coh))
      call check_true('accuracy: soil water within 10 atol of the tight run', worst_theta < 10.0_wp * WATER_ATOL,  &
                      worst_theta)
      call check_true('accuracy: tissue water within 1e-2 RWC of the tight run', worst_rwc < 1.0e-2_wp, worst_rwc)
      print '(a,i6)', '  info : sub-steps over the day at the default tolerances: ', isub
   end subroutine test_accuracy

   !=======================================================================================!
   !----- Two stores a BCI regrowth stand reaches: a cohort with 1e-22 kgC of leaves,              !
   !      whose leaf store holds 1e-26 kg/m2 against its plant's 1e-9 tolerance, and a pond holding a !
   !      trace of water over a dry top layer in light rain, which drains tens of millimetres an hour. !
   !      A day of them: no sub-step forced, the leaf never above saturation, the pond never below  !
   !      empty. --------------------------------------------------------------------------------!
   subroutine test_vanishing_stores()
      type(spac_solver_t) :: spac
      real(wp) :: theta(soil%n_active), leaf_w(3), wood_w(3), w0, flows, rwc_max, pond_min
      integer(ik) :: step, nfail
      theta = 0.16_wp
      call make_stand(spac, 3_ik, theta, 0.0_wp, 0.9_wp, 0.85_wp, 0.0_wp)
      spac%bleaf(3) = 1.3e-22_wp
      leaf_w = spac%water_coh(LEAF, :) / spac%nplant ; wood_w = spac%water_coh(WOOD, :) / spac%nplant
      leaf_w(3) = 2.0_wp * spac%bleaf(3) * 0.9_wp
      call spac_set_state(spac, soil, hyd, theta, 3.0e-8_wp, 0.0_wp, leaf_w, wood_w)
      w0 = total_water(spac) ; flows = 0.0_wp ; nfail = 0_ik ; rwc_max = 0.0_wp ; pond_min = 0.0_wp
      do step = 1_ik, 96_ik
         spac%rain = merge(1.0e-6_wp, 1.2e-3_wp, step < 40)
         spac%e_demand(1:2) = 2.0e-5_wp ; spac%e_demand(3) = 1.0e-32_wp
         call spac_advance(spac, soil, sopts, forc, hyd, 900.0_wp)
         flows = flows + spac%rain * 900.0_wp - spac%drain - spac%evap - spac%spill - sum(spac%transp)
         nfail = nfail + spac%nfail
         rwc_max  = max(rwc_max, spac%tau_coh(LEAF, 3))
         pond_min = min(pond_min, spac%water_soil(POND))
      end do
      call check_int('vanishing stores: no sub-step forced through', nfail, 0_ik)
      call check_true('vanishing stores: the near-leafless leaf stays at or below saturation', rwc_max <= 1.0_wp, rwc_max)
      call check_true('vanishing stores: the pond never holds less than nothing', pond_min >= 0.0_wp, pond_min)
      call check_close('vanishing stores: water conserved', total_water(spac) - w0, flows, 1.0e-10_wp * w0)
   end subroutine test_vanishing_stores

   !=======================================================================================!
   !----- Two transpiring days without rain on a Campbell clay (psi_sat = -0.26 m, b = 5.65), on   !
   !      the van Genuchten clay's layers. The column carries one parameter pair whose meaning the !
   !      family sets, so the network must read (psi_sat, b) through it: each layer's potential is !
   !      the Campbell curve's at its water content, finite and below zero, the stand draws water  !
   !      from the soil, and every step closes. No layer reaches saturation, where the curve reads !
   !      zero and the network a head above air entry. The fixture's soil is swapped for the test  !
   !      and put back after it. ------------------------------------------------------------!
   subroutine test_campbell()
      type(spac_solver_t)  :: spac
      type(soil_profile_t) :: vg_soil
      real(wp), parameter  :: THETA_SAT = 0.44_wp, THETA_RES = 0.0_wp, PSI_SAT = -0.26_wp, B = 5.65_wp
      real(wp) :: bottoms(n_soil_layer_max), theta(soil%n_active), psi_curve(soil%n_active)
      real(wp) :: w0, w1, dt, hour, flows, worst, worst_psi, uptake
      integer(ik) :: step, nfail, nsl
      logical :: below_zero
      vg_soil = soil
      nsl = vg_soil%n_active
      call geometric_layer_bottoms(1.44_wp, nsl, 3.0_wp, bottoms)
      call build_profile_hydraulics(nsl, SOIL_RETENTION_CAMPBELL, bottoms, THETA_SAT, THETA_RES, 4.54e-6_wp,    &
                                  PSI_SAT, B, -3.37_wp, soil)
      theta = 0.30_wp
      call make_stand(spac, 4_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp)
      dt = 900.0_wp ; worst = 0.0_wp ; worst_psi = 0.0_wp ; uptake = 0.0_wp ; nfail = 0_ik
      below_zero = .true.
      do step = 1_ik, 192
         hour = mod(real(step, wp) * dt / 3600.0_wp, 24.0_wp)
         spac%e_demand = 5.0e-5_wp * max(0.0_wp, sin(pi * (hour - 6.0_wp) / 12.0_wp))                        &
                        * spac%nplant * spac%bleaf / sum(spac%nplant * spac%bleaf)
         spac%rain = 0.0_wp
         w0 = total_water(spac)
         call spac_advance(spac, soil, sopts, forc, hyd, dt)
         w1 = total_water(spac)
         flows = -spac%drain - spac%evap - spac%spill - sum(spac%transp)
         worst = max(worst, abs((w1 - w0) - flows) / max(w0, 1.0_wp))
         nfail = nfail + spac%nfail
         uptake = uptake + sum(spac%uptake)
         !----- Each layer's potential against the Campbell curve at the water it holds, which the  !
         !      sub-step commits from its flows: the two agree to the Newton tolerance. ------------!
         theta = spac%water_soil(2:nsl+1) / (rho_h2o * soil%dz(1:nsl))
         psi_curve = grav_head * soil_theta_to_head(SOIL_RETENTION_CAMPBELL, theta, THETA_SAT, THETA_RES,    &
                                                     PSI_SAT, B)
         worst_psi = max(worst_psi, maxval(abs(spac%psi_soil(2:nsl+1) - psi_curve) / abs(psi_curve)))
         below_zero = below_zero .and. all(spac%psi_soil(2:nsl+1) < 0.0_wp)                                 &
                                 .and. all(spac%psi_soil(2:nsl+1) > -huge(1.0_wp))
      end do
      soil = vg_soil
      call check_true('Campbell: every soil potential is finite and below zero', below_zero)
      call check_true('Campbell: every soil potential is the Campbell curve''s', worst_psi < 1.0e-5_wp, worst_psi)
      call check_true('Campbell: the stand takes up water from the soil', uptake > 0.0_wp, uptake)
      call check_true('Campbell: every step closes to round-off', worst < 1.0e-13_wp, worst)
      call check_int('Campbell: no sub-step forced through', nfail, 0_ik)
   end subroutine test_campbell

   !=======================================================================================!
   !     Small helpers.                                                                     !
   !=======================================================================================!
   !----- A linear congruential generator, in 64-bit arithmetic so the product cannot overflow. ---!
   real(wp) function pseudo(seed)
      integer(ik), intent(inout) :: seed
      integer(int64) :: s
      s = mod(1103515245_int64 * int(seed, int64) + 12345_int64, 2147483648_int64)
      seed = int(s, ik)
      pseudo = real(s, wp) / 2147483648.0_wp - 0.5_wp
   end function pseudo

   subroutine dense_solve(a_in, b, x, n)
      integer(ik), intent(in)  :: n
      real(wp),    intent(in)  :: a_in(n, n), b(n)
      real(wp),    intent(out) :: x(n)
      real(wp) :: a(n, n+1), row(n+1), f
      integer(ik) :: i, j, p
      a(:, 1:n) = a_in ; a(:, n+1) = b
      do j = 1_ik, n
         p = j
         do i = j + 1_ik, n
            if (abs(a(i, j)) > abs(a(p, j))) p = i
         end do
         row = a(j, :) ; a(j, :) = a(p, :) ; a(p, :) = row
         do i = j + 1_ik, n
            f = a(i, j) / a(j, j)
            a(i, :) = a(i, :) - f * a(j, :)
         end do
      end do
      do i = n, 1_ik, -1_ik
         x(i) = (a(i, n+1) - sum(a(i, i+1:n) * x(i+1:n))) / a(i, i)
      end do
   end subroutine dense_solve

end program test_spac_solver
