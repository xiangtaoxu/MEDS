! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_water_network -- the soil-plant water network solver (MEDS_SPAC_WATER_NETWORK_PLAN.md P2). !
!                                                                                          !
!   1. GROUPED SOLVE   : the Schur solve equals dense elimination on a random grouped system.   !
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
!   9. DRY-LAYER STORM : the 1994-04-20 pattern -- a thin top layer near residual, a stand         !
!                        refilling its wood, then a storm: converges, no store leaves its range.   !
!  10. SATURATED       : rain on a saturated column ponds, spills, and conserves.                   !
!  11. ACCURACY        : the default tolerances track a tight-tolerance run.                       !
!==========================================================================================!
program test_water_network
   use meds_test_assert, only : check_true, check_close, check_int, test_report
   use, intrinsic :: iso_fortran_env, only : int64
   use meds_kinds,       only : wp, ik
   use meds_constants,   only : grav_head, rho_h2o, pi, tiny_num
   use meds_grouped_solve, only : grouped_system_t, grouped_allocate, grouped_factor, grouped_solve
   use meds_water_network, only : water_network_t, water_network_allocate, water_network_start,            &
                                  water_network_advance, water_network_jacobian_error,                     &
                                  WATER_NODE_POND, WATER_WOOD, WATER_LEAF
   use meds_water_paths,   only : xylem_max_conductance, root_conductance_geometry
   use meds_water_storage, only : tissue_store, soil_store
   use meds_water_retention, only : SOIL_RETENTION_VG, water_curve_t, rwc_from_psi, pv_psi_tlp,             &
                                    soil_psi_from_theta, soil_hydr_cond_from_theta
   use meds_hydr_lib,      only : build_hydro_table
   use meds_plant_types,   only : hydro_params_table_t
   use meds_plant_hydraulics, only : cohort_root_profile
   use meds_column_params, only : soil_params_t, build_soil_hydr_params
   use meds_soil_types,    only : chydro_forcing_t
   use meds_biophysics_opts, only : soil_opts_t, SOIL_BC_BEDROCK
   implicit none

   type(hydro_params_table_t) :: hyd
   type(soil_params_t)        :: soil
   type(soil_opts_t)          :: sopts
   type(chydro_forcing_t)     :: forc

   call make_pft(hyd)
   call make_soil(soil, 10_ik, 1.44_wp, 3.0_wp)
   forc%t_ground = 299.0_wp ; forc%q_air = 0.016_wp ; forc%rho_air = 1.15_wp ; forc%r_aero = 80.0_wp

   call test_grouped_solve()
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
   call test_report('test_water_network')

contains

   !=======================================================================================!
   !     Fixtures.                                                                          !
   !=======================================================================================!

   !----- One PFT with the shipped [hydraulics] defaults and a tropical wood conductivity. ----!
   subroutine make_pft(h)
      type(hydro_params_table_t), intent(out) :: h
      h%n_pft = 1_ik
      allocate(h%pft(1))
      h%pft(1)%leaf_curve = water_curve_t(pi0=-1.5_wp, elastic_mod=12.0_wp, apoplast_frac=0.30_wp, water_sat=2.0_wp)
      h%pft(1)%wood_curve = water_curve_t(pi0=-1.0_wp, elastic_mod=8.0_wp,  apoplast_frac=0.20_wp, water_sat=1.0_wp)
      h%pft(1)%wood_psi50  = -2.0_wp ; h%pft(1)%wood_kexp = 2.0_wp
      h%pft(1)%wood_kmax   = 3.3_wp  ; h%pft(1)%vessel_curl = 1.5_wp
      call build_hydro_table(h%pft(1)%vuln_table, h%pft(1)%wood_kexp)
   end subroutine make_pft

   !----- BCI's clay (example 04), on a grid whose top layer is 2.6 cm. ----------------------!
   subroutine make_soil(s, nsl, depth, growth)
      type(soil_params_t), intent(out) :: s
      integer(ik),         intent(in)  :: nsl
      real(wp),            intent(in)  :: depth, growth
      call build_soil_hydr_params(nsl, SOIL_RETENTION_VG, depth, growth, 0.4314_wp, 0.0981_wp, 4.54e-6_wp,    &
                                  0.7607_wp, 1.3554_wp, -3.37_wp, s)
   end subroutine make_soil

   !----- A stand of n cohorts, a canopy tree first and smaller ones after, built as             !
   !      build_column_frozen builds them: rooting depth and profile from height, root length    !
   !      from fine-root carbon, root spacing from the patch's root-length density. --------------!
   subroutine make_stand(net, n, theta, w_pond, leaf_rwc, wood_rwc, e_total, leafless)
      type(water_network_t), intent(inout) :: net
      integer(ik),           intent(in)    :: n
      real(wp),              intent(in)    :: theta(:), w_pond, leaf_rwc, wood_rwc, e_total
      logical, optional,     intent(in)    :: leafless
      real(wp) :: height(n), nplant(n), bleaf(n), bsap(n), broot(n), sap_area(n)
      real(wp) :: frac(soil%n_active), root_len(soil%n_active, n), half(soil%n_active)
      real(wp) :: leaf_w(n), wood_w(n), lai_share(n)
      integer(ik) :: i, k, nsl
      nsl = soil%n_active
      call water_network_allocate(net, nsl, n)
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
         net%pft(i) = 1_ik ; net%nplant(i) = nplant(i) ; net%bleaf(i) = bleaf(i)
         net%bwood(i) = bsap(i) + broot(i)
         net%grav(i)  = grav_head * height(i)
         net%k_xylem(i) = nplant(i) * xylem_max_conductance(hyd%pft(1)%wood_kmax, sap_area(i), height(i),     &
                                                            hyd%pft(1)%vessel_curl)
         net%e_demand(i) = e_total * lai_share(i)
         net%r_close(i)  = rwc_from_psi(2.0_wp * pv_psi_tlp(hyd%pft(1)%leaf_curve%pi0,                        &
                                        hyd%pft(1)%leaf_curve%elastic_mod), hyd%pft(1)%leaf_curve%pi0,         &
                                        hyd%pft(1)%leaf_curve%elastic_mod)
         call cohort_root_profile(height(i), hyd%pft(1), -soil%soil_layer_z(1:nsl), -soil%soil_layer_z(2:nsl+1), frac)
         root_len(:, i) = broot(i) * frac * hyd%pft(1)%specific_root_length
         leaf_w(i) = bleaf(i) * hyd%pft(1)%leaf_curve%water_sat                                              &
                     * ((1.0_wp - hyd%pft(1)%leaf_curve%apoplast_frac) * leaf_rwc + hyd%pft(1)%leaf_curve%apoplast_frac)
         wood_w(i) = net%bwood(i) * hyd%pft(1)%wood_curve%water_sat                                          &
                     * ((1.0_wp - hyd%pft(1)%wood_curve%apoplast_frac) * wood_rwc + hyd%pft(1)%wood_curve%apoplast_frac)
      end do
      do k = 1_ik, nsl
         half(k) = 1.0_wp / sqrt(pi * max(sum(nplant * root_len(k, :)) / soil%dz(k), tiny_num))
         do i = 1_ik, n
            net%g_root(k, i) = root_conductance_geometry(nplant(i), root_len(k, i), half(k),                  &
                                                         hyd%pft(1)%fine_root_radius)
         end do
      end do
      net%rtol = sopts%rtol ; net%atol = sopts%atol ; net%w_pond_max = sopts%w_pond_max
      net%t_soil = 299.0_wp ; net%t_rain = 297.0_wp
      call water_network_start(net, soil, hyd, theta, w_pond, 0.0_wp, leaf_w, wood_w)
   end subroutine make_stand

   !----- All the water the network holds [kg/m2]. ---------------------------------------!
   real(wp) function total_water(net)
      type(water_network_t), intent(in) :: net
      total_water = sum(net%w_core) + sum(net%w_group)
   end function total_water

   !=======================================================================================!
   subroutine test_grouped_solve()
      type(grouped_system_t) :: s
      integer(ik), parameter :: nc = 7, ng = 5, mg = 2, nt = nc + ng*mg
      real(wp) :: full(nt, nt), rhs(nt), x_ref(nt), r_c(nc), r_g(mg, ng), x_c(nc), x_g(mg, ng), v
      real(wp) :: add_c(nc), add_g(mg, ng), a_kept(nc, nc)
      integer(ik) :: i, j, g, l, seed
      logical :: ok
      call grouped_allocate(s, nc, ng, mg)
      seed = 12345_ik ; full = 0.0_wp
      do j = 1_ik, nc
         do i = 1_ik, nc
            v = pseudo(seed) ; s%a(i, j) = v ; full(i, j) = v
         end do
         s%a(j, j) = s%a(j, j) + 8.0_wp ; full(j, j) = full(j, j) + 8.0_wp
      end do
      do g = 1_ik, ng
         do i = 1_ik, nc            ! each group meets the core through its first unknown
            v = pseudo(seed) ; s%b(i, g) = v ; full(i, nc + (g-1)*mg + 1) = v
            v = pseudo(seed) ; s%c(i, g) = v ; full(nc + (g-1)*mg + 1, i) = v
         end do
         do l = 1_ik, mg
            do j = 1_ik, mg
               v = pseudo(seed) + merge(5.0_wp, 0.0_wp, j == l)
               s%d(l, j, g) = v ; full(nc + (g-1)*mg + l, nc + (g-1)*mg + j) = v
            end do
         end do
      end do
      !----- Additions to the blocks' diagonals, as a Newton matrix's storage terms are added. --------!
      do i = 1_ik, nc
         add_c(i) = 1.5_wp + 0.1_wp * real(i, wp) ; full(i, i) = full(i, i) + add_c(i)
      end do
      do g = 1_ik, ng
         do l = 1_ik, mg
            add_g(l, g) = 0.7_wp ; full(nc + (g-1)*mg + l, nc + (g-1)*mg + l) = full(nc + (g-1)*mg + l, nc + (g-1)*mg + l) + 0.7_wp
         end do
      end do
      do i = 1_ik, nt
         rhs(i) = pseudo(seed)
      end do
      r_c = rhs(1:nc) ; r_g = reshape(rhs(nc+1:nt), [mg, ng])
      call dense_solve(full, rhs, x_ref, nt)
      a_kept = s%a
      call grouped_factor(s, add_c, add_g, ok)
      call check_true('grouped solve: factorisation succeeds', ok)
      call check_true('grouped solve: factoring leaves the blocks as they were', all(s%a == a_kept))
      call grouped_solve(s, r_c, r_g, x_c, x_g)
      call check_true('grouped solve: equals dense elimination',                                             &
                      maxval(abs([x_c, reshape(x_g, [ng*mg])] - x_ref)) < 1.0e-12_wp * maxval(abs(x_ref)))
   end subroutine test_grouped_solve

   !=======================================================================================!
   subroutine test_jacobian()
      type(water_network_t) :: net
      real(wp) :: theta(soil%n_active), err
      theta = 0.36_wp
      call make_stand(net, 3_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 4.0e-5_wp)
      err = water_network_jacobian_error(net, soil, sopts, forc, hyd, 300.0_wp)
      call check_true('jacobian: wet soil, transpiring stand', err < 1.0e-5_wp, err)
      theta = 0.12_wp ; theta(1) = 0.1 + 0.003_wp
      call make_stand(net, 3_ik, theta, 0.0_wp, 0.4_wp, 0.35_wp, 2.0e-5_wp)
      err = water_network_jacobian_error(net, soil, sopts, forc, hyd, 300.0_wp)
      call check_true('jacobian: dry soil, stressed stand', err < 1.0e-5_wp, err)
      theta = 0.4314_wp - 1.0e-3_wp ; theta(1) = 0.4314_wp + 1.0e-5_wp
      call make_stand(net, 2_ik, theta, 12.0_wp, 0.98_wp, 0.97_wp, 0.0_wp)
      err = water_network_jacobian_error(net, soil, sopts, forc, hyd, 300.0_wp)
      call check_true('jacobian: ponded, supersaturated top layer', err < 1.0e-5_wp, err)
      theta = 0.3_wp
      call make_stand(net, 3_ik, theta, 0.0_wp, 0.9_wp, 0.85_wp, 1.0e-5_wp, leafless=.true.)
      err = water_network_jacobian_error(net, soil, sopts, forc, hyd, 300.0_wp)
      call check_true('jacobian: a leafless cohort', err < 1.0e-5_wp, err)
   end subroutine test_jacobian

   !=======================================================================================!
   subroutine test_conservation()
      type(water_network_t) :: net
      real(wp) :: theta(soil%n_active), w0, w1, dt, hour, resid, worst, flows
      integer(ik) :: step, nfail
      theta = 0.33_wp
      call make_stand(net, 4_ik, theta, 0.0_wp, 0.95_wp, 0.9_wp, 0.0_wp)
      dt = 900.0_wp ; worst = 0.0_wp ; nfail = 0_ik
      do step = 1_ik, 192
         hour = mod(real(step, wp) * dt / 3600.0_wp, 24.0_wp)
         net%e_demand = 5.0e-5_wp * max(0.0_wp, sin(pi * (hour - 6.0_wp) / 12.0_wp))                        &
                        * net%nplant * net%bleaf / sum(net%nplant * net%bleaf)
         net%rain = 0.0_wp
         if (step >= 60 .and. step < 68) net%rain = 2.0e-2_wp      ! 72 mm/h for two hours
         if (step >= 150 .and. step < 160) net%rain = 1.0e-3_wp    ! a long drizzle
         w0 = total_water(net)
         call water_network_advance(net, soil, sopts, forc, hyd, dt)
         w1 = total_water(net)
         flows = net%rain * dt - net%drain - net%evap - net%spill - sum(net%transp)
         resid = abs((w1 - w0) - flows)
         worst = max(worst, resid / max(w0, 1.0_wp))
         nfail = nfail + net%nfail
      end do
      call check_true('conservation: every step closes to round-off', worst < 1.0e-13_wp, worst)
      call check_int('conservation: no sub-step forced through', nfail, 0_ik)
   end subroutine test_conservation

   !=======================================================================================!
   !----- The soil held at its potential (prescribed core nodes: a plant-only network). Under a      !
   !      steady demand the plant reaches its steady state, where every path carries the            !
   !      transpiration: root uptake = sapflow = transpiration, and the stores stop changing. --------!
   subroutine test_plant_only()
      type(water_network_t) :: net
      real(wp) :: theta(soil%n_active), w_plant0, w_plant1, e, u, sap
      integer(ik) :: step
      theta = 0.28_wp
      call make_stand(net, 1_ik, theta, 0.0_wp, 0.95_wp, 0.92_wp, 3.0e-5_wp)
      net%core_fixed = .true.
      net%rtol = 1.0e-8_wp ; net%atol = 1.0e-9_wp
      do step = 1_ik, 96_ik
         w_plant0 = sum(net%w_group)
         call water_network_advance(net, soil, sopts, forc, hyd, 900.0_wp)
         w_plant1 = sum(net%w_group)
      end do
      e = net%transp(1) ; u = sum(net%uptake(:, 1)) ; sap = net%sapflow(1)
      call check_close('plant only: at steady state root uptake = transpiration', u, e, 1.0e-6_wp * e)
      call check_close('plant only: at steady state sapflow = transpiration', sap, e, 1.0e-6_wp * e)
      call check_close('plant only: the plant water stops changing', w_plant1 - w_plant0, 0.0_wp, 1.0e-8_wp * w_plant1)
      call check_true('plant only: the leaf sits below the wood by more than the crown head',                   &
                      net%psi_group(WATER_LEAF, 1) < net%psi_group(WATER_WOOD, 1) - net%grav(1),                &
                      net%psi_group(WATER_LEAF, 1))
   end subroutine test_plant_only

   !=======================================================================================!
   !----- No cohorts, a sealed bottom, no rain and no evaporation: the column relaxes to hydrostatic !
   !      equilibrium, the total head (matric head plus elevation) the same in every layer. ---------!
   subroutine test_soil_only()
      type(water_network_t)  :: net
      type(soil_opts_t)      :: sealed
      type(chydro_forcing_t) :: still
      real(wp) :: theta(soil%n_active), head(soil%n_active), w0
      integer(ik) :: step, nsl
      nsl = soil%n_active
      theta = 0.40_wp ; theta(1:3) = 0.425_wp          ! wet: BCI's clay settles in days, not years
      sealed = sopts ; sealed%bottom_bc = SOIL_BC_BEDROCK
      still = forc ; still%q_air = 0.05_wp                 ! air wetter than any surface: no evaporation
      call make_stand(net, 0_ik, theta, 0.0_wp, 0.9_wp, 0.9_wp, 0.0_wp)
      w0 = total_water(net)
      do step = 1_ik, 96_ik * 30_ik
         call water_network_advance(net, soil, sealed, still, hyd, 900.0_wp)
      end do
      head = net%psi_core(2:nsl+1) / grav_head + soil%z_node(1:nsl)
      call check_true('soil only: hydrostatic equilibrium (total head uniform to 1 cm)',                        &
                      maxval(head) - minval(head) < 1.0e-2_wp, maxval(head) - minval(head))
      call check_close('soil only: sealed column keeps its water', total_water(net), w0, 1.0e-11_wp * w0)
   end subroutine test_soil_only

   !=======================================================================================!
   subroutine test_stiff_exchange()
      type(water_network_t) :: net
      real(wp) :: theta(soil%n_active)
      integer(ik) :: step, nsub
      theta = 0.42_wp
      call make_stand(net, 4_ik, theta, 0.0_wp, 0.97_wp, 0.95_wp, 0.0_wp)
      net%g_root = 50.0_wp * net%g_root                ! a fifty-fold root density: a sub-second exchange
      nsub = 0_ik
      do step = 1_ik, 8_ik
         net%e_demand = merge(6.0e-5_wp, 0.0_wp, mod(step, 2_ik) == 0_ik) * net%nplant / sum(net%nplant)
         call water_network_advance(net, soil, sopts, forc, hyd, 900.0_wp)
         nsub = nsub + net%nsub + net%nreject
      end do
      call check_true('stiff exchange: demand jumps cost at most 6 sub-steps a step', nsub <= 48_ik, real(nsub, wp))
   end subroutine test_stiff_exchange

   !=======================================================================================!
   subroutine test_infeasible_demand()
      type(water_network_t) :: net
      real(wp) :: theta(soil%n_active), demand
      theta = 0.15_wp
      call make_stand(net, 1_ik, theta, 0.0_wp, 0.6_wp, 0.5_wp, 0.0_wp)
      net%bleaf(1) = 1.0e-3_wp                          ! a few grams of leaf
      net%k_xylem(1) = 1.0e-3_wp * net%k_xylem(1)        ! and a xylem that can hardly feed it
      call water_network_start(net, soil, hyd, theta, 0.0_wp, 0.0_wp,                                      &
                               [1.0e-3_wp * 2.0_wp * 0.75_wp], [net%w_group(WATER_WOOD, 1) / net%nplant(1)])
      demand = 1.0e-4_wp
      net%e_demand(1) = demand
      call water_network_advance(net, soil, sopts, forc, hyd, 900.0_wp)
      call check_int('infeasible demand: no sub-step forced through', net%nfail, 0_ik)
      call check_true('infeasible demand: the leaf stays above empty', net%tau_group(WATER_LEAF, 1) > 0.0_wp,     &
                      net%tau_group(WATER_LEAF, 1))
      call check_true('infeasible demand: transpiration falls short of the demand',                          &
                      net%transp(1) < 0.5_wp * demand * 900.0_wp, net%transp(1) / (demand * 900.0_wp))
      call check_true('infeasible demand: leaf potential is finite', abs(net%psi_group(WATER_LEAF, 1)) < 1.0e6_wp,  &
                      net%psi_group(WATER_LEAF, 1))
   end subroutine test_infeasible_demand

   !=======================================================================================!
   subroutine test_leafless()
      type(water_network_t) :: net
      real(wp) :: theta(soil%n_active)
      theta = 0.30_wp
      call make_stand(net, 3_ik, theta, 0.0_wp, 0.9_wp, 0.8_wp, 2.0e-5_wp, leafless=.true.)
      net%e_demand(3) = 0.0_wp
      call water_network_advance(net, soil, sopts, forc, hyd, 900.0_wp)
      call check_close('leafless: leaf at the wood less the crown head', net%psi_group(WATER_LEAF, 3),            &
                       net%psi_group(WATER_WOOD, 3) - net%grav(3), 1.0e-6_wp)
      call check_close('leafless: its empty leaf store stays empty', net%w_group(WATER_LEAF, 3), 0.0_wp, 1.0e-12_wp)
   end subroutine test_leafless

   !=======================================================================================!
   !----- The BCI regrowth crash (1994-04-20): the top 2.6 cm a hair above residual, a stand    !
   !      whose wood is half empty and refilling, then 30 mm/h of rain for two steps. ----------!
   subroutine test_dry_layer_storm()
      type(water_network_t) :: net
      real(wp) :: theta(soil%n_active), w0, flows, theta_min
      integer(ik) :: step, nsl
      logical :: finite
      nsl = soil%n_active
      theta = 0.16_wp ; theta(1) = 0.0981_wp + 2.0e-3_wp
      call make_stand(net, 6_ik, theta, 0.0_wp, 0.45_wp, 0.40_wp, 0.0_wp)
      w0 = total_water(net) ; flows = 0.0_wp
      finite = .true. ; theta_min = 1.0_wp
      do step = 1_ik, 6_ik
         net%rain = merge(30.0_wp / 3600.0_wp, 0.0_wp, step == 3 .or. step == 4)
         call water_network_advance(net, soil, sopts, forc, hyd, 900.0_wp)
         flows = flows + net%rain * 900.0_wp - net%drain - net%evap - net%spill - sum(net%transp)
         finite = finite .and. all(net%w_core == net%w_core) .and. all(net%psi_group == net%psi_group)
         theta_min = min(theta_min, minval(net%w_core(2:nsl+1) / (rho_h2o * soil%dz(1:nsl)) - soil%theta_res(1:nsl)))
         call check_int('dry-layer storm: no sub-step forced through', net%nfail, 0_ik)
      end do
      call check_true('dry-layer storm: every state finite', finite)
      call check_true('dry-layer storm: no layer below residual', theta_min > 0.0_wp, theta_min)
      call check_close('dry-layer storm: water conserved', total_water(net) - w0, flows, 1.0e-10_wp * w0)
      call check_true('dry-layer storm: the wood refilled', net%tau_group(WATER_WOOD, 1) > 0.40_wp,                &
                      net%tau_group(WATER_WOOD, 1))
   end subroutine test_dry_layer_storm

   !=======================================================================================!
   subroutine test_saturated()
      type(water_network_t) :: net
      real(wp) :: theta(soil%n_active), w0, flows, spill, over
      integer(ik) :: step, nsl
      nsl = soil%n_active
      theta = 0.4314_wp - 1.0e-3_wp
      call make_stand(net, 2_ik, theta, 0.0_wp, 0.97_wp, 0.95_wp, 1.0e-5_wp)
      w0 = total_water(net) ; flows = 0.0_wp ; spill = 0.0_wp ; over = 0.0_wp
      do step = 1_ik, 8_ik
         net%rain = merge(50.0_wp / 3600.0_wp, 0.0_wp, step <= 4)
         call water_network_advance(net, soil, sopts, forc, hyd, 900.0_wp)
         flows = flows + net%rain * 900.0_wp - net%drain - net%evap - net%spill - sum(net%transp)
         spill = spill + net%spill
         over = max(over, maxval(net%w_core(2:nsl+1) / (rho_h2o * soil%dz(1:nsl)) - soil%theta_sat(1:nsl)))
         call check_int('saturated: no sub-step forced through', net%nfail, 0_ik)
      end do
      call check_true('saturated: the pond spilled', spill > 0.0_wp, spill)
      call check_true('saturated: no layer more than 1e-3 above saturation', over < 1.0e-3_wp, over)
      call check_true('saturated: pond at or below its capacity', net%w_core(WATER_NODE_POND) <= sopts%w_pond_max,  &
                      net%w_core(WATER_NODE_POND))
      call check_close('saturated: water conserved', total_water(net) - w0, flows, 1.0e-10_wp * w0)
   end subroutine test_saturated

   !=======================================================================================!
   !----- A transpiring day with a midday storm at the default tolerances, against the same day at  !
   !      tolerances a thousand times tighter. ------------------------------------------------!
   subroutine test_accuracy()
      type(water_network_t) :: net, ref
      real(wp) :: theta(soil%n_active), hour, worst_theta, worst_rwc
      integer(ik) :: step, nsl, isub
      nsl = soil%n_active
      theta = 0.30_wp
      call make_stand(net, 3_ik, theta, 0.0_wp, 0.9_wp, 0.85_wp, 0.0_wp)
      call make_stand(ref, 3_ik, theta, 0.0_wp, 0.9_wp, 0.85_wp, 0.0_wp)
      ref%rtol = 1.0e-3_wp * net%rtol ; ref%atol = 1.0e-3_wp * net%atol
      isub = 0_ik
      do step = 1_ik, 96_ik
         hour = real(step, wp) * 0.25_wp
         net%e_demand = 5.0e-5_wp * max(0.0_wp, sin(pi * (hour - 6.0_wp) / 12.0_wp))                        &
                        * net%nplant * net%bleaf / sum(net%nplant * net%bleaf)
         ref%e_demand = net%e_demand
         net%rain = merge(1.0e-2_wp, 0.0_wp, step >= 50 .and. step < 54) ; ref%rain = net%rain
         call water_network_advance(net, soil, sopts, forc, hyd, 900.0_wp)
         call water_network_advance(ref, soil, sopts, forc, hyd, 900.0_wp)
         isub = isub + net%nsub
      end do
      worst_theta = maxval(abs(net%w_core(2:nsl+1) - ref%w_core(2:nsl+1)) / (rho_h2o * soil%dz(1:nsl)))
      worst_rwc   = maxval(abs(net%tau_group - ref%tau_group))
      call check_true('accuracy: soil water within 10 atol of the tight run', worst_theta < 10.0_wp * sopts%atol,  &
                      worst_theta)
      call check_true('accuracy: tissue water within 1e-2 RWC of the tight run', worst_rwc < 1.0e-2_wp, worst_rwc)
      print '(a,i6)', '  info : sub-steps over the day at the default tolerances: ', isub
   end subroutine test_accuracy

   !=======================================================================================!
   !----- Two stores the BCI regrowth stand (1994-04-20) showed: a cohort with 1e-22 kgC of leaves,  !
   !      whose leaf store holds 1e-26 kg/m2 against its plant's 1e-9 tolerance, and a pond holding a !
   !      trace of water over a dry top layer in light rain, which drains tens of millimetres an hour. !
   !      A day of them: no sub-step forced, the leaf never above saturation, the pond never below  !
   !      empty. --------------------------------------------------------------------------------!
   subroutine test_vanishing_stores()
      type(water_network_t) :: net
      real(wp) :: theta(soil%n_active), leaf_w(3), wood_w(3), w0, flows, rwc_max, pond_min
      integer(ik) :: step, nfail
      theta = 0.16_wp
      call make_stand(net, 3_ik, theta, 0.0_wp, 0.9_wp, 0.85_wp, 0.0_wp)
      net%bleaf(3) = 1.3e-22_wp
      leaf_w = net%w_group(WATER_LEAF, :) / net%nplant ; wood_w = net%w_group(WATER_WOOD, :) / net%nplant
      leaf_w(3) = 2.0_wp * net%bleaf(3) * 0.9_wp
      call water_network_start(net, soil, hyd, theta, 3.0e-8_wp, 0.0_wp, leaf_w, wood_w)
      w0 = total_water(net) ; flows = 0.0_wp ; nfail = 0_ik ; rwc_max = 0.0_wp ; pond_min = 0.0_wp
      do step = 1_ik, 96_ik
         net%rain = merge(1.0e-6_wp, 1.2e-3_wp, step < 40)
         net%e_demand(1:2) = 2.0e-5_wp ; net%e_demand(3) = 1.0e-32_wp
         call water_network_advance(net, soil, sopts, forc, hyd, 900.0_wp)
         flows = flows + net%rain * 900.0_wp - net%drain - net%evap - net%spill - sum(net%transp)
         nfail = nfail + net%nfail
         rwc_max  = max(rwc_max, net%tau_group(WATER_LEAF, 3))
         pond_min = min(pond_min, net%w_core(WATER_NODE_POND))
      end do
      call check_int('vanishing stores: no sub-step forced through', nfail, 0_ik)
      call check_true('vanishing stores: the near-leafless leaf stays at or below saturation', rwc_max <= 1.0_wp, rwc_max)
      call check_true('vanishing stores: the pond never holds less than nothing', pond_min >= 0.0_wp, pond_min)
      call check_close('vanishing stores: water conserved', total_water(net) - w0, flows, 1.0e-10_wp * w0)
   end subroutine test_vanishing_stores

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

end program test_water_network
