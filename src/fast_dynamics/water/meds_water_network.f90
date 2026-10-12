! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_water_network -- one patch's soil and plant water, solved together as one network     !
! (MEDS_SPAC_WATER_NETWORK_PLAN.md §4, §12).                                                  !
!                                                                                          !
! THE NETWORK. Water sits in stores joined by paths:                                         !
!   * the shared core: the pond (node 1) and the soil layers (node 1 + k);                   !
!   * one group per cohort: its wood and its leaves, which touch only each other and the     !
!     soil layers its roots reach;                                                           !
!   * paths: Darcy flow between adjacent layers and between the pond and the top layer; root  !
!     uptake from each layer into each cohort's wood; sapflow up each cohort's xylem;          !
!   * boundaries: rain into the pond, ground evaporation from the top layer, drainage from the  !
!     bottom one, and each cohort's transpiration, limited by what its leaves hold.             !
! What each store holds and how fast each path moves water are the process library's          !
! (meds_water_storage, meds_water_paths, and the soil's boundary processes in meds_soil_water);  !
! this module only joins them.                                                                    !
!                                                                                          !
! THE SOLVE. Over one fast step, implicit Euler on sub-steps. Each sub-step is a Newton solve    !
! of the water balance of every store at the sub-step's end,                                     !
!     (W(tau) - W_start) / h - (net inflow at tau) = 0,                                           !
! in each store's own unknown tau (meds_water_storage). Its linear systems are solved by          !
! eliminating each cohort onto the soil (meds_grouped_solve). Newton steps are bounded so no       !
! store is driven past empty, and backtracked until the balance improves. A sub-step's error is    !
! the implicit-Euler local error, estimated from the change in the inflows across the sub-step     !
! and filtered twice through the Newton matrix so a stiff exchange that has already settled does   !
! not read as error (§12.2 M1); it sets the next sub-step.                                          !
!                                                                                          !
! WHAT IS KEPT EXACT. Each store's water is advanced by the flows themselves, W += h * inflow,     !
! never re-read from its curve, so every kilogram that leaves one store arrives in another or      !
! crosses a boundary: the network conserves water to round-off whatever Newton's residual. The     !
! amounts each path moved over the step are reported for the energy integrator and the ledgers.    !
!==========================================================================================!
module meds_water_network
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : grav_head, rho_h2o, tiny_num
   use meds_water_storage,   only : soil_curve_t, soil_curve_setup, soil_store, soil_unknown_from_theta,   &
                                    tissue_store, tissue_rwc_from_water, tissue_unknown_from_rwc,          &
                                    pond_store, pond_unknown_from_water
   use meds_water_paths,     only : darcy_flow, root_uptake_profile, xylem_flow, supply_limited_transpiration
   use meds_water_retention, only : water_curve_t, rwc_from_psi
   use meds_grouped_solve,   only : grouped_system_t, grouped_allocate, grouped_factor, grouped_solve
   use meds_soil_water,      only : ground_evaporation_rate, bottom_flux
   use meds_soil_types,      only : chydro_forcing_t
   use meds_column_params,   only : soil_params_t, curve_a, curve_n
   use meds_biophysics_opts, only : soil_opts_t
   use meds_plant_types,     only : hydro_params_table_t
   use meds_therm_lib,       only : internal_energy_liquid, temp_of_liquid_enthalpy
   implicit none
   private

   public :: water_network_t, water_network_allocate, water_network_start, water_network_advance
   public :: water_network_jacobian_error
   public :: WATER_NODE_POND, WATER_WOOD, WATER_LEAF

   integer(ik), parameter :: WATER_NODE_POND = 1_ik   !< the pond is core node 1; soil layer k is 1 + k
   integer(ik), parameter :: WATER_WOOD = 1_ik        !< a cohort's wood, first in its group
   integer(ik), parameter :: WATER_LEAF = 2_ik        !< its leaves, second
   integer(ik), parameter :: GROUP_SIZE = 2_ik

   !----- Newton. The balance counts as solved when every store's imbalance over the sub-step is  !
   !      below NEWTON_TOL of the step's absolute water tolerance -- far inside the error the     !
   !      sub-step itself is allowed. A Newton step may take a soil layer or a tissue at most        !
   !      STEP_MAX_SHARE of the way to empty (fraction to the boundary) and move its saturation or   !
   !      relative water content by at most STEP_MAX_DRY. --------------------------------------!
   real(wp),    parameter :: NEWTON_TOL     = 1.0e-2_wp
   integer(ik), parameter :: NEWTON_MAX     = 25_ik
   integer(ik), parameter :: BACKTRACK_MAX  = 12_ik
   real(wp),    parameter :: STEP_MAX_SHARE = 0.9_wp
   real(wp),    parameter :: STEP_MAX_DRY   = 0.2_wp
   !----- Sub-steps: the error controller of an order-one method (exponent -1/2), its growth and   !
   !      shrink limits, the factor a sub-step whose Newton failed is cut by, and the shortest      !
   !      sub-step, as a share of the fast step, below which a sub-step is taken as it stands. ----!
   real(wp),    parameter :: SAFETY = 0.9_wp, GROW_MAX = 4.0_wp, SHRINK_MIN = 0.2_wp
   real(wp),    parameter :: NEWTON_FAIL_CUT = 0.25_wp
   real(wp),    parameter :: H_MIN_SHARE = 1.0e-5_wp
   integer(ik), parameter :: TRIALS_MAX = 2000_ik
   !----- [kg/m2] Below this a pond holds no water to carry a temperature (as meds_soil_water). ---!
   real(wp),    parameter :: POND_TINY = 1.0e-9_wp

   !----- Everything one evaluation of the network gives: each store's potential and water and  !
   !      their slopes, the net inflow into each store, and each path's flow. -------------------!
   type :: net_eval_t
      real(wp), allocatable :: psi_c(:), dpsi_c(:), w_c(:), dw_c(:), k_c(:), dk_c(:), f_c(:)
      real(wp), allocatable :: psi_g(:,:), dpsi_g(:,:), w_g(:,:), dw_g(:,:), f_g(:,:)
      real(wp), allocatable :: q_face(:), u_root(:,:), s_xyl(:), e_leaf(:)
      real(wp) :: q_pond = 0.0_wp, q_drain = 0.0_wp, e_soil = 0.0_wp
   end type net_eval_t

   !----- One patch's network: what it is built from, its state, and what it moved. ------------!
   type :: water_network_t
      integer(ik) :: nsl = 0_ik           !< soil layers
      integer(ik) :: n   = 0_ik           !< cohorts
      integer(ik) :: nc  = 0_ik           !< core nodes, 1 + nsl
      !----- The soil and the pond. -------------------------------------------------------!
      type(soil_curve_t), allocatable :: soil(:)       !< (nsl) each layer's curve
      real(wp), allocatable :: z_node(:)               !< (nsl) [m, < 0] node depth
      real(wp), allocatable :: dz_node(:)              !< (nsl) [m] distance to the node below
      real(wp), allocatable :: face_weight(:)          !< (nsl) lower layer's weight in the face conductivity
      real(wp), allocatable :: t_soil(:)               !< (nsl) [K] temperature of water a layer gives up
      real(wp) :: rain       = 0.0_wp                  !< [kg/m2/s] water reaching the pond
      real(wp) :: t_rain     = 298.15_wp               !< [K] the temperature that values it
      real(wp) :: w_pond_max = 5.0_wp                  !< [kg/m2] pond capacity; above it, runoff
      real(wp) :: rtol = 1.0e-3_wp, atol = 1.0e-4_wp   !< sub-step error tolerances (plan D8)
      !----- The cohorts. ----------------------------------------------------------------!
      integer(ik), allocatable :: pft(:)               !< (n) PFT, selecting curves and vulnerability
      real(wp), allocatable :: nplant(:)               !< (n) [plants/m2]
      real(wp), allocatable :: bleaf(:), bwood(:)      !< (n) [kgC/plant] leaf; sapwood + fine root
      real(wp), allocatable :: grav(:)                 !< (n) [MPa] hydrostatic head from wood to leaves
      real(wp), allocatable :: k_xylem(:)              !< (n) [kg/s/MPa per m2] intact xylem conductance
      real(wp), allocatable :: e_demand(:)             !< (n) [kg/m2/s] the stomata's transpiration demand
      real(wp), allocatable :: r_close(:)              !< (n) [-] leaf RWC at which the stomata are shut
      real(wp), allocatable :: g_root(:,:)             !< (nsl, n) [kg/s/MPa per m/s per m2] root geometry
      !----- Core nodes held at their potential: a plant-only network prescribes the soil. -----!
      logical,  allocatable :: core_fixed(:)           !< (nc)
      !----- State: each store's unknown and the water it holds [kg/m2]. -----------------------!
      real(wp), allocatable :: tau_core(:), tau_group(:,:)   !< (nc), (2, n)
      real(wp), allocatable :: w_core(:),   w_group(:,:)     !< (nc), (2, n)
      real(wp) :: e_pond = 0.0_wp                            !< [J/m2] the pond's enthalpy
      !----- What moved over the step [kg/m2, J/m2]. --------------------------------------!
      real(wp), allocatable :: face(:)                 !< (nsl) down through the face below layer k
      real(wp), allocatable :: uptake(:,:)             !< (nsl, n) from each layer into each cohort
      real(wp), allocatable :: sapflow(:), transp(:)   !< (n)
      real(wp) :: infil = 0.0_wp, infil_enth = 0.0_wp  !< pond -> top layer, and its enthalpy
      real(wp) :: exfil = 0.0_wp, exfil_enth = 0.0_wp  !< top layer -> pond
      real(wp) :: drain = 0.0_wp, evap = 0.0_wp        !< out of the bottom; out of the top layer to the air
      real(wp) :: spill = 0.0_wp, spill_enth = 0.0_wp  !< pond runoff over its capacity
      real(wp) :: rain_enth = 0.0_wp                   !< enthalpy the rain brought
      !----- End-of-step potentials [MPa]. ----------------------------------------------!
      real(wp), allocatable :: psi_core(:), psi_group(:,:)
      !----- Work done: accepted sub-steps, Newton iterations, rejected sub-steps, and sub-steps  !
      !      taken although they did not converge or meet the tolerance (counted, not hidden). ---!
      integer(ik) :: nsub = 0_ik, nnewton = 0_ik, nreject = 0_ik, nfail = 0_ik
      !----- Scratch, kept from step to step so the arrays are allocated again only when the shape  !
      !      changes. ev(i_start) is the network at a sub-step's start; Newton evaluates its iterates  !
      !      into the other. sys holds the flows' Jacobian of the latest evaluation (its factorisation !
      !      is kept apart). A sub-step starts where the last accepted one ended, so while            !
      !      start_ready holds, its start's evaluation and Jacobian are already there. ---------------!
      type(grouped_system_t) :: sys
      type(net_eval_t)       :: ev(2)
      integer(ik)            :: i_start = 1_ik
      logical                :: start_ready = .false.
   end type water_network_t

contains

   !----- Size a network for nsl soil layers and n cohorts, keeping the storage when the shape is   !
   !      unchanged, and reset what the caller fills. --------------------------------------------!
   subroutine water_network_allocate(net, nsl, n)
      type(water_network_t), intent(inout) :: net
      integer(ik),           intent(in)    :: nsl, n
      if (.not. allocated(net%soil) .or. net%nsl /= nsl .or. net%n /= n) then
         net%nsl = nsl ; net%n = n ; net%nc = 1_ik + nsl
         if (allocated(net%soil)) deallocate(net%soil, net%z_node, net%dz_node, net%face_weight, net%t_soil,  &
              net%pft, net%nplant, net%bleaf, net%bwood, net%grav, net%k_xylem, net%e_demand, net%r_close,     &
              net%g_root, net%core_fixed, net%tau_core, net%tau_group, net%w_core, net%w_group, net%face,       &
              net%uptake, net%sapflow, net%transp, net%psi_core, net%psi_group)
         allocate(net%soil(nsl), net%z_node(nsl), net%dz_node(nsl), net%face_weight(nsl), net%t_soil(nsl))
         allocate(net%pft(n), net%nplant(n), net%bleaf(n), net%bwood(n), net%grav(n), net%k_xylem(n),        &
                  net%e_demand(n), net%r_close(n), net%g_root(nsl, n))
         allocate(net%core_fixed(net%nc), net%tau_core(net%nc), net%tau_group(GROUP_SIZE, n),               &
                  net%w_core(net%nc), net%w_group(GROUP_SIZE, n))
         allocate(net%face(nsl), net%uptake(nsl, n), net%sapflow(n), net%transp(n),                         &
                  net%psi_core(net%nc), net%psi_group(GROUP_SIZE, n))
         call eval_allocate(net%ev(1), net) ; call eval_allocate(net%ev(2), net)
      end if
      net%t_soil = 298.15_wp ; net%core_fixed = .false. ; net%g_root = 0.0_wp
      net%e_demand = 0.0_wp ; net%r_close = 0.5_wp
   end subroutine water_network_allocate

   !----- Set the network's state from the model's: soil water content, the pond's water and       !
   !      enthalpy, and each cohort's leaf and wood water [kg per plant]. The cohort inputs          !
   !      (pft ... g_root) must already be filled. Each store's unknown is read from its water; a    !
   !      leafless cohort's leaves hold none, so their unknown starts at the wood's potential less   !
   !      the head to the crown, where an empty leaf comes to rest. ----------------------------!
   subroutine water_network_start(net, params, hyd, theta, w_pond, e_pond, leaf_water, wood_water)
      type(water_network_t),      intent(inout) :: net
      type(soil_params_t),        intent(in)    :: params
      type(hydro_params_table_t), intent(in)    :: hyd
      real(wp),                   intent(in)    :: theta(:)         !< (nsl) [m3/m3]
      real(wp),                   intent(in)    :: w_pond, e_pond   !< [kg/m2], [J/m2]
      real(wp),                   intent(in)    :: leaf_water(:), wood_water(:)   !< (n) [kg/plant]
      integer(ik) :: k, i
      real(wp)    :: psi, dpsi, w, dw, kk, dk, head1
      do k = 1_ik, net%nsl
         call soil_curve_setup(params%retention, params%theta_sat(k), params%theta_res(k),                   &
                               curve_a(params, k), curve_n(params, k), params%ksat(k), params%dz(k), net%soil(k))
         net%z_node(k) = params%z_node(k)
         if (k < net%nsl) then
            net%dz_node(k)     = params%dz_node(k)
            net%face_weight(k) = params%dz(k) / (params%dz(k) + params%dz(k+1))
         end if
         net%tau_core(1_ik + k) = soil_unknown_from_theta(net%soil(k), theta(k))
         net%w_core(1_ik + k)   = rho_h2o * params%dz(k) * theta(k)
      end do
      !----- The pond starts at its water's head, but no higher than the head at which it would     !
      !      neither give nor take water from the top layer. Over dry soil the pond face drains tens of !
      !      millimetres an hour, so a pond holding a trace of water empties within the step: starting !
      !      it at that trace's head would put the step's first Newton iterate at a flow thousands of   !
      !      times what the pond holds, far outside the linear model's reach. ----------------------!
      net%w_core(WATER_NODE_POND) = w_pond
      net%e_pond = e_pond
      call soil_store(net%soil(1), net%tau_core(2), psi, dpsi, w, dw, kk, dk)
      head1 = psi / grav_head
      net%tau_core(WATER_NODE_POND) = min(pond_unknown_from_water(w_pond), head1 + net%z_node(1))
      do i = 1_ik, net%n
         associate (hp => hyd%pft(net%pft(i)))
            net%tau_group(WATER_WOOD, i) = tissue_unknown_from_rwc(                                          &
                                              tissue_rwc_from_water(hp%wood_curve, net%bwood(i), wood_water(i)))
            net%w_group(WATER_WOOD, i)   = wood_water(i) * net%nplant(i)
            if (hp%leaf_curve%water_sat * net%bleaf(i) > tiny_num) then
               net%tau_group(WATER_LEAF, i) = tissue_unknown_from_rwc(                                       &
                                                 tissue_rwc_from_water(hp%leaf_curve, net%bleaf(i), leaf_water(i)))
               net%w_group(WATER_LEAF, i) = leaf_water(i) * net%nplant(i)
            else
               call tissue_store(hp%wood_curve, net%bwood(i), net%nplant(i), net%tau_group(WATER_WOOD, i),     &
                                 psi, dpsi, w, dw)
               net%tau_group(WATER_LEAF, i) = tissue_unknown_from_rwc(rwc_from_psi(psi - net%grav(i),           &
                                                 hp%leaf_curve%pi0, hp%leaf_curve%elastic_mod))
               net%w_group(WATER_LEAF, i) = 0.0_wp
            end if
         end associate
      end do
   end subroutine water_network_start

   !----- Advance the network over dt. On return the state holds the end-of-step water and          !
   !      unknowns, and the amounts each path moved. -----------------------------------------!
   subroutine water_network_advance(net, params, opts, forcing, hyd, dt)
      type(water_network_t),      intent(inout) :: net
      type(soil_params_t),        intent(in)    :: params
      type(soil_opts_t),          intent(in)    :: opts
      type(chydro_forcing_t),     intent(in)    :: forcing
      type(hydro_params_table_t), intent(in)    :: hyd
      real(wp),                   intent(in)    :: dt

      real(wp) :: t, h, h_min, err, scale_c(net%nc), scale_g(GROUP_SIZE, net%n)
      real(wp) :: tau_c(net%nc), tau_g(GROUP_SIZE, net%n)
      integer(ik) :: trials, iters, i_sol
      logical  :: converged, last, spilled

      call clear_amounts(net)
      call grouped_allocate(net%sys, net%nc, net%n, GROUP_SIZE)
      call water_scales(net, hyd, scale_c, scale_g)
      net%start_ready = .false.

      t = 0.0_wp ; h = dt ; h_min = H_MIN_SHARE * dt ; trials = 0_ik
      do while (t < dt * (1.0_wp - 1.0e-12_wp))
         h = min(h, dt - t)
         last = (h <= h_min) .or. (trials >= TRIALS_MAX)
         trials = trials + 1_ik
         tau_c = net%tau_core ; tau_g = net%tau_group
         call newton_substep(net, params, opts, forcing, hyd, h, scale_c, scale_g, tau_c, tau_g, converged, iters,  &
                             i_sol)
         net%nnewton = net%nnewton + iters
         if (.not. converged .and. .not. last) then
            net%nreject = net%nreject + 1_ik
            h = max(h * NEWTON_FAIL_CUT, h_min)
            cycle
         end if
         if (converged) then
            err = substep_error(net, h, net%ev(net%i_start), net%ev(i_sol), scale_c, scale_g)
         else
            err = huge(1.0_wp)
         end if
         if (err > 1.0_wp .and. .not. last) then
            net%nreject = net%nreject + 1_ik
            h = max(h * max(SHRINK_MIN, SAFETY * err ** (-0.5_wp)), h_min)
            cycle
         end if
         if (err > 1.0_wp) net%nfail = net%nfail + 1_ik
         call commit_substep(net, h, tau_c, tau_g, net%ev(i_sol), spilled)
         !----- The next sub-step starts where this one ended: its start is this one's solution, whose  !
         !      Jacobian sys still holds -- unless the spill moved the pond. ----------------------------!
         net%i_start = i_sol
         net%start_ready = .not. spilled
         t = t + h
         net%nsub = net%nsub + 1_ik
         if (err <= 1.0_wp)                                                                                   &
            h = h * min(GROW_MAX, max(SHRINK_MIN, SAFETY * max(err, 1.0e-10_wp) ** (-0.5_wp))) ! clamp-ok: controller floor
      end do
      net%psi_core  = net%ev(net%i_start)%psi_c
      net%psi_group = net%ev(net%i_start)%psi_g
   end subroutine water_network_advance

   !----- How far the assembled Jacobian is from the residual's own derivative, at the network's   !
   !      current state for a sub-step of h: the largest difference between J v and a central       !
   !      difference of the residual along v, relative to the largest entry of either, with v a     !
   !      fixed pattern of small changes in every unknown. For the test suite: it checks the whole  !
   !      network's Jacobian against the flows it is built from. ---------------------------------!
   real(wp) function water_network_jacobian_error(net, params, opts, forcing, hyd, h) result(err)
      type(water_network_t),      intent(inout) :: net
      type(soil_params_t),        intent(in)    :: params
      type(soil_opts_t),          intent(in)    :: opts
      type(chydro_forcing_t),     intent(in)    :: forcing
      type(hydro_params_table_t), intent(in)    :: hyd
      real(wp),                   intent(in)    :: h
      type(net_eval_t) :: ev
      real(wp) :: v_c(net%nc), v_g(GROUP_SIZE, net%n), jv_c(net%nc), jv_g(GROUP_SIZE, net%n)
      real(wp) :: rp_c(net%nc), rp_g(GROUP_SIZE, net%n), rm_c(net%nc), rm_g(GROUP_SIZE, net%n)
      real(wp) :: eps, big, add_c(net%nc), add_g(GROUP_SIZE, net%n)
      integer(ik) :: c, i, j, l
      call eval_allocate(ev, net)
      call grouped_allocate(net%sys, net%nc, net%n, GROUP_SIZE)
      do c = 1_ik, net%nc
         v_c(c) = 1.0e-2_wp * (1.0_wp + 0.1_wp * real(c, wp)) * max(abs(net%tau_core(c)), 1.0e-3_wp) ! clamp-ok: test step scale
      end do
      do i = 1_ik, net%n
         do j = 1_ik, GROUP_SIZE
            v_g(j, i) = -1.0e-2_wp * (1.0_wp + 0.07_wp * real(i + j, wp)) * net%tau_group(j, i)
         end do
      end do
      call evaluate(net, params, opts, forcing, hyd, net%tau_core, net%tau_group, ev, .true.)
      call storage_diagonal(net, h, ev, add_c, add_g)
      jv_c = matmul(net%sys%a, v_c) + add_c * v_c ; jv_g = 0.0_wp
      do i = 1_ik, net%n
         jv_c = jv_c + net%sys%b(:, i) * v_g(WATER_WOOD, i)
         do l = 1_ik, GROUP_SIZE
            jv_g(l, i) = sum(net%sys%d(l, :, i) * v_g(:, i)) + add_g(l, i) * v_g(l, i)
         end do
         jv_g(WATER_WOOD, i) = jv_g(WATER_WOOD, i) + sum(net%sys%c(:, i) * v_c)
      end do
      eps = 1.0e-5_wp
      call evaluate(net, params, opts, forcing, hyd, net%tau_core + eps * v_c, net%tau_group + eps * v_g, ev, .false.)
      call residual(net, h, ev, rp_c, rp_g)
      call evaluate(net, params, opts, forcing, hyd, net%tau_core - eps * v_c, net%tau_group - eps * v_g, ev, .false.)
      call residual(net, h, ev, rm_c, rm_g)
      rp_c = (rp_c - rm_c) / (2.0_wp * eps) ; rp_g = (rp_g - rm_g) / (2.0_wp * eps)
      where (net%core_fixed) jv_c = 0.0_wp
      big = max(maxval(abs(jv_c)), maxval(abs(rp_c)), tiny_num)
      err = maxval(abs(jv_c - rp_c)) / big
      if (net%n > 0_ik) then
         big = max(maxval(abs(jv_g)), maxval(abs(rp_g)), tiny_num)
         err = max(err, maxval(abs(jv_g - rp_g)) / big)
      end if
   end function water_network_jacobian_error

   !=======================================================================================!
   !     One sub-step: Newton on the end-of-sub-step water balance.                          !
   !=======================================================================================!
   subroutine newton_substep(net, params, opts, forcing, hyd, h, scale_c, scale_g, tau_c, tau_g, converged, iters, &
                             i_sol)
      type(water_network_t),      intent(inout) :: net
      type(soil_params_t),        intent(in)    :: params
      type(soil_opts_t),          intent(in)    :: opts
      type(chydro_forcing_t),     intent(in)    :: forcing
      type(hydro_params_table_t), intent(in)    :: hyd
      real(wp),                   intent(in)    :: h, scale_c(:), scale_g(:,:)
      real(wp),                   intent(inout) :: tau_c(:), tau_g(:,:)   !< start of sub-step in, solution out
      logical,                    intent(out)   :: converged
      integer(ik),                intent(out)   :: iters
      integer(ik),                intent(out)   :: i_sol   !< which of net%ev holds the solution
      real(wp) :: r_c(net%nc), r_g(GROUP_SIZE, net%n), d_c(net%nc), d_g(GROUP_SIZE, net%n)
      real(wp) :: try_c(net%nc), try_g(GROUP_SIZE, net%n), r2_c(net%nc), r2_g(GROUP_SIZE, net%n)
      real(wp) :: merit, merit_try, worst, lam, lam_max
      integer(ik) :: it, ib, i_try
      logical  :: ok

      converged = .false. ; iters = 0_ik
      !----- The network at the start, evaluated unless it already is; the iterates are evaluated into  !
      !      the other slot, and the solution is the start until one is taken. ------------------------!
      i_sol = net%i_start ; i_try = 3_ik - net%i_start
      if (.not. net%start_ready) then
         call evaluate(net, params, opts, forcing, hyd, tau_c, tau_g, net%ev(i_sol), .true.)
         net%start_ready = .true.
      end if
      call residual(net, h, net%ev(i_sol), r_c, r_g)
      merit = residual_merit(net, h, r_c, r_g, scale_c, scale_g)
      worst = residual_worst(net, h, r_c, r_g, scale_c, scale_g)
      if (merit /= merit) return
      do it = 1_ik, NEWTON_MAX
         if (worst <= 1.0_wp) exit
         iters = it
         call newton_factor(net, h, net%ev(i_sol), ok)
         if (.not. ok) return
         call grouped_solve(net%sys, -r_c, -r_g, d_c, d_g)
         lam_max = step_limit(net, tau_c, tau_g, d_c, d_g)
         lam = lam_max
         !----- A trial's evaluation replaces the Jacobian sys held, so the start's is gone. ---------!
         net%start_ready = .false.
         do ib = 1_ik, BACKTRACK_MAX
            try_c = tau_c + lam * d_c ; try_g = tau_g + lam * d_g
            call evaluate(net, params, opts, forcing, hyd, try_c, try_g, net%ev(i_try), .true.)
            call residual(net, h, net%ev(i_try), r2_c, r2_g)
            merit_try = residual_merit(net, h, r2_c, r2_g, scale_c, scale_g)
            if (merit_try == merit_try .and. merit_try <= (1.0_wp - 1.0e-4_wp * lam) * merit) exit
            lam = 0.5_wp * lam
         end do
         if (merit_try /= merit_try) return
         tau_c = try_c ; tau_g = try_g ; r_c = r2_c ; r_g = r2_g ; merit = merit_try ; i_sol = i_try
         worst = residual_worst(net, h, r_c, r_g, scale_c, scale_g)
      end do
      !----- Factor the Newton matrix at the solution: the error estimate solves with it. -----------!
      call newton_factor(net, h, net%ev(i_sol), ok)
      converged = ok .and. worst <= 1.0_wp
   end subroutine newton_substep

   !----- The water-balance residual of every store [kg/m2/s]: (W(tau) - W_start)/h - inflow.  !
   !      A prescribed core node has none. -------------------------------------------------!
   subroutine residual(net, h, ev, r_c, r_g)
      type(water_network_t), intent(in)  :: net
      real(wp),              intent(in)  :: h
      type(net_eval_t),      intent(in)  :: ev
      real(wp),              intent(out) :: r_c(:), r_g(:,:)
      r_c = (ev%w_c - net%w_core) / h - ev%f_c
      r_g = (ev%w_g - net%w_group) / h - ev%f_g
      where (net%core_fixed) r_c = 0.0_wp
   end subroutine residual

   !----- Each store's imbalance over the sub-step against NEWTON_TOL of its absolute water        !
   !      tolerance. residual_worst, below one, means the balance is solved; residual_merit, the sum of !
   !      squares, is what the line search asks to fall -- a Newton step always points down it,      !
   !      which is not true of the worst store alone. ---------------------------------------------!
   pure real(wp) function residual_worst(net, h, r_c, r_g, scale_c, scale_g) result(m)
      type(water_network_t), intent(in) :: net
      real(wp),              intent(in) :: h, r_c(:), r_g(:,:), scale_c(:), scale_g(:,:)
      real(wp) :: tol
      tol = NEWTON_TOL * net%atol
      m = max(maxval(abs(r_c) * h / (tol * scale_c)), 0.0_wp)
      if (net%n > 0_ik) m = max(m, maxval(abs(r_g) * h / (tol * scale_g)))
   end function residual_worst

   pure real(wp) function residual_merit(net, h, r_c, r_g, scale_c, scale_g) result(m)
      type(water_network_t), intent(in) :: net
      real(wp),              intent(in) :: h, r_c(:), r_g(:,:), scale_c(:), scale_g(:,:)
      real(wp) :: tol
      tol = NEWTON_TOL * net%atol
      m = sum((r_c * h / (tol * scale_c)) ** 2)
      if (net%n > 0_ik) m = m + sum((r_g * h / (tol * scale_g)) ** 2)
   end function residual_merit

   !----- The largest share of a Newton step every store can take: none is driven more than        !
   !      STEP_MAX_SHARE of the way to empty, and none moves its saturation or relative water content  !
   !      by more than STEP_MAX_DRY on the dry side. The pond and a wet soil layer, whose unknowns are  !
   !      potentials, are not bounded. ---------------------------------------------------------!
   pure real(wp) function step_limit(net, tau_c, tau_g, d_c, d_g) result(lam)
      type(water_network_t), intent(in) :: net
      real(wp),              intent(in) :: tau_c(:), tau_g(:,:), d_c(:), d_g(:,:)
      integer(ik) :: k, i, j
      lam = 1.0_wp
      do k = 1_ik, net%nsl
         if (tau_c(1_ik + k) <= net%soil(k)%se_switch) call bound(tau_c(1_ik + k), d_c(1_ik + k), lam)
      end do
      do i = 1_ik, net%n
         do j = 1_ik, GROUP_SIZE
            call bound(tau_g(j, i), d_g(j, i), lam)
         end do
      end do
   contains
      !----- The cap on the change applies below full saturation (an RWC or Se of one), where the     !
      !      potential steepens towards empty; above it the curve is gentle and Newton moves freely. ---!
      pure subroutine bound(tau, d, lam)
         real(wp), intent(in)    :: tau, d
         real(wp), intent(inout) :: lam
         if (d < 0.0_wp) lam = min(lam, STEP_MAX_SHARE * tau / (-d))
         if (tau < 1.0_wp .and. abs(d) > STEP_MAX_DRY) lam = min(lam, STEP_MAX_DRY / abs(d))
      end subroutine bound
   end function step_limit

   !=======================================================================================!
   !     The network at one state: storage, flows, and (optionally) the Jacobian.            !
   !=======================================================================================!
   subroutine evaluate(net, params, opts, forcing, hyd, tau_c, tau_g, ev, with_jacobian)
      type(water_network_t),      intent(inout) :: net
      type(soil_params_t),        intent(in)    :: params
      type(soil_opts_t),          intent(in)    :: opts
      type(chydro_forcing_t),     intent(in)    :: forcing
      type(hydro_params_table_t), intent(in)    :: hyd
      real(wp),                   intent(in)    :: tau_c(:), tau_g(:,:)
      type(net_eval_t),           intent(inout) :: ev
      logical,                    intent(in)    :: with_jacobian
      integer(ik) :: k, i, ca, cb, c
      real(wp) :: q, dq_pa, dq_pb, dq_ka, dq_kb, dqa, dqb, u, s, ds_w, ds_l
      real(wp) :: e, de_dtheta, de_dpsi, de, qbot, dqbot_dpsi, dqbot_dk, theta1, head
      real(wp) :: psi_w, dpsi_w, u_sum, dqb_sum, du_ps_k(net%nsl), du_k_k(net%nsl)

      ev%f_c = 0.0_wp ; ev%f_g = 0.0_wp
      !----- The Jacobian's entries are written afresh: the core block's band, which the flows add  !
      !      into, is cleared; every other entry an evaluation touches is set outright below, and the  !
      !      rest stay zero from grouped_allocate. ---------------------------------------------------!
      if (with_jacobian) then
         do c = 1_ik, net%nc
            net%sys%a(c, c) = 0.0_wp
            if (c < net%nc) then
               net%sys%a(c, c + 1_ik) = 0.0_wp ; net%sys%a(c + 1_ik, c) = 0.0_wp
            end if
         end do
      end if

      !----- Storage. ------------------------------------------------------------------!
      call pond_store(tau_c(WATER_NODE_POND), ev%psi_c(1), ev%dpsi_c(1), ev%w_c(1), ev%dw_c(1))
      ev%k_c(1) = params%ksat(1) ; ev%dk_c(1) = 0.0_wp      ! the pond conducts as saturated soil
      do k = 1_ik, net%nsl
         c = 1_ik + k
         call soil_store(net%soil(k), tau_c(c), ev%psi_c(c), ev%dpsi_c(c), ev%w_c(c), ev%dw_c(c),          &
                         ev%k_c(c), ev%dk_c(c))
      end do
      do i = 1_ik, net%n
         associate (hp => hyd%pft(net%pft(i)))
            call tissue_store(hp%wood_curve, net%bwood(i), net%nplant(i), tau_g(WATER_WOOD, i),               &
                              ev%psi_g(WATER_WOOD, i), ev%dpsi_g(WATER_WOOD, i), ev%w_g(WATER_WOOD, i),       &
                              ev%dw_g(WATER_WOOD, i))
            call tissue_store(hp%leaf_curve, net%bleaf(i), net%nplant(i), tau_g(WATER_LEAF, i),               &
                              ev%psi_g(WATER_LEAF, i), ev%dpsi_g(WATER_LEAF, i), ev%w_g(WATER_LEAF, i),       &
                              ev%dw_g(WATER_LEAF, i))
         end associate
      end do

      !----- The pond and the top layer: Darcy, the pond a saturated layer as thick as the top one. !
      ca = WATER_NODE_POND ; cb = 2_ik
      call darcy_flow(ev%k_c(ca), ev%k_c(cb), 0.5_wp, ev%psi_c(ca), ev%psi_c(cb), -net%z_node(1),          &
                      -net%z_node(1), q, dq_pa, dq_pb, dq_ka, dq_kb)
      ev%q_pond = q
      call core_edge(ca, cb, q, dq_pa * ev%dpsi_c(ca), dq_pb * ev%dpsi_c(cb) + dq_kb * ev%dk_c(cb))

      !----- Between soil layers. --------------------------------------------------------!
      do k = 1_ik, net%nsl - 1_ik
         ca = 1_ik + k ; cb = 2_ik + k
         call darcy_flow(ev%k_c(ca), ev%k_c(cb), net%face_weight(k), ev%psi_c(ca), ev%psi_c(cb),            &
                         net%dz_node(k), net%dz_node(k), q, dq_pa, dq_pb, dq_ka, dq_kb)
         ev%q_face(k) = q
         dqa = dq_pa * ev%dpsi_c(ca) + dq_ka * ev%dk_c(ca)
         dqb = dq_pb * ev%dpsi_c(cb) + dq_kb * ev%dk_c(cb)
         call core_edge(ca, cb, q, dqa, dqb)
      end do

      !----- Rain into the pond; evaporation from the top layer; drainage from the bottom one. ----!
      ev%f_c(WATER_NODE_POND) = ev%f_c(WATER_NODE_POND) + net%rain
      c = 2_ik
      theta1 = ev%w_c(c) / (rho_h2o * net%soil(1)%dz)
      call ground_evaporation_rate(theta1, ev%psi_c(c) / grav_head, params, forcing, opts, e, de_dtheta, de_dpsi)
      ev%e_soil = e
      ev%f_c(c) = ev%f_c(c) - e
      if (with_jacobian) net%sys%a(c, c) = net%sys%a(c, c) + de_dtheta * ev%dw_c(c) / (rho_h2o * net%soil(1)%dz) &
                                                         + de_dpsi * ev%dpsi_c(c) / grav_head
      c = 1_ik + net%nsl
      head = ev%psi_c(c) / grav_head
      call bottom_flux(params, opts, net%nsl, head, ev%k_c(c), qbot, dqbot_dpsi, dqbot_dk)
      ev%q_drain = rho_h2o * qbot
      ev%f_c(c) = ev%f_c(c) - ev%q_drain
      if (with_jacobian) net%sys%a(c, c) = net%sys%a(c, c)                                                   &
                         + rho_h2o * (dqbot_dpsi * ev%dpsi_c(c) / grav_head + dqbot_dk * ev%dk_c(c))

      !----- Each cohort: uptake from its layers, sapflow to its leaves, and transpiration. The      !
      !      cohort's uptake comes from its whole root profile at once; its own sums are kept in       !
      !      scalars over its layers and written once. ------------------------------------------!
      do i = 1_ik, net%n
         psi_w = ev%psi_g(WATER_WOOD, i) ; dpsi_w = ev%dpsi_g(WATER_WOOD, i)
         call root_uptake_profile(net%nsl, net%g_root(1:net%nsl, i), ev%k_c(2:net%nc), ev%psi_c(2:net%nc),     &
                                  net%z_node(1:net%nsl), psi_w, ev%u_root(1:net%nsl, i), du_ps_k, du_k_k)
         u_sum = 0.0_wp ; dqb_sum = 0.0_wp
         do k = 1_ik, net%nsl
            c = 1_ik + k
            if (net%g_root(k, i) <= 0.0_wp) then
               if (with_jacobian) then
                  net%sys%b(c, i) = 0.0_wp ; net%sys%c(c, i) = 0.0_wp
               end if
               cycle
            end if
            u = ev%u_root(k, i)
            ev%f_c(c) = ev%f_c(c) - u
            u_sum = u_sum + u
            if (with_jacobian) then
               dqa = du_ps_k(k) * ev%dpsi_c(c) + du_k_k(k) * ev%dk_c(c)
               dqb = -du_ps_k(k) * dpsi_w
               net%sys%a(c, c) = net%sys%a(c, c) + dqa
               net%sys%b(c, i) = dqb
               net%sys%c(c, i) = -dqa
               dqb_sum = dqb_sum - dqb
            end if
         end do
         ev%f_g(WATER_WOOD, i) = ev%f_g(WATER_WOOD, i) + u_sum
         associate (hp => hyd%pft(net%pft(i)))
            call xylem_flow(net%k_xylem(i), ev%psi_g(WATER_WOOD, i), ev%psi_g(WATER_LEAF, i), net%grav(i),        &
                            hp%wood_psi50, hp%wood_kexp, hp%vuln_table, s, ds_w, ds_l)
         end associate
         ev%s_xyl(i) = s
         ev%f_g(WATER_WOOD, i) = ev%f_g(WATER_WOOD, i) - s
         ev%f_g(WATER_LEAF, i) = ev%f_g(WATER_LEAF, i) + s
         call supply_limited_transpiration(net%e_demand(i), tau_g(WATER_LEAF, i), net%r_close(i), e, de)
         ev%e_leaf(i) = e
         ev%f_g(WATER_LEAF, i) = ev%f_g(WATER_LEAF, i) - e
         if (with_jacobian) then
            dqa = ds_w * ev%dpsi_g(WATER_WOOD, i)
            dqb = ds_l * ev%dpsi_g(WATER_LEAF, i)
            net%sys%d(WATER_WOOD, WATER_WOOD, i) = dqb_sum + dqa
            net%sys%d(WATER_WOOD, WATER_LEAF, i) = dqb
            net%sys%d(WATER_LEAF, WATER_WOOD, i) = -dqa
            net%sys%d(WATER_LEAF, WATER_LEAF, i) = -dqb + de
         end if
      end do

   contains

      !----- A flow q from core node a to core node b, with dq/dtau_a and dq/dtau_b: the inflows  !
      !      and, when wanted, the Jacobian of the residual (W - W0)/h - inflow. ----------------!
      subroutine core_edge(a, b, q, dqa, dqb)
         integer(ik), intent(in) :: a, b
         real(wp),    intent(in) :: q, dqa, dqb
         ev%f_c(a) = ev%f_c(a) - q
         ev%f_c(b) = ev%f_c(b) + q
         if (.not. with_jacobian) return
         net%sys%a(a, a) = net%sys%a(a, a) + dqa
         net%sys%a(a, b) = net%sys%a(a, b) + dqb
         net%sys%a(b, a) = net%sys%a(b, a) - dqa
         net%sys%a(b, b) = net%sys%a(b, b) - dqb
      end subroutine core_edge

   end subroutine evaluate

   !----- The storage term dW/dtau / h each store adds to the Newton matrix's diagonal, and none    !
   !      for a prescribed core node. ----------------------------------------------------------!
   subroutine storage_diagonal(net, h, ev, add_c, add_g)
      type(water_network_t), intent(in)  :: net
      real(wp),              intent(in)  :: h
      type(net_eval_t),      intent(in)  :: ev
      real(wp),              intent(out) :: add_c(:), add_g(:,:)
      add_c = ev%dw_c / h
      add_g = ev%dw_g / h
      where (net%core_fixed) add_c = 0.0_wp
   end subroutine storage_diagonal

   !----- Factor the Newton matrix: the flows' Jacobian of the latest evaluation (sys) plus the     !
   !      storage terms of ev. A prescribed core node's row is the identity, so its unknown stays     !
   !      put; making it so edits the assembled rows, which the next evaluation writes afresh. ----!
   subroutine newton_factor(net, h, ev, ok)
      type(water_network_t), intent(inout) :: net
      real(wp),              intent(in)    :: h
      type(net_eval_t),      intent(in)    :: ev
      logical,               intent(out)   :: ok
      real(wp)    :: add_c(net%nc), add_g(GROUP_SIZE, net%n)
      integer(ik) :: c
      do c = 1_ik, net%nc
         if (net%core_fixed(c)) then
            net%sys%a(c, :) = 0.0_wp ; net%sys%a(c, c) = 1.0_wp
            net%sys%b(c, :) = 0.0_wp
         end if
      end do
      call storage_diagonal(net, h, ev, add_c, add_g)
      call grouped_factor(net%sys, add_c, add_g, ok)
   end subroutine newton_factor

   !=======================================================================================!
   !     Error, commit, and bookkeeping.                                                    !
   !=======================================================================================!

   !----- The sub-step's local error relative to its tolerance (<= 1 to accept). Implicit Euler's   !
   !      error is half the change in the inflows across the sub-step, h (F_end - F_start)/2. Taken !
   !      raw that over-states a stiff exchange that has already settled -- after rain, the 0.4 s  !
   !      soil-to-wood exchange -- by a factor of its stiffness, so it is filtered twice through    !
   !      (I - h dF/dW)^-1, which leaves slow modes alone and divides a stiff one's error by its     !
   !      stiffness squared: what step doubling would measure, at the price of two back-             !
   !      substitutions with the last Newton factorisation (§12.2 M1). Each store's error is        !
   !      weighed against atol of its scale plus rtol of its water (plan D8). --------------------!
   real(wp) function substep_error(net, h, ev0, ev1, scale_c, scale_g) result(err)
      type(water_network_t), intent(in) :: net
      real(wp),              intent(in) :: h
      type(net_eval_t),      intent(in) :: ev0, ev1
      real(wp),              intent(in) :: scale_c(:), scale_g(:,:)
      real(wp) :: v_c(net%nc), v_g(GROUP_SIZE, net%n), y_c(net%nc), y_g(GROUP_SIZE, net%n)
      integer(ik) :: pass
      v_c = h * (ev1%f_c - ev0%f_c)
      v_g = h * (ev1%f_g - ev0%f_g)
      where (net%core_fixed) v_c = 0.0_wp
      do pass = 1, 2
         call grouped_solve(net%sys, v_c / h, v_g / h, y_c, y_g)
         v_c = ev1%dw_c * y_c
         v_g = ev1%dw_g * y_g
         where (net%core_fixed) v_c = 0.0_wp
      end do
      err = maxval(0.5_wp * abs(v_c) / (net%atol * scale_c + net%rtol * abs(ev1%w_c)))
      if (net%n > 0_ik) err = max(err, maxval(0.5_wp * abs(v_g) / (net%atol * scale_g + net%rtol * abs(ev1%w_g))))
      if (err /= err) err = huge(1.0_wp)
   end function substep_error

   !----- Accept a sub-step: advance each store's water by the flows themselves, so water is        !
   !      conserved to round-off; spill the pond above its capacity; carry the pond's enthalpy; and  !
   !      add what each path moved to the step's amounts.                                             !
   !                                                                                                   !
   !      Newton balances every store only to a tolerance, and the flows carry that slack into the     !
   !      water they move. For a store hanging off the network by one path -- a cohort's leaves on     !
   !      its xylem, the pond on its face to the top layer -- the slack can be many times what the     !
   !      store holds (a near-leafless cohort's leaves hold 1e-26 kg/m2, its plant's tolerance is     !
   !      1e-9). So such a store takes its water from its own curve, and its one path carries the     !
   !      flow that closes its balance exactly; the store across that path (the wood, the top layer)   !
   !      absorbs the slack, which is small against what it holds. Every kilogram is still counted     !
   !      once. -----------------------------------------------------------------------------------!
   subroutine commit_substep(net, h, tau_c, tau_g, ev, spilled)
      type(water_network_t), intent(inout) :: net
      real(wp),              intent(in)    :: h, tau_c(:), tau_g(:,:)
      type(net_eval_t),      intent(in)    :: ev
      logical,               intent(out)   :: spilled   !< the pond spilled, so its unknown moved
      real(wp) :: w_mix, e_mix, t_mix, u_mix, amount, spill, t_pond, q_pond, sap(net%n)
      integer(ik) :: c, i

      !----- The flows that close the pond's and each leaf's balance at their curves' water. -------!
      q_pond = (net%w_core(WATER_NODE_POND) + h * net%rain - ev%w_c(WATER_NODE_POND)) / h
      do i = 1_ik, net%n
         sap(i) = (ev%w_g(WATER_LEAF, i) - net%w_group(WATER_LEAF, i)) / h + ev%e_leaf(i)
      end do

      !----- The pond's enthalpy, in the order its water moves: the rain joins it first, so what    !
      !      infiltrates leaves at the mixture's temperature; water the top layer pushes up arrives   !
      !      at that layer's temperature; then the spill leaves at the pond's own. --------------!
      w_mix = net%w_core(WATER_NODE_POND) + h * net%rain
      e_mix = net%e_pond + h * net%rain * internal_energy_liquid(net%t_rain)
      net%rain_enth = net%rain_enth + h * net%rain * internal_energy_liquid(net%t_rain)
      t_mix = net%t_rain
      if (w_mix > POND_TINY) t_mix = temp_of_liquid_enthalpy(e_mix / w_mix)
      u_mix = internal_energy_liquid(t_mix)
      if (q_pond >= 0.0_wp) then
         amount = h * q_pond
         net%infil = net%infil + amount ; net%infil_enth = net%infil_enth + amount * u_mix
         net%e_pond = e_mix - amount * u_mix
      else
         amount = -h * q_pond
         net%exfil = net%exfil + amount
         net%exfil_enth = net%exfil_enth + amount * internal_energy_liquid(net%t_soil(1))
         net%e_pond = e_mix + amount * internal_energy_liquid(net%t_soil(1))
      end if

      !----- Water, by the flows; the pond and the leaves at their curves, the top layer and the woods !
      !      across the paths that closed them. -----------------------------------------------!
      net%w_core  = net%w_core  + h * ev%f_c
      net%w_group = net%w_group + h * ev%f_g
      net%w_core(WATER_NODE_POND) = ev%w_c(WATER_NODE_POND)
      net%w_core(2) = net%w_core(2) + h * (q_pond - ev%q_pond)
      do i = 1_ik, net%n
         net%w_group(WATER_LEAF, i) = ev%w_g(WATER_LEAF, i)
         net%w_group(WATER_WOOD, i) = net%w_group(WATER_WOOD, i) - h * (sap(i) - ev%s_xyl(i))
      end do
      net%tau_core = tau_c ; net%tau_group = tau_g
      do c = 1_ik, net%nc
         if (net%core_fixed(c)) net%w_core(c) = ev%w_c(c)
      end do

      !----- Spill: water above the pond's capacity runs off at the pond's temperature. --------!
      spill = max(0.0_wp, net%w_core(WATER_NODE_POND) - net%w_pond_max)
      spilled = spill > 0.0_wp
      if (spilled) then
         t_pond = temp_of_liquid_enthalpy(net%e_pond / net%w_core(WATER_NODE_POND))
         net%spill = net%spill + spill
         net%spill_enth = net%spill_enth + spill * internal_energy_liquid(t_pond)
         net%e_pond = net%e_pond - spill * internal_energy_liquid(t_pond)
         net%w_core(WATER_NODE_POND) = net%w_pond_max
         net%tau_core(WATER_NODE_POND) = pond_unknown_from_water(net%w_pond_max)
      end if
      if (net%w_core(WATER_NODE_POND) <= POND_TINY) then
         net%e_pond = 0.0_wp
      end if

      !----- What each path moved. --------------------------------------------------------!
      net%face(1:net%nsl-1) = net%face(1:net%nsl-1) + h * ev%q_face(1:net%nsl-1)
      net%drain = net%drain + h * ev%q_drain
      net%evap  = net%evap  + h * ev%e_soil
      net%uptake  = net%uptake  + h * ev%u_root
      net%sapflow = net%sapflow + h * sap
      net%transp  = net%transp  + h * ev%e_leaf
   end subroutine commit_substep

   !----- Each store's water scale [kg/m2], against which its tolerances are set: a soil layer's or  !
   !      the pond's is a metre of water per metre of the (top) layer, so its tolerances read as     !
   !      water content; a cohort's is the water its wood and leaf cells hold at full turgor, shared  !
   !      by both stores so a leafless cohort's empty leaf is weighed against its plant. ---------!
   subroutine water_scales(net, hyd, scale_c, scale_g)
      type(water_network_t),      intent(in)  :: net
      type(hydro_params_table_t), intent(in)  :: hyd
      real(wp),                   intent(out) :: scale_c(:), scale_g(:,:)
      integer(ik) :: k, i
      real(wp)    :: s
      scale_c(WATER_NODE_POND) = rho_h2o * net%soil(1)%dz
      do k = 1_ik, net%nsl
         scale_c(1_ik + k) = rho_h2o * net%soil(k)%dz
      end do
      do i = 1_ik, net%n
         associate (lc => hyd%pft(net%pft(i))%leaf_curve, wc => hyd%pft(net%pft(i))%wood_curve)
            s = net%nplant(i) * ((1.0_wp - lc%apoplast_frac) * lc%water_sat * net%bleaf(i)                    &
                                + (1.0_wp - wc%apoplast_frac) * wc%water_sat * net%bwood(i))
         end associate
         scale_g(:, i) = max(s, tiny_num)
      end do
   end subroutine water_scales

   subroutine clear_amounts(net)
      type(water_network_t), intent(inout) :: net
      net%face = 0.0_wp ; net%uptake = 0.0_wp ; net%sapflow = 0.0_wp ; net%transp = 0.0_wp
      net%infil = 0.0_wp ; net%infil_enth = 0.0_wp ; net%exfil = 0.0_wp ; net%exfil_enth = 0.0_wp
      net%drain = 0.0_wp ; net%evap = 0.0_wp ; net%spill = 0.0_wp ; net%spill_enth = 0.0_wp
      net%rain_enth = 0.0_wp
      net%nsub = 0_ik ; net%nnewton = 0_ik ; net%nreject = 0_ik ; net%nfail = 0_ik
   end subroutine clear_amounts

   subroutine eval_allocate(ev, net)
      type(net_eval_t),      intent(inout) :: ev
      type(water_network_t), intent(in)    :: net
      if (allocated(ev%psi_c)) deallocate(ev%psi_c, ev%dpsi_c, ev%w_c, ev%dw_c, ev%k_c, ev%dk_c, ev%f_c,          &
                                          ev%psi_g, ev%dpsi_g, ev%w_g, ev%dw_g, ev%f_g, ev%q_face, ev%u_root,     &
                                          ev%s_xyl, ev%e_leaf)
      allocate(ev%psi_c(net%nc), ev%dpsi_c(net%nc), ev%w_c(net%nc), ev%dw_c(net%nc), ev%k_c(net%nc),       &
               ev%dk_c(net%nc), ev%f_c(net%nc))
      allocate(ev%psi_g(GROUP_SIZE, net%n), ev%dpsi_g(GROUP_SIZE, net%n), ev%w_g(GROUP_SIZE, net%n),        &
               ev%dw_g(GROUP_SIZE, net%n), ev%f_g(GROUP_SIZE, net%n))
      allocate(ev%q_face(net%nsl), ev%u_root(net%nsl, net%n), ev%s_xyl(net%n), ev%e_leaf(net%n))
      ev%q_face = 0.0_wp
   end subroutine eval_allocate

end module meds_water_network
