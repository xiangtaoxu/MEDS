! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_spac_solver -- one patch's soil and plant water over one fast step, solved together as    !
! one network (docs/science/soil_plant_water.md).                                                !
!                                                                                                !
! THE NETWORK (docs/science/soil_plant_water.md §1). Water sits in stores joined by paths:       !
!   * the soil nodes: the pond (node 1) and the soil layers (node 1 + k);                        !
!   * two stores per cohort, its wood and its leaves, which touch only each other and the        !
!     soil layers its roots reach;                                                               !
!   * paths: soil water flux between adjacent layers and between the pond and the top layer;     !
!     root water uptake from each layer into each cohort's wood; sapflow up each cohort's        !
!     xylem;                                                                                     !
!   * boundaries: rain into the pond, ground evaporation from the top layer, drainage from       !
!     the bottom one, and each cohort's transpiration, limited by what its leaves hold           !
!     (docs/science/soil_plant_water.md eq. (3)).                                                !
! How much each store holds is meds_soil_hydr_lib's and meds_plant_hydr_lib's (the               !
! *_water_relations); how fast water moves along each path is meds_plant_hydr_lib's              !
! (rhizosphere, xylem, transpiration) and meds_soil_water's (soil water flux, evaporation,       !
! drainage). This module only joins them.                                                        !
!                                                                                                !
! THE SOLVE (docs/science/soil_plant_water.md §3, eq. (2)). Implicit Euler on sub-steps. Each    !
! sub-step is a Newton solve of the water balance of every store at the sub-step's end,          !
!     (W(tau) - W_start) / h - (net inflow at tau) = 0,                                          !
! in each store's own unknown tau: a soil layer's saturation (dry) or potential (wet), a         !
! tissue's relative water content, the pond's head. Its linear systems are block arrowheads --   !
! the soil nodes bordered by the cohorts -- solved by eliminating each cohort onto the soil      !
! (meds_linear_algebra). Newton steps are bounded so no store is driven past empty, and         !
! backtracked until the balance improves. A sub-step's error is the implicit-Euler local         !
! error, estimated from the change in the inflows across the sub-step and filtered twice         !
! through the Newton matrix so a stiff exchange that has already settled does not read as        !
! error; it sets the next sub-step.                                                              !
!                                                                                                !
! WHAT IS KEPT EXACT. Each store's water is advanced by the flows themselves, W += h * inflow,   !
! never re-read from its curve, so every kilogram that leaves one store arrives in another or    !
! crosses a boundary: the solve conserves water to round-off whatever Newton's residual. The     !
! amounts each path moved over the step are reported for the energy integrator and the ledgers.  !
!==========================================================================================!
module meds_spac_solver
   use meds_kinds,           only : wp, ik
   use meds_constants,       only : grav_head, rho_h2o, tiny_num
   use meds_soil_hydr_lib,   only : soil_layer_curve_t, build_soil_layer_curve, soil_water_relations,  &
                                    soil_theta_to_unknown, pond_water_relations, pond_water_to_head
   use meds_plant_hydr_lib,  only : tissue_water_relations, tissue_water_to_rwc, tissue_rwc_to_unknown,  &
                                    tissue_psi_to_rwc, root_water_uptake_profile, xylem_sapflow
   use meds_plant_hydraulics, only : supply_limited_transpiration
   use meds_soil_water,      only : soil_water_flux, ground_evaporation_rate, bottom_flux
   use meds_linear_algebra,  only : arrowhead_system_t, arrowhead_allocate, arrowhead_factor, arrowhead_solve
   use meds_soil_types,      only : chydro_forcing_t
   use meds_soil_column_config, only : soil_params_t
   use meds_biophysics_opts, only : soil_opts_t
   use meds_plant_types,     only : hydro_args_table_t
   use meds_numerics,        only : sufficient_decrease, fraction_to_boundary
   use meds_therm_lib,       only : internal_energy_liquid, liquid_energy_to_temp
   implicit none
   private

   public :: spac_solver_t, spac_allocate, spac_set_state, spac_advance
   public :: spac_jacobian_error
   public :: POND, WOOD, LEAF

   integer(ik), parameter :: POND = 1_ik            !< the pond is soil node 1; soil layer k is node 1 + k
   integer(ik), parameter :: WOOD = 1_ik            !< a cohort's wood, its first store
   integer(ik), parameter :: LEAF = 2_ik            !< its leaves, its second
   integer(ik), parameter :: COHORT_STORES = 2_ik   !< stores per cohort

   !----- Newton. The balance counts as solved when every store's imbalance over the sub-step is  !
   !      below NEWTON_TOL of the step's absolute water tolerance -- far inside the error the     !
   !      sub-step itself is allowed. A Newton step may take a soil layer or a tissue at most      !
   !      STEP_MAX_SHARE of the way to empty (fraction to the boundary) and move its saturation or !
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

   !----- The network at one state: each store's potential [MPa] and water [kg/m2] with their       !
   !      slopes in its unknown, each soil node's conductivity [m/s], the net inflow into each store !
   !      [kg/m2/s], and the flow along each path [kg/m2/s]. ---------------------------------------!
   type :: spac_eval_t
      real(wp), allocatable :: psi_soil(:), dpsi_soil(:), water_soil(:), dwater_soil(:)   !< (nsoil)
      real(wp), allocatable :: k_soil(:), dk_soil(:), inflow_soil(:)                       !< (nsoil)
      real(wp), allocatable :: psi_coh(:,:), dpsi_coh(:,:), water_coh(:,:), dwater_coh(:,:) !< (2, ncoh)
      real(wp), allocatable :: inflow_coh(:,:)                                             !< (2, ncoh)
      real(wp), allocatable :: face_flux(:)        !< (nsl) down through the face below layer k
      real(wp), allocatable :: uptake(:,:)         !< (nsl, ncoh) from each layer into each cohort
      real(wp), allocatable :: sapflow(:), transp(:)   !< (ncoh)
      real(wp) :: pond_flux = 0.0_wp               !< pond -> top layer
      real(wp) :: drain = 0.0_wp, evap = 0.0_wp    !< out of the bottom; out of the top layer to the air
   end type spac_eval_t

   !----- One patch's soil-plant water: what it is built from, its state, and what it moved. -----!
   type :: spac_solver_t
      integer(ik) :: nsl   = 0_ik         !< soil layers
      integer(ik) :: ncoh  = 0_ik         !< cohorts
      integer(ik) :: nsoil = 0_ik         !< soil nodes: the pond and the nsl layers
      !----- The soil and the pond. -------------------------------------------------------!
      type(soil_layer_curve_t), allocatable :: layer(:)   !< (nsl) each layer's curve
      real(wp), allocatable :: z_node(:)        !< (nsl) [m, < 0] node depth
      real(wp), allocatable :: dz_node(:)       !< (nsl) [m] distance to the node below
      real(wp), allocatable :: face_weight(:)   !< (nsl) lower layer's weight in the face conductivity
      real(wp), allocatable :: t_soil(:)        !< (nsl) [K] temperature of water a layer gives up
      real(wp) :: rain       = 0.0_wp           !< [kg/m2/s] water reaching the pond
      real(wp) :: t_rain     = 298.15_wp        !< [K] the temperature that values it
      real(wp) :: w_pond_max = 5.0_wp           !< [kg/m2] pond capacity; above it, runoff
      real(wp) :: rtol = 1.0e-3_wp, atol = 1.0e-4_wp   !< sub-step error tolerances ([soil] rtol, atol)
      !----- The cohorts. ----------------------------------------------------------------!
      integer(ik), allocatable :: pft(:)        !< (ncoh) PFT, selecting curves and vulnerability
      real(wp), allocatable :: nplant(:)        !< (ncoh) [plants/m2]
      real(wp), allocatable :: bleaf(:), bwood(:)   !< (ncoh) [kgC/plant] leaf; sapwood + fine root
      real(wp), allocatable :: grav(:)          !< (ncoh) [MPa] hydrostatic head from wood to leaves
      real(wp), allocatable :: k_xylem(:)       !< (ncoh) [kg/s/MPa per m2] intact xylem conductance
      real(wp), allocatable :: e_demand(:)      !< (ncoh) [kg/m2/s] the stomata's transpiration demand
      real(wp), allocatable :: r_close(:)       !< (ncoh) [-] leaf RWC at which the stomata are shut
      real(wp), allocatable :: g_rhizo(:,:)     !< (nsl, ncoh) rhizosphere conductance per unit soil
                                                !<   conductivity [kg/s/MPa per m2, per m/s]
      !----- Soil nodes held at their potential: a plant-only solve prescribes the soil. -------!
      logical,  allocatable :: soil_fixed(:)    !< (nsoil)
      !----- State: each store's unknown and the water it holds [kg/m2]. -----------------------!
      real(wp), allocatable :: tau_soil(:), tau_coh(:,:)       !< (nsoil), (2, ncoh)
      real(wp), allocatable :: water_soil(:), water_coh(:,:)   !< (nsoil), (2, ncoh)
      real(wp) :: e_pond = 0.0_wp               !< [J/m2] the pond's enthalpy
      !----- What moved over the step [kg/m2, J/m2]. --------------------------------------!
      real(wp), allocatable :: face(:)          !< (nsl) down through the face below layer k
      real(wp), allocatable :: uptake(:,:)      !< (nsl, ncoh) from each layer into each cohort
      real(wp), allocatable :: sapflow(:), transp(:)   !< (ncoh)
      real(wp) :: infil = 0.0_wp, infil_enth = 0.0_wp  !< pond -> top layer, and its enthalpy
      real(wp) :: exfil = 0.0_wp, exfil_enth = 0.0_wp  !< top layer -> pond
      real(wp) :: drain = 0.0_wp, evap = 0.0_wp        !< out of the bottom; out of the top layer to the air
      real(wp) :: spill = 0.0_wp, spill_enth = 0.0_wp  !< pond runoff over its capacity
      real(wp) :: rain_enth = 0.0_wp                   !< enthalpy the rain brought
      !----- End-of-step potentials [MPa]. ----------------------------------------------!
      real(wp), allocatable :: psi_soil(:), psi_coh(:,:)   !< (nsoil), (2, ncoh)
      !----- Work done: accepted sub-steps, Newton iterations, rejected sub-steps, and sub-steps  !
      !      taken although they did not converge or meet the tolerance (counted, not hidden). ---!
      integer(ik) :: nsub = 0_ik, nnewton = 0_ik, nreject = 0_ik, nfail = 0_ik
      !----- Scratch, kept from step to step so the arrays are allocated again only when the shape  !
      !      changes. ev(i_start) is the network at a sub-step's start; Newton evaluates its iterates  !
      !      into the other. sys holds the flows' Jacobian of the latest evaluation (its factorisation !
      !      is kept apart). A sub-step starts where the last accepted one ended, so while            !
      !      start_ready holds, its start's evaluation and Jacobian are already there. ---------------!
      type(arrowhead_system_t) :: sys
      type(spac_eval_t)        :: ev(2)
      integer(ik)              :: i_start = 1_ik
      logical                  :: start_ready = .false.
   end type spac_solver_t

contains

   !----- Size the solver for nsl soil layers and ncoh cohorts, keeping the storage when the shape  !
   !      is unchanged, and reset what the caller fills. ------------------------------------------!
   subroutine spac_allocate(spac, nsl, ncoh)
      type(spac_solver_t), intent(inout) :: spac
      integer(ik),         intent(in)    :: nsl, ncoh
      if (.not. allocated(spac%layer) .or. spac%nsl /= nsl .or. spac%ncoh /= ncoh) then
         spac%nsl = nsl ; spac%ncoh = ncoh ; spac%nsoil = 1_ik + nsl
         if (allocated(spac%layer)) deallocate(spac%layer, spac%z_node, spac%dz_node, spac%face_weight,          &
              spac%t_soil, spac%pft, spac%nplant, spac%bleaf, spac%bwood, spac%grav, spac%k_xylem,               &
              spac%e_demand, spac%r_close, spac%g_rhizo, spac%soil_fixed, spac%tau_soil, spac%tau_coh,           &
              spac%water_soil, spac%water_coh, spac%face, spac%uptake, spac%sapflow, spac%transp,                &
              spac%psi_soil, spac%psi_coh)
         allocate(spac%layer(nsl), spac%z_node(nsl), spac%dz_node(nsl), spac%face_weight(nsl), spac%t_soil(nsl))
         allocate(spac%pft(ncoh), spac%nplant(ncoh), spac%bleaf(ncoh), spac%bwood(ncoh), spac%grav(ncoh),      &
                  spac%k_xylem(ncoh), spac%e_demand(ncoh), spac%r_close(ncoh), spac%g_rhizo(nsl, ncoh))
         allocate(spac%soil_fixed(spac%nsoil), spac%tau_soil(spac%nsoil), spac%tau_coh(COHORT_STORES, ncoh),    &
                  spac%water_soil(spac%nsoil), spac%water_coh(COHORT_STORES, ncoh))
         allocate(spac%face(nsl), spac%uptake(nsl, ncoh), spac%sapflow(ncoh), spac%transp(ncoh),                &
                  spac%psi_soil(spac%nsoil), spac%psi_coh(COHORT_STORES, ncoh))
         call eval_allocate(spac%ev(1), spac) ; call eval_allocate(spac%ev(2), spac)
      end if
      spac%t_soil = 298.15_wp ; spac%soil_fixed = .false. ; spac%g_rhizo = 0.0_wp
      spac%e_demand = 0.0_wp ; spac%r_close = 0.5_wp
   end subroutine spac_allocate

   !----- Set the solver's state from the model's: soil water content, the pond's water and       !
   !      enthalpy, and each cohort's leaf and wood water [kg per plant]. The cohort inputs         !
   !      (pft ... g_rhizo) must already be filled. Each store's unknown is read from its water; a  !
   !      leafless cohort's leaves hold none, so their unknown starts at the wood's potential less  !
   !      the head to the crown, where an empty leaf comes to rest. ----------------------------!
   subroutine spac_set_state(spac, params, hyd, theta, w_pond, e_pond, leaf_water, wood_water)
      type(spac_solver_t),        intent(inout) :: spac
      type(soil_params_t),        intent(in)    :: params
      type(hydro_args_table_t),   intent(in)    :: hyd
      real(wp),                   intent(in)    :: theta(:)                       !< (nsl) [m3/m3]
      real(wp),                   intent(in)    :: w_pond, e_pond                 !< [kg/m2], [J/m2]
      real(wp),                   intent(in)    :: leaf_water(:), wood_water(:)   !< (ncoh) [kg/plant]
      integer(ik) :: k, i
      real(wp)    :: psi, dpsi, w, dw, kk, dk, head1
      do k = 1_ik, spac%nsl
         call build_soil_layer_curve(params%retention, params%theta_sat(k), params%theta_res(k),              &
                                     params%curve_par_a(k), params%curve_par_n(k), params%ksat(k), params%dz(k), &
                                     spac%layer(k))
         spac%z_node(k) = params%z_node(k)
         if (k < spac%nsl) then
            spac%dz_node(k)     = params%dz_node(k)
            spac%face_weight(k) = params%dz(k) / (params%dz(k) + params%dz(k+1))
         end if
         spac%tau_soil(1_ik + k)   = soil_theta_to_unknown(spac%layer(k), theta(k))
         spac%water_soil(1_ik + k) = rho_h2o * params%dz(k) * theta(k)
      end do
      !----- The pond starts at its water's head, but no higher than the head at which it would     !
      !      neither give nor take water from the top layer. Over dry soil the pond face drains tens of !
      !      millimetres an hour, so a pond holding a trace of water empties within the step: starting !
      !      it at that trace's head would put the step's first Newton iterate at a flow thousands of   !
      !      times what the pond holds, far outside the linear model's reach. ----------------------!
      spac%water_soil(POND) = w_pond
      spac%e_pond = e_pond
      call soil_water_relations(spac%layer(1), spac%tau_soil(2), psi, dpsi, w, dw, kk, dk)
      head1 = psi / grav_head
      spac%tau_soil(POND) = min(pond_water_to_head(w_pond), head1 + spac%z_node(1))
      do i = 1_ik, spac%ncoh
         associate (hp => hyd%pft(spac%pft(i)))
            spac%tau_coh(WOOD, i)   = tissue_rwc_to_unknown(                                               &
                                         tissue_water_to_rwc(hp%wood_curve, spac%bwood(i), wood_water(i)))
            spac%water_coh(WOOD, i) = wood_water(i) * spac%nplant(i)
            if (hp%leaf_curve%water_sat * spac%bleaf(i) > tiny_num) then
               spac%tau_coh(LEAF, i)   = tissue_rwc_to_unknown(                                            &
                                            tissue_water_to_rwc(hp%leaf_curve, spac%bleaf(i), leaf_water(i)))
               spac%water_coh(LEAF, i) = leaf_water(i) * spac%nplant(i)
            else
               call tissue_water_relations(hp%wood_curve, spac%bwood(i), spac%nplant(i), spac%tau_coh(WOOD, i), &
                                           psi, dpsi, w, dw)
               spac%tau_coh(LEAF, i)   = tissue_rwc_to_unknown(tissue_psi_to_rwc(psi - spac%grav(i),             &
                                            hp%leaf_curve%pi0, hp%leaf_curve%elastic_mod))
               spac%water_coh(LEAF, i) = 0.0_wp
            end if
         end associate
      end do
   end subroutine spac_set_state

   !----- Advance over dt. On return the state holds the end-of-step water and unknowns, and the     !
   !      amounts each path moved. -----------------------------------------------------------!
   subroutine spac_advance(spac, params, opts, forcing, hyd, dt)
      type(spac_solver_t),        intent(inout) :: spac
      type(soil_params_t),        intent(in)    :: params
      type(soil_opts_t),          intent(in)    :: opts
      type(chydro_forcing_t),     intent(in)    :: forcing
      type(hydro_args_table_t),   intent(in)    :: hyd
      real(wp),                   intent(in)    :: dt

      real(wp) :: t, h, h_min, err, scale_soil(spac%nsoil), scale_coh(COHORT_STORES, spac%ncoh)
      real(wp) :: tau_soil(spac%nsoil), tau_coh(COHORT_STORES, spac%ncoh)
      integer(ik) :: trials, iters, i_sol
      logical  :: converged, last, spilled

      call clear_amounts(spac)
      call arrowhead_allocate(spac%sys, spac%nsoil, spac%ncoh, COHORT_STORES)
      call water_scales(spac, hyd, scale_soil, scale_coh)
      spac%start_ready = .false.

      t = 0.0_wp ; h = dt ; h_min = H_MIN_SHARE * dt ; trials = 0_ik
      do while (t < dt * (1.0_wp - 1.0e-12_wp))
         h = min(h, dt - t)
         last = (h <= h_min) .or. (trials >= TRIALS_MAX)
         trials = trials + 1_ik
         tau_soil = spac%tau_soil ; tau_coh = spac%tau_coh
         call newton_substep(spac, params, opts, forcing, hyd, h, scale_soil, scale_coh, tau_soil, tau_coh,     &
                             converged, iters, i_sol)
         spac%nnewton = spac%nnewton + iters
         if (.not. converged .and. .not. last) then
            spac%nreject = spac%nreject + 1_ik
            h = max(h * NEWTON_FAIL_CUT, h_min)
            cycle
         end if
         if (converged) then
            err = substep_error(spac, h, spac%ev(spac%i_start), spac%ev(i_sol), scale_soil, scale_coh)
         else
            err = huge(1.0_wp)
         end if
         if (err > 1.0_wp .and. .not. last) then
            spac%nreject = spac%nreject + 1_ik
            h = max(h * max(SHRINK_MIN, SAFETY * err ** (-0.5_wp)), h_min)
            cycle
         end if
         if (err > 1.0_wp) spac%nfail = spac%nfail + 1_ik
         call commit_substep(spac, h, tau_soil, tau_coh, spac%ev(i_sol), spilled)
         !----- The next sub-step starts where this one ended: its start is this one's solution, whose  !
         !      Jacobian sys still holds -- unless the spill moved the pond. ----------------------------!
         spac%i_start = i_sol
         spac%start_ready = .not. spilled
         t = t + h
         spac%nsub = spac%nsub + 1_ik
         if (err <= 1.0_wp)                                                                                   &
            h = h * min(GROW_MAX, max(SHRINK_MIN, SAFETY * max(err, 1.0e-10_wp) ** (-0.5_wp))) ! clamp-ok: controller floor
      end do
      spac%psi_soil = spac%ev(spac%i_start)%psi_soil
      spac%psi_coh  = spac%ev(spac%i_start)%psi_coh
   end subroutine spac_advance

   !----- How far the assembled Jacobian is from the residual's own derivative, at the current     !
   !      state for a sub-step of h: the largest difference between J v and a central difference   !
   !      of the residual along v, relative to the largest entry of either, with v a fixed pattern !
   !      of small changes in every unknown. For the test suite: it checks the whole Jacobian     !
   !      against the flows it is built from. -----------------------------------------------------!
   real(wp) function spac_jacobian_error(spac, params, opts, forcing, hyd, h) result(err)
      type(spac_solver_t),        intent(inout) :: spac
      type(soil_params_t),        intent(in)    :: params
      type(soil_opts_t),          intent(in)    :: opts
      type(chydro_forcing_t),     intent(in)    :: forcing
      type(hydro_args_table_t),   intent(in)    :: hyd
      real(wp),                   intent(in)    :: h
      type(spac_eval_t) :: ev
      real(wp) :: v_soil(spac%nsoil), v_coh(COHORT_STORES, spac%ncoh)
      real(wp) :: jv_soil(spac%nsoil), jv_coh(COHORT_STORES, spac%ncoh)
      real(wp) :: rp_soil(spac%nsoil), rp_coh(COHORT_STORES, spac%ncoh)
      real(wp) :: rm_soil(spac%nsoil), rm_coh(COHORT_STORES, spac%ncoh)
      real(wp) :: add_soil(spac%nsoil), add_coh(COHORT_STORES, spac%ncoh)
      real(wp) :: eps, big
      integer(ik) :: c, i, j, l
      call eval_allocate(ev, spac)
      call arrowhead_allocate(spac%sys, spac%nsoil, spac%ncoh, COHORT_STORES)
      do c = 1_ik, spac%nsoil
         v_soil(c) = 1.0e-2_wp * (1.0_wp + 0.1_wp * real(c, wp)) * max(abs(spac%tau_soil(c)), 1.0e-3_wp) ! clamp-ok: test step
      end do
      do i = 1_ik, spac%ncoh
         do j = 1_ik, COHORT_STORES
            v_coh(j, i) = -1.0e-2_wp * (1.0_wp + 0.07_wp * real(i + j, wp)) * spac%tau_coh(j, i)
         end do
      end do
      call evaluate(spac, params, opts, forcing, hyd, spac%tau_soil, spac%tau_coh, ev, .true.)
      call storage_diagonal(spac, h, ev, add_soil, add_coh)
      jv_soil = matmul(spac%sys%a, v_soil) + add_soil * v_soil ; jv_coh = 0.0_wp
      do i = 1_ik, spac%ncoh
         jv_soil = jv_soil + spac%sys%b(:, i) * v_coh(WOOD, i)
         do l = 1_ik, COHORT_STORES
            jv_coh(l, i) = sum(spac%sys%d(l, :, i) * v_coh(:, i)) + add_coh(l, i) * v_coh(l, i)
         end do
         jv_coh(WOOD, i) = jv_coh(WOOD, i) + sum(spac%sys%c(:, i) * v_soil)
      end do
      eps = 1.0e-5_wp
      call evaluate(spac, params, opts, forcing, hyd, spac%tau_soil + eps * v_soil, spac%tau_coh + eps * v_coh,   &
                    ev, .false.)
      call residual(spac, h, ev, rp_soil, rp_coh)
      call evaluate(spac, params, opts, forcing, hyd, spac%tau_soil - eps * v_soil, spac%tau_coh - eps * v_coh,   &
                    ev, .false.)
      call residual(spac, h, ev, rm_soil, rm_coh)
      rp_soil = (rp_soil - rm_soil) / (2.0_wp * eps) ; rp_coh = (rp_coh - rm_coh) / (2.0_wp * eps)
      where (spac%soil_fixed) jv_soil = 0.0_wp
      big = max(maxval(abs(jv_soil)), maxval(abs(rp_soil)), tiny_num)
      err = maxval(abs(jv_soil - rp_soil)) / big
      if (spac%ncoh > 0_ik) then
         big = max(maxval(abs(jv_coh)), maxval(abs(rp_coh)), tiny_num)
         err = max(err, maxval(abs(jv_coh - rp_coh)) / big)
      end if
   end function spac_jacobian_error

   !=======================================================================================!
   !     One sub-step: Newton on the end-of-sub-step water balance.                          !
   !=======================================================================================!
   subroutine newton_substep(spac, params, opts, forcing, hyd, h, scale_soil, scale_coh, tau_soil, tau_coh,       &
                             converged, iters, i_sol)
      type(spac_solver_t),        intent(inout) :: spac
      type(soil_params_t),        intent(in)    :: params
      type(soil_opts_t),          intent(in)    :: opts
      type(chydro_forcing_t),     intent(in)    :: forcing
      type(hydro_args_table_t),   intent(in)    :: hyd
      real(wp),                   intent(in)    :: h, scale_soil(:), scale_coh(:,:)
      real(wp),                   intent(inout) :: tau_soil(:), tau_coh(:,:)   !< start of sub-step in, solution out
      logical,                    intent(out)   :: converged
      integer(ik),                intent(out)   :: iters
      integer(ik),                intent(out)   :: i_sol   !< which of spac%ev holds the solution
      real(wp) :: res_soil(spac%nsoil), res_coh(COHORT_STORES, spac%ncoh)
      real(wp) :: step_soil(spac%nsoil), step_coh(COHORT_STORES, spac%ncoh)
      real(wp) :: try_soil(spac%nsoil), try_coh(COHORT_STORES, spac%ncoh)
      real(wp) :: res2_soil(spac%nsoil), res2_coh(COHORT_STORES, spac%ncoh)
      real(wp) :: merit, merit_try, worst, lam, lam_max
      integer(ik) :: it, ib, i_try
      logical  :: ok

      converged = .false. ; iters = 0_ik
      !----- The network at the start, evaluated unless it already is; the iterates are evaluated into  !
      !      the other slot, and the solution is the start until one is taken. ------------------------!
      i_sol = spac%i_start ; i_try = 3_ik - spac%i_start
      if (.not. spac%start_ready) then
         call evaluate(spac, params, opts, forcing, hyd, tau_soil, tau_coh, spac%ev(i_sol), .true.)
         spac%start_ready = .true.
      end if
      call residual(spac, h, spac%ev(i_sol), res_soil, res_coh)
      merit = residual_merit(spac, h, res_soil, res_coh, scale_soil, scale_coh)
      worst = residual_worst(spac, h, res_soil, res_coh, scale_soil, scale_coh)
      if (merit /= merit) return
      do it = 1_ik, NEWTON_MAX
         if (worst <= 1.0_wp) exit
         iters = it
         call newton_factor(spac, h, spac%ev(i_sol), ok)
         if (.not. ok) return
         call arrowhead_solve(spac%sys, -res_soil, -res_coh, step_soil, step_coh)
         lam_max = step_limit(spac, tau_soil, tau_coh, step_soil, step_coh)
         lam = lam_max
         !----- A trial's evaluation replaces the Jacobian sys held, so the start's is gone. ---------!
         spac%start_ready = .false.
         do ib = 1_ik, BACKTRACK_MAX
            try_soil = tau_soil + lam * step_soil ; try_coh = tau_coh + lam * step_coh
            call evaluate(spac, params, opts, forcing, hyd, try_soil, try_coh, spac%ev(i_try), .true.)
            call residual(spac, h, spac%ev(i_try), res2_soil, res2_coh)
            merit_try = residual_merit(spac, h, res2_soil, res2_coh, scale_soil, scale_coh)
            if (merit_try == merit_try .and. sufficient_decrease(merit_try, merit, lam)) exit
            lam = 0.5_wp * lam
         end do
         if (merit_try /= merit_try) return
         tau_soil = try_soil ; tau_coh = try_coh ; res_soil = res2_soil ; res_coh = res2_coh
         merit = merit_try ; i_sol = i_try
         worst = residual_worst(spac, h, res_soil, res_coh, scale_soil, scale_coh)
      end do
      !----- Factor the Newton matrix at the solution: the error estimate solves with it. -----------!
      call newton_factor(spac, h, spac%ev(i_sol), ok)
      converged = ok .and. worst <= 1.0_wp
   end subroutine newton_substep

   !----- The water-balance residual of every store [kg/m2/s]: (W(tau) - W_start)/h - inflow.  !
   !      A prescribed soil node has none. -------------------------------------------------!
   subroutine residual(spac, h, ev, res_soil, res_coh)
      type(spac_solver_t), intent(in)  :: spac
      real(wp),            intent(in)  :: h
      type(spac_eval_t),   intent(in)  :: ev
      real(wp),            intent(out) :: res_soil(:), res_coh(:,:)
      res_soil = (ev%water_soil - spac%water_soil) / h - ev%inflow_soil
      res_coh  = (ev%water_coh - spac%water_coh) / h - ev%inflow_coh
      where (spac%soil_fixed) res_soil = 0.0_wp
   end subroutine residual

   !----- Each store's imbalance over the sub-step against NEWTON_TOL of its absolute water        !
   !      tolerance. residual_worst, below one, means the balance is solved; residual_merit, the sum of !
   !      squares, is what the line search asks to fall -- a Newton step always points down it,      !
   !      which is not true of the worst store alone. ---------------------------------------------!
   pure real(wp) function residual_worst(spac, h, res_soil, res_coh, scale_soil, scale_coh) result(m)
      type(spac_solver_t), intent(in) :: spac
      real(wp),            intent(in) :: h, res_soil(:), res_coh(:,:), scale_soil(:), scale_coh(:,:)
      real(wp) :: tol
      tol = NEWTON_TOL * spac%atol
      m = max(maxval(abs(res_soil) * h / (tol * scale_soil)), 0.0_wp)
      if (spac%ncoh > 0_ik) m = max(m, maxval(abs(res_coh) * h / (tol * scale_coh)))
   end function residual_worst

   pure real(wp) function residual_merit(spac, h, res_soil, res_coh, scale_soil, scale_coh) result(m)
      type(spac_solver_t), intent(in) :: spac
      real(wp),            intent(in) :: h, res_soil(:), res_coh(:,:), scale_soil(:), scale_coh(:,:)
      real(wp) :: tol
      tol = NEWTON_TOL * spac%atol
      m = sum((res_soil * h / (tol * scale_soil)) ** 2)
      if (spac%ncoh > 0_ik) m = m + sum((res_coh * h / (tol * scale_coh)) ** 2)
   end function residual_merit

   !----- The largest share of a Newton step every store can take: no soil layer or tissue is      !
   !      driven more than STEP_MAX_SHARE of the way to empty, and none on the dry side moves its      !
   !      saturation or relative water content by more than STEP_MAX_DRY. A wet soil layer's unknown   !
   !      is its potential, so only the first applies to it: it keeps a step that would carry it past  !
   !      the branch point and below empty inside the curve. The pond's head is not bounded. -------!
   pure real(wp) function step_limit(spac, tau_soil, tau_coh, step_soil, step_coh) result(lam)
      type(spac_solver_t), intent(in) :: spac
      real(wp),            intent(in) :: tau_soil(:), tau_coh(:,:), step_soil(:), step_coh(:,:)
      integer(ik) :: k, i, j
      lam = 1.0_wp
      do k = 1_ik, spac%nsl
         lam = min(lam, fraction_to_boundary(tau_soil(1_ik + k), step_soil(1_ik + k), STEP_MAX_SHARE))
         if (tau_soil(1_ik + k) <= spac%layer(k)%se_switch) call cap_dry_change(tau_soil(1_ik + k), step_soil(1_ik + k), lam)
      end do
      do i = 1_ik, spac%ncoh
         do j = 1_ik, COHORT_STORES
            lam = min(lam, fraction_to_boundary(tau_coh(j, i), step_coh(j, i), STEP_MAX_SHARE))
            call cap_dry_change(tau_coh(j, i), step_coh(j, i), lam)
         end do
      end do
   contains
      !----- The cap on the change applies below full saturation (an RWC or Se of one), where the     !
      !      potential steepens towards empty; above it the curve is gentle and Newton moves freely. ---!
      pure subroutine cap_dry_change(tau, d, lam)
         real(wp), intent(in)    :: tau, d
         real(wp), intent(inout) :: lam
         if (tau < 1.0_wp .and. abs(d) > STEP_MAX_DRY) lam = min(lam, STEP_MAX_DRY / abs(d))
      end subroutine cap_dry_change
   end function step_limit

   !=======================================================================================!
   !     The network at one state: storage, flows, and (optionally) the Jacobian.            !
   !=======================================================================================!
   subroutine evaluate(spac, params, opts, forcing, hyd, tau_soil, tau_coh, ev, with_jacobian)
      type(spac_solver_t),        intent(inout) :: spac
      type(soil_params_t),        intent(in)    :: params
      type(soil_opts_t),          intent(in)    :: opts
      type(chydro_forcing_t),     intent(in)    :: forcing
      type(hydro_args_table_t),   intent(in)    :: hyd
      real(wp),                   intent(in)    :: tau_soil(:), tau_coh(:,:)
      type(spac_eval_t),          intent(inout) :: ev
      logical,                    intent(in)    :: with_jacobian
      integer(ik) :: k, i, ca, cb, c
      real(wp) :: q, dq_pa, dq_pb, dq_ka, dq_kb, dqa, dqb, u, s, ds_w, ds_l
      real(wp) :: e, de_dtheta, de_dpsi, de, qbot, dqbot_dpsi, dqbot_dk, theta1, head
      real(wp) :: psi_w, dpsi_w, u_sum, dqb_sum, du_dpsi(spac%nsl), du_dk(spac%nsl)

      ev%inflow_soil = 0.0_wp ; ev%inflow_coh = 0.0_wp
      !----- The Jacobian's entries are written afresh: the soil block's band, which the flows add  !
      !      into, is cleared; every other entry an evaluation touches is set outright below, and the  !
      !      rest stay zero from arrowhead_allocate. ---------------------------------------------------!
      if (with_jacobian) then
         do c = 1_ik, spac%nsoil
            spac%sys%a(c, c) = 0.0_wp
            if (c < spac%nsoil) then
               spac%sys%a(c, c + 1_ik) = 0.0_wp ; spac%sys%a(c + 1_ik, c) = 0.0_wp
            end if
         end do
      end if

      !----- Storage. ------------------------------------------------------------------!
      call pond_water_relations(tau_soil(POND), ev%psi_soil(1), ev%dpsi_soil(1), ev%water_soil(1),             &
                                ev%dwater_soil(1))
      ev%k_soil(1) = params%ksat(1) ; ev%dk_soil(1) = 0.0_wp      ! the pond conducts as saturated soil
      do k = 1_ik, spac%nsl
         c = 1_ik + k
         call soil_water_relations(spac%layer(k), tau_soil(c), ev%psi_soil(c), ev%dpsi_soil(c), ev%water_soil(c), &
                                   ev%dwater_soil(c), ev%k_soil(c), ev%dk_soil(c))
      end do
      do i = 1_ik, spac%ncoh
         associate (hp => hyd%pft(spac%pft(i)))
            call tissue_water_relations(hp%wood_curve, spac%bwood(i), spac%nplant(i), tau_coh(WOOD, i),          &
                                        ev%psi_coh(WOOD, i), ev%dpsi_coh(WOOD, i), ev%water_coh(WOOD, i),        &
                                        ev%dwater_coh(WOOD, i))
            call tissue_water_relations(hp%leaf_curve, spac%bleaf(i), spac%nplant(i), tau_coh(LEAF, i),          &
                                        ev%psi_coh(LEAF, i), ev%dpsi_coh(LEAF, i), ev%water_coh(LEAF, i),        &
                                        ev%dwater_coh(LEAF, i))
         end associate
      end do

      !----- The pond and the top layer: Darcy, the pond a saturated layer as thick as the top one. !
      ca = POND ; cb = 2_ik
      call soil_water_flux(ev%k_soil(ca), ev%k_soil(cb), 0.5_wp, ev%psi_soil(ca), ev%psi_soil(cb), -spac%z_node(1),  &
                           -spac%z_node(1), q, dq_pa, dq_pb, dq_ka, dq_kb)
      ev%pond_flux = q
      call add_soil_flow(ca, cb, q, dq_pa * ev%dpsi_soil(ca), dq_pb * ev%dpsi_soil(cb) + dq_kb * ev%dk_soil(cb))

      !----- Between soil layers. --------------------------------------------------------!
      do k = 1_ik, spac%nsl - 1_ik
         ca = 1_ik + k ; cb = 2_ik + k
         call soil_water_flux(ev%k_soil(ca), ev%k_soil(cb), spac%face_weight(k), ev%psi_soil(ca), ev%psi_soil(cb), &
                              spac%dz_node(k), spac%dz_node(k), q, dq_pa, dq_pb, dq_ka, dq_kb)
         ev%face_flux(k) = q
         dqa = dq_pa * ev%dpsi_soil(ca) + dq_ka * ev%dk_soil(ca)
         dqb = dq_pb * ev%dpsi_soil(cb) + dq_kb * ev%dk_soil(cb)
         call add_soil_flow(ca, cb, q, dqa, dqb)
      end do

      !----- Rain into the pond; evaporation from the top layer; drainage from the bottom one. ----!
      ev%inflow_soil(POND) = ev%inflow_soil(POND) + spac%rain
      c = 2_ik
      theta1 = ev%water_soil(c) / (rho_h2o * spac%layer(1)%dz)
      call ground_evaporation_rate(theta1, ev%psi_soil(c) / grav_head, params, forcing, opts, e, de_dtheta, de_dpsi)
      ev%evap = e
      ev%inflow_soil(c) = ev%inflow_soil(c) - e
      if (with_jacobian) spac%sys%a(c, c) = spac%sys%a(c, c)                                                     &
                         + de_dtheta * ev%dwater_soil(c) / (rho_h2o * spac%layer(1)%dz)                        &
                         + de_dpsi * ev%dpsi_soil(c) / grav_head
      c = 1_ik + spac%nsl
      head = ev%psi_soil(c) / grav_head
      call bottom_flux(params, opts, spac%nsl, head, ev%k_soil(c), qbot, dqbot_dpsi, dqbot_dk)
      ev%drain = rho_h2o * qbot
      ev%inflow_soil(c) = ev%inflow_soil(c) - ev%drain
      if (with_jacobian) spac%sys%a(c, c) = spac%sys%a(c, c)                                                     &
                         + rho_h2o * (dqbot_dpsi * ev%dpsi_soil(c) / grav_head + dqbot_dk * ev%dk_soil(c))

      !----- Each cohort: uptake from its layers, sapflow to its leaves, and transpiration. The      !
      !      cohort's uptake comes from its whole root profile at once; its own sums are kept in       !
      !      scalars over its layers and written once. ------------------------------------------!
      do i = 1_ik, spac%ncoh
         psi_w = ev%psi_coh(WOOD, i) ; dpsi_w = ev%dpsi_coh(WOOD, i)
         call root_water_uptake_profile(spac%nsl, spac%g_rhizo(1:spac%nsl, i), ev%k_soil(2:spac%nsoil),             &
                                        ev%psi_soil(2:spac%nsoil), spac%z_node(1:spac%nsl), psi_w,                &
                                        ev%uptake(1:spac%nsl, i), du_dpsi, du_dk)
         u_sum = 0.0_wp ; dqb_sum = 0.0_wp
         do k = 1_ik, spac%nsl
            c = 1_ik + k
            if (spac%g_rhizo(k, i) <= 0.0_wp) then
               if (with_jacobian) then
                  spac%sys%b(c, i) = 0.0_wp ; spac%sys%c(c, i) = 0.0_wp
               end if
               cycle
            end if
            u = ev%uptake(k, i)
            ev%inflow_soil(c) = ev%inflow_soil(c) - u
            u_sum = u_sum + u
            if (with_jacobian) then
               dqa = du_dpsi(k) * ev%dpsi_soil(c) + du_dk(k) * ev%dk_soil(c)
               dqb = -du_dpsi(k) * dpsi_w
               spac%sys%a(c, c) = spac%sys%a(c, c) + dqa
               spac%sys%b(c, i) = dqb
               spac%sys%c(c, i) = -dqa
               dqb_sum = dqb_sum - dqb
            end if
         end do
         ev%inflow_coh(WOOD, i) = ev%inflow_coh(WOOD, i) + u_sum
         associate (hp => hyd%pft(spac%pft(i)))
            call xylem_sapflow(spac%k_xylem(i), ev%psi_coh(WOOD, i), ev%psi_coh(LEAF, i), spac%grav(i),           &
                               hp%wood_psi50, hp%wood_kexp, hp%vuln_table, s, ds_w, ds_l)
         end associate
         ev%sapflow(i) = s
         ev%inflow_coh(WOOD, i) = ev%inflow_coh(WOOD, i) - s
         ev%inflow_coh(LEAF, i) = ev%inflow_coh(LEAF, i) + s
         call supply_limited_transpiration(spac%e_demand(i), tau_coh(LEAF, i), spac%r_close(i), e, de)
         ev%transp(i) = e
         ev%inflow_coh(LEAF, i) = ev%inflow_coh(LEAF, i) - e
         if (with_jacobian) then
            dqa = ds_w * ev%dpsi_coh(WOOD, i)
            dqb = ds_l * ev%dpsi_coh(LEAF, i)
            spac%sys%d(WOOD, WOOD, i) = dqb_sum + dqa
            spac%sys%d(WOOD, LEAF, i) = dqb
            spac%sys%d(LEAF, WOOD, i) = -dqa
            spac%sys%d(LEAF, LEAF, i) = -dqb + de
         end if
      end do

   contains

      !----- A flow q from soil node a to soil node b, with dq/dtau_a and dq/dtau_b: the inflows   !
      !      and, when wanted, the Jacobian of the residual (W - W0)/h - inflow. ----------------!
      subroutine add_soil_flow(a, b, q, dqa, dqb)
         integer(ik), intent(in) :: a, b
         real(wp),    intent(in) :: q, dqa, dqb
         ev%inflow_soil(a) = ev%inflow_soil(a) - q
         ev%inflow_soil(b) = ev%inflow_soil(b) + q
         if (.not. with_jacobian) return
         spac%sys%a(a, a) = spac%sys%a(a, a) + dqa
         spac%sys%a(a, b) = spac%sys%a(a, b) + dqb
         spac%sys%a(b, a) = spac%sys%a(b, a) - dqa
         spac%sys%a(b, b) = spac%sys%a(b, b) - dqb
      end subroutine add_soil_flow

   end subroutine evaluate

   !----- The storage term dW/dtau / h each store adds to the Newton matrix's diagonal, and none    !
   !      for a prescribed soil node. ----------------------------------------------------------!
   subroutine storage_diagonal(spac, h, ev, add_soil, add_coh)
      type(spac_solver_t), intent(in)  :: spac
      real(wp),            intent(in)  :: h
      type(spac_eval_t),   intent(in)  :: ev
      real(wp),            intent(out) :: add_soil(:), add_coh(:,:)
      add_soil = ev%dwater_soil / h
      add_coh  = ev%dwater_coh / h
      where (spac%soil_fixed) add_soil = 0.0_wp
   end subroutine storage_diagonal

   !----- Factor the Newton matrix: the flows' Jacobian of the latest evaluation (sys) plus the     !
   !      storage terms of ev. A prescribed soil node's row is the identity, so its unknown stays     !
   !      put; making it so edits the assembled rows, which the next evaluation writes afresh. ----!
   subroutine newton_factor(spac, h, ev, ok)
      type(spac_solver_t), intent(inout) :: spac
      real(wp),            intent(in)    :: h
      type(spac_eval_t),   intent(in)    :: ev
      logical,             intent(out)   :: ok
      real(wp)    :: add_soil(spac%nsoil), add_coh(COHORT_STORES, spac%ncoh)
      integer(ik) :: c
      do c = 1_ik, spac%nsoil
         if (spac%soil_fixed(c)) then
            spac%sys%a(c, :) = 0.0_wp ; spac%sys%a(c, c) = 1.0_wp
            spac%sys%b(c, :) = 0.0_wp
         end if
      end do
      call storage_diagonal(spac, h, ev, add_soil, add_coh)
      call arrowhead_factor(spac%sys, add_soil, add_coh, ok)
   end subroutine newton_factor

   !=======================================================================================!
   !     Error, commit, and bookkeeping.                                                    !
   !=======================================================================================!

   !----- The sub-step's local error relative to its tolerance (<= 1 to accept). Implicit Euler's   !
   !      error is half the change in the inflows across the sub-step, h (F_end - F_start)/2. Taken !
   !      raw that over-states a stiff exchange that has already settled -- after rain, the 0.4 s   !
   !      soil-to-wood exchange -- by a factor of its stiffness, so it is filtered twice through    !
   !      (I - h dF/dW)^-1, which leaves slow modes alone and divides a stiff one's error by its    !
   !      stiffness squared: what step doubling would measure, at the price of two back-            !
   !      substitutions with the last Newton factorisation (docs/science/soil_plant_water.md §3).   !
   !      Each store's error is weighed against atol of its scale plus rtol of its water. ----------!
   real(wp) function substep_error(spac, h, ev0, ev1, scale_soil, scale_coh) result(err)
      type(spac_solver_t), intent(in) :: spac
      real(wp),            intent(in) :: h
      type(spac_eval_t),   intent(in) :: ev0, ev1   !< the network at the sub-step's start and end
      real(wp),            intent(in) :: scale_soil(:), scale_coh(:,:)
      real(wp) :: v_soil(spac%nsoil), v_coh(COHORT_STORES, spac%ncoh)
      real(wp) :: y_soil(spac%nsoil), y_coh(COHORT_STORES, spac%ncoh)
      integer(ik) :: pass
      v_soil = h * (ev1%inflow_soil - ev0%inflow_soil)
      v_coh  = h * (ev1%inflow_coh - ev0%inflow_coh)
      where (spac%soil_fixed) v_soil = 0.0_wp
      do pass = 1, 2
         call arrowhead_solve(spac%sys, v_soil / h, v_coh / h, y_soil, y_coh)
         v_soil = ev1%dwater_soil * y_soil
         v_coh  = ev1%dwater_coh * y_coh
         where (spac%soil_fixed) v_soil = 0.0_wp
      end do
      err = maxval(0.5_wp * abs(v_soil) / (spac%atol * scale_soil + spac%rtol * abs(ev1%water_soil)))
      if (spac%ncoh > 0_ik)                                                                                     &
         err = max(err, maxval(0.5_wp * abs(v_coh) / (spac%atol * scale_coh + spac%rtol * abs(ev1%water_coh))))
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
   subroutine commit_substep(spac, h, tau_soil, tau_coh, ev, spilled)
      type(spac_solver_t), intent(inout) :: spac
      real(wp),            intent(in)    :: h, tau_soil(:), tau_coh(:,:)
      type(spac_eval_t),   intent(in)    :: ev
      logical,             intent(out)   :: spilled   !< the pond spilled, so its unknown moved
      real(wp) :: w_mix, e_mix, t_mix, u_mix, amount, spill, t_pond, pond_flux, sap(spac%ncoh)
      integer(ik) :: c, i

      !----- The flows that close the pond's and each leaf's balance at their curves' water. -------!
      pond_flux = (spac%water_soil(POND) + h * spac%rain - ev%water_soil(POND)) / h
      do i = 1_ik, spac%ncoh
         sap(i) = (ev%water_coh(LEAF, i) - spac%water_coh(LEAF, i)) / h + ev%transp(i)
      end do

      !----- The pond's enthalpy, in the order its water moves: the rain joins it first, so what    !
      !      infiltrates leaves at the mixture's temperature; water the top layer pushes up arrives   !
      !      at that layer's temperature; then the spill leaves at the pond's own. --------------!
      w_mix = spac%water_soil(POND) + h * spac%rain
      e_mix = spac%e_pond + h * spac%rain * internal_energy_liquid(spac%t_rain)
      spac%rain_enth = spac%rain_enth + h * spac%rain * internal_energy_liquid(spac%t_rain)
      t_mix = spac%t_rain
      if (w_mix > POND_TINY) t_mix = liquid_energy_to_temp(e_mix / w_mix)
      u_mix = internal_energy_liquid(t_mix)
      if (pond_flux >= 0.0_wp) then
         amount = h * pond_flux
         spac%infil = spac%infil + amount ; spac%infil_enth = spac%infil_enth + amount * u_mix
         spac%e_pond = e_mix - amount * u_mix
      else
         amount = -h * pond_flux
         spac%exfil = spac%exfil + amount
         spac%exfil_enth = spac%exfil_enth + amount * internal_energy_liquid(spac%t_soil(1))
         spac%e_pond = e_mix + amount * internal_energy_liquid(spac%t_soil(1))
      end if

      !----- Water, by the flows; the pond and the leaves at their curves, the top layer and the woods !
      !      across the paths that closed them. -----------------------------------------------!
      spac%water_soil = spac%water_soil + h * ev%inflow_soil
      spac%water_coh  = spac%water_coh + h * ev%inflow_coh
      spac%water_soil(POND) = ev%water_soil(POND)
      spac%water_soil(2) = spac%water_soil(2) + h * (pond_flux - ev%pond_flux)
      do i = 1_ik, spac%ncoh
         spac%water_coh(LEAF, i) = ev%water_coh(LEAF, i)
         spac%water_coh(WOOD, i) = spac%water_coh(WOOD, i) - h * (sap(i) - ev%sapflow(i))
      end do
      spac%tau_soil = tau_soil ; spac%tau_coh = tau_coh
      do c = 1_ik, spac%nsoil
         if (spac%soil_fixed(c)) spac%water_soil(c) = ev%water_soil(c)
      end do

      !----- Spill: water above the pond's capacity runs off at the pond's temperature. --------!
      spill = max(0.0_wp, spac%water_soil(POND) - spac%w_pond_max)
      spilled = spill > 0.0_wp
      if (spilled) then
         t_pond = liquid_energy_to_temp(spac%e_pond / spac%water_soil(POND))
         spac%spill = spac%spill + spill
         spac%spill_enth = spac%spill_enth + spill * internal_energy_liquid(t_pond)
         spac%e_pond = spac%e_pond - spill * internal_energy_liquid(t_pond)
         spac%water_soil(POND) = spac%w_pond_max
         spac%tau_soil(POND) = pond_water_to_head(spac%w_pond_max)
      end if
      if (spac%water_soil(POND) <= POND_TINY) then
         spac%e_pond = 0.0_wp
      end if

      !----- What each path moved. --------------------------------------------------------!
      spac%face(1:spac%nsl-1) = spac%face(1:spac%nsl-1) + h * ev%face_flux(1:spac%nsl-1)
      spac%drain   = spac%drain   + h * ev%drain
      spac%evap    = spac%evap    + h * ev%evap
      spac%uptake  = spac%uptake  + h * ev%uptake
      spac%sapflow = spac%sapflow + h * sap
      spac%transp  = spac%transp  + h * ev%transp
   end subroutine commit_substep

   !----- Each store's water scale [kg/m2], against which its tolerances are set: a soil layer's or  !
   !      the pond's is a metre of water per metre of the (top) layer, so its tolerances read as     !
   !      water content; a cohort's is the water its wood and leaf cells hold at full turgor, shared  !
   !      by both stores so a leafless cohort's empty leaf is weighed against its plant. ---------!
   subroutine water_scales(spac, hyd, scale_soil, scale_coh)
      type(spac_solver_t),        intent(in)  :: spac
      type(hydro_args_table_t),   intent(in)  :: hyd
      real(wp),                   intent(out) :: scale_soil(:), scale_coh(:,:)
      integer(ik) :: k, i
      real(wp)    :: s
      scale_soil(POND) = rho_h2o * spac%layer(1)%dz
      do k = 1_ik, spac%nsl
         scale_soil(1_ik + k) = rho_h2o * spac%layer(k)%dz
      end do
      do i = 1_ik, spac%ncoh
         associate (lc => hyd%pft(spac%pft(i))%leaf_curve, wc => hyd%pft(spac%pft(i))%wood_curve)
            s = spac%nplant(i) * ((1.0_wp - lc%apoplast_frac) * lc%water_sat * spac%bleaf(i)                    &
                                + (1.0_wp - wc%apoplast_frac) * wc%water_sat * spac%bwood(i))
         end associate
         scale_coh(:, i) = max(s, tiny_num)
      end do
   end subroutine water_scales

   subroutine clear_amounts(spac)
      type(spac_solver_t), intent(inout) :: spac
      spac%face = 0.0_wp ; spac%uptake = 0.0_wp ; spac%sapflow = 0.0_wp ; spac%transp = 0.0_wp
      spac%infil = 0.0_wp ; spac%infil_enth = 0.0_wp ; spac%exfil = 0.0_wp ; spac%exfil_enth = 0.0_wp
      spac%drain = 0.0_wp ; spac%evap = 0.0_wp ; spac%spill = 0.0_wp ; spac%spill_enth = 0.0_wp
      spac%rain_enth = 0.0_wp
      spac%nsub = 0_ik ; spac%nnewton = 0_ik ; spac%nreject = 0_ik ; spac%nfail = 0_ik
   end subroutine clear_amounts

   subroutine eval_allocate(ev, spac)
      type(spac_eval_t),   intent(inout) :: ev
      type(spac_solver_t), intent(in)    :: spac
      integer(ik) :: ns, nc
      ns = spac%nsoil ; nc = spac%ncoh
      if (allocated(ev%psi_soil)) deallocate(ev%psi_soil, ev%dpsi_soil, ev%water_soil, ev%dwater_soil, ev%k_soil, &
                                             ev%dk_soil, ev%inflow_soil, ev%psi_coh, ev%dpsi_coh, ev%water_coh,   &
                                             ev%dwater_coh, ev%inflow_coh, ev%face_flux, ev%uptake, ev%sapflow,  &
                                             ev%transp)
      allocate(ev%psi_soil(ns), ev%dpsi_soil(ns), ev%water_soil(ns), ev%dwater_soil(ns), ev%k_soil(ns),          &
               ev%dk_soil(ns), ev%inflow_soil(ns))
      allocate(ev%psi_coh(COHORT_STORES, nc), ev%dpsi_coh(COHORT_STORES, nc), ev%water_coh(COHORT_STORES, nc),    &
               ev%dwater_coh(COHORT_STORES, nc), ev%inflow_coh(COHORT_STORES, nc))
      allocate(ev%face_flux(spac%nsl), ev%uptake(spac%nsl, nc), ev%sapflow(nc), ev%transp(nc))
      ev%face_flux = 0.0_wp
   end subroutine eval_allocate

end module meds_spac_solver
