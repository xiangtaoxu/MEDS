! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_ark -- the default fast-loop integrator, an L-stable two-stage ESDIRK2 with    !
! the ARS(2,2,2) tableau (config name "ark"; its explicit part is empty, so it is not an   !
! IMEX scheme). column_fast_step_ark builds the frozen record, marches one dt_fast,        !
! commits the water network's soil water and pond, and closes the seven conservation      !
! ledgers; ark2_column_step is one step and adaptive_ark_march its embedded-error march   !
! (docs/science/numerical_scheme.md §3 `ark`, eqs. (1)-(4), §5b).                          !
!==========================================================================================!
module meds_fast_ark
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, rho_h2o
   use meds_config, only : meds_config_t
   use meds_site_state_vars, only : site_surface_t
   use meds_fast_control, only : state_wrms_grouped, step_control_factor, build_error_control
   use meds_fast_types, only : fast_dynamic_state_t, fast_static_state_t, fast_forcing_t, fast_output_t,          &
                               column_state_t, column_frozen_t, stage_bflux_t, column_bflux_t,                   &
                               error_control_t, ark_stages_t, fast_workspace_t, WORK_INTEG_STEPS,               &
                               WORK_INTEG_REJECT, WORK_CLAMP_STAGE
   use meds_column_state_ops, only : state_init, state_extrap, state_err_diff, state_sub, bflux_zero,             &
                                     bflux_add, bflux_bweight, clamp_cas, soil_water_store, soil_energy_store,    &
                                     plant_water_store, canopy_film_store, deposit_condensate,                    &
                                     clamp_canopy_film, unpack_column_state, diagnose_soil_temps,                 &
                                     transp_per_leaf_area, commit_tissue_energy
   use meds_fast_be_stage, only : column_be_stage, advance_plant_water, advance_surf_water_full, water_step
   use meds_fast_frozen, only : build_column_frozen, record_water_flows
   use meds_diag_vars,       only : CO_TRANSP
   use meds_therm_lib, only : internal_energy_liquid
   use meds_budget_check, only : budget_check, BUDGET_ENERGY_RATE_FLOOR, BUDGET_WATER_RATE_FLOOR, BUDGET_CO2_RATE_FLOOR
   implicit none
   private

   public :: column_fast_step_ark
   public :: ark2_column_step, adaptive_ark_march   ! used only by the tests

contains

   !---------------------------------------------------------------------------------------!
   ! ark2_column_step -- one ESDIRK2 step, docs/science/numerical_scheme.md eqs. (1)-(2): the  !
   ! ARS(2,2,2) tableau (Ascher, Ruuth & Spiteri 1997), two backward-Euler stages over gamma*dt !
   ! (column_be_stage), stiffly accurate, so Y3 is the step. Soil water moves at the network's  !
   ! rate in the stages; plant water and the films are advanced once over dt on the stages'    !
   ! b-weighted fluxes, so neither adds to the embedded error y_err. --------------------------!
   subroutine ark2_column_step(y, frozen, site, ec, n, nsl, dt, y_out, y_err, stages, bf, clamp_n)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      type(site_surface_t),  intent(in)  :: site
      type(error_control_t), intent(in)  :: ec     !< its budget stop (column_be_stage)
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(inout) :: y_out, y_err   !< filled in place
      type(ark_stages_t),    intent(inout) :: stages    !< this thread's stage storage, filled in place
      type(column_bflux_t), optional, intent(out) :: bf   !< b-weighted boundary-flux AMOUNTS over dt (ledger)
      !----- Stage-clamp count (the report's WORK_CLAMP_STAGE): only the stage-3 base is clamped. ---!
      integer(ik), optional, intent(inout) :: clamp_n
      real(wp), parameter :: GAMMA = 0.2928932188134524_wp   ! 1 - 1/sqrt(2)
      real(wp), parameter :: BETA  = 2.4142135623730951_wp   ! (1-gamma)/gamma = 1 + sqrt(2)
      type(stage_bflux_t)   :: bf2, bf3
      real(wp)              :: transp_bw(n)
      real(wp)              :: film_evap_leaf_bw(n), film_evap_wood_bw(n)
      associate (Y2 => stages%y2, base3 => stages%base3, Y3 => stages%y3, sf2 => stages%sf2, sf3 => stages%sf3)
         !----- Stage 2: a gamma*dt backward-Euler stage from y_n. --------------------------------------!
         call column_be_stage(y, frozen, site, ec, n, nsl, GAMMA*dt, Y2, sf2, bf=bf2)
         !----- Stage 3 from the extrapolated base. Extrapolating by BETA can carry the canopy-air      !
         !      enthalpy to a temperature where qsat overflows, so the base's canopy air is clamped to a  !
         !      physical range: an identity in range, and an oversized step is rejected anyway. Soil     !
         !      water needs no clamp: it moves at a steady rate, so the base stays between its ends. ----!
         call state_extrap(y, BETA, Y2, n, nsl, base3)
         call clamp_cas(base3, nfire=clamp_n)
         call column_be_stage(base3, frozen, site, ec, n, nsl, GAMMA*dt, Y3, sf3, bf=bf3)
         call state_init(Y3, n, nsl, y_out)
         !----- Plant water over the full dt at the network's rates, debited the same b-weighted         !
         !      transpiration the canopy air's vapour balance was credited. -----------------------------!
         transp_bw(1:n) = (1.0_wp - GAMMA)*sf2%transp_c(1:n) + GAMMA*sf3%transp_c(1:n)
         call advance_plant_water(y, frozen, n, dt, transp_bw, y_out)
         !----- The films the same way (a no-op with canopy_water_on off). ------------------------------!
         film_evap_leaf_bw(1:n) = (1.0_wp - GAMMA)*sf2%film_evap_leaf(1:n) + GAMMA*sf3%film_evap_leaf(1:n)
         film_evap_wood_bw(1:n) = (1.0_wp - GAMMA)*sf2%film_evap_wood(1:n) + GAMMA*sf3%film_evap_wood(1:n)
         call advance_surf_water_full(y, frozen, n, dt, film_evap_leaf_bw, film_evap_wood_bw, y_out)
         !----- Embedded first-order error, eq. (2). ------------------------------------------------------!
         call state_err_diff(Y3, base3, Y2, y, n, nsl, y_err)
         !----- b-weighted boundary-flux amounts, b = (0, 1-gamma, gamma), which telescope exactly. -----!
         if (present(bf)) then
            call bflux_bweight(bf, bf2, bf3, dt, GAMMA)
            !----- ...the tissue-temperature time integrals the tissues' energy is committed from. -----!
            allocate(bf%tissue_leaf_int(n), bf%tissue_wood_int(n))
            bf%tissue_leaf_int(1:n) = dt * ((1.0_wp - GAMMA)*sf2%leaf_temp(1:n) + GAMMA*sf3%leaf_temp(1:n))
            bf%tissue_wood_int(1:n) = dt * ((1.0_wp - GAMMA)*sf2%wood_temp(1:n) + GAMMA*sf3%wood_temp(1:n))
            !----- ...and each cohort's transpiration, the transp_bw the leaf water was debited. -------!
            allocate(bf%transp_coh(n))
            bf%transp_coh(1:n) = dt * transp_bw(1:n)
         end if
      end associate
   end subroutine ark2_column_step

   !---------------------------------------------------------------------------------------!
   ! adaptive_ark_march -- march to t_end, accepting a step when the grouped norm of its       !
   ! embedded error is within tolerance (docs/science/numerical_scheme.md eq. (3)) and sizing  !
   ! the next by step_control_factor (eq. (4), p = 1). Returns the step and reject counts.     !
   !---------------------------------------------------------------------------------------!
   subroutine adaptive_ark_march(y0, frozen, site, n, nsl, t_end, ec, dt_init, y_out, work, nsteps, nrej, acc,        &
                                 dt_warm_out, clamp_n)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: frozen
      type(site_surface_t),  intent(in)  :: site
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: t_end, dt_init
      type(error_control_t), intent(in)  :: ec       !< tolerances and controller (meds_fast_control)
      type(column_state_t),  intent(inout) :: y_out  !< filled in place
      type(fast_workspace_t), intent(inout) :: work  !< this thread's march and stage storage
      integer(ik),           intent(out) :: nsteps, nrej
      !----- Warm start (docs/science/numerical_scheme.md §3 `ark`): the last accepted step that the     !
      !      interval's end did not cut short, to seed the next call. -------------------------------!
      real(wp),    optional, intent(out) :: dt_warm_out
      type(column_bflux_t), optional, intent(out) :: acc  !< accumulated boundary-flux amounts (ledger)
      !----- Stage-clamp count over the march, rejected trials included: a clamp still says the        !
      !      controller probed too far. Accumulated, so the caller zeroes it. ------------------------!
      integer(ik), optional, intent(inout) :: clamp_n

      type(column_bflux_t) :: bfsub
      real(wp)             :: t, dt, err, err_prev, fac, dt_floor
      real(wp)             :: dt_try, dt_warm
      logical              :: clamped
      associate (y => work%y_march, y_new => work%y_new, y_err => work%y_err, y_lo => work%y_lo)
         if (present(acc)) call bflux_zero(acc, n)
         !----- Sub-step floor t_end/64, at least 1e-2 s. The backward-Euler stages are L-stable, so a  !
         !      floor step stays bounded even when its error is above tolerance; the floor caps a stiff  !
         !      transient at 64 sub-steps and surfaces a non-finite state promptly. ----------------------!
         dt_floor = max(1.0e-2_wp, t_end / 64.0_wp)   ! clamp-ok: the smallest sub-step [s] (above)

         call state_init(y0, n, nsl, y)
         t = 0.0_wp
         dt = min(dt_init, t_end)
         nsteps = 0_ik
         nrej = 0_ik
         err_prev = -1.0_wp                                          ! < 0 => first step uses the I-controller
         dt_warm = dt                                                ! fallback if every step is end-clamped
         do
            if (t >= t_end - tiny_num) exit
            dt_try = dt
            dt = min(dt, t_end - t)
            clamped = dt < dt_try - tiny_num
            call ark2_column_step(y, frozen, site, ec, n, nsl, dt, y_new, y_err, work%stages, bf=bfsub,      &
                                  clamp_n=clamp_n)
            call state_sub(y_new, y_err, n, nsl, y_lo)               ! the 1st-order embedded solution
            !----- Grouped norm over the whole column state, eq. (3). ARK's theta and plant-water terms are  !
            !      structurally zero and dilute it (docs/science/numerical_scheme.md §7, item 3). ----------!
            err = state_wrms_grouped(y_new, y_lo, y, n, nsl, ec%tols)
            !----- A non-finite error means the step was too big (typically the stage-3 base overshot to   !
            !      where qsat overflows): reject it and shrink by fmin, since the NaN poisons the factor.   !
            !      A non-finite floor step is committed, so the driver's state guard reports it. ---------!
            if (err /= err .or. dt /= dt) then
               if (dt <= dt_floor) then
                  call state_init(y_new, n, nsl, y)
                  t = t + dt_floor
                  nsteps = nsteps + 1_ik
                  exit
               end if
               nrej = nrej + 1_ik
               dt = max(dt * ec%fmin, dt_floor)
               cycle
            end if
            fac = step_control_factor(err, err_prev, ec)             ! I (default) or PI (Gustafsson) controller
            if (err <= 1.0_wp .or. dt <= dt_floor) then
               call state_init(y_new, n, nsl, y)
               if (present(acc)) call bflux_add(acc, bfsub)          ! accumulate only accepted sub-steps
               t = t + dt
               nsteps = nsteps + 1_ik
               err_prev = err                                        ! remember for the PI controller
               !----- The warm-start seed is the step just accepted, before the growth factor: a grown   !
               !      proposal over-estimates (docs/science/numerical_scheme.md §3 `ark`). ---------------!
               if (.not. clamped) dt_warm = dt
               dt = dt * fac
            else
               nrej = nrej + 1_ik
               dt = dt * fac
            end if
            if (nsteps + nrej > 4096_ik) exit                        ! hard backstop (should never trigger)
         end do
         call state_init(y, n, nsl, y_out)
         if (present(dt_warm_out)) dt_warm_out = dt_warm
      end associate
   end subroutine adaptive_ark_march
   !=======================================================================================!
   !  column_fast_step_ark -- the INTEG_ARK fast step (docs/science/numerical_scheme.md §3 `ark`):  !
   !  build the frozen record and step the water network into it, advance one dt_fast, commit the   !
   !  network's water and the tissue energy, deposit the condensate, and close the seven ledgers.   !
   !=======================================================================================!
   subroutine column_fast_step_ark(dt_fast, cfg, site, static, forcing, dynamic, out, work)
      real(wp),                   intent(in)    :: dt_fast
      type(meds_config_t),        intent(in)    :: cfg
      type(site_surface_t),       intent(in)    :: site
      type(fast_static_state_t),  intent(in)    :: static
      type(fast_forcing_t),       intent(in)    :: forcing
      type(fast_dynamic_state_t), intent(inout) :: dynamic
      type(fast_output_t),        intent(inout) :: out
      type(fast_workspace_t),     intent(inout) :: work    !< this thread's reusable step storage

      type(column_frozen_t)  :: frozen
      type(column_state_t)   :: y, y_out, ycur, ytmp, yerr
      type(column_bflux_t)   :: acc, bfsub
      real(wp)    :: dt0, cas_mass_capacity, cas_molar_capacity, enth0, shv0, co20, enth1, shv1, co21, e_soil0, e_soil1, &
           w_soil0, w_soil1
      real(wp)    :: e_pond0, e_pond1   !< [J/m2] pond enthalpy, start and end of the step
      real(wp)    :: w_plant0, w_plant1
      real(wp)    :: w_surface0, dt_warm_next
      !----- Canopy-film ledger scratch. --------------------------------------------------------!
      real(wp)    :: surf_water0, surf_water1, surf_enth0, surf_enth1
      real(wp)    :: surf_overflow, surf_deficit, intercept_total
      !----- Condensate deposited into soil layer 1 [kg/m2] and its enthalpy [J/m2], kept so the     !
      !      ledgers book the same numbers the state update applied. --------------------------------!
      real(wp)    :: cond_dep_mass, cond_dep_enth
      real(wp)    :: tissue_store0, tissue_store1
      real(wp)    :: t_leaf_avg(static%n), t_wood_avg(static%n)
      type(error_control_t) :: ec
      integer(ik) :: n, nsl, isub, nsub, nsteps, nrej
      logical     :: halt_budgets     !< hard-stop on a non-closing budget ([fast].debug_error)

      n = static%n
      nsl = site%soil_profile%n_active
      !----- The stage-clamp count accumulates down the call chain, so this step's tally starts here. !
      !      ARK has no explicit stage, so theta_ood_max stays 0. ------------------------------------!
      out%report%work(WORK_CLAMP_STAGE) = 0_ik
      out%report%theta_ood_max = 0.0_wp

      !----- The run's error control: tolerances and controller, eqs. (3)-(4), and the budget stops. --!
      ec = build_error_control(cfg)
      halt_budgets = ec%debug_error

      !----- State^n pond mass, for the whole-water ledger. ------------------------------------------!
      w_surface0 = dynamic%soil_water%w_surface

      call build_column_frozen(dt_fast, cfg, site, static, forcing, dynamic, out, work, frozen, y)
      call water_step(dt_fast, static, site%soil_profile, site%soil, cfg%pft%hydraulics, frozen%water_drive,     &
                      ec%tols, work%water, out%report%work)
      call record_water_flows(dt_fast, cfg, site, static, dynamic, work%water, out, frozen)

      !----- State^n pond enthalpy, taken before the pond commit below overwrites it                  !
      !      (build_column_frozen writes nothing to dynamic%soil_water). -----------------------------!
      e_pond0 = dynamic%soil_water%w_surface_energy

      !----- Advance one dt_fast: adaptive, or fixed equal sub-steps. ---------------------------------!
      if (cfg%ark_adaptive) then
         !----- Warm start (docs/science/numerical_scheme.md §3 `ark`): seed from the step this patch    !
         !      settled on last dt_fast; a first call starts from dt_fast. ------------------------------!
         dt0 = dt_fast
         if (dynamic%adapt_dt_last > tiny_num) dt0 = min(dynamic%adapt_dt_last, dt_fast)
         call adaptive_ark_march(y, frozen, site, n, nsl, dt_fast, ec, dt0, y_out, work, nsteps, nrej,         &
                                 acc=acc, dt_warm_out=dt_warm_next,                                 &
                                 clamp_n=out%report%work(WORK_CLAMP_STAGE))
         dynamic%adapt_dt_last = dt_warm_next
      else
         nsub = max(1_ik, cfg%ark_fixed_substep)
         nrej = 0_ik
         ycur = y
         call bflux_zero(acc, n)
         do isub = 1_ik, nsub
            call ark2_column_step(ycur, frozen, site, ec, n, nsl, dt_fast/real(nsub, wp), ytmp, yerr, work%stages, &
                                  bf=bfsub, clamp_n=out%report%work(WORK_CLAMP_STAGE))
            call bflux_add(acc, bfsub)
            ycur = ytmp
         end do
         y_out = ycur
         nsteps = nsub
      end if
      !----- Work counters; the water network's own work is counted in water_step. -----------------!
      out%report%work(WORK_INTEG_STEPS) = nsteps
      out%report%work(WORK_INTEG_REJECT) = nrej

      !----- Soil water is the water network's: the stages moved theta at its steady rate, and its   !
      !      end-of-step theta is committed once, here, for the state, the soil temperatures and the  !
      !      ledgers. The pond likewise, through the state vector. ----------------------------------!
      y_out%theta(1:nsl) = frozen%hydrology%theta1(1:nsl)
      y_out%w_surface      = frozen%hydrology%w_surface1
      y_out%w_surface_energy = frozen%hydrology%w_surface_energy1     ! paired with the mass

      !----- Canopy films: capacity clamp and overflow/deficit at the store's boundary, after the      !
      !      march. A film evaporating at its frozen conductance can overdraw below 0, which would     !
      !      create water in the next step's interception, so it is floored at 0 and the deficit is    !
      !      booked as negative outflow. --------------------------------------------------------------!
      surf_overflow = 0.0_wp
      surf_deficit = 0.0_wp
      if (cfg%canopy_water_on) then
         call clamp_canopy_film(y_out, static%lai, static%wai,                                  &
                                cfg%pft%leaf_surf_water_max(static%pft(1:n)),                    &
                                cfg%pft%wood_surf_water_max(static%pft(1:n)), n, surf_overflow, surf_deficit)
      end if

      !----- Unpack into dynamic and diagnose the soil temperatures. -------------------------------!
      call unpack_column_state(y_out, n, nsl, dynamic)
      dynamic%soil_water%w_surface      = y_out%w_surface
      dynamic%soil_water%w_surface_energy = y_out%w_surface_energy
      call diagnose_soil_temps(y_out, site%soil_profile%dry_heat_capacity, nsl, dynamic%soil_energy%soil_temp, &
            dynamic%soil_energy%soil_fliq)
      !----- The tissue heat store is solved in the stages, jointly with the canopy air               !
      !      (docs/science/vegetation_energy_dynamics.md §One formulation). Each tissue's energy is     !
      !      committed from its time-average temperature over the march, the complement of what the   !
      !      canopy air received (docs/science/numerical_scheme.md §5b). ------------------------------!
      t_leaf_avg(1:n) = acc%tissue_leaf_int(1:n) / dt_fast
      t_wood_avg(1:n) = acc%tissue_wood_int(1:n) / dt_fast
      call commit_tissue_energy(frozen%tissue, static%nplant, dt_fast, t_leaf_avg, t_wood_avg, n, dynamic, &
                                tissue_store0, tissue_store1)
      !----- Each cohort's transpiration over the step, as the leaf water was debited. -------------!
      if (out%capture .and. allocated(acc%transp_coh))                                              &
         out%cdiag(CO_TRANSP, 1:n) = transp_per_leaf_area(acc%transp_coh(1:n), dt_fast, static%nplant(1:n), &
                                                          static%leaf_area(1:n))

      !----- The conservation ledgers, from the b-weighted boundary-flux amounts accumulated over the  !
      !      sub-steps (acc), close to round-off (docs/science/numerical_scheme.md §5b). First deposit  !
      !      the condensate, as on RK45: dew and fog land inside the column, so it goes, mass and      !
      !      enthalpy paired, into soil layer 1 (§7, item 4), before the store totals below. ----------!
      cond_dep_mass = 0.0_wp
      cond_dep_enth = 0.0_wp
      if (acc%whole_cond > 0.0_wp) then
         cond_dep_mass = acc%whole_cond
         cond_dep_enth = acc%whole_cond_enth       ! b-weighted at the stage CAS temperatures (see column_bflux_t)
         call deposit_condensate(y_out, site%soil_profile%dz(1), cond_dep_mass, cond_dep_enth)
         dynamic%soil_water%theta(1)       = y_out%theta(1)
         dynamic%soil_energy%soil_energy(1) = y_out%soil_energy(1)
      end if

      !----- The wood's absorbed radiation is already inside coh_rnet and acc%whole_enth_in. -------!

      cas_mass_capacity = frozen%cas%box%air_mass_capacity
      cas_molar_capacity = frozen%cas%box%air_molar_capacity
      enth0 = y%cas_enthalpy
      shv0 = y%cas_shv
      co20 = y%cas_co2
      enth1 = y_out%cas_enthalpy
      shv1 = y_out%cas_shv
      co21 = y_out%cas_co2
      e_pond1 = frozen%hydrology%w_surface_energy1
      e_soil0 = soil_energy_store(y%soil_energy,     site%soil_profile%dz, nsl)
      e_soil1 = soil_energy_store(y_out%soil_energy, site%soil_profile%dz, nsl)
      w_soil0 = soil_water_store(y%theta,     site%soil_profile%dz, nsl)
      w_soil1 = soil_water_store(y_out%theta, site%soil_profile%dz, nsl)
      !----- Plant water is a store: it absorbs the step's mismatch between transpiration and uptake  !
      !      (docs/science/plant_hydraulics.md eqs. (1)-(2)). ----------------------------------------!
      w_plant0 = plant_water_store(static%nplant, y%leaf_water_mass,     y%wood_water_mass,     n)
      w_plant1 = plant_water_store(static%nplant, y_out%leaf_water_mass, y_out%wood_water_mass, n)
      !----- The films, per m2 of ground, valued at the liquid enthalpy the intercepted water arrived  !
      !      with, as surface_derivs charges the tissue (0 with canopy_water_on off). ----------------!
      surf_water0 = canopy_film_store(y%leaf_surf_water,     y%wood_surf_water,     n)
      surf_water1 = canopy_film_store(y_out%leaf_surf_water, y_out%wood_surf_water, n)
      surf_enth0  = surf_water0 * internal_energy_liquid(frozen%hydrology%t_film_valuation)
      surf_enth1  = surf_water1 * internal_energy_liquid(frozen%hydrology%t_film_valuation)
      intercept_total = sum(frozen%film%intercept_leaf(1:n) + frozen%film%intercept_wood(1:n))
      !----- [fast].debug_error turns a non-closing budget into an error stop; the tolerance scales   !
      !      with the boundary flux, not the store (docs/science/numerical_scheme.md §5b). ------------!
      out%report%exports%atm_heat_export = acc%atm_heat_out
      out%report%exports%atm_vap_export  = acc%atm_vap_out
      call budget_check(out%report%ledgers%cas_energy,                                                    &
                        cas_mass_capacity*enth0, cas_mass_capacity*enth1, acc%cas_enth_in, acc%cas_enth_out, &
                        dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'cas_energy (ark)', halt_budgets)
      call budget_check(out%report%ledgers%cas_water,                                                     &
                        cas_mass_capacity*shv0,  cas_mass_capacity*shv1,  acc%cas_vap_in,  acc%cas_vap_out,  &
                        dt_fast, BUDGET_WATER_RATE_FLOOR, 'cas_water (ark)', halt_budgets)
      call budget_check(out%report%ledgers%cas_co2,                                                       &
                        cas_molar_capacity*co20, cas_molar_capacity*co21, acc%cas_co2_in,  acc%cas_co2_out,  &
                        dt_fast, BUDGET_CO2_RATE_FLOOR, 'cas_co2 (ark)', halt_budgets)
      !----- The condensate deposit is an inflow to the soil budgets, though the whole-column ledger   !
      !      sees it as a transfer between two of its stores. RK45 keeps no per-store budgets: every   !
      !      store rides the same right-hand side, so its whole-column ledgers are the per-store ones. -!
      call budget_check(out%report%ledgers%soil_energy, e_soil0, e_soil1, acc%soil_enth_in + cond_dep_enth,            &
                        acc%soil_enth_out, dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'soil_energy (ark)', halt_budgets)
      !----- Soil water: theta^n -> theta1, in through the top face net of evaporation and            !
      !      exfiltration (q_top), out through the bottom and to the roots, plus the condensate. ------!
      call budget_check(out%report%ledgers%soil_water,  w_soil0, w_soil1,                                              &
                        frozen%hydrology%q_top*rho_h2o*dt_fast + cond_dep_mass,                            &
                        (frozen%hydrology%drainage + frozen%roots%uptake)*dt_fast,                         &
                        dt_fast, BUDGET_WATER_RATE_FLOOR, 'soil_water (ark)', halt_budgets)
      !----- Whole water: rainfall, snowfall and shed water in; the CAS vapour export, drainage and    !
      !      runoff out; soil, pond, plant water, films, snow pack and canopy-air vapour stored.       !
      !      Sublimation leaves through the vapour export and meltwater moves pack -> pond, so        !
      !      neither needs a term. -------------------------------------------------------------------!
      call budget_check(out%report%ledgers%whole_water,                                                                &
                        w_soil0 + cas_mass_capacity*shv0 + w_surface0 + w_plant0 + surf_water0 + frozen%snow%swe0,  &
                        w_soil1 + cas_mass_capacity*shv1 + frozen%hydrology%w_surface1 + w_plant1 + surf_water1 &
                        + frozen%snow%swe1, &
                        acc%whole_wat_in + (forcing%rainfall + forcing%snowfall + static%shed_water_rate)*dt_fast,       &
                        acc%whole_wat_out                                                                   &
                        + (frozen%hydrology%runoff_surf + frozen%hydrology%drainage)*dt_fast                &
                                          + surf_overflow - surf_deficit,                                 &
                        dt_fast, BUDGET_WATER_RATE_FLOOR, 'whole_water (ark)', halt_budgets)
      call budget_check(out%report%ledgers%whole_energy,                                                               &
                             !----- Meltwater passes pack -> pond, both stores, so it telescopes. ----------!
                        e_soil0                           + cas_mass_capacity*enth0 + surf_enth0                     &
                        + frozen%snow%enth0 + e_pond0 + tissue_store0,                               &
                        e_soil1 + cas_mass_capacity*enth1 + surf_enth1 + frozen%snow%enth1 + e_pond1                &
                        + tissue_store1,                                                                 &
                        acc%whole_enth_in + intercept_total*dt_fast*internal_energy_liquid(frozen%hydrology%t_film_valuation) &
                                          + frozen%snow%acc_enth                                       &
                                          + frozen%hydrology%ground_inflow_enth*dt_fast,                              &
                        acc%whole_enth_out + (surf_overflow - surf_deficit)                                 &
                        * internal_energy_liquid(frozen%hydrology%t_film_valuation) &
                                           + frozen%hydrology%runoff_enth*dt_fast,                                    &
                        dt_fast, BUDGET_ENERGY_RATE_FLOOR, 'whole_energy (ark)', halt_budgets)
   end subroutine column_fast_step_ark

end module meds_fast_ark
