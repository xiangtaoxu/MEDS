! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_column_hydrology -- the soil column's water: the retention curves, canopy           !
! interception, the soil's boundary processes, and the column itself, solved as the        !
! soil-plant water network with no cohorts (docs/science/soil_plant_water.md §1; the       !
! processes are docs/science/soil_biophysics.md §The soil water column, §Boundary          !
! conditions and §Canopy interception).                                                    !
!                                                                                          !
!   1. CONSTITUTIVE  : van Genuchten and Campbell round trips and moisture capacity.           !
!   2. INTERCEPTION  : the canopy film's capture, capacity and drip.                            !
!   3. BEDROCK       : a sealed column keeps its water and drains nothing.                       !
!   4. FREE DRAIN    : a wet column drains, and its water closes against the drainage.          !
!   4b. LEAKY        : drainage is f K(theta_N): f = 1 is free drainage and f = 0 bedrock, bit   !
!                      for bit; the flux scales by f; a leaky column drains less than a free     !
!                      one, and its water closes against the drainage.                           !
!   5. DOWNPOUR      : rain beyond a wet clay's infiltration ponds, and runs off past the pond.  !
!   6. REWETTING     : a dried top layer re-wets within 12 h of rain, with little runoff.        !
!   7. SUB-FREEZING  : rain below the triple point: the pond's enthalpy closes, and what         !
!                      infiltrates carries its own enthalpy.                                     !
!   8. AQUIFER       : the head-driven bottom is two-way and settles to hydrostatic equilibrium. !
!   9. EVAPORATION   : area-weighted by the snow-free share; rises monotonically with moisture.   !
!==========================================================================================!
program test_column_hydrology
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, grav_head, n_soil_layer_max
   use meds_soil_types, only : chydro_forcing_t
   use meds_soil_params,      only : soil_profile_t, soil_params_t, SOIL_BC_DRAIN, SOIL_BC_AQUIFER
   use meds_config,             only : build_soil_hydr_params
   use meds_test_soil_grid, only : geometric_layer_bottoms
   use meds_soil_hydr_lib,   only : SOIL_RETENTION_VG, SOIL_RETENTION_CAMPBELL, soil_head_to_theta,  &
                                    soil_theta_to_head
   use meds_soil_water,       only : ground_evaporation_rate, bottom_flux
   use meds_spac_solver,    only : spac_solver_t, spac_allocate, spac_set_state,     &
                                     spac_advance, POND
   use meds_plant_types,      only : hydro_args_table_t
   use meds_therm_lib,        only : internal_energy_liquid, liquid_energy_to_temp
   use meds_plant_biophysics, only : intercept_canopy_layer
   implicit none

   type(hydro_args_table_t)   :: no_plants

   allocate(no_plants%pft(1))
   call test_constitutive()
   call test_interception()
   call test_bedrock_conserve()
   call test_free_drain()
   call test_leaky_drain()
   call test_downpour()
   call test_dried_surface_rewets()
   call test_pond_subfreezing_inflow()
   call test_aquifer_head_bc()
   call test_snow_free_evap()
   call test_evap_moisture_response()
   call test_report('test_column_hydrology')

contains

   !----- A 10-layer loam column (2 m). --------------------------------------------------------!
   subroutine loam_params(params)
      type(soil_profile_t), intent(out) :: params
      real(wp) :: bottoms(n_soil_layer_max)
      call geometric_layer_bottoms(2.0_wp, 10_ik, 3.0_wp, bottoms)
      call build_soil_hydr_params(10_ik, SOIL_RETENTION_VG, bottoms, 0.43_wp, 0.078_wp,                 &
           2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, params)
   end subroutine loam_params

   !----- The soil-only network for a column at theta with a pond of w_pond, and air that keeps    !
   !      the ground from evaporating unless the caller says otherwise. ---------------------------!
   subroutine soil_network(spac, params, opts, theta, w_pond, e_pond)
      type(spac_solver_t), intent(inout) :: spac
      type(soil_profile_t),  intent(in)    :: params
      type(soil_params_t),   intent(in)    :: opts
      real(wp),              intent(in)    :: theta(:), w_pond, e_pond
      real(wp) :: dum(0)
      call spac_allocate(spac, params%n_active, 0_ik)
      !----- [fast].rtol and water_atol at their defaults. -----------------------------------!
      spac%rtol = 1.0e-3_wp ; spac%atol = 1.0e-4_wp ; spac%w_pond_max = opts%w_pond_max
      spac%t_soil = 290.0_wp ; spac%t_rain = 290.0_wp
      call spac_set_state(spac, params, no_plants, theta, w_pond, e_pond, dum, dum)
   end subroutine soil_network

   real(wp) function column_water(spac)
      type(spac_solver_t), intent(in) :: spac
      column_water = sum(spac%water_soil)
   end function column_water

   function humid_air() result(f)
      type(chydro_forcing_t) :: f
      f%t_ground = 290.0_wp ; f%q_air = 0.05_wp ; f%rho_air = 1.2_wp ; f%r_aero = 100.0_wp   ! no evaporation
   end function humid_air

   !=======================================================================================!
   subroutine test_constitutive()
      real(wp) :: th, ps, th_rt
      print '(a)', 'test_constitutive:'
      !----- van Genuchten loam round-trip. -----!
      th    = 0.30_wp
      ps    = soil_theta_to_head(SOIL_RETENTION_VG, th, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp)
      th_rt = soil_head_to_theta(SOIL_RETENTION_VG, ps, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp)
      call check('vG theta round-trip', th_rt, th, 1.0e-9_wp)
      !----- Campbell loam round-trip. -----!
      th    = 0.30_wp
      ps    = soil_theta_to_head(SOIL_RETENTION_CAMPBELL, th, 0.44_wp, 0.0_wp, -0.26_wp, 5.65_wp)
      th_rt = soil_head_to_theta(SOIL_RETENTION_CAMPBELL, ps, 0.44_wp, 0.0_wp, -0.26_wp, 5.65_wp)
      call check('Campbell theta round-trip', th_rt, th, 1.0e-9_wp)
   end subroutine test_constitutive

   !=======================================================================================!
   subroutine test_interception()
      real(wp) :: lw, tf, dr, sw, bal, rain, dt
      print '(a)', 'test_interception:'
      dt = 600.0_wp ; rain = 1.0e-4_wp
      !----- Dry canopy fills below capacity: no drip, exact balance. -----!
      lw = 0.0_wp
      call intercept_canopy_layer(lw, rain, 2.0_wp, 0.5_wp, 0.0_wp, dt, 0.1_wp, 0.1_wp, 0.5_wp, 1.0_wp,&
                                  tf, dr, sw)
      bal = tf + lw / dt + 0.0_wp                       ! throughfall + storage-rate + evap
      call check('interception water balance', bal, rain, 1.0e-12_wp)
      call check_true('no drip below capacity', dr < 1.0e-30_wp, dr)
      call check_true('leaf_water within capacity', lw <= 0.1_wp * 2.5_wp + 1.0e-12_wp, lw)
      !----- Saturated canopy overflows: drip appears, storage capped. -----!
      lw = 0.25_wp                                       ! = 0.1*LAI + 0.1*SAI (full)
      call intercept_canopy_layer(lw, rain, 2.0_wp, 0.5_wp, 0.0_wp, dt, 0.1_wp, 0.1_wp, 0.5_wp, 1.0_wp,&
                                  tf, dr, sw)
      call check_true('drip when full', dr > 0.0_wp, dr)
      call check_true('capacity respected', lw <= 0.25_wp + 1.0e-12_wp, lw)
      call check_true('sigma_w saturates to 1', abs(sw - 1.0_wp) < 1.0e-9_wp, sw)
      !----- A film already ABOVE capacity (leaf area lost under a full film) drips the excess:       !
      !      storage capped, and throughfall + storage change = rain -- nothing discarded. ---------!
      lw = 0.30_wp                                       ! capacity is 0.25: 0.05 too much
      call intercept_canopy_layer(lw, rain, 2.0_wp, 0.5_wp, 0.0_wp, dt, 0.1_wp, 0.1_wp, 0.5_wp, 1.0_wp,&
                                  tf, dr, sw)
      call check_true('over capacity: storage capped', lw <= 0.25_wp + 1.0e-12_wp, lw)
      bal = tf + (lw - 0.30_wp) / dt
      call check('over capacity: the excess drips, the water balance holds', bal, rain, 1.0e-12_wp)
   end subroutine test_interception

   !=======================================================================================!
   subroutine test_bedrock_conserve()
      type(soil_profile_t)  :: params
      type(soil_params_t)   :: opts
      type(spac_solver_t) :: spac
      real(wp) :: w_before, theta(10)
      integer(ik) :: step
      print '(a)', 'test_bedrock_conserve:'
      call loam_params(params) ; theta = 0.30_wp ; theta(1:3) = 0.40_wp
      opts%water_bottom_bc = SOIL_BC_DRAIN ; opts%drain_factor = 0.0_wp       ! sealed bedrock
      call soil_network(spac, params, opts, theta, 0.0_wp, 0.0_wp)
      w_before = column_water(spac)
      do step = 1_ik, 10_ik
         call spac_advance(spac, params, opts, humid_air(), no_plants, 600.0_wp)
      end do
      call check('sealed column conserves water', column_water(spac), w_before, 1.0e-10_wp * w_before)
      call check_true('bedrock drainage = 0', abs(spac%drain) < 1.0e-30_wp, spac%drain)
   end subroutine test_bedrock_conserve

   !=======================================================================================!
   subroutine test_free_drain()
      type(soil_profile_t)  :: params
      type(soil_params_t)   :: opts
      type(spac_solver_t) :: spac
      real(wp) :: w0, w1, drained, theta(10)
      integer(ik) :: step
      print '(a)', 'test_free_drain:'
      call loam_params(params) ; theta = 0.40_wp
      opts%water_bottom_bc = SOIL_BC_DRAIN ; opts%drain_factor = 1.0_wp       ! free drainage
      call soil_network(spac, params, opts, theta, 0.0_wp, 0.0_wp)
      drained = 0.0_wp ; w0 = column_water(spac)
      do step = 1_ik, 20_ik
         call spac_advance(spac, params, opts, humid_air(), no_plants, 600.0_wp)
         drained = drained + spac%drain
      end do
      w1 = column_water(spac)
      call check_true('a wet column drains', drained > 0.0_wp, drained)
      call check('its water closes against the drainage', w0 - w1, drained, 1.0e-10_wp * w0)
   end subroutine test_free_drain

   !=======================================================================================!
   subroutine test_leaky_drain()
      type(soil_profile_t)  :: params
      type(soil_params_t)   :: opts
      type(spac_solver_t) :: spac
      real(wp), parameter :: F_LEAK = 0.05_wp          ! the ED soil maps' factor
      real(wp), parameter :: K_N = 1.234e-7_wp, PSI_N = -0.37_wp   ! [m/s], [m] a deepest layer's state
      real(wp) :: q, dq_dpsi, dq_dk, w0, w1, drained, drained_free, theta(10)
      integer(ik) :: step
      print '(a)', 'test_leaky_drain:'
      call loam_params(params)
      opts%water_bottom_bc = SOIL_BC_DRAIN
      !----- At one state: the two ends of the factor are free drainage and bedrock, exactly. -----!
      opts%drain_factor = 1.0_wp
      call bottom_flux(params, opts, 10_ik, PSI_N, K_N, q, dq_dpsi, dq_dk)
      call check('f = 1 is free drainage: q = K(theta_N), bit for bit', q, K_N, 0.0_wp)
      call check('f = 1: dq/dK = 1 and dq/dpsi = 0', abs(dq_dk - 1.0_wp) + abs(dq_dpsi), 0.0_wp, 0.0_wp)
      opts%drain_factor = 0.0_wp
      call bottom_flux(params, opts, 10_ik, PSI_N, K_N, q, dq_dpsi, dq_dk)
      call check('f = 0 is sealed bedrock: q = 0 and no derivative, bit for bit',                     &
                 abs(q) + abs(dq_dk) + abs(dq_dpsi), 0.0_wp, 0.0_wp)
      opts%drain_factor = F_LEAK
      call bottom_flux(params, opts, 10_ik, PSI_N, K_N, q, dq_dpsi, dq_dk)
      call check('a leaky base drains f K(theta_N)', q, F_LEAK * K_N, 0.0_wp)
      call check('and its flux moves with K(theta_N) by f', dq_dk, F_LEAK, 0.0_wp)
      !----- A wet column, free and leaky. -------------------------------------------------------!
      drained_free = 0.0_wp
      opts%drain_factor = 1.0_wp
      call loam_params(params) ; theta = 0.40_wp
      call soil_network(spac, params, opts, theta, 0.0_wp, 0.0_wp)
      do step = 1_ik, 20_ik
         call spac_advance(spac, params, opts, humid_air(), no_plants, 600.0_wp)
         drained_free = drained_free + spac%drain
      end do
      drained = 0.0_wp
      opts%drain_factor = F_LEAK
      call soil_network(spac, params, opts, theta, 0.0_wp, 0.0_wp)
      w0 = column_water(spac)
      do step = 1_ik, 20_ik
         call spac_advance(spac, params, opts, humid_air(), no_plants, 600.0_wp)
         drained = drained + spac%drain
      end do
      w1 = column_water(spac)
      call check_true('a leaky column drains, by about f of what a free one does',                    &
                      drained > F_LEAK * drained_free .and. drained < 2.0_wp * F_LEAK * drained_free,  &
                      drained / drained_free)
      call check('its water closes against the drainage', w0 - w1, drained, 1.0e-10_wp * w0)
   end subroutine test_leaky_drain

   !=======================================================================================!
   !----- WET clay (low Ksat, small suction) under a downpour: infiltration is limited by the     !
   !      clay's conductivity, so the excess ponds and, past the pond's capacity, runs off       !
   !      (Hortonian). A bone-dry clay would instead take it up by suction (Green-Ampt). --------!
   subroutine test_downpour()
      type(soil_profile_t)   :: params
      type(soil_params_t)    :: opts
      type(spac_solver_t)  :: spac
      real(wp) :: theta(10), w0, rain, infil, runoff, drained, bottoms(n_soil_layer_max)
      integer(ik) :: step
      print '(a)', 'test_downpour:'
      call geometric_layer_bottoms(2.0_wp, 10_ik, 3.0_wp, bottoms)
      call build_soil_hydr_params(10_ik, SOIL_RETENTION_VG, bottoms, 0.38_wp, 0.068_wp,            &
           5.6e-7_wp, 0.8_wp, 1.09_wp, -3.37_wp, params)
      theta = 0.36_wp
      call soil_network(spac, params, opts, theta, 0.0_wp, 0.0_wp)
      w0 = column_water(spac) ; rain = 0.0_wp ; infil = 0.0_wp ; runoff = 0.0_wp ; drained = 0.0_wp
      do step = 1_ik, 4_ik
         spac%rain = 1.0e-2_wp                                     ! 36 mm/hr
         call spac_advance(spac, params, opts, humid_air(), no_plants, 900.0_wp)
         rain    = rain + spac%rain * 900.0_wp
         infil   = infil + spac%infil - spac%exfil
         runoff  = runoff + spac%spill
         drained = drained + spac%drain
      end do
      call check_true('infiltration < rain (conductivity-limited)', infil < rain, infil / rain)
      call check_true('the pond filled to its capacity',                                            &
                      spac%water_soil(POND) >= opts%w_pond_max - 1.0e-9_wp, spac%water_soil(POND))
      call check_true('and the rest ran off', runoff > 0.0_wp, runoff)
      call check('water closes: rain = change in store + runoff + drainage', column_water(spac) - w0,        &
                 rain - runoff - drained, 1.0e-10_wp * rain)
   end subroutine test_downpour

   !=======================================================================================!
   !----- A top layer dried to residual takes up 2 mm/h of rain: within 12 h it has re-wetted     !
   !      and almost none of the rain has run off. The pond face's conductivity is the geometric  !
   !      mean of the saturated pond's and the top layer's, and the dry layer's suction drives it;  !
   !      the top layer's conductivity alone would seal a dried surface against rain. -----------!
   subroutine test_dried_surface_rewets()
      type(soil_profile_t)  :: params
      type(soil_params_t)   :: opts
      type(spac_solver_t) :: spac
      type(chydro_forcing_t) :: f
      real(wp) :: theta(10), rain, runoff
      integer(ik) :: step
      print '(a)', 'test_dried_surface_rewets:'
      call loam_params(params)
      theta = 0.25_wp ; theta(1) = 0.08_wp
      call soil_network(spac, params, opts, theta, 0.0_wp, 0.0_wp)
      f = humid_air() ; f%t_ground = 298.0_wp ; f%q_air = 0.02_wp
      rain = 0.0_wp ; runoff = 0.0_wp
      do step = 1_ik, 48_ik
         spac%rain = 2.0_wp / 3600.0_wp                            ! 2 mm/h
         call spac_advance(spac, params, opts, f, no_plants, 900.0_wp)
         rain = rain + spac%rain * 900.0_wp ; runoff = runoff + spac%spill
      end do
      call check_true('the dried top layer re-wets within 12 h of rain',                              &
                      spac%water_soil(2) / (rho_h2o * params%dz(1)) > 0.20_wp, spac%water_soil(2) / (rho_h2o * params%dz(1)))
      call check_true('and under 5 % of the rain runs off', runoff < 0.05_wp * rain, runoff / rain)
   end subroutine test_dried_surface_rewets

   !=======================================================================================!
   !----- Rain (or sub-threshold snow routed to the ground) below the triple point joins the pond:  !
   !      what the pond received in enthalpy is what it gave the soil, kept, and lost as runoff,    !
   !      and the water that infiltrated carried its own enthalpy -- an effective temperature below  !
   !      the melt plateau, not u_liq at the plateau, which would create energy all winter. -------!
   subroutine test_pond_subfreezing_inflow()
      type(soil_profile_t)  :: params
      type(soil_params_t)   :: opts
      type(spac_solver_t) :: spac
      type(chydro_forcing_t) :: f
      real(wp) :: theta(10), dt, e_in, t_infil
      print '(a)', 'test_pond_subfreezing_inflow:'
      call loam_params(params) ; theta = 0.25_wp
      call soil_network(spac, params, opts, theta, 0.0_wp, 0.0_wp)
      spac%rain = 5.0e-6_wp ; spac%t_rain = 270.0_wp                ! 18 mm/day at 3.16 K below the triple point
      f = humid_air() ; f%t_ground = 271.0_wp ; f%q_air = 0.003_wp
      dt = 150.0_wp
      call spac_advance(spac, params, opts, f, no_plants, dt)
      e_in = spac%rain * dt * internal_energy_liquid(270.0_wp)
      call check_true('sub-freezing inflow: all of it infiltrated (pond empty after the step)',          &
                      spac%water_soil(POND) <= 1.0e-9_wp .and. spac%infil > 0.0_wp, spac%water_soil(POND))
      call check('sub-freezing inflow: enthalpy to soil + left in pond + runoff == enthalpy received',    &
                 spac%infil_enth - spac%exfil_enth + spac%e_pond + spac%spill_enth, e_in, 1.0e-9_wp * abs(e_in))
      call check('sub-freezing inflow: the rain brought exactly that enthalpy', spac%rain_enth, e_in, 0.0_wp)
      t_infil = temp_of_liquid(spac%infil_enth / spac%infil)
      call check_true('sub-freezing inflow: infiltrating water carries its own (sub-plateau) enthalpy',    &
                      t_infil < 273.0_wp, t_infil)
   end subroutine test_pond_subfreezing_inflow

   real(wp) function temp_of_liquid(u)
      real(wp), intent(in) :: u
      temp_of_liquid = liquid_energy_to_temp(u)
   end function temp_of_liquid

   !=======================================================================================!
   !----- The aquifer bottom is head-driven against a saturated zone at the column base: a dry     !
   !      column draws water UP (free drainage never can), a wet one still drains down, the water   !
   !      closes with the flux running upward, and the column relaxes to hydrostatic equilibrium    !
   !      with the base (the bottom node's head -> minus half its thickness). -------------------!
   subroutine test_aquifer_head_bc()
      type(soil_profile_t)  :: params
      type(soil_params_t)   :: opts
      type(spac_solver_t) :: spac
      type(chydro_forcing_t) :: f
      real(wp) :: theta(10), drain_dry, drain_wet, w0, drained, theta_bot0, head_n, delta
      integer(ik) :: step, n
      print '(a)', 'test_aquifer_head_bc:'
      n = 10_ik
      call loam_params(params)
      opts%water_bottom_bc = SOIL_BC_AQUIFER
      f = humid_air() ; f%r_aero = 1.0e6_wp                       ! suppress soil evaporation
      theta = 0.12_wp
      call soil_network(spac, params, opts, theta, 0.0_wp, 0.0_wp)
      call spac_advance(spac, params, opts, f, no_plants, 600.0_wp)
      drain_dry = spac%drain
      call check_true('AQUIFER: dry column draws water UP from the saturated base (drainage < 0)',        &
                      drain_dry < 0.0_wp, drain_dry)
      theta = 0.40_wp
      call soil_network(spac, params, opts, theta, 0.0_wp, 0.0_wp)
      call spac_advance(spac, params, opts, f, no_plants, 600.0_wp)
      drain_wet = spac%drain
      call check_true('AQUIFER: wet column still drains DOWNWARD (drainage > 0)', drain_wet > 0.0_wp, drain_wet)
      theta = 0.15_wp
      call soil_network(spac, params, opts, theta, 0.0_wp, 0.0_wp)
      w0 = column_water(spac) ; drained = 0.0_wp ; theta_bot0 = 0.15_wp
      do step = 1_ik, 40_ik
         call spac_advance(spac, params, opts, f, no_plants, 600.0_wp)
         drained = drained + spac%drain
      end do
      call check('AQUIFER: water closes while drawing upward', column_water(spac) - w0, -drained, 1.0e-10_wp * w0)
      call check_true('AQUIFER: the dry column WETS from below',                                           &
                      spac%water_soil(1+n) / (rho_h2o * params%dz(n)) > theta_bot0, spac%water_soil(1+n) / (rho_h2o * params%dz(n)))
      do step = 1_ik, 1200_ik
         call spac_advance(spac, params, opts, f, no_plants, 3600.0_wp)
      end do
      head_n = spac%psi_soil(1+n) / grav_head
      delta  = 0.5_wp * params%dz(n)
      call check_true('AQUIFER: bottom layer reaches hydrostatic equilibrium (head_n -> -Delta)',           &
                      abs(head_n + delta) < 0.2_wp, head_n + delta)
      call check_true('AQUIFER: and the residual bottom flux has died', abs(spac%drain / 3600.0_wp) < 1.0e-5_wp, &
                      spac%drain / 3600.0_wp)
   end subroutine test_aquifer_head_bc

   !=======================================================================================!
   !----- Ground evaporation is an AREA-weighted tile flux, E = snow_free_frac rho dq / (r_aero +    !
   !      r_soil): linear in the snow-free share and exactly 0 under full cover, with a dry surface   !
   !      layer present (theta_1 = 0.30 is below dsl_theta_init * theta_sat = 0.344), where dividing   !
   !      r_aero by (1 - snowfac) instead would not be linear. --------------------------------------!
   subroutine test_snow_free_evap()
      type(soil_profile_t)   :: params
      type(soil_params_t)    :: opts
      type(chydro_forcing_t) :: f
      real(wp) :: e_full, e_half, e_none, de1, de2, head
      print '(a)', 'test_snow_free_evap:'
      call loam_params(params)
      f%t_ground = 305.0_wp ; f%q_air = 0.001_wp ; f%rho_air = 1.2_wp ; f%r_aero = 50.0_wp
      head = soil_theta_to_head(SOIL_RETENTION_VG, 0.30_wp, 0.43_wp, 0.078_wp, 3.6_wp, 1.56_wp)
      f%snow_free_frac = 1.0_wp
      call ground_evaporation_rate(0.30_wp, head, params, f, opts, e_full, de1, de2)
      f%snow_free_frac = 0.5_wp
      call ground_evaporation_rate(0.30_wp, head, params, f, opts, e_half, de1, de2)
      f%snow_free_frac = 0.0_wp
      call ground_evaporation_rate(0.30_wp, head, params, f, opts, e_none, de1, de2)
      call check_true('bare ground evaporates (non-trivial case)', e_full > 1.0e-8_wp, e_full)
      call check('half snow cover halves E exactly (area weighting)', e_half, 0.5_wp * e_full,             &
                 1.0e-14_wp * max(e_full, 1.0e-12_wp))
      call check('full snow cover gives exactly zero soil evaporation', e_none, 0.0_wp, 1.0e-30_wp)
   end subroutine test_snow_free_evap

   !=======================================================================================!
   !----- A wetter soil evaporates MORE. The dry-surface-layer resistance d_sl/(D_v tau) is built on  !
   !      the air-dry water content's tortuosity, so it falls as the layer wets and the response is   !
   !      monotone, approaching the no-DSL limit continuously at theta_init. Built from the bulk       !
   !      layer's air-filled porosity instead, the tortuosity would fall faster than the layer thins  !
   !      and a wetter soil would evaporate less. ----------------------------------------------!
   subroutine test_evap_moisture_response()
      type(soil_profile_t)   :: params
      type(soil_params_t)    :: opts
      type(chydro_forcing_t) :: f
      real(wp) :: theta, e_prev, e_now, e_sat, e_below, theta_init
      integer  :: i
      logical  :: monotone
      print '(a)', 'test_evap_moisture_response:'
      call loam_params(params)
      f%t_ground = 300.0_wp ; f%q_air = 0.005_wp ; f%rho_air = 1.2_wp ; f%r_aero = 100.0_wp
      theta_init = opts%dsl_theta_init * params%theta_sat(1)
      monotone = .true. ; e_prev = -1.0_wp
      do i = 1, 40
         theta = params%theta_res(1) + real(i, wp) / 40.0_wp * (theta_init - params%theta_res(1))
         e_now = loam_evap(theta, params, f, opts)
         if (e_now < e_prev - 1.0e-14_wp) monotone = .false.
         e_prev = e_now
      end do
      call check_true('ground evaporation increases monotonically with soil moisture', monotone, 0.0_wp)
      e_sat   = loam_evap(theta_init, params, f, opts)
      e_below = loam_evap(theta_init - 1.0e-3_wp * (theta_init - params%theta_res(1)), params, f, opts)
      call check_true('approaches the DSL-free limit continuously at theta_init',                       &
                      e_below > 0.8_wp * e_sat, e_below / max(e_sat, 1.0e-30_wp))
   end subroutine test_evap_moisture_response

   !----- Ground evaporation of the loam's top layer at water content th. ----------------------!
   real(wp) function loam_evap(th, params, f, opts) result(e)
      real(wp),               intent(in) :: th
      type(soil_profile_t),   intent(in) :: params
      type(chydro_forcing_t), intent(in) :: f
      type(soil_params_t),    intent(in) :: opts
      real(wp) :: de1, de2
      call ground_evaporation_rate(th, soil_theta_to_head(SOIL_RETENTION_VG, th, 0.43_wp, 0.078_wp,     &
                                   3.6_wp, 1.56_wp), params, f, opts, e, de1, de2)
   end function loam_evap

end program test_column_hydrology
