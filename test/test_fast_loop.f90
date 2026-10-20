! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_fast_loop -- the SITE-LEVEL fast-biophysics driver over the state-hub-owned per-patch     !
! reservoirs. Builds a two-patch site (one vegetated, one bare ground), seeds the reservoirs,     !
! and runs the fast loop:                                                                        !
!   1. CONSERVATION: every patch's whole-column energy + water budgets close (n_fail == 0) --     !
!      the fast loop conserves over the reservoirs the state hub owns, for a VEGETATED and a      !
!      BARE (zero-cohort) patch alike.                                                            !
!   2. ACTIVITY: the reservoirs evolve (CAS warms under the absorbed shortwave; soil responds).   !
!   3. STEPPER WIRING + GATE: advance_one_step runs the fast loop only when fast_biophysics_on    !
!      AND a fast context is supplied; with the gate off the reservoirs are untouched.            !
!==========================================================================================!
program test_fast_loop
   use meds_kinds,               only : wp, ik
   use meds_config,              only : meds_config_t
   use meds_cohort_state_vars,   only : set_cohort_size_from_carbon, hold_cohort_tissue_temps
   use meds_polygon_state_vars,  only : polygon_t
   use meds_init,                only : init_bare_ground, add_cohort, finalize_init
   use meds_test_soil_grid, only : geometric_layer_bottoms
   use meds_soil_column_config, only : build_soil_hydr_params, build_soil_therm_params, top_soil_psi
   use meds_soil_hydr_lib,   only : soil_retention_vg
   use meds_plant_hydr_lib,  only : tissue_psi_to_water
   use meds_slow_dynamics,   only : settle_tissue_water
   use meds_fast_dynamics,       only : fast_context_t, init_fast_reservoirs, fast_dynamics, &
                                        build_fast_context
   use meds_fast_types,          only : apply_hydraulics_config, apply_canopy_film_config
   use meds_fast_config, only : build_leaf_photo_table, build_integrator_opts
   use meds_stepper,             only : advance_one_step
   use meds_test_support, only : banner, build_test_config, check, check_close
   use meds_time,                only : meds_time_t, time_advance_seconds
   use meds_forcing_config,      only : MET_BACKEND_ED_DEFAULT, SWPART_CLEARIDX, METAVG_END
   use meds_forcing_types,       only : met_source_t, met_cursor_t, met_forcing_t
   use meds_met_driver,          only : met_open, met_cursor_init, met_close, met_advance, met_instant
   use meds_output_types,        only : output_buffers_t
   use meds_diag_vars,           only : PY_TAIR, PY_QAIR, PY_PAR_BEAM, PY_COSZ, PY_SW_IN, PY_LWDOWN, PY_PAR_DIFFUSE, PY_RHO_AIR, &
                                        PY_WIND, PA_CAS_TEMP, PA_SKIN_TEMP, PA_LW_UP
   use meds_constants,           only : stefan_boltzmann, n_soil_layer_max
   use meds_netcdf_c
   use iso_c_binding,            only : c_int, c_size_t, c_double
   implicit none
   real(wp) :: test_bottoms(n_soil_layer_max)   ! the geometric layers these checks were written on

   integer(ik), parameter :: nsl = 10_ik
   type(meds_config_t)  :: cfg
   type(polygon_t)      :: poly
   type(fast_context_t) :: ctx
   real(wp)    :: we, ww, t_cas0, t_cas1, theta0_1, theta1_1, mass0_leaf, mass1_leaf, mass2_leaf
   real(wp)    :: t_leaf0, t_wood0
   real(wp)    :: cbal0, cbal1
   integer(ik) :: n_budget_fail

   call banner('site-level fast-biophysics loop (state-hub reservoirs)')

   cfg = build_test_config()
   cfg%fast_biophysics_on = .true.
   cfg%dt_fast            = 150.0_wp
   cfg%n_fast_per_slow    = 48_ik           ! 48 x 150 s = 2 h of fast integration per slow step

   !----- Build the column config + reference met inside the fast context. -----------------!
   call geometric_layer_bottoms(2.0_wp, nsl, 3.0_wp, test_bottoms)
   call build_soil_hydr_params(nsl, soil_retention_vg, test_bottoms, 0.43_wp, 0.078_wp,            &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, -3.37_wp, ctx%col_config%soil)
   call build_soil_therm_params(nsl, 3.0_wp, 0.15_wp, 2.0e6_wp, ctx%col_config%soil_thermal)
   call apply_hydraulics_config(cfg%pft, ctx%col_config%hydraulics_table)
   call apply_canopy_film_config(cfg%pft, ctx%col_config%leaf_surf_water_max, ctx%col_config%wood_surf_water_max)
   call build_leaf_photo_table(cfg, ctx%col_config%leaf_photo)
   ctx%col_config%integrator = build_integrator_opts(cfg)
   ctx%air_temp = 290.0_wp ; ctx%rad_sw_top = 500.0_wp ; ctx%rad_sw_ground = 75.0_wp
   ctx%soil_temp_init = 288.0_wp

   !----- Two patches: patch 1 vegetated, patch 2 BARE (zero cohorts). ---------------------!
   call init_bare_ground(poly, cfg, 2_ik)
   call use_test_column(poly)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)     ! patch 1: one tree cohort
   call finalize_init(poly)
   call init_fast_reservoirs(poly, ctx)

   t_cas0    = poly%patch%cas(1)%can_temp
   theta0_1  = poly%patch%soil_water(1)%theta(1)
   !----- init_fast_reservoirs starts the cohort in its patch's state: holding the water its tissues !
   !      hold at the top soil layer's potential, so the fast loop's first touch reads a real store,  !
   !      not a sentinel to seed. --------------------------------------------------------------------!
   mass0_leaf = poly%cohort%leaf_water_mass(1)
   t_leaf0    = poly%cohort%leaf_temp(1) ; t_wood0 = poly%cohort%wood_temp(1)

   !=== 1+2. Run the fast loop directly; conservation + activity + per-cohort persistence. =!
   call fast_dynamics(poly, ctx, cfg, worst_energy=we, worst_water=ww, n_budget_fail=n_budget_fail)
   t_cas1    = poly%patch%cas(1)%can_temp
   theta1_1  = poly%patch%soil_water(1)%theta(1)
   mass1_leaf = poly%cohort%leaf_water_mass(1)       ! lazy-init seeded THIS call, then evolved once

   call check(n_budget_fail == 0_ik, 'whole-column budgets closed on every patch (n_fail == 0)')
   call check(we < 1.0e-3_wp, 'whole-column energy residual tiny')
   call check(ww < 1.0e-8_wp, 'whole-column water residual tiny')
   call check(abs(t_cas1 - t_cas0) > 0.05_wp, 'CAS temperature evolved under the fast loop')
   call check(poly%patch%cas(2)%can_temp > 200.0_wp, 'bare (zero-cohort) patch fast step stayed physical')
   call check_close(mass0_leaf, tissue_psi_to_water(top_soil_psi(ctx%col_config%soil, theta0_1),           &
                    ctx%col_config%hydraulics_table%pft(1)%leaf_curve, poly%cohort%leaf_carbon(1)),     &
                    1.0e-12_wp, 'a cohort starts with its leaf water at its patch''s top-soil potential')
   !----- To round-off: the temperatures are what the energies the cohort was born with give. ------!
   call check(abs(t_leaf0 - t_cas0) < 1.0e-9_wp .and. abs(t_wood0 - t_cas0) < 1.0e-9_wp,            &
              'a cohort starts at its patch''s canopy-air temperature')
   call check(mass1_leaf > 0.0_wp .and. mass1_leaf /= mass0_leaf,                                 &
              'the fast loop evolved the birth water and wrote it back')
   !----- A SECOND call starts from mass1_leaf, so its own evolution is a clean round-trip proof:   !
   !      the fast loop READ leaf_water_mass from the cohort block, evolved it, and WROTE it back   !
   !      (persist) onto the SoA, not just a local/scratch copy. -----------------------------------!
   call fast_dynamics(poly, ctx, cfg, worst_energy=we, worst_water=ww, n_budget_fail=n_budget_fail)
   mass2_leaf = poly%cohort%leaf_water_mass(1)
   call check(abs(mass2_leaf - mass1_leaf) > 1.0e-9_wp,                                           &
              'per-cohort leaf water mass evolved + persisted on the cohort block (2nd fast_dynamics call)')
   !----- Fast->slow carbon bridge: the vegetated cohort accumulated positive GROSS GPP. -------!
   call check(poly%cohort%gpp_accum(1) > 0.0_wp, 'fast loop accumulated positive gross GPP (fast->slow bridge)')

   !=== 3. Stepper wiring + gate. The gate is cfg%fast_biophysics_on; when it is ON a fast      !
   !    context is REQUIRED (advance_one_step error-stops otherwise -- not exercised here,      !
   !    since asserting an abort would kill the process). ======================================!
   block
      real(wp) :: t_before
      !----- Gate OFF (fast_biophysics_on = .false.): the fast loop must NOT run, even with a ctx. !
      cfg%fast_biophysics_on = .false.
      call init_fast_reservoirs(poly, ctx)
      t_before = poly%patch%cas(1)%can_temp
      call advance_one_step(poly, cfg, ctx)     ! gate off -> fast loop skipped
      call check_close(poly%patch%cas(1)%can_temp, t_before, 1.0e-12_wp, &
                       'fast loop must NOT run when fast_biophysics_on is off')
      !----- Gate ON (fast_biophysics_on = .true. + context): the hook fires, reservoirs change. -!
      cfg%fast_biophysics_on = .true.
      call init_fast_reservoirs(poly, ctx)
      t_before = poly%patch%cas(1)%can_temp
      call advance_one_step(poly, cfg, ctx)
      call check(abs(poly%patch%cas(1)%can_temp - t_before) > 0.05_wp, &
                 'advance_one_step ran the fast loop when gated on')

      !----- CAS DEPTH tracks the stand: refresh_canopy_depth sets it to the tallest cohort plus a   !
      !      freeboard (docs/science/canopy_air_space_biophysics.md, The three CAS twins).           !
      !      cas_mass_capacity IS the canopy air's heat capacity, so a fixed depth would give a 1 m  !
      !      gap and a 35 m canopy the same response time. Asserted end-to-end through the stepper,  !
      !      because the slow loop owns it. ---------------------------------------------------------!
      block
         real(wp)    :: h_top, expect
         integer(ik) :: i
         h_top = 0.0_wp
         do i = 1_ik, poly%patch%cohort_count(1)
            h_top = max(h_top, poly%cohort%height(poly%patch%cohort_offset(1) + i - 1_ik))
         end do
         expect = max(cfg%aero%min_canopy_depth, h_top + cfg%aero%canopy_freeboard)
         call check_close(poly%patch%cas(1)%can_depth, expect, 1.0e-12_wp,                        &
                          'CAS depth = tallest cohort + freeboard')
         call check(abs(poly%patch%cas(1)%can_depth - 20.0_wp) > 1.0e-6_wp,                       &
                    'CAS depth is not pinned at a fixed 20 m')
      end block
   end block

   !=== 4. END-TO-END fast->slow handoff via the stepper. The SAME carbon-mode site is advanced  !
   !    with the fast loop OFF then ON, with gpp_ref = 0 so the stub adds nothing. Turnover is    !
   !    identical between the two, so more carbon under fast-ON is EXACTLY the fast GPP flowing    !
   !    through gpp_accum into carbon growth -- the fast->slow bridge, isolated from turnover.  ===!
   cfg%gpp_ref       = 0.0_wp

   cfg%fast_biophysics_on = .false.                    ! stub GPP = 0; carbon change is pure turnover
   call init_bare_ground(poly, cfg, 1_ik)
   call use_test_column(poly)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
   call finalize_init(poly)
   call advance_one_step(poly, cfg, ctx)
   cbal0 = poly%cohort%leaf_carbon(1) + poly%cohort%fineroot_carbon(1)                          &
         + poly%cohort%wood_carbon(1) + poly%cohort%nonstructural_carbon(1)

   cfg%fast_biophysics_on = .true.                     ! fast loop -> gpp_accum -> carbon growth
   call init_bare_ground(poly, cfg, 1_ik)
   call use_test_column(poly)
   call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
   call finalize_init(poly)
   call init_fast_reservoirs(poly, ctx)
   call advance_one_step(poly, cfg, ctx)
   cbal1 = poly%cohort%leaf_carbon(1) + poly%cohort%fineroot_carbon(1)                          &
         + poly%cohort%wood_carbon(1) + poly%cohort%nonstructural_carbon(1)

   call check(cbal1 > cbal0, 'fast GPP raised carbon vs fast-off (gpp_ref=0): the fast->slow bridge is live')
   !----- The fast loop accumulates POSITIVE maintenance respiration (leaf Rd + stem + root),      !
   !      which compute_carbon_allocation nets out of GPP                                          !
   !      (docs/science/plant_carbon_allocation.md eq. (1)). --------------------------------------!
   call check(poly%cohort%leaf_resp_accum(1) > 0.0_wp .and. poly%cohort%stem_resp_accum(1) > 0.0_wp &
              .and. poly%cohort%root_resp_accum(1) > 0.0_wp,                                        &
              'fast loop accumulated positive leaf/stem/root maintenance respiration')

   !=== 5. FORCING-DRIVEN diurnal cycle, end to end. Drive the SAME fast loop from a real             !
   !    MEDS forcing NetCDF (SW peaks ~17 UTC at Ithaca). A 2 h NIGHT window (02-04 UTC, SW=0) must   !
   !    accumulate ~0 GPP; a 2 h DAY window (15-17 UTC, high SW) must accumulate clearly more -- i.e. !
   !    the forcing's diurnal shortwave really propagates through met_instant -> the fast loop.  =====!
   block
      type(met_source_t) :: drv
      type(met_cursor_t) :: cur
      real(wp)           :: gpp_night, gpp_day, tcas_n, light_night, light_day
      real(wp)           :: g1_ref, g2_ref, g3_ref   ! 1-thread multi-patch GPP, for the thread-invariance check
      cfg%gpp_ref = 0.0_wp     ! isolate the FORCING-driven GPP
      call build_fast_context(cfg, ctx)                          ! rebuild ctx WITH the RT optics table (rad_opt)
      call write_diurnal_forcing('test_fast_loop_forcing.nc')
      cfg%forcing%forcing_on   = .true.
      cfg%forcing%backend      = MET_BACKEND_ED_DEFAULT
      cfg%forcing%path         = 'test_fast_loop_forcing.nc'
      cfg%forcing%grid_index   = 1_ik ; cfg%forcing%dt_forcing = 3600.0_wp
      cfg%forcing%avg_convention = METAVG_END ; cfg%forcing%sw_partition = SWPART_CLEARIDX
      cfg%forcing%latitude_deg = 42.44_wp ; cfg%forcing%longitude_deg = -76.50_wp
      cfg%forcing%recycle = .false.
      cfg%forcing%co2_const = 415.0_wp          ! CO2 comes from config, not the file (docs/science/forcing.md §12)

      call init_bare_ground(poly, cfg, 1_ik)

      call use_test_column(poly)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
      call finalize_init(poly)
      call met_open(drv, cfg%forcing)
      call met_cursor_init(drv, cur, 1_ik, cfg%forcing%latitude_deg,      &
                           cfg%forcing%longitude_deg,                            &
                           cfg%forcing%elevation_m)

      call init_fast_reservoirs(poly, ctx)
      call fast_dynamics(poly, ctx, cfg, met_src=drv, met_cur=cur,               &
                               step_start=meds_time_t(2020_ik,7_ik,1_ik,2_ik))    ! 02-04 UTC (night)
      gpp_night = poly%cohort%gpp_accum(1)
      light_night = poly%cohort%light_hours_accum(1)
      tcas_n    = poly%patch%cas(1)%can_temp      ! night canopy-air temp after the SW=0 window

      call init_fast_reservoirs(poly, ctx)
      call fast_dynamics(poly, ctx, cfg, met_src=drv, met_cur=cur,               &
                               step_start=meds_time_t(2020_ik,7_ik,1_ik,15_ik))   ! 15-17 UTC (day)
      gpp_day = poly%cohort%gpp_accum(1)
      light_day = poly%cohort%light_hours_accum(1)

      !=== The FORCING ECHO (MEDS_FORCING_DESIGN.md §6.7): with the polygon diagnostics and the     !
      !    fast tier on, a step must record exactly the forcing the reader hands the fast loop --    !
      !    the dt-weighted mean of the sub-step samples in the polygon block (poly%diag, once per    !
      !    polygon), and each sample itself in the fast staging. A second cursor on the same source  !
      !    replays the reader's samples. ---------------------------------------------------------!
      block
         type(output_buffers_t) :: ob
         type(met_cursor_t)     :: cur2
         type(met_forcing_t)    :: m
         type(meds_time_t)      :: t0, ts
         real(wp)    :: sum_t, sum_q, sum_pb, sum_cz, sum_sw, w
         integer(ik) :: isub
         logical     :: exact
         t0 = meds_time_t(2020_ik,7_ik,1_ik,15_ik)
         poly%diag%active = .true.
         ob%fast_on = .true.
         call init_fast_reservoirs(poly, ctx)
         call fast_dynamics(poly, ctx, cfg, met_src=drv, met_cur=cur, step_start=t0, out_bufs=ob)
         call met_cursor_init(drv, cur2, 1_ik, cfg%forcing%latitude_deg, cfg%forcing%longitude_deg, &
                              cfg%forcing%elevation_m)
         sum_t = 0.0_wp ; sum_q = 0.0_wp ; sum_pb = 0.0_wp ; sum_cz = 0.0_wp ; sum_sw = 0.0_wp
         exact = .true.
         do isub = 1_ik, cfg%n_fast_per_slow
            ts = time_advance_seconds(t0, (real(isub, wp) - 1.0_wp + cfg%forcing_sample_frac) * cfg%dt_fast)
            call met_advance(drv, cur2, ts) ; m = met_instant(drv, cur2, ts)
            sum_t = sum_t + m%tair_k * cfg%dt_fast ; sum_q = sum_q + m%qair * cfg%dt_fast
            sum_pb = sum_pb + m%par_beam * cfg%dt_fast ; sum_cz = sum_cz + m%cosz * cfg%dt_fast
            sum_sw = sum_sw + m%swdown() * cfg%dt_fast
            associate (f => ob%fast_forcing(:, isub))
               exact = exact .and. f(PY_QAIR) == m%qair .and. f(PY_LWDOWN) == m%lwdown               &
                       .and. f(PY_PAR_DIFFUSE) == m%par_diffuse .and. f(PY_COSZ) == m%cosz           &
                       .and. f(PY_RHO_AIR) == m%rho_air .and. f(PY_WIND) == m%wind                   &
                       .and. f(PY_TAIR) == m%tair_k .and. f(PY_SW_IN) == m%swdown()
            end associate
         end do
         w = poly%diag%w
         call check(ob%n_fast_sub == cfg%n_fast_per_slow .and. exact,                              &
                    'the fast tier stages each sub-step''s forcing sample exactly')
         !----- The soil slabs fold the ACTIVE layers only, so the inactive tail reads back as the   !
         !      fill value, as on the coarse tiers. The fixture needs such a tail. ------------------!
         call check(ctx%col_config%soil%n_active < n_soil_layer_max,                               &
                    'fixture has an inactive soil tail to be worth anything')
         call check(ob%fast_n_site == 1_ik .and. count(ob%fast_soil_live(1:n_soil_layer_max))       &
                    == ctx%col_config%soil%n_active, 'the fast tier folds only the active soil layers')
         !----- The patch axis: each patch's staged row is its own, unweighted -- the last           !
         !      sub-step's is the state the step wrote back -- and the site mean is their area-      !
         !      weighted sum, in patch order (docs/science/diagnostics.md §4, The FAST tier by       !
         !      patch). -----------------------------------------------------------------------------!
         block
            integer(ik) :: ip, nsub, na, kb
            real(wp)    :: site_cas
            logical     :: own, folded, skin
            nsub = cfg%n_fast_per_slow ; na = ctx%col_config%soil%n_active
            own = ob%fast_n_patch == poly%patch%n ; site_cas = 0.0_wp ; skin = .true.
            do ip = 1_ik, poly%patch%n
               kb = (ip - 1_ik) * n_soil_layer_max
               own = own .and. ob%fast_patch(PA_CAS_TEMP, ip, nsub) == poly%patch%cas(ip)%can_temp     &
                         .and. all(ob%fast_soil_temp_patch(kb+1_ik:kb+na, nsub)                        &
                                   == poly%patch%soil_energy(ip)%soil_temp(1:na))
               site_cas = site_cas + poly%patch%area(ip) * ob%fast_patch(PA_CAS_TEMP, ip, nsub)
               !----- The skin temperature is the black-body temperature of the longwave up        !
               !      (docs/science/diagnostics.md §4, The skin temperature). ---------------------!
               skin = skin .and. abs(stefan_boltzmann * ob%fast_patch(PA_SKIN_TEMP, ip, nsub)**4        &
                                     - ob%fast_patch(PA_LW_UP, ip, nsub)) <= 1.0e-9_wp * ob%fast_patch(PA_LW_UP, ip, nsub)
            end do
            call check(skin, 'the skin temperature radiates the patch''s longwave up as a black body')
            folded = ob%fast_site(PA_CAS_TEMP, nsub) == site_cas
            call check(own, 'the fast tier stages each patch''s own sample and soil column')
            call check(folded, 'the fast tier''s site mean is the area-weighted sum of the patch rows')
         end block
         call check_close(w, cfg%n_fast_per_slow * cfg%dt_fast, 1.0e-9_wp,                         &
                          'the polygon block weighs each sub-step once (not once per patch)')
         call check_close(poly%diag%v(PY_TAIR) / w, sum_t / w, 1.0e-12_wp,                         &
                          'the polygon block holds the step-mean air temperature the loop used')
         call check_close(poly%diag%v(PY_QAIR) / w, sum_q / w, 1.0e-15_wp,                         &
                          'the polygon block holds the step-mean humidity')
         call check_close(poly%diag%v(PY_PAR_BEAM) / w, sum_pb / w, 1.0e-10_wp,                    &
                          'the polygon block holds the step-mean direct-beam PAR')
         call check_close(poly%diag%v(PY_SW_IN) / w, sum_sw / w, 1.0e-10_wp,                       &
                          'the polygon block holds the step-mean total shortwave')
         call check(poly%diag%v(PY_COSZ) / w > 0.5_wp,                                             &
                    'a 15-17 UTC summer window at Ithaca has the sun high (mean cos zenith > 0.5)')
         poly%diag%active = .false.
      end block
      call met_close(drv)

      call check(gpp_night < 1.0e-9_wp, 'night forcing window -> ~zero GPP (SW=0 propagated through met_instant)')
      call check(gpp_day > 1.0e-7_wp,   'day forcing window -> positive GPP')
      call check(gpp_day > 10.0_wp * max(gpp_night, 1.0e-30_wp), 'GPP tracks the diurnal shortwave (day >> night)')
      !----- The phenology light count: no step of the night window has PAR above par_min, every   !
      !      step of the mid-afternoon window does, so the cohort counts 0 h and the window's hours. !
      call check(light_night == 0.0_wp, 'night window -> no hours of light counted')
      call check(abs(light_day - cfg%n_fast_per_slow * cfg%dt_fast / 3600.0_wp) < 1.0e-9_wp, &
                 'day window -> every fast step counted as light (the window''s hours)')

      !=== MULTI-PATCH x FORCING. ===========================================================================!
      !    The forced block above runs ONE patch, so it cannot see the patch loop's interaction with the   !
      !    met reader at all: met_advance/met_instant are called once per sub-step, outside the patch        !
      !    loop, and a single-patch fixture is order-identical whether or not they are, so it proves         !
      !    nothing about it. This runs the SAME window over THREE patches carrying identical cohorts and     !
      !    asserts they get identical GPP: the met stream is site-uniform, so any per-patch divergence       !
      !    means the reader's state leaked across the patch loop. -------------------------------------------!
      block
         type(met_source_t) :: drv3
         type(met_cursor_t) :: cur3
         real(wp) :: g1, g2, g3
         call init_bare_ground(poly, cfg, 3_ik)
         call use_test_column(poly)
         call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
         call add_cohort(poly, cfg, 2_ik, 1_ik, 0.3_wp, 16.0_wp)
         call add_cohort(poly, cfg, 3_ik, 1_ik, 0.3_wp, 16.0_wp)
         call finalize_init(poly)
         call met_open(drv3, cfg%forcing)
         call met_cursor_init(drv3, cur3, 1_ik, cfg%forcing%latitude_deg,      &
                              cfg%forcing%longitude_deg,                            &
                              cfg%forcing%elevation_m)
         call init_fast_reservoirs(poly, ctx)
         call fast_dynamics(poly, ctx, cfg, met_src=drv3, met_cur=cur3,             &
                            step_start=meds_time_t(2020_ik,7_ik,1_ik,15_ik))
         call met_close(drv3)
         g1 = poly%cohort%gpp_accum(poly%patch%cohort_offset(1))
         g2 = poly%cohort%gpp_accum(poly%patch%cohort_offset(2))
         g3 = poly%cohort%gpp_accum(poly%patch%cohort_offset(3))
         call check(g1 > 1.0e-7_wp, 'multi-patch forced run produces GPP at all')
         !----- The reader-hoist assertion: patches 2 and 3 are BOTH downstream of patch 1 in the loop, !
         !      so if the met reader were advanced per patch they would see different                   !
         !      streams. They must agree EXACTLY. ------------------------------------------------------!
         call check(abs(g3 - g2) <= 1.0e-14_wp * max(abs(g2), 1.0e-30_wp),                   &
                    'patches 2 and 3 see the SAME met stream (reader not advanced per patch)')
         !----- adapt_dt_last (the adaptive controller's warm start) has per-patch storage and rides the !
         !      patch lockstep, so identical patches must agree EXACTLY. Held on scratch shared across  !
         !      the patch loop it would be loop-carried -- patch 1 cold-starts and patches 2..N inherit !
         !      their predecessor's step, ~1e-8 relative apart on IDENTICAL cohorts and forcing -- and  !
         !      the answer would depend on patch ORDER serially and on THREAD SCHEDULING in parallel    !
         !      (docs/science/numerical_scheme.md §6a). ------------------------------------------------!
         call check(abs(g1 - g2) <= 1.0e-14_wp * max(abs(g2), 1.0e-30_wp),                   &
                    'patch 1 == patch 2 exactly: no loop-carried controller warm start')
         g1_ref = g1 ; g2_ref = g2 ; g3_ref = g3     ! the 1-thread reference for the thread-invariance check below
      end block

      !=== THREAD INVARIANCE (docs/science/numerical_scheme.md §6a). =====================================!
      !    Re-run the identical multi-patch forced window at [run].n_threads = 4 and require the answer  !
      !    to be BYTE-identical to the 1-thread run above. In a build without OpenMP flags the directives !
      !    are comments and this passes trivially -- but it is the assertion that fails the moment the    !
      !    reductions stop being folded in patch order, or a per-thread buffer starts being shared, and    !
      !    those are exactly the two failure modes that produce plausible wrong numbers rather than a       !
      !    crash. (The two traps behind them: ifx builds `private`/`firstprivate` copies of a derived     !
      !    type through a STATIC mold that every thread writes, and the kernel libraries' local arrays    !
      !    sit in static storage without -auto; see docs/building.md.) -----------------------------------!
      block
         type(met_source_t)  :: drv4
         type(met_cursor_t)  :: cur4
         type(fast_context_t) :: ctx_mt
         real(wp) :: h1, h2, h3
         call init_bare_ground(poly, cfg, 3_ik)
         call use_test_column(poly)
         call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
         call add_cohort(poly, cfg, 2_ik, 1_ik, 0.3_wp, 16.0_wp)
         call add_cohort(poly, cfg, 3_ik, 1_ik, 0.3_wp, 16.0_wp)
         call finalize_init(poly)
         call met_open(drv4, cfg%forcing)
         call met_cursor_init(drv4, cur4, 1_ik, cfg%forcing%latitude_deg,      &
                              cfg%forcing%longitude_deg,                            &
                              cfg%forcing%elevation_m)
         call init_fast_reservoirs(poly, ctx)
         ctx_mt = ctx ; ctx_mt%patch_threads = 4_ik
         call fast_dynamics(poly, ctx_mt, cfg, met_src=drv4, met_cur=cur4,             &
                            step_start=meds_time_t(2020_ik,7_ik,1_ik,15_ik))
         call met_close(drv4)
         h1 = poly%cohort%gpp_accum(poly%patch%cohort_offset(1))
         h2 = poly%cohort%gpp_accum(poly%patch%cohort_offset(2))
         h3 = poly%cohort%gpp_accum(poly%patch%cohort_offset(3))
         call check(h1 == g1_ref .and. h2 == g2_ref .and. h3 == g3_ref,                      &
                    'n_threads = 4 reproduces n_threads = 1 EXACTLY')
      end block
      !----- Net-longwave wiring: at night (SW=0) the sky (LWdown=340) is far colder than the surface !
      !      (~sigma*T^4 ~ 385 W/m2), so the two-stream NET longwave (leaf + ground) is a radiative    !
      !      LOSS that drives the canopy air below its ~288 K start and the 288.5-291.5 K atmosphere. !
      !----- The threshold, measured 2x2 (2 h night window; LW off = abs_lw/abs_lw_ground at 0;     !
      !      store off = the leaf and wood heat capacities at 0):                                     !
      !                        LW on     LW off    LW signal                                        !
      !          store off     285.95    288.39    2.44 K                                           !
      !          store on      286.80    288.41    1.60 K                                           !
      !      The tissue heat store warms the night canopy air by 0.85 K, which is REAL: leaf and wood !
      !      carry a heat capacity comparable to the canopy air's, so tissue shedding its excess      !
      !      lifts the air. It also eats part of the LW signal, so the threshold has to sit between   !
      !      the store-on pair, not the store-off pair. 288.0 K discriminates in BOTH configurations  !
      !      (margins 1.20/0.41 K with the store on, 2.05/0.39 K with it off). Re-measure both        !
      !      controls if either is touched. ----------------------------------------------------------!
      call check(tcas_n < 288.0_wp, 'net longwave wired: canopy air radiatively cools (net-LW loss to the cold sky)')
      write(*,'(a,es10.3,a,es10.3,a)') '   (forcing GPP: night=', gpp_night, '  day=', gpp_day, ' kgC/plant)'
      write(*,'(a,f7.3,a)') '   (night CAS temp=', tcas_n, ' K, radiatively cooled below the ~288 K start)'
   end block

   !=== 6. FORCING + canopy-RT cohort ORDER. Two same-PFT cohorts -- a TALL overstory and              !
   !    a SHORT understory -- share one patch, with EQUAL per-cohort LAI (so per-leaf PAR is a clean  !
   !    overstory-vs-understory discriminator, not diluted by differing leaf area). The two-stream    !
   !    canopy RT + the wind cascade must hand the TOP cohort more absorbed PAR AND more in-canopy     !
   !    wind (=> higher boundary-layer gb) than the shaded understory below it, so with the leaf's     !
   !    hydraulic water-stress NEUTRALIZED (below) its GPP per leaf area is higher. A cohort order     !
   !    fed to the RT or the aerodynamics the wrong way round (top<->bottom inverted) flips this, so   !
   !    it guards BOTH the apply_rt_forcing `perm` reversal AND the aero_bottom_to_top reversal.       !
   !    (Left hydraulically limited, the taller cohort's more negative psi_leaf would legitimately     !
   !    close its stomata and INVERT the GPP ranking -- correct physics, but it masks the wiring.) ===!
   block
      type(met_source_t) :: drv
      type(met_cursor_t) :: cur
      integer(ik)        :: itall, ishort, ii
      real(wp)           :: gpp_la_top, gpp_la_under
      real(wp), parameter :: LAI_EACH = 2.0_wp
      ! forcing cfg + the diurnal NetCDF are still in place from block 5.
      call build_fast_context(cfg, ctx)                          ! ctx WITH the RT optics table (rad_opt)
      cfg%pft%wstress_psi_open  = -50.0_wp                        ! neutralize hydraulic stress (beta==1),
      cfg%pft%wstress_psi_close = -100.0_wp                       !   so the LIGHT/gb ordering drives GPP
      call init_bare_ground(poly, cfg, 1_ik)
      call use_test_column(poly)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 1.0_wp, 45.0_wp)     ! tall overstory  (big dbh -> tall)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 1.0_wp,  6.0_wp)     ! short understory (shaded below it)
      call finalize_init(poly)
      !----- equalize per-cohort LAI (leaf_area is per-plant, set from dbh, independent of nplant),  !
      !      so the RT's vertical light gradient is the ONLY thing separating the two per-leaf PARs.  !
      do ii = 1_ik, poly%cohort%n
         poly%cohort%nplant(ii) = LAI_EACH / max(poly%cohort%leaf_area(ii), 1.0e-9_wp)
      end do
      call met_open(drv, cfg%forcing)
      call met_cursor_init(drv, cur, 1_ik, cfg%forcing%latitude_deg,      &
                           cfg%forcing%longitude_deg,                            &
                           cfg%forcing%elevation_m)

      call init_fast_reservoirs(poly, ctx)
      call fast_dynamics(poly, ctx, cfg, met_src=drv, met_cur=cur,               &
                               step_start=meds_time_t(2020_ik,7_ik,1_ik,15_ik))   ! 15-17 UTC (day)
      call met_close(drv)

      itall = 1_ik ; ishort = 1_ik                               ! locate the tallest / shortest cohort
      do ii = 2_ik, poly%cohort%n
         if (poly%cohort%height(ii) > poly%cohort%height(itall))  itall  = ii
         if (poly%cohort%height(ii) < poly%cohort%height(ishort)) ishort = ii
      end do
      gpp_la_top   = poly%cohort%gpp_accum(itall)  / max(poly%cohort%leaf_area(itall),  1.0e-9_wp)
      gpp_la_under = poly%cohort%gpp_accum(ishort) / max(poly%cohort%leaf_area(ishort), 1.0e-9_wp)

      call check(itall /= ishort, 'two cohorts of distinct height were seeded (overstory + understory)')
      call check(gpp_la_top > gpp_la_under,                                                        &
                 'canopy-RT shading: top cohort gets more PAR/leaf -> higher GPP/leaf than understory')
      write(*,'(a,es10.3,a,es10.3,a)') '   (day GPP per leaf area: overstory=', gpp_la_top,        &
                                       '  understory=', gpp_la_under, ' kgC/m2)'
   end block

   !=== 7. MULTI-PATCH forcing: the per-patch forcing buffers must RESIZE to each patch's cohort   !
   !    count. A site whose FIRST-processed patch has FEWER cohorts than a later one (patch1=1,     !
   !    patch2=3; finalize_init keeps patch insertion order) would, with an allocate-once buffer,   !
   !    write patch 2's cohorts out of bounds. Under the Debug (-check all / -Mbounds) build this   !
   !    traps; here we assert every cohort gets a finite positive day GPP. The trailing BARE patch  !
   !    (patch3, 0 cohorts) also exercises apply_rt_forcing's empty-canopy path under forcing.  ====!
   block
      type(met_source_t) :: drv
      type(met_cursor_t) :: cur
      integer(ik)        :: ii, ncoh_all
      real(wp)           :: gpp_min
      call build_fast_context(cfg, ctx)                          ! ctx WITH the RT optics table
      call init_bare_ground(poly, cfg, 3_ik)                     ! THREE patches (patch3 stays bare)
      call use_test_column(poly)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.20_wp, 20.0_wp)   ! patch 1: 1 cohort (processed first)
      call add_cohort(poly, cfg, 2_ik, 1_ik, 0.20_wp, 25.0_wp)   ! patch 2: 3 cohorts (larger -> would overflow)
      call add_cohort(poly, cfg, 2_ik, 1_ik, 0.20_wp, 12.0_wp)
      call add_cohort(poly, cfg, 2_ik, 1_ik, 0.20_wp,  5.0_wp)
      call finalize_init(poly)                                   ! patch 3 keeps 0 cohorts (bare)
      call met_open(drv, cfg%forcing)
      call met_cursor_init(drv, cur, 1_ik, cfg%forcing%latitude_deg,      &
                           cfg%forcing%longitude_deg,                            &
                           cfg%forcing%elevation_m)
      call init_fast_reservoirs(poly, ctx)
      call fast_dynamics(poly, ctx, cfg, met_src=drv, met_cur=cur,               &
                               step_start=meds_time_t(2020_ik,7_ik,1_ik,15_ik))   ! 15-17 UTC (day)
      call met_close(drv)
      ncoh_all = poly%cohort%n
      gpp_min  = huge(1.0_wp)
      do ii = 1_ik, ncoh_all
         gpp_min = min(gpp_min, poly%cohort%gpp_accum(ii))
      end do
      call check(ncoh_all == 4_ik, 'multi-patch: all 4 cohorts (1 in patch1 + 3 in patch2) present')
      call check(gpp_min > 0.0_wp, 'multi-patch: every cohort got positive day GPP (forcing buffers resized)')
      write(*,'(a,i0,a,es10.3,a)') '   (multi-patch cohorts=', ncoh_all, '  min day GPP=', gpp_min, ' kgC/plant)'
   end block

   !=== 8. SEAM CONTINUITY: a forced leaf-carbon collapse between two fast_dynamics calls (an        !
   !    operator that drops leaf carbon without its water) must not reach the fast loop as a          !
   !    supersaturated store: settle_tissue_water, which the stepper runs before the fast loop, sheds  !
   !    the water above the new ceiling to the ground, and the fast loop then closes its budgets. ====!
   block
      real(wp) :: mass_before, leaf_carbon_before

      cfg%forcing%forcing_on = .false.     ! plain constant-forcing path (isolate from blocks 5-7)
      call init_bare_ground(poly, cfg, 1_ik)
      call use_test_column(poly)
      call add_cohort(poly, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
      call finalize_init(poly)
      call init_fast_reservoirs(poly, ctx)

      !----- One normal call evolves the birth water (the round trip block 1 exercises). ----------!
      call fast_dynamics(poly, ctx, cfg, worst_energy=we, worst_water=ww, n_budget_fail=n_budget_fail)
      mass_before        = poly%cohort%leaf_water_mass(1)
      leaf_carbon_before = poly%cohort%leaf_carbon(1)
      call check(mass_before > 0.0_wp, 'seam test: the leaves hold water before the forced collapse')

      !----- Force a snap-to-bare-sized leaf-carbon collapse directly on the cohort: bleaf drops to  !
      !      1% of its former value, so the OLD mass is now ~100x above the new ceiling               !
      !      water_sat*bleaf -- exactly the state a real dormant-canopy snap would leave behind.  ----!
      poly%cohort%leaf_carbon(1) = leaf_carbon_before * 0.01_wp
      call hold_cohort_tissue_temps(poly%cohort, 1_ik)   ! at the leaf's temperature, as an operator would

      call settle_tissue_water(poly, cfg)
      call check_close(poly%cohort%leaf_water_mass(1),                                              &
                       ctx%col_config%hydraulics_table%pft(1)%leaf_curve%water_sat                  &
                       * poly%cohort%leaf_carbon(1), 1.0e-12_wp,                                    &
                       'seam test: the settle capped the leaf water at its new saturation')
      call check(poly%patch%shed_water_rate(1) > 0.0_wp, 'seam test: the excess went to the ground')
      call fast_dynamics(poly, ctx, cfg, worst_energy=we, worst_water=ww, n_budget_fail=n_budget_fail)
      call check(n_budget_fail == 0_ik, 'seam test: whole-column budgets still close after the forced collapse')
      call check(poly%cohort%leaf_water_mass(1) < 0.5_wp*mass_before,                               &
                 'seam test: leaf water mass dropped sharply (the excess was shed), not left at the old value')
      write(*,'(a,es10.3,a,es10.3,a)') '   (seam clamp: leaf water mass ', mass_before, ' -> ',      &
                                       poly%cohort%leaf_water_mass(1), ' kg/plant after a 100x bleaf collapse)'

      !----- A BARE cohort (leaf_carbon = 0 -> leaf_water_mass = 0, its physical state) keeps its    !
      !      wood water: nothing re-seeds a store, so a dormant cohort carries its deficit through    !
      !      the winter. Perturb the wood store to a value a seed cannot produce, run a              !
      !      day with no leaves (no transpiration, so the wood store barely moves), and assert it    !
      !      carried through. ----------------------------------------------------------------------!
      block
         real(wp)    :: wood_set, wood_after, dt_fast_save
         integer(ik) :: nsub_save
         poly%cohort%leaf_carbon(1)     = 0.0_wp
         poly%cohort%leaf_water_mass(1) = 0.0_wp
         wood_set = 0.98_wp * poly%cohort%wood_water_mass(1)  ! a 2% deficit a seed cannot produce
         poly%cohort%wood_water_mass(1) = wood_set
         call hold_cohort_tissue_temps(poly%cohort, 1_ik)
         !----- ONE short sub-step with a SMALL deficit: near saturation the PV curve gives a small    !
         !      driving potential, so root uptake moves the store by << 0.5% in 60 s (a 30% deficit,   !
         !      by contrast, refills a third of the way in that minute). A re-seed from the soil puts  !
         !      the store back to 100% -- a 2% jump this check sees. ----------------------------------!
         !----- ...and with the plant-hydraulics process MASKED (frozen tissue water), because the    !
         !      fixture's leafless cohort refills even a 2% wood deficit within the minute: the store   !
         !      must come out exactly as it went in, so the only thing that can move it is the seed. --!
         dt_fast_save = cfg%dt_fast ; nsub_save = cfg%n_fast_per_slow
         cfg%dt_fast = 60.0_wp ; cfg%n_fast_per_slow = 1_ik
         ctx%col_config%mask%hydraulics = .false.        ! the fast loop reads the mask from the context, not cfg
         call fast_dynamics(poly, ctx, cfg, worst_energy=we, worst_water=ww, n_budget_fail=n_budget_fail)
         cfg%dt_fast = dt_fast_save ; cfg%n_fast_per_slow = nsub_save
         ctx%col_config%mask%hydraulics = .true.
         wood_after = poly%cohort%wood_water_mass(1)
         write(*,'(a,es10.3,a,es10.3,a)') '   (dormant seam: wood water set ', wood_set, ' -> ', wood_after, &
               ' kg/plant after one 60 s step)'
         call check(abs(wood_after - wood_set) < 1.0e-9_wp * wood_set,                               &
                    'dormant seam: bare cohort keeps its integrated wood water (not re-seeded)')
         call check(poly%cohort%leaf_water_mass(1) == 0.0_wp,                                        &
                    'dormant seam: a bare cohort holds no leaf water')
      end block
   end block

   !=== 9. SPLIT INVARIANCE and the LEAFLESS LEAF (docs/science/vegetation_energy_dynamics.md, "A      !
   !    thin tissue is a small tissue" and "A tissue with no area"). A thin stand -- 20 stems per      !
   !    hectare of a 20 cm tree, wood area index ~0.01 -- as one cohort and as two cohorts at half the !
   !    density, through the whole fast loop (two-stream, aerodynamics, water network, ARK) over a     !
   !    forced day window: the tissue temperatures, the canopy air and the per-plant GPP must agree.   !
   !    A floor in ground units on a tissue's coupling or heat capacity, which only the half-density   !
   !    pair would cross, fails this. The step controller is held off (rtol_all and atol_scale large)  !
   !    so both stands take the same steps: its error norm (docs/science/numerical_scheme.md eq. (3))  !
   !    counts the pair's states twice, each against the same atol, so under a stiff ground coupling   !
   !    the two sub-step differently, by ~1e-3 K -- the integrator's error, not the physics'. Then the !
   !    same tree without leaves: its leaf temperature is its wood's (a bud sits on its twig), and the !
   !    budgets close. ================================================================================!
   block
      type(met_source_t) :: drv
      type(met_cursor_t) :: cur
      integer(ik)        :: icase, nfail(3)
      real(wp)           :: t_leaf(3), t_wood(3), t_cas(3), gpp(3), d_split
      real(wp), parameter :: NPLANT_THIN = 0.002_wp, DBH_THIN = 20.0_wp
      real(wp), parameter :: HOLD_OFF = 1.0e6_wp     ! tolerances no step error reaches
      real(wp)           :: rtol_all0, atol_scale0
      cfg%forcing%forcing_on = .true.
      rtol_all0 = cfg%rtol_all ; atol_scale0 = cfg%atol_scale
      cfg%rtol_all = HOLD_OFF ; cfg%atol_scale = HOLD_OFF
      call build_fast_context(cfg, ctx)
      do icase = 1_ik, 3_ik
         call init_bare_ground(poly, cfg, 1_ik)
         call use_test_column(poly)
         if (icase == 2_ik) then
            call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp*NPLANT_THIN, DBH_THIN)
            call add_cohort(poly, cfg, 1_ik, 1_ik, 0.5_wp*NPLANT_THIN, DBH_THIN)
         else
            call add_cohort(poly, cfg, 1_ik, 1_ik, NPLANT_THIN, DBH_THIN)
         end if
         call finalize_init(poly)
         if (icase == 3_ik) then                     ! leafless: no leaf carbon, so no leaf area
            poly%cohort%leaf_carbon(1) = 0.0_wp
            poly%cohort%leaf_water_mass(1) = 0.0_wp
            call set_cohort_size_from_carbon(poly%cohort, 1_ik)
            call hold_cohort_tissue_temps(poly%cohort, 1_ik)
         end if
         call met_open(drv, cfg%forcing)
         call met_cursor_init(drv, cur, 1_ik, cfg%forcing%latitude_deg, cfg%forcing%longitude_deg,        &
                              cfg%forcing%elevation_m)
         call init_fast_reservoirs(poly, ctx)
         call fast_dynamics(poly, ctx, cfg, met_src=drv, met_cur=cur, step_start=meds_time_t(2020_ik,7_ik,1_ik,15_ik), &
                            n_budget_fail=nfail(icase))
         call met_close(drv)
         t_leaf(icase) = poly%cohort%leaf_temp(poly%cohort%n) ; t_wood(icase) = poly%cohort%wood_temp(poly%cohort%n)
         t_cas(icase)  = poly%patch%cas(1)%can_temp          ; gpp(icase)    = poly%cohort%gpp_accum(poly%cohort%n)
         if (icase == 2_ik) call check(poly%cohort%leaf_temp(1) == poly%cohort%leaf_temp(2) .and.              &
                                       poly%cohort%wood_temp(1) == poly%cohort%wood_temp(2),                   &
                                       'split: the two half-density cohorts are identical')
      end do
      d_split = max(abs(t_leaf(2) - t_leaf(1)), abs(t_wood(2) - t_wood(1)), abs(t_cas(2) - t_cas(1)))
      write(*,'(a,es10.3,a,es10.3,a,f8.3,a,f8.3,a)') '   (split: worst temperature difference ', d_split,         &
            ' K, per-plant GPP ', abs(gpp(2) - gpp(1)) / gpp(1), ' relative; wood ', t_wood(1), ' K, canopy air ', &
            t_cas(1), ' K)'
      call check(d_split < 1.0e-6_wp, 'split: one cohort and two at half density reach the same temperatures')
      call check(abs(gpp(2) - gpp(1)) <= 1.0e-6_wp * gpp(1), 'split: the same GPP per plant')
      call check(all(nfail(1:3) == 0_ik), 'split and leafless: the whole-column budgets close')
      call check(poly%cohort%leaf_area(1) == 0.0_wp, 'leafless: the cohort has no leaf area')
      call check(t_leaf(3) == t_wood(3), 'leafless: the leaf takes its wood''s temperature')
      call check(abs(t_wood(3) - t_cas(3)) > 1.0e-3_wp, 'leafless: the wood is not at the canopy air''s temperature')
      cfg%rtol_all = rtol_all0 ; cfg%atol_scale = atol_scale0
   end block

   write(*,'(a)')          '   PASS'
   write(*,'(a,f7.2,a,f7.2,a)') '   (patch-1 CAS temp ', t_cas0, ' -> ', t_cas1, ' K)'
   write(*,'(a,es10.3,a,es10.3,a)') '   (worst whole-column resid: energy=', we, ' J/m2  water=', ww, ' kg/m2)'

contains

   !----- The test's hand-built soil column is the site's, which the fast loop runs on, starting  !
   !      at 0.30 in every layer. ------------------------------------------------------------!
   subroutine use_test_column(p)
      type(polygon_t), intent(inout) :: p
      p%site%surface(1)%soil         = ctx%col_config%soil
      p%site%surface(1)%soil_thermal = ctx%col_config%soil_thermal
      p%site%surface(1)%theta_start  = 0.30_wp
   end subroutine use_test_column

   !----- Write a 1-day hourly (time=25, grid=1) MEDS forcing NetCDF with a diurnal shortwave    !
   !      (peaks ~17 UTC at Ithaca) so the fast loop sees a real day/night cycle.                 !
   subroutine write_diurnal_forcing(path)
      character(len=*), intent(in) :: path
      integer, parameter :: NT = 25
      integer(c_int) :: st, ncid, td, gd, vt, vla, vlo, vv(7), dims2(2), dims1(1)
      integer(c_size_t) :: s1(1), c1(1), s2(2), c2(2)
      real(c_double) :: tsec(NT), dat(1, NT)
      integer :: it, k
      real(wp) :: hh, sw
      character(len=8), parameter :: nm(7) = ['Tair    ','Qair    ','PSurf   ','Wind    ',        &
                                              'Rainf   ','SWdown  ','LWdown  ']
      st = nc_create_f(path, NC_NETCDF4, ncid) ; call nc_check(st,'c')
      st = nc_def_dim_f(ncid,'time',int(NT,c_size_t),td) ; call nc_check(st,'t')
      st = nc_def_dim_f(ncid,'grid',1_c_size_t,gd) ; call nc_check(st,'g')
      dims1(1)=td ; st = nc_def_var_f(ncid,'time',NC_DOUBLE,1,dims1,vt) ; call nc_check(st,'tv')
      st = nc_put_att_text_f(ncid,vt,'units',int(len_trim('seconds since 2020-07-01 00:00:00'),c_size_t), &
                             'seconds since 2020-07-01 00:00:00') ; call nc_check(st,'u')
      dims1(1)=gd ; st = nc_def_var_f(ncid,'latitude',NC_DOUBLE,1,dims1,vla) ; call nc_check(st,'la')
      st = nc_def_var_f(ncid,'longitude',NC_DOUBLE,1,dims1,vlo) ; call nc_check(st,'lo')
      dims2(1)=td ; dims2(2)=gd
      do k=1,7 ; st = nc_def_var_f(ncid,trim(nm(k)),NC_DOUBLE,2,dims2,vv(k)) ; call nc_check(st,'v') ; end do
      st = nc_put_att_text_f(ncid,NC_GLOBAL,'time_zone',3_c_size_t,'UTC') ; call nc_check(st,'tz')
      st = nc_enddef(ncid) ; call nc_check(st,'e')
      do it=1,NT ; tsec(it)=real(it-1,c_double)*3600.0_c_double ; end do
      s1=0_c_size_t ; c1(1)=int(NT,c_size_t) ; st=nc_put_vara_double(ncid,vt,s1,c1,tsec) ; call nc_check(st,'tw')
      c1(1)=1_c_size_t
      block ; real(c_double) :: la(1),lo(1) ; la=42.44_c_double ; lo=-76.50_c_double
         st=nc_put_vara_double(ncid,vla,s1,c1,la) ; call nc_check(st,'law')
         st=nc_put_vara_double(ncid,vlo,s1,c1,lo) ; call nc_check(st,'low') ; end block
      s2=0_c_size_t ; c2=[int(NT,c_size_t),1_c_size_t]
      do k=1,7
         do it=1,NT
            hh=real(it-1,wp)
            sw = max(0.0_wp, 850.0_wp*sin(3.14159265_wp*(hh-11.0_wp)/12.0_wp))
            if (hh<11.0_wp .or. hh>23.0_wp) sw=0.0_wp
            select case (k)
            case(1) ; dat(1,it)=290.0_wp + 6.0_wp*sin(2.0_wp*3.14159265_wp*(hh-15.0_wp)/24.0_wp)
            case(2) ; dat(1,it)=0.010_wp
            case(3) ; dat(1,it)=97500.0_wp
            case(4) ; dat(1,it)=2.5_wp
            case(5) ; dat(1,it)=0.0_wp
            case(6) ; dat(1,it)=sw
            case(7) ; dat(1,it)=340.0_wp
            end select
         end do
         st=nc_put_vara_double(ncid,vv(k),s2,c2,dat) ; call nc_check(st,'vw')
      end do
      st = nc_close(ncid) ; call nc_check(st,'cl')
   end subroutine write_diurnal_forcing

end program test_fast_loop
