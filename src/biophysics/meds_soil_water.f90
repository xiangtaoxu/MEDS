!==========================================================================================!
! meds_soil_water -- THE stateless seam of the 1-D soil-water column (design                   !
! MEDS_COLUMN_HYDROLOGY_DESIGN.md). It advances one patch's prognostic soil moisture + ponding !
! (+ a lumped aquifer) over a fast step dt: canopy-throughfall infiltration (conductivity-       !
! limited), ground evaporation, the multi-layer soil-water balance (implicit backward-Euler       !
! Thomas), and a water-potential-limited root-uptake sink; it returns the boundary fluxes plus     !
! the per-layer matric potential psi_soil that closes the plant-hydraulics boundary condition.      !
! Exposes the store in two forms (like meds_soil_energy):                                            !
!   * soil_water_time_deriv    -- the EXPLICIT Richards RHS dtheta_k/dt [1/s] (for the IMEX-ARK       !
!                                 integrator; pure, commits nothing).                                 !
!   * soil_water_step_implicit -- one IMPLICIT backward-Euler / Celia-Picard sub-step over dt.         !
! The seam column_hydrology_flux orchestrates the adaptive substepping + surface BCs around the         !
! implicit step. Canopy interception lives in meds_vegetation_biophysics (a per-cohort film).            !
!==========================================================================================!
module meds_soil_water
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, grav, r_wv, tiny_num, grav_head, p_std
   use meds_biophysics_types, only : soil_column_t, chydro_forcing_t, soil_params_t,          &
                                     soil_opts_t, chydro_flux_t, n_soil_layer_max,             &
                                     SOIL_BC_BEDROCK, SOIL_BC_AQUIFER, SOIL_RETENTION_CAMPBELL, &
                                     SOIL_LIN_PICARD, SOIL_SUBSTEP_FIXED
   use meds_hydr_lib,         only : soil_psi_from_theta, soil_theta_from_psi, soil_hydr_cond_from_theta, &
                                     soil_moist_cap_from_psi
   use meds_numerics,         only : thomas_solve
   use meds_therm_lib,        only : sat_specific_humidity, internal_energy_liquid, uext_to_temp
   implicit none
   private

   public :: column_hydrology_flux, soil_water_time_deriv, soil_water_step_implicit

   !----- Unconfined-aquifer specific yield (fraction) for the z_wt <-> w_aquifer diagnosis. -!
   real(wp), parameter :: SPECIFIC_YIELD = 0.2_wp
   !----- [kg/m2] below this the pond is treated as empty: it has no meaningful temperature, and any
   !      enthalpy residue left in a drained store would read as heat in an empty pond next step.
   real(wp), parameter :: POND_TINY = 1.0e-9_wp

contains

   !---------------------------------------------------------------------------------------!
   ! The seam: advance one patch soil column over dt. `col` (theta + stores) is the only      !
   ! mutated data; `flux` carries boundary fluxes, exported psi_soil, and the mass-budget       !
   ! residual (asserts to ~round-off).                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine column_hydrology_flux(col, forcing, params, opts, dt, flux)
      type(soil_column_t),    intent(inout) :: col
      type(chydro_forcing_t), intent(in)    :: forcing
      type(soil_params_t),    intent(in)    :: params
      type(soil_opts_t),      intent(in)    :: opts
      real(wp),               intent(in)    :: dt
      type(chydro_flux_t),    intent(out)   :: flux

      integer(ik) :: n, k, rc, nsub
      real(wp), dimension(n_soil_layer_max) :: theta0, theta1, psi_e, clip_l, floor_l
      real(wp) :: z_wt, z_bottom, psi1, kn1, e_soil, q_inf_max, q_avail, infl, q_top
      real(wp) :: f_sat, q_over, q_liq, drain_amt, uptake_amt, clip_ex, deficit, want, give
      real(wp) :: q_drai, recharge_amt, baseflow_amt, site_drain, wsurf, runoff, w0, w1
      real(wp) :: w_surf0, w_aqf0
      real(wp) :: e_surf0, esurf, t_pond, fliq_pond, over_mass   ! pond enthalpy (#78 item 4)
      real(wp) :: face_resid, f_in, f_out, f_sink
      logical  :: ok

      n  = params%n_active
      rc = params%retention
      z_bottom = params%soil_layer_z(n+1)                                 ! deepest interface (<= 0)
      theta0(1:n) = col%theta(1:n)
      w_surf0 = col%w_surface                                             ! initial stores (for the budget)
      e_surf0 = col%w_surface_enth                                        ! pond ENTHALPY at entry (#78 item 4)
      w_aqf0  = col%w_aquifer

      !----- Water-table elevation z_wt (<= 0): diagnosed for the aquifer BC, else the column   !
      !      bottom (deep -> ZD/Dunne inert), per the design.                                    !
      if (opts%bottom_bc == SOIL_BC_AQUIFER) then
         z_wt = col%z_wt
      else
         z_wt = z_bottom
      end if

      !----- Equilibrium potential for Zeng-Decker (retention-integral; else unused). -------!
      psi_e = 0.0_wp
      if (opts%zeng_decker) call compute_psi_e(params, rc, z_wt, n, psi_e)

      !----- Top-layer psi/K for evaporation + infiltration capacity. -----------------------!
      psi1 = soil_psi_from_theta(rc, theta0(1), params%theta_sat(1), params%theta_res(1),        &
                               curve_a(params,1), curve_n(params,1))
      kn1  = soil_hydr_cond_from_theta(rc, theta0(1), params%theta_sat(1), params%theta_res(1),            &
                            curve_a(params,1), curve_n(params,1), params%ksat(1))

      !----- Ground evaporation (Philip alpha_soil + Swenson-Lawrence DSL), capped. ---------!
      e_soil = ground_evaporation(theta0(1), psi1, params, forcing, opts)
      e_soil = min(e_soil, max(0.0_wp, (theta0(1) - params%theta_res(1)) * params%dz(1)         &
                               * rho_h2o / dt))

      !----- Dunne saturation-excess runoff, split off BEFORE infiltration (design 3d/3e). It is   !
      !      only meaningful with a GENUINELY DIAGNOSED water table (the SIMTOP aquifer BC). Under   !
      !      FREE_DRAIN/BEDROCK there is no water table -- z_wt is pinned at the geometric column     !
      !      bottom, so the old unconditional f_max*exp(..z_bottom) shed a fixed fraction of ALL rain !
      !      from an unsaturated (even bone-dry) column, worse for a shallow column (BUG7). Make it   !
      !      inert there: only Horton infiltration-excess (ponding overflow, below) runs off.  ------!
      if (opts%bottom_bc == SOIL_BC_AQUIFER) then
         f_sat = opts%f_max * exp(0.5_wp * opts%f_over * z_wt)             ! z_wt <= 0 => <= f_max
      else
         f_sat = 0.0_wp
      end if
      q_over = f_sat * forcing%precip_ground
      q_liq  = forcing%precip_ground - q_over

      !----- Conductivity-limited infiltration + ponding partition (design 3d). -------------!
      q_inf_max = kn1 * (1.0_wp + (-psi1) / (-params%z_node(1)))           ! Darcy velocity [m/s]
      q_avail   = col%w_surface / dt + q_liq                              ! [kg/m2/s]
      infl      = min(q_avail, q_inf_max * rho_h2o)
      q_top     = (infl - e_soil) / rho_h2o                              ! [m/s down], constant over substeps

      !----- Advance the interior (adaptive/fixed substepping of the implicit solve). --------!
      theta1(1:n) = theta0(1:n)
      call soil_water_advance(theta1, params, opts, rc, n, dt, q_top, psi_e,                   &
                              forcing%root_uptake, drain_amt, uptake_amt, flux%w_flux, nsub, ok)

      !----- FACE-FLUX CONSISTENCY (the contract the soil ENERGY column depends on). ---------!
      !      soil_energy_step_implicit advects liquid enthalpy on wflux_out, so those face fluxes    !
      !      must reproduce the mass that actually moved -- otherwise enthalpy is transported for      !
      !      water that did not move (or vice versa), and because internal_energy_liquid carries the   !
      !      tsupercool_liq datum (~4.2e5 J/kg in absolute terms) even a small inconsistency becomes    !
      !      a huge spurious heat flux. The column mass_resid below CANNOT see this: it checks the      !
      !      column against its BOUNDARY fluxes, and any interior face error cancels in that sum.       !
      !      Checked pre-clip, since the clip/give-back corrections are separate bookkeeping applied     !
      !      after the solve. --------------------------------------------------------------------------!
      face_resid = 0.0_wp
      do k = 1_ik, n
         if (k == 1_ik) then
            f_in = (infl - e_soil) * dt                                  ! top face [kg/m2]
         else
            f_in = flux%w_flux(k-1_ik) * rho_h2o * dt
         end if
         if (k == n) then
            f_out = drain_amt                                            ! bottom face [kg/m2]
         else
            f_out = flux%w_flux(k) * rho_h2o * dt
         end if
         f_sink = forcing%root_uptake(k) * dt
         face_resid = face_resid + abs((theta1(k) - theta0(k)) * params%dz(k) * rho_h2o           &
                                       - (f_in - f_out - f_sink))
      end do

      !----- Bottom boundary: aquifer store + baseflow, or direct site drainage. -------------!
      if (opts%bottom_bc == SOIL_BC_AQUIFER) then
         recharge_amt = drain_amt                                         ! soil -> aquifer [kg/m2]
         q_drai       = opts%q_drai_max * exp(opts%f_drai * z_wt)         ! baseflow [kg/m2/s], z_wt<=0
         baseflow_amt = q_drai * dt
         col%w_aquifer = max(0.0_wp, col%w_aquifer + recharge_amt - baseflow_amt)
         col%z_wt      = min(0.0_wp, z_bottom + col%w_aquifer / (rho_h2o * SPECIFIC_YIELD))
         site_drain    = baseflow_amt / dt                               ! only baseflow LEAVES the site
      else
         site_drain    = drain_amt / dt                                  ! free-drain/bedrock: leaves directly
      end if

      !----- Post-solve clip + theta_wp cap (design 5.1 step 5), all bookkept. --------------!
      !      Each correction below moves water with NO face, so the soil ENERGY column (which advects   !
      !      enthalpy on the faces but consumes the CORRECTED theta) cannot see it. The two that are     !
      !      genuine boundary crossings with a well-defined temperature -- the saturation clip (water     !
      !      leaves layer k for the pond) and the theta_res hard floor (water is created in layer k) --   !
      !      are recorded PER LAYER so the caller can move the matching enthalpy at layer k's own         !
      !      temperature. The two GIVE-BACKs are deliberately NOT recorded: they undo part of the root    !
      !      sink, whose enthalpy counterpart is keyed to transpiration at LEAF temperature rather than   !
      !      to the realized per-layer mass, so compensating them belongs with that deferred rework       !
      !      (see meds_fast_split.f90's root_heat_sink note) and cannot be done here alone. --------------!
      clip_ex = 0.0_wp
      deficit = 0.0_wp
      clip_l(1:n)  = 0.0_wp
      floor_l(1:n) = 0.0_wp
      do k = 1_ik, n
         if (theta1(k) > params%theta_sat(k)) then
            clip_l(k) = (theta1(k) - params%theta_sat(k)) * params%dz(k) * rho_h2o
            clip_ex   = clip_ex + clip_l(k)
            theta1(k) = params%theta_sat(k)
         end if
         if (theta1(k) < params%theta_wp(k) .and. forcing%root_uptake(k) > 0.0_wp) then
            want      = (params%theta_wp(k) - theta1(k)) * params%dz(k) * rho_h2o
            give      = min(want, max(0.0_wp, uptake_amt - deficit))     ! total give back <= total extracted
            theta1(k) = theta1(k) + give / (params%dz(k) * rho_h2o)
            deficit   = deficit + give
         end if
         !----- Residual floor: theta cannot fall below theta_res (the constitutive curves need    !
         !      Se >= 0). A numerical undershoot is water the solver over-removed; SOURCE the        !
         !      correction from the site's extracted root water (deficit), same as the give-back     !
         !      above, so it is BOOKKEPT (reduces uptake_total) rather than silently created. Any     !
         !      remainder once the extracted water is exhausted is hard-floored below and shows up in !
         !      the mass-budget residual (Debug error-stop).  --------------------------------------!
         if (theta1(k) < params%theta_res(k)) then
            want      = (params%theta_res(k) - theta1(k)) * params%dz(k) * rho_h2o
            give      = min(want, max(0.0_wp, uptake_amt - deficit))
            theta1(k) = theta1(k) + give / (params%dz(k) * rho_h2o)
            deficit   = deficit + give
         end if
         floor_l(k) = max(0.0_wp, (params%theta_res(k) - theta1(k)) * params%dz(k) * rho_h2o)
         theta1(k) = max(theta1(k), params%theta_res(k))                 ! hard residual floor (last resort)
      end do

      !----- Commit soil moisture; export psi_soil [MPa]. ----------------------------------!
      col%theta(1:n) = theta1(1:n)
      do k = 1_ik, n
         flux%psi_soil(k) = grav_head * soil_psi_from_theta(rc, theta1(k), params%theta_sat(k),  &
                            params%theta_res(k), curve_a(params,k), curve_n(params,k))
      end do

      !----- Ponding store: PAIRED mass + enthalpy, in the order the transfers physically happen.      !
      !                                                                                                !
      !      The pond used to be mass-only, so water crossing into it shed its enthalpy and that energy  !
      !      left the whole-column ledger while the water sat on the surface still holding its heat      !
      !      (issue #78 item 4). It now carries col%w_surface_enth, and every seam below moves mass and  !
      !      enthalpy together. Enthalpy is owned HERE, with the mass, rather than in the callers: a     !
      !      store whose two halves are updated in different places drifts, which is exactly how the     !
      !      ARK condensate deposit and the RK45 double-clip (both fixed in PR #81) happened.            !
      !                                                                                                !
      !      ORDER MATTERS and follows the physics of the step:                                         !
      !        1. rain joins the pond at t_precip -- it arrives before anything is drawn from it;        !
      !        2. infiltration draws from that MIXTURE, so it leaves at the mixed temperature. This is   !
      !           why the soil's top-face advection must use flux%t_infil and not rain_temp: the water   !
      !           entering layer 1 came out of the pond, not out of the sky. With a dry pond the mixture !
      !           IS the rain, so the common case is unchanged;                                          !
      !        3. the saturation clip joins afterwards, at each layer's OWN temperature -- it is water   !
      !           the soil rejected at the END of the solve, so it cannot have been available to         !
      !           infiltrate earlier in the same step;                                                   !
      !        4. overflow leaves at the final pond temperature.                                         !
      !                                                                                                  !
      !      Temperature is uext_to_temp with dry_hcap = 0, the same read-off the snow pack uses, so a    !
      !      freezing pond hits the melt plateau instead of going unphysically cold. An EMPTY pond has no !
      !      temperature to speak of: guard on the mass and fall back to t_precip, so a dry column        !
      !      reproduces the old behaviour exactly. -------------------------------------------------------!
      wsurf = w_surf0 + q_liq * dt
      esurf = e_surf0 + q_liq * dt * internal_energy_liquid(forcing%t_precip)
      !----- 2. infiltration leaves at the MIXED pond temperature. -----------------------------------!
      t_pond = forcing%t_precip
      if (wsurf > POND_TINY) call uext_to_temp(esurf, wsurf, 0.0_wp, t_pond, fliq_pond)
      flux%t_infil = t_pond
      wsurf = wsurf - infl * dt
      esurf = esurf - infl * dt * internal_energy_liquid(t_pond)
      !----- 3. the clip joins, each layer's water at that layer's temperature. ----------------------!
      do k = 1_ik, n
         esurf = esurf + clip_l(k) * internal_energy_liquid(forcing%soil_temp(k))
      end do
      wsurf = wsurf + clip_ex
      !----- 4. overflow (Horton) at the final pond temperature, plus the Dunne share, which never     !
      !      entered the pond and so leaves at t_precip. ---------------------------------------------!
      t_pond = forcing%t_precip
      if (wsurf > POND_TINY) call uext_to_temp(esurf, wsurf, 0.0_wp, t_pond, fliq_pond)
      over_mass = max(0.0_wp, wsurf - opts%w_pond_max)
      runoff = q_over + over_mass / dt
      flux%runoff_enth = q_over * internal_energy_liquid(forcing%t_precip)                            &
                       + over_mass / dt * internal_energy_liquid(t_pond)
      wsurf  = wsurf - over_mass
      esurf  = esurf - over_mass * internal_energy_liquid(t_pond)
      !----- A pond that has drained to nothing must carry no enthalpy either, or the residue reads as !
      !      heat in an empty store on the next step. -------------------------------------------------!
      if (wsurf <= POND_TINY) then
         wsurf = max(0.0_wp, wsurf) ; esurf = 0.0_wp
      end if
      col%w_surface      = wsurf
      col%w_surface_enth = esurf

      !----- Mass-budget closure check (total stores vs boundary fluxes). -------------------!
      w0 = w_surf0 + w_aqf0
      w1 = col%w_surface + col%w_aquifer
      do k = 1_ik, n
         w0 = w0 + theta0(k) * params%dz(k) * rho_h2o
         w1 = w1 + theta1(k) * params%dz(k) * rho_h2o
      end do

      flux%infiltration   = infl
      flux%drainage       = site_drain
      flux%runoff_surf    = runoff
      flux%soil_evap      = e_soil
      flux%uptake_total   = max(0.0_wp, uptake_amt - deficit) / dt
      flux%uptake_deficit = deficit / dt
      flux%clip_excess    = clip_ex / dt
      flux%clip_layer(1:n)  = clip_l(1:n)  / dt
      flux%floor_layer(1:n) = floor_l(1:n) / dt
      flux%mass_resid     = (w1 - w0) - dt * (forcing%precip_ground - e_soil - site_drain       &
                            - flux%uptake_total - runoff)
      flux%face_mass_resid = face_resid
      flux%nsub      = nsub
      flux%converged = ok
      if (opts%debug_error .and. abs(flux%mass_resid) > opts%atol) then
         error stop 'column_hydrology_flux: mass budget did not close'
      end if
      if (opts%debug_error .and. face_resid > opts%atol) then
         error stop 'column_hydrology_flux: per-face mass budget did not close (w_flux inconsistent)'
      end if
   end subroutine column_hydrology_flux

   !=======================================================================================!
   !  Interior solver: adaptive step-doubling wrapper + one implicit (BE/Picard) sub-step.  !
   !=======================================================================================!

   !----- Advance theta over dt with substepping; accumulate drainage + uptake amounts       !
   !      [kg/m2 over dt] and the sub-step count. q_top (top-face Darcy velocity) is held      !
   !      constant across the step; psi_e is the frozen Zeng-Decker reference.                 !
   !---------------------------------------------------------------------------------------!
   subroutine soil_water_advance(theta, params, opts, rc, n, dt, q_top, psi_e, root_uptake,   &
                                 drainage_tot, uptake_tot, wflux_out, nsub, ok)
      real(wp),            intent(inout) :: theta(n_soil_layer_max)
      type(soil_params_t), intent(in)    :: params
      type(soil_opts_t),   intent(in)    :: opts
      integer(ik),         intent(in)    :: rc, n
      real(wp),            intent(in)    :: dt, q_top, psi_e(n_soil_layer_max)
      real(wp),            intent(in)    :: root_uptake(n_soil_layer_max)
      real(wp),            intent(out)   :: drainage_tot, uptake_tot
      real(wp),            intent(out)   :: wflux_out(n_soil_layer_max)   !< [m/s] time-mean face flux (k=1..n-1)
      integer(ik),         intent(out)   :: nsub
      logical,             intent(out)   :: ok

      real(wp), dimension(n_soil_layer_max) :: th_big, th_h1, th_two, wf_b, wf_1, wf_2, wface_tot
      real(wp)    :: h, t, hmin, err, dr_b, up_b, dr_1, up_1, dr_2, up_2
      integer(ik) :: k, nfix
      real(wp), parameter :: safety = 0.9_wp, fmin = 0.25_wp, fmax = 4.0_wp
      logical :: dum

      drainage_tot = 0.0_wp
      uptake_tot = 0.0_wp
      nsub = 0_ik
      ok = .true.
      wface_tot = 0.0_wp

      if (opts%substep == SOIL_SUBSTEP_FIXED) then
         !----- Warp-uniform fixed count (GPU-friendly). ----------------------------------!
         nfix = max(1_ik, min(opts%max_substep, nint(dt / max(opts%h_init, tiny_num), ik)))
         h    = dt / real(nfix, wp)
         do k = 1_ik, nfix
            call soil_water_step_implicit(theta, params, opts, rc, n, h, q_top, psi_e, root_uptake, &
                                          th_big, dr_b, up_b, wf_b, dum)
            ok = ok .and. dum                          ! aggregate inner convergence (was discarded)
            theta(1:n) = th_big(1:n)
            drainage_tot = drainage_tot + dr_b
            uptake_tot = uptake_tot + up_b
            wface_tot(1:n) = wface_tot(1:n) + wf_b(1:n)
         end do
         nsub = nfix
         wflux_out(1:n) = wface_tot(1:n) / dt
         return
      end if

      !----- Adaptive step-doubling (BE is L-stable; substep only for accuracy). ------------!
      t = 0.0_wp
      h = min(opts%h_init, dt)
      hmin = dt / real(opts%max_substep, wp)
      do while (t < dt - 1.0e-9_wp * dt .and. nsub < opts%max_substep)
         h = min(h, dt - t)
         call soil_water_step_implicit(theta, params, opts, rc, n, h,        q_top, psi_e,       &
                                       root_uptake, th_big, dr_b, up_b, wf_b, dum)
         call soil_water_step_implicit(theta, params, opts, rc, n, 0.5_wp*h, q_top, psi_e,       &
                                       root_uptake, th_h1, dr_1, up_1, wf_1, dum)
         call soil_water_step_implicit(th_h1, params, opts, rc, n, 0.5_wp*h, q_top, psi_e,       &
                                       root_uptake, th_two, dr_2, up_2, wf_2, dum)
         err = 1.0e-12_wp
         do k = 1_ik, n
            err = max(err, abs(th_two(k) - th_big(k)) / (opts%atol + opts%rtol * abs(th_two(k))))
         end do
         !----- NON-FINITE GUARD. Both loop-exit variables (t, nsub) advance ONLY on the accept branch,  !
         !      so a trial that can never be accepted spins forever. That is exactly what a non-finite   !
         !      err does: NaN fails `err <= 1`, and once the reject branch turns h into NaN it fails      !
         !      `h <= hmin` too -- an infinite loop that hangs the whole model inside one fast step, with !
         !      no error and no output (diagnosed from a 30-yr run that froze mid-day at 100% CPU). Bail  !
         !      out flagged unconverged instead, leaving theta at the last good value, so the caller sees !
         !      ok = .false. rather than the model wedging. -----------------------------------------------!
         if (err /= err .or. h /= h .or. any(th_two(1:n) /= th_two(1:n))) then
            ok = .false. ; exit
         end if
         if (err <= 1.0_wp .or. h <= hmin * 1.0001_wp) then
            theta(1:n)   = th_two(1:n)                         ! the more accurate (two half-steps)
            drainage_tot = drainage_tot + dr_1 + dr_2
            uptake_tot   = uptake_tot + up_1 + up_2
            wface_tot(1:n) = wface_tot(1:n) + wf_1(1:n) + wf_2(1:n)
            t    = t + h
            nsub = nsub + 1_ik
            h    = h * min(fmax, safety * err ** (-0.5_wp))
         else
            h    = h * max(fmin, safety * err ** (-0.5_wp))    ! reject, shrink, retry
         end if
      end do
      if (t < dt - 1.0e-9_wp * dt) ok = .false.               ! ran out of sub-steps
      wflux_out(1:n) = wface_tot(1:n) / dt
   end subroutine soil_water_advance

   !----- One implicit backward-Euler sub-step of length h. Frozen-coefficient (1 iterate) or  !
   !      Celia (1990) modified-Picard (up to max_picard). Upstream K; retention-integral       !
   !      Zeng-Decker gravity via psi_e; conservative flux-divergence theta update. Returns      !
   !      drainage + uptake AMOUNTS [kg/m2 over h].                                              !
   !---------------------------------------------------------------------------------------!
   subroutine soil_water_step_implicit(theta_in, params, opts, rc, n, h, q_top, psi_e, root_uptake, &
                                       theta_out, drainage_amt, uptake_amt, wface_amt, ok)
      real(wp),            intent(in)  :: theta_in(n_soil_layer_max)
      type(soil_params_t), intent(in)  :: params
      type(soil_opts_t),   intent(in)  :: opts
      integer(ik),         intent(in)  :: rc, n
      real(wp),            intent(in)  :: h, q_top, psi_e(n_soil_layer_max)
      real(wp),            intent(in)  :: root_uptake(n_soil_layer_max)
      real(wp),            intent(out) :: theta_out(n_soil_layer_max)
      real(wp),            intent(out) :: drainage_amt, uptake_amt
      real(wp),            intent(out) :: wface_amt(n_soil_layer_max)   !< [m] downward water depth across face k
      logical,             intent(out) :: ok

      real(wp), dimension(n_soil_layer_max) :: psi_m, theta_m, kk, cc, kface, gface, qface
      real(wp), dimension(n_soil_layer_max) :: a, b, c, rhs, dpsi, sk, dsk
      integer(ik) :: k, iter, maxit
      real(wp)    :: qbot, in_k, out_k, err, theta_prev

      theta_m(1:n) = theta_in(1:n)
      do k = 1_ik, n
         psi_m(k) = soil_psi_from_theta(rc, theta_m(k), params%theta_sat(k), params%theta_res(k),&
                                      curve_a(params,k), curve_n(params,k))
      end do
      maxit = 1_ik
      if (opts%linearize == SOIL_LIN_PICARD) maxit = max(1_ik, opts%max_picard)
      ok = .false.

      do iter = 1_ik, maxit
         call face_and_sink(params, opts, rc, n, psi_m, theta_m, psi_e, root_uptake,           &
                            kk, cc, kface, gface, sk, dsk)
         qbot = 0.0_wp
         if (opts%bottom_bc /= SOIL_BC_BEDROCK) qbot = kk(n)             ! free-drain/aquifer: K(theta_N)
         do k = 1_ik, n - 1_ik
            qface(k) = kface(k) * ((psi_m(k) - psi_m(k+1)) / params%dz_node(k) + gface(k))
         end do
         do k = 1_ik, n
            if (k >= 2_ik) then
               a(k) = -kface(k-1) / params%dz_node(k-1)
            else
               a(k) = 0.0_wp
            end if
            if (k <= n-1_ik) then
               c(k) = -kface(k) / params%dz_node(k)
            else
               c(k) = 0.0_wp
            end if
            b(k)  = cc(k) * params%dz(k) / h - a(k) - c(k) + dsk(k) * params%dz(k)
            in_k  = q_top
            if (k >= 2_ik)   in_k  = qface(k-1)
            out_k = qbot
            if (k <= n-1_ik) out_k = qface(k)
            rhs(k) = in_k - out_k - sk(k) * params%dz(k)                                        &
                     - (theta_m(k) - theta_in(k)) * params%dz(k) / h      ! Celia mass term
         end do
         call thomas_solve(a, b, c, rhs, dpsi, n)
         !----- Picard convergence on the MOISTURE increment |dtheta| [m3/m3], dimensionally      !
         !      consistent with opts%atol (a moisture tolerance) -- the old test compared the       !
         !      matric-head increment |dpsi| [m] against the same atol (mixed units).  -------------!
         err = 0.0_wp
         do k = 1_ik, n
            psi_m(k)   = psi_m(k) + dpsi(k)
            theta_prev = theta_m(k)
            theta_m(k) = soil_theta_from_psi_l(rc, params, k, psi_m(k))
            err = max(err, abs(theta_m(k) - theta_prev))
         end do
         if (maxit == 1_ik .or. err < opts%atol) then
            ok = .true.
            exit
         end if
      end do

      !----- Conservative theta update from the converged interior fluxes (telescoping). ----!
      call face_and_sink(params, opts, rc, n, psi_m, theta_m, psi_e, root_uptake,              &
                         kk, cc, kface, gface, sk, dsk)
      qbot = 0.0_wp
      if (opts%bottom_bc /= SOIL_BC_BEDROCK) qbot = kk(n)
      wface_amt = 0.0_wp
      do k = 1_ik, n - 1_ik
         qface(k)     = kface(k) * ((psi_m(k) - psi_m(k+1)) / params%dz_node(k) + gface(k))
         wface_amt(k) = qface(k) * h                             ! [m] downward depth crossing face k over h
      end do
      uptake_amt = 0.0_wp
      do k = 1_ik, n
         in_k  = q_top
         if (k >= 2_ik)   in_k  = qface(k-1)
         out_k = qbot
         if (k <= n-1_ik) out_k = qface(k)
         theta_out(k) = theta_in(k) + h / params%dz(k) * (in_k - out_k - sk(k) * params%dz(k))
         uptake_amt   = uptake_amt + sk(k) * params%dz(k) * rho_h2o * h
      end do
      drainage_amt = qbot * rho_h2o * h
   end subroutine soil_water_step_implicit

   !---------------------------------------------------------------------------------------!
   ! soil_water_time_deriv -- the EXPLICIT Richards RHS dtheta_k/dt [1/s] at the current state, for !
   ! the IMEX-ARK fast integrator (docs/dev_plans/MEDS_IMEX_ARK_DESIGN.md). Same flux-divergence + root    !
   ! sink form as soil_water_step_implicit's conservative update, but evaluated ONCE at the current !
   ! theta (no Celia/frozen BE solve): reuses face_and_sink so upstream K, the Zeng-Decker gravity  !
   ! factor, and the psi-limited sink are IDENTICAL to the split -- no re-derivation. The top flux    !
   ! q_top, equilibrium psi_e, and root_uptake are the frozen surface BCs. Commits nothing.          !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_water_time_deriv(theta, params, opts, n, q_top, psi_e, root_uptake,      &
                                         dtheta_dt, drainage_rate, uptake_rate)
      real(wp),            intent(in)  :: theta(n_soil_layer_max)
      type(soil_params_t), intent(in)  :: params
      type(soil_opts_t),   intent(in)  :: opts
      integer(ik),         intent(in)  :: n
      real(wp),            intent(in)  :: q_top, psi_e(n_soil_layer_max), root_uptake(n_soil_layer_max)
      real(wp),            intent(out) :: dtheta_dt(n_soil_layer_max)   !< [1/s]     dtheta/dt per layer (0 for k>n)
      real(wp),            intent(out) :: drainage_rate                 !< [kg/m2/s] bottom drainage
      real(wp),            intent(out) :: uptake_rate                   !< [kg/m2/s] total psi-limited root uptake

      real(wp), dimension(n_soil_layer_max) :: psi_m, kk, cc, kface, gface, qface, sk, dsk
      real(wp)    :: qbot, in_k, out_k
      integer(ik) :: k, rc

      dtheta_dt = 0.0_wp
      rc = params%retention

      do k = 1_ik, n
         psi_m(k) = soil_psi_from_theta(rc, theta(k), params%theta_sat(k), params%theta_res(k),   &
                                      curve_a(params, k), curve_n(params, k))
      end do
      call face_and_sink(params, opts, rc, n, psi_m, theta, psi_e, root_uptake,                  &
                         kk, cc, kface, gface, sk, dsk)
      qbot = 0.0_wp
      if (opts%bottom_bc /= SOIL_BC_BEDROCK) qbot = kk(n)
      do k = 1_ik, n - 1_ik
         qface(k) = kface(k) * ((psi_m(k) - psi_m(k+1)) / params%dz_node(k) + gface(k))
      end do

      uptake_rate = 0.0_wp
      do k = 1_ik, n
         in_k  = q_top
         if (k >= 2_ik)   in_k  = qface(k-1)
         out_k = qbot
         if (k <= n-1_ik) out_k = qface(k)
         dtheta_dt(k) = (in_k - out_k) / params%dz(k) - sk(k)
         uptake_rate  = uptake_rate + sk(k) * params%dz(k) * rho_h2o
      end do
      drainage_rate = qbot * rho_h2o
   end subroutine soil_water_time_deriv

   !----- Fill node K/C, upstream face K, Zeng-Decker gravity factor, and the psi-limited sink !
   !      + its psi-derivative, all at the current iterate (psi_m, theta_m).                    !
   !---------------------------------------------------------------------------------------!
   pure subroutine face_and_sink(params, opts, rc, n, psi_m, theta_m, psi_e, root_uptake,      &
                                 kk, cc, kface, gface, sk, dsk)
      type(soil_params_t), intent(in)  :: params
      type(soil_opts_t),   intent(in)  :: opts
      integer(ik),         intent(in)  :: rc, n
      real(wp),            intent(in)  :: psi_m(n_soil_layer_max), theta_m(n_soil_layer_max)
      real(wp),            intent(in)  :: psi_e(n_soil_layer_max), root_uptake(n_soil_layer_max)
      real(wp),            intent(out) :: kk(n_soil_layer_max), cc(n_soil_layer_max)
      real(wp),            intent(out) :: kface(n_soil_layer_max), gface(n_soil_layer_max)
      real(wp),            intent(out) :: sk(n_soil_layer_max), dsk(n_soil_layer_max)
      integer(ik) :: k
      real(wp)    :: dh, fwilt, dfwilt
      do k = 1_ik, n
         kk(k) = soil_hydr_cond_from_theta(rc, theta_m(k), params%theta_sat(k), params%theta_res(k),      &
                                curve_a(params,k), curve_n(params,k), params%ksat(k))
         cc(k) = soil_moist_cap_from_psi(rc, psi_m(k), params%theta_sat(k), params%theta_res(k),        &
                                curve_a(params,k), curve_n(params,k))
         call f_wilt_ramp(psi_m(k), opts%psi_wilt, opts%psi_open, fwilt, dfwilt)
         sk(k)  = root_uptake(k) / (rho_h2o * params%dz(k)) * fwilt
         dsk(k) = root_uptake(k) / (rho_h2o * params%dz(k)) * dfwilt
      end do
      do k = 1_ik, n - 1_ik
         !----- Upstream K by the total-head gradient (down if dh >= 0). ------------------!
         dh = (psi_m(k) - psi_m(k+1)) + params%dz_node(k)
         if (dh >= 0.0_wp) then
            kface(k) = kk(k)
         else
            kface(k) = kk(k+1)
         end if
         if (opts%zeng_decker) then
            gface(k) = -(psi_e(k) - psi_e(k+1)) / params%dz_node(k)      ! = 1 for a linear-midpoint psi_e
         else
            gface(k) = 1.0_wp
         end if
      end do
   end subroutine face_and_sink

   !----- Retention-integral equilibrium potential psi_E,k = psi( layer-mean theta_eq ), with   !
   !      theta_eq(z) = theta( z_wt - z ) the hydrostatic-equilibrium moisture (design 3e).      !
   !---------------------------------------------------------------------------------------!
   subroutine compute_psi_e(params, rc, z_wt, n, psi_e)
      type(soil_params_t), intent(in)  :: params
      integer(ik),         intent(in)  :: rc, n
      real(wp),            intent(in)  :: z_wt
      real(wp),            intent(out) :: psi_e(n_soil_layer_max)
      integer(ik), parameter :: nq = 5_ik
      integer(ik) :: k, j
      real(wp)    :: theta_e, zc, psi_eq
      psi_e = 0.0_wp
      do k = 1_ik, n
         theta_e = 0.0_wp
         do j = 1_ik, nq
            zc     = params%soil_layer_z(k+1) + (real(j,wp) - 0.5_wp) / real(nq,wp) * params%dz(k)
            psi_eq = z_wt - zc
            theta_e = theta_e + soil_theta_from_psi_l(rc, params, k, psi_eq) / real(nq,wp)
         end do
         psi_e(k) = soil_psi_from_theta(rc, theta_e, params%theta_sat(k), params%theta_res(k),   &
                                      curve_a(params,k), curve_n(params,k))
      end do
   end subroutine compute_psi_e

   !=======================================================================================!
   !  Small helpers                                                                         !
   !=======================================================================================!

   !----- theta(psi) with the layer-k curve params (thin wrapper). ------------------------!
   pure function soil_theta_from_psi_l(rc, params, k, psi) result(theta)
      integer(ik),         intent(in) :: rc, k
      type(soil_params_t), intent(in) :: params
      real(wp),            intent(in) :: psi
      real(wp)                        :: theta
      theta = soil_theta_from_psi(rc, psi, params%theta_sat(k), params%theta_res(k),             &
                                curve_a(params,k), curve_n(params,k))
   end function soil_theta_from_psi_l

   !----- Curve parameter accessors: alpha/n for vG, psi_sat/b for Campbell. --------------!
   pure function curve_a(params, k) result(a)
      type(soil_params_t), intent(in) :: params
      integer(ik),         intent(in) :: k
      real(wp)                        :: a
      if (params%retention == SOIL_RETENTION_CAMPBELL) then
         a = params%psi_sat(k)
      else
         a = params%vg_alpha(k)
      end if
   end function curve_a

   pure function curve_n(params, k) result(nn)
      type(soil_params_t), intent(in) :: params
      integer(ik),         intent(in) :: k
      real(wp)                        :: nn
      if (params%retention == SOIL_RETENTION_CAMPBELL) then
         nn = params%b_camp(k)
      else
         nn = params%vg_n(k)
      end if
   end function curve_n

   !----- Smooth wilting ramp f_wilt(psi) in [0,1] and its derivative. ---------------------!
   pure subroutine f_wilt_ramp(psi, psi_wilt, psi_open, f, df)
      real(wp), intent(in)  :: psi, psi_wilt, psi_open
      real(wp), intent(out) :: f, df
      real(wp) :: span
      span = psi_open - psi_wilt                          ! > 0 (psi_open the less negative)
      if (psi >= psi_open) then
         f = 1.0_wp
         df = 0.0_wp
      else if (psi <= psi_wilt) then
         f = 0.0_wp
         df = 0.0_wp
      else
         f = (psi - psi_wilt) / span
         df = 1.0_wp / span
      end if
   end subroutine f_wilt_ramp

   !----- Ground evaporation: pore-space RH (Philip) + Swenson-Lawrence DSL resistance. ----!
   pure function ground_evaporation(theta1, psi1, params, forcing, opts) result(e_soil)
      real(wp),               intent(in) :: theta1, psi1
      type(soil_params_t),    intent(in) :: params
      type(chydro_forcing_t), intent(in) :: forcing
      type(soil_opts_t),      intent(in) :: opts
      real(wp) :: e_soil, alpha_soil, q_g, theta_init, dsl, dvap, phi, phi_air, tau, r_soil
      alpha_soil = exp(max(-40.0_wp, psi1 * grav / (r_wv * forcing%t_ground)))
      q_g        = alpha_soil * sat_specific_humidity(forcing%t_ground, p_std)
      theta_init = opts%dsl_theta_init * params%theta_sat(1)
      if (theta1 < theta_init) then
         dsl = opts%dsl_dmax * (theta_init - theta1) / max(theta_init, tiny_num)
      else
         dsl = 0.0_wp
      end if
      dvap    = 2.12e-5_wp * (forcing%t_ground / 273.15_wp) ** 1.75_wp
      phi     = params%theta_sat(1)
      phi_air = max(phi - theta1, tiny_num)
      tau     = phi_air ** (10.0_wp / 3.0_wp) / max(phi * phi, tiny_num)     ! Millington-Quirk
      r_soil  = dsl / max(dvap * tau, tiny_num)
      !----- AREA-weighted tile flux: only the snow-free fraction of the ground evaporates, and it does  !
      !      so at the BARE-SOIL resistance (r_aero = 1/ggnet, r_soil from the dry surface layer). The     !
      !      alternative of dividing r_aero by (1-snowfac) throttles only the aerodynamic leg, which is    !
      !      the same thing ONLY when r_soil = 0 -- see chydro_forcing_t%snow_free_frac. snow_free_frac    !
      !      defaults to 1, so a caller that does not model snow is unaffected. -------------------------!
      e_soil  = forcing%snow_free_frac * forcing%rho_air * (q_g - forcing%q_air)                        &
                / (forcing%r_aero + r_soil)
      e_soil  = max(0.0_wp, e_soil)                        ! no dew in v1
   end function ground_evaporation

end module meds_soil_water
