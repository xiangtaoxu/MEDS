!==========================================================================================!
! meds_column_hydrology -- THE stateless seam of the 1-D soil-water column (design            !
! MEDS_COLUMN_HYDROLOGY_DESIGN.md). It advances one patch's prognostic soil moisture + ponding !
! (+ a lumped aquifer) over a fast step dt: canopy-throughfall infiltration (conductivity-       !
! limited), ground evaporation, the multi-layer soil-water balance (implicit backward-Euler       !
! Thomas), and a water-potential-limited root-uptake sink; it returns the boundary fluxes plus     !
! the per-layer matric potential psi_soil that closes the plant-hydraulics boundary condition.      !
! Canopy interception is a separate per-cohort routine the caller sweeps top->bottom.               !
!                                                                                          !
! P2 numerics/physics (this file):                                                                  !
!   * upstream-weighted interface conductivity;                                                      !
!   * retention-integral Zeng-Decker equilibrium correction on interior faces (opts%zeng_decker);     !
!   * Celia (1990) modified-Picard OR frozen-coefficient linearization (opts%linearize);               !
!   * adaptive step-doubling substep control OR a fixed count (opts%substep);                           !
!   * free-drainage / bedrock / SIMTOP-aquifer bottom BC (+ diagnosed water table z_wt);                 !
!   * Dunne saturation-excess (f_sat) runoff.                                                             !
! State-free and device-eligibility is deferred to P5 (the inner step takes derived types for now).       !
! Deferred within P2: the Neumann->Dirichlet ponded-surface top-BC switch (noted in 3d).                    !
!==========================================================================================!
module meds_column_hydrology
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, grav, r_wv, tiny_num, grav_head, p_std
   use meds_biophysics_types, only : soil_column_t, chydro_forcing_t, soil_params_t,          &
                                     soil_opts_t, chydro_flux_t, n_soil_layer_max,             &
                                     SOIL_BC_BEDROCK, SOIL_BC_AQUIFER, SOIL_RETENTION_CAMPBELL, &
                                     SOIL_LIN_PICARD, SOIL_SUBSTEP_FIXED
   use meds_soil_parameters,  only : soil_psi_of_theta, soil_theta_of_psi, soil_hydr_cond,     &
                                     soil_moist_cap
   use meds_numerics,         only : thomas_solve
   use meds_thermo,           only : sat_specific_humidity
   implicit none
   private

   public :: column_hydrology_flux, soil_water_tendency, intercept_canopy_layer

   !----- Unconfined-aquifer specific yield (fraction) for the z_wt <-> w_aquifer diagnosis. -!
   real(wp), parameter :: SPECIFIC_YIELD = 0.2_wp

contains

   !---------------------------------------------------------------------------------------!
   ! Per-cohort canopy interception (design 3c). One canopy layer: capacity-limited bucket    !
   ! with a Beer interception fraction; the caller sweeps the height-sorted cohorts top->bottom !
   ! feeding each `throughfall` as the next cohort's `rain_above`. `leaf_water` is the per-cohort !
   ! prognostic film [kg/m2 ground]; `sigma_w` (wetted fraction) is exported for the future        !
   ! canopy-air-space evaporation. Not `elemental`: the cascade is a sequential recurrence.         !
   !---------------------------------------------------------------------------------------!
   pure subroutine intercept_canopy_layer(leaf_water, rain_above, lai, sai, e_canopy, dt,     &
                                          dewmx, k_int, alpha_pi, throughfall, drip, sigma_w)
      real(wp), intent(inout) :: leaf_water
      real(wp), intent(in)    :: rain_above, lai, sai, e_canopy, dt, dewmx, k_int, alpha_pi
      real(wp), intent(out)   :: throughfall, drip, sigma_w
      real(wp) :: pai, f_pi, w_max, q_grab, room, q_intr
      pai    = lai + sai
      f_pi   = alpha_pi * (1.0_wp - exp(-k_int * pai))
      w_max  = dewmx * pai
      q_grab = f_pi * rain_above
      room   = max(0.0_wp, w_max - leaf_water) / dt
      q_intr = min(q_grab, room + e_canopy)                      ! bounded by capacity + evap headroom
      leaf_water  = min(max(leaf_water + (q_intr - e_canopy) * dt, 0.0_wp), w_max)
      drip        = max(0.0_wp, q_grab - q_intr)
      throughfall = (rain_above - q_grab) + drip                 ! gap-throughfall + drip -> below
      if (w_max > tiny_num) then
         sigma_w = min(1.0_wp, (leaf_water / w_max) ** (2.0_wp / 3.0_wp))
      else
         sigma_w = 0.0_wp
      end if
   end subroutine intercept_canopy_layer

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
      real(wp), dimension(n_soil_layer_max) :: theta0, theta1, psi_e
      real(wp) :: z_wt, z_bottom, psi1, kn1, e_soil, q_inf_max, q_avail, infl, q_top
      real(wp) :: f_sat, q_over, q_liq, drain_amt, uptake_amt, clip_ex, deficit, want, give
      real(wp) :: q_drai, recharge_amt, baseflow_amt, site_drain, wsurf, runoff, w0, w1
      real(wp) :: w_surf0, w_aqf0
      logical  :: ok

      n  = params%n_active
      rc = params%retention
      z_bottom = params%soil_layer_z(n+1)                                 ! deepest interface (<= 0)
      theta0(1:n) = col%theta(1:n)
      w_surf0 = col%w_surface                                             ! initial stores (for the budget)
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
      psi1 = soil_psi_of_theta(rc, theta0(1), params%theta_sat(1), params%theta_res(1),        &
                               curve_a(params,1), curve_n(params,1))
      kn1  = soil_hydr_cond(rc, theta0(1), params%theta_sat(1), params%theta_res(1),            &
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
      clip_ex = 0.0_wp ; deficit = 0.0_wp
      do k = 1_ik, n
         if (theta1(k) > params%theta_sat(k)) then
            clip_ex   = clip_ex + (theta1(k) - params%theta_sat(k)) * params%dz(k) * rho_h2o
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
         theta1(k) = max(theta1(k), params%theta_res(k))                 ! hard residual floor (last resort)
      end do

      !----- Commit soil moisture; export psi_soil [MPa]. ----------------------------------!
      col%theta(1:n) = theta1(1:n)
      do k = 1_ik, n
         flux%psi_soil(k) = grav_head * soil_psi_of_theta(rc, theta1(k), params%theta_sat(k),  &
                            params%theta_res(k), curve_a(params,k), curve_n(params,k))
      end do

      !----- Ponding store (+ theta-clip overflow) and surface runoff. ----------------------!
      wsurf  = w_surf0 + (q_liq - infl) * dt + clip_ex
      runoff = q_over + max(0.0_wp, wsurf - opts%w_pond_max) / dt
      wsurf  = min(wsurf, opts%w_pond_max)
      col%w_surface = wsurf

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
      flux%mass_resid     = (w1 - w0) - dt * (forcing%precip_ground - e_soil - site_drain       &
                            - flux%uptake_total - runoff)
      flux%nsub      = nsub
      flux%converged = ok
      if (opts%debug_error .and. abs(flux%mass_resid) > opts%atol) then
         error stop 'column_hydrology_flux: mass budget did not close'
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

      drainage_tot = 0.0_wp ; uptake_tot = 0.0_wp ; nsub = 0_ik ; ok = .true.
      wface_tot = 0.0_wp

      if (opts%substep == SOIL_SUBSTEP_FIXED) then
         !----- Warp-uniform fixed count (GPU-friendly). ----------------------------------!
         nfix = max(1_ik, min(opts%max_substep, nint(dt / max(opts%h_init, tiny_num), ik)))
         h    = dt / real(nfix, wp)
         do k = 1_ik, nfix
            call soil_be_single_step(theta, params, opts, rc, n, h, q_top, psi_e, root_uptake, &
                                     th_big, dr_b, up_b, wf_b, dum)
            ok = ok .and. dum                          ! aggregate inner convergence (was discarded)
            theta(1:n) = th_big(1:n)
            drainage_tot = drainage_tot + dr_b ; uptake_tot = uptake_tot + up_b
            wface_tot(1:n) = wface_tot(1:n) + wf_b(1:n)
         end do
         nsub = nfix
         wflux_out(1:n) = wface_tot(1:n) / dt
         return
      end if

      !----- Adaptive step-doubling (BE is L-stable; substep only for accuracy). ------------!
      t = 0.0_wp ; h = min(opts%h_init, dt) ; hmin = dt / real(opts%max_substep, wp)
      do while (t < dt - 1.0e-9_wp * dt .and. nsub < opts%max_substep)
         h = min(h, dt - t)
         call soil_be_single_step(theta, params, opts, rc, n, h,        q_top, psi_e,          &
                                  root_uptake, th_big, dr_b, up_b, wf_b, dum)
         call soil_be_single_step(theta, params, opts, rc, n, 0.5_wp*h, q_top, psi_e,          &
                                  root_uptake, th_h1, dr_1, up_1, wf_1, dum)
         call soil_be_single_step(th_h1, params, opts, rc, n, 0.5_wp*h, q_top, psi_e,          &
                                  root_uptake, th_two, dr_2, up_2, wf_2, dum)
         err = 1.0e-12_wp
         do k = 1_ik, n
            err = max(err, abs(th_two(k) - th_big(k)) / (opts%atol + opts%rtol * abs(th_two(k))))
         end do
         if (err <= 1.0_wp .or. h <= hmin * 1.0001_wp) then
            theta(1:n)   = th_two(1:n)                         ! the more accurate (two half-steps)
            drainage_tot = drainage_tot + dr_1 + dr_2
            uptake_tot   = uptake_tot + up_1 + up_2
            wface_tot(1:n) = wface_tot(1:n) + wf_1(1:n) + wf_2(1:n)
            t    = t + h ; nsub = nsub + 1_ik
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
   subroutine soil_be_single_step(theta_in, params, opts, rc, n, h, q_top, psi_e, root_uptake, &
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
         psi_m(k) = soil_psi_of_theta(rc, theta_m(k), params%theta_sat(k), params%theta_res(k),&
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
            if (k >= 2_ik) then ; a(k) = -kface(k-1) / params%dz_node(k-1) ; else ; a(k) = 0.0_wp ; end if
            if (k <= n-1_ik) then ; c(k) = -kface(k) / params%dz_node(k) ; else ; c(k) = 0.0_wp ; end if
            b(k)  = cc(k) * params%dz(k) / h - a(k) - c(k) + dsk(k) * params%dz(k)
            in_k  = q_top ; if (k >= 2_ik)   in_k  = qface(k-1)
            out_k = qbot  ; if (k <= n-1_ik) out_k = qface(k)
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
            theta_m(k) = soil_theta_of_psi_l(rc, params, k, psi_m(k))
            err = max(err, abs(theta_m(k) - theta_prev))
         end do
         if (maxit == 1_ik .or. err < opts%atol) then ; ok = .true. ; exit ; end if
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
         in_k  = q_top ; if (k >= 2_ik)   in_k  = qface(k-1)
         out_k = qbot  ; if (k <= n-1_ik) out_k = qface(k)
         theta_out(k) = theta_in(k) + h / params%dz(k) * (in_k - out_k - sk(k) * params%dz(k))
         uptake_amt   = uptake_amt + sk(k) * params%dz(k) * rho_h2o * h
      end do
      drainage_amt = qbot * rho_h2o * h
   end subroutine soil_be_single_step

   !---------------------------------------------------------------------------------------!
   ! soil_water_tendency -- the EXPLICIT Richards RHS dtheta_k/dt [1/s] at the current state, for  !
   ! the IMEX-ARK fast integrator (archive/MEDS_IMEX_ARK_DESIGN.md). Same flux-divergence + root    !
   ! sink form as soil_be_single_step's conservative update (:362-379), but evaluated ONCE at the   !
   ! current theta (no Celia/frozen BE solve): reuses face_and_sink so upstream K, the Zeng-Decker  !
   ! gravity factor, and the psi-limited sink are IDENTICAL to the split -- no re-derivation. The    !
   ! top flux q_top, equilibrium psi_e, and root_uptake are the frozen surface BCs (explicit part    !
   ! of the ARK split), passed in as in soil_be_single_step. Commits nothing.                        !
   !---------------------------------------------------------------------------------------!
   pure subroutine soil_water_tendency(theta, params, opts, n, q_top, psi_e, root_uptake,        &
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
         psi_m(k) = soil_psi_of_theta(rc, theta(k), params%theta_sat(k), params%theta_res(k),   &
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
         in_k  = q_top ; if (k >= 2_ik)   in_k  = qface(k-1)
         out_k = qbot  ; if (k <= n-1_ik) out_k = qface(k)
         dtheta_dt(k) = (in_k - out_k) / params%dz(k) - sk(k)
         uptake_rate  = uptake_rate + sk(k) * params%dz(k) * rho_h2o
      end do
      drainage_rate = qbot * rho_h2o
   end subroutine soil_water_tendency

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
         kk(k) = soil_hydr_cond(rc, theta_m(k), params%theta_sat(k), params%theta_res(k),      &
                                curve_a(params,k), curve_n(params,k), params%ksat(k))
         cc(k) = soil_moist_cap(rc, psi_m(k), params%theta_sat(k), params%theta_res(k),        &
                                curve_a(params,k), curve_n(params,k))
         call f_wilt_ramp(psi_m(k), opts%psi_wilt, opts%psi_open, fwilt, dfwilt)
         sk(k)  = root_uptake(k) / (rho_h2o * params%dz(k)) * fwilt
         dsk(k) = root_uptake(k) / (rho_h2o * params%dz(k)) * dfwilt
      end do
      do k = 1_ik, n - 1_ik
         !----- Upstream K by the total-head gradient (down if dh >= 0). ------------------!
         dh = (psi_m(k) - psi_m(k+1)) + params%dz_node(k)
         if (dh >= 0.0_wp) then ; kface(k) = kk(k) ; else ; kface(k) = kk(k+1) ; end if
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
            theta_e = theta_e + soil_theta_of_psi_l(rc, params, k, psi_eq) / real(nq,wp)
         end do
         psi_e(k) = soil_psi_of_theta(rc, theta_e, params%theta_sat(k), params%theta_res(k),   &
                                      curve_a(params,k), curve_n(params,k))
      end do
   end subroutine compute_psi_e

   !=======================================================================================!
   !  Small helpers                                                                         !
   !=======================================================================================!

   !----- theta(psi) with the layer-k curve params (thin wrapper). ------------------------!
   pure function soil_theta_of_psi_l(rc, params, k, psi) result(theta)
      integer(ik),         intent(in) :: rc, k
      type(soil_params_t), intent(in) :: params
      real(wp),            intent(in) :: psi
      real(wp)                        :: theta
      theta = soil_theta_of_psi(rc, psi, params%theta_sat(k), params%theta_res(k),             &
                                curve_a(params,k), curve_n(params,k))
   end function soil_theta_of_psi_l

   !----- Curve parameter accessors: alpha/n for vG, psi_sat/b for Campbell. --------------!
   pure function curve_a(params, k) result(a)
      type(soil_params_t), intent(in) :: params
      integer(ik),         intent(in) :: k
      real(wp)                        :: a
      if (params%retention == SOIL_RETENTION_CAMPBELL) then ; a = params%psi_sat(k)
      else ; a = params%vg_alpha(k) ; end if
   end function curve_a

   pure function curve_n(params, k) result(nn)
      type(soil_params_t), intent(in) :: params
      integer(ik),         intent(in) :: k
      real(wp)                        :: nn
      if (params%retention == SOIL_RETENTION_CAMPBELL) then ; nn = params%b_camp(k)
      else ; nn = params%vg_n(k) ; end if
   end function curve_n

   !----- Smooth wilting ramp f_wilt(psi) in [0,1] and its derivative. ---------------------!
   pure subroutine f_wilt_ramp(psi, psi_wilt, psi_open, f, df)
      real(wp), intent(in)  :: psi, psi_wilt, psi_open
      real(wp), intent(out) :: f, df
      real(wp) :: span
      span = psi_open - psi_wilt                          ! > 0 (psi_open the less negative)
      if (psi >= psi_open) then
         f = 1.0_wp ; df = 0.0_wp
      else if (psi <= psi_wilt) then
         f = 0.0_wp ; df = 0.0_wp
      else
         f = (psi - psi_wilt) / span ; df = 1.0_wp / span
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
      e_soil  = forcing%rho_air * (q_g - forcing%q_air) / (forcing%r_aero + r_soil)
      e_soil  = max(0.0_wp, e_soil)                        ! no dew in v1
   end function ground_evaporation

end module meds_column_hydrology
