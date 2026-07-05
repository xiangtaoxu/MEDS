!==========================================================================================!
! meds_vertical_hydrology -- THE stateless seam of the 1-D soil-water column (design            !
! MEDS_VERTICAL_HYDROLOGY_DESIGN.md). It advances one patch's prognostic soil moisture + ponding !
! over a fast step dt: canopy-throughfall infiltration (conductivity-limited), ground            !
! evaporation, the multi-layer soil-water balance (implicit backward-Euler Thomas, plain-gravity !
! flux, free-drainage / bedrock bottom BC), and a water-potential-limited root-uptake sink; it    !
! returns the boundary fluxes plus the per-layer matric potential psi_soil that closes the plant- !
! hydraulics boundary condition. Canopy interception is a separate per-cohort routine the caller  !
! sweeps top->bottom.                                                                             !
!                                                                                          !
! State-free and device-eligible: it takes plain value types (never a site_t). This is the P0/P1  !
! MVP -- frozen-coefficient single implicit solve per call; adaptive step-doubling, Celia Picard,  !
! Zeng-Decker, the aquifer BC and van-Genuchten-only extras are P2 (see the design).                !
!==========================================================================================!
module meds_vertical_hydrology
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, grav, r_wv, tiny_num, grav_head, p_std
   use meds_biophysics_types, only : soil_column_t, vhydro_forcing_t, soil_params_t,          &
                                     soil_opts_t, vhydro_flux_t, n_soil_layer_max,             &
                                     SOIL_BC_BEDROCK, SOIL_RETENTION_CAMPBELL
   use meds_soil_parameters,  only : soil_psi_of_theta, soil_hydr_cond, soil_moist_cap
   use meds_soil_solver,      only : thomas_solve
   implicit none
   private

   public :: vertical_hydrology_flux, intercept_canopy_layer

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
   ! The seam: advance one patch soil column over dt. `col%theta(:)` + stores are the only    !
   ! mutated data; `flux` carries the boundary fluxes, exported psi_soil, and the mass-budget   !
   ! residual (asserts to ~round-off).                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine vertical_hydrology_flux(col, forcing, params, opts, dt, flux)
      type(soil_column_t),    intent(inout) :: col
      type(vhydro_forcing_t), intent(in)    :: forcing
      type(soil_params_t),    intent(in)    :: params
      type(soil_opts_t),      intent(in)    :: opts
      real(wp),               intent(in)    :: dt
      type(vhydro_flux_t),    intent(out)   :: flux

      integer(ik) :: n, k, rc
      real(wp), dimension(n_soil_layer_max) :: theta0, psi_n, kn, cn, sk, dsk, fw, dfw
      real(wp), dimension(n_soil_layer_max) :: a, b, c, rhs, dpsi, psi1, kface, qint, sl, theta1
      real(wp) :: q_top, q_bot, infl, e_soil, q_avail, q_inf_max
      real(wp) :: w0, w1, uptake_lin, uptake_real, deficit, clip_ex, runoff, wsurf, fwilt, dfwilt
      real(wp) :: give, want, in_k, out_k

      n  = params%n_active
      rc = params%retention

      !----- Snapshot + node-wise constitutive quantities at theta^n. ----------------------!
      theta0(1:n) = col%theta(1:n)
      do k = 1_ik, n
         psi_n(k) = soil_psi_of_theta(rc, theta0(k), params%theta_sat(k), params%theta_res(k),&
                                      curve_a(params,k), curve_n(params,k))
         kn(k)    = soil_hydr_cond(rc, theta0(k), params%theta_sat(k), params%theta_res(k),   &
                                   curve_a(params,k), curve_n(params,k), params%ksat(k))
         cn(k)    = soil_moist_cap(rc, psi_n(k), params%theta_sat(k), params%theta_res(k),    &
                                   curve_a(params,k), curve_n(params,k))
      end do

      !----- Interior face conductivities (arithmetic mean; upstream weighting is a P2       !
      !      refinement). kface(k) sits at the k<->k+1 interface, k = 1..n-1.                 !
      do k = 1_ik, n - 1_ik
         kface(k) = 0.5_wp * (kn(k) + kn(k+1))
      end do

      !----- Ground evaporation (Philip alpha_soil + Swenson-Lawrence DSL resistance). ------!
      e_soil = ground_evaporation(theta0(1), psi_n(1), params, forcing, opts)
      e_soil = min(e_soil, max(0.0_wp, (theta0(1) - params%theta_res(1)) * params%dz(1)       &
                               * rho_h2o / dt))

      !----- Conductivity-limited infiltration + ponding partition (design 3d). -------------!
      q_inf_max = kn(1) * (1.0_wp + (-psi_n(1)) / (-params%z_node(1)))       ! Darcy velocity [m/s]
      q_avail   = col%w_surface / dt + forcing%precip_ground                 ! [kg/m2/s]
      infl      = min(q_avail, q_inf_max * rho_h2o)                          ! [kg/m2/s]

      !----- Top face (Neumann): infiltration minus ground evaporation, as a Darcy velocity. -!
      q_top = (infl - e_soil) / rho_h2o                                      ! [m/s, positive down]

      !----- Bottom face BC. ---------------------------------------------------------------!
      if (opts%bottom_bc == SOIL_BC_BEDROCK) then
         q_bot = 0.0_wp
      else
         q_bot = kn(n)                                                       ! free drainage: unit gradient
      end if

      !----- Root-uptake sink S_k [m3/m3/s] + its psi-derivative (design 3f). ---------------!
      do k = 1_ik, n
         call f_wilt_ramp(psi_n(k), opts%psi_wilt, opts%psi_open, fwilt, dfwilt)
         fw(k)  = fwilt ; dfw(k) = dfwilt
         sk(k)  = forcing%root_uptake(k) / (rho_h2o * params%dz(k)) * fw(k)
         dsk(k) = forcing%root_uptake(k) / (rho_h2o * params%dz(k)) * dfw(k)
      end do

      !----- Interior face fluxes at theta^n (plain-gravity form q = K(dpsi/dz + 1)). --------!
      do k = 1_ik, n - 1_ik
         qint(k) = kface(k) * ((psi_n(k) - psi_n(k+1)) / params%dz_node(k) + 1.0_wp)
      end do

      !----- Assemble the tridiagonal for the increment dpsi (design 5.2). ------------------!
      do k = 1_ik, n
         if (k >= 2_ik) then
            a(k) = -kface(k-1) / params%dz_node(k-1)
         else
            a(k) = 0.0_wp
         end if
         if (k <= n - 1_ik) then
            c(k) = -kface(k) / params%dz_node(k)
         else
            c(k) = 0.0_wp
         end if
         b(k) = cn(k) * params%dz(k) / dt - a(k) - c(k) + dsk(k) * params%dz(k)
         in_k  = q_top                                     ! default for k = 1
         if (k >= 2_ik)     in_k  = qint(k-1)
         out_k = q_bot                                     ! default for k = n
         if (k <= n - 1_ik) out_k = qint(k)
         rhs(k) = in_k - out_k - sk(k) * params%dz(k)
      end do

      call thomas_solve(a, b, c, rhs, dpsi, n)
      psi1(1:n) = psi_n(1:n) + dpsi(1:n)

      !----- Conservative theta update from fluxes at psi^{n+1} (frozen K). Interior faces    !
      !      are computed once and telescope, guaranteeing mass conservation.                 !
      do k = 1_ik, n - 1_ik
         qint(k) = kface(k) * ((psi1(k) - psi1(k+1)) / params%dz_node(k) + 1.0_wp)
      end do
      uptake_lin = 0.0_wp
      do k = 1_ik, n
         sl(k) = sk(k) + dsk(k) * dpsi(k)                  ! linearized sink
         in_k  = q_top ; if (k >= 2_ik)     in_k  = qint(k-1)
         out_k = q_bot ; if (k <= n - 1_ik) out_k = qint(k)
         theta1(k)  = theta0(k) + dt / params%dz(k) * (in_k - out_k - sl(k) * params%dz(k))
         uptake_lin = uptake_lin + sl(k) * params%dz(k) * rho_h2o
      end do

      !----- Post-solve clip + theta_wp cap (design 5.1 step 5), all bookkept. --------------!
      clip_ex = 0.0_wp ; deficit = 0.0_wp
      do k = 1_ik, n
         if (theta1(k) > params%theta_sat(k)) then
            clip_ex   = clip_ex + (theta1(k) - params%theta_sat(k)) * params%dz(k) * rho_h2o / dt
            theta1(k) = params%theta_sat(k)
         end if
         if (theta1(k) < params%theta_wp(k) .and. sl(k) > 0.0_wp) then
            want      = (params%theta_wp(k) - theta1(k)) * params%dz(k) * rho_h2o / dt
            give      = min(want, sl(k) * params%dz(k) * rho_h2o)       ! give back <= extracted
            theta1(k) = theta1(k) + give * dt / (params%dz(k) * rho_h2o)
            deficit   = deficit + give
         end if
         theta1(k) = max(theta1(k), params%theta_res(k))               ! hard residual floor
      end do
      uptake_real = max(0.0_wp, uptake_lin - deficit)

      !----- Commit soil moisture; export psi_soil [MPa]. ----------------------------------!
      w0 = col%w_surface
      col%theta(1:n) = theta1(1:n)
      do k = 1_ik, n
         flux%psi_soil(k) = grav_head * soil_psi_of_theta(rc, theta1(k), params%theta_sat(k), &
                            params%theta_res(k), curve_a(params,k), curve_n(params,k))
      end do

      !----- Ponding store + surface runoff. -----------------------------------------------!
      wsurf  = w0 + (forcing%precip_ground - infl) * dt + clip_ex * dt
      runoff = max(0.0_wp, wsurf - opts%w_pond_max) / dt
      wsurf  = min(wsurf, opts%w_pond_max)
      col%w_surface = wsurf

      !----- Mass-budget closure check (should be ~round-off). ------------------------------!
      w1 = col%w_surface
      do k = 1_ik, n
         w0 = w0 + theta0(k) * params%dz(k) * rho_h2o
         w1 = w1 + theta1(k) * params%dz(k) * rho_h2o
      end do
      flux%infiltration   = infl
      flux%drainage       = q_bot * rho_h2o
      flux%runoff_surf    = runoff
      flux%soil_evap      = e_soil
      flux%uptake_total   = uptake_real
      flux%uptake_deficit = deficit
      flux%clip_excess    = clip_ex
      flux%mass_resid     = (w1 - w0) - dt * (forcing%precip_ground - e_soil - flux%drainage  &
                            - uptake_real - runoff)
      flux%nsub      = 1_ik
      flux%converged = .true.
      if (opts%debug_error .and. abs(flux%mass_resid) > opts%atol) then
         error stop 'vertical_hydrology_flux: mass budget did not close'
      end if
   end subroutine vertical_hydrology_flux

   !=======================================================================================!
   !  Private helpers                                                                       !
   !=======================================================================================!

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
      type(vhydro_forcing_t), intent(in) :: forcing
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

   !----- Saturation specific humidity [kg/kg] (Bolton 1980). ------------------------------!
   elemental function sat_specific_humidity(t_k, p_pa) result(qs)
      real(wp), intent(in) :: t_k, p_pa
      real(wp)             :: qs, esat, tc
      tc   = t_k - 273.15_wp
      esat = 611.2_wp * exp(17.67_wp * tc / (tc + 243.5_wp))
      qs   = 0.622_wp * esat / max(p_pa - 0.378_wp * esat, tiny_num)
   end function sat_specific_humidity

end module meds_vertical_hydrology
