!==========================================================================================!
! meds_ark_stepper -- time-integrators for the fast-loop column state over the pure RHS            !
! meds_column_derivs (design archive/MEDS_IMEX_ARK_DESIGN.md). The overhaul replaces the split +    !
! Picard fast step with ONE additive Runge-Kutta advance; this module is where the tableaux live.   !
!                                                                                          !
! PHASE P1 (here): the fully-EXPLICIT classical RK4 reference integrator `rk4_column_step`. It is    !
! the plan's verification ORACLE (roadmap INTEG_RK4): it shares no code with the split's backward-   !
! Euler machinery, so agreement between it and any implicit integrator rules out a shared-bug false  !
! pass. It is NOT a production integrator -- the fast loop is stiff (stiffness ratio ~6.5e5), so RK4  !
! is stable only for dt below ~2.785*tau_fast (the ~17 s plant-hydraulic mode among the INTEGRATED   !
! reservoirs; leaf energy is diagnostic, not integrated). Used as an oracle at small dt.             !
!                                                                                          !
! DEFERRED here: the L-stable IMEX-ARK stage solve (P2 arrowhead + Thomas + 2x2) and the adaptive    !
! embedded-error controller (P3). This module is their home; `rk4_column_step` establishes the       !
! state-arithmetic + stage-evaluation scaffold they reuse.                                           !
!==========================================================================================!
module meds_ark_stepper
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, tiny_num
   use meds_numerics,         only : adaptive_step_update
   use meds_thermo,           only : uext_to_temp, cas_temp_of_enthalpy, sat_specific_humidity, &
                                     internal_energy_liquid, cas_enthalpy_of_temp
   use meds_biophysics_types, only : n_soil_layer_max, soil_energy_column_t, energy_forcing_t, energy_flux_t
   use meds_plant_types,      only : N_HYDRO, NODE_LEAF, NODE_WOOD, hydro_env_t, hydro_flux_t
   use meds_column_energy,    only : soil_energy_flux
   use meds_plant_hydraulics, only : solve_plant_water
   use meds_column_derivs,    only : column_state_t, column_frozen_t, column_tend_t, column_derivs, &
                                     surface_state_t, surface_frozen_t, surface_tend_t, surface_derivs, &
                                     stage_bflux_t, column_bflux_t
   implicit none
   private

   public :: rk4_column_step, imex_euler_column_step, adaptive_imex_march
   public :: ark2_column_step, adaptive_ark_march, bflux_zero, bflux_add

   !----- default step-doubling error-scale tolerances (WRMS: |dy| <= atol + rtol*|y|). ----------!
   real(wp), parameter :: ATOL_ENTH = 5.0e1_wp    !< [J/kg]   (~0.05 K in enthalpy)
   real(wp), parameter :: ATOL_SHV  = 1.0e-6_wp   !< [kg/kg]
   real(wp), parameter :: ATOL_CO2  = 1.0e-1_wp   !< [umol/mol]
   real(wp), parameter :: ATOL_SE   = 1.0e3_wp    !< [J/m3]
   real(wp), parameter :: ATOL_TH   = 1.0e-5_wp   !< [m3/m3]
   real(wp), parameter :: ATOL_PSI  = 1.0e-3_wp   !< [MPa]

contains

   !---------------------------------------------------------------------------------------!
   ! adaptive_imex_march -- integrate the column from t=0 to t_end with a STEP-DOUBLING adaptive     !
   ! controller over the coupled IMEX-Euler step (the plan's P3 adaptive machinery, using the shipped !
   ! adaptive_step_update primitive; p=1 -> exponent -1/(p+1) = -1/2). Each trial compares one step   !
   ! of dt against two of dt/2; the local extrapolation (the two-half-step result) is committed and   !
   ! the WRMS of their difference drives accept/reject + the next dt. Reports the step + reject count. !
   ! The embedded-estimate ARK tableau (2nd order, one solve/stage) is the follow-on that replaces    !
   ! step-doubling's 3x cost; this establishes the adaptive-controller contract.                      !
   !---------------------------------------------------------------------------------------!
   subroutine adaptive_imex_march(y0, fro, n, nsl, t_end, rtol, dt_init, y_out, nsteps, nrej)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: t_end, rtol, dt_init
      type(column_state_t),  intent(out) :: y_out
      integer(ik),           intent(out) :: nsteps, nrej

      type(column_state_t) :: y, y_big, y_h, y_small
      real(wp)             :: t, dt, err, fac
      real(wp), parameter  :: DT_FLOOR = 1.0e-2_wp, SAFETY = 0.9_wp, FMIN = 0.2_wp, FMAX = 5.0_wp
      integer(ik), parameter :: NP = 8_ik
      real(wp),    parameter :: RLX = 0.6_wp

      call state_init(y0, n, nsl, y)
      t = 0.0_wp ; dt = min(dt_init, t_end) ; nsteps = 0_ik ; nrej = 0_ik

      do
         if (t >= t_end - tiny_num) exit
         dt = min(dt, t_end - t)
         call imex_euler_column_step(y,   fro, n, nsl, dt,          y_big,   niter=NP, relax=RLX)
         call imex_euler_column_step(y,   fro, n, nsl, 0.5_wp*dt,   y_h,     niter=NP, relax=RLX)
         call imex_euler_column_step(y_h, fro, n, nsl, 0.5_wp*dt,   y_small, niter=NP, relax=RLX)
         err = state_wrms(y_big, y_small, y, n, nsl, rtol)
         fac = adaptive_step_update(max(err, tiny_num), SAFETY, FMIN, FMAX)
         if (err <= 1.0_wp .or. dt <= DT_FLOOR) then
            call state_init(y_small, n, nsl, y)       ! local extrapolation: commit the finer result
            t = t + dt ; nsteps = nsteps + 1_ik
            dt = dt * fac
         else
            nrej = nrej + 1_ik
            dt = dt * fac                        ! reject + shrink (fac < 1 since err > 1)
         end if
      end do
      call state_init(y, n, nsl, y_out)
   end subroutine adaptive_imex_march

   !----- WRMS error scale of (a - b), normalized per reservoir by atol + rtol*|y_ref|. -----------!
   pure function state_wrms(a, b, y_ref, n, nsl, rtol) result(err)
      type(column_state_t), intent(in) :: a, b, y_ref
      integer(ik),          intent(in) :: n, nsl
      real(wp),             intent(in) :: rtol
      real(wp)    :: err, s
      integer(ik) :: k, i, cnt
      s = 0.0_wp ; cnt = 0_ik
      s = s + ((a%cas_enthalpy - b%cas_enthalpy) / (ATOL_ENTH + rtol*abs(y_ref%cas_enthalpy)))**2 ; cnt = cnt + 1_ik
      s = s + ((a%cas_shv - b%cas_shv) / (ATOL_SHV + rtol*abs(y_ref%cas_shv)))**2 ; cnt = cnt + 1_ik
      s = s + ((a%cas_co2 - b%cas_co2) / (ATOL_CO2 + rtol*abs(y_ref%cas_co2)))**2 ; cnt = cnt + 1_ik
      do k = 1_ik, nsl
         s = s + ((a%soil_energy(k) - b%soil_energy(k)) / (ATOL_SE + rtol*abs(y_ref%soil_energy(k))))**2
         cnt = cnt + 1_ik    ! theta is operator-split (frozen across ESDIRK stages) -> its diff is identically
      end do                 ! zero; excluding it keeps the WRMS from being diluted by nsl no-op terms

      do i = 1_ik, n
         s = s + ((a%psi(NODE_LEAF,i) - b%psi(NODE_LEAF,i)) / (ATOL_PSI + rtol*abs(y_ref%psi(NODE_LEAF,i))))**2
         s = s + ((a%psi(NODE_WOOD,i) - b%psi(NODE_WOOD,i)) / (ATOL_PSI + rtol*abs(y_ref%psi(NODE_WOOD,i))))**2
         cnt = cnt + 2_ik
      end do
      err = sqrt(s / real(cnt, wp))
   end function state_wrms

   !---------------------------------------------------------------------------------------!
   ! column_be_stage -- ONE backward-Euler stage of the STIFF, backward-Euler-integrable block only:  !
   ! the CAS twins (BE-in-atm; np>1 -> the coupled 2x2 Newton arrowhead) + soil heat + soil water     !
   ! (BE-Thomas), driven by the frozen surface sources. Plant hydraulics is DELIBERATELY EXCLUDED --  !
   ! solve_plant_water is an EXACT matrix exponential, not a backward-Euler stage, so it cannot ride   !
   ! an ESDIRK accumulation (it would drop the order + overshoot psi); psi is PASSED THROUGH here and  !
   ! advanced separately by advance_hydraulics_full over the full step. This is the reusable ESDIRK    !
   ! stage primitive: for the CAS+soil block it solves Y = base + dt*f_I(Y), so both imex_euler_column_!
   ! step (gamma=1) and each ark2 stage (gamma*dt) are just a column_be_stage call. Reuses the         !
   ! validated production kernels -- no new numerics.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine column_be_stage(y, fro, n, nsl, dt, y_out, niter, relax, bf)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      integer(ik), optional, intent(in)  :: niter    !< 1 = uncoupled BE baseline; >1 = coupled leaf<->CAS Newton
      real(wp),    optional, intent(in)  :: relax    !< vestigial (kept for API symmetry; ignored)
      type(stage_bflux_t), optional, intent(out) :: bf  !< per-stage boundary-flux RATES for the ARK ledger

      type(surface_state_t)      :: ys
      type(surface_frozen_t)     :: fs
      type(surface_tend_t)       :: sf
      type(soil_energy_column_t) :: se
      type(energy_forcing_t)     :: eforc
      type(energy_flux_t)        :: eflux
      real(wp)    :: t_ground, fliq1, wmass1, wcap, ccap, gah, gaw, gac
      real(wp)    :: enth1, shv1, e_infil, e_runof, e_drain, t_cas1
      integer(ik) :: k, np, nfeval
      logical     :: ok

      np = 1_ik ; if (present(niter)) np = max(1_ik, niter)

      call state_init(y, n, nsl, y_out)

      !----- diagnose the soil-top temperature so the ground skin sees the current state. --------!
      wmass1 = y%theta(1) * rho_h2o
      call uext_to_temp(y%soil_energy(1), wmass1, fro%therm%soil_dry_heat_capacity(1), t_ground, fliq1)

      wcap = fro%surf%wcap ; ccap = fro%surf%ccap
      gah  = fro%surf%gah  ; gaw  = fro%surf%gaw ; gac = fro%surf%gac
      fs = fro%surf ; fs%t_ground = t_ground

      !----- CAS enthalpy + humidity. np==1: the uncoupled single-BE-pass baseline. np>1: a DIRECT 2x2  !
      !      Newton solve of the coupled backward-Euler surface block (the arrowhead). The FINAL sf     !
      !      drives the soil sinks (single-flux-per-interface).                                         !
      if (np <= 1_ik) then
         ys%cas_enthalpy = y%cas_enthalpy ; ys%cas_shv = y%cas_shv ; ys%cas_co2 = y%cas_co2
         call surface_derivs(ys, fs, n, sf)
         enth1 = (wcap*y%cas_enthalpy + dt*(sf%src_enth + gah*fro%surf%enth_atm)) / (wcap + dt*gah)
         shv1  = (wcap*y%cas_shv      + dt*(sf%src_vap  + gaw*fro%surf%shv_atm )) / (wcap + dt*gaw)
      else
         call newton_surface_solve(y, fs, n, dt, wcap, gah, gaw, enth1, shv1, sf, nfeval, ok)
      end if
      y_out%cas_enthalpy = enth1
      y_out%cas_shv      = shv1
      y_out%cas_co2      = (ccap*y%cas_co2 + dt*(fro%surf%nee_biotic + gac*fro%surf%co2_atm)) / (ccap + dt*gac)

      !----- soil-heat column: implicit BE-Thomas (soil_energy_flux). ---------------------------!
      se%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      eforc%g_top = sf%g_top ; eforc%geothermal = fro%geothermal
      do k = 1_ik, nsl
         eforc%soil_water(k)     = y%theta(k)
         eforc%root_heat_sink(k) = sf%coh_qsoil * fro%soil%root_frac(k)
         eforc%w_flux(k)         = 0.0_wp
      end do
      !----- boundary water-enthalpy advection (guard-lift): infiltrating throughfall carries its liquid !
      !      enthalpy INTO the soil top, runoff/drainage carry it OUT (state^n temps, frozen -> the fixed !
      !      source telescopes exactly under the b-weighted ledger). Matches the split :436-439. root_    !
      !      heat_sink is a SINK, so q_src = -sink/dz: subtract an inflow, add an outflow. --------------!
      e_infil = fro%infiltration * internal_energy_liquid(fro%rain_temp)
      e_runof = fro%runoff_surf  * internal_energy_liquid(fro%surf%t_ground)   ! t_ground frozen @ state^n
      e_drain = fro%drainage     * internal_energy_liquid(fro%t_bot)
      eforc%root_heat_sink(1)   = eforc%root_heat_sink(1)   - e_infil + e_runof
      eforc%root_heat_sink(nsl) = eforc%root_heat_sink(nsl) + e_drain
      call soil_energy_flux(se, eforc, fro%therm, fro%soil, fro%energy_opts, dt, eflux)
      y_out%soil_energy(1:nsl) = se%soil_energy(1:nsl)

      !----- soil water is OPERATOR-SPLIT OUT of the ESDIRK stages: theta is PASSED THROUGH (held at the   !
      !      stage input = theta^n) and the AUTHORITATIVE end-of-step theta is committed once, from the     !
      !      scratch column_hydrology_flux (fro%theta1), in column_fast_step_ark. Re-solving it here with a !
      !      relief-free single-BE Richards drifted to saturation over long wet runs (no ponding/runoff),   !
      !      then hung the next scratch solve; the robust ponding/runoff/free-drain solve is the SOLE       !
      !      soil-water authority now (the ED2 "single soil-water authority" principle). theta feeds only   !
      !      the t_ground diagnosis + the soil-energy thermal property above, both correctly at theta^n. ---!
      y_out%theta(1:nsl) = y%theta(1:nsl)

      !----- plant hydraulics PASSED THROUGH (advanced by advance_hydraulics_full, not here). ----!
      y_out%psi(:, 1:n) = y%psi(:, 1:n)

      !----- emit this stage's boundary-flux RATES for the conservation ledger (§2.3). The b-weight   !
      !      + cross-substep accumulation happens in ark2_column_step / adaptive_ark_march. Every      !
      !      quantity is the committed-state flux, so the accumulated amounts telescope to closure.    !
      if (present(bf)) then
         t_cas1 = cas_temp_of_enthalpy(enth1, shv1)         ! committed CAS temp for the dew liquid enthalpy
         associate (fs2 => fro%surf)
            bf%cas_enth_in  = sf%src_enth + gah*fs2%enth_atm    ; bf%cas_enth_out = gah*enth1
            bf%cas_vap_in   = sf%src_vap  + gaw*fs2%shv_atm     ; bf%cas_vap_out  = gaw*shv1
            bf%cas_co2_in   = fs2%nee_biotic + gac*fs2%co2_atm  ; bf%cas_co2_out  = gac*y_out%cas_co2
            bf%soil_enth_in = sf%g_top + fro%geothermal + e_infil
            bf%soil_enth_out= sf%coh_qsoil * sum(fro%soil%root_frac(1:nsl)) + e_runof + e_drain
            !----- soil water is out of the ARK: its storage delta + q_top/drainage/uptake fluxes are     !
            !      re-sourced once/step from the frozen hflux in column_fast_step_ark, so the per-stage    !
            !      bf carries ONLY the CAS-vapour exchange (drainage/runoff/precip are frozen fast-step).  !
            bf%soil_wat_in  = 0.0_wp                            ; bf%soil_wat_out = 0.0_wp
            !----- condensation (dew) leaves the CAS as liquid at Tcas -> a whole-column water + liquid-   !
            !      enthalpy OUTPUT (the CAS-side loss is already in src_vap/src_enth, so cas_water/energy   !
            !      close automatically). -----------------------------------------------------------------!
            bf%whole_enth_in= sf%coh_rnet + fs2%abs_sw_ground + fs2%abs_lw_ground + e_infil
            bf%whole_enth_out= gah*(enth1 - fs2%enth_atm) + e_runof + e_drain                             &
                             + sf%cond*internal_energy_liquid(t_cas1)
            bf%whole_wat_in = 0.0_wp                            ; bf%whole_wat_out = gaw*(shv1 - fs2%shv_atm) + sf%cond
         end associate
      end if
   end subroutine column_be_stage

   !---------------------------------------------------------------------------------------!
   ! imex_euler_column_step -- the COMPLETE gamma=1, first-order IMEX-Euler column integrator = one    !
   ! column_be_stage (CAS + soil) composed with the exact-exponential plant-hydraulics operator split  !
   ! over the full dt (advance_hydraulics_full). It is L-stable at the full dt_fast where the explicit  !
   ! RK4 oracle blows up, and is the gamma=1 member of the ARK family / the P2 baseline. Hydraulics is  !
   ! an EXPLICIT operator split (not folded into the BE stage) because solve_plant_water is an exact    !
   ! matrix exponential, not a backward-Euler stage. ark2_column_step composes the same two pieces at   !
   ! 2nd order (two column_be_stage calls + one hydraulics split), so there is ONE hydraulics path.     !
   !---------------------------------------------------------------------------------------!
   subroutine imex_euler_column_step(y, fro, n, nsl, dt, y_out, niter, relax)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      integer(ik), optional, intent(in)  :: niter
      real(wp),    optional, intent(in)  :: relax

      call column_be_stage(y, fro, n, nsl, dt, y_out, niter, relax)         ! CAS + soil (psi passed through)
      call advance_hydraulics_full(y, fro, n, nsl, dt, y_out)               ! exact-exp hydraulics over full dt
   end subroutine imex_euler_column_step

   !---------------------------------------------------------------------------------------!
   ! newton_surface_solve -- the ARROWHEAD: a direct 2x2 Newton solve of the coupled backward-Euler   !
   ! CAS surface block { R_H, R_q } = 0 for (cas_enthalpy H1, cas_shv q1), where the surface sources   !
   ! src_enth/src_vap depend nonlinearly on (H1,q1) through tcas, qcas and qsat. Replaces the leaf<->  !
   ! CAS Picard iteration: quadratic convergence (~1 step) + robust near saturation. Numerical         !
   ! Jacobian (finite-difference surface_derivs) -- captures the strong VPD self-limiting d src_vap/dq  !
   ! with no derivation risk. Singular-Jacobian guard + line search + supersaturation clamp + eval cap; !
   ! never error stops (GPU-safe). Commits the CAS via the FLUX form so budgets close for ANY sf.      !
   !---------------------------------------------------------------------------------------!
   subroutine newton_surface_solve(y, fs, n, dt, wcap, gah, gaw, enth1, shv1, sf, nfeval, ok)
      type(column_state_t),   intent(in)    :: y
      type(surface_frozen_t), intent(in)    :: fs
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt, wcap, gah, gaw
      real(wp),               intent(out)   :: enth1, shv1
      type(surface_tend_t),   intent(out)   :: sf
      integer(ik),            intent(out)   :: nfeval
      logical,                intent(out)   :: ok

      type(surface_state_t) :: ys
      real(wp)    :: H0, q0, Hk, qk, R_H, R_q, J11, J12, J21, J22, detJ, delH, delq
      real(wp)    :: lam, rn0, Ht, qt, RHt, Rqt
      integer(ik) :: it, ls
      real(wp),    parameter :: RTOL_N = 1.0e-7_wp, ATOL_H = 5.0e1_wp, ATOL_Q = 1.0e-6_wp
      real(wp),    parameter :: DETEPS = 1.0e-30_wp
      integer(ik), parameter :: NEWT_MAX = 4_ik, LS_MAX = 6_ik, FEVAL_CAP = 24_ik

      H0 = y%cas_enthalpy ; q0 = y%cas_shv
      Hk = H0 ; qk = q0 ; nfeval = 0_ik ; ok = .false.
      ys%cas_co2 = y%cas_co2
      ys%cas_enthalpy = Hk ; ys%cas_shv = qk
      call surface_derivs(ys, fs, n, sf) ; nfeval = nfeval + 1_ik
      R_H = wcap*(Hk - H0)/dt - sf%src_enth - gah*(fs%enth_atm - Hk)
      R_q = wcap*(qk - q0)/dt - sf%src_vap  - gaw*(fs%shv_atm  - qk)

      do it = 1_ik, NEWT_MAX
         if ( abs(R_H)*dt/wcap <= ATOL_H + RTOL_N*abs(Hk) .and.                                  &
              abs(R_q)*dt/wcap <= ATOL_Q + RTOL_N*abs(qk) ) then
            ok = .true. ; exit
         end if
         call jac_surface(Hk, qk, y%cas_co2, fs, sf, n, wcap, gah, gaw, dt, J11, J12, J21, J22, nfeval)
         detJ = J11*J22 - J12*J21
         if (detJ <= DETEPS*abs(J11*J22) .or. detJ <= 0.0_wp) then       ! singular / sign-flipped guard
            delH = -R_H / max(J11, tiny_num)                             ! damped-diagonal (Picard-like) fallback
            delq = -R_q / max(J22, tiny_num)
         else
            delH = (-R_H*J22 + R_q*J12)/detJ
            delq = (-R_q*J11 + R_H*J21)/detJ
         end if
         lam = 1.0_wp ; rn0 = R_H*R_H + R_q*R_q ; Ht = Hk ; qt = qk ; RHt = R_H ; Rqt = R_q
         do ls = 1_ik, LS_MAX
            Ht = Hk + lam*delH ; qt = qk + lam*delq
            !----- NO supersaturation STATE clamp: ED2 never clamps can_shv to SAT*qsat -- doing so is a   !
            !      harsh state discontinuity that inflates the ARK2 embedded error in cas_shv AND cas_     !
            !      enthalpy near RH=1 and thrashes the adaptive controller. Like ED2 we TOLERATE transient !
            !      supersaturation; the smooth condensation SINK in surface_derivs relaxes it physically.  !
            ys%cas_enthalpy = Ht ; ys%cas_shv = qt
            call surface_derivs(ys, fs, n, sf) ; nfeval = nfeval + 1_ik
            RHt = wcap*(Ht - H0)/dt - sf%src_enth - gah*(fs%enth_atm - Ht)
            Rqt = wcap*(qt - q0)/dt - sf%src_vap  - gaw*(fs%shv_atm  - qt)
            if (RHt*RHt + Rqt*Rqt <= (1.0_wp - 1.0e-4_wp*lam)*rn0) exit          ! Armijo
            lam = 0.5_wp*lam
         end do
         Hk = Ht ; qk = qt ; R_H = RHt ; R_q = Rqt
         if (nfeval >= FEVAL_CAP) exit
      end do

      !----- authoritative final eval + flux-form commit (conservation holds for ANY sf). -----------!
      ys%cas_enthalpy = Hk ; ys%cas_shv = qk
      call surface_derivs(ys, fs, n, sf) ; nfeval = nfeval + 1_ik
      enth1 = (wcap*H0 + dt*(sf%src_enth + gah*fs%enth_atm)) / (wcap + dt*gah)
      shv1  = (wcap*q0 + dt*(sf%src_vap  + gaw*fs%shv_atm )) / (wcap + dt*gaw)
   end subroutine newton_surface_solve

   !----- 2x2 numerical Jacobian of (R_H, R_q) w.r.t. (H, q) by forward-differencing surface_derivs. --!
   subroutine jac_surface(Hk, qk, co2, fs, sf, n, wcap, gah, gaw, dt, J11, J12, J21, J22, nfeval)
      real(wp),               intent(in)    :: Hk, qk, co2, wcap, gah, gaw, dt
      type(surface_frozen_t), intent(in)    :: fs
      type(surface_tend_t),   intent(in)    :: sf         ! base eval at (Hk,qk)
      integer(ik),            intent(in)    :: n
      real(wp),               intent(out)   :: J11, J12, J21, J22
      integer(ik),            intent(inout) :: nfeval
      type(surface_state_t) :: ys
      type(surface_tend_t)  :: sfp
      real(wp) :: dH, dq, dse_dH, dsv_dH, dse_dq, dsv_dq
      real(wp), parameter :: SQEPS = 1.4901161e-8_wp, HSCALE = 1.0e4_wp, QSCALE = 1.0e-3_wp
      dH = SQEPS * max(abs(Hk), HSCALE)
      dq = SQEPS * max(abs(qk), QSCALE)
      ys%cas_co2 = co2
      ys%cas_enthalpy = Hk + dH ; ys%cas_shv = qk
      call surface_derivs(ys, fs, n, sfp) ; nfeval = nfeval + 1_ik
      dse_dH = (sfp%src_enth - sf%src_enth)/dH ; dsv_dH = (sfp%src_vap - sf%src_vap)/dH
      ys%cas_enthalpy = Hk ; ys%cas_shv = qk + dq
      call surface_derivs(ys, fs, n, sfp) ; nfeval = nfeval + 1_ik
      dse_dq = (sfp%src_enth - sf%src_enth)/dq ; dsv_dq = (sfp%src_vap - sf%src_vap)/dq
      J11 = wcap/dt + gah - dse_dH ; J12 =              - dse_dq
      J21 =              - dsv_dH  ; J22 = wcap/dt + gaw - dsv_dq
   end subroutine jac_surface

   !---------------------------------------------------------------------------------------!
   ! rk4_column_step -- one classical 4th-order Runge-Kutta step of the whole column state over the  !
   ! pure RHS column_derivs, with the frozen forcing `fro` held constant across the four stages (the  !
   ! explicit part of the additive split). Commits into y_out; y is unchanged.                        !
   !---------------------------------------------------------------------------------------!
   subroutine rk4_column_step(y, fro, n, nsl, dt, y_out, freeze_theta)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      logical, optional,     intent(in)  :: freeze_theta   !< zero the soil-water tendency (theta held fixed):
                                                           !< makes the oracle solve the SAME reduced system the
                                                           !< ARK stepper does (soil water operator-split OUT).

      type(column_tend_t)  :: k1, k2, k3, k4
      type(column_state_t) :: ys
      logical              :: frz

      frz = .false. ; if (present(freeze_theta)) frz = freeze_theta

      call column_derivs(y, fro, n, nsl, k1) ; if (frz) k1%dtheta_dt = 0.0_wp
      call state_axpy(y, 0.5_wp * dt, k1, n, nsl, ys) ; call column_derivs(ys, fro, n, nsl, k2) ; if (frz) k2%dtheta_dt = 0.0_wp
      call state_axpy(y, 0.5_wp * dt, k2, n, nsl, ys) ; call column_derivs(ys, fro, n, nsl, k3) ; if (frz) k3%dtheta_dt = 0.0_wp
      call state_axpy(y,          dt, k3, n, nsl, ys) ; call column_derivs(ys, fro, n, nsl, k4) ; if (frz) k4%dtheta_dt = 0.0_wp

      !----- y_out = y + dt/6 (k1 + 2 k2 + 2 k3 + k4). -------------------------------------!
      call state_init(y, n, nsl, y_out)
      call state_accum(y_out, dt / 6.0_wp, k1, n, nsl)
      call state_accum(y_out, dt / 3.0_wp, k2, n, nsl)
      call state_accum(y_out, dt / 3.0_wp, k3, n, nsl)
      call state_accum(y_out, dt / 6.0_wp, k4, n, nsl)
   end subroutine rk4_column_step

   !----- ys = y + a*k  (state + a * tendency). ------------------------------------------!
   pure subroutine state_axpy(y, a, k, n, nsl, ys)
      type(column_state_t), intent(in)  :: y
      real(wp),             intent(in)  :: a
      type(column_tend_t),  intent(in)  :: k
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: ys
      integer(ik) :: j, i
      ys%cas_enthalpy = y%cas_enthalpy + a * k%d_cas_enthalpy
      ys%cas_shv      = y%cas_shv      + a * k%d_cas_shv
      ys%cas_co2      = y%cas_co2      + a * k%d_cas_co2
      ys%soil_energy  = y%soil_energy
      ys%theta        = y%theta
      do j = 1_ik, nsl
         ys%soil_energy(j) = y%soil_energy(j) + a * k%dedt(j)
         ys%theta(j)       = y%theta(j)       + a * k%dtheta_dt(j)
      end do
      allocate(ys%psi(N_HYDRO, n))
      do i = 1_ik, n
         ys%psi(:, i) = y%psi(:, i) + a * k%dpsi_dt(:, i)
      end do
   end subroutine state_axpy

   !----- copy the prognostic state (used to seed the RK combination). --------------------!
   pure subroutine state_init(y, n, nsl, ys)
      type(column_state_t), intent(in)  :: y
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: ys
      ys%cas_enthalpy = y%cas_enthalpy ; ys%cas_shv = y%cas_shv ; ys%cas_co2 = y%cas_co2
      ys%soil_energy  = y%soil_energy  ; ys%theta   = y%theta
      allocate(ys%psi(N_HYDRO, n))
      ys%psi(:, 1:n) = y%psi(:, 1:n)
   end subroutine state_init

   !----- ys += a*k  (accumulate a weighted tendency into a state). -----------------------!
   pure subroutine state_accum(ys, a, k, n, nsl)
      type(column_state_t), intent(inout) :: ys
      real(wp),             intent(in)    :: a
      type(column_tend_t),  intent(in)    :: k
      integer(ik),          intent(in)    :: n, nsl
      integer(ik) :: j, i
      ys%cas_enthalpy = ys%cas_enthalpy + a * k%d_cas_enthalpy
      ys%cas_shv      = ys%cas_shv      + a * k%d_cas_shv
      ys%cas_co2      = ys%cas_co2      + a * k%d_cas_co2
      do j = 1_ik, nsl
         ys%soil_energy(j) = ys%soil_energy(j) + a * k%dedt(j)
         ys%theta(j)       = ys%theta(j)       + a * k%dtheta_dt(j)
      end do
      do i = 1_ik, n
         ys%psi(:, i) = ys%psi(:, i) + a * k%dpsi_dt(:, i)
      end do
   end subroutine state_accum

   !---------------------------------------------------------------------------------------!
   ! ark2_column_step -- one 2nd-order L-stable IMEX step via the ARS(2,2,2) additive Runge-Kutta      !
   ! (Ascher-Ruuth-Spiteri 1997, Appl.Numer.Math. 25:151; identical gamma in Giraldo et al. 2013     !
   ! "ARK2"). Stiffly accurate on both tableaux, so the two ESDIRK stages map onto imex_euler's       !
   ! "last BE solve = committed state" structure with dt -> gamma*dt. The biotic CO2 source is folded !
   ! IMPLICIT (stays in the CO2 BE numerator), so f_E == 0 and the scheme reduces to a clean 2-solve  !
   ! ESDIRK2. PLANT HYDRAULICS IS OPERATOR-SPLIT OUT of the tableau (solve_plant_water is an EXACT     !
   ! matrix exponential, NOT a backward-Euler stage -- putting it in the ESDIRK accumulation silently  !
   ! drops the order and overshoots psi): psi is frozen through the stages, then advanced once over    !
   ! the full dt, and excluded from the embedded error. y_err = (Y3-base3)-(Y2-y_n) is the free        !
   ! embedded 1st-order estimate for the adaptive controller (2 solves/step vs step-doubling's 3).     !
   !---------------------------------------------------------------------------------------!
   subroutine ark2_column_step(y, fro, n, nsl, dt, y_out, y_err, niter, relax, bf)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out, y_err
      integer(ik), optional, intent(in)  :: niter
      real(wp),    optional, intent(in)  :: relax
      type(column_bflux_t), optional, intent(out) :: bf   !< b-weighted boundary-flux AMOUNTS over dt (ledger)
      real(wp), parameter :: GAMMA = 0.2928932188134524_wp   ! 1 - 1/sqrt(2)
      real(wp), parameter :: BETA  = 2.4142135623730951_wp   ! (1-gamma)/gamma = 1 + sqrt(2)
      type(column_state_t) :: Y2, base3, Y3
      type(stage_bflux_t)  :: bf2, bf3
      integer(ik) :: np
      real(wp)    :: rlx
      np = 1_ik ; if (present(niter)) np = max(1_ik, niter)
      rlx = 1.0_wp ; if (present(relax)) rlx = relax

      !----- Stage 2: gamma*dt BE stage from y_n (CAS+soil only; psi frozen -- it is split out). -----!
      call column_be_stage(y, fro, n, nsl, GAMMA*dt, Y2, niter=np, relax=rlx, bf=bf2)
      !----- Stage 3: extrapolated base. The BETA=2.414 extrapolation can overshoot BOTH the vG theta   !
      !      range AND the CAS enthalpy into a wild temperature where qsat(T) overflows to NaN; clamp     !
      !      both to physical ranges so the stage stays FINITE. This only bites on a genuinely oversized  !
      !      step (which the adaptive controller then rejects and shrinks normally) -- in range it is an  !
      !      identity, so no accuracy cost. Without the CAS clamp a big transient poisons the whole march !
      !      with NaN. --------------------------------------------------------------------------------!
      call state_extrap(y, BETA, Y2, n, nsl, base3)
      call clamp_theta(base3, fro, nsl)
      call clamp_cas(base3)
      call column_be_stage(base3, fro, n, nsl, GAMMA*dt, Y3, niter=np, relax=rlx, bf=bf3)
      call state_init(Y3, n, nsl, y_out)
      !----- operator-split hydraulics: exact 2x2 over the FULL dt from y_n, endpoint transp. -------!
      call advance_hydraulics_full(y, fro, n, nsl, dt, y_out)
      !----- embedded 1st-order error estimate (psi is split out -> zeroed). -------------------------!
      call state_err_diff(Y3, base3, Y2, y, n, nsl, y_err)
      !----- b-weighted boundary-flux amounts: b^I = (0, 1-gamma, gamma) -> exact telescoping.         !
      !      (Water closure is exact only when clamp_theta is inactive; it barely moves over gamma*dt.) !
      if (present(bf)) call bflux_bweight(bf, bf2, bf3, dt, GAMMA)
   end subroutine ark2_column_step

   !----- ledger helpers: b-weight two stage RATE structs into accumulated AMOUNTS over dt (weights   !
   !      b^I = (1-gamma, gamma) times dt); zero an accumulator; add one substep's amounts. -----------!
   pure subroutine bflux_bweight(acc, s2, s3, dt, gam)
      type(column_bflux_t), intent(out) :: acc
      type(stage_bflux_t),  intent(in)  :: s2, s3
      real(wp),             intent(in)  :: dt, gam
      real(wp) :: b2, b3
      b2 = (1.0_wp - gam) * dt ; b3 = gam * dt
      acc%cas_enth_in   = b2*s2%cas_enth_in   + b3*s3%cas_enth_in
      acc%cas_enth_out  = b2*s2%cas_enth_out  + b3*s3%cas_enth_out
      acc%cas_vap_in    = b2*s2%cas_vap_in    + b3*s3%cas_vap_in
      acc%cas_vap_out   = b2*s2%cas_vap_out   + b3*s3%cas_vap_out
      acc%cas_co2_in    = b2*s2%cas_co2_in    + b3*s3%cas_co2_in
      acc%cas_co2_out   = b2*s2%cas_co2_out   + b3*s3%cas_co2_out
      acc%soil_enth_in  = b2*s2%soil_enth_in  + b3*s3%soil_enth_in
      acc%soil_enth_out = b2*s2%soil_enth_out + b3*s3%soil_enth_out
      acc%soil_wat_in   = b2*s2%soil_wat_in   + b3*s3%soil_wat_in
      acc%soil_wat_out  = b2*s2%soil_wat_out  + b3*s3%soil_wat_out
      acc%whole_enth_in = b2*s2%whole_enth_in + b3*s3%whole_enth_in
      acc%whole_enth_out= b2*s2%whole_enth_out+ b3*s3%whole_enth_out
      acc%whole_wat_in  = b2*s2%whole_wat_in  + b3*s3%whole_wat_in
      acc%whole_wat_out = b2*s2%whole_wat_out + b3*s3%whole_wat_out
   end subroutine bflux_bweight

   pure subroutine bflux_zero(acc)
      type(column_bflux_t), intent(out) :: acc
      acc = column_bflux_t()
   end subroutine bflux_zero

   pure subroutine bflux_add(acc, s)
      type(column_bflux_t), intent(inout) :: acc
      type(column_bflux_t), intent(in)    :: s
      acc%cas_enth_in   = acc%cas_enth_in   + s%cas_enth_in
      acc%cas_enth_out  = acc%cas_enth_out  + s%cas_enth_out
      acc%cas_vap_in    = acc%cas_vap_in    + s%cas_vap_in
      acc%cas_vap_out   = acc%cas_vap_out   + s%cas_vap_out
      acc%cas_co2_in    = acc%cas_co2_in    + s%cas_co2_in
      acc%cas_co2_out   = acc%cas_co2_out   + s%cas_co2_out
      acc%soil_enth_in  = acc%soil_enth_in  + s%soil_enth_in
      acc%soil_enth_out = acc%soil_enth_out + s%soil_enth_out
      acc%soil_wat_in   = acc%soil_wat_in   + s%soil_wat_in
      acc%soil_wat_out  = acc%soil_wat_out  + s%soil_wat_out
      acc%whole_enth_in = acc%whole_enth_in + s%whole_enth_in
      acc%whole_enth_out= acc%whole_enth_out+ s%whole_enth_out
      acc%whole_wat_in  = acc%whole_wat_in  + s%whole_wat_in
      acc%whole_wat_out = acc%whole_wat_out + s%whole_wat_out
   end subroutine bflux_add

   !----- out = (1-b)*y + b*Y2  (the ARS stage-3 extrapolation base). --------------------------!
   pure subroutine state_extrap(y, b, Y2, n, nsl, out)
      type(column_state_t), intent(in)  :: y, Y2
      real(wp),             intent(in)  :: b
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: out
      real(wp)    :: a
      integer(ik) :: k, i
      a = 1.0_wp - b
      out%cas_enthalpy = a*y%cas_enthalpy + b*Y2%cas_enthalpy
      out%cas_shv      = a*y%cas_shv      + b*Y2%cas_shv
      out%cas_co2      = a*y%cas_co2      + b*Y2%cas_co2
      out%soil_energy = y%soil_energy ; out%theta = y%theta
      do k = 1_ik, nsl
         out%soil_energy(k) = a*y%soil_energy(k) + b*Y2%soil_energy(k)
         out%theta(k)       = a*y%theta(k)       + b*Y2%theta(k)
      end do
      allocate(out%psi(N_HYDRO, n))
      do i = 1_ik, n
         out%psi(:, i) = a*y%psi(:, i) + b*Y2%psi(:, i)     ! == y%psi (psi frozen in the stages)
      end do
   end subroutine state_extrap

   !----- clamp the extrapolated CAS enthalpy + humidity into a wide PHYSICAL range so a BETA=2.414   !
   !      overshoot cannot drive cas_temp_of_enthalpy to a wild T where qsat(T) overflows to NaN. Only  !
   !      active on a pathological overshoot (then the step is rejected); an in-range base3 is untouched.!
   pure subroutine clamp_cas(s)
      type(column_state_t), intent(inout) :: s
      real(wp) :: t, shv_c
      real(wp), parameter :: T_LO = 180.0_wp, T_HI = 350.0_wp, SHV_LO = 1.0e-8_wp, SHV_HI = 0.06_wp
      shv_c = min(max(s%cas_shv, SHV_LO), SHV_HI)
      t     = cas_temp_of_enthalpy(s%cas_enthalpy, shv_c)
      t     = min(max(t, T_LO), T_HI)
      s%cas_shv      = shv_c
      s%cas_enthalpy = cas_enthalpy_of_temp(t, shv_c)
   end subroutine clamp_cas

   !----- clamp the extrapolated theta into [theta_res, theta_sat] (van Genuchten domain). -----!
   pure subroutine clamp_theta(s, fro, nsl)
      type(column_state_t),  intent(inout) :: s
      type(column_frozen_t), intent(in)    :: fro
      integer(ik),           intent(in)    :: nsl
      integer(ik) :: k
      do k = 1_ik, nsl
         s%theta(k) = min(max(s%theta(k), fro%soil%theta_res(k)), fro%soil%theta_sat(k))
      end do
   end subroutine clamp_theta

   !----- err = (Y3 - base3) - (Y2 - y)  (the embedded 2nd-1st order difference); psi zeroed. ---!
   pure subroutine state_err_diff(Y3, base3, Y2, y, n, nsl, err)
      type(column_state_t), intent(in)  :: Y3, base3, Y2, y
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: err
      integer(ik) :: k
      err%cas_enthalpy = (Y3%cas_enthalpy - base3%cas_enthalpy) - (Y2%cas_enthalpy - y%cas_enthalpy)
      err%cas_shv      = (Y3%cas_shv      - base3%cas_shv)      - (Y2%cas_shv      - y%cas_shv)
      err%cas_co2      = (Y3%cas_co2      - base3%cas_co2)      - (Y2%cas_co2      - y%cas_co2)
      err%soil_energy = 0.0_wp ; err%theta = 0.0_wp
      do k = 1_ik, nsl
         err%soil_energy(k) = (Y3%soil_energy(k) - base3%soil_energy(k)) - (Y2%soil_energy(k) - y%soil_energy(k))
         err%theta(k)       = (Y3%theta(k)       - base3%theta(k))       - (Y2%theta(k)       - y%theta(k))
      end do
      allocate(err%psi(N_HYDRO, n))
      err%psi(:, 1:n) = 0.0_wp
   end subroutine state_err_diff

   !----- operator-split plant hydraulics: exact 2x2 matrix-exp over the FULL dt from y%psi, driven   !
   !      by the ENDPOINT transpiration demand (surface_derivs at the committed CAS). ---------------!
   subroutine advance_hydraulics_full(y, fro, n, nsl, dt, y_out)
      type(column_state_t),  intent(in)    :: y
      type(column_frozen_t), intent(in)    :: fro
      integer(ik),           intent(in)    :: n, nsl
      real(wp),              intent(in)    :: dt
      type(column_state_t),  intent(inout) :: y_out
      type(surface_state_t)  :: ys
      type(surface_frozen_t) :: fs
      type(surface_tend_t)   :: sf
      type(hydro_env_t)      :: henv
      type(hydro_flux_t)     :: hfx
      real(wp)    :: tg, fl, psi_i(N_HYDRO)
      integer(ik) :: i
      call uext_to_temp(y_out%soil_energy(1), y_out%theta(1)*rho_h2o,                             &
                        fro%therm%soil_dry_heat_capacity(1), tg, fl)
      fs = fro%surf ; fs%t_ground = tg
      ys%cas_enthalpy = y_out%cas_enthalpy ; ys%cas_shv = y_out%cas_shv ; ys%cas_co2 = y_out%cas_co2
      call surface_derivs(ys, fs, n, sf)
      do i = 1_ik, n
         henv%transp     = sf%transp_c(i) * fro%surf%src_frac / max(fro%nplant(i), tiny_num)
         henv%soil_psi   = fro%soil_psi_root ; henv%rhizo_cond = fro%rhizo_cond
         henv%bleaf      = fro%bleaf(i) ; henv%bsap = fro%bsap(i) ; henv%broot = fro%broot(i)
         henv%sap_area   = fro%sap_area(i) ; henv%height = fro%height(i) ; henv%leaf_area = fro%leaf_area(i)
         psi_i = y%psi(:, i)
         call solve_plant_water(henv, fro%hydro_p, fro%hydro_o, dt, psi_i, hfx)
         y_out%psi(:, i) = psi_i
      end do
   end subroutine advance_hydraulics_full

   !----- out = a - b  (state difference; used to form the low-order embedded solution). --------!
   pure subroutine state_sub(a, b, n, nsl, out)
      type(column_state_t), intent(in)  :: a, b
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: out
      integer(ik) :: k
      out%cas_enthalpy = a%cas_enthalpy - b%cas_enthalpy
      out%cas_shv      = a%cas_shv      - b%cas_shv
      out%cas_co2      = a%cas_co2      - b%cas_co2
      out%soil_energy = a%soil_energy ; out%theta = a%theta
      do k = 1_ik, nsl
         out%soil_energy(k) = a%soil_energy(k) - b%soil_energy(k)
         out%theta(k)       = a%theta(k)       - b%theta(k)
      end do
      allocate(out%psi(N_HYDRO, n))
      out%psi(:, 1:n) = a%psi(:, 1:n) - b%psi(:, 1:n)
   end subroutine state_sub

   !---------------------------------------------------------------------------------------!
   ! adaptive_ark_march -- integrate to t_end with the ARK2 embedded error estimate driving the       !
   ! step controller (2 solves/step, no step-doubling). y_err (already the 2nd-1st order difference)   !
   ! is the local error; the WRMS of it vs tolerance drives accept/reject via adaptive_step_update     !
   ! (p=1 embedded -> exponent -1/2). Reports the step + reject count.                                 !
   !---------------------------------------------------------------------------------------!
   subroutine adaptive_ark_march(y0, fro, n, nsl, t_end, rtol, dt_init, y_out, nsteps, nrej, niter, relax, acc)
      type(column_state_t),  intent(in)  :: y0
      type(column_frozen_t), intent(in)  :: fro
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: t_end, rtol, dt_init
      type(column_state_t),  intent(out) :: y_out
      integer(ik),           intent(out) :: nsteps, nrej
      integer(ik), optional, intent(in)  :: niter    !< coupled leaf<->CAS Newton cap (default 8)
      real(wp),    optional, intent(in)  :: relax    !< under-relaxation (default 0.6)
      type(column_bflux_t), optional, intent(out) :: acc  !< accumulated boundary-flux amounts (ledger)

      type(column_state_t) :: y, y_new, y_err, y_lo
      type(column_bflux_t) :: bfsub
      real(wp)             :: t, dt, err, fac, rlx, dt_floor
      real(wp), parameter  :: SAFETY = 0.9_wp, FMIN = 0.2_wp, FMAX = 5.0_wp
      integer(ik)          :: np

      np = 8_ik ; if (present(niter)) np = max(1_ik, niter)
      rlx = 0.6_wp ; if (present(relax)) rlx = relax
      if (present(acc)) call bflux_zero(acc)
      !----- substep FLOOR: bound the worst case to ~t_end/DT_FLOOR sub-steps. The ARK2 BE stages are    !
      !      L-stable, so a floor step is STABLE (bounded) even when the embedded error stays above tol   !
      !      -- e.g. a stiff transient the tolerance can't resolve. A tiny absolute floor (the old 1e-2s) !
      !      let a pathological step balloon to ~1.8e5 sub-steps and stall the march; t_end/64 caps it at !
      !      64 and degrades gracefully. (Also surfaces a genuine non-finite state promptly rather than   !
      !      grinding at the floor forever.) -----------------------------------------------------------!
      dt_floor = max(1.0e-2_wp, t_end / 64.0_wp)

      call state_init(y0, n, nsl, y)
      t = 0.0_wp ; dt = min(dt_init, t_end) ; nsteps = 0_ik ; nrej = 0_ik
      do
         if (t >= t_end - tiny_num) exit
         dt = min(dt, t_end - t)
         call ark2_column_step(y, fro, n, nsl, dt, y_new, y_err, niter=np, relax=rlx, bf=bfsub)
         call state_sub(y_new, y_err, n, nsl, y_lo)               ! the 1st-order embedded solution
         err = state_wrms(y_new, y_lo, y, n, nsl, rtol)           ! = WRMS(y_err)
         !----- ROBUSTNESS: a non-finite err (a stage -- typically the BETA=2.414 stage-3 extrapolation    !
         !      base3 -- overshot the CAS enthalpy into a region where qsat(T) overflows) is a step that   !
         !      is simply TOO BIG: REJECT it and shrink dt deterministically (the NaN poisons the adaptive !
         !      fac, so use FMIN directly). At a smaller dt, Y2 ~ y and base3 no longer overshoots, so the !
         !      step becomes finite and the march recovers -- the correct adaptive response, not a force-  !
         !      accept. Only if even a floor-sized step is non-finite do we commit + bail so meds_main's    !
         !      has_nan check reports it cleanly instead of the march hanging.                              !
         if (err /= err .or. dt /= dt) then
            if (dt <= dt_floor) then
               call state_init(y_new, n, nsl, y) ; t = t + dt_floor ; nsteps = nsteps + 1_ik ; exit
            end if
            nrej = nrej + 1_ik ; dt = max(dt * FMIN, dt_floor) ; cycle
         end if
         fac = adaptive_step_update(max(err, tiny_num), SAFETY, FMIN, FMAX)
         if (err <= 1.0_wp .or. dt <= dt_floor) then
            call state_init(y_new, n, nsl, y)
            if (present(acc)) call bflux_add(acc, bfsub)          ! accumulate ONLY accepted substeps
            t = t + dt ; nsteps = nsteps + 1_ik
            dt = dt * fac
         else
            nrej = nrej + 1_ik
            dt = dt * fac
         end if
         if (nsteps + nrej > 4096_ik) exit                        ! hard backstop (should never trigger)
      end do
      call state_init(y, n, nsl, y_out)
   end subroutine adaptive_ark_march

end module meds_ark_stepper
