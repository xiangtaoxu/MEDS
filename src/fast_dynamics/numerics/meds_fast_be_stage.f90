!==========================================================================================!
! meds_fast_be_stage -- the IMPLICIT STAGE of the fast-loop march, and the water-mass        !
! advance that rides with it.                                                                !
!                                                                                          !
! One ESDIRK stage is a backward-Euler solve of the coupled surface system: the canopy air's !
! three twins against the leaf, wood and ground temperatures, closed by a small Newton        !
! iteration over an analytic Jacobian. That machinery is independent of WHICH tableau calls    !
! it, which is the point of pulling it out: `meds_fast_ark` keeps the tableau and the march,   !
! and the RK4 oracle reaches the same stage solver directly instead of importing it from ARK.  !
!                                                                                          !
! Contents: `column_be_stage` (the stage), `newton_surface_solve` + `jac_surface` +            !
! `cas_box_commit` (its inner solve), and `advance_water_mass_full` / `advance_surf_water_full` !
! (the mass-form water advance a stage commits alongside the energy solve).                     !
!==========================================================================================!
module meds_fast_be_stage
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, cp_air, rho_h2o
   use meds_plant_hydraulics, only : solve_plant_water_batch
   use meds_hydr_lib, only : psi_from_water_content, water_content
   use meds_config, only : INTEG_ARK, CTRL_L2_STRICT
   use meds_fast_control, only : step_control_factor
   use meds_soil_types, only : energy_forcing_t, energy_flux_t, snow_melt_t
   use meds_column_state_types, only : soil_energy_column_t
   use meds_fast_time_derivs, only : surface_derivs, cas_conductances
   use meds_cas_biophysics, only : cas_column_step_implicit, cas_column_t, cas_source_t
   use meds_column_state_ops, only : state_init, state_sub, bflux_zero, bflux_add, bflux_bweight, clamp_cas, clamp_theta, &
                                     clamp_soil_energy, soil_water_store, soil_energy_store, plant_water_store, &
                                     canopy_film_store, deposit_condensate, clamp_canopy_film, unpack_column_state, &
                                     diagnose_soil_temps, assemble_soil_energy_forcing, apply_process_mask
   use meds_fast_types, only : column_budget_t, alloc_column_cohort, column_state_t, column_frozen_t, surface_state_t, &
                               cas_boundary_t, surface_tend_t, stage_bflux_t, column_bflux_t, error_control_t, column_tend_t, &
                               mask_is_full
   use meds_soil_energy, only : soil_energy_step_implicit
   use meds_cas_biophysics, only : cas_column_t, cas_source_t, cas_column_step_implicit
   use meds_ground_biophysics, only : snow_accumulate, snow_drain_meltwater, snow_cover_fraction, ground_surface_fluxes
   use meds_plant_types, only : N_HYDRO, NODE_LEAF, NODE_WOOD
   use meds_plant_hydraulics, only : solve_plant_water_batch
   use meds_therm_lib, only : cas_temp_of_enthalpy, internal_energy_liquid, internal_energy_to_temp, internal_energy_ice, &
                              temp_of_liquid_enthalpy
   use meds_budget_check, only : budget_check, budget_energy_rate_floor, budget_water_rate_floor, budget_co2_rate_floor
   implicit none
   private

   public :: column_be_stage, advance_water_mass_full, advance_surf_water_full

contains

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
   subroutine column_be_stage(y, frozen, n, nsl, dt, y_out, niter, bf, sf_out)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      type(column_state_t),  intent(out) :: y_out
      integer(ik), optional, intent(in)  :: niter    !< 1 = uncoupled BE baseline; >1 = coupled leaf<->CAS Newton
      type(stage_bflux_t), optional, intent(out) :: bf  !< per-stage boundary-flux RATES for the ARK ledger
      !----- this stage's OWN surface tendencies (incl. transp_c), for the caller to b-weight into      !
      !      the plant water-mass update (MEDS_ED2_RK45_DESIGN.md sec 1/4/5, P2) -- the SAME surf_tend the      !
      !      stage's own bf/CAS-source already used, so the mass debit and the CAS credit agree. --------!
      type(surface_tend_t), optional, intent(out) :: sf_out

      type(surface_state_t)      :: y_stage
      type(cas_boundary_t)       :: cas_stage   !< the CAS boundary with THIS stage's conductances
      type(surface_tend_t)       :: surf_tend
      type(soil_energy_column_t) :: se
      type(energy_forcing_t)     :: eforc
      type(energy_flux_t)        :: eflux
      real(wp)    :: t_ground, fliq1, wmass1, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2
      real(wp)    :: enth1, shv1, e_infil, e_drain, e_clip, e_floor, t_cas1, qloss_total
      real(wp)    :: co21, enth_unused, shv_unused
      integer(ik) :: k, np, nfeval
      logical     :: ok

      np = 1_ik ; if (present(niter)) np = max(1_ik, niter)

      call state_init(y, n, nsl, y_out)

      !----- diagnose the soil-top temperature so the ground skin sees the current state. --------!
      wmass1 = y%theta(1) * rho_h2o
      call internal_energy_to_temp(y%soil_energy(1), wmass1, frozen%params%therm%soil_dry_heat_capacity(1), t_ground, fliq1)

      cas_mass_capacity = frozen%cas%cas_mass_capacity ; cas_molar_capacity = frozen%cas%cas_molar_capacity
      cas_stage = frozen%cas              ! plain scalars: a cheap copy, overridden with this stage's conductances below

      !----- Re-solve the Monin-Obukhov surface layer at THIS STAGE's canopy-air state, so the    !
      !      ventilation the stage is charged for is the ventilation its own temperature earns.    !
      !      Everything downstream (the BE commit, the Newton, and the boundary-flux ledger) reads !
      !      the LOCAL g_atm_heat/g_atm_vapour/g_atm_co2, so refreshing them here keeps "one flux, both sides"           !
      !      automatically -- the state update and the ledger cannot disagree.                     !
      !                                                                                          !
      !      Mirror into cas_stage and CLEAR ITS mo_live: surface_derivs would otherwise re-solve the !
      !      surface layer on every one of the Newton's residual evaluations (up to 24 per stage)  !
      !      to fill a CAS tendency this scheme does not even read -- it commits the CAS through   !
      !      its own backward-Euler denominator.  So ARK pays for exactly ONE solve per stage. ----!
      call cas_conductances(frozen%cas, y%cas_enthalpy, y%cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      cas_stage%g_atm_heat = g_atm_heat ; cas_stage%g_atm_vapour = g_atm_vapour ; cas_stage%g_atm_co2 = g_atm_co2
      cas_stage%mo_live = .false.

      !----- CAS enthalpy + humidity. np==1: the uncoupled single-BE-pass baseline. np>1: a DIRECT 2x2  !
      !      Newton solve of the coupled backward-Euler surface block (the arrowhead). The FINAL surf_tend     !
      !      drives the soil sinks (single-flux-per-interface).                                         !
      if (np <= 1_ik) then
         y_stage%cas_enthalpy = y%cas_enthalpy ; y_stage%cas_shv = y%cas_shv ; y_stage%cas_co2 = y%cas_co2
         call surface_derivs(y_stage, cas_stage, frozen%tissue, frozen%film, frozen%ground, frozen%snow,    &
                             t_ground, n, surf_tend)
         call cas_box_commit(y%cas_enthalpy, y%cas_shv, y%cas_co2, surf_tend, cas_stage, cas_mass_capacity, cas_molar_capacity, &
                             g_atm_heat, g_atm_vapour, g_atm_co2, &
                             dt, enth1, shv1, co21)
      else
         call newton_surface_solve(y, cas_stage, frozen, t_ground, n, dt, cas_mass_capacity, g_atm_heat, g_atm_vapour, enth1, &
                                   shv1, surf_tend,   &
                                   nfeval, ok)
      end if
      y_out%cas_enthalpy = enth1
      y_out%cas_shv      = shv1
      if (np > 1_ik) call cas_box_commit(y%cas_enthalpy, y%cas_shv, y%cas_co2, surf_tend, cas_stage, cas_mass_capacity, &
          cas_molar_capacity,     &
                                         g_atm_heat, g_atm_vapour, g_atm_co2, dt, enth_unused, shv_unused, &
                                              co21)   ! CO2 rides the same box
      y_out%cas_co2      = co21
      if (present(sf_out)) sf_out = surf_tend

      !----- soil-heat column: implicit BE-Thomas (soil_energy_step_implicit). ---------------------------!
      se%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      !----- Root heat sink = qloss_total (uptake's advected enthalpy, sec 2/6, P2), the SAME sink     !
      !      column_derivs uses (meds_fast_time_derivs.f90), distributed by the static root_share       !
      !      profile: the soil pays once for the water the roots extract, the leaf/wood side gains it   !
      !      via qwflux_wl/q_wood_net, and the leaf pays the full vapour enthalpy of what it transpires. !
      !      (The old coh_qsoil proxy charged the soil a second time for that vapour's liquid part;      !
      !      2026-09 review, item 1A #10.) qloss_total sums to 0 when the P2 wiring is unset.            !
      !                                                                                                 !
      !      The ARK commits the SCRATCH hydrology's theta verbatim, so its faces (w_flux_frozen), its   !
      !      drainage and its two UNFACED post-solve mass corrections are the right numbers here: the    !
      !      clip (water leaving layer k for the pond, valued at the layer's state^n temperature) ADDS   !
      !      to the sink and the theta_res floor (water created in layer k) SUBTRACTS. Both are 0 unless !
      !      the hydrology actually corrected that layer. --------------------------------------------!
      qloss_total = sum(frozen%roots%qloss_frozen(1:n))
      e_infil = frozen%hydrology%infiltration * internal_energy_liquid(frozen%hydrology%t_infil)
      e_drain = frozen%hydrology%drainage     * internal_energy_liquid(frozen%hydrology%t_bot)
      e_clip  = sum(frozen%hydrology%clip_enth(1:nsl))
      e_floor = sum(frozen%hydrology%floor_enth(1:nsl))
      call assemble_soil_energy_forcing(eforc, nsl, surf_tend%g_top, frozen%hydrology%geothermal, y%theta,          &
                                        frozen%roots%root_share, qloss_total, frozen%hydrology%w_flux_frozen,             &
                                        frozen%hydrology%infiltration, frozen%hydrology%t_infil, e_drain,                     &
                                        sink_add=frozen%hydrology%clip_enth, sink_sub=frozen%hydrology%floor_enth)
      call soil_energy_step_implicit(se, eforc, frozen%params%therm, frozen%params%soil, frozen%params%energy_opts, dt, eflux)
      y_out%soil_energy(1:nsl) = se%soil_energy(1:nsl)

      !----- soil water is OPERATOR-SPLIT OUT of the ESDIRK stages: theta is PASSED THROUGH (held at the   !
      !      stage input = theta^n) and the AUTHORITATIVE end-of-step theta is committed once, from the     !
      !      scratch advance_soil_water_column (frozen%hydrology%theta1), in column_fast_step_ark. Re-solving it here with a !
      !      relief-free single-BE Richards drifted to saturation over long wet runs (no ponding/runoff),   !
      !      then hung the next scratch solve; the robust ponding/runoff/free-drain solve is the SOLE       !
      !      soil-water authority now (the ED2 "single soil-water authority" principle). theta feeds only   !
      !      the t_ground diagnosis + the soil-energy thermal property above, both correctly at theta^n. ---!
      y_out%theta(1:nsl) = y%theta(1:nsl)

      !----- pond PASSED THROUGH. y_out is intent(out), so without this it would default-initialise !
      !      to 0 rather than carry state^n -- harmless today (nothing reads it in a stage and it is   !
      !      excluded from the error norm) but wrong the moment #93 Phase 1 gives it a stage RHS. ----!
      y_out%w_surface      = y%w_surface
      y_out%w_surface_enth = y%w_surface_enth

      !----- plant water MASS PASSED THROUGH (advanced by advance_water_mass_full, not here). ----!
      y_out%leaf_water_mass(1:n) = y%leaf_water_mass(1:n)
      y_out%wood_water_mass(1:n) = y%wood_water_mass(1:n)

      !----- emit this stage's boundary-flux RATES for the conservation ledger (§2.3). The b-weight   !
      !      + cross-substep accumulation happens in ark2_column_step / adaptive_ark_march. Every      !
      !      quantity is the committed-state flux, so the accumulated amounts telescope to closure.    !
      if (present(bf)) then
         t_cas1 = cas_temp_of_enthalpy(enth1, shv1)         ! committed CAS temp for the dew liquid enthalpy
            bf%cas_enth_in  = surf_tend%src_enth + g_atm_heat*frozen%cas%enthalpy_atm    ; bf%cas_enth_out = g_atm_heat*enth1
            bf%cas_vap_in   = surf_tend%src_vap  + g_atm_vapour*frozen%cas%shv_atm     ; bf%cas_vap_out  = g_atm_vapour*shv1
            bf%cas_co2_in   = frozen%cas%nee_biotic + g_atm_co2*frozen%cas%co2_atm  ; bf%cas_co2_out  = g_atm_co2*y_out%cas_co2
            bf%soil_enth_in = surf_tend%g_top + frozen%hydrology%geothermal + e_infil + e_floor
            bf%soil_enth_out= qloss_total * sum(frozen%params%soil%root_frac(1:nsl)) + e_drain + e_clip
            !----- soil water is out of the ARK: its storage delta + q_top/drainage/uptake fluxes are     !
            !      re-sourced once/step from the frozen hflux in column_fast_step_ark, so the per-stage    !
            !      bf carries ONLY the CAS-vapour exchange (drainage/runoff/rainfall are frozen fast-step).  !
            bf%soil_wat_in  = 0.0_wp                            ; bf%soil_wat_out = 0.0_wp
            !----- condensation (dew) leaves the CAS as liquid at Tcas -> a whole-column water + liquid-   !
            !      enthalpy OUTPUT (the CAS-side loss is already in src_vap/src_enth, so cas_water/energy   !
            !      close automatically). -----------------------------------------------------------------!
            !----- ground_rad is the snowfac-BLENDED radiative input (= abs_sw+abs_lw when bare). ---!
            !----- #78 item 4: e_infil (pond -> soil) and e_clip (soil -> pond) are now transfers    !
            !      between two TRACKED stores, so they telescope and must NOT be boundary terms.     !
            !      The boundary rainfall input and the runoff output are added once at the outer level. !
            bf%whole_enth_in= surf_tend%coh_rnet + frozen%snow%ground_rad + e_floor
            !----- row 1b: surf_tend%cond's enthalpy is NO LONGER a boundary loss -- the condensate is        !
            !      deposited into soil layer 1 by the caller, carrying this same u_liq(t_cas1). ------!
            bf%whole_enth_out= g_atm_heat*(enth1 - frozen%cas%enthalpy_atm) + e_drain
            bf%whole_wat_in = 0.0_wp                            ; bf%whole_wat_out = g_atm_vapour*(shv1 - frozen%cas%shv_atm)
            bf%whole_cond   = surf_tend%cond                     ! row 1b: deposited into a store, not lost
            bf%whole_cond_enth = surf_tend%cond_enth   ! EXACTLY what surface_derivs debited from the CAS (one number, both sides)
            bf%atm_heat_out = g_atm_heat*cp_air*(t_cas1 - frozen%cas%mo_theta_atm)   ! the reported H, on the ledger's basis
            bf%atm_vap_out  = g_atm_vapour*(shv1  - frozen%cas%shv_atm)
      end if
   end subroutine column_be_stage

   !---------------------------------------------------------------------------------------!
   ! newton_surface_solve -- the ARROWHEAD: a direct 2x2 Newton solve of the coupled backward-Euler   !
   ! CAS surface block { R_H, R_q } = 0 for (cas_enthalpy H1, cas_shv q1), where the surface sources   !
   ! src_enth/src_vap depend nonlinearly on (H1,q1) through tcas, qcas and qsat. Replaces the leaf<->  !
   ! CAS Picard iteration: quadratic convergence (~1 step) + robust near saturation. Numerical         !
   ! Jacobian (finite-difference surface_derivs) -- captures the strong VPD self-limiting d src_vap/dq  !
   ! with no derivation risk. Singular-Jacobian guard + line search + supersaturation clamp + eval cap; !
   ! never error stops (GPU-safe). Commits the CAS via the FLUX form so budgets close for ANY surf_tend.      !
   !---------------------------------------------------------------------------------------!
   subroutine newton_surface_solve(y, cas, frozen, t_ground, n, dt, cas_mass_capacity, g_atm_heat, g_atm_vapour, enth1, shv1, &
                                   surf_tend, nfeval, ok)
      type(column_state_t),   intent(in)    :: y
      type(cas_boundary_t),   intent(in)    :: cas      !< the stage's CAS boundary (its own conductances)
      type(column_frozen_t),  intent(in)    :: frozen   !< tissue / film / ground / snow for surface_derivs
      real(wp),               intent(in)    :: t_ground
      real(wp) :: co2_unused
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt, cas_mass_capacity, g_atm_heat, g_atm_vapour
      real(wp),               intent(out)   :: enth1, shv1
      type(surface_tend_t),   intent(out)   :: surf_tend
      integer(ik),            intent(out)   :: nfeval
      logical,                intent(out)   :: ok

      type(surface_state_t) :: y_stage
      real(wp)    :: H0, q0, Hk, qk, R_H, R_q, J11, J12, J21, J22, detJ, delH, delq
      real(wp)    :: lam, rn0, Ht, qt, RHt, Rqt
      integer(ik) :: it, ls
      real(wp),    parameter :: RTOL_N = 1.0e-7_wp, ATOL_H = 5.0e1_wp, ATOL_Q = 1.0e-6_wp
      real(wp),    parameter :: DETEPS = 1.0e-30_wp
      integer(ik), parameter :: NEWT_MAX = 4_ik, LS_MAX = 6_ik, FEVAL_CAP = 24_ik

      H0 = y%cas_enthalpy ; q0 = y%cas_shv
      Hk = H0 ; qk = q0 ; nfeval = 0_ik ; ok = .false.
      y_stage%cas_co2 = y%cas_co2
      y_stage%cas_enthalpy = Hk ; y_stage%cas_shv = qk
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,               &
         t_ground, n, surf_tend) ; nfeval = nfeval + 1_ik
      R_H = cas_mass_capacity*(Hk - H0)/dt - surf_tend%src_enth - g_atm_heat*(cas%enthalpy_atm - Hk)
      R_q = cas_mass_capacity*(qk - q0)/dt - surf_tend%src_vap  - g_atm_vapour*(cas%shv_atm  - qk)

      do it = 1_ik, NEWT_MAX
         if ( abs(R_H)*dt/cas_mass_capacity <= ATOL_H + RTOL_N*abs(Hk) .and.                                  &
              abs(R_q)*dt/cas_mass_capacity <= ATOL_Q + RTOL_N*abs(qk) ) then
            ok = .true. ; exit
         end if
         call jac_surface(Hk, qk, y%cas_co2, cas, frozen, t_ground, surf_tend, n, cas_mass_capacity, g_atm_heat, g_atm_vapour, dt, &
                          J11, J12, J21, J22, nfeval)
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
            y_stage%cas_enthalpy = Ht ; y_stage%cas_shv = qt
            call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,         &
               t_ground, n, surf_tend) ; nfeval = nfeval + 1_ik
            RHt = cas_mass_capacity*(Ht - H0)/dt - surf_tend%src_enth - g_atm_heat*(cas%enthalpy_atm - Ht)
            Rqt = cas_mass_capacity*(qt - q0)/dt - surf_tend%src_vap  - g_atm_vapour*(cas%shv_atm  - qt)
            if (RHt*RHt + Rqt*Rqt <= (1.0_wp - 1.0e-4_wp*lam)*rn0) exit          ! Armijo
            lam = 0.5_wp*lam
         end do
         Hk = Ht ; qk = qt ; R_H = RHt ; R_q = Rqt
         if (nfeval >= FEVAL_CAP) exit
      end do

      !----- authoritative final eval + flux-form commit (conservation holds for ANY surf_tend). -----------!
      y_stage%cas_enthalpy = Hk ; y_stage%cas_shv = qk
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,               &
         t_ground, n, surf_tend) ; nfeval = nfeval + 1_ik
      call cas_box_commit(H0, q0, 0.0_wp, surf_tend, cas, cas_mass_capacity, 1.0_wp, g_atm_heat, g_atm_vapour, 0.0_wp, dt, enth1, &
                          shv1, co2_unused)
   end subroutine newton_surface_solve

   !---------------------------------------------------------------------------------------!
   ! The backward-Euler canopy-air box commit, routed through the SHARED kernel                      !
   ! meds_cas_biophysics%cas_column_step_implicit (which was exported but had no caller while this    !
   ! module re-implemented its three formulas inline). One implementation, both schemes' box.         !
   !---------------------------------------------------------------------------------------!
   pure subroutine cas_box_commit(h0, q0, c0, surf_tend, cas, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, &
                                  g_atm_co2, dt, h1, q1, c1)
      real(wp),               intent(in)  :: h0, q0, c0, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, &
           g_atm_co2, dt
      type(surface_tend_t),   intent(in)  :: surf_tend
      type(cas_boundary_t),   intent(in)  :: cas
      real(wp),               intent(out) :: h1, q1, c1
      type(cas_column_t) :: box
      type(cas_source_t) :: src
      box%air_mass_capacity        = cas_mass_capacity ; box%air_molar_capacity     = cas_molar_capacity
      box%atm_conductance_enthalpy = g_atm_heat  ; box%atm_conductance_vapor  = g_atm_vapour ; box%atm_conductance_co2 = g_atm_co2
      box%atm_enthalpy             = cas%enthalpy_atm ; box%atm_specific_humidity = cas%shv_atm ; box%atm_co2 = cas%co2_atm
      src%surface_enthalpy_source  = surf_tend%src_enth ; src%surface_vapor_source = surf_tend%src_vap
      src%biotic_co2_source        = cas%nee_biotic
      call cas_column_step_implicit(h0, q0, c0, src, box, dt, h1, q1, c1)
   end subroutine cas_box_commit

   !----- 2x2 numerical Jacobian of (R_H, R_q) w.r.t. (H, q) by forward-differencing surface_derivs. --!
   subroutine jac_surface(Hk, qk, co2, cas, frozen, t_ground, surf_tend, n, cas_mass_capacity, g_atm_heat, g_atm_vapour, dt, J11, &
                          J12, J21, J22, &
                          nfeval)
      real(wp),               intent(in)    :: Hk, qk, co2, cas_mass_capacity, g_atm_heat, g_atm_vapour, dt
      type(cas_boundary_t),   intent(in)    :: cas
      type(column_frozen_t),  intent(in)    :: frozen
      real(wp),               intent(in)    :: t_ground
      type(surface_tend_t),   intent(in)    :: surf_tend         ! base eval at (Hk,qk)
      integer(ik),            intent(in)    :: n
      real(wp),               intent(out)   :: J11, J12, J21, J22
      integer(ik),            intent(inout) :: nfeval
      type(surface_state_t) :: y_stage
      type(surface_tend_t)  :: sfp
      real(wp) :: dH, dq, dse_dH, dsv_dH, dse_dq, dsv_dq
      real(wp), parameter :: SQEPS = 1.4901161e-8_wp, HSCALE = 1.0e4_wp, QSCALE = 1.0e-3_wp
      dH = SQEPS * max(abs(Hk), HSCALE)
      dq = SQEPS * max(abs(qk), QSCALE)
      y_stage%cas_co2 = co2
      y_stage%cas_enthalpy = Hk + dH ; y_stage%cas_shv = qk
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,               &
         t_ground, n, sfp) ; nfeval = nfeval + 1_ik
      dse_dH = (sfp%src_enth - surf_tend%src_enth)/dH ; dsv_dH = (sfp%src_vap - surf_tend%src_vap)/dH
      y_stage%cas_enthalpy = Hk ; y_stage%cas_shv = qk + dq
      call surface_derivs(y_stage, cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,               &
         t_ground, n, sfp) ; nfeval = nfeval + 1_ik
      dse_dq = (sfp%src_enth - surf_tend%src_enth)/dq ; dsv_dq = (sfp%src_vap - surf_tend%src_vap)/dq
      J11 = cas_mass_capacity/dt + g_atm_heat - dse_dH ; J12 =              - dse_dq
      J21 =              - dsv_dH  ; J22 = cas_mass_capacity/dt + g_atm_vapour - dsv_dq
   end subroutine jac_surface

   !---------------------------------------------------------------------------------------!
   ! advance_water_mass_full -- operator-split plant water MASS over the FULL dt from y%*_water_mass,   !
   ! driven by the CALLER-SUPPLIED transp_c_bw -- the SAME per-cohort transpiration (already b-weighted  !
   ! across the stage(s) that make up this step, sec below) the CAS's own vapour balance used, so the     !
   ! mass debit and the CAS credit are consistent to within the tableau's own stage algebra (not a         !
   ! separate, later evaluation at a possibly-different point).                                             !
   !                                                                                          !
   ! THE UPDATE IS A CORRECTOR, not the bare Euler step it used to be. The reason is exact algebra, not     !
   ! a tolerance judgement. solve_plant_water defines its own reported flux FROM its own storage change,    !
   !     flux%sapflow = dw_l/dt + transp,   flux%root_uptake = (dw_l + dw_w)/dt + transp,                   !
   ! so `W + dt*(sapflow - transp)` reproduces the kernel's matrix-exponential endpoint EXACTLY -- to        !
   ! machine precision, at any dt -- PROVIDED the transp used to price the sapflow is the transp actually    !
   ! debited. The pre-pass priced sapflow_frozen at transp_pp (the state^n FULL demand); this routine        !
   ! debits transp_c_bw (the b-weighted REALISED demand). The gap between the two is therefore not a          !
   ! discretisation error of the mass ODE at all -- it is exactly dt*(transp_pp - transp_bw), a pure FLUX     !
   ! INCONSISTENCY that grows linearly in dt and never converges away within a step. Measured on the 3 h      !
   ! midday fixture it was the WHOLE of the psi_leaf dt-divergence: 1.45 MPa of error at dt_fast = 900 s.     !
   !                                                                                          !
   ! Re-solving the kernel here on the REALISED transpiration -- same Category-0 frozen conductances and      !
   ! soil potentials as the pre-pass, so this is one semi-discretisation and not a second linearisation      !
   ! point -- makes the credit and the debit the same number again. Measured effect at dt_fast = 900 s:       !
   ! psi_leaf error 1.454 -> 0.0046 MPa (314x), against a common dt = 5 s reference that BOTH versions        !
   ! agree on to 8e-4 MPa (i.e. this changes the discretisation, NOT the dt->0 limit). T_cas differs by        !
   ! 6e-5 K at dt = 25 s and 0.007 K at 900 s -- also vanishing with dt, as a consistent scheme must.          !
   ! (MEDS_ED2_RK45_DESIGN.md sec 1/4/5 described the pre-corrector Euler form.) --------------------------!
   subroutine advance_water_mass_full(y, frozen, n, nsl, dt, transp_c_bw, y_out)
      type(column_state_t),  intent(in)    :: y
      type(column_frozen_t), intent(in)    :: frozen
      integer(ik),           intent(in)    :: n, nsl
      real(wp),              intent(in)    :: dt
      real(wp),              intent(in)    :: transp_c_bw(n)  !< [kg/m2 ground/s] per-cohort, b-weighted over the step
      type(column_state_t),  intent(inout) :: y_out
      real(wp)    :: transp_i
      integer(ik) :: i
      !----- CORRECTOR locals (only touched on the re-solve path). --------------------------------!
      real(wp)    :: transp_pp(n), psi_c(N_HYDRO, n)
      real(wp)    :: sapflow_c(n), uptake_c(n), uptake_layer_c(nsl, n)
      real(wp)    :: psi_leaf_c(n), psi_wood_c(n), plc_c(n)
      integer(ik) :: nsub_c(n)
      logical     :: conv_c(n)
      real(wp)    :: sap_use(n), upt_use(n)

      !----- THE CORRECTOR (see the header). The pre-pass solved the kernel on transp_pp = the       !
      !      state^n FULL demand; the step actually debits transp_c_bw, the b-weighted REALISED       !
      !      demand. Since the kernel defines sapflow = dw_l/dt + transp, an Euler step on the mass   !
      !      reproduces the kernel's own dw EXACTLY when the two transpirations agree -- and differs   !
      !      by exactly dt*(transp_pp - transp_bw) when they do not. That residual is a pure FLUX      !
      !      INCONSISTENCY (credit priced at one transpiration, debit taken at another), it grows      !
      !      linearly in dt, and it is the whole of the psi_leaf dt-divergence. Re-solving here on     !
      !      the realised transpiration -- same FROZEN conductances and soil potentials, so this is    !
      !      still one Category-0 semi-discretisation, not a second linearisation point -- makes the   !
      !      credit and the debit the same number again.                                              !
      !                                                                                               !
      !      Guarded on allocation: the RK4/IMEX-Euler oracle hand-builds its own `frozen` and never      !
      !      populates the hydraulics boundary inputs, so it keeps the frozen fluxes (same rule as     !
      !      surface_frozen_t%mo_live). --------------------------------------------------------------!
      sap_use(1:n) = frozen%plant%sapflow_frozen(1:n)
      upt_use(1:n) = frozen%plant%uptake_frozen(1:n)
      if (allocated(frozen%roots%rhizo_cond) .and. allocated(frozen%roots%psi_soil_pre) .and. dt > tiny_num) then
         do i = 1_ik, n
            transp_pp(i) = transp_c_bw(i) / max(frozen%plant%nplant(i), tiny_num)
         end do
         psi_c(NODE_LEAF, 1:n) = psi_from_water_content(y%leaf_water_mass(1:n),                      &
              frozen%params%hydraulics_params%leaf_pi0, frozen%params%hydraulics_params%leaf_elastic_mod, &
                   frozen%params%hydraulics_params%leaf_apoplast_frac,    &
              frozen%params%hydraulics_params%leaf_water_sat, frozen%plant%bleaf(1:n))
         psi_c(NODE_WOOD, 1:n) = psi_from_water_content(y%wood_water_mass(1:n),                      &
              frozen%params%hydraulics_params%wood_pi0, frozen%params%hydraulics_params%wood_elastic_mod, &
                   frozen%params%hydraulics_params%wood_apoplast_frac,    &
              frozen%params%hydraulics_params%wood_water_sat, frozen%plant%bsap(1:n) + frozen%plant%broot(1:n))
         call solve_plant_water_batch(n, nsl, transp_pp(1:n), frozen%plant%bleaf(1:n), frozen%plant%bsap(1:n),         &
              frozen%plant%broot(1:n), frozen%plant%sap_area(1:n), frozen%plant%height(1:n),                  &
                 frozen%plant%leaf_area(1:n),                &
              frozen%roots%psi_soil_pre(1:nsl), frozen%params%soil%z_node(1:nsl), frozen%roots%rhizo_cond(1:nsl, 1:n),           &
              frozen%params%hydraulics_params, frozen%params%hydraulics_opts, dt, psi_c(:, 1:n), sapflow_c(1:n), uptake_c(1:n), &
              uptake_layer_c(1:nsl, 1:n), psi_leaf_c(1:n), psi_wood_c(1:n), plc_c(1:n),              &
              nsub_c(1:n), conv_c(1:n))
         !----- TAKE THE CORRECTOR'S SAPFLOW ONLY; the wood<->soil interface KEEPS uptake_frozen, which  !
         !      is the SAME number the soil column already committed as its root sink (frozen%roots%uptake,        !
         !      rescaled by the soil's own fwilt-limited supply). That keeps "one flux, both sides"       !
         !      across an interface spanning two OPERATOR-SPLIT solves, so the whole-column water         !
         !      ledger closes exactly.                                                                    !
         !                                                                                          !
         !      MEASURED CONSEQUENCE, and why this is the right asymmetry rather than a shortcut.         !
         !      uptake_frozen is priced at the pre-pass transpiration, so at long dt_fast the plant       !
         !      receives ~1.3% less water than it should: plant total 1.9139 vs 1.9392 kg/m2 at 900 s     !
         !      against a 5 s reference. That TOTAL error is fixed by the committed uptake and is         !
         !      INDEPENDENT of the root boundary condition -- a specified-flux (Neumann) root BC was      !
         !      built and measured, and it produced the IDENTICAL total (1.9139). All the BC chooses is   !
         !      WHERE the deficit is parked:                                                              !
         !        potential BC (this code): all of it into the wood -- psi_leaf 4.6e-3, psi_wood 1.74e-1  !
         !        specified-flux BC:        split by capacitance -- psi_leaf 1.61e-1, psi_wood 1.60e-1    !
         !      The leaf is the consequential node (it drives the stomatal water-stress feedback, hence   !
         !      transpiration and GPP), while psi_wood at these potentials is far from wood_psi50 and     !
         !      drives almost no PLC. Parking the error in the wood is therefore the better trade, and    !
         !      the Neumann variant was NOT landed. Closing the deficit itself requires re-pricing the    !
         !      SOIL side (the pre-pass would have to see the realised transpiration); no plant-side      !
         !      boundary condition can do it. See <scratchpad>/neumann_root_bc.patch.                     !
         sap_use(1:n) = sapflow_c(1:n)
      end if

      do i = 1_ik, n
         transp_i = transp_c_bw(i) / max(frozen%plant%nplant(i), tiny_num)
         !----- KNOWN DEFERRED EDGE CASE: unlike psi (whose PV-curve capacitance self-limits as        !
         !      tissue dries, dw/dpsi -> 0 in the flaccid tail), the mass ODE is a plain linear Euler    !
         !      step with no such restoring force -- sapflow_frozen/uptake_frozen are the STATE-n        !
         !      average, but transp_i here is the ENDPOINT-refreshed demand, so a large intra-step       !
         !      CAS swing could in principle debit more water than is actually in storage. Floor at a    !
         !      tiny positive value (not 0, so a subsequent psi_from_water_content diagnosis never        !
         !      divides by an exact 0 rwc) rather than let it go negative -- a bookkept boundary term      !
         !      for this clamp is deferred (mirrors the P1 surf_overflow precedent, not yet needed        !
         !      here: unobserved in this pass's test scenarios, see MEDS_ED2_RK45_DESIGN.md P2 notes). ---!
         y_out%leaf_water_mass(i) = max(y%leaf_water_mass(i) + dt*(sap_use(i) - transp_i), tiny_num)
         y_out%wood_water_mass(i) = max(y%wood_water_mass(i)                                          &
                                   + dt*(upt_use(i) - sap_use(i)), tiny_num)
      end do
   end subroutine advance_water_mass_full

   !---------------------------------------------------------------------------------------!
   ! advance_surf_water_full -- operator-split canopy-SURFACE water (sec 3.4, P2c) over the FULL dt   !
   ! from y%*_surf_water, mirroring advance_water_mass_full exactly: frozen%film%intercept_leaf/wood is the     !
   ! FROZEN Act-1 inflow, film_evap_*_bw is the b-weighted per-stage outflow (the SAME number the        !
   ! caller's sf2/sf3 credited to the CAS, sec 1/3/4/5's "one flux, both sides" discipline). Floors at    !
   ! exactly 0.0 (not tiny_num like the mass version -- surf_water feeds no psi/rwc division downstream).  !
   ! CAPACITY clamping + overflow bookkeeping (mirroring the split path's post-hoc treatment, sec 9)         !
   ! happens LATER, in column_fast_step_ark, not here -- this is the raw, unclamped Euler step only. --------!
   pure subroutine advance_surf_water_full(y, frozen, n, dt, film_evap_leaf_bw, film_evap_wood_bw, y_out)
      type(column_state_t),  intent(in)    :: y
      type(column_frozen_t), intent(in)    :: frozen
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt
      real(wp),               intent(in)    :: film_evap_leaf_bw(n), film_evap_wood_bw(n)
      type(column_state_t),  intent(inout) :: y_out
      integer(ik) :: i
      !----- NOT floored here (unlike an earlier version of this routine): the caller's post-hoc         !
      !      overflow/deficit bookkeeping (column_fast_step_ark, mirroring column_fast_step_rk45's own    !
      !      treatment) needs to SEE a possibly-negative value to correctly credit the shortfall back to    !
      !      the whole-column ledger's outflow -- clamping prematurely here would silently discard it,        !
      !      fabricating mass exactly like an unclamped store corrupts the NEXT step's frozen interception.  !
      do i = 1_ik, n
         y_out%leaf_surf_water(i) = y%leaf_surf_water(i) + dt*(frozen%film%intercept_leaf(i) - film_evap_leaf_bw(i))
         y_out%wood_surf_water(i) = y%wood_surf_water(i) + dt*(frozen%film%intercept_wood(i) - film_evap_wood_bw(i))
      end do
   end subroutine advance_surf_water_full

end module meds_fast_be_stage
