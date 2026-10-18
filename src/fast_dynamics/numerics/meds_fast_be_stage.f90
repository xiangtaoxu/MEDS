! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_be_stage -- the IMPLICIT STAGE of the fast-loop march, and the water-mass      !
! advance that rides with it.                                                              !
!                                                                                          !
! One ESDIRK stage (docs/science/numerical_scheme.md §3 `ark`, eq. (1)) is a               !
! backward-Euler solve of the coupled surface system: the canopy air's enthalpy and        !
! humidity against the leaf, wood and ground temperatures, closed by a small Newton        !
! iteration over a forward-difference Jacobian, with CO2 riding the same backward-Euler    !
! box. That machinery is independent of WHICH tableau calls it, which is the point of      !
! pulling it out: `meds_fast_ark` keeps the tableau and the march, and the RHS test        !
! (test_column_derivs) reaches the same stage solver directly.                             !
!                                                                                          !
! Contents: `column_be_stage` (the stage), `newton_surface_solve` + `jac_surface` +        !
! `cas_box_commit` (its inner solve), and `advance_plant_water` /                          !
! `advance_surf_water_full` (the mass-form water advance a stage commits alongside the     !
! energy solve).                                                                           !
!==========================================================================================!
module meds_fast_be_stage
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, cp_air, rho_h2o
   use meds_config, only : INTEG_ARK, CTRL_L2_STRICT
   use meds_fast_control, only : step_control_factor
   use meds_soil_types, only : energy_forcing_t, energy_flux_t, snow_melt_t
   use meds_patch_state_vars,   only : soil_energy_column_t
   use meds_fast_time_derivs, only : surface_derivs, cas_conductances
   use meds_cas_biophysics, only : cas_column_step_implicit, cas_column_t, cas_source_t
   use meds_column_state_ops, only : state_init, state_sub, bflux_zero, bflux_add, bflux_bweight, clamp_cas, &
                                     soil_water_store, soil_energy_store, plant_water_store,                    &
                                     canopy_film_store, deposit_condensate, clamp_canopy_film, unpack_column_state, &
                                     diagnose_soil_temps, assemble_soil_energy_forcing, apply_process_mask
   use meds_fast_types, only : column_config_t, column_budget_t, alloc_column_cohort, column_state_t, &
                               column_frozen_t, surface_state_t, cas_boundary_t, surface_tend_t, stage_bflux_t, &
                               column_bflux_t, error_control_t, column_tend_t, &
                               mask_is_full
   use meds_soil_energy, only : soil_energy_step_implicit
   use meds_cas_biophysics, only : cas_column_t, cas_source_t, cas_column_step_implicit
   use meds_ground_biophysics, only : snow_accumulate, snow_drain_meltwater, snow_cover_fraction, ground_surface_fluxes
   use meds_therm_lib, only : cas_temp_of_enthalpy, internal_energy_liquid, internal_energy_to_temp, internal_energy_ice, &
                              temp_of_liquid_enthalpy
   use meds_budget_check, only : budget_check, budget_energy_rate_floor, budget_water_rate_floor, budget_co2_rate_floor
   implicit none
   private

   public :: column_be_stage, advance_plant_water, advance_surf_water_full

contains

   !---------------------------------------------------------------------------------------!
   ! column_be_stage -- ONE backward-Euler stage of the STIFF block: the CAS twins (implicit in the     !
   ! atmospheric exchange; np > 1 -> the coupled 2x2 Newton) and soil heat (implicit, Thomas solve),    !
   ! driven by the frozen surface sources. Water is the soil-plant water network's: soil water moves    !
   ! at its steady rate and plant water is passed through, then advanced once over the step             !
   ! (advance_plant_water). It solves Y = base + dt*f(Y), so each ESDIRK2 stage (dt = gamma*h) is       !
   ! one call.                                                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine column_be_stage(y, frozen, col_config, n, nsl, dt, y_out, surf_tend, niter, bf)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      type(column_config_t), intent(in)  :: col_config   !< the column's parameters (soil, thermal, hydraulics)
      integer(ik),           intent(in)  :: n, nsl
      real(wp),              intent(in)  :: dt
      !----- The stage's state and its OWN surface tendencies (incl. transp_c), for the caller to        !
      !      b-weight into the plant water-mass update -- the SAME surf_tend the stage's own bf/CAS      !
      !      source used, so the leaves' water debit and the canopy air's vapour credit agree            !
      !      (docs/science/numerical_scheme.md §5a). Both are filled in place (the caller's              !
      !      workspace keeps their arrays). -------------------------------------------------------------!
      type(column_state_t),  intent(inout) :: y_out
      type(surface_tend_t),  intent(inout) :: surf_tend
      integer(ik), optional, intent(in)  :: niter    !< 1 = uncoupled BE baseline; >1 = coupled leaf<->CAS Newton
      type(stage_bflux_t), optional, intent(out) :: bf  !< per-stage boundary-flux RATES for the ARK ledger

      type(surface_state_t)      :: y_stage
      type(cas_boundary_t)       :: cas_stage   !< the CAS boundary with THIS stage's conductances
      type(soil_energy_column_t) :: se
      type(energy_forcing_t)     :: eforc
      type(energy_flux_t)        :: eflux
      real(wp)    :: t_ground, fliq1, wmass1, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2
      real(wp)    :: enth1, shv1, e_infil, e_drain, t_cas1, water_heat(nsl)
      real(wp)    :: co21, enth_unused, shv_unused
      integer(ik) :: k, np, nfeval
      logical     :: ok

      np = 1_ik ; if (present(niter)) np = max(1_ik, niter)

      call state_init(y, n, nsl, y_out)

      !----- diagnose the soil-top temperature so the ground skin sees the current state. --------!
      wmass1 = y%theta(1) * rho_h2o
      call internal_energy_to_temp(y%soil_energy(1), wmass1, col_config%soil_thermal%soil_dry_heat_capacity(1), t_ground, fliq1)

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

      !----- soil-heat column: implicit BE-Thomas (soil_energy_step_implicit). ---------------------------!
      se%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      !----- The water the soil-plant network moved over the step carries its enthalpy: across the      !
      !      faces at the network's time-mean flows, in from the pond at the pond's temperature, out with !
      !      root uptake and the top layer's exfiltration at each layer's start-of-step temperature (the  !
      !      numbers the wood and the pond were credited), and out of the bottom with drainage. The       !
      !      leaf pays the full vapour enthalpy of what it transpires (surface_derivs). The stages move    !
      !      theta at the network's steady rate (below), so the water the faces carry and the water the   !
      !      layers hold stay in step. --------------------------------------------------------------!
      water_heat(1:nsl) = frozen%roots%root_heat(1:nsl)
      water_heat(1)     = water_heat(1) + frozen%hydrology%exfil_enth
      e_infil = frozen%hydrology%infiltration * internal_energy_liquid(frozen%hydrology%t_infil)
      e_drain = frozen%hydrology%drainage     * internal_energy_liquid(frozen%hydrology%t_bot)
      call assemble_soil_energy_forcing(eforc, nsl, surf_tend%g_top, frozen%hydrology%geothermal, y%theta,          &
                                        water_heat, frozen%hydrology%w_flux_frozen,                                &
                                        frozen%hydrology%infiltration, frozen%hydrology%t_infil, e_drain)
      call soil_energy_step_implicit(se, eforc, col_config%soil_thermal, col_config%soil, col_config%energy, dt, eflux)
      y_out%soil_energy(1:nsl) = se%soil_energy(1:nsl)

      !----- Soil water is the water network's: the end-of-step theta is committed once, from the          !
      !      network (frozen%hydrology%theta1), in column_fast_step_ark -- one authority for soil          !
      !      water (docs/science/numerical_scheme.md §5b, "Every soil-water flow comes from one            !
      !      solve").                                                                                      !
      !                                                                                                    !
      !      The stages move theta at the network's steady rate rather than holding it at theta^n. The     !
      !      energy step above carries the enthalpy of the water the faces move, at the same steady        !
      !      rates; with theta held back, the next stage would read that enthalpy as heat in a too-dry     !
      !      layer: in heavy rain the top layer reads tens of kelvin too warm late in a step (up to        !
      !      ~50 K at BCI), and the next stage's ground skin sees it. A steady rate is integrated          !
      !      exactly by both the step and its error estimate, so theta adds nothing to the error, and      !
      !      the step ends at theta1.                                                                      !
      y_out%theta(1:nsl) = y%theta(1:nsl) + dt * frozen%hydrology%theta_rate(1:nsl)

      !----- pond PASSED THROUGH at state^n (state_init above already copied it): the network's        !
      !      end-of-step pond is committed once per dt_fast in column_fast_step_ark, and the pond      !
      !      is outside the error norm. ---------------------------------------------------------------!
      y_out%w_surface      = y%w_surface
      y_out%w_surface_enth = y%w_surface_enth

      !----- plant water MASS PASSED THROUGH (advanced by advance_plant_water, not here). ----!
      y_out%leaf_water_mass(1:n) = y%leaf_water_mass(1:n)
      y_out%wood_water_mass(1:n) = y%wood_water_mass(1:n)

      !----- emit this stage's boundary-flux RATES for the conservation ledger. The b-weight           !
      !      + cross-substep accumulation happens in ark2_column_step / adaptive_ark_march. Every      !
      !      quantity is the committed-state flux, so the accumulated amounts telescope to closure.    !
      if (present(bf)) then
         t_cas1 = cas_temp_of_enthalpy(enth1, shv1)         ! committed CAS temperature, for the reported H
            bf%cas_enth_in  = surf_tend%src_enth + g_atm_heat*frozen%cas%enthalpy_atm    ; bf%cas_enth_out = g_atm_heat*enth1
            bf%cas_vap_in   = surf_tend%src_vap  + g_atm_vapour*frozen%cas%shv_atm     ; bf%cas_vap_out  = g_atm_vapour*shv1
            bf%cas_co2_in   = frozen%cas%nee_biotic + g_atm_co2*frozen%cas%co2_atm  ; bf%cas_co2_out  = g_atm_co2*y_out%cas_co2
            !----- The bottom face is the flux the soil solve committed: zero on the adiabatic base,      !
            !      the conduction to deep_temp under the Dirichlet anchor. It crosses the boundary of     !
            !      both books, the soil's and the whole column's. ----------------------------------------!
            bf%soil_enth_in = surf_tend%g_top - eflux%bottom_heat + e_infil
            bf%soil_enth_out= sum(water_heat(1:nsl)) + e_drain
            !----- Soil water's storage change and its q_top/drainage/uptake flows are the water           !
            !      network's, booked once per step in column_fast_step_ark, so the per-stage soil-water    !
            !      amounts are zero and the whole-column water carries only the CAS-vapour exchange        !
            !      (drainage, runoff and rainfall are fixed for the step). --------------------------------!
            bf%soil_wat_in  = 0.0_wp                            ; bf%soil_wat_out = 0.0_wp
            !----- ground_rad is the snowfac-BLENDED radiative input (= abs_sw+abs_lw when bare). ---!
            !----- e_infil (pond -> soil) and the exfiltration (soil -> pond) move water between two       !
            !      TRACKED stores, so they telescope and are not boundary terms; the root uptake's         !
            !      enthalpy moves soil -> tissue, also tracked. The rainfall input and the runoff          !
            !      output are added once at the outer level (docs/science/numerical_scheme.md §5b,         !
            !      "The pond holds heat"). ----------------------------------------------------------------!
            bf%whole_enth_in= surf_tend%coh_rnet + frozen%snow%ground_rad - eflux%bottom_heat
            !----- The condensate is not a boundary loss: the caller deposits it into soil layer 1         !
            !      with the enthalpy surface_derivs debited from the CAS (whole_cond,                      !
            !      whole_cond_enth); the CAS-side loss is already in src_vap/src_enth, so the CAS          !
            !      water and energy budgets close as they are. --------------------------------------------!
            bf%whole_enth_out= g_atm_heat*(enth1 - frozen%cas%enthalpy_atm) + e_drain
            bf%whole_wat_in = 0.0_wp                            ; bf%whole_wat_out = g_atm_vapour*(shv1 - frozen%cas%shv_atm)
            bf%whole_cond   = surf_tend%cond                     ! deposited into soil layer 1, not lost
            bf%whole_cond_enth = surf_tend%cond_enth   ! EXACTLY what surface_derivs debited from the CAS (one number, both sides)
            !----- The reported H, against theta_atm (T at zref), the ledger's basis. ----------------------!
            bf%atm_heat_out = g_atm_heat*cp_air*(t_cas1 - frozen%cas%mo_theta_atm)
            bf%atm_vap_out  = g_atm_vapour*(shv1  - frozen%cas%shv_atm)
      end if
   end subroutine column_be_stage

   !---------------------------------------------------------------------------------------!
   ! newton_surface_solve -- the ARROWHEAD: a direct 2x2 Newton solve of the coupled backward-Euler     !
   ! CAS surface block { R_H, R_q } = 0 for (cas_enthalpy H1, cas_shv q1), where the surface sources    !
   ! src_enth/src_vap depend nonlinearly on (H1,q1) through tcas, qcas and qsat. Converges              !
   ! quadratically (~1 step) and stays robust near saturation, where a leaf<->CAS Picard iteration      !
   ! does not. Numerical Jacobian (forward-difference surface_derivs, jac_surface) -- captures the      !
   ! strong VPD self-limiting d src_vap/dq with no derivation risk. Singular-Jacobian guard + line      !
   ! search + eval cap, and no supersaturation clamp (below); never error stops (GPU-safe). Commits     !
   ! the CAS via the FLUX form so budgets close for ANY surf_tend.                                      !
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
      type(surface_tend_t),   intent(inout) :: surf_tend   !< filled in place by every evaluation
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
   ! The backward-Euler canopy-air box commit, through the SHARED box kernel                           !
   ! meds_cas_biophysics%cas_column_step_implicit (docs/science/canopy_air_space_biophysics.md         !
   ! eq. (2)): one implementation of the box, not a copy of its formulas here.                         !
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
   ! advance_plant_water -- each cohort's leaf and wood water over an ARK step of length dt, at the   !
   ! soil-plant water network's time-mean rates: sapflow and root uptake in, sapflow and the          !
   ! transpiration the stages charged the canopy air (transp_c_bw, b-weighted over the step) out. The !
   ! stages hold each cohort's transpiration at the network's (surface_derivs), so a step of dt_fast  !
   ! ends exactly where the network did, and the canopy air received exactly the water the leaves     !
   ! lost. Nothing is floored: the network keeps every store inside its curve.                       !
   !---------------------------------------------------------------------------------------!
   pure subroutine advance_plant_water(y, frozen, n, dt, transp_c_bw, y_out)
      type(column_state_t),  intent(in)    :: y
      type(column_frozen_t), intent(in)    :: frozen
      integer(ik),           intent(in)    :: n
      real(wp),              intent(in)    :: dt
      real(wp),              intent(in)    :: transp_c_bw(n)  !< [kg/m2 ground/s] per cohort, b-weighted over the step
      type(column_state_t),  intent(inout) :: y_out
      integer(ik) :: i
      do i = 1_ik, n
         y_out%leaf_water_mass(i) = y%leaf_water_mass(i)                                                  &
              + dt * (frozen%plant%sapflow_frozen(i) - transp_c_bw(i) / max(frozen%plant%nplant(i), tiny_num))
         y_out%wood_water_mass(i) = y%wood_water_mass(i)                                                  &
              + dt * (frozen%plant%uptake_frozen(i) - frozen%plant%sapflow_frozen(i))
      end do
   end subroutine advance_plant_water

   !---------------------------------------------------------------------------------------!
   ! advance_surf_water_full -- operator-split canopy-SURFACE water over the FULL dt from               !
   ! y%*_surf_water, the film analogue of advance_plant_water: frozen%film%intercept_leaf/wood is       !
   ! the frozen interception inflow, film_evap_*_bw the b-weighted per-stage outflow (the SAME          !
   ! number the caller's sf2/sf3 credited to the CAS). The raw, unclamped Euler step only:              !
   ! capacity clamping and the overflow/deficit bookkeeping happen in column_fast_step_ark. ------------!
   pure subroutine advance_surf_water_full(y, frozen, n, dt, film_evap_leaf_bw, film_evap_wood_bw, y_out)
      type(column_state_t),  intent(in)    :: y
      type(column_frozen_t), intent(in)    :: frozen
      integer(ik),            intent(in)    :: n
      real(wp),               intent(in)    :: dt
      real(wp),               intent(in)    :: film_evap_leaf_bw(n), film_evap_wood_bw(n)
      type(column_state_t),  intent(inout) :: y_out
      integer(ik) :: i
      !----- Not floored here: the caller's overflow/deficit bookkeeping (column_fast_step_ark, as       !
      !      column_fast_step_rk45 does for its own films) must SEE a negative value to credit the       !
      !      shortfall back to the whole-column outflow; clamping here would silently create that        !
      !      water. -------------------------------------------------------------------------------------!
      do i = 1_ik, n
         y_out%leaf_surf_water(i) = y%leaf_surf_water(i) + dt*(frozen%film%intercept_leaf(i) - film_evap_leaf_bw(i))
         y_out%wood_surf_water(i) = y%wood_surf_water(i) + dt*(frozen%film%intercept_wood(i) - film_evap_wood_bw(i))
      end do
   end subroutine advance_surf_water_full

end module meds_fast_be_stage
