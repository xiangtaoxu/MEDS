! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_time_derivs -- the pure whole-column tendency (RHS) the fast-loop integrators share      !
! (docs/science/numerical_scheme.md §3): a side-effect-free f(y) = dy/dt for the whole column, with  !
! no flux hidden inside a backward-Euler denominator and no store mutated, so any tableau can        !
! combine its evaluations.                                                                           !
!                                                                                                    !
! TWO ENTRY POINTS:                                                                                  !
!   * surface_derivs -- the CAS surface block (leaf/wood energy balance + ground skin + the three    !
!     CAS twins).                                                                                    !
!   * column_derivs  -- the WHOLE column: surface_derivs + the soil-heat column (the tendency-       !
!     exposing soil_energy_time_deriv, a side-effect-free sibling of its BE kernel) + soil and plant !
!     water, which move at the soil-plant water network's steady rates over the step.                !
!                                                                                                    !
! Faithfulness: the soil-heat tendency is the EXPLICIT RHS whose backward-Euler advance reproduces   !
! its BE kernel as dt -> 0, and the CAS twins' exactly (implicit in the atmosphere). Verified per-   !
! reservoir in test_column_derivs. The water network's flows and the leaf gas-exchange pre-pass are  !
! the frozen part of the step (Category 0, docs/science/numerical_scheme.md §2), held constant       !
! across the dt_fast.                                                                                !
!==========================================================================================!
module meds_fast_time_derivs
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, rho_h2o, n_soil_layer_max
   use meds_therm_lib, only : cas_molar_density, cas_enthalpy_to_temp, sat_specific_humidity,             &
                              sat_specific_humidity_temp_deriv, enthalpy_vapor, internal_energy_to_temp,  &
                              internal_energy_liquid
   use meds_soil_types, only : energy_forcing_t
   use meds_patch_state_vars,   only : soil_energy_column_t
   use meds_soil_energy, only : soil_energy_time_deriv
   use meds_cas_biophysics, only : cas_column_t, cas_source_t, cas_column_time_deriv
   use meds_ground_biophysics, only : ground_surface_fluxes
   use meds_canopy_aerodynamics, only : mo_surface_layer, cas_atm_conductances
   use meds_plant_biophysics, only : veg_energy_balance, lw_emission_slope
   use meds_column_state_ops, only : assemble_soil_energy_forcing
   use meds_fast_types, only : column_config_t, snow_stage_t, surface_state_t, surface_tend_t, cas_boundary_t, &
                               tissue_coefficients_t, canopy_film_capacity_t, ground_boundary_t, column_state_t, &
                               column_frozen_t, column_tend_t, &
                               stage_bflux_t, column_bflux_t
   implicit none
   private

   public :: surface_derivs, column_derivs, cas_conductances
   !----- public so any other user of the condensation timescale reads this one rather than     !
   !      keeping a copy that could drift out of step with it. ---------------------------------!
   public :: TAU_COND

   !----- CAS supersaturation relaxation timescale [s] for the smooth condensation sink: q relaxes to    !
   !      qsat with this e-folding time (dq/dt|_cond = -(q-qsat)/TAU_COND). Short vs dt_fast=1800s so a    !
   !      supersaturated canopy air condenses within the step, but SMOOTH (no clamp) so the ARK2 embedded  !
   !      error stays small. Physical dew/fog forms fast; the exact value is not critical (any strong     !
   !      restoring rate pins q near qsat). ----------------------------------------------------------!
   real(wp), parameter :: TAU_COND = 300.0_wp

contains

   !---------------------------------------------------------------------------------------!
   ! refresh_cas_conductances -- re-solve ONLY the bulk Monin-Obukhov surface layer at a live         !
   ! canopy-air state and return the three CAS<->atmosphere conductances                              !
   ! (docs/science/canopy_aerodynamics.md eq. (5)). Deliberately NOT a full canopy_aerodynamics       !
   ! call: the per-cohort boundary layers (leaf_gbw/leaf_gbh, and therefore h_coeff_leaf/             !
   ! g_transp_leaf) stay frozen, because each carries <= 0.7% of the lag gain while they are what     !
   ! makes a full refresh expensive (docs/science/numerical_scheme.md §2). Canopy geometry            !
   ! (displacement, roughness) and the reference-level state are frozen inputs -- only (T_cas, q_cas) !
   ! is live, which is exactly the loop that oscillates when frozen.                                  !
   !                                                                                          !
   ! temp2 == temp1 in canopy_aerodynamics (z0q = z0h), so vapour reuses the heat transfer factor;      !
   ! the CO2 conductance rides the same factor on the MOLAR capacity, computed at the live humidity.    !
   !---------------------------------------------------------------------------------------!
   pure subroutine refresh_cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      type(cas_boundary_t),   intent(in)  :: cas
      real(wp),               intent(in)  :: cas_enthalpy, cas_shv
      real(wp),               intent(out) :: g_atm_heat, g_atm_vapour, g_atm_co2
      real(wp) :: tcas, ustar, temp1, zeta, rib, obu, can_dmol
      tcas = cas_enthalpy_to_temp(cas_enthalpy, cas_shv)
      !----- column_prepass sets aenv%can_theta = tcas, so passing tcas here reproduces exactly    !
      !      what the frozen pre-pass fed this same kernel -- the ONLY difference is the state. ---!
      call mo_surface_layer(cas%aero_cfg, cas%mo_u_ref, cas%mo_zref, cas%mo_displace, cas%mo_rough, &
                            cas%mo_theta_atm, cas%mo_shv_atm, tcas, cas_shv,                      &
                            ustar, temp1, zeta, rib, obu)
      can_dmol = cas_molar_density(cas%mo_rho, cas_shv)
      call cas_atm_conductances(cas%mo_rho, can_dmol, ustar, temp1, temp1, g_atm_heat, g_atm_vapour, g_atm_co2)
   end subroutine refresh_cas_conductances

   !---------------------------------------------------------------------------------------!
   ! cas_conductances -- the ONE place that decides whether a surface-layer re-solve is needed.       !
   ! Returns the CAS<->atmosphere conductances for the canopy-air state (cas_enthalpy, cas_shv): a     !
   ! fresh solve when the bundle carries live MO inputs, otherwise the bundle's own values verbatim.   !
   ! Single-sourcing the rule is what keeps a tendency and the boundary flux charged against it from   !
   ! ever disagreeing -- see stage_bnd (meds_fast_rk45) for the failure this prevents.                 !
   !---------------------------------------------------------------------------------------!
   pure subroutine cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      type(cas_boundary_t),   intent(in)  :: cas
      real(wp),               intent(in)  :: cas_enthalpy, cas_shv
      real(wp),               intent(out) :: g_atm_heat, g_atm_vapour, g_atm_co2
      if (cas%mo_live) then
         call refresh_cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      else
         g_atm_heat = cas%g_atm_heat ; g_atm_vapour = cas%g_atm_vapour ; g_atm_co2 = cas%g_atm_co2
      end if
   end subroutine cas_conductances



   !---------------------------------------------------------------------------------------!
   ! surface_derivs -- the CAS surface-block RHS (leaf/wood energy balance + ground skin + the three !
   ! CAS twins), side-effect-free. Every tissue flux is linearized about the current tcas, which is  !
   ! also the longwave emission base (docs/science/vegetation_energy_dynamics.md eqs. (1)-(2)); the  !
   ! twins' tendencies are docs/science/canopy_air_space_biophysics.md eq. (1).                      !
   !---------------------------------------------------------------------------------------!
   pure subroutine surface_derivs(y, cas, tissue, film, ground, snow, t_ground, n, f)
      type(surface_state_t),        intent(in)  :: y
      type(cas_boundary_t),         intent(in)  :: cas      !< CAS <-> atmosphere boundary (a stage may pass its own conductances)
      type(tissue_coefficients_t),  intent(in)  :: tissue   !< per-cohort leaf/wood energy coefficients
      type(canopy_film_capacity_t), intent(in)  :: film     !< interception film
      type(ground_boundary_t),      intent(in)  :: ground   !< bare-ground boundary
      type(snow_stage_t),           intent(in)  :: snow     !< the snow stage's outcome
      real(wp),               intent(in)  :: t_ground   !< [K] soil-top temperature at THIS evaluation (a live input, not frozen)
      integer(ik),            intent(in)  :: n
      !----- Filled in place: every field is assigned below, and the per-cohort arrays are allocated  !
      !      only when missing or of another length, so a Newton solve that evaluates this many times !
      !      on one record allocates them once. ------------------------------------------------------!
      type(surface_tend_t),   intent(inout) :: f

      real(wp)    :: tcas, qcas, qsat_c, dqdt, esat
      real(wp)    :: lw_slope, le_slope, le_ref, dtl, tl, transp_i, dh, drnet
      real(wp)    :: lw_slope_w, dtw, tw, transp_w, dh_w, drnet_w    !< WOOD balance (own store)
      real(wp)    :: coh_h, coh_qw, coh_transp, coh_rnet, coh_film_evap
      real(wp)    :: h_evap_l, h_film_l, h_evap_w, h_film_w   !< [J/kg] energy per kg evaporated (leaf/wood; transp/film)
      real(wp)    :: h_bare, le_soil   !< bare-soil half of the snowfac blend
      !----- Canopy-SURFACE water: the wetted-fraction film-evap latent terms, using the FROZEN     !
      !      conductance (film%g_film_leaf/w, from leaf_film_coeff in the pre-pass) but the state-  !
      !      dependent dqdt/qsat_c-qcas -- mirroring le_slope/le_ref's own frozen-conductance/      !
      !      live-state split for the dry pathway. Harmless when canopy_water_on is off:            !
      !      film%f_wet_c(i) stays 0.0 there, which makes veg_energy_balance's wet pathway vanish   !
      !      identically, so these are pure no-ops. ------------------------------------------------!
      real(wp)    :: le_slope_wet, le_ref_wet, le_slope_wet_w, le_ref_wet_w
      real(wp)    :: gah_l, gaw_l, gac_l      !< the conductances this evaluation uses
      type(cas_source_t) :: cas_src
      type(cas_column_t) :: cas_col
      integer(ik) :: i

      call ensure_surface_tend(f, n)

      tcas   = cas_enthalpy_to_temp(y%cas_enthalpy, y%cas_shv)
      qcas   = y%cas_shv
      qsat_c = sat_specific_humidity(tcas, cas%press)
      dqdt   = sat_specific_humidity_temp_deriv(tcas, cas%press)

      coh_h = 0.0_wp ; coh_qw = 0.0_wp ; coh_transp = 0.0_wp ; coh_rnet = 0.0_wp
      coh_film_evap = 0.0_wp
      do i = 1_ik, n
         lw_slope = lw_emission_slope(tissue%leaf_emiss(i), tcas, tissue%lai(i))
         !----- The leaf pays the FULL specific enthalpy of the vapour it sheds, h_evap = enthalpy_vapor !
         !      at the canopy-air temperature the balance is linearized around (the same reference       !
         !      qsat_c/dqdt use; the cp_vap*(t_leaf - t_cas) difference is ~0.2% of h_evap), and the     !
         !      CAS receives the SAME number -- one flux, both sides. Charging the leaf latent_heat_vap  !
         !      only, with the liquid part (h_evap - L) charged to the soil, would charge the soil twice !
         !      beside the explicit soil -> wood -> leaf advected-enthalpy chain (qloss, qwflux_wl) and  !
         !      heat the leaf for free (~30 W/m2 at 3 mm/day). ED2's rk4_derivs charges the leaf         !
         !      tq2enthalpy(T_leaf) and the soil uint_water once, as this does. The FILM pays h_evap     !
         !      minus the liquid enthalpy the intercepted water arrived with (film_liquid_enthalpy), so  !
         !      film store + leaf + CAS close with no slack term                                         !
         !      (docs/science/vegetation_energy_dynamics.md eq. (2)). -----------------------------------!
         h_evap_l = enthalpy_vapor(tcas)
         h_film_l = h_evap_l - film%film_liquid_enthalpy
         if (tissue%transp_held) then
            !----- The water network has fixed this cohort's transpiration for the step: the leaf    !
            !      balance carries it as a fixed latent term on the dry share of the leaf, so the     !
            !      transpiration it reports is exactly the water the leaves lost. -------------------!
            le_slope = 0.0_wp
            le_ref   = h_evap_l * tissue%transp_net(i) / max(1.0_wp - film%f_wet_c(i), tiny_num)
         else
            le_slope = h_evap_l * cas%rho * tissue%g_transp_leaf(i) * dqdt
            le_ref   = h_evap_l * cas%rho * tissue%g_transp_leaf(i) * (qsat_c - qcas)
         end if
         le_slope_wet = h_film_l * cas%rho * film%g_film_leaf(i) * dqdt
         le_ref_wet   = h_film_l * cas%rho * film%g_film_leaf(i) * (qsat_c - qcas)
         !----- WOOD = the le_slope = le_ref = 0 case of the same kernel (no transpiration), solved     !
         !      first because a leaf with no area takes its wood's temperature (a bud sits on its      !
         !      twig), and wood with no area takes the canopy air's                                    !
         !      (docs/science/vegetation_energy_dynamics.md, "A tissue with no area"). q_wood_net      !
         !      (qloss - qwflux_wl) folds in via q_extra the same way qwflux_wl does for the leaf      !
         !      below (kept out of drnet). ------------------------------------------------------------!
         lw_slope_w = lw_emission_slope(tissue%wood_emiss(i), tcas, tissue%wai(i))
         h_evap_w = h_evap_l
         h_film_w = h_evap_w - film%film_liquid_enthalpy
         le_slope_wet_w = h_film_w * cas%rho * film%g_film_w(i) * dqdt
         le_ref_wet_w   = h_film_w * cas%rho * film%g_film_w(i) * (qsat_c - qcas)
         call veg_energy_balance(tissue%abs_sw_wood(i), tissue%abs_lw_wood(i), tissue%h_coeff_w(i),          &
                                    0.0_wp, lw_slope_w, 0.0_wp, tcas, tcas, tissue%wood_hcap_per_dt(i),        &
                                    tissue%t_wood0(i),                                                &
                                    dtw, tw, transp_w, dh_w, drnet_w, q_extra=tissue%q_wood_net(i),   &
                                    f_wet=film%f_wet_c(i), le_slope_wet=le_slope_wet_w,             &
                                    le_ref_wet=le_ref_wet_w, film_evap=f%film_evap_wood(i),        &
                                    h_evap=h_evap_w, h_evap_wet=h_film_w)
         f%wood_temp(i) = tw
         !----- LEAF: linearized about t_cas, relaxing with its own heat store (leaf_hcap_per_dt,        !
         !      t_leaf0; docs/science/vegetation_energy_dynamics.md eqs. (3)-(6)). qwflux_wl (sapflow's  !
         !      advected enthalpy) folds in via q_extra -- it shifts the equilibrium temperature (and    !
         !      hence dh/transp) like any other absorbed energy, but is kept OUT of drnet (an internal   !
         !      soil<->leaf transfer, not a boundary radiative input; see veg_energy_balance's own       !
         !      doc-comment). ---------------------------------------------------------------------------!
         call veg_energy_balance(tissue%abs_sw(i), tissue%abs_lw(i), tissue%h_coeff_leaf(i), le_slope,          &
                                    lw_slope, le_ref, tcas, tw, tissue%leaf_hcap_per_dt(i), tissue%t_leaf0(i),    &
                                    dtl, tl, transp_i, dh, drnet, q_extra=tissue%qwflux_wl(i),        &
                                    f_wet=film%f_wet_c(i), le_slope_wet=le_slope_wet,               &
                                    le_ref_wet=le_ref_wet, film_evap=f%film_evap_leaf(i),          &
                                    h_evap=h_evap_l, h_evap_wet=h_film_l)
         f%leaf_temp(i) = tl
         f%transp_c(i)  = transp_i                                          ! per-cohort transpiration this stage charges
         coh_h      = coh_h      + dh
         coh_qw     = coh_qw     + (transp_i + f%film_evap_leaf(i)) * h_evap_l   ! what the leaf paid
         coh_transp = coh_transp + transp_i
         coh_film_evap = coh_film_evap + f%film_evap_leaf(i)
         coh_rnet   = coh_rnet   + drnet
         !----- The wood's sensible heat, film evaporation and net longwave join coh_h, coh_qw and   !
         !      coh_rnet after the leaf's. ---------------------------------------------------------!
         coh_h    = coh_h    + dh_w
         coh_qw   = coh_qw   + f%film_evap_wood(i) * h_evap_w   ! wood film-evap -> CAS, what the wood paid
         coh_film_evap = coh_film_evap + f%film_evap_wood(i)
         coh_rnet = coh_rnet + drnet_w
      end do

      !----- GROUND SURFACE = snowfac-blended snow + (1-snowfac) bare soil. The snow terms come from !
      !      the shared pre-column stage (meds_fast_frozen%advance_snow_stage) and are ALREADY       !
      !      snowfac-weighted; the bare-soil sensible is scaled by (1-snowfac). soil_evap already    !
      !      carries its own (1-snowfac) AREA factor via chydro_forcing_t%snow_free_frac, so le_soil !
      !      is an area-integrated tile flux and is added whole -- do NOT scale it again.            !
      !                                                                                              !
      !      With no snow every added term is 0 and snowfac = 0, so this reduces EXACTLY to          !
      !        h_ground = h_bare, le_ground = le_soil, g_top = rad - h_bare - le_soil                !
      !      and a snow-free column needs no branch. snow%ground_rad is seeded to abs_sw_ground +    !
      !      abs_lw_ground by advance_snow_stage for the same reason. -------------------------------!
      call ground_surface_fluxes(t_ground, tcas, ground%ggnet, cas%rho, ground%soil_evap, h_bare, le_soil)
      f%h_ground  = snow%h_snow  + (1.0_wp - snow%snowfac) * h_bare
      f%le_ground = snow%le_snow + le_soil
      f%g_top     = snow%g_base                                                                   &
                    + (1.0_wp - snow%snowfac) * (ground%abs_sw_ground + ground%abs_lw_ground - h_bare)       &
                    - le_soil

      f%src_enth   = coh_h + coh_qw + f%h_ground + f%le_ground
      f%src_vap    = coh_transp + coh_film_evap + ground%soil_evap + snow%subl_rate
      f%coh_rnet   = coh_rnet

      !----- SMOOTH condensation sink (dew/fog): when the CAS is supersaturated, condense the excess as a   !
      !      relaxation of q toward qsat with timescale TAU_COND -- continuous, with only a kink           !
      !      in slope at RH=1, where a HARD state clamp would jump. The condensed vapour leaves as liquid   !
      !      at Tcas: the CAS loses vapour mass (src_vap) AND only the LIQUID enthalpy (src_enth -=         !
      !      cond*u_liq); the latent heat (enthalpy_vapor - u_liq) STAYS in the air (warming it, as real    !
      !      condensation does). Being smooth keeps the ARK2 embedded error smooth through saturation.      !
      !      cas%cas_condensation ([fast].cas_condensation, default on) switches the sink off.              !
      if (cas%cas_condensation) then
         f%cond     = (cas%cas_mass_capacity / TAU_COND) * max(0.0_wp, y%cas_shv - qsat_c)
      else
         f%cond     = 0.0_wp
      end if
      f%src_vap  = f%src_vap  - f%cond
      f%cond_enth = f%cond * internal_energy_liquid(tcas)
      f%src_enth = f%src_enth - f%cond_enth

      !----- CAS box tendencies (shared kernel; the condensation adjustment above is folded into    !
      !      src_enth/src_vap, so the box math is scheme-agnostic). ---------------------------------!
      cas_src%surface_enthalpy_source = f%src_enth
      cas_src%surface_vapor_source    = f%src_vap
      cas_src%biotic_co2_source       = cas%nee_biotic
      cas_col%air_mass_capacity        = cas%cas_mass_capacity
      cas_col%air_molar_capacity       = cas%cas_molar_capacity
      !----- The CAS<->atm conductances follow the LIVE canopy-air state.  This is the RK45 path's     !
      !      whole story -- its CAS tendency is built here, so this call is what unfreezes the           !
      !      Monin-Obukhov feedback for that scheme.  ARK arrives with mo_live already cleared on its   !
      !      stage copy, so it reuses the value it solved once for this stage rather than re-solving on  !
      !      every one of the Newton's residual evaluations. ---------------------------------------------!
      call cas_conductances(cas, y%cas_enthalpy, y%cas_shv, gah_l, gaw_l, gac_l)
      cas_col%atm_conductance_enthalpy = gah_l
      cas_col%atm_conductance_vapor    = gaw_l
      cas_col%atm_conductance_co2      = gac_l
      cas_col%atm_enthalpy             = cas%enthalpy_atm
      cas_col%atm_specific_humidity    = cas%shv_atm
      cas_col%atm_co2                  = cas%co2_atm
      call cas_column_time_deriv(y%cas_enthalpy, y%cas_shv, y%cas_co2, cas_src, cas_col,        &
                                 f%d_cas_enthalpy, f%d_cas_shv, f%d_cas_co2)
   end subroutine surface_derivs

   !----- Allocate a surface-tendency record's per-cohort arrays for n cohorts, unless they are. -!
   pure subroutine ensure_surface_tend(f, n)
      type(surface_tend_t), intent(inout) :: f
      integer(ik),          intent(in)    :: n
      if (allocated(f%leaf_temp) .and. allocated(f%wood_temp) .and. allocated(f%transp_c) .and.   &
          allocated(f%film_evap_leaf) .and. allocated(f%film_evap_wood)) then
         if (size(f%leaf_temp) == n .and. size(f%wood_temp) == n .and. size(f%transp_c) == n .and. &
             size(f%film_evap_leaf) == n .and. size(f%film_evap_wood) == n) return
      end if
      if (allocated(f%leaf_temp))      deallocate(f%leaf_temp)
      if (allocated(f%wood_temp))      deallocate(f%wood_temp)
      if (allocated(f%transp_c))       deallocate(f%transp_c)
      if (allocated(f%film_evap_leaf)) deallocate(f%film_evap_leaf)
      if (allocated(f%film_evap_wood)) deallocate(f%film_evap_wood)
      allocate(f%leaf_temp(n), f%wood_temp(n), f%transp_c(n), f%film_evap_leaf(n), f%film_evap_wood(n))
   end subroutine ensure_surface_tend

   !---------------------------------------------------------------------------------------!
   ! column_derivs -- the WHOLE-column RHS. Diagnoses the soil-top temperature from the state (so    !
   ! the ground skin couples to the current soil-top energy), runs the surface block, then assembles !
   ! the soil-heat, soil-water, and per-cohort plant-water-MASS tendencies from the frozen forcing +   !
   ! the surface couplings. Commits nothing. It is the tendency an RK45 stage evaluates.             !
   !                                                                                                 !
   ! WATER is the soil-plant water network's (docs/science/numerical_scheme.md §5b,                  !
   ! docs/science/soil_plant_water.md §4): soil water, sapflow and root uptake move at the network's !
   ! time-mean rates, and the leaves transpire what the network gave them (surface_derivs holds it), !
   ! so every stage sees the same water and the step ends at the network's state. Only the energy    !
   ! and the canopy air are integrated here. --------------------------------------------------------!
   pure subroutine column_derivs(y, frozen, col_config, n, nsl, f, sf_out)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      type(column_config_t), intent(in)  :: col_config   !< the column's parameters (soil, thermal, hydraulics)
      integer(ik),           intent(in)  :: n      !< number of cohorts
      integer(ik),           intent(in)  :: nsl    !< number of active soil layers
      type(column_tend_t),   intent(out) :: f
      !----- this stage's OWN surface tendencies (coh_rnet, src_enth/src_vap, ...), for a fully-    !
      !      explicit multi-stage integrator (RK45) to b-weight into its OWN whole-column boundary- !
      !      flux ledger -- mirrors column_be_stage's surf_tend output (meds_fast_be_stage). -------!
      type(surface_tend_t),  intent(out), optional :: sf_out

      type(surface_state_t)      :: y_stage
      type(surface_tend_t)       :: surf_tend
      type(soil_energy_column_t) :: soil_energy
      type(energy_forcing_t)     :: eforc
      real(wp)                   :: t_ground, fliq1, wmass1, transp_i, e_drain, water_heat(nsl)
      integer(ik)                :: i

      allocate(f%d_leaf_water_mass(n), f%d_wood_water_mass(n), f%leaf_temp(n))
      allocate(f%d_leaf_surf_water(n), f%d_wood_surf_water(n))

      !----- Diagnose the soil-top temperature from the current state so the ground skin sees the   !
      !      prognostic soil-top energy (the coupling the surface block needs). ---------------------!
      wmass1   = y%theta(1) * rho_h2o
      call internal_energy_to_temp(y%soil_energy(1), wmass1, col_config%soil_thermal%soil_dry_heat_capacity(1), t_ground, fliq1)

      !----- 1. Surface block (leaf + ground + CAS twins). ------------------------------------!
      y_stage%cas_enthalpy = y%cas_enthalpy ; y_stage%cas_shv = y%cas_shv ; y_stage%cas_co2 = y%cas_co2
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,    &
                          t_ground, n, surf_tend)
      f%d_cas_enthalpy = surf_tend%d_cas_enthalpy ; f%d_cas_shv = surf_tend%d_cas_shv ; f%d_cas_co2 = surf_tend%d_cas_co2
      f%g_top = surf_tend%g_top ; f%leaf_temp(1:n) = surf_tend%leaf_temp(1:n)

      !----- 2. Soil water moves at the soil-plant water network's steady rate over the step: the     !
      !      network is the one authority for soil water on every scheme, as it is on ARK, so RK45's    !
      !      theta ends the step exactly where the network's did (a steady rate integrates exactly     !
      !      under any consistent tableau). Its drainage and root uptake are the network's too. ------!
      f%dtheta_dt(1:nsl) = frozen%hydrology%theta_rate(1:nsl)
      f%drainage_rate    = frozen%hydrology%drainage
      f%uptake_rate      = frozen%roots%uptake

      !----- 3. Soil-heat column: g_top from the surface; the water's enthalpy moves with the         !
      !      network's flows -- its time-mean faces, infiltration at the pond's temperature, root     !
      !      uptake and the top layer's exfiltration at each layer's start-of-step temperature, and   !
      !      drainage at the bottom layer's -- the same numbers the ARK stage uses, and the ones the  !
      !      committed theta actually moved (docs/science/numerical_scheme.md §5b). The leaf pays the !
      !      full vapour enthalpy of what it transpires in surface_derivs, as ED2 does. --------------!
      soil_energy%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      water_heat(1:nsl) = frozen%roots%root_heat(1:nsl)
      water_heat(1)     = water_heat(1) + frozen%hydrology%exfil_enth
      e_drain = f%drainage_rate * internal_energy_liquid(frozen%hydrology%t_bot)
      call assemble_soil_energy_forcing(eforc, nsl, surf_tend%g_top, frozen%hydrology%geothermal, y%theta,          &
                                        water_heat, frozen%hydrology%w_flux_frozen,                                &
                                        frozen%hydrology%infiltration, frozen%hydrology%t_infil, e_drain)
      !----- Record the faces the ENERGY equation was just given, so the stepper can check them     !
      !      against the theta it actually commits. Taken from eforc -- the assembled forcing -- so !
      !      the check keeps its meaning if the argument above changes. ----------------------------!
      f%soil_face(1:nsl) = -eforc%w_flux(1:nsl)      ! eforc negates; store DOWNWARD
      call soil_energy_time_deriv(soil_energy, eforc, col_config%soil_thermal, col_config%soil, col_config%energy, &
                                  f%dedt, bottom_face=f%g_bottom)

      !----- 4. Per-cohort plant WATER MASS: the network's sapflow and uptake in, and out the         !
      !      transpiration this stage charged the canopy air (surf_tend%transp_c, per m2 of ground;     !
      !      per plant divides by nplant) -- the network's own when the stages hold it, so the step     !
      !      ends at the network's water. ---------------------------------------------------------!
      do i = 1_ik, n
         transp_i = surf_tend%transp_c(i) / max(frozen%plant%nplant(i), tiny_num)
         f%d_leaf_water_mass(i) = frozen%plant%sapflow_frozen(i) - transp_i
         f%d_wood_water_mass(i) = frozen%plant%uptake_frozen(i)  - frozen%plant%sapflow_frozen(i)
      end do

      !----- 5. Canopy-SURFACE water: frozen interception in, REFRESHED per-stage film evaporation out  !
      !      -- the surface-film analogue of the internal mass ODE just above. Already per m2 of ground !
      !      on both sides (unlike leaf_water_mass, no nplant division: frozen%film%intercept_leaf/wood !
      !      and surf_tend%film_evap_leaf/wood share that convention). Zero when canopy_water_on is     !
      !      off, so this is a no-op unless build_column_frozen populates it. --------------------------!
      f%d_leaf_surf_water(1:n) = frozen%film%intercept_leaf(1:n) - surf_tend%film_evap_leaf(1:n)
      f%d_wood_surf_water(1:n) = frozen%film%intercept_wood(1:n) - surf_tend%film_evap_wood(1:n)
      if (present(sf_out)) sf_out = surf_tend
   end subroutine column_derivs

end module meds_fast_time_derivs
