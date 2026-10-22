! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_time_derivs -- the side-effect-free whole-column tendency f(y) = dy/dt the fast-loop     !
! integrators share (docs/science/numerical_scheme.md §3): no flux hidden in a backward-Euler        !
! denominator and no store mutated, so any tableau can combine its evaluations.                      !
!   * surface_derivs   -- leaf and wood energy balance, ground skin, the three canopy-air twins;     !
!   * column_derivs    -- surface_derivs plus the soil-heat column and the soil and plant water,     !
!     which move at the soil-plant water network's frozen rates (§2);                                !
!   * cas_conductances -- the CAS <-> atmosphere conductances at a live canopy-air state.           !
!==========================================================================================!
module meds_fast_time_derivs
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, rho_h2o
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
   use meds_fast_types, only : snow_stage_t, surface_state_t, surface_tend_t, cas_boundary_t,                     &
                               tissue_coefficients_t, canopy_film_capacity_t, ground_boundary_t,                  &
                               column_state_t, column_frozen_t, column_tend_t, column_site_t
   implicit none
   private

   public :: surface_derivs, column_derivs, cas_conductances

   !----- [s] e-folding time of the smooth condensation sink, dq/dt = -(q - qsat)/TAU_COND: short       !
   !      against dt_fast, so supersaturated canopy air condenses within the step, and smooth (no clamp), !
   !      so the embedded error stays small. Any strong restoring rate pins q near qsat. ----------------!
   real(wp), parameter :: TAU_COND = 300.0_wp

contains

   !---------------------------------------------------------------------------------------!
   ! refresh_cas_conductances -- re-solve only the bulk Monin-Obukhov surface layer at a live canopy-  !
   ! air state for the three CAS <-> atmosphere conductances (docs/science/canopy_aerodynamics.md      !
   ! eq. (5)). The per-cohort boundary layers stay frozen (docs/science/numerical_scheme.md §2), as do !
   ! the geometry and the reference-level state: only (T_cas, q_cas) is live. Vapour and CO2 reuse the !
   ! heat transfer factor (z0q = z0h), CO2 on the molar capacity at the live humidity.                 !
   !---------------------------------------------------------------------------------------!
   pure subroutine refresh_cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      type(cas_boundary_t),   intent(in)  :: cas
      real(wp),               intent(in)  :: cas_enthalpy, cas_shv
      real(wp),               intent(out) :: g_atm_heat, g_atm_vapour, g_atm_co2
      real(wp) :: tcas, ustar, temp1, zeta, rib, obu, can_dmol
      tcas = cas_enthalpy_to_temp(cas_enthalpy, cas_shv)
      !----- The pre-pass also passes can_theta = tcas, so only the state differs from its solve. -!
      call mo_surface_layer(cas%aero_cfg, cas%mo_u_ref, cas%mo_zref, cas%mo_displace, cas%mo_rough, &
                            cas%mo_theta_atm, cas%mo_shv_atm, tcas, cas_shv,                      &
                            ustar, temp1, zeta, rib, obu)
      can_dmol = cas_molar_density(cas%mo_rho, cas_shv)
      call cas_atm_conductances(cas%mo_rho, can_dmol, ustar, temp1, temp1, g_atm_heat, g_atm_vapour, g_atm_co2)
   end subroutine refresh_cas_conductances

   !---------------------------------------------------------------------------------------!
   ! cas_conductances -- the one place that decides whether to re-solve the surface layer: a fresh    !
   ! solve when the bundle carries live MO inputs, else the bundle's own values. One rule keeps a      !
   ! tendency and the boundary flux charged against it in agreement (stage_bnd, meds_fast_rk45).       !
   !---------------------------------------------------------------------------------------!
   pure subroutine cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      type(cas_boundary_t),   intent(in)  :: cas
      real(wp),               intent(in)  :: cas_enthalpy, cas_shv
      real(wp),               intent(out) :: g_atm_heat, g_atm_vapour, g_atm_co2
      if (cas%mo_live) then
         call refresh_cas_conductances(cas, cas_enthalpy, cas_shv, g_atm_heat, g_atm_vapour, g_atm_co2)
      else
         g_atm_heat = cas%g_atm_heat
         g_atm_vapour = cas%g_atm_vapour
         g_atm_co2 = cas%g_atm_co2
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
      real(wp),               intent(in)  :: t_ground   !< [K] soil-top temperature at this evaluation (live, not frozen)
      integer(ik),            intent(in)  :: n
      !----- Filled in place; the per-cohort arrays are allocated only when missing or resized, so a  !
      !      Newton solve that evaluates this many times on one record allocates them once. ---------!
      type(surface_tend_t),   intent(inout) :: f

      real(wp)    :: tcas, qcas, qsat_c, dqdt
      real(wp)    :: lw_slope, le_slope, le_ref, dtl, tl, transp_i, dh, drnet
      real(wp)    :: lw_slope_w, dtw, tw, transp_w, dh_w, drnet_w    !< wood balance (own store)
      real(wp)    :: coh_h, coh_qw, coh_transp, coh_rnet, coh_film_evap
      real(wp)    :: h_evap_l, h_film_l, h_evap_w, h_film_w   !< [J/kg] energy per kg evaporated (leaf/wood; transp/film)
      real(wp)    :: h_bare, le_soil   !< bare-soil half of the snowfac blend
      !----- Film-evaporation latent terms: frozen conductance, live dqdt and qsat_c - qcas, as for the !
      !      dry pathway. With canopy_water_on off, f_wet_c = 0 and the wet pathway vanishes. --------!
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

      coh_h = 0.0_wp
      coh_qw = 0.0_wp
      coh_transp = 0.0_wp
      coh_rnet = 0.0_wp
      coh_film_evap = 0.0_wp
      do i = 1_ik, n
         lw_slope = lw_emission_slope(tissue%leaf_emiss(i), tcas, tissue%lai(i))
         !----- The leaf pays the full specific enthalpy of the vapour it sheds, enthalpy_vapor at tcas   !
         !      (the linearization base), and the CAS receives the same number; the film pays it less the !
         !      liquid enthalpy its water arrived with, so film, leaf and CAS close with no slack         !
         !      (docs/science/vegetation_energy_dynamics.md eq. (2)). -----------------------------------!
         h_evap_l = enthalpy_vapor(tcas)
         h_film_l = h_evap_l - film%film_liquid_enthalpy
         if (tissue%transp_held) then
            !----- Transpiration held at the network's: a fixed latent term on the leaf's dry share, so  !
            !      the transpiration it reports is the water the leaves lost (soil_plant_water.md §4). -!
            le_slope = 0.0_wp
            le_ref   = h_evap_l * tissue%transp_net(i) / max(1.0_wp - film%f_wet_c(i), tiny_num)
         else
            le_slope = h_evap_l * cas%rho * tissue%g_transp_leaf(i) * dqdt
            le_ref   = h_evap_l * cas%rho * tissue%g_transp_leaf(i) * (qsat_c - qcas)
         end if
         le_slope_wet = h_film_l * cas%rho * film%g_film_leaf(i) * dqdt
         le_ref_wet   = h_film_l * cas%rho * film%g_film_leaf(i) * (qsat_c - qcas)
         !----- Wood: the same kernel with no transpiration, solved first because a leaf with no area    !
         !      takes its wood's temperature (docs/science/vegetation_energy_dynamics.md, "A tissue with !
         !      no area"). q_wood_net enters as q_extra, as qwflux_wl does for the leaf. -----------------!
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
         !----- Leaf: linearized about tcas, relaxing with its own heat store                            !
         !      (docs/science/vegetation_energy_dynamics.md eqs. (3)-(6)). qwflux_wl enters as q_extra:  !
         !      it shifts the equilibrium like absorbed energy but is an internal transfer, not in drnet. !
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

      !----- Ground surface = snowfac-weighted snow + (1 - snowfac) bare soil. The snow terms arrive    !
      !      snowfac-weighted (advance_snow_stage); soil_evap already carries its snow-free area factor   !
      !      (chydro_forcing_t%snow_free_frac), so le_soil is added whole. With no snow this is exactly  !
      !      bare ground, with no branch. --------------------------------------------------------------!
      call ground_surface_fluxes(t_ground, tcas, ground%ggnet, cas%rho, ground%soil_evap, h_bare, le_soil)
      f%h_ground  = snow%h_snow  + (1.0_wp - snow%snowfac) * h_bare
      f%le_ground = snow%le_snow + le_soil
      f%g_top     = snow%g_base                                                                   &
                    + (1.0_wp - snow%snowfac) * (ground%abs_sw_ground + ground%abs_lw_ground - h_bare)       &
                    - le_soil

      f%src_enth   = coh_h + coh_qw + f%h_ground + f%le_ground
      f%src_vap    = coh_transp + coh_film_evap + ground%soil_evap + snow%subl_rate
      f%coh_rnet   = coh_rnet

      !----- Smooth condensation sink (dew, fog): supersaturated canopy air relaxes to qsat over TAU_COND. !
      !      The condensate leaves as liquid at tcas, taking only its liquid enthalpy; the latent heat    !
      !      stays in the air, as in real condensation. ----------------------------------------------!
      f%cond     = (cas%cas_mass_capacity / TAU_COND) * max(0.0_wp, y%cas_shv - qsat_c)
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
      !----- The CAS <-> atmosphere conductances follow the live canopy-air state: RK45 re-solves them   !
      !      here. ARK arrives with mo_live cleared on its stage copy, having solved once for the stage,  !
      !      so its Newton's residual evaluations reuse that value. -----------------------------------!
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
   ! column_derivs -- the whole-column RHS an RK45 stage evaluates; commits nothing. Diagnoses the    !
   ! soil-top temperature from the state, runs the surface block, and assembles the soil-heat, soil-  !
   ! water and plant-water-mass tendencies. Water moves at the soil-plant water network's time-mean   !
   ! rates and the leaves transpire what the network gave them, so every stage sees the same water    !
   ! (docs/science/soil_plant_water.md §4, docs/science/numerical_scheme.md §5b).                     !
   !---------------------------------------------------------------------------------------!
   pure subroutine column_derivs(y, frozen, col_site, n, nsl, f, sf_out)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      type(column_site_t),   intent(in)  :: col_site
      integer(ik),           intent(in)  :: n      !< number of cohorts
      integer(ik),           intent(in)  :: nsl    !< number of active soil layers
      type(column_tend_t),   intent(out) :: f
      !----- This stage's own surface tendencies, for RK45's whole-column boundary-flux ledger (as   !
      !      column_be_stage's surf_tend serves ARK's). ----------------------------------------------!
      type(surface_tend_t),  intent(out), optional :: sf_out

      type(surface_state_t)      :: y_stage
      type(surface_tend_t)       :: surf_tend
      type(soil_energy_column_t) :: soil_energy
      type(energy_forcing_t)     :: eforc
      real(wp)                   :: t_ground, fliq1, wmass1, transp_i, e_drain, water_heat(nsl)
      integer(ik)                :: i

      allocate(f%d_leaf_water_mass(n), f%d_wood_water_mass(n), f%leaf_temp(n))
      allocate(f%d_leaf_surf_water(n), f%d_wood_surf_water(n))

      !----- The soil-top temperature from the current state, for the ground skin. -----------------!
      wmass1   = y%theta(1) * rho_h2o
      call internal_energy_to_temp(y%soil_energy(1), wmass1, col_site%soil_profile%dry_heat_capacity(1), t_ground, fliq1)

      !----- 1. Surface block (leaf + ground + CAS twins). ------------------------------------!
      y_stage%cas_enthalpy = y%cas_enthalpy
      y_stage%cas_shv = y%cas_shv
      y_stage%cas_co2 = y%cas_co2
      call surface_derivs(y_stage, frozen%cas, frozen%tissue, frozen%film, frozen%ground, frozen%snow,    &
                          t_ground, n, surf_tend)
      f%d_cas_enthalpy = surf_tend%d_cas_enthalpy
      f%d_cas_shv = surf_tend%d_cas_shv
      f%d_cas_co2 = surf_tend%d_cas_co2
      f%g_top = surf_tend%g_top
      f%leaf_temp(1:n) = surf_tend%leaf_temp(1:n)

      !----- 2. Soil water at the network's steady rate, which any consistent tableau integrates      !
      !      exactly, so RK45's theta ends where the network's does; so do drainage and uptake. ----!
      f%dtheta_dt(1:nsl) = frozen%hydrology%theta_rate(1:nsl)
      f%drainage_rate    = frozen%hydrology%drainage
      f%uptake_rate      = frozen%roots%uptake

      !----- 3. Soil heat: g_top from the surface, and the enthalpy of the network's flows, the numbers !
      !      the ARK stage uses (docs/science/numerical_scheme.md §5b, soil_plant_water.md §5). The   !
      !      leaf pays its transpiration's full vapour enthalpy in surface_derivs. ------------------!
      soil_energy%soil_energy(1:nsl) = y%soil_energy(1:nsl)
      water_heat(1:nsl) = frozen%roots%root_heat(1:nsl)
      water_heat(1)     = water_heat(1) + frozen%hydrology%exfil_enth
      e_drain = f%drainage_rate * internal_energy_liquid(frozen%hydrology%t_bot)
      call assemble_soil_energy_forcing(eforc, nsl, surf_tend%g_top, y%theta,          &
                                        water_heat, frozen%hydrology%w_flux_frozen,                                &
                                        frozen%hydrology%infiltration, frozen%hydrology%t_infil, e_drain)
      !----- The faces the energy equation was given, read from eforc, for the stepper to check      !
      !      against the theta it commits. ---------------------------------------------------------!
      f%soil_face(1:nsl) = -eforc%w_flux(1:nsl)      ! eforc negates; store downward
      call soil_energy_time_deriv(soil_energy, eforc, col_site%soil_profile, col_site%energy, &
                                  f%dedt, bottom_face=f%g_bottom)

      !----- 4. Plant water mass per plant: the network's sapflow and uptake in, and out the          !
      !      transpiration this stage charged the canopy air, the network's own when the stages hold it. !
      do i = 1_ik, n
         transp_i = surf_tend%transp_c(i) / max(frozen%plant%nplant(i), tiny_num)
         f%d_leaf_water_mass(i) = frozen%plant%sapflow_frozen(i) - transp_i
         f%d_wood_water_mass(i) = frozen%plant%uptake_frozen(i)  - frozen%plant%sapflow_frozen(i)
      end do

      !----- 5. Canopy film: frozen interception in, per-stage film evaporation out, both per m2 of   !
      !      ground. Zero with canopy_water_on off. -------------------------------------------------!
      f%d_leaf_surf_water(1:n) = frozen%film%intercept_leaf(1:n) - surf_tend%film_evap_leaf(1:n)
      f%d_wood_surf_water(1:n) = frozen%film%intercept_wood(1:n) - surf_tend%film_evap_wood(1:n)
      if (present(sf_out)) sf_out = surf_tend
   end subroutine column_derivs

end module meds_fast_time_derivs
