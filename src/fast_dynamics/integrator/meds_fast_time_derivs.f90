! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_time_derivs -- the side-effect-free whole-column tendency f(y) = dy/dt the fast-loop     !
! integrators share (docs/science/numerical_scheme.md §3): no flux hidden in a backward-Euler        !
! denominator and no store mutated, so any tableau can combine its evaluations.                      !
!   * surface_derivs   -- leaf and wood energy balance, ground skin, the three canopy-air twins;     !
!   * column_derivs    -- surface_derivs plus the soil-heat column and the soil and plant water,     !
!     which move at the soil-plant water network's frozen rates (§2);                                !
!   * cas_conductances -- the CAS <-> atmosphere conductances at a live canopy-air state;          !
!   * water_derivs     -- the soil-plant water network's flows, inflows and Jacobian at one state,  !
!     and set_water_state its state from the model's (docs/science/soil_plant_water.md §1-§2).      !
!==========================================================================================!
module meds_fast_time_derivs
   use meds_kinds, only : wp, ik
   use meds_constants, only : tiny_num, rho_h2o, grav_head
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
                               column_state_t, column_frozen_t, column_tend_t, fast_static_state_t,               &
                               water_drive_t, water_tend_t, water_work_t, ensure_water_capacity, POND, WOOD, LEAF
   use meds_site_state_vars, only : site_surface_t
   use meds_soil_hydr_lib, only : soil_water_relations, soil_theta_to_unknown, pond_water_relations,              &
                                  pond_water_to_head
   use meds_plant_hydr_lib, only : tissue_water_relations, tissue_water_to_rwc, tissue_rwc_to_unknown,            &
                                   tissue_psi_to_rwc, root_water_uptake_profile, xylem_sapflow
   use meds_plant_hydraulics, only : supply_limited_transpiration
   use meds_soil_water, only : soil_water_flux, ground_evaporation_rate, bottom_flux
   use meds_linear_algebra, only : arrowhead_system_t
   use meds_soil_params, only : soil_profile_t, soil_params_t
   use meds_pft_params, only : hydraulics_table_t
   implicit none
   private

   public :: surface_derivs, column_derivs, cas_conductances
   public :: set_water_state, water_derivs

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
         g_atm_heat = cas%box%atm_conductance_enthalpy
         g_atm_vapour = cas%box%atm_conductance_vapor
         g_atm_co2 = cas%box%atm_conductance_co2
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
      !      (ground_evap_env_t%snow_free_frac), so le_soil is added whole. With no snow this is exactly  !
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
      f%cond     = (cas%box%air_mass_capacity / TAU_COND) * max(0.0_wp, y%cas_shv - qsat_c)
      f%src_vap  = f%src_vap  - f%cond
      f%cond_enth = f%cond * internal_energy_liquid(tcas)
      f%src_enth = f%src_enth - f%cond_enth

      !----- CAS box tendencies (shared kernel; the condensation adjustment above is folded into    !
      !      src_enth/src_vap, so the box math is scheme-agnostic). ---------------------------------!
      cas_src%surface_enthalpy_source = f%src_enth
      cas_src%surface_vapor_source    = f%src_vap
      cas_src%biotic_co2_source       = cas%nee_biotic
      !----- The CAS <-> atmosphere conductances follow the live canopy-air state: RK45 re-solves them   !
      !      here. ARK arrives with mo_live cleared on its stage copy, having solved once for the stage,  !
      !      so its Newton's residual evaluations reuse that value. -----------------------------------!
      cas_col = cas%box
      call cas_conductances(cas, y%cas_enthalpy, y%cas_shv, gah_l, gaw_l, gac_l)
      cas_col%atm_conductance_enthalpy = gah_l
      cas_col%atm_conductance_vapor    = gaw_l
      cas_col%atm_conductance_co2      = gac_l
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
   pure subroutine column_derivs(y, frozen, site, n, nsl, f, sf_out)
      type(column_state_t),  intent(in)  :: y
      type(column_frozen_t), intent(in)  :: frozen
      type(site_surface_t),  intent(in)  :: site
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
      call internal_energy_to_temp(y%soil_energy(1), wmass1, site%soil_profile%dry_heat_capacity(1), t_ground, fliq1)

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
      call soil_energy_time_deriv(soil_energy, eforc, site%soil_profile, site%soil,   &
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

   !=======================================================================================!
   !  The soil-plant water network (docs/science/soil_plant_water.md §1-§2): the pond, the   !
   !  soil layers, and each cohort's wood and leaves. The stores' curves and the paths'      !
   !  flows are the process libraries'; these join them.                                    !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! set_water_state -- the network's state from the model's: each soil layer's water       !
   ! content, the pond's water and enthalpy, and each cohort's leaf and wood water [kg per   !
   ! plant]; each store's unknown is read from its water. A leafless cohort's leaves hold    !
   ! none, so their unknown starts at the wood's potential less the head to the crown, where !
   ! an empty leaf comes to rest.                                                           !
   !---------------------------------------------------------------------------------------!
   subroutine set_water_state(static, profile, hyd, theta, w_pond, e_pond, leaf_water, wood_water, water)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(hydraulics_table_t),  intent(in)    :: hyd
      real(wp),                  intent(in)    :: theta(:)                       !< (nsl) [m3/m3]
      real(wp),                  intent(in)    :: w_pond, e_pond                 !< [kg/m2], [J/m2]
      real(wp),                  intent(in)    :: leaf_water(:), wood_water(:)   !< (ncoh) [kg/plant]
      type(water_work_t),        intent(inout) :: water
      integer(ik) :: k, i
      real(wp)    :: psi, dpsi, w, dw, kk, dk, head1
      call ensure_water_capacity(water, profile%n_active, static%n)
      associate (st => water%state)
         do k = 1_ik, profile%n_active
            st%tau_soil(1_ik + k)   = soil_theta_to_unknown(profile%layer_curve(k), theta(k))
            st%water_soil(1_ik + k) = rho_h2o * profile%dz(k) * theta(k)
         end do
         !----- The pond starts at its water's head, but no higher than the head at which it neither   !
         !      gives nor takes water from the top layer: over dry soil a trace pond empties within the  !
         !      step, and its own head would start Newton at a flow thousands of times what it holds. --!
         st%water_soil(POND) = w_pond
         st%e_pond = e_pond
         call soil_water_relations(profile%layer_curve(1), st%tau_soil(2), psi, dpsi, w, dw, kk, dk)
         head1 = psi / grav_head
         st%tau_soil(POND) = min(pond_water_to_head(w_pond), head1 + profile%z_node(1))
         do i = 1_ik, static%n
            associate (hp => hyd%pft(static%pft(i)))
               st%tau_coh(WOOD, i)   = tissue_rwc_to_unknown(                                                  &
                                          tissue_water_to_rwc(hp%wood_curve, static%bwood_water(i), wood_water(i)))
               st%water_coh(WOOD, i) = wood_water(i) * static%nplant(i)
               if (hp%leaf_curve%water_sat * static%bleaf(i) > tiny_num) then
                  st%tau_coh(LEAF, i)   = tissue_rwc_to_unknown(                                               &
                                             tissue_water_to_rwc(hp%leaf_curve, static%bleaf(i), leaf_water(i)))
                  st%water_coh(LEAF, i) = leaf_water(i) * static%nplant(i)
               else
                  call tissue_water_relations(hp%wood_curve, static%bwood_water(i), static%nplant(i),            &
                                              st%tau_coh(WOOD, i), psi, dpsi, w, dw)
                  st%tau_coh(LEAF, i)   = tissue_rwc_to_unknown(tissue_psi_to_rwc(psi - static%crown_head(i),   &
                                             hp%leaf_curve%pi0, hp%leaf_curve%elastic_mod))
                  st%water_coh(LEAF, i) = 0.0_wp
               end if
            end associate
         end do
      end associate
   end subroutine set_water_state

   !---------------------------------------------------------------------------------------!
   ! water_derivs -- the network at one state: each store's potential and water from its     !
   ! curve; the flow along every path (pond to top layer, between layers, evaporation,       !
   ! drainage, each cohort's uptake, sapflow and transpiration); the net inflow into each    !
   ! store; and, when asked, the flows' Jacobian in the stores' unknowns (sys).              !
   !---------------------------------------------------------------------------------------!
   subroutine water_derivs(static, profile, soil, hyd, drive, tau_soil, tau_coh, ev, sys, with_jacobian)
      type(fast_static_state_t), intent(in)    :: static
      type(soil_profile_t),      intent(in)    :: profile
      type(soil_params_t),       intent(in)    :: soil
      type(hydraulics_table_t),  intent(in)    :: hyd
      type(water_drive_t),       intent(in)    :: drive
      real(wp),                  intent(in)    :: tau_soil(:), tau_coh(:,:)
      type(water_tend_t),        intent(inout) :: ev
      type(arrowhead_system_t),  intent(inout) :: sys   !< the flows' Jacobian, written when with_jacobian
      logical,                   intent(in)    :: with_jacobian
      integer(ik) :: k, i, ca, cb, c, nsl, nsoil
      real(wp) :: q, dq_pa, dq_pb, dq_ka, dq_kb, dqa, dqb, u, s, ds_w, ds_l
      real(wp) :: e, de_dtheta, de_dpsi, de, qbot, dqbot_dpsi, dqbot_dk, theta1, head
      real(wp) :: psi_w, dpsi_w, u_sum, dqb_sum, du_dpsi(profile%n_active), du_dk(profile%n_active)

      nsl   = profile%n_active
      nsoil = 1_ik + nsl

      ev%inflow_soil = 0.0_wp
      ev%inflow_coh = 0.0_wp
      !----- The Jacobian's entries are written afresh: the soil block's band, which the flows add  !
      !      into, is cleared; every other entry an evaluation touches is set outright below, and the  !
      !      rest stay zero from arrowhead_allocate. ---------------------------------------------------!
      if (with_jacobian) then
         do c = 1_ik, nsoil
            sys%a(c, c) = 0.0_wp
            if (c < nsoil) then
               sys%a(c, c + 1_ik) = 0.0_wp
               sys%a(c + 1_ik, c) = 0.0_wp
            end if
         end do
      end if

      !----- Storage. ------------------------------------------------------------------!
      call pond_water_relations(tau_soil(POND), ev%psi_soil(1), ev%dpsi_soil(1), ev%water_soil(1),             &
                                ev%dwater_soil(1))
      ev%k_soil(1) = profile%ksat(1)
      ev%dk_soil(1) = 0.0_wp      ! the pond conducts as saturated soil
      do k = 1_ik, nsl
         c = 1_ik + k
         call soil_water_relations(profile%layer_curve(k), tau_soil(c), ev%psi_soil(c), ev%dpsi_soil(c),          &
                                   ev%water_soil(c), ev%dwater_soil(c), ev%k_soil(c), ev%dk_soil(c))
      end do
      do i = 1_ik, static%n
         associate (hp => hyd%pft(static%pft(i)))
            call tissue_water_relations(hp%wood_curve, static%bwood_water(i), static%nplant(i), tau_coh(WOOD, i), &
                                        ev%psi_coh(WOOD, i), ev%dpsi_coh(WOOD, i), ev%water_coh(WOOD, i),        &
                                        ev%dwater_coh(WOOD, i))
            call tissue_water_relations(hp%leaf_curve, static%bleaf(i), static%nplant(i), tau_coh(LEAF, i),       &
                                        ev%psi_coh(LEAF, i), ev%dpsi_coh(LEAF, i), ev%water_coh(LEAF, i),        &
                                        ev%dwater_coh(LEAF, i))
         end associate
      end do

      !----- The pond and the top layer: Darcy, the pond a saturated layer as thick as the top one. !
      ca = POND
      cb = 2_ik
      call soil_water_flux(ev%k_soil(ca), ev%k_soil(cb), 0.5_wp, ev%psi_soil(ca), ev%psi_soil(cb),                &
                           -profile%z_node(1), -profile%z_node(1), q, dq_pa, dq_pb, dq_ka, dq_kb)
      ev%pond_flux = q
      call add_soil_flow(ca, cb, q, dq_pa * ev%dpsi_soil(ca), dq_pb * ev%dpsi_soil(cb) + dq_kb * ev%dk_soil(cb))

      !----- Between soil layers. --------------------------------------------------------!
      do k = 1_ik, nsl - 1_ik
         ca = 1_ik + k
         cb = 2_ik + k
         call soil_water_flux(ev%k_soil(ca), ev%k_soil(cb), profile%face_weight(k), ev%psi_soil(ca),              &
                              ev%psi_soil(cb), profile%dz_node(k), profile%dz_node(k), q, dq_pa, dq_pb, dq_ka, dq_kb)
         ev%face_flux(k) = q
         dqa = dq_pa * ev%dpsi_soil(ca) + dq_ka * ev%dk_soil(ca)
         dqb = dq_pb * ev%dpsi_soil(cb) + dq_kb * ev%dk_soil(cb)
         call add_soil_flow(ca, cb, q, dqa, dqb)
      end do

      !----- Rain into the pond; evaporation from the top layer; drainage from the bottom one. ----!
      ev%inflow_soil(POND) = ev%inflow_soil(POND) + drive%rain
      c = 2_ik
      theta1 = ev%water_soil(c) / (rho_h2o * profile%dz(1))
      call ground_evaporation_rate(theta1, ev%psi_soil(c) / grav_head, profile, drive%ground, soil, e, de_dtheta,  &
                                   de_dpsi)
      ev%evap = e
      ev%inflow_soil(c) = ev%inflow_soil(c) - e
      if (with_jacobian) sys%a(c, c) = sys%a(c, c)                                                              &
                         + de_dtheta * ev%dwater_soil(c) / (rho_h2o * profile%dz(1))                               &
                         + de_dpsi * ev%dpsi_soil(c) / grav_head
      c = 1_ik + nsl
      head = ev%psi_soil(c) / grav_head
      call bottom_flux(profile, soil, nsl, head, ev%k_soil(c), qbot, dqbot_dpsi, dqbot_dk)
      ev%drain = rho_h2o * qbot
      ev%inflow_soil(c) = ev%inflow_soil(c) - ev%drain
      if (with_jacobian) sys%a(c, c) = sys%a(c, c)                                                              &
                         + rho_h2o * (dqbot_dpsi * ev%dpsi_soil(c) / grav_head + dqbot_dk * ev%dk_soil(c))

      !----- Each cohort: uptake from its layers, sapflow to its leaves, and transpiration. The      !
      !      cohort's uptake comes from its whole root profile at once; its own sums are kept in       !
      !      scalars over its layers and written once. ------------------------------------------!
      do i = 1_ik, static%n
         psi_w = ev%psi_coh(WOOD, i)
         dpsi_w = ev%dpsi_coh(WOOD, i)
         call root_water_uptake_profile(nsl, static%g_rhizo(1:nsl, i), ev%k_soil(2:nsoil), ev%psi_soil(2:nsoil), &
                                        profile%z_node(1:nsl), psi_w, ev%uptake(1:nsl, i), du_dpsi, du_dk)
         u_sum = 0.0_wp
         dqb_sum = 0.0_wp
         do k = 1_ik, nsl
            c = 1_ik + k
            if (static%g_rhizo(k, i) <= 0.0_wp) then
               if (with_jacobian) then
                  sys%b(c, i) = 0.0_wp
                  sys%c(c, i) = 0.0_wp
               end if
               cycle
            end if
            u = ev%uptake(k, i)
            ev%inflow_soil(c) = ev%inflow_soil(c) - u
            u_sum = u_sum + u
            if (with_jacobian) then
               dqa = du_dpsi(k) * ev%dpsi_soil(c) + du_dk(k) * ev%dk_soil(c)
               dqb = -du_dpsi(k) * dpsi_w
               sys%a(c, c) = sys%a(c, c) + dqa
               sys%b(c, i) = dqb
               sys%c(c, i) = -dqa
               dqb_sum = dqb_sum - dqb
            end if
         end do
         ev%inflow_coh(WOOD, i) = ev%inflow_coh(WOOD, i) + u_sum
         associate (hp => hyd%pft(static%pft(i)))
            call xylem_sapflow(static%k_xylem(i), ev%psi_coh(WOOD, i), ev%psi_coh(LEAF, i), static%crown_head(i),  &
                               hp%wood_psi50, hp%wood_kexp, hp%vuln_table, s, ds_w, ds_l)
         end associate
         ev%sapflow(i) = s
         ev%inflow_coh(WOOD, i) = ev%inflow_coh(WOOD, i) - s
         ev%inflow_coh(LEAF, i) = ev%inflow_coh(LEAF, i) + s
         call supply_limited_transpiration(drive%e_demand(i), tau_coh(LEAF, i), hyd%pft(static%pft(i))%rwc_close,   &
                                           e, de)
         ev%transp(i) = e
         ev%inflow_coh(LEAF, i) = ev%inflow_coh(LEAF, i) - e
         if (with_jacobian) then
            dqa = ds_w * ev%dpsi_coh(WOOD, i)
            dqb = ds_l * ev%dpsi_coh(LEAF, i)
            sys%d(WOOD, WOOD, i) = dqb_sum + dqa
            sys%d(WOOD, LEAF, i) = dqb
            sys%d(LEAF, WOOD, i) = -dqa
            sys%d(LEAF, LEAF, i) = -dqb + de
         end if
      end do

   contains

      !----- A flow q from soil node a to soil node b, with dq/dtau_a and dq/dtau_b: the inflows   !
      !      and, when wanted, the Jacobian of the residual (W - W0)/h - inflow. ----------------!
      subroutine add_soil_flow(a, b, q, dqa, dqb)
         integer(ik), intent(in) :: a, b
         real(wp),    intent(in) :: q, dqa, dqb
         ev%inflow_soil(a) = ev%inflow_soil(a) - q
         ev%inflow_soil(b) = ev%inflow_soil(b) + q
         if (.not. with_jacobian) return
         sys%a(a, a) = sys%a(a, a) + dqa
         sys%a(a, b) = sys%a(a, b) + dqb
         sys%a(b, a) = sys%a(b, a) - dqa
         sys%a(b, b) = sys%a(b, b) - dqb
      end subroutine add_soil_flow

   end subroutine water_derivs

end module meds_fast_time_derivs
