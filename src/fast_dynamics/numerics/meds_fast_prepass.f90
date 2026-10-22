! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_fast_prepass -- the once-per-dt_fast pre-pass of the fast loop (Category 0 of                !
! docs/science/numerical_scheme.md §2): canopy aerodynamics from the current canopy-air state,      !
! leaf gas exchange (GPP, gs, Rd) with the frozen leaf-energy coefficients, stem and fine-root      !
! maintenance respiration, heterotrophic respiration, and the canopy-air capacities and            !
! atmosphere conductances. column_prepass runs them in order and assembles the biotic CO2 flux;     !
! both integrators reach it through build_column_frozen, so ARK and RK45 see the same GPP.         !
!==========================================================================================!
module meds_fast_prepass
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, cp_air, r_gas, kgcday_2_umols
   use meds_config,           only : top_soil_psi
   use meds_diag_vars,        only : CO_ANET, CO_AGROSS, CO_GSW, CO_GBW, CO_CI, CO_CS, CO_RD, CO_BETA_STOM, CO_BETA_NONSTOM, &
                                     CO_LEAF_TEMP, CO_WOOD_TEMP, CO_LEAF_VPD, CO_PSI_LEAF, CO_ABS_PAR, CO_ABS_SW, CO_ABS_LW, &
                                     CO_WIND, CO_LEAF_WATER, CO_WOOD_WATER, CO_GPP_RATE, CO_LEAF_PAR, CO_GB_MOL, CO_PSI_PREDAWN, &
                                     CO_CA, CO_PRESSURE
   use meds_plant_hydraulics, only : leaf_node_psi
   use meds_canopy_types, only : aero_env_t, aero_geom_t, aero_out_t
   use meds_plant_types, only : leaf_photo_table_t, hydro_args_table_t
   use meds_fast_types, only : patch_biophys_t, column_config_t, column_cohort_t, column_forcing_t,               &
                               column_budget_t, column_site_t
   use meds_soil_params,      only : soil_profile_t
   use meds_state_defaults, only : is_unset_negative
   use meds_biophysics_params, only : aero_cfg_t
   use meds_patch_state_vars,   only : cas_state_t, soil_carbon_t
   use meds_biogeochem_types, only : N_SOIL_POOL
   use meds_canopy_aerodynamics, only : canopy_aerodynamics, cas_atm_conductances
   use meds_plant_biophysics, only : sensible_heat_coeff, leaf_transp_coeff, EFFAREA_HEAT, EFFAREA_TRANSP
   use meds_leaf_gas_exchange, only : leaf_gas_exchange_batch
   use meds_plant_respiration, only : stem_maintenance_respiration, fine_root_maintenance_respiration, &
                                     root_zone_temp_scale
   use meds_soil_biogeochem,  only : heterotrophic_respiration_matrix, &
                                     assemble_env_scalar, assemble_transfer_matrix
   use meds_therm_lib,        only : cas_molar_density, cas_enthalpy_to_temp, sat_vapor_pressure
   use meds_numerics,         only : weighted_mean, ascending_order
   implicit none
   private

   public :: column_prepass
   public :: cohort_leaf_par, aero_bottom_to_top   ! used only by the tests

contains

   !---------------------------------------------------------------------------------------!
   ! column_prepass -- run the five processes in order and assemble the biotic CO2 flux            !
   ! NEE_biotic = Rd + stem + root maintenance + Rh - GPP + slow_co2_rate [umol/m2/s], positive to  !
   ! the canopy air (docs/science/plant_respiration.md eq. (7)). biophys is read-only: a caller     !
   ! that persists the CAS temperature writes biophys%cas%can_temp = tcas itself.                   !
   !---------------------------------------------------------------------------------------!
   subroutine column_prepass(col_config, col_site, aenv, ageom, col_cohort, forc, biophys, aero, budget, root_frac, &
                             tcas, qcas, press, rho, t_ground, h_coeff_leaf, g_transp_leaf,                     &
                             cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2, nee_biotic, &
                             gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)
      type(column_config_t),   intent(in)    :: col_config
      type(column_site_t),     intent(in)    :: col_site
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: col_cohort
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(in)    :: biophys
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budget
      real(wp),                intent(in)    :: root_frac(:)   !< [-] the patch's fine roots by layer (sum 1)
      real(wp),                intent(out)   :: tcas, qcas, press, rho, t_ground
      real(wp),                intent(out)   :: h_coeff_leaf(:), g_transp_leaf(:)
      real(wp),                intent(out)   :: cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2, &
           nee_biotic
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      !----- Per-cohort diagnostic capture (docs/science/diagnostics.md §1, stage [2]); absent, the    !
      !      extra fields are never requested from the batch kernel. ----------------------------!
      real(wp), optional,      intent(inout) :: cdiag(:,:)   !< (N_CO_FAST_ROWS, ncoh) instantaneous values

      real(wp) :: gpp, ra_leaf, ra_stem, ra_root, rh, soil_temp_root, theta_mean

      !----- 1. aerodynamics from the current canopy-air state. ------------------------------------!
      call refresh_canopy_aerodynamics(col_config%aero, aenv, ageom, col_cohort, biophys%cas,           &
                                       biophys%soil_energy%soil_temp(1), biophys%leaf_temp,                 &
                                       biophys%wood_temp, aero, tcas, qcas, press, rho, t_ground)

      !----- 2. root-zone environment (root + heterotrophic respiration drivers). ------------------!
      call root_zone_environment(col_site%soil_profile, root_frac, biophys%soil_energy%soil_temp, biophys%soil_water%theta, &
                                 soil_temp_root, theta_mean)

      !----- 3. leaf gas exchange + the frozen leaf-energy coefficients. ---------------------------!
      call canopy_leaf_gas_exchange(col_site%leaf_photo, col_config%hydraulics_table, col_site%soil_profile,       &
                                    col_cohort, forc, aero, biophys,                                       &
                                    qcas, press, rho, gpp, ra_leaf, h_coeff_leaf, g_transp_leaf,                  &
                                    gpp_coh, leaf_resp_coh, cdiag)

      !----- 4. stem + fine-root maintenance respiration. -----------------------------------------!
      call canopy_maintenance_respiration(col_config%ea_rd, col_config%hd_rd, col_config%ds_rd, col_site%soil_profile,   &
                                          root_frac, col_cohort,                                                 &
                                          biophys%wood_temp, biophys%soil_energy%soil_temp, ra_stem, ra_root, &
                                          stem_resp_coh, root_resp_coh)

      !----- 5. heterotrophic respiration. -------------------------------------------------------!
      call patch_heterotrophic_respiration(col_config, col_site, biophys%soil_carbon, t_ground,    &
                                           soil_temp_root, theta_mean, rh, budget)

      !----- NEE (docs/science/plant_respiration.md eq. (7)). slow_co2_rate is the autotrophic part this !
      !      loop cannot compute, the allocator's daily growth respiration on realized growth, handed    !
      !      down as a frozen rate with the opposite-signed starvation correction (§5-§6 there). ------!
      nee_biotic = ra_leaf + ra_stem + ra_root + rh - gpp + biophys%slow_co2_rate
      budget%gpp_last = gpp
      budget%nee_last = nee_biotic

      !----- 6. CAS capacities + atm-exchange conductances (frozen across the macro-step). --------!
      call cas_capacities_and_conductances(rho, qcas, biophys%cas%can_depth, aero%ustar, aero%temp1,    &
                                           aero%temp2, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2)
   end subroutine column_prepass

   !---------------------------------------------------------------------------------------!
   ! refresh_canopy_aerodynamics -- the canopy-air scalars diagnosed from the CAS state, written    !
   ! into the aerodynamic environment, and the whole-canopy aerodynamics solve (friction velocity,  !
   ! scalar profile factors, ground conductance, per-cohort boundary layers;                        !
   ! docs/science/canopy_aerodynamics.md §2-§4).                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine refresh_canopy_aerodynamics(aero_cfg, aenv, ageom, col_cohort, cas, soil_temp_top,      &
                                          leaf_temp, wood_temp, aero, tcas, qcas, press, rho, t_ground)
      type(aero_cfg_t),      intent(in)    :: aero_cfg
      type(aero_env_t),      intent(inout) :: aenv
      type(aero_geom_t),     intent(in)    :: ageom
      type(column_cohort_t), intent(in)    :: col_cohort
      type(cas_state_t),     intent(in)    :: cas
      real(wp),              intent(in)    :: soil_temp_top    !< [K] top soil-node temperature (ground skin)
      real(wp),              intent(in)    :: leaf_temp(:)     !< [K] per-cohort leaf temperature
      real(wp),              intent(in)    :: wood_temp(:)     !< [K] per-cohort wood temperature
      type(aero_out_t),      intent(inout) :: aero
      real(wp),              intent(out)   :: tcas, qcas, press, rho, t_ground
      tcas = cas_enthalpy_to_temp(cas%can_enthalpy, cas%can_shv)
      qcas = cas%can_shv
      press = aenv%press
      rho = aenv%rho_air
      t_ground = soil_temp_top
      aenv%can_temp = tcas
      aenv%can_theta = tcas
      aenv%can_shv = qcas
      aenv%can_co2 = cas%can_co2
      aenv%t_ground = t_ground
      call aero_bottom_to_top(aero_cfg, aenv, ageom, col_cohort%n, col_cohort%height, col_cohort%lai,      &
                              col_cohort%leaf_width, col_cohort%branch_diam, leaf_temp, wood_temp, aero)
   end subroutine refresh_canopy_aerodynamics

   !---------------------------------------------------------------------------------------!
   ! root_zone_environment -- root-weighted soil temperature and the column-mean soil moisture     !
   ! that drive fine-root and heterotrophic respiration.                                             !
   !---------------------------------------------------------------------------------------!
   pure subroutine root_zone_environment(soil, root_frac, soil_temp, theta, soil_temp_root, theta_mean)
      type(soil_profile_t), intent(in)  :: soil
      real(wp),            intent(in)  :: root_frac(:)      !< [-]     the patch's fine roots by layer
      real(wp),            intent(in)  :: soil_temp(:)      !< [K]     per-layer soil temperature
      real(wp),            intent(in)  :: theta(:)          !< [m3/m3] per-layer soil moisture
      real(wp),            intent(out) :: soil_temp_root    !< [K]     root-fraction-weighted mean
      real(wp),            intent(out) :: theta_mean        !< [m3/m3] depth-weighted column mean
      integer(ik) :: k, nsl
      nsl = soil%n_active
      soil_temp_root = weighted_mean(soil_temp(1:nsl), root_frac, nsl)
      theta_mean = 0.0_wp
      do k = 1_ik, nsl
         theta_mean     = theta_mean     + theta(k) * soil%dz(k)
      end do
      theta_mean = theta_mean / max(-soil%soil_layer_z(nsl+1_ik), tiny_num)
   end subroutine root_zone_environment

   !---------------------------------------------------------------------------------------!
   ! canopy_leaf_gas_exchange -- per-cohort leaf gas exchange (GPP, gs, Rd) through the batch kernel,  !
   ! the patch GPP and leaf-respiration totals, and the frozen leaf-energy coefficients h_coeff_leaf   !
   ! and g_transp_leaf (docs/science/leaf_gas_exchange.md). psi_leaf is read once per dt_fast from the !
   ! state^n leaf water. The stomatal limb (eq. (10), §4.2) takes the cohort's predawn leaf potential, !
   ! dmax_psi_leaf, through the kernel's psi keyword; a cohort with no history yet (a recruit, a run's !
   ! first step) starts from its patch's top-layer soil potential rather than at fully turgid.        !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_leaf_gas_exchange(leaf_photo, hyd_table, soil, col_cohort,                         &
                                       forc, aero,                                                       &
                                       biophys, qcas, press, rho, gpp, ra_leaf, h_coeff_leaf, g_transp_leaf,      &
                                       gpp_coh, leaf_resp_coh, cdiag)
      type(leaf_photo_table_t),   intent(in)  :: leaf_photo   !< per-PFT leaf parameters (once per run)
      type(hydro_args_table_t),   intent(in)  :: hyd_table              !< per-PFT leaf PV curve
      type(soil_profile_t),       intent(in)  :: soil         !< surface-layer retention (dmax_psi seed)
      type(column_cohort_t),      intent(in)  :: col_cohort
      type(column_forcing_t),     intent(in)  :: forc         !< absorbed PAR per cohort
      type(aero_out_t),           intent(in)  :: aero         !< leaf boundary-layer conductances, wind
      type(patch_biophys_t),      intent(in)  :: biophys      !< leaf/wood temperature + water, CAS CO2, theta(1)
      real(wp),                   intent(in)  :: qcas, press, rho
      real(wp),                   intent(out) :: gpp, ra_leaf          !< [umol/m2 ground/s] patch totals
      real(wp),                   intent(out) :: h_coeff_leaf(:), g_transp_leaf(:)
      real(wp), optional,         intent(out) :: gpp_coh(:), leaf_resp_coh(:)   !< [umol/plant/s]
      real(wp), optional,         intent(inout) :: cdiag(:,:)

      !----- Bare-array batch I/O for the per-cohort physiology kernel. ----------------------------!
      real(wp) :: par_arr(col_cohort%n), vpd_arr(col_cohort%n), gb_arr(col_cohort%n), rho_mol_arr(col_cohort%n), &
            psi_leaf_arr(col_cohort%n)
      real(wp) :: dmax_psi_arr(col_cohort%n), dmax_psi_seed
      real(wp) :: a_gross_arr(col_cohort%n), gs_arr(col_cohort%n), rd_arr(col_cohort%n)
      real(wp) :: a_net_arr(col_cohort%n), ci_arr(col_cohort%n), cs_arr(col_cohort%n)
      real(wp) :: bstom_arr(col_cohort%n), bnstom_arr(col_cohort%n)
      real(wp) :: e_air, gsw_ms
      integer(ik) :: i, n

      n = col_cohort%n
      gpp = 0.0_wp
      ra_leaf = 0.0_wp
      if (present(gpp_coh))       gpp_coh(1:n)       = 0.0_wp
      if (present(leaf_resp_coh)) leaf_resp_coh(1:n) = 0.0_wp
      e_air = qcas * press / (0.622_wp + 0.378_wp * qcas)          ! loop-invariant
      do i = 1_ik, n
         rho_mol_arr(i)  = press / (r_gas * biophys%leaf_temp(i))
         par_arr(i)      = cohort_leaf_par(forc%abs_par(i), col_cohort%lai(i), forc%par_per_w)
         vpd_arr(i)      = max(sat_vapor_pressure(biophys%leaf_temp(i)) - e_air, 0.0_wp)
         gb_arr(i)       = aero%leaf_gbw(i) * rho_mol_arr(i)
         psi_leaf_arr(i) = leaf_node_psi(biophys%leaf_water_mass(i), col_cohort%bleaf(i),                 &
               biophys%wood_water_mass(i), col_cohort%bsap(i) + col_cohort%broot(i), col_cohort%height(i),   &
               hyd_table%pft(col_cohort%pft(i)))
      end do
      dmax_psi_seed = top_soil_psi(soil, biophys%soil_water%theta(1))
      do i = 1_ik, n
         if (is_unset_negative(col_cohort%dmax_psi_leaf(i))) then
            dmax_psi_arr(i) = dmax_psi_seed                    ! unset: seed from the soil
         else
            dmax_psi_arr(i) = col_cohort%dmax_psi_leaf(i)
         end if
      end do
      if (present(cdiag)) then
         call leaf_gas_exchange_batch(n, par_arr, biophys%leaf_temp(1:n), vpd_arr, biophys%cas%can_co2, press, &
                                      psi_leaf_arr, gb_arr, leaf_photo, col_cohort%pft(1:n),                  &
                                      col_cohort%vcmax25(1:n), col_cohort%rd25(1:n), a_gross_arr, gs_arr, rd_arr, &
                                      psi=dmax_psi_arr(1:n),                                            &
                                      a_net=a_net_arr, ci=ci_arr, cs=cs_arr,                            &
                                      beta_stom=bstom_arr, beta_nonstom=bnstom_arr)
         do i = 1_ik, n
            cdiag(CO_ANET,         i) = a_net_arr(i)
            cdiag(CO_AGROSS,       i) = a_gross_arr(i)
            cdiag(CO_GSW,          i) = gs_arr(i)
            cdiag(CO_GBW,          i) = aero%leaf_gbw(i)
            cdiag(CO_CI,           i) = ci_arr(i)
            cdiag(CO_CS,           i) = cs_arr(i)
            cdiag(CO_RD,           i) = rd_arr(i)
            cdiag(CO_BETA_STOM,    i) = bstom_arr(i)
            cdiag(CO_BETA_NONSTOM, i) = bnstom_arr(i)
            cdiag(CO_LEAF_TEMP,    i) = biophys%leaf_temp(i)
            cdiag(CO_WOOD_TEMP,    i) = biophys%wood_temp(i)
            cdiag(CO_LEAF_VPD,     i) = vpd_arr(i)
            cdiag(CO_PSI_LEAF,     i) = psi_leaf_arr(i)
            cdiag(CO_ABS_PAR,      i) = forc%abs_par(i)
            cdiag(CO_ABS_SW,       i) = forc%abs_sw(i)
            cdiag(CO_ABS_LW,       i) = forc%abs_lw(i)
            cdiag(CO_WIND,         i) = aero%wind(i)
            cdiag(CO_LEAF_WATER,   i) = biophys%leaf_water_mass(i)
            cdiag(CO_WOOD_WATER,   i) = biophys%wood_water_mass(i)
            cdiag(CO_LEAF_PAR,     i) = par_arr(i)
            cdiag(CO_GB_MOL,       i) = gb_arr(i)
            cdiag(CO_PSI_PREDAWN,  i) = dmax_psi_arr(i)
            cdiag(CO_CA,           i) = biophys%cas%can_co2
            cdiag(CO_PRESSURE,     i) = press
         end do
      else
         call leaf_gas_exchange_batch(n, par_arr, biophys%leaf_temp(1:n), vpd_arr, biophys%cas%can_co2, press, &
                                      psi_leaf_arr, gb_arr, leaf_photo, col_cohort%pft(1:n),                  &
                                      col_cohort%vcmax25(1:n), col_cohort%rd25(1:n), a_gross_arr, gs_arr, rd_arr, &
                                      psi=dmax_psi_arr(1:n))
      end if
      do i = 1_ik, n
         gsw_ms  = gs_arr(i) / max(rho_mol_arr(i), tiny_num)
         !----- GPP is gross A, and Rd is charged below either way, so the canopy air receives the !
         !      leaf's net exchange A_g - Rd whatever the water-stress controls did to A_g. -------!
         gpp     = gpp     + a_gross_arr(i) * col_cohort%leaf_area(i) * col_cohort%nplant(i)
         if (present(gpp_coh)) gpp_coh(i) = a_gross_arr(i) * col_cohort%leaf_area(i)
         if (present(cdiag))   cdiag(CO_GPP_RATE, i) = a_gross_arr(i) * col_cohort%leaf_area(i)
         ra_leaf = ra_leaf + rd_arr(i)      * col_cohort%leaf_area(i) * col_cohort%nplant(i)
         if (present(leaf_resp_coh)) leaf_resp_coh(i) = rd_arr(i) * col_cohort%leaf_area(i)
         h_coeff_leaf(i) = sensible_heat_coeff(EFFAREA_HEAT * col_cohort%lai(i), aero%leaf_gbh(i), rho, cp_air)
         g_transp_leaf(i)    = leaf_transp_coeff(EFFAREA_TRANSP, col_cohort%lai(i), aero%leaf_gbw(i), gsw_ms)
      end do
   end subroutine canopy_leaf_gas_exchange

   !---------------------------------------------------------------------------------------!
   ! A cohort's leaf PAR: the PAR its leaves absorb per unit of its OWN leaf area. The guard only   !
   ! keeps a leafless cohort, which absorbs nothing, at 0; a floor on the leaf area would spread a  !
   ! thin cohort's light over leaf area it does not have.                                           !
   !---------------------------------------------------------------------------------------!
   elemental real(wp) function cohort_leaf_par(abs_par, lai, par_per_w) result(par)
      real(wp), intent(in) :: abs_par    !< [W/m2 ground] PAR the cohort's leaves absorb
      real(wp), intent(in) :: lai        !< [m2 leaf/m2 ground] the cohort's leaf area index
      real(wp), intent(in) :: par_per_w  !< [umol/J] PAR photons per joule
      par = abs_par / max(lai, tiny_num) * par_per_w
   end function cohort_leaf_par

   !---------------------------------------------------------------------------------------!
   ! canopy_maintenance_respiration -- stem and fine-root maintenance respiration, per cohort       !
   ! [umol/plant/s] and as patch totals [umol/m2 ground/s] (docs/science/plant_respiration.md       !
   ! eqs. (3)-(4)), both on leaf Rd's temperature response. Elemental kernels: the array actuals    !
   ! drive the broadcast; the response and the patch-uniform soil temperature scale broadcast.      !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_maintenance_respiration(ea_rd, hd_rd, ds_rd, soil, root_frac, col_cohort, wood_temp, soil_temp, &
                                             ra_stem, ra_root, stem_resp_coh, root_resp_coh)
      real(wp),              intent(in)  :: ea_rd, hd_rd, ds_rd !< leaf Rd's peaked response (stem and roots)
      real(wp),              intent(in)  :: root_frac(:)        !< [-] the patch's fine roots by layer
      type(soil_profile_t),  intent(in)  :: soil                !< active layer count
      type(column_cohort_t), intent(in)  :: col_cohort
      real(wp),              intent(in)  :: wood_temp(:)        !< [K] per-cohort wood temperature
      real(wp),              intent(in)  :: soil_temp(:)        !< [K] PER-LAYER soil temperature
      real(wp),              intent(out) :: ra_stem, ra_root    !< [umol/m2 ground/s]
      real(wp), optional,    intent(out) :: stem_resp_coh(:), root_resp_coh(:)   !< [umol/plant/s]
      real(wp) :: stem_resp_arr(col_cohort%n), root_resp_arr(col_cohort%n), tscale_root
      integer(ik) :: i, n
      n = col_cohort%n
      !----- The temperature response, summed over layers (docs/science/plant_respiration.md eq. (4)), !
      !      is patch-uniform: evaluated once and broadcast over the cohorts. ---------------------!
      tscale_root = root_zone_temp_scale(soil_temp, root_frac, soil%n_active, ea_rd, hd_rd, ds_rd)
      ra_stem = 0.0_wp
      ra_root = 0.0_wp
      if (present(stem_resp_coh)) stem_resp_coh(1:n) = 0.0_wp
      if (present(root_resp_coh)) root_resp_coh(1:n) = 0.0_wp
      call stem_maintenance_respiration(wood_temp(1:n), col_cohort%dbh(1:n), col_cohort%height(1:n),   &
                                   col_cohort%wai(1:n), col_cohort%nplant(1:n),                          &
                                   col_cohort%aboveground_frac(1:n), col_cohort%is_woody(1:n),          &
                                   col_cohort%stem_resp_factor25(1:n), col_cohort%sap_volume(1:n),      &
                                   col_cohort%stem_resp_sapwood25(1:n), ea_rd, hd_rd, ds_rd,          &
                                   stem_resp_arr(1:n))
      call fine_root_maintenance_respiration(tscale_root, col_cohort%broot(1:n),                        &
                                   col_cohort%root_resp_factor25(1:n), root_resp_arr(1:n))
      do i = 1_ik, n
         ra_stem = ra_stem + stem_resp_arr(i) * col_cohort%nplant(i)
         ra_root = ra_root + root_resp_arr(i) * col_cohort%nplant(i)
         if (present(stem_resp_coh)) stem_resp_coh(i) = stem_resp_arr(i)
         if (present(root_resp_coh)) root_resp_coh(i) = root_resp_arr(i)
      end do
   end subroutine canopy_maintenance_respiration

   !---------------------------------------------------------------------------------------!
   ! patch_heterotrophic_respiration -- Rh [umol/m2/s] from the CENTURY matrix on the patch's frozen  !
   ! pools (docs/science/soil_carbon.md §7, eq. (2)), so the day's fast Rh equals the daily step's    !
   ! pool debit by construction; xi goes to budget%xi_step for the caller to accumulate. Zero with     !
   ! soil carbon off.                                                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine patch_heterotrophic_respiration(col_config, col_site, soil_carbon, t_ground, soil_temp_root,  &
                                              theta_mean, rh, budget)
      type(column_config_t), intent(in)    :: col_config
      type(column_site_t),   intent(in)    :: col_site
      type(soil_carbon_t),   intent(in)    :: soil_carbon
      real(wp),              intent(in)    :: t_ground, soil_temp_root, theta_mean
      real(wp),              intent(out)   :: rh              !< [umol/m2/s]
      type(column_budget_t), intent(inout) :: budget
      real(wp) :: xi(N_SOIL_POOL), a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL)
      real(wp) :: rh_kgc                                     !< [kgC/m2/day] the matrix's own currency
      if (col_config%soil_carbon%soil_carbon_on) then
         call assemble_env_scalar(t_ground, soil_temp_root, theta_mean, col_site%soil_profile%theta_res(1),        &
                                  col_site%soil_profile%theta_sat(1), soil_carbon, col_config%soil_carbon, xi)
         call assemble_transfer_matrix(soil_carbon, col_config%soil_carbon, a_mat, k_diag, er)
         !----- The matrix returns kgC/m2/day, the pools' currency, which budget%rh_matrix_step keeps for  !
         !      the seam check; the canopy air gets umol/m2/s (docs/science/soil_carbon.md §7). ----------!
         rh_kgc = heterotrophic_respiration_matrix(a_mat, k_diag, xi, soil_carbon)
         budget%xi_step = xi
         budget%rh_matrix_step = rh_kgc      ! [kgC/m2/day] -- the seam's currency
         rh = rh_kgc * kgcday_2_umols                              ! [umol/m2/s]  -- what the CAS gets
      else
         !----- No soil carbon, no soil respiration: the pools are empty and take no litter, so any Rh   !
         !      here would be carbon from nothing (docs/science/soil_carbon.md §Parameters). ----------!
         rh = 0.0_wp
      end if
   end subroutine patch_heterotrophic_respiration

   !---------------------------------------------------------------------------------------!
   ! cas_capacities_and_conductances -- canopy-air-space mass and molar capacities [kg/m2],         !
   ! [mol/m2] and the CAS <-> atmosphere conductances for enthalpy, vapour and CO2                  !
   ! (docs/science/canopy_aerodynamics.md eq. (5)). Heat rides the temp1 profile factor, vapour and !
   ! CO2 ride temp2 (equal in canopy_aerodynamics, z0q = z0h).                                      !
   !---------------------------------------------------------------------------------------!
   pure subroutine cas_capacities_and_conductances(rho, qcas, can_depth, ustar, temp1, temp2,          &
                                                   cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2)
      real(wp), intent(in)  :: rho, qcas, can_depth, ustar, temp1, temp2
      real(wp), intent(out) :: cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2
      real(wp) :: can_dmol
      can_dmol = cas_molar_density(rho, qcas)
      cas_mass_capacity = rho      * can_depth
      cas_molar_capacity = can_dmol * can_depth
      call cas_atm_conductances(rho, can_dmol, ustar, temp1, temp2, g_atm_heat, g_atm_vapour, g_atm_co2)
   end subroutine cas_capacities_and_conductances

   !---------------------------------------------------------------------------------------!
   ! aero_bottom_to_top -- canopy_aerodynamics takes its cohorts bottom (1) to top (n); the cohort     !
   ! arrays are tallest first. Only the wind cascade and the per-cohort boundary layers depend on     !
   ! order, so the inputs are permuted by ord (ord(k) = the k-th cohort from the bottom,              !
   ! ascending_order) and the per-cohort outputs scattered back. An identity for n <= 1.              !
   !---------------------------------------------------------------------------------------!
   subroutine aero_bottom_to_top(acfg, aenv, ageom, n, height, lai, leaf_width, branch_diam,           &
                                 leaf_temp, wood_temp, aero)
      type(aero_cfg_t),      intent(in)    :: acfg
      type(aero_env_t),      intent(in)    :: aenv
      type(aero_geom_t),     intent(in)    :: ageom
      integer(ik),           intent(in)    :: n
      real(wp),              intent(in)    :: height(:)       !< [m]  cohort height (gather order, top first)
      real(wp),              intent(in)    :: lai(:)          !< [m2/m2] leaf area index
      real(wp),              intent(in)    :: leaf_width(:)   !< [m]  leaf width (boundary layer)
      real(wp),              intent(in)    :: branch_diam(:)  !< [m]  branch diameter (wood boundary layer)
      real(wp),              intent(in)    :: leaf_temp(:)    !< [K]  leaf temperature
      !----- The wood's own temperature drives its boundary layer's free convection; it lags the leaves. !
      real(wp),              intent(in)    :: wood_temp(:)    !< [K]  wood temperature
      type(aero_out_t),      intent(inout) :: aero
      integer(ik) :: ord(n), k, i
      real(wp)    :: h_bt(n), lai_bt(n), lt_bt(n), wt_bt(n), lw_bt(n), bd_bt(n)
      real(wp)    :: wind_bt(n), lgbh_bt(n), lgbw_bt(n), wgbh_bt(n), wgbw_bt(n)

      call ascending_order(height, n, ord)
      do k = 1_ik, n
         i = ord(k)
         h_bt(k)   = height(i)
         lai_bt(k) = lai(i)
         lt_bt(k)  = leaf_temp(i)
         wt_bt(k)  = wood_temp(i)
         lw_bt(k)  = leaf_width(i)
         bd_bt(k)  = branch_diam(i)
      end do

      call canopy_aerodynamics(acfg, aenv, ageom, n, h_bt, lai_bt, lt_bt, wt_bt, lw_bt, bd_bt, aero)

      !----- aero%*(k) is now bottom->top; copy out, then scatter back to gather order. ----------!
      do k = 1_ik, n
         wind_bt(k) = aero%wind(k)
         lgbh_bt(k) = aero%leaf_gbh(k)
         lgbw_bt(k) = aero%leaf_gbw(k)
         wgbh_bt(k) = aero%wood_gbh(k)
         wgbw_bt(k) = aero%wood_gbw(k)
      end do
      do k = 1_ik, n
         aero%wind(ord(k))     = wind_bt(k)
         aero%leaf_gbh(ord(k)) = lgbh_bt(k)
         aero%leaf_gbw(ord(k)) = lgbw_bt(k)
         aero%wood_gbh(ord(k)) = wgbh_bt(k)
         aero%wood_gbw(ord(k)) = wgbw_bt(k)
      end do
   end subroutine aero_bottom_to_top

end module meds_fast_prepass
