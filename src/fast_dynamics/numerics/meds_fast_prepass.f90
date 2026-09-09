!==========================================================================================!
! meds_fast_prepass -- the once-per-dt_fast PRE-PASS of the fast biophysics loop, split into the  !
! five processes it fuses (2026-09 review, item 4 #4): canopy aerodynamics from the current       !
! canopy-air state, leaf gas exchange (GPP / stomatal conductance / dark respiration) with the      !
! frozen leaf-energy coefficients, stem + fine-root maintenance respiration, heterotrophic          !
! respiration, and the canopy-air-space capacities + atmosphere-exchange conductances.              !
! `column_prepass` is the thin orchestrator that runs them in order and assembles the biotic CO2   !
! flux; both integrators reach it through build_column_frozen (meds_fast_ark), so it is the ONE     !
! authority that keeps ARK and RK45 GPP bit-for-bit. ED2 freezes gs and hydraulics per DTLSM in the !
! same way.                                                                                          !
!                                                                                                    !
! The per-cohort arithmetic keeps the SAME i = 1..n accumulation order the fused routine used, and   !
! each patch total is summed independently of the others, so splitting the one loop into a leaf     !
! loop and a respiration loop leaves every output bit-identical.                                     !
!==========================================================================================!
module meds_fast_prepass
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, cp_air, r_gas, grav_head
   use meds_config,           only : meds_config_t
   use meds_site_diag_types,  only : CD_ANET, CD_AGROSS, CD_GSW, CD_GBW, CD_CI, CD_CS, CD_RD,      &
                                     CD_TRANSP, CD_BETA_STOM, CD_BETA_NONSTOM, CD_LEAF_TEMP,       &
                                     CD_WOOD_TEMP, CD_LEAF_VPD, CD_PSI_LEAF, CD_ABS_PAR, CD_ABS_SW, &
                                     CD_ABS_LW, CD_WIND, CD_LEAF_WATER, CD_WOOD_WATER, CD_GPP_RATE
   use meds_hydr_lib,         only : soil_psi_from_theta, psi_from_water_content
   use meds_biophysics_types, only : aero_env_t, aero_geom_t, aero_out_t, patch_biophys_t, veg_thermal_params_t
   use meds_column_params, only : soil_params_t
   use meds_biophysics_opts, only : aero_cfg_t
   use meds_column_reservoirs, only : cas_state_t, soil_carbon_t
   use meds_plant_types, only : leaf_photo_table_t, hydro_params_t, wood_params_t, root_params_t
   use meds_biogeochem_types, only : co2_opts_t, n_soil_pool
   use meds_fast_types,       only : column_config_t, column_cohort_t, column_forcing_t, column_budget_t
   use meds_canopy_aerodynamics, only : canopy_aerodynamics, cas_atm_conductances
   use meds_vegetation_biophysics, only : sensible_heat_coeff, leaf_transp_coeff
   use meds_leaf_gas_exchange, only : leaf_gas_exchange_batch
   use meds_plant_respiration, only : stem_maintenance_respiration, fine_root_maintenance_respiration
   use meds_soil_biogeochem,  only : heterotrophic_respiration_flux, heterotrophic_respiration_matrix, &
                                     assemble_env_scalar, assemble_transfer_matrix
   use meds_therm_lib,        only : cas_molar_density, cas_temp_of_enthalpy, sat_vapor_pressure
   use meds_numerics,         only : weighted_mean
   implicit none
   private

   public :: column_prepass
   public :: refresh_canopy_aerodynamics, canopy_leaf_gas_exchange, canopy_maintenance_respiration
   public :: root_zone_environment, patch_heterotrophic_respiration, cas_capacities_and_conductances
   public :: aero_bottom_to_top

contains

   !---------------------------------------------------------------------------------------!
   ! column_prepass -- run the five processes in order and assemble the biotic CO2 flux.          !
   !   NEE_biotic = leaf Rd + stem + root maintenance respiration + heterotrophic Rh - GPP          !
   ! [umol/m2/s], positive to the canopy air. `biophys` is intent(in): callers that need the CAS     !
   ! temperature persisted write biophys%cas%can_temp = tcas themselves right after the call.         !
   !---------------------------------------------------------------------------------------!
   subroutine column_prepass(cfg, col_config, aenv, ageom, col_cohort, forc, biophys, aero, budget,   &
                             tcas, qcas, press, rho, t_ground, h_coeff_f, g_tr_f,                     &
                             cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2, nee_biotic, &
                             gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, cdiag)
      type(meds_config_t),     intent(in)    :: cfg
      type(column_config_t),   intent(in)    :: col_config
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: col_cohort
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(in)    :: biophys
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budget
      real(wp),                intent(out)   :: tcas, qcas, press, rho, t_ground
      real(wp),                intent(out)   :: h_coeff_f(:), g_tr_f(:)
      real(wp),                intent(out)   :: cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2, &
           nee_biotic
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      !----- OPTIONAL per-cohort DIAGNOSTIC capture (MEDS_IO_V01_PLAN.md section 3.4). Present only    !
      !      when the run reports per-cohort ecophysiology; absent, the extra leaf_flux_t fields are    !
      !      never even requested from the batch kernel, so this costs nothing.  ---------------------!
      real(wp), optional,      intent(inout) :: cdiag(:,:)   !< (N_CDIAG, ncoh) INSTANTANEOUS values

      real(wp) :: gpp, ra_leaf, ra_stem, ra_root, rh, soil_temp_root, theta_mean

      !----- 1. aerodynamics from the current canopy-air state. ------------------------------------!
      call refresh_canopy_aerodynamics(col_config%aero, aenv, ageom, col_cohort, biophys%cas,           &
                                       biophys%soil_e%soil_temp(1), biophys%leaf_temp, aero,           &
                                       tcas, qcas, press, rho, t_ground)

      !----- 2. root-zone environment (root + heterotrophic respiration drivers). ------------------!
      call root_zone_environment(col_config%soil, biophys%soil_e%soil_temp, biophys%soil_w%theta,       &
                                 soil_temp_root, theta_mean)

      !----- 3. leaf gas exchange + the frozen leaf-energy coefficients. ---------------------------!
      call canopy_leaf_gas_exchange(col_config%leaf_photo, col_config%hydro_p, col_config%veg_thermal,  &
                                    col_config%soil, col_cohort, forc, aero, biophys,                   &
                                    qcas, press, rho, gpp, ra_leaf, h_coeff_f, g_tr_f,                  &
                                    gpp_coh, leaf_resp_coh, cdiag)

      !----- 4. stem + fine-root maintenance respiration. -----------------------------------------!
      call canopy_maintenance_respiration(col_config%wood, col_config%root, col_cohort,                 &
                                          biophys%wood_temp, soil_temp_root, ra_stem, ra_root,          &
                                          stem_resp_coh, root_resp_coh)

      !----- 5. heterotrophic respiration. -------------------------------------------------------!
      call patch_heterotrophic_respiration(cfg, col_config, biophys%soil_carbon, t_ground,              &
                                           soil_temp_root, theta_mean, rh, budget)

      !----- NEE = autotrophic (leaf Rd + stem + root) + heterotrophic Rh - GPP. ------------------!
      nee_biotic = ra_leaf + ra_stem + ra_root + rh - gpp
      budget%gpp_last = gpp ; budget%nee_last = nee_biotic

      !----- 6. CAS capacities + atm-exchange conductances (frozen across the macro-step). --------!
      call cas_capacities_and_conductances(rho, qcas, biophys%cas%can_depth, aero%ustar, aero%temp1,    &
                                           aero%temp2, cas_mass_capacity, cas_molar_capacity, g_atm_heat, g_atm_vapour, g_atm_co2)
   end subroutine column_prepass

   !---------------------------------------------------------------------------------------!
   ! refresh_canopy_aerodynamics -- the canopy-air scalars diagnosed from the CAS state, written    !
   ! into the aerodynamic environment, and the whole-canopy aerodynamics solve (friction velocity,   !
   ! scalar profile factors, ground conductance, per-cohort boundary layers).                        !
   !---------------------------------------------------------------------------------------!
   subroutine refresh_canopy_aerodynamics(aero_cfg, aenv, ageom, col_cohort, cas, soil_temp_top,      &
                                          leaf_temp, aero, tcas, qcas, press, rho, t_ground)
      type(aero_cfg_t),      intent(in)    :: aero_cfg
      type(aero_env_t),      intent(inout) :: aenv
      type(aero_geom_t),     intent(in)    :: ageom
      type(column_cohort_t), intent(in)    :: col_cohort
      type(cas_state_t),     intent(in)    :: cas
      real(wp),              intent(in)    :: soil_temp_top    !< [K] top soil-node temperature (ground skin)
      real(wp),              intent(in)    :: leaf_temp(:)     !< [K] per-cohort leaf temperature
      type(aero_out_t),      intent(inout) :: aero
      real(wp),              intent(out)   :: tcas, qcas, press, rho, t_ground
      tcas = cas_temp_of_enthalpy(cas%can_enthalpy, cas%can_shv)
      qcas = cas%can_shv ; press = aenv%press ; rho = aenv%rho_air ; t_ground = soil_temp_top
      aenv%can_temp = tcas ; aenv%can_theta = tcas ; aenv%can_shv = qcas ; aenv%can_co2 = cas%can_co2
      aenv%t_ground = t_ground
      call aero_bottom_to_top(aero_cfg, aenv, ageom, col_cohort%n, col_cohort%height, col_cohort%lai,      &
                              col_cohort%crown, col_cohort%leaf_width, col_cohort%branch_diam, leaf_temp, aero)
   end subroutine refresh_canopy_aerodynamics

   !---------------------------------------------------------------------------------------!
   ! root_zone_environment -- root-weighted soil temperature and the column-mean soil moisture     !
   ! that drive fine-root and heterotrophic respiration.                                             !
   !---------------------------------------------------------------------------------------!
   pure subroutine root_zone_environment(soil, soil_temp, theta, soil_temp_root, theta_mean)
      type(soil_params_t), intent(in)  :: soil
      real(wp),            intent(in)  :: soil_temp(:)      !< [K]     per-layer soil temperature
      real(wp),            intent(in)  :: theta(:)          !< [m3/m3] per-layer soil moisture
      real(wp),            intent(out) :: soil_temp_root    !< [K]     root-fraction-weighted mean
      real(wp),            intent(out) :: theta_mean        !< [m3/m3] depth-weighted column mean
      integer(ik) :: k, nsl
      nsl = soil%n_active
      soil_temp_root = weighted_mean(soil_temp(1:nsl), soil%root_frac, nsl)
      theta_mean = 0.0_wp
      do k = 1_ik, nsl
         theta_mean     = theta_mean     + theta(k) * soil%dz(k)
      end do
      theta_mean = theta_mean / max(-soil%soil_layer_z(nsl+1_ik), tiny_num)
   end subroutine root_zone_environment

   !---------------------------------------------------------------------------------------!
   ! canopy_leaf_gas_exchange -- per-cohort leaf gas exchange (GPP / gs / Rd) over the whole      !
   ! patch at once through the bare-array batch kernel, the patch GPP and leaf-respiration totals, !
   ! and the FROZEN leaf-energy coefficients h_coeff_f / g_tr_f that the stage kernels consume.     !
   !                                                                                                !
   ! psi_leaf for gs stays FROZEN (Category-0, ED2-faithful): diagnosed ONCE per dt_fast from the   !
   ! prognostic leaf_water_mass^n -- never refreshed per stage.                                     !
   !                                                                                                !
   ! STOMATAL WATER STRESS (issue #95): beta_stomata = min(1, exp(sref*psi)) is driven by           !
   ! YESTERDAY's daily-maximum leaf water potential -- the model's predawn potential.               !
   ! DMAX_PSI_LEAF_UNSET (positive, so unmistakable) means the cohort has no history yet: a recruit, !
   ! or the first step of a run. It is seeded from the SURFACE-LAYER soil potential so it starts at  !
   ! its patch's actual water status rather than at 0, which would read as fully turgid.             !
   !                                                                                                !
   ! NOTE ON THE `psi=` KEYWORD: what is passed is dmax_psi_leaf -- the cohort's own predawn LEAF   !
   ! potential, NOT a soil potential. The leaf kernel's dummy is called psi because the Sabot        !
   ! stomatal limb is conventionally keyed on soil/predawn potential, and predawn leaf psi IS the    !
   ! plant's overnight equilibration with the soil in WET soil. They do NOT coincide under drought,  !
   ! which is exactly the regime this feedback exists for. Renaming the kernel dummy is deferred     !
   ! because `psi` is a published Python keyword (meds.plant.leaf) -- see issue #99.                 !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_leaf_gas_exchange(leaf_photo, hydro_p, veg_thermal, soil, col_cohort, forc, aero,  &
                                       biophys, qcas, press, rho, gpp, ra_leaf, h_coeff_f, g_tr_f,      &
                                       gpp_coh, leaf_resp_coh, cdiag)
      type(leaf_photo_table_t),   intent(in)  :: leaf_photo   !< per-PFT leaf parameters (once per run)
      type(hydro_params_t),       intent(in)  :: hydro_p      !< leaf PV curve (psi_leaf from water content)
      type(veg_thermal_params_t), intent(in)  :: veg_thermal  !< effective exchange areas
      type(soil_params_t),        intent(in)  :: soil         !< surface-layer retention (dmax_psi seed)
      type(column_cohort_t),      intent(in)  :: col_cohort
      type(column_forcing_t),     intent(in)  :: forc         !< absorbed PAR per cohort
      type(aero_out_t),           intent(in)  :: aero         !< leaf boundary-layer conductances, wind
      type(patch_biophys_t),      intent(in)  :: biophys      !< leaf/wood temperature + water, CAS CO2, theta(1)
      real(wp),                   intent(in)  :: qcas, press, rho
      real(wp),                   intent(out) :: gpp, ra_leaf          !< [umol/m2 ground/s] patch totals
      real(wp),                   intent(out) :: h_coeff_f(:), g_tr_f(:)
      real(wp), optional,         intent(out) :: gpp_coh(:), leaf_resp_coh(:)   !< [umol/plant/s]
      real(wp), optional,         intent(inout) :: cdiag(:,:)

      !----- Bare-array batch I/O for the per-cohort physiology kernel (MEDS_NUMERICS_SCOPING.md).   !
      real(wp) :: par_arr(col_cohort%n), vpd_arr(col_cohort%n), gb_arr(col_cohort%n), rho_mol_arr(col_cohort%n), &
            psi_leaf_arr(col_cohort%n)
      real(wp) :: dmax_psi_arr(col_cohort%n), dmax_psi_seed
      real(wp) :: a_gross_arr(col_cohort%n), gs_arr(col_cohort%n), rd_arr(col_cohort%n)
      real(wp) :: a_net_arr(col_cohort%n), ci_arr(col_cohort%n), cs_arr(col_cohort%n), transp_arr(col_cohort%n)
      real(wp) :: bstom_arr(col_cohort%n), bnstom_arr(col_cohort%n)
      real(wp) :: e_air, gsw_ms
      integer(ik) :: i, n

      n = col_cohort%n
      gpp = 0.0_wp ; ra_leaf = 0.0_wp
      if (present(gpp_coh))       gpp_coh(1:n)       = 0.0_wp
      if (present(leaf_resp_coh)) leaf_resp_coh(1:n) = 0.0_wp
      e_air = qcas * press / (0.622_wp + 0.378_wp * qcas)          ! loop-invariant
      do i = 1_ik, n
         rho_mol_arr(i)  = press / (r_gas * biophys%leaf_temp(i))
         par_arr(i)      = forc%abs_par(i) / max(col_cohort%lai(i), 0.1_wp) * forc%par_per_w
         vpd_arr(i)      = max(sat_vapor_pressure(biophys%leaf_temp(i)) - e_air, 0.0_wp)
         gb_arr(i)       = aero%leaf_gbw(i) * rho_mol_arr(i)
         psi_leaf_arr(i) = psi_from_water_content(biophys%leaf_water_mass(i), hydro_p%leaf_pi0,      &
              hydro_p%leaf_elastic_mod, hydro_p%leaf_apoplast_frac,                                  &
              hydro_p%leaf_water_sat, col_cohort%bleaf(i))
      end do
      dmax_psi_seed = grav_head * soil_psi_from_theta(soil%retention, biophys%soil_w%theta(1),        &
                    soil%theta_sat(1), soil%theta_res(1), soil%vg_alpha(1), soil%vg_n(1))
      do i = 1_ik, n
         if (col_cohort%dmax_psi_leaf(i) > 0.0_wp) then
            dmax_psi_arr(i) = dmax_psi_seed                    ! UNSET sentinel -> seed from the soil
         else
            dmax_psi_arr(i) = col_cohort%dmax_psi_leaf(i)
         end if
      end do
      if (present(cdiag)) then
         call leaf_gas_exchange_batch(n, par_arr, biophys%leaf_temp(1:n), vpd_arr, biophys%cas%can_co2, press, &
                                      psi_leaf_arr, gb_arr, leaf_photo, col_cohort%pft(1:n),                  &
                                      col_cohort%vcmax25(1:n), col_cohort%rd25(1:n), a_gross_arr, gs_arr, rd_arr, &
                                      psi=dmax_psi_arr(1:n),                                            &
                                      a_net=a_net_arr, ci=ci_arr, cs=cs_arr, transp=transp_arr,         &
                                      beta_stom=bstom_arr, beta_nonstom=bnstom_arr)
         do i = 1_ik, n
            cdiag(CD_ANET,         i) = a_net_arr(i)
            cdiag(CD_AGROSS,       i) = a_gross_arr(i)
            cdiag(CD_GSW,          i) = gs_arr(i)
            cdiag(CD_GBW,          i) = aero%leaf_gbw(i)
            cdiag(CD_CI,           i) = ci_arr(i)
            cdiag(CD_CS,           i) = cs_arr(i)
            cdiag(CD_RD,           i) = rd_arr(i)
            cdiag(CD_TRANSP,       i) = transp_arr(i)
            cdiag(CD_BETA_STOM,    i) = bstom_arr(i)
            cdiag(CD_BETA_NONSTOM, i) = bnstom_arr(i)
            cdiag(CD_LEAF_TEMP,    i) = biophys%leaf_temp(i)
            cdiag(CD_WOOD_TEMP,    i) = biophys%wood_temp(i)
            cdiag(CD_LEAF_VPD,     i) = vpd_arr(i)
            cdiag(CD_PSI_LEAF,     i) = psi_leaf_arr(i)
            cdiag(CD_ABS_PAR,      i) = forc%abs_par(i)
            cdiag(CD_ABS_SW,       i) = forc%abs_sw(i)
            cdiag(CD_ABS_LW,       i) = forc%abs_lw(i)
            cdiag(CD_WIND,         i) = aero%wind(i)
            cdiag(CD_LEAF_WATER,   i) = biophys%leaf_water_mass(i)
            cdiag(CD_WOOD_WATER,   i) = biophys%wood_water_mass(i)
         end do
      else
         call leaf_gas_exchange_batch(n, par_arr, biophys%leaf_temp(1:n), vpd_arr, biophys%cas%can_co2, press, &
                                      psi_leaf_arr, gb_arr, leaf_photo, col_cohort%pft(1:n),                  &
                                      col_cohort%vcmax25(1:n), col_cohort%rd25(1:n), a_gross_arr, gs_arr, rd_arr, &
                                      psi=dmax_psi_arr(1:n))
      end if
      do i = 1_ik, n
         gsw_ms  = gs_arr(i) / max(rho_mol_arr(i), tiny_num)
         gpp     = gpp     + a_gross_arr(i) * col_cohort%leaf_area(i) * col_cohort%nplant(i)
         if (present(gpp_coh)) gpp_coh(i) = a_gross_arr(i) * col_cohort%leaf_area(i)
         if (present(cdiag))   cdiag(CD_GPP_RATE, i) = a_gross_arr(i) * col_cohort%leaf_area(i)
         ra_leaf = ra_leaf + rd_arr(i)      * col_cohort%leaf_area(i) * col_cohort%nplant(i)
         if (present(leaf_resp_coh)) leaf_resp_coh(i) = rd_arr(i) * col_cohort%leaf_area(i)
         h_coeff_f(i) = sensible_heat_coeff(veg_thermal%effarea_heat * col_cohort%lai(i), aero%leaf_gbh(i), rho, cp_air)
         g_tr_f(i)    = leaf_transp_coeff(veg_thermal%effarea_transp, col_cohort%lai(i), aero%leaf_gbw(i), gsw_ms)
      end do
   end subroutine canopy_leaf_gas_exchange

   !---------------------------------------------------------------------------------------!
   ! canopy_maintenance_respiration -- stem and fine-root maintenance respiration, per cohort     !
   ! [umol/plant/s] and as patch totals [umol/m2 ground/s]. Elemental kernels: the array actuals    !
   ! drive the broadcast; the parameter records and the patch-uniform soil temperature broadcast.   !
   !---------------------------------------------------------------------------------------!
   subroutine canopy_maintenance_respiration(wood, root, col_cohort, wood_temp, soil_temp_root,        &
                                             ra_stem, ra_root, stem_resp_coh, root_resp_coh)
      type(wood_params_t),   intent(in)  :: wood
      type(root_params_t),   intent(in)  :: root
      type(column_cohort_t), intent(in)  :: col_cohort
      real(wp),              intent(in)  :: wood_temp(:)        !< [K] per-cohort wood temperature
      real(wp),              intent(in)  :: soil_temp_root      !< [K] root-weighted soil temperature
      real(wp),              intent(out) :: ra_stem, ra_root    !< [umol/m2 ground/s]
      real(wp), optional,    intent(out) :: stem_resp_coh(:), root_resp_coh(:)   !< [umol/plant/s]
      real(wp) :: stem_resp_arr(col_cohort%n), root_resp_arr(col_cohort%n)
      integer(ik) :: i, n
      n = col_cohort%n
      ra_stem = 0.0_wp ; ra_root = 0.0_wp
      if (present(stem_resp_coh)) stem_resp_coh(1:n) = 0.0_wp
      if (present(root_resp_coh)) root_resp_coh(1:n) = 0.0_wp
      call stem_maintenance_respiration(wood_temp(1:n), col_cohort%dbh(1:n), col_cohort%height(1:n),   &
                                   col_cohort%wai(1:n), col_cohort%nplant(1:n), wood, stem_resp_arr(1:n))
      call fine_root_maintenance_respiration(soil_temp_root, col_cohort%broot(1:n), root, root_resp_arr(1:n))
      do i = 1_ik, n
         ra_stem = ra_stem + stem_resp_arr(i) * col_cohort%nplant(i)
         ra_root = ra_root + root_resp_arr(i) * col_cohort%nplant(i)
         if (present(stem_resp_coh)) stem_resp_coh(i) = stem_resp_arr(i)
         if (present(root_resp_coh)) root_resp_coh(i) = root_resp_arr(i)
      end do
   end subroutine canopy_maintenance_respiration

   !---------------------------------------------------------------------------------------!
   ! patch_heterotrophic_respiration -- Rh [umol/m2/s] from EITHER the constant-pool scalar form   !
   ! (soil_carbon_on = .false.) OR the CENTURY matrix over the FROZEN per-patch pool held constant   !
   ! across today's sub-steps (B2, MEDS_SLOW_DYNAMICS_DESIGN.md Part II section 9): the day's total  !
   ! fast Rh then equals the daily soil_carbon_step's pool debit BY CONSTRUCTION, since both read    !
   ! the same frozen pool + the same per-pool env scalar xi (written to budget%xi_step for the caller !
   ! to integrate into xi_int).                                                                       !
   !---------------------------------------------------------------------------------------!
   subroutine patch_heterotrophic_respiration(cfg, col_config, soil_carbon, t_ground, soil_temp_root,  &
                                              theta_mean, rh, budget)
      type(meds_config_t),   intent(in)    :: cfg
      type(column_config_t), intent(in)    :: col_config
      type(soil_carbon_t),   intent(in)    :: soil_carbon
      real(wp),              intent(in)    :: t_ground, soil_temp_root, theta_mean
      real(wp),              intent(out)   :: rh              !< [umol/m2/s]
      type(column_budget_t), intent(inout) :: budget
      real(wp) :: xi(n_soil_pool), a_mat(n_soil_pool, n_soil_pool), k_diag(n_soil_pool), er(n_soil_pool)
      if (cfg%soil_carbon_on) then
         call assemble_env_scalar(t_ground, soil_temp_root, theta_mean, col_config%soil%theta_res(1),      &
                                  col_config%soil%theta_sat(1), soil_carbon, cfg%soil_carbon, xi)
         call assemble_transfer_matrix(soil_carbon, cfg%soil_carbon, a_mat, k_diag, er)
         rh = heterotrophic_respiration_matrix(a_mat, k_diag, xi, soil_carbon)
         budget%xi_step = xi ; budget%rh_matrix_step = rh
      else
         rh = heterotrophic_respiration_flux(col_config%fast_soil_carbon, soil_temp_root, theta_mean,      &
                                             col_config%soil%theta_res(1), col_config%soil%theta_sat(1), col_config%co2)
      end if
   end subroutine patch_heterotrophic_respiration

   !---------------------------------------------------------------------------------------!
   ! cas_capacities_and_conductances -- canopy-air-space mass and molar capacities [kg/m2],        !
   ! [mol/m2] and the CAS <-> atmosphere conductances for enthalpy, vapour and CO2. Heat rides the  !
   ! temp1 profile factor, vapour and CO2 ride temp2 (equal in canopy_aerodynamics, z0q = z0h).     !
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
   ! aero_bottom_to_top -- solve canopy aerodynamics with the cohort order it CONTRACTS for --      !
   ! BOTTOM(1)->TOP(n) -- from height-DESCENDING per-cohort arrays. Only the wind cascade + the       !
   ! per-cohort boundary layers depend on order; the whole-canopy scalars (ustar/temp1/temp2/uh) do !
   ! not. An ascending-height permutation `ord` reverses the per-cohort inputs; the per-cohort wind  !
   ! and leaf/wood conductance outputs are scattered back to gather order. Identity for n<=1, so     !
   ! single-cohort behaviour is bit-unchanged.                                                       !
   !                                                                                                 !
   ! ord(k) = gather index of the k-th cohort counting from the canopy BOTTOM. This used to be an   !
   ! O(n^2) selection sort -- the ONLY superlinear term in the whole fast loop (0.97 s at n = 2000,  !
   ! 16.5 s at n = 8000 over a 6-day run) -- and it is REDUNDANT in the normal case: sort_cohorts    !
   ! leaves the cohort block height-DESCENDING and the gather preserves that order, so bottom-to-top !
   ! is simply the reverse, ord(k) = n-k+1. Detect that in O(n) and take the reverse; fall back to   !
   ! the sort otherwise, because unit tests construct cohorts in arbitrary order.                    !
   !                                                                                                 !
   ! TIE-BREAK, and why the two branches agree exactly: the sort's `<=` keeps the LAST index         !
   ! achieving the running minimum, so among equal heights it emits the largest index first. In a   !
   ! descending array equal heights are consecutive and the largest remaining index is always       !
   ! minimal, so the reverse produces the identical permutation -- ties included.                    !
   !---------------------------------------------------------------------------------------!
   subroutine aero_bottom_to_top(acfg, aenv, ageom, n, height, lai, crown, leaf_width, branch_diam,    &
                                 leaf_temp, aero)
      type(aero_cfg_t),      intent(in)    :: acfg
      type(aero_env_t),      intent(in)    :: aenv
      type(aero_geom_t),     intent(in)    :: ageom
      integer(ik),           intent(in)    :: n
      real(wp),              intent(in)    :: height(:)       !< [m]  cohort height (gather order, top first)
      real(wp),              intent(in)    :: lai(:)          !< [m2/m2] leaf area index
      real(wp),              intent(in)    :: crown(:)        !< [-]  crown fraction
      real(wp),              intent(in)    :: leaf_width(:)   !< [m]  leaf width (boundary layer)
      real(wp),              intent(in)    :: branch_diam(:)  !< [m]  branch diameter (wood boundary layer)
      real(wp),              intent(in)    :: leaf_temp(:)    !< [K]  leaf temperature
      type(aero_out_t),      intent(inout) :: aero
      integer(ik) :: ord(n), k, j, imin
      real(wp)    :: hmin
      logical     :: used(n), descending
      real(wp)    :: h_bt(n), lai_bt(n), cr_bt(n), lt_bt(n), lw_bt(n), bd_bt(n)
      real(wp)    :: wind_bt(n), lgbh_bt(n), lgbw_bt(n), wgbh_bt(n), wgbw_bt(n)

      descending = .true.
      do j = 1_ik, n - 1_ik
         if (height(j) < height(j+1_ik)) then ; descending = .false. ; exit ; end if
      end do

      used = .false.
      do k = 1_ik, n
         if (descending) then
            imin = n - k + 1_ik                                  ! O(n) fast path
         else
            imin = 0_ik ; hmin = huge(1.0_wp)                    ! O(n^2) fallback (unsorted input)
            do j = 1_ik, n
               if (.not. used(j) .and. height(j) <= hmin) then ; hmin = height(j) ; imin = j ; end if
            end do
         end if
         ord(k)    = imin ; used(imin) = .true.
         h_bt(k)   = height(imin)     ; lai_bt(k) = lai(imin)
         cr_bt(k)  = crown(imin)      ; lt_bt(k)  = leaf_temp(imin)
         lw_bt(k)  = leaf_width(imin) ; bd_bt(k)  = branch_diam(imin)
      end do

      call canopy_aerodynamics(acfg, aenv, ageom, n, h_bt, lai_bt, cr_bt, lt_bt, lt_bt, lw_bt, bd_bt, aero)

      !----- aero%*(k) is now bottom->top; copy out, then scatter back to gather order. ----------!
      do k = 1_ik, n
         wind_bt(k) = aero%wind(k)     ; lgbh_bt(k) = aero%leaf_gbh(k) ; lgbw_bt(k) = aero%leaf_gbw(k)
         wgbh_bt(k) = aero%wood_gbh(k) ; wgbw_bt(k) = aero%wood_gbw(k)
      end do
      do k = 1_ik, n
         aero%wind(ord(k))     = wind_bt(k)
         aero%leaf_gbh(ord(k)) = lgbh_bt(k) ; aero%leaf_gbw(ord(k)) = lgbw_bt(k)
         aero%wood_gbh(ord(k)) = wgbh_bt(k) ; aero%wood_gbw(ord(k)) = wgbw_bt(k)
      end do
   end subroutine aero_bottom_to_top

end module meds_fast_prepass
