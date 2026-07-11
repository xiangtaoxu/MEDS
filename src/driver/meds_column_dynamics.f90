!==========================================================================================!
! meds_column_dynamics -- the fast-timescale (dt_fast) integrator that couples the fast          !
! biophysics processes into one sub-daily operator-split sweep (design Part II). The canopy-air-  !
! space (CAS) twins are the shared coupling reservoir: each dt_fast the aerodynamics kernel sets   !
! the conductances; the leaf and GROUND surface fluxes feed the CAS; the soil WATER column         !
! (Richards) and soil THERMAL column (heat) are advanced; and the three CAS twins (enthalpy /       !
! specific humidity / CO2) are advanced IMPLICITLY in the atmospheric-exchange term, gated by u*.    !
!                                                                                          !
! Coupled so far: aerodynamics -> {diagnostic leaf balance, ground balance} -> soil WATER column     !
! (infiltration / DSL soil-evap / root uptake / drainage) -> CAS enthalpy/vapour/CO2 (implicit) ->    !
! soil THERMAL column (implicit BE-Thomas), which reads the just-updated soil moisture. Carries the    !
! §3.5 fix (atm<->CAS conductance = profile-factored rho*ustar*temp1/temp2) AND the §3.6 GROUND-        !
! evaporation single-authority: the hydrology kernel's DSL/alpha_soil soil_evap is THE ground latent    !
! flux -- it drives the CAS vapour twin AND the ground energy balance's LE (no double-count).           !
!                                                                                          !
! Leaf temperature is DIAGNOSED from a linearized steady-state balance. CARBON is now REAL: leaf gas     !
! exchange gives GPP + stomatal conductance (driving transpiration) + leaf Rd; stem + fine-root           !
! maintenance respiration (autotrophic) and heterotrophic Rh assemble NEE = (Rd+stem+root) + Rh - GPP     !
! for the CAS CO2 twin. Plant HYDRAULICS is now REAL too: each cohort carries prognostic node water       !
! potentials psi(NODE_LEAF/WOOD) in patch_biophys_t; plant_water_flux advances them from the realized     !
! (supply-limited) transpiration demand and the root-weighted soil psi, and the updated psi_leaf feeds     !
! NEXT step's leaf gas exchange -- the soil -> plant -> stomata drought feedback (lagged one dt_fast).     !
! Still prescribed (next layers): absorbed radiation, precip. Canopy interception (leaf film), the         !
! stepper hook, and cross-demography persistence remain to wire.                                            !
!                                                                                          !
! WHOLE-COLUMN CONSERVATION -- verified by budg%whole_energy / budg%whole_water (Δ of ALL stores vs the     !
! true boundary fluxes; these CATCH cross-seam leaks the per-kernel budgets miss). Water-borne enthalpy    !
! is transported consistently: the CAS latent uses enthalpy_vapor(tl) (matching the CAS inverter + ground); !
! the soil sheds the transpiration water's liquid enthalpy via root_heat_sink; infiltration/drainage water  !
! carry internal_energy_liquid across the soil boundaries. INTER-LAYER advective heat: the hydrology kernel  !
! now EXPOSES the time-mean per-face Darcy flux (hflux%w_flux), and soil_energy_flux can advect the liquid   !
! enthalpy on it -- an OPT-IN coupling (cfg%advect_soil_heat, default OFF). It reconciles the soil moisture   !
! <-> energy coupling that the standalone kernels never exercised (their unit tests forced moisture constant):!
! with the internal-energy soil store referenced at tsupercool_liq (~193 K), a per-step Δθ carries a large    !
! liquid enthalpy, so turning the advection ON shifts the coupled soil/CAS temperatures NOTABLY. It still      !
! conserves the whole-column energy to machine precision and keeps soil temps bounded, but the MAGNITUDE      !
! awaits validation against a reference (ED2) with real met forcing -- hence gated OFF by default. A          !
! SUPPLY-limited leaf is not re-solved for the extra sensible (small drought residual). The stepper hook +    !
! cross-demography persistence remain to wire.                                                                !
!==========================================================================================!
module meds_column_dynamics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mmdry, tiny_num, cp_air, stefan, latent_heat_vap, rho_h2o, r_gas
   use meds_config,           only : meds_config_t, SCHEME_SPLIT_SEQUENTIAL, SCHEME_PICARD_COUPLED, &
                                     INTEG_SPLIT, INTEG_ARK
   use meds_biophysics_types, only : aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t,          &
                                     alloc_aero_out, veg_thermal_params_t, patch_biophys_t,    &
                                     soil_params_t, soil_thermal_params_t, soil_opts_t,        &
                                     energy_forcing_t, energy_opts_t, energy_flux_t,           &
                                     soil_column_t, soil_energy_column_t, chydro_forcing_t, chydro_flux_t, &
                                     SOIL_BC_FREE_DRAIN
   use meds_column_derivs,    only : column_state_t, column_frozen_t, surface_state_t,         &
                                     surface_frozen_t, surface_tend_t, surface_derivs
   use meds_ark_stepper,      only : ark2_column_step, adaptive_ark_march
   use meds_canopy_aerodynamics, only : canopy_aerodynamics
   use meds_column_energy,    only : soil_energy_flux
   use meds_column_hydrology, only : column_hydrology_flux
   use meds_plant_interface,  only : leaf_env_t, leaf_flux_t, leaf_gas_exchange,               &
                                     wood_env_t, wood_params_t, wood_flux_t, stem_maintenance_respiration, &
                                     root_env_t, root_params_t, root_flux_t, fine_root_maintenance_respiration, &
                                     hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t,  &
                                     plant_water_flux, N_HYDRO, NODE_LEAF, NODE_WOOD
   use meds_column_co2,       only : heterotrophic_respiration_flux
   use meds_biogeochem_types, only : co2_opts_t
   use meds_thermo,           only : cas_temp_of_enthalpy, sat_specific_humidity,             &
                                     d_sat_vapor_pressure_dt, enthalpy_vapor, internal_energy_liquid,  &
                                     sat_vapor_pressure, uext_to_temp
   use meds_budget_check,     only : budget_t, budget_accumulate, closure_ok
   implicit none
   private

   public :: column_config_t, column_cohort_t, column_forcing_t, column_budget_t
   public :: alloc_column_cohort, column_fast_step, aero_bottom_to_top

   !----- Static per-run column configuration (built once; constant across dt_fast steps). ----!
   type :: column_config_t
      type(aero_cfg_t)            :: aero            !< aerodynamics constants
      type(veg_thermal_params_t)  :: veg_thermal    !< leaf/wood thermal params
      type(soil_params_t)         :: soil           !< soil geometry + texture (n_active layers)
      type(soil_thermal_params_t) :: soil_thermal   !< soil thermal texture
      type(energy_opts_t)         :: energy         !< soil-thermal solver options
      type(soil_opts_t)           :: hydro          !< soil-water (Richards) solver options
      type(wood_params_t)         :: wood           !< stem-respiration parameters
      type(root_params_t)         :: root           !< fine-root-respiration parameters
      type(co2_opts_t)            :: co2            !< heterotrophic-respiration options
      type(hydro_params_t)        :: hydro_p        !< plant-hydraulics parameters (PV curves, vulnerability)
      type(hydro_opts_t)          :: hydro_o        !< plant-hydraulics solver options
      real(wp)                    :: fast_soil_carbon = 5.0_wp   !< [kgC/m2] decomposable soil-C pool (prescribed, MVP)
      real(wp)                    :: rhizo_cond       = 5.0e-4_wp !< [kg/s/MPa] soil->root conductance (prescribed, MVP)
      logical                     :: advect_soil_heat = .false.  !< opt-in: advect liquid enthalpy on the interior
                                                                 !< per-layer Darcy flux (moisture<->energy coupling)
      !----- P3 coupled-surface (Picard) solver knobs; only consulted under SCHEME_PICARD_COUPLED. !
      integer(ik) :: picard_max_iter = 20_ik        !< outer-iteration cap
      real(wp)    :: picard_tol_temp = 1.0e-3_wp     !< [K]     temperature convergence tolerance
      real(wp)    :: picard_tol_shv  = 1.0e-6_wp     !< [kg/kg] CAS specific-humidity convergence tolerance
      real(wp)    :: picard_relax    = 0.5_wp        !< [-]     under-relaxation of the next-pass seed. The CAS<->ground
                                                     !<         sensible coupling gives the fixed-point map a slope ~ -1
                                                     !<         (oscillatory); 0.5 makes the relaxed map a contraction.
      logical     :: picard_fixed_iter = .false.     !< run a uniform pass count (GPU warp-uniform; no early exit)
      integer(ik) :: leaf_energy_model  = 0_ik       !< LEAFEN_DIAGNOSTIC (0) | LEAFEN_PROGNOSTIC (1)
      integer(ik) :: soil_water_coupling = 0_ik      !< SOILH2O_LAGGED (0) | SOILH2O_COUPLED (1)
   end type column_config_t

   !----- Leaf thermal model + soil-water-in-loop selector codes (P3). -----------------------!
   integer(ik), parameter, public :: LEAFEN_DIAGNOSTIC = 0_ik  !< steady-state leaf (tl = tcas + dtl)
   integer(ik), parameter, public :: LEAFEN_PROGNOSTIC = 1_ik  !< prognostic leaf_energy via veg_energy_balance (P3e)
   integer(ik), parameter, public :: SOILH2O_LAGGED    = 0_ik  !< soil water/hydraulics frozen per sub-step
   integer(ik), parameter, public :: SOILH2O_COUPLED   = 1_ik  !< soil water re-solved inside the Picard loop (P3f)

   !----- Per-patch cohort state (SoA; the demographic slice the fast loop consumes). ---------!
   type :: column_cohort_t
      integer(ik)              :: n = 0_ik
      integer(ik), allocatable :: pft(:)                       !< PFT index (into cfg%pft)
      real(wp),    allocatable :: lai(:), wai(:), height(:), crown(:)
      real(wp),    allocatable :: leaf_width(:), branch_diam(:)
      real(wp),    allocatable :: leaf_area(:), nplant(:), dbh(:), broot(:)   !< [m2/plant],[plant/m2],[cm],[kgC/plant]
      real(wp),    allocatable :: bleaf(:), bsap(:), sap_area(:)              !< [kgC/plant],[kgC/plant],[m2] (hydraulics)
   end type column_cohort_t

   !----- Prescribed per-step forcing the higher layers (RT, met) supply; photosynthesis/    !
   !      respiration/NEE are now computed from the plant kernels (no longer prescribed).     !
   type :: column_forcing_t
      real(wp)              :: enthalpy_atm  = 0.0_wp   !< [J/kg]     reference-level specific enthalpy
      real(wp)              :: shv_atm       = 0.0_wp   !< [kg/kg]    reference-level specific humidity
      real(wp)              :: co2_atm       = 400.0_wp !< [umol/mol] free-atmosphere CO2
      real(wp)              :: abs_sw_ground = 0.0_wp   !< [W/m2] shortwave reaching the ground
      real(wp)              :: abs_lw_ground = 0.0_wp   !< [W/m2] net longwave at the ground
      real(wp)              :: precip        = 0.0_wp   !< [kg/m2/s] ground-reaching rainfall (interception deferred)
      real(wp)              :: par_per_w     = 2.1_wp   !< [umol photon / (W absorbed)] absorbed->PAR-photon factor
      real(wp), allocatable :: abs_sw(:), abs_lw(:)     !< [W/m2] absorbed SW (VIS+NIR) / net LW per cohort (leaf ENERGY)
      real(wp), allocatable :: abs_par(:)               !< [W/m2] INCIDENT-equiv PAR (VIS) per cohort; the leaf
                                                        !< model re-applies leaf_absorptance internally (PHOTOSYNTHESIS)
   end type column_forcing_t

   !----- The per-patch conservation budgets (one place; the driver accumulates the closed resids).!
   !      The per-kernel budgets close BY CONSTRUCTION; whole_energy/whole_water are the CROSS-      !
   !      seam column totals (Δ all stores vs the true boundary fluxes) that actually catch leaks.   !
   type :: column_budget_t
      type(budget_t) :: cas_energy, cas_water, cas_co2, soil_energy, soil_water
      type(budget_t) :: whole_energy, whole_water
      real(wp)       :: gpp_last = 0.0_wp, nee_last = 0.0_wp   !< [umol/m2/s] last-step diagnostics
      !----- P3 Picard diagnostics (reporting only; not conserved state). --------------------!
      integer(ik)    :: picard_iters       = 0_ik    !< worst outer-iteration count over the sub-steps
      integer(ik)    :: picard_nonconv     = 0_ik    !< number of sub-steps that hit picard_max_iter unconverged
      real(wp)       :: picard_worst_resid = 0.0_wp  !< [K] worst residual temperature at exit
   end type column_budget_t

contains

   !----- Allocate a column_cohort_t (the per-patch cohort SoA the fast loop consumes). ------!
   subroutine alloc_column_cohort(coh, n)
      type(column_cohort_t), intent(out) :: coh
      integer(ik),           intent(in)  :: n
      coh%n = n
      allocate(coh%pft(n), coh%lai(n), coh%wai(n), coh%height(n), coh%crown(n),                &
               coh%leaf_width(n), coh%branch_diam(n), coh%leaf_area(n), coh%nplant(n),         &
               coh%dbh(n), coh%broot(n), coh%bleaf(n), coh%bsap(n), coh%sap_area(n))
      coh%pft = 1_ik
      coh%lai = 0.0_wp ; coh%wai = 0.0_wp ; coh%height = 0.0_wp ; coh%crown = 1.0_wp
      coh%leaf_width = 0.04_wp ; coh%branch_diam = 0.02_wp
      coh%leaf_area = 0.0_wp ; coh%nplant = 0.0_wp ; coh%dbh = 0.0_wp ; coh%broot = 0.0_wp
      coh%bleaf = 0.0_wp ; coh%bsap = 0.0_wp ; coh%sap_area = 0.0_wp
   end subroutine alloc_column_cohort

   !=======================================================================================!
   !  One fast (dt_fast) operator-split sweep for a single patch. Cohort arrays BOTTOM(1)->TOP. !
   !  Leaf gas exchange (real GPP + stomata + leaf Rd), stem/root maintenance respiration and    !
   !  heterotrophic Rh feed a physically-decomposed NEE = (Rd_leaf + stem + root) + Rh - GPP.    !
   !=======================================================================================!
   subroutine column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, gpp_coh, &
                               leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters)
      real(wp),                intent(in)    :: dt_fast
      type(meds_config_t),     intent(in)    :: cfg          !< PFT traits for leaf gas exchange
      type(column_config_t),   intent(in)    :: ccfg
      type(aero_env_t),        intent(inout) :: aenv         !< can_* fields refreshed from CAS state
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: coh
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(inout) :: bio
      type(aero_out_t),        intent(inout) :: aero         !< preallocated (alloc_aero_out)
      type(column_budget_t),   intent(inout) :: budg
      real(wp), optional,      intent(out)   :: gpp_coh(:)   !< [umol CO2/plant/s] per-cohort GROSS GPP (fast->slow)
      real(wp), optional,      intent(out)   :: leaf_resp_coh(:) !< [umol CO2/plant/s] leaf dark respiration (fast->slow)
      real(wp), optional,      intent(out)   :: stem_resp_coh(:) !< [umol CO2/plant/s] stem maintenance resp (fast->slow)
      real(wp), optional,      intent(out)   :: root_resp_coh(:) !< [umol CO2/plant/s] fine-root maint. resp (fast->slow)
      logical,     optional,   intent(out)   :: converged    !< Picard converged this sub-step (true for split)
      integer(ik), optional,   intent(out)   :: iters        !< outer-iteration count taken (1 for split)

      type(chydro_forcing_t) :: hforc
      type(chydro_flux_t)    :: hflux
      type(energy_forcing_t) :: eforc
      type(energy_flux_t)    :: sflux
      type(leaf_env_t)       :: lenv
      type(leaf_flux_t)      :: lf
      type(wood_env_t)       :: wenv
      type(wood_flux_t)      :: wf
      type(root_env_t)       :: renv
      type(root_flux_t)      :: rf
      type(hydro_env_t)      :: henv
      type(hydro_flux_t)     :: hfx
      real(wp)               :: transp_c(coh%n)     !< [kg/m2 ground/s] per-cohort transpiration demand (automatic)
      real(wp)               :: h_coeff_f(coh%n), g_tr_f(coh%n), leaf_in(coh%n)   !< frozen coeffs + prev-iterate leaf temp
      real(wp)               :: t_emit(coh%n)      !< LW emission base (start leaf_temp; matches the RT tcan_bt, P3c)
      real(wp)    :: te
      type(soil_column_t)        :: soil_w_n        !< snapshot of the soil-water column at state^n (Picard reset)
      type(soil_energy_column_t) :: soil_e_n        !< snapshot of the soil thermal column at state^n (Picard reset)
      real(wp), allocatable  :: psi_n(:,:)          !< snapshot of the per-cohort plant water potentials at state^n
      real(wp)    :: soil_psi_root, tground_in, t_ground_dia, t_bot_dia
      real(wp)    :: tcas, qcas, press, rho, le_slope, lw_slope, qsat_c, dqdt
      real(wp)    :: le_ref, dtl, tl, transp_i, gsw_ms, e_air, rho_mol
      real(wp)    :: resid_T, tcas_in, qcas_in, tcas_new
      real(wp), parameter :: LAI_SLAVE_MIN = 1.0e-2_wp    !< [m2/m2] below this a cohort is slaved to tcas (Picard)
      integer(ik) :: iter, niter, niter_taken
      logical     :: picard, nconv
      real(wp)    :: coh_h, coh_qw, coh_qsoil, coh_transp, coh_rnet
      real(wp)    :: gpp, ra_leaf, ra_stem, ra_root, rh, nee_biotic, soil_temp_root, theta_mean
      real(wp)    :: t_ground, t_bot, g_top, h_ground, le_ground, soil_evap, rain_temp
      real(wp)    :: gah, gaw, gac, wcap, ccap, can_dmol, src_enth, src_vap, src_frac
      real(wp)    :: enth0, shv0, co20, enth1, shv1, co21
      real(wp)    :: e_soil0, e_soil1, w_soil0, w_soil1, e_in, e_out, w_in, w_out
      integer(ik) :: i, n, nsl, k

      n = coh%n ; nsl = ccfg%soil%n_active

      !----- TIME-INTEGRATOR dispatch (inserted BEFORE the first bio mutation, so the split path       !
      !      below is byte-for-byte unentered -- the golden anchor is preserved structurally). The     !
      !      coupled IMEX-ARK path is opt-in ([fast].time_integrator="ark"); default is the split.     !
      if (cfg%time_integrator == INTEG_ARK) then
         call column_fast_step_ark(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg,   &
                                   gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters)
         return
      end if

      picard = (cfg%integration_scheme == SCHEME_PICARD_COUPLED)
      niter  = 1_ik ; if (picard) niter = max(1_ik, ccfg%picard_max_iter)
      !----- The prognostic-leaf option (leaf_energy SoA + veg_energy_balance) is DEFERRED (P3e);   !
      !      guard it rather than silently running the diagnostic leaf. The diagnostic leaf is the   !
      !      correct default at dt_fast ~ 900 s (leaf thermal inertia negligible).                    !
      if (ccfg%leaf_energy_model == LEAFEN_PROGNOSTIC) &
         error stop 'column_fast_step: leaf_energy_model="prognostic" not yet implemented (P3e); use "diagnostic"'
      !----- The soil-water coupling selector: LAGGED and COUPLED both currently re-solve the soil    !
      !      water from state^n each Picard pass (required for conservation while the leaf demand      !
      !      iterates). A true frozen/lagged optimization (thermal-only, cheaper) is deferred (P3f).   !

      !----- Snapshot start-of-step SOIL stores (for the whole-column budgets). --------------!
      e_soil0 = 0.0_wp ; w_soil0 = bio%soil_w%w_surface
      do k = 1_ik, nsl
         e_soil0 = e_soil0 + bio%soil_e%soil_energy(k) * ccfg%soil%dz(k)
         w_soil0 = w_soil0 + bio%soil_w%theta(k) * ccfg%soil%dz(k) * rho_h2o
      end do

      !----- 1. Refresh the aerodynamics env from the current CAS state, then solve. ---------!
      bio%cas%can_temp = cas_temp_of_enthalpy(bio%cas%can_enthalpy, bio%cas%can_shv)
      tcas = bio%cas%can_temp ; qcas = bio%cas%can_shv ; press = aenv%press ; rho = aenv%rho_air
      t_ground = bio%soil_e%soil_temp(1) ; t_bot = bio%soil_e%soil_temp(nsl) ; rain_temp = tcas
      aenv%can_temp = tcas ; aenv%can_theta = tcas ; aenv%can_shv = qcas ; aenv%can_co2 = bio%cas%can_co2
      aenv%t_ground = t_ground
      !      canopy_aerodynamics expects cohorts ordered BOTTOM(1)->TOP(n) (its wind cascade walks !
      !      top->bottom), but the column buffer is gathered height-DESCENDING (index 1 = top). Feed !
      !      it the bottom->top order and scatter the per-cohort conductances back to gather order.   !
      call aero_bottom_to_top(ccfg%aero, aenv, ageom, n, coh, bio%leaf_temp, aero)

      !----- Root-weighted soil temperature + column-mean moisture (root / heterotrophic resp). !
      soil_temp_root = 0.0_wp ; theta_mean = 0.0_wp
      do k = 1_ik, nsl
         soil_temp_root = soil_temp_root + bio%soil_e%soil_temp(k) * ccfg%soil%root_frac(k)
         theta_mean     = theta_mean     + bio%soil_w%theta(k) * ccfg%soil%dz(k)
      end do
      theta_mean = theta_mean / max(-ccfg%soil%soil_layer_z(nsl+1_ik), tiny_num)

      !----- 2. PRE-PASS (once per sub-step; ED2 freezes gs/hydraulics per DTLSM): LEAF gas      !
      !         exchange (GPP/gs/Rd), the FROZEN per-cohort leaf-energy coefficients h_coeff_f/    !
      !         g_tr_f, and stem/root maintenance respiration. These do NOT change across the      !
      !         Picard passes (they use the lagged, start-of-sub-step leaf_temp).                  !
      gpp = 0.0_wp ; ra_leaf = 0.0_wp ; ra_stem = 0.0_wp ; ra_root = 0.0_wp
      if (present(gpp_coh))       gpp_coh(1:n)       = 0.0_wp
      if (present(leaf_resp_coh)) leaf_resp_coh(1:n) = 0.0_wp
      if (present(stem_resp_coh)) stem_resp_coh(1:n) = 0.0_wp
      if (present(root_resp_coh)) root_resp_coh(1:n) = 0.0_wp
      do i = 1_ik, n
         rho_mol       = press / (r_gas * bio%leaf_temp(i))                     ! [mol/m3] molar air density
         e_air         = qcas * press / (0.622_wp + 0.378_wp * qcas)            ! [Pa] canopy-air vapour pressure
         lenv%par      = forc%abs_par(i) / max(coh%lai(i), 0.1_wp) * forc%par_per_w   ! absorbed PAR (VIS), not total SW
         lenv%leaf_temp = bio%leaf_temp(i)
         t_emit(i)      = bio%leaf_temp(i)     ! start-of-sub-step leaf temp = the RT LW emission base (P3c)
         lenv%vpd      = max(sat_vapor_pressure(bio%leaf_temp(i)) - e_air, 0.0_wp)
         lenv%ca       = bio%cas%can_co2
         lenv%pressure = press
         lenv%psi_leaf = bio%psi(NODE_LEAF, i)                                  ! lagged plant water status (hydraulics)
         lenv%gb       = aero%leaf_gbw(i) * rho_mol                             ! m/s -> mol H2O/m2/s
         call leaf_gas_exchange(lenv, cfg, coh%pft(i), lf)
         gsw_ms        = lf%gs / max(rho_mol, tiny_num)                         ! mol/m2/s -> m/s
         gpp           = gpp     + lf%a_gross * coh%leaf_area(i) * coh%nplant(i)
         if (present(gpp_coh)) gpp_coh(i) = lf%a_gross * coh%leaf_area(i)          ! [umol/plant/s] per-plant gross
         ra_leaf       = ra_leaf + lf%rd      * coh%leaf_area(i) * coh%nplant(i)
         if (present(leaf_resp_coh)) leaf_resp_coh(i) = lf%rd * coh%leaf_area(i)   ! [umol/plant/s] per-plant leaf Rd
         !----- frozen leaf-energy coefficients (sensible + transpiration series conductance). --!
         h_coeff_f(i) = ccfg%veg_thermal%effarea_heat * coh%lai(i) * aero%leaf_gbh(i) * rho * cp_air
         g_tr_f(i)    = 0.0_wp
         if (aero%leaf_gbw(i) + gsw_ms > tiny_num) then
            g_tr_f(i) = ccfg%veg_thermal%effarea_transp * coh%lai(i)                            &
                        * aero%leaf_gbw(i) * gsw_ms / (aero%leaf_gbw(i) + gsw_ms)
         end if
         !----- Autotrophic maintenance respiration: stem + fine root (per plant -> per m2). ---!
         wenv%wood_temp = bio%leaf_temp(i) ; wenv%dbh = coh%dbh(i) ; wenv%height = coh%height(i)
         wenv%wai = coh%wai(i) ; wenv%nplant = coh%nplant(i)
         call stem_maintenance_respiration(wenv, ccfg%wood, wf)
         renv%soil_temp = soil_temp_root ; renv%broot = coh%broot(i)
         call fine_root_maintenance_respiration(renv, ccfg%root, rf)
         ra_stem = ra_stem + wf%stem_resp * coh%nplant(i)
         ra_root = ra_root + rf%root_resp * coh%nplant(i)
         if (present(stem_resp_coh)) stem_resp_coh(i) = wf%stem_resp   ! [umol/plant/s] already per-plant
         if (present(root_resp_coh)) root_resp_coh(i) = rf%root_resp   ! [umol/plant/s] already per-plant
      end do

      !----- NEE = autotrophic (leaf Rd + stem + root) + heterotrophic Rh - GPP (thermally      !
      !      passive: frozen across passes; the CO2 twin is solved once after convergence). ----!
      rh = heterotrophic_respiration_flux(ccfg%fast_soil_carbon, soil_temp_root, theta_mean,   &
                                          ccfg%soil%theta_res(1), ccfg%soil%theta_sat(1), ccfg%co2)
      nee_biotic = ra_leaf + ra_stem + ra_root + rh - gpp
      budg%gpp_last = gpp ; budg%nee_last = nee_biotic

      !----- CAS capacities + atm-exchange conductances (frozen across passes, §3.5). ---------!
      can_dmol = rho * (1.0_wp - qcas) / mmdry
      wcap     = rho      * bio%cas%can_depth
      ccap     = can_dmol * bio%cas%can_depth
      gah      = rho      * aero%ustar * aero%temp1
      gaw      = rho      * aero%ustar * aero%temp2
      gac      = can_dmol * aero%ustar * aero%temp2

      !----- Snapshot state^n once: the Picard passes re-solve the SAME backward-Euler steps FROM  !
      !      this base each iteration (only the source, re-evaluated at the iterate, changes). The   !
      !      soil-water column + plant-psi are prognostic and column_hydrology_flux / plant_water_flux !
      !      ADVANCE them, so they must be reset to state^n before each re-solve or they double-step.  !
      enth0 = bio%cas%can_enthalpy ; shv0 = bio%cas%can_shv ; co20 = bio%cas%can_co2
      soil_w_n = bio%soil_w ; soil_e_n = bio%soil_e ; psi_n = bio%psi

      !======================================================================================!
      !  3. Outer PICARD fixed point over { leaf energy -> soil water/src_frac -> CAS twins }.   !
      !     niter = 1 under SCHEME_SPLIT_SEQUENTIAL reproduces the operator-split sweep EXACTLY   !
      !     (one pass, no convergence test, soil water solved that pass). Under PICARD the block   !
      !     iterates at the current tcas/qcas until the store temperatures converge; the soil       !
      !     water/src_frac/hydraulics are frozen after pass 1 (SOILH2O_LAGGED) unless coupled.      !
      !======================================================================================!
      src_frac = 1.0_wp ; soil_evap = 0.0_wp ; nconv = .false. ; resid_T = 0.0_wp
      niter_taken = 0_ik
      do iter = 1_ik, niter
         niter_taken = iter
         tcas_in = tcas ; qcas_in = qcas ; leaf_in(1:n) = bio%leaf_temp(1:n)

         !----- 3a. Leaf energy balance (diagnostic) at the CURRENT tcas/qcas, frozen coeffs. --!
         qsat_c = sat_specific_humidity(tcas, press)
         dqdt   = 0.622_wp * press / max((press - 0.378_wp * sat_vapor_pressure(tcas))**2, tiny_num) &
                  * d_sat_vapor_pressure_dt(tcas)
         coh_h = 0.0_wp ; coh_qw = 0.0_wp ; coh_qsoil = 0.0_wp ; coh_transp = 0.0_wp ; coh_rnet = 0.0_wp
         do i = 1_ik, n
            if (picard .and. coh%lai(i) < LAI_SLAVE_MIN) then    ! near-zero LAI: slave to CAS, no exchange
               bio%leaf_temp(i) = tcas ; transp_c(i) = 0.0_wp ; cycle
            end if
            !----- LW emission linearized around T_emit: tcas for SPLIT (reduces to the current    !
            !      form, so split stays bit-identical) or the start leaf_temp for PICARD (which the   !
            !      two-stream also emits at via tcan_bt, P3c) -> leaf emission consistent at leaf_temp. !
            te = tcas ; if (picard) te = t_emit(i)
            lw_slope = 4.0_wp * ccfg%veg_thermal%leaf_emiss * stefan * te**3 * coh%lai(i)
            le_slope = latent_heat_vap * rho * g_tr_f(i) * dqdt
            le_ref   = latent_heat_vap * rho * g_tr_f(i) * (qsat_c - qcas)
            dtl = (forc%abs_sw(i) + forc%abs_lw(i) - le_ref - lw_slope * (tcas - te))                &
                  / max(h_coeff_f(i) + le_slope + lw_slope, tiny_num)
            tl  = tcas + dtl
            bio%leaf_temp(i) = tl
            transp_i    = (le_ref + le_slope * dtl) / latent_heat_vap
            transp_c(i) = transp_i                                                       ! per-cohort demand (hydraulics)
            coh_h      = coh_h      + h_coeff_f(i) * dtl
            coh_qw     = coh_qw     + transp_i * enthalpy_vapor(tl)                       ! CAS latent (vapour enthalpy)
            coh_qsoil  = coh_qsoil  + transp_i * (enthalpy_vapor(tl) - latent_heat_vap)   ! liquid enthalpy soil sheds
            coh_transp = coh_transp + transp_i
            !----- NET leaf radiation. Write (tl - te) as ((tcas - te) + dtl): for SPLIT (te = tcas) !
            !      the first term is EXACTLY 0.0, so this is bit-identical to the old lw_slope*dtl     !
            !      (whereas (tcas+dtl)-tcas would carry a rounding ulp); identical value for PICARD.   !
            coh_rnet   = coh_rnet   + (forc%abs_sw(i) + forc%abs_lw(i) - lw_slope * ((tcas - te) + dtl))
         end do

         !----- 3b. Soil WATER column + supply limiter + plant hydraulics, RE-SOLVED FROM state^n  !
         !          each pass so the realized transpiration (= root uptake) stays consistent with    !
         !          the iterated leaf demand -- freezing it while the demand iterates would leak     !
         !          water/enthalpy. pass 1 is already at state^n; later passes reset the prognostic   !
         !          water + psi to the snapshot before re-advancing.                                  !
         if (iter > 1_ik) then
            bio%soil_w = soil_w_n ; bio%psi = psi_n
         end if
         hforc%precip_ground      = forc%precip
         hforc%root_uptake(1:nsl) = coh_transp * ccfg%soil%root_frac(1:nsl)
         hforc%t_ground           = t_ground
         hforc%q_air              = qcas
         hforc%rho_air            = rho
         hforc%r_aero             = 1.0_wp / max(aero%ggnet, tiny_num)   ! §3.6: r_aero = 1/ggnet
         call column_hydrology_flux(bio%soil_w, hforc, ccfg%soil, ccfg%hydro, dt_fast, hflux)
         soil_evap = hflux%soil_evap                                     ! §3.6: THE ground latent authority
         src_frac  = 1.0_wp
         if (coh_transp > tiny_num) src_frac = min(1.0_wp, hflux%uptake_total / coh_transp)
         soil_psi_root = 0.0_wp
         do k = 1_ik, nsl
            soil_psi_root = soil_psi_root + hflux%psi_soil(k) * ccfg%soil%root_frac(k)
         end do
         do i = 1_ik, n
            henv%transp     = transp_c(i) * src_frac / max(coh%nplant(i), tiny_num)   ! [kg/plant/s]
            henv%soil_psi   = soil_psi_root
            henv%rhizo_cond = ccfg%rhizo_cond
            henv%bleaf      = coh%bleaf(i) ; henv%bsap = coh%bsap(i) ; henv%broot = coh%broot(i)
            henv%sap_area   = coh%sap_area(i) ; henv%height = coh%height(i) ; henv%leaf_area = coh%leaf_area(i)
            call plant_water_flux(henv, ccfg%hydro_p, ccfg%hydro_o, dt_fast, bio%psi(:,i), hfx)
         end do
         coh_qw     = coh_qw     * src_frac      ! the CAS gains only the water the soil gave up
         coh_qsoil  = coh_qsoil  * src_frac
         coh_transp = coh_transp * src_frac

         !----- 3c. GROUND surface: sensible at the current tcas + hydrology latent (frozen). ---!
         h_ground  = aero%ggnet * rho * cp_air * (t_ground - tcas)
         le_ground = soil_evap * enthalpy_vapor(t_ground)
         g_top     = forc%abs_sw_ground + forc%abs_lw_ground - h_ground - le_ground

         !----- 3d. CAS three-twin update: IMPLICIT atm exchange, FROM the state^n snapshot. ----!
         src_enth = coh_h + coh_qw + h_ground + le_ground                 ! [W/m2]  sensible + latent
         src_vap  = coh_transp + soil_evap                               ! [kg/m2/s] leaf transp + soil evap
         enth1 = (wcap*enth0 + dt_fast*(src_enth + gah*forc%enthalpy_atm)) / (wcap + dt_fast*gah)
         shv1  = (wcap*shv0  + dt_fast*(src_vap  + gaw*forc%shv_atm))       / (wcap + dt_fast*gaw)
         tcas_new = cas_temp_of_enthalpy(enth1, shv1)

         !----- 3d'. SOIL THERMAL step (P3b): reset to state^n, apply the water-enthalpy boundary  !
         !          advection (at THIS pass's t_ground/t_bot, so the budget uses the same values),   !
         !          then the BE heat diffusion with g_top; DIAGNOSE the surface/bottom temps for the  !
         !          NEXT pass. On the converged exit t_ground/t_bot stay at the values the advection   !
         !          + budget used (so split@niter=1 is bit-identical to the old post-loop step).       !
         if (iter > 1_ik) bio%soil_e = soil_e_n
         bio%soil_e%soil_energy(1)   = bio%soil_e%soil_energy(1)                                    &
              + (hflux%infiltration * internal_energy_liquid(rain_temp)                             &
                 - hflux%runoff_surf * internal_energy_liquid(t_ground)) * dt_fast / ccfg%soil%dz(1)
         bio%soil_e%soil_energy(nsl) = bio%soil_e%soil_energy(nsl)                                  &
              - hflux%drainage * internal_energy_liquid(t_bot) * dt_fast / ccfg%soil%dz(nsl)
         eforc%g_top = g_top ; eforc%geothermal = 0.0_wp
         eforc%soil_water(1:nsl) = bio%soil_w%theta(1:nsl)
         if (ccfg%advect_soil_heat) then
            eforc%w_flux(1:nsl) = -hflux%w_flux(1:nsl)      ! down-positive (hydro) -> up-positive (energy)
         else
            eforc%w_flux(1:nsl) = 0.0_wp                    ! interior advection lumped (validated baseline)
         end if
         eforc%root_heat_sink(1:nsl) = coh_qsoil * ccfg%soil%root_frac(1:nsl)   ! shed transpiration-water enthalpy
         call soil_energy_flux(bio%soil_e, eforc, ccfg%soil_thermal, ccfg%soil, ccfg%energy, dt_fast, sflux)
         t_ground_dia = bio%soil_e%soil_temp(1) ; t_bot_dia = bio%soil_e%soil_temp(nsl)

         !----- 3e. Convergence: inter-iterate temperature (CAS + leaf + ground) + CAS humidity. -!
         resid_T = abs(tcas_new - tcas_in)
         resid_T = max(resid_T, abs(t_ground_dia - t_ground))
         do i = 1_ik, n
            if (picard .and. coh%lai(i) < LAI_SLAVE_MIN) cycle    ! slaved cohorts excluded from the norm
            resid_T = max(resid_T, abs(bio%leaf_temp(i) - leaf_in(i)))
         end do
         if (.not. picard) then
            nconv = .true. ; exit                                 ! split: single pass is the answer
         else if (.not. ccfg%picard_fixed_iter .and. resid_T < ccfg%picard_tol_temp                &
                  .and. abs(shv1 - qcas_in) < ccfg%picard_tol_shv) then
            nconv = .true. ; exit                                 ! converged: keep this pass's t_ground/t_bot
         end if
         !----- Seed the NEXT pass only (under-relax the CAS seed; committed enth1/shv1 stay exact;  !
         !      t_ground/t_bot take the freshly diagnosed values). Skipped on the LAST pass of a      !
         !      non-converged / picard_fixed_iter run so the post-loop budget + the last advection    !
         !      reference the SAME t_ground/t_bot (else a small energy asymmetry on those exits).      !
         if (iter < niter) then
            tcas = ccfg%picard_relax * tcas_new + (1.0_wp - ccfg%picard_relax) * tcas_in
            qcas = ccfg%picard_relax * shv1     + (1.0_wp - ccfg%picard_relax) * qcas_in
            t_ground = t_ground_dia ; t_bot = t_bot_dia
         end if
      end do
      if (picard .and. ccfg%picard_fixed_iter) nconv = .true.    ! fixed-count run: accept the last iterate

      !----- Commit the CAS (enth1/shv1 = exact BE solution at convergence) + the passive CO2 twin. !
      co21 = (ccap*co20 + dt_fast*(nee_biotic + gac*forc%co2_atm)) / (ccap + dt_fast*gac)
      bio%cas%can_enthalpy = enth1 ; bio%cas%can_shv = shv1 ; bio%cas%can_co2 = co21
      bio%cas%can_temp     = cas_temp_of_enthalpy(enth1, shv1)

      !----- Picard diagnostics + non-convergence contract (clamp = last iterate, never partial). !
      budg%picard_iters       = max(budg%picard_iters, niter_taken)
      budg%picard_worst_resid = max(budg%picard_worst_resid, resid_T)
      if (picard .and. .not. nconv) then
         budg%picard_nonconv = budg%picard_nonconv + 1_ik
         if (ccfg%energy%debug_error) error stop 'column_fast_step: Picard did not converge (debug_error)'
      end if
      if (present(converged)) converged = nconv
      if (present(iters))     iters     = niter_taken

      !----- The soil thermal step now runs INSIDE the Picard loop (§3d'); bio%soil_e, sflux, the   !
      !      converged g_top, and t_ground/t_bot (the values the last pass's advection used) are all  !
      !      final here, so the budgets below close against the consistent boundary fluxes.           !

      !----- 7. Per-kernel closed budgets (each closes by construction). ----------------------!
      call budget_accumulate(budg%cas_energy, wcap*enth0, wcap*enth1, src_enth + gah*forc%enthalpy_atm, &
                             gah*enth1, dt_fast, abs(wcap*enth1), 1.0e-8_wp, 1.0e-3_wp)
      call budget_accumulate(budg%cas_water,  wcap*shv0, wcap*shv1, src_vap + gaw*forc%shv_atm,        &
                             gaw*shv1, dt_fast, max(abs(wcap*shv1), 1.0e-6_wp), 1.0e-8_wp, 1.0e-10_wp)
      call budget_accumulate(budg%cas_co2,    ccap*co20, ccap*co21, nee_biotic + gac*forc%co2_atm,&
                             gac*co21, dt_fast, abs(ccap*co21), 1.0e-6_wp, 1.0e-3_wp)
      call track_resid(budg%soil_energy, sflux%energy_resid, abs(g_top)*dt_fast + 1.0_wp, 1.0e-6_wp, 1.0e-3_wp)
      call track_resid(budg%soil_water,  hflux%mass_resid,   1.0_wp,             1.0e-6_wp, 1.0e-4_wp)

      !----- 7b. WHOLE-COLUMN budgets: Δ(all stores) vs the TRUE boundary fluxes (catches leaks). !
      e_soil1 = 0.0_wp ; w_soil1 = bio%soil_w%w_surface
      do k = 1_ik, nsl
         e_soil1 = e_soil1 + bio%soil_e%soil_energy(k) * ccfg%soil%dz(k)
         w_soil1 = w_soil1 + bio%soil_w%theta(k) * ccfg%soil%dz(k) * rho_h2o
      end do
      !----- Water: precip IN; drainage + runoff + atm-vapour OUT. --------------------------!
      w_in  = forc%precip
      w_out = hflux%drainage + hflux%runoff_surf + gaw * (shv1 - forc%shv_atm)
      call budget_accumulate(budg%whole_water, w_soil0 + wcap*shv0, w_soil1 + wcap*shv1, w_in, w_out, &
                             dt_fast, max(w_soil1 + wcap*shv1, 1.0_wp), 1.0e-6_wp, 1.0e-4_wp)
      !----- Energy: net radiation + rain enthalpy IN; atm exchange + drainage/runoff enthalpy OUT.!
      e_in  = coh_rnet + forc%abs_sw_ground + forc%abs_lw_ground                                     &
              + hflux%infiltration * internal_energy_liquid(rain_temp)   ! energy that reached the SOIL store
      e_out = gah * (enth1 - forc%enthalpy_atm) + hflux%drainage * internal_energy_liquid(t_bot)      &
              + hflux%runoff_surf * internal_energy_liquid(t_ground)
      call budget_accumulate(budg%whole_energy, e_soil0 + wcap*enth0, e_soil1 + wcap*enth1, e_in, e_out, &
                             dt_fast, abs(e_soil1 + wcap*enth1), 1.0e-6_wp, 1.0e0_wp)
   end subroutine column_fast_step

   !=======================================================================================!
   !  INTEG_ARK path: the coupled IMEX-ARK fast step (archive/MEDS_IMEX_ARK_DESIGN.md). Shares the   !
   !  split's frozen pre-pass (build_column_frozen), packs the state into the pure column vector,     !
   !  advances one dt_fast with the ARK stepper, then unpacks. MVP restriction: inert hydrology       !
   !  (free-drain, precip==0, no Zeng-Decker) -- column_state_t carries no ponding/aquifer/water-      !
   !  table state, and the ARK omits the split's soil-boundary water-enthalpy advection.              !
   !=======================================================================================!
   subroutine column_fast_step_ark(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg,  &
                                   gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh, converged, iters)
      real(wp),                intent(in)    :: dt_fast
      type(meds_config_t),     intent(in)    :: cfg
      type(column_config_t),   intent(in)    :: ccfg
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: coh
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(inout) :: bio
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budg
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      logical,     optional,   intent(out)   :: converged
      integer(ik), optional,   intent(out)   :: iters

      type(column_frozen_t)  :: fro
      type(column_state_t)   :: y, y_out, ycur, ytmp, yerr
      type(surface_state_t)  :: ys
      type(surface_frozen_t) :: fs
      type(surface_tend_t)   :: sf
      real(wp)    :: tg, fl, dt0
      integer(ik) :: n, nsl, k, isub, nsub, nsteps, nrej

      n = coh%n ; nsl = ccfg%soil%n_active

      !----- inert-hydrology guard (see header). -------------------------------------------------!
      if (ccfg%hydro%zeng_decker .or. ccfg%hydro%bottom_bc /= SOIL_BC_FREE_DRAIN .or.             &
          forc%precip > tiny_num)                                                                 &
         error stop 'column_fast_step_ark: INTEG_ARK requires inert hydrology (free-drain, precip==0, no zeng_decker)'

      call build_column_frozen(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, n, nsl, &
                               fro, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh)

      !----- advance one dt_fast: adaptive (embedded-error) or GPU-warp-uniform fixed substeps. ----!
      if (cfg%ark_adaptive) then
         dt0 = dt_fast ; if (cfg%ark_dt_init > tiny_num) dt0 = min(cfg%ark_dt_init, dt_fast)
         call adaptive_ark_march(y, fro, n, nsl, dt_fast, cfg%ark_rtol, dt0, y_out, nsteps, nrej)
      else
         nsub = max(1_ik, cfg%ark_fixed_substep) ; nrej = 0_ik ; ycur = y
         do isub = 1_ik, nsub
            call ark2_column_step(ycur, fro, n, nsl, dt_fast/real(nsub, wp), ytmp, yerr,          &
                                  niter=cfg%ark_niter, relax=cfg%ark_relax)
            ycur = ytmp
         end do
         y_out = ycur ; nsteps = nsub
      end if

      !----- unpack into bio + re-derive the diagnostic soil temperatures + leaf temperatures. -----!
      bio%cas%can_enthalpy = y_out%cas_enthalpy ; bio%cas%can_shv = y_out%cas_shv ; bio%cas%can_co2 = y_out%cas_co2
      bio%cas%can_temp = cas_temp_of_enthalpy(y_out%cas_enthalpy, y_out%cas_shv)
      bio%soil_e%soil_energy(1:nsl) = y_out%soil_energy(1:nsl)
      bio%soil_w%theta(1:nsl)       = y_out%theta(1:nsl)
      bio%psi(:, 1:n)               = y_out%psi(:, 1:n)
      do k = 1_ik, nsl
         call uext_to_temp(y_out%soil_energy(k), y_out%theta(k)*rho_h2o,                          &
                           ccfg%soil_thermal%soil_dry_heat_capacity(k), bio%soil_e%soil_temp(k), bio%soil_e%soil_fliq(k))
      end do
      call uext_to_temp(y_out%soil_energy(1), y_out%theta(1)*rho_h2o,                             &
                        ccfg%soil_thermal%soil_dry_heat_capacity(1), tg, fl)
      fs = fro%surf ; fs%t_ground = tg
      ys%cas_enthalpy = y_out%cas_enthalpy ; ys%cas_shv = y_out%cas_shv ; ys%cas_co2 = y_out%cas_co2
      call surface_derivs(ys, fs, n, sf)
      bio%leaf_temp(1:n) = sf%leaf_temp(1:n)

      !----- gpp_last/nee_last were set by build_column_frozen. The whole-column conservation ledger  !
      !      for the multi-substep ARK (b-weighted boundary fluxes) is the deferred follow-on, so the  !
      !      per-store budget accumulators are left untouched here (n_check = 0).                      !
      if (present(converged)) converged = (nrej == 0_ik)
      if (present(iters))     iters     = nsteps
   end subroutine column_fast_step_ark

   !----- Build the frozen pre-pass (leaf gas exchange / respiration / CAS caps / aero) + the frozen  !
   !      hydrology BCs into a column_frozen_t, and pack the prognostic state into a column_state_t.   !
   !      The pre-pass loop is a VERBATIM copy of column_fast_step's :237-314 (writing struct fields   !
   !      instead of locals), so gpp_coh/resp are bit-identical to the split; the split's inline       !
   !      pre-pass is untouched (the golden anchor stays byte-for-byte).                               !
   subroutine build_column_frozen(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg, n, nsl, &
                                  fro, y, gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh)
      real(wp),                intent(in)    :: dt_fast
      type(meds_config_t),     intent(in)    :: cfg
      type(column_config_t),   intent(in)    :: ccfg
      type(aero_env_t),        intent(inout) :: aenv
      type(aero_geom_t),       intent(in)    :: ageom
      type(column_cohort_t),   intent(in)    :: coh
      type(column_forcing_t),  intent(in)    :: forc
      type(patch_biophys_t),   intent(in)    :: bio
      type(aero_out_t),        intent(inout) :: aero
      type(column_budget_t),   intent(inout) :: budg
      integer(ik),             intent(in)    :: n, nsl
      type(column_frozen_t),   intent(out)   :: fro
      type(column_state_t),    intent(out)   :: y
      real(wp), optional,      intent(out)   :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)

      type(leaf_env_t)       :: lenv ; type(leaf_flux_t) :: lf
      type(wood_env_t)       :: wenv ; type(wood_flux_t) :: wf
      type(root_env_t)       :: renv ; type(root_flux_t) :: rf
      type(chydro_forcing_t) :: hforc ; type(chydro_flux_t) :: hflux
      type(soil_column_t)    :: soil_w_scratch
      type(surface_state_t)  :: ys ; type(surface_tend_t) :: sf0
      real(wp) :: tcas, qcas, press, rho, t_ground, rho_mol, e_air, gsw_ms, can_dmol
      real(wp) :: gpp, ra_leaf, ra_stem, ra_root, rh, nee_biotic, soil_temp_root, theta_mean
      integer(ik) :: i, k

      allocate(fro%surf%h_coeff_f(n), fro%surf%g_tr_f(n), fro%surf%abs_sw(n), fro%surf%abs_lw(n), fro%surf%lai(n))
      allocate(fro%psi_e(nsl), fro%nplant(n), fro%bleaf(n), fro%bsap(n), fro%broot(n),            &
               fro%sap_area(n), fro%height(n), fro%leaf_area(n))
      allocate(y%psi(N_HYDRO, n))

      !----- aerodynamics from the current CAS state. -------------------------------------------!
      tcas = cas_temp_of_enthalpy(bio%cas%can_enthalpy, bio%cas%can_shv)
      qcas = bio%cas%can_shv ; press = aenv%press ; rho = aenv%rho_air ; t_ground = bio%soil_e%soil_temp(1)
      aenv%can_temp = tcas ; aenv%can_theta = tcas ; aenv%can_shv = qcas ; aenv%can_co2 = bio%cas%can_co2
      aenv%t_ground = t_ground
      call aero_bottom_to_top(ccfg%aero, aenv, ageom, n, coh, bio%leaf_temp, aero)

      soil_temp_root = 0.0_wp ; theta_mean = 0.0_wp
      do k = 1_ik, nsl
         soil_temp_root = soil_temp_root + bio%soil_e%soil_temp(k) * ccfg%soil%root_frac(k)
         theta_mean     = theta_mean     + bio%soil_w%theta(k) * ccfg%soil%dz(k)
      end do
      theta_mean = theta_mean / max(-ccfg%soil%soil_layer_z(nsl+1_ik), tiny_num)

      !----- leaf gas exchange + frozen leaf-energy coefficients + maintenance respiration. -----!
      gpp = 0.0_wp ; ra_leaf = 0.0_wp ; ra_stem = 0.0_wp ; ra_root = 0.0_wp
      if (present(gpp_coh))       gpp_coh(1:n)       = 0.0_wp
      if (present(leaf_resp_coh)) leaf_resp_coh(1:n) = 0.0_wp
      if (present(stem_resp_coh)) stem_resp_coh(1:n) = 0.0_wp
      if (present(root_resp_coh)) root_resp_coh(1:n) = 0.0_wp
      do i = 1_ik, n
         rho_mol        = press / (r_gas * bio%leaf_temp(i))
         e_air          = qcas * press / (0.622_wp + 0.378_wp * qcas)
         lenv%par       = forc%abs_par(i) / max(coh%lai(i), 0.1_wp) * forc%par_per_w
         lenv%leaf_temp = bio%leaf_temp(i)
         lenv%vpd       = max(sat_vapor_pressure(bio%leaf_temp(i)) - e_air, 0.0_wp)
         lenv%ca        = bio%cas%can_co2 ; lenv%pressure = press
         lenv%psi_leaf  = bio%psi(NODE_LEAF, i)
         lenv%gb        = aero%leaf_gbw(i) * rho_mol
         call leaf_gas_exchange(lenv, cfg, coh%pft(i), lf)
         gsw_ms         = lf%gs / max(rho_mol, tiny_num)
         gpp            = gpp     + lf%a_gross * coh%leaf_area(i) * coh%nplant(i)
         if (present(gpp_coh)) gpp_coh(i) = lf%a_gross * coh%leaf_area(i)
         ra_leaf        = ra_leaf + lf%rd      * coh%leaf_area(i) * coh%nplant(i)
         if (present(leaf_resp_coh)) leaf_resp_coh(i) = lf%rd * coh%leaf_area(i)
         fro%surf%h_coeff_f(i) = ccfg%veg_thermal%effarea_heat * coh%lai(i) * aero%leaf_gbh(i) * rho * cp_air
         fro%surf%g_tr_f(i)    = 0.0_wp
         if (aero%leaf_gbw(i) + gsw_ms > tiny_num) then
            fro%surf%g_tr_f(i) = ccfg%veg_thermal%effarea_transp * coh%lai(i)                     &
                                 * aero%leaf_gbw(i) * gsw_ms / (aero%leaf_gbw(i) + gsw_ms)
         end if
         wenv%wood_temp = bio%leaf_temp(i) ; wenv%dbh = coh%dbh(i) ; wenv%height = coh%height(i)
         wenv%wai = coh%wai(i) ; wenv%nplant = coh%nplant(i)
         call stem_maintenance_respiration(wenv, ccfg%wood, wf)
         renv%soil_temp = soil_temp_root ; renv%broot = coh%broot(i)
         call fine_root_maintenance_respiration(renv, ccfg%root, rf)
         ra_stem = ra_stem + wf%stem_resp * coh%nplant(i)
         ra_root = ra_root + rf%root_resp * coh%nplant(i)
         if (present(stem_resp_coh)) stem_resp_coh(i) = wf%stem_resp
         if (present(root_resp_coh)) root_resp_coh(i) = rf%root_resp
         !----- per-cohort geometry + radiation the ARK frozen inputs need. --------------------!
         fro%surf%lai(i)    = coh%lai(i)
         fro%surf%abs_sw(i) = forc%abs_sw(i) ; fro%surf%abs_lw(i) = forc%abs_lw(i)
         fro%nplant(i)   = coh%nplant(i)  ; fro%bleaf(i)    = coh%bleaf(i)  ; fro%bsap(i) = coh%bsap(i)
         fro%broot(i)    = coh%broot(i)   ; fro%sap_area(i) = coh%sap_area(i)
         fro%height(i)   = coh%height(i)  ; fro%leaf_area(i) = coh%leaf_area(i)
      end do
      rh = heterotrophic_respiration_flux(ccfg%fast_soil_carbon, soil_temp_root, theta_mean,      &
                                          ccfg%soil%theta_res(1), ccfg%soil%theta_sat(1), ccfg%co2)
      nee_biotic = ra_leaf + ra_stem + ra_root + rh - gpp
      budg%gpp_last = gpp ; budg%nee_last = nee_biotic

      !----- CAS capacities + atm conductances + the rest of the frozen surface inputs. ---------!
      can_dmol = rho * (1.0_wp - qcas) / mmdry
      fro%surf%leaf_emiss = ccfg%veg_thermal%leaf_emiss
      fro%surf%wcap = rho * bio%cas%can_depth ; fro%surf%ccap = can_dmol * bio%cas%can_depth
      fro%surf%gah  = rho * aero%ustar * aero%temp1
      fro%surf%gaw  = rho * aero%ustar * aero%temp2
      fro%surf%gac  = can_dmol * aero%ustar * aero%temp2
      fro%surf%enth_atm = forc%enthalpy_atm ; fro%surf%shv_atm = forc%shv_atm ; fro%surf%co2_atm = forc%co2_atm
      fro%surf%nee_biotic = nee_biotic
      fro%surf%abs_sw_ground = forc%abs_sw_ground ; fro%surf%abs_lw_ground = forc%abs_lw_ground
      fro%surf%ggnet = aero%ggnet ; fro%surf%rho = rho ; fro%surf%press = press
      fro%surf%src_frac = 1.0_wp ; fro%surf%t_ground = t_ground

      !----- params + hydraulics BCs. -----------------------------------------------------------!
      fro%soil = ccfg%soil ; fro%therm = ccfg%soil_thermal ; fro%energy_opts = ccfg%energy
      fro%hydro_opts = ccfg%hydro ; fro%hydro_p = ccfg%hydro_p ; fro%hydro_o = ccfg%hydro_o
      fro%geothermal = 0.0_wp ; fro%rhizo_cond = ccfg%rhizo_cond ; fro%psi_e(1:nsl) = 0.0_wp

      !----- FROZEN hydrology BCs: total transp demand (surface_derivs @ state^n, src_frac=1), then a  !
      !      SCRATCH column_hydrology_flux for soil_evap / infiltration / psi_soil / uptake_total. ----!
      ys%cas_enthalpy = bio%cas%can_enthalpy ; ys%cas_shv = bio%cas%can_shv ; ys%cas_co2 = bio%cas%can_co2
      call surface_derivs(ys, fro%surf, n, sf0)
      hforc%precip_ground      = forc%precip
      hforc%root_uptake(1:nsl) = sf0%coh_transp * ccfg%soil%root_frac(1:nsl)
      hforc%t_ground           = t_ground ; hforc%q_air = qcas ; hforc%rho_air = rho
      hforc%r_aero             = 1.0_wp / max(aero%ggnet, tiny_num)
      soil_w_scratch = bio%soil_w
      call column_hydrology_flux(soil_w_scratch, hforc, ccfg%soil, ccfg%hydro, dt_fast, hflux)
      fro%surf%soil_evap = hflux%soil_evap
      fro%q_top          = (hflux%infiltration - hflux%soil_evap) / rho_h2o
      fro%soil_psi_root  = 0.0_wp
      do k = 1_ik, nsl
         fro%soil_psi_root = fro%soil_psi_root + hflux%psi_soil(k) * ccfg%soil%root_frac(k)
      end do
      if (sf0%coh_transp > tiny_num) fro%surf%src_frac = min(1.0_wp, hflux%uptake_total / sf0%coh_transp)

      !----- pack the prognostic state. ---------------------------------------------------------!
      y%cas_enthalpy = bio%cas%can_enthalpy ; y%cas_shv = bio%cas%can_shv ; y%cas_co2 = bio%cas%can_co2
      y%soil_energy(1:nsl) = bio%soil_e%soil_energy(1:nsl)
      y%theta(1:nsl)       = bio%soil_w%theta(1:nsl)
      y%psi(:, 1:n)        = bio%psi(:, 1:n)
   end subroutine build_column_frozen

   !----- Solve canopy aerodynamics with the cohort order it CONTRACTS for -- BOTTOM(1)->TOP(n)  !
   !      -- from the height-DESCENDING column buffer. Only the wind cascade + the per-cohort       !
   !      boundary layers depend on order; the whole-canopy scalars (ustar/temp1/temp2/uh) do not.  !
   !      An ascending-height permutation `ord` reverses the per-cohort inputs; the per-cohort wind  !
   !      and leaf/wood conductance outputs are scattered back to gather order. Identity for n<=1,   !
   !      so single-cohort behaviour is bit-unchanged.                                               !
   subroutine aero_bottom_to_top(acfg, aenv, ageom, n, coh, leaf_temp, aero)
      type(aero_cfg_t),      intent(in)    :: acfg
      type(aero_env_t),      intent(in)    :: aenv
      type(aero_geom_t),     intent(in)    :: ageom
      integer(ik),           intent(in)    :: n
      type(column_cohort_t), intent(in)    :: coh
      real(wp),              intent(in)    :: leaf_temp(:)
      type(aero_out_t),      intent(inout) :: aero
      integer(ik) :: ord(n), k, j, imin
      real(wp)    :: hmin
      logical     :: used(n)
      real(wp)    :: h_bt(n), lai_bt(n), cr_bt(n), lt_bt(n), lw_bt(n), bd_bt(n)
      real(wp)    :: wind_bt(n), lgbh_bt(n), lgbw_bt(n), wgbh_bt(n), wgbw_bt(n)

      !----- ord(k) = gather index of the k-th cohort counting from the canopy BOTTOM. ----------!
      used = .false.
      do k = 1_ik, n
         imin = 0_ik ; hmin = huge(1.0_wp)
         do j = 1_ik, n
            if (.not. used(j) .and. coh%height(j) <= hmin) then ; hmin = coh%height(j) ; imin = j ; end if
         end do
         ord(k)    = imin ; used(imin) = .true.
         h_bt(k)   = coh%height(imin)     ; lai_bt(k) = coh%lai(imin)
         cr_bt(k)  = coh%crown(imin)      ; lt_bt(k)  = leaf_temp(imin)
         lw_bt(k)  = coh%leaf_width(imin) ; bd_bt(k)  = coh%branch_diam(imin)
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

   !----- Track a kernel's own closed-budget residual into a budget_t (worst + fail count). ---!
   pure subroutine track_resid(b, resid, scale, rtol, atol)
      type(budget_t), intent(inout) :: b
      real(wp),       intent(in)    :: resid, scale, rtol, atol
      b%n_check = b%n_check + 1_ik
      b%resid   = resid
      b%worst   = max(b%worst, abs(resid))
      if (.not. closure_ok(resid, scale, rtol, atol)) b%n_fail = b%n_fail + 1_ik
   end subroutine track_resid

end module meds_column_dynamics
