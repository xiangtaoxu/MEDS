!==========================================================================================!
! meds_fast_loop -- the per-SITE fast-biophysics driver (the fast-cadence analogue of          !
! meds_vegetation_dynamics). It owns the ORCHESTRATION only: for each patch it gathers the       !
! demographic cohort slice into the column_cohort_t buffer, assembles a patch_biophys_t working  !
! bundle from the state-hub-owned per-patch reservoirs (cas/soil_e/soil_w), runs n_fast_per_slow  !
! operator-split sweeps of the per-patch kernel column_fast_step, and writes the evolved          !
! reservoirs back to the site. The static column_config_t + base met arrive via fast_context_t    !
! (the caller builds them -- no model parameters are hard-coded here); per-cohort leaf_temp/psi    !
! are RESEEDED each slow step (they relax within hours, so cross-slow-step memory is not carried   !
! -- the soil + CAS reservoirs hold the genuine multi-step memory).                                !
!                                                                                          !
! The stepper calls run_fast_biophysics before the slow loop when cfg%fast_biophysics_on. This is  !
! the fast->slow seam's fast half; the daily-GPP handoff into carbon growth lands in a later step. !
!==========================================================================================!
module meds_fast_loop
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, rho_h2o, umol_2_kgC, grav, cp_air
   use meds_config,           only : meds_config_t, SCHEME_PICARD_COUPLED
   use meds_thermo,           only : cas_enthalpy_of_temp, cas_temp_of_enthalpy, temp_to_uext
   use meds_demography_types, only : site_t
   use meds_biophysics_types, only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out,       &
                                     patch_biophys_t, alloc_patch_biophys, SOIL_RETENTION_VG
   use meds_soil_parameters,  only : build_soil_params
   use meds_soil_thermal,     only : build_soil_thermal
   use meds_column_dynamics,  only : column_config_t, column_cohort_t, column_forcing_t,        &
                                     column_budget_t, alloc_column_cohort, column_fast_step
   implicit none
   private

   public :: fast_context_t, init_fast_reservoirs, run_fast_biophysics, build_fast_context

   !----- Everything the fast driver needs beyond the site + cfg: the static column config plus !
   !      the reference met + initial soil state. The CALLER builds this (from TOML in the       !
   !      production path; the MVP holds constant, horizontally-uniform boundary conditions).    !
   type :: fast_context_t
      type(column_config_t) :: ccfg                 !< soil/thermal/hydro/aero/resp column config
      real(wp) :: u_ref   = 2.0_wp, zref = 30.0_wp  !< [m/s],[m] reference wind + height
      real(wp) :: press   = 101325.0_wp             !< [Pa]
      real(wp) :: rho_air = 1.2_wp                  !< [kg/m3]
      real(wp) :: air_temp = 288.0_wp               !< [K]        reference-level air temperature
      real(wp) :: shv_atm  = 0.008_wp               !< [kg/kg]    reference-level specific humidity
      real(wp) :: co2_atm  = 400.0_wp               !< [umol/mol] free-atmosphere CO2
      real(wp) :: rad_sw_top    = 400.0_wp          !< [W/m2] shortwave into the canopy (leaves)
      real(wp) :: rad_sw_ground = 60.0_wp           !< [W/m2] shortwave reaching the ground
      real(wp) :: precip        = 0.0_wp            !< [kg/m2/s] ground-reaching rainfall
      real(wp) :: theta_init      = 0.30_wp         !< [m3/m3] initial soil moisture (all layers)
      real(wp) :: soil_temp_init  = 288.0_wp        !< [K]     initial soil + CAS temperature
      real(wp) :: veg_height_bare = 1.0_wp          !< [m] canopy height for a cohort-free patch
   end type fast_context_t

contains

   !=======================================================================================!
   !  Build the fast-loop context for the production run. The scalar reference met + reservoir  !
   !  seeds keep the fast_context_t defaults (constant, horizontally-uniform MVP boundary          !
   !  conditions). The static column config `ccfg` is assembled here from DOCUMENTED MVP            !
   !  PLACEHOLDER parameters (soil texture/thermal, stem/root maintenance-respiration factors,      !
   !  heterotrophic-Rh and plant-hydraulics constants) -- there is no per-site column_config TOML   !
   !  loader yet, so these mirror the validated test_fast_loop configuration (they pass the fast-   !
   !  loop energy/water/CO2 budget checks). The fast loop is OPT-IN (cfg%fast_biophysics_on,         !
   !  default .false.), so this leaves the default demographic run unchanged. Follow-up: source      !
   !  these from a [column]/[soil] TOML block + per-PFT respiration traits.                          !
   !=======================================================================================!
   subroutine build_fast_context(cfg, ctx)
      type(meds_config_t),  intent(in)  :: cfg
      type(fast_context_t), intent(out) :: ctx
      integer(ik), parameter :: NSL_MVP = 10_ik      ! MVP soil-layer count (matches test_fast_loop)
      !----- Soil geometry/texture + thermal (van Genuchten; loam-ish MVP placeholders). -------!
      call build_soil_params(NSL_MVP, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,      &
                             2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, ctx%ccfg%soil)
      call build_soil_thermal(NSL_MVP, 3.0_wp, 0.15_wp, 2.0e6_wp, ctx%ccfg%soil_thermal)
      !----- Autotrophic maintenance-respiration + heterotrophic-Rh + prescribed soil-C pool. ---!
      ctx%ccfg%wood%is_woody = .true.
      ctx%ccfg%wood%stem_resp_factor25 = 0.06_wp ; ctx%ccfg%wood%agf_bs = 0.7_wp
      ctx%ccfg%root%root_resp_factor25 = 0.30_wp
      ctx%ccfg%co2%rh_k_base           = 0.01_wp
      ctx%ccfg%fast_soil_carbon        = 5.0_wp
      !----- Plant-hydraulics pressure-volume + vulnerability curves (MVP placeholders). --------!
      ctx%ccfg%hydro_p%leaf_pi0 = -1.5_wp ; ctx%ccfg%hydro_p%leaf_eps = 12.0_wp
      ctx%ccfg%hydro_p%leaf_af  = 0.30_wp ; ctx%ccfg%hydro_p%leaf_water_sat = 2.0_wp
      ctx%ccfg%hydro_p%wood_pi0 = -1.0_wp ; ctx%ccfg%hydro_p%wood_eps = 8.0_wp
      ctx%ccfg%hydro_p%wood_af  = 0.20_wp ; ctx%ccfg%hydro_p%wood_water_sat = 1.0_wp
      ctx%ccfg%hydro_p%wood_psi50 = -2.0_wp ; ctx%ccfg%hydro_p%wood_kexp = 2.0_wp
      ctx%ccfg%hydro_p%k_plant_max = 6.0e-4_wp ; ctx%ccfg%hydro_p%wood_kmax = 8.0_wp
      ctx%ccfg%hydro_p%vessel_curl = 1.5_wp
      ctx%ccfg%rhizo_cond = 5.0e-4_wp
   end subroutine build_fast_context

   !----- Seed every patch's fast reservoirs to a horizontally-uniform equilibrium. Called once !
   !      after the community is built (production) or before a fast-loop test. --------------!
   subroutine init_fast_reservoirs(site, ctx)
      type(site_t),         intent(inout) :: site
      type(fast_context_t), intent(in)    :: ctx
      integer(ik) :: ip, k, nsl
      nsl = ctx%ccfg%soil%n_active
      do ip = 1_ik, site%patch%n
         associate (cas => site%patch%cas(ip), se => site%patch%soil_e(ip), sw => site%patch%soil_w(ip))
            sw%theta(1:nsl)  = ctx%theta_init ; sw%w_surface = 0.0_wp ; sw%w_aquifer = 0.0_wp ; sw%z_wt = 0.0_wp
            do k = 1_ik, nsl
               se%soil_energy(k) = temp_to_uext(ctx%ccfg%soil_thermal%soil_dry_heat_capacity(k),    &
                                   ctx%theta_init * rho_h2o, ctx%soil_temp_init, 1.0_wp)
               se%soil_temp(k)   = ctx%soil_temp_init ; se%soil_fliq(k) = 1.0_wp
            end do
            cas%can_shv      = ctx%shv_atm ; cas%can_co2 = ctx%co2_atm
            cas%can_enthalpy = cas_enthalpy_of_temp(ctx%air_temp, ctx%shv_atm)
            cas%can_temp     = ctx%air_temp
         end associate
      end do
   end subroutine init_fast_reservoirs

   !=======================================================================================!
   !  Advance every patch of the site by one slow-step's worth of fast biophysics: n_fast_per_ !
   !  slow operator-split sweeps over the state-hub reservoirs, on constant (MVP) forcing.     !
   !  Optional out-args report the worst whole-column budget residuals + the fail count so a    !
   !  caller/test can assert conservation.                                                     !
   !=======================================================================================!
   subroutine run_fast_biophysics(site, ctx, cfg, worst_energy, worst_water, n_budget_fail)
      type(site_t),         intent(inout) :: site
      type(fast_context_t), intent(in)    :: ctx
      type(meds_config_t),  intent(in)    :: cfg
      real(wp),    optional, intent(out)  :: worst_energy, worst_water
      integer(ik), optional, intent(out)  :: n_budget_fail

      type(column_cohort_t)  :: coh
      type(column_forcing_t) :: forc
      type(aero_env_t)       :: aenv
      type(aero_geom_t)      :: ageom
      type(aero_out_t)       :: aero
      type(patch_biophys_t)  :: bio
      type(column_budget_t)  :: budg
      real(wp), allocatable  :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      real(wp)    :: we, ww, sum_lai
      integer(ik) :: ip, isub, j, i, i0, ncoh, nfail

      if (cfg%integration_scheme == SCHEME_PICARD_COUPLED) then
         error stop 'run_fast_biophysics: PICARD scheme not implemented (only split-sequential)'
      end if

      we = 0.0_wp ; ww = 0.0_wp ; nfail = 0_ik
      !----- Reset the fast->slow carbon accumulators (gross GPP + maintenance-resp losses) BEFORE !
      !      the fast window (carbon_growth reads them after; it has site intent(in), cannot reset).!
      site%cohort%gpp_accum(1:site%cohort%n)       = 0.0_wp
      site%cohort%leaf_resp_accum(1:site%cohort%n) = 0.0_wp
      site%cohort%stem_resp_accum(1:site%cohort%n) = 0.0_wp
      site%cohort%root_resp_accum(1:site%cohort%n) = 0.0_wp

      do ip = 1_ik, site%patch%n
         ncoh = site%patch%cohort_count(ip)
         i0   = site%patch%cohort_offset(ip)

         !----- Gather the patch's cohort slice into the column buffer (+ MVP derived inputs). !
         call alloc_column_cohort(coh, ncoh)
         if (allocated(gpp_coh)) deallocate(gpp_coh, leaf_resp_coh, stem_resp_coh, root_resp_coh)
         allocate(gpp_coh(ncoh), leaf_resp_coh(ncoh), stem_resp_coh(ncoh), root_resp_coh(ncoh))
         sum_lai = 0.0_wp
         do j = 1_ik, ncoh
            i = i0 + j - 1_ik
            coh%pft(j)       = site%cohort%pft(i)
            coh%nplant(j)    = site%cohort%nplant(i)
            coh%dbh(j)       = site%cohort%dbh(i)
            coh%height(j)    = site%cohort%height(i)
            coh%leaf_area(j) = site%cohort%leaf_area(i)
            coh%lai(j)       = site%cohort%nplant(i) * site%cohort%leaf_area(i)
            coh%bleaf(j)     = site%cohort%leaf_carbon(i)
            coh%broot(j)     = site%cohort%fineroot_carbon(i)
            !----- MVP derived geometry (proper allometry + PFT traits land with the RT/config step). !
            coh%wai(j)       = 0.20_wp * coh%lai(j)
            coh%bsap(j)      = 0.10_wp * site%cohort%wood_carbon(i)
            coh%sap_area(j)  = 0.05_wp * site%cohort%basal_area(i)
            sum_lai          = sum_lai + coh%lai(j)
         end do

         !----- Per-patch canopy geometry + constant forcing. -----------------------------!
         ageom%veg_height   = ctx%veg_height_bare
         do j = 1_ik, ncoh
            ageom%veg_height = max(ageom%veg_height, coh%height(j))
         end do
         ageom%opencan_frac = 0.0_wp ; ageom%snowfac = 0.0_wp

         call build_forcing(forc, coh, ctx, sum_lai)

         !----- Assemble the working bundle: adopt the owned per-patch reservoirs AND the        !
         !      PERSISTED per-cohort leaf_temp/psi carried on the cohort block (no reseeding).    !
         call alloc_patch_biophys(bio, ncoh, ctx%air_temp, ctx%shv_atm, ctx%co2_atm, ctx%air_temp)
         bio%cas    = site%patch%cas(ip)
         bio%soil_e = site%patch%soil_e(ip)
         bio%soil_w = site%patch%soil_w(ip)
         do j = 1_ik, ncoh
            i = i0 + j - 1_ik
            bio%leaf_temp(j) = site%cohort%leaf_temp(i)
            bio%psi(:,j)     = site%cohort%psi(:,i)
         end do

         call alloc_aero_out(aero, ncoh)
         call fill_aenv(aenv, bio, ctx)

         !----- n_fast_per_slow operator-split sweeps (budgets accumulate across the sweeps). !
         budg = column_budget_t()
         do isub = 1_ik, cfg%n_fast_per_slow
            call column_fast_step(cfg%dt_fast, cfg, ctx%ccfg, aenv, ageom, coh, forc, bio, aero, budg, &
                                  gpp_coh=gpp_coh, leaf_resp_coh=leaf_resp_coh,                        &
                                  stem_resp_coh=stem_resp_coh, root_resp_coh=root_resp_coh)
            !----- Integrate GROSS GPP + maintenance-resp losses [umol/plant/s] -> [kgC/plant].  !
            !      Keep gross and loss terms SEPARATE (carbon_growth nets them; mirrors ED2).     !
            do j = 1_ik, ncoh
               i = i0 + j - 1_ik
               site%cohort%gpp_accum(i)       = site%cohort%gpp_accum(i)       + gpp_coh(j)       * cfg%dt_fast * umol_2_kgC
               site%cohort%leaf_resp_accum(i) = site%cohort%leaf_resp_accum(i) + leaf_resp_coh(j) * cfg%dt_fast * umol_2_kgC
               site%cohort%stem_resp_accum(i) = site%cohort%stem_resp_accum(i) + stem_resp_coh(j) * cfg%dt_fast * umol_2_kgC
               site%cohort%root_resp_accum(i) = site%cohort%root_resp_accum(i) + root_resp_coh(j) * cfg%dt_fast * umol_2_kgC
            end do
            call fill_aenv(aenv, bio, ctx)          ! refresh CAS/ground state for the next sweep
         end do

         !----- Write the evolved state back to the site: per-patch reservoirs + per-cohort psi. !
         site%patch%cas(ip)    = bio%cas
         site%patch%soil_e(ip) = bio%soil_e
         site%patch%soil_w(ip) = bio%soil_w
         do j = 1_ik, ncoh
            i = i0 + j - 1_ik
            site%cohort%leaf_temp(i) = bio%leaf_temp(j)
            site%cohort%psi(:,i)     = bio%psi(:,j)
         end do

         we    = max(we, budg%whole_energy%worst) ; ww = max(ww, budg%whole_water%worst)
         nfail = nfail + budg%whole_energy%n_fail + budg%whole_water%n_fail
      end do

      if (present(worst_energy))  worst_energy  = we
      if (present(worst_water))   worst_water   = ww
      if (present(n_budget_fail)) n_budget_fail = nfail
   end subroutine run_fast_biophysics

   !----- Build the per-patch prescribed forcing from the reference met (constant MVP). ------!
   subroutine build_forcing(forc, coh, ctx, sum_lai)
      type(column_forcing_t), intent(out) :: forc
      type(column_cohort_t),  intent(in)  :: coh
      type(fast_context_t),   intent(in)  :: ctx
      real(wp),               intent(in)  :: sum_lai
      integer(ik) :: j
      allocate(forc%abs_sw(coh%n), forc%abs_lw(coh%n))
      forc%enthalpy_atm  = cas_enthalpy_of_temp(ctx%air_temp, ctx%shv_atm)
      forc%shv_atm       = ctx%shv_atm
      forc%co2_atm       = ctx%co2_atm
      forc%abs_sw_ground = ctx%rad_sw_ground
      forc%abs_lw_ground = 0.0_wp
      forc%precip        = ctx%precip
      !----- Split the canopy-top shortwave across cohorts by LAI share (MVP; real per-cohort !
      !      absorbed PAR comes from canopy RT at the config/RT step).                         !
      do j = 1_ik, coh%n
         if (sum_lai > tiny_num) then
            forc%abs_sw(j) = ctx%rad_sw_top * coh%lai(j) / sum_lai
         else
            forc%abs_sw(j) = 0.0_wp
         end if
         forc%abs_lw(j) = 0.0_wp
      end do
   end subroutine build_forcing

   !----- Fill the aerodynamics env from the reference met + the patch's current CAS/ground. -!
   subroutine fill_aenv(aenv, bio, ctx)
      type(aero_env_t),     intent(inout) :: aenv
      type(patch_biophys_t), intent(in)   :: bio
      type(fast_context_t), intent(in)    :: ctx
      aenv%u_ref = ctx%u_ref ; aenv%zref = ctx%zref ; aenv%press = ctx%press ; aenv%rho_air = ctx%rho_air
      !----- The aerodynamics buoyancy uses POTENTIAL temperatures (theta*(1+0.61 q)); convert the  !
      !      reference-level ACTUAL temperature with the shallow-layer dry-adiabatic form theta =    !
      !      T + (g/cp)*z. The CAS is the near-surface reference (z~0). Approximation: ignores the   !
      !      displacement height (a proper met driver will use zref - displace).  ------------------!
      aenv%theta_atm = ctx%air_temp + (grav / cp_air) * ctx%zref
      aenv%shv_atm = ctx%shv_atm ; aenv%co2_atm = ctx%co2_atm
      aenv%can_theta = bio%cas%can_temp ; aenv%can_temp = bio%cas%can_temp
      aenv%can_shv   = bio%cas%can_shv  ; aenv%can_co2  = bio%cas%can_co2
      aenv%t_ground  = bio%soil_e%soil_temp(1)
   end subroutine fill_aenv

end module meds_fast_loop
