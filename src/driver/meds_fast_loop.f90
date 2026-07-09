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
   use meds_time,             only : meds_time_t, time_advance_seconds
   use meds_forcing_types,    only : met_driver_t, met_forcing_t
   use meds_met_driver,       only : met_advance, met_instant
   use meds_demography_types, only : site_t
   use meds_biophysics_types, only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out,       &
                                     patch_biophys_t, alloc_patch_biophys, SOIL_RETENTION_VG,    &
                                     rad_pft_optics_t, rad_forcing_t, rad_flux_t,                &
                                     alloc_rad_forcing, N_RAD_BAND_DEFAULT, RAD_VIS, RAD_NIR, RAD_LW
   use meds_optics,           only : derive_rad_optics, ground_optics, surface_state_t,          &
                                     beta_params_from_mean
   use meds_canopy_radiation, only : canopy_radiation
   use meds_soil_parameters,  only : build_soil_params
   use meds_soil_thermal,     only : build_soil_thermal
   use meds_column_dynamics,  only : column_config_t, column_cohort_t, column_forcing_t,        &
                                     column_budget_t, alloc_column_cohort, column_fast_step
   implicit none
   private

   public :: fast_context_t, init_fast_reservoirs, run_fast_biophysics, build_fast_context

   !----- Absorbed-PAR (VIS) energy -> photon-flux conversion [umol photon / J], the 400-700 nm !
   !      value (~4.57). Used ONLY on the RT path (true absorbed PAR); the const path keeps the    !
   !      2.1 total-SW blend (column_forcing_t default) -- see run_fast_biophysics.                 !
   real(wp), parameter :: PAR_W_2_UMOL = 4.6_wp

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
      !----- Canopy-RT optics (per-PFT spectral/angle table + ground surface), for the RT-driven   !
      !      per-cohort absorbed SW/PAR on the forcing path (§6.3). MVP placeholders, built once.    !
      type(rad_pft_optics_t) :: rad_opt             !< per-PFT, per-band leaf/wood scattering + LIDF
      real(wp) :: soil_albedo(3) = [0.15_wp, 0.30_wp, 0.0_wp]  !< [VIS,NIR,LW] ground albedo (LW=0)
      real(wp) :: soil_emiss     = 0.95_wp          !< [-] ground longwave emissivity
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
      !----- P3 coupled-surface (Picard) solver knobs + option selectors, from the [fast] block. --!
      ctx%ccfg%picard_max_iter    = cfg%picard_max_iter
      ctx%ccfg%picard_tol_temp    = cfg%picard_tol_temp
      ctx%ccfg%picard_tol_shv     = cfg%picard_tol_shv
      ctx%ccfg%picard_relax       = cfg%picard_relax
      ctx%ccfg%picard_fixed_iter  = cfg%picard_fixed_iter
      ctx%ccfg%leaf_energy_model  = cfg%leaf_energy_model
      ctx%ccfg%soil_water_coupling = cfg%soil_water_coupling

      !----- Canopy-RT optics table (MVP placeholders; PFT-UNIFORM -- optics do not vary by PFT   !
      !      yet, that is the Phase-2 [radiation] PFT-TOML block). Values mirror                    !
      !      test_canopy_radiation.f90 build_optics. Built ONCE; read-only downstream.               !
      block
         integer(ik), parameter :: NB = N_RAD_BAND_DEFAULT
         integer(ik) :: np
         real(wp), allocatable :: rl(:,:), tl(:,:), rw(:,:), tw(:,:), cl(:), cw(:), bp(:), bq(:)
         logical  :: hb(NB), he(NB)
         real(wp) :: bpp, bqq
         np = cfg%pft%n
         allocate(rl(NB,np), tl(NB,np), rw(NB,np), tw(NB,np), cl(np), cw(np), bp(np), bq(np))
         rl(RAD_VIS,:) = 0.10_wp ; tl(RAD_VIS,:) = 0.05_wp        ! leaf VIS reflect/transmit
         rl(RAD_NIR,:) = 0.45_wp ; tl(RAD_NIR,:) = 0.25_wp        ! leaf NIR
         rl(RAD_LW,:)  = 0.03_wp ; tl(RAD_LW,:)  = 0.0_wp         ! leaf_emiss = 0.97
         rw(RAD_VIS,:) = 0.11_wp ; tw(RAD_VIS,:) = 0.001_wp       ! wood VIS (near-opaque)
         rw(RAD_NIR,:) = 0.25_wp ; tw(RAD_NIR,:) = 0.001_wp       ! wood NIR
         rw(RAD_LW,:)  = 0.10_wp ; tw(RAD_LW,:)  = 0.0_wp         ! wood_emiss = 0.90
         cl = 0.80_wp ; cw = 0.50_wp                              ! leaf/wood clumping
         call beta_params_from_mean(45.0_wp, 20.0_wp, bpp, bqq)   ! mean 45deg, std 20deg leaf angle
         bp = bpp ; bq = bqq
         hb = [.true.,  .true.,  .false.]                         ! VIS/NIR have a beam; LW does not
         he = [.false., .false., .true. ]                         ! only LW emits
         call derive_rad_optics(NB, np, rl, tl, rw, tw, cl, cw, bp, bq, hb, he, ctx%rad_opt)
      end block
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
   subroutine run_fast_biophysics(site, ctx, cfg, met_drv, step_start, worst_energy, worst_water, n_budget_fail)
      type(site_t),         intent(inout) :: site
      type(fast_context_t), intent(in)    :: ctx
      type(meds_config_t),  intent(in)    :: cfg
      type(met_driver_t), optional, intent(inout) :: met_drv       !< live forcing reader (per-sub-step met)
      type(meds_time_t),  optional, intent(in)    :: step_start    !< calendar time at the START of this slow step
      real(wp),    optional, intent(out)  :: worst_energy, worst_water
      integer(ik), optional, intent(out)  :: n_budget_fail

      type(column_cohort_t)  :: coh
      type(column_forcing_t) :: forc
      type(aero_env_t)       :: aenv
      type(aero_geom_t)      :: ageom
      type(aero_out_t)       :: aero
      type(patch_biophys_t)  :: bio
      type(column_budget_t)  :: budg
      type(fast_context_t)   :: ctx_now                            !< per-sub-step met overlay on ctx
      type(met_forcing_t)    :: met
      type(meds_time_t)      :: t_sub
      real(wp), allocatable  :: gpp_coh(:), leaf_resp_coh(:), stem_resp_coh(:), root_resp_coh(:)
      real(wp)    :: we, ww, sum_lai, f_ground
      integer(ik) :: ip, isub, j, i, i0, ncoh, nfail
      logical     :: do_forcing

      !----- Live forcing drives the fast loop only when it is ON and a reader + step time are    !
      !      supplied; otherwise ctx_now stays == ctx and the loop runs the CONSTANT-forcing MVP    !
      !      bit-identically (the diurnal cycle lives INSIDE the sub-step loop, design §1.1/§6.2).  !
      do_forcing = cfg%forcing%forcing_on .and. present(met_drv) .and. present(step_start)
      ctx_now  = ctx
      f_ground = ctx%rad_sw_ground / max(ctx%rad_sw_top, tiny_num)   ! ground/canopy-top SW transmittance

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

         call alloc_forcing(forc, ncoh)

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

         !----- n_fast_per_slow operator-split sweeps. Forcing is re-evaluated PER SUB-STEP (the   !
         !      diurnal cycle lives here): refresh the met overlay ctx_now, then fill_forcing +     !
         !      fill_aenv from it. CONSTANT path (do_forcing=.false.): ctx_now==ctx, so these        !
         !      reproduce the old build_forcing-once + fill_aenv sequence bit-identically.           !
         budg = column_budget_t()
         do isub = 1_ik, cfg%n_fast_per_slow
            if (do_forcing) then
               t_sub = time_advance_seconds(step_start, (real(isub, wp) - 0.5_wp) * cfg%dt_fast)  ! sub-interval midpoint
               call met_advance(met_drv, t_sub)
               met = met_instant(met_drv, t_sub)
               call apply_met_to_ctx(ctx_now, met, f_ground)
            end if
            call fill_forcing(forc, coh, ctx_now, sum_lai)
            !----- RT join (§6.3): when forcing is on, REPLACE the LAI-share SW split with real     !
            !      per-cohort absorbed SW/PAR from the two-stream canopy radiation (ctx%rad_opt read !
            !      directly -- not the ctx_now overlay -- so the allocatable table is not deep-copied). !
            if (do_forcing) call apply_rt_forcing(forc, coh, bio, ctx, met, cfg)
            call fill_aenv(aenv, bio, ctx_now)
            call column_fast_step(cfg%dt_fast, cfg, ctx_now%ccfg, aenv, ageom, coh, forc, bio, aero, budg, &
                                  gpp_coh=gpp_coh, leaf_resp_coh=leaf_resp_coh,                            &
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

   !----- Allocate the per-patch forcing buffers ONCE; fill_forcing then updates VALUES each   !
   !      sub-step (so time-varying met can drive them without re-allocating). ----------------!
   subroutine alloc_forcing(forc, ncoh)
      type(column_forcing_t), intent(inout) :: forc
      integer(ik),            intent(in)    :: ncoh
      !----- Resize (not just allocate-once): `forc` is reused across the patch loop, and a later    !
      !      patch can hold MORE cohorts than the first. The sibling per-patch buffers (coh/bio/aero) !
      !      resize via intent(out); the three per-cohort forcing arrays must match, or fill_forcing  !
      !      / apply_rt_forcing would write out of bounds on a larger downstream patch.               !
      if (allocated(forc%abs_sw)) then
         if (size(forc%abs_sw) /= ncoh) deallocate(forc%abs_sw, forc%abs_lw, forc%abs_par)
      end if
      if (.not. allocated(forc%abs_sw)) allocate(forc%abs_sw(ncoh), forc%abs_lw(ncoh), forc%abs_par(ncoh))
   end subroutine alloc_forcing

   !----- Fill the per-patch prescribed forcing from the (possibly per-sub-step) reference met. !
   subroutine fill_forcing(forc, coh, ctx, sum_lai)
      type(column_forcing_t), intent(inout) :: forc
      type(column_cohort_t),  intent(in)    :: coh
      type(fast_context_t),   intent(in)    :: ctx
      real(wp),               intent(in)    :: sum_lai
      integer(ik) :: j
      forc%enthalpy_atm  = cas_enthalpy_of_temp(ctx%air_temp, ctx%shv_atm)
      forc%shv_atm       = ctx%shv_atm
      forc%co2_atm       = ctx%co2_atm
      forc%abs_sw_ground = ctx%rad_sw_ground
      forc%abs_lw_ground = 0.0_wp
      forc%precip        = ctx%precip
      forc%par_per_w     = 2.1_wp                    ! LAI-split path: total-SW->PAR blend (abs_par == abs_sw)
      !----- Split the canopy-top shortwave across cohorts by LAI share (MVP; the RT join (§6.3) !
      !      replaces this with real per-cohort absorbed SW/PAR when forcing is on).             !
      do j = 1_ik, coh%n
         if (sum_lai > tiny_num) then
            forc%abs_sw(j) = ctx%rad_sw_top * coh%lai(j) / sum_lai
         else
            forc%abs_sw(j) = 0.0_wp
         end if
         forc%abs_par(j) = forc%abs_sw(j)            ! no PAR/NIR split in the LAI path -> PAR==SW (biased high)
         forc%abs_lw(j)  = 0.0_wp
      end do
   end subroutine fill_forcing

   !----- The met-source shim (design §6.2/§6.5, retired at P1): copy the instantaneous          !
   !      met_forcing_t's raw scalars into fast_context_t's met fields, so fill_forcing/fill_aenv  !
   !      keep their present logic unchanged. cosz/leaf_temp/ustar/can_co2 are NOT overwritten     !
   !      (they are prognostic or derived). Ground SW scales with canopy-top SW at the reference    !
   !      transmittance f_ground (reproduces the ad-hoc rad_sw_ground=60 at swdown=400).            !
   subroutine apply_met_to_ctx(ctx, met, f_ground)
      type(fast_context_t), intent(inout) :: ctx
      type(met_forcing_t),  intent(in)    :: met
      real(wp),             intent(in)    :: f_ground
      ctx%air_temp      = met%tair_k
      ctx%shv_atm       = met%qair
      ctx%press         = met%psurf_pa
      ctx%rho_air       = met%rho_air
      ctx%co2_atm       = met%co2
      ctx%u_ref         = met%wind
      ctx%precip        = met%rainf
      ctx%rad_sw_top    = met%swdown()
      ctx%rad_sw_ground = f_ground * met%swdown()
   end subroutine apply_met_to_ctx

   !----- RT join (§6.3): run the two-stream canopy radiation for this patch and OVERWRITE the      !
   !      LAI-share SW split in `forc` with real per-cohort absorbed SW (leaf energy) + PAR           !
   !      (photosynthesis) + NET longwave (leaf + ground) + below-canopy ground SW.                   !
   !      Cohorts are gathered height-DESCENDING (top=1) but the two-stream wants BOTTOM(1)->TOP, so  !
   !      we build an ascending-height permutation `perm` and inverse-scatter the outputs.            !
   !      NOTE on PAR: the leaf model's `env%par` is INCIDENT PAR (it re-applies cfg%leaf_absorptance !
   !      internally for electron transport), so we divide the two-stream ABSORBED VIS by that same   !
   !      absorptance to hand back an incident-equivalent PAR -- otherwise leaf absorptance would be  !
   !      applied twice. abs_sw stays true ABSORBED SW (the leaf energy balance wants absorbed).      !
   subroutine apply_rt_forcing(forc, coh, bio, ctx, met, cfg)
      type(column_forcing_t), intent(inout) :: forc
      type(column_cohort_t),  intent(in)    :: coh
      type(patch_biophys_t),  intent(in)    :: bio
      type(fast_context_t),   intent(in)    :: ctx
      type(met_forcing_t),    intent(in)    :: met
      type(meds_config_t),    intent(in)    :: cfg
      integer(ik) :: ncoh, j, k, ig, imin
      integer(ik) :: perm(coh%n), pft_bt(coh%n)
      real(wp)    :: lai_bt(coh%n), wai_bt(coh%n), tcan_bt(coh%n)
      logical     :: used(coh%n)
      real(wp)    :: hmin, tcas
      type(rad_forcing_t)   :: rf
      type(rad_flux_t)      :: flux
      type(surface_state_t) :: surf
      logical :: he(N_RAD_BAND_DEFAULT)

      ncoh = coh%n
      !----- A bare patch (ncoh == 0) is NOT special-cased: the zero-trip perm/scatter loops fall     !
      !      through and canopy_radiation's own empty-canopy branch returns the correct NET ground SW  !
      !      (incident * (1 - soil albedo)), so a patch shedding its last cohort stays continuous.     !

      !----- LW emission base = the CAS temperature. The diagnostic leaf energy balance linearizes    !
      !      leaf LW emission around tcas (lw_slope*dtl, dtl=tl-tcas), so it needs abs_lw = NET LW AT   !
      !      tcas; feeding the two-stream tcas as the canopy emission temp makes abs_leaf(LW) exactly   !
      !      that. (Prognostic CAS enthalpy is always valid here, unlike the lagged bio%cas%can_temp.)  !
      tcas = cas_temp_of_enthalpy(bio%cas%can_enthalpy, bio%cas%can_shv)

      !----- perm: gather-indices in ASCENDING height (bottom -> top). Selection sort (ncoh small). !
      used = .false.
      do j = 1_ik, ncoh
         imin = 0_ik ; hmin = huge(1.0_wp)
         do k = 1_ik, ncoh
            if (.not. used(k) .and. coh%height(k) <= hmin) then ; hmin = coh%height(k) ; imin = k ; end if
         end do
         perm(j) = imin ; used(imin) = .true.
         pft_bt(j) = coh%pft(imin) ; lai_bt(j) = coh%lai(imin)
         !----- LW emission temperature: the diagnostic (split) leaf balance linearizes around tcas, !
         !      so feed tcas; the PICARD balance re-bases emission to the cohort's own leaf_temp (P3c), !
         !      so feed leaf_temp -- the leaf then emits LW at leaf_temp, consistent within the sub-step.!
         wai_bt(j) = coh%wai(imin)
         if (cfg%integration_scheme == SCHEME_PICARD_COUPLED) then
            tcan_bt(j) = bio%leaf_temp(imin)
         else
            tcan_bt(j) = tcas
         end if
      end do

      !----- rad_forcing_t from met (§6.3 mapping table; all W/m2, direct assignment). -----------!
      call alloc_rad_forcing(rf, N_RAD_BAND_DEFAULT)
      rf%cosz = met%cosz
      rf%incid_beam(RAD_VIS) = met%par_beam ; rf%incid_diff(RAD_VIS) = met%par_diffuse
      rf%incid_beam(RAD_NIR) = met%nir_beam ; rf%incid_diff(RAD_NIR) = met%nir_diffuse
      rf%incid_beam(RAD_LW)  = 0.0_wp       ; rf%incid_diff(RAD_LW)  = met%lwdown

      !----- ground optics (MVP soil albedo/emiss from ctx; ground skin temp from the soil column). !
      surf%n_band = N_RAD_BAND_DEFAULT
      allocate(surf%soil_albedo(N_RAD_BAND_DEFAULT))
      surf%soil_albedo = ctx%soil_albedo ; surf%soil_emiss = ctx%soil_emiss
      surf%soil_temp   = bio%soil_e%soil_temp(1)
      he = [.false., .false., .true.]
      call ground_optics(surf, N_RAD_BAND_DEFAULT, he, rf%grnd_refl, rf%grnd_emiss)

      call canopy_radiation(ctx%rad_opt, rf, ncoh, pft_bt, lai_bt, wai_bt, tcan_bt, flux)

      !----- inverse-scatter: RT index j (bottom->top) maps to gather index perm(j). --------------!
      do j = 1_ik, ncoh
         ig = perm(j)
         forc%abs_sw(ig)  = flux%abs_leaf(RAD_VIS, j) + flux%abs_leaf(RAD_NIR, j)   ! total ABSORBED leaf SW (energy)
         forc%abs_par(ig) = flux%abs_leaf(RAD_VIS, j) / max(cfg%leaf_absorptance, tiny_num)  ! -> INCIDENT-equiv PAR
         forc%abs_lw(ig)  = flux%abs_leaf(RAD_LW, j)                                ! NET leaf LW at tcas (emission incl.)
      end do
      forc%abs_sw_ground = (flux%dn_ground(RAD_VIS) - flux%up_ground(RAD_VIS))                      &
                         + (flux%dn_ground(RAD_NIR) - flux%up_ground(RAD_NIR))
      forc%abs_lw_ground = flux%dn_ground(RAD_LW) - flux%up_ground(RAD_LW)          ! NET ground LW (soil emission incl.)
      forc%par_per_w     = PAR_W_2_UMOL                        ! true VIS absorbed -> photon flux
   end subroutine apply_rt_forcing

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
