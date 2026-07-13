!==========================================================================================!
! test_column_dynamics -- integration test for the fast-loop coupling core: aerodynamics ->     !
! {LEAF GAS EXCHANGE (real GPP + stomata + Rd) + PLANT HYDRAULICS, ground balance} -> soil WATER  !
! column -> CAS three-twin update -> soil THERMAL column, with autotrophic (leaf+stem+root) +      !
! heterotrophic respiration assembling the CAS-CO2 NEE. Driven by a diurnal cycle + morning rain.   !
!   1. CONSERVATION: every fast step closes the CAS energy/water/CO2, soil thermal, soil water,     !
!      AND the WHOLE-COLUMN energy + water budgets (meds_budget_check n_fail == 0).                 !
!   2. PHYSICAL SANITY: CAS temp tracks a diurnal cycle; leaf warms; the soil-surface swing damps   !
!      with depth; soil moisture responds to rain; real photosynthesis draws the CAS CO2 down in    !
!      daylight; the leaf water potential is under tension and more negative at midday.             !
!   3. OPT-IN INTER-LAYER ADVECTION: re-running with cfg%advect_soil_heat = .true. (the hydrology    !
!      per-face Darcy flux advecting liquid enthalpy) STILL closes the whole-column budgets and      !
!      keeps soil temperatures bounded -- the moisture<->energy coupling conserves.                  !
!==========================================================================================!
program test_column_dynamics
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : rho_h2o
   use meds_config,              only : meds_config_t
   use meds_time,                only : meds_time_t, solar_cosz
   use meds_thermo,              only : cas_enthalpy_of_temp, temp_to_uext
   use meds_biophysics_types,    only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out,    &
                                        patch_biophys_t, alloc_patch_biophys, SOIL_RETENTION_VG
   use meds_soil_parameters,     only : build_soil_params
   use meds_soil_thermal,        only : build_soil_thermal
   use meds_column_dynamics,     only : column_config_t, column_cohort_t, column_forcing_t,     &
                                        column_budget_t, alloc_column_cohort, column_fast_step,  &
                                        aero_bottom_to_top
   use meds_plant_interface,     only : NODE_LEAF
   use meds_test_support,        only : build_test_config
   implicit none

   integer(ik), parameter :: n = 1_ik, nsl = 10_ik, nstep = 96_ik    ! 96 x 900 s = 24 h
   real(wp),    parameter :: dt_fast = 900.0_wp, lat = 40.0_wp, t0 = 288.0_wp, theta0 = 0.30_wp
   type(meds_config_t)    :: cfg
   type(column_config_t)  :: ccfg
   type(column_cohort_t)  :: coh
   type(aero_env_t)       :: aenv
   type(aero_geom_t)      :: ageom
   type(aero_out_t)       :: aero
   type(patch_biophys_t)  :: bio
   type(column_forcing_t) :: forc
   type(column_budget_t)  :: budg
   type(meds_time_t)      :: sim_date
   real(wp) :: ct_night, ct_noon, co2_night, co2_noon, tleaf_noon, tleaf_night
   real(wp) :: ss_min, ss_max, sd_min, sd_max, th_min, th_max, gpp_noon, nee_noon
   real(wp) :: psileaf_noon, psileaf_night
   integer(ik) :: nfail

   nfail = 0_ik
   sim_date = meds_time_t(2001_ik, 6_ik, 21_ik)

   !----- Full model config (PFT traits for leaf gas exchange) + geometry. ----------------!
   cfg = build_test_config()
   ageom%veg_height = 18.0_wp ; ageom%opencan_frac = 0.0_wp ; ageom%snowfac = 0.0_wp
   aenv%u_ref = 2.0_wp ; aenv%zref = 30.0_wp ; aenv%press = 101325.0_wp ; aenv%rho_air = 1.2_wp

   !----- Cohort SoA (single cohort): LAI = nplant*leaf_area = 0.3*10 = 3. -----------------!
   call alloc_column_cohort(coh, n)
   coh%pft(1) = 1_ik ; coh%lai(1) = 3.0_wp ; coh%wai(1) = 0.5_wp
   coh%height(1) = 16.0_wp ; coh%crown(1) = 0.9_wp
   coh%leaf_width(1) = 0.04_wp ; coh%branch_diam(1) = 0.02_wp
   coh%leaf_area(1) = 10.0_wp ; coh%nplant(1) = 0.3_wp ; coh%dbh(1) = 20.0_wp ; coh%broot(1) = 0.5_wp
   coh%bleaf(1) = 0.5_wp ; coh%bsap(1) = 5.0_wp ; coh%sap_area(1) = 0.01_wp     ! for hydraulic capacitance

   !----- Static column config: soil column + respiration parameters. ---------------------!
   call build_soil_params(nsl, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,           &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, ccfg%soil)
   call build_soil_thermal(nsl, 3.0_wp, 0.15_wp, 2.0e6_wp, ccfg%soil_thermal)
   ccfg%wood%is_woody = .true. ; ccfg%wood%stem_resp_factor25 = 0.06_wp ; ccfg%wood%agf_bs = 0.7_wp
   ccfg%root%root_resp_factor25 = 0.30_wp
   ccfg%co2%rh_k_base = 0.01_wp                        ! nonzero decomposition rate so Rh > 0
   ccfg%fast_soil_carbon = 5.0_wp

   !----- Plant hydraulics parameters (Xu-2016-style capacitance + vulnerability). --------!
   ccfg%hydro_p%leaf_pi0 = -1.5_wp ; ccfg%hydro_p%leaf_eps = 12.0_wp
   ccfg%hydro_p%leaf_af  = 0.30_wp ; ccfg%hydro_p%leaf_water_sat = 2.0_wp
   ccfg%hydro_p%wood_pi0 = -1.0_wp ; ccfg%hydro_p%wood_eps = 8.0_wp
   ccfg%hydro_p%wood_af  = 0.20_wp ; ccfg%hydro_p%wood_water_sat = 1.0_wp
   ccfg%hydro_p%wood_psi50 = -2.0_wp ; ccfg%hydro_p%wood_kexp = 2.0_wp
   ccfg%hydro_p%k_plant_max = 6.0e-4_wp ; ccfg%hydro_p%wood_kmax = 8.0_wp
   ccfg%hydro_p%vessel_curl = 1.5_wp
   ccfg%rhizo_cond = 5.0e-4_wp                         ! soil->root rhizosphere conductance

   call alloc_aero_out(aero, n)
   allocate(forc%abs_sw(n), forc%abs_lw(n), forc%abs_par(n), forc%abs_sw_wood(n), forc%abs_lw_wood(n))
   forc%abs_sw_wood = 0.0_wp ; forc%abs_lw_wood = 0.0_wp

   !=====================================================================================!
   !  RUN 1 -- default coupling (advect_soil_heat = .false.): the full physical-sanity suite. !
   !=====================================================================================!
   call integrate_day(.false.)

   !----- 1. Conservation: all seven budgets closed every step. ----------------------------!
   call ck(budg%cas_energy%n_fail    == 0_ik, 'CAS energy budget closed',   real(budg%cas_energy%n_fail, wp))
   call ck(budg%cas_water%n_fail     == 0_ik, 'CAS water budget closed',    real(budg%cas_water%n_fail, wp))
   call ck(budg%cas_co2%n_fail       == 0_ik, 'CAS CO2 budget closed',      real(budg%cas_co2%n_fail, wp))
   call ck(budg%soil_energy%n_fail   == 0_ik, 'soil thermal budget closed', real(budg%soil_energy%n_fail, wp))
   call ck(budg%soil_water%n_fail    == 0_ik, 'soil water budget closed',   real(budg%soil_water%n_fail, wp))
   call ck(budg%whole_water%n_fail   == 0_ik, 'WHOLE-COLUMN water budget closes',  real(budg%whole_water%n_fail, wp))
   call ck(budg%whole_energy%n_fail  == 0_ik, 'WHOLE-COLUMN energy budget closes', real(budg%whole_energy%n_fail, wp))

   !----- 2. Physical sanity. -------------------------------------------------------------!
   call ck(ct_noon > ct_night, 'CAS warmer near solar noon than at night', ct_noon - ct_night)
   call ck(tleaf_noon > tleaf_night, 'leaf warms under absorbed shortwave', tleaf_noon - tleaf_night)
   call ck((ss_max - ss_min) > (sd_max - sd_min), 'soil diurnal swing damped with depth', &
           (ss_max - ss_min) - (sd_max - sd_min))
   call ck(th_max - th_min > 1.0e-4_wp, 'soil moisture responds to the rain pulse', th_max - th_min)
   call ck(th_min > 0.05_wp .and. th_max < 0.43_wp, 'soil moisture stays physical', th_max)
   call ck(gpp_noon > 1.0_wp, 'daytime GPP is active (real photosynthesis)', gpp_noon)
   call ck(nee_noon < 0.0_wp, 'daytime NEE is net uptake (GPP > respiration)', nee_noon)
   call ck(co2_noon < co2_night, 'CAS CO2 lower at midday than at night', co2_night - co2_noon)
   call ck(psileaf_noon < 0.0_wp, 'leaf water potential under tension in daylight', psileaf_noon)
   call ck(psileaf_noon < psileaf_night, 'leaf more tensioned at midday than at night',  &
           psileaf_noon - psileaf_night)

   if (nfail == 0_ik) then
      print '(a)', 'test_column_dynamics: RUN 1 (advect_soil_heat=F) PASSED'
      print '(a,f7.2,a,f7.2,a)', '   (CAS temp night=', ct_night, ' K  noon=', ct_noon, ' K)'
      print '(a,f7.2,a,f7.2,a)', '   (CAS CO2  night=', co2_night, '     noon=', co2_noon, ' umol/mol)'
      print '(a,f7.2,a,f7.2,a)', '   (noon GPP=', gpp_noon, ' umol/m2/s  NEE=', nee_noon, ' umol/m2/s)'
      print '(a,f7.3,a,f7.3,a)', '   (leaf psi night=', psileaf_night, ' MPa  noon=', psileaf_noon, ' MPa)'
      print '(a,es10.3,a,es10.3,a)', '   (whole-column worst resid: energy=', budg%whole_energy%worst,      &
                                     ' J/m2  water=', budg%whole_water%worst, ' kg/m2)'
   end if

   !=====================================================================================!
   !  RUN 2 -- opt-in inter-layer advection (advect_soil_heat = .true.): conserve + bounded.  !
   !=====================================================================================!
   call integrate_day(.true.)
   call ck(budg%whole_energy%n_fail == 0_ik, 'ADVECT: whole-column energy still closes', &
           real(budg%whole_energy%n_fail, wp))
   call ck(budg%whole_water%n_fail  == 0_ik, 'ADVECT: whole-column water still closes',  &
           real(budg%whole_water%n_fail, wp))
   call ck(ss_min > 270.0_wp .and. ss_max < 330.0_wp, 'ADVECT: soil surface temp stays bounded', ss_max)
   call ck(sd_min > 270.0_wp .and. sd_max < 330.0_wp, 'ADVECT: deep soil temp stays bounded',    sd_max)
   call ck(gpp_noon > 1.0_wp, 'ADVECT: daytime GPP still active', gpp_noon)

   !=====================================================================================!
   !  RUN 3 -- caller-side cohort ORDER: aero_bottom_to_top must respect the wind cascade.  !
   !=====================================================================================!
   call test_aero_order()

   if (nfail == 0_ik) then
      print '(a)', 'test_column_dynamics: RUN 2 (advect_soil_heat=T) PASSED'
      print '(a,f7.2,a,f7.2,a)', '   (CAS noon=', ct_noon, ' K  soil surf max=', ss_max, ' K)'
      print '(a,es10.3,a,es10.3,a)', '   (whole-column worst resid: energy=', budg%whole_energy%worst,       &
                                     ' J/m2  water=', budg%whole_water%worst, ' kg/m2)'
      print '(a)', 'test_column_dynamics: ALL PASSED'
   else
      print '(a,i0,a)', 'test_column_dynamics: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   !----- One 24 h diurnal integration from a freshly seeded column state. Fills the host    !
   !      diagnostic variables (min/max ranges, noon/night captures) + budg. ----------------!
   subroutine integrate_day(advect_heat)
      logical, intent(in) :: advect_heat
      real(wp)    :: t_sec, cosz, t_air
      integer(ik) :: istep, k

      ccfg%advect_soil_heat = advect_heat

      !----- (Re)seed the prognostic column state. --------------------------------------!
      if (allocated(bio%leaf_temp)) deallocate(bio%leaf_temp)
      if (allocated(bio%psi))       deallocate(bio%psi)
      call alloc_patch_biophys(bio, n, t0, 0.008_wp, 400.0_wp, t0)
      budg = column_budget_t()
      bio%soil_w%theta(1:nsl) = theta0
      do k = 1_ik, nsl
         bio%soil_e%soil_energy(k) = temp_to_uext(ccfg%soil_thermal%soil_dry_heat_capacity(k),  &
                                     theta0 * rho_h2o, t0, 1.0_wp)
         bio%soil_e%soil_temp(k)   = t0
      end do

      ss_min = 1.0e9_wp ; ss_max = -1.0e9_wp ; sd_min = 1.0e9_wp ; sd_max = -1.0e9_wp
      th_min = 1.0e9_wp ; th_max = -1.0e9_wp

      do istep = 1_ik, nstep
         t_sec = (real(istep, wp) - 0.5_wp) * dt_fast
         cosz  = solar_cosz(sim_date, t_sec, lat)
         t_air = 288.0_wp + 6.0_wp * (cosz - 0.3_wp)

         forc%abs_sw   = 500.0_wp * cosz                            ! leaf-absorbed shortwave [W/m2]
         forc%abs_par  = forc%abs_sw                                 ! MVP: absorbed PAR == absorbed SW (par_per_w=2.1)
         forc%abs_lw   = 0.0_wp
         forc%abs_sw_ground = 75.0_wp * cosz
         forc%abs_lw_ground = 0.0_wp
         forc%precip   = 0.0_wp
         if (istep >= 12_ik .and. istep <= 28_ik) forc%precip = 1.5e-4_wp     ! morning rain pulse
         forc%enthalpy_atm = cas_enthalpy_of_temp(t_air, 0.008_wp)
         forc%shv_atm      = 0.008_wp
         forc%co2_atm      = 400.0_wp

         call column_fast_step(dt_fast, cfg, ccfg, aenv, ageom, coh, forc, bio, aero, budg)

         ss_min = min(ss_min, bio%soil_e%soil_temp(1))   ; ss_max = max(ss_max, bio%soil_e%soil_temp(1))
         sd_min = min(sd_min, bio%soil_e%soil_temp(nsl)) ; sd_max = max(sd_max, bio%soil_e%soil_temp(nsl))
         th_min = min(th_min, bio%soil_w%theta(1))       ; th_max = max(th_max, bio%soil_w%theta(1))
         if (istep == 54_ik) then
            ct_noon = bio%cas%can_temp ; tleaf_noon = bio%leaf_temp(1) ; co2_noon = bio%cas%can_co2
            gpp_noon = budg%gpp_last ; nee_noon = budg%nee_last ; psileaf_noon = bio%psi(NODE_LEAF, 1)
         end if
         if (istep == 2_ik) then
            ct_night = bio%cas%can_temp ; tleaf_night = bio%leaf_temp(1) ; co2_night = bio%cas%can_co2
            psileaf_night = bio%psi(NODE_LEAF, 1)
         end if
      end do
   end subroutine integrate_day

   subroutine ck(cond, name, val)
      logical,          intent(in) :: cond
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: val
      if (.not. cond) then
         print '(a,a,a,es14.6)', '  FAIL ', name, ': val = ', val
         nfail = nfail + 1_ik
      end if
   end subroutine ck

   !----- Guard the caller-side cohort-order fix. aero_bottom_to_top feeds canopy_aerodynamics the !
   !      BOTTOM(1)->TOP(n) order it contracts for, starting from a height-DESCENDING column buffer  !
   !      (index 1 = tallest). So the TALL cohort (gather index 1 = canopy top) must come back with   !
   !      the top-of-canopy wind (=> higher boundary-layer gb); the short understory gets the         !
   !      attenuated wind. The old direct call (gather order straight into the kernel) inverts this,   !
   !      and reverting the fix must fail here.                                                        !
   subroutine test_aero_order()
      type(column_cohort_t) :: c2
      type(aero_out_t)      :: a2
      type(aero_env_t)      :: e2      ! defaults are fine (u_ref/zref/press/rho_air/can_* all sensible)
      type(aero_geom_t)     :: g2
      real(wp)              :: lt(2)
      call alloc_column_cohort(c2, 2_ik)
      c2%height = [18.0_wp, 6.0_wp] ; c2%lai = [2.0_wp, 2.0_wp] ; c2%crown = [0.9_wp, 0.8_wp]   ! DESCENDING
      c2%leaf_width = [0.04_wp, 0.04_wp] ; c2%branch_diam = [0.02_wp, 0.02_wp]
      lt = [300.0_wp, 300.0_wp]
      e2%u_ref = 3.0_wp
      g2%veg_height = 18.0_wp ; g2%opencan_frac = 0.0_wp ; g2%snowfac = 0.0_wp
      call alloc_aero_out(a2, 2_ik)
      call aero_bottom_to_top(ccfg%aero, e2, g2, 2_ik, c2, lt, a2)
      call ck(a2%wind(1) > a2%wind(2), 'aero order: tall cohort (gather idx1=top) gets more wind',    &
              a2%wind(1) - a2%wind(2))
      call ck(a2%leaf_gbw(1) > a2%leaf_gbw(2), 'aero order: tall cohort gets higher leaf gb',         &
              a2%leaf_gbw(1) - a2%leaf_gbw(2))
   end subroutine test_aero_order

end program test_column_dynamics
