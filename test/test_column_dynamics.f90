!==========================================================================================!
! test_column_dynamics -- integration test for the fast-loop coupling core: aerodynamics ->     !
! {leaf balance, ground balance} -> soil WATER column -> canopy-air-space (CAS) three-twin       !
! update -> soil THERMAL column, driven by a diurnal cycle with a morning rain pulse.            !
!   1. CONSERVATION: every fast step closes the CAS energy / water / CO2 budgets AND the soil     !
!      thermal AND soil water budgets to machine precision (meds_budget_check n_fail == 0).        !
!   2. PHYSICAL SANITY: CAS tracks a diurnal cycle; leaf warms under shortwave; CAS CO2 responds  !
!      to the biotic flux; the soil-surface diurnal swing DAMPS with depth; and soil moisture      !
!      responds to the rain event while staying within [theta_res, theta_sat].                     !
!==========================================================================================!
program test_column_dynamics
   use meds_kinds,               only : wp, ik
   use meds_constants,           only : rho_h2o
   use meds_time,                only : meds_time_t, solar_cosz
   use meds_thermo,              only : cas_enthalpy_of_temp, temp_to_uext
   use meds_biophysics_types,    only : aero_env_t, aero_geom_t, aero_out_t, alloc_aero_out,    &
                                        patch_biophys_t, alloc_patch_biophys, SOIL_RETENTION_VG
   use meds_soil_parameters,     only : build_soil_params
   use meds_soil_thermal,        only : build_soil_thermal
   use meds_column_dynamics,     only : column_config_t, column_forcing_t, column_budget_t,     &
                                        column_fast_step
   implicit none

   integer(ik), parameter :: n = 1_ik, nsl = 10_ik, nstep = 96_ik    ! 96 x 900 s = 24 h
   real(wp),    parameter :: dt_fast = 900.0_wp, lat = 40.0_wp, t0 = 288.0_wp, theta0 = 0.30_wp
   type(column_config_t)  :: ccfg
   type(aero_env_t)       :: aenv
   type(aero_geom_t)      :: ageom
   type(aero_out_t)       :: aero
   type(patch_biophys_t)  :: bio
   type(column_forcing_t) :: forc
   type(column_budget_t)  :: budg
   type(meds_time_t)      :: sim_date
   real(wp) :: lai(n), wai(n), height(n), crown(n), lwid(n), bdia(n)
   real(wp) :: t_sec, cosz, t_air
   real(wp) :: ct_night, ct_noon, co2_min, co2_max, tleaf_noon, tleaf_night
   real(wp) :: ss_min, ss_max, sd_min, sd_max, th_min, th_max
   integer(ik) :: istep, k, nfail

   nfail = 0_ik
   sim_date = meds_time_t(2001_ik, 6_ik, 21_ik)

   ageom%veg_height = 18.0_wp ; ageom%opencan_frac = 0.0_wp ; ageom%snowfac = 0.0_wp
   aenv%u_ref = 2.0_wp ; aenv%zref = 30.0_wp ; aenv%press = 101325.0_wp ; aenv%rho_air = 1.2_wp
   lai(1) = 3.0_wp ; wai(1) = 0.5_wp ; height(1) = 16.0_wp ; crown(1) = 0.9_wp
   lwid(1) = 0.04_wp ; bdia(1) = 0.02_wp

   !----- Static column config: aero + veg-thermal + soil-opts defaults; build the soil column.!
   call build_soil_params(nsl, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,           &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, ccfg%soil)
   call build_soil_thermal(nsl, 3.0_wp, 0.15_wp, 2.0e6_wp, ccfg%soil_thermal)

   call alloc_aero_out(aero, n)
   call alloc_patch_biophys(bio, n, t0, 0.008_wp, 400.0_wp, t0)
   allocate(forc%abs_sw(n), forc%abs_lw(n), forc%gsw(n), forc%fs_open(n))

   !----- Seed the soil water + thermal columns. ------------------------------------------!
   bio%soil_w%theta(1:nsl) = theta0
   do k = 1_ik, nsl
      bio%soil_e%soil_energy(k) = temp_to_uext(ccfg%soil_thermal%soil_dry_heat_capacity(k),     &
                                  theta0 * rho_h2o, t0, 1.0_wp)
      bio%soil_e%soil_temp(k)   = t0
   end do

   co2_min = 1.0e9_wp ; co2_max = -1.0e9_wp
   ss_min = 1.0e9_wp ; ss_max = -1.0e9_wp ; sd_min = 1.0e9_wp ; sd_max = -1.0e9_wp
   th_min = 1.0e9_wp ; th_max = -1.0e9_wp

   do istep = 1_ik, nstep
      t_sec = (real(istep, wp) - 0.5_wp) * dt_fast
      cosz  = solar_cosz(sim_date, t_sec, lat)
      t_air = 288.0_wp + 6.0_wp * (cosz - 0.3_wp)

      forc%abs_sw   = 500.0_wp * cosz
      forc%abs_lw   = 0.0_wp
      forc%gsw      = 0.003_wp * cosz + 1.0e-4_wp
      forc%fs_open  = 1.0_wp
      forc%abs_sw_ground = 75.0_wp * cosz
      forc%abs_lw_ground = 0.0_wp
      forc%precip   = 0.0_wp
      if (istep >= 12_ik .and. istep <= 28_ik) forc%precip = 1.5e-4_wp     ! morning rain pulse
      forc%enthalpy_atm = cas_enthalpy_of_temp(t_air, 0.008_wp)
      forc%shv_atm      = 0.008_wp
      forc%co2_atm      = 400.0_wp
      forc%nee_biotic   = 4.0_wp - 12.0_wp * cosz

      call column_fast_step(dt_fast, ccfg, aenv, ageom, n, lai, wai, height, crown, lwid, bdia, &
                            forc, bio, aero, budg)

      co2_min = min(co2_min, bio%cas%can_co2) ; co2_max = max(co2_max, bio%cas%can_co2)
      ss_min  = min(ss_min, bio%soil_e%soil_temp(1))   ; ss_max = max(ss_max, bio%soil_e%soil_temp(1))
      sd_min  = min(sd_min, bio%soil_e%soil_temp(nsl)) ; sd_max = max(sd_max, bio%soil_e%soil_temp(nsl))
      th_min  = min(th_min, bio%soil_w%theta(1))       ; th_max = max(th_max, bio%soil_w%theta(1))
      if (istep == 54_ik) then ; ct_noon  = bio%cas%can_temp ; tleaf_noon  = bio%leaf_temp(1) ; end if
      if (istep ==  2_ik) then ; ct_night = bio%cas%can_temp ; tleaf_night = bio%leaf_temp(1) ; end if
   end do

   !----- 1. Conservation: all five budgets closed every step. ----------------------------!
   call ck(budg%cas_energy%n_fail  == 0_ik, 'CAS energy budget closed every step',   real(budg%cas_energy%n_fail, wp))
   call ck(budg%cas_water%n_fail   == 0_ik, 'CAS water  budget closed every step',   real(budg%cas_water%n_fail, wp))
   call ck(budg%cas_co2%n_fail     == 0_ik, 'CAS CO2    budget closed every step',   real(budg%cas_co2%n_fail, wp))
   call ck(budg%soil_energy%n_fail == 0_ik, 'soil thermal budget closed every step', real(budg%soil_energy%n_fail, wp))
   call ck(budg%soil_water%n_fail  == 0_ik, 'soil water budget closed every step',   real(budg%soil_water%n_fail, wp))

   !----- 2. Physical sanity. -------------------------------------------------------------!
   call ck(ct_noon > ct_night, 'CAS warmer near solar noon than at night', ct_noon - ct_night)
   call ck(tleaf_noon > tleaf_night, 'leaf warms under absorbed shortwave', tleaf_noon - tleaf_night)
   call ck(co2_max > co2_min + 1.0_wp, 'CAS CO2 responds to day/night biotic flux', co2_max - co2_min)
   call ck(ss_max > ss_min, 'soil surface warms during the day', ss_max - ss_min)
   call ck((ss_max - ss_min) > (sd_max - sd_min), 'diurnal swing damped with depth', &
           (ss_max - ss_min) - (sd_max - sd_min))
   call ck(th_max - th_min > 1.0e-4_wp, 'soil moisture responds to the rain pulse', th_max - th_min)
   call ck(th_min > 0.05_wp .and. th_max < 0.43_wp, 'soil moisture stays physical', th_max)

   if (nfail == 0_ik) then
      print '(a)', 'test_column_dynamics: ALL PASSED'
      print '(a,f7.2,a,f7.2,a)', '   (CAS night=', ct_night, ' K  noon=', ct_noon, ' K)'
      print '(a,f6.3,a,f6.3,a)', '   (soil swing: surface=', ss_max-ss_min, ' K  deep=', sd_max-sd_min, ' K)'
      print '(a,f6.4,a,f6.4,a)', '   (soil theta(1): min=', th_min, '  max=', th_max, ')'
   else
      print '(a,i0,a)', 'test_column_dynamics: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   subroutine ck(cond, name, val)
      logical,          intent(in) :: cond
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: val
      if (.not. cond) then
         print '(a,a,a,es14.6)', '  FAIL ', name, ': val = ', val
         nfail = nfail + 1_ik
      end if
   end subroutine ck

end program test_column_dynamics
