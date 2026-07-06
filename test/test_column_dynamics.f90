!==========================================================================================!
! test_column_dynamics -- integration test for the fast-loop coupling core: aerodynamics ->     !
! leaf/wood energy -> canopy-air-space (CAS) three-twin update, driven by a diurnal cycle.      !
!   1. CONSERVATION: every fast step closes the CAS energy / water / CO2 budgets to machine       !
!      precision (meds_budget_check n_fail == 0 across a full simulated day).                     !
!   2. PHYSICAL SANITY: CAS temperature stays physical and tracks a diurnal cycle (warmer near     !
!      solar noon than at night); the leaf warms under absorbed shortwave; CAS CO2 responds to     !
!      the day/night net biotic flux.                                                              !
!==========================================================================================!
program test_column_dynamics
   use meds_kinds,               only : wp, ik
   use meds_time,                only : meds_time_t, solar_cosz
   use meds_thermo,              only : cas_enthalpy_of_temp
   use meds_biophysics_types,    only : aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t,        &
                                        alloc_aero_out, veg_thermal_params_t
   use meds_budget_check,        only : budget_t
   use meds_column_dynamics,     only : patch_biophys_t, alloc_patch_biophys,                   &
                                        column_forcing_t, column_fast_step
   implicit none

   integer(ik), parameter :: n = 1_ik, nstep = 96_ik            ! 96 x 900 s = 24 h
   real(wp),    parameter :: dt_fast = 900.0_wp, lat = 40.0_wp
   type(aero_cfg_t)           :: acfg
   type(aero_env_t)           :: aenv
   type(aero_geom_t)          :: ageom
   type(aero_out_t)           :: aero
   type(veg_thermal_params_t) :: tparams
   type(patch_biophys_t)      :: bio
   type(column_forcing_t)     :: forc
   type(budget_t)             :: be_e, be_w, be_c
   type(meds_time_t)          :: sim_date
   real(wp) :: lai(n), wai(n), height(n), crown(n), lwid(n), bdia(n)
   real(wp) :: t_sec, cosz, t_air
   real(wp) :: ct_min, ct_max, ct_night, ct_noon, co2_min, co2_max, tleaf_noon, tleaf_night
   integer(ik) :: istep, nfail

   nfail = 0_ik
   sim_date  = meds_time_t(2001_ik, 6_ik, 21_ik)                     ! N summer solstice

   !----- Canopy + cohort setup (single cohort). -----------------------------------------!
   ageom%veg_height = 18.0_wp ; ageom%opencan_frac = 0.0_wp ; ageom%snowfac = 0.0_wp
   aenv%u_ref = 2.0_wp ; aenv%zref = 30.0_wp ; aenv%press = 101325.0_wp ; aenv%rho_air = 1.2_wp
   lai(1) = 3.0_wp ; wai(1) = 0.5_wp ; height(1) = 16.0_wp ; crown(1) = 0.9_wp
   lwid(1) = 0.04_wp ; bdia(1) = 0.02_wp

   call alloc_aero_out(aero, n)
   call alloc_patch_biophys(bio, n, 288.0_wp, 0.008_wp, 400.0_wp, 288.0_wp)
   allocate(forc%abs_sw(n), forc%abs_lw(n), forc%gsw(n), forc%fs_open(n))

   ct_min = 1.0e9_wp ; ct_max = -1.0e9_wp ; co2_min = 1.0e9_wp ; co2_max = -1.0e9_wp

   !----- Diurnal loop. ------------------------------------------------------------------!
   do istep = 1_ik, nstep
      t_sec = (real(istep, wp) - 0.5_wp) * dt_fast               ! mid-step, seconds into the day
      cosz  = solar_cosz(sim_date, t_sec, lat)
      t_air = 288.0_wp + 6.0_wp * (cosz - 0.3_wp)                ! diurnal air temperature

      forc%abs_sw   = 500.0_wp * cosz                            ! absorbed shortwave (leaf) [W/m2]
      forc%abs_lw   = 0.0_wp
      forc%gsw      = 0.003_wp * cosz + 1.0e-4_wp                ! stomata open in daylight [m/s]
      forc%fs_open  = 1.0_wp
      forc%enthalpy_atm = cas_enthalpy_of_temp(t_air, 0.008_wp)
      forc%shv_atm      = 0.008_wp
      forc%co2_atm      = 400.0_wp
      forc%nee_biotic   = 4.0_wp - 12.0_wp * cosz                ! Reco 4, GPP ~12*cosz [umol/m2/s]

      call column_fast_step(dt_fast, acfg, aenv, ageom, n, lai, wai, height, crown, lwid, bdia, &
                            tparams, forc, bio, aero, be_e, be_w, be_c)

      ct_min = min(ct_min, bio%cas%can_temp) ; ct_max = max(ct_max, bio%cas%can_temp)
      co2_min = min(co2_min, bio%cas%can_co2) ; co2_max = max(co2_max, bio%cas%can_co2)
      if (istep == 48_ik) then ; ct_noon  = bio%cas%can_temp ; tleaf_noon  = bio%leaf_temp(1) ; end if
      if (istep ==  2_ik) then ; ct_night = bio%cas%can_temp ; tleaf_night = bio%leaf_temp(1) ; end if
   end do

   !----- 1. Conservation: budgets closed every step. -------------------------------------!
   call ck(be_e%n_fail == 0_ik, 'CAS energy budget closed every step', real(be_e%n_fail, wp))
   call ck(be_w%n_fail == 0_ik, 'CAS water  budget closed every step', real(be_w%n_fail, wp))
   call ck(be_c%n_fail == 0_ik, 'CAS CO2    budget closed every step', real(be_c%n_fail, wp))

   !----- 2. Physical sanity. -------------------------------------------------------------!
   call ck(ct_min > 270.0_wp .and. ct_max < 330.0_wp, 'CAS temp stays physical', ct_max)
   call ck(ct_noon > ct_night, 'CAS warmer near solar noon than at night', ct_noon - ct_night)
   call ck(tleaf_noon > tleaf_night, 'leaf warms under absorbed shortwave', tleaf_noon - tleaf_night)
   call ck(co2_max > co2_min + 1.0_wp, 'CAS CO2 responds to day/night biotic flux', co2_max - co2_min)

   if (nfail == 0_ik) then
      print '(a)', 'test_column_dynamics: ALL PASSED'
      print '(a,f7.2,a,f7.2,a)', '   (CAS temp night=', ct_night, ' K  noon=', ct_noon, ' K)'
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
