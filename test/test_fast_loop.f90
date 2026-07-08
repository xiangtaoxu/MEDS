!==========================================================================================!
! test_fast_loop -- the SITE-LEVEL fast-biophysics driver over the state-hub-owned per-patch     !
! reservoirs. Builds a two-patch site (one vegetated, one bare ground), seeds the reservoirs,     !
! and runs the fast loop:                                                                        !
!   1. CONSERVATION: every patch's whole-column energy + water budgets close (n_fail == 0) --     !
!      the fast loop conserves over the reservoirs the state hub owns, for a VEGETATED and a      !
!      BARE (zero-cohort) patch alike.                                                            !
!   2. ACTIVITY: the reservoirs evolve (CAS warms under the absorbed shortwave; soil responds).   !
!   3. STEPPER WIRING + GATE: advance_one_step runs the fast loop only when fast_biophysics_on    !
!      AND a fast context is supplied; with the gate off the reservoirs are untouched.            !
!==========================================================================================!
program test_fast_loop
   use meds_kinds,               only : wp, ik
   use meds_config,              only : meds_config_t, GS_CARBON, GS_EMPIRICAL
   use meds_demography_types,    only : site_t
   use meds_init,                only : init_bare_ground, add_cohort, finalize_init
   use meds_soil_parameters,     only : build_soil_params
   use meds_soil_thermal,        only : build_soil_thermal
   use meds_biophysics_types,    only : SOIL_RETENTION_VG
   use meds_fast_loop,           only : fast_context_t, init_fast_reservoirs, run_fast_biophysics
   use meds_stepper,             only : advance_one_step
   use meds_test_support,        only : build_test_config, check, check_close, banner
   use meds_time,                only : meds_time_t
   use meds_forcing_config,      only : MET_BACKEND_NETCDF, SWPART_CLEARIDX, METAVG_END
   use meds_forcing_types,       only : met_driver_t
   use meds_met_driver,          only : met_open, met_close
   use meds_netcdf_c
   use iso_c_binding,            only : c_int, c_size_t, c_double
   implicit none

   integer(ik), parameter :: nsl = 10_ik
   type(meds_config_t)  :: cfg
   type(site_t)         :: site
   type(fast_context_t) :: ctx
   real(wp)    :: we, ww, t_cas0, t_cas1, theta0_1, theta1_1, psi0_leaf, psi1_leaf
   real(wp)    :: cbal0, cbal1
   integer(ik) :: nfail

   call banner('site-level fast-biophysics loop (state-hub reservoirs)')

   cfg = build_test_config()
   cfg%fast_biophysics_on = .true.
   cfg%dt_fast            = 900.0_wp
   cfg%n_fast_per_slow    = 8_ik            ! 8 x 900 s = 2 h of fast integration per slow step

   !----- Build the (MVP) column config + reference met inside the fast context. -----------!
   call build_soil_params(nsl, SOIL_RETENTION_VG, 2.0_wp, 3.0_wp, 0.43_wp, 0.078_wp,            &
                          2.89e-6_wp, 3.6_wp, 1.56_wp, 2.0_wp, -3.37_wp, ctx%ccfg%soil)
   call build_soil_thermal(nsl, 3.0_wp, 0.15_wp, 2.0e6_wp, ctx%ccfg%soil_thermal)
   ctx%ccfg%wood%is_woody = .true. ; ctx%ccfg%wood%stem_resp_factor25 = 0.06_wp ; ctx%ccfg%wood%agf_bs = 0.7_wp
   ctx%ccfg%root%root_resp_factor25 = 0.30_wp
   ctx%ccfg%co2%rh_k_base = 0.01_wp
   ctx%ccfg%fast_soil_carbon = 5.0_wp
   ctx%ccfg%hydro_p%leaf_pi0 = -1.5_wp ; ctx%ccfg%hydro_p%leaf_eps = 12.0_wp
   ctx%ccfg%hydro_p%leaf_af  = 0.30_wp ; ctx%ccfg%hydro_p%leaf_water_sat = 2.0_wp
   ctx%ccfg%hydro_p%wood_pi0 = -1.0_wp ; ctx%ccfg%hydro_p%wood_eps = 8.0_wp
   ctx%ccfg%hydro_p%wood_af  = 0.20_wp ; ctx%ccfg%hydro_p%wood_water_sat = 1.0_wp
   ctx%ccfg%hydro_p%wood_psi50 = -2.0_wp ; ctx%ccfg%hydro_p%wood_kexp = 2.0_wp
   ctx%ccfg%hydro_p%k_plant_max = 6.0e-4_wp ; ctx%ccfg%hydro_p%wood_kmax = 8.0_wp
   ctx%ccfg%hydro_p%vessel_curl = 1.5_wp
   ctx%ccfg%rhizo_cond = 5.0e-4_wp
   ctx%air_temp = 290.0_wp ; ctx%rad_sw_top = 500.0_wp ; ctx%rad_sw_ground = 75.0_wp
   ctx%theta_init = 0.30_wp ; ctx%soil_temp_init = 288.0_wp

   !----- Two patches: patch 1 vegetated, patch 2 BARE (zero cohorts). ---------------------!
   call init_bare_ground(site, cfg, 2_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)     ! patch 1: one tree cohort
   call finalize_init(site)
   call init_fast_reservoirs(site, ctx)

   t_cas0    = site%patch%cas(1)%can_temp
   theta0_1  = site%patch%soil_w(1)%theta(1)
   psi0_leaf = site%cohort%psi(1, 1)                 ! patch-1 cohort-1 leaf-node psi (== PSI_INIT)

   !=== 1+2. Run the fast loop directly; conservation + activity + per-cohort persistence. =!
   call run_fast_biophysics(site, ctx, cfg, worst_energy=we, worst_water=ww, n_budget_fail=nfail)
   t_cas1    = site%patch%cas(1)%can_temp
   theta1_1  = site%patch%soil_w(1)%theta(1)
   psi1_leaf = site%cohort%psi(1, 1)

   call check(nfail == 0_ik, 'whole-column budgets closed on every patch (n_fail == 0)')
   call check(we < 1.0e-3_wp, 'whole-column energy residual tiny')
   call check(ww < 1.0e-8_wp, 'whole-column water residual tiny')
   call check(abs(t_cas1 - t_cas0) > 0.05_wp, 'CAS temperature evolved under the fast loop')
   call check(site%patch%cas(2)%can_temp > 200.0_wp, 'bare (zero-cohort) patch fast step stayed physical')
   !----- The fast loop READ psi from the cohort block, evolved it, and WROTE it back (persist). !
   call check(psi1_leaf < psi0_leaf - 1.0e-4_wp, 'per-cohort leaf psi evolved + persisted on the cohort block')
   !----- Fast->slow carbon bridge: the vegetated cohort accumulated positive GROSS GPP. -------!
   call check(site%cohort%gpp_accum(1) > 0.0_wp, 'fast loop accumulated positive gross GPP (fast->slow bridge)')

   !=== 3. Stepper wiring + gate. The gate is cfg%fast_biophysics_on; when it is ON a fast      !
   !    context is REQUIRED (advance_one_step error-stops otherwise -- the hardened BUG1 guard,  !
   !    not exercised here since asserting an abort would kill the process). =====================!
   block
      real(wp) :: t_before
      !----- Gate OFF (fast_biophysics_on = .false.): the fast loop must NOT run, even with a ctx. !
      cfg%fast_biophysics_on = .false.
      call init_fast_reservoirs(site, ctx)
      t_before = site%patch%cas(1)%can_temp
      call advance_one_step(site, cfg, .false., .false., ctx)     ! gate off -> fast loop skipped
      call check_close(site%patch%cas(1)%can_temp, t_before, 1.0e-12_wp, &
                       'fast loop must NOT run when fast_biophysics_on is off')
      !----- Gate ON (fast_biophysics_on = .true. + context): the hook fires, reservoirs change. -!
      cfg%fast_biophysics_on = .true.
      call init_fast_reservoirs(site, ctx)
      t_before = site%patch%cas(1)%can_temp
      call advance_one_step(site, cfg, .false., .false., ctx)
      call check(abs(site%patch%cas(1)%can_temp - t_before) > 0.05_wp, &
                 'advance_one_step ran the fast loop when gated on')
   end block

   !=== 4. END-TO-END fast->slow handoff via the stepper. The SAME carbon-mode site is advanced  !
   !    with the fast loop OFF then ON, with gpp_ref = 0 so the stub adds nothing. Turnover is    !
   !    identical between the two, so more carbon under fast-ON is EXACTLY the fast GPP flowing    !
   !    through gpp_accum into carbon growth -- the fast->slow bridge, isolated from turnover.  ===!
   cfg%growth_source = GS_CARBON
   cfg%gpp_ref       = 0.0_wp

   cfg%fast_biophysics_on = .false.                    ! stub GPP = 0; carbon change is pure turnover
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
   call finalize_init(site)
   call advance_one_step(site, cfg, .false., .false., ctx)
   cbal0 = site%cohort%leaf_carbon(1) + site%cohort%fineroot_carbon(1)                          &
         + site%cohort%wood_carbon(1) + site%cohort%nonstructural_carbon(1)

   cfg%fast_biophysics_on = .true.                     ! fast loop -> gpp_accum -> carbon growth
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
   call finalize_init(site)
   call init_fast_reservoirs(site, ctx)
   call advance_one_step(site, cfg, .false., .false., ctx)
   cbal1 = site%cohort%leaf_carbon(1) + site%cohort%fineroot_carbon(1)                          &
         + site%cohort%wood_carbon(1) + site%cohort%nonstructural_carbon(1)

   call check(cbal1 > cbal0, 'fast GPP raised carbon vs fast-off (gpp_ref=0): the fast->slow bridge is live')
   !----- BUG3: the fast loop accumulated POSITIVE maintenance respiration (leaf Rd + stem + root), !
   !      which carbon_growth now nets out of GPP. Before the fix these were never tracked. --------!
   call check(site%cohort%leaf_resp_accum(1) > 0.0_wp .and. site%cohort%stem_resp_accum(1) > 0.0_wp &
              .and. site%cohort%root_resp_accum(1) > 0.0_wp,                                        &
              'fast loop accumulated positive leaf/stem/root maintenance respiration (BUG3)')

   !=== 5. FORCING-DRIVEN diurnal cycle (Phase B end-to-end). Drive the SAME fast loop from a real  !
   !    MEDS forcing NetCDF (SW peaks ~17 UTC at Ithaca). A 2 h NIGHT window (02-04 UTC, SW=0) must   !
   !    accumulate ~0 GPP; a 2 h DAY window (15-17 UTC, high SW) must accumulate clearly more -- i.e. !
   !    the forcing's diurnal shortwave really propagates through met_instant -> the fast loop.  =====!
   block
      type(met_driver_t) :: drv
      real(wp)           :: gpp_night, gpp_day
      cfg%growth_source = GS_EMPIRICAL ; cfg%gpp_ref = 0.0_wp     ! isolate the FORCING-driven GPP
      call write_diurnal_forcing('test_fast_loop_forcing.nc')
      cfg%forcing%forcing_on   = .true.
      cfg%forcing%backend      = MET_BACKEND_NETCDF
      cfg%forcing%path         = 'test_fast_loop_forcing.nc'
      cfg%forcing%grid_index   = 1_ik ; cfg%forcing%dt_forcing = 3600.0_wp
      cfg%forcing%avg_convention = METAVG_END ; cfg%forcing%sw_partition = SWPART_CLEARIDX
      cfg%forcing%latitude_deg = 42.44_wp ; cfg%forcing%longitude_deg = -76.50_wp
      cfg%forcing%utc_offset_h = 0.0_wp ; cfg%forcing%apply_solar_longitude = .true.
      cfg%forcing%recycle = .false.

      call init_bare_ground(site, cfg, 1_ik)
      call add_cohort(site, cfg, 1_ik, 1_ik, 0.3_wp, 16.0_wp)
      call finalize_init(site)
      call met_open(drv, cfg%forcing)

      call init_fast_reservoirs(site, ctx)
      call run_fast_biophysics(site, ctx, cfg, met_drv=drv,                                       &
                               step_start=meds_time_t(2020_ik,7_ik,1_ik,2_ik))    ! 02-04 UTC (night)
      gpp_night = site%cohort%gpp_accum(1)

      call init_fast_reservoirs(site, ctx)
      call run_fast_biophysics(site, ctx, cfg, met_drv=drv,                                       &
                               step_start=meds_time_t(2020_ik,7_ik,1_ik,15_ik))   ! 15-17 UTC (day)
      gpp_day = site%cohort%gpp_accum(1)
      call met_close(drv)

      call check(gpp_night < 1.0e-9_wp, 'night forcing window -> ~zero GPP (SW=0 propagated through met_instant)')
      call check(gpp_day > 1.0e-7_wp,   'day forcing window -> positive GPP')
      call check(gpp_day > 10.0_wp * max(gpp_night, 1.0e-30_wp), 'GPP tracks the diurnal shortwave (day >> night)')
      write(*,'(a,es10.3,a,es10.3,a)') '   (forcing GPP: night=', gpp_night, '  day=', gpp_day, ' kgC/plant)'
   end block

   write(*,'(a)')          '   PASS'
   write(*,'(a,f7.2,a,f7.2,a)') '   (patch-1 CAS temp ', t_cas0, ' -> ', t_cas1, ' K)'
   write(*,'(a,es10.3,a,es10.3,a)') '   (worst whole-column resid: energy=', we, ' J/m2  water=', ww, ' kg/m2)'

contains

   !----- Write a 1-day hourly (time=25, grid=1) MEDS forcing NetCDF with a diurnal shortwave    !
   !      (peaks ~17 UTC at Ithaca) so the fast loop sees a real day/night cycle.                 !
   subroutine write_diurnal_forcing(path)
      character(len=*), intent(in) :: path
      integer, parameter :: NT = 25
      integer(c_int) :: st, ncid, td, gd, vt, vla, vlo, vv(8), dims2(2), dims1(1)
      integer(c_size_t) :: s1(1), c1(1), s2(2), c2(2)
      real(c_double) :: tsec(NT), dat(1, NT)
      integer :: it, k
      real(wp) :: hh, sw
      character(len=8), parameter :: nm(8) = ['Tair    ','Qair    ','PSurf   ','Wind    ',        &
                                              'Rainf   ','SWdown  ','LWdown  ','CO2air  ']
      st = nc_create_f(path, NC_NETCDF4, ncid) ; call nc_check(st,'c')
      st = nc_def_dim_f(ncid,'time',int(NT,c_size_t),td) ; call nc_check(st,'t')
      st = nc_def_dim_f(ncid,'grid',1_c_size_t,gd) ; call nc_check(st,'g')
      dims1(1)=td ; st = nc_def_var_f(ncid,'time',NC_DOUBLE,1,dims1,vt) ; call nc_check(st,'tv')
      st = nc_put_att_text_f(ncid,vt,'units',int(len_trim('seconds since 2020-07-01 00:00:00'),c_size_t), &
                             'seconds since 2020-07-01 00:00:00') ; call nc_check(st,'u')
      dims1(1)=gd ; st = nc_def_var_f(ncid,'latitude',NC_DOUBLE,1,dims1,vla) ; call nc_check(st,'la')
      st = nc_def_var_f(ncid,'longitude',NC_DOUBLE,1,dims1,vlo) ; call nc_check(st,'lo')
      dims2(1)=td ; dims2(2)=gd
      do k=1,8 ; st = nc_def_var_f(ncid,trim(nm(k)),NC_DOUBLE,2,dims2,vv(k)) ; call nc_check(st,'v') ; end do
      st = nc_enddef(ncid) ; call nc_check(st,'e')
      do it=1,NT ; tsec(it)=real(it-1,c_double)*3600.0_c_double ; end do
      s1=0_c_size_t ; c1(1)=int(NT,c_size_t) ; st=nc_put_vara_double(ncid,vt,s1,c1,tsec) ; call nc_check(st,'tw')
      c1(1)=1_c_size_t
      block ; real(c_double) :: la(1),lo(1) ; la=42.44_c_double ; lo=-76.50_c_double
         st=nc_put_vara_double(ncid,vla,s1,c1,la) ; call nc_check(st,'law')
         st=nc_put_vara_double(ncid,vlo,s1,c1,lo) ; call nc_check(st,'low') ; end block
      s2=0_c_size_t ; c2=[int(NT,c_size_t),1_c_size_t]
      do k=1,8
         do it=1,NT
            hh=real(it-1,wp)
            sw = max(0.0_wp, 850.0_wp*sin(3.14159265_wp*(hh-11.0_wp)/12.0_wp))
            if (hh<11.0_wp .or. hh>23.0_wp) sw=0.0_wp
            select case (k)
            case(1) ; dat(1,it)=290.0_wp + 6.0_wp*sin(2.0_wp*3.14159265_wp*(hh-15.0_wp)/24.0_wp)
            case(2) ; dat(1,it)=0.010_wp
            case(3) ; dat(1,it)=97500.0_wp
            case(4) ; dat(1,it)=2.5_wp
            case(5) ; dat(1,it)=0.0_wp
            case(6) ; dat(1,it)=sw
            case(7) ; dat(1,it)=340.0_wp
            case(8) ; dat(1,it)=415.0_wp
            end select
         end do
         st=nc_put_vara_double(ncid,vv(k),s2,c2,dat) ; call nc_check(st,'vw')
      end do
      st = nc_close(ncid) ; call nc_check(st,'cl')
   end subroutine write_diurnal_forcing

end program test_fast_loop
