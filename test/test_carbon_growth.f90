!----- Carbon-driven growth: growth-respiration fraction, PARTEH carbon closure (incl. --------!
!----- reproduction), and a carbon-mode step that grows wood_carbon -> dbh (the geometry flip). !
program test_carbon_growth
   use meds_kinds,                  only : wp, ik
   use meds_config,                 only : meds_config_t, GS_CARBON
   use meds_demography_interface,   only : site_t
   use meds_plant_interface,        only : get_plant_flux_slow, growth_respiration,             &
                                           carbon_env_t, carbon_demand_t, carbon_npp_t, PHEN_ON
   use meds_init,                   only : init_bare_ground, add_cohort, finalize_init
   use meds_stepper,                only : advance_one_step
   use meds_demography_diagnostics, only : has_nan
   use meds_test_support,           only : build_test_config, check, check_close, banner
   implicit none

   type(meds_config_t)   :: cfg
   type(site_t)          :: site
   type(carbon_env_t)    :: env
   type(carbon_demand_t) :: demand
   type(carbon_npp_t)    :: npp
   real(wp)    :: a_gpp, net, sum_npp, wood0, dbh0
   integer(ik) :: istep

   call banner('carbon-driven growth')
   cfg = build_test_config()
   cfg%pft%growth_resp_factor = [ 0.66_wp, 0.66_wp, 0.66_wp ]   ! only a third of GPP becomes NPP

   !=== 1. growth respiration: with factor 0.66, NPP = GPP - 0.66*GPP = 0.34*GPP. =========!
   a_gpp = 1.0_wp
   net   = a_gpp - growth_respiration(a_gpp, 0.66_wp)
   call check_close(net, 0.34_wp, 1.0e-12_wp, 'NPP should be 0.34*GPP at growth_resp_factor=0.66')

   !=== 2. get_plant_flux_slow CLOSES carbon. With no standing tissue there is no turnover, ===!
   !    so the allocation sums exactly to net_carbon; a mature cohort banks a reproduction share.
   env%net_carbon = net ; env%nonstructural = 0.0_wp
   env%leaf_carbon = 0.0_wp ; env%fineroot_carbon = 0.0_wp        ! -> zero turnover loss
   env%tissue_temp = 298.15_wp ; env%dt_yr = 1.0_wp ; env%phenology_status = PHEN_ON
   demand%leaf = 0.05_wp ; demand%fineroot = 0.05_wp ; demand%storage = 0.02_wp
   demand%wood = 1.0e6_wp ; demand%reproduction_fraction = 0.3_wp
   call get_plant_flux_slow(env, cfg, 3_ik, demand, npp)
   sum_npp = npp%leaf + npp%fineroot + npp%wood + npp%nonstructural + npp%repro
   call check_close(sum_npp, net, 1.0e-12_wp, 'carbon not conserved: sum(npp) /= net_carbon')
   call check(.not. npp%starving,       'positive NPP must not starve')
   call check(npp%repro > 0.0_wp,       'mature cohort should allocate carbon to reproduction')
   call check(npp%wood  > 0.0_wp,       'surplus carbon should grow wood (residual sink)')

   !=== 3. NEGATIVE net carbon: no growth, storage drains, and it starves once storage is 0. ==!
   env%net_carbon = -0.1_wp ; env%nonstructural = 0.03_wp
   call get_plant_flux_slow(env, cfg, 3_ik, demand, npp)
   call check(npp%wood <= 0.0_wp, 'no wood growth on negative NPP')
   call check(npp%starving,       'storage below the debt must set starving')
   call check_close(npp%deficit, 0.07_wp, 1.0e-12_wp, 'deficit should be debt - storage')

   !=== 4. A carbon-mode step grows wood_carbon -> dbh, and leaf_area stays leaf_carbon*sla. ==!
   cfg%growth_source = GS_CARBON ; cfg%gpp_ref = 0.5_wp
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 3_ik, 0.1_wp, 20.0_wp)       ! climax cohort
   call finalize_init(site)
   wood0 = site%cohort%wood_carbon(1) ; dbh0 = site%cohort%dbh(1)
   do istep = 1_ik, 30_ik
      call advance_one_step(site, cfg, .false., .false.)         ! growth + mortality only
   end do
   call check(site%cohort%wood_carbon(1) > wood0, 'carbon growth did not add wood_carbon')
   call check(site%cohort%dbh(1) > dbh0,          'dbh did not increase with wood_carbon (the flip)')
   call check_close(site%cohort%leaf_area(1), site%cohort%leaf_carbon(1) * cfg%pft%sla(3_ik),   &
                    1.0e-9_wp, 'carbon-mode leaf_area /= leaf_carbon*sla')
   call check(.not. has_nan(site), 'carbon-mode step produced NaN')

   write(*,'(a)') '   PASS'
end program test_carbon_growth
