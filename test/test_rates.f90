!----- External rate application: growth (+cap), mortality survivorship, recruitment, swap. -!
program test_rates
   use meds_kinds,         only : wp, ik
   use meds_constants,     only : yr_day
   use meds_config,        only : meds_config_t, build_config
   use meds_demography_types,         only : community
   use meds_setup,         only : init_bare_ground, add_cohort, finalize_init
   use meds_demography_kernels,       only : compute_competition, growth_step, mortality_accumulate,  &
                                  mortality_apply_month
   use meds_recruitment,   only : recruitment_month
   use meds_step,          only : advance_one_step
   use meds_diagnostics,   only : mean_dbh
   use meds_test_support,  only : check, check_close, banner, const_provider_t
   implicit none

   type(meds_config_t)    :: cfg
   type(community)        :: comm
   type(const_provider_t) :: prov
   integer(ik)            :: istep, nday, n0
   real(wp)               :: dexp, mA, mB

   call banner('rate application: growth/mortality/recruitment/swap')
   cfg = build_config()
   cfg%veget_dyn_on = .false.       ! isolate the daily rate kernels

   !=== Growth: constant rate over a year advances DBH by g*t. ============================!
   call init_bare_ground(comm, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call add_cohort(comm, cfg, 1_ik, 3_ik, 0.1_wp, 10.0_wp)   ! climax, dbh_crit=120
   call finalize_init(comm)
   prov = const_provider_t(g = 2.0_wp, m = 0.0_wp, r = 0.0_wp)
   nday = 365_ik
   do istep = 1_ik, nday
      call compute_competition(comm)
      call growth_step(comm, prov, cfg%dt_years)
   end do
   dexp = 10.0_wp + 2.0_wp * real(nday, wp) / yr_day
   call check_close(comm%coh%dbh(1), dexp, 1.0e-6_wp, 'constant growth did not advance DBH')

   !=== Growth cap: DBH clamps at dbh_crit, never overshoots. =============================!
   call init_bare_ground(comm, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call add_cohort(comm, cfg, 1_ik, 3_ik, 0.1_wp, 119.0_wp)
   call finalize_init(comm)
   prov = const_provider_t(g = 100.0_wp, m = 0.0_wp, r = 0.0_wp)
   do istep = 1_ik, 30_ik
      call growth_step(comm, prov, cfg%dt_years)
   end do
   call check(comm%coh%dbh(1) <= 120.0_wp + 1.0e-12_wp, 'DBH overshot dbh_crit')
   call check_close(comm%coh%dbh(1), 120.0_wp, 1.0e-9_wp, 'DBH did not clamp at dbh_crit')

   !=== Mortality: 30 daily accumulations then monthly apply -> exp survivorship. ==========!
   call init_bare_ground(comm, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call add_cohort(comm, cfg, 1_ik, 1_ik, 1.0_wp, 10.0_wp)
   call finalize_init(comm)
   prov = const_provider_t(g = 0.0_wp, m = 0.1_wp, r = 0.0_wp)
   do istep = 1_ik, 30_ik
      call mortality_accumulate(comm, prov, cfg%dt_years)
   end do
   call mortality_apply_month(comm, cfg%negligible_nplant)
   call check_close(comm%coh%nplant(1), exp(-0.1_wp * 30.0_wp / yr_day), 1.0e-9_wp,        &
                    'survivorship /= exp(-m*dt)')

   !=== Recruitment: spawns above threshold, blocked when the temperature gate fails. ======!
   call init_bare_ground(comm, cfg, 1_ik, 298.15_wp, 295.15_wp)
   call finalize_init(comm)
   prov = const_provider_t(g = 0.0_wp, m = 0.0_wp, r = 0.02_wp, rmin = -1.0e30_wp)
   call recruitment_month(comm, cfg, prov)         ! pool 0.02 >= min_recruit_size 0.01
   call check(comm%coh%n == comm%n_pft, 'recruitment should spawn one cohort per PFT')
   n0 = comm%coh%n
   prov = const_provider_t(g = 0.0_wp, m = 0.0_wp, r = 0.02_wp, rmin = 400.0_wp)  ! gate shut
   call recruitment_month(comm, cfg, prov)
   call check(comm%coh%n == n0, 'temperature gate failed to block recruitment')

   !=== Provider swap: same engine binary, two providers -> different outcomes. ============!
   mA = grow_for(cfg, 0.5_wp)
   mB = grow_for(cfg, 2.0_wp)
   call check(mB > mA + 1.0_wp, 'provider swap did not change the result via dispatch')

   write(*,'(a)') '   PASS'

contains

   real(wp) function grow_for(cfg_in, g) result(dm)
      type(meds_config_t), intent(in) :: cfg_in
      real(wp),            intent(in) :: g
      type(community)        :: cc
      type(const_provider_t) :: pp
      integer(ik)            :: s
      call init_bare_ground(cc, cfg_in, 1_ik, 298.15_wp, 295.15_wp)
      call add_cohort(cc, cfg_in, 1_ik, 3_ik, 0.1_wp, 10.0_wp)
      call finalize_init(cc)
      pp = const_provider_t(g = g, m = 0.0_wp, r = 0.0_wp)
      do s = 1_ik, 730_ik          ! ~2 years: Delta-g of 1.5 cm/yr => ~3 cm divergence
         call advance_one_step(cc, pp, cfg_in, .false., .false.)
      end do
      dm = mean_dbh(cc)
   end function grow_for

end program test_rates
