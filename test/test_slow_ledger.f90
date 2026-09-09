!==========================================================================================!
! test_slow_ledger -- the site conservation ledger's ARITHMETIC and its STORE COVERAGE.       !
!                                                                                          !
! The ledger reports numbers nobody can check by inspection, so what it must never do is be   !
! quietly BLIND. Two ways it can be, and one test for each:                                   !
!                                                                                          !
!   1. THE STORE MISSES A RESERVOIR. If slow_site_store forgets a field, every operator that   !
!      moves that field reads as perfectly conservative -- the ledger prints zeros and the      !
!      defect it exists to find is the one thing it cannot see. Fortran cannot enumerate a      !
!      derived type's components, so this is done the way test_state_combinators does it:      !
!      perturb each reservoir in turn and assert the total MOVES, by the right amount and in    !
!      the right currency. Adding a state field means adding a row here.                        !
!                                                                                          !
!   2. THE DECLARATION ARITHMETIC IS WRONG. A declared term is subtracted from the residual, so  !
!      a sign error would let a real leak cancel against a real source. Asserted directly.       !
!                                                                                          !
! The third test binds the canopy-air open-volume terms to the store they are supposed to        !
! explain: resize the control volume, and the three declared terms must equal the three store    !
! changes EXACTLY. That is the one declaration in the skeleton that does real arithmetic rather   !
! than passing a quantity through, and the run it was written against closes it to 1e-17.         !
!==========================================================================================!
program test_slow_ledger
   use meds_kinds,              only : wp, ik
   use meds_config,             only : meds_config_t
   use meds_site_state_types,   only : site_t
   use meds_column_state_types, only : cas_set_depth
   use meds_column_params,      only : soil_params_t, build_soil_hydr_params
   use meds_init,               only : init_bare_ground, add_cohort, finalize_init
   use meds_slow_ledger,        only : slow_store_t, slow_ledger_t, slow_site_store,               &
                                       slow_ledger_open, slow_ledger_declare, slow_ledger_mark,    &
                                       slow_fast_carbon_handover, SLOW_PHASE_GROW
   use meds_test_support,       only : build_test_config, check, check_close, banner
   implicit none

   real(wp), parameter :: RHO = 1.2_wp        !< the density the canopy-air store is valued at
   type(meds_config_t) :: cfg
   type(site_t)        :: site
   type(soil_params_t) :: soil
   type(slow_store_t)  :: s0, s1
   type(slow_ledger_t) :: led
   real(wp)            :: de, dw, dc, expect, hand

   call banner('slow ledger: store coverage and declaration arithmetic')

   cfg = build_test_config()
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 2_ik, 0.3_wp, 18.0_wp)
   call add_cohort(site, cfg, 1_ik, 1_ik, 0.8_wp,  6.0_wp)
   call finalize_init(site)
   call check(site%cohort%n == 2_ik .and. site%patch%n == 1_ik, 'fixture: 2 cohorts in 1 patch')

   associate (sc => cfg%soil_column)
      call build_soil_hydr_params(sc%n_layer, sc%retention, sc%depth, sc%grid_growth, sc%theta_sat, &
                                  sc%theta_res, sc%ksat, sc%curve_par_a, sc%curve_par_n,            &
                                  sc%root_beta, sc%psi_fc, soil)
   end associate

   !----- Give every reservoir a non-zero value, so a perturbation below is a change to a store !
   !      that already exists rather than the creation of one. --------------------------------!
   site%patch%soil_w(1)%theta(1:soil%n_active)      = 0.30_wp
   site%patch%soil_e(1)%soil_energy(1:soil%n_active) = 5.0e7_wp
   site%patch%soil_w(1)%w_surface      = 2.0_wp
   site%patch%soil_w(1)%w_surface_enth = 8.0e5_wp
   site%patch%snow(1)%swe(1)           = 5.0_wp
   site%patch%snow(1)%snow_energy(1)   = -1.0e6_wp
   site%patch%cas(1)%can_shv      = 0.008_wp
   site%patch%cas(1)%can_enthalpy = 2.9e5_wp
   site%patch%cas(1)%can_co2      = 400.0_wp
   site%patch%cas(1)%can_depth    = 20.0_wp
   site%patch%soil_carbon(1)%slow_carbon = 3.0_wp
   site%patch%recruit_pool(:,1)          = 0.0_wp
   site%cohort%leaf_water_mass(1:2) = 0.02_wp
   site%cohort%wood_water_mass(1:2) = 0.05_wp
   site%cohort%leaf_surf_water(1:2) = 0.01_wp
   site%cohort%wood_surf_water(1:2) = 0.03_wp

   !=== 1. STORE COVERAGE: every reservoir must reach the total, in its own currency. =========!
   call moves_water ('soil theta',        perturb_theta,        0.01_wp * soil%dz(1) * 1000.0_wp)
   call moves_water ('pond w_surface',    perturb_pond,         1.0_wp)
   call moves_water ('snow swe',          perturb_swe,          1.0_wp)
   call moves_water ('canopy vapour',     perturb_shv,          RHO * 20.0_wp * 0.001_wp)
   call moves_water ('tissue water',      perturb_tissue_water, site%cohort%nplant(1) * 0.01_wp)
   call moves_water ('interception film', perturb_film,         0.01_wp)
   call moves_energy('soil energy',       perturb_soil_e,       1.0e6_wp * soil%dz(1))
   call moves_energy('pond enthalpy',     perturb_pond_enth,    1.0e5_wp)
   call moves_energy('snow enthalpy',     perturb_snow_enth,    1.0e5_wp)
   call moves_energy('canopy enthalpy',   perturb_cas_enth,     RHO * 20.0_wp * 1.0e4_wp)
   call moves_carbon('live pools',        perturb_live_c,       site%cohort%nplant(1) * 0.5_wp)
   call moves_carbon('CENTURY pools',     perturb_necromass,    1.0_wp)
   !----- Tissue HEAT and canopy CO2 have no closed-form expectation worth restating here (they  !
   !      are the store's own formulae); assert only that they reach the total, which is the     !
   !      blindness this test is about. ----------------------------------------------------!
   call moves('tissue heat',   perturb_tissue_temp, 'energy')
   call moves('canopy CO2',    perturb_cas_co2,     'carbon')
   call moves('recruit pool',  perturb_recruit,     'carbon')

   !=== 2. DECLARATION ARITHMETIC: a declared term is subtracted from the residual. ===========!
   led%active = .true.
   call slow_ledger_open(led, site, cfg, RHO)
   call slow_ledger_declare(led, carbon_in = 1.0_wp, water_out = 2.0_wp, energy_in = 3.0_wp)
   call slow_ledger_mark(led, site, cfg, SLOW_PHASE_GROW)   ! nothing moved: residual = -(in - out)
   call check_close(led%resid_sum(SLOW_PHASE_GROW)%carbon, -1.0_wp, 1.0e-12_wp, 'declared carbon_in is subtracted')
   call check_close(led%resid_sum(SLOW_PHASE_GROW)%water,   2.0_wp, 1.0e-12_wp, 'declared water_out is added back')
   call check_close(led%resid_sum(SLOW_PHASE_GROW)%energy, -3.0_wp, 1.0e-12_wp, 'declared energy_in is subtracted')

   !----- An UNDECLARED, UNCHANGED step must be exactly zero, or every phase carries a bias. ----!
   led = slow_ledger_t()
   led%active = .true.
   call slow_ledger_open(led, site, cfg, RHO)
   call slow_ledger_mark(led, site, cfg, SLOW_PHASE_GROW)
   call check(led%resid_sum(SLOW_PHASE_GROW)%carbon == 0.0_wp .and.                                &
              led%resid_sum(SLOW_PHASE_GROW)%water  == 0.0_wp .and.                                &
              led%resid_sum(SLOW_PHASE_GROW)%energy == 0.0_wp, 'a step that changes nothing closes at exactly zero')

   !=== 3. The canopy-air OPEN-VOLUME terms must equal the store change they explain. ==========!
   s0 = slow_site_store(site, cfg, soil, RHO)
   call cas_set_depth(site%patch%cas(1), 8.0_wp, rho_air=RHO, de_open=de, dw_open=dw, dc_open=dc)
   s1 = slow_site_store(site, cfg, soil, RHO)
   call check_close(s1%energy - s0%energy, de,             1.0e-6_wp, 'de_open explains the canopy energy change')
   call check_close(s1%water  - s0%water,  dw,             1.0e-12_wp, 'dw_open explains the canopy vapour change')
   call check_close(s1%carbon - s0%carbon, dc * 1.2e-8_wp, 1.0e-14_wp, 'dc_open explains the canopy CO2 change')
   call check(de < 0.0_wp, 'shrinking the canopy DETRAINS (negative)')

   !=== 4. The fast->slow handover is nplant- and area-weighted. ==============================!
   cfg%fast_biophysics_on = .true.
   site%cohort%gpp_accum(1:2)       = 1.0e-3_wp
   site%cohort%leaf_resp_accum(1:2) = 2.0e-4_wp
   site%cohort%stem_resp_accum(1:2) = 1.0e-4_wp
   site%cohort%root_resp_accum(1:2) = 1.0e-4_wp
   expect = site%patch%area(1) * (site%cohort%nplant(1) + site%cohort%nplant(2)) * 6.0e-4_wp
   hand   = slow_fast_carbon_handover(site, cfg)
   call check_close(hand, expect, 1.0e-14_wp, 'handover = area * sum(nplant * (gpp - maintenance resp))')

   print '(a)', 'test_slow_ledger: ALL PASSED'

contains

   !----- Assert a perturbation moves the named currency by `amount` (and moves it at all). ----!
   subroutine moves_water(name, perturb, amount)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: amount
      interface
         subroutine perturb()
         end subroutine perturb
      end interface
      type(slow_store_t) :: a, b
      a = slow_site_store(site, cfg, soil, RHO)
      call perturb()
      b = slow_site_store(site, cfg, soil, RHO)
      call check_close(b%water - a%water, amount, 1.0e-9_wp * max(1.0_wp, abs(amount)),            &
                       'store: '//name//' reaches the water total')
   end subroutine moves_water

   subroutine moves_energy(name, perturb, amount)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: amount
      interface
         subroutine perturb()
         end subroutine perturb
      end interface
      type(slow_store_t) :: a, b
      a = slow_site_store(site, cfg, soil, RHO)
      call perturb()
      b = slow_site_store(site, cfg, soil, RHO)
      call check_close(b%energy - a%energy, amount, 1.0e-6_wp * max(1.0_wp, abs(amount)),          &
                       'store: '//name//' reaches the energy total')
   end subroutine moves_energy

   subroutine moves_carbon(name, perturb, amount)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: amount
      interface
         subroutine perturb()
         end subroutine perturb
      end interface
      type(slow_store_t) :: a, b
      a = slow_site_store(site, cfg, soil, RHO)
      call perturb()
      b = slow_site_store(site, cfg, soil, RHO)
      call check_close(b%carbon - a%carbon, amount, 1.0e-9_wp * max(1.0_wp, abs(amount)),          &
                       'store: '//name//' reaches the carbon total')
   end subroutine moves_carbon

   subroutine moves(name, perturb, currency)
      character(len=*), intent(in) :: name, currency
      interface
         subroutine perturb()
         end subroutine perturb
      end interface
      type(slow_store_t) :: a, b
      real(wp) :: d
      a = slow_site_store(site, cfg, soil, RHO)
      call perturb()
      b = slow_site_store(site, cfg, soil, RHO)
      select case (currency)
      case ('carbon') ; d = b%carbon - a%carbon
      case ('water')  ; d = b%water  - a%water
      case default    ; d = b%energy - a%energy
      end select
      call check(abs(d) > 0.0_wp, 'store: '//name//' reaches the '//currency//' total')
   end subroutine moves

   subroutine perturb_theta()        ; site%patch%soil_w(1)%theta(1) = site%patch%soil_w(1)%theta(1) + 0.01_wp
   end subroutine perturb_theta
   subroutine perturb_pond()         ; site%patch%soil_w(1)%w_surface = site%patch%soil_w(1)%w_surface + 1.0_wp
   end subroutine perturb_pond
   subroutine perturb_swe()          ; site%patch%snow(1)%swe(1) = site%patch%snow(1)%swe(1) + 1.0_wp
   end subroutine perturb_swe
   subroutine perturb_shv()          ; site%patch%cas(1)%can_shv = site%patch%cas(1)%can_shv + 0.001_wp
   end subroutine perturb_shv
   subroutine perturb_tissue_water() ; site%cohort%leaf_water_mass(1) = site%cohort%leaf_water_mass(1) + 0.01_wp
   end subroutine perturb_tissue_water
   subroutine perturb_film()         ; site%cohort%leaf_surf_water(1) = site%cohort%leaf_surf_water(1) + 0.01_wp
   end subroutine perturb_film
   subroutine perturb_soil_e()       ; site%patch%soil_e(1)%soil_energy(1) = site%patch%soil_e(1)%soil_energy(1) + 1.0e6_wp
   end subroutine perturb_soil_e
   subroutine perturb_pond_enth()    ; site%patch%soil_w(1)%w_surface_enth = site%patch%soil_w(1)%w_surface_enth + 1.0e5_wp
   end subroutine perturb_pond_enth
   subroutine perturb_snow_enth()    ; site%patch%snow(1)%snow_energy(1) = site%patch%snow(1)%snow_energy(1) + 1.0e5_wp
   end subroutine perturb_snow_enth
   subroutine perturb_cas_enth()     ; site%patch%cas(1)%can_enthalpy = site%patch%cas(1)%can_enthalpy + 1.0e4_wp
   end subroutine perturb_cas_enth
   subroutine perturb_cas_co2()      ; site%patch%cas(1)%can_co2 = site%patch%cas(1)%can_co2 + 50.0_wp
   end subroutine perturb_cas_co2
   subroutine perturb_live_c()       ; site%cohort%wood_carbon(1) = site%cohort%wood_carbon(1) + 0.5_wp
   end subroutine perturb_live_c
   subroutine perturb_necromass()    ; site%patch%soil_carbon(1)%slow_carbon = site%patch%soil_carbon(1)%slow_carbon + 1.0_wp
   end subroutine perturb_necromass
   subroutine perturb_tissue_temp()  ; site%cohort%leaf_temp(1) = site%cohort%leaf_temp(1) + 5.0_wp
   end subroutine perturb_tissue_temp
   subroutine perturb_recruit()      ; site%patch%recruit_pool(1,1) = site%patch%recruit_pool(1,1) + 0.5_wp
   end subroutine perturb_recruit

end program test_slow_ledger
