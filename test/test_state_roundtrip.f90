!----- State restart round-trip: the per-cohort PLASTIC leaf traits (sla/vcmax25/rd25/llspan) ----!
!----- must be written to and recovered from a state checkpoint, not reset to top-of-canopy.  ----!
program test_state_roundtrip
   use meds_kinds,            only : wp, ik
   use meds_config,           only : meds_config_t
   use meds_core_interface,   only : site_t
   use meds_init,             only : init_bare_ground, add_cohort, finalize_init
   use meds_io,               only : io_write_state, io_read_state
   use meds_time,             only : meds_time_t
   use meds_test_support,     only : build_test_config, check, check_close, banner
   implicit none

   type(meds_config_t) :: cfg
   type(site_t)        :: site, site2
   type(meds_time_t)   :: now, rt
   logical             :: found, fast_ok
   real(wp)            :: sla_set, vc_set, rd_set, ll_set
   integer(ik)         :: i

   call banner('state roundtrip: plastic traits')
   cfg = build_test_config()
   call init_bare_ground(site, cfg, 1_ik)
   call add_cohort(site, cfg, 1_ik, 3_ik, 0.2_wp, 20.0_wp)
   call finalize_init(site)
   call check(site%cohort%n >= 1_ik, 'have at least one cohort to checkpoint')

   !----- Stamp distinct, NON-top-of-canopy plastic traits on every cohort. ----------------!
   sla_set = 99.0_wp ; vc_set = 42.0_wp ; rd_set = 3.5_wp ; ll_set = 7.0_wp
   do i = 1_ik, site%cohort%n
      site%cohort%sla(i)     = sla_set ; site%cohort%vcmax25(i) = vc_set
      site%cohort%rd25(i)    = rd_set  ; site%cohort%llspan(i)  = ll_set
      !----- P6: per-cohort hydraulics/temperature, stamped to distinct NON-sentinel values so a       !
      !      round-trip that fails to persist them (leaving the lazy-init 0 / LEAF_TEMP_INIT) is caught. !
      site%cohort%leaf_water_mass(i) = 0.0123_wp + 0.001_wp * real(i, wp)
      site%cohort%wood_water_mass(i) = 0.0456_wp + 0.001_wp * real(i, wp)
      site%cohort%leaf_temp(i)       = 291.5_wp + real(i, wp)
      site%cohort%wood_temp(i)       = 289.5_wp + real(i, wp)
   end do

   !----- P5 (MEDS_ED2_RK45_DESIGN.md): stamp distinct, non-default FAST reservoir values on every  !
   !      patch -- a bare-ground site never calls init_fast_reservoirs (that is meds_main's own,      !
   !      fast_biophysics_on-gated step), so these still sit at patch_alloc's raw defaults (theta=0,   !
   !      soil_energy=0 -> 0 K) unless stamped here, exactly the discontinuity the fix closes. --------!
   call check(site%patch%n >= 1_ik, 'have at least one patch to checkpoint')
   do i = 1_ik, site%patch%n
      site%patch%cas(i)%can_enthalpy      = 3.1e5_wp ; site%patch%cas(i)%can_shv = 0.011_wp
      site%patch%cas(i)%can_co2           = 415.0_wp ; site%patch%cas(i)%can_temp = 291.0_wp
      site%patch%soil_w(i)%theta(1:3)     = [0.28_wp, 0.31_wp, 0.34_wp]
      site%patch%soil_e(i)%soil_energy(1:3) = [2.1e8_wp, 2.2e8_wp, 2.3e8_wp]
      site%patch%soil_e(i)%soil_temp(1:3)   = [284.0_wp, 285.0_wp, 286.0_wp]
      site%patch%soil_e(i)%soil_fliq(1:3)   = [0.9_wp, 0.95_wp, 1.0_wp]
      site%patch%snow(i)%swe(1)           = 4.0_wp   ; site%patch%snow(i)%snow_energy(1) = -1.0e5_wp
      site%patch%snow(i)%snow_depth(1)    = 0.02_wp  ; site%patch%snow(i)%snow_temp(1)   = 270.0_wp
      site%patch%snow(i)%snow_fliq(1)     = 0.0_wp   ; site%patch%snow(i)%nlayer         = 1_ik
   end do

   now = meds_time_t(year=2000, month=1, day=1, hour=0, minute=0, second=0)
   call io_write_state(site, cfg, '.', 'test_sr', now)
   call io_read_state(site2, cfg, './test_sr-S-20000101000000.nc', rt, found, fast_found=fast_ok)

   call check(found, 'state file read back')
   call check(fast_ok, 'fast reservoirs found in state file just written')
   call check(site2%cohort%n == site%cohort%n, 'cohort count preserved')
   do i = 1_ik, site2%cohort%n
      call check_close(site2%cohort%sla(i),     sla_set, 1.0e-9_wp, 'sla recovered from state')
      call check_close(site2%cohort%vcmax25(i), vc_set,  1.0e-9_wp, 'vcmax25 recovered from state')
      call check_close(site2%cohort%rd25(i),    rd_set,  1.0e-9_wp, 'rd25 recovered from state')
      call check_close(site2%cohort%llspan(i),  ll_set,  1.0e-9_wp, 'llspan recovered from state')
      call check_close(site2%cohort%leaf_water_mass(i), site%cohort%leaf_water_mass(i), 1.0e-12_wp, &
                        'leaf_water_mass recovered from state (not re-seeded)')
      call check_close(site2%cohort%wood_water_mass(i), site%cohort%wood_water_mass(i), 1.0e-12_wp, &
                        'wood_water_mass recovered from state (not re-seeded)')
      call check_close(site2%cohort%leaf_temp(i), site%cohort%leaf_temp(i), 1.0e-9_wp, &
                        'leaf_temp recovered from state (not reset to LEAF_TEMP_INIT)')
      call check_close(site2%cohort%wood_temp(i), site%cohort%wood_temp(i), 1.0e-9_wp, &
                        'wood_temp recovered from state')
   end do
   do i = 1_ik, site2%patch%n
      call check_close(site2%patch%cas(i)%can_enthalpy, site%patch%cas(i)%can_enthalpy, 1.0e-6_wp, &
                        'cas can_enthalpy recovered from state')
      call check_close(site2%patch%cas(i)%can_shv,      site%patch%cas(i)%can_shv,      1.0e-9_wp, &
                        'cas can_shv recovered from state')
      call check_close(site2%patch%cas(i)%can_temp,     site%patch%cas(i)%can_temp,     1.0e-9_wp, &
                        'cas can_temp recovered from state (not re-seeded to a generic value)')
      call check_close(site2%patch%soil_w(i)%theta(2),  site%patch%soil_w(i)%theta(2),  1.0e-9_wp, &
                        'soil theta recovered from state (not reset to 0 = bone dry)')
      call check_close(site2%patch%soil_e(i)%soil_temp(2), site%patch%soil_e(i)%soil_temp(2), 1.0e-9_wp, &
                        'soil temp recovered from state (not reset to 0 K)')
      call check_close(site2%patch%snow(i)%swe(1),      site%patch%snow(i)%swe(1),      1.0e-9_wp, &
                        'snow swe recovered from state')
      call check(site2%patch%snow(i)%nlayer == site%patch%snow(i)%nlayer, 'snow nlayer recovered from state')
   end do

   write(*,'(a)') '   PASS'
end program test_state_roundtrip
