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
   logical             :: found
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
   end do

   now = meds_time_t(year=2000, month=1, day=1, hour=0, minute=0, second=0)
   call io_write_state(site, cfg, '.', 'test_sr', now)
   call io_read_state(site2, cfg, './test_sr-S-20000101000000.nc', rt, found)

   call check(found, 'state file read back')
   call check(site2%cohort%n == site%cohort%n, 'cohort count preserved')
   do i = 1_ik, site2%cohort%n
      call check_close(site2%cohort%sla(i),     sla_set, 1.0e-9_wp, 'sla recovered from state')
      call check_close(site2%cohort%vcmax25(i), vc_set,  1.0e-9_wp, 'vcmax25 recovered from state')
      call check_close(site2%cohort%rd25(i),    rd_set,  1.0e-9_wp, 'rd25 recovered from state')
      call check_close(site2%cohort%llspan(i),  ll_set,  1.0e-9_wp, 'llspan recovered from state')
   end do

   write(*,'(a)') '   PASS'
end program test_state_roundtrip
