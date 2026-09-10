!==========================================================================================!
! test_capi_demography -- COVERAGE FOR THE DEMOGRAPHY C-API SHIM.                             !
!                                                                                          !
! WHY THIS TEST EXISTS, in the sharpest possible terms: when it was written, the shim it covers  !
! DID NOT COMPILE against main. PR #137 changed `apply_recruitment`'s signature and the shim      !
! still passed the old argument, and 42/42 tests passed on both back ends for a day because the   !
! only thing that built this file was the optional `-DMEDS_BUILD_PYLIB=ON` library. That is the    !
! same hole as #95 -> #100, one subsystem over, and structure-plan §7.6 #4 asked for exactly this   !
! target to close it.                                                                              !
!                                                                                          !
! Two things are therefore load-bearing here and neither is an assertion:                          !
!   * the shim is COMPILED by a mandatory target, so a signature or struct change is a build break;!
!   * the test runs from the SOURCE directory against the SHIPPED example config, so that config    !
!     has a build-time consumer. It had none, and had drifted 37 required keys behind the schema    !
!     -- the example could not load at all. A config nothing loads rots exactly like a shim         !
!     nothing compiles.                                                                             !
!                                                                                          !
! The assertions cover what neither of those can: that the opaque-handle registry hands out         !
! distinct handles, releases them, and that a site actually advances behind one.                     !
!==========================================================================================!
program test_capi_demography
   use, intrinsic :: iso_c_binding, only : c_double, c_int, c_char, c_long
   use meds_kinds,            only : wp
   use meds_capi_demography,  only : meds_config_load, meds_config_n_pft, meds_config_dt_years,   &
                                     meds_site_create, meds_site_init_bare, meds_site_free,       &
                                     meds_advance_slow, meds_site_n_patch, meds_site_n_cohort,    &
                                     meds_site_total_agb, meds_site_generation
   use meds_test_support,     only : check, check_close, banner
   implicit none

   character(len=*), parameter :: CFG = 'examples/example_demography/example_config_main.toml'
   integer(c_int) :: ch, sh, sh2, npft, npatch
   real(c_double) :: dt_yr, agb0, agb1
   integer(c_long) :: gen0, gen1
   integer :: k
   character(kind=c_char), allocatable :: cpath(:)

   call banner('demography C-API shim')

   !=== 1. The SHIPPED example config loads. ================================================!
   !      This is the assertion that would have caught the 37-key schema drift. It is a real  !
   !      consumer of examples/example_demography/, which previously had none.                !
   cpath = to_c(CFG)
   ch = meds_config_load(cpath, int(len(CFG), c_int))
   call check(ch > 0_c_int, 'the shipped example config loads through the C-API')
   npft  = meds_config_n_pft(ch)
   dt_yr = meds_config_dt_years(ch)
   call check(npft >= 1_c_int, 'config reports a sane PFT count')
   call check(dt_yr > 0.0_wp .and. dt_yr <= 1.0_wp, 'config reports a sane dt_years')

   !=== 2. The opaque-handle registry hands out DISTINCT handles and releases them. ==========!
   !      This is the shim's own state -- the one thing no kernel test can cover -- and the    !
   !      reason this file is not split further (see its header).                              !
   sh  = meds_site_create()
   sh2 = meds_site_create()
   call check(sh > 0_c_int .and. sh2 > 0_c_int, 'site_create returns valid handles')
   call check(sh /= sh2, 'two live sites get DISTINCT handles')
   call meds_site_free(sh2)
   sh2 = meds_site_create()
   call check(sh2 /= sh, 'a handle freed and re-taken does not collide with the live one')
   call meds_site_free(sh2)

   !=== 3. A site behind a handle initialises and ADVANCES. =================================!
   call meds_site_init_bare(sh, ch, 2_c_int)
   npatch = meds_site_n_patch(sh)
   call check(npatch == 2_c_int, 'init_bare honours the requested patch count')
   call check(meds_site_n_cohort(sh) >= 0_c_int, 'a bare site reports a cohort count')

   gen0 = meds_site_generation(sh)
   agb0 = meds_site_total_agb(sh)
   !----- A month of daily steps: enough to cross the monthly fuse/fission cadence once, which  !
   !      is where the recruitment path this shim drives actually does something.  -------------!
   do k = 1, 31
      call meds_advance_slow(sh, ch, merge(1_c_int, 0_c_int, k == 1), merge(1_c_int, 0_c_int, k == 1))
   end do
   gen1 = meds_site_generation(sh)
   agb1 = meds_site_total_agb(sh)
   call check(gen1 > gen0, 'advance_slow bumps the site generation counter')
   call check(agb1 >= 0.0_wp, 'total_agb is readable after advancing')

   call meds_site_free(sh)
   print '(a)', 'test_capi_demography: ALL PASSED'

contains

   !----- A NUL-free c_char array; the shim takes (buffer, length), not a C string. ----------!
   function to_c(s) result(a)
      character(len=*), intent(in)        :: s
      character(kind=c_char), allocatable :: a(:)
      integer :: i
      allocate(a(len(s)))
      do i = 1, len(s)
         a(i) = s(i:i)
      end do
   end function to_c

end program test_capi_demography
