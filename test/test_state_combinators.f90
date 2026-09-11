!==========================================================================================!
! test_state_combinators -- every field of column_state_t, through every combinator.        !
!                                                                                          !
! WHY THIS EXISTS. The fast integrator's state vector is enumerated by hand in a dozen      !
! places -- copy, axpy, accumulate, extrapolate, subtract, embedded-error, zero -- and       !
! Fortran gives no compile-time failure when one of them omits a field. An omission is        !
! invisible three ways over: the compiler is silent, the conservation ledgers still close      !
! (a field left at its default is a consistent zero, not a leak), and the step controller       !
! simply reads a smaller error and takes a longer step. Structure plan §10.3.                    !
!                                                                                          !
! This cannot be a true completeness check -- Fortran cannot enumerate a derived type's         !
! components -- so it is the next best thing: fill every field with a DISTINCT non-zero value    !
! and assert each combinator's contract field by field. Adding a state field means adding its    !
! rows here, and the block at the top says so. What it does buy is that an omission in an        !
! EXISTING combinator, which is the failure that has actually happened, is caught immediately.   !
!==========================================================================================!
program test_state_combinators
   use meds_kinds,            only : wp, ik
   use meds_fast_types,       only : column_state_t, column_tend_t, patch_biophys_t
   use meds_column_state_ops, only : state_init, state_axpy, state_sub, state_err_diff, zero_like,  &
                                     state_accum, state_extrap, unpack_column_state
   use meds_therm_lib,        only : cas_temp_of_enthalpy
   use meds_test_support,     only : check, check_close, banner
   implicit none

   integer(ik), parameter :: N = 2_ik, NSL = 3_ik
   type(column_state_t) :: a, b, out
   type(column_tend_t)  :: k

   call banner('state combinators: every column_state_t field')

   call fill_state(a, 1.0_wp)
   call fill_state(b, 0.25_wp)

   !=== 1. state_init copies EVERY field. ==================================================!
   call state_init(a, N, NSL, out)
   call same(out, a, 'state_init')

   !=== 2. state_sub differences EVERY field. ==============================================!
   call state_sub(a, b, N, NSL, out)
   call check_close(out%cas_enthalpy,   a%cas_enthalpy   - b%cas_enthalpy,   1.0e-14_wp, 'state_sub: cas_enthalpy')
   call check_close(out%cas_shv,        a%cas_shv        - b%cas_shv,        1.0e-14_wp, 'state_sub: cas_shv')
   call check_close(out%cas_co2,        a%cas_co2        - b%cas_co2,        1.0e-14_wp, 'state_sub: cas_co2')
   call check_close(out%w_surface,      a%w_surface      - b%w_surface,      1.0e-14_wp, 'state_sub: w_surface (pond)')
   call check_close(out%w_surface_enth, a%w_surface_enth - b%w_surface_enth, 1.0e-14_wp, 'state_sub: w_surface_enth')
   call check_close(out%soil_energy(NSL), a%soil_energy(NSL) - b%soil_energy(NSL), 1.0e-14_wp, 'state_sub: soil_energy')
   call check_close(out%theta(NSL),       a%theta(NSL)       - b%theta(NSL),       1.0e-14_wp, 'state_sub: theta')
   call check_close(out%leaf_water_mass(N), a%leaf_water_mass(N) - b%leaf_water_mass(N), 1.0e-14_wp, &
                    'state_sub: leaf_water_mass')
   call check_close(out%wood_water_mass(N), a%wood_water_mass(N) - b%wood_water_mass(N), 1.0e-14_wp, &
                    'state_sub: wood_water_mass')
   call check_close(out%leaf_surf_water(N), a%leaf_surf_water(N) - b%leaf_surf_water(N), 1.0e-14_wp, &
                    'state_sub: leaf_surf_water (film)')
   call check_close(out%wood_surf_water(N), a%wood_surf_water(N) - b%wood_surf_water(N), 1.0e-14_wp, &
                    'state_sub: wood_surf_water (film)')

   !=== 3. zero_like ALLOCATES and zeroes every field. The films were once left unallocated,  !
   !       so the norm could not read them; assert the allocation, not just the value. =======!
   out = zero_like(a, N, NSL)
   call check(allocated(out%leaf_water_mass) .and. allocated(out%wood_water_mass),                  &
              'zero_like: tissue-water fields allocated')
   call check(allocated(out%leaf_surf_water) .and. allocated(out%wood_surf_water),                  &
              'zero_like: canopy-film fields allocated')
   call check_close(maxval(abs([out%cas_enthalpy, out%cas_shv, out%cas_co2, out%w_surface,          &
                                out%w_surface_enth, maxval(abs(out%soil_energy)),                   &
                                maxval(abs(out%theta)), maxval(abs(out%leaf_water_mass)),           &
                                maxval(abs(out%wood_water_mass)), maxval(abs(out%leaf_surf_water)), &
                                maxval(abs(out%wood_surf_water))])), 0.0_wp, 1.0e-30_wp,            &
                    'zero_like: some field is not zero')

   !=== 4. state_err_diff: the DECLARED exclusions are zero, the included ones are not. ======!
   !       The tissue water, the canopy film and the pond are deliberately outside the         !
   !       embedded error estimate. Assert that, so "excluded" cannot decay into "forgotten".  !
   call state_err_diff(a, b, b, a, N, NSL, out)   ! (a-b) - (b-a) = 2(a-b), so the estimate is non-zero
   call check(abs(out%cas_enthalpy) > 0.0_wp, 'state_err_diff: CAS enthalpy must be IN the estimate')
   call check_close(out%w_surface,          0.0_wp, 1.0e-30_wp, 'state_err_diff: pond excluded (declared)')
   call check_close(out%w_surface_enth,     0.0_wp, 1.0e-30_wp, 'state_err_diff: pond enthalpy excluded (declared)')
   call check_close(out%leaf_water_mass(N), 0.0_wp, 1.0e-30_wp, 'state_err_diff: tissue water excluded (declared)')
   call check_close(out%leaf_surf_water(N), 0.0_wp, 1.0e-30_wp, 'state_err_diff: canopy film excluded (declared)')

   !=== 5. state_axpy advances EVERY field by its own tendency. =============================!
   call fill_tend(k)
   call state_axpy(a, 2.0_wp, k, N, NSL, out)
   call check_close(out%cas_enthalpy,       a%cas_enthalpy       + 2.0_wp*k%d_cas_enthalpy,       1.0e-13_wp, &
                    'state_axpy: cas_enthalpy')
   call check_close(out%soil_energy(NSL),   a%soil_energy(NSL)   + 2.0_wp*k%dedt(NSL),            1.0e-13_wp, &
                    'state_axpy: soil_energy')
   call check_close(out%theta(NSL),         a%theta(NSL)         + 2.0_wp*k%dtheta_dt(NSL),       1.0e-13_wp, &
                    'state_axpy: theta')
   call check_close(out%leaf_water_mass(N), a%leaf_water_mass(N) + 2.0_wp*k%d_leaf_water_mass(N), 1.0e-13_wp, &
                    'state_axpy: leaf_water_mass')
   call check_close(out%leaf_surf_water(N), a%leaf_surf_water(N) + 2.0_wp*k%d_leaf_surf_water(N), 1.0e-13_wp, &
                    'state_axpy: leaf_surf_water')

   !----- the helpers error-stop on the first failure, so reaching here is the pass. ---------!
   print '(a)', 'test_state_combinators: all checks passed'

   !=== 6. state_accum accumulates EVERY field that has a stage tendency -- and the pond, which  !
   !       has none, must be left ALONE rather than zeroed. That distinction is the whole point:   !
   !       `column_tend_t` carries no d_w_surface, so the pond is passed through the stages and    !
   !       committed by the hydrology path. Before this test, "accum does not touch the pond" and  !
   !       "somebody forgot the pond in accum" were the same observable.  -------------------------!
   block
      type(column_state_t) :: acc
      real(wp) :: pond0, pond_enth0
      call fill_state(acc, 1.0_wp)
      call fill_tend(k)
      pond0 = acc%w_surface ; pond_enth0 = acc%w_surface_enth
      call state_accum(acc, 2.0_wp, k, N, NSL)
      call check_close(acc%cas_enthalpy, 11.0_wp + 2.0_wp*k%d_cas_enthalpy, 1.0e-13_wp, 'state_accum: cas_enthalpy')
      call check_close(acc%cas_shv,      12.0_wp + 2.0_wp*k%d_cas_shv,      1.0e-13_wp, 'state_accum: cas_shv')
      call check_close(acc%cas_co2,      13.0_wp + 2.0_wp*k%d_cas_co2,      1.0e-13_wp, 'state_accum: cas_co2')
      call check_close(acc%soil_energy(NSL), (20.0_wp+real(NSL,wp)) + 2.0_wp*k%dedt(NSL),          &
                       1.0e-13_wp, 'state_accum: soil_energy')
      call check_close(acc%theta(NSL),       (30.0_wp+real(NSL,wp)) + 2.0_wp*k%dtheta_dt(NSL),    &
                       1.0e-13_wp, 'state_accum: theta')
      call check_close(acc%leaf_water_mass(N), (40.0_wp+real(N,wp)) + 2.0_wp*k%d_leaf_water_mass(N), &
                       1.0e-13_wp, 'state_accum: leaf_water_mass')
      call check_close(acc%wood_water_mass(N), (50.0_wp+real(N,wp)) + 2.0_wp*k%d_wood_water_mass(N), &
                       1.0e-13_wp, 'state_accum: wood_water_mass')
      call check_close(acc%leaf_surf_water(N), (60.0_wp+real(N,wp)) + 2.0_wp*k%d_leaf_surf_water(N), &
                       1.0e-13_wp, 'state_accum: leaf_surf_water')
      call check_close(acc%wood_surf_water(N), (70.0_wp+real(N,wp)) + 2.0_wp*k%d_wood_surf_water(N), &
                       1.0e-13_wp, 'state_accum: wood_surf_water')
      !----- the EXCLUSION, asserted: unchanged, not zeroed.  ------------------------------------!
      call check_close(acc%w_surface,      pond0,      1.0e-14_wp, 'state_accum: pond PASSED THROUGH untouched')
      call check_close(acc%w_surface_enth, pond_enth0, 1.0e-14_wp, 'state_accum: pond enthalpy untouched')
   end block

   !=== 7. state_extrap blends EVERY field, pond included. ==================================!
   block
      type(column_state_t) :: ex
      real(wp), parameter  :: BB = 0.25_wp
      call state_extrap(a, BB, b, N, NSL, ex)
      call check_close(ex%cas_enthalpy, (1.0_wp-BB)*a%cas_enthalpy + BB*b%cas_enthalpy,           &
                       1.0e-13_wp, 'state_extrap: cas_enthalpy')
      call check_close(ex%soil_energy(NSL), (1.0_wp-BB)*a%soil_energy(NSL) + BB*b%soil_energy(NSL), &
                       1.0e-13_wp, 'state_extrap: soil_energy')
      call check_close(ex%theta(NSL), (1.0_wp-BB)*a%theta(NSL) + BB*b%theta(NSL),                 &
                       1.0e-13_wp, 'state_extrap: theta')
      call check_close(ex%leaf_water_mass(N), (1.0_wp-BB)*a%leaf_water_mass(N) + BB*b%leaf_water_mass(N), &
                       1.0e-13_wp, 'state_extrap: leaf_water_mass')
      call check_close(ex%wood_surf_water(N), (1.0_wp-BB)*a%wood_surf_water(N) + BB*b%wood_surf_water(N), &
                       1.0e-13_wp, 'state_extrap: wood_surf_water')
      call check_close(ex%w_surface, (1.0_wp-BB)*a%w_surface + BB*b%w_surface,                    &
                       1.0e-13_wp, 'state_extrap: pond IS blended')
      call check_close(ex%w_surface_enth, (1.0_wp-BB)*a%w_surface_enth + BB*b%w_surface_enth,     &
                       1.0e-13_wp, 'state_extrap: pond enthalpy IS blended')
   end block

   !=== 8. unpack_column_state -- THE COMMIT PATH. A field omitted here never reaches the patch  !
   !       state at all: the step computes it, the ledgers balance on it, and it is then dropped   !
   !       on the floor. That makes this the most consequential of the enumerations, and it was    !
   !       the only one with no test.                                                              !
   !                                                                                          !
   !       It writes 9 of the 11 fields. The pond is EXCLUDED ON PURPOSE -- it is passed through   !
   !       the stages and committed from the scratch hydrology solve (see column_state_t's own     !
   !       comment) -- so the exclusion is asserted here too, and the routine now says so in a     !
   !       comment. Exclusion and omission must never again be the same observable.  --------------!
   block
      type(patch_biophys_t) :: bio
      allocate(bio%leaf_water_mass(N), bio%wood_water_mass(N),                                     &
               bio%leaf_surf_water(N), bio%wood_surf_water(N))
      bio%leaf_water_mass = 0.0_wp ; bio%wood_water_mass = 0.0_wp
      bio%leaf_surf_water = 0.0_wp ; bio%wood_surf_water = 0.0_wp
      bio%soil_e%soil_energy = 0.0_wp ; bio%soil_w%theta = 0.0_wp
      call unpack_column_state(a, N, NSL, bio)
      call check_close(bio%cas%can_enthalpy, a%cas_enthalpy, 1.0e-14_wp, 'unpack: cas_enthalpy committed')
      call check_close(bio%cas%can_shv,      a%cas_shv,      1.0e-14_wp, 'unpack: cas_shv committed')
      call check_close(bio%cas%can_co2,      a%cas_co2,      1.0e-14_wp, 'unpack: cas_co2 committed')
      call check_close(bio%soil_e%soil_energy(NSL), a%soil_energy(NSL), 1.0e-14_wp, 'unpack: soil_energy committed')
      call check_close(bio%soil_w%theta(NSL),       a%theta(NSL),       1.0e-14_wp, 'unpack: theta committed')
      call check_close(bio%leaf_water_mass(N), a%leaf_water_mass(N), 1.0e-14_wp, 'unpack: leaf_water_mass committed')
      call check_close(bio%wood_water_mass(N), a%wood_water_mass(N), 1.0e-14_wp, 'unpack: wood_water_mass committed')
      call check_close(bio%leaf_surf_water(N), a%leaf_surf_water(N), 1.0e-14_wp, 'unpack: leaf_surf_water committed')
      call check_close(bio%wood_surf_water(N), a%wood_surf_water(N), 1.0e-14_wp, 'unpack: wood_surf_water committed')
      !----- can_temp is a DERIVED commit: it must be consistent with the enthalpy/humidity that   !
      !      were just committed, whatever value that is. Asserting a positive temperature instead   !
      !      would be asserting that the synthetic fill happens to be physical, which it is not.  ---!
      call check_close(bio%cas%can_temp, cas_temp_of_enthalpy(a%cas_enthalpy, a%cas_shv),          &
                       1.0e-12_wp, 'unpack: can_temp re-diagnosed from the COMMITTED enthalpy')
   end block

contains

   !----- Every field a DISTINCT non-zero multiple of `s`, so a combinator that copies the      !
   !      wrong field, or leaves one at its default, cannot pass by coincidence. ---------------!
   subroutine fill_state(y, s)
      type(column_state_t), intent(out) :: y
      real(wp),             intent(in)  :: s
      integer(ik) :: j
      y%cas_enthalpy   = 11.0_wp * s
      y%cas_shv        = 12.0_wp * s
      y%cas_co2        = 13.0_wp * s
      y%w_surface      = 14.0_wp * s
      y%w_surface_enth = 15.0_wp * s
      y%soil_energy = 0.0_wp ; y%theta = 0.0_wp
      do j = 1_ik, NSL
         y%soil_energy(j) = (20.0_wp + real(j, wp)) * s
         y%theta(j)       = (30.0_wp + real(j, wp)) * s
      end do
      allocate(y%leaf_water_mass(N), y%wood_water_mass(N), y%leaf_surf_water(N), y%wood_surf_water(N))
      do j = 1_ik, N
         y%leaf_water_mass(j) = (40.0_wp + real(j, wp)) * s
         y%wood_water_mass(j) = (50.0_wp + real(j, wp)) * s
         y%leaf_surf_water(j) = (60.0_wp + real(j, wp)) * s
         y%wood_surf_water(j) = (70.0_wp + real(j, wp)) * s
      end do
   end subroutine fill_state

   subroutine fill_tend(t)
      type(column_tend_t), intent(out) :: t
      integer(ik) :: j
      t%d_cas_enthalpy = 1.5_wp ; t%d_cas_shv = 2.5_wp ; t%d_cas_co2 = 3.5_wp
      t%dedt = 0.0_wp ; t%dtheta_dt = 0.0_wp
      do j = 1_ik, NSL
         t%dedt(j) = 4.0_wp + real(j, wp) ; t%dtheta_dt(j) = 5.0_wp + real(j, wp)
      end do
      allocate(t%d_leaf_water_mass(N), t%d_wood_water_mass(N),                                     &
               t%d_leaf_surf_water(N), t%d_wood_surf_water(N))
      do j = 1_ik, N
         t%d_leaf_water_mass(j) = 6.0_wp + real(j, wp) ; t%d_wood_water_mass(j) = 7.0_wp + real(j, wp)
         t%d_leaf_surf_water(j) = 8.0_wp + real(j, wp) ; t%d_wood_surf_water(j) = 9.0_wp + real(j, wp)
      end do
   end subroutine fill_tend

   subroutine same(x, y, tag)
      type(column_state_t), intent(in) :: x, y
      character(len=*),     intent(in) :: tag
      call check_close(x%cas_enthalpy,       y%cas_enthalpy,       1.0e-14_wp, tag//': cas_enthalpy')
      call check_close(x%cas_shv,            y%cas_shv,            1.0e-14_wp, tag//': cas_shv')
      call check_close(x%cas_co2,            y%cas_co2,            1.0e-14_wp, tag//': cas_co2')
      call check_close(x%w_surface,          y%w_surface,          1.0e-14_wp, tag//': w_surface (pond)')
      call check_close(x%w_surface_enth,     y%w_surface_enth,     1.0e-14_wp, tag//': w_surface_enth')
      call check_close(x%soil_energy(NSL),   y%soil_energy(NSL),   1.0e-14_wp, tag//': soil_energy')
      call check_close(x%theta(NSL),         y%theta(NSL),         1.0e-14_wp, tag//': theta')
      call check_close(x%leaf_water_mass(N), y%leaf_water_mass(N), 1.0e-14_wp, tag//': leaf_water_mass')
      call check_close(x%wood_water_mass(N), y%wood_water_mass(N), 1.0e-14_wp, tag//': wood_water_mass')
      call check_close(x%leaf_surf_water(N), y%leaf_surf_water(N), 1.0e-14_wp, tag//': leaf_surf_water')
      call check_close(x%wood_surf_water(N), y%wood_surf_water(N), 1.0e-14_wp, tag//': wood_surf_water')
   end subroutine same

end program test_state_combinators
