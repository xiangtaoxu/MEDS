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
   use meds_fast_types,       only : column_state_t, column_tend_t
   use meds_column_state_ops, only : state_init, state_axpy, state_sub, state_err_diff, zero_like
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
