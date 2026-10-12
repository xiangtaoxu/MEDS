! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_grouped_solve -- a linear solve for systems whose unknowns split into a shared core and  !
! groups that each touch only the core and themselves:                                      !
!                                                                                          !
!      [ A    B_1  ...  B_G ] [x_0]   [r_0]                                                 !
!      [ C_1  D_1           ] [x_1] = [r_1]                                                 !
!      [ ...       ...      ] [...]   [...]                                                 !
!      [ C_G            D_G ] [x_G]   [r_G]                                                 !
!                                                                                          !
! Each group's small block is eliminated onto the core (a Schur complement); the core system  !
! that remains is solved by LU with partial pivoting; the groups are recovered from the core.  !
! The cost grows as groups x core^2 + core^3, not as (core + groups)^3. The soil-plant water     !
! network is one such system -- soil layers and the pond are the core, each cohort a group --    !
! but nothing here knows that (MEDS_SPAC_WATER_NETWORK_PLAN.md §4.3).                              !
!                                                                                          !
! Factor once, then solve for as many right-hand sides as needed: the network's error estimate  !
! reuses the last Newton factorisation.                                                          !
!==========================================================================================!
module meds_grouped_solve
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: grouped_system_t, grouped_allocate, grouped_factor, grouped_solve

   !----- One system: the blocks, filled by the caller, and the factorisation that replaces them. -!
   type :: grouped_system_t
      integer(ik) :: nc = 0_ik            !< core unknowns
      integer(ik) :: ng = 0_ik            !< groups
      integer(ik) :: mg = 0_ik            !< unknowns per group
      real(wp),    allocatable :: a(:,:)    !< (nc, nc)     core block; on factor, the LU of the Schur complement
      real(wp),    allocatable :: b(:,:,:)  !< (nc, mg, ng) core rows, group columns
      real(wp),    allocatable :: c(:,:,:)  !< (mg, nc, ng) group rows, core columns; on factor, D^-1 C
      real(wp),    allocatable :: d(:,:,:)  !< (mg, mg, ng) group blocks; on factor, their inverses
      integer(ik), allocatable :: piv(:)    !< (nc) the core LU's row interchanges
   end type grouped_system_t

contains

   !----- Size a system for nc core unknowns and ng groups of mg, keeping the storage when the    !
   !      shape is unchanged. The blocks are zeroed. ---------------------------------------------!
   pure subroutine grouped_allocate(s, nc, ng, mg)
      type(grouped_system_t), intent(inout) :: s
      integer(ik),            intent(in)    :: nc, ng, mg
      if (s%nc /= nc .or. s%ng /= ng .or. s%mg /= mg .or. .not. allocated(s%a)) then
         if (allocated(s%a))   deallocate(s%a)
         if (allocated(s%b))   deallocate(s%b)
         if (allocated(s%c))   deallocate(s%c)
         if (allocated(s%d))   deallocate(s%d)
         if (allocated(s%piv)) deallocate(s%piv)
         allocate(s%a(nc, nc), s%b(nc, mg, ng), s%c(mg, nc, ng), s%d(mg, mg, ng), s%piv(nc))
         s%nc = nc ; s%ng = ng ; s%mg = mg
      end if
      s%a = 0.0_wp ; s%b = 0.0_wp ; s%c = 0.0_wp ; s%d = 0.0_wp
   end subroutine grouped_allocate

   !----- Factor in place: invert each group block, fold the groups into the core, and LU-factor   !
   !      the core. ok is false if a group block or the core is singular. -----------------------!
   pure subroutine grouped_factor(s, ok)
      type(grouped_system_t), intent(inout) :: s
      logical,                intent(out)   :: ok
      integer(ik) :: g, i, j, l
      real(wp)    :: dc(s%mg, s%nc)
      ok = .true.
      do g = 1_ik, s%ng
         call invert_small(s%d(:, :, g), s%mg, ok)
         if (.not. ok) return
         !----- D^-1 C, kept for the back-substitution, and the Schur update A -= B D^-1 C. -------!
         do j = 1_ik, s%nc
            do i = 1_ik, s%mg
               dc(i, j) = 0.0_wp
               do l = 1_ik, s%mg
                  dc(i, j) = dc(i, j) + s%d(i, l, g) * s%c(l, j, g)
               end do
            end do
         end do
         s%c(:, :, g) = dc
         do j = 1_ik, s%nc
            do l = 1_ik, s%mg
               if (dc(l, j) == 0.0_wp) cycle
               do i = 1_ik, s%nc
                  s%a(i, j) = s%a(i, j) - s%b(i, l, g) * dc(l, j)
               end do
            end do
         end do
      end do
      call lu_factor(s%a, s%piv, s%nc, ok)
   end subroutine grouped_factor

   !----- Solve the factored system for one right-hand side (r_core, r_group) -> (x_core, x_group). !
   pure subroutine grouped_solve(s, r_core, r_group, x_core, x_group)
      type(grouped_system_t), intent(in)  :: s
      real(wp),               intent(in)  :: r_core(:), r_group(:, :)   !< (nc), (mg, ng)
      real(wp),               intent(out) :: x_core(:), x_group(:, :)
      real(wp)    :: y(s%mg), rc(s%nc)
      integer(ik) :: g, i, l
      rc = r_core(1:s%nc)
      !----- Fold each group's right-hand side onto the core: r_0 - B_g D_g^-1 r_g. -------------!
      do g = 1_ik, s%ng
         do i = 1_ik, s%mg
            y(i) = 0.0_wp
            do l = 1_ik, s%mg
               y(i) = y(i) + s%d(i, l, g) * r_group(l, g)
            end do
         end do
         x_group(1:s%mg, g) = y
         do l = 1_ik, s%mg
            if (y(l) == 0.0_wp) cycle
            do i = 1_ik, s%nc
               rc(i) = rc(i) - s%b(i, l, g) * y(l)
            end do
         end do
      end do
      call lu_solve(s%a, s%piv, s%nc, rc)
      x_core(1:s%nc) = rc
      !----- Back to the groups: x_g = D_g^-1 r_g - D_g^-1 C_g x_0. -----------------------------!
      do g = 1_ik, s%ng
         do i = 1_ik, s%mg
            do l = 1_ik, s%nc
               x_group(i, g) = x_group(i, g) - s%c(i, l, g) * rc(l)
            end do
         end do
      end do
   end subroutine grouped_solve

   !----- Invert a small dense block in place by Gauss-Jordan with partial pivoting. ------------!
   pure subroutine invert_small(m, n, ok)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: m(n, n)
      logical,     intent(out)   :: ok
      real(wp)    :: aug(n, 2*n), piv_val, f, row(2*n)
      integer(ik) :: i, j, p
      aug = 0.0_wp
      aug(:, 1:n) = m
      do i = 1_ik, n
         aug(i, n + i) = 1.0_wp
      end do
      ok = .true.
      do j = 1_ik, n
         p = j
         do i = j + 1_ik, n
            if (abs(aug(i, j)) > abs(aug(p, j))) p = i
         end do
         piv_val = aug(p, j)
         if (piv_val == 0.0_wp .or. piv_val /= piv_val) then
            ok = .false. ; return
         end if
         if (p /= j) then
            row = aug(j, :) ; aug(j, :) = aug(p, :) ; aug(p, :) = row
         end if
         aug(j, :) = aug(j, :) / piv_val
         do i = 1_ik, n
            if (i == j) cycle
            f = aug(i, j)
            if (f /= 0.0_wp) aug(i, :) = aug(i, :) - f * aug(j, :)
         end do
      end do
      m = aug(:, n+1:2*n)
   end subroutine invert_small

   !----- Dense LU with partial pivoting, in place (Doolittle; L unit-lower, U upper). ----------!
   pure subroutine lu_factor(a, piv, n, ok)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: a(n, n)
      integer(ik), intent(out)   :: piv(n)
      logical,     intent(inout) :: ok
      integer(ik) :: i, j, k, p
      real(wp)    :: row(n), f
      do k = 1_ik, n
         p = k
         do i = k + 1_ik, n
            if (abs(a(i, k)) > abs(a(p, k))) p = i
         end do
         piv(k) = p
         if (a(p, k) == 0.0_wp .or. a(p, k) /= a(p, k)) then
            ok = .false. ; return
         end if
         if (p /= k) then
            row = a(k, :) ; a(k, :) = a(p, :) ; a(p, :) = row
         end if
         do i = k + 1_ik, n
            f = a(i, k) / a(k, k)
            a(i, k) = f
            if (f == 0.0_wp) cycle
            do j = k + 1_ik, n
               a(i, j) = a(i, j) - f * a(k, j)
            end do
         end do
      end do
   end subroutine lu_factor

   !----- Solve with lu_factor's result, in place. ------------------------------------------!
   pure subroutine lu_solve(a, piv, n, x)
      integer(ik), intent(in)    :: n
      real(wp),    intent(in)    :: a(n, n)
      integer(ik), intent(in)    :: piv(n)
      real(wp),    intent(inout) :: x(n)
      integer(ik) :: i, j
      real(wp)    :: t
      do i = 1_ik, n
         if (piv(i) /= i) then
            t = x(i) ; x(i) = x(piv(i)) ; x(piv(i)) = t
         end if
      end do
      do i = 2_ik, n
         do j = 1_ik, i - 1_ik
            x(i) = x(i) - a(i, j) * x(j)
         end do
      end do
      do i = n, 1_ik, -1_ik
         do j = i + 1_ik, n
            x(i) = x(i) - a(i, j) * x(j)
         end do
         x(i) = x(i) / a(i, i)
      end do
   end subroutine lu_solve

end module meds_grouped_solve
