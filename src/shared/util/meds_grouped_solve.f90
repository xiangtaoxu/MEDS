! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_grouped_solve -- a linear solve for systems whose unknowns split into a shared core and  !
! groups that each touch only themselves and the core, the core through their FIRST unknown:  !
!                                                                                          !
!      [ A        b_1 e_1'  ...  b_G e_1' ] [x_0]   [r_0]                                  !
!      [ e_1 c_1'   D_1                   ] [x_1] = [r_1]                                  !
!      [ ...               ...            ] [...]   [...]                                  !
!      [ e_1 c_G'                    D_G  ] [x_G]   [r_G]                                  !
!                                                                                          !
! b_g is the core's column for group g's first unknown, c_g that unknown's row over the core,    !
! and e_1 picks the group's first unknown. Each group's block is eliminated onto the core: the   !
! Schur update is the rank-1 product (D_g^-1)_11 b_g c_g', so the cost grows as groups x core^2  !
! + core^3, not as (core + groups)^3. The core system that remains is solved by LU with partial  !
! pivoting, and the groups are recovered from it.                                                 !
!                                                                                          !
! The soil-plant water network is one such system -- soil layers and the pond are the core, each !
! cohort a group whose wood meets the soil and whose leaves meet only the wood (a root node, when  !
! it comes, will be a group's first unknown in the wood's place) -- but nothing here knows that    !
! (MEDS_SPAC_WATER_NETWORK_PLAN.md §4.3). Factor once, then solve for as many right-hand sides as  !
! needed: the network's error estimate reuses the last Newton factorisation.                       !
!==========================================================================================!
module meds_grouped_solve
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: grouped_system_t, grouped_allocate, grouped_copy, grouped_factor, grouped_solve

   !----- One system: the blocks, filled by the caller, and the factorisation that replaces them. -!
   type :: grouped_system_t
      integer(ik) :: nc = 0_ik            !< core unknowns
      integer(ik) :: ng = 0_ik            !< groups
      integer(ik) :: mg = 0_ik            !< unknowns per group
      real(wp),    allocatable :: a(:,:)    !< (nc, nc)     core block; on factor, the LU of the Schur complement
      real(wp),    allocatable :: b(:,:)    !< (nc, ng)     core rows, each group's first column
      real(wp),    allocatable :: c(:,:)    !< (nc, ng)     each group's first row, over the core columns
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
         allocate(s%a(nc, nc), s%b(nc, ng), s%c(nc, ng), s%d(mg, mg, ng), s%piv(nc))
         s%nc = nc ; s%ng = ng ; s%mg = mg
      end if
      s%a = 0.0_wp ; s%b = 0.0_wp ; s%c = 0.0_wp ; s%d = 0.0_wp
   end subroutine grouped_allocate

   !----- Copy one system's blocks into another of the same shape, in place. --------------------!
   pure subroutine grouped_copy(dst, src)
      type(grouped_system_t), intent(inout) :: dst
      type(grouped_system_t), intent(in)    :: src
      dst%a = src%a ; dst%b = src%b ; dst%c = src%c ; dst%d = src%d
   end subroutine grouped_copy

   !----- Factor in place: invert each group block, fold the groups into the core by their rank-1  !
   !      updates, and LU-factor the core. ok is false if a group block or the core is singular. ---!
   pure subroutine grouped_factor(s, ok)
      type(grouped_system_t), intent(inout) :: s
      logical,                intent(out)   :: ok
      integer(ik) :: g, i, j
      real(wp)    :: e11, cj
      ok = .true.
      do g = 1_ik, s%ng
         call invert_small(s%d(:, :, g), s%mg, ok)
         if (.not. ok) return
         e11 = s%d(1, 1, g)
         do j = 1_ik, s%nc
            cj = e11 * s%c(j, g)
            if (cj == 0.0_wp) cycle
            do i = 1_ik, s%nc
               s%a(i, j) = s%a(i, j) - s%b(i, g) * cj
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
      real(wp)    :: rc(s%nc), cx
      integer(ik) :: g, i, l
      rc = r_core(1:s%nc)
      !----- Each group's own solution, D_g^-1 r_g, folded onto the core through its first unknown. -!
      do g = 1_ik, s%ng
         do i = 1_ik, s%mg
            x_group(i, g) = 0.0_wp
            do l = 1_ik, s%mg
               x_group(i, g) = x_group(i, g) + s%d(i, l, g) * r_group(l, g)
            end do
         end do
         if (x_group(1, g) /= 0.0_wp) rc = rc - s%b(:, g) * x_group(1, g)
      end do
      call lu_solve(s%a, s%piv, s%nc, rc)
      x_core(1:s%nc) = rc
      !----- Back to the groups: x_g = D_g^-1 r_g - D_g^-1 e_1 (c_g' x_0). --------------------------!
      do g = 1_ik, s%ng
         cx = sum(s%c(:, g) * rc)
         x_group(1:s%mg, g) = x_group(1:s%mg, g) - s%d(1:s%mg, 1, g) * cx
      end do
   end subroutine grouped_solve

   !----- Invert a small dense block in place: closed form for 1x1 and 2x2, Gauss-Jordan with      !
   !      partial pivoting otherwise. ok is false if it is singular. -------------------------------!
   pure subroutine invert_small(m, n, ok)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: m(n, n)
      logical,     intent(out)   :: ok
      real(wp)    :: det, m11
      ok = .true.
      if (n == 1_ik) then
         if (m(1, 1) == 0.0_wp .or. m(1, 1) /= m(1, 1)) then
            ok = .false. ; return
         end if
         m(1, 1) = 1.0_wp / m(1, 1)
      else if (n == 2_ik) then
         det = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)
         if (det == 0.0_wp .or. det /= det) then
            ok = .false. ; return
         end if
         m11 = m(1, 1)
         m(1, 1) =  m(2, 2) / det ; m(2, 2) = m11 / det
         m(1, 2) = -m(1, 2) / det ; m(2, 1) = -m(2, 1) / det
      else
         call invert_gauss_jordan(m, n, ok)
      end if
   end subroutine invert_small

   pure subroutine invert_gauss_jordan(m, n, ok)
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
   end subroutine invert_gauss_jordan

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
