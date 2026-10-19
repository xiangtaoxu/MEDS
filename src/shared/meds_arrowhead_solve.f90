! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_arrowhead_solve -- a linear solve for a block arrowhead system: a shared core bordered by  !
! diagonal blocks that each touch only themselves and the core, the core through their FIRST    !
! unknown:                                                                                 !
!                                                                                          !
!      [ A        b_1 e_1'  ...  b_B e_1' ] [x_0]   [r_0]                                  !
!      [ e_1 c_1'   D_1                   ] [x_1] = [r_1]                                  !
!      [ ...               ...            ] [...]   [...]                                  !
!      [ e_1 c_B'                    D_B  ] [x_B]   [r_B]                                  !
!                                                                                          !
! b_j is the core's column for block j's first unknown, c_j that unknown's row over the core,    !
! and e_1 picks the block's first unknown. Each block is eliminated onto the core (a Schur       !
! complement): its update is the rank-1 product (D_j^-1)_11 b_j c_j', so the cost grows as        !
! blocks x core^2 + core^3, not as (core + blocks)^3. The core system that remains is solved by  !
! LU with partial pivoting, and the blocks are recovered from it. The factorisation is kept      !
! apart from the blocks, which it leaves as they were, and takes additions to the diagonals (a   !
! Newton matrix's storage terms), so the same blocks can be factored again with other additions. !
!                                                                                          !
! The soil-plant water solve (meds_spac_solver; docs/science/soil_plant_water.md §3) is one     !
! such system -- the pond and the soil layers are the core, each cohort a block whose wood      !
! meets the soil and whose leaves meet only the wood (a root store would be a block's first     !
! unknown in the wood's place; see docs/ROADMAP.md) -- but nothing here knows that. Factor      !
! once, then solve for as many right-hand sides as needed: the solver's error estimate          !
! reuses the last factorisation.                                                                !
!==========================================================================================!
module meds_arrowhead_solve
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: arrowhead_system_t, arrowhead_allocate, arrowhead_factor, arrowhead_solve

   !----- One system: the blocks, filled by the caller, and their factorisation. ---------------!
   type :: arrowhead_system_t
      integer(ik) :: nc = 0_ik            !< core unknowns
      integer(ik) :: nb = 0_ik            !< blocks
      integer(ik) :: mb = 0_ik            !< unknowns per block
      real(wp),    allocatable :: a(:,:)      !< (nc, nc)     core block
      real(wp),    allocatable :: b(:,:)      !< (nc, nb)     core rows, each block's first column
      real(wp),    allocatable :: c(:,:)      !< (nc, nb)     each block's first row, over the core columns
      real(wp),    allocatable :: d(:,:,:)    !< (mb, mb, nb) the diagonal blocks
      real(wp),    allocatable :: lu(:,:)     !< (nc, nc)     the LU of the Schur complement
      real(wp),    allocatable :: dinv(:,:,:) !< (mb, mb, nb) their inverses
      integer(ik), allocatable :: piv(:)      !< (nc) the core LU's row interchanges
   end type arrowhead_system_t

contains

   !----- Size a system for nc core unknowns and nb blocks of mb, keeping the storage when the    !
   !      shape is unchanged. The blocks are zeroed. ---------------------------------------------!
   pure subroutine arrowhead_allocate(s, nc, nb, mb)
      type(arrowhead_system_t), intent(inout) :: s
      integer(ik),              intent(in)    :: nc, nb, mb
      if (s%nc /= nc .or. s%nb /= nb .or. s%mb /= mb .or. .not. allocated(s%a)) then
         if (allocated(s%a))   deallocate(s%a, s%b, s%c, s%d, s%lu, s%dinv, s%piv)
         allocate(s%a(nc, nc), s%b(nc, nb), s%c(nc, nb), s%d(mb, mb, nb), s%lu(nc, nc), s%dinv(mb, mb, nb),  &
                  s%piv(nc))
         s%nc = nc ; s%nb = nb ; s%mb = mb
      end if
      s%a = 0.0_wp ; s%b = 0.0_wp ; s%c = 0.0_wp ; s%d = 0.0_wp
   end subroutine arrowhead_allocate

   !----- Factor the system whose blocks' diagonals carry the additions add_core and add_block:     !
   !      invert each block, fold the blocks into the core by their rank-1 updates, and            !
   !      LU-factor the core. The blocks are left as they were. ok is false if a block or the        !
   !      core is singular. ----------------------------------------------------------------------!
   pure subroutine arrowhead_factor(s, add_core, add_block, ok)
      type(arrowhead_system_t), intent(inout) :: s
      real(wp),                 intent(in)    :: add_core(:), add_block(:, :)   !< (nc), (mb, nb)
      logical,                  intent(out)   :: ok
      integer(ik) :: ib, i, j
      real(wp)    :: e11, cj
      ok = .true.
      s%lu = s%a
      do i = 1_ik, s%nc
         s%lu(i, i) = s%lu(i, i) + add_core(i)
      end do
      do ib = 1_ik, s%nb
         s%dinv(:, :, ib) = s%d(:, :, ib)
         do i = 1_ik, s%mb
            s%dinv(i, i, ib) = s%dinv(i, i, ib) + add_block(i, ib)
         end do
         call invert_small(s%dinv(:, :, ib), s%mb, ok)
         if (.not. ok) return
         e11 = s%dinv(1, 1, ib)
         do j = 1_ik, s%nc
            cj = e11 * s%c(j, ib)
            if (cj == 0.0_wp) cycle
            do i = 1_ik, s%nc
               s%lu(i, j) = s%lu(i, j) - s%b(i, ib) * cj
            end do
         end do
      end do
      call lu_factor(s%lu, s%piv, s%nc, ok)
   end subroutine arrowhead_factor

   !----- Solve the factored system for one right-hand side (r_core, r_block) -> (x_core, x_block). !
   pure subroutine arrowhead_solve(s, r_core, r_block, x_core, x_block)
      type(arrowhead_system_t), intent(in)  :: s
      real(wp),                 intent(in)  :: r_core(:), r_block(:, :)   !< (nc), (mb, nb)
      real(wp),                 intent(out) :: x_core(:), x_block(:, :)
      real(wp)    :: rc(s%nc), cx
      integer(ik) :: ib, i, l
      rc = r_core(1:s%nc)
      !----- Each block's own solution, D_j^-1 r_j, folded onto the core through its first unknown. -!
      do ib = 1_ik, s%nb
         do i = 1_ik, s%mb
            x_block(i, ib) = 0.0_wp
            do l = 1_ik, s%mb
               x_block(i, ib) = x_block(i, ib) + s%dinv(i, l, ib) * r_block(l, ib)
            end do
         end do
         if (x_block(1, ib) /= 0.0_wp) rc = rc - s%b(:, ib) * x_block(1, ib)
      end do
      call lu_solve(s%lu, s%piv, s%nc, rc)
      x_core(1:s%nc) = rc
      !----- Back to the blocks: x_j = D_j^-1 r_j - D_j^-1 e_1 (c_j' x_0). --------------------------!
      do ib = 1_ik, s%nb
         cx = sum(s%c(:, ib) * rc)
         x_block(1:s%mb, ib) = x_block(1:s%mb, ib) - s%dinv(1:s%mb, 1, ib) * cx
      end do
   end subroutine arrowhead_solve

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

end module meds_arrowhead_solve
