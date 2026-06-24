!==========================================================================================!
! meds_sort -- ordering of cohorts (within each patch) and patches (across the site).      !
!                                                                                          !
! Cohort order is LOAD-BEARING: cohorts are kept tallest-first (height desc, DBH desc tie)  !
! so the overtopping-basal-area competition sweep and the fusion receptor/donor scan are    !
! well defined. Patch order is age-descending (oldest first), as in ED2 sort_patches.       !
!                                                                                          !
! Both sorts build a permutation and route it through the centralized lockstep machinery    !
! in meds_demography_types, so no per-field array is ever forgotten.                                  !
!==========================================================================================!
module meds_sort
   use meds_kinds, only : wp, ik
   use meds_demography_types, only : site, cohort_reorder, rebuild_csr
   implicit none
   private

   public :: sort_cohorts, sort_patches

contains

   !---------------------------------------------------------------------------------------!
   ! Sort cohorts within every patch slice: height descending, DBH descending on ties.      !
   !---------------------------------------------------------------------------------------!
   subroutine sort_cohorts(comm)
      type(site), intent(inout) :: comm
      integer(ik), allocatable :: perm(:)
      integer(ik) :: ip, i0, i1, i, j, key
      integer(ik) :: n

      n = comm%coh%n
      if (n < 1_ik) return
      allocate(perm(n))
      do i = 1_ik, n
         perm(i) = i
      end do
      associate (c => comm%coh, p => comm%pat)
         do ip = 1_ik, p%n
            i0 = p%coff(ip)
            i1 = i0 + p%ccount(ip) - 1_ik
            !----- Insertion sort of the index window perm(i0:i1), descending. ------------!
            do i = i0 + 1_ik, i1
               key = perm(i)
               j   = i - 1_ik
               do while (j >= i0)
                  if (higher(c%hite(perm(j)), c%dbh(perm(j)), c%hite(key), c%dbh(key))) exit
                  perm(j + 1_ik) = perm(j)
                  j = j - 1_ik
               end do
               perm(j + 1_ik) = key
            end do
         end do
      end associate
      call cohort_reorder(comm%coh, perm, n)
   end subroutine sort_cohorts

   !----- True if (h1,d1) ranks at least as high as (h2,d2) under the sort order. ----------!
   pure logical function higher(h1, d1, h2, d2)
      real(wp), intent(in) :: h1, d1, h2, d2
      if (h1 > h2) then
         higher = .true.
      else if (h1 < h2) then
         higher = .false.
      else
         higher = (d1 >= d2)
      end if
   end function higher

   !---------------------------------------------------------------------------------------!
   ! Sort patches by age descending; remap cohort ownership and rebuild the CSR map.        !
   !---------------------------------------------------------------------------------------!
   subroutine sort_patches(comm)
      type(site), intent(inout) :: comm
      integer(ik), allocatable :: pperm(:), inv(:)
      integer(ik) :: np, i, j, key, k
      np = comm%pat%n
      if (np < 2_ik) return
      allocate(pperm(np), inv(np))
      do i = 1_ik, np
         pperm(i) = i
      end do
      associate (p => comm%pat, c => comm%coh)
         !----- Insertion sort patch index by age descending. -----------------------------!
         do i = 2_ik, np
            key = pperm(i)
            j   = i - 1_ik
            do while (j >= 1_ik)
               if (p%age(pperm(j)) >= p%age(key)) exit
               pperm(j + 1_ik) = pperm(j)
               j = j - 1_ik
            end do
            pperm(j + 1_ik) = key
         end do
         !----- Apply permutation to patch arrays (RHS vector subscript -> safe temp). -----!
         p%area(1:np)           = p%area(pperm(1:np))
         p%age(1:np)            = p%age(pperm(1:np))
         p%dist_type(1:np)      = p%dist_type(pperm(1:np))
         p%avg_daily_temp(1:np) = p%avg_daily_temp(pperm(1:np))
         p%min_month_temp(1:np) = p%min_month_temp(pperm(1:np))
         p%recruit_pool(:,1:np) = p%recruit_pool(:,pperm(1:np))
         !----- Remap owner_patch: old index -> new position. -----------------------------!
         do k = 1_ik, np
            inv(pperm(k)) = k
         end do
         do i = 1_ik, c%n
            c%owner_patch(i) = inv(c%owner_patch(i))
         end do
      end associate
      call rebuild_csr(comm)
   end subroutine sort_patches

end module meds_sort
