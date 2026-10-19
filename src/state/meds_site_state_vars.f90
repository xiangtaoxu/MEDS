! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_site_state_vars -- the sites of a polygon: each site's area and its slice of the     !
! patch block.                                                                             !
!                                                                                          !
! A site is the part of a polygon with one soil and one terrain. Its patches are one        !
! contiguous slice of the polygon's patch block, which is sorted by site and then by age,   !
! and the site block maps each site to its slice with an offset and a count, as the patch   !
! block maps each patch to its cohorts. A patch's area is a fraction of its site, and a     !
! site's area a fraction of the polygon.                                                    !
!==========================================================================================!
module meds_site_state_vars
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: site_block, alloc_site_block, free_site_block

   !----- All sites of a polygon: an array per number. --------------------------------------!
   type :: site_block
      !==========================================================================================!
      ! 0. BLOCK RECORD                                                                          !
      !==========================================================================================!
      integer(ik) :: n   = 0_ik       !< sites in use
      integer(ik) :: cap = 0_ik       !< allocated capacity
      !==========================================================================================!
      ! 1. AREA                                                                                  !
      !==========================================================================================!
      real(wp),    allocatable :: area(:)           !< fraction of the polygon, sum = 1
      !==========================================================================================!
      ! 2. PATCH MAP: a CSR map into the patch block (rebuild_patch_map derives it)              !
      !==========================================================================================!
      integer(ik), allocatable :: patch_offset(:)   !< first patch of the site (1-based)
      integer(ik), allocatable :: patch_count(:)    !< patches in the site
   end type site_block

contains

   !=======================================================================================!
   !  THE SITE LOCKSTEP: the routines that list every site field, in section order         !
   !=======================================================================================!

   subroutine alloc_site_block(site, cap)
      type(site_block), intent(inout) :: site
      integer(ik),      intent(in)    :: cap
      site%cap = cap
      site%n = 0_ik
      !----- 1. Area ----------------------------------------------------------------------------!
      allocate(site%area(cap))
      site%area = 0.0_wp
      !----- 2. Patch map -----------------------------------------------------------------------!
      allocate(site%patch_offset(cap), site%patch_count(cap))
      site%patch_offset = 0_ik
      site%patch_count = 0_ik
   end subroutine alloc_site_block

   !----- Free every per-site array. ------------------------------------------------------!
   subroutine free_site_block(site)
      type(site_block), intent(inout) :: site
      site%n = 0_ik
      site%cap = 0_ik
      if (allocated(site%area)) deallocate(site%area, site%patch_offset, site%patch_count)
   end subroutine free_site_block

end module meds_site_state_vars
