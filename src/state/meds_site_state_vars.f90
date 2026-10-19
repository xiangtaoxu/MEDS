! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_site_state_vars -- the sites of a polygon: each site's area and its slice of the     !
! patch block.                                                                             !
!                                                                                          !
! A site is the part of a polygon with one soil and one terrain. Its patches are one        !
! contiguous slice of the polygon's patch block, which is sorted by site and then by age,   !
! and the site block maps each site to its slice with an offset and a count, as the patch   !
! block maps each patch to its cohorts. A patch's area is a fraction of its site, and a     !
! site's area a fraction of the polygon. Each site has its own surface: its soil column and !
! the parameters derived from it, the column's bottom, and its terrain.                     !
!==========================================================================================!
module meds_site_state_vars
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : n_soil_layer_max
   use meds_config,             only : meds_config_t
   use meds_soil_column_config, only : soil_column_config_t, soil_params_t, soil_thermal_params_t,     &
                                       soil_column_params, soil_column_start_theta
   implicit none
   private

   public :: site_block, site_surface_t, alloc_site_block, free_site_block, make_site_surface

   !----- What makes a site: its soil and its terrain. It is configuration held on the site:   !
   !      never integrated, blended or reordered, and its parameters are derived once, when the !
   !      site is made (make_site_surface). ---------------------------------------------------!
   type :: site_surface_t
      !----- The soil column, from [soil] or a cell's soil maps: its layers and texture, and the   !
      !      parameters derived from them. ----------------------------------------------------!
      type(soil_column_config_t)  :: soil_column
      type(soil_params_t)         :: soil
      type(soil_thermal_params_t) :: soil_thermal
      real(wp)    :: theta_start(n_soil_layer_max) = 0.0_wp   !< [m3/m3] each layer's starting water
      !----- The column's bottom: the water boundary (SOIL_BC_*) with drainage's factor f of        !
      !      q = f K, and the temperature at the deep depth below the column. -------------------!
      integer(ik) :: water_bottom_bc = 0_ik
      real(wp)    :: drain_factor    = 1.0_wp      !< [-]
      real(wp)    :: deep_temp       = -1.0_wp     !< [K]; < 0 not set
      real(wp)    :: deep_depth      = 0.0_wp      !< [m], positive
      !----- The terrain: what corrects the polygon's forcing to the site. ---------------------!
      real(wp)    :: elevation_m = 0.0_wp          !< [m]
      real(wp)    :: slope_deg   = 0.0_wp          !< [degrees]
      real(wp)    :: aspect_deg  = 0.0_wp          !< [degrees from north, clockwise]
   end type site_surface_t

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
      !==========================================================================================!
      ! 3. SURFACE                                                                               !
      !==========================================================================================!
      type(site_surface_t), allocatable :: surface(:)
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
      !----- 3. Surface ------------------------------------------------------------------------!
      allocate(site%surface(cap))                                ! default-initialised by its type
   end subroutine alloc_site_block

   !----- Free every per-site array. ------------------------------------------------------!
   subroutine free_site_block(site)
      type(site_block), intent(inout) :: site
      site%n = 0_ik
      site%cap = 0_ik
      if (allocated(site%area)) deallocate(site%area, site%patch_offset, site%patch_count, site%surface)
   end subroutine free_site_block

   !=======================================================================================!
   !  A SITE'S SURFACE                                                                     !
   !=======================================================================================!

   !----- A site's surface from a configuration: its soil column, from [soil] or the soil     !
   !      maps of its cell (apply_map_inputs), and the column's bottom. The terrain is the      !
   !      caller's: it knows where the site is. --------------------------------------------!
   subroutine make_site_surface(surf, cfg)
      type(site_surface_t), intent(inout) :: surf
      type(meds_config_t),  intent(in)    :: cfg
      surf%soil_column = cfg%soil_column
      call soil_column_params(cfg%soil_column, surf%soil, surf%soil_thermal)
      call soil_column_start_theta(cfg%soil_column, surf%soil, cfg%init_soil_theta, surf%theta_start)
      surf%water_bottom_bc = cfg%soil%bottom_bc
      surf%drain_factor    = cfg%soil%drain_factor
      surf%deep_temp       = cfg%energy%deep_temp
      surf%deep_depth      = cfg%energy%deep_depth
   end subroutine make_site_surface

end module meds_site_state_vars
