! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_init_bare_ground -- init_mode 0, a near-bare-ground start, and the builders every   !
! start that makes cohorts shares (docs/science/order_of_processes.md §1, setup step 2;    !
! docs/configuration.md, How a run starts). init_bare_ground (one site) and                !
! init_bare_sites (several) lay out empty patches; seed_bare_ground plants each PFT's      !
! nbg_init_density; add_cohort appends a cohort and finalize_init maps and sorts them.     !
!==========================================================================================!
module meds_init_bare_ground
   use meds_kinds,      only : wp, ik
   use meds_config,     only : meds_config_t, growth_window_steps
   use meds_allometry_lib, only : height_to_dbh
   use meds_cohort_state_vars,     only : expand_cohort_block, init_cohort
   use meds_site_state_vars,       only : make_site_surface
   use meds_polygon_state_vars,    only : polygon_t, alloc_polygon, rebuild_patch_map, rebuild_cohort_map,      &
                                          assign_cohort_id, assign_patch_id
   use meds_patch_state_vars,      only : patch_physical_state, DIST_PRIMARY
   use meds_demography_cohort_fusefiss, only : sort_cohorts
   implicit none
   private

   public :: init_bare_ground, init_bare_sites, seed_bare_ground, add_cohort, finalize_init
   public :: N_PATCH_INIT

   integer(ik), parameter :: N_PATCH_INIT = 6_ik     !< a bare-ground start's patches per site

contains

   !---------------------------------------------------------------------------------------!
   ! A bare polygon of one site: n_patch identical empty patches sharing the site's area.  !
   !---------------------------------------------------------------------------------------!
   subroutine init_bare_ground(poly, cfg, n_patch)
      type(polygon_t),     intent(out) :: poly
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: n_patch
      call init_bare_sites(poly, cfg, [n_patch], [cfg], [1.0_wp])
   end subroutine init_bare_ground

   !-------------------------------------------------------------------------------------------!
   ! A near-bare-ground polygon of several sites: site `is` takes its share `site_area(is)` of !
   ! the polygon, its surface from `site_cfg(is)` (its soil, from [soil] or its own maps) and  !
   ! n_patch(is) identical empty patches sharing its area. The patches are numbered site by    !
   ! site, so the block is sorted by site. Everything else comes from `cfg`, the run's.        !
   !-------------------------------------------------------------------------------------------!
   subroutine init_bare_sites(poly, cfg, n_patch, site_cfg, site_area)
      type(polygon_t),     intent(out) :: poly
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: n_patch(:)     !< patches per site
      type(meds_config_t), intent(in)  :: site_cfg(:)
      real(wp),            intent(in)  :: site_area(:)   !< [-] fractions of the polygon
      integer(ik) :: ip, is, k, ns
      ns = int(size(n_patch), ik)
      call alloc_polygon(poly, cfg%pft%n, coh_cap = 64_ik, pat_cap = max(sum(n_patch), 1_ik),  &
                      n_growth_window = growth_window_steps(cfg), pft = cfg%pft,                    &
                      n_site = ns)
      poly%patch%n = sum(n_patch)
      ip = 0_ik
      do is = 1_ik, ns
         call make_site_surface(poly%site%surface(is), site_cfg(is))
         poly%site%area(is) = site_area(is)
         do k = 1_ik, n_patch(is)
            ip = ip + 1_ik
            poly%patch%area(ip)       = 1.0_wp / real(n_patch(is), wp)
            poly%patch%age(ip)        = 0.0_wp
            poly%patch%dist_type(ip)  = DIST_PRIMARY
            poly%patch%owner_site(ip) = is
            call assign_patch_id(poly, ip)
         end do
      end do
      poly%patch%recruit_pool = 0.0_wp
      poly%cohort%n = 0_ik
      call rebuild_patch_map(poly)
      call rebuild_cohort_map(poly)
   end subroutine init_bare_sites

   !---------------------------------------------------------------------------------------!
   ! A near-bare-ground start's seedlings: every patch gets one cohort of each PFT with    !
   ! nbg_init_density > 0, that many plants per m2, at min_cohort_height on the PFT's own  !
   ! height curve, where recruits are born. All densities 0 leave the patches empty.       !
   !---------------------------------------------------------------------------------------!
   subroutine seed_bare_ground(poly, cfg)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp)    :: dbh
      integer(ik) :: ip, ipft
      if (.not. any(cfg%pft%nbg_init_density > 0.0_wp)) return
      do ip = 1_ik, poly%patch%n
         do ipft = 1_ik, cfg%pft%n
            if (cfg%pft%nbg_init_density(ipft) <= 0.0_wp) cycle
            dbh = height_to_dbh(poly%cohort%pft_table%allom(ipft), cfg%pft%min_cohort_height)
            call add_cohort(poly, cfg, ip, ipft, cfg%pft%nbg_init_density(ipft), dbh)
         end do
      end do
      call finalize_init(poly)
   end subroutine seed_bare_ground

   !---------------------------------------------------------------------------------------!
   ! Append one cohort to patch ip (caller invokes finalize_init afterwards).               !
   !---------------------------------------------------------------------------------------!
   subroutine add_cohort(poly, cfg, ip, ipft, nplant, dbh)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      integer(ik),         intent(in)    :: ip, ipft
      real(wp),            intent(in)    :: nplant, dbh
      integer(ik) :: m
      real(wp)    :: t_birth, psi_birth
      call expand_cohort_block(poly%cohort, poly%cohort%n + 1_ik)
      m = poly%cohort%n + 1_ik
      !----- Born in its patch's state (patch_physical_state); at setup the patch is not seeded yet, so !
      !      that is the run's starting state, and init_fast_reservoirs re-seeds it from the patch. -!
      call patch_physical_state(poly%patch%cas(ip), poly%patch%soil_water(ip), cfg,                   &
                                poly%site%surface(poly%patch%owner_site(ip)), t_birth, psi_birth)
      call init_cohort(poly%cohort, m, cfg%pft, ipft, ip, nplant, dbh, t_birth, psi_birth)
      call assign_cohort_id(poly, m)
      poly%cohort%n = m
   end subroutine add_cohort

   !----- The cohorts' map onto their patches, then each patch's cohorts sorted by height. --!
   subroutine finalize_init(poly)
      type(polygon_t), intent(inout) :: poly
      call rebuild_cohort_map(poly)
      call sort_cohorts(poly)
   end subroutine finalize_init

end module meds_init_bare_ground
