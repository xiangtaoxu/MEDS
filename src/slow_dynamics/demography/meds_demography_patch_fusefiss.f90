! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_demography_patch_fusefiss -- the patch operators (docs/science/demography.md §4).   !
! sort_patches keeps the patches by site, then oldest first. apply_patch_disturbance opens !
! one treefall gap per site; new_fuse_patches fuses patches of one site whose light        !
! profiles (eq. (10)) are alike; terminate_patches merges a patch too small to track into  !
! its site's largest. What each conserves and hands to the slow ledger is demography.md    !
! §5. They re-sort the cohorts with meds_demography_cohort_fusefiss.                       !
!==========================================================================================!
module meds_demography_patch_fusefiss
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : tiny_num, almost_one, yr_sec
   use meds_config,     only : meds_config_t
   use meds_patch_state_vars, only : DIST_TREEFALL
   use meds_cohort_state_vars, only : compact_cohort_block, expand_cohort_block, copy_cohort,  &
                                      cohort_tissue_water, cohort_film_water, cohort_water_enthalpy
   use meds_polygon_state_vars, only : polygon_t, rebuild_patch_map, rebuild_cohort_map, assign_cohort_id,      &
                                       assign_patch_id
   use meds_diag_vars,        only : blend_patch_diag, inherit_patch_diag, PA_DISTURB_AREA, PA_MORT_C_DISTURB
   use meds_demography_cohort_fusefiss, only : sort_cohorts
   use meds_patch_state_vars,   only : expand_patch_block, reorder_patch_block, blend_cas, blend_soil_water,      &
                                       blend_soil_energy, blend_snow, snow_state_t, blend_soil_carbon,          &
                                       blend_litter_input, blend_decomp_xi
   use meds_litter_partition, only : necromass_to_litter
   implicit none
   private

   public :: sort_patches, apply_patch_disturbance
   public :: new_fuse_patches, fuse_2_patches, terminate_patches, patch_light_profile

contains

   !---------------------------------------------------------------------------------------!
   ! Sort by site, then oldest first; remap owner_patch and rebuild both maps.             !
   !---------------------------------------------------------------------------------------!
   subroutine sort_patches(poly)
      type(polygon_t), intent(inout) :: poly
      integer(ik), allocatable :: pperm(:), inv(:)
      integer(ik) :: np, i, j, key, k
      np = poly%patch%n
      if (np < 2_ik) return
      allocate(pperm(np), inv(np))
      do i = 1_ik, np
         pperm(i) = i
      end do
      associate (patch => poly%patch, cohort => poly%cohort)
         !----- Insertion sort patch index by site, then by age descending. ----------------!
         do i = 2_ik, np
            key = pperm(i)
            j   = i - 1_ik
            do while (j >= 1_ik)
               if (patch%owner_site(pperm(j)) < patch%owner_site(key)) exit
               if (patch%owner_site(pperm(j)) == patch%owner_site(key) .and.                        &
                   patch%age(pperm(j)) >= patch%age(key)) exit
               pperm(j + 1_ik) = pperm(j)
               j = j - 1_ik
            end do
            pperm(j + 1_ik) = key
         end do
         !----- Apply the permutation to every patch field, the diagnostics included. ---------!
         call reorder_patch_block(patch, pperm, np)
         !----- Remap owner_patch: old index -> new position. -----------------------------!
         do k = 1_ik, np
            inv(pperm(k)) = k
         end do
         do i = 1_ik, cohort%n
            cohort%owner_patch(i) = inv(cohort%owner_patch(i))
         end do
      end associate
      call rebuild_patch_map(poly)
      call rebuild_cohort_map(poly)
   end subroutine sort_patches

   !---------------------------------------------------------------------------------------!
   ! A patch's light profile, eq. (10) of demography.md §4: the fraction of full sunlight  !
   ! left at each height layer, exp(-k x the LAI at and above it), with k =                !
   ! light_ext_coeff_fusefiss. Layer 1 is the ground, layer n_height_layers the top.       !
   !---------------------------------------------------------------------------------------!
   subroutine patch_light_profile(poly, cfg, ip, light)
      type(polygon_t),       intent(in)  :: poly
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: ip
      real(wp),            intent(out) :: light(:)         ! size cfg%n_height_layers
      call light_profile_of(poly, cfg, ip, poly%patch%cohort_offset(ip),                                   &
                            poly%patch%cohort_offset(ip) + poly%patch%cohort_count(ip) - 1_ik, light)
   end subroutine patch_light_profile

   !----- The light profile of patch ip's cohorts among positions i0..i1 of the cohort block: its !
   !      slice, or in a fusion pass the slice widened to the patches already fused into it. -----!
   subroutine light_profile_of(poly, cfg, ip, i0, i1, light)
      type(polygon_t),     intent(in)  :: poly
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: ip, i0, i1
      real(wp),            intent(out) :: light(:)         ! size cfg%n_height_layers
      real(wp)    :: cum_lai(cfg%n_height_layers)
      integer(ik) :: i, ihgt, nl
      nl = cfg%n_height_layers
      cum_lai = 0.0_wp
      do i = i0, i1
         if (poly%cohort%owner_patch(i) /= ip) cycle
         ihgt = min(nl, max(1_ik, count(cfg%height_edges < poly%cohort%height(i)) + 1_ik))
         cum_lai(ihgt) = cum_lai(ihgt) + poly%cohort%nplant(i) * poly%cohort%leaf_area(i)
      end do
      !----- Integrate top-down: cum_lai(ihgt) = LAI in this layer and every layer above. --!
      do ihgt = nl - 1_ik, 1_ik, -1_ik
         cum_lai(ihgt) = cum_lai(ihgt) + cum_lai(ihgt + 1_ik)
      end do
      do ihgt = 1_ik, nl
         light(ihgt) = exp(-cfg%light_ext_coeff_fusefiss * cum_lai(ihgt))
      end do
   end subroutine light_profile_of

   !---------------------------------------------------------------------------------------!
   ! Patch fusion (demography.md §4) in passes, the tolerance growing from                 !
   ! patch_light_tol to patch_light_tol_max, until no site holds more than max_patch       !
   ! patches. Patches less alike than the last tolerance stay apart even then, unless      !
   ! max_patch < 0 forces fusion. Only patches of one site fuse.                           !
   !---------------------------------------------------------------------------------------!
   subroutine new_fuse_patches(poly, cfg)
      type(polygon_t),  intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp)    :: tol
      integer(ik) :: it, maxp
      logical     :: force

      if (cfg%max_patch == 0_ik) return
      force = (cfg%max_patch < 0_ik)
      maxp  = abs(cfg%max_patch)
      tol   = cfg%patch_light_tol

      do it = 1_ik, cfg%n_patch_fusion_iter
         call patch_fuse_pass(poly, cfg, tol, force)
         if (.not. force .and. maxval(poly%site%patch_count(1:poly%site%n)) <= maxp) exit
         tol = tol * cfg%patch_light_tol_mult
      end do
   end subroutine new_fuse_patches

   !----- One fusion pass at a fixed tolerance: a receptor only absorbs patches after it, and an !
   !      absorbed patch absorbs nothing, so the cohort block can stay as it is until the pass   !
   !      ends. patch_compact then regroups it once. --------------------------------------------!
   subroutine patch_fuse_pass(poly, cfg, tol, force)
      type(polygon_t),  intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: tol
      logical,             intent(in)    :: force
      logical,     allocatable :: alive(:)
      real(wp),    allocatable :: prof(:,:), lr(:), ld(:)
      integer(ik), allocatable :: last(:), ncoh(:)   !< each patch's last cohort position and its count
      real(wp)    :: davg, dmax, dif
      integer(ik) :: np, recp, donp, ihgt, npop, nl
      logical     :: similar, both_empty

      np = poly%patch%n
      if (np < 2_ik) return
      nl = cfg%n_height_layers
      allocate(alive(np), prof(nl, np), lr(nl), ld(nl), last(np), ncoh(np))
      alive = .true.
      do recp = 1_ik, np
         ncoh(recp) = poly%patch%cohort_count(recp)
         last(recp) = poly%patch%cohort_offset(recp) + ncoh(recp) - 1_ik
         call light_profile_of(poly, cfg, recp, poly%patch%cohort_offset(recp), last(recp), prof(:, recp))
      end do

      do recp = 1_ik, np - 1_ik
         if (.not. alive(recp)) cycle
         do donp = recp + 1_ik, np
            if (.not. alive(donp)) cycle
            if (poly%patch%owner_site(donp) /= poly%patch%owner_site(recp)) cycle
            if (poly%patch%dist_type(donp) /= poly%patch%dist_type(recp)) cycle
            both_empty = (ncoh(recp) == 0_ik .and. ncoh(donp) == 0_ik)
            if (both_empty) then
               similar = .true.
            else
               lr = prof(:, recp) ; ld = prof(:, donp)
               davg = 0.0_wp ; dmax = 0.0_wp ; npop = 0_ik
               do ihgt = 1_ik, nl
                  !----- Skip layers above both canopies (both at full light). -----------!
                  if (lr(ihgt) >= almost_one .and. ld(ihgt) >= almost_one) cycle
                  dif  = abs(lr(ihgt) - ld(ihgt))
                  davg = davg + dif
                  dmax = max(dmax, dif)
                  npop = npop + 1_ik
               end do
               if (npop > 0_ik) davg = davg / real(npop, wp)
               similar = (davg <= tol) .and. (dmax <= tol * cfg%patch_light_maxdev_factor)
            end if
            if (force .or. similar) then
               call fuse_patch_into(poly, recp, donp, poly%patch%cohort_offset(recp), last(recp))
               last(recp) = max(last(recp), last(donp))
               ncoh(recp) = ncoh(recp) + ncoh(donp)
               call light_profile_of(poly, cfg, recp, poly%patch%cohort_offset(recp), last(recp), prof(:, recp))
               alive(donp) = .false.
            end if
         end do
      end do

      call patch_compact(poly, alive)
      call sort_cohorts(poly)
   end subroutine patch_fuse_pass

   !---------------------------------------------------------------------------------------!
   ! Merge patch donp into recp: the stores blend on area and the densities scale by each  !
   ! patch's share of the fused area, which conserves plant number (demography.md §4).     !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_patches(poly, recp, donp)
      type(polygon_t), intent(inout) :: poly
      integer(ik),     intent(in)    :: recp, donp
      call fuse_patch_into(poly, recp, donp, poly%patch%cohort_offset(recp),                               &
                           poly%patch%cohort_offset(recp) + poly%patch%cohort_count(recp) - 1_ik)
   end subroutine fuse_2_patches

   !----- fuse_2_patches, with the receptor's cohorts those it owns among positions r0..r1 (its slice, or a  !
   !      fusion pass's widened one); the donor has absorbed nothing, so its own slice holds its cohorts. ---!
   subroutine fuse_patch_into(poly, recp, donp, r0, r1)
      type(polygon_t), intent(inout) :: poly
      integer(ik),     intent(in)    :: recp, donp, r0, r1
      real(wp)    :: ar, ad, anew, rawgt, dawgt
      integer(ik) :: i, i0, i1

      associate (patch => poly%patch, cohort => poly%cohort)
         ar = patch%area(recp) ; ad = patch%area(donp) ; anew = ar + ad
         if (anew <= tiny_num) then                     ! degenerate: still hand donor cohorts to recp
            i0 = patch%cohort_offset(donp) ; i1 = i0 + patch%cohort_count(donp) - 1_ik
            do i = i0, i1
               cohort%owner_patch(i) = recp
            end do
            return
         end if
         rawgt = ar / anew ; dawgt = ad / anew
         !----- The fields in the patch file's section order; the area, the weight, is set last. ---!
         !----- 1. The age, area-weighted. ---------------------------------------------------------!
         patch%age(recp)            = rawgt * patch%age(recp)            + dawgt * patch%age(donp)
         !----- 3. The recruit bank. ------------------------------------------------------!
         patch%recruit_pool(:,recp) = rawgt * patch%recruit_pool(:,recp) + dawgt * patch%recruit_pool(:,donp)
         !----- 4-6. The canopy air, snow, soil water and soil heat, area-weighted, which conserves !
         !      each store; temperatures and liquid fractions are re-diagnosed in the fast loop. ---!
         patch%cas(recp)    = blend_cas(rawgt,    patch%cas(recp),    dawgt, patch%cas(donp))
         patch%snow(recp)   = blend_snow(rawgt,   patch%snow(recp),   dawgt, patch%snow(donp))
         patch%soil_water(recp) = blend_soil_water(rawgt, patch%soil_water(recp), dawgt, patch%soil_water(donp))
         patch%soil_energy(recp) = blend_soil_energy(rawgt, patch%soil_energy(recp), dawgt, patch%soil_energy(donp))
         !----- 7. The soil carbon, area-weighted. ----------------------------------------!
         patch%soil_carbon(recp) = blend_soil_carbon(rawgt, patch%soil_carbon(recp), dawgt, patch%soil_carbon(donp))
         !----- 8. The shed water, its enthalpy and the CO2 owed are per-area rates: area-weighted. ---!
         patch%shed_water_rate(recp) = rawgt*patch%shed_water_rate(recp) + dawgt*patch%shed_water_rate(donp)
         patch%shed_water_enth(recp) = rawgt*patch%shed_water_enth(recp) + dawgt*patch%shed_water_enth(donp)
         patch%slow_co2_rate(recp)   = rawgt*patch%slow_co2_rate(recp)   + dawgt*patch%slow_co2_rate(donp)
         !----- 9. The day's decomposition integral and litter are per-area amounts: area-weighted. ---!
         patch%decomp_xi(recp)    = blend_decomp_xi(rawgt, patch%decomp_xi(recp), dawgt, patch%decomp_xi(donp))
         patch%litter_in(recp)   = blend_litter_input(rawgt, patch%litter_in(recp), dawgt, patch%litter_in(donp))
         !----- 10. adapt_dt_last only seeds the step controller; area-weighted to be deterministic. ---!
         patch%adapt_dt_last(recp)   = rawgt*patch%adapt_dt_last(recp)   + dawgt*patch%adapt_dt_last(donp)
         !----- fast_cost orders the fast loop only: the fused patch carries both patches' work. --------!
         patch%fast_cost(recp)       = patch%fast_cost(recp) + patch%fast_cost(donp)
         !----- 11. The diagnostic rows, each per m2 of its own patch's ground: area-weighted. ---!
         call blend_patch_diag(patch%diag, recp, donp, ar, ad)
         !----- Scale the receptor's densities by its share of the fused area (per-plant amounts stay). ---!
         do i = r0, r1
            if (cohort%owner_patch(i) /= recp) cycle
            cohort%nplant(i) = cohort%nplant(i) * rawgt
         end do
         !----- Scale the donor's densities likewise and hand its cohorts to the receptor. ---!
         i0 = patch%cohort_offset(donp) ; i1 = i0 + patch%cohort_count(donp) - 1_ik
         do i = i0, i1
            cohort%nplant(i)      = cohort%nplant(i) * dawgt
            cohort%owner_patch(i) = recp
         end do
         patch%area(recp) = anew
         patch%area(donp) = 0.0_wp
      end associate
   end subroutine fuse_patch_into

   !---------------------------------------------------------------------------------------!
   ! Termination: merge each patch below min_patch_area into its site's largest surviving  !
   ! patch, then renormalize each site's areas to sum to 1 (demography.md §4).             !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_patches(poly, cfg)
      type(polygon_t),  intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      logical, allocatable :: pkeep(:)
      integer(ik)          :: np, ip, imax, is, p0, p1
      real(wp)             :: s

      np = poly%patch%n
      if (np < 1_ik) return
      allocate(pkeep(np))
      pkeep = poly%patch%area(1:np) >= cfg%min_patch_area
      !----- Every site keeps at least its largest patch. -------------------------------------!
      do is = 1_ik, poly%site%n
         p0 = poly%site%patch_offset(is)
         p1 = p0 + poly%site%patch_count(is) - 1_ik
         if (p1 < p0) cycle
         if (.not. any(pkeep(p0:p1))) then
            imax = p0 - 1_ik + maxloc(poly%patch%area(p0:p1), dim=1)
            pkeep(imax) = .true.
         end if
      end do

      if (.not. all(pkeep)) then
         !----- Merged, not deleted: dropping a sliver and renormalizing would redistribute every conserved !
         !      quantity across the site, while merging conserves exactly, so the renormalization below     !
         !      corrects round-off only (demography.md §4). The survivor is in the sliver's own site. ------!
         do is = 1_ik, poly%site%n
            p0 = poly%site%patch_offset(is)
            p1 = p0 + poly%site%patch_count(is) - 1_ik
            if (p1 < p0) cycle
            if (all(pkeep(p0:p1))) cycle
            imax = p0 - 1_ik + maxloc(poly%patch%area(p0:p1), dim=1)
            if (.not. pkeep(imax)) imax = p0 - 1_ik + maxloc(poly%patch%area(p0:p1), dim=1, mask=pkeep(p0:p1))
            do ip = p0, p1
               if (pkeep(ip) .or. ip == imax) cycle
               call fuse_2_patches(poly, imax, ip)
               call rebuild_cohort_map(poly) ! regroup after each merge: the next one reads the survivor's slice
            end do
         end do
         call patch_compact(poly, pkeep)
      end if

      !----- Renormalize each site's areas to sum to 1; after the merges, round-off only. ---!
      do is = 1_ik, poly%site%n
         p0 = poly%site%patch_offset(is)
         p1 = p0 + poly%site%patch_count(is) - 1_ik
         if (p1 < p0) cycle
         s = sum(poly%patch%area(p0:p1))
         if (s > tiny_num) then
            do ip = p0, p1
               poly%patch%area(ip) = poly%patch%area(ip) / s
            end do
         end if
      end do
   end subroutine terminate_patches

   !---------------------------------------------------------------------------------------!
   ! Drop the patches pkeep leaves out, with their cohorts; remap owner_patch to the new   !
   ! patch indices and rebuild both maps.                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine patch_compact(poly, pkeep)
      type(polygon_t), intent(inout) :: poly
      logical,         intent(in)    :: pkeep(:)
      integer(ik), allocatable :: newidx(:), keep_idx(:)
      logical,     allocatable :: ckeep(:)
      integer(ik) :: np, ip, k, i, nc

      np = poly%patch%n
      allocate(newidx(np), keep_idx(max(np, 1_ik)))
      k = 0_ik
      do ip = 1_ik, np
         if (pkeep(ip)) then
            k = k + 1_ik
            newidx(ip) = k
            keep_idx(k) = ip
         else
            newidx(ip) = 0_ik
         end if
      end do

      !----- Cohorts: keep only those whose patch survives; remap owner. ------------------!
      nc = poly%cohort%n
      allocate(ckeep(nc))
      do i = 1_ik, nc
         ckeep(i) = pkeep(poly%cohort%owner_patch(i))
         if (ckeep(i)) poly%cohort%owner_patch(i) = newidx(poly%cohort%owner_patch(i))
      end do
      call compact_cohort_block(poly%cohort, ckeep)

      !----- Patches: keep the survivors, in order; their indices are reorder_patch_block's permutation. ---!
      call reorder_patch_block(poly%patch, keep_idx, k)
      call rebuild_patch_map(poly)
      call rebuild_cohort_map(poly)
   end subroutine patch_compact

   !---------------------------------------------------------------------------------------!
   ! Treefall gaps (demography.md §4): a share f = 1 - exp(-rate x dt) of every patch's    !
   ! area becomes one age-0 gap per site. Cohorts at least disturbance_survive_height tall !
   ! die in it, shorter ones survive into it, and each donor keeps (1-f) of its area. Site !
   ! area is conserved; plant number is not, since the killed canopy is the disturbance.   !
   !---------------------------------------------------------------------------------------!
   subroutine apply_patch_disturbance(poly, cfg, dt_yr, water_shed)
      type(polygon_t),       intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: dt_yr
      !----- What the killed canopy handed on, for the caller to declare to the slow ledger. ---!
      real(wp), optional,  intent(out)   :: water_shed  !< [kg/m2 polygon] tissue + film water -> the gaps' ground
      integer(ik) :: np0, newp, d, i, i0, i1, m, m0, nsurv, pf, is, p0, p1
      real(wp)    :: frac, new_area, atot, wd
      real(wp), allocatable :: donor_area(:)
      real(wp)    :: lost_density, lab_g, lab_s, str_g, str_s, lig_g, lig_s
      real(wp)    :: kill_frac, w_kill, w_tot

      w_tot = 0.0_wp
      if (present(water_shed)) water_shed = 0.0_wp

      np0 = poly%patch%n
      if (np0 < 1_ik .or. cfg%patch_disturbance_rate <= 0.0_wp) return

      !----- Disturbed area fraction this step, and the aggregate new-patch area. ----------!
      frac = 1.0_wp - exp(-cfg%patch_disturbance_rate * dt_yr)
      if (frac <= tiny_num) return
      new_area = frac * sum(poly%patch%area(1:np0))
      if (new_area <= tiny_num) return

      !----- The disturbed area, written on the donors before the gaps are added; each gap inherits !
      !      its donors' rows (inherit_patch_diag, below). The rate frac/dt_yr [1/yr] times the     !
      !      weight dt_yr x yr_sec [s] is frac x yr_sec (diagnostics.md §4, "The slow rows"). ------!
      if (poly%patch%diag%active) then
         do i = 1_ik, np0
            poly%patch%diag%v(PA_DISTURB_AREA, i) = poly%patch%diag%v(PA_DISTURB_AREA, i)          &
                                                  + frac * yr_sec
         end do
      end if

      !----- Count the survivors. The other branch writes the disturbance's mortality carbon on the    !
      !      donors, per m2 of their original ground: the shrunk donors and the gap that inherits      !
      !      their rows weigh it by that area, so the site total is exact (dt_yr cancels as above). ---!
      nsurv = 0_ik
      do i = 1_ik, poly%cohort%n
         if (poly%cohort%height(i) < cfg%disturbance_survive_height) then
            nsurv = nsurv + 1_ik
         else if (poly%patch%diag%active) then
            d = poly%cohort%owner_patch(i)
            poly%patch%diag%v(PA_MORT_C_DISTURB, d) = poly%patch%diag%v(PA_MORT_C_DISTURB, d)      &
               + poly%cohort%nplant(i) * frac * yr_sec                                             &
                 * (poly%cohort%leaf_carbon(i)     + poly%cohort%fineroot_carbon(i)                &
                  + poly%cohort%wood_carbon(i)     + poly%cohort%nonstructural_carbon(i))
         end if
      end do

      call expand_patch_block(poly%patch, np0 + poly%site%n, poly%n_pft)
      call expand_cohort_block(poly%cohort, poly%cohort%n + nsurv)
      m0 = poly%cohort%n                         ! cohorts before the survivors are moved in
      allocate(donor_area(np0))

      !----- One gap per site, cut from that site's patches: their areas are fractions of it. The  !
      !      gaps go at the end of the block, and sort_patches below puts each in its site. -------!
      newp = np0
      do is = 1_ik, poly%site%n
         p0 = poly%site%patch_offset(is)
         p1 = p0 + poly%site%patch_count(is) - 1_ik
         if (p1 < p0) cycle
         new_area = frac * sum(poly%patch%area(p0:p1))
         if (new_area <= tiny_num) cycle
         newp = newp + 1_ik

         associate (cohort => poly%cohort, patch => poly%patch)
            !----- A new age-0 treefall gap. ----------------------------------------------!
            patch%area(newp)           = new_area
            patch%age(newp)            = 0.0_wp
            patch%dist_type(newp)      = DIST_TREEFALL
            patch%owner_site(newp)     = is
            !----- A fresh gap has no cohorts of its own yet, so none of its own rates (demography.md §4). ---!
            patch%shed_water_rate(newp) = 0.0_wp
            patch%shed_water_enth(newp) = 0.0_wp
            patch%slow_co2_rate(newp)   = 0.0_wp   ! fresh gap: no cohort of its own respired today
            patch%adapt_dt_last(newp)   = 0.0_wp   ! a fresh gap cold-starts the controller
            patch%fast_cost(newp)       = 0.0_wp   ! and is not yet measured
            patch%n = newp

            !----- The gap inherits the area-weighted mean of its donors' stores, the recruit pool among !
            !      them: a gap does not sterilise the ground it opens (demography.md §4). ---------------!
            atot = sum(patch%area(p0:p1))                       ! donors still hold pre-disturbance area here
            if (atot > tiny_num) then
               patch%recruit_pool(:,newp) = (patch%area(p0)/atot) * patch%recruit_pool(:,p0)
               patch%cas(newp)    = blend_cas(   patch%area(p0)/atot, patch%cas(p0),    0.0_wp, patch%cas(p0))
               patch%snow(newp)   = blend_snow(  patch%area(p0)/atot, patch%snow(p0),   0.0_wp, patch%snow(p0))
               patch%soil_water(newp) = blend_soil_water(patch%area(p0)/atot, patch%soil_water(p0), 0.0_wp,   &
                                                         patch%soil_water(p0))
               patch%soil_energy(newp) = blend_soil_energy(patch%area(p0)/atot, patch%soil_energy(p0), 0.0_wp, &
                                                           patch%soil_energy(p0))
               patch%soil_carbon(newp) = blend_soil_carbon(patch%area(p0)/atot, patch%soil_carbon(p0), &
                                                           0.0_wp, patch%soil_carbon(p0))
               patch%decomp_xi(newp) = blend_decomp_xi(patch%area(p0)/atot, patch%decomp_xi(p0), 0.0_wp,       &
                                                     patch%decomp_xi(p0))
               patch%litter_in(newp) = blend_litter_input(patch%area(p0)/atot, patch%litter_in(p0),     &
                                                          0.0_wp, patch%litter_in(p0))
               do d = p0 + 1_ik, p1
                  wd = patch%area(d) / atot
                  patch%recruit_pool(:,newp) = patch%recruit_pool(:,newp) + wd * patch%recruit_pool(:,d)
                  patch%cas(newp)    = blend_cas(   1.0_wp, patch%cas(newp),    wd, patch%cas(d))
                  patch%snow(newp)   = blend_snow(  1.0_wp, patch%snow(newp),   wd, patch%snow(d))  ! conserve snow into the gap
                  patch%soil_water(newp) = blend_soil_water(1.0_wp, patch%soil_water(newp), wd, patch%soil_water(d))
                  patch%soil_energy(newp) = blend_soil_energy(1.0_wp, patch%soil_energy(newp), wd, patch%soil_energy(d))
                  patch%soil_carbon(newp) = blend_soil_carbon(1.0_wp, patch%soil_carbon(newp), wd, patch%soil_carbon(d))
                  patch%decomp_xi(newp)    = blend_decomp_xi(1.0_wp, patch%decomp_xi(newp), wd, patch%decomp_xi(d))
                  patch%litter_in(newp)   = blend_litter_input(1.0_wp, patch%litter_in(newp), wd,        &
                                                               patch%litter_in(d))
               end do
            end if

            !----- Survivors move in at the density the disturbed share of their patch gives over the gap's !
            !      area; the killed canopy, at that density, becomes litter in the gap when soil carbon is  !
            !      on (soil_carbon.md §5). -----------------------------------------------------------------!
            m = cohort%n
            do d = p0, p1
               i0 = patch%cohort_offset(d) ; i1 = i0 + patch%cohort_count(d) - 1_ik
               do i = i0, i1
                  if (cohort%height(i) >= cfg%disturbance_survive_height) then   ! canopy dies in gap
                     !----- The killed canopy's water goes to the gap's ground down the channel turnover !
                     !      shedding uses (conservation.md §3); its tissue heat is a thermal-mass change !
                     !      the caller declares (conservation.md §4). kill_frac is the share of the      !
                     !      cohort's area density the gap takes, the factor the survivor branch uses. ---!
                     kill_frac = frac * patch%area(d)
                     w_kill = kill_frac * (cohort_tissue_water(cohort, i) + cohort_film_water(cohort, i))
                     patch%shed_water_rate(newp) = patch%shed_water_rate(newp)                         &
                                                 + w_kill / new_area / max(cfg%dt_slow, tiny_num)
                     patch%shed_water_enth(newp) = patch%shed_water_enth(newp)                         &
                                                 + kill_frac * cohort%nplant(i) * cohort_water_enthalpy(cohort, i) &
                                                   / new_area / max(cfg%dt_slow, tiny_num)
                     w_tot = w_tot + poly%site%area(is) * w_kill   ! w_kill is per m2 of the site
                     if (cfg%soil_carbon%soil_carbon_on) then
                        lost_density = cohort%nplant(i) * (frac * patch%area(d) / new_area)
                        pf = cohort%pft(i)
                        call necromass_to_litter(lost_density * cohort%leaf_carbon(i),                  &
                                 lost_density * cohort%fineroot_carbon(i),                               &
                                 lost_density * cohort%wood_carbon(i),                                   &
                                 lost_density * cohort%nonstructural_carbon(i),                          &
                                 cfg%pft%f_labile_leaf(pf), cfg%pft%f_labile_stem(pf),                   &
                                 cohort%pft_table%allom(pf)%aboveground_frac,                           &
                                 cfg%pft%struct_lignin_frac(pf),                                         &
                                 lab_g, lab_s, str_g, str_s, lig_g, lig_s)
                        patch%soil_carbon(newp)%fast_grnd_carbon   = patch%soil_carbon(newp)%fast_grnd_carbon   + lab_g
                        patch%soil_carbon(newp)%fast_soil_carbon   = patch%soil_carbon(newp)%fast_soil_carbon   + lab_s
                        patch%soil_carbon(newp)%struct_grnd_carbon = patch%soil_carbon(newp)%struct_grnd_carbon + str_g
                        patch%soil_carbon(newp)%struct_soil_carbon = patch%soil_carbon(newp)%struct_soil_carbon + str_s
                        patch%soil_carbon(newp)%struct_grnd_lignin  = patch%soil_carbon(newp)%struct_grnd_lignin  + lig_g
                        patch%soil_carbon(newp)%struct_soil_lignin  = patch%soil_carbon(newp)%struct_soil_lignin  + lig_s
                     end if
                     cycle
                  end if
                  m = m + 1_ik
                  call copy_cohort(cohort, m, i)
                  cohort%nplant(m)      = cohort%nplant(i) * (frac * patch%area(d) / new_area)
                  cohort%owner_patch(m) = newp
               end do
            end do
            cohort%n = m

            !----- Donors keep the undisturbed remainder (all their cohorts intact). ----------!
            do d = p0, p1
               patch%area(d) = (1.0_wp - frac) * patch%area(d)
            end do
         end associate

         !----- The gap takes a new id, and its donors' diagnostic rows blended on area, so the    !
         !      stand's Sum(area x value) is unchanged by the carving (diagnostics.md §4). Every   !
         !      donor lost the same share, so their current areas blend as their original ones. ---!
         call assign_patch_id(poly, newp)
         donor_area = 0.0_wp
         donor_area(p0:p1) = poly%patch%area(p0:p1)
         call inherit_patch_diag(poly%patch%diag, newp, donor_area)
      end do

      if (present(water_shed)) water_shed = w_tot

      !----- The survivors copied into the gaps take new ids; the cohorts they came from keep theirs. ---!
      do i = m0 + 1_ik, poly%cohort%n
         call assign_cohort_id(poly, i)
      end do

      call sort_patches(poly)                    ! each gap into its site's slice; both maps rebuilt
      call sort_cohorts(poly)
   end subroutine apply_patch_disturbance

end module meds_demography_patch_fusefiss
