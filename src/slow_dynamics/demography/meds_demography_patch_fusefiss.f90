! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_demography_patch_fusefiss -- PATCH numerical-resolution control for the adaptive patch          !
! discretization (the patch half of the ED2 fuse_fiss_utils analogue) plus patch DISTURBANCE.    !
! Driven by ARITHMETIC (fusion/termination) and by ECOLOGY (treefall disturbance), it keeps the  !
! patch count bounded and conserves site-level plant number / area.                              !
!                                                                                          !
!   Sorting        -- patches oldest-first (age descending); remap cohort ownership.            !
!   Patch fusion   -- patches with a similar vertical light structure (cumulative-LAI light     !
!                     profile, see patch_light_profile) and the same disturbance type fuse,     !
!                     conserving site-level plant number via area-fraction rescaling.           !
!   Termination    -- drop patches below the area floor; renormalize areas to sum to 1.         !
!   Disturbance    -- apply_patch_disturbance: treefall aggregates a fraction of every patch's   !
!                     area into ONE new age-0 gap (canopy dies, understorey survives).          !
! Depends on the cohort sibling meds_demography_cohort_fusefiss (patch fusion + disturbance re-sort    !
! the cohorts they touch).                                                                       !
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
                                       blend_soil_energy, blend_snow, snow_column_t, blend_soil_carbon,          &
                                       blend_litter_input, blend_xi_accum
   use meds_litter_partition, only : necromass_to_litter
   implicit none
   private

   public :: sort_patches, apply_patch_disturbance
   public :: new_fuse_patches, fuse_2_patches, terminate_patches, patch_light_profile

contains

   !---------------------------------------------------------------------------------------!
   ! Sort patches by site, then by age descending; remap cohort ownership and rebuild both !
   ! maps.                                                                                  !
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
   ! Cumulative-LAI LIGHT profile for one patch (ED2 patch-fusion scheme): the light         !
   ! fraction reaching the top of each height layer, light(ihgt) = exp(-k LAI_at-or-above), !
   ! with k = [demography].light_ext_coeff_fusefiss. Layer 1 is the ground (most overlying   !
   ! LAI -> darkest), layer n_height_layers is the canopy top. Two patches with a similar     !
   ! light profile (and the same disturbance type) are fused. Mirrors ED2                     !
   ! patch_pft_size_profile + new_fuse_patches.                                               !
   !---------------------------------------------------------------------------------------!
   subroutine patch_light_profile(poly, cfg, ip, light)
      type(polygon_t),       intent(in)  :: poly
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: ip
      real(wp),            intent(out) :: light(:)         ! size cfg%n_height_layers
      real(wp)    :: cum_lai(cfg%n_height_layers)
      integer(ik) :: i, i0, i1, ihgt, nl
      nl = cfg%n_height_layers
      cum_lai = 0.0_wp
      i0 = poly%patch%cohort_offset(ip)
      i1 = i0 + poly%patch%cohort_count(ip) - 1_ik
      do i = i0, i1
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
   end subroutine patch_light_profile

   !---------------------------------------------------------------------------------------!
   ! Patch fusion with geometric tolerance relaxation, from patch_light_tol to             !
   ! patch_light_tol_max over n_patch_fusion_iter passes. It stops at the first pass that  !
   ! brings every site's count to max_patch. The ceiling is hard: patches more different   !
   ! than it stay apart even when the count is still above max_patch (forced fusion        !
   ! aside). Two patches fuse only within one site: their areas are fractions of it.       !
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

   !----- One patch-fusion sweep at a fixed tolerance, comparing light profiles. ----------!
   subroutine patch_fuse_pass(poly, cfg, tol, force)
      type(polygon_t),  intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: tol
      logical,             intent(in)    :: force
      logical,     allocatable :: alive(:)
      real(wp),    allocatable :: prof(:,:), lr(:), ld(:)
      real(wp)    :: davg, dmax, dif
      integer(ik) :: np, recp, donp, ihgt, npop, nl
      logical     :: similar, both_empty

      np = poly%patch%n
      if (np < 2_ik) return
      nl = cfg%n_height_layers
      allocate(alive(np), prof(nl, np), lr(nl), ld(nl))
      alive = .true.
      do recp = 1_ik, np
         call patch_light_profile(poly, cfg, recp, prof(:, recp))
      end do

      do recp = 1_ik, np - 1_ik
         if (.not. alive(recp)) cycle
         do donp = recp + 1_ik, np
            if (.not. alive(donp)) cycle
            if (poly%patch%owner_site(donp) /= poly%patch%owner_site(recp)) cycle
            if (poly%patch%dist_type(donp) /= poly%patch%dist_type(recp)) cycle
            both_empty = (poly%patch%cohort_count(recp) == 0_ik .and. poly%patch%cohort_count(donp) == 0_ik)
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
               call fuse_2_patches(poly, recp, donp)
               call rebuild_cohort_map(poly)          ! recp slice now holds merged cohorts
               call patch_light_profile(poly, cfg, recp, prof(:, recp))
               alive(donp) = .false.
            end if
         end do
      end do

      call patch_compact(poly, alive)
      call sort_cohorts(poly)
   end subroutine patch_fuse_pass

   !---------------------------------------------------------------------------------------!
   ! Merge donor patch into receptor, conserving the polygon's plant number via area weights.  !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_patches(poly, recp, donp)
      type(polygon_t), intent(inout) :: poly
      integer(ik),     intent(in)    :: recp, donp
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
         !----- The patch's fields, in the patch file's section order. The area is set last,      !
         !      below: it is the weight. ----------------------------------------------------------!
         !----- 1. Identity and area: area-weighted age. -------------------------------------------!
         patch%age(recp)            = rawgt * patch%age(recp)            + dawgt * patch%age(donp)
         !----- 3. Recruit bank. --------------------------------------------------------------------!
         patch%recruit_pool(:,recp) = rawgt * patch%recruit_pool(:,recp) + dawgt * patch%recruit_pool(:,donp)
         !----- 4-6. Area-weighted fast-biophysics reservoirs (conserves the area-extensive store);   !
         !      temp/fliq are re-diagnosed in the fast loop. -------------------------------------------!
         patch%cas(recp)    = blend_cas(rawgt,    patch%cas(recp),    dawgt, patch%cas(donp))
         patch%snow(recp)   = blend_snow(rawgt,   patch%snow(recp),   dawgt, patch%snow(donp))
         patch%soil_water(recp) = blend_soil_water(rawgt, patch%soil_water(recp), dawgt, patch%soil_water(donp))
         patch%soil_energy(recp) = blend_soil_energy(rawgt, patch%soil_energy(recp), dawgt, patch%soil_energy(donp))
         !----- 7. Area-weighted slow soil-carbon reservoir (conserves site-wide soil carbon). -------!
         patch%soil_carbon(recp) = blend_soil_carbon(rawgt, patch%soil_carbon(recp), dawgt, patch%soil_carbon(donp))
         !----- 8. shed_water_rate is a per-area RATE (like age): area-weighted, not nplant-weighted. -!
         patch%shed_water_rate(recp) = rawgt*patch%shed_water_rate(recp) + dawgt*patch%shed_water_rate(donp)
         patch%shed_water_enth(recp) = rawgt*patch%shed_water_enth(recp) + dawgt*patch%shed_water_enth(donp)
         patch%slow_co2_rate(recp)   = rawgt*patch%slow_co2_rate(recp)   + dawgt*patch%slow_co2_rate(donp)
         !----- 9. Today's litter is a per-area rate like the rest: area-weight it, or the fused patch !
         !      loses whatever the donor's cohorts shed this step.  ---------------------------------!
         patch%xi_accum(recp)    = blend_xi_accum(rawgt, patch%xi_accum(recp), dawgt, patch%xi_accum(donp))
         patch%litter_in(recp)   = blend_litter_input(rawgt, patch%litter_in(recp), dawgt, patch%litter_in(donp))
         !----- 10. adapt_dt_last is a controller SEED, not a conserved amount -- any value is valid and   !
         !      the controller re-adapts within a step. Area-weight it like its neighbours purely so    !
         !      the result is DETERMINISTIC and order-independent. ------------------------------------!
         patch%adapt_dt_last(recp)   = rawgt*patch%adapt_dt_last(recp)   + dawgt*patch%adapt_dt_last(donp)
         !----- 11. Diagnostic accumulators: every patch diagnostic is per m2 of ITS OWN ground, so!
         !      they all area-weight, exactly like age and shed_water_rate above.  ---------------!
         call blend_patch_diag(patch%diag, recp, donp, ar, ad)
         !----- Rescale receptor cohort densities (slice currently holds all recp cohorts). Their !
         !      per-plant amounts need nothing. --------------------------------------------------!
         i0 = patch%cohort_offset(recp) ; i1 = i0 + patch%cohort_count(recp) - 1_ik
         do i = i0, i1
            cohort%nplant(i) = cohort%nplant(i) * rawgt
         end do
         !----- Rescale and reassign donor cohorts to the receptor. -----------------------!
         i0 = patch%cohort_offset(donp) ; i1 = i0 + patch%cohort_count(donp) - 1_ik
         do i = i0, i1
            cohort%nplant(i)      = cohort%nplant(i) * dawgt
            cohort%owner_patch(i) = recp
         end do
         patch%area(recp) = anew
         patch%area(donp) = 0.0_wp
      end associate
   end subroutine fuse_2_patches

   !---------------------------------------------------------------------------------------!
   ! Remove patches with area below the floor, site by site; renormalize each site's areas !
   ! to sum to 1.                                                                          !
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
         !----- A patch too small to track is MERGED INTO the largest survivor, not deleted.    !
         !                                                                                     !
         !      Dropping it and renormalizing the survivors' areas back to 1 would silently     !
         !      redistribute the whole site: every conserved quantity would change by           !
         !      (a/(1-a)).sum_kept(a_i.X_i) - a.X_dropped, which vanishes only if the doomed    !
         !      patch held the survivors' mean -- and a doomed patch is atypical by             !
         !      construction, usually a fresh disturbance gap.                                  !
         !                                                                                     !
         !      Fusing instead conserves exactly, reuses the operator patch fusion already      !
         !      relies on, and says something defensible: a sliver of ground joins its biggest  !
         !      neighbour rather than evaporating. fuse_2_patches zeroes the donor's area, so   !
         !      the renormalisation below becomes a round-off no-op instead of a redistribution.!
         !      The neighbour is in the sliver's own site: the areas are fractions of it.       !
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
               call rebuild_cohort_map(poly) ! the survivor's CSR slice now holds the merged cohorts,
            end do                        ! and the NEXT fusion reads that slice -- as patch_fuse_pass does
         end do
         call patch_compact(poly, pkeep)
      end if

      !----- Renormalize each site's areas to sum to 1. After the merge above this corrects      !
      !      round-off only. --------------------------------------------------------------------!
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
   ! Compact the patch list by keep-mask: drop dropped patches AND their cohorts, remap      !
   ! owner_patch to the new patch indices, and rebuild both maps.                            !
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

      !----- Patches: keep the survivors, in order. The surviving-index list is the permutation   !
      !      reorder_patch_block takes, exactly as compact_cohort_block builds it for the cohort axis. --------!
      call reorder_patch_block(poly%patch, keep_idx, k)
      call rebuild_patch_map(poly)
      call rebuild_cohort_map(poly)
   end subroutine patch_compact

   !---------------------------------------------------------------------------------------!
   ! Treefall patch disturbance (ED2 analogue). A fraction f = 1 - exp(-rate*dt) of EVERY     !
   ! patch's area is disturbed and aggregated into one new age-0 gap patch per site           !
   ! (DIST_TREEFALL), cut from that site's patches: their areas are fractions of it.          !
   ! Survivorship follows ED2 treefall: tall canopy cohorts (height >= survive_height) die in !
   ! the gap; the short understory (height < threshold) survives into it. Each donor keeps    !
   ! area (1-f)*area with ALL its cohorts intact. Site area is conserved; plant number is     !
   ! conserved for survivors and reduced by f for the killed canopy -- that loss IS the        !
   ! disturbance. The gap is later consolidated by patch and cohort fusion.                   !
   !---------------------------------------------------------------------------------------!
   subroutine apply_patch_disturbance(poly, cfg, dt_yr, water_shed)
      type(polygon_t),       intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: dt_yr
      !----- What the KILLED CANOPY handed on, for the caller to declare to the slow ledger. The   !
      !      engine applies and reports; the driver accounts. -------------------------------------!
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

      !----- Report the disturbed-area flux, WRITTEN ON THE DONORS and written HERE, before          !
      !      patch%n grows: every patch loses the same fraction `frac`, and at this point the diag   !
      !      block's slots still line up with the donor patches. The gap then inherits the           !
      !      donors' blended slot (inherit_patch_diag, at the end of this routine), so               !
      !      Sum(area * value) counts the rate on all of the ground it was disturbed on.             !
      !      UNITS: the block's weight is Sum(dt) in SECONDS, so the contribution is the             !
      !      declared [1/yr] rate times this step's weight in seconds. `frac` is the fraction        !
      !      disturbed over dt_yr, so the rate is frac/dt_yr and this step's weight is               !
      !      dt_yr*yr_sec -- their product is frac*yr_sec, with dt_yr cancelling (the bare `frac`    !
      !      would be a per-SECOND rate under the per-year label).                                   !
      if (poly%patch%diag%active) then
         do i = 1_ik, np0
            poly%patch%diag%v(PA_DISTURB_AREA, i) = poly%patch%diag%v(PA_DISTURB_AREA, i)          &
                                                  + frac * yr_sec
         end do
      end if

      !----- Count understorey survivors (all current cohorts live in donor patches). ------!
      !                                                                                          !
      !      The same sweep records the DISTURBANCE pathway's mortality carbon off the other         !
      !      branch of the same test, so the split and the kill can never drift apart. Written on   !
      !      the DONORS and written HERE, before patch%n grows -- the same rule PA_DISTURB_AREA     !
      !      above follows.                                                                        !
      !                                                                                          !
      !      Valued on the donor's ORIGINAL ground. Every patch diagnostic is per m2 of its OWN      !
      !      ground and is aggregated to the site as Sum(area * value) at READ time. The donors are  !
      !      shrunk to (1-frac)*area by then, and the gap -- the frac part the carbon stood on --    !
      !      inherits the donors' slot, so together they weigh the value by the donors' original     !
      !      area and the site total is exact.                                                       !
      !      Rate [kgC/m2/yr] x weight in seconds is the block's contract; dt_yr cancels             !
      !      against the weight dt_yr*yr_sec, as it does for PA_DISTURB_AREA.                        !
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
            !----- New age-0 gap patch (treefall). --------------------------------------------!
            patch%area(newp)           = new_area
            patch%age(newp)            = 0.0_wp
            patch%dist_type(newp)      = DIST_TREEFALL
            patch%owner_site(newp)     = is
            !----- Fresh gap, like age: it has no cohorts of its own that contributed to today's shed  !
            !      rate (unlike soil_carbon/xi_accum, which are genuine inherited material). ----------!
            patch%shed_water_rate(newp) = 0.0_wp
            patch%shed_water_enth(newp) = 0.0_wp
            patch%slow_co2_rate(newp)   = 0.0_wp   ! fresh gap: no cohort of its own respired today
            patch%adapt_dt_last(newp)   = 0.0_wp   ! a fresh gap cold-starts the controller
            patch%n = newp

            !----- Seed the gap's fast reservoirs = area-weighted donor mean (the soil column and !
            !      canopy air are INHERITED from the disturbed area, not created bare). ----------!
            atot = sum(patch%area(p0:p1))                       ! donors still hold pre-disturbance area here
            if (atot > tiny_num) then
               !----- The SEED BANK is inherited material too. A treefall gap does not sterilise      !
               !      the ground it opens: the carry-forward recruit pool sits in the soil with the   !
               !      litter and the CENTURY carbon, both of which the gap inherits. Zeroing it would !
               !      destroy the pool's carbon on the disturbed fraction. ---------------------------!
               patch%recruit_pool(:,newp) = (patch%area(p0)/atot) * patch%recruit_pool(:,p0)
               patch%cas(newp)    = blend_cas(   patch%area(p0)/atot, patch%cas(p0),    0.0_wp, patch%cas(p0))
               patch%snow(newp)   = blend_snow(  patch%area(p0)/atot, patch%snow(p0),   0.0_wp, patch%snow(p0))
               patch%soil_water(newp) = blend_soil_water(patch%area(p0)/atot, patch%soil_water(p0), 0.0_wp,   &
                                                         patch%soil_water(p0))
               patch%soil_energy(newp) = blend_soil_energy(patch%area(p0)/atot, patch%soil_energy(p0), 0.0_wp, &
                                                           patch%soil_energy(p0))
               patch%soil_carbon(newp) = blend_soil_carbon(patch%area(p0)/atot, patch%soil_carbon(p0), &
                                                           0.0_wp, patch%soil_carbon(p0))
               patch%xi_accum(newp) = blend_xi_accum(patch%area(p0)/atot, patch%xi_accum(p0), 0.0_wp,       &
                                                     patch%xi_accum(p0))
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
                  patch%xi_accum(newp)    = blend_xi_accum(1.0_wp, patch%xi_accum(newp), wd, patch%xi_accum(d))
                  patch%litter_in(newp)   = blend_litter_input(1.0_wp, patch%litter_in(newp), wd,        &
                                                               patch%litter_in(d))
               end do
            end if

            !----- Move understorey survivors into the gap at area-weighted density; the killed           !
            !      canopy's carbon becomes litter in the SAME gap patch when soil carbon is on            !
            !      ([soil_carbon].soil_carbon_on, default on; docs/science/soil_carbon.md §5): the        !
            !      density it would have carried into the gap had it survived (the conversion factor      !
            !      the survivor branch below uses), times its per-plant carbon pools. Added directly onto !
            !      soil_carbon(newp), since this module cannot link biogeochemistry (necromass_to_litter  !
            !      takes plain scalars). -----------------------------------------------------------------!
            m = cohort%n
            do d = p0, p1
               i0 = patch%cohort_offset(d) ; i1 = i0 + patch%cohort_count(d) - 1_ik
               do i = i0, i1
                  if (cohort%height(i) >= cfg%disturbance_survive_height) then   ! canopy dies in gap
                     !----- The killed canopy's WATER goes to the gap's ground down the same channel   !
                     !      turnover shedding uses; this is its largest single use -- a whole canopy   !
                     !      dies at once. Its tissue HEAT is a thermal-mass change, declared with the  !
                     !      rest of the phase's by the caller. `kill_frac` is the fraction of the      !
                     !      cohort's area-density the gap takes -- the SAME factor the survivor branch !
                     !      uses. ---------------------------------------------------------------------!
                     kill_frac = frac * patch%area(d)
                     w_kill = kill_frac * (cohort_tissue_water(cohort, i) + cohort_film_water(cohort, i))
                     patch%shed_water_rate(newp) = patch%shed_water_rate(newp)                         &
                                                 + w_kill / new_area / max(cfg%dt_slow, tiny_num)
                     patch%shed_water_enth(newp) = patch%shed_water_enth(newp)                         &
                                                 + kill_frac * cohort%nplant(i) * cohort_water_enthalpy(cohort, i) &
                                                   / new_area / max(cfg%dt_slow, tiny_num)
                     w_tot = w_tot + poly%site%area(is) * w_kill   ! w_kill is per m2 of the site
                     if (cfg%soil_carbon_on) then
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

         !----- Stamp the new gap patch with a fresh global id (the gap fragments are new            !
         !      entities). Its diagnostics: it inherits its donors' slots blended by area, as its    !
         !      reservoirs above inherit their state, so the stand total Sum(area * v) is unchanged  !
         !      by the carving. At a calendar boundary the donors hold only that boundary's events   !
         !      (the step's sums were read and reset), and the gap keeps its share of them; inside a !
         !      step (the C API) they hold the step so far. Every donor lost the same fraction, so   !
         !      their current areas give the same blend as their original ones. (expand_patch_block  !
         !      above grew the block with the patch arrays.) ----------------------------------------!
         call assign_patch_id(poly, newp)
         donor_area = 0.0_wp
         donor_area(p0:p1) = poly%patch%area(p0:p1)
         call inherit_patch_diag(poly%patch%diag, newp, donor_area)
      end do

      if (present(water_shed)) water_shed = w_tot

      !----- The moved-in survivor cohorts take fresh global ids (their donor cohorts keep their own). !
      do i = m0 + 1_ik, poly%cohort%n
         call assign_cohort_id(poly, i)
      end do

      call sort_patches(poly)                    ! each gap into its site's slice; both maps rebuilt
      call sort_cohorts(poly)
   end subroutine apply_patch_disturbance

end module meds_demography_patch_fusefiss
