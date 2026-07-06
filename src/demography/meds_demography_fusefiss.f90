!==========================================================================================!
! meds_demography_fusefiss -- NUMERICAL RESOLUTION control for the adaptive cohort/patch       !
! discretization (the ED2 fuse_fiss_utils analogue). Driven by ARITHMETIC, not ecology: it     !
! keeps the cohort/patch counts bounded and the representation conservative WITHOUT changing    !
! what the ecosystem does. The ecological EVENTS (growth, mortality, recruitment, disturbance)  !
! live in meds_demography_dynamics; this is their housekeeping counterpart -- sort/fuse/split/cull.!
!                                                                                          !
!   Sorting        -- cohorts tallest-first within each patch (height desc, DBH desc tie);    !
!                     patches oldest-first. The order is load-bearing for the overtopping     !
!                     sweeps and the receptor/donor fusion scans.                             !
!   Cohort fusion  -- receptor (taller) absorbs same-PFT donors of similar HEIGHT unless the  !
!   /fission          merged cohort LAI exceeds the cap; split a cohort whose LAI exceeds the  !
!                     cap. The conserved invariant is TOTAL ABOVEGROUND BIOMASS (carbon):      !
!                     fused/split diameters are re-derived from conserved AGB, never averaged. !
!   Termination    -- cull cohorts/patches below the density/area floors: a negligible entry    !
!                     to drop (also the numerical tail of ecological extinction).               !
!   Patch fusion   -- patches with a similar vertical light structure (cumulative-LAI light    !
!                     profile, see patch_light_profile) and the same disturbance type fuse,    !
!                     conserving site_t-level plant number via area-fraction rescaling.          !
!==========================================================================================!
module meds_demography_fusefiss
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : tiny_num, almost_one
   use meds_allometry,  only : agb_to_dbh, agb_c2, b2Ht, light_ext
   use meds_config,     only : meds_config_t
   use meds_demography_types, only : site_t, cohort_reorder, rebuild_csr, cohort_compact,        &
                                     cohort_ensure_capacity, copy_cohort_slot, set_cohort_size,  &
                                     assign_cohort_id
   use meds_column_state_types, only : blend_cas, blend_soil_w, blend_soil_e
   implicit none
   private

   public :: sort_cohorts, sort_patches
   public :: new_fuse_cohorts, fuse_2_cohorts, split_cohorts, terminate_cohorts, max_cohort_count
   public :: new_fuse_patches, fuse_2_patches, terminate_patches, patch_light_profile

   real(wp), parameter :: prof_relax = 1.5_wp     !< per-iteration patch-fusion tolerance growth

contains

   !=======================================================================================!
   !  Sorting                                                                               !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Sort cohorts within every patch slice: height descending, DBH descending on ties.      !
   !---------------------------------------------------------------------------------------!
   subroutine sort_cohorts(site)
      type(site_t), intent(inout) :: site
      integer(ik), allocatable :: perm(:)
      integer(ik) :: ip, i0, i1, i, j, key
      integer(ik) :: n

      n = site%cohort%n
      if (n < 1_ik) return
      allocate(perm(n))
      do i = 1_ik, n
         perm(i) = i
      end do
      associate (cohort => site%cohort, patch => site%patch)
         do ip = 1_ik, patch%n
            i0 = patch%cohort_offset(ip)
            i1 = i0 + patch%cohort_count(ip) - 1_ik
            !----- Insertion sort of the index window perm(i0:i1), descending. ------------!
            do i = i0 + 1_ik, i1
               key = perm(i)
               j   = i - 1_ik
               do while (j >= i0)
                  if (higher(cohort%height(perm(j)), cohort%dbh(perm(j)), cohort%height(key), cohort%dbh(key))) exit
                  perm(j + 1_ik) = perm(j)
                  j = j - 1_ik
               end do
               perm(j + 1_ik) = key
            end do
         end do
      end associate
      call cohort_reorder(site%cohort, perm, n)
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
   subroutine sort_patches(site)
      type(site_t), intent(inout) :: site
      integer(ik), allocatable :: pperm(:), inv(:)
      integer(ik) :: np, i, j, key, k
      np = site%patch%n
      if (np < 2_ik) return
      allocate(pperm(np), inv(np))
      do i = 1_ik, np
         pperm(i) = i
      end do
      associate (patch => site%patch, cohort => site%cohort)
         !----- Insertion sort patch index by age descending. -----------------------------!
         do i = 2_ik, np
            key = pperm(i)
            j   = i - 1_ik
            do while (j >= 1_ik)
               if (patch%age(pperm(j)) >= patch%age(key)) exit
               pperm(j + 1_ik) = pperm(j)
               j = j - 1_ik
            end do
            pperm(j + 1_ik) = key
         end do
         !----- Apply permutation to patch arrays (RHS vector subscript -> safe temp). -----!
         patch%area(1:np)           = patch%area(pperm(1:np))
         patch%age(1:np)            = patch%age(pperm(1:np))
         patch%dist_type(1:np)      = patch%dist_type(pperm(1:np))
         patch%recruit_pool(:,1:np) = patch%recruit_pool(:,pperm(1:np))
         patch%global_id(1:np)      = patch%global_id(pperm(1:np))
         patch%cas(1:np)            = patch%cas(pperm(1:np))
         patch%soil_e(1:np)         = patch%soil_e(pperm(1:np))
         patch%soil_w(1:np)         = patch%soil_w(pperm(1:np))
         !----- Remap owner_patch: old index -> new position. -----------------------------!
         do k = 1_ik, np
            inv(pperm(k)) = k
         end do
         do i = 1_ik, cohort%n
            cohort%owner_patch(i) = inv(cohort%owner_patch(i))
         end do
      end associate
      call rebuild_csr(site)
   end subroutine sort_patches

   !=======================================================================================!
   !  Cohort fusion / fission / termination                                                 !
   !=======================================================================================!

   integer(ik) function max_cohort_count(site)
      type(site_t), intent(in) :: site
      if (site%patch%n < 1_ik) then
         max_cohort_count = 0_ik
      else
         max_cohort_count = maxval(site%patch%cohort_count(1:site%patch%n))
      end if
   end function max_cohort_count

   !---------------------------------------------------------------------------------------!
   ! Cohort fusion with geometric tolerance relaxation.                                     !
   !---------------------------------------------------------------------------------------!
   subroutine new_fuse_cohorts(site, cfg)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp)    :: tol
      integer(ik) :: ifus, maxc
      logical     :: force

      if (cfg%max_cohort == 0_ik) return
      force = (cfg%max_cohort < 0_ik)
      maxc  = abs(cfg%max_cohort)
      tol   = cfg%cohort_size_tol_min

      do ifus = 1_ik, cfg%n_cohort_fusion_iter
         call fuse_pass(site, cfg, tol, force)
         if (.not. force) then
            if (max_cohort_count(site) <= maxc) exit
         end if
         tol = tol * cfg%cohort_size_tol_mult
      end do
   end subroutine new_fuse_cohorts

   !----- One fusion sweep over all patches at a fixed tolerance. --------------------------!
   subroutine fuse_pass(site, cfg, tol, force)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: tol
      logical,             intent(in)    :: force
      logical, allocatable :: alive(:)
      integer(ik) :: ip, i0, i1, recc, donc
      real(wp)    :: diff_height, lai_comb
      logical     :: did_fuse

      did_fuse = .false.
      allocate(alive(site%cohort%n))
      alive = .true.
      associate (cohort => site%cohort, patch => site%patch)
         do ip = 1_ik, patch%n
            i0 = patch%cohort_offset(ip)
            i1 = i0 + patch%cohort_count(ip) - 1_ik
            do recc = i0, i1 - 1_ik
               if (.not. alive(recc)) cycle
               do donc = recc + 1_ik, i1
                  if (.not. alive(donc)) cycle
                  if (cohort%pft(donc) /= cohort%pft(recc)) cycle
                  lai_comb = cohort%nplant(recc) * cohort%leaf_area(recc) + cohort%nplant(donc) * cohort%leaf_area(donc)
                  if (.not. force .and. lai_comb >= cfg%cohort_lai_cap) cycle
                  diff_height = abs(cohort%height(donc) - cohort%height(recc))
                  !----- Same-PFT fusion (checked above), so recc and donc share hgt_max. -------!
                  if (force .or. diff_height < cohort%p_hgt_max(recc) * tol) then
                     call fuse_2_cohorts(site, recc, donc, cfg%conservation_tol)
                     alive(donc) = .false.
                     did_fuse    = .true.
                  end if
               end do
            end do
         end do
      end associate
      if (did_fuse) then
         call cohort_compact(site%cohort, alive)
         call rebuild_csr(site)
         call sort_cohorts(site)
      end if
   end subroutine fuse_pass

   !---------------------------------------------------------------------------------------!
   ! Merge donor cohort into receptor, conserving plant number and total aboveground         !
   ! biomass (carbon). The fused per-plant AGB is the nplant-weighted mean; the diameter is   !
   ! re-derived from it (agb_to_dbh) and the rest of the geometry refreshed.                  !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_cohorts(site, recc, donc, conservation_tol)
      type(site_t), intent(inout) :: site
      integer(ik),     intent(in)    :: recc, donc
      real(wp),        intent(in)    :: conservation_tol
      real(wp) :: nr, nd, ntot, agb_tot, agb_new

      associate (cohort => site%cohort)
         nr      = cohort%nplant(recc)
         nd      = cohort%nplant(donc)
         ntot    = nr + nd
         agb_tot = nr * cohort%agb(recc) + nd * cohort%agb(donc)     ! [kgC/m2] conserved
         !----- The survivor keeps its own moving-average growth history (ring buffer + accum  !
         !      + count + growth_avg are left untouched); the donor's is discarded with it. ---!
         !----- Conserve plant number and AGB; re-derive size from carbon. ----------------!
         cohort%nplant(recc) = ntot
         cohort%agb(recc)    = agb_tot / ntot                   ! [kgC/plant]
         cohort%dbh(recc)    = agb_to_dbh(cohort%agb(recc), cohort%p_wood_density(recc), cohort%p_hgt_max(recc))
         call set_cohort_size(site%cohort, recc)              ! refresh height/BA/agb/leaf_area
         !----- Conservation guard (AGB density). -----------------------------------------!
         agb_new = cohort%nplant(recc) * cohort%agb(recc)
         if (abs(agb_new - agb_tot) > conservation_tol * max(agb_tot, tiny_num))             &
            error stop 'fuse_2_cohorts: AGB conservation violated'
      end associate
   end subroutine fuse_2_cohorts

   !---------------------------------------------------------------------------------------!
   ! Split every cohort whose LAI density exceeds the cap, iterating until none do.          !
   !---------------------------------------------------------------------------------------!
   subroutine split_cohorts(site, cfg)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      integer(ik) :: iter, i, n0, m, nsplit
      real(wp)    :: d0, eps, renorm, p_agb, agb_before, agb_after

      if (.not. cfg%enable_cohort_fission) return
      eps = cfg%split_eps
      !----- AGB ~ dbh^p_agb (uncapped); scale the two daughters' diameters so the +/-eps   !
      !      split conserves AGB exactly: 0.5*[(1+eps)^p + (1-eps)^p] * renorm^p = 1.        !
      p_agb  = agb_c2 * (2.0_wp + b2Ht)
      renorm = (0.5_wp * ((1.0_wp + eps) ** p_agb + (1.0_wp - eps) ** p_agb)) ** (-1.0_wp / p_agb)

      do iter = 1_ik, cfg%n_cohort_fusion_iter
         n0 = site%cohort%n
         nsplit = 0_ik
         do i = 1_ik, n0
            if (site%cohort%nplant(i) * site%cohort%leaf_area(i) > cfg%cohort_lai_cap) nsplit = nsplit + 1_ik
         end do
         if (nsplit == 0_ik) exit
         agb_before = sum(site%cohort%nplant(1:n0) * site%cohort%agb(1:n0))
         call cohort_ensure_capacity(site%cohort, n0 + nsplit)
         m = n0
         associate (cohort => site%cohort)
            do i = 1_ik, n0
               if (cohort%nplant(i) * cohort%leaf_area(i) <= cfg%cohort_lai_cap) cycle
               d0          = cohort%dbh(i)
               cohort%nplant(i) = 0.5_wp * cohort%nplant(i)
               !----- '+eps' half stays in slot i. ----------------------------------------!
               cohort%dbh(i) = d0 * (1.0_wp + eps) * renorm
               call set_cohort_size(cohort, i)
               !----- '-eps' half appended in slot m+1. -----------------------------------!
               m = m + 1_ik
               call copy_cohort_slot(cohort, m, i)              ! copies halved nplant + params
               cohort%dbh(m) = d0 * (1.0_wp - eps) * renorm
               call set_cohort_size(cohort, m)
            end do
            cohort%n = m
            agb_after = sum(cohort%nplant(1:m) * cohort%agb(1:m))
         end associate
         !----- The '-eps' daughters (slots n0+1..m) are NEW cohorts -> fresh global ids; the !
         !      '+eps' half kept slot i and its parent id (the continuation).                 !
         do i = n0 + 1_ik, m
            call assign_cohort_id(site, i)
         end do
         if (abs(agb_after - agb_before) > cfg%conservation_tol * max(agb_before, tiny_num))       &
            error stop 'split_cohorts: AGB conservation violated'
         call rebuild_csr(site)
         call sort_cohorts(site)
      end do
   end subroutine split_cohorts

   !---------------------------------------------------------------------------------------!
   ! Cull cohorts below the AGB-density floor or the absolute density floor.                 !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_cohorts(site, cfg)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical, allocatable :: keep(:)
      integer(ik)          :: i, n

      n = site%cohort%n
      if (n < 1_ik) return
      allocate(keep(n))
      do i = 1_ik, n
         keep(i) = (site%cohort%nplant(i) * site%cohort%agb(i) >= cfg%min_cohort_agb) .and.        &
                   (site%cohort%nplant(i) >= cfg%negligible_nplant)
      end do
      if (all(keep)) return
      call cohort_compact(site%cohort, keep)
      call rebuild_csr(site)
   end subroutine terminate_cohorts

   !=======================================================================================!
   !  Patch fusion / termination                                                            !
   !=======================================================================================!

   !---------------------------------------------------------------------------------------!
   ! Cumulative-LAI LIGHT profile for one patch (ED2 patch-fusion scheme): the light         !
   ! fraction reaching the top of each height layer, light(ihgt) = exp(-light_ext *          !
   ! LAI_at-or-above(ihgt)). Layer 1 is the ground (most overlying LAI -> darkest), layer     !
   ! n_height_layers is the canopy top. Two patches with a similar light profile (and the     !
   ! same disturbance type) are fused. Mirrors ED2 patch_pft_size_profile + new_fuse_patches. !
   !---------------------------------------------------------------------------------------!
   subroutine patch_light_profile(site, cfg, ip, light)
      type(site_t),          intent(in)  :: site
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: ip
      real(wp),            intent(out) :: light(:)         ! size cfg%n_height_layers
      real(wp)    :: cum_lai(cfg%n_height_layers)
      integer(ik) :: i, i0, i1, ihgt, nl
      nl = cfg%n_height_layers
      cum_lai = 0.0_wp
      i0 = site%patch%cohort_offset(ip)
      i1 = i0 + site%patch%cohort_count(ip) - 1_ik
      do i = i0, i1
         ihgt = min(nl, max(1_ik, count(cfg%height_edges < site%cohort%height(i)) + 1_ik))
         cum_lai(ihgt) = cum_lai(ihgt) + site%cohort%nplant(i) * site%cohort%leaf_area(i)
      end do
      !----- Integrate top-down: cum_lai(ihgt) = LAI in this layer and every layer above. --!
      do ihgt = nl - 1_ik, 1_ik, -1_ik
         cum_lai(ihgt) = cum_lai(ihgt) + cum_lai(ihgt + 1_ik)
      end do
      do ihgt = 1_ik, nl
         light(ihgt) = exp(-light_ext * cum_lai(ihgt))
      end do
   end subroutine patch_light_profile

   !---------------------------------------------------------------------------------------!
   ! Patch fusion with geometric tolerance relaxation.                                      !
   !---------------------------------------------------------------------------------------!
   subroutine new_fuse_patches(site, cfg)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp)    :: tol
      integer(ik) :: it, maxp
      logical     :: force

      if (cfg%max_patch == 0_ik) return
      force = (cfg%max_patch < 0_ik)
      maxp  = abs(cfg%max_patch)
      tol   = cfg%patch_light_tol

      do it = 1_ik, cfg%n_patch_fusion_iter
         call patch_fuse_pass(site, cfg, tol, force)
         if (.not. force .and. site%patch%n <= maxp) exit
         tol = tol * prof_relax
      end do
   end subroutine new_fuse_patches

   !----- One patch-fusion sweep at a fixed tolerance, comparing light profiles. ----------!
   subroutine patch_fuse_pass(site, cfg, tol, force)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: tol
      logical,             intent(in)    :: force
      logical,     allocatable :: alive(:)
      real(wp),    allocatable :: prof(:,:), lr(:), ld(:)
      real(wp)    :: davg, dmax, dif
      integer(ik) :: np, recp, donp, ihgt, npop, nl
      logical     :: similar, both_empty

      np = site%patch%n
      if (np < 2_ik) return
      nl = cfg%n_height_layers
      allocate(alive(np), prof(nl, np), lr(nl), ld(nl))
      alive = .true.
      do recp = 1_ik, np
         call patch_light_profile(site, cfg, recp, prof(:, recp))
      end do

      do recp = 1_ik, np - 1_ik
         if (.not. alive(recp)) cycle
         do donp = recp + 1_ik, np
            if (.not. alive(donp)) cycle
            if (site%patch%dist_type(donp) /= site%patch%dist_type(recp)) cycle
            both_empty = (site%patch%cohort_count(recp) == 0_ik .and. site%patch%cohort_count(donp) == 0_ik)
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
               call fuse_2_patches(site, recp, donp)
               call rebuild_csr(site)                 ! recp slice now holds merged cohorts
               call patch_light_profile(site, cfg, recp, prof(:, recp))
               alive(donp) = .false.
            end if
         end do
      end do

      call patch_compact(site, alive)
      call sort_cohorts(site)
   end subroutine patch_fuse_pass

   !---------------------------------------------------------------------------------------!
   ! Merge donor patch into receptor, conserving site_t-level plant number via area weights.   !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_patches(site, recp, donp)
      type(site_t), intent(inout) :: site
      integer(ik),     intent(in)    :: recp, donp
      real(wp)    :: ar, ad, anew, rawgt, dawgt
      integer(ik) :: i, i0, i1

      associate (patch => site%patch, cohort => site%cohort)
         ar = patch%area(recp) ; ad = patch%area(donp) ; anew = ar + ad
         if (anew <= tiny_num) return
         rawgt = ar / anew ; dawgt = ad / anew
         !----- Area-weighted patch scalars. ----------------------------------------------!
         patch%age(recp)            = rawgt * patch%age(recp)            + dawgt * patch%age(donp)
         patch%recruit_pool(:,recp) = rawgt * patch%recruit_pool(:,recp) + dawgt * patch%recruit_pool(:,donp)
         !----- Area-weighted fast-biophysics reservoirs (conserves the area-extensive store). !
         patch%cas(recp)    = blend_cas(rawgt,    patch%cas(recp),    dawgt, patch%cas(donp))
         patch%soil_e(recp) = blend_soil_e(rawgt, patch%soil_e(recp), dawgt, patch%soil_e(donp))
         patch%soil_w(recp) = blend_soil_w(rawgt, patch%soil_w(recp), dawgt, patch%soil_w(donp))
         !----- Rescale receptor cohort densities (slice currently holds all recp cohorts). !
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
   ! Remove patches with area below the floor; renormalize remaining areas to sum to 1.     !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_patches(site, cfg)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical, allocatable :: pkeep(:)
      integer(ik)          :: np, ip, imax
      real(wp)             :: s

      np = site%patch%n
      if (np < 1_ik) return
      allocate(pkeep(np))
      pkeep = site%patch%area(1:np) >= cfg%min_patch_area
      if (.not. any(pkeep)) then
         imax = maxloc(site%patch%area(1:np), dim=1)
         pkeep(imax) = .true.
      end if
      if (.not. all(pkeep)) call patch_compact(site, pkeep)

      !----- Renormalize areas to sum to 1. -----------------------------------------------!
      s = sum(site%patch%area(1:site%patch%n))
      if (s > tiny_num) then
         do ip = 1_ik, site%patch%n
            site%patch%area(ip) = site%patch%area(ip) / s
         end do
      end if
   end subroutine terminate_patches

   !---------------------------------------------------------------------------------------!
   ! Compact the patch list by keep-mask: drop dropped patches AND their cohorts, remap      !
   ! owner_patch to the new patch indices, and rebuild the CSR map.                          !
   !---------------------------------------------------------------------------------------!
   subroutine patch_compact(site, pkeep)
      type(site_t), intent(inout) :: site
      logical,         intent(in)    :: pkeep(:)
      integer(ik), allocatable :: newidx(:)
      logical,     allocatable :: ckeep(:)
      integer(ik) :: np, ip, k, i, nc

      np = site%patch%n
      allocate(newidx(np))
      k = 0_ik
      do ip = 1_ik, np
         if (pkeep(ip)) then
            k = k + 1_ik
            newidx(ip) = k
         else
            newidx(ip) = 0_ik
         end if
      end do

      !----- Cohorts: keep only those whose patch survives; remap owner. ------------------!
      nc = site%cohort%n
      allocate(ckeep(nc))
      do i = 1_ik, nc
         ckeep(i) = pkeep(site%cohort%owner_patch(i))
         if (ckeep(i)) site%cohort%owner_patch(i) = newidx(site%cohort%owner_patch(i))
      end do
      call cohort_compact(site%cohort, ckeep)

      !----- Patches: compact every patch array by pkeep (RHS pack -> safe temp). ---------!
      associate (patch => site%patch)
         patch%area(1:k)           = pack(patch%area(1:np),           pkeep)
         patch%age(1:k)            = pack(patch%age(1:np),            pkeep)
         patch%dist_type(1:k)      = pack(patch%dist_type(1:np),      pkeep)
         patch%global_id(1:k)      = pack(patch%global_id(1:np),      pkeep)
         patch%cas(1:k)            = pack(patch%cas(1:np),            pkeep)
         patch%soil_e(1:k)         = pack(patch%soil_e(1:np),         pkeep)
         patch%soil_w(1:k)         = pack(patch%soil_w(1:np),         pkeep)
         block
            integer(ik) :: jp
            do jp = 1_ik, site%n_pft
               patch%recruit_pool(jp, 1:k) = pack(patch%recruit_pool(jp, 1:np), pkeep)
            end do
         end block
         patch%n = k
      end associate
      call rebuild_csr(site)
   end subroutine patch_compact

end module meds_demography_fusefiss
