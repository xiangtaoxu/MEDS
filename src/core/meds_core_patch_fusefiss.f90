!==========================================================================================!
! meds_core_patch_fusefiss -- PATCH numerical-resolution control for the adaptive patch          !
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
! Depends on the cohort sibling meds_core_cohort_fusefiss (patch fusion + disturbance re-sort    !
! the cohorts they touch).                                                                       !
!==========================================================================================!
module meds_core_patch_fusefiss
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : tiny_num, almost_one
   use meds_allometry,  only : light_ext
   use meds_config,     only : meds_config_t, DIST_TREEFALL
   use meds_core_state_types, only : site_t, rebuild_csr, cohort_compact,                        &
                                      cohort_ensure_capacity, copy_cohort_slot,                    &
                                      patch_ensure_capacity, assign_cohort_id, assign_patch_id
   use meds_core_cohort_fusefiss, only : sort_cohorts
   use meds_column_state_types, only : blend_cas, blend_soil_w, blend_soil_e, blend_snow, snow_column_t, &
                                      blend_soil_carbon, necromass_to_litter, blend_xi_accum
   implicit none
   private

   public :: sort_patches, apply_patch_disturbance
   public :: new_fuse_patches, fuse_2_patches, terminate_patches, patch_light_profile

   real(wp), parameter :: prof_relax = 1.5_wp     !< per-iteration patch-fusion tolerance growth

contains

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
         patch%snow(1:np)           = patch%snow(pperm(1:np))
         patch%soil_carbon(1:np)    = patch%soil_carbon(pperm(1:np))
         patch%xi_accum(1:np)       = patch%xi_accum(pperm(1:np))
         patch%shed_water_rate(1:np) = patch%shed_water_rate(pperm(1:np))
         patch%adapt_dt_last(1:np)   = patch%adapt_dt_last(pperm(1:np))
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
         if (anew <= tiny_num) then                     ! degenerate: still hand donor cohorts to recp
            i0 = patch%cohort_offset(donp) ; i1 = i0 + patch%cohort_count(donp) - 1_ik
            do i = i0, i1
               cohort%owner_patch(i) = recp
            end do
            return
         end if
         rawgt = ar / anew ; dawgt = ad / anew
         !----- Area-weighted patch scalars. ----------------------------------------------!
         patch%age(recp)            = rawgt * patch%age(recp)            + dawgt * patch%age(donp)
         patch%recruit_pool(:,recp) = rawgt * patch%recruit_pool(:,recp) + dawgt * patch%recruit_pool(:,donp)
         !----- Area-weighted fast-biophysics reservoirs (conserves the area-extensive store). !
         patch%cas(recp)    = blend_cas(rawgt,    patch%cas(recp),    dawgt, patch%cas(donp))
         patch%soil_e(recp) = blend_soil_e(rawgt, patch%soil_e(recp), dawgt, patch%soil_e(donp))
         patch%soil_w(recp) = blend_soil_w(rawgt, patch%soil_w(recp), dawgt, patch%soil_w(donp))
         patch%snow(recp)   = blend_snow(rawgt,   patch%snow(recp),   dawgt, patch%snow(donp))  ! temp/fliq re-diagnosed in the fast loop
         !----- Area-weighted slow soil-carbon reservoir (conserves site-wide soil carbon). -----!
         patch%soil_carbon(recp) = blend_soil_carbon(rawgt, patch%soil_carbon(recp), dawgt, patch%soil_carbon(donp))
         patch%xi_accum(recp)    = blend_xi_accum(rawgt, patch%xi_accum(recp), dawgt, patch%xi_accum(donp))
         !----- shed_water_rate is a per-area RATE (like age): area-weighted, not nplant-weighted. -----!
         patch%shed_water_rate(recp) = rawgt*patch%shed_water_rate(recp) + dawgt*patch%shed_water_rate(donp)
         !----- adapt_dt_last is a controller SEED, not a conserved amount -- any value is valid and   !
         !      the controller re-adapts within a step. Area-weight it like its neighbours purely so    !
         !      the result is DETERMINISTIC and order-independent (issue #106's whole point). ---------!
         patch%adapt_dt_last(recp)   = rawgt*patch%adapt_dt_last(recp)   + dawgt*patch%adapt_dt_last(donp)
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
         patch%snow(1:k)           = pack(patch%snow(1:np),           pkeep)
         patch%soil_carbon(1:k)    = pack(patch%soil_carbon(1:np),    pkeep)
         patch%xi_accum(1:k)       = pack(patch%xi_accum(1:np),       pkeep)
         patch%shed_water_rate(1:k) = pack(patch%shed_water_rate(1:np), pkeep)
         patch%adapt_dt_last(1:k)   = pack(patch%adapt_dt_last(1:np),   pkeep)
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

   !---------------------------------------------------------------------------------------!
   ! Treefall patch disturbance (ED2 analogue). A fraction f = 1 - exp(-rate*dt) of EVERY     !
   ! patch's area is disturbed and aggregated into ONE new age-0 gap patch (DIST_TREEFALL).   !
   ! Survivorship follows ED2 treefall: tall canopy cohorts (height >= survive_height) die in !
   ! the gap; the short understory (height < threshold) survives into it. Each donor keeps    !
   ! area (1-f)*area with ALL its cohorts intact. Site area is conserved; plant number is     !
   ! conserved for survivors and reduced by f for the killed canopy -- that loss IS the        !
   ! disturbance. The gap is later consolidated by patch and cohort fusion.                   !
   !---------------------------------------------------------------------------------------!
   subroutine apply_patch_disturbance(site, cfg, dt_yr)
      type(site_t),          intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: dt_yr
      integer(ik) :: np0, newp, d, i, i0, i1, m, m0, nsurv, pf
      real(wp)    :: frac, new_area, atot, wd
      real(wp)    :: lost_density, lab_g, lab_s, str_g, str_s, lig_g, lig_s

      np0 = site%patch%n
      if (np0 < 1_ik .or. cfg%patch_disturbance_rate <= 0.0_wp) return

      !----- Disturbed area fraction this step, and the aggregate new-patch area. ----------!
      frac = 1.0_wp - exp(-cfg%patch_disturbance_rate * dt_yr)
      if (frac <= tiny_num) return
      new_area = frac * sum(site%patch%area(1:np0))
      if (new_area <= tiny_num) return

      !----- Count understorey survivors (all current cohorts live in donor patches). ------!
      nsurv = 0_ik
      do i = 1_ik, site%cohort%n
         if (site%cohort%height(i) < cfg%disturbance_survive_height) nsurv = nsurv + 1_ik
      end do

      call patch_ensure_capacity(site%patch, np0 + 1_ik, site%n_pft)
      call cohort_ensure_capacity(site%cohort, site%cohort%n + nsurv)
      newp = np0 + 1_ik
      m0   = site%cohort%n                       ! cohorts before the survivors are moved in

      associate (cohort => site%cohort, patch => site%patch)
         !----- New age-0 gap patch (treefall). --------------------------------------------!
         patch%area(newp)           = new_area
         patch%age(newp)            = 0.0_wp
         patch%dist_type(newp)      = DIST_TREEFALL
         patch%recruit_pool(:,newp) = 0.0_wp
         !----- Fresh gap, like age: it has no cohorts of its own that contributed to today's shed  !
         !      rate (unlike soil_carbon/xi_accum, which are genuine inherited material). ----------!
         patch%shed_water_rate(newp) = 0.0_wp
         patch%adapt_dt_last(newp)   = 0.0_wp   ! a fresh gap cold-starts the controller
         patch%n = newp

         !----- Seed the gap's fast reservoirs = area-weighted donor mean (the soil column and !
         !      canopy air are INHERITED from the disturbed area, not created bare). ----------!
         atot = sum(patch%area(1:np0))                       ! donors still hold pre-disturbance area here
         if (atot > tiny_num) then
            patch%cas(newp)    = blend_cas(   patch%area(1)/atot, patch%cas(1),    0.0_wp, patch%cas(1))
            patch%soil_e(newp) = blend_soil_e(patch%area(1)/atot, patch%soil_e(1), 0.0_wp, patch%soil_e(1))
            patch%soil_w(newp) = blend_soil_w(patch%area(1)/atot, patch%soil_w(1), 0.0_wp, patch%soil_w(1))
            patch%snow(newp)   = blend_snow(  patch%area(1)/atot, patch%snow(1),   0.0_wp, patch%snow(1))
            patch%soil_carbon(newp) = blend_soil_carbon(patch%area(1)/atot, patch%soil_carbon(1), &
                                                        0.0_wp, patch%soil_carbon(1))
            patch%xi_accum(newp) = blend_xi_accum(patch%area(1)/atot, patch%xi_accum(1), 0.0_wp, patch%xi_accum(1))
            do d = 2_ik, np0
               wd = patch%area(d) / atot
               patch%cas(newp)    = blend_cas(   1.0_wp, patch%cas(newp),    wd, patch%cas(d))
               patch%soil_e(newp) = blend_soil_e(1.0_wp, patch%soil_e(newp), wd, patch%soil_e(d))
               patch%soil_w(newp) = blend_soil_w(1.0_wp, patch%soil_w(newp), wd, patch%soil_w(d))
               patch%snow(newp)   = blend_snow(  1.0_wp, patch%snow(newp),   wd, patch%snow(d))  ! conserve snow into the gap
               patch%soil_carbon(newp) = blend_soil_carbon(1.0_wp, patch%soil_carbon(newp), wd, patch%soil_carbon(d))
               patch%xi_accum(newp)    = blend_xi_accum(1.0_wp, patch%xi_accum(newp), wd, patch%xi_accum(d))
            end do
         end if

         !----- Move understorey survivors into the gap at area-weighted density; the killed       !
         !      canopy's carbon becomes litter into the SAME gap patch (B1, MEDS_SLOW_DYNAMICS_     !
         !      DESIGN.md Part II; OPT-IN [soil_carbon].soil_carbon_on -- default .false. keeps       !
         !      this bit-identical) -- the density it would have carried into the gap had it            !
         !      survived (the same conversion factor line 394 uses for survivors), times its per-      !
         !      plant carbon pools. Added directly onto soil_carbon(newp) since this module cannot     !
         !      link biogeochemistry (necromass_to_litter is DAG-safe: plain scalars). -----------------!
         m = cohort%n
         do d = 1_ik, np0
            i0 = patch%cohort_offset(d) ; i1 = i0 + patch%cohort_count(d) - 1_ik
            do i = i0, i1
               if (cohort%height(i) >= cfg%disturbance_survive_height) then   ! canopy dies in gap
                  if (cfg%soil_carbon_on) then
                     lost_density = cohort%nplant(i) * (frac * patch%area(d) / new_area)
                     pf = cohort%pft(i)
                     call necromass_to_litter(lost_density * cohort%leaf_carbon(i),                  &
                              lost_density * cohort%fineroot_carbon(i),                               &
                              lost_density * cohort%wood_carbon(i),                                   &
                              lost_density * cohort%nonstructural_carbon(i),                          &
                              cfg%pft%f_labile_leaf(pf), cfg%pft%f_labile_stem(pf),                   &
                              cfg%pft%aboveground_frac(pf), cfg%pft%struct_lignin_frac(pf),           &
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
               call copy_cohort_slot(cohort, m, i)
               cohort%nplant(m)      = cohort%nplant(i) * (frac * patch%area(d) / new_area)
               cohort%owner_patch(m) = newp
            end do
         end do
         cohort%n = m

         !----- Donors keep the undisturbed remainder (all their cohorts intact). ----------!
         do d = 1_ik, np0
            patch%area(d) = (1.0_wp - frac) * patch%area(d)
         end do
      end associate

      !----- Stamp the new gap patch and the moved-in survivor cohorts with fresh global ids !
      !      (the gap fragments are new entities; their donor cohorts keep their own ids).    !
      call assign_patch_id(site, newp)
      do i = m0 + 1_ik, site%cohort%n
         call assign_cohort_id(site, i)
      end do

      call rebuild_csr(site)
      call sort_cohorts(site)
   end subroutine apply_patch_disturbance

end module meds_core_patch_fusefiss
