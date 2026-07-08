!==========================================================================================!
! meds_demography_dynamics -- the per-step demographic EVENTS: what the ECOSYSTEM does to      !
! itself each step -- grow (growth_step | apply_carbon_npp), die (mortality_step), be born      !
! (apply_recruitment), be knocked down (apply_patch_disturbance). Driven by ECOLOGY; they       !
! change the state and the counts. The numerical housekeeping that keeps the adaptive           !
! discretization bounded (sort/fuse/split/cull) is the counterpart module meds_demography_fusefiss.!
!                                                                                          !
! growth_step / mortality_step are PURE ARRAY kernels: branch-light arithmetic over plain     !
! 1-D arrays (no derived types), each wrapped in an explicit OpenMP `target` region with       !
! `map` clauses, so the data is copied to/from the device per call and the host keeps ALL      !
! state in normal memory. (This deliberately avoids `-stdpar=gpu`'s global CUDA-managed        !
! allocator, whose deep-copy/finalize of the allocatable-component `site_t` type double-frees.)  !
! One directive serves three back ends via the build flag: nvfortran `-mp=gpu` -> GPU, `-mp`   !
! -> CPU threads, no flag -> the `!$omp` lines are comments and the loop runs serially.        !
!                                                                                          !
! apply_carbon_npp / apply_recruitment / apply_patch_disturbance are HOST operations (they call !
! the allometric inverse, or change the cohort/patch count); they sit with the offload kernels  !
! because all are the per-step ecological processes.                                            !
!==========================================================================================!
module meds_demography_dynamics
   use meds_kinds,     only : wp, ik
   use meds_constants, only : pio4, tiny_num, lnexp_min, lnexp_max, mon_per_yr
   use meds_allometry, only : height_to_dbh
   use meds_config,    only : meds_config_t, DIST_TREEFALL
   use meds_demography_types,     only : site_t, cohort_block, carbon_flux_block,                 &
                                         patch_ensure_capacity, cohort_ensure_capacity,          &
                                         copy_cohort_slot, rebuild_csr, assign_cohort_id,         &
                                         assign_patch_id, set_cohort_size,                        &
                                         set_cohort_size_from_carbon, GROWTH_AVG_UNSET
   use meds_demography_fusefiss, only : sort_cohorts
   use meds_column_state_types,  only : blend_cas, blend_soil_w, blend_soil_e, LEAF_TEMP_INIT, PSI_INIT
   implicit none
   private

   public :: growth_step, mortality_step, apply_carbon_npp, apply_patch_disturbance, apply_recruitment

contains

   !---------------------------------------------------------------------------------------!
   ! Daily growth: advance DBH by the supplied per-cohort growth rate (capped at dbh_critical)   !
   ! and re-derive the cached geometry -- height, basal area, AGB, leaf area AND the four        !
   ! on-allometry carbon pools -- from the allometry. The allometry is INLINED (not a call)      !
   ! so the loop offloads cleanly; it must mirror meds_allometry / set_cohort_size exactly.   !
   ! It also advances the SLIDING SIMPLE MOVING AVERAGE of the requested growth RATE that      !
   ! predicts growth-dependent mortality. Each cohort keeps a ring buffer (growth_hist column)  !
   ! of the last n_window samples (n_window = memory-window / step); all cohorts write the same  !
   ! site-global slot hist_pos this step. While the buffer is filling (growth_count < n_window)  !
   ! the new sample is added; once full, it replaces (evicts) the n_window-old sample at         !
   ! hist_pos. growth_avg = running sum / count is the windowed mean (GROWTH_AVG_UNSET while     !
   ! count is 0). Each cohort writes consecutive slots, so the slot at hist_pos when it is full  !
   ! is exactly its oldest sample -- correct eviction regardless of when it was born. Order 1:n. !
   !---------------------------------------------------------------------------------------!
   subroutine growth_step(n, dbh, height, basal_area, agb, leaf_area,                             &
                          leaf_carbon, fineroot_carbon, wood_carbon, nonstructural_carbon,         &
                          growth_avg, growth_accum,                                                &
                          growth_count, growth_hist, dbh_critical, wood_density, hgt_max,           &
                          sla, aboveground_frac, root_to_leaf_ratio, storage_cushion,              &
                          growth, dt_yr, n_window, hist_pos, b1Ht, b2Ht, agb_c1, agb_c2, lai_b1, lai_b2)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: dbh(:), height(:), basal_area(:), agb(:), leaf_area(:)
      real(wp),    intent(inout) :: leaf_carbon(:), fineroot_carbon(:), wood_carbon(:), nonstructural_carbon(:)
      real(wp),    intent(inout) :: growth_avg(:), growth_accum(:), growth_hist(:,:)
      integer(ik), intent(inout) :: growth_count(:)
      real(wp),    intent(in)    :: dbh_critical(:), wood_density(:), hgt_max(:)
      real(wp),    intent(in)    :: sla(:), aboveground_frac(:), root_to_leaf_ratio(:), storage_cushion(:)
      real(wp),    intent(in)    :: growth(:)
      real(wp),    intent(in)    :: dt_yr
      integer(ik), intent(in)    :: n_window, hist_pos
      !----- Global allometry coefficients as firstprivate scalars (the kernel can't read host  !
      !       module state on the device); the per-PFT height cap arrives as the hgt_max array.  !
      real(wp),    intent(in)    :: b1Ht, b2Ht, agb_c1, agb_c2, lai_b1, lai_b2
      integer(ik) :: i
      real(wp)    :: size_var

      !----- Read-modify-write state stays tofrom (whole array, tail preserved); the derived        !
      !       fields are WRITE-ONLY here (device recomputes them), so map(from: 1:n) skips the         !
      !       host->device copy AND avoids clobbering the dead tail slots n+1:cap on the GPU. --------!
      !$omp target teams distribute parallel do simd                                        &
      !$omp&        map(to: dbh_critical, wood_density, hgt_max, growth,                          &
      !$omp&               sla, aboveground_frac, root_to_leaf_ratio, storage_cushion)            &
      !$omp&        map(tofrom: dbh, growth_accum, growth_count, growth_hist)                     &
      !$omp&        map(from: height(1:n), basal_area(1:n), agb(1:n), leaf_area(1:n),             &
      !$omp&                  leaf_carbon(1:n), fineroot_carbon(1:n), wood_carbon(1:n),           &
      !$omp&                  nonstructural_carbon(1:n), growth_avg(1:n)) private(size_var)
      do i = 1_ik, n
         dbh(i)        = min(dbh(i) + growth(i) * dt_yr, dbh_critical(i))
         height(i)     = min(exp(b1Ht + b2Ht * log(dbh(i))), hgt_max(i))
         basal_area(i) = pio4 * dbh(i) * dbh(i)
         size_var      = dbh(i) * dbh(i) * height(i)
         agb(i)        = agb_c1 * wood_density(i) ** agb_c2 * size_var ** agb_c2
         leaf_area(i)  = lai_b1 * size_var ** lai_b2
         !----- On-allometry carbon pools (mirror set_cohort_size; reuse agb & leaf_area). --------!
         leaf_carbon(i)          = leaf_area(i) / max(sla(i), tiny_num)
         wood_carbon(i)          = agb(i) / max(aboveground_frac(i), tiny_num)
         fineroot_carbon(i)      = root_to_leaf_ratio(i) * leaf_carbon(i)
         nonstructural_carbon(i) = storage_cushion(i) * leaf_carbon(i)
         !----- Sliding moving average: add (filling) or evict-and-replace (full). ------------!
         if (growth_count(i) < n_window) then
            growth_accum(i) = growth_accum(i) + growth(i)
            growth_count(i) = growth_count(i) + 1_ik
         else
            growth_accum(i) = growth_accum(i) + growth(i) - growth_hist(hist_pos, i)
         end if
         growth_hist(hist_pos, i) = growth(i)
         growth_avg(i) = growth_accum(i) / real(growth_count(i), wp)
      end do
   end subroutine growth_step

   !---------------------------------------------------------------------------------------!
   ! Per-step mortality: apply the supplied per-cohort total mortality rate directly to       !
   ! nplant for one step of length dt_yr, floored (downward only) so density never rises and  !
   ! never drops below `negligible`.                                                          !
   !---------------------------------------------------------------------------------------!
   subroutine mortality_step(n, nplant, mortality, dt_yr, negligible)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: nplant(:)
      real(wp),    intent(in)    :: mortality(:)
      real(wp),    intent(in)    :: dt_yr, negligible
      integer(ik) :: i
      real(wp)    :: dln, floor_dln

      !$omp target teams distribute parallel do simd                                        &
      !$omp&        map(to: mortality) map(tofrom: nplant) private(dln, floor_dln)
      do i = 1_ik, n
         dln         = -mortality(i) * dt_yr
         floor_dln   = min(0.0_wp, log(max(negligible, tiny_num) / max(nplant(i), tiny_num)))
         dln         = max(dln, floor_dln)
         nplant(i)   = nplant(i) * exp(min(max(dln, lnexp_min), lnexp_max))
      end do
   end subroutine mortality_step

   !---------------------------------------------------------------------------------------!
   ! CARBON-mode growth application (host; the carbon analogue of growth_step). Adds this        !
   ! step's per-cohort NPP to the four prognostic carbon pools, re-derives dbh from the (now     !
   ! anchor) wood_carbon and the cached geometry via set_cohort_size_from_carbon (the FLIP), and !
   ! feeds the resulting dbh-increment RATE into the SAME moving-average machinery that drives   !
   ! Camac mortality. Carbon starvation is physical: npp_wood ~ 0 -> dbh flat -> growth_avg      !
   ! falls -> low-growth hazard rises. Not offloaded (it calls the allometric wood_to_dbh).      !
   !---------------------------------------------------------------------------------------!
   subroutine apply_carbon_npp(cohort, npp, dt_yr, n_window, hist_pos)
      type(cohort_block),      intent(inout) :: cohort
      type(carbon_flux_block), intent(in)    :: npp
      real(wp),                intent(in)    :: dt_yr
      integer(ik),             intent(in)    :: n_window, hist_pos
      integer(ik) :: i
      real(wp)    :: dbh_old, dbh_rate

      do i = 1_ik, cohort%n
         !----- Add this step's NPP to the prognostic pools (never let a pool go negative). ---!
         cohort%leaf_carbon(i)          = max(cohort%leaf_carbon(i)          + npp%leaf(i),          0.0_wp)
         cohort%fineroot_carbon(i)      = max(cohort%fineroot_carbon(i)      + npp%fineroot(i),      0.0_wp)
         cohort%wood_carbon(i)          = max(cohort%wood_carbon(i)          + npp%wood(i),          0.0_wp)
         cohort%nonstructural_carbon(i) = max(cohort%nonstructural_carbon(i) + npp%nonstructural(i), 0.0_wp)
         !----- FLIP the geometry: dbh from wood_carbon, leaf_area from leaf_carbon. ----------!
         dbh_old = cohort%dbh(i)
         call set_cohort_size_from_carbon(cohort, i)
         !----- Feed the derived dbh-increment RATE [cm/yr] into the moving-average ring buffer !
         !      (identical eviction logic to growth_step) so Camac mortality sees carbon growth. !
         dbh_rate = (cohort%dbh(i) - dbh_old) / dt_yr
         if (cohort%growth_count(i) < n_window) then
            cohort%growth_accum(i) = cohort%growth_accum(i) + dbh_rate
            cohort%growth_count(i) = cohort%growth_count(i) + 1_ik
         else
            cohort%growth_accum(i) = cohort%growth_accum(i) + dbh_rate - cohort%growth_hist(hist_pos, i)
         end if
         cohort%growth_hist(hist_pos, i) = dbh_rate
         cohort%growth_avg(i) = cohort%growth_accum(i) / real(cohort%growth_count(i), wp)
      end do
   end subroutine apply_carbon_npp

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
      integer(ik) :: np0, newp, d, i, i0, i1, m, m0, nsurv
      real(wp)    :: frac, new_area, atot, wd

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
         patch%n = newp

         !----- Seed the gap's fast reservoirs = area-weighted donor mean (the soil column and !
         !      canopy air are INHERITED from the disturbed area, not created bare). ----------!
         atot = sum(patch%area(1:np0))                       ! donors still hold pre-disturbance area here
         if (atot > tiny_num) then
            patch%cas(newp)    = blend_cas(   patch%area(1)/atot, patch%cas(1),    0.0_wp, patch%cas(1))
            patch%soil_e(newp) = blend_soil_e(patch%area(1)/atot, patch%soil_e(1), 0.0_wp, patch%soil_e(1))
            patch%soil_w(newp) = blend_soil_w(patch%area(1)/atot, patch%soil_w(1), 0.0_wp, patch%soil_w(1))
            do d = 2_ik, np0
               wd = patch%area(d) / atot
               patch%cas(newp)    = blend_cas(   1.0_wp, patch%cas(newp),    wd, patch%cas(d))
               patch%soil_e(newp) = blend_soil_e(1.0_wp, patch%soil_e(newp), wd, patch%soil_e(d))
               patch%soil_w(newp) = blend_soil_w(1.0_wp, patch%soil_w(newp), wd, patch%soil_w(d))
            end do
         end if

         !----- Move understorey survivors into the gap at area-weighted density. ----------!
         m = cohort%n
         do d = 1_ik, np0
            i0 = patch%cohort_offset(d) ; i1 = i0 + patch%cohort_count(d) - 1_ik
            do i = i0, i1
               if (cohort%height(i) >= cfg%disturbance_survive_height) cycle   ! canopy dies in gap
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

   !---------------------------------------------------------------------------------------!
   ! Apply the supplied per-(PFT,patch) recruitment rate: accumulate it into a carry-forward !
   ! pool and, when a pool reaches `min_recruit_size`, spawn ONE new cohort at the shared      !
   ! minimum cohort height (the smallest tracked size; the pool is reset, otherwise it carries  !
   ! over so rare recruiters still establish eventually). A HOST structural process -- it changes !
   ! the cohort count -- so it lives next to apply_patch_disturbance; the recruitment RATE comes  !
   ! from the vegetation-dynamics driver (meds_vegetation_dynamics) via the update_demography seam.!
   !---------------------------------------------------------------------------------------!
   subroutine apply_recruitment(site, cfg, recruitment)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: recruitment(:,:)  !< [plant/m2/yr] (pft, patch)
      integer(ik) :: ip, pf, np, m, nspawn, n_before, k
      real(wp)    :: recruit_dbh

      np = site%patch%n
      if (np < 1_ik) return

      !----- All PFTs recruit at the smallest tracked size -> the same diameter. ----------!
      recruit_dbh = height_to_dbh(cfg%pft%min_cohort_height)

      !----- Accumulate one month of the supplied per-YEAR recruit density into the carry- !
      !      forward pool (this routine runs once a month, so add rate / mon_per_yr). -------!
      do ip = 1_ik, np
         do pf = 1_ik, site%n_pft
            site%patch%recruit_pool(pf, ip) = site%patch%recruit_pool(pf, ip)                &
                                              + recruitment(pf, ip) / mon_per_yr
         end do
      end do

      !----- Count pools that have reached the spawn threshold. ---------------------------!
      nspawn = 0_ik
      do ip = 1_ik, np
         do pf = 1_ik, site%n_pft
            if (site%patch%recruit_pool(pf, ip) >= cfg%min_recruit_size) nspawn = nspawn + 1_ik
         end do
      end do
      if (nspawn == 0_ik) return

      call cohort_ensure_capacity(site%cohort, site%cohort%n + nspawn)
      n_before = site%cohort%n
      m = site%cohort%n
      associate (cohort => site%cohort, patch => site%patch, pft => cfg%pft)
         do ip = 1_ik, np
            do pf = 1_ik, site%n_pft
               if (patch%recruit_pool(pf, ip) < cfg%min_recruit_size) cycle
               m = m + 1_ik
               cohort%pft(m)            = pf
               cohort%owner_patch(m)    = ip
               cohort%nplant(m)         = patch%recruit_pool(pf, ip)
               cohort%dbh(m)            = recruit_dbh
               cohort%growth_avg(m)     = GROWTH_AVG_UNSET   ! set on its first growth step
               cohort%growth_accum(m)   = 0.0_wp
               cohort%growth_count(m)   = 0_ik
               cohort%p_dbh_critical(m)     = pft%dbh_critical(pf)
               cohort%p_wood_density(m) = pft%wood_density(pf)
               cohort%p_hgt_max(m)      = pft%hgt_max(pf)
               cohort%p_sla(m)                = pft%sla(pf)
               cohort%p_aboveground_frac(m)   = pft%aboveground_frac(pf)
               cohort%p_root_to_leaf_ratio(m) = pft%root_to_leaf_ratio(pf)
               cohort%p_storage_cushion(m)    = pft%storage_cushion(pf)
               cohort%leaf_temp(m)      = LEAF_TEMP_INIT   ! fresh fast state (slot may be a reused, stale cull)
               cohort%psi(:,m)          = PSI_INIT
               call set_cohort_size(cohort, m)         ! height/basal_area/agb/leaf_area from dbh
               patch%recruit_pool(pf, ip) = 0.0_wp
            end do
         end do
         cohort%n = m
      end associate

      !----- Stamp each freshly spawned cohort with a persistent global id. ----------------!
      do k = n_before + 1_ik, site%cohort%n
         call assign_cohort_id(site, k)
      end do

      call rebuild_csr(site)
      call sort_cohorts(site)
   end subroutine apply_recruitment

end module meds_demography_dynamics
