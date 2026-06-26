!==========================================================================================!
! meds_demography_dynamics -- the per-step demographic PROCESSES that change cohort/patch     !
! state: the dominant per-cohort array kernels (growth, mortality) and treefall patch          !
! disturbance. (Renamed from meds_demography_kernels and merged with meds_disturbance.)        !
!                                                                                          !
! growth_step / mortality_step are PURE ARRAY kernels: branch-light arithmetic over plain     !
! 1-D arrays (no derived types), each wrapped in an explicit OpenMP `target` region with       !
! `map` clauses, so the data is copied to/from the device per call and the host keeps ALL      !
! state in normal memory. (This deliberately avoids `-stdpar=gpu`'s global CUDA-managed        !
! allocator, whose deep-copy/finalize of the allocatable-component `site_t` type double-frees.)  !
! One directive serves three back ends via the build flag: nvfortran `-mp=gpu` -> GPU, `-mp`   !
! -> CPU threads, no flag -> the `!$omp` lines are comments and the loop runs serially.        !
!                                                                                          !
! apply_patch_disturbance is a HOST structural operation (it adds a patch); it shares the      !
! module with the kernels because all three are the per-step state-changing processes.         !
!==========================================================================================!
module meds_demography_dynamics
   use meds_kinds,     only : wp, ik
   use meds_constants, only : pio4, tiny_num, lnexp_min, lnexp_max
   use meds_allometry, only : b1Ht, b2Ht, height_max, agb_c1, agb_c2, lai_b1, lai_b2, height_to_dbh
   use meds_config,    only : meds_config_t, DIST_TREEFALL
   use meds_demography_types,     only : site_t, patch_ensure_capacity, cohort_ensure_capacity,  &
                                         copy_cohort_slot, rebuild_csr, assign_cohort_id,         &
                                         assign_patch_id, set_cohort_size, GROWTH_AVG_UNSET
   use meds_demography_structure, only : sort_cohorts
   implicit none
   private

   public :: growth_step, mortality_step, apply_patch_disturbance, apply_recruitment

contains

   !---------------------------------------------------------------------------------------!
   ! Daily growth: advance DBH by the supplied per-cohort growth rate (capped at dbh_critical)   !
   ! and re-derive the cached geometry -- height, basal area, AGB (carbon) and per-stem leaf  !
   ! area -- from the pan-tropical allometry. The allometry is INLINED (not a procedure call) !
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
   subroutine growth_step(n, dbh, height, basal_area, agb, leaf_area, growth_avg, growth_accum,     &
                          growth_count, growth_hist, dbh_critical, wood_density, growth, dt_yr,      &
                          n_window, hist_pos)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: dbh(:), height(:), basal_area(:), agb(:), leaf_area(:)
      real(wp),    intent(inout) :: growth_avg(:), growth_accum(:), growth_hist(:,:)
      integer(ik), intent(inout) :: growth_count(:)
      real(wp),    intent(in)    :: dbh_critical(:), wood_density(:)
      real(wp),    intent(in)    :: growth(:)
      real(wp),    intent(in)    :: dt_yr
      integer(ik), intent(in)    :: n_window, hist_pos
      integer(ik) :: i
      real(wp)    :: size_var

      !$omp target teams distribute parallel do simd                                        &
      !$omp&        map(to: dbh_critical, wood_density, growth)                                  &
      !$omp&        map(tofrom: dbh, height, basal_area, agb, leaf_area, growth_avg, growth_accum, &
      !$omp&            growth_count, growth_hist) private(size_var)
      do i = 1_ik, n
         dbh(i)        = min(dbh(i) + growth(i) * dt_yr, dbh_critical(i))
         height(i)     = min(exp(b1Ht + b2Ht * log(dbh(i))), height_max)
         basal_area(i) = pio4 * dbh(i) * dbh(i)
         size_var      = dbh(i) * dbh(i) * height(i)
         agb(i)        = agb_c1 * wood_density(i) ** agb_c2 * size_var ** agb_c2
         leaf_area(i)  = lai_b1 * size_var ** lai_b2
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
      real(wp)    :: frac, new_area

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
   ! from the physiology layer (meds_phenomenological_vital_rates) via the update_demography seam.!
   !---------------------------------------------------------------------------------------!
   subroutine apply_recruitment(site, cfg, recruitment)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: recruitment(:,:)  !< [plant/m2/month] (pft, patch)
      integer(ik) :: ip, pf, np, m, nspawn, n_before, k
      real(wp)    :: recruit_dbh

      np = site%patch%n
      if (np < 1_ik) return

      !----- All PFTs recruit at the smallest tracked size -> the same diameter. ----------!
      recruit_dbh = height_to_dbh(cfg%pft%min_cohort_height)

      !----- Accumulate the supplied monthly recruit density into the carry-forward pool. -!
      do ip = 1_ik, np
         do pf = 1_ik, site%n_pft
            site%patch%recruit_pool(pf, ip) = site%patch%recruit_pool(pf, ip) + recruitment(pf, ip)
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
