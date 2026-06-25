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
   use meds_allometry, only : b1Ht, b2Ht, height_max, agb_c1, agb_c2, lai_b1, lai_b2
   use meds_config,    only : meds_config_t, DIST_TREEFALL
   use meds_demography_types,     only : site_t, patch_ensure_capacity, cohort_ensure_capacity,  &
                                         copy_cohort_slot, rebuild_csr
   use meds_demography_structure, only : sort_cohorts
   implicit none
   private

   public :: growth_step, mortality_step, apply_patch_disturbance

contains

   !---------------------------------------------------------------------------------------!
   ! Daily growth: advance DBH by the supplied per-cohort growth rate (capped at dbh_critical)   !
   ! and re-derive the cached geometry -- height, basal area, AGB (carbon) and per-stem leaf  !
   ! area -- from the pan-tropical allometry. The allometry is INLINED (not a procedure call) !
   ! so the loop offloads cleanly; it must mirror meds_allometry / set_cohort_size exactly.   !
   ! Arrays are indexed in the current cohort order (1:n).                                    !
   !---------------------------------------------------------------------------------------!
   subroutine growth_step(n, dbh, height, basal_area, agb, leaf_area, dbh_critical, wood_density,    &
                          growth, dt_yr)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: dbh(:), height(:), basal_area(:), agb(:), leaf_area(:)
      real(wp),    intent(in)    :: dbh_critical(:), wood_density(:)
      real(wp),    intent(in)    :: growth(:)
      real(wp),    intent(in)    :: dt_yr
      integer(ik) :: i
      real(wp)    :: size_var

      !$omp target teams distribute parallel do simd                                        &
      !$omp&        map(to: dbh_critical, wood_density, growth)                                  &
      !$omp&        map(tofrom: dbh, height, basal_area, agb, leaf_area) private(size_var)
      do i = 1_ik, n
         dbh(i)        = min(dbh(i) + growth(i) * dt_yr, dbh_critical(i))
         height(i)     = min(exp(b1Ht + b2Ht * log(dbh(i))), height_max)
         basal_area(i) = pio4 * dbh(i) * dbh(i)
         size_var      = dbh(i) * dbh(i) * height(i)
         agb(i)        = agb_c1 * wood_density(i) ** agb_c2 * size_var ** agb_c2
         leaf_area(i)  = lai_b1 * size_var ** lai_b2
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
      integer(ik) :: np0, newp, d, i, i0, i1, m, nsurv
      real(wp)    :: frac, new_area, wsum, tavg, tmin, w

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

      associate (cohort => site%cohort, patch => site%patch)
         !----- New age-0 gap patch: area-weighted temperatures over the donors. -----------!
         wsum = 0.0_wp ; tavg = 0.0_wp ; tmin = huge(1.0_wp)
         do d = 1_ik, np0
            w    = frac * patch%area(d)
            wsum = wsum + w
            tavg = tavg + w * patch%avg_daily_temp(d)
            tmin = min(tmin, patch%min_month_temp(d))
         end do
         tavg = tavg / max(wsum, tiny_num)
         patch%area(newp)           = new_area
         patch%age(newp)            = 0.0_wp
         patch%dist_type(newp)      = DIST_TREEFALL
         patch%avg_daily_temp(newp) = tavg
         patch%min_month_temp(newp) = tmin
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

      call rebuild_csr(site)
      call sort_cohorts(site)
   end subroutine apply_patch_disturbance

end module meds_demography_dynamics
