! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_site_state_vars -- the site: one cohort block and one patch block, and the          !
! routines that need both (the CSR map, the persistent ids, the step diagnostic reset).    !
!                                                                                          !
! A polygon holds exactly one site_t, so this is the polygon's ecological state too;       !
! meds_polygon_t (src/main) wraps it with the run's forcing, budgets and ledgers.          !
!==========================================================================================!
module meds_site_state_vars
   use meds_diag_vars,          only : polygon_diag_block, cohort_diag_reset, patch_diag_reset,          &
                                       polygon_diag_reset
   use meds_kinds,              only : wp, ik
   use meds_allometry,          only : allometry_t
   use meds_pft_params,         only : pft_table_t
   use meds_cohort_state_vars,  only : cohort_block, cohort_deriv_block, alloc_cohort_block, free_cohort_block, reorder_cohort_block
   use meds_patch_state_vars,   only : patch_block, alloc_patch_block, free_patch_block
   implicit none
   private

   public :: site_t, site_alloc, site_free, rebuild_csr
   public :: assign_cohort_id, assign_patch_id, reset_step_diagnostics

   type :: site_t
      !==========================================================================================!
      ! 0. THE BLOCKS IT HOLDS                                                                   !
      !==========================================================================================!
      type(cohort_block) :: cohort
      type(patch_block)  :: patch
      !==========================================================================================!
      ! 1. DIMENSIONS AND ID COUNTERS                                                            !
      !==========================================================================================!
      integer(ik)        :: n_pft     = 0_ik
      !----- Monotonic counters for the persistent global ids (never reused within a run). ---!
      integer(ik)        :: next_cohort_id = 1_ik
      integer(ik)        :: next_patch_id  = 1_ik
      !----- Unused: only the restart writes and reads it (one site per polygon, patch areas sum to 1). !
      real(wp)           :: site_area = 1.0_wp
      !==========================================================================================!
      ! 2. SHARED COHORT CLOCK                                                                   !
      !==========================================================================================!
      !----- Global ring-buffer write slot for the moving-average growth (same for all cohorts !
      !      since they all take one sample per step); advanced by update_demography.           !
      integer(ik)        :: growth_hist_pos = 0_ik
      !==========================================================================================!
      ! 3. CLIMATE MEMORY                                                                        !
      !==========================================================================================!
      !----- GROWTH TEMPERATURE for thermal acclimation (docs/science/leaf_gas_exchange.md eq.       !
      !      (4)): an exponential running mean of the daily-mean air temperature, window             !
      !      [leaf_physiology].acclim_window_days (~30 d), which is what Kattge & Knorr (2007)       !
      !      fitted dS against. PROGNOSTIC and slow, so it is written to the state file:             !
      !      rebuilding a month of memory after every restart would make the acclimated optimum      !
      !      jump. Negative means "not yet seeded" -- the first slow step adopts that day's mean     !
      !      outright rather than relaxing from an arbitrary origin. --------------------------------!
      real(wp)           :: t_growth_avg = -1.0_wp       !< [K] running-mean growth temperature
      !==========================================================================================!
      ! 4. FAST-TO-SLOW ACCUMULATORS (the model reads them)                                      !
      !==========================================================================================!
      !----- Site daily-mean air-temperature accumulator for phenology: the fast loop sums the   !
      !      met air temperature over the slow step (reset                                        !
      !      each step), and the slow-loop phenology driver reads the mean (sum/n) then. Transient !
      !      (never restarted). Air temperature is site-uniform (single-site forcing).            !
      real(wp)           :: pheno_tair_sum = 0.0_wp
      integer(ik)        :: pheno_tair_n   = 0_ik
      !==========================================================================================!
      ! 5. TENDENCIES (scratch)                                                                  !
      !==========================================================================================!
      !----- TRANSIENT slow-loop scratch: the per-cohort tendency bundle the driver fills each     !
      !      step and update_cohort_states applies. Site-carried so it allocates once (resize-on-   !
      !      grow) instead of per step; refilled every step, so it is NOT lockstep-reordered and    !
      !      NOT serialized (a scratch buffer, not prognostic state).                               !
      type(cohort_deriv_block) :: deriv
      !==========================================================================================!
      ! 6. OUTPUT ACCUMULATORS (only output reads them; never restarted)                         !
      !==========================================================================================!
      !----- TRANSIENT fast-loop diagnostics that are the same across the polygon -- its forcing   !
      !      (PY_* in meds_diag_vars): one slot, the cohort and patch blocks' polygon twin. ---      !
      type(polygon_diag_block) :: diag
      !----- Site evapotranspiration accumulator [kg/m2 = mm] (site-uniform, single-site): the fast   !
      !      loop sums the area-weighted canopy-air -> atmosphere water-vapour flux * dt_fast over the  !
      !      slow step (reset each step, mirrors the gpp_accum/pheno lifecycle); read as a diagnostic   !
      !      (total_et) at the output tick. Transient (never restarted).                                !
      real(wp)           :: et_accum       = 0.0_wp
      !----- And the water the soil column lost the same way [kg/m2 = mm]: the pond's spill (surface   !
      !      runoff) and the bottom flux (drainage; negative when an aquifer feeds the column). --------!
      real(wp)           :: runoff_surf_accum = 0.0_wp
      real(wp)           :: drainage_accum    = 0.0_wp
      !----- Integrator WORK accumulators: summed over patches and fast sub-steps within a slow      !
      !      step, reset alongside et_accum. These are the COST axis of a scheme comparison --       !
      !      wall-clock is too coarse and too machine-dependent to rank integrators, whereas         !
      !      step/sub-solver counts are exact and reproducible. Transient (diagnostic only,          !
      !      never restarted). Real-valued so the output integrators can average them over a         !
      !      period like any other diagnostic. ------------------------------------------------------!
      real(wp)           :: work_integ_steps = 0.0_wp   !< accepted integrator sub-steps
      real(wp)           :: work_integ_rej   = 0.0_wp   !< rejected integrator steps
      !----- The soil-plant water network's work. A nonzero work_water_nfail is also a correctness     !
      !      flag: a sub-step taken at the shortest length although Newton had not converged keeps      !
      !      water exactly but leaves the stores off their curves. -------------------------------------!
      real(wp)           :: work_water_nsub   = 0.0_wp   !< accepted sub-steps
      real(wp)           :: work_water_newton = 0.0_wp   !< Newton iterations
      real(wp)           :: work_water_reject = 0.0_wp   !< rejected sub-steps
      real(wp)           :: work_water_nfail  = 0.0_wp   !< sub-steps taken although unconverged
      !----- INTEGRATOR-HEALTH accumulators. Same cadence and area weighting as the work             !
      !      counters above, but these measure CORRECTNESS rather than cost: how often a scheme      !
      !      had to leave its own model to stay standing, and how much mass/energy that cost         !
      !      outside the conservation ledger. Without them a run cannot report whether it            !
      !      actually ran the integrator it was configured for (docs/science/numerical_scheme.md     !
      !      §3, `rk45`). ---------------------------------------------------------------------------!
      real(wp)           :: work_rk45_rescue = 0.0_wp   !< dt_fast steps RK45 abandoned and redid with the ARK
      real(wp)           :: work_clamp_stage = 0.0_wp   !< stage-input clamp activations (throwaway states)
   end type site_t

contains

   !=======================================================================================!
   !  LIFECYCLE                                                                            !
   !=======================================================================================!

   !=======================================================================================!
   !  Allocation                                                                            !
   !=======================================================================================!
   subroutine site_alloc(site, n_pft, coh_cap, pat_cap, n_growth_window, allom, pft)
      type(site_t),      intent(out) :: site
      integer(ik),       intent(in)  :: n_pft, coh_cap, pat_cap, n_growth_window
      type(allometry_t), intent(in)  :: allom   !< the site's allometry (cfg%allom of its config)
      type(pft_table_t), intent(in)  :: pft     !< the site's PFT trait table (cfg%pft of its config)
      site%n_pft = n_pft
      site%cohort%allom     = allom
      site%cohort%pft_table = pft
      call alloc_cohort_block(site%cohort, max(coh_cap, 1_ik), max(n_growth_window, 1_ik))
      call alloc_patch_block(site%patch, max(pat_cap, 1_ik), n_pft)
   end subroutine site_alloc

   !----- Each level frees its own fields: the cohort and patch files list them. -------------!
   subroutine site_free(site)
      type(site_t), intent(inout) :: site
      call free_cohort_block(site%cohort)
      call free_patch_block(site%patch)
   end subroutine site_free

   !=======================================================================================!
   !  STRUCTURE                                                                            !
   !=======================================================================================!

   !=======================================================================================!
   !  Rebuild the CSR patch map by stably regrouping cohorts by owner_patch (counting sort).!
   !=======================================================================================!
   subroutine rebuild_csr(site)
      type(site_t), intent(inout) :: site
      integer(ik), allocatable :: perm(:), pos(:), slot(:)
      integer(ik)              :: i, ip, np, nc
      np = site%patch%n
      nc = site%cohort%n
      if (np < 1_ik) then                                  ! no patches: nothing to group
         if (nc > 0_ik) error stop 'rebuild_csr: cohorts present with no patches'
         return
      end if
      site%patch%cohort_count(1:np) = 0_ik
      do i = 1_ik, nc
         ip = site%cohort%owner_patch(i)
         if (ip < 1_ik .or. ip > np) error stop 'rebuild_csr: owner_patch out of range'
         site%patch%cohort_count(ip) = site%patch%cohort_count(ip) + 1_ik
      end do
      site%patch%cohort_offset(1) = 1_ik
      do ip = 2_ik, np
         site%patch%cohort_offset(ip) = site%patch%cohort_offset(ip-1) + site%patch%cohort_count(ip-1)
      end do
      !----- Stable placement: preserve within-patch order. -------------------------------!
      allocate(pos(np), slot(nc), perm(nc))
      pos(1:np) = site%patch%cohort_offset(1:np)
      do i = 1_ik, nc
         ip = site%cohort%owner_patch(i)
         slot(i) = pos(ip)
         pos(ip) = pos(ip) + 1_ik
      end do
      do i = 1_ik, nc
         perm(slot(i)) = i
      end do
      call reorder_cohort_block(site%cohort, perm, nc)
   end subroutine rebuild_csr

   !=======================================================================================!
   !  IDENTITY                                                                             !
   !=======================================================================================!

   !=======================================================================================!
   !  Persistent identity: stamp a slot with the next unused global id and bump the         !
   !  monotonic counter. Ids are never reused within a run, so a terminated/fused-away id    !
   !  simply never reappears -- exactly what an external tracker needs.                      !
   !=======================================================================================!
   subroutine assign_cohort_id(site, slot)
      type(site_t), intent(inout) :: site
      integer(ik),  intent(in)    :: slot
      site%cohort%global_id(slot) = site%next_cohort_id
      site%next_cohort_id = site%next_cohort_id + 1_ik
   end subroutine assign_cohort_id

   subroutine assign_patch_id(site, slot)
      type(site_t), intent(inout) :: site
      integer(ik),  intent(in)    :: slot
      site%patch%global_id(slot) = site%next_patch_id
      site%next_patch_id = site%next_patch_id + 1_ik
   end subroutine assign_patch_id

   !=======================================================================================!
   !  PER STEP                                                                             !
   !=======================================================================================!

   !----- Zero the step diagnostic blocks -- per cohort, per patch, per polygon -- after the       !
   !      output has read them. They accumulate from one output tick to the next: the fast loop    !
   !      adds each sub-step, the slow step its rows, and a calendar boundary between two steps    !
   !      its events (disturbance, culls), which therefore reach the record of the step that       !
   !      follows. The slow cohort block is SET each step instead and needs no reset here. --------!
   subroutine reset_step_diagnostics(site)
      type(site_t), intent(inout) :: site
      call cohort_diag_reset(site%cohort%diag)
      call patch_diag_reset(site%patch%diag)
      call polygon_diag_reset(site%diag)
   end subroutine reset_step_diagnostics

end module meds_site_state_vars
