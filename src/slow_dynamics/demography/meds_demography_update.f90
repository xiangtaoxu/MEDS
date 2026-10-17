! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_demography_update -- APPLY/refresh of (derived) cohort + patch state.                     !
!                                                                                          !
! update_cohort_states is the PURE APPLIER: it advances the cohort's PROGNOSTIC state -- the four  !
! carbon pools (state += rate*dt, round-off below zero snapped to 0) and nplant (nplant *=         !
! exp(dln*dt), log-space, stopped at the floor) -- by the per-cohort TENDENCY bundle               !
! (cohort_deriv_block) the slow-loop driver computed this step, and then DERIVES the geometry      !
! (dbh, height, basal area, AGB, leaf area and the wood caches) from the committed carbon through  !
! the allometry. Geometry is never integrated: it is a function of the carbon, so it cannot drift  !
! from it or leave its range (a leaf area integrated by its own rate can go negative on a bare     !
! cohort whose leaf carbon sits at 0). It is MODE-AGNOSTIC: whether the tendencies came from the   !
! carbon flip or the empirical growth law is the driver's concern, invisible here. The derived-type !
! SEAM (host) unpacks the SoA and calls the bare-array KERNEL, which is wrapped in an OpenMP       !
! `target` region so it offloads cleanly while the host keeps ALL state in normal memory.          !
!                                                                                          !
! fill_cohort_deriv is the shared per-cohort TENDENCY BUILDER used by BOTH the carbon (driver) and !
! empirical (C-API) computers: given the new carbon pools (however the caller derived them) it     !
! backs out their time-derivatives and advances the growth-average ring buffer with the caller's   !
! sample -- so the field list and the fill/evict logic live in ONE place.                          !
!                                                                                          !
! update_patch_states advances patch age; the soil-carbon pools, the other slow per-patch state,   !
! advance in meds_biogeochem_dynamics. update_overtopping_lai fills the stored competition         !
! context from the sorted stand.                                                                   !
!==========================================================================================!
module meds_demography_update
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num, lnexp_min, lnexp_max
   use meds_site_state_types, only : site_t, patch_block, cohort_block, cohort_deriv_block,      &
                                     set_cohort_size_from_carbon
   implicit none
   private

   public :: update_cohort_states, fill_cohort_deriv, update_patch_states, update_overtopping_lai

   !----- How far past its bound a quantity may land from round-off alone, in epsilons of the      !
   !      magnitudes it touches. A pool's tendency is (new - old)/dt, so the applier reaches a pool  !
   !      the driver emptied through four roundings (subtract, divide, multiply, add), each within   !
   !      half an epsilon: 4 bounds that with a factor of two to spare. A density decline aimed at   !
   !      the floor is as many (the log of a quotient, a rate divided and multiplied by dt). ------!
   real(wp), parameter :: ROUNDOFF_EPSILONS = 4.0_wp

contains

   !---------------------------------------------------------------------------------------!
   ! Build the per-cohort tendency for slot i from the supplied NEW carbon pools (the caller     !
   ! computed them from NPP or from the empirical forward allometry). Each pool's rate is        !
   ! (new - old)/dt; mortality is the log-space rate -mortality. It ALSO advances                 !
   ! the moving-average ring buffer with `ring_sample` (carbon: the realized dbh rate; empirical: !
   ! the supplied growth rate) -- the fill-while-filling / evict-when-full logic that both callers !
   ! share. Mutates cohort's ring buffer (growth_*); fills deriv(i).                             !
   !---------------------------------------------------------------------------------------!
   subroutine fill_cohort_deriv(cohort, i, deriv, dt_yr, mortality, ring_sample,                   &
                                leaf_carbon_new, fineroot_carbon_new, wood_carbon_new,              &
                                nonstructural_carbon_new, n_window, hist_pos)
      type(cohort_block),       intent(inout) :: cohort
      integer(ik),              intent(in)    :: i
      type(cohort_deriv_block), intent(inout) :: deriv
      real(wp),                 intent(in)    :: dt_yr, mortality, ring_sample
      real(wp),                 intent(in)    :: leaf_carbon_new, fineroot_carbon_new, wood_carbon_new
      real(wp),                 intent(in)    :: nonstructural_carbon_new
      integer(ik),              intent(in)    :: n_window, hist_pos

      !----- Back out each prognostic pool's time-derivative for the applier. ---------------!
      deriv%d_leaf_carbon_dt(i)          = (leaf_carbon_new          - cohort%leaf_carbon(i))          / dt_yr
      deriv%d_fineroot_carbon_dt(i)      = (fineroot_carbon_new      - cohort%fineroot_carbon(i))      / dt_yr
      deriv%d_wood_carbon_dt(i)          = (wood_carbon_new          - cohort%wood_carbon(i))          / dt_yr
      deriv%d_nonstructural_carbon_dt(i) = (nonstructural_carbon_new - cohort%nonstructural_carbon(i)) / dt_yr
      deriv%dln_nplant_dt(i)             = -mortality
      !----- Advance the moving-average ring buffer with the caller's sample (fill or evict), and   !
      !      take the window's sum from the buffer itself. A running sum (add the new sample, take  !
      !      off the evicted one) carries the round-off of every sample it has seen, so a window of  !
      !      zeros after a large sample could read slightly negative. The slots this cohort          !
      !      has not written yet are zero: init_cohort clears the column. -----------------------!
      if (cohort%growth_count(i) < n_window) cohort%growth_count(i) = cohort%growth_count(i) + 1_ik
      cohort%growth_hist(hist_pos, i) = ring_sample
      cohort%growth_accum(i) = sum(cohort%growth_hist(1:n_window, i))
      cohort%growth_avg(i) = cohort%growth_accum(i) / real(cohort%growth_count(i), wp)
   end subroutine fill_cohort_deriv

   !---------------------------------------------------------------------------------------!
   ! Advance every cohort's prognostic + cached-geometry state by the supplied per-cohort     !
   ! TIME-DERIVATIVES for one step of length dt_yr. Derived-type SEAM (host): unpack the SoA + !
   ! the bundle and call the bare-array offload kernel below.                                  !
   !---------------------------------------------------------------------------------------!
   subroutine update_cohort_states(cohort, deriv, dt_yr, negligible)
      type(cohort_block),       intent(inout) :: cohort
      type(cohort_deriv_block), intent(in)    :: deriv
      real(wp),                 intent(in)    :: dt_yr, negligible
      integer(ik) :: i
      call update_cohort_states_kernel(cohort%n, cohort%leaf_carbon, cohort%fineroot_carbon,      &
                                       cohort%wood_carbon, cohort%nonstructural_carbon,            &
                                       cohort%nplant, deriv%d_leaf_carbon_dt,                      &
                                       deriv%d_fineroot_carbon_dt, deriv%d_wood_carbon_dt,          &
                                       deriv%d_nonstructural_carbon_dt, deriv%dln_nplant_dt,        &
                                       dt_yr, negligible)
      !----- The geometry follows the committed carbon: dbh from wood carbon, leaf area from leaf  !
      !      carbon x SLA, and the wood caches from that dbh -- the single flip fusion and fission  !
      !      already use. Host-side and outside the offload region: it is branchy allometry, and    !
      !      the kernel above is arithmetic-only by design.                                         !
      do i = 1_ik, cohort%n
         call set_cohort_size_from_carbon(cohort, i)
      end do
   end subroutine update_cohort_states

   !---------------------------------------------------------------------------------------!
   ! Bare-array OFFLOAD kernel: the carbon pools are additive (state += rate*dt, advance_pool);  !
   ! nplant follows the log-space mortality rate (nplant *= exp(dln*dt)), stopped at the floor   !
   ! `negligible` (advance_density). Arithmetic-only over plain arrays, so the `map` clauses are  !
   ! clean and the host keeps all state in normal memory.                                        !
   !---------------------------------------------------------------------------------------!
   subroutine update_cohort_states_kernel(n, leaf_carbon, fineroot_carbon, wood_carbon,           &
                                   nonstructural_carbon, nplant, d_leaf_carbon_dt,                  &
                                   d_fineroot_carbon_dt, d_wood_carbon_dt, d_nonstructural_carbon_dt,&
                                   dln_nplant_dt, dt_yr, negligible)
      integer(ik), intent(in)    :: n
      real(wp),    intent(inout) :: leaf_carbon(:), fineroot_carbon(:), wood_carbon(:), nonstructural_carbon(:)
      real(wp),    intent(inout) :: nplant(:)
      real(wp),    intent(in)    :: d_leaf_carbon_dt(:), d_fineroot_carbon_dt(:), d_wood_carbon_dt(:)
      real(wp),    intent(in)    :: d_nonstructural_carbon_dt(:), dln_nplant_dt(:)
      real(wp),    intent(in)    :: dt_yr, negligible
      integer(ik) :: i

      !$omp target teams distribute parallel do simd                                        &
      !$omp&        map(to: d_leaf_carbon_dt, d_fineroot_carbon_dt, d_wood_carbon_dt,             &
      !$omp&               d_nonstructural_carbon_dt, dln_nplant_dt)                              &
      !$omp&        map(tofrom: leaf_carbon, fineroot_carbon, wood_carbon, nonstructural_carbon,  &
      !$omp&                    nplant)
      do i = 1_ik, n
         leaf_carbon(i)          = advance_pool(leaf_carbon(i),          d_leaf_carbon_dt(i)          * dt_yr)
         fineroot_carbon(i)      = advance_pool(fineroot_carbon(i),      d_fineroot_carbon_dt(i)      * dt_yr)
         wood_carbon(i)          = advance_pool(wood_carbon(i),          d_wood_carbon_dt(i)          * dt_yr)
         nonstructural_carbon(i) = advance_pool(nonstructural_carbon(i), d_nonstructural_carbon_dt(i) * dt_yr)
         nplant(i)               = advance_density(nplant(i), dln_nplant_dt(i) * dt_yr, negligible)
      end do
   end subroutine update_cohort_states_kernel

   !---------------------------------------------------------------------------------------!
   ! A cohort's density after the log-space change `dln` over the step (nplant *= exp(dln)),   !
   ! stopped AT the density floor. A decline that reaches the floor, to within the round-off   !
   ! of comparing the two in log space, ends exactly on it: nplant*exp(log(floor/nplant)) is   !
   ! the floor only to an epsilon either side, and the monthly cull, which removes a cohort at !
   ! the floor, must find it there, not an epsilon above it for ever. A cohort already at or   !
   ! below the floor is not taken lower. The `1 +` counts the quotient's own rounding, half an !
   ! epsilon in its log whatever the log's size.                                               !
   !---------------------------------------------------------------------------------------!
   elemental pure real(wp) function advance_density(nplant, dln, floor) result(nplant_new)
      real(wp), intent(in) :: nplant   !< [plant/m2] the density
      real(wp), intent(in) :: dln      !< [-] its change in log space over the step, -mortality*dt
      real(wp), intent(in) :: floor    !< [plant/m2] the density floor, negligible_nplant
      real(wp) :: floor_dln
      floor_dln = min(0.0_wp, log(max(floor, tiny_num) / max(nplant, tiny_num)))
      if (floor_dln < 0.0_wp .and. dln - floor_dln <= ROUNDOFF_EPSILONS * epsilon(dln)              &
                                                     * (1.0_wp + abs(floor_dln) + abs(dln))) then
         nplant_new = floor
      else
         nplant_new = nplant * exp(min(max(max(dln, floor_dln), lnexp_min), lnexp_max))
      end if
   end function advance_density

   !---------------------------------------------------------------------------------------!
   ! A carbon pool advanced by its change over the step. A result below zero by no more than    !
   ! ROUNDOFF_EPSILONS epsilons of the pool and its change is a pool the driver emptied, carried !
   ! past zero by rounding: it is 0. Anything further below zero is an overdraw -- an operator   !
   ! that took more than the pool held, which is a bug -- and is left as it is, so the bounds    !
   ! guard (cohort_bad_field) stops the run on this step, naming the pool, instead of a clamp    !
   ! creating the carbon the overdraw took.                                                       !
   !---------------------------------------------------------------------------------------!
   elemental pure real(wp) function advance_pool(pool, change) result(pool_new)
      real(wp), intent(in) :: pool, change    !< [kgC/plant] the pool, and its change over the step
      pool_new = pool + change
      if (pool_new < 0.0_wp .and.                                                                  &
          -pool_new <= ROUNDOFF_EPSILONS * epsilon(pool_new) * (abs(pool) + abs(change))) pool_new = 0.0_wp
   end function advance_pool

   !---------------------------------------------------------------------------------------!
   ! Advance the slow per-PATCH state by one step of length dt_yr: patch AGEING, the             !
   ! patch-level analogue of update_cohort_states. The CENTURY soil-carbon pools, the other      !
   ! slow per-patch state, are advanced by meds_biogeochem_dynamics (soil_carbon_step is         !
   ! their sole writer).                                                                         !
   !---------------------------------------------------------------------------------------!
   subroutine update_patch_states(patch, dt_yr)
      type(patch_block), intent(inout) :: patch
      real(wp),          intent(in)    :: dt_yr
      integer(ik) :: ip
      do ip = 1_ik, patch%n
         patch%age(ip) = patch%age(ip) + dt_yr
      end do
   end subroutine update_patch_states

   !---------------------------------------------------------------------------------------!
   ! Fill cohort%overtopping_lai for every cohort from the current (height-sorted) stand:    !
   ! for each patch, a single top-down pass accumulates the overtopping leaf-area index (the !
   ! cumulative LAI of all TALLER cohorts), with equal-height cohorts sharing one layer      !
   ! (co-dominance). A pure GEOMETRIC reduction -- NO rate law, NO plant physiology.          !
   !---------------------------------------------------------------------------------------!
   subroutine update_overtopping_lai(site)
      type(site_t), intent(inout) :: site
      integer(ik) :: ip, i0, i1, i, j, k
      real(wp)    :: cum, over_lai, layer_lai

      associate (cohort => site%cohort, patch => site%patch)
         do ip = 1_ik, patch%n
            i0  = patch%cohort_offset(ip)
            i1  = i0 + patch%cohort_count(ip) - 1_ik
            cum = 0.0_wp                                   ! overtopping LAI accumulator
            i   = i0
            do while (i <= i1)
               !----- Equal-height cohorts are CO-DOMINANT here: they share the overtopping LAI  !
               !      and do not shade each other, so the result does not depend on the arbitrary   !
               !      within-layer order. Exact-equality is the right test -- recruits are born from !
               !      one shared `recruit_dbh`, so they tie bit-exactly.                             !
               !                                                                                  !
               !      The light capture agrees: canopy_radiation puts exactly tied cohorts in one   !
               !      radiative-transfer layer. What this sweep drives is the light-plastic         !
               !      traits and the demography's competition index.                                !
               k = i
               do while (k < i1)
                  if (cohort%height(k + 1_ik) /= cohort%height(i)) exit
                  k = k + 1_ik
               end do
               over_lai  = cum
               layer_lai = 0.0_wp
               do j = i, k
                  cohort%overtopping_lai(j) = over_lai
                  layer_lai = layer_lai + cohort%nplant(j) * cohort%leaf_area(j)
               end do
               cum = cum + layer_lai                        ! whole layer shades those below
               i   = k + 1_ik
            end do
         end do
      end associate
   end subroutine update_overtopping_lai

end module meds_demography_update
