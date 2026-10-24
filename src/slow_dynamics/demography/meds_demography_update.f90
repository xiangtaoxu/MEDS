! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_demography_update -- applies the step's rates (docs/science/demography.md §2).      !
! update_cohort_states advances the four carbon pools (advance_pool, eq. (8)) and the      !
! density (advance_density, in log space) by the slow driver's tendencies, then derives    !
! the geometry from the committed carbon. fill_cohort_deriv builds the tendencies and      !
! advances the growth history (eq. (2)), for the carbon and empirical growth paths alike.  !
! update_patch_states ages the patches, and update_overtopping_lai computes eq. (11).      !
!==========================================================================================!
module meds_demography_update
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num, lnexp_min, lnexp_max
   use meds_numerics,  only : equal_height_layers
   use meds_cohort_state_vars, only : cohort_block, cohort_deriv_block, set_cohort_size_from_carbon, &
                                      hold_cohort_tissue_temps
   use meds_patch_state_vars, only : patch_block
   use meds_polygon_state_vars, only : polygon_t
   implicit none
   private

   public :: update_cohort_states, fill_cohort_deriv, advance_pool, update_patch_states, update_overtopping_lai

   !----- How far past its bound a quantity may land from round-off alone, in epsilons of the magnitudes !
   !      it touches (demography.md §2). An emptied pool, or a density decline aimed at the floor, gets  !
   !      there through four roundings of half an epsilon each, so 4 leaves a factor of two to spare. ---!
   real(wp), parameter :: ROUNDOFF_EPSILONS = 4.0_wp

contains

   !---------------------------------------------------------------------------------------!
   ! Cohort i's tendency from the new carbon pools the caller computed (carbon allocation  !
   ! or empirical growth law): each pool's rate is (new - old)/dt, the density's is        !
   ! -mortality in log space. It also adds ring_sample to the growth history (eq. (2)).    !
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

      deriv%d_leaf_carbon_dt(i)          = (leaf_carbon_new          - cohort%leaf_carbon(i))          / dt_yr
      deriv%d_fineroot_carbon_dt(i)      = (fineroot_carbon_new      - cohort%fineroot_carbon(i))      / dt_yr
      deriv%d_wood_carbon_dt(i)          = (wood_carbon_new          - cohort%wood_carbon(i))          / dt_yr
      deriv%d_nonstructural_carbon_dt(i) = (nonstructural_carbon_new - cohort%nonstructural_carbon(i)) / dt_yr
      deriv%dln_nplant_dt(i)             = -mortality
      !----- Write the sample over the oldest and sum the window from the buffer itself, not as !
      !      a running total (demography.md §1); the slots not yet written are 0, cleared by    !
      !      init_cohort. ----------------------------------------------------------------------!
      if (cohort%growth_count(i) < n_window) cohort%growth_count(i) = cohort%growth_count(i) + 1_ik
      cohort%growth_hist(hist_pos, i) = ring_sample
      cohort%growth_accum(i) = sum(cohort%growth_hist(1:n_window, i))
      cohort%growth_avg(i) = cohort%growth_accum(i) / real(cohort%growth_count(i), wp)
   end subroutine fill_cohort_deriv

   !---------------------------------------------------------------------------------------!
   ! Advance every cohort's pools and density by the tendencies over dt_yr, then derive    !
   ! the geometry from the committed carbon (demography.md §2).                            !
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
      !----- The geometry from the committed carbon, as fusion and fission derive it; the tissue !
      !      made or lost came or went at the tissue's temperature, which its energy is held to  !
      !      (conservation.md §4). --------------------------------------------------------------!
      do i = 1_ik, cohort%n
         call set_cohort_size_from_carbon(cohort, i)
         call hold_cohort_tissue_temps(cohort, i)
      end do
   end subroutine update_cohort_states

   !---------------------------------------------------------------------------------------!
   ! Each pool moves by its change over the step (advance_pool) and the density by its     !
   ! log-space change, stopped at the floor `negligible` (advance_density).                !
   ! The loop is an OpenMP target region.                                                  !
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
   ! The density after the log-space change dln, stopped at the floor (demography.md §2).  !
   ! A decline that reaches the floor, to within the round-off of comparing the two in log !
   ! space, ends exactly on it, so the cull finds it there; a cohort already at or below   !
   ! the floor is not taken lower. The `1 +` counts the quotient's own rounding.           !
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
   ! A carbon pool advanced by its change over the step, eq. (8). A result below zero by   !
   ! no more than ROUNDOFF_EPSILONS of |pool| + |change| is round-off and set to 0. An     !
   ! overdraw further below is left as it is, so the bounds check (cohort_bad_field) stops !
   ! the run and names the pool (demography.md §2).                                        !
   !---------------------------------------------------------------------------------------!
   elemental pure real(wp) function advance_pool(pool, change) result(pool_new)
      real(wp), intent(in) :: pool, change    !< [kgC/plant] the pool, and its change over the step
      pool_new = pool + change
      if (pool_new < 0.0_wp .and.                                                                  &
          -pool_new <= ROUNDOFF_EPSILONS * epsilon(pool_new) * (abs(pool) + abs(change))) pool_new = 0.0_wp
   end function advance_pool

   !---------------------------------------------------------------------------------------!
   ! Age every patch by dt_yr. The soil-carbon pools, the other slow per-patch state,      !
   ! advance in meds_biogeochem_dynamics (soil_carbon_step is their only writer).          !
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
   ! Each cohort's overtopping LAI, eq. (11) of demography.md §6: the LAI of the taller    !
   ! cohorts in its patch, in one top-down pass. Equal-height cohorts share a layer and do !
   ! not shade each other, as in the radiation (equal_height_layers).                      !
   !---------------------------------------------------------------------------------------!
   subroutine update_overtopping_lai(poly)
      type(polygon_t), intent(inout) :: poly
      integer(ik) :: ip, i0, i1, i, j, k, il, n_layer
      real(wp)    :: cum, over_lai, layer_lai

      associate (cohort => poly%cohort, patch => poly%patch)
         do ip = 1_ik, patch%n
            i0  = patch%cohort_offset(ip)
            i1  = i0 + patch%cohort_count(ip) - 1_ik
            block
               integer(ik) :: layer(i0:i1)
               call equal_height_layers(cohort%height(i0:i1), i1 - i0 + 1_ik, layer, n_layer)
               cum = 0.0_wp                                ! overtopping LAI accumulator
               i   = i0
               do il = 1_ik, n_layer
                  k = i
                  do while (k < i1)
                     if (layer(k + 1_ik) /= il) exit
                     k = k + 1_ik
                  end do
                  over_lai  = cum
                  layer_lai = 0.0_wp
                  do j = i, k
                     cohort%overtopping_lai(j) = over_lai
                     layer_lai = layer_lai + cohort%nplant(j) * cohort%leaf_area(j)
                  end do
                  cum = cum + layer_lai                     ! whole layer shades those below
                  i   = k + 1_ik
               end do
            end block
         end do
      end associate
   end subroutine update_overtopping_lai

end module meds_demography_update
