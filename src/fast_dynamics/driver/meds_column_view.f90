! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_column_view -- the one way a column_cohort_t is built: the read-only per-patch view !
! of the cohort block that the fast loop reads. copy_column_cohort copies one patch's      !
! cohorts; column_cohort_init builds a test fixture through the birth path (init_cohort),  !
! then copies it, so a fixture tree is on-allometry by construction and has no field left  !
! unset. A new per-cohort input to column physics is one field on the cohort block, one on !
! column_cohort_t and one line here.                                                       !
!==========================================================================================!
module meds_column_view
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : cm2_to_m2, c2b_wood, t_3ple
   use meds_cohort_state_vars, only : cohort_block, alloc_cohort_block, init_cohort
   use meds_pft_params,       only : pft_table_t
   use meds_fast_types,       only : column_cohort_t, ensure_column_cohort_capacity
   implicit none
   private

   public :: copy_column_cohort
   public :: column_cohort_init   ! used only by the tests

contains

   !---------------------------------------------------------------------------------------!
   ! Copy one patch's contiguous cohort section [i0, i0+ncoh-1] into the view. Each line   !
   ! is a copy or a per-ground index; the geometry is cached on the cohort block           !
   ! (set_cohort_size_from_carbon), and the indices are formed here, so none carries a     !
   ! stale plant density.                                                                  !
   !---------------------------------------------------------------------------------------!
   subroutine copy_column_cohort(cc, cohort, i0, ncoh)
      type(column_cohort_t), intent(inout) :: cc
      type(cohort_block),    intent(in)    :: cohort
      integer(ik),           intent(in)    :: i0      !< first cohort of the patch (CSR offset)
      integer(ik),           intent(in)    :: ncoh    !< cohorts in the patch
      integer(ik) :: i, j, pf
      call ensure_column_cohort_capacity(cc, ncoh)
      associate (pft => cohort%pft_table)
      do j = 1_ik, ncoh
         i = i0 + j - 1_ik
         pf = cohort%pft(i)
         !----- Demographic inputs. -----------------------------------------------------!
         cc%pft(j)           = cohort%pft(i)
         cc%nplant(j)        = cohort%nplant(i)
         cc%dbh(j)           = cohort%dbh(i)
         cc%height(j)        = cohort%height(i)
         cc%leaf_area(j)     = cohort%leaf_area(i)
         !----- Per-ground indices: nplant * per-plant area, formed here, never stored. ---!
         cc%lai(j)           = cohort%nplant(i) * cohort%leaf_area(i)
         cc%wai(j)           = cohort%nplant(i) * cohort%wood_area(i)
         !----- Carbon pools the fast loop reads (capacitance + thermal stores). ----------!
         cc%bleaf(j)         = cohort%leaf_carbon(i)
         cc%broot(j)         = cohort%fineroot_carbon(i)
         cc%bsap(j)          = cohort%sapwood_carbon(i)   ! sapwood ring -> hydraulics
         cc%bwood(j)         = cohort%wood_carbon(i)      ! all wood     -> thermal store
         cc%sap_area(j)      = cohort%sapwood_area(i) * cm2_to_m2   ! [cm2] -> [m2] (hydraulics)
         !----- Plastic leaf capacities + yesterday's daily-max leaf potential. -----------!
         cc%vcmax25(j)       = cohort%vcmax25(i)
         cc%rd25(j)          = cohort%rd25(i)
         cc%dmax_psi_leaf(j) = cohort%dmax_psi_leaf(i)
         !----- Canopy-element geometry and respiration traits (the PFT's, from the block's table). !
         cc%leaf_width(j)    = pft%leaf_width(pf)
         cc%branch_diam(j)   = pft%branch_diameter(pf)
         cc%aboveground_frac(j) = pft%allom(pf)%aboveground_frac
         cc%is_woody(j)           = pft%is_woody(pf)
         cc%stem_resp_factor25(j) = pft%stem_resp_factor25(pf)
         cc%stem_resp_sapwood25(j) = pft%stem_resp_sapwood25(pf)
         !----- Sapwood volume [m3/plant]: its carbon as dry mass over the PFT's wood density   !
         !      [g/cm3 = 1000 kg/m3]; biomass grows as density^0.976, so it barely depends on it. !
         cc%sap_volume(j)    = cohort%sapwood_carbon(i) * c2b_wood / (1.0e3_wp * pft%wood_density(pf))
         cc%root_resp_factor25(j) = pft%root_resp_factor25(pf)
      end do
      end associate
   end subroutine copy_column_cohort

   !----- A fixture view of n cohorts, born through init_cohort as a real recruit is and gathered  !
   !      by copy_column_cohort: the caller gives only each cohort's PFT, diameter and density. ---!
   subroutine column_cohort_init(cc, pft, ipft, dbh, nplant)
      type(column_cohort_t), intent(inout) :: cc
      type(pft_table_t),     intent(in)    :: pft
      integer(ik),           intent(in)    :: ipft(:)    !< (n) PFT index per cohort
      real(wp),              intent(in)    :: dbh(:)     !< (n) [cm]       diameter
      real(wp),              intent(in)    :: nplant(:)  !< (n) [plant/m2] density
      type(cohort_block) :: block
      integer(ik)        :: j, n
      n = int(size(ipft), ik)
      call alloc_cohort_block(block, n, 1_ik)   ! growth ring buffer irrelevant to a fixture
      block%pft_table = pft               ! as alloc_polygon sets it
      block%n = n
      do j = 1_ik, n
         !----- The view carries no tissue water or temperature, so the birth state reaches nothing: !
         !      saturated (0 MPa) at the triple point. --------------------------------------------!
         call init_cohort(block, j, pft, ipft(j), 1_ik, nplant(j), dbh(j), t_3ple, 0.0_wp)
      end do
      call copy_column_cohort(cc, block, 1_ik, n)
   end subroutine column_cohort_init

end module meds_column_view
