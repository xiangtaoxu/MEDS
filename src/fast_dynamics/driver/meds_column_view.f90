! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_column_view -- the ONE way a column_cohort_t gets built.                            !
!                                                                                          !
! `column_cohort_t` is the per-patch demographic slice the fast loop reads: a read-only    !
! SoA of one patch's cohorts. It is not state -- every field is a copy of, or a per-ground !
! index over, a field the cohort block already owns -- so the only question it raises is   !
! whether two different fillers can disagree.                                              !
!                                                                                          !
! WHY ONE FILLER. A view assembled field by field need not be allometrically consistent    !
! (dbh 20 cm with a 16 m height, a leaf area of 10 m2/plant and an LAI of 3 that does not  !
! equal nplant*leaf_area), and a field it forgets keeps a zero or stale value: an unset    !
! `bwood` gives the tree a wrong wood heat capacity. A fixture that cannot represent a     !
! real tree cannot falsify anything about one.                                             !
!                                                                                          !
! So there are exactly two entry points, and both go through the cohort block:             !
!   * copy_column_cohort -- production: copy one patch's contiguous CSR section into the   !
!                           view.                                                          !
!   * column_cohort_init -- tests and probes: INITIALIZE a cohort block through the        !
!                           canonical birth path (init_cohort, which sizes it with         !
!                           set_cohort_size), then copy that, so a fixture tree is         !
!                           on-allometry by construction rather than by the author         !
!                           remembering to make it so.                                     !
!                                                                                          !
! Adding a per-cohort INPUT to column physics is therefore: one field on the cohort block, !
! one on column_cohort_t, one line here. Tests inherit it with no edit.                    !
!==========================================================================================!
module meds_column_view
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : cm2_to_m2, C2B_WOOD, t_3ple
   use meds_cohort_state_vars, only : cohort_block, alloc_cohort_block, init_cohort, set_cohort_size
   use meds_pft_params,       only : pft_table_t
   use meds_config,           only : hydraulics_config_t
   use meds_fast_types,       only : column_cohort_t, ensure_column_cohort_capacity
   implicit none
   private

   public :: copy_column_cohort, column_cohort_init

contains

   !---------------------------------------------------------------------------------------!
   ! Copy one patch's contiguous CSR cohort section [i0, i0+ncoh-1] into the fast loop's      !
   ! per-patch view. Every line is a plain copy or a per-ground index: nothing is DERIVED     !
   ! here, because the derived geometry is cached on the cohort block beside height/agb and   !
   ! refreshed whenever the size changes (set_cohort_wood_allometry).                         !
   !                                                                                          !
   ! The per-ground indices are formed HERE rather than cached, so no stored value can carry   !
   ! a stale plant density: mortality changes nplant every step without touching geometry.     !
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
         cc%bsap(j)          = cohort%sapwood_carbon(i)   ! sapwood ring -> HYDRAULICS
         cc%bwood(j)         = cohort%wood_carbon(i)      ! ALL wood      -> THERMAL store
         cc%sap_area(j)      = cohort%sapwood_area(i) * cm2_to_m2   ! [cm2] -> [m2] (hydraulics)
         !----- Plastic leaf capacities + yesterday's daily-max leaf potential. -----------!
         cc%vcmax25(j)       = cohort%vcmax25(i)
         cc%rd25(j)          = cohort%rd25(i)
         cc%dmax_psi_leaf(j) = cohort%dmax_psi_leaf(i)
         !----- Canopy-element geometry and respiration traits (the PFT's, from the block's table). !
         cc%leaf_width(j)    = pft%leaf_width(pf)
         cc%branch_diam(j)   = pft%branch_diameter(pf)
         cc%aboveground_frac(j) = pft%aboveground_frac(pf)
         cc%is_woody(j)           = pft%is_woody(pf)
         cc%stem_resp_factor25(j) = pft%stem_resp_factor25(pf)
         cc%stem_resp_sapwood25(j) = pft%stem_resp_sapwood25(pf)
         !----- Sapwood volume [m3/plant]: its carbon as dry mass, over the PFT's wood density   !
         !      [g/cm3 = 1000 kg/m3]. Biomass grows as density^0.976, so the volume barely      !
         !      depends on the density: it is the stem's geometry. ---------------------------!
         cc%sap_volume(j)    = cohort%sapwood_carbon(i) * C2B_WOOD / (1.0e3_wp * pft%wood_density(pf))
         cc%root_resp_factor25(j) = pft%root_resp_factor25(pf)
      end do
      end associate
   end subroutine copy_column_cohort

   !---------------------------------------------------------------------------------------!
   ! A fixture view of `n` cohorts, on-allometry BY CONSTRUCTION: born through init_cohort    !
   ! and sized by set_cohort_size, exactly as a real recruit is, then gathered by the          !
   ! production path above. The caller supplies only what a real cohort is born from --        !
   ! its PFT, its diameter and its density -- and cannot invent a tree that does not exist.    !
   !---------------------------------------------------------------------------------------!
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
      block%pft_table = pft               ! as site_alloc sets it
      block%n = n
      do j = 1_ik, n
         !----- The view carries no tissue water or temperature (the column's stores are the         !
         !      caller's), so the [hydraulics] defaults for the curves and any birth state reach       !
         !      nothing: the fixture is born saturated (0 MPa) at the triple point. ------------------!
         call init_cohort(block, j, pft, ipft(j), 1_ik, nplant(j), dbh(j), hydraulics_config_t(),      &
                          t_3ple, 0.0_wp)
      end do
      call copy_column_cohort(cc, block, 1_ik, n)
   end subroutine column_cohort_init

end module meds_column_view
