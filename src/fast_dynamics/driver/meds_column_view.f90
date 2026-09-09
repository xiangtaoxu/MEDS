!==========================================================================================!
! meds_column_view -- the ONE way a column_cohort_t gets built.                              !
!                                                                                          !
! `column_cohort_t` is the per-patch demographic slice the fast loop reads: a read-only SoA !
! of one patch's cohorts. It is not state -- every field is a copy of, or a per-ground index  !
! over, a field the cohort block already owns -- so the only question it raises is whether     !
! two different fillers can disagree. They could, and they did.                                !
!                                                                                          !
! THE DEFECT THIS MODULE CLOSES. The production filler copied a contiguous CSR section, but  !
! every column test hand-assembled its own view field by field. Those hand-built trees were   !
! not allometrically consistent -- dbh 20 cm with a 16 m height, a leaf area of 10 m2/plant    !
! and an LAI of 3 that does not equal nplant*leaf_area -- and none of them set `bwood` at all, !
! which is allocated and never initialized, so the wood heat capacity in those tests ran on    !
! uninitialized memory. A fixture that cannot represent a real tree cannot falsify anything     !
! about one.                                                                                     !
!                                                                                          !
! So there are exactly two entry points, and both go through the cohort block:                  !
!   * copy_column_cohort  -- production: copy one patch's contiguous CSR section into the view. !
!   * column_cohort_init  -- tests and probes: INITIALIZE a cohort block through the canonical  !
!                            birth path (init_cohort + set_cohort_size), then copy that, so a    !
!                            fixture tree is on-allometry by construction rather than by the     !
!                            author remembering to make it so.                                    !
!                                                                                          !
! Adding a per-cohort INPUT to column physics is therefore: one field on the cohort block, one   !
! on column_cohort_t, one line here. Tests inherit it with no edit.                               !
!==========================================================================================!
module meds_column_view
   use meds_kinds,            only : wp, ik
   use meds_site_state_types, only : cohort_block, cohort_alloc, init_cohort, set_cohort_size
   use meds_pft_params,       only : pft_table_t
   use meds_fast_types,       only : column_cohort_t, ensure_column_cohort_capacity
   implicit none
   private

   public :: copy_column_cohort, column_cohort_init

contains

   !---------------------------------------------------------------------------------------!
   ! Copy one patch's contiguous CSR cohort section [i0, i0+ncoh-1] into the fast loop's      !
   ! per-patch view. Every line is a plain copy or a per-ground index: nothing is DERIVED     !
   ! here, because the derived geometry is cached on the cohort block beside height/agb and   !
   ! refreshed whenever the size changes (set_cohort_wood_geometry).                          !
   !                                                                                          !
   ! The per-ground indices are formed HERE rather than cached, so no stored value can carry   !
   ! a stale plant density: mortality changes nplant every step without touching geometry.     !
   !---------------------------------------------------------------------------------------!
   subroutine copy_column_cohort(cc, cohort, i0, ncoh)
      type(column_cohort_t), intent(inout) :: cc
      type(cohort_block),    intent(in)    :: cohort
      integer(ik),           intent(in)    :: i0      !< first cohort of the patch (CSR offset)
      integer(ik),           intent(in)    :: ncoh    !< cohorts in the patch
      integer(ik) :: i, j
      call ensure_column_cohort_capacity(cc, ncoh)
      do j = 1_ik, ncoh
         i = i0 + j - 1_ik
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
         cc%sap_area(j)      = cohort%sapwood_area(i)
         !----- Plastic leaf capacities + yesterday's daily-max leaf potential. -----------!
         cc%vcmax25(j)       = cohort%vcmax25(i)
         cc%rd25(j)          = cohort%rd25(i)
         cc%dmax_psi_leaf(j) = cohort%dmax_psi_leaf(i)
         !----- Canopy-element geometry (per-PFT traits, gathered on the block). ----------!
         cc%crown(j)         = cohort%p_crown_area_frac(i)
         cc%leaf_width(j)    = cohort%p_leaf_width(i)
         cc%branch_diam(j)   = cohort%p_branch_diameter(i)
         cc%aboveground_frac(j) = cohort%p_aboveground_frac(i)
      end do
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
      call cohort_alloc(block, n, 1_ik)   ! growth ring buffer irrelevant to a fixture
      block%n = n
      do j = 1_ik, n
         call init_cohort(block, j, pft, ipft(j), 1_ik, nplant(j), dbh(j))
      end do
      call copy_column_cohort(cc, block, 1_ik, n)
   end subroutine column_cohort_init

end module meds_column_view
