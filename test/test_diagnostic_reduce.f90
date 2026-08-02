!==========================================================================================!
! test_diagnostic_reduce -- the WEIGHTED scale-aggregation layer (MEDS_IO_V01_PLAN.md section 3.3   !
! and test plan items 2, 3, 3b, 4, 6).                                                        !
!                                                                                          !
! What these tests are FOR. The reducer's whole job is to apply the right weight to the right   !
! quantity, and the two ways that goes wrong are both silent: an EXTENSIVE quantity aggregated    !
! as a mean (or an INTENSIVE one as a sum) produces a plausible number with no error anywhere,     !
! and an empty bin reported as 0 is indistinguishable in the output file from a real measurement   !
! of zero. So the assertions below are deliberately chosen so that a wrong weight cannot pass:     !
! the fixture uses DIFFERENT nplant per cohort and different values per PFT and size class, which  !
! makes the sum, the stem-weighted mean and the basal-area-weighted mean three distinct numbers.   !
!==========================================================================================!
program test_diagnostic_reduce
   use meds_kinds,              only : wp, ik
   use meds_core_state_types,   only : site_t, site_alloc, site_free
   use meds_column_state_types, only : n_soil_layer_max
   use meds_diagnostic_reduce,  only : reduce_cohort_to_site, reduce_cohort_to_patch,            &
                                       reduce_cohort_to_pft, reduce_cohort_to_size,             &
                                       reduce_patch_to_site, reduce_patch_column_to_site,       &
                                       gather_patch_columns,                                     &
                                       W_NONE, W_NPLANT, W_LEAF_AREA, W_BASAL_AREA,             &
                                       total_agb, total_lai, total_nplant
   use meds_diagnostic_kernels, only : dbh_class_index, soil_wetness, cohort_gsc, cohort_wue,    &
                                       safe_ratio
   use meds_core_diag_types,    only : cohort_diag_block, cohort_diag_alloc, cohort_diag_reorder, &
                                       cohort_diag_copy_slot, cohort_diag_fuse,                  &
                                       cohort_diag_clear_slot, cohort_diag_value,                &
                                       CDIAG_FUSE, N_CDIAG, CD_LEAF_TEMP, CD_SAPFLOW, CD_ABS_SW
   use meds_test_support,       only : check, check_close, banner
   implicit none

   call banner('diagnostic_reduce')
   call test_extensive_vs_intensive()
   call test_pft_and_size_closure()
   call test_empty_sets()
   call test_class_edges()
   call test_soil_columns()
   call test_kernels()
   call test_diag_lockstep()
   write(*,'(a)') 'test_diagnostic_reduce: ALL PASSED'

contains

   !----- Two patches, three cohorts, deliberately unequal nplant / basal_area so that a SUM, a  !
   !      stem-weighted MEAN and a basal-area-weighted MEAN are three different numbers.  -------!
   subroutine build_fixture(site)
      type(site_t), intent(out) :: site
      call site_alloc(site, 3_ik, 16_ik, 4_ik, 4_ik)
      site%n_pft = 3_ik
      !----- patch 1 (area 0.75): cohorts 1,2 ; patch 2 (area 0.25): cohort 3. -----!
      site%patch%n = 2_ik
      site%patch%area(1) = 0.75_wp ; site%patch%area(2) = 0.25_wp
      site%patch%cohort_offset(1) = 1_ik ; site%patch%cohort_count(1) = 2_ik
      site%patch%cohort_offset(2) = 3_ik ; site%patch%cohort_count(2) = 1_ik
      site%cohort%n = 3_ik
      site%cohort%owner_patch(1:3) = [1_ik, 1_ik, 2_ik]
      site%cohort%pft(1:3)         = [1_ik, 2_ik, 3_ik]
      site%cohort%nplant(1:3)      = [2.0_wp,  1.0_wp,  4.0_wp]
      site%cohort%agb(1:3)         = [10.0_wp, 20.0_wp, 30.0_wp]
      site%cohort%leaf_area(1:3)   = [ 1.0_wp,  3.0_wp,  5.0_wp]
      site%cohort%basal_area(1:3)  = [ 5.0_wp,  1.0_wp,  2.0_wp]
      site%cohort%dbh(1:3)         = [ 5.0_wp, 15.0_wp, 45.0_wp]   ! -> classes 1, 2, 4
   end subroutine build_fixture

   !=======================================================================================!
   !  EXTENSIVE vs INTENSIVE. The same field, reduced with the two different contracts, must   !
   !  give the two different right answers -- and swapping the weight must CHANGE the answer,   !
   !  which is what makes the assertion able to fail.                                           !
   !=======================================================================================!
   subroutine test_extensive_vs_intensive()
      type(site_t) :: site
      real(wp) :: v, v2, out(8)
      logical  :: ok, vld(8)
      integer(ik) :: n
      call build_fixture(site)

      !----- EXTENSIVE: agb_site = SUM area*nplant*agb                                        !
      !      = 0.75*(2*10 + 1*20) + 0.25*(4*30) = 0.75*40 + 0.25*120 = 30 + 30 = 60.  --------!
      call reduce_cohort_to_site(site, site%cohort%agb, W_NPLANT, .false., v, ok)
      call check_close(v, 60.0_wp, 1.0e-12_wp, 'agb_site = weighted SUM over nplant')
      call check(ok, 'extensive reduction is always valid')
      call check_close(total_agb(site), 60.0_wp, 1.0e-12_wp, 'total_agb wrapper agrees')

      !----- The same call with mean=.true. is a DIFFERENT number (mean per-plant agb):        !
      !      60 / (0.75*3 + 0.25*4) = 60/3.25.  ---------------------------------------------!
      call reduce_cohort_to_site(site, site%cohort%agb, W_NPLANT, .true., v2, ok)
      call check_close(v2, 60.0_wp/3.25_wp, 1.0e-12_wp, 'mean=.true. gives the per-plant mean')
      call check(abs(v - v2) > 1.0_wp, 'sum and mean are genuinely different (a wrong flag cannot pass)')

      !----- INTENSIVE with a WEIGHT that matters: basal-area-weighted mean dbh                !
      !      w = nplant*basal_area = [10, 1, 8]; area-weighted:                                 !
      !      num = 0.75*(10*5 + 1*15) + 0.25*(8*45) = 0.75*65 + 0.25*360 = 48.75 + 90 = 138.75  !
      !      den = 0.75*(10 + 1)      + 0.25*8      = 8.25 + 2 = 10.25                          !
      call reduce_cohort_to_site(site, site%cohort%dbh, W_BASAL_AREA, .true., v, ok)
      call check_close(v, 138.75_wp/10.25_wp, 1.0e-12_wp, 'BA-weighted mean dbh')
      !----- The STEM-weighted mean of the same field is a different number, so the weight kind  !
      !      is doing real work rather than being decorative.  -------------------------------!
      call reduce_cohort_to_site(site, site%cohort%dbh, W_NPLANT, .true., v2, ok)
      call check(abs(v - v2) > 1.0_wp, 'BA-weighted /= stem-weighted mean dbh')

      !----- LEAF-AREA weighting, the contract canopy temperatures and conductances use. -----!
      call reduce_cohort_to_site(site, site%cohort%dbh, W_LEAF_AREA, .true., v, ok)
      call check(abs(v - v2) > 0.1_wp, 'leaf-area-weighted /= stem-weighted mean dbh')

      !----- W_NONE + SUM is the stem-density contract: 0.75*(2+1) + 0.25*4 = 2.25 + 1 = 3.25. !
      call reduce_cohort_to_site(site, site%cohort%nplant, W_NONE, .false., v, ok)
      call check_close(v, 3.25_wp, 1.0e-12_wp, 'nplant_site = plain area-weighted sum')
      call check_close(total_nplant(site), 3.25_wp, 1.0e-12_wp, 'total_nplant wrapper agrees')

      !----- PATCH axis carries NO area factor (a patch value is per m2 of its OWN ground):    !
      !      patch 1 agb = 2*10 + 1*20 = 40 ; patch 2 = 4*30 = 120.  -------------------------!
      call reduce_cohort_to_patch(site, site%cohort%agb, W_NPLANT, .false., out, vld, n)
      call check(n == 2_ik, 'patch reduction length')
      call check_close(out(1),  40.0_wp, 1.0e-12_wp, 'agb_patch(1) has no area factor')
      call check_close(out(2), 120.0_wp, 1.0e-12_wp, 'agb_patch(2) has no area factor')

      call site_free(site)
   end subroutine test_extensive_vs_intensive

   !=======================================================================================!
   !  CLOSURE. sum over PFTs and sum over size classes must both reproduce the site total     !
   !  EXACTLY -- that identity is the entire reason those axes are trustworthy, and it is the   !
   !  one property a reader will assume without checking.                                      !
   !=======================================================================================!
   subroutine test_pft_and_size_closure()
      type(site_t) :: site
      real(wp) :: site_val, pft_out(8), size_out(8), edges(8)
      logical  :: ok, vld(8)
      integer(ik) :: n
      call build_fixture(site)
      edges(1:5) = [0.0_wp, 10.0_wp, 20.0_wp, 30.0_wp, 100.0_wp]     ! 4 classes

      call reduce_cohort_to_site(site, site%cohort%agb, W_NPLANT, .false., site_val, ok)
      call reduce_cohort_to_pft (site, site%cohort%agb, W_NPLANT, .false., 3_ik, pft_out, vld, n)
      call check(n == 3_ik, 'pft axis length = n_pft')
      call check_close(sum(pft_out(1:3)), site_val, 1.0e-12_wp, 'sum over PFTs == agb_site')

      call reduce_cohort_to_size(site, site%cohort%agb, W_NPLANT, .false., edges, 4_ik,          &
                                 size_out, vld, n)
      call check(n == 4_ik, 'size axis length = n_class')
      call check_close(sum(size_out(1:4)), site_val, 1.0e-12_wp, 'sum over DBH classes == agb_site')

      !----- The stem count must be conserved by the binning too: no cohort split, none dropped. !
      call reduce_cohort_to_site(site, site%cohort%nplant, W_NONE, .false., site_val, ok)
      call reduce_cohort_to_size(site, site%cohort%nplant, W_NONE, .false., edges, 4_ik,         &
                                 size_out, vld, n)
      call check_close(sum(size_out(1:4)), site_val, 1.0e-12_wp, 'sum over classes == nplant_site')
      !----- Each cohort landed in exactly ONE bin (dbh 5, 15, 45 -> classes 1, 2, 4).  -------!
      call check(size_out(3) == 0.0_wp, 'the empty class-3 bin holds nothing')

      !----- LAI closes on both axes as well (a second field, so the test is not agb-specific). !
      call reduce_cohort_to_site(site, site%cohort%leaf_area, W_NPLANT, .false., site_val, ok)
      call reduce_cohort_to_pft (site, site%cohort%leaf_area, W_NPLANT, .false., 3_ik, pft_out, vld, n)
      call check_close(sum(pft_out(1:3)), site_val, 1.0e-12_wp, 'sum over PFTs == lai_site')
      call check_close(total_lai(site), site_val, 1.0e-12_wp, 'total_lai wrapper agrees')

      call site_free(site)
   end subroutine test_pft_and_size_closure

   !=======================================================================================!
   !  EMPTY SETS. A patch with no cohorts, and a MEAN over a PFT with no members, must report   !
   !  invalid (-> _FillValue), never 0/0 and never a bare 0 that a reader would take for data.   !
   !  But an empty-PFT SUM is a true 0, because reporting fill there would break the closure     !
   !  asserted above the moment a PFT went locally extinct. The two conventions are different    !
   !  on purpose and both are pinned here.                                                       !
   !=======================================================================================!
   subroutine test_empty_sets()
      type(site_t) :: site
      real(wp) :: out(8), v
      logical  :: vld(8), ok
      integer(ik) :: n
      call build_fixture(site)
      !----- Make patch 2 empty (its cohort moves to patch 1). -----!
      site%patch%cohort_count(1) = 3_ik
      site%patch%cohort_count(2) = 0_ik
      site%patch%cohort_offset(2) = 4_ik
      site%cohort%owner_patch(3) = 1_ik

      call reduce_cohort_to_patch(site, site%cohort%agb, W_NPLANT, .false., out, vld, n)
      call check(n == 2_ik,      'empty patch still occupies a slot')
      call check(vld(1),         'occupied patch is valid')
      call check(.not. vld(2),   'EMPTY patch is invalid -> _FillValue, not 0')

      !----- A MEAN over a PFT with no members is invalid; a SUM over it is a true 0. ---------!
      site%n_pft = 4_ik                                   ! PFT 4 exists but has no cohorts
      call reduce_cohort_to_pft(site, site%cohort%dbh, W_NPLANT, .true., 4_ik, out, vld, n)
      call check(.not. vld(4), 'MEAN over an empty PFT is invalid')
      call reduce_cohort_to_pft(site, site%cohort%agb, W_NPLANT, .false., 4_ik, out, vld, n)
      call check(vld(4),            'SUM over an empty PFT is VALID (closure would break otherwise)')
      call check(out(4) == 0.0_wp,  'SUM over an empty PFT is exactly 0')

      !----- A site with no cohorts at all: the extensive sum is 0, the mean is invalid. ------!
      site%cohort%n = 0_ik ; site%patch%cohort_count(1:2) = 0_ik
      call reduce_cohort_to_site(site, site%cohort%agb, W_NPLANT, .false., v, ok)
      call check(v == 0.0_wp, 'bare ground: extensive site total is 0')
      call reduce_cohort_to_site(site, site%cohort%dbh, W_BASAL_AREA, .true., v, ok)
      call check(.not. ok,    'bare ground: weighted MEAN is invalid, not 0/0')

      call site_free(site)
   end subroutine test_empty_sets

   !=======================================================================================!
   !  SIZE-CLASS EDGES. Half-open [lo, hi) except the last class, which is CLOSED at the top    !
   !  so the largest tree in the stand is never dropped; below the first edge clamps into bin 1  !
   !  for the same reason. Both are closure-preserving choices, so both are pinned.              !
   !=======================================================================================!
   subroutine test_class_edges()
      real(wp) :: edges(5)
      edges = [0.0_wp, 10.0_wp, 20.0_wp, 30.0_wp, 100.0_wp]
      call check(dbh_class_index( 0.0_wp, edges, 4_ik) == 1_ik, 'exactly the first edge -> bin 1')
      call check(dbh_class_index( 9.9_wp, edges, 4_ik) == 1_ik, 'just below an edge -> lower bin')
      call check(dbh_class_index(10.0_wp, edges, 4_ik) == 2_ik, 'exactly on an edge -> UPPER bin')
      call check(dbh_class_index(29.9_wp, edges, 4_ik) == 3_ik, 'interior bin')
      call check(dbh_class_index(30.0_wp, edges, 4_ik) == 4_ik, 'exactly the last lower edge -> last bin')
      call check(dbh_class_index(1.0e4_wp, edges, 4_ik) == 4_ik, 'above the top edge -> last bin (never dropped)')
      call check(dbh_class_index(-1.0_wp, edges, 4_ik) == 1_ik, 'below the bottom edge -> bin 1 (never dropped)')
   end subroutine test_class_edges

   !=======================================================================================!
   !  SOIL COLUMNS. The area-weighted site column, and the FLATTENED 2-D (layer, patch) slab   !
   !  whose stride is the COMPILE-TIME layer ceiling. The closure below -- area-weighted mean    !
   !  over the patch axis of the 2-D slab == the site column -- is the P1 acceptance test, and   !
   !  it is what proves the flattening index is right. Varying the ACTIVE layer count while the  !
   !  stride stays fixed is the specific mis-mapping it has to rule out.                         !
   !=======================================================================================!
   subroutine test_soil_columns()
      type(site_t) :: site
      real(wp) :: col(n_soil_layer_max, 2), sitecol(n_soil_layer_max)
      real(wp) :: flat(2*n_soil_layer_max), wmean
      logical  :: vld(2*n_soil_layer_max)
      integer(ik) :: n, k, nlayer
      call build_fixture(site)
      nlayer = 5_ik                                     ! ACTIVE layers < the ceiling, on purpose
      do k = 1_ik, nlayer
         col(k, 1) = 0.10_wp * real(k, wp)              ! patch 1 profile
         col(k, 2) = 0.50_wp * real(k, wp)              ! patch 2 profile (deliberately different)
      end do

      call reduce_patch_column_to_site(site, col, nlayer, sitecol, n)
      call check(n == nlayer, 'site column length = active layers')
      !----- area-weighted: 0.75*0.1k + 0.25*0.5k = 0.2k  ------------------------------------!
      do k = 1_ik, nlayer
         call check_close(sitecol(k), 0.2_wp*real(k,wp), 1.0e-12_wp, 'area-weighted site soil column')
      end do

      call gather_patch_columns(site, col, nlayer, flat, vld, n)
      call check(n == 2_ik*n_soil_layer_max, 'flat slab length = n_patch * n_soil_layer_max')
      !----- The STRIDE is the compile-time ceiling, NOT nlayer: patch 2 layer 1 sits at         !
      !      n_soil_layer_max + 1, not nlayer + 1. Striding by the live count is the mis-mapping  !
      !      that would produce plausible shifted profiles instead of an error.  ---------------!
      call check_close(flat(1),                       col(1, 1), 1.0e-12_wp, 'patch 1 layer 1')
      call check_close(flat(nlayer),                  col(nlayer, 1), 1.0e-12_wp, 'patch 1 last active layer')
      call check_close(flat(n_soil_layer_max + 1_ik), col(1, 2), 1.0e-12_wp, 'patch 2 starts at the CEILING stride')
      call check(.not. vld(nlayer + 1_ik), 'inactive layers of patch 1 are invalid -> _FillValue')
      call check(vld(n_soil_layer_max + 1_ik), 'patch 2 active layers are valid')

      !----- The P1 acceptance identity, computed the way a reader would.  --------------------!
      do k = 1_ik, nlayer
         wmean = site%patch%area(1) * flat(k)                                                    &
               + site%patch%area(2) * flat(n_soil_layer_max + k)
         call check_close(wmean, sitecol(k), 1.0e-12_wp, 'area-weighted 2-D slab == site column')
      end do

      call site_free(site)
   end subroutine test_soil_columns

   !=======================================================================================!
   !  The derived-quantity kernels, including the guarded division that keeps a night-time      !
   !  zero denominator from poisoning a period mean with a NaN or an Inf.                       !
   !=======================================================================================!
   subroutine test_kernels()
      call check_close(cohort_gsc(1.6_wp), 1.0_wp, 1.0e-12_wp, 'gsc = gsw / 1.6')
      call check_close(soil_wetness(0.3_wp, 0.1_wp, 0.5_wp), 0.5_wp, 1.0e-12_wp, 'relative saturation')
      call check_close(soil_wetness(0.6_wp, 0.1_wp, 0.5_wp), 1.0_wp, 1.0e-12_wp, 'oversaturated clamps to 1')
      call check_close(soil_wetness(0.0_wp, 0.1_wp, 0.5_wp), 0.0_wp, 1.0e-12_wp, 'below residual clamps to 0')
      call check_close(safe_ratio(4.0_wp, 2.0_wp), 2.0_wp, 1.0e-12_wp, 'ordinary ratio')
      call check(safe_ratio(4.0_wp, 0.0_wp) == 0.0_wp, 'zero denominator -> 0, not NaN/Inf')
      call check(cohort_wue(5.0_wp, 0.0_wp) == 0.0_wp, 'night-time WUE (E=0) -> 0, not a spike')
   end subroutine test_kernels

   !=======================================================================================!
   !  THE LOCKSTEP + FUSION CONTRACT of the fast diagnostic block. This is the highest-risk    !
   !  part of the whole subsystem: the block is filled by the fast loop and read at the output  !
   !  tick, with restructuring in between, so slot i must keep meaning cohort i across every     !
   !  permutation -- and a fused cohort's diagnostics must blend by each field's own kind.        !
   !                                                                                          !
   !  The three fields exercised below are deliberately one of EACH kind, so a blend that used    !
   !  one rule for all of them cannot pass: leaf_temp is INTENSIVE (leaf-area-weighted),           !
   !  sapflow is EXTENSIVE (nplant-weighted), abs_sw is GROUND-referenced (summed).                !
   !=======================================================================================!
   subroutine test_diag_lockstep()
      type(cohort_diag_block) :: d
      real(wp)    :: x(8)
      integer(ik) :: perm(4), n

      call cohort_diag_alloc(d, 8_ik, .true.)
      d%n = 4_ik
      d%v(CD_LEAF_TEMP, 1:4) = [290.0_wp, 300.0_wp, 310.0_wp, 320.0_wp]
      d%v(CD_SAPFLOW,   1:4) = [  1.0_wp,   2.0_wp,   3.0_wp,   4.0_wp]
      d%v(CD_ABS_SW,    1:4) = [ 10.0_wp,  20.0_wp,  30.0_wp,  40.0_wp]
      d%w(1:4)               = 1.0_wp

      !----- REORDER: one statement must permute EVERY field, not the ones someone remembered. !
      perm = [4_ik, 3_ik, 2_ik, 1_ik]
      call cohort_diag_reorder(d, perm, 4_ik)
      call cohort_diag_value(d, CD_LEAF_TEMP, x, n)
      call check(n == 4_ik, 'reorder keeps the live count')
      call check_close(x(1), 320.0_wp, 1.0e-12_wp, 'reorder permuted leaf_temp')
      call cohort_diag_value(d, CD_SAPFLOW, x, n)
      call check_close(x(1), 4.0_wp, 1.0e-12_wp, 'reorder permuted sapflow in lockstep')
      call cohort_diag_value(d, CD_ABS_SW, x, n)
      call check_close(x(1), 40.0_wp, 1.0e-12_wp, 'reorder permuted abs_sw in lockstep')

      !----- COPY SLOT (a cohort split): the daughter is a full copy. -----!
      call cohort_diag_copy_slot(d, 5_ik, 1_ik)
      call check_close(d%v(CD_SAPFLOW, 5_ik), d%v(CD_SAPFLOW, 1_ik), 1.0e-12_wp, 'copy_slot copies every field')
      call check_close(d%w(5_ik), d%w(1_ik), 1.0e-12_wp, 'copy_slot carries the dt weight')

      !----- CLEAR SLOT (a recruit into a reused, stale cull slot). -----!
      call cohort_diag_clear_slot(d, 5_ik)
      call check(d%v(CD_SAPFLOW, 5_ik) == 0.0_wp, 'clear_slot zeroes a reused slot')
      call check(d%w(5_ik) == 0.0_wp,             'clear_slot zeroes the weight')

      !----- FUSE slot 2 into slot 1 with leaf-area weights (3, 1) and nplant weights (2, 6).   !
      !      After the reorder above: slot1 = (320, 4, 40), slot2 = (310, 3, 30).                !
      !        INTENSIVE  leaf_temp = (3*320 + 1*310)/4      = 317.5                              !
      !        EXTENSIVE  sapflow   = (2*4   + 6*3  )/8      = 3.25                               !
      !        GROUND     abs_sw    =  40 + 30                = 70                                 !
      call cohort_diag_fuse(d, 1_ik, 2_ik, 3.0_wp, 1.0_wp, 2.0_wp, 6.0_wp, CDIAG_FUSE)
      call cohort_diag_value(d, CD_LEAF_TEMP, x, n)
      call check_close(x(1), 317.5_wp, 1.0e-12_wp, 'INTENSIVE field fuses leaf-area-weighted')
      call cohort_diag_value(d, CD_SAPFLOW, x, n)
      call check_close(x(1), 3.25_wp, 1.0e-12_wp, 'EXTENSIVE field fuses nplant-weighted')
      call cohort_diag_value(d, CD_ABS_SW, x, n)
      call check_close(x(1), 70.0_wp, 1.0e-12_wp, 'GROUND-referenced field fuses by SUM')
      !----- The three answers are mutually distinct, so a blend that applied one rule to all    !
      !      of them would fail at least two of the assertions above.  -----------------------!

      !----- An INACTIVE block is a total no-op: nothing allocated, every entry point safe. ---!
      block
         type(cohort_diag_block) :: dead
         call cohort_diag_alloc(dead, 8_ik, .false.)
         call cohort_diag_reorder(dead, perm, 4_ik)      ! must not touch anything
         call cohort_diag_value(dead, CD_SAPFLOW, x, n)
         call check(n == 0_ik, 'inactive block reports no data and does not crash')
      end block
   end subroutine test_diag_lockstep

end program test_diagnostic_reduce
