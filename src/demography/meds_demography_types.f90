!==========================================================================================!
! meds_demography_types -- the demographic state, as a flat site_t-wide Structure-of-Arrays.            !
!                                                                                          !
! ALL cohorts of the WHOLE site_t live in one contiguous set of 1-D arrays (`cohort_block`),  !
! so the dominant daily kernels are a single unit-stride sweep, ideal for SIMD and GPU.     !
! Patch membership is a CSR map (`cohort_offset`,`cohort_count`) over the flat cohort arrays. Cohorts of   !
! a patch occupy a contiguous slice, so per-patch operations (sort, fusion) work on slices.  !
!                                                                                          !
! Every structural change goes through ONE centralized routine, `cohort_reorder`, which     !
! permutes EVERY per-cohort array in lockstep -- this is the single place to update when a   !
! field is added, eliminating the ED2 "forgot to reallocate an array" class of bug.         !
!==========================================================================================!
module meds_demography_types
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : pio4
   use meds_pft_params, only : pft_table_t
   use meds_allometry,  only : dbh_to_height, dbh_to_agb, dbh_to_leaf_area
   implicit none
   private

   public :: cohort_block, patch_index, site_t
   public :: site_alloc, site_free
   public :: cohort_ensure_capacity, cohort_reorder, cohort_compact, gather_pft_params
   public :: patch_ensure_capacity, rebuild_csr, copy_cohort_slot, set_cohort_size
   public :: assign_cohort_id, assign_patch_id
   public :: GROWTH_AVG_UNSET

   !----- Sentinel for a not-yet-sampled moving-average growth: a freshly created cohort holds  !
   !      this until its first growth step adds a sample (negative => "unset"; a real growth     !
   !      average is >= 0). Used by the growth kernel and the rate evaluator.                    !
   real(wp), parameter :: GROWTH_AVG_UNSET = -1.0_wp

   !---------------------------------------------------------------------------------------!
   ! All cohorts of the site_t (contiguous SoA). `dbh` is the prognostic size axis; `height`,   !
   ! `basal_area`, `agb` and `leaf_area` are derived but stored so the hot kernels stay         !
   ! arithmetic-only. `agb` (per-plant carbon) is the quantity conserved across fusion/        !
   ! fission; `leaf_area` (per-plant) gives the cohort's LAI contribution (nplant*leaf_area).   !
   !---------------------------------------------------------------------------------------!
   type :: cohort_block
      integer(ik) :: n   = 0_ik       !< cohorts in use
      integer(ik) :: cap = 0_ik       !< allocated capacity
      !----- Prognostic / primary state. --------------------------------------------------!
      integer(ik), allocatable :: pft(:)            !< PFT index
      real(wp),    allocatable :: nplant(:)         !< [plant/m2] PRIMARY demographic state
      real(wp),    allocatable :: dbh(:)            !< [cm]       prognostic size
      real(wp),    allocatable :: height(:)           !< [m]        = dbh_to_height(dbh), cached
      real(wp),    allocatable :: basal_area(:)        !< [cm2/plant]= pio4*dbh^2, cached
      real(wp),    allocatable :: agb(:)             !< [kgC/plant] conserved carbon, cached
      real(wp),    allocatable :: leaf_area(:)           !< [m2/plant]  leaf area, cached (LAI=nplant*leaf_area)
      !----- Sliding simple-moving-average growth rate [cm/yr] driving growth-dependent          !
      !      mortality. growth_hist is the per-cohort ring buffer of the last growth_window steps !
      !      (column per cohort; the write slot is the site-global growth_hist_pos), growth_accum !
      !      its running sum and growth_count its fill level (<= window); growth_avg =            !
      !      accum/count is the SMA (GROWTH_AVG_UNSET until the first sample). A PROGNOSTIC        !
      !      history (not derivable from dbh).                                                    !
      real(wp),    allocatable :: growth_avg(:)
      real(wp),    allocatable :: growth_accum(:)   !< [cm/yr] running sum of the samples in the window
      integer(ik), allocatable :: growth_count(:)   !< samples currently in the window (<= growth_window)
      real(wp),    allocatable :: growth_hist(:,:)  !< (growth_window, cohort) ring buffer of growth samples
      !----- Per-cohort gathered PFT params (kernels never index the PFT table). ----------!
      real(wp),    allocatable :: p_dbh_critical(:)
      real(wp),    allocatable :: p_wood_density(:)
      !----- Host-only back-index used to regroup the flat array by patch. ----------------!
      integer(ik), allocatable :: owner_patch(:)
      !----- Persistent identity: a global id stamped at creation and carried (in lockstep   !
      !      with every other field) through sorts/fusion/compaction, so an external tracker  !
      !      can follow one cohort across output records until it fuses away or is culled.     !
      integer(ik), allocatable :: global_id(:)
   end type cohort_block

   !---------------------------------------------------------------------------------------!
   ! Patches of the site_t (SoA) with a CSR map into the cohort block.                        !
   !---------------------------------------------------------------------------------------!
   type :: patch_index
      integer(ik) :: n   = 0_ik
      integer(ik) :: cap = 0_ik
      real(wp),    allocatable :: area(:)           !< fraction, sum = 1 (kept real64)
      real(wp),    allocatable :: age(:)            !< [yr]
      integer(ik), allocatable :: dist_type(:)      !< disturbance/land-use class
      integer(ik), allocatable :: cohort_offset(:)           !< CSR: first cohort of patch (1-based)
      integer(ik), allocatable :: cohort_count(:)         !< cohorts in patch
      real(wp),    allocatable :: recruit_pool(:,:)  !< (pft, patch) carry-forward density
      !----- Persistent identity (as for cohorts): stamped at creation, carried through       !
      !      patch sort/fusion/compaction so a tracker can follow one patch across records.    !
      integer(ik), allocatable :: global_id(:)
   end type patch_index

   type :: site_t
      type(cohort_block) :: cohort
      type(patch_index)  :: patch
      real(wp)           :: site_area = 1.0_wp
      integer(ik)        :: n_pft     = 0_ik
      !----- Monotonic counters for the persistent global ids (never reused within a run). ---!
      integer(ik)        :: next_cohort_id = 1_ik
      integer(ik)        :: next_patch_id  = 1_ik
      !----- Global ring-buffer write slot for the moving-average growth (same for all cohorts !
      !      since they all take one sample per step); advanced by update_demography.           !
      integer(ik)        :: growth_hist_pos = 0_ik
   end type site_t

contains

   !=======================================================================================!
   !  Allocation                                                                            !
   !=======================================================================================!
   subroutine site_alloc(site, n_pft, coh_cap, pat_cap, n_growth_window)
      type(site_t), intent(out) :: site
      integer(ik),     intent(in)  :: n_pft, coh_cap, pat_cap, n_growth_window
      site%n_pft = n_pft
      call cohort_alloc(site%cohort, max(coh_cap, 1_ik), max(n_growth_window, 1_ik))
      call patch_alloc(site%patch, max(pat_cap, 1_ik), n_pft)
   end subroutine site_alloc

   subroutine site_free(site)
      type(site_t), intent(inout) :: site
      site%cohort%n = 0_ik ; site%cohort%cap = 0_ik
      site%patch%n = 0_ik ; site%patch%cap = 0_ik
      if (allocated(site%cohort%nplant)) deallocate(site%cohort%pft, site%cohort%nplant, site%cohort%dbh, &
         site%cohort%height, site%cohort%basal_area, site%cohort%agb, site%cohort%leaf_area,                       &
         site%cohort%growth_avg, site%cohort%growth_accum, site%cohort%growth_count,                     &
         site%cohort%growth_hist, site%cohort%p_dbh_critical, site%cohort%p_wood_density,                &
         site%cohort%owner_patch, site%cohort%global_id)
      if (allocated(site%patch%area)) deallocate(site%patch%area, site%patch%age, site%patch%dist_type, &
         site%patch%cohort_offset, site%patch%cohort_count, site%patch%recruit_pool, site%patch%global_id)
   end subroutine site_free

   subroutine cohort_alloc(cohort, cap, nwin)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: cap, nwin
      cohort%cap = cap ; cohort%n = 0_ik
      allocate(cohort%pft(cap), cohort%owner_patch(cap), cohort%global_id(cap))
      allocate(cohort%nplant(cap), cohort%dbh(cap), cohort%height(cap), cohort%basal_area(cap),            &
               cohort%agb(cap), cohort%leaf_area(cap), cohort%growth_avg(cap),                            &
               cohort%growth_accum(cap), cohort%growth_count(cap), cohort%growth_hist(nwin, cap))
      allocate(cohort%p_dbh_critical(cap), cohort%p_wood_density(cap))
      cohort%pft = 0_ik ; cohort%owner_patch = 0_ik ; cohort%global_id = 0_ik
      cohort%nplant = 0.0_wp ; cohort%dbh = 0.0_wp ; cohort%height = 0.0_wp ; cohort%basal_area = 0.0_wp
      cohort%agb = 0.0_wp ; cohort%leaf_area = 0.0_wp ; cohort%growth_avg = GROWTH_AVG_UNSET
      cohort%growth_accum = 0.0_wp ; cohort%growth_count = 0_ik ; cohort%growth_hist = 0.0_wp
      cohort%p_dbh_critical = 0.0_wp ; cohort%p_wood_density = 0.0_wp
   end subroutine cohort_alloc

   subroutine patch_alloc(patch, cap, n_pft)
      type(patch_index), intent(inout) :: patch
      integer(ik),       intent(in)    :: cap, n_pft
      patch%cap = cap ; patch%n = 0_ik
      allocate(patch%area(cap), patch%age(cap), patch%dist_type(cap), patch%global_id(cap))
      allocate(patch%cohort_offset(cap), patch%cohort_count(cap), patch%recruit_pool(n_pft, cap))
      patch%area = 0.0_wp ; patch%age = 0.0_wp ; patch%dist_type = 1_ik ; patch%global_id = 0_ik
      patch%cohort_offset = 0_ik ; patch%cohort_count = 0_ik ; patch%recruit_pool = 0.0_wp
   end subroutine patch_alloc

   !=======================================================================================!
   !  Capacity growth (1.5x, copy old into new).                                            !
   !=======================================================================================!
   subroutine cohort_ensure_capacity(cohort, need)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: need
      type(cohort_block)                :: tmp
      integer(ik)                       :: newcap, m
      if (need <= cohort%cap) return
      newcap = max(need, (cohort%cap * 3_ik) / 2_ik + 1_ik)
      call cohort_alloc(tmp, newcap, int(size(cohort%growth_hist, 1), ik))  ! keep the window size
      m = cohort%n
      tmp%n = m
      tmp%pft(1:m)            = cohort%pft(1:m)
      tmp%owner_patch(1:m)    = cohort%owner_patch(1:m)
      tmp%nplant(1:m)         = cohort%nplant(1:m)
      tmp%dbh(1:m)            = cohort%dbh(1:m)
      tmp%height(1:m)           = cohort%height(1:m)
      tmp%basal_area(1:m)        = cohort%basal_area(1:m)
      tmp%agb(1:m)            = cohort%agb(1:m)
      tmp%leaf_area(1:m)          = cohort%leaf_area(1:m)
      tmp%growth_avg(1:m)     = cohort%growth_avg(1:m)
      tmp%growth_accum(1:m)   = cohort%growth_accum(1:m)
      tmp%growth_count(1:m)   = cohort%growth_count(1:m)
      tmp%growth_hist(:,1:m)  = cohort%growth_hist(:,1:m)
      tmp%p_dbh_critical(1:m)     = cohort%p_dbh_critical(1:m)
      tmp%p_wood_density(1:m) = cohort%p_wood_density(1:m)
      tmp%global_id(1:m)      = cohort%global_id(1:m)
      call move_alloc_block(tmp, cohort)
   end subroutine cohort_ensure_capacity

   !----- Move every allocatable from src into dst (src is reset). ------------------------!
   subroutine move_alloc_block(src, dst)
      type(cohort_block), intent(inout) :: src
      type(cohort_block), intent(inout) :: dst
      dst%n = src%n ; dst%cap = src%cap
      call move_alloc(src%pft, dst%pft)
      call move_alloc(src%owner_patch, dst%owner_patch)
      call move_alloc(src%nplant, dst%nplant)
      call move_alloc(src%dbh, dst%dbh)
      call move_alloc(src%height, dst%height)
      call move_alloc(src%basal_area, dst%basal_area)
      call move_alloc(src%agb, dst%agb)
      call move_alloc(src%leaf_area, dst%leaf_area)
      call move_alloc(src%growth_avg, dst%growth_avg)
      call move_alloc(src%growth_accum, dst%growth_accum)
      call move_alloc(src%growth_count, dst%growth_count)
      call move_alloc(src%growth_hist, dst%growth_hist)
      call move_alloc(src%p_dbh_critical, dst%p_dbh_critical)
      call move_alloc(src%p_wood_density, dst%p_wood_density)
      call move_alloc(src%global_id, dst%global_id)
   end subroutine move_alloc_block

   subroutine patch_ensure_capacity(patch, need, n_pft)
      type(patch_index), intent(inout) :: patch
      integer(ik),       intent(in)    :: need, n_pft
      type(patch_index)                :: tmp
      integer(ik)                       :: newcap, m
      if (need <= patch%cap) return
      newcap = max(need, (patch%cap * 3_ik) / 2_ik + 1_ik)
      call patch_alloc(tmp, newcap, n_pft)
      m = patch%n ; tmp%n = m
      tmp%area(1:m)           = patch%area(1:m)
      tmp%age(1:m)            = patch%age(1:m)
      tmp%dist_type(1:m)      = patch%dist_type(1:m)
      tmp%cohort_offset(1:m)           = patch%cohort_offset(1:m)
      tmp%cohort_count(1:m)         = patch%cohort_count(1:m)
      tmp%recruit_pool(:,1:m) = patch%recruit_pool(:,1:m)
      tmp%global_id(1:m)      = patch%global_id(1:m)
      patch%n = tmp%n ; patch%cap = tmp%cap
      call move_alloc(tmp%area, patch%area)             ; call move_alloc(tmp%age, patch%age)
      call move_alloc(tmp%dist_type, patch%dist_type)
      call move_alloc(tmp%cohort_offset, patch%cohort_offset)             ; call move_alloc(tmp%cohort_count, patch%cohort_count)
      call move_alloc(tmp%recruit_pool, patch%recruit_pool)
      call move_alloc(tmp%global_id, patch%global_id)
   end subroutine patch_ensure_capacity

   !=======================================================================================!
   !  Centralized lockstep reorder: place old index perm(k) at new position k, for k=1..m.  !
   !  Array assignment with a vector subscript on the RHS is evaluated to a temporary       !
   !  first, so self-permutation is safe. THIS is the one place that knows every field.     !
   !=======================================================================================!
   subroutine cohort_reorder(cohort, perm, m)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: perm(:)
      integer(ik),        intent(in)    :: m
      cohort%pft(1:m)            = cohort%pft(perm(1:m))
      cohort%owner_patch(1:m)    = cohort%owner_patch(perm(1:m))
      cohort%nplant(1:m)         = cohort%nplant(perm(1:m))
      cohort%dbh(1:m)            = cohort%dbh(perm(1:m))
      cohort%height(1:m)           = cohort%height(perm(1:m))
      cohort%basal_area(1:m)        = cohort%basal_area(perm(1:m))
      cohort%agb(1:m)            = cohort%agb(perm(1:m))
      cohort%leaf_area(1:m)          = cohort%leaf_area(perm(1:m))
      cohort%growth_avg(1:m)     = cohort%growth_avg(perm(1:m))
      cohort%growth_accum(1:m)   = cohort%growth_accum(perm(1:m))
      cohort%growth_count(1:m)   = cohort%growth_count(perm(1:m))
      cohort%growth_hist(:,1:m)  = cohort%growth_hist(:,perm(1:m))
      cohort%p_dbh_critical(1:m)     = cohort%p_dbh_critical(perm(1:m))
      cohort%p_wood_density(1:m) = cohort%p_wood_density(perm(1:m))
      cohort%global_id(1:m)      = cohort%global_id(perm(1:m))
      cohort%n = m
   end subroutine cohort_reorder

   !----- Compact by keep-mask (length n): drop entries where keep is .false. -------------!
   subroutine cohort_compact(cohort, keep)
      type(cohort_block), intent(inout) :: cohort
      logical,            intent(in)    :: keep(:)
      integer(ik), allocatable :: perm(:)
      integer(ik)              :: i, m
      allocate(perm(cohort%n)) ; m = 0_ik
      do i = 1_ik, cohort%n
         if (keep(i)) then
            m = m + 1_ik
            perm(m) = i
         end if
      end do
      call cohort_reorder(cohort, perm, m)
   end subroutine cohort_compact

   !----- Copy every per-cohort field from slot src to slot dst (centralized). ------------!
   subroutine copy_cohort_slot(cohort, dst, src)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: dst, src
      cohort%pft(dst)            = cohort%pft(src)
      cohort%owner_patch(dst)    = cohort%owner_patch(src)
      cohort%nplant(dst)         = cohort%nplant(src)
      cohort%dbh(dst)            = cohort%dbh(src)
      cohort%height(dst)           = cohort%height(src)
      cohort%basal_area(dst)        = cohort%basal_area(src)
      cohort%agb(dst)            = cohort%agb(src)
      cohort%leaf_area(dst)          = cohort%leaf_area(src)
      cohort%growth_avg(dst)     = cohort%growth_avg(src)
      cohort%growth_accum(dst)   = cohort%growth_accum(src)
      cohort%growth_count(dst)   = cohort%growth_count(src)
      cohort%growth_hist(:,dst)  = cohort%growth_hist(:,src)
      cohort%p_dbh_critical(dst)     = cohort%p_dbh_critical(src)
      cohort%p_wood_density(dst) = cohort%p_wood_density(src)
      cohort%global_id(dst)      = cohort%global_id(src)
   end subroutine copy_cohort_slot

   !----- Fill the gathered per-cohort PFT params from the trait table. -------------------!
   subroutine gather_pft_params(cohort, pft)
      type(cohort_block), intent(inout) :: cohort
      type(pft_table_t),  intent(in)    :: pft
      integer(ik) :: i, p
      do i = 1_ik, cohort%n
         p = cohort%pft(i)
         cohort%p_dbh_critical(i)     = pft%dbh_critical(p)
         cohort%p_wood_density(i) = pft%wood_density(p)
      end do
   end subroutine gather_pft_params

   !---------------------------------------------------------------------------------------!
   ! Re-derive the cached geometry (height, basal area, AGB, leaf area) of ONE cohort slot   !
   ! from its prognostic diameter and gathered wood density. Single-slot host analogue of    !
   ! the array math in growth_step; used by recruitment, setup, fusion and fission so the     !
   ! allometry lives in exactly one (shared) place.                                          !
   !---------------------------------------------------------------------------------------!
   subroutine set_cohort_size(cohort, i)
      type(cohort_block), intent(inout) :: cohort
      integer(ik),        intent(in)    :: i
      cohort%height(i)      = dbh_to_height(cohort%dbh(i))
      cohort%basal_area(i)  = pio4 * cohort%dbh(i) * cohort%dbh(i)
      cohort%agb(i)         = dbh_to_agb(cohort%dbh(i), cohort%height(i), cohort%p_wood_density(i))
      cohort%leaf_area(i)       = dbh_to_leaf_area(cohort%dbh(i), cohort%height(i))
   end subroutine set_cohort_size

   !=======================================================================================!
   !  Rebuild the CSR patch map by stably regrouping cohorts by owner_patch (counting sort).!
   !=======================================================================================!
   subroutine rebuild_csr(site)
      type(site_t), intent(inout) :: site
      integer(ik), allocatable :: perm(:), pos(:), slot(:)
      integer(ik)              :: i, ip, np, nc
      np = site%patch%n ; nc = site%cohort%n
      site%patch%cohort_count(1:np) = 0_ik
      do i = 1_ik, nc
         ip = site%cohort%owner_patch(i)
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
      call cohort_reorder(site%cohort, perm, nc)
   end subroutine rebuild_csr

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

end module meds_demography_types
