!==========================================================================================!
! meds_demography_types -- the demographic state, as a flat site-wide Structure-of-Arrays.            !
!                                                                                          !
! ALL cohorts of the WHOLE site live in one contiguous set of 1-D arrays (`cohort_block`),  !
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

   public :: cohort_block, patch_index, site
   public :: site_alloc, site_free
   public :: cohort_ensure_capacity, cohort_reorder, cohort_compact, gather_pft_params
   public :: patch_ensure_capacity, rebuild_csr, copy_cohort_slot, set_cohort_size

   !---------------------------------------------------------------------------------------!
   ! All cohorts of the site (contiguous SoA). `dbh` is the prognostic size axis; `height`,   !
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
      !----- Per-cohort gathered PFT params (kernels never index the PFT table). ----------!
      real(wp),    allocatable :: p_dbh_crit(:)
      real(wp),    allocatable :: p_wood_density(:)
      !----- Host-only back-index used to regroup the flat array by patch. ----------------!
      integer(ik), allocatable :: owner_patch(:)
   end type cohort_block

   !---------------------------------------------------------------------------------------!
   ! Patches of the site (SoA) with a CSR map into the cohort block.                        !
   !---------------------------------------------------------------------------------------!
   type :: patch_index
      integer(ik) :: n   = 0_ik
      integer(ik) :: cap = 0_ik
      real(wp),    allocatable :: area(:)           !< fraction, sum = 1 (kept real64)
      real(wp),    allocatable :: age(:)            !< [yr]
      integer(ik), allocatable :: dist_type(:)      !< disturbance/land-use class
      real(wp),    allocatable :: avg_daily_temp(:) !< [K] frost driver
      real(wp),    allocatable :: min_month_temp(:) !< [K] recruit eligibility
      integer(ik), allocatable :: cohort_offset(:)           !< CSR: first cohort of patch (1-based)
      integer(ik), allocatable :: cohort_count(:)         !< cohorts in patch
      real(wp),    allocatable :: recruit_pool(:,:)  !< (pft, patch) carry-forward density
   end type patch_index

   type :: site
      type(cohort_block) :: coh
      type(patch_index)  :: pat
      real(wp)           :: site_area = 1.0_wp
      integer(ik)        :: n_pft     = 0_ik
   end type site

contains

   !=======================================================================================!
   !  Allocation                                                                            !
   !=======================================================================================!
   subroutine site_alloc(comm, n_pft, coh_cap, pat_cap)
      type(site), intent(out) :: comm
      integer(ik),     intent(in)  :: n_pft, coh_cap, pat_cap
      comm%n_pft = n_pft
      call cohort_alloc(comm%coh, max(coh_cap, 1_ik))
      call patch_alloc(comm%pat, max(pat_cap, 1_ik), n_pft)
   end subroutine site_alloc

   subroutine site_free(comm)
      type(site), intent(inout) :: comm
      comm%coh%n = 0_ik ; comm%coh%cap = 0_ik
      comm%pat%n = 0_ik ; comm%pat%cap = 0_ik
      if (allocated(comm%coh%nplant)) deallocate(comm%coh%pft, comm%coh%nplant, comm%coh%dbh, &
         comm%coh%height, comm%coh%basal_area, comm%coh%agb, comm%coh%leaf_area,                       &
         comm%coh%p_dbh_crit, comm%coh%p_wood_density, comm%coh%owner_patch)
      if (allocated(comm%pat%area)) deallocate(comm%pat%area, comm%pat%age, comm%pat%dist_type, &
         comm%pat%avg_daily_temp, comm%pat%min_month_temp, comm%pat%cohort_offset, comm%pat%cohort_count,    &
         comm%pat%recruit_pool)
   end subroutine site_free

   subroutine cohort_alloc(coh, cap)
      type(cohort_block), intent(inout) :: coh
      integer(ik),        intent(in)    :: cap
      coh%cap = cap ; coh%n = 0_ik
      allocate(coh%pft(cap), coh%owner_patch(cap))
      allocate(coh%nplant(cap), coh%dbh(cap), coh%height(cap), coh%basal_area(cap),            &
               coh%agb(cap), coh%leaf_area(cap))
      allocate(coh%p_dbh_crit(cap), coh%p_wood_density(cap))
      coh%pft = 0_ik ; coh%owner_patch = 0_ik
      coh%nplant = 0.0_wp ; coh%dbh = 0.0_wp ; coh%height = 0.0_wp ; coh%basal_area = 0.0_wp
      coh%agb = 0.0_wp ; coh%leaf_area = 0.0_wp
      coh%p_dbh_crit = 0.0_wp ; coh%p_wood_density = 0.0_wp
   end subroutine cohort_alloc

   subroutine patch_alloc(pat, cap, n_pft)
      type(patch_index), intent(inout) :: pat
      integer(ik),       intent(in)    :: cap, n_pft
      pat%cap = cap ; pat%n = 0_ik
      allocate(pat%area(cap), pat%age(cap), pat%dist_type(cap))
      allocate(pat%avg_daily_temp(cap), pat%min_month_temp(cap))
      allocate(pat%cohort_offset(cap), pat%cohort_count(cap), pat%recruit_pool(n_pft, cap))
      pat%area = 0.0_wp ; pat%age = 0.0_wp ; pat%dist_type = 1_ik
      pat%avg_daily_temp = 0.0_wp ; pat%min_month_temp = 0.0_wp
      pat%cohort_offset = 0_ik ; pat%cohort_count = 0_ik ; pat%recruit_pool = 0.0_wp
   end subroutine patch_alloc

   !=======================================================================================!
   !  Capacity growth (1.5x, copy old into new).                                            !
   !=======================================================================================!
   subroutine cohort_ensure_capacity(coh, need)
      type(cohort_block), intent(inout) :: coh
      integer(ik),        intent(in)    :: need
      type(cohort_block)                :: tmp
      integer(ik)                       :: newcap, m
      if (need <= coh%cap) return
      newcap = max(need, (coh%cap * 3_ik) / 2_ik + 1_ik)
      call cohort_alloc(tmp, newcap)
      m = coh%n
      tmp%n = m
      tmp%pft(1:m)            = coh%pft(1:m)
      tmp%owner_patch(1:m)    = coh%owner_patch(1:m)
      tmp%nplant(1:m)         = coh%nplant(1:m)
      tmp%dbh(1:m)            = coh%dbh(1:m)
      tmp%height(1:m)           = coh%height(1:m)
      tmp%basal_area(1:m)        = coh%basal_area(1:m)
      tmp%agb(1:m)            = coh%agb(1:m)
      tmp%leaf_area(1:m)          = coh%leaf_area(1:m)
      tmp%p_dbh_crit(1:m)     = coh%p_dbh_crit(1:m)
      tmp%p_wood_density(1:m) = coh%p_wood_density(1:m)
      call move_alloc_block(tmp, coh)
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
      call move_alloc(src%p_dbh_crit, dst%p_dbh_crit)
      call move_alloc(src%p_wood_density, dst%p_wood_density)
   end subroutine move_alloc_block

   subroutine patch_ensure_capacity(pat, need, n_pft)
      type(patch_index), intent(inout) :: pat
      integer(ik),       intent(in)    :: need, n_pft
      type(patch_index)                :: tmp
      integer(ik)                       :: newcap, m
      if (need <= pat%cap) return
      newcap = max(need, (pat%cap * 3_ik) / 2_ik + 1_ik)
      call patch_alloc(tmp, newcap, n_pft)
      m = pat%n ; tmp%n = m
      tmp%area(1:m)           = pat%area(1:m)
      tmp%age(1:m)            = pat%age(1:m)
      tmp%dist_type(1:m)      = pat%dist_type(1:m)
      tmp%avg_daily_temp(1:m) = pat%avg_daily_temp(1:m)
      tmp%min_month_temp(1:m) = pat%min_month_temp(1:m)
      tmp%cohort_offset(1:m)           = pat%cohort_offset(1:m)
      tmp%cohort_count(1:m)         = pat%cohort_count(1:m)
      tmp%recruit_pool(:,1:m) = pat%recruit_pool(:,1:m)
      pat%n = tmp%n ; pat%cap = tmp%cap
      call move_alloc(tmp%area, pat%area)             ; call move_alloc(tmp%age, pat%age)
      call move_alloc(tmp%dist_type, pat%dist_type)
      call move_alloc(tmp%avg_daily_temp, pat%avg_daily_temp)
      call move_alloc(tmp%min_month_temp, pat%min_month_temp)
      call move_alloc(tmp%cohort_offset, pat%cohort_offset)             ; call move_alloc(tmp%cohort_count, pat%cohort_count)
      call move_alloc(tmp%recruit_pool, pat%recruit_pool)
   end subroutine patch_ensure_capacity

   !=======================================================================================!
   !  Centralized lockstep reorder: place old index perm(k) at new position k, for k=1..m.  !
   !  Array assignment with a vector subscript on the RHS is evaluated to a temporary       !
   !  first, so self-permutation is safe. THIS is the one place that knows every field.     !
   !=======================================================================================!
   subroutine cohort_reorder(coh, perm, m)
      type(cohort_block), intent(inout) :: coh
      integer(ik),        intent(in)    :: perm(:)
      integer(ik),        intent(in)    :: m
      coh%pft(1:m)            = coh%pft(perm(1:m))
      coh%owner_patch(1:m)    = coh%owner_patch(perm(1:m))
      coh%nplant(1:m)         = coh%nplant(perm(1:m))
      coh%dbh(1:m)            = coh%dbh(perm(1:m))
      coh%height(1:m)           = coh%height(perm(1:m))
      coh%basal_area(1:m)        = coh%basal_area(perm(1:m))
      coh%agb(1:m)            = coh%agb(perm(1:m))
      coh%leaf_area(1:m)          = coh%leaf_area(perm(1:m))
      coh%p_dbh_crit(1:m)     = coh%p_dbh_crit(perm(1:m))
      coh%p_wood_density(1:m) = coh%p_wood_density(perm(1:m))
      coh%n = m
   end subroutine cohort_reorder

   !----- Compact by keep-mask (length n): drop entries where keep is .false. -------------!
   subroutine cohort_compact(coh, keep)
      type(cohort_block), intent(inout) :: coh
      logical,            intent(in)    :: keep(:)
      integer(ik), allocatable :: perm(:)
      integer(ik)              :: i, m
      allocate(perm(coh%n)) ; m = 0_ik
      do i = 1_ik, coh%n
         if (keep(i)) then
            m = m + 1_ik
            perm(m) = i
         end if
      end do
      call cohort_reorder(coh, perm, m)
   end subroutine cohort_compact

   !----- Copy every per-cohort field from slot src to slot dst (centralized). ------------!
   subroutine copy_cohort_slot(coh, dst, src)
      type(cohort_block), intent(inout) :: coh
      integer(ik),        intent(in)    :: dst, src
      coh%pft(dst)            = coh%pft(src)
      coh%owner_patch(dst)    = coh%owner_patch(src)
      coh%nplant(dst)         = coh%nplant(src)
      coh%dbh(dst)            = coh%dbh(src)
      coh%height(dst)           = coh%height(src)
      coh%basal_area(dst)        = coh%basal_area(src)
      coh%agb(dst)            = coh%agb(src)
      coh%leaf_area(dst)          = coh%leaf_area(src)
      coh%p_dbh_crit(dst)     = coh%p_dbh_crit(src)
      coh%p_wood_density(dst) = coh%p_wood_density(src)
   end subroutine copy_cohort_slot

   !----- Fill the gathered per-cohort PFT params from the trait table. -------------------!
   subroutine gather_pft_params(coh, pft)
      type(cohort_block), intent(inout) :: coh
      type(pft_table_t),  intent(in)    :: pft
      integer(ik) :: i, p
      do i = 1_ik, coh%n
         p = coh%pft(i)
         coh%p_dbh_crit(i)     = pft%dbh_crit(p)
         coh%p_wood_density(i) = pft%wood_density(p)
      end do
   end subroutine gather_pft_params

   !---------------------------------------------------------------------------------------!
   ! Re-derive the cached geometry (height, basal area, AGB, leaf area) of ONE cohort slot   !
   ! from its prognostic diameter and gathered wood density. Single-slot host analogue of    !
   ! the array math in growth_step; used by recruitment, setup, fusion and fission so the     !
   ! allometry lives in exactly one (shared) place.                                          !
   !---------------------------------------------------------------------------------------!
   subroutine set_cohort_size(coh, i)
      type(cohort_block), intent(inout) :: coh
      integer(ik),        intent(in)    :: i
      coh%height(i)      = dbh_to_height(coh%dbh(i))
      coh%basal_area(i)  = pio4 * coh%dbh(i) * coh%dbh(i)
      coh%agb(i)         = dbh_to_agb(coh%dbh(i), coh%height(i), coh%p_wood_density(i))
      coh%leaf_area(i)       = dbh_to_leaf_area(coh%dbh(i), coh%height(i))
   end subroutine set_cohort_size

   !=======================================================================================!
   !  Rebuild the CSR patch map by stably regrouping cohorts by owner_patch (counting sort).!
   !=======================================================================================!
   subroutine rebuild_csr(comm)
      type(site), intent(inout) :: comm
      integer(ik), allocatable :: perm(:), pos(:), slot(:)
      integer(ik)              :: i, ip, np, nc
      np = comm%pat%n ; nc = comm%coh%n
      comm%pat%cohort_count(1:np) = 0_ik
      do i = 1_ik, nc
         ip = comm%coh%owner_patch(i)
         comm%pat%cohort_count(ip) = comm%pat%cohort_count(ip) + 1_ik
      end do
      comm%pat%cohort_offset(1) = 1_ik
      do ip = 2_ik, np
         comm%pat%cohort_offset(ip) = comm%pat%cohort_offset(ip-1) + comm%pat%cohort_count(ip-1)
      end do
      !----- Stable placement: preserve within-patch order. -------------------------------!
      allocate(pos(np), slot(nc), perm(nc))
      pos(1:np) = comm%pat%cohort_offset(1:np)
      do i = 1_ik, nc
         ip = comm%coh%owner_patch(i)
         slot(i) = pos(ip)
         pos(ip) = pos(ip) + 1_ik
      end do
      do i = 1_ik, nc
         perm(slot(i)) = i
      end do
      call cohort_reorder(comm%coh, perm, nc)
   end subroutine rebuild_csr

end module meds_demography_types
