!==========================================================================================!
! meds_cohort_dynamics -- cohort fusion, fission (split), and termination.                 !
!                                                                                          !
! Diameter/size-distribution analogue of ED2 fuse_fiss_utils, with LAI/leaf-area criteria   !
! replaced by DBH + height similarity and a per-cohort basal-area density cap. The          !
! conserved invariant on every merge/split is TOTAL BASAL AREA (the diameter analogue of    !
! ED2's biomass): fused diameters are re-derived from conserved BA, never averaged.         !
!                                                                                          !
!   * new_fuse_cohorts -- geometric tolerance-relaxation loop; receptor (taller) absorbs    !
!                         donors of the same PFT when |dDBH| < dbh_crit*tol and              !
!                         |dheight| < height_max*tol, unless the merged BA density exceeds the    !
!                         cap (hysteresis vs split). Stops when every patch has              !
!                         cohort_count <= |max_cohort|.                                            !
!   * fuse_2_cohorts -- merge donor into receptor: nplant summed, BA summed, per-plant       !
!                       fields nplant-weighted, DBH re-derived from BA; 1% BA check.         !
!   * split_cohorts -- split any cohort whose BA density exceeds the cap into two halves      !
!                      (halved nplant, +/-eps DBH), conserving nplant and BA.               !
!   * terminate_cohorts -- cull cohorts below the size/density floor.                       !
!==========================================================================================!
module meds_cohort_dynamics
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : pio4, tiny_num
   use meds_pft_params, only : dbh_to_height
   use meds_config,     only : meds_config_t
   use meds_demography_types,      only : site, cohort_compact, cohort_ensure_capacity,          &
                               copy_cohort_slot, rebuild_csr
   use meds_sort,       only : sort_cohorts
   implicit none
   private

   public :: new_fuse_cohorts, fuse_2_cohorts, split_cohorts, terminate_cohorts, max_cohort_count

contains

   !---------------------------------------------------------------------------------------!
   integer(ik) function max_cohort_count(comm)
      type(site), intent(in) :: comm
      if (comm%pat%n < 1_ik) then
         max_cohort_count = 0_ik
      else
         max_cohort_count = maxval(comm%pat%cohort_count(1:comm%pat%n))
      end if
   end function max_cohort_count

   !---------------------------------------------------------------------------------------!
   ! Cohort fusion with geometric tolerance relaxation.                                     !
   !---------------------------------------------------------------------------------------!
   subroutine new_fuse_cohorts(comm, cfg)
      type(site),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      real(wp)    :: tol
      integer(ik) :: ifus, maxc
      logical     :: force

      if (cfg%max_cohort == 0_ik) return
      force = (cfg%max_cohort < 0_ik)
      maxc  = abs(cfg%max_cohort)
      tol   = cfg%cohort_size_tol_min

      do ifus = 1_ik, cfg%n_cohort_fusion_iter
         call fuse_pass(comm, cfg, tol, force)
         if (.not. force) then
            if (max_cohort_count(comm) <= maxc) exit
         end if
         tol = tol * cfg%cohort_size_tol_mult
      end do
   end subroutine new_fuse_cohorts

   !----- One fusion sweep over all patches at a fixed tolerance. --------------------------!
   subroutine fuse_pass(comm, cfg, tol, force)
      type(site),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: tol
      logical,             intent(in)    :: force
      logical, allocatable :: alive(:)
      integer(ik) :: ip, i0, i1, recc, donc
      real(wp)    :: diff_dbh, diff_height, height_max_r, ba_comb
      logical     :: did_fuse

      did_fuse = .false.
      allocate(alive(comm%coh%n))
      alive = .true.
      associate (c => comm%coh, p => comm%pat)
         do ip = 1_ik, p%n
            i0 = p%cohort_offset(ip)
            i1 = i0 + p%cohort_count(ip) - 1_ik
            do recc = i0, i1 - 1_ik
               if (.not. alive(recc)) cycle
               do donc = recc + 1_ik, i1
                  if (.not. alive(donc)) cycle
                  if (c%pft(donc) /= c%pft(recc)) cycle
                  ba_comb = c%nplant(recc) * c%basal_area(recc) + c%nplant(donc) * c%basal_area(donc)
                  if (.not. force .and. ba_comb >= cfg%basal_area_bin_cap) cycle
                  diff_dbh  = abs(c%dbh(donc)  - c%dbh(recc))
                  diff_height  = abs(c%height(donc) - c%height(recc))
                  height_max_r = c%p_height_min(recc) + c%p_b1_height(recc)
                  if (force .or. (diff_dbh < c%p_dbh_crit(recc) * tol .and.                 &
                                  diff_height < height_max_r * tol)) then
                     call fuse_2_cohorts(comm, recc, donc, cfg%conservation_tol)
                     alive(donc) = .false.
                     did_fuse    = .true.
                  end if
               end do
            end do
         end do
      end associate
      if (did_fuse) then
         call cohort_compact(comm%coh, alive)
         call rebuild_csr(comm)
         call sort_cohorts(comm)
      end if
   end subroutine fuse_pass

   !---------------------------------------------------------------------------------------!
   ! Merge donor cohort into receptor, conserving plant number and total basal area.        !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_cohorts(comm, recc, donc, conservation_tol)
      type(site), intent(inout) :: comm
      integer(ik),     intent(in)    :: recc, donc
      real(wp),        intent(in)    :: conservation_tol
      real(wp) :: nr, nd, ntot, ba_tot, ba_new

      associate (c => comm%coh)
         nr     = c%nplant(recc)
         nd     = c%nplant(donc)
         ntot   = nr + nd
         ba_tot = nr * c%basal_area(recc) + nd * c%basal_area(donc)     ! [cm2/m2] conserved
         !----- Conserve plant number and basal area; re-derive size. ---------------------!
         c%nplant(recc)  = ntot
         c%basal_area(recc) = ba_tot / ntot                         ! [cm2/plant]
         c%dbh(recc)     = sqrt(c%basal_area(recc) / pio4)
         c%height(recc)    = dbh_to_height(c%p_height_min(recc), c%p_b1_height(recc), c%p_b2_height(recc), c%dbh(recc))
         !----- Conservation guard. -------------------------------------------------------!
         ba_new = c%nplant(recc) * c%basal_area(recc)
         if (abs(ba_new - ba_tot) > conservation_tol * max(ba_tot, tiny_num))                       &
            error stop 'fuse_2_cohorts: basal-area conservation violated'
      end associate
   end subroutine fuse_2_cohorts

   !---------------------------------------------------------------------------------------!
   ! Split every cohort whose basal-area density exceeds the cap, iterating until none do.   !
   !---------------------------------------------------------------------------------------!
   subroutine split_cohorts(comm, cfg)
      type(site),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      integer(ik) :: iter, i, n0, m, nsplit
      real(wp)    :: d0, eps, renorm, ba_before, ba_after

      if (.not. cfg%enable_cohort_fission) return
      eps = cfg%split_eps
      !----- Scale the two daughters' diameters so (1+-eps) splits conserve BA exactly:    !
      !      0.5*[(1+eps)^2+(1-eps)^2]/(1+eps^2) = 1.                                       !
      renorm = 1.0_wp / sqrt(1.0_wp + eps * eps)

      do iter = 1_ik, cfg%n_cohort_fusion_iter
         n0 = comm%coh%n
         nsplit = 0_ik
         do i = 1_ik, n0
            if (comm%coh%nplant(i) * comm%coh%basal_area(i) > cfg%basal_area_bin_cap) nsplit = nsplit + 1_ik
         end do
         if (nsplit == 0_ik) exit
         ba_before = sum(comm%coh%nplant(1:n0) * comm%coh%basal_area(1:n0))
         call cohort_ensure_capacity(comm%coh, n0 + nsplit)
         m = n0
         associate (c => comm%coh)
            do i = 1_ik, n0
               if (c%nplant(i) * c%basal_area(i) <= cfg%basal_area_bin_cap) cycle
               d0          = c%dbh(i)
               c%nplant(i) = 0.5_wp * c%nplant(i)
               !----- '+eps' half stays in slot i. ----------------------------------------!
               c%dbh(i)     = d0 * (1.0_wp + eps) * renorm
               c%basal_area(i) = pio4 * c%dbh(i) * c%dbh(i)
               c%height(i)    = dbh_to_height(c%p_height_min(i), c%p_b1_height(i), c%p_b2_height(i), c%dbh(i))
               !----- '-eps' half appended in slot m+1. -----------------------------------!
               m = m + 1_ik
               call copy_cohort_slot(c, m, i)              ! copies halved nplant + params
               c%dbh(m)     = d0 * (1.0_wp - eps) * renorm
               c%basal_area(m) = pio4 * c%dbh(m) * c%dbh(m)
               c%height(m)    = dbh_to_height(c%p_height_min(m), c%p_b1_height(m), c%p_b2_height(m), c%dbh(m))
            end do
            c%n = m
            ba_after = sum(c%nplant(1:m) * c%basal_area(1:m))
         end associate
         if (abs(ba_after - ba_before) > cfg%conservation_tol * max(ba_before, tiny_num))          &
            error stop 'split_cohorts: basal-area conservation violated'
         call rebuild_csr(comm)
         call sort_cohorts(comm)
      end do
   end subroutine split_cohorts

   !---------------------------------------------------------------------------------------!
   ! Cull cohorts below the basal-area-density floor or the absolute density floor.         !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_cohorts(comm, cfg)
      type(site),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      logical, allocatable :: keep(:)
      integer(ik)          :: i, n

      n = comm%coh%n
      if (n < 1_ik) return
      allocate(keep(n))
      do i = 1_ik, n
         keep(i) = (comm%coh%nplant(i) * comm%coh%basal_area(i) >= cfg%min_cohort_size) .and.  &
                   (comm%coh%nplant(i) >= cfg%negligible_nplant)
      end do
      if (all(keep)) return
      call cohort_compact(comm%coh, keep)
      call rebuild_csr(comm)
   end subroutine terminate_cohorts

end module meds_cohort_dynamics
