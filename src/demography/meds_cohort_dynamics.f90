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
!                         |dHite| < hgt_max*tol, unless the merged BA density exceeds the    !
!                         cap (hysteresis vs split). Stops when every patch has              !
!                         ccount <= |maxcohort|.                                            !
!   * fuse_2_cohorts -- merge donor into receptor: nplant summed, BA summed, per-plant       !
!                       fields nplant-weighted, DBH re-derived from BA; 1% BA check.         !
!   * split_cohorts -- split any cohort whose BA density exceeds the cap into two halves      !
!                      (halved nplant, +/-eps DBH), conserving nplant and BA.               !
!   * terminate_cohorts -- cull cohorts below the size/density floor.                       !
!==========================================================================================!
module meds_cohort_dynamics
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : pio4, tiny_num
   use meds_pft_params, only : dbh2h
   use meds_config,     only : meds_config_t
   use meds_demography_types,      only : community, cohort_compact, cohort_ensure_capacity,          &
                               copy_cohort_slot, rebuild_csr
   use meds_sort,       only : sort_cohorts
   implicit none
   private

   public :: new_fuse_cohorts, fuse_2_cohorts, split_cohorts, terminate_cohorts, max_ccount

contains

   !---------------------------------------------------------------------------------------!
   integer(ik) function max_ccount(comm)
      type(community), intent(in) :: comm
      if (comm%pat%n < 1_ik) then
         max_ccount = 0_ik
      else
         max_ccount = maxval(comm%pat%ccount(1:comm%pat%n))
      end if
   end function max_ccount

   !---------------------------------------------------------------------------------------!
   ! Cohort fusion with geometric tolerance relaxation.                                     !
   !---------------------------------------------------------------------------------------!
   subroutine new_fuse_cohorts(comm, cfg)
      type(community),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      real(wp)    :: tol
      integer(ik) :: ifus, maxc
      logical     :: force

      if (cfg%maxcohort == 0_ik) return
      force = (cfg%maxcohort < 0_ik)
      maxc  = abs(cfg%maxcohort)
      tol   = cfg%coh_size_tol_min

      do ifus = 1_ik, cfg%niter_cohfus
         call fuse_pass(comm, cfg, tol, force)
         if (.not. force) then
            if (max_ccount(comm) <= maxc) exit
         end if
         tol = tol * cfg%coh_size_tol_mult
      end do
   end subroutine new_fuse_cohorts

   !----- One fusion sweep over all patches at a fixed tolerance. --------------------------!
   subroutine fuse_pass(comm, cfg, tol, force)
      type(community),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: tol
      logical,             intent(in)    :: force
      logical, allocatable :: alive(:)
      integer(ik) :: ip, i0, i1, recc, donc
      real(wp)    :: diff_dbh, diff_hgt, hgt_max_r, ba_comb
      logical     :: did_fuse

      did_fuse = .false.
      allocate(alive(comm%coh%n))
      alive = .true.
      associate (c => comm%coh, p => comm%pat)
         do ip = 1_ik, p%n
            i0 = p%coff(ip)
            i1 = i0 + p%ccount(ip) - 1_ik
            do recc = i0, i1 - 1_ik
               if (.not. alive(recc)) cycle
               do donc = recc + 1_ik, i1
                  if (.not. alive(donc)) cycle
                  if (c%pft(donc) /= c%pft(recc)) cycle
                  ba_comb = c%nplant(recc) * c%basarea(recc) + c%nplant(donc) * c%basarea(donc)
                  if (.not. force .and. ba_comb >= cfg%ba_bin_cap) cycle
                  diff_dbh  = abs(c%dbh(donc)  - c%dbh(recc))
                  diff_hgt  = abs(c%hite(donc) - c%hite(recc))
                  hgt_max_r = c%p_hgt_min(recc) + c%p_b1ht(recc)
                  if (force .or. (diff_dbh < c%p_dbh_crit(recc) * tol .and.                 &
                                  diff_hgt < hgt_max_r * tol)) then
                     call fuse_2_cohorts(comm, recc, donc, cfg%cons_tol)
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
   subroutine fuse_2_cohorts(comm, recc, donc, cons_tol)
      type(community), intent(inout) :: comm
      integer(ik),     intent(in)    :: recc, donc
      real(wp),        intent(in)    :: cons_tol
      real(wp) :: nr, nd, ntot, ba_tot, rw, dw, ba_new

      associate (c => comm%coh)
         nr     = c%nplant(recc)
         nd     = c%nplant(donc)
         ntot   = nr + nd
         ba_tot = nr * c%basarea(recc) + nd * c%basarea(donc)     ! [cm2/m2] conserved
         rw     = nr / ntot
         dw     = nd / ntot
         !----- Per-plant fields: nplant-weighted means. ----------------------------------!
         c%monthly_dlnndt(recc) = rw * c%monthly_dlnndt(recc) + dw * c%monthly_dlnndt(donc)
         c%ddbh_dt(recc)        = rw * c%ddbh_dt(recc)        + dw * c%ddbh_dt(donc)
         !----- Conserve plant number and basal area; re-derive size. ---------------------!
         c%nplant(recc)  = ntot
         c%basarea(recc) = ba_tot / ntot                         ! [cm2/plant]
         c%dbh(recc)     = sqrt(c%basarea(recc) / pio4)
         c%hite(recc)    = dbh2h(c%p_hgt_min(recc), c%p_b1ht(recc), c%p_b2ht(recc), c%dbh(recc))
         !----- Conservation guard. -------------------------------------------------------!
         ba_new = c%nplant(recc) * c%basarea(recc)
         if (abs(ba_new - ba_tot) > cons_tol * max(ba_tot, tiny_num))                       &
            error stop 'fuse_2_cohorts: basal-area conservation violated'
      end associate
   end subroutine fuse_2_cohorts

   !---------------------------------------------------------------------------------------!
   ! Split every cohort whose basal-area density exceeds the cap, iterating until none do.   !
   !---------------------------------------------------------------------------------------!
   subroutine split_cohorts(comm, cfg)
      type(community),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      integer(ik) :: iter, i, n0, m, nsplit
      real(wp)    :: d0, eps, renorm, ba_before, ba_after

      if (.not. cfg%enable_cohort_fission) return
      eps = cfg%split_eps
      !----- Scale the two daughters' diameters so (1+-eps) splits conserve BA exactly:    !
      !      0.5*[(1+eps)^2+(1-eps)^2]/(1+eps^2) = 1.                                       !
      renorm = 1.0_wp / sqrt(1.0_wp + eps * eps)

      do iter = 1_ik, cfg%niter_cohfus
         n0 = comm%coh%n
         nsplit = 0_ik
         do i = 1_ik, n0
            if (comm%coh%nplant(i) * comm%coh%basarea(i) > cfg%ba_bin_cap) nsplit = nsplit + 1_ik
         end do
         if (nsplit == 0_ik) exit
         ba_before = sum(comm%coh%nplant(1:n0) * comm%coh%basarea(1:n0))
         call cohort_ensure_capacity(comm%coh, n0 + nsplit)
         m = n0
         associate (c => comm%coh)
            do i = 1_ik, n0
               if (c%nplant(i) * c%basarea(i) <= cfg%ba_bin_cap) cycle
               d0          = c%dbh(i)
               c%nplant(i) = 0.5_wp * c%nplant(i)
               !----- '+eps' half stays in slot i. ----------------------------------------!
               c%dbh(i)     = d0 * (1.0_wp + eps) * renorm
               c%basarea(i) = pio4 * c%dbh(i) * c%dbh(i)
               c%hite(i)    = dbh2h(c%p_hgt_min(i), c%p_b1ht(i), c%p_b2ht(i), c%dbh(i))
               !----- '-eps' half appended in slot m+1. -----------------------------------!
               m = m + 1_ik
               call copy_cohort_slot(c, m, i)              ! copies halved nplant + params
               c%dbh(m)     = d0 * (1.0_wp - eps) * renorm
               c%basarea(m) = pio4 * c%dbh(m) * c%dbh(m)
               c%hite(m)    = dbh2h(c%p_hgt_min(m), c%p_b1ht(m), c%p_b2ht(m), c%dbh(m))
            end do
            c%n = m
            ba_after = sum(c%nplant(1:m) * c%basarea(1:m))
         end associate
         if (abs(ba_after - ba_before) > cfg%cons_tol * max(ba_before, tiny_num))          &
            error stop 'split_cohorts: basal-area conservation violated'
         call rebuild_csr(comm)
         call sort_cohorts(comm)
      end do
   end subroutine split_cohorts

   !---------------------------------------------------------------------------------------!
   ! Cull cohorts below the basal-area-density floor or the absolute density floor.         !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_cohorts(comm, cfg)
      type(community),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      logical, allocatable :: keep(:)
      integer(ik)          :: i, n

      n = comm%coh%n
      if (n < 1_ik) return
      allocate(keep(n))
      do i = 1_ik, n
         keep(i) = (comm%coh%nplant(i) * comm%coh%basarea(i) >= cfg%min_cohort_size) .and.  &
                   (comm%coh%nplant(i) >= cfg%negligible_nplant)
      end do
      if (all(keep)) return
      call cohort_compact(comm%coh, keep)
      call rebuild_csr(comm)
   end subroutine terminate_cohorts

end module meds_cohort_dynamics
