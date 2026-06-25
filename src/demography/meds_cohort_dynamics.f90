!==========================================================================================!
! meds_cohort_dynamics -- cohort fusion, fission (split), and termination.                 !
!                                                                                          !
! Height/LAI analogue of ED2 fuse_fiss_utils. Cohorts are merged when they are similar in    !
! HEIGHT and split when their per-cohort LAI exceeds a cap; the conserved invariant on every  !
! merge/split is TOTAL ABOVEGROUND BIOMASS (carbon): fused diameters are re-derived from the   !
! conserved AGB via agb_to_dbh, never averaged.                                              !
!                                                                                          !
!   * new_fuse_cohorts -- geometric tolerance-relaxation loop; receptor (taller) absorbs    !
!                         donors of the same PFT when |dheight| < hgt_max*tol, unless the      !
!                         merged cohort LAI exceeds the cap (hysteresis vs split). Stops when  !
!                         every patch has cohort_count <= |max_cohort|.                        !
!   * fuse_2_cohorts -- merge donor into receptor: nplant summed, AGB summed, per-plant AGB   !
!                       averaged, DBH re-derived from AGB, geometry refreshed; 1% AGB check.   !
!   * split_cohorts -- split any cohort whose LAI density exceeds the cap into two halves      !
!                      (halved nplant, +/-eps DBH), conserving nplant and AGB.               !
!   * terminate_cohorts -- cull cohorts below the AGB/density floor.                          !
!==========================================================================================!
module meds_cohort_dynamics
   use meds_kinds,      only : wp, ik
   use meds_constants,  only : tiny_num
   use meds_allometry,  only : agb_to_dbh, agb_c2, b2Ht, hgt_max
   use meds_config,     only : meds_config_t
   use meds_demography_types,      only : site, cohort_compact, cohort_ensure_capacity,          &
                               copy_cohort_slot, rebuild_csr, set_cohort_size
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
      real(wp)    :: diff_height, lai_comb
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
                  lai_comb = c%nplant(recc) * c%larea(recc) + c%nplant(donc) * c%larea(donc)
                  if (.not. force .and. lai_comb >= cfg%cohort_lai_cap) cycle
                  diff_height = abs(c%height(donc) - c%height(recc))
                  if (force .or. diff_height < hgt_max * tol) then
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
   ! Merge donor cohort into receptor, conserving plant number and total aboveground         !
   ! biomass (carbon). The fused per-plant AGB is the nplant-weighted mean; the diameter is   !
   ! re-derived from it (agb_to_dbh) and the rest of the geometry refreshed.                  !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_cohorts(comm, recc, donc, conservation_tol)
      type(site), intent(inout) :: comm
      integer(ik),     intent(in)    :: recc, donc
      real(wp),        intent(in)    :: conservation_tol
      real(wp) :: nr, nd, ntot, agb_tot, agb_new

      associate (c => comm%coh)
         nr      = c%nplant(recc)
         nd      = c%nplant(donc)
         ntot    = nr + nd
         agb_tot = nr * c%agb(recc) + nd * c%agb(donc)     ! [kgC/m2] conserved
         !----- Conserve plant number and AGB; re-derive size from carbon. ----------------!
         c%nplant(recc) = ntot
         c%agb(recc)    = agb_tot / ntot                   ! [kgC/plant]
         c%dbh(recc)    = agb_to_dbh(c%agb(recc), c%p_wood_density(recc))
         call set_cohort_size(comm%coh, recc)              ! refresh height/BA/agb/larea
         !----- Conservation guard (AGB density). -----------------------------------------!
         agb_new = c%nplant(recc) * c%agb(recc)
         if (abs(agb_new - agb_tot) > conservation_tol * max(agb_tot, tiny_num))             &
            error stop 'fuse_2_cohorts: AGB conservation violated'
      end associate
   end subroutine fuse_2_cohorts

   !---------------------------------------------------------------------------------------!
   ! Split every cohort whose basal-area density exceeds the cap, iterating until none do.   !
   !---------------------------------------------------------------------------------------!
   subroutine split_cohorts(comm, cfg)
      type(site),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      integer(ik) :: iter, i, n0, m, nsplit
      real(wp)    :: d0, eps, renorm, p_agb, agb_before, agb_after

      if (.not. cfg%enable_cohort_fission) return
      eps = cfg%split_eps
      !----- AGB ~ dbh^p_agb (uncapped); scale the two daughters' diameters so the +/-eps   !
      !      split conserves AGB exactly: 0.5*[(1+eps)^p + (1-eps)^p] * renorm^p = 1.        !
      p_agb  = agb_c2 * (2.0_wp + b2Ht)
      renorm = (0.5_wp * ((1.0_wp + eps) ** p_agb + (1.0_wp - eps) ** p_agb)) ** (-1.0_wp / p_agb)

      do iter = 1_ik, cfg%n_cohort_fusion_iter
         n0 = comm%coh%n
         nsplit = 0_ik
         do i = 1_ik, n0
            if (comm%coh%nplant(i) * comm%coh%larea(i) > cfg%cohort_lai_cap) nsplit = nsplit + 1_ik
         end do
         if (nsplit == 0_ik) exit
         agb_before = sum(comm%coh%nplant(1:n0) * comm%coh%agb(1:n0))
         call cohort_ensure_capacity(comm%coh, n0 + nsplit)
         m = n0
         associate (c => comm%coh)
            do i = 1_ik, n0
               if (c%nplant(i) * c%larea(i) <= cfg%cohort_lai_cap) cycle
               d0          = c%dbh(i)
               c%nplant(i) = 0.5_wp * c%nplant(i)
               !----- '+eps' half stays in slot i. ----------------------------------------!
               c%dbh(i) = d0 * (1.0_wp + eps) * renorm
               call set_cohort_size(c, i)
               !----- '-eps' half appended in slot m+1. -----------------------------------!
               m = m + 1_ik
               call copy_cohort_slot(c, m, i)              ! copies halved nplant + params
               c%dbh(m) = d0 * (1.0_wp - eps) * renorm
               call set_cohort_size(c, m)
            end do
            c%n = m
            agb_after = sum(c%nplant(1:m) * c%agb(1:m))
         end associate
         if (abs(agb_after - agb_before) > cfg%conservation_tol * max(agb_before, tiny_num))       &
            error stop 'split_cohorts: AGB conservation violated'
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
         keep(i) = (comm%coh%nplant(i) * comm%coh%agb(i) >= cfg%min_cohort_agb) .and.        &
                   (comm%coh%nplant(i) >= cfg%negligible_nplant)
      end do
      if (all(keep)) return
      call cohort_compact(comm%coh, keep)
      call rebuild_csr(comm)
   end subroutine terminate_cohorts

end module meds_cohort_dynamics
