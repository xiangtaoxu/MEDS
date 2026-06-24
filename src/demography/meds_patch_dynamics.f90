!==========================================================================================!
! meds_patch_dynamics -- patch size-distribution profiles, patch fusion, and termination.   !
!                                                                                          !
! ED2's patch fusion compares patches by their cumulative-LAI light profile; MEDS replaces  !
! that with the cumulative BASAL-AREA size distribution (top-down by DBH class), normalized  !
! by the site reference so the avg/max profile distances behave like ED2's [0,1] light       !
! differences. Two patches fuse when avg distance <= pat_prof_tol AND max distance <=         !
! pat_prof_tol*pat_prof_mxd_fac (and they share a disturbance type), with geometric          !
! tolerance relaxation until npatch <= |maxpatch|.                                           !
!                                                                                          !
! fuse_2_patches conserves SITE-LEVEL plant number: cohort densities are rescaled by area    !
! fractions when two patches of areas Ar, Ad combine into Ar+Ad, so Ar*dens_r is preserved.  !
!==========================================================================================!
module meds_patch_dynamics
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num, area_tol
   use meds_config,    only : meds_config_t
   use meds_demography_types,     only : site, cohort_compact, rebuild_csr
   use meds_sort,      only : sort_cohorts
   implicit none
   private

   public :: new_fuse_patches, fuse_2_patches, terminate_patches, patch_profile

   real(wp), parameter :: cm2_to_m2  = 1.0e-4_wp
   real(wp), parameter :: prof_relax = 1.5_wp     !< per-iteration tolerance growth

contains

   !---------------------------------------------------------------------------------------!
   ! DBH bin index for a diameter: number of interior edges below it, plus one (1..ndbh).   !
   !---------------------------------------------------------------------------------------!
   pure integer(ik) function bin_index(dbh, edges) result(ib)
      real(wp), intent(in) :: dbh, edges(:)
      integer(ik)          :: k
      ib = 1_ik
      do k = 1_ik, size(edges, kind=ik)
         if (dbh >= edges(k)) ib = ib + 1_ik
      end do
   end function bin_index

   !---------------------------------------------------------------------------------------!
   ! Cumulative (top-down) basal-area profile [m2 m-2] for one patch, by DBH class.         !
   !---------------------------------------------------------------------------------------!
   subroutine patch_profile(comm, cfg, ip, ba_cum)
      type(site),     intent(in)  :: comm
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: ip
      real(wp),            intent(out) :: ba_cum(:)        ! size cfg%ff_ndbh
      integer(ik) :: i, i0, i1, ib
      ba_cum = 0.0_wp
      i0 = comm%pat%coff(ip)
      i1 = i0 + comm%pat%ccount(ip) - 1_ik
      do i = i0, i1
         ib = bin_index(comm%coh%dbh(i), cfg%dbh_edges)
         ba_cum(ib) = ba_cum(ib) + comm%coh%nplant(i) * comm%coh%basarea(i) * cm2_to_m2
      end do
      !----- Accumulate from the largest class downward. ----------------------------------!
      do ib = cfg%ff_ndbh - 1_ik, 1_ik, -1_ik
         ba_cum(ib) = ba_cum(ib) + ba_cum(ib + 1_ik)
      end do
   end subroutine patch_profile

   !---------------------------------------------------------------------------------------!
   ! Patch fusion with geometric tolerance relaxation.                                      !
   !---------------------------------------------------------------------------------------!
   subroutine new_fuse_patches(comm, cfg)
      type(site),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      real(wp)    :: tol
      integer(ik) :: it, maxp
      logical     :: force

      if (cfg%maxpatch == 0_ik) return
      force = (cfg%maxpatch < 0_ik)
      maxp  = abs(cfg%maxpatch)
      tol   = cfg%pat_prof_tol

      do it = 1_ik, cfg%niter_patfus
         call patch_fuse_pass(comm, cfg, tol, force)
         if (.not. force .and. comm%pat%n <= maxp) exit
         tol = tol * prof_relax
      end do
   end subroutine new_fuse_patches

   !----- One patch-fusion sweep at a fixed tolerance. ------------------------------------!
   subroutine patch_fuse_pass(comm, cfg, tol, force)
      type(site),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      real(wp),            intent(in)    :: tol
      logical,             intent(in)    :: force
      logical,     allocatable :: alive(:)
      real(wp),    allocatable :: prof(:,:), pr(:), pd(:)
      real(wp)    :: ba_ref, davg, dmax, dif
      integer(ik) :: np, recp, donp, ib, npop
      logical     :: similar, both_empty

      np = comm%pat%n
      if (np < 2_ik) return
      allocate(alive(np), prof(cfg%ff_ndbh, np), pr(cfg%ff_ndbh), pd(cfg%ff_ndbh))
      alive = .true.
      ba_ref = tiny_num
      do recp = 1_ik, np
         call patch_profile(comm, cfg, recp, prof(:, recp))
         ba_ref = max(ba_ref, prof(1, recp))
      end do

      do recp = 1_ik, np - 1_ik
         if (.not. alive(recp)) cycle
         do donp = recp + 1_ik, np
            if (.not. alive(donp)) cycle
            if (comm%pat%dist_type(donp) /= comm%pat%dist_type(recp)) cycle
            both_empty = (comm%pat%ccount(recp) == 0_ik .and. comm%pat%ccount(donp) == 0_ik)
            if (both_empty) then
               similar = .true.
            else
               pr = prof(:, recp) ; pd = prof(:, donp)
               davg = 0.0_wp ; dmax = 0.0_wp ; npop = 0_ik
               do ib = 1_ik, cfg%ff_ndbh
                  if (pr(ib) <= tiny_num .and. pd(ib) <= tiny_num) cycle  ! top-down early skip
                  dif  = abs(pr(ib) - pd(ib)) / ba_ref
                  davg = davg + dif
                  dmax = max(dmax, dif)
                  npop = npop + 1_ik
               end do
               if (npop > 0_ik) davg = davg / real(npop, wp)
               similar = (davg <= tol) .and. (dmax <= tol * cfg%pat_prof_mxd_fac)
            end if
            if (force .or. similar) then
               call fuse_2_patches(comm, recp, donp)
               call rebuild_csr(comm)                 ! recp slice now holds merged cohorts
               call patch_profile(comm, cfg, recp, prof(:, recp))
               alive(donp) = .false.
            end if
         end do
      end do

      call patch_compact(comm, alive)
      call sort_cohorts(comm)
   end subroutine patch_fuse_pass

   !---------------------------------------------------------------------------------------!
   ! Merge donor patch into receptor, conserving site-level plant number via area weights.   !
   !---------------------------------------------------------------------------------------!
   subroutine fuse_2_patches(comm, recp, donp)
      type(site), intent(inout) :: comm
      integer(ik),     intent(in)    :: recp, donp
      real(wp)    :: ar, ad, anew, rawgt, dawgt
      integer(ik) :: i, i0, i1

      associate (p => comm%pat, c => comm%coh)
         ar = p%area(recp) ; ad = p%area(donp) ; anew = ar + ad
         if (anew <= tiny_num) return
         rawgt = ar / anew ; dawgt = ad / anew
         !----- Area-weighted patch scalars. ----------------------------------------------!
         p%age(recp)            = rawgt * p%age(recp)            + dawgt * p%age(donp)
         p%avg_daily_temp(recp) = rawgt * p%avg_daily_temp(recp) + dawgt * p%avg_daily_temp(donp)
         p%min_month_temp(recp) = min(p%min_month_temp(recp), p%min_month_temp(donp))
         p%recruit_pool(:,recp) = rawgt * p%recruit_pool(:,recp) + dawgt * p%recruit_pool(:,donp)
         !----- Rescale receptor cohort densities (slice currently holds all recp cohorts). !
         i0 = p%coff(recp) ; i1 = i0 + p%ccount(recp) - 1_ik
         do i = i0, i1
            c%nplant(i) = c%nplant(i) * rawgt
         end do
         !----- Rescale and reassign donor cohorts to the receptor. -----------------------!
         i0 = p%coff(donp) ; i1 = i0 + p%ccount(donp) - 1_ik
         do i = i0, i1
            c%nplant(i)      = c%nplant(i) * dawgt
            c%owner_patch(i) = recp
         end do
         p%area(recp) = anew
         p%area(donp) = 0.0_wp
      end associate
   end subroutine fuse_2_patches

   !---------------------------------------------------------------------------------------!
   ! Remove patches with area below the floor; renormalize remaining areas to sum to 1.     !
   !---------------------------------------------------------------------------------------!
   subroutine terminate_patches(comm, cfg)
      type(site),     intent(inout) :: comm
      type(meds_config_t), intent(in)    :: cfg
      logical, allocatable :: pkeep(:)
      integer(ik)          :: np, ip, imax
      real(wp)             :: s

      np = comm%pat%n
      if (np < 1_ik) return
      allocate(pkeep(np))
      pkeep = comm%pat%area(1:np) >= cfg%min_patch_area
      if (.not. any(pkeep)) then
         imax = maxloc(comm%pat%area(1:np), dim=1)
         pkeep(imax) = .true.
      end if
      if (.not. all(pkeep)) call patch_compact(comm, pkeep)

      !----- Renormalize areas to sum to 1. -----------------------------------------------!
      s = sum(comm%pat%area(1:comm%pat%n))
      if (s > tiny_num) then
         do ip = 1_ik, comm%pat%n
            comm%pat%area(ip) = comm%pat%area(ip) / s
         end do
      end if
   end subroutine terminate_patches

   !---------------------------------------------------------------------------------------!
   ! Compact the patch list by keep-mask: drop dropped patches AND their cohorts, remap      !
   ! owner_patch to the new patch indices, and rebuild the CSR map.                          !
   !---------------------------------------------------------------------------------------!
   subroutine patch_compact(comm, pkeep)
      type(site), intent(inout) :: comm
      logical,         intent(in)    :: pkeep(:)
      integer(ik), allocatable :: newidx(:)
      logical,     allocatable :: ckeep(:)
      integer(ik) :: np, ip, k, i, nc

      np = comm%pat%n
      allocate(newidx(np))
      k = 0_ik
      do ip = 1_ik, np
         if (pkeep(ip)) then
            k = k + 1_ik
            newidx(ip) = k
         else
            newidx(ip) = 0_ik
         end if
      end do

      !----- Cohorts: keep only those whose patch survives; remap owner. ------------------!
      nc = comm%coh%n
      allocate(ckeep(nc))
      do i = 1_ik, nc
         ckeep(i) = pkeep(comm%coh%owner_patch(i))
         if (ckeep(i)) comm%coh%owner_patch(i) = newidx(comm%coh%owner_patch(i))
      end do
      call cohort_compact(comm%coh, ckeep)

      !----- Patches: compact every patch array by pkeep (RHS pack -> safe temp). ---------!
      associate (p => comm%pat)
         p%area(1:k)           = pack(p%area(1:np),           pkeep)
         p%age(1:k)            = pack(p%age(1:np),            pkeep)
         p%dist_type(1:k)      = pack(p%dist_type(1:np),      pkeep)
         p%avg_daily_temp(1:k) = pack(p%avg_daily_temp(1:np), pkeep)
         p%min_month_temp(1:k) = pack(p%min_month_temp(1:np), pkeep)
         block
            integer(ik) :: jp
            do jp = 1_ik, comm%n_pft
               p%recruit_pool(jp, 1:k) = pack(p%recruit_pool(jp, 1:np), pkeep)
            end do
         end block
         p%n = k
      end associate
      call rebuild_csr(comm)
   end subroutine patch_compact

end module meds_patch_dynamics
