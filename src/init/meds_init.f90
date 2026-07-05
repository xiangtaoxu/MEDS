!==========================================================================================!
! meds_init -- helpers to build the initial community a run starts from.                    !
!                                                                                          !
! Two entry points exist: a near-bare-ground site (for spin-up demos and the test suite),  !
! and a restart from a cohort CENSUS file (init_from_census) -- a CSV produced by a previous !
! MEDS run or a field inventory. add_cohort/finalize_init are the low-level builders both    !
! the demo/tests and the census reader share.                                              !
!==========================================================================================!
module meds_init
   use meds_kinds,      only : wp, ik
   use meds_config,     only : meds_config_t, DIST_PRIMARY, growth_window_steps
   use meds_demography_types,      only : site_t, site_alloc, cohort_ensure_capacity, rebuild_csr,  &
                                          set_cohort_size, assign_cohort_id, assign_patch_id,        &
                                          GROWTH_AVG_UNSET
   use meds_demography_fusefiss, only : sort_cohorts
   implicit none
   private

   public :: init_bare_ground, add_cohort, finalize_init, init_from_census

contains

   !---------------------------------------------------------------------------------------!
   ! Near-bare-ground site_t: n_patch identical empty patches sharing the site_t area.     !
   !---------------------------------------------------------------------------------------!
   subroutine init_bare_ground(site, cfg, n_patch)
      type(site_t),        intent(out) :: site
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: n_patch
      integer(ik) :: ip
      call site_alloc(site, cfg%pft%n, coh_cap = 64_ik, pat_cap = max(n_patch, 1_ik),          &
                      n_growth_window = growth_window_steps(cfg))
      site%patch%n = n_patch
      do ip = 1_ik, n_patch
         site%patch%area(ip)      = 1.0_wp / real(n_patch, wp)
         site%patch%age(ip)       = 0.0_wp
         site%patch%dist_type(ip) = DIST_PRIMARY
         call assign_patch_id(site, ip)
      end do
      site%patch%recruit_pool = 0.0_wp
      site%cohort%n = 0_ik
      call rebuild_csr(site)
   end subroutine init_bare_ground

   !---------------------------------------------------------------------------------------!
   ! Append one cohort to patch ip (caller invokes finalize_init afterwards).               !
   !---------------------------------------------------------------------------------------!
   subroutine add_cohort(site, cfg, ip, ipft, nplant, dbh)
      type(site_t),     intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      integer(ik),         intent(in)    :: ip, ipft
      real(wp),            intent(in)    :: nplant, dbh
      integer(ik) :: m
      call cohort_ensure_capacity(site%cohort, site%cohort%n + 1_ik)
      m = site%cohort%n + 1_ik
      associate (cohort => site%cohort, pft => cfg%pft)
         cohort%pft(m)            = ipft
         cohort%owner_patch(m)    = ip
         cohort%nplant(m)         = nplant
         cohort%dbh(m)            = dbh
         cohort%growth_avg(m)     = GROWTH_AVG_UNSET    ! set on its first growth step
         cohort%growth_accum(m)   = 0.0_wp
         cohort%growth_count(m)   = 0_ik
         cohort%p_dbh_critical(m) = pft%dbh_critical(ipft)
         cohort%p_wood_density(m) = pft%wood_density(ipft)
         cohort%p_hgt_max(m)      = pft%hgt_max(ipft)
         cohort%p_sla(m)                = pft%sla(ipft)
         cohort%p_aboveground_frac(m)   = pft%aboveground_frac(ipft)
         cohort%p_root_to_leaf_ratio(m) = pft%root_to_leaf_ratio(ipft)
         cohort%p_storage_cushion(m)    = pft%storage_cushion(ipft)
         call set_cohort_size(cohort, m)            ! height/basal_area/agb/leaf_area + carbon pools from dbh
      end associate
      call assign_cohort_id(site, m)
      site%cohort%n = m
   end subroutine add_cohort

   subroutine finalize_init(site)
      type(site_t), intent(inout) :: site
      call rebuild_csr(site)
      call sort_cohorts(site)
   end subroutine finalize_init

   !---------------------------------------------------------------------------------------!
   ! Initialize a site from a cohort CENSUS file: a CSV with a header line and one row per   !
   ! cohort, columns:                                                                        !
   !                                                                                          !
   !     site_id, patch_id, cohort_id, dbh, height, pft, nplant                              !
   !                                                                                          !
   ! Such a file is produced by a previous MEDS run or a field inventory. One site is built:  !
   ! rows are filtered to a single site_id (`use_site_id` if present, else the first site_id  !
   ! encountered); each distinct patch_id then becomes an equal-area DIST_PRIMARY patch (the   !
   ! census carries no area or age), and each row becomes a cohort. `dbh` drives the           !
   ! allometry (height/basal_area/agb/leaf_area via set_cohort_size, the model's invariant),   !
   ! so the file's `height` and `cohort_id` are read for provenance / sanity-checking only.    !
   !                                                                                          !
   ! Blank lines and lines beginning with '#' are ignored; the first remaining line is taken  !
   ! as the column header and skipped. found=.false. (the site is left unbuilt) if the file    !
   ! cannot be opened or holds no row for the selected site -- the caller can then fall back    !
   ! to bare ground. The file is read twice (count patches, then add cohorts) so no row buffer  !
   ! is needed; cohort capacity grows on demand.                                              !
   !---------------------------------------------------------------------------------------!
   subroutine init_from_census(site, cfg, path, found, use_site_id)
      type(site_t),        intent(out) :: site
      type(meds_config_t), intent(in)  :: cfg
      character(len=*),    intent(in)  :: path
      logical,             intent(out) :: found
      integer(ik), intent(in), optional :: use_site_id

      integer,     parameter :: LINELEN = 1024
      character(len=LINELEN) :: line
      integer(ik)            :: sid, pid, cid, ipft, target, npatch, ip
      real(wp)               :: dbh, height, nplant
      integer(ik), allocatable :: patch_id_of(:)
      logical                :: ok_row, seen_header, target_known
      integer                :: u, ios

      found = .false.
      open(newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return

      !----- Pass 1: pick the target site and collect its distinct patch ids (first-seen). -!
      target_known = present(use_site_id)
      if (target_known) target = use_site_id
      allocate(patch_id_of(0))
      seen_header = .false.
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (.not. is_data_line(line)) cycle
         if (.not. seen_header) then ; seen_header = .true. ; cycle ; end if   ! column header
         call parse_census_row(line, sid, pid, cid, dbh, height, ipft, nplant, ok_row)
         if (.not. ok_row) error stop 'meds_init: malformed census row in '//trim(path)
         if (.not. target_known) then ; target = sid ; target_known = .true. ; end if
         if (sid /= target) cycle
         if (.not. any(patch_id_of == pid)) patch_id_of = [patch_id_of, pid]
      end do

      npatch = int(size(patch_id_of), ik)
      if (npatch == 0_ik) then ; close(u) ; return ; end if   ! no rows for the selected site

      !----- Build the patches (equal area, primary). --------------------------------------!
      call init_bare_ground(site, cfg, npatch)

      !----- Pass 2: add every cohort of the target site to its patch. ---------------------!
      rewind(u)
      seen_header = .false.
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (.not. is_data_line(line)) cycle
         if (.not. seen_header) then ; seen_header = .true. ; cycle ; end if
         call parse_census_row(line, sid, pid, cid, dbh, height, ipft, nplant, ok_row)
         if (.not. ok_row) cycle
         if (sid /= target) cycle
         if (ipft < 1_ik .or. ipft > cfg%pft%n)                                              &
            error stop 'meds_init: census pft index out of range'
         if (dbh <= 0.0_wp .or. height <= 0.0_wp .or. nplant <= 0.0_wp)                      &
            error stop 'meds_init: census row with non-positive dbh/height/nplant'
         ip = int(findloc(patch_id_of, pid, dim=1), ik)
         call add_cohort(site, cfg, ip, ipft, nplant, dbh)
      end do
      close(u)

      call finalize_init(site)
      found = .true.
   end subroutine init_from_census

   !----- A line carries data unless it is blank or a '#' comment. ------------------------!
   logical function is_data_line(line)
      character(len=*), intent(in) :: line
      character(len=len(line))     :: s
      s = adjustl(line)
      is_data_line = (len_trim(s) > 0) .and. (s(1:1) /= '#')
   end function is_data_line

   !----- Parse one comma-separated census row (commas -> blanks, then list-directed read). !
   subroutine parse_census_row(line, sid, pid, cid, dbh, height, ipft, nplant, ok)
      character(len=*), intent(in)  :: line
      integer(ik),      intent(out) :: sid, pid, cid, ipft
      real(wp),         intent(out) :: dbh, height, nplant
      logical,          intent(out) :: ok
      character(len=len(line)) :: buf
      integer :: k, ios
      buf = line
      do k = 1, len(buf)
         if (buf(k:k) == ',') buf(k:k) = ' '
      end do
      read(buf, *, iostat=ios) sid, pid, cid, dbh, height, ipft, nplant
      ok = (ios == 0)
   end subroutine parse_census_row

end module meds_init
