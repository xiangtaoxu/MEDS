! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_init -- helpers to build the initial community a run starts from.                    !
!                                                                                          !
! Two entry points exist: a near-bare-ground site (for spin-up demos and the test suite),  !
! and a start from a cohort CENSUS file (init_from_census) -- a CSV produced by a previous   !
! MEDS run or a field inventory. add_cohort/finalize_init are the low-level builders both    !
! the demo/tests and the census reader share.                                              !
!==========================================================================================!
module meds_init
   use meds_kinds,      only : wp, ik
   use meds_config,     only : meds_config_t, DIST_PRIMARY, growth_window_steps
   use meds_site_state_types,      only : site_t, site_alloc, cohort_ensure_capacity, rebuild_csr,  &
                                          assign_cohort_id, assign_patch_id, init_cohort
   use meds_demography_cohort_fusefiss, only : sort_cohorts
   implicit none
   private

   public :: init_bare_ground, add_cohort, finalize_init, init_from_census

   !----- The census file's column layout: each index is the column's position, 0 if absent. -!
   integer, parameter :: CENSUS_MAX_COLS = 32, CENSUS_FIELD_LEN = 64
   type :: census_columns_t
      logical :: header = .false.
      integer :: n = 0
      integer :: site = 0, patch = 0, cohort = 0, dbh = 0, height = 0, pft = 0, nplant = 0, area = 0
   end type census_columns_t

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
      call init_cohort(site%cohort, m, cfg%pft, ipft, ip, nplant, dbh)
      call assign_cohort_id(site, m)
      site%cohort%n = m
   end subroutine add_cohort

   subroutine finalize_init(site)
      type(site_t), intent(inout) :: site
      call rebuild_csr(site)
      call sort_cohorts(site)
   end subroutine finalize_init

   !---------------------------------------------------------------------------------------!
   ! Initialize a site from a cohort CENSUS file: a CSV with one row per cohort, produced by   !
   ! a previous MEDS run or from a field inventory. One site is built.                        !
   !                                                                                          !
   ! COLUMNS ARE MATCHED BY NAME from the header line, in any order:                           !
   !   required  patch_id, dbh [cm], pft, nplant [plants per m2 of the patch]                  !
   !   optional  patch_area [any unit], site_id, cohort_id, height [m]                         !
   ! Each distinct patch_id is one primary patch of age 0. With patch_area the patches take     !
   ! their areas normalized to fractions of the site; without it they share the site equally.  !
   ! Rows are filtered to one site_id (`use_site_id` if present, else the first one seen); a    !
   ! file without the column is one site. `dbh` drives the allometry (set_cohort_size), so      !
   ! `height` and `cohort_id` are provenance: height, when present, must be positive.          !
   !                                                                                          !
   ! A file whose first data line is seven numbers is read positionally as                      !
   ! `site_id, patch_id, cohort_id, dbh, height, pft, nplant`, which is how census files        !
   ! without a header were written.                                                            !
   !                                                                                          !
   ! Blank lines and lines beginning with '#' are ignored. found=.false. (the site is left      !
   ! unbuilt) if the file cannot be opened or holds no row for the selected site, so the caller !
   ! can fall back to bare ground. An unknown or duplicated column, a missing required one, a  !
   ! malformed row, a non-positive value, or a patch whose rows disagree on patch_area stops    !
   ! the run. The file is read twice (patches, then cohorts), so no row buffer is needed.      !
   !---------------------------------------------------------------------------------------!
   subroutine init_from_census(site, cfg, path, found, use_site_id)
      type(site_t),        intent(out) :: site
      type(meds_config_t), intent(in)  :: cfg
      character(len=*),    intent(in)  :: path
      logical,             intent(out) :: found
      integer(ik), intent(in), optional :: use_site_id

      integer,     parameter :: LINELEN = 1024
      character(len=LINELEN) :: line
      type(census_columns_t) :: cols
      integer(ik)            :: sid, pid, ipft, sel_site, npatch, ip
      real(wp)               :: dbh, height, nplant, area
      integer(ik), allocatable :: patch_id_of(:)
      real(wp),    allocatable :: patch_area_of(:)
      logical                :: target_known, have_layout
      integer                :: u, ios, ipass
      character(len=24)      :: pid_str

      found = .false.
      open(newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return

      target_known = present(use_site_id)
      if (target_known) sel_site = use_site_id
      allocate(patch_id_of(0), patch_area_of(0))
      npatch = 0_ik

      do ipass = 1, 2
         if (ipass == 2) then
            npatch = int(size(patch_id_of), ik)
            if (npatch == 0_ik) exit                         ! no rows for the selected site
            !----- The patches: equal areas, or the file's areas normalized to the site. --------!
            call init_bare_ground(site, cfg, npatch)
            if (cols%area > 0) site%patch%area(1:npatch) = patch_area_of / sum(patch_area_of)
            rewind(u)
         end if
         have_layout = .false.
         do
            read(u, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (.not. is_data_line(line)) cycle
            if (.not. have_layout) then
               call census_layout(line, path, cols)
               have_layout = .true.
               if (cols%header) cycle
            end if
            call parse_census_line(line, cols, path, sid, pid, dbh, height, ipft, nplant, area)
            if (.not. target_known) then ; sel_site = sid ; target_known = .true. ; end if
            if (sid /= sel_site) cycle
            if (ipass == 1) then
               !----- Pass 1: the distinct patches, first seen first, and their areas. ----------!
               ip = int(findloc(patch_id_of, pid, dim=1), ik)
               if (ip == 0_ik) then
                  patch_id_of   = [patch_id_of, pid]
                  patch_area_of = [patch_area_of, area]
               else if (cols%area > 0) then
                  if (abs(area - patch_area_of(ip)) > 1.0e-9_wp * patch_area_of(ip)) then
                     write(pid_str, '(i0)') pid
                     error stop 'meds_init: census patch_id '//trim(pid_str)//                     &
                                ' has rows with two different patch_area values in '//trim(path)
                  end if
               end if
            else
               !----- Pass 2: every cohort of the selected site, into its patch. ------------------!
               if (ipft < 1_ik .or. ipft > cfg%pft%n)                                              &
                  error stop 'meds_init: census pft index out of range in '//trim(path)
               ip = int(findloc(patch_id_of, pid, dim=1), ik)
               call add_cohort(site, cfg, ip, ipft, nplant, dbh)
            end if
         end do
      end do
      close(u)
      if (npatch == 0_ik) return

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

   !---------------------------------------------------------------------------------------!
   ! The file's layout, from its first data line: seven numbers are a header-less row read   !
   ! positionally; anything else is the header, and each name is mapped to its column.       !
   !---------------------------------------------------------------------------------------!
   subroutine census_layout(line, path, cols)
      character(len=*),       intent(in)  :: line, path
      type(census_columns_t), intent(out) :: cols
      character(len=CENSUS_FIELD_LEN) :: names(CENSUS_MAX_COLS)
      integer(ik) :: sid, pid, cid, ipft
      real(wp)    :: dbh, height, nplant
      logical     :: numeric
      integer     :: k

      call parse_census_row(line, sid, pid, cid, dbh, height, ipft, nplant, numeric)
      if (numeric) then
         cols%header = .false.
         cols%n = 7 ; cols%site = 1 ; cols%patch = 2 ; cols%cohort = 3 ; cols%dbh = 4
         cols%height = 5 ; cols%pft = 6 ; cols%nplant = 7 ; cols%area = 0
         return
      end if
      cols%header = .true.
      call split_csv(line, names, cols%n)
      if (cols%n > CENSUS_MAX_COLS) error stop 'meds_init: census header has too many columns in '//trim(path)
      do k = 1, cols%n
         select case (lower(trim(names(k))))
         case ('site_id')    ; call claim(cols%site,   k, 'site_id')
         case ('patch_id')   ; call claim(cols%patch,  k, 'patch_id')
         case ('cohort_id')  ; call claim(cols%cohort, k, 'cohort_id')
         case ('dbh')        ; call claim(cols%dbh,    k, 'dbh')
         case ('height')     ; call claim(cols%height, k, 'height')
         case ('pft')        ; call claim(cols%pft,    k, 'pft')
         case ('nplant')     ; call claim(cols%nplant, k, 'nplant')
         case ('patch_area') ; call claim(cols%area,   k, 'patch_area')
         case default
            error stop 'meds_init: census column "'//trim(names(k))//'" is unknown in '//trim(path)
         end select
      end do
      if (cols%patch  == 0) error stop 'meds_init: census has no patch_id column: '//trim(path)
      if (cols%dbh    == 0) error stop 'meds_init: census has no dbh column: '//trim(path)
      if (cols%pft    == 0) error stop 'meds_init: census has no pft column: '//trim(path)
      if (cols%nplant == 0) error stop 'meds_init: census has no nplant column: '//trim(path)
   contains
      subroutine claim(slot, k, name)
         integer,          intent(inout) :: slot
         integer,          intent(in)    :: k
         character(len=*), intent(in)    :: name
         if (slot /= 0) error stop 'meds_init: census column "'//name//'" appears twice in '//trim(path)
         slot = k
      end subroutine claim
   end subroutine census_layout

   !----- One data row under the layout. Absent optional columns: site 1, no area, height 1. -!
   subroutine parse_census_line(line, cols, path, sid, pid, dbh, height, ipft, nplant, area)
      character(len=*),       intent(in)  :: line, path
      type(census_columns_t), intent(in)  :: cols
      integer(ik),            intent(out) :: sid, pid, ipft
      real(wp),               intent(out) :: dbh, height, nplant, area
      character(len=CENSUS_FIELD_LEN) :: f(CENSUS_MAX_COLS)
      integer(ik) :: cid
      integer     :: nf, ios
      logical     :: ok

      area = 1.0_wp
      if (.not. cols%header) then
         call parse_census_row(line, sid, pid, cid, dbh, height, ipft, nplant, ok)
         if (.not. ok) error stop 'meds_init: malformed census row in '//trim(path)
      else
         call split_csv(line, f, nf)
         if (nf < cols%n) error stop 'meds_init: census row with too few columns in '//trim(path)
         sid = 1_ik ; height = 1.0_wp ; ios = 0
         if (cols%site   > 0 .and. ios == 0) read(f(cols%site),   *, iostat=ios) sid
         if (ios == 0) read(f(cols%patch),  *, iostat=ios) pid
         if (ios == 0) read(f(cols%dbh),    *, iostat=ios) dbh
         if (ios == 0) read(f(cols%pft),    *, iostat=ios) ipft
         if (ios == 0) read(f(cols%nplant), *, iostat=ios) nplant
         if (cols%height > 0 .and. ios == 0) read(f(cols%height), *, iostat=ios) height
         if (cols%area   > 0 .and. ios == 0) read(f(cols%area),   *, iostat=ios) area
         if (ios /= 0) error stop 'meds_init: malformed census row in '//trim(path)
      end if
      if (dbh <= 0.0_wp .or. height <= 0.0_wp .or. nplant <= 0.0_wp .or. area <= 0.0_wp)         &
         error stop 'meds_init: census row with non-positive dbh/height/nplant/patch_area in '//trim(path)
   end subroutine parse_census_line

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

   !----- Split a comma-separated line into its fields, each left-adjusted. nf counts them  !
   !      all, so a caller can tell a line with more fields than `fields` holds.  ----------!
   subroutine split_csv(line, fields, nf)
      character(len=*), intent(in)  :: line
      character(len=*), intent(out) :: fields(:)
      integer,          intent(out) :: nf
      integer :: k, i0, n
      n = len_trim(line) ; nf = 0 ; i0 = 1
      fields = ' '
      do k = 1, n + 1
         if (k > n) then
            nf = nf + 1
            if (nf <= size(fields)) fields(nf) = adjustl(line(i0:n))
         else if (line(k:k) == ',') then
            nf = nf + 1
            if (nf <= size(fields)) fields(nf) = adjustl(line(i0:k-1))
            i0 = k + 1
         end if
      end do
   end subroutine split_csv

   !----- Lower-case copy of an ASCII string (column names are matched case-blind). -------!
   pure function lower(s) result(t)
      character(len=*), intent(in) :: s
      character(len=len(s))        :: t
      integer :: k, c
      t = s
      do k = 1, len(s)
         c = iachar(s(k:k))
         if (c >= iachar('A') .and. c <= iachar('Z')) t(k:k) = achar(c + 32)
      end do
   end function lower

end module meds_init
