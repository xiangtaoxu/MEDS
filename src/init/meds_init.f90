! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_init -- helpers to build the initial community a run starts from.                    !
!                                                                                          !
! Two entry points exist: a near-bare-ground site (for spin-up demos and the test suite),  !
! which seed_bare_ground plants with each PFT's nbg_init_density, and a start from a cohort !
! CENSUS file (init_from_census) -- a CSV produced by a previous MEDS run or a field        !
! inventory. add_cohort/finalize_init are the low-level builders both share.               !
!==========================================================================================!
module meds_init
   use meds_kinds,      only : wp, ik
   use meds_config,     only : meds_config_t, growth_window_steps
   use meds_allometry_lib, only : height_to_dbh
   use meds_cohort_state_vars,     only : expand_cohort_block, init_cohort
   use meds_site_state_vars,       only : make_site_surface
   use meds_polygon_state_vars,    only : polygon_t, alloc_polygon, rebuild_patch_map, rebuild_cohort_map,      &
                                          assign_cohort_id, assign_patch_id
   use meds_patch_state_vars,      only : patch_physical_state, DIST_PRIMARY
   use meds_demography_cohort_fusefiss, only : sort_cohorts, new_fuse_cohorts, terminate_cohorts,     &
                                               split_cohorts
   use meds_demography_patch_fusefiss,  only : sort_patches, new_fuse_patches, terminate_patches
   implicit none
   private

   public :: init_bare_ground, init_bare_sites, seed_bare_ground, add_cohort, finalize_init, init_from_census
   public :: restructure_census_stand

   !----- The census file's column layout: each index is the column's position, 0 if absent. -!
   integer, parameter :: CENSUS_MAX_COLS = 32, CENSUS_FIELD_LEN = 64
   !----- Two rows of one patch agree on its area within this fraction of it. ------------------!
   real(wp), parameter :: PATCH_AREA_RTOL = 1.0e-9_wp
   type :: census_columns_t
      logical :: header = .false.
      integer :: n = 0
      integer :: site = 0, patch = 0, cohort = 0, dbh = 0, height = 0, pft = 0, nplant = 0, area = 0
   end type census_columns_t

contains

   !---------------------------------------------------------------------------------------!
   ! A bare polygon of one site: n_patch identical empty patches sharing the site's area.  !
   !---------------------------------------------------------------------------------------!
   subroutine init_bare_ground(poly, cfg, n_patch)
      type(polygon_t),     intent(out) :: poly
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: n_patch
      call init_bare_sites(poly, cfg, [n_patch], [cfg], [1.0_wp])
   end subroutine init_bare_ground

   !-------------------------------------------------------------------------------------------!
   ! A near-bare-ground polygon of several sites: site `is` takes its share `site_area(is)` of !
   ! the polygon, its surface from `site_cfg(is)` (its soil, from [soil] or its own maps) and  !
   ! n_patch(is) identical empty patches sharing its area. The patches are numbered site by    !
   ! site, so the block is sorted by site. Everything else comes from `cfg`, the run's.        !
   !-------------------------------------------------------------------------------------------!
   subroutine init_bare_sites(poly, cfg, n_patch, site_cfg, site_area)
      type(polygon_t),     intent(out) :: poly
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: n_patch(:)     !< patches per site
      type(meds_config_t), intent(in)  :: site_cfg(:)
      real(wp),            intent(in)  :: site_area(:)   !< [-] fractions of the polygon
      integer(ik) :: ip, is, k, ns
      ns = int(size(n_patch), ik)
      call alloc_polygon(poly, cfg%pft%n, coh_cap = 64_ik, pat_cap = max(sum(n_patch), 1_ik),  &
                      n_growth_window = growth_window_steps(cfg), pft = cfg%pft,                    &
                      n_site = ns)
      poly%patch%n = sum(n_patch)
      ip = 0_ik
      do is = 1_ik, ns
         call make_site_surface(poly%site%surface(is), site_cfg(is))
         poly%site%area(is) = site_area(is)
         do k = 1_ik, n_patch(is)
            ip = ip + 1_ik
            poly%patch%area(ip)       = 1.0_wp / real(n_patch(is), wp)
            poly%patch%age(ip)        = 0.0_wp
            poly%patch%dist_type(ip)  = DIST_PRIMARY
            poly%patch%owner_site(ip) = is
            call assign_patch_id(poly, ip)
         end do
      end do
      poly%patch%recruit_pool = 0.0_wp
      poly%cohort%n = 0_ik
      call rebuild_patch_map(poly)
      call rebuild_cohort_map(poly)
   end subroutine init_bare_sites

   !---------------------------------------------------------------------------------------!
   ! A near-bare-ground start's seedlings (ED2's init_density): every patch gets one cohort of !
   ! each PFT with nbg_init_density > 0, that many plants per m2, at min_cohort_height on the  !
   ! PFT's own height curve, where recruits are born. All densities 0 leave the site empty.    !
   !---------------------------------------------------------------------------------------!
   subroutine seed_bare_ground(poly, cfg)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      real(wp)    :: dbh
      integer(ik) :: ip, ipft
      if (.not. any(cfg%pft%nbg_init_density > 0.0_wp)) return
      do ip = 1_ik, poly%patch%n
         do ipft = 1_ik, cfg%pft%n
            if (cfg%pft%nbg_init_density(ipft) <= 0.0_wp) cycle
            dbh = height_to_dbh(poly%cohort%pft_table%allom(ipft), cfg%pft%min_cohort_height)
            call add_cohort(poly, cfg, ip, ipft, cfg%pft%nbg_init_density(ipft), dbh)
         end do
      end do
      call finalize_init(poly)
   end subroutine seed_bare_ground

   !---------------------------------------------------------------------------------------!
   ! Append one cohort to patch ip (caller invokes finalize_init afterwards).               !
   !---------------------------------------------------------------------------------------!
   subroutine add_cohort(poly, cfg, ip, ipft, nplant, dbh)
      type(polygon_t),  intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      integer(ik),         intent(in)    :: ip, ipft
      real(wp),            intent(in)    :: nplant, dbh
      integer(ik) :: m
      real(wp)    :: t_birth, psi_birth
      call expand_cohort_block(poly%cohort, poly%cohort%n + 1_ik)
      m = poly%cohort%n + 1_ik
      !----- Born in its patch's state (patch_physical_state); at setup the patch is not seeded yet, so !
      !      that is the run's starting state, and init_fast_reservoirs re-seeds it from the patch. -!
      call patch_physical_state(poly%patch%cas(ip), poly%patch%soil_water(ip), cfg,                   &
                                poly%site%surface(poly%patch%owner_site(ip)), t_birth, psi_birth)
      call init_cohort(poly%cohort, m, cfg%pft, ipft, ip, nplant, dbh, t_birth, psi_birth)
      call assign_cohort_id(poly, m)
      poly%cohort%n = m
   end subroutine add_cohort

   subroutine finalize_init(poly)
      type(polygon_t), intent(inout) :: poly
      call rebuild_cohort_map(poly)
      call sort_cohorts(poly)
   end subroutine finalize_init

   !---------------------------------------------------------------------------------------!
   ! Initialize a polygon from a cohort CENSUS file: a CSV with one row per cohort, produced    !
   ! by a previous MEDS run or from a field inventory.                                          !
   !                                                                                          !
   ! COLUMNS ARE MATCHED BY NAME from the header line, in any order:                           !
   !   required  patch_id, dbh [cm], pft, nplant [plants per m2 of the patch]                  !
   !   optional  patch_area [any unit], site_id, cohort_id, height [m]                         !
   ! Each distinct (site_id, patch_id) is one primary patch of age 0. With patch_area a site's  !
   ! patches take their areas normalized to fractions of the site; without it they share it     !
   ! equally. One site is built from the rows of one site_id (`use_site_id` if present, else    !
   ! the first one seen); with `site_ids`, site k is built from the rows of site_ids(k), with   !
   ! its surface from site_cfg(k) and its share site_area(k) of the polygon. A file without     !
   ! the column is one site. `dbh` drives the allometry (set_cohort_size), so `height` and      !
   ! `cohort_id` are provenance: height, when present, must be positive.                        !
   !                                                                                          !
   ! A file whose first data line is seven numbers is read positionally as                    !
   ! `site_id, patch_id, cohort_id, dbh, height, pft, nplant`, the layout of a census file    !
   ! without a header.                                                                        !
   !                                                                                          !
   ! Blank lines and lines beginning with '#' are ignored. found=.false. (the site is left      !
   ! unbuilt) if the file cannot be opened or holds no row for the selected site, and `why`     !
   ! says which, so the caller can stop and say so. An unknown or duplicated column, a missing  !
   ! required one, a malformed row, a non-positive value, or a patch whose rows disagree on     !
   ! patch_area stops the run. The file is read twice (patches, then cohorts), so no row buffer !
   ! is needed.                                                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine init_from_census(poly, cfg, path, found, use_site_id, why, site_ids, site_cfg, site_area)
      type(polygon_t),     intent(out) :: poly
      type(meds_config_t), intent(in)  :: cfg
      character(len=*),    intent(in)  :: path
      logical,             intent(out) :: found
      integer(ik), intent(in), optional :: use_site_id
      character(len=:), allocatable, intent(out), optional :: why   !< why found = .false.; '' when found
      integer(ik),         intent(in), optional :: site_ids(:)      !< several sites: each one's census site_id,
      type(meds_config_t), intent(in), optional :: site_cfg(:)      !< configuration
      real(wp),            intent(in), optional :: site_area(:)     !< and share of the polygon

      integer,     parameter :: LINELEN = 1024
      character(len=LINELEN) :: line
      type(census_columns_t) :: cols
      integer(ik)            :: sid, pid, ipft, npatch, ip, k, i, j
      real(wp)               :: dbh, height, nplant, area
      integer(ik), allocatable :: sel(:), site_of(:), patch_id_of(:), slot_of(:), n_patch(:)
      real(wp),    allocatable :: patch_area_of(:), site_patch_area(:)
      logical                :: target_known, have_layout
      integer                :: u, ios, ipass
      character(len=24)      :: pid_str

      found = .false.
      if (present(why)) why = 'cannot be opened'
      open(newunit=u, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      if (present(why)) why = 'holds no cohort for the site'

      !----- The census site of each site built: the given ones, or the first one the file names. -!
      target_known = .true.
      if (present(site_ids)) then
         sel = site_ids
      else if (present(use_site_id)) then
         sel = [use_site_id]
      else
         allocate(sel(0))
         target_known = .false.
      end if
      allocate(site_of(0), patch_id_of(0), patch_area_of(0))
      npatch = 0_ik

      do ipass = 1, 2
         if (ipass == 2) then
            npatch = int(size(patch_id_of), ik)
            if (npatch == 0_ik) exit                         ! no rows for the selected sites
            !----- The patches, site by site, first seen first within a site: equal areas, or   !
            !      the file's areas normalized to their site. ----------------------------------!
            allocate(slot_of(npatch), n_patch(size(sel)))
            j = 0_ik
            do k = 1_ik, int(size(sel), ik)
               n_patch(k) = int(count(site_of == k), ik)
               do i = 1_ik, npatch
                  if (site_of(i) /= k) cycle
                  j = j + 1_ik
                  slot_of(i) = j
               end do
            end do
            if (present(site_ids)) then
               call init_bare_sites(poly, cfg, n_patch, site_cfg, site_area)
            else
               call init_bare_ground(poly, cfg, npatch)
            end if
            if (cols%area > 0) then
               do k = 1_ik, int(size(sel), ik)
                  site_patch_area = pack(patch_area_of, site_of == k)
                  site_patch_area = site_patch_area / sum(site_patch_area)
                  j = 0_ik
                  do i = 1_ik, npatch
                     if (site_of(i) /= k) cycle
                     j = j + 1_ik
                     poly%patch%area(slot_of(i)) = site_patch_area(j)
                  end do
               end do
            end if
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
            if (.not. target_known) then
               sel = [sid]
               target_known = .true.
            end if
            k = int(findloc(sel, sid, dim=1), ik)
            if (k == 0_ik) cycle
            ip = census_patch(k, pid)
            if (ipass == 1) then
               !----- Pass 1: the distinct patches, first seen first, and their areas. ----------!
               if (ip == 0_ik) then
                  site_of       = [site_of, k]
                  patch_id_of   = [patch_id_of, pid]
                  patch_area_of = [patch_area_of, area]
               else if (cols%area > 0) then
                  if (abs(area - patch_area_of(ip)) > PATCH_AREA_RTOL * patch_area_of(ip)) then
                     write(pid_str, '(i0)') pid
                     error stop 'meds_init: census patch_id '//trim(pid_str)//                     &
                                ' has rows with two different patch_area values in '//trim(path)
                  end if
               end if
            else
               !----- Pass 2: every cohort of the selected sites, into its patch. -----------------!
               if (ipft < 1_ik .or. ipft > cfg%pft%n)                                              &
                  error stop 'meds_init: census pft index out of range in '//trim(path)
               call add_cohort(poly, cfg, slot_of(ip), ipft, nplant, dbh)
            end if
         end do
      end do
      close(u)
      if (npatch == 0_ik) return

      call finalize_init(poly)
      found = .true.
      if (present(why)) why = ''

   contains

      !----- The patch of site k with census patch_id pid, in first-seen order (0 if new). ----!
      integer(ik) function census_patch(k, pid) result(ip)
         integer(ik), intent(in) :: k, pid
         integer(ik) :: i
         ip = 0_ik
         do i = 1_ik, int(size(patch_id_of), ik)
            if (site_of(i) == k .and. patch_id_of(i) == pid) then
               ip = i
               return
            end if
         end do
      end function census_patch
   end subroutine init_from_census

   !---------------------------------------------------------------------------------------!
   ! Restructure a census stand before the first step, with the slow step's own operators:    !
   ! its monthly cohort block, then its yearly patch block, less recruitment and disturbance.  !
   ! A census arrives with a cohort per measured size and a patch per plot cell, far more of   !
   ! both than a step should carry, and nothing else fuses them before the first month or year !
   ! boundary. Each block runs under the same switches as in the slow step.                    !
   !---------------------------------------------------------------------------------------!
   subroutine restructure_census_stand(poly, cfg)
      type(polygon_t),     intent(inout) :: poly
      type(meds_config_t), intent(in)    :: cfg
      if (cfg%demography_on .and. cfg%do_cohort_fissfuse) then
         call new_fuse_cohorts(poly, cfg)
         call terminate_cohorts(poly, cfg)
         call split_cohorts(poly, cfg)
         call sort_cohorts(poly)
      end if
      if (cfg%demography_on .and. cfg%do_patch_fissfuse) then
         call sort_patches(poly)
         call new_fuse_patches(poly, cfg)
         call terminate_patches(poly, cfg)
         call new_fuse_cohorts(poly, cfg)
         call terminate_cohorts(poly, cfg)
         call sort_cohorts(poly)
      end if
   end subroutine restructure_census_stand

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
