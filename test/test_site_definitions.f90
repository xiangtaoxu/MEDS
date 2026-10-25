! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_site_definitions -- a run's sites from [[sites]] (docs/configuration.md, Several    !
! sites).                                                                                  !
!                                                                                          !
! The configs are derived from example 04's meds_config_eval.toml (the test runs from its  !
! directory): an override block is prepended -- meds_toml returns the FIRST match of a key !
! -- that cuts the run to three days from bare ground, without the fast loop or forcing,   !
! and writes into the work directory (the first argument).                                 !
!                                                                                          !
! 1. Two sites, 0.25 of the polygon on its own 0.6 m soil and terrain, 0.75 on [soil]: the !
!    config holds them, the polygon is built with them, its files carry the site axis and  !
!    its state file the sites, and a restart from that file reads them back.               !
! 2. A census is read by census_site_id, each site from its own rows.                      !
! 3. meds_main (the second argument) refuses shares that do not sum to 1, a table without  !
!    area_fraction, a key no table takes, and a census run without census_site_id.         !
!==========================================================================================!
program test_site_definitions
   use iso_c_binding,     only : c_int, c_size_t, c_double
   use meds_kinds,        only : wp, ik
   use meds_config,       only : meds_config_t
   use meds_config_io,    only : load_meds_config
   use meds_region,       only : meds_region_t, region_open, region_step, region_done, region_finalize, &
                                 region_free, RUN_OK
   use meds_netcdf_c
   use meds_test_support, only : banner, check, check_close, test_report
   implicit none

   character(len=*), parameter :: TEMPLATE = 'meds_config_eval.toml'
   character(len=512) :: work, meds_main
   character(len=:), allocatable :: cfg_path, nl
   type(meds_config_t) :: cfg
   type(meds_region_t) :: reg
   logical :: ok
   integer(ik) :: status, i
   real(wp) :: biggest

   call banner('sites from [[sites]]')
   call get_command_argument(1, work)
   call get_command_argument(2, meds_main)
   call execute_command_line('mkdir -p "'//trim(work)//'"')
   nl = new_line('a')
   cfg_path = trim(work)//'/sites.toml'

   !=== 1. Two sites: the config, the polygon, its files and its state. =====================!
   call derive(cfg_path, run_block('2024-07-01', '2024-07-04', 0_ik, '')//                      &
                         '[output]'//nl//'state_write = true'//nl//two_sites())
   call load_meds_config(cfg_path, cfg)
   call check(allocated(cfg%sites), 'the config holds the [[sites]]')
   if (allocated(cfg%sites)) then
      call check(size(cfg%sites) == 2, 'two [[sites]] tables are two sites')
      call check(cfg%sites(1)%own_ground .and. .not. cfg%sites(2)%own_ground,                     &
                 'a table with [soil] keys has its own ground, one without takes the run''s')
      call check_close(cfg%sites(1)%soil%depth, 0.6_wp, 1.0e-12_wp, 'the first site is 0.6 m deep')
      call check_close(cfg%sites(1)%soil%theta_sat, cfg%soil%theta_sat, 0.0_wp,                    &
                       'a key the table leaves out takes [soil]''s value')
      call check(cfg%sites(1)%has_elevation .and. .not. cfg%sites(2)%has_elevation,               &
                 'the first site has its own elevation, the second the polygon''s')
   end if

   call region_open(cfg_path, reg, ok, verbose=.false.)
   call check(ok, 'the run of two sites opens')
   associate (poly => reg%poly(1))
      call check(poly%site%n == 2_ik, 'the polygon has the two sites')
      call check_close(poly%site%area(1), 0.25_wp, 1.0e-15_wp, 'the first site is a quarter of it')
      call check_close(poly%site%area(2), 0.75_wp, 1.0e-15_wp, 'the second the rest')
      call check(poly%site%surface(1)%soil_profile%n_active < poly%site%surface(2)%soil_profile%n_active,          &
                 'the shallow site has fewer soil layers')
      call check_close(poly%site%surface(1)%elevation_m, 400.0_wp, 0.0_wp, 'the first site''s elevation')
      call check_close(poly%site%surface(2)%elevation_m, poly%elevation_m, 0.0_wp,                &
                       'the second site sits at the polygon''s elevation')
      call check(poly%site%surface(1)%slope_deg == 12.0_wp .and. poly%site%surface(1)%aspect_deg == 135.0_wp, &
                 'the first site''s slope and aspect')
      call check(all(poly%site%patch_count(1:2) > 0_ik), 'each site has bare-ground patches')
      call check(all(poly%patch%owner_site(1:poly%patch%n - 1_ik) <= poly%patch%owner_site(2:poly%patch%n)), &
                 'the patches are sorted by site')
   end associate
   do while (.not. region_done(reg))
      call region_step(reg, status)
      if (status /= RUN_OK) exit
   end do
   call check(status == RUN_OK, 'the two sites run three days')
   call region_finalize(reg, status)
   call check(status == RUN_OK, 'the run finishes with its area conserved')
   call region_free(reg)
   call check_output(trim(work)//'/sites-D-202407.nc')   ! a cohort tier's files are monthly

   !----- A restart from the state file reads the two sites back. ------------------------------!
   call derive(cfg_path, run_block('2024-07-04', '2024-07-05', 2_ik,                            &
                                   trim(work)//'/sites-S-20240704000000.nc')//two_sites())
   call region_open(cfg_path, reg, ok, verbose=.false.)
   call check(ok, 'the run restarts from its state file')
   if (ok) then
      call check(reg%poly(1)%site%n == 2_ik .and. reg%poly(1)%site%area(1) == 0.25_wp,            &
                 'the restart reads the two sites and their shares')
      call check(reg%poly(1)%site%surface(1)%soil_profile%n_active < reg%poly(1)%site%surface(2)%soil_profile%n_active, &
                 'each site keeps its soil column')
   end if
   call region_free(reg)

   !=== 2. A census read by census_site_id. ==================================================!
   call write_census(trim(work)//'/sites_census.csv')
   call derive(cfg_path, run_block('2024-07-01', '2024-07-02', 1_ik, '')//                      &
                         '[init]'//nl//'census_file = "'//trim(work)//'/sites_census.csv"'//nl//     &
                         '[[sites]]'//nl//'area_fraction = 0.4'//nl//'census_site_id = 9'//nl//       &
                         '[[sites]]'//nl//'area_fraction = 0.6'//nl//'census_site_id = 7'//nl)
   call region_open(cfg_path, reg, ok, verbose=.false.)
   call check(ok, 'a census run of two sites opens')
   if (ok) then
      associate (poly => reg%poly(1))
         call check(poly%site%n == 2_ik .and. all(poly%site%patch_count(1:2) > 0_ik),            &
                    'each site has its census patches')
         biggest = 0.0_wp
         do i = 1_ik, poly%cohort%n
            if (poly%patch%owner_site(poly%cohort%owner_patch(i)) == 1_ik) biggest = max(biggest, poly%cohort%dbh(i))
         end do
         call check(biggest > 40.0_wp, 'the first site holds census site 9''s big tree')
      end associate
   end if
   call region_free(reg)

   !=== 3. What meds_main refuses. ===========================================================!
   call refused('shares that do not sum to 1',                                                  &
                '[[sites]]'//nl//'area_fraction = 0.5'//nl//'[[sites]]'//nl//'area_fraction = 0.4'//nl, &
                'must sum to 1')
   call refused('a table without area_fraction',                                                &
                '[[sites]]'//nl//'slope = 5.0'//nl, 'sites.area_fraction')
   call refused('a key no [[sites]] table takes',                                               &
                '[[sites]]'//nl//'area_fraction = 1.0'//nl//'w_pond_max = 1.0'//nl,              &
                'sites.1.w_pond_max')
   call refused('a census run without census_site_id',                                          &
                '[init]'//nl//'init_mode = 1'//nl//'[[sites]]'//nl//'area_fraction = 1.0'//nl,    &
                'census_site_id')

   call test_report('test_site_definitions')

contains

   !----- The override block: a short run from `mode` (0 bare ground, 1 census, 2 restart), no    !
   !      fast loop or forcing, daily water and structure output into the work directory. -------!
   function run_block(start, end, mode, restart) result(s)
      character(len=*), intent(in) :: start, end, restart
      integer(ik),      intent(in) :: mode
      character(len=:), allocatable :: s
      character(len=8) :: m
      write(m, '(i0)') mode
      s = '[run]'//nl//'start_time = "'//start//'"'//nl//'end_time = "'//end//'"'//nl//            &
          'dt_slow = "1d"'//nl//'[fast]'//nl//'fast_biophysics_on = false'//nl//                   &
          'gpp_ref = 0.3'//nl//                                                                    &
          '[forcing]'//nl//'forcing_on = false'//nl//                                              &
          '[init]'//nl//'init_mode = '//trim(m)//nl//'restart_file = "'//restart//'"'//nl//         &
          '[output]'//nl//'state_file_dir = "'//trim(work)//'"'//nl//'state_file_prefix = "sites"'//nl// &
          'dir = "'//trim(work)//'"'//nl//'prefix = "sites"'//nl//                                  &
          'io_config = ""'//nl//                                                                    &
          '[output.fast]'//nl//'enabled = false'//nl//'[output.daily]'//nl//'enabled = true'//nl//  &
          'file_chunk = "run"'//nl
   end function run_block

   !----- Two sites: a quarter of the polygon on 0.6 m of soil, at 400 m on a 12-degree slope      !
   !      facing south-east; the rest on [soil]. -----------------------------------------------!
   function two_sites() result(s)
      character(len=:), allocatable :: s
      s = '[[sites]]'//nl//'area_fraction = 0.25'//nl//'depth = 0.6'//nl//'elevation = 400.0'//nl//  &
          'slope = 12.0'//nl//'aspect = 135.0'//nl//'[[sites]]'//nl//'area_fraction = 0.75'//nl
   end function two_sites

   !----- The override block, then the template. -----------------------------------------------!
   subroutine derive(path, overrides)
      character(len=*), intent(in) :: path, overrides
      character(len=1024) :: line
      integer :: uin, uout, ios
      open(newunit=uout, file=path, status='replace', action='write')
      write(uout, '(a)') overrides
      open(newunit=uin, file=TEMPLATE, status='old', action='read')
      do
         read(uin, '(a)', iostat=ios) line
         if (ios /= 0) exit
         write(uout, '(a)') trim(line)
      end do
      close(uin)
      close(uout)
   end subroutine derive

   !----- Census sites 7 and 9, two patches each; site 9 holds the one tree over 40 cm. ---------!
   subroutine write_census(path)
      character(len=*), intent(in) :: path
      integer :: u
      open(newunit=u, file=path, status='replace', action='write')
      write(u, '(a)') 'site_id,patch_id,dbh,pft,nplant'
      write(u, '(a)') '7,1,12.0,1,0.02'
      write(u, '(a)') '7,2,8.0,1,0.05'
      write(u, '(a)') '9,1,45.0,1,0.005'
      write(u, '(a)') '9,2,5.0,1,0.10'
      close(u)
   end subroutine write_census

   !----- The daily file: a site axis of two, each site's share, and the shallow site's soil      !
   !      coordinate filled below its bottom layer. --------------------------------------------!
   subroutine check_output(path)
      character(len=*), intent(in) :: path
      integer(c_int)    :: ncid, vid
      integer(c_size_t) :: n
      real(c_double)    :: area(2), z(1)
      call check(nc_open_f(path, NC_NOWRITE, ncid) == NC_NOERR, 'the daily file is written')
      call check(nc_inq_dimlen_f(ncid, 'site', n) == NC_NOERR .and. n == 2_c_size_t,              &
                 'its site axis holds the two sites')
      call nc_check(nc_inq_varid_f(ncid, 'area_si', vid), 'area_si')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t], [1_c_size_t, 2_c_size_t], area), &
                    'get site_area')
      call check(area(1) == 0.25_c_double .and. area(2) == 0.75_c_double, 'it holds each site''s share')
      call nc_check(nc_inq_varid_f(ncid, 'soil_z', vid), 'soil_z')
      call nc_check(nc_get_vara_double(ncid, vid, [0_c_size_t, 0_c_size_t, 6_c_size_t],            &
                    [1_c_size_t, 1_c_size_t, 1_c_size_t], z), 'get soil_z')
      call check(z(1) > 1.0e30_c_double, 'the shallow site''s soil coordinate stops at its bottom layer')
      call nc_check(nc_close(ncid), 'close')
   end subroutine check_output

   !----- meds_main must stop on the config, saying `message`. The case's keys come first, so they   !
   !      win over the run block's. ------------------------------------------------------------------!
   subroutine refused(what, overrides, message)
      character(len=*), intent(in) :: what, overrides, message
      character(len=:), allocatable :: path, log
      character(len=1024) :: line
      integer :: es, u, ios
      logical :: said
      path = trim(work)//'/refused.toml'
      log = trim(work)//'/refused.log'
      call derive(path, overrides//run_block('2024-07-01', '2024-07-02', 0_ik, ''))
      call execute_command_line('"'//trim(meds_main)//'" "'//path//'" > "'//log//'" 2>&1', exitstat=es)
      said = .false.
      open(newunit=u, file=log, status='old', action='read')
      do
         read(u, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (index(line, message) > 0) said = .true.
      end do
      close(u)
      call check(es /= 0 .and. said, 'meds_main refuses '//what)
   end subroutine refused

end program test_site_definitions
