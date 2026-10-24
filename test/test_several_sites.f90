! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_several_sites -- a polygon of several sites.                                        !
!                                                                                          !
! 1. Two identical sites of areas 0.3 and 0.7 each evolve as a polygon of one site does,   !
!    over 60 days of fast and slow steps on a short calendar -- a month boundary every 15  !
!    days, a year boundary every 30, each with a disturbance and patch fusion: every       !
!    site's stand is the one-site stand, and the polygon totals are the one site's to      !
!    round-off.                                                                            !
! 2. Two sites on different soils each evolve as a one-site polygon on its own soil, and   !
!    their soil columns differ.                                                            !
! 3. The patch loop gives the same polygon on one thread and on four.                      !
! 4. A state file keeps the sites: their areas and soils, and each patch's site.           !
!==========================================================================================!
program test_several_sites
   use meds_kinds,                only : wp, ik
   use meds_constants,            only : rho_h2o
   use meds_config,               only : meds_config_t
   use meds_polygon_state_vars,   only : polygon_t
   use meds_init_bare_ground, only : init_bare_ground, init_bare_sites, add_cohort, finalize_init
   use meds_init,                 only : init_fast_state
   use meds_stepper,              only : advance_one_step, run_ledgers_t
   use meds_output_reduce,    only : total_agb, total_nplant, total_area
   use meds_init_restart,                   only : restart_file_t, restart_file_open, restart_file_write, restart_file_close, &
                                         read_restart_file
   use meds_time,                 only : meds_time_t
   use meds_test_support,         only : banner, build_test_config, check, check_close, test_report
   implicit none

   integer(ik), parameter :: N_DAY = 60_ik            ! two short years of two 15-day months
   real(wp),    parameter :: RTOL  = 1.0e-12_wp   !< relative: check_close scales it by the expected value
   type(meds_config_t)  :: cfg, cfg_shallow
   type(polygon_t)      :: one, two, two_mt, one_shallow, mixed, back
   type(meds_time_t)    :: now, t_read
   logical              :: found
   integer(ik)          :: is

   call banner('several sites per polygon')
   cfg = build_test_config()
   cfg%fast_biophysics_on = .true.
   cfg_shallow = cfg
   cfg_shallow%soil%depth = 0.6_wp

   !=== 1. Two identical sites, 0.3 and 0.7 of the polygon. ================================!
   call stand(one, cfg, [cfg], [1.0_wp])
   call stand(two, cfg, [cfg, cfg], [0.3_wp, 0.7_wp])
   call stand(two_mt, cfg, [cfg, cfg], [0.3_wp, 0.7_wp])
   call run(one, 1_ik)
   call run(two, 1_ik)
   call run(two_mt, 4_ik)
   call check(one%patch%n > 2_ik, 'the year boundary added a disturbance gap to the one site')
   do is = 1_ik, 2_ik
      call check(two%site%patch_count(is) == one%patch%n,                                            &
                 'each identical site has the one site''s patches after the disturbance and fusion')
      call check_close(site_carbon(two, is), site_carbon(one, 1_ik), RTOL,                             &
                       'each identical site holds the one site''s plant carbon')
      call check_close(site_water(two, is), site_water(one, 1_ik), RTOL,                              &
                       'each identical site holds the one site''s soil water')
   end do
   call check(all(two%patch%owner_site(1:two%patch%n) == [(1_ik, is = 1_ik, two%site%patch_count(1)),  &
                                                            (2_ik, is = 1_ik, two%site%patch_count(2))]), &
              'the patch block stays sorted by site through the disturbance and fusion')
   call check_close(total_agb(two), total_agb(one), RTOL,                                           &
                    'the polygon''s AGB is the one site''s to round-off')
   call check_close(total_nplant(two), total_nplant(one), RTOL,                                     &
                    'the polygon''s stem density is the one site''s to round-off')
   call check_close(total_area(two), 1.0_wp, 1.0e-12_wp, 'the patches still cover the polygon')
   do is = 1_ik, 2_ik
      associate (o => two%site%patch_offset(is), c => two%site%patch_count(is))
         call check_close(sum(two%patch%area(o:o + c - 1_ik)), 1.0_wp, 1.0e-12_wp,                  &
                          'each site''s patch areas still sum to 1')
      end associate
   end do

   !=== 3. One thread or four: the same polygon. ===========================================!
   call check(all(two_mt%patch%area(1:two%patch%n) == two%patch%area(1:two%patch%n)) .and.        &
              site_carbon(two_mt, 1_ik) == site_carbon(two, 1_ik) .and.                             &
              site_carbon(two_mt, 2_ik) == site_carbon(two, 2_ik) .and.                             &
              site_water(two_mt, 2_ik) == site_water(two, 2_ik),                                    &
              'four patch threads give the polygon one thread gives, bit for bit')

   !=== 2. Two soils: each site evolves on its own. =========================================!
   call stand(one_shallow, cfg_shallow, [cfg_shallow], [1.0_wp])
   call stand(mixed, cfg, [cfg, cfg_shallow], [0.5_wp, 0.5_wp])
   call run(one_shallow, 1_ik)
   call run(mixed, 1_ik)
   call check(mixed%site%surface(1)%soil_profile%n_active /= mixed%site%surface(2)%soil_profile%n_active,           &
              'the two sites have different soil columns')
   call check(abs(site_water(mixed, 1_ik) - site_water(mixed, 2_ik)) > 1.0_wp,                      &
              'the two soils hold different water')
   call check_close(site_carbon(mixed, 1_ik), site_carbon(one, 1_ik), RTOL,                         &
                    'the deep site grows the one-site stand on the deep soil')
   call check_close(site_carbon(mixed, 2_ik), site_carbon(one_shallow, 1_ik), RTOL,                 &
                    'the shallow site grows the one-site stand on its soil')
   call check_close(site_water(mixed, 2_ik), site_water(one_shallow, 1_ik), RTOL,                   &
                    'the shallow site keeps the one-site water on its soil')

   !=== 4. The mixed polygon through a state file. ==========================================!
   now = meds_time_t(year=2000, month=3, day=1)
   call write_state(mixed, now)
   call read_restart_file(back, cfg, './test_several_sites-S-20000301000000.nc', t_read, found,       &
                      polygon=1_ik, site_cfg=[cfg, cfg_shallow])
   call check(found .and. back%site%n == 2_ik, 'the state file holds the two sites')
   call check(all(back%site%area(1:2) == mixed%site%area(1:2)), 'each site keeps its area')
   call check(back%site%surface(2)%soil_profile%n_active == mixed%site%surface(2)%soil_profile%n_active,        &
              'the shallow site keeps its soil column')
   call check(all(back%patch%owner_site(1:back%patch%n) == mixed%patch%owner_site(1:mixed%patch%n)), &
              'each patch keeps its site')
   do is = 1_ik, 2_ik
      call check(back%site%patch_count(is) == mixed%site%patch_count(is), 'each site keeps its patches')
      call check_close(site_carbon(back, is), site_carbon(mixed, is), RTOL,                         &
                       'each site keeps its plant carbon')
      call check(site_water(back, is) == site_water(mixed, is), 'each site keeps its soil water')
   end do

   call test_report('test_several_sites')

contains

   !----- A two-patch stand in each site: a canopy tree and an understorey in one patch, a     !
   !      sapling in the other. ---------------------------------------------------------------!
   subroutine stand(poly, run_cfg, site_cfg, site_area)
      type(polygon_t),     intent(out) :: poly
      type(meds_config_t), intent(in)  :: run_cfg, site_cfg(:)
      real(wp),            intent(in)  :: site_area(:)
      integer(ik) :: k, ns
      ns = int(size(site_cfg), ik)
      call init_bare_sites(poly, run_cfg, [(2_ik, k = 1_ik, ns)], site_cfg, site_area)
      do k = 0_ik, ns - 1_ik
         call add_cohort(poly, run_cfg, 2_ik*k + 1_ik, 1_ik, 0.05_wp, 30.0_wp)
         call add_cohort(poly, run_cfg, 2_ik*k + 1_ik, 2_ik, 0.20_wp,  4.0_wp)
         call add_cohort(poly, run_cfg, 2_ik*k + 2_ik, 3_ik, 0.50_wp,  2.0_wp)
      end do
      call finalize_init(poly)
      call init_fast_state(poly, run_cfg)
   end subroutine stand

   !----- N_DAY daily steps, a month's restructuring owed every 15 days and a year's every 30,   !
   !      which the next step performs first, as a run does. -------------------------------------!
   subroutine run(poly, patch_threads)
      type(polygon_t), intent(inout) :: poly
      integer(ik),     intent(in)    :: patch_threads
      type(run_ledgers_t) :: ledgers
      integer(ik) :: day
      do day = 1_ik, N_DAY
         call advance_one_step(poly, cfg, ledgers, patch_threads=patch_threads)
         poly%restructure_pending  = mod(day, 15_ik) == 0_ik
         poly%restructure_new_year = mod(day, 30_ik) == 0_ik
      end do
   end subroutine run

   !----- The polygon as the one polygon of a state file (id 0, at 0 N 0 E). ------------------!
   subroutine write_state(poly, t)
      type(polygon_t),   intent(in) :: poly
      type(meds_time_t), intent(in) :: t
      type(restart_file_t) :: sf
      call restart_file_open(sf, '.', 'test_several_sites', t, poly%n_pft,                          &
                           int(size(poly%cohort%growth_hist, 1), ik), poly%cohort%n, poly%patch%n,    &
                           poly%site%n, polygon_id=[0_ik], polygon_lat=[0.0_wp], polygon_lon=[0.0_wp])
      call restart_file_write(sf, poly, t, 1_ik)
      call restart_file_close(sf)
   end subroutine write_state

   !----- Plant carbon [kgC/m2 of the site]: its patches' cohorts, weighted by patch area. ----!
   real(wp) function site_carbon(poly, is) result(c)
      type(polygon_t), intent(in) :: poly
      integer(ik),     intent(in) :: is
      integer(ik) :: ip, i
      c = 0.0_wp
      do ip = poly%site%patch_offset(is), poly%site%patch_offset(is) + poly%site%patch_count(is) - 1_ik
         do i = poly%patch%cohort_offset(ip), poly%patch%cohort_offset(ip) + poly%patch%cohort_count(ip) - 1_ik
            c = c + poly%patch%area(ip) * poly%cohort%nplant(i) * (poly%cohort%leaf_carbon(i)           &
                  + poly%cohort%fineroot_carbon(i) + poly%cohort%wood_carbon(i)                       &
                  + poly%cohort%nonstructural_carbon(i))
         end do
      end do
   end function site_carbon

   !----- Soil water [kg/m2 of the site] over its patches' active layers. ---------------------!
   real(wp) function site_water(poly, is) result(w)
      type(polygon_t), intent(in) :: poly
      integer(ik),     intent(in) :: is
      integer(ik) :: ip, k
      w = 0.0_wp
      associate (soil => poly%site%surface(is)%soil_profile)
         do ip = poly%site%patch_offset(is), poly%site%patch_offset(is) + poly%site%patch_count(is) - 1_ik
            do k = 1_ik, soil%n_active
               w = w + poly%patch%area(ip) * poly%patch%soil_water(ip)%theta(k) * soil%dz(k) * rho_h2o
            end do
         end do
      end associate
   end function site_water

end program test_several_sites
