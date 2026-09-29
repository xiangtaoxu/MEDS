! SPDX-License-Identifier: Apache-2.0
!----- A census the reader must refuse: run with `area` (a patch whose rows disagree on         !
!      patch_area) or `column` (an unknown column). CTest passes it only on the refusal's message; !
!      reaching the end means the file was accepted, which is the failure.                          !
program test_init_census_refuses
   use meds_kinds,            only : ik
   use meds_config,           only : meds_config_t
   use meds_site_state_types, only : site_t
   use meds_init,             only : init_from_census
   use meds_test_support,     only : build_test_config
   implicit none
   type(meds_config_t) :: cfg
   type(site_t)        :: site
   character(len=16)   :: which
   character(len=64)   :: path
   integer             :: u
   logical             :: found

   cfg = build_test_config()
   call get_command_argument(1, which)
   path = 'test_census_refuse_'//trim(which)//'.csv'
   open(newunit=u, file=trim(path), status='replace', action='write')
   select case (trim(which))
   case ('area')
      write(u,'(a)') 'patch_id,patch_area,dbh,pft,nplant'
      write(u,'(a)') '1,400.0,10.0,1,0.01'
      write(u,'(a)') '1,300.0,20.0,1,0.01'
   case ('column')
      write(u,'(a)') 'patch_id,basal,dbh,pft,nplant'
      write(u,'(a)') '1,0.5,10.0,1,0.01'
   end select
   close(u)
   call init_from_census(site, cfg, trim(path), found)
   write(*,'(a)') 'NOT REFUSED: the census was accepted'
   error stop 1
end program test_init_census_refuses
