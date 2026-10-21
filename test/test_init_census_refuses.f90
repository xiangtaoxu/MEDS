! SPDX-License-Identifier: Apache-2.0
!----- A census the reader must refuse: run with `area` (a patch whose rows disagree on         !
!      patch_area) or `column` (an unknown column). The reader returns its reason, which CTest      !
!      matches; a census accepted is the failure.                                                    !
program test_init_census_refuses
   use meds_kinds,            only : ik
   use meds_config,           only : meds_config_t
   use meds_polygon_state_vars, only : polygon_t
   use meds_init_census, only : init_from_census
   use meds_test_support,     only : build_test_config
   implicit none
   type(meds_config_t) :: cfg
   type(polygon_t)     :: poly
   character(len=16)   :: which
   character(len=64)   :: path
   integer             :: u
   logical             :: found
   character(len=:), allocatable :: why

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
   call init_from_census(poly, cfg, trim(path), found, why=why)
   if (found) then
      write(*,'(a)') 'NOT REFUSED: the census was accepted'
      error stop 1
   end if
   write(*,'(4a)') 'refused: the census file ', trim(path), ' ', why
end program test_init_census_refuses
