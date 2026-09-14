!==========================================================================================!
! test_io_config_example -- the shipped meds_io_config.toml must list EVERY registered variable   !
! (#241).                                                                                        !
!                                                                                          !
! That file is a GENERATED artifact (`meds_main --dump-io-config`) which is nonetheless TRACKED,   !
! because it is the shipped example of the per-variable override surface that                      !
! docs/science/diagnostics.md points users at. It is a pure function of the registry, so it can    !
! only ever be right or stale -- and it went stale across five pull requests, ending up short by   !
! twenty variables, because nothing compared the two.                                              !
!                                                                                          !
! It degrades QUIETLY, which is why: an unknown name in the file is a hard error, but a MISSING    !
! name is not an error at all -- that variable simply keeps its registry-default streams. A user   !
! who copies the example as a starting point and then switches everything off would silently fail  !
! to switch off the variables the file forgot.                                                     !
!                                                                                          !
! The check is name presence, not a byte diff: regenerating rewrites every stream mask and comment !
! line, so a byte comparison would fail on cosmetic churn and get suppressed. What matters is that !
! the surface is complete.                                                                          !
!==========================================================================================!
program test_io_config_example
   use meds_kinds,           only : ik
   use meds_config,          only : meds_config_t
   use meds_output_types,    only : output_registry_t
   use meds_output_registry, only : build_output_registry
   use meds_test_support,    only : banner, build_test_config, check, test_report
   implicit none

   type(output_registry_t)       :: reg
   type(meds_config_t)           :: cfg
   character(len=:), allocatable :: path, body
   character(len=512)            :: arg
   integer(ik)                   :: k, n_missing
   logical                       :: ok

   call banner('the shipped meds_io_config.toml example is complete (#241)')

   !----- CMake passes the SOURCE-tree path; ctest runs from the build tree. -------------!
   call get_command_argument(1, arg)
   path = trim(arg)
   call check(len_trim(path) > 0, 'the example file path was passed on the command line')

   call slurp(path, body, ok)
   call check(ok, 'the shipped meds_io_config.toml exists and is readable')
   if (.not. ok) then
      call test_report('test_io_config_example')
      stop
   end if

   cfg = build_test_config()
   call build_output_registry(reg, cfg)
   call check(reg%nvar > 0_ik, 'the registry has variables to compare against')

   n_missing = 0_ik
   do k = 1_ik, reg%nvar
      !----- The dump writes one commented assignment per variable: "# <name> = ...". Matching on  !
      !      "# name = " rather than the bare name avoids a short name matching inside a longer     !
      !      one (gpp_site inside gpp_site_something) or inside a long_name or units string.        !
      if (index(body, '# '//trim(reg%var(k)%name)//' = ') == 0) then
         n_missing = n_missing + 1_ik
         if (n_missing <= 10_ik) write(*,'(a)') '   MISSING: '//trim(reg%var(k)%name)
      end if
   end do
   if (n_missing > 10_ik) write(*,'(a,i0,a)') '   ... and ', n_missing - 10_ik, ' more'

   call check(n_missing == 0_ik,                                                                  &
              'every registered variable appears in the shipped example (regenerate with '//       &
              'meds_main --dump-io-config)')
   write(*,'(a,i0,a)') '   (', reg%nvar, ' registered variables, all present)'

   call test_report('test_io_config_example')

contains

   subroutine slurp(p, s, good)
      character(len=*),              intent(in)  :: p
      character(len=:), allocatable, intent(out) :: s
      logical,                       intent(out) :: good
      integer(ik)        :: u, ios
      character(len=1024) :: line
      good = .false. ; s = ''
      open(newunit=u, file=p, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      do
         read(u,'(a)', iostat=ios) line
         if (ios /= 0) exit
         s = s//trim(line)//new_line('a')
      end do
      close(u)
      good = .true.
   end subroutine slurp

end program test_io_config_example
