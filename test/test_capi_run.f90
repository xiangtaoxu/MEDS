!==========================================================================================!
! test_capi_run -- COVERAGE FOR THE FULL-COUPLED-MODEL C-API SHIM (`meds.model`).            !
!                                                                                          !
! Two things here are load-bearing and neither is an assertion:                                !
!                                                                                          !
!   * `src/capi/meds_capi_run.f90` is COMPILED by a MANDATORY target, so a change to           !
!     `meds_driver`'s open/step/finalize seam is a build break rather than a silent rot in an   !
!     optional library. That is the §7.6 #4 rule, and it exists because `meds_capi_demography`  !
!     spent a day not compiling while 42/42 stayed green.                                        !
!                                                                                          !
!   * this test runs from `examples/example_biophysics/` and LOADS BOTH of that example's        !
!     shipped configs. Before this, nothing in the build touched them -- and the sibling          !
!     example's config had drifted 37 required keys behind the schema exactly that way. A config   !
!     nothing loads rots like a shim nothing compiles.                                              !
!                                                                                          !
! The assertions cover the shim's own logic: the invalid-handle guard (every accessor has to      !
! return a sentinel rather than index a stale registry slot), distinct handles, and a real         !
! open -> step -> finalize -> free cycle through the C boundary.                                    !
!                                                                                          !
! The cycle runs on a config DERIVED from the shipped spin-up one: fast biophysics and live met     !
! forcing off (the forcing netCDF is not tracked in git), bare ground, three days. The derivation   !
! prepends an override block and then copies the original verbatim -- `meds_toml`'s lookup returns  !
! the FIRST match for a key, so the prepended values win. Deriving it, rather than writing a config !
! from scratch, means this test also proves the shipped file is complete enough to INITIALISE, not  !
! merely to parse.                                                                                   !
!==========================================================================================!
program test_capi_run
   use, intrinsic :: iso_c_binding, only : c_double, c_int, c_char
   use meds_kinds,        only : wp
   use meds_config,       only : meds_config_t
   use meds_config_io,    only : load_meds_config
   use meds_capi_run,     only : meds_run_open, meds_run_step, meds_run_is_done,                 &
                                 meds_run_finalize, meds_run_free, meds_run_year, meds_run_day,  &
                                 meds_run_istep, meds_run_n_patch, meds_run_n_cohort,            &
                                 meds_run_total_agb, meds_run_total_lai, meds_run_soil_carbon
   use meds_test_support, only : check, banner
   implicit none

   character(len=*), parameter :: CFG_SPINUP = 'meds_config_spinup.toml'
   character(len=*), parameter :: CFG_JULY   = 'meds_config_july.toml'
   character(len=*), parameter :: CFG_DERIVED = 'test_capi_run_derived.toml'
   integer(c_int)      :: h, h2, status, d0, d1, nsteps
   real(c_double)      :: agb, lai, soilc
   type(meds_config_t) :: cfg

   call banner('full-model C-API shim (meds.model)')

   !=== 1. Both SHIPPED example_biophysics configs load. ====================================!
   !      The schema-drift guard. This example is driven by a shell script and three plotting !
   !      scripts, none of which the build knows about, so until now nothing here would have  !
   !      noticed the configs falling behind the config reader.                               !
   call load_meds_config(CFG_SPINUP, cfg)
   call check(cfg%pft%n >= 1_c_int, 'the shipped spin-up config loads and carries PFTs')
   call check(cfg%fast_biophysics_on, 'the spin-up config still requests fast biophysics')
   call check(cfg%soil_carbon_on, 'the spin-up config still requests SOIL CARBON')
   call load_meds_config(CFG_JULY, cfg)
   call check(cfg%pft%n >= 1_c_int, 'the shipped July config loads and carries PFTs')
   call check(cfg%soil_carbon_on, 'the July config still requests SOIL CARBON')
   !----- Both stages MUST agree on soil carbon. They did not, which is how the example came to !
   !      report Rh = 0 and discard its litter: with the feature off there is no soil pool to    !
   !      receive it, and the slow ledger flagged the discarded carbon as an undeclared residual.!
   !      A July stage restarting from a spin-up that never built the pools respires nothing     !
   !      whatever this flag says, so the two have to move together.  ---------------------------!

   !=== 2. The invalid-handle guard. ========================================================!
   !      Every accessor takes a bare integer from the caller. One that was never opened, or  !
   !      already freed, must return a sentinel -- not index a stale registry slot.            !
   call check(meds_run_step(99_c_int)      == -1_c_int, 'step on a bogus handle returns -1')
   call check(meds_run_is_done(99_c_int)   ==  1_c_int, 'is_done on a bogus handle reports done')
   call check(meds_run_finalize(0_c_int)   == -1_c_int, 'finalize on handle 0 returns -1')
   call check(meds_run_year(-3_c_int)      == -1_c_int, 'a negative handle yields the year sentinel')
   call check(meds_run_total_agb(99_c_int) == 0.0_wp,   'aggregates on a bogus handle read 0')

   !=== 3. A real open -> step -> finalize -> free cycle across the C boundary. ==============!
   call write_derived_config()
   h = meds_run_open(to_c(CFG_DERIVED), int(len(CFG_DERIVED), c_int), 0_c_int)
   call check(h > 0_c_int, 'a run opens from a config through the C-API')
   call check(meds_run_is_done(h) == 0_c_int, 'a freshly opened run is not already done')
   call check(meds_run_n_patch(h) >= 1_c_int, 'the opened run reports its patch count')

   h2 = meds_run_open(to_c(CFG_DERIVED), int(len(CFG_DERIVED), c_int), 0_c_int)
   call check(h2 > 0_c_int .and. h2 /= h, 'two live runs get DISTINCT handles')
   call meds_run_free(h2)
   call check(meds_run_is_done(h2) == 1_c_int, 'a freed handle reads as done, not as a live run')

   !----- Step to the end. The CALLER owns the loop -- that is the whole point of this shim. ---!
   d0 = meds_run_day(h)
   nsteps = 0_c_int
   do while (meds_run_is_done(h) == 0_c_int)
      status = meds_run_step(h)
      call check(status == 0_c_int, 'each step reports DRIVER_OK')
      nsteps = nsteps + 1_c_int
      if (nsteps > 10_c_int) exit                     ! the derived config is a 3-day window
   end do
   d1 = meds_run_day(h)
   call check(nsteps == 3_c_int, 'a 3-day window takes exactly 3 daily steps')
   call check(d1 /= d0, 'the calendar advanced')
   call check(meds_run_istep(h) == nsteps, 'the driver istep counter matches the steps taken')

   !----- Stepping past the end is a STATUS, not a crash: the Python loop relies on it. -------!
   call check(meds_run_step(h) == 1_c_int, 'a step past end_time returns DRIVER_FINISHED')

   !=== 3b. A REUSED registry slot must not inherit the previous run's counters. =============!
   !      The registry is module-`save`, so a freed handle is handed out again -- and the run       !
   !      object in that slot still holds whatever the last run left in it. `meds_run_t`'s default  !
   !      initialisers do NOT help: they apply to a fresh variable, not to a slot being reopened.   !
   !                                                                                          !
   !      This is not hypothetical. examples/example_biophysics opens the spin-up, closes it, and   !
   !      opens the July stage in the SAME process; before driver_open reset these, stage 2         !
   !      inherited stage 1's step counters, both whole-column budgets, the slow ledger AND the     !
   !      soil-carbon seam maximum -- which is how it was caught, both stages reporting an          !
   !      identical worst gap to four significant figures.  ---------------------------------------!
   call meds_run_free(h)
   h = meds_run_open(to_c(CFG_DERIVED), int(len(CFG_DERIVED), c_int), 0_c_int)
   call check(h > 0_c_int, 'a second run opens after the first is freed')
   call check(meds_run_istep(h) == 0_c_int,                                                        &
              'a reopened registry slot starts at step 0, not the previous run''s count')
   call check(meds_run_is_done(h) == 0_c_int, 'and is not already done')
   status = meds_run_step(h)
   call check(meds_run_istep(h) == 1_c_int, 'the second run counts its own steps from 1')

   !=== 4. State is readable, and readable AFTER finalize (that is why free is separate). ====!
   agb   = meds_run_total_agb(h)
   lai   = meds_run_total_lai(h)
   soilc = meds_run_soil_carbon(h)
   call check(agb >= 0.0_wp .and. lai >= 0.0_wp, 'site aggregates read back non-negative')
   call check(soilc >= 0.0_wp, 'site soil carbon reads back non-negative')
   call check(meds_run_n_cohort(h) >= 0_c_int, 'the cohort count is readable')

   status = meds_run_finalize(h)
   call check(status == 0_c_int, 'finalize reports no area-conservation failure')
   call check(meds_run_total_agb(h) == agb, 'the state is still readable after finalize')
   call meds_run_free(h)
   call check(meds_run_total_agb(h) == 0.0_wp, 'a freed handle no longer reads state')

   call delete_file(CFG_DERIVED)
   print '(a)', 'test_capi_run: ALL PASSED'

contains

   !----- Prepend an override block to the shipped spin-up config. `meds_toml` returns the FIRST !
   !      match for a key, so these win over the same keys later in the file. What is turned off: !
   !      the fast loop and live met forcing (the forcing netCDF is untracked), the output        !
   !      streams and the state checkpoint (a test must not litter the source tree), and the      !
   !      window is cut to three days from bare ground.  -----------------------------------------!
   subroutine write_derived_config()
      integer :: uin, uout, ios
      character(len=512) :: line
      open(newunit=uout, file=CFG_DERIVED, status='replace', action='write')
      write(uout,'(a)') '# GENERATED by test_capi_run -- deleted when the test passes.'
      write(uout,'(a)') '[run]'
      write(uout,'(a)') 'start_time = "2024-07-01"'
      write(uout,'(a)') 'end_time   = "2024-07-04"'
      write(uout,'(a)') 'dt_slow    = "1d"'
      write(uout,'(a)') '[fast]'
      write(uout,'(a)') 'fast_biophysics_on = false'
      write(uout,'(a)') '[forcing]'
      write(uout,'(a)') 'forcing_on = false'
      write(uout,'(a)') '[init]'
      write(uout,'(a)') 'init_mode = 0'
      write(uout,'(a)') '[io]'
      write(uout,'(a)') 'write_state = false'
      write(uout,'(a)') '[output]'
      write(uout,'(a)') 'enabled = false'
      open(newunit=uin, file=CFG_SPINUP, status='old', action='read')
      do
         read(uin,'(a)',iostat=ios) line
         if (ios /= 0) exit
         write(uout,'(a)') trim(line)
      end do
      close(uin)
      close(uout)
   end subroutine write_derived_config

   subroutine delete_file(path)
      character(len=*), intent(in) :: path
      integer :: u
      logical :: there
      inquire(file=path, exist=there)
      if (.not. there) return
      open(newunit=u, file=path, status='old')
      close(u, status='delete')
   end subroutine delete_file

   !----- A NUL-free c_char array; the shim takes (buffer, length), not a C string. ----------!
   function to_c(s) result(a)
      character(len=*), intent(in)        :: s
      character(kind=c_char), allocatable :: a(:)
      integer :: i
      allocate(a(len(s)))
      do i = 1, len(s)
         a(i) = s(i:i)
      end do
   end function to_c

end program test_capi_run
