!==========================================================================================!
! meds_capi_run -- the ISO_C_BINDING shim that exposes the FULL COUPLED MODEL (`meds.model`).  !
!                                                                                          !
! The other three shims expose pieces: `meds_capi_leaf` and `meds_capi_phenology` are stateless !
! kernels, and `meds_capi_demography` drives the slow loop alone (no biophysics, no forcing, no !
! output streams). This one drives what the `meds_main` EXECUTABLE drives -- coupled fast + slow, !
! live met forcing, netCDF output, both conservation ledgers -- by holding a `meds_run_t` from    !
! `meds_driver` and handing the caller the OPEN / STEP / FINALIZE seam.                            !
!                                                                                          !
! The caller owns the time loop. That is the whole point: `examples/example_biophysics` used to    !
! be a bash script that exec'd the binary twice and then ran three plotting scripts over the        !
! netCDF it left behind; it is now one Python process that opens a run, steps it, reads the site    !
! state AS IT GOES, and plots. Nothing here re-implements the model -- `driver_step` is the same    !
! call the executable makes, so a Python-driven run and a `meds_main` run of the same config        !
! produce byte-identical output. The example's README records that comparison.                       !
!                                                                                          !
! Conventions match the other shims: `site_t` and friends never cross the boundary. A run is an     !
! opaque integer handle into a module-`save` registry; scalars return by value; arrays COPY into a  !
! caller-allocated buffer. The bind(c) entry points are also PUBLIC Fortran so `test_capi_run`      !
! can call them -- an untestable shim is exactly how `meds_capi_demography` came to stop compiling  !
! without anyone noticing.                                                                            !
!                                                                                          !
! `driver_step` RETURNS a NaN status rather than `error stop`ping, which is what makes this safe to  !
! call from an interpreter: a bad state surfaces as a Python exception, not a dead process.           !
!==========================================================================================!
module meds_capi_run
   use iso_c_binding
   use meds_kinds,             only : wp, ik
   use meds_driver,            only : meds_run_t, driver_open, driver_step, driver_finalize,   &
                                      driver_free, driver_done
   use meds_diagnostic_reduce, only : total_agb, total_lai, total_nplant, total_basal_area,    &
                                      count_cohorts
   implicit none
   private

   public :: meds_run_open, meds_run_step, meds_run_is_done, meds_run_finalize, meds_run_free
   public :: meds_run_year, meds_run_month, meds_run_day, meds_run_istep, meds_run_iyear
   public :: meds_run_n_patch, meds_run_n_cohort
   public :: meds_run_total_agb, meds_run_total_lai, meds_run_total_nplant
   public :: meds_run_total_basal_area, meds_run_soil_carbon, meds_run_get_real, meds_run_get_int

   !----- Small registry: a handful of concurrent runs is plenty (the example opens two, one per   !
   !      stage, and could hold both at once). Fixed-size and module-`save` for the same reason    !
   !      meds_capi_demography's is: the handle has to survive across ctypes calls.                !
   integer, parameter :: MAXR = 4
   type(meds_run_t), target, save :: g_run(MAXR)
   logical,                  save :: run_used(MAXR) = .false.

contains

   !----- Build a Fortran string from an explicit-length C char buffer (no strlen needed). -------!
   pure function f_string(cbuf, n) result(s)
      character(kind=c_char), intent(in) :: cbuf(*)
      integer(c_int),         intent(in) :: n
      character(len=n) :: s
      integer :: i
      do i = 1, n
         s(i:i) = cbuf(i)
      end do
   end function f_string

   !----- Guard every accessor: a handle the caller never opened (or already freed) must return a  !
   !      sentinel, not index a stale registry slot.  --------------------------------------------!
   pure logical function live(h)
      integer(c_int), intent(in) :: h
      live = (h >= 1_c_int .and. h <= int(MAXR, c_int))
      if (live) live = run_used(h)
   end function live

   !=======================================================================================!
   !  Lifecycle                                                                             !
   !=======================================================================================!

   !----- Open a run from a TOML config -> opaque handle (>=1); -1 registry full, -2 open failed. -!
   !      `verbose` /= 0 keeps the progress lines meds_main prints on stdout.                      !
   function meds_run_open(path, path_len, verbose) result(h) bind(c, name="meds_run_open")
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value,  intent(in) :: path_len, verbose
      integer(c_int)                     :: h
      integer :: i
      logical :: ok
      h = -1_c_int
      do i = 1, MAXR
         if (.not. run_used(i)) then ; h = int(i, c_int) ; exit ; end if
      end do
      if (h < 0_c_int) return
      call driver_open(f_string(path, path_len), g_run(h), ok, verbose=(verbose /= 0_c_int))
      if (.not. ok) then ; h = -2_c_int ; return ; end if
      run_used(h) = .true.
   end function meds_run_open

   !----- Advance ONE slow step. Returns the driver status: 0 stepped, 1 already finished, 2 NaN. -!
   function meds_run_step(h) result(status) bind(c, name="meds_run_step")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: status
      integer(ik) :: st
      if (.not. live(h)) then ; status = -1_c_int ; return ; end if
      call driver_step(g_run(h), st)
      status = int(st, c_int)
   end function meds_run_step

   !----- 1 when the calendar has reached end_time. The loop predicate, so the CALLER writes the   !
   !      `while` -- Python owning the time loop is the point of this shim.  ---------------------!
   function meds_run_is_done(h) result(d) bind(c, name="meds_run_is_done")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: d
      if (.not. live(h)) then ; d = 1_c_int ; return ; end if
      d = merge(1_c_int, 0_c_int, driver_done(g_run(h)))
   end function meds_run_is_done

   !----- Terminal checkpoint + conservation reports + close the streams. The state stays readable !
   !      afterwards (that is why free is separate).  --------------------------------------------!
   function meds_run_finalize(h) result(status) bind(c, name="meds_run_finalize")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: status
      integer(ik) :: st
      if (.not. live(h)) then ; status = -1_c_int ; return ; end if
      call driver_finalize(g_run(h), st)
      status = int(st, c_int)
   end function meds_run_finalize

   subroutine meds_run_free(h) bind(c, name="meds_run_free")
      integer(c_int), value, intent(in) :: h
      if (.not. live(h)) return
      call driver_free(g_run(h))
      run_used(h) = .false.
   end subroutine meds_run_free

   !=======================================================================================!
   !  Where the calendar is                                                                 !
   !=======================================================================================!
   function meds_run_year(h) result(v) bind(c, name="meds_run_year")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int ; if (live(h)) v = int(g_run(h)%now%year, c_int)
   end function meds_run_year

   function meds_run_month(h) result(v) bind(c, name="meds_run_month")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int ; if (live(h)) v = int(g_run(h)%now%month, c_int)
   end function meds_run_month

   function meds_run_day(h) result(v) bind(c, name="meds_run_day")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int ; if (live(h)) v = int(g_run(h)%now%day, c_int)
   end function meds_run_day

   function meds_run_istep(h) result(v) bind(c, name="meds_run_istep")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int ; if (live(h)) v = int(g_run(h)%istep, c_int)
   end function meds_run_istep

   function meds_run_iyear(h) result(v) bind(c, name="meds_run_iyear")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int ; if (live(h)) v = int(g_run(h)%iyear, c_int)
   end function meds_run_iyear

   !=======================================================================================!
   !  Site aggregates -- the same reducers print_summary uses, so a Python progress line and !
   !  the executable's own summary line report the same numbers.                             !
   !=======================================================================================!
   function meds_run_n_patch(h) result(v) bind(c, name="meds_run_n_patch")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int ; if (live(h)) v = int(g_run(h)%site%patch%n, c_int)
   end function meds_run_n_patch

   function meds_run_n_cohort(h) result(v) bind(c, name="meds_run_n_cohort")
      integer(c_int), value, intent(in) :: h
      integer(c_int)                    :: v
      v = -1_c_int ; if (live(h)) v = int(count_cohorts(g_run(h)%site), c_int)
   end function meds_run_n_cohort

   function meds_run_total_agb(h) result(v) bind(c, name="meds_run_total_agb")
      integer(c_int), value, intent(in) :: h
      real(c_double)                    :: v
      v = 0.0_c_double ; if (live(h)) v = real(total_agb(g_run(h)%site), c_double)
   end function meds_run_total_agb

   function meds_run_total_lai(h) result(v) bind(c, name="meds_run_total_lai")
      integer(c_int), value, intent(in) :: h
      real(c_double)                    :: v
      v = 0.0_c_double ; if (live(h)) v = real(total_lai(g_run(h)%site), c_double)
   end function meds_run_total_lai

   function meds_run_total_nplant(h) result(v) bind(c, name="meds_run_total_nplant")
      integer(c_int), value, intent(in) :: h
      real(c_double)                    :: v
      v = 0.0_c_double ; if (live(h)) v = real(total_nplant(g_run(h)%site), c_double)
   end function meds_run_total_nplant

   function meds_run_total_basal_area(h) result(v) bind(c, name="meds_run_total_basal_area")
      integer(c_int), value, intent(in) :: h
      real(c_double)                    :: v
      v = 0.0_c_double ; if (live(h)) v = real(total_basal_area(g_run(h)%site), c_double)
   end function meds_run_total_basal_area

   !----- Site total soil carbon [kgC/m2], area-weighted over patches: all seven CENTURY pools.    !
   !      Identically zero unless [soil_carbon].soil_carbon_on -- which is exactly what makes it    !
   !      worth watching from the driver loop, since with the feature off litter is DISCARDED       !
   !      rather than stored, and this number stays flat at zero while the stand grows.             !
   function meds_run_soil_carbon(h) result(v) bind(c, name="meds_run_soil_carbon")
      integer(c_int), value, intent(in) :: h
      real(c_double)                    :: v
      real(wp)    :: tot
      integer(ik) :: ip
      v = 0.0_c_double
      if (.not. live(h)) return
      tot = 0.0_wp
      associate (p => g_run(h)%site%patch)
         do ip = 1_ik, p%n
            tot = tot + p%area(ip) * (p%soil_carbon(ip)%fast_grnd_carbon                        &
                                    + p%soil_carbon(ip)%fast_soil_carbon                        &
                                    + p%soil_carbon(ip)%struct_grnd_carbon                      &
                                    + p%soil_carbon(ip)%struct_soil_carbon                      &
                                    + p%soil_carbon(ip)%microbial_carbon                        &
                                    + p%soil_carbon(ip)%slow_carbon                             &
                                    + p%soil_carbon(ip)%passive_carbon)
         end do
      end associate
      v = real(tot, c_double)
   end function meds_run_soil_carbon

   !=======================================================================================!
   !  Copy-out per-cohort getters (caller allocates buf of length n_cohort). Field ids match !
   !  meds_capi_demography's, so the two shims read the same way from Python.                !
   !=======================================================================================!
   subroutine meds_run_get_real(h, field_id, buf) bind(c, name="meds_run_get_real")
      integer(c_int), value, intent(in)  :: h, field_id
      real(c_double),        intent(out) :: buf(*)
      integer(ik) :: n
      if (.not. live(h)) return
      associate (c => g_run(h)%site%cohort)
         n = c%n
         if (n < 1_ik) return
         select case (field_id)
         case (0_c_int) ; buf(1:n) = real(c%dbh(1:n),             c_double)
         case (1_c_int) ; buf(1:n) = real(c%height(1:n),          c_double)
         case (2_c_int) ; buf(1:n) = real(c%nplant(1:n),          c_double)
         case (3_c_int) ; buf(1:n) = real(c%agb(1:n),             c_double)
         case (4_c_int) ; buf(1:n) = real(c%leaf_area(1:n),       c_double)
         case (5_c_int) ; buf(1:n) = real(c%overtopping_lai(1:n), c_double)
         case (6_c_int) ; buf(1:n) = real(c%growth_avg(1:n),      c_double)
         case (7_c_int) ; buf(1:n) = real(c%wood_carbon(1:n),     c_double)
         end select
      end associate
   end subroutine meds_run_get_real

   subroutine meds_run_get_int(h, field_id, buf) bind(c, name="meds_run_get_int")
      integer(c_int), value, intent(in)  :: h, field_id
      integer(c_int),        intent(out) :: buf(*)
      integer(ik) :: n
      if (.not. live(h)) return
      associate (c => g_run(h)%site%cohort)
         n = c%n
         if (n < 1_ik) return
         select case (field_id)
         case (0_c_int) ; buf(1:n) = int(c%pft(1:n),         c_int)
         case (1_c_int) ; buf(1:n) = int(c%owner_patch(1:n), c_int)
         case (2_c_int) ; buf(1:n) = int(c%global_id(1:n),   c_int)
         end select
      end associate
   end subroutine meds_run_get_int

end module meds_capi_run
