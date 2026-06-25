!==========================================================================================!
! meds_toml -- a minimal TOML reader for MEDS run configuration.                           !
!                                                                                          !
! A full TOML library (toml-f) is not usable here: like netCDF-Fortran, its compiled `.mod` !
! is compiler-specific and unavailable for ifx/nvfortran on this machine. MEDS configs are   !
! a flat set of `key = value` lines under `[section]` headers, so we parse that subset        !
! directly: `# comments`, `[section]` headers, and scalar int/real/bool/string values plus    !
! `[a, b, c]` real arrays. Keys are stored fully qualified as `section.key`. The typed        !
! getters return a supplied default when the key is absent or unparseable, so a partial (or    !
! missing) config file simply leaves the defaults in place.                                   !
!==========================================================================================!
module meds_toml
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: toml_table_t, toml_parse_file
   public :: toml_has, toml_int, toml_real, toml_logical, toml_string, toml_real_array

   integer, parameter :: KEYLEN = 64, VALLEN = 256, MAXKEYS = 512

   type :: toml_table_t
      integer(ik)                     :: n = 0_ik
      character(len=KEYLEN), allocatable :: key(:)
      character(len=VALLEN), allocatable :: val(:)
   end type toml_table_t

contains

   !---------------------------------------------------------------------------------------!
   ! Parse a TOML file into the (section.key -> raw-value) table. ok=.false. if it cannot    !
   ! be opened (the caller then keeps all defaults).                                         !
   !---------------------------------------------------------------------------------------!
   subroutine toml_parse_file(path, t, ok)
      character(len=*),    intent(in)  :: path
      type(toml_table_t),  intent(out) :: t
      logical,             intent(out) :: ok
      integer                :: unit, ios, ieq, nb
      character(len=VALLEN)  :: line
      character(len=KEYLEN)  :: section
      character(len=VALLEN)  :: raw

      ok = .false.
      allocate(t%key(MAXKEYS), t%val(MAXKEYS))
      t%n = 0_ik
      section = ''

      open(newunit=unit, file=path, status='old', action='read', iostat=ios)
      if (ios /= 0) return
      ok = .true.

      do
         read(unit, '(a)', iostat=ios) raw
         if (ios /= 0) exit
         line = strip_comment(raw)
         line = adjustl(line)
         if (len_trim(line) == 0) cycle

         if (line(1:1) == '[') then                       ! [section] header
            nb = index(line, ']')
            if (nb > 2) section = adjustl(line(2:nb-1))
            cycle
         end if

         ieq = index(line, '=')
         if (ieq < 2) cycle                               ! not a key = value line
         if (t%n >= MAXKEYS) exit
         t%n = t%n + 1_ik
         if (len_trim(section) > 0) then
            t%key(t%n) = trim(section) // '.' // trim(adjustl(line(1:ieq-1)))
         else
            t%key(t%n) = trim(adjustl(line(1:ieq-1)))
         end if
         t%val(t%n) = adjustl(line(ieq+1:))
      end do
      close(unit)
   end subroutine toml_parse_file

   !----- Strip a trailing `# comment` that is not inside a quoted string. ----------------!
   pure function strip_comment(s) result(out)
      character(len=*), intent(in) :: s
      character(len=len(s))        :: out
      integer :: i, nq
      out = s
      nq = 0
      do i = 1, len(s)
         if (s(i:i) == '"') nq = nq + 1
         if (s(i:i) == '#' .and. mod(nq, 2) == 0) then
            out = s(1:i-1)
            return
         end if
      end do
   end function strip_comment

   !----- Index of a fully-qualified key, or 0 if absent. ---------------------------------!
   pure integer(ik) function find_key(t, key) result(idx)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      integer(ik) :: i
      idx = 0_ik
      do i = 1_ik, t%n
         if (trim(t%key(i)) == trim(key)) then
            idx = i
            return
         end if
      end do
   end function find_key

   logical function toml_has(t, key) result(yes)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      yes = find_key(t, key) > 0_ik
   end function toml_has

   integer(ik) function toml_int(t, key, default) result(v)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      integer(ik),        intent(in) :: default
      integer(ik) :: idx, ios
      v = default
      idx = find_key(t, key)
      if (idx > 0_ik) read(t%val(idx), *, iostat=ios) v
      if (idx > 0_ik .and. ios /= 0) v = default
   end function toml_int

   real(wp) function toml_real(t, key, default) result(v)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      real(wp),           intent(in) :: default
      integer(ik) :: idx, ios
      v = default
      idx = find_key(t, key)
      if (idx > 0_ik) read(t%val(idx), *, iostat=ios) v
      if (idx > 0_ik .and. ios /= 0) v = default
   end function toml_real

   logical function toml_logical(t, key, default) result(v)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      logical,            intent(in) :: default
      integer(ik) :: idx
      character(len=VALLEN) :: s
      v = default
      idx = find_key(t, key)
      if (idx == 0_ik) return
      s = adjustl(t%val(idx))
      if (s(1:4) == 'true')  v = .true.
      if (s(1:5) == 'false') v = .false.
   end function toml_logical

   !----- String value with surrounding double quotes stripped. ---------------------------!
   function toml_string(t, key, default) result(v)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key, default
      character(len=VALLEN)          :: v
      integer(ik) :: idx
      character(len=VALLEN) :: s
      v = default
      idx = find_key(t, key)
      if (idx == 0_ik) return
      s = adjustl(t%val(idx))
      if (len_trim(s) >= 2 .and. s(1:1) == '"') then
         v = s(2:index(s(2:), '"'))
      else
         v = trim(s)
      end if
   end function toml_string

   !----- Parse `[a, b, c]` into out(:) reals; nout = count read (0 if absent). ------------!
   subroutine toml_real_array(t, key, out, nout)
      type(toml_table_t), intent(in)  :: t
      character(len=*),   intent(in)  :: key
      real(wp),           intent(out) :: out(:)
      integer(ik),        intent(out) :: nout
      integer(ik) :: idx, i, lb, rb, ios
      character(len=VALLEN) :: s
      nout = 0_ik
      idx = find_key(t, key)
      if (idx == 0_ik) return
      s  = t%val(idx)
      lb = index(s, '[') ; rb = index(s, ']', back=.true.)
      if (lb == 0 .or. rb <= lb) return
      s = s(lb+1:rb-1)
      do i = 1, len_trim(s)                         ! commas -> spaces for list-directed read
         if (s(i:i) == ',') s(i:i) = ' '
      end do
      block
         real(wp) :: tmp(size(out))
         integer(ik) :: k
         tmp = 0.0_wp
         read(s, *, iostat=ios) (tmp(k), k = 1_ik, size(out, kind=ik))
         ! list-directed read stops at end-of-record; count the leading non-defaulted entries
         do k = 1_ik, size(out, kind=ik)
            out(k) = tmp(k)
         end do
      end block
      nout = count_tokens(s)
   end subroutine toml_real_array

   !----- Number of whitespace-separated tokens in a string. ------------------------------!
   pure integer(ik) function count_tokens(s) result(nt)
      character(len=*), intent(in) :: s
      integer(ik) :: i
      logical     :: in_tok
      nt = 0_ik ; in_tok = .false.
      do i = 1_ik, len_trim(s)
         if (s(i:i) /= ' ') then
            if (.not. in_tok) nt = nt + 1_ik
            in_tok = .true.
         else
            in_tok = .false.
         end if
      end do
   end function count_tokens

end module meds_toml
