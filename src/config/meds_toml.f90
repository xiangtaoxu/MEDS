!==========================================================================================!
! meds_toml -- a minimal TOML reader for MEDS run configuration.                           !
!                                                                                          !
! A full TOML library (toml-f) is not usable here: like netCDF-Fortran, its compiled `.mod` !
! is compiler-specific and unavailable for ifx/nvfortran on this machine. MEDS configs are   !
! a flat set of `key = value` lines under `[section]` headers, so we parse that subset        !
! directly: `# comments`, `[section]` headers, and scalar int/real/bool/string values plus    !
! `[a, b, c]` real arrays. Keys are stored fully qualified as `section.key`. The typed        !
! getters return the supplied default ONLY when the key is ABSENT; a key that is PRESENT but    !
! unparseable is a HARD ERROR (error stop with the offending key + raw text) -- silently zero-  !
! filling malformed values would let a typo'd REQUIRED parameter pass validation (see the       !
! meds_config_io 'every parameter required, no defaults' contract).                             !
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

   !----- One place for the 'present-but-malformed value' hard error: prints the offending key  !
   !      + raw text then aborts (mirrors the missing-key error-stop style in meds_config_io). --!
   subroutine toml_parse_error(key, raw, what)
      character(len=*), intent(in) :: key, raw, what
      write(*,'(5a)') ' meds_toml: cannot parse ', trim(what), ' for key "', trim(key), '"'
      write(*,'(3a)') '   raw value: "', trim(raw), '"'
      error stop 'meds_toml: malformed configuration value'
   end subroutine toml_parse_error

   integer(ik) function toml_int(t, key, default) result(v)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      integer(ik),        intent(in) :: default
      integer(ik) :: idx, ios
      v = default
      idx = find_key(t, key)
      if (idx == 0_ik) return                             ! absent -> default
      read(t%val(idx), *, iostat=ios) v
      if (ios /= 0) call toml_parse_error(key, t%val(idx), 'integer')
   end function toml_int

   real(wp) function toml_real(t, key, default) result(v)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      real(wp),           intent(in) :: default
      integer(ik) :: idx, ios
      v = default
      idx = find_key(t, key)
      if (idx == 0_ik) return                             ! absent -> default
      read(t%val(idx), *, iostat=ios) v
      if (ios /= 0) call toml_parse_error(key, t%val(idx), 'real')
   end function toml_real

   logical function toml_logical(t, key, default) result(v)
      type(toml_table_t), intent(in) :: t
      character(len=*),   intent(in) :: key
      logical,            intent(in) :: default
      integer(ik) :: idx
      character(len=VALLEN) :: s
      v = default
      idx = find_key(t, key)
      if (idx == 0_ik) return                             ! absent -> default
      s = adjustl(t%val(idx))
      !----- Exact match (documented lowercase true/false + common Fortran/case variants).      !
      !      Kills the old fixed-position prefix bug ('trueish' -> .true., '.true.' -> default). !
      select case (trim(s))
      case ('true', '.true.', 'True', 'TRUE')
         v = .true.
      case ('false', '.false.', 'False', 'FALSE')
         v = .false.
      case default
         call toml_parse_error(key, t%val(idx), 'logical')
      end select
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
      if (idx == 0_ik) return                       ! absent -> nout = 0 (caller records missing key)
      s  = t%val(idx)
      lb = index(s, '[') ; rb = index(s, ']', back=.true.)
      if (lb == 0 .or. rb <= lb) call toml_parse_error(key, t%val(idx), 'real array (expected [..] literal)')
      s = s(lb+1:rb-1)
      do i = 1, len_trim(s)                         ! commas -> spaces for list-directed read
         if (s(i:i) == ',') s(i:i) = ' '
      end do
      block
         real(wp)    :: tmp(size(out))
         integer(ik) :: k, ntok
         ntok = count_tokens(s)
         if (ntok == 0_ik) return                   ! empty [] -> nout = 0
         !----- Bound by the receiving buffer BEFORE reading (prevents an out-of-bounds read on  !
         !      an over-length list, e.g. > MAXPFT trait entries into a MAXPFT-sized buffer). ----!
         if (ntok > size(out, kind=ik)) call toml_parse_error(key, t%val(idx), 'real array (too many entries)')
         tmp = 0.0_wp
         read(s, *, iostat=ios) (tmp(k), k = 1_ik, ntok)
         if (ios /= 0) call toml_parse_error(key, t%val(idx), 'real array')   ! any unparseable token -> hard error
         out(1:ntok) = tmp(1:ntok)
         nout = ntok                                ! TRUE parsed count (never over-reported)
      end block
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
