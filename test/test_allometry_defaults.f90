! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_allometry_defaults -- the shipped meds_config_pft.toml [allometry] block must equal !
! the meds_allometry_lib initializers; the biomass law must be Chave et al. (2014) eq. 4   !
! in carbon, and the leaf-area scale ED2's BAAD fit divided by C2B.                        !
!                                                                                          !
! The initializers exist for the Python shared library, which reaches the allometry without a    !
! config; a model run takes every value from the PFT config instead. The two copies can drift     !
! apart silently, and nothing else compares them.                                                 !
!                                                                                          !
! Both fits are dry mass or its equivalent, and MEDS keeps carbon: Chave's 0.0673 kg dry mass     !
! becomes 0.0673/C2B kgC, and ED2's size2bl uses c14f15_bl_xx(1)/(SLA*C2B). Taken without the      !
! division, AGB is dry mass labelled as carbon and leaf area is twice ED2's; this test guards     !
! against both.                                                                                   !
!==========================================================================================!
program test_allometry_defaults
   use meds_kinds,         only : wp, ik
   use meds_toml,          only : toml_table_t, toml_parse_file, toml_has, toml_real, toml_string, toml_real_array
   use meds_allometry_lib, only : allometry_t, HEIGHT_POWER
   use meds_test_support,  only : banner, check, check_close
   implicit none

   real(wp), parameter :: C2B              = 2.0_wp          ! ED2 C2B, carbon -> biomass
   real(wp), parameter :: CHAVE14_SCALE    = 0.0673_wp       ! Chave et al. (2014) eq. 4 [kg dry mass]
   real(wp), parameter :: CHAVE14_EXPONENT = 0.976_wp        ! Chave et al. (2014) eq. 4
   real(wp), parameter :: ED2_BL_XX_SCALE  = 0.46769540_wp   ! ED2 c14f15_bl_xx(1)
   real(wp), parameter :: RTOL             = 1.0e-12_wp
   type(toml_table_t)  :: t
   character(len=512)  :: arg
   logical             :: ok
   type(allometry_t)   :: d          ! a record at its initializers

   call banner('shipped allometry defaults')

   !----- CMake passes the SOURCE-tree path; ctest runs from the build tree. -------------!
   call get_command_argument(1, arg)
   call check(len_trim(arg) > 0, 'the meds_config_pft.toml path was passed on the command line')
   call toml_parse_file(trim(arg), t, ok)
   call check(ok, 'the shipped meds_config_pft.toml parses')

   !----- The Fortran initializers: Chave 2014 in carbon, and ED2's leaf fit over C2B. ----!
   call check_close(d%agb_c1, CHAVE14_SCALE / C2B, RTOL, 'agb_c1 initializer = Chave 2014 eq. 4 scale / C2B')
   call check_close(d%agb_c2, CHAVE14_EXPONENT, RTOL, 'agb_c2 initializer = Chave 2014 eq. 4 exponent')
   call check_close(d%lai_b1, ED2_BL_XX_SCALE / C2B, RTOL, 'lai_b1 initializer = ED2 c14f15_bl_xx(1)/C2B')

   !----- The shipped config equals the initializers, key by key. -----------------------!
   call same('b1Ht', d%b1Ht)
   call same('b2Ht', d%b2Ht)
   call same('agb_c1', d%agb_c1)
   call same('agb_c2', d%agb_c2)
   call same('lai_b1', d%lai_b1)
   call same('lai_b2', d%lai_b2)
   !----- The keys the shipped config lists per PFT: every entry equals the initializer. --!
   call same_per_pft('hgt_max',          d%hgt_max)
   call same_per_pft('wai_b1',           d%wai_b1)
   call same_per_pft('wai_b2',           d%wai_b2)
   call same_per_pft('sapwood_area_b1',  d%sapwood_area_b1)
   call same_per_pft('sapwood_area_b2',  d%sapwood_area_b2)
   call same_per_pft('aboveground_frac', d%aboveground_frac)
   !----- The log-quadratic initializers are Chave et al. (2014) eq. 6a. -----------------!
   call check_close(d%chave_a, 0.893_wp, RTOL, 'chave_a initializer = Chave 2014 eq. 6a')
   call check_close(d%chave_b, 0.760_wp, RTOL, 'chave_b initializer = Chave 2014 eq. 6a')
   call check_close(d%chave_c, -0.0340_wp, RTOL, 'chave_c initializer = Chave 2014 eq. 6a')
   call check_close(d%chave_rse, 0.243_wp, RTOL, 'chave_rse initializer = Chave 2014 eq. 6a')

   !----- The height curve: the shipped config and the initializer are both the power law. ---!
   call check(trim(toml_string(t, 'allometry.height_allometry', '')) == 'power',                    &
              'meds_config_pft.toml sets [allometry].height_allometry = "power"')
   call check(d%height_allometry == HEIGHT_POWER, 'the height-curve initializer is the power law')

   print '(a)', 'test_allometry_defaults: ALL PASSED'

contains

   subroutine same(key, initializer)
      character(len=*), intent(in) :: key
      real(wp),         intent(in) :: initializer
      call check(toml_has(t, 'allometry.'//key), 'meds_config_pft.toml sets [allometry].'//key)
      call check_close(toml_real(t, 'allometry.'//key, -huge(1.0_wp)), initializer, RTOL,         &
                       'meds_config_pft.toml [allometry].'//key//' equals the meds_allometry_lib initializer')
   end subroutine same

   subroutine same_per_pft(key, initializer)
      character(len=*), intent(in) :: key
      real(wp),         intent(in) :: initializer
      real(wp)    :: buf(64)
      integer(ik) :: n
      call toml_real_array(t, 'allometry.'//key, buf, n)
      call check(n > 0, 'meds_config_pft.toml lists [allometry].'//key//' per PFT')
      call check(all(abs(buf(1:n) - initializer) <= RTOL * abs(initializer)),                         &
                 'meds_config_pft.toml [allometry].'//key//' equals the meds_allometry_lib initializer for every PFT')
   end subroutine same_per_pft

end program test_allometry_defaults
