! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_allometry_defaults -- the shipped meds_config_pft.toml [allometry] block must equal the  !
! meds_allometry initializers; the biomass law must be Chave et al. (2014) eq. 4 in carbon, and  !
! the leaf-area scale ED2's BAAD fit divided by C2B.                                              !
!                                                                                          !
! The initializers exist for the Python shared library, which reaches the allometry without a    !
! config; a model run takes every value from the PFT config instead. The two copies can drift     !
! apart silently, and nothing else compares them.                                                 !
!                                                                                          !
! Both fits are dry mass or its equivalent, and MEDS keeps carbon: Chave's 0.0673 kg dry mass     !
! becomes 0.0673/C2B kgC, and ED2's size2bl uses c14f15_bl_xx(1)/(SLA*C2B). Taken without the      !
! division, AGB is dry mass labelled as carbon and leaf area is twice ED2's; the shipped values    !
! carried exactly that until the fix.                                                              !
!==========================================================================================!
program test_allometry_defaults
   use meds_kinds,        only : wp
   use meds_toml,         only : toml_table_t, toml_parse_file, toml_has, toml_real, toml_string
   use meds_allometry,    only : b1Ht, b2Ht, agb_c1, agb_c2, ca_b1, ca_b2, lai_b1, lai_b2,      &
                                 light_ext, height_form, gmm_a, gmm_b, gmm_k, HEIGHT_POWER
   use meds_test_support, only : banner, check, check_close
   implicit none

   real(wp), parameter :: c2b              = 2.0_wp          ! ED2 C2B, carbon -> biomass
   real(wp), parameter :: chave14_scale    = 0.0673_wp       ! Chave et al. (2014) eq. 4 [kg dry mass]
   real(wp), parameter :: chave14_exponent = 0.976_wp        ! Chave et al. (2014) eq. 4
   real(wp), parameter :: ed2_bl_xx_scale  = 0.46769540_wp   ! ED2 c14f15_bl_xx(1)
   real(wp), parameter :: rtol             = 1.0e-12_wp
   type(toml_table_t)  :: t
   character(len=512)  :: arg
   logical             :: ok

   call banner('shipped allometry defaults')

   !----- CMake passes the SOURCE-tree path; ctest runs from the build tree. -------------!
   call get_command_argument(1, arg)
   call check(len_trim(arg) > 0, 'the meds_config_pft.toml path was passed on the command line')
   call toml_parse_file(trim(arg), t, ok)
   call check(ok, 'the shipped meds_config_pft.toml parses')

   !----- The Fortran initializers: Chave 2014 in carbon, and ED2's leaf fit over C2B. ----!
   call check_close(agb_c1, chave14_scale / c2b, rtol, 'agb_c1 initializer = Chave 2014 eq. 4 scale / C2B')
   call check_close(agb_c2, chave14_exponent, rtol, 'agb_c2 initializer = Chave 2014 eq. 4 exponent')
   call check_close(lai_b1, ed2_bl_xx_scale / c2b, rtol, 'lai_b1 initializer = ED2 c14f15_bl_xx(1)/C2B')

   !----- The shipped config equals the initializers, key by key. -----------------------!
   call same('b1Ht', b1Ht)
   call same('b2Ht', b2Ht)
   call same('agb_c1', agb_c1)
   call same('agb_c2', agb_c2)
   call same('ca_b1', ca_b1)
   call same('ca_b2', ca_b2)
   call same('lai_b1', lai_b1)
   call same('lai_b2', lai_b2)
   call same('light_ext', light_ext)

   !----- The height curve: the shipped config and the initializer are both the power law, and  !
   !      the gMM initializers are Cano et al. (2019) eq. 7, the values the config lists. -------!
   call check(trim(toml_string(t, 'allometry.height_allometry', '')) == 'power',                    &
              'meds_config_pft.toml sets [allometry].height_allometry = "power"')
   call check(height_form == HEIGHT_POWER, 'the height-curve initializer is the power law')
   call check_close(gmm_a, 58.0_wp, rtol, 'gmm_a initializer = Cano et al. (2019) eq. 7')
   call check_close(gmm_b, 0.73_wp, rtol, 'gmm_b initializer = Cano et al. (2019) eq. 7')
   call check_close(gmm_k, 21.8_wp, rtol, 'gmm_k initializer = Cano et al. (2019) eq. 7')

   print '(a)', 'test_allometry_defaults: ALL PASSED'

contains

   subroutine same(key, initializer)
      character(len=*), intent(in) :: key
      real(wp),         intent(in) :: initializer
      call check(toml_has(t, 'allometry.'//key), 'meds_config_pft.toml sets [allometry].'//key)
      call check_close(toml_real(t, 'allometry.'//key, -huge(1.0_wp)), initializer, rtol,         &
                       'meds_config_pft.toml [allometry].'//key//' equals the meds_allometry initializer')
   end subroutine same

end program test_allometry_defaults
