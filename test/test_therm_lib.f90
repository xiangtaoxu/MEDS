! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_therm_lib -- the thermodynamics the stores share (meds_therm_lib).                  !
!   1. TEMPERATURE RESPONSE (arrhenius_scale, peaked_arrhenius_scale,                      !
!      kattge_knorr_entropy; docs/science/plant_respiration.md eq. (5),                    !
!      docs/science/leaf_gas_exchange.md eq. (4)): each form is k25 at 25 C; Arrhenius     !
!      rises with temperature and the peaked form has an interior optimum; a warm-grown    !
!      plant's entropy shift keeps more capacity in heat and moves the optimum up.         !
!==========================================================================================!
program test_therm_lib
   use meds_test_assert, only : check, check_true, test_report
   use meds_kinds,       only : wp, ik
   use meds_constants,   only : t_kelvin
   use meds_therm_lib,   only : arrhenius_scale, peaked_arrhenius_scale, kattge_knorr_entropy
   implicit none

   call test_temperature_response()

   call test_report('test_therm_lib')

contains

   !----- 1. The temperature response of enzyme rates, at Vcmax's activation, deactivation and     !
   !      entropy (60 at 25 C; Ea 65330 J/mol, Hd 200000 J/mol, dS 650 J/mol/K). Kattge & Knorr    !
   !      (2007)'s acclimation shifts dS by the growth temperature; what it is for is moving where  !
   !      the response peaks, measured here from the function by a scan. --------------------------!
   subroutine test_temperature_response()
      real(wp), parameter :: VC25 = 60.0_wp, EA = 65330.0_wp, HD = 200000.0_wp, DS = 650.0_wp
      real(wp), parameter :: T_COLD = t_kelvin + 10.0_wp, T_WARM = t_kelvin + 30.0_wp
      real(wp) :: k25, k32, k45, ds_cold, ds_warm, k_cold, k_warm, topt_cold, topt_warm
      print '(a)', 'test_temperature_response:'
      call check('Arrhenius is k25 at 25 C', arrhenius_scale(40.49_wp, 79430.0_wp, t_kelvin + 25.0_wp), 40.49_wp,  &
                 1.0e-12_wp * 40.49_wp)
      call check_true('Arrhenius increases with temperature',                                                  &
                      arrhenius_scale(40.49_wp, 79430.0_wp, t_kelvin + 5.0_wp) <                              &
                      arrhenius_scale(40.49_wp, 79430.0_wp, t_kelvin + 35.0_wp))
      k25 = peaked_arrhenius_scale(VC25, EA, HD, DS, t_kelvin + 25.0_wp)
      k32 = peaked_arrhenius_scale(VC25, EA, HD, DS, t_kelvin + 32.0_wp)
      k45 = peaked_arrhenius_scale(VC25, EA, HD, DS, t_kelvin + 45.0_wp)
      call check('peaked is k25 at 25 C', k25, VC25, 1.0e-12_wp * VC25)
      call check_true('peaked has an interior thermal optimum', k32 > k25 .and. k32 > k45, k32)
      !----- At a hot leaf temperature the warm-acclimated plant keeps more capacity than the      !
      !      cold-acclimated one, and its optimum sits higher. ------------------------------------!
      ds_cold = kattge_knorr_entropy(668.39_wp, 1.07_wp, T_COLD)
      ds_warm = kattge_knorr_entropy(668.39_wp, 1.07_wp, T_WARM)
      k_cold  = peaked_arrhenius_scale(VC25, EA, HD, ds_cold, t_kelvin + 38.0_wp)
      k_warm  = peaked_arrhenius_scale(VC25, EA, HD, ds_warm, t_kelvin + 38.0_wp)
      call check_true('warm acclimation retains more Vcmax at 38 C than cold acclimation', k_warm > k_cold,   &
                      k_warm - k_cold)
      topt_cold = peak_temperature(ds_cold)
      topt_warm = peak_temperature(ds_warm)
      call check_true('the thermal optimum moves up with growth temperature (> 1 K over 20 K)',                &
                      topt_warm > topt_cold + 1.0_wp, topt_warm - topt_cold)
   end subroutine test_temperature_response

   !----- The peak of the peaked response [K] by a scan at 0.05 K, so the optimum is measured from   !
   !      the function rather than derived from a formula the test would check against itself. ---!
   pure real(wp) function peak_temperature(ds) result(t_opt)
      real(wp), intent(in) :: ds   !< [J/mol/K] the entropy term
      real(wp)    :: t, k, k_best
      integer(ik) :: j
      k_best = -1.0_wp
      t_opt  = 0.0_wp
      do j = 0_ik, 600_ik
         t = t_kelvin + 5.0_wp + 0.05_wp * real(j, wp)
         k = peaked_arrhenius_scale(60.0_wp, 65330.0_wp, 200000.0_wp, ds, t)
         if (k > k_best) then
            k_best = k
            t_opt  = t
         end if
      end do
   end function peak_temperature

end program test_therm_lib
