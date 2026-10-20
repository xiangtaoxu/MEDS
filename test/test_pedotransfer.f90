! SPDX-License-Identifier: Apache-2.0
!----- The soil pedotransfer function (meds_soil_hydr_lib): Hodnett & Tomasella (2002)'s        !
!      worked example, Ottoni et al. (2019)'s Ks of it, coarse fragments, the residual-water    !
!      floor, and the elemental call over a profile against the call per layer, and CLM5's      !
!      thermal properties of a soil (soil_texture_to_thermal, in meds_therm_lib);               !
!      docs/science/soil_biophysics.md eqs. (6)-(9). -------------------------------------------!
program test_pedotransfer
   use meds_kinds,         only : wp
   use meds_constants,     only : grav_head
   use meds_soil_hydr_lib, only : soil_texture_t, soil_texture_to_hydraulics
   use meds_therm_lib,     only : soil_texture_to_thermal
   use meds_test_assert,   only : banner, check_close, check_true, test_report
   implicit none

   type(soil_texture_t)   :: example, stony, bare_sand, profile(3)
   real(wp)               :: theta_sat, theta_res, vg_alpha, vg_n, ksat
   real(wp)               :: theta_sat_stony, theta_res_stony, alpha_stony, n_stony, ksat_stony
   real(wp), dimension(3) :: theta_sat_p, theta_res_p, alpha_p, n_p, ksat_p
   real(wp)               :: theta_sat_k, theta_res_k, alpha_k, n_k, ksat_k
   real(wp)               :: k_solid, k_dry, dry_cvol
   logical                :: same
   integer                :: k

   call banner('soil pedotransfer function')

   !=== Hodnett & Tomasella (2002) Table 8's worked example, as Marthews et al. (2014, GMD 7:711)  !
   !    repeat it: alpha 0.2452 per kPa, n 1.3674, theta_sat 0.4993, theta_res 0.2134. ============!
   example = soil_texture_t(sand=39.2_wp, silt=24.2_wp, clay=36.7_wp, org_carbon=14.0_wp, cec=15.9_wp,  &
                            ph=5.85_wp, bulk_density=1.2_wp, coarse_frag=0.0_wp)
   call soil_texture_to_hydraulics(example, theta_sat, theta_res, vg_alpha, vg_n, ksat)
   call check_close('worked example: theta_sat', theta_sat, 0.4993_wp, 1.0e-4_wp)
   call check_close('worked example: theta_res', theta_res, 0.2134_wp, 1.0e-4_wp)
   call check_close('worked example: alpha per kPa', vg_alpha / (grav_head*1000.0_wp), 0.2452_wp, 1.0e-4_wp)
   call check_close('worked example: n', vg_n, 1.3674_wp, 1.0e-4_wp)

   !=== Ottoni et al. (2019) for the example: 19.31 m/day x (theta_sat - theta(-33 kPa))^1.948,    !
   !    with theta(-33 kPa) = 0.34409; 5.9343e-6 m/s (21.4 mm/h), computed apart from MEDS. ======!
   call check_close('worked example: Ks [m/s]', ksat, 5.934313e-6_wp, 1.0e-11_wp)

   !=== 30 % coarse fragments by volume at a fine-earth bulk density of 1.2: the curve holds 0.7 of !
   !    the fine earth's water; the fragments are 0.48624 of the mass, so Ks keeps 0.60870 of the  !
   !    fine earth's (Saxton & Rawls 2006 eqs. 19 and 22). The curve's shape does not change. =====!
   stony = example
   stony%coarse_frag = 30.0_wp
   call soil_texture_to_hydraulics(stony, theta_sat_stony, theta_res_stony, alpha_stony, n_stony, ksat_stony)
   call check_close('coarse fragments: theta_sat', theta_sat_stony, 0.7_wp*theta_sat, 1.0e-12_wp)
   call check_close('coarse fragments: theta_res', theta_res_stony, 0.7_wp*theta_res, 1.0e-12_wp)
   call check_close('coarse fragments: Ks ratio', ksat_stony/ksat, 0.6086956522_wp, 1.0e-9_wp)
   call check_true('coarse fragments: alpha and n unchanged', alpha_stony == vg_alpha .and. n_stony == vg_n)

   !=== A sandy, alkaline soil of little CEC: the residual-water regression is negative there   !
   !    (-0.0086), and the floor holds it at zero. ================================================!
   bare_sand = soil_texture_t(sand=98.0_wp, silt=1.0_wp, clay=1.0_wp, org_carbon=1.0_wp, cec=0.5_wp,  &
                              ph=9.5_wp, bulk_density=1.6_wp, coarse_frag=0.0_wp)
   call soil_texture_to_hydraulics(bare_sand, theta_sat_k, theta_res_k, alpha_k, n_k, ksat_k)
   call check_close('sandy alkaline soil: theta_res floored at 0', theta_res_k, 0.0_wp, 0.0_wp)
   call check_true('sandy alkaline soil: a curve that drains (theta_sat > theta_res, n > 1, Ks > 0)',  &
                   theta_sat_k > theta_res_k .and. n_k > 1.0_wp .and. ksat_k > 0.0_wp)

   !=== One elemental call over a profile gives exactly the per-layer results. ==================!
   profile = [example, stony, bare_sand]
   call soil_texture_to_hydraulics(profile, theta_sat_p, theta_res_p, alpha_p, n_p, ksat_p)
   same = .true.
   do k = 1, 3
      call soil_texture_to_hydraulics(profile(k), theta_sat_k, theta_res_k, alpha_k, n_k, ksat_k)
      same = same .and. theta_sat_p(k) == theta_sat_k .and. theta_res_p(k) == theta_res_k .and.     &
             alpha_p(k) == alpha_k .and. n_p(k) == n_k .and. ksat_p(k) == ksat_k
   end do
   call check_true('elemental call over a profile equals the per-layer calls', same)

   !=== CLM5's mineral soil at Palo Verde's column (sand 33 %, clay 40 %, porosity 0.4456), computed  !
   !    apart from MEDS: solids (8.80 x 33 + 2.92 x 40) / 73 = 5.57808 W/m/K; dry bulk density 1496.88 !
   !    kg/m3, so dry soil (0.135 x 1496.88 + 64.7) / (2700 - 0.947 x 1496.88) = 0.208022 W/m/K; and   !
   !    the solids' heat capacity (2.128 x 33 + 2.385 x 40) / 73 x 1e6 x 0.5544 = 1.257835e6 J/m3/K. ==!
   call soil_texture_to_thermal(33.0_wp, 40.0_wp, 0.4456_wp, k_solid, k_dry, dry_cvol)
   call check_close('CLM5 solids conductivity', k_solid, 5.5780822_wp, 1.0e-6_wp)
   call check_close('CLM5 dry conductivity', k_dry, 0.20802202_wp, 1.0e-7_wp)
   call check_close('CLM5 solids heat capacity', dry_cvol, 1.2578349e6_wp, 1.0_wp)

   call test_report('test_pedotransfer')
end program test_pedotransfer
