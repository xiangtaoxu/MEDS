! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! test_pft_optics_config -- the [pft] canopy optics reach the two-stream's per-PFT table,    !
! and DIFFERENTIATE the PFTs.                                                                !
!                                                                                          !
! These fourteen traits were literals broadcast over every PFT, which meant two plant        !
! functional types could not differ in how they intercept or scatter light -- in a model      !
! whose two-stream solver exists to resolve that competition. Moving them to the table is      !
! only worth anything if the values arrive per PFT, and that wiring fails silently: a          !
! mis-threaded index leaves every PFT on the first one's optics, the radiation still solves,   !
! the energy budget still closes, and the canopy is quietly one optical type again. So this    !
! asserts the values arrive AND that PFTs differ, not merely that the code compiles.           !
!==========================================================================================!
program test_pft_optics_config
   use meds_kinds,          only : wp, ik
   use meds_config,         only : meds_config_t
   use meds_config_io,      only : write_pft_params_csv
   use meds_plant_types,    only : hydro_params_table_t
   use meds_fast_types,     only : apply_hydraulics_config
   use meds_hydr_lib,       only : plc_retained
   use meds_fast_dynamics,  only : fast_context_t, build_fast_context
   use meds_canopy_types,   only : RAD_VIS, RAD_NIR, RAD_LW
   use meds_test_support, only : banner, build_test_config, check, check_close
   implicit none

   type(meds_config_t)  :: cfg
   type(fast_context_t) :: ctx
   real(wp) :: om1, om2

   call banner('[pft] canopy optics: per PFT, and reaching the two-stream')
   cfg = build_test_config()

   !=== 1. The DEFAULTS are the literals the fast driver used to carry. =====================!
   call check_close(cfg%pft%leaf_reflect_vis(1),  0.10_wp,  1.0e-14_wp, 'default leaf VIS reflectance changed')
   call check_close(cfg%pft%leaf_transmit_nir(1), 0.25_wp,  1.0e-14_wp, 'default leaf NIR transmittance changed')
   call check_close(cfg%pft%leaf_emissivity(1),   0.97_wp,  1.0e-14_wp, 'default leaf emissivity changed')
   call check_close(cfg%pft%wood_emissivity(1),   0.90_wp,  1.0e-14_wp, 'default wood emissivity changed')
   call check_close(cfg%pft%leaf_clumping(1),     0.80_wp,  1.0e-14_wp, 'default leaf clumping changed')
   call check_close(cfg%pft%leaf_angle_mean(1),   45.0_wp,  1.0e-14_wp, 'default leaf-angle mean changed')

   !=== 2. Uniform traits => a uniform table. The control for test 3. =======================!
   call build_fast_context(cfg, ctx)
   call check(ctx%rad_opt%n_pft >= 2_ik, 'fixture must have at least two PFTs to differentiate')
   call check_close(ctx%rad_opt%omega_leaf(RAD_VIS,1), ctx%rad_opt%omega_leaf(RAD_VIS,2), 1.0e-14_wp,       &
                    'identical traits must give identical optics')

   !=== 3. Differentiate ONE PFT and the table must follow it -- only that one. ==============!
   !       Single-scatter albedo is reflectance + transmittance, so raising the reflectance    !
   !       of PFT 2 must raise ITS omega and leave PFT 1 alone.                                !
   om1 = ctx%rad_opt%omega_leaf(RAD_VIS,1)
   cfg%pft%leaf_reflect_vis(2) = 0.30_wp
   call build_fast_context(cfg, ctx)
   om2 = ctx%rad_opt%omega_leaf(RAD_VIS,2)
   call check_close(ctx%rad_opt%omega_leaf(RAD_VIS,1), om1, 1.0e-14_wp,                                    &
                    'changing PFT 2 must not move PFT 1 (per-PFT indexing)')
   call check_close(om2, 0.30_wp + cfg%pft%leaf_transmit_vis(2), 1.0e-12_wp,                               &
                    'leaf VIS omega must be reflectance + transmittance of THAT PFT')
   call check(om2 > om1, 'a more reflective PFT must have the larger single-scatter albedo')

   !=== 4. Longwave is configured as EMISSIVITY: omega = 1 - emissivity, transmittance 0. ====!
   cfg%pft%leaf_emissivity(2) = 0.80_wp
   call build_fast_context(cfg, ctx)
   call check_close(ctx%rad_opt%omega_leaf(RAD_LW,2), 1.0_wp - 0.80_wp, 1.0e-12_wp,                        &
                    'leaf LW omega must be 1 - emissivity (a leaf is opaque in the thermal band)')

   !=== 5. Clumping and the leaf-angle distribution are per PFT too. ========================!
   cfg%pft%leaf_clumping(2)   = 0.40_wp
   cfg%pft%leaf_angle_mean(2) = 20.0_wp        ! erectophile against PFT 1's spherical 45 deg
   call build_fast_context(cfg, ctx)
   call check_close(ctx%rad_opt%clumping_leaf(2), 0.40_wp, 1.0e-14_wp, 'leaf clumping did not reach the table')
   call check_close(ctx%rad_opt%clumping_leaf(1), 0.80_wp, 1.0e-14_wp, 'PFT 1 clumping must be untouched')
   call check(abs(ctx%rad_opt%bf(2) - ctx%rad_opt%bf(1)) > 1.0e-6_wp,                            &
              'a different leaf-angle mean must give a different angle distribution')

   !=== 5b. PER-PFT HYDRAULIC TRAITS (#179). Same shape as the optics checks above, and for the =!
   !        same reason: until now every PFT shared ONE hydraulic parameter set, so wood density
   !        was the only axis on which PFTs could differ hydraulically -- in a model whose point is
   !        that plant strategies differ. Assert that the table SEPARATES them and that the
   !        separation CHANGES the answer, not merely that it compiles.
   block
      type(hydro_params_table_t) :: tab
      real(wp) :: plc1, plc2
      cfg%pft%hyd_wood_psi50(2)  = -0.8_wp        ! vs the shared -2.0
      cfg%pft%hyd_leaf_pi0(2)    = -2.5_wp        ! vs the shared -1.5
      cfg%pft%hyd_k_plant_max(2) =  3.0e-4_wp     ! vs the shared 6.0e-4
      call apply_hydraulics_config(cfg%hydraulics, cfg%pft, tab)
      call check_close(tab%pft(1)%wood_psi50, cfg%hydraulics%wood_psi50, 1.0e-12_wp,              &
                       'PFT 1 with no override takes the shared [hydraulics] value')
      call check_close(tab%pft(2)%wood_psi50, -0.8_wp,   1.0e-12_wp, 'PFT 2 takes its own psi50')
      call check_close(tab%pft(2)%leaf_pi0,   -2.5_wp,   1.0e-12_wp, 'PFT 2 takes its own leaf pi0')
      call check_close(tab%pft(2)%k_plant_max, 3.0e-4_wp, 1.0e-15_wp, 'PFT 2 takes its own conductance')
      !----- One override must not require restating the other twelve. -------------------------!
      call check_close(tab%pft(2)%wood_kmax, cfg%hydraulics%wood_kmax, 1.0e-12_wp,                &
                       'an unset trait still falls back to the shared value, per PFT')
      !----- And it MATTERS: at the same potential the more vulnerable xylem has lost more. -----!
      plc1 = 1.0_wp - plc_retained(-1.5_wp, tab%pft(1)%wood_psi50, tab%pft(1)%wood_kexp)
      plc2 = 1.0_wp - plc_retained(-1.5_wp, tab%pft(2)%wood_psi50, tab%pft(2)%wood_kexp)
      call check(plc2 > plc1, 'a more vulnerable PFT loses more conductance at the same psi')
   end block

   !=== 6. The PFT-parameter CSV dump has as many values as it has column headers. ==========!
   !        A Fortran format SHORTER than its output list does not fail: it REVERTS to the last  !
   !        repeat group and keeps going, so an integer edit descriptor silently receives a real  !
   !        and prints its bit pattern, and the trailing columns disappear. That is exactly what  !
   !        adding two PFT traits did (#118), and nothing in the suite noticed, because the run   !
   !        completed and every budget still closed. Counting the two lines is the whole test.    !
   block
      character(len=*), parameter :: CSVFILE = 'test_pft_params_tmp.csv'
      character(len=8192) :: head_line, val_line
      integer(ik)         :: u, ios, nhead, nval
      call build_fast_context(cfg, ctx)          ! keep cfg consistent with the edits above
      call write_pft_params_csv(cfg, CSVFILE)
      open(newunit=u, file=CSVFILE, status='old', action='read', iostat=ios)
      call check(ios == 0, 'PFT parameter CSV must be written')
      read(u,'(a)') head_line
      read(u,'(a)') val_line
      close(u, status='delete')
      nhead = count_fields(head_line)
      nval  = count_fields(val_line)
      call check(nhead == nval, 'PFT CSV: header column count must equal the value count')
      if (nhead /= nval) print '(a,i0,a,i0)', '   header fields = ', nhead, ', value fields = ', nval
      !----- And no field may be a bit pattern: every value is either an integer or es15.8, so a  !
      !      19-digit run of digits is the format-reversion signature, not a real number. --------!
      call check(index(val_line, '4605380978949069210') == 0,                                     &
                 'PFT CSV: a real printed through an integer descriptor (format reversion)')
   end block

   print '(a)', 'test_pft_optics_config: all checks passed'

contains

   pure integer(ik) function count_fields(line) result(n)
      character(len=*), intent(in) :: line
      integer(ik) :: i
      n = 1_ik
      do i = 1_ik, len_trim(line)
         if (line(i:i) == ',') n = n + 1_ik
      end do
   end function count_fields

end program test_pft_optics_config
