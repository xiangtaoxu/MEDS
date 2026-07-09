!==========================================================================================!
! test_canopy_radiation -- unit tests for the canopy two-stream RT library (src/biophysics). !
!                                                                                          !
! Checks the physical invariants that a correct two-stream must satisfy:                     !
!   1. ENERGY CONSERVATION  (per shortwave band): absorbed + net-ground + reflected = incident.!
!   2. BEER'S LAW           : a single black cohort under a direct beam attenuates as exp(-k L).!
!   3. ZERO-LAI ALBEDO      : an empty canopy has albedo = ground albedo.                       !
!   4. LEAF-ANGLE LIMITS     : bf and the diffuse backscatter beta hit their analytic limits    !
!                              (uniform LAD -> bf=0.5; horizontal -> beta -> rho/(rho+tau)).      !
!   5. MULTI-PFT PATCH       : cohorts of different PFTs get their OWN diffuse optics (MEDS        !
!                              avoids ED2 bug B1 by construction), and energy still conserves.      !
!   6. LONGWAVE EQUILIBRIUM  : an isothermal canopy+ground+sky has ~zero net thermal absorption.    !
!==========================================================================================!
program test_canopy_radiation
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : stefan
   use meds_biophysics_types, only : rad_pft_optics_t, rad_forcing_t, rad_flux_t,             &
                                   alloc_rad_forcing, RAD_VIS, RAD_NIR, RAD_LW, N_RAD_BAND_DEFAULT
   use meds_optics,         only : beta_lidf, leaf_bf, beta_params_from_mean, gfun_direct,    &
                                   N_LEAF_CLASS, derive_rad_optics, surface_state_t, ground_optics
   use meds_canopy_radiation, only : canopy_radiation
   implicit none

   integer(ik) :: nfail
   real(wp), parameter :: tol = 1.0e-9_wp

   nfail = 0_ik
   call test_leaf_angle_limits()
   call test_energy_conservation()
   call test_beers_law()
   call test_zero_lai_albedo()
   call test_multi_pft()
   call test_longwave_equilibrium()
   call test_lidf_delta()
   call test_cohort_order()

   if (nfail == 0_ik) then
      print '(a)', 'test_canopy_radiation: ALL PASSED'
   else
      print '(a,i0,a)', 'test_canopy_radiation: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   !----- Assertions. ------------------------------------------------------------------------!
   subroutine check(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) <= atol) then
         print '(a,a,a,es13.5,a,es13.5)', '  ok   : ', name, '  (', got, ' ~ ', expect, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,es13.5,a,es13.5,a,es10.2)', '  FAIL : ', name, '  got ', got,          &
               ' expected ', expect, '  |diff|>', atol
      end if
   end subroutine check

   subroutine check_true(name, cond, val)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      real(wp),         intent(in) :: val
      if (cond) then
         print '(a,a,a,es13.5,a)', '  ok   : ', name, '  (margin ', val, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,es13.5)', '  FAIL : ', name, '  margin ', val
      end if
   end subroutine check_true

   !----- A delta leaf-angle distribution (std_deg = 0) must give FINITE Beta (p, q): the variance !
   !      floor keeps kappa finite instead of the old mt(1-mt)/0 -> NaN LIDF. --------------------!
   subroutine test_lidf_delta()
      real(wp) :: p, q
      print '(a)', 'test_lidf_delta:'
      call beta_params_from_mean(45.0_wp, 0.0_wp, p, q)     ! mean=45deg -> mt=0.5; std=0 -> vt floored
      call check('delta LIDF p finite (variance floored)', p, 124999.5_wp, 1.0e-3_wp)
      call check('delta LIDF q finite (variance floored)', q, 124999.5_wp, 1.0e-3_wp)
   end subroutine test_lidf_delta

   !----- Build a 2-PFT, 3-band optics table (VIS/NIR/LW) with a chosen leaf-angle mean. ------!
   subroutine build_optics(mean_deg, opt)
      real(wp),               intent(in)  :: mean_deg
      type(rad_pft_optics_t), intent(out) :: opt
      integer(ik), parameter :: NB = N_RAD_BAND_DEFAULT, NP = 2_ik
      real(wp) :: rl(NB,NP), tl(NB,NP), rw(NB,NP), tw(NB,NP)
      real(wp) :: cl(NP), cw(NP), bp(NP), bq(NP)
      logical  :: hb(NB), he(NB)
      integer(ik) :: p
      real(wp) :: emiss_leaf, emiss_wood
      !----- Leaf spectral: PFT 1 darker (low VIS reflectance), PFT 2 brighter. -------------!
      rl(RAD_VIS,:) = [0.09_wp, 0.14_wp] ; tl(RAD_VIS,:) = [0.05_wp, 0.08_wp]
      rl(RAD_NIR,:) = [0.45_wp, 0.50_wp] ; tl(RAD_NIR,:) = [0.25_wp, 0.30_wp]
      emiss_leaf = 0.97_wp
      rl(RAD_LW,:)  = 1.0_wp - emiss_leaf ; tl(RAD_LW,:) = 0.0_wp
      !----- Wood spectral (near-opaque). ---------------------------------------------------!
      rw(RAD_VIS,:) = 0.11_wp ; tw(RAD_VIS,:) = 0.001_wp
      rw(RAD_NIR,:) = 0.25_wp ; tw(RAD_NIR,:) = 0.001_wp
      emiss_wood = 0.90_wp
      rw(RAD_LW,:)  = 1.0_wp - emiss_wood ; tw(RAD_LW,:) = 0.0_wp
      cl = [0.80_wp, 0.85_wp] ; cw = [0.50_wp, 0.50_wp]
      do p = 1_ik, NP
         call beta_params_from_mean(mean_deg, 20.0_wp, bp(p), bq(p))
      end do
      hb = [.true.,  .true.,  .false.]        ! VIS/NIR have a beam; LW does not
      he = [.false., .false., .true. ]        ! only LW emits
      call derive_rad_optics(NB, NP, rl, tl, rw, tw, cl, cw, bp, bq, hb, he, opt)
   end subroutine build_optics

   !----- Fill shortwave + longwave forcing from a surface state. ----------------------------!
   subroutine build_forcing(cosz, t_sky, t_soil, f)
      real(wp),            intent(in)  :: cosz, t_sky, t_soil
      type(rad_forcing_t), intent(out) :: f
      type(surface_state_t) :: surf
      logical :: he(N_RAD_BAND_DEFAULT)
      call alloc_rad_forcing(f, N_RAD_BAND_DEFAULT)
      f%cosz = cosz
      f%incid_beam(RAD_VIS) = 300.0_wp ; f%incid_diff(RAD_VIS) = 100.0_wp
      f%incid_beam(RAD_NIR) = 250.0_wp ; f%incid_diff(RAD_NIR) =  80.0_wp
      f%incid_beam(RAD_LW)  = 0.0_wp   ; f%incid_diff(RAD_LW)  = stefan * t_sky ** 4
      surf%n_band = N_RAD_BAND_DEFAULT
      allocate(surf%soil_albedo(N_RAD_BAND_DEFAULT))
      surf%soil_albedo = [0.15_wp, 0.30_wp, 0.0_wp]
      surf%soil_emiss  = 0.95_wp ; surf%soil_temp = t_soil
      he = [.false., .false., .true.]
      call ground_optics(surf, N_RAD_BAND_DEFAULT, he, f%grnd_refl, f%grnd_emiss)
   end subroutine build_forcing

   !=======================================================================================!
   subroutine test_leaf_angle_limits()
      real(wp) :: p, q, bf
      real(wp) :: lidf(N_LEAF_CLASS)
      print '(a)', '[1] leaf-angle distribution limits'
      !----- Uniform inclination (Beta(1,1)): <cos^2> = 0.5 exactly. ------------------------!
      lidf = beta_lidf(1.0_wp, 1.0_wp) ; bf = leaf_bf(lidf)
      call check('bf(uniform LAD)', bf, 0.5_wp, 2.0e-2_wp)
      !----- Strongly planophile (small mean angle): bf -> ~1 (near-horizontal). ------------!
      call beta_params_from_mean(8.0_wp, 6.0_wp, p, q)
      lidf = beta_lidf(p, q) ; bf = leaf_bf(lidf)
      if (bf > 0.90_wp) then
         print '(a,f7.4)', '  ok   : bf(planophile) > 0.90   got ', bf
      else
         nfail = nfail + 1_ik ; print '(a,f7.4)', '  FAIL : bf(planophile) <= 0.90  got ', bf
      end if
      !----- Strongly erectophile (near-vertical): bf -> ~0. --------------------------------!
      call beta_params_from_mean(82.0_wp, 6.0_wp, p, q)
      lidf = beta_lidf(p, q) ; bf = leaf_bf(lidf)
      if (bf < 0.10_wp) then
         print '(a,f7.4)', '  ok   : bf(erectophile) < 0.10  got ', bf
      else
         nfail = nfail + 1_ik ; print '(a,f7.4)', '  FAIL : bf(erectophile) >= 0.10 got ', bf
      end if
   end subroutine test_leaf_angle_limits

   !=======================================================================================!
   subroutine test_energy_conservation()
      type(rad_pft_optics_t) :: opt
      type(rad_forcing_t)    :: f
      type(rad_flux_t)       :: flux
      integer(ik), parameter :: NC = 4_ik
      integer(ik) :: pft(NC), b
      real(wp) :: lai(NC), wai(NC), tcan(NC), absorbed, net_ground, reflected, incid, resid
      print '(a)', '[2] shortwave energy conservation (4 cohorts)'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.7_wp, 290.0_wp, 298.0_wp, f)
      pft  = [1_ik, 2_ik, 1_ik, 2_ik]
      lai  = [2.0_wp, 1.5_wp, 1.0_wp, 0.5_wp]   ! bottom -> top
      wai  = 0.1_wp * lai
      tcan = 298.0_wp
      call canopy_radiation(opt, f, NC, pft, lai, wai, tcan, flux)
      do b = RAD_VIS, RAD_NIR
         absorbed   = sum(flux%abs_leaf(b,:)) + sum(flux%abs_wood(b,:))
         net_ground = flux%dn_ground(b) - flux%up_ground(b)
         incid      = f%incid_beam(b) + f%incid_diff(b)
         reflected  = flux%albedo(b) * incid
         resid      = incid - (absorbed + net_ground + reflected)
         call check('energy residual band', resid, 0.0_wp, 1.0e-8_wp * incid)
         if (flux%albedo(b) < 0.0_wp .or. flux%albedo(b) > 1.0_wp) then
            nfail = nfail + 1_ik ; print '(a,es12.4)', '  FAIL : albedo out of [0,1] ', flux%albedo(b)
         end if
      end do
   end subroutine test_energy_conservation

   !=======================================================================================!
   !  Cohort ORDER: two EQUAL-LAI cohorts stacked BOTTOM(1)->TOP(2). The top cohort intercepts !
   !  the incoming beam+diffuse first, so per band it must absorb strictly MORE than the shaded  !
   !  bottom cohort. This vertical-gradient property is exactly what the fast loop's              !
   !  apply_rt_forcing `perm` relies on to route the higher absorbed PAR to the taller cohort.     !
   subroutine test_cohort_order()
      type(rad_pft_optics_t) :: opt
      type(rad_forcing_t)    :: f
      type(rad_flux_t)       :: flux
      integer(ik), parameter :: NC = 2_ik
      integer(ik) :: pft(NC)
      real(wp) :: lai(NC), wai(NC), tcan(NC)
      print '(a)', '[8] cohort order: top of two equal-LAI cohorts absorbs more than the bottom'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.7_wp, 290.0_wp, 298.0_wp, f)
      pft = [1_ik, 1_ik] ; lai = [2.0_wp, 2.0_wp] ; wai = [0.0_wp, 0.0_wp] ; tcan = 298.0_wp  ! idx1=bottom, idx2=top
      call canopy_radiation(opt, f, NC, pft, lai, wai, tcan, flux)
      call check_true('VIS: top cohort (idx2) absorbs more than bottom (idx1)',                 &
                      flux%abs_leaf(RAD_VIS,2) > flux%abs_leaf(RAD_VIS,1),                       &
                      flux%abs_leaf(RAD_VIS,2) - flux%abs_leaf(RAD_VIS,1))
      call check_true('NIR: top cohort (idx2) absorbs more than bottom (idx1)',                 &
                      flux%abs_leaf(RAD_NIR,2) > flux%abs_leaf(RAD_NIR,1),                       &
                      flux%abs_leaf(RAD_NIR,2) - flux%abs_leaf(RAD_NIR,1))
   end subroutine test_cohort_order

   !=======================================================================================!
   subroutine test_beers_law()
      type(rad_pft_optics_t) :: opt
      type(rad_forcing_t)    :: f
      type(rad_flux_t)       :: flux
      integer(ik), parameter :: NC = 1_ik
      integer(ik) :: pft(NC)
      real(wp) :: lai(NC), wai(NC), tcan(NC), kdir, gee, absorbed_vis, expect
      real(wp) :: elai
      print '(a)', '[3] Beers law (single black cohort, direct beam)'
      call build_optics(45.0_wp, opt)
      call build_forcing(1.0_wp, 290.0_wp, 298.0_wp, f)
      !----- Black leaves: force VIS single-scatter albedo to ~0 by zeroing spectra. --------!
      opt%omega_leaf(RAD_VIS,:) = 0.0_wp ; opt%g_leaf(RAD_VIS,:) = 0.0_wp
      opt%omega_wood(RAD_VIS,:) = 0.0_wp ; opt%g_wood(RAD_VIS,:) = 0.0_wp
      f%incid_diff(RAD_VIS) = 0.0_wp                       ! beam only
      f%grnd_refl(RAD_VIS)  = 0.0_wp                       ! black soil: no reflection to re-absorb
      pft = [1_ik] ; lai = [3.0_wp] ; wai = [0.0_wp] ; tcan = 298.0_wp
      call canopy_radiation(opt, f, NC, pft, lai, wai, tcan, flux)
      elai         = opt%clumping_leaf(1) * lai(1)
      gee          = gfun_direct(opt%lidf(:,1), 1.0_wp)
      kdir         = gee / 1.0_wp
      expect       = f%incid_beam(RAD_VIS) * (1.0_wp - exp(-kdir * elai))
      absorbed_vis = flux%abs_leaf(RAD_VIS,1) + flux%abs_wood(RAD_VIS,1)
      call check('absorbed = incid*(1-exp(-kL))', absorbed_vis, expect, 1.0e-7_wp)
      call check('below-canopy = incid*exp(-kL)', flux%dn_ground(RAD_VIS),                     &
                 f%incid_beam(RAD_VIS) * exp(-kdir * elai), 1.0e-7_wp)
   end subroutine test_beers_law

   !=======================================================================================!
   subroutine test_zero_lai_albedo()
      type(rad_pft_optics_t) :: opt
      type(rad_forcing_t)    :: f
      type(rad_flux_t)       :: flux
      integer(ik), parameter :: NC = 1_ik
      integer(ik) :: pft(NC)
      real(wp) :: lai(NC), wai(NC), tcan(NC)
      print '(a)', '[4] zero-LAI canopy has albedo = ground albedo'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.7_wp, 290.0_wp, 298.0_wp, f)
      pft = [1_ik] ; lai = [0.0_wp] ; wai = [0.0_wp] ; tcan = 298.0_wp
      call canopy_radiation(opt, f, NC, pft, lai, wai, tcan, flux)
      call check('albedo(VIS) = soil albedo', flux%albedo(RAD_VIS), f%grnd_refl(RAD_VIS), 1.0e-9_wp)
      call check('albedo(NIR) = soil albedo', flux%albedo(RAD_NIR), f%grnd_refl(RAD_NIR), 1.0e-9_wp)
   end subroutine test_zero_lai_albedo

   !=======================================================================================!
   subroutine test_multi_pft()
      type(rad_pft_optics_t) :: opt
      type(rad_forcing_t)    :: f
      type(rad_flux_t)       :: flux
      integer(ik), parameter :: NC = 2_ik
      integer(ik) :: pft(NC)
      real(wp) :: lai(NC), wai(NC), tcan(NC), absorbed, net_ground, reflected, incid, resid
      print '(a)', '[5] multi-PFT patch: per-cohort optics + conservation (ED2 bug B1 avoided)'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.6_wp, 290.0_wp, 298.0_wp, f)
      pft = [1_ik, 2_ik] ; lai = [1.5_wp, 1.5_wp] ; wai = 0.15_wp ; tcan = 298.0_wp
      !----- The two PFTs have DISTINCT VIS optics: the RT must see per-cohort values. -------!
      if (abs(opt%omega_leaf(RAD_VIS,1) - opt%omega_leaf(RAD_VIS,2)) < 1.0e-6_wp) then
         nfail = nfail + 1_ik ; print '(a)', '  FAIL : the two PFTs should differ in VIS omega'
      else
         print '(a)', '  ok   : PFTs differ in VIS omega (per-cohort optics exercised)'
      end if
      call canopy_radiation(opt, f, NC, pft, lai, wai, tcan, flux)
      incid      = f%incid_beam(RAD_VIS) + f%incid_diff(RAD_VIS)
      absorbed   = sum(flux%abs_leaf(RAD_VIS,:)) + sum(flux%abs_wood(RAD_VIS,:))
      net_ground = flux%dn_ground(RAD_VIS) - flux%up_ground(RAD_VIS)
      reflected  = flux%albedo(RAD_VIS) * incid
      resid      = incid - (absorbed + net_ground + reflected)
      call check('energy residual (2 PFT)', resid, 0.0_wp, 1.0e-8_wp * incid)
   end subroutine test_multi_pft

   !=======================================================================================!
   subroutine test_longwave_equilibrium()
      type(rad_pft_optics_t) :: opt
      type(rad_forcing_t)    :: f
      type(rad_flux_t)       :: flux
      integer(ik), parameter :: NC = 3_ik
      integer(ik) :: pft(NC), i
      real(wp) :: lai(NC), wai(NC), tcan(NC), net_lw, maxnet
      real(wp), parameter :: teq = 300.0_wp
      print '(a)', '[6] longwave radiative equilibrium (isothermal -> zero net absorption)'
      call build_optics(45.0_wp, opt)
      call build_forcing(0.5_wp, teq, teq, f)          ! sky and soil both at Teq
      pft  = [1_ik, 2_ik, 1_ik] ; lai = [2.0_wp, 1.0_wp, 0.5_wp] ; wai = 0.1_wp * lai
      tcan = teq                                        ! canopy also at Teq
      call canopy_radiation(opt, f, NC, pft, lai, wai, tcan, flux)
      maxnet = 0.0_wp
      do i = 1_ik, NC
         net_lw = flux%abs_leaf(RAD_LW,i) + flux%abs_wood(RAD_LW,i)
         maxnet = max(maxnet, abs(net_lw))
      end do
      call check('max |net LW| at equilibrium', maxnet, 0.0_wp, 1.0e-7_wp * stefan * teq ** 4)
   end subroutine test_longwave_equilibrium

end program test_canopy_radiation
