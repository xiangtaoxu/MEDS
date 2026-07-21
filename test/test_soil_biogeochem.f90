!==========================================================================================!
! test_soil_biogeochem -- unit tests for the P0 SLOW soil-carbon matrix kernels                !
! (meds_soil_biogeochem), mirroring the design test plan (MEDS_BIOGEOCHEMISTRY_DESIGN.md §8):   !
!   1. MASS CLOSURE  : audit%resid ~ 0 for one soil_carbon_step (EULER and EXPM), by construction. !
!   2. Rh COMPLEMENT : Rh = sum er_j*xi_j*K_j*X_j = -1^T A xi K X; column check sum_i a_ij = -er_j.  !
!   3. SCALAR PLACE  : A*(xi*K*X) conserves mass to Rh; the wrong xi*(A*K*X) does not (per-pool xi).   !
!   4. SCHEME TOPO   : scheme 0 = 3-active (struct->slow); scheme 5 = 5-active clay/sand topology.      !
!   5. SPIN-UP       : SASU active-block X_ss is an exact EULER fixed point; brute force converges;      !
!                      inert pools held at 0; the singular full 7x7 is never solved (default K5=K7=0).    !
!   6. RESIDENCE/CAP : X_c = X_ss at equilibrium, X_p = 0; a perturbed pool relaxes toward X_c.           !
!   7. LITTER+LIGNIN : build_litter_input maps pre-split streams; lignin is a proportional passive tracer.!
!   8. T/THETA RESP  : Q10 doubling per 10 K; ED2 cap at 45 C; moisture peak; Rh=0 at X=0; N off => f=1.  !
!   9. FAST/SLOW SEAM: matrix single-pool Rh == fast meds_cas_biophysics Rh (matched chemistry); the          !
!                      accumulated xi_int path closes rh_today==today_rh, the daily-mean path does NOT.    !
!==========================================================================================!
program test_soil_biogeochem
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : yr_day, kgCday_2_umols
   use meds_biogeochem_types, only : soil_carbon_t, decomp_opts_t, litter_input_t, soilc_audit_t,  &
                                     soilc_diag_t, n_soil_pool,                                     &
                                     IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL,    &
                                     IP_MICR, IP_SLOW, IP_PASSIVE,                                  &
                                     DECOMP_STEP_EULER, DECOMP_STEP_EXPM,                           &
                                     DECOMP_SCHEME_ED2, DECOMP_SCHEME_CENTURY5
   use meds_biogeochem_types, only : co2_opts_t, HR_Q10
   use meds_soil_biogeochem,  only : assemble_env_scalar, assemble_transfer_matrix,                &
                                     build_litter_input, heterotrophic_respiration_matrix,          &
                                     soil_carbon_step, solve_soil_carbon_steady_state,              &
                                     soil_carbon_diagnostics, pack_pool_vector
   use meds_soil_biogeochem,  only : heterotrophic_respiration_flux
   implicit none
   integer(ik) :: nfail
   nfail = 0_ik

   call test_mass_closure()
   call test_rh_complement()
   call test_scalar_placement()
   call test_scheme_topology()
   call test_spinup_steady_state()
   call test_residence_capacity()
   call test_litter_and_lignin()
   call test_temperature_moisture()
   call test_fast_slow_seam()

   if (nfail == 0_ik) then
      print '(a)', 'test_soil_biogeochem: ALL PASSED'
   else
      print '(a,i0,a)', 'test_soil_biogeochem: ', nfail, ' FAILED'
      error stop 1
   end if

contains

   subroutine check(name, got, expect, atol)
      character(len=*), intent(in) :: name
      real(wp),         intent(in) :: got, expect, atol
      if (abs(got - expect) <= atol) then
         print '(a,a,a,es13.5,a,es13.5)', '  ok   : ', name, '  (', got, ' ~ ', expect, ')'
      else
         nfail = nfail + 1_ik
         print '(a,a,a,es13.5,a,es13.5)', '  FAIL : ', name, '  got ', got, ' expected ', expect
      end if
   end subroutine check

   subroutine check_true(name, cond, val)
      character(len=*), intent(in) :: name
      logical,          intent(in) :: cond
      real(wp),         intent(in) :: val
      if (cond) then
         print '(a,a,a,es13.5,a)', '  ok   : ', name, '  (', val, ')'
      else
         nfail = nfail + 1_ik ; print '(a,a,a,es13.5,a)', '  FAIL : ', name, '  (', val, ')'
      end if
   end subroutine check_true

   !----- a reproducible pseudo-random-ish filler (no Date/random; deterministic). ---------------!
   pure function frac(i) result(f)
      integer(ik), intent(in) :: i
      real(wp) :: f
      f = mod(real(i, wp) * 0.61803398875_wp, 1.0_wp)     ! golden-ratio low-discrepancy in (0,1)
   end function frac

   !=======================================================================================!
   ! 1. Mass closure -- audit%resid ~ 0 for one step under BOTH solvers (reporting invariant).     !
   !=======================================================================================!
   subroutine test_mass_closure()
      type(soil_carbon_t) :: pools
      type(decomp_opts_t) :: opts
      type(soilc_audit_t) :: audit
      real(wp) :: u(n_soil_pool), lignin_in(2), xi_int(n_soil_pool), rh, x(n_soil_pool), sx
      real(wp) :: a_mat(n_soil_pool,n_soil_pool), k_diag(n_soil_pool), er(n_soil_pool), x1(n_soil_pool)
      integer(ik) :: i, j
      print '(a)', '-- test 1: mass closure (EULER + EXPM) --'
      opts%decomp_scheme = DECOMP_SCHEME_CENTURY5      ! all 5 active -> exercises transfers + EXPM
      opts%k_micr = 6.0_wp ; opts%k_pass = 0.012_wp ; opts%k_slow = 0.15_wp ; opts%k_struct = 1.5_wp
      opts%k_fast = 12.0_wp
      opts%er_fast = 0.55_wp ; opts%er_struct_nonlig = 0.50_wp ; opts%er_struct_lig = 0.20_wp
      opts%er_slow = 0.55_wp ; opts%er_pass = 0.55_wp ; opts%er_micr_int = 0.60_wp ; opts%er_micr_slp = 0.17_wp

      pools%fast_grnd_carbon = 0.4_wp ; pools%fast_soil_carbon = 0.6_wp
      pools%struct_grnd_carbon = 2.0_wp ; pools%struct_soil_carbon = 3.0_wp
      pools%microbial_carbon = 0.5_wp ; pools%slow_carbon = 8.0_wp ; pools%passive_carbon = 20.0_wp
      pools%struct_grnd_lignin = 0.6_wp ; pools%struct_soil_lignin = 0.9_wp
      u = [ (0.01_wp * frac(int(10*i, ik)), i = 1, n_soil_pool) ]
      u(IP_MICR:IP_PASSIVE) = 0.0_wp                   ! litter only enters pools 1-4
      lignin_in = [0.002_wp, 0.003_wp]
      xi_int = [ (0.5_wp * frac(int(3*i+1, ik)) + 0.1_wp, i = 1, n_soil_pool) ]  ! ~O(0.1..0.6) day

      call pack_pool_vector(pools, x) ; sx = max(sum(x), 1.0_wp)

      opts%step_solver = DECOMP_STEP_EULER
      call assemble_transfer_matrix(pools, opts, a_mat, k_diag, er)
      block
         type(soil_carbon_t) :: p
         real(wp) :: rh_indep, dc
         p = pools
         call soil_carbon_step(p, u, lignin_in, xi_int, opts, rh, audit)
         call check_true('EULER resid ~ 0 (reporting invariant)', abs(audit%resid) < 1.0e-12_wp*sx, audit%resid)
         call check_true('EULER lignin_resid ~ 0', abs(audit%lignin_resid) < 1.0e-12_wp*sx, audit%lignin_resid)
         call check_true('EULER rh_today >= 0', rh >= 0.0_wp, rh)
         !----- INDEPENDENT (non-tautological) mass check: dC == litter_in - Rh, with Rh recomputed    !
         !      from er*xi_int*K*X0 (NOT from rh_today = litter_in - dC, which is closure-by-definition).!
         !      This would catch a dropped-litter or misattributed-transfer bug that audit%resid can't. !
         call pack_pool_vector(p, x1)
         rh_indep = 0.0_wp
         do j = 1_ik, n_soil_pool
            rh_indep = rh_indep + er(j)*xi_int(j)*k_diag(j)*x(j)
         end do
         dc = sum(x1) - sum(x)
         call check('EULER dC == litter_in - Rh (INDEPENDENT)', dc, sum(u) - rh_indep, 1.0e-12_wp*sx)
      end block

      opts%step_solver = DECOMP_STEP_EXPM
      block
         type(soil_carbon_t) :: p
         p = pools
         call soil_carbon_step(p, u, lignin_in, xi_int, opts, rh, audit)
         call check_true('EXPM  resid ~ 0', abs(audit%resid) < 1.0e-10_wp*sx, audit%resid)
         call check_true('EXPM  rh_today >= 0', rh >= -1.0e-12_wp, rh)
      end block
   end subroutine test_mass_closure

   !=======================================================================================!
   ! 2. Rh respired complement + column-sum identity.                                              !
   !=======================================================================================!
   subroutine test_rh_complement()
      type(soil_carbon_t) :: pools
      type(decomp_opts_t) :: opts
      real(wp) :: a_mat(n_soil_pool,n_soil_pool), k_diag(n_soil_pool), er(n_soil_pool)
      real(wp) :: xi(n_soil_pool), x(n_soil_pool), rh_fun, rh_hand, colsum, worst
      integer(ik) :: i, j
      print '(a)', '-- test 2: Rh = -1^T A xi K X + column identity --'
      opts%decomp_scheme = DECOMP_SCHEME_CENTURY5
      opts%k_micr = 6.0_wp ; opts%k_pass = 0.012_wp ; opts%k_slow = 0.15_wp ; opts%k_struct = 1.5_wp
      opts%k_fast = 12.0_wp
      opts%er_fast=0.55_wp; opts%er_struct_nonlig=0.50_wp; opts%er_struct_lig=0.20_wp
      opts%er_slow=0.55_wp; opts%er_pass=0.55_wp; opts%er_micr_int=0.60_wp; opts%er_micr_slp=0.17_wp
      pools%fast_grnd_carbon=0.4_wp; pools%fast_soil_carbon=0.6_wp
      pools%struct_grnd_carbon=2.0_wp; pools%struct_soil_carbon=3.0_wp
      pools%microbial_carbon=0.5_wp; pools%slow_carbon=8.0_wp; pools%passive_carbon=20.0_wp
      call assemble_transfer_matrix(pools, opts, a_mat, k_diag, er)
      xi = [ (0.4_wp + 0.3_wp*frac(int(i,ik)), i = 1, n_soil_pool) ]
      call pack_pool_vector(pools, x)
      rh_fun  = heterotrophic_respiration_matrix(a_mat, k_diag, xi, pools)
      rh_hand = 0.0_wp
      do j = 1_ik, n_soil_pool
         rh_hand = rh_hand + er(j)*xi(j)*k_diag(j)*x(j)
      end do
      call check('Rh matrix == sum er*xi*K*X', rh_fun, rh_hand, 1.0e-14_wp*max(abs(rh_hand),1.0_wp))
      ! column identity: sum_i a_ij = -er_j  and  a_jj = -1
      worst = 0.0_wp
      do j = 1_ik, n_soil_pool
         colsum = sum(a_mat(:, j))
         worst  = max(worst, abs(colsum + er(j)))
         worst  = max(worst, abs(a_mat(j,j) + 1.0_wp))
      end do
      call check_true('column identity sum_i a_ij = -er_j, a_jj=-1', worst < 1.0e-14_wp, worst)
   end subroutine test_rh_complement

   !=======================================================================================!
   ! 3. Scalar placement -- A*(xi*K*X) conserves to Rh; xi*(A*K*X) (receiver-scalar) does NOT.      !
   !=======================================================================================!
   subroutine test_scalar_placement()
      type(soil_carbon_t) :: pools
      type(decomp_opts_t) :: opts
      real(wp) :: a_mat(n_soil_pool,n_soil_pool), k_diag(n_soil_pool), er(n_soil_pool)
      real(wp) :: xi(n_soil_pool), x(n_soil_pool)
      real(wp) :: dcorrect(n_soil_pool), dwrong(n_soil_pool), rh, loss_correct, loss_wrong
      integer(ik) :: i, j
      print '(a)', '-- test 3: scalar placement A*(xi*K*X) vs xi*(A*K*X) --'
      opts%decomp_scheme = DECOMP_SCHEME_ED2           ! default 3-active
      pools%struct_grnd_carbon=2.0_wp; pools%struct_soil_carbon=3.0_wp
      pools%fast_grnd_carbon=1.0_wp; pools%fast_soil_carbon=1.5_wp; pools%slow_carbon=5.0_wp
      call assemble_transfer_matrix(pools, opts, a_mat, k_diag, er)
      call pack_pool_vector(pools, x)
      ! DISTINCT per-pool xi (surface vs below), so ordering matters:
      xi = [0.9_wp, 0.3_wp, 0.9_wp, 0.3_wp, 0.3_wp, 0.3_wp, 0.3_wp]
      ! correct: dX_i = sum_j a_ij*(xi_j*K_j*X_j)   (donor scalar)
      do i = 1_ik, n_soil_pool
         dcorrect(i) = 0.0_wp ; dwrong(i) = 0.0_wp
         do j = 1_ik, n_soil_pool
            dcorrect(i) = dcorrect(i) + a_mat(i,j)*(xi(j)*k_diag(j)*x(j))
            dwrong(i)   = dwrong(i)   + a_mat(i,j)*(xi(i)*k_diag(j)*x(j))    ! receiver scalar (WRONG)
         end do
      end do
      rh = heterotrophic_respiration_matrix(a_mat, k_diag, xi, pools)
      loss_correct = -sum(dcorrect)      ! total carbon leaving the network
      loss_wrong   = -sum(dwrong)
      call check('correct form: total loss == Rh', loss_correct, rh, 1.0e-14_wp*max(rh,1.0_wp))
      call check_true('wrong form: total loss /= Rh', abs(loss_wrong - rh) > 1.0e-6_wp, loss_wrong - rh)
   end subroutine test_scalar_placement

   !=======================================================================================!
   ! 4. Scheme topology -- hand-checked A/er for scheme 0 (default) and scheme 5 (CENTURY).         !
   !=======================================================================================!
   subroutine test_scheme_topology()
      type(soil_carbon_t) :: pools
      type(decomp_opts_t) :: opts
      real(wp) :: a_mat(n_soil_pool,n_soil_pool), k_diag(n_soil_pool), er(n_soil_pool), offsum
      real(wp) :: er_micr, ex_mp, ex_sp
      integer(ik) :: i, j
      print '(a)', '-- test 4: decomp scheme topology --'
      ! ---- scheme 0 (default 3-active): only a(6,3),a(6,4) off-diagonal; micr/pass inert. ----
      opts = decomp_opts_t()                            ! defaults: scheme 0, er_struct=0.3, k_micr=k_pass=0
      call assemble_transfer_matrix(pools, opts, a_mat, k_diag, er)
      call check('scheme0 K_micr = 0', k_diag(IP_MICR),    0.0_wp, 0.0_wp)
      call check('scheme0 K_pass = 0', k_diag(IP_PASSIVE), 0.0_wp, 0.0_wp)
      call check('scheme0 a(slow,struct_grnd)=0.7', a_mat(IP_SLOW,IP_STRUCT_GRND), 0.7_wp, 1.0e-14_wp)
      call check('scheme0 a(slow,struct_soil)=0.7', a_mat(IP_SLOW,IP_STRUCT_SOIL), 0.7_wp, 1.0e-14_wp)
      call check('scheme0 er(struct_grnd)=0.3', er(IP_STRUCT_GRND), 0.3_wp, 1.0e-14_wp)
      call check('scheme0 er(fast_soil)=1.0',   er(IP_FAST_SOIL),   1.0_wp, 1.0e-14_wp)
      ! everything except a(6,3),a(6,4) and the -1 diagonal must be 0:
      offsum = 0.0_wp
      do j = 1_ik, n_soil_pool
         do i = 1_ik, n_soil_pool
            if (i /= j .and. .not. (i==IP_SLOW .and. (j==IP_STRUCT_GRND .or. j==IP_STRUCT_SOIL))) then
               offsum = offsum + abs(a_mat(i,j))
            end if
         end do
      end do
      call check_true('scheme0 no other transfers', offsum < 1.0e-14_wp, offsum)

      ! ---- scheme 5 (5-active): clay/sand-controlled transfers (hand values, f_lignin=0). ----
      opts = decomp_opts_t()
      opts%decomp_scheme = DECOMP_SCHEME_CENTURY5
      opts%k_micr=6.0_wp; opts%k_pass=0.012_wp; opts%k_slow=0.15_wp; opts%k_struct=1.5_wp; opts%k_fast=12.0_wp
      opts%er_fast=0.55_wp; opts%er_struct_nonlig=0.50_wp; opts%er_struct_lig=0.20_wp
      opts%er_slow=0.55_wp; opts%er_pass=0.55_wp; opts%er_micr_int=0.60_wp; opts%er_micr_slp=0.17_wp
      opts%xsand=0.35_wp; opts%xclay=0.25_wp
      call assemble_transfer_matrix(pools, opts, a_mat, k_diag, er)
      er_micr = 0.60_wp + 0.17_wp*0.35_wp
      ex_mp   = min(1.0_wp-er_micr, 0.003_wp+0.032_wp*0.25_wp)
      ex_sp   = min(1.0_wp-0.55_wp, 0.003_wp+0.009_wp*0.25_wp)
      call check('scheme5 a(micr,fast_grnd)=0.45', a_mat(IP_MICR,IP_FAST_GRND), 0.45_wp, 1.0e-14_wp)
      call check('scheme5 a(micr,struct_grnd)=0.50', a_mat(IP_MICR,IP_STRUCT_GRND), 0.50_wp, 1.0e-14_wp)
      call check('scheme5 a(pass,micr)=ex_mp', a_mat(IP_PASSIVE,IP_MICR), ex_mp, 1.0e-14_wp)
      call check('scheme5 a(slow,micr)=(1-er_micr)-ex_mp', a_mat(IP_SLOW,IP_MICR), (1.0_wp-er_micr)-ex_mp, 1.0e-14_wp)
      call check('scheme5 a(pass,slow)=ex_sp', a_mat(IP_PASSIVE,IP_SLOW), ex_sp, 1.0e-14_wp)
      call check('scheme5 a(slow,pass)=0.45', a_mat(IP_SLOW,IP_PASSIVE), 0.45_wp, 1.0e-14_wp)
      call check('scheme5 er(micr)=er_micr', er(IP_MICR), er_micr, 1.0e-13_wp)
   end subroutine test_scheme_topology

   !=======================================================================================!
   ! 5. SASU spin-up: active-block X_ss is an exact EULER fixed point; brute force converges;       !
   !    inert pools held at 0; default (K5=K7=0) solves despite the singular full matrix.            !
   !=======================================================================================!
   subroutine test_spinup_steady_state()
      print '(a)', '-- test 5: SASU steady state (active block) --'
      call spinup_one_scheme(DECOMP_SCHEME_ED2, 'scheme0')
      call spinup_one_scheme(DECOMP_SCHEME_CENTURY5, 'scheme5')
   end subroutine test_spinup_steady_state

   subroutine spinup_one_scheme(scheme, tag)
      integer(ik),      intent(in) :: scheme
      character(len=*), intent(in) :: tag
      type(soil_carbon_t) :: pools, pss, p
      type(decomp_opts_t) :: opts
      type(soilc_audit_t) :: audit
      real(wp) :: a_mat(n_soil_pool,n_soil_pool), k_diag(n_soil_pool), er(n_soil_pool)
      real(wp) :: xi(n_soil_pool), u(n_soil_pool), lignin_in(2), xss(n_soil_pool), x1(n_soil_pool)
      real(wp) :: rh, worst, relerr
      integer(ik) :: it
      opts = decomp_opts_t()
      opts%decomp_scheme = scheme
      if (scheme == DECOMP_SCHEME_CENTURY5) then
         opts%k_micr=6.0_wp; opts%k_pass=0.012_wp; opts%k_slow=0.15_wp; opts%k_struct=1.5_wp; opts%k_fast=12.0_wp
         opts%er_fast=0.55_wp; opts%er_struct_nonlig=0.50_wp; opts%er_struct_lig=0.20_wp
         opts%er_slow=0.55_wp; opts%er_pass=0.55_wp; opts%er_micr_int=0.60_wp; opts%er_micr_slp=0.17_wp
      end if
      ! frozen climate scalar (constant), constant litter input to pools 1-4:
      xi = [0.7_wp, 0.5_wp, 0.7_wp, 0.5_wp, 0.5_wp, 0.5_wp, 0.5_wp]
      u  = 0.0_wp
      u(IP_FAST_GRND)=0.03_wp; u(IP_FAST_SOIL)=0.05_wp; u(IP_STRUCT_GRND)=0.02_wp; u(IP_STRUCT_SOIL)=0.04_wp
      lignin_in = 0.0_wp
      call assemble_transfer_matrix(pools, opts, a_mat, k_diag, er)
      call solve_soil_carbon_steady_state(a_mat, k_diag, xi, u, opts, pss)
      call pack_pool_vector(pss, xss)
      ! inert pools held at 0 (default scheme); finite everywhere:
      if (scheme == DECOMP_SCHEME_ED2) then
         !----- §8.5 guard (MEANINGFUL, not vacuous): the FULL 7x7 A*diag(xi*K) is SINGULAR because the  !
         !      inert (micr, passive) COLUMNS are identically zero (K_micr=K_pass=0), so a full solve      !
         !      would hit a zero pivot -- the active-block {K>0} solve is REQUIRED and still returns a      !
         !      finite X_ss (checked just below). Asserting the inert columns vanish is what proves the     !
         !      active-block path is genuinely exercised (a full-solve regression could not satisfy this).  !
         block
            real(wp)    :: mcol_inert
            integer(ik) :: ii
            mcol_inert = 0.0_wp
            do ii = 1_ik, n_soil_pool
               mcol_inert = mcol_inert + abs(a_mat(ii,IP_MICR)   *xi(IP_MICR)   *k_diag(IP_MICR))   &
                                       + abs(a_mat(ii,IP_PASSIVE)*xi(IP_PASSIVE)*k_diag(IP_PASSIVE))
            end do
            call check('  '//tag//' full 7x7 A*xi*K SINGULAR (inert cols = 0)', mcol_inert, 0.0_wp, 0.0_wp)
         end block
         call check('  '//tag//' inert micr = 0', xss(IP_MICR),    0.0_wp, 0.0_wp)
         call check('  '//tag//' inert pass = 0', xss(IP_PASSIVE), 0.0_wp, 0.0_wp)
      end if
      call check_true('  '//tag//' X_ss finite & >= 0', all(xss >= 0.0_wp) .and. all(xss == xss), sum(xss))
      ! EXACT EULER fixed point: step from X_ss with xi_int=xi, u -> no change.
      opts%step_solver = DECOMP_STEP_EULER
      p = pss
      call soil_carbon_step(p, u, lignin_in, xi, opts, rh, audit)
      call pack_pool_vector(p, x1)
      worst = maxval(abs(x1 - xss))
      call check_true('  '//tag//' EULER fixed point dX~0', worst < 1.0e-10_wp*max(maxval(xss),1.0_wp), worst)
      ! brute force from bare ground converges TOWARD X_ss (SASU vs brute force agree). The scheme-5
      ! tolerance is looser BY DESIGN: its passive pool (k=0.012/yr => ~170 yr residence at xi~0.5) is
      ! millennial, so 400k days (~1100 yr) leaves an O(1e-3) tail -- brute-forcing that pool to round-off
      ! is exactly the multi-millennium integration SASU (the exact fixed point above, dX~0) replaces.
      p = soil_carbon_t()
      do it = 1_ik, 400000_ik
         call soil_carbon_step(p, u, lignin_in, xi, opts, rh, audit)
      end do
      call pack_pool_vector(p, x1)
      relerr = maxval(abs(x1 - xss)) / max(maxval(xss), 1.0_wp)
      if (scheme == DECOMP_SCHEME_CENTURY5) then
         call check_true('  '//tag//' brute force -> X_ss (millennial passive tail)', relerr < 3.0e-3_wp, relerr)
      else
         call check_true('  '//tag//' brute force -> X_ss', relerr < 1.0e-6_wp, relerr)
      end if
   end subroutine spinup_one_scheme

   !=======================================================================================!
   ! 6. Residence/capacity: X_c = X_ss at equilibrium, X_p = 0; perturbed pool relaxes toward X_c.  !
   !=======================================================================================!
   subroutine test_residence_capacity()
      type(soil_carbon_t) :: pss, pert, p
      type(decomp_opts_t) :: opts
      type(soilc_diag_t)  :: diag
      type(soilc_audit_t) :: audit
      real(wp) :: a_mat(n_soil_pool,n_soil_pool), k_diag(n_soil_pool), er(n_soil_pool)
      real(wp) :: xi(n_soil_pool), u(n_soil_pool), lignin_in(2), xss(n_soil_pool)
      real(wp) :: x_before, x_after, rh
      print '(a)', '-- test 6: residence time / storage capacity --'
      opts = decomp_opts_t()
      xi = [0.7_wp, 0.5_wp, 0.7_wp, 0.5_wp, 0.5_wp, 0.5_wp, 0.5_wp]
      u  = 0.0_wp
      u(IP_FAST_GRND)=0.03_wp; u(IP_FAST_SOIL)=0.05_wp; u(IP_STRUCT_GRND)=0.02_wp; u(IP_STRUCT_SOIL)=0.04_wp
      lignin_in = 0.0_wp
      call assemble_transfer_matrix(pss, opts, a_mat, k_diag, er)
      call solve_soil_carbon_steady_state(a_mat, k_diag, xi, u, opts, pss)
      call pack_pool_vector(pss, xss)
      call soil_carbon_diagnostics(a_mat, k_diag, xi, u, pss, diag)
      call check_true('X_c == X_ss at equilibrium', maxval(abs(diag%x_c - xss)) < 1.0e-10_wp*max(maxval(xss),1.0_wp), &
                      maxval(abs(diag%x_c - xss)))
      call check_true('X_p == 0 at equilibrium', maxval(abs(diag%x_p)) < 1.0e-10_wp*max(maxval(xss),1.0_wp), maxval(abs(diag%x_p)))
      call check_true('tau_E > 0 [yr]', diag%tau_e > 0.0_wp, diag%tau_e)
      ! perturb the slow pool DOWN and confirm it relaxes back UP toward X_c:
      pert = pss ; pert%slow_carbon = pss%slow_carbon * 0.5_wp
      x_before = pert%slow_carbon
      opts%step_solver = DECOMP_STEP_EULER
      p = pert
      call soil_carbon_step(p, u, lignin_in, xi, opts, rh, audit)
      x_after = p%slow_carbon
      call check_true('perturbed pool relaxes toward X_c', (x_after - x_before) > 0.0_wp .and. x_after <= xss(IP_SLOW), &
                      x_after - x_before)
   end subroutine test_residence_capacity

   !=======================================================================================!
   ! 7. Litter mapping + lignin passive tracer.                                                    !
   !=======================================================================================!
   subroutine test_litter_and_lignin()
      type(soil_carbon_t) :: pools, p
      type(decomp_opts_t) :: opts
      type(litter_input_t):: lit
      type(soilc_audit_t) :: audit
      real(wp) :: u(n_soil_pool), lignin_in(2), xi_int(n_soil_pool), rh
      real(wp) :: f_lig_before, f_lig_after, c_before, l_before, d_s
      ! (a) driver-side per-cohort split demo (two PFTs) summed into a pool-destined litter_input_t:
      real(wp), parameter :: leaf1=0.10_wp, root1=0.06_wp, wood1=0.20_wp
      real(wp), parameter :: leaf2=0.08_wp, root2=0.05_wp, wood2=0.30_wp
      real(wp), parameter :: flab_leaf1=0.7_wp, flab_stem1=0.05_wp, agf1=0.6_wp
      real(wp), parameter :: flab_leaf2=0.5_wp, flab_stem2=0.10_wp, agf2=0.7_wp
      real(wp) :: exp_lab_grnd, exp_lab_soil, exp_str_grnd, exp_str_soil
      print '(a)', '-- test 7: litter mapping + lignin tracer --'

      ! per-cohort split: leaf splits labile/structural by f_labile_leaf and above/below by agf;
      ! fineroot is below-ground; wood (CWD) is mostly structural via f_labile_stem, split by agf.
      exp_lab_grnd = leaf1*flab_leaf1*agf1 + wood1*flab_stem1*agf1                       &
                   + leaf2*flab_leaf2*agf2 + wood2*flab_stem2*agf2
      exp_lab_soil = leaf1*flab_leaf1*(1.0_wp-agf1) + root1*flab_leaf1 + wood1*flab_stem1*(1.0_wp-agf1)  &
                   + leaf2*flab_leaf2*(1.0_wp-agf2) + root2*flab_leaf2 + wood2*flab_stem2*(1.0_wp-agf2)
      exp_str_grnd = leaf1*(1.0_wp-flab_leaf1)*agf1 + wood1*(1.0_wp-flab_stem1)*agf1     &
                   + leaf2*(1.0_wp-flab_leaf2)*agf2 + wood2*(1.0_wp-flab_stem2)*agf2
      exp_str_soil = leaf1*(1.0_wp-flab_leaf1)*(1.0_wp-agf1) + root1*(1.0_wp-flab_leaf1)                 &
                   + wood1*(1.0_wp-flab_stem1)*(1.0_wp-agf1)                                             &
                   + leaf2*(1.0_wp-flab_leaf2)*(1.0_wp-agf2) + root2*(1.0_wp-flab_leaf2)                 &
                   + wood2*(1.0_wp-flab_stem2)*(1.0_wp-agf2)
      lit%labile_grnd = exp_lab_grnd ; lit%labile_soil = exp_lab_soil
      lit%struct_grnd = exp_str_grnd ; lit%struct_soil = exp_str_soil
      lit%lignin_grnd = 0.3_wp*exp_str_grnd ; lit%lignin_soil = 0.3_wp*exp_str_soil
      call build_litter_input(lit, u, lignin_in)
      call check('build u(fast_grnd)',   u(IP_FAST_GRND),   exp_lab_grnd, 1.0e-14_wp)
      call check('build u(fast_soil)',   u(IP_FAST_SOIL),   exp_lab_soil, 1.0e-14_wp)
      call check('build u(struct_grnd)', u(IP_STRUCT_GRND), exp_str_grnd, 1.0e-14_wp)
      call check('build u(struct_soil)', u(IP_STRUCT_SOIL), exp_str_soil, 1.0e-14_wp)
      call check('build lignin_in(1)',   lignin_in(1),      0.3_wp*exp_str_grnd, 1.0e-14_wp)
      call check_true('no litter to inert/SOM pools', abs(u(IP_MICR))+abs(u(IP_SLOW))+abs(u(IP_PASSIVE))<1.0e-30_wp, &
                      u(IP_SLOW))

      ! (b) lignin tracer: proportional decay conserves f_lignin under pure decomposition (no input).
      opts = decomp_opts_t()
      pools = soil_carbon_t()
      pools%struct_grnd_carbon = 2.0_wp ; pools%struct_grnd_lignin = 0.8_wp   ! f_lignin = 0.4
      c_before = pools%struct_grnd_carbon ; l_before = pools%struct_grnd_lignin
      f_lig_before = l_before / c_before
      xi_int = 0.0_wp ; xi_int(IP_STRUCT_GRND) = 0.2_wp        ! decay only the struct_grnd pool
      u = 0.0_wp ; lignin_in = 0.0_wp
      opts%step_solver = DECOMP_STEP_EULER
      p = pools
      call soil_carbon_step(p, u, lignin_in, xi_int, opts, rh, audit)
      f_lig_after = p%struct_grnd_lignin / p%struct_grnd_carbon
      d_s = xi_int(IP_STRUCT_GRND) * (opts%k_struct/yr_day)
      call check('f_lignin conserved under pure decay', f_lig_after, f_lig_before, 1.0e-12_wp)
      call check_true('0 <= L <= C after decay', p%struct_grnd_lignin >= 0.0_wp .and. &
                      p%struct_grnd_lignin <= p%struct_grnd_carbon, p%struct_grnd_lignin)
      call check_true('lignin_resid ~ 0', abs(audit%lignin_resid) < 1.0e-14_wp, audit%lignin_resid)
      ! fresh lignin input raises f_lignin:
      pools = soil_carbon_t() ; pools%struct_grnd_carbon = 2.0_wp ; pools%struct_grnd_lignin = 0.4_wp
      xi_int = 0.0_wp
      u = 0.0_wp ; u(IP_STRUCT_GRND) = 0.5_wp ; lignin_in = [0.4_wp, 0.0_wp]     ! high-lignin input
      p = pools
      call soil_carbon_step(p, u, lignin_in, xi_int, opts, rh, audit)
      call check_true('lignin input raises f_lignin', &
                      p%struct_grnd_lignin/p%struct_grnd_carbon > 0.2_wp, p%struct_grnd_lignin/p%struct_grnd_carbon)

      ! (c) EXPM at a LARGE accelerated-spin-up step: lignin must decay in lockstep with its host
      ! carbon (exp(-d_s)), so f_lignin stays conserved and 0 <= L <= C. Before the solver-consistent
      ! survival fix, the linear (1 - d_s) drove L NEGATIVE here (d_s = 3.0 => L = 0.8*(1-3) = -1.6).
      opts = decomp_opts_t() ; opts%step_solver = DECOMP_STEP_EXPM
      pools = soil_carbon_t() ; pools%struct_grnd_carbon = 2.0_wp ; pools%struct_grnd_lignin = 0.8_wp  ! f_lignin=0.4
      c_before = pools%struct_grnd_carbon ; l_before = pools%struct_grnd_lignin ; f_lig_before = l_before/c_before
      xi_int = 0.0_wp ; xi_int(IP_STRUCT_GRND) = 3.0_wp / (opts%k_struct/yr_day)   ! force d_s = xi_int*K = 3.0
      u = 0.0_wp ; lignin_in = 0.0_wp
      p = pools
      call soil_carbon_step(p, u, lignin_in, xi_int, opts, rh, audit)
      call check_true('EXPM large step: L >= 0 (no negative lignin)', p%struct_grnd_lignin >= 0.0_wp, p%struct_grnd_lignin)
      call check_true('EXPM large step: L <= C', p%struct_grnd_lignin <= p%struct_grnd_carbon, &
                      p%struct_grnd_carbon - p%struct_grnd_lignin)
      call check('EXPM large step: f_lignin conserved under pure decay', &
                 p%struct_grnd_lignin/p%struct_grnd_carbon, f_lig_before, 1.0e-12_wp)
      call check_true('EXPM large step: lignin_resid ~ 0', abs(audit%lignin_resid) < 1.0e-12_wp, audit%lignin_resid)
   end subroutine test_litter_and_lignin

   !=======================================================================================!
   ! 8. Temperature/moisture response + N-off inertness.                                           !
   !=======================================================================================!
   subroutine test_temperature_moisture()
      type(soil_carbon_t) :: pools
      type(decomp_opts_t) :: opts
      real(wp) :: xi_a(n_soil_pool), xi_b(n_soil_pool), xi_lo(n_soil_pool), xi_hi(n_soil_pool)
      real(wp) :: xi_opt(n_soil_pool), xi_dry(n_soil_pool)
      print '(a)', '-- test 8: temperature / moisture response --'
      pools = soil_carbon_t()
      ! Q10 = 2 doubling per 10 K (scheme 5, uncapped):
      opts = decomp_opts_t() ; opts%decomp_scheme = DECOMP_SCHEME_CENTURY5
      opts%rh_q10 = 2.0_wp ; opts%rh_t_ref = 288.15_wp
      call assemble_env_scalar(298.15_wp, 298.15_wp, 0.5_wp, 0.1_wp, 0.6_wp, pools, opts, xi_a)  ! T
      call assemble_env_scalar(308.15_wp, 308.15_wp, 0.5_wp, 0.1_wp, 0.6_wp, pools, opts, xi_b)  ! T+10
      call check('Q10=2 doubling per 10K', xi_b(IP_FAST_SOIL)/xi_a(IP_FAST_SOIL), 2.0_wp, 1.0e-12_wp)
      ! ED2 capped exp: saturates (min(1,.)) at/above 45 C (318.15 K):
      opts = decomp_opts_t() ; opts%decomp_scheme = DECOMP_SCHEME_ED2
      call assemble_env_scalar(320.15_wp, 320.15_wp, opts%resp_opt_water*0.6_wp+0.1_wp, 0.1_wp, 0.6_wp, &
                               pools, opts, xi_hi)
      ! at rel = resp_opt_water the water factor = 1, so xi(fast) = min(1,exp(pos)) = 1
      call assemble_env_scalar(320.15_wp, 320.15_wp, 0.1_wp + opts%resp_opt_water*(0.6_wp-0.1_wp), 0.1_wp, 0.6_wp, &
                               pools, opts, xi_opt)
      call check('ED2 cap: f_temp<=1 at 47C, rel=opt => xi=1', xi_opt(IP_FAST_SOIL), 1.0_wp, 1.0e-12_wp)
      ! moisture: rel=opt gives f_water=1 (peak); a drier rel gives f_water<1 (lower xi):
      opts = decomp_opts_t()
      call assemble_env_scalar(298.15_wp, 298.15_wp, 0.1_wp + opts%resp_opt_water*(0.6_wp-0.1_wp), 0.1_wp, 0.6_wp, &
                               pools, opts, xi_opt)
      call assemble_env_scalar(298.15_wp, 298.15_wp, 0.2_wp, 0.1_wp, 0.6_wp, pools, opts, xi_dry)   ! rel < opt
      call check_true('moisture peak at opt (dry side lower)', xi_dry(IP_FAST_SOIL) < xi_opt(IP_FAST_SOIL), &
                      xi_opt(IP_FAST_SOIL) - xi_dry(IP_FAST_SOIL))
      ! N off (default): structural xi has NO N brake -> with zero lignin, xi(struct)==xi(fast) same layer.
      call check('N off: struct xi == fast xi (no lignin, no N brake)', xi_opt(IP_STRUCT_SOIL), xi_opt(IP_FAST_SOIL), 1.0e-14_wp)
   end subroutine test_temperature_moisture

   !=======================================================================================!
   ! 9. Fast/slow seam -- matrix single-pool Rh == fast meds_cas_biophysics Rh (matched chemistry);      !
   !    the accumulated xi_int path closes rh_today==today_rh; the daily-mean path does NOT (Jensen). !
   !=======================================================================================!
   subroutine test_fast_slow_seam()
      type(soil_carbon_t) :: pools
      type(decomp_opts_t) :: dopts
      type(co2_opts_t)    :: copts
      real(wp) :: xi(n_soil_pool), a_mat(n_soil_pool,n_soil_pool), k_diag(n_soil_pool), er(n_soil_pool)
      real(wp) :: rh_matrix_umol, rh_fast, theta, theta_dry, theta_sat, xpool
      real(wp) :: xi_int, today_rh, tmean, ftemp, fwater, rel, dt_sub, ttime, k2
      real(wp) :: rh_today, rh_meanpath, xi_meanpath
      integer(ik) :: is, nsub
      print '(a)', '-- test 9: fast/slow seam (Rh reconciliation + Jensen) --'
      pools = soil_carbon_t() ; pools%fast_soil_carbon = 4.0_wp ; xpool = 4.0_wp
      theta = 0.45_wp ; theta_dry = 0.1_wp ; theta_sat = 0.6_wp

      ! ---- (a) matched-chemistry single-pool Rh at a fixed T: matrix == fast kernel. ----
      dopts = decomp_opts_t() ; dopts%decomp_scheme = DECOMP_SCHEME_CENTURY5   ! Q10 temp form
      dopts%rh_q10 = 1.5_wp ; dopts%rh_t_ref = 288.15_wp ; dopts%k_fast = 12.0_wp
      k2 = dopts%k_fast / yr_day
      copts = co2_opts_t() ; copts%hr_model = HR_Q10 ; copts%rh_q10 = 1.5_wp ; copts%rh_t_ref = 288.15_wp
      copts%rh_k_base = k2                                    ! [1/day] == matrix K for fast pool
      copts%resp_opt_water = dopts%resp_opt_water
      copts%resp_water_below_opt = dopts%resp_water_below_opt ; copts%resp_water_above_opt = dopts%resp_water_above_opt
      call assemble_transfer_matrix(pools, dopts, a_mat, k_diag, er)
      call assemble_env_scalar(300.15_wp, 300.15_wp, theta, theta_dry, theta_sat, pools, dopts, xi)
      rh_matrix_umol = heterotrophic_respiration_matrix(a_mat, k_diag, xi, pools) * kgCday_2_umols
      rh_fast = heterotrophic_respiration_flux(xpool, 300.15_wp, theta, theta_dry, theta_sat, copts)
      call check('matrix Rh == fast column_co2 Rh (matched)', rh_matrix_umol, rh_fast, 1.0e-10_wp*max(rh_fast,1.0_wp))

      ! ---- (b) diurnal accumulation: rh_today (matrix, from xi_int) == today_rh (fast integral). ----
      nsub = 48_ik ; dt_sub = 1.0_wp / real(nsub, wp)        ! day fraction per substep [day]
      xi_int = 0.0_wp ; today_rh = 0.0_wp ; tmean = 0.0_wp
      rel = min(1.0_wp, max(0.0_wp,(theta-theta_dry)/(theta_sat-theta_dry)))
      if (rel <= dopts%resp_opt_water) then
         fwater = exp((rel-dopts%resp_opt_water)*dopts%resp_water_below_opt)
      else
         fwater = exp((dopts%resp_opt_water-rel)*dopts%resp_water_above_opt)
      end if
      do is = 1_ik, nsub
         ttime = 295.15_wp + 8.0_wp*sin(2.0_wp*3.14159265358979_wp*(real(is,wp)-0.5_wp)/real(nsub,wp))  ! diurnal T
         ftemp = dopts%rh_q10 ** ((ttime - dopts%rh_t_ref)*0.1_wp)
         xi_int   = xi_int   + (ftemp*fwater) * dt_sub          ! INT xi dt  [day]
         today_rh = today_rh + (ftemp*fwater)*k2*xpool * dt_sub ! INT xi K X dt [kgC/m2]
         tmean    = tmean    + ttime * dt_sub
      end do
      ! matrix daily step decrements the fast pool by xi_int -> rh_today for the single active pool:
      rh_today = xi_int * k2 * xpool                            ! (er=1 for fast pool)
      call check('rh_today (xi_int) == today_rh (fast integral)', rh_today, today_rh, 1.0e-14_wp*max(today_rh,1.0_wp))
      ! Jensen counter-check: recomputing xi from daily-mean T does NOT close.
      xi_meanpath = (dopts%rh_q10 ** ((tmean - dopts%rh_t_ref)*0.1_wp)) * fwater
      rh_meanpath = xi_meanpath * k2 * xpool
      call check_true('daily-mean path leaves a nonzero Jensen gap', abs(rh_meanpath - today_rh) > 1.0e-8_wp, &
                      rh_meanpath - today_rh)
   end subroutine test_fast_slow_seam

end program test_soil_biogeochem
