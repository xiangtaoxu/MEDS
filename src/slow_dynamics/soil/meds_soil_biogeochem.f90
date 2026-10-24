! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_soil_biogeochem -- the soil carbon pools' daily step (docs/science/soil_carbon.md). !
! The per-patch decomposable soil carbon obeys the CENTURY-family carbon matrix ODE (eq. (1))    !
!                                                                                          !
!         dX/dt = B*I(t) + A*xi(t)*K*X(t)        [kgC/m2/day]                                       !
!                                                                                          !
! -- litter/CWD inputs `u = B*I` feed a multi-pool state X; each donor decays at its own rate         !
! xi_j*K_j; the decayed carbon is split (fraction er_j respired to CO2, fractions a_ij transferred      !
! to downstream pools). This is ED2's `update_C_and_N_pools` CENTURY, REORGANIZED as a matrix so a      !
! semi-analytic accelerated spin-up (SASU) and traceability diagnostics fall out of the same             !
! assembled (A, xi, K).                                                                                   !
! The decay and respiration chemistry it integrates (xi, A, K, Rh) is meds_respiration's.  !
!                                                                                          !
! LOAD-BEARING CORRECTNESS RULES (docs/science/soil_carbon.md §3, §4, §6, §8):                            !
!   * Scalar placement is A*(xi*K*X), NOT xi*(A*K*X): each donor decays at ITS OWN xi_j*K_j, and ONLY     !
!     THEN is the loss split. xi and K are diagonal (commute), but neither commutes with A.               !
!   * a_jj = -1 UNIFORMLY for every pool (all decayed carbon leaves the donor). er_j = 1 - sum_{i!=j}a_ij.!
!     Rh = sum_j er_j*xi_j*K_j*X_j = -1^T * A*xi*K*X  (the respired complement A does NOT transfer).       !
!   * The steady-state / diagnostic solve runs on the ACTIVE sub-block {j: K_j>0} ONLY -- the full 7x7    !
!     A*xi*K is SINGULAR when inert pools have K=0 (the DEFAULT scheme). Inert pools are held at 0.        !
!   * soil_carbon_step decrements each donor by the fast loop's ACCUMULATED per-pool loss integral         !
!     xi_int_j = INT_day xi_j dt [day] (NOT xi recomputed from daily-mean T/theta) -- so the pool loss     !
!     equals the day's integrated atmospheric Rh BY CONSTRUCTION (closes the Jensen gap, §6).              !
!                                                                                          !
! All compute kernels take bare scalars/arrays + small derived types and link the shared,  !
! state and config libraries and meds_respiration; they are `pure` where the arithmetic    !
! allows and allocation/LAPACK-free (the active-block solve is meds_linear_algebra's dense !
! LU), so they are GPU-eligible.                                                           !
!==========================================================================================!
module meds_soil_biogeochem
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : yr_day, tiny_num
   use meds_linear_algebra,   only : matrix_exp, lu_factor, lu_solve
   use meds_biogeochem_types, only : litter_input_t, soilc_audit_t, soilc_diag_t
   use meds_patch_state_vars,   only : soil_carbon_t, pack_pool_vector, unpack_pool_vector
   use meds_respiration,        only : decomposition_transfer_matrix, heterotrophic_respiration
   use meds_biogeochem_params, only : soil_carbon_params_t, DECOMP_STEP_EULER, DECOMP_STEP_EXPM,                            &
                                      DECOMP_SCHEME_CENTURY5, N_SOIL_POOL, IP_FAST_GRND, IP_FAST_SOIL,                      &
                                      IP_STRUCT_GRND, IP_STRUCT_SOIL, IP_MICR, IP_SLOW, IP_PASSIVE
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
   implicit none
   private

   public :: build_litter_input, soil_carbon_step
   public :: solve_soil_carbon_steady_state, soil_carbon_diagnostics
   public :: soil_carbon_bad_pool, soil_carbon_pool_name, SOILC_MAX_PLAUSIBLE, SOILC_NEG_TOL

   !=======================================================================================!
   !  PLAUSIBILITY, as distinct from CONSERVATION.                                           !
   !                                                                                          !
   !  A closed budget proves bookkeeping, not physics: a runaway that CONSERVES carbon closes    !
   !  every phase and every currency of the slow ledger while moving absurd amounts of it        !
   !  (docs/science/soil_carbon.md §9).                                                          !
   !                                                                                          !
   !  These two numbers are what "physically impossible" means for a CENTURY pool:                !
   !    * a carbon pool is a MASS. It is never negative. The tolerance is round-off only --       !
   !      the daily Euler step's own arithmetic, not a modelling allowance.                       !
   !    * 1e4 kgC/m2 is ~500x the richest real soil profile (a temperate forest carries 20-40).   !
   !      Nothing between "plausible" and this ceiling is asserted on; the ceiling exists to      !
   !      catch divergence, not to police science.                                                !
   !=======================================================================================!
   real(wp), parameter :: SOILC_NEG_TOL       = 1.0e-9_wp   !< [kgC/m2] round-off allowance below zero
   real(wp), parameter :: SOILC_MAX_PLAUSIBLE = 1.0e4_wp    !< [kgC/m2] per-pool divergence ceiling

contains

   !=======================================================================================!
   !  Map the ALREADY-PARTITIONED, pool-destined litter fluxes onto the input vector u = B*I. The   !
   !  per-PFT labile/structural (f_labile) and above/below (agf) split is applied PER-COHORT by the   !
   !  driver BEFORE summing to the patch (docs/science/soil_carbon.md §5); this kernel does NO          !
   !  f_labile arithmetic -- it just places the four pre-split streams (+ lignin) into the pool         !
   !  vector. CWD -> structural is part of the driver split (wood uses f_labile_stem).                  !
   !=======================================================================================!
   pure subroutine build_litter_input(litter, u, lignin_in)
      type(litter_input_t), intent(in)  :: litter
      real(wp),             intent(out) :: u(N_SOIL_POOL)   !< [kgC/m2/day] input to each pool (1-4 nonzero)
      real(wp),             intent(out) :: lignin_in(2)     !< [kgC/m2/day] lignin flux to struct_grnd/struct_soil
      u = 0.0_wp
      u(IP_FAST_GRND)   = litter%labile_grnd
      u(IP_FAST_SOIL)   = litter%labile_soil
      u(IP_STRUCT_GRND) = litter%struct_grnd
      u(IP_STRUCT_SOIL) = litter%struct_soil
      lignin_in(1) = litter%lignin_grnd
      lignin_in(2) = litter%lignin_soil
   end subroutine build_litter_input

   !=======================================================================================!
   !  Advance the pools ONE DAILY step: X_new = X + u_day + A*diag(xi_int*K)*X, decrementing each      !
   !  donor by the TRUE sub-daily loss integral. xi_int_j = INT_day xi_j dt [day] (the fast loop's       !
   !  accumulator), NOT xi recomputed from daily-mean T/theta -- so the pool loss equals the day's        !
   !  integrated atmospheric Rh by construction (docs/science/soil_carbon.md eqs. (5)-(6)). Two            !
   !  frozen-operator solvers:                                                                             !
   !    DECOMP_STEP_EULER -- ED2-faithful forward Euler (one step; daily k*dt << 1, already stable).       !
   !    DECOMP_STEP_EXPM  -- exact matrix exponential of the frozen operator (EXACTNESS for LARGE           !
   !                         accelerated/spin-up steps; via an augmented (n+1) exponential so no A inverse).!
   !  rh_today is reported SOLVER-CONSISTENTLY as the actual pool loss to atmosphere, rh_today =            !
   !  litter_in - dC_pool (inter-pool transfers cancel in the total), so it matches the pool decrement      !
   !  for BOTH solvers and equals the fast loop's accumulated Rh. It is an AUDIT value -- NOT re-fed     !
   !  to CAS.                                                                                            !
   !  Structural lignin decays with its host (L_s *= 1 - d_s, or exp(-d_s) under EXPM) + lignin_in          !
   !  (a passive tracer, its OWN balance, outside the carbon mass vector). dt_day is one day (u in          !
   !  [kgC/m2/day], xi_int in [day]).                                                                       !
   !=======================================================================================!
   subroutine soil_carbon_step(pools, u, lignin_in, xi_int, opts, rh_today, audit)
      type(soil_carbon_t),        intent(inout) :: pools
      real(wp),                   intent(in)    :: u(N_SOIL_POOL), lignin_in(2)
      real(wp),                   intent(in)    :: xi_int(N_SOIL_POOL)   !< [day] INT_day xi_j dt (fast accumulator)
      type(soil_carbon_params_t), intent(in)    :: opts
      real(wp),                   intent(out)   :: rh_today              !< [kgC/m2/day] actual pool loss = fast loop's Rh
      type(soilc_audit_t),        intent(out)   :: audit
      real(wp) :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL), er(N_SOIL_POOL)
      real(wp) :: x0(N_SOIL_POOL), x1(N_SOIL_POOL), dvec(N_SOIL_POOL), loss(N_SOIL_POOL)
      real(wp) :: lig0(2), lig1(2), d_s(2), surv(2), litter_in, dc_pool, lig_resid
      integer(ik) :: i, j

      call decomposition_transfer_matrix(pools, opts, a_mat, k_diag, er)
      call pack_pool_vector(pools, x0)
      lig0(1) = pools%struct_grnd_lignin
      lig0(2) = pools%struct_soil_lignin

      !----- per-pool day-integrated decay FRACTION d_j = xi_int_j * K_j [-] (dimensionless). ------!
      do j = 1_ik, N_SOIL_POOL
         dvec(j) = xi_int(j) * k_diag(j)
      end do

      select case (opts%step_solver)
      case (DECOMP_STEP_EXPM)
         call step_expm(a_mat, dvec, u, x0, x1)
      case default   ! DECOMP_STEP_EULER
         do j = 1_ik, N_SOIL_POOL
            loss(j) = dvec(j) * x0(j)                 ! [kgC/m2] day-integrated loss of donor j
         end do
         do i = 1_ik, N_SOIL_POOL
            x1(i) = x0(i) + u(i)
            do j = 1_ik, N_SOIL_POOL
               x1(i) = x1(i) + a_mat(i, j) * loss(j)  ! a_ii=-1 removes donor i's own loss; a_ij adds transfers
            end do
         end do
      end select

      !----- lignin passive tracer: it decays in LOCKSTEP with its host structural CARBON, so         !
      !      f_lignin is conserved under pure decomposition and 0 <= L <= C -- which requires the      !
      !      SOLVER-CONSISTENT survival fraction. EULER carbon decays by the linear (1 - d_s) (the      !
      !      docs/science/soil_carbon.md §6 law); EXPM decays a pure-donor structural pool by the      !
      !      exact exp(-d_s), so lignin must too. Using (1 - d_s) under EXPM would drift f_lignin      !
      !      and drive L NEGATIVE for the large d_s of accelerated spin-up steps, which the daily      !
      !      EULER step never shows.                                                                   !
      d_s(1) = dvec(IP_STRUCT_GRND)
      d_s(2) = dvec(IP_STRUCT_SOIL)
      select case (opts%step_solver)
      case (DECOMP_STEP_EXPM) ; surv(1) = exp(-d_s(1))   ; surv(2) = exp(-d_s(2))
      case default            ; surv(1) = 1.0_wp - d_s(1) ; surv(2) = 1.0_wp - d_s(2)
      end select
      lig1(1) = lig0(1) * surv(1) + lignin_in(1)
      lig1(2) = lig0(2) * surv(2) + lignin_in(2)

      !----- commit the new state. ----------------------------------------------------------------!
      call unpack_pool_vector(x1, pools)
      pools%struct_grnd_lignin = lig1(1)
      pools%struct_soil_lignin = lig1(2)

      !----- carbon-mass audit (both solvers): rh_today = litter_in - dC_pool (transfers cancel). --!
      litter_in = sum(u)
      dc_pool   = sum(x1) - sum(x0)
      rh_today  = litter_in - dc_pool
      audit%litter_in_matrix = litter_in
      audit%rh_out       = rh_today
      audit%dC_pool      = dc_pool
      !----- No pool-mass residual is reported: with rh_today DEFINED just above as                   !
      !      `litter_in - dC_pool`, `dC_pool - (litter_in - rh_today)` is identically zero. The       !
      !      independent check (recompute Rh from er*xi_int*K*X0 and compare) lives in                !
      !      test_soil_biogeochem and is the one that can catch a dropped-litter or                   !
      !      misattributed-transfer bug (docs/science/soil_carbon.md §9).                             !
      !----- THE FREEZE NUMBER: the largest fraction of any pool this step withdraws. See its comment !
      !      in soilc_audit_t for why this is the number that says whether freezing the pool across   !
      !      the slow step is sound, and how it differs from rh_seam_gap.  ---------------------------!
      audit%lambda_max   = 0.0_wp
      audit%lambda_pool  = 0_ik
      do i = 1_ik, N_SOIL_POOL
         if (dvec(i) > audit%lambda_max) then
            audit%lambda_max  = dvec(i)
            audit%lambda_pool = i
         end if
      end do
      audit%rh_fast_accum= rh_today                                  ! the driver overwrites it with the fast loop's Rh
      audit%rh_seam_gap  = 0.0_wp                                    ! rh_out - rh_fast_accum (set by the driver)
      !----- lignin passive-tracer balance: dL_s == lignin_in_s - (1-surv_s)*L_s, 0 by construction   !
      !      for EITHER solver (EULER loss (1-surv)=d_s recovers d_s*L_s). ----------------------------!
      lig_resid = 0.0_wp
      do i = 1_ik, 2_ik
         lig_resid = max(lig_resid, abs( (lig1(i) - lig0(i)) - (lignin_in(i) - (1.0_wp - surv(i)) * lig0(i)) ))
      end do
      audit%lignin_resid = lig_resid
   end subroutine soil_carbon_step

   !=======================================================================================!
   !  SEMI-ANALYTIC SPIN-UP (SASU; docs/science/soil_carbon.md eq. (7)). Steady state under           !
   !  climatological-mean drivers is a direct linear                                                  !
   !  solve, X_ss = -(A*xi_bar*K)^{-1} * B*I_bar, over the ACTIVE sub-block {j: K_j>0} ONLY -- the full !
   !  N_SOIL_POOL system is SINGULAR (inert pools have K=0). Inert pools are held at their conserved      !
   !  value (0: no litter input, no transfer into them in the default scheme). The reduced block is       !
   !  solved by LU factorisation with partial pivoting (meds_linear_algebra), which handles the            !
   !  non-triangular scheme-5 back-transfers as well as the default block.                                 !
   !=======================================================================================!
   pure subroutine solve_soil_carbon_steady_state(a_mat, k_diag, xi_bar, u_bar, opts, pools_ss)
      real(wp),                   intent(in)  :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL)
      real(wp),                   intent(in)  :: xi_bar(N_SOIL_POOL), u_bar(N_SOIL_POOL)
      type(soil_carbon_params_t), intent(in)  :: opts
      type(soil_carbon_t),        intent(out) :: pools_ss
      real(wp)    :: x_ss(N_SOIL_POOL)
      !----- opts is not read: the active set {K_j>0} is derived inside steady_state_vector,           !
      !      which already encodes the scheme via k_micr/k_pass = 0. It stays in the interface for     !
      !      an iterated scheme-5 steady state. --------------------------------------------------------!
      call steady_state_vector(a_mat, k_diag, xi_bar, u_bar, x_ss)
      pools_ss = soil_carbon_t()                       ! zero all (inert pools + lignin stay 0)
      call unpack_pool_vector(x_ss, pools_ss)
   end subroutine solve_soil_carbon_steady_state

   !----- Solve M[active]*x = -u[active] on the active sub-block; inert indices held at 0. --------!
   pure subroutine steady_state_vector(a_mat, k_diag, xi_bar, u_bar, x_ss)
      real(wp), intent(in)  :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL)
      real(wp), intent(in)  :: xi_bar(N_SOIL_POOL), u_bar(N_SOIL_POOL)
      real(wp), intent(out) :: x_ss(N_SOIL_POOL)
      real(wp)    :: xik(N_SOIL_POOL), m_red(N_SOIL_POOL, N_SOIL_POOL), b_red(N_SOIL_POOL)
      real(wp)    :: x_red(N_SOIL_POOL)
      integer(ik) :: idx(N_SOIL_POOL), piv(N_SOIL_POOL), na, i, j
      logical     :: ok

      xik   = xi_bar * k_diag                             ! diag(xi*K)
      m_red = 0.0_wp
      b_red = 0.0_wp
      idx   = 0_ik
      !----- active index set {j: K_j > 0}. -------------------------------------------------------!
      na = 0_ik
      do j = 1_ik, N_SOIL_POOL
         if (k_diag(j) > tiny_num) then
            na = na + 1_ik
            idx(na) = j
         end if
      end do
      !----- reduced block M = (A*diag(xi*K))[active,active], RHS b = -(u)[active]. ----------------!
      do j = 1_ik, na
         b_red(j) = -u_bar(idx(j))
         do i = 1_ik, na
            m_red(i, j) = a_mat(idx(i), idx(j)) * xik(idx(j))
         end do
      end do
      x_red = b_red
      call lu_factor(m_red(1:na, 1:na), piv(1:na), na, ok)
      if (.not. ok) error stop 'meds_soil_biogeochem: the steady-state soil carbon system is singular'
      call lu_solve(m_red(1:na, 1:na), piv(1:na), na, x_red(1:na))
      !----- scatter back (inert pools stay 0). ---------------------------------------------------!
      x_ss = 0.0_wp
      do j = 1_ik, na
         x_ss(idx(j)) = x_red(j)
      end do
   end subroutine steady_state_vector

   !=======================================================================================!
   !  TRACEABILITY DIAGNOSTICS (pure post-processing off the assembled matrices;                     !
   !  docs/science/soil_carbon.md eq. (8)):                                                          !
   !    X_c = -(A*xi*K)^{-1}*u    storage CAPACITY (equilibrium chased under current climate+input)     !
   !    X_p = X_c - X             storage POTENTIAL (remaining sink; <0 = a source)                     !
   !    tau_E = sum(X_c)/sum(u)   ecosystem residence time [yr] (stock/throughput; u in /day -> /yr).    !
   !  All on the ACTIVE sub-block (inert pools carry no residence time; held at 0).                      !
   !=======================================================================================!
   pure subroutine soil_carbon_diagnostics(a_mat, k_diag, xi, u, pools, diag)
      real(wp),            intent(in)  :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL)
      real(wp),            intent(in)  :: xi(N_SOIL_POOL), u(N_SOIL_POOL)
      type(soil_carbon_t), intent(in)  :: pools
      type(soilc_diag_t),  intent(out) :: diag
      real(wp) :: x_c(N_SOIL_POOL), x(N_SOIL_POOL), throughput
      call steady_state_vector(a_mat, k_diag, xi, u, x_c)        ! X_c = -(A xi K)^{-1} u
      call pack_pool_vector(pools, x)
      diag%x_c  = x_c
      diag%x_p  = x_c - x
      throughput = sum(u)
      if (throughput > tiny_num) then
         diag%tau_e = (sum(x_c) / throughput) / yr_day           ! [kgC/m2] / [kgC/m2/day] = day -> yr
      else
         diag%tau_e = 0.0_wp
      end if
      diag%rh = heterotrophic_respiration(a_mat, k_diag, xi, pools)
   end subroutine soil_carbon_diagnostics


   !=======================================================================================!
   !  EXPM daily step via an AUGMENTED matrix exponential. For constant operator M = A*diag(xi_int*K)  !
   !  and constant day-input u, the exact evolution of [X;1] is exp(Aug)*[X0;1] with                    !
   !  Aug = [[M, u],[0, 0]] -- this handles the constant input WITHOUT inverting M (safe for the         !
   !  singular inert-pool rows). exp is meds_linear_algebra's matrix_exp (scaling and squaring +           !
   !  Taylor series) on the (n+1)x(n+1) augmented matrix (fixed size, allocation-free, GPU-safe).         !
   !  Mass closes by construction since rh_today is reported as litter_in - dC_pool.                      !
   !  docs/science/soil_carbon.md eq. (6).                                                                !
   !=======================================================================================!
   pure subroutine step_expm(a_mat, dvec, u, x0, x1)
      real(wp), intent(in)  :: a_mat(N_SOIL_POOL, N_SOIL_POOL)
      real(wp), intent(in)  :: dvec(N_SOIL_POOL)          !< day-integrated decay fraction xi_int*K [-]
      real(wp), intent(in)  :: u(N_SOIL_POOL), x0(N_SOIL_POOL)
      real(wp), intent(out) :: x1(N_SOIL_POOL)
      integer(ik), parameter :: M = N_SOIL_POOL + 1_ik
      real(wp)    :: aug(M, M), e(M, M), y0(M), y1(M)
      integer(ik) :: i, j

      !----- assemble Aug = [[A*diag(dvec), u],[0,0]]. ---------------------------------------------!
      aug = 0.0_wp
      do j = 1_ik, N_SOIL_POOL
         do i = 1_ik, N_SOIL_POOL
            aug(i, j) = a_mat(i, j) * dvec(j)
         end do
      end do
      do i = 1_ik, N_SOIL_POOL
         aug(i, M) = u(i)
      end do
      call matrix_exp(aug, e, M)
      !----- y1 = e * [x0; 1]. --------------------------------------------------------------------!
      y0(1:N_SOIL_POOL) = x0
      y0(M)             = 1.0_wp
      do i = 1_ik, M
         y1(i) = 0.0_wp
         do j = 1_ik, M
            y1(i) = y1(i) + e(i, j) * y0(j)
         end do
      end do
      x1 = y1(1:N_SOIL_POOL)
   end subroutine step_expm

   !---------------------------------------------------------------------------------------!
   ! soil_carbon_bad_pool -- 0 when every CENTURY pool is physically possible, else the index of  !
   ! the FIRST pool that is not. Checks, in order: not-a-number, negative beyond round-off, and    !
   ! above the divergence ceiling. Cheap enough (7 comparisons) to call at every operator boundary. !
   !                                                                                          !
   ! It deliberately reports an INDEX rather than a logical: the caller's message can then name    !
   ! the pool, and which pool went bad is the first thing that narrows down which operator did it. !
   !---------------------------------------------------------------------------------------!
   pure function soil_carbon_bad_pool(pools) result(k)
      type(soil_carbon_t), intent(in) :: pools
      integer(ik)                     :: k
      real(wp)    :: x(N_SOIL_POOL)
      integer(ik) :: j
      call pack_pool_vector(pools, x)
      k = 0_ik
      do j = 1_ik, N_SOIL_POOL
         !----- NaN first, and alone, by classification. A comparison with a NaN can raise IEEE      !
         !      invalid -- ordered ones always, and ifx compiles even `x /= x` that way at -O0 --     !
         !      which a trapping build (-fpe0, -ffpe-trap=invalid) stops on, and `.or.` does not    !
         !      short-circuit. ieee_is_nan is elemental, so this stays `pure`. ----------------------!
         if (ieee_is_nan(x(j))) then
            k = j
            return
         end if
         if (x(j) < -SOILC_NEG_TOL .or. x(j) > SOILC_MAX_PLAUSIBLE) then
            k = j
            return
         end if
      end do
   end function soil_carbon_bad_pool

   !----- Name for the index soil_carbon_bad_pool returns, so a failure message is readable. -----!
   pure function soil_carbon_pool_name(k) result(nm)
      integer(ik), intent(in) :: k
      character(len=16)       :: nm
      select case (k)
      case (IP_FAST_GRND)   ; nm = 'fast_grnd'
      case (IP_FAST_SOIL)   ; nm = 'fast_soil'
      case (IP_STRUCT_GRND) ; nm = 'struct_grnd'
      case (IP_STRUCT_SOIL) ; nm = 'struct_soil'
      case (IP_MICR)        ; nm = 'microbial'
      case (IP_SLOW)        ; nm = 'slow'
      case (IP_PASSIVE)     ; nm = 'passive'
      case default          ; nm = '(none)'
      end select
   end function soil_carbon_pool_name

end module meds_soil_biogeochem
