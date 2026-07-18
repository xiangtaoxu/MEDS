!==========================================================================================!
! meds_plant_capi -- a thin ISO_C_BINDING shim that exposes the leaf-physiology model to C /  !
! Python (ctypes). It is NOT part of libmeds_plant; it is compiled only into the    !
! optional shared library libmeds_leaf_c (CMake -DMEDS_BUILD_PYLIB=ON), so the core model      !
! stays foreign-call-free.                                                                     !
!                                                                                          !
! Three C-interoperable derived types mirror the model's leaf_env_t / leaf_photo_params_t /    !
! leaf_flux_t (same field ORDER; real64 -> c_double, int32 -> c_int, the flux `converged`       !
! logical -> c_int 0/1). Three exported procedures wrap the model:                              !
!   * meds_leaf_solve         -- the coupled A-gs-Ci solver (solve_leaf_gas_exchange).           !
!   * meds_assimilation_demand_c3    -- the raw C3 FvCB demand kernel at a prescribed Ci (assimilation_demand_c3):!
!                                gross A and the Ac/Aj/Ap limitation rates, NO temperature scaling. !
!   * meds_electron_transport_j -- J from Jmax and light (the non-rectangular hyperbola).          !
!   * meds_peaked_arrhenius   -- the peaked temperature-response function.                        !
!   * meds_arrhenius          -- the plain Arrhenius temperature-response function.               !
! Composing the last four (kinetics via meds_arrhenius, J via meds_electron_transport_j, then       !
! meds_assimilation_demand_c3) draws an A-Ci demand curve from Python using Vcmax/Jmax DIRECTLY (no        !
! capacity temperature-correction). The Python side (python/meds/leaf/_ffi.py) mirrors the structs. !
!==========================================================================================!
module meds_plant_capi
   use iso_c_binding,            only : c_double, c_int
   use meds_kinds,              only : wp, ik
   use meds_plant_types,         only : leaf_env_t, leaf_photo_params_t, leaf_flux_t
   use meds_leaf_gas_exchange,        only : solve_leaf_gas_exchange
   use meds_leaf_gas_exchange,only : assimilation_demand_c3, electron_transport_j
   use meds_temp_response, only : peaked_arrhenius_scale, arrhenius_scale
   implicit none
   private

   public :: leaf_env_c, leaf_params_c, leaf_flux_c, leaf_c3_demand_c
   public :: meds_leaf_solve, meds_assimilation_demand_c3, meds_electron_transport_j
   public :: meds_peaked_arrhenius, meds_arrhenius

   !----- C-interoperable mirror of leaf_env_t (8 doubles). --------------------------------!
   type, bind(c) :: leaf_env_c
      real(c_double) :: par, leaf_temp, vpd, ca, pressure, psi_leaf, gb, psi_soil
   end type leaf_env_c

   !----- C-interoperable mirror of leaf_flux_t (7 doubles + 2 ints; converged 0/1). --------!
   type, bind(c) :: leaf_flux_c
      real(c_double) :: A_net, A_gross, gs, ci, cs, transpiration, rd
      integer(c_int) :: limitation
      integer(c_int) :: converged
   end type leaf_flux_c

   !----- C-interoperable mirror of leaf_photo_params_t (1 int + 35 doubles, same order). ---!
   type, bind(c) :: leaf_params_c
      integer(c_int) :: pathway
      real(c_double) :: vcmax25, jmax25, tpu25, rd25, kp25
      real(c_double) :: g0, g1, d0, quantum_yield, theta_j, theta_cj, theta_ic
      real(c_double) :: lambda25, psi_open, psi_close, lambda_psi_exp, sref_stomata
      real(c_double) :: kc25, ko25, gstar25
      real(c_double) :: ea_kc, ea_ko, ea_gstar, ea_vcmax, ea_jmax, ea_rd
      real(c_double) :: hd_vcmax, hd_jmax, hd_rd, ds_vcmax, ds_jmax, ds_rd
      real(c_double) :: o2_mol_frac, absorptance, phi_psii
   end type leaf_params_c

   !----- C-interoperable C3 demand rates at a prescribed Ci (4 doubles). -------------------!
   type, bind(c) :: leaf_c3_demand_c
      real(c_double) :: A_gross, Ac, Aj, Ap
   end type leaf_c3_demand_c

contains

   !---------------------------------------------------------------------------------------!
   ! Coupled leaf gas-exchange: unpack the C structs into the model types, solve, pack back. !
   ! sm/tresp/colim are the SM_*/TRESP_*/COLIM_* integer codes; use_boundary_layer is 0/1.               !
   !---------------------------------------------------------------------------------------!
   subroutine meds_leaf_solve(env_c, p_c, sm, tresp, colim, use_boundary_layer, flux_c) bind(c, name="meds_leaf_solve")
      type(leaf_env_c),    intent(in)  :: env_c
      type(leaf_params_c), intent(in)  :: p_c
      integer(c_int), value, intent(in) :: sm, tresp, colim, use_boundary_layer
      type(leaf_flux_c),   intent(out) :: flux_c
      type(leaf_env_t)          :: env
      type(leaf_photo_params_t) :: p
      type(leaf_flux_t)         :: flux

      env = to_env(env_c) ; p = to_params(p_c)

      call solve_leaf_gas_exchange(env, p, int(sm), int(tresp), int(colim), use_boundary_layer /= 0_c_int, flux)

      flux_c%A_net = flux%A_net ; flux_c%A_gross = flux%A_gross ; flux_c%gs = flux%gs
      flux_c%ci = flux%ci ; flux_c%cs = flux%cs ; flux_c%transpiration = flux%transpiration
      flux_c%rd = flux%rd ; flux_c%limitation = int(flux%limitation, c_int)
      flux_c%converged = merge(1_c_int, 0_c_int, flux%converged)
   end subroutine meds_leaf_solve

   !---------------------------------------------------------------------------------------!
   ! Raw C3 FvCB demand at a PRESCRIBED intercellular CO2 (assimilation_demand_c3), stomata bypassed !
   ! and NO temperature scaling: the caller passes already-in-situ values -- vcmax, j (the      !
   ! electron-transport RATE, from meds_electron_transport_j), tpu, and the mole-fraction        !
   ! kinetics gstar/kc/ko/o2 [umol/mol]. Returns gross A and the Ac/Aj/Ap limitation rates      !
   ! (net = gross - Rd is the caller's business). colim is a COLIM_* code; theta is the C3        !
   ! co-limitation curvature (use COLIM_MIN for a sharp min(Ac,Aj,Ap) envelope).                 !
   !---------------------------------------------------------------------------------------!
   subroutine meds_assimilation_demand_c3(ci, vcmax, j, tpu, gstar, kc, ko, o2, colim, theta, dem_c)  &
                                   bind(c, name="meds_assimilation_demand_c3")
      real(c_double), value, intent(in) :: ci, vcmax, j, tpu, gstar, kc, ko, o2, theta
      integer(c_int), value, intent(in) :: colim
      type(leaf_c3_demand_c), intent(out) :: dem_c
      real(wp) :: a_gross, ac, aj, ap
      call assimilation_demand_c3(real(ci, wp), real(vcmax, wp), real(j, wp), real(tpu, wp),           &
                           real(gstar, wp), real(kc, wp), real(ko, wp), real(o2, wp),           &
                           int(colim), real(theta, wp), a_gross, ac, aj, ap)
      dem_c%A_gross = a_gross ; dem_c%Ac = ac ; dem_c%Aj = aj ; dem_c%Ap = ap
   end subroutine meds_assimilation_demand_c3

   !---------------------------------------------------------------------------------------!
   ! Electron-transport rate J from Jmax and incident PAR (non-rectangular hyperbola).      !
   !---------------------------------------------------------------------------------------!
   function meds_electron_transport_j(par, absorptance, phi_psii, jmax, theta) result(j)       &
                                      bind(c, name="meds_electron_transport_j")
      real(c_double), value, intent(in) :: par, absorptance, phi_psii, jmax, theta
      real(c_double)                    :: j
      j = electron_transport_j(par, absorptance, phi_psii, jmax, theta)
   end function meds_electron_transport_j

   !---------------------------------------------------------------------------------------!
   ! Peaked-Arrhenius temperature response (k25 anchored; Ea/Hd/dS in J/mol, T in K).       !
   !---------------------------------------------------------------------------------------!
   function meds_peaked_arrhenius(k25, ea, hd, ds, t_leaf) result(y) bind(c, name="meds_peaked_arrhenius")
      real(c_double), value, intent(in) :: k25, ea, hd, ds, t_leaf
      real(c_double)                    :: y
      y = peaked_arrhenius_scale(k25, ea, hd, ds, t_leaf)
   end function meds_peaked_arrhenius

   !---------------------------------------------------------------------------------------!
   ! Plain Arrhenius temperature response.                                                 !
   !---------------------------------------------------------------------------------------!
   function meds_arrhenius(k25, ea, t_leaf) result(y) bind(c, name="meds_arrhenius")
      real(c_double), value, intent(in) :: k25, ea, t_leaf
      real(c_double)                    :: y
      y = arrhenius_scale(k25, ea, t_leaf)
   end function meds_arrhenius

   !---------------------------------------------------------------------------------------!
   ! Unpack the C mirror structs into the model's derived types (used by meds_leaf_solve).   !
   !---------------------------------------------------------------------------------------!
   pure function to_env(env_c) result(env)
      type(leaf_env_c), intent(in) :: env_c
      type(leaf_env_t)             :: env
      env = leaf_env_t(env_c%par, env_c%leaf_temp, env_c%vpd, env_c%ca, env_c%pressure,        &
                       env_c%psi_leaf, env_c%gb, env_c%psi_soil)
   end function to_env

   pure function to_params(p_c) result(p)
      type(leaf_params_c), intent(in) :: p_c
      type(leaf_photo_params_t)       :: p
      p = leaf_photo_params_t(p_c%pathway, p_c%vcmax25, p_c%jmax25, p_c%tpu25, p_c%rd25,       &
             p_c%kp25, p_c%g0, p_c%g1, p_c%d0, p_c%quantum_yield, p_c%theta_j, p_c%theta_cj,   &
             p_c%theta_ic, p_c%lambda25, p_c%psi_open, p_c%psi_close, p_c%lambda_psi_exp,      &
             p_c%sref_stomata,                                                                 &
             p_c%kc25, p_c%ko25, p_c%gstar25, p_c%ea_kc, p_c%ea_ko, p_c%ea_gstar,             &
             p_c%ea_vcmax, p_c%ea_jmax, p_c%ea_rd, p_c%hd_vcmax, p_c%hd_jmax, p_c%hd_rd,       &
             p_c%ds_vcmax, p_c%ds_jmax, p_c%ds_rd, p_c%o2_mol_frac, p_c%absorptance, p_c%phi_psii)
   end function to_params

end module meds_plant_capi
