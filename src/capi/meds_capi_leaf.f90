!==========================================================================================!
! meds_capi_leaf -- the C-API shim for LEAF GAS EXCHANGE (`meds.plant.leaf`).                 !
!                                                                                          !
! One shim per subsystem, mirroring the Fortran tree (structure plan §7.6 #3). Split out of  !
! the former `meds_plant_capi`, which carried leaf and phenology together and so had to be    !
! rebuilt, re-reviewed and re-tested as a unit whenever either moved.                          !
!                                                                                          !
! EVERY bind(c) struct here is an ABI CONTRACT with `python/meds/plant/_ffi.py`: the field     !
! ORDER must match, member for member. It is compiled by a mandatory ctest target              !
! (`test_capi_leaf`) precisely so a mismatch is a BUILD failure -- issue #95 -> #100 was a      !
! component inserted mid-type in `leaf_photo_params_t` that broke the C API while the whole     !
! suite stayed green, because this file was compiled only by an optional target.                 !
!==========================================================================================!
module meds_capi_leaf
   use iso_c_binding,          only : c_double, c_int
   use meds_kinds,             only : wp, ik
   use meds_plant_types,       only : leaf_env_t, leaf_photo_params_t, leaf_flux_t
   use meds_leaf_gas_exchange, only : solve_leaf_gas_exchange
   use meds_leaf_gas_exchange, only : assimilation_demand_c3, electron_transport_j
   use meds_temp_response,     only : peaked_arrhenius_scale, arrhenius_scale
   implicit none
   private

   public :: leaf_env_c, leaf_params_c, leaf_flux_c, leaf_c3_demand_c
   public :: meds_leaf_solve, meds_assimilation_demand_c3, meds_electron_transport_j
   public :: meds_peaked_arrhenius, meds_arrhenius

   !----- C-interoperable mirror of leaf_env_t (8 doubles). --------------------------------!
   type, bind(c) :: leaf_env_c
      real(c_double) :: par, leaf_temp, vpd, ca, pressure, psi_leaf, gb, psi
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
      !----- KEYWORD form, deliberately. A positional structure constructor silently re-binds every  !
      !      value when a component is inserted mid-type, and components WITH defaults cannot be      !
      !      skipped positionally at all -- which is exactly how this file stopped compiling when      !
      !      #95 added psi_tlp/stress_arrestor to leaf_photo_params_t (see to_params below). ----------!
      env = leaf_env_t(par=env_c%par, leaf_temp=env_c%leaf_temp, vpd=env_c%vpd, ca=env_c%ca,      &
                       pressure=env_c%pressure, psi_leaf=env_c%psi_leaf, gb=env_c%gb,             &
                       psi=env_c%psi)
   end function to_env

   pure function to_params(p_c) result(p)
      type(leaf_params_c), intent(in) :: p_c
      type(leaf_photo_params_t)       :: p
      !----- KEYWORD form (see to_env). This constructor was POSITIONAL and broke the moment           !
      !      leaf_photo_params_t grew components in the middle: #95 inserted psi_tlp and                 !
      !      stress_arrestor after sref_stomata, so the trailing values shifted by three and ifx         !
      !      rejected the file with "Omitted component is not initialized" for o2_mol_frac /             !
      !      absorptance / phi_psii. The pylib is NOT part of the default build or of ctest, so nothing  !
      !      caught it -- see issue #100.                                                                !
      !                                                                                          !
      !      psi_tlp, stress_arrestor and wstress_nonstomatal are absent from leaf_params_c and take     !
      !      their type defaults (-2.0 MPa, ARREST_GS_CLAMP, .false.). Exposing them across the C ABI    !
      !      is a deliberate follow-up, not an oversight: adding a field to a bind(c) struct is an ABI   !
      !      break for any existing caller. --------------------------------------------------------!
      p = leaf_photo_params_t(pathway=p_c%pathway, vcmax25=p_c%vcmax25, jmax25=p_c%jmax25,        &
             tpu25=p_c%tpu25, rd25=p_c%rd25, kp25=p_c%kp25, g0=p_c%g0, g1=p_c%g1, d0=p_c%d0,      &
             quantum_yield=p_c%quantum_yield, theta_j=p_c%theta_j, theta_cj=p_c%theta_cj,         &
             theta_ic=p_c%theta_ic, lambda25=p_c%lambda25, psi_open=p_c%psi_open,                 &
             psi_close=p_c%psi_close, lambda_psi_exp=p_c%lambda_psi_exp,                          &
             sref_stomata=p_c%sref_stomata,                                                       &
             kc25=p_c%kc25, ko25=p_c%ko25, gstar25=p_c%gstar25, ea_kc=p_c%ea_kc,                  &
             ea_ko=p_c%ea_ko, ea_gstar=p_c%ea_gstar, ea_vcmax=p_c%ea_vcmax,                       &
             ea_jmax=p_c%ea_jmax, ea_rd=p_c%ea_rd, hd_vcmax=p_c%hd_vcmax, hd_jmax=p_c%hd_jmax,    &
             hd_rd=p_c%hd_rd, ds_vcmax=p_c%ds_vcmax, ds_jmax=p_c%ds_jmax, ds_rd=p_c%ds_rd,        &
             o2_mol_frac=p_c%o2_mol_frac, absorptance=p_c%absorptance, phi_psii=p_c%phi_psii)
   end function to_params

end module meds_capi_leaf
