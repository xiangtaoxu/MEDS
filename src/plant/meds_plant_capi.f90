!==========================================================================================!
! meds_plant_capi -- a thin ISO_C_BINDING shim that exposes the plant-ecophysiology model to  !
! C / Python (ctypes): leaf gas exchange AND leaf phenology. It is NOT part of libmeds_plant;  !
! it is compiled only into the optional shared library libmeds_plant_c (CMake                  !
! -DMEDS_BUILD_PYLIB=ON), so the core model stays foreign-call-free.                            !
!                                                                                          !
! LEAF gas exchange -- C-interoperable mirrors of leaf_env_t / leaf_photo_params_t / leaf_flux_t !
! (same field ORDER; real64 -> c_double, int32 -> c_int, the flux `converged` logical -> 0/1):  !
!   * meds_leaf_solve         -- the coupled A-gs-Ci solver (solve_leaf_gas_exchange).           !
!   * meds_assimilation_demand_c3 -- the raw C3 FvCB demand kernel at a prescribed Ci: gross A    !
!                                and the Ac/Aj/Ap limitation rates, NO temperature scaling.       !
!   * meds_electron_transport_j -- J from Jmax and light (the non-rectangular hyperbola).          !
!   * meds_peaked_arrhenius / meds_arrhenius -- the temperature-response functions.                !
! Composing the kinetics + J + demand draws an A-Ci curve from Vcmax/Jmax DIRECTLY.               !
!                                                                                          !
! LEAF PHENOLOGY -- C mirrors of pheno_env_t / pheno_params_t / pheno_state_t / pheno_out_t      !
! (the two logicals hemis_north / water_use_potential -> c_int 0/1):                             !
!   * meds_phenology_step     -- advance phenology_kernel ONE daily step (state advanced in place).!
! The Python side (python/meds/leaf/_ffi.py + meds/leaf/pheno.py) mirrors all of these structs.  !
!==========================================================================================!
module meds_plant_capi
   use iso_c_binding,            only : c_double, c_int
   use meds_kinds,              only : wp, ik
   use meds_plant_types,         only : leaf_env_t, leaf_photo_params_t, leaf_flux_t
   use meds_plant_types,         only : pheno_env_t, pheno_params_t, pheno_state_t, pheno_out_t
   use meds_leaf_gas_exchange,        only : solve_leaf_gas_exchange
   use meds_leaf_gas_exchange,only : assimilation_demand_c3, electron_transport_j
   use meds_temp_response, only : peaked_arrhenius_scale, arrhenius_scale
   use meds_phenology,     only : phenology_kernel
   implicit none
   private

   public :: leaf_env_c, leaf_params_c, leaf_flux_c, leaf_c3_demand_c
   public :: meds_leaf_solve, meds_assimilation_demand_c3, meds_electron_transport_j
   public :: meds_peaked_arrhenius, meds_arrhenius
   public :: pheno_env_c, pheno_params_c, pheno_state_c, pheno_out_c, meds_phenology_step

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

   !----- PHENOLOGY: C mirror of pheno_env_t (6 doubles + 2 ints; hemis_north 0/1). ---------!
   type, bind(c) :: pheno_env_c
      real(c_double) :: temp_day, soil_temp, avail_water, dmax_leaf_psi, rad, daylength
      integer(c_int) :: doy, hemis_north
   end type pheno_env_c

   !----- C mirror of pheno_params_t (masks + rate scales + all cue params; water_use_pot 0/1). !
   type, bind(c) :: pheno_params_c
      integer(c_int) :: flush_cue_mask, shed_cue_mask
      real(c_double) :: cue_sharpness, k_flush_max, k_shed_max, tau_flush, tau_shed
      real(c_double) :: gdd_base_temp, chill_base_temp, phen_a, phen_b, phen_c
      real(c_double) :: cold_drop_daylength, cold_drop_soiltemp1, cold_drop_soiltemp2
      integer(c_int) :: water_use_potential
      real(c_double) :: water_off_threshold, water_on_threshold, water_window, water_width
      real(c_double) :: leaf_psi_tlp, low_psi_threshold, high_psi_threshold
      real(c_double) :: photo_crit, photo_slope
      real(c_double) :: light_on_threshold, light_width, light_window
   end type pheno_params_c

   !----- C mirror of pheno_state_t (8 doubles; the prognostic memory, in/out). -------------!
   type, bind(c) :: pheno_state_c
      real(c_double) :: flush_drive, shed_drive, gdd, chill, water_avg, low_psi_days, &
                        high_psi_days, light_avg
   end type pheno_state_c

   !----- C mirror of pheno_out_t (2 doubles + 1 int). --------------------------------------!
   type, bind(c) :: pheno_out_c
      real(c_double) :: leaf_flush_rate, leaf_shed_rate
      integer(c_int) :: cue_limiting
   end type pheno_out_c

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
   ! Advance one cohort's leaf phenology ONE step: unpack the C structs, call phenology_kernel !
   ! (state advanced in place), pack the advanced state + the two relative rates back out.     !
   ! `state_c` is intent(inout): the caller keeps it across days (the phenological memory).    !
   !---------------------------------------------------------------------------------------!
   subroutine meds_phenology_step(env_c, p_c, dt, state_c, out_c) bind(c, name="meds_phenology_step")
      type(pheno_env_c),    intent(in)    :: env_c
      type(pheno_params_c), intent(in)    :: p_c
      real(c_double), value, intent(in)   :: dt
      type(pheno_state_c),  intent(inout) :: state_c
      type(pheno_out_c),    intent(out)   :: out_c
      type(pheno_env_t)    :: env
      type(pheno_params_t) :: p
      type(pheno_state_t)  :: state
      type(pheno_out_t)    :: out

      !----- Unpack the environment. -----------------------------------------------------!
      env%temp_day      = env_c%temp_day
      env%soil_temp     = env_c%soil_temp
      env%avail_water   = env_c%avail_water
      env%dmax_leaf_psi = env_c%dmax_leaf_psi
      env%rad           = env_c%rad
      env%daylength     = env_c%daylength
      env%doy           = int(env_c%doy, ik)
      env%hemis_north   = env_c%hemis_north /= 0_c_int

      !----- Unpack the parameters. ------------------------------------------------------!
      p%flush_cue_mask      = int(p_c%flush_cue_mask, ik)
      p%shed_cue_mask       = int(p_c%shed_cue_mask, ik)
      p%cue_sharpness       = p_c%cue_sharpness
      p%k_flush_max         = p_c%k_flush_max
      p%k_shed_max          = p_c%k_shed_max
      p%tau_flush           = p_c%tau_flush
      p%tau_shed            = p_c%tau_shed
      p%gdd_base_temp       = p_c%gdd_base_temp
      p%chill_base_temp     = p_c%chill_base_temp
      p%phen_a              = p_c%phen_a
      p%phen_b              = p_c%phen_b
      p%phen_c              = p_c%phen_c
      p%cold_drop_daylength = p_c%cold_drop_daylength
      p%cold_drop_soiltemp1 = p_c%cold_drop_soiltemp1
      p%cold_drop_soiltemp2 = p_c%cold_drop_soiltemp2
      p%water_use_potential = p_c%water_use_potential /= 0_c_int
      p%water_off_threshold = p_c%water_off_threshold
      p%water_on_threshold  = p_c%water_on_threshold
      p%water_window        = p_c%water_window
      p%water_width         = p_c%water_width
      p%leaf_psi_tlp        = p_c%leaf_psi_tlp
      p%low_psi_threshold   = p_c%low_psi_threshold
      p%high_psi_threshold  = p_c%high_psi_threshold
      p%photo_crit          = p_c%photo_crit
      p%photo_slope         = p_c%photo_slope
      p%light_on_threshold  = p_c%light_on_threshold
      p%light_width         = p_c%light_width
      p%light_window        = p_c%light_window

      !----- Unpack the prognostic state (the caller keeps it across days). ---------------!
      state%flush_drive   = state_c%flush_drive
      state%shed_drive    = state_c%shed_drive
      state%gdd           = state_c%gdd
      state%chill         = state_c%chill
      state%water_avg     = state_c%water_avg
      state%low_psi_days  = state_c%low_psi_days
      state%high_psi_days = state_c%high_psi_days
      state%light_avg     = state_c%light_avg

      call phenology_kernel(env, p, real(dt, wp), state, out)

      !----- Pack the advanced state + outputs back. -------------------------------------!
      state_c%flush_drive   = state%flush_drive
      state_c%shed_drive    = state%shed_drive
      state_c%gdd           = state%gdd
      state_c%chill         = state%chill
      state_c%water_avg     = state%water_avg
      state_c%low_psi_days  = state%low_psi_days
      state_c%high_psi_days = state%high_psi_days
      state_c%light_avg     = state%light_avg
      out_c%leaf_flush_rate = out%leaf_flush_rate
      out_c%leaf_shed_rate  = out%leaf_shed_rate
      out_c%cue_limiting    = int(out%cue_limiting, c_int)
   end subroutine meds_phenology_step

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

end module meds_plant_capi
