!==========================================================================================!
! meds_leaf_demo -- a standalone driver that exercises the leaf-physiology module and writes  !
! canonical leaf gas-exchange response curves as CSV (for post_proc/plot_leaf_response.py).  !
!                                                                                          !
! It reads the SAME TOML configuration as meds_main (so the per-PFT photosynthesis traits and  !
! the model selectors are the live ones), picks one PFT, and sweeps four diagnostics for a      !
! reference leaf:                                                                              !
!   * <prefix>_aci.csv   -- A-Ci demand curve at fixed light/temperature (Ac/Aj/Ap envelope),  !
!                           computed from the demand functions directly (stomata bypassed).     !
!   * <prefix>_apar.csv  -- light response: the fully coupled solve over PAR.                   !
!   * <prefix>_atemp.csv -- temperature response: the coupled solve over leaf temperature.      !
!   * <prefix>_gsvpd.csv -- stomatal closure: gs(VPD) for all three stomatal models.            !
!                                                                                          !
! Usage:  meds_leaf_demo [main.toml] [pft_index] [out_prefix]   (defaults: meds_config_main.toml 1 leafdemo)!
!==========================================================================================!
program meds_leaf_demo
   use meds_kinds,              only : wp, ik
   use meds_constants,          only : t_kelvin, p_std
   use meds_config,             only : meds_config_t, SM_LEUNING, SM_MEDLYN, SM_KATUL
   use meds_pft_params,         only : PATH_C4
   use meds_config_io,          only : load_meds_config
   use meds_leaf_temp_response, only : temp_response, arrhenius_scale
   use meds_leaf_photosynthesis,only : assim_demand_c3, assim_demand_c4, electron_transport_j
   use meds_leaf_physiology,    only : leaf_env_t, leaf_flux_t, leaf_gas_exchange
   implicit none

   type(meds_config_t) :: cfg
   type(leaf_env_t)    :: env
   type(leaf_flux_t)   :: flux, fl, fm, fk
   character(len=256)  :: path, prefix, arg
   integer(ik)         :: ipft, i, u
   real(wp)            :: ci, ag, ac, aj, ap, par, tleaf, vpd

   !----- Command line: config path, PFT index, output prefix. -------------------------!
   path = 'meds_config_main.toml' ; ipft = 1_ik ; prefix = 'leafdemo'
   if (command_argument_count() >= 1_ik) call get_command_argument(1, path)
   if (command_argument_count() >= 2_ik) then
      call get_command_argument(2, arg) ; read(arg,*) ipft
   end if
   if (command_argument_count() >= 3_ik) call get_command_argument(3, prefix)

   call load_meds_config(trim(path), cfg)
   if (ipft < 1_ik .or. ipft > cfg%pft%n) error stop 'meds_leaf_demo: PFT index out of range'
   write(*,'(2a)')      ' config : ', trim(path)
   write(*,'(a,i0,a,i0)') ' PFT    : ', ipft, ' of ', cfg%pft%n
   write(*,'(a,i0)')    ' stomatal_model (gs-VPD sweeps all three) : ', cfg%stomatal_model

   !----- A reference midday tropical leaf. ---------------------------------------------!
   env%par = 1500.0_wp ; env%leaf_temp = t_kelvin + 25.0_wp ; env%vpd = 1500.0_wp
   env%ca = 400.0_wp ; env%pressure = p_std ; env%psi_leaf = 0.0_wp ; env%gb = 0.0_wp

   !----- 1. A-Ci demand curve (stomata bypassed): vary Ci, report the three limits. The     !
   !         pathway comment lets the plotter label the third limit correctly (C3 TPU vs C4   !
   !         PEPcase). -------------------------------------------------------------------!
   open(newunit=u, file=trim(prefix)//'_aci.csv', status='replace', action='write')
   if (cfg%pft%photosynthetic_pathway(ipft) == PATH_C4) then
      write(u,'(a)') '# pathway = c4'
   else
      write(u,'(a)') '# pathway = c3'
   end if
   write(u,'(a)') 'ci,a_net,ac,aj,ap'
   do i = 0_ik, 60_ik
      ci = 20.0_wp + 10.0_wp * real(i, wp)             ! 20 .. 620 ppm
      call demand_at_ci(cfg, ipft, env, ci, ag, ac, aj, ap)
      write(u,'(f0.3,4(",",es15.8))') ci, ag - leaf_rd(cfg, ipft, env%leaf_temp), ac, aj, ap
   end do
   close(u)

   !----- 2. A-PAR light response (fully coupled solve). ---------------------------------!
   open(newunit=u, file=trim(prefix)//'_apar.csv', status='replace', action='write')
   write(u,'(a)') 'par,a_net,gs,ci'
   do i = 0_ik, 40_ik
      env%par = 50.0_wp * real(i, wp)                  ! 0 .. 2000 umol/m2/s
      call leaf_gas_exchange(env, cfg, ipft, flux)
      write(u,'(f0.2,3(",",es15.8))') env%par, flux%a_net, flux%gs, flux%ci
   end do
   env%par = 1500.0_wp

   !----- 3. A-leaf-temperature response (fully coupled solve). --------------------------!
   open(newunit=u, file=trim(prefix)//'_atemp.csv', status='replace', action='write')
   write(u,'(a)') 'tleaf_c,a_net,gs,ci'
   do i = 0_ik, 40_ik
      tleaf = 5.0_wp + real(i, wp)                     ! 5 .. 45 degC
      env%leaf_temp = t_kelvin + tleaf
      call leaf_gas_exchange(env, cfg, ipft, flux)
      write(u,'(f0.2,3(",",es15.8))') tleaf, flux%a_net, flux%gs, flux%ci
   end do
   env%leaf_temp = t_kelvin + 25.0_wp

   !----- 4. gs(VPD) for all three stomatal models. --------------------------------------!
   open(newunit=u, file=trim(prefix)//'_gsvpd.csv', status='replace', action='write')
   write(u,'(a)') 'vpd_pa,gs_leuning,gs_medlyn,gs_katul'
   do i = 1_ik, 40_ik
      vpd = 100.0_wp * real(i, wp)                     ! 100 .. 4000 Pa
      env%vpd = vpd
      cfg%stomatal_model = SM_LEUNING ; call leaf_gas_exchange(env, cfg, ipft, fl)
      cfg%stomatal_model = SM_MEDLYN  ; call leaf_gas_exchange(env, cfg, ipft, fm)
      cfg%stomatal_model = SM_KATUL   ; call leaf_gas_exchange(env, cfg, ipft, fk)
      write(u,'(f0.2,3(",",es15.8))') vpd, fl%gs, fm%gs, fk%gs
   end do
   close(u)

   write(*,'(4a)') ' wrote  : ', trim(prefix), '_{aci,apar,atemp,gsvpd}.csv'

contains

   !----- Leaf respiration at temperature T (mirrors the solver's T-scaling). ------------!
   function leaf_rd(cfg, ipft, t_leaf) result(rd)
      type(meds_config_t), intent(in) :: cfg
      integer(ik),         intent(in) :: ipft
      real(wp),            intent(in) :: t_leaf
      real(wp)                        :: rd
      rd = temp_response(cfg%temp_response_form, cfg%pft%rd25(ipft), cfg%ea_rd, cfg%hd_rd,     &
                         cfg%ds_rd, t_leaf)
   end function leaf_rd

   !----- Gross + raw limitation rates at a prescribed Ci (mirrors the solver's demand).    !
   !       Uses UNSTRESSED capacity (no beta downregulation); the demo's reference leaf is    !
   !       well-watered (psi_leaf = 0), so this matches the coupled A-PAR/A-T panels.  -------!
   subroutine demand_at_ci(cfg, ipft, env, ci, ag, ac, aj, ap)
      type(meds_config_t), intent(in)  :: cfg
      integer(ik),         intent(in)  :: ipft
      type(leaf_env_t),    intent(in)  :: env
      real(wp),            intent(in)  :: ci
      real(wp),            intent(out) :: ag, ac, aj, ap
      real(wp) :: vcmax, jmax, tpu, kc_ppm, ko_ppm, gstar_ppm, o2_ppm, jrate, aj_light, kp_eff
      associate (t => cfg%pft, tl => env%leaf_temp, pr => env%pressure)
         vcmax = temp_response(cfg%temp_response_form, t%vcmax25(ipft), cfg%ea_vcmax,          &
                               cfg%hd_vcmax, cfg%ds_vcmax, tl)
         o2_ppm = cfg%o2_mol_frac * 1.0e6_wp
         if (t%photosynthetic_pathway(ipft) == PATH_C4) then
            aj_light = t%quantum_yield_c4(ipft) * cfg%leaf_absorptance * env%par
            kp_eff   = t%kp25(ipft) * pr / p_std
            call assim_demand_c4(ci, vcmax, aj_light, kp_eff, cfg%colimitation,                &
                                 t%theta_cj_c4(ipft), t%theta_ic_c4(ipft), ag, ac, aj, ap)
         else
            jmax = temp_response(cfg%temp_response_form, t%jmax25(ipft), cfg%ea_jmax,          &
                                 cfg%hd_jmax, cfg%ds_jmax, tl)
            tpu  = temp_response(cfg%temp_response_form, t%tpu25(ipft), cfg%ea_vcmax,          &
                                 cfg%hd_vcmax, cfg%ds_vcmax, tl)
            kc_ppm    = arrhenius_scale(cfg%kc25,    cfg%ea_kc,    tl) / pr * 1.0e6_wp
            ko_ppm    = arrhenius_scale(cfg%ko25,    cfg%ea_ko,    tl) / pr * 1.0e6_wp
            gstar_ppm = arrhenius_scale(cfg%gstar25, cfg%ea_gstar, tl) / pr * 1.0e6_wp
            jrate     = electron_transport_j(env%par, cfg%leaf_absorptance, cfg%phi_psii,      &
                                             jmax, t%theta_j(ipft))
            call assim_demand_c3(ci, vcmax, jrate, tpu, gstar_ppm, kc_ppm, ko_ppm, o2_ppm,     &
                                 cfg%colimitation, t%theta_j(ipft), ag, ac, aj, ap)
         end if
      end associate
   end subroutine demand_at_ci

end program meds_leaf_demo
