!==========================================================================================!
! meds_column_co2 -- THE stateless column CO2-balance kernels (design MEDS_COLUMN_CO2_BALANCE_ !
! DESIGN.md, section 4). The canopy-air-space CO2 mixing ratio `can_co2 [umol/mol]` is the      !
! THIRD prognostic CAS twin beside can_enthalpy / can_shv: it is advanced each fast substep from  !
! the net biotic source (respiration - GPP) and the turbulent exchange with the free atmosphere,   !
! closing a machine-precision CO2 budget (resid ~ 0), exactly mirroring meds_column_energy's         !
! canopy_air_update -- with one unit change: CO2 is a MOLAR mixing ratio, so the CAS capacity is     !
! the dry-air MOLAR column ccapcan = can_dmol*can_depth [mol/m2] (vs the mass wcapcan [kg/m2]).       !
!                                                                                          !
!   * canopy_air_co2_update        -- the CAS storage + atmospheric flux (implicit-atm; closed resid). !
!   * aggregate_cohort_co2_fluxes  -- reduce per-cohort GPP/respiration to [umol CO2 / m2 ground / s].   !
!   * heterotrophic_respiration_flux -- the MVP soil Rh: Q10 or ED2 capped-exp x moisture, on a frozen    !
!                                       (prescribed at P0) soil-carbon pool. (DAMM lands at P1.)            !
!   * column_co2_step              -- the host-side assembler: aggregate -> Rh -> advance twin -> guard.    !
!                                                                                          !
! All compute kernels are `pure` and take bare scalars/arrays (the shipped canopy_air_update idiom),     !
! so they are GPU-eligible and link src/shared ONLY. can_co2 rides in cas_state_t (biophysics); the       !
! driver passes it here as a bare scalar -- no biogeochem->biophysics library edge.                        !
!==========================================================================================!
module meds_column_co2
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : mmdry, tiny_num, kgCday_2_umols, r_gas_kj, o2_air_frac, damm_flux_factor
   use meds_biogeochem_types, only : co2_opts_t, column_co2_budget_t, cohort_co2_flux_t,       &
                                     damm_params_t, HR_Q10, HR_EXP_ED2, HR_DAMM
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   private

   public :: canopy_air_co2_update, aggregate_cohort_co2_fluxes
   public :: heterotrophic_respiration_flux, heterotrophic_respiration_damm, column_co2_step

contains

   !---------------------------------------------------------------------------------------!
   ! Canopy-air-space CO2 update -- the THIRD prognostic twin. Advances the dry-air CO2      !
   ! mixing ratio one step from the net biotic source (respiration - GPP) and the atmospheric !
   ! exchange (IMPLICIT in the atm term, like canopy_air_update's enthalpy branch). Molar      !
   ! capacity because can_co2 is a mole fraction. Fills the closed CO2 budget (resid ~ 0).      !
   !---------------------------------------------------------------------------------------!
   pure subroutine canopy_air_co2_update(can_co2, can_depth, can_shv,                        &
                                         gross_primary_prod, plant_respiration,              &
                                         heterotrophic_respiration,                          &
                                         ustar, temp2, co2_atm, rho_air, dt, budget)
      real(wp), intent(inout) :: can_co2                    !< [umol/mol]  prognostic third twin
      real(wp), intent(in)    :: can_depth                  !< [m]         CAS depth (fixed within a step at MVP)
      real(wp), intent(in)    :: can_shv                    !< [kg/kg]     CAS humidity (-> dry-air molar density)
      real(wp), intent(in)    :: gross_primary_prod         !< [umol/m2/s] patch GPP (uptake, >= 0)
      real(wp), intent(in)    :: plant_respiration          !< [umol/m2/s] leaf+stem+root+growth+storage (source)
      real(wp), intent(in)    :: heterotrophic_respiration  !< [umol/m2/s] soil Rh (source, >= 0)
      real(wp), intent(in)    :: ustar                      !< [m/s]       friction velocity (shared twin conductance)
      real(wp), intent(in)    :: temp2                      !< [--]        M-O scalar-transfer coeff c3 (aero%temp2)
      real(wp), intent(in)    :: co2_atm                    !< [umol/mol]  free-atmosphere reference
      real(wp), intent(in)    :: rho_air                    !< [kg/m3]     moist CAS air density
      real(wp), intent(in)    :: dt                         !< [s]
      type(column_co2_budget_t), intent(out) :: budget      !< nee/nep/loss2atm/storage/resid
      real(wp) :: can_dmol, ccapcan, cci, f_bio, gatm_co2, co2_new

      can_dmol = rho_air * (1.0_wp - can_shv) / mmdry        ! [mol_dryair/m3] DRY-air molar density (exact)
      ccapcan  = can_dmol * can_depth                        ! [mol_air/m2]    MOLAR CAS capacity (cf wcapcan)
      cci      = 1.0_wp / max(ccapcan, tiny_num)             ! [m2/mol]
      f_bio    = heterotrophic_respiration + plant_respiration - gross_primary_prod  ! [umol/m2/s] Reco - GPP
      !----- atm<->CAS MOLAR conductance = can_dmol*ustar*c3, with c3 (temp2) the M-O scalar-        !
      !      transfer coefficient (aero_out%temp2). Sharing ustar AND temp2 with the energy/vapour   !
      !      twins keeps all three CAS twins on one turbulence basis; dropping it (c3=1) over-couples !
      !      the CAS to the free atmosphere and damps nocturnal sub-canopy CO2 build-up (BUG8). ------!
      gatm_co2 = can_dmol * ustar * temp2                    ! [mol_air/m2/s]  atm<->CAS molar exchange
      !----- Implicit in the atmospheric-exchange term (L-stable); biotic source explicit. --------!
      co2_new  = (can_co2 + dt*cci*(f_bio + gatm_co2*co2_atm)) / (1.0_wp + dt*cci*gatm_co2)
      !----- Diagnostics + closed budget (resid = 0 by substitution). -----------------------------!
      budget%resid    = ccapcan*(co2_new - can_co2) - dt*(f_bio + gatm_co2*(co2_atm - co2_new))  ! = 0
      budget%nee      = f_bio
      budget%nep      = -f_bio
      budget%loss2atm = gatm_co2*(co2_new - co2_atm)
      budget%storage  = ccapcan*co2_new
      can_co2  = co2_new
   end subroutine canopy_air_co2_update

   !---------------------------------------------------------------------------------------!
   ! Reduce the patch's per-cohort carbon fluxes to patch-ground [umol CO2 / m2 ground / s].  !
   ! Leaf GPP / respiration are per m2 LEAF (x leaf_area x nplant = LAI weighting); stem / root !
   ! respiration are per PLANT (x nplant). Explicit do-loop (not sum(...)) for clarity and to     !
   ! sidestep the nvfortran array-temporary trap (CLAUDE.md issue #7). growth/storage resp are     !
   ! set by the caller (committed daily amounts amortized), not here; both are 0 at MVP.            !
   !---------------------------------------------------------------------------------------!
   pure subroutine aggregate_cohort_co2_fluxes(n, nplant, leaf_area, a_gross, rd,             &
                                               stem_resp, root_resp, out)
      integer(ik), intent(in) :: n
      real(wp), intent(in)  :: nplant(n), leaf_area(n)        !< [plant/m2], [m2 leaf/plant]
      real(wp), intent(in)  :: a_gross(n), rd(n)              !< [umol CO2/m2 leaf/s]
      real(wp), intent(in)  :: stem_resp(n), root_resp(n)     !< [umol CO2/plant/s]
      type(cohort_co2_flux_t), intent(out) :: out
      integer(ik) :: i

      out%gross_primary_prod = 0.0_wp ; out%leaf_respiration = 0.0_wp
      out%stem_respiration   = 0.0_wp ; out%root_respiration = 0.0_wp
      out%growth_respiration = 0.0_wp ; out%storage_respiration = 0.0_wp
      do i = 1_ik, n
         out%gross_primary_prod = out%gross_primary_prod + a_gross(i)  * leaf_area(i) * nplant(i)
         out%leaf_respiration   = out%leaf_respiration   + rd(i)       * leaf_area(i) * nplant(i)
         out%stem_respiration   = out%stem_respiration   + stem_resp(i)              * nplant(i)
         out%root_respiration   = out%root_respiration   + root_resp(i)              * nplant(i)
      end do
   end subroutine aggregate_cohort_co2_fluxes

   !---------------------------------------------------------------------------------------!
   ! Heterotrophic respiration [umol/m2/s] -- a pluggable soil-gas-flux dispatcher over        !
   ! opts%hr_model. HR_DAMM computes the ABSOLUTE mechanistic rate (own moisture physics; it     !
   ! does NOT use rh_k_base / theta_dry / the empirical f_water). HR_Q10 / HR_EXP_ED2 share the    !
   ! empirical `frozen pool x rate x f_temp(T) x f_water(theta)` form. The soil-carbon pool is     !
   ! read-only in the fast loop (prescribed at P0); this is the fast coupling flux to the CAS twin. !
   !---------------------------------------------------------------------------------------!
   pure function heterotrophic_respiration_flux(fast_soil_carbon, soil_temp, theta,          &
                                                theta_dry, theta_sat, opts) result(rh)
      real(wp), intent(in) :: fast_soil_carbon      !< [kgC/m2]  decomposable pool (frozen in the fast loop)
      real(wp), intent(in) :: soil_temp             !< [K]       representative (root-weighted) soil temperature
      real(wp), intent(in) :: theta, theta_dry, theta_sat   !< [m3/m3] moisture, air-dry floor, porosity
      type(co2_opts_t), intent(in) :: opts
      real(wp) :: rh                                !< [umol/m2/s]
      real(wp) :: f_temp, f_water, rel

      select case (opts%hr_model)
      case (HR_DAMM)                                ! mechanistic (theta_dry / rh_k_base UNUSED)
         rh = heterotrophic_respiration_damm(fast_soil_carbon, soil_temp, theta, theta_sat, opts%damm)
      case (HR_Q10, HR_EXP_ED2)                     ! empirical: frozen pool x f_temp(T) x f_water(theta)
         rel     = min(1.0_wp, max(0.0_wp, (theta - theta_dry) / max(theta_sat - theta_dry, tiny_num)))
         f_water = water_modifier(rel, opts)
         if (opts%hr_model == HR_EXP_ED2) then      ! ED2 scheme-0: min(1, exp(a*(T - T_sat)))
            f_temp = min(1.0_wp, exp(opts%resp_temp_increase * (soil_temp - opts%resp_temp_ref)))
         else                                       ! HR_Q10 (Collatz/K13): q10**((T - T_ref)/10)
            f_temp = opts%rh_q10 ** ((soil_temp - opts%rh_t_ref) * 0.1_wp)
         end if
         rh = fast_soil_carbon * opts%rh_k_base * f_temp * f_water * kgCday_2_umols   ! kgC/m2/day -> umol/m2/s
      case default                                  ! unrecognized selector: fail loud (pure -> NaN, not a
         rh = ieee_value(rh, ieee_quiet_nan)        ! silent Q10 fallback; caught by the has_nan/budget guard)
      end select
   end function heterotrophic_respiration_flux

   !---------------------------------------------------------------------------------------!
   ! DAMM heterotrophic respiration [umol CO2 m-2 s-1] (Davidson et al. 2012 GCB 18:371-384). !
   ! An Arrhenius maximum velocity x TWO Michaelis-Menten factors: soluble-C substrate (rises    !
   ! with moisture via a liquid-film theta^3 diffusion term) and O2 (falls with moisture via the   !
   ! air-filled-porosity (theta_sat - theta)^(4/3) gas term). The unimodal moisture response is     !
   ! EMERGENT -- no empirical f_water. Needs only soil_temp[K], theta, theta_sat, and the frozen     !
   ! soil-C pool. The respiring depth appears twice (SOC->conc AND flux depth-integral) and does     !
   ! NOT cancel (Sx enters the nonlinear MM term). Pure/bare-scalar/GPU-eligible.                     !
   !---------------------------------------------------------------------------------------!
   pure function heterotrophic_respiration_damm(fast_soil_carbon, soil_temp, theta,          &
                                                theta_sat, damm) result(rh)
      real(wp), intent(in) :: fast_soil_carbon      !< [kgC/m2]  frozen decomposable pool (= DAMM total soil C)
      real(wp), intent(in) :: soil_temp             !< [K]       soil temperature (NATIVE Kelvin -> Arrhenius)
      real(wp), intent(in) :: theta                 !< [m3/m3]   volumetric soil moisture (native; no conversion)
      real(wp), intent(in) :: theta_sat             !< [m3/m3]   porosity = air-filled-porosity ceiling
      type(damm_params_t), intent(in) :: damm
      real(wp) :: rh                                !< [umol/m2/s]
      real(wp) :: vmax, sx_total, sx, a_air, o2, mm_sx, mm_o2, r_sx

      !----- Arrhenius max velocity (Ea & R BOTH in kJ/mol). --------------------------------------!
      vmax     = damm%alpha_sx * exp( -damm%ea_sx / (r_gas_kj * soil_temp) )        ! [mgC cm-3 h-1]
      !----- Soluble-C substrate: column pool -> volumetric conc, soluble fraction, liquid diff (^3).!
      sx_total = fast_soil_carbon * 0.1_wp / damm%depth_cm                          ! [gC cm-3] (kgC/m2 over depth)
      sx       = damm%p_soluble * sx_total * damm%d_liq * theta**3                  ! [gC cm-3]
      !----- Oxygen: air-filled porosity CLAMPED >= 0 before the 4/3 power (else NaN when theta>sat).!
      a_air    = max(theta_sat - theta, 0.0_wp)                                     ! [m3/m3]
      o2       = damm%d_gas * o2_air_frac * a_air ** (4.0_wp / 3.0_wp)              ! [cm3 O2 cm-3 air]
      !----- Dual Michaelis-Menten (both self-bounded in [0,1); kM > 0 => safe at conc = 0). ------!
      mm_sx    = sx / (damm%km_sx + sx)
      mm_o2    = o2 / (damm%km_o2 + o2)
      r_sx     = vmax * mm_sx * mm_o2                                               ! [mgC cm-3 h-1]
      !----- Depth-integrate + convert to MEDS units (231.269 folds cm, h, mgC -> umol). ----------!
      rh       = r_sx * damm%depth_cm * damm_flux_factor                           ! [umol CO2 m-2 s-1]
   end function heterotrophic_respiration_damm

   !----- One-sided exponential moisture modifier about resp_opt_water (shared by both empirical  !
   !      Rh cases). rel = relative saturation in [0,1]; f_water in (0,1]. --------------------------!
   pure function water_modifier(rel, opts) result(f_water)
      real(wp),         intent(in) :: rel
      type(co2_opts_t), intent(in) :: opts
      real(wp) :: f_water
      if (rel <= opts%resp_opt_water) then
         f_water = exp((rel - opts%resp_opt_water) * opts%resp_water_below_opt)
      else
         f_water = exp((opts%resp_opt_water - rel) * opts%resp_water_above_opt)
      end if
   end function water_modifier

   !---------------------------------------------------------------------------------------!
   ! The host-side per-patch column CO2 step: aggregate cohorts, attach committed growth /    !
   ! storage respiration, compute Rh, sum autotrophic respiration, advance the twin, and       !
   ! enforce the closed-budget guard in Debug (the uniform biophysics discipline). NOT pure     !
   ! (it error-stops on a numerical fault; the algebra guarantees the physics).                  !
   !---------------------------------------------------------------------------------------!
   subroutine column_co2_step(cas_can_co2, can_depth, can_shv, ustar, temp2, co2_atm, rho_air, dt, &
                              n, nplant, leaf_area, a_gross, rd, stem_resp, root_resp,         &
                              growth_resp_committed, storage_resp_committed,                   &
                              fast_soil_carbon, soil_temp, theta, theta_dry, theta_sat,        &
                              opts, budget)
      real(wp),    intent(inout) :: cas_can_co2            !< [umol/mol] = cas%can_co2 (passed by reference)
      real(wp),    intent(in)    :: can_depth, can_shv, ustar, temp2, co2_atm, rho_air, dt
      integer(ik), intent(in)    :: n
      real(wp),    intent(in)    :: nplant(n), leaf_area(n), a_gross(n), rd(n), stem_resp(n), root_resp(n)
      real(wp),    intent(in)    :: growth_resp_committed, storage_resp_committed  ! [umol/m2/s] MVP = 0
      real(wp),    intent(in)    :: fast_soil_carbon, soil_temp, theta, theta_dry, theta_sat
      type(co2_opts_t),          intent(in)  :: opts
      type(column_co2_budget_t), intent(out) :: budget
      type(cohort_co2_flux_t) :: coh
      real(wp) :: hetero, plant_resp, scale

      call aggregate_cohort_co2_fluxes(n, nplant, leaf_area, a_gross, rd, stem_resp, root_resp, coh)
      coh%growth_respiration  = growth_resp_committed
      coh%storage_respiration = storage_resp_committed
      hetero     = heterotrophic_respiration_flux(fast_soil_carbon, soil_temp, theta,          &
                                                  theta_dry, theta_sat, opts)
      plant_resp = coh%leaf_respiration + coh%stem_respiration + coh%root_respiration           &
                 + coh%growth_respiration + coh%storage_respiration
      call canopy_air_co2_update(cas_can_co2, can_depth, can_shv, coh%gross_primary_prod,       &
                                 plant_resp, hetero, ustar, temp2, co2_atm, rho_air, dt, budget)
      !----- Closed-budget guard (mixed rtol/atol form; cf. soil_energy_flux). --------------------!
      scale = max(abs(budget%storage), 1.0_wp)
      if (opts%debug_error .and. abs(budget%resid) > opts%rtol * scale + opts%atol) then
         error stop 'column_co2_step: CO2 budget did not close'
      end if
   end subroutine column_co2_step

end module meds_column_co2
