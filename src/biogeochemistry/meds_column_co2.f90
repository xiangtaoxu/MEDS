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
   use meds_constants,        only : mmdry, tiny_num, kgCday_2_umols
   use meds_biogeochem_types, only : co2_opts_t, column_co2_budget_t, cohort_co2_flux_t, HR_EXP_ED2
   implicit none
   private

   public :: canopy_air_co2_update, aggregate_cohort_co2_fluxes
   public :: heterotrophic_respiration_flux, column_co2_step

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
                                         ustar, co2_atm, rho_air, dt, budget)
      real(wp), intent(inout) :: can_co2                    !< [umol/mol]  prognostic third twin
      real(wp), intent(in)    :: can_depth                  !< [m]         CAS depth (fixed within a step at MVP)
      real(wp), intent(in)    :: can_shv                    !< [kg/kg]     CAS humidity (-> dry-air molar density)
      real(wp), intent(in)    :: gross_primary_prod         !< [umol/m2/s] patch GPP (uptake, >= 0)
      real(wp), intent(in)    :: plant_respiration          !< [umol/m2/s] leaf+stem+root+growth+storage (source)
      real(wp), intent(in)    :: heterotrophic_respiration  !< [umol/m2/s] soil Rh (source, >= 0)
      real(wp), intent(in)    :: ustar                      !< [m/s]       friction velocity (shared twin conductance)
      real(wp), intent(in)    :: co2_atm                    !< [umol/mol]  free-atmosphere reference
      real(wp), intent(in)    :: rho_air                    !< [kg/m3]     moist CAS air density
      real(wp), intent(in)    :: dt                         !< [s]
      type(column_co2_budget_t), intent(out) :: budget      !< nee/nep/loss2atm/storage/resid
      real(wp) :: can_dmol, ccapcan, cci, f_bio, gatm_co2, co2_new

      can_dmol = rho_air * (1.0_wp - can_shv) / mmdry        ! [mol_dryair/m3] DRY-air molar density (exact)
      ccapcan  = can_dmol * can_depth                        ! [mol_air/m2]    MOLAR CAS capacity (cf wcapcan)
      cci      = 1.0_wp / max(ccapcan, tiny_num)             ! [m2/mol]
      f_bio    = heterotrophic_respiration + plant_respiration - gross_primary_prod  ! [umol/m2/s] Reco - GPP
      gatm_co2 = can_dmol * ustar                            ! [mol_air/m2/s]  atm<->CAS molar exchange (c3 = 1)
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
   ! MVP heterotrophic respiration [umol/m2/s]: a frozen decomposable pool x an intrinsic     !
   ! decay rate x a dimensionless temperature modifier x a dimensionless moisture modifier.    !
   ! ED2 default scheme-0 (capped exponential) OR Collatz Q10, selectable via opts%hr_model.    !
   ! The soil-carbon pool is read-only in the fast loop (prescribed at P0); this is the fast     !
   ! coupling flux from the slow pool into the CAS twin. (DAMM lands at P1 as HR_DAMM.)           !
   !---------------------------------------------------------------------------------------!
   pure function heterotrophic_respiration_flux(fast_soil_carbon, soil_temp, theta,          &
                                                theta_dry, theta_sat, opts) result(rh)
      real(wp), intent(in) :: fast_soil_carbon      !< [kgC/m2]  decomposable pool (frozen in the fast loop)
      real(wp), intent(in) :: soil_temp             !< [K]       representative (root-weighted) soil temperature
      real(wp), intent(in) :: theta, theta_dry, theta_sat   !< [m3/m3] moisture, air-dry floor, porosity
      type(co2_opts_t), intent(in) :: opts
      real(wp) :: rh                                !< [umol/m2/s]
      real(wp) :: f_temp, f_water, rel

      rel     = min(1.0_wp, max(0.0_wp, (theta - theta_dry) / max(theta_sat - theta_dry, tiny_num)))
      f_water = water_modifier(rel, opts)
      select case (opts%hr_model)
      case (HR_EXP_ED2)                             ! ED2 scheme-0: min(1, exp(a*(T - T_sat)))
         f_temp = min(1.0_wp, exp(opts%resp_temp_increase * (soil_temp - opts%resp_temp_ref)))
      case default                                  ! HR_Q10 (Collatz/K13): q10**((T - T_ref)/10)
         f_temp = opts%rh_q10 ** ((soil_temp - opts%rh_t_ref) * 0.1_wp)
      end select
      rh = fast_soil_carbon * opts%rh_k_base * f_temp * f_water * kgCday_2_umols   ! kgC/m2/day -> umol/m2/s
   end function heterotrophic_respiration_flux

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
   subroutine column_co2_step(cas_can_co2, can_depth, can_shv, ustar, co2_atm, rho_air, dt,   &
                              n, nplant, leaf_area, a_gross, rd, stem_resp, root_resp,         &
                              growth_resp_committed, storage_resp_committed,                   &
                              fast_soil_carbon, soil_temp, theta, theta_dry, theta_sat,        &
                              opts, budget)
      real(wp),    intent(inout) :: cas_can_co2            !< [umol/mol] = cas%can_co2 (passed by reference)
      real(wp),    intent(in)    :: can_depth, can_shv, ustar, co2_atm, rho_air, dt
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
                                 plant_resp, hetero, ustar, co2_atm, rho_air, dt, budget)
      !----- Closed-budget guard (mixed rtol/atol form; cf. soil_energy_flux). --------------------!
      scale = max(abs(budget%storage), 1.0_wp)
      if (opts%debug_error .and. abs(budget%resid) > opts%rtol * scale + opts%atol) then
         error stop 'column_co2_step: CO2 budget did not close'
      end if
   end subroutine column_co2_step

end module meds_column_co2
