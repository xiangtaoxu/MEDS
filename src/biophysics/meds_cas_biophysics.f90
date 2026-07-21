!==========================================================================================!
! meds_cas_biophysics -- the CANOPY-AIR-SPACE (CAS) balance: the well-mixed sub-canopy air box    !
! and its three prognostic twins (specific enthalpy, specific humidity, CO2 mixing ratio), all    !
! sharing one capacity and the ustar-based atm<->CAS conductance. A fast diffusion/venting         !
! exchange, so it lives with the biophysics kernels (not the slow soil-carbon pools).              !
!                                                                                          !
!   * canopy_air_update      -- advance the enthalpy + humidity twins one step (implicit atm term). !
!   * canopy_air_co2_update  -- advance the CO2 twin (molar capacity; the THIRD twin).              !
!   * aggregate_cohort_co2_fluxes -- reduce per-cohort GPP/respiration to the patch-ground CO2 flux. !
! (Heterotrophic soil respiration -- a carbon-decomposition process -- lives in meds_soil_biogeochem.)!
!==========================================================================================!
module meds_cas_biophysics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, mmdry
   use meds_therm_lib,        only : cas_temp_of_enthalpy
   use meds_biophysics_types, only : column_co2_budget_t, cohort_co2_flux_t
   implicit none
   private

   public :: canopy_air_update
   public :: canopy_air_co2_update, aggregate_cohort_co2_fluxes

contains

   !---------------------------------------------------------------------------------------!
   ! Canopy-air-space enthalpy + humidity update (design 4b). Advances BOTH prognostic twins  !
   ! one step from the summed cohort/ground fluxes and the atmospheric exchange (implicit in    !
   ! the atm term for L-stability). Returns the closed-budget residual (~0).                    !
   !---------------------------------------------------------------------------------------!
   pure subroutine canopy_air_update(cas_enthalpy, cas_shv, cas_temp, can_depth,               &
                                     coh_h_flux, coh_qw_flux, coh_w_flux, coh_transp,           &
                                     ground_h_flux, ground_qw_flux, ground_w_flux, dew,         &
                                     ustar, temp1, enthalpy_atm, w_flux_ac, rho_air, dt, resid)
      real(wp), intent(inout) :: cas_enthalpy, cas_shv, cas_temp
      real(wp), intent(in)    :: can_depth
      real(wp), intent(in)    :: coh_h_flux, coh_qw_flux, coh_w_flux, coh_transp
      real(wp), intent(in)    :: ground_h_flux, ground_qw_flux, ground_w_flux, dew
      real(wp), intent(in)    :: ustar, temp1, enthalpy_atm, w_flux_ac, rho_air, dt
      real(wp), intent(out)   :: resid
      real(wp) :: wcapcan, wci, f_sens, gatm, enth_new, shv_new
      wcapcan = rho_air * can_depth                                         ! [kg/m2] CAS air mass per ground area
      wci     = 1.0_wp / max(wcapcan, tiny_num)
      f_sens  = coh_h_flux + coh_qw_flux + ground_h_flux + ground_qw_flux   ! [W/m2] into CAS from surfaces
      !----- atm<->CAS scalar conductance = rho*ustar*c3, where c3 (temp1) is the dimensionless   !
      !      Monin-Obukhov scalar-transfer coefficient from the aerodynamics solver (aero_out%      !
      !      temp1). Dropping it (c3=1) over-couples the CAS to the free atmosphere ~1/c3 (BUG8).   !
      !      The vapour twin's temp2 (== temp1 while z0q==z0h) rides in the caller-formed w_flux_ac !
      !      (= rho*ustar*temp2*(shv_atm - can_shv)). -----------------------------------------------!
      gatm    = rho_air * ustar * temp1                                    ! [kg/m2/s] atm<->CAS exchange
      !----- Enthalpy: implicit in the atmospheric-exchange term. --------------------------!
      enth_new = (cas_enthalpy + dt * wci * (f_sens + gatm * enthalpy_atm)) / (1.0_wp + dt * wci * gatm)
      !----- Specific humidity twin (explicit in the caller-formed atm vapour flux w_flux_ac;    !
      !      the live driver uses the implicit twin -- see column_fast_step). Clamp non-negative   !
      !      so a large sink cannot drive the CAS humidity below zero.  --------------------------!
      shv_new  = max(0.0_wp, cas_shv + dt * wci * (coh_w_flux + coh_transp + ground_w_flux - dew + w_flux_ac))
      resid    = wcapcan * (enth_new - cas_enthalpy)                                            &
                 - dt * (f_sens + gatm * (enthalpy_atm - enth_new))         ! = 0 by construction
      cas_enthalpy = enth_new
      cas_shv      = shv_new
      cas_temp     = cas_temp_of_enthalpy(enth_new, shv_new)
   end subroutine canopy_air_update

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

      out%gross_primary_prod = 0.0_wp
      out%leaf_respiration   = 0.0_wp
      out%stem_respiration   = 0.0_wp
      out%root_respiration   = 0.0_wp
      out%growth_respiration = 0.0_wp
      out%storage_respiration = 0.0_wp
      do i = 1_ik, n
         out%gross_primary_prod = out%gross_primary_prod + a_gross(i)  * leaf_area(i) * nplant(i)
         out%leaf_respiration   = out%leaf_respiration   + rd(i)       * leaf_area(i) * nplant(i)
         out%stem_respiration   = out%stem_respiration   + stem_resp(i)              * nplant(i)
         out%root_respiration   = out%root_respiration   + root_resp(i)              * nplant(i)
      end do
   end subroutine aggregate_cohort_co2_fluxes


end module meds_cas_biophysics
