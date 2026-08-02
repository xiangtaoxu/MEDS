!==========================================================================================!
! meds_diagnostic_kernels -- pure/elemental DERIVED diagnostic quantities.                  !
!                                                                                          !
! Stage [1] of the diagnostic wall (MEDS_IO_V01_PLAN.md section 3.1): quantities that are not      !
! stored anywhere in state but are a closed-form function of things that are. Plain arguments      !
! only -- no site_t, no config aggregator -- so every routine here is directly unit-testable.      !
!                                                                                          !
! CONTENT RULE. Anything a physics library already owns is CALLED, never re-implemented:           !
! soil_psi_from_theta / psi_from_water_content live in meds_hydr_lib, uext_to_temp and the moist-  !
! air thermodynamics in meds_therm_lib, the allometry in meds_allometry. This module holds only     !
! diagnostics that no physics kernel owns. Where a thin re-export is more readable than making      !
! every caller `use` two modules, the wrapper is marked as such and simply forwards.                !
!                                                                                          !
! PLACEMENT. src/io/ (the meds_io_prep target, which links meds_core only -- no netCDF), NOT        !
! src/shared/functions/. Accepted consequence: nothing BELOW io in the library DAG can call these.  !
! That is correct rather than merely tolerated -- a diagnostic is by definition something no        !
! physics kernel needs, and putting it below the wall would be the first step toward physics        !
! depending on its own reporting layer. If a quantity here ever turns out to be needed by physics,  !
! that is the signal it was never a diagnostic and belongs in the owning physics library instead.   !
!==========================================================================================!
module meds_diagnostic_kernels
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num, gsw_2_gsc, grav_head, mmdry, r_gas, r_wv
   use meds_therm_lib, only : sat_vapor_pressure
   use meds_hydr_lib,  only : soil_psi_from_theta
   implicit none
   private

   public :: cohort_lai, cohort_npp_per_plant, cohort_gsc, cohort_wue, cohort_ci_ca
   public :: soil_wetness, soil_matric_potential, air_vpd, specific_humidity_to_vpd
   public :: bowen_ratio, safe_ratio
   public :: dbh_class_index

   !----- Dry-air / water-vapour gas-constant ratio eps = R_d/R_v (~0.622), the mixing-ratio    !
   !      constant in q <-> e. Derived from the shared constants rather than re-stated, so it     !
   !      can never drift from the thermodynamics the model actually integrates.                  !
   real(wp), parameter :: EPS_MOL = (r_gas / mmdry) / r_wv

contains

   !=======================================================================================!
   !  Structural / carbon                                                                   !
   !=======================================================================================!

   !----- Per-cohort leaf area index [m2 leaf / m2 ground of its own patch]. ----------------!
   pure elemental real(wp) function cohort_lai(nplant, leaf_area) result(lai)
      real(wp), intent(in) :: nplant     !< [plant/m2]
      real(wp), intent(in) :: leaf_area  !< [m2/plant]
      lai = nplant * leaf_area
   end function cohort_lai

   !----- Per-plant NPP [kgC/plant] over whatever period the accumulators span: GROSS GPP    !
   !      minus the three autotrophic MAINTENANCE respiration terms. Growth respiration is     !
   !      NOT included here -- it is charged inside the slow-loop allocation step, and adding   !
   !      it would double-count against npp_*_site. Mirrors total_npp's definition exactly.     !
   pure elemental real(wp) function cohort_npp_per_plant(gpp, leaf_resp, stem_resp, root_resp) &
                                    result(npp)
      real(wp), intent(in) :: gpp, leaf_resp, stem_resp, root_resp   !< [kgC/plant]
      npp = gpp - leaf_resp - stem_resp - root_resp
   end function cohort_npp_per_plant

   !=======================================================================================!
   !  Leaf gas exchange                                                                     !
   !=======================================================================================!

   !----- Stomatal conductance to CO2 from the conductance to water vapour (Fick: CO2        !
   !      diffuses 1.6x slower than H2O through the stomatal pore).                            !
   pure elemental real(wp) function cohort_gsc(gsw) result(gsc)
      real(wp), intent(in) :: gsw   !< [mol H2O/m2 leaf/s]
      gsc = gsw / gsw_2_gsc
   end function cohort_gsc

   !----- Instantaneous water-use efficiency [umol CO2 / mmol H2O]. The transpiration input   !
   !      is the leaf kernel's [mol H2O/m2/s], so the 1e3 converts mol -> mmol. Guarded: at     !
   !      night E -> 0 and the ratio is meaningless, so it returns 0 rather than a huge spike.  !
   pure elemental real(wp) function cohort_wue(a_net, transpiration) result(wue)
      real(wp), intent(in) :: a_net          !< [umol CO2/m2 leaf/s]
      real(wp), intent(in) :: transpiration  !< [mol H2O/m2 leaf/s]
      wue = safe_ratio(a_net, 1.0e3_wp * transpiration)
   end function cohort_wue

   !----- Ci/Ca ratio [-], the standard water-use-strategy index. ---------------------------!
   pure elemental real(wp) function cohort_ci_ca(ci, ca) result(r)
      real(wp), intent(in) :: ci, ca   !< [umol/mol]
      r = safe_ratio(ci, ca)
   end function cohort_ci_ca

   !=======================================================================================!
   !  Soil                                                                                  !
   !=======================================================================================!

   !----- Relative saturation (ED2 `soil_wetness`) [-]: where the layer sits between residual !
   !      and saturation. Clamped to [0,1] so a solver excursion outside the constitutive       !
   !      domain reports as exactly dry/saturated rather than as an out-of-range number.        !
   pure elemental real(wp) function soil_wetness(theta, theta_res, theta_sat) result(w)
      real(wp), intent(in) :: theta, theta_res, theta_sat
      w = (theta - theta_res) / max(theta_sat - theta_res, tiny_num)
      w = min(max(w, 0.0_wp), 1.0_wp)
   end function soil_wetness

   !----- Soil matric potential [MPa] from volumetric moisture. A thin re-export of            !
   !      meds_hydr_lib's retention curve (which returns metres of head) with the unit          !
   !      conversion the output layer wants, so a caller needs one `use`, not two.              !
   !      head [m] -> MPa via the shared `grav_head` (rho_w*g), the SAME constant plant           !
   !      hydraulics uses -- so a soil psi reported here and a psi_soil the roots saw are on one   !
   !      scale by construction, not by two independent unit conversions agreeing.                 !
   pure elemental real(wp) function soil_matric_potential(retention, theta, theta_sat,         &
                                    theta_res, par_a, par_n) result(psi_mpa)
      integer(ik), intent(in) :: retention
      real(wp),    intent(in) :: theta, theta_sat, theta_res, par_a, par_n
      psi_mpa = grav_head * soil_psi_from_theta(retention, theta, theta_sat, theta_res,       &
                                                par_a, par_n)
   end function soil_matric_potential

   !=======================================================================================!
   !  Atmosphere / canopy air                                                               !
   !=======================================================================================!

   !----- Vapour-pressure deficit [Pa] from temperature and ACTUAL vapour pressure. ---------!
   pure elemental real(wp) function air_vpd(temp, e_vap) result(vpd)
      real(wp), intent(in) :: temp    !< [K]
      real(wp), intent(in) :: e_vap   !< [Pa] actual vapour pressure
      vpd = max(0.0_wp, sat_vapor_pressure(temp) - e_vap)
   end function air_vpd

   !----- Vapour-pressure deficit [Pa] from the canopy-air-space prognostic twins            !
   !      (temperature + SPECIFIC humidity) at a given pressure. This is the form the CAS      !
   !      diagnostics need: e = q*p / (eps + (1-eps)*q), eps = mmh2o/mmdry.                    !
   pure elemental real(wp) function specific_humidity_to_vpd(temp, shv, pressure) result(vpd)
      real(wp), intent(in) :: temp      !< [K]
      real(wp), intent(in) :: shv       !< [kg/kg] specific humidity
      real(wp), intent(in) :: pressure  !< [Pa]
      real(wp) :: e_vap
      e_vap = shv * pressure / max(EPS_MOL + (1.0_wp - EPS_MOL) * shv, tiny_num)
      vpd   = air_vpd(temp, e_vap)
   end function specific_humidity_to_vpd

   !=======================================================================================!
   !  Generic                                                                               !
   !=======================================================================================!

   !----- Bowen ratio H/LE [-]; 0 when the latent flux is negligible (night, frozen). -------!
   pure elemental real(wp) function bowen_ratio(h_flux, le_flux) result(b)
      real(wp), intent(in) :: h_flux, le_flux   !< [W/m2]
      b = safe_ratio(h_flux, le_flux)
   end function bowen_ratio

   !----- x/y, returning 0 when |y| is negligible. The ONE guarded division in this module:  !
   !      every ratio diagnostic routes through it so a zero denominator can never produce a   !
   !      NaN or an Inf that would then poison a period mean (and be indistinguishable from a  !
   !      real value in the output file). A diagnostic that cannot be formed is 0 here, and     !
   !      the aggregation layer's own empty-period guard emits _FillValue for a period that      !
   !      genuinely had no samples -- the two mechanisms answer different questions.             !
   pure elemental real(wp) function safe_ratio(x, y) result(r)
      real(wp), intent(in) :: x, y
      if (abs(y) > tiny_num) then ; r = x / y ; else ; r = 0.0_wp ; end if
   end function safe_ratio

   !=======================================================================================!
   !  Size-class binning (DIM_SIZE)                                                          !
   !=======================================================================================!

   !----- 1-based DBH class of one plant, ED2-style (MEDS_IO_V01_PLAN.md section 8 D4).       !
   !                                                                                          !
   !      `edges` holds n_class+1 ASCENDING boundaries. Bin k is the half-open interval          !
   !      [edges(k), edges(k+1)), EXCEPT the last, which is closed at the top so the largest      !
   !      tree in the stand is never dropped. A dbh below edges(1) is clamped into bin 1 rather   !
   !      than discarded -- silently losing the smallest cohorts would break the closure           !
   !      sum(nplant_size) == nplant_site that makes this axis trustworthy.                       !
   !                                                                                          !
   !      A cohort is NEVER split across bins: it lands whole in the bin of its mean dbh. That     !
   !      is the standard ED2 treatment, and it is what makes nplant_size a genuine stem-density   !
   !      distribution comparable to a forest inventory.                                          !
   pure integer(ik) function dbh_class_index(dbh, edges, n_class) result(k)
      real(wp),    intent(in) :: dbh
      real(wp),    intent(in) :: edges(:)   !< (n_class+1) ascending
      integer(ik), intent(in) :: n_class
      integer(ik) :: j
      k = 1_ik
      if (n_class <= 0_ik) return
      do j = n_class, 1_ik, -1_ik
         if (dbh >= edges(j)) then ; k = j ; return ; end if
      end do
   end function dbh_class_index

end module meds_diagnostic_kernels
