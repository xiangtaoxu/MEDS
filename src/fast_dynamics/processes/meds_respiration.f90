! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_respiration -- the patch's respiration besides the leaf's: the plants' maintenance  !
! respiration of stems and fine roots, per plant, at rates referenced to 25 degC           !
! (docs/science/plant_respiration.md), and the soil's heterotrophic respiration from the   !
! decomposing carbon pools (docs/science/soil_carbon.md eqs. (2)-(4)):                     !
!   * stem_maintenance_respiration      -- stem surface and sapwood terms, eq. (3);        !
!   * fine_root_maintenance_respiration -- per unit fine-root carbon, eq. (4), with the    !
!                                          layer-summed scale root_zone_temp_scale;        !
!   * decomposition_env_scalar          -- each pool's temperature, moisture and lignin    !
!                                          scalar xi;                                      !
!   * decomposition_transfer_matrix     -- the pools' decay rates K, the routing A of the  !
!                                          decayed carbon and the respired fractions er;   !
!   * heterotrophic_respiration         -- Rh = sum_j er_j xi_j K_j X_j.                   !
! Leaf Rd is the leaf solver's, growth respiration the daily allocator's, and the daily    !
! step of the pools meds_soil_biogeochem's.                                                !
!==========================================================================================!
module meds_respiration
   use meds_kinds,             only : wp, ik
   use meds_constants,         only : pi, tiny_num, yr_day
   use meds_therm_lib,         only : peaked_arrhenius_scale
   use meds_patch_state_vars,  only : soil_carbon_t, pack_pool_vector
   use meds_biogeochem_params, only : soil_carbon_params_t, DECOMP_SCHEME_CENTURY5, N_SOIL_POOL, IP_FAST_GRND,      &
                                      IP_FAST_SOIL, IP_STRUCT_GRND, IP_STRUCT_SOIL, IP_MICR, IP_SLOW, IP_PASSIVE
   implicit none
   private

   public :: stem_maintenance_respiration, fine_root_maintenance_respiration, root_zone_temp_scale
   public :: decomposition_env_scalar, decomposition_transfer_matrix, heterotrophic_respiration

contains

   !========== Plant (autotrophic) maintenance respiration (plant_respiration.md) =========!

   !---------------------------------------------------------------------------------------!
   ! Stem maintenance respiration [umol CO2/plant/s], docs/science/plant_respiration.md    !
   ! eq. (3). Elemental: a scalar call does one cohort, an array call a patch, with params !
   ! broadcast. A grass (is_woody = .false.) respires nothing here.                        !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine stem_maintenance_respiration(wood_temp, dbh, height, wai, nplant,          &
                                                         aboveground_frac, is_woody, resp_factor25,     &
                                                         sap_volume, resp_sapwood25, ea, hd, ds, stem_resp)
      real(wp),            intent(in)  :: wood_temp   !< [K]  woody-tissue temperature
      real(wp),            intent(in)  :: dbh         !< [cm] stem diameter at breast height
      real(wp),            intent(in)  :: height      !< [m]  cohort height
      real(wp),            intent(in)  :: wai         !< [m2/m2 ground] wood area index
      real(wp),            intent(in)  :: nplant      !< [plant/m2] stem density
      real(wp),            intent(in)  :: aboveground_frac  !< [--] cohort PFT's aboveground fraction
      logical,             intent(in)  :: is_woody          !< cohort PFT is woody (grass => 0)
      real(wp),            intent(in)  :: resp_factor25     !< [umol CO2/m2 stem/s @25C] cohort PFT's surface rate
      !----- The sapwood term: the living wood respires by its volume, coarse roots included (§2). ---!
      real(wp),            intent(in)  :: sap_volume        !< [m3/plant] sapwood volume
      real(wp),            intent(in)  :: resp_sapwood25    !< [umol CO2/m3 sapwood/s @25C] cohort PFT's sapwood rate
      !----- Leaf Rd's temperature parameters, always in the peaked form (plant_respiration.md §4). -!
      real(wp),            intent(in)  :: ea          !< [J/mol]   activation energy
      real(wp),            intent(in)  :: hd          !< [J/mol]   deactivation energy
      real(wp),            intent(in)  :: ds          !< [J/mol/K] entropy term
      real(wp),            intent(out) :: stem_resp   !< [umol CO2 / plant / s]
      real(wp) :: tscale, stem_area

      if (.not. is_woody) then
         stem_resp = 0.0_wp
         return
      end if

      !----- Peaked temperature response (= 1 at 25 degC) on leaf Rd's parameters. -----------!
      tscale = peaked_arrhenius_scale(1.0_wp, ea, hd, ds, wood_temp)
      !----- Per-plant stem surface: the cylinder plus the WAI branch term (WAI is per ground, !
      !      hence /nplant), over the aboveground fraction to cover the whole woody plant.      !
      stem_area = ( pi * (dbh * 1.0e-2_wp) * height                                             &
                  + pi * wai / max(nplant, tiny_num) ) / max(aboveground_frac, tiny_num)
      !----- The surface and sapwood terms share the temperature response (eq. (3)); a PFT with   !
      !      no sapwood rate respires by its surface alone. --------------------------------------!
      stem_resp = resp_factor25 * tscale * stem_area + resp_sapwood25 * tscale * sap_volume
   end subroutine stem_maintenance_respiration

   !---------------------------------------------------------------------------------------!
   ! Fine-root maintenance respiration [umol CO2/plant/s], docs/science/plant_respiration.md !
   ! eq. (4). Elemental; the patch-uniform tscale_root broadcasts over the cohorts.        !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine fine_root_maintenance_respiration(tscale_root, broot, resp_factor25,      &
                                                              root_resp)
      real(wp),            intent(in)  :: tscale_root !< [-] root-weighted temperature scale (root_zone_temp_scale)
      real(wp),            intent(in)  :: broot       !< [kgC/plant] fine-root biomass
      real(wp),            intent(in)  :: resp_factor25 !< [umol CO2/kgC root/s @25C] cohort PFT's baseline
      real(wp),            intent(out) :: root_resp   !< [umol CO2 / plant / s]
      root_resp = resp_factor25 * tscale_root * broot
   end subroutine fine_root_maintenance_respiration

   !---------------------------------------------------------------------------------------!
   ! The root-weighted temperature scale sum_k root_frac_k * f(T_k): summed over layers,   !
   ! not f at the mean temperature, whose error flips sign with the season                 !
   ! (docs/science/plant_respiration.md §3). Patch-uniform: evaluated once per patch.      !
   !---------------------------------------------------------------------------------------!
   pure function root_zone_temp_scale(soil_temp, root_frac, n, ea, hd, ds) result(tscale)
      real(wp),            intent(in) :: soil_temp(:)  !< [K] per-layer soil temperature
      real(wp),            intent(in) :: root_frac(:)  !< [-] per-layer root fraction (sums to 1)
      integer(ik),         intent(in) :: n             !< active layers
      real(wp),            intent(in) :: ea, hd, ds    !< [J/mol], [J/mol], [J/mol/K] leaf Rd's response
      real(wp)    :: tscale
      integer(ik) :: k
      tscale = 0.0_wp
      do k = 1_ik, n
         tscale = tscale + root_frac(k)                                                          &
                * peaked_arrhenius_scale(1.0_wp, ea, hd, ds, soil_temp(k))
      end do
   end function root_zone_temp_scale

   !========== Soil (heterotrophic) respiration (soil_carbon.md) ==========================!

   !=======================================================================================!
   !  ENVIRONMENTAL SCALAR xi = temperature x moisture x oxygen, per pool. Grnd pools (1,3) use   !
   !  the SURFACE temperature (ED2 A_decomp); below pools (2,4,5,6,7) the ACTIVE-layer temperature  !
   !  (B_decomp). Structural pools (3,4) additionally carry the lignin brake L = exp(-e_lignin*      !
   !  f_lignin). Temperature form: scheme 5 = Collatz Q10 (UNCAPPED; xi can exceed 1), else ED2         !
   !  scheme-0 capped exp min(1, exp(a*(T-T_sat))).                                                      !
   !  Moisture is the one-sided exponential about resp_opt_water (dry side) / anoxia (wet side).         !
   !  The fast loop's sub-daily Rh calls this same routine, so xi and the fast Rh reconcile at           !
   !  matched inputs. Invariant: xi > 0 (scheme-0 additionally <= 1). For spin-up, evaluate at           !
   !  climatological-mean T/theta to build xi_bar. docs/science/soil_carbon.md eqs. (3)-(4).             !
   !=======================================================================================!
   pure subroutine decomposition_env_scalar(soil_temp_surf, soil_temp_active, theta_active, theta_dry, &
                                       theta_sat, pools, params, xi)
      real(wp),                   intent(in)  :: soil_temp_surf, soil_temp_active     !< [K]
      real(wp),                   intent(in)  :: theta_active, theta_dry, theta_sat   !< [m3/m3]
      type(soil_carbon_t),        intent(in)  :: pools                               !< for f_lignin
      type(soil_carbon_params_t), intent(in)  :: params
      real(wp),                   intent(out) :: xi(N_SOIL_POOL)                     !< [-] diagonal scalar, > 0
      real(wp) :: f_temp_surf, f_temp_active, f_water, rel, f_lig_g, f_lig_s, brake_g, brake_s

      f_temp_surf   = decomp_temperature_factor(soil_temp_surf,   params)
      f_temp_active = decomp_temperature_factor(soil_temp_active, params)
      rel     = min(1.0_wp, max(0.0_wp, (theta_active - theta_dry) / max(theta_sat - theta_dry, tiny_num)))
      f_water = decomp_water_factor(rel, params)

      f_lig_g = lignin_fraction(pools%struct_grnd_lignin, pools%struct_grnd_carbon)
      f_lig_s = lignin_fraction(pools%struct_soil_lignin, pools%struct_soil_carbon)
      brake_g = exp(-params%e_lignin * f_lig_g)
      brake_s = exp(-params%e_lignin * f_lig_s)

      xi(IP_FAST_GRND)   = f_temp_surf   * f_water
      xi(IP_FAST_SOIL)   = f_temp_active * f_water
      xi(IP_STRUCT_GRND) = f_temp_surf   * f_water * brake_g
      xi(IP_STRUCT_SOIL) = f_temp_active * f_water * brake_s
      xi(IP_MICR)        = f_temp_active * f_water
      xi(IP_SLOW)        = f_temp_active * f_water
      xi(IP_PASSIVE)     = f_temp_active * f_water
   end subroutine decomposition_env_scalar

   !----- Decomposition temperature factor: scheme-5 UNCAPPED Q10; else ED2 capped exp. ----------!
   pure function decomp_temperature_factor(soil_temp, params) result(f_temp)
      real(wp),                   intent(in) :: soil_temp
      type(soil_carbon_params_t), intent(in) :: params
      real(wp) :: f_temp
      if (params%decomp_scheme == DECOMP_SCHEME_CENTURY5) then         ! Collatz/K13 Q10 (uncapped, > 1 above ref)
         f_temp = params%rh_q10 ** ((soil_temp - params%rh_t_ref) * 0.1_wp)
      else                                                           ! ED2 scheme-0 capped exp, min(1, .)
         f_temp = min(1.0_wp, exp(params%resp_temp_increase * (soil_temp - params%resp_temp_ref)))
      end if
   end function decomp_temperature_factor

   !----- One-sided moisture/oxygen modifier about resp_opt_water (SAME form as the fast Rh). ----!
   pure function decomp_water_factor(rel, params) result(f_water)
      real(wp),                   intent(in) :: rel
      type(soil_carbon_params_t), intent(in) :: params
      real(wp) :: f_water
      if (rel <= params%resp_opt_water) then
         f_water = exp((rel - params%resp_opt_water) * params%resp_water_below_opt)      ! dry side
      else
         f_water = exp((params%resp_opt_water - rel) * params%resp_water_above_opt)      ! wet side (anoxia)
      end if
   end function decomp_water_factor

   pure function lignin_fraction(lignin, carbon) result(f_lignin)
      real(wp), intent(in) :: lignin, carbon
      real(wp) :: f_lignin
      if (carbon > tiny_num) then
         f_lignin = min(1.0_wp, max(0.0_wp, lignin / carbon))
      else
         f_lignin = 0.0_wp                                            ! C=0 => f_lignin undefined, read as 0
      end if
   end function lignin_fraction

   !=======================================================================================!
   !  TRANSFER MATRIX A + decay diagonal K + respired fractions er. A(i,j) is the fraction of pool  !
   !  j's decayed carbon entering pool i, a_jj = -1 UNIFORMLY, er_j = 1 - sum_{i!=j} a_ij. K is in    !
   !  [1/day] (params carries [1/yr]; converted by yr_day here). Default (scheme 0-4): metabolic litter  !
   !  respires 100%, structural sends (1-er_struct) to slow (er_struct lignin-weighted), slow respires  !
   !  100%, microbial/passive INERT (K=0). Scheme 5 (Bolker/Koven): fast->micr, struct non-lignin->micr  !
   !  / lignin->slow, micr<->slow and slow<->passive back-transfers with clay/sand-controlled routing.    !
   !=======================================================================================!
   pure subroutine decomposition_transfer_matrix(pools, params, a_mat, k_diag, er)
      type(soil_carbon_t),        intent(in)  :: pools
      type(soil_carbon_params_t), intent(in)  :: params
      real(wp),                   intent(out) :: a_mat(N_SOIL_POOL, N_SOIL_POOL)   !< [-] a_ii=-1, a_ij transfer frac
      real(wp),                   intent(out) :: k_diag(N_SOIL_POOL)               !< [1/day] baseline decay rates
      real(wp),                   intent(out) :: er(N_SOIL_POOL)                   !< [-] respired fraction 1 - sum_i a_ij
      real(wp)    :: f_lig_g, f_lig_s, er_stg, er_sts, er_micr, ex_mp, ex_sp
      integer(ik) :: j

      !----- decay diagonal K [1/day] (params in /yr). -----------------------------------------------!
      k_diag(IP_FAST_GRND)   = params%k_fast   / yr_day
      k_diag(IP_FAST_SOIL)   = params%k_fast   / yr_day
      k_diag(IP_STRUCT_GRND) = params%k_struct / yr_day
      k_diag(IP_STRUCT_SOIL) = params%k_struct / yr_day
      k_diag(IP_MICR)        = params%k_micr   / yr_day
      k_diag(IP_SLOW)        = params%k_slow   / yr_day
      k_diag(IP_PASSIVE)     = params%k_pass   / yr_day

      !----- lignin-weighted structural respired fractions. ----------------------------------------!
      f_lig_g = lignin_fraction(pools%struct_grnd_lignin, pools%struct_grnd_carbon)
      f_lig_s = lignin_fraction(pools%struct_soil_lignin, pools%struct_soil_carbon)
      er_stg  = (1.0_wp - f_lig_g) * params%er_struct_nonlig + f_lig_g * params%er_struct_lig
      er_sts  = (1.0_wp - f_lig_s) * params%er_struct_nonlig + f_lig_s * params%er_struct_lig

      !----- uniform a_jj = -1 (all decayed carbon leaves the donor), off-diagonals zeroed. --------!
      a_mat = 0.0_wp
      do j = 1_ik, N_SOIL_POOL
         a_mat(j, j) = -1.0_wp
      end do

      if (params%decomp_scheme == DECOMP_SCHEME_CENTURY5) then
         !----- FULL CENTURY (5-active, Bolker/Koven): texture-controlled topology. ----------------!
         er_micr = params%er_micr_int + params%er_micr_slp * params%xsand             ! sand -> microbial respired frac
         ex_mp   = min(1.0_wp - er_micr, params%fx_micr_pass_int + params%fx_micr_pass_slp * params%xclay)
         ex_sp   = min(1.0_wp - params%er_slow, params%fx_slow_pass_int + params%fx_slow_pass_slp * params%xclay)
         !----- fast litter (1,2) -> microbial. ---------------------------------------------------!
         a_mat(IP_MICR, IP_FAST_GRND) = 1.0_wp - params%er_fast
         a_mat(IP_MICR, IP_FAST_SOIL) = 1.0_wp - params%er_fast
         !----- structural (3,4): non-lignin -> microbial, lignin -> slow. -------------------------!
         a_mat(IP_MICR, IP_STRUCT_GRND) = (1.0_wp - er_stg) * (1.0_wp - f_lig_g)
         a_mat(IP_MICR, IP_STRUCT_SOIL) = (1.0_wp - er_sts) * (1.0_wp - f_lig_s)
         a_mat(IP_SLOW, IP_STRUCT_GRND) = (1.0_wp - er_stg) * f_lig_g
         a_mat(IP_SLOW, IP_STRUCT_SOIL) = (1.0_wp - er_sts) * f_lig_s
         !----- microbial (5): -> passive (clay), remainder -> slow. -------------------------------!
         a_mat(IP_PASSIVE, IP_MICR) = ex_mp
         a_mat(IP_SLOW,    IP_MICR) = (1.0_wp - er_micr) - ex_mp
         !----- slow (6): -> passive (clay), remainder -> microbial. -------------------------------!
         a_mat(IP_PASSIVE, IP_SLOW) = ex_sp
         a_mat(IP_MICR,    IP_SLOW) = (1.0_wp - params%er_slow) - ex_sp
         !----- passive (7): -> slow. -------------------------------------------------------------!
         a_mat(IP_SLOW, IP_PASSIVE) = 1.0_wp - params%er_pass
      else
         !----- DEFAULT (schemes 0-4, 3-active): structural -> slow, everything else terminal. -----!
         a_mat(IP_SLOW, IP_STRUCT_GRND) = 1.0_wp - er_stg
         a_mat(IP_SLOW, IP_STRUCT_SOIL) = 1.0_wp - er_sts
      end if

      !----- respired fraction er_j = 1 - sum_{i!=j} a_ij  (= -sum_i a_ij, since a_jj=-1). ---------!
      do j = 1_ik, N_SOIL_POOL
         er(j) = 1.0_wp - (sum(a_mat(:, j)) + 1.0_wp)      ! sum_i a_ij = a_jj(-1) + sum_{i!=j}; er = 1 - sum_{i!=j}
      end do
   end subroutine decomposition_transfer_matrix

   !=======================================================================================!
   !  INSTANTANEOUS heterotrophic respiration Rh = sum_j er_j*xi_j*K_j*X_j = -1^T*A*xi*K*X            !
   !  [kgC/m2/day] -- the respired complement A does NOT transfer (the carbon that leaves the network). !
   !  This is ALSO what the fast loop respires into the canopy air, per sub-step, against the frozen    !
   !  daily pool -- which is what makes the day's atmospheric Rh equal the daily step's pool debit by    !
   !  construction (the rh_seam_gap contract). Callers crossing into the CAS must convert with           !
   !  kgcday_2_umols; patch_heterotrophic_respiration is where that happens.                             !
   !  docs/science/soil_carbon.md eq. (2) and §7.                                                        !
   !=======================================================================================!
   pure function heterotrophic_respiration(a_mat, k_diag, xi, pools) result(rh)
      real(wp),            intent(in) :: a_mat(N_SOIL_POOL, N_SOIL_POOL), k_diag(N_SOIL_POOL)
      real(wp),            intent(in) :: xi(N_SOIL_POOL)
      type(soil_carbon_t), intent(in) :: pools
      real(wp) :: rh, x(N_SOIL_POOL), er_j
      integer(ik) :: j
      call pack_pool_vector(pools, x)
      rh = 0.0_wp
      do j = 1_ik, N_SOIL_POOL
         er_j = -sum(a_mat(:, j))                    ! column sum of A is -er_j (a_jj=-1)
         rh   = rh + er_j * xi(j) * k_diag(j) * x(j)
      end do
   end function heterotrophic_respiration

end module meds_respiration
