!==========================================================================================!
! meds_canopy_aerodynamics -- THE stateless in-canopy aerodynamics kernel (design Part I). A     !
! pure function of (free-atmosphere forcing, canopy-air-space state, canopy geometry) that        !
! produces every turbulence/conductance quantity the coupled fast loop's other kernels take as    !
! forced inputs: friction velocity `ustar` + the scalar profile factors temp1/temp2 (which set    !
! the atm<->CAS conductances gah/gaw/gac for all three CAS twins), the per-cohort in-canopy wind   !
! and leaf/wood boundary-layer conductances, and the ground<->CAS conductance ggnet (=> the soil-  !
! evaporation resistance r_aero = 1/ggnet). It carries NO integrated state -- like ED2's           !
! canopy_turbulence8, it is recomputed each fast substep from the current CAS + vegetation state.  !
!                                                                                          !
! Physics (see the design's ED2<->CLM reconciliation): the surface layer uses CLM5's clean         !
! four-range Monin-Obukhov denominators with a fixed-iteration solve (no root-find, no I/O -> GPU  !
! portable); the leaf/wood boundary layers use ED2's forced+free convection Nusselt kernel (needed !
! for low-wind free convection); the in-canopy wind uses ED2's per-cohort crown-area extinction;   !
! the ground conductance uses the CLM4-like dense-canopy form. All `pure`; per-cohort outputs are  !
! written into caller-owned SoA slices (never an array-valued function result into a call, #7).    !
!==========================================================================================!
module meds_canopy_aerodynamics
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : grav, pi, tiny_num
   use meds_biophysics_types, only : aero_cfg_t, aero_env_t, aero_geom_t, aero_out_t
   implicit none
   private

   public :: canopy_aerodynamics             !< master per-patch seam
   public :: mo_surface_layer               !< CLM Monin-Obukhov solve (exposed for unit tests)
   public :: reduced_wind                    !< log-profile wind at a height (exposed for tests)
   public :: boundary_gbh_mos                !< Nusselt boundary-layer conductance (exposed for tests)

contains

   !=======================================================================================!
   !  Master seam: fill an aero_out_t for one patch. Cohort arrays are ordered BOTTOM(1) ->    !
   !  TOP(n) (the canopy-radiation convention); the wind cascade walks top -> bottom.          !
   !=======================================================================================!
   pure subroutine canopy_aerodynamics(cfg, env, geom, n, height, lai, crown_area,        &
                                            leaf_temp, wood_temp, leaf_width, branch_diam, out)
      type(aero_cfg_t),  intent(in)    :: cfg
      type(aero_env_t),  intent(in)    :: env
      type(aero_geom_t), intent(in)    :: geom
      integer(ik),       intent(in)    :: n
      real(wp),          intent(in)    :: height(n), lai(n), crown_area(n)
      real(wp),          intent(in)    :: leaf_temp(n), wood_temp(n), leaf_width(n), branch_diam(n)
      type(aero_out_t),  intent(inout) :: out

      real(wp)    :: rough, displace, zldis, uh, ca, ext_half, ext_full, stab, denom
      integer(ik) :: ico

      !----- 1. Roughness + displacement from canopy height (snow blend). ------------------!
      rough    = (1.0_wp - geom%snowfac) * cfg%z0m_ratio * geom%veg_height                     &
                 + geom%snowfac * cfg%snow_rough
      rough    = max(rough, cfg%snow_rough)
      displace = cfg%d_ratio * geom%veg_height * (1.0_wp - geom%snowfac)
      zldis    = max(env%zref - displace, 2.0_wp * rough)
      out%rough = rough
      out%displace = displace
      out%can_depth = max(cfg%min_canopy_depth, geom%veg_height)
      out%n_coh = n

      !----- 2. CLM Monin-Obukhov surface layer (canopy z0h = z0m). ------------------------!
      call mo_surface_layer(cfg, env%u_ref, env%zref, displace, rough, env%theta_atm,          &
                            env%shv_atm, env%can_theta, env%can_shv, out%ustar, out%temp1,      &
                            out%zeta, out%rib, out%obu)
      out%temp2 = out%temp1                                        ! z0q = z0h => scalar factors equal
      out%tstar = out%temp1 * (env%theta_atm - env%can_theta)
      out%qstar = out%temp2 * (env%shv_atm   - env%can_shv)
      out%cstar = out%temp2 * (env%co2_atm   - env%can_co2)

      !----- 3. Ground <-> CAS conductances (bare similarity + CLM4 dense-canopy blend). ---!
      out%ggbare = out%ustar * out%temp1                           ! [m/s] bare-ground scalar conductance
      denom      = (env%can_temp + env%t_ground) * out%ustar ** 2
      stab       = 0.0_wp
      if (denom > tiny_num) then
         stab = max(0.0_wp, min(10.0_wp, 2.0_wp * grav * geom%veg_height                       &
                    * (env%can_temp - env%t_ground) / denom))
      end if
      out%ggveg = cfg%cs_dense * out%ustar / (1.0_wp + cfg%gamma_g * stab)
      if (geom%opencan_frac > 0.999_wp .or. geom%snowfac >= 0.9_wp) then
         out%ggnet = out%ggbare
      else
         out%ggnet = out%ggbare * out%ggveg                                                     &
                     / (out%ggveg + (1.0_wp - geom%opencan_frac) * out%ggbare + tiny_num)
      end if

      !----- 4. In-canopy wind: top-of-canopy wind, then crown-area extinction top -> bottom. !
      uh     = reduced_wind(cfg, out%ustar, out%zeta, geom%veg_height, displace, rough, zldis)
      out%uh = uh
      do ico = n, 1_ik, -1_ik
         ca            = min(max(crown_area(ico), 0.01_wp), 1.0_wp)
         ext_half      = ca * exp(-0.25_wp * lai(ico) / ca) + (1.0_wp - ca)
         ext_full      = ca * exp(-0.50_wp * lai(ico) / ca) + (1.0_wp - ca)
         out%wind(ico) = max(cfg%ugbmin, uh * ext_half)          ! wind at crown mid-depth
         uh            = uh * ext_full                            ! attenuate for the cohort below
      end do

      !----- 5. Per-cohort leaf (flat plate) + wood (cylinder) boundary layers [m/s]. ------!
      do ico = 1_ik, n
         out%leaf_gbh(ico) = boundary_gbh_mos(cfg, out%wind(ico), leaf_temp(ico), env%can_temp, &
                                              leaf_width(ico), .true.)
         out%leaf_gbw(ico) = cfg%gbh_2_gbw * out%leaf_gbh(ico)
         out%wood_gbh(ico) = boundary_gbh_mos(cfg, out%wind(ico), wood_temp(ico), env%can_temp, &
                                              branch_diam(ico), .false.)
         out%wood_gbw(ico) = cfg%gbh_2_gbw * out%wood_gbh(ico)
      end do
   end subroutine canopy_aerodynamics

   !=======================================================================================!
   !  CLM5 Monin-Obukhov surface layer. Fixed n_iter_mo sweeps of the four-range denominators, !
   !  seeded by the analytic bulk-Richardson guess (MoninObukIni). Returns ustar + the scalar   !
   !  profile factor temp1 (= vonk/D_h; temp2 identical for z0q=z0h) + stability diagnostics.    !
   !  Neutral limit: ustar -> vonk*u/ln(zldis/z0m), temp1 -> vonk/ln(zldis/z0m). No division by  !
   !  zeta near neutral (zeta0m = z0m*zeta/zldis; zeta appears in a denominator ONLY in the very- !
   !  unstable/very-stable branches, where |zeta| > 1), so -fpe0 is safe.                        !
   !=======================================================================================!
   pure subroutine mo_surface_layer(cfg, u_ref, zref, displace, rough, theta_atm, shv_atm,     &
                                    can_theta, can_shv, ustar, temp1, zeta, rib, obu)
      type(aero_cfg_t), intent(in)  :: cfg
      real(wp),         intent(in)  :: u_ref, zref, displace, rough, theta_atm, shv_atm
      real(wp),         intent(in)  :: can_theta, can_shv
      real(wp),         intent(out) :: ustar, temp1, zeta, rib, obu
      real(wp)    :: uref, zldis, z0m, thv_atm, thv_can, dthv, dth, dqh, um, thvstar
      integer(ik) :: it

      uref    = max(u_ref, cfg%ubmin)
      z0m     = rough
      zldis   = max(zref - displace, 2.0_wp * z0m)
      thv_atm = theta_atm * (1.0_wp + 0.61_wp * shv_atm)
      thv_can = can_theta * (1.0_wp + 0.61_wp * can_shv)
      dthv    = thv_atm - thv_can
      dth     = theta_atm - can_theta
      dqh     = shv_atm   - can_shv

      !----- MoninObukIni: convective velocity if unstable, then analytic Rib -> zeta guess. -!
      if (dthv >= 0.0_wp) then
         um = max(uref, 0.1_wp)
      else
         um = sqrt(uref * uref + cfg%wc * cfg%wc)
      end if
      rib = grav * zldis * dthv / (thv_atm * um * um)
      if (rib >= 0.0_wp) then
         zeta = rib * log(zldis / z0m) / (1.0_wp - 5.0_wp * min(rib, 0.19_wp))
         zeta = min(cfg%zeta_max_stable, max(zeta, 0.01_wp))
      else
         zeta = rib * log(zldis / z0m)
         zeta = max(-100.0_wp, min(zeta, -0.01_wp))
      end if

      !----- Fixed-iteration solve (no data-dependent exit -> warp-uniform on GPU). --------!
      ustar = cfg%ustmin
      temp1 = cfg%vonk / log(zldis / z0m)
      do it = 1_ik, cfg%n_iter_mo
         ustar   = max(cfg%ustmin, cfg%vonk * um / d_mom(zeta, zldis, z0m, cfg))
         temp1   = cfg%vonk / d_heat(zeta, zldis, z0m, cfg)
         thvstar = (temp1 * dth) * (1.0_wp + 0.61_wp * shv_atm) + 0.61_wp * theta_atm * (temp1 * dqh)
         zeta    = zldis * cfg%vonk * grav * thvstar / (ustar * ustar * thv_atm)
      end do
      !----- Final ustar/temp1 consistent with the converged zeta. -------------------------!
      ustar = max(cfg%ustmin, cfg%vonk * um / d_mom(zeta, zldis, z0m, cfg))
      temp1 = cfg%vonk / d_heat(zeta, zldis, z0m, cfg)
      obu   = zldis / sign(max(abs(zeta), 1.0e-6_wp), zeta)
   end subroutine mo_surface_layer

   !----- Momentum profile denominator D_m (ustar = vonk*um/D_m), four CLM5 ranges. ----------!
   pure function d_mom(zeta, zldis, z0m, cfg) result(dm)
      real(wp),         intent(in) :: zeta, zldis, z0m
      type(aero_cfg_t), intent(in) :: cfg
      real(wp) :: dm, zeta0m
      zeta0m = z0m * zeta / zldis                                  ! = z0m / L  (no division by zeta)
      if (zeta < -cfg%zeta_m) then                                 ! very unstable
         dm = log(-cfg%zeta_m * zldis / (zeta * z0m)) - stabfunc1(-cfg%zeta_m) + stabfunc1(zeta0m) &
              + 1.14_wp * ((-zeta) ** 0.333_wp - cfg%zeta_m ** 0.333_wp)
      else if (zeta < 0.0_wp) then                                 ! unstable
         dm = log(zldis / z0m) - stabfunc1(zeta) + stabfunc1(zeta0m)
      else if (zeta <= 1.0_wp) then                                ! stable (Businger-Dyer)
         dm = log(zldis / z0m) + 5.0_wp * zeta - 5.0_wp * zeta0m
      else                                                         ! very stable (log cap)
         dm = log(zldis / (zeta * z0m)) + 5.0_wp - 5.0_wp * zeta0m + (5.0_wp * log(zeta) + zeta - 1.0_wp)
      end if
   end function d_mom

   !----- Heat/vapour profile denominator D_h (temp1 = vonk/D_h), four CLM5 ranges. ----------!
   pure function d_heat(zeta, zldis, z0h, cfg) result(dh)
      real(wp),         intent(in) :: zeta, zldis, z0h
      type(aero_cfg_t), intent(in) :: cfg
      real(wp) :: dh, zeta0h
      zeta0h = z0h * zeta / zldis
      if (zeta < -cfg%zeta_t) then                                 ! very unstable
         dh = log(-cfg%zeta_t * zldis / (zeta * z0h)) - stabfunc2(-cfg%zeta_t) + stabfunc2(zeta0h) &
              + 0.8_wp * (cfg%zeta_t ** (-0.333_wp) - (-zeta) ** (-0.333_wp))
      else if (zeta < 0.0_wp) then                                 ! unstable
         dh = log(zldis / z0h) - stabfunc2(zeta) + stabfunc2(zeta0h)
      else if (zeta <= 1.0_wp) then                                ! stable
         dh = log(zldis / z0h) + 5.0_wp * zeta - 5.0_wp * zeta0h
      else                                                         ! very stable
         dh = log(zldis / (zeta * z0h)) + 5.0_wp - 5.0_wp * zeta0h + (5.0_wp * log(zeta) + zeta - 1.0_wp)
      end if
   end function d_heat

   !----- CLM StabilityFunc1 (psi_m, UNSTABLE branch; call only with zeta < 0). --------------!
   pure function stabfunc1(zeta) result(psi)
      real(wp), intent(in) :: zeta
      real(wp) :: psi, x
      x   = (1.0_wp - 16.0_wp * zeta) ** 0.25_wp
      psi = 2.0_wp * log((1.0_wp + x) * 0.5_wp) + log((1.0_wp + x * x) * 0.5_wp)                &
            - 2.0_wp * atan(x) + pi * 0.5_wp
   end function stabfunc1

   !----- CLM StabilityFunc2 (psi_h, UNSTABLE branch; call only with zeta < 0). --------------!
   pure function stabfunc2(zeta) result(psi)
      real(wp), intent(in) :: zeta
      real(wp) :: psi, x2
      x2  = sqrt(1.0_wp - 16.0_wp * zeta)
      psi = 2.0_wp * log((1.0_wp + x2) * 0.5_wp)
   end function stabfunc2

   !----- Integrated momentum stability function for the WIND PROFILE (both signs). ----------!
   pure function psim(zeta) result(p)
      real(wp), intent(in) :: zeta
      real(wp) :: p
      if (zeta < 0.0_wp) then
         p = stabfunc1(zeta)
      else
         p = -5.0_wp * zeta                                        ! Businger-Dyer stable (profile use)
      end if
   end function psim

   !=======================================================================================!
   !  Log-profile wind at `height` from the surface-layer solution (floored at ugbmin). Used   !
   !  for the canopy-top wind; monotonically increasing with height in neutral/stable air.      !
   !=======================================================================================!
   pure function reduced_wind(cfg, ustar, zeta, height, displace, rough, zldis) result(u)
      type(aero_cfg_t), intent(in) :: cfg
      real(wp),         intent(in) :: ustar, zeta, height, displace, rough, zldis
      real(wp) :: u, hd, zeta_h, zeta_0
      hd     = max(height - displace, rough * 1.001_wp)
      zeta_h = zeta * hd    / zldis
      zeta_0 = zeta * rough / zldis
      u      = (ustar / cfg%vonk) * (log(hd / rough) - psim(zeta_h) + psim(zeta_0))
      u      = max(u, cfg%ugbmin)
   end function reduced_wind

   !=======================================================================================!
   !  ED2 leaf (flat-plate) / wood (cylinder) boundary-layer conductance gbh_mos [m/s], the     !
   !  parallel sum of forced (Reynolds) and free (Grashof) convection Nusselt numbers. This is   !
   !  the PRIMITIVE the energy kernel consumes directly (its gbh/gbw are m/s velocities); the     !
   !  driver forms gbw = gbh_2_gbw*gbh_mos and, for the leaf gas-exchange seam, gb_mol = gbw*rho. !
   !=======================================================================================!
   pure function boundary_gbh_mos(cfg, wind, elem_temp, can_temp, lchar, is_leaf) result(gbh_mos)
      type(aero_cfg_t), intent(in) :: cfg
      real(wp),         intent(in) :: wind, elem_temp, can_temp, lchar
      logical,          intent(in) :: is_leaf
      real(wp) :: gbh_mos, nu_v, al_th, re, gr, nu_forced, nu_free
      nu_v  = cfg%kin_visc0 * (1.0_wp + cfg%dkin_visc * (can_temp - cfg%t_ref_air))
      al_th = cfg%th_diff0  * (1.0_wp + cfg%dth_diff  * (can_temp - cfg%t_ref_air))
      re    = max(wind * lchar / al_th, tiny_num)
      gr    = max(grav / (can_temp * nu_v * nu_v) * abs(elem_temp - can_temp) * lchar ** 3, 0.0_wp)
      if (is_leaf) then                                            ! flat plate
         nu_forced = max(cfg%aflat_lami * re ** cfg%nflat_lami, cfg%aflat_turb * re ** cfg%nflat_turb)
         nu_free   = max(cfg%bflat_lami * gr ** cfg%mflat_lami, cfg%bflat_turb * gr ** cfg%mflat_turb)
      else                                                         ! cylinder
         nu_forced = max(cfg%ocyli_lami + cfg%acyli_lami * re ** cfg%ncyli_lami,                &
                         cfg%ocyli_turb + cfg%acyli_turb * re ** cfg%ncyli_turb)
         nu_free   = max(cfg%bcyli_lami * gr ** cfg%mcyli_lami, cfg%bcyli_turb * gr ** cfg%mcyli_turb)
      end if
      gbh_mos = max(cfg%gbhmos_min, (nu_forced + nu_free) * al_th / lchar)
   end function boundary_gbh_mos

end module meds_canopy_aerodynamics
