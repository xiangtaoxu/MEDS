!==========================================================================================!
! meds_plant_hydraulics -- the plant-hydraulics NETWORK SOLVER: the coupled matrix-exponential   !
! sub-step integrator (solve_plant_water, + its batch wrapper solve_plant_water_batch) that        !
! assembles the tissue CONSTITUTIVE curves (pressure-volume + Kirchhoff conductance, now in         !
! meds_hydr_lib) into a 2-node leaf<->wood ODE, plus the optional soil->root rhizosphere        !
! conductance helper. Every fast-loop integrator (split/ARK/RK45) now calls this ONCE per macro-  !
! step as the Act-1 pre-pass (MEDS_ED2_RK45_DESIGN.md sec 1/4/5) and freezes its returned time-     !
! averaged sapflow/root_uptake for the whole step -- the retired plant_water_tendency (a per-stage  !
! psi RHS) is no longer needed now that internal water MASS, not psi, is the fast-loop prognostic   !
! state (mass's own ODE has no self-feedback stiffness, so it needs no per-stage re-linearization). !
! solve_plant_water is re-exported through meds_plant_interface.                !
!==========================================================================================!
module meds_plant_hydraulics
   use meds_kinds,     only : wp, ik
   use meds_constants, only : pi, grav_head, safe_exp, tiny_num
   use meds_plant_types,      only : hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t,  &
                                     N_HYDRO, NODE_LEAF, NODE_WOOD, NROOT_MAX,                  &
                                     HYDRO_NODES_2, HYDRO_COND_SEGMENT, HYDRO_SUBSTEP_FIXED
   use meds_hydr_lib,      only : kirchhoff_edge, kirchhoff_edge_tab, capacitance,           &
                                     water_content, plc_retained
   use meds_numerics,         only : adaptive_step_update
   implicit none
   private

   public :: solve_plant_water, solve_plant_water_batch, rhizosphere_cond
   public :: root_fraction_profile, effective_root_boundary

   real(wp),    parameter :: c_floor   = 1.0e-12_wp  !< capacitance floor (linearization only)
   real(wp),    parameter :: k_floor   = 1.0e-15_wp  !< conductance floor (keeps M finite)
   real(wp),    parameter :: sinhc_eps = 1.0e-8_wp   !< sinhc series threshold
   real(wp),    parameter :: safety    = 0.9_wp
   real(wp),    parameter :: fmax      = 4.0_wp
   real(wp),    parameter :: fmin      = 0.2_wp
   real(wp),    parameter :: h_floor   = 1.0e-6_wp   !< accept regardless below this sub-step [s]


contains

   !----- Wood<->leaf edge conductance: the EXACT closed form for wood_kexp in {1,2}, else the     !
   !       precomputed lookup table (kirchhoff_edge_tab). At kexp in {1,2} the table is never        !
   !       consulted, so existing (kexp=2) runs are bit-identical; the general-kexp branch trades    !
   !       the 7-pt quadrature for a linear lookup. Pure => callable from the ARK RHS; the stale-     !
   !       table guard lives in solve_plant_water (non-pure, so it can error stop).                    !
   pure real(wp) function edge_cond(psi_up, psi_down, k_cond, p) result(keff)
      real(wp),             intent(in) :: psi_up, psi_down, k_cond
      type(hydro_params_t), intent(in) :: p
      if (abs(p%wood_kexp - 1.0_wp) < 1.0e-9_wp .or. abs(p%wood_kexp - 2.0_wp) < 1.0e-9_wp) then
         keff = kirchhoff_edge(psi_up, psi_down, k_cond, p%wood_psi50, p%wood_kexp)
      else
         keff = kirchhoff_edge_tab(psi_up, psi_down, k_cond, p%wood_psi50, p%vuln_table)
      end if
   end function edge_cond

   !----- Per-layer soil->root (rhizosphere) conductance [kg/s/MPa] per plant (Katul et al. 2003),   !
   !      from soil hydraulic conductivity, fine-root biomass, specific root area, and the layer's    !
   !      root fraction. This is ED2's `gw_cond`; call it per soil layer to fill env%rhizo_cond_layer  !
   !      for the multi-layer root boundary (MEDS_MULTILAYER_ROOTS_DESIGN).                            !
   pure real(wp) function rhizosphere_cond(soil_cond, broot, sra, root_frac, dz, nplant) result(gw)
      real(wp), intent(in) :: soil_cond   !< [kg/m/s/MPa] soil hydraulic conductivity (per MPa gradient)
      real(wp), intent(in) :: broot       !< [kgC] fine-root biomass (per plant)
      real(wp), intent(in) :: sra         !< [m2/kgC] specific root area
      real(wp), intent(in) :: root_frac   !< [-] fraction of roots in this layer
      real(wp), intent(in) :: dz          !< [m] layer thickness
      real(wp), intent(in) :: nplant      !< [pl/m2] plant density
      real(wp) :: rai
      rai = broot * sra * root_frac * nplant                     ! root area index [m2/m2]
      gw  = soil_cond * sqrt(max(rai, 0.0_wp)) / (pi * dz) / max(nplant, tiny(1.0_wp))
   end function rhizosphere_cond

   !----- ED2 cumulative-exponential root fraction in a soil layer spanning depths [z_top, z_bot]     !
   !      (both >= 0, below surface, z_bot > z_top): frac = beta^(z_top/D) - beta^(z_bot/D), with      !
   !      D = root_depth and beta in (0,1). Shallower layers get more roots; depths clamp to [0, D] so  !
   !      layers below the rooting depth contribute 0. Summed over [0, D] the profile telescopes to     !
   !      1 - beta (ED2 convention; only the RELATIVE distribution enters the conductance weights).     !
   pure real(wp) function root_fraction_profile(root_beta, root_depth, z_top, z_bot) result(frac)
      real(wp), intent(in) :: root_beta   !< [-]  root-profile decay (0,1); smaller => shallower
      real(wp), intent(in) :: root_depth  !< [m]  maximum rooting depth (> 0)
      real(wp), intent(in) :: z_top       !< [m]  depth of the layer top    (>= 0)
      real(wp), intent(in) :: z_bot       !< [m]  depth of the layer bottom (> z_top)
      real(wp) :: inv_d, a, b
      inv_d = 1.0_wp / max(root_depth, tiny_num)
      a = min(max(z_top, 0.0_wp), root_depth) * inv_d
      b = min(max(z_bot, 0.0_wp), root_depth) * inv_d
      frac = root_beta**a - root_beta**b
   end function root_fraction_profile

   !----- Reduce the root boundary to an effective (conductance, soil potential) at the wood node.    !
   !      Multi-layer (n_root_layer > 1): the per-layer conductances are in PARALLEL to the common     !
   !      wood node, so they collapse EXACTLY to G = sum_k g_k and psi_eff = sum_k g_k*(psi_k +         !
   !      grav_head*z_k) / G (conductance-weighted mean soil potential, gravity head to each layer).    !
   !      Single-layer (n_root_layer <= 1): return the scalar rhizo_cond/soil_psi BC verbatim, so the   !
   !      solver is bit-identical to the pre-multilayer path. Pure => usable from the ARK RHS; the       !
   !      k_floor clamp on the conductance stays at the call site (matching the legacy code).            !
   pure subroutine effective_root_boundary(env, rhizo_eff, soil_psi_eff)
      type(hydro_env_t), intent(in)  :: env
      real(wp),          intent(out) :: rhizo_eff, soil_psi_eff
      real(wp)    :: g_tot, gpsi
      integer(ik) :: k
      if (env%n_root_layer > 1_ik) then
         g_tot = 0.0_wp ; gpsi = 0.0_wp
         do k = 1_ik, env%n_root_layer
            g_tot = g_tot + env%rhizo_cond_layer(k)
            gpsi  = gpsi  + env%rhizo_cond_layer(k)                                                &
                          * (env%soil_psi_layer(k) + grav_head*env%root_z_layer(k))
         end do
         rhizo_eff    = g_tot
         soil_psi_eff = gpsi / max(g_tot, k_floor)
      else
         rhizo_eff    = env%rhizo_cond
         soil_psi_eff = env%soil_psi
      end if
   end subroutine effective_root_boundary


   !========== meds_hydro_solver.f90 ====================================================!

   !---------------------------------------------------------------------------------------!
   ! Advance the node potentials psi(1:2) = (leaf, wood) over dt. env/p/o are read-only; psi   !
   ! is inout; flux carries the outputs. Only 2-node (leaf+wood) is implemented for now.         !
   !---------------------------------------------------------------------------------------!
   subroutine solve_plant_water(env, p, o, dt, psi, flux)
      type(hydro_env_t),    intent(in)    :: env
      type(hydro_params_t), intent(in)    :: p
      type(hydro_opts_t),   intent(in)    :: o
      real(wp),             intent(in)    :: dt
      real(wp),             intent(inout) :: psi(N_HYDRO)
      type(hydro_flux_t),   intent(out)   :: flux
      !----- Locals frozen over dt (boundary conditions + geometry). ----------------------!
      real(wp)    :: e_transp, soil_psi, rhz, grav, k_cond, bwood, supply_tot
      real(wp)    :: psi_l0, psi_w0, psi_l, psi_w, dw_l, dw_w
      real(wp)    :: h, t_rem, xl, xw, ml, mw, xl2, xw2, err, errl, errw
      real(wp)    :: fm11, fm12, fm21, fm22, fps_l, fps_w   ! start-state frozen coeffs shared by full + half step
      integer(ik) :: nsub, i, nfix
      !------------------------------------------------------------------------------------!

      if (o%topology /= HYDRO_NODES_2) then
         call fatal_hydro('solve_plant_water: only HYDRO_NODES_2 is implemented (3-node is a follow-up).')
      end if

      !----- Multi-layer root boundary must fit the fixed NROOT_MAX arrays (else the aggregation and    !
      !       per-layer distribution loops read out of bounds). The caller fills n_root_layer. ----------!
      if (env%n_root_layer > NROOT_MAX) then
         call fatal_hydro('solve_plant_water: n_root_layer exceeds NROOT_MAX (grow NROOT_MAX in meds_plant_types).')
      end if

      !----- Stale/absent-table guard: when kexp is not in {1,2} the solver consults p%vuln_table,     !
      !       which must have been BUILT (build_hydro_table) for THIS wood_kexp. inv_h>0 iff built      !
      !       (it is NTAB/r_max); testing inv_h too closes the wood_kexp==0.0 default coincidence,       !
      !       where kexp matches an all-zero unbuilt table. Catches a post-derive wood_kexp mutation     !
      !       or a partially-initialised hydro_params_t (edge_cond is pure and cannot check).            !
      if (abs(p%wood_kexp - 1.0_wp) >= 1.0e-9_wp .and. abs(p%wood_kexp - 2.0_wp) >= 1.0e-9_wp) then
         if (p%vuln_table%inv_h <= 0.0_wp .or. abs(p%vuln_table%kexp - p%wood_kexp) > 1.0e-9_wp)      &
            call fatal_hydro('solve_plant_water: vuln_table not built for this wood_kexp (stale/absent table).')
      end if

      !----- Boundary conditions / geometry, held constant over dt. ------------------------!
      e_transp = env%transp
      call effective_root_boundary(env, rhz, soil_psi) ! multi-layer aggregate OR scalar BC (identical if <=1)
      rhz      = max(rhz, k_floor)                     ! ground the network (avoid singular M)
      bwood    = env%bsap + env%broot
      if (o%gravity_on) then ; grav = grav_head*env%height ; else ; grav = 0.0_wp ; end if
      k_cond   = cond_max()                            ! max (plc=1) internal conductance [kg/s/MPa]

      psi_l0 = psi(NODE_LEAF) ; psi_w0 = psi(NODE_WOOD)
      psi_l  = psi_l0         ; psi_w  = psi_w0
      nsub   = 0_ik

      if (o%substep_mode == HYDRO_SUBSTEP_FIXED) then
         !----- Fixed equal sub-steps (GPU lockstep). ------------------------------------!
         nfix = max(o%max_substep, 1_ik)
         h    = dt / real(nfix, wp)
         do i = 1_ik, nfix
            call exact_substep(psi_l, psi_w, h, xl, xw)
            psi_l = xl ; psi_w = xw ; nsub = nsub + 1_ik
         end do
         flux%converged = .true.
      else
         !----- Adaptive step-doubling. --------------------------------------------------!
         t_rem = dt
         if (o%h_init > 0.0_wp) then ; h = min(o%h_init, dt) ; else ; h = dt ; end if
         flux%converged = .true.
         do
            if (t_rem <= tiny_num) exit
            h = min(h, t_rem)
            !----- Full step and first half step share the SAME start state (psi_l, psi_w), so     !
            !       their frozen coefficients are identical -- compute once, reuse for both. -------!
            call freeze_coeffs(psi_l, psi_w, fm11, fm12, fm21, fm22, fps_l, fps_w)
            call advance_exact_linear(psi_l, psi_w, h,        fm11, fm12, fm21, fm22, fps_l, fps_w, xl,  xw )  ! full step
            call advance_exact_linear(psi_l, psi_w, 0.5_wp*h, fm11, fm12, fm21, fm22, fps_l, fps_w, ml,  mw )  ! first half
            call exact_substep(ml,    mw,    0.5_wp*h,   xl2, xw2)                                       ! second half
            errl = abs(xl2 - xl) / (o%atol + o%rtol*abs(xl2))
            errw = abs(xw2 - xw) / (o%atol + o%rtol*abs(xw2))
            err  = max(errl, errw, 1.0e-12_wp)
            if (err <= 1.0_wp .or. h <= h_floor) then
               if (err > 1.0_wp) flux%converged = .false.          ! floor-forced accept below tolerance
               psi_l = xl2 ; psi_w = xw2 ; t_rem = t_rem - h ; nsub = nsub + 1_ik
               h   = h * adaptive_step_update(err, safety, fmin, fmax)   ! accept: grow (p=1 controller)
            else
               h   = h * adaptive_step_update(err, safety, fmin, fmax)   ! reject: shrink (same clamp)
            end if
            if (nsub >= o%max_substep .and. t_rem > tiny_num) then  ! cap hit before finishing interval
               flux%converged = .false. ; exit
            end if
         end do
      end if

      !----- Mass-closing boundary fluxes from the converged storage change. ---------------!
      dw_l = water_content(psi_l, p%leaf_pi0, p%leaf_elastic_mod, p%leaf_apoplast_frac, p%leaf_water_sat, env%bleaf) &
           - water_content(psi_l0,p%leaf_pi0, p%leaf_elastic_mod, p%leaf_apoplast_frac, p%leaf_water_sat, env%bleaf)
      dw_w = water_content(psi_w, p%wood_pi0, p%wood_elastic_mod, p%wood_apoplast_frac, p%wood_water_sat, bwood)      &
           - water_content(psi_w0,p%wood_pi0, p%wood_elastic_mod, p%wood_apoplast_frac, p%wood_water_sat, bwood)

      flux%sapflow     = dw_l/dt + e_transp
      flux%root_uptake = (dw_l + dw_w)/dt + e_transp

      !----- Distribute the total root uptake to soil layers, proportional to each layer's POSITIVE     !
      !       supply max(g_k*(psi_soil_k + grav_head*z_k - psi_wood), 0), so the per-layer uptake sums    !
      !       EXACTLY to root_uptake and NO layer effluxes. A dry layer that would give a negative supply !
      !       (root->soil flux) is floored to 0 rather than effluxing: hydraulic redistribution (HR) is   !
      !       intentionally NOT enabled here -- see docs/dev_plans/MEDS_MULTILAYER_ROOTS_DESIGN.md (HR      !
      !       deferred to a future version). Single-layer path assigns it all to layer 1. ----------------!
      flux%root_uptake_layer = 0.0_wp
      if (env%n_root_layer > 1_ik) then
         supply_tot = 0.0_wp
         do i = 1_ik, env%n_root_layer
            supply_tot = supply_tot + max(env%rhizo_cond_layer(i)                                   &
                       * (env%soil_psi_layer(i) + grav_head*env%root_z_layer(i) - psi_w), 0.0_wp)
         end do
         if (supply_tot > tiny_num) then
            do i = 1_ik, env%n_root_layer
               flux%root_uptake_layer(i) = max(env%rhizo_cond_layer(i)                              &
                    * (env%soil_psi_layer(i) + grav_head*env%root_z_layer(i) - psi_w), 0.0_wp)      &
                    / supply_tot * flux%root_uptake
            end do
         end if
      else
         flux%root_uptake_layer(1) = flux%root_uptake
      end if

      flux%psi_leaf    = psi_l
      flux%psi_wood    = psi_w
      flux%plc         = 1.0_wp - plc_retained(psi_w, p%wood_psi50, p%wood_kexp)
      flux%nsub        = nsub

      psi(NODE_LEAF) = psi_l
      psi(NODE_WOOD) = psi_w

   contains

      !----- Maximum (plc=1) internal conductance, per plant [kg/s/MPa]. -------------------!
      real(wp) function cond_max() result(kc)
         if (o%cond_mode == HYDRO_COND_SEGMENT) then
            kc = p%wood_kmax * env%sap_area / max(env%height*p%vessel_curl, tiny_num)
         else
            kc = p%k_plant_max * env%leaf_area
         end if
         kc = max(kc, k_floor)
      end function cond_max

      !----- Frozen linear-system coefficients + Ohm's-law steady state at the sub-step start   !
      !       state (pl, pw). h-INDEPENDENT, so the adaptive full + first-half steps share them. !
      subroutine freeze_coeffs(pl, pw, m11, m12, m21, m22, ps_l, ps_w)
         real(wp), intent(in)  :: pl, pw
         real(wp), intent(out) :: m11, m12, m21, m22, ps_l, ps_w
         real(wp) :: cl, cw, keff, c1, c2, detm
         cl   = max(capacitance(pl, p%leaf_pi0, p%leaf_elastic_mod, p%leaf_apoplast_frac, p%leaf_water_sat, env%bleaf), c_floor)
         cw   = max(capacitance(pw, p%wood_pi0, p%wood_elastic_mod, p%wood_apoplast_frac, p%wood_water_sat, bwood),     c_floor)
         keff = max(edge_cond(pw, pl, k_cond, p), k_floor)
         !----- Linear system psi' = M psi + c. -------------------------------------------!
         m11 = -keff/cl        ; m12 =  keff/cl
         m21 =  keff/cw        ; m22 = -(keff + rhz)/cw
         c1  = (-keff*grav - e_transp)/cl
         c2  = ( keff*grav + rhz*soil_psi)/cw
         detm = m11*m22 - m12*m21                            ! = keff*rhz/(cl*cw) > 0
         ps_l = -( m22*c1 - m12*c2)/detm                     ! psi* = -M^{-1} c (Ohm's-law steady state)
         ps_w =  ( m21*c1 - m11*c2)/detm
      end subroutine freeze_coeffs

      !----- Advance one exact sub-step of length hs under pre-computed FROZEN coefficients.  !
      subroutine advance_exact_linear(pl, pw, hs, m11, m12, m21, m22, ps_l, ps_w, plo, pwo)
         real(wp), intent(in)  :: pl, pw, hs, m11, m12, m21, m22, ps_l, ps_w
         real(wp), intent(out) :: plo, pwo
         real(wp) :: n11, n12, n21, n22, mu, dl, ep, em, ch, sh, e11, e12, e21, e22, ddl, ddw
         !----- e^{Mh} via the underflow-safe sinhc form (real eigenvalues <= 0). ----------!
         n11 = m11*hs ; n12 = m12*hs ; n21 = m21*hs ; n22 = m22*hs
         mu  = 0.5_wp*(n11 + n22)
         dl  = sqrt(max(mu*mu - (n11*n22 - n12*n21), 0.0_wp))
         ep  = safe_exp(mu + dl) ; em = safe_exp(mu - dl)    ! both <= 1
         ch  = 0.5_wp*(ep + em)
         if (dl > sinhc_eps) then ; sh = 0.5_wp*(ep - em)/dl ; else ; sh = ch ; end if
         e11 = ch + sh*(n11 - mu) ; e12 = sh*n12
         e21 = sh*n21             ; e22 = ch + sh*(n22 - mu)
         ddl = pl - ps_l ; ddw = pw - ps_w
         plo = ps_l + e11*ddl + e12*ddw
         pwo = ps_w + e21*ddl + e22*ddw
      end subroutine advance_exact_linear

      !----- One exact frozen-coefficient 2x2 matrix-exponential sub-step (freeze + apply).  !
      subroutine exact_substep(pl, pw, hs, plo, pwo)
         real(wp), intent(in)  :: pl, pw, hs
         real(wp), intent(out) :: plo, pwo
         real(wp) :: m11, m12, m21, m22, ps_l, ps_w
         call freeze_coeffs(pl, pw, m11, m12, m21, m22, ps_l, ps_w)
         call advance_exact_linear(pl, pw, hs, m11, m12, m21, m22, ps_l, ps_w, plo, pwo)
      end subroutine exact_substep

   end subroutine solve_plant_water

   !---------------------------------------------------------------------------------------!
   ! solve_plant_water_batch -- BARE-ARRAY wrapper over solve_plant_water for a whole patch's n   !
   ! cohorts (MEDS_NUMERICS_SCOPING.md BB1 phase 2: the per-cohort inner loop of column_fast_step   !
   ! is genuinely independent across cohorts -- no cross-cohort coupling within one dt_fast sub-     !
   ! step -- so this is the natural offload-eligible seam. The per-cohort PHYSICS is completely      !
   ! UNCHANGED: this does not re-derive anything, it just calls the SAME validated solve_plant_water  !
   ! once per cohort inside a plain `do i=1,n` loop, so the result is bit-identical to the caller's   !
   ! old inline loop.                                                                                 !
   !                                                                                          !
   ! Every dummy argument is a BARE array (or a scalar broadcast), mirroring meds_core_state_update's  !
   ! proven `!$omp target`-eligible pattern (CLAUDE.md: "it takes bare arrays (no site_t, no derived   !
   ! types), so the map clauses are clean") -- NOT a derived-type bundle. The caller passes CONTIGUOUS  !
   ! slices of its own (possibly capacity-oversized, BB1 phase 1) backing arrays, e.g. coh%bleaf(1:n);   !
   ! every array here is declared to exactly the active extent (n cohorts / nsl layers), so passing an   !
   ! oversized backing array unsliced would be WRONG -- always slice at the call site.                   !
   !                                                                                          !
   ! Root boundary condition (mirrors henv%n_root_layer's two branches in the original inline loop):     !
   !   * multilayer_roots = .false. (default, bit-identical single-BC path): soil_psi_scalar/            !
   !     rhizo_cond_scalar broadcast to every cohort; soil_psi_layer/root_z_layer/rhizo_cond_layer        !
   !     are unused (the caller may pass them un-filled).                                                 !
   !   * multilayer_roots = .true.: soil_psi_layer(nsl)/root_z_layer(nsl) are PER-LAYER, the SAME for      !
   !     every cohort (soil state does not vary by cohort); rhizo_cond_layer(nsl,n) is genuinely           !
   !     per-(layer,cohort) (it depends on each cohort's broot/nplant) and the caller must precompute      !
   !     it (rhizosphere_cond is cheap + already `pure`, itself a candidate for the same treatment).       !
   !                                                                                          !
   ! NOT YET GPU-OFFLOADED: solve_plant_water is NOT `pure` (it can `error stop` on a stale vulnerability  !
   ! table -- meds_plant_hydraulics.f90's fatal_hydro), which blocks a literal `!$omp target` on this      !
   ! loop today; device code generally cannot execute error-termination/I-O statements. Landing this as    !
   ! a bare-array CPU-serial refactor first (this pass) separates that follow-up (making the stale-table   !
   ! guard a caller-side precondition check instead) from the data-layout change, matching BB1's own       !
   ! "land it bit-identical and still serial first" discipline.                                            !
   !---------------------------------------------------------------------------------------!
   subroutine solve_plant_water_batch(n, nsl, multilayer_roots, transp, bleaf, bsap, broot, sap_area, &
                                      height, leaf_area, soil_psi_scalar, rhizo_cond_scalar,          &
                                      soil_psi_layer, root_z_layer, rhizo_cond_layer, p, o, dt, psi,  &
                                      sapflow, root_uptake, root_uptake_layer, psi_leaf, psi_wood,    &
                                      plc, nsub, converged)
      integer(ik),          intent(in)    :: n, nsl
      logical,              intent(in)    :: multilayer_roots
      real(wp),             intent(in)    :: transp(n), bleaf(n), bsap(n), broot(n), sap_area(n)
      real(wp),             intent(in)    :: height(n), leaf_area(n)
      real(wp),             intent(in)    :: soil_psi_scalar, rhizo_cond_scalar    !< broadcast (single-BC path)
      real(wp),             intent(in)    :: soil_psi_layer(nsl), root_z_layer(nsl)      !< per-layer, all cohorts
      real(wp),             intent(in)    :: rhizo_cond_layer(nsl, n)                    !< per-(layer,cohort)
      type(hydro_params_t), intent(in)    :: p
      type(hydro_opts_t),   intent(in)    :: o
      real(wp),             intent(in)    :: dt
      real(wp),             intent(inout) :: psi(N_HYDRO, n)
      real(wp),             intent(out)   :: sapflow(n), root_uptake(n), root_uptake_layer(nsl, n)
      real(wp),             intent(out)   :: psi_leaf(n), psi_wood(n), plc(n)
      integer(ik),          intent(out)   :: nsub(n)
      logical,              intent(out)   :: converged(n)

      type(hydro_env_t)  :: env
      type(hydro_flux_t) :: flux
      integer(ik) :: i, k

      !----- Independent per-cohort solve (no cross-iteration dependency): the offload-eligible loop.  !
      do i = 1_ik, n
         env%transp    = transp(i)
         env%bleaf     = bleaf(i) ; env%bsap = bsap(i) ; env%broot = broot(i)
         env%sap_area  = sap_area(i) ; env%height = height(i) ; env%leaf_area = leaf_area(i)
         if (multilayer_roots) then
            env%n_root_layer = nsl
            do k = 1_ik, nsl
               env%soil_psi_layer(k)   = soil_psi_layer(k)
               env%root_z_layer(k)     = root_z_layer(k)
               env%rhizo_cond_layer(k) = rhizo_cond_layer(k, i)
            end do
         else
            env%n_root_layer = 0_ik
            env%soil_psi      = soil_psi_scalar
            env%rhizo_cond    = rhizo_cond_scalar
         end if
         call solve_plant_water(env, p, o, dt, psi(:, i), flux)
         sapflow(i)     = flux%sapflow
         root_uptake(i) = flux%root_uptake
         psi_leaf(i)    = flux%psi_leaf
         psi_wood(i)    = flux%psi_wood
         plc(i)         = flux%plc
         nsub(i)        = flux%nsub
         converged(i)   = flux%converged
         if (multilayer_roots) then
            do k = 1_ik, nsl
               root_uptake_layer(k, i) = flux%root_uptake_layer(k)
            end do
         end if
      end do
   end subroutine solve_plant_water_batch

   !----- Catchable fatal (error stop), matching the MEDS convention. ----------------------!
   subroutine fatal_hydro(msg)
      character(len=*), intent(in) :: msg
      write(*,'(2a)') 'meds_plant_hydraulics: ', msg
      error stop 1
   end subroutine fatal_hydro

end module meds_plant_hydraulics
