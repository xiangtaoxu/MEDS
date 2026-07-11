!==========================================================================================!
! meds_plant_hydraulics -- plant-hydraulics COMPUTE kernels, merged into one module: the       !
! nonlinear pressure-volume curves (Bartlett/Tyree-Hammel), the Kirchhoff-integrated whole-     !
! plant / segment conductance, and the coupled matrix-exponential sub-step solver               !
! (solve_plant_water). The public seam plant_water_flux lives in meds_plant_interface.           !
!==========================================================================================!
module meds_plant_hydraulics
   use meds_kinds,     only : wp, ik
   use meds_constants, only : pi, grav_head, safe_exp, tiny_num
   use meds_plant_types,      only : hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t,  &
                                     N_HYDRO, NODE_LEAF, NODE_WOOD,                             &
                                     HYDRO_NODES_2, HYDRO_COND_SEGMENT, HYDRO_SUBSTEP_FIXED
   implicit none
   private

   !----- from meds_hydro_conductance.f90 ------------------------------------------------!

   public :: plc_retained, dplc_dpsi, flux_potential, phi_inverse, kirchhoff_edge
   public :: rhizosphere_cond

   real(wp), parameter :: dpsi_eps = 1.0e-6_wp   !< |up-down| below which K_eff -> pointwise limit
   integer(ik), parameter :: NG = 7_ik           !< Gauss-Legendre points for general-exponent Phi

   !----- 7-point Gauss-Legendre nodes/weights on [-1,1]. ----------------------------------!
   real(wp), parameter :: gl_x(NG) = [ -0.9491079123427585_wp, -0.7415311855993945_wp,        &
                                       -0.4058451513773972_wp,  0.0000000000000000_wp,        &
                                        0.4058451513773972_wp,  0.7415311855993945_wp,        &
                                        0.9491079123427585_wp ]
   real(wp), parameter :: gl_w(NG) = [  0.1294849661688697_wp,  0.2797053914892766_wp,        &
                                        0.3818300505051189_wp,  0.4179591836734694_wp,        &
                                        0.3818300505051189_wp,  0.2797053914892766_wp,        &
                                        0.1294849661688697_wp ]

   !----- from meds_hydro_pv.f90 ---------------------------------------------------------!

   public :: pv_psi_tlp, pv_rwc_tlp, rwc_from_psi, psi_from_rwc
   public :: water_content, capacitance, pv_water_cap_from_traits

   real(wp), parameter :: rwc_floor = 1.0e-4_wp   !< keep R strictly positive in the flaccid tail

   !----- from meds_hydro_solver.f90 -----------------------------------------------------!

   public :: solve_plant_water, plant_water_tendency

   real(wp),    parameter :: c_floor   = 1.0e-12_wp  !< capacitance floor (linearization only)
   real(wp),    parameter :: k_floor   = 1.0e-15_wp  !< conductance floor (keeps M finite)
   real(wp),    parameter :: sinhc_eps = 1.0e-8_wp   !< sinhc series threshold
   real(wp),    parameter :: safety    = 0.9_wp
   real(wp),    parameter :: fmax      = 4.0_wp
   real(wp),    parameter :: fmin      = 0.2_wp
   real(wp),    parameter :: h_floor   = 1.0e-6_wp   !< accept regardless below this sub-step [s]


contains

   !========== meds_hydro_conductance.f90 ===============================================!

   !----- Fraction of conductance retained (1 - PLC). psi<0, psi50<0 => r>0; clamp for psi>0. -!
   elemental real(wp) function plc_retained(psi, psi50, kexp) result(f)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r = max(psi/psi50, 0.0_wp)
      f = 1.0_wp / (1.0_wp + r**kexp)
   end function plc_retained

   !----- d(plc_retained)/d(psi); finite at psi=0 for kexp>1 (avoids the -(a/psi)f(1-f) NaN). --!
   elemental real(wp) function dplc_dpsi(psi, psi50, kexp) result(df)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r
      r  = max(psi/psi50, 0.0_wp)
      df = -kexp * r**(kexp - 1.0_wp) / (psi50 * (1.0_wp + r**kexp)**2)
   end function dplc_dpsi

   !----- Kirchhoff (matric flux) potential Phi(psi) = integral_0^psi plc ds [MPa]. ----------!
   !      Closed form for kexp in {1,2}; fixed Gauss-Legendre on [0,r] otherwise. Phi(0)=0,     !
   !      Phi(psi<0)<0, strictly increasing in psi.                                             !
   pure real(wp) function flux_potential(psi, psi50, kexp) result(phi)
      real(wp), intent(in) :: psi, psi50, kexp
      real(wp) :: r, s, u, acc
      integer(ik) :: g
      r = max(psi/psi50, 0.0_wp)                 ! r >= 0
      if (abs(kexp - 1.0_wp) < 1.0e-9_wp) then
         phi = psi50 * log(1.0_wp + r)
      else if (abs(kexp - 2.0_wp) < 1.0e-9_wp) then
         phi = psi50 * atan(r)
      else
         !----- integral_0^r du/(1+u^kexp) by 7-pt Gauss-Legendre (map [-1,1] -> [0,r]). ----!
         acc = 0.0_wp
         do g = 1_ik, NG
            u   = 0.5_wp*r*(gl_x(g) + 1.0_wp)
            acc = acc + gl_w(g) / (1.0_wp + u**kexp)
         end do
         s   = 0.5_wp*r*acc
         phi = psi50 * s
      end if
   end function flux_potential

   !----- Inverse of Phi: find psi in [psi_lo, 0] with flux_potential(psi)=phi_target. --------!
   !      Monotone => robust bisection. Used by the vertical-profile diagnostic (3-node).       !
   pure real(wp) function phi_inverse(phi_target, psi50, kexp, psi_lo) result(psi)
      real(wp), intent(in) :: phi_target, psi50, kexp, psi_lo
      real(wp) :: lo, hi, mid, fmid
      integer(ik) :: it
      lo = psi_lo ; hi = 0.0_wp
      do it = 1_ik, 60_ik
         mid  = 0.5_wp*(lo + hi)
         fmid = flux_potential(mid, psi50, kexp)
         if (fmid < phi_target) then ; lo = mid ; else ; hi = mid ; end if
      end do
      psi = 0.5_wp*(lo + hi)
   end function phi_inverse

   !----- Kirchhoff edge conductance K_eff = k_cond * <plc> [kg/s/MPa]. k_cond is the maximum    !
   !      (plc=1) whole-plant/segment conductance already scaled to per-plant [kg/s/MPa]. The     !
   !      |dpsi|->0 limit is the pointwise value (L'Hopital of DeltaPhi/Deltapsi).                 !
   pure real(wp) function kirchhoff_edge(psi_up, psi_down, k_cond, psi50, kexp) result(keff)
      real(wp), intent(in) :: psi_up, psi_down, k_cond, psi50, kexp
      real(wp) :: dpsi
      dpsi = psi_up - psi_down
      if (abs(dpsi) > dpsi_eps) then
         keff = k_cond * ( flux_potential(psi_up,   psi50, kexp)                              &
                         - flux_potential(psi_down, psi50, kexp) ) / dpsi
      else
         keff = k_cond * plc_retained(0.5_wp*(psi_up + psi_down), psi50, kexp)
      end if
   end function kirchhoff_edge

   !----- Optional TEST helper: Katul-2003 rhizosphere (soil->root) conductance [kg/s/MPa] per   !
   !      plant, from soil conductance, fine-root biomass and specific root area. Production      !
   !      fills rhizo_cond from a future soil module; do not use on the hot path.                 !
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


   !========== meds_hydro_pv.f90 ========================================================!

   !----- Turgor loss point [MPa] (Bartlett eqn 1). ---------------------------------------!
   elemental real(wp) function pv_psi_tlp(pi0, eps) result(psi_tlp)
      real(wp), intent(in) :: pi0, eps
      psi_tlp = pi0*eps / (pi0 + eps)
   end function pv_psi_tlp

   !----- Symplastic relative water content at turgor loss (Bartlett eqn 2). --------------!
   elemental real(wp) function pv_rwc_tlp(pi0, eps) result(rwc_tlp)
      real(wp), intent(in) :: pi0, eps
      rwc_tlp = (pi0 + eps) / eps
   end function pv_rwc_tlp

   !----- Water potential [MPa] from symplastic RWC. --------------------------------------!
   elemental real(wp) function psi_from_rwc(rwc, pi0, eps) result(psi)
      real(wp), intent(in) :: rwc, pi0, eps
      real(wp) :: rwc_tlp, r
      rwc_tlp = (pi0 + eps) / eps
      r       = max(rwc, rwc_floor)
      if (r >= rwc_tlp) then
         psi = eps*(r - rwc_tlp) + pi0/r      ! turgor + osmotic
      else
         psi = pi0/r                          ! turgor lost
      end if
   end function psi_from_rwc

   !----- Symplastic RWC from water potential (closed-form turgid inverse). ----------------!
   elemental real(wp) function rwc_from_psi(psi, pi0, eps) result(rwc)
      real(wp), intent(in) :: psi, pi0, eps
      real(wp) :: psi_tlp, b, disc
      psi_tlp = pi0*eps / (pi0 + eps)
      if (psi >= psi_tlp) then                ! turgid (psi less negative than the TLP)
         !----- eps*R^2 - (psi+eps+pi0)*R + pi0 = 0; take the + root in [rwc_tlp, 1]. ------!
         b    = psi + eps + pi0
         disc = max(b*b - 4.0_wp*eps*pi0, 0.0_wp)
         rwc  = (b + sqrt(disc)) / (2.0_wp*eps)
      else                                    ! flaccid
         rwc  = pi0/psi
      end if
      rwc = max(rwc, rwc_floor)
   end function rwc_from_psi

   !----- Total tissue water [kg, per plant] at potential psi. -----------------------------!
   !      W = W_sat_sym*R + W_apoplast, with W_sat = water_sat*biomass, W_sat_sym =           !
   !      (1-af)*W_sat and W_apoplast = af*W_sat (a constant reservoir).                       !
   elemental real(wp) function water_content(psi, pi0, eps, af, water_sat, biomass) result(w)
      real(wp), intent(in) :: psi, pi0, eps, af, water_sat, biomass
      real(wp) :: w_sat, r
      w_sat = water_sat*biomass
      r     = rwc_from_psi(psi, pi0, eps)
      w     = (1.0_wp - af)*w_sat*r + af*w_sat
   end function water_content

   !----- Capacitance C = dW/dpsi [kg/MPa, per plant] at potential psi. ---------------------!
   elemental real(wp) function capacitance(psi, pi0, eps, af, water_sat, biomass) result(c)
      real(wp), intent(in) :: psi, pi0, eps, af, water_sat, biomass
      real(wp) :: psi_tlp, w_sat_sym, r, cr
      psi_tlp   = pi0*eps / (pi0 + eps)
      w_sat_sym = (1.0_wp - af)*water_sat*biomass
      r         = rwc_from_psi(psi, pi0, eps)
      if (psi >= psi_tlp) then
         cr = 1.0_wp / (eps - pi0/(r*r))      ! dR/dpsi, turgid  (-> 1/(eps+|pi0|) at R=1)
      else
         cr = -(r*r)/pi0                       ! dR/dpsi, flaccid (= |pi0|/psi^2)
      end if
      c = w_sat_sym*cr
   end function capacitance

   !----- Legacy-linear mapping: the constant capacitance per biomass [kg/kgC/MPa] that a       !
   !      pure-elastic (X16) reservoir would need to match this curve's full-turgor capacitance. !
   elemental real(wp) function pv_water_cap_from_traits(pi0, eps, af, water_sat) result(water_cap)
      real(wp), intent(in) :: pi0, eps, af, water_sat
      water_cap = (1.0_wp - af)*water_sat / (eps + abs(pi0))
   end function pv_water_cap_from_traits


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
      real(wp)    :: e_transp, soil_psi, rhz, grav, k_cond, bwood
      real(wp)    :: psi_l0, psi_w0, psi_l, psi_w, dw_l, dw_w
      real(wp)    :: h, t_rem, xl, xw, ml, mw, xl2, xw2, err, errl, errw, fac
      real(wp)    :: fm11, fm12, fm21, fm22, fps_l, fps_w   ! start-state frozen coeffs shared by full + half step
      integer(ik) :: nsub, i, nfix
      !------------------------------------------------------------------------------------!

      if (o%topology /= HYDRO_NODES_2) then
         call fatal_hydro('solve_plant_water: only HYDRO_NODES_2 is implemented (3-node is a follow-up).')
      end if

      !----- Boundary conditions / geometry, held constant over dt. ------------------------!
      e_transp = env%transp
      soil_psi = env%soil_psi
      rhz      = max(env%rhizo_cond, k_floor)          ! ground the network (avoid singular M)
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
            call expm_step(psi_l, psi_w, h, xl, xw)
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
            call apply_expm(psi_l, psi_w, h,        fm11, fm12, fm21, fm22, fps_l, fps_w, xl,  xw )  ! full step
            call apply_expm(psi_l, psi_w, 0.5_wp*h, fm11, fm12, fm21, fm22, fps_l, fps_w, ml,  mw )  ! first half
            call expm_step(ml,    mw,    0.5_wp*h,   xl2, xw2)                                       ! second half
            errl = abs(xl2 - xl) / (o%atol + o%rtol*abs(xl2))
            errw = abs(xw2 - xw) / (o%atol + o%rtol*abs(xw2))
            err  = max(errl, errw, 1.0e-12_wp)
            if (err <= 1.0_wp .or. h <= h_floor) then
               if (err > 1.0_wp) flux%converged = .false.          ! floor-forced accept below tolerance
               psi_l = xl2 ; psi_w = xw2 ; t_rem = t_rem - h ; nsub = nsub + 1_ik
               fac = safety * err**(-0.5_wp)                        ! p=1 => exponent -1/(p+1)
               h   = h * min(fmax, max(fmin, fac))
            else
               h   = h * max(fmin, safety*err**(-0.5_wp))          ! reject, shrink
            end if
            if (nsub >= o%max_substep .and. t_rem > tiny_num) then  ! cap hit before finishing interval
               flux%converged = .false. ; exit
            end if
         end do
      end if

      !----- Mass-closing boundary fluxes from the converged storage change. ---------------!
      dw_l = water_content(psi_l, p%leaf_pi0, p%leaf_eps, p%leaf_af, p%leaf_water_sat, env%bleaf) &
           - water_content(psi_l0,p%leaf_pi0, p%leaf_eps, p%leaf_af, p%leaf_water_sat, env%bleaf)
      dw_w = water_content(psi_w, p%wood_pi0, p%wood_eps, p%wood_af, p%wood_water_sat, bwood)      &
           - water_content(psi_w0,p%wood_pi0, p%wood_eps, p%wood_af, p%wood_water_sat, bwood)

      flux%sapflow     = dw_l/dt + e_transp
      flux%root_uptake = (dw_l + dw_w)/dt + e_transp
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
         cl   = max(capacitance(pl, p%leaf_pi0, p%leaf_eps, p%leaf_af, p%leaf_water_sat, env%bleaf), c_floor)
         cw   = max(capacitance(pw, p%wood_pi0, p%wood_eps, p%wood_af, p%wood_water_sat, bwood),     c_floor)
         keff = max(kirchhoff_edge(pw, pl, k_cond, p%wood_psi50, p%wood_kexp), k_floor)
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
      subroutine apply_expm(pl, pw, hs, m11, m12, m21, m22, ps_l, ps_w, plo, pwo)
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
      end subroutine apply_expm

      !----- One exact frozen-coefficient 2x2 matrix-exponential sub-step (freeze + apply).  !
      subroutine expm_step(pl, pw, hs, plo, pwo)
         real(wp), intent(in)  :: pl, pw, hs
         real(wp), intent(out) :: plo, pwo
         real(wp) :: m11, m12, m21, m22, ps_l, ps_w
         call freeze_coeffs(pl, pw, m11, m12, m21, m22, ps_l, ps_w)
         call apply_expm(pl, pw, hs, m11, m12, m21, m22, ps_l, ps_w, plo, pwo)
      end subroutine expm_step

   end subroutine solve_plant_water

   !---------------------------------------------------------------------------------------!
   ! plant_water_tendency -- the EXPLICIT 2-node hydraulics RHS dpsi/dt [MPa/s] at the current    !
   ! state, for the IMEX-ARK fast integrator (archive/MEDS_IMEX_ARK_DESIGN.md). It is dpsi/dt =    !
   ! M*psi + c with the SAME linear system solve_plant_water freezes each sub-step (the nested      !
   ! freeze_coeffs, :328-343); since apply_expm integrates exactly this system, (expm(psi,h)-psi)/h !
   ! -> this tendency as h -> 0. cond_max (:317-323) is mirrored here (it is a nested closure). The  !
   ! coefficients are frozen at the passed psi -- exactly what an ESDIRK stage evaluates. Pure.      !
   !---------------------------------------------------------------------------------------!
   pure subroutine plant_water_tendency(env, p, o, psi, dpsi_dt)
      type(hydro_env_t),    intent(in)  :: env
      type(hydro_params_t), intent(in)  :: p
      type(hydro_opts_t),   intent(in)  :: o
      real(wp),             intent(in)  :: psi(N_HYDRO)
      real(wp),             intent(out) :: dpsi_dt(N_HYDRO)

      real(wp) :: pl, pw, e_transp, soil_psi, rhz, grav, k_cond, bwood
      real(wp) :: cl, cw, keff, m11, m12, m21, m22, c1, c2

      dpsi_dt = 0.0_wp
      pl = psi(NODE_LEAF) ; pw = psi(NODE_WOOD)

      !----- Boundary conditions / geometry (matches solve_plant_water :246-252). -------------!
      e_transp = env%transp
      soil_psi = env%soil_psi
      rhz      = max(env%rhizo_cond, k_floor)
      bwood    = env%bsap + env%broot
      if (o%gravity_on) then ; grav = grav_head * env%height ; else ; grav = 0.0_wp ; end if
      !----- Maximum (plc=1) internal conductance (mirrors the nested cond_max :317-323). -----!
      if (o%cond_mode == HYDRO_COND_SEGMENT) then
         k_cond = p%wood_kmax * env%sap_area / max(env%height * p%vessel_curl, tiny_num)
      else
         k_cond = p%k_plant_max * env%leaf_area
      end if
      k_cond = max(k_cond, k_floor)

      !----- Linear system psi' = M psi + c, frozen at the current psi (mirrors freeze_coeffs). --!
      cl   = max(capacitance(pl, p%leaf_pi0, p%leaf_eps, p%leaf_af, p%leaf_water_sat, env%bleaf), c_floor)
      cw   = max(capacitance(pw, p%wood_pi0, p%wood_eps, p%wood_af, p%wood_water_sat, bwood),     c_floor)
      keff = max(kirchhoff_edge(pw, pl, k_cond, p%wood_psi50, p%wood_kexp), k_floor)
      m11 = -keff / cl              ; m12 =  keff / cl
      m21 =  keff / cw              ; m22 = -(keff + rhz) / cw
      c1  = (-keff * grav - e_transp) / cl
      c2  = ( keff * grav + rhz * soil_psi) / cw

      dpsi_dt(NODE_LEAF) = m11 * pl + m12 * pw + c1
      dpsi_dt(NODE_WOOD) = m21 * pl + m22 * pw + c2
   end subroutine plant_water_tendency

   !----- Catchable fatal (error stop), matching the MEDS convention. ----------------------!
   subroutine fatal_hydro(msg)
      character(len=*), intent(in) :: msg
      write(*,'(2a)') 'meds_plant_hydraulics: ', msg
      error stop 1
   end subroutine fatal_hydro

end module meds_plant_hydraulics
