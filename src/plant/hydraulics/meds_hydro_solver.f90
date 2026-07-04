!==========================================================================================!
! meds_hydro_solver -- advance the coupled plant water-transport network over one fast step.  !
!                                                                                          !
! The 2-node network (leaf L, lumped wood W) is the RC circuit                                 !
!     C_L dpsi_L/dt = K_LW (psi_W - psi_L - G) - E                                              !
!     C_W dpsi_W/dt = K_LW (psi_L - psi_W + G) + K_S (psi_soil - psi_W)                          !
! with G = grav_head*height the hydrostatic head, E the transpiration sink, K_S the given         !
! rhizosphere conductance, and K_LW the Kirchhoff-integrated whole-plant xylem conductance         !
! (meds_hydro_conductance). Capacitances C_L, C_W come from the nonlinear p-v curve                 !
! (meds_hydro_pv). Over a sub-step [t,t+h] the coefficients M,c are frozen and the linear system     !
! psi' = M psi + c is advanced by the EXACT 2x2 matrix exponential psi(t+h) = psi* + e^{Mh}(psi-psi*),!
! psi* = -M^{-1}c. Because M is similar to a symmetric matrix its eigenvalues are real and <=0, so     !
! the closed form is always valid and L-stable. Adaptive step-doubling controls the frozen-coeff       !
! (nonlinearity) error. Boundary fluxes are derived from the converged storage change so mass closes    !
! exactly: sapflow = dW_leaf/dt + E, root_uptake = (dW_leaf + dW_wood)/dt + E.                            !
!                                                                                          !
! Stateless: psi is passed in/out; nothing is stored. (3-node backward Euler is a planned follow-up.)     !
!==========================================================================================!
module meds_hydro_solver
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : grav_head, safe_exp, tiny_num
   use meds_hydro_types,      only : hydro_env_t, hydro_params_t, hydro_opts_t, hydro_flux_t,  &
                                     N_HYDRO, NODE_LEAF, NODE_WOOD,                             &
                                     HYDRO_NODES_2, HYDRO_COND_SEGMENT, HYDRO_SUBSTEP_FIXED
   use meds_hydro_pv,         only : capacitance, water_content
   use meds_hydro_conductance,only : kirchhoff_edge, plc_retained
   implicit none
   private

   public :: solve_plant_water

   real(wp),    parameter :: c_floor   = 1.0e-12_wp  !< capacitance floor (linearization only)
   real(wp),    parameter :: k_floor   = 1.0e-15_wp  !< conductance floor (keeps M finite)
   real(wp),    parameter :: sinhc_eps = 1.0e-8_wp   !< sinhc series threshold
   real(wp),    parameter :: safety    = 0.9_wp
   real(wp),    parameter :: fmax      = 4.0_wp
   real(wp),    parameter :: fmin      = 0.2_wp
   real(wp),    parameter :: h_floor   = 1.0e-6_wp   !< accept regardless below this sub-step [s]

contains

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
            call expm_step(psi_l, psi_w, h,          xl,  xw )      ! one full step
            call expm_step(psi_l, psi_w, 0.5_wp*h,   ml,  mw )      ! two half steps
            call expm_step(ml,    mw,    0.5_wp*h,   xl2, xw2)
            errl = abs(xl2 - xl) / (o%atol + o%rtol*abs(xl2))
            errw = abs(xw2 - xw) / (o%atol + o%rtol*abs(xw2))
            err  = max(errl, errw, 1.0e-12_wp)
            if (err <= 1.0_wp .or. h <= h_floor) then
               psi_l = xl2 ; psi_w = xw2 ; t_rem = t_rem - h ; nsub = nsub + 1_ik
               fac = safety * err**(-0.5_wp)                        ! p=1 => exponent -1/(p+1)
               h   = h * min(fmax, max(fmin, fac))
            else
               h   = h * max(fmin, safety*err**(-0.5_wp))          ! reject, shrink
            end if
            if (nsub >= o%max_substep) then
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

      !----- One exact frozen-coefficient 2x2 matrix-exponential sub-step. -----------------!
      subroutine expm_step(pl, pw, hs, plo, pwo)
         real(wp), intent(in)  :: pl, pw, hs
         real(wp), intent(out) :: plo, pwo
         real(wp) :: cl, cw, keff, m11, m12, m21, m22, c1, c2, detm, ps_l, ps_w
         real(wp) :: n11, n12, n21, n22, mu, dl, ep, em, ch, sh, e11, e12, e21, e22, ddl, ddw
         !----- Frozen coefficients at the sub-step start state (pl, pw). ------------------!
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
      end subroutine expm_step

   end subroutine solve_plant_water

   !----- Catchable fatal (error stop), matching the MEDS convention. ----------------------!
   subroutine fatal_hydro(msg)
      character(len=*), intent(in) :: msg
      write(*,'(2a)') 'meds_hydro_solver: ', msg
      error stop 1
   end subroutine fatal_hydro

end module meds_hydro_solver
