!==========================================================================================!
! meds_twostream_band -- THE unified single-band two-stream solver over a vertical stack of    !
! cohort layers. ONE routine serves every band (VIS, NIR, LW): shortwave bands carry a direct   !
! beam and no emission; the thermal band carries emission and no beam. The diffuse multiple-     !
! scattering operator is identical across bands -- only the source terms (beam scattering,        !
! thermal emission) and the incident/ground boundary conditions differ.                          !
!                                                                                          !
! Method: the ADDING (matrix-operator / layer-recursion) method -- an O(N), unconditionally      !
! stable direct solve of the block-tridiagonal two-stream system (the block-Thomas benchmarked    !
! in docs/dev_plans/radiative_transfer_design.md). Per cohort layer we form the analytic diffuse          !
! reflectance/transmittance (rdd, tdd) and the up/down source fluxes; a bottom-up sweep builds      !
! the reflectance & upward source of the growing stack, then a top-down sweep recovers the           !
! diffuse fluxes at every interface. The direct beam attenuates by Beer's law (crown fraction        !
! cai = 1) and injects its single-scattered fraction as a diffuse layer source.                      !
!                                                                                          !
! Cohorts are ordered BOTTOM (index 1) to TOP (index ncoh); a virtual sky interface sits above.      !
! Everything is absolute W m-2; total absorbed per layer is the beam + diffuse flux divergence, so     !
! energy is conserved by construction (in - out, telescoping).                                         !
!==========================================================================================!
module meds_twostream_band
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num, lnexp_max
   implicit none
   private

   public :: solve_band

contains

   !---------------------------------------------------------------------------------------!
   ! Solve ONE band. Inputs (per cohort, bottom->top): etai (effective total area index),     !
   ! omega/beta/beta0/kdir (blended optics), emission (blackbody sigma*T^4, 0 for shortwave).   !
   ! Scalars: incident beam & diffuse at canopy top, ground reflectance & emission, band flags.  !
   ! Outputs: absorbed(ncoh) [W/m2 per cohort, beam+diffuse], dn_ground/up_ground [W/m2] and     !
   ! albedo (upward fraction leaving the top).                                                    !
   !---------------------------------------------------------------------------------------!
   subroutine solve_band(ncoh, etai, omega, beta, beta0, kdir, emission, incid_beam,         &
                         incid_diff, grnd_refl, grnd_emiss, has_beam, has_emission,           &
                         absorbed, dn_ground, up_ground, albedo)
      integer(ik), intent(in)  :: ncoh
      real(wp),    intent(in)  :: etai(ncoh), omega(ncoh), beta(ncoh), beta0(ncoh), kdir(ncoh)
      real(wp),    intent(in)  :: emission(ncoh)
      real(wp),    intent(in)  :: incid_beam, incid_diff, grnd_refl, grnd_emiss
      logical,     intent(in)  :: has_beam, has_emission
      real(wp),    intent(out) :: absorbed(ncoh)
      real(wp),    intent(out) :: dn_ground, up_ground, albedo

      real(wp), dimension(ncoh)   :: rdd, tdd, sup, sdn          ! layer R/T + up/down sources
      real(wp), dimension(ncoh+1) :: down0                       ! direct-beam profile (interfaces)
      real(wp), dimension(0:ncoh) :: rbel, sbel                  ! reflectance & up-source of stack below
      real(wp), dimension(ncoh)   :: din, dbot, uup, utop        ! diffuse fluxes per layer
      real(wp) :: intercept, scat, emit_face, denom, incid_tot
      integer(ik) :: i

      !----- Direct-beam profile: Beer's law from the top down (cai = 1). -----------------!
      if (has_beam) then
         down0(ncoh+1) = incid_beam
         do i = ncoh, 1_ik, -1_ik
            down0(i) = down0(i+1) * exp(-min(kdir(i) * etai(i), lnexp_max))
         end do
      else
         down0(1:ncoh+1) = 0.0_wp                    ! no beam: zero the whole profile, incl. top face
      end if

      !----- Per-layer diffuse reflectance/transmittance and source fluxes. ----------------!
      do i = 1_ik, ncoh
         call layer_rt(omega(i), beta(i), etai(i), rdd(i), tdd(i))
         intercept = down0(i+1) - down0(i)                        ! beam intercepted by layer i (>=0)
         scat      = omega(i) * intercept                         ! beam single-scattered into diffuse
         emit_face = 0.0_wp
         if (has_emission) emit_face = (1.0_wp - rdd(i) - tdd(i)) * emission(i)
         sup(i) = beta0(i)          * scat + emit_face            ! up-source (top of layer)
         sdn(i) = (1.0_wp - beta0(i)) * scat + emit_face          ! down-source (bottom of layer)
      end do

      !----- Bottom-up sweep: reflectance & upward source of {ground + layers 1..i}. The       !
      !      ground's upward source is its thermal emission PLUS the reflection of the direct    !
      !      beam that reaches it (down0(1)); the diffuse reflection is handled by rbel(0).       !
      rbel(0) = grnd_refl
      sbel(0) = grnd_emiss + grnd_refl * down0(1)
      do i = 1_ik, ncoh
         denom   = 1.0_wp - rbel(i-1) * rdd(i)
         rbel(i) = rdd(i) + tdd(i) * tdd(i) * rbel(i-1) / denom
         sbel(i) = sup(i) + tdd(i) * (sbel(i-1) + rbel(i-1) * sdn(i)) / denom
      end do

      !----- Top-down sweep: recover diffuse fluxes at every interface. --------------------!
      din(ncoh) = incid_diff
      do i = ncoh, 1_ik, -1_ik
         denom   = 1.0_wp - rdd(i) * rbel(i-1)
         dbot(i) = (tdd(i) * din(i) + sdn(i) + rdd(i) * sbel(i-1)) / denom
         uup(i)  = sbel(i-1) + rbel(i-1) * dbot(i)                ! up flux at bottom of layer i
         utop(i) = rdd(i) * din(i) + tdd(i) * uup(i) + sup(i)     ! up flux at top of layer i
         if (i > 1_ik) din(i-1) = dbot(i)
      end do

      !----- Absorbed per layer = beam divergence + diffuse divergence (in - out). ---------!
      do i = 1_ik, ncoh
         absorbed(i) = (down0(i+1) - down0(i)) + (din(i) + uup(i)) - (utop(i) + dbot(i))
      end do

      !----- Ground and top-of-canopy diagnostics. ----------------------------------------!
      dn_ground = dbot(1) + down0(1)                              ! diffuse + beam reaching the ground
      up_ground = uup(1)                                          ! upwelling leaving the ground
      albedo    = 0.0_wp
      incid_tot = incid_beam + incid_diff
      if (incid_tot > tiny_num) albedo = utop(ncoh) / incid_tot
   end subroutine solve_band

   !---------------------------------------------------------------------------------------!
   ! Analytic two-stream diffuse reflectance (rdd) and transmittance (tdd) of a homogeneous    !
   ! layer of effective area index `etai`, single-scatter albedo `omega`, backscatter `beta`.   !
   ! Uses the classical Stokes formulas with rinf = (att-m)/sigb, e = exp(-m*etai).             !
   !---------------------------------------------------------------------------------------!
   pure subroutine layer_rt(omega, beta, etai, rdd, tdd)
      real(wp), intent(in)  :: omega, beta, etai
      real(wp), intent(out) :: rdd, tdd
      real(wp) :: w, sigb, sigf, att, m, rinf, e, e2, denom
      w = min(omega, 1.0_wp - 1.0e-6_wp)          ! guard the conservative-scattering singularity
      sigb = w * beta                              ! backward diffuse scattering
      sigf = w * (1.0_wp - beta)                   ! forward diffuse scattering
      att  = 1.0_wp - sigf
      if (sigb < tiny_num) then                    ! no backscatter -> no diffuse reflectance
         rdd = 0.0_wp
         tdd = exp(-min(att * etai, lnexp_max))
         return
      end if
      m    = sqrt(max(att * att - sigb * sigb, 0.0_wp))
      rinf = (att - m) / sigb
      e    = exp(-min(m * etai, lnexp_max))
      e2   = e * e
      denom = 1.0_wp - rinf * rinf * e2
      tdd  = (1.0_wp - rinf * rinf) * e / denom
      rdd  = rinf * (1.0_wp - e2) / denom
   end subroutine layer_rt

end module meds_twostream_band
