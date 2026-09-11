!==========================================================================================!
! meds_column_state_ops -- the algebra and bookkeeping helpers on the fast loop's column state     !
! (column_state_t / column_tend_t / stage_bflux_t / column_bflux_t), owned by the types they act   !
! on rather than by one integrator. Both schemes (meds_fast_ark, meds_fast_rk45) and the test-only  !
! oracle import from here, so RK45 no longer depends on the ARK module for non-ARK code.            !
!                                                                                          !
!   * state combinators  -- state_init / state_axpy / state_accum / state_extrap / state_sub /       !
!                           state_err_diff / zero_like  (every prognostic field, every time)        !
!   * boundary-flux ledger accumulators -- bflux_zero / bflux_add / bflux_bweight                   !
!   * stage-domain clamps -- clamp_cas / clamp_theta / clamp_soil_energy                             !
!   * whole-column store totals for the ledgers -- soil_water_store / soil_energy_store /            !
!                           plant_water_store / canopy_film_store  (one loop order, so the two       !
!                           schemes' ledgers are the same arithmetic)                                !
!   * post-march commits shared by both schemes -- clamp_canopy_film / deposit_condensate /          !
!                           unpack_column_state / diagnose_soil_temps                                !
!                                                                                          !
! ADDING A STATE FIELD: every combinator here must carry it. There is no compile-time check, so     !
! grep this file for an existing field (e.g. wood_surf_water) and mirror each occurrence.           !
!==========================================================================================!
module meds_column_state_ops
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, tiny_num
   use meds_therm_lib,        only : internal_energy_to_temp, temp_to_internal_energy, cas_temp_of_enthalpy, cas_enthalpy_of_temp
   use meds_fast_types,       only : column_state_t, column_tend_t, column_frozen_t,                     &
                                     stage_bflux_t, column_bflux_t, process_mask_t
   use meds_soil_types, only : energy_forcing_t
   use meds_fast_types, only : patch_biophys_t
   implicit none
   private

   public :: state_init, state_axpy, state_accum, state_extrap, state_sub, state_err_diff, zero_like
   public :: bflux_zero, bflux_add, bflux_bweight
   public :: assemble_soil_energy_forcing, apply_process_mask
   public :: clamp_cas, clamp_theta, clamp_soil_energy
   public :: soil_water_store, soil_energy_store, plant_water_store, canopy_film_store
   public :: clamp_canopy_film, deposit_condensate, unpack_column_state, diagnose_soil_temps

contains

   !----- copy the prognostic state (used to seed the RK combination). --------------------!
   pure subroutine state_init(y, n, nsl, y_stage)
      type(column_state_t), intent(in)  :: y
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: y_stage
      y_stage%cas_enthalpy = y%cas_enthalpy ; y_stage%cas_shv = y%cas_shv ; y_stage%cas_co2 = y%cas_co2
      y_stage%soil_energy  = y%soil_energy  ; y_stage%theta   = y%theta
      y_stage%w_surface    = y%w_surface    ; y_stage%w_surface_enth = y%w_surface_enth
      allocate(y_stage%leaf_water_mass(n), y_stage%wood_water_mass(n))
      y_stage%leaf_water_mass(1:n) = y%leaf_water_mass(1:n)
      y_stage%wood_water_mass(1:n) = y%wood_water_mass(1:n)
      allocate(y_stage%leaf_surf_water(n), y_stage%wood_surf_water(n))
      y_stage%leaf_surf_water(1:n) = y%leaf_surf_water(1:n)
      y_stage%wood_surf_water(1:n) = y%wood_surf_water(1:n)
   end subroutine state_init

   !----- y_stage = y + a*k  (state + a * tendency) -- the single-term combinator classical RK4's mid-  !
   !      point/endpoint stages use. Cash-Karp's later stages need a MULTI-term combination (each    !
   !      reads several prior k's), for which state_init + repeated state_accum is the pattern; both  !
   !      live here together as the ONE set of generic column_state_t/column_tend_t combinators        !
   !      every explicit fast-loop integrator (the RK4 oracle, RK45) builds its stages from. -----------!
   pure subroutine state_axpy(y, a, k, n, nsl, y_stage)
      type(column_state_t), intent(in)  :: y
      real(wp),             intent(in)  :: a
      type(column_tend_t),  intent(in)  :: k
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: y_stage
      integer(ik) :: j, i
      y_stage%cas_enthalpy = y%cas_enthalpy + a * k%d_cas_enthalpy
      y_stage%cas_shv      = y%cas_shv      + a * k%d_cas_shv
      y_stage%cas_co2      = y%cas_co2      + a * k%d_cas_co2
      y_stage%soil_energy  = y%soil_energy
      y_stage%theta        = y%theta
      !----- pond PASSED THROUGH (no stage tendency yet -- #93 Phase 1 gives it one). ----------!
      y_stage%w_surface    = y%w_surface ; y_stage%w_surface_enth = y%w_surface_enth
      do j = 1_ik, nsl
         y_stage%soil_energy(j) = y%soil_energy(j) + a * k%dedt(j)
         y_stage%theta(j)       = y%theta(j)       + a * k%dtheta_dt(j)
      end do
      allocate(y_stage%leaf_water_mass(n), y_stage%wood_water_mass(n))
      allocate(y_stage%leaf_surf_water(n), y_stage%wood_surf_water(n))
      do i = 1_ik, n
         y_stage%leaf_water_mass(i) = y%leaf_water_mass(i) + a * k%d_leaf_water_mass(i)
         y_stage%wood_water_mass(i) = y%wood_water_mass(i) + a * k%d_wood_water_mass(i)
         y_stage%leaf_surf_water(i) = y%leaf_surf_water(i) + a * k%d_leaf_surf_water(i)
         y_stage%wood_surf_water(i) = y%wood_surf_water(i) + a * k%d_wood_surf_water(i)
      end do
   end subroutine state_axpy

   !----- y_stage += a*k  (accumulate a weighted tendency into a state). -----------------------!
   pure subroutine state_accum(y_stage, a, k, n, nsl)
      type(column_state_t), intent(inout) :: y_stage
      real(wp),             intent(in)    :: a
      type(column_tend_t),  intent(in)    :: k
      integer(ik),          intent(in)    :: n, nsl
      integer(ik) :: j, i
      y_stage%cas_enthalpy = y_stage%cas_enthalpy + a * k%d_cas_enthalpy
      y_stage%cas_shv      = y_stage%cas_shv      + a * k%d_cas_shv
      y_stage%cas_co2      = y_stage%cas_co2      + a * k%d_cas_co2
      do j = 1_ik, nsl
         y_stage%soil_energy(j) = y_stage%soil_energy(j) + a * k%dedt(j)
         y_stage%theta(j)       = y_stage%theta(j)       + a * k%dtheta_dt(j)
      end do
      do i = 1_ik, n
         y_stage%leaf_water_mass(i) = y_stage%leaf_water_mass(i) + a * k%d_leaf_water_mass(i)
         y_stage%wood_water_mass(i) = y_stage%wood_water_mass(i) + a * k%d_wood_water_mass(i)
         y_stage%leaf_surf_water(i) = y_stage%leaf_surf_water(i) + a * k%d_leaf_surf_water(i)
         y_stage%wood_surf_water(i) = y_stage%wood_surf_water(i) + a * k%d_wood_surf_water(i)
      end do
      !----- w_surface / w_surface_enth are EXCLUDED, not forgotten. `column_tend_t` carries no      !
      !      d_w_surface: the pond has no stage RHS, so it is passed through the stages untouched    !
      !      and committed from the scratch hydrology solve. Left ALONE here (this is intent(inout), !
      !      so "leave alone" is the correct action, and zeroing would destroy it). Asserted in      !
      !      test_state_combinators, because otherwise this exclusion and an omission look the same. !
   end subroutine state_accum

   !----- ledger helpers: b-weight two stage RATE structs into accumulated AMOUNTS over dt (weights   !
   !      b^I = (1-gamma, gamma) times dt); zero an accumulator; add one substep's amounts. -----------!
   pure subroutine bflux_bweight(acc, s2, s3, dt, gam)
      type(column_bflux_t), intent(out) :: acc
      type(stage_bflux_t),  intent(in)  :: s2, s3
      real(wp),             intent(in)  :: dt, gam
      real(wp) :: b2, b3
      b2 = (1.0_wp - gam) * dt ; b3 = gam * dt
      acc%cas_enth_in   = b2*s2%cas_enth_in   + b3*s3%cas_enth_in
      acc%cas_enth_out  = b2*s2%cas_enth_out  + b3*s3%cas_enth_out
      acc%cas_vap_in    = b2*s2%cas_vap_in    + b3*s3%cas_vap_in
      acc%cas_vap_out   = b2*s2%cas_vap_out   + b3*s3%cas_vap_out
      acc%cas_co2_in    = b2*s2%cas_co2_in    + b3*s3%cas_co2_in
      acc%cas_co2_out   = b2*s2%cas_co2_out   + b3*s3%cas_co2_out
      acc%soil_enth_in  = b2*s2%soil_enth_in  + b3*s3%soil_enth_in
      acc%soil_enth_out = b2*s2%soil_enth_out + b3*s3%soil_enth_out
      acc%soil_wat_in   = b2*s2%soil_wat_in   + b3*s3%soil_wat_in
      acc%soil_wat_out  = b2*s2%soil_wat_out  + b3*s3%soil_wat_out
      acc%whole_enth_in = b2*s2%whole_enth_in + b3*s3%whole_enth_in
      acc%whole_enth_out= b2*s2%whole_enth_out+ b3*s3%whole_enth_out
      acc%whole_wat_in  = b2*s2%whole_wat_in  + b3*s3%whole_wat_in
      acc%whole_wat_out = b2*s2%whole_wat_out + b3*s3%whole_wat_out
      acc%whole_cond    = b2*s2%whole_cond    + b3*s3%whole_cond
      acc%whole_cond_enth = b2*s2%whole_cond_enth + b3*s3%whole_cond_enth
      acc%atm_heat_out  = b2*s2%atm_heat_out  + b3*s3%atm_heat_out
      acc%atm_vap_out   = b2*s2%atm_vap_out   + b3*s3%atm_vap_out
   end subroutine bflux_bweight

   pure subroutine bflux_zero(acc, n)
      type(column_bflux_t), intent(out) :: acc
      integer(ik), optional, intent(in) :: n   !< allocate + zero the per-cohort tissue integrals
      !----- intent(out) already default-initialises every scalar to zero and deallocates the tissue  !
      !      integrals on entry (F2018 8.5.10). An explicit `acc = column_bflux_t()` here is redundant !
      !      and nvfortran 25.11 rejects it in this module ("Empty structure constructor", F-0155). ---!
      if (present(n)) then
         allocate(acc%tissue_leaf_int(n), acc%tissue_wood_int(n))
         acc%tissue_leaf_int = 0.0_wp ; acc%tissue_wood_int = 0.0_wp
      end if
   end subroutine bflux_zero

   pure subroutine bflux_add(acc, s)
      type(column_bflux_t), intent(inout) :: acc
      type(column_bflux_t), intent(in)    :: s
      acc%cas_enth_in   = acc%cas_enth_in   + s%cas_enth_in
      acc%cas_enth_out  = acc%cas_enth_out  + s%cas_enth_out
      acc%cas_vap_in    = acc%cas_vap_in    + s%cas_vap_in
      acc%cas_vap_out   = acc%cas_vap_out   + s%cas_vap_out
      acc%cas_co2_in    = acc%cas_co2_in    + s%cas_co2_in
      acc%cas_co2_out   = acc%cas_co2_out   + s%cas_co2_out
      acc%soil_enth_in  = acc%soil_enth_in  + s%soil_enth_in
      acc%soil_enth_out = acc%soil_enth_out + s%soil_enth_out
      acc%soil_wat_in   = acc%soil_wat_in   + s%soil_wat_in
      acc%soil_wat_out  = acc%soil_wat_out  + s%soil_wat_out
      acc%whole_enth_in = acc%whole_enth_in + s%whole_enth_in
      acc%whole_enth_out= acc%whole_enth_out+ s%whole_enth_out
      acc%whole_wat_in  = acc%whole_wat_in  + s%whole_wat_in
      acc%whole_wat_out = acc%whole_wat_out + s%whole_wat_out
      acc%whole_cond    = acc%whole_cond    + s%whole_cond
      acc%whole_cond_enth = acc%whole_cond_enth + s%whole_cond_enth
      acc%atm_heat_out  = acc%atm_heat_out  + s%atm_heat_out
      acc%atm_vap_out   = acc%atm_vap_out   + s%atm_vap_out
      !----- Only ACCEPTED sub-steps reach here, so the tissue integrals accumulate over exactly the  !
      !      accepted march -- the same set of sub-steps every other amount above is summed over. -----!
      if (allocated(acc%tissue_leaf_int) .and. allocated(s%tissue_leaf_int)) then
         acc%tissue_leaf_int = acc%tissue_leaf_int + s%tissue_leaf_int
         acc%tissue_wood_int = acc%tissue_wood_int + s%tissue_wood_int
      end if
   end subroutine bflux_add

   !----- out = (1-b)*y + b*Y2  (the ARS stage-3 extrapolation base). --------------------------!
   pure subroutine state_extrap(y, b, Y2, n, nsl, out)
      type(column_state_t), intent(in)  :: y, Y2
      real(wp),             intent(in)  :: b
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: out
      real(wp)    :: a
      integer(ik) :: k, i
      a = 1.0_wp - b
      out%cas_enthalpy = a*y%cas_enthalpy + b*Y2%cas_enthalpy
      out%cas_shv      = a*y%cas_shv      + b*Y2%cas_shv
      out%cas_co2      = a*y%cas_co2      + b*Y2%cas_co2
      out%soil_energy = y%soil_energy ; out%theta = y%theta
      !----- == y%w_surface (the pond is frozen in the stages, like the mass stores). ----------!
      out%w_surface      = a*y%w_surface      + b*Y2%w_surface
      out%w_surface_enth = a*y%w_surface_enth + b*Y2%w_surface_enth
      do k = 1_ik, nsl
         out%soil_energy(k) = a*y%soil_energy(k) + b*Y2%soil_energy(k)
         out%theta(k)       = a*y%theta(k)       + b*Y2%theta(k)
      end do
      allocate(out%leaf_water_mass(n), out%wood_water_mass(n))
      allocate(out%leaf_surf_water(n), out%wood_surf_water(n))
      do i = 1_ik, n
         !----- == y%*_water_mass (mass is frozen in the stages, like psi was). ----------------!
         out%leaf_water_mass(i) = a*y%leaf_water_mass(i) + b*Y2%leaf_water_mass(i)
         out%wood_water_mass(i) = a*y%wood_water_mass(i) + b*Y2%wood_water_mass(i)
         !----- == y%*_surf_water (surface water is ALSO frozen/split out of the stages, sec 3.4/P2c). !
         out%leaf_surf_water(i) = a*y%leaf_surf_water(i) + b*Y2%leaf_surf_water(i)
         out%wood_surf_water(i) = a*y%wood_surf_water(i) + b*Y2%wood_surf_water(i)
      end do
   end subroutine state_extrap

   !----- clamp the extrapolated CAS enthalpy + humidity into a wide PHYSICAL range so a BETA=2.414   !
   !      overshoot cannot drive cas_temp_of_enthalpy to a wild T where qsat(T) overflows to NaN. Only  !
   !      active on a pathological overshoot (then the step is rejected); an in-range base3 is untouched.!
   !      `nfire` (optional) counts this call as an activation when the clamp actually moved the      !
   !      state -- see column_budget_t's clamp_* fields for why activations are tracked at all.       !
   pure subroutine clamp_cas(s, nfire)
      type(column_state_t), intent(inout) :: s
      integer(ik), optional, intent(inout) :: nfire
      real(wp) :: t, shv_c, enth_in
      real(wp), parameter :: T_LO = 180.0_wp, T_HI = 350.0_wp, SHV_LO = 1.0e-8_wp, SHV_HI = 0.06_wp
      enth_in = s%cas_enthalpy
      shv_c = min(max(s%cas_shv, SHV_LO), SHV_HI)
      t     = cas_temp_of_enthalpy(s%cas_enthalpy, shv_c)
      t     = min(max(t, T_LO), T_HI)
      s%cas_shv      = shv_c
      s%cas_enthalpy = cas_enthalpy_of_temp(t, shv_c)
      !----- an in-range state round-trips through cas_temp_of_enthalpy/cas_enthalpy_of_temp, so    !
      !      compare against the INPUT rather than testing the bounds -- that also catches a         !
      !      shv-only clamp, which moves enthalpy through the humidity term. ------------------------!
      if (present(nfire)) then
         if (s%cas_enthalpy /= enth_in) nfire = nfire + 1_ik
      end if
   end subroutine clamp_cas

   !----- clamp the extrapolated theta into [theta_res, theta_sat] (van Genuchten domain).       !
   !      dmass (optional) accumulates |water| moved, in kg/m2 of GROUND -- the mass this clamp   !
   !      creates or destroys with no ledger entry. -----------------------------------------!
   pure subroutine clamp_theta(s, frozen, nsl, nfire, dmass)
      type(column_state_t),  intent(inout) :: s
      type(column_frozen_t), intent(in)    :: frozen
      integer(ik),           intent(in)    :: nsl
      integer(ik), optional, intent(inout) :: nfire
      real(wp),    optional, intent(inout) :: dmass
      integer(ik) :: k
      real(wp)    :: th_in
      do k = 1_ik, nsl
         th_in      = s%theta(k)
         s%theta(k) = min(max(s%theta(k), frozen%params%soil%theta_res(k)), frozen%params%soil%theta_sat(k))
         if (s%theta(k) /= th_in) then
            if (present(nfire)) nfire = nfire + 1_ik
            if (present(dmass)) dmass = dmass + abs(s%theta(k) - th_in) * frozen%params%soil%dz(k) * rho_h2o
         end if
      end do
   end subroutine clamp_theta

   !----- clamp each soil layer's internal energy into a wide PHYSICAL temperature range, the      !
   !      soil-column analogue of clamp_cas above: an explicit-stage overshoot in soil_energy         !
   !      otherwise diagnoses (internal_energy_to_temp) a wild soil temperature that overflows ground_evaporation's  !
   !      fractional pow() (a negative base to a non-integer exponent is a domain error, not just an     !
   !      overflow) or qsat. Reconstructs soil_energy (temp_to_internal_energy) at the CLAMPED temperature and the    !
   !      SAME liquid fraction the (possibly wild) input diagnosed -- an in-range input is untouched, and  !
   !      the step is rejected normally by the adaptive controller when this bites. Call AFTER clamp_theta  !
   !      (uses the already-clamped theta for the water-mass term of the phase-change inverter). ----------!
   !      denergy (optional) accumulates |energy| moved, in J/m2 of GROUND -- the energy this clamp    !
   !      creates or destroys with no ledger entry. -------------------------------------------------!
   pure subroutine clamp_soil_energy(s, frozen, nsl, nfire, denergy)
      type(column_state_t),  intent(inout) :: s
      type(column_frozen_t), intent(in)    :: frozen
      integer(ik),           intent(in)    :: nsl
      integer(ik), optional, intent(inout) :: nfire
      real(wp),    optional, intent(inout) :: denergy
      real(wp), parameter :: T_LO = 180.0_wp, T_HI = 350.0_wp
      real(wp) :: temp, fliq, wmass, e_in
      integer(ik) :: k
      do k = 1_ik, nsl
         wmass = s%theta(k) * rho_h2o
         e_in  = s%soil_energy(k)
         call internal_energy_to_temp(s%soil_energy(k), wmass, frozen%params%therm%soil_dry_heat_capacity(k), temp, fliq)
         fliq  = min(max(fliq, 0.0_wp), 1.0_wp)
         temp  = min(max(temp, T_LO), T_HI)
         s%soil_energy(k) = temp_to_internal_energy(frozen%params%therm%soil_dry_heat_capacity(k), wmass, temp, fliq)
         !----- compare against the INPUT, not the T bounds: the internal_energy_to_temp/temp_to_internal_energy round trip   !
         !      is the identity only for an in-range state, so this also catches a clamp that bit       !
         !      through the liquid-fraction bound rather than the temperature bound. -------------------!
         if (s%soil_energy(k) /= e_in) then
            if (present(nfire))   nfire   = nfire   + 1_ik
            if (present(denergy)) denergy = denergy + abs(s%soil_energy(k) - e_in) * frozen%params%soil%dz(k)
         end if
      end do
   end subroutine clamp_soil_energy

   !----- err = (Y3 - base3) - (Y2 - y)  (the embedded 2nd-1st order difference); mass zeroed     !
   !      (like psi before it -- mass is frozen/operator-split through the ESDIRK stages, so       !
   !      Y3%*_water_mass == base3%*_water_mass == Y2%*_water_mass == y%*_water_mass exactly). -----!
   pure subroutine state_err_diff(Y3, base3, Y2, y, n, nsl, err)
      type(column_state_t), intent(in)  :: Y3, base3, Y2, y
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: err
      integer(ik) :: k
      err%cas_enthalpy = (Y3%cas_enthalpy - base3%cas_enthalpy) - (Y2%cas_enthalpy - y%cas_enthalpy)
      err%cas_shv      = (Y3%cas_shv      - base3%cas_shv)      - (Y2%cas_shv      - y%cas_shv)
      err%cas_co2      = (Y3%cas_co2      - base3%cas_co2)      - (Y2%cas_co2      - y%cas_co2)
      err%soil_energy = 0.0_wp ; err%theta = 0.0_wp
      do k = 1_ik, nsl
         err%soil_energy(k) = (Y3%soil_energy(k) - base3%soil_energy(k)) - (Y2%soil_energy(k) - y%soil_energy(k))
         err%theta(k)       = (Y3%theta(k)       - base3%theta(k))       - (Y2%theta(k)       - y%theta(k))
      end do
      allocate(err%leaf_water_mass(n), err%wood_water_mass(n))
      err%leaf_water_mass(1:n) = 0.0_wp
      err%wood_water_mass(1:n) = 0.0_wp
      !----- surface water is ALSO split out of the embedded estimate (sec 3.4/P2c), like mass above. --!
      allocate(err%leaf_surf_water(n), err%wood_surf_water(n))
      err%leaf_surf_water(1:n) = 0.0_wp
      err%wood_surf_water(1:n) = 0.0_wp
      !----- The POND is excluded too, and is written EXPLICITLY rather than left to the type's      !
      !      default initializer. The value is the same; what changes is that a reader can tell       !
      !      exclusion from omission. Every other field of column_state_t is assigned above, so an     !
      !      unassigned one here is indistinguishable from a field somebody forgot -- which is the      !
      !      whole silent-omission hazard (structure plan §10.3): the compiler cannot flag it, the      !
      !      conservation ledgers cannot see it, and the step controller just reads a smaller error.    !
      err%w_surface      = 0.0_wp
      err%w_surface_enth = 0.0_wp
   end subroutine state_err_diff

   !----- out = a - b  (state difference; used to form the low-order embedded solution). --------!
   pure subroutine state_sub(a, b, n, nsl, out)
      type(column_state_t), intent(in)  :: a, b
      integer(ik),          intent(in)  :: n, nsl
      type(column_state_t), intent(out) :: out
      integer(ik) :: k
      out%cas_enthalpy = a%cas_enthalpy - b%cas_enthalpy
      out%cas_shv      = a%cas_shv      - b%cas_shv
      out%cas_co2      = a%cas_co2      - b%cas_co2
      out%soil_energy = a%soil_energy ; out%theta = a%theta
      out%w_surface      = a%w_surface      - b%w_surface
      out%w_surface_enth = a%w_surface_enth - b%w_surface_enth
      do k = 1_ik, nsl
         out%soil_energy(k) = a%soil_energy(k) - b%soil_energy(k)
         out%theta(k)       = a%theta(k)       - b%theta(k)
      end do
      allocate(out%leaf_water_mass(n), out%wood_water_mass(n))
      out%leaf_water_mass(1:n) = a%leaf_water_mass(1:n) - b%leaf_water_mass(1:n)
      out%wood_water_mass(1:n) = a%wood_water_mass(1:n) - b%wood_water_mass(1:n)
      allocate(out%leaf_surf_water(n), out%wood_surf_water(n))
      out%leaf_surf_water(1:n) = a%leaf_surf_water(1:n) - b%leaf_surf_water(1:n)
      out%wood_surf_water(1:n) = a%wood_surf_water(1:n) - b%wood_surf_water(1:n)
   end subroutine state_sub

   !----- a column_state_t of the SAME shape as `ref`, every field zeroed -- lets state_wrms_grouped   !
   !      (which takes two STATES to difference) read a already-a-difference y_err directly, without    !
   !      a bespoke "WRMS of one state" variant. Trivial and allocation-only; not a hot path (once      !
   !      per accept/reject trial, not per stage). ------------------------------------------------------!
   pure function zero_like(ref, n, nsl) result(z)
      type(column_state_t), intent(in) :: ref
      integer(ik),           intent(in) :: n, nsl
      type(column_state_t) :: z
      z%cas_enthalpy = 0.0_wp ; z%cas_shv = 0.0_wp ; z%cas_co2 = 0.0_wp
      z%soil_energy = 0.0_wp ; z%theta = 0.0_wp
      allocate(z%leaf_water_mass(n), z%wood_water_mass(n), z%leaf_surf_water(n), z%wood_surf_water(n))
      z%leaf_water_mass = 0.0_wp ; z%wood_water_mass = 0.0_wp
      z%leaf_surf_water = 0.0_wp ; z%wood_surf_water = 0.0_wp   ! every field, so the norm may read any of them
      z%w_surface = 0.0_wp ; z%w_surface_enth = 0.0_wp
   end function zero_like

   !---------------------------------------------------------------------------------------!
   ! Whole-column STORE TOTALS for the conservation ledgers [per m2 ground]. Written as the         !
   ! ascending loops the two schemes used to carry inline, so a ledger built from these is the      !
   ! same arithmetic on both.                                                                       !
   !---------------------------------------------------------------------------------------!
   pure function soil_water_store(theta, dz, nsl) result(w)     ! [kg/m2]
      real(wp),    intent(in) :: theta(:), dz(:)
      integer(ik), intent(in) :: nsl
      real(wp) :: w
      integer(ik) :: k
      w = 0.0_wp
      do k = 1_ik, nsl
         w = w + theta(k) * dz(k) * rho_h2o
      end do
   end function soil_water_store

   pure function soil_energy_store(soil_energy, dz, nsl) result(e)   ! [J/m2] from [J/m3]
      real(wp),    intent(in) :: soil_energy(:), dz(:)
      integer(ik), intent(in) :: nsl
      real(wp) :: e
      integer(ik) :: k
      e = 0.0_wp
      do k = 1_ik, nsl
         e = e + soil_energy(k) * dz(k)
      end do
   end function soil_energy_store

   pure function plant_water_store(nplant, leaf_mass, wood_mass, n) result(w)   ! [kg/m2] from [kg/plant]
      real(wp),    intent(in) :: nplant(:), leaf_mass(:), wood_mass(:)
      integer(ik), intent(in) :: n
      real(wp) :: w
      w = sum(nplant(1:n) * (leaf_mass(1:n) + wood_mass(1:n)))
   end function plant_water_store

   pure function canopy_film_store(leaf_surf, wood_surf, n) result(w)   ! [kg/m2 ground], already ground-referenced
      real(wp),    intent(in) :: leaf_surf(:), wood_surf(:)
      integer(ik), intent(in) :: n
      real(wp) :: w
      w = sum(leaf_surf(1:n) + wood_surf(1:n))
   end function canopy_film_store

   !---------------------------------------------------------------------------------------!
   ! Canopy-film capacity clamp on the COMMITTED state: cap each cohort's leaf/wood film at        !
   ! dewmx*LAI / dewmx*WAI and floor it at 0; the clipped excess (overflow) and the floored          !
   ! shortfall (deficit) are returned so the caller books them in the whole-column ledger.          !
   !---------------------------------------------------------------------------------------!
   pure subroutine clamp_canopy_film(y_out, lai, wai, dewmx, n, surf_overflow, surf_deficit)
      type(column_state_t), intent(inout) :: y_out
      real(wp),             intent(in)    :: lai(:), wai(:), dewmx
      integer(ik),          intent(in)    :: n
      real(wp),             intent(inout) :: surf_overflow, surf_deficit   !< [kg/m2] accumulated
      real(wp)    :: leaf_cap_i, wood_cap_i
      integer(ik) :: i
      do i = 1_ik, n
         leaf_cap_i = dewmx * lai(i) ; wood_cap_i = dewmx * wai(i)
         surf_overflow = surf_overflow + max(0.0_wp, y_out%leaf_surf_water(i) - leaf_cap_i)          &
                                        + max(0.0_wp, y_out%wood_surf_water(i) - wood_cap_i)
         surf_deficit  = surf_deficit  + max(0.0_wp, -y_out%leaf_surf_water(i))                      &
                                        + max(0.0_wp, -y_out%wood_surf_water(i))
         y_out%leaf_surf_water(i) = min(max(y_out%leaf_surf_water(i), 0.0_wp), leaf_cap_i)
         y_out%wood_surf_water(i) = min(max(y_out%wood_surf_water(i), 0.0_wp), wood_cap_i)
      end do
   end subroutine clamp_canopy_film

   !---------------------------------------------------------------------------------------!
   ! Deposit condensed canopy-air vapour into soil layer 1 as paired mass [kg/m2] + liquid         !
   ! enthalpy [J/m2] (an internal CAS -> soil transfer; the ledger needs no boundary term).         !
   !---------------------------------------------------------------------------------------!
   pure subroutine deposit_condensate(y_out, dz1, mass, enth)
      type(column_state_t), intent(inout) :: y_out
      real(wp),             intent(in)    :: dz1, mass, enth
      y_out%theta(1)       = y_out%theta(1) + mass / (rho_h2o * dz1)
      y_out%soil_energy(1) = y_out%soil_energy(1) + enth / dz1
   end subroutine deposit_condensate

   !---------------------------------------------------------------------------------------!
   ! Unpack the committed column state into the patch biophysics record: CAS twins (+ diagnosed    !
   ! temperature), soil energy and moisture, plant internal water, canopy films. The ponding       !
   ! store is NOT unpacked here -- both schemes commit it separately under the soil-water mask.    !
   !---------------------------------------------------------------------------------------!
   pure subroutine unpack_column_state(y_out, n, nsl, biophys)
      type(column_state_t),  intent(in)    :: y_out
      integer(ik),           intent(in)    :: n, nsl
      type(patch_biophys_t), intent(inout) :: biophys
      biophys%cas%can_enthalpy = y_out%cas_enthalpy ; biophys%cas%can_shv = y_out%cas_shv ; biophys%cas%can_co2 = y_out%cas_co2
      biophys%cas%can_temp = cas_temp_of_enthalpy(y_out%cas_enthalpy, y_out%cas_shv)
      biophys%soil_e%soil_energy(1:nsl) = y_out%soil_energy(1:nsl)
      biophys%soil_w%theta(1:nsl)       = y_out%theta(1:nsl)
      biophys%leaf_water_mass(1:n) = y_out%leaf_water_mass(1:n)
      biophys%wood_water_mass(1:n) = y_out%wood_water_mass(1:n)
      biophys%leaf_surf_water(1:n) = y_out%leaf_surf_water(1:n)
      biophys%wood_surf_water(1:n) = y_out%wood_surf_water(1:n)
      !----- w_surface / w_surface_enth are EXCLUDED, not forgotten. The pond is committed from the  !
      !      scratch hydrology solve (column_fast_step_ark), exactly as theta's clip/floor           !
      !      corrections are, so writing it here would commit it twice. THIS IS THE COMMIT PATH: a   !
      !      field genuinely omitted here is computed, balanced by every ledger, and then dropped on !
      !      the floor -- which is why the exclusion is stated rather than left to be inferred, and  !
      !      asserted in test_state_combinators.  ----------------------------------------------------!
   end subroutine unpack_column_state

   !----- Re-diagnose every layer's temperature and liquid fraction from the committed energy + water. -!
   pure subroutine diagnose_soil_temps(y_out, dry_hcap, nsl, soil_temp, soil_fliq)
      type(column_state_t), intent(in)    :: y_out
      real(wp),             intent(in)    :: dry_hcap(:)
      integer(ik),          intent(in)    :: nsl
      real(wp),             intent(inout) :: soil_temp(:), soil_fliq(:)
      integer(ik) :: k
      do k = 1_ik, nsl
         call internal_energy_to_temp(y_out%soil_energy(k), y_out%theta(k)*rho_h2o, dry_hcap(k), soil_temp(k), soil_fliq(k))
      end do
   end subroutine diagnose_soil_temps

   !---------------------------------------------------------------------------------------!
   ! assemble_soil_energy_forcing -- the soil-heat column's boundary and volumetric forcing from  !
   ! the surface ground heat flux and the water fluxes that advect enthalpy through it. The ONE    !
   ! assembler for the ARK stage (meds_fast_ark%column_be_stage, faces and drainage from the frozen !
   ! scratch hydrology, plus the scratch solve's clip/floor corrections) and the whole-column RHS   !
   ! (meds_fast_time_derivs%column_derivs, faces and drainage from the stage's OWN water tendency).  !
   ! Which water trajectory the faces come from is the caller's decision -- borrowing another        !
   ! solve's faces while committing your own theta is the defect class the 2026-07 RK45 soil-surface !
   ! blow-up belonged to -- so this routine only lays the numbers out with one sign convention:      !
   !                                                                                                  !
   !   * face_flux(k) [m/s] is the hydrology's DOWN-positive water flux at layer k's face; the energy !
   !     kernel wants UP-positive, hence the sign flip;                                               !
   !   * root_heat_sink(k) [W/m2] is the enthalpy the roots extract (qloss_total by root_share), plus  !
   !     the optional per-layer corrections `sink_add - sink_sub` (ARK: the scratch clip valued at    !
   !     each layer's state^n temperature ADDS, the theta_res floor SUBTRACTS), plus the drainage      !
   !     enthalpy e_drain leaving through the bottom layer;                                           !
   !   * the TOP face carries infiltration at the pond temperature t_infil as a kernel term with the   !
   !     same upwind rule as the interior faces; the BOTTOM face is 0 here because the drainage        !
   !     enthalpy is charged through root_heat_sink(nsl) instead (the ledger's convention). There is   !
   !     deliberately NO runoff term: runoff leaves the POND, not soil layer 1.                        !
   !---------------------------------------------------------------------------------------!
   pure subroutine assemble_soil_energy_forcing(eforc, nsl, g_top, geothermal, theta, root_share,       &
                                                qloss_total, face_flux, infiltration, t_infil, e_drain, &
                                                sink_add, sink_sub)
      type(energy_forcing_t), intent(out) :: eforc
      integer(ik),            intent(in)  :: nsl
      real(wp),               intent(in)  :: g_top          !< [W/m2]  net ground heat flux into the soil top
      real(wp),               intent(in)  :: geothermal     !< [W/m2]  bottom boundary flux
      real(wp),               intent(in)  :: theta(:)       !< [m3/m3] soil moisture (thermal properties)
      real(wp),               intent(in)  :: root_share(:)  !< [-]     static root-uptake share per layer
      real(wp),               intent(in)  :: qloss_total    !< [W/m2]  enthalpy advected out with root uptake
      real(wp),               intent(in)  :: face_flux(:)   !< [m/s]   DOWN-positive water flux at layer faces
      real(wp),               intent(in)  :: infiltration   !< [kg/m2/s] top-face infiltration (pond -> soil)
      real(wp),               intent(in)  :: t_infil        !< [K]     temperature of the infiltrating water
      real(wp),               intent(in)  :: e_drain        !< [W/m2]  enthalpy leaving with bottom drainage
      real(wp), optional,     intent(in)  :: sink_add(:), sink_sub(:)   !< [W/m2] per-layer sink corrections
      integer(ik) :: k
      eforc%g_top = g_top ; eforc%geothermal = geothermal
      do k = 1_ik, nsl
         eforc%soil_water(k)     = theta(k)
         if (present(sink_add)) then
            eforc%root_heat_sink(k) = qloss_total * root_share(k) + sink_add(k) - sink_sub(k)
         else
            eforc%root_heat_sink(k) = qloss_total * root_share(k)
         end if
         eforc%w_flux(k)         = -face_flux(k)
      end do
      eforc%w_flux_top  = -infiltration / rho_h2o
      eforc%t_water_top = t_infil
      eforc%w_flux_bot  = 0.0_wp
      eforc%root_heat_sink(nsl) = eforc%root_heat_sink(nsl) + e_drain
   end subroutine assemble_soil_energy_forcing

   !---------------------------------------------------------------------------------------!
   ! apply_process_mask -- hold every masked-OFF process at its start-of-step state. The process    !
   ! mask reduces the column for diagnosis (a frozen store still exchanges with its neighbours, so a  !
   ! reduced column cannot conserve -- see mask_is_full); the march integrates everything and the     !
   ! masked fields are then restored here, once, on the committed state. ONE routine for both        !
   ! schemes: they used to carry their own copies, which had already diverged on the pond fields     !
   ! (2026-09 review, item 5 #3). The pond rides with soil water; on a scheme that passes it through  !
   ! the stages the restore is an identity.                                                           !
   !---------------------------------------------------------------------------------------!
   pure subroutine apply_process_mask(mask, y, y_out, n, nsl)
      type(process_mask_t), intent(in)    :: mask
      type(column_state_t), intent(in)    :: y        !< start-of-step state
      type(column_state_t), intent(inout) :: y_out    !< integrated state, masked fields restored
      integer(ik),          intent(in)    :: n, nsl
      if (.not. mask%cas_energy) y_out%cas_enthalpy        = y%cas_enthalpy
      if (.not. mask%cas_vapour) y_out%cas_shv             = y%cas_shv
      if (.not. mask%cas_co2)    y_out%cas_co2             = y%cas_co2
      if (.not. mask%soil_heat)  y_out%soil_energy(1:nsl)  = y%soil_energy(1:nsl)
      if (.not. mask%soil_water) then
         y_out%theta(1:nsl)   = y%theta(1:nsl)
         y_out%w_surface      = y%w_surface
         y_out%w_surface_enth = y%w_surface_enth
      end if
      if (.not. mask%hydraulics) then
         y_out%leaf_water_mass(1:n) = y%leaf_water_mass(1:n)
         y_out%wood_water_mass(1:n) = y%wood_water_mass(1:n)
         y_out%leaf_surf_water(1:n) = y%leaf_surf_water(1:n)
         y_out%wood_surf_water(1:n) = y%wood_surf_water(1:n)
      end if
   end subroutine apply_process_mask

end module meds_column_state_ops
