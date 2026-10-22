! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_column_state_ops -- algebra and bookkeeping on the fast loop's column state, shared by both !
! schemes (meds_fast_ark, meds_fast_rk45): the state's fields listed once (state_to_array,          !
! array_to_state, tend_to_array, state_entry_rules) and the combinators on the flat array; the      !
! boundary-flux ledger (bflux_*); the stage clamp clamp_cas; the store totals the ledgers read; the  !
! soil-energy forcing; and the post-march commits (clamp_canopy_film, deposit_condensate,            !
! unpack_column_state, diagnose_soil_temps, commit_tissue_energy). A new state                      !
! field goes in the four listing routines, the type, the frozen-record pack, unpack_column_state     !
! and test_state_combinators; the combinators never name a field.                                    !
!==========================================================================================!
module meds_column_state_ops
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : rho_h2o, tiny_num, mmh2o, cp_liq
   use meds_therm_lib,        only : internal_energy_to_temp, cas_enthalpy_to_temp, cas_temp_to_enthalpy,  &
                                     tissue_energy_to_temp
   use meds_fast_types,       only : column_state_t, column_tend_t, stage_bflux_t, column_bflux_t, GRP_ENTH,      &
                                     GRP_SHV, GRP_CO2, GRP_SE, GRP_THETA, GRP_LEAF_W, GRP_WOOD_W,                 &
                                     patch_biophys_t, tissue_coefficients_t
   use meds_soil_types, only : energy_forcing_t
   implicit none
   private

   public :: state_to_array, state_entry_rules
   public :: state_init, state_accum, state_extrap, state_sub, state_err_diff, zero_like
   public :: bflux_zero, bflux_add, bflux_bweight, transp_per_leaf_area
   public :: assemble_soil_energy_forcing, commit_tissue_energy
   public :: clamp_cas
   public :: soil_water_store, soil_energy_store, plant_water_store, canopy_film_store
   public :: clamp_canopy_film, deposit_condensate, unpack_column_state, diagnose_soil_temps
   public :: state_size, array_to_state, state_axpy   ! used only by the tests

contains

   !=======================================================================================!
   ! The state's fields, listed once. The integrator's bookkeeping works on one flat array per     !
   ! state; the three routines below are the only ones that name column_state_t's fields, and     !
   ! state_entry_rules the only one that says how each is treated. Fortran cannot check that a     !
   ! field is listed, so test_state_combinators runs every field through every combinator.         !
   ! Layout: CAS enthalpy, humidity, CO2 | soil energy (1:nsl) | soil water (1:nsl) |             !
   !   per cohort: leaf, wood water | per cohort: leaf, wood film | pond: water, enthalpy         !
   !=======================================================================================!
   pure integer(ik) function state_size(n, nsl) result(m)
      integer(ik), intent(in) :: n, nsl
      m = 5_ik + 2_ik*nsl + 4_ik*n
   end function state_size

   pure subroutine state_to_array(y, n, nsl, v)
      type(column_state_t), intent(in)  :: y
      integer(ik),          intent(in)  :: n, nsl
      real(wp),             intent(out) :: v(:)
      integer(ik) :: i, o
      v(1) = y%cas_enthalpy
      v(2) = y%cas_shv
      v(3) = y%cas_co2
      v(4:3+nsl)       = y%soil_energy(1:nsl)
      v(4+nsl:3+2*nsl) = y%theta(1:nsl)
      o = 3_ik + 2_ik*nsl
      do i = 1_ik, n
         v(o+2*i-1) = y%leaf_water_mass(i)
         v(o+2*i) = y%wood_water_mass(i)
      end do
      o = o + 2_ik*n
      do i = 1_ik, n
         v(o+2*i-1) = y%leaf_surf_water(i)
         v(o+2*i) = y%wood_surf_water(i)
      end do
      o = o + 2_ik*n
      v(o+1) = y%w_surface
      v(o+2) = y%w_surface_energy
   end subroutine state_to_array

   !----- The inverse. The per-cohort arrays are allocated only when they are missing or of another !
   !      length, so a state reused from call to call is filled in place. Soil layers past nsl are 0. !
   pure subroutine array_to_state(v, n, nsl, y)
      real(wp),             intent(in)    :: v(:)
      integer(ik),          intent(in)    :: n, nsl
      type(column_state_t), intent(inout) :: y
      integer(ik) :: i, o
      call ensure_cohort_arrays(y, n)
      y%cas_enthalpy = v(1)
      y%cas_shv = v(2)
      y%cas_co2 = v(3)
      y%soil_energy = 0.0_wp
      y%theta = 0.0_wp
      y%soil_energy(1:nsl) = v(4:3+nsl)
      y%theta(1:nsl)       = v(4+nsl:3+2*nsl)
      o = 3_ik + 2_ik*nsl
      do i = 1_ik, n
         y%leaf_water_mass(i) = v(o+2*i-1)
         y%wood_water_mass(i) = v(o+2*i)
      end do
      o = o + 2_ik*n
      do i = 1_ik, n
         y%leaf_surf_water(i) = v(o+2*i-1)
         y%wood_surf_water(i) = v(o+2*i)
      end do
      o = o + 2_ik*n
      y%w_surface = v(o+1)
      y%w_surface_energy = v(o+2)
   end subroutine array_to_state

   !----- A tendency in the same layout. The pond has none (it is committed from the water network),  !
   !      so its entries are 0 and y + a*k leaves it as it was. ---------------------------------!
   pure subroutine tend_to_array(k, n, nsl, v)
      type(column_tend_t), intent(in)  :: k
      integer(ik),         intent(in)  :: n, nsl
      real(wp),            intent(out) :: v(:)
      integer(ik) :: i, o
      v(1) = k%d_cas_enthalpy
      v(2) = k%d_cas_shv
      v(3) = k%d_cas_co2
      v(4:3+nsl)       = k%dedt(1:nsl)
      v(4+nsl:3+2*nsl) = k%dtheta_dt(1:nsl)
      o = 3_ik + 2_ik*nsl
      do i = 1_ik, n
         v(o+2*i-1) = k%d_leaf_water_mass(i)
         v(o+2*i) = k%d_wood_water_mass(i)
      end do
      o = o + 2_ik*n
      do i = 1_ik, n
         v(o+2*i-1) = k%d_leaf_surf_water(i)
         v(o+2*i) = k%d_wood_surf_water(i)
      end do
      o = o + 2_ik*n
      v(o+1) = 0.0_wp
      v(o+2) = 0.0_wp
   end subroutine tend_to_array

   !----- How entry j of the flat array is treated: norm_group, its group in the step-size error norm !
   !      (0: not counted), and in_embedded_error. The films and the pond are outside the norm; the    !
   !      plant water, films and pond ride operator-split maps outside the stages, so their embedded    !
   !      difference is set to 0 rather than left as the extrapolation's rounding noise. ------------!
   pure subroutine state_entry_rules(j, n, nsl, norm_group, in_embedded_error)
      integer(ik), intent(in)  :: j, n, nsl
      integer(ik), intent(out) :: norm_group
      logical,     intent(out) :: in_embedded_error
      integer(ik) :: o
      o = 3_ik + 2_ik*nsl
      if (j == 1_ik) then
         norm_group = GRP_ENTH
         in_embedded_error = .true.
      else if (j == 2_ik) then
         norm_group = GRP_SHV
         in_embedded_error = .true.
      else if (j == 3_ik) then
         norm_group = GRP_CO2
         in_embedded_error = .true.
      else if (j <= 3_ik + nsl) then
         norm_group = GRP_SE
         in_embedded_error = .true.
      else if (j <= o) then
         norm_group = GRP_THETA
         in_embedded_error = .true.
      else if (j <= o + 2_ik*n) then
         norm_group = merge(GRP_LEAF_W, GRP_WOOD_W, mod(j - o, 2_ik) == 1_ik)
         in_embedded_error = .false.
      else
         norm_group = 0_ik
         in_embedded_error = .false.         ! films and pond
      end if
   end subroutine state_entry_rules

   !----- Allocate a state's per-cohort arrays for n cohorts, unless they already are. ------------!
   pure subroutine ensure_cohort_arrays(y, n)
      type(column_state_t), intent(inout) :: y
      integer(ik),          intent(in)    :: n
      if (allocated(y%leaf_water_mass) .and. allocated(y%wood_water_mass) .and.                   &
          allocated(y%leaf_surf_water) .and. allocated(y%wood_surf_water)) then
         if (size(y%leaf_water_mass) == n .and. size(y%wood_water_mass) == n .and.               &
             size(y%leaf_surf_water) == n .and. size(y%wood_surf_water) == n) return
      end if
      if (allocated(y%leaf_water_mass)) deallocate(y%leaf_water_mass)
      if (allocated(y%wood_water_mass)) deallocate(y%wood_water_mass)
      if (allocated(y%leaf_surf_water)) deallocate(y%leaf_surf_water)
      if (allocated(y%wood_surf_water)) deallocate(y%wood_surf_water)
      allocate(y%leaf_water_mass(n), y%wood_water_mass(n), y%leaf_surf_water(n), y%wood_surf_water(n))
   end subroutine ensure_cohort_arrays

   !=======================================================================================!
   ! The combinators, on the flat array. Each fills its output in place.                       !
   !=======================================================================================!
   !----- copy: y_stage = y. ---------------------------------------------------------------------!
   pure subroutine state_init(y, n, nsl, y_stage)
      type(column_state_t), intent(in)    :: y
      integer(ik),          intent(in)    :: n, nsl
      type(column_state_t), intent(inout) :: y_stage
      real(wp) :: v(5_ik + 2_ik*nsl + 4_ik*n)
      call state_to_array(y, n, nsl, v)
      call array_to_state(v, n, nsl, y_stage)
   end subroutine state_init

   !----- y_stage = y + a*k; a stage of several k's is state_init and repeated state_accum. ------!
   pure subroutine state_axpy(y, a, k, n, nsl, y_stage)
      type(column_state_t), intent(in)    :: y
      real(wp),             intent(in)    :: a
      type(column_tend_t),  intent(in)    :: k
      integer(ik),          intent(in)    :: n, nsl
      type(column_state_t), intent(inout) :: y_stage
      real(wp) :: vy(5_ik + 2_ik*nsl + 4_ik*n), vk(5_ik + 2_ik*nsl + 4_ik*n)
      call state_to_array(y, n, nsl, vy)
      call tend_to_array(k, n, nsl, vk)
      vy = vy + a * vk
      call array_to_state(vy, n, nsl, y_stage)
   end subroutine state_axpy

   !----- y_stage = y_stage + a*k, in place. ------------------------------------------------------!
   pure subroutine state_accum(y_stage, a, k, n, nsl)
      type(column_state_t), intent(inout) :: y_stage
      real(wp),             intent(in)    :: a
      type(column_tend_t),  intent(in)    :: k
      integer(ik),          intent(in)    :: n, nsl
      real(wp) :: vy(5_ik + 2_ik*nsl + 4_ik*n), vk(5_ik + 2_ik*nsl + 4_ik*n)
      call state_to_array(y_stage, n, nsl, vy)
      call tend_to_array(k, n, nsl, vk)
      vy = vy + a * vk
      call array_to_state(vy, n, nsl, y_stage)
   end subroutine state_accum

   !----- Ledger helpers: b-weight two stage rates into amounts over dt (weights (1-gamma, gamma)      !
   !      times dt); zero an accumulator; add one sub-step's amounts. --------------------------------!
   pure subroutine bflux_bweight(acc, s2, s3, dt, gam)
      type(column_bflux_t), intent(out) :: acc
      type(stage_bflux_t),  intent(in)  :: s2, s3
      real(wp),             intent(in)  :: dt, gam
      real(wp) :: b2, b3
      b2 = (1.0_wp - gam) * dt
      b3 = gam * dt
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
      !----- intent(out) zeroes the scalars; no empty structure constructor (docs/building.md traps). -!
      if (present(n)) then
         allocate(acc%tissue_leaf_int(n), acc%tissue_wood_int(n), acc%transp_coh(n))
         acc%tissue_leaf_int = 0.0_wp
         acc%tissue_wood_int = 0.0_wp
         acc%transp_coh = 0.0_wp
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
      !----- Only accepted sub-steps reach here, so the tissue integrals cover the same sub-steps as  !
      !      every other amount. -------------------------------------------------------------------!
      if (allocated(acc%tissue_leaf_int) .and. allocated(s%tissue_leaf_int)) then
         acc%tissue_leaf_int = acc%tissue_leaf_int + s%tissue_leaf_int
         acc%tissue_wood_int = acc%tissue_wood_int + s%tissue_wood_int
      end if
      if (allocated(acc%transp_coh) .and. allocated(s%transp_coh)) acc%transp_coh = acc%transp_coh + s%transp_coh
   end subroutine bflux_add

   !----- A cohort's transpiration over a step [kg/m2 ground] as the rate per leaf area the            !
   !      diagnostics report [mol H2O/m2 leaf/s]; zero for a leafless cohort. -------------------------!
   elemental pure function transp_per_leaf_area(amount, dt, nplant, leaf_area) result(rate)
      real(wp), intent(in) :: amount      !< [kg/m2 ground] over dt
      real(wp), intent(in) :: dt          !< [s]
      real(wp), intent(in) :: nplant      !< [plant/m2]
      real(wp), intent(in) :: leaf_area   !< [m2 leaf/plant]
      real(wp) :: rate
      rate = 0.0_wp
      if (nplant * leaf_area > tiny_num) rate = amount / dt / (nplant * leaf_area) / mmh2o
   end function transp_per_leaf_area

   !----- out = (1-b)*y + b*Y2, the ARS stage-3 extrapolation base (b = (1-gamma)/gamma,            !
   !      docs/science/numerical_scheme.md eq. (1)). The plant water, the films and the pond are    !
   !      frozen through the stages, so for them this is y again (up to rounding). -----------------!
   pure subroutine state_extrap(y, b, Y2, n, nsl, out)
      type(column_state_t), intent(in)    :: y, Y2
      real(wp),             intent(in)    :: b
      integer(ik),          intent(in)    :: n, nsl
      type(column_state_t), intent(inout) :: out
      real(wp) :: vy(5_ik + 2_ik*nsl + 4_ik*n), v2(5_ik + 2_ik*nsl + 4_ik*n)
      call state_to_array(y, n, nsl, vy)
      call state_to_array(Y2, n, nsl, v2)
      vy = (1.0_wp - b)*vy + b*v2
      call array_to_state(vy, n, nsl, out)
   end subroutine state_extrap

   !----- Clamp the extrapolated CAS enthalpy and humidity into a wide physical range, so an overshoot  !
   !      of the ARS base (b = 1 + sqrt(2)) cannot reach a temperature where qsat overflows; in range it !
   !      is an identity. nfire counts the calls that moved the state (column_budget_t's clamp_*). ---!
   pure subroutine clamp_cas(s, nfire)
      type(column_state_t), intent(inout) :: s
      integer(ik), optional, intent(inout) :: nfire
      real(wp) :: t, shv_c, enth_in
      real(wp), parameter :: T_LO = 180.0_wp, T_HI = 350.0_wp, SHV_LO = 1.0e-8_wp, SHV_HI = 0.06_wp
      enth_in = s%cas_enthalpy
      shv_c = min(max(s%cas_shv, SHV_LO), SHV_HI)
      t     = cas_enthalpy_to_temp(s%cas_enthalpy, shv_c)
      t     = min(max(t, T_LO), T_HI)
      s%cas_shv      = shv_c
      s%cas_enthalpy = cas_temp_to_enthalpy(t, shv_c)
      !----- An in-range state round-trips, so compare with the input: that also catches a          !
      !      humidity-only clamp, which moves the enthalpy. ------------------------------------------!
      if (present(nfire)) then
         if (s%cas_enthalpy /= enth_in) nfire = nfire + 1_ik
      end if
   end subroutine clamp_cas

   !----- The embedded error estimate of the ARS(2,2,2) step, err = (Y3 - base3) - (Y2 - y)          !
   !      (docs/science/numerical_scheme.md eq. (2)), for the entries state_entry_rules counts, and  !
   !      0 for the rest. ---------------------------------------------------------------------------!
   pure subroutine state_err_diff(Y3, base3, Y2, y, n, nsl, err)
      type(column_state_t), intent(in)    :: Y3, base3, Y2, y
      integer(ik),          intent(in)    :: n, nsl
      type(column_state_t), intent(inout) :: err
      real(wp)    :: v3(5_ik + 2_ik*nsl + 4_ik*n), vb(5_ik + 2_ik*nsl + 4_ik*n)
      real(wp)    :: v2(5_ik + 2_ik*nsl + 4_ik*n), vy(5_ik + 2_ik*nsl + 4_ik*n)
      integer(ik) :: j, grp
      logical     :: counted
      call state_to_array(Y3, n, nsl, v3)
      call state_to_array(base3, n, nsl, vb)
      call state_to_array(Y2, n, nsl, v2)
      call state_to_array(y, n, nsl, vy)
      do j = 1_ik, size(v3, kind=ik)
         call state_entry_rules(j, n, nsl, grp, counted)
         if (counted) then
            v3(j) = (v3(j) - vb(j)) - (v2(j) - vy(j))
         else
            v3(j) = 0.0_wp
         end if
      end do
      call array_to_state(v3, n, nsl, err)
   end subroutine state_err_diff

   !----- out = a - b, the state difference that forms the low-order embedded solution. --------!
   pure subroutine state_sub(a, b, n, nsl, out)
      type(column_state_t), intent(in)    :: a, b
      integer(ik),          intent(in)    :: n, nsl
      type(column_state_t), intent(inout) :: out
      real(wp) :: va(5_ik + 2_ik*nsl + 4_ik*n), vb(5_ik + 2_ik*nsl + 4_ik*n)
      call state_to_array(a, n, nsl, va)
      call state_to_array(b, n, nsl, vb)
      va = va - vb
      call array_to_state(va, n, nsl, out)
   end subroutine state_sub

   !----- A state for n cohorts with every field zero, so the step-size norm can read an error that  !
   !      is already a difference (state_wrms_grouped takes two states). ---------------------------!
   pure function zero_like(ref, n, nsl) result(z)
      type(column_state_t), intent(in) :: ref
      integer(ik),          intent(in) :: n, nsl
      type(column_state_t) :: z
      real(wp) :: v(5_ik + 2_ik*nsl + 4_ik*n)
      v = 0.0_wp
      call array_to_state(v, n, nsl, z)
   end function zero_like

   !---------------------------------------------------------------------------------------!
   ! Whole-column store totals for the conservation ledgers [per m2 ground]. Plain ascending loops, !
   ! shared by both schemes, so a ledger built from these is the same arithmetic on both.           !
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
   ! Canopy-film capacity clamp on the committed state: cap each cohort's leaf/wood film at        !
   ! leaf_max*LAI / wood_max*WAI and floor it at 0; the clipped excess (overflow) and the floored    !
   ! shortfall (deficit) are returned so the caller books them in the whole-column ledger.          !
   !---------------------------------------------------------------------------------------!
   pure subroutine clamp_canopy_film(y_out, lai, wai, leaf_max, wood_max, n, surf_overflow, surf_deficit)
      type(column_state_t), intent(inout) :: y_out
      real(wp),             intent(in)    :: lai(:), wai(:)
      real(wp),             intent(in)    :: leaf_max(:), wood_max(:)   !< [kg/m2 leaf], [kg/m2 wood] per cohort
      integer(ik),          intent(in)    :: n
      real(wp),             intent(inout) :: surf_overflow, surf_deficit   !< [kg/m2] accumulated
      real(wp)    :: leaf_cap_i, wood_cap_i
      integer(ik) :: i
      do i = 1_ik, n
         leaf_cap_i = leaf_max(i) * lai(i)
         wood_cap_i = wood_max(i) * wai(i)
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
   ! Unpack the committed column state into the patch biophysics record: CAS twins (and their      !
   ! temperature), soil energy and moisture, plant internal water, canopy films; not the pond.    !
   !---------------------------------------------------------------------------------------!
   pure subroutine unpack_column_state(y_out, n, nsl, biophys)
      type(column_state_t),  intent(in)    :: y_out
      integer(ik),           intent(in)    :: n, nsl
      type(patch_biophys_t), intent(inout) :: biophys
      biophys%cas%can_enthalpy = y_out%cas_enthalpy
      biophys%cas%can_shv = y_out%cas_shv
      biophys%cas%can_co2 = y_out%cas_co2
      biophys%cas%can_temp = cas_enthalpy_to_temp(y_out%cas_enthalpy, y_out%cas_shv)
      biophys%soil_energy%soil_energy(1:nsl) = y_out%soil_energy(1:nsl)
      biophys%soil_water%theta(1:nsl)       = y_out%theta(1:nsl)
      biophys%leaf_water_mass(1:n) = y_out%leaf_water_mass(1:n)
      biophys%wood_water_mass(1:n) = y_out%wood_water_mass(1:n)
      biophys%leaf_surf_water(1:n) = y_out%leaf_surf_water(1:n)
      biophys%wood_surf_water(1:n) = y_out%wood_surf_water(1:n)
      !----- The pond (w_surface, w_surface_energy) is left out on purpose: each scheme commits it from  !
      !      the water network, so writing it here would commit it twice (test_state_combinators). ---!
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

   !------------------------------------------------------------------------------------------------!
   ! commit_tissue_energy -- each tissue's end-of-step energy and temperature, and the tissue store !
   ! the whole-column ledger counts (docs/science/vegetation_energy_dynamics.md eqs. (9)-(10)). The  !
   ! energy gains the heat cap*(t_avg - t0), t_avg being the march's time-average of where the      !
   ! tissue relaxed to, and the enthalpy its water gained (docs/science/soil_plant_water.md eq. (5)); !
   ! the temperature is read back with the end-of-step water. A masked-off process adds nothing; a  !
   ! leaf with no heat capacity takes its wood's temperature.                                        !
   !------------------------------------------------------------------------------------------------!
   pure subroutine commit_tissue_energy(tissue, nplant, dt, t_leaf_avg, t_wood_avg, n,                          &
                                        biophys, store0, store1)
      type(tissue_coefficients_t), intent(in)    :: tissue
      real(wp),                    intent(in)    :: nplant(:)                     !< [plant/m2]
      real(wp),                    intent(in)    :: dt                            !< [s] dt_fast
      real(wp),                    intent(in)    :: t_leaf_avg(:), t_wood_avg(:)  !< [K]
      integer(ik),                 intent(in)    :: n
      type(patch_biophys_t),       intent(inout) :: biophys
      real(wp),                    intent(out)   :: store0, store1                !< [J/m2] the tissue store
      real(wp)    :: heat_leaf, heat_wood, enth_leaf, enth_wood, w_leaf, w_wood
      integer(ik) :: i
      store0 = 0.0_wp
      store1 = 0.0_wp
      do i = 1_ik, n
         heat_leaf = tissue%leaf_hcap_per_dt(i) * dt * (t_leaf_avg(i) - tissue%t_leaf0(i))
         heat_wood = tissue%wood_hcap_per_dt(i) * dt * (t_wood_avg(i) - tissue%t_wood0(i))
         enth_leaf = tissue%leaf_water_enth(i) * dt
         enth_wood = tissue%wood_water_enth(i) * dt
         store0 = store0 + biophys%leaf_energy(i) + biophys%wood_energy(i)
         biophys%leaf_energy(i) = biophys%leaf_energy(i) + heat_leaf + enth_leaf
         biophys%wood_energy(i) = biophys%wood_energy(i) + heat_wood + enth_wood
         store1 = store1 + biophys%leaf_energy(i) + biophys%wood_energy(i)
         w_leaf = max(biophys%leaf_water_mass(i), 0.0_wp) * nplant(i)
         w_wood = max(biophys%wood_water_mass(i), 0.0_wp) * nplant(i)
         if (tissue%wood_dry_hcap(i) + w_wood * cp_liq > 0.0_wp)                                        &
            biophys%wood_temp(i) = tissue_energy_to_temp(biophys%wood_energy(i), w_wood, tissue%wood_dry_hcap(i))
         biophys%leaf_temp(i) = biophys%wood_temp(i)
         if (tissue%leaf_dry_hcap(i) + w_leaf * cp_liq > 0.0_wp)                                        &
            biophys%leaf_temp(i) = tissue_energy_to_temp(biophys%leaf_energy(i), w_leaf, tissue%leaf_dry_hcap(i))
      end do
   end subroutine commit_tissue_energy

   !---------------------------------------------------------------------------------------!
   ! assemble_soil_energy_forcing -- the soil-heat column's forcing from the ground heat flux and the  !
   ! water network's flows, for the ARK stage and column_derivs alike                                  !
   ! (docs/science/numerical_scheme.md §5b, docs/science/soil_plant_water.md §5). face_flux is down-   !
   ! positive (the kernel wants up, hence the flip); water_heat_out is the enthalpy leaving each layer !
   ! other than across a face, plus the drainage's in the bottom layer, whose face is therefore 0; the !
   ! top face carries infiltration at t_infil with the interior faces' upwind rule. No runoff term:    !
   ! runoff leaves the pond.                                                                           !
   !---------------------------------------------------------------------------------------!
   pure subroutine assemble_soil_energy_forcing(eforc, nsl, g_top, theta, water_heat_out,                  &
                                                face_flux, infiltration, t_infil, e_drain)
      type(energy_forcing_t), intent(out) :: eforc
      integer(ik),            intent(in)  :: nsl
      real(wp),               intent(in)  :: g_top          !< [W/m2]  net ground heat flux into the soil top
      real(wp),               intent(in)  :: theta(:)       !< [m3/m3] soil moisture (thermal properties)
      real(wp),               intent(in)  :: water_heat_out(:) !< [W/m2] enthalpy leaving each layer with water
      real(wp),               intent(in)  :: face_flux(:)   !< [m/s]   DOWN-positive water flux at layer faces
      real(wp),               intent(in)  :: infiltration   !< [kg/m2/s] top-face infiltration (pond -> soil)
      real(wp),               intent(in)  :: t_infil        !< [K]     temperature of the infiltrating water
      real(wp),               intent(in)  :: e_drain        !< [W/m2]  enthalpy leaving with bottom drainage
      integer(ik) :: k
      eforc%g_top = g_top
      do k = 1_ik, nsl
         eforc%soil_water(k)     = theta(k)
         eforc%root_heat_sink(k) = water_heat_out(k)
         eforc%w_flux(k)         = -face_flux(k)
      end do
      eforc%w_flux_top  = -infiltration / rho_h2o
      eforc%t_water_top = t_infil
      eforc%root_heat_sink(nsl) = eforc%root_heat_sink(nsl) + e_drain
   end subroutine assemble_soil_energy_forcing

end module meds_column_state_ops
