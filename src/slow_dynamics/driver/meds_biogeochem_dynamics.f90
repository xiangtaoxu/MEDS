!==========================================================================================!
! meds_biogeochem_dynamics -- the slow-loop BIOGEOCHEMISTRY driver (MEDS_SLOW_DYNAMICS_DESIGN.md !
! Part II, B2): the peer of meds_vegetation_dynamics under meds_slow_dynamics. Per patch, it       !
! advances the CENTURY soil-carbon matrix ONE daily step: the litter accumulated by the             !
! vegetation driver this step (`lit`, already routed through build_litter_input) plus the fast       !
! loop's day-integrated environmental scalar (`site%patch%xi_accum`, accumulated once per (patch,     !
! sub-step) by column_prepass over the SAME frozen pool the fast loop's heterotrophic_respiration_     !
! matrix respired against) drive soil_carbon_step, which commits the pool change once, at day-end.      !
!                                                                                          !
! THE DOUBLE-COUNTING CONTRACT (design section 9): the fast loop respires the FROZEN daily pool      !
! (biophys%soil_carbon, seeded once per patch per day in meds_fast_dynamics, never mutated there); this   !
! module is the SOLE writer of the real site%patch%soil_carbon. So the day's total fast Rh equals      !
! the pool debit BY CONSTRUCTION -- rh_seam_gap = soil_carbon_step's rh_today - the fast loop's         !
! accumulated rh_fast_accum is an ASSERTION GUARD (should be ~0), not a live correction.                 !
!==========================================================================================!
module meds_biogeochem_dynamics
   use meds_kinds,            only : wp, ik
   use meds_config,           only : meds_config_t
   use meds_site_state_types, only : site_t
   use meds_biogeochem_types, only : litter_input_t, soilc_audit_t, n_soil_pool, IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND, &
                                     IP_STRUCT_SOIL, IP_MICR, IP_SLOW, IP_PASSIVE
   use meds_soil_biogeochem,  only : build_litter_input, soil_carbon_step
   use meds_slow_ledger,      only : slow_ledger_t, slow_ledger_declare
   implicit none
   private

   public :: advance_biogeochem_dynamics

contains

   !---------------------------------------------------------------------------------------!
   ! Advance the slow soil-carbon matrix one daily step, per patch. `lit` is vegetation_       !
   ! dynamics's per-patch litter accumulator (turnover + continuous-mortality carbon; cull-     !
   ! termination and disturbance-kill litter were already added directly onto site%patch%        !
   ! soil_carbon by the core engine -- see meds_demography_cohort_fusefiss/meds_demography_patch_fusefiss).   !
   ! `worst_rh_seam_gap` (optional) reports the worst |rh_today - rh_fast_accum| across patches    !
   ! for a caller to assert on (mirrors fast_dynamics's worst_energy/worst_water pattern).          !
   !---------------------------------------------------------------------------------------!
   subroutine advance_biogeochem_dynamics(site, cfg, worst_rh_seam_gap, ledger)
      type(site_t),          intent(inout) :: site
      type(meds_config_t),   intent(in)    :: cfg
      !----- The litter arriving today is site%patch%litter_in -- patch state, in lockstep with the
      !      patch array. It used to be a caller-supplied array sized BEFORE patch disturbance ran,
      !      and the loops below (both of them, including the ledger declaration) indexed it by the
      !      post-disturbance patch count. Bounds checking reports it as `Subscript #1 of the array
      !      LIT has value 12 which is greater than the upper bound of 11`; without bounds checking
      !      it is undefined memory going straight into the CENTURY source term u.
      real(wp), optional,    intent(out)   :: worst_rh_seam_gap
      !----- The site ledger (plan §10.2). The litter the vegetation driver declared LEAVING the   !
      !      live pools arrives here; declaring the same quantity at both ends turns the litter    !
      !      seam into something the ledger TESTS rather than something it has to be told to       !
      !      ignore. Rh is not declared: it leaves the CENTURY pools into the canopy air, and BOTH  !
      !      are stores this ledger carries, so it is an internal transfer.  -----------------------!
      type(slow_ledger_t), optional, intent(inout) :: ledger
      real(wp)             :: u(n_soil_pool), lignin_in(2), xi_int(n_soil_pool), rh_today
      type(soilc_audit_t)  :: audit
      real(wp)             :: worst
      integer(ik)          :: ip

      worst = 0.0_wp
      if (present(ledger)) then
         do ip = 1_ik, site%patch%n
            call slow_ledger_declare(ledger, carbon_in = site%patch%area(ip)                       &
                     * (site%patch%litter_in(ip)%labile_grnd + site%patch%litter_in(ip)%labile_soil                                  &
                      + site%patch%litter_in(ip)%struct_grnd + site%patch%litter_in(ip)%struct_soil))
         end do
      end if
      do ip = 1_ik, site%patch%n
         call build_litter_input(site%patch%litter_in(ip), u, lignin_in)
         associate (xa => site%patch%xi_accum(ip))
            xi_int(IP_FAST_GRND)   = xa%fast_grnd
            xi_int(IP_FAST_SOIL)   = xa%fast_soil
            xi_int(IP_STRUCT_GRND) = xa%struct_grnd
            xi_int(IP_STRUCT_SOIL) = xa%struct_soil
            xi_int(IP_MICR)        = xa%microbial
            xi_int(IP_SLOW)        = xa%slow
            xi_int(IP_PASSIVE)     = xa%passive
            call soil_carbon_step(site%patch%soil_carbon(ip), u, lignin_in, xi_int, cfg%soil_carbon, &
                                  rh_today, audit)
            !----- The assertion guard (design section 9): with the frozen-pool identity, rh_today   !
            !      (the matrix's own pool debit) should equal the fast loop's independently            !
            !      accumulated Rh to within floating-point/adaptive-substep error. -------------------!
            audit%rh_fast_accum = xa%rh_fast_accum
            audit%rh_seam_gap   = rh_today - xa%rh_fast_accum
            worst = max(worst, abs(audit%rh_seam_gap))
            !----- LEDGER: heterotrophic respiration is the MIRROR of the GPP handover. The matrix  !
            !      debits the real pools HERE, inside the slow step; the fast loop already credited !
            !      that carbon to the canopy air yesterday, against the FROZEN copy of the pools    !
            !      (the double-counting contract in this module's header). One end of the transfer  !
            !      is inside this ledger's window and the other is outside it, so the debit must be !
            !      declared or it reads as the soil losing carbon to nowhere. Declaring rh_today --  !
            !      the matrix's own debit rather than the fast loop's accumulation -- is deliberate: !
            !      what is left in the residual is then rh_seam_gap, which is the quantity the       !
            !      contract actually asserts on.  ------------------------------------------------!
            if (present(ledger)) call slow_ledger_declare(ledger,                                  &
                                          carbon_out = site%patch%area(ip) * rh_today)
         end associate
      end do
      if (present(worst_rh_seam_gap)) worst_rh_seam_gap = worst
   end subroutine advance_biogeochem_dynamics

end module meds_biogeochem_dynamics
