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
   use meds_biogeochem_types, only : litter_input_t, soilc_audit_t, soilc_seam_t, n_soil_pool,       &
                                     IP_FAST_GRND, IP_FAST_SOIL, IP_STRUCT_GRND,                     &
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
   ! `seam` (optional) accumulates the per-RUN worst of the three soil-carbon seam diagnostics --   !
   ! the fast/slow Rh reconciliation gap, the lignin passive-tracer residual, and lambda, the        !
   ! largest fraction of any pool a single slow step withdraws (mirrors fast_dynamics's              !
   ! worst_energy/worst_water pattern).                                                              !
   !---------------------------------------------------------------------------------------!
   subroutine advance_biogeochem_dynamics(site, cfg, seam, ledger)
      type(site_t),          intent(inout) :: site
      type(meds_config_t),   intent(in)    :: cfg
      !----- The litter arriving today is site%patch%litter_in -- patch state, in lockstep with the
      !      patch array. It used to be a caller-supplied array sized BEFORE patch disturbance ran,
      !      and the loops below (both of them, including the ledger declaration) indexed it by the
      !      post-disturbance patch count. Bounds checking reports it as `Subscript #1 of the array
      !      LIT has value 12 which is greater than the upper bound of 11`; without bounds checking
      !      it is undefined memory going straight into the CENTURY source term u.
      !----- The seam diagnostics, accumulated as per-RUN worsts (intent(inout), reset by the driver !
      !      at open). One record rather than an optional real per number: worst_rh_seam_gap spent its !
      !      whole life unreported partly because adding the next one meant a second argument through  !
      !      three layers of driver.  ------------------------------------------------------------!
      type(soilc_seam_t), optional, intent(inout) :: seam
      !----- The site ledger (plan §10.2). The litter the vegetation driver declared LEAVING the   !
      !      live pools arrives here; declaring the same quantity at both ends turns the litter    !
      !      seam into something the ledger TESTS rather than something it has to be told to       !
      !      ignore. Rh is not declared: it leaves the CENTURY pools into the canopy air, and BOTH  !
      !      are stores this ledger carries, so it is an internal transfer.  -----------------------!
      type(slow_ledger_t), optional, intent(inout) :: ledger
      real(wp)             :: u(n_soil_pool), lignin_in(2), xi_int(n_soil_pool), rh_today
      type(soilc_audit_t)  :: audit
      integer(ik)          :: ip

      if (present(ledger)) then
         do ip = 1_ik, site%patch%n
            associate (lit => site%patch%litter_in(ip))
               call slow_ledger_declare(ledger, carbon_in = site%patch%area(ip)                    &
                        * (lit%labile_grnd + lit%labile_soil + lit%struct_grnd + lit%struct_soil))
            end associate
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
            if (present(seam)) then
               seam%worst_rh_gap = max(seam%worst_rh_gap, abs(audit%rh_seam_gap))
               seam%worst_lignin = max(seam%worst_lignin, abs(audit%lignin_resid))
               if (audit%lambda_max > seam%worst_lambda) then
                  seam%worst_lambda = audit%lambda_max
                  seam%lambda_pool  = audit%lambda_pool
               end if
            end if
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
   end subroutine advance_biogeochem_dynamics

end module meds_biogeochem_dynamics
