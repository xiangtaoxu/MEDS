!==========================================================================================!
! meds_slow_dynamics -- the THIN slow-tier coordinator (MEDS_SLOW_DYNAMICS_DESIGN.md Part II,   !
! section 10a). Sequences the slow-tier domain drivers as PEERS -- vegetation dynamics, then     !
! biogeochemistry -- and owns the shared per-patch slow-state application (update_patch_states,  !
! hoisted out of meds_vegetation_dynamics) between them. It stays a coordinator, not a numerical   !
! stepper: every actual computation lives in the domain drivers it calls.                          !
!==========================================================================================!
module meds_slow_dynamics
   use meds_kinds,                 only : wp, ik
   use meds_config,                only : meds_config_t
   use meds_core_interface,        only : site_t, update_patch_states
   use meds_vegetation_dynamics,   only : vegetation_dynamics
   use meds_biogeochem_dynamics,   only : advance_biogeochem_dynamics
   use meds_biogeochem_types,      only : litter_input_t
   implicit none
   private

   public :: advance_slow_dynamics

contains

   !---------------------------------------------------------------------------------------!
   ! Advance one slow step: vegetation dynamics (demographic rates + apply-primitives; returns  !
   ! this step's litter accumulator), the shared patch-ageing applier, then -- when soil_carbon_  !
   ! on -- the daily soil-carbon matrix step consuming that litter + the fast loop's day-          !
   ! integrated environmental scalar. Same signature as vegetation_dynamics (doy optional).        !
   !---------------------------------------------------------------------------------------!
   subroutine advance_slow_dynamics(site, cfg, is_new_month, is_new_year, doy, worst_rh_seam_gap)
      type(site_t),        intent(inout) :: site
      type(meds_config_t), intent(in)    :: cfg
      logical,             intent(in)    :: is_new_month, is_new_year
      integer(ik),         intent(in), optional :: doy
      real(wp),            intent(out),  optional :: worst_rh_seam_gap
      type(litter_input_t), allocatable :: lit(:)

      if (present(doy)) then
         call vegetation_dynamics(site, cfg, is_new_month, is_new_year, doy, lit)
      else
         call vegetation_dynamics(site, cfg, is_new_month, is_new_year, lit=lit)
      end if

      call update_patch_states(site%patch, cfg%dt_years)

      if (cfg%soil_carbon_on) call advance_biogeochem_dynamics(site, cfg, lit, worst_rh_seam_gap)
   end subroutine advance_slow_dynamics

end module meds_slow_dynamics
