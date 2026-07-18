!==========================================================================================!
! meds_plant_carbon_dynamics -- stateless per-cohort CARBON budget + allocation kernels.     !
!                                                                                          !
! The mechanistic replacement for the phenomenological growth engine: it turns the daily      !
! carbon gain (NPP) into tissue growth, covering the science of ED2's growth_balive +          !
! structural_growth but UNIFIED into one step (woody growth is a daily demand, not a separate   !
! monthly cadence). It follows FATES PARTEH Hypothesis-1 ("Allometrically Guided, Carbon         !
! Only"): every pool has an allometric target, the gain fills toward the targets in PRIORITY      !
! order, and the residual advances stature (wood). See docs/dev_plans/MEDS_PLANT_CARBON_DYNAMICS_DESIGN.md.!
!                                                                                          !
!   * tissue_turnover_rates  -- per-tissue turnover RATES [1/yr]. Constant baselines now, with   !
!       the evergreen cold-suppression factor; the seam exists so a future light / tropical-       !
!       phenology form can vary the leaf rate without touching allocation.                          !
!   * plant_carbon_allocation -- distribute the step's carbon (gains) net of tissue turnover        !
!       (losses) toward the allometric deficits (demands), gated by the phenology flag, returning    !
!       the NET carbon flux to each pool (npp). Pure distribution arithmetic: targets/deficits are    !
!       computed upstream (allometry) and passed in, so the kernel is state-free and GPU/SIMD-safe.    !
!                                                                                          !
! ALL quantities are CARBON [kgC/plant]; every biomass<->carbon conversion is done once at param    !
! initialization, so the kernels never convert. Amounts are integrated over one step (rate x pool x   !
! dt is the caller's job). The public seams are re-exported through meds_plant_interface.               !
!                                                                                          !
! Carbon closure (holds every call): (npp_leaf + npp_fineroot + npp_wood + npp_nonstructural + npp_repro)!
!                                    - deficit  ==  net_carbon - (turnover_leaf + turnover_fineroot).    !
! `deficit` (unpaid respiration once storage is exhausted) and `starving` flag a cohort the STATEFUL    !
! updater (src/demography, next PR) must resolve by destroying tissue -- this pure kernel never          !
! mutates a pool; it only reports the transfers.                                                        !
!==========================================================================================!
module meds_plant_carbon_dynamics
   use meds_kinds,       only : wp, ik
   use meds_constants,   only : safe_exp
   use meds_plant_types, only : turnover_env_t, turnover_params_t, turnover_rates_t,           &
                                carbon_gain_t, carbon_loss_t, carbon_demand_t, carbon_npp_t,   &
                                PHEN_ON
   implicit none
   private

   public :: tissue_turnover_rates, plant_carbon_allocation

   !----- Evergreen cold-suppression of turnover (ED2 form): factor = 1/(1+exp(k*(T0 - T))),  !
   !       ~1 when warm, ~0 in the cold, 0.5 at T0. -----------------------------------------!
   real(wp), parameter :: evg_ref_temp = 278.15_wp   !< [K]    5 degC reference
   real(wp), parameter :: evg_slope    = 0.4_wp      !< [1/K]  suppression sharpness

contains

   !---------------------------------------------------------------------------------------!
   ! Effective per-tissue turnover rates [1/yr] for one cohort. Baselines come straight from  !
   ! the PFT traits; evergreen PFTs cold-suppress turnover (dormant leaves persist), so their  !
   ! rates are scaled by the ED2 evergreen factor. Deciduous PFTs use the flat baseline (their  !
   ! leaf display is governed by phenology, not by a temperature response here).               !
   !---------------------------------------------------------------------------------------!
   pure subroutine tissue_turnover_rates(env, params, rates)
      type(turnover_env_t),    intent(in)  :: env
      type(turnover_params_t), intent(in)  :: params
      type(turnover_rates_t),  intent(out) :: rates
      real(wp) :: f_t

      f_t = 1.0_wp
      if (params%evergreen) f_t = 1.0_wp / (1.0_wp + safe_exp(evg_slope * (evg_ref_temp - env%tissue_temp)))
      rates%leaf     = params%leaf_turnover_rate     * f_t
      rates%fineroot = params%fineroot_turnover_rate * f_t
   end subroutine tissue_turnover_rates

   !---------------------------------------------------------------------------------------!
   ! Allocate the step's carbon among the pools (PARTEH priority ladder, wood as the residual !
   ! sink). Inputs are amounts over the step [kgC/plant]: `gain` (NPP after growth resp, plus   !
   ! the storage pool available to draw), `loss` (leaf/fine-root turnover -> litter), `demand`   !
   ! (allometric deficits toward each target; wood is the residual growth demand), and the       !
   ! phenology status flag. Output `npp` is the NET carbon flux to each pool (allocation minus    !
   ! that pool's turnover), signed.                                                              !
   !                                                                                          !
   ! Priority when net_carbon >= 0 (fund from NPP first, then storage where allowed):            !
   !   P1  replace leaf + fine-root turnover      (storage allowed: keep foliage/roots whole)     !
   !   P2  grow leaf + fine-root toward target    (storage allowed only when PHEN_ON: flush)       !
   !   P3  refill storage toward target           (NPP only)                                        !
   !   P4  reproduction = reproduction_fraction x residual   (NPP only; 0 below the maturity height)  !
   !   P5  grow wood                              (NPP only; residual sink -> npp_wood = 0 if none)  !
   !   leftover NPP banks to storage.                                                              !
   ! When net_carbon < 0 the respiration debt is paid from storage (no growth); if storage cannot  !
   ! cover it, `starving` is set and the shortfall is reported as `deficit` for the caller to        !
   ! resolve (tissue destruction lives in the stateful updater, not here).                           !
   !---------------------------------------------------------------------------------------!
   pure subroutine plant_carbon_allocation(gain, loss, demand, phenology_status, npp)
      type(carbon_gain_t),   intent(in)  :: gain
      type(carbon_loss_t),   intent(in)  :: loss
      type(carbon_demand_t), intent(in)  :: demand
      integer(ik),           intent(in)  :: phenology_status
      type(carbon_npp_t),    intent(out) :: npp
      real(wp) :: npp_left, store_left, draw, debt
      real(wp) :: a_leaf, a_fineroot, a_store, a_wood, a_repro
      logical  :: can_flush

      a_leaf = 0.0_wp ; a_fineroot = 0.0_wp ; a_store = 0.0_wp ; a_wood = 0.0_wp ; a_repro = 0.0_wp
      draw       = 0.0_wp
      store_left = max(gain%storage, 0.0_wp)
      npp%starving = .false.
      npp%deficit  = 0.0_wp

      if (gain%net_carbon >= 0.0_wp) then
         npp_left = gain%net_carbon
         !----- P1: replace leaf + fine-root turnover (draw storage if NPP short). ----------!
         call fund(loss%leaf,     .true., npp_left, store_left, draw, a_leaf)
         call fund(loss%fineroot, .true., npp_left, store_left, draw, a_fineroot)
         !----- P2: grow leaf + fine-root toward target; storage only during a PHEN_ON flush. !
         can_flush = (phenology_status == PHEN_ON)
         call fund(demand%leaf,     can_flush, npp_left, store_left, draw, a_leaf)
         call fund(demand%fineroot, can_flush, npp_left, store_left, draw, a_fineroot)
         !----- P3: refill storage toward target from remaining NPP only. -------------------!
         a_store  = min(max(demand%storage, 0.0_wp), max(npp_left, 0.0_wp))
         npp_left = npp_left - a_store
         !----- P4: reproduction -- a fraction of the post-storage residual (the fraction is    !
         !      0 below the maturity height, so immature cohorts reproduce nothing). ----------!
         a_repro  = max(demand%reproduction_fraction, 0.0_wp) * max(npp_left, 0.0_wp)
         npp_left = npp_left - a_repro
         !----- P5: grow wood from remaining NPP only (residual sink). ----------------------!
         a_wood   = min(max(demand%wood, 0.0_wp), max(npp_left, 0.0_wp))
         npp_left = npp_left - a_wood
         !----- Any leftover NPP banks to storage. -----------------------------------------!
         a_store  = a_store + max(npp_left, 0.0_wp)
      else
         !----- Negative NPP: pay the respiration debt from storage; no growth. -------------!
         debt = -gain%net_carbon
         draw = min(debt, store_left)
         if (draw < debt) then
            npp%starving = .true.
            npp%deficit  = debt - draw     ! unpaid resp -> caller destroys tissue (next PR)
         end if
      end if

      npp%leaf          = a_leaf     - loss%leaf
      npp%fineroot      = a_fineroot - loss%fineroot
      npp%wood          = a_wood
      npp%repro         = a_repro
      npp%nonstructural = a_store    - draw
   end subroutine plant_carbon_allocation

   !----- Fund a demand from NPP first, then (if allowed) from storage; accumulate what was    !
   !       obtained into `got` and the storage drawdown into `draw`. --------------------------!
   pure subroutine fund(need, allow_store, npp_left, store_left, draw, got)
      real(wp), intent(in)    :: need
      logical,  intent(in)    :: allow_store
      real(wp), intent(inout) :: npp_left, store_left, draw, got
      real(wp) :: want, take

      want = max(need, 0.0_wp)
      !----- from NPP first. --------------------------------------------------------------!
      take     = min(want, max(npp_left, 0.0_wp))
      got      = got + take
      npp_left = npp_left - take
      want     = want - take
      !----- then storage, if allowed and still needed. -----------------------------------!
      if (allow_store .and. want > 0.0_wp) then
         take       = min(want, store_left)
         got        = got + take
         store_left = store_left - take
         draw       = draw + take
      end if
   end subroutine fund

end module meds_plant_carbon_dynamics
