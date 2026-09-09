!==========================================================================================!
! meds_plant_carbon_allocation -- how a plant spends (or fails to pay) its daily carbon.    !
!                                                                                          !
! The stateless, ELEMENTAL kernel that distributes one cohort's daily carbon budget to tissue  !
! GROWTH. It is the mechanistic replacement for the phenomenological growth engine, covering    !
! the science of ED2's growth_balive + structural_growth UNIFIED into one daily step, and it    !
! follows FATES PARTEH Hypothesis-1 ("Allometrically Guided, Carbon Only"): every pool has an    !
! allometric target, the net gain fills toward the targets in PRIORITY order, and the residual   !
! advances stature (wood). See docs/science/plant_carbon_allocation.md +                          !
! docs/dev_plans/MEDS_PLANT_CARBON_ALLOCATION_REFACTOR_DESIGN.md.                                       !
!                                                                                          !
!   * plant_carbon_allocation -- MASTER. Net carbon = GPP - maintenance respiration; distribute  !
!       it (plus storage where allowed) down the priority ladder to GROWTH, charging GROWTH        !
!       respiration (1+g per unit of tissue built) on realized growth only. It returns the         !
!       per-pool GROWTH (>= 0) + the net storage change + growth respiration + starvation flags.   !
!       Tissue LOSS (turnover / shed) is applied UPSTREAM by the driver (update_biomass_turnover): !
!       the pools handed in are already net of this step's shed, and the driver forms the net      !
!       per-pool change (growth - shed) + the litter. So this kernel only ever ADDS carbon.        !
!   * growth_respiration -- construction respiration = g x (carbon allocated to NEW growth). It     !
!       lives here (with the growth it charges, not with maintenance respiration) and takes the     !
!       realized growth npp_growth, opening the door to a more mechanistic form later.              !
!                                                                                          !
! ALL quantities are CARBON [kgC/plant]; every biomass<->carbon conversion is done once at param   !
! initialization, so the kernel never converts. Everything is elemental + scalar arithmetic, so    !
! the driver maps it straight over the cohort Structure-of-Arrays (GPU/SIMD-safe, issue-#7 N/A).   !
!                                                                                          !
! Carbon closure (holds every call, growth side):                                               !
!   (growth_leaf + growth_fineroot + growth_wood + growth_repro + npp_store) - deficit             !
!        == (gpp - resp_maint) - growth_resp                                                       !
! The driver then subtracts this step's shed (litter) from leaf/fine-root to get the net pool      !
! change. `deficit` (unpaid maintenance once storage is exhausted) + `starving` flag a cohort the  !
! STATEFUL updater must resolve by destroying tissue; this pure kernel never mutates a pool.       !
!==========================================================================================!
module meds_plant_carbon_allocation
   use meds_kinds, only : wp
   implicit none
   private

   public :: plant_carbon_allocation, growth_respiration

contains

   !---------------------------------------------------------------------------------------!
   ! MASTER daily carbon allocation for one cohort. Net carbon = gpp - resp_maint; building one !
   ! unit of growth tissue costs (1+g) carbon (the fraction g respired as GROWTH respiration),   !
   ! while storage refill is 1:1. leaf/fine-root demands arrive already FLUSH-CAPPED (from the    !
   ! phenology flush rate) and computed against the POST-SHED pool, so a dormant canopy presents   !
   ! zero growth demand and a just-shed canopy is refilled here within the same step.             !
   !                                                                                          !
   ! Priority: maintenance debt (from storage, if net < 0) > leaf/fine-root growth (NPP then       !
   ! storage) > storage refill (NPP only) > reproduction > wood (residual). Storage funds leaf/    !
   ! fine-root growth EVEN WHEN net < 0 (spring leaf-out is storage-funded), after the maintenance !
   ! debt is paid; a plant that cannot even cover maintenance (starving) has no reserves to grow.  !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine plant_carbon_allocation(gpp, resp_maint, growth_resp_frac, storage, &
                                    leaf_demand, fineroot_demand, storage_demand, repro_frac,     &
                                    growth_leaf, growth_fineroot, growth_wood, npp_store,         &
                                    growth_repro, growth_resp, deficit, starving)
      real(wp), intent(in)  :: gpp              !< [kgC/plant] gross primary production this step (>= 0)
      real(wp), intent(in)  :: resp_maint       !< [kgC/plant] maintenance respiration this step (>= 0)
      real(wp), intent(in)  :: growth_resp_frac !< [--] construction-cost fraction g (PFT trait)
      real(wp), intent(in)  :: storage          !< [kgC/plant] current nonstructural carbon (available to draw)
      real(wp), intent(in)  :: leaf_demand      !< [kgC/plant] flush-capped leaf growth demand (post-shed; >= 0)
      real(wp), intent(in)  :: fineroot_demand  !< [kgC/plant] flush-capped fine-root growth demand (post-shed; >= 0)
      real(wp), intent(in)  :: storage_demand   !< [kgC/plant] deficit toward the storage target (>= 0)
      real(wp), intent(in)  :: repro_frac       !< [--] fraction of the post-storage residual -> reproduction
      real(wp), intent(out) :: growth_leaf      !< [kgC/plant] carbon grown into leaf (>= 0)
      real(wp), intent(out) :: growth_fineroot  !< [kgC/plant] carbon grown into fine root (>= 0)
      real(wp), intent(out) :: growth_wood      !< [kgC/plant] carbon grown into wood (residual; >= 0)
      real(wp), intent(out) :: npp_store        !< [kgC/plant] NET storage change (refill - drawdown; signed)
      real(wp), intent(out) :: growth_repro     !< [kgC/plant] carbon -> reproduction (>= 0; -> recruits)
      real(wp), intent(out) :: growth_resp      !< [kgC/plant] growth respiration on realized growth (>= 0)
      real(wp), intent(out) :: deficit          !< [kgC/plant] unpaid maintenance after storage exhausted (>= 0)
      logical,  intent(out) :: starving         !< .true. => storage could not cover the maintenance debt
      real(wp) :: net, g, cost, avail, store_left, draw, debt, a_store, c_repro, c_wood

      growth_leaf = 0.0_wp
      growth_fineroot = 0.0_wp
      growth_wood = 0.0_wp
      growth_repro = 0.0_wp
      a_store = 0.0_wp
      draw = 0.0_wp
      growth_resp = 0.0_wp
      deficit = 0.0_wp
      starving = .false.
      store_left = max(storage, 0.0_wp)
      net = gpp - resp_maint
      g   = max(growth_resp_frac, 0.0_wp)
      cost = 1.0_wp + g

      if (net >= 0.0_wp) then
         avail = net
      else
         !----- Pay the maintenance debt from storage FIRST; leftover reserves may still grow. --!
         debt = -net
         draw = min(debt, store_left)
         store_left = store_left - draw
         if (draw < debt) then
            starving = .true.
            deficit  = debt - draw
         end if
         avail = 0.0_wp
      end if

      !----- P1: grow leaf + fine-root toward the flush-capped demand -- NPP first, then storage  !
      !          (so a leafless canopy can flush from reserves even when net < 0). --------------!
      call fill_carbon_demand(leaf_demand,     g, cost, avail, store_left, draw, growth_leaf,     growth_resp)
      call fill_carbon_demand(fineroot_demand, g, cost, avail, store_left, draw, growth_fineroot, growth_resp)
      !----- P2: refill storage from remaining NPP only (no construction cost). -----------------!
      a_store = min(max(storage_demand, 0.0_wp), max(avail, 0.0_wp))
      avail   = avail - a_store
      !----- P3: reproduction -- a fraction of the post-storage residual, construction-charged. --!
      c_repro      = max(repro_frac, 0.0_wp) * max(avail, 0.0_wp)
      growth_repro = c_repro / cost
      growth_resp  = growth_resp + growth_respiration(growth_repro, g)
      avail        = avail - c_repro
      !----- P4: wood is the residual sink -- build everything left, construction-charged. -------!
      c_wood      = max(avail, 0.0_wp)
      growth_wood = c_wood / cost
      growth_resp = growth_resp + growth_respiration(growth_wood, g)

      npp_store = a_store - draw
   end subroutine plant_carbon_allocation

   !---------------------------------------------------------------------------------------!
   ! Fund a tissue demand [kgC] from NPP first, then (if any remains) from storage; each unit   !
   ! of tissue built consumes (1+g) carbon, of which g is GROWTH respiration. Accumulate the     !
   ! tissue into `built`, the storage carbon drawn into `draw`, and the respired carbon into     !
   ! `growth_resp`. elemental (called only from the elemental master).                          !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine fill_carbon_demand(demand, g, cost, avail, store_left, draw, built, growth_resp)
      real(wp), intent(in)    :: demand, g, cost
      real(wp), intent(inout) :: avail, store_left, draw, built, growth_resp
      real(wp) :: want, build
      want = max(demand, 0.0_wp)
      if (want <= 0.0_wp) return
      !----- from NPP first (carbon-limited: at most avail/cost units of tissue). --------------!
      build       = min(want, max(avail, 0.0_wp) / cost)
      avail       = avail - build * cost
      built       = built + build
      growth_resp = growth_resp + growth_respiration(build, g)
      want        = want - build
      !----- then from storage, if still needed and available. --------------------------------!
      if (want > 0.0_wp .and. store_left > 0.0_wp) then
         build       = min(want, store_left / cost)
         store_left  = store_left - build * cost
         draw        = draw + build * cost
         built       = built + build
         growth_resp = growth_resp + growth_respiration(build, g)
      end if
   end subroutine fill_carbon_demand

   !---------------------------------------------------------------------------------------!
   ! Growth (construction) respiration [kgC/plant] = g x the carbon allocated to NEW growth this !
   ! step. It lives here (with the growth it charges, not with maintenance respiration), and     !
   ! takes the REALIZED growth npp_growth -- not the pre-allocation balance -- so it charges only !
   ! what was actually built. A single, simple form today; the seam is where a more mechanistic   !
   ! (tissue- or temperature-dependent) construction cost would go.                              !
   !---------------------------------------------------------------------------------------!
   elemental pure function growth_respiration(npp_growth, growth_resp_factor) result(rg)
      real(wp), intent(in) :: npp_growth          !< [kgC/plant] carbon allocated to NEW growth this step
      real(wp), intent(in) :: growth_resp_factor  !< [--] per-PFT construction-cost fraction
      real(wp)             :: rg                   !< [kgC/plant] growth respiration
      rg = growth_resp_factor * max(0.0_wp, npp_growth)
   end function growth_respiration

end module meds_plant_carbon_allocation
