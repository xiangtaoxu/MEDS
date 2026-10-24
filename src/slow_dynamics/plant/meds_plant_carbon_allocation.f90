! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_carbon_allocation -- how a cohort spends (or fails to pay) its daily carbon.  !
! Every pool has an allometric target, the net carbon fills toward the targets in priority !
! order, and the residual builds wood (docs/science/plant_carbon_allocation.md §1-§3, §5). !
!   * plant_carbon_allocation -- the priority ladder: per-pool growth, the net storage     !
!       change, growth respiration, root exudate and the unpaid maintenance.               !
!   * growth_respiration      -- g x the carbon built (eq. (5)).                           !
! All amounts are carbon [kgC/plant]. The driver (cohort_carbon_demand) removes the step's !
! tissue loss before the call, so this kernel only adds carbon; it closes eq. (6).         !
!==========================================================================================!
module meds_plant_carbon_allocation
   use meds_kinds, only : wp
   implicit none
   private

   public :: plant_carbon_allocation, growth_respiration

contains

   !---------------------------------------------------------------------------------------!
   ! Daily allocation for one cohort down the ladder of plant_carbon_allocation.md §2:     !
   ! maintenance debt (from storage when net = gpp - resp_maint < 0) > leaf and fine-root  !
   ! growth (NPP, then storage) > storage refill (NPP only) > reproduction > wood (the     !
   ! residual), at (1+g) per unit of tissue and 1:1 for storage. The leaf and fine-root    !
   ! demands arrive flush-capped and net of the step's shed. Maintenance no pool can pay   !
   ! is `deficit`, which the driver nets out of the cohort's CO2 (plant_respiration.md §6). !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine plant_carbon_allocation(gpp, resp_maint, growth_resp_frac, storage, &
                                    leaf_demand, fineroot_demand, storage_demand, repro_frac,     &
                                    growth_leaf, growth_fineroot, growth_wood, npp_store,         &
                                    growth_repro, growth_resp, deficit, starving,                 &
                                    wood_growth_max, exudate)
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
      real(wp), intent(out) :: npp_store        !< [kgC/plant] net storage change (refill - drawdown; signed)
      real(wp), intent(out) :: growth_repro     !< [kgC/plant] carbon -> reproduction (>= 0; -> recruits)
      real(wp), intent(out) :: growth_resp      !< [kgC/plant] growth respiration on realized growth (>= 0)
      real(wp), intent(out) :: deficit          !< [kgC/plant] unpaid maintenance after storage exhausted (>= 0)
      logical,  intent(out) :: starving         !< .true. => storage could not cover the maintenance debt
      !----- The sink limit (meds_sink_limitation): pass both or neither; without them wood !
      !      takes the whole residual.                                                     !
      real(wp), intent(in),  optional :: wood_growth_max  !< [kgC/plant] sink-limited wood growth this step
      real(wp), intent(out), optional :: exudate          !< [kgC/plant] residual the sink cannot use (>= 0)
      real(wp) :: net, g, cost, avail, store_left, draw, debt, storage_refill, c_repro, c_wood, c_exudate

      growth_leaf = 0.0_wp
      growth_fineroot = 0.0_wp
      growth_wood = 0.0_wp
      growth_repro = 0.0_wp
      storage_refill = 0.0_wp
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
         !----- Pay the maintenance debt from storage first; leftover reserves may still grow. ---!
         debt = -net
         draw = min(debt, store_left)
         store_left = store_left - draw
         if (draw < debt) then
            starving = .true.
            deficit  = debt - draw
         end if
         avail = 0.0_wp
      end if

      !----- P1: leaf and fine-root growth toward the flush-capped demand, from NPP and then !
      !          storage, so a leafless canopy can flush on its reserves. -----------------!
      call fill_carbon_demand(leaf_demand,     g, cost, avail, store_left, draw, growth_leaf,     growth_resp)
      call fill_carbon_demand(fineroot_demand, g, cost, avail, store_left, draw, growth_fineroot, growth_resp)
      !----- P2: refill storage from remaining NPP only (no construction cost). -----------------!
      storage_refill = min(max(storage_demand, 0.0_wp), max(avail, 0.0_wp))
      avail   = avail - storage_refill
      !----- P3: reproduction -- a fraction of the post-storage residual, construction-charged. --!
      c_repro      = max(repro_frac, 0.0_wp) * max(avail, 0.0_wp)
      growth_repro = c_repro / cost
      growth_resp  = growth_resp + growth_respiration(growth_repro, g)
      avail        = avail - c_repro
      !----- P4: wood takes the residual up to the sink limit; the rest is root exudate at the !
      !          same (1+g) cost, so the limit leaves growth respiration unchanged. -------!
      c_wood      = max(avail, 0.0_wp)
      growth_wood = c_wood / cost
      c_exudate   = 0.0_wp
      if (present(wood_growth_max)) then
         if (growth_wood > wood_growth_max) then
            c_exudate   = growth_wood - max(wood_growth_max, 0.0_wp)
            growth_wood = max(wood_growth_max, 0.0_wp)
         end if
      end if
      growth_resp = growth_resp + growth_respiration(growth_wood + c_exudate, g)
      if (present(exudate)) exudate = c_exudate

      npp_store = storage_refill - draw
   end subroutine plant_carbon_allocation

   !---------------------------------------------------------------------------------------!
   ! Fund a tissue demand from NPP first, then from storage, at (1+g) carbon per unit built; !
   ! add the tissue to `built`, the storage drawn to `draw` and the respired carbon to     !
   ! `growth_resp` (docs/science/plant_carbon_allocation.md §3).                           !
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
   ! Growth (construction) respiration [kgC/plant] = g x the carbon actually built this    !
   ! step (docs/science/plant_carbon_allocation.md eq. (5)). A tissue- or                  !
   ! temperature-dependent construction cost would replace this one form.                  !
   !---------------------------------------------------------------------------------------!
   elemental pure function growth_respiration(npp_growth, growth_resp_factor) result(rg)
      real(wp), intent(in) :: npp_growth          !< [kgC/plant] carbon allocated to new growth this step
      real(wp), intent(in) :: growth_resp_factor  !< [--] per-PFT construction-cost fraction
      real(wp)             :: rg                   !< [kgC/plant] growth respiration
      rg = growth_resp_factor * max(0.0_wp, npp_growth)
   end function growth_respiration

end module meds_plant_carbon_allocation
