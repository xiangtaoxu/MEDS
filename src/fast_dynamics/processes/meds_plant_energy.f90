! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_plant_energy -- a cohort's leaf and wood energy                                     !
! (docs/science/vegetation_energy_dynamics.md): one tissue's balance relaxed exactly over  !
! the step (tissue_energy_balance, eqs. (1)-(7)), with its linearization coefficients and  !
! latent pathways; every cohort's wood and leaf at one canopy-air state                    !
! (canopy_tissue_energy); and each tissue's end-of-step energy and temperature             !
! (commit_tissue_energy). The wetted fraction its film pathway takes is the interception   !
! film's (meds_canopy_water).                                                              !
!==========================================================================================!
module meds_plant_energy
   use meds_kinds,            only : wp, ik
   use meds_constants,        only : tiny_num, stefan_boltzmann, latent_heat_vap, cp_liq
   use meds_therm_lib,        only : enthalpy_vapor, tissue_energy_to_temp
   use meds_canopy_water,     only : canopy_water_out_t
   implicit none
   private

   public :: plant_energy_env_t, latent_pathway_t
   public :: tissue_energy_balance, canopy_tissue_energy, commit_tissue_energy
   public :: sensible_heat_coeff, lw_emission_slope, leaf_transp_coeff, leaf_film_coeff, latent_pathway
   public :: tissue_hcap_per_dt
   public :: EFFAREA_HEAT, EFFAREA_EVAP, EFFAREA_TRANSP

   !----- The leaves' sidedness (docs/science/vegetation_energy_dynamics.md eq. (1)): a flat leaf    !
   !      exchanges heat from both sides, while its stomata and its film each use one. ----------!
   real(wp), parameter :: EFFAREA_HEAT   = 2.0_wp   !< [-] sensible heat
   real(wp), parameter :: EFFAREA_EVAP   = 1.0_wp   !< [-] film evaporation
   real(wp), parameter :: EFFAREA_TRANSP = 1.0_wp   !< [-] transpiration

   !----- Per-cohort leaf and wood energy-balance coefficients and the tissue heat store, held for  !
   !      the step. Each evaluation relaxes a tissue from t_leaf0/t_wood0 against hcap_per_dt =     !
   !      cap/dt, an algebraic closure rather than an unknown of the step, so it needs no error     !
   !      group and no Newton (docs/science/vegetation_energy_dynamics.md §One formulation: the     !
   !      tissue relaxes exactly). ------------------------------------------------------------!
   type :: plant_energy_env_t
      real(wp), allocatable :: h_coeff_leaf(:)   !< [W/m2/K]  frozen sensible coefficient
      real(wp), allocatable :: g_transp_leaf(:)      !< [m/s]     frozen leaf transpiration series conductance
      real(wp), allocatable :: abs_sw(:)      !< [W/m2]    absorbed shortwave (frozen source)
      real(wp), allocatable :: abs_lw(:)      !< [W/m2]    net longwave at the emission base (frozen source)
      real(wp), allocatable :: lai(:)         !< [m2/m2]   cohort leaf area index
      real(wp), allocatable :: h_coeff_w(:)   !< [W/m2/K]  frozen WOOD sensible coefficient (pi*wai*wood_gbh*rho*cp)
      real(wp), allocatable :: abs_sw_wood(:), abs_lw_wood(:) !< [W/m2] frozen absorbed SW / net LW on wood
      real(wp), allocatable :: wai(:)         !< [m2/m2]   cohort wood area index
      real(wp), allocatable :: leaf_hcap_per_dt(:), wood_hcap_per_dt(:)   !< [W/m2/K] cap/dt_fast
      real(wp), allocatable :: t_leaf0(:), t_wood0(:) !< [K]      start-of-step tissue temperatures
      !----- Heat the moving plant water brings each tissue, and the enthalpy each tissue's own water  !
      !      gains as its mass changes, both at start-of-step temperatures (docs/science/soil_plant_water.md !
      !      eqs. (4)-(5)). ----------------------------------------------------------------------------!
      real(wp), allocatable :: qwflux_wl(:)   !< [W/m2 ground] into the leaf, from sapflow (wood->leaf)
      real(wp), allocatable :: q_wood_net(:)  !< [W/m2 ground] into the wood, from root uptake less sapflow
      real(wp), allocatable :: leaf_water_enth(:)  !< [W/m2 ground] the leaf's water, eq. (5)
      real(wp), allocatable :: wood_water_enth(:)  !< [W/m2 ground] the wood's water, eq. (5)
      !----- Each cohort's transpiration over the step [kg/m2 ground/s], from the water network         !
      !      (docs/science/soil_plant_water.md §4). With transp_held the leaf energy balance holds it, so  !
      !      the canopy air gains exactly the water the leaves lose; without (the pre-pass's own call, a  !
      !      hand-built record) the demand follows the stomata and the canopy air's humidity. ----------!
      logical :: transp_held = .false.
      real(wp), allocatable :: transp_net(:)
      !----- Per-cohort longwave emissivities, the ones the radiation solver used, so each tissue's     !
      !      emission slope 4*eps*sigma*T^3 matches the longwave it was given. ----------------------!
      real(wp), allocatable :: leaf_emiss(:)  !< [-]       leaf LW emissivity
      real(wp), allocatable :: wood_emiss(:)  !< [-]       wood LW emissivity
      !----- The heat-capacity inputs behind leaf_hcap_per_dt/wood_hcap_per_dt = (dry hcap + water !
      !      mass*cp_liq)/dt (docs/science/vegetation_energy_dynamics.md eq. (8)). ----------------!
      real(wp), allocatable :: wood_dry_hcap(:)   !< [J/m2/K]  dry wood heat capacity
      real(wp), allocatable :: wood_wmass(:)      !< [kg/m2]   the wood store's water (sapwood and fine roots)
      real(wp), allocatable :: leaf_dry_hcap(:)   !< [J/m2/K]  dry leaf heat capacity
      real(wp), allocatable :: leaf_wmass(:)      !< [kg/m2]   internal (symplast) leaf water mass
   end type plant_energy_env_t

   !----- One latent pathway of a tissue's balance (docs/science/vegetation_energy_dynamics.md  !
   !      eq. (2)): the latent flux ref + slope*(T - T_cas), linearized about the canopy-air    !
   !      temperature, and h_vap, the energy per kg evaporated it was built with, which turns   !
   !      it back into mass. The default is no pathway. ----------------------------------------!
   type :: latent_pathway_t
      real(wp) :: slope = 0.0_wp            !< [W/m2/K]
      real(wp) :: ref   = 0.0_wp            !< [W/m2]   at the canopy-air temperature
      real(wp) :: h_vap = latent_heat_vap   !< [J/kg]
   end type latent_pathway_t

contains

   !---------------------------------------------------------------------------------------!
   ! Leaf or wood energy balance (docs/science/vegetation_energy_dynamics.md eqs.          !
   ! (1)-(7)): fluxes linearized about the canopy-air temperature, the tissue relaxed      !
   ! exactly over the step from t_store0 with heat capacity store_hcap_per_dt = cap/dt (0: !
   ! zero inertia). The latent flux takes two pathways: the dry share 1 - f_wet through    !
   ! the stomata (dry; none for wood), and the wetted share f_wet evaporating the film     !
   ! (wet). A tissue with no area takes t_parent: its wood's temperature for a leaf, the   !
   ! canopy air's for wood.                                                                !
   !---------------------------------------------------------------------------------------!
   elemental pure subroutine tissue_energy_balance(abs_sw, abs_lw, q_extra, h_coeff, lw_slope, dry, wet, f_wet, &
                                                   t_cas, t_parent, store_hcap_per_dt, t_store0,             &
                                                   dt_temp, t_store, transp, film_evap, dh, drnet)
      real(wp), intent(in)  :: abs_sw     !< [W/m2]   absorbed shortwave
      real(wp), intent(in)  :: abs_lw     !< [W/m2]   net longwave at the emission base
      !----- The sapflow's advected enthalpy, Q_adv of eq. (1): it shifts the balance like abs_sw but     !
      !      stays out of drnet, being an internal soil-to-plant transfer (eq. (6)). --------------------!
      real(wp), intent(in)  :: q_extra    !< [W/m2]
      real(wp), intent(in)  :: h_coeff    !< [W/m2/K] sensible-heat coefficient (sensible_heat_coeff)
      real(wp), intent(in)  :: lw_slope   !< [W/m2/K] longwave emission slope (lw_emission_slope)
      type(latent_pathway_t), intent(in) :: dry   !< the stomatal pathway (transpiration)
      type(latent_pathway_t), intent(in) :: wet   !< the film's pathway, through the boundary layer only
      real(wp), intent(in)  :: f_wet      !< [-]      wetted fraction
      real(wp), intent(in)  :: t_cas      !< [K]      canopy-air temperature, the linearization base
      real(wp), intent(in)  :: t_parent   !< [K]      the temperature a tissue with no area takes
      real(wp), intent(in)  :: store_hcap_per_dt   !< [W/m2/K] heat capacity over the step, cap/dt
      real(wp), intent(in)  :: t_store0   !< [K]      start-of-step temperature
      real(wp), intent(out) :: dt_temp    !< temperature offset from the CAS [K]
      real(wp), intent(out) :: t_store    !< diagnosed store temperature (t_cas + dt_temp) [K]
      real(wp), intent(out) :: transp     !< transpiration mass flux (0 for wood) [kg/m2/s]
      real(wp), intent(out) :: film_evap  !< film evaporation (dew if negative) [kg/m2/s]
      real(wp), intent(out) :: dh         !< sensible-flux contribution to coh_h [W/m2]
      real(wp), intent(out) :: drnet      !< net-radiation contribution to coh_rnet [W/m2]
      real(wp) :: les_dry, ler_dry, les_wet, ler_wet, denom
      real(wp) :: numer, imbalance, dt_prev, dt_avg, x, w_avg, w_flux
      les_wet = f_wet * wet%slope
      ler_wet = f_wet * wet%ref
      les_dry = (1.0_wp - f_wet) * dry%slope
      ler_dry = (1.0_wp - f_wet) * dry%ref

      !----- Coupling conductance (denom, without the tissue's own inertia) and forcing (numer) of eq. (1). !
      !      Both scale with the tissue's area, so splitting a cohort changes no temperature. ------------!
      denom   = h_coeff + les_dry + les_wet + lw_slope
      numer   = abs_sw + abs_lw + q_extra - (ler_dry + ler_wet)
      dt_prev = t_store0 - t_cas

      !=========================================================================================!
      ! Exact exponential relaxation, eqs. (3)-(7), with x = dt/tau = denom/store_hcap_per_dt:  !
      ! the endpoint (the committed state) takes the share w_avg = (1-exp(-x))/x of the         !
      ! start-of-step imbalance and the step average (every reported flux) w_flux = (1-w_avg)/x;!
      ! pairing them closes eq. (7). On the imbalance, it holds at x = 0 with no division by    !
      ! denom.                                                                                  !
      !=========================================================================================!
      if (store_hcap_per_dt > tiny_num) then
         x = denom / store_hcap_per_dt
         !----- Both weights without cancellation: their series at small x (Fortran has no expm1; a big  !
         !      cohort lives there), the 1/x tails at large x. ------------------------------------------!
         if (x < 1.0e-4_wp) then
            w_avg  = 1.0_wp - 0.5_wp*x*(1.0_wp - x/3.0_wp*(1.0_wp - 0.25_wp*x))
            w_flux = 0.5_wp*(1.0_wp - x/3.0_wp*(1.0_wp - 0.25_wp*x*(1.0_wp - 0.2_wp*x)))
         else if (x > 500.0_wp) then
            w_avg  = 1.0_wp / x
            w_flux = (1.0_wp - w_avg) / x
         else
            w_avg  = (1.0_wp - exp(-x)) / x
            w_flux = (1.0_wp - w_avg) / x
         end if
         imbalance = numer - denom * dt_prev
         dt_temp   = dt_prev + w_avg  * imbalance / store_hcap_per_dt   ! endpoint: the committed state
         dt_avg    = dt_prev + w_flux * imbalance / store_hcap_per_dt   ! average: every reported flux
      else if (denom > tiny_num) then
         dt_temp = numer / denom          ! zero heat capacity: the steady balance, at once
         dt_avg  = dt_temp
      else
         dt_temp = t_parent - t_cas       ! no area: no coupling and no heat (eq. (8))
         dt_avg  = dt_temp
      end if

      t_store   = t_cas + dt_temp
      transp    = (ler_dry + les_dry * dt_avg) / dry%h_vap
      dh        = h_coeff * dt_avg
      drnet     = abs_sw + abs_lw - lw_slope * dt_avg
      film_evap = (ler_wet + les_wet * dt_avg) / wet%h_vap
   end subroutine tissue_energy_balance

   !---------------------------------------------------------------------------------------!
   ! canopy_tissue_energy -- every cohort's wood and leaf balance at one canopy-air state, !
   ! and what the canopy as a whole gives the canopy air. Wood comes first, because a leaf !
   ! with no area takes its wood's temperature                                             !
   ! (docs/science/vegetation_energy_dynamics.md, "A tissue with no area"). A tissue pays  !
   ! the full specific enthalpy of the vapour it sheds, enthalpy_vapor at t_cas (the       !
   ! linearization base), and the canopy air receives the same number; the film pays it    !
   ! less the liquid enthalpy its water arrived with, so film, tissue and canopy air close !
   ! with no slack (eq. (2)). The latent terms take the live saturation deficit and its    !
   ! slope; the conductances are env's and film's.                                         !
   !---------------------------------------------------------------------------------------!
   pure subroutine canopy_tissue_energy(env, film, t_cas, q_cas, qsat_cas, dqsat_dt, rho, n, leaf_temp,       &
                                        wood_temp, transp, film_evap_leaf, film_evap_wood, h_canopy, qw_canopy, &
                                        transp_canopy, film_evap_canopy, rnet_canopy)
      type(plant_energy_env_t), intent(in)    :: env
      type(canopy_water_out_t), intent(in)    :: film
      real(wp),                 intent(in)    :: t_cas      !< [K]       canopy-air temperature
      real(wp),                 intent(in)    :: q_cas      !< [kg/kg]   canopy-air specific humidity
      real(wp),                 intent(in)    :: qsat_cas   !< [kg/kg]   saturation specific humidity at t_cas
      real(wp),                 intent(in)    :: dqsat_dt   !< [kg/kg/K] its temperature derivative
      real(wp),                 intent(in)    :: rho        !< [kg/m3]   canopy-air density
      integer(ik),              intent(in)    :: n
      real(wp),                 intent(inout) :: leaf_temp(:), wood_temp(:)          !< [K] each cohort's tissues
      real(wp),                 intent(inout) :: transp(:)                           !< [kg/m2 ground/s]
      real(wp),                 intent(inout) :: film_evap_leaf(:), film_evap_wood(:)   !< [kg/m2 ground/s]
      real(wp),                 intent(out)   :: h_canopy           !< [W/m2]     sensible heat into the canopy air
      real(wp),                 intent(out)   :: qw_canopy          !< [W/m2]     the vapour's enthalpy into it
      real(wp),                 intent(out)   :: transp_canopy      !< [kg/m2/s]  transpiration into it
      real(wp),                 intent(out)   :: film_evap_canopy   !< [kg/m2/s]  film evaporation into it
      real(wp),                 intent(out)   :: rnet_canopy        !< [W/m2]     net radiation the canopy absorbs
      type(latent_pathway_t) :: leaf_dry, leaf_wet, wood_dry, wood_wet
      real(wp)    :: q_def, h_evap, h_film, lw_slope, lw_slope_w
      real(wp)    :: dtl, transp_i, dh, drnet, dtw, transp_w, dh_w, drnet_w
      integer(ik) :: i
      q_def  = qsat_cas - q_cas
      h_evap = enthalpy_vapor(t_cas)
      h_film = h_evap - film%film_liquid_enthalpy
      h_canopy         = 0.0_wp
      qw_canopy        = 0.0_wp
      transp_canopy    = 0.0_wp
      rnet_canopy      = 0.0_wp
      film_evap_canopy = 0.0_wp
      do i = 1_ik, n
         if (env%transp_held) then
            !----- Transpiration held at the network's: a fixed latent term on the leaf's dry share, so  !
            !      the transpiration it reports is the water the leaves lost (soil_plant_water.md §4). -!
            leaf_dry = latent_pathway_t(slope=0.0_wp, ref=h_evap * env%transp_net(i)                    &
                                                          / max(1.0_wp - film%f_wet_c(i), tiny_num), h_vap=h_evap)
         else
            leaf_dry = latent_pathway(h_evap, rho, env%g_transp_leaf(i), dqsat_dt, q_def)
         end if
         leaf_wet = latent_pathway(h_film, rho, film%g_film_leaf(i), dqsat_dt, q_def)
         wood_dry = latent_pathway_t(h_vap=h_evap)                ! wood does not transpire
         wood_wet = latent_pathway(h_film, rho, film%g_film_w(i), dqsat_dt, q_def)
         lw_slope   = lw_emission_slope(env%leaf_emiss(i), t_cas, env%lai(i))
         lw_slope_w = lw_emission_slope(env%wood_emiss(i), t_cas, env%wai(i))
         !----- Wood: q_wood_net enters as q_extra, as qwflux_wl does for the leaf. ----------------!
         call tissue_energy_balance(env%abs_sw_wood(i), env%abs_lw_wood(i), env%q_wood_net(i), env%h_coeff_w(i), &
                                    lw_slope_w, wood_dry, wood_wet, film%f_wet_c(i), t_cas, t_cas,              &
                                    env%wood_hcap_per_dt(i), env%t_wood0(i), dtw, wood_temp(i), transp_w,       &
                                    film_evap_wood(i), dh_w, drnet_w)
         call tissue_energy_balance(env%abs_sw(i), env%abs_lw(i), env%qwflux_wl(i), env%h_coeff_leaf(i),         &
                                    lw_slope, leaf_dry, leaf_wet, film%f_wet_c(i), t_cas, wood_temp(i),          &
                                    env%leaf_hcap_per_dt(i), env%t_leaf0(i), dtl, leaf_temp(i), transp_i,       &
                                    film_evap_leaf(i), dh, drnet)
         transp(i) = transp_i
         h_canopy         = h_canopy         + dh
         qw_canopy        = qw_canopy        + (transp_i + film_evap_leaf(i)) * h_evap   ! what the leaf paid
         transp_canopy    = transp_canopy    + transp_i
         film_evap_canopy = film_evap_canopy + film_evap_leaf(i)
         rnet_canopy      = rnet_canopy      + drnet
         !----- The wood's sensible heat, film evaporation and net longwave join the leaf's. -------!
         h_canopy         = h_canopy         + dh_w
         qw_canopy        = qw_canopy        + film_evap_wood(i) * h_evap   ! what the wood paid
         film_evap_canopy = film_evap_canopy + film_evap_wood(i)
         rnet_canopy      = rnet_canopy      + drnet_w
      end do
   end subroutine canopy_tissue_energy

   !---------------------------------------------------------------------------------------!
   ! commit_tissue_energy -- each tissue's end-of-step energy and temperature, and the     !
   ! tissue store the whole-column ledger counts                                           !
   ! (docs/science/vegetation_energy_dynamics.md eqs. (9)-(10)). The energy gains the heat !
   ! cap*(t_avg - t0), t_avg being the step's time-average of where the tissue relaxed to, !
   ! and the enthalpy its water gained (docs/science/soil_plant_water.md eq. (5)); the     !
   ! temperature is read back with the end-of-step water. A masked-off process adds        !
   ! nothing; a leaf with no heat capacity takes its wood's temperature.                   !
   !---------------------------------------------------------------------------------------!
   pure subroutine commit_tissue_energy(env, nplant, dt, t_leaf_avg, t_wood_avg, n, leaf_water_mass,          &
                                        wood_water_mass, leaf_energy, wood_energy, leaf_temp, wood_temp,       &
                                        store0, store1)
      type(plant_energy_env_t), intent(in)    :: env
      real(wp),                 intent(in)    :: nplant(:)                          !< [plant/m2]
      real(wp),                 intent(in)    :: dt                                 !< [s]
      real(wp),                 intent(in)    :: t_leaf_avg(:), t_wood_avg(:)       !< [K]
      integer(ik),              intent(in)    :: n
      real(wp),                 intent(in)    :: leaf_water_mass(:), wood_water_mass(:)   !< [kg/plant] end of step
      real(wp),                 intent(inout) :: leaf_energy(:), wood_energy(:)     !< [J/m2]
      real(wp),                 intent(inout) :: leaf_temp(:), wood_temp(:)         !< [K]
      real(wp),                 intent(out)   :: store0, store1                     !< [J/m2] the tissue store
      real(wp)    :: heat_leaf, heat_wood, enth_leaf, enth_wood, w_leaf, w_wood
      integer(ik) :: i
      store0 = 0.0_wp
      store1 = 0.0_wp
      do i = 1_ik, n
         heat_leaf = env%leaf_hcap_per_dt(i) * dt * (t_leaf_avg(i) - env%t_leaf0(i))
         heat_wood = env%wood_hcap_per_dt(i) * dt * (t_wood_avg(i) - env%t_wood0(i))
         enth_leaf = env%leaf_water_enth(i) * dt
         enth_wood = env%wood_water_enth(i) * dt
         store0 = store0 + leaf_energy(i) + wood_energy(i)
         leaf_energy(i) = leaf_energy(i) + heat_leaf + enth_leaf
         wood_energy(i) = wood_energy(i) + heat_wood + enth_wood
         store1 = store1 + leaf_energy(i) + wood_energy(i)
         w_leaf = max(leaf_water_mass(i), 0.0_wp) * nplant(i)
         w_wood = max(wood_water_mass(i), 0.0_wp) * nplant(i)
         if (env%wood_dry_hcap(i) + w_wood * cp_liq > 0.0_wp)                                              &
            wood_temp(i) = tissue_energy_to_temp(wood_energy(i), w_wood, env%wood_dry_hcap(i))
         leaf_temp(i) = wood_temp(i)
         if (env%leaf_dry_hcap(i) + w_leaf * cp_liq > 0.0_wp)                                              &
            leaf_temp(i) = tissue_energy_to_temp(leaf_energy(i), w_leaf, env%leaf_dry_hcap(i))
      end do
   end subroutine commit_tissue_energy

   !---------------------------------------------------------------------------------------!
   ! The linearization coefficients of tissue_energy_balance (leaf or wood), one authority !
   ! each (docs/science/vegetation_energy_dynamics.md eqs. (1)-(2)).                       !
   !---------------------------------------------------------------------------------------!
   pure function sensible_heat_coeff(area_basis, gbh, rho, cp) result(h_coeff)
      real(wp), intent(in) :: area_basis   !< effarea_heat*LAI (leaf) or pi*WAI (wood) [m2/m2]
      real(wp), intent(in) :: gbh, rho, cp
      real(wp) :: h_coeff
      h_coeff = area_basis * gbh * rho * cp
   end function sensible_heat_coeff

   pure function leaf_transp_coeff(effarea_transp, lai, gbw, gsw) result(g_tr)
      real(wp), intent(in) :: effarea_transp, lai, gbw, gsw
      real(wp) :: g_tr
      g_tr = 0.0_wp
      if (gbw + gsw > tiny_num) g_tr = effarea_transp * lai * gbw * gsw / (gbw + gsw)
   end function leaf_transp_coeff

   !----- Film-evaporation conductance: the wet share's latent pathway, through the boundary layer   !
   !      only; same sidedness convention as leaf_transp_coeff, so leaf and wood share it. ---------!
   pure function leaf_film_coeff(effarea_evap, area_index, gbw) result(g_ev)
      real(wp), intent(in) :: effarea_evap, area_index, gbw
      real(wp) :: g_ev
      g_ev = effarea_evap * area_index * gbw
   end function leaf_film_coeff

   pure function lw_emission_slope(emiss, te, area_index) result(lw_slope)
      real(wp), intent(in) :: emiss, te, area_index
      real(wp) :: lw_slope
      lw_slope = 4.0_wp * emiss * stefan_boltzmann * te**3 * area_index
   end function lw_emission_slope

   !----- A latent pathway through conductance g (eq. (2)): h_vap*rho*g times the saturation deficit  !
   !      q_def at the canopy-air temperature (ref) and times its slope dqsat_dt (slope). -----------!
   pure function latent_pathway(h_vap, rho, g, dqsat_dt, q_def) result(path)
      real(wp), intent(in) :: h_vap      !< [J/kg]    energy per kg evaporated
      real(wp), intent(in) :: rho        !< [kg/m3]   canopy-air density
      real(wp), intent(in) :: g          !< [m/s]     the pathway's conductance
      real(wp), intent(in) :: dqsat_dt   !< [kg/kg/K]
      real(wp), intent(in) :: q_def      !< [kg/kg]   qsat(T_cas) - q_cas
      type(latent_pathway_t) :: path
      path%slope = h_vap * rho * g * dqsat_dt
      path%ref   = h_vap * rho * g * q_def
      path%h_vap = h_vap
   end function latent_pathway

   !----- A tissue's heat capacity over the step, (dry heat capacity + water mass * c_liq)/dt (eq.  !
   !      (8)). ----------------------------------------------------------------------------------!
   elemental pure function tissue_hcap_per_dt(dry_hcap, water_mass, dt) result(hcap_per_dt)
      real(wp), intent(in) :: dry_hcap     !< [J/m2/K]
      real(wp), intent(in) :: water_mass   !< [kg/m2]
      real(wp), intent(in) :: dt           !< [s]
      real(wp) :: hcap_per_dt
      hcap_per_dt = (dry_hcap + water_mass * cp_liq) / dt
   end function tissue_hcap_per_dt

end module meds_plant_energy
