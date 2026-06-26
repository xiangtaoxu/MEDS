!==========================================================================================!
! meds_pft_params -- plant functional type (PFT) trait table.                              !
!                                                                                          !
! Structure-of-arrays: one allocatable array per trait, indexed by PFT. Size allometry is    !
! pan-tropical and PFT-independent (see meds_allometry); the per-PFT allometric inputs are    !
! `dbh_critical` (the maximum diameter) and `wood_density` (rho, which enters AGB).            !
!                                                                                          !
! The PFTs carry the parameters of the PHENOMENOLOGICAL vital rates                            !
! (meds_phenomenological_vital_rates):                                                        !
!   * GROWTH: an intrinsic capped log-linear function of dbh (growth_dbh_slope/cap/max),       !
!     suppressed multiplicatively by neighbourhood competition (growth_lai_slope on overtopping!
!     LAI) and by reproductive allocation (reproduction_investment_fraction above maturity).   !
!   * MORTALITY: the Camac et al. (2018, PNAS) additive hazard, simplified to                  !
!     rate = mort_gamma + mort_alpha*exp(-mort_beta*growth_avg), with the three parameters     !
!     DERIVED from wood density (low rho => higher baseline & low-growth hazard).              !
!   * RECRUITMENT: a baseline external `seed_rain_recruits` plus a reproduction flux computed  !
!     in the rate module from the carbon diverted to reproduction.                             !
!                                                                                          !
! Two shared height thresholds: `min_cohort_height` (the smallest tracked cohort; recruits are !
! born here) and `min_reproduction_height` (the height a cohort must exceed to reproduce). The !
! default table seeds three strategies (pioneer/mid/climax) differing only in wood density and !
! maximum diameter.                                                                          !
!==========================================================================================!
module meds_pft_params
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: pft_table_t, init_default_pfts, derive_pft_rates

   !----- Wood-density -> mortality-hazard parameters (Camac et al. 2018 PNAS). A POWER LAW in --!
   !      wood density centred on mort_rho_ref: param = param_0 * (rho/rho_ref)^exp, equivalently!
   !      exp(log(param_0) + (log(rho) - log(rho_ref))*exp). The power-law form keeps every       !
   !      parameter strictly positive for any rho > 0 (no clamping needed). The hazard is         !
   !      rate = gamma + alpha*exp(-beta*growth_avg) [1/yr], growth_avg in cm/yr.                 !
   real(wp), parameter :: mort_rho_ref   = 0.6_wp
   real(wp), parameter :: mort_gamma_0   = 0.0094_wp,  mort_gamma_exp = -1.8392_wp
   real(wp), parameter :: mort_alpha_0   = 0.05_wp,    mort_alpha_exp = -1.1493_wp
   real(wp), parameter :: mort_beta_0    = 18.72_wp,   mort_beta_exp  =  0.2792_wp

   !---------------------------------------------------------------------------------------!
   ! PFT trait table (SoA).  Units in brackets.                                            !
   !---------------------------------------------------------------------------------------!
   type :: pft_table_t
      integer(ik) :: n = 0_ik
      !----- Size limits + the wood-density axis (allometry itself is global). ------------!
      real(wp), allocatable :: dbh_critical(:)   !< [cm]    maximum diameter (growth clamp)
      real(wp), allocatable :: wood_density(:)   !< [g/cm3] rho: AGB + mortality anchor
      !----- Intrinsic growth: capped log-linear in dbh. ----------------------------------!
      real(wp), allocatable :: growth_dbh_slope(:) !< [--]    slope of ln(growth) vs (ln cap - ln dbh)
      real(wp), allocatable :: growth_dbh_cap(:)   !< [cm]    dbh at/above which growth floors at the max
      real(wp), allocatable :: growth_dbh_max(:)   !< [cm/yr] intrinsic growth at/above the cap
      !----- Growth suppression: competition (overtopping LAI) + reproductive allocation. -!
      real(wp), allocatable :: growth_lai_slope(:) !< [(m2/m2)^-1] slope of ln(suppression) vs overtopping LAI
      real(wp), allocatable :: reproduction_investment_fraction(:) !< [--] growth fraction diverted to reproduction
      real(wp), allocatable :: repro_carbon_efficiency(:)          !< [--] reproduction carbon -> establishable recruits
      !----- Wood-density-derived mortality-hazard parameters (Camac 2018; power law in rho). !
      real(wp), allocatable :: mort_gamma(:)     !< [1/yr]  growth-independent baseline hazard
      real(wp), allocatable :: mort_alpha(:)     !< [1/yr]  low-growth hazard magnitude
      real(wp), allocatable :: mort_beta(:)      !< [yr/cm] growth sensitivity of the hazard
      !----- Recruitment. -----------------------------------------------------------------!
      real(wp),    allocatable :: seed_rain_recruits(:) !< [plant/m2/yr] baseline external seed rain
      integer(ik), allocatable :: include_pft(:)        !< 1 = PFT may recruit, 0 = excluded
      !----- Shared height thresholds. ----------------------------------------------------!
      real(wp) :: min_cohort_height       = 2.0_wp   !< [m] smallest tracked cohort; recruits born here
      real(wp) :: min_reproduction_height = 20.0_wp  !< [m] height a cohort must exceed to reproduce
   end type pft_table_t

contains

   !---------------------------------------------------------------------------------------!
   ! Allocate every trait array to n PFTs.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_pft_table(pft, n)
      type(pft_table_t), intent(inout) :: pft
      integer(ik),       intent(in)    :: n
      pft%n = n
      allocate(pft%dbh_critical(n), pft%wood_density(n))
      allocate(pft%growth_dbh_slope(n), pft%growth_dbh_cap(n), pft%growth_dbh_max(n))
      allocate(pft%growth_lai_slope(n), pft%reproduction_investment_fraction(n),              &
               pft%repro_carbon_efficiency(n))
      allocate(pft%mort_gamma(n), pft%mort_alpha(n), pft%mort_beta(n))
      allocate(pft%seed_rain_recruits(n), pft%include_pft(n))
   end subroutine alloc_pft_table

   !---------------------------------------------------------------------------------------!
   ! Seed three contrasting PFTs:  1 = pioneer, 2 = mid-successional, 3 = climax,           !
   ! distinguished by wood density (low->high) and maximum diameter. The growth/reproduction!
   ! parameters are PFT-specific but seeded identically for now (tune per PFT later).        !
   !---------------------------------------------------------------------------------------!
   subroutine init_default_pfts(pft)
      type(pft_table_t), intent(out) :: pft

      call alloc_pft_table(pft, 3_ik)

      !----- The wood-density axis + maximum diameter. ------------------------------------!
      pft%wood_density = [ 0.40_wp, 0.60_wp, 0.85_wp ]
      pft%dbh_critical = [ 40.0_wp, 80.0_wp, 120.0_wp ]

      !----- Intrinsic growth: shape uniform for now, but the max rate falls pioneer->climax. !
      pft%growth_dbh_slope = [ 0.25_wp,  0.25_wp,  0.25_wp ]
      pft%growth_dbh_cap   = [ 100.0_wp, 100.0_wp, 100.0_wp ]
      pft%growth_dbh_max   = [ 1.5_wp,   1.0_wp,   0.5_wp ]

      !----- Competition + reproduction suppression (same for all PFTs for now). ----------!
      pft%growth_lai_slope                 = [ -0.5_wp, -0.5_wp, -0.5_wp ]
      pft%reproduction_investment_fraction = [  0.3_wp,  0.3_wp,  0.3_wp ]
      pft%repro_carbon_efficiency          = [  1.0e-3_wp, 1.0e-3_wp, 1.0e-3_wp ]

      !----- Recruitment: baseline seed rain (identical across PFTs by default). -----------!
      pft%seed_rain_recruits = [ 0.01_wp, 0.01_wp, 0.01_wp ]
      pft%include_pft        = [ 1_ik, 1_ik, 1_ik ]

      pft%min_cohort_height       = 2.0_wp
      pft%min_reproduction_height = 20.0_wp

      !----- Derive the wood-density-dependent mortality-hazard parameters. ---------------!
      call derive_pft_rates(pft)
   end subroutine init_default_pfts

   !---------------------------------------------------------------------------------------!
   ! Derive the Camac (2018) mortality-hazard parameters from wood density as a power law      !
   ! centred on mort_rho_ref (param = param_0 * (rho/rho_ref)^exp). Always positive. Call after !
   ! wood_density is set or changed (e.g. by a meds_config.toml override).                     !
   !---------------------------------------------------------------------------------------!
   subroutine derive_pft_rates(pft)
      type(pft_table_t), intent(inout) :: pft
      pft%mort_gamma = mort_gamma_0 * (pft%wood_density / mort_rho_ref) ** mort_gamma_exp
      pft%mort_alpha = mort_alpha_0 * (pft%wood_density / mort_rho_ref) ** mort_alpha_exp
      pft%mort_beta  = mort_beta_0  * (pft%wood_density / mort_rho_ref) ** mort_beta_exp
   end subroutine derive_pft_rates

end module meds_pft_params
