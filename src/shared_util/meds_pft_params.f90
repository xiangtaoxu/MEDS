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

   public :: pft_table_t, alloc_pft_table, derive_pft_rates

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
      !----- Wood-density -> mortality-hazard DERIVATION coefficients (Camac et al. 2018 PNAS,  !
      !       a power law centred on mort_rho_ref: param = param_0 * (rho/rho_ref)^exp, always   !
      !       positive). Shared scalars for now (PFT-specific later). --------------------------!
      real(wp) :: mort_rho_ref            !< [g/cm3] wood-density reference for the power law
      real(wp) :: mort_gamma_0, mort_gamma_exp   !< baseline hazard scale + exponent
      real(wp) :: mort_alpha_0, mort_alpha_exp   !< low-growth hazard scale + exponent
      real(wp) :: mort_beta_0,  mort_beta_exp    !< growth-sensitivity scale + exponent
      !----- Wood-density-derived mortality-hazard parameters (DERIVED from the above). -------!
      real(wp), allocatable :: mort_gamma(:)     !< [1/yr]  growth-independent baseline hazard
      real(wp), allocatable :: mort_alpha(:)     !< [1/yr]  low-growth hazard magnitude
      real(wp), allocatable :: mort_beta(:)      !< [yr/cm] growth sensitivity of the hazard
      !----- Recruitment. -----------------------------------------------------------------!
      real(wp),    allocatable :: seed_rain_recruits(:) !< [plant/m2/yr] baseline external seed rain
      integer(ik), allocatable :: include_pft(:)        !< 1 = PFT may recruit, 0 = excluded
      !----- Shared height thresholds. ----------------------------------------------------!
      real(wp) :: min_cohort_height          !< [m] smallest tracked cohort; recruits born here
      real(wp) :: min_reproduction_height    !< [m] height a cohort must exceed to reproduce
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
   ! Derive the Camac (2018) mortality-hazard parameters from wood density as a power law      !
   ! centred on mort_rho_ref (param = param_0 * (rho/rho_ref)^exp). Always positive. Call after !
   ! wood_density and the Camac coefficients are set (e.g. from the PFT config file).           !
   !---------------------------------------------------------------------------------------!
   subroutine derive_pft_rates(pft)
      type(pft_table_t), intent(inout) :: pft
      pft%mort_gamma = pft%mort_gamma_0 * (pft%wood_density / pft%mort_rho_ref) ** pft%mort_gamma_exp
      pft%mort_alpha = pft%mort_alpha_0 * (pft%wood_density / pft%mort_rho_ref) ** pft%mort_alpha_exp
      pft%mort_beta  = pft%mort_beta_0  * (pft%wood_density / pft%mort_rho_ref) ** pft%mort_beta_exp
   end subroutine derive_pft_rates

end module meds_pft_params
