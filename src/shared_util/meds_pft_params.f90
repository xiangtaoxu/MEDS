!==========================================================================================!
! meds_pft_params -- plant functional type (PFT) trait table.                              !
!                                                                                          !
! Structure-of-arrays: one allocatable array per trait, indexed by PFT. Size allometry is    !
! pan-tropical and PFT-independent (see meds_allometry); the ONLY allometric per-PFT input   !
! is `dbh_critical` (the maximum diameter) and `wood_density` (rho, which enters AGB).           !
!                                                                                          !
! The growth/mortality contrast between PFTs collapses onto ONE physical axis: WOOD DENSITY. !
! Following the classic ED / Camac-2018 trade-off, low wood density => high maximum relative  !
! growth rate, but also high baseline mortality AND high extra mortality under shade; high    !
! wood density => slow but tolerant. The derived rates are powers of (rho_ref/rho):           !
!                                                                                          !
!   gr_max     = gr_ref  * (rho_ref/rho)^k_grow     [1/yr]  max relative DBH growth           !
!   mort_base  = m_ref   * (rho_ref/rho)^k_mort      [1/yr]  density-independent baseline       !
!   mort_shade = ms_ref  * (rho_ref/rho)^k_shade     [1/yr]  extra mortality at full shade      !
!                                                                                          !
! All PFTs share `min_reproduction_height` (recruits are born at this height) and the same    !
! `recruit_dens` (identical recruitment output), so PFTs differ ONLY through wood density and  !
! their maximum diameter. The default table seeds three strategies (pioneer/mid/climax).       !
!==========================================================================================!
module meds_pft_params
   use meds_kinds, only : wp, ik
   implicit none
   private

   public :: pft_table_t, init_default_pfts

   !---------------------------------------------------------------------------------------!
   ! PFT trait table (SoA).  Units in brackets.                                            !
   !---------------------------------------------------------------------------------------!
   type :: pft_table_t
      integer(ik) :: n = 0_ik
      !----- Size limits (allometry itself is global, see meds_allometry). ----------------!
      real(wp), allocatable :: dbh_critical(:)       !< [cm]    maximum diameter (growth clamp)
      real(wp), allocatable :: wood_density(:)   !< [g/cm3] rho: the growth-mortality anchor
      !----- Wood-density-derived vital-rate parameters. ----------------------------------!
      real(wp), allocatable :: gr_max(:)         !< [1/yr]  maximum relative DBH growth rate
      real(wp), allocatable :: mort_base(:)      !< [1/yr]  density-independent baseline mortality
      real(wp), allocatable :: mort_shade(:)     !< [1/yr]  extra mortality at full shade
      !----- Recruitment (identical across PFTs by design). -------------------------------!
      real(wp),    allocatable :: recruit_dens(:)     !< [plant/m2/month] potential recruit density
      real(wp),    allocatable :: recruit_min_temp(:) !< [K] minimum monthly temp to recruit
      integer(ik), allocatable :: include_pft(:)      !< 1 = PFT may recruit, 0 = excluded
      !----- Shared scalar: recruit birth height (same for all PFTs). ---------------------!
      real(wp) :: min_reproduction_height = 2.0_wp    !< [m] height at which recruits are born
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
      allocate(pft%gr_max(n), pft%mort_base(n), pft%mort_shade(n))
      allocate(pft%recruit_dens(n), pft%recruit_min_temp(n), pft%include_pft(n))
   end subroutine alloc_pft_table

   !---------------------------------------------------------------------------------------!
   ! Seed three contrasting PFTs:  1 = pioneer, 2 = mid-successional, 3 = climax,           !
   ! distinguished ONLY by wood density (low->high) and maximum diameter.                   !
   !---------------------------------------------------------------------------------------!
   subroutine init_default_pfts(pft)
      type(pft_table_t), intent(out) :: pft
      !----- Trade-off anchors (reference = mid-successional). ----------------------------!
      real(wp), parameter :: rho_ref = 0.60_wp
      real(wp), parameter :: gr_ref  = 0.15_wp,  k_grow  = 1.5_wp   ! [1/yr]
      real(wp), parameter :: m_ref   = 0.03_wp,  k_mort  = 1.5_wp   ! [1/yr]
      real(wp), parameter :: ms_ref  = 0.15_wp,  k_shade = 1.5_wp   ! [1/yr]

      call alloc_pft_table(pft, 3_ik)

      !----- The single free axis: wood density. ------------------------------------------!
      pft%wood_density = [ 0.40_wp, 0.60_wp, 0.85_wp ]
      pft%dbh_critical     = [ 40.0_wp, 80.0_wp, 120.0_wp ]

      !----- Growth-mortality trade-off, all powers of (rho_ref/rho). ---------------------!
      pft%gr_max     = gr_ref * (rho_ref / pft%wood_density) ** k_grow
      pft%mort_base  = m_ref  * (rho_ref / pft%wood_density) ** k_mort
      pft%mort_shade = ms_ref * (rho_ref / pft%wood_density) ** k_shade

      !----- Recruitment: identical output across PFTs (shared birth height in the type). -!
      pft%recruit_dens           = [ 0.015_wp, 0.015_wp, 0.015_wp ]
      pft%recruit_min_temp       = [ 283.15_wp, 283.15_wp, 283.15_wp ]
      pft%include_pft            = [ 1_ik, 1_ik, 1_ik ]
      pft%min_reproduction_height = 2.0_wp
   end subroutine init_default_pfts

end module meds_pft_params
