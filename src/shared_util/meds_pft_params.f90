!==========================================================================================!
! meds_pft_params -- plant functional type (PFT) trait table and allometry.                !
!                                                                                          !
! Structure-of-arrays: one allocatable array per trait, indexed by PFT. Only the traits    !
! relevant to PURE DEMOGRAPHY are kept (no carbon pools, no leaf area). Allometry is a      !
! monotone, analytically invertible diameter->height curve so that fused/derived cohorts    !
! recover a consistent height from diameter.                                               !
!                                                                                          !
!   h(dbh)  = height_min + b1_height * (1 - exp(-b2_height*dbh))          (asymptote height_min+b1_height)       !
!   dbh(h)  = -ln(1 - (h-height_min)/b1_height) / b2_height                                              !
!                                                                                          !
! The default table seeds three contrasting strategies (pioneer / mid / climax) so the     !
! demographic engine exhibits succession and self-thinning out of the box.                  !
!==========================================================================================!
module meds_pft_params
   use meds_kinds,     only : wp, ik
   use meds_constants, only : tiny_num
   implicit none
   private

   public :: pft_table_t, init_default_pfts, dbh_to_height, height_to_dbh

   !---------------------------------------------------------------------------------------!
   ! PFT trait table (SoA).  Units in brackets.                                            !
   !---------------------------------------------------------------------------------------!
   type :: pft_table_t
      integer(ik) :: n = 0_ik
      !----- Allometry / size limits. -----------------------------------------------------!
      real(wp), allocatable :: dbh_min(:)        !< [cm]  recruit size & lower clamp
      real(wp), allocatable :: dbh_crit(:)       !< [cm]  asymptotic max diameter (upper clamp)
      real(wp), allocatable :: height_min(:)        !< [m]   height at dbh=0
      real(wp), allocatable :: height_max(:)        !< [m]   asymptotic height
      real(wp), allocatable :: b1_height(:)           !< [m]   height scale (height_max-height_min)
      real(wp), allocatable :: b2_height(:)           !< [1/cm] height curvature (>0)
      !----- Empirical growth: g = g0*(1-dbh/dbh_crit)*exp(-gcomp*comp). -------------------!
      real(wp), allocatable :: g0(:)             !< [cm/yr] potential growth, small & open
      real(wp), allocatable :: gcomp(:)          !< [m2/m2]^-1 competition sensitivity
      !----- Empirical mortality (ED2 "scheme 2" analog). ---------------------------------!
      real(wp), allocatable :: mort1(:)          !< [1/yr] growth-dependent amplitude
      real(wp), allocatable :: mort2(:)          !< [yr/cm] growth-dependence decay
      real(wp), allocatable :: mort3(:)          !< [1/yr] density-independent background
      real(wp), allocatable :: frost_mort(:)     !< [1/yr] frost mortality scale
      real(wp), allocatable :: plant_min_temp(:) !< [K]    frost onset temperature
      !----- Size-dependent treefall mortality. -------------------------------------------!
      real(wp), allocatable :: treefall_dbh(:)   !< [cm]   threshold diameter
      real(wp), allocatable :: treefall_rate(:)  !< [1/yr] treefall mortality above threshold
      !----- Empirical recruitment. -------------------------------------------------------!
      real(wp),    allocatable :: recruit_dens(:)     !< [plant/m2/month] potential recruit density
      real(wp),    allocatable :: recruit_min_temp(:) !< [K] minimum monthly temp to recruit
      integer(ik), allocatable :: include_pft(:)      !< 1 = PFT may recruit, 0 = excluded
   end type pft_table_t

contains

   !---------------------------------------------------------------------------------------!
   ! Allocate every trait array to n PFTs.                                                 !
   !---------------------------------------------------------------------------------------!
   subroutine alloc_pft_table(t, n)
      type(pft_table_t), intent(inout) :: t
      integer(ik),       intent(in)    :: n
      t%n = n
      allocate(t%dbh_min(n), t%dbh_crit(n), t%height_min(n), t%height_max(n), t%b1_height(n), t%b2_height(n))
      allocate(t%g0(n), t%gcomp(n))
      allocate(t%mort1(n), t%mort2(n), t%mort3(n), t%frost_mort(n), t%plant_min_temp(n))
      allocate(t%treefall_dbh(n), t%treefall_rate(n))
      allocate(t%recruit_dens(n), t%recruit_min_temp(n), t%include_pft(n))
   end subroutine alloc_pft_table

   !---------------------------------------------------------------------------------------!
   ! Seed three contrasting PFTs:  1 = pioneer, 2 = mid-successional, 3 = climax.          !
   !---------------------------------------------------------------------------------------!
   subroutine init_default_pfts(t)
      type(pft_table_t), intent(out) :: t
      call alloc_pft_table(t, 3_ik)

      !----- Size / allometry. ------------------------------------------------------------!
      t%dbh_min  = [  1.0_wp,  1.0_wp,  1.0_wp ]
      t%dbh_crit = [ 40.0_wp, 80.0_wp,120.0_wp ]
      t%height_min  = [  0.5_wp,  0.5_wp,  0.5_wp ]
      t%height_max  = [ 18.0_wp, 30.0_wp, 40.0_wp ]
      t%b1_height     = t%height_max - t%height_min
      t%b2_height     = [ 0.060_wp, 0.040_wp, 0.030_wp ]

      !----- Growth: pioneers fast & competition-sensitive, climax slow & tolerant. -------!
      t%g0       = [ 2.5_wp, 1.4_wp, 0.8_wp ]
      t%gcomp    = [ 1.2_wp, 0.6_wp, 0.25_wp ]

      !----- Mortality: pioneers high background, climax low; all growth-dependent. -------!
      t%mort1    = [ 1.5_wp, 1.0_wp, 0.6_wp ]
      t%mort2    = [ 5.0_wp, 5.0_wp, 5.0_wp ]
      t%mort3    = [ 0.06_wp, 0.02_wp, 0.008_wp ]
      t%frost_mort     = [ 3.0_wp, 3.0_wp, 3.0_wp ]
      t%plant_min_temp = [ 278.15_wp, 278.15_wp, 278.15_wp ]

      !----- Treefall (large trees). ------------------------------------------------------!
      t%treefall_dbh  = [ 25.0_wp, 40.0_wp, 60.0_wp ]
      t%treefall_rate = [ 0.014_wp, 0.014_wp, 0.014_wp ]

      !----- Recruitment. -----------------------------------------------------------------!
      t%recruit_dens     = [ 0.30_wp, 0.15_wp, 0.06_wp ]
      t%recruit_min_temp = [ 283.15_wp, 283.15_wp, 283.15_wp ]
      t%include_pft      = [ 1_ik, 1_ik, 1_ik ]
   end subroutine init_default_pfts

   !---------------------------------------------------------------------------------------!
   ! Diameter -> height (elemental, pure).                                                 !
   !---------------------------------------------------------------------------------------!
   elemental pure function dbh_to_height(height_min, b1_height, b2_height, dbh) result(h)
      real(wp), intent(in) :: height_min, b1_height, b2_height, dbh
      real(wp)             :: h
      h = height_min + b1_height * (1.0_wp - exp(-b2_height * dbh))
   end function dbh_to_height

   !---------------------------------------------------------------------------------------!
   ! Height -> diameter (inverse of dbh_to_height, guarded against the asymptote).                 !
   !---------------------------------------------------------------------------------------!
   elemental pure function height_to_dbh(height_min, b1_height, b2_height, h) result(dbh)
      real(wp), intent(in) :: height_min, b1_height, b2_height, h
      real(wp)             :: dbh, frac
      frac = max(1.0_wp - (h - height_min) / b1_height, tiny_num)
      dbh  = -log(frac) / b2_height
   end function height_to_dbh

end module meds_pft_params
