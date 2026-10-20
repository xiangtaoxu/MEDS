! SPDX-License-Identifier: Apache-2.0
!==========================================================================================!
! meds_sites_opts -- the [[sites]] tables, as a low-level config leaf: one record per site  !
! of the polygon (docs/configuration.md, "Several sites"). A site's surface is built from  !
! it at setup (make_site_surface, through site_config); the record itself is never state.  !
!==========================================================================================!
module meds_sites_opts
   use meds_kinds,              only : wp, ik
   use meds_soil_column_config, only : soil_column_config_t
   implicit none
   private

   public :: site_spec_t

   !----- One site: its share of the polygon, its terrain, the census rows it reads and, when its !
   !      table sets any of [soil]'s column keys, its own ground, [soil]'s values for each key it  !
   !      leaves out. -----------------------------------------------------------------------------!
   type :: site_spec_t
      real(wp)    :: area_fraction  = 0.0_wp      !< [-] its share of the polygon
      logical     :: has_elevation  = .false.     !< .false.: the polygon's elevation
      real(wp)    :: elevation      = 0.0_wp      !< [m]
      real(wp)    :: slope          = 0.0_wp      !< [degrees] from horizontal
      real(wp)    :: aspect         = 0.0_wp      !< [degrees] clockwise from north
      integer(ik) :: census_site_id = 0_ik        !< the census rows it reads; 0: not given
      logical     :: own_ground     = .false.     !< .false.: the polygon's ground
      type(soil_column_config_t) :: soil_column
      integer(ik) :: water_bottom_bc = 0_ik       !< SOIL_BC_*
      real(wp)    :: drain_factor    = 1.0_wp     !< [-]
      real(wp)    :: deep_temp       = 0.0_wp     !< [K]
      real(wp)    :: deep_depth      = 0.0_wp     !< [m]
   end type site_spec_t

end module meds_sites_opts
