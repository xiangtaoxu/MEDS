!==========================================================================================!
! meds_demography_capi -- the ISO_C_BINDING shim that exposes the DEMOGRAPHIC model to C /     !
! Python (ctypes), compiled ONLY into the optional shared library libmeds_c (MEDS_BUILD_PYLIB). !
!                                                                                          !
! `site_t` (allocatable-component derived type -- not f2py-able) crosses as an OPAQUE HANDLE:   !
! a module-`save` registry of live sites indexed by a small integer the caller holds. Config    !
! is likewise an opaque handle loaded from the same TOML `meds_main` uses. Only flat arrays      !
! (copied into caller buffers) and scalars cross the boundary; getters always COPY.             !
!                                                                                          !
! This first cut drives the standalone CARBON model from Python: load config -> create site ->   !
! init bare -> advance_slow (the driver's carbon orchestration) -> read state -> free. A         !
! per-handle GENERATION counter (bumped every advance_slow, which reorders the SoA on fuse/fiss) !
! lets the caller detect a stale positional snapshot -- global_id is the only stable key.        !
!==========================================================================================!
module meds_demography_capi
   use iso_c_binding
   use meds_kinds,                  only : wp, ik
   use meds_constants,              only : pio4, tiny_num
   use meds_config,                 only : meds_config_t, growth_window_steps
   use meds_config_io,              only : load_meds_config
   use meds_core_state_types,        only : site_t, site_free, cohort_deriv_block, cohort_deriv_alloc
   use meds_init,                   only : init_bare_ground
   use meds_vegetation_dynamics,    only : vegetation_dynamics
   use meds_output_diagnostics, only : total_agb, total_lai, total_nplant, total_basal_area, count_cohorts
   use meds_allometry,              only : b1Ht, b2Ht, agb_c1, agb_c2, lai_b1, lai_b2
   use meds_core_interface,   only : update_cohort_states, apply_recruitment,               &
                                           apply_patch_disturbance, new_fuse_cohorts,             &
                                           terminate_cohorts, split_cohorts, new_fuse_patches,    &
                                           terminate_patches, sort_cohorts, sort_patches,         &
                                           update_overtopping_lai
   implicit none
   private

   integer, parameter :: MAXH = 8
   type(site_t),        target, save :: g_site(MAXH)
   type(meds_config_t), target, save :: g_cfg(MAXH)
   logical,             save :: site_used(MAXH) = .false.
   logical,             save :: cfg_used(MAXH)  = .false.
   integer(c_long),     save :: g_generation(MAXH) = 0_c_long

contains

   !----- Build a Fortran string from an explicit-length C char buffer (no strlen needed). ---!
   pure function f_string(cbuf, n) result(s)
      character(kind=c_char), intent(in) :: cbuf(*)
      integer(c_int),         intent(in) :: n
      character(len=n) :: s
      integer :: i
      do i = 1, n
         s(i:i) = cbuf(i)
      end do
   end function f_string

   !----- Load a config from a TOML path -> opaque cfg handle (>=1), or -1 if the registry is full. -!
   function meds_config_load(path, path_len) result(h) bind(c, name="meds_config_load")
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value, intent(in)  :: path_len
      integer(c_int)                     :: h
      integer :: i
      h = -1_c_int
      do i = 1, MAXH
         if (.not. cfg_used(i)) then ; h = int(i, c_int) ; exit ; end if
      end do
      if (h < 0_c_int) return
      call load_meds_config(f_string(path, path_len), g_cfg(h))
      cfg_used(h) = .true.
   end function meds_config_load

   !----- Create an empty site slot -> opaque site handle (>=1), or -1 if full. --------------!
   function meds_site_create() result(h) bind(c, name="meds_site_create")
      integer(c_int) :: h
      integer :: i
      h = -1_c_int
      do i = 1, MAXH
         if (.not. site_used(i)) then
            h = int(i, c_int) ; site_used(i) = .true. ; g_generation(i) = 0_c_long ; exit
         end if
      end do
   end function meds_site_create

   subroutine meds_site_init_bare(sh, ch, n_patch) bind(c, name="meds_site_init_bare")
      integer(c_int), value, intent(in) :: sh, ch, n_patch
      call init_bare_ground(g_site(sh), g_cfg(ch), int(n_patch, ik))
      g_generation(sh) = g_generation(sh) + 1_c_long
   end subroutine meds_site_init_bare

   !----- Advance one slow step (the driver's CARBON orchestration + cadence). Bumps the        !
   !       generation (the SoA is reordered by fuse/fission). is_new_month/is_new_year are 0/1.  !
   subroutine meds_advance_slow(sh, ch, is_new_month, is_new_year) bind(c, name="meds_advance_slow")
      integer(c_int), value, intent(in) :: sh, ch, is_new_month, is_new_year
      call vegetation_dynamics(g_site(sh), g_cfg(ch), is_new_month /= 0_c_int, is_new_year /= 0_c_int)
      g_generation(sh) = g_generation(sh) + 1_c_long
   end subroutine meds_advance_slow

   !---------------------------------------------------------------------------------------!
   ! Apply CALLER-SUPPLIED demographic rates (the empirical / Python path): grow by the       !
   ! per-cohort growth [cm/yr], die by the per-cohort mortality [1/yr], recruit from the       !
   ! per-(PFT,patch) recruitment [plant/m2/yr], then age patches + run the fuse/fission        !
   ! cadence -- the engine's law-free apply-primitives, sequenced exactly as the former         !
   ! empirical update_demography. `recr` is a flat column-major (PFT-fastest) npft*npatch array. !
   ! is_new_month / is_new_year are 0/1.                                                        !
   !---------------------------------------------------------------------------------------!
   subroutine meds_apply_rates(sh, ch, growth, mortality, recr, is_new_month, is_new_year) &
                               bind(c, name="meds_apply_rates")
      integer(c_int), value, intent(in) :: sh, ch, is_new_month, is_new_year
      real(c_double),        intent(in) :: growth(*), mortality(*), recr(*)
      integer(ik) :: n, np, npft, ip, pf, n_window
      real(wp), allocatable :: g(:), m(:), rec(:,:)
      type(cohort_deriv_block) :: deriv
      logical :: do_cohort_fissfuse, do_patch_disturbance, do_patch_fissfuse
      real(wp), parameter :: PATCH_DYNAMICS_INTERVAL = 1.0_wp

      associate (site => g_site(sh), cfg => g_cfg(ch))
         n = site%cohort%n ; np = site%patch%n ; npft = cfg%pft%n
         allocate(g(max(n,1)), m(max(n,1)), rec(npft, max(np,1)))
         if (n > 0_ik) then
            g(1:n) = real(growth(1:n),    wp)
            m(1:n) = real(mortality(1:n), wp)
         end if
         do ip = 1_ik, np
            do pf = 1_ik, npft
               rec(pf, ip) = real(recr((ip - 1_ik) * npft + pf), wp)
            end do
         end do

         do_cohort_fissfuse   = is_new_month /= 0_c_int .and. cfg%demography_on .and. cfg%do_cohort_fissfuse
         do_patch_disturbance = is_new_year  /= 0_c_int .and. cfg%demography_on .and. cfg%do_patch_disturbance
         do_patch_fissfuse    = is_new_year  /= 0_c_int .and. cfg%demography_on .and. cfg%do_patch_fissfuse

         n_window = growth_window_steps(cfg)
         site%growth_hist_pos = mod(site%growth_hist_pos, n_window) + 1_ik
         !----- Empirical growth + mortality via the TENDENCY seam: compute the per-cohort         !
         !      derivatives from the supplied rates (forward allometry + the supplied-rate ring      !
         !      buffer), then let the core engine's pure applier advance the state. --------------!
         call cohort_deriv_alloc(deriv, site%cohort%n)
         call compute_empirical_derivatives(site, g, m, cfg%dt_years, n_window, site%growth_hist_pos, deriv)
         call update_cohort_states(site%cohort%n, site%cohort%dbh, site%cohort%height,              &
                                   site%cohort%basal_area, site%cohort%agb, site%cohort%leaf_area,  &
                                   site%cohort%leaf_carbon, site%cohort%fineroot_carbon,            &
                                   site%cohort%wood_carbon, site%cohort%nonstructural_carbon,       &
                                   site%cohort%nplant, deriv%d_dbh_dt, deriv%d_height_dt,           &
                                   deriv%d_basal_area_dt, deriv%d_agb_dt, deriv%d_leaf_area_dt,     &
                                   deriv%d_leaf_carbon_dt, deriv%d_fineroot_carbon_dt,              &
                                   deriv%d_wood_carbon_dt, deriv%d_nonstructural_carbon_dt,         &
                                   deriv%dln_nplant_dt, cfg%dt_years, cfg%negligible_nplant)
         !----- NOTE: this deliberately mirrors the ORIGINAL empirical update_demography order      !
         !      (sort only on the monthly/annual fuse-fiss cadence, below) so the Python empirical   !
         !      example reproduces the historical golden. The go-forward CARBON model                !
         !      (meds_vegetation_dynamics) re-sorts every step -- the #2 correctness fix.             !
         do ip = 1_ik, site%patch%n
            site%patch%age(ip) = site%patch%age(ip) + cfg%dt_years
         end do
         if (do_cohort_fissfuse) then
            call apply_recruitment(site, cfg, rec)
            call new_fuse_cohorts(site, cfg) ; call terminate_cohorts(site, cfg)
            call split_cohorts(site, cfg)    ; call sort_cohorts(site)
         end if
         if (do_patch_disturbance) call apply_patch_disturbance(site, cfg, PATCH_DYNAMICS_INTERVAL)
         if (do_patch_fissfuse) then
            call sort_patches(site)      ; call new_fuse_patches(site, cfg)
            call terminate_patches(site, cfg) ; call new_fuse_cohorts(site, cfg)
            call terminate_cohorts(site, cfg) ; call sort_cohorts(site)
         end if
         call update_overtopping_lai(site)
      end associate
      g_generation(sh) = g_generation(sh) + 1_c_long
   end subroutine meds_apply_rates

   !---------------------------------------------------------------------------------------!
   ! EMPIRICAL-mode TENDENCY computer (the empirical analogue of the driver's carbon           !
   ! compute_slow_derivatives; the former growth_step, split so the ENGINE only applies). Per    !
   ! cohort it advances dbh by the SUPPLIED growth rate [cm/yr] (capped at dbh_critical), forward- !
   ! derives the new geometry + on-allometry pools from the inlined allometry (mirrors            !
   ! set_cohort_size), and backs out each field's time-derivative for the core applier. It also   !
   ! advances the moving-average ring buffer with the SUPPLIED growth rate (the empirical         !
   ! growth_step behaviour -- records the request even when dbh is capped) and records the log-    !
   ! space mortality rate dln_nplant_dt = -mortality.                                             !
   !---------------------------------------------------------------------------------------!
   subroutine compute_empirical_derivatives(site, growth, mortality, dt_yr, n_window, hist_pos, deriv)
      type(site_t),             intent(inout) :: site       ! inout: refreshes the ring buffer
      real(wp),                 intent(in)    :: growth(:), mortality(:)
      real(wp),                 intent(in)    :: dt_yr
      integer(ik),              intent(in)    :: n_window, hist_pos
      type(cohort_deriv_block), intent(inout) :: deriv
      integer(ik) :: i
      real(wp)    :: dbh_new, height_new, ba_new, agb_new, la_new, size_var
      real(wp)    :: lc_new, fc_new, wc_new, nc_new

      associate (cohort => site%cohort)
         do i = 1_ik, cohort%n
            !----- Forward allometry from the supplied dbh growth (inlined; mirrors set_cohort_size). !
            dbh_new    = min(cohort%dbh(i) + growth(i) * dt_yr, cohort%p_dbh_critical(i))
            height_new = min(exp(b1Ht + b2Ht * log(dbh_new)), cohort%p_hgt_max(i))
            ba_new     = pio4 * dbh_new * dbh_new
            size_var   = dbh_new * dbh_new * height_new
            agb_new    = agb_c1 * cohort%p_wood_density(i) ** agb_c2 * size_var ** agb_c2
            la_new     = lai_b1 * size_var ** lai_b2
            lc_new     = la_new / max(cohort%p_sla(i), tiny_num)
            wc_new     = agb_new / max(cohort%p_aboveground_frac(i), tiny_num)
            fc_new     = cohort%p_root_to_leaf_ratio(i) * lc_new
            nc_new     = cohort%p_storage_cushion(i) * lc_new
            !----- Back out every field's time-derivative for the applier. --------------------!
            deriv%d_dbh_dt(i)                  = (dbh_new    - cohort%dbh(i))                  / dt_yr
            deriv%d_height_dt(i)               = (height_new - cohort%height(i))               / dt_yr
            deriv%d_basal_area_dt(i)           = (ba_new     - cohort%basal_area(i))           / dt_yr
            deriv%d_agb_dt(i)                  = (agb_new    - cohort%agb(i))                  / dt_yr
            deriv%d_leaf_area_dt(i)            = (la_new     - cohort%leaf_area(i))            / dt_yr
            deriv%d_leaf_carbon_dt(i)          = (lc_new - cohort%leaf_carbon(i))              / dt_yr
            deriv%d_fineroot_carbon_dt(i)      = (fc_new - cohort%fineroot_carbon(i))          / dt_yr
            deriv%d_wood_carbon_dt(i)          = (wc_new - cohort%wood_carbon(i))              / dt_yr
            deriv%d_nonstructural_carbon_dt(i) = (nc_new - cohort%nonstructural_carbon(i))     / dt_yr
            deriv%dln_nplant_dt(i)             = -mortality(i)
            !----- Advance the ring buffer with the SUPPLIED growth rate (empirical behaviour). --!
            if (cohort%growth_count(i) < n_window) then
               cohort%growth_accum(i) = cohort%growth_accum(i) + growth(i)
               cohort%growth_count(i) = cohort%growth_count(i) + 1_ik
            else
               cohort%growth_accum(i) = cohort%growth_accum(i) + growth(i) - cohort%growth_hist(hist_pos, i)
            end if
            cohort%growth_hist(hist_pos, i) = growth(i)
            cohort%growth_avg(i) = cohort%growth_accum(i) / real(cohort%growth_count(i), wp)
         end do
      end associate
   end subroutine compute_empirical_derivatives

   function meds_site_n_patch(sh) result(np) bind(c, name="meds_site_n_patch")
      integer(c_int), value, intent(in) :: sh
      integer(c_int)                    :: np
      np = int(g_site(sh)%patch%n, c_int)
   end function meds_site_n_patch

   function meds_config_n_pft(ch) result(npft) bind(c, name="meds_config_n_pft")
      integer(c_int), value, intent(in) :: ch
      integer(c_int)                    :: npft
      npft = int(g_cfg(ch)%pft%n, c_int)
   end function meds_config_n_pft

   !=======================================================================================!
   !  Read accessors (scalars). All COPY out; none alias the SoA.                           !
   !=======================================================================================!
   function meds_config_dt_years(ch) result(v) bind(c, name="meds_config_dt_years")
      integer(c_int), value, intent(in) :: ch
      real(c_double)                    :: v
      v = real(g_cfg(ch)%dt_years, c_double)
   end function meds_config_dt_years

   function meds_site_generation(sh) result(g) bind(c, name="meds_site_generation")
      integer(c_int), value, intent(in) :: sh
      integer(c_long)                   :: g
      g = g_generation(sh)
   end function meds_site_generation

   function meds_site_n_cohort(sh) result(n) bind(c, name="meds_site_n_cohort")
      integer(c_int), value, intent(in) :: sh
      integer(c_int)                    :: n
      n = int(count_cohorts(g_site(sh)), c_int)
   end function meds_site_n_cohort

   function meds_site_total_agb(sh) result(v) bind(c, name="meds_site_total_agb")
      integer(c_int), value, intent(in) :: sh
      real(c_double)                    :: v
      v = real(total_agb(g_site(sh)), c_double)
   end function meds_site_total_agb

   function meds_site_total_lai(sh) result(v) bind(c, name="meds_site_total_lai")
      integer(c_int), value, intent(in) :: sh
      real(c_double)                    :: v
      v = real(total_lai(g_site(sh)), c_double)
   end function meds_site_total_lai

   function meds_site_total_nplant(sh) result(v) bind(c, name="meds_site_total_nplant")
      integer(c_int), value, intent(in) :: sh
      real(c_double)                    :: v
      v = real(total_nplant(g_site(sh)), c_double)
   end function meds_site_total_nplant

   function meds_site_total_basal_area(sh) result(v) bind(c, name="meds_site_total_basal_area")
      integer(c_int), value, intent(in) :: sh
      real(c_double)                    :: v
      v = real(total_basal_area(g_site(sh)), c_double)
   end function meds_site_total_basal_area

   !=======================================================================================!
   !  Copy-out per-cohort field getters (caller allocates buf of length n_cohort).          !
   !=======================================================================================!
   subroutine meds_site_get_real(sh, field_id, buf) bind(c, name="meds_site_get_real")
      integer(c_int), value, intent(in)  :: sh, field_id
      real(c_double),        intent(out) :: buf(*)
      integer(ik) :: n
      associate (c => g_site(sh)%cohort)
         n = c%n
         if (n < 1_ik) return
         select case (field_id)
         case (0_c_int) ; buf(1:n) = real(c%dbh(1:n),             c_double)
         case (1_c_int) ; buf(1:n) = real(c%height(1:n),          c_double)
         case (2_c_int) ; buf(1:n) = real(c%nplant(1:n),          c_double)
         case (3_c_int) ; buf(1:n) = real(c%agb(1:n),             c_double)
         case (4_c_int) ; buf(1:n) = real(c%leaf_area(1:n),       c_double)
         case (5_c_int) ; buf(1:n) = real(c%overtopping_lai(1:n), c_double)
         case (6_c_int) ; buf(1:n) = real(c%growth_avg(1:n),      c_double)
         case (7_c_int) ; buf(1:n) = real(c%wood_carbon(1:n),     c_double)
         end select
      end associate
   end subroutine meds_site_get_real

   subroutine meds_site_get_int(sh, field_id, buf) bind(c, name="meds_site_get_int")
      integer(c_int), value, intent(in)  :: sh, field_id
      integer(c_int),        intent(out) :: buf(*)
      integer(ik) :: n
      associate (c => g_site(sh)%cohort)
         n = c%n
         if (n < 1_ik) return
         select case (field_id)
         case (0_c_int) ; buf(1:n) = int(c%pft(1:n),         c_int)
         case (1_c_int) ; buf(1:n) = int(c%owner_patch(1:n), c_int)
         case (2_c_int) ; buf(1:n) = int(c%global_id(1:n),   c_int)
         end select
      end associate
   end subroutine meds_site_get_int

   subroutine meds_site_free(sh) bind(c, name="meds_site_free")
      integer(c_int), value, intent(in) :: sh
      call site_free(g_site(sh))
      site_used(sh) = .false.
   end subroutine meds_site_free

end module meds_demography_capi
