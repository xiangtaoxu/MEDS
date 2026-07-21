!==========================================================================================!
! meds_forcing_kernels -- the PURE/ELEMENTAL meteorological-forcing math (design MEDS_FORCING_ !
! DESIGN.md section 5): per-variable temporal interpolation, the interval-mean-conserving        !
! solar-zenith shortwave disaggregation, the total->4-stream shortwave partition, humidity        !
! from dewpoint / RH, and the precip phase split. Depends ONLY on meds_shared (meds_kinds,          !
! meds_constants, meds_therm_lib, meds_time) -- REUSES meds_therm_lib%sat_vapor_pressure (no re-invented   !
! esat) and meds_time%solar_cosz (no re-invented solar geometry). GPU-safe (leaf math over scalars). !
!==========================================================================================!
module meds_forcing_kernels
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : pi, t_3ple, tiny_num, p_std, grav, r_dry
   use meds_therm_lib,         only : sat_vapor_pressure
   use meds_time,           only : meds_time_t, solar_cosz, day_of_year
   use meds_forcing_config, only : INTERP_LINEAR, INTERP_STEP, INTERP_COSZ,                    &
                                   SWPART_PASSTHROUGH, SWPART_CLEARIDX, SWPART_WEISS_NORMAN
   implicit none
   private

   public :: interpolate_forcing, interpolate_wind_energy
   public :: apparent_solar_seconds, met_solar_cosz, cosz_reconstruct_factor, disaggregate_sw
   public :: partition_shortwave, dewpoint_to_specific_humidity, rh_to_specific_humidity
   public :: precip_phase
   public :: great_circle_distance, nearest_grid_index
   public :: wind_log_profile, lapse_air_temperature, lapse_pressure

   real(wp), parameter :: SOLAR_CONSTANT = 1361.0_wp   !< [W/m2] TOA normal irradiance (mean; eccentricity ~+-3% ignored)
   real(wp), parameter :: COSZ_MIN       = 1.0e-3_wp   !< [-] cosz floor: below this it is night (SW = 0)
   real(wp), parameter :: COSZ_BAR_MIN   = 1.0e-3_wp   !< [-] window-mean-cosz floor for a valid reconstruction
   real(wp), parameter :: FPAR_BEAM      = 0.43_wp     !< [-] PAR fraction of direct-beam SW energy
   real(wp), parameter :: FPAR_DIFFUSE   = 0.52_wp     !< [-] PAR fraction of diffuse SW energy (diffuse is PAR-enriched)
   real(wp), parameter :: PHASE_BAND_K   = 1.0_wp      !< [K] half-width of the rain/snow phase-transition band
   real(wp), parameter :: EARTH_RADIUS_M = 6.371e6_wp  !< [m] mean Earth radius (great-circle grid match; argmin only)

   !----- Weiss & Norman (1985) band-specific SW partition constants (port of ED2 -------------!
   !      short_bdown_weissnorman, radiate_utils.f90). Empirical/physical -> module parameters.  !
   real(wp), parameter :: WN_PAR_BEAM_EXPEXT = -0.185_wp   !< visible-beam atmospheric extinction (WN85 eq.1)
   real(wp), parameter :: WN_NIR_BEAM_EXPEXT = -0.060_wp   !< NIR-beam atmospheric extinction (WN85 eq.4)
   real(wp), parameter :: WN_PAR2DIFF_SUN    = 0.400_wp    !< visible beam->diffuse potential fraction (WN85 eq.3)
   real(wp), parameter :: WN_NIR2DIFF_SUN    = 0.600_wp    !< NIR beam->diffuse potential fraction (WN85 eq.5)
   real(wp), parameter :: WN_W10_A = -1.1950_wp, WN_W10_B = 0.4459_wp, WN_W10_C = -0.0345_wp  !< NIR water absorption (WN85 eq.6)
   real(wp), parameter :: WN_PAR_ACT_A = 0.90_wp, WN_PAR_ACT_B = 0.70_wp   !< visible actual-beam fraction (WN85 eq.11)
   real(wp), parameter :: WN_NIR_ACT_A = 0.88_wp, WN_NIR_ACT_B = 0.68_wp   !< NIR actual-beam fraction (WN85 eq.12)
   real(wp), parameter :: WN_FVIS_BEAM = 0.43_wp, WN_FNIR_BEAM = 0.57_wp   !< TOA visible/NIR beam energy fractions
   real(wp), parameter :: WN_FVIS_DIFF = 0.52_wp, WN_FNIR_DIFF = 0.48_wp   !< dawn all-diffuse visible/NIR fractions
   real(wp), parameter :: WN_TWOTHIRDS = 2.0_wp / 3.0_wp
   real(wp), parameter :: WN_COSZ_MIN  = 0.017452406437_wp  !< cos(89 deg): below this the sun is grazing -> all diffuse

contains

   !=======================================================================================!
   !  TEMPORAL INTERPOLATION between the two bracketing records (weight w_next in [0,1]).      !
   !=======================================================================================!
   elemental function interpolate_forcing(policy, prev, next, w_next) result(val)
      integer(ik), intent(in) :: policy        !< INTERP_LINEAR | INTERP_STEP
      real(wp),    intent(in) :: prev, next, w_next
      real(wp) :: val
      select case (policy)
      case (INTERP_STEP)                        ! step-constant: hold the previous value across the interval
         val = prev
      case default                              ! INTERP_LINEAR
         val = (1.0_wp - w_next) * prev + w_next * next
      end select
   end function interpolate_forcing

   !----- Wind interpolates as SQUARED quantities (kinetic-energy-conserving, ED2 vels). ----!
   elemental function interpolate_wind_energy(prev, next, w_next, u_min) result(u)
      real(wp), intent(in) :: prev, next, w_next, u_min
      real(wp) :: u
      u = sqrt(max((1.0_wp - w_next) * prev * prev + w_next * next * next, u_min * u_min))
   end function interpolate_wind_energy

   !=======================================================================================!
   !  SOLAR TIME + ZENITH. ERA5-Land is UTC, so local apparent solar seconds require the         !
   !  longitude + equation-of-time correction (design §5.1). The single transform:               !
   !     t_solar = t_utc + (longitude - 15*utc_offset)*240 s/deg + eot(doy)                        !
   !  reduces to t_utc + longitude*240 + eot for a UTC file (utc_offset = 0). A local-clock         !
   !  file (utc_offset /= 0) shares the same one line. apply_solar_longitude = .false. passes the    !
   !  clock seconds straight through (a file already in local apparent solar time).                   !
   !=======================================================================================!
   pure function equation_of_time(t) result(eot_sec)
      type(meds_time_t), intent(in) :: t
      real(wp) :: eot_sec, b
      b = 2.0_wp * pi * real(day_of_year(t) - 1_ik, wp) / 365.0_wp        ! Spencer 1971
      eot_sec = 229.18_wp * ( 0.000075_wp + 0.001868_wp*cos(b) - 0.032077_wp*sin(b)             &
                            - 0.014615_wp*cos(2.0_wp*b) - 0.040849_wp*sin(2.0_wp*b) ) * 60.0_wp
   end function equation_of_time

   pure function apparent_solar_seconds(t, sec_clock, longitude_deg, utc_offset_h, apply_lon) result(s)
      type(meds_time_t), intent(in) :: t              !< record/model time (for declination + EoT via doy)
      real(wp),          intent(in) :: sec_clock      !< [s] seconds into day on the FILE clock
      real(wp),          intent(in) :: longitude_deg, utc_offset_h
      logical,           intent(in) :: apply_lon
      real(wp) :: s
      if (apply_lon) then
         s = sec_clock + (longitude_deg - 15.0_wp * utc_offset_h) * 240.0_wp + equation_of_time(t)
      else
         s = sec_clock                                 ! file clock IS local apparent solar time
      end if
   end function apparent_solar_seconds

   !----- cos(solar zenith) at a file-clock sub-daily cursor (converts to solar time first). --!
   pure function met_solar_cosz(t, sec_clock, latitude_deg, longitude_deg, utc_offset_h, apply_lon) result(cosz)
      type(meds_time_t), intent(in) :: t
      real(wp),          intent(in) :: sec_clock, latitude_deg, longitude_deg, utc_offset_h
      logical,           intent(in) :: apply_lon
      real(wp) :: cosz, sec_solar
      sec_solar = apparent_solar_seconds(t, sec_clock, longitude_deg, utc_offset_h, apply_lon)
      cosz = solar_cosz(t, sec_solar, latitude_deg)
   end function met_solar_cosz

   !=======================================================================================!
   !  SHORTWAVE DISAGGREGATION (design §5.1, the load-bearing method). Reconstruct an           !
   !  instantaneous flux F(now) = F_avg*cosz(now)/<cosz>_win from an interval-mean F_avg so its   !
   !  time-average over the SAME window returns F_avg. The factor is 1/<cosz>_win -- the           !
   !  RECIPROCAL OF THE MEAN COSINE (night sub-samples clamped to 0), NOT <sec z> (which is         !
   !  biased high at low sun and does NOT conserve the interval mean). Fully-night window ->         !
   !  <cosz>_win = 0 and F_avg = 0, so factor is returned 0 and all SW routes to 0.                   !
   !=======================================================================================!
   pure function cosz_reconstruct_factor(t, win_start_sec, dt_sub, dt_win,                      &
                                         latitude_deg, longitude_deg, utc_offset_h, apply_lon) result(factor)
      type(meds_time_t), intent(in) :: t
      real(wp),          intent(in) :: win_start_sec, dt_sub, dt_win
      real(wp),          intent(in) :: latitude_deg, longitude_deg, utc_offset_h
      logical,           intent(in) :: apply_lon
      real(wp)    :: factor, cosz_sum, cosz_bar, sec, ci
      integer(ik) :: nsub, i
      nsub = max(1_ik, nint(dt_win / max(dt_sub, tiny_num), ik))
      cosz_sum = 0.0_wp
      do i = 1_ik, nsub                                   ! midpoint rule over the FULL window
         sec = win_start_sec + (real(i, wp) - 0.5_wp) * (dt_win / real(nsub, wp))
         ci  = met_solar_cosz(t, sec, latitude_deg, longitude_deg, utc_offset_h, apply_lon)
         cosz_sum = cosz_sum + max(ci, 0.0_wp)            ! night sub-samples contribute exactly 0
      end do
      cosz_bar = cosz_sum / real(nsub, wp)
      if (cosz_bar > COSZ_BAR_MIN) then
         factor = 1.0_wp / cosz_bar
      else
         factor = 0.0_wp                                  ! fully-night window; F_avg is 0 too
      end if
   end function cosz_reconstruct_factor

   !----- Instantaneous SW stream from its interval mean: F_now = F_avg*cosz(now)*factor. ------!
   !      Night (cosz(now) <= COSZ_MIN) -> 0 automatically; factor already 0 for a night window. !
   elemental function disaggregate_sw(f_avg, cosz_now, factor) result(f_now)
      real(wp), intent(in) :: f_avg, cosz_now, factor
      real(wp) :: f_now
      if (cosz_now > COSZ_MIN) then
         f_now = f_avg * cosz_now * factor
      else
         f_now = 0.0_wp
      end if
   end function disaggregate_sw

   !=======================================================================================!
   !  SHORTWAVE PARTITION: total SWdown -> the four (beam/diffuse) x (PAR/NIR) streams (§5.6).   !
   !  ERA5-Land ships only total SW, so this is P0-required. SWPART_CLEARIDX (default): the         !
   !  diffuse fraction from the Erbs (1982) clearness-index correlation (kt = SW/(S0*cosz)); PAR    !
   !  vs NIR by fixed energy fractions (diffuse is PAR-enriched). SWPART_WEISS_NORMAN: the           !
   !  band-specific, pressure-dependent Weiss & Norman (1985) split (ED2 short_bdown_weissnorman).   !
   !  BOTH conserve total energy exactly (the four streams sum to sw).                                !
   !=======================================================================================!
   pure subroutine partition_shortwave(swdown_total, cosz, psurf_pa, scheme,                   &
                                       par_beam, par_diffuse, nir_beam, nir_diffuse)
      real(wp),    intent(in)  :: swdown_total, cosz, psurf_pa
      integer(ik), intent(in)  :: scheme
      real(wp),    intent(out) :: par_beam, par_diffuse, nir_beam, nir_diffuse
      real(wp) :: sw, i_toa, kt, fdiff, beam, diff

      sw = max(swdown_total, 0.0_wp)
      if (sw <= tiny_num .or. cosz <= COSZ_MIN) then          ! night / no light
         par_beam = 0.0_wp ; par_diffuse = 0.0_wp ; nir_beam = 0.0_wp ; nir_diffuse = 0.0_wp
         return
      end if

      select case (scheme)
      case (SWPART_WEISS_NORMAN)
         call weiss_norman_partition(sw, cosz, psurf_pa, par_beam, par_diffuse, nir_beam, nir_diffuse)
      case default   ! SWPART_CLEARIDX (also SWPART_SIB / anything-else fallback)
         i_toa = SOLAR_CONSTANT * cosz                        ! TOA irradiance on a horizontal surface
         kt    = min(1.0_wp, max(0.0_wp, sw / max(i_toa, tiny_num)))
         fdiff = erbs_diffuse_fraction(kt)
         diff = fdiff * sw
         beam = sw - diff
         par_beam    = FPAR_BEAM    * beam
         nir_beam    = beam - par_beam
         par_diffuse = FPAR_DIFFUSE * diff
         nir_diffuse = diff - par_diffuse
      end select
   end subroutine partition_shortwave

   !=======================================================================================!
   !  WEISS & NORMAN (1985) band-specific SW partition (port of ED2 short_bdown_weissnorman,     !
   !  radiate_utils.f90). Potential band beam/diffuse fluxes from air-mass optical depth on a      !
   !  clear sky; the ratio of observed-to-potential total sets the actual beam fraction per band    !
   !  (cloudier -> more diffuse). Conserves energy EXACTLY: par_full + nir_full = ratio*(pot) = sw, !
   !  and beam + diffuse = full within each band. Grazing sun (cosz <= WN_COSZ_MIN) -> all diffuse.  !
   !=======================================================================================!
   pure subroutine weiss_norman_partition(sw, cosz, psurf_pa, par_beam, par_diffuse, nir_beam, nir_diffuse)
      real(wp), intent(in)  :: sw, cosz, psurf_pa
      real(wp), intent(out) :: par_beam, par_diffuse, nir_beam, nir_diffuse
      real(wp) :: secz, log10secz, prat, w10
      real(wp) :: par_beam_top, nir_beam_top, par_beam_pot, par_diff_pot, par_full_pot
      real(wp) :: nir_beam_pot, nir_diff_pot, nir_full_pot, ratio, par_full, nir_full
      real(wp) :: aux_par, aux_nir, fvis_beam, fnir_beam

      if (cosz <= WN_COSZ_MIN) then                           ! grazing sun: secant unstable -> all diffuse
         par_beam = 0.0_wp ; nir_beam = 0.0_wp
         par_diffuse = WN_FVIS_DIFF * sw ; nir_diffuse = WN_FNIR_DIFF * sw   ! 0.52+0.48 = 1 (conserves sw)
         return
      end if

      secz      = 1.0_wp / cosz
      log10secz = log10(secz)
      prat      = psurf_pa / p_std
      par_beam_top = WN_FVIS_BEAM * SOLAR_CONSTANT           ! TOA visible/NIR beam energy
      nir_beam_top = WN_FNIR_BEAM * SOLAR_CONSTANT

      !----- potential (clear-sky) visible beam + diffuse (WN85 eq.1, eq.3, eq.9). ------------!
      par_beam_pot = par_beam_top * exp(WN_PAR_BEAM_EXPEXT * prat * secz) * cosz
      par_diff_pot = WN_PAR2DIFF_SUN * (par_beam_top - par_beam_pot) * cosz
      par_full_pot = par_beam_pot + par_diff_pot

      !----- potential NIR: water-vapour absorption w10, then beam + diffuse (WN85 eq.6,4,5,10). !
      w10          = SOLAR_CONSTANT * 10.0_wp ** (WN_W10_A + log10secz * (WN_W10_B + WN_W10_C * log10secz))
      nir_beam_pot = (nir_beam_top * exp(WN_NIR_BEAM_EXPEXT * prat * secz) - w10) * cosz
      nir_diff_pot = WN_NIR2DIFF_SUN * (nir_beam_top - nir_beam_pot - w10) * cosz
      nir_full_pot = nir_beam_pot + nir_diff_pot

      !----- scale the potential to the OBSERVED total; beam fraction shrinks under cloud (eq.7/8). !
      !      NOTE (ED2-faithful edge): in a narrow band just above the grazing cutoff (cosz ~           !
      !      0.0175-0.029) the NIR water absorption w10 exceeds the attenuated NIR beam, so            !
      !      nir_full_pot < 0 and the scaled nir_diffuse can go slightly NEGATIVE (with PAR > sw)       !
      !      while the four streams still sum to sw EXACTLY. ED2's short_bdown_weissnorman does the     !
      !      same; the window is extreme twilight where SW is tiny, and WN is opt-in (Erbs is default). !
      ratio    = sw / max(par_full_pot + nir_full_pot, tiny_num)
      par_full = ratio * par_full_pot
      nir_full = ratio * nir_full_pot

      !----- actual beam fractions per band (WN85 eq.11/12; clipped to [0,1]). ----------------!
      aux_par   = min(WN_PAR_ACT_A, max(0.0_wp, ratio))
      aux_nir   = min(WN_NIR_ACT_A, max(0.0_wp, ratio))
      fvis_beam = min(1.0_wp, max(0.0_wp,                                                       &
                  par_beam_pot * (1.0_wp - ((WN_PAR_ACT_A - aux_par) / WN_PAR_ACT_B) ** WN_TWOTHIRDS) &
                  / max(par_full_pot, tiny_num)))
      fnir_beam = min(1.0_wp, max(0.0_wp,                                                       &
                  nir_beam_pot * (1.0_wp - ((WN_NIR_ACT_A - aux_nir) / WN_NIR_ACT_B) ** WN_TWOTHIRDS) &
                  / max(nir_full_pot, tiny_num)))

      par_beam    = fvis_beam * par_full ; par_diffuse = (1.0_wp - fvis_beam) * par_full
      nir_beam    = fnir_beam * nir_full ; nir_diffuse = (1.0_wp - fnir_beam) * nir_full
   end subroutine weiss_norman_partition

   !----- Erbs et al. (1982) diffuse fraction of hourly global radiation vs clearness index. ---!
   pure function erbs_diffuse_fraction(kt) result(fd)
      real(wp), intent(in) :: kt
      real(wp) :: fd
      if (kt <= 0.22_wp) then
         fd = 1.0_wp - 0.09_wp * kt
      else if (kt <= 0.80_wp) then
         fd = 0.9511_wp - 0.1604_wp*kt + 4.388_wp*kt**2 - 16.638_wp*kt**3 + 12.336_wp*kt**4
      else
         fd = 0.165_wp
      end if
      fd = min(1.0_wp, max(0.0_wp, fd))
   end function erbs_diffuse_fraction

   !=======================================================================================!
   !  HUMIDITY. Both forms reuse meds_therm_lib%sat_vapor_pressure (Bolton 1980) -- IDENTICAL to     !
   !  the formatter script, so file-built Qair reconciles with any reader-side humidity math.     !
   !=======================================================================================!
   !----- Specific humidity from dewpoint: actual e = e_sat AT the dewpoint. -----------------!
   elemental function dewpoint_to_specific_humidity(td_k, p_pa) result(q)
      real(wp), intent(in) :: td_k, p_pa
      real(wp) :: q, e
      e = sat_vapor_pressure(td_k)
      q = 0.622_wp * e / max(p_pa - 0.378_wp * e, tiny_num)
   end function dewpoint_to_specific_humidity

   !----- Specific humidity from relative humidity (fraction) + temperature + pressure. ------!
   elemental function rh_to_specific_humidity(rh, t_k, p_pa) result(q)
      real(wp), intent(in) :: rh, t_k, p_pa
      real(wp) :: q, e
      e = min(1.0_wp, max(0.0_wp, rh)) * sat_vapor_pressure(t_k)
      q = 0.622_wp * e / max(p_pa - 0.378_wp * e, tiny_num)
   end function rh_to_specific_humidity

   !=======================================================================================!
   !  PRECIP PHASE (ED2 Jin 1999 style, simplified): a linear rain/snow ramp across a band      !
   !  centered on the triple point. Mass-conserving (rainf + snowf = total).                     !
   !=======================================================================================!
   elemental subroutine precip_phase(precip_total, tair_k, rainf, snowf)
      real(wp), intent(in)  :: precip_total, tair_k
      real(wp), intent(out) :: rainf, snowf
      real(wp) :: frac_liq
      frac_liq = min(1.0_wp, max(0.0_wp, (tair_k - (t_3ple - PHASE_BAND_K)) / (2.0_wp * PHASE_BAND_K)))
      rainf = frac_liq * precip_total
      snowf = precip_total - rainf
   end subroutine precip_phase

   !=======================================================================================!
   !  GRID MATCH (multi-polygon P2 subset): great-circle distance + nearest-cell argmin.        !
   !  great_circle_distance ports ED2 dist_gc (great_circle.f90); nearest_grid_index ports the    !
   !  match_poly_grid argmin loop. Only the argmin matters for selection, so EARTH_RADIUS_M is     !
   !  immaterial; strict '<' keeps the lowest index on ties.                                        !
   !=======================================================================================!
   pure function great_circle_distance(lon1_deg, lat1_deg, lon2_deg, lat2_deg) result(dist_m)
      real(wp), intent(in) :: lon1_deg, lat1_deg, lon2_deg, lat2_deg
      real(wp) :: dist_m, d2r, lat1, lat2, dlon, x, y
      d2r  = pi / 180.0_wp
      lat1 = lat1_deg * d2r ; lat2 = lat2_deg * d2r
      dlon = (lon2_deg - lon1_deg) * d2r
      x = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(dlon)
      y = sqrt( (cos(lat2) * sin(dlon)) ** 2                                                    &
              + (cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon)) ** 2 )
      dist_m = EARTH_RADIUS_M * atan2(y, x)
   end function great_circle_distance

   pure function nearest_grid_index(site_lon, site_lat, lon, lat) result(idx)
      real(wp), intent(in) :: site_lon, site_lat, lon(:), lat(:)
      integer(ik) :: idx, i
      real(wp) :: dmin, d
      idx = 1_ik ; dmin = huge(1.0_wp)
      do i = 1_ik, int(size(lon), ik)
         d = great_circle_distance(site_lon, site_lat, lon(i), lat(i))
         if (d < dmin) then ; dmin = d ; idx = i ; end if
      end do
   end function nearest_grid_index

   !=======================================================================================!
   !  WIND-HEIGHT + ELEVATION-LAPSE corrections (design §5.2 / §10-Q2). Applied at ingest to      !
   !  the raw file record; both are OFF by default (forcing_config gates). Pure kernels.            !
   !=======================================================================================!
   !----- Neutral-log wind from the measurement height to the model reference height. The factor !
   !      is independent of u (commutes with the energy-form interpolation); degenerate z0 -> no-op. !
   elemental function wind_log_profile(u_meas, z_meas, z_ref, z0) result(u_ref)
      real(wp), intent(in) :: u_meas, z_meas, z_ref, z0
      real(wp) :: u_ref
      if (z0 <= 0.0_wp .or. z_meas <= z0 .or. z_ref <= z0) then
         u_ref = u_meas                                       ! degenerate: leave the wind untouched
      else
         u_ref = u_meas * log(z_ref / z0) / log(z_meas / z0)
      end if
   end function wind_log_profile

   !----- Linear environmental lapse of air temperature (ED2 calc_met_lapse). dz = site - grid,  !
   !      gamma > 0 cools with height, so a higher site is colder.                                  !
   elemental function lapse_air_temperature(tair_grid, dz, gamma) result(tair_site)
      real(wp), intent(in) :: tair_grid, dz, gamma
      real(wp) :: tair_site
      tair_site = tair_grid - gamma * dz
   end function lapse_air_temperature

   !----- Hydrostatic hypsometric pressure CONSISTENT with the same linear T(z): dP/dz=-Pg/(Rd T),  !
   !      T(z)=T_grid-gamma*z -> P_site = P_grid*(T_site/T_grid)^(g/(Rd*gamma)); isothermal limit    !
   !      (gamma -> 0) is the barometric exp(-g*dz/(Rd*T)). Keeps P and T ideal-gas-consistent.      !
   elemental function lapse_pressure(psurf_grid, tair_grid, dz, gamma) result(psurf_site)
      real(wp), intent(in) :: psurf_grid, tair_grid, dz, gamma
      real(wp) :: psurf_site, tair_site
      real(wp), parameter :: LAPSE_GAMMA_MIN = 1.0e-6_wp
      if (abs(gamma) > LAPSE_GAMMA_MIN) then
         tair_site  = tair_grid - gamma * dz
         psurf_site = psurf_grid * (tair_site / tair_grid) ** (grav / (r_dry * gamma))
      else
         psurf_site = psurf_grid * exp(-grav * dz / (r_dry * tair_grid))
      end if
   end function lapse_pressure

end module meds_forcing_kernels
