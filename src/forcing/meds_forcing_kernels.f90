!==========================================================================================!
! meds_forcing_kernels -- the PURE/ELEMENTAL meteorological-forcing math (design MEDS_FORCING_ !
! DESIGN.md section 5): per-variable temporal interpolation, the interval-mean-conserving        !
! solar-zenith shortwave disaggregation, the total->4-stream shortwave partition, humidity        !
! from dewpoint / RH, and the precip phase split. Depends ONLY on meds_shared (meds_kinds,          !
! meds_constants, meds_thermo, meds_time) -- REUSES meds_thermo%sat_vapor_pressure (no re-invented   !
! esat) and meds_time%solar_cosz (no re-invented solar geometry). GPU-safe (leaf math over scalars). !
!==========================================================================================!
module meds_forcing_kernels
   use meds_kinds,          only : wp, ik
   use meds_constants,      only : pi, t_3ple, tiny_num
   use meds_thermo,         only : sat_vapor_pressure
   use meds_time,           only : meds_time_t, solar_cosz, day_of_year
   use meds_forcing_config, only : INTERP_LINEAR, INTERP_STEP, INTERP_COSZ,                    &
                                   SWPART_PASSTHROUGH, SWPART_CLEARIDX, SWPART_WEISS_NORMAN
   implicit none
   private

   public :: interpolate_forcing, interpolate_wind_energy
   public :: apparent_solar_seconds, met_solar_cosz, cosz_reconstruct_factor, disaggregate_sw
   public :: partition_shortwave, dewpoint_to_specific_humidity, rh_to_specific_humidity
   public :: precip_phase

   real(wp), parameter :: SOLAR_CONSTANT = 1361.0_wp   !< [W/m2] TOA normal irradiance (mean; eccentricity ~+-3% ignored)
   real(wp), parameter :: COSZ_MIN       = 1.0e-3_wp   !< [-] cosz floor: below this it is night (SW = 0)
   real(wp), parameter :: COSZ_BAR_MIN   = 1.0e-3_wp   !< [-] window-mean-cosz floor for a valid reconstruction
   real(wp), parameter :: FPAR_BEAM      = 0.43_wp     !< [-] PAR fraction of direct-beam SW energy
   real(wp), parameter :: FPAR_DIFFUSE   = 0.52_wp     !< [-] PAR fraction of diffuse SW energy (diffuse is PAR-enriched)
   real(wp), parameter :: PHASE_BAND_K   = 1.0_wp      !< [K] half-width of the rain/snow phase-transition band

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
   !  ERA5-Land ships only total SW, so this is P0-required. SWPART_CLEARIDX: the diffuse           !
   !  fraction from the Erbs (1982) clearness-index correlation (kt = SW / (S0*cosz)); PAR vs NIR   !
   !  by fixed energy fractions (diffuse is PAR-enriched). SWPART_WEISS_NORMAN (band-specific) is    !
   !  P1 and currently falls back to CLEARIDX.                                                        !
   !=======================================================================================!
   pure subroutine partition_shortwave(swdown_total, cosz, scheme,                             &
                                       par_beam, par_diffuse, nir_beam, nir_diffuse)
      real(wp),    intent(in)  :: swdown_total, cosz
      integer(ik), intent(in)  :: scheme
      real(wp),    intent(out) :: par_beam, par_diffuse, nir_beam, nir_diffuse
      real(wp) :: sw, i_toa, kt, fdiff, beam, diff

      sw = max(swdown_total, 0.0_wp)
      if (sw <= tiny_num .or. cosz <= COSZ_MIN) then          ! night / no light
         par_beam = 0.0_wp ; par_diffuse = 0.0_wp ; nir_beam = 0.0_wp ; nir_diffuse = 0.0_wp
         return
      end if

      select case (scheme)
      case default   ! SWPART_CLEARIDX (also the P1 WEISS_NORMAN fallback)
         i_toa = SOLAR_CONSTANT * cosz                        ! TOA irradiance on a horizontal surface
         kt    = min(1.0_wp, max(0.0_wp, sw / max(i_toa, tiny_num)))
         fdiff = erbs_diffuse_fraction(kt)
      end select

      diff = fdiff * sw
      beam = sw - diff
      par_beam    = FPAR_BEAM    * beam
      nir_beam    = beam - par_beam
      par_diffuse = FPAR_DIFFUSE * diff
      nir_diffuse = diff - par_diffuse
   end subroutine partition_shortwave

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
   !  HUMIDITY. Both forms reuse meds_thermo%sat_vapor_pressure (Bolton 1980) -- IDENTICAL to     !
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

end module meds_forcing_kernels
