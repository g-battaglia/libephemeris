//! Rise set implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Resolve the event type from the `rsmi` low bits with the release
/// compatibility precedence `SET > RISE > ITRANSIT > MTRANSIT`. A
/// modifier-only mask or a missing event type defaults to
/// [`RiseSetEvent::Rise`].
pub fn resolve_event_kind(rsmi: i32) -> RiseSetEvent {
    if rsmi & CALC_SET != 0 {
        RiseSetEvent::Set
    } else if rsmi & CALC_RISE != 0 {
        RiseSetEvent::Rise
    } else if rsmi & CALC_ITRANSIT != 0 {
        RiseSetEvent::LowerTransit
    } else if rsmi & CALC_MTRANSIT != 0 {
        RiseSetEvent::UpperTransit
    } else {
        RiseSetEvent::Rise
    }
}

/// Validate a geographic observer triple: longitude/latitude in
/// degrees, altitude in meters. Non-finite components are
/// [`Error::InvalidObserver`]; the release `ValueError` for a short or
/// non-numeric *sequence* lives in the binding adapter, which owns
/// sequence decoding. Latitude/range semantics belong to the `azalt`
/// serving-layer binding.
pub fn validate_geopos(lon_deg: f64, lat_deg: f64, alt_m: f64) -> Result<()> {
    if !lon_deg.is_finite() {
        return Err(Error::InvalidObserver {
            field: "lon_deg",
            value: lon_deg,
        });
    }
    if !lat_deg.is_finite() {
        return Err(Error::InvalidObserver {
            field: "lat_deg",
            value: lat_deg,
        });
    }
    if !alt_m.is_finite() {
        return Err(Error::InvalidObserver {
            field: "alt_m",
            value: alt_m,
        });
    }
    Ok(())
}

/// Estimate atmospheric pressure from observer altitude with the release
/// barometric expression. An input pressure of 0 means "estimate", not "no
/// pressure" (the binding checks that); the temperature is used
/// (0 means 0 C, never an implicit standard atmosphere). The base goes
/// negative above ~44.3 km, where a fractional power would not be real; it
/// is clamped to 0 so extreme altitudes give pressure ~0 (refraction ~0)
/// instead of an error.
pub fn pressure_from_altitude_m(alt_m: f64) -> f64 {
    let base = 1.0 - 0.0065 * alt_m / 288.0;
    SEA_LEVEL_PRESSURE_HPA * base.max(0.0).powf(5.255)
}

/// Return the twilight override when `rsmi` carries a twilight bit and the
/// body is the Sun or a fixed star, else `None`.
pub fn twilight_override(rsmi: i32, is_sun_or_star: bool) -> Option<TwilightOverride> {
    if !is_sun_or_star {
        return None;
    }
    let horizon_alt_deg = if rsmi & BIT_CIVIL_TWILIGHT != 0 {
        CIVIL_TWILIGHT_ALT_DEG
    } else if rsmi & BIT_NAUTIC_TWILIGHT != 0 {
        NAUTIC_TWILIGHT_ALT_DEG
    } else if rsmi & BIT_ASTRO_TWILIGHT != 0 {
        ASTRO_TWILIGHT_ALT_DEG
    } else {
        return None;
    };
    Some(TwilightOverride {
        horizon_alt_deg,
        force_no_refraction: true,
        force_disc_center: true,
    })
}

/// Select the crossing disc point: fixed stars always use
/// the center; else `BIT_DISC_CENTER` outranks `BIT_DISC_BOTTOM`, which
/// outranks the upper-limb default.
pub fn select_disc_point(rsmi_eff: i32, is_fixed_star: bool) -> DiscPoint {
    if rsmi_eff & BIT_DISC_CENTER != 0 || is_fixed_star {
        DiscPoint::Center
    } else if rsmi_eff & BIT_DISC_BOTTOM != 0 {
        DiscPoint::LowerLimb
    } else {
        DiscPoint::UpperLimb
    }
}

/// Semidiameter altitude term for the selected disc point, degrees
/// `+/- asin(min(1, R / dist))`
/// from the actual distance, so the term tracks the Earth-Sun distance
/// through the year instead of a fixed 16 arcminutes. `Center` and
/// non-positive radii contribute 0. A non-finite or non-positive distance
/// is [`Error::NonFiniteSearchValue`]: the release divides by the
/// ephemeris distance, always positive in practice, and never guards this
/// path. The `BIT_FIXED_DISC_SIZE` pinned distances
/// ([`FIXED_DISC_SUN_DIST_AU`]/[`FIXED_DISC_MOON_DIST_AU`]) are applied by
/// the caller before this function.
pub fn disc_altitude_term_deg(point: DiscPoint, body_radius_au: f64, dist_au: f64) -> Result<f64> {
    let sign = match point {
        DiscPoint::UpperLimb => 1.0,
        DiscPoint::Center => return Ok(0.0),
        DiscPoint::LowerLimb => -1.0,
    };
    if body_radius_au <= 0.0 || !body_radius_au.is_finite() {
        return Ok(0.0);
    }
    if !dist_au.is_finite() || dist_au <= 0.0 {
        return Err(Error::NonFiniteSearchValue { x: dist_au });
    }
    Ok(sign * (body_radius_au / dist_au).min(1.0).asin().to_degrees())
}

/// Resolve the event horizon altitude, degrees. The
/// sentinel [`HORIZON_DIP_SENTINEL_DEG`] asks for the dip of the sea
/// horizon at the observer's altitude — computed by the `calc_dip`
/// horizon-dip kernel on the binding side and passed as `dip_deg` — plus
/// [`HORIZON_DIP_EPSILON_DEG`]. Any other value is used:
/// negative values are legitimate (a depressed horizon for an elevated
/// observer).
pub fn resolve_horizon_altitude_deg(horhgt_deg: f64, dip_deg: f64) -> f64 {
    if horhgt_deg == HORIZON_DIP_SENTINEL_DEG {
        dip_deg + HORIZON_DIP_EPSILON_DEG
    } else {
        horhgt_deg
    }
}

/// Build the 10-element `tret` compatibility tuple: only index 0 is
/// meaningful (the Julian Day of the event); the remaining slots are 0.0.
/// A circumpolar outcome carries 0.0 in slot 0.
pub fn make_tret(jd_event: f64) -> [f64; 10] {
    [jd_event, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
}

/// True-to-apparent altitude for rise/set events. The
/// refraction law is a caller input — `refraction_deg(app)` returns
/// `R(app)` in degrees — and
/// `dip_deg` is the observer's horizon dip. `x_hint` labels non-finite evaluations.
///
/// The apparent altitude solves `app = true + R(app)` by Newton iteration
/// with a numerical derivative (`eps = 1e-4`, at most 10 steps, stopping
/// below `1e-9`). When the refracted point would still sit below the dip,
/// the altitude is left unrefracted — the release convention for rise/set.
/// Non-finite refraction values are [`Error::NonFiniteSearchValue`].
pub fn true_to_apparent_altitude(
    true_alt_deg: f64,
    dip_deg: f64,
    refraction_deg: &mut impl FnMut(f64) -> std::result::Result<f64, Error>,
    x_hint: f64,
) -> Result<f64> {
    if !true_alt_deg.is_finite() || !dip_deg.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: x_hint });
    }
    let checked_r = &mut |app: f64| -> Result<f64> {
        let r = refraction_deg(app)?;
        if !r.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: app });
        }
        Ok(r)
    };
    let mut app = true_alt_deg + checked_r(true_alt_deg)?;
    if !app.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: x_hint });
    }
    for _ in 0..10 {
        let f = app - checked_r(app)? - true_alt_deg;
        let eps = 1e-4;
        let f2 = (app + eps) - checked_r(app + eps)? - true_alt_deg;
        let dfd = (f2 - f) / eps;
        let step = if dfd != 0.0 { f / dfd } else { f };
        app -= step;
        if !app.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: x_hint });
        }
        if step.abs() < 1e-9 {
            break;
        }
    }
    let refr = app - true_alt_deg;
    if true_alt_deg + refr < dip_deg {
        Ok(true_alt_deg)
    } else {
        Ok(true_alt_deg + refr)
    }
}

/// Compose the event height from its parts:
/// the disc-raised true altitude, optionally replaced by its refracted
/// apparent value, minus the horizon altitude. Positive means the
/// configured disc point stands above the configured horizon.
/// `apparent_alt_deg` is `Some` exactly when refraction applies.
pub fn event_height_deg(
    alt_point_true_deg: f64,
    apparent_alt_deg: Option<f64>,
    horizon_alt_deg: f64,
) -> f64 {
    apparent_alt_deg.unwrap_or(alt_point_true_deg) - horizon_alt_deg
}

/// Hour angle from apparent RA, GAST and observer longitude, degrees in
/// `(-180, 180]` ` first, so values at or below -180 wrap up — mirrored here
/// via `mod360_deg`, unlike the symmetric [`wrapped_residual_deg`]).
pub(super) fn transit_hour_angle_deg(ra_hours: f64, gast_hours: f64, lon_deg_east: f64) -> f64 {
    let ha = mod360_deg((gast_hours + lon_deg_east / 15.0) * 15.0 - ra_hours * 15.0);
    if ha > 180.0 { ha - 360.0 } else { ha }
}

/// Meridian transit time from backend-agnostic providers: `get_ra_dec(jd)` returns apparent `(RA hours,
/// Dec deg)`; `gast_hours(jd)` returns Greenwich *apparent* sidereal time
/// in hours — the RA is apparent (equinox of date with nutation), so GAST
/// is the consistent pairing.
///
/// `lower_transit` selects the anti-culmination (`CALC_ITRANSIT`
/// 180-degree shift with the same single-branch wrap). The seed advances
/// by mean-sidereal-day fraction from the current hour angle; each step
/// subtracts `ha / (360 / SIDEREAL_DAY)` days. Non-finite provider values
/// are [`Error::NonFiniteSearchValue`].
pub fn calculate_transit(
    jd_start: f64,
    lon_deg_east: f64,
    lower_transit: bool,
    get_ra_dec: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    gast_hours: &mut impl FnMut(f64) -> std::result::Result<f64, Error>,
) -> Result<TransitSolution> {
    if !jd_start.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: jd_start });
    }
    if !lon_deg_east.is_finite() {
        return Err(Error::InvalidObserver {
            field: "lon_deg",
            value: lon_deg_east,
        });
    }
    let hour_angle = &mut |jd: f64| -> Result<f64> {
        let (ra_hours, dec_deg) = get_ra_dec(jd)?;
        if !ra_hours.is_finite() || !dec_deg.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        let gast = gast_hours(jd)?;
        if !gast.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        let mut ha = transit_hour_angle_deg(ra_hours, gast, lon_deg_east);
        if lower_transit {
            ha = mod360_deg(ha - 180.0);
            if ha > 180.0 {
                ha -= 360.0;
            }
        }
        Ok(ha)
    };

    let ha = hour_angle(jd_start)?;
    let mut jd_guess = if ha > 0.0 {
        jd_start + (360.0 - ha) / 360.0 * SIDEREAL_DAY
    } else {
        jd_start + (-ha) / 360.0 * SIDEREAL_DAY
    };
    for _ in 0..TRANSIT_MAX_ITER {
        let ha = hour_angle(jd_guess)?;
        if ha.abs() < TRANSIT_TOL_HA_DEG {
            return Ok(TransitSolution {
                jd_ut: jd_guess,
                converged: true,
            });
        }
        jd_guess -= ha / (360.0 / SIDEREAL_DAY);
    }
    Ok(TransitSolution {
        jd_ut: jd_guess,
        converged: false,
    })
}

/// Find the next rise or set: first zero crossing of `event_height`
///
/// `event_height(jd)` is positive when the configured disc point stands
/// above the configured horizon (disc, refraction and horizon already
/// applied by the caller); `get_ra_dec(jd)` returns apparent `(RA hours,
/// Dec deg)` for the circumpolar pre-check — under the Hindu bit that
/// must be the *projected* declination the altitude solve uses (binding
/// side). `horizon_alt_deg` is the true center-altitude of the horizon
/// and `lift_max_deg` the maximum refraction + disc lift; both feed only
/// the analytic circumpolar pre-check and the step sizing. `rising`
/// selects the `h_prev < 0 <= h` (rise) versus `h_prev >= 0 > h` (set)
/// crossing polarity, including the equality sides.
///
/// The scan covers [`RISE_SET_MAX_SEARCH_DAYS`] days with a grazing-aware
/// step (10/20/60 minutes); the bracket is refined by
/// [`RISE_SET_BISECT_ITER`] bisection steps with the release rate-aware
/// exit: inside `|h| < 1e-4` the local slope (tight symmetric probe) pins
/// the *time* to ~0.03 s for near-tangential grazing, ordinary crossings
/// take one Newton step onto the true root, and the sub-polar refraction
/// "dip" discontinuity keeps the historic midpoint on the refracted branch
/// above the jump.
#[allow(clippy::too_many_arguments)]
pub fn calculate_rise_set(
    jd_start: f64,
    lat_deg: f64,
    rising: bool,
    event_height: &mut impl FnMut(f64) -> std::result::Result<f64, Error>,
    get_ra_dec: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    horizon_alt_deg: f64,
    lift_max_deg: f64,
) -> Result<RiseSetOutcome> {
    if !jd_start.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: jd_start });
    }
    if !lat_deg.is_finite() {
        return Err(Error::InvalidObserver {
            field: "lat_deg",
            value: lat_deg,
        });
    }
    if !horizon_alt_deg.is_finite() || !lift_max_deg.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: jd_start });
    }
    let height = &mut |jd: f64| -> Result<f64> {
        let h = event_height(jd)?;
        if !h.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        Ok(h)
    };
    let declination = &mut |jd: f64| -> Result<f64> {
        let (_, dec) = get_ra_dec(jd)?;
        if !dec.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        Ok(dec)
    };

    // Circumpolar pre-check from the declination extremes, widened by a
    // finite-difference declination-speed margin for fast movers (the Moon
    // moves ~13 deg/day, so its declination shifts inside the window).
    let dec = declination(jd_start)?;
    let dec2 = declination(jd_start + 1.0 / 24.0)?;
    let dec_speed = (dec2 - dec).abs() * 24.0;
    let margin = dec_speed * RISE_SET_MAX_SEARCH_DAYS;

    let max_alt = 90.0 - (lat_deg - dec).abs();
    let min_alt = if lat_deg >= 0.0 {
        dec - (90.0 - lat_deg)
    } else {
        -dec - (90.0 + lat_deg)
    };
    if min_alt - margin > horizon_alt_deg || max_alt + margin + lift_max_deg < horizon_alt_deg {
        return Ok(RiseSetOutcome::Circumpolar);
    }

    // Step size: near-grazing geometries can stay above/below the horizon
    // for less than an hour, so sample more densely there.
    let grazing_margin = (max_alt - horizon_alt_deg)
        .abs()
        .min((min_alt - horizon_alt_deg).abs());
    let search_step = if grazing_margin < 2.0 {
        10.0 / 1440.0
    } else if grazing_margin < 5.0 {
        20.0 / 1440.0
    } else {
        1.0 / 24.0
    };

    let mut h_prev = height(jd_start)?;
    let mut bracket: Option<(f64, f64)> = None;
    let n_steps = (RISE_SET_MAX_SEARCH_DAYS / search_step) as usize + 1;
    for i in 1..=n_steps {
        let jd = jd_start + i as f64 * search_step;
        let h = height(jd)?;
        let crossed = if rising {
            h_prev < 0.0 && h >= 0.0
        } else {
            h_prev >= 0.0 && h < 0.0
        };
        if crossed {
            bracket = Some((jd - search_step, jd));
            break;
        }
        h_prev = h;
    }
    let (mut lo, mut hi) = match bracket {
        Some(b) => b,
        None => return Ok(RiseSetOutcome::Circumpolar),
    };

    for _ in 0..RISE_SET_BISECT_ITER {
        let jd_mid = 0.5 * (lo + hi);
        let h_mid = height(jd_mid)?;
        if h_mid.abs() < RISE_SET_ALT_TOLERANCE_DEG {
            let span = hi - lo;
            let mut dt = RISE_SET_SLOPE_PROBE_S / 86400.0;
            if dt > 0.25 * span {
                dt = 0.25 * span;
            }
            let mut slope = 0.0;
            if dt > 0.0 {
                let slope_signed =
                    (height(jd_mid + dt)? - height(jd_mid - dt)?) / (2.0 * dt * 86400.0);
                slope = slope_signed.abs();
            }
            if slope <= 0.0
                || slope >= RISE_SET_GRAZE_SLOPE_DEG_S
                || h_mid.abs() < slope * RISE_SET_TIME_TOL_S
            {
                // Refine the actual crossing rather than its local slope.
                // Refraction's dip clamp is a visibility jump, whose visible
                // endpoint must be retained instead of solving a smooth zero.
                let jd_ut = if slope > 0.0 {
                    refine_horizon_crossing(
                        height,
                        lo,
                        hi,
                        rising,
                        slope >= RISE_SET_DIP_SLOPE_DEG_S,
                    )?
                } else {
                    jd_mid
                };
                return Ok(RiseSetOutcome::Found { jd_ut });
            }
        }
        if rising {
            if h_mid < 0.0 {
                lo = jd_mid;
            } else {
                hi = jd_mid;
            }
        } else if h_mid >= 0.0 {
            lo = jd_mid;
        } else {
            hi = jd_mid;
        }
    }

    Ok(RiseSetOutcome::Found {
        jd_ut: 0.5 * (lo + hi),
    })
}

/// Refine an already acquired horizon-crossing bracket.
///
/// Epochs and bounds are UT1 Julian days; `height` is the actual refracted
/// selected disc-point altitude relative to its horizon, in degrees. The
/// bracket's crossing direction is known. For a smooth crossing bisection
/// retains the evaluated point of smallest residual and stops at exact zero
/// or adjacent binary64 epochs. For the refraction-dip visibility jump,
/// `visible_boundary` retains the first visible rising epoch or last visible
/// setting epoch, even when the visible branch begins with a zero plateau.
/// The same sixty-step bracket refinement serves both physical cases without
/// changing the horizon model. Callback errors and nonfinite heights propagate.
pub(crate) fn refine_horizon_crossing(
    height: &mut impl FnMut(f64) -> Result<f64>,
    mut lo: f64,
    mut hi: f64,
    rising: bool,
    visible_boundary: bool,
) -> Result<f64> {
    let left = height(lo)?;
    let right = height(hi)?;
    if !left.is_finite() || !right.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: lo });
    }
    let mut best = if left.abs() <= right.abs() {
        (lo, left.abs())
    } else {
        (hi, right.abs())
    };
    if !visible_boundary && best.1 == 0.0 {
        return Ok(best.0);
    }
    for _ in 0..RISE_SET_BISECT_ITER {
        let mid = lo + (hi - lo) / 2.0;
        if mid == lo || mid == hi {
            break;
        }
        let value = height(mid)?;
        if !value.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: mid });
        }
        if value.abs() < best.1 {
            best = (mid, value.abs());
        }
        if value == 0.0 && !visible_boundary {
            return Ok(mid);
        }
        if (value < 0.0) == rising {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Ok(if visible_boundary {
        if rising { hi } else { lo }
    } else {
        best.0
    })
}

/// Retained horizon-frame angles from one independent TT epoch.
///
/// True obliquity and longitude nutation are degrees. A prepared value owns
/// no provider or clock callback, so callers can consume their actual Time
/// before acquiring the separate model clocks used for sidereal rotation.
#[derive(Clone, Copy, Debug)]
pub struct HorizonFrame {
    eps_true_deg: f64,
    dpsi_deg: f64,
}

impl HorizonFrame {
    /// Retain already evaluated true obliquity and longitude nutation.
    ///
    /// Both angles are degrees. This constructor performs no epoch, timescale,
    /// trigonometric or finite-input validation: it preserves the supplied IEEE
    /// values, including signed zero and NaN. It supports callers which have
    /// already consumed their own frame primitives and must retain their error
    /// and non-finite propagation order. Use [`Self::at_tt`] to evaluate a frame
    /// from a TT Julian Day with the checked numerical primitive domains.
    pub const fn from_degrees(eps_true_deg: f64, dpsi_deg: f64) -> Self {
        Self {
            eps_true_deg,
            dpsi_deg,
        }
    }

    /// Return the retained longitude nutation in degrees, including IEEE values.
    /// No frame, epoch, clock or provider is evaluated by this accessor.
    pub const fn nutation_longitude_degrees(&self) -> f64 {
        self.dpsi_deg
    }

    /// Evaluate Vondrák mean obliquity, then IAU 2006/2000A nutation at
    /// `jd_tt`, in that order. The existing primitive domains/errors apply.
    /// No UT conversion, user Delta-T or ephemeris is acquired here.
    pub fn at_tt(jd_tt: f64) -> Result<Self> {
        let eps_mean = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        let (dpsi, deps) = crate::erfa::nut06a(jd_tt)?;
        Ok(Self {
            eps_true_deg: eps_mean + deps.to_degrees(),
            dpsi_deg: dpsi.to_degrees(),
        })
    }

    /// Combine this independent frame with the configured UT1 sidereal clock.
    ///
    /// `jd_ut` is finite UT1 Julian days; `delta` retains checked Delta-T days
    /// and continuity offsets in degrees. The result is unnormalized Greenwich
    /// apparent sidereal hours, before any longitude division, addition or
    /// modulo. The mean-time radians round-trip and equation of equinoxes are
    /// evaluated in that order. Previously evaluated frame angles retain IEEE
    /// propagation; no provider, TT conversion or additional finite-angle check
    /// is performed. Errors are those of the checked mean sidereal primitive.
    pub fn greenwich_apparent_hours(
        &self,
        jd_ut: f64,
        delta: &crate::time::SiderealDeltaT,
    ) -> Result<f64> {
        crate::time::greenwich_apparent_sidereal_time_hours_from_angles(
            jd_ut,
            self.eps_true_deg,
            self.dpsi_deg,
            delta,
        )
    }
}

/// Project already evaluated equatorial angles onto the observer's horizon.
///
/// Hour angle, declination and geographic latitude are radians. The closed-form
/// result is `(azimuth, true_altitude)` in degrees; azimuth is South-westward
/// and normalized to `[0, 360)`. No clock, frame, atmosphere or provider is
/// acquired. With latitude `phi`, declination `delta` and hour angle `H`,
/// `sin(altitude) = sin(phi)sin(delta) + cos(phi)cos(delta)cos(H)`; azimuth is
/// `atan2(sin(H)cos(delta), cos(H)cos(delta)sin(phi) - sin(delta)cos(phi))`.
/// Inputs retain IEEE values: non-finite derived directions propagate NaN.
/// Callers requiring checked finite inputs validate before this primitive.
pub fn horizon_forward_geometry_radians(
    hour_angle_rad: f64,
    declination_rad: f64,
    latitude_rad: f64,
) -> (f64, f64) {
    let sin_alt = latitude_rad.sin() * declination_rad.sin()
        + latitude_rad.cos() * declination_rad.cos() * hour_angle_rad.cos();
    let true_alt_deg = sin_alt.clamp(-1.0, 1.0).asin().to_degrees();
    let y = hour_angle_rad.sin() * declination_rad.cos();
    let x = hour_angle_rad.cos() * declination_rad.cos() * latitude_rad.sin()
        - declination_rad.sin() * latitude_rad.cos();
    let azimuth_deg = y.atan2(x).to_degrees().rem_euclid(360.0);
    (azimuth_deg, true_alt_deg)
}

/// Recover equatorial geometry from already evaluated horizontal angles.
///
/// Azimuth, true altitude and geographic latitude are radians; azimuth uses
/// the South-westward convention. The closed-form result is `(hour_angle,
/// declination)` in degrees. With latitude `phi`, altitude `h` and azimuth `A`,
/// `sin(declination) = sin(phi)sin(h) - cos(phi)cos(h)cos(A)`. Declination uses
/// constant-first asin bounds and retains that bounded sine in the hour-angle
/// reconstruction. At a geographic
/// or celestial pole with cosine magnitude below `1e-10`, hour angle is zero.
/// No clock, frame, refraction or provider is acquired and no input is refused.
/// Unknown declination intermediates remain NaN. Callers own any prior
/// conversion or trigonometric domain errors.
pub fn horizon_reverse_geometry_radians(
    azimuth_rad: f64,
    true_altitude_rad: f64,
    latitude_rad: f64,
) -> (f64, f64) {
    let sin_dec = latitude_rad.sin() * true_altitude_rad.sin()
        - latitude_rad.cos() * true_altitude_rad.cos() * azimuth_rad.cos();
    let sin_dec = sin_dec.clamp(-1.0, 1.0);
    let dec_rad = sin_dec.asin();
    let dec_deg = dec_rad.to_degrees();
    let cos_dec = dec_rad.cos();
    let cos_lat = latitude_rad.cos();
    let sin_lat = latitude_rad.sin();
    let ha_deg = if cos_dec.abs() < 1e-10 || cos_lat.abs() < 1e-10 {
        0.0
    } else {
        let sin_h = azimuth_rad.sin() * true_altitude_rad.cos() / cos_dec;
        let cos_h = (true_altitude_rad.sin() - sin_lat * sin_dec) / (cos_dec * cos_lat);
        sin_h.atan2(cos_h).to_degrees()
    };
    (ha_deg, dec_deg)
}

/// Celestial-to-horizontal conversion with the same composition as
/// release `utils.azalt`.
///
/// Returns `(azimuth_deg, true_altitude_deg, apparent_altitude_deg)` in
/// degrees: the azimuth keeps the release convention (from South,
/// westward: 0 = South, 90 = West, 180 = North, 270 = East, wrapped to
/// `[0, 360)`), the true altitude is geometric (no refraction) and the
/// apparent altitude carries the refraction leg.
/// `tjdut` is the Julian Day in UT1; `flag` selects the input frame
/// (`ECL2HOR` ecliptic lon/lat, anything else equatorial RA/Dec);
/// `lon_deg_east`/`lat_deg`/`alt_m` are the observer triple;
/// `atpress_mbar <= 0.0` estimates the pressure from `alt_m`;
/// `attemp_c` travels into the refraction leg; `xin` is
/// `(x1, x2, dist)` in degrees (distance passes through the rotation
/// and is otherwise unused). `deltat_at` maps UT days to Delta-T days
/// the nutation/obliquity frame epoch rides the
/// bundled default-timescale TT instead, independent of the caller
/// model — the two legs differ by seconds, enough to flip ULP-level
/// search decisions downstream. `offset_lo_deg` /
/// `offset_hi_deg` are the sidereal continuity offsets for the local
/// apparent sidereal time; `lapse_rate_k_m` is the explicit lapse rate
/// for the `refrac_extended` forward branch. Non-finite inputs
/// (including a non-finite Delta-T word) are
/// [`Error::NonFiniteArgument`].
#[allow(clippy::too_many_arguments)]
pub fn azalt_deg(
    tjdut: f64,
    flag: i32,
    lon_deg_east: f64,
    lat_deg: f64,
    alt_m: f64,
    atpress_mbar: f64,
    attemp_c: f64,
    xin: (f64, f64, f64),
    deltat_at: &dyn Fn(f64) -> Result<f64>,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
    lapse_rate_k_m: f64,
) -> Result<(f64, f64, f64)> {
    azalt_with_frame(
        tjdut,
        flag,
        lon_deg_east,
        lat_deg,
        alt_m,
        atpress_mbar,
        attemp_c,
        xin,
        deltat_at,
        offset_lo_deg,
        offset_hi_deg,
        lapse_rate_k_m,
        None,
    )
}

/// Forward horizon conversion with a previously consumed independent frame.
/// Coordinate, epoch, model-clock and atmospheric input domains follow
/// [`azalt_deg`]. `frame` replaces bundled TT/frame acquisition; its previously
/// evaluated angles retain IEEE non-finite propagation instead of an additional
/// finite-angle check. The actual model Delta-T callback and sidereal offsets
/// remain independent and use the same numerical projection. Unordered altitude
/// clamps use constant-first min/max before asin, as for the ordinary route.
#[allow(clippy::too_many_arguments)]
pub fn azalt_deg_with_frame(
    tjdut: f64,
    flag: i32,
    lon_deg_east: f64,
    lat_deg: f64,
    alt_m: f64,
    atpress_mbar: f64,
    attemp_c: f64,
    xin: (f64, f64, f64),
    deltat_at: &dyn Fn(f64) -> Result<f64>,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
    lapse_rate_k_m: f64,
    frame: HorizonFrame,
) -> Result<(f64, f64, f64)> {
    azalt_with_frame(
        tjdut,
        flag,
        lon_deg_east,
        lat_deg,
        alt_m,
        atpress_mbar,
        attemp_c,
        xin,
        deltat_at,
        offset_lo_deg,
        offset_hi_deg,
        lapse_rate_k_m,
        Some(frame),
    )
}

/// Shared checked forward projection, selecting only frame acquisition.
#[allow(clippy::too_many_arguments)]
fn azalt_with_frame(
    tjdut: f64,
    flag: i32,
    lon_deg_east: f64,
    lat_deg: f64,
    alt_m: f64,
    atpress_mbar: f64,
    attemp_c: f64,
    xin: (f64, f64, f64),
    deltat_at: &dyn Fn(f64) -> Result<f64>,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
    lapse_rate_k_m: f64,
    frame: Option<HorizonFrame>,
) -> Result<(f64, f64, f64)> {
    use crate::time::{SiderealDeltaT, local_apparent_sidereal_time_hours};
    for v in [
        tjdut,
        lon_deg_east,
        lat_deg,
        alt_m,
        atpress_mbar,
        attemp_c,
        xin.0,
        xin.1,
        xin.2,
        offset_lo_deg,
        offset_hi_deg,
        lapse_rate_k_m,
    ] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let dt_days = deltat_at(tjdut)?;
    if !dt_days.is_finite() {
        return Err(Error::NonFiniteArgument { value: dt_days });
    }
    // Frame epoch for nutation/obliquity: the bundled default-timescale TT
    // not the caller-model TT.
    let supplied_frame = frame.is_some();
    let frame = match frame {
        Some(frame) => frame,
        None => {
            let jd_tt_frame = crate::time::default_timescale_tt_from_ut1(tjdut)?;
            let (dpsi, deps) = crate::erfa::nut06a(jd_tt_frame)?;
            HorizonFrame {
                eps_true_deg: crate::positions::vondrak_mean_obliquity_deg(jd_tt_frame)?
                    + deps.to_degrees(),
                dpsi_deg: dpsi.to_degrees(),
            }
        }
    };
    let eps_true_deg = frame.eps_true_deg;
    let (ra_deg, dec_deg) = if flag == ECL2HOR {
        let (ra, dec, _) = cotrans(xin.0, xin.1, xin.2, -eps_true_deg);
        (ra, dec)
    } else {
        (xin.0, xin.1)
    };
    let sdt = SiderealDeltaT::try_new(dt_days, offset_lo_deg, offset_hi_deg)?;
    let sidereal = if supplied_frame {
        crate::time::local_apparent_sidereal_time_hours_from_angles
    } else {
        local_apparent_sidereal_time_hours
    };
    let lst_hours = sidereal(tjdut, lon_deg_east, eps_true_deg, frame.dpsi_deg, &sdt)?;
    let ha_deg = (lst_hours * 15.0 - ra_deg).rem_euclid(360.0);
    let ha_rad = ha_deg.to_radians();
    let dec_rad = dec_deg.to_radians();
    let lat_rad = lat_deg.to_radians();
    let (azimuth_deg, true_alt_deg) = horizon_forward_geometry_radians(ha_rad, dec_rad, lat_rad);
    let pressure_mbar = if atpress_mbar <= 0.0 {
        crate::refraction::estimate_pressure_mbar(alt_m)
    } else {
        atpress_mbar
    };
    let (apparent_alt_deg, _) = crate::refraction::refrac_extended_true_to_app(
        true_alt_deg,
        alt_m,
        pressure_mbar,
        attemp_c,
        lapse_rate_k_m,
    );
    if !supplied_frame && !apparent_alt_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: apparent_alt_deg,
        });
    }
    Ok((azimuth_deg, true_alt_deg, apparent_alt_deg))
}

/// Horizontal-to-celestial conversion with the same composition as
/// release `utils.azalt_rev` (not a numerical inverse of [`azalt_deg`] beyond the
/// geometric legs — refraction applied forward is not removed here).
///
/// Returns the celestial `(x1, x2)` pair in degrees in the selected frame
/// (equatorial RA/Dec or ecliptic lon/lat); distance is not recovered.
/// `azimuth_deg` keeps the release input convention (from South, westward).
/// `flag` selects the output frame (`HOR2EQU` equatorial RA/Dec,
/// anything else ecliptic lon/lat via the positive-obliquity
/// `cotrans` rotation); `alt_m` is accepted for `geopos`-triple shape
/// the interface and ignored. At a celestial pole (`|cos(dec)| < 1e-10`) or a
/// geographic pole (`|cos(lat)| < 1e-10`) the hour angle is genuinely
/// undefined and reports `0.0`, so right ascension equals local
/// sidereal time there. All other legs (epoch, nutation, sidereal
/// time, Delta-T/offsets) follow [`azalt_deg`].
#[allow(clippy::too_many_arguments)]
pub fn azalt_rev_deg(
    tjdut: f64,
    flag: i32,
    lon_deg_east: f64,
    lat_deg: f64,
    alt_m: f64,
    azimuth_deg: f64,
    true_alt_deg: f64,
    deltat_at: &dyn Fn(f64) -> Result<f64>,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
) -> Result<(f64, f64)> {
    azalt_rev_with_frame(
        tjdut,
        flag,
        lon_deg_east,
        lat_deg,
        alt_m,
        azimuth_deg,
        true_alt_deg,
        deltat_at,
        offset_lo_deg,
        offset_hi_deg,
        None,
    )
}

/// Reverse horizon conversion with a previously consumed independent frame.
/// Units, pole behavior and output-frame selection follow [`azalt_rev_deg`].
/// `frame` owns only independent obliquity/nutation; the supplied model clock
/// still supplies UT1 sidereal time and its continuity corrections. Supplied
/// angle NaNs propagate through sidereal time and right ascension; they are not
/// refused after the caller has already evaluated its frame. The ordinary
/// [`azalt_rev_deg`] entry retains checked bundled frame acquisition.
#[allow(clippy::too_many_arguments)]
pub fn azalt_rev_deg_with_frame(
    tjdut: f64,
    flag: i32,
    lon_deg_east: f64,
    lat_deg: f64,
    alt_m: f64,
    azimuth_deg: f64,
    true_alt_deg: f64,
    deltat_at: &dyn Fn(f64) -> Result<f64>,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
    frame: HorizonFrame,
) -> Result<(f64, f64)> {
    azalt_rev_with_frame(
        tjdut,
        flag,
        lon_deg_east,
        lat_deg,
        alt_m,
        azimuth_deg,
        true_alt_deg,
        deltat_at,
        offset_lo_deg,
        offset_hi_deg,
        Some(frame),
    )
}

/// Shared checked reverse projection, selecting only frame acquisition.
#[allow(clippy::too_many_arguments)]
fn azalt_rev_with_frame(
    tjdut: f64,
    flag: i32,
    lon_deg_east: f64,
    lat_deg: f64,
    alt_m: f64,
    azimuth_deg: f64,
    true_alt_deg: f64,
    deltat_at: &dyn Fn(f64) -> Result<f64>,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
    frame: Option<HorizonFrame>,
) -> Result<(f64, f64)> {
    use crate::time::{SiderealDeltaT, local_apparent_sidereal_time_hours};
    for v in [
        tjdut,
        lon_deg_east,
        lat_deg,
        alt_m,
        azimuth_deg,
        true_alt_deg,
        offset_lo_deg,
        offset_hi_deg,
    ] {
        if !v.is_finite() {
            return Err(Error::NonFiniteArgument { value: v });
        }
    }
    let _ = alt_m;
    let dt_days = deltat_at(tjdut)?;
    if !dt_days.is_finite() {
        return Err(Error::NonFiniteArgument { value: dt_days });
    }
    // Frame epoch for nutation/obliquity: the bundled default-timescale TT
    // not the caller-model TT.
    let supplied_frame = frame.is_some();
    let frame = match frame {
        Some(frame) => frame,
        None => {
            let jd_tt_frame = crate::time::default_timescale_tt_from_ut1(tjdut)?;
            let (dpsi, deps) = crate::erfa::nut06a(jd_tt_frame)?;
            HorizonFrame {
                eps_true_deg: crate::positions::vondrak_mean_obliquity_deg(jd_tt_frame)?
                    + deps.to_degrees(),
                dpsi_deg: dpsi.to_degrees(),
            }
        }
    };
    let eps_true_deg = frame.eps_true_deg;
    let az_rad = azimuth_deg.to_radians();
    let alt_rad = true_alt_deg.to_radians();
    let lat_rad = lat_deg.to_radians();
    let (ha_deg, dec_deg) = horizon_reverse_geometry_radians(az_rad, alt_rad, lat_rad);
    let sdt = SiderealDeltaT::try_new(dt_days, offset_lo_deg, offset_hi_deg)?;
    let sidereal = if supplied_frame {
        crate::time::local_apparent_sidereal_time_hours_from_angles
    } else {
        local_apparent_sidereal_time_hours
    };
    let lst_hours = sidereal(tjdut, lon_deg_east, eps_true_deg, frame.dpsi_deg, &sdt)?;
    let ra_deg = (lst_hours * 15.0 - ha_deg).rem_euclid(360.0);
    if flag == HOR2EQU {
        Ok((ra_deg, dec_deg))
    } else {
        let (lon, lat, _) = cotrans(ra_deg, dec_deg, 1.0, eps_true_deg);
        Ok((lon, lat))
    }
}

/// Apparent disc radii in degrees from geocentric distances. Distances must be
/// positive AU words; non-finite words are refused.
pub fn eclipse_disc_angular_radii_deg(sun_dist_au: f64, moon_dist_au: f64) -> Result<(f64, f64)> {
    for value in [sun_dist_au, moon_dist_au] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    Ok((
        (SUN_SEMIDIAMETER_ARCSEC / 3600.0) / sun_dist_au,
        (MOON_SEMIDIAMETER_ARCSEC / 3600.0) * (MOON_SEMIDIAMETER_REF_AU / moon_dist_au),
    ))
}

impl<'a, P> SealedTransitServing<'a, P>
where
    P: StateProvider,
{
    /// Checked constructor: borrows the sealed serving pair and the
    /// Delta-T provider, resolving the sidereal-continuity offsets from
    /// the provider's own boundary words. No I/O, no globals.
    pub fn try_new(
        serving: &'a SealedServing<'a, P>,
        deltat_at: &'a dyn Fn(f64) -> Result<f64>,
        body_id: i32,
        calc_flags: i32,
    ) -> Result<Self> {
        use crate::time::{
            SIDEREAL_MODERN_END_JD, SIDEREAL_MODERN_START_JD, sidereal_continuity_offsets,
        };
        let dt_lo = deltat_at(SIDEREAL_MODERN_START_JD)?;
        if !dt_lo.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_lo });
        }
        let dt_hi = deltat_at(SIDEREAL_MODERN_END_JD)?;
        if !dt_hi.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_hi });
        }
        let (offset_lo_deg, offset_hi_deg) = sidereal_continuity_offsets(dt_lo, dt_hi)?;
        Ok(Self {
            serving,
            deltat_at,
            body_id,
            calc_flags,
            offset_lo_deg,
            offset_hi_deg,
        })
    }

    /// Bound sealed serving pair.
    pub fn serving(&self) -> &SealedServing<'a, P> {
        self.serving
    }

    /// Bound body id.
    pub fn body_id(&self) -> i32 {
        self.body_id
    }

    /// Bound base calculation flags (without the implied
    /// `FLG_EQUATORIAL | FLG_SPEED`).
    pub fn calc_flags(&self) -> i32 {
        self.calc_flags
    }

    /// TT epoch for a UT epoch through the bound Delta-T provider.
    pub fn jd_tt_for_ut(&self, jd_ut: f64) -> Result<f64> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt = (self.deltat_at)(jd_ut)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_ut + dt)
    }

    /// Apparent `(RA hours, Dec deg)` at a UT epoch through the sealed
    /// equatorial words.
    pub fn ra_dec_at_ut(&self, jd_ut: f64) -> Result<(f64, f64)> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (pos, _) = serving_compute_tt(
            self.serving,
            jd_tt,
            self.body_id,
            self.calc_flags | crate::constants::FLG_EQUATORIAL | crate::constants::FLG_SPEED,
        )?;
        Ok((pos[0] / 15.0, pos[1]))
    }

    /// Greenwich apparent sidereal time in hours `[0, 24)` at a UT epoch
    /// the
    /// native sidereal-time composition over the sealed nutation pair at the
    /// TT epoch and the Vondrak mean obliquity the frame pipeline uses.
    pub fn gast_hours_at_ut(&self, jd_ut: f64) -> Result<f64> {
        use crate::time::{SiderealDeltaT, greenwich_apparent_sidereal_time_hours};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt_days = (self.deltat_at)(jd_ut)?;
        if !dt_days.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_days });
        }
        let jd_tt = jd_ut + dt_days;
        let (dpsi_rad, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let dpsi_deg = dpsi_rad.to_degrees();
        let eps_mean_deg = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        let eps_true_deg = eps_mean_deg + deps_rad.to_degrees();
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        greenwich_apparent_sidereal_time_hours(jd_ut, eps_true_deg, dpsi_deg, &sdt)
    }

    /// Meridian transit over the sealed serving.
    /// `lower_transit` selects the anti-culmination (`CALC_ITRANSIT`
    /// 180-degree shift).
    pub fn transit_at_ut(
        &self,
        jd_start: f64,
        lon_deg_east: f64,
        lower_transit: bool,
    ) -> Result<TransitSolution> {
        if !jd_start.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd_start });
        }
        let mut get_ra_dec = |jd: f64| self.ra_dec_at_ut(jd);
        let mut gast_hours = |jd: f64| self.gast_hours_at_ut(jd);
        calculate_transit(
            jd_start,
            lon_deg_east,
            lower_transit,
            &mut get_ra_dec,
            &mut gast_hours,
        )
    }

    /// Local apparent sidereal time in hours `[0, 24)` at a UT epoch for
    /// an observer longitude. This is the
    /// `azalt` LST acquisition the sealed rise/set altitude engine uses.
    pub fn last_hours_at_ut(&self, jd_ut: f64, lon_deg_east: f64) -> Result<f64> {
        use crate::time::{SiderealDeltaT, local_apparent_sidereal_time_hours};
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        if !lon_deg_east.is_finite() {
            return Err(Error::InvalidObserver {
                field: "lon_deg",
                value: lon_deg_east,
            });
        }
        let dt_days = (self.deltat_at)(jd_ut)?;
        if !dt_days.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt_days });
        }
        let jd_tt = jd_ut + dt_days;
        let (dpsi_rad, deps_rad) = (self.serving.nutation_provider())(jd_tt)?;
        let dpsi_deg = dpsi_rad.to_degrees();
        let eps_mean_deg = crate::positions::vondrak_mean_obliquity_deg(jd_tt)?;
        let eps_true_deg = eps_mean_deg + deps_rad.to_degrees();
        let sdt = SiderealDeltaT::try_new(dt_days, self.offset_lo_deg, self.offset_hi_deg)?;
        local_apparent_sidereal_time_hours(jd_ut, lon_deg_east, eps_true_deg, dpsi_deg, &sdt)
    }
}

/// Meridian-transit estimate as a UT1 Julian day over the bound serving.
///
/// `lon_deg_east` is east-positive longitude in degrees. `lower_transit=false`
/// selects upper culmination; `true` selects anti-culmination. XYZ/RADIANS are
/// refused, and consumed provider/clock errors propagate. This convenience
/// entry discards [`TransitSolution::converged`]: reaching the iteration cap
/// can still return the last estimate. Use the serving's `transit_at_ut` method
/// when the convergence verdict is required.
pub fn transit_ut_serving<P: StateProvider>(
    serving: &SealedTransitServing<'_, P>,
    jd_start: f64,
    lon_deg_east: f64,
    lower_transit: bool,
) -> Result<f64> {
    check_angular_serving_flags(serving.calc_flags(), "transit_ut")?;
    Ok(serving
        .transit_at_ut(jd_start, lon_deg_east, lower_transit)?
        .jd_ut)
}

/// Physical disc radius of a classical body in AU, using the table's
/// volumetric-mean radii (with the separately defined Sun/Moon constants).
/// These are not uniformly equatorial radii. Bodies without a finite apparent disc
/// on this path map to `0.0` (no disc term) — the release lookup misses
/// every minor id, so Chiron and the asteroids are point sources here.
pub fn body_disc_radius_au(body_id: i32) -> f64 {
    use crate::constants::{
        JUPITER, MARS, MERCURY, MOON, NEPTUNE, PLUTO, SATURN, SUN, URANUS, VENUS,
    };
    if body_id == SUN {
        return SUN_RADIUS_AU;
    }
    if body_id == MOON {
        return MOON_RADIUS_AU;
    }
    let radius_km = match body_id {
        _ if body_id == MERCURY => 2439.4,
        _ if body_id == VENUS => 6051.8,
        _ if body_id == MARS => 3389.5,
        _ if body_id == JUPITER => 69911.0,
        _ if body_id == SATURN => 58232.0,
        _ if body_id == URANUS => 25362.0,
        _ if body_id == NEPTUNE => 24622.0,
        _ if body_id == PLUTO => 1188.3,
        _ => return 0.0,
    };
    radius_km / RISE_SET_AU_KM
}

/// Next rise/set over the sealed serving: `(res, tret)` with the release
/// `rise_trans_true_hor` contract — `res = 0` carries the event JD (UT)
/// in `tret[0]`, `res = -2` is the circumpolar/exhausted verdict with
/// `tret[0] = 0.0`. `transit` supplies the sealed ephemeris, Delta-T and
/// sidereal words (its `calc_flags` travel with `FLG_TOPOCTR | FLG_SPEED`
/// added, the release `calc_ut` acquisition spelling).
pub fn rise_set_ut_serving<P: StateProvider>(
    transit: &SealedTransitServing<'_, P>,
    req: &RiseSetRequest,
) -> Result<(i32, [f64; 10])> {
    use crate::constants::{
        BIT_FIXED_DISC_SIZE, BIT_NO_REFRACTION, EARTH, ECL_NUT, FLG_SPEED, FLG_TOPOCTR, MOON, SUN,
    };
    let body_id = transit.body_id();
    // FIX-C37 (follows FIX-C34): the altitude engine reads angular
    // RA/Dec words from this serving, so format-flagged servings are
    // refused here too — not only on the transit-culmination path.
    check_angular_serving_flags(transit.calc_flags(), "rise_set_ut")?;
    if body_id == EARTH || body_id == ECL_NUT {
        // Neither id denotes a direction in the sky (the observer itself /
        // a nutation carrier): the release fails explicitly rather than
        // manufacturing a location-dependent event from non-directional
        // data.
        return Err(Error::InvalidBody {
            body_id,
            operation: "rise_set_ut",
        });
    }
    let rising = match resolve_event_kind(req.rsmi) {
        RiseSetEvent::Rise => true,
        RiseSetEvent::Set => false,
        _ => {
            return Err(Error::InvalidSearchStep {
                value: req.rsmi as f64,
            });
        }
    };
    validate_geopos(req.lon_deg_east, req.lat_deg, req.height_m)?;
    if !req.jd_start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: req.jd_start_ut });
    }
    // `atpress == 0`: the release barometric estimate (never "no
    // pressure"); the temperature travels.
    let pressure = if req.atpress_mbar == 0.0 {
        pressure_from_altitude_m(req.height_m)
    } else {
        req.atpress_mbar
    };
    // The `-100` sentinel dips the sea horizon at the observer altitude
    // under the release lapse default, plus the 0.0001 epsilon; any other
    // horizon (including negatives) travels.
    let mut horizon_alt = resolve_horizon_altitude_deg(
        req.horhgt_deg,
        crate::refraction::calc_dip(
            req.height_m.max(0.0),
            RISE_SET_DIP_LAPSE_RATE_K_M,
            pressure,
            req.attemp_c,
        ),
    );
    // Sun twilight depressions force the geometric-center, no-refraction
    // combination (planets and the Moon keep their ordinary crossing).
    let mut rsmi_eff = req.rsmi;
    if let Some(twilight) = twilight_override(req.rsmi, body_id == SUN) {
        rsmi_eff |= crate::constants::BIT_NO_REFRACTION | crate::constants::BIT_DISC_CENTER;
        horizon_alt = twilight.horizon_alt_deg;
    }
    let use_refraction = rsmi_eff & BIT_NO_REFRACTION == 0;
    let disc = select_disc_point(rsmi_eff, false);
    let body_radius_au = body_disc_radius_au(body_id);
    let mut lift_max = 0.0;
    if use_refraction {
        lift_max += 0.75;
    }
    if disc != DiscPoint::Center {
        lift_max += 0.30;
    }
    // Hindu-rising convention (`BIT_GEOCTR_NO_ECL_LAT`, the
    // latitude-zeroing component of `BIT_HINDU_RISING`): the release
    // replaces the ordinary topocentric place with the body's geocentric
    // apparent place projected onto the ecliptic (ecliptic latitude set
    // to zero, longitude and distance kept), rotated to the observer's
    // horizon. This drops both the topocentric parallax and the ecliptic
    // latitude, so the event moves (e.g. the Moon shifts by tens of
    // minutes toward its ecliptic-longitude crossing). Meridian transits
    // ignore the bit (the slice-1 transit path returns before the
    // altitude engine), so only this rise/set path branches on it.
    let hindu = rsmi_eff & crate::constants::BIT_GEOCTR_NO_ECL_LAT != 0;
    let serving = transit.serving();
    let calc_flags = transit.calc_flags() | FLG_TOPOCTR | FLG_SPEED;
    // Geocentric acquisition spelling under the Hindu bit: no topocentric site, no observer
    // scope — the horizon rotation still uses the observer's location.
    let geo_flags = transit.calc_flags() | FLG_SPEED;
    let site = (req.lon_deg_east, req.lat_deg, req.height_m);
    // Native celestial-to-terrestrial matrix at each sample epoch with
    // `xp = yp = 0`, exactly the release fast-calc `c2t06a` call.
    let c2t_at = |tt: f64, ut1: f64| crate::erfa::c2t06a(tt, ut1, 0.0, 0.0);
    // True altitude plus distance from the apparent place rotated to the
    // observer's horizon.
    let alt_dist_at_ut = |jd_ut: f64| -> Result<(f64, f64)> {
        let jd_tt = transit.jd_tt_for_ut(jd_ut)?;
        let (pos, _) = if hindu {
            crate::positions::serving_compute_tt(serving, jd_tt, body_id, geo_flags)?
        } else {
            crate::positions::serving_compute_tt_topocentric(
                serving,
                jd_tt,
                jd_ut,
                body_id,
                calc_flags,
                site,
                &c2t_at,
                false,
                crate::positions::SiteEllipsoid::Wgs84,
            )?
        };
        let (_, deps_rad) = (serving.nutation_provider())(jd_tt)?;
        let eps_true_deg =
            crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
        let (ecl_lon, ecl_lat) = if hindu {
            (pos[0], 0.0)
        } else {
            (pos[0], pos[1])
        };
        let (ra_deg, dec_deg, dist) = cotrans(ecl_lon, ecl_lat, pos[2], -eps_true_deg);
        let last_hours = transit.last_hours_at_ut(jd_ut, req.lon_deg_east)?;
        let ha_deg = (last_hours * 15.0 - ra_deg).rem_euclid(360.0);
        let ha_rad = ha_deg.to_radians();
        let dec_rad = dec_deg.to_radians();
        let lat_rad = req.lat_deg.to_radians();
        let sin_alt = lat_rad.sin() * dec_rad.sin() + lat_rad.cos() * dec_rad.cos() * ha_rad.cos();
        Ok((sin_alt.clamp(-1.0, 1.0).asin().to_degrees(), dist))
    };
    // Height of the configured disc point above the configured horizon
    // .
    let dip_deg = crate::refraction::calc_dip(
        req.height_m.max(0.0),
        RISE_SET_DIP_LAPSE_RATE_K_M,
        pressure,
        req.attemp_c,
    );
    let mut event_height = |jd: f64| -> Result<f64> {
        let (alt_true, dist) = alt_dist_at_ut(jd)?;
        let mut cur = dist;
        if rsmi_eff & BIT_FIXED_DISC_SIZE != 0 {
            if body_id == SUN {
                cur = FIXED_DISC_SUN_DIST_AU;
            } else if body_id == MOON {
                cur = FIXED_DISC_MOON_DIST_AU;
            }
        }
        let mut alt_pt = alt_true + disc_altitude_term_deg(disc, body_radius_au, cur)?;
        if use_refraction {
            let mut refraction_deg = |app: f64| -> Result<f64> {
                Ok(crate::refraction::sinclair_refraction_deg(
                    app,
                    pressure,
                    req.attemp_c,
                ))
            };
            alt_pt = true_to_apparent_altitude(alt_pt, dip_deg, &mut refraction_deg, jd)?;
        }
        Ok(alt_pt - horizon_alt)
    };
    // The circumpolar pre-check sees the same declination the altitude
    // solve uses: the geocentric apparent declination (exactly the
    // transit serving's words) on the ordinary path, or the projected
    // declination of the latitude-zeroed geocentric place under the Hindu
    // bit.
    let mut get_ra_dec = |jd: f64| -> Result<(f64, f64)> {
        if hindu {
            let jd_tt = transit.jd_tt_for_ut(jd)?;
            let (pos, _) =
                crate::positions::serving_compute_tt(serving, jd_tt, body_id, geo_flags)?;
            let (_, deps_rad) = (serving.nutation_provider())(jd_tt)?;
            let eps_true_deg =
                crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
            let (ra_deg, dec_deg, _) = cotrans(pos[0], 0.0, pos[2], -eps_true_deg);
            Ok((ra_deg / 15.0, dec_deg))
        } else {
            transit.ra_dec_at_ut(jd)
        }
    };
    match calculate_rise_set(
        req.jd_start_ut,
        req.lat_deg,
        rising,
        &mut event_height,
        &mut get_ra_dec,
        horizon_alt,
        lift_max,
    )? {
        RiseSetOutcome::Found { jd_ut } => Ok((0, make_tret(jd_ut))),
        RiseSetOutcome::Circumpolar => Ok((-2, make_tret(0.0))),
    }
}

/// Next meridian transit JD (UT) of the bound star over the sealed star
/// serving: `lower_transit` `false` is the upper culmination
/// (`CALC_MTRANSIT`), `true` the anti-culmination (`CALC_ITRANSIT`).
pub fn star_transit_ut_serving<P>(
    serving: &SealedStarEventServing<'_, P>,
    jd_start: f64,
    lon_deg_east: f64,
    lower_transit: bool,
) -> Result<f64>
where
    P: crate::positions::StateProvider,
{
    Ok(serving
        .star_transit_at_ut(jd_start, lon_deg_east, lower_transit)?
        .jd_ut)
}

/// Next star rise/set over the sealed star serving: `(res, tret)` with
/// the release `rise_trans_true_hor` contract — `res = 0` carries the
/// event JD (UT) in `tret[0]`, `res = -2` is the circumpolar/exhausted
/// verdict with `tret[0] = 0.0`.
pub fn star_rise_set_ut_serving<P>(
    serving: &SealedStarEventServing<'_, P>,
    req: &StarRiseSetRequest,
) -> Result<(i32, [f64; 10])>
where
    P: crate::positions::StateProvider,
{
    use crate::constants::BIT_NO_REFRACTION;
    let rising = match resolve_event_kind(req.rsmi) {
        RiseSetEvent::Rise => true,
        RiseSetEvent::Set => false,
        _ => {
            return Err(Error::InvalidSearchStep {
                value: req.rsmi as f64,
            });
        }
    };
    validate_geopos(req.lon_deg_east, req.lat_deg, req.height_m)?;
    if !req.jd_start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: req.jd_start_ut });
    }
    // `atpress == 0`: the release barometric estimate (never "no
    // pressure"); the temperature travels.
    let pressure = if req.atpress_mbar == 0.0 {
        pressure_from_altitude_m(req.height_m)
    } else {
        req.atpress_mbar
    };
    // The `-100` sentinel dips the sea horizon at the observer altitude
    // under the release lapse default, plus the 0.0001 epsilon; any other
    // horizon (including negatives) travels.
    let mut horizon_alt = resolve_horizon_altitude_deg(
        req.horhgt_deg,
        crate::refraction::calc_dip(
            req.height_m.max(0.0),
            RISE_SET_DIP_LAPSE_RATE_K_M,
            pressure,
            req.attemp_c,
        ),
    );
    // Star twilight depressions force the geometric-center, no-refraction
    // combination.
    let mut rsmi_eff = req.rsmi;
    if let Some(twilight) = twilight_override(req.rsmi, true) {
        rsmi_eff |= crate::constants::BIT_NO_REFRACTION | crate::constants::BIT_DISC_CENTER;
        horizon_alt = twilight.horizon_alt_deg;
    }
    let use_refraction = rsmi_eff & BIT_NO_REFRACTION == 0;
    // Fixed stars always use the disc center; the disc modifier bits are accepted and ignored.
    let disc = select_disc_point(rsmi_eff, true);
    debug_assert!(disc == DiscPoint::Center);
    let mut lift_max = 0.0;
    if use_refraction {
        lift_max += 0.75;
    }
    // Hindu-rising convention for stars: the geocentric star place with
    // the ecliptic latitude zeroed, rotated to the observer's horizon
    // .
    let hindu = rsmi_eff & crate::constants::BIT_GEOCTR_NO_ECL_LAT != 0;
    // True altitude plus geocentric distance from the sealed star place
    // rotated to the observer's horizon.
    let alt_dist_at_ut = |jd_ut: f64| -> Result<(f64, f64)> {
        let triple = serving.star_ecliptic_at_ut(jd_ut)?;
        let jd_tt = serving.jd_tt_for_ut(jd_ut)?;
        let (_, deps_rad) = (serving.nutation_provider())(jd_tt)?;
        let eps_true_deg =
            crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
        let (ecl_lon, ecl_lat) = if hindu {
            (triple[0], 0.0)
        } else {
            (triple[0], triple[1])
        };
        let (ra_deg, dec_deg, dist) = cotrans(ecl_lon, ecl_lat, triple[2], -eps_true_deg);
        let last_hours = serving.last_hours_at_ut(jd_ut, req.lon_deg_east)?;
        let ha_deg = (last_hours * 15.0 - ra_deg).rem_euclid(360.0);
        let ha_rad = ha_deg.to_radians();
        let dec_rad = dec_deg.to_radians();
        let lat_rad = req.lat_deg.to_radians();
        let sin_alt = lat_rad.sin() * dec_rad.sin() + lat_rad.cos() * dec_rad.cos() * ha_rad.cos();
        Ok((sin_alt.clamp(-1.0, 1.0).asin().to_degrees(), dist))
    };
    // Height of the disc center above the configured horizon.
    let dip_deg = crate::refraction::calc_dip(
        req.height_m.max(0.0),
        RISE_SET_DIP_LAPSE_RATE_K_M,
        pressure,
        req.attemp_c,
    );
    let mut event_height = |jd: f64| -> Result<f64> {
        let (alt_true, dist) = alt_dist_at_ut(jd)?;
        let mut alt_pt = alt_true + disc_altitude_term_deg(disc, 0.0, dist)?;
        if use_refraction {
            let mut refraction_deg = |app: f64| -> Result<f64> {
                Ok(crate::refraction::sinclair_refraction_deg(
                    app,
                    pressure,
                    req.attemp_c,
                ))
            };
            alt_pt = true_to_apparent_altitude(alt_pt, dip_deg, &mut refraction_deg, jd)?;
        }
        Ok(alt_pt - horizon_alt)
    };
    // The circumpolar pre-check sees the same declination the altitude
    // solve uses: the star apparent declination on the ordinary path, or
    // the projected declination of the latitude-zeroed place under the
    // Hindu bit.
    let mut get_ra_dec = |jd: f64| -> Result<(f64, f64)> {
        if hindu {
            let triple = serving.star_ecliptic_at_ut(jd)?;
            let jd_tt = serving.jd_tt_for_ut(jd)?;
            let (_, deps_rad) = (serving.nutation_provider())(jd_tt)?;
            let eps_true_deg =
                crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
            let (ra_deg, dec_deg, _) = cotrans(triple[0], 0.0, triple[2], -eps_true_deg);
            Ok((ra_deg / 15.0, dec_deg))
        } else {
            serving.star_ra_dec_at_ut(jd)
        }
    };
    match calculate_rise_set(
        req.jd_start_ut,
        req.lat_deg,
        rising,
        &mut event_height,
        &mut get_ra_dec,
        horizon_alt,
        lift_max,
    )? {
        RiseSetOutcome::Found { jd_ut } => Ok((0, make_tret(jd_ut))),
        RiseSetOutcome::Circumpolar => Ok((-2, make_tret(0.0))),
    }
}

/// Next lunar-point meridian transit JD (UT) over the sealed serving:
/// `lower_transit` `false` is the upper culmination (`CALC_MTRANSIT`),
/// `true` the anti-culmination (`CALC_ITRANSIT`).
pub fn lunar_transit_ut_serving<P: crate::positions::StateProvider>(
    serving: &SealedLunarEventServing<'_, P>,
    jd_start: f64,
    lon_deg_east: f64,
    lower_transit: bool,
) -> Result<f64> {
    Ok(serving
        .lunar_transit_at_ut(jd_start, lon_deg_east, lower_transit)?
        .jd_ut)
}

/// Next lunar-point rise/set over the sealed serving: `(res, tret)` with
/// the release `rise_trans_true_hor` contract — `res = 0` carries the
/// event JD (UT) in `tret[0]`, `res = -2` is the circumpolar/exhausted
/// verdict with `tret[0] = 0.0`.
///
/// Mirrors [`rise_set_ut_serving`] arm for arm (geopos unpack, event-kind
/// gate, upfront placeability probe, pressure estimate, `-100` dip
/// horizon, disc selection with the zero-radius point-source term,
/// Sinclair-via-transit-core event height, `calculate_rise_set`)
/// over the geocentric lunar words. The twilight override is inert on
/// this path;
/// the Hindu `BIT_GEOCTR_NO_ECL_LAT` projection applies to the
/// geocentric words with its projected-declination circumpolar
/// pre-check, like the release.
pub fn lunar_rise_set_ut_serving<P: crate::positions::StateProvider>(
    serving: &SealedLunarEventServing<'_, P>,
    req: &RiseSetRequest,
) -> Result<(i32, [f64; 10])> {
    use crate::constants::{BIT_FIXED_DISC_SIZE, BIT_NO_REFRACTION, MOON, SUN};
    let rising = match resolve_event_kind(req.rsmi) {
        RiseSetEvent::Rise => true,
        RiseSetEvent::Set => false,
        _ => {
            return Err(Error::InvalidSearchStep {
                value: req.rsmi as f64,
            });
        }
    };
    validate_geopos(req.lon_deg_east, req.lat_deg, req.height_m)?;
    if !req.jd_start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: req.jd_start_ut });
    }
    // Upfront placeability probe: a coverage miss keeps its typed
    // range error, any other miss already refused at construction, so a
    // readable word here means the search runs on served words.
    let jd_probe = serving.jd_tt_for_ut(req.jd_start_ut)?;
    serving.lunar_ecliptic_at_tt(jd_probe)?;
    let pressure = if req.atpress_mbar == 0.0 {
        pressure_from_altitude_m(req.height_m)
    } else {
        req.atpress_mbar
    };
    let horizon_alt = resolve_horizon_altitude_deg(
        req.horhgt_deg,
        crate::refraction::calc_dip(
            req.height_m.max(0.0),
            RISE_SET_DIP_LAPSE_RATE_K_M,
            pressure,
            req.attemp_c,
        ),
    );
    let rsmi_eff = req.rsmi;
    let use_refraction = rsmi_eff & BIT_NO_REFRACTION == 0;
    let disc = select_disc_point(rsmi_eff, false);
    let body_radius_au = body_disc_radius_au(serving.body_id());
    debug_assert_eq!(body_radius_au, 0.0);
    let mut lift_max = 0.0;
    if use_refraction {
        lift_max += 0.75;
    }
    if disc != DiscPoint::Center {
        lift_max += 0.30;
    }
    let hindu = rsmi_eff & crate::constants::BIT_GEOCTR_NO_ECL_LAT != 0;
    // True altitude plus distance from the geocentric lunar place rotated
    // to the observer's horizon.
    let alt_dist_at_ut = |jd_ut: f64| -> Result<(f64, f64)> {
        let jd_tt = serving.jd_tt_for_ut(jd_ut)?;
        let (lon, lat, dist) = serving.lunar_ecliptic_at_tt(jd_tt)?;
        let (_, deps_rad) = (serving.serving.nutation_provider())(jd_tt)?;
        let eps_true_deg =
            crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
        let (ecl_lon, ecl_lat) = if hindu { (lon, 0.0) } else { (lon, lat) };
        let (ra_deg, dec_deg, dist) = cotrans(ecl_lon, ecl_lat, dist, -eps_true_deg);
        let last_hours = serving.last_hours_at_ut(jd_ut, req.lon_deg_east)?;
        let ha_deg = (last_hours * 15.0 - ra_deg).rem_euclid(360.0);
        let ha_rad = ha_deg.to_radians();
        let dec_rad = dec_deg.to_radians();
        let lat_rad = req.lat_deg.to_radians();
        let sin_alt = lat_rad.sin() * dec_rad.sin() + lat_rad.cos() * dec_rad.cos() * ha_rad.cos();
        Ok((sin_alt.clamp(-1.0, 1.0).asin().to_degrees(), dist))
    };
    // Height of the configured disc point above the configured horizon
    // .
    let dip_deg = crate::refraction::calc_dip(
        req.height_m.max(0.0),
        RISE_SET_DIP_LAPSE_RATE_K_M,
        pressure,
        req.attemp_c,
    );
    let mut event_height = |jd: f64| -> Result<f64> {
        let (alt_true, dist) = alt_dist_at_ut(jd)?;
        let mut cur = dist;
        if rsmi_eff & BIT_FIXED_DISC_SIZE != 0 {
            if serving.body_id() == SUN {
                cur = FIXED_DISC_SUN_DIST_AU;
            } else if serving.body_id() == MOON {
                cur = FIXED_DISC_MOON_DIST_AU;
            }
        }
        let mut alt_pt = alt_true + disc_altitude_term_deg(disc, body_radius_au, cur)?;
        if use_refraction {
            let mut refraction_deg = |app: f64| -> Result<f64> {
                Ok(crate::refraction::sinclair_refraction_deg(
                    app,
                    pressure,
                    req.attemp_c,
                ))
            };
            alt_pt = true_to_apparent_altitude(alt_pt, dip_deg, &mut refraction_deg, jd)?;
        }
        Ok(alt_pt - horizon_alt)
    };
    // The circumpolar pre-check sees the same declination the altitude
    // solve uses: the lunar apparent declination on the ordinary path, or
    // the projected declination of the latitude-zeroed place under the
    // Hindu bit.
    let mut get_ra_dec = |jd: f64| -> Result<(f64, f64)> {
        if hindu {
            let jd_tt = serving.jd_tt_for_ut(jd)?;
            let (lon, _, dist) = serving.lunar_ecliptic_at_tt(jd_tt)?;
            let (_, deps_rad) = (serving.serving.nutation_provider())(jd_tt)?;
            let eps_true_deg =
                crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
            let (ra_deg, dec_deg, _) = cotrans(lon, 0.0, dist, -eps_true_deg);
            Ok((ra_deg / 15.0, dec_deg))
        } else {
            serving.lunar_ra_dec_at_ut(jd)
        }
    };
    match calculate_rise_set(
        req.jd_start_ut,
        req.lat_deg,
        rising,
        &mut event_height,
        &mut get_ra_dec,
        horizon_alt,
        lift_max,
    )? {
        RiseSetOutcome::Found { jd_ut } => Ok((0, make_tret(jd_ut))),
        RiseSetOutcome::Circumpolar => Ok((-2, make_tret(0.0))),
    }
}

/// Next minor-body meridian transit JD (UT) over the sealed serving:
/// `lower_transit` `false` is the upper culmination (`CALC_MTRANSIT`),
/// `true` the anti-culmination (`CALC_ITRANSIT`). LEB-channeled bodies
/// run the shared transit serving; the channel-less remainder
/// runs the fallback transit over the Keplerian words.
pub fn minor_transit_ut_serving<P: crate::positions::StateProvider>(
    serving: &SealedMinorEventServing<'_, P>,
    jd_start: f64,
    lon_deg_east: f64,
    lower_transit: bool,
) -> Result<f64> {
    if serving.is_leb_channeled() {
        transit_ut_serving(
            &serving.inner_transit()?,
            jd_start,
            lon_deg_east,
            lower_transit,
        )
    } else {
        Ok(serving
            .kepler_transit_at_ut(jd_start, lon_deg_east, lower_transit)?
            .jd_ut)
    }
}

/// Next minor-body rise/set over the sealed serving: `(res, tret)` with
/// the release `rise_trans_true_hor` contract — `res = 0` carries the
/// event JD (UT) in `tret[0]`, `res = -2` is the circumpolar/exhausted
/// verdict with `tret[0] = 0.0`.
///
/// LEB-channeled bodies delegate to the shared [`rise_set_ut_serving`]
/// (identical acquisition: `calc_ut(FLG_TOPOCTR | FLG_SPEED)` words
/// rotated by `azalt` ECL2HOR, `calc_ut(FLG_EQUATORIAL | FLG_SPEED)`
/// transit words, point-source disc, Hindu projection arms). The
/// channel-less remainder mirrors it arm for arm over the geocentric
/// fallback words. The
/// twilight override is inert on both arms.
pub fn minor_rise_set_ut_serving<P: crate::positions::StateProvider>(
    serving: &SealedMinorEventServing<'_, P>,
    req: &RiseSetRequest,
) -> Result<(i32, [f64; 10])> {
    if serving.is_leb_channeled() {
        return rise_set_ut_serving(&serving.inner_transit()?, req);
    }
    kepler_rise_set_ut_serving(serving, req)
}

/// Next fallback-path rise/set over the geocentric Keplerian words:
/// mirrors [`rise_set_ut_serving`] arm for arm (geopos unpack,
/// event-kind gate, upfront placeability probe, pressure estimate,
/// `-100` dip horizon, disc selection with the zero-radius
/// point-source term, Sinclair-via-transit-core event height,
/// `calculate_rise_set`) with the topocentric composition replaced by
/// the geocentric fallback words, exactly like the release
/// `_calc_keplerian_fallback` path.
pub(super) fn kepler_rise_set_ut_serving<P: crate::positions::StateProvider>(
    serving: &SealedMinorEventServing<'_, P>,
    req: &RiseSetRequest,
) -> Result<(i32, [f64; 10])> {
    use crate::constants::{BIT_FIXED_DISC_SIZE, BIT_NO_REFRACTION, MOON, SUN};
    let rising = match resolve_event_kind(req.rsmi) {
        RiseSetEvent::Rise => true,
        RiseSetEvent::Set => false,
        _ => {
            return Err(Error::InvalidSearchStep {
                value: req.rsmi as f64,
            });
        }
    };
    validate_geopos(req.lon_deg_east, req.lat_deg, req.height_m)?;
    if !req.jd_start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: req.jd_start_ut });
    }
    // Upfront placeability probe: a coverage miss keeps its typed
    // range error (the escalated Earth/Sun miss naming the requested
    // body window); a readable word means the search runs on served
    // words.
    let jd_probe = serving.jd_tt_for_ut(req.jd_start_ut)?;
    serving.kepler_ecliptic_at_tt(jd_probe)?;
    let pressure = if req.atpress_mbar == 0.0 {
        pressure_from_altitude_m(req.height_m)
    } else {
        req.atpress_mbar
    };
    let horizon_alt = resolve_horizon_altitude_deg(
        req.horhgt_deg,
        crate::refraction::calc_dip(
            req.height_m.max(0.0),
            RISE_SET_DIP_LAPSE_RATE_K_M,
            pressure,
            req.attemp_c,
        ),
    );
    let rsmi_eff = req.rsmi;
    let use_refraction = rsmi_eff & BIT_NO_REFRACTION == 0;
    let disc = select_disc_point(rsmi_eff, false);
    let body_radius_au = body_disc_radius_au(serving.minor_id());
    debug_assert_eq!(body_radius_au, 0.0);
    let mut lift_max = 0.0;
    if use_refraction {
        lift_max += 0.75;
    }
    if disc != DiscPoint::Center {
        lift_max += 0.30;
    }
    let hindu = rsmi_eff & crate::constants::BIT_GEOCTR_NO_ECL_LAT != 0;
    // True altitude plus distance from the geocentric fallback place
    // rotated to the observer's horizon.
    let alt_dist_at_ut = |jd_ut: f64| -> Result<(f64, f64)> {
        let jd_tt = serving.jd_tt_for_ut(jd_ut)?;
        let (lon, lat, dist) = serving.kepler_ecliptic_at_tt(jd_tt)?;
        let (_, deps_rad) = (serving.serving.nutation_provider())(jd_tt)?;
        let eps_true_deg =
            crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
        let (ecl_lon, ecl_lat) = if hindu { (lon, 0.0) } else { (lon, lat) };
        let (ra_deg, dec_deg, dist) = cotrans(ecl_lon, ecl_lat, dist, -eps_true_deg);
        let last_hours = serving.last_hours_at_ut(jd_ut, req.lon_deg_east)?;
        let ha_deg = (last_hours * 15.0 - ra_deg).rem_euclid(360.0);
        let ha_rad = ha_deg.to_radians();
        let dec_rad = dec_deg.to_radians();
        let lat_rad = req.lat_deg.to_radians();
        let sin_alt = lat_rad.sin() * dec_rad.sin() + lat_rad.cos() * dec_rad.cos() * ha_rad.cos();
        Ok((sin_alt.clamp(-1.0, 1.0).asin().to_degrees(), dist))
    };
    // Height of the configured disc point above the configured horizon
    // .
    let dip_deg = crate::refraction::calc_dip(
        req.height_m.max(0.0),
        RISE_SET_DIP_LAPSE_RATE_K_M,
        pressure,
        req.attemp_c,
    );
    let mut event_height = |jd: f64| -> Result<f64> {
        let (alt_true, dist) = alt_dist_at_ut(jd)?;
        let mut cur = dist;
        if rsmi_eff & BIT_FIXED_DISC_SIZE != 0 {
            if serving.minor_id() == SUN {
                cur = FIXED_DISC_SUN_DIST_AU;
            } else if serving.minor_id() == MOON {
                cur = FIXED_DISC_MOON_DIST_AU;
            }
        }
        let mut alt_pt = alt_true + disc_altitude_term_deg(disc, body_radius_au, cur)?;
        if use_refraction {
            let mut refraction_deg = |app: f64| -> Result<f64> {
                Ok(crate::refraction::sinclair_refraction_deg(
                    app,
                    pressure,
                    req.attemp_c,
                ))
            };
            alt_pt = true_to_apparent_altitude(alt_pt, dip_deg, &mut refraction_deg, jd)?;
        }
        Ok(alt_pt - horizon_alt)
    };
    // The circumpolar pre-check sees the same declination the altitude
    // solve uses: the fallback apparent declination on the ordinary
    // path, or the projected declination of the latitude-zeroed place
    // under the Hindu bit.
    let mut get_ra_dec = |jd: f64| -> Result<(f64, f64)> {
        if hindu {
            let jd_tt = serving.jd_tt_for_ut(jd)?;
            let (lon, _, dist) = serving.kepler_ecliptic_at_tt(jd_tt)?;
            let (_, deps_rad) = (serving.serving.nutation_provider())(jd_tt)?;
            let eps_true_deg =
                crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
            let (ra_deg, dec_deg, _) = cotrans(lon, 0.0, dist, -eps_true_deg);
            Ok((ra_deg / 15.0, dec_deg))
        } else {
            serving.kepler_ra_dec_at_ut(jd)
        }
    };
    match calculate_rise_set(
        req.jd_start_ut,
        req.lat_deg,
        rising,
        &mut event_height,
        &mut get_ra_dec,
        horizon_alt,
        lift_max,
    )? {
        RiseSetOutcome::Found { jd_ut } => Ok((0, make_tret(jd_ut))),
        RiseSetOutcome::Circumpolar => Ok((-2, make_tret(0.0))),
    }
}

/// Next fictitious-body meridian transit JD (UT) over the sealed serving:
/// `lower_transit` `false` is the upper culmination (`CALC_MTRANSIT`),
/// `true` the anti-culmination (`CALC_ITRANSIT`).
pub fn fictitious_transit_ut_serving<P: crate::positions::StateProvider>(
    serving: &SealedFictitiousEventServing<'_, P>,
    jd_start: f64,
    lon_deg_east: f64,
    lower_transit: bool,
) -> Result<f64> {
    Ok(serving
        .fict_transit_at_ut(jd_start, lon_deg_east, lower_transit)?
        .jd_ut)
}

/// Next fictitious-body rise/set over the sealed serving: `(res, tret)`
/// with the release `rise_trans_true_hor` contract — `res = 0` carries
/// the event JD (UT) in `tret[0]`, `res = -2` is the circumpolar/exhausted
/// verdict with `tret[0] = 0.0`.
///
/// Mirrors [`minor_rise_set_ut_serving`] arm for arm (geopos unpack,
/// event-kind gate, upfront placeability probe, pressure estimate,
/// `-100` dip horizon, disc selection with the zero-radius point-source
/// term, Sinclair-via-transit-core event height,
/// `calculate_rise_set`) over the topocentric apparent words — the
/// release fictitious path honors `FLG_TOPOCTR` through
/// `_apply_hypothetical_topocentric`, unlike the 4b Keplerian fallback.
/// The Hindu arm uses the geocentric words with the ecliptic latitude
/// zeroed. The twilight override
/// is inert.
pub fn fictitious_rise_set_ut_serving<P: crate::positions::StateProvider>(
    serving: &SealedFictitiousEventServing<'_, P>,
    req: &RiseSetRequest,
) -> Result<(i32, [f64; 10])> {
    fict_rise_set_ut_serving(serving, req)
}

/// Next rise/set over the fictitious apparent words: mirrors
/// [`kepler_rise_set_ut_serving`] with the topocentric composition in
/// place of the geocentric fallback words.
pub(super) fn fict_rise_set_ut_serving<P: crate::positions::StateProvider>(
    serving: &SealedFictitiousEventServing<'_, P>,
    req: &RiseSetRequest,
) -> Result<(i32, [f64; 10])> {
    use crate::constants::{BIT_FIXED_DISC_SIZE, BIT_NO_REFRACTION, MOON, SUN};
    let rising = match resolve_event_kind(req.rsmi) {
        RiseSetEvent::Rise => true,
        RiseSetEvent::Set => false,
        _ => {
            return Err(Error::InvalidSearchStep {
                value: req.rsmi as f64,
            });
        }
    };
    validate_geopos(req.lon_deg_east, req.lat_deg, req.height_m)?;
    if !req.jd_start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: req.jd_start_ut });
    }
    // Upfront placeability probe: a coverage miss keeps its typed
    // range error; a readable word means the search runs on served
    // words. Model-less ids never reach here (refused at construction
    // with the probe's illegal-body mapping).
    let jd_probe = serving.jd_tt_for_ut(req.jd_start_ut)?;
    serving.fict_geo_ecliptic_at_tt(jd_probe)?;
    let pressure = if req.atpress_mbar == 0.0 {
        pressure_from_altitude_m(req.height_m)
    } else {
        req.atpress_mbar
    };
    let horizon_alt = resolve_horizon_altitude_deg(
        req.horhgt_deg,
        crate::refraction::calc_dip(
            req.height_m.max(0.0),
            RISE_SET_DIP_LAPSE_RATE_K_M,
            pressure,
            req.attemp_c,
        ),
    );
    let rsmi_eff = req.rsmi;
    let use_refraction = rsmi_eff & BIT_NO_REFRACTION == 0;
    let disc = select_disc_point(rsmi_eff, false);
    let body_radius_au = body_disc_radius_au(serving.body_id());
    debug_assert_eq!(body_radius_au, 0.0);
    let mut lift_max = 0.0;
    if use_refraction {
        lift_max += 0.75;
    }
    if disc != DiscPoint::Center {
        lift_max += 0.30;
    }
    let hindu = rsmi_eff & crate::constants::BIT_GEOCTR_NO_ECL_LAT != 0;
    // True altitude plus distance from the topocentric apparent place
    // rotated to the observer's horizon.
    let site = (req.lon_deg_east, req.lat_deg, req.height_m);
    let alt_dist_at_ut = |jd_ut: f64| -> Result<(f64, f64)> {
        let jd_tt = serving.jd_tt_for_ut(jd_ut)?;
        let (lon, lat, dist) = if hindu {
            let (lon, _, dist) = serving.fict_geo_ecliptic_at_tt(jd_tt)?;
            (lon, 0.0, dist)
        } else {
            serving.fict_topo_ecliptic_at_ut(jd_ut, site.0, site.1, site.2)?
        };
        let (_, deps_rad) = (serving.serving.nutation_provider())(jd_tt)?;
        let eps_true_deg =
            crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
        let (ra_deg, dec_deg, dist) = cotrans(lon, lat, dist, -eps_true_deg);
        let last_hours = serving.last_hours_at_ut(jd_ut, req.lon_deg_east)?;
        let ha_deg = (last_hours * 15.0 - ra_deg).rem_euclid(360.0);
        let ha_rad = ha_deg.to_radians();
        let dec_rad = dec_deg.to_radians();
        let lat_rad = req.lat_deg.to_radians();
        let sin_alt = lat_rad.sin() * dec_rad.sin() + lat_rad.cos() * dec_rad.cos() * ha_rad.cos();
        Ok((sin_alt.clamp(-1.0, 1.0).asin().to_degrees(), dist))
    };
    // Height of the configured disc point above the configured horizon
    // .
    let dip_deg = crate::refraction::calc_dip(
        req.height_m.max(0.0),
        RISE_SET_DIP_LAPSE_RATE_K_M,
        pressure,
        req.attemp_c,
    );
    let mut event_height = |jd: f64| -> Result<f64> {
        let (alt_true, dist) = alt_dist_at_ut(jd)?;
        let mut cur = dist;
        if rsmi_eff & BIT_FIXED_DISC_SIZE != 0 {
            if serving.body_id() == SUN {
                cur = FIXED_DISC_SUN_DIST_AU;
            } else if serving.body_id() == MOON {
                cur = FIXED_DISC_MOON_DIST_AU;
            }
        }
        let mut alt_pt = alt_true + disc_altitude_term_deg(disc, body_radius_au, cur)?;
        if use_refraction {
            let mut refraction_deg = |app: f64| -> Result<f64> {
                Ok(crate::refraction::sinclair_refraction_deg(
                    app,
                    pressure,
                    req.attemp_c,
                ))
            };
            alt_pt = true_to_apparent_altitude(alt_pt, dip_deg, &mut refraction_deg, jd)?;
        }
        Ok(alt_pt - horizon_alt)
    };
    // The circumpolar pre-check sees the same declination the altitude
    // solve uses: the topocentric apparent declination on the ordinary
    // path, or the projected declination of the latitude-zeroed
    // geocentric place under the Hindu bit.
    let mut get_ra_dec = |jd: f64| -> Result<(f64, f64)> {
        if hindu {
            let jd_tt = serving.jd_tt_for_ut(jd)?;
            let (lon, _, dist) = serving.fict_geo_ecliptic_at_tt(jd_tt)?;
            let (_, deps_rad) = (serving.serving.nutation_provider())(jd_tt)?;
            let eps_true_deg =
                crate::positions::vondrak_mean_obliquity_deg(jd_tt)? + deps_rad.to_degrees();
            let (ra_deg, dec_deg, _) = cotrans(lon, 0.0, dist, -eps_true_deg);
            Ok((ra_deg / 15.0, dec_deg))
        } else {
            serving.fict_ra_dec_at_ut(jd)
        }
    };
    match calculate_rise_set(
        req.jd_start_ut,
        req.lat_deg,
        rising,
        &mut event_height,
        &mut get_ra_dec,
        horizon_alt,
        lift_max,
    )? {
        RiseSetOutcome::Found { jd_ut } => Ok((0, make_tret(jd_ut))),
        RiseSetOutcome::Circumpolar => Ok((-2, make_tret(0.0))),
    }
}

#[cfg(test)]
mod acquired_horizon_frame_tests {
    use super::*;

    #[test]
    fn completed_radians_preserve_unknown_geometry_and_distance() {
        let (azimuth, altitude) = horizon_forward_geometry_radians(f64::NAN, 0.0, 0.0);
        assert!(azimuth.is_nan());
        assert!(altitude.is_nan());
        for (azimuth, altitude, latitude) in [
            (f64::NAN, 0.0, 0.0),
            (0.0, f64::NAN, 0.0),
            (0.0, 0.0, f64::NAN),
        ] {
            let (hour_angle, declination) =
                horizon_reverse_geometry_radians(azimuth, altitude, latitude);
            assert!(hour_angle.is_nan());
            assert!(declination.is_nan());
        }
        let distance = f64::from_bits(0x7ff8_0000_0000_1234);
        let (longitude, latitude, retained) =
            super::super::cotrans_radians(f64::NAN, 0.0, 0.0, distance);
        assert!(longitude.is_nan());
        assert!(latitude.is_nan());
        assert_eq!(retained.to_bits(), distance.to_bits());
    }

    #[test]
    fn greenwich_completed_frame_precedes_local_normalization() {
        let dt = crate::time::SiderealDeltaT::try_new(0.0, 0.0, 0.0).unwrap();
        let frame = HorizonFrame::from_degrees(0.0, 10_800.0);
        let greenwich = frame.greenwich_apparent_hours(2451545.0, &dt).unwrap();
        assert!(greenwich > 24.0);
        let local = crate::time::local_apparent_sidereal_time_hours_from_angles(
            2451545.0, 0.0, 0.0, 10_800.0, &dt,
        )
        .unwrap();
        assert_eq!((greenwich % 24.0).to_bits(), local.to_bits());
        assert!(
            HorizonFrame::from_degrees(f64::NAN, f64::NAN)
                .greenwich_apparent_hours(2451545.0, &dt)
                .unwrap()
                .is_nan()
        );
        assert!(matches!(
            frame.greenwich_apparent_hours(f64::NAN, &dt),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn prepared_default_frame_preserves_both_checked_projection_words() {
        for jd in [2415020.5, 2451545.0, 2460000.5] {
            let frame =
                HorizonFrame::at_tt(crate::time::default_timescale_tt_from_ut1(jd).unwrap())
                    .unwrap();
            let frame = HorizonFrame::from_degrees(frame.eps_true_deg, frame.dpsi_deg);
            for flag in [0, 1] {
                let legacy = azalt_deg(
                    jd,
                    flag,
                    12.0,
                    45.0,
                    20.0,
                    1013.25,
                    15.0,
                    (90.0, 23.5, 1.0),
                    &|_| Ok(0.25),
                    0.0,
                    0.0,
                    0.0065,
                )
                .unwrap();
                let supplied = azalt_deg_with_frame(
                    jd,
                    flag,
                    12.0,
                    45.0,
                    20.0,
                    1013.25,
                    15.0,
                    (90.0, 23.5, 1.0),
                    &|_| Ok(0.25),
                    0.0,
                    0.0,
                    0.0065,
                    frame,
                )
                .unwrap();
                assert_eq!(
                    [legacy.0.to_bits(), legacy.1.to_bits(), legacy.2.to_bits()],
                    [
                        supplied.0.to_bits(),
                        supplied.1.to_bits(),
                        supplied.2.to_bits()
                    ]
                );
                let legacy = azalt_rev_deg(
                    jd,
                    flag,
                    12.0,
                    45.0,
                    0.0,
                    90.0,
                    23.5,
                    &|_| Ok(0.25),
                    0.0,
                    0.0,
                )
                .unwrap();
                let supplied = azalt_rev_deg_with_frame(
                    jd,
                    flag,
                    12.0,
                    45.0,
                    0.0,
                    90.0,
                    23.5,
                    &|_| Ok(0.25),
                    0.0,
                    0.0,
                    frame,
                )
                .unwrap();
                assert_eq!(
                    [legacy.0.to_bits(), legacy.1.to_bits()],
                    [supplied.0.to_bits(), supplied.1.to_bits()]
                );
            }
        }
    }

    #[test]
    fn supplied_ieee_frame_preserves_checked_clock_domains() {
        use std::cell::Cell;
        let frame = HorizonFrame::from_degrees(f64::NAN, f64::NAN);
        let calls = Cell::new(0);
        let delta = |_| {
            calls.set(calls.get() + 1);
            Ok(0.0)
        };
        let (ra, dec) = azalt_rev_deg_with_frame(
            2451545.0, HOR2EQU, 0.0, 0.0, 0.0, 90.0, 45.0, &delta, 0.0, 0.0, frame,
        )
        .unwrap();
        assert!(ra.is_nan());
        assert!(dec.is_finite());
        assert_eq!(calls.get(), 1);
        let (longitude, latitude) = azalt_rev_deg_with_frame(
            2451545.0, 0, 0.0, 0.0, 0.0, 90.0, 45.0, &delta, 0.0, 0.0, frame,
        )
        .unwrap();
        assert!(longitude.is_nan());
        assert!(latitude.is_nan());
        let (azimuth, true_altitude, apparent_altitude) = azalt_deg_with_frame(
            2451545.0,
            1,
            0.0,
            0.0,
            0.0,
            1013.25,
            15.0,
            (0.0, 0.0, 1.0),
            &delta,
            0.0,
            0.0,
            0.0065,
            frame,
        )
        .unwrap();
        assert!(azimuth.is_nan());
        assert!(true_altitude.is_nan());
        assert!(apparent_altitude.is_nan());
        assert_eq!(calls.get(), 3);
        assert!(matches!(
            HorizonFrame::at_tt(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        let dt = crate::time::SiderealDeltaT::try_new(0.0, 0.0, 0.0).unwrap();
        assert!(matches!(
            crate::time::local_apparent_sidereal_time_hours(
                2451545.0,
                0.0,
                f64::NAN,
                f64::NAN,
                &dt,
            ),
            Err(Error::NonFiniteArgument { .. })
        ));
        for jd in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                azalt_rev_deg_with_frame(
                    jd, HOR2EQU, 0.0, 0.0, 0.0, 90.0, 45.0, &delta, 0.0, 0.0, frame,
                ),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
        assert_eq!(calls.get(), 3);
    }
}

#[cfg(test)]
mod smooth_horizon_refinement_tests {
    use super::*;

    #[test]
    fn nonlinear_shallow_horizons_use_the_actual_crossing_in_both_directions() {
        let epoch: f64 = 2460000.0;
        let root = epoch + 0.03141592653589793;
        for rising in [false, true] {
            let direction = if rising { 1.0 } else { -1.0 };
            let mut height = |t| {
                let dt = t - root;
                Ok(direction * (0.0001 * dt + 0.002 * dt * dt * dt))
            };
            let actual =
                refine_horizon_crossing(&mut height, epoch, epoch + 0.1, rising, false).unwrap();
            assert_eq!(actual, root);
        }
    }

    #[test]
    fn an_adjacent_horizon_bracket_keeps_the_endpoint_with_less_residual() {
        let lo: f64 = 2460000.0;
        let hi = lo.next_up();
        let mut height = |t| Ok((t - lo) - 0.75 * (hi - lo));
        assert_eq!(
            refine_horizon_crossing(&mut height, lo, hi, true, false).unwrap(),
            hi
        );
    }

    #[test]
    fn actual_height_errors_propagate_from_refinement() {
        let err = Error::NonFiniteArgument {
            value: f64::INFINITY,
        };
        let mut height = |_| Err(err.clone());
        assert_eq!(
            refine_horizon_crossing(&mut height, 0.0, 1.0, true, false).unwrap_err(),
            err
        );
    }

    #[test]
    fn refraction_visibility_steps_keep_the_boundary_in_both_directions() {
        let epoch: f64 = 2460000.0;
        let root = epoch + 0.03141592653589793;
        for rising in [false, true] {
            let mut height = |t| {
                let visible = if rising { t >= root } else { t <= root };
                Ok(if visible { 0.0 } else { -0.5 })
            };
            let actual =
                refine_horizon_crossing(&mut height, epoch, epoch + 0.1, rising, true).unwrap();
            assert_eq!(actual, root);
        }
    }
}
