//! Fixed-star horizon events over the public scalar-UT catalog place.
use super::rise::{
    EventParameters, default_gast_hours, native_star_event_callbacks, projected_ra_dec,
};
use super::*;
use crate::events::{self, RiseSetRequest};

impl Session {
    /// Next fixed-star rise, set or meridian transit from an original UT1 JD.
    /// This is [`Self::star_rise_trans_true_hor`] with horizon zero. The key uses
    /// V1 catalog lookup through [`Self::fixstar_ut`]; no lookup/data acquisition
    /// occurs until the selected event consumes its first public star place.
    /// Geographic coordinates are east longitude/latitude degrees and metres;
    /// pressure is mbar (zero estimates it), temperature is Celsius.
    /// Returns `(0, tret)` on success or `(-2, [0;10])` on no crossing. Only
    /// `tret[0]` carries a UT1 JD; the other nine words are positive zero.
    #[allow(clippy::too_many_arguments)]
    pub fn star_rise_trans(
        &self,
        key: &str,
        jd_ut: f64,
        rsmi: i32,
        geopos: [f64; 3],
        pressure: f64,
        temperature: f64,
        flags: i32,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.star_rise_trans_true_hor(key, jd_ut, rsmi, geopos, pressure, temperature, 0.0, flags)
    }

    /// Fixed-star event with lazy binding main-image configuration.
    /// Inputs, units, consumed-stage order and output follow [`Self::star_rise_trans`].
    /// The resolver runs only at the actual public scalar-UT star acquisition;
    /// setup/catalog/clock errors retain their precedence. One owned success is
    /// retained for all event samples, without mutating the raw caller context.
    #[allow(clippy::too_many_arguments)]
    pub fn star_rise_trans_with_position_context(
        &self,
        key: &str,
        jd_ut: f64,
        rsmi: i32,
        geopos: [f64; 3],
        pressure: f64,
        temperature: f64,
        flags: i32,
        resolve: &crate::positions::star_serving::StarPositionContextResolver<'_>,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.star_rise_trans_true_hor_with_position_context(
            key,
            jd_ut,
            rsmi,
            geopos,
            pressure,
            temperature,
            0.0,
            flags,
            resolve,
        )
    }

    /// Fixed-star horizon event with horizon degrees at original UT1.
    ///
    /// Event precedence is SET > RISE > ITRANSIT > MTRANSIT, defaulting to rise.
    /// Stars use their geometric centre regardless of disc-size/bottom bits;
    /// twilight selects -6/-12/-18 degrees and disables refraction. Horizon
    /// -100 requests the refracted sea-horizon dip plus 1e-4 degree. Pressure
    /// zero uses `1013.25*max(1-0.0065*height/288,0)^5.255`; height is metres.
    /// Refraction uses the shared Sinclair inverse and dip-floor branch.
    ///
    /// Every position sample is a genuine [`Self::fixstar_ut`] call: ordinary
    /// altitude and Hindu latitude-zero projection request SPEED, while transit
    /// and ordinary circumpolar precheck request EQUATORIAL|SPEED. Caller `flags`
    /// is not forwarded to these fixed acquisitions; its one-try/XYZ/RADIANS
    /// bits do not change the event geometry. The configured provider/image
    /// belongs to this Session. Explicit LEB still refuses JPL substitution;
    /// NativeJpl uses the strict catalog/state/range path. No temporary
    /// observer is installed or validated for a geocentric star calculation.
    ///
    /// [`Self::azalt`] supplies the spherical ECL2HOR true-altitude projection
    /// with its configured model sidereal clock and independent frame TT. Star
    /// UT1-to-scalar-TT/default observer clocks stay in `fixstar_ut`; meridian
    /// GAST and Hindu obliquity use their separately bundled original-UT1 clocks.
    /// A vector-altaz substitute or model-time replacement is not performed.
    ///
    /// Rise/set searches a two-day sample window and shares
    /// the bounded grazing-aware scan, 60 bisections, slope-aware 0.03-second
    /// exit and conditional in-bracket Newton correction. Transit shares the
    /// thirty-iteration hour-angle solver. These limits are algorithmic, not
    /// certified physical accuracy. There is no eager all-field finite check:
    /// pressure/horizon setup precedes the first catalog/place acquisition,
    /// and unused transit latitude is not consumed. Actual catalog, clock/frame/provider/support and
    /// sample-range errors propagate in acquisition order. Existing contexts,
    /// registry allocations and active snapshots are retained unchanged. All
    /// callbacks/key borrows are synchronous and end before this method returns.
    #[allow(clippy::too_many_arguments)]
    pub fn star_rise_trans_true_hor(
        &self,
        key: &str,
        jd_ut: f64,
        rsmi: i32,
        geopos: [f64; 3],
        pressure: f64,
        temperature: f64,
        horhgt: f64,
        flags: i32,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.star_rise_trans_true_hor_with_main_context(
            key,
            jd_ut,
            rsmi,
            geopos,
            pressure,
            temperature,
            horhgt,
            flags,
            None,
        )
    }

    /// Horizon event with a borrowed lazy main-configuration resolver.
    /// Semantics follow [`Self::star_rise_trans_true_hor`]; only the consumed
    /// stellar provider selects binding settings. The callback is synchronous,
    /// runs unlocked, and cannot replace nonempty supplied registry images.
    #[allow(clippy::too_many_arguments)]
    pub fn star_rise_trans_true_hor_with_position_context(
        &self,
        key: &str,
        jd_ut: f64,
        rsmi: i32,
        geopos: [f64; 3],
        pressure: f64,
        temperature: f64,
        horhgt: f64,
        flags: i32,
        resolve: &crate::positions::star_serving::StarPositionContextResolver<'_>,
    ) -> Result<(i32, [f64; 10]), Error> {
        self.star_rise_trans_true_hor_with_main_context(
            key,
            jd_ut,
            rsmi,
            geopos,
            pressure,
            temperature,
            horhgt,
            flags,
            Some(resolve),
        )
    }

    /// Shared event body; setup and every numerical expression remain in the
    /// same order. Optional configuration belongs to the scalar star callback,
    /// not the horizon/model/GAST clocks or eager event-input validation.
    #[allow(clippy::too_many_arguments)]
    fn star_rise_trans_true_hor_with_main_context(
        &self,
        key: &str,
        jd_ut: f64,
        rsmi: i32,
        geopos: [f64; 3],
        pressure: f64,
        temperature: f64,
        horhgt: f64,
        flags: i32,
        resolve: Option<&crate::positions::star_serving::StarPositionContextResolver<'_>>,
    ) -> Result<(i32, [f64; 10]), Error> {
        let req = RiseSetRequest {
            jd_start_ut: jd_ut,
            lon_deg_east: geopos[0],
            lat_deg: geopos[1],
            height_m: geopos[2],
            atpress_mbar: pressure,
            attemp_c: temperature,
            horhgt_deg: horhgt,
            rsmi,
        };
        let _ = flags;
        let parameters = star_event_parameters(&req)?;
        let calc = |ut, internal_flags| {
            star_position_callback(ut, internal_flags, &|ut, actual_flags| {
                match resolve {
                    Some(resolve) => self.fixstar_ut_captured_with_position_context(
                        key,
                        ut,
                        actual_flags,
                        resolve,
                    ),
                    None => self.fixstar_ut_captured(key, ut, actual_flags),
                }
                .into_catalog_result()
                .map(|star| star.values)
            })
        };
        let altitude = |ut, place: [f64; 3]| {
            self.azalt(
                ut,
                crate::events::ECL2HOR,
                geopos[0],
                geopos[1],
                geopos[2],
                0.0,
                0.0,
                place[0],
                place[1],
            )
            .map(|horizontal| horizontal.true_altitude)
        };
        native_star_event_callbacks(
            -1,
            &req,
            &parameters,
            &calc,
            &altitude,
            &projected_ra_dec,
            &default_gast_hours,
        )
    }
}

/// Map only the shared kernel's internal ordinary-altitude TOPOCTR request to
/// the fixed star's geocentric SPEED acquisition. EQUATORIAL|SPEED is retained
/// for transit/precheck. No caller flags, observer mutation, position values or
/// epoch conversion are introduced; the actual public call owns those roles.
/// `ut` remains scalar UT1 JD; output is the public six spherical components in
/// degrees/AU and daily rates. Lookup and acquisition errors pass through.
fn star_position_callback(
    ut: f64,
    internal_flags: i32,
    public_fixstar: &dyn Fn(f64, i32) -> Result<[f64; 6], Error>,
) -> Result<[f64; 6], Error> {
    public_fixstar(ut, internal_flags & !FLG_TOPOCTR)
}

/// Resolve the star's centre/refraction/horizon choices without star/data I/O.
/// Pressure is mbar, temperature Celsius and horizon/dip degrees. The -100 dip
/// calculation precedes twilight replacement.
/// Twilight priority is civil > nautical > astronomical. Fixed stars have no
/// disc-radius lift: only refraction contributes the 0.75-degree precheck
/// allowance. This preliminary stage performs no observer or clock acquisition.
pub(super) fn star_event_parameters(req: &RiseSetRequest) -> Result<EventParameters, Error> {
    let pressure = if req.atpress_mbar == 0.0 {
        star_pressure_from_altitude_m(req.height_m)?
    } else {
        req.atpress_mbar
    };
    let mut horizon = req.horhgt_deg;
    if horizon == -100.0 {
        let dip = crate::refraction::calc_dip(
            star_nonnegative(req.height_m),
            events::RISE_SET_DIP_LAPSE_RATE_K_M,
            pressure,
            req.attemp_c,
        );
        horizon = events::resolve_horizon_altitude_deg(horizon, dip);
    }
    let mut effective = req.rsmi;
    if let Some(twilight) = events::twilight_override(req.rsmi, true) {
        effective |= BIT_NO_REFRACTION | BIT_DISC_CENTER;
        horizon = twilight.horizon_alt_deg;
    }
    Ok(EventParameters {
        pressure,
        horizon,
        effective,
        disc: events::DiscPoint::Center,
        radius: 0.0,
        lift: if effective & BIT_NO_REFRACTION == 0 {
            0.75
        } else {
            0.0
        },
    })
}

/// Barometric pressure `1013.25 * max(1 - 0.0065*h/288, 0)^5.255`, in mbar.
/// `h` is height in metres. The ordered first-argument maximum retains NaN;
/// an infinite base can produce infinity without a power-overflow refusal.
/// Only a finite base whose positive fractional power overflows is a consumed
/// primitive failure. Subsequent floating multiplication is a separate stage
/// and can produce infinity without raising this error. Owned operand words
/// let an attached binding consume the genuine built-in power constructor;
/// no errno, exception arguments or backend graph is inferred from Display.
fn star_pressure_from_altitude_m(height_m: f64) -> Result<f64, Error> {
    let base = 1.0 - 0.0065 * height_m / 288.0;
    let clamped = star_nonnegative(base);
    let exponent = 5.255_f64;
    let powered = clamped.powf(exponent);
    if clamped.is_finite() && powered.is_infinite() {
        return Err(Error::NumericalOverflow {
            operation: "barometric pressure",
            message: "barometric pressure overflow".into(),
        });
    }
    Ok(1013.25 * powered)
}

/// First-argument maximum with zero. Ordered comparison leaves NaN and a tied
/// negative zero unchanged; negative finite/infinite values select positive
/// zero. This performs no domain admission, callback or source acquisition.
pub(super) fn star_nonnegative(value: f64) -> f64 {
    if value < 0.0 { 0.0 } else { value }
}

/// Meridian transit with source-consumed stage ordering. RA hours is acquired
/// before bundled GAST hours, then east longitude degrees forms the hour angle.
/// Latitude, height and atmospheric fields are unused by this solve. The mean
/// sidereal-day seed, lower-transit wrap and thirty iterations retain the shared
/// checked solver's arithmetic; there is no eager finite-field admission here.
/// With GAST `G` and RA `a` in hours and longitude `l` in degrees,
/// `H = mod360((G+l/15)*15-a*15)` is folded once to (-180,180]. Lower transit
/// shifts H by 180 degrees. The seed advances by the forward arc divided by
/// 360 times S=0.99726957 days; refinement subtracts H/(360/S) days, ending at
/// |H|<1/3600 degree or after thirty steps.
/// Position/clock callbacks own their real failures; cap exhaustion returns the
/// last UT1 estimate with status zero, without a physical-accuracy assertion.
/// Callback borrows are synchronous and no source/native resource is replaced.
pub(super) fn star_calculate_transit(
    jd_start: f64,
    lon_deg_east: f64,
    lower_transit: bool,
    get_ra_dec: &dyn Fn(f64) -> Result<(f64, f64), Error>,
    gast_hours: &dyn Fn(f64) -> Result<f64, Error>,
) -> Result<(i32, [f64; 10]), Error> {
    let hour_angle = |jd: f64| -> Result<f64, Error> {
        let (ra_hours, _) = get_ra_dec(jd)?;
        let gast = gast_hours(jd)?;
        let mut ha =
            events::normalize_target_deg((gast + lon_deg_east / 15.0) * 15.0 - ra_hours * 15.0);
        if ha > 180.0 {
            ha -= 360.0;
        }
        if lower_transit {
            ha = events::normalize_target_deg(ha - 180.0);
            if ha > 180.0 {
                ha -= 360.0;
            }
        }
        Ok(ha)
    };
    let ha = hour_angle(jd_start)?;
    let mut jd_guess = if ha > 0.0 {
        jd_start + (360.0 - ha) / 360.0 * events::SIDEREAL_DAY
    } else {
        jd_start + (-ha) / 360.0 * events::SIDEREAL_DAY
    };
    for _ in 0..events::TRANSIT_MAX_ITER {
        let ha = hour_angle(jd_guess)?;
        if ha.abs() < events::TRANSIT_TOL_HA_DEG {
            return Ok((0, events::make_tret(jd_guess)));
        }
        jd_guess -= ha / (360.0 / events::SIDEREAL_DAY);
    }
    Ok((0, events::make_tret(jd_guess)))
}

#[cfg(test)]
#[path = "star_rise_tests.rs"]
mod tests;
