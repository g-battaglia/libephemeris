//! Compact eclipse-module visibility models with direct native apparent geometry.
//!
//! The fixed-instant model returns eight details and uses a piecewise solar
//! background plus a simple lunar background. The event model searches a
//! 400-day sampled visibility window with its own twilight magnitude relation.
//! These are distinct project convenience models, not projections of the
//! ten-detail VISLIMIT retinal-contrast model. Apparent positions, photometry,
//! stellar catalog resolution and horizon transforms use the shared engine.
use super::*;
use crate::positions::vislim::NativeVisibilityObserver;
use std::fmt;

/// A disabled native-main capability or a calculation/data failure.
#[derive(Debug)]
pub enum EclipseVisibilityError {
    /// Forced local-LEB mode does not permit this direct-JPL model.
    NativeMainDisabled,
    /// An ordinary native input, image, clock or calculation error.
    Native(Error),
}
impl fmt::Display for EclipseVisibilityError {
    fn fmt(&self, output: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NativeMainDisabled => output.write_str(
                "this visibility model requires native JPL/SPK access, disabled in forced LEB mode",
            ),
            Self::Native(error) => error.fmt(output),
        }
    }
}
impl std::error::Error for EclipseVisibilityError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Native(error) => Some(error),
            Self::NativeMainDisabled => None,
        }
    }
}
impl From<Error> for EclipseVisibilityError {
    fn from(error: Error) -> Self {
        Self::Native(error)
    }
}
type VisibilityResult<T> = std::result::Result<T, EclipseVisibilityError>;

fn invalid(detail: impl Into<String>) -> Error {
    Error::InvalidEclipseSearch {
        entry: "eclipse_visibility",
        detail: detail.into(),
    }
}
fn planet(body: i32) -> bool {
    (SUN..=PLUTO).contains(&body) || body == EARTH
}
fn mathematical_failure(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::NonFiniteArgument { .. }
            | Error::MathDomain { .. }
            | Error::NumericalOverflow { .. }
            | Error::RationalDivisionByZero
    )
}
fn magnitude(session: &Session, jd: f64, body: i32, flags: i32) -> Result<f64, Error> {
    match session.pheno_ut(jd, body, flags) {
        Ok(result) => Ok(result.values[4]),
        Err(error) if mathematical_failure(&error) => Ok(0.0),
        Err(error) => Err(error),
    }
}

impl Session {
    fn eclipse_visibility_owner(&self) -> VisibilityResult<Session> {
        if self.ctx.mode() == ProviderMode::SealedLeb {
            return Err(EclipseVisibilityError::NativeMainDisabled);
        }
        let context = self.ctx.clone().with_mode(ProviderMode::NativeJpl);
        let owner = MainImageLifecycle::new().acquire(&context)?;
        let context = owner.native_provider_context()?;
        Ok(Session::from_context(context))
    }

    /// Fixed-instant compact visual model with eight diagnostic values.
    ///
    /// `jd` is UT1, geographic longitude/latitude degrees east/north, height
    /// metres, atmosphere pressure mbar/temperature Celsius/humidity percent/
    /// meteorological range km. Missing sequence values select 0-degree site,
    /// 1013.25 mbar/15 C/50%/0 km and age36/Snellen1. Objects are ordinary
    /// planet names/numeric IDs or native fixed-star catalog names. The native
    /// JPL main supplies apparent planet geometry at the bundled frame clock;
    /// photometry and catalog/horizon calls retain the configured clock.
    ///
    /// Returns vision status and `(limit,object_alt,object_az,Sun_alt,Sun_az,
    /// Moon_alt,Moon_az,apparent_object_mag)`. Planet azimuths start north;
    /// stellar azimuth uses the public south-westward horizon convention.
    /// A below-horizon object returns limit0 and its intrinsic magnitude.
    /// DARK/NOMOON replace the corresponding altitude by -90 degrees.
    /// Model limits: approximate solar/lunar backgrounds, constant V-band
    /// molecular/ozone extinction and linear age correction. Forced LEB is
    /// unavailable; data/provider and non-mathematical photometry failures
    /// propagate. It does not select the separate retinal-contrast model.
    #[allow(clippy::too_many_arguments)]
    pub fn eclipse_vis_limit_mag(
        &self,
        jd: f64,
        geopos: &[f64],
        atmo: &[f64],
        observer: &[f64],
        objname: &str,
        flags: i32,
    ) -> VisibilityResult<(f64, [f64; 8])> {
        if objname.is_empty() {
            return Err(invalid("objname cannot be empty").into());
        }
        let word = |values: &[f64], index: usize, default: f64| {
            values.get(index).copied().unwrap_or(default)
        };
        let site = [
            word(geopos, 0, 0.0),
            word(geopos, 1, 0.0),
            word(geopos, 2, 0.0),
        ];
        let pressure = word(atmo, 0, 1013.25);
        let temperature = word(atmo, 1, 15.0);
        let humidity = word(atmo, 2, 50.0);
        let met_range = word(atmo, 3, 0.0);
        let age = word(observer, 0, 36.0);
        let snellen = word(observer, 1, 1.0);
        let owner = self.eclipse_visibility_owner()?;
        // The same acquired main supplies public photometry and stellar calls.
        // Keep the caller's provider mode so Auto still uses its owned LEB
        // channels where available, with identical model/observer settings.
        let calculation = Session::from_context(owner.ctx.clone().with_mode(self.ctx.mode()));
        let sky = NativeVisibilityObserver::new(owner.ctx.clone(), jd, site)?;
        let sun = sky.observe(SUN)?;
        let moon = sky.observe(MOON)?;
        let numeric = objname.trim().parse::<i32>().ok();
        let named = match objname.trim().to_uppercase().as_str() {
            "SUN" => Some(SUN),
            "MOON" => Some(MOON),
            "MERCURY" => Some(MERCURY),
            "VENUS" => Some(VENUS),
            "MARS" => Some(MARS),
            "JUPITER" => Some(JUPITER),
            "SATURN" => Some(SATURN),
            "URANUS" => Some(URANUS),
            "NEPTUNE" => Some(NEPTUNE),
            "PLUTO" => Some(PLUTO),
            _ => None,
        };
        let (altitude, azimuth, object_magnitude) = if let Some(body) = numeric.or(named) {
            if !planet(body) {
                return Err(invalid(format!("illegal planet number {body}.")).into());
            }
            let place = sky.observe(body)?;
            (
                place.altitude_deg,
                place.azimuth_deg,
                magnitude(&calculation, jd, body, flags)?,
            )
        } else {
            let place = calculation.fixstar2_ut(objname, jd, flags & 255)?;
            let object_magnitude = crate::stars::fixstar2_mag(objname)?.0;
            let horizontal = self.azalt(
                jd,
                crate::events::ECL2HOR,
                site[0],
                site[1],
                site[2],
                pressure,
                temperature,
                place.values[0],
                place.values[1],
            )?;
            (
                horizontal.true_altitude,
                horizontal.azimuth,
                object_magnitude,
            )
        };
        let sun_alt = if flags & HELFLAG_VISLIM_DARK != 0 {
            -90.0
        } else {
            sun.altitude_deg
        };
        let moon_alt = if flags & HELFLAG_VISLIM_NOMOON != 0 {
            -90.0
        } else {
            moon.altitude_deg
        };
        let details = |limit, mag| {
            [
                limit,
                altitude,
                azimuth,
                sun_alt,
                sun.azimuth_deg,
                moon_alt,
                moon.azimuth_deg,
                mag,
            ]
        };
        if altitude < 0.0 {
            return Ok((HELFLAG_BELOW_HORIZON as f64, details(0.0, object_magnitude)));
        }
        let extinction = if met_range > 0.0 && met_range < 1.0 {
            met_range
        } else {
            let aerosol = if met_range >= 1.0 {
                3.912 / met_range - 0.106
            } else {
                0.1 + 0.2 * (humidity / 100.0)
            };
            crate::schaefer::RAYLEIGH_COEFF_SEA_LEVEL * (pressure / 1013.25)
                + aerosol
                + crate::schaefer::OZONE_COEFF
        };
        let airmass = if altitude <= 0.0 {
            40.0
        } else {
            let z = 90.0 - altitude;
            (1.0 / (z.to_radians().cos() + 0.50572 * (96.07995 - z).powf(-1.6364))).min(40.0)
        };
        let solar_background = if sun_alt >= 0.0 {
            1e10
        } else if sun_alt >= -6.0 {
            10.0_f64.powf(4.0 + sun_alt * 0.4)
        } else if sun_alt >= -12.0 {
            10.0_f64.powf(1.6 + (sun_alt + 6.0) * 0.2)
        } else if sun_alt >= -18.0 {
            10.0_f64.powf(0.4 + (sun_alt + 12.0) * 0.1)
        } else {
            1.0
        };
        let lunar_background = if moon_alt <= 0.0 {
            0.0
        } else {
            let illumination = match calculation.pheno_ut(jd, MOON, flags & 255) {
                Ok(value) => (1.0 + value.values[0].to_radians().cos()) / 2.0,
                Err(error) if mathematical_failure(&error) => 0.5,
                Err(error) => return Err(error.into()),
            };
            illumination * moon_alt.to_radians().sin() * 100.0
        };
        let background = solar_background + lunar_background + 1.0;
        let age_loss = ((age - 30.0) / 10.0).max(0.0) * 0.05;
        let snellen_loss = if snellen > 0.0 {
            -2.5 * snellen.log10()
        } else {
            0.0
        };
        let sky_loss = if background > 1.0 {
            2.5 * background.log10()
        } else {
            0.0
        };
        let limit = 6.5 - age_loss - snellen_loss - sky_loss;
        let status = if sun_alt >= -6.0 {
            HELFLAG_PHOTOPIC
        } else if sun_alt >= -12.0 {
            HELFLAG_MIXED
        } else {
            HELFLAG_SCOTOPIC
        };
        Ok((
            status as f64,
            details(limit, object_magnitude + extinction * airmass),
        ))
    }

    /// Compact heliacal event model with a 400-day sampled search.
    ///
    /// Time is UT1; site angles are degrees/metres and atmospheric operands
    /// mbar/Celsius/humidity fraction are retained in the public contract. This
    /// model's visibility gate uses the piecewise twilight magnitude, geometric
    /// altitude, catalog magnitude and elongation >=10 degrees; atmospheric
    /// operands do not enter that relation. Events1..4 are morning first,
    /// evening last, evening first and morning last. Sun/Moon/unknown IDs refuse.
    /// Sampling uses three hourly dawn/dusk points and a0.2-day bisection window
    /// with1e-6-day stopping width. A completed miss is `(0,-1)`; actual data or
    /// calculation errors propagate. Forced LEB cannot serve direct geometry.
    #[allow(clippy::too_many_arguments)]
    pub fn eclipse_heliacal_ut(
        &self,
        jd_start: f64,
        lat: f64,
        lon: f64,
        altitude: f64,
        pressure: f64,
        temperature: f64,
        humidity: f64,
        body: i32,
        event_type: i32,
        flags: i32,
    ) -> VisibilityResult<(f64, i32)> {
        if !(1..=4).contains(&event_type) {
            return Err(invalid(format!("invalid event_type: {event_type}")).into());
        }
        if body == SUN || body == MOON || !planet(body) {
            return Err(invalid(format!("invalid heliacal body: {body}")).into());
        }
        let _ = (pressure, temperature, humidity);
        let owner = self.eclipse_visibility_owner()?;
        let calculation = Session::from_context(owner.ctx.clone().with_mode(self.ctx.mode()));
        let search = CompactSearch {
            photometry: &calculation,
            owner: &owner,
            site: [lon, lat, altitude],
            body,
            flags,
        };
        let morning = event_type == 1 || event_type == 4;
        let hours = if morning {
            [4.0, 5.0, 6.0]
        } else {
            [18.0, 19.0, 20.0]
        };
        let mut first = 0.0;
        let mut last = 0.0;
        let mut invisible = false;
        for day in 0..400 {
            for hour in hours {
                let jd = jd_start + day as f64 + hour / 24.0;
                let (sun, body_place) = search.places(jd)?;
                if !(sun.altitude_deg > -15.0 && sun.altitude_deg < -3.0) {
                    continue;
                }
                if matches!(event_type, 1 | 2) && body_place.altitude_deg <= 1.0 {
                    continue;
                }
                let visible = search.visible(jd, sun.altitude_deg, body_place.altitude_deg)?;
                match event_type {
                    1 if visible => return Ok((search.refine(jd, true)?, event_type)),
                    2 if visible => {
                        first = jd;
                        break;
                    }
                    3 if !visible && body_place.altitude_deg < 5.0 => invisible = true,
                    3 if visible && invisible => {
                        return Ok((search.refine(jd, false)?, event_type));
                    }
                    4 if visible => last = jd,
                    4 if last > 0.0 => return Ok((search.refine(last, true)?, event_type)),
                    _ => {}
                }
            }
            if event_type == 2 && first > 0.0 {
                break;
            }
        }
        if event_type == 2 && first > 0.0 {
            last = first;
            let start_day = (first - jd_start).trunc() as usize;
            for day in start_day..start_day + 120 {
                let mut visible_day = false;
                for hour in hours {
                    let jd = jd_start + day as f64 + hour / 24.0;
                    let (sun, body_place) = search.places(jd)?;
                    if sun.altitude_deg > -15.0
                        && sun.altitude_deg < -3.0
                        && body_place.altitude_deg > 1.0
                        && search.visible(jd, sun.altitude_deg, body_place.altitude_deg)?
                    {
                        last = jd;
                        visible_day = true;
                        break;
                    }
                }
                if !visible_day && day > start_day {
                    return Ok((search.refine(last, false)?, event_type));
                }
            }
            return Ok((search.refine(last, false)?, event_type));
        }
        if event_type == 4 && last > 0.0 {
            Ok((last, event_type))
        } else {
            Ok((0.0, -1))
        }
    }
}

struct CompactSearch<'a> {
    photometry: &'a Session,
    owner: &'a Session,
    site: [f64; 3],
    body: i32,
    flags: i32,
}
impl CompactSearch<'_> {
    fn places(
        &self,
        jd: f64,
    ) -> Result<
        (
            crate::vislim::VisibilityPlace,
            crate::vislim::VisibilityPlace,
        ),
        Error,
    > {
        let sky = NativeVisibilityObserver::new(self.owner.ctx.clone(), jd, self.site)?;
        Ok((sky.observe(SUN)?, sky.observe(self.body)?))
    }
    fn visible(&self, jd: f64, sun: f64, altitude: f64) -> Result<bool, Error> {
        let elongation = match self.photometry.pheno_ut(jd, self.body, self.flags) {
            Ok(value) => value.values[2],
            Err(error) if mathematical_failure(&error) => {
                let sky = NativeVisibilityObserver::new(self.owner.ctx.clone(), jd, self.site)?;
                crate::positions::vislim::apparent_separation(
                    sky.geocentric_ray(SUN)?,
                    sky.geocentric_ray(self.body)?,
                )
            }
            Err(error) => return Err(error),
        };
        if altitude < 0.0 || sun > 0.0 {
            return Ok(false);
        }
        let limit = if sun >= 0.0 {
            -2.0
        } else if sun >= -6.0 {
            -2.0 + (sun + 6.0) / 6.0
        } else if sun >= -12.0 {
            3.0 + (sun + 12.0) * 2.0 / 6.0
        } else if sun >= -18.0 {
            5.0 + (sun + 18.0) * 1.5 / 6.0
        } else {
            6.5
        };
        Ok(magnitude(self.photometry, jd, self.body, self.flags)? <= limit && elongation >= 10.0)
    }
    fn refine(&self, jd: f64, morning: bool) -> Result<f64, Error> {
        let (mut lower, mut upper) = (jd - 0.1, jd + 0.1);
        for _ in 0..30 {
            let middle = (lower + upper) / 2.0;
            let (sun, body) = self.places(middle)?;
            let visible = self.visible(middle, sun.altitude_deg, body.altitude_deg)?;
            if visible == morning {
                upper = middle;
            } else {
                lower = middle;
            }
            if upper - lower < 1e-6 {
                break;
            }
        }
        Ok((lower + upper) / 2.0)
    }
}
