//! One independent bundled-timescale surface observer observer for visual sky samples.
use super::context_serving::{
    native_pheno_ecliptic_frame, native_pheno_frame_angles, native_pheno_itrs_position,
};
use super::heliacal_star::{CatalogHeliacalStar, NativeHeliacalStar};
use super::reference_frames::*;
use super::{LocalSpkStateProvider, SiteEllipsoid};
use crate::context::EphemerisContext;
use crate::errors::Result;
use crate::types::TwoPartEpoch;
use crate::vislim::VisibilityPlace;

/// Ecliptic-to-horizontal reduction for stored visibility geometry.
///
/// Input epoch is UT1 JD; site and ecliptic angles are degrees, height metres,
/// pressure mbar, temperature Celsius and lapse K/m. The configured Delta-T
/// callback supplies TT-minus-UT1 days for the sidereal model. Native frame,
/// refraction and sidereal kernels perform the reduction without acquiring a
/// body state or main ephemeris image; their errors propagate.
pub(crate) fn stored_visibility_horizontal(
    jd: f64,
    site: [f64; 3],
    atmo: [f64; 4],
    place: [f64; 3],
    lapse: f64,
    delta: &dyn Fn(f64) -> Result<f64>,
) -> Result<VisibilityPlace> {
    let tt = crate::time::default_timescale_tt_from_ut1(jd)?;
    let eps0 = super::vondrak_mean_obliquity_deg(tt)?;
    let (psi, deps) = crate::erfa::nut06a(tt)?;
    let eps = eps0 + deps.to_degrees();
    let (ra, dec) = super::cotrans(place[0], place[1], -eps)?;
    // The public _gmst06 returns radians, then _sidtime_internal degrees.
    let gmst = crate::time::public_sidereal::mean_deg(jd, delta)?
        .to_radians()
        .to_degrees()
        / 15.0;
    let equation = (psi.to_degrees() * eps.to_radians().cos()) / 15.0;
    let lst = ((gmst + equation) + site[0] / 15.0).rem_euclid(24.0);
    let hour = (lst * 15.0 - ra).rem_euclid(360.0).to_radians();
    let dec = dec.to_radians();
    let lat = site[1].to_radians();
    let sine = lat.sin() * dec.sin() + lat.cos() * dec.cos() * hour.cos();
    let sine = sine.clamp(-1.0, 1.0);
    let altitude = sine.asin().to_degrees();
    let y = hour.sin() * dec.cos();
    let x = hour.cos() * dec.cos() * lat.sin() - dec.sin() * lat.cos();
    let azimuth = y.atan2(x).to_degrees().rem_euclid(360.0);
    let pressure = if atmo[0] <= 0.0 {
        crate::refraction::estimate_pressure_mbar(site[2])
    } else {
        atmo[0]
    };
    // Evaluate atmospheric refraction while retaining true altitude for the visibility model.
    crate::refraction::refrac_extended_true_to_app(altitude, site[2], pressure, atmo[1], lapse);
    Ok(VisibilityPlace {
        azimuth_deg: (azimuth + 180.0) % 360.0,
        altitude_deg: altitude,
        ray: place,
    })
}

/// Owned input configuration and one original UT1 Time, WGS84 site and
/// latitude/longitude rotation. Each observation rereads its actual Earth
/// state before target light time; no caller model TT replaces the Time.
pub(crate) struct NativeVisibilityObserver {
    context: EphemerisContext,
    epochs: VectorRotationEpochs,
    scalar_tt: f64,
    site_itrs: [f64; 3],
    latlon: [[f64; 3]; 3],
    site: [f64; 3],
    heliacal_star: Option<NativeHeliacalStar>,
    held_heliacal_frame: Option<([[f64; 3]; 3], f64)>,
}

impl NativeVisibilityObserver {
    /// Construct the true WGS84 surface and resolve Earth/Sun/Moon aliases
    /// before the independent UT1 Time. The context owns immutable images;
    /// its clones share their bytes. East longitude/latitude are degrees,
    /// height metres. Site and raw provider/time errors remain in read order.
    pub(crate) fn new(context: EphemerisContext, jd_ut1: f64, site: [f64; 3]) -> Result<Self> {
        let site_itrs =
            native_pheno_itrs_position((site[0], site[1], site[2]), SiteEllipsoid::Wgs84)?;
        let provider = LocalSpkStateProvider::new(context.spk_registry(), context.mode())?;
        provider.validate_source_alias(crate::constants::EARTH)?;
        provider.validate_source_alias(crate::constants::SUN)?;
        provider.validate_source_alias(crate::constants::MOON)?;
        let (epochs, scalar_tt) = VectorRotationEpochs::from_ut1(jd_ut1)?;
        let phi = site[1] / 360.0 * FRAME_TAU;
        let zero = phi * 0.0;
        let one = zero + 1.0;
        let lat = [
            [-scalar_sin(phi), zero, scalar_cos(phi)],
            [zero, one, zero],
            [scalar_cos(phi), zero, scalar_sin(phi)],
        ];
        let latlon = matmul3(lat, rot_z(-(site[0] / 360.0 * FRAME_TAU)))?;
        Ok(Self {
            context,
            epochs,
            scalar_tt,
            site_itrs,
            latlon,
            site,
            heliacal_star: None,
            held_heliacal_frame: None,
        })
    }

    /// Heliacal acquisition is Sun, Earth, then the selected numeric target or
    /// direct catalog Star, before the WGS84 site and original UT1 Time. The
    /// Star keeps RA/declination/proper motion with omitted parallax/radial
    /// fields; no SPK star alias or public fixstar path is selected. Moon is
    /// acquired only by a consumed refinement, so twilight needs no Moon image.
    pub(crate) fn new_for_heliacal(
        context: EphemerisContext,
        jd_ut1: f64,
        site: [f64; 3],
        body: i32,
    ) -> Result<Self> {
        let catalog = admit_heliacal_catalog(&context, body)?;
        let heliacal_star = catalog
            .map(|catalog| NativeHeliacalStar::new(context.clone(), catalog))
            .transpose()?;
        let site_itrs =
            native_pheno_itrs_position((site[0], site[1], site[2]), SiteEllipsoid::Wgs84)?;
        let (epochs, scalar_tt) = VectorRotationEpochs::from_ut1(jd_ut1)?;
        let phi = site[1] / 360.0 * FRAME_TAU;
        let zero = phi * 0.0;
        let one = zero + 1.0;
        let lat = [
            [-scalar_sin(phi), zero, scalar_cos(phi)],
            [zero, one, zero],
            [scalar_cos(phi), zero, scalar_sin(phi)],
        ];
        let latlon = matmul3(lat, rot_z(-(site[0] / 360.0 * FRAME_TAU)))?;
        Ok(Self {
            context,
            epochs,
            scalar_tt,
            site_itrs,
            latlon,
            site,
            heliacal_star,
            held_heliacal_frame: None,
        })
    }

    /// Retained direct-Star identity, canonical name and magnitude, selected
    /// before site/Time construction. Numeric and fixed-vislim owners have no
    /// direct catalog carrier. This metadata query acquires no state or clock.
    pub(crate) fn heliacal_catalog(&self) -> Option<CatalogHeliacalStar> {
        self.heliacal_star.as_ref().map(NativeHeliacalStar::catalog)
    }

    /// Acquire and retain original epoch's equinox matrix and GAST, hours,
    /// before heliacal Sun geometry. Subsequent rotations/date reductions use
    /// that same record, including its actual GAST word. A failed acquisition
    /// installs nothing. No provider state, configured model clock or new Time
    /// construction is consumed; the existing default frame guards apply.
    pub(crate) fn acquire_heliacal_gast(&mut self) -> Result<f64> {
        if let Some((_, gast)) = self.held_heliacal_frame {
            return Ok(gast);
        }
        let frame = equinox_and_gast_at_parts(
            self.epochs.whole,
            self.epochs.tt_fraction,
            self.epochs.ut1_fraction,
            self.scalar_tt,
        )?;
        self.held_heliacal_frame = Some(frame);
        Ok(frame.1)
    }

    /// Reuse a successfully retained heliacal frame, or preserve the original
    /// independent-Time frame acquisition for other observer paths.
    fn original_equinox_and_gast(&self) -> Result<([[f64; 3]; 3], f64)> {
        match self.held_heliacal_frame {
            Some(frame) => Ok(frame),
            None => equinox_and_gast_at_parts(
                self.epochs.whole,
                self.epochs.tt_fraction,
                self.epochs.ut1_fraction,
                self.scalar_tt,
            ),
        }
    }

    /// Original-Time ITRS rotation. A retained heliacal record uses the exact
    /// existing Rz(-GAST)@M primitive; the fixed-vislim/uncached tail is the
    /// unchanged split-Time constructor, without a collapsed epoch substitute.
    fn original_rotation(&self) -> Result<[[f64; 3]; 3]> {
        match self.held_heliacal_frame {
            Some((equinox, gast)) => itrs_rotation_matrix(equinox, gast),
            None => itrs_rotation_at_parts(
                self.epochs.whole,
                self.epochs.tt_fraction,
                self.epochs.ut1_fraction,
                self.scalar_tt,
            ),
        }
    }

    /// Select only this heliacal entry's retained direct-Star target. Sun and
    /// Moon samples retain the existing numeric acquisition, and the separate
    /// fixed-vislim owner never enters the direct-Star channel.
    fn direct_star(&self, body: i32) -> Option<&NativeHeliacalStar> {
        self.heliacal_star
            .as_ref()
            .filter(|star| star.catalog().id() == body)
    }

    /// Apparent geocentric ICRS AU ray at the original epoch. Earth position
    /// and velocity precede target light time and apparent corrections. There
    /// is no surface observer surface or Earth-deflection channel in this observation.
    pub(crate) fn geocentric_ray(&self, body: i32) -> Result<[f64; 3]> {
        let provider =
            LocalSpkStateProvider::new(self.context.spk_registry(), self.context.mode())?;
        let time = TwoPartEpoch {
            hi: self.epochs.whole,
            lo: self.epochs.tt_fraction,
        };
        let earth = provider.star_source_time_state(crate::constants::EARTH, time)?;
        if let Some(star) = self.direct_star(body) {
            return star.apparent_geocentric_ray(time, earth.0, earth.1);
        }
        provider.eclipse_source_observed_vector(body, time, earth.0, earth.1, None)
    }

    /// Apparent longitude in the fixed IAU 1976/SPICE ECLIPJ2000 ecliptic,
    /// in degrees. The geocentric ray uses that stored rotation, with no moving
    /// frame or IAU 2006 frame-bias substitution.
    pub(crate) fn geocentric_longitude(&self, body: i32) -> Result<f64> {
        Ok(native_pheno_frame_angles(
            super::routing::iau1976_ecliptic_j2000_matrix(),
            self.geocentric_ray(body)?,
        )?
        .0)
    }

    /// Geometric target-to-site ray and north-through-east horizon angles.
    ///
    /// Barycentric target and Earth/site use the same split observation epoch;
    /// planetary geometry has no retardation, aberration or deflection. A stellar
    /// row uses its propagated astrometric ray and the true-equator-of-date frame,
    /// then GAST and the site's latitude reduce RA/declination to horizon angles.
    /// The reported ray is ICRS AU, before apparent corrections.
    pub(crate) fn geometric_place(&self, body: i32) -> Result<VisibilityPlace> {
        let provider =
            LocalSpkStateProvider::new(self.context.spk_registry(), self.context.mode())?;
        let time = TwoPartEpoch {
            hi: self.epochs.whole,
            lo: self.epochs.tt_fraction,
        };
        if let Some(star) = self.direct_star(body) {
            let earth = provider.star_source_time_state(crate::constants::EARTH, time)?;
            let rotation = self.original_rotation()?;
            let (surface, _) = itrs_to_gcrs(rotation, self.site_itrs)?;
            let observer = std::array::from_fn(|i| earth.0[i] + surface[i]);
            let astrometric = star.astrometric(time, observer);
            let (frame, gast) = self.original_equinox_and_gast()?;
            let (ra_hours, dec_deg) = date_equatorial(matvec3(frame, astrometric.ray)?);
            let (altitude_deg, azimuth_deg) = crate::heliacal::altaz_from_radec(
                ra_hours,
                dec_deg,
                self.site[0],
                self.site[1],
                gast,
            );
            return Ok(VisibilityPlace {
                azimuth_deg,
                altitude_deg,
                ray: astrometric.ray,
            });
        }
        let target = provider.star_source_time_state(body, time)?;
        let earth = provider.star_source_time_state(crate::constants::EARTH, time)?;
        let rotation = self.original_rotation()?;
        let (surface, _) = itrs_to_gcrs(rotation, self.site_itrs)?;
        let ray = std::array::from_fn(|i| target.0[i] - (earth.0[i] + surface[i]));
        let horizon = matmul3(self.latlon, rotation)?;
        let (azimuth_deg, altitude_deg) = native_pheno_frame_angles(horizon, ray)?;
        Ok(VisibilityPlace {
            azimuth_deg,
            altitude_deg,
            ray,
        })
    }

    /// Geometric geocentric altitude in degrees from of-date RA/Dec and GAST.
    /// The ray omits light time, aberration and surface parallax. The derived
    /// altitude sine is limited to `[-1, 1]` for inverse-trigonometric roundoff;
    /// NaN remains NaN. Stellar rows acquire and propagate their geocentric
    /// astrometry separately from the topocentric observation.
    pub(crate) fn geocentric_altitude(&self, body: i32) -> Result<f64> {
        let provider =
            LocalSpkStateProvider::new(self.context.spk_registry(), self.context.mode())?;
        let time = TwoPartEpoch {
            hi: self.epochs.whole,
            lo: self.epochs.tt_fraction,
        };
        if let Some(star) = self.direct_star(body) {
            let earth = provider.star_source_time_state(crate::constants::EARTH, time)?;
            let astrometric = star.astrometric(time, earth.0);
            let (frame, gast) = self.original_equinox_and_gast()?;
            let (ra_hours, dec_deg) = date_equatorial(matvec3(frame, astrometric.ray)?);
            let lst = gast + self.site[0] / 15.0;
            let hour_angle = ((lst - ra_hours) * 15.0).to_radians();
            let dec = dec_deg.to_radians();
            let lat = self.site[1].to_radians();
            let sine = lat.sin() * dec.sin() + lat.cos() * dec.cos() * hour_angle.cos();
            let sine = sine.clamp(-1.0, 1.0);
            return Ok(sine.asin().to_degrees());
        }
        let target = provider.star_source_time_state(body, time)?;
        let earth = provider.star_source_time_state(crate::constants::EARTH, time)?;
        let ray = std::array::from_fn(|i| target.0[i] - earth.0[i]);
        let (frame, gast) = self.original_equinox_and_gast()?;
        let [x, y, z] = matvec3(frame, ray)?;
        let ra_hours = floored_mod(y.atan2(x), FRAME_TAU) * 24.0 / FRAME_TAU;
        let dec_deg = z.atan2(x.hypot(y)) * 360.0 / FRAME_TAU;
        let lst = gast + self.site[0] / 15.0;
        let hour_angle = ((lst - ra_hours) * 15.0).to_radians();
        let dec = dec_deg.to_radians();
        let lat = self.site[1].to_radians();
        let sine = lat.sin() * dec.sin() + lat.cos() * dec.cos() * hour_angle.cos();
        Ok(sine.clamp(-1.0, 1.0).asin().to_degrees())
    }

    /// Retain one target's real apparent ICRS AU ray and unrefracted
    /// north-through-east horizontal coordinates. Earth and surface
    /// position/velocity are composed before the target's ten-step light-time
    /// and strict apparent corrections; Some(actual_surface) keeps the
    /// surface observer Earth-deflection channel. No public scalar TT envelope or output
    /// flags are applied. Split provider and frame errors propagate unchanged.
    pub(crate) fn observe(&self, body: i32) -> Result<VisibilityPlace> {
        let provider =
            LocalSpkStateProvider::new(self.context.spk_registry(), self.context.mode())?;
        let time = TwoPartEpoch {
            hi: self.epochs.whole,
            lo: self.epochs.tt_fraction,
        };
        let earth = provider.star_source_time_state(crate::constants::EARTH, time)?;
        let rotation = self.original_rotation()?;
        let (surface, site_velocity) = itrs_to_gcrs(rotation, self.site_itrs)?;
        let observer = std::array::from_fn(|i| earth.0[i] + surface[i]);
        let velocity = std::array::from_fn(|i| earth.1[i] + site_velocity[i]);
        if let Some(star) = self.direct_star(body) {
            return star.apparent_topocentric_place(time, observer, velocity, surface, || {
                matmul3(self.latlon, rotation)
            });
        }
        let ray = provider.eclipse_source_observed_vector(
            body,
            time,
            observer,
            velocity,
            Some(surface),
        )?;
        let horizon = matmul3(self.latlon, rotation)?;
        let (azimuth_deg, altitude_deg) = native_pheno_frame_angles(horizon, ray)?;
        Ok(VisibilityPlace {
            azimuth_deg,
            altitude_deg,
            ray,
        })
    }

    /// True-ecliptic-of-date longitude/latitude of an already retained ray,
    /// degrees, using this same epoch's independent frame clock.
    pub(crate) fn ecliptic(&self, ray: [f64; 3]) -> Result<(f64, f64)> {
        let frame = native_pheno_ecliptic_frame(TwoPartEpoch {
            hi: self.epochs.whole,
            lo: self.epochs.tt_fraction,
        })?;
        native_pheno_frame_angles(frame, ray)
    }
}

/// Resolve heliacal Sun, Earth and target aliases before observer optics and
/// atmospheric model preparation. This epoch-free admission retains the owned
/// registry, performs no sample/Time acquisition and deliberately skips Moon.
pub(crate) fn admit_heliacal_target(context: &EphemerisContext, body: i32) -> Result<()> {
    admit_heliacal_catalog(context, body).map(|_| ())
}

/// Resolve Sun/Earth before the selected catalog model or numeric target.
/// The catalog constructor owns missing-Star metadata and consumes no Time;
/// numeric aliases retain the same raw provider lookup and refusal order.
fn admit_heliacal_catalog(
    context: &EphemerisContext,
    body: i32,
) -> Result<Option<CatalogHeliacalStar>> {
    let provider = LocalSpkStateProvider::new(context.spk_registry(), context.mode())?;
    admit_heliacal_catalog_via(
        body,
        |id| provider.validate_source_alias(id),
        CatalogHeliacalStar::from_id,
    )
}

/// Ordered epoch-free constituent admission. Separate primitive producers
/// retain the actual alias/catalog errors and never consume an unused arm.
fn admit_heliacal_catalog_via(
    body: i32,
    mut alias_at: impl FnMut(i32) -> Result<()>,
    catalog_at: impl FnOnce(i32) -> Result<CatalogHeliacalStar>,
) -> Result<Option<CatalogHeliacalStar>> {
    alias_at(crate::constants::SUN)?;
    alias_at(crate::constants::EARTH)?;
    if crate::heliacal::is_fixed_star(body) {
        catalog_at(body).map(Some)
    } else {
        alias_at(body)?;
        Ok(None)
    }
}

/// Of-date equatorial ray to RA hours and signed declination degrees.
/// Preserve spherical atan2/hypot and Angle's multiply-then-divide conversion;
/// no radius normalization, post-arithmetic guard or full-turn reclassification
/// is introduced. The frame belongs to the caller's original split epoch.
fn date_equatorial([x, y, z]: [f64; 3]) -> (f64, f64) {
    (
        floored_mod(y.atan2(x), FRAME_TAU) * 24.0 / FRAME_TAU,
        z.atan2(x.hypot(y)) * 360.0 / FRAME_TAU,
    )
}

/// Kahan half-angle separation of two retained ICRS AU rays, degrees.
/// Scale each ray by the other's norm, then use 2 atan2(norm(a-b),norm(a+b)).
/// All norms/products keep left-associated three-component arithmetic.
/// Derived zero/nonfinite rays retain IEEE behavior; no clamp or physical
/// correction replaces the acquired apparent vectors.
pub(crate) fn apparent_separation(u: [f64; 3], v: [f64; 3]) -> f64 {
    let length = |p: [f64; 3]| (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
    let lu = length(u);
    let lv = length(v);
    let a: [f64; 3] = std::array::from_fn(|i| u[i] * lv);
    let b: [f64; 3] = std::array::from_fn(|i| v[i] * lu);
    let minus = std::array::from_fn(|i| a[i] - b[i]);
    let plus = std::array::from_fn(|i| a[i] + b[i]);
    (2.0 * length(minus).atan2(length(plus))).to_degrees()
}

#[cfg(test)]
mod heliacal_star_integration_tests {
    use super::*;
    use crate::constants::{EARTH, MERCURY, SIRIUS, SUN};
    use crate::errors::Error;
    use std::cell::{Cell, RefCell};

    #[test]
    fn stored_horizontal_consumes_only_modern_sidereal_clock() {
        let reads = RefCell::new(Vec::new());
        let value = stored_visibility_horizontal(
            2460400.25,
            [12.0, 42.0, 50.0],
            [1013.25, 15.0, 50.0, 0.0],
            [80.0, 5.0, 1.0],
            0.0065,
            &|jd| {
                reads.borrow_mut().push(jd);
                Ok(0.001)
            },
        )
        .unwrap();
        assert_eq!(*reads.borrow(), [2460400.25]);
        assert!(value.azimuth_deg.is_finite() && value.altitude_deg.is_finite());
        assert_eq!(value.ray, [80.0, 5.0, 1.0]);
    }

    #[test]
    fn stored_horizontal_keeps_its_actual_clock_refusal() {
        let marker = Error::ConfigurationMissing {
            missing_config: "actual consumed horizon clock",
            suggestion: "same owning clock refusal",
        };
        let result = stored_visibility_horizontal(
            2460400.25,
            [12.0, 42.0, 50.0],
            [1013.25, 15.0, 50.0, 0.0],
            [80.0, 5.0, 1.0],
            0.0065,
            &|jd| {
                assert_eq!(jd, 2460400.25);
                Err(marker.clone())
            },
        );
        assert!(matches!(result, Err(error) if error == marker));
    }

    #[test]
    fn stellar_admission_is_sun_earth_catalog_without_a_star_alias() {
        let reads = RefCell::new(Vec::new());
        let star = admit_heliacal_catalog_via(
            SIRIUS,
            |body| {
                assert_ne!(body, SIRIUS, "a catalog Star has no SPK target alias");
                reads.borrow_mut().push((false, body));
                Ok(())
            },
            |body| {
                reads.borrow_mut().push((true, body));
                CatalogHeliacalStar::from_id(body)
            },
        )
        .unwrap()
        .unwrap();
        assert_eq!(star.id(), SIRIUS);
        assert_eq!(
            *reads.borrow(),
            [(false, SUN), (false, EARTH), (true, SIRIUS)]
        );
    }

    #[test]
    fn numeric_admission_never_consumes_a_catalog_constructor() {
        let reads = RefCell::new(Vec::new());
        let target = admit_heliacal_catalog_via(
            MERCURY,
            |body| {
                reads.borrow_mut().push(body);
                Ok(())
            },
            |_| panic!("numeric acquisition consumed a catalog constructor"),
        )
        .unwrap();
        assert!(target.is_none());
        assert_eq!(*reads.borrow(), [SUN, EARTH, MERCURY]);
    }

    #[test]
    fn alias_refusal_precedes_catalog_and_retains_its_actual_error() {
        let reads = RefCell::new(Vec::new());
        let catalog_read = Cell::new(false);
        let raw = Error::SpkUnregisteredBody { body: EARTH };
        let actual = admit_heliacal_catalog_via(
            SIRIUS,
            |body| {
                reads.borrow_mut().push(body);
                if body == EARTH {
                    Err(raw.clone())
                } else {
                    Ok(())
                }
            },
            |body| {
                catalog_read.set(true);
                CatalogHeliacalStar::from_id(body)
            },
        );
        assert_eq!(actual.unwrap_err(), raw);
        assert_eq!(*reads.borrow(), [SUN, EARTH]);
        assert!(!catalog_read.get());
    }

    #[test]
    fn date_equatorial_keeps_hour_units_wrap_and_raw_declination() {
        assert_eq!(date_equatorial([0.0, 1.0, 0.0]), (6.0, 0.0));
        assert_eq!(date_equatorial([0.0, -1.0, 0.0]), (18.0, 0.0));
        assert_eq!(date_equatorial([0.0, 0.0, 1.0]), (0.0, 90.0));
        assert_eq!(date_equatorial([0.0, 0.0, -1.0]), (0.0, -90.0));
        let (ra, dec) = date_equatorial([f64::NAN, 0.0, 0.0]);
        assert!(ra.is_nan() && dec.is_nan());
    }

    #[test]
    fn held_gast_record_is_reused_without_an_unused_clock_acquisition() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let gast = 5.75_f64;
        // Invalid unconsumed epoch words make an accidental reacquisition
        // refuse. Production installs this record only after checked admission.
        let mut sky = NativeVisibilityObserver {
            context: EphemerisContext::new(),
            epochs: VectorRotationEpochs {
                whole: f64::NAN,
                tt_fraction: f64::NAN,
                ut1_fraction: f64::NAN,
            },
            scalar_tt: f64::NAN,
            site_itrs: [0.0; 3],
            latlon: identity,
            site: [0.0; 3],
            heliacal_star: None,
            held_heliacal_frame: Some((identity, gast)),
        };
        assert_eq!(
            sky.acquire_heliacal_gast().unwrap().to_bits(),
            gast.to_bits()
        );
        assert_eq!(sky.original_equinox_and_gast().unwrap(), (identity, gast));
        assert_eq!(
            sky.original_rotation().unwrap(),
            itrs_rotation_matrix(identity, gast).unwrap()
        );
        sky.held_heliacal_frame = None;
        assert!(sky.original_rotation().is_err());
        assert!(sky.acquire_heliacal_gast().is_err());
        assert!(sky.held_heliacal_frame.is_none());
    }
}
