//! Fixed-instant visual limiting magnitude with ordered external acquisitions.
use crate::errors::{Error, Result};
use crate::heliacal::{VislimBody, VislimGeometry, VislimModel, vis_limit_mag_words};

/// Optical model retained for one fixed-instant visibility calculation.
///
/// Geometry supplies unrefracted altitudes, north-zero eastward azimuths and
/// separations in degrees, visual-model magnitude and lunar illuminated fraction.
/// `finish` returns a numerical **status flag**, not the limiting magnitude,
/// followed by the ten values documented on [`crate::heliacal::vis_limit_mag_words`].
/// The ordinary model uses 0 for photopic, 1 for scotopic and −2 for below
/// horizon; limiting magnitude occupies diagnostic slot 0. Implementations own
/// their data and return typed calculation failures.
pub trait VisibilityOpticalModel: Send + Sync {
    /// Complete a visibility reduction using the supplied physical geometry.
    /// No foreign arithmetic, interpreter object or mutable callback state is
    /// required; the implementation remains native Rust.
    fn finish(
        &self,
        flags: i32,
        object: VislimBody,
        geometry: &VislimGeometry,
    ) -> Result<(f64, [f64; 10])>;
}

/// Owned apparent sky place: north-through-east azimuth and true altitude in
/// degrees, with a retained apparent ICRS AU ray or ecliptic star direction.
#[derive(Debug, Clone, Copy)]
pub struct VisibilityPlace {
    /// North-zero eastward azimuth, degrees.
    pub azimuth_deg: f64,
    /// Geometric altitude, degrees; negative is the below-horizon sentinel.
    pub altitude_deg: f64,
    /// Observer-relative apparent ICRS AU ray for planets. The prepared star
    /// carrier instead retains ecliptic longitude/latitude degrees and unit
    /// distance until its public horizontal/separation stages consume them.
    pub ray: [f64; 3],
}

/// Ordered acquisition provider for one limiting-magnitude call.
/// Main/Time/observer configuration is retained by the implementation; these
/// methods own each sample and never return a borrowed kernel or foreign value.
pub trait VisibilityFeed {
    /// Sun apparent true-altitude/azimuth acquisition.
    fn sun(&mut self) -> Result<VisibilityPlace>;
    /// Moon apparent true-altitude/azimuth, after the Sun.
    fn moon(&mut self) -> Result<VisibilityPlace>;
    /// Resolve the input object only after both sky-brightness bodies.
    fn object(&mut self) -> Result<VislimBody>;
    /// Apparent object position, after parsing/validation.
    fn object_place(&mut self, object: VislimBody) -> Result<VisibilityPlace>;
    /// Star horizon projection after its magnitude lookup. Planet providers
    /// already project in object_place; their implementation is not called here.
    fn star_horizontal(&mut self, place: VisibilityPlace) -> Result<VisibilityPlace>;
    /// Visual magnitude; only the defined arithmetic/value/type catch can
    /// choose the numeric fallback (zero for planets, two for stars).
    fn magnitude(&mut self, object: VislimBody) -> Result<f64>;
    /// Atmosphere and observer optics after the above-horizon decision.
    fn model(&mut self) -> Result<VislimModel>;
    /// Optional live model acquired only after the above-horizon decision.
    /// None preserves the existing standalone scalar model and its guards.
    fn acquired_model(&mut self) -> Result<Option<Box<dyn VisibilityOpticalModel>>> {
        Ok(None)
    }
    /// Geocentric Moon illuminated fraction; its error-family handling differs
    /// from object-magnitude handling.
    fn moon_phase(&mut self) -> Result<f64>;
    /// Acquire the fresh object ray or both Moon/Sun ecliptic frames before
    /// assigning either separation. A refused second frame leaves both angle
    /// defaults intact while preserving the already acquired Moon phase.
    fn prepare_separations(
        &mut self,
        object: VislimBody,
        place: VisibilityPlace,
        moon: VisibilityPlace,
        sun: VisibilityPlace,
    ) -> Result<()>;
    /// Object-Moon angular separation, degrees, after phase acquisition.
    fn moon_separation(
        &mut self,
        object: VislimBody,
        place: VisibilityPlace,
        moon: VisibilityPlace,
    ) -> Result<f64>;
    /// Object-Sun angular separation after the Moon separation.
    fn sun_separation(
        &mut self,
        object: VislimBody,
        place: VisibilityPlace,
        sun: VisibilityPlace,
    ) -> Result<f64>;
    /// True only for actually retained ValueError/TypeError/ArithmeticError
    /// constructors at the object-magnitude stage, not arbitrary library errors.
    fn catches_magnitude(&self, error: &Error) -> bool;
    /// True only for the retained ValueError/TypeError/AttributeError families
    /// at the later phase/separation stage.
    fn catches_phase(&self, error: &Error) -> bool;
}

/// Compose a fixed-instant limiting visual magnitude.
/// Sun and Moon precede target parsing, target position precedes magnitude,
/// and a negative target altitude returns flag -2 with (-100, nine zeros)
/// before atmosphere/phase/separation. Otherwise model construction precedes
/// phase, Moon separation and Sun separation. A caught late error retains
/// partial assignments already completed, rather than resetting all defaults.
/// Dark-sky/no-Moon/Moon-object branches and the ten-float output are supplied
/// by the shared Schaefer model composition. Return order is `(status, words)`,
/// with limiting magnitude in `words[0]`, not in the first tuple field; units
/// and layout follow [`crate::heliacal::vis_limit_mag_words`]. An acquired model
/// is owned for this call and prevents construction of the scalar model.
/// Uncaught acquisition/model errors propagate; no iteration or event search
/// occurs. The feed owns the observation epoch and observer configuration.
pub fn evaluate_visibility(feed: &mut impl VisibilityFeed, flags: i32) -> Result<(f64, [f64; 10])> {
    let sun = feed.sun()?;
    let moon = feed.moon()?;
    let object = feed.object()?;
    let mut place = feed.object_place(object)?;
    let magnitude = match feed.magnitude(object) {
        Ok(value) => value,
        Err(error) if feed.catches_magnitude(&error) => match object {
            VislimBody::Planet(_) => 0.0,
            VislimBody::Star(_) | VislimBody::UnresolvedStar => 2.0,
        },
        Err(error) => return Err(error),
    };
    if matches!(object, VislimBody::Star(_) | VislimBody::UnresolvedStar) {
        place = feed.star_horizontal(place)?;
    }
    if place.altitude_deg < 0.0 {
        let mut words = [0.0; 10];
        words[0] = -100.0;
        return Ok((crate::constants::HELFLAG_BELOW_HORIZON as f64, words));
    }
    let acquired = feed.acquired_model()?;
    let model = if acquired.is_none() {
        Some(feed.model()?)
    } else {
        None
    };
    let mut moon_phase = 0.0;
    let mut moon_angle = 180.0;
    let mut sun_angle = 180.0;
    let extras = (|| -> Result<()> {
        moon_phase = feed.moon_phase()?;
        feed.prepare_separations(object, place, moon, sun)?;
        moon_angle = feed.moon_separation(object, place, moon)?;
        sun_angle = feed.sun_separation(object, place, sun)?;
        Ok(())
    })();
    if let Err(error) = extras {
        if !feed.catches_phase(&error) {
            return Err(error);
        }
    }
    let geometry = VislimGeometry {
        obj_alt_deg: place.altitude_deg,
        obj_az_deg: place.azimuth_deg,
        obj_mag: magnitude,
        sun_alt_deg: sun.altitude_deg,
        sun_az_deg: sun.azimuth_deg,
        moon_alt_deg: moon.altitude_deg,
        moon_az_deg: moon.azimuth_deg,
        moon_phase,
        sun_obj_angle_deg: sun_angle,
        moon_obj_angle_deg: moon_angle,
    };
    match acquired {
        Some(model) => model.finish(flags, object, &geometry),
        None => {
            let (status, words) = vis_limit_mag_words(
                model.as_ref().expect("scalar model retained"),
                flags,
                object,
                &geometry,
            );
            Ok((status as f64, words))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Marked {
        trace: Vec<&'static str>,
        star: bool,
        below: bool,
        failed: Option<&'static str>,
        live_model: bool,
    }
    impl Marked {
        fn at(&mut self, name: &'static str) -> Result<()> {
            self.trace.push(name);
            if self.failed == Some(name) {
                return Err(Error::MathDomain { value: 2.0 });
            }
            Ok(())
        }
        fn place(&self) -> VisibilityPlace {
            VisibilityPlace {
                azimuth_deg: 80.0,
                altitude_deg: if self.below { -1.0 } else { 20.0 },
                ray: [1.0, 0.0, 0.0],
            }
        }
    }
    impl VisibilityFeed for Marked {
        fn sun(&mut self) -> Result<VisibilityPlace> {
            self.at("sun")?;
            Ok(self.place())
        }
        fn moon(&mut self) -> Result<VisibilityPlace> {
            self.at("moon")?;
            Ok(self.place())
        }
        fn object(&mut self) -> Result<VislimBody> {
            self.at("object")?;
            Ok(if self.star {
                VislimBody::Star(1000001)
            } else {
                VislimBody::Planet(crate::constants::VENUS)
            })
        }
        fn object_place(&mut self, _: VislimBody) -> Result<VisibilityPlace> {
            self.at("place")?;
            Ok(self.place())
        }
        fn star_horizontal(&mut self, place: VisibilityPlace) -> Result<VisibilityPlace> {
            self.at("star-horizon")?;
            Ok(place)
        }
        fn magnitude(&mut self, _: VislimBody) -> Result<f64> {
            self.at("magnitude")?;
            Ok(1.0)
        }
        fn model(&mut self) -> Result<VislimModel> {
            assert!(
                !self.live_model,
                "live model must not build the scalar fallback"
            );
            self.at("model")?;
            Ok(VislimModel::new(
                &crate::heliacal::VislimAtmosphere::standard(0.0, 0.0, 2451545.0)?,
                crate::heliacal::HeliacalObserver::naked_eye(),
                true,
            ))
        }
        fn acquired_model(&mut self) -> Result<Option<Box<dyn VisibilityOpticalModel>>> {
            if !self.live_model {
                return Ok(None);
            }
            self.at("live-model")?;
            struct Held;
            impl VisibilityOpticalModel for Held {
                fn finish(
                    &self,
                    _: i32,
                    _: VislimBody,
                    g: &VislimGeometry,
                ) -> Result<(f64, [f64; 10])> {
                    let mut words = [0.0; 10];
                    words[0] = g.moon_phase;
                    words[1] = g.obj_alt_deg;
                    Ok((1.0, words))
                }
            }
            Ok(Some(Box::new(Held)))
        }
        fn moon_phase(&mut self) -> Result<f64> {
            self.at("phase")?;
            Ok(0.75)
        }
        fn prepare_separations(
            &mut self,
            _: VislimBody,
            _: VisibilityPlace,
            _: VisibilityPlace,
            _: VisibilityPlace,
        ) -> Result<()> {
            self.at("separation-setup")
        }
        fn moon_separation(
            &mut self,
            _: VislimBody,
            _: VisibilityPlace,
            _: VisibilityPlace,
        ) -> Result<f64> {
            self.at("moon-separation")?;
            Ok(60.0)
        }
        fn sun_separation(
            &mut self,
            _: VislimBody,
            _: VisibilityPlace,
            _: VisibilityPlace,
        ) -> Result<f64> {
            self.at("sun-separation")?;
            Ok(90.0)
        }
        fn catches_magnitude(&self, _: &Error) -> bool {
            true
        }
        fn catches_phase(&self, _: &Error) -> bool {
            true
        }
    }
    fn feed() -> Marked {
        Marked {
            trace: vec![],
            star: false,
            below: false,
            failed: None,
            live_model: false,
        }
    }
    #[test]
    fn below_horizon_skips_atmosphere_phase_and_keeps_exact_sentinel() {
        let mut f = feed();
        f.below = true;
        let (flag, words) = evaluate_visibility(&mut f, 0).unwrap();
        assert_eq!(flag, -2.0);
        assert_eq!(words[0], -100.0);
        assert!(words[1..].iter().all(|x| x.to_bits() == 0));
        assert_eq!(f.trace, ["sun", "moon", "object", "place", "magnitude"]);
    }
    #[test]
    fn star_magnitude_precedes_horizon_and_below_horizon_decision() {
        let mut f = feed();
        f.star = true;
        f.below = true;
        evaluate_visibility(&mut f, 0).unwrap();
        assert_eq!(
            f.trace,
            [
                "sun",
                "moon",
                "object",
                "place",
                "magnitude",
                "star-horizon"
            ]
        );
    }
    #[test]
    fn atmosphere_precedes_phase_and_partial_assignments_survive_late_catch() {
        let mut f = feed();
        f.failed = Some("sun-separation");
        let (_, words) = evaluate_visibility(&mut f, 0).unwrap();
        assert!(words[0].is_finite());
        assert_eq!(
            f.trace,
            [
                "sun",
                "moon",
                "object",
                "place",
                "magnitude",
                "model",
                "phase",
                "separation-setup",
                "moon-separation",
                "sun-separation"
            ]
        );
        // Independently compose the successful partial values and still-default
        // Sun separation. No expected astronomical word is fitted or stored.
        let model = VislimModel::new(
            &crate::heliacal::VislimAtmosphere::standard(0.0, 0.0, 2451545.0).unwrap(),
            crate::heliacal::HeliacalObserver::naked_eye(),
            true,
        );
        let expected = vis_limit_mag_words(
            &model,
            0,
            VislimBody::Planet(crate::constants::VENUS),
            &VislimGeometry {
                obj_alt_deg: 20.0,
                obj_az_deg: 80.0,
                obj_mag: 1.0,
                sun_alt_deg: 20.0,
                sun_az_deg: 80.0,
                moon_alt_deg: 20.0,
                moon_az_deg: 80.0,
                moon_phase: 0.75,
                sun_obj_angle_deg: 180.0,
                moon_obj_angle_deg: 60.0,
            },
        );
        assert_eq!(words.map(f64::to_bits), expected.1.map(f64::to_bits));
    }
    #[test]
    fn first_failed_sun_never_reads_moon_or_object() {
        let mut f = feed();
        f.failed = Some("sun");
        let error = evaluate_visibility(&mut f, 0).unwrap_err();
        assert!(matches!(
            error.primary_error(),
            Error::MathDomain { value: 2.0 }
        ));
        assert_eq!(f.trace, ["sun"]);
    }
    #[test]
    fn live_model_is_lazy_and_retained_across_later_geometry() {
        let mut f = feed();
        f.live_model = true;
        f.below = true;
        assert_eq!(evaluate_visibility(&mut f, 0).unwrap().0, -2.0);
        assert!(!f.trace.contains(&"live-model"));
        f.below = false;
        f.trace.clear();
        let (status, words) = evaluate_visibility(&mut f, 0).unwrap();
        assert_eq!((status, words[0], words[1]), (1.0, 0.75, 20.0));
        assert_eq!(
            &f.trace[5..],
            &[
                "live-model",
                "phase",
                "separation-setup",
                "moon-separation",
                "sun-separation"
            ]
        );
    }
}
