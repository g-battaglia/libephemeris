//! Independent WGS84 apparent-place acquisition for a numeric occulted body.
//! Parent serving primitives own images, clocks, Earth state and corrections.
use super::*;

/// One retained native observer for mutual planetary occultations. Site setup
/// and Earth lookup happen once; every observed place owns its independent
/// bundled UT1 Time and the actual apparent ICRS ray in AU.
pub(crate) struct NativePlanetObserver {
    context: EphemerisContext,
    itrs: [f64; 3],
    latlon: [[f64; 3]; 3],
}

/// Numeric body or the local walk's linearly propagated bare star direction.
/// The latter consumes neither parallax nor radial velocity.
#[derive(Clone, Copy)]
pub(crate) enum PlanetObservedTarget {
    Body(i32),
    Star { ra: f64, dec: f64 },
}

/// One observation's actual ICRS ray, AU distance and retained site rotations.
pub(crate) struct PlanetObservedPlace {
    pub(crate) ray: [f64; 3],
    pub(crate) distance: f64,
    rotation: [[f64; 3]; 3],
    latlon: [[f64; 3]; 3],
}
impl PlanetObservedPlace {
    /// Consume the real horizon projection only when altitude is requested.
    /// Return true altitude and south-westward azimuth, both in degrees.
    pub(crate) fn altaz(&self) -> Result<(f64, f64)> {
        let horizon = positions::reference_frames::matmul3(self.latlon, self.rotation)?;
        let (az, alt) = native_pheno_frame_angles(horizon, self.ray)?;
        Ok((alt, (az + 180.0).rem_euclid(360.0)))
    }
}

impl NativePlanetObserver {
    /// Bind an owned main image, east longitude/north latitude degrees and
    /// WGS84 height metres. Earth lookup precedes site construction.
    pub(crate) fn new(context: &EphemerisContext, site: [f64; 3]) -> Result<Self> {
        use positions::reference_frames::*;
        let provider =
            positions::LocalSpkStateProvider::new(context.spk_registry(), context.mode())?;
        provider.validate_source_alias(crate::constants::EARTH)?;
        let itrs = native_pheno_itrs_position(
            (site[0], site[1], site[2]),
            positions::SiteEllipsoid::Wgs84,
        )?;
        let phi = site[1] / 360.0 * FRAME_TAU;
        let zero = phi * 0.0;
        let one = zero + 1.0;
        let latitude = [
            [-scalar_sin(phi), zero, scalar_cos(phi)],
            [zero, one, zero],
            [scalar_cos(phi), zero, scalar_sin(phi)],
        ];
        let latlon = matmul3(latitude, rot_z(-(site[0] / 360.0 * FRAME_TAU)))?;
        Ok(Self {
            context: context.clone(),
            itrs,
            latlon,
        })
    }

    /// Numeric target lookup precedes Time. Earth/sample/rotation/light-time
    /// errors propagate. Distance is measured before spherical projection;
    /// true altitude and south-westward azimuth use the actual site frame.
    pub(crate) fn place(
        &self,
        jd: f64,
        target: PlanetObservedTarget,
    ) -> Result<PlanetObservedPlace> {
        use positions::reference_frames::*;
        let provider = positions::LocalSpkStateProvider::new(
            self.context.spk_registry(),
            self.context.mode(),
        )?;
        if let PlanetObservedTarget::Body(body) = target {
            provider.validate_source_alias(body)?;
        }
        let (epochs, tt) = VectorRotationEpochs::from_ut1(jd)?;
        self.observe_at(epochs, tt, target)
    }

    /// Consume an already acquired independent Time, then the numeric target
    /// lookup. Target-altitude and maximum geometry use this order; foreground
    /// altitude keeps its separate lookup-before-Time entry.
    pub(crate) fn place_after_time(
        &self,
        epochs: positions::reference_frames::VectorRotationEpochs,
        tt: f64,
        target: PlanetObservedTarget,
    ) -> Result<PlanetObservedPlace> {
        let provider = positions::LocalSpkStateProvider::new(
            self.context.spk_registry(),
            self.context.mode(),
        )?;
        if let PlanetObservedTarget::Body(body) = target {
            provider.validate_source_alias(body)?;
        }
        self.observe_at(epochs, tt, target)
    }

    /// Earth/site apparent observation at the actual selected split epoch.
    /// Target lookup belongs to the caller's distinct acquisition role.
    fn observe_at(
        &self,
        epochs: positions::reference_frames::VectorRotationEpochs,
        tt: f64,
        target: PlanetObservedTarget,
    ) -> Result<PlanetObservedPlace> {
        use positions::reference_frames::*;
        let provider = positions::LocalSpkStateProvider::new(
            self.context.spk_registry(),
            self.context.mode(),
        )?;
        let time = TwoPartEpoch {
            hi: epochs.whole,
            lo: epochs.tt_fraction,
        };
        let earth = provider.star_source_time_state(crate::constants::EARTH, time)?;
        let rotation =
            itrs_rotation_at_parts(epochs.whole, epochs.tt_fraction, epochs.ut1_fraction, tt)?;
        let (surface, site_velocity) = itrs_to_gcrs(rotation, self.itrs)?;
        let observer = std::array::from_fn(|i| earth.0[i] + surface[i]);
        let velocity = std::array::from_fn(|i| earth.1[i] + site_velocity[i]);
        let ray = match target {
            PlanetObservedTarget::Body(body) => provider.eclipse_source_observed_vector(
                body,
                time,
                observer,
                velocity,
                Some(surface),
            )?,
            PlanetObservedTarget::Star { ra, dec } => provider
                .arbitrary_star_apparent_vector(time, ra, dec, observer, velocity, surface)?,
        };
        let distance = native_pheno_vector_length(ray);
        Ok(PlanetObservedPlace {
            ray,
            distance,
            rotation,
            latlon: self.latlon,
        })
    }
}

/// Observe a classical numeric target from the actual height-metre WGS84 site
/// at original UT1 Julian Day. Acquire main image, independent bundled Time,
/// Earth alias/state, site rotation, target light-time/apparent corrections and
/// true ecliptic-of-date projection in that order. The target is selected only
/// at the observation stage; caller scalar model TT and flags are not consumed.
///
/// Returns longitude degrees in [0,360), latitude degrees and distance AU.
/// Surface position is present for real surface observer Earth-deflection semantics.
/// Missing aliases/ranges/directories, frame and light-time errors propagate;
/// no no-event sentinel, output oracle or global observer mutation is used.
pub(crate) fn native_occultation_topocentric_place(
    ctx: &EphemerisContext,
    jd_ut1: f64,
    site: [f64; 3],
    body: i32,
) -> Result<[f64; 3]> {
    use positions::reference_frames::*;
    if ctx.mode() != ProviderMode::NativeJpl {
        return Err(Error::ConfigurationMissing {
            missing_config: "native_occultation_context",
            suggestion: "select NativeJpl with the intended main-image configuration",
        });
    }
    let acquired = if ctx.spk_registry().is_empty() {
        Some(positions::main_image::acquire_main_image(ctx)?)
    } else {
        None
    };
    let context = acquired.as_ref().map_or(ctx, |owned| &owned.context);
    let (epochs, tt) = VectorRotationEpochs::from_ut1(jd_ut1)?;
    let provider = positions::LocalSpkStateProvider::new(context.spk_registry(), context.mode())?;
    provider.validate_source_alias(crate::constants::EARTH)?;
    let itrs =
        native_pheno_itrs_position((site[0], site[1], site[2]), positions::SiteEllipsoid::Wgs84)?;
    let time = TwoPartEpoch {
        hi: epochs.whole,
        lo: epochs.tt_fraction,
    };
    let earth = provider.star_source_time_state(crate::constants::EARTH, time)?;
    let rotation =
        itrs_rotation_at_parts(epochs.whole, epochs.tt_fraction, epochs.ut1_fraction, tt)?;
    let (surface, site_velocity) = itrs_to_gcrs(rotation, itrs)?;
    let observer = std::array::from_fn(|i| earth.0[i] + surface[i]);
    let velocity = std::array::from_fn(|i| earth.1[i] + site_velocity[i]);
    let ray =
        provider.eclipse_source_observed_vector(body, time, observer, velocity, Some(surface))?;
    let frame = native_pheno_ecliptic_frame(time)?;
    let (lon, lat) = native_pheno_frame_angles(frame, ray)?;
    let distance = native_pheno_vector_length(matvec3(frame, ray)?);
    Ok([lon, lat, distance])
}
