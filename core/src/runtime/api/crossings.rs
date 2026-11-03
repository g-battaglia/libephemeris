//! Longitude and lunar-node searches over the session's explicit provider.
//! Every search consumes the Session's public positions, retaining its
//! observer, sidereal frame, data ownership and input clock. Stored searches
//! require spherical ecliptic degrees.
use super::*;
pub use crate::events::SearchDirection;
#[path = "native_crossings.rs"]
pub(super) mod native;

/// Longitude crossing time in the input time scale, with any advisory.
#[derive(Debug, Clone, PartialEq)]
pub struct LongitudeCrossing {
    /// Julian Day, UT or TT according to the called entry.
    pub time: f64,
    /// Session configuration advisories, or empty.
    pub diagnostic: String,
}
/// Moon's ecliptic latitude-zero crossing with its actual coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct NodeCrossing {
    /// Crossing Julian Day in the input time scale.
    pub time: f64,
    /// Moon component zero; normally ecliptic longitude in degrees.
    /// Native node-search coordinate flags retain the selected output units.
    pub longitude: f64,
    /// Moon component one; normally ecliptic latitude in degrees.
    /// Native node-search coordinate flags retain the selected output units.
    pub latitude: f64,
    /// Session configuration advisories, or empty.
    pub diagnostic: String,
}
/// Restrict stored searches to spherical ecliptic components in degrees.
/// This guard runs before stored provider acquisition. Native searches instead
/// preserve the public position flags and normalize only longitude targets.
fn check_crossing_shape(flags: i32) -> Result<(), Error> {
    if flags & (FLG_RADIANS | FLG_XYZ | FLG_EQUATORIAL) != 0 {
        return Err(Error::UnsupportedFlagCombination {
            entry: "api crossing angular shape",
            flags,
        });
    }
    Ok(())
}
/// Convert a scalar longitude target to degrees and remove its radians flag.
/// Other bits remain available to the public position calculation. The node
/// search consumes its own coordinate components and does not use this helper.
fn cross_target(target: f64, flags: i32) -> (f64, i32) {
    let target = if flags & FLG_RADIANS != 0 {
        target.to_degrees()
    } else {
        target
    };
    (native::normalize_target(target), flags & !FLG_RADIANS)
}

/// Station motion catch: acquired range/ValueError constructors first, then
/// known numerical ValueError and calculation/range producer roles. Search
/// failures, loader I/O/corruption, lookup and unclassified errors do not fold.
pub(super) fn station_position_error(error: &Error) -> bool {
    native::input_failure(error)
        || matches!(
            error.primary_error(),
            Error::EphemerisRange { .. }
                | Error::LebOutOfRange { .. }
                | Error::SpkOutOfRange { .. }
                | Error::EpochOutOfSegment { .. }
                | Error::PolarCircle { .. }
                | Error::ConvergenceFailed { .. }
        )
}
/// Public Error/RuntimeError station recovery. An acquired closed public class
/// determines its hierarchy. Without one, only known public constructor roles
/// are caught; raw I/O, malformed providers, model lookup, builtin arithmetic
/// and unclassified errors remain terminal. Unknown backend class correspondence
/// requires acquisition at its actual producer, not inference from Display.
pub(super) fn station_recoverable(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::EphemerisRange { .. }
            | Error::LebOutOfRange { .. }
            | Error::SpkOutOfRange { .. }
            | Error::EpochOutOfSegment { .. }
            | Error::ConvergenceFailed { .. }
            | Error::ConfigurationMissing { .. }
            | Error::UnknownBody { .. }
            | Error::LebBodyNotFound { .. }
            | Error::NetworkSealed { .. }
            | Error::NoStationInInterval { .. }
    )
}

impl Session {
    /// Next forward longitude crossing for a selected body at UT1.
    ///
    /// Native positions retain the scalar UT1 query and independent position
    /// clock at every sample; radians targets are converted to degrees before
    /// the search. The adaptive solver uses longitude and its daily speed,
    /// scans retrograde/station branches and refines the first crossing.
    /// Stored providers require degree spherical ecliptic flags. Native flags
    /// select the public coordinate components consumed by the search, so
    /// ordinary ecliptic flags should be used for physical longitude events.
    /// Entry and Newton calculation/range/value failures become `PublicCrossing`
    /// with stage text and their owned cause. Raw scan failures propagate. Only
    /// a Newton public range failure initiates range recovery, whose new range
    /// failure is retained. Other lookup, configuration and I/O failures retain
    /// their types. The captured initial iteration cap does not grow.
    pub fn cross_ut(
        &self,
        body: i32,
        target: f64,
        start: f64,
        flags: i32,
    ) -> Result<LongitudeCrossing, Error> {
        if self.ctx.mode() != ProviderMode::NativeJpl {
            check_crossing_shape(flags)?;
        }
        let (target, flags) = cross_target(target, flags);
        let mut sample = |jd, sample_flags| {
            self.calc_ut(jd, body, sample_flags)
                .map(|position| position.values)
        };
        let time = native::planet(&mut sample, body, target, start, flags)?;
        Ok(LongitudeCrossing {
            time,
            diagnostic: self.diagnostics.join("; "),
        })
    }

    /// Acquire a public six-component position in the search's input scale.
    /// The orchestration supplies exact sample flags, including SPEED only
    /// where consumed. UT1 is preserved on the UT branch; TT is passed unchanged
    /// on the TT branch. This acquisition helper does not wrap failures.
    fn crossing_position(
        &self,
        jd: f64,
        body: i32,
        flags: i32,
        is_ut: bool,
    ) -> Result<Position, Error> {
        if is_ut {
            self.calc_ut(jd, body, flags)
        } else {
            self.calc(jd, body, flags)
        }
    }
    /// Next Sun longitude crossing at UT, using the existing Newton/search
    /// contract. Native radians targets are converted to degrees; public
    /// scalar UT1 position samples supply longitude and daily speed. Stored
    /// RADIANS/XYZ/EQUATORIAL shapes remain unsupported. Native flags select
    /// the actual position components, with ordinary ecliptic flags appropriate
    /// for physical longitude events. Entry and Newton calculation/range/value
    /// failures become `PublicCrossing` with their stage and owned cause;
    /// unrelated failures propagate. Exhaustion performs no extra sample.
    pub fn solcross_ut(
        &self,
        target: f64,
        start: f64,
        flags: i32,
    ) -> Result<LongitudeCrossing, Error> {
        self.solcross_ut_direction(target, start, flags, SearchDirection::Forward)
    }
    /// Next Sun crossing at TT; TT is passed unchanged to the stored provider.
    pub fn solcross(
        &self,
        target: f64,
        start: f64,
        flags: i32,
    ) -> Result<LongitudeCrossing, Error> {
        self.solcross_direction(target, start, flags, SearchDirection::Forward)
    }
    /// Next Moon longitude crossing at UT, with the same units/error policy.
    pub fn mooncross_ut(
        &self,
        target: f64,
        start: f64,
        flags: i32,
    ) -> Result<LongitudeCrossing, Error> {
        self.mooncross_ut_direction(target, start, flags, SearchDirection::Forward)
    }
    /// Next Moon crossing at TT, without an approximate UT round trip.
    pub fn mooncross(
        &self,
        target: f64,
        start: f64,
        flags: i32,
    ) -> Result<LongitudeCrossing, Error> {
        self.mooncross_direction(target, start, flags, SearchDirection::Forward)
    }
    /// Sun longitude crossing in either direction at UT1. Native acquisitions
    /// retain each scalar UT1 and its independent position clocks. Backward
    /// searches exclude an at-start crossing; forward searches include it.
    pub fn solcross_ut_direction(
        &self,
        target: f64,
        start: f64,
        flags: i32,
        direction: SearchDirection,
    ) -> Result<LongitudeCrossing, Error> {
        self.longitude_crossing(SUN, target, start, flags, direction, true)
    }
    /// Sun component-zero crossing in either direction at TT.
    /// The input and every position sample remain scalar TT. Targets are degrees
    /// or radians according to RADIANS; the solver works in degrees after
    /// stripping that bit. Backward searches seed the preceding cycle at an
    /// at-start crossing. Fifty Newton samples use the component-three daily
    /// rate and a one-milliarcsecond residual tolerance. Selected position
    /// failures become `PublicCrossing` with an owned cause; unrelated failures
    /// remain unchanged. Stored searches retain their angular-shape guards.
    pub fn solcross_direction(
        &self,
        target: f64,
        start: f64,
        flags: i32,
        direction: SearchDirection,
    ) -> Result<LongitudeCrossing, Error> {
        self.longitude_crossing(SUN, target, start, flags, direction, false)
    }
    /// Moon component-zero crossing in either direction at TT, with scalar TT
    /// retained at every sample. Units, direction and error ownership follow
    /// `solcross_direction`; the Moon uses a 13.176 deg/day seed and 0.1 deg/day
    /// Newton-rate floor, with a 31-day divergence envelope about its seed.
    pub fn mooncross_direction(
        &self,
        target: f64,
        start: f64,
        flags: i32,
        direction: SearchDirection,
    ) -> Result<LongitudeCrossing, Error> {
        self.longitude_crossing(MOON, target, start, flags, direction, false)
    }
    /// Moon longitude crossing in either direction at UT1, retaining scalar UT1 samples.
    pub fn mooncross_ut_direction(
        &self,
        target: f64,
        start: f64,
        flags: i32,
        direction: SearchDirection,
    ) -> Result<LongitudeCrossing, Error> {
        self.longitude_crossing(MOON, target, start, flags, direction, true)
    }
    /// Next Moon latitude-zero crossing at UT, with its actual coordinates.
    /// Native UT uses independent bundled UT1-to-TT and TT-to-UT1 conversions
    /// around the TT search. Native coordinate flags select the retained
    /// components and their rate; ordinary ecliptic flags identify lunar nodes.
    /// Stored serving retains its configured-clock angular domain.
    pub fn mooncross_node_ut(&self, start: f64, flags: i32) -> Result<NodeCrossing, Error> {
        self.mooncross_node_ut_direction(start, flags, SearchDirection::Forward)
    }
    /// Next Moon latitude-zero crossing at TT, with the same coordinate layout.
    pub fn mooncross_node(&self, start: f64, flags: i32) -> Result<NodeCrossing, Error> {
        self.mooncross_node_direction(start, flags, SearchDirection::Forward)
    }
    /// Moon latitude-zero crossing in either direction at UT1. Independent bundled
    /// UT1/TT conversions surround the scalar TT node search.
    pub fn mooncross_node_ut_direction(
        &self,
        start: f64,
        flags: i32,
        direction: SearchDirection,
    ) -> Result<NodeCrossing, Error> {
        self.node_crossing(start, flags, direction, true)
    }
    /// Moon component-one zero in either direction at TT. The returned TT,
    /// component zero and component one are from the accepted sample. Native
    /// samples add SPEED and use component four directly as the daily latitude
    /// derivative. Two-day scanning and sixty midpoint refinements select the
    /// first root; fifty Newton acquisitions enforce the requested direction.
    /// At-start roots are excluded by a day-scale epsilon. Stored searches keep
    /// their distinct configured-clock and latitude-stencil contract. Selected
    /// position failures own their public crossing stage and original cause.
    pub fn mooncross_node_direction(
        &self,
        start: f64,
        flags: i32,
        direction: SearchDirection,
    ) -> Result<NodeCrossing, Error> {
        self.node_crossing(start, flags, direction, false)
    }
    /// Search from the scalar TT acquired by an independent Time object.
    ///
    /// Native positions retain that scalar's Time provenance in their public
    /// range/error constructors. The returned epoch is TT: the owner of the
    /// original Time performs any inverse conversion after a successful search.
    /// No configured Delta-T, UT1-to-TT refinement or inverse clock is consumed
    /// here. Flags, direction, sampling and solver limits are unchanged. Stored
    /// providers retain the ordinary TT node-search contract.
    pub fn mooncross_node_from_time_direction(
        &self,
        start_tt: f64,
        flags: i32,
        direction: SearchDirection,
    ) -> Result<NodeCrossing, Error> {
        if self.ctx.mode() == ProviderMode::NativeJpl {
            self.native_node_from_tt(start_tt, flags, direction, true)
        } else {
            self.node_crossing(start_tt, flags, direction, false)
        }
    }
    /// Heliocentric crossing at UT; negative `dir` searches backward, otherwise
    /// forward. Native targets follow the radians flag and public scalar UT1
    /// sampling; stored non-angular/radian shapes refuse. Mean-rate seeds and
    /// bracket refinement retain nearest-past/first-forward semantics. Selected
    /// entry and Newton failures own a public crossing stage and cause.
    pub fn helio_cross_ut(
        &self,
        body: i32,
        target: f64,
        start: f64,
        flags: i32,
        dir: i32,
    ) -> Result<LongitudeCrossing, Error> {
        self.helio_crossing(body, target, start, flags, dir, true)
    }
    /// Heliocentric crossing at TT; no model-clock conversion of input samples.
    pub fn helio_cross(
        &self,
        body: i32,
        target: f64,
        start: f64,
        flags: i32,
        dir: i32,
    ) -> Result<LongitudeCrossing, Error> {
        self.helio_crossing(body, target, start, flags, dir, false)
    }
    /// Search Sun or Moon component zero using its daily component-three rate.
    /// Native UT samples retain their original UT1 clock acquisition; native TT
    /// samples remain scalar TT. Radians targets become degrees before native
    /// sampling. Stored acquisition keeps its configured Delta-T/identity-TT
    /// serving and angular-shape guards. No coordinate result is cached or
    /// substituted. Selected entry/iteration failures become public crossing
    /// errors with owned causes; unlisted failures propagate unchanged.
    fn longitude_crossing(
        &self,
        body: i32,
        target: f64,
        start: f64,
        flags: i32,
        direction: SearchDirection,
        is_ut: bool,
    ) -> Result<LongitudeCrossing, Error> {
        if self.ctx.mode() != ProviderMode::NativeJpl {
            check_crossing_shape(flags)?;
        }
        let (target, flags) = cross_target(target, flags);
        let mut sample = |jd, sample_flags| {
            self.crossing_position(jd, body, sample_flags, is_ut)
                .map(|position| position.values)
        };
        let time = native::luminary(&mut sample, body, target, start, flags, direction)?;
        Ok(LongitudeCrossing {
            time,
            diagnostic: self.diagnostics.join("; "),
        })
    }
    /// Search the Moon's component-one zero and retain components zero and one
    /// at the resulting sample. Native UT converts once to independent bundled
    /// TT, searches through scalar public TT positions and converts the final
    /// epoch back to bundled UT1. The Newton derivative is the actual public
    /// component-four rate; output coordinate flags remain intact. Stored
    /// acquisition samples the public configured-clock positions directly. Domain, data and solver failures remain typed errors.
    fn node_crossing(
        &self,
        start: f64,
        flags: i32,
        direction: SearchDirection,
        is_ut: bool,
    ) -> Result<NodeCrossing, Error> {
        if self.ctx.mode() == ProviderMode::NativeJpl {
            let start_tt = if is_ut {
                crate::smh2016::source_time_scalar_tt_from_ut1(start)
            } else {
                start
            };
            let mut hit = self.native_node_from_tt(start_tt, flags, direction, is_ut)?;
            if is_ut {
                hit.time = crate::time::default_timescale_ut1_from_tt(hit.time)?;
            }
            return Ok(hit);
        }
        check_crossing_shape(flags)?;
        let mut sample = |jd, sample_flags| {
            self.crossing_position(jd, MOON, sample_flags, is_ut)
                .map(|position| position.values)
        };
        let hit = native::node(&mut sample, start, flags, direction)?;
        Ok(NodeCrossing {
            time: hit.time,
            longitude: hit.longitude,
            latitude: hit.latitude,
            diagnostic: self.diagnostics.join("; "),
        })
    }
    /// Native TT node loop with its consumed scalar-constructor origin.
    /// The sampling callback is the only distinction between a public float
    /// and a scalar acquired from Time. Search arithmetic, coordinate words,
    /// flags and error wrapping share the same native implementation.
    fn native_node_from_tt(
        &self,
        start_tt: f64,
        flags: i32,
        direction: SearchDirection,
        from_time: bool,
    ) -> Result<NodeCrossing, Error> {
        let mut sample = |jd, sample_flags| {
            let result = if from_time {
                self.calc_from_native_time_scalar(jd, MOON, sample_flags)
            } else {
                self.calc(jd, MOON, sample_flags)
            };
            result.map(|position| position.values)
        };
        let hit = native::node(&mut sample, start_tt, flags, direction)?;
        Ok(NodeCrossing {
            time: hit.time,
            longitude: hit.longitude,
            latitude: hit.latitude,
            diagnostic: self.diagnostics.join("; "),
        })
    }
    /// Search a heliocentric component-zero crossing in the input time scale.
    /// Native sampling adds HELCTR and SPEED to public scalar UT1/TT positions
    /// after converting radians targets to degrees. Direction is backward for
    /// negative `dir`, forward otherwise. The existing eccentric-orbit seed,
    /// nearest-past/first-forward repair and numerical stopping rules remain
    /// in Native orchestration. Stored serving keeps its strict angular-shape domain;
    /// selected entry/iteration calculation/range/value failures retain their
    /// public crossing stage and cause. Bracket recovery catches only its
    /// designated public Error branch; other failures remain visible.
    fn helio_crossing(
        &self,
        body: i32,
        target: f64,
        start: f64,
        flags: i32,
        dir: i32,
        is_ut: bool,
    ) -> Result<LongitudeCrossing, Error> {
        if self.ctx.mode() != ProviderMode::NativeJpl {
            check_crossing_shape(flags)?;
        }
        let direction = if dir < 0 {
            SearchDirection::Backward
        } else {
            SearchDirection::Forward
        };
        let (target, flags) = cross_target(target, flags);
        let mut sample = |jd, sample_flags| {
            self.crossing_position(jd, body, sample_flags, is_ut)
                .map(|position| position.values)
        };
        let time = native::helio(&mut sample, body, target, start, flags, direction)?;
        Ok(LongitudeCrossing {
            time,
            diagnostic: self.diagnostics.join("; "),
        })
    }
}

#[cfg(test)]
mod native_public_crossing_routes {
    use super::*;
    use std::fs::OpenOptions;
    use std::path::PathBuf;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static COUNTER: AtomicUsize = AtomicUsize::new(0);

    struct OwnedRootFile(PathBuf);

    impl OwnedRootFile {
        fn new() -> Self {
            let path = std::env::temp_dir().join(format!(
                "libephemeris-native-crossing-root-{}-{}",
                std::process::id(),
                COUNTER.fetch_add(1, Ordering::Relaxed)
            ));
            OpenOptions::new()
                .create_new(true)
                .write(true)
                .open(&path)
                .unwrap();
            Self(path)
        }
    }

    impl Drop for OwnedRootFile {
        fn drop(&mut self) {
            std::fs::remove_file(&self.0).unwrap();
        }
    }

    fn assert_actual_main_root_failure<T>(result: Result<T, Error>, root: &OwnedRootFile) {
        match result {
            Err(Error::EphemerisLoaderIo {
                path, kind, errno, ..
            }) => {
                assert_eq!(path, root.0.to_string_lossy());
                assert_eq!(kind, std::io::ErrorKind::AlreadyExists);
                assert!(errno.is_some());
            }
            Err(other) => panic!("unexpected error {other:?}"),
            Ok(_) => panic!("configured regular-file main root was accepted"),
        }
    }

    #[test]
    fn native_longitude_methods_enter_actual_owned_main_loader() {
        let root = OwnedRootFile::new();
        let ctx = EphemerisContext::new()
            .with_mode(ProviderMode::NativeJpl)
            .with_network_policy(Some(NetworkPolicy::Sealed))
            .with_jpl_data_dir(Some(root.0.clone()));
        let session = SessionBuilder::from_context(ctx).build().unwrap();
        let flags = FLG_SWIEPH | FLG_RADIANS;
        assert_actual_main_root_failure(session.solcross_ut(0.0, 2451545.125, flags), &root);
        assert_actual_main_root_failure(session.solcross(0.0, 2451545.125, flags), &root);
        assert_actual_main_root_failure(
            session.solcross_direction(0.0, 2451545.125, flags, SearchDirection::Backward),
            &root,
        );
        assert_actual_main_root_failure(session.mooncross_ut(0.0, 2451545.125, flags), &root);
        assert_actual_main_root_failure(session.mooncross(0.0, 2451545.125, flags), &root);
        assert_actual_main_root_failure(
            session.mooncross_direction(0.0, 2451545.125, flags, SearchDirection::Backward),
            &root,
        );
        assert_actual_main_root_failure(session.cross_ut(MARS, 0.0, 2451545.125, flags), &root);
        assert_actual_main_root_failure(
            session.helio_cross_ut(MARS, 0.0, 2451545.125, flags, -1),
            &root,
        );
        assert_actual_main_root_failure(
            session.helio_cross(MARS, 0.0, 2451545.125, flags, 1),
            &root,
        );
    }

    #[test]
    fn native_node_methods_enter_public_tt_main_loader() {
        let root = OwnedRootFile::new();
        let ctx = EphemerisContext::new()
            .with_mode(ProviderMode::NativeJpl)
            .with_network_policy(Some(NetworkPolicy::Sealed))
            .with_user_delta_t(Some(100.0))
            .with_jpl_data_dir(Some(root.0.clone()));
        let session = SessionBuilder::from_context(ctx).build().unwrap();
        assert_actual_main_root_failure(session.mooncross_node_ut(2451545.125, FLG_SWIEPH), &root);
        assert_actual_main_root_failure(session.mooncross_node(2451545.125, FLG_SWIEPH), &root);
        assert_actual_main_root_failure(
            session.mooncross_node_from_time_direction(
                2451545.125,
                FLG_SWIEPH,
                SearchDirection::Backward,
            ),
            &root,
        );
        assert_actual_main_root_failure(
            session.mooncross_node_direction(2451545.125, FLG_SWIEPH, SearchDirection::Backward),
            &root,
        );
    }

    #[test]
    fn stored_provider_angular_guards_remain_before_acquisition() {
        let session = SessionBuilder::new().build().unwrap();
        for flag in [FLG_RADIANS, FLG_XYZ, FLG_EQUATORIAL] {
            assert!(matches!(
                session.solcross_ut(0.0, 2451545.125, flag),
                Err(Error::UnsupportedFlagCombination { .. })
            ));
            assert!(matches!(
                session.cross_ut(MARS, 0.0, 2451545.125, flag),
                Err(Error::UnsupportedFlagCombination { .. })
            ));
            assert!(matches!(
                session.mooncross_node_ut(2451545.125, flag),
                Err(Error::UnsupportedFlagCombination { .. })
            ));
        }
    }
}
