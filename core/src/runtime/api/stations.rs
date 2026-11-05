//! Public station routes over actual scalar-UT position calls. Each operation
//! retains its Session and owned provider inputs until all samples complete.
use super::*;
pub use crate::events::{StationKind, StationMotionKind, StationRequest};
#[path = "native_stations.rs"]
mod native;

/// Nearest-station information, including the public partial-success branch.
/// A partial result owns the first caught search error without ambient state.
#[derive(Debug, Clone, PartialEq)]
pub enum StationInformation {
    /// Both current motion and a selected past/future station are available.
    Complete(crate::events::StationInfo),
    /// Both bounded searches failed with caught public errors. The requested
    /// epoch's motion remains valid; four station fields are absent.
    Unavailable {
        /// Current component-three sign at the original UT1 sample.
        currently_retrograde: bool,
        /// Current component-three rate, normally degrees/day.
        velocity: f64,
        /// First forward/backward caught failure; its owned constructor is
        /// retained when acquired. Missing backend presentation is not inferred.
        error: Error,
    },
}

/// Closed station request role. Flags remain original public position flags.
/// Dates and steps are Julian days; no checked solver domain is added here.
#[derive(Debug, Clone, Copy)]
pub enum StationQuery {
    /// Current retrograde predicate; Sun/Moon consume no position.
    Retrograde {
        /// Public body identifier; Sun and Moon return false without sampling.
        body: i32,
        /// Epoch as a UT1 Julian day.
        jd_ut: f64,
        /// Public position flags; SPEED is added when acquiring the rate.
        flags: i32,
    },
    /// Current motion type with a consumed one-day trend near stations.
    Motion {
        /// Public body identifier; Sun and Moon are classified as direct.
        body: i32,
        /// UT1 Julian day of the current-motion sample.
        jd_ut: f64,
        /// Public coordinate flags retained for the current and trend samples.
        flags: i32,
    },
    /// First forward station matching a normalized request.
    Find {
        /// Public body identifier; Sun and Moon have no station search.
        body: i32,
        /// Starting UT1 Julian day for the forward search.
        jd_ut: f64,
        /// Requested SR/SD transition, or either kind.
        request: StationRequest,
        /// Public position flags used throughout the search; SPEED is added.
        flags: i32,
    },
    /// Next complete SR/SD period, skipping an already active period.
    Next {
        /// Public body identifier; Sun and Moon are rejected.
        body: i32,
        /// Starting UT1 Julian day; an already active retrograde period is skipped.
        jd_ut: f64,
        /// Public coordinate flags retained for the constituent station searches.
        flags: i32,
    },
    /// Central difference of public components with caller flags unchanged.
    Velocity {
        /// Public body identifier sampled at both difference endpoints.
        body: i32,
        /// Central UT1 Julian day.
        jd_ut: f64,
        /// Original position flags, forwarded without adding SPEED.
        flags: i32,
        /// Difference half-step in days. Negative steps are accepted; zero
        /// refuses at division after both position acquisitions succeed.
        dt: f64,
    },
    /// Nearest past/future station or partial motion information.
    Information {
        /// Public body identifier; Sun and Moon are rejected.
        body: i32,
        /// Reference UT1 Julian day for current motion and nearest-station distance.
        jd_ut: f64,
        /// Public coordinate flags used for motion and the selected station position.
        flags: i32,
    },
}

/// Owned query result; no borrowed provider or raw pointer.
#[derive(Debug, Clone, PartialEq)]
pub enum StationOutput {
    /// Public negative-rate predicate. Sun/Moon identities and handled position
    /// failures yield false without establishing a sampled rate.
    Retrograde(bool),
    /// Direct, retrograde or near-station classification. Sun/Moon identities
    /// and handled position failures yield direct motion.
    Motion(StationMotionKind),
    /// Station UT1 Julian day followed by its SR/SD transition kind.
    Station(f64, StationKind),
    /// Start (SR) and end (SD) of the next complete period, as UT1 Julian days.
    Period(f64, f64),
    /// Differences of public components zero, one and two, divided by twice
    /// the step in days. Component zero receives one +/-360 fold, including
    /// under nonangular coordinate flags. Ordinary spherical flags give
    /// degrees/day, degrees/day and AU/day; other flags retain their components.
    Velocity(f64, f64, f64),
    /// Nearest-station details or current-motion information with caught search errors.
    Information(StationInformation),
}

/// Optional binding-owned configuration acquisition, never astronomical data.
/// Its synchronous borrow ends with the query; it preserves provider and clocks.
type PositionResolver<'a> = dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext, Error> + 'a;

/// Normalize accepted ASCII station names; all other strings select any.
/// The binding performs any user-object `.upper` operation before this
/// bounded string role. No unknown-spelling refusal is added to the SDK route.
pub fn station_request(value: &str) -> StationRequest {
    if value.eq_ignore_ascii_case("SR") || value.eq_ignore_ascii_case("RETROGRADE") {
        StationRequest::Retrograde
    } else if value.eq_ignore_ascii_case("SD") || value.eq_ignore_ascii_case("DIRECT") {
        StationRequest::Direct
    } else {
        StationRequest::Any
    }
}

/// Source-defined Sun/Moon identity before unused configuration acquisition.
/// Velocity still consumes both public positions. Search and summary failures
/// own their distinct builtin ValueError constructors; all other queries return
/// None and must use their genuine Session route.
pub fn station_identity(query: StationQuery) -> Option<Result<StationOutput, Error>> {
    match query {
        StationQuery::Retrograde {
            body: SUN | MOON, ..
        } => Some(Ok(StationOutput::Retrograde(false))),
        StationQuery::Motion {
            body: SUN | MOON, ..
        } => Some(Ok(StationOutput::Motion(StationMotionKind::Direct))),
        StationQuery::Find {
            body: body @ (SUN | MOON),
            ..
        }
        | StationQuery::Next {
            body: body @ (SUN | MOON),
            ..
        } => Some(native::check_body(body, false).map(|()| unreachable!("identity refuses"))),
        StationQuery::Information {
            body: body @ (SUN | MOON),
            ..
        } => Some(native::check_body(body, true).map(|()| unreachable!("identity refuses"))),
        _ => None,
    }
}

impl Session {
    /// True when the public component-three rate is negative. Sun/Moon return
    /// false without consuming date, flags, clocks or providers. Selected
    /// calculation/range/ValueError failures fold to false; other failures
    /// retain their types. SPEED is added at the actual scalar-UT sample.
    pub fn is_retrograde(&self, body: i32, jd_ut: f64, flags: i32) -> Result<bool, Error> {
        match self.station_query(StationQuery::Retrograde { body, jd_ut, flags })? {
            StationOutput::Retrograde(value) => Ok(value),
            _ => unreachable!("closed query role"),
        }
    }
    /// Current direct/retrograde/near-station motion type. The near threshold
    /// is 0.01 component units/day; minus-one then plus-one-day rates choose a
    /// decreasing/increasing trend. Sun/Moon are direct without acquisition.
    pub fn get_station_type(
        &self,
        body: i32,
        jd_ut: f64,
        flags: i32,
    ) -> Result<StationMotionKind, Error> {
        match self.station_query(StationQuery::Motion { body, jd_ut, flags })? {
            StationOutput::Motion(value) => Ok(value),
            _ => unreachable!("closed query role"),
        }
    }
    /// First requested forward station at UT1. Initial rate acquisition
    /// precedes four half-period attempts, fifty-interval scans and Brent
    /// refinement (rate tolerance 1e-6, time guard 1e-8 days, cap100). Search
    /// recovery covers scan/refinement/classification, not builtin arithmetic
    /// or ValueError. SR is positive-before/negative-after; all other trends
    /// are SD. The date and independent public position clocks are retained at
    /// every sample. Exhaustion owns a public base Error constructor. Sun/Moon
    /// own a ValueError before provider/time lookup. No physical certification
    /// follows from the solver tolerance or an unusual coordinate flag.
    pub fn find_station_ut(
        &self,
        body: i32,
        jd_ut: f64,
        request: StationRequest,
        flags: i32,
    ) -> Result<(f64, StationKind), Error> {
        match self.station_query(StationQuery::Find {
            body,
            jd_ut,
            request,
            flags,
        })? {
            StationOutput::Station(time, kind) => Ok((time, kind)),
            _ => unreachable!("closed query role"),
        }
    }
    /// Next complete retrograde period as SR/SD UT1 Julian days. Consume the
    /// current retrograde predicate first; an ongoing period is skipped by
    /// finding its SD before the following SR and SD. A subsequent search
    /// begins one day after the preceding root. Errors retain their stages.
    pub fn next_retrograde_ut(
        &self,
        body: i32,
        jd_ut: f64,
        flags: i32,
    ) -> Result<(f64, f64), Error> {
        match self.station_query(StationQuery::Next { body, jd_ut, flags })? {
            StationOutput::Period(sr, sd) => Ok((sr, sd)),
            _ => unreachable!("closed query role"),
        }
    }
    /// Difference public components at jd-dt then jd+dt and divide by 2*dt.
    /// Original flags are forwarded without SPEED. Component zero receives
    /// one degree-turn wrap; rates normally have deg/day,deg/day,AU/day units.
    /// Negative steps are accepted; zero division refuses after both position
    /// acquisitions. No checked positive-step guard is inherited from sealed
    /// low-level kernels. Coordinate-shaped inputs retain their actual values.
    pub fn calc_velocity_at_station(
        &self,
        body: i32,
        jd_ut: f64,
        flags: i32,
        dt: f64,
    ) -> Result<(f64, f64, f64), Error> {
        match self.station_query(StationQuery::Velocity {
            body,
            jd_ut,
            flags,
            dt,
        })? {
            StationOutput::Velocity(a, b, c) => Ok((a, b, c)),
            _ => unreachable!("closed query role"),
        }
    }
    /// Nearest station information. Current motion precedes forward then
    /// backward searches; a time-distance tie retains the forward candidate.
    /// The selected longitude uses original flags. Double caught search
    /// failure is a partial success with current motion and the first error.
    /// Sun/Moon own a distinct ValueError. Raw backend graph gaps stay explicit.
    pub fn get_station_info(
        &self,
        body: i32,
        jd_ut: f64,
        flags: i32,
    ) -> Result<StationInformation, Error> {
        match self.station_query(StationQuery::Information { body, jd_ut, flags })? {
            StationOutput::Information(value) => Ok(value),
            _ => unreachable!("closed query role"),
        }
    }
    /// Execute one closed station operation using public `calc_ut` at every
    /// scalar sample. Stored modes retain explicit LEB refusal; NativeJpl uses
    /// the current main-image, split samples, whole-target/range and per-segment
    /// AU contracts. No separate vector provider or clock is assembled.
    pub fn station_query(&self, query: StationQuery) -> Result<StationOutput, Error> {
        self.station_query_inner(query, None)
    }
    /// Binding configuration counterpart of [`Self::station_query`]. A resolver
    /// is consumed once at the first standard-body Native position with an empty
    /// registry, after original nonfinite/priority-resolved flag refusals. Other body IDs, supplied
    /// registries and provider modes skip it. The resolved immutable Session is
    /// retained through this operation; public `calc_ut` still owns each image
    /// acquisition and its independent position/model/frame clocks. The core
    /// reads no process environment. Resolver failures follow the station's
    /// source catch stage. No callback executes under a cache mutex, and no
    /// borrowed configuration outlives this call. Source-wide global mutable
    /// cache sharing is not introduced by this call-local configuration role.
    pub fn station_query_with_position_context(
        &self,
        query: StationQuery,
        resolve: &PositionResolver<'_>,
    ) -> Result<StationOutput, Error> {
        self.station_query_inner(query, Some(resolve))
    }
    /// Build one call-local position closure. A supplied registry is borrowed
    /// directly; an optional empty-registry configuration is retained once.
    /// Early position refusals delegate to the actual public entry before any
    /// unused configuration callback. The shared public flag validator resolves
    /// TOPOCTR center and JPLEPH/SWIEPH ephemeris priority before refusal. Raw
    /// overlapping bits alone never suppress a valid request's resolver.
    /// Success/errors are never cached as words.
    fn station_query_inner(
        &self,
        query: StationQuery,
        resolve: Option<&PositionResolver<'_>>,
    ) -> Result<StationOutput, Error> {
        let body = match query {
            StationQuery::Retrograde { body, .. }
            | StationQuery::Motion { body, .. }
            | StationQuery::Find { body, .. }
            | StationQuery::Next { body, .. }
            | StationQuery::Velocity { body, .. }
            | StationQuery::Information { body, .. } => body,
        };
        let mut selected = None;
        let mut sample = |jd: f64, flags: i32| {
            if !jd.is_finite()
                || crate::positions::validate_barycentric_moseph(flags, "calc_ut").is_err()
            {
                return self.calc_ut(jd, body, flags).map(|p| p.values);
            }
            if let Some(resolve) = resolve {
                if self.ctx.mode() == ProviderMode::NativeJpl
                    && matches!(body, SUN..=PLUTO | EARTH)
                    && self.ctx.spk_registry().is_empty()
                    && selected.is_none()
                {
                    let context = resolve(&self.ctx, body)?;
                    if context.mode() != ProviderMode::NativeJpl {
                        return Err(Error::ConfigurationMissing {
                            missing_config: "native_station_provider",
                            suggestion: "The position resolver must preserve NativeJpl mode",
                        });
                    }
                    selected = Some(Session::from_context(context));
                }
            }
            selected
                .as_ref()
                .unwrap_or(self)
                .calc_ut(jd, body, flags)
                .map(|p| p.values)
        };
        match query {
            StationQuery::Retrograde { jd_ut, flags, .. } => {
                native::retrograde(&mut sample, body, jd_ut, flags).map(StationOutput::Retrograde)
            }
            StationQuery::Motion { jd_ut, flags, .. } => {
                native::motion(&mut sample, body, jd_ut, flags).map(StationOutput::Motion)
            }
            StationQuery::Find {
                jd_ut,
                request,
                flags,
                ..
            } => native::find(&mut sample, body, jd_ut, request, flags)
                .map(|(a, b)| StationOutput::Station(a, b)),
            StationQuery::Next { jd_ut, flags, .. } => {
                native::next_period(&mut sample, body, jd_ut, flags)
                    .map(|(a, b)| StationOutput::Period(a, b))
            }
            StationQuery::Velocity {
                jd_ut, flags, dt, ..
            } => native::velocity(&mut sample, jd_ut, flags, dt)
                .map(|(a, b, c)| StationOutput::Velocity(a, b, c)),
            StationQuery::Information { jd_ut, flags, .. } => {
                native::info(&mut sample, body, jd_ut, flags).map(StationOutput::Information)
            }
        }
    }
}

#[cfg(test)]
#[path = "native_station_tests.rs"]
mod tests;

#[cfg(test)]
mod configuration_priority_tests {
    use super::*;
    use std::cell::Cell;

    fn native() -> Session {
        Session::from_context(
            EphemerisContext::default()
                .with_mode(ProviderMode::NativeJpl)
                .with_observer(crate::context::Observer::new(0.0, 0.0, 0.0).unwrap()),
        )
    }
    #[test]
    fn accepted_priority_flags_consume_native_configuration_once() {
        for flags in [
            FLG_BARYCTR | FLG_MOSEPH | FLG_SWIEPH,
            FLG_BARYCTR | FLG_MOSEPH | FLG_JPLEPH,
            FLG_BARYCTR | FLG_MOSEPH | FLG_TOPOCTR,
        ] {
            let called = Cell::new(0);
            let resolver = |_: &EphemerisContext, body| {
                assert_eq!(body, MARS);
                called.set(called.get() + 1);
                Err(Error::ConfigurationMissing {
                    missing_config: "station_resolver_witness",
                    suggestion: "owned test marker",
                })
            };
            let error = native()
                .station_query_with_position_context(
                    StationQuery::Velocity {
                        body: MARS,
                        jd_ut: 2451545.0,
                        flags,
                        dt: 1e-4,
                    },
                    &resolver,
                )
                .unwrap_err();
            assert_eq!(called.get(), 1);
            assert!(matches!(
                error.primary_error(),
                Error::ConfigurationMissing {
                    missing_config: "station_resolver_witness",
                    ..
                }
            ));
        }
    }
    #[test]
    fn actual_early_public_refusals_do_not_consume_configuration() {
        for (jd, flags, unsupported) in [
            (2451545.0, FLG_BARYCTR | FLG_MOSEPH, true),
            (f64::NAN, FLG_BARYCTR | FLG_MOSEPH, true),
            (f64::NAN, FLG_SWIEPH, false),
        ] {
            let called = Cell::new(0);
            let resolver = |_: &EphemerisContext, _| {
                called.set(called.get() + 1);
                Err(Error::ConfigurationMissing {
                    missing_config: "unused",
                    suggestion: "unused",
                })
            };
            let error = native()
                .station_query_with_position_context(
                    StationQuery::Velocity {
                        body: MARS,
                        jd_ut: jd,
                        flags,
                        dt: 1e-4,
                    },
                    &resolver,
                )
                .unwrap_err();
            assert_eq!(called.get(), 0);
            if unsupported {
                assert!(matches!(
                    error.primary_error(),
                    Error::UnsupportedFlagCombination { .. }
                ));
            } else {
                assert!(matches!(
                    error.primary_error(),
                    Error::NonFiniteArgument { .. }
                ));
            }
        }
    }
}
