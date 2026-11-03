//! Bounded scalar root/search primitives for the event solvers.
//!
//! Generic infrastructure for the event families (crossings, stations,
//! rise/set, eclipse contacts): bracket validation, bisection, Brent's
//! bracketed method, Newton refinement with a derivative floor, and a
//! uniform first-bracket scan. Serving-layer bindings connect these to actual ephemeris
//! residuals, including wrapped angular differences. The generic root functions
//! impose no body, frame or time semantics; astronomical servings define them.
//!
//! Quantity, units and domain: the caller supplies a scalar function
//! `f(x)` over `f64` and owns the meaning of both axes. `tolerance` is in
//! `|f|` (function) units; `min_step` is an `x`-resolution floor in `x`
//! units. The two must never be compared against each other.
//!
//! Source: the algorithms are the published methods — bisection, Newton
//! iteration, Brent, R. P. (1973) *Algorithms for Minimization without
//! Derivatives*, ch. 4 (Brent–Dekker zero finder) — orchestrated here with
//! the same structural contract as the release `crossing.py` solver it will
//! serve: endpoint-root pre-check, sign-change bracket gate, fallback
//! guards on the interpolation step, and a best-estimate report when the
//! iteration cap is reached. No release coefficient, ephemeris value, or
//! fitted constant is embedded.
//!
//! Numerical method, branch rules and stopping conditions:
//!
//! - [`check_bracket`] evaluates both endpoints once. A non-finite endpoint
//!   value is [`Error::NonFiniteSearchValue`]. An endpoint with
//!   `|f| <= tolerance` is [`BracketOutcome::EndpointRoot`]. A strict sign
//!   change is [`BracketOutcome::Bracketed`]; anything else (same sign,
//!   with neither endpoint a root) is [`BracketOutcome::NoBracket`].
//!   Tangent roots (even multiplicity, no sign change) are therefore
//!   reported as `NoBracket`, never silently refined: they need the
//!   derivative-aware [`newton`] path.
//! - [`bisection`] halves the bracket. It converges when `|f(mid)|` reaches
//!   the tolerance or when the bracket width reaches `min_step` (the root
//!   is then localized to `min_step`; the endpoint with the smaller `|f|`
//!   is reported).
//! - [`brent`] keeps the release solver's structure: the better estimate
//!   stays in `b`, inverse-quadratic interpolation or the secant step is
//!   attempted, and the bisection fallback triggers on the same five
//!   guards (interpolant outside the `(3a+b)/4..b` cone, insufficient
//!   step decrease with/without a previous interpolation, `x`-stagnation
//!   below `min_step` with/without a previous interpolation). Convergence
//!   is declared on `|f(b)| <= tolerance` only — never on `x`-stagnation
//!   alone. Exhausting `max_iter` yields
//!   [`RootOutcome::NonConverged`] carrying the best estimate, mirroring
//!   the release behavior of returning it, but labeled honestly.
//! - [`newton`] iterates `x -= f/df` from `x0`. A step from a derivative
//!   with `|df| < deriv_floor` is refused with `NonConverged` instead of
//!   dividing by a near-zero slope; the station-guard caller
//!   switches to a bracketed method there, as the release solver does near
//!   retrograde stations. A non-finite Newton step also yields
//!   `NonConverged` with the last finite point.
//! - [`scan_first_bracket`] samples `[lo, hi]` uniformly (`count`
//!   sub-intervals, endpoints evaluated exactly once each) and returns the
//!   first strict sign change, or [`ScanOutcome::EndpointRoot`] when a
//!   sample itself satisfies the tolerance. No sign change is
//!   [`ScanOutcome::NoRoot`]. Wrap-aware residuals (angular `±180°`
//!   discontinuities) are the caller's responsibility: this scan treats
//!   `f` as a plain scalar function.
//!
//! Generic root-solver semantics: checked constructors reject invalid
//! configurations, and nonfinite evaluations return `Err(Error)`; use those
//! constructors rather than supplying invalid public configuration fields.
//! Every legitimate search question is an
//! `Ok(outcome)` variant, including `NoBracket`, `NoRoot` and
//! `NonConverged`. An unresolved valid computation is never reported as a
//! certified no-event result: `NoRoot`/`NoBracket` describe the scanned
//! bracket only.
//!
//! Sealed ephemeris-bound serving layer.
//!
//! The remainder of this module binds the solvers above to actual ephemeris
//! words through *sealed servings* ([`SealedSolarGlobalServing`],
//! [`SealedSolarLocalServing`], [`SealedLunarGlobalServing`],
//! [`SealedLunarLocalServing`], [`SealedOccultGlobalServing`],
//! [`SealedOccultLocalServing`], [`SealedPlanetOccultGlobalServing`],
//! [`SealedPlanetOccultLocalServing`], [`SealedTransitServing`],
//! [`SealedStarEventServing`], [`SealedMinorEventServing`],
//! [`SealedFictitiousEventServing`], [`SealedLongitudeServing`]): each one
//! pairs the sealed position pipeline (Sun/Moon/planet geometric states)
//! with the native Delta-T chain and drives one provider-agnostic walk,
//! contact, path, circumstance or magnitude core. Public `..._sealed` and
//! `..._serving` entries expose the release search geometries: global and
//! local solar-eclipse walks with the 10-slot phase layout, central-line
//! and limit samplers, ground width and duration, C1–C4 contacts,
//! Besselian elements and rates, local circumstances and magnitudes,
//! lunar-eclipse walks, contacts, durations, magnitudes and gamma,
//! lunar- and planetary-occultation walks, fixed-instant where/how
//! geometries, rise/set/transit/station searches, the horizon-pair
//! transform and its inverse, and Gauquelin rise/set sectors.
//!
//! Search instants and contact timestamps use UT1 Julian days (`jd_ut`), not
//! UTC. Bound model Delta-T callbacks usually supply TT-minus-UT1 in days for
//! ephemeris evaluation, but do not replace every frame clock. For example,
//! global solar ground geometry and horizon transforms acquire their nutation/
//! obliquity through an independent bundled timescale while retaining the
//! caller's model for other sidereal terms. Read each serving's clock contract;
//! nominally identical Julian-day values do not make those routes interchangeable.
//!
//! Angles are degrees unless stated otherwise. Horizon `azalt` and solar
//! circumstance azimuths are south-zero westward, unlike the north-zero eastward
//! heliacal convention. Path widths are kilometres, Besselian x/y/l1/l2 Earth
//! radii, and durations minutes. Eclipse magnitudes are diameter fractions;
//! obscuration conventions are entry-specific and can be uncapped in totality.
//! Classifications are integer bitmasks with indexed time/attribute arrays.
//! A zero/empty slot can mean absent, inapplicable or unresolved geometry,
//! depending on the entry. It is not a universal proof of physical absence.
//!
//! Invalid checked inputs and consumed provider failures return typed errors,
//! subject to each entry's documented catches and zero-sentinel policies.
//! Exhausted eclipse walks normally return [`Error::EclipseSearchExhausted`].
//! Fixed-star rise/set and transits are supported through
//! [`SealedStarEventServing`], which binds a resolved catalogue star and its
//! separate stellar place/clock feeds. They are not ordinary planet-state reads.
//!
//! Method: Meeus, *Astronomical Algorithms* (1998), for the eclipse
//! geometry implemented here (Besselian elements, contacts, magnitudes,
//! obscuration); the horizon transform follows the release
//! apparent-altitude composition. Per-entry docs state the exact output
//! layout, gates and no-event shapes. Numerical precision and known
//! limitations are described in the maintained precision reference.
//!
//! # Example: distinguish a root estimate from a failed refinement
//!
//! ```
//! use libephemeris_core::events::{bisection, RootOutcome, SearchConfig};
//!
//! let config = SearchConfig::new(1e-10, 1e-12, 64).expect("valid tolerances");
//! let mut residual = |x: f64| x - 2.0;
//! let outcome = bisection(&mut residual, 0.0, 4.0, config).expect("finite samples");
//! match outcome {
//!     RootOutcome::Converged { x, .. } | RootOutcome::EndpointRoot { x, .. } =>
//!         println!("Root estimate: {x}"),
//!     other => println!("No accepted root: {other:?}"),
//! }
//! // Here x has caller-defined units; an astronomical serving supplies its clock.
//! ```

use crate::errors::{Error, Result};

/// Convergence configuration shared by the bracketed refiners.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SearchConfig {
    /// Convergence tolerance in `|f|` (function) units. Must be finite
    /// and positive.
    pub tolerance: f64,
    /// `x`-resolution floor in `x` units: bisection stops narrowing here
    /// and Brent treats smaller steps as stagnation. Must be finite and
    /// positive.
    pub min_step: f64,
    /// Maximum refinement iterations. Must be non-zero.
    pub max_iter: usize,
}

/// Newton configuration: function tolerance, derivative floor and cap.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NewtonConfig {
    /// Convergence tolerance in `|f|` units. Must be finite and positive.
    pub tolerance: f64,
    /// Steps from `|df| < deriv_floor` are refused (station guard).
    /// Must be finite and positive.
    pub deriv_floor: f64,
    /// Maximum iterations. Must be non-zero.
    pub max_iter: usize,
}

/// Outcome of bracket validation ([`check_bracket`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BracketOutcome {
    /// Strict sign change across `[a, b]`; both values finite.
    Bracketed {
        /// Lower endpoint in the caller's independent-variable units.
        a: f64,
        /// Upper endpoint in the caller's independent-variable units.
        b: f64,
        /// Signed finite function value at `a`, in function units.
        fa: f64,
        /// Signed finite function value at `b`, in function units.
        fb: f64,
    },
    /// An endpoint already satisfies the tolerance.
    EndpointRoot {
        /// Accepted endpoint in the caller's independent-variable units.
        x: f64,
        /// Signed function value whose absolute value meets the tolerance.
        fx: f64,
    },
    /// Endpoint checks did not admit a strict sign-change bracket.
    /// This does not establish that the interval contains no root.
    NoBracket {
        /// Signed finite function value at the supplied lower endpoint.
        fa: f64,
        /// Signed finite function value at the supplied upper endpoint.
        fb: f64,
    },
}

/// Outcome of a refinement ([`bisection`], [`brent`], [`newton`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RootOutcome {
    /// A solver stopping criterion was reached within the iteration cap.
    /// This is the function tolerance, or, for bisection, bracket width at
    /// most `min_step`; the latter does not guarantee a small residual.
    Converged {
        /// Returned root estimate in the caller's independent-variable units.
        x: f64,
        /// Signed residual in function units. Width-based bisection termination
        /// may leave its absolute value above the configured function tolerance.
        fx: f64,
        /// Refinement iterations reported by the solver, not callback-call count.
        iterations: usize,
    },
    /// An input endpoint already satisfied the tolerance; no iteration.
    EndpointRoot {
        /// Accepted input endpoint in the caller's independent-variable units.
        x: f64,
        /// Signed endpoint residual satisfying the configured function tolerance.
        fx: f64,
    },
    /// Bracket admission failed, so no bracketed refinement was attempted.
    NoBracket {
        /// Signed finite residual at the supplied lower endpoint.
        fa: f64,
        /// Signed finite residual at the supplied upper endpoint.
        fb: f64,
    },
    /// Refinement stopped without meeting its convergence criterion.
    /// The retained estimate is not a certified root or a global closest point;
    /// Newton retains its last finite evaluated point when it stops.
    NonConverged {
        /// Finite estimate retained at termination; it is not a certified root.
        best_x: f64,
        /// Signed function value at `best_x`, in function units.
        best_fx: f64,
        /// Iterations reported at the cap or earlier termination condition.
        iterations: usize,
    },
}

/// Outcome of a uniform bracket scan ([`scan_first_bracket`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ScanOutcome {
    /// First strict sign change found between two adjacent samples.
    Bracket {
        /// Earlier sample coordinate, in the caller's independent-variable units.
        a: f64,
        /// Adjacent later sample coordinate, in the same units as `a`.
        b: f64,
        /// Signed finite function value at `a`.
        fa: f64,
        /// Signed finite function value at `b`, with the opposite sign.
        fb: f64,
    },
    /// A sample itself satisfies the tolerance.
    EndpointRoot {
        /// Sample coordinate in the caller's independent-variable units.
        x: f64,
        /// Signed function value satisfying the configured tolerance.
        fx: f64,
    },
    /// No sign change and no within-tolerance sample on `[lo, hi]`.
    NoRoot,
}

// Longitude-crossing search layer.

/// Angular convergence tolerance: 0.001 arcsecond in degrees, shared by
/// the Sun, Moon and planet paths.
pub const CROSSING_TOLERANCE_DEG: f64 = 0.001 / 3600.0;

/// Speed below which a body counts as near a retrograde station and the
/// search switches to Brent's method.
pub const STATION_SPEED_THRESHOLD_DEG_DAY: f64 = 0.01;

/// Moon Newton derivative floor, degrees/day: below it the
/// Newton step uses the mean rate. This is NOT the shared 0.01 station
/// threshold the solar path uses.
pub const MOONCROSS_SPEED_FLOOR_DEG_DAY: f64 = 0.1;

/// Brent x-stagnation guard in days.
/// A TIME quantity: never compared against angular residuals.
pub const BRENT_TIME_EPS_DAYS: f64 = 1e-8;

/// Newton iteration caps.
pub const NR_MAX_ITER_SUN: usize = 50;
/// Moon cap.
pub const NR_MAX_ITER_MOON: usize = 50;
/// Generic-planet cap for normal speeds.
pub const NR_MAX_ITER_PLANET: usize = 60;
/// Slow-planet cap.
pub const NR_MAX_ITER_SLOW: usize = 80;
/// Near-station cap.
pub const NR_MAX_ITER_STATION: usize = 100;

/// Sun mean motion seed, degrees/day (Meeus ch. 25/31
/// root-finder seed, non-final).
pub const SUN_MEAN_RATE_DEG_DAY: f64 = 0.9856;
/// Moon mean motion seed, degrees/day (Meeus ch. 47
/// root-finder seed, non-final).
pub const MOON_MEAN_RATE_DEG_DAY: f64 = 13.176;
/// Newton divergence window for the Sun in days.
pub const SUN_DIVERGENCE_DAYS: f64 = 366.0;
/// Newton divergence window for the Moon in days.
pub const MOON_DIVERGENCE_DAYS: f64 = 31.0;

/// Search direction: forward-only (`t >= start`) or strictly-past
/// (`t < start`). The release spells this as `backwards` (bool/int/str via
/// `_coerce_backwards`); the binding layer will coerce to this enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SearchDirection {
    /// First crossing at or after the start epoch.
    Forward,
    /// Nearest crossing strictly before the start epoch.
    Backward,
}

// Station / retrograde / lunar-node search layer.

/// Station convergence tolerance: velocity within this of zero counts as a
/// station.
pub const STATION_VELOCITY_TOLERANCE_DEG_DAY: f64 = 1e-6;

/// Default bracket-scan resolution for one station search window.
pub const STATION_BRACKET_SAMPLES: usize = 50;

/// Forward station search attempts before giving up.
pub const STATION_SEARCH_ATTEMPTS: usize = 4;

/// Classification and advance step in days.
pub const STATION_TYPE_OFFSET_DAYS: f64 = 1.0;

/// Lunar-node latitude tolerance in degrees.
pub const MOON_NODE_TOLERANCE_DEG: f64 = CROSSING_TOLERANCE_DEG;

/// Lunar-node Newton cap.
pub const MOON_NODE_MAX_ITER: usize = NR_MAX_ITER_MOON;

/// Half nodal (draconic) month in days between successive lunar node
/// crossings (Meeus ch. 47 seed, non-final).
pub const HALF_NODAL_MONTH_DAYS: f64 = 13.6;

/// Lunar-node divergence window in days.
pub const MOON_NODE_DIVERGENCE_DAYS: f64 = 30.0;

/// At-node restart gate in days.
pub const MOON_NODE_AT_NODE_GATE_DAYS: f64 = 2.5e-7;

/// Station kind: what the velocity trend does across the station.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StationKind {
    /// Direct-to-retrograde.
    StationaryRetrograde,
    /// Retrograde-to-direct.
    StationaryDirect,
}

/// Which station to find.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StationRequest {
    /// Next station regardless of type.
    Any,
    /// Next direct-to-retrograde station.
    Retrograde,
    /// Next retrograde-to-direct station.
    Direct,
}

/// Nearest-station summary over an abstract speed provider (the
/// ephemeris-free core of release `get_station_info`; the longitude at the
/// station and the string-error partial dict belong to the serving-layer binding).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NearestStation {
    /// Julian Day of the nearest station (past or future).
    pub jd_station: f64,
    /// Its kind.
    pub kind: StationKind,
    /// Signed days from `start` (negative when the station is past).
    pub days_to_station: f64,
    /// Whether the body is retrograde at `start`.
    pub currently_retrograde: bool,
    /// Velocity at `start` in degrees/day.
    pub speed_now_deg_day: f64,
}

/// Lunar-node crossing.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MoonNodeCrossing {
    /// Julian Day of the crossing (same scale as the provider epochs).
    pub jd: f64,
    /// Ecliptic longitude at the crossing in degrees.
    pub longitude_deg: f64,
    /// Ecliptic latitude at the crossing in degrees (within tolerance).
    pub latitude_deg: f64,
}

// Rise/set/transit searches over abstract altitude providers.
use crate::constants::{
    BIT_ASTRO_TWILIGHT, BIT_CIVIL_TWILIGHT, BIT_DISC_BOTTOM, BIT_DISC_CENTER, BIT_NAUTIC_TWILIGHT,
    CALC_ITRANSIT, CALC_MTRANSIT, CALC_RISE, CALC_SET,
};

/// Maximum forward search span for a rise/set crossing, days.
pub const RISE_SET_MAX_SEARCH_DAYS: f64 = 2.0;

/// Altitude residual band for the bisection exit, degrees.
/// A fixed altitude residual maps to a time error of `|h| / |dh/dt|`; on
/// steep crossings this is far below the time tolerance, so the loop exits
/// exactly as before, while near-tangential grazing crossings keep
/// bisecting until the *time* is pinned (see [`RISE_SET_TIME_TOL_S`]).
pub const RISE_SET_ALT_TOLERANCE_DEG: f64 = 1e-4;

/// Bisection iterations over the first-crossing bracket.
pub const RISE_SET_BISECT_ITER: usize = 60;

/// Transit fixed-point iterations.
pub const TRANSIT_MAX_ITER: usize = 30;

/// Transit hour-angle tolerance, degrees.
pub const TRANSIT_TOL_HA_DEG: f64 = 1.0 / 3600.0;

/// Mean sidereal day in solar days, the transit seed and step scale
pub const SIDEREAL_DAY: f64 = 0.99726957;

/// Target time accuracy for near-tangential grazing crossings, seconds
pub const RISE_SET_TIME_TOL_S: f64 = 0.03;

/// Symmetric slope-probe half-width, seconds.
pub const RISE_SET_SLOPE_PROBE_S: f64 = 0.05;

/// Slope below which a crossing counts as near-tangential grazing,
/// degrees/second. Only this
/// population keeps the time-based exit; steeper crossings keep the
/// historic altitude exit plus a one-step Newton refinement.
pub const RISE_SET_GRAZE_SLOPE_DEG_S: f64 = 2.0 / 3600.0;

/// Slope above which the apparent-altitude jump is the sub-polar
/// refraction "dip" discontinuity rather than a genuine crossing,
/// degrees/second. The dip branch keeps the
/// historic midpoint so the search still stops on the refracted branch
/// above the jump.
pub const RISE_SET_DIP_SLOPE_DEG_S: f64 = 100.0 / 3600.0;

/// `horhgt` sentinel requesting the dip of the sea horizon at the
/// observer's altitude: not an angle.
pub const HORIZON_DIP_SENTINEL_DEG: f64 = -100.0;

/// Epsilon past the dip edge the sentinel resolves to, degrees: just past the refraction dead band, so the effect is
/// determinate.
pub const HORIZON_DIP_EPSILON_DEG: f64 = 1e-4;

/// Sea-level pressure for the barometric estimate, hPa.
pub const SEA_LEVEL_PRESSURE_HPA: f64 = 1013.25;

/// Twilight geometric-center altitudes, degrees.
pub const CIVIL_TWILIGHT_ALT_DEG: f64 = -6.0;
/// Nautical twilight center altitude, degrees.
pub const NAUTIC_TWILIGHT_ALT_DEG: f64 = -12.0;
/// Astronomical twilight center altitude, degrees.
pub const ASTRO_TWILIGHT_ALT_DEG: f64 = -18.0;

/// Pinned Sun distance for the fixed-disc-size convention, AU. Binding-side input to
/// [`disc_altitude_term_deg`].
pub const FIXED_DISC_SUN_DIST_AU: f64 = 1.0;
/// Pinned Moon distance for the fixed-disc-size convention, AU.
pub const FIXED_DISC_MOON_DIST_AU: f64 = 0.00257;

/// Maximum refraction lift used only by the circumpolar pre-check, degrees
pub const REFRACTION_LIFT_MAX_DEG: f64 = 0.75;
/// Maximum disc lift used only by the circumpolar pre-check, degrees
pub const DISC_LIFT_MAX_DEG: f64 = 0.30;

/// Rise/set/transit event selected from the `rsmi` bitmask.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RiseSetEvent {
    /// Body crossing the horizon going up (`CALC_RISE`).
    Rise,
    /// Body crossing the horizon going down (`CALC_SET`).
    Set,
    /// Upper meridian transit / culmination (`CALC_MTRANSIT`).
    UpperTransit,
    /// Lower meridian transit / anti-culmination (`CALC_ITRANSIT`).
    LowerTransit,
}

/// Twilight geometric-center override. Honored for the Sun and for fixed stars
/// only; planets and the Moon keep their ordinary horizon crossing. The
/// override forces the `BIT_NO_REFRACTION | BIT_DISC_CENTER` combination
/// and replaces the horizon altitude; civil has precedence over nautical
/// over astronomical.
pub struct TwilightOverride {
    /// Geometric-center altitude replacing `horhgt`, degrees.
    pub horizon_alt_deg: f64,
    /// Always true: twilight forces `BIT_NO_REFRACTION`.
    pub force_no_refraction: bool,
    /// Always true: twilight forces `BIT_DISC_CENTER`.
    pub force_disc_center: bool,
}

/// Which point of the disc crosses the horizon.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiscPoint {
    /// Upper limb (default): the true altitude is raised by the
    /// topocentric semidiameter, computed from the actual distance.
    UpperLimb,
    /// Disc center (`BIT_DISC_CENTER`, always for fixed stars).
    Center,
    /// Lower limb (`BIT_DISC_BOTTOM`): a full diameter away from the
    /// upper limb, so the event moves by roughly twice the semidiameter
    /// in time.
    LowerLimb,
}

/// Meridian-transit solution. The release always reports `res = 0`, so the
/// outcome carries whether the hour-angle tolerance was actually reached:
/// `converged == false` is the best estimate after [`TRANSIT_MAX_ITER`]
/// steps, never a certified transit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TransitSolution {
    /// Julian Day (UT) of the transit estimate.
    pub jd_ut: f64,
    /// Whether `|hour angle| < 1 arcsecond` was reached.
    pub converged: bool,
}

/// Rise/set search outcome. Both no-event paths are the release `res = -2`
/// with `tret[0] = 0.0`: the analytic never-above/never-below verdict and
/// an exhausted 2-day scan are observationally identical to the caller.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RiseSetOutcome {
    /// Event found: Julian Day (UT) of the rise/set.
    Found {
        /// Rise/set epoch as a UT1 Julian day, within the searched window.
        jd_ut: f64,
    },
    /// Circumpolar body, or no crossing in the search window.
    Circumpolar,
}

// Named shadow-geometry data and fixed-state eclipse helpers.
use crate::eclipse_tables::{
    LUNAR_INEX_REFERENCES, LUNAR_SAROS_REFERENCES, SOLAR_INEX_REFERENCES, SOLAR_SAROS_REFERENCES,
};

/// Minimum eclipse magnitude for reliable contact times.
pub const SHALLOW_ECLIPSE_MAG_THRESHOLD: f64 = 0.01;
/// Margin from the gamma limit marking near-miss edge handling.
pub const NEAR_MISS_GAMMA_MARGIN: f64 = 0.02;
/// Separation below which the lens branch is skipped.
pub const MINIMUM_SEPARATION_FOR_LENS_DEG: f64 = 1e-10;
/// Saros cycle period, days.
pub const SAROS_CYCLE_DAYS: f64 = 6585.3211;
/// Inex cycle period, days.
pub const INEX_CYCLE_DAYS: f64 = 10571.9509;
/// Half-width of the Besselian central-difference stencil, days.
pub const BESSELIAN_RATE_STEP_DAYS: f64 = 1.0 / 1440.0;
/// Penumbral lunar radius ratio.
pub const BESSELIAN_K_PENUMBRA: f64 = 0.2725076;
/// Umbral lunar radius ratio.
pub const BESSELIAN_K_UMBRA: f64 = 0.272281;
/// Earth equatorial radius in AU.
pub const ECL_REARTH_AU: f64 = 6378.140 / 149597870.700;
/// Solar radius in AU.
pub const ECL_RSUN_AU: f64 = 696000.0 / 149597870.700;
/// No-match sentinel for the saros lookup.
pub const SAROS_NO_MATCH: f64 = -99999999.0;
/// Residual acceptance window around the mean-cycle grid, days.
pub const SAROS_MATCH_WINDOW_DAYS: f64 = 2.0;
/// Partial-eclipse gamma visibility limit.
pub const ECLIPSE_GAMMA_LIMIT_PARTIAL: f64 = 1.55;
/// Eclipse-magnitude clamp ceiling.
pub const ECLIPSE_MAGNITUDE_MAX: f64 = 1.5;

/// Eclipse family selector for the cycle-number helpers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EclipseKind {
    /// Solar eclipses.
    Solar,
    /// Lunar eclipses.
    Lunar,
}

/// Besselian elements for a solar eclipse at a reference UT1 Julian day.
///
/// x/y and cone radii l1/l2 are Earth radii; declination d and Greenwich hour
/// angle mu are degrees. Rates are per **hour**, not per day: Earth radii/hour
/// or degrees/hour. Construction checks finiteness only, without normalizing
/// angular ranges or certifying physical shadow geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BesselianElements {
    t0: f64,
    x: f64,
    y: f64,
    d: f64,
    l1: f64,
    l2: f64,
    mu: f64,
    dx_dt: f64,
    dy_dt: f64,
    dd_dt: f64,
    dl1_dt: f64,
    dl2_dt: f64,
    dmu_dt: f64,
}

/// The nine `_besselian_core` words: fundamental-plane `(x, y)` in Earth
/// radii, axis declination/hour angle in degrees, cone radii and their
/// tangents, and the Moon's axial distance in Earth radii.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BesselianCore {
    x: f64,
    y: f64,
    d_deg: f64,
    mu_deg: f64,
    l1: f64,
    l2: f64,
    tan_f1: f64,
    tan_f2: f64,
    z_moon_er: f64,
}

/// Per-hour Besselian rates (see [`besselian_rates_via`]).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BesselianRates {
    dx_dt: f64,
    dy_dt: f64,
    dd_dt: f64,
    dl1_dt: f64,
    dl2_dt: f64,
    dmu_dt: f64,
}

// Solar/lunar where/how fixed-state cores
use crate::constants::{
    ECL_ANNULAR, ECL_CENTRAL, ECL_NONCENTRAL, ECL_PARTIAL, ECL_PENUMBRAL, ECL_TOTAL, ECL_VISIBLE,
};

/// Lunar radius in AU.
pub const ECL_RMOON_AU: f64 = 1738.15 / 149597870.700;
/// Earth flattening for the shadow-axis ellipsoid.
pub const ECL_EARTH_FLATTENING: f64 = 1.0 / 298.25642;
/// Astronomical unit in km.
pub const ECL_AU_KM: f64 = 149597870.700;
/// Danjon atmospheric scale on the lunar horizontal-parallax term
/// (Espenak & Meeus,
/// NASA/TP-2009-214173: scaling only the leading Earth-radius term
/// implements `Ru = 1.01*Pm - Ss + Ps`, not Chauvenet's 1/50 radius
/// convention).
pub const DANJON_MOON_PARALLAX_SCALE: f64 = 1.01;
/// Horizon-refraction constant in the solar visibility gate, arcminutes
pub const SOL_HOW_HMIN_REFRACTION_ARCMIN: f64 = 34.4556;
/// Horizon-dip coefficient in the solar visibility gate, arcminutes per
/// square-root metre of observer height.
pub const SOL_HOW_HMIN_DIP_COEF: f64 = 2.12;

/// Shadow-axis ground geometry from `_eclipse_where_core`: global
/// `retc` (0 when no eclipse is in progress anywhere — the carried
/// lon/lat then hold the closest-approach point — else `ECL_CENTRAL` or
/// `ECL_NONCENTRAL`, plus `ECL_PARTIAL` when only the penumbra touches,
/// plus `ECL_TOTAL`/`ECL_ANNULAR` for the non-partial case), the
/// shadow-center lon/lat in degrees (lon in `(-180, 180]`), and the
/// 7-slot shadow-geometry `dcore` layout padded to the release 10-tuple:
/// `[0]` core-shadow diameter at the ground point, km (negative when the
/// umbra reaches the surface), `[1]` penumbra diameter there, km,
/// `[2]` axis distance from the geocenter, km, `[3]` core diameter on the
/// fundamental plane, km, `[4]` penumbra diameter on the fundamental
/// plane, km, `[5]` cos of the core-cone half angle, `[6]` cos of the
/// penumbra-cone half angle, `[7..9]` reserved zeros.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WhereGeometry {
    retc: i32,
    lon_deg: f64,
    lat_deg: f64,
    dcore: [f64; 10],
}

/// Local solar-eclipse circumstances from `_sol_how_core`: local phase
/// `retc` (`ECL_TOTAL`/`ECL_ANNULAR`/`ECL_PARTIAL`, 0 when no eclipse is
/// in progress at this place and time) plus `ECL_VISIBLE` when part of
/// the eclipsed Sun can stand above the local horizon, and the 20-slot
/// `attr` layout — `[0]` magnitude as diameter fraction (negative when
/// the limbs do not overlap), `[1]` lunar/solar diameter ratio,
/// `[2]` obscuration (an uncapped lunar/solar area ratio in totality for
/// the solar convention, with `1.0` as the raw no-event sentinel),
/// `[3]` 0.0 here (the composition fills the core-shadow width),
/// `[4]` south-zero westward azimuth in degrees, `[5]` unrefracted true
/// altitude in degrees, `[6]` atmospherically apparent altitude in degrees of the
/// occulted body, `[7]` center separation in degrees, `[8]` NASA magnitude,
/// `[9]`/`[10]` saros series and member, `[11..19]` zeros.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolHowCircumstances {
    retc: i32,
    attr: [f64; 20],
}

/// Geocentric lunar-eclipse circumstances from `_lun_how_core`:
/// `retc` (`ECL_TOTAL`, `ECL_PARTIAL`, `ECL_PENUMBRAL` or 0), the
/// 20-slot `attr` layout — `[0]` umbral magnitude (0 for penumbral
/// eclipses), `[1]` penumbral magnitude (negative when the Moon stands
/// clear of the penumbra), `[2]`/`[3]` zeros, `[4..6]` zeros here (the
/// composition fills azimuth/true/apparent altitude of the Moon),
/// `[7]` Moon's distance from opposition in degrees (0 when no eclipse),
/// `[8]` `= attr[0]`, `[9]`/`[10]` lunar saros series and member,
/// `[11..19]` zeros — and the AU-based 10-slot `dcore`: `[0]` shadow-axis
/// distance from the selenocenter, `[1]` core diameter on the
/// fundamental plane, `[2]` penumbra diameter there, `[3]` cos of the
/// core-cone half angle, `[4]` cos of the penumbra-cone half angle,
/// `[5..9]` reserved zeros.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LunHowCircumstances {
    retc: i32,
    attr: [f64; 20],
    dcore: [f64; 10],
}

// Global eclipse searches and contact mapping
use crate::constants::{ECL_ALLTYPES_LUNAR, ECL_ANNULAR_TOTAL, ECL_ONE_TRY};

/// Finite safety backstop for the global eclipse lunation walks, years
/// The walks stop earlier at
/// the active ephemeris boundary; this only bounds a never-matching
/// filter on very long ephemerides.
pub const ECLIPSE_SEARCH_HORIZON_YEARS: f64 = 2000.0;
/// Lunations per year sizing the walk backstop.
pub const ECLIPSE_WALK_LUNATIONS_PER_YEAR: f64 = 12.4;
/// Maximum lunation steps per global search.
pub const ECLIPSE_WHEN_MAX_LUNATIONS: usize = 24_800;
/// Epoch-vs-maximum coincidence margin, days: a candidate maximum
/// within this margin of the search epoch counts as already reached, so
/// restarting a search from a returned maximum advances to the
/// neighbouring event.
pub const ECLIPSE_WHEN_EPOCH_MARGIN_DAYS: f64 = 100.0e-7;
/// Coarse new-moon node pre-filter, degrees. Only a pre-filter: the true eclipse test
/// is the classification at the refined maximum.
pub const SOLAR_NODE_PREFILTER_DEG: f64 = 19.0;
/// Coarse full-moon node pre-filter, degrees. Only a pre-filter: the true test is the
/// shadow-geometry classification at the refined maximum.
pub const LUNAR_NODE_PREFILTER_DEG: f64 = 20.0;
/// Bidirectional/forward pre-probe window, days:
/// the previous (solar: either) adjacent syzygy is probed first so an
/// in-progress eclipse straddling the search epoch is not skipped.
pub const ECLIPSE_PRE_PROBE_WINDOW_DAYS: f64 = 15.0;
/// Lunation-walk advance past a consumed conjunction, days.
pub const ECLIPSE_WALK_ADVANCE_DAYS: f64 = 25.0;
/// Earth equatorial radius in km for the eclipse shadow geometry
pub const ECLIPSE_EARTH_RADIUS_KM: f64 = 6378.137;
/// Solar photospheric radius in km for the penumbral/umbral cone limits
pub const ECLIPSE_SUN_RADIUS_KM: f64 = 696340.0;
/// Lunar radius in km for the penumbral/umbral cone limits.
pub const ECLIPSE_MOON_RADIUS_KM: f64 = 1737.4;
/// Atmospheric enlargement of the lunar shadow radii.
pub const LUNAR_SHADOW_ENLARGEMENT: f64 = 1.0 + 1.0 / 85.0;
/// Solar angular semi-diameter at 1 AU, arcsec.
pub const SUN_SEMIDIAMETER_ARCSEC: f64 = 959.63;
/// Lunar angular semi-diameter at the reference distance, arcsec
pub const MOON_SEMIDIAMETER_ARCSEC: f64 = 932.56;
/// Reference distance for the lunar semi-diameter, AU.
pub const MOON_SEMIDIAMETER_REF_AU: f64 = 0.002569;
/// Earth radius in AU for the lunar shadow geometry.
pub const LUNAR_EARTH_RADIUS_AU: f64 = 4.2635e-5;
/// Mean Sun-Moon relative longitude speed seeding the syzygy searches,
/// deg/day.
pub const SYZYGY_MEAN_RELATIVE_SPEED_DEG_DAY: f64 = 12.190749;
/// Syzygy Newton-step speed floor, deg/day: below it the
/// step uses [`SYZYGY_FALLBACK_SPEED_DEG_DAY`].
pub const SYZYGY_SPEED_FLOOR_DEG_DAY: f64 = 0.1;
/// Syzygy Newton-step fallback speed, deg/day.
pub const SYZYGY_FALLBACK_SPEED_DEG_DAY: f64 = 12.19;
/// Syzygy refinement convergence, degrees.
pub const SYZYGY_CONVERGENCE_DEG: f64 = 1e-5;
/// Syzygy Newton-Raphson refinement iterations.
pub const SYZYGY_REFINEMENT_ITERATIONS: usize = 20;
/// Golden-section search iterations.
pub const GOLDEN_SECTION_ITERATIONS: usize = 60;
/// Solar gamma-minimum search half-window, days.
pub const SOLAR_MAXIMUM_SEARCH_RANGE_DAYS: f64 = 0.125;
/// Solar gamma-minimum convergence, days.
pub const SOLAR_MAXIMUM_TOLERANCE_DAYS: f64 = 1e-8;
/// Generic golden-section convergence quoted for `_golden_min`
pub const GOLDEN_MIN_TOLERANCE_DAYS: f64 = 1e-7;
/// Lunar selenocentric-depth search half-window, days.
pub const LUNAR_MAXIMUM_SEARCH_RANGE_DAYS: f64 = 0.3;
/// Lunar phase-contact search half-window, days.
pub const LUNAR_CONTACT_SEARCH_RANGE_DAYS: f64 = 4.0 / 24.0;
/// Besselian contact bracket-expansion cap, days.
pub const BESSELIAN_CONTACT_MAX_WINDOW_DAYS: f64 = 0.35;
/// Besselian contact bracket-expansion seed step, days
/// with the release `1.6` geometric growth.
pub const BESSELIAN_CONTACT_STEP_DAYS: f64 = 0.02;
/// Contact bisection convergence: bracket width in days and gamma-residual gate.
pub const CONTACT_BISECT_TOLERANCE_DAYS: f64 = 1e-8;
/// Hybrid-probe gamma contact search range, days.
pub const HYBRID_CONTACT_SEARCH_RANGE_DAYS: f64 = 0.10;
/// Hybrid fast-reject band around zero fundamental-plane `l2`, Earth
/// radii.
pub const HYBRID_L2_FAST_REJECT_NEG: f64 = -0.001;
/// Hybrid fast-reject band upper edge: fundamental-plane `l2` above
/// `0.006` Earth radii skips the contact solve.
pub const HYBRID_L2_FAST_REJECT_POS: f64 = 0.006;

/// Search direction for the global solar search.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EclipseSearchDirection {
    /// Only eclipses with maximum after the epoch (reference-API direction).
    Forward,
    /// Only eclipses with maximum before the epoch.
    Backward,
    /// Probe ±15 days first so an in-progress eclipse is not missed;
    /// otherwise forward.
    Bidirectional,
}

/// A public `backwards` argument:
/// booleans pass through, integers use truthiness, strings
/// (case-insensitive, whitespace-stripped) map the backward/forward
/// spellings explicitly so a truthy `"forward"` never selects a backward
/// search. Rust has no dynamic-type fallback: unsupported static types
/// cannot occur, and unknown strings are `InvalidEclipseSearch` carrying
/// the release `TypeError` text.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BackwardsArg<'a> {
    /// Boolean direction flag.
    Bool(bool),
    /// Integer direction flag (truthiness).
    Int(i64),
    /// Direction string (`"backward"`/`"forward"` spellings).
    Text(&'a str),
}

/// Apparent Sun/Moon ecliptic state for the syzygy and classification
/// providers. The ephemeris
/// acquisition stays provider-side; this bundle only carries the words
/// the searches consume.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SunMoonState {
    sun_lon_deg: f64,
    sun_lat_deg: f64,
    sun_dist_au: f64,
    sun_speed_deg_day: f64,
    moon_lon_deg: f64,
    moon_lat_deg: f64,
    moon_dist_au: f64,
    moon_speed_deg_day: f64,
}

/// One golden-section decision row: the bracket and objective words at
/// decision time, the arm taken, and the width after the update.
/// Degrees/days follow [`golden_section_min`]; all words are the live
/// computation values, recorded without rounding or normalization.
#[derive(Clone, Copy, Debug)]
pub struct GoldenTraceRow {
    /// Zero-based decision index (the initial probe pair is level 0).
    pub level: usize,
    /// Left interior point word at decision time.
    pub jd_a: f64,
    /// Right interior point word at decision time.
    pub jd_b: f64,
    /// Objective word at `jd_a`.
    pub f_a: f64,
    /// Objective word at `jd_b`.
    pub f_b: f64,
    /// `true` when the `f_a < f_b` arm fired (`hi = b`), else the
    /// `lo = a` arm.
    pub arm_a_taken: bool,
    /// Bracket width after the update.
    pub width_after: f64,
}

/// Signed Besselian contact residuals on the Earth ellipsoid, km
/// each residual is negative while
/// the corresponding shadow feature touches the ellipsoid and positive
/// once cleared — `pen` for the penumbral edge (P1/P4), `umb` for the
/// umbral/antumbral edge (U1/U4), `cen` for the shadow axis itself
/// (center line). `where_dcore` is the [`eclipse_where_core_via`]
/// 10-slot tuple (`r0`, fundamental-plane core/penumbra diameters and
/// cone cosines live in slots 2..6). Contacts are cone-tangency instants:
/// with the polar axis stretched so the flattened Earth counts as a
/// sphere, a shadow cone of half-angle `f` and fundamental-plane radius
/// `L` grazes a sphere of radius `Re` whose center stands `r0` off the
/// axis when `(r0 − L)·cos f = Re`, to second order in `f` (the cone
/// depth a flat circle-on-the-plane test would miss is carried by the
/// `cos f` factor). The reduction follows the Explanatory Supplement to
/// the Astronomical Almanac, 3rd ed., ch. 11, and Meeus, Elements of
/// Solar Eclipses, ch. 4.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolarContactResiduals {
    /// Penumbral-edge residual, km (slot 0).
    pub pen_km: f64,
    /// Umbral/antumbral-edge residual, km (slot 1).
    pub umb_km: f64,
    /// Shadow-axis residual, km (slot 2).
    pub cen_km: f64,
}

/// Call-local memo for the `(pen, umb, cen)` residual triple. Keyed by `f64` bits; slot selects the
/// component. A non-boundary provider failure degrades to `None` (the
/// contact-solve sentinel path); a coverage miss propagates.
struct MemoTriple<'a, F>
where
    F: FnMut(f64) -> Result<(f64, f64, f64)>,
{
    inner: &'a mut F,
    memo: Vec<(u64, (f64, f64, f64))>,
}

/// Solar-eclipse classification at maximum: angular-size ratio from the apparent distances, gamma from
/// the caller (computed via [`gamma_from_vectors`]), shadow limits via
/// [`penumbra_limit_er`] / [`umbra_limit_er`], and the graze-boundary
/// gate from the ellipsoidal where-geometry (`where_retc` /
/// `where_magnitude`, i.e. `_eclipse_where_core` retc plus
/// `_where_eclipse_magnitude`). The centrality bit comes from the same
/// where geometry the compositions use — never a spherical gamma test.
///
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolarEclipseClassification {
    /// Type bitmask (`ECL_TOTAL` / `ECL_ANNULAR` / `ECL_PARTIAL` /
    /// `ECL_ANNULAR_TOTAL` with `ECL_CENTRAL` / `ECL_NONCENTRAL`), or `0`
    /// for no eclipse.
    pub eclipse_type: i32,
    /// Eclipse magnitude (Sun-diameter fraction, clamped to `[0, 1.5]`).
    pub magnitude: f64,
    /// Gamma at maximum, Earth radii.
    pub gamma: f64,
    /// Apparent Moon/Sun diameter ratio.
    pub moon_sun_ratio: f64,
}

/// Acquisition contract for a global solar search at UT1 Julian days.
///
/// Methods supply the stated geocentric apparent geometry and units; bound
/// servings own TT conversion, independent frame clocks and typed acquisition
/// errors. `&self` permits interior state but does not imply callback purity or
/// authorize replacement of repeated observations by cached values. Default
/// methods consume their delegated feeds; overriding one preserves the other
/// methods' separate acquisition roles.
pub trait SolarGlobalEphemeris {
    /// Moon longitude for the conjunction pre-filter, in ecliptic degrees.
    /// Default serving consumes its existing full state. Native public search
    /// overrides this with the one Moon call carrying caller flags plus SPEED.
    fn candidate_moon_longitude_deg(&self, jd_ut: f64) -> Result<f64> {
        self.sun_moon_state(jd_ut).map(|state| state.moon_lon_deg())
    }
    /// Apparent Sun/Moon ecliptic state (lon/lat/dist/speed words).
    fn sun_moon_state(&self, jd_ut: f64) -> Result<SunMoonState>;
    /// Mean lunar-node longitude, degrees.
    fn mean_node_lon_deg(&self, jd_ut: f64) -> Result<f64>;
    /// Gamma at `jd_ut`, Earth radii.
    fn gamma_er(&self, jd_ut: f64) -> Result<f64>;
    /// Fundamental-plane shadow radii `(l1, l2)` at `jd_ut`, Earth radii.
    fn shadow_limits_er(&self, jd_ut: f64) -> Result<(f64, f64)>;
    /// Where-core retc plus where-magnitude at `jd_ut` (the
    /// graze-boundary gate words).
    fn where_gate(&self, jd_ut: f64) -> Result<(i32, f64)>;
    /// Where-core classification bits, without consuming a local magnitude.
    /// Existing sealed providers retain their pair callback by default.
    fn where_retflag(&self, jd_ut: f64) -> Result<i32> {
        self.where_gate(jd_ut).map(|pair| pair.0)
    }
    /// Local maximum-site magnitude after one global geometry selection.
    /// Native providers override this role with the actual conditional calls.
    fn where_magnitude(&self, jd_ut: f64) -> Result<f64> {
        self.where_gate(jd_ut).map(|pair| pair.1)
    }

    /// `(pen, umb, cen)` contact residuals at `jd_ut`, km.
    fn contact_residuals_km(&self, jd_ut: f64) -> Result<(f64, f64, f64)>;
    /// Besselian-x at `jd_ut`, Earth radii.
    fn besselian_x_er(&self, jd_ut: f64) -> Result<f64>;
    /// Ground-point core-shadow width (`dcore[0]`) at `jd_ut`, km, for
    /// the hybrid probe.
    fn ground_core_width_km(&self, jd_ut: f64) -> Result<f64>;
}

/// Lunar-eclipse classification at maximum: Earth-shadow cone sizes at the
/// Moon's distance with Danjon's `1 + 1/85` atmospheric enlargement,
/// Moon distance from the shadow axis (anti-solar point with the
/// latitude-cosine longitude term), gamma in Earth radii, and the
/// penumbral/umbral magnitude arms with the `0.01` shallow gates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LunarEclipseClassification {
    /// Type bits (`ECL_TOTAL` / `ECL_PARTIAL` / `ECL_PENUMBRAL`), or `0`
    /// for no eclipse.
    pub eclipse_type: i32,
    /// Umbral magnitude (Moon diameter fraction in umbra).
    pub magnitude_umbral: f64,
    /// Penumbral magnitude.
    pub magnitude_penumbral: f64,
    /// Gamma: Moon distance from the shadow axis, Earth radii.
    pub gamma: f64,
    /// Penumbral shadow radius at the Moon's distance, degrees.
    pub penumbra_radius_deg: f64,
    /// Umbral shadow radius at the Moon's distance, degrees.
    pub umbra_radius_deg: f64,
}

/// Canonical lunar shadow words for the phase-contact solvers: the
/// [`lun_how_core_via`] AU tuple — Moon-axis offset `r0`, umbral and
/// penumbral diameters, and the two cone cosines.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LunarShadowDcore {
    r0_au: f64,
    umbral_diam_au: f64,
    penumbral_diam_au: f64,
    cos_f_umbral: f64,
    cos_f_penumbral: f64,
}

/// Ephemeris acquisition for the global lunar search. All methods take
/// `&self` (read-only search, as in the `StateProvider` state-provider trait). The
/// serving-layer binding implements this over the sealed/native
/// pipeline; tests implement it over analytic models.
/// Geocentric equatorial `(moon, sun)` AU vector pair shared by the
/// lunar global-search acquisition trait and its implementations.
pub type MoonSunXyzAu = ((f64, f64, f64), (f64, f64, f64));

/// Read-only ephemeris acquisition for the global lunar-eclipse search.
/// Every method samples the Sun/Moon geometry at a UT Julian day and
/// reports computation failures as errors; the sealed serving binding
/// implements this over the sealed/native pipeline, tests over analytic
/// models.
pub trait LunarGlobalEphemeris {
    /// Moon-only candidate longitude in degrees, sampled at UT1 with the
    /// caller's eclipse flags plus SPEED. Implementations with a distinct
    /// public scalar acquisition override this channel; the default preserves
    /// the shared apparent state of existing checked providers.
    fn candidate_moon_longitude_deg(&self, jd_ut: f64) -> Result<f64> {
        Ok(self.sun_moon_state(jd_ut)?.moon_lon_deg)
    }
    /// Apparent Sun/Moon ecliptic state (lon/lat/dist/speed words).
    fn sun_moon_state(&self, jd_ut: f64) -> Result<SunMoonState>;
    /// Mean lunar-node longitude, degrees.
    fn mean_node_lon_deg(&self, jd_ut: f64) -> Result<f64>;
    /// Geocentric equatorial `(moon, sun)` AU vectors.
    fn geocentric_xyz_au(&self, jd_ut: f64) -> Result<MoonSunXyzAu>;
    /// Canonical shadow classification plus shadow words.
    fn lun_shadow(&self, jd_ut: f64) -> Result<(i32, LunarShadowDcore)>;
}

// Local eclipse search/visibility and maxima.
use crate::constants::{
    ECL_1ST_VISIBLE, ECL_2ND_VISIBLE, ECL_3RD_VISIBLE, ECL_4TH_VISIBLE, ECL_MAX_VISIBLE,
    ECL_PARTBEG_VISIBLE, ECL_PARTEND_VISIBLE, ECL_PENUMBBEG_VISIBLE, ECL_PENUMBEND_VISIBLE,
    ECL_TOTBEG_VISIBLE, ECL_TOTEND_VISIBLE, FLG_JPLEPH, FLG_MOSEPH, FLG_SWIEPH,
};

/// Epoch-vs-local-contact coincidence window, days: location searches anchor
/// their candidate lunations on the geocentric maximum, but the
/// observer's local phases can still unfold after that instant, so the
/// candidate generator is anchored this far before/after the epoch and
/// the local-contact epoch gate below keeps or skips it.
pub const WHEN_LOC_EPOCH_WINDOW_DAYS: f64 = 0.35;
/// Inner-contact ("umbral") lunar radius in AU: second /
/// third contact reckon the disappearance/reappearance of the solar
/// disc against this slightly smaller disc while first/fourth use the
/// full [`ECL_RMOON_AU`] silhouette.
pub const ECL_RMOON_INNER_AU: f64 = 0.272281 * 6378.140 / 149597870.700;
/// Half-window of the local-maximum golden search, days.
pub const LOCAL_MAX_SEARCH_RANGE_DAYS: f64 = 3.0 / 24.0;
/// Golden-section iterations of the served local-maximum search
pub const LOCAL_MAX_ITERATIONS: usize = 40;
/// Golden-section width stop of the served local-maximum search,
/// days.
pub const LOCAL_MAX_TOLERANCE_DAYS: f64 = 1e-7;
/// Golden-section iterations of the legacy local-phase search
pub const LOCAL_PHASES_GOLDEN_ITERATIONS: usize = 30;
/// Golden-section width stop of the local searches, days.
pub const LOCAL_PHASES_TOLERANCE_DAYS: f64 = 1e-7;
/// Half-window of the outer-contact search around the local maximum,
/// days.
pub const LOCAL_CONTACT_SEARCH_DAYS: f64 = 2.0 / 24.0;
/// Half-window of the inner-contact search around the local maximum,
/// days.
pub const LOCAL_INNER_CONTACT_SEARCH_DAYS: f64 = 0.5 / 24.0;
/// Bisection iterations of the legacy local contact solve.
pub const LOCAL_CONTACT_BISECT_ITERATIONS: usize = 50;
/// Separation-match stop of the legacy local contact solve, degrees
pub const LOCAL_CONTACT_SEP_TOLERANCE_DEG: f64 = 1e-6;
/// Post-solve validation bound of the legacy local contact solve,
/// degrees.
pub const LOCAL_CONTACT_VALID_TOLERANCE_DEG: f64 = 0.001;
/// True-altitude gate below which the Sun (solar) or the eclipse is
/// treated as not visible, degrees.
pub const LOCAL_SUN_ALT_GATE_DEG: f64 = -1.0;
/// Maximum lunations examined by the solar local walk.
pub const SOL_LOC_MAX_ECLIPSES: usize = 120;
/// Maximum lunations examined by the lunar local walk.
pub const LUN_LOC_MAX_ITERATIONS: usize = 250;
/// Backward-walk collection reach, days.
pub const SOL_LOC_BACKWARD_COLLECT_DAYS: f64 = 400.0;

/// Ten-slot local eclipse circumstances: `[0]` local
/// maximum, `[1]` first contact, `[2]` second contact, `[3]` third
/// contact, `[4]` fourth contact, `[5]` magnitude, `[6]` fraction of
/// the solar diameter covered, `[7]` obscuration (area), `[8]` solar
/// azimuth at maximum, `[9]` solar altitude at maximum. All zeros
/// means no local eclipse (Sun down, no disc overlap).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocalEclipsePhases {
    /// Ten-slot circumstances in the release order.
    pub slots: [f64; 10],
}

/// Lower-limb horizon crossings of a body around a local eclipse
/// return codes plus the
/// `tret[0]` instants. `rc == -2` on either side marks the
/// circumpolar Moon/Sun that leaves the geometric phase times
/// untouched.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyHorizonCrossings {
    /// Return code of the rise search (`-2` circumpolar).
    pub rise_rc: i32,
    /// Rise instant (`tret[0]`), JD UT.
    pub rise_jd: f64,
    /// Return code of the set search (`-2` circumpolar).
    pub set_rc: i32,
    /// Set instant (`tret[0]`), JD UT.
    pub set_jd: f64,
}

/// Ephemeris acquisition for the served solar local search. All methods take
/// `&self` (read-only search): interior mutability covers future
/// caching. The serving-layer binding implements this over the
/// sealed/native pipeline; tests implement it over analytic models.
pub trait SolarLocalEphemeris {
    /// Next global solar maximum at or after `jd` (the
    /// global solar-eclipse search): `(type, tret10)`.
    /// A coverage miss is reported as an ephemeris-boundary error so
    /// the walk converts it to exhaustion.
    fn next_global_forward(&self, jd: f64) -> Result<(i32, [f64; 10])>;
    /// Topocentric Sun-Moon center separation, degrees.
    fn separation_deg(&self, jd: f64) -> Result<f64>;
    /// Apparent radii plus center separation `(rsun, rmoon, sep)`,
    /// degrees; `inner_moon` selects the smaller
    /// [`ECL_RMOON_INNER_AU`] lunar disc for second/third contact.
    fn radii_sep_deg(&self, jd: f64, inner_moon: bool) -> Result<(f64, f64, f64)>;
    /// Solar `(true_altitude_deg, azimuth_deg, apparent_altitude_deg)`.
    fn sun_altaz_deg(&self, jd: f64) -> Result<(f64, f64, f64)>;
    /// Topocentric apparent diameters `(moon_diam_deg, sun_diam_deg)`.
    fn topo_diameters_deg(&self, jd: f64) -> Result<(f64, f64)>;
    /// Lower-limb solar horizon crossings from `jd_from`.
    fn sun_horizon_crossings(&self, jd_from: f64) -> Result<BodyHorizonCrossings>;
    /// Solar local-circumstance core at `jd`: `(retc, attr20)` in the
    /// `_sol_how_core` layout.
    fn how_core(&self, jd: f64) -> Result<(i32, [f64; 20])>;
    /// Shadow-axis core at `jd`: `(retc, dcore10)`; `dcore[0]` is the
    /// ground core-shadow width in km.
    fn where_core(&self, jd: f64) -> Result<(i32, [f64; 10])>;
}

/// Ephemeris acquisition for the served lunar local search. All
/// methods take `&self` (read-only search). The serving-layer binding
/// implements this over the sealed/native pipeline; tests
/// implement it over analytic models.
pub trait LunarLocalEphemeris {
    /// Next global lunar eclipse `(retflag, tret10)` in the
    /// `lun_eclipse_when` layout from `jd` (`backwards` selects the
    /// direction).
    fn next_global(&self, jd: f64, backwards: bool) -> Result<(i32, [f64; 10])>;
    /// Lunar local-circumstance core at `jd`: `(retc, attr20)` in the
    /// `_lun_how_core` layout (`attr[6]` is the apparent altitude).
    fn how_attr(&self, jd: f64) -> Result<(i32, [f64; 20])>;
    /// Lower-limb lunar horizon crossings from `jd_from`.
    fn moon_horizon_crossings(&self, jd_from: f64) -> Result<BodyHorizonCrossings>;
}

// Lunar-occultation fixed-state APIs (finite vs point sources).
use crate::constants::{
    CERES, CHIRON, EARTH, JUNO, JUPITER, MARS, MERCURY, MOON, NEPTUNE, PALLAS, PHOLUS, PLUTO,
    SATURN, SUN, URANUS, VENUS, VESTA,
};

/// Ecliptic-latitude bound past which a fixed star can never be occulted,
/// degrees.
pub const OCC_STAR_LAT_LIMIT_DEG: f64 = 7.0;
/// Geocentric Moon-body latitude gap past which no occultation is possible
/// anywhere on Earth at a conjunction, degrees.
pub const OCC_LAT_GATE_DEG: f64 = 2.0;
/// Point-source angular-radius fallback, degrees.
pub const OCC_POINT_SOURCE_RADIUS_DEG: f64 = 0.0001;

/// Fixed-instant `lun_occult_where` outcome: the where-core character
/// bits (0 when no occultation is in progress anywhere — the carried
/// lon/lat then hold the closest-approach point — else `ECL_CENTRAL` or
/// `ECL_NONCENTRAL`, plus `ECL_PARTIAL` for penumbral-only touch, plus
/// `ECL_TOTAL`/`ECL_ANNULAR`; stars are point sources, always total once
/// in progress), the release 10-slot `geopos` (shadow-center lon/lat,
/// remaining slots zero) and the 20-slot local-circumstance `attr` with
/// `attr[3]` set to the core-shadow diameter at the center point in km
/// (negative for a total occultation — for a star, minus the lunar
/// diameter).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OccultWhereOutcome {
    retflag: i32,
    geopos: [f64; 10],
    attr: [f64; 20],
}

// Lunar-occultation global/local searches and contact mapping.
use crate::constants::{ECL_OCC_BEG_DAYLIGHT, ECL_OCC_END_DAYLIGHT};

/// Mean lunar sidereal motion seeding the conjunction Newton step,
/// deg/day.
pub const OCC_MOON_RATE_DEG_DAY: f64 = 13.0;
/// Conjunction convergence: 0.1 deg of longitude is ~11 minutes of lunar
/// motion, well inside the ±1 day refinement window used afterwards
pub const OCC_CONJ_TOL_DEG: f64 = 0.1;
/// Conjunction advance past a consumed candidate, days.
pub const OCC_CONJUNCTION_HOP_DAYS: f64 = 20.0;
/// Newton-iteration backstop for one conjunction seek.
pub const OCC_CONJUNCTION_MAX_STEPS: usize = 300;
/// Walk backstop sizing the occultation searches to the search horizon
pub const OCCULT_MAX_CONJUNCTIONS: usize = 27_000;
/// Earth equatorial radius in km for the occultation shadow-geometry
/// phase equations.
pub const OCC_EARTH_EQ_RADIUS_KM: f64 = 6378.140;
/// Maximum-refinement half-window around a conjunction, days.
pub const OCC_MAX_SEARCH_HALF_WINDOW_DAYS: f64 = 1.0;
/// Global outer-contact search half-window, days, with the `wide = 5.0 / 24.0` fallback probed when the narrow
/// window reports the `0.0` sentinel.
pub const OCC_GLOBAL_OUTER_WINDOW_DAYS: f64 = 2.0 / 24.0;
/// Global outer-contact fallback half-window, days.
pub const OCC_GLOBAL_OUTER_WIDE_DAYS: f64 = 5.0 / 24.0;
/// Local outer-contact half-window, days, with the doubled fallback probed when the narrow window
/// reports the `0.0` sentinel.
pub const OCC_LOCAL_OUTER_WINDOW_DAYS: f64 = 2.0 / 24.0;
/// Local inner-contact half-window, days.
pub const OCC_LOCAL_INNER_WINDOW_DAYS: f64 = 1.5 / 24.0;

/// Ephemeris acquisition for the global occultation search. All methods take `&self`
/// (read-only search, as in [`SolarGlobalEphemeris`]). The serving-layer
/// binding implements this over the sealed/native pipeline;
/// tests implement it over analytic models.
pub trait OccultGlobalEphemeris {
    /// Geocentric apparent ecliptic `(lon, lat)` of the occulted body,
    /// degrees (planets via `calc_ut`, stars via `fixstar_ut`).
    fn body_lonlat_deg(&self, jd_ut: f64) -> Result<(f64, f64)>;
    /// Geocentric apparent ecliptic `(lon, lat)` of the Moon, degrees.
    fn moon_lonlat_deg(&self, jd_ut: f64) -> Result<(f64, f64)>;
    /// Geocentric separation minus the sum of the apparent radii
    /// degrees.
    fn sep_minus_radii_deg(&self, jd_ut: f64) -> Result<f64>;
    /// Shadow-axis ground geometry at `jd_ut`.
    fn where_core(&self, jd_ut: f64) -> Result<WhereGeometry>;
    /// Wrapped body-minus-Moon right-ascension gap, degrees.
    fn ra_gap_deg(&self, jd_ut: f64) -> Result<f64>;
}

/// Ephemeris acquisition for the local occultation search. All methods take `&self`
/// (read-only search). The serving-layer binding implements this over
/// the sealed/native pipeline; tests implement it over analytic models.
pub trait OccultLocalEphemeris {
    /// Geocentric apparent ecliptic `(lon, lat)` of the occulted body,
    /// degrees (conjunction-seek source; stars geocentric).
    fn body_lonlat_deg(&self, jd_ut: f64) -> Result<(f64, f64)>;
    /// Geocentric apparent ecliptic `(lon, lat)` of the Moon, degrees.
    fn moon_lonlat_deg(&self, jd_ut: f64) -> Result<(f64, f64)>;
    /// Topocentric apparent `(body_radius, moon_radius, separation)` at
    /// `jd_ut`, degrees.
    fn topo_radii_sep_deg(&self, jd_ut: f64, inner: bool) -> Result<(f64, f64, f64)>;
    /// Local-circumstance altitude of the occulted body at `jd_ut`,
    /// degrees.
    fn body_alt_deg(&self, jd_ut: f64) -> Result<f64>;
    /// Local circumstances at the maximum.
    fn local_how_attr(&self, jd_ut: f64) -> Result<[f64; 20]>;
    /// Shadow-axis ground geometry at `jd_ut`.
    fn where_core(&self, jd_ut: f64) -> Result<WhereGeometry>;
    /// Body rise/set transits from `jd_from`: `Some(t)` carries
    /// `tr[0]` when `rc >= 0`, `None` when the body never rises/sets.
    fn body_rise_set(&self, jd_from: f64) -> Result<(Option<f64>, Option<f64>)>;
    /// Solar rise/set transits at `jd`: `Some(t)` carries `tr[0]` when `rc >=
    /// 0`, `None` otherwise.
    fn sun_rise_set(&self, jd: f64) -> Result<(Option<f64>, Option<f64>)>;
}

// Sealed lunar-occultation servings.

/// Occulted-body leg bound into the sealed occult servings: a numeric
/// planet id (planets and the Sun travel the `calc_ut` acquisition)
/// or a resolved star (star names and integer fixed-star ids travel
/// the `fixstar_ut` acquisition, geocentric throughout).
#[derive(Clone, Debug)]
pub enum OccultBodyLeg {
    /// Numeric planet id answered by `calc_ut` (Sun included).
    Planet(i32),
    /// Resolved star answered by `fixstar_ut` (name plus match).
    Star {
        /// Resolved star match (id plus canonical name).
        star: crate::stars::StarMatch,
        /// Search name as given (reused for star rise/set acquisition).
        name: String,
    },
}

/// Sealed global lunar-occultation serving: implements
/// [`OccultGlobalEphemeris`] over the sealed/native pipeline for one
/// bound body leg. Site-independent (the global walk is geocentric);
/// the observer triple enters only at the O-A binding for the local
/// slice.
pub struct SealedOccultGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    serving: &'a crate::positions::SealedServing<'a, P>,
    star_serving: &'a crate::stars::SealedStarServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    calc_flags: i32,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
    body: OccultBodyLeg,
}

// Public `azalt` / `azalt_rev` horizon entries

/// Ecliptic-to-horizontal flag for [`azalt_deg`].
pub const ECL2HOR: i32 = 0;

/// Equatorial-to-horizontal flag for [`azalt_deg`].
pub const EQU2HOR: i32 = 1;

/// Horizontal-to-ecliptic flag for [`azalt_rev_deg`].
pub const HOR2ECL: i32 = 0;

/// Horizontal-to-equatorial flag for [`azalt_rev_deg`].
pub const HOR2EQU: i32 = 1;

/// Sealed local lunar-occultation serving: implements
/// [`OccultLocalEphemeris`] over the sealed/native pipeline for one
/// bound body leg and one bound observer triple. The site travels
/// explicitly (explicit-site contract) — there is no global set/restore
/// so concurrent searches never share
/// observer state.
pub struct SealedOccultLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    serving: &'a crate::positions::SealedServing<'a, P>,
    star_serving: &'a crate::stars::SealedStarServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
    calc_flags: i32,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
    body: OccultBodyLeg,
    lon_deg_east: f64,
    lat_deg: f64,
    alt_m: f64,
}

// Planetary-occultation global/local searches and contact mapping.
use crate::positions::body_uses_jpl_ephemeris;

/// Search horizon for the planetary-occultation walk, years.
pub const PLANET_OCCULT_MAX_SEARCH_YEARS: f64 = 150.0;
/// Daily-sample iteration bound.
pub const PLANET_OCCULT_MAX_ITERATIONS: usize = 54_750;
/// Forward end-of-window break.
pub const PLANET_OCCULT_WINDOW_DAYS: f64 = 150.0 * 365.25;
/// Close-approach refinement gate, degrees.
pub const PLANET_OCCULT_CLOSE_SEP_DEG: f64 = 0.2;
/// Close-approach refinement half-window, days.
pub const PLANET_OCCULT_CLOSE_HALF_WINDOW_DAYS: f64 = 0.5;
/// Occultation-hit refinement half-window, days.
pub const PLANET_OCCULT_HIT_HALF_WINDOW_DAYS: f64 = 1.0;
/// Golden-section minimum iterations and width.
pub const PLANET_OCCULT_GOLDEN_ITERATIONS: usize = 50;
/// Golden-section stop width: a bracket narrower than `1e-6` days
/// (≈ 0.1 s) ends the refinement.
pub const PLANET_OCCULT_GOLDEN_TOL_DAYS: f64 = 1e-6;
/// Relative-speed stencil half-step, days.
pub const PLANET_OCCULT_SPEED_SAMPLE_DAYS: f64 = 1.0 / 24.0;
/// Stencil-speed floor and fallback, deg/day.
pub const PLANET_OCCULT_SPEED_FLOOR_DEG_DAY: f64 = 0.001;
/// Fallback relative speed used when the stencil speed falls below the
/// floor, deg/day.
pub const PLANET_OCCULT_SPEED_FALLBACK_DEG_DAY: f64 = 0.5;
/// Adaptive walk steps, days.
pub const PLANET_OCCULT_STEP_CLOSE_DAYS: f64 = 0.01;
/// Near-range adaptive walk step, days.
pub const PLANET_OCCULT_STEP_NEAR_DAYS: f64 = 0.1;
/// Mid-range adaptive walk step, days.
pub const PLANET_OCCULT_STEP_MID_DAYS: f64 = 0.5;
/// Far-range adaptive walk step, days.
pub const PLANET_OCCULT_STEP_FAR_DAYS: f64 = 1.0;
/// Adaptive-step separation gates, degrees.
pub const PLANET_OCCULT_STEP_CLOSE_SEP_DEG: f64 = 0.5;
/// Near-range separation gate selecting the near walk step, degrees
pub const PLANET_OCCULT_STEP_NEAR_SEP_DEG: f64 = 2.0;
/// Mid-range separation gate selecting the mid walk step, degrees
pub const PLANET_OCCULT_STEP_MID_SEP_DEG: f64 = 10.0;
/// Local walk global re-search bound.
pub const PLANET_OCCULT_MAX_GLOBAL_SEARCHES: usize = 100;
/// Local visibility floor, degrees.
pub const PLANET_OCCULT_MIN_ALT_DEG: f64 = -0.5;
/// Local visibility sampling.
pub const PLANET_OCCULT_LOCAL_SAMPLES: usize = 10;
/// Point-source diameter ratio.
pub const PLANET_OCCULT_POINT_SOURCE_RATIO: f64 = 999.0;

/// Ephemeris acquisition for the planetary-occultation walk. All methods take `&self` (read-only
/// search). The serving-layer binding implements this over
/// `calc_ut(FLG_EQUATORIAL | FLG_SPEED)` distances via
/// [`calc_planet_angular_radius_deg`] (planets) and
/// [`OCC_POINT_SOURCE_RADIUS_DEG`] (resolved stars); tests implement it
/// over analytic models.
pub trait PlanetOccultEphemeris {
    /// One occulting-body acquisition, including its distance-derived radius.
    /// Implementations may override this to retain the single scalar position
    /// read that precedes the target read. The default composes legacy getters.
    fn occulting_position(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        let (ra, dec) = self.occulting_ra_dec_deg(jd_ut)?;
        Ok((ra, dec, self.occulting_radius_deg(jd_ut)?))
    }
    /// One target acquisition and its angular radius, after the occulting
    /// place. Star resolution belongs at this consumed target boundary.
    fn target_position(&self, jd_ut: f64) -> Result<(f64, f64, f64)> {
        let (ra, dec) = self.target_ra_dec_deg(jd_ut)?;
        Ok((ra, dec, self.target_radius_deg(jd_ut)?))
    }
    /// Geocentric RA/Dec of the occulting planet, degrees.
    fn occulting_ra_dec_deg(&self, jd_ut: f64) -> Result<(f64, f64)>;
    /// Geocentric RA/Dec of the occulted body (planet or resolved star),
    /// degrees.
    fn target_ra_dec_deg(&self, jd_ut: f64) -> Result<(f64, f64)>;
    /// Angular radius of the occulting planet, degrees.
    fn occulting_radius_deg(&self, jd_ut: f64) -> Result<f64>;
    /// Angular radius of the occulted body, degrees
    /// ([`OCC_POINT_SOURCE_RADIUS_DEG`] for stars).
    fn target_radius_deg(&self, jd_ut: f64) -> Result<f64>;
}

/// Ephemeris acquisition for the local planetary-occultation search
/// Altitudes carry the release `-0.5`-deg visibility-floor semantics;
/// azimuths travel unmodified from the provider (the serving binding owns
/// the north-zero to south-zero azimuth conversion). The serving-layer binding
/// implements this over the topocentric pipeline; tests
/// implement it over analytic models.
pub trait PlanetOccultLocalEphemeris: PlanetOccultEphemeris {
    /// Apparent altitude and azimuth of the occulting planet, degrees.
    fn occulting_alt_az_deg(&self, jd_ut: f64) -> Result<(f64, f64)>;
    /// Apparent altitude and azimuth of the occulted body, degrees.
    fn target_alt_az_deg(&self, jd_ut: f64) -> Result<(f64, f64)>;
    /// Topocentric angular radii and separation at one epoch, degrees
    fn local_radii_sep_deg(&self, jd_ut: f64) -> Result<(f64, f64, f64)>;
}

// Sealed planetary-occultation servings.

/// Sealed global planetary-occultation serving: implements
/// [`PlanetOccultEphemeris`] over the sealed/native pipeline for one
/// bound occulting planet and one bound target (planet or resolved
/// star). Site-independent (the global walk is geocentric).
pub struct SealedPlanetOccultGlobalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    serving: &'a crate::positions::SealedServing<'a, P>,
    star_serving: &'a crate::stars::SealedStarServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    calc_flags: i32,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
    occulting: i32,
    target: OccultBodyLeg,
}

/// Sealed local planetary-occultation serving: implements
/// [`PlanetOccultEphemeris`] (by delegation to the bound global
/// serving, so geocentric acquisition exists once) and
/// [`PlanetOccultLocalEphemeris`] over the sealed/native pipeline for
/// one bound occulting planet, one bound target and one bound observer
/// triple. The site travels explicitly (explicit-site contract) — there is no
/// global set/restore, so concurrent
/// searches never share observer state.
pub struct SealedPlanetOccultLocalServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    global: SealedPlanetOccultGlobalServing<'a, P>,
    lon_deg_east: f64,
    lat_deg: f64,
    alt_m: f64,
}

// Path samplers, geographic input contracts and remaining event utilities.

/// Central-line surface gate, Earth radii.
pub const PATH_CENTRAL_GAMMA_LIMIT: f64 = 1.5;
/// Central-line nonzero-coordinates gate, degrees.
pub const PATH_CENTRAL_COORD_EPS_DEG: f64 = 0.0001;
/// WGS84 flattening used by the limit samplers.
pub const PATH_LIMIT_FLATTENING: f64 = 1.0 / 298.257223563;
/// Limit-sampler umbral-radius floor, Earth radii.
pub const PATH_LIMIT_UMBRA_FLOOR_ER: f64 = 0.001;
/// Limit-sampler `z`-factor gate.
pub const PATH_LIMIT_SURFACE_GAMMA: f64 = 0.9999;
/// Limit-sampler polar `cos_lat` floor.
pub const PATH_LIMIT_COS_LAT_FLOOR: f64 = 0.001;
/// Minutes per day.
pub const MINUTES_PER_DAY: f64 = 1440.0;
/// Path-width cone radii, km.
pub const PATH_WIDTH_SUN_RADIUS_KM: f64 = 696340.0;
/// Lunar radius used by the path-width cone geometry, km.
pub const PATH_WIDTH_MOON_RADIUS_KM: f64 = 1737.4;
/// Path-width grazing `sin_alt` floor and width cap, km
pub const PATH_WIDTH_SIN_ALT_FLOOR: f64 = 0.05;
/// Path-width cap applied after the grazing `sin_alt` floor, km
pub const PATH_WIDTH_MAX_KM: f64 = 1500.0;
/// How-details local-maximum search: half-range in days, golden iterations, convergence in days
pub const DETAILS_LOCAL_MAX_RANGE_DAYS: f64 = 3.0 / 24.0;
/// How-details local-maximum golden iterations.
pub const DETAILS_LOCAL_MAX_ITERATIONS: usize = 50;
/// How-details local-maximum convergence in days.
pub const DETAILS_LOCAL_MAX_TOL_DAYS: f64 = 1e-8;
/// How-details contact search: half-range in days, bisection iterations
/// separation match in degrees, bracket
/// width in days.
pub const DETAILS_CONTACT_RANGE_DAYS: f64 = 2.5 / 24.0;
/// How-details contact bisection iterations.
pub const DETAILS_CONTACT_ITERATIONS: usize = 60;
/// How-details contact separation match in degrees.
pub const DETAILS_CONTACT_MATCH_DEG: f64 = 1e-7;
/// How-details contact bracket width in days.
pub const DETAILS_CONTACT_TOL_DAYS: f64 = 1e-9;
/// `geopos` length-rule message of `sol_eclipse_when_loc`, carried
/// unchanged from the release.
pub const GEOPOS_SEQUENCE_MESSAGE: &str =
    "geopos must be a sequence of [longitude, latitude, altitude]";
/// `geopos` length-rule message of `_sol_eclipse_how_details_impl`,
/// carried unchanged from the release.
pub const GEOPOS_ELEMENTS_MESSAGE: &str = "geopos must have at least 3 elements: [lon, lat, alt]";

/// Central-line / path-width acquisition: shadow-center `sol_eclipse_where`
/// triple plus Besselian fundamental-plane coordinates.
pub trait EclipsePathEphemeris {
    /// `(retflag, lon_deg, lat_deg)` of the shadow-axis ground point.
    fn path_where(&self, jd_ut: f64) -> Result<(i32, f64, f64)>;
    /// Besselian `(x, y)` in Earth radii at `jd_ut`.
    fn path_besselian_xy(&self, jd_ut: f64) -> Result<(f64, f64)>;
}

/// Limit-sampler acquisition: Besselian `(x, y, d, mu, l2)` with angles in
/// degrees and radii in Earth radii.
pub trait EclipseLimitEphemeris {
    /// `(x_er, y_er, d_deg, mu_deg, l2_er)` at `jd_ut`.
    fn limit_besselian(&self, jd_ut: f64) -> Result<(f64, f64, f64, f64, f64)>;
}

/// Path-width acquisition: shadow-center `sol_eclipse_where` triple plus
/// one topocentric Sun/Moon sample.
pub struct PathWidthSample {
    /// Topocentric Sun/Moon distances in AU.
    pub sun_dist_au: f64,
    /// Topocentric Moon distance in AU.
    pub moon_dist_au: f64,
    /// True Sun altitude in degrees at the calculation site.
    pub sun_alt_true_deg: f64,
    /// Topocentric Sun-Moon angular separation in degrees.
    pub separation_deg: f64,
}

/// Path-width provider (see [`PathWidthSample`]).
pub trait EclipsePathWidthEphemeris {
    /// `(retflag, lon_deg, lat_deg)` of the shadow-axis ground point.
    fn path_width_where(&self, jd_ut: f64) -> Result<(i32, f64, f64)>;
    /// Topocentric Sun/Moon sample at `(lon_deg, lat_deg)`, altitude 0.
    fn path_width_sample(&self, jd_ut: f64, lon_deg: f64, lat_deg: f64) -> Result<PathWidthSample>;
}

/// Solar-duration acquisition: shadow-center `sol_eclipse_where` triple
/// plus local central contacts.
pub trait SolarDurationEphemeris {
    /// `(retflag, lon_deg, lat_deg)` of the shadow-axis ground point.
    fn duration_where(&self, jd_ut: f64) -> Result<(i32, f64, f64)>;
    /// Local `(second_contact_jd, third_contact_jd)` at the site
    /// (`0.0` slots where the contact does not exist).
    fn duration_local_contacts(&self, jd_ut: f64, lat_deg: f64, lon_deg: f64)
    -> Result<(f64, f64)>;
}

/// Completed local planetary event with binary64 contact times and attributes.
#[derive(Debug, Clone)]
pub struct PlanetOccultLocalResult {
    /// Eclipse character flags.
    pub flags: i32,
    /// Ten UT1 Julian-day contact slots; absent contacts are floating zeros.
    pub times: [f64; 10],
    /// Twenty attributes with the ordinary local layout and reserved zeros.
    pub attributes: [f64; 20],
}

/// Local solar-eclipse circumstances at the requested instant and nearby maximum.
///
/// Contact/maximum dates are UT1 Julian days, with zero for unavailable contacts.
/// Angular values are degrees; solar azimuths are south-zero westward and solar
/// altitudes are unrefracted. Magnitudes are solar-diameter fractions and
/// obscuration fields use the producing entry's area-ratio convention, which
/// can exceed one in totality; percentage fields multiply those values by 100.
/// Instantaneous and local-maximum fields are distinct and need not describe
/// the same epoch. Widths are signed kilometres and durations minutes.
#[derive(Debug, Clone)]
pub struct EclipseHowDetails {
    /// Eclipse-type bitmask (`ECL_*`).
    pub eclipse_type: i32,
    /// Visibility/class predicates from the initial type.
    pub is_visible: bool,
    /// Updated from the local maximum (central-only arms).
    pub is_total: bool,
    /// Updated from the local maximum (central-only arms).
    pub is_annular: bool,
    /// Updated from the local maximum (false when central).
    pub is_partial: bool,
    /// Contact times in JD UT (`0.0` where not applicable).
    pub jd_c1: f64,
    /// Second contact (`0.0` unless central).
    pub jd_c2: f64,
    /// Local-maximum time.
    pub jd_max: f64,
    /// Third contact (`0.0` unless central).
    pub jd_c3: f64,
    /// Fourth contact.
    pub jd_c4: f64,
    /// Instant magnitude/obscuration from the initial `sol_eclipse_how`.
    pub magnitude: f64,
    /// Local-maximum magnitude.
    pub max_magnitude: f64,
    /// Instant obscuration from the initial call.
    pub obscuration: f64,
    /// Local-maximum obscuration.
    pub max_obscuration: f64,
    /// Obscuration as a percentage.
    pub obscuration_percent: f64,
    /// Local-maximum obscuration as a percentage.
    pub max_obscuration_percent: f64,
    /// Lunar/solar diameter ratio (updated at the maximum).
    pub ratio: f64,
    /// Core shadow width in km from the initial call.
    pub shadow_width_km: f64,
    /// Position angles in degrees (Moon relative to Sun).
    pub position_angle_c1: f64,
    /// Second-contact position angle (`0.0` unless central).
    pub position_angle_c2: f64,
    /// Third-contact position angle (`0.0` unless central).
    pub position_angle_c3: f64,
    /// Fourth-contact position angle.
    pub position_angle_c4: f64,
    /// True Sun altitude/azimuth at each phase.
    pub sun_alt_c1: f64,
    /// Sun azimuth at first contact.
    pub sun_az_c1: f64,
    /// True Sun altitude at second contact.
    pub sun_alt_c2: f64,
    /// Sun azimuth at second contact.
    pub sun_az_c2: f64,
    /// True Sun altitude at the local maximum.
    pub sun_alt_max: f64,
    /// Sun azimuth at the local maximum.
    pub sun_az_max: f64,
    /// True Sun altitude at third contact.
    pub sun_alt_c3: f64,
    /// Sun azimuth at third contact.
    pub sun_az_c3: f64,
    /// True Sun altitude at fourth contact.
    pub sun_alt_c4: f64,
    /// Sun azimuth at fourth contact.
    pub sun_az_c4: f64,
    /// Partial-phase span in minutes.
    pub duration_partial_minutes: f64,
    /// Central-phase span in minutes.
    pub duration_total_minutes: f64,
    /// Apparent angular radii in degrees.
    pub sun_angular_radius: f64,
    /// Apparent lunar radius in degrees.
    pub moon_angular_radius: f64,
    /// Instant Sun-Moon separation from the initial call.
    pub separation: f64,
}

/// How-details acquisition: the initial `sol_eclipse_how` pair plus the
/// site providers.
pub trait EclipseHowDetailsEphemeris {
    /// `(eclipse_type, attr20)` at the requested instant.
    fn details_initial(&self, jd_ut: f64) -> Result<(i32, [f64; 20])>;
    /// Topocentric Sun-Moon separation in degrees.
    fn details_separation_deg(&self, jd_ut: f64) -> Result<f64>;
    /// True Sun `(altitude, azimuth)` in degrees.
    fn details_sun_alt_az_deg(&self, jd_ut: f64) -> Result<(f64, f64)>;
    /// Apparent `(sun_radius, moon_radius)` in degrees.
    fn details_disc_radii_deg(&self, jd_ut: f64) -> Result<(f64, f64)>;
    /// Equatorial `(sun_ra, sun_dec, moon_ra, moon_dec)` in degrees for
    /// position angles.
    fn details_equatorial_deg(&self, jd_ut: f64) -> Result<(f64, f64, f64, f64)>;
    /// Position angle in degrees, retaining the original angular representation.
    /// The default consumes equatorial degrees; direct native ICRS observers
    /// override it to avoid a radians/degrees round trip.
    fn details_position_angle_deg(&self, jd_ut: f64) -> Result<f64> {
        let (sra, sdec, mra, mdec) = self.details_equatorial_deg(jd_ut)?;
        position_angle_deg(sra, sdec, mra, mdec)
    }
}

/// One solved how-details contact plus its optional tail words: the
/// contact JD always travels with the block; `alt_az`/`position_angle_deg` are `None`
/// when their acquisition degraded. A `None` block means the contact solve itself
/// degraded and the contact stays `0.0`.
struct DetailsContactBlock {
    jd: f64,
    alt_az: Option<(f64, f64)>,
    position_angle_deg: Option<f64>,
}

// Sealed serving binding for longitude searches
use crate::positions::{SealedServing, StateProvider, serving_compute_tt};

/// Sealed longitude serving for one body: UT-to-TT Delta-T plus the
/// sealed apparent ecliptic longitude/speed words.
///
/// `deltat_at` maps a UT1 Julian day to TT-minus-UT1 in days; it is an injected
/// model, not an implicit selection of a particular Delta-T formula.
/// `calc_flags` are retained and SPEED is added at acquisition. Direct component
/// methods preserve those frame/output flags; degree-target search wrappers
/// reject XYZ/RADIANS rather than interpreting their components as degrees.
/// Provider and callback references are borrowed for the serving's lifetime.
/// Construction acquires nothing; consumed feeds own any I/O and errors.
pub struct SealedLongitudeServing<'a, P> {
    serving: &'a SealedServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    body_id: i32,
    calc_flags: i32,
}

/// Station-motion kind reported by release `get_station_type`
/// (`"direct"` / `"retrograde"` / `"stationary_retrograde"` /
/// `"stationary_direct"`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StationMotionKind {
    /// Prograde motion, or any station query that cannot be served
    /// (Sun/Moon identity, provider failure — the release folds both to
    /// `"direct"`).
    Direct,
    /// Retrograde motion (negative longitude speed).
    Retrograde,
    /// Within the near-station band with decreasing speed (about to go
    /// retrograde).
    StationaryRetrograde,
    /// Within the near-station band with non-decreasing speed (about to
    /// resume direct motion).
    StationaryDirect,
}

/// Fixed-instant station summary over the sealed serving (same composition
/// as release `get_station_info`: the nearest station in
/// either direction via [`nearest_station`] — forward search paired with
/// the backward search, nearest by absolute time distance, signed
/// `days_to_station` — plus the ecliptic longitude at the station and
/// the current motion words. Sun/Moon refuse via [`check_station_body`];
/// the release string-error partial dict is a binding-layer shape, not
/// reproduced here — a double search failure surfaces the first error).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StationInfo {
    /// Julian Day of the nearest station (past or future).
    pub jd_station: f64,
    /// Its kind.
    pub kind: StationKind,
    /// Signed days from `start` (negative when the station is past).
    pub days_to_station: f64,
    /// Ecliptic longitude at the station in degrees.
    pub longitude_at_station_deg: f64,
    /// Whether the body is retrograde at `start`.
    pub currently_retrograde: bool,
    /// Longitude speed at `start` in degrees/day.
    pub speed_now_deg_day: f64,
}

// Sealed meridian-transit serving
/// Borrowed position and Delta-T inputs for planetary meridian-transit searches.
///
/// Search epochs are UT1 Julian days. The clock callback maps a UT1 Julian day
/// to TT-minus-UT1 in days for TT position evaluation; retained sidereal offsets
/// are degrees. The serving and callback must outlive this value. Construction
/// and evaluation methods propagate their clock and position-provider errors;
/// this type neither owns an ephemeris image nor acquires process configuration.
pub struct SealedTransitServing<'a, P> {
    serving: &'a SealedServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    body_id: i32,
    calc_flags: i32,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
}

// Sealed rise/set serving

/// Astronomical unit in km.
pub const RISE_SET_AU_KM: f64 = 149597870.700;
/// Solar equatorial radius in AU.
pub const SUN_RADIUS_AU: f64 = 696000.0 / RISE_SET_AU_KM;
/// Lunar mean radius in AU.
pub const MOON_RADIUS_AU: f64 = 1738.15 / RISE_SET_AU_KM;
/// Default tropospheric lapse rate in K/m for the `-100` dip horizon
pub const RISE_SET_DIP_LAPSE_RATE_K_M: f64 = 0.0065;

/// Sealed rise/set request: one `rise_trans_true_hor` call worth of
/// inputs. `rsmi` carries the event kind (`CALC_RISE` / `CALC_SET` take
/// effect here; transit kinds are handled by the transit serving) plus the modifier bits
/// (`BIT_DISC_CENTER` / `BIT_DISC_BOTTOM`, `BIT_NO_REFRACTION`,
/// `BIT_FIXED_DISC_SIZE`, twilight bits). `atpress_mbar == 0.0` means
/// "estimate from the observer altitude" (never "no pressure");
/// `horhgt_deg == -100.0` requests the dipped sea horizon.
pub struct RiseSetRequest {
    /// Search start as a UT Julian Day.
    pub jd_start_ut: f64,
    /// Observer longitude East, degrees.
    pub lon_deg_east: f64,
    /// Observer latitude, degrees.
    pub lat_deg: f64,
    /// Observer altitude above sea level, metres.
    pub height_m: f64,
    /// Atmospheric pressure at the observer, mbar (`0.0` estimates).
    pub atpress_mbar: f64,
    /// Atmospheric temperature at the observer, C (used ).
    pub attemp_c: f64,
    /// Custom horizon altitude, degrees (`-100.0` dips the sea horizon).
    pub horhgt_deg: f64,
    /// Event kind plus modifier bits (see above).
    pub rsmi: i32,
}

// Sealed fixed-star rise/set + transit serving

/// Sealed fixed-star event serving for one catalog star: the v1-resolved
/// [`crate::stars::StarMatch`], the sealed star-catalog serving, the UT-to-TT
/// Delta-T provider, the sealed nutation pair for GAST/LST, and the
/// sidereal-continuity offsets resolved once at construction from the
/// same Delta-T provider. Fields are
/// private with accessors (construction-invariant layer); there is no global state and no I/O
/// on this path.
pub struct SealedStarEventServing<'a, P>
where
    P: crate::positions::StateProvider,
{
    star_serving: &'a crate::stars::SealedStarServing<'a, P>,
    star: crate::stars::StarMatch,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
}

/// Sealed fixed-star rise/set request: one star `rise_trans_true_hor`
/// call worth of inputs. The star is bound in the serving (resolved once
/// at construction); `rsmi` carries the event kind (`CALC_RISE` /
/// `CALC_SET` take effect here; transit kinds belong to the star transit
/// entry) plus the modifier bits (`BIT_NO_REFRACTION`, twilight bits,
/// the Hindu `BIT_GEOCTR_NO_ECL_LAT`). Disc bits are accepted and
/// ignored — fixed stars always use the disc center, like the release.
/// `atpress_mbar == 0.0` means "estimate from the observer altitude";
/// `horhgt_deg == -100.0` requests the dipped sea horizon.
pub struct StarRiseSetRequest {
    /// Search start as a UT Julian Day.
    pub jd_start_ut: f64,
    /// Observer longitude East, degrees.
    pub lon_deg_east: f64,
    /// Observer latitude, degrees.
    pub lat_deg: f64,
    /// Observer altitude above sea level, metres.
    pub height_m: f64,
    /// Atmospheric pressure at the observer, mbar (`0.0` estimates).
    pub atpress_mbar: f64,
    /// Atmospheric temperature at the observer, C (passes through unchanged).
    pub attemp_c: f64,
    /// Custom horizon altitude, degrees (`-100.0` dips the sea horizon).
    pub horhgt_deg: f64,
    /// Event kind plus modifier bits (see above).
    pub rsmi: i32,
}

/// Sealed lunar-point event serving for one out-of-map lunar id: the
/// sealed Pipeline-A serving pair (Moon/Earth geometric states plus the
/// LEB nutation pair), the UT-to-TT Delta-T provider, and the
/// sidereal-continuity offsets resolved once at construction from the
/// same provider.
///
/// `body_id` is the requested id (a south alias resolves to its north
/// source with the antipode applied on every word, the release
/// `calc_ut` interception rule). `calc_flags` carries the caller's
/// sealed ephemeris bits; sidereal, J2000, NONUT, XYZ, radians and
/// helio/bary selections are refused at construction — the served words
/// stay on the tropical ecliptic/equatorial basis served by this item.
pub struct SealedLunarEventServing<'a, P> {
    serving: &'a crate::positions::SealedServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    body_id: i32,
    north_id: i32,
    south: bool,
    calc_flags: i32,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
}

/// Sealed minor-body event serving for one out-of-map table id: the
/// sealed Pipeline-A serving pair (composite core+asteroids+exotics
/// provider plus the LEB nutation pair), the UT-to-TT Delta-T provider, and the
/// sidereal-continuity offsets resolved once at construction from the
/// same provider.
///
/// `body_id` is the requested id (an `AST_OFFSET + N` alias resolves
/// through the release `_remap_ast_offset` rule before the gate).
/// `calc_flags` carries the caller's sealed ephemeris bits; sidereal,
/// J2000, NONUT, XYZ, radians and helio/bary selections are refused at
/// construction — the served words stay on the tropical
/// ecliptic/equatorial basis served by this item.
pub struct SealedMinorEventServing<'a, P> {
    serving: &'a crate::positions::SealedServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    body_id: i32,
    minor_id: i32,
    calc_flags: i32,
    leb_channeled: bool,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
}

// Sealed fictitious/Uranian rise/set + transit serving --

/// Fictitious composition family driving the per-family apparent-place
/// reduction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FictFamily {
    /// Hamburg-School Uranian planets 40-47 (same-jd heliocentric
    /// differencing, Bradley aberration, no target light-time, no PPN).
    Uranian,
    /// Transpluto/Isis 48 (3-iteration geocentric-distance light-time,
    /// first-order aberration, no PPN).
    Transpluto,
    /// Vulcan/Proserpina/predicted planets 50-55/57 (3-iteration
    /// light-time, shared PPN deflection, first-order aberration).
    PredictedHelio,
    /// Selena/White Moon 56 (geocentric symbolic model, first-order
    /// aberration, Δψ tail without a precession step).
    WhiteMoon,
    /// Waldemath 58.
    Waldemath,
}

/// Sealed fictitious-body event serving for one analytical-model id: the
/// sealed Pipeline-A serving pair (core/asteroids/exotics composite
/// provider for the Earth/Sun/deflector vectors plus the LEB nutation
/// pair), the UT-to-TT Delta-T provider, and the sidereal-continuity
/// offsets resolved once at construction.
///
/// `calc_flags` carries the caller's sealed ephemeris bits; sidereal,
/// J2000, NONUT, XYZ, radians and helio/bary selections are refused at
/// construction — the served words stay on the tropical
/// ecliptic/equatorial basis served by this item.
pub struct SealedFictitiousEventServing<'a, P> {
    serving: &'a crate::positions::SealedServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    body_id: i32,
    calc_flags: i32,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
}

// Sealed global-solar-eclipse walk serving

/// Sun/Moon Cartesian pair behind the sealed global-solar walk: two
/// geocentric apparent equatorial AU vectors in `(sun, moon)` order,
/// named once so every walk word shares the pair shape.
type SunMoonXyz = ((f64, f64, f64), (f64, f64, f64));

/// Sealed global-solar-eclipse serving: the sealed Pipeline-A Sun/Moon
/// acquisitions plus Delta-T behind the [`SolarGlobalEphemeris`] walk
/// contract.
///
/// `calc_flags` carries the caller's sealed ephemeris bits (normally
/// `FLG_SWIEPH`); the per-word shaping bits (`FLG_SPEED`,
/// `FLG_EQUATORIAL`, `FLG_J2000`, `FLG_TOPOCTR`) are added per
/// acquisition exactly at the release call sites above, so the served
/// words stay on the verified tropical basis.
pub struct SealedSolarGlobalServing<'a, P> {
    serving: &'a crate::positions::SealedServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    calc_flags: i32,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
}

// Global lunar-eclipse serving over owned Sun/Moon/Earth data.
// UT1 is converted through the bound Delta-T provider. Opposition seeds,
// shadow depth, contact roots and type filters share those resources.
// Unbracketed contacts retain their documented zero sentinel; provider
// failures propagate instead of fabricating an event.

/// Sealed global-lunar-eclipse serving: the sealed Pipeline-A Sun/Moon
/// acquisitions plus Delta-T behind the [`LunarGlobalEphemeris`] walk
/// contract.
///
/// `calc_flags` carries the caller's sealed ephemeris bits (normally
/// `FLG_SWIEPH`); the per-word shaping bits (`FLG_SPEED`,
/// `FLG_EQUATORIAL`) are added per acquisition exactly at the release
/// call sites above, so the served words stay on the verified tropical
/// basis.
pub struct SealedLunarGlobalServing<'a, P> {
    serving: &'a crate::positions::SealedServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    calc_flags: i32,
}

// Sealed local-solar-eclipse walk serving

/// Sealed local-solar-eclipse serving: the sealed Pipeline-A Sun/Moon
/// acquisitions, Delta-T and the observer triple behind the
/// [`SolarLocalEphemeris`] walk contract (plus the
/// [`SolarGlobalEphemeris`] facet feeding its candidate lunations).
///
/// `calc_flags` carries the caller's sealed ephemeris bits (normally
/// `FLG_SWIEPH`); the per-word shaping bits (`FLG_SPEED`,
/// `FLG_EQUATORIAL`, `FLG_J2000`, `FLG_TOPOCTR`) are added per
/// acquisition exactly at the release call sites above, so the served
/// words stay on the verified tropical basis. The observer triple is
/// `(longitude_deg_east, latitude_deg, altitude_m)` in release
/// `geopos` order (longitude first).
pub struct SealedSolarLocalServing<'a, P> {
    serving: &'a crate::positions::SealedServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    calc_flags: i32,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
    lon_deg_east: f64,
    lat_deg: f64,
    alt_m: f64,
}

// Sealed local-lunar-eclipse walk serving

/// Sealed local-lunar-eclipse serving: the sealed Pipeline-A Sun/Moon
/// acquisitions, Delta-T and the observer triple behind the
/// [`LunarLocalEphemeris`] walk contract (plus the
/// [`LunarGlobalEphemeris`] facet feeding its candidate lunations).
///
/// `calc_flags` carries the caller's sealed ephemeris bits (normally
/// `FLG_SWIEPH`); the per-word shaping bits (`FLG_SPEED`,
/// `FLG_EQUATORIAL`, `FLG_TOPOCTR`) are added per acquisition exactly
/// at the release call sites above, so the served words stay on the
/// verified tropical basis. The observer triple is
/// `(longitude_deg_east, latitude_deg, altitude_m)` in release
/// `geopos` order (longitude first).
pub struct SealedLunarLocalServing<'a, P> {
    serving: &'a crate::positions::SealedServing<'a, P>,
    deltat_at: &'a dyn Fn(f64) -> Result<f64>,
    calc_flags: i32,
    offset_lo_deg: f64,
    offset_hi_deg: f64,
    lon_deg_east: f64,
    lat_deg: f64,
    alt_m: f64,
}

// Private implementation boundaries; public logical paths remain unchanged.
mod search;
pub use search::bisection;
pub use search::brent;
pub use search::check_bracket;
pub use search::newton;
pub use search::scan_first_bracket;
mod crossings;
pub use crossings::adaptive_max_iterations;
pub use crossings::backward_start_exclusion_days;
pub use crossings::brent_crossing;
pub use crossings::cross_max_range_days;
pub use crossings::find_crossing_bracket;
pub use crossings::first_forward_crossing_adaptive;
pub use crossings::first_forward_helio_crossing;
pub use crossings::forward_first_crossing;
pub use crossings::guarded_first_crossing;
pub use crossings::helio_brent;
pub use crossings::helio_cross_bracketed;
pub use crossings::helio_crossing;
pub use crossings::helio_crossing_serving;
pub use crossings::helio_speed_known;
pub use crossings::helio_typical_speed_deg_day;
pub use crossings::is_near_station;
use crossings::mod360_deg;
pub use crossings::moon_crossing;
pub use crossings::mooncross_node_ut_serving;
pub use crossings::mooncross_ut_serving;
pub use crossings::nearest_past_helio_crossing;
pub use crossings::needs_first_crossing_guard;
pub use crossings::newton_crossing;
pub use crossings::newton_residual_deg;
pub use crossings::normalize_target_deg;
pub use crossings::planet_crossing;
pub use crossings::retro_arc_deg;
pub use crossings::seed_diff_deg;
pub use crossings::slow_max_rate_deg_day;
pub use crossings::sun_crossing;
pub use crossings::synodic_days;
pub use crossings::typical_speed_deg_day;
pub use crossings::wrapped_residual_deg;
mod stations;
pub use stations::brent_station;
pub use stations::check_station_body;
#[cfg(test)]
use stations::checked_speed;
pub use stations::classify_station;
pub use stations::find_previous_station;
pub use stations::find_previous_station_ut_serving;
pub use stations::find_station;
pub use stations::find_station_bracket;
pub use stations::find_station_ut_serving;
pub use stations::is_retrograde_speed;
pub use stations::is_retrograde_ut_serving;
pub use stations::moon_node_crossing;
pub use stations::nearest_station;
pub use stations::nearest_station_ut_serving;
pub use stations::next_retrograde;
pub use stations::next_retrograde_ut_serving;
pub use stations::station_central_velocity;
pub use stations::station_info_ut_serving;
pub use stations::station_matches_request;
pub use stations::station_type_ut_serving;
pub use stations::velocity_at_station_ut_serving;
mod rise_set;
pub use rise_set::HorizonFrame;
pub use rise_set::azalt_deg;
pub use rise_set::azalt_deg_with_frame;
pub use rise_set::azalt_rev_deg;
pub use rise_set::azalt_rev_deg_with_frame;
pub use rise_set::body_disc_radius_au;
pub use rise_set::calculate_rise_set;
pub use rise_set::calculate_transit;
pub use rise_set::disc_altitude_term_deg;
pub use rise_set::eclipse_disc_angular_radii_deg;
pub use rise_set::event_height_deg;
pub use rise_set::fictitious_rise_set_ut_serving;
pub use rise_set::fictitious_transit_ut_serving;
pub use rise_set::horizon_forward_geometry_radians;
pub use rise_set::horizon_reverse_geometry_radians;
pub use rise_set::lunar_rise_set_ut_serving;
pub use rise_set::lunar_transit_ut_serving;
pub use rise_set::make_tret;
pub use rise_set::minor_rise_set_ut_serving;
pub use rise_set::minor_transit_ut_serving;
pub use rise_set::pressure_from_altitude_m;
pub(crate) use rise_set::refine_horizon_crossing;
pub use rise_set::resolve_event_kind;
pub use rise_set::resolve_horizon_altitude_deg;
pub use rise_set::rise_set_ut_serving;
pub use rise_set::select_disc_point;
pub use rise_set::star_rise_set_ut_serving;
pub use rise_set::star_transit_ut_serving;
pub use rise_set::transit_ut_serving;
pub use rise_set::true_to_apparent_altitude;
pub use rise_set::twilight_override;
pub use rise_set::validate_geopos;
mod geometry;
use geometry::angular_separation_ecliptic;
pub use geometry::besselian_contact_residuals;
pub(crate) use geometry::besselian_core_lazy;
pub use geometry::besselian_core_via;
pub use geometry::besselian_elements_sealed;
pub use geometry::besselian_rates_via;
pub use geometry::calc_besselian_d_sealed;
pub use geometry::calc_besselian_dd_dt_sealed;
pub use geometry::calc_besselian_dl1_dt_sealed;
pub use geometry::calc_besselian_dl2_dt_sealed;
pub use geometry::calc_besselian_dmu_dt_sealed;
pub use geometry::calc_besselian_dx_dt_sealed;
pub use geometry::calc_besselian_dy_dt_sealed;
pub use geometry::calc_besselian_l1_sealed;
pub use geometry::calc_besselian_l2_sealed;
pub use geometry::calc_besselian_mu_sealed;
pub use geometry::calc_besselian_x_sealed;
pub use geometry::calc_besselian_y_sealed;
pub use geometry::calculate_obscuration;
pub use geometry::calculate_solar_phases_besselian;
pub use geometry::details_max_obscuration;
use geometry::dist3;
pub use geometry::find_contact_time_gamma;
pub use geometry::gamma_from_vectors;
pub use geometry::interpolate_besselian_elements;
pub use geometry::planet_angular_separation_deg;
pub use geometry::safe_acos;
pub use geometry::safe_sqrt;
pub use geometry::sol_eclipse_obscuration_at_loc_sealed;
pub use geometry::sol_eclipse_obscuration_at_loc_via;
pub(crate) use geometry::solar_gamma_reader_vectors_at_tt;
pub(crate) use geometry::solar_gamma_vectors_at_tt;
pub use geometry::solve_besselian_contact;
pub use geometry::vincenty_separation_deg;
mod solar_eclipses;
mod solar_path;
pub(crate) use solar_path::{
    SolarPathWidthSource, solar_limit_source_via, solar_path_width_source_via,
};
mod solar_local_search;
use solar_eclipses::altaz_of_topo_ecliptic_impl;
pub use solar_eclipses::bisect_sign_change;
#[cfg(test)]
use solar_eclipses::bradley_aberration_xyz;
pub use solar_eclipses::calc_local_eclipse_max_time_via;
pub use solar_eclipses::calc_planet_angular_radius_deg;
pub use solar_eclipses::calculate_eclipse_magnitude;
pub use solar_eclipses::calculate_eclipse_magnitude_with_limit;
pub use solar_eclipses::calculate_local_eclipse_phases_via;
use solar_eclipses::check_angular_serving_flags;
use solar_eclipses::check_full_moon_for_eclipse;
pub use solar_eclipses::classify_solar_eclipse;
#[cfg(test)]
pub(crate) use solar_eclipses::classify_solar_eclipse_consumed;
pub use solar_eclipses::coerce_backwards;
pub use solar_eclipses::collapse_point_source_contacts;
use solar_eclipses::compensated_sum3;
use solar_eclipses::contact_provider_value;
pub use solar_eclipses::cotrans;
pub use solar_eclipses::cotrans_radians;
pub use solar_eclipses::cross_ut_serving;
pub use solar_eclipses::details_contact_time;
pub use solar_eclipses::details_max_magnitude;
pub use solar_eclipses::determine_local_eclipse_type;
pub use solar_eclipses::ecl_eph_flags;
#[cfg(test)]
use solar_eclipses::ecl_xyz_to_spherical;
pub use solar_eclipses::eclipse_central_line_via;
pub use solar_eclipses::eclipse_how_details_via;
pub use solar_eclipses::eclipse_limit_point_deg;
pub use solar_eclipses::eclipse_max_time_via;
pub use solar_eclipses::eclipse_northern_limit_via;
pub use solar_eclipses::eclipse_path_width_via;
pub use solar_eclipses::eclipse_phase_duration_min;
pub use solar_eclipses::eclipse_sampling_step_days;
pub use solar_eclipses::eclipse_southern_limit_via;
pub use solar_eclipses::eclipse_where_core_via;
pub use solar_eclipses::find_local_contact_time;
pub use solar_eclipses::find_local_maximum_via;
pub use solar_eclipses::find_next_full_moon;
pub use solar_eclipses::find_next_new_moon;
pub use solar_eclipses::find_previous_full_moon;
pub use solar_eclipses::find_previous_new_moon;
#[cfg(test)]
use solar_eclipses::first_order_aberration_xyz;
use solar_eclipses::fold_180_deg;
pub use solar_eclipses::fold_pseudo_body_to_sun;
pub use solar_eclipses::golden_section_min;
pub(crate) use solar_eclipses::golden_section_min_source;
pub use solar_eclipses::golden_section_min_traced;
pub use solar_eclipses::inex_number;
pub use solar_eclipses::is_ephemeris_boundary;
pub use solar_eclipses::is_fictitious_event_body;
pub use solar_eclipses::is_hybrid_solar_eclipse;
pub use solar_eclipses::is_minor_event_body;
pub use solar_eclipses::is_near_miss_eclipse;
pub use solar_eclipses::is_near_miss_eclipse_with_limit;
pub use solar_eclipses::is_shallow_eclipse;
pub use solar_eclipses::kernel_jd_bounds;
pub use solar_eclipses::moon_node_distance_deg;
pub use solar_eclipses::normalize_path_longitude_deg;
pub use solar_eclipses::occ_body_radius_au;
pub use solar_eclipses::parse_eclipse_kind;
pub use solar_eclipses::parse_geopos;
pub use solar_eclipses::parse_search_direction;
pub use solar_eclipses::penumbra_limit_er;
pub use solar_eclipses::position_angle_deg;
pub use solar_eclipses::refine_solar_eclipse_maximum;
#[cfg(test)]
use solar_eclipses::rotate_vector_to_of_date;
#[cfg(test)]
use solar_eclipses::round_half_even;
pub use solar_eclipses::saros_info;
pub use solar_eclipses::saros_number;
pub use solar_eclipses::seek_moon_conjunction;
pub use solar_eclipses::selenocentric_depth_deg;
pub use solar_eclipses::sol_contact_c1_sealed;
pub use solar_eclipses::sol_contact_c2_sealed;
pub use solar_eclipses::sol_contact_c3_sealed;
pub use solar_eclipses::sol_contact_c4_sealed;
pub use solar_eclipses::sol_eclipse_central_line_sealed;
pub use solar_eclipses::sol_eclipse_duration_min_sealed;
pub use solar_eclipses::sol_eclipse_how_details_sealed;
pub use solar_eclipses::sol_eclipse_how_parts;
pub use solar_eclipses::sol_eclipse_magnitude_at_loc_sealed;
pub use solar_eclipses::sol_eclipse_magnitude_at_loc_via;
pub use solar_eclipses::sol_eclipse_northern_limit_sealed;
pub use solar_eclipses::sol_eclipse_path_width_sealed;
pub use solar_eclipses::sol_eclipse_southern_limit_sealed;
pub use solar_eclipses::sol_eclipse_when_glob_sealed;
pub(crate) use solar_eclipses::sol_eclipse_when_glob_source_via;
pub use solar_eclipses::sol_eclipse_when_glob_via;
pub use solar_eclipses::sol_eclipse_when_loc_sealed;
pub use solar_eclipses::sol_eclipse_when_loc_via;
pub use solar_eclipses::sol_eclipse_where_parts;
pub use solar_eclipses::sol_eclipse_where_sealed;
pub use solar_eclipses::sol_glob_accepts;
pub use solar_eclipses::sol_glob_reject_impossible;
pub use solar_eclipses::sol_global_max_time_sealed;
pub use solar_eclipses::sol_how_core_via;
pub use solar_eclipses::sol_how_sealed;
pub use solar_eclipses::sol_local_max_time_sealed;
pub use solar_eclipses::solar_contact_c1_via;
pub use solar_eclipses::solar_contact_c2_via;
pub use solar_eclipses::solar_contact_c3_via;
pub use solar_eclipses::solar_contact_c4_via;
pub use solar_eclipses::solar_eclipse_duration_min_via;
pub use solar_eclipses::solcross_ut_serving;
use solar_eclipses::solve_contact_opt;
pub use solar_eclipses::strip_one_try_bit;
pub use solar_eclipses::umbra_limit_er;
pub use solar_eclipses::validate_contact_time;
pub use solar_eclipses::wrap_ra_gap_deg;
pub(crate) use solar_eclipses::{eclipse_where_core_lazy, sol_how_core_lazy};
pub use solar_local_search::sol_eclipse_when_loc_source_via;
mod lunar_eclipses;
mod lunar_global_search;
mod lunar_local_search;
pub use lunar_eclipses::classify_lunar_eclipse;
pub use lunar_eclipses::find_next_syzygy;
pub use lunar_eclipses::find_previous_syzygy;
pub use lunar_eclipses::is_lunar_event_body;
pub use lunar_eclipses::lun_contact_p1_sealed;
pub use lunar_eclipses::lun_contact_p4_sealed;
pub use lunar_eclipses::lun_contact_u1_sealed;
pub use lunar_eclipses::lun_contact_u2_sealed;
pub use lunar_eclipses::lun_contact_u3_sealed;
pub use lunar_eclipses::lun_contact_u4_sealed;
pub use lunar_eclipses::lun_contact_via_core;
pub use lunar_eclipses::lun_eclipse_how_parts;
pub use lunar_eclipses::lun_eclipse_max_time_sealed;
pub use lunar_eclipses::lun_eclipse_phase_times_via;
pub use lunar_eclipses::lun_eclipse_when_glob_sealed;
pub use lunar_eclipses::lun_eclipse_when_loc_sealed;
pub use lunar_eclipses::lun_eclipse_when_loc_via;
pub use lunar_eclipses::lun_eclipse_when_via;
pub use lunar_eclipses::lun_gamma_ecliptic_via;
pub use lunar_eclipses::lun_gamma_sealed;
pub use lunar_eclipses::lun_gamma_via;
pub use lunar_eclipses::lun_how_core_via;
pub use lunar_eclipses::lun_how_sealed;
pub use lunar_eclipses::lun_penumbral_magnitude_sealed;
pub use lunar_eclipses::lun_penumbral_magnitude_via;
pub use lunar_eclipses::lun_total_duration_sealed;
pub use lunar_eclipses::lun_umbral_duration_sealed;
pub use lunar_eclipses::lun_umbral_magnitude_sealed;
pub use lunar_eclipses::lun_umbral_magnitude_via;
pub use lunar_eclipses::lunar_penumbral_contact_p1_via;
pub use lunar_eclipses::lunar_penumbral_contact_p4_via;
pub use lunar_eclipses::lunar_total_duration_min_via;
pub use lunar_eclipses::lunar_umbral_contact_u1_via;
pub use lunar_eclipses::lunar_umbral_contact_u2_via;
pub use lunar_eclipses::lunar_umbral_contact_u3_via;
pub use lunar_eclipses::lunar_umbral_contact_u4_via;
pub use lunar_eclipses::lunar_umbral_duration_min_via;
pub use lunar_eclipses::refine_lunar_eclipse_maximum;
pub use lunar_global_search::lun_eclipse_when_source_via;
pub use lunar_global_search::normalize_lunar_search_mask;
pub use lunar_local_search::{LunarLocalSearchResult, lun_eclipse_when_loc_source_via};
mod occultations;
pub use occultations::check_planet_occultation_via;
pub use occultations::classify_occultation_phase;
pub use occultations::classify_planet_occultation;
pub use occultations::compose_occult_where_outcome;
pub use occultations::conjunction_latitude_admits_occultation;
pub use occultations::finish_occult_local_attr;
pub use occultations::is_occultation_in_progress;
pub use occultations::lun_occult_when_glob_sealed;
pub use occultations::lun_occult_when_glob_via;
pub use occultations::lun_occult_when_loc_sealed;
pub use occultations::lun_occult_when_loc_via;
pub use occultations::normalize_occultation_filter;
pub use occultations::occult_global_central_residual_km;
pub use occultations::occult_global_outer_residual_km;
pub use occultations::occult_global_total_residual_km;
pub use occultations::occultation_inner_residual_deg;
pub use occultations::occultation_outer_residual_deg;
pub use occultations::occultation_radii_deg;
pub use occultations::planet_occult_apparent_altitude_deg;
pub use occultations::planet_occult_contact_times_via;
pub use occultations::planet_occult_minimum_via;
pub use occultations::planet_occult_separation_via;
pub use occultations::planet_occult_threshold_met;
pub use occultations::planet_occult_when_glob_sealed;
pub use occultations::planet_occult_when_glob_via;
pub use occultations::planet_occult_when_loc_sealed;
pub use occultations::planet_occult_when_loc_via;
pub use occultations::reject_moon_self_occultation;
pub use occultations::star_never_occulted;
pub use occultations::validate_planet_occult_request;
pub use occultations::{planet_occult_local_attr, planet_occult_when_loc_result_via};

#[cfg(test)]
#[path = "tests/evt13_b28_occult_serving_tests.rs"]
mod evt13_b28_occult_serving_tests;

#[cfg(test)]
#[path = "tests/evt13_b29_planet_occult_serving_tests.rs"]
mod evt13_b29_planet_occult_serving_tests;

#[cfg(test)]
#[path = "tests/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "tests/evt13_serving_tests.rs"]
mod evt13_serving_tests;

#[cfg(test)]
#[path = "tests/evt14_slice20_solcontact_tests.rs"]
mod evt14_slice20_solcontact_tests;

#[cfg(test)]
#[path = "tests/evt14_slice9_solpath_tests.rs"]
mod evt14_slice9_solpath_tests;

#[cfg(test)]
#[path = "tests/evt14_slice5a_solglob_tests.rs"]
mod evt14_slice5a_solglob_tests;

#[cfg(test)]
#[path = "tests/evt14_slice25_lunmaxtime_tests.rs"]
mod evt14_slice25_lunmaxtime_tests;

#[cfg(test)]
#[path = "tests/evt14_slice17_luncontact_tests.rs"]
mod evt14_slice17_luncontact_tests;

#[cfg(test)]
#[path = "tests/evt14_slice18_lunhow_tests.rs"]
mod evt14_slice18_lunhow_tests;

#[cfg(test)]
#[path = "tests/evt14_slice19_solhow_tests.rs"]
mod evt14_slice19_solhow_tests;

#[cfg(test)]
#[path = "tests/evt14_slice5b_lunglob_tests.rs"]
mod evt14_slice5b_lunglob_tests;

#[cfg(test)]
#[path = "tests/evt14_slice21_solmagloc_tests.rs"]
mod evt14_slice21_solmagloc_tests;

#[cfg(test)]
#[path = "tests/evt14_slice22_besselian_tests.rs"]
mod evt14_slice22_besselian_tests;

#[cfg(test)]
#[path = "tests/evt14_slice26_solwhere_tests.rs"]
mod evt14_slice26_solwhere_tests;

#[cfg(test)]
#[path = "tests/evt14_slice23_maxtime_tests.rs"]
mod evt14_slice23_maxtime_tests;

// Sealed lunar-occultation where serving

#[cfg(test)]
#[path = "tests/evt14_slice24_occwhere_tests.rs"]
mod evt14_slice24_occwhere_tests;

#[cfg(test)]
#[path = "tests/evt14_slice10_sealed_solwidth_tests.rs"]
mod evt14_slice10_sealed_solwidth_tests;

#[cfg(test)]
#[path = "tests/evt14_slice6_solloc_tests.rs"]
mod evt14_slice6_solloc_tests;

#[cfg(test)]
#[path = "tests/evt14_slice7_lunloc_tests.rs"]
mod evt14_slice7_lunloc_tests;

// Lunar-side contact contact invariants

#[cfg(test)]
#[path = "tests/evt14_lunpar_contact_tests.rs"]
mod evt14_lunpar_contact_tests;

pub(crate) use solar_path::{solar_path_width_captured_source_via, solar_path_width_geometry};

mod numeric_occultation;
pub(crate) use numeric_occultation::{
    NumericOccultWhereSource, numeric_occult_how_via, numeric_occult_where_via,
};

mod numeric_occult_global;
pub(crate) use numeric_occult_global::numeric_occult_separation_source;
pub(crate) use occultations::lun_occult_when_glob_source_via;
pub(crate) use occultations::lun_occult_when_glob_star_source_via;

pub(crate) use geometry::vincenty_separation_deg as numeric_occult_local_separation;
pub(crate) use occultations::lun_occult_when_loc_source_via;
pub(crate) use occultations::lun_occult_when_loc_star_source_via;
