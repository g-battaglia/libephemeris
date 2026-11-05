//! Typed failures and non-fatal advisories for the native numerical engine.
//!
//! Errors own their data and can cross worker boundaries. Numerical fields use
//! the units documented on each variant. [`Error::category`] groups failures
//! for clients; [`std::error::Error::source`] exposes an actual native cause.
//! Diagnostic wording is not a machine-readable interface. Match variants or
//! categories instead of parsing [`std::fmt::Display`] text. [`Error`] is
//! non-exhaustive: downstream matches need a catch-all arm.
//!
//! A returned `Err` is distinct from a panic. Caller-supplied callbacks can
//! panic, and their owning API documents whether worker panics propagate;
//! these failures are not automatically converted into an [`Error`].

use std::fmt;

/// A successful value or an owned native [`Error`].
///
/// Use `?` to propagate failures from a function returning this type. An API
/// returning `Result<Option<T>>` uses `Ok(None)` for its documented absence or
/// uncovered-data case, not for every possible failure.
pub type Result<T> = std::result::Result<T, Error>;

/// Native input, data, configuration and calculation failures.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum Error {
    /// A Chebyshev coefficient slice was empty.
    EmptyCoefficients,
    /// An argument was NaN or infinite.
    NonFiniteArgument {
        /// Rejected NaN or infinite value, in the calling operation's units.
        value: f64,
    },
    /// A normalized Chebyshev argument was outside `[-1, 1]`.
    NormalizedTimeOutOfRange {
        /// Dimensionless Chebyshev coordinate that lies outside `[-1, 1]`.
        x: f64,
    },
    /// The requested derivative order exceeds the supported maximum.
    InvalidDerivativeOrder {
        /// Derivative order requested by the caller.
        order: usize,
        /// Largest derivative order supported by this evaluator.
        max: usize,
    },
    /// Ephemeris segment bounds were degenerate or non-finite.
    InvalidSegment {
        /// Declared segment start, in the evaluator's time-coordinate units.
        begin: f64,
        /// Declared segment end, in the same units as `begin`.
        end: f64,
    },
    /// The evaluation epoch was outside its segment, in the caller's time units.
    EpochOutOfSegment {
        /// Requested epoch, in the same time-coordinate units as the segment bounds.
        t: f64,
        /// Start of the segment assigned to the evaluation.
        begin: f64,
        /// End of that segment, in the same units as `t`.
        end: f64,
    },
    /// The zero vector has no normalization.
    ZeroVectorNormalize,
    /// The requested 3×3 matrix inverse does not exist.
    SingularMatrix,
    /// The zero vector has no spherical direction.
    ZeroVectorToSpherical,
    /// Integer arithmetic or an integer-valued conversion exceeded its range.
    ///
    /// The destination depends on the operation: exact rational kernels use
    /// `i128`, while calendar/time fields and offsets can use `i32` or `i64`.
    IntegerOverflow,
    /// A rational denominator was zero.
    RationalDivisionByZero,
    /// An interval's lower bound exceeds its upper bound.
    IntervalInvertedBounds,
    /// An interval divisor contains zero.
    IntervalSpansZero,
    /// An interval operation was outside its mathematical domain.
    IntervalDomainViolation,
    /// The calendar flag does not identify the Julian or Gregorian calendar.
    InvalidCalendar {
        /// Unrecognized integer calendar selector supplied by the caller.
        value: i32,
    },
    /// The calendar code is neither `j` nor `g` (case-insensitive).
    InvalidCalendarCode {
        /// Unrecognized character supplied instead of a Julian or Gregorian selector.
        code: char,
    },
    /// The calendar date is invalid.
    InvalidDate {
        /// Calendar year supplied for the rejected date.
        year: i32,
        /// Calendar month supplied for the rejected date.
        month: i32,
        /// Day of the month supplied for the rejected date.
        day: i32,
    },
    /// The clock time is invalid; seconds permit an actual leap-second label.
    InvalidTime {
        /// Hour component of the rejected clock time.
        hour: i32,
        /// Minute component of the rejected clock time.
        minute: i32,
        /// Second component, including any fractional seconds, of the rejected time.
        second: f64,
    },
    /// The date has no leap second matching the supplied `23:59:60.x` label.
    NoLeapSecond {
        /// Hour component accompanying the unrecognized leap-second label.
        hour: i32,
        /// Minute component accompanying the unrecognized leap-second label.
        minute: i32,
        /// Supplied leap-second label in seconds, including its fractional part.
        second: f64,
    },
    /// An observer coordinate is outside its domain (angles in degrees).
    InvalidObserver {
        /// Name of the observer or coordinate field whose domain check failed.
        field: &'static str,
        /// Rejected value, in the units associated with `field`.
        value: f64,
    },
    /// A segment index cannot be determined without arithmetic overflow.
    SegmentOwnershipUnresolved {
        /// Epoch whose grid index could not be established, in the grid's time units.
        t: f64,
        /// Grid origin, in the same time-coordinate units as `t`.
        begin: f64,
        /// Segment duration, in the same units as `t` and `begin`.
        length: f64,
    },
    /// The segment derivative scale is not finite.
    DerivativeScaleOverflow {
        /// Nonfinite multiplier used to convert the segment derivative to time units.
        scale: f64,
    },
    /// The requested LEB file could not be opened.
    LebFileNotFound {
        /// Requested LEB path retained for the diagnostic.
        path: String,
    },
    /// The file identification bytes do not identify a supported LEB format.
    LebUnknownFormat {
        /// Four identification bytes read from the file header.
        magic: [u8; 4],
    },
    /// The LEB format version is unsupported.
    LebUnsupportedVersion {
        /// Format version declared by the LEB header.
        found: u32,
    },
    /// The LEB compression scheme is unsupported.
    LebUnsupportedCompression {
        /// Compression identifier declared by the LEB2 header.
        found: u32,
    },
    /// LEB metadata or offsets are invalid; `offset` is a byte offset.
    LebCorrupt {
        /// Description of the invalid metadata, shape or byte range.
        what: String,
        /// Byte offset associated with the structural failure in the LEB image.
        offset: usize,
    },
    /// A declared resource limit was exceeded before the associated work.
    /// Sizes are bytes unless `resource` explicitly names a count.
    ResourceLimit {
        /// Resource label; it identifies whether the counts are bytes or another unit.
        resource: &'static str,
        /// Amount required by the refused operation, in the resource's units.
        required: usize,
        /// Configured maximum, in the same units as `required`.
        limit: usize,
    },
    /// Indexed storage cannot provide a borrowed whole-file byte slice.
    LebBytesUnavailable,
    /// The requested LEB reader or its resource configuration is invalid.
    InvalidLebReader {
        /// Rejected reader selection or resource-configuration value.
        value: String,
    },
    /// The requested body is absent from the LEB image.
    LebBodyNotFound {
        /// Requested body or channel identifier absent from the image.
        body_id: i32,
    },
    /// The requested catalog index is absent from the LEB image.
    LebStarNotFound {
        /// Requested stored-star identifier absent from the image.
        star_id: i32,
    },
    /// The requested operation requires an open LEB reader.
    LebReaderClosed,
    /// The Julian day is outside the stored coverage for this body.
    LebOutOfRange {
        /// Stored body or channel identifier associated with the coverage failure.
        body_id: i32,
        /// Requested Julian day, in the stored channel's time scale.
        jd: f64,
        /// First covered Julian day, in the same time scale as `jd`.
        start: f64,
        /// Last covered Julian day, in the same time scale as `jd`.
        end: f64,
    },
    /// A required asset identity is absent from the data manifest.
    DataManifestMissing {
        /// Asset name whose expected identity was requested from the manifest.
        name: String,
    },
    /// The selected Delta-T model is unavailable in this configuration.
    DeltaTModelUnavailable {
        /// Selected Delta-T model that cannot be served in this configuration.
        model: &'static str,
    },
    /// The requested Delta-T model name is unsupported.
    InvalidDeltaTModel {
        /// Unrecognized Delta-T model name supplied to the selector.
        value: String,
    },
    /// The requested nutation model name is unsupported.
    InvalidNutationModel {
        /// Unrecognized nutation model name supplied to the selector.
        value: String,
    },
    /// The requested ephemeris mode is unsupported.
    InvalidCalcMode {
        /// Rejected ephemeris calculation-mode value.
        value: String,
        /// Setting or entry point from which the mode value was read.
        source: &'static str,
    },
    /// The requested data tier is unsupported.
    InvalidTier {
        /// Rejected precision-tier name.
        value: String,
    },
    /// The network policy is invalid; `source` identifies its setting.
    InvalidNetworkPolicy {
        /// Rejected network-policy value.
        value: String,
        /// Setting or entry point from which the policy value was read.
        source: &'static str,
    },
    /// The configured network policy forbids the named operation.
    NetworkSealed {
        /// Description of the network operation refused by the configured policy.
        purpose: String,
    },
    /// The HTTP timeout in seconds is not finite and positive.
    InvalidTimeout {
        /// Rejected HTTP timeout, in seconds.
        seconds: f64,
    },
    /// The HTTP response-size limit in bytes is zero.
    InvalidBodyLimit {
        /// Rejected maximum HTTP response size, in bytes.
        limit_bytes: usize,
    },
    /// SPK coverage padding in days is not finite and nonnegative.
    InvalidPadding {
        /// Name of the request or setting whose coverage padding was checked.
        entry: &'static str,
        /// Rejected SPK coverage padding, in days.
        value: f64,
    },
    /// The request URL has an unsupported scheme or no host.
    HttpInvalidUrl {
        /// URL whose scheme or host failed validation.
        url: String,
    },
    /// The named HTTP request exceeded its timeout in seconds.
    HttpTimeout {
        /// Description of the HTTP operation that timed out.
        purpose: String,
        /// Timeout applied to that operation, in seconds.
        seconds: f64,
    },
    /// The named request received an unsuccessful HTTP status.
    HttpStatus {
        /// Unsuccessful HTTP response status code.
        status: u16,
        /// Description of the request that received the status.
        purpose: String,
    },
    /// The named request failed at its connection or protocol layer.
    HttpTransport {
        /// Description of the request whose transport failed.
        purpose: String,
        /// Connection or protocol failure diagnostic supplied by the transport.
        detail: String,
    },
    /// The response exceeded the configured byte limit.
    HttpResponseTooLarge {
        /// Description of the request whose response exceeded the limit.
        purpose: String,
        /// Maximum permitted response length, in bytes.
        limit_bytes: usize,
    },
    /// The selected ephemeris filename or loader setting is invalid.
    EphemerisLoaderValue {
        /// Diagnostic explaining the invalid filename or loader selection.
        message: String,
    },
    /// A local loader path failed an OS operation, retaining its native error metadata.
    EphemerisLoaderIo {
        /// Local path on which the loader's OS operation failed.
        path: String,
        /// Portable classification of the underlying I/O error.
        kind: std::io::ErrorKind,
        /// Native OS error number when supplied by the I/O error; otherwise `None`.
        errno: Option<i32>,
        /// Diagnostic retained from the underlying I/O failure.
        message: String,
    },
    /// Horizons refused the requested command.
    HorizonsApiError {
        /// Horizons target command associated with the rejected request.
        command: String,
        /// Rejection diagnostic returned by Horizons.
        message: String,
    },
    /// A Horizons response failed structural validation.
    HorizonsMalformed {
        /// Horizons target command whose response was being validated.
        command: String,
        /// Structural requirement that the response failed to satisfy.
        reason: &'static str,
    },
    /// A Horizons command exhausted its request budget.
    HorizonsRequestFailed {
        /// Horizons target command for the failed acquisition.
        command: String,
        /// Number of request attempts made before returning the failure.
        attempts: u32,
        /// Diagnostic retained from the failed request sequence.
        detail: String,
    },
    /// The supplied asset SHA-256 differs from the expected digest.
    HashMismatch {
        /// Expected SHA-256 digest string for the selected asset.
        expected: String,
        /// SHA-256 digest string computed from the supplied asset.
        actual: String,
    },
    /// An asset failed structural validation before installation.
    AssetRejected {
        /// Selected asset filename whose contents failed validation.
        filename: String,
        /// Structural validation requirement that the asset failed.
        reason: &'static str,
    },
    /// The file identification bytes do not identify a binary SPK.
    SpkUnknownFormat {
        /// Eight identification bytes read from the DAF/SPK header.
        magic: [u8; 8],
    },
    /// The DAF descriptor shape is not the SPK shape (ND=2, NI=6).
    SpkUnsupportedNdNi {
        /// Declared number of floating-point summary components; SPK requires two.
        nd: i32,
        /// Declared number of integer summary components; SPK requires six.
        ni: i32,
    },
    /// SPK metadata or addresses are invalid; `offset` is a byte offset.
    SpkCorrupt {
        /// Description of the invalid SPK metadata, address or record.
        what: String,
        /// Byte offset associated with the failure in the SPK image.
        offset: usize,
    },
    /// The SPK segment type has no implemented evaluator.
    SpkUnsupportedType {
        /// DAF descriptor's SPK data-type identifier lacking an evaluator.
        seg_type: i32,
    },
    /// The requested SPK file could not be opened.
    SpkFileNotFound {
        /// Requested SPK path retained for the diagnostic.
        path: String,
    },
    /// The requested NAIF target is absent from the SPK image.
    SpkBodyNotFound {
        /// Requested NAIF target identifier absent from the image.
        target: i32,
    },
    /// The image contains no planetary-moon target chainable to the solar-system barycenter.
    SpkNoValidMoons {
        /// SPK path whose targets provided no supported barycentric moon chain.
        path: String,
    },
    /// The target epoch is outside the stored coverage, in ephemeris seconds past J2000.
    SpkOutOfRange {
        /// NAIF target identifier of the segment being evaluated.
        target: i32,
        /// Requested TDB epoch, in seconds since J2000.0.
        et: f64,
        /// Segment coverage start, in TDB seconds since J2000.0.
        start: f64,
        /// Segment coverage end, in TDB seconds since J2000.0.
        end: f64,
    },
    /// Segments for one target declare more than one center of motion.
    SpkAmbiguousCenter {
        /// NAIF target whose segment descriptors disagree on the center identifier.
        target: i32,
    },
    /// The SPK frame is unsupported; frame 1 denotes J2000/ICRS.
    SpkUnsupportedFrame {
        /// Reference-frame identifier declared by the SPK segment.
        frame: i32,
    },
    /// The library body has no registered SPK target.
    SpkUnregisteredBody {
        /// Library body identifier with no registered SPK target mapping.
        body: i32,
    },
    /// A search interval has non-finite endpoints or is not increasing.
    InvalidSearchInterval {
        /// Supplied lower endpoint, in the calling search's coordinate units.
        a: f64,
        /// Supplied upper endpoint, in the same units as `a`.
        b: f64,
    },
    /// The search residual tolerance is not finite and positive.
    InvalidSearchTolerance {
        /// Rejected residual tolerance, in the searched function's units.
        value: f64,
    },
    /// A search step or derivative-magnitude floor is not finite and positive.
    InvalidSearchStep {
        /// Rejected step or derivative-magnitude floor, in the calling algorithm's units.
        value: f64,
    },
    /// The search iteration limit is zero.
    InvalidSearchIterations {
        /// Rejected maximum number of search iterations.
        max_iter: usize,
    },
    /// The search scan sample count is zero.
    InvalidScanSamples {
        /// Rejected number of bracket-scan samples.
        count: usize,
    },
    /// The search function returned a non-finite value at `x`.
    NonFiniteSearchValue {
        /// Search coordinate where the function evaluation was nonfinite.
        x: f64,
    },
    /// An SPK range failure retaining the segment metadata and its immutable image owner.
    BackendRange {
        /// Primary typed failure exposed to the caller.
        error: Box<Error>,
        /// Original SPK range metadata and its retained immutable image owner.
        raw: std::sync::Arc<crate::spk::OwnedSpkRange>,
    },
    /// A binary64 numerical operation overflowed; `operation` identifies the calculation.
    NumericalOverflow {
        /// Name of the calculation in which the binary64 overflow occurred.
        operation: &'static str,
        /// Diagnostic describing the failed numerical operation.
        message: String,
    },
    /// A native failure with its directly owned cause.
    Caused {
        /// Failure exposed to the caller.
        error: Box<Error>,
        /// Underlying failure exposed through [`std::error::Error::source`].
        cause: Box<Error>,
    },
    /// A crossing calculation failed, optionally retaining the native calculation cause.
    PublicCrossing {
        /// Diagnostic exposed for the failed public crossing calculation.
        message: String,
        /// Actual native calculation failure, when retained; `None` if no cause is attached.
        cause: Option<Box<Error>>,
    },
    /// No longitude crossing was bracketed in the requested Julian-day interval.
    NoCrossingInInterval {
        /// Lower Julian-day bound, in the invoking entry's time scale.
        lo: f64,
        /// Upper Julian-day bound, in the same time scale as `lo`.
        hi: f64,
        /// Requested longitude crossing, in degrees.
        target: f64,
    },
    /// A crossing search left its admissible Julian-day interval.
    CrossingDiverged {
        /// Initial search Julian day, in the invoking entry's time scale.
        start: f64,
        /// Divergent Julian-day iterate, in the same time scale as `start`.
        estimate: f64,
    },
    /// A longitude crossing search exhausted its iterations; residual is in degrees.
    CrossingNotConverged {
        /// Retained search estimate as a Julian day in the invoking entry's time scale.
        best: f64,
        /// Longitude residual at the retained estimate, in degrees.
        residual: f64,
    },
    /// No velocity sign change was found in the requested Julian-day interval.
    NoStationInInterval {
        /// Lower Julian-day search bound, in the invoking entry's time scale.
        lo: f64,
        /// Upper Julian-day search bound, in the same time scale as `lo`.
        hi: f64,
        /// Library body identifier whose longitude-speed sign change was sought.
        planet: i32,
    },
    /// A lunar-node search left its admissible Julian-day interval.
    MoonNodeDiverged {
        /// Initial lunar-node search Julian day, in the invoking entry's time scale.
        start: f64,
        /// Divergent Julian-day iterate, in the same time scale as `start`.
        estimate: f64,
    },
    /// A lunar-node search did not converge; latitude residual is in degrees.
    MoonNodeNotConverged {
        /// Retained node-search Julian day, in the invoking entry's time scale.
        best: f64,
        /// Lunar ecliptic latitude at the retained estimate, in degrees.
        latitude: f64,
    },
    /// The requested body has no retrograde stations.
    BodyHasNoStations {
        /// Library body identifier for which retrograde-station searches are unsupported.
        body: i32,
    },
    /// The requested library body identifier is unknown.
    UnknownBody {
        /// Unrecognized library body identifier.
        body_id: i32,
    },
    /// The known body is unsupported by the named operation.
    InvalidBody {
        /// Known library body identifier rejected by the operation.
        body_id: i32,
        /// Name of the operation that does not support the body.
        operation: &'static str,
    },
    /// The catalog query found no star; `search_type` identifies the query form when known.
    StarNotFound {
        /// Catalog query or identifier string that did not resolve to a star.
        star_id: String,
        /// Query form when known; `None` if the failing lookup did not attach one.
        search_type: Option<&'static str>,
    },
    /// A resolved star has no stored magnitude.
    StarMagnitudeUnavailable {
        /// Resolved catalog identifier whose record has no magnitude.
        star_id: i32,
    },
    /// Strict precision requires an SPK for the requested body.
    SpkRequired {
        /// Library body identifier for which strict precision requires SPK data.
        body: i32,
    },
    /// A registered SPK excludes the requested Julian day and strict precision forbids approximation.
    SpkRegisteredOutOfCoverage {
        /// Library body identifier whose registered SPK lacks the requested epoch.
        body: i32,
        /// Requested TT Julian day.
        requested_jd: f64,
    },
    /// The requested Julian day is outside the available ephemeris coverage.
    EphemerisRange {
        /// Requested Julian day in the public entry's time scale.
        requested_jd: f64,
        /// Available coverage start as a Julian day in the reported time scale.
        start_jd: f64,
        /// Available coverage end as a Julian day in the reported time scale.
        end_jd: f64,
        /// Associated library body identifier, or `None` when no body is attached.
        body: Option<i32>,
    },
    /// The selected house system is undefined at this latitude; angles are in degrees.
    PolarCircle {
        /// Observer latitude in degrees.
        latitude: f64,
        /// Polar latitude threshold applied by the selected house system, in degrees.
        threshold: f64,
        /// Obliquity used to determine the polar threshold, in degrees.
        obliquity: f64,
        /// House-system character whose latitude restriction was exceeded.
        house_system: char,
    },
    /// The house-system selector is not a single character.
    InvalidHouseSystem {
        /// Supplied house-system string that was not a single character.
        value: String,
    },
    /// A house cusp index is outside the inclusive range 1 through 12.
    InvalidHouseCuspIndex {
        /// Requested cusp index outside the supported one-based range.
        index: i64,
    },
    /// A sequential star index is not a decimal digit string.
    InvalidStarSearch {
        /// Supplied sequential-star selector that was not a decimal digit string.
        input: String,
    },
    /// A required configuration setting is absent.
    ConfigurationMissing {
        /// Name of the configuration required by the operation.
        missing_config: &'static str,
        /// Diagnostic guidance for supplying the missing configuration.
        suggestion: &'static str,
    },
    /// An iterative calculation exhausted its limit; tolerance and last value use algorithm units.
    ConvergenceFailed {
        /// Name of the iterative algorithm that failed to converge.
        algorithm: &'static str,
        /// Number of iterations completed before the failure.
        iterations: usize,
        /// Iteration limit applied by the algorithm.
        max_iterations: usize,
        /// Convergence threshold, in the named algorithm's units.
        tolerance: f64,
        /// Last iterate or diagnostic value, in the named algorithm's units.
        last_value: f64,
    },
    /// The named calculation does not support this flag combination.
    UnsupportedFlagCombination {
        /// Calculation entry that rejected the flags.
        entry: &'static str,
        /// Rejected calculation-flag bitmask as supplied to the check.
        flags: i32,
    },
    /// The selected body channel has no evaluator for the requested coordinate type.
    SealedChannelNotServed {
        /// Library body identifier stored in the unsupported channel.
        body_id: i32,
        /// Stored coordinate-type identifier for which this binding has no evaluator.
        coord_type: u32,
    },
    /// The selected sidereal mode requires an unavailable star anchor.
    AyanamshaAnchorUnserved {
        /// Sidereal mode identifier whose required live anchor is unavailable.
        mode: i32,
    },
    /// The selected dispatcher does not implement this house system.
    HouseSystemUnserved {
        /// Requested house-system character lacking a served cusp kernel.
        house_system: char,
    },
    /// The Gauquelin method requires an unavailable rise/set provider.
    GauquelinMethodUnserved {
        /// Requested Gauquelin-sector method lacking its required event provider.
        method: i32,
    },
    /// A required Gauquelin rise or set event could not be found.
    GauquelinRiseSetNotFound {
        /// Gauquelin-sector method for which a necessary rise or set was not found.
        method: i32,
    },
    /// The requested planetary-moon theory is unavailable.
    MoonTheoryUnserved {
        /// Requested planetary-moon theory name.
        system: &'static str,
    },
    /// The named table generator received an invalid grid, shape, order or step.
    InvalidGeneratorInput {
        /// Name of the generator that rejected its inputs.
        entry: &'static str,
        /// Description of the invalid grid, shape, order or step.
        detail: &'static str,
    },
    /// The eclipse type is neither solar nor lunar.
    InvalidEclipseType {
        /// Rejected eclipse-type selector.
        value: String,
    },
    /// An eclipse search input is invalid.
    InvalidEclipseSearch {
        /// Eclipse-search entry that rejected the request.
        entry: &'static str,
        /// Diagnostic explaining the invalid search input.
        detail: String,
    },
    /// An eclipse search reached its available coverage without finding an event.
    EclipseSearchExhausted {
        /// Eclipse-search entry whose available search coverage was exhausted.
        entry: &'static str,
        /// Diagnostic describing the exhausted search or unavailable event.
        detail: String,
    },
    /// The named orbital-elements calculation does not support this body.
    InvalidOrbitalObject {
        /// Orbital-elements entry that rejected the body selection.
        entry: &'static str,
        /// Diagnostic describing the unsupported orbital object.
        detail: String,
    },
    /// An orbital-elements input or file is invalid.
    InvalidOrbitElements {
        /// Orbital-elements entry that rejected an input or file.
        entry: &'static str,
        /// Diagnostic identifying the invalid element data or file.
        detail: String,
    },
    /// The input lies outside the domain of a numerical function.
    MathDomain {
        /// Value at which the numerical function's domain check failed.
        value: f64,
    },
    /// The N-body integrator name is unsupported.
    UnknownIntegrator {
        /// Unrecognized N-body integrator name.
        name: String,
    },
    /// The requested optional N-body backend or its data are unavailable.
    NBodyBackendUnavailable {
        /// Requested optional backend or backend-data configuration that is unavailable.
        backend: &'static str,
    },
    /// A direct REBOUND or ASSIST operation failed with the supplied stage and backend status.
    NBodyNative {
        /// Operation stage: 0 input, 1 data, 2 initialization, 3 orbit, 4 integration, 5 output.
        stage: i32,
        /// Backend status code; codes 101–106 identify invalid orbital elements.
        code: i32,
        /// Diagnostic supplied by the native backend.
        diagnostic: String,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyCoefficients => write!(f, "empty Chebyshev coefficient slice"),
            Self::NonFiniteArgument { value } => {
                write!(f, "non-finite evaluation argument: {value}")
            }
            Self::NormalizedTimeOutOfRange { x } => {
                write!(f, "normalized Chebyshev argument out of [-1, 1]: {x}")
            }
            Self::InvalidDerivativeOrder { order, max } => {
                write!(
                    f,
                    "derivative order {order} exceeds supported maximum {max}"
                )
            }
            Self::InvalidSegment { begin, end } => {
                write!(f, "invalid segment bounds: [{begin}, {end}]")
            }
            Self::EpochOutOfSegment { t, begin, end } => {
                write!(f, "epoch {t} outside owning segment [{begin}, {end}]")
            }
            Self::ZeroVectorNormalize => write!(f, "cannot normalize the zero vector"),
            Self::SingularMatrix => write!(f, "cannot invert a singular matrix"),
            Self::ZeroVectorToSpherical => {
                write!(f, "spherical angles of the zero vector are undefined")
            }
            Self::IntegerOverflow => write!(f, "exact integer arithmetic overflowed"),
            Self::RationalDivisionByZero => write!(f, "rational division by zero"),
            Self::IntervalInvertedBounds => {
                write!(f, "interval lower bound exceeds upper bound")
            }
            Self::IntervalSpansZero => {
                write!(f, "interval division by a divisor spanning zero")
            }
            Self::IntervalDomainViolation => {
                write!(f, "interval operation outside its domain")
            }
            Self::InvalidCalendar { value } => {
                write!(f, "invalid calendar flag: {value}")
            }
            Self::InvalidCalendarCode { code } => {
                write!(f, "calendar must be 'j' or 'g', got: '{code}'")
            }
            Self::InvalidDate { year, month, day } => {
                write!(
                    f,
                    "invalid date: year = {year}, month = {month}, day = {day}"
                )
            }
            Self::InvalidTime {
                hour,
                minute,
                second,
            } => {
                write!(f, "invalid time: {hour}:{minute}:{second:.2}")
            }
            Self::NoLeapSecond {
                hour,
                minute,
                second,
            } => {
                write!(
                    f,
                    "invalid time (no leap second!): {hour}:{minute}:{second:.2}"
                )
            }
            Self::InvalidObserver { field, value } => {
                write!(f, "invalid observer {field}: {value}")
            }
            Self::SegmentOwnershipUnresolved { t, begin, length } => {
                write!(
                    f,
                    "cannot establish ownership of epoch {t} on grid [{begin}, +{length}]: exact comparison overflowed"
                )
            }
            Self::DerivativeScaleOverflow { scale } => {
                write!(f, "non-finite segment derivative scale: {scale}")
            }
            Self::LebFileNotFound { path } => {
                write!(f, "LEB file not found: {path}")
            }
            Self::LebUnknownFormat { magic } => {
                write!(
                    f,
                    "unknown LEB format magic: 0x{:02x}{:02x}{:02x}{:02x}",
                    magic[0], magic[1], magic[2], magic[3]
                )
            }
            Self::LebUnsupportedVersion { found } => {
                write!(f, "unsupported LEB version: {found}")
            }
            Self::LebUnsupportedCompression { found } => {
                write!(f, "unsupported LEB2 compression scheme: {found}")
            }
            Self::LebCorrupt { what, offset } => {
                write!(f, "corrupt LEB {what} at byte {offset}")
            }
            Self::ResourceLimit {
                resource,
                required,
                limit,
            } => {
                write!(
                    f,
                    "{resource} requires {required}; configured limit is {limit}"
                )
            }
            Self::InvalidLebReader { value } => {
                write!(f, "invalid LEB reader configuration: {value}")
            }
            Self::LebBytesUnavailable => {
                write!(f, "indexed LEB storage has no borrowed byte image")
            }
            Self::LebBodyNotFound { body_id } => {
                write!(f, "body {body_id} not in LEB file")
            }
            Self::LebStarNotFound { star_id } => {
                write!(f, "star {star_id} not in LEB file")
            }
            Self::LebReaderClosed => {
                write!(f, "LEB reader is closed")
            }
            Self::LebOutOfRange {
                body_id,
                jd,
                start,
                end,
            } => {
                write!(
                    f,
                    "JD {jd} outside range [{start}, {end}] for body {body_id}"
                )
            }
            Self::DataManifestMissing { name } => {
                write!(f, "Reviewed manifest pin missing for {name}")
            }
            Self::DeltaTModelUnavailable { model } => {
                write!(f, "Delta-T model '{model}' has no usable backend")
            }
            Self::InvalidNutationModel { value } => {
                write!(f, "Unknown nutation model {value:?}; use nut00a or nut06a")
            }
            Self::InvalidDeltaTModel { value } => {
                write!(
                    f,
                    "Unknown delta_t model {value:?}; valid models: smh2016, espenak_meeus"
                )
            }
            Self::InvalidCalcMode { value, source } => {
                write!(
                    f,
                    "Invalid calculation mode from {source}: {value:?}. Must be one of: [\"auto\", \"skyfield\", \"leb\", \"horizons\"]"
                )
            }
            Self::InvalidTier { value } => {
                write!(
                    f,
                    "Invalid tier: {value:?}. Must be one of: [\"base\", \"medium\", \"extended\"]"
                )
            }
            Self::InvalidNetworkPolicy { value, source } => {
                write!(
                    f,
                    "Invalid network policy from {source}: {value:?}. Must be one of: [\"auto\", \"allow\", \"sealed\"]"
                )
            }
            Self::NetworkSealed { purpose } => {
                write!(
                    f,
                    "Network access is sealed; blocked operation: {purpose}. Provision required files explicitly with 'libephemeris download ...' before starting the service."
                )
            }
            Self::InvalidTimeout { seconds } => {
                write!(
                    f,
                    "invalid network timeout: {seconds} s (must be finite and positive)"
                )
            }
            Self::InvalidBodyLimit { limit_bytes } => {
                write!(
                    f,
                    "invalid response body limit: {limit_bytes} bytes (must be at least 1)"
                )
            }
            Self::InvalidPadding { entry, value } => {
                write!(
                    f,
                    "invalid window padding in {entry}: {value} days (must be finite and non-negative)"
                )
            }
            Self::HttpInvalidUrl { url } => {
                write!(f, "invalid request URL (need http(s)://host/...): {url}")
            }
            Self::HttpTimeout { purpose, seconds } => {
                write!(f, "network request timed out after {seconds} s: {purpose}")
            }
            Self::HttpStatus { status, purpose } => {
                write!(f, "HTTP status {status} for {purpose}")
            }
            Self::HttpTransport { purpose, detail } => {
                write!(f, "network transport failure for {purpose}: {detail}")
            }
            Self::EphemerisLoaderIo { path, message, .. } => {
                write!(f, "ephemeris loader I/O failure for {path}: {message}")
            }
            Self::HttpResponseTooLarge {
                purpose,
                limit_bytes,
            } => {
                write!(
                    f,
                    "response body for {purpose} exceeded the {limit_bytes}-byte in-memory bound"
                )
            }
            Self::HorizonsApiError { command, message } => {
                write!(f, "Horizons API error for {command}: {message}")
            }
            Self::HorizonsMalformed { command, reason } => {
                write!(f, "malformed Horizons response for {command}: {reason}")
            }
            Self::HorizonsRequestFailed {
                command,
                attempts,
                detail,
            } => {
                write!(
                    f,
                    "Horizons API request failed after {attempts} attempts for {command}: {detail}"
                )
            }
            Self::HashMismatch { expected, actual } => {
                write!(
                    f,
                    "downloaded asset failed SHA-256 verification: expected {expected}, got {actual}"
                )
            }
            Self::AssetRejected { filename, reason } => {
                write!(f, "downloaded asset {filename} failed validation: {reason}")
            }
            Self::SpkUnknownFormat { magic } => {
                write!(
                    f,
                    "unknown SPK format magic: 0x{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}{:02x}",
                    magic[0], magic[1], magic[2], magic[3], magic[4], magic[5], magic[6], magic[7]
                )
            }
            Self::SpkUnsupportedNdNi { nd, ni } => {
                write!(
                    f,
                    "unsupported DAF summary shape for SPK: ND = {nd}, NI = {ni} (need 2 and 6)"
                )
            }
            Self::SpkCorrupt { what, offset } => {
                write!(f, "corrupt SPK {what} at byte {offset}")
            }
            Self::SpkUnsupportedType { seg_type } => {
                write!(f, "unsupported SPK data type: {seg_type}")
            }
            Self::SpkFileNotFound { path } => {
                write!(f, "SPK file not found: {path}")
            }
            Self::SpkBodyNotFound { target } => {
                write!(f, "target {target} not in SPK kernel")
            }
            Self::SpkNoValidMoons { path } => {
                write!(f, "no valid planetary moons found in SPK kernel: {path}")
            }
            Self::SpkAmbiguousCenter { target } => {
                write!(
                    f,
                    "target {target} declares more than one SPK center; a single center per target is required"
                )
            }
            Self::SpkUnsupportedFrame { frame } => {
                write!(f, "unsupported SPK reference frame: {frame}")
            }
            Self::SpkUnregisteredBody { body } => {
                write!(f, "body {body} has no registered SPK kernel")
            }
            Self::InvalidSearchInterval { a, b } => {
                write!(f, "invalid search interval: [{a}, {b}]")
            }
            Self::InvalidSearchTolerance { value } => {
                write!(f, "invalid search tolerance: {value}")
            }
            Self::InvalidSearchStep { value } => {
                write!(f, "invalid search x-step floor: {value}")
            }
            Self::InvalidSearchIterations { max_iter } => {
                write!(f, "invalid search iteration cap: {max_iter}")
            }
            Self::InvalidScanSamples { count } => {
                write!(f, "invalid bracket-scan sample count: {count}")
            }
            Self::NonFiniteSearchValue { x } => {
                write!(f, "search function returned a non-finite value at x = {x}")
            }
            Self::BackendRange { error, .. } | Self::Caused { error, .. } => {
                fmt::Display::fmt(error.as_ref(), f)
            }
            Self::NumericalOverflow { message, .. } => f.write_str(message),
            Self::PublicCrossing { message, .. } => f.write_str(message),
            Self::NoCrossingInInterval { lo, hi, target } => {
                write!(f, "no crossing of {target} in search interval [{lo}, {hi}]")
            }
            Self::CrossingDiverged { start, estimate } => {
                write!(f, "crossing search from {start} diverged to {estimate}")
            }
            Self::CrossingNotConverged { best, residual } => {
                write!(
                    f,
                    "crossing search did not converge near {best} (residual {residual})"
                )
            }
            Self::NoStationInInterval { lo, hi, planet } => {
                write!(f, "no station found for planet {planet} in [{lo}, {hi}]")
            }
            Self::MoonNodeDiverged { start, estimate } => {
                write!(f, "lunar node search from {start} diverged to {estimate}")
            }
            Self::MoonNodeNotConverged { best, latitude } => {
                write!(
                    f,
                    "lunar node search did not converge near {best} (latitude {latitude})"
                )
            }
            Self::BodyHasNoStations { body } => {
                write!(f, "body {body} has no retrograde stations")
            }
            Self::SpkOutOfRange {
                target,
                et,
                start,
                end,
            } => {
                write!(
                    f,
                    "epoch {et} outside range [{start}, {end}] for target {target}"
                )
            }
            Self::UnknownBody { body_id } => {
                write!(f, "unknown body ID: {body_id}")
            }
            Self::InvalidBody { body_id, operation } => {
                write!(f, "body {body_id} is not valid for operation '{operation}'")
            }
            Self::StarNotFound {
                star_id,
                search_type,
            } => match search_type {
                Some(kind) => write!(
                    f,
                    "fixed star not found: '{star_id}' (search type '{kind}')"
                ),
                None => write!(f, "fixed star not found: '{star_id}'"),
            },
            Self::StarMagnitudeUnavailable { star_id } => {
                write!(f, "Magnitude not available for star ID {star_id}")
            }
            Self::SpkRequired { body } => {
                write!(
                    f,
                    "SPK kernel required for body {body} (strict precision mode)"
                )
            }
            Self::SpkRegisteredOutOfCoverage { body, requested_jd } => {
                write!(
                    f,
                    "SPK kernel registered for body {body} does not cover JD {requested_jd:.1} TT (strict precision mode)"
                )
            }
            Self::EphemerisRange {
                requested_jd,
                start_jd,
                end_jd,
                body,
            } => match body {
                Some(id) => write!(
                    f,
                    "JD {requested_jd} outside ephemeris range [{start_jd}, {end_jd}] for body {id}"
                ),
                None => write!(
                    f,
                    "JD {requested_jd} outside ephemeris range [{start_jd}, {end_jd}]"
                ),
            },
            Self::PolarCircle {
                latitude,
                threshold,
                obliquity,
                house_system,
            } => {
                write!(
                    f,
                    "house system '{house_system}' cannot be calculated at latitude {latitude} (polar threshold {threshold} at obliquity {obliquity})"
                )
            }
            Self::EphemerisLoaderValue { message } => f.write_str(message),
            Self::InvalidHouseCuspIndex { index } => {
                write!(f, "house cusp index {index} is outside 1..=12")
            }
            Self::InvalidHouseSystem { value } => {
                write!(
                    f,
                    "house-system identifier must be a single character, not {value:?}"
                )
            }
            Self::InvalidStarSearch { input } => {
                write!(
                    f,
                    "invalid sequential star number {input:?}: not a decimal digit string"
                )
            }
            Self::ConfigurationMissing {
                missing_config,
                suggestion,
            } => {
                write!(f, "missing configuration '{missing_config}': {suggestion}")
            }
            Self::ConvergenceFailed {
                algorithm,
                iterations,
                max_iterations,
                tolerance,
                last_value,
            } => {
                write!(
                    f,
                    "algorithm '{algorithm}' did not converge after {iterations}/{max_iterations} iterations (tolerance {tolerance}, last value {last_value})"
                )
            }
            Self::UnsupportedFlagCombination { entry, flags } => {
                write!(
                    f,
                    "{entry}: unsupported calculation flag combination (flags {flags})"
                )
            }
            Self::SealedChannelNotServed {
                body_id,
                coord_type,
            } => {
                write!(
                    f,
                    "body {body_id} is stored in LEB channel {coord_type}, which the sealed Pipeline-A binding does not serve"
                )
            }
            Self::AyanamshaAnchorUnserved { mode } => {
                write!(
                    f,
                    "sidereal mode {mode} needs live star/galactic anchor geometry, which no provider serves (STAR-04 scope)"
                )
            }
            Self::HouseSystemUnserved { house_system } => {
                write!(
                    f,
                    "house system '{house_system}' needs a cusp kernel this entry does not serve yet (HOU-04/05/06 scope)"
                )
            }
            Self::GauquelinMethodUnserved { method } => {
                write!(
                    f,
                    "gauquelin_sector method {method} needs rise/set event times, which this provider does not serve (EVT-14 scope)"
                )
            }
            Self::GauquelinRiseSetNotFound { method } => {
                write!(
                    f,
                    "gauquelin_sector method {method}: rise or set not found (circumpolar body, no method-0 fallback)"
                )
            }
            Self::MoonTheoryUnserved { system } => {
                write!(
                    f,
                    "moon theory '{system}' is not served by this build (BODY-05 scope)"
                )
            }
            Self::InvalidGeneratorInput { entry, detail } => {
                write!(f, "{entry}: invalid generator input ({detail})")
            }
            Self::InvalidEclipseType { value } => {
                write!(f, "eclipse_type must be 'solar' or 'lunar', got {value:?}")
            }
            Self::InvalidEclipseSearch { entry, detail } => {
                write!(f, "{entry}: {detail}")
            }
            Self::EclipseSearchExhausted { entry, detail } => {
                write!(f, "{entry}: {detail}")
            }
            Self::InvalidOrbitalObject { entry, detail } => {
                write!(f, "{entry}: {detail}")
            }
            Self::InvalidOrbitElements { entry, detail } => {
                write!(f, "{entry}: {detail}")
            }
            Self::MathDomain { value } => {
                write!(f, "math domain error for value: {value}")
            }
            Self::UnknownIntegrator { name } => {
                write!(f, "unknown N-body integrator: {name}")
            }
            Self::NBodyNative {
                stage,
                code,
                diagnostic,
                ..
            } => {
                write!(f, "REBOUND stage {stage}, code {code}: {diagnostic}")
            }
            Self::NBodyBackendUnavailable { backend } => {
                write!(
                    f,
                    "N-body backend '{backend}' is unavailable in this configuration"
                )
            }
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.cause().map(|cause| cause as &dyn std::error::Error)
    }
}

/// Broad failure class, independent of diagnostic wording.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ErrorCategory {
    /// Caller input is invalid or outside a numerical domain.
    InputValidation,
    /// Required data is absent or unusable.
    DataNotFound,
    /// A numerical calculation failed or exceeded its coverage.
    Calculation,
    /// Required configuration or session state is absent.
    Configuration,
}

/// A non-fatal numerical or configuration advisory, returned alongside a result.
#[derive(Debug, Clone, PartialEq)]
pub enum Warning {
    /// Lunar polynomial outside its optimal range; the epoch is in Julian centuries from J2000.
    MeeusPolynomial {
        /// Evaluation epoch in Julian centuries of 36,525 days from J2000.0.
        centuries_from_j2000: f64,
        /// Whether the epoch exceeds the validity range rather than only the optimal range.
        outside_valid_range: bool,
    },
    /// Unknown sidereal mode replaced by Fagan/Bradley (mode 0).
    UnknownSiderealMode {
        /// Unrecognized sidereal mode for which mode zero was selected.
        mode: i32,
    },
    /// Unsupported projection bits omitted from the sidereal mode.
    UnsupportedSiderealProjection {
        /// Unsupported projection-bit subset omitted from the retained sidereal selection.
        bits: i32,
    },
    /// A malformed configuration file was ignored.
    MalformedConfig {
        /// Configuration file path whose contents were ignored.
        path: String,
        /// Diagnostic explaining why the configuration could not be parsed.
        detail: String,
    },
    /// One ASSIST asteroid self-force was excluded because it coincides with the target.
    ReboundSelfInteraction {
        /// Name of the ASSIST asteroid whose force on the identical target was excluded.
        perturber: String,
    },
}

impl Warning {
    /// Stable machine-readable category for diagnostic filtering.
    pub fn category(&self) -> &'static str {
        match self {
            Self::MeeusPolynomial { .. } => "meeus_polynomial",
            Self::UnknownSiderealMode { .. } => "unknown_sidereal_mode",
            Self::UnsupportedSiderealProjection { .. } => "unsupported_sidereal_projection",
            Self::MalformedConfig { .. } => "malformed_config",
            Self::ReboundSelfInteraction { .. } => "rebound_self_interaction",
        }
    }
}

impl fmt::Display for Warning {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MeeusPolynomial {
                centuries_from_j2000,
                outside_valid_range,
            } => {
                if *outside_valid_range {
                    write!(
                        f,
                        "date is {centuries_from_j2000:.1} centuries from J2000.0, outside the valid range (+-20 centuries): error may exceed 1 degree"
                    )
                } else {
                    write!(
                        f,
                        "date is {centuries_from_j2000:.1} centuries from J2000.0, outside the optimal range (+-10 centuries): precision is degraded"
                    )
                }
            }
            Self::UnknownSiderealMode { mode } => {
                write!(
                    f,
                    "unknown sidereal mode {mode} is not recognized; falling back to Fagan/Bradley (mode 0)"
                )
            }
            Self::UnsupportedSiderealProjection { bits } => {
                write!(
                    f,
                    "SIDBIT projection flags 0x{bits:04X} are not supported in this version and were ignored"
                )
            }
            Self::MalformedConfig { path, detail } => {
                write!(f, "ignoring malformed libephemeris config {path}: {detail}")
            }
            Self::ReboundSelfInteraction { perturber } => {
                write!(
                    f,
                    "propagate: the body being propagated coincides with the ASSIST asteroid perturber '{perturber}', so that perturber's self-force was excluded from this integration; all other configured asteroid forces remain active."
                )
            }
        }
    }
}

impl Error {
    /// Attach an actual SPK range failure and retain its immutable image owner.
    /// Replacing this metadata does not change the primary failure or its cause.
    pub fn with_backend_range(self, raw: crate::spk::OwnedSpkRange) -> Self {
        let raw = std::sync::Arc::new(raw);
        match self {
            Self::BackendRange { error, .. } => Self::BackendRange { error, raw },
            error => Self::BackendRange {
                error: Box::new(error),
                raw,
            },
        }
    }

    /// Borrow SPK segment metadata without traversing a different failure's cause.
    pub fn backend_range(&self) -> Option<&crate::spk::OwnedSpkRange> {
        match self {
            Self::BackendRange { raw, .. } => Some(raw.as_ref()),
            Self::Caused { error, .. } => error.backend_range(),
            _ => None,
        }
    }

    /// Attach the native failure that caused this error.
    /// The primary category and diagnostic remain unchanged.
    pub fn with_cause(self, cause: Self) -> Self {
        Self::Caused {
            error: Box::new(self),
            cause: Box::new(cause),
        }
    }

    /// Borrow the primary typed failure beneath native metadata and cause wrappers.
    pub fn primary_error(&self) -> &Self {
        match self {
            Self::BackendRange { error, .. } | Self::Caused { error, .. } => error.primary_error(),
            _ => self,
        }
    }

    /// Borrow the directly owned native cause, also exposed by standard `source`.
    pub fn cause(&self) -> Option<&Self> {
        match self {
            Self::BackendRange { error, .. } => error.cause(),
            Self::Caused { cause, .. } => Some(cause.as_ref()),
            Self::PublicCrossing { cause, .. } => cause.as_deref(),
            _ => None,
        }
    }

    /// Broad native failure category.
    pub fn category(&self) -> ErrorCategory {
        match self {
            Self::BackendRange { error, .. } | Self::Caused { error, .. } => error.category(),
            Self::NumericalOverflow { .. }
            | Self::EmptyCoefficients
            | Self::NonFiniteArgument { .. }
            | Self::NormalizedTimeOutOfRange { .. }
            | Self::InvalidDerivativeOrder { .. }
            | Self::InvalidSegment { .. }
            | Self::ZeroVectorNormalize
            | Self::SingularMatrix
            | Self::ZeroVectorToSpherical
            | Self::IntegerOverflow
            | Self::RationalDivisionByZero
            | Self::IntervalInvertedBounds
            | Self::IntervalSpansZero
            | Self::IntervalDomainViolation
            | Self::InvalidCalendar { .. }
            | Self::InvalidCalendarCode { .. }
            | Self::InvalidDate { .. }
            | Self::InvalidTime { .. }
            | Self::NoLeapSecond { .. }
            | Self::InvalidObserver { .. }
            | Self::SegmentOwnershipUnresolved { .. }
            | Self::DerivativeScaleOverflow { .. }
            | Self::InvalidDeltaTModel { .. }
            | Self::InvalidNutationModel { .. }
            | Self::InvalidCalcMode { .. }
            | Self::InvalidTier { .. }
            | Self::InvalidNetworkPolicy { .. }
            | Self::InvalidTimeout { .. }
            | Self::InvalidBodyLimit { .. }
            | Self::InvalidPadding { .. }
            | Self::HttpInvalidUrl { .. }
            | Self::InvalidSearchInterval { .. }
            | Self::InvalidSearchTolerance { .. }
            | Self::InvalidSearchStep { .. }
            | Self::InvalidSearchIterations { .. }
            | Self::InvalidScanSamples { .. }
            | Self::NonFiniteSearchValue { .. }
            | Self::InvalidBody { .. }
            | Self::BodyHasNoStations { .. }
            | Self::UnsupportedFlagCombination { .. }
            | Self::InvalidGeneratorInput { .. }
            | Self::InvalidEclipseType { .. }
            | Self::InvalidEclipseSearch { .. }
            | Self::InvalidOrbitalObject { .. }
            | Self::InvalidOrbitElements { .. }
            | Self::MathDomain { .. }
            | Self::EphemerisLoaderValue { .. }
            | Self::InvalidHouseSystem { .. }
            | Self::InvalidHouseCuspIndex { .. }
            | Self::InvalidStarSearch { .. }
            | Self::UnknownIntegrator { .. } => ErrorCategory::InputValidation,
            Self::NBodyNative { stage: 0, .. } => ErrorCategory::InputValidation,
            Self::NBodyNative { stage: 1, .. } => ErrorCategory::DataNotFound,
            Self::NBodyNative { stage: 3, code, .. } if (101..=106).contains(code) => {
                ErrorCategory::InputValidation
            }
            Self::NBodyNative { .. } => ErrorCategory::Calculation,
            Self::LebFileNotFound { .. }
            | Self::LebUnknownFormat { .. }
            | Self::LebUnsupportedVersion { .. }
            | Self::LebUnsupportedCompression { .. }
            | Self::LebCorrupt { .. }
            | Self::LebBodyNotFound { .. }
            | Self::LebStarNotFound { .. }
            | Self::DataManifestMissing { .. }
            | Self::DeltaTModelUnavailable { .. }
            | Self::SpkUnknownFormat { .. }
            | Self::SpkUnsupportedNdNi { .. }
            | Self::SpkCorrupt { .. }
            | Self::SpkUnsupportedType { .. }
            | Self::SpkFileNotFound { .. }
            | Self::EphemerisLoaderIo { .. }
            | Self::SpkBodyNotFound { .. }
            | Self::SpkNoValidMoons { .. }
            | Self::SpkAmbiguousCenter { .. }
            | Self::SpkUnsupportedFrame { .. }
            | Self::SpkUnregisteredBody { .. }
            | Self::UnknownBody { .. }
            | Self::StarNotFound { .. }
            | Self::StarMagnitudeUnavailable { .. }
            | Self::SpkRequired { .. }
            | Self::SpkRegisteredOutOfCoverage { .. }
            | Self::SealedChannelNotServed { .. }
            | Self::AyanamshaAnchorUnserved { .. }
            | Self::HouseSystemUnserved { .. }
            | Self::GauquelinMethodUnserved { .. }
            | Self::GauquelinRiseSetNotFound { .. }
            | Self::MoonTheoryUnserved { .. }
            | Self::NBodyBackendUnavailable { .. }
            | Self::HashMismatch { .. }
            | Self::AssetRejected { .. }
            | Self::HorizonsApiError { .. } => ErrorCategory::DataNotFound,
            Self::EpochOutOfSegment { .. }
            | Self::LebOutOfRange { .. }
            | Self::SpkOutOfRange { .. }
            | Self::PublicCrossing { .. }
            | Self::NoCrossingInInterval { .. }
            | Self::CrossingDiverged { .. }
            | Self::CrossingNotConverged { .. }
            | Self::NoStationInInterval { .. }
            | Self::MoonNodeDiverged { .. }
            | Self::MoonNodeNotConverged { .. }
            | Self::EphemerisRange { .. }
            | Self::PolarCircle { .. }
            | Self::EclipseSearchExhausted { .. }
            | Self::ConvergenceFailed { .. } => ErrorCategory::Calculation,
            Self::LebReaderClosed
            | Self::LebBytesUnavailable
            | Self::InvalidLebReader { .. }
            | Self::ResourceLimit { .. }
            | Self::NetworkSealed { .. }
            | Self::HttpTimeout { .. }
            | Self::HttpStatus { .. }
            | Self::HttpTransport { .. }
            | Self::HttpResponseTooLarge { .. }
            | Self::HorizonsMalformed { .. }
            | Self::HorizonsRequestFailed { .. }
            | Self::ConfigurationMissing { .. } => ErrorCategory::Configuration,
        }
    }

    /// Whether required data or network provisioning is absent or unusable.
    /// Range misses, caller input and numerical search outcomes are not provisioning failures.
    pub fn is_infrastructure(&self) -> bool {
        if let Self::BackendRange { error, .. } | Self::Caused { error, .. } = self {
            return error.is_infrastructure();
        }
        matches!(
            self,
            Self::LebFileNotFound { .. }
                | Self::LebUnknownFormat { .. }
                | Self::LebUnsupportedVersion { .. }
                | Self::LebUnsupportedCompression { .. }
                | Self::LebCorrupt { .. }
                | Self::DataManifestMissing { .. }
                | Self::DeltaTModelUnavailable { .. }
                | Self::NetworkSealed { .. }
                | Self::HttpTimeout { .. }
                | Self::HttpStatus { .. }
                | Self::HttpTransport { .. }
                | Self::EphemerisLoaderIo { .. }
                | Self::HttpResponseTooLarge { .. }
                | Self::HorizonsMalformed { .. }
                | Self::HorizonsRequestFailed { .. }
                | Self::HashMismatch { .. }
                | Self::AssetRejected { .. }
                | Self::SpkUnknownFormat { .. }
                | Self::SpkUnsupportedNdNi { .. }
                | Self::SpkCorrupt { .. }
                | Self::SpkUnsupportedType { .. }
                | Self::SpkFileNotFound { .. }
                | Self::SpkRequired { .. }
                | Self::SpkRegisteredOutOfCoverage { .. }
                | Self::NBodyNative { stage: 1, .. }
        )
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error as _;

    #[test]
    fn resource_refusals_are_configuration_not_corruption_or_provisioning() {
        for error in [
            Error::ResourceLimit {
                resource: "LEB decoded frame",
                required: 1024,
                limit: 512,
            },
            Error::InvalidLebReader {
                value: "typo".into(),
            },
            Error::LebBytesUnavailable,
        ] {
            assert_eq!(error.category(), ErrorCategory::Configuration);
            assert!(!error.is_infrastructure());
            assert!(!error.to_string().is_empty());
        }
    }

    #[test]
    fn numerical_cause_retains_primary_and_standard_source() {
        let primary = Error::ConvergenceFailed {
            algorithm: "contact refinement",
            iterations: 80,
            max_iterations: 80,
            tolerance: 1e-9,
            last_value: 0.002,
        };
        let cause = Error::NonFiniteSearchValue { x: 2_451_545.0 };
        let diagnostic = primary.to_string();
        let error = primary.clone().with_cause(cause.clone());
        assert_eq!(error.primary_error(), &primary);
        assert_eq!(error.category(), ErrorCategory::Calculation);
        assert_eq!(error.to_string(), diagnostic);
        assert_eq!(error.cause(), Some(&cause));
        assert_eq!(
            error.source().unwrap().downcast_ref::<Error>(),
            Some(&cause)
        );
        assert!(cause.source().is_none());
    }

    #[test]
    fn range_miss_is_not_corrupt_data_or_missing_provisioning() {
        let miss = Error::SpkOutOfRange {
            target: 499,
            et: 1.0,
            start: 2.0,
            end: 3.0,
        };
        assert_eq!(miss.category(), ErrorCategory::Calculation);
        assert!(!miss.is_infrastructure());
        let corrupt = Error::SpkCorrupt {
            what: "segment records".into(),
            offset: 1024,
        };
        assert_eq!(corrupt.category(), ErrorCategory::DataNotFound);
        assert!(corrupt.is_infrastructure());
        assert!(miss.clone().with_cause(corrupt.clone()).source().is_some());
        assert!(!miss.with_cause(corrupt).is_infrastructure());
    }

    #[test]
    fn crossing_reports_the_actual_calculation_cause() {
        let cause = Error::SpkBodyNotFound { target: 499 };
        let error = Error::PublicCrossing {
            message: "crossing position failed".into(),
            cause: Some(Box::new(cause.clone())),
        };
        assert_eq!(
            error.source().unwrap().downcast_ref::<Error>(),
            Some(&cause)
        );
        assert_eq!(error.category(), ErrorCategory::Calculation);
    }

    #[test]
    fn backend_stage_distinguishes_data_input_and_calculation_failures() {
        let error = |stage, code| Error::NBodyNative {
            stage,
            code,
            diagnostic: "native failure".into(),
        };
        assert_eq!(error(0, 1).category(), ErrorCategory::InputValidation);
        assert_eq!(error(1, 1).category(), ErrorCategory::DataNotFound);
        assert!(error(1, 1).is_infrastructure());
        assert_eq!(error(3, 102).category(), ErrorCategory::InputValidation);
        assert_eq!(error(4, 3).category(), ErrorCategory::Calculation);
    }

    #[test]
    fn meeus_bands_retain_distinct_numerical_advisories() {
        let degraded = Warning::MeeusPolynomial {
            centuries_from_j2000: 12.0,
            outside_valid_range: false,
        };
        let outside = Warning::MeeusPolynomial {
            centuries_from_j2000: 25.0,
            outside_valid_range: true,
        };
        assert_eq!(degraded.category(), outside.category());
        assert!(degraded.to_string().contains("precision is degraded"));
        assert!(outside.to_string().contains("error may exceed 1 degree"));
    }
}
