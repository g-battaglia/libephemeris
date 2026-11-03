//! Native DAF/SPK segment lookup, evaluation, and kernel registry.
//!
//! Quantity decoded: the DAF file record, summary-record chain, and SPK
//! segment descriptors of consumed JPL kernels — identity, descriptor
//! shape, summary-chain head, binary format, per-target segment tables
//! with closed-interval coverage, and the type 2/3 Chebyshev plus type 21
//! Extended Modified Difference Array payloads. Raw decoders/evaluators consume
//! caller-held bytes; [`SpkRegistry`] separately owns file acquisition, immutable
//! image snapshots, registrations and lazily retained decoded segments.
//!
//! DAF/SPK layout follows the NAIF file/summary/name records, closed
//! descriptor priority and type2/3 Chebyshev and type21 difference-line
//! structures. The type21 recurrence is Krogh's Modified Difference Array
//! method. Hand-built little/big-endian fixtures check domains, values,
//! analytic velocities, malformed metadata and immutable resource ownership.
//!
//! # Ordered layers
//!
//! * Container decode: [`FileRecord`] identity and [`DafEndian`] binary
//!   format, the summary-record chain walk, and [`SpkIndex`] segment
//!   tables with [`SpkIndex::select`] latest-wins resolution.
//! * Type 2/3 evaluation: [`decode_chebyshev`] directory/record decode
//!   and [`eval_chebyshev`] value/velocity evaluation with analytic
//!   type-2 velocity.
//! * Type 21 evaluation: [`decode_emda`] segment decode,
//!   [`emda_record_index`] strict-`>` record selection, and [`eval_emda`]
//!   Krogh-recurrence evaluation.
//! * Chain composition: [`resolve_chain`] target/center chains over
//!   [`CenterKind`], [`eval_relative`] center-relative states, and
//!   [`eval_ssb`] composed barycentric states behind the J2000-only
//!   [`check_frame`] gate with the [`icrs_to_ecliptic_j2000`] rotation.
//! * Registration and lifetime: [`SpkRegistry`] kernel images with
//!   target/center resolution, plus local discovery
//!   ([`discover_local_spks`]) and path resolution
//!   ([`resolve_spk_path`]).
//! * Discovery and automatic acquisition: NAIF-ID deduction, the Horizons SPK
//!   request/download composition, day-granular [`jd_to_iso_date`] /
//!   [`iso_to_jd`] converters, tier date windows, and [`AutoSpkConfig`]
//!   descriptors. Path/name planning is deterministic. Explicit download
//!   operations enforce the supplied network policy and structurally validate
//!   transported SPK bytes before publication. Registry file reads retain
//!   acquired images; process-environment acquisition stays context-side.
//!
//! # Frames, units, epochs
//!
//! Words are 8-byte double slots; addresses are 1-based word numbers.
//! Segment epochs are ephemeris seconds past J2000 (JD 2451545.0 TDB);
//! descriptor coverage is closed: `start <= et <= end`. Actual record
//! directories and polynomial domains can further refuse an advertised epoch.
//! Priority: for one target, the segment stored later in the file wins
//! over earlier overlapping segments ("latest supersedes"). Epoch
//! converters ([`jd_tdb_to_et`] / [`et_to_jd_tdb`]) bridge JD TDB and
//! ephemeris seconds. Evaluated states are km/km/s with
//! [`ChebState::position_au`] / [`ChebState::velocity_au_per_day`] AU
//! converters; the served frame is J2000 only. Calendar helpers are
//! day-granular proleptic-Gregorian (Meeus order, floor semantics).
//!
//! # Refusals
//!
//! Non-finite epochs are [`Error::NonFiniteArgument`]; truncated, over-long, or inconsistent
//! metadata is [`Error::SpkCorrupt`], never a silent mis-selection.
//! Segments naming any frame but J2000 and types outside 2/3 + 21 are
//! refused; missing targets keep their typed errors, and out-of-range
//! epochs are never silently extrapolated.
use crate::constants::AUNIT_TO_KM;
use crate::errors::{Error, Result};
use std::cell::Cell;
use std::collections::HashMap;
use std::path::Path;
use std::sync::{Arc, OnceLock};

// ---------------------------------------------------------------------------
// Layout constants
// ---------------------------------------------------------------------------

/// DAF physical record size in bytes (128 double words).
pub const RECORD_BYTES: usize = 1024;
/// Words per DAF record.
pub const WORDS_PER_RECORD: usize = 128;
/// File-record identification word of a binary SPK kernel.
pub const SPK_FILE_ID: &[u8; 8] = b"DAF/SPK ";
/// `ND` (double components per summary) required for SPK descriptors.
pub const SPK_ND: i32 = 2;
/// `NI` (integer components per summary) required for SPK descriptors.
pub const SPK_NI: i32 = 6;
/// Summary size in double words: `ND + (NI + 1) / 2`.
pub const SPK_SUMMARY_WORDS: usize = 5;
/// Segment-name size in bytes: `8 * SUMMARY_WORDS`.
pub const SPK_NAME_BYTES: usize = 40;
/// Maximum summaries per summary record: `125 / SUMMARY_WORDS`.
pub const MAX_SUMMARIES_PER_RECORD: usize = 25;
/// J2000 epoch: segment epochs count ephemeris seconds from this JD (TDB).
pub const J2000_JD_TDB: f64 = 2451545.0;
/// Seconds per Julian day, for epoch conversion.
pub const SECONDS_PER_DAY: f64 = 86400.0;

// File-record byte offsets (NAIF DAF Required Reading).
const OFF_ID: usize = 0;
const OFF_ND: usize = 8;
const OFF_NI: usize = 12;
const OFF_FWARD: usize = 76;
const OFF_BWARD: usize = 80;
const OFF_FORMAT: usize = 88;

// ---------------------------------------------------------------------------
// Binary format (endianness)
// ---------------------------------------------------------------------------

/// Numeric binary format declared by the DAF file record (`LOCFMT`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DafEndian {
    /// `"LTL-IEEE"`: least-significant byte first.
    Little,
    /// `"BIG-IEEE"`: most-significant byte first.
    Big,
}

// ---------------------------------------------------------------------------
// File record
// ---------------------------------------------------------------------------

/// Decoded DAF file record (record 1): identity, descriptor shape,
/// summary-chain head, and binary format.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileRecord {
    /// Summary-chain head: record number of the first summary record
    /// (`FWARD`; 0 when the file holds no segments).
    pub first_summary_record: u64,
    /// Record number of the final summary record (`BWARD`).
    pub last_summary_record: u64,
    /// Numeric binary format of the file (`LOCFMT`).
    pub endian: DafEndian,
}

// ---------------------------------------------------------------------------
// Segment descriptors
// ---------------------------------------------------------------------------

/// Unpacked SPK segment descriptor: coverage in ephemeris seconds past
/// J2000 plus the six integer components in NAIF order
/// (target, center, frame, type, initial address, final address).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Descriptor {
    /// Initial epoch of coverage (seconds past J2000, TDB).
    pub start_et: f64,
    /// Final epoch of coverage (seconds past J2000, TDB).
    pub end_et: f64,
    /// NAIF integer code of the target body.
    pub target: i32,
    /// NAIF integer code of the center of motion.
    pub center: i32,
    /// NAIF integer code of the reference frame.
    pub frame: i32,
    /// SPK data-type code (2/3 Chebyshev and 21 EMDA in Phase 1 scope).
    pub seg_type: i32,
    /// 1-based word address of the first array element.
    pub initial_address: i32,
    /// 1-based word address of the final array element.
    pub final_address: i32,
}

/// One SPK segment: descriptor plus its 40-byte Fortran name.
#[derive(Debug, Clone, PartialEq)]
pub struct Segment {
    /// Unpacked descriptor (coverage, codes, array addresses).
    pub descriptor: Descriptor,
    /// Raw segment name (space-padded Fortran string).
    pub name: [u8; SPK_NAME_BYTES],
}

// ---------------------------------------------------------------------------
// Kernel index
// ---------------------------------------------------------------------------

/// All segments of one SPK kernel image, in file order.
///
/// File order is priority order: [`SpkIndex::select`] returns the last
/// covering segment for a target (NAIF "latest supersedes").
#[derive(Debug, Clone, PartialEq)]
pub struct SpkIndex {
    segments: Vec<Segment>,
}

// ---------------------------------------------------------------------------
// Type 2/3 Chebyshev evaluation
// ---------------------------------------------------------------------------

/// SPK data type 2: Chebyshev polynomials, position only (velocity by
/// differentiation). Used for planet barycenters and integrated
/// satellites (DE440/DE441, planet-center subsets).
pub const TYPE_CHEBYSHEV_POS: i32 = 2;
/// SPK data type 3: Chebyshev polynomials, position and velocity sets.
/// Used for analytically computed satellites.
pub const TYPE_CHEBYSHEV_POS_VEL: i32 = 3;

/// Decoded type 2/3 segment: directory plus owned record words.
///
/// Layout (NAIF SPK Required Reading): `N` records of `RSIZE` words
/// followed by the 4-word directory `[INIT, INTLEN, RSIZE, N]`.
/// Each record holds `MID`, `RADIUS`, then 3 (type 2) or 6 (type 3)
/// coefficient sets of `(RSIZE - 2) / 3` or `/ 6` words. Epochs are
/// ephemeris seconds past J2000 (TDB); positions are km, velocities
/// km/s.
#[derive(Debug, Clone, PartialEq)]
pub struct ChebSegment {
    seg_type: i32,
    init_et: f64,
    int_len: f64,
    record_words: usize,
    n_records: usize,
    words: Vec<f64>,
    data_offset: usize,
}

/// Geometric state from a type 2/3 record, NAIF-native units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChebState {
    /// Position relative to the segment center, km (J2000/ICRS frame
    /// carried by the descriptor; frame conversion belongs to the
    /// chain-composition layer).
    pub position_km: [f64; 3],
    /// Velocity relative to the segment center, km/s.
    pub velocity_km_s: [f64; 3],
}

// ---------------------------------------------------------------------------
// Type 21 Extended Modified Difference Arrays
// ---------------------------------------------------------------------------

/// SPK data type 21: Extended Modified Difference Arrays, the format JPL
/// Horizons returns for most asteroid/comet kernels.
pub const TYPE_EMDA: i32 = 21;
/// NAIF `spk21.inc` bound on the difference-table dimension (`MAXTRM`),
/// also recorded in the release inventory.
pub const EMDA_MAXTRM: usize = 25;

/// Decoded type 21 segment: owned records plus the final-epoch table.
///
/// Layout (NAIF SPK Required Reading): `N` records of `DLSIZE =
/// `4*MAXDIM + 11` words, the `N` record final epochs, `N/100`
/// directory epochs (every 100th), then `MAXDIM`, then `N`. `MAXDIM`
/// (not the documented `DLSIZE`) sits at the penultimate word of
/// actual data — the layout the release reader also implements.
/// Epochs are ephemeris seconds past J2000 (TDB); states are km/km/s.
#[derive(Debug, Clone, PartialEq)]
pub struct EmdaSegment {
    maxdim: usize,
    n_records: usize,
    dlsize: usize,
    records: Vec<f64>,
    epochs: Vec<f64>,
    data_offset: usize,
}

// ---------------------------------------------------------------------------
// Target/center chain composition and frame handling
// ---------------------------------------------------------------------------

/// NAIF frame code served by this module: J2000, the native ICRS-aligned
/// frame of DE and Horizons kernels in the release inventory
/// (`docs/data.md`). Segments naming any other frame are
/// refused by [`check_frame`]; inventing a rotation for them would
/// corrupt every downstream pipeline (the apparent-pipeline layers own
/// frame transforms).
pub const FRAME_J2000: i32 = 1;
/// NAIF center code of the solar-system barycenter. A segment declaring
/// this center is already SSB-relative: its state is served unchanged
/// and no base-ephemeris lookup runs.
pub const CENTER_SSB: i32 = 0;
/// NAIF center code of the Sun. Resolved in the base ephemeris like any
/// other non-barycentric center. Sun and SSB links are exempt from the
/// registration reach-table check; evaluation still needs the Sun state.
pub const CENTER_SUN: i32 = 10;
/// Fixed J2000 obliquity in degrees used by the ICRS-to-ecliptic-J2000
/// rotation. This scalar frame convention is distinct from a date-dependent
/// true or mean obliquity model.
pub const OBLIQUITY_J2000_DEG: f64 = 23.4392911;

/// Classification of a segment-declared center: the barycenter needs no
/// lookup; the Sun and every other center resolve in the base ephemeris
/// through a [`CenterProvider`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CenterKind {
    /// Solar-system barycenter (code 0): state identically zero.
    Barycenter,
    /// Sun (code 10): base-ephemeris lookup.
    Sun,
    /// Any other NAIF center id: base-ephemeris lookup.
    Other(i32),
}

/// Registration-time chain for one target: the target id, the single
/// center its segments declare, and the verified reference frame.
///
/// Built by [`resolve_chain`]; the registry layer stores it at
/// registration so an
/// unknown target, an ambiguous center, or an unsupported frame fails
/// at registration, never per-evaluation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TargetChain {
    /// NAIF integer code of the target body.
    pub target: i32,
    /// NAIF integer code of the single declared center of motion.
    pub center: i32,
    /// Verified reference frame (always [`FRAME_J2000`] for now).
    pub frame: i32,
}

/// Base-ephemeris source for non-barycentric segment centers.
///
/// Returns the center's SSB-relative ICRS/J2000 state in km and km/s at
/// ephemeris seconds past J2000. The epoch reaching the provider already
/// passed the caller's coverage selection and is finite; providers must
/// still answer honestly (typed error, never invented data) when their
/// own coverage ends. Target descriptor coverage does not establish center
/// coverage. Implementations own acquisition/caching and may perform I/O;
/// synchronous calls do not imply purity, and provider panics are not caught.
pub trait CenterProvider {
    /// SSB-relative ICRS/J2000 state of NAIF `center` at TDB seconds past J2000.
    /// Return km and km/s, or a typed error for missing/uncovered data.
    fn center_state(&self, center: i32, et: f64) -> Result<ChebState>;

    /// Center state at a two-part TDB Julian date, in km and km/s.
    ///
    /// The default preserves compatibility with scalar-ET implementations by
    /// collapsing ephemeris seconds. Fraction-preserving providers override
    /// this method; both date parts must be finite. Registration/coverage and
    /// frame faults remain the provider's responsibility.
    fn center_state_two_part_jd(&self, center: i32, jd: f64, fraction: f64) -> Result<ChebState> {
        let (_, _, hi, lo) = evaluation::jd_parts_seconds(jd, fraction)?;
        self.center_state(center, hi + lo)
    }
}

// Kernel registration, selection, and lifetime

/// One loaded kernel image: owned bytes plus the decoded summary index.
///
/// Owning the bytes here is the lifetime contract: segment evaluation
/// reads words out of `bytes` through descriptor addresses, so the image
/// must outlive every evaluation. The registry (or a standalone owner)
/// holds the `SpkKernel`; evaluations borrow it (`&self`), and
/// [`SpkRegistry::remove_kernel`] / [`SpkRegistry::clear`] drop images
/// whose dependents are released with them (body registrations cascade).
/// Clones share immutable image bytes and the decoded index. Registry mutation
/// replaces or drops its own handles; active snapshots and captured refusals
/// retain their image without copying the full kernel.
///
/// `PartialEq` compares the full stored image by value. Contexts own independent
/// registration maps while their immutable image allocations may be shared.
#[derive(Clone)]
pub struct SpkKernel {
    path: String,
    bytes: Arc<[u8]>,
    endian: DafEndian,
    index: Arc<SpkIndex>,
    source_file: Option<Arc<std::fs::File>>,
    file_image_owner: Option<Arc<file_image_cache::FileImageOwner>>,
}

/// Equality retains the image-value contract independently of whether its
/// original read-only acquisition descriptor remains owned for presentation.
impl PartialEq for SpkKernel {
    fn eq(&self, other: &Self) -> bool {
        self.path == other.path
            && self.bytes == other.bytes
            && self.endian == other.endian
            && self.index == other.index
    }
}

/// Registration of one library body id onto a kernel NAIF target.
///
/// Mirrors one release `_SPK_BODY_MAP` entry `(spk_file, naif_id)`: `body`
/// is the caller's library body key (e.g. Chiron), `naif_id` the NAIF
/// target served from the kernel stored under `kernel`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BodyRegistration {
    /// Library body key.
    pub body: i32,
    /// Registry key of the serving kernel image.
    pub kernel: String,
    /// NAIF target id served from that kernel.
    pub naif_id: i32,
}

/// Owned SPK kernel store plus the body-to-kernel map.
///
/// Selection has two levels: the body map picks
/// exactly one `(kernel, NAIF target)` per library body (no cross-kernel
/// fallback, no Keplerian fall-through — an unregistered body is
/// [`Error::SpkUnregisteredBody`]), and inside the kernel the latest-wins
/// covering segment serves the epoch ([`SpkIndex::select`]).
///
/// Re-registering a body replaces its mapping; replacing a kernel installs a
/// fresh image/cache without revalidating dependent body mappings. Later
/// evaluation reports any missing target, ambiguous center, unsupported frame
/// or unsupported type through its corresponding typed error.
///
/// File registration retains a whole immutable byte image. Registry evaluation
/// lazily decodes selected segment payloads, retaining both successful decoded
/// segments and decode failures for that image. This storage has no LEB-style
/// configurable coefficient budget or per-request eviction: it grows with the
/// distinct segments decoded, not just the number of loaded kernel bytes.
///
/// Cloning shares image/index/decoded storage while copying registration maps.
/// Removing/replacing kernels or clearing releases this registry's handles;
/// snapshots and other clones can keep prior allocations and decode errors
/// alive until their last owner is dropped. Unregistering only a body does not
/// evict its kernel or decoded segments.
#[derive(Debug, Clone, Default)]
pub struct SpkRegistry {
    kernels: HashMap<String, SpkKernel>,
    bodies: HashMap<i32, (String, i32)>,
    decoded: HashMap<String, Arc<DecodedSpkCache>>,
    /// Automatically selected center target → (owned image key, system barycenter).
    automatic_centers: HashMap<i32, (String, i32)>,
}

/// Immutable decoded record storage belonging to one installed kernel image.
/// Cells initialize independently at first use, preserving record-error timing
/// while sharing results between readers and exact-image registry clones.
#[derive(Default)]
struct DecodedSpkCache {
    entries: Vec<OnceLock<Result<DecodedSpkSegment>>>,
}

impl std::fmt::Debug for DecodedSpkCache {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DecodedSpkCache")
            .field("segments", &self.entries.len())
            .finish()
    }
}

/// Validated payload of one descriptor, independent of subsequent evaluation.
#[derive(Debug)]
enum DecodedSpkSegment {
    Chebyshev(ChebSegment),
    Emda(EmdaSegment),
}

impl PartialEq for SpkRegistry {
    fn eq(&self, other: &Self) -> bool {
        // Cache initialization is not part of the stored-image value contract.
        self.kernels == other.kernels
            && self.bodies == other.bodies
            && self.automatic_centers == other.automatic_centers
    }
}

/// One local-discovery result: the display body name plus the status string (`"registered"`,
/// `"already_registered"`, or `"error: <detail>"`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpkDiscoveryEntry {
    /// Display body name.
    pub body_name: String,
    /// `"registered"`, `"already_registered"`, or `"error: <detail>"`.
    pub status: String,
}

/// Recursive center provider over one [`SpkRegistry`] whose bodies are
/// NAIF ids ([`register_all_naif_targets`]).
///
/// `center_state(center, et)` evaluates body `center` SSB-relative
/// through the same provider, so multi-level chains (Moon → EMB → SSB)
/// resolve without the caller unrolling them. The barycenter short-
/// circuits inside [`eval_ssb`] and never reaches this provider; an
/// unregistered center surfaces [`Error::SpkUnregisteredBody`], never
/// invented data. A center cycle (impossible in a valid NAIF kernel)
/// fails with [`Error::SpkCorrupt`] after [`NaifBaseProvider::MAX_DEPTH`]
/// levels instead of overflowing the stack.
///
/// Single-threaded tooling use: the depth counter is a [`Cell`], so this
/// provider is `!Sync` by design.
pub struct NaifBaseProvider<'a> {
    registry: &'a SpkRegistry,
    depth: Cell<usize>,
}

// Exported SPK/minor-body source helpers and automatic

/// JPL Horizons API endpoint.
pub const HORIZONS_API_URL: &str = "https://ssd.jpl.nasa.gov/api/horizons.api";

/// Legacy/default asteroid offset.
pub const NAIF_ASTEROID_OFFSET_LEGACY: i64 = 2_000_000;

/// JPL Horizons SPK offset.
pub const NAIF_ASTEROID_OFFSET_HORIZONS: i64 = 20_000_000;

/// Coverage edge tolerances: the reported coverage is the
/// *usable* span — its start sits one light-time inside the stored data
/// and Horizons reports boundaries at sub-day precision.
pub const SPK_COVERAGE_START_TOLERANCE_DAYS: f64 = 1.5;
/// See [`SPK_COVERAGE_START_TOLERANCE_DAYS`].
pub const SPK_COVERAGE_END_TOLERANCE_DAYS: f64 = 0.5;

/// Center assumed for a type-21 kernel whose reader exposes no segment
/// summaries: JPL Horizons small-body kernels are heliocentric
pub const TYPE21_DEFAULT_CENTER: i32 = 10;

/// Upper bound (AU) on the Earth-to-center distance for the Sun or the
/// barycenter.
pub const EARTH_REACH_FROM_SUN_OR_SSB_AU: f64 = 1.1;

/// Per-planetary-system aphelion bounds (AU), rounded up to the next
/// 0.1 AU, keyed by system number.
pub const PLANET_SYSTEM_APHELION_AU: &[(i32, f64)] = &[
    (1, 0.5),
    (2, 0.8),
    (3, 1.1),
    (4, 1.7),
    (5, 5.5),
    (6, 10.2),
    (7, 20.2),
    (8, 30.4),
    (9, 49.4),
];

/// Offset (days) of the light-time distance probe inside the first
/// segment start.
pub const TYPE21_PROBE_OFFSET_DAYS: f64 = 1.0e-6;

/// Half-step (days) of the legacy `FLG_SPEED` central difference in the
/// direct type-21 path.
pub const TYPE21_LEGACY_SPEED_HALF_STEP_DAYS: f64 = 1.0 / 86400.0;

/// Speed of light in AU/day.
pub const LIGHT_SPEED_AU_PER_DAY: f64 = 173.1446326846693;

/// JPL Horizons SPK generation window. The request-range clamp that consumes them is
/// download-flow work and stays in the download layer.
pub const HORIZONS_SPK_DATE_MIN: &str = "1600-01-01";
/// See [`HORIZONS_SPK_DATE_MIN`].
pub const HORIZONS_SPK_DATE_MAX: &str = "2500-01-01";
/// See [`HORIZONS_SPK_DATE_MIN`].
pub const HORIZONS_SPK_JD_MIN: f64 = 2_305_448.5;
/// See [`HORIZONS_SPK_DATE_MIN`].
pub const HORIZONS_SPK_JD_MAX: f64 = 2_634_167.5;

/// Default auto-download SPK window.
pub const DEFAULT_AUTO_SPK_START_DATE: &str = "2000-01-01";
/// See [`DEFAULT_AUTO_SPK_START_DATE`].
pub const DEFAULT_AUTO_SPK_END_DATE: &str = "2100-01-01";

/// Per-tier SPK date windows:
/// `(tier name, start, end)`.
pub const TIER_SPK_DATE_RANGES: &[(&str, &str, &str)] = &[
    ("base", "1850-01-01", "2150-01-01"),
    ("medium", "1900-01-01", "2100-01-01"),
    ("extended", "1600-01-01", "2500-01-01"),
];

/// Environment variable selecting automatic SPK downloads
pub const AUTO_SPK_ENV_VAR: &str = "LIBEPHEMERIS_AUTO_SPK";

/// An SPK-file request to the JPL Horizons API: `center`
/// defaults to `"500@0"` (solar-system barycenter, for JPL/DE
/// compatibility).
pub struct HorizonsSpkRequest<'a> {
    /// Target body: number, name, or combined form (e.g. `"2060 Chiron"`).
    pub body: &'a str,
    /// Start date (`YYYY-MM-DD`).
    pub start: &'a str,
    /// End date (`YYYY-MM-DD`).
    pub end: &'a str,
    /// Reference center (default `"500@0"`).
    pub center: &'a str,
}

// ---------------------------------------------------------------------------
// Horizons SPK download composition (API-07 U16g; release `spk.download_spk`
// Network transport and file publication are separate owned operations.
// ---------------------------------------------------------------------------

/// Release `User-Agent` for Horizons SPK requests (`spk.download_spk`).
pub const HORIZONS_SPK_USER_AGENT: &str = "libephemeris/0.1";

/// Exchange attempts for the Horizons SPK retry loop. Retries are
/// immediate: the release 1 s inter-attempt pause is an interactive
/// courtesy, not a protocol requirement, and an engine must not sleep.
pub const HORIZONS_SPK_ATTEMPTS: u32 = 3;

/// One tier-download plan row (pure half of release
/// `spk_auto.ensure_all_ephemerides` step 2): the library body key, the
/// Horizons/NAIF pair from [`crate::constants::SPK_BODY_NAME_MAP`], the
/// JPL block reason (`None` when downloadable), and the cache filename
/// under the tier window. No filesystem or socket work.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnsurePlanEntry {
    /// Library body key.
    pub ipl: i32,
    /// Horizons query id from the map.
    pub horizons_id: &'static str,
    /// NAIF id from the map.
    pub naif_id: i32,
    /// JPL block reason from [`crate::constants::SPK_AUTO_DOWNLOAD_BLOCKED`].
    pub blocked: Option<&'static str>,
    /// Cache filename under the tier window.
    pub filename: String,
}

/// Per-body verdict of the ensure sweep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnsureVerdict {
    /// Already registered (or the tier kernel already on disk).
    Cached,
    /// Fetched and registered by this call.
    Downloaded,
    /// JPL-blocked or failed: the body keeps its Keplerian fallback.
    Fallback,
}

/// One body's ensure outcome: the release display name plus the
/// machine-readable verdict.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnsureEntryOutcome {
    /// Display name.
    pub body_name: String,
    /// Machine-readable verdict.
    pub verdict: EnsureVerdict,
    /// Release `"fallback (...)"` reason, else `""`.
    pub detail: String,
}

/// Verified ensure report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnsureSummary {
    /// Tier the sweep ran for.
    pub tier: String,
    /// Tier kernel filename.
    pub ephemeris_file: String,
    /// `"cached"` when the kernel sits in `data_dir`, else
    /// `"missing: <file>"`.
    pub ephemeris_status: String,
    /// Per-body outcomes in map order.
    pub entries: Vec<EnsureEntryOutcome>,
    /// Bodies already registered.
    pub cached: u32,
    /// Bodies fetched and registered by this call.
    pub downloaded: u32,
    /// JPL-blocked or failed bodies.
    pub fallback: u32,
    /// Missed tier-kernel errors.
    pub errors: u32,
}

/// Configuration for automatic SPK download of one minor body
/// `get_cache_dir` creation and the actual
/// MD5 `{stem}_{hash[:8]}.bsp` filename stay in the download layer,
/// which owns both the hash dependency and the directory writes.
pub struct AutoSpkConfig {
    /// LibEphemeris body ID (e.g. `CHIRON`).
    pub ipl: i32,
    /// JPL Horizons identifier (e.g. `"2060"`).
    pub body_id: String,
    /// NAIF ID for SPK lookup (`None` until resolved from the file).
    pub naif_id: Option<i32>,
    /// Coverage start (`YYYY-MM-DD`).
    pub start: String,
    /// Coverage end (`YYYY-MM-DD`).
    pub end: String,
    /// Cache directory override (`None` defers to the SPK cache
    /// resolution).
    pub cache_dir: Option<String>,
    /// Whether auto-download is enabled for this body.
    pub enabled: bool,
    /// Downloaded SPK path (`None` until downloaded).
    pub spk_path: Option<String>,
}

/// Major-asteroid SPK download info:
/// `(body id, asteroid number, Horizons id, NAIF id, display name)`.
/// Bodies formerly in JPL's major-body index use the `"Name;"` syntax
/// to bypass the Horizons pre-computed-object restriction.
pub const MAJOR_ASTEROID_SPK_INFO: &[(i32, i32, &str, i32, &str)] = &[
    (
        crate::constants::CERES,
        1,
        "Ceres;",
        crate::constants::NAIF_CERES,
        "Ceres",
    ),
    (
        crate::constants::PALLAS,
        2,
        "Pallas;",
        crate::constants::NAIF_PALLAS,
        "Pallas",
    ),
    (
        crate::constants::JUNO,
        3,
        "Juno;",
        crate::constants::NAIF_JUNO,
        "Juno",
    ),
    (
        crate::constants::VESTA,
        4,
        "Vesta;",
        crate::constants::NAIF_VESTA,
        "Vesta",
    ),
    (
        crate::constants::CHIRON,
        2060,
        "2060",
        crate::constants::NAIF_CHIRON,
        "Chiron",
    ),
];

/// Human-readable SPK diagnostic names.
pub const SPK_BODY_NAMES: &[(i32, &str)] = &[
    (15, "Chiron"),
    (16, "Pholus"),
    (17, "Ceres"),
    (18, "Pallas"),
    (19, "Juno"),
    (20, "Vesta"),
    (146199, "Eris"),
    (100377, "Sedna"),
    (146108, "Haumea"),
    (146472, "Makemake"),
    (38978, "Ixion"),
    (100482, "Orcus"),
    (60000, "Quaoar"),
    (17066, "Nessus"),
    (18405, "Asbolus"),
    (20199, "Chariklo"),
    (235088, "Gonggong"),
    (109942, "Apophis"),
    (10010, "Hygiea"),
    (10433, "Eros"),
];

/// Hand-built kernel-image pieces (test-only, module scope so
/// sibling-module tests can reuse them through
/// [`moon_test_kernel_bytes`]).
#[cfg(test)]
#[derive(Clone, Copy)]
pub(crate) struct SegSpec {
    start: f64,
    end: f64,
    target: i32,
    center: i32,
    frame: i32,
    seg_type: i32,
    name: &'static str,
}

// Private implementation boundaries; public logical paths remain unchanged.
mod acquisition;
pub use acquisition::build_horizons_url;
pub use acquisition::download_and_register_spk_to;
pub use acquisition::download_and_register_spk_with_policy;
pub use acquisition::download_horizons_spk_to;
pub use acquisition::ensure_all_ephemerides_to;
pub use acquisition::ensure_all_plan;
pub use acquisition::extract_json_string_field;
pub use acquisition::has_json_key;
pub use acquisition::spk_download_dest;
pub use acquisition::spk_download_filename;
use acquisition::text_char_at;
mod decode;
pub use decode::check_frame;
use decode::declared_endian;
pub use decode::decode_base64_standard;
pub use decode::decode_chebyshev;
pub use decode::decode_emda;
use decode::decode_json_string;
pub use decode::parse_horizons_spk_response;
mod range_presentation;
pub use range_presentation::{OwnedSpkRange, RangeImage};
mod file_image_cache;
mod reader;
pub use reader::Type21ReaderCache;
#[cfg(test)]
use reader::build_image;
pub use reader::candidate_spk_paths;
pub use reader::classify_center;
pub use reader::declared_center;
pub use reader::deduce_naif_id;
pub use reader::detect_naif_id_from_targets;
pub use reader::discover_local_spks;
pub use reader::et_to_jd_tdb;
pub use reader::extract_asteroid_number;
pub use reader::find_covering_spk;
pub use reader::find_naif_id_for_asteroid;
pub use reader::get_major_asteroid_info;
pub use reader::icrs_to_ecliptic_j2000;
pub use reader::is_major_asteroid;
pub use reader::is_spk_downloadable;
pub use reader::iso_to_jd;
pub use reader::jd_tdb_to_et;
pub use reader::jd_to_iso_date;
#[cfg(test)]
use reader::put_f64;
#[cfg(test)]
use reader::put_i32;
use reader::record_offset;
pub use reader::register_all_naif_targets;
#[cfg(test)]
use reader::request_label_to_et;
pub use reader::resolve_auto_spk_download;
pub use reader::resolve_spk_path;
pub use reader::sanitize_filename;
pub use reader::spk_body_name;
pub use reader::spk_targets_of;
pub use reader::type21_center_reach_ok;
pub use reader::type21_default_stencil_days;
pub use reader::type21_earth_reach_au;
pub use reader::type21_request_stencil_days;
pub use reader::type21_usable_span;
mod index;
pub use index::auto_spk_cache_stem;
pub use index::generate_spk_cache_filename;
pub use index::is_in_jpl_major_body_index;
#[cfg(test)]
pub(crate) use index::moon_chain_test_kernel_bytes;
#[cfg(test)]
pub(crate) use index::moon_test_kernel_bytes;
pub use index::resolve_chain;
pub use index::resolve_spk_cache_dir;
pub use index::spk_coverage_covers_request;
pub use index::spk_date_range_for_tier;
pub use index::spk_date_range_for_tier_name;
#[cfg(feature = "nbody")]
mod positional;
#[cfg(feature = "nbody")]
pub(crate) use positional::PositionalChebyshev;
mod evaluation;
// Selected force-file ownership is bounded for ordinary N-body serving.
#[cfg(any(test, feature = "nbody"))]
mod bounded;
#[cfg(any(test, feature = "nbody"))]
pub(crate) use bounded::ChebyshevFile;
#[cfg(all(test, feature = "nbody"))]
pub(crate) use bounded::tests::image as force_test_image;
mod geometric_source;
mod subset;
pub use evaluation::emda_record_index;
pub use evaluation::emda_record_reference;
pub use evaluation::eval_chebyshev;
pub use evaluation::eval_chebyshev_two_part_jd;
pub use evaluation::eval_emda;
pub use evaluation::eval_relative;
pub use evaluation::eval_relative_two_part_jd;
pub use evaluation::eval_ssb;
pub use evaluation::eval_ssb_two_part_jd;
pub use geometric_source::GeometricSpkSource;
pub use subset::{SpkSubsetInput, SpkSubsetReport, write_spk_subset};

#[cfg(test)]
#[path = "tests/tests.rs"]
mod tests;
#[cfg(test)]
pub(crate) use tests::geometric_minor_test_image;
#[cfg(test)]
pub(crate) use tests::range_capture_chain_image;
