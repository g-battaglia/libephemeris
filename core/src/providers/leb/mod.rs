//! LEB/LEB2 container decoding, evaluation, and coverage policy.
//!
//! Reads project-native LEB1/LEB2 metadata, decodes compressed coefficient
//! chunks and evaluates Chebyshev body/nutation channels and auxiliary tables.
//! Wire integers are little-endian and coefficient words binary64; offsets
//! are absolute bytes from the container start. Layouts are described by the
//! header/entry constants here and `docs/data-format.md`.
//!
//! [`LebFile`] binds an index and cache to shared immutable backing. Its safe
//! owned constructor snapshots bytes; mapped and positional constructors have
//! explicit immutable-inode preconditions. Raw [`FileIndex`] methods instead
//! accept caller-held bytes, which must match the image used for decoding.
//! Public caches also require the image-affinity contract of [`LebChunkCache`].
//!
//! Checked byte arithmetic detects truncation, inflated counts and invalid
//! offsets. Numerical evaluation uses the same segment-selection core for
//! scalar and two-part TT epochs. Decompression and polynomial evaluation are
//! lazy, so structural admission is not complete payload validation.
//!
//! # Ordered layers
//!
//! * Container decode: [`LebFormat`] detection, [`FileHeader`] metadata,
//!   and the [`SectionEntry`] directory with last-wins resolution.
//! * Body and chunk index: [`BodyEntry`] table validation with
//!   [`BodyEntry::expected_uncompressed_size`], and [`ChunkEntry`] chunk
//!   lookup with [`find_chunk`] / [`select_chunk`].
//! * Blob decoding: shuffle lanes, coefficient-major reorder, and bounded
//!   zstd decode over caller-held slices.
//! * Scalar evaluation: [`FileIndex::eval_body`] for monolithic bodies,
//!   [`FileIndex::eval_leb2_body`] over decoded blobs, and the two-part
//!   entry shapes sharing the [`select_segment`] core with [`NormalizedTime`]
//!   `tau`.
//! * Auxiliary reads: [`FileIndex::eval_nutation`], [`FileIndex::delta_t`],
//!   and [`FileIndex::get_star`] over the [`NutationHeader`],
//!   [`DeltaTHeader`], and [`StarEntry`] sections.
//! * Groups and composition: [`CompositeIndex`] / [`TieredIndex`] dispatch
//!   with body-group tables, and system-barycenter composition.
//! * Coverage and tier policy: [`CoverageReport`] date-aware coverage,
//!   [`SealedOutcome`] sealed classification, and tier eligibility with
//!   runtime asset checks.
//! * Inspection and verification: index inspection, container verification
//!   ([`LebVerifyReport`]), and deterministic LEB2-vs-LEB1 comparison
//!   ([`Leb2CompareReport`]).
//! * Lifecycle: [`LebFile`] owned images with advisory warm/cool ranges,
//!   opened via [`open_leb`].
//!
//! # Frames, units, epochs
//!
//! Integers are little-endian and floats binary64; offsets are absolute
//! byte positions. Epochs are Julian Days (TT at the evaluation entries);
//! segment geometry is `jd_start` plus `interval_days` per segment.
//! Returned states carry one entry per stored component: AU position with
//! AU/day velocity for ICRS bodies, and degrees/degrees/AU with
//! (degrees/day, degrees/day, AU/day) for ecliptic bodies. Ecliptic longitudes
//! wrap with `rem_euclid(360.0)`.
//!
//! # Refusals
//!
//! Truncated, inflated, or out-of-bounds bytes are [`Error::LebCorrupt`];
//! unknown magic, versions, and compression schemes keep their distinct
//! variants. Body/nutation channels refuse uncovered epochs rather than
//! extrapolating; stored Delta-T interpolation instead clamps to endpoint
//! values. Missing bodies, tier mismatches and sealed refusals retain typed
//! errors. Checked numerical entries reject nonfinite inputs/data as documented;
//! the raw [`eval_delta_t`] helper does not validate caller-supplied table words.
//!

use crate::errors::{Error, Result};
use crate::math::chebyshev::{chebyshev_deriv_coeffs, chebyshev_value};
use crate::types::{NormalizedTime, TwoPartEpoch};
use std::collections::HashMap;
use std::sync::Arc;

mod convert;
pub(crate) use convert::validate_leb2_conversion_target;
pub use convert::{Leb2ConversionReport, Leb2ConvertOptions, convert_leb1_to_leb2};

#[cfg(feature = "bundled-data")]
mod bundled;
#[cfg(feature = "bundled-data")]
pub use bundled::bundled_base_core;

/// Common header size in bytes (`<4sIIIdddI20s`).
pub const HEADER_SIZE: usize = 64;
/// Section directory entry size in bytes (`<IIQQ`).
pub const SECTION_DIR_SIZE: usize = 24;
/// LEB1 body index entry size in bytes (`<iIIdddIIQ`).
pub const LEB1_BODY_ENTRY_SIZE: usize = 52;
/// LEB2 body index entry size in bytes (`<iIIdddIIQQQ`).
pub const LEB2_BODY_ENTRY_SIZE: usize = 68;
/// LEB2 v2 chunk index header size in bytes (`<IIfI`: count, reserved,
/// nominal chunk interval in days as `f32`, reserved).
pub const CHUNK_INDEX_HEADER_SIZE: usize = 16;
/// LEB2 v2 chunk index entry size in bytes
/// (`<ddQQQII`: JD start/end, blob offset, compressed/uncompressed sizes,
/// first segment, segment count).
pub const CHUNK_ENTRY_SIZE: usize = 48;
/// Byte-shuffle lane width: one lane per byte of a binary64 word.
const SHUFFLE_LANES: usize = 8;
/// Nutation section header size in bytes (`<dddIIII`).
pub const NUTATION_HEADER_SIZE: usize = 40;
/// Delta-T section header size in bytes (`<II`).
pub const DELTA_T_HEADER_SIZE: usize = 8;
/// Delta-T table entry size in bytes (`<dd`: JD plus Delta-T in days).
pub const DELTA_T_ENTRY_SIZE: usize = 16;
/// Fixed-star catalog entry size in bytes (`<iddddddd4s`).
pub const STAR_ENTRY_SIZE: usize = 64;

/// Magic of monolithic raw-coefficient files.
pub const MAGIC_LEB1: [u8; 4] = *b"LEB1";
/// Magic of compressed per-body files.
pub const MAGIC_LEB2: [u8; 4] = *b"LEB2";
/// Only accepted LEB1 version.
pub const VERSION_LEB1: u32 = 1;
/// Legacy LEB2 version: one zstd blob per body.
pub const VERSION_LEB2_V1: u32 = 1;
/// Chunked LEB2 version: one zstd blob per ~10-year chunk.
pub const VERSION_LEB2_V2: u32 = 2;
/// Only accepted LEB2 header `flags`: mantissa truncation, coefficient
/// reorder, byte shuffle, then zstd.
pub const COMPRESSION_ZSTD_TRUNC_SHUFFLE: u32 = 1;

/// Section id of the body index (both formats).
pub const SECTION_BODY_INDEX: u32 = 0;
/// Section id of raw Chebyshev data (LEB1 only).
pub const SECTION_CHEBYSHEV: u32 = 1;
/// Section id of the nutation table.
pub const SECTION_NUTATION: u32 = 2;
/// Section id of the Delta-T table.
pub const SECTION_DELTA_T: u32 = 3;
/// Section id of the fixed-star catalog.
pub const SECTION_STARS: u32 = 4;
/// Section id of orbital elements (reserved, absent at the pin).
pub const SECTION_ORBITAL_ELEMENTS: u32 = 5;
/// Section id of per-body compressed blobs (LEB2 only).
pub const SECTION_COMPRESSED_CHEBYSHEV: u32 = 6;

/// Coordinate type: ICRS barycentric planet center, AU xyz.
pub const COORD_ICRS_BARY: u32 = 0;
/// Coordinate type: ecliptic of date, deg/deg/AU lon/lat/dist.
pub const COORD_ECLIPTIC: u32 = 1;
/// Coordinate type 2 is retired (pre-3.1.0 uranians). It still decodes —
/// legacy files in the wild encode it — but the calculation path never
/// serves it. Never reuse the value.
pub const COORD_HELIO_ECL_RETIRED: u32 = 2;
/// Coordinate type 3 is reserved and unused.
/// Reserved: never written or served at the pin (retrograde-cusp design note lives with the algorithms record, §16).
pub const COORD_GEO_ECLIPTIC: u32 = 3;
/// Coordinate type: ICRS system barycenter, AU xyz (runtime COB
/// correction is POS-pipeline work, not reader work).
pub const COORD_ICRS_BARY_SYSTEM: u32 = 4;

/// LEB container flavor established by magic, version, and flags.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LebFormat {
    /// `LEB1` version 1: monolithic raw-f64 coefficient data.
    Leb1,
    /// `LEB2` version 1: one zstd blob per body.
    Leb2V1,
    /// `LEB2` version 2: chunked per-body zstd blobs.
    Leb2V2,
}

/// Decoded 64-byte common file header.
#[derive(Debug, Clone, PartialEq)]
pub struct FileHeader {
    magic: [u8; 4],
    version: u32,
    section_count: u32,
    body_count: u32,
    jd_start: f64,
    jd_end: f64,
    generation_epoch: f64,
    flags: u32,
    format: LebFormat,
}

/// Decoded 24-byte section directory entry: absolute `offset` and `size`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SectionEntry {
    section_id: u32,
    offset: u64,
    size: u64,
}

/// Decoded body index entry, unified across flavors.
///
/// `data_offset` is the absolute offset of the first raw coefficient
/// (LEB1) or of the compressed blob / chunk index (LEB2); the
/// `compressed_size` pair is `None` for LEB1.
#[derive(Debug, Clone, PartialEq)]
pub struct BodyEntry {
    body_id: i32,
    coord_type: u32,
    segment_count: u32,
    jd_start: f64,
    jd_end: f64,
    interval_days: f64,
    degree: u32,
    components: u32,
    data_offset: u64,
    compressed_size: Option<u64>,
    uncompressed_size: Option<u64>,
}

/// Header, section directory, and unique body index of one container.
///
/// Repeated body ids resolve to the last physical entry while retaining their
/// first insertion position. The header's body count remains the physical
/// record count. Composite and tier dispatch retain their separate reader
/// precedence rules.
///
/// The index owns metadata, not the coefficient image. Every raw byte argument
/// must come from the same immutable container used to decode this index.
/// Cached evaluation additionally requires a cache dedicated to that image;
/// [`LebChunkCache`] explains why clearing is not an image reset. [`LebFile`]
/// maintains this pairing for file-backed callers.
#[derive(Debug, Clone, PartialEq)]
pub struct FileIndex {
    header: FileHeader,
    sections: Vec<SectionEntry>,
    bodies: Vec<BodyEntry>,
}

/// Decoded 40-byte nutation section header (`<dddIIII`).
///
/// `segment_count == 0` (the `--skip-aux` generator mode) means the table
/// is absent and must not be evaluated. A zero `interval_days` with a
/// nonzero count is an evaluation-time failure at the body-evaluation layer, not a decode
/// failure, matching the release `eval_nutation` order.
#[derive(Debug, Clone, PartialEq)]
pub struct NutationHeader {
    jd_start: f64,
    jd_end: f64,
    interval_days: f64,
    degree: u32,
    components: u32,
    segment_count: u32,
}

/// Decoded 8-byte Delta-T section header (`<II`: entry count, reserved).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DeltaTHeader {
    n_entries: u32,
}

/// Decoded 64-byte fixed-star catalog entry
/// (`<iddddddd4s`; the trailing 4 bytes are reserved).
#[derive(Debug, Clone, PartialEq)]
pub struct StarEntry {
    star_id: i32,
    ra_j2000_deg: f64,
    dec_j2000_deg: f64,
    pm_ra_deg_per_yr: f64,
    pm_dec_deg_per_yr: f64,
    parallax_arcsec: f64,
    rv_km_per_s: f64,
    magnitude: f64,
}

/// Decoded 16-byte LEB2 v2 chunk index header (`<IIfI`).
///
/// `chunk_count` is the entry count guarded by [`validate_entry_count`];
/// `chunk_interval_days` is the nominal chunk length (10 years at
/// generation) carried as `f32`, advisory only — the segment index owns
/// the epoch, never this rounded value.
#[derive(Debug, Clone, PartialEq)]
pub struct ChunkIndexHeader {
    chunk_count: u32,
    chunk_interval_days: f32,
}

/// Decoded 48-byte LEB2 v2 chunk index entry (`<ddQQQII`).
///
/// Each chunk covers a temporal slice of one body's coefficients and
/// decodes independently through the §1.5 pipeline. `blob_offset` is
/// absolute; `segment_start`/`segment_count` address the body's global
/// segment table. Chunk JD bounds are advisory (write-time rounded):
/// lookup binary-searches them, then ownership is verified against the
/// segment table (see [`select_chunk`]).
#[derive(Debug, Clone, PartialEq)]
pub struct ChunkEntry {
    jd_start: f64,
    jd_end: f64,
    blob_offset: u64,
    compressed_size: u64,
    uncompressed_size: u64,
    segment_start: u32,
    segment_count: u32,
}

// COB correction do not exist yet natively ; the composition below

/// LEB1 merged-main partition, in generation + merge order
/// (`leb_groups.LEB1_GROUPS`).
/// Single source of truth for every consumer (generators, merge, downloader, wizard, diagnostics); generator-dict keys must match this tuple.
pub const LEB1_GROUPS: [&str; 4] = ["planets", "asteroids", "exotics", "analytical"];
/// Standalone LEB1 companions generated and verified but never merged.
/// Empty since 3.1.0: the `uranians` companion was retired and fictitious
/// bodies are always served by runtime analytical models.
pub const LEB1_COMPANION_GROUPS: [&str; 0] = [];
/// LEB2 distribution partition: `{tier}_{group}.leb2` companion files
/// (`leb_groups.LEB2_GROUPS`).
/// `uranians` was a fifth group until 3.1.0; legacy files with that name are still recognized but never consulted.
pub const LEB2_GROUPS: [&str; 4] = ["core", "asteroids", "exotics", "apogee"];
/// Filename group suffixes recognized by the composite heuristics,
/// including the legacy-only `uranians` name (classified, never opened).
pub const GROUP_SUFFIXES: [&str; 5] = ["core", "asteroids", "apogee", "uranians", "exotics"];
/// Retired group suffixes: recognized for classification but never
/// attached — pre-3.1.0 installs may still carry such companions on disk,
/// and opening one would resurface retired `source = "LEB"` coverage.
pub const RETIRED_GROUP_SUFFIXES: [&str; 1] = ["uranians"];

/// Reviewed LEB precision tier. Selection priority is always
/// base -> medium -> extended, independent of installation order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LebTier {
    /// Short-DE440-kernel solution, modern interval.
    Base,
    /// Extended DE440 solution.
    Medium,
    /// DE441 solution outside DE440 coverage.
    Extended,
}

/// First-wins composite over several decoded container indexes with
/// their borrowed container images.
///
/// Body dispatch is first-wins in
/// construction order; auxiliary channels (nutation, Delta-T, stars) come
/// from the first constituent that carries usable data; `jd_range` is the
/// widest file-level range. Support reads (Sun/Earth states consumed by
/// modelled bodies) use the same dispatch, so a missing support channel
/// surfaces [`Error::LebBodyNotFound`] naming the support id — never a
/// silent substitution.
///
/// Raw-index construction borrows its byte slices. File-backed construction
/// retains shared immutable image owners; file I/O, mmap and warm/cool lifecycle
/// remain with those owners. Retaining an index does not copy its buffers.
///
/// Each constituent uses an image-specific [`LebChunkCache`]: LEB2 values
/// remain bitwise-identical to direct decoding; LEB1 bypasses the memo.
/// [`Self::new`] owns private bounded caches; [`Self::from_files`] shares
/// the immutable files' bounded caches across provider rebuilds. Equal body
/// and chunk keys from different images never share entries. Tiered views
/// inherit the same caches through their composites.
pub struct CompositeIndex<'a> {
    entries: CompositeSources<'a>,
    caches: Vec<Arc<LebChunkCache>>,
    body_map: HashMap<i32, BodyLocation>,
    // One fixed allocation at construction keeps provider enum storage small.
    builtin_bodies: Box<[Option<BodyLocation>; BUILTIN_BODY_COUNT]>,
    nutation_at: Option<usize>,
    delta_t_at: Option<usize>,
    stars_at: Option<usize>,
}

/// Fixed declared Sun-through-interpolated-perigee identifier range. Arbitrary
/// dataset IDs never determine allocation size or change randomized lookup.
const BUILTIN_BODY_COUNT: usize = crate::constants::INTP_PERG as usize + 1;

/// Original source and resolved row within its immutable decoded body index.
#[derive(Clone, Copy)]
struct BodyLocation {
    /// Position of the first serving source in the constituent array.
    source_position: usize,
    /// Resolved row within that source, after duplicate-record handling.
    body_row: usize,
}

/// Raw views borrow bytes; shared views retain complete immutable image owners.
/// Keeping the owner, index and coefficient cache together avoids self references
/// and permits a reusable provider without copying file data or decoded words.
enum CompositeSources<'a> {
    Raw(Vec<(FileIndex, LebSource<'a>)>),
    Shared(Vec<LebFile>),
}

/// Date-aware selection over per-tier composites, independently for each body.
/// Priority is base, then medium, then extended, regardless of installation
/// order. If no stored interval covers an epoch, the first available reader
/// returns its range failure; the serving dispatcher owns any permitted
/// analytical-model route. No JPL/SPK image is opened by this index.
pub struct TieredIndex<'a> {
    tiers: Vec<(LebTier, CompositeIndex<'a>)>,
    candidates: HashMap<i32, Vec<TierCandidate>>,
}

/// Immutable dispatch metadata copied from the original serving body row.
/// Bounds keep their original words; no date selection or state is retained.
#[derive(Clone, Copy)]
struct TierCandidate {
    /// Position in the priority-ordered tier array.
    position: usize,
    /// Original per-body Julian-day bounds, absent if the reader has none.
    bounds: Option<(f64, f64)>,
}

/// Bodies whose LEB channels store the system barycenter
/// (`COORD_ICRS_BARY_SYSTEM`), needing a runtime center-of-body offset
/// to reach the planet center.
pub const SYSTEM_BARY_BODY_IDS: [i32; 5] = [5, 6, 7, 8, 9];

/// Centered-difference step for the center-offset velocity, in days.
/// Source: release `fast_calc._VEL_H = 0.0005`.
/// Shared numerical-differentiation step reused across the release pipeline (not a fitted astronomical constant): small enough to bound truncation, large enough to contain amplified input noise.
pub const CENTER_OFFSET_VELOCITY_STEP_DAYS: f64 = 0.0005;

// Retarded-epoch pipeline ; `FLG_TRUEPOS` works from the boundary

/// Library body id of the Sun, the LEB support state consulted by the
/// sealed companion-escalation rule.
pub const SUN_BODY_ID: i32 = 0;

/// First fictitious body id: ids `40..=58` are served by declared runtime
/// analytical models, so no LEB file — not even a hand-selected legacy
/// companion — provides coverage for them.
pub const FICTITIOUS_BODY_MIN: i32 = 40;
/// Last fictitious body id (see [`FICTITIOUS_BODY_MIN`]).
pub const FICTITIOUS_BODY_MAX: i32 = 58;

/// First Julian Day (TT) of the packaged apogee residual grid.
/// Source: generated `lunar_apse_model.APOGEE_CT_JD_START` @ the pin.
pub const APOGEE_MODEL_JD_START: f64 = 2286820.5;
/// Apogee residual-grid sampling step in days.
/// Source: generated `lunar_apse_model.APOGEE_CT_STEP_DAYS` @ the pin.
pub const APOGEE_MODEL_STEP_DAYS: f64 = 10.0;
/// Apogee residual-grid sample count.
/// Source: generated `lunar_apse_model.APOGEE_CT_COUNT` @ the pin.
pub const APOGEE_MODEL_COUNT: u64 = 40250;
/// Last Julian Day (TT) of the packaged apogee residual grid:
/// `START + STEP * (COUNT - 1)`.
pub const APOGEE_MODEL_JD_END: f64 =
    APOGEE_MODEL_JD_START + APOGEE_MODEL_STEP_DAYS * ((APOGEE_MODEL_COUNT - 1) as f64);
/// First Julian Day (TT) of the packaged perigee residual grid.
/// Source: generated `lunar_apse_model.PERIGEE_CT_JD_START` @ the pin.
pub const PERIGEE_MODEL_JD_START: f64 = 2286820.5;
/// Perigee residual-grid sampling step in days.
/// Source: generated `lunar_apse_model.PERIGEE_CT_STEP_DAYS` @ the pin.
pub const PERIGEE_MODEL_STEP_DAYS: f64 = 2.0;
/// Perigee residual-grid sample count.
/// Source: generated `lunar_apse_model.PERIGEE_CT_COUNT` @ the pin.
pub const PERIGEE_MODEL_COUNT: u64 = 201249;
/// Last Julian Day (TT) of the packaged perigee residual grid:
/// `START + STEP * (COUNT - 1)`.
pub const PERIGEE_MODEL_JD_END: f64 =
    PERIGEE_MODEL_JD_START + PERIGEE_MODEL_STEP_DAYS * ((PERIGEE_MODEL_COUNT - 1) as f64);

/// Coverage answer for one body, optionally resolved for a target date.
///
/// For a date-aware tiered reader a covered `jd` reports the concrete
/// selected tier's interval (`selected` names its priority position); a
/// date outside every stored interval reports the union coverage with
/// `selected` unset — never pretending one file served the request

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoverageReport {
    /// Stored interval start (JD).
    pub start: f64,
    /// Stored interval end (JD).
    pub end: f64,
    /// Priority position of the tier selected for the query date, when
    /// one covers it. `None` for date-blind answers and for union
    /// reports outside every stored interval.
    /// Positions index tier order (base, medium, extended).
    pub selected: Option<usize>,
}

/// Sealed-mode disposition of one body/date coverage question.
///
/// `Covered` is the only disposition that evaluates from LEB. Every other
/// disposition is decided without opening another provider:
/// `LocalModelFallthrough` and `NoLebChannel` leave the declared-model
/// routing decision to the provider-routing/event pipeline), while the
/// `Refused*` dispositions carry the typed error that must surface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SealedOutcome {
    /// The stored interval covers `jd`: evaluate from LEB.
    Covered,
    /// A curated companion outside its own window while LEB support
    /// (Sun) still covers the date: the declared local model may run.
    /// Never a silent substitution — the dispatcher records the source.
    LocalModelFallthrough,
    /// No LEB channel serves this body here (fictitious, derived, or
    /// companion ids without a stored interval, or no support coverage
    /// to classify against): not a range miss, so no range error is
    /// manufactured. Geometric evaluation surfaces
    /// [`Error::LebBodyNotFound`]; declared-model routing decides.
    NoLebChannel {
        /// Requested body id.
        body_id: i32,
    },
    /// A stored-state core body absent from the active reader
    /// (custom/partial file): a declared missing-body condition, not a
    /// date-range one. Surfaces [`Error::LebBodyNotFound`].
    MissingCoreChannel {
        /// Requested body id.
        body_id: i32,
    },
    /// The date is outside the applicable stored or fitted window: fail
    /// closed naming the requested body. Surfaces
    /// [`Error::LebOutOfRange`].
    RefusedRange {
        /// Requested body id.
        body_id: i32,
        /// Requested Julian Day (caller's own scale).
        jd: f64,
        /// Applicable window start (stored interval, or the fitted model
        /// window for 21/22, or the support window when the requested
        /// body has no channel of its own).
        start: f64,
        /// Applicable window end.
        end: f64,
    },
}

// Stored body coverage, cumulative tier requirements and inspection records.

/// Numerical origin stored in a selected LEB artifact, with the release
/// `_leb_precision_class` spellings (`"ephemeris"`, `"analytical"`,
/// `"numerical-model"`, `"mixed"`, `"unverified-local"`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PrecisionClass {
    /// Integrated ephemeris Chebyshev state (planets, asteroids, Moon).
    Ephemeris,
    /// Runtime analytical model (lunar nodes/apsides, fictitious bodies).
    Analytical,
    /// Extended-tier numerical trajectory of an N-body exotic.
    NumericalModel,
    /// Date-less union over a direct-SPK base/medium interval and the
    /// extended numerical trajectory: no single class applies.
    Mixed,
    /// The serving file is not a reviewed manifest artifact.
    UnverifiedLocal,
}

/// Lunar node/apse ids served by the runtime analytical model
/// (mean/true node, mean/osculating/interpolated apogee, interpolated perigee).
pub const ANALYTICAL_BODY_IDS: [i32; 6] = [10, 11, 12, 13, 21, 22];

/// Extended-tier bodies flown as massless N-body trajectories:
/// the extended exotic set minus the assist perturbers, which are windowed exactly like the classic asteroids.
/// The perturber set is built from its SPK window exactly like the
/// classic asteroids, so only the remainder takes the `"mixed"` /
/// `"numerical-model"` branch of [`precision_class`].
pub const NBODY_MODEL_BODY_IDS: [i32; 17] = [
    16, 10055, 10080, 10944, 11181, 17066, 18405, 20199, 30000, 38978, 60000, 100377, 100482,
    146108, 146199, 146472, 235088,
];

/// Coverage contract for one body from the active calculation source: `source` is always `"LEB"` (see
/// [`BodyCoverage::source`]), fictitious bodies (40-58) never receive a
/// record (every constructor returns `None` for them, like the release
/// early return), and the reported interval is the stored coefficient
/// range — an apparent-place request at the exact lower boundary still
/// needs up to one light-time of retardation inside it (the light-time layer).
#[derive(Debug, Clone, PartialEq)]
pub struct BodyCoverage {
    body_id: i32,
    precision_class: PrecisionClass,
    jd_start: f64,
    jd_end: f64,
    data_file: Option<String>,
    group: Option<String>,
    reviewed: bool,
}

/// File identity of one constituent of a composite or tiered reader: the
/// release `.path` / `._manifest_verified` attributes, carried
/// explicitly because native indexes own no file-system handles.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ServingFile<'s> {
    /// Serving file path (`None` labels a date-less union, never a file).
    pub path: Option<&'s str>,
    /// Whether the file is a reviewed manifest artifact.
    pub reviewed: bool,
}

/// Stored per-body interval of one file: the `bodies` row of release
/// `_reader_file_info`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BodyRange {
    body_id: i32,
    jd_start: f64,
    jd_end: f64,
}

/// File identity and body/range metadata of one open LEB container: `name` is the base name of `path`,
/// `group` follows [`group_from_path`], `size_bytes` is the caller-observed
/// file length (`None` when the release `os.path.getsize` fails), and
/// `bodies` lists every indexed body in id order with its stored bounds.
#[derive(Debug, Clone, PartialEq)]
pub struct FileInfo {
    name: String,
    path: String,
    group: Option<String>,
    size_bytes: Option<u64>,
    reviewed: bool,
    bodies: Vec<BodyRange>,
}

/// One body's row of a [`LebVerifyReport`]: the stored interval plus the
/// outcome of the deterministic evaluation probes. `error` is `None` on a
/// passing body and carries the first probe failure otherwise; a body is
/// never skipped silently.
#[derive(Debug, Clone, PartialEq)]
pub struct LebBodyVerify {
    body_id: i32,
    jd_start: f64,
    jd_end: f64,
    probes: u32,
    error: Option<String>,
}

/// Self-contained structural verification of one decoded LEB container:
/// the container flavor, one row per indexed body in id order, and the
/// aggregate verdict. This is the engine half of the `leb verify` bridge
/// (developer tooling): it proves the file decodes and every indexed body evaluates
/// over its stored interval with the native reader. Value verification
/// against an external reference is a
/// separate step, never implied by this report.
#[derive(Debug, Clone, PartialEq)]
pub struct LebVerifyReport {
    format: LebFormat,
    bodies: Vec<LebBodyVerify>,
}

/// One compared body's row of a [`Leb2CompareReport`]: the deterministic
/// sample count and the worst absolute component difference observed
/// between candidate and reference over the shared coverage. Units are the
/// body's native stored components (AU xyz / AU-day for ICRS bodies,
/// deg/deg/AU and (deg/day, deg/day, AU/day) for ecliptic bodies) — the same
/// units the release `leb2 verify` reports in.
#[derive(Debug, Clone, PartialEq)]
pub struct Leb2BodyComparison {
    body_id: i32,
    samples: u32,
    max_abs_err: f64,
}

/// One skipped body's row of a [`Leb2CompareReport`]: bodies are skipped,
/// never silently dropped, with the reason recorded.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Leb2BodySkipped {
    body_id: i32,
    reason: String,
}

/// LEB2-vs-LEB1 sampled comparison of two decoded container images: the
/// engine half of the `leb2 verify` bridge (developer tooling). Bodies present in
/// both containers are evaluated at a deterministic uniform grid over
/// their shared coverage (endpoints included, so segment-clamp
/// boundaries are always probed) through each container's own flavor
/// path, and the worst absolute native-component difference is reported
/// per body. Bodies that cannot be compared become [`Leb2BodySkipped`]
/// rows. Either image failing to decode is `Err`; a per-body evaluation
/// failure aborts the comparison with `Err` as well — unlike structural
/// verification, a comparison with unevaluated samples would be evidence
/// of nothing.
#[derive(Debug, Clone, PartialEq)]
pub struct Leb2CompareReport {
    compared: Vec<Leb2BodyComparison>,
    skipped: Vec<Leb2BodySkipped>,
}

/// Upper bound for the per-body sample grid: 4096 epochs × two
/// evaluations stays a sub-second structural check on real companions
/// while resolving every segment of any shipped body several times over.
pub const LEB2_COMPARE_MAX_SAMPLES: usize = 4096;

/// One immutable file required by a sealed LEB runtime tier
/// (`name`, `kind`, `group`, `path`,
/// `sha256`); `kind` is always `"leb2"` (see
/// [`RuntimeDataRequirement::kind`]).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RuntimeDataRequirement {
    name: String,
    group: String,
    path: String,
    sha256: String,
}

// By `download_file`, which stays work).

/// Caller-observed state of one `{tier}_{group}.leb2` candidate: the
/// release `check_data_status` row (`exists`/`size`) plus the bytes needed
/// for the header parse when present.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ObservedFile<'s> {
    /// Manifest file name (`{tier}_{group}.leb2`).
    pub name: &'s str,
    /// File bytes when present (`None` reports the file missing).
    pub bytes: Option<&'s [u8]>,
    /// Caller-observed file length in bytes, if known.
    pub size_bytes: Option<u64>,
}

/// Validation outcome of one required asset: the release
/// `check_data_status` exists/size reporting crossed with the
/// `_is_valid_leb` parse gate.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AssetState {
    /// Present and header-parseable (reviewed in the asset-validation sense).
    Ready {
        /// Caller-observed file length in bytes, if known.
        size_bytes: Option<u64>,
        /// Number of indexed bodies.
        body_count: usize,
    },
    /// No bytes were observed for this requirement.
    Missing,
    /// Bytes were observed but [`FileIndex::decode`] refused them (the
    /// release reader-open rejection, with its message).
    Unreadable {
        /// Decode failure detail.
        what: String,
    },
}

/// Per-file validation row for one tier requirement.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssetStatus {
    name: String,
    group: String,
    path: String,
    state: AssetState,
}

/// Tier validation report over the cumulative requirement set: `ready` is
/// true only when every required asset is [`AssetState::Ready`] — a
/// partial tier never reports a full pass. `missing` / `unreadable` name
/// the failing requirements explicitly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TierAssets {
    tier: LebTier,
    files: Vec<AssetStatus>,
}

// `SpkRegistry`), for three documented reasons: (1) warm/cool are

/// Advisory page size for range merging, in bytes.
///
/// Source: release `leb_format._PAGE_SIZE` (`os.sysconf("SC_PAGE_SIZE")`
/// with a 4096 fallback; the supported macOS/Linux hosts report 4096 or
/// a multiple of it, and the merge rule only widens the advised span —
/// it never changes a decoded value).
pub const ADVISE_PAGE_SIZE: u64 = 4096;

/// A shared LEB image with owned/static/mapped or positional backing, its
/// decoded [`FileIndex`] and lifecycle (`open`, `warm`, `cool`, `close`).
///
/// Raw composite construction borrows byte slices; file-backed composites
/// retain shared image handles and caches. Both keep their images alive for
/// evaluation under their public borrow/ownership lifetimes. Positional reads
/// share the numerical kernels under [`IndexedLimits`] and the immutable-inode
/// contract; neither a read-only descriptor nor a cached digest proves that
/// another process cannot modify the inode.
///
/// `Clone` shares the immutable bytes/index and bounded decoded cache. Only
/// handle configuration (including closed state) is independent. `PartialEq`
/// compares the path and stored image, never cache counters or residency.
/// Indexed clones compare by shared owner; independent indexed opens do not
/// scan or hash file contents to establish equality.
/// [`Self::open`] keeps the safe whole-file snapshot contract. For immutable
/// release assets, [`Self::open_mapped`] avoids a full heap copy under its explicit
/// safety contract. Mapped file length is not resident or total process memory.
#[derive(Debug, Clone)]
pub struct LebFile {
    path: String,
    image: Option<Arc<LebImage>>,
}

/// Immutable container snapshot and its image-specific decoded-word memo.
/// Data and cache move together, preventing identical body/chunk keys from
/// different snapshots (even under the same path label) from sharing words.
struct LebImage {
    bytes: storage::ImageStorage,
    index: FileIndex,
    cache: Arc<LebChunkCache>,
    sha256: std::sync::OnceLock<String>,
}

impl std::fmt::Debug for LebImage {
    /// Show metadata, not an entire potentially multi-GiB byte buffer.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LebImage")
            .field("byte_len", &self.bytes.source().len())
            .field("mapped_file_bytes", &self.bytes.mapped_len())
            .field("index", &self.index)
            .finish_non_exhaustive()
    }
}

impl PartialEq for LebFile {
    /// Compare semantic container data; warming/cooling does not change equality.
    fn eq(&self, other: &Self) -> bool {
        self.path == other.path
            && match (&self.image, &other.image) {
                (None, None) => true,
                (Some(a), Some(b)) => {
                    // Comparing distinct images is intentionally a full-content
                    // operation and may fault their mapped pages. Clones do not scan.
                    Arc::ptr_eq(a, b)
                        || (a.index == b.index
                            && match (a.bytes.source(), b.bytes.source()) {
                                (LebSource::Slice(a), LebSource::Slice(b)) => a == b,
                                _ => false,
                            })
                }
                _ => false,
            }
    }
}

impl LebFile {
    /// Compare retained owners without reading image bytes or cache contents.
    pub(crate) fn shares_image_with(&self, other: &Self) -> bool {
        self.path == other.path
            && match (&self.image, &other.image) {
                (None, None) => true,
                (Some(a), Some(b)) => Arc::ptr_eq(a, b),
                _ => false,
            }
    }
}

// Private implementation boundaries; public logical paths remain unchanged.
mod access;
mod access_api;
mod packing;
mod resources;
pub(crate) use access::LebSource;
pub use access::{IndexedLimits, LebReadPolicy};
pub(crate) use resources::IndexedResources;
use resources::{CoefficientWords, TrackedVec};
mod storage;
pub(crate) use storage::HashCacheAdvice;
pub(crate) use storage::open_regular_file;
pub use storage::release_data_cache;
pub mod open_diagnostics;
mod reader;
pub use reader::apply_center_offset;
use reader::coeff_bytes;
pub use reader::compare_leb2_against_leb1;
pub use reader::decompress_body;
pub use reader::decompress_chunk;
pub use reader::decompress_monolith;
pub use reader::file_group_suffix;
pub use reader::find_chunk;
use reader::finite_offset;
pub use reader::group_from_path;
pub use reader::group_prefix;
pub use reader::has_fitted_model_window;
pub use reader::has_pinned_siblings;
pub use reader::inspect_bytes;
pub use reader::inspect_file;
pub use reader::interval_contains;
pub use reader::is_core_body;
use reader::is_ecliptic_coord;
pub use reader::is_fictitious_body;
pub use reader::is_nbody_model_body;
pub use reader::is_retired_group_file;
pub use reader::is_stored_state_body;
pub use reader::is_system_barycenter_channel;
pub use reader::merge_ranges;
pub use reader::open_leb;
pub use reader::precision_class;
use reader::prefault;
pub use reader::reorder_to_coeff_major;
use reader::reorder_words;
pub use reader::retard_epoch;
pub use reader::runtime_data_requirements;
pub use reader::sealed_outcome;
pub use reader::shuffle_bytes;
use reader::slice_at;
pub use reader::star_count;
pub use reader::system_barycenter_name;
pub use reader::unshuffle_bytes;
pub use reader::verify_container;
mod decode;
pub use decode::check_asset;
pub use decode::check_single_tier;
pub use decode::check_tier_subset;
pub use decode::decode_chunk_index;
use decode::decode_delta_t_table;
pub use decode::decode_section_table;
pub use decode::read_blob;
use decode::read_f32;
use decode::read_f64;
use decode::read_i32;
use decode::read_segment_coefficients;
use decode::read_u32;
use decode::read_u64;
pub use decode::validate_entry_count;
pub use decode::validate_tier_assets;
use decode::with_table_words;
mod image_evaluation;
mod index;
pub(crate) use image_evaluation::ImageBinding;
pub use index::eligible_tiers;
pub use index::file_tier;
pub use index::inspect_index;
pub use index::reorder_to_segment_major;
pub use index::select_chunk;
pub use index::select_segment;
pub use index::select_segment_two_part;
pub use index::{DEFAULT_CHUNK_CACHE_BYTES, LebChunkCache};
mod acceleration;
mod evaluation;
mod relative;
mod segment_domain;
use evaluation::batch_eval;
use evaluation::batch_eval_epochs;
pub use evaluation::eval_body_center;
pub use evaluation::eval_delta_t;
use evaluation::eval_nutation_words;
use evaluation::eval_verify_probe;
pub use evaluation::interpolated_apse_model_window;
use evaluation::{eval_one_segment, eval_one_segment_inline};
pub use relative::CartesianState;

#[cfg(test)]
#[path = "tests/tests.rs"]
mod tests;
