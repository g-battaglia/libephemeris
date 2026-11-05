//! Context-owned IERS row acquisition for model Delta-T.
//!
//! Queries use UT1 Julian days, with row coordinates in MJD and offsets in
//! seconds. Covered monthly Delta-T wins; otherwise finals rows supply
//! `TAI-UTC + 32.184 - (UT1-UTC)` under the row kernel's adjacent-day policy.
//! This layer does not change independent frame/observer timescales.
//!
//! Configuration copies share their existing acquisition cache, even while
//! unloaded. Changing a path does not invalidate rows. Clear starts a new lazy
//! cache; explicit reload and supplied-row replacement return new row owners.

use super::{
    DELTA_T_CACHE_FILENAME, DeltaTDataPoint, FINALS_CACHE_FILENAME, IersDataPoint,
    IersLoadedCounts, LEAP_SECONDS_CACHE_FILENAME, delta_t_iers_for_jd, parse_delta_t,
    parse_finals_2000a,
};
use crate::context::{
    EffectiveNetworkPolicy, NetworkPolicy, ProviderMode, effective_network_policy,
};
use crate::errors::{Error, Result};
use crate::time::{LeapSecondEntry, hardcoded_leap_seconds, parse_leap_seconds};
use std::collections::HashMap;
use std::fmt;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

/// Immutable owned finals, leap-second and observed Delta-T tables.
///
/// Parsing preserves the existing row kernels and their domains. Missing or
/// unreadable finals/monthly files give empty tables; missing or empty leap
/// inputs use the packaged leap table. The nonempty leap fallback is retained
/// even when both Earth-orientation tables are empty. Monthly interpolation
/// stops at its endpoints; finals queries retain the row kernel's single-sided
/// adjacent-day values, not an unrestricted extrapolation.
/// Separate file views retain permissive finals/leap values for language
/// boundaries; they never change checked row constructors or SDK interpolation.
#[derive(Debug, Clone, PartialEq)]
pub struct IersDataSnapshot {
    finals: Vec<IersDataPoint>,
    source_finals: Vec<(f64, f64)>,
    leaps: Vec<LeapSecondEntry>,
    observed: Vec<DeltaTDataPoint>,
    source_leap_fields: Vec<(f64, f64)>,
    source_observed_text: Option<String>,
    source_views: OnceLock<SourceFileViews>,
}

/// Completed public file fields after the boundary's stable leap-key sort.
/// Both tables belong to the same retained acquisition generation.
#[derive(Debug, Clone, PartialEq)]
struct SourceFileViews {
    leaps: Vec<(f64, f64)>,
    observed: Vec<DeltaTDataPoint>,
}

impl IersDataSnapshot {
    /// Parse optional cache contents in finals, leap-second, observed order.
    ///
    /// Text may have been decoded from Latin-1 bytes by the filesystem owner.
    /// The monthly C 04 layout consumes this snapshot's parsed leap table.
    /// Parsers sort their rows and preserve their existing duplicate policy.
    /// All supplied text is parsed/copied into owned storage; no input borrow,
    /// file handle or network client survives construction.
    pub fn from_texts(finals: Option<&str>, leaps: Option<&str>, observed: Option<&str>) -> Self {
        let source_finals = finals.map(parse_finals_file_view).unwrap_or_default();
        let finals = finals.map(parse_finals_2000a).unwrap_or_default();
        let source_leap_fields = leaps
            .map(parse_leap_file_view)
            .unwrap_or_else(leap_fallback_view);
        let source_observed_text = observed.map(str::to_owned);
        let leaps = leaps
            .map(parse_leap_seconds)
            .unwrap_or_else(hardcoded_leap_seconds);
        let observed = observed
            .map(|text| parse_delta_t(text, &leaps))
            .unwrap_or_default();
        Self {
            finals,
            source_finals,
            leaps,
            observed,
            source_leap_fields,
            source_observed_text,
            source_views: OnceLock::new(),
        }
    }

    /// Snapshot the three named local files, in their load order.
    ///
    /// File errors are local misses. Every byte has a Latin-1 scalar mapping;
    /// a non-UTF 8 comment must not discard an otherwise usable cache file.
    /// Directory creation/download policy belongs to the acquisition owner.
    pub fn read_in(directory: &Path) -> Self {
        fn text(path: PathBuf) -> Option<String> {
            std::fs::read(path)
                .ok()
                .map(|raw| raw.into_iter().map(char::from).collect())
        }
        let finals_text = text(directory.join(FINALS_CACHE_FILENAME));
        let source_finals = finals_text
            .as_deref()
            .map(parse_finals_file_view)
            .unwrap_or_default();
        let finals = finals_text
            .as_deref()
            .map(parse_finals_2000a)
            .unwrap_or_default();
        let leap_text = text(directory.join(LEAP_SECONDS_CACHE_FILENAME));
        let source_leap_fields = leap_text
            .as_deref()
            .map(parse_leap_file_view)
            .unwrap_or_else(leap_fallback_view);
        let leaps = leap_text
            .as_deref()
            .map(parse_leap_seconds)
            .unwrap_or_else(hardcoded_leap_seconds);
        let observed_text = text(directory.join(DELTA_T_CACHE_FILENAME));
        let observed = observed_text
            .as_deref()
            .map(|raw| parse_delta_t(raw, &leaps))
            .unwrap_or_default();
        Self {
            finals,
            source_finals,
            leaps,
            observed,
            source_leap_fields,
            source_observed_text: observed_text,
            source_views: OnceLock::new(),
        }
    }

    /// Reload the three local inputs in finals, leap, monthly order.
    ///
    /// An absent finals file retains preceding finals rows, as explicit reload
    /// does; a present unreadable file replaces them with an empty table. Missing
    /// leaps use the packaged fallback and missing monthly data becomes empty.
    /// Latin-1 decoding and existing parsers preserve their units and policies.
    fn reload_in(directory: &Path, previous: Option<&Self>) -> Self {
        Self::reload_paths(
            &super::iers_cache_filenames().map(|name| directory.join(name)),
            previous,
        )
    }

    /// Reparse the independently retained finals, leap and monthly paths.
    /// Reentrant boundary callbacks may change later paths; absent finals retain
    /// the current preceding rows rather than a pre-callback snapshot.
    fn reload_paths(paths: &[PathBuf; 3], previous: Option<&Self>) -> Self {
        fn text(path: &Path) -> Option<String> {
            std::fs::read(path)
                .ok()
                .map(|raw| raw.into_iter().map(char::from).collect())
        }
        let finals_path = &paths[0];
        let (finals, source_finals) = if finals_path.exists() {
            let raw = text(finals_path);
            (
                raw.as_deref().map(parse_finals_2000a).unwrap_or_default(),
                raw.as_deref()
                    .map(parse_finals_file_view)
                    .unwrap_or_default(),
            )
        } else {
            previous
                .map(|rows| (rows.finals.clone(), rows.source_finals.clone()))
                .unwrap_or_default()
        };
        let leap_text = text(&paths[1]);
        let source_leap_fields = leap_text
            .as_deref()
            .map(parse_leap_file_view)
            .unwrap_or_else(leap_fallback_view);
        let leaps = leap_text
            .as_deref()
            .map(parse_leap_seconds)
            .unwrap_or_else(hardcoded_leap_seconds);
        let observed_text = text(&paths[2]);
        let observed = observed_text
            .as_deref()
            .map(|raw| parse_delta_t(raw, &leaps))
            .unwrap_or_default();
        Self {
            finals,
            source_finals,
            leaps,
            observed,
            source_leap_fields,
            source_observed_text: observed_text,
            source_views: OnceLock::new(),
        }
    }

    /// Inclusive loaded finals envelope in JD UT1, or no finals rows.
    /// Ordered comparisons preserve the first row on equal values and do not
    /// replace an unordered first endpoint; no finite-epoch guard is imposed.
    fn finals_range(&self) -> Option<(f64, f64)> {
        let first = self.finals.first()?.mjd;
        let (mut lo, mut hi) = (first, first);
        for row in &self.finals[1..] {
            if row.mjd < lo {
                lo = row.mjd;
            }
            if row.mjd > hi {
                hi = row.mjd;
            }
        }
        Some((lo + 2_400_000.5, hi + 2_400_000.5))
    }

    /// Monthly envelope in JD UT1 using the existing sorted table endpoints.
    fn observed_range(&self) -> Option<(f64, f64)> {
        Some((
            self.observed.first()?.mjd + 2_400_000.5,
            self.observed.last()?.mjd + 2_400_000.5,
        ))
    }

    /// Loaded row counts; these describe inputs, not coverage at an epoch.
    pub fn counts(&self) -> IersLoadedCounts {
        IersLoadedCounts {
            finals: self.finals.len(),
            leap_seconds: self.leaps.len(),
            delta_t: self.observed.len(),
        }
    }

    /// Parsed leap-second rows in their retained source order. The slice is
    /// borrowed from this immutable snapshot; it cannot outlive its owner.
    /// Missing input has already selected the packaged nonempty fallback.
    pub fn leap_seconds(&self) -> &[LeapSecondEntry] {
        &self.leaps
    }

    /// Daily finals rows retained by this immutable generation, in parser order.
    /// Values use MJD and UT1−UTC seconds; this slice cannot outlive its owner.
    pub fn finals_rows(&self) -> &[IersDataPoint] {
        &self.finals
    }

    /// Finals file fields for permissive public queries: `(MJD, UT1−UTC s)`.
    /// These values may be nonfinite and are separate from checked SDK rows.
    /// Distinct NaN keys survive; equal keys keep their initial position and
    /// last value, including the single shared key for signed zero.
    pub fn source_finals_rows(&self) -> &[(f64, f64)] {
        &self.source_finals
    }

    /// Loaded finals file envelope in JD UT1, including nonfinite endpoints.
    /// Ordered comparisons begin with the first retained key; an unordered
    /// first key therefore stays the endpoint. Empty input returns None.
    pub fn source_finals_range(&self) -> Option<(f64, f64)> {
        let (mut lo, _) = *self.source_finals.first()?;
        let mut hi = lo;
        for &(mjd, _) in &self.source_finals[1..] {
            if mjd < lo {
                lo = mjd;
            }
            if mjd > hi {
                hi = mjd;
            }
        }
        Some((lo + 2_400_000.5, hi + 2_400_000.5))
    }

    /// Permissive leap file `(MJD, TAI−UTC s)` fields before boundary sorting.
    /// Nonfinite values remain data. Missing/unparsable files use the packaged
    /// fallback after the original format-selection and alternate-parser stages.
    pub fn source_leap_fields(&self) -> &[(f64, f64)] {
        &self.source_leap_fields
    }

    /// Whether this retained generation has its boundary-ordered file views.
    pub fn source_views_ready(&self) -> bool {
        self.source_views.get().is_some()
    }

    /// Publish an actual stable sort permutation once for this generation.
    /// Returns false if indices are not a permutation of the parsed leap rows;
    /// no publication then occurs. The supplied order must be the language's
    /// stable MJD-key sort, preserving its unordered-NaN comparison semantics.
    /// Monthly C 04 conversion uses that leap view and its unchanged sanity band.
    pub fn prepare_source_views(&self, order: &[usize]) -> bool {
        if self.source_views_ready() {
            return true;
        }
        if order.len() != self.source_leap_fields.len() {
            return false;
        }
        let mut seen = vec![false; order.len()];
        for &index in order {
            let Some(slot) = seen.get_mut(index) else {
                return false;
            };
            if *slot {
                return false;
            }
            *slot = true;
        }
        let leaps = order
            .iter()
            .map(|&index| self.source_leap_fields[index])
            .collect::<Vec<_>>();
        let observed = self
            .source_observed_text
            .as_deref()
            .map(|text| parse_observed_file_view(text, &leaps))
            .unwrap_or_default();
        let _ = self.source_views.set(SourceFileViews { leaps, observed });
        true
    }

    /// Boundary-ordered leap fields after [`Self::prepare_source_views`].
    /// Before preparation this is the retained file order, with no coercion.
    pub fn source_leap_seconds(&self) -> &[(f64, f64)] {
        self.source_views
            .get()
            .map(|views| views.leaps.as_slice())
            .unwrap_or(&self.source_leap_fields)
    }

    /// Monthly rows derived using the prepared permissive leap file view.
    /// Before boundary preparation the checked monthly rows remain available;
    /// public language queries prepare before consuming these fields.
    pub fn source_observed_rows(&self) -> &[DeltaTDataPoint] {
        self.source_views
            .get()
            .map(|views| views.observed.as_slice())
            .unwrap_or(&self.observed)
    }

    /// Public TAI−UTC seconds using the prepared permissive leap view.
    /// The pre-1972 rate model remains independent of loaded leap fields.
    /// Modern ordered comparisons preserve nonfinite keys and offsets.
    pub fn source_tai_utc(&self, mjd: f64) -> f64 {
        file_tai_utc(mjd, self.source_leap_seconds())
    }

    /// Monthly observed Delta-T rows retained by this immutable generation.
    /// Epochs are ascending MJD and values are TT−UT1 seconds. Duplicate epochs
    /// retain the parser's policy; the borrowed slice cannot outlive its owner.
    pub fn observed_rows(&self) -> &[DeltaTDataPoint] {
        &self.observed
    }

    /// Delta-T seconds at the original UT1 Julian Day, or no covering rows.
    ///
    /// Nonfinite epochs retain the checked row kernel's typed refusal. The
    /// model callback maps ordinary row errors to a miss at its fallback seam.
    pub fn delta_t_seconds(&self, jd: f64) -> Result<Option<f64>> {
        delta_t_iers_for_jd(jd, &self.finals, &self.observed, &self.leaps)
    }
}

/// Whitespace stripping for Latin-1 fields, including the separator controls.
fn file_trim(field: &str) -> &str {
    field.trim_matches(|c: char| c.is_whitespace() || ('\u{1c}'..='\u{1f}').contains(&c))
}

/// Float field grammar for Latin-1 file columns, retaining nan/inf and valid
/// decimal digit separators without admitting malformed underscore placement.
fn file_float(field: &str) -> Option<f64> {
    let field = file_trim(field);
    // A decimal underscore must be surrounded by digits. Removing only
    // those separators preserves the file float grammar, including nan/inf.
    let bytes = field.as_bytes();
    if bytes.iter().enumerate().any(|(i, byte)| {
        *byte == b'_'
            && (i == 0
                || i + 1 == bytes.len()
                || !bytes[i - 1].is_ascii_digit()
                || !bytes[i + 1].is_ascii_digit())
    }) {
        return None;
    }
    field.replace('_', "").parse().ok()
}

/// Parse permissive finals fields without weakening checked row types.
/// Latin-1 columns are character positions. Date fields must parse as integers;
/// missing or malformed MJD/UT1 fields skip the line. Error-column failures do
/// not discard a usable UT1 value. Values retain ordinary float nonfinite values.
fn parse_finals_file_view(text: &str) -> Vec<(f64, f64)> {
    let normalized = text.replace("\r\n", "\n").replace('\r', "\n");
    let mut rows: Vec<(f64, f64)> = Vec::new();
    let mut positions: HashMap<u64, usize> = HashMap::new();
    for line in normalized.split_inclusive('\n') {
        let columns: Vec<char> = line.chars().collect();
        if columns.len() < 68 {
            continue;
        }
        let field = |lo, hi| columns[lo..hi].iter().collect::<String>();
        if [(0, 2), (2, 4), (4, 6)]
            .iter()
            .any(|&(lo, hi)| file_trim(&field(lo, hi)).parse::<i32>().is_err())
        {
            continue;
        }
        let (Some(mjd), Some(ut1)) = (file_float(&field(7, 15)), file_float(&field(58, 68))) else {
            continue;
        };
        if mjd.is_nan() {
            rows.push((mjd, ut1));
            continue;
        }
        let key = if mjd == 0.0 { 0 } else { mjd.to_bits() };
        if let Some(&index) = positions.get(&key) {
            rows[index].1 = ut1;
        } else {
            positions.insert(key, rows.len());
            rows.push((mjd, ut1));
        }
    }
    rows
}

/// Packaged leap fields for the same missing/malformed-file fallback policy.
fn leap_fallback_view() -> Vec<(f64, f64)> {
    hardcoded_leap_seconds()
        .into_iter()
        .map(|row| (row.mjd(), row.tai_utc()))
        .collect()
}

/// Preserve IERS numeric/textual leap fields without checked finite admission.
/// Format detection and alternate/fallback stages precede boundary sorting;
/// the returned rows are in file order, including every unordered MJD.
fn parse_leap_file_view(text: &str) -> Vec<(f64, f64)> {
    fn iers(text: &str) -> Vec<(f64, f64)> {
        const MONTHS: [&str; 12] = [
            "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec",
        ];
        let mut rows = Vec::new();
        for line in text.lines().map(file_trim) {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts = line.split_whitespace().collect::<Vec<_>>();
            if parts.len() < 5 {
                continue;
            }
            let numeric = parts[1..4]
                .iter()
                .all(|word| !word.is_empty() && word.bytes().all(|b| b.is_ascii_digit()));
            let parsed = if numeric {
                file_float(parts[0]).zip(file_float(parts[4]))
            } else {
                let month = MONTHS
                    .iter()
                    .position(|&month| parts[1].to_ascii_lowercase().get(..3) == Some(month));
                match (
                    parts[0].parse::<i32>(),
                    month,
                    parts[2].parse::<i32>(),
                    file_float(parts[4].trim_end_matches('s')),
                ) {
                    (Ok(year), Some(month), Ok(day), Some(offset)) => Some((
                        crate::time::calendar_to_mjd(year, month as i32 + 1, day),
                        offset,
                    )),
                    _ => None,
                }
            };
            if let Some(row) = parsed {
                rows.push(row);
            }
        }
        rows
    }
    let ietf = || {
        crate::time::parse_leap_seconds_ietf(text)
            .into_iter()
            .map(|row| (row.mjd(), row.tai_utc()))
            .collect::<Vec<_>>()
    };
    let is_ietf = ["NTP", "#$", "#@", "leap-seconds.list"]
        .iter()
        .any(|marker| text.contains(marker));
    let first = if is_ietf { ietf() } else { iers(text) };
    let rows = if first.is_empty() {
        if is_ietf { iers(text) } else { ietf() }
    } else {
        first
    };
    if rows.is_empty() {
        leap_fallback_view()
    } else {
        rows
    }
}

/// Actual ordered leap lookup for file views; offsets need not be finite.
fn file_tai_utc(mjd: f64, rows: &[(f64, f64)]) -> f64 {
    if mjd < 41317.0 {
        return crate::time::pre_1972_tai_utc(mjd);
    }
    if rows.is_empty() {
        return if mjd >= 57754.0 { 37.0 } else { 32.0 };
    }
    let mut offset = 0.0;
    for &(start, value) in rows {
        if mjd >= start {
            offset = value;
        } else {
            break;
        }
    }
    offset
}

/// Monthly file parser with the same date/layout/sanity filters, using the
/// permissive loaded leap view only at the C 04 conversion's consumed lookup.
/// The accepted rows remain finite checked monthly points, stably MJD sorted.
fn parse_observed_file_view(text: &str, leaps: &[(f64, f64)]) -> Vec<DeltaTDataPoint> {
    let mut rows = Vec::new();
    for line in text.lines().map(file_trim) {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts = line.split_whitespace().collect::<Vec<_>>();
        if parts.len() < 4 {
            continue;
        }
        let parsed = (|| {
            let year = parts[0].parse::<i32>().ok()?;
            let month = parts[1].parse::<i32>().ok()?;
            let day = parts[2].parse::<i32>().ok()?;
            if !(1900..=2100).contains(&year)
                || !(1..=12).contains(&month)
                || !(1..=31).contains(&day)
            {
                return None;
            }
            let mjd = crate::time::calendar_to_mjd(year, month, day);
            let col4 = file_float(parts[3])?;
            let delta = if parts.len() >= 7 && (col4 - mjd).abs() < 1.0 {
                file_tai_utc(mjd, leaps) + 32.184 - file_float(parts[6])?
            } else {
                col4
            };
            if !(-500.0..=1000.0).contains(&delta) {
                return None;
            }
            DeltaTDataPoint::new(mjd, year, month, day, delta).ok()
        })();
        if let Some(row) = parsed {
            rows.push(row);
        }
    }
    rows.sort_by(|a, b| a.mjd().total_cmp(&b.mjd()));
    rows
}

/// Stored acquisition inputs, without environment/home or filesystem reads.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct IersAcquisitionRequest {
    /// Explicit cache directory. SDK resolution makes relative paths absolute
    /// against the current directory at acquisition, without home expansion.
    pub cache_dir: Option<PathBuf>,
    /// Caller-resolved common data root; its `iers_cache` child is the default.
    pub data_dir: Option<PathBuf>,
    /// Explicit automatic-download selection; unset means disabled in the SDK.
    pub auto_download: Option<bool>,
    /// Explicit network policy; unset means mode-dependent SDK policy.
    pub network_policy: Option<NetworkPolicy>,
    /// Calculation mode for the automatic network-policy gate.
    pub mode: ProviderMode,
}

/// Fully resolved inputs for an actually consumed acquisition.
#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedIersAcquisition {
    /// Absolute selected cache directory; local acquisition may create it.
    pub directory: PathBuf,
    /// Whether missing files may enter the existing IERS transport owner.
    pub auto_download: bool,
    /// Effective allow/sealed gate or an ordinary network-policy refusal.
    ///
    /// A policy refusal belongs to the handled download boundary: local rows
    /// still load and may be retained. It is not a directory-resolution failure.
    /// Sealed requests and refused policies never open a socket.
    pub network: Result<EffectiveNetworkPolicy>,
}

/// Boundary callback that resolves process settings only on an IERS lookup.
///
/// It receives immutable programmatic selections and returns owned resolved
/// values. Ordinary errors cause model fallback without a retained failed
/// acquisition. It must not substitute frame-timescale rows or caller outcomes,
/// or recursively query the same owner while resolving; acquisition is serialized.
pub type IersBoundaryResolver = fn(&IersAcquisitionRequest) -> Result<ResolvedIersAcquisition>;

/// Resolve stored SDK selections without reading process configuration.
///
/// An explicit cache directory wins over the supplied data root's
/// `iers_cache` child. With neither selection, this returns
/// [`Error::ConfigurationMissing`]. Relative paths become lexical absolute
/// paths using the current directory; failure to obtain that directory returns
/// the same error variant. No home expansion, environment, TOML, filesystem
/// creation or row acquisition occurs here.
///
/// The returned owned settings use the stored auto-download choice (default
/// false) and effective network policy for the selected provider mode. Directory
/// creation and handled download failures belong to the consumed acquisition.
fn resolve_sdk(request: &IersAcquisitionRequest) -> Result<ResolvedIersAcquisition> {
    let directory = request
        .cache_dir
        .clone()
        .or_else(|| {
            request
                .data_dir
                .as_ref()
                .map(|root| root.join("iers_cache"))
        })
        .ok_or(Error::ConfigurationMissing {
            missing_config: "iers_cache_directory",
            suggestion: "Set the IERS cache directory or attach parsed IERS rows.",
        })?;
    let directory = std::path::absolute(&directory).map_err(|_| Error::ConfigurationMissing {
        missing_config: "current_directory",
        suggestion: "Supply an absolute IERS cache directory.",
    })?;
    Ok(ResolvedIersAcquisition {
        directory,
        auto_download: request.auto_download.unwrap_or(false),
        network: Ok(effective_network_policy(
            request.network_policy.unwrap_or(NetworkPolicy::Auto),
            request.mode,
        )),
    })
}

/// Acquire existing local rows, attempting only genuinely missing cache files.
///
/// Directory/path errors remain acquisition failures and are retried on a later
/// lookup. An unsuccessful ordinary download falls through to the local files:
/// completed files remain installed; later downloads in that attempt are not
/// entered. Existing files are never refreshed merely because they are old.
fn acquire(settings: ResolvedIersAcquisition) -> Result<IersDataSnapshot> {
    acquire_reload(settings, false, None)
}

/// Acquire an explicit reload without discarding the caller's prior owner.
/// With automatic acquisition enabled, force selects all three downloads with
/// their force flag; otherwise only absent files are selected. A handled
/// transport failure leaves local completed files available for ordered reads.
/// Directory failure is returned; only the caller chooses model fallback.
fn acquire_reload(
    settings: ResolvedIersAcquisition,
    force_download: bool,
    previous: Option<&IersDataSnapshot>,
) -> Result<IersDataSnapshot> {
    if !settings.directory.exists() {
        std::fs::create_dir_all(&settings.directory).map_err(|_| Error::AssetRejected {
            filename: settings.directory.to_string_lossy().into_owned(),
            reason: "cannot create IERS cache directory",
        })?;
    }
    if settings.auto_download {
        let missing: Vec<_> = crate::download::IERS_FILES
            .iter()
            .filter(|file| force_download || !settings.directory.join(file.name).exists())
            .cloned()
            .collect();
        if !missing.is_empty() {
            // Policy ValueError is caught by the file-download boundary before
            // local loading, so it retains the same successful local snapshot.
            let attempt = settings.network.and_then(|network| {
                crate::download::refresh_iers_files(
                    &settings.directory,
                    &missing,
                    force_download,
                    network,
                    &crate::net::TransportConfig::default(),
                )
            });
            if let Err(error) = attempt {
                download_failure(error)?;
            }
        }
    }
    Ok(IersDataSnapshot::reload_in(&settings.directory, previous))
}

/// Classify only the download boundary's known caught exception roles.
/// URL/HTTP/timeout/I/O publication, invalid policy/mode ValueErrors and sealed
/// RuntimeError correspond to a failed attempt. Unknown typed failures propagate
/// unchanged; a caller cannot silently install rows after an unrelated refusal.
fn download_failure(error: Error) -> Result<bool> {
    match error {
        Error::InvalidNetworkPolicy { .. }
        | Error::InvalidCalcMode { .. }
        | Error::NetworkSealed { .. }
        | Error::HttpInvalidUrl { .. }
        | Error::HttpTimeout { .. }
        | Error::HttpStatus { .. }
        | Error::HttpTransport { .. }
        | Error::HttpResponseTooLarge { .. }
        | Error::InvalidTimeout { .. }
        | Error::InvalidBodyLimit { .. }
        | Error::HashMismatch { .. }
        | Error::AssetRejected { .. } => Ok(false),
        error => Err(error),
    }
}

/// Record an ordinary lookup failure and return the model-fallback marker.
///
/// Diagnostic emission uses the existing configured debug sink. It retains
/// the calculation's fallback behavior; unrelated panics are not intercepted.
fn fallback(error: Error) -> Option<f64> {
    crate::logging::emit(
        crate::logging::LogLevel::Debug,
        &format!("IERS Delta T lookup failed ({error}); falling back to Skyfield."),
    );
    None
}

/// Retained Julian-Day bounds for finals and observed Delta-T tables.
/// The first tuple element is the finals interval, the second the monthly
/// interval. Each optional pair contains inclusive first/last row epochs;
/// an empty or unloaded table has no interval. These bounds do not extrapolate.
pub type IersLoadedRanges = (Option<(f64, f64)>, Option<(f64, f64)>);

type RowCache = Arc<Mutex<Option<Arc<IersDataSnapshot>>>>;

/// Retain one successful immutable snapshot under the acquisition lock.
/// Failed resolution/acquisition leaves the cache unloaded for a retry. The
/// returned Arc bounds query lifetime independently of later clear/reload.
fn ensure_rows(
    request: &IersAcquisitionRequest,
    boundary: Option<IersBoundaryResolver>,
    cache: &RowCache,
) -> Result<Arc<IersDataSnapshot>> {
    let mut loaded = cache.lock().unwrap_or_else(|poison| poison.into_inner());
    if loaded.is_none() {
        let settings = boundary.unwrap_or(resolve_sdk)(request)?;
        let rows = acquire(settings)?;
        *loaded = Some(Arc::new(rows));
    }
    Ok(loaded
        .as_ref()
        .expect("successful row acquisition installs a snapshot")
        .clone())
}

/// Lazy context-local provider retaining owned row snapshots through clones.
///
/// A successful acquisition retains even empty finals/monthly tables, with the
/// packaged nonempty leap fallback. Resolution/directory failure leaves the
/// cache unloaded for retry; existing files are not refreshed merely for age.
/// Default SDK construction supplies no cache directory and enables no automatic
/// download. Configure a path or attach rows before requesting local acquisition.
///
/// `Clone`, [`Self::with_request`] and [`Self::with_boundary`] share the same
/// cache, even before loading. The first successful acquisition through any
/// sibling installs rows for all siblings, using that caller's stored settings.
/// A settings copy therefore is not an independent data generation. Use
/// [`Self::cleared`], [`Self::with_rows`] or successful [`Self::loaded`] to obtain
/// a new row owner; preceding serving clones and returned row Arcs remain valid.
/// No process-global row store or implicit SDK home lookup exists.
///
/// ```
/// use libephemeris_core::iers::{IersContextSource, IersDataSnapshot};
/// use libephemeris_core::time::{Calendar, julday};
///
/// let rows = IersDataSnapshot::from_texts(None, None,
///     Some("2020 1 1 69.0\n2020 2 1 70.0\n"));
/// let source = IersContextSource::default().with_rows(rows);
/// let jd_ut1 = julday(2020, 1, 1, 0.0, Calendar::Gregorian)?;
/// assert_eq!(source.provider()(jd_ut1), Some(69.0)); // seconds; no file or network
/// let previous = source.row_snapshot()?;
/// let cleared = source.cleared(); // no acquisition
/// assert_eq!(cleared.loaded_counts().finals, 0);
/// assert_eq!(previous.observed_rows().len(), 2); // retained generation survives
/// # Ok::<(), libephemeris_core::errors::Error>(())
/// ```
#[derive(Clone)]
pub struct IersContextSource {
    request: IersAcquisitionRequest,
    boundary: Option<IersBoundaryResolver>,
    cache: RowCache,
    lookup: Arc<dyn Fn(f64) -> Option<f64> + Send + Sync>,
}

impl fmt::Debug for IersContextSource {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        out.debug_struct("IersContextSource")
            .field("request", &self.request)
            .field("boundary_configured", &self.boundary.is_some())
            .field(
                "loaded",
                &self.cache.lock().ok().is_some_and(|rows| rows.is_some()),
            )
            .finish()
    }
}

impl PartialEq for IersContextSource {
    /// Compare acquisition configuration and already acquired snapshot values.
    fn eq(&self, other: &Self) -> bool {
        let same_boundary = match (self.boundary, other.boundary) {
            (None, None) => true,
            (Some(a), Some(b)) => std::ptr::fn_addr_eq(a, b),
            _ => false,
        };
        if self.request != other.request || !same_boundary {
            return false;
        }
        if Arc::ptr_eq(&self.cache, &other.cache) {
            return true;
        }
        let own = self
            .cache
            .lock()
            .unwrap_or_else(|poison| poison.into_inner())
            .clone();
        let other = other
            .cache
            .lock()
            .unwrap_or_else(|poison| poison.into_inner())
            .clone();
        own == other
    }
}

impl Default for IersContextSource {
    fn default() -> Self {
        Self::build(
            IersAcquisitionRequest::default(),
            None,
            Arc::new(Mutex::new(None)),
        )
    }
}

impl IersContextSource {
    /// Build a settings-only owner and its lazy original-UT1 lookup callback.
    ///
    /// The callback captures an owned request and a shared row-cache owner.
    /// Its first lookup serializes resolver and local acquisition under the
    /// cache lock, installing only a successful immutable snapshot. Resolution
    /// or directory failures return no IERS value with a diagnostic, leave the
    /// cache empty and are retried by a later lookup. A completed load, including
    /// a nonempty packaged leap fallback without finals/monthly coverage, stays
    /// retained until an explicit clear creates a new cache owner.
    ///
    /// After acquisition, the lock is released before querying the owned rows
    /// at UT1 Julian days; results are optional Delta-T seconds. Ordinary query
    /// errors return no value while preserving those acquired rows. Unrelated
    /// panics propagate; a later call can recover the poisoned lock. A resolver
    /// must not recursively query this same owner during serialized acquisition.
    ///
    /// The provider reference is borrowed from the returned source, and cannot
    /// outlive it. Clones keep the shared row allocation alive independently of
    /// caller files or later configuration stores. Construction performs no
    /// path resolution, I/O, network access or model/frame-time conversion.
    fn build(
        request: IersAcquisitionRequest,
        boundary: Option<IersBoundaryResolver>,
        cache: RowCache,
    ) -> Self {
        let held = cache.clone();
        let selected = request.clone();
        let lookup = Arc::new(move |jd| {
            let rows = match ensure_rows(&selected, boundary, &held) {
                Ok(rows) => rows,
                Err(error) => return fallback(error),
            };
            match rows.delta_t_seconds(jd) {
                Ok(value) => value,
                Err(error) => fallback(error),
            }
        });
        Self {
            request,
            boundary,
            cache,
            lookup,
        }
    }

    /// Return new stored settings sharing this owner's acquisition cache.
    ///
    /// No resolution, I/O or invalidation occurs. Sharing includes an unloaded
    /// cache: the first successful load through either sibling supplies both.
    /// Use [`Self::cleared`] on the returned value for independent acquisition.
    pub fn with_request(&self, request: IersAcquisitionRequest) -> Self {
        Self::build(request, self.boundary, self.cache.clone())
    }

    /// Return a resolver selection sharing this owner's loaded or unloaded cache.
    /// No callback runs here; first-successful-load semantics match
    /// [`Self::with_request`]. A new resolver does not replace retained rows.
    pub fn with_boundary(&self, boundary: Option<IersBoundaryResolver>) -> Self {
        Self::build(self.request.clone(), boundary, self.cache.clone())
    }

    /// Move parsed rows into a new independently owned acquisition generation.
    /// The returned source is already loaded and performs no acquisition for
    /// ordinary lookups. Existing clones retain their prior cache and rows.
    pub fn with_rows(&self, rows: IersDataSnapshot) -> Self {
        Self::build(
            self.request.clone(),
            self.boundary,
            Arc::new(Mutex::new(Some(Arc::new(rows)))),
        )
    }

    /// Start a fresh unloaded cache for this context; do not invoke acquisition.
    pub fn cleared(&self) -> Self {
        Self::build(
            self.request.clone(),
            self.boundary,
            Arc::new(Mutex::new(None)),
        )
    }

    /// Counts of presently retained rows without resolving settings or loading.
    /// An unloaded owner reports three zero counts; a leap-only completed load
    /// reports its actual fallback entries, independently of epoch coverage.
    pub fn loaded_counts(&self) -> IersLoadedCounts {
        self.cache
            .lock()
            .unwrap_or_else(|poison| poison.into_inner())
            .as_ref()
            .map(|rows| rows.counts())
            .unwrap_or(IersLoadedCounts {
                finals: 0,
                leap_seconds: 0,
                delta_t: 0,
            })
    }

    /// Present permissive file counts, with monthly C 04 rows taken from the
    /// prepared leap view. No settings resolution, parsing or load occurs.
    /// Checked SDK counts remain available through [`Self::loaded_counts`].
    pub fn loaded_source_counts(&self) -> IersLoadedCounts {
        self.cache
            .lock()
            .unwrap_or_else(|poison| poison.into_inner())
            .as_ref()
            .map(|rows| IersLoadedCounts {
                finals: rows.source_finals.len(),
                leap_seconds: rows.source_leap_fields.len(),
                delta_t: rows.source_observed_rows().len(),
            })
            .unwrap_or(IersLoadedCounts {
                finals: 0,
                leap_seconds: 0,
                delta_t: 0,
            })
    }

    /// Retain a currently loaded generation without resolving or loading files.
    /// The cache lock is released before a language boundary sorts its fields.
    pub fn loaded_row_snapshot(&self) -> Option<Arc<IersDataSnapshot>> {
        self.cache.lock().unwrap_or_else(|p| p.into_inner()).clone()
    }

    /// Ensure the actual row generation and retain its leap-second table.
    /// Lazy directory, policy and file acquisition follows the same owning
    /// source as observed Delta-T, even when only leaps are consumed. The
    /// returned immutable small table survives clear/reload of the originating
    /// context. No date guard, model fallback or ephemeris image is acquired.
    pub fn leap_second_snapshot(&self) -> Result<Arc<[LeapSecondEntry]>> {
        let rows = ensure_rows(&self.request, self.boundary, &self.cache)?;
        Ok(Arc::from(rows.leap_seconds()))
    }

    /// Acquire and retain the complete immutable IERS row generation.
    ///
    /// This uses the same serialized lazy path, policy and file acquisition as
    /// the existing leap/Delta-T consumers. The returned Arc survives a later
    /// clear/reload and exposes no mutable table or borrowed cache lock. No
    /// epoch, ephemeris image, frame clock or Delta-T model is consumed here.
    pub fn row_snapshot(&self) -> Result<Arc<IersDataSnapshot>> {
        ensure_rows(&self.request, self.boundary, &self.cache)
    }

    /// Explicitly reparse local data, returning a new owner on success.
    /// Every call reloads, including a warm owner and force=false. The force
    /// argument controls re-download only while automatic download is enabled;
    /// it never forces network access past that setting/policy. Prior clones
    /// and in-flight row Arcs retain their preceding immutable snapshot. A
    /// resolution/directory error leaves the original owner unchanged.
    pub fn loaded(&self, force_download: bool) -> Result<Self> {
        let previous = self
            .cache
            .lock()
            .unwrap_or_else(|poison| poison.into_inner())
            .clone();
        let settings = self.boundary.unwrap_or(resolve_sdk)(&self.request)?;
        let rows = acquire_reload(settings, force_download, previous.as_deref())?;
        Ok(self.with_rows(rows))
    }

    /// Resolve and create the selected local cache directory, without transport.
    /// Path/configuration errors are returned before any force callback. This
    /// boundary is repeated for each source path getter; no row load occurs.
    pub fn prepare_directory(&self) -> Result<PathBuf> {
        let settings = self.boundary.unwrap_or(resolve_sdk)(&self.request)?;
        if !settings.directory.exists() {
            std::fs::create_dir_all(&settings.directory).map_err(|_| Error::AssetRejected {
                filename: settings.directory.to_string_lossy().into_owned(),
                reason: "cannot create IERS cache directory",
            })?;
        }
        Ok(settings.directory)
    }

    /// Reparse the three captured paths against this owner's current rows.
    /// No force/configuration callback or transport is entered here. Returning
    /// a new owner leaves already serving Arcs valid. Paths use finals, leap,
    /// monthly order, Latin-1 decoding and the unchanged existing row parsers.
    pub fn reloaded_paths(&self, paths: &[PathBuf; 3]) -> Self {
        let previous = self.cache.lock().unwrap_or_else(|p| p.into_inner()).clone();
        self.with_rows(IersDataSnapshot::reload_paths(paths, previous.as_deref()))
    }

    /// Download one selected file after the caller has consumed its force/age
    /// branch. No second freshness decision is made. Known policy, URL, HTTP,
    /// timeout and asset publication failures map to the handled false result;
    /// unrelated typed errors remain terminal. False aborts the remaining file
    /// attempts at the outer load boundary, which still reparses local paths.
    /// This bool projection does not expose a complete acquisition error graph.
    pub fn download_selected_in(
        directory: &Path,
        kind: super::IersFileKind,
        policy: Result<EffectiveNetworkPolicy>,
    ) -> Result<bool> {
        let attempt = policy.and_then(|policy| {
            crate::download::refresh_iers_files(
                directory,
                &[crate::download::IERS_FILES[kind.index()]],
                true,
                policy,
                &crate::net::TransportConfig::default(),
            )
        });
        match attempt {
            Ok(_) => Ok(true),
            Err(error) => download_failure(error),
        }
    }

    /// Consume exactly one prepared primary/backup endpoint, without retry.
    /// The caller prepared the file before its policy callback and retains the
    /// force/path decision. Known native download failures return false; unknown
    /// failures propagate unchanged. Policy callback exceptions are handled by
    /// the language boundary, never by this owned-file primitive.
    pub fn finish_selected_attempt(
        prepared: crate::download::PreparedUnpinnedDownload,
        kind: super::IersFileKind,
        backup: bool,
        policy: Result<EffectiveNetworkPolicy>,
    ) -> Result<bool> {
        let file = kind.file();
        let url = if backup {
            file.backup_url
        } else {
            file.primary_url
        };
        let purpose = format!("download {}", file.name);
        match policy.and_then(|policy| {
            prepared.finish(
                url,
                &purpose,
                policy,
                &crate::net::TransportConfig::default(),
            )
        }) {
            Ok(_) => Ok(true),
            Err(error) => download_failure(error),
        }
    }

    /// Retained loaded ranges without configuration, I/O or a lazy reload.
    /// The owned row Arc bounds both envelope reads; an unloaded owner has no
    /// ranges. The language facade separately performs its actual lazy load at
    /// the consumed range getter, after file metadata and pre-load counts.
    pub fn loaded_ranges(&self) -> IersLoadedRanges {
        let rows = self.cache.lock().unwrap_or_else(|p| p.into_inner()).clone();
        rows.map(|rows| (rows.finals_range(), rows.observed_range()))
            .unwrap_or((None, None))
    }

    /// Retained permissive finals and prepared monthly file envelopes.
    /// An unloaded generation returns no ranges without acquiring files.
    pub fn loaded_source_ranges(&self) -> IersLoadedRanges {
        let rows = self.cache.lock().unwrap_or_else(|p| p.into_inner()).clone();
        rows.map(|rows| {
            let table = rows.source_observed_rows();
            let range = table
                .first()
                .zip(table.last())
                .map(|(first, last)| (first.mjd() + 2400000.5, last.mjd() + 2400000.5));
            (rows.source_finals_range(), range)
        })
        .unwrap_or((None, None))
    }

    /// Counts and file ages only, without consuming lazy range acquisition.
    /// Directory resolution, counts and per-existing-file clocks retain their
    /// order. A language boundary can then invoke its genuine public loader,
    /// including consumed policy callbacks, before reading loaded ranges.
    pub fn file_info_with_clock(
        &self,
        clock: impl FnMut() -> Result<f64>,
    ) -> Result<super::IersCacheInfo> {
        let directory = self.prepare_directory()?;
        let counts = self.loaded_counts();
        let mut info = super::iers_cache_info_with_clock(&directory, clock)?;
        info.data_points = counts.finals;
        info.leap_second_entries = counts.leap_seconds;
        info.delta_t_entries = counts.delta_t;
        Ok(info)
    }

    /// File metadata with prepared public file counts sampled before ages.
    /// Checked SDK metadata remains in [`Self::file_info_with_clock`]; the
    /// original clock-acquisition order is retained.
    pub fn source_file_info_with_clock(
        &self,
        clock: impl FnMut() -> Result<f64>,
    ) -> Result<super::IersCacheInfo> {
        let directory = self.prepare_directory()?;
        let counts = self.loaded_source_counts();
        let mut info = super::iers_cache_info_with_clock(&directory, clock)?;
        info.data_points = counts.finals;
        info.leap_second_entries = counts.leap_seconds;
        info.delta_t_entries = counts.delta_t;
        Ok(info)
    }

    /// Snapshot metadata with a fresh Unix-second clock at each existing file.
    /// Counts precede ages and lazy ranges. Timestamp errors propagate through
    /// the typed filesystem boundary; prior row Arcs remain retained on refusal.
    pub fn cache_info_with_clock(
        &self,
        clock: impl FnMut() -> Result<f64>,
    ) -> Result<super::IersCacheInfo> {
        let directory = self.prepare_directory()?;
        let counts = self.loaded_counts();
        let mut info = super::iers_cache_info_with_clock(&directory, clock)?;
        info.data_points = counts.finals;
        info.leap_second_entries = counts.leap_seconds;
        info.delta_t_entries = counts.delta_t;
        let rows = ensure_rows(&self.request, self.boundary, &self.cache)?;
        info.data_range = rows.finals_range();
        info.delta_t_range = rows.observed_range();
        Ok(info)
    }

    /// Cache metadata with counts captured before the lazy range acquisition.
    /// The selected directory is resolved first. File presence and Unix ages
    /// use the existing explicit-time filesystem kernel; then range getters
    /// acquire rows only if cold. Thus a first call can report zero counts and
    /// nonempty ranges, matching the observable load order. A later call reports
    /// retained counts. Days/JD UT1 and Unix seconds stay separate, with no model
    /// TT or independent frame conversion. This typed SDK result does not carry
    /// the complete source filesystem exception graph.
    pub fn cache_info(&self, now_unix_seconds: f64) -> Result<super::IersCacheInfo> {
        let settings = self.boundary.unwrap_or(resolve_sdk)(&self.request)?;
        let counts = self.loaded_counts();
        let mut info = super::iers_cache_info_in(&settings.directory, now_unix_seconds);
        info.data_points = counts.finals;
        info.leap_second_entries = counts.leap_seconds;
        info.delta_t_entries = counts.delta_t;
        let rows = ensure_rows(&self.request, self.boundary, &self.cache)?;
        info.data_range = rows.finals_range();
        info.delta_t_range = rows.observed_range();
        Ok(info)
    }

    /// Remove only the three named files, then start an unloaded memory owner.
    /// Existing row Arcs remain readable. Missing files and individual removal
    /// errors do not count as deletions. Directory resolution/creation errors
    /// leave the memory owner unchanged and are returned. The caller controls
    /// which writable cache is selected; this operation is never implicit.
    pub fn deleted_files(&self) -> Result<(Self, usize)> {
        let settings = self.boundary.unwrap_or(resolve_sdk)(&self.request)?;
        if !settings.directory.exists() {
            std::fs::create_dir_all(&settings.directory).map_err(|_| Error::AssetRejected {
                filename: settings.directory.to_string_lossy().into_owned(),
                reason: "cannot create IERS cache directory",
            })?;
        }
        let count = super::delete_iers_cache_files_in(&settings.directory);
        Ok((self.cleared(), count))
    }

    /// One native value/rate query over the same acquired immutable row owner.
    /// Acquisition and scalar query errors retain the ordinary model-fallback
    /// seam; the derivative never retries acquisition or reads a later snapshot.
    pub(crate) fn delta_t_seconds_and_rate(&self, jd: f64) -> Option<(f64, f64)> {
        let rows = match ensure_rows(&self.request, self.boundary, &self.cache) {
            Ok(rows) => rows,
            Err(error) => {
                fallback(error);
                return None;
            }
        };
        let value = match rows.delta_t_seconds(jd) {
            Ok(Some(value)) => value,
            Ok(None) => return None,
            Err(error) => {
                fallback(error);
                return None;
            }
        };
        match super::delta_t_iers_rate_for_jd(jd, &rows.finals, &rows.observed, &rows.leaps) {
            Ok(Some(rate)) => Some((value, rate)),
            Ok(None) => None,
            Err(error) => {
                fallback(error);
                None
            }
        }
    }

    /// Borrow a synchronous `UT1 JD -> Option<TT−UT1 seconds>` callback.
    ///
    /// A cold call resolves settings and may create a directory, read files and
    /// attempt policy-controlled downloads. Successful acquisition, including
    /// empty Earth-orientation tables, stays shared by this cache's siblings.
    /// After acquisition, the row query runs without the acquisition lock.
    ///
    /// Uncovered epochs and ordinary acquisition/query errors return `None`
    /// for the outer Delta-T model fallback. Failed resolution/directory
    /// acquisition is retried later; query errors retain the loaded generation.
    /// The callback borrows this source and is not a promise of pure arithmetic.
    ///
    /// # Panics
    /// Boundary callback panics propagate. Resolution executes under the shared
    /// acquisition lock, so a resolver must not recursively query the same
    /// cache; doing so can deadlock. Poisoned locks are recovered on later calls.
    pub fn provider(&self) -> &(dyn Fn(f64) -> Option<f64> + Send + Sync) {
        &*self.lookup
    }
}

#[cfg(test)]
#[path = "iers_context_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "iers_lifecycle_tests.rs"]
mod lifecycle_tests;
