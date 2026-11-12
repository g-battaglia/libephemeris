//! IERS Earth-orientation rows, interpolation and local cache metadata.
//!
//! Daily finals rows carry UT1–UTC and monthly Delta-T rows carry TT–UT1, both
//! in seconds. Epochs are Modified Julian Days (`MJD = JD − 2400000.5`). Parsers
//! accept fixed-width finals records and whitespace Delta-T/C 04 layouts; malformed
//! lines are skipped. Query JDs use UT1 and are converted to MJD without an
//! additional UTC/UT1 correction. Finals lookup uses the query's floor day and
//! following day, including a single-sided value when only one is present.
//! With both rows it interpolates UT1–TAI across leap seconds. Monthly observed
//! Delta-T has priority within its closed first-to-last span;
//! otherwise finals data and the leap table give `TAI−UTC + 32.184 − (UT1−UTC)`.
//!
//! Pure row kernels perform no I/O. Cache helpers use an explicit directory;
//! [`IersContextSource`] owns lazy loading and policy-controlled acquisition.
//! Resetting a context does not invalidate snapshots retained by active calls.
//! Uncovered queries return `None`, and checked nonfinite queries return typed
//! errors. Prediction flags describe the supplied data rather than extrapolation
//! performed by these interpolation routines.

use crate::errors::{Error, Result};
use crate::time::{LeapSecondEntry, jd_to_mjd, mjd_to_jd, tai_utc_for_mjd};

#[path = "iers_context.rs"]
mod context_rows;
pub use context_rows::{
    IersAcquisitionRequest, IersBoundaryResolver, IersContextSource, IersDataSnapshot,
    IersLoadedRanges, ResolvedIersAcquisition,
};

/// Primary IERS finals URL (`finals2000A.data`).
///
/// URL used by the missing-file acquisition owner; row kernels never fetch.
pub const IERS_FINALS_URL: &str = "https://datacenter.iers.org/data/latestVersion/finals2000A.data";
/// Backup IERS finals URL.
pub const IERS_FINALS_URL_BACKUP: &str = "https://maia.usno.navy.mil/ser7/finals2000A.data";
/// Primary leap-second URL.
pub const IERS_LEAP_SECONDS_URL: &str = "https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat";
/// Backup leap-second URL.
pub const IERS_LEAP_SECONDS_URL_BACKUP: &str =
    "https://www.ietf.org/timezones/data/leap-seconds.list";
/// Primary observed Delta-T URL.
pub const IERS_DELTA_T_URL: &str = "https://maia.usno.navy.mil/ser7/deltat.data";
/// Backup observed Delta-T URL (C 04 one-file layout).
pub const IERS_DELTA_T_URL_BACKUP: &str =
    "https://datacenter.iers.org/data/latestVersion/EOP_14_C04_IAU1980_one_file_1962-now.txt";

/// Cache filenames within the IERS cache directory.
///
/// basenames (`finals2000A.data`,
/// `leap_seconds.dat`, `deltat.data`).
pub const FINALS_CACHE_FILENAME: &str = "finals2000A.data";
/// Cache filename for the leap-second asset. See [`FINALS_CACHE_FILENAME`].
pub const LEAP_SECONDS_CACHE_FILENAME: &str = "leap_seconds.dat";
/// Cache filename for the observed Delta-T asset. See [`FINALS_CACHE_FILENAME`].
pub const DELTA_T_CACHE_FILENAME: &str = "deltat.data";

/// Advisory refresh age in days for explicit download operations.
/// Lazy [`IersContextSource`] acquisition does not refresh an existing file
/// merely because its age exceeds this value.
pub const DEFAULT_MAX_AGE_DAYS: f64 = 30.0;

/// Stateless IERS cache-file management over an explicit directory.
///
/// [`delete_iers_cache_files_in`] removes only [`FINALS_CACHE_FILENAME`],
/// [`LEAP_SECONDS_CACHE_FILENAME`] and [`DELTA_T_CACHE_FILENAME`], counting
/// successful deletions and ignoring individual removal failures. It does
/// not clear an owned row generation; a context or binding lifecycle owner
/// performs that separate memory operation.
///
/// [`iers_cache_info_in`] reports file presence and one-decimal ages in days
/// against supplied Unix seconds. Its loaded counts are zero and ranges are
/// `None` because this stateless helper has no row-owner argument. Context,
/// Session and binding lifecycle routes overlay their actual retained counts
/// and lazy UT1-JD ranges. A caller-supplied directory is distinct from each
/// boundary's effective cache-path resolution.
///
/// These two stateless helpers attempt to create an absent directory. They
/// swallow creation/read failures into zero deletion counts or absent/unknown
/// readings; checked lifecycle acquisition helpers report their own errors.
use std::path::{Path, PathBuf};

/// Cache basenames in finals, leap-second, observed-Delta-T order.
pub fn iers_cache_filenames() -> [&'static str; 3] {
    [
        FINALS_CACHE_FILENAME,
        LEAP_SECONDS_CACHE_FILENAME,
        DELTA_T_CACHE_FILENAME,
    ]
}

/// Delete the three named IERS files inside `dir`, returning successful removals.
///
/// Missing paths, dangling links and individual removal failures are ignored.
/// An absent directory is created best-effort before deletion. This changes
/// filesystem state but does not clear any retained [`IersDataSnapshot`].
pub fn delete_iers_cache_files_in(dir: &Path) -> usize {
    let _ = std::fs::create_dir_all(dir);
    let mut deleted = 0usize;
    for name in iers_cache_filenames() {
        let path = dir.join(name);
        if path.exists() && std::fs::remove_file(&path).is_ok() {
            deleted += 1;
        }
    }
    deleted
}

/// Presence plus rounded age of one cache file.
///
/// Mirrors one `get_iers_cache_info` file triple: `exists` is the
/// file gate, and `age_days` is the model `round(age, 1)` readout in
/// days (`None` when the file is absent or its timestamp is unreadable).
/// Snapshot only: no invariant is enforced, so the fields are public.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IersCacheFileStatus {
    /// Whether the cache file exists.
    pub exists: bool,
    /// Rounded age in days, or `None` when absent/unreadable.
    pub age_days: Option<f64>,
}

/// Filesystem half of `get_iers_cache_info` dictionary.
///
/// File presence/age and an optional owned-row metadata projection.
/// The filesystem-only constructor fills zero counts/absent ranges. A context
/// owner overlays its pre-acquisition counts and lazy JD UT1 ranges. This
/// owned value contains no file handle, row borrow or process-global reference.
#[derive(Debug, Clone, PartialEq)]
pub struct IersCacheInfo {
    /// The directory this snapshot was taken over (as supplied).
    pub cache_dir: PathBuf,
    /// `finals2000A.data` presence + age.
    pub finals: IersCacheFileStatus,
    /// `leap_seconds.dat` presence + age.
    pub leap_seconds: IersCacheFileStatus,
    /// `deltat.data` presence + age.
    pub delta_t: IersCacheFileStatus,
    /// Loaded finals row count; zero in the filesystem-only constructor.
    pub data_points: usize,
    /// Loaded leap row count; zero in the filesystem-only constructor.
    pub leap_second_entries: usize,
    /// Loaded monthly row count; zero in the filesystem-only constructor.
    pub delta_t_entries: usize,
    /// Optional inclusive finals envelope in JD UT1.
    pub data_range: Option<(f64, f64)>,
    /// Optional inclusive monthly envelope in JD UT1.
    pub delta_t_range: Option<(f64, f64)>,
}

/// Round the exact binary64 age to one decimal place, ties to even.
/// A finite value is m*2^e; multiplying its integer significand by five and
/// shifting by e+1 evaluates the decimal-grid quotient without binary scaling
/// error. The selected integer tenth is converted once by decimal parsing,
/// avoiding an integer-to-float/division double rounding. Signed zero and
/// nonfinite values are preserved; large already integral values need no change.
fn round_age_days_tenth(age_days: f64) -> f64 {
    let bits = age_days.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    if exponent == 0x7ff {
        return age_days;
    }
    let fraction = bits & ((1_u64 << 52) - 1);
    let (mantissa, shift) = if exponent == 0 {
        (fraction, -1073)
    } else {
        (fraction | (1_u64 << 52), exponent - 1074)
    };
    if shift >= 0 {
        return age_days;
    }
    let scaled = mantissa * 5;
    let right = (-shift) as u32;
    let rounded = if right >= 64 {
        0
    } else {
        let integral = scaled >> right;
        let remainder = scaled & ((1_u64 << right) - 1);
        let half = 1_u64 << (right - 1);
        integral + u64::from(remainder > half || (remainder == half && integral & 1 != 0))
    };
    let decimal = format!("{}.{:01}", rounded / 10, rounded % 10);
    let value: f64 = decimal.parse().expect("finite integer-tenth decimal");
    value.copysign(age_days)
}

/// Signed Unix seconds, including valid modification times before 1970.
/// The absolute duration is converted once; no date/model timescale is involved.
pub fn iers_unix_seconds(time: std::time::SystemTime) -> f64 {
    match time.duration_since(std::time::UNIX_EPOCH) {
        Ok(duration) => duration.as_secs_f64(),
        Err(error) => -error.duration().as_secs_f64(),
    }
}

/// Strict age read after the caller's fresh clock sample, in days.
/// A missing/unreadable timestamp is a typed filesystem refusal. The time
/// difference and division retain their order; future files have negative age.
pub fn iers_file_age_at(path: &Path, now: f64) -> crate::errors::Result<f64> {
    let modified = std::fs::metadata(path)
        .and_then(|m| m.modified())
        .map_err(|_| crate::errors::Error::AssetRejected {
            filename: path.to_string_lossy().into_owned(),
            reason: "cannot read IERS cache modification time",
        })?;
    Ok((now - iers_unix_seconds(modified)) / 86_400.0)
}

/// Source-ordered file metadata with caller-owned fresh Unix-second samples.
/// All three existence gates precede clock/timestamp reads. Only existing files
/// consume a clock, in finals/leap/monthly order; then round their exact age to
/// one decimal. Directory and timestamp errors propagate, without partial info.
/// Counts/ranges remain zero here and belong to the row-owning overlay.
pub fn iers_cache_info_with_clock(
    dir: &Path,
    mut clock: impl FnMut() -> crate::errors::Result<f64>,
) -> crate::errors::Result<IersCacheInfo> {
    if !dir.exists() {
        std::fs::create_dir_all(dir).map_err(|_| crate::errors::Error::AssetRejected {
            filename: dir.to_string_lossy().into_owned(),
            reason: "cannot create IERS cache directory",
        })?;
    }
    let paths = iers_cache_filenames().map(|name| dir.join(name));
    let exists = paths.each_ref().map(|p| p.exists());
    let mut status = [IersCacheFileStatus {
        exists: false,
        age_days: None,
    }; 3];
    for i in 0..3 {
        status[i].exists = exists[i];
        if exists[i] {
            let now = clock()?;
            status[i].age_days = Some(round_age_days_tenth(iers_file_age_at(&paths[i], now)?));
        }
    }
    Ok(IersCacheInfo {
        cache_dir: dir.to_path_buf(),
        finals: status[0],
        leap_seconds: status[1],
        delta_t: status[2],
        data_points: 0,
        leap_second_entries: 0,
        delta_t_entries: 0,
        data_range: None,
        delta_t_range: None,
    })
}

/// Snapshot the IERS cache directory `dir` at Unix epoch `now_unix_secs`.
///
/// Mirrors the filesystem half of `get_iers_cache_info`
/// (`iers_data.py`): per-file existence plus the rounded age in days.
/// The supplied clock value avoids an internal wall-clock read, but this is
/// still a filesystem operation: it best-effort creates a missing directory
/// and reads file metadata. Creation/read failures become absent/unknown
/// readings. A nonfinite `now_unix_secs` gives `None` ages. No rows are loaded
/// and no network request is made.
pub fn iers_cache_info_in(dir: &Path, now_unix_secs: f64) -> IersCacheInfo {
    let _ = std::fs::create_dir_all(dir);
    let finite_now = now_unix_secs.is_finite();
    let status = |name: &str| {
        let path = dir.join(name);
        let metadata = std::fs::metadata(&path).ok();
        let age_days = metadata.as_ref().and_then(|md| {
            let mtime = md.modified().ok()?;
            let mtime_secs = iers_unix_seconds(mtime);
            if !finite_now {
                return None;
            }
            Some(round_age_days_tenth(
                (now_unix_secs - mtime_secs) / 86_400.0,
            ))
        });
        IersCacheFileStatus {
            exists: metadata.is_some(),
            age_days,
        }
    };
    // Destructure order follows the model field order.
    let filenames = iers_cache_filenames();
    IersCacheInfo {
        cache_dir: dir.to_path_buf(),
        finals: status(filenames[0]),
        leap_seconds: status(filenames[1]),
        delta_t: status(filenames[2]),
        data_points: 0,
        leap_second_entries: 0,
        delta_t_entries: 0,
        data_range: None,
        delta_t_range: None,
    }
}

/// MJD of 2000-01-01, the epoch of the model two-digit-year pivot:
/// finals rows on/after it decode to 2000+, earlier rows to 1900+.
const YEAR_PIVOT_MJD: f64 = 51544.0;

/// TT-TAI offset in seconds (`TT = TAI + 32.184 s`), applied exactly as in
/// `get_delta_t_iers` / C 04 derivation.
const TT_TAI_SECONDS: f64 = 32.184;

/// Plausibility band for ingested Delta-T rows in seconds.
///
/// Mirrors the model sanity band: observed Delta T is tens of seconds, so
/// anything outside `[-500, 1000]` (e.g. an MJD misread as Delta T) is
/// rejected and can never poison the table.
const DELTA_T_SANITY_MIN: f64 = -500.0;
/// Upper end of the Delta-T plausibility band. See `DELTA_T_SANITY_MIN`.
const DELTA_T_SANITY_MAX: f64 = 1000.0;

/// One IERS Earth-orientation row: daily UT1-UTC with its observed/
/// predicted provenance.
///
/// Fields are private so a consumer can
/// never observe a non-finite MJD/UT1-UTC (see [`IersDataPoint::new`]); the
/// year/month/day labels are carried from the asset unchanged, exactly as the
/// source stores them without validation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IersDataPoint {
    mjd: f64,
    year: i32,
    month: i32,
    day: i32,
    ut1_utc: f64,
    ut1_utc_err: Option<f64>,
    is_prediction: bool,
}

impl IersDataPoint {
    /// Build a row, rejecting non-finite MJD, UT1-UTC, or error values.
    ///
    /// The source's `float` would admit `nan`/`inf` literals; those rows
    /// are skipped by the parser instead (documented strictness).
    pub fn new(
        mjd: f64,
        year: i32,
        month: i32,
        day: i32,
        ut1_utc: f64,
        ut1_utc_err: Option<f64>,
        is_prediction: bool,
    ) -> Result<Self> {
        if !mjd.is_finite() {
            return Err(Error::NonFiniteArgument { value: mjd });
        }
        if !ut1_utc.is_finite() {
            return Err(Error::NonFiniteArgument { value: ut1_utc });
        }
        if let Some(err) = ut1_utc_err {
            if !err.is_finite() {
                return Err(Error::NonFiniteArgument { value: err });
            }
        }
        Ok(Self {
            mjd,
            year,
            month,
            day,
            ut1_utc,
            ut1_utc_err,
            is_prediction,
        })
    }

    /// Modified Julian Date of this row.
    pub fn mjd(self) -> f64 {
        self.mjd
    }

    /// Calendar year label (four digits via the MJD-epoch pivot).
    pub fn year(self) -> i32 {
        self.year
    }

    /// Calendar month label.
    pub fn month(self) -> i32 {
        self.month
    }

    /// Calendar day label.
    pub fn day(self) -> i32 {
        self.day
    }

    /// UT1-UTC in seconds.
    pub fn ut1_utc(self) -> f64 {
        self.ut1_utc
    }

    /// UT1-UTC error in seconds, when the asset carries it.
    pub fn ut1_utc_err(self) -> Option<f64> {
        self.ut1_utc_err
    }

    /// `true` when `P` (prediction) flag marks this row;
    /// `false` for observed (`I`) rows.
    pub fn is_prediction(self) -> bool {
        self.is_prediction
    }
}

/// One observed Delta-T row: monthly TT-minus-UT1 in seconds.
///
/// Fields are private so a consumer can
/// never observe a non-finite MJD/Delta-T (see [`DeltaTDataPoint::new`]);
/// the year/month/day labels are carried from the asset unchanged.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DeltaTDataPoint {
    mjd: f64,
    year: i32,
    month: i32,
    day: i32,
    delta_t: f64,
}

impl DeltaTDataPoint {
    /// Build a row, rejecting non-finite MJD or Delta-T values.
    pub fn new(mjd: f64, year: i32, month: i32, day: i32, delta_t: f64) -> Result<Self> {
        if !mjd.is_finite() {
            return Err(Error::NonFiniteArgument { value: mjd });
        }
        if !delta_t.is_finite() {
            return Err(Error::NonFiniteArgument { value: delta_t });
        }
        Ok(Self {
            mjd,
            year,
            month,
            day,
            delta_t,
        })
    }

    /// Modified Julian Date of this row.
    pub fn mjd(self) -> f64 {
        self.mjd
    }

    /// Calendar year label.
    pub fn year(self) -> i32 {
        self.year
    }

    /// Calendar month label.
    pub fn month(self) -> i32 {
        self.month
    }

    /// Calendar day label.
    pub fn day(self) -> i32 {
        self.day
    }

    /// Delta T (TT minus UT1) in seconds.
    pub fn delta_t(self) -> f64 {
        self.delta_t
    }
}

/// Parse `finals2000A.data` content into daily UT1-UTC rows.
///
/// Fixed-width columns, exactly as in the:
/// year (last 2 digits, bytes 0-2), month (2-4), day (4-6), MJD (7-15),
/// UT1-UTC flag `I`/`P` (byte 57, i.e. column 58 one-indexed), UT1-UTC
/// seconds (58-68), error seconds (68-78, when the line is that long).
/// Lines shorter than 68 bytes, blank fields, and malformed numbers are
/// skipped, exactly like the source's `try/except: continue`.
///
/// The two-digit year decodes against the MJD epoch (`MJD 51544.0` =
/// 2000-01-01): rows on/after it are 2000+, earlier rows 1900+ — never a
/// fixed two-digit pivot. The prediction flag compares ASCII
/// case-insensitively (`P` == predicted, anything else observed). A repeated
/// MJD keeps the last occurrence (the source dict overwrites); the result is
/// sorted ascending by MJD (stable).
///
/// Polar-motion columns (including the column-17 flag named in the model
/// docstring) are not consumed.
pub fn parse_finals_2000a(text: &str) -> Vec<IersDataPoint> {
    let mut rows: Vec<IersDataPoint> = Vec::new();
    for raw_line in text.lines() {
        // The asset is ASCII in practice; a line that is not valid UTF-8
        // fails conversion and is skipped (same class as malformed).
        let bytes = raw_line.as_bytes();
        if bytes.len() < 68 {
            continue;
        }
        let field =
            |lo: usize, hi: usize| -> Option<&str> { std::str::from_utf8(bytes.get(lo..hi)?).ok() };
        let row = (|| -> Option<IersDataPoint> {
            let year_2digit: i32 = field(0, 2)?.trim().parse().ok()?;
            let month: i32 = field(2, 4)?.trim().parse().ok()?;
            let day: i32 = field(4, 6)?.trim().parse().ok()?;
            let mjd_text = field(7, 15)?.trim();
            if mjd_text.is_empty() {
                return None;
            }
            let mjd: f64 = mjd_text.parse().ok()?;
            let year = if mjd >= YEAR_PIVOT_MJD {
                2000 + year_2digit
            } else {
                1900 + year_2digit
            };
            let ut1_text = field(58, 68)?.trim();
            if ut1_text.is_empty() {
                return None;
            }
            let ut1_utc: f64 = ut1_text.parse().ok()?;
            // Byte 57 always exists under the 68-byte guard (the source's
            // `len(line) > 57` check is statically true for the same reason).
            let is_prediction = bytes[57].eq_ignore_ascii_case(&b'P');
            let ut1_utc_err = if bytes.len() >= 78 {
                let err_text = field(68, 78)?.trim();
                if err_text.is_empty() {
                    None
                } else {
                    Some(err_text.parse::<f64>().ok()?)
                }
            } else {
                None
            };
            IersDataPoint::new(mjd, year, month, day, ut1_utc, ut1_utc_err, is_prediction).ok()
        })();
        if let Some(point) = row {
            rows.push(point);
        }
    }
    // Stable ascending sort, then keep the last occurrence of a repeated
    // MJD (source-dict overwrite semantics): reverse, dedup (keeps the
    // first == last occurrence), reverse back.
    rows.sort_by(|a, b| a.mjd.total_cmp(&b.mjd));
    rows.reverse();
    let mut seen = std::collections::HashSet::new();
    rows.retain(|row| seen.insert(row.mjd.to_bits()));
    rows.reverse();
    rows
}

/// Parse observed Delta-T content into monthly rows.
///
/// Two on-disk layouts, distinguished exactly as in the:
/// - Primary `deltat.data`: `YEAR MONTH DAY DELTA_T` with Delta T
///   (TT minus UT1) already in seconds.
/// - Backup C 04 one-file: `year month day MJD x y UT1-UTC ...`, where the
///   4 th column is the MJD (matching this row's date within a day) and
///   Delta T derives as `TAI-UTC + 32.184 - (UT1-UTC)` against `leap_table`
///
/// Year/month/day ranges are validated (`1900..=2100`, `1..=12`, `1..=31`)
/// and rows outside the `DELTA_T_SANITY_MIN..=DELTA_T_SANITY_MAX`
/// plausibility band are rejected, so a mis-detected column (e.g. an MJD
/// read as Delta T) can never poison the table. Comment (`#`), blank, short
/// (`< 4` columns), and malformed lines are skipped. The result is sorted
/// ascending by MJD (stable).
pub fn parse_delta_t(text: &str, leap_table: &[LeapSecondEntry]) -> Vec<DeltaTDataPoint> {
    let mut rows: Vec<DeltaTDataPoint> = Vec::new();
    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let parsed = (|| -> Option<DeltaTDataPoint> {
            let year: i32 = parts[0].parse().ok()?;
            let month: i32 = parts[1].parse().ok()?;
            let day: i32 = parts[2].parse().ok()?;
            if !(1900..=2100).contains(&year) {
                return None;
            }
            if !(1..=12).contains(&month) {
                return None;
            }
            if !(1..=31).contains(&day) {
                return None;
            }
            // `calendar_to_mjd` on small ints is finite by construction.
            let mjd = crate::time::calendar_to_mjd(year, month, day);
            if !mjd.is_finite() {
                return None;
            }
            let col4: f64 = parts[3].parse().ok()?;
            let is_c04 = parts.len() >= 7 && (col4 - mjd).abs() < 1.0;
            let delta_t = if is_c04 {
                let ut1_utc: f64 = parts[6].parse().ok()?;
                // `mjd` is finite, so only an (impossible here) table error
                // can fail; such a row is skipped rather than invented.
                let tai_utc = tai_utc_for_mjd(mjd, leap_table).ok()?;
                tai_utc + TT_TAI_SECONDS - ut1_utc
            } else {
                col4
            };
            // NaN fails the band comparison and is rejected with it (the
            // source's chained comparison rejects NaN the same way).
            if !(DELTA_T_SANITY_MIN..=DELTA_T_SANITY_MAX).contains(&delta_t) {
                return None;
            }
            DeltaTDataPoint::new(mjd, year, month, day, delta_t).ok()
        })();
        if let Some(point) = parsed {
            rows.push(point);
        }
    }
    // All MJDs are finite by construction; `total_cmp` keeps the sort
    // stable for repeated MJDs.
    rows.sort_by(|a, b| a.mjd.total_cmp(&b.mjd));
    rows
}

/// UT1-UTC in seconds for a Modified Julian Date against `rows`.
///
/// The query coordinate is UT1 MJD. Searches for rows exactly at `floor(mjd)`
/// and `floor(mjd) + 1`, without requiring the input slice to be sorted. The
/// last matching row wins. An exact floor-day hit returns its value; with two
/// rows, interpolate UT1−TAI and restore TAI−UTC at the query, avoiding a
/// smeared leap-second step. `leap_table` must follow the ascending-order
/// contract of [`tai_utc_for_mjd`].
///
/// A single adjacent row returns its value even outside the first-to-last row
/// envelope; no two adjacent rows means `Ok(None)`. Thus this is not general
/// interpolation across arbitrary gaps. Nonfinite MJD returns
/// [`Error::NonFiniteArgument`]; leap-table lookup errors propagate.
pub fn ut1_utc_for_mjd(
    mjd: f64,
    rows: &[IersDataPoint],
    leap_table: &[LeapSecondEntry],
) -> Result<Option<f64>> {
    if !mjd.is_finite() {
        return Err(Error::NonFiniteArgument { value: mjd });
    }
    let floor = mjd.floor();
    let ceil = floor + 1.0;
    let frac = mjd - floor;
    let mut floor_hit: Option<f64> = None;
    let mut ceil_hit: Option<f64> = None;
    for row in rows {
        if row.mjd == floor {
            floor_hit = Some(row.ut1_utc);
        }
        if row.mjd == ceil {
            ceil_hit = Some(row.ut1_utc);
        }
    }
    if frac == 0.0 {
        if let Some(value) = floor_hit {
            return Ok(Some(value));
        }
    }
    match (floor_hit, ceil_hit) {
        (Some(lo), Some(hi)) => {
            let tai_floor = tai_utc_for_mjd(floor, leap_table)?;
            let tai_ceil = tai_utc_for_mjd(ceil, leap_table)?;
            let tai_now = tai_utc_for_mjd(mjd, leap_table)?;
            Ok(Some(
                (lo - tai_floor) + frac * ((hi - tai_ceil) - (lo - tai_floor)) + tai_now,
            ))
        }
        (Some(value), None) | (None, Some(value)) => Ok(Some(value)),
        (None, None) => Ok(None),
    }
}

/// Observed Delta T in seconds for a UT1 Julian Day against `rows`.
///
/// outside the first-to-last span
/// the result is `None`; inside, the bracketing monthly rows interpolate
/// linearly (binary search; `rows` must ascend by MJD, as the parsers
/// produce). A zero-width bracket returns its value directly.
pub fn observed_delta_t_for_jd(jd: f64, rows: &[DeltaTDataPoint]) -> Result<Option<f64>> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    if rows.is_empty() {
        return Ok(None);
    }
    let mjd = jd_to_mjd(jd)?;
    if mjd < rows[0].mjd || mjd > rows[rows.len() - 1].mjd {
        return Ok(None);
    }
    // A single row spans exactly its own instant.
    // Without this guard `right - 1` underflows the index type.
    if rows.len() == 1 {
        return Ok(Some(rows[0].delta_t));
    }
    let mut left = 0usize;
    let mut right = rows.len() - 1;
    while left < right - 1 {
        let mid = (left + right) / 2;
        if rows[mid].mjd <= mjd {
            left = mid;
        } else {
            right = mid;
        }
    }
    let p1 = rows[left];
    let p2 = rows[right];
    if p1.mjd == p2.mjd {
        return Ok(Some(p1.delta_t));
    }
    let frac = (mjd - p1.mjd) / (p2.mjd - p1.mjd);
    Ok(Some(p1.delta_t + frac * (p2.delta_t - p1.delta_t)))
}

/// Inclusive UT1-JD first-to-last envelope of ascending observed rows.
/// Returns `None` for empty input. No sorting, acquisition or interpolation is
/// performed; input ordering is the same precondition as the query kernel.
pub fn observed_delta_t_range_jd(rows: &[DeltaTDataPoint]) -> Option<(f64, f64)> {
    if rows.is_empty() {
        return None;
    }
    // Parser output ascends, so the ends bound the span; `mjd_to_jd` on
    // finite MJDs is infallible in practice — a non-finite row (impossible
    // via the checked constructor) surfaces as `None`, never a wrong span.
    let start = mjd_to_jd(rows[0].mjd).ok()?;
    let end = mjd_to_jd(rows[rows.len() - 1].mjd).ok()?;
    Some((start, end))
}

/// Whether observed Delta-T rows cover a UT1 Julian Day.
///
/// (inclusive ends).
pub fn is_observed_delta_t_available(jd: f64, rows: &[DeltaTDataPoint]) -> Result<bool> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    if rows.is_empty() {
        return Ok(false);
    }
    let mjd = jd_to_mjd(jd)?;
    Ok(rows[0].mjd <= mjd && mjd <= rows[rows.len() - 1].mjd)
}

/// TT−UT1 seconds at UT1 JD `jd`, using observed rows before finals.
///
/// The ascending monthly table wins within its closed span. Otherwise the
/// finals adjacent-day lookup supplies `TAI-UTC + 32.184 - (UT1-UTC)` under
/// [`ut1_utc_for_mjd`]'s single-sided edge policy. `Ok(None)` means neither
/// query can answer; this function does not select a predictive Delta-T model.
/// Nonfinite requests and consumed row/table errors propagate.
pub fn delta_t_iers_for_jd(
    jd: f64,
    finals_rows: &[IersDataPoint],
    delta_t_rows: &[DeltaTDataPoint],
    leap_table: &[LeapSecondEntry],
) -> Result<Option<f64>> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    if let Some(observed) = observed_delta_t_for_jd(jd, delta_t_rows)? {
        return Ok(Some(observed));
    }
    let mjd = jd_to_mjd(jd)?;
    let Some(ut1_utc) = ut1_utc_for_mjd(mjd, finals_rows, leap_table)? else {
        return Ok(None);
    };
    Ok(Some(
        tai_utc_for_mjd(mjd, leap_table)? + TT_TAI_SECONDS - ut1_utc,
    ))
}

/// Derivative of the selected observed/finals Delta-T piece, seconds/day.
/// Monthly brackets match the value binary search, including duplicate rows.
/// Finals interpolate UT1-TAI, so leap offsets cancel from Delta T. A lone
/// covered daily row is locally constant. None retains missing coverage.
pub(crate) fn delta_t_iers_rate_for_jd(
    jd: f64,
    finals: &[IersDataPoint],
    observed: &[DeltaTDataPoint],
    leaps: &[LeapSecondEntry],
) -> Result<Option<f64>> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    let mjd = jd_to_mjd(jd)?;
    if observed_delta_t_for_jd(jd, observed)?.is_some() {
        if observed.len() < 2 {
            return Ok(Some(0.0));
        }
        let mut left = 0;
        let mut right = observed.len() - 1;
        while left < right - 1 {
            let mid = (left + right) / 2;
            if observed[mid].mjd <= mjd {
                left = mid
            } else {
                right = mid
            }
        }
        let (a, b) = (observed[left], observed[right]);
        return Ok(Some(if a.mjd == b.mjd {
            0.0
        } else {
            (b.delta_t - a.delta_t) / (b.mjd - a.mjd)
        }));
    }
    let floor = mjd.floor();
    let ceil = floor + 1.0;
    let a = finals.iter().rfind(|r| r.mjd == floor);
    let b = finals.iter().rfind(|r| r.mjd == ceil);
    match (a, b) {
        (Some(a), Some(b)) => Ok(Some(
            -((b.ut1_utc - tai_utc_for_mjd(ceil, leaps)?)
                - (a.ut1_utc - tai_utc_for_mjd(floor, leaps)?)),
        )),
        (Some(_), None) | (None, Some(_)) => Ok(Some(0.0)),
        (None, None) => Ok(None),
    }
}

/// Adapt caller-owned IERS row slices to the
/// [`crate::time::DeltaTConfig::iers_lookup_seconds`] provider seam (no I/O, no globals, no network):
/// observed monthly rows win where they cover the epoch, otherwise the
/// finals-derived value, `None` where neither covers it (the chain then
/// falls through to its model arm — never a silent substitution). The
/// leap table is an explicit caller input, like every other query in
/// this module. Epochs are UT1 Julian days; values are TT−UT1 seconds.
/// Unlike [`delta_t_iers_for_jd`], this convenience callback maps typed query
/// errors, including nonfinite requests, to `None`. All row slices remain
/// borrowed for the callback's lifetime; it performs no lazy acquisition.
pub fn iers_delta_t_provider<'a>(
    finals_rows: &'a [IersDataPoint],
    delta_t_rows: &'a [DeltaTDataPoint],
    leap_table: &'a [LeapSecondEntry],
) -> impl Fn(f64) -> Option<f64> + 'a {
    move |tjdut| {
        delta_t_iers_for_jd(tjdut, finals_rows, delta_t_rows, leap_table)
            .ok()
            .flatten()
    }
}

/// UT1-JD minimum/maximum row envelope, or `None` for empty input.
/// This is metadata, not a guarantee of continuous coverage through gaps.
/// Actual availability follows [`ut1_utc_for_mjd`], whose single-sided policy
/// can also answer adjacent-day queries outside this envelope.
pub fn iers_data_range_jd(rows: &[IersDataPoint]) -> Option<(f64, f64)> {
    if rows.is_empty() {
        return None;
    }
    let mut min = rows[0].mjd;
    let mut max = rows[0].mjd;
    for row in &rows[1..] {
        if row.mjd < min {
            min = row.mjd;
        }
        if row.mjd > max {
            max = row.mjd;
        }
    }
    let start = mjd_to_jd(min).ok()?;
    let end = mjd_to_jd(max).ok()?;
    Some((start, end))
}

/// Whether finals rows can answer a UT1 Julian Day.
///
/// exactly, including its
/// floor/`floor+1` bracketing: `true` exactly when [`ut1_utc_for_mjd`]
/// would return a non-`None` value.
pub fn is_iers_data_available(jd: f64, rows: &[IersDataPoint]) -> Result<bool> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    let mjd_floor = jd_to_mjd(jd)?.floor();
    Ok(rows
        .iter()
        .any(|row| row.mjd == mjd_floor || row.mjd == mjd_floor + 1.0))
}

/// Exact-MJD finals row lookup.
///
/// Returns the row stored for `mjd` itself (repeated MJDs keep the last
/// parsed occurrence, so at most one row matches). This is how callers read
/// the observed/predicted provenance ([`IersDataPoint::is_prediction`]) of
/// the day-boundary values the interpolation blends.
pub fn find_finals_row(mjd: f64, rows: &[IersDataPoint]) -> Option<IersDataPoint> {
    if !mjd.is_finite() {
        return None;
    }
    rows.iter().find(|row| row.mjd == mjd).copied()
}

/// Which IERS cache file a download/load leg addresses.
///
/// Order mirrors [`crate::download::IERS_FILES`] (finals, leap seconds,
/// observed Delta T) and `download_iers_finals` /
/// `download_leap_seconds` / `download_delta_t_data` triple.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IersFileKind {
    /// `finals2000A.data`.
    Finals,
    /// `leap_seconds.dat`.
    LeapSeconds,
    /// `deltat.data`.
    DeltaT,
}

impl IersFileKind {
    /// Position in [`crate::download::IERS_FILES`].
    pub fn index(self) -> usize {
        match self {
            IersFileKind::Finals => 0,
            IersFileKind::LeapSeconds => 1,
            IersFileKind::DeltaT => 2,
        }
    }

    /// The refreshable-file row this leg fetches.
    pub fn file(self) -> &'static crate::download::IersFile {
        &crate::download::IERS_FILES[self.index()]
    }

    /// Cache filename this leg reads when loading.
    pub fn filename(self) -> &'static str {
        match self {
            IersFileKind::Finals => FINALS_CACHE_FILENAME,
            IersFileKind::LeapSeconds => LEAP_SECONDS_CACHE_FILENAME,
            IersFileKind::DeltaT => DELTA_T_CACHE_FILENAME,
        }
    }
}

/// Counts of the retained finals, leap-second and observed Delta-T rows.
///
/// Counts describe memory, not current file presence or temporal coverage.
/// An owned reload can retain preceding finals when the new file is missing;
/// missing monthly data instead empties that table. Missing or unreadable leap
/// data uses the packaged fallback. Clearing an owner resets its counts while
/// previously retained generations remain valid.
///
/// [`load_iers_counts_in`] supplies a separate stateless, parse-only snapshot.
/// It performs no download or installation into an owned or default cache.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IersLoadedCounts {
    /// Retained finals rows, including preceding rows preserved by a reload.
    pub finals: usize,
    /// Retained leap-second entries, including the packaged fallback.
    pub leap_seconds: usize,
    /// Retained observed Delta-T rows; missing monthly data leaves this empty.
    pub delta_t: usize,
}

impl IersLoadedCounts {
    /// `load_iers_data` return: true when any table holds rows.
    /// The hardcoded leap fallback keeps this true whenever parsing
    /// itself works, exactly like the model.
    pub fn verdict(self) -> bool {
        self.finals > 0 || self.leap_seconds > 0 || self.delta_t > 0
    }
}

/// Parse the three IERS cache files inside `dir` and count the rows.
///
/// This stateless helper reads each present UTF-8 file without downloading or
/// installing rows into an owned/default cache. Missing or unreadable finals
/// and monthly files give zero counts; a missing/unreadable leap file supplies
/// the packaged leap table. It does not retain finals from an earlier call or
/// provide the lifecycle semantics of [`IersContextSource::loaded`].
pub fn load_iers_counts_in(dir: &Path) -> IersLoadedCounts {
    let read = |name: &str| std::fs::read_to_string(dir.join(name)).ok();
    let finals = read(FINALS_CACHE_FILENAME)
        .map(|text| parse_finals_2000a(&text).len())
        .unwrap_or(0);
    let leap_seconds = read(LEAP_SECONDS_CACHE_FILENAME)
        .map(|text| crate::time::parse_leap_seconds(&text).len())
        .unwrap_or_else(|| crate::time::hardcoded_leap_seconds().len());
    let leap_table = read(LEAP_SECONDS_CACHE_FILENAME)
        .map(|text| crate::time::parse_leap_seconds(&text))
        .unwrap_or_else(crate::time::hardcoded_leap_seconds);
    let delta_t = read(DELTA_T_CACHE_FILENAME)
        .map(|text| parse_delta_t(&text, &leap_table).len())
        .unwrap_or(0);
    IersLoadedCounts {
        finals,
        leap_seconds,
        delta_t,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time::hardcoded_leap_seconds;

    /// Build one synthetic fixed-width finals line with the model column
    /// layout: `YYMMDD MJD flag UT1-UTC error`.
    fn finals_line(
        year_2digit: u32,
        month: u32,
        day: u32,
        mjd: f64,
        flag: char,
        ut1_utc: f64,
        err: Option<f64>,
    ) -> String {
        let mut line = vec![b' '; 100];
        let head = format!("{year_2digit:02}{month:02}{day:02} {mjd:<8}");
        line[..head.len()].copy_from_slice(head.as_bytes());
        line[57] = flag as u8;
        let value = format!("{ut1_utc:>10.7}");
        line[58..58 + value.len()].copy_from_slice(value.as_bytes());
        if let Some(err) = err {
            let err_text = format!("{err:>10.7}");
            line[68..68 + err_text.len()].copy_from_slice(err_text.as_bytes());
        }
        String::from_utf8(line).expect("synthetic finals line must be UTF-8")
    }

    #[test]
    fn finals_numeric_row_parses_with_provenance() {
        let text = finals_line(0, 1, 1, 51544.0, 'I', 0.1234567, Some(0.0001234));
        let rows = parse_finals_2000a(&text);
        assert_eq!(rows.len(), 1);
        let row = rows[0];
        assert_eq!((row.year(), row.month(), row.day()), (2000, 1, 1));
        assert_eq!(row.mjd(), 51544.0);
        assert!((row.ut1_utc() - 0.1234567).abs() < 1e-9);
        assert!(!row.is_prediction());
        assert!(row.ut1_utc_err().expect("error must parse") - 0.0001234 < 1e-9);
    }

    #[test]
    fn finals_prediction_flag_and_pivot() {
        // MJD below the 2000-01-01 epoch decodes to the 1900 s; `P` (either
        // case) marks predictions.
        let text = [
            finals_line(73, 1, 2, 41683.0, 'P', -0.5, None),
            finals_line(25, 6, 15, 60845.0, 'p', 0.25, None),
        ]
        .join("\n");
        let rows = parse_finals_2000a(&text);
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].year(), 1973);
        assert!(rows[0].is_prediction());
        assert_eq!(rows[1].year(), 2025);
        assert!(rows[1].is_prediction());
    }

    #[test]
    fn finals_skips_short_malformed_and_bad_rows() {
        let text = [
            "too short".to_string(),
            finals_line(0, 1, 1, 51544.0, 'I', 0.1, None),
            "991332 XX".to_string() + &" ".repeat(70),
            finals_line(0, 1, 3, f64::NAN, 'I', 0.1, None),
            "# a comment line padded out".to_string() + &" ".repeat(60),
        ]
        .join("\n");
        let rows = parse_finals_2000a(&text);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].mjd(), 51544.0);
    }

    #[test]
    fn finals_repeated_mjd_keeps_last_and_sorts() {
        let text = [
            finals_line(0, 1, 2, 51545.0, 'I', 0.2, None),
            finals_line(0, 1, 1, 51544.0, 'I', 0.1, None),
            finals_line(0, 1, 1, 51544.0, 'P', 0.3, None),
        ]
        .join("\n");
        let rows = parse_finals_2000a(&text);
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].mjd(), 51544.0);
        assert!((rows[0].ut1_utc() - 0.3).abs() < 1e-12);
        assert!(rows[0].is_prediction());
    }

    #[test]
    fn ut1_interpolates_midday_on_non_leap_day() {
        let table = hardcoded_leap_seconds();
        // 2025-06-15/16: no leap boundary between them, so the TAI-UTC
        // terms cancel and midday is the plain mean.
        let (lo, hi) = (60840.0, 60841.0);
        let text = [
            finals_line(25, 6, 15, lo, 'I', 0.10, None),
            finals_line(25, 6, 16, hi, 'I', 0.12, None),
        ]
        .join("\n");
        let rows = parse_finals_2000a(&text);
        let mid = ut1_utc_for_mjd(lo + 0.5, &rows, &table)
            .expect("finite query must not error")
            .expect("bracketed query must hit");
        assert!((mid - 0.11).abs() < 1e-12, "midday must average, got {mid}");
        // Exact-boundary hit returns the row itself.
        assert_eq!(ut1_utc_for_mjd(lo, &rows, &table).unwrap(), Some(0.10));
    }

    #[test]
    fn ut1_crosses_leap_day_without_midday_dip() {
        let table = hardcoded_leap_seconds();
        // Rows bracketing the 2017-01-01 leap second (MJD 57754): UT1-TAI
        // is near-continuous, so UT1-UTC steps by ~+1 s with TAI-UTC (the
        // added leap second holds UTC back while UT1 runs on).
        // Interpolating UT1-UTC naively would dip ~0.5 s at midday; the
        // UT1-TAI path must not.
        let (lo, hi) = (57753.0, 57754.0);
        let text = [
            finals_line(16, 12, 31, lo, 'I', -0.4143, None),
            finals_line(17, 1, 1, hi, 'I', -0.4143 + 1.0 + 0.001, None),
        ]
        .join("\n");
        let rows = parse_finals_2000a(&text);
        let mid = ut1_utc_for_mjd(lo + 0.5, &rows, &table).unwrap().unwrap();
        // Recompute the model formula independently from the parsed rows.
        let tai_lo = tai_utc_for_mjd(lo, &table).unwrap();
        let tai_hi = tai_utc_for_mjd(hi, &table).unwrap();
        let tai_mid = tai_utc_for_mjd(lo + 0.5, &table).unwrap();
        assert_eq!((tai_lo, tai_hi), (36.0, 37.0));
        let expected = (rows[0].ut1_utc() - tai_lo)
            + 0.5 * ((rows[1].ut1_utc() - tai_hi) - (rows[0].ut1_utc() - tai_lo))
            + tai_mid;
        assert!((mid - expected).abs() < 1e-12);
        assert!(
            (mid - rows[0].ut1_utc()).abs() < 0.01,
            "midday must stay near the boundary values, got {mid}"
        );
    }

    #[test]
    fn ut1_missing_rows_are_none_with_edge_policy() {
        let table = hardcoded_leap_seconds();
        let rows = parse_finals_2000a(&finals_line(0, 1, 1, 51544.0, 'I', 0.1, None));
        // Same-day queries hit the single row (floor) or its ceiling edge.
        assert!(ut1_utc_for_mjd(51544.0, &rows, &table).unwrap().is_some());
        assert!(ut1_utc_for_mjd(51544.7, &rows, &table).unwrap().is_some());
        // Far outside, no bracket exists.
        assert_eq!(ut1_utc_for_mjd(52000.0, &rows, &table).unwrap(), None);
        assert_eq!(ut1_utc_for_mjd(51544.0, &[], &table).unwrap(), None);
        // Availability takes Julian Days (it converts via jd_to_mjd
        // internally); MJD 51544.5 is JD 2451545.0.
        assert!(is_iers_data_available(2451545.0, &rows).unwrap());
        assert!(!is_iers_data_available(2451546.5, &rows).unwrap());
        assert!(iers_data_range_jd(&rows).is_some());
        assert!(iers_data_range_jd(&[]).is_none());
    }

    #[test]
    fn delta_t_primary_layout_parses_and_interpolates() {
        let table = hardcoded_leap_seconds();
        let text = "1973  2  1  43.4724\n1973  3  1  43.5648\n";
        let rows = parse_delta_t(text, &table);
        assert_eq!(rows.len(), 2);
        let mjd_feb = crate::time::calendar_to_mjd(1973, 2, 1);
        assert_eq!(rows[0].mjd(), mjd_feb);
        assert!((rows[0].delta_t() - 43.4724).abs() < 1e-12);
        // Midpoint between monthly rows interpolates linearly.
        let mjd_mar = crate::time::calendar_to_mjd(1973, 3, 1);
        let jd_mid = mjd_to_jd((mjd_feb + mjd_mar) / 2.0).unwrap();
        let mid = observed_delta_t_for_jd(jd_mid, &rows).unwrap().unwrap();
        assert!((mid - (43.4724 + 43.5648) / 2.0).abs() < 1e-9);
        // Outside the span there is no extrapolation.
        let before = mjd_to_jd(mjd_feb - 1.0).unwrap();
        assert_eq!(observed_delta_t_for_jd(before, &rows).unwrap(), None);
        assert!(!is_observed_delta_t_available(before, &rows).unwrap());
        assert!(is_observed_delta_t_available(jd_mid, &rows).unwrap());
        assert!(observed_delta_t_range_jd(&rows).is_some());
        assert!(observed_delta_t_range_jd(&[]).is_none());
    }

    #[test]
    fn single_monthly_row_answers_its_own_instant_only() {
        // Regression: with one row `right - 1` underflowed the index type
        // and panicked (exposed by the model-agreement sweep over these
        // parsers). the model
        // skips its search loop with `left == right == 0` and returns the
        // row through the zero-width arm; off-instant queries miss the
        // degenerate span exactly like any other out-of-range instant.
        let table = hardcoded_leap_seconds();
        let rows = parse_delta_t("2021  1  1  69.4000\n", &table);
        assert_eq!(rows.len(), 1);
        let mjd = crate::time::calendar_to_mjd(2021, 1, 1);
        let jd = mjd_to_jd(mjd).unwrap();
        assert_eq!(observed_delta_t_for_jd(jd, &rows).unwrap(), Some(69.4));
        assert_eq!(
            observed_delta_t_for_jd(mjd_to_jd(mjd - 0.5).unwrap(), &rows).unwrap(),
            None
        );
        assert_eq!(
            observed_delta_t_for_jd(mjd_to_jd(mjd + 0.5).unwrap(), &rows).unwrap(),
            None
        );
        assert!(is_observed_delta_t_available(jd, &rows).unwrap());
    }

    #[test]
    fn iers_provider_threads_rows_through_the_delta_t_chain() {
        use crate::time::{DeltaTConfig, DeltaTModel, deltat_with_config};
        let table = hardcoded_leap_seconds();
        // Two finals rows (MJD 60000/60001) plus two monthly observed rows
        // spanning 60000.0-60030.0.
        let finals = parse_finals_2000a(
            &[
                finals_line(23, 2, 25, 60000.0, 'I', -0.1012345, Some(0.0000123)),
                finals_line(23, 2, 26, 60001.0, 'I', -0.1023456, Some(0.0000123)),
                finals_line(23, 4, 6, 60040.0, 'I', -0.11, Some(0.0000123)),
            ]
            .join("\n"),
        );
        assert_eq!(finals.len(), 3);
        let monthly = parse_delta_t("2023  2  25  69.3000\n2023  3  27  69.4000\n", &table);
        assert_eq!(monthly.len(), 2);
        let provider = iers_delta_t_provider(&finals, &monthly, &table);
        let mjd_mid = 60000.5;
        let jd_mid = mjd_to_jd(mjd_mid).unwrap();
        // Observed rows cover the epoch: linear interpolation answers.
        assert_eq!(provider(jd_mid), Some(69.3 + 0.5 / 30.0 * 0.1));
        let cfg = DeltaTConfig {
            model: DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: crate::constants::TIDAL_DEFAULT,
            iers_enabled: true,
            iers_lookup_seconds: Some(&provider),
            smh2016_lookup_seconds: None,
        };
        let observed = 69.3 + 0.5 / 30.0 * 0.1;
        let adjust =
            crate::time::tid_acc_adjustment_seconds(jd_mid, crate::constants::TIDAL_DEFAULT);
        assert_eq!(
            deltat_with_config(jd_mid, &cfg),
            Ok((observed + adjust) / 86400.0)
        );
        // Past the monthly span but on a finals row the finals-derived
        // value answers instead (single-sided table edge).
        let jd_later = mjd_to_jd(60040.5).unwrap();
        assert!(provider(jd_later).is_some());
        // Past every table the provider misses and the chain falls back
        // to the model arm.
        let jd_far = mjd_to_jd(61000.0).unwrap();
        assert_eq!(provider(jd_far), None);
        let model_cfg = DeltaTConfig {
            model: DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: crate::constants::TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        };
        assert_eq!(
            deltat_with_config(jd_far, &cfg),
            deltat_with_config(jd_far, &model_cfg)
        );
    }

    #[test]
    fn delta_t_c04_layout_derives_from_ut1() {
        let table = hardcoded_leap_seconds();
        // C 04 row: `year month day MJD x y UT1-UTC ...`; the 4 th column
        // matches the row MJD so the backup layout triggers.
        let mjd = crate::time::calendar_to_mjd(2025, 6, 15);
        let text = format!("2025 06 15 {mjd:.1} 0.05 0.09 0.10 0.0 0.0\n");
        let rows = parse_delta_t(&text, &table);
        assert_eq!(rows.len(), 1);
        let tai = tai_utc_for_mjd(mjd, &table).unwrap();
        assert!((rows[0].delta_t() - (tai + TT_TAI_SECONDS - 0.10)).abs() < 1e-9);
    }

    #[test]
    fn delta_t_rejects_implausible_and_out_of_range_rows() {
        let table = hardcoded_leap_seconds();
        let mjd = crate::time::calendar_to_mjd(2025, 6, 15);
        let text = [
            // MJD-scale value misread as Delta T: rejected by the band.
            format!("2025 06 15 {mjd:.1}"),
            // Year/month/day outside the validated ranges.
            "1899 06 15 50.0".to_string(),
            "2025 13 15 50.0".to_string(),
            "2025 06 32 50.0".to_string(),
            "# comment".to_string(),
            String::new(),
            "2025 06 15".to_string(),
        ]
        .join("\n");
        assert!(parse_delta_t(&text, &table).is_empty());
    }

    #[test]
    fn delta_t_iers_prefers_observed_then_computes() {
        let table = hardcoded_leap_seconds();
        // Observed monthly table covers 1973-02..03.
        let delta_rows = parse_delta_t("1973  2  1  43.4724\n1973  3  1  43.5648\n", &table);
        // Finals rows cover the same span with constant UT1-UTC.
        let mjd_feb = crate::time::calendar_to_mjd(1973, 2, 1);
        let mjd_mar = crate::time::calendar_to_mjd(1973, 3, 1);
        let finals = [
            finals_line(73, 2, 1, mjd_feb, 'I', 0.2, None),
            finals_line(73, 3, 1, mjd_mar, 'I', 0.2, None),
        ]
        .join("\n");
        let finals_rows = parse_finals_2000a(&finals);
        // Inside the observed span the monthly value wins.
        let jd_mid = mjd_to_jd((mjd_feb + mjd_mar) / 2.0).unwrap();
        let got = delta_t_iers_for_jd(jd_mid, &finals_rows, &delta_rows, &table)
            .unwrap()
            .unwrap();
        let want = observed_delta_t_for_jd(jd_mid, &delta_rows)
            .unwrap()
            .unwrap();
        assert!((got - want).abs() < 1e-12);
        // Past the observed span the finals-derived value applies. The
        // finals getter needs daily rows bracketing the query, so add the
        // two day-boundary rows around it (monthly rows alone bracket
        // nothing — the model returns None there too).
        let mjd_apr = crate::time::calendar_to_mjd(1973, 4, 1);
        let extra = finals_line(73, 4, 1, mjd_apr, 'I', 0.2, None);
        let daily_lo = finals_line(73, 3, 11, mjd_mar + 10.0, 'I', 0.2, None);
        let daily_hi = finals_line(73, 3, 12, mjd_mar + 11.0, 'I', 0.2, None);
        let finals_wide =
            parse_finals_2000a(&[finals.as_str(), &extra, &daily_lo, &daily_hi].join("\n"));
        let jd_apr = mjd_to_jd(mjd_mar + 10.0).unwrap();
        let derived = delta_t_iers_for_jd(jd_apr, &finals_wide, &delta_rows, &table)
            .unwrap()
            .unwrap();
        let tai = tai_utc_for_mjd(mjd_mar + 10.0, &table).unwrap();
        assert!((derived - (tai + TT_TAI_SECONDS - 0.2)).abs() < 1e-9);
        // With neither source covering the epoch the answer is missing.
        assert_eq!(delta_t_iers_for_jd(jd_apr, &[], &[], &table).unwrap(), None);
    }

    #[test]
    fn non_finite_epochs_are_typed_errors() {
        let table = hardcoded_leap_seconds();
        let finals = parse_finals_2000a(&finals_line(0, 1, 1, 51544.0, 'I', 0.1, None));
        let delta = parse_delta_t("1973  2  1  43.4724\n", &table);
        for bad in [f64::NAN, f64::INFINITY] {
            assert!(matches!(
                ut1_utc_for_mjd(bad, &finals, &table),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                observed_delta_t_for_jd(bad, &delta),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                delta_t_iers_for_jd(bad, &finals, &delta, &table),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                is_iers_data_available(bad, &finals),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                is_observed_delta_t_available(bad, &delta),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
        assert!(IersDataPoint::new(f64::NAN, 2000, 1, 1, 0.0, None, false).is_err());
        assert!(DeltaTDataPoint::new(51544.0, 2000, 1, 1, f64::INFINITY).is_err());
        assert!(find_finals_row(f64::NAN, &finals).is_none());
    }

    #[test]
    fn provenance_lookup_reads_the_flag() {
        let rows = parse_finals_2000a(
            &[
                finals_line(0, 1, 1, 51544.0, 'I', 0.1, None),
                finals_line(0, 1, 2, 51545.0, 'P', 0.2, None),
            ]
            .join("\n"),
        );
        assert!(!find_finals_row(51544.0, &rows).unwrap().is_prediction());
        assert!(find_finals_row(51545.0, &rows).unwrap().is_prediction());
        assert!(find_finals_row(51546.0, &rows).is_none());
    }

    /// Read-only spot check against the installed IERS assets in `iers_cache/`.
    ///
    /// Never writes: it resolves the directory from `LIBEPHEMERIS_DATA_DIR`
    /// (with an `iers_cache` child) or the default `~/.libephemeris`, reads
    /// at most the two text files, and passes vacuously when they are
    /// absent (extended-tier-style gap, not a failure).
    #[test]
    fn real_assets_parse_with_expected_shape() {
        let cache_dir = std::env::var("LIBEPHEMERIS_DATA_DIR")
            .ok()
            .map(std::path::PathBuf::from)
            .map(|base| base.join("iers_cache"))
            .filter(|dir| dir.is_dir())
            .or_else(|| {
                std::env::var("HOME")
                    .ok()
                    .map(std::path::PathBuf::from)
                    .map(|home| home.join(".libephemeris").join("iers_cache"))
                    .filter(|dir| dir.is_dir())
            });
        let Some(cache_dir) = cache_dir else {
            return;
        };
        let finals_path = cache_dir.join(FINALS_CACHE_FILENAME);
        let delta_path = cache_dir.join(DELTA_T_CACHE_FILENAME);
        if !finals_path.is_file() || !delta_path.is_file() {
            return;
        }
        let finals_text = std::fs::read_to_string(&finals_path).expect("finals asset must read");
        let delta_text = std::fs::read_to_string(&delta_path).expect("delta-t asset must read");
        let table = hardcoded_leap_seconds();

        let finals_rows = parse_finals_2000a(&finals_text);
        assert!(
            finals_rows.len() > 10_000,
            "finals asset must hold daily rows, got {}",
            finals_rows.len()
        );
        assert!(
            finals_rows
                .windows(2)
                .all(|pair| pair[0].mjd() <= pair[1].mjd())
        );
        assert!(finals_rows.iter().any(|row| row.is_prediction()));

        let delta_rows = parse_delta_t(&delta_text, &table);
        assert!(
            delta_rows.len() > 500,
            "delta-t asset must hold monthly rows, got {}",
            delta_rows.len()
        );
        let (start_jd, end_jd) = observed_delta_t_range_jd(&delta_rows).unwrap();
        // Monthly observed Delta T starts in 1973 (JD ~2441700).
        assert!(start_jd > 2_441_000.0 && start_jd < 2_442_000.0);
        assert!(end_jd > start_jd);

        // A year-2000 epoch inside both spans answers from both paths and
        // the two Delta-T routes agree to interpolation tolerance.
        let jd_2000 = crate::time::julday(2000, 6, 15, 12.0, crate::time::Calendar::Gregorian)
            .expect("fixed test epoch must convert");
        assert!(is_iers_data_available(jd_2000, &finals_rows).unwrap());
        let from_iers = delta_t_iers_for_jd(jd_2000, &finals_rows, &delta_rows, &table)
            .unwrap()
            .unwrap();
        let from_monthly = observed_delta_t_for_jd(jd_2000, &delta_rows)
            .unwrap()
            .unwrap();
        assert!(
            (from_iers - from_monthly).abs() < 0.05,
            "observed and finals-derived Delta T must agree, got {from_iers} vs {from_monthly}"
        );
        assert!((from_iers - 64.0).abs() < 5.0, "year-2000 Delta T is ~64 s");
    }

    /// Hermetic scratch directory for the cache-file tests below
    /// (machine-local `temp_dir`, per-process name, removed first so a
    /// stale tree can never leak values between runs).
    fn cache_scratch_dir(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "libephemeris-iers-cache-test-{name}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("scratch dir");
        dir
    }

    #[test]
    fn cache_filenames_match_release_order() {
        // Basenames tails, in
        // the model deletion/info order.
        assert_eq!(
            iers_cache_filenames(),
            ["finals2000A.data", "leap_seconds.dat", "deltat.data",]
        );
    }

    #[test]
    fn delete_counts_cache_files_and_ignores_foreign_entries() {
        let dir = cache_scratch_dir("delete");
        for name in iers_cache_filenames() {
            std::fs::write(dir.join(name), b"stale").expect("fixture file");
        }
        std::fs::write(dir.join("other.txt"), b"keep").expect("foreign file");
        // A directory masquerading as a cache file is skipped, never
        // counted.
        std::fs::remove_file(dir.join(LEAP_SECONDS_CACHE_FILENAME)).expect("remove leap file");
        std::fs::create_dir_all(dir.join(LEAP_SECONDS_CACHE_FILENAME)).expect("masquerade dir");

        // Two real files + one masquerading directory → 2 deletions.
        assert_eq!(delete_iers_cache_files_in(&dir), 2);
        assert!(dir.join("other.txt").is_file());
        assert!(dir.join(LEAP_SECONDS_CACHE_FILENAME).is_dir());
        // Idempotent: second run deletes nothing.
        assert_eq!(delete_iers_cache_files_in(&dir), 0);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn delete_creates_missing_dir_and_reads_zero() {
        let dir = std::env::temp_dir().join(format!(
            "libephemeris-iers-cache-test-missing-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        assert_eq!(delete_iers_cache_files_in(&dir), 0);
        // side effect.
        assert!(dir.is_dir());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn info_snapshot_reports_presence_and_zero_load_halves() {
        let dir = cache_scratch_dir("info");
        let now = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("test clock must be post-epoch")
            .as_secs_f64();
        std::fs::write(dir.join(FINALS_CACHE_FILENAME), b"rows").expect("finals file");
        let info = iers_cache_info_in(&dir, now);
        assert_eq!(info.cache_dir, dir);
        assert!(info.finals.exists);
        // Fresh file: the tiny clock skew rounds to exactly 0.0.
        assert_eq!(info.finals.age_days, Some(0.0));
        assert!(!info.leap_seconds.exists);
        assert_eq!(info.leap_seconds.age_days, None);
        assert!(!info.delta_t.exists);
        assert_eq!(info.delta_t.age_days, None);
        // No process-global store: loaded halves always read empty.
        assert_eq!(info.data_points, 0);
        assert_eq!(info.leap_second_entries, 0);
        assert_eq!(info.delta_t_entries, 0);
        assert_eq!(info.data_range, None);
        assert_eq!(info.delta_t_range, None);
        // Non-finite epochs guard the snapshot instead of propagating.
        let guarded = iers_cache_info_in(&dir, f64::NAN);
        assert!(guarded.finals.exists);
        assert_eq!(guarded.finals.age_days, None);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn round_age_days_tenth_is_ties_even() {
        // 2.25 is exact in binary: halfway between 2.2 and 2.3 rounds
        // to the even tenth, matching the model `round(age, 1)`.
        assert_eq!(round_age_days_tenth(2.25), 2.2);
        assert_eq!(round_age_days_tenth(2.75), 2.8);
        assert_eq!(round_age_days_tenth(0.0), 0.0);
        assert_eq!(round_age_days_tenth(3.24), 3.2);
        assert_eq!(round_age_days_tenth(3.26), 3.3);
    }

    #[test]
    fn file_kind_rows_mirror_release_names_and_urls() {
        // Kinds index the shared `IERS_FILES` triple in the model
        // download order (finals, leap seconds, Delta T).
        assert_eq!(IersFileKind::Finals.index(), 0);
        assert_eq!(IersFileKind::LeapSeconds.index(), 1);
        assert_eq!(IersFileKind::DeltaT.index(), 2);
        assert_eq!(IersFileKind::Finals.filename(), "finals2000A.data");
        assert_eq!(IersFileKind::LeapSeconds.filename(), "leap_seconds.dat");
        assert_eq!(IersFileKind::DeltaT.filename(), "deltat.data");
        assert_eq!(
            IersFileKind::Finals.file().name,
            IersFileKind::Finals.filename()
        );
        assert_eq!(
            IersFileKind::LeapSeconds.file().name,
            IersFileKind::LeapSeconds.filename()
        );
        assert_eq!(
            IersFileKind::DeltaT.file().name,
            IersFileKind::DeltaT.filename()
        );
        for kind in [
            IersFileKind::Finals,
            IersFileKind::LeapSeconds,
            IersFileKind::DeltaT,
        ] {
            assert!(kind.file().primary_url.starts_with("https://"));
            assert!(kind.file().backup_url.starts_with("https://"));
        }
    }

    #[test]
    fn load_counts_parse_present_files_and_verdict() {
        let dir = cache_scratch_dir("load");
        let finals = [
            finals_line(25, 1, 1, 60676.0, 'I', 0.05, None),
            finals_line(25, 1, 2, 60677.0, 'I', 0.06, None),
        ]
        .join("\n");
        std::fs::write(dir.join(FINALS_CACHE_FILENAME), finals).expect("finals fixture");
        // IETF-shaped leap fixture: NTP seconds + TAI-UTC pairs.
        let leap = "#$ 3948576000\n3948576000 33\n3950000000 34\n";
        std::fs::write(dir.join(LEAP_SECONDS_CACHE_FILENAME), leap).expect("leap fixture");
        let delta_t = "2025 01 01 69.2\n2025 02 01 69.3\n";
        std::fs::write(dir.join(DELTA_T_CACHE_FILENAME), delta_t).expect("delta-t fixture");
        let counts = load_iers_counts_in(&dir);
        assert_eq!(counts.finals, 2);
        assert_eq!(counts.leap_seconds, 2);
        assert_eq!(counts.delta_t, 2);
        assert!(counts.verdict());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn load_counts_empty_dir_reads_leap_fallback_only() {
        let dir = cache_scratch_dir("load-empty");
        let counts = load_iers_counts_in(&dir);
        assert_eq!(counts.finals, 0);
        assert_eq!(counts.delta_t, 0);
        // `load_iers_data` falls back to its hardcoded leap
        // table when the file is missing, so the leap half (and the
        // verdict) stays populated even with an empty directory.
        assert_eq!(
            counts.leap_seconds,
            crate::time::hardcoded_leap_seconds().len()
        );
        assert!(counts.verdict());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn load_counts_garbage_files_parse_empty_without_panicking() {
        let dir = cache_scratch_dir("load-garbage");
        std::fs::write(dir.join(FINALS_CACHE_FILENAME), b"\xff\xfe binary").expect("finals");
        std::fs::write(dir.join(DELTA_T_CACHE_FILENAME), b"not rows at all").expect("deltat");
        let counts = load_iers_counts_in(&dir);
        assert_eq!(counts.finals, 0);
        assert_eq!(counts.delta_t, 0);
        assert_eq!(
            counts.leap_seconds,
            crate::time::hardcoded_leap_seconds().len()
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
}

#[cfg(test)]
mod live_iers_slope_tests {
    use super::*;
    #[test]
    fn live_iers_finals_slope_cancels_leap_step_and_preserves_single_row() {
        let leaps = crate::time::hardcoded_leap_seconds();
        let day = crate::time::calendar_to_mjd(2016, 12, 31);
        let rows = [
            IersDataPoint::new(day, 2016, 12, 31, -0.4, None, false).unwrap(),
            IersDataPoint::new(day + 1.0, 2017, 1, 1, 0.5, None, false).unwrap(),
        ];
        for fraction in [0.0, 0.5, 0.999] {
            let rate = delta_t_iers_rate_for_jd(day + 2400000.5 + fraction, &rows, &[], &leaps)
                .unwrap()
                .unwrap();
            assert!((rate - 0.1).abs() < 1e-14);
        }
        assert_eq!(
            delta_t_iers_rate_for_jd(day + 2400000.75, &rows[..1], &[], &leaps).unwrap(),
            Some(0.0)
        );
        assert_eq!(
            delta_t_iers_rate_for_jd(day + 2399999.0, &rows, &[], &leaps).unwrap(),
            None
        );
    }
}
