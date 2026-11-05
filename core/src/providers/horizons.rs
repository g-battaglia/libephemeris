//! NASA JPL Horizons vector requests, parsing, caching and bounded batch work.
//!
//! [`HorizonsClient`] acquires uncorrected ICRF vectors in AU and AU/day.
//! Requests identify center and time scale explicitly; a barycentric request
//! uses center `@0`. The position pipeline owns observer translation,
//! light-time, apparent corrections and output-frame conversion. Supplied JDs
//! are transmitted in the selected [`HorizonsTimeType`] without local clock
//! conversion; in particular, a TDB request does not convert a TT input.
//!
//! The response parser validates the JSON envelope, Horizons error payload,
//! `$$SOE`/`$$EOE` delimiters and vector fields. API refusals are terminal;
//! retryable transport failures use bounded exponential backoff.
//!
//! The owned client cache keys epoch, command, center and time scale and
//! uses bounded LRU storage. Batch work limits parallel native requests and
//! retains per-request errors. Sealed network policy refuses before sockets.
//! Wire escaping and decimal epoch formatting follow the Horizons protocol;
//! no astronomical values are embedded as request expectations.

use std::collections::{HashMap, VecDeque};
use std::num::NonZeroUsize;
use std::sync::Mutex;
use std::time::Duration;

use crate::constants::{
    FICT_OFFSET, HARRINGTON, NEPTUNE_ADAMS, NEPTUNE_LEVERRIER, PLUTO_LOWELL, PLUTO_PICKERING,
    PROSERPINA, VULCAN, WALDEMATH, WHITE_MOON,
};
use crate::context::EffectiveNetworkPolicy;
use crate::errors::{Error, Result};
use crate::net::{HttpRequest, TransportConfig, fetch};
use crate::spk::HORIZONS_API_URL;

/// User-Agent for Horizons vector requests, overriding the transport default.
pub const HORIZONS_USER_AGENT: &str = "libephemeris/1.0";

/// Default Horizons query timeout in seconds.
///
/// Mirrors `HorizonsClient(timeout: int = 30)`; SPK generation (120 s) and
/// the ASSIST bootstrap (300 s) are separate callers with their own
/// budgets.
pub const HORIZONS_DEFAULT_TIMEOUT_SECS: f64 = 30.0;

/// Default Horizons query center.
///
/// `"@0"` is the solar-system barycenter. A different requested center produces
/// vectors relative to that center; the client does not translate them to SSB.
pub const HORIZONS_DEFAULT_CENTER: &str = "@0";

/// Retry budget after the initial exchange: at most three attempts per vector.
pub const HORIZONS_MAX_RETRIES: u32 = 2;

/// Default LRU vector-cache capacity in entries.
///
/// Mirrors `HorizonsClient(max_cache_size: int = 4096)`.
pub const HORIZONS_DEFAULT_MAX_CACHE_SIZE: usize = 4096;

/// Default batch worker bound.
///
/// Mirrors `HorizonsClient(max_workers: int = 8)`.
pub const HORIZONS_DEFAULT_MAX_WORKERS: usize = 8;

/// Backoff base in seconds.
///
/// Mirrors `time.sleep(0.5 * (2 ** attempt))` between attempts.
pub const HORIZONS_RETRY_BASE_DELAY_SECS: f64 = 0.5;

/// Delay before retry number `retry_index` (0-based), in seconds.
///
/// `0.5 * 2^min(retry_index, 30)`, in seconds. The normal retry budget uses
/// indices 0 and 1 (0.5 s and 1.0 s). This helper only computes a delay;
/// [`HorizonsClient::fetch_vectors`] performs the blocking sleep.
pub fn retry_delay_secs(retry_index: u32) -> f64 {
    HORIZONS_RETRY_BASE_DELAY_SECS * 2f64.powi(retry_index.min(30) as i32)
}

/// Whether a failed exchange may be retried.
///
/// Retries HTTP timeout/transport/status errors and malformed Horizons
/// payloads. API refusal payloads, sealed policy, invalid URL/configuration,
/// unknown bodies and nonfinite epochs are terminal. An exhausted-budget
/// [`Error::HorizonsRequestFailed`] is itself terminal.
pub fn is_retryable_horizons_error(err: &Error) -> bool {
    matches!(
        err,
        Error::HttpTimeout { .. }
            | Error::HttpTransport { .. }
            | Error::HttpStatus { .. }
            | Error::HorizonsMalformed { .. }
    )
}

/// Horizons `COMMAND` string for a library body ID.
///
/// Supports major planets/Earth, Chiron and the four main asteroids through
/// the explicit command identifiers below. Pholus is not mapped. Unknown
/// library IDs return [`Error::UnknownBody`]; this function performs no request.
pub fn horizons_command(body_id: i32) -> Result<&'static str> {
    match body_id {
        0 => Ok("10"),
        1 => Ok("301"),
        2 => Ok("199"),
        3 => Ok("299"),
        4 => Ok("499"),
        5 => Ok("599"),
        6 => Ok("699"),
        7 => Ok("799"),
        8 => Ok("899"),
        9 => Ok("999"),
        14 => Ok("399"),
        15 => Ok("2060;"),
        17 => Ok("1;"),
        18 => Ok("2;"),
        19 => Ok("3;"),
        20 => Ok("4;"),
        _ => Err(Error::UnknownBody { body_id }),
    }
}

/// Whether a body is served analytically with no HTTP fetch.
///
/// The `_ANALYTICAL_BODIES` set with the same members, `{10, 12}` (Mean Node, Mean Apogee).
pub fn is_analytical_body(body_id: i32) -> bool {
    matches!(body_id, 10 | 12)
}

/// Whether a body is derived from one Moon state vector.
///
/// The `_MOON_DERIVED_BODIES` set with the same members, `{11, 13, 21, 22}` (True Node, Osculating
/// Apogee, Interpolated Apogee/Perigee).
pub fn is_moon_derived_body(body_id: i32) -> bool {
    matches!(body_id, 11 | 13 | 21 | 22)
}

/// Whether a body ID is a local hypothetical model.
///
/// The `_HYPOTHETICAL_BODIES` set with the same members, `frozenset(range(FICT_OFFSET,
/// WALDEMATH + 1))`, i.e. `40..=58` at the release pin.
pub fn is_hypothetical_body(body_id: i32) -> bool {
    (FICT_OFFSET..=WALDEMATH).contains(&body_id)
}

/// Whether a hypothetical body is natively geocentric.
///
/// The `_GEOCENTRIC_HYPOTHETICAL_BODIES` set with the same members, `{WHITE_MOON, WALDEMATH}`.
pub fn is_geocentric_hypothetical_body(body_id: i32) -> bool {
    body_id == WHITE_MOON || body_id == WALDEMATH
}

/// Whether a hypothetical body is defined in the mean ecliptic of date.
///
/// The `_OF_DATE_HYPOTHETICAL_BODIES` set with the same members, `{VULCAN, WHITE_MOON, PROSERPINA,
/// WALDEMATH}`.
pub fn is_of_date_hypothetical_body(body_id: i32) -> bool {
    body_id == VULCAN || body_id == WHITE_MOON || body_id == PROSERPINA || body_id == WALDEMATH
}

/// Whether a hypothetical body carries a heliocentric predicted model served
/// through the native JPL fallback in every center mode.
///
/// The release `planets._calc_body` `VULCAN`/`PROSERPINA`/`_FICT_HELIO_IDS`
/// arm with the same members: Weston Vulcan (55), Proserpina (57) and the
/// five classical predictions (Harrington, LeVerrier, Adams, Lowell,
/// Pickering). Nibiru (49) is model-less and keeps its existing refusal;
/// Selena/Waldemath stay with the geocentric-symbolic arm and the Uranian
/// bodies with theirs.
pub fn is_predicted_hypothetical_body(body_id: i32) -> bool {
    body_id == VULCAN
        || body_id == PROSERPINA
        || matches!(
            body_id,
            HARRINGTON | NEPTUNE_LEVERRIER | NEPTUNE_ADAMS | PLUTO_LOWELL | PLUTO_PICKERING
        )
}

/// Horizons command for a gravitational-deflector body.
///
/// The `_DEFLECTOR_NAIF` map with the same entries: Sun (`0 -> "10"`), Jupiter barycenter
/// (`5 -> "5"`), Saturn barycenter (`6 -> "6"`). `None` for any other body,
/// mirroring the release `KeyError` the shared deflection routine skips
/// gracefully.
pub fn deflector_command(body_id: i32) -> Option<&'static str> {
    match body_id {
        0 => Some("10"),
        5 => Some("5"),
        6 => Some("6"),
        _ => None,
    }
}

/// Horizons epoch time type for a vectors query.
///
/// Selects the remote `TIME_TYPE` parameter, not a local scale conversion.
/// With [`Self::Tdb`], the supplied Julian date must already be TDB; use the
/// [`crate::time`] conversions when starting from TT. [`Self::Ut`] sends the
/// literal Horizons `UT` option: the service uses UT1 before 1962 and UTC from
/// 1962 onward. For future dates it holds the last known leap-second correction,
/// rather than predicting future leap seconds. The client applies no local
/// UTC/UT1 conversion or Earth-orientation model to the supplied date.
///
/// Source: the [JPL Horizons manual's time-scale conventions](https://ssd.jpl.nasa.gov/horizons/manual.html).
/// The manual discusses TDB, TT and UT vector tables, while the
/// [HTTP API's `TIME_TYPE` description](https://ssd-api.jpl.nasa.gov/doc/horizons.html)
/// lists TDB/UT for vectors. This Rust enum exposes only TDB and UT; it does not
/// promise an additional server option or a TT-to-TDB conversion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HorizonsTimeType {
    /// Barycentric Dynamical Time; input Julian dates are TDB.
    Tdb,
    /// Horizons UT1 before 1962 and UTC from 1962 onward, interpreted remotely.
    /// Future dates retain the service's last known leap-second correction;
    /// see the enum-level JPL references.
    Ut,
}

impl HorizonsTimeType {
    /// Wire spelling (`'TDB'` / `'UT'` before quoting).
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Tdb => "TDB",
            Self::Ut => "UT",
        }
    }
}

/// Percent-encode one vectors-query value.
///
/// Release rule `quote(str(v), safe='')`: bytes outside `[A-Za-z0-9_.~-]`
/// become uppercase `%HH` escapes over the UTF-8 encoding; space becomes
/// `%20` (never `+`). In particular the wrapping single quotes become
/// `%27` and the small-body `';'` becomes `%3B` — a bare `;` would let the
/// server split the query and answer HTTP 400 on every asteroid fetch.
pub fn quote_strict(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for b in value.as_bytes() {
        match b {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' | b'_' | b'.' | b'-' | b'~' => {
                out.push(*b as char);
            }
            _ => {
                out.push('%');
                out.push_str(&format!("{b:02X}"));
            }
        }
    }
    out
}

/// Format a Julian Day for `TLIST`.
///
/// Uses Rust's binary64 debug formatting, retaining a round-trip decimal
/// representation (ordinary integral JDs include `.0`). This does not change
/// the time scale or round the epoch to the cache grid. Nonfinite epochs return
/// [`Error::NonFiniteArgument`] and are not sent.
pub fn format_tlist(jd: f64) -> Result<String> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    Ok(format!("{jd:?}"))
}

/// Build the Horizons vectors request URL.
///
/// Exact `fetch_state_vector` contract: fixed parameter order (`format`,
/// `COMMAND`, `EPHEM_TYPE`, `CENTER`, `TLIST`, `TLIST_TYPE`, `VEC_TABLE`,
/// `OUT_UNITS`, `VEC_CORR`, `CSV_FORMAT`, `REF_SYSTEM`, `REF_PLANE`,
/// `TIME_TYPE`), single-quoted values (except the unquoted `TLIST`), every
/// value [`quote_strict`]-encoded. `REF_PLANE='FRAME'` is load-bearing: it
/// selects the ICRF x-y plane (Earth mean equator); without it Horizons
/// defaults to the ecliptic and the whole pipeline would be misframed by
/// ~23.4°.
///
/// `command` and `center` are Horizons identifier strings; [`horizons_command`]
/// supplies supported library-body mappings. Arbitrary strings, including
/// empties, are encoded without local identifier validation. `jd` is a Julian
/// date already in `time_type`; no TT/TDB/UTC/UT1 conversion occurs. Nonfinite
/// epochs return [`Error::NonFiniteArgument`]. This constructs a URL only and
/// performs no I/O.
pub fn build_vectors_url(
    command: &str,
    jd: f64,
    center: &str,
    time_type: HorizonsTimeType,
) -> Result<String> {
    let params = [
        ("format", "json".to_string()),
        ("COMMAND", format!("'{command}'")),
        ("EPHEM_TYPE", "'VECTORS'".to_string()),
        ("CENTER", format!("'{center}'")),
        ("TLIST", format_tlist(jd)?),
        ("TLIST_TYPE", "'JD'".to_string()),
        ("VEC_TABLE", "'2'".to_string()),
        ("OUT_UNITS", "'AU-D'".to_string()),
        ("VEC_CORR", "'NONE'".to_string()),
        ("CSV_FORMAT", "'YES'".to_string()),
        ("REF_SYSTEM", "'ICRF'".to_string()),
        ("REF_PLANE", "'FRAME'".to_string()),
        ("TIME_TYPE", format!("'{}'", time_type.as_str())),
    ];
    let query: Vec<String> = params
        .iter()
        .map(|(k, v)| format!("{k}={}", quote_strict(v)))
        .collect();
    Ok(format!("{HORIZONS_API_URL}?{}", query.join("&")))
}

/// Geometric ICRF state relative to the requested Horizons center.
///
/// Position is in AU and velocity in AU/day, with `VEC_CORR='NONE'`. Only
/// requests using the SSB center (`"@0"`) are barycentric. No apparent, observer
/// or frame correction is applied here. This value stores finite components,
/// not its target, center or epoch metadata; retain the request alongside it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StateVector {
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
}

impl StateVector {
    /// Checked constructor: every component must be finite.
    ///
    /// A state vector with a non-finite component is malformed on arrival;
    /// infinity/NaN must never flow into the POS pipeline as geometry.
    pub fn new(x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) -> Result<Self> {
        for value in [x, y, z, vx, vy, vz] {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        Ok(Self {
            x,
            y,
            z,
            vx,
            vy,
            vz,
        })
    }

    /// Position components in AU (ICRS, requested center).
    pub fn pos(&self) -> (f64, f64, f64) {
        (self.x, self.y, self.z)
    }

    /// Velocity components in AU/day (ICRS, requested center).
    pub fn vel(&self) -> (f64, f64, f64) {
        (self.vx, self.vy, self.vz)
    }

    /// Position x in AU.
    pub fn x(&self) -> f64 {
        self.x
    }

    /// Position y in AU.
    pub fn y(&self) -> f64 {
        self.y
    }

    /// Position z in AU.
    pub fn z(&self) -> f64 {
        self.z
    }

    /// Velocity x in AU/day.
    pub fn vx(&self) -> f64 {
        self.vx
    }

    /// Velocity y in AU/day.
    pub fn vy(&self) -> f64 {
        self.vy
    }

    /// Velocity z in AU/day.
    pub fn vz(&self) -> f64 {
        self.vz
    }
}

/// Locate a top-level `"key":` field in a flat JSON envelope.
///
/// Returns the byte offset just past the colon, or `None` when the key is
/// absent. This is an envelope scanner, not a JSON parser: it only serves
/// the Horizons vectors envelope (`{"error"?, "result"?,...}`), where both
/// fields of interest are top-level.
fn top_field_value_start(body: &str, key: &str) -> Option<usize> {
    let needle = format!("\"{key}\"");
    let bytes = body.as_bytes();
    let mut search_from = 0;
    while let Some(rel) = body[search_from..].find(&needle) {
        let mut idx = search_from + rel + needle.len();
        while idx < bytes.len() && bytes[idx].is_ascii_whitespace() {
            idx += 1;
        }
        if idx < bytes.len() && bytes[idx] == b':' {
            return Some(idx + 1);
        }
        search_from = search_from + rel + 1;
    }
    None
}

/// Decode one JSON string literal starting at `start` (the opening quote).
///
/// Returns the decoded text and the byte offset just past the closing
/// quote. Handles all standard escapes including `\uXXXX` with surrogate
/// pairs; rejects control characters, bad escapes, lone surrogates and
/// truncation with a static reason. Never lossy: undecodable input is an
/// error, never replacement characters.
fn decode_json_string(
    body: &str,
    start: usize,
) -> std::result::Result<(String, usize), &'static str> {
    let bytes = body.as_bytes();
    if bytes.get(start) != Some(&b'"') {
        return Err("expected string value");
    }
    let mut out = String::new();
    let mut idx = start + 1;
    while idx < bytes.len() {
        let b = bytes[idx];
        match b {
            b'"' => return Ok((out, idx + 1)),
            b'\\' => {
                idx += 1;
                let esc = *bytes.get(idx).ok_or("truncated string escape")?;
                match esc {
                    b'"' => out.push('"'),
                    b'\\' => out.push('\\'),
                    b'/' => out.push('/'),
                    b'b' => out.push('\u{0008}'),
                    b'f' => out.push('\u{000C}'),
                    b'n' => out.push('\n'),
                    b'r' => out.push('\r'),
                    b't' => out.push('\t'),
                    b'u' => {
                        if idx + 4 >= bytes.len() {
                            return Err("truncated unicode escape");
                        }
                        let hex = &body[idx + 1..idx + 5];
                        let unit =
                            u32::from_str_radix(hex, 16).map_err(|_| "invalid unicode escape")?;
                        idx += 4;
                        let ch = if (0xD800..0xDC00).contains(&unit) {
                            // High surrogate: a low surrogate must follow.
                            if body[idx + 1..].starts_with("\\u") && idx + 7 < bytes.len() {
                                let low_hex = &body[idx + 3..idx + 7];
                                let low = u32::from_str_radix(low_hex, 16)
                                    .map_err(|_| "invalid low surrogate")?;
                                if !(0xDC00..0xE000).contains(&low) {
                                    return Err("lone high surrogate");
                                }
                                idx += 6;
                                let scalar = 0x10000 + ((unit - 0xD800) << 10) + (low - 0xDC00);
                                char::from_u32(scalar).ok_or("invalid code point")?
                            } else {
                                return Err("lone high surrogate");
                            }
                        } else if (0xDC00..0xE000).contains(&unit) {
                            return Err("lone low surrogate");
                        } else {
                            char::from_u32(unit).ok_or("invalid code point")?
                        };
                        out.push(ch);
                    }
                    _ => return Err("invalid string escape"),
                }
                idx += 1;
            }
            0x00..=0x1F => return Err("unescaped control character"),
            _ => {
                // Copy one full UTF-8 scalar (literal non-ASCII is legal JSON).
                let rest = &body[idx..];
                let ch = rest.chars().next().ok_or("truncated string")?;
                out.push(ch);
                idx += ch.len_utf8();
            }
        }
    }
    Err("unterminated string value")
}

/// Extract a top-level string field from the Horizons envelope.
///
/// `Ok(None)` when the key is absent; `Ok(Some(_))` with the decoded text
/// when present; `Err(reason)` when the value is truncated or undecodable.
/// A present-but-non-string value (never sent by Horizons) reports
/// `Err("expected string value")` so it cannot be mistaken for absence.
fn top_string_field(body: &str, key: &str) -> std::result::Result<Option<String>, &'static str> {
    let Some(mut idx) = top_field_value_start(body, key) else {
        return Ok(None);
    };
    let bytes = body.as_bytes();
    while idx < bytes.len() && bytes[idx].is_ascii_whitespace() {
        idx += 1;
    }
    if bytes.get(idx) != Some(&b'"') {
        return Err("expected string value");
    }
    decode_json_string(body, idx).map(|(text, _)| Some(text))
}

/// Parse one Horizons vectors response body.
///
/// Exact `_parse_response` contract over the response text:
///
/// 1. A top-level `"error"` field is [`Error::HorizonsApiError`] (never
///    retried, mirroring the release `raise KeyError` before parsing).
/// 2. Otherwise the `"result"` text must hold `$$SOE`/`$$EOE` markers
///    around one `VEC_TABLE='2'` + `CSV_FORMAT='YES'` record shaped
///    `JDTDB, Calendar Date, X, Y, Z, VX, VY, VZ`, parsed positionally
///    from the first non-empty line: fields 2–8 are the state vector.
///    Marker/record/field/numeric failures are
///    [`Error::HorizonsMalformed`] (retried, mirroring the release
///    `ValueError` path).
///
/// `command` names the requested body for error attribution only.
pub fn parse_vectors_response(body: &str, command: &str) -> Result<StateVector> {
    if let Some(error_text) =
        top_string_field(body, "error").map_err(|reason| Error::HorizonsMalformed {
            command: command.to_string(),
            reason,
        })?
    {
        return Err(Error::HorizonsApiError {
            command: command.to_string(),
            message: error_text,
        });
    }
    let result_text =
        top_string_field(body, "result").map_err(|reason| Error::HorizonsMalformed {
            command: command.to_string(),
            reason,
        })?;
    let Some(result_text) = result_text else {
        return Err(Error::HorizonsMalformed {
            command: command.to_string(),
            reason: "missing result field",
        });
    };
    let soe = result_text.find("$$SOE").ok_or(Error::HorizonsMalformed {
        command: command.to_string(),
        reason: "missing $$SOE marker",
    })?;
    let eoe = result_text.find("$$EOE").ok_or(Error::HorizonsMalformed {
        command: command.to_string(),
        reason: "missing $$EOE marker",
    })?;
    if eoe < soe {
        // Cannot happen from the server (markers are always ordered); the
        // release slices an empty block here and reports it empty below.
        // The saturating slice keeps that exact observable outcome.
        return Err(Error::HorizonsMalformed {
            command: command.to_string(),
            reason: "empty data block",
        });
    }
    let block = result_text[soe + "$$SOE".len()..eoe].trim();
    let line = block
        .lines()
        .map(str::trim)
        .find(|line| !line.is_empty())
        .ok_or(Error::HorizonsMalformed {
            command: command.to_string(),
            reason: "empty data block",
        })?;
    // Positional CSV parse of the single-epoch record: accumulating every
    // parseable float across the block would let a stray numeric token
    // shift the coordinate indices silently.
    let fields: Vec<&str> = line.split(',').map(str::trim).collect();
    if fields.len() < 8 {
        return Err(Error::HorizonsMalformed {
            command: command.to_string(),
            reason: "short data record",
        });
    }
    let mut components = [0.0f64; 6];
    for (slot, field) in components.iter_mut().zip(fields.iter().skip(2).take(6)) {
        *slot = field.parse::<f64>().map_err(|_| Error::HorizonsMalformed {
            command: command.to_string(),
            reason: "non-numeric state-vector field",
        })?;
        if !slot.is_finite() {
            return Err(Error::HorizonsMalformed {
                command: command.to_string(),
                reason: "non-finite state-vector field",
            });
        }
    }
    StateVector::new(
        components[0],
        components[1],
        components[2],
        components[3],
        components[4],
        components[5],
    )
}

/// Normalize a finite Julian Day to twelve decimal places for cache indexing.
/// The key also includes command, center and time scale. Binary64 spacing
/// at astronomical JD magnitudes is coarser than 1e-12 day, so this operation
/// does not promise to combine distinct representable epochs. Non-finite
/// epochs return [`Error::NonFiniteArgument`] and are never cached.
pub fn round_jd12(jd: f64) -> Result<f64> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    Ok((jd * 1e12).round_ties_even() / 1e12)
}

/// Identity of one cached Horizons vector.
///
/// The full release cache key `(round(jd, 12), command, center, time_type)`:
/// a cached vector is served only to its exact command, center and time
/// type at a key-equal epoch. The epoch is stored as the bits of the
/// [`round_jd12`] normalization so the key is `Hash + Eq` (raw `f64`
/// cannot key a map: NaN never equals itself).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct HorizonsCacheKey {
    jd_key_bits: u64,
    command: String,
    center: String,
    time_type: HorizonsTimeType,
}

impl HorizonsCacheKey {
    /// Build the cache key for one query.
    ///
    /// Command and center spellings are kept as given (empties pass
    /// through to the server like the release; they simply key their own
    /// entry and can only ever hit on an identical re-query).
    pub fn new(command: &str, jd: f64, center: &str, time_type: HorizonsTimeType) -> Result<Self> {
        Ok(Self {
            jd_key_bits: round_jd12(jd)?.to_bits(),
            command: command.to_string(),
            center: center.to_string(),
            time_type,
        })
    }

    /// Cache-normalized epoch this key serves.
    pub fn jd_key(&self) -> f64 {
        f64::from_bits(self.jd_key_bits)
    }

    /// Horizons command this key serves.
    pub fn command(&self) -> &str {
        &self.command
    }

    /// Horizons center this key serves.
    pub fn center(&self) -> &str {
        &self.center
    }

    /// Epoch time type this key serves.
    pub fn time_type(&self) -> HorizonsTimeType {
        self.time_type
    }
}

/// Bounded least-recently-used Horizons vector store.
///
/// Pure `OrderedDict` equivalent: `lookup` refreshes recency, `insert` evicts the oldest entry past capacity.
/// Thread safety lives in the owning [`HorizonsClient`] mutex; this type
/// itself is a plain value so the eviction order stays unit-testable
/// without sockets.
#[derive(Debug, Clone, PartialEq)]
pub struct HorizonsCache {
    capacity: usize,
    entries: HashMap<HorizonsCacheKey, StateVector>,
    order: VecDeque<HorizonsCacheKey>,
}

impl HorizonsCache {
    /// Empty cache holding at most `capacity` entries (`capacity == 0`
    /// admits no entry; clients reject that shape at construction via
    /// `NonZeroUsize`).
    pub fn new(capacity: usize) -> Self {
        Self {
            capacity,
            entries: HashMap::new(),
            order: VecDeque::new(),
        }
    }

    /// Capacity in entries.
    pub fn capacity(&self) -> usize {
        self.capacity
    }

    /// Number of stored vectors.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the cache holds no vector.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Whether `key` is stored (no recency change: use [`HorizonsCache::lookup`]
    /// for the release `move_to_end` hit path).
    pub fn contains(&self, key: &HorizonsCacheKey) -> bool {
        self.entries.contains_key(key)
    }

    /// Serve the cached vector and refresh its recency.
    pub fn lookup(&mut self, key: &HorizonsCacheKey) -> Option<StateVector> {
        let hit = self.entries.get(key).copied()?;
        if let Some(pos) = self.order.iter().position(|member| member == key) {
            self.order.remove(pos);
        }
        self.order.push_back(key.clone());
        Some(hit)
    }

    /// Store `vector` under `key`, evicting oldest entries past capacity
    pub fn insert(&mut self, key: HorizonsCacheKey, vector: StateVector) {
        if self.entries.contains_key(&key) {
            self.entries.insert(key.clone(), vector);
            if let Some(pos) = self.order.iter().position(|member| member == &key) {
                self.order.remove(pos);
            }
            self.order.push_back(key);
            return;
        }
        self.entries.insert(key.clone(), vector);
        self.order.push_back(key);
        while self.entries.len() > self.capacity {
            if let Some(oldest) = self.order.pop_front() {
                self.entries.remove(&oldest);
            } else {
                break;
            }
        }
    }

    /// Drop every entry.
    pub fn clear(&mut self) {
        self.entries.clear();
        self.order.clear();
    }

    /// Snapshot of the stored entries (clone support for the client).
    fn snapshot(&self) -> Self {
        self.clone()
    }
}

/// One batch vector request: command, epoch and center.
///
/// [`HorizonsClient::fetch_batch`] interprets `jd` as TDB Julian days and
/// hardcodes [`HorizonsTimeType::Tdb`]. It does not convert TT or accept an
/// alternative time scale. The center is a Horizons identifier, not a library
/// body ID. Construction validates epoch finiteness only.
#[derive(Debug, Clone, PartialEq)]
pub struct VectorsRequest {
    command: String,
    jd: f64,
    center: String,
}

impl VectorsRequest {
    /// Checked constructor: the epoch must be finite (a non-finite epoch
    /// can never be fetched nor cached, so it is refused up front rather
    /// than sent to the server).
    pub fn new(command: &str, jd: f64, center: &str) -> Result<Self> {
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        Ok(Self {
            command: command.to_string(),
            jd,
            center: center.to_string(),
        })
    }

    /// Horizons command to fetch.
    pub fn command(&self) -> &str {
        &self.command
    }

    /// Requested TDB Julian date, before cache-key normalization.
    pub fn jd(&self) -> f64 {
        self.jd
    }

    /// Horizons center to fetch.
    pub fn center(&self) -> &str {
        &self.center
    }

    /// Build a key using the supplied time-type tag, without converting `jd`.
    /// Batch execution always uses TDB regardless of separately constructed keys.
    pub fn cache_key(&self, time_type: HorizonsTimeType) -> Result<HorizonsCacheKey> {
        HorizonsCacheKey::new(&self.command, self.jd, &self.center, time_type)
    }

    /// Result-map key: the original request identity.
    pub fn batch_key(&self) -> BatchKey {
        BatchKey {
            command: self.command.clone(),
            jd_bits: self.jd.to_bits(),
            center: self.center.clone(),
        }
    }
}

/// Identity of one batch result row.
///
/// The requested `(command, jd, center)` triple with the epoch as bits
/// (raw `f64` is not `Hash`; every request epoch is finite by
/// [`VectorsRequest::new`], so bits round-trip exactly).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct BatchKey {
    command: String,
    jd_bits: u64,
    center: String,
}

impl BatchKey {
    /// Requested Horizons command.
    pub fn command(&self) -> &str {
        &self.command
    }

    /// Requested Julian Day (original spelling).
    pub fn jd(&self) -> f64 {
        f64::from_bits(self.jd_bits)
    }

    /// Requested Horizons center.
    pub fn center(&self) -> &str {
        &self.center
    }
}

/// Owned Horizons network configuration: timeout, cache capacity and batch
/// worker bound.
///
/// `timeout_secs` is the per-exchange timeout, not a whole-batch deadline;
/// retries and backoff can extend total blocking time. `max_cache_size` bounds
/// retained vector entries, not process memory. `max_workers` bounds threads
/// within one batch call, not all concurrent client calls. The retry budget
/// stays fixed at [`HORIZONS_MAX_RETRIES`].
///
/// Construction validates eagerly: a non-finite or non-positive timeout is
/// [`Error::InvalidTimeout`]; zero capacities are unrepresentable
/// (`NonZeroUsize`), so a zero-entry cache or a zero-worker pool is a type
/// error rather than a silently useless client. Building a client from this
/// value performs no I/O (the cache starts empty), so attaching a config to
/// a context is observably network-free; every exchange still builds a
/// per-call agent in [`crate::net`], so the cache holds result vectors
/// only, never pooled connection state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HorizonsConfig {
    timeout_secs: f64,
    max_cache_size: NonZeroUsize,
    max_workers: NonZeroUsize,
}

impl HorizonsConfig {
    /// Explicit configuration with release constructor semantics.
    ///
    /// Non-finite or non-positive `timeout_secs` is
    /// [`Error::InvalidTimeout`]; nothing is silently clamped.
    pub fn new(
        timeout_secs: f64,
        max_cache_size: NonZeroUsize,
        max_workers: NonZeroUsize,
    ) -> Result<Self> {
        if !timeout_secs.is_finite() || timeout_secs <= 0.0 {
            return Err(Error::InvalidTimeout {
                seconds: timeout_secs,
            });
        }
        Ok(Self {
            timeout_secs,
            max_cache_size,
            max_workers,
        })
    }

    /// Query timeout in seconds.
    pub fn timeout_secs(&self) -> f64 {
        self.timeout_secs
    }

    /// LRU cache capacity in entries.
    pub fn max_cache_size(&self) -> NonZeroUsize {
        self.max_cache_size
    }

    /// Batch worker bound.
    pub fn max_workers(&self) -> NonZeroUsize {
        self.max_workers
    }
}

impl Default for HorizonsConfig {
    /// Release defaults: 30 s timeout, 4096-entry cache, 8 batch workers
    fn default() -> Self {
        Self::new(
            HORIZONS_DEFAULT_TIMEOUT_SECS,
            NonZeroUsize::new(HORIZONS_DEFAULT_MAX_CACHE_SIZE)
                .expect("default Horizons cache size is nonzero"),
            NonZeroUsize::new(HORIZONS_DEFAULT_MAX_WORKERS)
                .expect("default Horizons worker count is nonzero"),
        )
        .expect("release-default Horizons config is valid")
    }
}

/// Owned Horizons vectors client: timeout, retry, cache and batch policy.
///
/// Sharing a client shares its mutex-protected result cache. Cloning snapshots
/// configuration and existing entries into an independent cache; subsequent
/// inserts and clears are not shared. Exchanges build per-call transport
/// agents, so the cache holds vectors rather than pooled connections.
///
/// Calls are synchronous and can block on transport, retries and backoff.
/// Concurrent misses are not coalesced, and each batch has its own worker
/// bound. [`Self::shutdown`] clears entries but does not cancel or drain work.
pub struct HorizonsClient {
    timeout_secs: f64,
    max_retries: u32,
    max_cache_size: usize,
    max_workers: usize,
    cache: Mutex<HorizonsCache>,
}

impl std::fmt::Debug for HorizonsClient {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("HorizonsClient")
            .field("timeout_secs", &self.timeout_secs)
            .field("max_retries", &self.max_retries)
            .field("max_cache_size", &self.max_cache_size)
            .field("max_workers", &self.max_workers)
            .field("cache_len", &self.cache_len())
            .finish_non_exhaustive()
    }
}

impl Clone for HorizonsClient {
    /// Snapshot clone: configuration plus the currently cached vectors.
    /// Later inserts on either clone stay local (no shared entry can be
    /// served across clones after the split).
    fn clone(&self) -> Self {
        Self {
            timeout_secs: self.timeout_secs,
            max_retries: self.max_retries,
            max_cache_size: self.max_cache_size,
            max_workers: self.max_workers,
            cache: Mutex::new(self.lock_cache().snapshot()),
        }
    }
}

impl PartialEq for HorizonsClient {
    /// Configuration plus current cache contents (a clone compares equal
    /// until either side caches or clears).
    fn eq(&self, other: &Self) -> bool {
        if std::ptr::eq(self, other) {
            return true;
        }
        if !(self.timeout_secs == other.timeout_secs
            && self.max_retries == other.max_retries
            && self.max_cache_size == other.max_cache_size
            && self.max_workers == other.max_workers)
        {
            return false;
        }
        // Shared context owners can refer to this exact client. Distinct
        // clients lock in address order, so opposite concurrent comparisons
        // cannot acquire the two cache locks in reverse order.
        let (first, second) = if std::ptr::from_ref(self) < std::ptr::from_ref(other) {
            (self, other)
        } else {
            (other, self)
        };
        *first.lock_cache() == *second.lock_cache()
    }
}

impl HorizonsClient {
    /// Client with the release defaults: 30 s timeout, 2 retries, a 4096-entry
    /// LRU cache and 8 batch workers.
    pub fn new() -> Result<Self> {
        Self::with_timeout_secs(HORIZONS_DEFAULT_TIMEOUT_SECS)
    }

    /// Client with an explicit timeout (seconds) and the release retry budget,
    /// cache size and worker count.
    ///
    /// Non-finite or non-positive values are [`Error::InvalidTimeout`];
    /// nothing is silently clamped.
    pub fn with_timeout_secs(timeout_secs: f64) -> Result<Self> {
        Self::with_config(
            timeout_secs,
            NonZeroUsize::new(HORIZONS_DEFAULT_MAX_CACHE_SIZE)
                .expect("default Horizons cache size is nonzero"),
            NonZeroUsize::new(HORIZONS_DEFAULT_MAX_WORKERS)
                .expect("default Horizons worker count is nonzero"),
        )
    }

    /// Client with explicit timeout, cache capacity and batch worker count.
    ///
    /// `NonZeroUsize` bounds carry the release "no validation" gap
    /// honestly: a zero-entry cache or a zero-worker pool is a type error
    /// here rather than a silently useless client.
    pub fn with_config(
        timeout_secs: f64,
        max_cache_size: NonZeroUsize,
        max_workers: NonZeroUsize,
    ) -> Result<Self> {
        if !timeout_secs.is_finite() || timeout_secs <= 0.0 {
            return Err(Error::InvalidTimeout {
                seconds: timeout_secs,
            });
        }
        Ok(Self {
            timeout_secs,
            max_retries: HORIZONS_MAX_RETRIES,
            max_cache_size: max_cache_size.get(),
            max_workers: max_workers.get(),
            cache: Mutex::new(HorizonsCache::new(max_cache_size.get())),
        })
    }

    /// Client from an explicit [`HorizonsConfig`] value.
    ///
    /// Same validation as [`HorizonsClient::with_config`]; the retry budget
    /// stays fixed. Construction performs no I/O (the cache starts empty).
    pub fn from_config(config: HorizonsConfig) -> Result<Self> {
        Self::with_config(
            config.timeout_secs,
            config.max_cache_size,
            config.max_workers,
        )
    }

    /// Query timeout in seconds.
    pub fn timeout_secs(&self) -> f64 {
        self.timeout_secs
    }

    /// Retry budget after the first attempt, fixed at [`HORIZONS_MAX_RETRIES`].
    pub fn max_retries(&self) -> u32 {
        self.max_retries
    }

    /// LRU cache capacity in entries.
    pub fn max_cache_size(&self) -> usize {
        self.max_cache_size
    }

    /// Batch worker bound.
    pub fn max_workers(&self) -> usize {
        self.max_workers
    }

    /// Number of currently cached vectors.
    pub fn cache_len(&self) -> usize {
        self.lock_cache().len()
    }

    /// Whether the exact query is currently cached (no fetch performed).
    pub fn cache_contains(
        &self,
        command: &str,
        jd: f64,
        center: &str,
        time_type: HorizonsTimeType,
    ) -> Result<bool> {
        Ok(self
            .lock_cache()
            .contains(&HorizonsCacheKey::new(command, jd, center, time_type)?))
    }

    /// Lock the result cache (poison-tolerant: a poisoned mutex still
    /// yields its guarded cache, since entries are plain data).
    fn lock_cache(&self) -> std::sync::MutexGuard<'_, HorizonsCache> {
        self.cache.lock().unwrap_or_else(|err| err.into_inner())
    }

    /// Drop every cached vector.
    ///
    /// Configuration and the client remain usable. This is not cancellation,
    /// a worker drain or an invalidation barrier: in-flight successful fetches
    /// may insert entries after the clear.
    pub fn clear_cache(&self) {
        self.lock_cache().clear();
    }

    /// Clear cached vectors without permanently closing the client.
    ///
    /// Equivalent to [`Self::clear_cache`]. It neither cancels nor waits for
    /// requests, and in-flight work may repopulate the cache. Further fetches
    /// are allowed; callers own shutdown coordination.
    pub fn shutdown(&self) {
        self.clear_cache();
    }

    /// Build the vectors URL for one query.
    ///
    /// Thin wrapper over [`build_vectors_url`] with the default center
    /// ([`HORIZONS_DEFAULT_CENTER`]) selection left to the caller.
    pub fn vectors_url(
        &self,
        command: &str,
        jd: f64,
        center: &str,
        time_type: HorizonsTimeType,
    ) -> Result<String> {
        build_vectors_url(command, jd, center, time_type)
    }

    /// Run one vectors exchange without retry.
    ///
    /// The sealed-policy gate fires before any socket work (via
    /// [`crate::net`]); HTTP 4xx/5xx surface as [`Error::HttpStatus`],
    /// API refusals as [`Error::HorizonsApiError`], bad records as
    /// [`Error::HorizonsMalformed`]. The exchange carries the release
    /// purpose string (`Horizons state vector for {command}`) with the
    /// caller `User-Agent` header; TLS uses the [`crate::net`] per-call
    /// agent rather than the release module-level cached SSL context.
    /// Callers needing the release retry schedule use
    /// [`HorizonsClient::fetch_vectors`].
    pub fn fetch_once(
        &self,
        url: &str,
        command: &str,
        policy: EffectiveNetworkPolicy,
    ) -> Result<StateVector> {
        let purpose = format!("Horizons state vector for {command}");
        let request = HttpRequest::get(url)?.header("User-Agent", HORIZONS_USER_AGENT);
        let config = TransportConfig::with_timeout_secs(self.timeout_secs)?;
        let response = fetch(&request, &purpose, policy, &config)?;
        parse_vectors_response(response.body_text()?, command)
    }

    /// Fetch a requested-center geometric state with bounded blocking retries.
    ///
    /// `jd` must already use `time_type`. The returned ICRF components are AU
    /// and AU/day; neither clock conversion nor center translation is applied.
    /// Up to three attempts use 0.5 s / 1.0 s backoff. API refusals and caller
    /// errors return immediately; exhausted retryable failures become
    /// [`Error::HorizonsRequestFailed`] carrying the last failure's text.
    ///
    /// This primitive does not consult or populate the cache. An invalid input
    /// or sealed policy can refuse before sockets. Use [`Self::fetch_vectors_cached`]
    /// for cache-first execution.
    pub fn fetch_vectors(
        &self,
        command: &str,
        jd: f64,
        center: &str,
        time_type: HorizonsTimeType,
        policy: EffectiveNetworkPolicy,
    ) -> Result<StateVector> {
        self.fetch_uncached(command, jd, center, time_type, policy)
    }

    /// Uncached URL fetch with the retry schedule, shared by the cached
    /// single-fetch and batch paths below.
    fn fetch_uncached(
        &self,
        command: &str,
        jd: f64,
        center: &str,
        time_type: HorizonsTimeType,
        policy: EffectiveNetworkPolicy,
    ) -> Result<StateVector> {
        let url = self.vectors_url(command, jd, center, time_type)?;
        self.fetch_loop(&url, command, policy)
    }

    /// Fetch one state vector through the LRU cache.
    ///
    /// Input clocks, requested-center units and errors follow [`Self::fetch_vectors`].
    /// A hit performs no socket work and bypasses the network-policy gate;
    /// sealed callers can therefore receive an already cached result. Misses
    /// fetch and store successes only. Separate clients are required when
    /// cached provenance must not cross application policy boundaries.
    /// Concurrent misses may fetch the same key independently.
    pub fn fetch_vectors_cached(
        &self,
        command: &str,
        jd: f64,
        center: &str,
        time_type: HorizonsTimeType,
        policy: EffectiveNetworkPolicy,
    ) -> Result<StateVector> {
        let key = HorizonsCacheKey::new(command, jd, center, time_type)?;
        if let Some(hit) = self.lock_cache().lookup(&key) {
            return Ok(hit);
        }
        let vector = self.fetch_uncached(command, jd, center, time_type, policy)?;
        self.lock_cache().insert(key, vector);
        Ok(vector)
    }

    /// Cache-aware single fetch over a caller-supplied fetcher.
    ///
    /// `fetch` runs synchronously on a miss, outside the cache lock; successes
    /// are cached and errors are not. The closure owns I/O, clock/center
    /// consistency and network policy. This adapter performs no direct network
    /// request and does not catch callback panics. Concurrent misses can each
    /// invoke their own closure.
    ///
    /// ```
    /// use libephemeris_core::horizons::{HorizonsClient, HorizonsTimeType, StateVector};
    ///
    /// let client = HorizonsClient::new()?;
    /// client.get_or_fetch_with("499", 2_451_545.0, "@0", HorizonsTimeType::Tdb,
    ///     || StateVector::new(1.0, 0.0, 0.0, 0.0, 0.01, 0.0))?;
    /// assert_eq!(client.cache_len(), 1);
    /// client.shutdown(); // clears entries; does not close the client
    /// assert_eq!(client.cache_len(), 0);
    /// # Ok::<(), libephemeris_core::errors::Error>(())
    /// ```
    pub fn get_or_fetch_with(
        &self,
        command: &str,
        jd: f64,
        center: &str,
        time_type: HorizonsTimeType,
        fetch: impl FnOnce() -> Result<StateVector>,
    ) -> Result<StateVector> {
        let key = HorizonsCacheKey::new(command, jd, center, time_type)?;
        if let Some(hit) = self.lock_cache().lookup(&key) {
            return Ok(hit);
        }
        let vector = fetch()?;
        self.lock_cache().insert(key, vector);
        Ok(vector)
    }

    /// Fetch multiple TDB state vectors with cache-first parallelism.
    ///
    /// Every request epoch must already be TDB; there is no TT conversion.
    /// Cache hits bypass sockets and the policy gate. Misses run on at most
    /// [`Self::max_workers`] scoped threads per call; the call waits for all
    /// workers. Results and typed failures occupy separate maps keyed by the
    /// original `(command, jd bits, center)`, not the normalized cache epoch.
    /// Empty input returns two empty maps. Duplicate request keys are not
    /// separate result rows and misses are not deduplicated; completion order
    /// is unspecified. Returned states retain each request's center.
    ///
    /// # Panics
    /// Scoped worker panics propagate to the caller rather than becoming typed
    /// errors in the result map. Thread creation can also panic.
    pub fn fetch_batch(
        &self,
        requests: &[VectorsRequest],
        policy: EffectiveNetworkPolicy,
    ) -> (HashMap<BatchKey, StateVector>, HashMap<BatchKey, Error>) {
        self.fetch_batch_with(requests, |command, jd, center| {
            self.fetch_uncached(command, jd, center, HorizonsTimeType::Tdb, policy)
        })
    }

    /// Batch orchestration over a caller-supplied fetcher.
    ///
    /// Uses [`Self::fetch_batch`]'s TDB cache keys, duplicate handling and
    /// per-call worker bound. `fetch(command, jd_tdb, center)` can run
    /// concurrently outside the cache lock and must return matching geometric
    /// ICRF AU/AU-day states. It owns all I/O and policy decisions; this adapter
    /// neither enforces network policy nor assumes the closure is pure.
    /// Typed errors are collected per key and successes are cached.
    ///
    /// # Panics
    /// Callback and other scoped-worker panics propagate after scoped workers
    /// have been joined; they are not converted into the returned error map.
    /// Thread creation can also panic.
    pub fn fetch_batch_with<F>(
        &self,
        requests: &[VectorsRequest],
        fetch: F,
    ) -> (HashMap<BatchKey, StateVector>, HashMap<BatchKey, Error>)
    where
        F: Fn(&str, f64, &str) -> Result<StateVector> + Sync,
    {
        // Partition under one short lock: hits are served, misses queued.
        // The lock is never held across a fetch.
        let mut results: HashMap<BatchKey, StateVector> = HashMap::new();
        let mut misses: Vec<(BatchKey, HorizonsCacheKey, &VectorsRequest)> = Vec::new();
        {
            let mut cache = self.lock_cache();
            for request in requests {
                // `VectorsRequest::new` guarantees a finite epoch, so the
                // TDB cache key below is infallible; a defensive miss keeps
                // a corrupt request from poisoning the batch.
                let Ok(key) = request.cache_key(HorizonsTimeType::Tdb) else {
                    continue;
                };
                match cache.lookup(&key) {
                    Some(hit) => {
                        results.insert(request.batch_key(), hit);
                    }
                    None => misses.push((request.batch_key(), key, request)),
                }
            }
        }
        // Fetch the misses on at most `max_workers` scoped threads.
        let fetched: Mutex<Vec<(BatchKey, HorizonsCacheKey, Result<StateVector>)>> =
            Mutex::new(Vec::with_capacity(misses.len()));
        if !misses.is_empty() {
            let chunks = misses.len().div_ceil(self.max_workers).max(1);
            let fetch_ref = &fetch;
            let fetched_ref = &fetched;
            std::thread::scope(|scope| {
                for chunk in misses.chunks(chunks) {
                    scope.spawn(move || {
                        let mut local = Vec::with_capacity(chunk.len());
                        for (batch_key, cache_key, request) in chunk {
                            local.push((
                                batch_key.clone(),
                                cache_key.clone(),
                                fetch_ref(request.command(), request.jd(), request.center()),
                            ));
                        }
                        fetched_ref
                            .lock()
                            .unwrap_or_else(|err| err.into_inner())
                            .extend(local);
                    });
                }
            });
        }
        // Merge under one short lock: successes are stored and reported; failures are captured per
        // request and never cached.
        let mut errors: HashMap<BatchKey, Error> = HashMap::new();
        {
            let mut cache = self.lock_cache();
            for (batch_key, cache_key, outcome) in
                fetched.into_inner().unwrap_or_else(|err| err.into_inner())
            {
                match outcome {
                    Ok(vector) => {
                        cache.insert(cache_key, vector);
                        results.insert(batch_key, vector);
                    }
                    Err(err) => {
                        errors.insert(batch_key, err);
                    }
                }
            }
        }
        (results, errors)
    }

    /// Retry loop over an explicit URL, shared by [`HorizonsClient::fetch_vectors`]
    /// and the hermetic loopback tests (which cannot use the release endpoint).
    fn fetch_loop(
        &self,
        url: &str,
        command: &str,
        policy: EffectiveNetworkPolicy,
    ) -> Result<StateVector> {
        let mut attempt: u32 = 0;
        loop {
            match self.fetch_once(url, command, policy) {
                Ok(vector) => return Ok(vector),
                Err(err) => {
                    if !is_retryable_horizons_error(&err) || attempt >= self.max_retries {
                        if !is_retryable_horizons_error(&err) {
                            return Err(err);
                        }
                        return Err(Error::HorizonsRequestFailed {
                            command: command.to_string(),
                            attempts: attempt + 1,
                            detail: err.to_string(),
                        });
                    }
                    std::thread::sleep(Duration::from_secs_f64(retry_delay_secs(attempt)));
                    attempt += 1;
                }
            }
        }
    }
}

impl Default for HorizonsClient {
    fn default() -> Self {
        Self::new().expect("default Horizons timeout is valid")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::context::{EphemerisContext, NetworkPolicy, ProviderMode};
    use std::io::Read as _;
    use std::io::Write as _;
    use std::net::TcpListener;
    use std::time::Duration;

    fn allow() -> EffectiveNetworkPolicy {
        EffectiveNetworkPolicy::Allow
    }

    #[test]
    fn shared_and_opposite_client_comparisons_do_not_relock_one_cache() {
        let first = HorizonsClient::new().unwrap();
        let second = first.clone();
        assert_eq!(first, first);
        std::thread::scope(|scope| {
            for reverse in [false, true] {
                let (a, b) = if reverse {
                    (&second, &first)
                } else {
                    (&first, &second)
                };
                scope.spawn(move || {
                    for _ in 0..1000 {
                        assert_eq!(a, b);
                    }
                });
            }
        });
    }

    #[test]
    fn command_map_matches_release_pins() {
        // Every release _HORIZONS_COMMAND entry,.
        let pins = [
            (0, "10"),
            (1, "301"),
            (2, "199"),
            (3, "299"),
            (4, "499"),
            (5, "599"),
            (6, "699"),
            (7, "799"),
            (8, "899"),
            (9, "999"),
            (14, "399"),
            (15, "2060;"),
            (17, "1;"),
            (18, "2;"),
            (19, "3;"),
            (20, "4;"),
        ];
        for (body, command) in pins {
            assert_eq!(horizons_command(body).unwrap(), command, "body {body}");
        }
        // No other body resolves: Pholus (16) stays absent like the release.
        for body in [16, 21, 22, 40, -1, 999] {
            assert!(
                matches!(
                    horizons_command(body),
                    Err(Error::UnknownBody { body_id }) if body_id == body
                ),
                "body {body} must be UnknownBody"
            );
        }
    }

    #[test]
    fn body_class_sets_match_release() {
        assert!(is_analytical_body(10) && is_analytical_body(12));
        assert!(!is_analytical_body(11) && !is_analytical_body(0));
        for body in [11, 13, 21, 22] {
            assert!(is_moon_derived_body(body), "body {body}");
        }
        assert!(!is_moon_derived_body(1) && !is_moon_derived_body(10));
        // Hypothetical range 40..=58 with the release subsets.
        assert!(!is_hypothetical_body(39) && !is_hypothetical_body(59));
        for body in FICT_OFFSET..=WALDEMATH {
            assert!(is_hypothetical_body(body), "body {body}");
        }
        assert!(is_geocentric_hypothetical_body(WHITE_MOON));
        assert!(is_geocentric_hypothetical_body(WALDEMATH));
        assert!(!is_geocentric_hypothetical_body(VULCAN));
        for body in [VULCAN, WHITE_MOON, PROSERPINA, WALDEMATH] {
            assert!(is_of_date_hypothetical_body(body), "body {body}");
        }
        assert!(!is_of_date_hypothetical_body(40));
        // Deflectors: Sun + Jupiter/Saturn barycenters only.
        assert_eq!(deflector_command(0), Some("10"));
        assert_eq!(deflector_command(5), Some("5"));
        assert_eq!(deflector_command(6), Some("6"));
        assert_eq!(deflector_command(4), None);
        assert_eq!(deflector_command(14), None);
    }

    #[test]
    fn strict_quoting_matches_urllib_safe_empty() {
        // Unreserved bytes pass; quotes/semicolons/spaces do not.
        assert_eq!(quote_strict("499"), "499");
        assert_eq!(quote_strict("'499'"), "%27499%27");
        assert_eq!(quote_strict("'2060;'"), "%272060%3B%27");
        assert_eq!(quote_strict("'500@0'"), "%27500%400%27");
        assert_eq!(quote_strict("'AU-D'"), "%27AU-D%27");
        assert_eq!(quote_strict("a b"), "a%20b");
        assert_eq!(quote_strict("~._-"), "~._-");
        // UTF-8 multibyte values encode per byte, uppercase hex.
        assert_eq!(quote_strict("é"), "%C3%A9");
    }

    #[test]
    fn tlist_formatting_is_python_compatible() {
        assert_eq!(format_tlist(2451545.0).unwrap(), "2451545.0");
        assert_eq!(format_tlist(2460000.5).unwrap(), "2460000.5");
        assert_eq!(format_tlist(2378496.123456).unwrap(), "2378496.123456");
        // Round-trips exactly.
        for jd in [2400000.5, 2451545.0, 2460676.7083333, 2689310.5] {
            let text = format_tlist(jd).unwrap();
            assert_eq!(text.parse::<f64>().unwrap().to_bits(), jd.to_bits());
        }
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                matches!(format_tlist(bad), Err(Error::NonFiniteArgument { .. })),
                "{bad} must be rejected"
            );
        }
    }

    #[test]
    fn vectors_url_carries_the_exact_release_contract() {
        let url = build_vectors_url("499", 2451545.0, "@0", HorizonsTimeType::Tdb).unwrap();
        assert_eq!(
            url,
            concat!(
                "https://ssd.jpl.nasa.gov/api/horizons.api",
                "?format=json",
                "&COMMAND=%27499%27",
                "&EPHEM_TYPE=%27VECTORS%27",
                "&CENTER=%27%400%27",
                "&TLIST=2451545.0",
                "&TLIST_TYPE=%27JD%27",
                "&VEC_TABLE=%272%27",
                "&OUT_UNITS=%27AU-D%27",
                "&VEC_CORR=%27NONE%27",
                "&CSV_FORMAT=%27YES%27",
                "&REF_SYSTEM=%27ICRF%27",
                "&REF_PLANE=%27FRAME%27",
                "&TIME_TYPE=%27TDB%27"
            )
        );
        // Small-body semicolon is encoded; UT spelling passes through.
        let chic = build_vectors_url("2060;", 2460000.5, "@0", HorizonsTimeType::Ut).unwrap();
        assert!(chic.contains("COMMAND=%272060%3B%27"), "{chic}");
        assert!(chic.contains("TLIST=2460000.5"), "{chic}");
        assert!(chic.contains("TIME_TYPE=%27UT%27"), "{chic}");
        // Non-default center is single-quoted and encoded like the release.
        let c = build_vectors_url("399", 2451545.0, "500@399", HorizonsTimeType::Tdb).unwrap();
        assert!(c.contains("CENTER=%27500%40399%27"), "{c}");
    }

    #[test]
    fn state_vector_constructor_rejects_non_finite_components() {
        let good = StateVector::new(1.0, 2.0, 3.0, 0.1, 0.2, 0.3).unwrap();
        assert_eq!(good.pos(), (1.0, 2.0, 3.0));
        assert_eq!(good.vel(), (0.1, 0.2, 0.3));
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(StateVector::new(bad, 0.0, 0.0, 0.0, 0.0, 0.0).is_err());
            assert!(StateVector::new(0.0, 0.0, 0.0, 0.0, 0.0, bad).is_err());
        }
    }

    /// Minimal Horizons vectors envelope around one CSV record.
    fn envelope(record: &str) -> String {
        format!(
            "{{\"signature\":{{\"version\":\"1.0\"}},\"result\":\"Dynamical point-mass ephemeris\\n$$SOE\\n{record}\\n$$EOE\\n\"}}"
        )
    }

    #[test]
    fn parse_accepts_the_release_record_shape() {
        let body = envelope(
            "2451545.000000000, A.D. 2000-Jan-01 12:00:00.0000, -1.234567890123456E-01, 2.0, 3.0, 4.0E-03, -5.0, 6.0",
        );
        let sv = parse_vectors_response(&body, "499").unwrap();
        assert_eq!(
            sv.pos(),
            (-0.1234567890123456, 2.0, 3.0),
            "fields 2-4 are the position"
        );
        assert_eq!(sv.vel(), (0.004, -5.0, 6.0), "fields 5-7 are the velocity");
    }

    #[test]
    fn parse_rejects_every_malformed_shape() {
        let good_record = "2451545.0, date, 1.0, 2.0, 3.0, 0.1, 0.2, 0.3";
        // API refusal payload (never retried): exact message preserved.
        let api = "{\"error\":\"unknown target '9999'\"}";
        match parse_vectors_response(api, "9999") {
            Err(Error::HorizonsApiError { command, message }) => {
                assert_eq!(command, "9999");
                assert_eq!(message, "unknown target '9999'");
            }
            other => panic!("API payload must be HorizonsApiError, got {other:?}"),
        }
        // Structural failures, each with its static reason.
        let cases = [
            ("{}", "missing result field"),
            ("{\"result\":123}", "expected string value"),
            ("{\"result\":\"no markers here\"}", "missing $$SOE marker"),
            ("{\"result\":\"$$EOE\\n$$SOE\\n1.0\"}", "empty data block"),
            (
                "{\"result\":\"prefix $$EOE middle $$SOE suffix\"}",
                "empty data block",
            ),
            ("{\"result\":\"$$SOE\\n   \\n$$EOE\"}", "empty data block"),
            (&envelope("2451545.0, date, 1.0, 2.0"), "short data record"),
            (
                &envelope("2451545.0, date, 1.0, 2.0, 3.0, 0.1, oops, 0.3"),
                "non-numeric state-vector field",
            ),
            (
                &envelope("2451545.0, date, 1.0, 2.0, 3.0, 0.1, NaN, 0.3"),
                "non-finite state-vector field",
            ),
        ];
        for (body, reason) in cases {
            match parse_vectors_response(body, "499") {
                Err(Error::HorizonsMalformed {
                    command,
                    reason: got,
                }) => {
                    assert_eq!(command, "499");
                    assert_eq!(got, reason, "wrong reason for {body}");
                }
                other => panic!("{body} must be HorizonsMalformed, got {other:?}"),
            }
        }
        // Stray numerics outside the first record never shift the indices:
        // the header line is skipped, the first data line wins.
        let stray = format!(
            "{{\"result\":\"epoch 2451545.0 scale 3\\n$$SOE\\n{good_record}\\n9.9, 9.9\\n$$EOE\"}}"
        );
        let sv = parse_vectors_response(&stray, "499").unwrap();
        assert_eq!(sv.pos(), (1.0, 2.0, 3.0));
        // Escaped quotes inside the envelope decode before marker search.
        let escaped = envelope(good_record).replace("date", "A.D. 2000-\\\"Jan\\\"-01");
        assert!(parse_vectors_response(&escaped, "499").is_ok());
        // Truncated JSON is malformed, never a partial vector.
        assert!(matches!(
            parse_vectors_response("{\"result\":\"$$SOE", "499"),
            Err(Error::HorizonsMalformed { .. })
        ));
    }

    #[test]
    fn retry_policy_matches_release_structure() {
        assert_eq!(retry_delay_secs(0), 0.5);
        assert_eq!(retry_delay_secs(1), 1.0);
        assert_eq!(HORIZONS_MAX_RETRIES, 2);
        // Retried: transport/timeout/status + malformed payloads.
        for err in [
            Error::HttpTimeout {
                purpose: "p".to_string(),
                seconds: 30.0,
            },
            Error::HttpTransport {
                purpose: "p".to_string(),
                detail: "d".to_string(),
            },
            Error::HttpStatus {
                status: 500,
                purpose: "p".to_string(),
            },
            Error::HorizonsMalformed {
                command: "499".to_string(),
                reason: "empty data block",
            },
        ] {
            assert!(is_retryable_horizons_error(&err), "{err:?} must retry");
        }
        // Never retried: API refusals and caller-side rejections.
        for err in [
            Error::HorizonsApiError {
                command: "499".to_string(),
                message: "nope".to_string(),
            },
            Error::HorizonsRequestFailed {
                command: "499".to_string(),
                attempts: 3,
                detail: "d".to_string(),
            },
            Error::NetworkSealed {
                purpose: "p".to_string(),
            },
            Error::HttpInvalidUrl {
                url: "x".to_string(),
            },
            Error::UnknownBody { body_id: 16 },
            Error::NonFiniteArgument { value: f64::NAN },
        ] {
            assert!(!is_retryable_horizons_error(&err), "{err:?} must not retry");
        }
    }

    #[test]
    fn client_config_is_checked() {
        let client = HorizonsClient::new().unwrap();
        assert_eq!(client.timeout_secs(), HORIZONS_DEFAULT_TIMEOUT_SECS);
        assert_eq!(client.max_retries(), HORIZONS_MAX_RETRIES);
        assert_eq!(HorizonsClient::default(), client);
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    HorizonsClient::with_timeout_secs(bad),
                    Err(Error::InvalidTimeout { .. })
                ),
                "{bad} must be InvalidTimeout"
            );
        }
    }

    #[test]
    fn horizons_config_carries_release_defaults() {
        let config = HorizonsConfig::default();
        assert_eq!(config.timeout_secs(), HORIZONS_DEFAULT_TIMEOUT_SECS);
        assert_eq!(
            config.max_cache_size().get(),
            HORIZONS_DEFAULT_MAX_CACHE_SIZE
        );
        assert_eq!(config.max_workers().get(), HORIZONS_DEFAULT_MAX_WORKERS);
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                matches!(
                    HorizonsConfig::new(
                        bad,
                        NonZeroUsize::new(16).unwrap(),
                        NonZeroUsize::new(2).unwrap()
                    ),
                    Err(Error::InvalidTimeout { .. })
                ),
                "{bad} must be InvalidTimeout"
            );
        }
    }

    #[test]
    fn client_from_config_round_trips_knobs() {
        let config = HorizonsConfig::new(
            5.0,
            NonZeroUsize::new(16).unwrap(),
            NonZeroUsize::new(2).unwrap(),
        )
        .unwrap();
        let client = HorizonsClient::from_config(config).unwrap();
        assert_eq!(client.timeout_secs(), 5.0);
        assert_eq!(client.max_cache_size(), 16);
        assert_eq!(client.max_workers(), 2);
        assert_eq!(client.max_retries(), HORIZONS_MAX_RETRIES);
        let defaulted = HorizonsClient::from_config(HorizonsConfig::default()).unwrap();
        assert_eq!(defaulted, HorizonsClient::new().unwrap());
    }

    /// Multi-shot loopback Horizons server: serves `responses` in order,
    /// records every raw request head for assertions, then exits.
    struct Fixture {
        url: String,
        received: std::sync::mpsc::Receiver<Vec<u8>>,
        _thread: std::thread::JoinHandle<()>,
    }

    fn serve_horizons(responses: Vec<(&'static str, String)>) -> Fixture {
        let listener = TcpListener::bind("127.0.0.1:0").expect("bind loopback");
        let addr = listener.local_addr().expect("loopback addr");
        let (tx, rx) = std::sync::mpsc::channel();
        let thread = std::thread::spawn(move || {
            for (status_line, body) in responses {
                let (mut stream, _) = listener.accept().expect("accept");
                stream
                    .set_read_timeout(Some(Duration::from_secs(10)))
                    .expect("read timeout");
                let mut raw = Vec::new();
                let mut one = [0u8; 1];
                loop {
                    match stream.read(&mut one) {
                        Ok(0) => break,
                        Ok(_) => {
                            raw.push(one[0]);
                            if raw.ends_with(b"\r\n\r\n") {
                                break;
                            }
                            if raw.len() > 65536 {
                                break;
                            }
                        }
                        Err(_) => break,
                    }
                }
                let _ = tx.send(raw);
                let head = format!(
                    "{status_line}\r\nContent-Length: {}\r\nContent-Type: application/json\r\n\r\n",
                    body.len()
                );
                let _ = stream.write_all(head.as_bytes());
                let _ = stream.write_all(body.as_bytes());
                let _ = stream.flush();
            }
        });
        Fixture {
            url: format!("http://{addr}/vectors"),
            received: rx,
            _thread: thread,
        }
    }

    fn recv_head(fixture: &Fixture) -> String {
        let raw = fixture
            .received
            .recv_timeout(Duration::from_secs(10))
            .expect("server saw the request");
        String::from_utf8_lossy(&raw).into_owned()
    }

    fn local_url(fixture: &Fixture, command: &str) -> String {
        // Same contract as build_vectors_url, pointed at loopback.
        let params = [
            ("format", "json".to_string()),
            ("COMMAND", format!("'{command}'")),
            ("TLIST", "2451545.0".to_string()),
        ];
        let query: Vec<String> = params
            .iter()
            .map(|(k, v)| format!("{k}={}", quote_strict(v)))
            .collect();
        format!("{}?{}", fixture.url, query.join("&"))
    }

    #[test]
    fn fetch_once_sends_release_headers_and_parses() {
        let record =
            "2451545.000000000, A.D. 2000-Jan-01 12:00:00.0000, 1.0, 2.0, 3.0, 0.1, 0.2, 0.3";
        let fixture = serve_horizons(vec![("HTTP/1.1 200 OK", envelope(record))]);
        let client = HorizonsClient::new().unwrap();
        let url = local_url(&fixture, "499");
        let sv = client
            .fetch_once(&url, "499", allow())
            .expect("loopback fetch");
        assert_eq!(sv.pos(), (1.0, 2.0, 3.0));
        assert_eq!(sv.vel(), (0.1, 0.2, 0.3));
        let head = recv_head(&fixture);
        assert!(head.starts_with("GET "), "{head}");
        assert!(head.contains("user-agent: libephemeris/1.0\r\n"), "{head}");
        assert!(head.contains("COMMAND=%27499%27"), "{head}");
    }

    #[test]
    fn api_refusal_returns_without_retry() {
        // Two servable responses: if the client retried, the second request
        // would consume the decoy. It must not.
        let fixture = serve_horizons(vec![
            (
                "HTTP/1.1 200 OK",
                "{\"error\":\"unknown target\"}".to_string(),
            ),
            ("HTTP/1.1 200 OK", envelope("1,2,3,4,5,6,7,8")),
        ]);
        let client = HorizonsClient::new().unwrap();
        let url = local_url(&fixture, "9999");
        match client.fetch_loop(&url, "9999", allow()) {
            Err(Error::HorizonsApiError { .. }) => {}
            other => panic!("API refusal must surface, got {other:?}"),
        }
        // Only one exchange happened: drain the first request head, then
        // the decoy response must stay untouched (no retry).
        let head = recv_head(&fixture);
        assert!(head.contains("COMMAND=%279999%27"), "{head}");
        assert!(
            fixture
                .received
                .recv_timeout(Duration::from_millis(300))
                .is_err(),
            "API refusal must not be retried"
        );
    }

    #[test]
    fn transient_status_is_retried_then_succeeds() {
        let record = "2451545.0, date, 1.0, 2.0, 3.0, 0.1, 0.2, 0.3";
        let fixture = serve_horizons(vec![
            ("HTTP/1.1 500 Internal Server Error", "boom".to_string()),
            ("HTTP/1.1 200 OK", envelope(record)),
        ]);
        let client = HorizonsClient::new().unwrap();
        let url = local_url(&fixture, "499");
        // fetch_once reports the raw status (transport never retries).
        assert!(matches!(
            client.fetch_once(&url, "499", allow()),
            Err(Error::HttpStatus { status: 500, .. })
        ));
        // fetch_vectors retries once and parses the second response.
        let sv = client
            .fetch_loop(&url, "499", allow())
            .expect("retry then parse");
        assert_eq!(sv.pos(), (1.0, 2.0, 3.0));
    }

    #[test]
    fn exhausted_budget_reports_attempts_and_cause() {
        let fixture = serve_horizons(vec![
            ("HTTP/1.1 500 Internal Server Error", "a".to_string()),
            ("HTTP/1.1 500 Internal Server Error", "b".to_string()),
            ("HTTP/1.1 500 Internal Server Error", "c".to_string()),
        ]);
        let client = HorizonsClient::new().unwrap();
        let url = local_url(&fixture, "499");
        match client.fetch_loop(&url, "499", allow()) {
            Err(Error::HorizonsRequestFailed {
                command,
                attempts,
                detail,
            }) => {
                assert_eq!(command, "499");
                assert_eq!(attempts, 3);
                assert!(detail.contains("500"), "{detail}");
            }
            other => panic!("exhausted budget must be HorizonsRequestFailed, got {other:?}"),
        }
    }

    #[test]
    fn sealed_policy_refuses_before_socket_work() {
        let client = HorizonsClient::new().unwrap();
        let sealed = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        assert_eq!(
            sealed.effective_network_policy(),
            EffectiveNetworkPolicy::Sealed
        );
        assert!(matches!(
            client.fetch_vectors(
                "499",
                2451545.0,
                "@0",
                HorizonsTimeType::Tdb,
                sealed.effective_network_policy()
            ),
            Err(Error::NetworkSealed { .. })
        ));
        let allowed = EphemerisContext::new().with_network_policy(Some(NetworkPolicy::Allow));
        // An unroutable port proves an allowed policy reaches the socket
        // (no live call: the loop never leaves loopback).
        match client.fetch_loop(
            "http://127.0.0.1:9/vectors",
            "499",
            allowed.effective_network_policy(),
        ) {
            Err(Error::NetworkSealed { .. }) => panic!("allow context must reach the socket"),
            Err(_) => {}
            Ok(_) => panic!("unroutable port must not succeed"),
        }
    }

    fn unit_vector(x: f64) -> StateVector {
        StateVector::new(x, 0.0, 0.0, 0.0, 0.0, 0.0).expect("finite fixture vector")
    }

    #[test]
    fn jd12_normalization_coalesces_and_rejects() {
        // At JD scale (~2.4e6) the f64 spacing (~4.7e-10) already exceeds
        // the 1e-12 grid, so every distinct f64 epoch keys distinctly and
        // the normalization is the identity there (verified below). The
        // 12-decimal coalescing itself is exercised at unit magnitude,
        // where sub-grid f64 steps are representable.
        let base = 2451545.0_f64;
        assert_eq!(round_jd12(base).unwrap(), base);
        assert_ne!(
            round_jd12(base).unwrap().to_bits(),
            round_jd12(base + 1e-9).unwrap().to_bits()
        );
        assert_eq!(
            round_jd12(1.0).unwrap().to_bits(),
            round_jd12(1.0 + 4e-13).unwrap().to_bits()
        );
        assert_ne!(
            round_jd12(1.0).unwrap().to_bits(),
            round_jd12(1.0 + 2e-12).unwrap().to_bits()
        );
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(
                matches!(round_jd12(bad), Err(Error::NonFiniteArgument { .. })),
                "{bad} must never be cached"
            );
        }
    }

    #[test]
    fn cache_key_identity_matches_release_tuple() {
        let jd = 2451545.0_f64;
        let key = HorizonsCacheKey::new("499", jd, "@0", HorizonsTimeType::Tdb).unwrap();
        // At JD scale every distinct f64 epoch keys distinctly (see the
        // normalization test); command/center/time-type never coalesce.
        assert_eq!(
            key,
            HorizonsCacheKey::new("499", jd, "@0", HorizonsTimeType::Tdb).unwrap()
        );
        assert_ne!(
            key,
            HorizonsCacheKey::new("499", jd + 1e-9, "@0", HorizonsTimeType::Tdb).unwrap()
        );
        assert_ne!(
            key,
            HorizonsCacheKey::new("599", jd, "@0", HorizonsTimeType::Tdb).unwrap()
        );
        assert_ne!(
            key,
            HorizonsCacheKey::new("499", jd, "@399", HorizonsTimeType::Tdb).unwrap()
        );
        assert_ne!(
            key,
            HorizonsCacheKey::new("499", jd, "@0", HorizonsTimeType::Ut).unwrap()
        );
        assert_eq!(key.jd_key().to_bits(), round_jd12(jd).unwrap().to_bits());
        assert!(HorizonsCacheKey::new("499", f64::NAN, "@0", HorizonsTimeType::Tdb).is_err());
    }

    #[test]
    fn lru_cache_refreshes_recency_and_evicts_oldest() {
        let mut cache = HorizonsCache::new(2);
        let jd = 2451545.0_f64;
        let ka = HorizonsCacheKey::new("499", jd, "@0", HorizonsTimeType::Tdb).unwrap();
        let kb = HorizonsCacheKey::new("599", jd, "@0", HorizonsTimeType::Tdb).unwrap();
        let kc = HorizonsCacheKey::new("699", jd, "@0", HorizonsTimeType::Tdb).unwrap();
        cache.insert(ka.clone(), unit_vector(1.0));
        cache.insert(kb.clone(), unit_vector(2.0));
        // Hit on `ka` refreshes recency: `kb` evicts next.
        assert_eq!(cache.lookup(&ka), Some(unit_vector(1.0)));
        cache.insert(kc.clone(), unit_vector(3.0));
        assert!(cache.contains(&ka));
        assert!(cache.contains(&kc));
        assert!(!cache.contains(&kb));
        assert_eq!(cache.len(), 2);
        // Overwrite refreshes recency and value without growing.
        cache.insert(kb.clone(), unit_vector(20.0));
        cache.insert(ka.clone(), unit_vector(10.0));
        assert_eq!(cache.len(), 2);
        assert_eq!(cache.lookup(&ka), Some(unit_vector(10.0)));
        assert!(cache.contains(&kb));
        assert!(!cache.contains(&kc));
        // Zero capacity admits no entry (clients reject this shape via NonZeroUsize).
        let mut zero = HorizonsCache::new(0);
        zero.insert(ka.clone(), unit_vector(1.0));
        assert!(zero.is_empty());
        // contains performs no recency change: covered by the eviction order above.
    }

    #[test]
    fn get_or_fetch_caches_success_never_failure() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let client = HorizonsClient::new().unwrap();
        let jd = 2451545.0_f64;
        let calls = AtomicUsize::new(0);
        let first = client
            .get_or_fetch_with("499", jd, "@0", HorizonsTimeType::Tdb, || {
                calls.fetch_add(1, Ordering::SeqCst);
                Ok(unit_vector(1.0))
            })
            .unwrap();
        assert_eq!(first, unit_vector(1.0));
        // Hit performs no second fetch (same epoch spelling).
        let second = client
            .get_or_fetch_with("499", jd, "@0", HorizonsTimeType::Tdb, || {
                calls.fetch_add(1, Ordering::SeqCst);
                Ok(unit_vector(2.0))
            })
            .unwrap();
        assert_eq!(second, unit_vector(1.0));
        assert_eq!(calls.load(Ordering::SeqCst), 1);
        // Failures are never cached: the retry runs the fetcher again.
        let fails = AtomicUsize::new(0);
        assert!(matches!(
            client.get_or_fetch_with("599", jd, "@0", HorizonsTimeType::Tdb, || {
                fails.fetch_add(1, Ordering::SeqCst);
                Err::<StateVector, _>(Error::HorizonsMalformed {
                    command: "599".to_string(),
                    reason: "boom",
                })
            }),
            Err(Error::HorizonsMalformed { .. })
        ));
        assert!(
            client
                .get_or_fetch_with("599", jd, "@0", HorizonsTimeType::Tdb, || {
                    fails.fetch_add(1, Ordering::SeqCst);
                    Ok(unit_vector(5.0))
                })
                .is_ok()
        );
        assert_eq!(fails.load(Ordering::SeqCst), 2);
        assert!(
            client
                .cache_contains("599", jd, "@0", HorizonsTimeType::Tdb)
                .unwrap()
        );
        assert!(
            !client
                .cache_contains("599", jd, "@0", HorizonsTimeType::Ut)
                .unwrap()
        );
    }

    #[test]
    fn cached_hit_bypasses_sealed_gate_like_release() {
        let client = HorizonsClient::new().unwrap();
        let jd = 2451545.0_f64;
        client
            .get_or_fetch_with("499", jd, "@0", HorizonsTimeType::Tdb, || {
                Ok(unit_vector(1.0))
            })
            .unwrap();
        let sealed = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        assert_eq!(
            sealed.effective_network_policy(),
            EffectiveNetworkPolicy::Sealed
        );
        // A cached hit serves without socket work, so the
        // sealed gate (which fires inside the transport) is not reached.
        // Sharing a fetching client with a sealed context therefore serves
        // stale vectors — callers needing sealed purity must not share.
        let hit = client
            .fetch_vectors_cached(
                "499",
                jd,
                "@0",
                HorizonsTimeType::Tdb,
                sealed.effective_network_policy(),
            )
            .expect("cached hit bypasses the gate");
        assert_eq!(hit, unit_vector(1.0));
        // A miss still refuses before any socket work.
        assert!(matches!(
            client.fetch_vectors_cached(
                "599",
                jd,
                "@0",
                HorizonsTimeType::Tdb,
                sealed.effective_network_policy()
            ),
            Err(Error::NetworkSealed { .. })
        ));
    }

    #[test]
    fn batch_serves_cache_first_captures_errors_and_stores_success() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let client = HorizonsClient::new().unwrap();
        let jd_hit = 2451545.0_f64;
        let jd_ok = 2451546.0_f64;
        let jd_bad = 2451547.0_f64;
        client
            .get_or_fetch_with("499", jd_hit, "@0", HorizonsTimeType::Tdb, || {
                Ok(unit_vector(1.0))
            })
            .unwrap();
        let requests = [
            VectorsRequest::new("499", jd_hit, "@0").unwrap(),
            VectorsRequest::new("599", jd_ok, "@0").unwrap(),
            VectorsRequest::new("699", jd_bad, "@0").unwrap(),
        ];
        let calls = AtomicUsize::new(0);
        let (results, errors) = client.fetch_batch_with(&requests, |command, jd, _| {
            calls.fetch_add(1, Ordering::SeqCst);
            assert_eq!(jd, if command == "599" { jd_ok } else { jd_bad });
            if command == "599" {
                Ok(unit_vector(2.0))
            } else {
                Err(Error::HorizonsMalformed {
                    command: command.to_string(),
                    reason: "bad record",
                })
            }
        });
        // Cached entry served without fetch; results keyed by the original
        // request triple.
        assert_eq!(calls.load(Ordering::SeqCst), 2);
        assert_eq!(results.len(), 2);
        assert_eq!(
            results[&requests[0].batch_key()],
            unit_vector(1.0),
            "cache-first hit"
        );
        assert_eq!(results[&requests[1].batch_key()], unit_vector(2.0));
        assert_eq!(errors.len(), 1);
        assert!(matches!(
            errors[&requests[2].batch_key()],
            Error::HorizonsMalformed { .. }
        ));
        // Successes are stored; failures never are.
        assert!(
            client
                .cache_contains("599", jd_ok, "@0", HorizonsTimeType::Tdb)
                .unwrap()
        );
        assert!(
            !client
                .cache_contains("699", jd_bad, "@0", HorizonsTimeType::Tdb)
                .unwrap()
        );
        // A failed request retries its fetch on the next batch.
        let (retry_results, retry_errors) =
            client.fetch_batch_with(&requests[2..], |_, _, _| Ok(unit_vector(3.0)));
        assert_eq!(retry_results[&requests[2].batch_key()], unit_vector(3.0));
        assert!(retry_errors.is_empty());
        assert!(
            client
                .cache_contains("699", jd_bad, "@0", HorizonsTimeType::Tdb)
                .unwrap()
        );
    }

    #[test]
    fn batch_keys_keep_original_spelling_and_tdb_only() {
        let client = HorizonsClient::new().unwrap();
        // Results report under the requested spelling: two distinct epochs each get their own row and cache entry.
        let jd_a = 2451545.0_f64;
        let jd_b = 2451545.0_f64 + 1e-9;
        let requests = [
            VectorsRequest::new("499", jd_a, "@0").unwrap(),
            VectorsRequest::new("499", jd_b, "@0").unwrap(),
        ];
        assert_ne!(requests[0].batch_key(), requests[1].batch_key());
        assert_ne!(
            requests[0].cache_key(HorizonsTimeType::Tdb).unwrap(),
            requests[1].cache_key(HorizonsTimeType::Tdb).unwrap()
        );
        let (results, errors) = client.fetch_batch_with(&requests, |_, _, _| Ok(unit_vector(7.0)));
        assert!(errors.is_empty());
        assert_eq!(results.len(), 2);
        assert_eq!(results[&requests[0].batch_key()], unit_vector(7.0));
        assert_eq!(results[&requests[1].batch_key()], unit_vector(7.0));
        assert_eq!(results[&requests[0].batch_key()].pos(), (7.0, 0.0, 0.0));
        assert_eq!(requests[0].batch_key().jd().to_bits(), jd_a.to_bits());
        // TDB-only: the batch wrote the TDB key, never the UT key.
        assert!(
            client
                .cache_contains("499", jd_a, "@0", HorizonsTimeType::Tdb)
                .unwrap()
        );
        assert!(
            !client
                .cache_contains("499", jd_a, "@0", HorizonsTimeType::Ut)
                .unwrap()
        );
        assert_eq!(requests[0].command(), "499");
        assert_eq!(requests[0].center(), "@0");
        assert!(VectorsRequest::new("499", f64::NAN, "@0").is_err());
    }

    #[test]
    fn client_defaults_cache_refresh_and_worker_bound() {
        use std::num::NonZeroUsize;
        let client = HorizonsClient::new().unwrap();
        assert_eq!(client.max_cache_size(), HORIZONS_DEFAULT_MAX_CACHE_SIZE);
        assert_eq!(client.max_cache_size(), 4096);
        assert_eq!(client.max_workers(), HORIZONS_DEFAULT_MAX_WORKERS);
        assert_eq!(client.max_workers(), 8);
        assert_eq!(client.cache_len(), 0);
        client
            .get_or_fetch_with("499", 2451545.0, "@0", HorizonsTimeType::Tdb, || {
                Ok(unit_vector(1.0))
            })
            .unwrap();
        assert_eq!(client.cache_len(), 1);
        // Snapshot clone carries the cached vectors; later inserts stay local.
        let snapshot = client.clone();
        assert_eq!(snapshot, client);
        client
            .get_or_fetch_with("599", 2451545.0, "@0", HorizonsTimeType::Tdb, || {
                Ok(unit_vector(2.0))
            })
            .unwrap();
        assert_ne!(snapshot, client);
        assert_eq!(snapshot.cache_len(), 1);
        // Refresh drops every entry but keeps the client usable.
        client.clear_cache();
        assert_eq!(client.cache_len(), 0);
        client.shutdown();
        client
            .get_or_fetch_with("499", 2451545.0, "@0", HorizonsTimeType::Tdb, || {
                Ok(unit_vector(9.0))
            })
            .unwrap();
        assert_eq!(client.cache_len(), 1);
        // Explicit small cache evicts past capacity.
        let tiny = HorizonsClient::with_config(
            HORIZONS_DEFAULT_TIMEOUT_SECS,
            NonZeroUsize::new(1).unwrap(),
            NonZeroUsize::new(2).unwrap(),
        )
        .unwrap();
        assert_eq!(tiny.max_workers(), 2);
        tiny.get_or_fetch_with("499", 2451545.0, "@0", HorizonsTimeType::Tdb, || {
            Ok(unit_vector(1.0))
        })
        .unwrap();
        tiny.get_or_fetch_with("599", 2451545.0, "@0", HorizonsTimeType::Tdb, || {
            Ok(unit_vector(2.0))
        })
        .unwrap();
        assert_eq!(tiny.cache_len(), 1);
        assert!(
            !tiny
                .cache_contains("499", 2451545.0, "@0", HorizonsTimeType::Tdb)
                .unwrap()
        );
        assert!(
            tiny.cache_contains("599", 2451545.0, "@0", HorizonsTimeType::Tdb)
                .unwrap()
        );
    }

    #[test]
    fn unroutable_host_is_typed_transport_failure_and_never_cached() {
        // Port 9 (discard) on loopback refuses connections, exercising the
        // real single-exchange path against unavailable network: the
        // failure must surface as HttpTransport (the URLError surface the
        // release callers catch), and failures — including sealed-policy
        // refusals through the cached path — must never populate the LRU.
        let client = HorizonsClient::new().unwrap();
        match client.fetch_once("http://127.0.0.1:9/api/horizons.api", "499", allow()) {
            Err(Error::HttpTransport { purpose, .. }) => {
                assert_eq!(purpose, "Horizons state vector for 499");
            }
            other => panic!("refused connection must be HttpTransport, got {other:?}"),
        }
        assert_eq!(client.cache_len(), 0);
        let jd = 2451545.0_f64;
        let sealed = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        assert_eq!(
            sealed.effective_network_policy(),
            EffectiveNetworkPolicy::Sealed
        );
        assert!(matches!(
            client.fetch_vectors_cached(
                "599",
                jd,
                "@0",
                HorizonsTimeType::Tdb,
                sealed.effective_network_policy()
            ),
            Err(Error::NetworkSealed { .. })
        ));
        assert_eq!(client.cache_len(), 0, "sealed refusal must not be cached");
        // Repeating the refused query still attempts (still refuses):
        // there is no negative caching to go stale.
        assert!(matches!(
            client.fetch_vectors_cached(
                "599",
                jd,
                "@0",
                HorizonsTimeType::Tdb,
                sealed.effective_network_policy()
            ),
            Err(Error::NetworkSealed { .. })
        ));
        assert_eq!(client.cache_len(), 0);
    }
}
